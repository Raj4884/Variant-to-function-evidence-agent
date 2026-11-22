"""Deterministic evidence aggregation. No LLM in this module.

The ranking and confidence rules are heuristics. They are written down here so they
can be evaluated and changed; they are not calibrated probabilities. The eval harness
(eval/run_eval.py) is where they get tested.
"""
from __future__ import annotations

from collections import defaultdict

from .models import Confidence, Direction, EvidenceItem, GeneHypothesis

HIGH_PIP = 0.9
HIGH_EFFECT = 0.2
DOMINANCE_RATIO = 2.0  # top gene must beat the runner-up by this factor for HIGH
# If the runner-up scores at least this fraction of the top gene, the locus is treated as multi-gene:
# the top gene is not singled out and confidence is capped at LOW. Heuristic, not calibrated.
AMBIGUITY_RATIO = 0.75
# Open Targets L2G counts as "preferring another gene" only if its score is at least this high.
L2G_DISAGREE_MIN = 0.5
# If the other L2G gene beats ours by less than this, the sources are near-tied: treat as multi-gene (cap LOW).
L2G_TIE_MARGIN = 0.1
# Best fine-mapping PIP below this: the variant is probably a tag for another causal variant (cap LOW).
LOW_PIP = 0.1


def _direction(lfc: float) -> Direction:
    if lfc > 0:
        return Direction.UP
    if lfc < 0:
        return Direction.DOWN
    return Direction.UNKNOWN


def best_pip(items: list[EvidenceItem]) -> tuple[float | None, list[str]]:
    pips = [(i.value, i.id) for i in items if i.kind == "credible_set" and i.value is not None]
    if not pips:
        return None, []
    top = max(p for p, _ in pips)
    return top, [i for p, i in pips if p == top]


def rank_genes(items: list[EvidenceItem]) -> list[GeneHypothesis]:
    effects = [i for i in items if i.kind == "variant_effect" and i.gene and i.value is not None]
    called = [i for i in effects if i.extras.get("_called", 0.0) == 1.0]
    by_gene: dict[str, list[EvidenceItem]] = defaultdict(list)
    for i in called:
        by_gene[i.gene].append(i)

    pip, pip_ids = best_pip(items)
    l2g = [i for i in items if i.kind == "l2g" and i.gene and i.value is not None]
    l2g_best = max(l2g, key=lambda i: i.value) if l2g else None
    hyps: list[GeneHypothesis] = []
    for gene, evs in by_gene.items():
        evs.sort(key=lambda e: abs(e.value), reverse=True)
        top = evs[0]
        # Score: strongest single-condition effect. Breadth is reported, not rewarded,
        # because the paper finds most real calls are cell-type specific.
        score = abs(top.value)
        # Direction comes from the strongest condition; flag if conditions disagree.
        signs = {_direction(e.value) for e in evs}
        notes = []
        if len(signs) > 1:
            same = sum(1 for e in evs if _direction(e.value) == _direction(top.value))
            label = "Direction is mixed" if same / len(evs) < 0.8 else "Minority of conditions disagree on direction"
            notes.append(
                f"{label}: {same} of {len(evs)} called conditions share the top direction "
                f"({_direction(top.value).value}); the rest are opposite."
            )
        hyps.append(
            GeneHypothesis(
                gene=gene,
                score=score,
                direction=_direction(top.value),
                top_cell_types=[e.cell_type for e in evs[:3] if e.cell_type],
                n_conditions_called=len(evs),
                confidence=Confidence.LOW,
                evidence_ids=[e.id for e in evs[:5]] + pip_ids,
                notes=notes,
            )
        )
    hyps.sort(key=lambda h: h.score, reverse=True)

    l2g_other = bool(hyps and l2g_best is not None and l2g_best.gene != hyps[0].gene and l2g_best.value >= L2G_DISAGREE_MIN)
    top_l2g = max((i.value for i in l2g if hyps and i.gene == hyps[0].gene), default=0.0)
    l2g_clear = l2g_other and (l2g_best.value - top_l2g) >= L2G_TIE_MARGIN   # other gene clearly preferred
    l2g_tied = l2g_other and not l2g_clear                                   # near-tie between two genes
    for idx, h in enumerate(hyps):
        runner_up = hyps[1].score if idx == 0 and len(hyps) > 1 else None
        h.confidence = _confidence(idx, h, pip, runner_up)
    if hyps:
        top = hyps[0]
        if l2g_clear and top.confidence == Confidence.HIGH:
            top.confidence = Confidence.MEDIUM   # two sources disagree on the gene: never HIGH
        if l2g_tied:
            top.confidence = Confidence.LOW
            top.notes.append(
                f"Near-tie between {top.gene} and {l2g_best.gene} in Open Targets L2G "
                f"({top_l2g:.2f} vs {l2g_best.value:.2f}): treat as a multi-gene locus."
            )
        if pip is not None and pip < LOW_PIP:
            top.confidence = Confidence.LOW
            top.notes.append(
                f"Fine-mapping gives this variant a low posterior probability (best {pip:.2f}): "
                "it may only tag a different causal variant, so the variant-level claim is weak."
            )
    if len(hyps) > 1 and hyps[1].score >= AMBIGUITY_RATIO * hyps[0].score:
        close = [h for h in hyps[1:] if h.score >= AMBIGUITY_RATIO * hyps[0].score]
        hyps[0].confidence = Confidence.LOW
        hyps[0].notes.append(
            "Multi-gene locus: " + ", ".join(f"{h.gene} ({h.score:.2f})" for h in close)
            + f" score within {AMBIGUITY_RATIO:.0%} of {hyps[0].gene} ({hyps[0].score:.2f}); "
            "effect size alone does not identify the causal gene."
        )
    if hyps and l2g_best is not None:
        top = hyps[0]
        top_l2g = max((i.value for i in l2g if i.gene == top.gene), default=None)
        if l2g_best.gene == top.gene:
            top.notes.append(f"Open Targets L2G also ranks {top.gene} first ({l2g_best.value:.2f}).")
        else:
            seen = f" (L2G score for {top.gene}: {top_l2g:.2f})" if top_l2g is not None else ""
            top.notes.append(f"Open Targets L2G prefers {l2g_best.gene} ({l2g_best.value:.2f}), a different gene{seen}.")
        top.evidence_ids.append(l2g_best.id)
    return hyps


def _confidence(idx: int, h: GeneHypothesis, pip: float | None, runner_up: float | None) -> Confidence:
    if idx != 0:
        return Confidence.LOW
    dominant = runner_up is None or h.score >= DOMINANCE_RATIO * runner_up
    if pip is not None and pip >= HIGH_PIP and h.score >= HIGH_EFFECT and dominant:
        return Confidence.HIGH
    if h.score >= HIGH_EFFECT or (pip is not None and pip >= HIGH_PIP):
        return Confidence.MEDIUM
    return Confidence.LOW


def standing_caveats(items: list[EvidenceItem], hyps: list[GeneHypothesis]) -> list[str]:
    caveats = [
        "AlphaGenome 2 recovers only a minority of fine-mapped causal variants at high precision; "
        "a missing call is not evidence of no effect.",
        "Predictions are for healthy, steady-state conditions and are not validated for clinical use.",
    ]
    if not any(i.kind == "credible_set" for i in items):
        caveats.append("No fine-mapping (credible set) evidence was found, so causality is untested here.")
    if not hyps:
        caveats.append("No gene passed AG2's call thresholds in any condition.")
    return caveats
