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
            notes.append("Predicted direction differs between conditions (opposite effects by cell type).")
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

    for idx, h in enumerate(hyps):
        runner_up = hyps[1].score if idx == 0 and len(hyps) > 1 else None
        h.confidence = _confidence(idx, h, pip, runner_up)
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
