"""Orchestration: fetch evidence -> rank deterministically -> synthesise claims -> verify.

The LLM only writes prose claims over evidence that was already retrieved and ranked.
It never decides the ranking, and anything it writes is checked by the verifier.
"""
from __future__ import annotations

import json
from typing import Callable, Optional

from .aggregator import rank_genes, standing_caveats
from .models import Claim, EvidenceItem, GeneHypothesis, Report, VariantQuery
from .providers.base import EvidenceProvider
from .verifier import verify

# A synthesizer takes (hypotheses, evidence) and returns candidate claims.
Synthesizer = Callable[[list[GeneHypothesis], list[EvidenceItem]], list[Claim]]


def template_synthesizer(hyps: list[GeneHypothesis], items: list[EvidenceItem]) -> list[Claim]:
    """No-LLM baseline. Also the control the LLM synthesizer is compared against."""
    index = {i.id: i for i in items}
    claims: list[Claim] = []
    for h in hyps[:3]:
        top_ev = next((index[e] for e in h.evidence_ids if index[e].kind == "variant_effect"), None)
        if top_ev is None:
            continue
        claims.append(Claim(text=top_ev.detail, evidence_ids=[top_ev.id]))
        for e in h.evidence_ids:
            if index[e].kind == "credible_set":
                claims.append(Claim(text=index[e].detail, evidence_ids=[e]))
                break
    return claims


SYSTEM_PROMPT = """You write short factual claims about a genetic variant using ONLY the evidence records given.
Rules:
- Return JSON: {"claims": [{"text": str, "evidence_ids": [str]}]}.
- Every claim must cite one or more evidence ids from the list.
- Any number you write must appear in a cited record. Do not round differently or compute new numbers.
- Do not mention genes, cell types or mechanisms that are not in the records.
- Say 'predicted' for model outputs. Do not claim causation.
"""


def llm_synthesizer(client, model: str) -> Synthesizer:
    """client: an anthropic.Anthropic instance (or anything with .messages.create)."""

    def run(hyps: list[GeneHypothesis], items: list[EvidenceItem]) -> list[Claim]:
        payload = {
            "ranked_genes": [h.model_dump(mode="json") for h in hyps[:3]],
            "evidence": [i.model_dump(mode="json") for i in items],
        }
        msg = client.messages.create(
            model=model,
            max_tokens=1200,
            system=SYSTEM_PROMPT,
            messages=[{"role": "user", "content": json.dumps(payload)}],
        )
        text = "".join(b.text for b in msg.content if getattr(b, "type", "") == "text")
        try:
            data = json.loads(text[text.index("{") : text.rindex("}") + 1])
            return [Claim(**c) for c in data["claims"]]
        except (ValueError, KeyError, TypeError):
            return []

    return run


def run_agent(
    query: VariantQuery,
    providers: list[EvidenceProvider],
    synthesizer: Optional[Synthesizer] = None,
    context: Optional[list[str]] = None,
) -> Report:
    """context: optional cell-type / tissue terms (e.g. ["hepatocyte"]). Variant-effect evidence is
    restricted to conditions whose name contains one of them, so ranking answers "which gene moves in
    the relevant cell type" instead of "in any of ~500 conditions"."""
    items: list[EvidenceItem] = []
    for p in providers:
        items += p.fetch(query)
    extra_caveats: list[str] = []
    if context:
        terms = [t.lower() for t in context]
        before = sum(1 for i in items if i.kind == "variant_effect")
        items = [i for i in items if i.kind != "variant_effect" or any(t in (i.cell_type or "").lower() for t in terms)]
        extra_caveats.append(f"Ranking restricted to conditions matching {context} "
                             f"({sum(1 for i in items if i.kind == 'variant_effect')} of {before} variant-effect records).")
    if not any(i.kind == "variant_effect" for i in items):
        for p in providers:
            st = getattr(p, "status", lambda r: "unknown")(query.rsid)
            if st == "not_in_panel":
                extra_caveats.append("Variant is not in the AG2 variant panel: it was never scored.")
            elif st == "scored_no_effect":
                extra_caveats.append("Variant was scored by AG2 but no protein-coding gene effect exceeded the release filter (|score| > 0.05).")
    hyps = rank_genes(items)
    candidates = (synthesizer or template_synthesizer)(hyps, items)
    kept, dropped = verify(candidates, items)
    return Report(
        query=query,
        hypotheses=hyps,
        claims=kept,
        caveats=standing_caveats(items, hyps) + extra_caveats,
        dropped_claims=len(dropped),
        evidence=items,
    )
