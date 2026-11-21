"""Simple baselines the agent has to beat to be worth anything."""
from __future__ import annotations

from .models import EvidenceItem


def l2g_top_gene(items: list[EvidenceItem]) -> tuple[str, float] | None:
    """Open Targets' own best gene for this variant: highest locus-to-gene score over its credible sets.

    Not trait-specific: a variant in credible sets for several traits gets the best score across them."""
    l2g = [i for i in items if i.kind == "l2g" and i.gene and i.value is not None]
    if not l2g:
        return None
    best = max(l2g, key=lambda i: i.value)
    return best.gene, best.value
