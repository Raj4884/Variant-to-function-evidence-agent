from __future__ import annotations

from typing import Protocol

from ..models import EvidenceItem, VariantQuery


class EvidenceProvider(Protocol):
    name: str

    def fetch(self, query: VariantQuery) -> list[EvidenceItem]:
        """Return evidence for the variant. Must not raise for 'no data'; return []."""
        ...
