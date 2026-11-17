"""Data models shared by providers, the aggregator, the verifier and the agent."""
from __future__ import annotations

from enum import Enum
from typing import Optional

from pydantic import BaseModel, Field


class Direction(str, Enum):
    UP = "up"
    DOWN = "down"
    UNKNOWN = "unknown"


class Confidence(str, Enum):
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


class VariantQuery(BaseModel):
    rsid: str


class EvidenceItem(BaseModel):
    """One retrievable fact. Every claim in a report must point at one or more of these."""

    id: str
    source: str  # "ag2_precomputed" | "open_targets" | ...
    kind: str  # "variant_effect" | "credible_set"
    gene: Optional[str] = None
    cell_type: Optional[str] = None
    # Meaning depends on kind: log fold change for variant_effect, PIP for credible_set.
    value: Optional[float] = None
    detail: str = ""
    # Secondary numbers that a claim may legitimately quote (e.g. predicted expression).
    extras: dict[str, float] = Field(default_factory=dict)


class GeneHypothesis(BaseModel):
    gene: str
    score: float
    direction: Direction
    top_cell_types: list[str]
    n_conditions_called: int
    confidence: Confidence
    evidence_ids: list[str]
    notes: list[str] = Field(default_factory=list)


class Claim(BaseModel):
    text: str
    evidence_ids: list[str]


class Report(BaseModel):
    query: VariantQuery
    hypotheses: list[GeneHypothesis]
    claims: list[Claim] = Field(default_factory=list)
    caveats: list[str] = Field(default_factory=list)
    dropped_claims: int = 0
    evidence: list[EvidenceItem] = Field(default_factory=list)
