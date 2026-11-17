"""Claim verifier: the part that makes the report trustworthy.

A claim survives only if:
  1. it cites at least one evidence id,
  2. every cited id exists,
  3. it is either a verbatim copy of a cited record's `detail`, or every standalone
     number in it equals a structured value (value / extras) of a cited record when
     rounded to the number of decimals the claim itself uses.

Numbers are never harvested from free text in the evidence: boilerplate such as the
"+1" in log(CPM+1) would otherwise whitelist unrelated numbers. Dropped claims are
counted and surfaced, never silently patched.
"""
from __future__ import annotations

import re

from .models import Claim, EvidenceItem

_PSEUDOCOUNT = re.compile(r"\(?\s*CPM\s*\+\s*1\s*\)?", re.IGNORECASE)
_NUM = re.compile(r"(?<![A-Za-z0-9_.])[-+]?\d+(?:\.\d+)?(?![A-Za-z0-9_])")


def _numbers(text: str) -> list[str]:
    text = _PSEUDOCOUNT.sub(" ", text)
    return [m.group() for m in _NUM.finditer(text)]


def _allowed(items: list[EvidenceItem]) -> list[float]:
    nums: list[float] = []
    for it in items:
        if it.value is not None:
            nums += [it.value, abs(it.value)]
        nums += [v for k, v in it.extras.items() if not k.startswith("_")]  # '_' = internal flag, not quotable
    return nums


def _matches(token: str, allowed: list[float]) -> bool:
    decimals = len(token.split(".")[1]) if "." in token else 0
    n = float(token)
    return any(round(a, decimals) == n or round(abs(a), decimals) == abs(n) for a in allowed)


def _norm(s: str) -> str:
    return " ".join(s.split())


def check_claim(claim: Claim, evidence: dict[str, EvidenceItem]) -> tuple[bool, str]:
    if not claim.evidence_ids:
        return False, "no evidence cited"
    unknown = [e for e in claim.evidence_ids if e not in evidence]
    if unknown:
        return False, f"unknown evidence id(s): {unknown}"
    cited = [evidence[e] for e in claim.evidence_ids]
    if any(_norm(claim.text) == _norm(c.detail) for c in cited):
        return True, "verbatim"
    allowed = _allowed(cited)
    for tok in _numbers(claim.text):
        if not _matches(tok, allowed):
            return False, f"number {tok} not found in cited evidence"
    return True, "ok"


def verify(claims: list[Claim], items: list[EvidenceItem]) -> tuple[list[Claim], list[tuple[Claim, str]]]:
    index = {i.id: i for i in items}
    kept: list[Claim] = []
    dropped: list[tuple[Claim, str]] = []
    for c in claims:
        ok, why = check_claim(c, index)
        if ok:
            kept.append(c)
        else:
            dropped.append((c, why))
    return kept, dropped
