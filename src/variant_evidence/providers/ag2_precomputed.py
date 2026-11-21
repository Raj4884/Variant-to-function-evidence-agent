"""Reader for AlphaGenome 2 precomputed variant scores.

The AG2 paper reports a precomputed table of single-cell expression variant effects
(Supplementary Table 15 / downloads page). I have not inspected the real download
format, so this reader expects a normalised CSV and the column mapping is explicit.
Once the real file is downloaded, write a small converter to this schema (or change
COLUMN_MAP). Nothing else in the project depends on the raw format.

Normalised columns:
    variant    rsID, e.g. rs2836882
    gene       gene symbol
    condition  cell type / tissue condition name
    lfc        predicted natural-log fold change (ALT vs REF) on gene expression
    log_cpm    log(predicted CPM + 1) of the gene in that condition
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from ..models import EvidenceItem, VariantQuery

COLUMN_MAP = {
    "variant": "variant",
    "gene": "gene",
    "condition": "condition",
    "lfc": "lfc",
    "log_cpm": "log_cpm",
}

# "Called" thresholds the AG2 paper calibrated on TenK10K sc-eQTLs (Methods, Definition of called variants).
LFC_THRESHOLD = 0.1
LOG_CPM_THRESHOLD = 1.5


class AG2PrecomputedProvider:
    name = "ag2_precomputed"

    def __init__(self, csv_path: str | Path, scored_rsids: set[str] | None = None):
        """scored_rsids: rsIDs confirmed present in the AG2 variant panel (optional). Lets the agent tell
        "not scored at all" from "scored, no gene effect above the release filter"."""
        self.scored_rsids = scored_rsids
        self.csv_path = Path(csv_path)
        self._index: dict[str, list[dict]] = defaultdict(list)
        self._load()

    def _load(self) -> None:
        with self.csv_path.open(newline="") as fh:
            reader = csv.DictReader(fh)
            missing = set(COLUMN_MAP.values()) - set(reader.fieldnames or [])
            if missing:
                raise ValueError(f"{self.csv_path} is missing columns: {sorted(missing)}")
            for row in reader:
                self._index[row[COLUMN_MAP["variant"]]].append(row)

    def status(self, rsid: str) -> str:
        if rsid in self._index:
            return "has_rows"
        if self.scored_rsids is None:
            return "unknown"
        return "scored_no_effect" if rsid in self.scored_rsids else "not_in_panel"

    def fetch(self, query: VariantQuery) -> list[EvidenceItem]:
        items: list[EvidenceItem] = []
        for row in self._index.get(query.rsid, []):
            lfc = float(row[COLUMN_MAP["lfc"]])
            log_cpm = float(row[COLUMN_MAP["log_cpm"]])
            called = abs(lfc) > LFC_THRESHOLD and log_cpm > LOG_CPM_THRESHOLD
            gene = row[COLUMN_MAP["gene"]]
            cond = row[COLUMN_MAP["condition"]]
            extras = {"log_cpm": log_cpm, "_called": 1.0 if called else 0.0}
            chan = ""
            if row.get("n_channels"):
                n_ch, n_ag = int(float(row["n_channels"])), int(float(row["n_agree"]))
                extras["n_channels"], extras["n_agree"] = float(n_ch), float(n_ag)
                chan = f" Median over {n_ch} assay channels; {n_ag} agree in direction."
            items.append(
                EvidenceItem(
                    id=f"ag2:{query.rsid}:{gene}:{cond}",
                    source=self.name,
                    kind="variant_effect",
                    gene=gene,
                    cell_type=cond,
                    value=lfc,
                    detail=(
                        f"AG2 predicts {query.rsid} changes {gene} expression in {cond} "
                        f"by log fold change {lfc:+.2f} (predicted log(CPM+1) {log_cpm:.2f}); "
                        f"{'passes' if called else 'does not pass'} the paper's call thresholds "
                        f"(|LFC|>{LFC_THRESHOLD}, log(CPM+1)>{LOG_CPM_THRESHOLD})." + chan
                    ),
                    extras=extras,
                )
            )
        return items
