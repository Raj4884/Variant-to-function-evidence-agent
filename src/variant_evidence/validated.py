"""Held-out evaluation set from the Alsheikh et al. 2022 review of experimentally validated non-coding variants
(BMC Med Genomics 15:74, PMC8973751, Additional file 3).

SELECTION RULES (fixed before any AG2 / Open Targets output was seen for these variants):
  R1  rsID matches ^rs\\d+$
  R2  Variant Class is 'Cis Regulatory Element' or 'Promoter' (non-coding-RNA / miRNA-mechanism rows excluded:
      the target gene is then post-transcriptional, not an expression effect AG2 can predict)
  R3  Experimental_Validation includes 'Genome Editing' or 'In-vivo/Animal models'
  R4  target gene symbol present; every gene listed for the rsID after R2-R3 counts as correct
  R5  rsID not already in eval/gold.yaml
  R6  (applied in Colab) the variant is a single-alt SNV in the AG2 variant panel at the listed hg38 position
Strata by distance from the variant to the target gene's TSS (column Closest_Target_Gene_TSS_Dist):
  proximal <= 10 kb, mid 10-100 kb, distal > 100 kb. Genes at proximal distance are mostly what nearest-gene picks,
  so the distal stratum is where a method has to add something.
Known limits: the review is biased toward nearby genes; validation quality varies; the table is a compilation
(the original papers were not re-read); direction of effect is not available, so only the gene is scored.
"""
from __future__ import annotations

import pandas as pd

OK_CLASSES = {"Cis Regulatory Element", "Promoter"}
STRONG = ("Genome Editing", "In-vivo/Animal models")


def stratum(dist) -> str:
    if pd.isna(dist):
        return "unknown"
    d = abs(float(dist))
    return "proximal" if d <= 10_000 else ("mid" if d <= 100_000 else "distal")


def select_validated(df: pd.DataFrame, exclude_rsids: set[str]) -> pd.DataFrame:
    d = df.copy()
    d = d[d["Validated Variant rsID"].astype(str).str.match(r"^rs\d+$")]
    d = d[d["Variant Class"].isin(OK_CLASSES)]
    d = d[d["Experimental_Validation"].astype(str).apply(lambda v: any(s in v for s in STRONG))]
    d = d[d["Target Gene"].notna() & ~d["Target Gene"].astype(str).isin({"-", ""})]
    d = d[~d["Validated Variant rsID"].isin(exclude_rsids)]
    d["stratum"] = d["Closest_Target_Gene_TSS_Dist"].apply(stratum)
    g = d.groupby("Validated Variant rsID").agg(
        genes=("Target Gene", lambda s: sorted(set(map(str, s)))),
        chr=("Chr", "first"),
        pos_hg38=("Loc (hg38)", "first"),
        pmids=("PMID", lambda s: sorted(set(map(str, s)))),
        trait=("GWAS_DISEASE.TRAIT", "first"),
        stratum=("stratum", "first"),
        class_=("Variant Class", "first"),
    ).reset_index().rename(columns={"Validated Variant rsID": "rsid"})
    return g
