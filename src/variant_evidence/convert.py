"""Convert a raw AlphaGenome 2 score table into the normalised CSV the reader expects.

The raw download format is unverified (see providers/ag2_precomputed.py), so the column
mapping is an argument, not a constant. Typical use from the Colab notebook:

    rows = con.execute("SELECT * FROM read_parquet(...) WHERE <variant col> IN (...)").fetchall()
    cols = [d[0] for d in con.description]
    normalise_rows((dict(zip(cols, r)) for r in rows), rename={...}, rsids={...})

The caller does the heavy filtering (DuckDB) so this module never needs the whole atlas in memory.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable, Iterator

NORMALISED = ["variant", "gene", "condition", "lfc", "log_cpm"]


def normalise_rows(
    rows: Iterable[dict],
    rename: dict[str, str],
    rsids: set[str] | None = None,
    lfc_scale: float = 1.0,
) -> Iterator[dict]:
    """rename maps normalised name -> raw column name.

    lfc_scale multiplies the raw effect. The paper's expression scores are natural-log fold
    changes. If the raw file is log2, pass lfc_scale=math.log(2). CHECK the units in the real file.
    """
    missing_keys = set(NORMALISED) - set(rename)
    if missing_keys:
        raise ValueError(f"rename is missing keys: {sorted(missing_keys)}")
    checked = False
    for row in rows:
        if not checked:
            absent = [raw for raw in rename.values() if raw not in row]
            if absent:
                raise ValueError(f"raw columns not found: {absent}; available: {sorted(row)}")
            checked = True
        variant = str(row[rename["variant"]])
        if rsids is not None and variant not in rsids:
            continue
        try:
            lfc = float(row[rename["lfc"]]) * lfc_scale
            log_cpm = float(row[rename["log_cpm"]])
        except (TypeError, ValueError):
            continue
        if math.isnan(lfc) or math.isnan(log_cpm):
            continue
        yield {
            "variant": variant,
            "gene": str(row[rename["gene"]]),
            "condition": str(row[rename["condition"]]),
            "lfc": lfc,
            "log_cpm": log_cpm,
        }


def write_normalised_csv(rows: Iterable[dict], out_path: str | Path) -> int:
    n = 0
    with Path(out_path).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=NORMALISED)
        w.writeheader()
        for r in rows:
            w.writerow(r)
            n += 1
    return n
