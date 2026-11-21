"""Evaluation harness.

Metrics per locus (reported separately for paper_example and heldout_candidate):
  top_gene_hit      the top-ranked gene equals the gold gene
  direction_hit     predicted direction of that gene equals gold (only scored when the gene is right)
  cell_type_hit     a gold cell type substring appears in the top gene's top cell types
  abstained         no gene passed the call thresholds (counted separately, not as a miss)
Report-level metrics:
  claim_drop_rate   fraction of candidate claims removed by the verifier

IMPORTANT: run against the bundled fixtures this only tests the plumbing. The fixtures are
synthetic. A real result needs the real AG2 table converted to the normalised schema.

usage: python eval/run_eval.py --ag2 data/fixtures/ag2_scores.csv --ot data/fixtures/ot_credible_sets.json
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from variant_evidence.agent import run_agent  # noqa: E402
from variant_evidence.models import VariantQuery  # noqa: E402
from variant_evidence.providers.ag2_precomputed import AG2PrecomputedProvider  # noqa: E402
from variant_evidence.providers.open_targets import OpenTargetsProvider, json_file_transport  # noqa: E402


def evaluate(gold_path: str, ag2_csv: str, ot_json: str | None, synthesizer=None, use_context: bool = False) -> dict:
    gold = yaml.safe_load(Path(gold_path).read_text())["loci"]
    providers = [AG2PrecomputedProvider(ag2_csv)]
    if ot_json:
        providers.append(OpenTargetsProvider(json_file_transport(ot_json)))

    tally: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    rows = []
    for g in gold:
        rep = run_agent(VariantQuery(rsid=g["rsid"]), providers, synthesizer,
                        context=g.get("context") if use_context else None)
        t = tally[g["split"]]
        t["n"] += 1
        t["claims_kept"] += len(rep.claims)
        t["claims_dropped"] += rep.dropped_claims
        if not rep.hypotheses:
            t["abstained"] += 1
            rows.append((g["rsid"], "abstain", "-", "-"))
            continue
        top = rep.hypotheses[0]
        gene_ok = top.gene == g["gene"]
        dir_known = g.get("direction") is not None
        dir_ok = dir_known and gene_ok and top.direction.value == g["direction"]
        t["direction_n"] += dir_known
        ct_ok = any(c.lower() in " ".join(top.top_cell_types).lower() for c in g["cell_types"])
        t["top_gene_hit"] += gene_ok
        t["direction_hit"] += dir_ok
        t["cell_type_hit"] += ct_ok
        rows.append((g["rsid"], top.gene, top.direction.value, top.confidence.value))
    return {"tally": {k: dict(v) for k, v in tally.items()}, "rows": rows}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gold", default=str(Path(__file__).with_name("gold.yaml")))
    ap.add_argument("--ag2", required=True)
    ap.add_argument("--ot")
    a = ap.parse_args()
    for label, ctx in (("TISSUE-AGNOSTIC", False), ("WITH TRAIT CONTEXT (loci that define one)", True)):
        res = evaluate(a.gold, a.ag2, a.ot, use_context=ctx)
        print(f"=== {label} ===")
        for r in res["rows"]:
            print("{:<12} {:<8} {:<8} {}".format(*r))
        print()
        for split, t in res["tally"].items():
            n = t["n"]
            total_claims = t["claims_kept"] + t["claims_dropped"]
            drop = t["claims_dropped"] / total_claims if total_claims else 0.0
            print(f"[{split}] n={n} top_gene={t.get('top_gene_hit', 0)}/{n} direction={t.get('direction_hit', 0)}/{t.get('direction_n', 0)} "
                  f"cell_type={t.get('cell_type_hit', 0)}/{n} abstained={t.get('abstained', 0)} claim_drop_rate={drop:.2f}")
        print()


if __name__ == "__main__":
    main()
