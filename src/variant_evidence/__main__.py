"""CLI: python -m variant_evidence rs2836882 --ag2 scores.csv [--ot ot.json] [--llm]"""
from __future__ import annotations

import argparse
import os

from .agent import llm_synthesizer, run_agent
from .models import VariantQuery
from .providers.ag2_precomputed import AG2PrecomputedProvider
from .providers.open_targets import OpenTargetsProvider, json_file_transport


def main() -> None:
    ap = argparse.ArgumentParser(prog="variant_evidence")
    ap.add_argument("rsid")
    ap.add_argument("--ag2", required=True, help="normalised AG2 score CSV")
    ap.add_argument("--ot", help="offline Open Targets JSON")
    ap.add_argument("--llm", action="store_true", help="use Claude for claim wording (needs ANTHROPIC_API_KEY)")
    ap.add_argument("--model", default="claude-sonnet-5-5")
    ap.add_argument("--context", nargs="+", help="restrict to cell types/tissues containing these terms")
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()

    providers = [AG2PrecomputedProvider(a.ag2)]
    if a.ot:
        providers.append(OpenTargetsProvider(json_file_transport(a.ot)))
    synth = None
    if a.llm:
        import anthropic

        if not os.environ.get("ANTHROPIC_API_KEY"):
            raise SystemExit("ANTHROPIC_API_KEY is not set")
        synth = llm_synthesizer(anthropic.Anthropic(), a.model)

    report = run_agent(VariantQuery(rsid=a.rsid), providers, synth, context=a.context)
    if a.json:
        print(report.model_dump_json(indent=2))
        return
    print(f"Variant {report.query.rsid}")
    for h in report.hypotheses[:3]:
        print(f"  {h.gene:10s} {h.direction.value:5s} conf={h.confidence.value:6s} "
              f"score={h.score:.2f} cell types: {', '.join(h.top_cell_types)}")
    print("Claims:")
    for c in report.claims:
        print(f"  - {c.text}  [{', '.join(c.evidence_ids)}]")
    if report.dropped_claims:
        print(f"({report.dropped_claims} claim(s) dropped by the verifier)")
    print("Caveats:")
    for c in report.caveats:
        print(f"  - {c}")


if __name__ == "__main__":
    main()
