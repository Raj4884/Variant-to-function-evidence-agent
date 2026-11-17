# variant-evidence

Given a non-coding variant (rsID), this tool reports which gene it most likely regulates, in which cell type, in which direction, and how confident that call is. Every sentence in the report points at a source record, and a verifier removes any sentence it cannot trace.

It sits on top of the precomputed predictions from AlphaGenome 2 (Gayoso et al., Google DeepMind, 2026) and Open Targets fine-mapping. It does not retrain or reimplement the model.

## Status

Early. What exists and is tested:

- AG2 score reader with the paper's call thresholds (`|LFC| > 0.1`, `log(CPM+1) > 1.5`)
- Open Targets credible-set parser (offline fixture only)
- Deterministic gene ranking and confidence rules
- Claim verifier (citation check plus number check)
- Evaluation harness with a gold set split into paper examples and held-out candidates

What does not exist yet:

- The real AG2 download has not been converted to the normalised schema. The format is unverified.
- The live Open Targets query has not been run. Field names are unverified.
- Everything in `data/fixtures/` is synthetic. Running the harness on it tests plumbing only and says nothing about accuracy.

## Run

```
pip install -e ".[dev,llm]"
pytest
python -m variant_evidence rs12913832 --ag2 data/fixtures/ag2_scores.csv --ot data/fixtures/ot_credible_sets.json
python eval/run_eval.py --ag2 data/fixtures/ag2_scores.csv --ot data/fixtures/ot_credible_sets.json
```

Add `--llm` (with `ANTHROPIC_API_KEY` set) to let Claude word the claims. The ranking never depends on the LLM.

## Design in one paragraph

Providers return `EvidenceItem` records. The aggregator ranks genes from those records with fixed rules. A synthesizer (template or LLM) writes claims over the ranked evidence. The verifier keeps a claim only if it cites existing evidence ids and every number in it appears in the cited records. See `docs/design.md`.

## Limits inherited from the model

AG2 recovers a minority of causal variants at high precision, is weak outside accessible chromatin, was trained on healthy steady-state samples, and the paper's atlas covers SNVs only. A missing call is not evidence of no effect. This tool is for research, not clinical use.
