# Design notes

## Goal

Turn a variant ID into a short, checkable statement: likely target gene, cell types, direction of effect, confidence, and the evidence behind each part.

## Decisions

**1. The LLM never ranks.** Gene ranking and confidence come from `aggregator.py`. The LLM only words claims over records that were already retrieved. Reason: ranking has to be testable and reproducible, and an LLM ranking cannot be audited.

**2. Verify, don't trust.** `verifier.py` drops any claim that cites nothing, cites an unknown id, or contains a number absent from its cited evidence. Dropped claims are counted in the report. Reason: the typical failure of a retrieval agent in this domain is a fluent sentence with an invented number or mechanism.

**3. Abstain explicitly.** If no gene passes AG2's call thresholds the report says so. The paper states recall at high precision is only about 26-30%, so silence is the common case, not an error.

**4. Providers are swappable.** Anything that returns `EvidenceItem` works. The real AG2 file format is isolated in one reader.

## Confidence rules (heuristic, not calibrated)

- HIGH: best PIP >= 0.9, top effect |LFC| >= 0.2, and top gene's score >= 2x the runner-up
- MEDIUM: top effect >= 0.2 or PIP >= 0.9
- LOW: otherwise, and every gene other than the top one

These thresholds are my starting choices. The next step is to calibrate them against the gold set and report a reliability curve, not to defend them.

## Evaluation plan

1. Convert the real AG2 table. Re-run the harness on `paper_example` loci (expect high scores; this is not generalisation).
2. Build a larger held-out set from sources the paper did not use (for example CRISPR-validated or MPRA-validated variants) and report the same metrics there.
3. Compare three synthesizers: template, LLM, LLM without verifier. Metrics: claim drop rate, and a manual audit of 50 claims for unsupported statements.
4. Calibration: bucket HIGH/MEDIUM/LOW and report how often the top gene is correct in each bucket.

## Known gaps

- Direction-of-effect check against rare-variant LoF burden (the paper's concordance analysis) is not implemented.
- Cell-type names are free text. Mapping to Cell Ontology terms is needed before cell-type matching is reliable.
- Indels are out of scope because the source atlas has SNVs only.
- No ancestry stratification. Training and benchmark data in the paper skew European; checking this is a candidate follow-up.
