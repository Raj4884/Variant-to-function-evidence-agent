from pathlib import Path

import pytest

from variant_evidence.agent import run_agent
from variant_evidence.aggregator import rank_genes
from variant_evidence.models import Claim, Confidence, Direction, EvidenceItem, VariantQuery
from variant_evidence.providers.ag2_precomputed import AG2PrecomputedProvider
from variant_evidence.providers.open_targets import OpenTargetsProvider, json_file_transport, parse_credible_sets
from variant_evidence.verifier import verify

FIX = Path(__file__).resolve().parents[1] / "data" / "fixtures"
AG2 = FIX / "ag2_scores.csv"
OT = FIX / "ot_credible_sets.json"


def providers():
    return [AG2PrecomputedProvider(AG2), OpenTargetsProvider(json_file_transport(str(OT)))]


# ---- providers ----
def test_ag2_applies_paper_call_thresholds():
    items = AG2PrecomputedProvider(AG2).fetch(VariantQuery(rsid="rs13066062"))
    by_gene = {i.gene: i for i in items}
    assert by_gene["CXCR6"].extras["_called"] == 1.0   # |lfc| 0.13 > 0.1
    assert by_gene["FYCO1"].extras["_called"] == 0.0   # |lfc| 0.06 < 0.1


def test_ag2_rejects_wrong_schema(tmp_path):
    bad = tmp_path / "bad.csv"
    bad.write_text("variant,gene\nrs1,A\n")
    with pytest.raises(ValueError, match="missing columns"):
        AG2PrecomputedProvider(bad)


def test_unknown_variant_returns_empty():
    assert AG2PrecomputedProvider(AG2).fetch(VariantQuery(rsid="rs0")) == []


def test_open_targets_parser():
    items = parse_credible_sets("rs1", [{"study_id": "S1", "trait": "t", "pip": "0.5"}])
    assert items[0].kind == "credible_set" and items[0].value == 0.5


# ---- aggregator ----
def test_top_gene_beats_decoy_and_direction_is_right():
    rep = run_agent(VariantQuery(rsid="rs12913832"), providers())
    top = rep.hypotheses[0]
    assert (top.gene, top.direction) == ("OCA2", Direction.DOWN)


def test_high_confidence_needs_pip_effect_and_dominance():
    rep = run_agent(VariantQuery(rsid="rs12913832"), providers())
    assert rep.hypotheses[0].confidence == Confidence.HIGH


def test_near_tie_is_not_high_confidence():
    # IRX3 0.12 vs IRX5 0.11: weak and not dominant
    rep = run_agent(VariantQuery(rsid="rs1421085"), providers())
    assert rep.hypotheses[0].confidence == Confidence.LOW


def test_no_call_means_abstain_with_caveat():
    rep = run_agent(VariantQuery(rsid="rs12740374"), providers())
    assert rep.hypotheses == []
    assert any("No gene passed" in c for c in rep.caveats)


def test_opposite_effects_across_conditions_are_flagged():
    items = [
        EvidenceItem(id="a", source="s", kind="variant_effect", gene="G", cell_type="epithelial",
                     value=-0.3, extras={"_called": 1.0}),
        EvidenceItem(id="b", source="s", kind="variant_effect", gene="G", cell_type="endothelial",
                     value=0.2, extras={"_called": 1.0}),
    ]
    h = rank_genes(items)[0]
    assert any("opposite" in n.lower() for n in h.notes)


# ---- verifier ----
def _items():
    return AG2PrecomputedProvider(AG2).fetch(VariantQuery(rsid="rs2836882"))


def test_verifier_keeps_grounded_claim():
    it = _items()[0]
    kept, dropped = verify([Claim(text=f"Predicted LFC is {it.value}", evidence_ids=[it.id])], _items())
    assert len(kept) == 1 and not dropped


def test_verifier_drops_uncited_unknown_and_wrong_number():
    it = _items()[0]
    claims = [
        Claim(text="ETS2 is the causal gene", evidence_ids=[]),
        Claim(text="Effect is large", evidence_ids=["ag2:nope"]),
        Claim(text="Predicted LFC is 0.99", evidence_ids=[it.id]),
    ]
    kept, dropped = verify(claims, _items())
    assert kept == [] and len(dropped) == 3


def test_verifier_number_matching_is_rounding_level_not_fuzzy():
    it = _items()[0]  # value 0.31, log_cpm 2.4
    ok = lambda text: bool(verify([Claim(text=text, evidence_ids=[it.id])], _items())[0])
    assert ok("Predicted LFC is 0.31")
    assert ok("Predicted LFC is 0.3")      # same value at 1 decimal
    assert ok("Predicted LFC is +0.31 and log(CPM+1) is 2.40")
    assert not ok("Predicted LFC is 0.99")
    assert not ok("Predicted LFC is 0.32")  # near but wrong
    assert not ok("Predicted LFC is 1")     # the old bug: '+1' in log(CPM+1) must not whitelist 1


def test_verifier_ignores_digits_inside_identifiers():
    it = _items()[0]
    kept, _ = verify([Claim(text="rs2836882 affects ETS2 expression", evidence_ids=[it.id])], _items())
    assert len(kept) == 1


# ---- end to end with a misbehaving LLM ----
def test_hallucinating_synthesizer_is_contained():
    def bad(hyps, items):
        good = items[0]
        return [
            Claim(text=good.detail, evidence_ids=[good.id]),
            Claim(text="Predicted LFC is 7.5 in neurons", evidence_ids=[good.id]),
            Claim(text="This variant causes disease", evidence_ids=[]),
        ]

    rep = run_agent(VariantQuery(rsid="rs2836882"), providers(), bad)
    assert len(rep.claims) == 1 and rep.dropped_claims == 2
