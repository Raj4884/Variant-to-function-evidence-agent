import math

import pytest

from variant_evidence.models import EvidenceItem
from variant_evidence.convert import normalise_rows, write_normalised_csv
from variant_evidence.models import VariantQuery
from variant_evidence.providers.ag2_precomputed import AG2PrecomputedProvider

RENAME = {"variant": "rsid", "gene": "symbol", "condition": "cell", "lfc": "effect", "log_cpm": "expr"}


def raw():
    return [
        {"rsid": "rs1", "symbol": "A", "cell": "T cell", "effect": "0.5", "expr": "2.0"},
        {"rsid": "rs2", "symbol": "B", "cell": "NK", "effect": "-0.2", "expr": "3.0"},
        {"rsid": "rs1", "symbol": "C", "cell": "T cell", "effect": "nan", "expr": "2.0"},
        {"rsid": "rs1", "symbol": "D", "cell": "T cell", "effect": "", "expr": "2.0"},
    ]


def test_filters_variants_and_drops_bad_numbers():
    out = list(normalise_rows(raw(), RENAME, rsids={"rs1"}))
    assert [r["gene"] for r in out] == ["A"]


def test_lfc_scale_converts_log2_to_natural_log():
    out = list(normalise_rows(raw(), RENAME, lfc_scale=math.log(2)))
    assert out[0]["lfc"] == pytest.approx(0.5 * math.log(2))


def test_missing_raw_column_gives_helpful_error():
    bad = dict(RENAME, lfc="nope")
    with pytest.raises(ValueError, match="available"):
        list(normalise_rows(raw(), bad))


def test_roundtrip_through_provider(tmp_path):
    p = tmp_path / "n.csv"
    n = write_normalised_csv(normalise_rows(raw(), RENAME), p)
    assert n == 2
    items = AG2PrecomputedProvider(p).fetch(VariantQuery(rsid="rs1"))
    assert items[0].gene == "A" and items[0].extras["_called"] == 1.0


def test_ag2_expression_rows_match_on_coordinates():
    from variant_evidence.convert import normalise_ag2_expression

    c2r = {("6", 100, "A", "G"): "rs9"}
    base = dict(chromosome="chr6", position=100, ref="A", alt="G", gene_name="G1",
                biosample_name="T cell", channel_name="chA",
                sc_gene_expression_score=0.4, sc_gene_expression_active_score=2.1)
    rows = [
        base,
        dict(base, position=101),                                   # different variant: skipped
        dict(base, channel_name="denoised", sc_gene_expression_active_score=float("nan")),  # no active score
    ]
    out = list(normalise_ag2_expression(rows, c2r))
    assert out == [{"variant": "rs9", "gene": "G1", "condition": "T cell", "channel": "chA", "lfc": 0.4, "log_cpm": 2.1}]


def test_parse_ensembl_snv_grch38_only():
    from variant_evidence.resolve import parse_ensembl_variation

    payload = {"rs9": {"mappings": [
        {"assembly_name": "GRCh38", "seq_region_name": "6", "start": 100, "end": 100, "allele_string": "A/G/T"},
        {"assembly_name": "GRCh37", "seq_region_name": "6", "start": 50, "end": 50, "allele_string": "A/G"},
        {"assembly_name": "GRCh38", "seq_region_name": "6", "start": 7, "end": 8, "allele_string": "AT/A"},
    ]}}
    assert parse_ensembl_variation(payload)["rs9"] == [("6", 100, "A", "G"), ("6", 100, "A", "T")]


def test_collapse_channels_uses_median_and_counts_agreement(tmp_path):
    from variant_evidence.convert import collapse_channels

    rows = [dict(variant="rs9", gene="G", condition="T cell", channel=c, lfc=l, log_cpm=2.0)
            for c, l in [("a", 0.30), ("b", 0.20), ("c", -0.10), ("d", 0.25)]]
    (out,) = collapse_channels(rows)
    assert out["n_channels"] == 4 and out["n_agree"] == 3
    assert out["lfc"] == 0.225


def test_provider_reports_channels_and_panel_status(tmp_path):
    p = tmp_path / "n.csv"
    rows = [dict(variant="rs1", gene="G", condition="T cell", channel=c, lfc=0.4, log_cpm=2.5) for c in "abc"]
    from variant_evidence.convert import collapse_channels
    write_normalised_csv(collapse_channels(rows), p)
    prov = AG2PrecomputedProvider(p, scored_rsids={"rs1", "rs2"})
    (item,) = prov.fetch(VariantQuery(rsid="rs1"))
    assert item.extras["n_channels"] == 3.0 and "3 agree" in item.detail
    assert prov.status("rs1") == "has_rows"
    assert prov.status("rs2") == "scored_no_effect"
    assert prov.status("rs3") == "not_in_panel"


def test_multigene_locus_is_not_singled_out_and_context_filters(tmp_path):
    from variant_evidence.agent import run_agent
    from variant_evidence.models import Confidence

    rows = []
    for g, l in [("CELSR2", 1.0), ("PSRC1", 0.9), ("SORT1", 0.4)]:
        rows.append(dict(variant="rs7", gene=g, condition="liver hepatocyte", lfc=l, log_cpm=3.0))
    rows.append(dict(variant="rs7", gene="OTHER", condition="lung cell", lfc=1.5, log_cpm=3.0))
    p = tmp_path / "n.csv"
    write_normalised_csv(rows, p)
    prov = AG2PrecomputedProvider(p)
    rep = run_agent(VariantQuery(rsid="rs7"), [prov])
    assert rep.hypotheses[0].gene == "OTHER"                     # tissue-agnostic ranking
    rep = run_agent(VariantQuery(rsid="rs7"), [prov], context=["hepatocyte"])
    top = rep.hypotheses[0]
    assert top.gene == "CELSR2" and top.confidence == Confidence.LOW
    assert any("Multi-gene locus" in n and "PSRC1" in n for n in top.notes)
    assert any("restricted" in c for c in rep.caveats)


def test_open_targets_live_response_uses_this_variants_pip_and_l2g():
    from variant_evidence.baselines import l2g_top_gene
    from variant_evidence.providers.open_targets import parse_credible_sets, rows_from_response

    vid = "1_100_A_G"
    data = {"variant": {"credibleSets": {"count": 1, "rows": [{
        "studyLocusId": "L1", "studyId": "S1", "studyType": "gwas", "finemappingMethod": "SuSiE",
        "study": {"traitFromSource": "LDL"},
        "locus": {"count": 2, "rows": [
            {"posteriorProbability": 0.99, "variant": {"id": "1_999_C_T"}},   # another variant in the set
            {"posteriorProbability": 0.40, "variant": {"id": vid}},
        ]},
        "l2GPredictions": {"rows": [{"score": 0.8, "target": {"id": "E1", "approvedSymbol": "SORT1"}},
                                    {"score": 0.1, "target": {"id": "E2", "approvedSymbol": "PSRC1"}}]},
    }]}}}
    rows = rows_from_response(data, vid)
    assert rows[0]["pip"] == 0.40                     # not the set maximum
    items = parse_credible_sets("rs7", rows)
    assert l2g_top_gene(items) == ("SORT1", 0.8)
    assert {i.kind for i in items} == {"credible_set", "l2g"}


def test_l2g_agreement_note_is_added():
    from variant_evidence.aggregator import rank_genes
    from variant_evidence.providers.open_targets import parse_credible_sets

    items = parse_credible_sets("rs7", [{"study_id": "S", "trait": "t", "pip": 0.5,
                                         "l2g": [{"gene": "PSRC1", "score": 0.7}]}])
    items.append(EvidenceItem(id="a", source="ag2_precomputed", kind="variant_effect", gene="SORT1",
                              cell_type="liver", value=0.4, extras={"_called": 1.0, "log_cpm": 3.0}))
    (h,) = rank_genes(items)
    assert any("prefers PSRC1" in n for n in h.notes)
