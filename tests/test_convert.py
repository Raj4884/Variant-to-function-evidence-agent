import math

import pytest

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
    assert out == [{"variant": "rs9", "gene": "G1", "condition": "T cell (chA)", "lfc": 0.4, "log_cpm": 2.1}]


def test_parse_ensembl_snv_grch38_only():
    from variant_evidence.resolve import parse_ensembl_variation

    payload = {"rs9": {"mappings": [
        {"assembly_name": "GRCh38", "seq_region_name": "6", "start": 100, "end": 100, "allele_string": "A/G/T"},
        {"assembly_name": "GRCh37", "seq_region_name": "6", "start": 50, "end": 50, "allele_string": "A/G"},
        {"assembly_name": "GRCh38", "seq_region_name": "6", "start": 7, "end": 8, "allele_string": "AT/A"},
    ]}}
    assert parse_ensembl_variation(payload)["rs9"] == [("6", 100, "A", "G"), ("6", 100, "A", "T")]
