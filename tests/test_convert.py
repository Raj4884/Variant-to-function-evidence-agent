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
