"""Open Targets Platform credible-set evidence.

STATUS: the parser is tested against a fixture; the live GraphQL query has NOT been
run. The sandbox that built this project could not reach the API, and the exact field
names below are from memory of the Platform schema. Before relying on it, run
the Colab notebook against the live API; see the note above CREDIBLE_SET_QUERY for what is confirmed.

The transport is injected so tests and offline runs never touch the network.
"""
from __future__ import annotations

from typing import Callable

from ..models import EvidenceItem, VariantQuery

ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"

# Confirmed against the live schema: credibleSets takes studyTypes (enum incl. gwas, eqtl, pqtl, ...) and
# page {index, size}; field names of CredibleSets/CredibleSet/Loci/Locus/L2GPredictions/L2GPrediction/Study.
# Only GWAS credible sets are fetched: L2G predictions exist for those, and QTL sets are a different evidence type.
# NOT confirmed: Target.approvedSymbol, the maximum allowed page size, and that page.index is 0-based.
PAGE_SIZE = 100
MAX_PAGES = 30
CREDIBLE_SET_QUERY = """
query ($variantId: String!, $idx: Int!) {
  variant(variantId: $variantId) {
    credibleSets(studyTypes: [gwas], page: {index: $idx, size: 100}) {
      count
      rows {
        studyLocusId studyId studyType finemappingMethod
        study { traitFromSource }
        locus { count rows { posteriorProbability variant { id } } }
        l2GPredictions { rows { score target { id approvedSymbol } } }
      }
    }
  }
}
"""

Transport = Callable[[VariantQuery], list[dict]]


def parse_credible_sets(rsid: str, rows: list[dict]) -> list[EvidenceItem]:
    """rows: [{"study_id", "trait", "pip", "method", optional "study_locus_id", "l2g": [{"gene", "score"}]}].

    One credible_set item per row (pip = this variant's posterior inclusion probability) and up to
    three "l2g" items per row (Open Targets' locus-to-gene scores for that credible set).
    """
    items = []
    for r in rows:
        pip = float(r["pip"])
        key = r.get("study_locus_id") or r["study_id"]
        items.append(
            EvidenceItem(
                id=f"ot:{rsid}:{key}",
                source="open_targets",
                kind="credible_set",
                value=pip,
                detail=(
                    f"{rsid} is in a {r.get('method', 'fine-mapped')} credible set for "
                    f"'{r['trait']}' (study {r['study_id']}) with posterior inclusion probability {pip:.2f}."
                ),
                extras={"pip": pip},
            )
        )
        for g in sorted(r.get("l2g", []), key=lambda g: -g["score"])[:3]:
            score = float(g["score"])
            items.append(
                EvidenceItem(
                    id=f"ot:l2g:{rsid}:{key}:{g['gene']}",
                    source="open_targets",
                    kind="l2g",
                    gene=g["gene"],
                    value=score,
                    detail=(
                        f"Open Targets locus-to-gene model scores {g['gene']} {score:.2f} for the credible set "
                        f"of '{r['trait']}' (study {r['study_id']}) that contains {rsid}."
                    ),
                    extras={"l2g_score": score},
                )
            )
    return items


class OpenTargetsProvider:
    name = "open_targets"

    def __init__(self, transport: Transport):
        self.transport = transport

    def fetch(self, query: VariantQuery) -> list[EvidenceItem]:
        return parse_credible_sets(query.rsid, self.transport(query))


# ---- Live access (ALL UNVERIFIED: written without network access to the API) ----------------
# Workflow in Colab: run introspect() first, read the real field names, then fix the three
# GraphQL documents above/below and rows_from_response() to match. Errors are loud on purpose.

SEARCH_QUERY = """
query ($q: String!) { search(queryString: $q, entityNames: ["variant"]) { hits { id entity } } }
"""

INTROSPECT_QUERY = """
query ($name: String!) { __type(name: $name) { name fields { name type { name kind ofType { name kind } } } } }
"""


def _post(query: str, variables: dict, timeout: float = 30.0) -> dict:
    import httpx

    r = httpx.post(ENDPOINT, json={"query": query, "variables": variables}, timeout=timeout)
    r.raise_for_status()
    body = r.json()
    if body.get("errors"):
        raise RuntimeError(f"Open Targets GraphQL errors: {body['errors']}")
    return body["data"]


def introspect(type_name: str = "Variant") -> list[str]:
    """List the field names of a GraphQL type, e.g. introspect('Variant'), introspect('CredibleSet')."""
    t = _post(INTROSPECT_QUERY, {"name": type_name})["__type"]
    if t is None:
        raise RuntimeError(f"type {type_name!r} not found in the schema")
    return [f["name"] for f in t["fields"]]


_TYPE_FIELDS = """
query ($name: String!) { __type(name: $name) { fields { name type { name kind ofType { name kind
  ofType { name kind ofType { name kind ofType { name kind } } } } } } } }
"""


def _base_type(t: dict) -> str:
    while t and t.get("ofType"):
        t = t["ofType"]
    return (t or {}).get("name") or "?"


def field_types(type_name: str) -> dict[str, str]:
    """{field name: underlying type name} for a GraphQL type. Used to walk the live schema."""
    t = _post(_TYPE_FIELDS, {"name": type_name})["__type"]
    if t is None:
        raise RuntimeError(f"type {type_name!r} not found")
    return {f["name"]: _base_type(f["type"]) for f in t["fields"]}


_FIELD_ARGS = """
query ($name: String!) { __type(name: $name) { fields { name args { name type { name kind
  ofType { name kind ofType { name kind ofType { name kind } } } } } } } }
"""


def field_args(type_name: str, field: str) -> dict[str, str]:
    """{argument name: underlying type name} for one field, e.g. field_args('Variant', 'credibleSets')."""
    t = _post(_FIELD_ARGS, {"name": type_name})["__type"]
    for f in (t or {}).get("fields", []):
        if f["name"] == field:
            return {a["name"]: _base_type(a["type"]) for a in f["args"]}
    raise RuntimeError(f"{type_name}.{field} not found")


def resolve_variant_id(rsid: str) -> str:
    """Open Targets keys variants by chr_pos_ref_alt, not rsID."""
    hits = _post(SEARCH_QUERY, {"q": rsid})["search"]["hits"]
    variants = [h["id"] for h in hits if h.get("entity") == "variant"]
    if not variants:
        raise LookupError(f"no Open Targets variant found for {rsid}")
    return variants[0]


def _convert_rows(cs_rows: list[dict], variant_id: str) -> list[dict]:
    out = []
    for s in cs_rows:
        mine = [r["posteriorProbability"] for r in s["locus"]["rows"] if r["variant"]["id"] == variant_id]
        if not mine or mine[0] is None:
            continue
        out.append({
            "study_locus_id": s["studyLocusId"],
            "study_id": s["studyId"],
            "trait": (s.get("study") or {}).get("traitFromSource") or s.get("studyType") or "unknown trait",
            "pip": mine[0],
            "method": s.get("finemappingMethod") or "fine-mapped",
            "study_type": s.get("studyType"),
            "l2g": [
                {"gene": p["target"]["approvedSymbol"], "score": p["score"]}
                for p in (s.get("l2GPredictions") or {}).get("rows", [])
                if p.get("target") and p["target"].get("approvedSymbol") and p.get("score") is not None
            ],
        })
    return out


def rows_from_response(data: dict, variant_id: str) -> list[dict]:
    """Map ONE page of the live response to normalised rows. pip is THIS variant's posterior inclusion
    probability (the locus lists every variant in the credible set, so it is picked out by id)."""
    import warnings

    try:
        cs = data["variant"]["credibleSets"]
        if cs["count"] > len(cs["rows"]):
            warnings.warn(f"Open Targets returned {len(cs['rows'])} of {cs['count']} credible sets for {variant_id}")
        return _convert_rows(cs["rows"], variant_id)
    except (KeyError, TypeError) as exc:
        raise RuntimeError(
            f"Open Targets response did not match the expected shape ({exc!r}). "
            "Update CREDIBLE_SET_QUERY and rows_from_response()."
        ) from exc


def fetch_all_credible_sets(variant_id: str, post=None) -> list[dict]:
    """Page through every GWAS credible set for a variant. `post` is injectable for tests."""
    import warnings

    post = post or _post
    out: list[dict] = []
    fetched = 0
    count = None
    for idx in range(MAX_PAGES):
        try:
            cs = post(CREDIBLE_SET_QUERY, {"variantId": variant_id, "idx": idx})["variant"]["credibleSets"]
            count = cs["count"]
            page = cs["rows"]
            out += _convert_rows(page, variant_id)
        except (KeyError, TypeError) as exc:
            raise RuntimeError(f"Open Targets response did not match the expected shape ({exc!r}).") from exc
        fetched += len(page)
        if not page or fetched >= count:
            break
    if count is not None and fetched < count:
        warnings.warn(f"Open Targets: fetched {fetched} of {count} GWAS credible sets for {variant_id} (MAX_PAGES reached)")
    return out


def http_transport() -> Transport:
    """Live transport. Variant effect for one rsID costs two calls (resolve + fetch)."""

    def run(q: VariantQuery) -> list[dict]:
        vid = resolve_variant_id(q.rsid)
        return fetch_all_credible_sets(vid)

    return run


def json_file_transport(path: str) -> Transport:
    """Offline transport: a JSON file mapping rsID -> list of normalised rows."""
    import json
    from pathlib import Path

    data = json.loads(Path(path).read_text())
    return lambda q: data.get(q.rsid, [])
