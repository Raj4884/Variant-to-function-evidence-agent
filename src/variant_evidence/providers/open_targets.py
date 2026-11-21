"""Open Targets Platform credible-set evidence.

STATUS: the parser is tested against a fixture; the live GraphQL query has NOT been
run. The sandbox that built this project could not reach the API, and the exact field
names below are from memory of the Platform schema. Before relying on it, run
`scripts/check_open_targets.py` (to be written) against the live API and fix
CREDIBLE_SET_QUERY and `parse_credible_sets` to match the real schema.

The transport is injected so tests and offline runs never touch the network.
"""
from __future__ import annotations

from typing import Callable

from ..models import EvidenceItem, VariantQuery

ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"

# UNVERIFIED against the live schema. See module docstring.
CREDIBLE_SET_QUERY = """
query ($variantId: String!) {
  variant(variantId: $variantId) {
    credibleSets { rows { studyId studyLocusId finemappingMethod
      locus { rows { posteriorProbability } } study { traitFromSource } } }
  }
}
"""

Transport = Callable[[VariantQuery], list[dict]]


def parse_credible_sets(rsid: str, rows: list[dict]) -> list[EvidenceItem]:
    """rows: [{"study_id", "trait", "pip", "method"}] (normalised shape)."""
    items = []
    for r in rows:
        pip = float(r["pip"])
        items.append(
            EvidenceItem(
                id=f"ot:{rsid}:{r['study_id']}",
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


def resolve_variant_id(rsid: str) -> str:
    """Open Targets keys variants by chr_pos_ref_alt, not rsID."""
    hits = _post(SEARCH_QUERY, {"q": rsid})["search"]["hits"]
    variants = [h["id"] for h in hits if h.get("entity") == "variant"]
    if not variants:
        raise LookupError(f"no Open Targets variant found for {rsid}")
    return variants[0]


def rows_from_response(data: dict) -> list[dict]:
    """Map the (guessed) response shape to normalised rows. Fix this after introspecting."""
    try:
        sets = data["variant"]["credibleSets"]["rows"]
        out = []
        for s in sets:
            pips = [r["posteriorProbability"] for r in s["locus"]["rows"]]
            out.append({
                "study_id": s["studyId"],
                "trait": s["study"]["traitFromSource"],
                "pip": max(pips),
                "method": s.get("finemappingMethod", "fine-mapped"),
            })
        return out
    except (KeyError, TypeError) as exc:
        raise RuntimeError(
            f"Open Targets response did not match the expected shape ({exc!r}). "
            "Run introspect() and update CREDIBLE_SET_QUERY and rows_from_response()."
        ) from exc


def http_transport() -> Transport:
    """Live transport. Variant effect for one rsID costs two calls (resolve + fetch)."""

    def run(q: VariantQuery) -> list[dict]:
        vid = resolve_variant_id(q.rsid)
        return rows_from_response(_post(CREDIBLE_SET_QUERY, {"variantId": vid}))

    return run


def json_file_transport(path: str) -> Transport:
    """Offline transport: a JSON file mapping rsID -> list of normalised rows."""
    import json
    from pathlib import Path

    data = json.loads(Path(path).read_text())
    return lambda q: data.get(q.rsid, [])
