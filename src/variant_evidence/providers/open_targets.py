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


def json_file_transport(path: str) -> Transport:
    """Offline transport: a JSON file mapping rsID -> list of normalised rows."""
    import json
    from pathlib import Path

    data = json.loads(Path(path).read_text())
    return lambda q: data.get(q.rsid, [])
