"""rsID -> hg38 coordinates, needed because the AG2 download is keyed by
(chromosome, position, ref, alt) and contains no rsIDs.

Uses the Ensembl REST batch endpoint. UNVERIFIED against the live API (written without network
access): the parser follows the documented `mappings` shape; check it on the first Colab run.
`parse_ensembl_variation` is pure so it can be tested offline.
"""
from __future__ import annotations

Coord = tuple[str, int, str, str]  # (chromosome without 'chr', 1-based hg38 position, ref, alt)

ENSEMBL_URL = "https://rest.ensembl.org/variation/human"


def norm_chrom(c: str) -> str:
    c = str(c)
    return c[3:] if c.lower().startswith("chr") else c


def parse_ensembl_variation(payload: dict) -> dict[str, list[Coord]]:
    """payload: {rsid: {"mappings": [{"assembly_name", "seq_region_name", "start", "end", "allele_string"}]}}.

    SNVs only (AG2 panel is SNV-only). One coordinate per alt allele, GRCh38 mappings only.
    """
    out: dict[str, list[Coord]] = {}
    for rsid, rec in payload.items():
        coords: list[Coord] = []
        for m in (rec or {}).get("mappings", []):
            if m.get("assembly_name") != "GRCh38" or m.get("start") != m.get("end"):
                continue
            alleles = str(m.get("allele_string", "")).split("/")
            if len(alleles) < 2:
                continue
            ref, alts = alleles[0], alleles[1:]
            if len(ref) != 1:
                continue
            for alt in alts:
                if len(alt) == 1 and alt != ref:
                    coords.append((norm_chrom(m["seq_region_name"]), int(m["start"]), ref, alt))
        out[rsid] = coords
    return out


def resolve_rsids(rsids: list[str]) -> dict[str, list[Coord]]:
    import httpx

    r = httpx.post(
        ENSEMBL_URL,
        json={"ids": list(rsids)},
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        timeout=60.0,
    )
    r.raise_for_status()
    result = parse_ensembl_variation(r.json())
    missing = [x for x in rsids if not result.get(x)]
    if missing:
        raise LookupError(f"no GRCh38 SNV coordinates found for: {missing}")
    return result
