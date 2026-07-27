from __future__ import annotations

from dotplot_tool.core.models import GenomeCandidate
from dotplot_tool.ncbi.client import NcbiClient


def _infer_source_from_accession(accession: str) -> str:
    accession = accession.strip().upper()
    if accession.startswith("GCF_"):
        return "RefSeq"
    if accession.startswith("GCA_"):
        return "GenBank"
    return "NA"


def _pick(d: dict, *keys: str, default: str = "") -> str:
    for k in keys:
        if k in d and isinstance(d[k], str) and d[k].strip():
            return d[k].strip()

        k2 = k.lower()
        for kk, vv in d.items():
            if (
                isinstance(kk, str)
                and kk.lower() == k2
                and isinstance(vv, str)
                and vv.strip()
            ):
                return vv.strip()

    return default


def search_assemblies(
    client: NcbiClient,
    query: str,
    retmax: int = 100,
) -> list[GenomeCandidate]:
    query = query.strip()
    if not query:
        return []

    ids = client.esearch(db="assembly", term=query, retmax=retmax)
    if not ids:
        return []

    summary = client.esummary(db="assembly", ids=ids)
    uids = summary.get("uids", [])

    results: list[GenomeCandidate] = []

    for uid in uids:
        rec = summary.get(uid)
        if not isinstance(rec, dict):
            continue

        accession = _pick(
            rec,
            "AssemblyAccession",
            "assemblyaccession",
            "assembly_accession",
            default="",
        )
        organism = _pick(
            rec,
            "Organism",
            "organism",
            "OrganismName",
            "organismname",
            default="Unknown organism",
        )
        level = _pick(
            rec,
            "AssemblyLevel",
            "assemblylevel",
            default="",
        )

        if not accession:
            continue

        source = _infer_source_from_accession(accession)

        results.append(
            GenomeCandidate(
                assembly_accession=accession,
                organism_name=organism,
                assembly_level=level,
                refseq_category=source,
            )
        )

    results.sort(
        key=lambda g: (
            0 if g.assembly_accession.startswith("GCF_") else 1,
            g.organism_name.lower(),
            g.assembly_accession,
        )
    )

    return results[:retmax]