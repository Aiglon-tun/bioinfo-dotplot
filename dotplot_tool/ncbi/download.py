from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import gzip
import json
import shutil
import tempfile

import requests

from dotplot_tool.ncbi.client import NcbiClient


class DownloadError(RuntimeError):
    """Erreur métier liée au téléchargement et à la préparation d'un génome."""
    pass


@dataclass(frozen=True)
class DownloadedGenome:
    assembly_accession: str
    organism_name: str
    assembly_name: str
    ftp_path: str
    proteins_faa: Path
    feature_table: Path
    metadata_json: Path
    root_dir: Path


def _pick(d: dict, *keys: str, default: str = "") -> str:
    for k in keys:
        if k in d and isinstance(d[k], str) and d[k].strip():
            return d[k].strip()

        lk = k.lower()
        for kk, vv in d.items():
            if (
                isinstance(kk, str)
                and kk.lower() == lk
                and isinstance(vv, str)
                and vv.strip()
            ):
                return vv.strip()

    return default


def _normalize_accession(acc: str) -> str:
    return acc.strip().split()[0]


def _find_assembly_record(client: NcbiClient, assembly_accession: str) -> dict:
    query = _normalize_accession(assembly_accession)

    ids = client.esearch(db="assembly", term=query, retmax=10)
    if not ids:
        raise DownloadError(f"Aucun assembly trouvé pour {query}")

    summary = client.esummary(db="assembly", ids=ids)
    uids = summary.get("uids", [])

    target = query.upper()
    fallback = None

    for uid in uids:
        rec = summary.get(uid)
        if not isinstance(rec, dict):
            continue

        acc = _pick(
            rec,
            "AssemblyAccession",
            "assemblyaccession",
            "assembly_accession",
            default="",
        )
        if not acc:
            continue

        if fallback is None:
            fallback = rec

        if acc.upper() == target:
            return rec

    if fallback is not None:
        return fallback

    raise DownloadError(f"Impossible de résoudre l'assembly {query}")


def _resolve_ftp_path(rec: dict) -> str:
    ftp_path = _pick(
        rec,
        "FtpPath_RefSeq",
        "ftppath_refseq",
        "FtpPath_GenBank",
        "ftppath_genbank",
        default="",
    )

    if not ftp_path:
        raise DownloadError(
            "Aucun chemin FTP RefSeq/GenBank trouvé pour cet assembly."
        )

    return ftp_path.rstrip("/")


def _as_download_base_url(path: str) -> str:
    path = path.strip().rstrip("/")
    if path.startswith("ftp://ftp.ncbi.nlm.nih.gov/"):
        return "https://ftp.ncbi.nlm.nih.gov/" + path[len("ftp://ftp.ncbi.nlm.nih.gov/") :]
    return path


def _download_file(url: str, dest: Path, timeout: int = 120) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)

    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()

        with tempfile.NamedTemporaryFile(delete=False, dir=str(dest.parent)) as tmp:
            tmp_path = Path(tmp.name)
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    tmp.write(chunk)

    tmp_path.replace(dest)
    return dest


def _gunzip_file(src_gz: Path, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(src_gz, "rb") as fin, open(dest, "wb") as fout:
        shutil.copyfileobj(fin, fout)

    return dest


def _write_metadata(path: Path, payload: dict) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    return path


def _validate_pipeline_files(
    assembly_accession: str,
    proteins_faa: Path,
    feature_table: Path,
) -> None:
    missing: list[str] = []

    if not proteins_faa.exists():
        missing.append("protein.faa")
    if not feature_table.exists():
        missing.append("feature_table.txt")

    if missing:
        if assembly_accession.startswith("GCF_"):
            hint = "Cet assembly RefSeq existe mais ne fournit pas tous les fichiers attendus."
        elif assembly_accession.startswith("GCA_"):
            hint = (
                "Cet assembly GenBank n'est probablement pas suffisamment annoté pour ce pipeline. "
                "Essaie de préférence un assembly RefSeq (GCF)."
            )
        else:
            hint = "Essaie de préférence un assembly annoté, idéalement RefSeq (GCF)."

        raise DownloadError(
            f"Assembly non exploitable par l'application : fichiers manquants "
            f"{', '.join(missing)} pour {assembly_accession}. {hint}"
        )


def download_genome_package(
    client: NcbiClient,
    assembly_accession: str,
    out_dir: str | Path,
) -> DownloadedGenome:
    out_dir = Path(out_dir)
    rec = _find_assembly_record(client, assembly_accession)

    accession = _pick(
        rec,
        "AssemblyAccession",
        "assemblyaccession",
        "assembly_accession",
        default=assembly_accession,
    )
    organism = _pick(
        rec,
        "Organism",
        "organism",
        "OrganismName",
        "organismname",
        default="Unknown organism",
    )
    assembly_name = _pick(
        rec,
        "AssemblyName",
        "assemblyname",
        default="assembly",
    )

    ftp_path = _as_download_base_url(_resolve_ftp_path(rec))
    prefix = ftp_path.split("/")[-1]

    genome_dir = out_dir / accession
    raw_dir = genome_dir / "raw"
    meta_dir = genome_dir / "meta"

    proteins_gz = raw_dir / f"{prefix}_protein.faa.gz"
    feature_gz = raw_dir / f"{prefix}_feature_table.txt.gz"
    proteins_faa = genome_dir / "protein.faa"
    feature_table = genome_dir / "feature_table.txt"
    metadata_json = meta_dir / "metadata.json"

    protein_url = f"{ftp_path}/{prefix}_protein.faa.gz"
    feature_url = f"{ftp_path}/{prefix}_feature_table.txt.gz"

    if not proteins_faa.exists():
        if not proteins_gz.exists():
            try:
                _download_file(protein_url, proteins_gz)
            except requests.HTTPError as e:
                raise DownloadError(
                    f"Téléchargement protéines impossible pour {accession}. "
                    f"Fichier absent ou inaccessible : {protein_url}. "
                    "Cet assembly n'est peut-être pas exploitable pour le pipeline."
                ) from e
        _gunzip_file(proteins_gz, proteins_faa)

    if not feature_table.exists():
        if not feature_gz.exists():
            try:
                _download_file(feature_url, feature_gz)
            except requests.HTTPError as e:
                raise DownloadError(
                    f"Téléchargement feature table impossible pour {accession}. "
                    f"Fichier absent ou inaccessible : {feature_url}. "
                    "Cet assembly n'est peut-être pas suffisamment annoté pour le pipeline."
                ) from e
        _gunzip_file(feature_gz, feature_table)

    _validate_pipeline_files(accession, proteins_faa, feature_table)

    _write_metadata(
        metadata_json,
        {
            "assembly_accession": accession,
            "organism_name": organism,
            "assembly_name": assembly_name,
            "ftp_path": ftp_path,
            "protein_url": protein_url,
            "feature_table_url": feature_url,
        },
    )

    return DownloadedGenome(
        assembly_accession=accession,
        organism_name=organism,
        assembly_name=assembly_name,
        ftp_path=ftp_path,
        proteins_faa=proteins_faa,
        feature_table=feature_table,
        metadata_json=metadata_json,
        root_dir=genome_dir,
    )