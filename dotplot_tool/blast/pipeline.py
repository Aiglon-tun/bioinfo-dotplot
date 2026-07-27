from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess


class BlastPipelineError(RuntimeError):
    pass


@dataclass(frozen=True)
class BlastHit:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float


@dataclass(frozen=True)
class BlastRunResult:
    hits: list[BlastHit]
    blast_out: Path
    used_cache: bool


def ensure_blast_installed() -> None:
    missing = [exe for exe in ("makeblastdb", "blastp") if shutil.which(exe) is None]
    if missing:
        raise BlastPipelineError(
            "BLAST+ introuvable. Commandes manquantes : " + ", ".join(missing)
        )


def _run_command(cmd: list[str], cwd: str | Path | None = None) -> None:
    try:
        subprocess.run(
            cmd,
            cwd=str(cwd) if cwd is not None else None,
            check=True,
            text=True,
            capture_output=True,
        )
    except subprocess.CalledProcessError as e:
        raise BlastPipelineError(
            f"Commande échouée:\n{' '.join(cmd)}\n\nSTDOUT:\n{e.stdout}\n\nSTDERR:\n{e.stderr}"
        ) from e


def _format_param_token(value: float) -> str:
    text = f"{value:.0e}" if value != 0 and (abs(value) < 1e-3 or abs(value) >= 1e4) else f"{value:g}"
    return text.replace("+", "").replace(".", "p")


def _build_blast_output_path(
    results_dir: Path,
    query_faa: Path,
    subject_faa: Path,
    evalue_threshold: float,
    max_target_seqs: int,
) -> Path:
    evalue_token = _format_param_token(evalue_threshold)
    return results_dir / (
        f"{query_faa.stem}_vs_{subject_faa.stem}"
        f"__e{evalue_token}__mts{max_target_seqs}.tsv"
    )


def make_protein_blast_db(
    subject_faa: str | Path,
    db_dir: str | Path,
    db_name: str | None = None,
    force: bool = False,
) -> Path:
    ensure_blast_installed()

    subject_faa = Path(subject_faa)
    db_dir = Path(db_dir)

    if not subject_faa.exists():
        raise FileNotFoundError(f"FASTA sujet introuvable : {subject_faa}")

    db_dir.mkdir(parents=True, exist_ok=True)

    if db_name is None:
        db_name = subject_faa.stem

    db_prefix = db_dir / db_name

    expected = [db_prefix.with_suffix(ext) for ext in (".phr", ".pin", ".psq")]
    if not force and all(p.exists() for p in expected):
        return db_prefix

    cmd = [
        "makeblastdb",
        "-in",
        str(subject_faa),
        "-dbtype",
        "prot",
        "-out",
        str(db_prefix),
    ]
    _run_command(cmd)
    return db_prefix


def run_blastp(
    query_faa: str | Path,
    db_prefix: str | Path,
    out_path: str | Path,
    evalue: float = 1e-20,
    max_target_seqs: int = 10,
    num_threads: int = 1,
    force: bool = False,
) -> tuple[Path, bool]:
    ensure_blast_installed()

    query_faa = Path(query_faa)
    db_prefix = Path(db_prefix)
    out_path = Path(out_path)

    if not query_faa.exists():
        raise FileNotFoundError(f"FASTA requête introuvable : {query_faa}")

    db_files = [db_prefix.with_suffix(ext) for ext in (".phr", ".pin", ".psq")]
    if not all(p.exists() for p in db_files):
        raise FileNotFoundError(
            f"Base BLAST absente ou incomplète pour le préfixe : {db_prefix}"
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.exists() and not force:
        return out_path, True

    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

    cmd = [
        "blastp",
        "-query",
        str(query_faa),
        "-db",
        str(db_prefix),
        "-out",
        str(out_path),
        "-evalue",
        str(evalue),
        "-max_target_seqs",
        str(max_target_seqs),
        "-num_threads",
        str(num_threads),
        "-outfmt",
        outfmt,
    ]
    _run_command(cmd)
    return out_path, False


def parse_blast_tabular(blast_path: str | Path) -> list[BlastHit]:
    blast_path = Path(blast_path)

    if not blast_path.exists():
        raise FileNotFoundError(f"Sortie BLAST introuvable : {blast_path}")

    hits: list[BlastHit] = []

    with blast_path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) != 12:
                raise ValueError(
                    f"Ligne BLAST inattendue ({len(parts)} colonnes au lieu de 12) : {line}"
                )

            hits.append(
                BlastHit(
                    qseqid=parts[0],
                    sseqid=parts[1],
                    pident=float(parts[2]),
                    length=int(parts[3]),
                    mismatch=int(parts[4]),
                    gapopen=int(parts[5]),
                    qstart=int(parts[6]),
                    qend=int(parts[7]),
                    sstart=int(parts[8]),
                    send=int(parts[9]),
                    evalue=float(parts[10]),
                    bitscore=float(parts[11]),
                )
            )

    return hits


def filter_hits_by_evalue(
    hits: list[BlastHit],
    threshold: float = 1e-20,
) -> list[BlastHit]:
    return [hit for hit in hits if hit.evalue < threshold]


def keep_best_hit_per_query(hits: list[BlastHit]) -> list[BlastHit]:
    best_by_query: dict[str, BlastHit] = {}

    for hit in hits:
        current = best_by_query.get(hit.qseqid)
        if current is None:
            best_by_query[hit.qseqid] = hit
            continue

        if (hit.evalue < current.evalue) or (
            hit.evalue == current.evalue and hit.bitscore > current.bitscore
        ):
            best_by_query[hit.qseqid] = hit

    return list(best_by_query.values())


def run_blastp_pipeline(
    query_faa: str | Path,
    subject_faa: str | Path,
    work_dir: str | Path,
    evalue_threshold: float = 1e-20,
    max_target_seqs: int = 10,
    num_threads: int = 1,
    keep_best_only: bool = False,
    force: bool = False,
) -> BlastRunResult:
    query_faa = Path(query_faa)
    subject_faa = Path(subject_faa)
    work_dir = Path(work_dir)

    work_dir.mkdir(parents=True, exist_ok=True)

    db_dir = work_dir / "db"
    results_dir = work_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    db_prefix = make_protein_blast_db(
        subject_faa=subject_faa,
        db_dir=db_dir,
        db_name=subject_faa.stem,
        force=force,
    )

    blast_out = _build_blast_output_path(
        results_dir=results_dir,
        query_faa=query_faa,
        subject_faa=subject_faa,
        evalue_threshold=evalue_threshold,
        max_target_seqs=max_target_seqs,
    )

    blast_out, used_cache = run_blastp(
        query_faa=query_faa,
        db_prefix=db_prefix,
        out_path=blast_out,
        evalue=evalue_threshold,
        max_target_seqs=max_target_seqs,
        num_threads=num_threads,
        force=force,
    )

    hits = parse_blast_tabular(blast_out)
    hits = filter_hits_by_evalue(hits, threshold=evalue_threshold)

    if keep_best_only:
        hits = keep_best_hit_per_query(hits)

    return BlastRunResult(
        hits=hits,
        blast_out=blast_out,
        used_cache=used_cache,
    )