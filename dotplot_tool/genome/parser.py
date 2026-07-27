from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv


@dataclass(frozen=True)
class GeneRecord:
    protein_id: str
    gene_index: int
    genomic_accession: str
    start: int
    end: int
    strand: str
    locus_tag: str
    symbol: str
    name: str
    feature: str
    assembly: str
    assembly_unit: str
    seq_type: str
    chromosome: str
    gene_id: str
    attributes: str


@dataclass(frozen=True)
class ProteinRecord : 
    protein_id : str
    fasta_header : str
    sequence : str



def _clean_str(value: str | None) -> str:
    if value is None:
        return ""
    return value.strip()


def _parse_int(value: str | None) -> int:
    value = _clean_str(value)
    if not value:
        return 0
    return int(value)


def parse_feature_table(feature_table_path: str | Path) -> list[GeneRecord]:
    feature_table_path = Path(feature_table_path)

    if not feature_table_path.exists():
        raise FileNotFoundError(f"Feature table introuvable : {feature_table_path}")

    records_raw: list[dict] = []

    with feature_table_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Feature table sans en-tête détectable.")

        required_columns = {
            "# feature",
            "genomic_accession",
            "start",
            "end",
            "strand",
            "product_accession",
        }
        missing = required_columns.difference(reader.fieldnames)
        if missing:
            raise ValueError(
                f"Colonnes obligatoires manquantes dans {feature_table_path}: {sorted(missing)}"
            )

        for row in reader:
            feature = _clean_str(row.get("# feature")).lower()
            if feature != "cds":
                continue

            protein_id = _clean_str(row.get("product_accession"))
            genomic_accession = _clean_str(row.get("genomic_accession"))

            if not protein_id or not genomic_accession:
                continue

            start = _parse_int(row.get("start"))
            end = _parse_int(row.get("end"))
            if start <= 0 or end <= 0:
                continue

            records_raw.append(
                {
                    "protein_id": protein_id,
                    "genomic_accession": genomic_accession,
                    "start": start,
                    "end": end,
                    "strand": _clean_str(row.get("strand")),
                    "locus_tag": _clean_str(row.get("locus_tag")),
                    "symbol": _clean_str(row.get("symbol")),
                    "name": _clean_str(row.get("name")),
                    "feature": _clean_str(row.get("# feature")),
                    "assembly": _clean_str(row.get("assembly")),
                    "assembly_unit": _clean_str(row.get("assembly_unit")),
                    "seq_type": _clean_str(row.get("seq_type")),
                    "chromosome": _clean_str(row.get("chromosome")),
                    "gene_id": _clean_str(row.get("GeneID")),
                    "attributes": _clean_str(row.get("attributes")),
                }
            )

    records_raw.sort(
        key=lambda r: (
            r["genomic_accession"],
            min(r["start"], r["end"]),
            max(r["start"], r["end"]),
            r["protein_id"],
        )
    )

    records: list[GeneRecord] = []
    for gene_index, row in enumerate(records_raw):
        records.append(
            GeneRecord(
                protein_id=row["protein_id"],
                gene_index=gene_index,
                genomic_accession=row["genomic_accession"],
                start=row["start"],
                end=row["end"],
                strand=row["strand"],
                locus_tag=row["locus_tag"],
                symbol=row["symbol"],
                name=row["name"],
                feature=row["feature"],
                assembly=row["assembly"],
                assembly_unit=row["assembly_unit"],
                seq_type=row["seq_type"],
                chromosome=row["chromosome"],
                gene_id=row["gene_id"],
                attributes=row["attributes"],
            )
        )

    return records


def build_gene_index_by_protein_id(records: list[GeneRecord]) -> dict[str, GeneRecord]:
    return {record.protein_id: record for record in records}


def parse_feature_table_to_index(
    feature_table_path: str | Path,
) -> dict[str, GeneRecord]:
    records = parse_feature_table(feature_table_path)
    return build_gene_index_by_protein_id(records)


def parse_faa_ids(faa_path: str | Path) -> dict[str, ProteinRecord]:
    faa_path = Path(faa_path)

    if not faa_path.exists():
        raise FileNotFoundError(f"Fichier FASTA introuvable : {faa_path}")

    records: dict[str, ProteinRecord] = {}

    current_header: str | None = None
    current_seq_parts: list[str] = []

    def flush_record() -> None:
        nonlocal current_header, current_seq_parts #pour éviter les arguments inutiles et garder la logique simple

        if current_header is None:
            return

        header = current_header[1:] if current_header.startswith(">") else current_header
        protein_id = header.split()[0].strip()
        sequence = "".join(current_seq_parts).strip()

        if protein_id:
            records[protein_id] = ProteinRecord(
                protein_id=protein_id,
                fasta_header=header,
                sequence=sequence,
            )

    with faa_path.open("r", encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                flush_record()
                current_header = line
                current_seq_parts = []
            else:
                current_seq_parts.append(line)

    flush_record()
    return records