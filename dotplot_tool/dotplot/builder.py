from __future__ import annotations

from dataclasses import dataclass

from dotplot_tool.blast.pipeline import BlastHit
from dotplot_tool.genome.parser import GeneRecord


@dataclass(frozen=True)
class DotplotPoint:
    x: int
    y: int
    qseqid: str
    sseqid: str
    evalue: float
    bitscore: float
    q_locus_tag: str
    s_locus_tag: str
    q_name: str
    s_name: str


def build_dotplot_points(
    hits: list[BlastHit],
    query_index: dict[str, GeneRecord],
    subject_index: dict[str, GeneRecord],
) -> list[DotplotPoint]:
    points: list[DotplotPoint] = []

    for hit in hits:
        qrec = query_index.get(hit.qseqid)
        srec = subject_index.get(hit.sseqid)

        if qrec is None or srec is None:
            continue

        points.append(
            DotplotPoint(
                x=qrec.gene_index,
                y=srec.gene_index,
                qseqid=hit.qseqid,
                sseqid=hit.sseqid,
                evalue=hit.evalue,
                bitscore=hit.bitscore,
                q_locus_tag=qrec.locus_tag,
                s_locus_tag=srec.locus_tag,
                q_name=qrec.name,
                s_name=srec.name,
            )
        )

    return points

