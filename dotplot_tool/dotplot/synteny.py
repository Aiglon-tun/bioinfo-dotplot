from __future__ import annotations

from dataclasses import dataclass

from dotplot_tool.dotplot.builder import DotplotPoint


@dataclass(frozen=True)
class SyntenySegment:
    points: tuple[DotplotPoint, ...]

    @property
    def size(self) -> int:
        return len(self.points)


def detect_synteny_segments(
    points: list[DotplotPoint],
    min_points: int = 5,
    max_gene_gap: int = 8,
    diagonal_tolerance: int = 2,
) -> list[SyntenySegment]:
    if not points:
        return []

    sorted_points = sorted(points, key=lambda p: (p.x, p.y))
    used_indices: set[int] = set()
    segments: list[SyntenySegment] = []

    for i, start_point in enumerate(sorted_points):
        if i in used_indices:
            continue

        chain = [start_point]
        local_indices = {i}
        current = start_point

        while True:
            best_j = None
            best_score = None

            for j in range(i + 1, len(sorted_points)):
                if j in used_indices or j in local_indices:
                    continue

                candidate = sorted_points[j]
                dx = candidate.x - current.x
                dy = candidate.y - current.y

                if dx <= 0 or dy <= 0:
                    continue

                if dx > max_gene_gap or dy > max_gene_gap:
                    continue

                if abs(dx - dy) > diagonal_tolerance:
                    continue

                score = (abs(dx - dy), dx + dy)
                if best_score is None or score < best_score:
                    best_score = score
                    best_j = j

            if best_j is None:
                break

            current = sorted_points[best_j]
            chain.append(current)
            local_indices.add(best_j)

        if len(chain) >= min_points:
            segments.append(SyntenySegment(points=tuple(chain)))
            used_indices.update(local_indices)

    segments.sort(key=lambda s: s.size, reverse=True)
    return segments