from __future__ import annotations
from pathlib import Path
import matplotlib.pyplot as plt 
from dotplot_tool.dotplot.builder import DotplotPoint


def plot_dotplot(
    points : list[DotplotPoint],
    n_genes_a : int,
    n_genes_b : int,
    out_path : str | Path | None = None,
    title : str | None = None,
    point_size : float = 4.0,
    point_alpha : float = 0.8,

):
    fig, ax = plt.subplots(figsize=(8,8), dpi=150)

    if points : 
        xs = [p.x for p in points]
        ys = [p.y for p in points]

        ax.scatter(
            xs,ys,
            s=point_size,
            c="black",
            alpha=point_alpha,
            linewidths=0,
            marker="s",
        )

    ax.set_xlim(-1, max(n_genes_a, 1))
    ax.set_ylim(-1, max(n_genes_b, 1))
    ax.set_xlabel("Gènes du génome A")
    ax.set_ylabel("Gènes du génome B")

    if title : 
        ax.set_title(title)

    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")

    ax.grid(False)
    ax.set_aspect("auto")

    for spine in ax.spines.values() : 
        spine.set_linewidth(1.0)

    plt.tight_layout()

    if out_path is not None :
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, bbox_inches="tight")

    return fig, ax

