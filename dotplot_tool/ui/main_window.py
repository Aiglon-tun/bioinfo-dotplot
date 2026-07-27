from __future__ import annotations

from pathlib import Path


from PySide6.QtCore import Qt
from PySide6.QtCore import QThreadPool
from PySide6.QtWidgets import (
    QAbstractItemView,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QFileDialog,
    QMenuBar,
    QPushButton,
    QPlainTextEdit,
    QTextEdit,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from dotplot_tool.core.models import GenomeCandidate
from dotplot_tool.ncbi.client import NcbiClient
from dotplot_tool.ncbi.search import search_assemblies
from dotplot_tool.ncbi.download import download_genome_package
from dotplot_tool.utils.threads import Worker
from dotplot_tool.genome.parser import parse_feature_table, parse_feature_table_to_index
from dotplot_tool.blast.pipeline import run_blastp_pipeline
from dotplot_tool.dotplot.builder import build_dotplot_points
from dotplot_tool.dotplot.synteny import detect_synteny_segments

import csv
from datetime import datetime

# Fonctions utilitaires pour les tâches en arrière-plan

SEARCH_LIMITS = 100  #nombre maximum de résultats à récupérer lors de la recherche NCBI, pour éviter de surcharger l'interface utilisateur et les ressources du système avec trop de données.


#dowload_two_genomes() sert à télécharger les données de deux génomes sélectionnés par l'utilisateur, en utilisant le client NCBI et en stockant les fichiers dans un répertoire de cache local. 
#Elle retourne les packages de données téléchargés pour les deux génomes, qui contiennent les chemins vers les fichiers nécessaires pour la suite du pipeline (FAA et feature table).
def download_two_genomes(
    client: NcbiClient,
    genome_a: GenomeCandidate,
    genome_b: GenomeCandidate,
    cache_dir: str | Path,
):
    cache_dir = Path(cache_dir)
    pkg_a = download_genome_package(client, genome_a.assembly_accession, cache_dir)
    pkg_b = download_genome_package(client, genome_b.assembly_accession, cache_dir)
    return pkg_a, pkg_b


#compute_dotplot_data() est une fonction qui exécute l'ensemble du pipeline de calcul du dotplot pour deux génomes donnés.
#Elle prend en entrée les chemins vers les fichiers FAA et feature table des deux génomes, ainsi que les paramètres pour le BLASTP, et retourne un dictionnaire contenant les points 
#du dotplot, le nombre de gènes dans chaque génome, et le nombre de hits retenus après le filtrage BLASTP.

def compute_dotplot_data(
    faa_a: str,
    faa_b: str,
    feature_a: str,
    feature_b: str,
    work_dir: str,
    evalue_threshold: float = 1e-20,
    max_target_seqs: int = 10,
    num_threads: int = 2,
    keep_best_only: bool = True,
):
    records_a = parse_feature_table(feature_a)
    records_b = parse_feature_table(feature_b)

    index_a = parse_feature_table_to_index(feature_a)
    index_b = parse_feature_table_to_index(feature_b)

    blast_result = run_blastp_pipeline(
        query_faa=faa_a,
        subject_faa=faa_b,
        work_dir=work_dir,
        evalue_threshold=evalue_threshold,
        max_target_seqs=max_target_seqs,
        num_threads=num_threads,
        keep_best_only=keep_best_only,
    )

    points = build_dotplot_points(blast_result.hits, index_a, index_b)
    
    synteny_segments = detect_synteny_segments(
        points,
        min_points = 5,
        max_gene_gap = 8,
        diagonal_tolerance = 2,
    )

    return {
        "points": points,
        "n_genes_a": len(records_a),
        "n_genes_b": len(records_b),
        "n_hits": len(blast_result.hits),
        "blast_cache_used": blast_result.used_cache,
        "blast_result_path": str(blast_result.blast_out),
        "synteny_segments" : synteny_segments,
        "n_synteny_segments" : len(synteny_segments),
    }

class DotplotWindow(QMainWindow):
    def __init__(self, dotplot_data: dict, title: str) -> None:
        super().__init__()

        self.setWindowTitle(title)
        self.resize(1400, 900)

        self._title = title
        self._dotplot_points = list(dotplot_data["points"])
        self._n_genes_a = int(dotplot_data["n_genes_a"])
        self._n_genes_b = int(dotplot_data["n_genes_b"])
        self._n_hits = int(dotplot_data["n_hits"])
        self.synteny_segments = dotplot_data.get("synteny_segments", [])

        self._dotplot_pixel_positions = []
        self._dotplot_ax = None
        self._dotplot_annotation = None
        self._last_hover_index = None
        self._selected_point_index = None
        self._selected_marker = None
        self._show_synteny = False

        root = QWidget()
        self.setCentralWidget(root)
        layout = QVBoxLayout(root)

        self.hover_label = QLabel("Survol : aucune paire de gènes.")
        self.hover_label.setWordWrap(True)
        self.hover_label.setFixedHeight(80)
        layout.addWidget(self.hover_label)

        self.selected_label = QLabel("Sélection : aucune paire de gènes figée.")
        layout.addWidget(self.selected_label)

        self.details_panel = QPlainTextEdit()
        self.details_panel.setReadOnly(True)
        self.details_panel.setFixedHeight(100)
        self.details_panel.setPlainText(
            "Clique sur un point du dotplot pour figer la paire de gènes."
        )
        layout.addWidget(self.details_panel)

        self.figure = Figure(figsize=(8, 8), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        self.btn_toggle_synteny = QPushButton("Synténies")
        self.btn_toggle_synteny.setCheckable(True)
        self.btn_toggle_synteny.setFixedWidth(90)
        self.btn_toggle_synteny.setFixedHeight(26)
        self.btn_toggle_synteny.setToolTip("Afficher / Masquer les régions de synténi en rouge")
        self.btn_toggle_synteny.setStyleSheet("""
        QPushButton {
        background-color: #f4f4f4;
        color: #333333;
        border: 1px solid #c8c8c8;
        border-radius: 5px;
        padding: 2px 8px;
        font-size: 11px;
    }
        QPushButton:hover:!disabled { background-color: #ececec; }
        QPushButton:checked {
        background-color: #f3d6d6;
        color: #7a1f1f;
        border: 1px solid #cc9d9d;
    }
    """)
        self.btn_toggle_synteny.clicked.connect(self._toggle_synteny_display)
        toolbar_row = QHBoxLayout()
        toolbar_row.addWidget(self.toolbar)
        toolbar_row.addWidget(self.btn_toggle_synteny)
        self.btn_toggle_synteny = QPushButton("Synténies")
        self.btn_toggle_synteny.setCheckable(True)
        self.btn_toggle_synteny.setFixedWidth(90)
        self.btn_toggle_synteny.setFixedHeight(26)
        self.btn_toggle_synteny.setToolTip("Afficher / masquer les régions de synténie en rouge")
        self.btn_toggle_synteny.setStyleSheet("""
        QPushButton {
        background-color: #f4f4f4;
        color: #333333;
        border: 1px solid #c8c8c8;
        border-radius: 5px;
        padding: 2px 8px;
        font-size: 11px;
        }
        QPushButton:hover:!disabled { background-color: #ececec; }
        QPushButton:checked {
        background-color: #f3d6d6;
        color: #7a1f1f;
        border: 1px solid #cc9d9d;
        }
        """)
        self.btn_toggle_synteny.clicked.connect(self._toggle_synteny_display)

        toolbar_row = QHBoxLayout()
        toolbar_row.addWidget(self.toolbar)
        toolbar_row.addWidget(self.btn_toggle_synteny)
        layout.addLayout(toolbar_row)
        layout.addWidget(self.canvas)


    

        self.canvas.mpl_connect("motion_notify_event", self._on_plot_hover)
        self.canvas.mpl_connect("figure_leave_event", self._on_plot_leave)
        self.canvas.mpl_connect("draw_event", self._on_canvas_draw)
        self.canvas.mpl_connect("button_press_event", self._on_plot_click)

        self._draw_plot()
        
    def _toggle_synteny_display(self) -> None :
        self._show_synteny = self.btn_toggle_synteny.isChecked()
        self._draw_plot()

    def _draw_plot(self) -> None:
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        self._dotplot_ax = ax

        if self._dotplot_points:
            xs = [p.x for p in self._dotplot_points]
            ys = [p.y for p in self._dotplot_points]
            ax.scatter(
                xs,
                ys,
                s=4,
                c="black",
                alpha=0.8,
                linewidths=0,
                marker="s",
            )
        if self._show_synteny :
            for segment in self.synteny_segments :
                seg_x = [p.x for p in segment.points]
                seg_y = [p.y for p in segment.points]
                ax.plot(seg_x, seg_y, color="red", linewidth=0.9, alpha=0.95, zorder=3)
                ax.scatter(seg_x, seg_y, s=10, c="red", alpha=0.95, linewidths=0, zorder=4)

        ax.set_xlim(-1, max(self._n_genes_a, 1))
        ax.set_ylim(-1, max(self._n_genes_b, 1))
        ax.set_xlabel("Gènes du génome A")
        ax.set_ylabel("Gènes du génome B")
        ax.set_title(f"{self._title} — Hits retenus : {self._n_hits}")
        ax.grid(False)

        self._dotplot_annotation = ax.annotate(
            "",
            xy=(0, 0),
            xytext=(12, 12),
            textcoords="offset points",
            bbox=dict(boxstyle="round", fc="white", ec="black", alpha=0.95),
            arrowprops=dict(arrowstyle="->", color="black"),
        )
        self._dotplot_annotation.set_visible(False)

        self.figure.tight_layout()
        self.canvas.draw()
        self._update_hover_cache()

    def _format_hover_text(self, point) -> str:
        q_locus = point.q_locus_tag or "-"
        s_locus = point.s_locus_tag or "-"
        q_name = point.q_name or "-"
        s_name = point.s_name or "-"

        return (
            f"A : qseqid={point.qseqid} | locus_tag={q_locus} | protéine={q_name}\n"
            f"B : sseqid={point.sseqid} | locus_tag={s_locus} | protéine={s_name}\n"
            f"BLAST : E-value={point.evalue:.2e} | bitscore={point.bitscore:.1f}"
        )

    def _format_selected_text(self, point) -> str:
        q_locus = point.q_locus_tag or "-"
        s_locus = point.s_locus_tag or "-"
        q_name = point.q_name or "-"
        s_name = point.s_name or "-"

        return (
            "Paire de gènes sélectionnée\n"
            "--------------------------\n"
            f"Génome A\n"
            f"  - qseqid      : {point.qseqid}\n"
            f"  - locus_tag   : {q_locus}\n"
            f"  - protéine    : {q_name}\n"
            f"\n"
            f"Génome B\n"
            f"  - sseqid      : {point.sseqid}\n"
            f"  - locus_tag   : {s_locus}\n"
            f"  - protéine    : {s_name}\n"
            f"\n"
            f"BLAST\n"
            f"  - E-value     : {point.evalue:.2e}\n"
            f"  - bitscore    : {point.bitscore:.1f}\n"
            f"\n"
            f"Coordonnées dotplot\n"
            f"  - x           : {point.x}\n"
            f"  - y           : {point.y}"
        )

    def _update_hover_cache(self) -> None:
        self._dotplot_pixel_positions = []

        if self._dotplot_ax is None or not self._dotplot_points:
            return

        transformed = self._dotplot_ax.transData.transform(
            [(p.x, p.y) for p in self._dotplot_points]
        )
        self._dotplot_pixel_positions = [
            (float(px), float(py)) for px, py in transformed
        ]

    def _on_canvas_draw(self, event) -> None:
        self._update_hover_cache()

    def _find_hovered_point_index(self, event, max_pixel_distance: float = 8.0):
        if self._dotplot_ax is None:
            return None
        if event.inaxes is not self._dotplot_ax:
            return None
        if event.x is None or event.y is None:
            return None
        if not self._dotplot_points or not self._dotplot_pixel_positions:
            return None

        mouse_x = float(event.x)
        mouse_y = float(event.y)

        best_index = None
        best_d2 = max_pixel_distance * max_pixel_distance

        for i, (px, py) in enumerate(self._dotplot_pixel_positions):
            d2 = (px - mouse_x) ** 2 + (py - mouse_y) ** 2
            if d2 <= best_d2:
                best_d2 = d2
                best_index = i

        return best_index

    def _clear_hover_feedback(self) -> None:
        self._last_hover_index = None
        self.hover_label.setText("Survol : aucune paire de gènes.")

        if self._dotplot_annotation is not None:
            self._dotplot_annotation.set_visible(False)

        self.canvas.draw_idle()

    def _on_plot_hover(self, event) -> None:
        point_index = self._find_hovered_point_index(event)

        if point_index is None:
            self._clear_hover_feedback()
            return

        if point_index == self._last_hover_index:
            return

        self._last_hover_index = point_index
        point = self._dotplot_points[point_index]

        self.hover_label.setText(self._format_hover_text(point))

        if self._dotplot_annotation is not None:
            q_short = point.q_locus_tag or point.qseqid
            s_short = point.s_locus_tag or point.sseqid
            self._dotplot_annotation.xy = (point.x, point.y)
            self._dotplot_annotation.set_text(f"{q_short} ↔ {s_short}")
            self._dotplot_annotation.set_visible(True)

        self.canvas.draw_idle()

    def _on_plot_leave(self, event) -> None:
        self._clear_hover_feedback()

    def _draw_selected_marker(self, point) -> None:
        if self._dotplot_ax is None:
            return

        if self._selected_marker is not None:
            try:
                self._selected_marker.remove()
            except ValueError:
                pass
            self._selected_marker = None

        self._selected_marker = self._dotplot_ax.scatter(
            [point.x],
            [point.y],
            s=90,
            facecolors="none",
            edgecolors="red",
            linewidths=1.8,
            zorder=5,
        )

    def _show_selected_point(self, point_index: int) -> None:
        if point_index < 0 or point_index >= len(self._dotplot_points):
            return

        self._selected_point_index = point_index
        point = self._dotplot_points[point_index]

        short_q = point.q_locus_tag or point.qseqid
        short_s = point.s_locus_tag or point.sseqid

        self.selected_label.setText(f"Sélection : {short_q} ↔ {short_s}")
        self.details_panel.setPlainText(self._format_selected_text(point))
        self._draw_selected_marker(point)
        self.canvas.draw_idle()

    def _clear_selected_point(self) -> None:
        self._selected_point_index = None
        self.selected_label.setText("Sélection : aucune paire de gènes figée.")
        self.details_panel.setPlainText(
            "Clique sur un point du dotplot pour figer la paire de gènes."
        )

        if self._selected_marker is not None:
            try:
                self._selected_marker.remove()
            except ValueError:
                pass
            self._selected_marker = None

        self.canvas.draw_idle()

    def _on_plot_click(self, event) -> None:
        if self._dotplot_ax is None:
            return
        if event.inaxes is not self._dotplot_ax:
            return

        if event.button == 3:
            self._clear_selected_point()
            return

        if event.button != 1:
            return

        point_index = self._find_hovered_point_index(event, max_pixel_distance=10.0)
        if point_index is None:
            return

        self._show_selected_point(point_index)



class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        # ── Menu Fichier (haut gauche, barre native Qt) ─────────────
        
        self.setWindowTitle("DotPlot - Recherche NCBI - Malek Louiz 2026 UM4BM748")
        self.resize(1100, 700)
        menu_bar = self.menuBar()
        file_menu = menu_bar.addMenu("Fichier")

        self._action_export_hits = file_menu.addAction("Exporter les hits BLASTP (.csv)")
        self._action_export_hits.setEnabled(False)
        self._action_export_hits.triggered.connect(self._export_hits_csv)

        self._action_export_synteny = file_menu.addAction("Exporter les segments de synténie (.csv)")
        self._action_export_synteny.setEnabled(False)
        self._action_export_synteny.triggered.connect(self._export_synteny_csv)

        file_menu.addSeparator()

        self._action_export_summary = file_menu.addAction("Sauvegarder le résumé de l'analyse (.txt)")
        self._action_export_summary.setEnabled(False)
        self._action_export_summary.triggered.connect(self._export_summary_txt)
        # ────────────────────────────────────────────────────────────
        self.selected_a: GenomeCandidate | None = None #permet de stocker le génome sélectionné comme "génome A" par l'utilisateur, qui sera utilisé dans la suite du pipeline pour le téléchargement et le calcul du dotplot.
        self.selected_b: GenomeCandidate | None = None

        self.downloaded_a = None #permet de stocker les données téléchargées pour le génome A après que l'utilisateur ait lancé le téléchargement, ce qui inclut les chemins vers les fichiers FAA et feature table nécessaires pour le calcul du dotplot.
        self.downloaded_b = None

        self.cache_dir = Path("cache/genome")
        self.threadpool = QThreadPool.globalInstance() #permet de gérer les tâches en arrière-plan (comme les recherches NCBI, les téléchargements, et le calcul du dotplot) sans bloquer l'interface utilisateur, en utilisant un pool de threads global.
        self.ncbi = NcbiClient(email="louiz_m@yahoo.fr") #ici, on crée une instance du client NCBI qui sera utilisée pour effectuer les requêtes de recherche et de téléchargement des génomes. L'email est requis par NCBI pour identifier l'utilisateur lors des requêtes.

        self._search_worker = None #permet de stocker la référence au worker en cours d'exécution pour la recherche NCBI, afin de pouvoir gérer son état et ses résultats.
        self._download_worker = None
        self._dotplot_worker = None
        
        self._last_dotplot_data = None
        self._plot_window = None

        root = QWidget()
        self.setCentralWidget(root) #permet de définir le widget central de la fenêtre principale, qui contiendra tous les éléments de l'interface utilisateur (comme les champs de recherche, les boutons, le tableau des résultats, et le canvas pour le dotplot).
        layout = QVBoxLayout(root) 

        #on va construire l'interface utilisateur en ajoutant d'abord une barre de recherche en haut, puis un tableau pour afficher les résultats de la recherche NCBI,
        # ensuite des boutons pour définir les génomes A et B, un bouton pour lancer le téléchargement, un autre pour lancer le pipeline du dotplot, 
        #et enfin un canvas pour afficher le dotplot lui-même.

        top = QHBoxLayout()
        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText("Ex : Escherichia coli K-12")
        self.btn_search = QPushButton("Rechercher")
        top.addWidget(self.search_input)
        top.addWidget(self.btn_search)
        layout.addLayout(top)

        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(["Assembly", "Organism", "Level", "Source"])
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.setSortingEnabled(True)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        layout.addWidget(self.table)

        bottom = QHBoxLayout()
        self.label_a = QLabel("Génome A : (aucun)")
        self.label_b = QLabel("Génome B : (aucun)")
        self.btn_set_a = QPushButton("Définir A")
        self.btn_set_b = QPushButton("Définir B")
        bottom.addWidget(self.label_a)
        bottom.addWidget(self.btn_set_a)
        bottom.addSpacing(30)
        bottom.addWidget(self.label_b)
        bottom.addWidget(self.btn_set_b)
        layout.addLayout(bottom)

        actions = QHBoxLayout()
        self.btn_download = QPushButton("Télécharger les génomes sélectionnés")
        self.btn_download.setEnabled(False)
        self.btn_download.setFixedHeight(30)
        self.btn_download.setFixedWidth(300)
        self.status_label = QLabel("Prêt.")
        actions.addWidget(self.btn_download)
        actions.addWidget(self.status_label)
        layout.addLayout(actions)

        filter_row = QHBoxLayout()
        self.evalue_label = QLabel("Seuil E-value : ")
        self.evalue_input = QLineEdit("1e-20")
        self.evalue_input.setMaximumWidth(120)
        self.evalue_input.setToolTip(
            "Deux gènes sont considérés homologues si E-value < seuil"
        )

        filter_row.addWidget(self.evalue_label)
        filter_row.addWidget(self.evalue_input)
        filter_row.addStretch()

        layout.addLayout(filter_row)
        
        
        

        self.btn_plot = QPushButton("Lancer BLASTP + Dotplot")
        
        self.btn_plot.setEnabled(False)
        self.btn_plot.move(300, 200)
        self.btn_plot.setFixedWidth(400)
        
        self.btn_plot.setStyleSheet("""
                                    QPushButton {
                                        background-color: #4f7c82;
                                        color: white;
                                        border : 1px solid #3f666b;
                                        border-radius : 6px;
                                        padding : 6px 12px;
                                        font-weight : 600;
                                    }
                                    QPushButton:hover:!disabled {
                                        background-color : #5b8b91;
                                        
                                    }
                                    QPushButton:pressed:!disabled {
                                        background-color : #446c71;
                                    }
                                    QPushButton:disabled {
                                        background-color: #cfd8da;
                                        color : #7a8688;
                                        border : 1px solid #bcc6c8;
                                    }
                                    """)
        
        
        
        self.plot_status = QLabel("Aucun dotplot affiché.")
        layout.addWidget(self.btn_plot, alignment=Qt.AlignHCenter)
        layout.addWidget(self.plot_status)
        
        

        
        

        self.hover_label = QLabel("Aucune paire de gènes sélectionnée.")
        self.hover_label.setWordWrap(True)
        self.hover_label.setFixedHeight(70)
        self.hover_label.setStyleSheet("QLabel { background-color : #f0f0f0; padding: 4px; border: 1px solid #ccc; }")
        layout.addWidget(self.hover_label)
        
        self.selected_label = QLabel("Sélection : aucune paire de gènes figée.")
        layout.addWidget(self.selected_label)

        self.details_panel = QTextEdit()
        self.details_panel.setReadOnly(True)
        self.details_panel.setMinimumHeight(70)
        self.details_panel.setPlainText(
            "Cliquez sur un point du dotplot pour figer la paire de gènes."
        )
        layout.addWidget(self.details_panel)

        self._selected_point_index = None
        self._selected_marker = None
        
        

        

        self._last_dotplot_data = None
        self._plot_window = None
        self._show_synteny = False

        self.btn_toggle_synteny = QPushButton("Synténies")
        self.btn_toggle_synteny.setCheckable(True)
        self.btn_toggle_synteny.setEnabled(False)
        self.btn_toggle_synteny.setFixedWidth(90)
        self.btn_toggle_synteny.setFixedHeight(26)
        self.btn_toggle_synteny.setToolTip("Afficher / masquer les régions de synténie en rouge")
        self.btn_toggle_synteny.setStyleSheet("""
            QPushButton {
                background-color: #f4f4f4;
                color: #333333;
                border: 1px solid #c8c8c8;
                border-radius: 5px;
                padding: 2px 8px;
                font-size: 11px;
            }
            QPushButton:hover:!disabled {
                background-color: #ececec;
            }
            QPushButton:checked {
                background-color: #f3d6d6;
                color: #7a1f1f;
                border: 1px solid #cc9d9d;
            }
            QPushButton:disabled {
                color: #9a9a9a;
                background-color: #f7f7f7;
                border: 1px solid #dddddd;
            }
        """)

        self.btn_open_plot_window = QPushButton("Agrandir")
        self.btn_open_plot_window.setEnabled(False)
        self.btn_open_plot_window.setFixedWidth(90)
        self.btn_open_plot_window.setFixedHeight(26)

        plot_header = QHBoxLayout()
        plot_header.addWidget(self.btn_toggle_synteny)
        plot_header.addStretch()
        plot_header.addWidget(self.btn_open_plot_window)
        layout.addLayout(plot_header)
        
        
        self.figure = Figure(figsize=(5, 5), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setFixedSize(500, 500)
        self.toolbar = NavigationToolbar(self.canvas, self)
        
        canvas_row=QHBoxLayout()
        canvas_row.addStretch()
        canvas_row.addWidget(self.canvas)
        canvas_row.addStretch()
        layout.addLayout(canvas_row)
        layout.addWidget(self.toolbar)
        
        
        
        self.canvas.mpl_connect("button_press_event", self._on_plot_click)

        self._dotplot_points = []
        self._dotplot_pixel_positions = []
        self._dotplot_ax = None
        self._dotplot_annotation = None
        self._last_hover_index = None

        self.canvas.mpl_connect("motion_notify_event", self._on_plot_hover)
        self.canvas.mpl_connect("figure_leave_event", self._on_plot_leave)
        self.canvas.mpl_connect("draw_event", self._on_canvas_draw)
        

        

        self.btn_search.clicked.connect(self.launch_search) #permet de connecter le clic sur le bouton de recherche à la méthode launch_search, qui va lancer la recherche NCBI en arrière-plan.
        self.search_input.returnPressed.connect(self.launch_search) #permet de connecter la validation du champ de recherche (en appuyant sur Entrée) à la même méthode launch_search, pour que l'utilisateur puisse lancer la recherche soit en cliquant sur le bouton, soit en appuyant sur Entrée.
        self.btn_set_a.clicked.connect(self.set_as_a)
        self.btn_set_b.clicked.connect(self.set_as_b)
        self.btn_download.clicked.connect(self.download_selected_genomes) #permet de connecter le clic sur le bouton de téléchargement à la méthode download_selected_genomes, qui va lancer le téléchargement des génomes sélectionnés en arrière-plan.
        self.btn_plot.clicked.connect(self.launch_dotplot_pipeline)
        self.btn_toggle_synteny.clicked.connect(self.toggle_synteny_display)
        self.btn_open_plot_window.clicked.connect(self.open_dotplot_window)

    def launch_search(self) -> None: #cette fonction permet de lancer la recherche NCBI en arrière-plan lorsque l'utilisateur clique sur le bouton de recherche ou appuie sur Entrée dans le champ de recherche. Elle récupère la requête saisie par l'utilisateur, désactive le bouton de recherche pour éviter les clics multiples, et crée un worker pour exécuter la fonction search_assemblies avec les paramètres appropriés. Les signaux du worker sont connectés à des méthodes pour gérer les résultats, les erreurs, et la fin de la tâche.
        query = self.search_input.text().strip()
        if not query:
            return

        self.btn_search.setEnabled(False)
        self.btn_search.setText("Recherche en cours...")

        self._search_worker = Worker(search_assemblies, self.ncbi, query, SEARCH_LIMITS)
        self._search_worker.signals.result.connect(self._on_search_results)
        self._search_worker.signals.error.connect(self._on_search_error)
        self._search_worker.signals.finished.connect(self._on_search_done)
        self.threadpool.start(self._search_worker)

    def _on_search_results(self, results: list[GenomeCandidate]) -> None: #cette fonction est appelée lorsque la recherche NCBI en arrière-plan a terminé et a renvoyé les résultats. Elle prend en entrée une liste de GenomeCandidate correspondant aux résultats de la recherche, et met à jour le tableau de l'interface utilisateur pour afficher ces résultats. Si aucun résultat n'est trouvé, elle affiche un message d'information. Enfin, elle réactive le tri du tableau et ajuste la taille des lignes.
        self.table.setSortingEnabled(False)  #permet de désactiver temporairement le tri du tableau pendant que nous mettons à jour les données, pour éviter des comportements inattendus.
        self.table.setRowCount(0)

        if not results:
            QMessageBox.information(self, "NCBI", "Aucun résultat trouvé pour cette requête.")
            self.table.setSortingEnabled(True) #permet de réactiver le tri du tableau après avoir mis à jour les données, pour que l'utilisateur puisse trier les résultats par n'importe quelle colonne.
            return

        for g in results: #pour chaque résultat de la recherche (qui est un GenomeCandidate), on ajoute une nouvelle ligne au tableau et on remplit les cellules avec les informations correspondantes (assembly accession, organism name, assembly level, refseq category).
            r = self.table.rowCount()
            self.table.insertRow(r)
            self.table.setItem(r, 0, QTableWidgetItem(g.assembly_accession))
            self.table.setItem(r, 1, QTableWidgetItem(g.organism_name))
            self.table.setItem(r, 2, QTableWidgetItem(g.assembly_level))
            self.table.setItem(r, 3, QTableWidgetItem(g.refseq_category))

        self.table.setSortingEnabled(True)
        self.table.resizeRowsToContents()

    def _on_search_error(self, msg: str) -> None: #cette fonction est appelée lorsque la recherche NCBI en arrière-plan rencontre une erreur. Elle prend en entrée un message d'erreur et affiche une boîte de dialogue critique pour informer l'utilisateur.
        QMessageBox.critical(self, "Erreur NCBI", msg)

    def _on_search_done(self) -> None: #cette fonction est appelée lorsque la recherche NCBI en arrière-plan est terminée, qu'elle ait réussi ou échoué. Elle réactive le bouton de recherche et remet son texte à la valeur par défaut.
        self.btn_search.setEnabled(True)
        self.btn_search.setText("Rechercher")

    def _current_candidate(self) -> GenomeCandidate | None: #cette fonction permet de récupérer le GenomeCandidate correspondant à la ligne actuellement sélectionnée dans le tableau. Elle vérifie d'abord si une ligne est sélectionnée, puis elle lit les valeurs des cellules de cette ligne pour construire et retourner un GenomeCandidate. Si aucune ligne n'est sélectionnée ou si les cellules sont vides, elle retourne None.
        r = self.table.currentRow()
        if r < 0:
            return None
        
        #ici on va récupérer les éléments des cellules de la ligne sélectionnée pour les colonnes correspondantes à l'assembly accession, l'organism name, l'assembly level, et le refseq category. Si l'un de ces éléments est manquant (None), on retourne None. Sinon, on construit un GenomeCandidate avec les valeurs extraites des cellules et on le retourne.
        item0 = self.table.item(r, 0)
        item1 = self.table.item(r, 1)
        item2 = self.table.item(r, 2)
        item3 = self.table.item(r, 3)

        if item0 is None or item1 is None or item2 is None or item3 is None:
            return None

        return GenomeCandidate(
            assembly_accession=item0.text(),
            organism_name=item1.text(),
            assembly_level=item2.text(),
            refseq_category=item3.text(),
        )  

    def _refresh_actions_state(self) -> None:
        ready = self.selected_a is not None and self.selected_b is not None
        self.btn_download.setEnabled(ready)

    def set_as_a(self) -> None:
        cand = self._current_candidate()
        if not cand:
            QMessageBox.warning(self, "Sélection", "Sélectionne d'abord une ligne dans le tableau.")
            return

        self.selected_a = cand
        self.label_a.setText(f"Génome A : {cand.assembly_accession} — {cand.organism_name[:60]}")
        self._refresh_actions_state()

    def set_as_b(self) -> None:
        cand = self._current_candidate()
        if not cand:
            QMessageBox.warning(self, "Sélection", "Sélectionne d'abord une ligne dans le tableau.")
            return

        self.selected_b = cand
        self.label_b.setText(f"Génome B : {cand.assembly_accession} — {cand.organism_name[:60]}")
        self._refresh_actions_state()

    def download_selected_genomes(self) -> None: #Cette fonction est appelée lorsque l'utilisateur clique sur le bouton de téléchargement pour lancer le téléchargement des génomes sélectionnés.
        # Elle vérifie d'abord que les génomes A et B ont été sélectionnés, puis elle désactive les boutons d'action pour éviter les interactions pendant le téléchargement, et met à jour le statut pour informer l'utilisateur que le téléchargement est en cours. 
        # Ensuite, elle crée un worker pour exécuter la fonction download_two_genomes avec les paramètres appropriés, et connecte les signaux du worker à des méthodes pour gérer les résultats, les erreurs, et la fin de la tâche.
        if self.selected_a is None or self.selected_b is None:
            QMessageBox.warning(self, "Téléchargement", "Sélectionne d'abord les génomes A et B.")
            return

        self.btn_download.setEnabled(False)
        self.btn_search.setEnabled(False)
        self.btn_set_a.setEnabled(False)
        self.btn_set_b.setEnabled(False)
        self.btn_plot.setEnabled(False)
        self.status_label.setText("Téléchargement en cours...")

        self._download_worker = Worker(
            download_two_genomes,
            self.ncbi,
            self.selected_a,
            self.selected_b,
            self.cache_dir,
        )
        self._download_worker.signals.result.connect(self._on_download_results)
        self._download_worker.signals.error.connect(self._on_download_error)
        self._download_worker.signals.finished.connect(self._on_download_done)
        self.threadpool.start(self._download_worker)

    def _on_download_results(self, result) -> None:
        #Cette fonction est appelée lorsque le téléchargement des génomes en arrière-plan a terminé avec succès. 
        # Elle prend en entrée le résultat du téléchargement, qui contient les packages de données pour les génomes A et B. 
        # Elle stocke ces packages dans les attributs self.downloaded_a et self.downloaded_b, puis met à jour le statut pour informer l'utilisateur que les génomes ont été téléchargés, et active le bouton pour lancer le pipeline du dotplot.
        pkg_a, pkg_b = result
        self.downloaded_a = pkg_a 
        self.downloaded_b = pkg_b

        self.status_label.setText(
            f"Téléchargés : {pkg_a.assembly_accession} et {pkg_b.assembly_accession}"
        )
        self.btn_plot.setEnabled(True)

        QMessageBox.information(
            self,
            "Téléchargement terminé",
            "Les génomes ont été téléchargés avec succès.",
        )

    def _on_download_error(self, msg: str) -> None:
        self.status_label.setText("Erreur de téléchargement.")
        QMessageBox.critical(self, "Téléchargement", msg)

    def _on_download_done(self) -> None:
        self.btn_search.setEnabled(True)
        self.btn_set_a.setEnabled(True)
        self.btn_set_b.setEnabled(True)
        self._refresh_actions_state()

    def launch_dotplot_pipeline(self) -> None:
        #Cette fonction est appelée lorsque l'utilisateur clique sur le bouton pour lancer le pipeline du dotplot.
        # Elle vérifie d'abord que les données des génomes A et B ont été téléchargées, puis elle désactive le bouton de lancement du pipeline et met à jour le statut pour informer
        # l'utilisateur que le BLASTP et la construction du dotplot sont en cours. Ensuite, elle crée un worker pour exécuter la fonction compute_dotplot_data avec les paramètres appropriés, et connecte les signaux du worker à des méthodes pour gérer les résultats, les erreurs, et la fin de la tâche.
        
        
        
        
        
        if self.downloaded_a is None or self.downloaded_b is None:
            QMessageBox.warning(self, "Dotplot", "Téléchargez d'abord les deux génomes.")
            return

        raw_evalue = self.evalue_input.text().strip()
        try :
            evalue_treshold = float(raw_evalue)
            if evalue_treshold <= 0:
                raise ValueError
        except ValueError :
            QMessageBox.warning(self, "Paramètre invalide", "Le seuil E-value doit être un nombre positif valide.")
            return


        self._show_synteny = False
        self.btn_toggle_synteny.setChecked(False)
        self.btn_toggle_synteny.setEnabled(False)
        
        self._action_export_hits.setEnabled(False)
        self._action_export_synteny.setEnabled(False)
        self._action_export_summary.setEnabled(False)
        
        
        self.btn_plot.setEnabled(False)
        self.plot_status.setText("BLASTP et construction du dotplot en cours...")
     
        work_dir = Path("cache/blast") / (
            f"{self.downloaded_a.assembly_accession}_vs_{self.downloaded_b.assembly_accession}"
        )

        self._dotplot_worker = Worker(
            compute_dotplot_data,
            str(self.downloaded_a.proteins_faa),
            str(self.downloaded_b.proteins_faa),
            str(self.downloaded_a.feature_table),
            str(self.downloaded_b.feature_table),
            str(work_dir),
            evalue_treshold,
            10,
            2,
            True,
        )
        self._dotplot_worker.signals.result.connect(self._on_dotplot_ready)
        self._dotplot_worker.signals.error.connect(self._on_dotplot_error)
        self._dotplot_worker.signals.finished.connect(self._on_dotplot_done)
        self.threadpool.start(self._dotplot_worker)

    def _on_dotplot_ready(self, result: dict) -> None:
        self._last_dotplot_data = result
        self.btn_open_plot_window.setEnabled(True)
        self.btn_toggle_synteny.setEnabled(True)
        self._action_export_hits.setEnabled(True)
        self._action_export_synteny.setEnabled(True)
        self._action_export_summary.setEnabled(True)
        self._render_dotplot(result)

    def toggle_synteny_display(self) -> None:
        self._show_synteny = self.btn_toggle_synteny.isChecked()
        if self._last_dotplot_data is not None:
            self._render_dotplot(self._last_dotplot_data)

    def _render_dotplot(self, result: dict) -> None:
        points = result["points"]
        n_genes_a = result["n_genes_a"]
        n_genes_b = result["n_genes_b"]
        n_hits = result["n_hits"]
        synteny_segments = result.get("synteny_segments", [])

        self._dotplot_points = points
        self._selected_point_index = None
        self.selected_label.setText("Sélection : aucune paire de gènes figée.")
        self.details_panel.setPlainText(
            "Sélectionnez un point du dotplot pour figer la paire de gènes."
        )
        self._dotplot_pixel_positions = []
        self._last_hover_index = None
        self.hover_label.setText("Survol : aucune paire de gènes.")
        self._selected_marker = None

        self.figure.clear()
        ax = self.figure.add_subplot(111)
        self._dotplot_ax = ax

        if points:
            xs = [p.x for p in points]
            ys = [p.y for p in points]
            ax.scatter(
                xs,
                ys,
                s=4,
                c="black",
                alpha=0.8,
                linewidths=0,
                marker="s",
            )

        if self._show_synteny:
            for segment in synteny_segments:
                seg_x = [p.x for p in segment.points]
                seg_y = [p.y for p in segment.points]
                ax.plot(seg_x, seg_y, color="red", linewidth=0.9, alpha=0.95, zorder=3)
                ax.scatter(seg_x, seg_y, s=10, c="red", alpha=0.95, linewidths=0, zorder=4)

        ax.set_xlim(-1, max(n_genes_a, 1))
        ax.set_ylim(-1, max(n_genes_b, 1))
        ax.set_xlabel("Gènes du génome A")
        ax.set_ylabel("Gènes du génome B")
        ax.set_title(
            f"Dotplot {self.downloaded_a.assembly_accession} vs {self.downloaded_b.assembly_accession}"
        )
        ax.grid(False)

        self._dotplot_annotation = ax.annotate(
            "",
            xy=(0, 0),
            xytext=(12, 12),
            textcoords="offset points",
            bbox=dict(boxstyle="round", fc="white", ec="black", alpha=0.95),
            arrowprops=dict(arrowstyle="->", color="black"),
        )
        self._dotplot_annotation.set_visible(False)

        self.figure.tight_layout()
        self.canvas.draw()
        self._update_hover_cache()

        if self._show_synteny:
            n_seg = result.get("n_synteny_segments", len(synteny_segments))
            self.plot_status.setText(
                f"Dotplot affiché avec succès. Hits retenus : {n_hits} | {n_seg} segment(s) de synténie visibles"
            )
        else:
            self.plot_status.setText(
                f"Dotplot affiché avec succès. Hits retenus : {n_hits}"
            )
    def _export_hits_csv(self) -> None:
        if self._last_dotplot_data is None:
            return
        points = self._last_dotplot_data["points"]
        default_name = (
            f"hits_{self.downloaded_a.assembly_accession}_vs_"
            f"{self.downloaded_b.assembly_accession}.csv"
        )
        path, _ = QFileDialog.getSaveFileName(
            self, "Exporter les hits BLASTP", default_name, "CSV (*.csv)"
        )
        if not path:
            return
        try:
            with open(path, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow([
                    "x_genome_a", "y_genome_b",
                    "qseqid", "sseqid",
                    "q_locus_tag", "s_locus_tag",
                    "q_name", "s_name",
                    "evalue", "bitscore",
                ])
                for p in points:
                    writer.writerow([
                        p.x, p.y,
                        p.qseqid, p.sseqid,
                        p.q_locus_tag or "", p.s_locus_tag or "",
                        p.q_name or "", p.s_name or "",
                        p.evalue, p.bitscore,
                    ])
            QMessageBox.information(self, "Export réussi", f"Hits exportés dans :\n{path}")
        except Exception as e:
            QMessageBox.critical(self, "Erreur export", str(e))

    def _export_synteny_csv(self) -> None:
        if self._last_dotplot_data is None:
            return
        segments = self._last_dotplot_data.get("synteny_segments", [])
        default_name = (
            f"syntenie_{self.downloaded_a.assembly_accession}_vs_"
            f"{self.downloaded_b.assembly_accession}.csv"
        )
        path, _ = QFileDialog.getSaveFileName(
            self, "Exporter les segments de synténie", default_name, "CSV (*.csv)"
        )
        if not path:
            return
        try:
            with open(path, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow([
                    "segment_id", "n_points",
                    "x_start", "y_start", "x_end", "y_end",
                    "point_index", "x", "y",
                    "qseqid", "sseqid",
                    "q_locus_tag", "s_locus_tag",
                    "q_name", "s_name",
                    "evalue", "bitscore",
                ])
                for seg_id, segment in enumerate(segments, start=1):
                    pts = segment.points
                    x_start = pts[0].x if pts else ""
                    y_start = pts[0].y if pts else ""
                    x_end = pts[-1].x if pts else ""
                    y_end = pts[-1].y if pts else ""
                    for pt_idx, p in enumerate(pts):
                        writer.writerow([
                            seg_id, len(pts),
                            x_start, y_start, x_end, y_end,
                            pt_idx,
                            p.x, p.y,
                            p.qseqid, p.sseqid,
                            p.q_locus_tag or "", p.s_locus_tag or "",
                            p.q_name or "", p.s_name or "",
                            p.evalue, p.bitscore,
                        ])
            QMessageBox.information(self, "Export réussi", f"Synténies exportées dans :\n{path}")
        except Exception as e:
            QMessageBox.critical(self, "Erreur export", str(e))

    def _export_summary_txt(self) -> None:
        if self._last_dotplot_data is None:
            return
        result = self._last_dotplot_data
        acc_a = self.downloaded_a.assembly_accession
        acc_b = self.downloaded_b.assembly_accession
        org_a = getattr(self.downloaded_a, "organism_name", "—")
        org_b = getattr(self.downloaded_b, "organism_name", "—")
        evalue = self.evalue_input.text().strip()
        n_hits = result["n_hits"]
        n_genes_a = result["n_genes_a"]
        n_genes_b = result["n_genes_b"]
        n_seg = result.get("n_synteny_segments", len(result.get("synteny_segments", [])))
        cache_used = result.get("blast_cache_used", False)
        blast_path = result.get("blast_result_path", "—")
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        lines = [
            "=" * 60,
            "  RÉSUMÉ DE L'ANALYSE DOTPLOT",
            "=" * 60,
            f"  Date                : {now}",
            "",
            "  GÉNOMES COMPARÉS",
            f"  Génome A            : {acc_a}",
            f"  Organisme A         : {org_a}",
            f"  Génome B            : {acc_b}",
            f"  Organisme B         : {org_b}",
            "",
            "  PARAMÈTRES",
            f"  Seuil E-value       : {evalue}",
            f"  Fichier BLAST       : {blast_path}",
            f"  Cache BLAST utilisé : {'Oui' if cache_used else 'Non'}",
            "",
            "  RÉSULTATS",
            f"  Gènes génome A      : {n_genes_a}",
            f"  Gènes génome B      : {n_genes_b}",
            f"  Hits retenus        : {n_hits}",
            f"  Segments de synténie: {n_seg}",
            "=" * 60,
        ]
        summary_text = "\n".join(lines)

        default_name = f"resume_{acc_a}_vs_{acc_b}.txt"
        path, _ = QFileDialog.getSaveFileName(
            self, "Sauvegarder le résumé", default_name, "Texte (*.txt)"
        )
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8") as f:
                f.write(summary_text)
            QMessageBox.information(self, "Sauvegarde réussie", f"Résumé sauvegardé dans :\n{path}")
        except Exception as e:
            QMessageBox.critical(self, "Erreur sauvegarde", str(e))
    
    def _on_dotplot_error(self, msg: str) -> None:
        self.plot_status.setText("Erreur lors du calcul du dotplot.")
        QMessageBox.critical(self, "Dotplot", msg)

    def _on_dotplot_done(self) -> None:
        self.btn_plot.setEnabled(True)
        
    def _format_hover_text(self, point) -> str:
        q_locus = point.q_locus_tag or "-"
        s_locus = point.s_locus_tag or "-"
        q_name = point.q_name or "-"
        s_name = point.s_name or "-"

        return (
            f"A : qseqid={point.qseqid} | locus_tag={q_locus} | protéine={q_name}\n"
            f"B : sseqid={point.sseqid} | locus_tag={s_locus} | protéine={s_name}\n"
            f"BLAST : E-value={point.evalue:.2e} | bitscore={point.bitscore:.1f}"
        )

    def _update_hover_cache(self) -> None:
        self._dotplot_pixel_positions = []

        if self._dotplot_ax is None or not self._dotplot_points:
            return

        transformed = self._dotplot_ax.transData.transform(
            [(p.x, p.y) for p in self._dotplot_points]
        )
        self._dotplot_pixel_positions = [
            (float(px), float(py)) for px, py in transformed
        ]

    def _on_canvas_draw(self, event) -> None:
        self._update_hover_cache()

    def _find_hovered_point_index(self, event, max_pixel_distance: float = 8.0):
        if self._dotplot_ax is None:
            return None
        if event.inaxes is not self._dotplot_ax:
            return None
        if event.x is None or event.y is None:
            return None
        if not self._dotplot_points or not self._dotplot_pixel_positions:
            return None

        mouse_x = float(event.x)
        mouse_y = float(event.y)

        best_index = None
        best_d2 = max_pixel_distance * max_pixel_distance

        for i, (px, py) in enumerate(self._dotplot_pixel_positions):
            d2 = (px - mouse_x) ** 2 + (py - mouse_y) ** 2
            if d2 <= best_d2:
                best_d2 = d2
                best_index = i

        return best_index

    def _clear_hover_feedback(self) -> None:
        if self._last_hover_index is None and (
            self._dotplot_annotation is None
            or not self._dotplot_annotation.get_visible()
        ):
            return

        self._last_hover_index = None
        self.hover_label.setText("Survol : aucune paire de gènes.")

        if self._dotplot_annotation is not None:
            self._dotplot_annotation.set_visible(False)

        self.canvas.draw_idle()

    def _on_plot_hover(self, event) -> None:
        point_index = self._find_hovered_point_index(event)

        if point_index is None:
            self._clear_hover_feedback()
            return

        if point_index == self._last_hover_index:
            return

        self._last_hover_index = point_index
        point = self._dotplot_points[point_index]

        self.hover_label.setText(self._format_hover_text(point))

        if self._dotplot_annotation is not None:
            q_short = point.q_locus_tag or point.qseqid
            s_short = point.s_locus_tag or point.sseqid
            self._dotplot_annotation.xy = (point.x, point.y)
            self._dotplot_annotation.set_text(f"{q_short} ↔ {s_short}")
            self._dotplot_annotation.set_visible(True)

        self.canvas.draw_idle()

    def _on_plot_leave(self, event) -> None:
        self._clear_hover_feedback()
        
    def _format_selected_text(self, point) -> str:
        
        q_locus = point.q_locus_tag or "-"
        s_locus = point.s_locus_tag or "-"
        q_name = point.q_name or "-"
        s_name = point.s_name or "-"

        return (
            "Paire de gènes sélectionnée\n"
            "--------------------------\n"
            f"Génome A\n"
            f"  - qseqid      : {point.qseqid}\n"
            f"  - locus_tag   : {q_locus}\n"
            f"  - protéine    : {q_name}\n"
            f"\n"
            f"Génome B\n"
            f"  - sseqid      : {point.sseqid}\n"
            f"  - locus_tag   : {s_locus}\n"
            f"  - protéine    : {s_name}\n"
            f"\n"
            f"BLAST\n"
            f"  - E-value     : {point.evalue:.2e}\n"
            f"  - bitscore    : {point.bitscore:.1f}\n"
            f"\n"
            f"Coordonnées dotplot\n"
            f"  - x           : {point.x}\n"
            f"  - y           : {point.y}"
        )

    def _draw_selected_marker(self, point) -> None:
        if self._dotplot_ax is None:
            return

        if self._selected_marker is not None:
            try:
                self._selected_marker.remove()
            except ValueError:
                pass
            self._selected_marker = None

        self._selected_marker = self._dotplot_ax.scatter(
            [point.x],
            [point.y],
            s=80,
            facecolors="none",
            edgecolors="red",
            linewidths=1.5,
            zorder=5,
        )

    def _show_selected_point(self, point_index: int) -> None:
        if point_index < 0 or point_index >= len(self._dotplot_points):
            return

        self._selected_point_index = point_index
        point = self._dotplot_points[point_index]

        short_q = point.q_locus_tag or point.qseqid
        short_s = point.s_locus_tag or point.sseqid

        self.selected_label.setText(
            f"Sélection : {short_q} ↔ {short_s}"
        )
        self.details_panel.setPlainText(self._format_selected_text(point))
        self._draw_selected_marker(point)
        self.canvas.draw_idle()

    def _clear_selected_point(self) -> None:
        self._selected_point_index = None
        self.selected_label.setText("Sélection : aucune paire de gènes figée.")
        self.details_panel.setPlainText(
            "Clique sur un point du dotplot pour figer la paire de gènes."
        )

        if self._selected_marker is not None:
            try:
                self._selected_marker.remove()
            except ValueError:
                pass
            self._selected_marker = None

        self.canvas.draw_idle()

    def _on_plot_click(self, event) -> None:
        if self._dotplot_ax is None:
            return

        if event.inaxes is not self._dotplot_ax:
            return

        if event.button == 3:
            self._clear_selected_point()
            return

        if event.button != 1:
            return

        point_index = self._find_hovered_point_index(event, max_pixel_distance=10.0)
        if point_index is None:
            return

        self._show_selected_point(point_index)
        
    def open_dotplot_window(self) -> None:
        if not self._last_dotplot_data :
            QMessageBox.information(self, "Dotplot", "Aucun dotplot à afficher.")
            return
        
        self._plot_window = DotplotWindow(
            dotplot_data = self._last_dotplot_data,
            title=f"Dotplot {self.downloaded_a.assembly_accession} vs {self.downloaded_b.assembly_accession}"
        )
        self._plot_window.showMaximized()
        