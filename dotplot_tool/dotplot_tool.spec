# -*- mode: python ; coding: utf-8 -*-
# ============================================================
#  dotplot_tool.spec  —  PyInstaller (Linux x86-64 / WSL2)
#  Projet : DotPlot NCBI  —  Malek Louiz 2026 UM4BM748
# ============================================================
#
#  Place ce fichier dans :
#      bioinfo_dotplot/src/dotplot_tool/dotplot_tool.spec
#
#  Lance ensuite, depuis ce dossier :
#      pyinstaller dotplot_tool.spec
#
#  Résultat :
#      bioinfo_dotplot/src/dotplot_tool/dist/dotplot_tool
#
#  Notes :
#  - app.py est dans le même dossier que ce .spec
#  - les sous-modules du projet sont : blast, core, dotplot, genome, ncbi, ui, utils
#  - Pillow (PIL) est inclus car Matplotlib en dépend à l'exécution
#  - tkinter reste exclu
# ============================================================

from PyInstaller.utils.hooks import collect_all

mpl_datas, mpl_bins, mpl_hidden = collect_all("matplotlib")
pyside_datas, pyside_bins, pyside_hidden = collect_all("PySide6")
try:
    pil_datas, pil_bins, pil_hidden = collect_all("PIL")
except Exception:
    pil_datas, pil_bins, pil_hidden = [], [], []


a = Analysis(
    ["app.py"],
    pathex=["."],

    binaries=(
        mpl_bins
        + pyside_bins
        + pil_bins
        # Décommentez si vous voulez embarquer BLAST+ dans l'exécutable :
        # + [("/usr/bin/blastp", "bin")]
        # + [("/usr/bin/makeblastdb", "bin")]
    ),

    datas=(
        mpl_datas
        + pyside_datas
        + pil_datas
        + [
            ("blast", "blast"),
            ("core", "core"),
            ("dotplot", "dotplot"),
            ("genome", "genome"),
            ("ncbi", "ncbi"),
            ("ui", "ui"),
            ("utils", "utils"),
        ]
    ),

    hiddenimports=(
        mpl_hidden
        + pyside_hidden
        + pil_hidden
        + [
            # Matplotlib
            "matplotlib.backends.backend_qtagg",
            "matplotlib.backends.backend_agg",
            "matplotlib.backends.backend_svg",
            "matplotlib.figure",

            # PySide6
            "PySide6.QtCore",
            "PySide6.QtWidgets",
            "PySide6.QtGui",
            "PySide6.QtSvg",
            "PySide6.QtXml",
            "PySide6.QtOpenGL",
            "PySide6.QtOpenGLWidgets",
            "PySide6.QtPrintSupport",

            # Pillow
            "PIL",
            "PIL.Image",
            "PIL.ImageFile",
            "PIL.ImageOps",

            # Sous-modules du projet
            "blast",
            "blast.pipeline",
            "core",
            "core.models",
            "dotplot",
            "dotplot.builder",
            "dotplot.synteny",
            "dotplot.render",
            "genome",
            "genome.parser",
            "ncbi",
            "ncbi.client",
            "ncbi.search",
            "ncbi.download",
            "ui",
            "ui.main_window",
            "utils",
            "utils.threads",

            # Stdlib utilisée
            "csv",
            "pathlib",
            "datetime",
            "subprocess",
            "shutil",
        ]
    ),

    excludes=[
        "tkinter", "_tkinter",
        "PyQt5", "PyQt6",
        "wx",
        "IPython", "notebook",
        "scipy", "sklearn",
        "pandas", "cv2",
    ],

    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    noarchive=False,
    optimize=0,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name="dotplot_tool",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=["libpyside6*", "libshiboken*"],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
