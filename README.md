## Graphical Tool for Genome Comparison by Dotplot
Overview
This application allows the comparison of two prokaryotic genomes based on their proteomes. It searches for and downloads the data from NCBI, runs a BLASTP comparison, and then displays the result as an interactive dotplot.

The project specification specifically requires the selection of two prokaryotic genomes, their download from NCBI, the detection of homologous gene pairs using BLASTP, and the display of a dotplot. It also recommends several ergonomic improvements, notably adjusting the E-value threshold, detecting synteny blocks, and interacting with the graph; these features are supported by the application.

Main Features
Search for prokaryotic assemblies in the NCBI database.

Select two genomes to compare.

Download the files required for the analysis.

Run a BLASTP comparison between the genomes.

Build and display the dotplot of homologous gene pairs.

Adjust the E-value threshold before computation.

Optionally display synteny regions.

Show detailed information about a gene pair on hover or click on the dotplot.

Export BLASTP hits in CSV format.

Export synteny segments in CSV format.

Save a textual summary of the analysis.

Required Environment
Development and testing were carried out under Linux. The application requires Python, the project’s Python dependencies, as well as BLAST+ for running blastp and makeblastdb.

Project Structure
text

Malek_Louiz_Dotplot_Project/
├── README.md
├── fiche_resume.md
├── requirements.txt
├── src/
│   └── dotplot_tool/
│       ├── app.py
│       ├── blast/
│       ├── core/
│       ├── dotplot/
│       ├── genome/
│       ├── ncbi/
│       ├── ui/
│       └── utils/
    └── app.py
└── examples/
    └── interface_screenshot.png
Minimum Dependencies
Python 3

PySide6

matplotlib

requests

BLAST+

Example of BLAST+ installation on Debian/Ubuntu:

bash

sudo apt install ncbi-blast+
Installation
From the project root:

bash

python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
Launch
bash

cd src
python app.py
Usage
Open the application.

Search for a first genome in NCBI.

Select it as genome A.

Search for a second genome.

Select it as genome B.

Download the required data.

Adjust the E-value threshold if needed.

Run the comparison pipeline.

Examine the generated dotplot.

Enable or disable the display of synteny regions.

Hover over or click on a point to view information about a gene pair.

Export the results if needed.

Generated Results
The application can produce several useful outputs for evaluation:

an interactive dotplot;

a CSV export of BLASTP hits;

a CSV export of synteny segments;

a textual summary of the analysis.

Possible Evaluation Points
The software makes it possible to directly verify the core elements expected in the project:

selection and download of two prokaryotic genomes;

comparison of their proteins using BLASTP;

visualization of homologs in a dotplot;

addition of useful features in a research context, such as threshold adjustment, synteny highlighting, and interactive point exploration.

Notes
The project statement reminds that the dotplot is particularly useful for visualizing syntenic regions and suggesting rearrangements such as duplications, inversions, or translocations. The tool was designed with this logic of visual exploration and usability in mind, for simple and user-friendly operation.

Author
Malek Louiz
21108469 MU5BM748 — Advanced Python
Academic Year 2025–2026
