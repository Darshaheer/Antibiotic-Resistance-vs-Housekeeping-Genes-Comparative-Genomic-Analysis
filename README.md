# Antibiotic-Resistance-vs-Housekeeping-Genes-Comparative-Genomic-Analysis
This project integrates Python and R for sequence analysis. Python extracts sequences, performs multiple alignments, and calculates GC content, while R visualizes GC distribution and nucleotide bias, providing insights into sequence composition and base-level variation across samples.

## Repository Structure
Sequence-Alignment-and-Visualization/
│── AMR.py                       # Python script for sequence parsing and alignment
│── AMR.R                        # R script for data visualization
│── results/                     # Output folder
│   └── (alignment data, GC content CSVs, and plots)
│── README.md                    # Project documentation
│── LICENSE                      # License file (if included)
│── .gitignore                   # Ignore unnecessary files

## Features

Reads and processes FASTA files
  Performs multiple sequence alignment
  Calculates GC content and nucleotide frequency per sequence
  Exports alignment statistics and summary files
  Visualizes GC content distribution and nucleotide bias using R

## Requirements

  Python 3.10+
  Biopython
  R 4.3+
  R packages: ggplot2, dplyr, readr

**Install Python dependencies with:**

  pip install biopython matplotlib

**Install R dependencies with:**

  install.packages(c("ggplot2", "dplyr", "readr"))

## Usage

**Prepare Input**

Run Python Script
python AMR.py
This script performs alignment, calculates GC content, and exports results to results/.

**Run R Script for Visualization**
source("AMR.R")
This generates GC content and nucleotide bias plots, saved in the results/ folder.

## Example Output
  Alignment summary text file
  GC content statistics CSV
  Nucleotide frequency table

**Plots:**
  GC Content Distribution
  Nucleotide Bias Plot

## Notes
Replace the sample Genbank Id's with your own Ids.
This project can be extended for comparative genomic studies, motif analysis, or evolutionary conservation visualization.

## License
MIT License
This project is open-source. You are free to use and modify it for research or educational purposes.
