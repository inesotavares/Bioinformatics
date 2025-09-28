# Yeast ORF Analysis
## Project Overview
This project analyses the first 30 kb of Saccharomyces cerevisiae chromosome I to identify open reading frames (ORFs) and compare them with reference genome annotations.

The workflow includes:

- Computing sequence statistics (length, nucleotide composition, GC content, codon usage).

- Detecting ORFs in both DNA strands with minimum length filtering.

- Translating ORFs into protein sequences (FASTA format).

- Saving ORF genomic coordinates in a structured text file.

- Comparing predicted ORFs against annotations (GTF file) and calculating overlap percentages.

The final deliverable is a Python script (yeast_orfs_chr1.py) that automates the entire analysis.
