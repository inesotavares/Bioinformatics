# Phylogenetic Analysis of the Hemoglobin Gene
## Project Overview
This project explores the evolutionary relationships of the hemoglobin beta subunit (HBB_HUMAN, UniProt P68871) across different species using a computational phylogenetics pipeline.

The analysis follows a standard bioinformatics workflow:

- Retrieving protein sequence data from UniProt

- Performing BLAST to identify homologous sequences from other species

- Selecting representative sequences (10 non-human species + human)

- Conducting multiple sequence alignment (MSA) with Clustal Omega

- Building and visualizing phylogenetic trees

- Comparing evolutionary patterns and drawing conclusions

The final deliverable is a Python script (build_phylotree.py) that automates the pipeline from sequence retrieval to tree construction.
