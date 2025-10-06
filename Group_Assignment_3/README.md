# Protein Sequence Classification
## Project Overview

- This project explores the use of machine learning to classify protein sequences based solely on their amino acid composition. Specifically, it focuses on two protein families: Globins and Zinc Fingers, using a Feature Frequency Profile (FFP) representation derived from the 2-mer (dipeptide) composition of each sequence.
- The main objective is to build and evaluate predictive models capable of distinguishing between these two protein families using only sequence data.

## Objectives

- Parse and represent protein sequences from FASTA files.
- Compute 2-mer frequency profiles (400 features) for each sequence.
- Construct a labelled dataset with class identifiers (Globin = 0, Zinc Finger = 1).
- Apply and compare multiple machine learning classifiers (e.g., SVM, Random Forest, Naive Bayes).
- Use cross-validation to evaluate models based on metrics such as precision, recall, and F1-score.
- Identify and report the best-performing model.
