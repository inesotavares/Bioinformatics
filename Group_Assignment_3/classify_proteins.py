### HEADER

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@   Bioinformatics --- Group Assignment 3  @@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@ By:   -> Bruna Rocha (up201906417 M:CC)          @@@
# @@@       -> Inês Tavares (up201706579 M:DS)         @@@
# @@@       -> Nuno Oliveira (up201703852 M:BBC)       @@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import itertools
from sklearn.dummy import DummyClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB

def read_fasta(seqfile: str) -> dict[str, str]:
    """Reads a multi-sequence FASTA file to a dictionary."""
    sequences = {}
    with open(seqfile, 'r') as fasta_file:
        sequence_id = ''
        sequence = ''
        for i in fasta_file:
            i=i.strip()
            if i.startswith('>'):
                if sequence_id:
                    sequences[sequence_id] =sequence
                    sequence = ''
                sequence_id = i.split('|')[1]                
            else:
                sequence += i
        if sequence_id:
            sequences[sequence_id] = sequence
    return sequences


def ffp(sequence: str) -> dict[str, float]:
    """Compute the FFP representation of one sequence as a dictionary
    (KEY: 2-mer | VALUE: '2-mer freq / total 2-mer count')."""
    total_number_two_mers = len(sequence) / 2
    two_mers_dict = {two_mer: 0 for two_mer in two_mers}
    for i in range(0, len(sequence)-1, 2):
        two_mer = sequence[i:i+2]
        if two_mer in two_mers_dict:
            two_mers_dict[two_mer] += 1
    result = {two_mer: count / total_number_two_mers for two_mer, count in two_mers_dict.items()}
    return result


def table_format(dicts: list[dict[str, str]]):
    """Converts FASTA data dictionaries to a dataframe."""
    # zincfinger = 0, globin = 1
    global df
    for i, protein_family in enumerate(dicts):
        for seq_id, sequence in protein_family.items():
            two_mers_count = ffp(sequence)
            two_mers_count['protFamily'] = i
            temp_df = pd.DataFrame([two_mers_count], index=[seq_id])
            df = pd.concat([df, temp_df], ignore_index=False)


def get_cv_metrics(model, X, y, cv_iter):
    """Computes metrics from a given model evaluated using cross-validation."""
    acc_scores = cross_val_score(model, X=X, y=y, cv=cv_iter, scoring='accuracy')
    prec_scores = cross_val_score(model, X=X, y=y, cv=cv_iter, scoring='precision')
    rec_scores = cross_val_score(model, X, y, cv=cv_iter, scoring='recall')
    f1_scores = cross_val_score(model, X, y, cv=cv_iter, scoring='f1')
    print("| (Shown as ==> Avg +- StdDev)")
    mean_acc = acc_scores.mean()
    acc_stdev = acc_scores.std()
    print(f"| Accuracy: {mean_acc:.3f} +- {acc_stdev:.3f}")
    mean_prec = prec_scores.mean()
    prec_stdev = prec_scores.std()
    print(f"| Precision: {mean_prec:.3f} +- {prec_stdev:.3f}")
    mean_rec = rec_scores.mean()
    rec_stdev = rec_scores.std()
    print(f"| Recall: {mean_rec:.3f} +- {rec_stdev:.3f}")
    mean_f1 = f1_scores.mean()
    f1_stdev = f1_scores.std()
    print(f"| F1-score: {mean_f1:.3f} +- {f1_stdev:.3f}")
    print()


def main():
    # Command line in Windows|Linux: "python|python3 classify_proteins.py -a zincfinger.fasta -b globin.fasta -k 2"
    parser = argparse.ArgumentParser(description="Classify proteins based on k-mers.")
    parser.add_argument('-a', type=str, required=True, help='Input FASTA file for protein family A')  # zincfinger.fasta
    parser.add_argument('-b', type=str, required=True, help='Input FASTA file for protein family B')  # globin.fasta
    parser.add_argument('-k', type=int, required=True, help='Length of the k-mer')  # 2
    args = parser.parse_args()
    
    ### Exercise 1
    zincfinger = read_fasta(args.a)  # 284 entries
    globin = read_fasta(args.b)  # 1161 entries
    
    ### Exercises 2 / 4
    amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                   'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    combinations = itertools.product(amino_acids, repeat=2)
    global two_mers, df
    two_mers = [''.join(comb) for comb in combinations]
    cols = two_mers + ['protFamily']
    df = pd.DataFrame(columns=cols)
    
    ### Exercises 3 / 4
    table_format([zincfinger, globin])
    print("> Preview of the dataframe ==> SeqIDs as indexes:")
    print(df)
    print()
    # Create a CSV of the dataframe for classification
    # | Necessary! |
    # | -> 'StratifiedKFold' does not read column 'dtypes' correctly from the 'df' variable! |
    df.to_csv("df.csv", index=False)
    
    ### Exercise 5
    data = pd.read_csv("df.csv")
    # Sets for classification    
    X = data.iloc[:, :-1]
    y = data['protFamily']
    # Heatmap overview of the data
    # | Used for visual analysis of data distribution |
    _, ax = plt.subplots(layout="constrained")
    sns.heatmap(X, vmin=0.0, vmax=0.02, cmap='coolwarm', robust=True,
                xticklabels=False, yticklabels=False, ax=ax)
    ax.set_title("FFPs across Zinc finger and Globin protein sequence families", fontsize=10, weight='bold')
    ax.set_xlabel("2-mer FFPs", fontsize=10)
    ax.set_ylabel("Sequences", fontsize=10)
    plt.savefig("df_heatmap.png", dpi=190.0)
    # Baseline - "Dummy" classifier
    # (Always chooses the most frequent class)
    dummy_clf = DummyClassifier(strategy="most_frequent")
    dummy_clf.fit(X, y)
    bline_acc = dummy_clf.score(X, y)
    print(f"> Baseline Accuracy ==> {bline_acc:.3f}")
    print()
    # 3 ML algorithms -> SVM | Random Forests | Naive Bayes
    svm_clf = svm.SVC(kernel='linear', C=1, random_state=42)
    rf_clf = RandomForestClassifier(random_state=42)
    gnb = GaussianNB()
    # Apply 10-fold cross-validation using Stratified k-fold
    cross_val = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
    # Get metrics for the 3 models:
    print("> SVM classifier metrics across the 10 folds:")
    get_cv_metrics(model=svm_clf, X=X, y=y, cv_iter=cross_val)
    print("> Random Forest classifier metrics across the 10 folds:")
    get_cv_metrics(model=rf_clf, X=X, y=y, cv_iter=cross_val)
    print("> Gaussian naive Bayes classifier metrics across the 10 folds:")
    get_cv_metrics(model=gnb, X=X, y=y, cv_iter=cross_val)


if __name__ == "__main__":
    main()
    
