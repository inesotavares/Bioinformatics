### HEADER

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@   Bioinformatics --- Group Assignment 1  @@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@ By:   -> Bruna Rocha (up201906417 M:CC)          @@@
# @@@       -> Inês Tavares (up...)                    @@@
# @@@       -> Nuno Oliveira (up201703852)             @@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

import sys 
import pandas as pd
import csv

### A. Get statistics

def read_fasta_file(file_path):
    try:
        with open(file_path, 'r') as file:
            # Skip the header lines starting with '>'
            lines = [line.strip() for line in file if not line.startswith('>')]
            # Concatenate the lines to form a single string
            fasta_string = ''.join(lines)
            return fasta_string
    except FileNotFoundError:
        print("File not found.")


def validate_dna(dna_seq: str) -> bool:
    """Checks if DNA sequence is valid. Returns True is sequence is valid, False otherwise."""
    seqm = dna_seq.upper()
    valid = seqm.count("A") + seqm.count("C") + seqm.count("G") + seqm.count("T")
    return valid == len(seqm)


def frequency(seq: str) -> dict:
    """Calculates the frequency of each symbol in the sequence. Returns a dictionary."""
    dic = {}
    for s in seq.upper():
        dic.setdefault(s, 0)
        dic[s] += 1
    return dic


def gc_content(dna_seq: str) -> float:
    """Returns the percentage of G and C nucleotides in a DNA sequence."""
    gc_count = 0
    for s in dna_seq:
        if s in "GCgc":
            gc_count += 1
    return gc_count / len(dna_seq) * 100


def reverse_complement(dna_seq: str) -> str:
    """Computes the reverse complement of the inputted DNA sequence."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    rev_comp = ""
    for nuc in dna_seq:
        if nuc == "A":
            rev_comp += "T"
        elif nuc == "T":
            rev_comp += "A"
        elif nuc == "C":
            rev_comp += "G"
        elif nuc == "G":
            rev_comp += "C"
    
    return rev_comp


def number_start(dna_seq: str):
    """ Returns the number of start codons in a DNA sequence."""
    start_count = 0
    for i in range(len(dna_seq) - 2):
        if dna_seq[i:i+3] == 'ATG':
            start_count += 1
    return start_count

def number_stop(dna_seq: str):
   """ Returns the number of stop codons in a DNA sequence.""" 
   stop_count = 0
   for i in range(len(dna_seq) - 2):
       if dna_seq[i:i+3] in ['TAA','TAG','TGA']:
           stop_count += 1
   return stop_count

from collections import Counter

def find_codon_freq(dna_seq):
    """Returns a dictionary with codons as keys and their frequencies in the DNA sequence."""
    codon_freq = Counter()
    
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        codon_freq[codon] += 1
    
    return codon_freq

def most_least_freq_codons(dna_seq):
    """Returns the most and least frequent codons in the DNA sequence."""
    codon_freq = find_codon_freq(dna_seq)
    
    most_freq_codon = max(codon_freq, key=codon_freq.get)
    least_freq_codon = min(codon_freq, key=codon_freq.get)
    return most_freq_codon, least_freq_codon

most_least_freq_codons(seq)


# Function to print results
def basic_stats():
    # 1. Show length
    size = len(pos_str)
    print("1. Sequence length:")
    print(f"-> Both the (+) and (-) strands are {size} nucleotides long.")
    print()

    # 2. Show nucleotide frequencies in %
    print("2. Nucleotide Frequencies:")
    freqs = [frequency(pos_str), frequency(neg_str)]
    
    for i in range(2):
        if i == 0:
            print("-> (+) strand:")
        else:
            print("-> (-) strand:")
        # Using lambda notation to sort the frequency dictionary
        ord_freqs = sorted(freqs[i].items(), key = lambda x: x[1], reverse = True)
        for k, v in ord_freqs:
            freq_str = f"-> Nucleotide: {k} | Freq: {v / size * 100:.2f} %"
            print(freq_str)
        print()
    
    # 3. Show GC content
    gc_cont = gc_content(pos_str)
    print("3. GC content:")
    print(f"-> Both the (+) and (-) strands have a GC content of {gc_cont:.2f} %.")
    print()

    # Aesthetic prints for 4 / 5 / 6
    # TODO


### B. Get ORFs
orf_count = 1

def find_ORF(fasta_string):
    global orf_count  # Declare orf_count as global inside the function

    sequence = fasta_string.upper()

    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]

    orfs = {}             # Dictionary to store ORF substrings
    #orfs_coords = {}       # Dictionary to store ORF coordinates

    for frame in range(3):
        start_index = sequence.find(start_codon, frame)
        #print(start_index)
        while start_index != -1:
            # j possible stop_condon
            for j in range(start_index + 3, len(sequence), 3):
                    codon = sequence[j:j+3]
                    if codon in stop_codons:
                        if(len(sequence[start_index:j+3]) >= 2):
                            orfs[f">Protein_ORF{orf_count}"] = [sequence[start_index:j+3], (start_index, j+2)]
                            #orfs[f">Protein_ORF{orf_count}"] = sequence[start_index:j+3]
                            #orfs_coords[f"ORF{orf_count}"] = (start_index, j+2) #2 because want the 0 position
                            orf_count += 1
                            break
                        
            # Find the next start codon
            #print(j)
            start_index = sequence.find(start_codon, j)
            #print(start_index)

    return orfs

def negative_strand(fasta_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_comp_seq = ''.join(complement[base] for base in reversed(fasta_sequence))

    return reverse_comp_seq

def remove_duplicates(orfs1, orfs2):
    combined_orfs = {**orfs1, **orfs2}  # Combine the two dictionaries
    print(combined_orfs)
    unique_orfs = {}  # New dictionary to store unique key-value pairs

    for key, value in combined_orfs.items():
        if value not in unique_orfs.values():  # Check if value is unique
            unique_orfs[key] = value  # Add key-value pair to unique_orfs

    return unique_orfs

def write_fasta(unique_orfs, output_file):
    with open(output_file, 'w') as fasta_file:
        count = 1
        for key, value in unique_orfs.items():
            fasta_file.write(f">Protein_ORF{count}\n")
            fasta_file.write(f"{value[0]}\n")  # Write the first element of the list
            count += 1

    fasta_file.close()
    
def write_coord(unique_orfs, output_file):
    count = 1
    with open(output_file, 'w') as coord_file:
        for key, value in unique_orfs.items():
            coord_file.write(f"{value[1][0]}, {value[1][1]}, ORF{count}\n")
            count += 1
            
    coord_file.close()

### C. Overlap with annotation
def gtf_to_list(filepath: str) -> list:
    """Returns a list of lists where each sublist is a row of a GTF file."""
    with open(filepath, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        gtf_rows = []
        for row in reader:
            gtf_rows.append(row)
    return gtf_rows


def get_exon_coords(gtf_rows: list) -> dict:
    """Returns a dictionary with the coordinates of each exon from a GTF file, by gene_ID.
    Requires the rows in list format (passing your file through "gtf_to_list")."""
    exons = []
    for i in range(len(gtf_rows)):
        if gtf_rows[i][2] == "exon":
            exons.append(gtf_rows[i])
    
    out = {}
    gene_id = ""
    coords = []
    for row in exons:
        raw_gene_id = row[8]
        p1, sep, p2 = raw_gene_id.partition("\"")
        for char in p2:
            if char != "\"":
                gene_id += char
            else:
                break
        
        coords.append(int(row[3]))
        coords.append(int(row[4]))
        out.setdefault(gene_id, tuple(coords))
        gene_id = ""
        coords = []
    return out

def gtf_to_list(filepath: str) -> list:
    """Returns a list of lists where each sublist is a row of a GTF file."""
    dataset = []
    with open(filepath, 'r') as gtf_file:
        gtf_reader = csv.reader(gtf_file, delimiter='\t')
        for row in gtf_reader:
            if not row[0].startswith('#'):  # Skip comment lines
                record = {
                    'seqname': row[0],
                    'source': row[1],
                    'feature': row[2],
                    'start': int(row[3]),
                    'end': int(row[4]),
                    'score': row[5],
                    'strand': row[6],
                    'frame': row[7],
                    'attributes': row[8]
                }
                dataset.append(record)
    return dataset

def get_exon_coords(dataset: list) -> dict:
    """Returns a dictionary with the coordinates of each exon from a GTF file, by gene_ID.
    Requires the rows in list format (passing your file through "gtf_to_list")."""
    exons_dict = {}
    for record in dataset:
        if record['feature'] == 'exon':
            gene_id_index = record['attributes'].find('gene_id')
            gene_id_start_index = gene_id_index + len('gene_id') + 2
            gene_id_end_index = record['attributes'].find('"', gene_id_start_index)
            gene_id = record['attributes'][gene_id_start_index:gene_id_end_index]

            exons_dict[gene_id] = (record['start'], record['end'])
    return exons_dict

def overlap(coord1, coord2):
    start1, end1 = coord1
    start2, end2 = coord2
    
    # Parse coordinates as integers
    start1, end1 = int(start1), int(end1)
    start2, end2 = int(start2), int(end2)

    #print("coord1")
    #print(start1)
    #print(end1)
    #print("coord2")
    #print(start2)
    #print(end2)

    # Calculate intersection
    intersection_start = max(start1, start2)
    intersection_end = min(end1, end2)
    
    # Calculate intersection length
    if(max(start1, start2) <= min(end1, end2)):
        intersection_length = abs(intersection_end - intersection_start)
    else:
        intersection_length = 0
    #print("intersection_length", intersection_length)

    # Min between A and B
    length1 = end1 - start1
    length2 = end2 - start2

    min_length = min(length1, length2)
    #print("min: ", min_length)

    overlap = intersection_length / min_length
    overlap_percentage = overlap * 100
    #print(overlap_percentage)
    
    return overlap_percentage

# What runs when calling the script goes here
if __name__ == "__main__":
    # Process interaction with the arguments
    if len(sys.argv) != 3:
        print("""USAGE: "python yeast_orfs_chr1.py <PATH_TO_FASTA_FILE> <PATH_TO_GTF_FILE>" """)
        print("WARNING: The arguments have to be passed in this order!")
        sys.exit()
    else:
        fasta = sys.argv[1]
        gtf = sys.argv[2]
        pos_str = read_fasta_file(fasta)[:30000]
        neg_str = reverse_complement(pos_str)

    # Outputs from part A
    basic_stats()

    # Ouputs from part B
    fasta_sequence = read_fasta_file(fasta_path)
    #fasta_sequence = "ATGATGTAAATGATGTAA"
    #print(fasta_sequence)
    #orf_positive, orf_postive_coord = find_ORF(fasta_sequence)
    orf_positive = find_ORF(fasta_sequence)
    print(orf_positive)
    #print(orf_postive_coord)
    fasta_sequence_negative = negative_strand(fasta_sequence)
    #orf_negative, orf_negative_coord = find_ORF(fasta_sequence_negative)
    orf_negative = find_ORF(fasta_sequence_negative)
    print(orf_negative)
    #print(orf_negative_coord)
    unique_orfs = remove_duplicates(orf_positive, orf_negative)
    #orfs1 = {'ORF1': 'ATGCTA', 'ORF2': 'TAGCTA'}
    #orfs2 = {'ORF3': 'ATGCTA', 'ORF4': 'TAGCTG'}
    #unique_orfs = remove_duplicates(orfs1, orfs2)
    print("UNIQUE")
    print(unique_orfs)
    write_fasta(unique_orfs, 'all_potential_proteins.fasta')
    write_coord(unique_orfs, 'orf_coordinates.txt')

    # Ouputs from part C
    gtf_data = gtf_to_list(gtf_path)
    exons = get_exon_coords(gtf_data)

    # Extract coordinates from unique_orfs
    coordinates_unique_orfs = [value[1] for value in unique_orfs.values()]

    best_result = -1
    gene_id = ""
    count = 1
    for coord_unique_orfs in coordinates_unique_orfs:
        for key_exons, coord_exons in exons.items():
            #print("Coordinates from exons:")
            overlap_percentage = overlap(coord_unique_orfs, coord_exons)
            if(overlap_percentage > best_result):
                gene_id = key_exons
                best_result = overlap_percentage
            
        print(f"ORF{count} {gene_id}    {best_result}%")
        count += 1
        gene_id = ""
        best_result = -1

