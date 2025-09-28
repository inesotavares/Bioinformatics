### HEADER

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@   Bioinformatics --- Group Assignment 1  @@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@ By:   -> Bruna Rocha (up201906417 M:CC)          @@@
# @@@       -> Inês Tavares (up201706579 M:DS)         @@@
# @@@       -> Nuno Oliveira (up201703852 M:BBC)       @@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

import sys
import csv

# The genetic code
gen_code = {
        "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
        "TGT":"C", "TGC":"C",
        "GAT":"D", "GAC":"D",
        "GAA":"E", "GAG":"E",
        "TTT":"F", "TTC":"F",
        "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
        "CAT":"H", "CAC":"H",
        "ATA":"I", "ATT":"I", "ATC":"I",
        "AAA":"K", "AAG":"K",
        "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
        "ATG":"M",
        "AAT":"N", "AAC":"N",
        "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
        "CAA":"Q", "CAG":"Q",
        "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
        "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
        "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
        "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
        "TGG":"W",
        "TAT":"Y", "TAC":"Y",
        "TAA":"_", "TAG":"_", "TGA":"_"
    }


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
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_comp = ''.join(complement[base] for base in reversed(dna_seq))
    return reverse_comp


def n_start_codons(dna_seq: str) -> int:
    """Returns the number of start codons in a DNA sequence."""
    start_count = 0
    for i in range(len(dna_seq) - 2):
        if dna_seq[i:i+3] == 'ATG':
            start_count += 1
    return start_count


def n_stop_codons(dna_seq: str) -> int:
   """Returns the number of stop codons in a DNA sequence.""" 
   stop_count = 0
   for i in range(len(dna_seq) - 2):
       if dna_seq[i:i+3] in ['TAA', 'TAG', 'TGA']:
           stop_count += 1
   return stop_count


def find_codon_freq(dna_seq: str) -> dict:
    """Returns a dictionary with codons as keys and their frequencies as values, in a DNA sequence."""
    codon_freqs = {}
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        if codon in gen_code.keys():
            codon_freqs.setdefault(codon, 0)
            codon_freqs[codon] += 1
    return codon_freqs


def most_least_freq_codons(dna_seq: str) -> tuple:
    """Returns a tuple with the following info on a DNA sequence:
    -> The most frequent codons -> list on index 0 / count on index 1;
    -> The least frequent codons -> list on index 2 / count on index 3."""
    codon_freqs = find_codon_freq(dna_seq)
    most_freq_v = max(codon_freqs.values())
    least_freq_v = min(codon_freqs.values())
    most_freq = {}
    least_freq = {}
    for k in codon_freqs.keys():
        if codon_freqs[k] == most_freq_v:
            most_freq.setdefault(k, most_freq_v)
        elif codon_freqs[k] == least_freq_v:
            least_freq.setdefault(k, least_freq_v)
    most_freq_cods = list(most_freq.keys())
    least_freq_cods = list(least_freq.keys())
    return most_freq_cods, most_freq_v, least_freq_cods, least_freq_v


def write_cod_strs(cod_data: tuple) -> tuple:
    """Returns a tuple of strings with a text version of the "most_least_freq_codons" output.
    -> To be used for pretty printing only."""
    # For the (+) strand
    if len(cod_data[0]) > 1:
        cods_pos = "/".join(cod_data[0])
        cods_pos_str = f"-> The most frequent codons are: {cods_pos} --> They appear {cod_data[1]} times."
    else:
        cods_pos = cod_data[0].pop()
        cods_pos_str = f"-> The most frequent codon is: {cods_pos} --> It appears {cod_data[1]} times."
    # For the (-) strand
    if len(cod_data[0]) > 1:
        cods_neg = "/".join(cod_data[2])
        cods_neg_str = f"-> The least frequent codons are: {cods_neg} --> They appear {cod_data[3]} times."
    else:
        cods_neg = cod_data[2].pop()
        cods_neg_str = f"-> The least frequent codon is: {cods_neg} --> It appears {cod_data[3]} times."
    return cods_pos_str, cods_neg_str


### B. Get ORFs

def translate_seq(dna_seq: str) -> str:
    """Translates a DNA sequence to an aminoacid sequence in the one-letter code format.
    -> It follows the genetic code dictionary."""
    assert validate_dna(dna_seq), "Invalid DNA sequence"
    aa_seq = ""
    for i in range(0, len(dna_seq) - 2, 3):
        cod = dna_seq[i:i+3]
        if cod in gen_code.keys():
            aa_seq += gen_code[cod]
    return aa_seq


def orf_finder(pos_str: str, neg_str: str) -> dict:
    """Function to create the output files and return the full coordinates dictionary."""
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = {}   # Dictionary to store ORF substrings

    # For the (+) strand
    for i in range(3):
        while i+3 < len(pos_str):
            codon = pos_str[i:i+3]
            if codon == start_codon:
                for j in range(i+3, len(pos_str), 3):
                    if pos_str[j:j+3] in stop_codons:
                        orfs[pos_str[i:j+3]] = (i, j+2)
                        i = j
                        break
            i += 3

    # For the (-) strand
    for i in range(3):
        while i+3 < len(neg_str):
            codon = neg_str[i:i+3]
            if codon == start_codon:
                for j in range(i+3, len(neg_str), 3):
                    if neg_str[j:j+3] in stop_codons:
                        orfs[neg_str[i:j+3]] = (i, j+2)
                        i = j
                        break
            i += 3
  
    # Verify for the ones with length >= 150
    filtered_dict = {}

    for key, value in orfs.items():
        if len(key) >= 150:
            filtered_dict[key] = value
    
    # 7. Write file "all_potential_proteins.fasta" with the protein sequences
    with open("all_potential_proteins.fasta", "w") as file_protein:
        counter = 1
        for key in filtered_dict.keys():
            file_protein.write(f">Protein_ORF{counter} \n")
            file_protein.write(translate_seq(key) + " \n")
            counter += 1

    # 8. Write file "orf_coordinates.txt" with the coordinates
    coords = {}  # Dictionary for the coordinates
    with open("orf_coordinates.txt", 'w') as coord_file:
        counter = 1
        for value in filtered_dict.values():
            coord_file.write(f"{value[0]}, {value[1]}, ORF{counter} \n")
            coords[f"ORF{counter}"] = value
            counter += 1

    return coords


### C. Overlap with annotation

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


def overlap(coords: dict, annotation: dict):
    """Calculates the overlap of experimental coordinates to annotated ones.
    -> Prints the best scores with labels to stdout."""
    for key_geneid, value_geneid in annotation.items():
        exon_size = value_geneid[1] - value_geneid[0]  # length of A
        best_overlap = -1
        for key, value in coords.items():
            intersec = min(value_geneid[1], value[1]) - max(value_geneid[0], value[0])
            if intersec > 0:
                p_overlap = intersec / exon_size * 100
            else:
                p_overlap = 0
            
            if p_overlap > best_overlap:
                best_overlap = round(p_overlap, 2)
                orf = key
        
        print(f"-> {orf} {key_geneid}   {best_overlap} %")
        best_overlap = -1


if __name__ == "__main__":
    # Process interaction with the arguments
    if len(sys.argv) != 3 or (".fasta" not in sys.argv[1] or ".gtf" not in sys.argv[2]):
        print("USAGE: python yeast_orfs_chr1.py <PATH_TO_FASTA_FILE> <PATH_TO_GTF_FILE>")
        sys.exit()
    else:
        fasta_path = sys.argv[1]
        pos_str = read_fasta_file(fasta_path)[:30000]
        neg_str = reverse_complement(pos_str)
        
        ### Outputs for part A
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
                freq_str = f"-> Nucleotide: {k} | Freq: {v / size * 100:.2f}%"
                print(freq_str)
        print()

        # 3. Show GC content
        gc_cont = gc_content(pos_str)
        print("3. GC content:")
        print(f"-> Both the (+) and (-) strands have a GC content of {gc_cont:.2f} %.")
        print()
    
        # 4. Show number of start codons found
        n_start_cods_pos = n_start_codons(pos_str)
        n_start_cods_neg = n_start_codons(neg_str)
        print("4. Number of start codons:")
        print(f"-> The (+) strand has {n_start_cods_pos} start (AUG) codons.")
        print(f"-> The (-) strand has {n_start_cods_neg} start (AUG) codons.")
        print()

        # 5. Show number of stop codons found
        n_stop_cods_pos = n_stop_codons(pos_str)
        n_stop_cods_neg = n_stop_codons(neg_str)
        print("5. Number of stop codons:")
        print(f"-> The (+) strand has {n_stop_cods_pos} stop (UAA, UAG, and UGA) codons.")
        print(f"-> The (-) strand has {n_stop_cods_neg} stop (UAA, UAG, and UGA) codons.")
        print()

        # 6. Show the most and least frequent codons
        most_least_freq_cods_pos = most_least_freq_codons(pos_str)
        cods_strs_pos = write_cod_strs(most_least_freq_cods_pos)
        most_least_freq_cods_neg = most_least_freq_codons(neg_str)
        cods_strs_neg = write_cod_strs(most_least_freq_cods_neg)
        print("6. Most and least frequent codons:")
        print("-> (+) strand:")
        print(cods_strs_pos[0])
        print(cods_strs_pos[1])
        print("-> (-) strand:")
        print(cods_strs_neg[0])
        print(cods_strs_neg[1])
        print()

        ### Messages for part B
        print("""7. Generating the "all_potential_proteins.fasta" file...""")
        print()
        print("""8. Generating the "orf_coordinates.txt" file...""")
        print()

        coords = orf_finder(pos_str, neg_str)
        gtf_path = sys.argv[2]
        gtf_data = gtf_to_list(gtf_path)
        exons = get_exon_coords(gtf_data)

        ### Outputs for part C
        print("9. Overlaps with annotation:")
        overlap(coords, exons)

