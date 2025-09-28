import sys 
import csv

orf_count = 1

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

def find_ORF(fasta_string):
    global orf_count  # Declare orf_count as global inside the function

    sequence = fasta_string.upper()

    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]

    orfs = {}             # Dictionary to store ORF substrings
    #orfs_coords = {}       # Dictionary to store ORF coordinates
    
    i = 0
    print(len(sequence))
    while(i+3 < len(sequence)):
        print("value of i", i)
        codon = sequence[i:i+3]
        if(codon == start_codon):
            for j in range(i+3, len(sequence), 3):
                print(j)
                print(sequence[j:j+3])
                if sequence[j:j+3] in stop_codons:
                    if(len(sequence[i:j+3]) >= 2):
                        orfs[sequence[i:j+3]] = (i, j+2)
                        i = j
                        break
        i += 3
        
    return orfs

def negative_strand(fasta_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_comp_seq = ''.join(complement[base] for base in reversed(fasta_sequence))

    return reverse_comp_seq

def remove_duplicates(orfs1, orfs2):
    combined_orfs = {**orfs1, **orfs2}  # Combine the two dictionaries
    #print(combined_orfs)
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

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <sequence_chr1.fasta_path> <genes_chr1.gtf_path>")
    else:
        fasta_path = sys.argv[1]
        gtf_path = sys.argv[2]

    fasta_sequence = read_fasta_file(fasta_path)
    #fasta_sequence = "ATGATGTAAATGATGTAA"
    #print(fasta_sequence)
    #orf_positive, orf_postive_coord = find_ORF(fasta_sequence)
    orf_positive = find_ORF(fasta_sequence)
    #print(orf_positive)
    #print(orf_postive_coord)
    fasta_sequence_negative = negative_strand(fasta_sequence)
    #orf_negative, orf_negative_coord = find_ORF(fasta_sequence_negative)
    orf_negative = find_ORF(fasta_sequence_negative)
    #print(orf_negative)
    #print(orf_negative_coord)
    unique_orfs = remove_duplicates(orf_positive, orf_negative)
    #orfs1 = {'ORF1': 'ATGCTA', 'ORF2': 'TAGCTA'}
    #orfs2 = {'ORF3': 'ATGCTA', 'ORF4': 'TAGCTG'}
    #unique_orfs = remove_duplicates(orfs1, orfs2)
    #print("UNIQUE")
    #print(unique_orfs)
    write_fasta(unique_orfs, 'all_potential_proteins.fasta')
    write_coord(unique_orfs, 'orf_coordinates.txt')
    #unique_coords

    #print(orf_negative)
    
    gtf_data = gtf_to_list(gtf_path)
    exons = get_exon_coords(gtf_data)
    
    #overlap(unique_orfs, exons)

    #unique_orfs = {'key1': ['value1', ('2', '20')], 'key2': ['value2', ('15', '25')]}
    #exons = {'key1': ('10', '20'), 'key2': ('15', '25')}
    
    # Extract coordinates from unique_orfs and exons
    coordinates_unique_orfs = [value[1] for value in unique_orfs.values()]
    #coordinates_exons = list(exons.values())

    #print("Coordinates from unique_orfs:", coordinates_unique_orfs)
    #print("Coordinates from exons:", coordinates_exons)
    
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
