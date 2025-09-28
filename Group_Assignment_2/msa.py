from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess

files_name = ["sequence.fasta", "sequence.fasta"]
sequences = []

def fasta_to_txt(fasta_file, txt_file):
    with open(fasta_file, 'r') as fasta:
        with open(txt_file, 'w') as txt:
            for line in fasta:
                if line.startswith('>'):  # Header line
                    txt.write(line.strip() + '\n')
                else:  # Sequence line
                    txt.write(line.strip())

def msa(sequences):
    seq_records = [SeqRecord(Seq(seq), id=f"seq{i+1}") for i, seq in enumerate(sequences)]

    input_file = "teste.fasta"
    SeqIO.write(seq_records, input_file, "fasta")

    # run clustal omega to perform MSA
    output_file = "output.fasta"
    txt_file = "alignment.txt"
    clustalomega_cline = ClustalOmegaCommandline(infile=input_file, outfile=output_file, verbose=True, auto=True)
    #stdout, stderr = clustalomega_cline()
    #
    #if stderr:
    #    print("An error occurred:", stderr)
    #else:
    #    print("Alignment complete. Alignment saved to:", output_file)

    # Clustal Omega via subprocess
    try:
        subprocess.run(str(clustalomega_cline), shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print("An error occurred:", e)

    alignment = AlignIO.read(output_file, "fasta")
    print(alignment)

    fasta_to_txt(output_file, txt_file)

def files():
    sequences = []
    for file in files_name:
        with open(file, "r") as fasta_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                sequences.append(str(record.seq))
    
    print("Sequence data:", sequences)
    #sequences = ["PHWAS", "HWASW", "HPHWA"]
    #sequences = ["ATAGC", "ATGAC", "AACG", "AATCG"]
    return sequences

if __name__ == "__main__":
    sequences = files()
    msa(sequences)

