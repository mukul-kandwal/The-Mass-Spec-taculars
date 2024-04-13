from Bio import SeqIO
from pyteomics import parser
import csv

def fasta_to_csv(filename, outputfile):
    list_sequence = []
    for record in SeqIO.parse("UP000000425_2024_01_18.fasta", "fasta"):
        sequence = str(record.seq)
        list_sequence.append(sequence)
    with open(outputfile, "w+") as csv_out:
        writer = csv.DictWriter(
            csv_out, fieldnames= ["Sequence"], lineterminator="\n"
        )
        writer.writeheader()
        for seq in list_sequence:
            writer.write({"Sequence": seq})

inputFasta = "C:\Users\Liv Caspersson\Documents\Proteomics\The-Mass-Spec-taculars\uniprotkb_proteome_UP000000425_2024_01_22.fasta"
outfile_path = "C:\Users\Liv Caspersson\Documents\Proteomics\The-Mass-Spec-taculars\fasta_to_csv output\output.csv"
fasta_to_csv(inputFasta, outfile_path)