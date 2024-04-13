from Bio import SeqIO
import csv

def fasta_to_csv(filename, outputfile):
    list_sequence = []
    for record in SeqIO.parse(filename, "fasta"):
        sequence = str(record.seq)
        list_sequence.append(sequence)
    with open(outputfile, "w") as csv_out:
        writer = csv.writer(
            csv_out, lineterminator="\n"
        )
        for seq in list_sequence:
            writer.writerow([seq])

inputFasta = "C:\Users\Liv Caspersson\Documents\Proteomics\The-Mass-Spec-taculars\uniprotkb_proteome_UP000000425_2024_01_22.fasta"
outfile_path = "C:\Users\Liv Caspersson\Documents\Proteomics\The-Mass-Spec-taculars/fasta_to_csv output\output.csv"
fasta_to_csv(inputFasta, outfile_path)