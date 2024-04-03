from Bio import SeqIO
from pyteomics import parser

def digest_peptides(filename):
    list_peptides = []
    set_peptides = set()
    dict_peptides = dict()
    for record in SeqIO.parse("UP000000425_2024_01_18.fasta", "fasta"):
        protein_id = str(record.id)
        sequence = str(record.seq)
    # use pyteomics to generate digested peptides
        digested_peptides = parser.cleave(filename, digest)
        for pep in digested_peptides:
            list_peptides.append(pep)
            set_peptides.add(pep)
            if pep in dict_peptides.keys():
                dict_peptides[pep] += 1
            else:
                dict_peptides[pep] = 1

print("number of peptides overall:", len(list_peptides))
print("number of peptide sequences", len(set_peptides))
print("number of peptide sequences in dict", len(dict_peptides.keys()))

counter = 0
for pep in dict_peptides.keys():
    if dict_peptides[pep] > 1:
        counter += 1

print(counter)