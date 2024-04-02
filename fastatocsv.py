from Bio import SeqIO
from pyteomics import mass
from pyteomics import electrochem 
from pyteomics import parser
import csv



def protein_prop(filename):
    protein_id = [] 
    mol_weight_list = [] #list of molecular weights of all proteins 
    pi_list = [] #list of protein pI
    charge_list = [] #list of charge
    gravy_list = [] #list of hydrophobicity values

    for record in SeqIO.parse(filename, "fasta"):
        proteinID = str(record.id)
        protein_id.append(proteinID)
        sequence = str(record.seq)
        #calculate molecular weight
        mw = mass.calculate_mass(sequence)
        mol_weight_list.append(mw)
        parsed_seq = parser.parse(record.seq, show_unmodified_termini=True)
        #calculate pI of protein and add to list
        pi = electrochem.pI(parsed_seq)
        pi_list.append(pi)
        #calculate charge and add to list
        charge = electrochem.charge(parsed_seq, 7)
        charge_list.append(charge)
    
    # #digesting protein
    # tryptic_peptides = parser.cleave(sequence, digest) #return set of peptide digested by trypsin for each protein
    # #create loop to add each digested peptide into dictionary 
    # for pep in tryptic_peptides:
    #     peptides.add(pep)
    #     if pep in dict_peptides.keys():
    #         dict_peptides[pep] += 1
    #     else:
    #         dict_peptides[pep] = 1
    #     #calculate retention time and add to list 
    #     retention_ls.append(achrom.calculate_RT(pep, achrom.RCs_krokhin_100A_tfa)) 
		#writing to a csv file
    with open("protein_prop.csv", "w") as csv_out:
# How the csv looks like
        writer = csv.DictWriter(
            csv_out, fieldnames=["Protein ID", "Molecular Weight", "pI", "Charge"], lineterminator="\n"
        )
        writer.writeheader()
    # Iterating over the protein entries
        for n, prot in enumerate(protein_id):
            writer.writerow({"Protein ID": prot, "Molecular Weight": mol_weight_list[n], "pI": pi_list[n] , "Charge": charge_list[n]})

protein_prop("uniprotkb_proteome_UP000036027_2024_04_01.fasta")