from Bio import SeqIO
from pyteomics import mass
from pyteomics import electrochem
from pyteomics import parser
import csv
import pathlib
import PySimpleGUI as sg
sg.theme('Dark Purple 2')


def protein_prop(filepath):
    protein_id = []
    mol_weight_list = []  # list of molecular weights of all proteins
    pi_list = []  # list of protein pI
    charge_list = []  # list of charge
    gravy_list = []  # list of hydrophobicity values

    for record in SeqIO.parse(filepath, "fasta"):
        proteinID = str(record.id)
        protein_id.append(proteinID)
        sequence = str(record.seq)
        # calculate molecular weight
        mw = mass.calculate_mass(sequence)
        mol_weight_list.append(mw)
        parsed_seq = parser.parse(record.seq, show_unmodified_termini=True)
        # calculate pI of protein and add to list
        pi = electrochem.pI(parsed_seq)
        pi_list.append(pi)
        # calculate charge and add to list
        charge = electrochem.charge(parsed_seq, 7)
        charge_list.append(charge)

    # writing to a csv file
    path = pathlib.Path(filepath)
    writepath = str(path.parent) + "/protein_prop.csv"
    with open(writepath, "w") as csv_out:
        # How the csv looks like
        writer = csv.DictWriter(
            csv_out, fieldnames=["Protein ID", "Molecular Weight", "pI", "Charge at pH 7"], lineterminator="\n"
        )
        writer.writeheader()
        # Iterating over the protein entries
        for n, prot in enumerate(protein_id):
            writer.writerow({"Protein ID": prot, "Molecular Weight": mol_weight_list[n], "pI": pi_list[n],
                             "Charge at pH 7": charge_list[n]})


layout = [[sg.Text("Choose Proteome File:"), sg.FileBrowse(key='-filepath-')], [sg.Button("Go!")]]
window = sg.Window("Calculate Protein Properties", layout)

while True:
    event, values = window.read()
    # End program if user closes window or
    # presses the OK button
    if event == sg.WIN_CLOSED:
        exit()
    if event == "Go!":
        protein_prop(values["-filepath-"])
        break

window.close()


