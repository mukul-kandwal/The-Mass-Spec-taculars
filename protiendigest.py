import csv
from pyteomics import parser


def writer(list, output_file_path):
    # Call this function like this: writer(<the list you want to save>, <the file path>)
    with open(output_file_path, "w") as csv_out:
        writer = csv.writer(
            csv_out, lineterminator="\n"
        )
        for seq in list:
            writer.writerow([seq])

def digest(CSV_file, method, outputfile):
    list_peptides = []
    dict_peptides = dict()
    set_peptides = set()
    with open(CSV_file, 'r') as read_obj:
        # Return a reader object which will
        # iterate over lines in the given csvfile
        csv_reader = csv.reader(read_obj)

        # convert string to list
        list_of_csv = list(csv_reader)
        input_list = []
        for row in list_of_csv:
            input_list.append(row[0])
        for seq in input_list:
            digested_peptides = parser.cleave(seq, method)
            for pep in digested_peptides:
                list_peptides.append(pep)
                set_peptides.add(pep)
                if pep in dict_peptides.keys():
                    dict_peptides[pep] += 1
                else:
                    dict_peptides[pep] = 1

    writer(list_peptides, outputfile)



outputfile = 'C:/Users\mukul\Desktop\Spring24\Proteomics\Group Project\output2.csv'
inputfile = 'C:/Users\mukul\Desktop\Spring24\Proteomics\Group Project\output.csv'
method = 'trypsin'
digest(inputfile, method, outputfile)