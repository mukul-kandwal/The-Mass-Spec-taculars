import csv
from pyteomics import mass
from pyteomics import electrochem 
from pyteomics import parser
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import plotly as py

def read(input_file_path):
    with open(input_file_path, 'r') as read_obj:
        csv_reader = csv.reader(read_obj)
        list_of_csv = list(csv_reader)
        input_list = []
        for row in list_of_csv:
            input_list.append(row[0])
        return input_list

input_list = read('output.csv')

def calc_pI(input_list):
    pI_list = []
    for seq in input_list:
        for string in seq: 
            parsed_seq = parser.parse(string, show_unmodified_termini=True)
            pi = electrochem.pI(parsed_seq)
            pI_list.append(pi)
    return pI_list

pI_list = calc_pI(input_list)

def plot_isoelectric(pI_list):
    labdict = {"value": "pI", 'variable' : ""}
    fig = px.strip(pd.DataFrame(pI_list),orientation='h', range_x=(0,14), labels = labdict)
    fig.show()

plot_isoelectric(pI_list)
