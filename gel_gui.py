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

def calc_MW(input_list):
    MW_list = []
    for seq in input_list: 
        mw = mass.calculate_mass(seq)
        MW_list.append(mw)
    return MW_list

MW_list = calc_MW(input_list)

def plot_gel(MW_list):
    labdict = {"value": "Molecular Weight", 'variable' : ""}
    fig = px.strip(pd.DataFrame(MW_list), labels = labdict)
    fig.show()

plot_gel(MW_list)
