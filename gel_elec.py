import csv
from pyteomics import mass
from pyteomics import electrochem 
from pyteomics import parser
import plotly.graph_objects as go

def read(input_file_path):
    with open(input_file_path, 'r') as read_obj:
        csv_reader = csv.reader(read_obj)
        list_of_csv = list(csv_reader)
        input_list = []
        for row in list_of_csv:
            input_list.append(row[0])
        return input_list

input_list = read('output.csv')

#isoelectric point 
def calc_pI(input_list):
    pI_list = []
    for seq in input_list:
        for string in seq: 
            parsed_seq = parser.parse(string, show_unmodified_termini=True)
            pi = electrochem.pI(parsed_seq)
            pI_list.append(pi)
    return pI_list

pI_list = calc_pI(input_list)

def isoelectric(pI_list):
    y_axis = ['string' + str(i) for i in range(len(pI_list))] #create an arbitrary string as y-axis to simulate the graph
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=pI_list, y=y_axis, mode='markers', name='pI Values'))
    fig.update_layout(title='Isoelectric Point (pI) of sequences',
                  xaxis_title='pH',
                  yaxis_visible=False, 
                  yaxis_showticklabels=False)
    fig.show()

isoelectric(pI_list)

def calc_MW(input_list):
    MW_list = []
    for seq in input_list: 
        mw = mass.calculate_mass(seq)
        MW_list.append(mw)
    return MW_list

MW_list = calc_MW(input_list)

def gel(MW_list):
    y_axis = ['string' + str(i) for i in range(len(MW_list))] #create an arbitrary string as y-axis to simulate the graph 
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=y_axis, y=MW_list, mode='markers', name='MW Values'))
    fig.update_layout(title='Molecular weight of sequences',
                  yaxis_title='Molecular Weight',
                  xaxis_visible=False, 
                  xaxis_showticklabels=False)
    fig.show()

gel(MW_list)


