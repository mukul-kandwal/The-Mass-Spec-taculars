from pyteomics import mass
from pyteomics import electrochem 
from pyteomics import parser
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import csv
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objs as go
import plotly as py


def csvToList(csv_file):
    data = []
    with open(csv_file, newline='') as csvfile:
        csv_reader = csv.reader(csvfile)
        for row in csv_reader:
            data.append(row)
    return data 

seq_list = csvToList('output.csv')

def MW(seq_list): #boxplot of molecular weight
    MW_list = []
    columns = ['MW']
    for seq in seq_list: 
        mw = mass.calculate_mass(seq)
        MW_list.append([mw])
    
    df = pd.DataFrame(MW_list, columns = columns) #convert list into dataframe for plotly 
    print(df)
    
    fig = px.box(df, y = 'MW', title='Distribution of molecular weight', points='all')
    fig.show()
#box plot, show all the datapoints right next to the main box plot
MW(seq_list) #you can create if user choose this, execute this 

def pI(seq_list):
    pI_list = []
    columns = ['pI']
    for seq in seq_list:
        for string in seq: 
            parsed_seq = parser.parse(string, show_unmodified_termini=True)
            pi = electrochem.pI(parsed_seq)
            pI_list.append([pi])
    df = pd.DataFrame(pI_list, columns = columns)
    fig = px.histogram(df, x="pI", nbins=20)
    fig.show()

pI(seq_list)

def hydrophobicity(seq_list):
    hydro_list = []
    columns = ['Hydrophobicity']
    for seq in seq_list:
        for string in seq:
            properties = ProteinAnalysis(string)
            hydro_list.append([properties.gravy()])
    df = pd.DataFrame(hydro_list, columns= columns)

    fig = px.histogram(df, x="Hydrophobicity", nbins=20)
    fig.show()

hydrophobicity(seq_list)

