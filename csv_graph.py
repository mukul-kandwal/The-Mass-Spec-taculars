from Bio import SeqIO
import csv
import pandas as pd 
import numpy as np
import plotly.express as px
import plotly.graph_objs as go
import plotly as py


def csvTodf(csv_file):
    return pd.read_csv(csv_file)


df = csvTodf('protein_properties.csv')

def MW(dataframe):
    #plot molecular weight results on boxplot 
    fig = px.box(dataframe, y = 'MW', title='Distribution of molecular weight of proteins within proteome')
    fig.show()

#MW(df) #you can create if user choose this, execute this 

def pI(dataframe):
    # create the bins
    counts, bins = np.histogram(dataframe.pI, bins=range(0, 13, 1))
    bins = 0.5 * (bins[:-1] + bins[1:])
    fig = px.bar(x=bins, y=counts, labels={'x':'isoelectric point (pI)', 'y':'Total number of proteins'})
    fig.show()
    # fig = px.histogram(dataframe, x="pI", nbins=20)
    # fig.show()

#pI(df)

def hydrophobicity(dataframe):
    fig = px.histogram(dataframe, x="Hydrophobicity", nbins=20)
    fig.show()

#hydrophobicity(df)

def type(dataframe):
    fig = px.histogram(dataframe, x= "Category")
    fig.show()

type(df)