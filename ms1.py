#random function (0-255) 
#nested list for color hex code 
#Input: csv with protein sequence and output: spectrum of m/z values and %intensity: assume 1 for now. Assign a random number for the 
#of isotopes and calculate distance of the peaks from each other using charge. List of hex code (for color) to indicate different proteins in the graph. 
#Print out something that indicate same color = isotopes
#user will the exact charge and calculate the charge for the whole proteins  
#range of random isotopes 
import csv 
import random 
from pyteomics import mass
from colordict import ColorDict
import plotly.graph_objs as go
import matplotlib as plt
#import pymzml


def read(input_file_path):
    # use this function like this: input_list = read(<input_file_path>)
    with open(input_file_path, 'r') as read_obj:
        csv_reader = csv.reader(read_obj)
        list_of_csv = list(csv_reader)
        input_list = []
        for row in list_of_csv:
            input_list.append(row[0])
        return input_list

input_list = read('output.csv')

#colors = ColorDict() 

#put in a dictionary to store all sequence and m/z 
#select 1 m/z and then calculate the isotopes for that 
#p.add for each protein in the dictionary with different color and legend? 
#calculate the m/z for isotopes and add to graph at the same time 
#different intensity of isotope= 4/4
#create dataframe with seq, isotope mz, intensity 

def calc_mz_isotope(input_list, charge):
    mz_isotope_dict = {}
    for seq in input_list: 
        num_isotope = 3
        #iterate and calculate individual isotope m/z value and also intensity value
        tuple_list = []
        for i in range(num_isotope): 
            mz = mass.calculate_mass(sequence=seq, ion_type="M", charge=int(charge))
            isotope = mz + i/int(charge)
            intensity = (num_isotope-i)/num_isotope
            tuple_data = (isotope, intensity)
            tuple_list.append(tuple_data) #add tuple to list of tuple
        if seq not in mz_isotope_dict:
            mz_isotope_dict[seq] = tuple_list
    return mz_isotope_dict


user_input = input("Enter charge of protein: ")
    
mz_dict = calc_mz_isotope(input_list, user_input)

# # Accessing the values
# for protein, mz_values in data.items():
#     print(f"Protein: {protein}")
#     for mz, intensity in mz_values:
#         print(f"m/z: {mz}, Intensity: {intensity}")

#using dictionary for plotly 

def plot_ms(mz_dict):
    # p = pymzml.plot.Factory()
    # p.new_plot()
    traces = []

    for protein, mz_values in mz_dict.items():
        x_values = [item[0] for item in mz_values]
        y_values = [item[1] for item in mz_values]
        if not x_values or not y_values:
            print(f"No data points for category: {protein}")
            continue
        trace = go.Bar(x = x_values, y = y_values)
        traces.append(trace)
    if not traces:
        print("No data to plot")
    else:   
        layout = go.Layout(
        title='m/z Values of Proteins',
        xaxis=dict(title='m/z'),
        yaxis=dict(title='Intensity'))

#Create figure
        fig = go.Figure(data=traces, layout=layout)

    # Plot figure
        fig.show()
plot_ms(mz_dict)

# def plot_ms(input):
#     categories = list(input.keys())
#     years = [pair[0] for pair in input[categories[0]]]  # Assuming all categories have the same years
#     values = {category: [pair[1] for pair in input[category]] for category in categories}

#     # Plotting
#     bar_width = 0.35
#     plt.figure()

#     for i, category in enumerate(categories):
#         plt.subplot(len(categories), 1, i+1)
#         plt.bar([x + i * bar_width for x in years], values[category], bar_width, label=category)
#         plt.xlabel('Year')
#         plt.ylabel('Value')
#         plt.title(category)
#         plt.xticks([x + bar_width / 2 for x in years], years)

#     plt.tight_layout()
#     plt.show()
    
# plot_ms(mz_dict)


#         for i in range(len(peaks)):
#             p.add(peaks[i], color = (51, 255, 246), style = 'sticks', name = "b & y ion peaks")
# # Example dictionary with list of tuples as values
# data = {
#     "protein1": [(1000, 200), (1500, 300), (2000, 400)],
#     "protein2": [(1200, 250), (1600, 350), (2200, 450)],
#     "protein3": [(1100, 220), (1400, 320), (1800, 420)]
# }

# Create traces for each protein
# traces = []
# for protein, mz_values in data.items():
#     mz_values_sorted = sorted(mz_values)  # Sort mz values
#     mz, intensity = zip(*mz_values_sorted)  # Unzip sorted mz values
#     trace = go.Scatter(x=mz, y=intensity, mode='markers', name=protein)
#     traces.append(trace)

# # Create layout
# layout = go.Layout(
#     title='m/z Values of Proteins',
#     xaxis=dict(title='m/z'),
#     yaxis=dict(title='Intensity')
# )

# # Create figure
# fig = go.Figure(data=traces, layout=layout)

# # Plot figure
# fig.show()