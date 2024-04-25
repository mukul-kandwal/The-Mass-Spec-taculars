import csv 
from pyteomics import mass
import plotly.graph_objs as go


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

def plot_ms(mz_dict):
    traces = []
    for protein, mz_values in mz_dict.items():
        x_values = []
        y_values = []
        for mz, i in mz_values:
            x_values.extend([mz - 0.1, mz, mz, mz, mz + 0.1])
            y_values.extend([0, 0, i, 0, 0])
        if not x_values or not y_values:
            print(f"No data points for category: {protein}")
            continue
        trace = go.Scatter(
            {
                "x": x_values,
                "y": y_values,
                "visible": True,
                "mode": "lines",
                "line": {
                    "width": 1,
                },
            }
        )
        traces.append(trace)
    if not traces:
        print("No data to plot")
    else:
        layout = go.Layout(
            title="m/z Values of Proteins",
            xaxis=dict(title="m/z"),
            yaxis=dict(title="Intensity"),
        )
        fig = go.Figure(data=traces, layout=layout)
        fig.show()


plot_ms(mz_dict)
