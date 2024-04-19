from pyteomics import mass
import plotly.graph_objects as go
import csv 

sequence = "MNYFPIFANLAGRP"

def calc_ions_mz(sequence,charge):
    mz_dict = {'b_ion': [], 'y_ion': []}
    b_ion_ls = []
    y_ion_ls = []
    for i in range(1, len(sequence)):
            b_ion = sequence[:i]
            b_ion_ls.append(mass.calculate_mass(sequence=b_ion, ion_type="b", charge=charge))

            y_ion = sequence[i:]
            y_ion_ls.append(mass.calculate_mass(sequence=y_ion, ion_type="y", charge=charge))

    mz_dict['b_ion']= (b_ion_ls)
    mz_dict['y_ion']=(y_ion_ls)
    return mz_dict, b_ion_ls, y_ion_ls
user_charge = int(input("what charge?"))

mz_dict, b_ion_ls, y_ion_ls = calc_ions_mz(sequence, user_charge)

def plot(mz_dict, b_ion_ls, y_ion_ls):
   traces = []
   for ion, mz_values in mz_dict.items():
        lower_xbound = min(mz_values)
        upper_xbound = max(mz_values)

        x_values = []
        y_values = []

        for mz in mz_values:
            x_values.extend([mz - 0.1, mz, mz, mz, mz + 0.1])
            y_values.extend([0, 0, 1, 0, 0])
            
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

        b_annotations = [dict(
            x=b_ion_ls[i],
            y=1,
            xref="x",
            yref="y",
            text=f"b{i+1}",
            showarrow=True,
            arrowhead=7,
            ax=0,
            ay=-40
        )for i in range(len(b_ion_ls))]

        y_annotations = [dict(
            x=y_ion_ls[i],
            y=1,
            xref="x",
            yref="y",
            text=f"y{i+1}",
            showarrow=True,
            arrowhead=7,
            ax=0,
            ay=-25 #you can change this to change the height of line for the annotation
        )for i in range(len(y_ion_ls))]

        traces.append(trace)

        
        layout = go.Layout(
            title="m/z Values of B and Y ions",
            xaxis=dict(title="m/z", range=[lower_xbound -50 , upper_xbound + 50]),
            yaxis=dict(title="Intensity"),
            annotations= b_annotations + y_annotations,
            showlegend =  False)

        # Create figure
        fig = go.Figure(data=traces, layout=layout)
        # Plot figure
        fig.show()

plot(mz_dict, b_ion_ls, y_ion_ls)