# from pyteomics import parser
# from pyteomics import electrochem
from Bio import SeqIO
import csv
import PySimpleGUI as sg
import os.path
from pyteomics import mass
from pyteomics import achrom
import pyteomics
import pandas as pd
import plotly.express as px
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import plotly.graph_objs as go
sg.theme('Dark Brown 5')

# General purpose functions ------------------------------------------------------------------------
def read(input_file_path):
    # use this function like this: input_list = read(<input_file_path>)
    with open(input_file_path, 'r') as read_obj:
        csv_reader = csv.reader(read_obj)
        list_of_csv = list(csv_reader)
        input_list = []
        for row in list_of_csv:
            input_list.append(row[0])
        return input_list
def writer(list_to_write, output_file_path):
    # Call this function like this: writer(<the list you want to save>, <the file path>)
    with open(output_file_path, "w") as csv_out:
        writer = csv.writer(
            csv_out, lineterminator="\n"
        )
        for seq in list_to_write:
            writer.writerow([seq])

def mainscreen():
    layout1 = [[sg.Text("General Purpose Functions")],
               [sg.Button("Convert Fasta to CSV")],
               [sg.Button("Enzymatic Digest")],
               [sg.Button("Graph Peptide Properties")]]

    layout2 = [[sg.Text("Protein Separation Methods")],
               [sg.Button("Isoelectric Focusing")],
               [sg.Button("Gel Electrophoresis")],
               [sg.Button("Chromatography")]]

    layout3 = [[sg.Text("Simulate Mass Spec!")],
               [sg.Button("MALDI-TOF")],
               [sg.Button("Electrospray-Orbitrap")],
               [sg.HSeparator()],
               [sg.Button("Exit")]]

    layout = [[sg.Column(layout1), sg.VSeperator(), sg.Column(layout2), sg.VSeperator(), sg.Column(layout3)]]

    window = sg.Window("Mass-Spectacular 3000", layout)
    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED or event == "Exit":
            exit()
        if event == "Convert Fasta to CSV":
            window.close()
            InitalRead()
        if event == "Enzymatic Digest":
            window.close()
            Digest()
        if event == "Graph Peptide Properties":
            window.close()
            Properties()
        if event == "Isoelectric Focusing":
            window.close()
            Isoelectric()
        if event == "Gel Electrophoresis":
            window.close()
            Electrophoresis()
        if event == "Chromatography":
            window.close()
            Chromatography()
        if event == "MALDI-TOF":
            window.close()
            MALDISpec()
        if event == "Electrospray-Orbitrap":
            window.close()
            OrbiSpec()

# functions for specific jobs ------------------------------------------------------------------------
def InitalRead():
    layout1 = [[sg.Text("Choose Fasta File:"), sg.FileBrowse(key='-FASTA-')],
               [sg.Text("Choose Output Folder:"), sg.FolderBrowse(key='-Folder-')],
               [sg.Button("<--"), sg.Button("Convert!")]]
    window = sg.Window("Convert FASTA", layout1)
    while True:
        event, values = window.read()
        exception = False
        if event == sg.WIN_CLOSED or event == "<--":
            window.close()
            mainscreen()
            break
        if event == "Convert!":
            file = values["-FASTA-"]
            folder = values["-Folder-"]
            if os.path.basename(file)[-6:] != ".fasta":
                sg.popup("Input file was not in FASTA format")
                event = ""
            if folder == "":
                sg.popup("No Output Folder was Selected")
                event = ""
            else:
                filename = (os.path.basename(file)[:-6] + ".csv")
                inputFasta = file
                outfile_path = os.path.join(folder, filename)
                list_sequence = []
                for record in SeqIO.parse(inputFasta, "fasta"):
                    sequence = str(record.seq)
                    list_sequence.append(sequence)
                writer(list_sequence, outfile_path)

                sg.popup("Finished")
    window.close()
def Digest():
    Lysis_Methods = ['arg-c', 'asp-n', 'bnps-skatole', 'caspase 1', 'caspase 2', 'caspase 3', 'caspase 4', 'caspase 5', 'caspase 6', 'caspase 7', 'caspase 8', 'caspase 9', 'caspase 10', 'chymotrypsin high specificity', 'chymotrypsin low specificity', 'clostripain', 'cnbr', 'enterokinase', 'factor xa', 'formic acid', 'glutamyl endopeptidase', 'granzyme b', 'hydroxylamine', 'iodosobenzoic acid', 'lysc', 'ntcb', 'pepsin ph1.3', 'pepsin ph2.0', 'proline endopeptidase', 'proteinase k', 'staphylococcal peptidase i', 'thermolysin', 'thrombin', 'trypsin', 'trypsin_exception']
    layout1 = [[sg.Text("Choose CSV File:"), sg.FileBrowse(key="-CSV-")],
               [sg.Text("Enzyme used for digestion:"), sg.Combo(Lysis_Methods, font=('Arial Bold', 8), expand_x=True, enable_events=True, readonly=False, key='-COMBO_Enzymes-')],
               [sg.Text("Choose Output Folder:"), sg.FolderBrowse(key='-FOLDER-')],
               [sg.Button("<--"), sg.Button("Digest!")]]
    window = sg.Window("Enzymatic Digest", layout1)
    while True:
        event, values = window.read()
        exception = False
        if event == sg.WIN_CLOSED or event == "<--":
            window.close()
            mainscreen()
            break
        if event == "Digest!":
            file = values["-CSV-"]
            folder = values["-FOLDER-"]
            if os.path.basename(file)[-4:] != ".csv":
                sg.popup("Input file was not in CSV format")
                event = ""
            elif values['-COMBO_Enzymes-'] not in Lysis_Methods:
                sg.popup("Incorrect Enzyme was Selected")
                event = ""
            elif folder == "":
                sg.popup("No Output Folder was Selected")
                event = ""
            else:
                filename = (os.path.basename(file)[:-4] + "_Digested.csv")
                outfile_path = os.path.join(folder, filename)
                input_list = read(file)
                list_peptides = []
                dict_peptides = dict()
                set_peptides = set()
                for seq in input_list:
                    digested_peptides = pyteomics.parser.cleave(seq, values["-COMBO_Enzymes-"])
                    for pep in digested_peptides:
                        list_peptides.append(pep)
                        set_peptides.add(pep)
                        if pep in dict_peptides.keys():
                            dict_peptides[pep] += 1
                        else:
                            dict_peptides[pep] = 1

                writer(list_peptides, outfile_path)
                sg.popup("I'm Full! (Digestion Finished)", "Total number of digested peptides: " + str(len(list_peptides)), "Number of prototypic peptides:    " + str(len(dict_peptides)))

    window.close()
def Properties():
    layout1 = [[sg.Text("Choose CSV File:"), sg.FileBrowse(key="-CSV-")],
               [sg.Checkbox("Molecular Weight Distribution", key='MW'), sg.Checkbox("Isoelectric Point Distribution", key='pI'), sg.Checkbox("Hydrophobicity Distribution", key='Hydro')],
               [sg.Button("<--"), sg.Button("Make Graphs!")]]
    window = sg.Window("Graph Peptide Properties", layout1)
    while True:
        event, values = window.read()
        exception = False
        if event == sg.WIN_CLOSED or event == "<--":
            window.close()
            mainscreen()
            break
        if event == "Make Graphs!":
            file = values["-CSV-"]
            if os.path.basename(file)[-4:] != ".csv":
                sg.popup("Input file was not in CSV format")
                event = ""
            else:
                input_list = read(file)
                if values["MW"]:
                    MW_list = []
                    columns = ['MW']
                    for seq in input_list:
                        mw = mass.calculate_mass(seq)
                        MW_list.append([mw])
                    df = pd.DataFrame(MW_list, columns=columns)  # convert list into dataframe for plotly
                    fig = px.box(df, y='MW', title='Distribution of molecular weight', points='all')
                    fig.show()
                if values["pI"]:
                    pI_list = []
                    columns = ['pI']
                    for seq in input_list:
                        for string in seq:
                            parsed_seq = pyteomics.parser.parse(string, show_unmodified_termini=True)
                            pi = pyteomics.electrochem.pI(parsed_seq)
                            pI_list.append([pi])
                    df = pd.DataFrame(pI_list, columns=columns)
                    fig = px.histogram(df, x="pI", nbins=20)
                    fig.show()
                if values["Hydro"]:
                    hydro_list = []
                    columns = ['Hydrophobicity']
                    for seq in input_list:
                        for string in seq:
                            properties = ProteinAnalysis(string)
                            hydro_list.append([properties.gravy()])
                    df = pd.DataFrame(hydro_list, columns=columns)
                    fig = px.histogram(df, x="Hydrophobicity", nbins=20)
                    fig.show()
                sg.Popup("Graphing Complete!")

    window.close()
def Isoelectric():
    layout = [[sg.Text("Choose CSV File:"), sg.FileBrowse(key="-CSV-")],
              [sg.Button("Show Isoelectric focusing gel")],
              [sg.HSeparator()],
              [sg.Text("Choose Output Folder:"), sg.FolderBrowse(key='-FOLDER-')],
              [sg.Text("Collect Peptides between pH: "), sg.Input(key='-Low-', s=4), sg.Text("and pH: "), sg.Input(key='-High-', s=4)],
              [sg.Button("<--"), sg.Button("Collect Peptides")]]
    window = sg.Window("Isoelectric Focusing", layout)
    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED or event == "<--":
            window.close()
            mainscreen()
            break
        if event == "Show Isoelectric focusing gel":
            file = values["-CSV-"]
            if os.path.basename(file)[-4:] != ".csv":
                sg.popup("Input file was not in CSV format")
                event = ""
            else:
                input_list = read(file)
                pI_list = []
                for seq in input_list:
                    parsed_seq = pyteomics.parser.parse(seq, show_unmodified_termini=True)
                    pi = pyteomics.electrochem.pI(parsed_seq)
                    pI_list.append(pi)
                # Graph the gel
                # y_axis = ['string' + str(i) for i in range(len(pI_list))]
                # fig = go.Figure()
                # fig.add_trace(go.Scatter(x=pI_list, y=y_axis, mode='markers', name='pI Values'))
                # fig.update_layout(title='Isoelectric Point (pI) of sequences',
                #                   xaxis_title='pH',
                #                   yaxis_visible=False,
                #                   yaxis_showticklabels=False)
                # fig.show()
                labdict = {"value": "Molecular Weight", "variable": ""}
                fig = px.strip(pd.DataFrame(pI_list), orientation='h', range_x=(0, 14), labels=labdict)
                fig.show()
        if event == "Collect Peptides":
            file = values["-CSV-"]
            filename = (os.path.basename(file)[:-4] + "_IsoelectricFocused.csv")
            folder = values["-FOLDER-"]
            input_list = read(file)
            if os.path.basename(file)[-4:] != ".csv":
                sg.popup("Input file was not in CSV format")
                event = ""
            elif folder == "":
                sg.popup("Output folder was not given")
                event = ""
            else:
                outfile_path = os.path.join(folder, filename)
                try:
                    int(values['-Low-'])
                except:
                    sg.popup("Lower bound of pH for Peptide collection was not valid")
                    event = ""
                else:
                    try:
                        int(values['-High-'])
                    except:
                        sg.popup("Upper bound of pH for Peptide collection was not valid")
                        event = ""
                    else:
                        if int(values['-High-']) < int(values['-Low-']):
                            sg.Popup("Upper bound of pH is smaller than Lower bound")
                            event = ""
                        else:
                            pI_list = []
                            for seq in input_list:
                                parsed_seq = pyteomics.parser.parse(seq, show_unmodified_termini=True)
                                pi = pyteomics.electrochem.pI(parsed_seq)
                                pI_list.append(pi)

                            user_seq = []
                            for i in range(len(pI_list)):
                                if int(values['-Low-']) <= pI_list[i] <= int(values['-High-']):
                                    user_seq.append(input_list[i])

                            writer(user_seq, outfile_path)
                            sg.popup("Done!")
def Electrophoresis():
    layout = [[sg.Text("Choose CSV File:"), sg.FileBrowse(key="-CSV-")],
              [sg.Button("Show Electrophoresis gel")],
              [sg.HSeparator()],
              [sg.Text("Choose Output Folder:"), sg.FolderBrowse(key='-FOLDER-')],
              [sg.Text("Collect Peptides between molecular weight: "), sg.Input(key='-Low-', s=10), sg.Text("and: "), sg.Input(key='-High-', s=10)],
              [sg.Button("<--"), sg.Button("Collect Peptides")]]
    window = sg.Window("Gel Electrophoresis", layout)
    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED or event == "<--":
            window.close()
            mainscreen()
            break
        if event == "Show Electrophoresis gel":
            file = values["-CSV-"]
            if os.path.basename(file)[-4:] != ".csv":
                sg.popup("Input file was not in CSV format")
                event = ""
            else:
                input_list = read(file)
                MW_list = []
                for seq in input_list:
                    mw = mass.calculate_mass(seq)
                    MW_list.append(mw)
                # Graph the gel
                labdict = {"value": "Molecular Weight", "variable": ""}
                fig = px.strip(pd.DataFrame(MW_list), labels=labdict)
                fig.show()
        if event == "Collect Peptides":
            file = values["-CSV-"]
            filename = (os.path.basename(file)[:-4] + "_ElectrophoresisSeperated.csv")
            folder = values["-FOLDER-"]
            input_list = read(file)
            if os.path.basename(file)[-4:] != ".csv":
                sg.popup("Input file was not in CSV format")
                event = ""
            elif folder == "":
                sg.popup("Output folder was not given")
                event = ""
            else:
                outfile_path = os.path.join(folder, filename)
                try:
                    int(values['-Low-'])
                except:
                    sg.popup("Lower bound of molecular weight for Peptide collection was not valid")
                    event = ""
                else:
                    try:
                        int(values['-High-'])
                    except:
                        sg.popup("Upper bound of molecular weight for Peptide collection was not valid")
                        event = ""
                    else:
                        if int(values['-High-']) < int(values['-Low-']):
                            sg.Popup("Upper bound of molecular weight is smaller than Lower bound")
                            event = ""
                        else:
                            MW_list = []
                            for seq in input_list:
                                mw = mass.calculate_mass(seq)
                                MW_list.append(mw)

                            user_seq = []
                            for i in range(len(MW_list)):
                                if int(values['-Low-']) <= MW_list[i] <= int(values['-High-']):
                                    user_seq.append(input_list[i])

                            writer(user_seq, outfile_path)
                            sg.popup("Done!")
def Chromatography():
    methodsTranslation = {'Reverse Phase at pH 2.0': achrom.RCs_guo_ph2_0,
               'Reverse Phase at pH 7.0': achrom.RCs_guo_ph7_0,
               # '': achrom.RCs_meek_ph2_1,
               # '': achrom.RCs_meek_ph7_4,
               # '': achrom.RCs_browne_tfa,
               # '': achrom.RCs_browne_hfba,
               # '': achrom.RCs_palmblad,
               'Normal phase chromatography': achrom.RCs_yoshida,
               # '': achrom.RCs_yoshida_lc,
               # '': achrom.RCs_zubarev,
               'Hydrophilic interaction liquid chromatography at pH 3.0': achrom.RCs_gilar_atlantis_ph3_0,
               'Hydrophilic interaction liquid chromatography at pH 4.5': achrom.RCs_gilar_atlantis_ph4_5,
               'Hydrophilic interaction liquid chromatography at pH 10.0': achrom.RCs_gilar_atlantis_ph10_0,
               # '': achrom.RCs_gilar_beh,
               # '': achrom.RCs_gilar_beh_amide,
               # '': achrom.RCs_gilar_rp,
               # '': achrom.RCs_krokhin_100A_fa,
               # '': achrom.RCs_krokhin_100A_tfa,
               }
    methods = ['Reverse Phase at pH 2.0', 'Reverse Phase at pH 7.0', 'Normal phase chromatography', 'Hydrophilic interaction liquid chromatography at pH 3.0', 'Hydrophilic interaction liquid chromatography at pH 4.5',  'Hydrophilic interaction liquid chromatography at pH 10.0']
    layout = [[sg.Text("Choose CSV File:"), sg.FileBrowse(key="-CSV-")],
              [sg.Text("Type of Chromatography:"), sg.Combo(methods, font=('Arial Bold', 8), expand_x=True, enable_events=True, readonly=False, key='-methods-')],
              [sg.Text("Chromatography run time"), sg.Input(key='-RunTime-', s=4), sg.Button("Show Elution Profile")],
              [sg.HSeparator()],
              [sg.Text("Choose Output Folder:"), sg.FolderBrowse(key='-FOLDER-')],
              [sg.Text("Fraction Collection Time Start: "), sg.Input(key='-Low-', s=4)],
              [sg.Text("Fraction Collection Time End: "), sg.Input(key='-High-', s=4)],
              [sg.Button("<--"), sg.Button("Collect Fraction")]]
    window = sg.Window("Chromatography", layout)
    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED or event == "<--":
            window.close()
            mainscreen()
            break
        if event == "Show Elution Profile":
            # code exeption handling
            file = values["-CSV-"]
            if os.path.basename(file)[-4:] != ".csv":
                sg.popup("Input file was not in CSV format")
                event = ""
            elif values['-methods-'] not in methods:
                sg.popup("Chromatography method given was invalid")
                event = ""
            else:
                try:
                    int(values['-RunTime-'])
                except:
                    sg.popup("Chromatography run time given was not an integer")
                    event = ""
                else:
                    input_list = read(file)
                    retention_ls = []
                    for seq in input_list:
                        retention_ls.append(pyteomics.achrom.calculate_RT(seq, methodsTranslation[values['-methods-']]))
                    normalized_time = []
                    for time in retention_ls:
                        ((time - min(retention_ls)) / (max(retention_ls) - min(retention_ls)))
                        normalized_time.append(((time - min(retention_ls)) / (max(retention_ls) - min(retention_ls))) * int(values['-RunTime-']))
                    fig2 = go.Figure()
                    fig2.add_trace(go.Histogram(x=normalized_time, marker=dict(color=[0, 255, 255])))
                    fig2.update_xaxes(title_text="Retention Time (min)")
                    fig2.update_yaxes(title_text="Number of sequences")
                    fig2.show()
        if event == "Collect Fraction":
            file = values["-CSV-"]
            folder = values["-FOLDER-"]
            if os.path.basename(file)[-4:] != ".csv":
                sg.popup("Input file was not in CSV format")
                event = ""
            elif folder == "":
                sg.popup("Output folder was not specified")
                event = ""
            elif values['-methods-'] not in methods:
                sg.popup("Chromatography method given was invalid")
                event = ""
            else:
                filename = (os.path.basename(file)[:-4] + "_Fraction.csv")
                outfile_path = os.path.join(folder, filename)
                input_list = read(file)
                try:
                    int(values['-RunTime-'])
                except:
                    sg.popup("Chromatography run time given was not an integer")
                    event = ""
                else:
                    try:
                        int(values['-Low-'])
                    except:
                        sg.popup("Fraction collection start time given was not an integer")
                        event = ""
                    else:
                        try:
                            int(values['-High-'])
                        except:
                            sg.popup("Fraction collection end time given was not an integer")
                            event = ""
                        else:
                            if int(values['-High-']) < int(values['-Low-']):
                                sg.Popup("Fraction start time is after Fraction end time!")
                                event = ""
                            else:
                                retention_ls = []
                                for seq in input_list:
                                    retention_ls.append(pyteomics.achrom.calculate_RT(seq, methodsTranslation[values['-methods-']]))
                                normalized_time = []
                                for time in retention_ls:
                                    normalized_time.append((((time / max(retention_ls)) - min(retention_ls)) * int(values['-RunTime-'])))
                                user_seq = []
                                for i in range(len(normalized_time)):
                                    if int(values['-Low-']) <= normalized_time[i] <= int(values['-High-']):
                                        user_seq.append(input_list[i])

                                writer(user_seq, outfile_path)
                                sg.popup("Done!")
def MALDISpec():
    Charges = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
    isos = [1, 2, 3, 4, 5]
    layout = [[sg.Text("Choose CSV File with Peptides in Matrix:"), sg.FileBrowse(key="-CSV-")],
              [sg.Text("Assumed Charge of Peptides: "), sg.Combo(Charges, font=('Arial Bold', 8), expand_x=True, default_value = 2, enable_events=True, readonly=False, key='-Charge-')],
              [sg.Text("Assumed Number of Isoforms: "), sg.Combo(isos, font=('Arial Bold', 8), expand_x=True, default_value=3, enable_events=True, readonly=False, key='-Isoforms-')],
              [sg.Button("<--"), sg.Button("Generate Spectra!"), sg.Button("MS2 Analysis")]]
    window = sg.Window("Simulating MALDI-TOF Spectrometry", layout)
    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED or event == "<--":
            window.close()
            mainscreen()
            break
        if event == "Generate Spectra!":
            # code exeption handling
            file = values["-CSV-"]
            if os.path.basename(file)[-4:] != ".csv":
                sg.popup("Input file was not in CSV format")
                event = ""
            elif values['-Charge-'] not in Charges:
                sg.popup("Charge given is not valid")
                event = ""
            elif values['-Isoforms-'] not in isos:
                sg.popup("Isoform number given is not valid")
                event = ""
            else:
                input_list = read(file)
                mz_isotope_dict = {}
                charge = values['-Charge-']
                lowerx = 10000000000000000000000000000000000000000000000
                upperx = 0
                for seq in input_list:
                    num_isotope = values['-Isoforms-']
                    # iterate and calculate individual isotope m/z value and also intensity value
                    tuple_list = []
                    mz = mass.calculate_mass(sequence=seq, ion_type="M", charge=int(charge))
                    if lowerx > mz:
                        lowerx = mz
                    if upperx < mz:
                        upperx = mz
                    for i in range(num_isotope):
                        isotope = mz + (i / int(charge))
                        intensity = (num_isotope - i) / num_isotope
                        tuple_data = (isotope, intensity)
                        tuple_list.append(tuple_data)  # add tuple to list of tuple
                    if seq not in mz_isotope_dict:
                        mz_isotope_dict[seq] = tuple_list

                traces = []
                for protein, mz_values in mz_isotope_dict.items():
                    x_values = []
                    y_values = []
                    for mz, i in mz_values:
                        x_values.extend([mz - 0.1, mz, mz, mz, mz + 0.1])
                        y_values.extend([0, 0, i, 0, 0])

                    if not x_values or not y_values:
                        print(f"No data points for category: {protein}")
                        continue

                    trace = go.Scatter({"x": x_values, "y": y_values, "visible": True, "mode": "lines", "line": {"width": 1, }})
                    traces.append(trace)

                if not traces:
                    print("No data to plot")
                else:
                    layout = go.Layout(
                        title="m/z Values of Proteins",
                        xaxis=dict(title="m/z", range=[0.9*lowerx, 1.1*upperx]),
                        yaxis=dict(title="Intensity"),
                    )
                    fig = go.Figure(data=traces, layout=layout)
                    fig.show()

        if event == "MS2 Analysis":
            # code exeption handling
            file = values["-CSV-"]
            if os.path.basename(file)[-4:] != ".csv":
                sg.popup("Input file was not in CSV format")
                event = ""
            elif values['-Charge-'] not in Charges:
                sg.popup("Charge given is not valid")
                event = ""
            else:
                input_list = read(file)
                mz_dict = {}
                mz_list = []
                charge = values['-Charge-']
                for seq in input_list:
                    # iterate and calculate individual isotope m/z value
                    mz = mass.calculate_mass(sequence=seq, ion_type="M", charge=int(charge))
                    if mz not in mz_dict:
                        mz_dict[mz] = seq
                        mz_list.append(mz)
                ms2(mz_dict, mz_list)
                event = ""
def OrbiSpec():
    methodsTranslation = {'Reverse Phase at pH 2.0': achrom.RCs_guo_ph2_0,
                          'Reverse Phase at pH 7.0': achrom.RCs_guo_ph7_0,
                          # '': achrom.RCs_meek_ph2_1,
                          # '': achrom.RCs_meek_ph7_4,
                          # '': achrom.RCs_browne_tfa,
                          # '': achrom.RCs_browne_hfba,
                          # '': achrom.RCs_palmblad,
                          'Normal phase chromatography': achrom.RCs_yoshida,
                          # '': achrom.RCs_yoshida_lc,
                          # '': achrom.RCs_zubarev,
                          'Hydrophilic interaction liquid chromatography at pH 3.0': achrom.RCs_gilar_atlantis_ph3_0,
                          'Hydrophilic interaction liquid chromatography at pH 4.5': achrom.RCs_gilar_atlantis_ph4_5,
                          'Hydrophilic interaction liquid chromatography at pH 10.0': achrom.RCs_gilar_atlantis_ph10_0,
                          # '': achrom.RCs_gilar_beh,
                          # '': achrom.RCs_gilar_beh_amide,
                          # '': achrom.RCs_gilar_rp,
                          # '': achrom.RCs_krokhin_100A_fa,
                          # '': achrom.RCs_krokhin_100A_tfa,
                          }
    methods = ['Reverse Phase at pH 2.0', 'Reverse Phase at pH 7.0', 'Normal phase chromatography',
               'Hydrophilic interaction liquid chromatography at pH 3.0',
               'Hydrophilic interaction liquid chromatography at pH 4.5',
               'Hydrophilic interaction liquid chromatography at pH 10.0']
    Charges = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
    isos = [1, 2, 3, 4, 5]
    # layout1 = [[sg.Text("Electrospray Column Settings")],
    #            [sg.Text("Choose CSV File:"), sg.FileBrowse(key="-CSV-")],
    #            [sg.Text("Type of Chromatography:"), sg.Combo(methods, font=('Arial Bold', 8), expand_x=True, enable_events=True, readonly=False, key='-methods-')],
    #            [sg.Text("Chromatography run time"), sg.Input(key='-RunTime-', s=4), sg.Button("Show Elution Profile")],
    #            [sg.Button("<--")]]
    # layout2 = [[sg.Text("Mass Spec Settings")],
    #            [sg.Text("Assumed Charge of Peptides: "), sg.Combo(Charges, font=('Arial Bold', 8), expand_x=True, default_value = 2, enable_events=True, readonly=False, key='-Charge-')],
    #            [sg.Text("Assumed Number of Isoforms: "), sg.Combo(isos, font=('Arial Bold', 8), expand_x=True, default_value=3, enable_events=True, readonly=False, key='-Isoforms-')],
    #            [sg.Text("Select Mass Spec Scan Time"), sg.Input(key='-ScanTime-', s=4)],
    #            [sg.Button("Generate Spectra!")]]
    # layout = [[sg.Column(layout1), sg.VSeperator(), sg.Column(layout2)]]
    layout3 = [[sg.Text("Electrospray Column Settings")],
               [sg.Text("Choose CSV File:"), sg.FileBrowse(key="-CSV-")],
               [sg.Text("Type of Chromatography:"), sg.Combo(methods, font=('Arial Bold', 8), expand_x=True, enable_events=True, readonly=False, key='-methods-')],
               [sg.Text("Chromatography run time"), sg.Input(key='-RunTime-', s=4), sg.Button("Show Elution Profile")],
               [sg.HSeparator()],
               [sg.Text("Mass Spec Settings")],
               [sg.Text("Assumed Charge of Peptides: "), sg.Combo(Charges, font=('Arial Bold', 8), expand_x=True, default_value=2, enable_events=True,readonly=False, key='-Charge-')],
               [sg.Text("Assumed Number of Isoforms: "), sg.Combo(isos, font=('Arial Bold', 8), expand_x=True, default_value=3, enable_events=True,readonly=False, key='-Isoforms-')],
               [sg.Text("Select Mass Spec Scan Time"), sg.Input(key='-ScanTime-', s=4)],
               [sg.Button("<--"), sg.Button("Generate Spectra!"), sg.Button("MS2 Analysis")]]
    window = sg.Window("Simulating Oribitrap Spectrometry", layout3)
    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED or event == "<--":
            window.close()
            mainscreen()
            break
        if event == "Show Elution Profile":
            # code exeption handling
            file = values["-CSV-"]
            if os.path.basename(file)[-4:] != ".csv":
                sg.popup("Input file was not in CSV format")
                event = ""
            elif values['-methods-'] not in methods:
                sg.popup("Chromatography method given was invalid")
                event = ""
            else:
                try:
                    int(values['-RunTime-'])
                except:
                    sg.popup("Chromatography run time given was not an integer")
                    event = ""
                else:
                    input_list = read(file)
                    retention_ls = []
                    for seq in input_list:
                        retention_ls.append(pyteomics.achrom.calculate_RT(seq, methodsTranslation[values['-methods-']]))
                    normalized_time = []
                    for time in retention_ls:
                        normalized_time.append((((time / max(retention_ls)) - min(retention_ls)) * int(values['-RunTime-'])))
                    fig2 = go.Figure()
                    fig2.add_trace(go.Histogram(x=normalized_time, marker=dict(color=[0, 255, 255])))
                    fig2.update_xaxes(title_text="Retention Time (min)")
                    fig2.update_yaxes(title_text="Number of sequences")
                    fig2.show()
        if event == "Generate Spectra!":
            # code exeption handling
            file = values["-CSV-"]
            if os.path.basename(file)[-4:] != ".csv":
                sg.popup("Input file was not in CSV format")
                event = ""
            elif values['-methods-'] not in methods:
                sg.popup("Chromatography method given was invalid")
                event = ""
            elif values['-Charge-'] not in Charges:
                sg.popup("Charge given is not valid")
                event = ""
            elif values['-Isoforms-'] not in isos:
                sg.popup("Isoform number given is not valid")
                event = ""
            else:
                try:
                    int(values['-RunTime-'])
                except:
                    sg.popup("Chromatography run time given was not an integer")
                    event = ""
                else:
                    try:
                        int(values['-ScanTime-'])
                    except:
                        sg.popup("Scan time time given was not an integer")
                        event = ""
                    else:
                        if int(values['-ScanTime-']) > int(values['-RunTime-']):
                            sg.popup("Scan time time given was not within the electrospray runtime")
                            event = ""
                        else:
                            input_list = read(file)
                            retention_ls = []
                            for seq in input_list:
                                retention_ls.append(pyteomics.achrom.calculate_RT(seq, methodsTranslation[values['-methods-']]))
                            normalized_time = []

                            for time in retention_ls:
                                normalized_time.append(normalized_time.append(((time - min(retention_ls)) / (max(retention_ls) - min(retention_ls))) * int(values['-RunTime-'])))

                            user_seq = []
                            for i in range(len(normalized_time)):
                                if int(values['-ScanTime-'])-1 <= normalized_time[i] <= int(values['-ScanTime-']):
                                    user_seq.append(input_list[i])

                            mz_isotope_dict = {}
                            charge = values['-Charge-']
                            lowerx = 10000000000000000000000000000000000000000000000
                            upperx = 0
                            for seq in user_seq:
                                num_isotope = values['-Isoforms-']
                                # iterate and calculate individual isotope m/z value and also intensity value
                                tuple_list = []
                                mz = mass.calculate_mass(sequence=seq, ion_type="M", charge=int(charge))
                                if lowerx > mz:
                                    lowerx = mz
                                if upperx < mz:
                                    upperx = mz
                                for i in range(num_isotope):
                                    isotope = mz + (i / int(charge))
                                    intensity = (num_isotope - i) / num_isotope
                                    tuple_data = (isotope, intensity)
                                    tuple_list.append(tuple_data)  # add tuple to list of tuple
                                if seq not in mz_isotope_dict:
                                    mz_isotope_dict[seq] = tuple_list

                            traces = []
                            for protein, mz_values in mz_isotope_dict.items():
                                x_values = []
                                y_values = []
                                for mz, i in mz_values:
                                    x_values.extend([mz - 0.1, mz, mz, mz, mz + 0.1])
                                    y_values.extend([0, 0, i, 0, 0])

                                if not x_values or not y_values:
                                    print(f"No data points for category: {protein}")
                                    continue

                                trace = go.Scatter({"x": x_values, "y": y_values, "visible": True, "mode": "lines", "line": {"width": 1, }})
                                traces.append(trace)

                            if not traces:
                                print("No data to plot")
                            else:
                                layout = go.Layout(
                                    title="m/z Values of Proteins",
                                    xaxis=dict(title="m/z", range=[0.9*lowerx, 1.1*upperx]),
                                    yaxis=dict(title="Intensity"),
                                )
                                fig = go.Figure(data=traces, layout=layout)
                                fig.show()
        if event == "MS2 Analysis":
            file = values["-CSV-"]
            if os.path.basename(file)[-4:] != ".csv":
                sg.popup("Input file was not in CSV format")
                event = ""
            elif values['-methods-'] not in methods:
                sg.popup("Chromatography method given was invalid")
                event = ""
            elif values['-Charge-'] not in Charges:
                sg.popup("Charge given is not valid")
                event = ""
            elif values['-Isoforms-'] not in isos:
                sg.popup("Isoform number given is not valid")
                event = ""
            else:
                try:
                    int(values['-RunTime-'])
                except:
                    sg.popup("Chromatography run time given was not an integer")
                    event = ""
                else:
                    try:
                        int(values['-ScanTime-'])
                    except:
                        sg.popup("Scan time time given was not an integer")
                        event = ""
                    else:
                        if int(values['-ScanTime-']) > int(values['-RunTime-']):
                            sg.popup("Scan time time given was not within the electrospray runtime")
                            event = ""
                        else:
                            input_list = read(file)
                            retention_ls = []
                            for seq in input_list:
                                retention_ls.append(pyteomics.achrom.calculate_RT(seq, methodsTranslation[values['-methods-']]))
                            normalized_time = []

                            for time in retention_ls:
                                normalized_time.append((((time / max(retention_ls)) - min(retention_ls)) * int(values['-RunTime-'])))

                            user_seq = []
                            for i in range(len(normalized_time)):
                                if int(values['-ScanTime-']) - 1 <= normalized_time[i] <= int(values['-ScanTime-']):
                                    user_seq.append(input_list[i])
                            mz_list = []
                            mz_dict = {}
                            charge = values['-Charge-']
                            for seq in user_seq:
                                # iterate and calculate individual isotope m/z value
                                mz = mass.calculate_mass(sequence=seq, ion_type="M", charge=int(charge))
                                if mz not in mz_dict:
                                    mz_dict[mz] = seq
                                    mz_list.append(mz)
                            ms2(mz_dict, mz_list)
                            event = ""
def ms2(mz_dict, mz_list):
    # pprint(mz_dict)
    # pprint(mz_list)
    Charges = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
    if len(mz_list) == 0:
        sg.popup("There are no peptides in the MS1 spectrum")
    else:
        layout = [[sg.Text("Choose Precursor MZ for MS2: "), sg.Combo(mz_list, font=('Arial Bold', 8), expand_x=True, enable_events=True, readonly=False, key='-MZ-')],
                  [sg.Text("Assumed Charge of Peptides: "), sg.Combo(Charges, font=('Arial Bold', 8), expand_x=True, enable_events=True, readonly=False, key='-charge-')],
                  [sg.Button("Run MS2!")],
                  [sg.HSeparator()],
                  [sg.Button("Add Post Translational Modifications")],
                  [sg.Text("PTM Position: "), sg.Text("0", key='-PTMPos-', enable_events=True)],
                  [sg.Text("PTM Mass: "), sg.Text("0", key='-PTMMass-', enable_events=True)]]
        window = sg.Window("Simulating Oribitrap Spectrometry", layout)
        while True:
            event, values = window.read()
            if event == sg.WIN_CLOSED:
                window.close()
                event = ""
                break
            if event == "-MZ-":
                window["-PTMPos-"].update(0)
                window["-PTMMass-"].update(0)
            if event == "Add Post Translational Modifications":
                if values['-MZ-'] not in mz_dict.keys():
                    sg.popup("Precursor chosen is not in MS1 spectrum")
                    event = ""
                else:
                    sequence = mz_dict[values['-MZ-']]
                    PTMmass, pos = PTM(sequence)
                    window["-PTMPos-"].update(pos)
                    window["-PTMMass-"].update(PTMmass)
            if event == "Run MS2!":
                if values['-MZ-'] not in mz_dict.keys():
                    sg.popup("Precursor chosen is not in MS1 spectrum")
                    event = ""
                elif values['-charge-'] not in Charges:
                    sg.popup("Charge given is not valid")
                    event = ""
                else:
                    sequence = mz_dict[values['-MZ-']]
                    charge = values['-charge-']
                    upper_xbound = 0
                    lower_xbound = 99999999999999999999999999999999999999999999999999999999999
                    PTM_Mass = int(window["-PTMMass-"].DisplayText)
                    PTM_Pos = int(window["-PTMPos-"].DisplayText)
                    sg.Popup(PTM_Mass)
                    # window["-PTMPos-"].DisplayText
                    ms2_mz_dict = {'b_ion': [], 'y_ion': []}
                    b_ion_ls = []
                    y_ion_ls = []
                    for i in range(1, len(sequence)):
                        b_ion = sequence[:i]
                        y_ion = sequence[i:]

                        if PTM_Pos == 0 or PTM_Mass == 0:
                            b_ion_ls.append(mass.calculate_mass(sequence=b_ion, ion_type="b", charge=charge))
                            y_ion_ls.append(mass.calculate_mass(sequence=y_ion, ion_type="y", charge=charge))
                        else:
                            if i >= PTM_Pos:
                                b_ion_ls.append(mass.calculate_mass(sequence=b_ion, ion_type="b", charge=charge) +
                                                (PTM_Mass / abs(charge)))
                                y_ion_ls.append(mass.calculate_mass(sequence=y_ion, ion_type="y", charge=charge))
                            else:
                                b_ion_ls.append((mass.calculate_mass(sequence=b_ion, ion_type="b", charge=charge)))
                                y_ion_ls.append(mass.calculate_mass(sequence=y_ion, ion_type="y", charge=charge) +
                                                (PTM_Mass / abs(charge)))

                    ms2_mz_dict['b_ion'] = b_ion_ls
                    ms2_mz_dict['y_ion'] = y_ion_ls

                    traces = []
                    for ion, mz_values in ms2_mz_dict.items():
                        if lower_xbound > min(mz_values):
                            lower_xbound = min(mz_values)
                        if upper_xbound < max(mz_values):
                            upper_xbound = max(mz_values)

                        x_values = []
                        y_values = []

                        for mz in mz_values:
                            x_values.extend([mz - 0.1, mz, mz, mz, mz + 0.1])
                            y_values.extend([0, 0, 1, 0, 0])

                        trace = go.Scatter({
                                "x": x_values,
                                "y": y_values,
                                "visible": True,
                                "mode": "lines",
                                "line": {"width": 1, },
                            })
                        traces.append(trace)

                    b_annotations = [dict(x=b_ion_ls[i], y=1, xref="x", yref="y", text=f"b{i + 1}", showarrow=True,
                                          arrowhead=7, ax=0, ay=-40) for i in range(len(b_ion_ls))]

                    y_annotations = [dict(x=y_ion_ls[i], y=1, xref="x", yref="y", text=f"y{len(y_ion_ls) - i}", showarrow=True,
                                          arrowhead=7, ax=0, ay=-25) for i in range(len(y_ion_ls))]
                    layout = go.Layout(
                        title="m/z Values of B and Y ions",
                        xaxis=dict(title="m/z", range=[0.9*lower_xbound, 1.1*upper_xbound]),
                        yaxis=dict(title="Intensity"),
                        annotations=b_annotations + y_annotations,
                        showlegend=False)

                    fig = go.Figure(data=traces, layout=layout)
                    fig.show()
def PTM(seq):
    seqlen = len(seq)
    layout4 = [[sg.Text("Select PTM Position between 1 and " + str(seqlen) + ": "), sg.Input(key='-Pos-', s=5)],
              # [sg.Text(seq)],
              [sg.Text("Select PTM Mass: "), sg.Input(key='-Mass-', s=5)],
              [sg.Button("Cancel"), sg.Button("Add PTM")]]
    window = sg.Window("Add PTM", layout4)
    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED or event == "Cancel":
            window.close()
            return 0, 0
        if event == "Add PTM":
            try:
                pos = int(values['-Pos-'])
            except:
                sg.Popup("PTM position given is not an integer")
                event = ""
            else:
                try:
                    mass = int(values['-Mass-'])
                except:
                    sg.Popup("PTM mass given is not an integer")
                    event = ""
                else:
                    if 0 > pos or seqlen < pos:
                        sg.Popup("PTM Position entered is not within range")
                        event = ""
                    elif mass < 0:
                        sg.Popup("PTM Mass entered cannot be negative")
                    else:
                        window.close()
                        return mass, pos

mainscreen()
