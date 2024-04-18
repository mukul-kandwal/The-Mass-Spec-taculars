import csv 
from pyteomics import achrom
import plotly.graph_objects as go

def read(input_file_path):
    # use this function like this: input_list = read(<input_file_path>)
    with open(input_file_path, 'r') as read_obj:
        csv_reader = csv.reader(read_obj)
        list_of_csv = list(csv_reader)
        input_list = []
        for row in list_of_csv:
            input_list.append(row[0])
        return input_list

input_list = read('output.csv') #read in csv file and output a list of protein sequences

def calc_time(input_list):
    retention_ls = []
    
    for seq in input_list:
        retention_ls.append(achrom.calculate_RT(seq, achrom.RCs_krokhin_100A_tfa))
    return retention_ls

retention_ls = calc_time(input_list) #calculate retention list using protein sequences and output list of calculated retention time 

def normalization(retention_ls, max_time):
    normalized_time = []
    for time in retention_ls:
        normalized_time.append((time / max(retention_ls) * max_time)) 

    fig2 = go.Figure()
    fig2.add_trace(go.Histogram(x=normalized_time, marker=dict(color=[0, 255, 255])))
    fig2.update_xaxes(title_text="Retention Time (min)")
    fig2.update_yaxes(title_text="Number of sequences")
    fig2.show()
    return normalized_time
    
max_time = int(input("How long would you like to run the chromatography for? "))
lower_range = int(input("What is the lower range of time? "))
upper_range = int(input("What is the higher range? "))
normalized_time_ls = normalization(retention_ls, max_time)

def output_csv(input_list, normalized_time, lower_range, upper_range, output_file_path):
    user_seq = []
    user_time_range = []
    for i in range(len(normalized_time)):
        if normalized_time[i] >= lower_range and normalized_time[i]<= upper_range:
            user_time_range.append(normalized_time[i])
            user_seq.append(input_list[i])
    with open(output_file_path, "w") as csv_out:
        writer = csv.writer(
            csv_out, lineterminator="\n"
        )
        for seq in user_seq:
            writer.writerow([seq])

output_csv(input_list, normalized_time_ls, lower_range, upper_range, "/Users/trangphan/Documents/Master/spring_2024/proteomics/project/The-Mass-Spec-taculars/retention_time.csv")
