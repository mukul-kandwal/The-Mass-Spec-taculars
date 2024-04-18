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

input_list = read('output.csv')

#calculate retention time
def calc_time(input_list):
    retention_ls = []
    
    for seq in input_list:
        retention_ls.append(achrom.calculate_RT(seq, achrom.RCs_krokhin_100A_tfa))
    return retention_ls

retention_ls = calc_time(input_list)

#plot the distribution of proteins eluting over range of time 
def plot_totalTime(retention_ls, max_time):
    times = []
    for time in retention_ls:
        if time < max_time:
            times.append(time)
    fig2 = go.Figure()
    fig2.add_trace(go.Histogram(x=times, marker=dict(color=[0, 255, 255])))
    fig2.update_xaxes(title_text="Retention Time (min)")
    fig2.update_yaxes(title_text="Number of sequences")
    fig2.show()

def normalize_time(input_list, retention_ls, max_time, lower_range, upper_range):
    normalized_times = []
    for time in retention_ls:
        if time < max_time:
            normalized_times.append(time/max_time)

    user_time_range = [] #list of retention time based on user defined time range
    user_seq = [] #list of sequences based on user defined time range
    low = lower_range/max_time
    high = upper_range/max_time
    for i in range(len(normalized_times)):
        if normalized_times[i] >= low or normalized_times[i]<= high:
            user_time_range.append(normalized_times[i])
            user_seq.append(input_list[i])

    return user_time_range, user_seq

def writer(user_sequence_ls, user_retention_ls, output_file_path):
    with open("retention_time.csv", "w") as csv_out:
        writer = csv.DictWriter(
            csv_out, fieldnames=["Sequences", "Normalized retention time"], lineterminator="\n")
        writer.writeheader()
        for n, prot in enumerate(user_sequence_ls):
            writer.writerow({"Sequences": prot[n], "Normalized retention time": user_retention_ls[n]})

max_time = int(input("How long would you like to run the chromatography for? "))
lower_range = int(input("What is the lower range of time? "))
upper_range = int(input("What is the higher range? "))
plot_totalTime(retention_ls, max_time)

user_time, user_seq = normalize_time(retention_ls, max_time, lower_range, upper_range)
writer(user_seq, user_time, "/Users/trangphan/Documents/Master/spring_2024/proteomics/project/The-Mass-Spec-taculars")

#normalization