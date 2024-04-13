import csv
def writer(list, output_file_path):
    # Call this function like this: writer(<the list you want to save>, <the file path>)
    with open(output_file_path, "w") as csv_out:
        writer = csv.writer(
            csv_out, lineterminator="\n"
        )
        for seq in list:
            writer.writerow([seq])