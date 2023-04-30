import os
import csv

# get the path to the current working directory
folder_path = os.getcwd()

# dictionary to store the number of lines for each csv file
csv_lines = {}

# loop through each file in the folder
for filename in os.listdir(folder_path):
    # check if the file has the extension csv
    if filename.endswith(".csv"):
        # read the file and count the number of lines
        with open(os.path.join(folder_path, filename), "r") as file:
            lines = file.readlines()
            num_lines = len(lines)
            # add the number of lines to the csv_lines dictionary
            csv_lines[filename] = num_lines

# dictionary to store the file numbers for each number of lines
lines_files = {}

# loop through the csv_lines dictionary and group the files by number of lines
for filename, num_lines in csv_lines.items():
    if num_lines not in lines_files:
        lines_files[num_lines] = []
    file_number = filename.split(".")[0]
    lines_files[num_lines].append(file_number)

# list to store deleted file numbers
deleted_files = []

# loop through the lines_files dictionary and delete the files with 1 line
for num_lines, file_numbers in lines_files.items():
    if num_lines == 1:
        for file_number in file_numbers:
            # delete the csv file
            csv_file = os.path.join(folder_path, f"{file_number}.csv")
            if os.path.exists(csv_file):
                os.remove(csv_file)
                deleted_files.append(file_number)
            # delete the gl2 file
            gl2_file = os.path.join(folder_path, f"{file_number}.gl2")
            if os.path.exists(gl2_file):
                os.remove(gl2_file)
            # delete the frq file
            frq_file = os.path.join(folder_path, f"{file_number}.frq")
            if os.path.exists(frq_file):
                os.remove(frq_file)

# create a CSV file with the deleted file numbers
with open("deleted_files.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Deleted File Numbers"])
    writer.writerow([",".join(set(deleted_files))])

print(f"{len(deleted_files)} files have been deleted.")

