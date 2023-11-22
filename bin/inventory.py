#! /usr/bin/env python

import csv
import os
import sys

nlines = 10
input_file_extension = ".glo"
output_file_extension = ".csv"

blue = "\033[94m"
cyan = "\033[96m"
green = "\033[92m"
red = "\033[91m"
yellow = "\033[33m"
bold = "\033[1m"
reset_format = "\033[0m"

if len(sys.argv) > 1:
    if os.path.isdir(sys.argv[1]):
        os.chdir(sys.argv[1])
    else:
        print(f"{red}Directory {sys.argv[1]} does not exist{reset_format}")
        exit()

folder_dict = {}
folder_dict["Independent"] = 1
folder_dict["Shuffle"] = 0
folder_dict["Discrete"] = 0
folder_dict["DeathRate"] = -7

current_dir = os.getcwd().split("/")[-1]
if "_shuffle" in current_dir:
    folder_dict["Shuffle"] = 1
if "noImimic" in current_dir:
    folder_dict["Independent"] = 0
if "_d" in current_dir:
    folder_dict["DeathRate"] = -3

mechanisms = [f for f in os.listdir(os.getcwd()) if os.path.isdir(f)]
#mechanisms.sort(key = lambda x: os.path.getmtime(x))
mechanisms.sort()
for mechanism in mechanisms:
    if "p" in mechanism:
        folder_dict["PartnerChoice"] = 1
    else:
        folder_dict["PartnerChoice"] = 0
    if "i" in mechanism:
        folder_dict["Reciprocity"] = 1
        folder_dict["IndirectR"] = 1
    elif "r" in mechanism:
        folder_dict["Reciprocity"] = 1
        folder_dict["IndirectR"] = 0
    else:
        folder_dict["Reciprocity"] = 0
        folder_dict["IndirectR"] = 0
    if "l" in mechanism:
        folder_dict["Language"] = 1
    else:
        folder_dict["Language"] = 0
    if "_8" in current_dir or "_8" in mechanism:
        folder_dict["GroupSize"] = 3
    else:
        folder_dict["GroupSize"] = 2
    # get the two-digit number after the string "cost" in "mechanism" and store it in variable "cost"

    givens = [f for f in os.listdir(mechanism) if os.path.isdir(os.path.join(mechanism, f))]
    if len(givens) == 0:
        print(f"{red}empty{reset_format}")
        continue
    givens.sort()
    for given in givens:
        folder_dict["Given"] = float(given[-3:]) / 100
        given_path = os.path.join(mechanism, given)
        input_files = [f for f in os.listdir(given_path) if f.endswith(input_file_extension)]
        print(f"{cyan}{mechanism}{reset_format}\t{blue}{given}{reset_format}", end = "  ")
        if len(input_files) == 0:
            print(f"{red}no {input_file_extension[1:]} files{reset_format}")
            continue
        data_dict = {}
        with open(os.path.join(given_path, input_files[0]), "r") as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                key, value = row
                if value.isdigit():
                    data_dict[key] = int(value)
                if key == "Given" or key == "DeathRate" or "Cost" or "a2" in key:
                    data_dict[key] = float(value)
        pass_params = True
        for key, folder_value in folder_dict.items():
            if data_dict[key] != folder_value:
                print(f"{red}{key} {data_dict[key]}{reset_format}", end = " ")
        f_smaller_nlines = 0
        f_equal_nlines = 0
        f_larger_nlines = 0
        for f in os.listdir(given_path):
            if f.endswith(output_file_extension):
                output_file = os.path.join(given_path, f)
                with open(output_file, "r") as output:
                    lines = output.readlines()
                    if len(lines) < nlines:
                        f_smaller_nlines += 1
                    elif len(lines) == nlines:
                        f_equal_nlines += 1
                    elif len(lines) > nlines:
                        f_larger_nlines += 1
        if f_equal_nlines == len(input_files):
            print(f"{green}{f_equal_nlines}{reset_format}")
        else:
            notstarted = len(input_files) - f_smaller_nlines - f_equal_nlines - f_larger_nlines
            if notstarted:
                print(f"{red}{notstarted}{reset_format}", end = " ")
            if f_smaller_nlines:
                print(f"{yellow}{f_smaller_nlines}{reset_format}", end = " ")
            if f_equal_nlines:
                print(f"{green}{f_equal_nlines}{reset_format}", end = " ")
            if f_larger_nlines:
                print(f"{red}{f_larger_nlines}{reset_format} > {nlines} lines", end = "")
            print()
