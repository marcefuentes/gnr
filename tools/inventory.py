#! /usr/bin/env python

import csv
import os
import sys

nlines = 10
input_file_extension = ".glo"
output_file_extension = ".csv"

blue = "\033[94m"
cyan = "\033[96m"
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

Independent = 1
Shuffle = 1
Discrete = 0
DeathRate = -7

current_dir = os.getcwd().split("/")[-1]
if current_dir[0] == "d":
    Discrete = 1
if "noshuffle" in current_dir:
    Shuffle = 0
if "noImimic" in current_dir:
    Independent = 0
if "_d" in current_dir:
    DeathRate = -3

mechanisms = [f for f in os.listdir(os.getcwd()) if os.path.isdir(f)]
mechanisms.sort()
for mechanism in mechanisms:
    print(f"{cyan}{mechanism} {reset_format}", end = "")
    if "p" in mechanism:
        PartnerChoice = 1
    else:
        PartnerChoice = 0
    if "i" in mechanism:
        Reciprocity = 1
        IndirectR = 1
    elif "r" in mechanism:
        Reciprocity = 1
        IndirectR = 0
    else:
        Reciprocity = 0
        IndirectR = 0
    if "l" in mechanism:
        Language = 1
    else:
        Language = 0
    if "8" in mechanism:
        GroupSize = 3
    else:
        GroupSize = 2
    givens = [f for f in os.listdir(mechanism) if os.path.isdir(os.path.join(mechanism, f))]
    if len(givens) == 0:
        print(f"{red}empty{reset_format}")
        continue
    givens.sort()
    for given in givens:
        Given = float(given[-3:]) / 100
        given_path = os.path.join(mechanism, given)
        input_files = [f for f in os.listdir(given_path) if f.endswith(input_file_extension)]
        if given != givens[0]:
            print(" " * (len(mechanism) + 1), end = "")
        print(f"{given}: ", end = "")
        if len(input_files) == 0:
            print(f"{red}no {input_file_extension[1:]} files{reset_format}")
            continue
        pass_params = True
        data_dict = {}
        with open(os.path.join(given_path, input_files[0]), "r") as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                key, value = row
                if value.isdigit():
                    data_dict[key] = int(value)
                if key == "Given" or key == "DeathRate":
                    data_dict[key] = float(value)
        if data_dict["Given"] != Given:
            print(f"{red}Given = {data_dict['Given']}{reset_format}", end = " ")
            pass_params = False
        if data_dict["DeathRate"] != DeathRate:
            print(f"{red}DeathRate = {data_dict['DeathRate']}{reset_format}", end = " ")
            pass_params = False
        if data_dict["GroupSize"] != GroupSize:
            print(f"{red}GroupSize = {data_dict['GroupSize']}{reset_format}", end = " ")
            pass_params = False
        if data_dict["PartnerChoice"] != PartnerChoice:
            print(f"{red}PartnerChoice = {data_dict['PartnerChoice']}{reset_format}", end = " ")
            pass_params = False
        if data_dict["Reciprocity"] != Reciprocity:
            print(f"{red}Reciprocity = {data_dict['Reciprocity']}{reset_format}", end = " ")
            pass_params = False
        if data_dict["IndirectR"] != IndirectR:
            print(f"{red}IndirectR = {data_dict['IndirectR']}{reset_format}", end = " ")
            pass_params = False
        if data_dict["Independent"] != Independent:
            print(f"{red}Independent = {data_dict['Independent']}{reset_format}", end = " ")
            pass_params = False
        if data_dict["Language"] != Language:
            print(f"{red}Language = {data_dict['Language']}{reset_format}", end = " ")
            pass_params = False
        if data_dict["Shuffle"] != Shuffle:
            print(f"{red}Shuffle = {data_dict['Shuffle']}{reset_format}", end = " ")
            pass_params = False
        if data_dict["Discrete"] != Discrete:
            print(f"{red}Discrete = {data_dict['Discrete']}{reset_format}", end = " ")
            pass_params = False

        if pass_params:
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
                print(f"{yellow}ok{reset_format}")
            else:
                first = True
                if f_equal_nlines:
                    print(f"{red}{f_equal_nlines}{reset_format} finished", end = "")
                    first = False
                if f_smaller_nlines:
                    if not first:
                        print(", ", end = "")
                    print(f"{red}{f_smaller_nlines}{reset_format} incomplete", end = "")
                    first = False
                notstarted = len(input_files) - f_smaller_nlines - f_equal_nlines - f_larger_nlines
                if notstarted:
                    if not first:
                        print(", ", end = "")
                    print(f"{red}{notstarted}{reset_format} not started", end = "")
                    first = False
                if f_larger_nlines:
                    if not first:
                        print(", ", end = "")
                    print(f"{red}{f_larger_nlines}{reset_format} > {nlines} lines", end = "")
                print()
        else:
            print() 
