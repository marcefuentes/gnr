#! /usr/bin/env python

import os

nlines = 10
input_file_extension = ".glo"
output_file_extension = ".csv"

blue = "\033[94m"
cyan = "\033[96m"
red = "\033[91m"
yellow = "\033[33m"
bold = "\033[1m"
reset_format = "\033[0m"

mechanisms = [f for f in os.listdir(os.getcwd()) if os.path.isdir(f)]
mechanisms.sort()
for mechanism in mechanisms:
    print(f"{cyan}{mechanism} {reset_format}", end = "")
    givens = [f for f in os.listdir(mechanism) if os.path.isdir(os.path.join(mechanism, f))]
    if len(givens) == 0:
        print(f"{red}empty{reset_format}")
        continue
    givens.sort()
    for given in givens:
        if given != givens[0]:
            print(" " * (len(mechanism) + 1), end = "")
        print(f"{given} ", end = "")
        given_path = os.path.join(mechanism, given)
        input_files = [f for f in os.listdir(given_path) if f.endswith(input_file_extension)]
        if len(input_files) == 0:
            print(f"{red}no {input_file_extension[1:]} files{reset_format}")
            continue
        f_equal_nlines = 0
        f_larger_nlines = 0
        f_smaller_nlines = 0
        for f in os.listdir(given_path):
            if f.endswith(output_file_extension):
                output_file = os.path.join(given_path, f)
                with open(output_file, "r") as output:
                    lines = output.readlines()
                    if len(lines) == nlines:
                        f_equal_nlines += 1
                    elif len(lines) > nlines:
                        f_larger_nlines += 1
                    elif len(lines) < nlines:
                        f_smaller_nlines += 1
        if f_equal_nlines == len(input_files):
            print(f"{yellow}ok{reset_format}")
        else:
            if f_equal_nlines + f_larger_nlines + f_smaller_nlines == 0:
                print(f"{red}0{reset_format} {output_file_extension[1:]} files")
            else:
                if f_larger_nlines:
                    print(f"{red}{f_larger_nlines}{reset_format} {output_file_extension[1:]} files with > {nlines} lines")
                if f_smaller_nlines:
                    if f_larger_nlines:
                        print(f"{' ' * (len(mechanism) + len(given) + 2)}{red}", end = "")
                    print(f"{red}{f_smaller_nlines}{reset_format} incomplete {output_file_extension[1:]} files")
                missing = len(input_files) - f_equal_nlines - f_larger_nlines - f_smaller_nlines
                if missing:
                    if f_larger_nlines or f_smaller_nlines:
                        print(f"{' ' * (len(mechanism) + len(given) + 2)}{red}", end = "")
                    print(f"{red}{missing}{reset_format} missing {output_file_extension[1:]} files")
            
