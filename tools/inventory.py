#! /usr/bin/env python

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
    os.chdir(sys.argv[1])

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
            
