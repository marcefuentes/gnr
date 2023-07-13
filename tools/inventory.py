#! /usr/bin/env python

import os

nlines = 10

blue = "\033[94m"
cyan = "\033[96m"
red = "\033[91m"
yellow = "\033[33m"
bold = "\033[1m"
reset_format = "\033[0m"

mechanisms = [f for f in os.listdir(os.getcwd()) if os.path.isdir(f)]
mechanisms.sort()
for mechanism in mechanisms:
    print(f"{cyan}{mechanism}{reset_format}")
    givens = [f for f in os.listdir(mechanism) if os.path.isdir(os.path.join(mechanism, f))]
    if len(givens) == 0:
        print(f"  {red}empty{reset_format}")
        continue
    givens.sort()
    for given in givens:
        given_path = os.path.join(mechanism, given)
        glo_files = [f for f in os.listdir(given_path) if f.endswith(".glo")]
        if len(glo_files) == 0:
            print(f"  {given}: {red}no .glo files{reset_format}")
            continue
        f_equal_nlines = 0
        f_larger_nlines = 0
        f_smaller_nlines = 0
        for f in os.listdir(given_path):
            if f.endswith(".csv"):
                csv_file = os.path.join(given_path, f)
                with open(csv_file, "r") as csv:
                    lines = csv.readlines()
                    if len(lines) == nlines:
                        f_equal_nlines += 1
                    elif len(lines) > nlines:
                        f_larger_nlines += 1
                    elif len(lines) < nlines:
                        f_smaller_nlines += 1
        if f_equal_nlines == len(glo_files):
            print(f"  {given}: {yellow}ok{reset_format}")
        else:
            if f_equal_nlines + f_larger_nlines + f_smaller_nlines == 0:
                print(f"  {given}: {red}no .csv files{reset_format}")
            else:
                print(f"  {given}:")
                if f_larger_nlines > 0:
                    print(f"            {red}{f_larger_nlines} files with > {nlines} lines{reset_format}")
                if f_smaller_nlines > 0:
                    print(f"            {red}{f_smaller_nlines} files with < {nlines} lines{reset_format}")
                missing = len(glo_files) - f_equal_nlines - f_larger_nlines - f_smaller_nlines
                if missing > 0:
                    print(f"            {red}{missing} missing .csv files{reset_format}")
            
