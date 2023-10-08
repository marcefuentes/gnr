#! /usr/bin/env python

from glob import glob
import os
import time

import pandas as pd

start_time = time.perf_counter()
this_file = os.path.basename(__file__)
file_name = this_file.split(".")[0]

# Options

columns = ["alpha",
           "logES",
           "Given",
           "Time",
           "wmean",
           "a2Seenmean",
           "ChooseGrainmean",
           "MimicGrainmean",
           "ImimicGrainmean"]

# The object of this script is to create a new csv file from the data in each basename subdirectory of the current directory. The new csv file will be named after the basename subdirectory.

# Get the names of the subdirectories in the current directory
subdirectories = glob("*/")

# Loop through the subdirectories
for subdirectory in subdirectories:
    # cd into the subdirectory
    os.chdir(subdirectory)
    # Get the names of the subdirectories in this subdirectory
    mechanisms = glob("*/")
    # Loop through the subdirectories
    for mechanism in mechanisms:
        os.chdir(mechanism)
        # Get the names of the subdirectories in this subdirectory
        givens = glob("*/")
        # Loop through the subdirectories
        for given in givens:
            # Create a dataframe from the last row of each csv file in given
            df = pd.concat([pd.read_csv(f) for f in glob(f"{given}*.csv")])
            # the dataframe must contain only the columns listed in the columns list
            df = df[columns]
            # write the dataframe to a csv file named after the subdirectory
            df.to_csv(f"{given[:-1]}.csv", index=False)
            # delete the subdirectory and its contents
            os.system(f"rm -rf {given}")
        # cd out of the subdirectory
        os.chdir("../")
    # cd out of the subdirectory
    os.chdir("../")

end_time = time.perf_counter()
print(f"\nTime elapsed: {(end_time - start_time):.2f} seconds")
