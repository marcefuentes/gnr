#! /usr/bin/env python

from glob import glob
import os
import pandas as pd
import mymodule as my

private_folder = "private"
social_folder = "social"
output_folder = "none"

for i in range(101, 542):

    private_file_path = os.path.join(private_folder, f"{i}.csv")
    dfprivate = my.read_file(private_file_path, False)
    a2low = dfprivate["a2Seenmean"].iloc[0]
    alpha = dfprivate["alpha"].iloc[0]
    logES = dfprivate["logES"].iloc[0]
    Given = dfprivate["Given"].iloc[0]

    social_file_path = os.path.join(social_folder, f"{i}.csv")
    dfsocial = my.read_file(social_file_path, False)
    a2high = dfsocial["a2Seenmean"].iloc[0]

    output_file_path = os.path.join(output_folder, f"{i}.glo")
    
    with open(output_file_path, "w") as f:
        f.write("Seed,1\n")
        f.write("N,12\n")
        f.write("Runs,30\n")
        f.write("Time,20\n")
        f.write("Periods,3\n")
        f.write("a1Max,1.0\n")
        f.write("a2Max,1.0\n")
        f.write("R1,2\n")
        f.write("R2,2\n")
        f.write("a2Init,0.5\n")
        f.write("ChooseGrainInit,1.0\n")
        f.write("MimicGrainInit,1.0\n")
        f.write("a2MutationSize,-6\n")
        f.write("GrainMutationSize,-6\n")
        f.write("DeathRate,-7\n")
        f.write("GroupSize,2\n")
        f.write("ChooseCost,-14\n")
        f.write("MimicCost,-14\n")
        f.write("PartnerChoice,0\n")
        f.write("Reciprocity,0\n")
        f.write("Discrete,1\n")
        f.write(f"a2low,{str(a2low)}\n")
        f.write(f"a2high,{str(a2high)}\n")
        f.write("IndirectR,0\n")
        f.write(f"alpha,{str(alpha)}\n")
        f.write(f"logES,{str(logES)}\n")
        f.write(f"Given,{str(Given)}\n")
