#! /usr/bin/env python

import numpy as np
import os
import sys

# the script accepts exactly three arguments
if len(sys.argv) != 4:
    print("Usage: python create.py <variant> <mechanism> <given>")
    sys.exit()

variant = sys.argv[1]
mechanism = sys.argv[2]
given = sys.argv[3]

# create subfolder path variant/mechanism/given in current directory
path = variant + "/" + mechanism + "/" + given
os.makedirs(path, exist_ok=True)

num = 21
alphas = np.linspace(0.1, 0.9, num)
logess = np.linspace(-5.0, 5.0, num)
c = 101

for alpha in alphas:
    for loges in logess:
        filename = path + "/" + str(c) + ".glo"
        f = open(filename, "w")

        f.write("Seed,1\n")
        f.write("N,12\n")
        # standard Runs,30
        f.write("Runs,10\n")
        # standard Time,21
        f.write("Time,20\n")
        f.write("Periods,3\n")
        f.write("a1Max,1.0\n")
        f.write("a2Max,1.0\n")
        f.write("R1,2\n")
        f.write("R2,2\n")
        f.write("a2Init,0.1\n")
        f.write("ChooseGrainInit,1.0\n")
        f.write("MimicGrainInit,1.0\n")
        f.write("ImimicGrainInit,1.0\n")
        f.write("a2MutationSize,-6\n")
        f.write("GrainMutationSize,-6\n")
        if "_d" in variant:
            f.write("DeathRate,-3\n")
        else:
            f.write("DeathRate,-7\n")
        if "8" in variant or "8" in mechanism:
            f.write("GroupSize,3\n")
        else:
            f.write("GroupSize,2\n")
        # standard costs,-14
        f.write("ChooseCost,-7\n")
        f.write("MimicCost,-7\n")
        f.write("ImimicCost,-7\n")
        if "p" in mechanism:
            f.write("PartnerChoice,1\n")
        else:
            f.write("PartnerChoice,0\n")
        if "i" in mechanism:
            f.write("Reciprocity,1\n")
            f.write("IndirectR,1\n")
        elif "r" in mechanism:
            f.write("Reciprocity,1\n")
            f.write("IndirectR,0\n")
        else:
            f.write("Reciprocity,0\n")
            f.write("IndirectR,0\n")
        if "noImimic" in variant:
            f.write("Independent,0\n")
        else:
            f.write("Independent,1\n")
        if "l" in mechanism:
            f.write("Language,1\n")
        else:
            f.write("Language,0\n")
        if "_shuffle" in variant:
            f.write("Shuffle,1\n")
        else:
            f.write("Shuffle,0\n")
        if "d_s" in variant or "d_n" in variant:
            f.write("Discrete,1\n")
        else:
            f.write("Discrete,0\n")
        f.write(f"a2low,0.0\n")
        f.write(f"a2high,0.0\n")
        f.write(f"alpha,{alpha:.6}\n")
        f.write(f"logES,{loges}\n")
        Given = float(given[-3:]) / 100
        f.write(f"Given,{Given}\n")

        f.close()
        c = c + 1
