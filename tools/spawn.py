#! /usr/bin/env python

import numpy as np

given = 0.95
num = 21
alphas = np.linspace(0.1, 0.9, num)
logess = np.linspace(-5.0, 5.0, num)
c = 101

for alpha in alphas:
    for loges in logess:
        filename = str(c) + ".glo"
        f = open(filename, "w")

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
        f.write("Discrete,0\n")
        f.write(f"a2low,0.0\n")
        f.write(f"a2high,0.0\n")
        f.write("IndirectR,0\n")
        f.write(f"alpha,{alpha:.6}\n")
        f.write(f"logES,{loges}\n")
        f.write(f"Given,{given}\n")

        f.close()
        c = c + 1
