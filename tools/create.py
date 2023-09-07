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

if "_d" in variant:
    deathrate = -3
else:
    deathrate = -7

if "8" in variant or "8" in mechanism:
    groupsize = 3
else:
    groupsize = 2

if "cost" in variant:
    cost = -8
else:
    cost = -14

if "p" in mechanism:
    partnerchoice = 1
else:
    partnerchoice = 0

if "i" in mechanism:
    reciprocity = 1
    indirectr = 1
elif "r" in mechanism:
    reciprocity = 1
    indirectr = 0
else:
    reciprocity = 0
    indirectr = 0

if "noImimic" in variant:
    independent = 0
else:
    independent = 1

if "l" in mechanism:
    language = 1
else:
    language = 0

if "_shuffle" in variant:
    shuffle = 1
else:
    shuffle = 0

if "d_s" in variant or "d_n" in variant:
    discrete = 1
    a2init = 0.6
    a2high = 0.8
else:
    discrete = 0
    a2init = 0.1
    a2high = 0.1

# create subfolder path variant/mechanism/given in current directory
path = variant + "/" + mechanism + "/" + given
os.makedirs(path, exist_ok=True)

num = 21
alphas = np.linspace(0.1, 0.9, num)
logess = np.linspace(-5.0, 5.0, num)
Given = float(given[-3:]) / 100
c = 101

for alpha in alphas:
    for loges in logess:
        filename = path + "/" + str(c) + ".glo"
        f = open(filename, "w")
        f.write("Seed,1\n")
        f.write("N,12\n")
        # standard Runs,30
        f.write("Runs,30\n")
        # standard Time,21
        if discrete:
            f.write("Time,20\n")
        else:
            f.write("Time,21\n")
        f.write("Periods,3\n")
        f.write("a1Max,1.0\n")
        f.write("a2Max,1.0\n")
        f.write("R1,2\n")
        f.write("R2,2\n")
        f.write(f"a2Init,{a2init}\n")
        f.write("ChooseGrainInit,1.0\n")
        f.write("MimicGrainInit,1.0\n")
        f.write("ImimicGrainInit,1.0\n")
        f.write("a2MutationSize,-6\n")
        f.write("GrainMutationSize,-6\n")
        f.write(f"DeathRate,{deathrate}\n")
        f.write(f"GroupSize,{groupsize}\n")
        f.write(f"ChooseCost,{cost}\n")
        f.write(f"MimicCost,{cost}\n")
        f.write(f"ImimicCost,{cost}\n")
        f.write(f"PartnerChoice,{partnerchoice}\n")
        f.write(f"Reciprocity,{reciprocity}\n")
        f.write(f"IndirectR,{indirectr}\n")
        f.write(f"Independent,{independent}\n")
        f.write(f"Language,{language}\n")
        f.write(f"Shuffle,{shuffle}\n")
        f.write(f"Discrete,{discrete}\n")
        f.write(f"a2low,{a2init}\n")
        f.write(f"a2high,{a2high}\n")
        f.write(f"alpha,{alpha:.6}\n")
        f.write(f"logES,{loges}\n")
        f.write(f"Given,{Given}\n")

        f.close()
        c = c + 1
