#! /usr/bin/env python

import numpy as np

num = 21
givens = np.linspace(0.0, 1.0, num)
alphas = np.linspace(0.1, 0.9, num)
loges = 1.0
c = 101

for alpha in alphas:
    for given in givens:
        filename = str(c) + '.glo'
        f = open(filename, 'w')

        f.write('Seed,1\n')
        f.write('N,12\n')
        f.write('Runs,30\n')
        f.write('Time,20\n')
        f.write('Periods,3\n')
        f.write('a1Max,1.0\n')
        f.write('a2Max,1.0\n')
        f.write('R1,2\n')
        f.write('R2,2\n')
        f.write('a2Init,0.5\n')
        f.write('ChooseGrainInit,1.0\n')
        f.write('MimicGrainInit,1.0\n')
        f.write('a2MutationSize,-6\n')
        f.write('GrainMutationSize,-6\n')
        f.write('DeathRate,-7\n')
        f.write('GroupSize,2\n')
        f.write('ChooseCost,-14\n')
        f.write('MimicCost,-14\n')
        f.write('PartnerChoice,0\n')
        f.write('Reciprocity,0\n')
        f.write('Discrete,1\n')
        f.write('IndirectR,0\n')
        f.write(f'alpha,{alpha}\n')
        f.write(f'ES,{loges}\n')
        f.write(f'given,{given}\n')

        f.close()
        c = c + 1
