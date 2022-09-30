#! /usr/bin/env python

import numpy as np

num = 21
ess = np.linspace(-5.0, 5.0, num)
givens = np.linspace(0.0, 1.0, num)
c = 101

for given in givens:
    for es in ess:
        filename = str(c) + '.glo'
        f = open(filename, 'w')

        f.write('Seed,1\n')
        f.write('N,12\n')
        f.write('Runs,30\n')
        f.write('Time,21\n')
        f.write('Periods,2\n')
        f.write('FFunction,c\n')
        f.write('a1Max,1.0\n')
        f.write('a2Max,1.0\n')
        f.write('R1,2\n')
        f.write('R2,2\n')
        f.write('a2Init,0.1\n')
        f.write('a2High,0.8\n')
        f.write('ChooseGrainInit,1.0\n')
        f.write('MimicGrainInit,1.0\n')
        f.write('a2MutationSize,-7\n')
        f.write('GrainMutationSize,-6\n')
        f.write('DeathRate,-8\n')
        f.write('c1,0.0\n')
        f.write('c2,0.0\n')
        f.write('alpha,0.5\n')
        f.write('ES,0\n')
        f.write('GroupSize,2\n')
        f.write('ChooseCost,-14\n')
        f.write('MimicCost,-14\n')
        f.write('Given,0.0\n')
        f.write('PartnerChoice,1\n')
        f.write('Reciprocity,1\n')
        f.write('Optimal,0\n')
        f.write('Macromutation,-8\n')
        f.write('IndirectR,0\n')
        f.write('factorName1,Given\n')
        f.write(f'fFirst1,{given}\n')
        f.write(f'fLast1,{given}\n')
        f.write('fLevels1,1\n')
        f.write('factorName2,ES\n')
        f.write(f'fFirst2,{es}\n')
        f.write(f'fLast2,{es}\n')
        f.write('fLevels2,1\n')
        f.write('factorName3,N\n')
        f.write('fFirst3,12\n')
        f.write('fLast3,12\n')
        f.write('fLevels3,1')

        f.close()
        c = c + 1