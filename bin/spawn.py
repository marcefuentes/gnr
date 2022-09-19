#! /usr/bin/env python

ess = range(-5, 6)
givens = range(0, 11)
c = 101

for given in givens:
    givenstring = str('{:.1f}'.format(given/10.0))
    for es in ess:
        esstring = str(es)

        filename = str(c) + '.glo'
        f = open(filename, 'w')

        f.write('Seed,1\n')
        f.write('N,12\n')
        f.write('Runs,30\n')
        f.write('Time,20\n')
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
        f.write('fFirst1,' + givenstring + '\n')
        f.write('fLast1,' + givenstring + '\n')
        f.write('fLevels1,1\n')
        f.write('factorName2,ES\n')
        f.write('fFirst2,' + esstring + '\n')
        f.write('fLast2,' + esstring + '\n')
        f.write('fLevels2,1\n')
        f.write('factorName3,N\n')
        f.write('fFirst3,12\n')
        f.write('fLast3,12\n')
        f.write('fLevels3,1')

        f.close()
        c = c + 1
