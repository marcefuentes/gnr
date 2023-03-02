#!/usr/bin/python

import numpy as np
import os
import pandas as pd

oldcolname = 'ES'
newcolname = 'logES'
#newcolvalue = 0.95

for root, dirs, files in os.walk('.'):
    for filename in files:
        if ('csv' in filename) or ('frq' in filename):
            #print(filename)
            df = pd.read_csv(os.path.join(root, filename))
            #newcolvalues = [newcolvalue]*len(df.index)
            #df[newcolname] = newcolvalues
            if oldcolname in df:
                df[newcolname] = np.log2(df[oldcolname])
                df.drop(oldcolname, axis=1, inplace=True)
                df.insert(2, newcolname, df.pop(newcolname))
                df.to_csv(os.path.join(root, filename), float_format='%.6f', index=False)