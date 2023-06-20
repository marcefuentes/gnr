#!/usr/bin/python

import os
import pandas as pd

newcolname = "a2low"
newcolvalue = 0.25

for root, dirs, files in os.walk("."):
    for filename in files:
        if ("csv" in filename) or ("frq" in filename):
            #print(filename)
            df = pd.read_csv(os.path.join(root, filename))
            newcolvalues = [newcolvalue]*len(df.index)
            df[newcolname] = newcolvalues
            df.insert(3, newcolname, df.pop(newcolname))
            df.to_csv(os.path.join(root, filename), float_format="%.6f", index=False)
