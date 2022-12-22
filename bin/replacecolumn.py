#!/usr/bin/python

import os
import pandas as pd

alpha = 0.5

for root, dirs, files in os.walk('.'):
    for filename in files:
        if '.csv' in filename:
            df = pd.read_csv(os.path.join(root, filename))
            newcol = [alpha]*len(df.index)
            df['alpha'] = newcol
            df.drop(['N'], inplace=True, axis=1)
            df.insert(2, 'alpha', df.pop('alpha'))
            df.to_csv(os.path.join(root, filename), index=False)
