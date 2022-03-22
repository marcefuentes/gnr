#!/usr/bin/python

import glob
import os
import sys
import pandas as pd

dirs = ['', 'drift/']

dfs = []

for d in dirs:
    files = glob.glob(os.path.join(d, '*.csv'))
    for f in files:
        df = pd.read_csv(f)
        df["GroupSize"] = sys.argv[1]
        df.to_csv(f, index=False)
