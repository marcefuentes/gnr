#!/usr/bin/python

import glob
import os
import pandas as pd

dirs = ["", "drift/"]

dfs = []

for d in dirs:
    files = glob.glob(os.path.join(d, "*.csv"))
    for f in files:
        df = pd.read_csv(f)
        gs = df.pop("GroupSize")
        df.insert(2, "GroupSize", gs)
        df.to_csv(f, index=False)
