#!/usr/bin/python

import sys
import pandas as pd

df = pd.read_csv(sys.argv[1])
df.drop(sys.argv[2], inplace=True, axis=1)
df[sys.argv[2]] = sys.argv[3]
df.to_csv(sys.argv[1], index=False)
