#!/usr/bin/python

import os
import pandas as pd

for dirs, files in os.walk():
    for subdir in dirs:
        print(subdir)
