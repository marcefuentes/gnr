#! /usr/bin/env python

from glob import glob
import os
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.ensemble import RandomForestRegressor
from sklearn.tree import plot_tree

import mymodule as my

start_time = time.perf_counter()
this_file = os.path.basename(__file__)
filename = this_file.split('.')[0]

folders = ['given100', 'given095', 'given050']
subfolders = ['none', 'p', 'r']

# Data

filelist = glob('given000/none/*.csv')
df = my.read_files(filelist, False)
df = df.sort_values(['alpha', 'logES'])
high = df['a2Seenmean'].values
high = np.tile(high, (len(folders), 1)).flatten()

dfs = []
for folder in folders:
    df_list = []
    for subfolder in subfolders:
        filelist = glob(os.path.join(folder, subfolder, '*.csv'))
        df_list.append(my.read_files(filelist, False)
    dfs.append(df_list)

dfns = []
for i, subfolder in enumerate(subfolders):
    df_concat = pd.concat([dfs[j][i] for j in range(len(folders))])
    dfns.append(df_concat)

for df in dfns:
    df = df.sort_values(['Given', 'alpha', 'logES'])

df = dfns[0]
low = df['a2Seenmean'].values
given = df['Given'].values
alpha = df['alpha'].values
loges = df['logES'].values
rho = 1.0 - 1.0/pow(2.0, loges)

T = my.fitness(high, low, given, alpha, rho)
R = my.fitness(high, high, given, alpha, rho)
P = my.fitness(low, low, given, alpha, rho)
S = my.fitness(low, high, given, alpha, rho)

dl = my.deadlock(T, R, P, S)
h = my.harmony(T, R, P, S)
pd = my.prisoner(T, R, P, S)
sd = my.snowdrift(T, R, P, S)
dl = dl.astype(int)
h = h.astype(int)
pd = pd.astype(int)
sd = sd.astype(int)

X = np.stack((T - R, R - P, R - S, P - S),axis=1)

CG = dfns[1]['ChooseGrainmean'].values
MG = dfns[2]['MimicGrainmean'].values
CG = 1.0 - CG
MG = 1.0 - MG

mask = pd 
X = X[mask]
CG = CG[mask]

# Split data into training and testing sets
X_train, X_test, CG_train, CG_test = train_test_split(X, CG, test_size=0.2, random_state=42)

# Choose a linear regression model
#model = LinearRegression()

# Create a random forest regressor object
rf = RandomForestRegressor(n_estimators=100, random_state=42)

# Fit the regressor to the training data
rf.fit(X_train, CG_train)

# Make predictions on the test data
CG_pred = rf.predict(X_test)

# Evaluate the model using mean squared error
mse = mean_squared_error(CG_test, CG_pred)
print(f"Mean Squared Error: {mse:.2f}")

# Get the individual decision trees in the random forest
trees = rf.estimators_

# Print the first tree
print(trees[0])

# Get the feature importances from the trained model
importances = rf.feature_importances_

# Print the feature importances
print(importances)

tree = rf.estimators_[0]

plt.figure(figsize=(20, 10))
plot_tree(tree, filled=True, fontsize=12)

plt.savefig('tree.png')

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
