#! /usr/bin/env python

from glob import glob
import os
import time

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score

import mymodule as my

start_time = time.perf_counter()
thisscript = os.path.basename(__file__)
filename = thisscript.split('.')[0]

folders = ['given100', 'given95', 'given50']
subfolders = ['none', 'p', 'r']

# Get data

def read_file(file, alltimes):
    df = pd.read_csv(file)
    if not alltimes:
        df = df.tail(1)
    return df

filelist = glob('given00/none/*.csv')
df = [read_file(file, False) for file in filelist]
df = pd.concat(df, ignore_index=True)
df = df.sort_values(['alpha', 'logES'])
high = df['a2Seenmean'].values
high = np.tile(high, (len(folders), 1)).flatten()

dfs = []
for folder in folders:
    df_list = []
    for subfolder in subfolders:
        filelist = glob(os.path.join(folder, subfolder, '*.csv'))
        df = [read_file(file, False) for file in filelist]
        df_concat = pd.concat(df, ignore_index=True)
        df_list.append(df_concat)
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

X = np.stack((T, R, P, S), axis=1)

CG = dfns[1]['ChooseGrainmean'].values
MG = dfns[2]['MimicGrainmean'].values
CG = 1.0 - CG
MG = 1.0 - MG

# Split data into training and testing sets
X_train, X_test, CG_train, CG_test = train_test_split(X, CG, test_size=0.2, random_state=42)

# Choose a linear regression model
model = LinearRegression()

# Train the model on the training data
model.fit(X_train, CG_train)

# Make predictions on the testing data
CG_pred = model.predict(X_test)

# Evaluate the performance of the model
mse = mean_squared_error(CG_test, CG_pred)
r2 = r2_score(CG_test, CG_pred)
print('Mean Squared Error:', mse)
print('R-squared:', r2)

# Print the regression equation
coef = model.coef_
intercept = model.intercept_
eqn = 'CG = {:.2f} + {:.2f}*T + {:.2f}*R + {:.2f}*P + {:.2f}*S'.format(intercept, coef[0], coef[1], coef[2], coef[3])
print('Regression Equation:', eqn)

# Split data into training and testing sets
X_train, X_test, MG_train, MG_test = train_test_split(X, MG, test_size=0.2, random_state=42)

# Choose a linear regression model
model = LinearRegression()

# Train the model on the training data
model.fit(X_train, MG_train)

# Make predictions on the testing data
MG_pred = model.predict(X_test)

# Evaluate the performance of the model
mse = mean_squared_error(MG_test, MG_pred)
r2 = r2_score(MG_test, MG_pred)
print('Mean Squared Error:', mse)
print('R-squared:', r2)

# Print the regression equation
coef = model.coef_
intercept = model.intercept_
eqn = 'MG = {:.2f} + {:.2f}*T + {:.2f}*R + {:.2f}*P + {:.2f}*S'.format(intercept, coef[0], coef[1], coef[2], coef[3])
print('Regression Equation:', eqn)

end_time = time.perf_counter()
print(f'\nTime elapsed: {(end_time - start_time):.2f} seconds')
