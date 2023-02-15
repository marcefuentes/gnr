#! /usr/bin/env python

import numpy as np

R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0
Rq = R2/R1
b = a2max/a1max
alphamin = 0.1
alphamax = 0.9
logesmin = -5.0
logesmax = 5.0

colormap = {
    'white' :       [1.0, 1.0, 1.0, 1.0],
    'red' :         [1.0, 0.0, 0.0, 1.0],
    'harmonyTS' :   [1.0, 0.0, 0.0, 1.0],
    'snowdrift' :   [0.7, 0.0, 0.7, 1.0],
    'snowdriftTS' : [1.0, 0.3, 1.0, 1.0],
    'prisoner' :    [0.0, 0.1, 0.3, 1.0],
    'prisonerTS' :  [0.3, 0.5, 0.8, 1.0],
    'deadlockTS' :  [0.6, 1.0, 0.6, 1.0],
}

def fitness(x, y, given, alpha, rho):
    q1 = (a2max - y)*R1/b
    q2 = y*R2*(1.0 - given) + x*R2*given
    w = q1*q2
    mask = (w > 0.0) & (rho == 0.0)
    w[mask] = pow(q1[mask], 1.0 - alpha[mask])*pow(q2[mask], alpha[mask])
    mask = (w > 0.0) & (rho < 0.0)
    w[mask] = (1.0 - alpha[mask])*pow(q1[mask], rho[mask]) + alpha[mask]*pow(q2[mask], rho[mask])
    mask = (w > 0.0) & (rho < 0.0)
    w[mask] = pow(w[mask], 1.0/rho[mask])
    mask = (rho > 0.0)
    w[mask] = pow((1.0 - alpha[mask])*pow(q1[mask], rho[mask]) + alpha[mask]*pow(q2[mask], rho[mask]), 1.0/rho[mask])
    return w

def harmony(T, R, P, S):
    mask = (T < R) & (R > P) & (P < S) 
    return mask

def deadlock(T, R, P, S):
    mask = (T > R) & (R < P) & (P > S)
    return mask

def prisoner(T, R, P, S):
    mask = (T > R) & (R > P) & (P > S)
    return mask

def snowdrift(T, R, P, S):
    mask = (T >= R) & (R > S) & (S >= P)
    return mask

def TS(mask, T, R, S):
    mask = (mask & (2.0*R < T + S))
    return mask

def diagonal(R, P):
    mask = np.isclose(R, P)
    return mask

def gamecolors(T, R, P, S, Z):
    harmonycolors(T, R, P, S, Z)
    deadlockcolors(T, R, P, S, Z)
    prisonercolors(T, R, P, S, Z)
    snowdriftcolors(T, R, P, S, Z)
    pass

def harmonycolors(T, R, P, S, Z):
    mask = harmony(T, R, P, S)
    Z[mask] = colormap['white']
    Z[diagonal(R, P)] = colormap['white']

def deadlockcolors(T, R, P, S, Z):
    mask = deadlock(T, R, P, S)
    Z[mask] = colormap['white']
    Z[TS(mask, T, P, S)] = colormap['deadlockTS']
    Z[diagonal(R, P)] = colormap['white']
    pass

def prisonercolors(T, R, P, S, Z):
    mask = prisoner(T, R, P, S)
    Z[mask] = colormap['prisoner']
    Z[TS(mask, T, R, S)] = colormap['prisonerTS']
    Z[diagonal(R, P)] = colormap['white']
    pass

def prisonerTScolors(T, R, P, S, Z):
    mask = prisoner(T, R, P, S)
    Z[TS(mask, T, R, S)] = colormap['prisonerTS']
    Z[diagonal(R, P)] = colormap['white']
    pass

def snowdriftcolors(T, R, P, S, Z):
    mask = snowdrift(T, R, P, S)
    Z[mask] = colormap['snowdrift']
    Z[TS(mask, T, R, S)] = colormap['snowdriftTS']
    Z[diagonal(R, P)] = colormap['white']
    pass

def equilibrium(T, R, P, S, low, high, a2eq, weq):

    # Harmony
    mask = (R > T) & (S > P) 
    a2eq[mask] = high[mask]
    weq[mask] = R[mask]

    # Deadlock or prisoner's dilemma
    mask = (T >= R) & (P >= S) 
    a2eq[mask] = low[mask]
    weq[mask] = P[mask]

    # Snowdrift (chicken)
    mask = (T > R) & (S > P)
    weq[mask] = (P[mask] - S[mask])/(R[mask] - S[mask] - T[mask] + P[mask])
    a2eq[mask] = high[mask]*weq[mask] + low[mask]*(1.0 - weq[mask])
    weq[mask] = (T[mask] + S[mask])*weq[mask]*(1.0 - weq[mask]) + R[mask]*weq[mask]*weq[mask] + P[mask]*(1.0 - weq[mask])*(1.0 - weq[mask])

    pass

