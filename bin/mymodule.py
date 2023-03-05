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
cost = pow(2, -14)
deathrate = pow(2, -7)
repeats = 1.0/(1.0 - pow(1.0 - deathrate, 2.0))

colormap = {
    'transparent':  [1.0, 1.0, 1.0, 0.0],
    'white' :       [1.0, 1.0, 1.0, 1.0],
    'grey' :        [0.8, 0.8, 0.8, 1.0],
    'red' :         [1.0, 0.0, 0.0, 1.0],
    'harmonyTS' :   [1.0, 0.0, 0.0, 1.0],
    'snowdrift' :   [0.7, 0.0, 0.7, 1.0],
    'snowdriftTS' : [1.0, 0.3, 1.0, 1.0],
    'leader' :      [1.0, 0.3, 0.5, 1.0],
    'prisoner' :    [0.0, 0.1, 0.3, 1.0],
    'prisonerTS' :  [0.3, 0.5, 0.8, 1.0],
    'deadlockTS' :  [0.6, 1.0, 0.6, 1.0],
    'deadlock' :    [0.9, 1.0, 0.9, 1.0],
}

def borderline(T, R, P, S):
    mask = (T == R) | (R == P) | (P == S)
    return mask

def harmony(T, R, P, S):
    mask = ((R > T) & (T >= S) & (S > P)) | ((R >= S) & (S >= T) & (T >= P))
    return mask

def deadlock(T, R, P, S):
    mask = (T > R) & (R < P) & (P > S)
    return mask

def prisoner(T, R, P, S):
    mask = (T > R) & (R > P) & (P > S)
    return mask

def snowdrift(T, R, P, S):
    mask = (T > R) & (R > S) & (S > P)
    return mask

def leader(T, R, P, S):
    mask = (T > S) & (S > R) & (R > P)
    return mask

def dilemma(T, R, P, S):
    mask = prisoner(T, R, P, S) | snowdrift(T, R, P, S) | (deadlock(T, R, P, S) & (2.0*P < T + S)) | leader(T, R, P, S)
    return mask

def TS(mask, T, R, S):
    mask = (mask & (2.0*R < T + S))
    return mask

def nodilemma(T, R, P, S):
    mask = harmony(T, R, P, S) | (deadlock(T, R, P, S) & (2.0*P >= T + S)) | diagonal(T, R, P, S)
    return mask

def diagonal(T, R, P, S):
    mask = (T == R) & (R == P) & (P == S)
    return mask

def borderlinecolors(T, R, P, S, Z):
    mask = borderline(T, R, P, S)
    Z[mask] = colormap['grey']
    Z[diagonal(T, R, P, S)] = colormap['white']

def harmonycolors(T, R, P, S, Z):
    mask = harmony(T, R, P, S)
    Z[mask] = colormap['white']
    Z[diagonal(T, R, P, S)] = colormap['white']

def deadlockcolors(T, R, P, S, Z):
    mask = deadlock(T, R, P, S)
    Z[mask] = colormap['white']
    Z[TS(mask, T, P, S)] = colormap['deadlockTS']
    Z[diagonal(T, R, P, S)] = colormap['white']
    pass

def prisonercolors(T, R, P, S, Z):
    mask = prisoner(T, R, P, S)
    Z[mask] = colormap['prisoner']
    Z[TS(mask, T, R, S)] = colormap['prisonerTS']
    Z[diagonal(T, R, P, S)] = colormap['white']
    pass

def prisonerTScolors(T, R, P, S, Z):
    mask = prisoner(T, R, P, S)
    Z[TS(mask, T, R, S)] = colormap['prisonerTS']
    Z[diagonal(T, R, P, S)] = colormap['white']
    pass

def snowdriftcolors(T, R, P, S, Z):
    mask = snowdrift(T, R, P, S)
    Z[mask] = colormap['snowdrift']
    Z[TS(mask, T, R, S)] = colormap['snowdriftTS']
    Z[diagonal(T, R, P, S)] = colormap['white']
    pass

def leadercolors(T, R, P, S, Z):
    mask = leader(T, R, P, S)
    Z[mask] = colormap['leader']
    Z[diagonal(T, R, P, S)] = colormap['white']

def nodilemmacolors(T, R, P, S):
    Z = np.full([*T.shape, 4], colormap['transparent'])
    mask = nodilemma(T, R, P, S)
    Z[mask] = colormap['white']
    return Z

def nodilemmacolorsg(T, R, P, S):
    Z = np.full([*T.shape, 4], colormap['transparent'])
    mask = nodilemma(T, R, P, S)
    Z[mask] = colormap['grey']
    return Z

def gamecolors(T, R, P, S):
    Z = np.full([*T.shape, 4], colormap['red'])
    borderlinecolors(T, R, P, S, Z)
    harmonycolors(T, R, P, S, Z)
    deadlockcolors(T, R, P, S, Z)
    prisonercolors(T, R, P, S, Z)
    snowdriftcolors(T, R, P, S, Z)
    leadercolors(T, R, P, S, Z)
    return Z

def equilibrium(T, R, P, S, low, high, a2eq, weq):

    # Harmony
    mask = (R > T) & (S > P) 
    a2eq[mask] = high[mask]
    weq[mask] = R[mask]

    # Deadlock or prisoner's dilemma
    mask = (T >= R) & (P >= S) 
    a2eq[mask] = low[mask]
    weq[mask] = P[mask]

    # Snowdrift (chicken) or leader
    mask = (T > R) & (S > P)
    weq[mask] = (P[mask] - S[mask])/(R[mask] - S[mask] - T[mask] + P[mask])
    a2eq[mask] = high[mask]*weq[mask] + low[mask]*(1.0 - weq[mask])
    weq[mask] = (T[mask] + S[mask])*weq[mask]*(1.0 - weq[mask]) + R[mask]*weq[mask]*weq[mask] + P[mask]*(1.0 - weq[mask])*(1.0 - weq[mask])

    pass

def fitness(x, y, given, alpha, rho):
    q1 = (a2max - y)*R1/b
    q2 = y*R2*(1.0 - given) + x*R2*given
    w = q1*q2
    if np.isscalar(alpha):
        alpha = np.full([len(q1)], alpha)
    if np.isscalar(rho):
        rho = np.full([len(q1)], rho)
    mask = (w > 0.0) & (rho == 0.0)
    w[mask] = pow(q1[mask], 1.0 - alpha[mask])*pow(q2[mask], alpha[mask])
    mask = (w > 0.0) & (rho < 0.0)
    w[mask] = (1.0 - alpha[mask])*pow(q1[mask], rho[mask]) + alpha[mask]*pow(q2[mask], rho[mask])
    mask = (w > 0.0) & (rho < 0.0)
    w[mask] = pow(w[mask], 1.0/rho[mask])
    mask = (rho > 0.0)
    w[mask] = pow((1.0 - alpha[mask])*pow(q1[mask], rho[mask]) + alpha[mask]*pow(q2[mask], rho[mask]), 1.0/rho[mask])
    return w

def a2eq(given, alpha, rho):
    if given < 1.0:
        MRT = b*Rq*(1.0 - given)
        Q = Rq*pow(MRT*alpha/(1.0 - alpha), 1.0/(rho - 1.0))
        a2 = a2max/(1.0 + Q*b)
    else:
        a2 = alpha*0.0
    return a2

def indifference(q, w, alpha, rho):
    if rho == 0.0:
        q2 = pow(w/pow(q, 1.0 - alpha), 1.0/alpha)
    else:
        A = pow(w, rho)
        B = (1.0 - alpha)*pow(q, rho)
        if A <= B:
            if rho < 0.0:
                q2 = 1000.0
            else:
                q2 = -0.1
        else:
            q2 = pow((A - B)/alpha, 1.0/rho)
    return q2

