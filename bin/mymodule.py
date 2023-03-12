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
    'grey' :        [0.7, 0.7, 0.7, 1.0],
    'greyTS' :      [0.9, 0.9, 0.9, 1.0],
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

def grey(T, R, P, S):
    m = (P == S)
    return m

def harmony(T, R, P, S):
    m = (T < R) & (R > P) & (P < S)
    return m

def deadlock(T, R, P, S):
    m = (T > R) & (R < P) & (P > S)
    return m

def prisoner(T, R, P, S):
    m = (T > R) & (R > P) & (P > S)
    return m

def snowdrift(T, R, P, S):
    m = (T > R) & (R > P) & (P < S)
    return m

def leader(T, R, P, S):
    m = (T > S) & (S > R) & (R > P)
    return m

def dilemma(T, R, P, S):
    m = prisoner(T, R, P, S) | snowdrift(T, R, P, S) | (deadlock(T, R, P, S) & (2.0*P < T + S)) | leader(T, R, P, S)
    return m

def TS(m, T, R, S):
    m = (m & (2.0*R < T + S))
    return m

def nodilemma(T, R, P, S):
    m = harmony(T, R, P, S) | (deadlock(T, R, P, S) & (2.0*P >= T + S)) | diagonal(T, R, P, S)
    return m

def diagonal(T, R, P, S):
    m = (T == R) & (R == P) & (P == S)
    return m

def greycolors(T, R, P, S, Z):
    m = grey(T, R, P, S)
    Z[m] = colormap['grey']
    Z[TS(m, T, R, S)] = colormap['greyTS']
    Z[diagonal(T, R, P, S)] = colormap['white']

def harmonycolors(T, R, P, S, Z):
    m = harmony(T, R, P, S)
    Z[m] = colormap['white']
    Z[diagonal(T, R, P, S)] = colormap['white']

def deadlockcolors(T, R, P, S, Z):
    m = deadlock(T, R, P, S)
    Z[m] = colormap['white']
    Z[TS(m, T, P, S)] = colormap['deadlockTS']
    Z[diagonal(T, R, P, S)] = colormap['white']
    pass

def prisonercolors(T, R, P, S, Z):
    m = prisoner(T, R, P, S)
    Z[m] = colormap['prisoner']
    Z[TS(m, T, R, S)] = colormap['prisonerTS']
    Z[diagonal(T, R, P, S)] = colormap['white']
    pass

def prisonerTScolors(T, R, P, S, Z):
    m = prisoner(T, R, P, S)
    Z[TS(m, T, R, S)] = colormap['prisonerTS']
    Z[diagonal(T, R, P, S)] = colormap['white']
    pass

def snowdriftcolors(T, R, P, S, Z):
    m = snowdrift(T, R, P, S)
    Z[m] = colormap['snowdrift']
    Z[TS(m, T, R, S)] = colormap['snowdriftTS']
    Z[diagonal(T, R, P, S)] = colormap['white']
    pass

def leadercolors(T, R, P, S, Z):
    m = leader(T, R, P, S)
    Z[m] = colormap['leader']
    Z[diagonal(T, R, P, S)] = colormap['white']

def nodilemmacolors(T, R, P, S):
    Z = np.full([*T.shape, 4], colormap['transparent'])
    m = nodilemma(T, R, P, S)
    Z[m] = colormap['white']
    return Z

def nodilemmacolorsg(T, R, P, S):
    Z = np.full([*T.shape, 4], colormap['transparent'])
    m = nodilemma(T, R, P, S)
    Z[m] = colormap['grey']
    return Z

def gamecolors(T, R, P, S):
    Z = np.full([*T.shape, 4], colormap['red'])
    greycolors(T, R, P, S, Z)
    harmonycolors(T, R, P, S, Z)
    deadlockcolors(T, R, P, S, Z)
    prisonercolors(T, R, P, S, Z)
    snowdriftcolors(T, R, P, S, Z)
    leadercolors(T, R, P, S, Z)
    return Z

def equilibrium(T, R, P, S, low, high, a2, weq):

    # Harmony
    m = (R > T) & (S > P) 
    a2[m] = high[m]
    weq[m] = R[m]

    # Deadlock or prisoner's dilemma
    m = (T >= R) & (P >= S) 
    a2[m] = low[m]
    weq[m] = P[m]

    # Snowdrift (chicken) or leader
    m = (T > R) & (S > P)
    weq[m] = (P[m] - S[m])/(R[m] - S[m] - T[m] + P[m])
    a2[m] = high[m]*weq[m] + low[m]*(1.0 - weq[m])
    weq[m] = (T[m] + S[m])*weq[m]*(1.0 - weq[m]) + R[m]*weq[m]*weq[m] + P[m]*(1.0 - weq[m])*(1.0 - weq[m])

    pass

def fitness(x, y, given, alpha, rho):
    q1 = (a2max - y)*R1/b
    q2 = y*R2*(1.0 - given) + x*R2*given
    w = q1*q2
    if np.isscalar(alpha):
        alpha = np.full(np.shape(q1), alpha)
    if np.isscalar(rho):
        rho = np.full(np.shape(q1), rho)
    m = (w > 0.0) & (rho == 0.0)
    w[m] = pow(q1[m], 1.0 - alpha[m])*pow(q2[m], alpha[m])
    m = ((w > 0.0) & (rho < 0.0)) | (rho > 0.0)
    w[m] = pow((1.0 - alpha[m])*pow(q1[m], rho[m]) + alpha[m]*pow(q2[m], rho[m]), 1.0/rho[m])
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

