#! /usr/bin/env python

colormap = {
    'default' :     [1.0, 1.0, 1.0, 1.0],
    'deadlock' :    [0.93, 0.96, 0.96, 1.0],
    'deadlockTS' :  [0.6, 1.0, 1.0, 1.0],
    'harmony' :     [1.0, 0.96, 0.98, 1.0],
    'prisoner' :    [0.2, 0.1, 0.4, 1.0],
    'prisonerTS' :  [0.0, 1.0, 1.0, 1.0],
    'snowdrift' :   [1.0, 0.0, 1.0, 1.0],
    'snowdriftTS' : [1.0, 0.0, 1.0, 1.0],
    'white' :       [1.0, 1.0, 1.0, 1.0],
    'black' :       [0.0, 0.0, 0.0, 1.0]
}

R1 = 2.0
R2 = 2.0
a1max = 1.0
a2max = 1.0
a2low = 0.25
a2high = 0.75
Rq = R2/R1
b = a2max/a1max

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

def gametypes(T, R, P, S, Z):

    # Harmony
    mask = ((R > T) & (T > S) & (S > P)) | ((R > S) & (S > T) & (S > P)) 
    Z[mask] = colormap['harmony']

    # Deadlock
    mask = ((T >= P) & (P > R) & (R >= S)) 
    Z[mask] = colormap['deadlock']

    mask = (mask & (2.0*P < T + S))
    Z[mask] = colormap['deadlockTS']

    # Prisoner's dilemma
    mask = (T > R) & (R >= P) & (P > S)
    Z[mask] = colormap['prisoner']

    mask = (mask & (2.0*R < T + S))
    Z[mask] = colormap['prisonerTS']

    # Snowdrift (chicken)
    mask = (T > R) & (R > S) & (S > P)
    Z[mask] = colormap['snowdrift']

    mask = (mask & (2.0*R < T + S))
    Z[mask] = colormap['snowdriftTS']

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

