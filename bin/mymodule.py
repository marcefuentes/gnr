#! /usr/bin/env python

colormap = {
    'default' :     [0.0, 1.0, 0.0, 1.0],
    'deadlock' :    [0.95, 0.95, 0.95, 1.0],
    'deadlockTS' :  [1.0, 0.7, 0.0, 1.0],
    'harmony' :     [1.0, 1.0, 0.0, 1.0],
    'harmonyTS' :   [0.8, 0.8, 0.0, 1.0],
    'prisoner' :    [0.5, 0.5, 0.5, 1.0],
    'prisonerTS' :  [1.0, 0.3, 0.0, 1.0],
    'snowdrift' :   [0.0, 1.0, 1.0, 1.0],
    'snowdriftTS' : [0.0, 0.5, 1.0, 1.0],
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

def gametypes(T, R, P, S, low, high, Z, TS, a2eq, xeq, weq):

    # Harmony
    mask = ((R > T) & (T > S) & (S > P)) | ((R > S) & (S > T) & (S > P)) 
    Z[mask] = colormap['harmony']
    a2eq[mask] = high[mask]
    weq[mask] = R[mask]

    mask = (mask & (2.0*R <= T + S))
    Z[mask] = colormap['harmonyTS']

    # Deadlock
    mask = ((T > P) & (P > R) & (R > S)) 
    Z[mask] = colormap['deadlock']
    a2eq[mask] = low[mask]
    weq[mask] = P[mask]

    mask = (mask & (2.0*P <= T + S))
    Z[mask] = colormap['deadlockTS']

    # Prisoner's dilemma
    mask = (T > R) & (R > P) & (P > S)
    Z[mask] = colormap['prisoner']
    TS[mask] = 1.0 + T[mask] + S[mask] - 2.0*R[mask]
    a2eq[mask] = low[mask]
    weq[mask] = P[mask]

    mask = (mask & (2.0*R <= T + S))
    Z[mask] = colormap['prisonerTS']
    TS[mask] = 1.0 + T[mask] + S[mask] - 2.0*R[mask]

    # Snowdrift or chicken
    mask = (T >= R) & (P <= S)
    Z[mask] = colormap['snowdrift']
    xeq[mask] = (P[mask] - S[mask])/(R[mask] - S[mask] - T[mask] + P[mask])
    a2eq[mask] = high[mask]*xeq[mask] + low[mask]*(1.0 - xeq[mask])
    weq[mask] = (T[mask] + S[mask])*xeq[mask]*(1.0 - xeq[mask]) + R[mask]*xeq[mask]*xeq[mask] + P[mask]*(1.0 - xeq[mask])*(1.0 - xeq[mask])

    mask = (mask & (2.0*R <= T + S))
    Z[mask] = colormap['snowdriftTS']

    pass

