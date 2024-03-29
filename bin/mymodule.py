import numpy as np
import pandas as pd

RA = 2.0
RB = 2.0
aAmax = 1.0
a2max = 1.0
wmax = 2.0 # For RA = 2.0, RB = 2.0, aAmax = 1.0, a2max = 1.0,
           # given = 1.0 and any values of alpha and rho
Rq = RB/RA
b = a2max/aAmax
alphamin = 0.1
alphamax = 0.9
logesmin = -5.0
logesmax = 5.0
cost = pow(2, -14)
deathrate = pow(2, -7)
repeats = 1.0/(1.0 - pow(1.0 - deathrate, 2.0))

colormap = {
    "transparent":  (1.0, 1.0, 1.0, 0.0),
    "red":          (1.0, 0.0, 0.0, 1.0),
    "white" :       (1.0, 1.0, 1.0, 1.0),
    "grey" :        (0.7, 0.7, 0.7, 1.0),
    "greyTS" :      (0.9, 0.9, 0.9, 1.0),
    "harmony" :     (1.0, 0.9, 1.0, 1.0),
    "harmonyTS" :   (1.0, 0.8, 0.9, 1.0),
    "snowdrift" :   (0.7, 0.0, 0.7, 1.0),
    "snowdriftTS" : (1.0, 0.3, 1.0, 1.0),
    "drift" :       (0.7, 0.5, 0.7, 1.0),
    "driftTS" :     (1.0, 0.7, 1.0, 1.0),
    "leader" :      (1.0, 0.3, 0.5, 1.0),
    "prisoner" :    (0.0, 0.1, 0.3, 1.0),
    "prisonerTS" :  (0.3, 0.5, 0.8, 1.0),
    "deadlockTS" :  (0.6, 1.0, 0.6, 1.0),
    "deadlock" :    (0.9, 1.0, 0.9, 1.0),
}

def harmony(T, R, P, S):
    m = ((T < R) & (R > P) & (P <= S)) | (R == P)
    return m

def deadlock(T, R, P, S):
    m = (T > R) & (R < P) & (P > S)
    return m

def prisoner(T, R, P, S):
    m = ((T > R) & (R > P) & (P >= S)) 
    return m

def snowdrift(T, R, P, S):
    m = (T > R) & (R > P) & (P < S)
    return m

def drift(T, R, P, S):
    m = (T == R) & (R > P) & (P == S)
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

def harmonycolors(T, R, P, S, Z):
    m = harmony(T, R, P, S)
    Z[m] = colormap["harmony"]
    Z[TS(m, T, R, S)] = colormap["harmonyTS"]
    Z[diagonal(T, R, P, S)] = colormap["white"]
    return Z

def deadlockcolors(T, R, P, S, Z):
    m = deadlock(T, R, P, S)
    Z[m] = colormap["deadlock"]
    Z[TS(m, T, P, S)] = colormap["deadlockTS"]
    Z[diagonal(T, R, P, S)] = colormap["white"]
    return Z

def prisonercolors(T, R, P, S, Z):
    m = prisoner(T, R, P, S)
    Z[m] = colormap["prisoner"]
    Z[TS(m, T, R, S)] = colormap["prisonerTS"]
    Z[diagonal(T, R, P, S)] = colormap["white"]
    return Z

def prisonerTScolors(T, R, P, S, Z):
    m = prisoner(T, R, P, S)
    Z[TS(m, T, R, S)] = colormap["prisonerTS"]
    Z[diagonal(T, R, P, S)] = colormap["white"]
    return Z

def snowdriftcolors(T, R, P, S, Z):
    m = snowdrift(T, R, P, S)
    Z[m] = colormap["snowdrift"]
    Z[TS(m, T, R, S)] = colormap["snowdriftTS"]
    Z[diagonal(T, R, P, S)] = colormap["white"]
    return Z

def driftcolors(T, R, P, S, Z):
    m = drift(T, R, P, S)
    Z[m] = colormap["drift"]
    Z[TS(m, T, R, S)] = colormap["driftTS"]
    Z[diagonal(T, R, P, S)] = colormap["white"]
    return Z

def leadercolors(T, R, P, S, Z):
    m = leader(T, R, P, S)
    Z[m] = colormap["leader"]
    Z[diagonal(T, R, P, S)] = colormap["white"]
    return Z

def nodilemmacolors(T, R, P, S):
    Z = np.full([*T.shape, 4], colormap["transparent"])
    m = nodilemma(T, R, P, S)
    Z[m] = colormap["white"]
    return Z

def nodilemmacolorsg(T, R, P, S):
    Z = np.full([*T.shape, 4], colormap["transparent"])
    m = nodilemma(T, R, P, S)
    Z[m] = colormap["grey"]
    return Z

def gamecolors(T, R, P, S):
    Z = np.full([*T.shape, 4], colormap["red"])
    Z = harmonycolors(T, R, P, S, Z)
    Z = deadlockcolors(T, R, P, S, Z)
    Z = prisonercolors(T, R, P, S, Z)
    Z = snowdriftcolors(T, R, P, S, Z)
    Z = driftcolors(T, R, P, S, Z)
    Z = leadercolors(T, R, P, S, Z)
    return Z

def equilibrium(T, R, P, S, low, high, a2, weq):

    # Harmony
    m = (R > T) & (S > P) 
    a2[m] = high[m]
    weq[m] = R[m]

    # Deadlock or prisoner"s dilemma
    m = (T >= R) & (P >= S) 
    a2[m] = low[m]
    weq[m] = P[m]

    # Snowdrift (chicken) or leader
    m = (T > R) & (S > P)
    weq[m] = (P[m] - S[m])/(R[m] - S[m] - T[m] + P[m])
    a2[m] = high[m]*weq[m] + low[m]*(1.0 - weq[m])
    weq[m] = (T[m] + S[m])*weq[m]*(1.0 - weq[m]) + R[m]*weq[m]*weq[m] + P[m]*(1.0 - weq[m])*(1.0 - weq[m])

    pass

def eqw(T, R, P, S):

    weq = np.full(T.shape, 0.0)

    # Harmony
    m = (R > T) & (S > P) 
    weq[m] = R[m]

    # Deadlock or prisoner"s dilemma
    m = (T >= R) & (P >= S) 
    weq[m] = P[m]

    # Snowdrift (chicken) or leader
    m = (T > R) & (S > P)
    weq[m] = (P[m] - S[m])/(R[m] - S[m] - T[m] + P[m])
    weq[m] = (T[m] + S[m])*weq[m]*(1.0 - weq[m]) + R[m]*weq[m]*weq[m] + P[m]*(1.0 - weq[m])*(1.0 - weq[m])

    return weq

def fitness(x, y, given, alpha, rho):
    qA = (a2max - y)*RA/b
    qB = y*RB*(1.0 - given) + x*RB*given
    w = qA*qB
    if not isinstance(qA, np.ndarray):
        qA = np.array([qA])
    if not isinstance(qB, np.ndarray):
        qB = np.array([qB])
    if not isinstance(w, np.ndarray):
        w = np.array([w])
    if not isinstance(alpha, np.ndarray):
        alpha = np.full(w.shape, alpha)
    if not isinstance(rho, np.ndarray):
        rho = np.full(w.shape, rho)
    m = (w > 0.0) & (rho == 0.0)
    w[m] = pow(qA[m], 1.0 - alpha[m])*pow(qB[m], alpha[m])
    m = ((w > 0.0) & (rho < 0.0)) | (rho > 0.0)
    w[m] = pow((1.0 - alpha[m])*pow(qA[m], rho[m]) + alpha[m]*pow(qB[m], rho[m]), 1.0/rho[m])
    return w

def a2eq(given, alpha, rho):
    if given < 1.0:
        MRT = b*Rq*(1.0 - given)
        Q = Rq*pow(MRT*alpha/(1.0 - alpha), 1.0/(rho - 1.0))
        a2 = a2max/(1.0 + Q*b)
    else:
        a2 = alpha*0.0
    return a2

def indifference(qs, w, alpha, rho):
    qB = np.full(qs.shape, 1000.0)
    for i, q in enumerate(qs):
        if rho == 0.0:
            qB[i] = pow(w/pow(q, 1.0 - alpha), 1.0/alpha)
        else:
            A = pow(w, rho)
            B = (1.0 - alpha)*pow(q, rho)
            if A <= B:
                if rho < 0.0:
                    qB[i] = 1000.0
                else:
                    qB[i] = -0.1
            else:
                qB[i] = pow((A - B)/alpha, 1.0/rho)
    return qB

def read_files(filelist, alltimes):
    df_list = [None] * len(filelist)
    for i, file in enumerate(filelist):
        df = pd.read_csv(file)
        if not alltimes:
            df = df.tail(1)
        df_list[i] = df
    dfc = pd.concat(df_list, ignore_index=True)
    return dfc

def get_Z(t, df, trait):
    m = df.Time == t
    Z = pd.pivot(df.loc[m],
                 values=trait,
                 index="alpha",
                 columns="logES")
    Z = Z.sort_index(axis=0, ascending=False)
    Z = Z.to_numpy()
    return Z

def get_Zd(t, df, alpha, loges, trait):
    m = (df.Time == t) & (df.alpha == alpha) & (df.logES == loges)
    Z = pd.pivot(df.loc[m],
                 values=trait,
                 index="a2high",
                 columns="a2low")
    Z = Z.sort_index(axis=0, ascending=False)
    Z = Z.to_numpy()
    Z = np.vstack((Z, np.full(Z.shape[1], np.nan)))
    Z = np.hstack((Z, np.full((Z.shape[0], 1), np.nan)))
    return Z

