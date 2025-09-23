# Definitions of the supported input shapers
#
# Copyright (C) 2020-2021  Dmitry Butyugin <dmbutyugin@google.com>
#
# This file may be distributed under the terms of the GNU GPLv3 license.
import collections, math
import numpy as np # FLSUN Changes

SHAPER_VIBRATION_REDUCTION=20.
DEFAULT_DAMPING_RATIO = 0.1

InputShaperCfg = collections.namedtuple(
        'InputShaperCfg', ('name', 'init_func', 'min_freq'))

# Start FLSUN Changes
def step(shaper_freq, damping_ratio,t):
    wn = 2. * math.pi * shaper_freq
    df = math.sqrt(1. - damping_ratio**2)
    wd = wn * df
    Beta = math.acos(damping_ratio)
    e = math.exp(-wn * damping_ratio * t)
    sin1 = math.sin(wd * t + Beta)
    cos1 = math.cos(wd * t + Beta)
    s = t + (wd * e * cos1 + damping_ratio * wn * e * sin1)/((damping_ratio **2 * wn **2 + wd **2)*df)
    return s

def generate_matrix_C_zv(shaper_freq, damping_ratio,n,Tc):
    wn = 2. * math.pi * shaper_freq
    df = math.sqrt(1. - damping_ratio**2)
    wd = wn * df
    matrix_C = np.zeros((4,n))
    for i in range(1,n+1):
        matrix_C[0][i-1] = np.exp(wn * damping_ratio * (i-1) * Tc) * np.cos(wd * (i-1) * Tc)
        matrix_C[1][i-1] = np.exp(wn * damping_ratio * (i-1) * Tc) * np.sin(wd * (i-1) * Tc)
        matrix_C[2][i-1] = 1
        matrix_C[3][i-1] = i-1
    return matrix_C

def generate_matrix_C_zvd(shaper_freq, damping_ratio,n,Tc):
    wn = 2. * math.pi * shaper_freq
    df = math.sqrt(1. - damping_ratio**2)
    wd = wn * df
    matrix_C = np.zeros((6,n))
    for i in range(1,n+1):
        matrix_C[0][i-1] = np.exp(wn * damping_ratio * (i-1) * Tc) * np.cos(wd * (i-1) * Tc)
        matrix_C[1][i-1] = np.exp(wn * damping_ratio * (i-1) * Tc) * np.sin(wd * (i-1) * Tc)
        matrix_C[2][i-1] = 1
        matrix_C[3][i-1] = (i-1) * np.exp(wn * damping_ratio * (i-1) * Tc) * np.cos(wd * (i-1) * Tc)
        matrix_C[4][i-1] = (i-1) * np.exp(wn * damping_ratio * (i-1) * Tc) * np.sin(wd * (i-1) * Tc)
        matrix_C[5][i-1] = i-1
    return matrix_C

def generate_vector_b_zv(shaper_freq, damping_ratio,m,Tc):
    wn = 2. * math.pi * shaper_freq
    vector_b = np.array([0,0,1,m - 2 * damping_ratio / (wn * Tc)]).reshape(4, 1)
    return vector_b

def generate_vector_b_zvd(shaper_freq, damping_ratio,m,Tc):
    wn = 2. * math.pi * shaper_freq
    vector_b = np.array([0,0,1,0,0,m - 2 * damping_ratio / (wn * Tc)]).reshape(6, 1)
    return vector_b

def generate_matrix_Q(n,qi):
    matrix_Q = np.zeros((n,n))
    for i in range(n) :
        matrix_Q[i][i] = qi
    return matrix_Q

def generate_matrix_Psi(shaper_freq, damping_ratio,n,Tc):
    wn = 2. * math.pi * shaper_freq
    df = math.sqrt(1. - damping_ratio**2)
    wd = wn * df
    Beta = math.acos(damping_ratio)
    matrix_Psi = np.zeros((n,n))
    for i in range(1,n+1) :
        for j in range(1,n+1) :
            a = (i-1)*Tc
            b = (j-1)*Tc
            t_min = max(a,b)
            t_max = (n-1)*Tc
            t = np.linspace(t_min,t_max)
            e1 = np.exp(-wn * damping_ratio * (t-a))
            e2 = np.exp(-wn * damping_ratio * (t-b))
            sin1 = np.sin(wd * (t-a) + Beta)
            sin2 = np.sin(wd * (t-b) + Beta)
            ft = (1 - e1*sin1/df)*(1-e2*sin2/df)
            matrix_Psi[i-1][j-1] = np.trapz(ft, t)
    return matrix_Psi

def generate_matrix_H(shaper_freq, damping_ratio,n,Tc):
    wn = 2. * math.pi * shaper_freq
    df = math.sqrt(1. - damping_ratio**2)
    wd = wn * df
    Beta = math.acos(damping_ratio)
    matrix_H = np.zeros((n,n))
    for i in range(1,n+1) :
        for j in range(1,n+1) :
            a = (i-1)*Tc
            b = (j-1)*Tc
            t_min = (n-1)*Tc
            t_max = 2*(n-1)*Tc
            t = np.linspace(t_min,t_max)
            e1 = np.exp(-wn * damping_ratio * (t-a))
            e2 = np.exp(-wn * damping_ratio * (t-b))
            sin1 = np.sin(wd * (t-a) + Beta)
            sin2 = np.sin(wd * (t-b) + Beta)
            ft = (1 - e1*sin1/df)*(1-e2*sin2/df)
            matrix_H[i-1][j-1] = np.trapz(ft, t)
    return matrix_H

def generate_vector_theta(shaper_freq, damping_ratio,n,Tc,m):
    vector_theta = np.zeros((n,1))
    for i in range(1,n+1):
        t_min = max(m-i-1, 0) * Tc
        vector_theta[i-1][0] = step(shaper_freq, damping_ratio,(n-i-2)*Tc) - step(shaper_freq, damping_ratio,t_min)
    return vector_theta

def generate_vector_g(shaper_freq, damping_ratio,n,Tc):
    vector_g = np.zeros((n,1))
    for i in range(1,n+1):
        t_min = max(n-i-2, 0) * Tc
        vector_g[i-1][0] = step(shaper_freq, damping_ratio,(2*n-i-1)*Tc) - step(shaper_freq, damping_ratio,t_min)
    return vector_g
# End FLSUN Changes

def get_none_shaper():
    return ([], [])

def get_zv_shaper(shaper_freq, damping_ratio):
    df = math.sqrt(1. - damping_ratio**2)
    K = math.exp(-damping_ratio * math.pi / df)
    t_d = 1. / (shaper_freq * df)
    A = [1., K]
    T = [0., .5*t_d]
    return (A, T)

def get_zvd_shaper(shaper_freq, damping_ratio):
    df = math.sqrt(1. - damping_ratio**2)
    K = math.exp(-damping_ratio * math.pi / df)
    t_d = 1. / (shaper_freq * df)
    A = [1., 2.*K, K**2]
    T = [0., .5*t_d, t_d]
    return (A, T)

def get_mzv_shaper(shaper_freq, damping_ratio):
    df = math.sqrt(1. - damping_ratio**2)
    K = math.exp(-.75 * damping_ratio * math.pi / df)
    t_d = 1. / (shaper_freq * df)

    a1 = 1. - 1. / math.sqrt(2.)
    a2 = (math.sqrt(2.) - 1.) * K
    a3 = a1 * K * K

    A = [a1, a2, a3]
    T = [0., .375*t_d, .75*t_d]
    return (A, T)

def get_ei_shaper(shaper_freq, damping_ratio):
    v_tol = 1. / SHAPER_VIBRATION_REDUCTION # vibration tolerance
    df = math.sqrt(1. - damping_ratio**2)
    t_d = 1. / (shaper_freq * df)
    dr = damping_ratio

    a1 = (0.24968 + 0.24961 * v_tol) + (( 0.80008 + 1.23328 * v_tol) +
                                        ( 0.49599 + 3.17316 * v_tol) * dr) * dr
    a3 = (0.25149 + 0.21474 * v_tol) + ((-0.83249 + 1.41498 * v_tol) +
                                        ( 0.85181 - 4.90094 * v_tol) * dr) * dr
    a2 = 1. - a1 - a3

    t2 = 0.4999 + ((( 0.46159 + 8.57843 * v_tol) * v_tol) +
                   (((4.26169 - 108.644 * v_tol) * v_tol) +
                    ((1.75601 + 336.989 * v_tol) * v_tol) * dr) * dr) * dr

    A = [a1, a2, a3]
    T = [0., t2 * t_d, t_d]
    return (A, T)

def _get_shaper_from_expansion_coeffs(shaper_freq, damping_ratio, t, a):
    tau = 1. / shaper_freq
    T = []
    A = []
    n = len(a)
    k = len(a[0])
    for i in range(n):
        u = t[i][k-1]
        v = a[i][k-1]
        for j in range(k-1):
            u = u * damping_ratio + t[i][k-j-2]
            v = v * damping_ratio + a[i][k-j-2]
        T.append(u * tau)
        A.append(v)
    return (A, T)

def get_2hump_ei_shaper(shaper_freq, damping_ratio):
    t = [[0., 0., 0., 0.],
         [0.49890,  0.16270, -0.54262, 6.16180],
         [0.99748,  0.18382, -1.58270, 8.17120],
         [1.49920, -0.09297, -0.28338, 1.85710]]
    a = [[0.16054,  0.76699,  2.26560, -1.22750],
         [0.33911,  0.45081, -2.58080,  1.73650],
         [0.34089, -0.61533, -0.68765,  0.42261],
         [0.15997, -0.60246,  1.00280, -0.93145]]
    return _get_shaper_from_expansion_coeffs(shaper_freq, damping_ratio, t, a)

def get_3hump_ei_shaper(shaper_freq, damping_ratio):
    t = [[0., 0., 0., 0.],
         [0.49974,  0.23834,  0.44559, 12.4720],
         [0.99849,  0.29808, -2.36460, 23.3990],
         [1.49870,  0.10306, -2.01390, 17.0320],
         [1.99960, -0.28231,  0.61536, 5.40450]]
    a = [[0.11275,  0.76632,  3.29160 -1.44380],
         [0.23698,  0.61164, -2.57850,  4.85220],
         [0.30008, -0.19062, -2.14560,  0.13744],
         [0.23775, -0.73297,  0.46885, -2.08650],
         [0.11244, -0.45439,  0.96382, -1.46000]]
    return _get_shaper_from_expansion_coeffs(shaper_freq, damping_ratio, t, a)

# Start FLSUN Changes
def get_zv_zvd_shaper(shaper_freq, damping_ratio):
    df = math.sqrt(1. - damping_ratio**2)
    K = math.exp(-damping_ratio * math.pi / df)
    t_d = 1. / (shaper_freq * df)
    A = [1., 3.*K, 3*K**2, K**3]
    T = [0., .5*t_d, t_d, 1.5*t_d]
    return (A, T)

def get_zero_zv_shaper(shaper_freq, damping_ratio):
    Tc = 0.001
    m = 15
    qi = 0.2
    k1 = 1
    k2 = 10
    b = generate_vector_b_zv(shaper_freq, damping_ratio,m,Tc)
    for n in range(6,700):
        A = np.zeros(n)
        T = np.zeros(n)
        for i in range(n):
            T[i] = i * Tc - m * Tc
        C = generate_matrix_C_zv(shaper_freq, damping_ratio,n,Tc)
        P1 = generate_matrix_Q(n,qi) + k1 * generate_matrix_Psi(shaper_freq, damping_ratio,n,Tc) + k2 * generate_matrix_H(shaper_freq, damping_ratio,n,Tc)
        P2 = k1 * generate_vector_theta(shaper_freq, damping_ratio,n,Tc,m) + k2 * generate_vector_g(shaper_freq, damping_ratio,n,Tc)
        P1_ = np.linalg.inv(P1)
        CT = np.transpose(C)
        CP1_CT = np.dot(np.dot(C,P1_),CT)
        CP1_CT_ = np.linalg.inv(CP1_CT)
        CP1_P2 = np.dot(np.dot(C,P1_),P2)
        CTCP1_CT_ = np.dot(CT,CP1_CT_)
        A1 = np.dot(CTCP1_CT_,(CP1_P2-b))
        A2 = P2 - A1
        A = np.dot(P1_,A2)
        max_A = np.max(A)
        min_A = np.min(A)
        diff = np.diff(np.transpose(A))
        max_diff = np.max(diff)
        min_diff = np.min(diff)
        if (max_A < 0.2) and (min_A > -0.2) and (max_diff < 0.15) and (min_diff > - 0.15):
            break
    A = A.flatten().tolist()
    T = T.flatten().tolist()
    return (A, T)
# End FLSUN Changes

# min_freq for each shaper is chosen to have projected max_accel ~= 1500
INPUT_SHAPERS = [
    InputShaperCfg('zv', get_zv_shaper, min_freq=21.),
    InputShaperCfg('mzv', get_mzv_shaper, min_freq=23.),
    InputShaperCfg('zvd', get_zvd_shaper, min_freq=29.),
    InputShaperCfg('ei', get_ei_shaper, min_freq=29.),
    InputShaperCfg('2hump_ei', get_2hump_ei_shaper, min_freq=39.),
    InputShaperCfg('3hump_ei', get_3hump_ei_shaper, min_freq=48.),
    # Start FLSUN Changes
    InputShaperCfg('zv_zvd', get_zv_zvd_shaper, min_freq=21.),
    InputShaperCfg('zero_zv', get_zero_zv_shaper, min_freq=21.),
    # End FLSUN Changes
]
