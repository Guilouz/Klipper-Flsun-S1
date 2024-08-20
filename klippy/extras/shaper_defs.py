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
    K = math.exp(-damping_ratio * math.pi / df)
    t_d = 1. / (shaper_freq * df)

    a1 = .25 * (1. + v_tol)
    a2 = .5 * (1. - v_tol) * K
    a3 = a1 * K * K

    A = [a1, a2, a3]
    T = [0., .5*t_d, t_d]
    return (A, T)

def get_2hump_ei_shaper(shaper_freq, damping_ratio):
    v_tol = 1. / SHAPER_VIBRATION_REDUCTION # vibration tolerance
    df = math.sqrt(1. - damping_ratio**2)
    K = math.exp(-damping_ratio * math.pi / df)
    t_d = 1. / (shaper_freq * df)

    V2 = v_tol**2
    X = pow(V2 * (math.sqrt(1. - V2) + 1.), 1./3.)
    a1 = (3.*X*X + 2.*X + 3.*V2) / (16.*X)
    a2 = (.5 - a1) * K
    a3 = a2 * K
    a4 = a1 * K * K * K

    A = [a1, a2, a3, a4]
    T = [0., .5*t_d, t_d, 1.5*t_d]
    return (A, T)

def get_3hump_ei_shaper(shaper_freq, damping_ratio):
    v_tol = 1. / SHAPER_VIBRATION_REDUCTION # vibration tolerance
    df = math.sqrt(1. - damping_ratio**2)
    K = math.exp(-damping_ratio * math.pi / df)
    t_d = 1. / (shaper_freq * df)

    K2 = K*K
    a1 = 0.0625 * (1. + 3. * v_tol + 2. * math.sqrt(2. * (v_tol + 1.) * v_tol))
    a2 = 0.25 * (1. - v_tol) * K
    a3 = (0.5 * (1. + v_tol) - 2. * a1) * K2
    a4 = a2 * K2
    a5 = a1 * K2 * K2

    A = [a1, a2, a3, a4, a5]
    T = [0., .5*t_d, t_d, 1.5*t_d, 2.*t_d]
    return (A, T)

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
