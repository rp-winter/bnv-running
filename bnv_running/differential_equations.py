import numpy as np
import math

# DE for the couplings
def couplings_de(yu, yd, g1, g2, g3, NJKL):
    Nf = 6 # flavor number
    N = NJKL[0]
    J = NJKL[1]
    K = NJKL[2]
    L = NJKL[3]
    
    # beta functions for the couplings
    sum_g1_sq = np.sum(g1**2)
    sum_g2_sq = np.sum(g2**2)
    sum_g3_sq = np.sum(g3**2)

    sum_yu_sq = np.sum(yu[1:]**2)
    sum_yd_sq = np.sum(yd[1:]**2)

    Yu = np.array([
        [1.3e-5, 0, 0],
        [0, 7.3e-3, 0],
        [0, 0, 1.0]
    ])

    Yd = np.array([
        [2.7e-5, 0, 0],
        [0, 5.5e-4, 0],
        [0, 0, 0.024]
    ])

    T = np.trace(3*np.matmul(Yd.T, Yd) + 3*np.matmul(Yu.T, Yu))

    beta1 = (-(2/3)*Nf - (N+1)/10)*g1**2
    beta2 = g2**2 * (16/3 - (2/3)*Nf - (N+1)/6) + 2*sum_g2_sq
    beta3 = g3**2 * (8 - 2*Nf/3) + 3*sum_g3_sq

    beta_yu = (3/2)*(Yu[0, 0]**2 - Yd[0, 0]**2) + T + ((9/2)*sum_yu_sq + (3/2)*sum_yd_sq) - ( (17/20)*sum_g1_sq + (9/4)*sum_g2_sq + 8*sum_g3_sq )
    beta_yd = (3/2)*(Yd[0, 0]**2 - Yu[0, 0]**2) + T + ((9/2)*sum_yd_sq + (3/2)*sum_yu_sq) - ( (1/4)*sum_g1_sq + (9/4)*sum_g2_sq + 8*sum_g3_sq )

    # return the derivatives
    # yu0, yu1, ..., yuN, yd0, yd1, ..., ydN, g10, g11, ..., g1J, g20, g21, ..., g2K, g30, g31, ..., g3L
    dcouplings_dt = []
    sixteen_pi_sq = 16*math.pi**2
    for i in range(N+1):
        dcouplings_dt.append(beta_yu*yu[i]/sixteen_pi_sq)
    for i in range(N+1):
        dcouplings_dt.append(beta_yd*yd[i]/sixteen_pi_sq)
    dcouplings_dt.extend(-beta1*g1/sixteen_pi_sq)
    dcouplings_dt.extend(-beta2*g2/sixteen_pi_sq)
    dcouplings_dt.extend(-beta3*g3/sixteen_pi_sq)
    return dcouplings_dt

# DE for the Wilson coefficients
def wc_de(t, wc, yu, yd, g1, g2, g3, NJKL, gamma_matrix_fn, wc_de_kwargs):
    Gamma = gamma_matrix_fn(t, yu, yd, g1, g2, g3, NJKL, **wc_de_kwargs)
    return np.dot(np.transpose(Gamma), wc)/(16*math.pi**2)
 