"""
Add IFN-gammm inhibition on virus infecting healthy cells.
Date: 2023.05.02
"""

import numpy as np
from scipy.integrate import odeint
from Equation import func
from Type_Characterization2 import Type_Characterization, WithinPhys
from E_Calculation import R0, E, E1, E2, E_kill, E_clear
import warnings

warnings.filterwarnings('error')


def model2(x0, time, Para, k_IFNg):
    # change k_infect
    IFNg = x0[29]
    parameters = np.copy(Para)
    parameters[25] = parameters[25]*(1-k_IFNg*IFNg/(50+IFNg))
    return func(x0, time, parameters)


def get_mode(Para):
    Time = np.arange(0, 80, 0.1)
    initial = [0.01, 0, Para[52] / Para[62], 0, 0, Para[53] / Para[65], 0, 0,
               Para[159], 0, 0, 0, 0, 0, 0,
               Para[54] / Para[70], 0, Para[160], 0, 0, 0, 0, 0, 0,
               Para[77] / Para[101], Para[82] / Para[102],
               Para[84] / Para[103],
               (Para[88] + Para[54] / Para[70] * Para[90]) / Para[104],
               Para[91] / Para[105],
               Para[95] / Para[106], 0, 0]

    k_IFNg = Para[-1]
    result = odeint(model2, initial, Time, args=(Para[:-1], k_IFNg))
    wp, constraints = WithinPhys(result)

    # compute other indexes
    sol = np.copy(result)
    v = sol[:, 0]
    v[v < 10**(-4)] = 10**(-4)
    sol[:, 0] = np.log10(v)

    R = R0(result, Para)
    e = E(result, Para)
    ei = E1(result, Para)
    ea = E2(result, Para)
    ek = E_kill(result, Para)
    ec = E_clear(result, Para)
    CD8T = np.sum(result[:, 17:21], axis=1)
    CD4T = np.sum(result[:, 8:14], axis=1) + result[:, 16]
    sol = np.vstack((sol.T, R, e, ei, ea, ek, ec, CD8T, CD4T)).T

    if wp:
        mode = Type_Characterization(result)
        return mode, sol, constraints
    else:
        return 100, sol, constraints
