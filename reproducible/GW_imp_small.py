#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 12:05:16 2024

@author: shuochieh
"""

#%%
# make sure python is operating at the directory where WI_util is in

import WI_util as WI
import numpy as np

#%%

p = 6

#%%

for pct in [2, 3, 4, 5]:
    Kalman = np.loadtxt("./real_data/small_series" + str(pct) + "_KS.csv", delimiter = ",", skiprows = 0)
    Kalman = Kalman.T

    n = Kalman.shape[0]
    d = Kalman.shape[1]

    raw = np.genfromtxt("./real_data/small_series" + str(pct) + "_masked.csv", delimiter = ",")
    raw = raw.T

    res_Kalman = np.zeros(shape = (n, d))
    kres_Kalman = np.zeros(shape = (n, d))

    alpha = n / (4 * p)

    for j in range(d):
        if j == 0:
            idx_obs = [[i for i in range(n) if ~np.isnan(raw[i, j])]]
        else:
            idx_obs.append([i for i in range(n) if ~np.isnan(raw[i, j])])

    res_Kalman = WI.WI_core_ordinary(np.copy(Kalman), 123, p, idx_obs, WI.solver_ordinary,
                                     alpha = alpha, WI_max_iter = 100)
    kres_Kalman = WI.kWI(np.copy(Kalman), WI.WI_core_ordinary, [231, 154, 77],
                         idx_obs = idx_obs, p = p, solver = WI.solver_ordinary, alpha = alpha,
                         WI_max_iter = 30, verbose = False)
    
    np.savetxt("./real_data/WI_Kalman_small" + str(pct) + ".csv", res_Kalman, delimiter = ",")
    np.savetxt("./real_data/kWI_Kalman_small" + str(pct) + ".csv", kres_Kalman, delimiter = ",")











    