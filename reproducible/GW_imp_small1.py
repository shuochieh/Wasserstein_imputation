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

p = 6

#%%

for pct in [1, 2, 3, 4, 5, 6]:

    init = np.loadtxt("./real_data/small_series" + str(pct) + "_LPlinimp.csv", delimiter = ",", skiprows = 0)
    init = init.T

    n = init.shape[0]
    d = init.shape[1]

    n1 = int(np.floor(n / 3))
    n2 = int(np.floor(2 * n / 3))

    raw = np.genfromtxt("./real_data/small_series" + str(pct) + "_LPresd.csv", delimiter = ",")
    raw = raw.T

    res = np.zeros(shape = (n, d))
    kres = np.zeros(shape = (n, d))

    alpha = n / (4 * 3 * p)

    for j in range(d):
        print("\n\n Series" + str(j))

        #idx_obs = [[i for i in range(n) if ~np.isnan(raw[i, j])]]

        idx_obs = [[i for i in range(n) if ~np.isnan(raw[i, j])]]

        # Impute by three sections
        # temp1 = WI.WI_core_ordinary(np.copy(init[:n1,j]).reshape(-1, 1), int(np.floor(0.4 * n1)), p, idx_obs1, WI.solver_ordinary,
        #                             alpha = alpha, WI_max_iter = 100, verbose = False)
        # temp2 = WI.WI_core_ordinary(np.copy(init[n1:n2,j]).reshape(-1, 1), int(np.floor(0.4 * n1)), p, idx_obs2, WI.solver_ordinary,
        #                             alpha = alpha, WI_max_iter = 100, verbose = False)
        # temp3 = WI.WI_core_ordinary(np.copy(init[n2:,j]).reshape(-1, 1), int(np.floor(0.4 * n1)), p, idx_obs3, WI.solver_ordinary,
        #                             alpha = alpha, WI_max_iter = 100, verbose = False)
        # res[:,j] = np.concatenate([temp1.reshape(-1), temp2.reshape(-1), temp3.reshape(-1)])
        temp = WI.WI_core_ordinary(np.copy(init[:,j]).reshape(-1, 1), int(np.floor(0.4 * n)), p, idx_obs, WI.solver_ordinary,
                                   alpha = alpha, WI_max_iter = 100, verbose = False)
        res[:,j] = temp.reshape(-1)

        temp = WI.kWI(np.copy(init[:,j]).reshape(-1, 1), WI.WI_core_ordinary, [92, 215, 154],
                      idx_obs = idx_obs, p = p, solver = WI.solver_ordinary,
                      alpha = alpha, WI_max_iter = 30, verbose = True)
        kres[:,j] = temp.reshape(-1)

        # Codes below have not been finished
        #temp1 = WI.kWI(np.copy(init[:n1,j]).reshape(-1,1), WI.WI_core_ordinary, [92, 215, 154],
        #              idx_obs = idx_obs, p = p, solver = WI.solver_ordinary,
        #              alpha = alpha, WI_max_iter = 30, verbose = False)
        #kres[:,j] = temp.reshape(-1)

    np.savetxt("./real_data/WI_small" + str(pct) + ".csv", res, delimiter = ",")
    np.savetxt("./real_data/kWI_small" + str(pct) + ".csv", kres, delimiter = ",")

