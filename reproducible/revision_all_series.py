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

init = np.loadtxt("./real_data/all_series/all_series_LPlinimp.csv", delimiter = ",", skiprows = 0)
init = init.T

n = init.shape[0]
d = init.shape[1]

raw = np.genfromtxt("./real_data/all_series/all_series_LPresd.csv", delimiter = ",")
raw = raw.T

res = np.zeros(shape = (n, d))
kres = np.zeros(shape = (n, d))

alpha = n / (4 * 3 * p)

for j in range(d):
    print("\n\n Series" + str(j))

    idx_obs = [[i for i in range(n) if ~np.isnan(raw[i, j])]]

    # Univariate impute
    temp = WI.WI_core_ordinary(np.copy(init[:,j]).reshape(-1, 1), int(np.floor(0.4 * n)), p, idx_obs, WI.solver_ordinary,
                               alpha = alpha, WI_max_iter = 100, verbose = False)
    res[:,j] = temp.reshape(-1)

    temp = WI.kWI(np.copy(init[:,j]).reshape(-1, 1), WI.WI_core_ordinary, [92, 215, 154],
                  idx_obs = idx_obs, p = p, solver = WI.solver_ordinary,
                  alpha = alpha, WI_max_iter = 30, verbose = True)
    kres[:,j] = temp.reshape(-1)

#%%
# multivariate impute
idx_obs = []
for j in range(d):
    temp_index = [i for i in range(n)]
    idx_obs.append(temp_index)

mres = WI.WI_core_ordinary(np.copy(init), int(np.floor(0.4 * n)), p, idx_obs, WI.solver_ordinary,
                           alpha = alpha, WI_max_iter = 100, verbose = False)

kmres = WI.kWI(np.copy(init), WI.WI_core_ordinary, [92, 215, 154], idx_obs = idx_obs, p = p, solver = WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, verbose = False)

#%%

np.savetxt("./real_data/all_series/WI_all_series.csv", res, delimiter = ",")
np.savetxt("./real_data/all_series/kWI_all_series.csv", kres, delimiter = ",")
np.savetxt("./real_data/all_series/WI_all_series_multi.csv", mres, delimiter = ",")
np.savetxt("./real_data/all_series/kWI_all_series_multi.csv", kmres, delimiter = ",")

