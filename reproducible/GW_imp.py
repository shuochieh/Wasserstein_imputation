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
linimp = np.loadtxt("./real_data/lin_imp.csv", delimiter = ",", skiprows = 0)
Kalman = np.loadtxt("./real_data/KS_imp.csv", delimiter = ",", skiprows = 0)

n = linimp.shape[0]
d = linimp.shape[1]

#%%

raw = np.genfromtxt("./real_data/GW_clip.csv", delimiter = ",")

#%%
res_lin = np.zeros(shape = (n, d))
kres_lin = np.zeros(shape = (n, d))
res_Kalman = np.zeros(shape = (n, d))
kres_Kalman = np.zeros(shape = (n, d))

alpha = n / (4 * p)

for j in range(d):    
    idx_obs = [[i for i in range(n) if ~np.isnan(raw[j, i])]]
    
    # linear initialization for WI
    print(f"\n === Station {j}: WI from linear interpolation ===\n")
    x_obs = np.copy(linimp[:,j]).reshape(-1, 1)
    temp = WI.WI_core_ordinary(x_obs, 160, p, idx_obs, WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_lin[:,j] = temp.squeeze()
    
    # linear initialization for kWI
    print(f"\n === Station {j}: kWI from linear interpolation ===\n")
    x_obs = np.copy(linimp[:,j]).reshape(-1, 1)
    temp = WI.kWI(x_obs, WI.WI_core_ordinary, [160, 225, 110], p = p, idx_obs = idx_obs, solver = WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    kres_lin[:,j] = temp.squeeze()

    # Kalman initialization for WI
    print(f"\n === Station {j}: WI from Kalman smoothing ===\n")
    x_obs = np.copy(Kalman[:,j]).reshape(-1, 1)
    temp = WI.WI_core_ordinary(x_obs, 160, p, idx_obs, WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_Kalman[:,j] = temp.squeeze()
    
    # linear initialization for kWI
    print(f"\n === Station {j}: kWI from Kalman smoothing ===\n")
    x_obs = np.copy(Kalman[:,j]).reshape(-1, 1)
    temp = WI.kWI(x_obs, WI.WI_core_ordinary, [160, 225, 110], p = p, idx_obs = idx_obs, solver = WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    kres_Kalman[:,j] = temp.squeeze()
    
#%%

np.savetxt("./real_data/WI_lin.csv", res_lin, delimiter = ",")    
np.savetxt("./real_data/kWI_lin.csv", kres_lin, delimiter = ",")    
np.savetxt("./real_data/WI_Kalman.csv", res_Kalman, delimiter = ",")    
np.savetxt("./real_data/kWI_Kalman.csv", kres_Kalman, delimiter = ",")    











    