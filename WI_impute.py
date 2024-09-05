#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 22:23:10 2024

@author: shuochieh
"""

#%%
# make sure python is operating at the directory where WI_util is in

import WI_util as WI
import numpy as np
import matplotlib.pyplot as plt


#%%

p = 6 # dimension of marginal distribution
n = 3000 # sample size (for our simulation)

#%% 

# First, conduct imputation for univariate time series
# AR(1) imputation

dta = np.loadtxt("./sim_data/AR/AR_miss1.csv", delimiter = ",")
benchmarks = np.loadtxt("./sim_data/AR/benchmarks_miss1.csv", delimiter = ",", skiprows = 1)
dta_lin = benchmarks[:,:1000]
dta_Kalman = benchmarks[:,2000:3000]

res_ord_lin = np.zeros(shape = (n, 1000))
res_ord_Kalman = np.zeros(shape = (n, 1000))
res_kWI_lin = np.zeros(shape = (n, 1000))
res_kWI_Kalman = np.zeros(shape = (n, 1000))

for sim in range(1000):
    idx_obs = [[i for i in range(n) if ~np.isnan(dta[i, sim])]]
    alpha = n / (4 * p) # learning rate
    
    print(f"\n === Simulation {sim}: WI from linear interpolation ===\n")
    x_obs = dta_lin[:,sim].reshape(-1, 1)
    temp = WI.WI_core_ordinary(x_obs, 1200, p, idx_obs, WI.solver_ordinary, alpha = alpha, WI_max_iter = 20, verbose = True)
    res_ord_lin[:,sim] = temp.squeeze()
    
    print("\n === WI from Kalman interpolation ===\n")
    
    
    


#%%
benchmarks = np.loadtxt("./sim_data/AR/benchmarks_miss2.csv", delimiter = ",", skiprows = 1)
dta_lin = benchmarks[:,:1000]
x_obs = dta_lin[:,0].reshape(-1,1)
plt.plot(np.linspace(0, n - 1, n), x_obs, 'c-', linewidth = 0.8)
plt.show()

temp = WI.WI_core_ordinary(x_obs, 1500, 2, idx_obs, WI.solver_ordinary, alpha = alpha, WI_max_iter = 20, verbose = True)

plt.plot(np.linspace(0, n - 1, n), temp, 'c-', linewidth = 0.8)
plt.show()




















