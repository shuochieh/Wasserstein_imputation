#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 23:22:55 2024

@author: shuochieh
"""

#%%
# make sure python is operating at the directory where WI_util is in

import WI_util as WI
import numpy as np
import matplotlib.pyplot as plt


#%%

p = 6 # dimension of marginal distribution
n = 1000 # sample size (for our simulation)

#%% 

# AR(1) imputation: Missing pattern I

dta = np.loadtxt("./sim_data/AR/AR_miss1.csv", delimiter = ",")
benchmarks = np.loadtxt("./sim_data/AR/benchmarks_miss1.csv", delimiter = ",", skiprows = 1)
dta_Kalman = benchmarks[:,2000:3000]

res_ord_Kalman = np.zeros(shape = (n, 1000))
res_kWI_Kalman = np.zeros(shape = (n, 1000))

for sim in range(1000):
    idx_obs = [[i for i in range(n) if ~np.isnan(dta[i, sim])]]
    alpha = n / (4 * p) # learning rate
    
    print(f"\n === Simulation {sim}: WI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim]).reshape(-1, 1)
    temp = WI.WI_core_ordinary(x_obs, 400, p, idx_obs, WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_ord_Kalman[:,sim] = temp.squeeze()
        
    print("\n === kWI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim]).reshape(-1, 1)
    temp = WI.kWI(x_obs, WI.WI_core_ordinary, [500, 250, 750], p = p, idx_obs = idx_obs, solver = WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_kWI_Kalman[:,sim] = temp.squeeze()
        
    
np.savetxt("./sim_data/AR/AR_WI_Kalman_miss1.csv", np.column_stack((res_ord_Kalman, res_kWI_Kalman)), delimiter = ",")    

# AR(1) imputation: Missing pattern II

dta = np.loadtxt("./sim_data/AR/AR_miss2.csv", delimiter = ",")
benchmarks = np.loadtxt("./sim_data/AR/benchmarks_miss2.csv", delimiter = ",", skiprows = 1)
dta_Kalman = benchmarks[:,2000:3000]

res_ord_Kalman = np.zeros(shape = (n, 1000))
res_kWI_Kalman = np.zeros(shape = (n, 1000))

for sim in range(1000):
    idx_obs = [[i for i in range(n) if ~np.isnan(dta[i, sim])]]
    alpha = n / (4 * p) # learning rate
    
    print(f"\n === Simulation {sim}: WI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim]).reshape(-1, 1)
    temp = WI.WI_core_ordinary(x_obs, 400, p, idx_obs, WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_ord_Kalman[:,sim] = temp.squeeze()
    
    print("\n === kWI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim]).reshape(-1, 1)
    temp = WI.kWI(x_obs, WI.WI_core_ordinary, [500, 250, 750], p = p, idx_obs = idx_obs, solver = WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_kWI_Kalman[:,sim] = temp.squeeze()
    
    
np.savetxt("./sim_data/AR/AR_WI_Kalman_miss2.csv", np.column_stack((res_ord_Kalman, res_kWI_Kalman)), delimiter = ",")    


# ARMA(1,1) imputation: Missing pattern I

dta = np.loadtxt("./sim_data/ARMA/ARMA_miss1.csv", delimiter = ",")
benchmarks = np.loadtxt("./sim_data/ARMA/benchmarks_miss1.csv", delimiter = ",", skiprows = 1)
dta_Kalman = benchmarks[:,2000:3000]

res_ord_Kalman = np.zeros(shape = (n, 1000))
res_kWI_Kalman = np.zeros(shape = (n, 1000))

for sim in range(1000):
    idx_obs = [[i for i in range(n) if ~np.isnan(dta[i, sim])]]
    alpha = n / (4 * p) # learning rate
    
    print(f"\n === Simulation {sim}: WI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim]).reshape(-1, 1)
    temp = WI.WI_core_ordinary(x_obs, 400, p, idx_obs, WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_ord_Kalman[:,sim] = temp.squeeze()
    
    print("\n === kWI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim]).reshape(-1, 1)
    temp = WI.kWI(x_obs, WI.WI_core_ordinary, [500, 250, 750], p = p, idx_obs = idx_obs, solver = WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_kWI_Kalman[:,sim] = temp.squeeze()
        
    
np.savetxt("./sim_data/ARMA/ARMA_WI_Kalman_miss1.csv", np.column_stack((res_ord_Kalman, res_kWI_Kalman)), delimiter = ",")    

# ARMA(1,1) imputation: Missing pattern II

dta = np.loadtxt("./sim_data/ARMA/ARMA_miss2.csv", delimiter = ",")
benchmarks = np.loadtxt("./sim_data/ARMA/benchmarks_miss2.csv", delimiter = ",", skiprows = 1)
dta_Kalman = benchmarks[:,2000:3000]

res_ord_Kalman = np.zeros(shape = (n, 1000))
res_kWI_Kalman = np.zeros(shape = (n, 1000))

for sim in range(1000):
    idx_obs = [[i for i in range(n) if ~np.isnan(dta[i, sim])]]
    alpha = n / (4 * p) # learning rate
    
    print(f"\n === Simulation {sim}: WI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim]).reshape(-1, 1)
    temp = WI.WI_core_ordinary(x_obs, 400, p, idx_obs, WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_ord_Kalman[:,sim] = temp.squeeze()
        
    print("\n === kWI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim]).reshape(-1, 1)
    temp = WI.kWI(x_obs, WI.WI_core_ordinary, [500, 250, 750], p = p, idx_obs = idx_obs, solver = WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_kWI_Kalman[:,sim] = temp.squeeze()
        
    
np.savetxt("./sim_data/ARMA/ARMA_WI_Kalman_miss2.csv", np.column_stack((res_ord_Kalman, res_kWI_Kalman)), delimiter = ",")    

# TAR imputation: Missing pattern I

dta = np.loadtxt("./sim_data/TAR/TAR_miss1.csv", delimiter = ",")
benchmarks = np.loadtxt("./sim_data/TAR/benchmarks_miss1.csv", delimiter = ",", skiprows = 1)
dta_Kalman = benchmarks[:,2000:3000]

res_ord_Kalman = np.zeros(shape = (n, 1000))
res_kWI_Kalman = np.zeros(shape = (n, 1000))

for sim in range(1000):
    idx_obs = [[i for i in range(n) if ~np.isnan(dta[i, sim])]]
    alpha = n / (4 * p) # learning rate
    
    print(f"\n === Simulation {sim}: WI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim]).reshape(-1, 1)
    temp = WI.WI_core_ordinary(x_obs, 400, p, idx_obs, WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_ord_Kalman[:,sim] = temp.squeeze()
        
    print("\n === kWI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim]).reshape(-1, 1)
    temp = WI.kWI(x_obs, WI.WI_core_ordinary, [500, 250, 750], p = p, idx_obs = idx_obs, solver = WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_kWI_Kalman[:,sim] = temp.squeeze()
        
    
np.savetxt("./sim_data/TAR/TAR_WI_Kalman_miss1.csv", np.column_stack((res_ord_Kalman, res_kWI_Kalman)), delimiter = ",")    

# TAR imputation: Missing pattern II

dta = np.loadtxt("./sim_data/TAR/TAR_miss2.csv", delimiter = ",")
benchmarks = np.loadtxt("./sim_data/TAR/benchmarks_miss2.csv", delimiter = ",", skiprows = 1)
dta_Kalman = benchmarks[:,2000:3000]

res_ord_Kalman = np.zeros(shape = (n, 1000))
res_kWI_Kalman = np.zeros(shape = (n, 1000))

for sim in range(1000):
    idx_obs = [[i for i in range(n) if ~np.isnan(dta[i, sim])]]
    alpha = n / (4 * p) # learning rate
    
    print(f"\n === Simulation {sim}: WI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim]).reshape(-1, 1)
    temp = WI.WI_core_ordinary(x_obs, 400, p, idx_obs, WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_ord_Kalman[:,sim] = temp.squeeze()
        
    print("\n === kWI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim]).reshape(-1, 1)
    temp = WI.kWI(x_obs, WI.WI_core_ordinary, [500, 250, 750], p = p, idx_obs = idx_obs, solver = WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_kWI_Kalman[:,sim] = temp.squeeze()
        
    
np.savetxt("./sim_data/TAR/TAR_WI_Kalman_miss2.csv", np.column_stack((res_ord_Kalman, res_kWI_Kalman)), delimiter = ",")    
    
# GARCH imputation: Missing pattern I

dta = np.loadtxt("./sim_data/GARCH/GARCH_miss1.csv", delimiter = ",")
benchmarks = np.loadtxt("./sim_data/GARCH/benchmarks_miss1.csv", delimiter = ",", skiprows = 1)
dta_Kalman = benchmarks[:,2000:3000]

res_ord_Kalman = np.zeros(shape = (n, 1000))
res_kWI_Kalman = np.zeros(shape = (n, 1000))

for sim in range(1000):
    idx_obs = [[i for i in range(n) if ~np.isnan(dta[i, sim])]]
    alpha = n / (4 * p) # learning rate
    
    print(f"\n === Simulation {sim}: WI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim]).reshape(-1, 1)
    temp = WI.WI_core_ordinary(x_obs, 400, p, idx_obs, WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_ord_Kalman[:,sim] = temp.squeeze()
        
    print("\n === kWI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim]).reshape(-1, 1)
    temp = WI.kWI(x_obs, WI.WI_core_ordinary, [500, 250, 750], p = p, idx_obs = idx_obs, solver = WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_kWI_Kalman[:,sim] = temp.squeeze()
        
    
np.savetxt("./sim_data/GARCH/GARCH_WI_Kalman_miss1.csv", np.column_stack((res_ord_Kalman, res_kWI_Kalman)), delimiter = ",")    

# GARCH imputation: Missing pattern II

dta = np.loadtxt("./sim_data/GARCH/GARCH_miss2.csv", delimiter = ",")
benchmarks = np.loadtxt("./sim_data/GARCH/benchmarks_miss2.csv", delimiter = ",", skiprows = 1)
dta_Kalman = benchmarks[:,2000:3000]

res_ord_Kalman = np.zeros(shape = (n, 1000))
res_kWI_Kalman = np.zeros(shape = (n, 1000))

for sim in range(1000):
    idx_obs = [[i for i in range(n) if ~np.isnan(dta[i, sim])]]
    alpha = n / (4 * p) # learning rate
    
    print(f"\n === Simulation {sim}: WI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim]).reshape(-1, 1)
    temp = WI.WI_core_ordinary(x_obs, 400, p, idx_obs, WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_ord_Kalman[:,sim] = temp.squeeze()
        
    print("\n === kWI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim]).reshape(-1, 1)
    temp = WI.kWI(x_obs, WI.WI_core_ordinary, [500, 250, 750], p = p, idx_obs = idx_obs, solver = WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_kWI_Kalman[:,sim] = temp.squeeze()
        
    
np.savetxt("./sim_data/GARCH/GARCH_WI_Kalman_miss2.csv", np.column_stack((res_ord_Kalman, res_kWI_Kalman)), delimiter = ",")    

# Cyc imputation: Missing pattern I

dta = np.loadtxt("./sim_data/Cyc/Cyc_miss1.csv", delimiter = ",")
benchmarks = np.loadtxt("./sim_data/Cyc/benchmarks_miss1.csv", delimiter = ",", skiprows = 1)
dta_Kalman = benchmarks[:,2000:3000]

res_ord_Kalman = np.zeros(shape = (n, 1000))
res_kWI_Kalman = np.zeros(shape = (n, 1000))

for sim in range(1000):
    idx_obs = [[i for i in range(n) if ~np.isnan(dta[i, sim])]]
    alpha = n / (4 * p) # learning rate
    
    print(f"\n === Simulation {sim}: WI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim]).reshape(-1, 1)
    temp = WI.WI_core_ordinary(x_obs, 400, p, idx_obs, WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_ord_Kalman[:,sim] = temp.squeeze()
    
    print("\n === kWI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim]).reshape(-1, 1)
    temp = WI.kWI(x_obs, WI.WI_core_ordinary, [500, 250, 750], p = p, idx_obs = idx_obs, solver = WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_kWI_Kalman[:,sim] = temp.squeeze()
        
    
np.savetxt("./sim_data/Cyc/Cyc_WI_Kalman_miss1.csv", np.column_stack((res_ord_Kalman, res_kWI_Kalman)), delimiter = ",")    

# Cyc imputation: Missing pattern II

dta = np.loadtxt("./sim_data/Cyc/Cyc_miss2.csv", delimiter = ",")
benchmarks = np.loadtxt("./sim_data/Cyc/benchmarks_miss2.csv", delimiter = ",", skiprows = 1)
dta_Kalman = benchmarks[:,2000:3000]

res_ord_Kalman = np.zeros(shape = (n, 1000))
res_kWI_Kalman = np.zeros(shape = (n, 1000))

for sim in range(1000):
    idx_obs = [[i for i in range(n) if ~np.isnan(dta[i, sim])]]
    alpha = n / (4 * p) # learning rate
    
    print(f"\n === Simulation {sim}: WI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim]).reshape(-1, 1)
    temp = WI.WI_core_ordinary(x_obs, 400, p, idx_obs, WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_ord_Kalman[:,sim] = temp.squeeze()
        
    print("\n === kWI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim]).reshape(-1, 1)
    temp = WI.kWI(x_obs, WI.WI_core_ordinary, [500, 250, 750], p = p, idx_obs = idx_obs, solver = WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    res_kWI_Kalman[:,sim] = temp.squeeze()
    
    
np.savetxt("./sim_data/Cyc/Cyc_WI_Kalman_miss2.csv", np.column_stack((res_ord_Kalman, res_kWI_Kalman)), delimiter = ",")    

# ALR imputation: Missing pattern I

dta = np.loadtxt("./sim_data/ALR/ALR_miss1.csv", delimiter = ",")
benchmarks = np.loadtxt("./sim_data/ALR/benchmarks_miss1.csv", delimiter = ",", skiprows = 1)
dta_Kalman = benchmarks[:,6000:9000]

res_ord_Kalman = np.zeros(shape = (n, 3000))
res_kWI_Kalman = np.zeros(shape = (n, 3000))


for sim in range(1000):
    sub_dta = dta[:n,(3 * sim):(3 * (sim + 1))]
    idx_obs = [[i for i in range(n) if ~np.isnan(sub_dta[i, 0])]] 
    
    # Constructing linear equations that delineates the admissible set
    K = np.zeros((len(idx_obs[0]) * 3 + (n - len(idx_obs[0])), 3 * n))
    b = np.zeros((len(idx_obs[0]) * 3 + (n - len(idx_obs[0])), 1))
    counter = 0
    for i in range(n):
        if i in idx_obs[0]:
            for j in range(3):
                K[counter, i + j * n] = 1
                b[counter, 0] = sub_dta[i,j]
                counter += 1
        else:
            for j in range(3):
                K[counter, i + (j * n)] = 1
            b[counter, 0] = 1
            counter += 1            
    
    print(f"\n === Simulation {sim}: WI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,(sim * 3):((sim + 1) * 3)])
    temp = WI.WI_core_exact(x_obs, 400, p, K, b, Lambda = 1e-4, WI_max_iter = 10, WI_tol = 1e-4, verbose = True)
    res_ord_Kalman[:,(sim * 3):((sim + 1) * 3)] = temp
        
    print("\n === kWI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,(sim * 3):((sim + 1) * 3)])
    temp = WI.kWI(x_obs, WI.WI_core_exact, [500, 250, 750], p = p, K = K, b = b, Lambda = 1e-5, WI_max_iter = 10, WI_tol = 1e-4, verbose = True)
    res_kWI_Kalman[:,(sim * 3):((sim + 1) * 3)] = temp
        
    
np.savetxt("./sim_data/ALR/ALR_WI_Kalman_miss1.csv", np.column_stack((res_ord_Kalman, res_kWI_Kalman)), delimiter = ",")    

# ALR imputation: Missing pattern II

dta = np.loadtxt("./sim_data/ALR/ALR_miss2.csv", delimiter = ",")
benchmarks = np.loadtxt("./sim_data/ALR/benchmarks_miss2.csv", delimiter = ",", skiprows = 1)
dta_Kalman = benchmarks[:,6000:9000]

res_ord_Kalman = np.zeros(shape = (n, 3000))
res_kWI_Kalman = np.zeros(shape = (n, 3000))


for sim in range(1000):
    sub_dta = dta[:n,(3 * sim):(3 * (sim + 1))]
    idx_obs = [[i for i in range(n) if ~np.isnan(sub_dta[i, 0])]] 
    
    # Constructing linear equations that delineates the admissible set
    K = np.zeros((len(idx_obs[0]) * 3 + (n - len(idx_obs[0])), 3 * n))
    b = np.zeros((len(idx_obs[0]) * 3 + (n - len(idx_obs[0])), 1))
    counter = 0
    for i in range(n):
        if i in idx_obs[0]:
            for j in range(3):
                K[counter, i + j * n] = 1
                b[counter, 0] = sub_dta[i,j]
                counter += 1
        else:
            for j in range(3):
                K[counter, i + (j * n)] = 1
            b[counter, 0] = 1
            counter += 1            
    
    print(f"\n === Simulation {sim}: WI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,(sim * 3):((sim + 1) * 3)])
    temp = WI.WI_core_exact(x_obs, 400, p, K, b, Lambda = 1e-4, WI_max_iter = 10, WI_tol = 1e-4, verbose = True)
    res_ord_Kalman[:,(sim * 3):((sim + 1) * 3)] = temp
        
    print("\n === kWI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,(sim * 3):((sim + 1) * 3)])
    temp = WI.kWI(x_obs, WI.WI_core_exact, [500, 250, 750], p = p, K = K, b = b, Lambda = 1e-5, WI_max_iter = 10, WI_tol = 1e-4, verbose = True)
    res_kWI_Kalman[:,(sim * 3):((sim + 1) * 3)] = temp
        
    
np.savetxt("./sim_data/ALR/ALR_WI_Kalman_miss2.csv", np.column_stack((res_ord_Kalman, res_kWI_Kalman)), delimiter = ",")    

# NLVAR imputation: Missing pattern I

dta = np.loadtxt("./sim_data/NLVAR/NLVAR_miss1.csv", delimiter = ",")
benchmarks = np.loadtxt("./sim_data/NLVAR/benchmarks_miss1.csv", delimiter = ",", skiprows = 1)
dta_Kalman = benchmarks[:,4000:6000]

res_ord_Kalman = np.zeros(shape = (n, 2000))
res_kWI_Kalman = np.zeros(shape = (n, 2000))


for sim in range(1000):
    sub_dta = dta[:n,(2 * sim):(2 * (sim + 1))]
    idx_obs = [[i for i in range(n) if ~np.isnan(sub_dta[i, 0])]]
    idx_obs = [idx_obs, idx_obs]
    
    alpha = n / (p * 4)
    print(f"\n === Simulation {sim}: WI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,(sim * 2):((sim + 1) * 2)])
    temp = WI.WI_core_ordinary(x_obs, 400, p, idx_obs, WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)    
    res_ord_Kalman[:,(sim * 2):((sim + 1) * 2)] = temp
        
    print("\n === kWI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,(sim * 2):((sim + 1) * 2)])
    temp = WI.kWI(x_obs, WI.WI_core_ordinary, [500, 250, 750], idx_obs = idx_obs, p = p, solver = WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)    
    res_kWI_Kalman[:,(sim * 2):((sim + 1) * 2)] = temp
        
    
np.savetxt("./sim_data/NLVAR/NLVAR_WI_Kalman_miss1.csv", np.column_stack((res_ord_Kalman, res_kWI_Kalman)), delimiter = ",")    

# NLVAR imputation: Missing pattern II

dta = np.loadtxt("./sim_data/NLVAR/NLVAR_miss2.csv", delimiter = ",")
benchmarks = np.loadtxt("./sim_data/NLVAR/benchmarks_miss2.csv", delimiter = ",", skiprows = 1)
dta_Kalman = benchmarks[:,4000:6000]

res_ord_Kalman = np.zeros(shape = (n, 2000))
res_kWI_Kalman = np.zeros(shape = (n, 2000))


for sim in range(1000):
    sub_dta = dta[:n,(2 * sim):(2 * (sim + 1))]
    idx_obs = [[i for i in range(n) if ~np.isnan(sub_dta[i, 0])]]
    idx_obs = [idx_obs, idx_obs]
    
    alpha = n / (p * 4)
    print(f"\n === Simulation {sim}: WI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,(sim * 2):((sim + 1) * 2)])
    temp = WI.WI_core_ordinary(x_obs, 400, p, idx_obs, WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)    
    res_ord_Kalman[:,(sim * 2):((sim + 1) * 2)] = temp
        
    print("\n === kWI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,(sim * 2):((sim + 1) * 2)])
    temp = WI.kWI(x_obs, WI.WI_core_ordinary, [500, 250, 750], idx_obs = idx_obs, p = p, solver = WI.solver_ordinary, alpha = alpha, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)    
    res_kWI_Kalman[:,(sim * 2):((sim + 1) * 2)] = temp
        
    
np.savetxt("./sim_data/NLVAR/NLVAR_WI_Kalman_miss2.csv", np.column_stack((res_ord_Kalman, res_kWI_Kalman)), delimiter = ",")    

# ARI imputation: Missing pattern I

dta = np.loadtxt("./sim_data/ARI/ARI_miss1.csv", delimiter = ",")
benchmarks = np.loadtxt("./sim_data/ARI/benchmarks_miss1.csv", delimiter = ",", skiprows = 1)
dta_Kalman = benchmarks[:,2000:3000]

res_ord_Kalman = np.zeros(shape = (n, 1000))
res_kWI_Kalman = np.zeros(shape = (n, 1000))


for sim in range(1000):
    sub_dta = dta[:n,sim]
    idx_obs = [[i for i in range(n) if ~np.isnan(sub_dta[i])]]
        
    K = np.zeros((len(idx_obs[0]) - 1, n - 1))
    b = np.zeros((len(idx_obs[0]) - 1, 1))
    counter = 0
    for i in range(1, n):
        if i in idx_obs[0]:
            K[counter,:i] = 1
            b[counter,0] = sub_dta[i] - sub_dta[0]
            counter += 1        
    
    print(f"\n === Simulation {sim}: WI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim])
    x_obs = x_obs[1:] - x_obs[:(n - 1)]
    x_obs = x_obs.reshape(-1, 1)
    temp = WI.WI_core_exact(x_obs, 400, p, K, b, Lambda = 1e-4, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    temp = np.cumsum(np.insert(temp.squeeze(), 0, sub_dta[0]))
    res_ord_Kalman[:,sim] = temp
        
    print("\n === kWI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim])
    x_obs = x_obs[1:] - x_obs[:(n - 1)]
    x_obs = x_obs.reshape(-1, 1)
    temp = WI.kWI(x_obs, WI.WI_core_exact, [500, 750, 250], p = p, K = K, b = b, Lambda = 1e-4, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    temp = np.cumsum(np.insert(temp.squeeze(), 0, sub_dta[0]))
    res_kWI_Kalman[:,sim] = temp
        
    
np.savetxt("./sim_data/ARI/ARI_WI_Kalman_miss1.csv", np.column_stack((res_ord_Kalman, res_kWI_Kalman)), delimiter = ",")    

# ARI imputation: Missing pattern II

dta = np.loadtxt("./sim_data/ARI/ARI_miss2.csv", delimiter = ",")
benchmarks = np.loadtxt("./sim_data/ARI/benchmarks_miss2.csv", delimiter = ",", skiprows = 1)
dta_Kalman = benchmarks[:,2000:3000]

res_ord_Kalman = np.zeros(shape = (n, 1000))
res_kWI_Kalman = np.zeros(shape = (n, 1000))


for sim in range(1000):
    sub_dta = dta[:n,sim]
    idx_obs = [[i for i in range(n) if ~np.isnan(sub_dta[i])]]
        
    K = np.zeros((len(idx_obs[0]) - 1, n - 1))
    b = np.zeros((len(idx_obs[0]) - 1, 1))
    counter = 0
    for i in range(1, n):
        if i in idx_obs[0]:
            K[counter,:i] = 1
            b[counter,0] = sub_dta[i] - sub_dta[0]
            counter += 1        
    
    print(f"\n === Simulation {sim}: WI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim])
    x_obs = x_obs[1:] - x_obs[:(n - 1)]
    x_obs = x_obs.reshape(-1, 1)
    temp = WI.WI_core_exact(x_obs, 400, p, K, b, Lambda = 1e-4, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    temp = np.cumsum(np.insert(temp.squeeze(), 0, sub_dta[0]))
    res_ord_Kalman[:,sim] = temp
        
    print("\n === kWI from linear interpolation ===\n")
    x_obs = np.copy(dta_Kalman[:,sim])
    x_obs = x_obs[1:] - x_obs[:(n - 1)]
    x_obs = x_obs.reshape(-1, 1)
    temp = WI.kWI(x_obs, WI.WI_core_exact, [500, 750, 250], p = p, K = K, b = b, Lambda = 1e-4, WI_max_iter = 30, WI_tol = 5e-4, verbose = True)
    temp = np.cumsum(np.insert(temp.squeeze(), 0, sub_dta[0]))
    res_kWI_Kalman[:,sim] = temp
        
    
np.savetxt("./sim_data/ARI/ARI_WI_Kalman_miss2.csv", np.column_stack((res_ord_Kalman, res_kWI_Kalman)), delimiter = ",")    
