#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: shuochieh
"""

#%%
import numpy as np
import random
import math
import matplotlib.pyplot as plt

#%%
def arma_sim(n, phi, psi):
    e = np.random.normal(size = n + 100, scale = 1)
    res = np.zeros((n + 100, 1))
    for i in range(1, n + 100):
        res[i] = phi * res[i - 1] + e[i] - psi * e[i - 1]
        
    return res[100:]

def TAR_sim(n, phi1, phi2, thres, s1 = 1, s2 = 1):
    res = np.zeros((n + 100, 1))
    for i in range(1, n + 100):
        if res[i - 1,0] <= thres:
            res[i] = phi1 * res[i - 1] + np.random.normal(size = 1, scale = s1)
        else:
            res[i] = phi2 * res[i - 1] + np.random.normal(size = 1, scale = s2)

    return res[100:]

def sigmoid(x):
    return 1 / (1 + np.exp(-x)) - 0.5

def NLVAR_sim(n, B, s1 = 1, s2 = 1):
    res = np.zeros((n + 100, 2))
    for i in range(1, n + 100):
        res[i, 0] = B[0, 0] * res[i - 1, 0] + B[0, 1] * sigmoid(3 * res[i - 1, 1]) + np.random.normal(size = 1, scale = s1)
        res[i, 1] = B[1, 0] * res[i - 1, 0] + B[1, 1] * res[i - 1, 1] + np.random.normal(size = 1, scale = s2)
        
    return res[100:]

def garch_sim(n, omega, alpha, beta):
    res = np.zeros((n + 100, 1))
    z = np.random.normal(size = n + 100)
    sigma2 = 0
    for i in range(2, n + 100):
        sigma2 = omega + alpha * np.power(res[i - 1,0], 2) + beta * sigma2
        res[i] = np.power(sigma2, 0.5) * z[i]
        
    return res[100:]

def alr_composition(n, k, A):
    """
    Create a compositional time series with k categories
    
    Inputs:
        n: sample size
        k: num of categories
        A: (k-1) by (k-1) matrix for the transformed VAR model
        
    Outputs:
        res: resulting compositional time series
        x: latent VAR realizations
    """
    res = np.zeros((n + 100, k))
    x = np.zeros((n + 100, k - 1))
    for i in range(4, n + 100):
        x[i,:] = np.array([0.1] + [0.1] * (k - 2)) + np.dot(A, x[i - 1,:]) + 0.15 * x[i - 4,:] + np.random.normal(size = k - 1, scale = 0.1)
        res[i, k - 1] = 1 / (1 + np.sum(np.exp(x[i,:])))
        for j in range(k - 1):
            res[i, j] = res[i,k - 1] * np.exp(x[i, j])
            
    return res[100:], x[100:]
        
def ARI_sim(n, alpha):
    res = np.zeros((n + 10, 1))
    dx = np.zeros((n + 10, 1))
    for i in range(2, n + 10):
        dx[i, 0] = alpha * dx[i - 1, 0] + np.random.normal(scale = 0.5)
        res[i, 0] = res[i - 1, 0] + dx[i, 0] + np.random.normal(scale = 1)
        
    return res[10:], dx[10:]

def cyc_sim(n):
    res = np.zeros((n, 1))
    for i in range(n):
        res[i, 0] = 10 * np.cos(i * (0.23) * np.pi) + 6 * np.cos(i * 0.17 * np.pi) + np.random.normal(scale = 0.3)
    return res

#%% 

# Sanity checks
n = 300

x = arma_sim(n, 0.8, 0)
plt.plot(x[:(n - 1)], x[1:], 'co')
plt.title("AR(1): original series")
plt.show()


x = arma_sim(n, 0.8, 0.7)
plt.plot(x[:(n - 1)], x[1:], 'co')
plt.title("ARMA(1,1): original series")
plt.show()
    
x = TAR_sim(n, -2, 0.7, 1, 1, 0.5)
plt.plot(x[:(n - 1)], x[1:], 'co')
plt.title("TAR(1): original series")
plt.show()

B = np.array([[0.3, 8], [0, 0.4]])
x = NLVAR_sim(n, B, s1 = 0.25, s2 = 3)
plt.plot(x[:(n - 1),1], x[1:,0], 'co')
plt.title("NL-VAR: original series")
plt.show()

x = garch_sim(n, 0.5, 0.8, 0.1)
plt.plot(x[:(n - 1)], x[1:], 'co')
plt.title("GARCH(1,1): original series")
plt.show()

A = np.array([[0.75, -0.1], [0, 0.5]])
x, a = alr_composition(n, 3, A)
plt.plot(np.linspace(0, n - 1, n), x[:,0], 'r-')
plt.plot(np.linspace(0, n - 1, n), x[:,1], 'g-')
plt.plot(np.linspace(0, n - 1, n), x[:,2], 'b-')
plt.title("ALR: original series")
plt.show()

x, dx = ARI_sim(n, -0.7)
plt.plot(np.linspace(0, n - 1, n), x[:,0], 'r-')
plt.title("ARI(1,1): original series")
plt.show()
plt.plot(x[:(n - 1)], x[1:], 'co')
plt.title("ARI(1,1): original series")
plt.show()
plt.plot(dx[:(n - 1)], dx[1:], 'co')
plt.title("ARI(1,1): original series")
plt.show()

x = cyc_sim(n)
plt.plot(np.linspace(0, n - 1, n), x, 'c-', linewidth = 0.8)
plt.title("Cyclical time series: original series")
plt.show()
plt.plot(x[:(n - 1),0], x[1:,0], 'co')
plt.show()
plt.plot(x[:(n - 4),0], x[4:,0], 'co')
plt.show()
plt.plot(x[:(n - 6),0], x[6:,0], 'co')
plt.show()
plt.plot(x[:(n - 12),0], x[12:,0], 'co')
plt.show()

#%%

### Data generation

n = 3000 # Sample size
n_exper = 1000 # Number of Monte Carlo experiment

# AR(1) Model
res_original = np.zeros((n, n_exper))
res_miss1 = np.zeros((n, n_exper))
res_miss2 = np.zeros((n, n_exper))

for b in range(n_exper):
    temp = arma_sim(n, 0.8, 0).squeeze()
    res_original[:,b] = temp
    m1 = random.sample(range(1, n - 1), k = math.floor(n * 0.3))
    m2 = list(range(math.floor(0.2 * n), math.ceil(0.3 * n))) + list(range(math.floor(0.7 * n), math.ceil(0.9 * n)))
    miss1 = np.copy(temp)
    miss1[m1] = np.nan
    res_miss1[:,b] = miss1
    
    miss2 = np.copy(temp)
    miss2[m2] = np.nan
    res_miss2[:,b] = miss2
    if (b+1) % 100 == 0 and b > 98:
        print(f"AR data generation: {b+1} / {n_exper}")


np.savetxt("./sim_data/AR_original.csv", res_original, delimiter = ",")
np.savetxt("./sim_data/AR_miss1.csv", res_miss1, delimiter = ",")
np.savetxt("./sim_data/AR_miss2.csv", res_miss2, delimiter = ",")

#%%

# ARMA(1,1) Model
res_original = np.zeros((n, n_exper))
res_miss1 = np.zeros((n, n_exper))
res_miss2 = np.zeros((n, n_exper))

for b in range(n_exper):
    temp = arma_sim(n, 0.8, 0.6).squeeze()
    res_original[:,b] = temp
    m1 = random.sample(range(1, n - 1), k = math.floor(n * 0.3))
    m2 = list(range(math.floor(0.2 * n), math.ceil(0.3 * n))) + list(range(math.floor(0.7 * n), math.ceil(0.9 * n)))
    miss1 = np.copy(temp)
    miss1[m1] = np.nan
    res_miss1[:,b] = miss1
    
    miss2 = np.copy(temp)
    miss2[m2] = np.nan
    res_miss2[:,b] = miss2
    if (b+1) % 100 == 0 and b > 98:
        print(f"ARMA data generation: {b+1} / {n_exper}")


np.savetxt("./sim_data/ARMA_original.csv", res_original, delimiter = ",")
np.savetxt("./sim_data/ARMA_miss1.csv", res_miss1, delimiter = ",")
np.savetxt("./sim_data/ARMA_miss2.csv", res_miss2, delimiter = ",")

#%%

# TAR(1) Model
res_original = np.zeros((n, n_exper))
res_miss1 = np.zeros((n, n_exper))
res_miss2 = np.zeros((n, n_exper))

for b in range(n_exper):
    temp = TAR_sim(n, -2, 0.7, 1, 1, 0.5).squeeze()
    res_original[:,b] = temp
    m1 = random.sample(range(1, n - 1), k = math.floor(n * 0.3))
    m2 = list(range(math.floor(0.2 * n), math.ceil(0.3 * n))) + list(range(math.floor(0.7 * n), math.ceil(0.9 * n)))
    miss1 = np.copy(temp)
    miss1[m1] = np.nan
    res_miss1[:,b] = miss1
    
    miss2 = np.copy(temp)
    miss2[m2] = np.nan
    res_miss2[:,b] = miss2
    if (b+1) % 100 == 0 and b > 98:
        print(f"TAR data generation: {b+1} / {n_exper}")


np.savetxt("./sim_data/TAR_original.csv", res_original, delimiter = ",")
np.savetxt("./sim_data/TAR_miss1.csv", res_miss1, delimiter = ",")
np.savetxt("./sim_data/TAR_miss2.csv", res_miss2, delimiter = ",")


#%%

# Nonlinear VAR Model
res_original = np.zeros((n, n_exper * 2))
res_miss1 = np.zeros((n, n_exper * 2))
res_miss2 = np.zeros((n, n_exper * 2))

B = np.array([[0.3, 8], [0, 0.4]])

for b in range(n_exper):
    temp = NLVAR_sim(n, B, s1 = 0.25, s2 = 3)
    res_original[:,2 * b] = temp[:,0]
    res_original[:,2 * b + 1] = temp[:,1]
    m1 = random.sample(range(1, n - 1), k = math.floor(n * 0.3))
    m2 = list(range(math.floor(0.2 * n), math.ceil(0.3 * n))) + list(range(math.floor(0.7 * n), math.ceil(0.9 * n)))
    miss1 = np.copy(temp)
    miss1[m1] = np.nan
    res_miss1[:,2 * b] = miss1[:,0]
    res_miss1[:,2 * b + 1] = miss1[:,1]
    
    miss2 = np.copy(temp)
    miss2[m2] = np.nan
    res_miss2[:,2 * b] = miss2[:,0]
    res_miss2[:,2 * b + 1] = miss2[:,1]
    if (b+1) % 100 == 0 and b > 98:
        print(f"Nonlinear VAR data generation: {b+1} / {n_exper}")

np.savetxt("./sim_data/NLVAR_original.csv", res_original, delimiter = ",")
np.savetxt("./sim_data/NLVAR_miss1.csv", res_miss1, delimiter = ",")
np.savetxt("./sim_data/NLVAR_miss2.csv", res_miss2, delimiter = ",")


#%%

# GARCH(1, 1) Model
res_original = np.zeros((n, n_exper))
res_miss1 = np.zeros((n, n_exper))
res_miss2 = np.zeros((n, n_exper))

for b in range(n_exper):
    temp = garch_sim(n, 0.5, 0.8, 0.1).squeeze()
    res_original[:,b] = temp
    m1 = random.sample(range(1, n - 1), k = math.floor(n * 0.3))
    m2 = list(range(math.floor(0.2 * n), math.ceil(0.3 * n))) + list(range(math.floor(0.7 * n), math.ceil(0.9 * n)))
    miss1 = np.copy(temp)
    miss1[m1] = np.nan
    res_miss1[:,b] = miss1
    
    miss2 = np.copy(temp)
    miss2[m2] = np.nan
    res_miss2[:,b] = miss2
    if (b+1) % 100 == 0 and b > 98:
        print(f"GARCH data generation: {b+1} / {n_exper}")


np.savetxt("./sim_data/GARCH_original.csv", res_original, delimiter = ",")
np.savetxt("./sim_data/GARCH_miss1.csv", res_miss1, delimiter = ",")
np.savetxt("./sim_data/GARCH_miss2.csv", res_miss2, delimiter = ",")

#%%

# ALR-composition Model
res_original = np.zeros((n, n_exper * 3))
res_miss1 = np.zeros((n, n_exper * 3))
res_miss2 = np.zeros((n, n_exper * 3))

A = np.array([[0.75, -0.1], [0, 0.5]])

for b in range(n_exper):
    temp, _ = alr_composition(n, 3, A)
    res_original[:,3 * b] = temp[:,0]
    res_original[:,3 * b + 1] = temp[:,1]
    res_original[:,3 * b + 2] = temp[:,2]
    m1 = random.sample(range(1, n - 1), k = math.floor(n * 0.3))
    m2 = list(range(math.floor(0.2 * n), math.ceil(0.3 * n))) + list(range(math.floor(0.7 * n), math.ceil(0.9 * n)))
    miss1 = np.copy(temp)
    miss1[m1] = np.nan
    res_miss1[:,3 * b] = miss1[:,0]
    res_miss1[:,3 * b + 1] = miss1[:,1]
    res_miss1[:,3 * b + 2] = miss1[:,2]
    
    miss2 = np.copy(temp)
    miss2[m2] = np.nan
    res_miss2[:,3 * b] = miss2[:,0]
    res_miss2[:,3 * b + 1] = miss2[:,1]
    res_miss2[:,3 * b + 2] = miss2[:,2]
    if (b+1) % 100 == 0 and b > 98:
        print(f"ALR data generation: {b+1} / {n_exper}")

np.savetxt("./sim_data/ALR_original.csv", res_original, delimiter = ",")
np.savetxt("./sim_data/ALR_miss1.csv", res_miss1, delimiter = ",")
np.savetxt("./sim_data/ALR_miss2.csv", res_miss2, delimiter = ",")

#%%

# ARI(1,1) Model
res_original = np.zeros((n, n_exper))
res_miss1 = np.zeros((n, n_exper))
res_miss2 = np.zeros((n, n_exper))


for b in range(n_exper):
    temp, _ = ARI_sim(n, -0.7)
    res_original[:,b] = temp.squeeze()
    m1 = random.sample(range(1, n - 1), k = math.floor(n * 0.3))
    m2 = list(range(math.floor(0.2 * n), math.ceil(0.3 * n))) + list(range(math.floor(0.7 * n), math.ceil(0.9 * n)))
    miss1 = np.copy(temp)
    miss1[m1] = np.nan
    res_miss1[:,b] = miss1.squeeze()
    
    miss2 = np.copy(temp)
    miss2[m2] = np.nan
    res_miss2[:,b] = miss2.squeeze()
    if (b+1) % 100 == 0 and b > 98:
        print(f"ARI data generation: {b+1} / {n_exper}")

np.savetxt("./sim_data/ARI_original.csv", res_original, delimiter = ",")
np.savetxt("./sim_data/ARI_miss1.csv", res_miss1, delimiter = ",")
np.savetxt("./sim_data/ARI_miss2.csv", res_miss2, delimiter = ",")


#%%

# Cyclic time series Model
res_original = np.zeros((n, n_exper))
res_miss1 = np.zeros((n, n_exper))
res_miss2 = np.zeros((n, n_exper))


for b in range(n_exper):
    temp = cyc_sim(n).squeeze()
    res_original[:,b] = temp
    m1 = random.sample(range(1, n - 1), k = math.floor(n * 0.3))
    m2 = list(range(math.floor(0.2 * n), math.ceil(0.3 * n))) + list(range(math.floor(0.7 * n), math.ceil(0.9 * n)))
    miss1 = np.copy(temp)
    miss1[m1] = np.nan
    res_miss1[:,b] = miss1
    
    miss2 = np.copy(temp)
    miss2[m2] = np.nan
    res_miss2[:,b] = miss2
    if (b+1) % 100 == 0 and b > 98:
        print(f"Cyclic data generation: {b+1} / {n_exper}")

np.savetxt("./sim_data/Cyc_original.csv", res_original, delimiter = ",")
np.savetxt("./sim_data/Cyc_miss1.csv", res_miss1, delimiter = ",")
np.savetxt("./sim_data/Cyc_miss2.csv", res_miss2, delimiter = ",")








