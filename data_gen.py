#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: shuochieh
"""

#%%
import numpy as np

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
        