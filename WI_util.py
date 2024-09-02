#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: shuochieh
"""

#%%

import ot
import numpy as np


#%%

def ts2stack(x, p = 2):
    """
    Return an numpy array consisting of lagged realizations of x
    
    Parameters
    ----------
    x : (n, 1) numpy array of time series
    p : dimension of the lagged time series; must be >= 2

    Returns
    -------
    A numpy array of size (n - p + 1, p)
    """
    
    n = len(x) - p + 1
    return np.lib.stride_tricks.as_strided(x, shape = (n, p), 
                                           strides = (x.itemsize, x.itemsize))

def slide_sum(A, d):
    """
    Parameters
    ----------
    A : A square 2-dim array
    d : dimension of the output array
    """
    
    res = np.zeros((d, d))
    k = A.shape[0]
    
    for i in range(d - k + 1):
        res[i:(i + k),i:(i + k)] += A
    return res

def solver_ordinary(x_obs, idx_obs, H, alpha, max_iter = 200, tol = 1e-4, verbose = False):
    """
    - Solve step (b) of the alternating minimization algorithm when there is no 
    side information
    - Euclidean distance is used as ground metric

    Parameters
    ----------
    x_obs : 
        (n, 1) array of observed time series (with initialized values)
    idx_obs : 
        a list of time indices corresponding to the observed values
    H : 
        (n, n) array of coefficients in the quadratic objective
    alpha : 
        learning rate of the proximal gradient descent algorithm (PDG)
    max_iter : 
        maximum number of iterations for PDG The default is 200.
    tol : 
        optimization toleration. The default is 1e-4.
    verbose : 
        whether to be verbose The default is False.

    Returns
    -------
    an array of size (n, 1) 
    """
    x_out = np.copy(x_obs)
    for _ in range(max_iter):
        temp = (x_out - alpha * np.dot(H, x_out)).reshape(-1, 1)
        temp[idx_obs, 0] = x_obs[idx_obs, 0]
        
        criterion = np.linalg.norm(x_out - temp) / np.sqrt(len(x_obs))
        x_out = temp
        if criterion < tol:
            return x_out
        
    if verbose:
        print(f"    solver: proximal gradient descent may have not converged; criterion: {criterion:.4f}")
        
    return x_out

def WI_core_ordinary(x, n1, p, idx_obs, solver, WI_max_iter = 1000, 
                     WI_tol = 1e-4, verbose = False, **kwargs):
    """
    - Implements Wasserstein imputation when there is no side information; 
    The target time series x may be multivariate. Certain speed-ups are applied 
    utilizing the purely missing structure
    - Step (b) is solved with the provided solver (w_solver)
 
    Parameters
    ----------
    x : 
        (n, d) array of time series, with missing entries already initialized
    n1 : 
        user-defined cutoff point in time
    p : 
        dimension of the marginal distribution to match
    idx_obs : 
        a list of lists of indices corresponding to the observed entries (each
        corresponds to one variable)
    solver : 
        generic solver for step (b) of the alternating minization algorithm; 
        must have the first few arguments as (x_obs, idx_obs, H). 
    WI_max_iter : 
        maximum number of iterations. The default is 1000.
    WI_tol : 
        optimization tolerance. The default is 1e-4.
    verbose : 
        whether to be verbose. The default is False.
    **kwargs : 
        Additional arguments for solver

    Returns
    -------
    A (n, d) array imputed by Wasserstein imputation
    """
    n = x.shape[0]
    d = x.shape[1]
    h1, h2 = np.ones((n1 - p + 2,)) / (n1 - p + 2), np.ones((n - n1 - 1,)) / (n - n1 - 1)
    h1 = h1 / np.sum(h1)
    h2 = h2 / np.sum(h2)
    loss = np.zeros(WI_max_iter)
    
    for i in range(WI_max_iter):
        for j in range(d):
            if j == 0:
                x_lags = ts2stack(x[:,j], p)
            else: 
                x_lags = np.column_stack((x_lags, ts2stack(x[:,j], p)))
        
        x_pre = x_lags[:(n1 - p + 2),:]
        x_post = x_lags[(n1 - p + 2):,:]
        
        
#%%



















