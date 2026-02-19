# -*- coding: utf-8 -*-
"""
Created on Fri Dec 26 14:29:48 2025

@author: whisk
"""
## imports
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import pandas as pd

## plot settings
fsize=22 ## font for xy label
nsize=20 ## number for tics
ffsize=26
lsize=19
MS=10
AP=0.7
mpl.rcParams['font.size']=fsize
mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['font.family']='serif'
mpl.rc('text', usetex=True)




## basic parameters
beta=1/0.592
convert= 25.7*beta

## skipping index
index = list(range(22))
skip = [2, 4, 6, 8]
index = [v for i, v in enumerate(index) if i not in skip]

## functions
def gradient_with_uncertainty(x, y, sigma_y):
    x = np.asarray(x); y = np.asarray(y); s = np.asarray(sigma_y)
    n = len(x)
    g  = np.empty(n)
    sg = np.empty(n)

    # interior points: 3-point central, non-uniform
    for i in range(1, n-1):
        h1 = x[i]   - x[i-1]
        h2 = x[i+1] - x[i]
        a_m1 = -h2/(h1*(h1+h2))
        a_0  =  (h2-h1)/(h1*h2)
        a_p1 =  h1/(h2*(h1+h2))
        g[i]  = a_m1*y[i-1] + a_0*y[i] + a_p1*y[i+1]
        sg[i] = np.sqrt((a_m1*s[i-1])**2 + (a_0*s[i])**2 + (a_p1*s[i+1])**2)

    # left edge: forward difference
    h = x[1] - x[0]
    g[0]  = (y[1]-y[0])/h
    sg[0] = np.sqrt(s[1]**2 + s[0]**2)/abs(h)

    # right edge: backward difference
    h = x[-1] - x[-2]
    g[-1]  = (y[-1]-y[-2])/h
    sg[-1] = np.sqrt(s[-1]**2 + s[-2]**2)/abs(h)

    return g, sg

def divisionE(A, eA, B, eB): ## return error of A/B
    val = A/B
    EA = eA/A
    EB = eB/B
    return val * np.sqrt( EA**2 + EB**2)

