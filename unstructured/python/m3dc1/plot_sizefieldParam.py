#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Apr  13 2020

@author: Andreas Kleiner
"""
import numpy as np
import matplotlib.pyplot as plt


def plot_sizefieldParam(filename='sizefieldParam',fignum=None):
    # Read sizefieldParam
    with open(filename, 'r') as f:
        for line in f:
            data_str = np.array(line.split())
            data = data_str.astype(np.float)
            break
    print(data)
    a1   = data[0]
    a2   = data[1]
    a3   = data[2]
    a4p  = data[3]
    a4v  = data[4]
    a5p  = data[5]
    a5v  = data[6]
    a6   = data[7]
    a7   = data[8]
    lc1  = data[9]
    lc2  = data[10]
    Wc   = data[11]
    psic = data[12]
    
    psip = np.linspace(0.0,a1,101)
    psiv = np.linspace(a1,3.0,101)
    
    # Normal length
    h1p = 1./(1./(a4p*(1 - np.exp(-np.abs(psip/a1 - 1)**a2)) + a7) + 1/lc1*(1./(1 + ((psip - psic)/Wc)**2)))
    h1v = 1./(1./(a4v*(1 - np.exp(-np.abs(psiv/a1 - 1)**a3)) + a7) + 1/lc1*(1./(1 + ((psiv - psic)/Wc)**2)))
    
    # Tangential length
    h2p = 1./(1./(a5p*(1 - np.exp(-np.abs(psip/a1 - 1)**a2)) + a6) + 1/lc2*(1./(1 + ((psip - psic)/Wc)**2)))
    h2v = 1./(1./(a5v*(1 - np.exp(-np.abs(psiv/a1 - 1)**a3)) + a6) + 1/lc2*(1./(1 + ((psiv - psic)/Wc)**2)))
    
    plt.figure(fignum)
    plt.plot(psip,h1p*1e2,c='C0',label='normal plasma')
    plt.plot(psiv,h1v*1e2,c='C0',ls='--',label='normal plasma')
    plt.plot(psip,h2p*1e2,c='C3',label='tangential plasma')
    plt.plot(psiv,h2v*1e2,c='C3',ls='--',label='tangential vacuum')
    plt.grid(True)
    plt.xlabel(r'$\psi_N$')
    plt.ylabel(r'Mesh Size (cm)')
    plt.legend(loc=0)
    return
