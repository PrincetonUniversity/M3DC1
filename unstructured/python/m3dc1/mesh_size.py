#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Apr  13 2020

@author: Andreas Kleiner
"""
import numpy as np
import matplotlib.pyplot as plt
import itertools as it
import m3dc1.fpylib as fpyl

def mesh_size(filename='sizefieldParam',params=[],fignum=None):
    if filename is not None:
        try:
            # Read sizefieldParam
            with open(filename, 'r') as f:
                for line in f:
                    data_str = np.array(line.split())
                    data = data_str.astype(np.float)
                    break
            print(data)

        except:
            fpyl.printerr('ERROR: Could not read '+ filename)
    else:
        if type(params)==str:
            data = np.asarray(params.split(),dtype=float)
        elif isinstance(params,(list,tuple,np.ndarray)):
            if len(params)==13:
                data = params
            else:
                raise ValueError("params should have a length of 13")
    
    # Define parameters
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
    
    #if grey:
    #    col_norm = 
    #    col_tang = 
    #else:
    #    col_norm = 'C0'
    #    col_tang = 'C3'

    plt.figure(fignum)
    ax = plt.gca()
    #print(ax.lines)
    #print(len(ax.lines))
    start = len(ax.lines)
    colors = it.cycle(['C0','C0','C3','C3',(0.5, 0.5, 0.5),(0.5, 0.5, 0.5),(0.5, 0.5, 0.5),(0.5, 0.5, 0.5), (0.4, 0.4, 0.4),(0.4, 0.4, 0.4),(0.4, 0.4, 0.4),(0.4, 0.4, 0.4), (0.3, 0.3, 0.3),(0.3, 0.3, 0.3),(0.3, 0.3, 0.3),(0.3, 0.3, 0.3), (0.2, 0.2, 0.2),(0.2, 0.2, 0.2),(0.2, 0.2, 0.2),(0.2, 0.2, 0.2), (0.1, 0.1, 0.1),(0.1, 0.1, 0.1),(0.1, 0.1, 0.1),(0.1, 0.1, 0.1)])
    colors = it.islice(colors, start, None)
    
    plt.plot(psip,h1p*1e2,c=next(colors),lw=2,label='normal plasma')
    plt.plot(psiv,h1v*1e2,c=next(colors),lw=1,label='normal plasma')
    plt.plot(psip,h2p*1e2,c=next(colors),lw=2,label='tangential plasma')
    plt.plot(psiv,h2v*1e2,c=next(colors),lw=1,label='tangential vacuum')
    plt.grid(True)
    plt.xlabel(r'$\psi_N$')
    plt.ylabel(r'Mesh Size (cm)')
    plt.legend(loc=0)
    return
