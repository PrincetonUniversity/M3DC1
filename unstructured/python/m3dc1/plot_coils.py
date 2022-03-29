#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Feb  14 2022

@author: Andreas Kleiner
"""
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.axes as mplax

def plot_coils(filename='coil.dat',angleUnits='deg',cycle_col=False,ax=None,fignum=None,pub=False):
    f = open(filename, 'r')
    data = f.readlines()
    
    # Color cycler:
    def col_cycler(cols):
        count = 0
        while True:
            yield cols[count]
            count = (count + 1)%len(cols)
    cols = col_cycler(['C1','C2','C3','C4','C5','C6','C7','C8','C9','k','C0'])
    line_col = next(cols)
    
    # Set font sizes and plot style parameters
    if pub:
        axlblfs = 20
        ticklblfs = 18
        linew = 1
    else:
        axlblfs = 12
        ticklblfs = 12
        linew = 1
    
    newfig=False
    if not isinstance(ax, (np.ndarray, mplax._axes.Axes)):
        fig = plt.figure(num=fignum)
        ax = plt.gca()
        newfig=True
    
    #Plot coils using same color for each coil winding of the same coil
    for line in data:
        coil = line.split()
        coil_format = len(coil)
        if coil_format>11:
            rc=float(coil[3])
            zc=float(coil[4])
            dr=float(coil[5])
            dz=float(coil[6])
            a1=float(coil[7])
            a2=float(coil[8])
        else:
            rc=float(coil[0])
            zc=float(coil[1])
            dr=float(coil[2])
            dz=float(coil[3])
            a1=float(coil[4])
            a2=float(coil[5])
        
        #Convert angle from degrees to radian
        if angleUnits == 'deg':
            a1 = a1/180*math.pi
            a2 = a2/180*math.pi
        if a2!=0.0:
            a2 = a2 + math.pi/2
        rcoord = [rc-(dr-dz*np.tan(a2))/2, rc+(dr+dz*np.tan(a2))/2, rc+(dr-dz*np.tan(a2))/2, rc-(dr+dz*np.tan(a2))/2, rc-(dr-dz*np.tan(a2))/2]
        zcoord = [zc-(dz+dr*np.tan(a1))/2, zc-(dz-dr*np.tan(a1))/2, zc+(dz+dr*np.tan(a1))/2, zc+(dz-dr*np.tan(a1))/2, zc-(dz+dr*np.tan(a1))/2]
        #if abs(rc/0.23560780-1)<0.1:
        #    #print(rc,zc,a1,np.tan(a1),a2,np.tan(a2))
        axarray = np.atleast_1d(ax)
        for ax in axarray:
            ax.plot(rcoord,zcoord,c=line_col,lw=linew,zorder=15)
        if cycle_col:
            line_col = next(cols)
        
        if newfig:
            plt.grid(True)
            ax.set_aspect('equal',adjustable='box')
            plt.xlabel(r'$R$',fontsize=axlblfs)
            plt.ylabel(r'$Z$',fontsize=axlblfs)
            ax.tick_params(labelsize=ticklblfs)
            plt.tight_layout()
    return
