#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Feb  14 2022

@author: Andreas Kleiner
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as mplax
import m3dc1.fpylib as fpyl



def plot_mag_probes(filename='',cycle_col=False,ax=None,fignum=None,pub=False,quiet=True):
    """
    Plots a field in the R,Z plane.
    
    Arguments:

    **filename**
    Name of the Slurm log file with magnetic probe information.

    **cycle_col**
    If True, use matplotlib color cycle in coil plot. If False, all coils are plotted in same color.

    **ax**
    axis object to plot coils into.

    **fignum**
    Number of figure for the plot.

    **pub**
    If True, format figure for publication (larger labels and thicker lines)

    **quiet**
    If True, suppress output to terminal.
    """

    slurmfile = fpyl.get_input_parameter_file()

    mag_probe_phi = np.asarray(fpyl.get_parameter_from_ascii('mag_probe_phi',slurmfile,quiet=quiet),dtype=np.float32)
    mag_probe_x = np.asarray(fpyl.get_parameter_from_ascii('mag_probe_x',slurmfile,quiet=quiet),dtype=np.float32)
    mag_probe_z = np.asarray(fpyl.get_parameter_from_ascii('mag_probe_z',slurmfile,quiet=quiet),dtype=np.float32)

    ind_phizero = np.where(mag_probe_phi==0)[0]
    R_mag_probe = np.take(mag_probe_x,ind_phizero)
    Z_mag_probe = np.take(mag_probe_z,ind_phizero)

    non_zero_elements = np.logical_or(R_mag_probe!=0, Z_mag_probe!=0) #Create mask of elements that have R or Z non-zero
    R_mag_probe = R_mag_probe[non_zero_elements] #Remove R,Z = 0,0 points from array
    Z_mag_probe = Z_mag_probe[non_zero_elements] #Remove R,Z = 0,0 points from array
    #print(mag_probe_phi)
    #print(ind_phizero)
    #print(R_mag_probe)
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
        #linew = 1
    else:
        axlblfs = 12
        ticklblfs = 12
        #linew = 1
    
    newfig=False
    if not isinstance(ax, (np.ndarray, mplax._axes.Axes)):
        fig = plt.figure(num=fignum)
        ax = plt.gca()
        newfig=True
    
    
    #Plot magnetic probes
    for i in range(len(R_mag_probe)):
        axarray = np.atleast_1d(ax)
        for ax in axarray:
            ax.plot(R_mag_probe[i],Z_mag_probe[i],c=line_col,lw=0,marker='.',zorder=15)
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
