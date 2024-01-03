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
from labellines import labelLines

def plot_coils(filename='coil.dat',angleUnits='deg',cycle_col=False,ax=None,print_coil_pos=False,show_labels=False,pick=False,fignum=None,pub=False):
    """
    Plots a field in the R,Z plane.
    
    Arguments:

    **filename**
    Name of the file with coil information, default is 'coil.dat'

    **angleUnits**
    'deg'/'rad'
    Units angles are given in.

    **cycle_col**
    If True, use matplotlib color cycle in coil plot. If False, all coils are plotted in same color.

    **ax**
    axis object to plot coils into.

    **print_coil_pos**
    If True, print (R,Z) values that define coil contour.

    **show_labels**
    If True, show name of each structure in plot. Only works when this information is present
    in coil file.

    **pick**
    If True, make points that define boundaries clickable, so the coordinates can directly be shown.

    **fignum**
    Number of figure for the plot.

    **pub**
    If True, format figure for publication (larger labels and thicker lines)
    """
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
    
    Rpts = []
    Zpts = []
    
    #Plot coils using same color for each coil winding of the same coil
    for line in data:
        coil = line.split()
        coil_format = len(coil)
        detailed = False #This variable is set to True, if the coil names are available in the input file.
        if coil_format>11:
            if coil[0].replace('.0','',1).isdigit() or coil[0].replace('.1','',1).isdigit():
                #general MAST-U format
                rc=float(coil[1])
                zc=float(coil[2])
                dr=float(coil[3])
                dz=float(coil[4])
                a1=float(coil[5])
                a2=float(coil[6])
            else:
                detailed = True
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
        
        Rpts.append(rcoord)
        Zpts.append(zcoord)
        
        if print_coil_pos:
            if detailed:
                print(coil[0])
            else:
                print('-----------------')
            print("\n".join("{} {}".format(R,Z) for R,Z in zip(rcoord, zcoord)))
        
        #if abs(rc/0.23560780-1)<0.1:
        #    #print(rc,zc,a1,np.tan(a1),a2,np.tan(a2))
        axarray = np.atleast_1d(ax)
        if show_labels and detailed:
            leglbl = coil[0]
        else:
            leglbl = None
        for ax in axarray:
            ax.plot(rcoord,zcoord,c=line_col,lw=linew,zorder=15,label=leglbl)
        if cycle_col:
            line_col = next(cols)
        
        if show_labels and detailed:
            plt.text(rc, zc,coil[0], ha='center', va='center', fontsize=10)
            #labelLines(ax.get_lines(), zorder=2.5)
    
    
    if pick:
        def line_picker(line, mouseevent):
            """
            find the points within a certain distance from the mouseclick in
            data coords and attach some extra attributes, pickx and picky
            which are the data points that were picked
            """
            if mouseevent.xdata is None:
                return False, dict()
            xdata = line.get_xdata()
            ydata = line.get_ydata()
            maxd = 0.001
            d = np.sqrt((xdata - mouseevent.xdata)**2. + (ydata - mouseevent.ydata)**2.)

            ind = np.nonzero(np.less_equal(d, maxd))
            if len(ind):
                pickx = np.take(xdata, ind)
                picky = np.take(ydata, ind)
                props = dict(ind=ind, pickx=pickx, picky=picky)
                return True, props
            else:
                return False, dict()

        def onpick2(event):
            print(event.pickx[0][0], event.picky[0][0])
        
        #Need to print markers separately outside of above loop, because otherwise it will loop through click events
        for ax in axarray:
            line, = ax.plot(np.asarray(Rpts).flatten(), np.asarray(Zpts).flatten(), 'o', picker=line_picker)
        fig.canvas.mpl_connect('pick_event', onpick2)
    
    
    if newfig:
        plt.grid(True)
        ax.set_aspect('equal',adjustable='box')
        plt.xlabel(r'$R$',fontsize=axlblfs)
        plt.ylabel(r'$Z$',fontsize=axlblfs)
        ax.tick_params(labelsize=ticklblfs)
        plt.tight_layout()
    
    return
