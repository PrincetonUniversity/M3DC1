#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Apr  13 2020

@author: Andreas Kleiner
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import itertools as it
import m3dc1.fpylib as fpyl


def get_h(psip,psiv,a1,a2,a3,a4p,a4v,a5p,a5v,a6,a7,lc1,lc2,Wc,psic):
    # Normal length
    h1p = 1./(1./(a4p*(1 - np.exp(-np.abs(psip/a1 - 1)**a2)) + a7) + 1/lc1*(1./(1 + ((psip - psic)/Wc)**2)))
    h1v = 1./(1./(a4v*(1 - np.exp(-np.abs(psiv/a1 - 1)**a3)) + a7) + 1/lc1*(1./(1 + ((psiv - psic)/Wc)**2)))
    
    # Tangential length
    h2p = 1./(1./(a5p*(1 - np.exp(-np.abs(psip/a1 - 1)**a2)) + a6) + 1/lc2*(1./(1 + ((psip - psic)/Wc)**2)))
    h2v = 1./(1./(a5v*(1 - np.exp(-np.abs(psiv/a1 - 1)**a3)) + a6) + 1/lc2*(1./(1 + ((psiv - psic)/Wc)**2)))
    return h1p, h1v, h2p, h2v



def mesh_size(filename='sizefieldParam',params=[],outfile='sizefieldParam.new'):
    """
    Visual and interactive tool to adjust the parameters in sizefieldParam.
    
    Arguments:

    **filename**
    Name of input sizefieldParam file. Should be None if instead a list of
    parameters is used as input.

    **params**
    List of 13 parameters (a1,a2,a3,a4p,a4v,a5p,a5v,a6,a7,lc1,lc2,Wc,psic)
    used as input.

    **outfile**
    Name of output file. Be careful to not overwrite existing files.
    """
    if filename is not None:
        try:
            # Read sizefieldParam
            with open(filename, 'r') as f:
                for line in f:
                    data_str = np.array(line.split())
                    data = data_str.astype(float)
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
    
    h1p, h1v, h2p, h2v = get_h(psip,psiv,a1,a2,a3,a4p,a4v,a5p,a5v,a6,a7,lc1,lc2,Wc,psic)
    
    
    fig, ax = plt.subplots(1, 1,figsize=(8,9))
    
    start = len(ax.lines)
    colors = it.cycle(['C0','C0','C3','C3',(0.5, 0.5, 0.5),(0.5, 0.5, 0.5),(0.5, 0.5, 0.5),(0.5, 0.5, 0.5), (0.4, 0.4, 0.4),(0.4, 0.4, 0.4),(0.4, 0.4, 0.4),(0.4, 0.4, 0.4), (0.3, 0.3, 0.3),(0.3, 0.3, 0.3),(0.3, 0.3, 0.3),(0.3, 0.3, 0.3), (0.2, 0.2, 0.2),(0.2, 0.2, 0.2),(0.2, 0.2, 0.2),(0.2, 0.2, 0.2), (0.1, 0.1, 0.1),(0.1, 0.1, 0.1),(0.1, 0.1, 0.1),(0.1, 0.1, 0.1)])
    colors = it.islice(colors, start, None)
    
    # Plot lines in color. These lines are replotted when sliders are moved
    line1, = ax.plot(psip,h1p*1e2,c=next(colors),lw=2,label='normal plasma',zorder=6)
    line2, = ax.plot(psiv,h1v*1e2,c=next(colors),lw=1,label='normal plasma',zorder=6)
    line3, = ax.plot(psip,h2p*1e2,c=next(colors),lw=2,label='tangential plasma',zorder=6)
    line4, = ax.plot(psiv,h2v*1e2,c=next(colors),lw=1,label='tangential vacuum',zorder=6)
    
    # Plot same lines in gray to indicate the original mesh size after sliders have been moved
    line01, = ax.plot(psip,h1p*1e2,c=next(colors),lw=2,zorder=1)
    line02, = ax.plot(psiv,h1v*1e2,c=next(colors),lw=1,zorder=1)
    line03, = ax.plot(psip,h2p*1e2,c=next(colors),lw=2,zorder=1)
    line04, = ax.plot(psiv,h2v*1e2,c=next(colors),lw=1,zorder=1)
    
    
    plt.grid(True)
    ax.set_xlabel(r'$\psi_N$')
    ax.set_ylabel(r'Mesh Size (cm)')
    plt.legend(loc=0)
    
    
    fig.subplots_adjust(top=0.98,bottom=0.5)
    
    #Parameters for slider positioning
    sl_x = 0.1
    sl_wid = 0.8
    sl_height = 0.02
    bot = 0.03
    delta_h = 0.04
    
    #Slider definition
    axa1 = fig.add_axes([sl_x, bot+9*delta_h, sl_wid, sl_height])
    a1_sl = Slider(ax=axa1, label='a1', valmin=0, valmax=5, valinit=a1)
    
    axa2 = fig.add_axes([sl_x, bot+8*delta_h, sl_wid, sl_height])
    a2_sl = Slider(ax=axa2, label='a2', valmin=0, valmax=5, valinit=a2)
    
    axa3 = fig.add_axes([sl_x, bot+7*delta_h, sl_wid, sl_height])
    a3_sl = Slider(ax=axa3, label='a3', valmin=0, valmax=5, valinit=a3)
    
    axa4p = fig.add_axes([sl_x, bot+6*delta_h, sl_wid, sl_height])
    a4p_sl = Slider(ax=axa4p, label='a4p', valmin=0, valmax=a4p*10, valinit=a4p)

    axa4v = fig.add_axes([sl_x, bot+5*delta_h, sl_wid, sl_height])
    a4v_sl = Slider(ax=axa4v, label="a4v", valmin=0, valmax=a4v*10, valinit=a4v)

    axa5p = fig.add_axes([sl_x, bot+4*delta_h, sl_wid, sl_height])
    a5p_sl = Slider(ax=axa5p, label='a5p', valmin=0, valmax=a5p*10, valinit=a5p)

    axa5v = fig.add_axes([sl_x, bot+3*delta_h, sl_wid, sl_height])
    a5v_sl = Slider(ax=axa5v, label="a5v", valmin=0, valmax=a5v*10, valinit=a5v)

    axa6 = fig.add_axes([sl_x, bot+2*delta_h, sl_wid, sl_height])
    a6_sl = Slider(ax=axa6, label='a6', valmin=0, valmax=a6*10, valinit=a6)

    axa7 = fig.add_axes([sl_x, bot+1*delta_h, sl_wid, sl_height])
    a7_sl = Slider(ax=axa7, label='a7', valmin=0, valmax=a7*10, valinit=a7)


    # The function to be called upon slider movement
    def update(val):
        h1p, h1v, h2p, h2v = get_h(psip,psiv,a1_sl.val,a2_sl.val,a3_sl.val,a4p_sl.val,a4v_sl.val,a5p_sl.val,a5v_sl.val,a6_sl.val,a7_sl.val,lc1,lc2,Wc,psic)
        line1.set_ydata(h1p*1e2)
        line2.set_ydata(h1v*1e2)
        line3.set_ydata(h2p*1e2)
        line4.set_ydata(h2v*1e2)
        fig.canvas.draw_idle()


    # Update plot upon slider movement
    a1_sl.on_changed(update)
    a2_sl.on_changed(update)
    a3_sl.on_changed(update)
    a4p_sl.on_changed(update)
    a4v_sl.on_changed(update)
    a5p_sl.on_changed(update)
    a5v_sl.on_changed(update)
    a6_sl.on_changed(update)
    a7_sl.on_changed(update)
    
    # Buttons to reset values and export sizefieldParam
    resetax = fig.add_axes([0.65, 0.01, 0.1, 0.03])
    reset_button = Button(resetax, 'Reset', hovercolor='0.975')

    exportax = fig.add_axes([0.8, 0.01, 0.1, 0.03])
    export_button = Button(exportax, 'Export', hovercolor='0.975')

    def reset(event):
        a1_sl.reset()
        a2_sl.reset()
        a3_sl.reset()
        a4p_sl.reset()
        a4v_sl.reset()
        a5p_sl.reset()
        a5v_sl.reset()
        a6_sl.reset()
        a7_sl.reset()
    
    reset_button.on_clicked(reset)

    def write_sizefieldparam(event):
        print(a1_sl.val,a2_sl.val,a3_sl.val,a4p_sl.val,a4v_sl.val,a5p_sl.val,a5v_sl.val,a6_sl.val,a7_sl.val,lc1,lc2,Wc,psic)
        os.system('cp '+filename + ' ' + filename+'.original')
        format_str = "{0:4.2f} {1:4.2f} {2:4.2f} {3:3.3f} {4:3.3f} {5:3.3f} {6:3.3f} {7:3.3f} {8:3.3f} {9:5.1f} {10:5.1f} {11:4.2f} {12:3.2f}"
        with open(outfile, 'w') as f:
            f.write(format_str.format(a1_sl.val,a2_sl.val,a3_sl.val,a4p_sl.val,a4v_sl.val,a5p_sl.val,a5v_sl.val,a6_sl.val,a7_sl.val,lc1,lc2,Wc,psic) + '\n')
            f.write('1.75 0 0.\n')
        print('Parameters written to file '+outfile)

    export_button.on_clicked(write_sizefieldparam)
    # Return buttons, because they would be garbage-collected otherwise
    return reset_button, export_button
