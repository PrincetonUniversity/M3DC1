#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 14:24:12 2019

@author: Andreas Kleiner
"""

import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from matplotlib import rc
from labellines import labelLines

import fpy
import m3dc1.fpylib as fpyl
from m3dc1.get_timetrace import get_timetrace
from m3dc1.unit_conv import unit_conv
from m3dc1.read_h5 import readParameter
from m3dc1.get_time_of_slice import get_time_of_slice
from m3dc1.time_trace_fast import create_plot_time_trace_fast


rc('text', usetex=True)
plt.rcParams.update({'figure.max_open_warning': 40})




def plot_time_trace_fast(trace,units='mks',millisec=False,sim=None,filename='C1.h5',diff=False,
                         growth=False,renorm=False,yscale='linear',unitlabel=None,fac=1,
                         show_legend=False,leglbl=None,in_plot_txt=None,time_marks=[],ts_marks=[],
                         ts_marks_all=False,rescale=False,save=False,savedir=None,pub=False,
                         fignum=None,figsize=None,drop_time_steps=None,skip_n0=False,export=False,txtname=None):
    """
    Plots a scalar quantity vs. time. All available
    time traces can be found in the M3D-C1 documentation.
    
    Arguments:

    **trace**
    String containing the trace to be plotted

    **units**
    The units in which the trace will be plotted

    **millisec**
    True/False. If True and units='mks' plot will be in terms of milliseconds, instead of seconds.

    **sim**
    fpy simulation object.

    **filename**
    Contains the path to the C1.h5 file.

    **diff**
    True / False. Plot derivative of scalar.

    **growth**
    Determines wether to calculate the derivative.
    True/False

    **renorm**
    Removes spikes that are caused by renormalizations
    in linear stability calculations. Interpolates at
    the locations of the spike. Should only be used if
    growth=True.

    **yscale**
    Scale of y axis, e.g. linear, log

    **unitlabel**
    Specify custom unitlabel. If not specified, default label will be used.

    **fac**
    Factor to apply to field when using mks units. fac=1.0E-3 converts to kilo..., fac=1.0E-6 to Mega...

    **show_legend**
    If True, show plot legend.

    **leglbl**
    Legend label for plot curve.

    **in_plot_txt**
    Overlay text to be shown inside of plot window.

    **time_marks**
    List of times that will be marked in terms of a vertical line in plot.

    **ts_marks**
    List of time slice numbers. A vertical line will be added for the time corresponding
    to each time slice listed.

    **ts_marks_all**
    True / False. If True, indicate times of all time slices.

    **rescale**
    Rescale y-axis such that noise in the beginning of the simulation is not considered for axis limits,
    i.e. plot is zoomed in to show relevant values.

    **save**
    If True, save plot to file

    **savedir**
    Directory where plot will be saved

    **pub**
    If True, format figure for publication (larger labels and thicker lines)

    **fignum**
    Figure number.

    **figsize**
    If True, format figure for publication (larger labels and thicker lines)

    **drop_time_steps**
    Number of time steps at the end of time trace to drop from plot.

    **skip_n0**
    When plotting energy spectrum, do not plot energy for n=0 mode.

    """
    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename)
    time,y_axis,label, unitlabel = get_timetrace(trace,sim=sim,units=units,growth=growth,diff=diff,renorm=renorm,unitlabel=unitlabel,fac=fac)
    if drop_time_steps is not None:
        time = time[:-drop_time_steps]
        y_axis = y_axis[:-drop_time_steps]
    if unitlabel is None or unitlabel == '':
        ylbl = label
    else:
        ylbl = label+' ('+unitlabel+')'
    
    filename = sim.filename
    create_plot_time_trace_fast(time,y_axis,trace,units=units,millisec=millisec,sim=sim,filename=filename,growth=growth,diff=diff,show_legend=show_legend,
                                leglbl=leglbl,yscale=yscale,rescale=rescale,save=save,savedir=savedir,pub=pub,in_plot_txt=in_plot_txt,
                                time_marks=time_marks,ts_marks=ts_marks,ts_marks_all=ts_marks_all,ylbl=ylbl,fignum=fignum,figsize=figsize,skip_n0=skip_n0,
                                export=export,txtname=txtname)
    return




def create_plot_time_trace_fast(time,scalar,trace,units='mks',millisec=False,sim=None,filename='C1.h5',
                                growth=False,diff=False,yscale='linear',rescale=False,save=False,in_plot_txt=None,
                                time_marks=[],ts_marks=[],ts_marks_all=False,savedir=None,title=False,pub=False,
                                ylbl=None,show_legend=False,leglbl=None,fignum=None,figsize=None,skip_n0=False,
                                export=False,txtname=None):

    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename)
        
    # If one array has new data but the other one doesn't 
    # plot only previous data
    if scalar.shape != time.shape:
        ymax   = np.amax(scalar.shape)
        tmax   = np.amax(time.shape)
        maxidx = np.amin([ymax,tmax])
        time   = time[0:maxidx]
        scalar = scalar[0:maxidx]
    
    ntor = readParameter('ntor',sim=sim)
    
    plt.figure(num=fignum,figsize=figsize)
    
    # Set font sizes and plot style parameters
    if pub:
        axlblfs = 20
        titlefs = 18
        ticklblfs = 16
        linew = 2
        legfs = 14
        inplttxtfs = 20
    else:
        axlblfs = None
        titlefs = None
        ticklblfs = None
        linew = 1
        legfs = None
        inplttxtfs = 16
    
    if units=='mks':
        #time = unit_conv(time,filename=filename,time=1)
        if millisec:
            time = time*1000
            plt.xlabel(r'time $(ms)$',fontsize=axlblfs)
        else:
            plt.xlabel(r'time $(s)$',fontsize=axlblfs)
    else:
        plt.xlabel(r'time $(\tau_A)$',fontsize=axlblfs)
    
    if diff:
        plt.ylabel('d '+ylbl+' / d t',fontsize=axlblfs)
    else:
        if trace == 'ke':
            if units=='mks':
                #scalar = unit_conv(scalar, arr_dim='M3DC1', filename=filename, energy=1)
                if growth:
                    plt.ylabel(r'$\gamma$ $[s^{-1}]$')
                else:
                    plt.ylabel(r'Kinetic energy $[J]$')
            elif units.lower()=='m3dc1':
                if growth:
                    plt.ylabel(r'$\gamma/\omega_A$')
                else:
                    plt.ylabel(r'Kinetic energy (M3DC1 units)')
        else:
            plt.ylabel(ylbl,fontsize=axlblfs)
    
    if export:
        plot_data = np.column_stack([time])
    
    if len(scalar.shape)>1:
        leglabels = ["n={:2}".format(n) for n in range(scalar.shape[1])]
        for i in range(scalar.shape[1]):
            if not (skip_n0 and i==0):
                plt.plot(time,scalar[:,i],lw=linew,label=leglabels[i])
        
        ncol = 2 if scalar.shape[1] > 5 else 1
        if show_legend:
            plt.legend(loc=0, ncol=ncol,fontsize=legfs)
        else:
            labelLines(plt.gca().get_lines(), zorder=2.5)
    else:
        plt.plot(time,scalar,lw=linew,label=leglbl)
        if show_legend and leglbl is not None:
            plt.legend(loc=0,fontsize=legfs)
            
    if export:
        plot_data = np.column_stack([plot_data, scalar])
    plt.grid(True)
    #Determine y-axis limits
    if rescale:
        if np.amax(scalar[1:]) < scalar[0]:
            start_time=250
            if units=='mks':
                start_time = unit_conv(start_time,arr_dim='m3dc1',sim=sim,time=1)
            start_ind = int(fpyl.find_nearest(time,start_time))
            top_lim=1.1*np.amax(scalar[start_ind:])
            plt.ylim([0,top_lim])
    
    #Plot vertical lines to mark certain points in time specified by time_marks.
    if len(time_marks)>0:
        for t in time_marks:
            plt.axvline(x=t,c='m')
    
    if ts_marks_all:
        slices = glob.glob('time*.h5')
        slices.sort(key=os.path.getmtime)
        temp = [s.replace('time_', '') for s in slices]
        temp = [s.replace('.h5', '') for s in temp]
        ts_marks = np.array(temp,dtype=int)
    
    if len(ts_marks)>0:
        for ts in ts_marks:
            ts_marks_time = get_time_of_slice(ts,filename=filename,units=units,millisec=millisec)
            plt.axvline(x=ts_marks_time,c='m')
    
    ax = plt.gca()
    if in_plot_txt is not None:
        plt.text(0.03, 0.95,in_plot_txt, ha='left', va='top', transform=ax.transAxes,fontsize=inplttxtfs)
    
    plt.yscale(yscale)
    ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
    if pub:
        yaxis_magn = ax.yaxis.get_offset_text()
        yaxis_magn.set_size(ticklblfs)
    if yscale == 'linear':
        plt.ticklabel_format( axis='y', style='sci',useOffset=False)
    if title:
        plt.title('n='+str(ntor),fontsize=titlefs)
    plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    plt.show()
    
    if save:
        tracestr = trace
        if growth:
            tracestr = tracestr + '-growth'
        if savedir is not None:
            tracestr = savedir + tracestr
        plt.savefig(tracestr+'_n'+"{:d}".format(ntor)+'.pdf', format='pdf',bbox_inches='tight')
        
    if export:
        print(plot_data)
        np.savetxt(txtname,plot_data,delimiter='   ')
    return
