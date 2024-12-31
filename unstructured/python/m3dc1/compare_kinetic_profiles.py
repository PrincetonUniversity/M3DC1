#!/usr/bin/env python3
#
# compare_kinetic_profiles: plots total p, electron p and ion p
#
# Coded on 08/27/2019 by:
# Andreas Kleiner:    akleiner@pppl.gov
# Ralf Mackenbach:    rmackenb@pppl.gov
# Chris Smiet    :    csmiet@pppl.gov

import numpy as np
import fpy
import matplotlib.pyplot as plt
#from matplotlib import rc
import m3dc1.fpylib as fpyl
from m3dc1.flux_average import flux_average
from m3dc1.flux_coordinates import flux_coordinates
#rc('text', usetex=True)

def compare_kinetic_profiles(sim=None, fcoords='pest', deriv=0, points=200, filename='C1.h5', time=0, units='mks',
                             fac=1, phit=0, pub=False, ls='-', xlimits=[None,None], ylimits=[None,None],
                             show_legend=True, leglbl=None, fignum=None, figsize=None):
    """
    Plots flux average of total pressure, electron pressure and ion pressure
    
    Arguments:

    **sim**
    simulation sim_data objects. If none is provided, plot_field will read a file and create
    an object.

    **fcoords**
    Name of desired flux coordinate system : 'pest', 'boozer', 'hamada', canonical, ''

    **deriv**
    If 1, calculate and return derivative of flux-averaged quantity dy/dpsin; if 2, calculate derivate w.r.t. to psi

    **points**
    Number of flux surfaces between psin = 0 to 1, where flux average is calculated.

    **filename**
    File name which will be read, i.e. "../C1.h5"
    Can also be a list of two filepaths when used for diff

    **time**
    The time-slice which will be used for the flux average

    **phit**
    The toroidal cross-section coordinate.

    **units**
    Units in which the result will be calculated

    **fac**
    Factor to apply to field when using mks units. fac=1.0E-3 converts to kilo..., fac=1.0E-6 to Mega...

    **pub**
    If True, format figure for publication (larger labels and thicker lines)

    **ls**
    line style

    **show_legend**
    True/False. Display plot legend.

    **leglbl**
    Legend label.

    **xlimits**
    x-axis limits, array of length 2, e.g. [0,1]

    **ylimits**
    y-axis limits, array of length 2, e.g. [0,100.0]

    **fignum**
    Number of figure to plot into

    **figsize**
    Array with length 2, width and height of figure window, e.g. [4.8,2.4]
    """
    
    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename=filename,time=time)
    else:
        if fcoords is None and sim.fc is not None:
            fcoords = sim.fc.fcoords
    
    # Calculate flux coodinates if it was not calculated yet or a different flux coordinate system than sim.fc.fcoords is desired
    if (not isinstance(sim.fc,fpy.flux_coordinates)) or (fcoords is not None and (sim.fc.fcoords!=fcoords)) or (sim.fc.points!=points):
        if fcoords is None:
            fcoords = ''
            print("FCOORDS SET TO NOTHING")
        sim = flux_coordinates(sim=sim, fcoords=fcoords, filename=filename, time=time, points=points, phit=phit,psin_range=None)
    
    flux_ptot, ptot = flux_average('p',coord='scalar',sim=sim, deriv=deriv, points=points, phit=phit, filename=filename, time=time, fcoords=fcoords, units=units)
    flux_ne, ne = flux_average('ne',coord='scalar',sim=sim, deriv=deriv, points=points, phit=phit, filename=filename, time=time, fcoords=fcoords, units=units)
    flux_te, te = flux_average('te',coord='scalar',sim=sim, deriv=deriv, points=points, phit=phit, filename=filename, time=time, fcoords=fcoords, units=units)
    flux_ni, ni = flux_average('ni',coord='scalar',sim=sim, deriv=deriv, points=points, phit=phit, filename=filename, time=time, fcoords=fcoords, units=units)
    flux_ti, ti = flux_average('ti',coord='scalar',sim=sim, deriv=deriv, points=points, phit=phit, filename=filename, time=time, fcoords=fcoords, units=units)
    
    if not np.all(flux_ne==flux_te):
        fpyl.printwarn('WARNING: fluxes for ne and te are not identical!')
    if not np.all(flux_ni==flux_ti):
        fpyl.printwarn('WARNING: fluxes for ni and ti are not identical!')
    
    norm = 11604.51812*1.380649E-23 if units=='mks' else 1.0
    pe = ne*te*norm
    pi = ni*ti*norm
    
    # Set font sizes and plot style parameters
    if pub:
        axlblfs = 20
        #titlefs = 18
        ticklblfs = 18
        if ls!=':':
            linew = 3
        else:
            linew = 4
        legfs = 14
        leghandlen = 3.0
    else:
        axlblfs = 12
        #titlefs = 12
        ticklblfs = 12
        linew = 1
        legfs = None
        leghandlen = 2.0
    
    # Plot pressure
    plt.figure(num=fignum,figsize=figsize)
    plt.plot(flux_ptot, ptot*fac, lw=linew,ls=ls,label='total pressure')
    plt.plot(flux_ne, pe*fac, lw=linew,ls=ls,label='electron pressure')
    plt.plot(flux_ni, pi*fac, lw=linew,ls=ls,label='ion pressure')
    ax = plt.gca()
    plt.grid(True)
    plt.xlabel(r'$\psi_N$',fontsize=axlblfs)
    
    ax.set_xlim(left=xlimits[0],right=xlimits[1])
    ax.set_ylim(bottom=ylimits[0],top=ylimits[1])
    
    ylbl = 'pressure'
    if units=='mks':
        ylbl = ylbl + ' (Pa)'
    plt.ylabel(ylbl,fontsize=axlblfs)
    #plt.ylabel(ylbl,fontsize=axlblfs)
    ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
    
    if show_legend:
        plt.legend(loc=0,fontsize=legfs,handlelength=leghandlen)
    plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    
    
    fig2num = plt.gcf().number+1 if fignum is None else fignum+1
    
    # Plot density
    plt.figure(num=fig2num,figsize=figsize)
    #plt.plot(flux_ni, ne+ni, lw=linew,ls=ls,label='total density')
    plt.plot(flux_ne, ne, lw=linew,ls=ls,label='electron density')
    plt.plot(flux_ni, ni, lw=linew,ls=ls,label='ion density')
    ax = plt.gca()
    plt.grid(True)
    plt.xlabel(r'$\psi_N$',fontsize=axlblfs)
    ylbl = 'density'
    plt.ylabel(ylbl,fontsize=axlblfs)
    ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
    if show_legend:
        plt.legend(loc=0,fontsize=legfs,handlelength=leghandlen)
    plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    
    # Plot temperature
    plt.figure(num=fig2num+1,figsize=figsize)
    plt.plot(flux_te, te, lw=linew,ls=ls,label='electron temperature')
    plt.plot(flux_ti, ti, lw=linew,ls=ls,label='ion temperature')
    ax = plt.gca()
    plt.grid(True)
    plt.xlabel(r'$\psi_N$',fontsize=axlblfs)
    ylbl = 'temperature'
    plt.ylabel(ylbl,fontsize=axlblfs)
    if units=='mks':
        ylbl = ylbl + ' (eV)'
    ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
    if show_legend:
        plt.legend(loc=0,fontsize=legfs,handlelength=leghandlen)
    plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    
    return
