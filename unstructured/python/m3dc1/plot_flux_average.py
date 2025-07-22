#!/usr/bin/env python3
#
# plot_flux_average: plots a flux averaged quantity
#
# Coded on 08/27/2019 by:
# Andreas Kleiner:    akleiner@pppl.gov
# Ralf Mackenbach:    rmackenb@pppl.gov
# Chris Smiet    :    csmiet@pppl.gov

import fpy
import matplotlib.pyplot as plt
#from matplotlib import rc
import m3dc1.fpylib as fpyl
from m3dc1.flux_average import flux_average
#rc('text', usetex=True)

#ToDo: Add rms
def plot_flux_average(field, sim=None, coord='scalar', fcoords='pest', deriv=0, points=200, filename='C1.h5', time=0, device=None, units='m3dc1',
                      fac=1, phit=0, rms=False,pub=False,c=None,ls='-',xlimits=[None,None],ylimits=[None,None],show_legend=False,leglbl=None,
                      shortlbl=False,fignum=None,figsize=None):
    """
    Plots flux surfaces
    
    Arguments:

    **field**
    Name of the field to flux average

    **sim**
    simulation sim_data objects. If none is provided, plot_field will read a file and create
    an object.

    **coord**
    For vector fields, component of field to flux average, e.g. R, phi, Z

    **fcoords**
    Name of desired flux coordinate system : 'pest', 'boozer', 'hamada', canonical, ''

    **sim**
    fpy simulation object

    **deriv**
    If 1, calculate and return derivative of flux-averaged quantity dy/dpsin; if 2, calculate derivate w.r.t. to psi

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

    **c**
    Line color

    **ls**
    line style

    **xlimits**
    x-axis limits, array of length 2, e.g. [0,1]

    **ylimits**
    y-axis limits, array of length 2, e.g. [0,100.0]

    **shortlbl**
    If True, use short y-axis label, e.g. 'T_e' instead of 'electron temperature'

    **fignum**
    Number of figure to plot into

    **figsize**
    Array with length 2, width and height of figure window, e.g. [4.8,2.4]
    """
    
    flux, fa = flux_average(field,coord=coord,sim=sim, deriv=deriv, points=points, phit=phit, filename=filename, time=time, fcoords=fcoords, units=units, device=device)
    
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
    
    plt.figure(num=fignum,figsize=figsize)
    plt.plot(flux, fa*fac, c=c, lw=linew,ls=ls,label=leglbl)
    ax = plt.gca()
    plt.grid(True)
    plt.xlabel(r'$\psi_N$',fontsize=axlblfs)
    
    ax.set_xlim(left=xlimits[0],right=xlimits[1])
    ax.set_ylim(bottom=ylimits[0],top=ylimits[1])
    
    fieldlabel,unitlabel = fpyl.get_fieldlabel(units,field,fac=fac,shortlbl=shortlbl)
    ylbl = fieldlabel
    if field not in ['q','safety factor','ne/ng']:
            if not (field in ['alpha','shear'] and units=='m3dc1') or field not in ['S','lundquist']:
                ylbl = ylbl + ' (' + unitlabel+')'
        
    plt.ylabel(ylbl,fontsize=axlblfs)
    ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
    
    if show_legend:
        if leglbl is not None:
            plt.legend(loc=0,fontsize=legfs,handlelength=leghandlen)
    plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    
    
    return
