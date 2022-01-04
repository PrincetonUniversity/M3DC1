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
def plot_flux_average(field, sim=None, coord='scalar', fcoords='', deriv=0, points=200, filename='C1.h5', time=0, units='m3dc1', phit=0, rms=False,pub=False,xlimit=0.0,ylimit=0.0, shortlbl=False,fignum=None):
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

    **pub**
    If True, format figure for publication (larger labels and thicker lines)

    **shortlbl**
    If True, use short y-axis label, e.g. 'T_e' instead of 'electron temperature'

    **fignum**
    Number of figure to plot into
    """
    
    flux, fa = flux_average(field,coord=coord,sim=sim, deriv=deriv, points=points, phit=phit, filename=filename, time=time, fcoords=fcoords, units=units)
    
    # Set font sizes and plot style parameters
    if pub:
        axlblfs = 20
        #titlefs = 18
        ticklblfs = 18
        linew = 2
    else:
        axlblfs = 12
        #titlefs = 12
        ticklblfs = 12
        linew = 1
    
    plt.figure(num=fignum)
    plt.plot(flux, fa, lw=linew)
    ax = plt.gca()
    plt.grid(True)
    plt.xlabel(r'$\psi_N$',fontsize=axlblfs)
    if xlimit>0:
        plt.xlim(left=xlimit,right=1.0)
    if ylimit>0:
        plt.ylim(top=ylimit)
    fieldlabel,unitlabel = fpyl.get_fieldlabel(units,field,shortlbl=shortlbl)
    ylbl = fieldlabel + ' (' + unitlabel+')'
    plt.ylabel(ylbl,fontsize=axlblfs)
    ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
    plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    
    
    return
