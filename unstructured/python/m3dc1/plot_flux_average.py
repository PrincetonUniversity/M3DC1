#!/usr/bin/env python3
#
# plot_shape: plots shapes of flux surfaces.
#
# Coded on 08/27/2019 by:
# Andreas Kleiner:    akleiner@pppl.gov
# Ralf Mackenbach:    rmackenb@pppl.gov
# Chris Smiet    :    csmiet@pppl.gov

import fpy
import matplotlib.pyplot as plt
from matplotlib import rc
import m3dc1.fpylib as fpyl
from m3dc1.flux_average import flux_average
#rc('text', usetex=True)

#ToDo: Add rms
def plot_flux_average(field, coord='scalar', sim=None, fcoords='', units='m3dc1', filename='C1.h5', time=0, phit=0, rms=False,pub=False):
    """
    Plots flux surfaces
    
    Arguments:

    **sim**
    simulation sim_data objects. If none is provided, plot_field will read a file and create
    an object.

    **filename**
    File name which will be read, i.e. "../C1.h5"
    Can also be a list of two filepaths when used for diff

    **time**
    The time-slice which will be used for the field plot

    **phit**
    The toroidal cross-section coordinate.
    """
    
    flux, fa = flux_average(field,coord=coord,sim=sim, phit=phit, filename=filename, time=time, fcoords=fcoords, units=units)
    
    # Set font sizes and plot style parameters
    if pub==False:
        axlblfs = 12
        titlefs = 12
        cbarlblfs = None
        cbarticklblfs = None
        ticklblfs = 12
        linew = 1
    elif pub==True:
        axlblfs = 20
        titlefs = 18
        cbarlblfs = 14
        cbarticklblfs = 14
        ticklblfs = 18
        linew = 2
    
    plt.figure()
    plt.plot(flux, fa, lw=linew)
    ax = plt.gca()
    plt.grid(True)
    plt.xlabel(r'$\psi_N$',fontsize=axlblfs)
    fieldlabel,unitlabel = fpyl.get_fieldlabel(units,field)
    ylbl = fieldlabel + ' (' + unitlabel+')'
    plt.ylabel(ylbl,fontsize=axlblfs)
    ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
    plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    
    
    return
