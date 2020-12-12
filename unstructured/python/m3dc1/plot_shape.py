#!/usr/bin/env python3
#
# plot_shape: plots shape of flux surfaces.
#
# Coded on 09/11/2019 by:
# Andreas Kleiner:    akleiner@pppl.gov
# Ralf Mackenbach:    rmackenb@pppl.gov
# Chris Smiet    :    csmiet@pppl.gov

import fpy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import m3dc1.fpylib as fpyl
from m3dc1.plot_mesh import plot_mesh
from m3dc1.eval_field import eval_field
#rc('text', usetex=True)



def plot_shape(sim=None, filename='C1.h5', time=None, phi=0,
               Nlvl_in=10, Nlvl_out=1, mesh=False, bound=False, lcfs=False,
               ax=None, pub=False, quiet=False):
    """
    Plot flux surfaces in poloidal plane
    
    Arguments:

    **sim**
    simulation sim_data object. If none is provided, plot_shape will read a file and create
    an object.

    **filename**
    File name which will be read, i.e. "../C1.h5"
    Can also be a list of two filepaths when used for diff

    **time**
    The time-slice which will be used for the shape plot

    **phi**
    The toroidal cross-section coordinate.

    **mesh**
    Overlay mesh on top of the plot. 
    True/False

    **bound**
    Overlay boundary (mesh boundary and walls) without other mesh lines on top of plot.

    **lcfs**
    Overlay last closed flux surface on top of the plot.

    **ax**
    matplotlib axes object to plot into. If None, new figure and axes are created.
    
    **pub**
    If True, plot will be formatted for publication
    """
    if sim is None:
        sim = fpy.sim_data(filename)
    
    psi_axis = sim.get_time_trace('psimin').values[0]
    psi_lcfs = sim.get_time_trace('psi_lcfs').values[0]
    if not quiet:
        print("Psi at magnetic axis: "+str(psi_axis))
        print("Psi at LCFS: "+str(psi_lcfs))
    
    
    fac = 1.0
    elms = sim.get_mesh(time=time,quiet=quiet)
    mp = elms.elements
    R = mp[:,4]*fac
    Z = mp[:,5]*fac
    
    vecpotfield = sim.get_field('A',time)
    Aphi = np.zeros_like(R)
    for idx in range(len(R)):
        Aphi[idx] = vecpotfield.evaluate((R[idx],phi,Z[idx]))[1]
    psifield = R*Aphi
    
    levels = (np.arange(Nlvl_in+Nlvl_out+1.)/Nlvl_in)
    levels = levels*(psi_lcfs-psi_axis)+psi_axis
    if psi_lcfs < psi_axis:
        levels = np.flipud(levels)
    
    # Set font sizes and plot style parameters
    if pub:
        axlblfs = 20
        ticklblfs = 18
        linew = 1
        bdlw = 3
    else:
        axlblfs = 12
        ticklblfs = 12
        linew = 1
        bdlw = 2
    
    if ax is None:
        fig, axs = plt.subplots(1, 1)
        fig.set_figheight(7)
    else:
        axs = ax
    
    if mesh:
        meshplt = plot_mesh(elms,boundary=False,ax=axs)
    elif bound:
        meshplt = plot_mesh(elms,boundary=bound,ax=axs)
    
    cont = axs.tricontour(R, Z, psifield,levels,zorder=0,linewidths=linew)
    
    if lcfs:
        cont = axs.tricontour(R, Z, psifield,[psi_lcfs],colors='magenta',linewidths=bdlw,linestyles='--',zorder=10)
    if ax is None:
        axs.set_xlim(fpyl.get_axlim(np.amin(R),'min'),fpyl.get_axlim(np.amax(R),'max'))
        axs.set_ylim(fpyl.get_axlim(np.amin(Z),'min'),fpyl.get_axlim(np.amax(Z),'max'))
        axs.set_aspect('equal')
        axs.set_xlabel(r'$R$',fontsize=axlblfs)
        axs.set_ylabel(r'$Z$',fontsize=axlblfs)
        axs.tick_params(labelsize=ticklblfs)
        axs.grid(True,zorder=10,alpha=0.5) #There seems to be a bug in matplotlib that ignores the zorder of the grid
        
        plt.rcParams["axes.axisbelow"] = False #Allows mesh to be on top of contour plot. This option conflicts with zorder (matplotlib
        plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
        plt.show()
    return axs
