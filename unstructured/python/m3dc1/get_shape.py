#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on July 07 2021

@author: Andreas Kleiner
"""
import numpy as np
import matplotlib.pyplot as plt
import fpy
from matplotlib import path
from m3dc1.eval_field import eval_field

def get_shape(sim,res=250,quiet=False):
    """
    Returns shaping parameters such as minor and major radius. The calculations is based on the countour of the last
    closed flux surface.
    
    Arguments:

    **sim**
    simulation sim_data object. Can also be list of such objects. If None is provided, plot_shape will read a file and create
    an object.

    **res**
    Resolution in R and Z direction.

    **quiet**
    If true, suppress some output to terminal.
    """
    lcfslw=1
    
    mesh_pts      = sim.get_mesh(quiet=quiet).elements
    R_mesh,Z_mesh = (mesh_pts[:,4],mesh_pts[:,5])
    R_range      = [np.amin(R_mesh),np.amax(R_mesh)]
    Z_range      = [np.amin(Z_mesh),np.amax(Z_mesh)]
    R_linspace   = np.linspace(R_range[0], R_range[1], res, endpoint=True)
    phi_linspace = np.linspace(0,         (360+0), 1, endpoint=False)
    Z_linspace   = np.linspace(Z_range[0], Z_range[1], res, endpoint=True)
    R, phi, Z    = np.meshgrid(R_linspace, phi_linspace,Z_linspace)
    R_ave = np.average(R, 0)
    Z_ave = np.average(Z, 0)

    psi_lcfs = sim.get_time_trace('psi_lcfs').values[0] #ToDo: BUG: get values of psi at actual time, c.f. flux coordinates
    R_magax = sim.get_time_trace('xmag').values[0] #ToDo: BUG: get values of psi at actual time, c.f. flux coordinates
    Z_magax = sim.get_time_trace('zmag').values[0] #ToDo: BUG: get values of psi at actual time, c.f. flux coordinates
    
    psifield = eval_field('psi', R, phi, Z, coord='scalar', sim=sim)
    print("Psi at LCFS: "+str(psi_lcfs))
    plt.figure(1729)
    cont = plt.contour(R_ave, Z_ave, np.average(psifield,0),[psi_lcfs],colors='magenta',linewidths=lcfslw,zorder=10)
    plt.close(1729)
    paths = cont.collections[0].get_paths()
    n_points = len(paths)
    for i in range(n_points):
        pp = paths[i]
        vert = pp.vertices
        v = path.Path(pp.vertices)
        if v.contains_point(point=(R_magax,Z_magax)) == True:
            #plt.plot(vert[:,0],vert[:,1],lw=2)
            Rmax = np.amax(vert[:,0])
            Rmin = np.amin(vert[:,0])
            print(Rmax,Rmin)
    a=(Rmax-Rmin)/2
    R0 = (Rmax+Rmin)/2
    if not quiet:
        print('a = '+str(a))
        print('R0 = '+str(R0))
        print('A=R0/a = '+str(R0/a))
    return a,R0
