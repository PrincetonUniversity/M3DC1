#!/usr/bin/env python3
#
# Coded on August 16th 2019 by:
# Andreas Kleiner:    akleiner@pppl.gov
# Ralf Mackenbach:    rmackenb@pppl.gov
# Chris Smiet    :    csmiet@pppl.gov

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.ticker as ticker
from m3dc1.unit_conv import unit_conv
from m3dc1.eval_field import eval_field
from m3dc1.plot_mesh import plot_mesh
import m3dc1.fpylib as fpyl
rc('text', usetex=True)


def plot_field_basic(field, coord='scalar', filename='C1.h5', sim=None, time=None,
                     phi=0, mesh=False, linear=False, diff=False, tor_av=1,
                     bound=False, units='mks', res=250, quiet=False):
    """
    Plots the field of a file. 
    
    Arguments:

    **field**
    The field that is to be plotted, i.e. 'B', 'j', etc.

    **coord**
    The chosen part of a vector field to be plotted, options are:
    'phi', 'R', 'Z', 'poloidal', 'radial', 'scalar'.
    Scalar is reserved for scalar fields.
    Poloidal takes the magnetic axis at time=0 as the origin, and
    defines anti-clockwise as the positive direction.
    Radial also takes the magnetic axis as the origin.

    **filename**
    File name which will be read, i.e. "../../C1.h5"
    Can also be a list of two filepaths when used for diff

    **time**
    The time-slice which will be used for the field plot

    **phi**
    The toroidal cross-section coordinate.

    **mesh**
    Overlay mesh on top of plot. 
    True/False

    **bound**
    Only plot boundary. Only works when mesh is true
    True/False

    **linear**
    Plot the linear part of the field (so the equilibrium is subtracted).
    True/False
    
    **diff**
    Plot the difference of two fields. 
    This could be the difference of two files (filename=['a/C1.h5','b/C1.h5']),
    or the difference between two time-slices (time=[t1,t2])
    If list for both time and filename are given file1 will be evaluated at time1,
    and file2 at time2

    **tor_av**
    Calculates the average field over tor_av number of toroidal planes

    **res**
    Resolution in R and Z direction.
    """
    
    sim, time = fpyl.setup_sims(sim,filename,time,linear,diff)
    field_idx = fpyl.get_field_idx(coord)

    # Make 3D grid based on the mesh points
    mesh_ob      = sim[0].get_mesh(quiet=quiet)
    mesh_pts     = mesh_ob.elements
    R_mesh       = mesh_pts[:,4]
    Z_mesh       = mesh_pts[:,5]
    R_range      = [np.amin(R_mesh),np.amax(R_mesh)]
    Z_range      = [np.amin(Z_mesh),np.amax(Z_mesh)]
    R_linspace   = np.linspace(R_range[0], R_range[1], res, endpoint=True)
    phi_linspace = np.linspace(phi,         (360+phi), tor_av, endpoint=False)
    Z_linspace   = np.linspace(Z_range[0], Z_range[1], res, endpoint=True)
    R, phi, Z    = np.meshgrid(R_linspace, phi_linspace,Z_linspace)


   # Get the magnetic axis at time zero, which will be used for poloidal coordinates
    if coord in ['poloidal', 'radial']:
        R_mag =  sim[0].get_time_trace('xmag').values[0]
        Z_mag =  sim[0].get_time_trace('zmag').values[0]


    # Evaluate usual vector components
    if coord not in ['poloidal', 'radial']:
        # Evaluate field
        if not quiet:
            print('Evaluating field... ')
        field1 = eval_field(field, R, phi, Z, coord=coord, sim=sim[0], time=time[0])
    
        # Evaluate second field and calculate difference between two if linear or diff is True
        if diff or linear:
            if not quiet:
                print('Evaluating second field... ')
            field2 = eval_field(field, R, phi, Z, coord=coord, sim=sim[1], time=time[1])
            field1 = field1 - field2


    # Evaluate poloidal component
    if coord == 'poloidal':
        if not quiet:
            print('Evaluating poloidal field... ')
        field1R, field1phi, field1Z  = eval_field(field, R, phi, Z, coord='vector', sim=sim[0], time=time[0])
        theta                        = np.arctan2(Z-Z_mag,R-R_mag)
        field1                       = -np.sin(theta)*field1R + np.cos(theta)*field1Z
        
        # Evaluate second field and calculate difference between two if linear or diff is True
        if diff or linear:
            if not quiet:
                print('Evaluating second poloidal field... ')
            field2R, field2phi, field2Z  = eval_field(field, R, phi, Z, coord='vector', sim=sim[1], time=time[1])
            field2                       = -np.sin(theta)*field2R + np.cos(theta)*field2Z
            field1                       = field1 - field2


    # Evaluate radial component
    if coord == 'radial':
        if not quiet:
            print('Evaluating radial field... ')
        field1R, field1phi, field1Z  = eval_field(field, R, phi, Z, coord='vector', sim=sim[0], time=time[0])
        theta                        = np.arctan2(Z-Z_mag,R-R_mag)
        field1                       = np.cos(theta)*field1R + np.sin(theta)*field1Z
        
        # Evaluate second field and calculate difference between two if linear or diff is True
        if diff or linear:
            if not quiet:
                print('Evaluating second poloidal field... ')
            field2R, field2phi, field2Z  = eval_field(field, R, phi, Z, coord='vector', sim=sim[1], time=time[1])
            field2                       = np.cos(theta)*field2R + np.sin(theta)*field2Z
            field1                       = field1 - field2


    # Calculate average over phi of the field. For a single slice the average is equal to itself
    field1_ave     = np.average(field1, 0)
    R_ave          = np.average(R, 0)
    Z_ave          = np.average(Z, 0)


    if units.lower()=='m3dc1':
        field1_ave = fpyl.get_conv_field(units,field,field1_ave,sim=sim[0])

    fieldlabel,unitlabel = fpyl.get_fieldlabel(units,field,shortlbl=shortlbl)
    if units.lower()=='m3dc1':
        unitlabel = fieldlabel + ' (' + unitlabel + ')'


    # Plot routines
    def fmt(x, pos):
        a, b = '{:.2e}'.format(x).split('e')
        b = int(b)
        return r'${} \cdot 10^{{{}}}$'.format(a, b)


    fig, axs = plt.subplots(1, 1)
    cont = axs.contourf(R_ave, Z_ave, field1_ave,100,cmap='viridis')
    cbar = fig.colorbar(cont,format=ticker.FuncFormatter(fmt),ax=axs)
    cbar.set_label(label)


    if mesh:
        meshplt = plot_mesh(mesh_ob,boundary=bound,ax=axs)
        plt.rcParams["axes.axisbelow"] = False
    
    plt.xlabel(r'$R$')
    plt.ylabel(r'$Z$')
    plt.title('Field plot')
    
    plt.show()
