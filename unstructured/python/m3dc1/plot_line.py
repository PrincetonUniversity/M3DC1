#!/usr/bin/env python3
#
# Coded on August 28th 2019 by:
# Andreas Kleiner:    akleiner@pppl.gov
# Ralf Mackenbach:    rmackenb@pppl.gov
# Chris Smiet    :    csmiet@pppl.gov

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from m3dc1.unit_conv  import unit_conv
from m3dc1.eval_field import eval_field
from m3dc1.plot_mesh import plot_mesh
import m3dc1.fpylib as fpyl
rc('text', usetex=True)




def plot_line(field, coord='scalar', angle=0, Zoff=0, dist_from_magax=False, filename='C1.h5',
              sim=None, time=None, phi=0, linear=False, diff=False, tor_av=1, units='mks',
              c=None,ls='-',shortlbl=False, show_legend=False, leglbl=None, quiet=False, fignum=None, pub=False):
    """
    Plots the values of a field on a line.
    
    Arguments:

    **field**
    The field that is to be plotted, i.e. 'B', or 'j'.

    **coord**
    The chosen part of a vector field to be plotted, options are:
    'phi', 'R', 'Z', 'poloidal', 'radial', 'scalar'.
    Scalar is reserved for scalar fields.
    Poloidal takes the magnetic axis at time=0 as the origin, and
    defines anti-clockwise as the positive direction.
    Radial also takes the magnetic axis as the origin.

    **angle**
    Angle in degrees between the midplane and the line on which the field
    is evaluated. Positive angles are anti-clockwise.

    **Zoff**
    Offset in Z coordinate. By default the field will be plotted along a
    line with angle=angle through the magnetic axis. If Zoff is different
    from zero, this line will be shifted by Zoff.

    **dist_from_magax**
    Determines radial variable for plot. Default: False
    If True, field will be plotted vs. signed distance from magnetic axis.
    If False, field will plotted vs. major radius.

    **sim**
    Simulation object.

    **filename**
    File name which will be read, i.e. "../../C1.h5"
    Can also be a list of two filepaths when used for diff

    **time**
    The time-slice which will be used for the field plot

    **phi**
    The toroidal cross-section coordinate.

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

    **units**
    units in which the result will be calculated

    **shortlbl**
    True/False. If True, axis labels will be shortened.

    **quiet**
    If True, suppress output to terminal.

    **fignum**
    Figure number.

    **pub**
    If True, format figure for publication (larger labels and thicker lines)
    """

    angle = np.radians(angle)

    sim, time = fpyl.setup_sims(sim,filename,time,linear,diff)
    field_idx = fpyl.get_field_idx(coord)

    # Make 3D grid based on the mesh points
    mesh_ob      = sim[0].get_mesh(quiet=quiet)
    mesh_pts     = mesh_ob.elements
    R_mag        = sim[0].get_time_trace('xmag').values[0]
    Z_mag        = sim[0].get_time_trace('zmag').values[0]
    R_range      = [np.amin(mesh_pts[:,4]), np.amax(mesh_pts[:,4])]
    Z_range      = [np.amin(mesh_pts[:,5]), np.amax(mesh_pts[:,5])]
    R_pts        = np.linspace(R_range[0], R_range[1], 10000, endpoint=True)
    R_straight   = R_pts-R_mag
    R_linspace   = (R_pts-R_mag)*np.cos(angle) + R_mag
    Z_linspace   = Zoff + (np.linspace(Z_range[0], Z_range[1], 10000, endpoint=True)-Z_mag)*np.sin(angle) + Z_mag
    phi_linspace = np.linspace(phi,      (360+phi), tor_av, endpoint=False)
    R, phi       = np.meshgrid(R_linspace, phi_linspace)
    Z, phi       = np.meshgrid(Z_linspace, phi_linspace)
    if dist_from_magax:
        xvar = np.sign(R_straight)*np.sqrt((R_linspace-R_mag)**2 + (Z_linspace-Z_mag)**2)
    else:
        xvar = R_linspace


    # Evaluate usual vector components
    if coord not in ['poloidal', 'radial']:
        # Evaluate field
        field1 = eval_field(field, R, phi, Z, coord=coord, sim=sim[0], time=time[0],quiet=quiet)
    
        # Evaluate second field and calculate difference between two if linear or diff is True
        if diff or linear:
            field2 = eval_field(field, R, phi, Z, coord=coord, sim=sim[1], time=time[1],quiet=quiet)
            field1 = field1 - field2


    # Evaluate poloidal component
    if coord == 'poloidal':
        field1R, field1phi, field1Z  = eval_field(field, R, phi, Z, coord='vector', sim=sim[0], time=time[0],quiet=quiet)
        theta                        = np.arctan2(Z-Z_mag,R-R_mag)
        field1                       = -np.sin(theta)*field1R + np.cos(theta)*field1Z
        
        # Evaluate second field and calculate difference between two if linear or diff is True
        if diff or linear:
            field2R, field2phi, field2Z  = eval_field(field, R, phi, Z, coord='vector', sim=sim[1], time=time[1],quiet=quiet)
            field2                       = -np.sin(theta)*field2R + np.cos(theta)*field2Z
            field1                       = field1 - field2


    # Evaluate radial component
    if coord == 'radial':
        field1R, field1phi, field1Z  = eval_field(field, R, phi, Z, coord='vector', sim=sim[0], time=time[0],quiet=quiet)
        theta                        = np.arctan2(Z-Z_mag,R-R_mag)
        field1                       = np.cos(theta)*field1R + np.sin(theta)*field1Z
        
        # Evaluate second field and calculate difference between two if linear or diff is True
        if diff or linear:
            field2R, field2phi, field2Z  = eval_field(field, R, phi, Z, coord='vector', sim=sim[1], time=time[1],quiet=quiet)
            field2                       = np.cos(theta)*field2R + np.sin(theta)*field2Z
            field1                       = field1 - field2


    # Calculate average over phi of the field. For a single slice the average is equal to itself
    field1_ave = np.average(field1, 0)
    R_ave = np.average(R, 0)
    Z_ave = np.average(Z, 0)


    if units.lower()=='m3dc1':
        field1_ave = fpyl.get_conv_field(units,field,field1_ave,sim=sim[0])

    fieldlabel,unitlabel = fpyl.get_fieldlabel(units,field,shortlbl=shortlbl)
    ylbl = fieldlabel + ' (' + unitlabel + ')'

    if pub:
        axlblfs = 20
        #titlefs = 20
        ticklblfs = 18
        linew = 2
        legfs = 12
    else:
        axlblfs = None
        #titlefs = None
        ticklblfs = None
        linew = 1
        legfs = None

    # Plot routines
    plt.figure(num=fignum)
    if ls==':' and linew is not None:
        linew = linew+1
    plt.plot(xvar, field1_ave,lw=linew,label=leglbl,c=c,ls=ls)
    if dist_from_magax:
        xlbl = 'Signed distance from magnetic axis [m]'
    else:
        xlbl = 'R / m'
    plt.xlabel(xlbl,fontsize=axlblfs)
    plt.ylabel(ylbl,fontsize=axlblfs)
    if show_legend and leglbl is not None:
        plt.legend(loc=0, fontsize=legfs)
    ax = plt.gca()
    ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
    plt.grid(True)
    plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    plt.show()

    return
