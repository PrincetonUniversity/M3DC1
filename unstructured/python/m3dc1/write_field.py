#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  26 18:31:00 2021

@author: Andreas Kleiner
"""
#import os
#import glob
import numpy as np
import h5py
import csv
#import math
import m3dc1.fpylib as fpyl
from m3dc1.get_field import get_field


def write_field(field, coord='scalar', row=1, sim=None, filename='C1.h5', time=None, phi=0, linear=False,
               diff=False, tor_av=1, units='mks', res=250, quiet=False, phys=False, suffix='', filetype='text'):
    """
    Writes field values and corresponding grid values (R,phi,Z) to a file. File type can be chosen.
    
    Arguments:

    **field**
    The field that is to be plotted, i.e. 'B', 'j', etc.

    **coord**
    The chosen part of a vector field to be plotted, options are:
    'phi', 'R', 'Z', 'poloidal', 'radial', 'scalar', 'vector', 'tensor'.
    Scalar is reserved for scalar fields.
    Poloidal takes the magnetic axis at time=0 as the origin, and
    defines anti-clockwise as the positive direction.
    Radial also takes the magnetic axis as the origin.

    **row**
    For tensor fields only, row of tensor to plot. Possibilities: 1,2,3

    **sim**
    simulation sim_data object or list of sim_data objects. If none is provided, the object will be created.

    **filename**
    File name which will be read, i.e. "../../C1.h5"
    Can also be a list of two filepaths when used for diff

    **time**
    The time-slice which will be used for the field plot. If time='last', the last time slice will be used.

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
    Units in which the field will be returned.

    **res**
    Resolution in R and Z direction.

    **phys**
    Use True for plotting in physical (stellarator) geometry

    **quiet**
    If True, suppress output to terminal.

    **suffix**
    Suffix added to file name

    **filetype**
    Type of output file. Currently supported: txt, csv, hdf5
    """
    
    sim, time, mesh_ob, R, phi_list, Z, R_mesh, Z_mesh, R_ave, Z_ave, field1_ave = get_field(field=field, coord=coord, row=row, sim=sim, filename=filename, time=time, phi=phi, linear=linear,
               diff=diff, tor_av=tor_av, units=units, res=res, quiet=quiet, phys=phys)
    
    if len(time)>1:
        ts = np.amax(time)
    else:
        ts = sim[0].timeslice
    outfile = 'field_'+str(field)
    if coord not in ['scalar','vector']:
        outfile = outfile + '_' + coord
    outfile = outfile+'_t'+str(np.amax(ts))+suffix

    # Write field to file:
    if filetype in ['text','txt','dat']:
        #for i in zip(columns):
        #    print(i)
        with open(outfile+'.txt', 'w') as out:
            for (R,p,Z,fi) in zip(R_ave.flatten(), phi_list[0].flatten(), Z_ave.flatten(), field1_ave[0].flatten()):
                #print(R,p,Z,fi)
                out.write(" {0:8.6f}   {1:8.6f}   {2:8.6f}   {3:g}\n".format(R,p,Z,fi))
    elif filetype == 'csv':
        with open(outfile+'.csv', 'w') as f:
            out = csv.writer(f, delimiter=';')
            out.writerows(zip(R_ave.flatten(), phi_list[0].flatten(), Z_ave.flatten(), field1_ave[0].flatten()))
    elif filetype in ['hdf5','h5']:
        hf = h5py.File(outfile+'.h5','w')
        hf.create_dataset('R', data=R_ave)
        hf.create_dataset('phi', data=phi_list[0])
        hf.create_dataset('Z', data=Z_ave)
        hf.create_dataset('field', data=field1_ave[0])
    else:
        fpyl.printerr('ERROR: file type '+ str(filetype) +' not yet supported!')
    return
