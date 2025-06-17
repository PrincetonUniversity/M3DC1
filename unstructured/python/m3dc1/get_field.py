#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  01 2022

@author: Andreas Kleiner
"""
import numpy as np
import m3dc1.fpylib as fpyl
from m3dc1.eval_field import eval_field
from fpy import fio_py


def get_field(field, coord='scalar', row=1, sim=None, filename='C1.h5', time=None, phi=0, linear=False,
               diff=False, scale=1.0, tor_av=1, units='mks', res=250, quiet=False, phys=False, R_range=None, Z_range=None):
    """
    Returns field values on an equidistant rectangular R,phi,Z grid.
    
    Arguments:

    **field**
    Name of the field to be returned, i.e. 'B', 'j', etc.

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
    """
    
    sim, time = fpyl.setup_sims(sim,filename,time,linear,diff)
    #field_idx = fpyl.get_field_idx(coord)
    if phys:
        fio_py.set_quiet_option(True)

    # Make 3D grid based on the mesh points
    mesh_ob      = sim[0].get_mesh(quiet=quiet)
    mesh_pts     = mesh_ob.elements
    R_mesh       = mesh_pts[:,4]
    Z_mesh       = mesh_pts[:,5]
    phi0 = phi*1

    if (R_range is None) or (Z_range is None):
        if phys:
            # rst = eval_field('rst', R, phi, Z, sim=sim, filename=filename, time=time)
            # zst = eval_field('zst', R, phi, Z, sim=sim, filename=filename, time=time)
            R_mesh = eval_field('rst', mesh_pts[:,4], phi0*np.ones_like(mesh_pts[:,4]), mesh_pts[:,5], sim=sim[0], filename=filename, time=-1, quiet=quiet)

        R_range      = [np.nanmin(R_mesh),np.nanmax(R_mesh)]
        if phys and not quiet:
            print("R_range=",R_range)
    if (R_range is None) or (Z_range is None):
        if phys:
            Z_mesh = eval_field('zst', mesh_pts[:,4], phi0*np.ones_like(mesh_pts[:,4]), mesh_pts[:,5], sim=sim[0], filename=filename, time=-1, quiet=quiet)
        Z_range      = [np.nanmin(Z_mesh),np.nanmax(Z_mesh)]
        if phys and not quiet:
            print("Z_range=",Z_range)
     # R_range=[1.1,1.3]
    # Z_range=[-0.35,0.35]
    R_linspace   = np.linspace(R_range[0], R_range[1], res, endpoint=True)
    phi_linspace = np.linspace(phi,         (360+phi), tor_av, endpoint=False)
    Z_linspace   = np.linspace(Z_range[0], Z_range[1], res, endpoint=True)
    R, phi, Z    = np.meshgrid(R_linspace, phi_linspace,Z_linspace)
    # if phys:
    #     rst = eval_field('rst', R, phi, Z, sim=sim, filename=filename, time=-1, quiet=quiet)
    #     zst = eval_field('zst', R, phi, Z, sim=sim, filename=filename, time=-1, quiet=quiet)
    #     R_mesh = eval_field('rst', mesh_pts[:,4], phi0*np.ones_like(mesh_pts[:,4]), mesh_pts[:,5], sim=sim, filename=filename, time=time)
    #     Z_mesh = eval_field('zst', mesh_pts[:,4], phi0*np.ones_like(mesh_pts[:,4]), mesh_pts[:,5], sim=sim, filename=filename, time=time)


    # Get the magnetic axis at time zero, which will be used for poloidal coordinates
    if coord in ['poloidal', 'radial']:
        R_mag =  sim[0].get_time_trace('xmag').values[0]
        Z_mag =  sim[0].get_time_trace('zmag').values[0]


    # Evaluate usual vector components
    if coord not in ['poloidal', 'radial', 'vector', 'tensor']:
        # Evaluate field
        field1 = eval_field(field, R, phi, Z, coord=coord, sim=sim[0], time=time[0],quiet=quiet)
        # Evaluate second field and calculate difference between two if linear or diff is True
        if diff or linear:
            field2 = eval_field(field, R, phi, Z, coord=coord, sim=sim[1], time=time[1],quiet=quiet)
            field1 = field1 - field2
            if not quiet:
                print('[DONE]')


    # Evaluate poloidal/radial field components or all field components (coord='vector')
    if coord in ['poloidal', 'radial', 'vector']:
        field1R, field1phi, field1Z  = eval_field(field, R, phi, Z, coord='vector', sim=sim[0], time=time[0],quiet=quiet)
        if diff or linear:
            field2R, field2phi, field2Z  = eval_field(field, R, phi, Z, coord='vector', sim=sim[1], time=time[1],quiet=quiet)
            if not quiet:
                print('[DONE]')
    elif coord =='tensor':
        field1RR, field1phiR, field1ZR, field1Rphi, field1phiphi, field1Zphi, field1RZ, field1phiZ, field1ZZ  = \
            eval_field(field, R, phi, Z, coord='tensor', sim=sim[0], time=time[0],quiet=quiet)
        if row == 1:
            field1R = field1RR
            field1phi = field1phiR
            field1Z = field1ZR
        elif row == 2:
            field1R = field1Rphi
            field1phi = field1phiphi
            field1Z = field1Zphi
        elif row == 3:
            field1R = field1RZ
            field1phi = field1phiZ
            field1Z = field1ZZ
        if diff or linear:
            field2RR, field2phiR, field2ZR, field2Rphi, field2phiphi, field2Zphi, field2RZ, field2phiZ, field2ZZ  = \
                eval_field(field, R, phi, Z, coord='tensor', sim=sim[1], time=time[1],quiet=quiet)
            if row == 1:
                field2R = field2RR
                field2phi = field2phiR
                field2Z = field2ZR
            elif row == 2:
                field2R = field2Rphi
                field2phi = field2phiphi
                field2Z = field2Zphi
            elif row == 3:
                field2R = field2RZ
                field2phi = field2phiZ
                field2Z = field2ZZ

    # Evaluate poloidal component
    if coord == 'poloidal':
        theta  = np.arctan2(Z-Z_mag,R-R_mag)
        field1 = -np.sin(theta)*field1R + np.cos(theta)*field1Z
        
        # Evaluate second field and calculate difference between two if linear or diff is True
        if diff or linear:
            field2 = -np.sin(theta)*field2R + np.cos(theta)*field2Z


    # Evaluate radial component
    if coord == 'radial':
        theta  = np.arctan2(Z-Z_mag,R-R_mag)
        field1 = np.cos(theta)*field1R + np.sin(theta)*field1Z
        
        # Evaluate second field and calculate difference between two if linear or diff is True
        if diff or linear:
            field2 = np.cos(theta)*field2R + np.sin(theta)*field2Z


    # Calculate difference between two if linear or diff is True
    if diff or linear:
        if coord in ['poloidal', 'radial']:
            field1 = field1 - field2
        elif coord in ['vector', 'tensor']:
            field1R   = field1R - field2R
            field1phi = field1phi - field2phi
            field1Z   = field1Z - field2Z

    
    # Calculate average over phi of the field. For a single slice the average is equal to itself
    if coord in ['vector', 'tensor']:
        field1R_ave   = np.average(field1R, 0)
        field1phi_ave = np.average(field1phi, 0)
        field1Z_ave   = np.average(field1Z, 0)
        field1_ave    = [field1R_ave,field1phi_ave,field1Z_ave]
    else:
        field1_ave    = [np.average(field1, 0)]
    R_ave = np.average(R, 0)
    Z_ave = np.average(Z, 0)
    # if phys:
    #     rst_ave    = np.average(rst, 0)
    #     zst_ave    = np.average(zst, 0)
    #     R_ave = np.where(np.isnan(rst_ave), R_ave, rst_ave)
    #     Z_ave = np.where(np.isnan(zst_ave), Z_ave, zst_ave)
    
    #if units.lower()=='m3dc1':
    #    field1_ave = fpyl.get_conv_field(units,field,field1_ave,sim=sim[0])
    #All fields listed in sim.available_fields are returned in MKS units by fusion-io. Other (scalar) fields are returned in m3dc1 units.
    if (field in sim[0].available_fields and units.lower()=='m3dc1') or (field not in sim[0].available_fields and units.lower()=='mks'):
        field1_ave = fpyl.get_conv_field(units,field,field1_ave,sim=sim[0])

    if phys:
        fio_py.set_quiet_option(False)
    
    return sim, time, mesh_ob, R, phi, Z, R_mesh, Z_mesh, R_ave, Z_ave, np.asarray(field1_ave)




def get_field_vs_phi(field, coord='scalar', row=1, sim=None, filename='C1.h5', time=None, linear=False,
               diff=False, cutr=None, cutz=None, res=250, phi_res=100, units='mks', quiet=False, phys=False):
    """
    Returns field values on an equidistant rectangular R,phi,Z grid.
    
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

    **linear**
    Plot the linear part of the field (so the equilibrium is subtracted).
    True/False

    **diff**
    Plot the difference of two fields. 
    This could be the difference of two files (filename=['a/C1.h5','b/C1.h5']),
    or the difference between two time-slices (time=[t1,t2])
    If list for both time and filename are given file1 will be evaluated at time1,
    and file2 at time2

    **cutr**
    Value of R, where field is plotted as function of phi and Z.

    **cutz**
    Value of Z, where field is plotted as function of R and phi.

    **res**
    Resolution in R and Z direction.

    **phi_res**
    Resolution in phi direction.

    **units**
    Units in which the field will be returned.

    **phys**
    Use True for plotting in physical (stellarator) geometry

    **quiet**
    If True, suppress output to terminal.
    """
    
    sim, time = fpyl.setup_sims(sim,filename,time,linear,diff)
    #field_idx = fpyl.get_field_idx(coord)

    # Make 3D grid based on the mesh points
    mesh_ob      = sim[0].get_mesh(quiet=quiet)
    mesh_pts     = mesh_ob.elements
    R_mesh       = mesh_pts[:,4]
    Z_mesh       = mesh_pts[:,5]

    R_range      = [np.amin(R_mesh),np.amax(R_mesh)]
    Z_range      = [np.amin(Z_mesh),np.amax(Z_mesh)]
    
    phi_linspace = np.linspace(0,2*np.pi, phi_res, endpoint=True)
    
    if cutz is not None:
        R_linspace   = np.linspace(R_range[0], R_range[1], res, endpoint=True)
    else:
        if cutr is not None:
            R_linspace   = np.linspace(cutr, cutr, 1, endpoint=False)
    
    if cutr is not None:
        Z_linspace   = np.linspace(Z_range[0], Z_range[1], res, endpoint=True)
    else:
        if cutz is not None:
            Z_linspace   = np.linspace(cutz, cutz, 1, endpoint=False)
    
    R, phi, Z    = np.meshgrid(R_linspace, phi_linspace,Z_linspace)
    
    #This needs to be rewritten:
    #if phys:
    #    rst = eval_field('rst', R, phi, Z, sim=sim, filename=filename, time=time)
    #    zst = eval_field('zst', R, phi, Z, sim=sim, filename=filename, time=time)
    #    R_mesh = eval_field('rst', mesh_pts[:,4], phi0*np.ones_like(mesh_pts[:,4]), mesh_pts[:,5], sim=sim, filename=filename, time=time)
    #    Z_mesh = eval_field('zst', mesh_pts[:,4], phi0*np.ones_like(mesh_pts[:,4]), mesh_pts[:,5], sim=sim, filename=filename, time=time)


    # Get the magnetic axis at time zero, which will be used for poloidal coordinates
    if coord in ['poloidal', 'radial']:
        R_mag =  sim[0].get_time_trace('xmag').values[0]
        Z_mag =  sim[0].get_time_trace('zmag').values[0]


    # Evaluate usual vector components
    if coord not in ['poloidal', 'radial', 'vector', 'tensor']:
        # Evaluate field
        field1 = eval_field(field, R, phi, Z, coord=coord, sim=sim[0], time=time[0],quiet=quiet)
        # Evaluate second field and calculate difference between two if linear or diff is True
        if diff or linear:
            field2 = eval_field(field, R, phi, Z, coord=coord, sim=sim[1], time=time[1],quiet=quiet)
            field1 = field1 - field2
            if not quiet:
                print('[DONE]')


    # Evaluate poloidal/radial field components or all field components (coord='vector')
    if coord in ['poloidal', 'radial', 'vector']:
        field1R, field1phi, field1Z  = eval_field(field, R, phi, Z, coord='vector', sim=sim[0], time=time[0],quiet=quiet)
        if diff or linear:
            field2R, field2phi, field2Z  = eval_field(field, R, phi, Z, coord='vector', sim=sim[1], time=time[1],quiet=quiet)
            if not quiet:
                print('[DONE]')
    elif coord =='tensor':
        field1RR, field1phiR, field1ZR, field1Rphi, field1phiphi, field1Zphi, field1RZ, field1phiZ, field1ZZ  = \
            eval_field(field, R, phi, Z, coord='tensor', sim=sim[0], time=time[0],quiet=quiet)
        if row == 1:
            field1R = field1RR
            field1phi = field1phiR
            field1Z = field1ZR
        elif row == 2:
            field1R = field1Rphi
            field1phi = field1phiphi
            field1Z = field1Zphi
        elif row == 3:
            field1R = field1RZ
            field1phi = field1phiZ
            field1Z = field1ZZ
        if diff or linear:
            field2RR, field2phiR, field2ZR, field2Rphi, field2phiphi, field2Zphi, field2RZ, field2phiZ, field2ZZ  = \
                eval_field(field, R, phi, Z, coord='tensor', sim=sim[1], time=time[1],quiet=quiet)
            if row == 1:
                field2R = field2RR
                field2phi = field2phiR
                field2Z = field2ZR
            elif row == 2:
                field2R = field2Rphi
                field2phi = field2phiphi
                field2Z = field2Zphi
            elif row == 3:
                field2R = field2RZ
                field2phi = field2phiZ
                field2Z = field2ZZ

    # Evaluate poloidal component
    if coord == 'poloidal':
        #theta  = np.arctan2(Z-Z_mag,R-R_mag)
        #field1 = -np.sin(theta)*field1R + np.cos(theta)*field1Z
        field1 = np.sqrt(field1R**2 + field1Z**2)
        
        # Evaluate second field and calculate difference between two if linear or diff is True
        if diff or linear:
            field2 = -np.sin(theta)*field2R + np.cos(theta)*field2Z


    # Evaluate radial component
    if coord == 'radial':
        theta  = np.arctan2(Z-Z_mag,R-R_mag)
        field1 = np.cos(theta)*field1R + np.sin(theta)*field1Z
        
        # Evaluate second field and calculate difference between two if linear or diff is True
        if diff or linear:
            field2 = np.cos(theta)*field2R + np.sin(theta)*field2Z


    # Calculate difference between two if linear or diff is True
    if diff or linear:
        if coord in ['poloidal', 'radial']:
            field1 = field1 - field2
        elif coord in ['vector', 'tensor']:
            field1R   = field1R - field2R
            field1phi = field1phi - field2phi
            field1Z   = field1Z - field2Z

    if cutr is not None:
        avg_ax = 1
    elif cutz is not None:
        avg_ax = 2
    
    # Calculate average over phi of the field. For a single slice the average is equal to itself
    if coord in ['vector', 'tensor']:
        field1R_ave   = np.average(field1R, avg_ax)
        field1phi_ave = np.average(field1phi, avg_ax)
        field1Z_ave   = np.average(field1Z, avg_ax)
        field1_ave    = [field1R_ave,field1phi_ave,field1Z_ave]
    else:
        field1_ave    = [np.average(field1, avg_ax)]
    R_ave = np.average(R, avg_ax)
    Z_ave = np.average(Z, avg_ax)
    phi_ave = np.average(phi, avg_ax)
    if phys:
        rst_ave    = np.average(rst, avg_ax)
        zst_ave    = np.average(zst, avg_ax)
        R_ave = np.where(np.isnan(rst_ave), R_ave, rst_ave)
        Z_ave = np.where(np.isnan(zst_ave), Z_ave, zst_ave)
    
    #All fields listed in sim.available_fields are returned in MKS units by fusion-io. Other (scalar) fields are returned in m3dc1 units.
    if (field in sim[0].available_fields and units.lower()=='m3dc1') or (field not in sim[0].available_fields and units.lower()=='mks'):
        field1_ave = fpyl.get_conv_field(units,field,field1_ave,sim=sim[0])
    
    
    return sim, time, mesh_ob, R, phi, Z, R_mesh, Z_mesh, R_ave, Z_ave, phi_ave, np.asarray(field1_ave)
