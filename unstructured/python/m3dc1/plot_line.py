#!/usr/bin/env python3
#
# Coded on August 28th 2019 by:
# Andreas Kleiner:    akleiner@pppl.gov
# Ralf Mackenbach:    rmackenb@pppl.gov
# Chris Smiet    :    csmiet@pppl.gov

import fpy
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.ticker as ticker
from m3dc1.unit_conv  import unit_conv
from m3dc1.eval_field import eval_field
from m3dc1.plot_mesh import plot_mesh
rc('text', usetex=True)




def plot_line(field, coord='scalar', angle=0, file_name='C1.h5', time=0, phi=0, linear=False, diff=False, tor_av=1, units='mks'):
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
    Angle in degrees between the Z=0 line going through the magnetic axis,
    and the line on which the field is evaluated. Anti-clockwise
    are positive angles.

    **file_name**
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
    This could be the difference of two files (file_name=['a/C1.h5','b/C1.h5']),
    or the difference between two time-slices (time=[t1,t2])
    If list for both time and file_name are given file1 will be evaluated at time1,
    and file2 at time2

    **tor_av**
    Calculates the average field over tor_av number of toroidal planes
    """
    angle = np.radians(angle)


    # make file name iterable if it is a string and not a list of strings
    file_name = (file_name,) if not isinstance(file_name, (tuple, list)) else file_name
    
    # make time iterable if it is a single int and not if it is list of ints
    time = (time,) if not isinstance(time, (tuple, list)) else time
    

    ### Input error handling ###
    if diff==True:
        if len(time)!=2:
            raise Exception('Please input two times for differences')
    
    if diff==True and linear==True:
        raise Exception('Please choose diff or linear, not both.')

    if diff==False:
        if (len(file_name)>1 or len(time)>1):
            raise Exception('Multiple file/time slices detected. Please set diff=True or input single slices')

    if (coord == 'R' or coord == 'scalar'):
        field_idx = 0
    elif coord == 'phi':
        field_idx = 1
    elif coord == 'Z':
        field_idx = 2
    elif coord=='poloidal':
        field_idx = None
    elif coord=='radial':
        field_idx = None
    else:
        raise Exception('Please enter valid coordinate. Accepted: \'R\', \'phi\', \'Z\', \'scalar\'')
    ### End input error handling ###
    


    
    # If linear, the difference needs to be taken with time 0
    if linear==True:
        time = [time[0],0]

    
    
    # If either file_name or time is a list, we will convert both of them to lists of length two.
    if (len(file_name)==2 and len(time)==1):
        time = [time, time]
    if (len(file_name)==1 and len(time)==2):
        file_name = [file_name[0],file_name[0]]
    

    
    
    
    # Make 3D grid based on the mesh points
    mesh_ob      = fpy.sim_data(file_name[0]).get_mesh(time=0)
    mesh_pts     = mesh_ob.elements
    R_mag        =  fpy.sim_data(file_name[0]).get_time_traces('xmag').values[0]
    Z_mag        =  fpy.sim_data(file_name[0]).get_time_traces('zmag').values[0]
    R_range      = [np.amin(mesh_pts[:,4]), np.amax(mesh_pts[:,4])]
    Z_range      = [np.amin(mesh_pts[:,5]), np.amax(mesh_pts[:,5])]
    R_straight   = (np.linspace(R_range[0], R_range[1], 10000, endpoint=True)-R_mag) 
    R_linspace   = (np.linspace(R_range[0], R_range[1], 10000, endpoint=True)-R_mag)*np.cos(angle) + R_mag   
    Z_linspace   = (np.linspace(Z_range[0], Z_range[1], 10000, endpoint=True)-Z_mag)*np.sin(angle) + Z_mag
    phi_linspace = np.linspace(phi,      (360+phi), tor_av, endpoint=False)
    R, phi       = np.meshgrid(R_linspace, phi_linspace)
    Z, phi       = np.meshgrid(Z_linspace, phi_linspace)
    L_signed     = np.sign(R_straight)*np.sqrt((R_linspace-R_mag)**2 + (Z_linspace-Z_mag)**2)
    
    



    # Evaluate usual vector components
    if coord != 'poloidal' and coord !='radial':
        # Evaluate field
        print('Evaluating field... ')
        field1 = eval_field(field, R, phi, Z, coord=coord, file_name=file_name[0], time=time[0])
    
        # Evaluate second field and calculate difference between two if linear or diff is True
        if (diff == True or linear==True):
            print('Evaluating second field... ')
            field2 = eval_field(field, R, phi, Z, coord=coord, file_name=file_name[1], time=time[1])
            field1 = field1 - field2



    # Evaluate poloidal component
    if coord == 'poloidal':
        print('Evaluating poloidal field... ')
        field1R, field1phi, field1Z  = eval_field(field, R, phi, Z, coord='vector', file_name=file_name[0], time=time[0])
        theta                        = np.arctan2(Z-Z_mag,R-R_mag)
        field1                       = -np.sin(theta)*field1R + np.cos(theta)*field1Z
        
        # Evaluate second field and calculate difference between two if linear or diff is True
        if (diff == True or linear==True):
            print('Evaluating second poloidal field... ')
            field2R, field2phi, field2Z  = eval_field(field, R, phi, Z, coord='vector', file_name=file_name[1], time=time[1])
            field2                       = -np.sin(theta)*field2R + np.cos(theta)*field2Z
            field1                       = field1 - field2


    # Evaluate radial component
    if coord == 'radial':
        print('Evaluating radial field... ')
        field1R, field1phi, field1Z  = eval_field(field, R, phi, Z, coord='vector', file_name=file_name[0], time=time[0])
        theta                        = np.arctan2(Z-Z_mag,R-R_mag)
        field1                       = np.cos(theta)*field1R + np.sin(theta)*field1Z
        
        # Evaluate second field and calculate difference between two if linear or diff is True
        if (diff == True or linear==True):
            print('Evaluating second poloidal field... ')
            field2R, field2phi, field2Z  = eval_field(field, R, phi, Z, coord='vector', file_name=file_name[1], time=time[1])
            field2                       = np.cos(theta)*field2R + np.sin(theta)*field2Z
            field1                       = field1 - field2



    
    # Calculate average over phi of the field. For a single slice the average is equal to itself
    field1_ave     = np.average(field1, 0)
    R_ave          = np.average(R, 0)
    Z_ave          = np.average(Z, 0)
    
 



    if units=='M3DC1':
        if field == 'j':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',current_density=1)
            label = 'current density (M3DC1 units)'
        if field == 'ni':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',particles=1,length=-3)
            label = 'number density (M3DC1 units)'
        if field == 'ne':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',particles=1,length=-3)
            label = 'number density (M3DC1 units)'
        if field == 'v':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',velocity=1)
            label = 'velocity (M3DC1 units)'
        if field == 'B':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',magnetic_field=1)
            label = 'magnetic field strength (M3DC1 units)'
        if field == 'p':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',pressure=1)
            label = 'pressure (M3DC1 units)'
        if field == 'pi':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',pressure=1)
            label = 'pressure (M3DC1 units)'
        if field == 'pe':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',pressure=1)
            label = 'pressure (M3DC1 units)'
        if field == 'ti':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',temperature=1)
            label = 'temperature (M3DC1 units)'
        if field == 'te':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',temperature=1)
            label = 'temperature (M3DC1 units)'
        if field == 'A':
            field1_ave = unit_conv(field1_ave,arr_dim='mks',magnetic_field=1,length=1)
            label = 'vector potential (M3DC1 units)'
    
    if units=='mks':
        if field == 'j':
            label = '$A/m^2$'
        if field == 'ni':
            label = 'particles/$m^3$'
        if field == 'ne':
            label = 'particles/$m^3$'
        if field == 'v':
            label = '$m/s$'
        if field == 'B':
            label = '$T$'
        if field == 'p':
            label = '$Pa$'
        if field == 'pi':
            label = '$Pa$'
        if field == 'pe':
            label = '$Pa$'
        if field == 'ti':
            label = '$eV$'
        if field == 'te':
            label = '$eV$'
        if field == 'A':
            label = '$Tesla \cdot m$'


    # Plot routines
    plt.plot(L_signed, field1_ave)
    plt.xlabel('Signed distance from magnetic axis [m]')
    plt.ylabel(label)
    
    plt.show()
