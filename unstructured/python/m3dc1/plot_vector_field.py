#!/usr/bin/env python3
#
# Coded on August 16th 2019 by:
# Andreas Kleiner:    akleiner@pppl.gov
# Ralf Mackenbach:    rmackenb@pppl.gov
# Chris Smiet    :    csmiet@pppl.gov

import fpy
import numpy as np
from m3dc1.eval_field import eval_field
from mayavi import mlab



def plot_vector_field(field, file_name='C1.h5', time=0, linear=False, diff=False, R_res=100, phi_res=8, Z_res=100):
    """
    Plots the field of a file
    
    Arguments:

    **field**
    The field that is to be plotted, i.e. 'magnetic field'

    **file_name**
    File name which will be read, i.e. "../C1.h5"
    Can also be a list of two filepaths when used for diff

    **time**
    The time-slice which will be used for the field plot

    **linear**
    Plot the linear part of the field (so the equilibrium is subtracted).
    True/False
    
    **diff**
    Plot the difference of two fields. 
    This could be the difference of two files (file_name=['a/C1.h5','b/C1.h5']),
    or the difference between two time-slices (time=[t1,t2])
    If list for both time and file_name are given file1 will be evaluated at time1,
    and file2 at time2

    **R_res**
    Number of points sampled in R direction

    **phi_res**
    Number of points sampled in phi direction

    **Z_res**
    Number of points sampled in Z direction
    """
    # make file name iterable if it is a string and not a list of strings
    file_name = (file_name,) if not isinstance(file_name, (tuple, list)) else file_name
    
    # make time iterable if it is a single int and not if it is list of ints
    time = (time,) if not isinstance(time, (tuple, list)) else time
    

    # Input error handling
    if diff==True:
        if (len(file_name)==1 and len(time)==1):
            raise Exception('Please input two files and/or timeslices for diff')
        elif(len(file_name)==2 and len(time)==1):
            raise Exception('Please choose the times at which the different files are to be evaluated. Single times are not allowed.')
    
    if diff==True and linear==True:
        raise Exception('Please choose diff or linear, not both.')

    if diff==False:
        if (len(file_name)>1 or len(time)>1):
            raise Exception('Multiple file/time slices detected. Please set diff=True or input single slices')
    
    


    
    # If linear, the difference needs to be taken with time 0
    if linear==True:
        time = [time[0],0]

    
    
    # If either file_name or time is a list, we will convert both of them to lists of length two.
    if (len(file_name)==2 and len(time)==1):
        time = [time, time]
    if (len(file_name)==1 and len(time)==2):
        file_name = [file_name[0],file_name[0]]
    

    
    
    
    # Make 3D grid based on max and min values of mesh points
    mesh_field  = fpy.sim_data(file_name[0]).get_mesh(time=0).elements
    

    # Check simulation bounds
    R_range = [np.amin(mesh_field[:,4]),np.amax(mesh_field[:,4])]
    Z_range = [np.amin(mesh_field[:,5]),np.amax(mesh_field[:,5])]

    # Make grid based on simulation bounds
    R_linspace   = np.linspace(R_range[0], R_range[1], R_res,   endpoint=True)
    phi_linspace = np.linspace(0,          360,        phi_res, endpoint=False)
    Z_linspace   = np.linspace(Z_range[0], Z_range[1], Z_res,   endpoint=True)
    R, phi, Z    = np.meshgrid(R_linspace, phi_linspace, Z_linspace)
    
    # Evaluate field
    print('Evaluating field... ')
    fields1 = eval_field(field, R, phi, Z, coord='vector', file_name=file_name[0], time=time[0])
    
    # Evaluate second field and calculate difference between two if linear or diff is True
    if (diff == True or linear==True):
        print('Evaluating second field... ')
        fields2 = eval_field(field, R, phi, Z, coord='vector', file_name=file_name[1], time=time[1])
        fields1 = [fields1[0]-fields2[0], fields1[1]-fields2[1], fields1[2]-fields2[2]]





    # Plot routines
    mlab.clf()
    mlab.quiver3d(R.flatten(), Z.flatten(), phi.flatten(), fields1[0].flatten(), fields1[2].flatten(), fields1[1].flatten())
    mlab.draw()
    mlab.show()
