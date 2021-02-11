#!/usr/bin/env python3
# This routine returns the field at locations given by arrays
#
# Coded on August 16th 2019 by:
# Andreas Kleiner:    akleiner@pppl.gov
# Ralf Mackenbach:    rmackenb@pppl.gov
# Chris Smiet    :    csmiet@pppl.gov

import fpy
import numpy as np



def eval_field(field_name=None, R=None, phi=None, Z=None, field=None, coord='scalar', sim=None, filename='C1.h5', time=None):
    """
    Evaluates the field at the locations specified by the 
    R, Z, phi arrays. The output will be array/arrays of the same size.
    
    Arguments:

    **field_name**
    The field that is to be evaluated, i.e. 'B' or 'j', etc..

    **coord**
    The chosen part of a field to be plotted, options are:
    'phi', 'R', 'Z', 'scalar', 'vector'. 'vector' will return three arrays.

    **filename**
    File name which will be read, i.e. "../C1.h5"

    **time**
    The time-slice which will be used for the field plot

    **elements**
    Input this if you want to accelarate evaluation by
    supplying the elements. The code will check if it is inside
    the convex hull.
    """
    
    # First, let's get the field and mesh from the simulation output
    if field  is None:
        if sim is None:
            sim = fpy.sim_data(filename)
        field = sim.get_field(field_name,time)
    
    # We check if the field is a scalar or vector field
    check_coord = (R.flatten()[0], phi.flatten()[0], Z.flatten()[0]) 
    length = len(field.evaluate(check_coord))
    if length == 1 and coord != 'scalar':
        raise Exception('You\'re trying to evaluate a component of a scalar field! Please set coord to \'scalar\' if you want to evaluate a scalar field')
    if length == 3 and coord == 'scalar':
        raise Exception('You\'re trying to evaluate a vector field as a scalar! Please set coord to \'R\', \'phi\', \'Z\', or \'vector\'.')
        
    # Set coord to corresponding index value
    if coord in ['R', 'scalar']:
        field_idx = 0
    elif coord == 'phi':
        field_idx = 1
    elif coord == 'Z':
        field_idx = 2
    elif coord in ['vector','tensor']:
        field_idx = -1
    else:
        raise Exception('Please enter valid coordinate. Accepted: \'R\', \'phi\', \'Z\', \'scalar\', \'vector\'.')

    
    # We create the output array, and evaluate the field
    if coord == 'vector':
        field_array_R   = np.zeros_like(R)
        field_array_phi = np.zeros_like(phi)
        field_array_Z   = np.zeros_like(Z)
        for (idx, r) in np.ndenumerate(R):
            field_tuple = field.evaluate((R[idx],phi[idx],Z[idx]))
            #print(len(field_tuple))
            field_array_R[idx]   = field_tuple[0]
            field_array_phi[idx] = field_tuple[1]
            field_array_Z[idx]   = field_tuple[2]
        return [field_array_R, field_array_phi, field_array_Z]
    elif coord == 'tensor':
        field_array_RR     = np.zeros_like(R)
        field_array_phiR   = np.zeros_like(phi)
        field_array_ZR     = np.zeros_like(Z)
        field_array_Rphi   = np.zeros_like(R)
        field_array_phiphi = np.zeros_like(phi)
        field_array_Zphi   = np.zeros_like(Z)
        field_array_RZ     = np.zeros_like(R)
        field_array_phiZ   = np.zeros_like(phi)
        field_array_ZZ     = np.zeros_like(Z)
        for (idx, r) in np.ndenumerate(R):
            field_tuple = field.evaluate((R[idx],phi[idx],Z[idx]))
            field_array_RR[idx]   = field_tuple[0]
            field_array_phiR[idx] = field_tuple[1]
            field_array_ZR[idx]   = field_tuple[2]
            field_array_Rphi[idx]   = field_tuple[3]
            field_array_phiphi[idx] = field_tuple[4]
            field_array_Zphi[idx]   = field_tuple[5]
            field_array_RZ[idx]   = field_tuple[6]
            field_array_phiZ[idx] = field_tuple[7]
            field_array_ZZ[idx]   = field_tuple[8]
        return [field_array_RR, field_array_phiR, field_array_ZR, field_array_Rphi, field_array_phiphi, field_array_Zphi, field_array_RZ, field_array_phiZ, field_array_ZZ]
    else:
        field_array = np.zeros_like(R)
        for (idx, r) in np.ndenumerate(R):
            field_array[idx] = field.evaluate((R[idx],phi[idx],Z[idx]))[field_idx]
        return field_array
