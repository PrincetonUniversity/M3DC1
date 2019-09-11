# Routine which returns the field at locations given by arrays
#
# Coded on August 16th 2019 by:
# Andreas Kleiner:    akleiner@pppl.gov
# Ralf Mackenbach:    rmackenb@pppl.gov
# Chris Smiet    :    csmiet@pppl.gov

import fpy
import numpy as np



def eval_field(field_name, R, phi, Z, coord='scalar', file_name='C1.h5', time=0):
    """
    Evaluates the field at the locations specified by the 
    R, Z, phi arrays. The output will be array/arrays of the same size.
    
    Arguments:

    **field_name**
    The field that is to be evaluated, i.e. 'B' or 'j', etc..

    **coord**
    The chosen part of a field to be plotted, options are:
    'phi', 'R', 'Z', 'scalar', 'vector'. 'vector' will return three arrays.

    **file_name**
    File name which will be read, i.e. "../C1.h5"

    **time**
    The time-slice which will be used for the field plot

    **elements**
    Input this if you want to accelarate evaluation by
    supplying the elements. The code will check if it is inside
    the convex hull.
    """
    
    # First, let's get the field and mesh from the simulation output
    field = fpy.sim_data(file_name).get_field(field_name,time)
    
    # We check if the field is a scalar or vector field
    check_coord = (R.flatten()[0], phi.flatten()[0], Z.flatten()[0]) 
    length = len(field.evaluate(check_coord))
    if length == 1 and coord != 'scalar':
        raise Exception('You\'re trying to evaluate a component of a scalar field! Please set coord to \'scalar\' if you want to evaluate a scalar field')
    if length == 3 and coord == 'scalar':
        raise Exception('You\'re trying to evaluate a vector field as a scalar! Please set coord to \'R\', \'phi\', \'Z\', or \'vector\'.')
        
    # Set coord to corresponding index value
    if (coord == 'R' or coord == 'scalar'):
        field_idx = 0
    elif coord == 'phi':
        field_idx = 1
    elif coord == 'Z':
        field_idx = 2
    elif coord == 'vector':
        field_idx = -1
    else:
        raise Exception('Please enter valid coordinate. Accepted: \'R\', \'phi\', \'Z\', \'scalar\', \'vector\'.')

    
    # We create the output array, and evaluate the field
    if coord != 'vector':
        field_array = np.zeros_like(R)
        for (idx, r) in np.ndenumerate(R):
            field_array[idx] = field.evaluate((R[idx],phi[idx],Z[idx]))[field_idx]
        
        return field_array

    if coord == 'vector':
        field_array_R   = np.zeros_like(R)
        field_array_phi = np.zeros_like(phi)
        field_array_Z   = np.zeros_like(Z)
        for (idx, r) in np.ndenumerate(R):
            field_tuple = field.evaluate((R[idx],phi[idx],Z[idx]))
            field_array_R[idx]   = field_tuple[0]
            field_array_phi[idx] = field_tuple[1]
            field_array_Z[idx]   = field_tuple[2]
            
        return [field_array_R, field_array_phi, field_array_Z]
