#!/usr/bin/env python3
# This routine returns the field at locations given by arrays
#
# Coded on August 16th 2019 by:
# Andreas Kleiner:    akleiner@pppl.gov
# Ralf Mackenbach:    rmackenb@pppl.gov
# Chris Smiet    :    csmiet@pppl.gov

import fpy
import matplotlib.pyplot as plt
from matplotlib import path
import numpy as np
import m3dc1.fpylib as fpyl
from m3dc1.get_time_of_slice import get_time_of_slice
from m3dc1.unit_conv  import unit_conv
from m3dc1.read_h5 import readParameter


def eval_field(field_name, R, phi, Z, coord='scalar', sim=None, filename='C1.h5', time=None,quiet=False):
    """
    Evaluates the field at the locations specified by the 
    R, Z, phi arrays. The output will be array/arrays of the same size.
    
    Arguments:

    **field_name**
    The field that is to be evaluated, i.e. 'B' or 'j', etc..

    **R, phi, Z**
    Scalars or arrays representing points where the field will be evaluated;
    e.g. R=1.2, phi=0.0, Z=0.0
    or
    R, phi, Z = np.meshgrid(R_linspace, phi_linspace, Z_linspace)

    **coord**
    The chosen part of a field to be plotted, options are:
    'phi', 'R', 'Z', 'scalar', 'vector'. 'vector' will return three arrays.

    **filename**
    File name which will be read, i.e. "../C1.h5"

    **time**
    The time-slice which will be used for the field plot
    """
    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename,time=time)
    if not quiet:
        print("Evaluating field '"+ field_name +"'... ", end=' ', flush=True)
    #Convert coordinates to numpy array if they are not provided as arrays
    R = np.asarray(R,dtype=np.float64)
    phi = np.asarray(phi,dtype=np.float64)
    Z = np.asarray(Z,dtype=np.float64)
    
    if field_name in sim.typedict and sim.typedict[field_name][3]=='simple':#Fields that are directly read with fusion-io, below are composite fields
        field_array = eval_m3dc1_field(field_name, R, phi, Z, coord, sim, filename, time)
    elif field_name=='psi':
        Aphi = eval_m3dc1_field('A', R=R, phi=phi, Z=Z, coord='phi', sim=sim, filename=filename, time=time)
        field_array = np.asarray(R)*Aphi
    elif field_name=='psin':
        psi = eval_field('psi', R=R, phi=phi, Z=Z, coord='scalar', sim=sim, filename=filename, time=time)
        ts_time_ind = fpyl.get_ind_near_val(sim.get_time_trace('psi_lcfs').values, get_time_of_slice(time,sim=sim,units='m3dc1'))
        print(sim.get_time_trace('psi_lcfs').values)
        print(get_time_of_slice(time,sim=sim,units='m3dc1'))
        print(ts_time_ind)
        psi_lcfs = unit_conv(sim.get_time_trace('psi_lcfs').values[ts_time_ind],arr_dim='m3dc1',sim=sim,magnetic_flux=1)
        psi_0 = unit_conv(sim.get_time_trace('psimin').values[ts_time_ind],arr_dim='m3dc1',sim=sim,magnetic_flux=1)
        field_array = (psi - psi_0) / (psi_lcfs - psi_0)
    elif field_name in ['f','I']:
        Bphi = eval_m3dc1_field('B', R=R, phi=phi, Z=Z, coord='phi', sim=sim, filename=filename, time=time)
        field_array = np.asarray(R)*Bphi
    elif field_name=='fp':
        dBphidphi = eval_m3dc1_field_deriv('B', R=R, phi=phi, Z=Z, sim=sim, filename=filename, time=time)[4]
        field_array = np.asarray(R)*dBphidphi
    elif field_name=='|B|':
        B = eval_m3dc1_field('B', R=R, phi=phi, Z=Z, coord='vector', sim=sim, filename=filename, time=time)
        #print(B.shape)
        #print(B)
        #print(B[0])
        #r = np.sqrt(R**2+Z**2)
        normB = np.sqrt(B[0]**2+B[1]**2+B[2]**2)
        field_array = normB
    elif field_name=='B2':
        B = eval_m3dc1_field('B', R=R, phi=phi, Z=Z, coord='vector', sim=sim, filename=filename, time=time)
        #r = np.sqrt(R**2+Z**2)
        B2 = B[0]**2+B[1]**2+B[2]**2
        field_array = B2
        
    elif field_name=='1/B2':
        B = eval_m3dc1_field('B', R=R, phi=phi, Z=Z, coord='vector', sim=sim, filename=filename, time=time)
        #r = np.sqrt(R**2+Z**2)
        B2 = B[0]**2+B[1]**2+B[2]**2
        field_array = 1.0/B2
    elif field_name=='va':
        if not isinstance(sim,fpy.sim_data):
            sim = fpy.sim_data(filename=filename)
        h5file = sim._all_attrs
        mu0 = 4.0E-7*np.pi
        mi = readParameter('ion_mass',h5file=h5file)
        mp = 1.67E-27
        
        B = eval_field('|B|', R=R, phi=phi, Z=Z, coord='scalar', sim=sim, filename=filename, time=time)
        den = eval_m3dc1_field('ne', R=R, phi=phi, Z=Z, coord='scalar', sim=sim, filename=filename, time=time)
        
        va = B/np.sqrt(mu0*den*mi*mp)
        field_array = va
    elif field_name in ['S','lundquist']:
        eta = eval_m3dc1_field('eta', R=R, phi=phi, Z=Z, coord='scalar', sim=sim, filename=filename, time=time)#m3dc1 units
        eta = fpyl.get_conv_field('mks','eta',eta,sim=sim)
        va = eval_field('va', R=R, phi=phi, Z=Z, coord='scalar', sim=sim, filename=filename, time=time)#mks units
        mu0 = 4.0E-7*np.pi
        shape = get_shape(sim=sim,quiet=True)
        a,R0 = (shape["a"],shape["R0"])
        L=a
        field_array = mu0 * L * va / eta
    else:
        field_array = eval_m3dc1_field(field_name, R=R, phi=phi, Z=Z, coord='scalar', sim=sim, filename=filename, time=time)
        # fpyl.printerr('ERROR: Field not supported!')
    if not quiet:
        print('[DONE]')
    return field_array



def eval_m3dc1_field(field_name, R, phi, Z, coord='scalar', sim=None, filename='C1.h5', time=None):
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
    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename)
    field = sim.get_field(field_name,time)
    
    
    # We check if the field is a scalar or vector field
    check_coord = (R.flatten()[0], phi.flatten()[0], Z.flatten()[0]) 
    length = len(field.evaluate(check_coord))
    if length == 1 and coord != 'scalar':
        #raise Exception('You\'re trying to evaluate a component of a scalar field! Please set coord to \'scalar\' if you want to evaluate a scalar field')
        coord = 'scalar'
        print(field_name + " is a scalar field. Setting coord='scalar'...")
    if length == 3 and coord == 'scalar':
        coord == 'vector'
        fpyl.printwarn("WARNING: You are trying to evaluate a vector field as a scalar! Setting coord=\'vector\'. In the future please choose coord=\'R\', \'phi\', \'Z\', or \'vector\'. ")
        #raise Exception('You are trying to evaluate a vector field as a scalar! Please set coord to \'R\', \'phi\', \'Z\', or \'vector\'.')
        
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
        return np.asarray([field_array_R, field_array_phi, field_array_Z])
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
        return np.asarray([field_array_RR, field_array_phiR, field_array_ZR, field_array_Rphi, field_array_phiphi, field_array_Zphi, field_array_RZ, field_array_phiZ, field_array_ZZ])
    else:
        field_array = np.zeros_like(R)
        for (idx, r) in np.ndenumerate(R):
            field_array[idx] = field.evaluate((R[idx],phi[idx],Z[idx]))[field_idx]
        return field_array



def eval_m3dc1_field_deriv(field_name, R, phi, Z, sim=None, filename='C1.h5', time=None):
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
    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename)
    field = sim.get_field(field_name,time)
    
    
    # We check if the field is a scalar or vector field
    check_coord = (R.flatten()[0], phi.flatten()[0], Z.flatten()[0]) 
    length = len(field.evaluate(check_coord))
    
    # We create the output array, and evaluate the field
    if length == 1:#scalar field
        field_array_R   = np.zeros_like(R)
        field_array_phi = np.zeros_like(phi)
        field_array_Z   = np.zeros_like(Z)
        for (idx, r) in np.ndenumerate(R):
            field_tuple = field.evaluate_deriv((R[idx],phi[idx],Z[idx]))
            #print(len(field_tuple))
            field_array_R[idx]   = field_tuple[0]
            field_array_phi[idx] = field_tuple[1]
            field_array_Z[idx]   = field_tuple[2]
        return np.asarray([field_array_R, field_array_phi, field_array_Z])
    elif length == 3:#vector field
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
            field_tuple = field.evaluate_deriv((R[idx],phi[idx],Z[idx]))
            field_array_RR[idx]   = field_tuple[0]
            field_array_phiR[idx] = field_tuple[1]
            field_array_ZR[idx]   = field_tuple[2]
            field_array_Rphi[idx]   = field_tuple[3]
            field_array_phiphi[idx] = field_tuple[4]
            field_array_Zphi[idx]   = field_tuple[5]
            field_array_RZ[idx]   = field_tuple[6]
            field_array_phiZ[idx] = field_tuple[7]
            field_array_ZZ[idx]   = field_tuple[8]
        return np.asarray([field_array_RR, field_array_phiR, field_array_ZR, field_array_Rphi, field_array_phiphi, field_array_Zphi, field_array_RZ, field_array_phiZ, field_array_ZZ])



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
    #paths = cont.collections[0].get_paths()
    paths = []
    for artist in cont.axes.artists + cont.axes.collections:
        if hasattr(artist, 'get_paths'):
            paths.extend(artist.get_paths())
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
    
    Z_max = np.amax(vert[:,1])
    Z_max_ind = fpyl.get_ind_at_val(vert[:,1], Z_max, unique=True)
    Z_min = np.amin(vert[:,1])
    Z_min_ind = fpyl.get_ind_at_val(vert[:,1], Z_min, unique=True)
    
    z0 = (Z_max + Z_min)/2.0
    b = (Z_max - Z_min)/2.0
    
    R_upper = vert[:,0][Z_max_ind]
    R_lower = vert[:,0][Z_min_ind]
    
    kappa = (Z_max-Z_min)/(2*a)
    delta_u = (R0-R_upper)/a
    delta_l = (R0-R_lower)/a
    delta = (delta_u + delta_l) / 2
    
    shape = {"a":a, "R0":R0, "A":R0/a, "z0":z0, "b":b,"kappa":kappa, "delta_u":delta_u, "delta_l":delta_l, "delta":delta}
    
    if not quiet:
        print('R0 = '+str(R0))
        print('a = '+str(a))
        print('A=R0/a = '+str(R0/a))
        print('z0 = '+str(z0))
        print('b = '+str(b))
        print('kappa     =  '+str(kappa))
        print('delta     =  '+str(delta))
        print('delta_u   =  '+str(delta_u))
        print('delta_l   =  '+str(delta_l))
    return shape
