#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on February 27 2020

@author: Andreas Kleiner
"""
import math
import numpy as np
import matplotlib.pyplot as plt
import fpy
import m3dc1.fpylib as fpyl
from m3dc1.eval_field import eval_field
from m3dc1.flux_coordinates import flux_coordinates
from m3dc1.unit_conv import unit_conv

#ToDo: Implement unit conversion
#ToDo: Allow for linear option, i.e. flux averaging of a difference between two time slides
def flux_average(field,coord='scalar',sim=None, linear=False, phit=0.0, filename='C1.h5', time=0, fcoords=None, psin_range=None, units='m3dc1'):
    """
    Calculates the flux average of a quantity
    
    Arguments:

    **field**
    Name of the field to flux average

    **coord**
    For vector fields, component of field to flux average, e.g. R, phi, Z

    **coord**
    fpy simulation object
    
    **fcoords**
    Name of desired flux coordinate system : 'pest', 'boozer', 'hamada', canonical, ''
    
    **psin_range**
    Range of normalized flux where a flux surface average will be performed. If None, the flux
    average will be done from the magnetic axis to the last closed flux surface.
    """
    
    if isinstance(sim,fpy.sim_data)==False:
        sim = fpy.sim_data(filename,time=time)
    # Calculate flux coodinates if it was not calculated yet or a different flux coordinate system than sim.fc.fcoords is desired
    if isinstance(sim.fc,fpy.flux_coordinates)==False or (fcoords!=None and (sim.fc.fcoords!=fcoords)):
        if fcoords==None:
            fcoords = ''
        sim = flux_coordinates(sim=sim, fcoords=fcoords, phit=phit,psin_range=psin_range)
    
    
    fc = sim.fc
    flux = fc.psi
    nflux = fc.psi_norm
    
    if field == 'q' or field.lower() == 'safety factor':
        fa = np.abs(fc.q)
    elif field=='rho':
        fa = fc.rho
    elif field=='flux_t':
        fa = fc.flux_tor
    elif field=='flux_p':
        fa = fc.flux_pol
    elif field=='V' or field.lower()=='volume':
        fa = fc.V
    elif field=='alpha':
        torphi = np.zeros_like(fc.rpath)
        p = eval_field('p', fc.rpath, torphi, fc.zpath, coord='scalar', sim=sim)
        pavg = flux_average_field(p,fc.j,fc.n,units,'p')
        
        pp = fpyl.deriv(pavg,fc.psi)
        dV = fc.dV_dchi / fc.dpsi_dchi
        alpha = -dV/(2.0*(math.pi)**2) * np.sqrt(fc.V/(2.0*(math.pi)**2*fc.r0)) * pp
        fa = alpha
        #fig = plt.figure()
        #plt.plot(nflux,alpha)
        #cont = plt.contourf(fc.rpath,fc.zpath,p,50)
        #ax = plt.gca()
        #fig.colorbar(cont,ax=ax)
        #plt.axis('equal')
    elif field=='shear':
        q = np.abs(fc.q)
        dqdV = deriv(q, fc.V)
        fa = 2.0*fc.V*dqdV/q
    elif field=='elongation':
        fa = None
    elif field=='dqdrho':
        q = np.abs(fc.q)
        fa = deriv(q, fc.rho)
    elif field=='lambda':
        fa = None
    elif field=='beta_pol':
        fa = None
    elif field=='alpha2':
        fa = None
    elif field=='kappa_implied':
        fa = None
    elif field=='lambda':
        fa = None
    elif field=='amu_implied':
        fa = None
    else:
        #Evaluate field via fusion io
        torphi = np.zeros_like(fc.rpath)
        field_val = eval_field(field, fc.rpath, torphi, fc.zpath, coord=coord, sim=sim)
        ffa = flux_average_field(field_val,fc.j,fc.n,units,field)
        fa = ffa
    #fig = plt.figure()
    #plt.plot(nflux,fa)
    #plt.grid(True)
    
    return nflux, fa



def flux_average_field(field,jac,n,units,fieldname):
    """
    Calculates and return the flux average of a field
    
    Arguments:

    **field**
    2D array containing field evaluated on points within the flux surfaces 
    
    **jac**
    2D array containing the Jacobian
    
    **n**
    Number of points in poloidal direction
    
    **units**
    Units for flux averaged field
    
    **fieldname**
    String containing field name
    """
    fa = np.zeros(n)
    if units.lower()=='m3dc1':
        field = fpyl.get_conv_field(units,fieldname,field)
    for i in range(n):
        fa[i] = np.sum(field[:,i]*jac[:,i])/np.sum(jac[:,i])
    return fa
