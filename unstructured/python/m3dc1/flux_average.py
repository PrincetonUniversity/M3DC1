#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on February 27 2020

@author: Andreas Kleiner
"""
import numpy as np
import matplotlib.pyplot as plt
import fpy
import m3dc1.fpylib as fpyl
from m3dc1.eval_field import eval_field
from m3dc1.flux_coordinates import flux_coordinates
from m3dc1.unit_conv import unit_conv

#ToDo: Implement unit conversion
#ToDo: Allow for linear option, i.e. flux averaging of a difference between two time slides
def flux_average(field,coord='scalar',sim=None, fcoords=None, linear=False, deriv=0, points=200, phit=0.0, filename='C1.h5', time=0, psin_range=None, device='nstx', units='m3dc1'):
    """
    Calculates the flux average of a quantity
    
    Arguments:

    **field**
    Name of the field to flux average

    **coord**
    For vector fields, component of field to flux average, e.g. R, phi, Z

    **sim**
    fpy simulation object

    **fcoords**
    Name of desired flux coordinate system : 'pest', 'boozer', 'hamada', canonical, ''

    **deriv**
    If 1, calculate and return derivative of flux-averaged quantity dy/dpsin; if 2, calculate derivate w.r.t. to psi

    **phit**
    Toroidal angle where flux average is calculated

    **filename**
    File name which will be read, i.e. "../C1.h5".

    **time**
    The time-slice which will be used for the flux average

    **psin_range**
    Range of normalized flux where a flux surface average will be performed. If None, the flux
    average will be done from the magnetic axis to the last closed flux surface.

    **units**
    Units in which the result will be calculated
    """
    
    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename=filename,time=time)
    else:
        if fcoords is None:
            fcoords = sim.fc.fcoords
    # Calculate flux coodinates if it was not calculated yet or a different flux coordinate system than sim.fc.fcoords is desired
    if isinstance(sim.fc,fpy.flux_coordinates)==False or (fcoords!=None and (sim.fc.fcoords!=fcoords)) or (sim.fc.points!=points):
        if fcoords is None:
            fcoords = ''
            print("FCOORDS SET TO NOTHING")
        sim = flux_coordinates(sim=sim, fcoords=fcoords, filename=filename, time=time, points=points, phit=phit,psin_range=psin_range)
        
    
    
    fc = sim.fc
    nflux = np.asarray(fc.psi_norm)
    
    if field.lower() in ['q','safety factor']:
        fa = np.abs(fc.q)
    elif field=='rho':
        fa = fc.rho
    elif field=='flux_t':
        fa = fc.flux_tor
    elif field=='flux_p':
        fa = fc.flux_pol
    elif field=='psi':
        fa = fc.psi
    elif field in ['current','I']:
        fa = fc.current
        # ToDo: Check if this is correct
        if units=='mks':
            fa = unit_conv(fc.current,arr_dim='m3dc1',filename=filename,current=1)
    elif field =='jav':
        fa = flux_average('current', coord='scalar', sim=sim, fcoords=fcoords, points=points, units=units)[1]/flux_average('polarea', coord='scalar', sim=sim, fcoords=fcoords, points=points, units=units)[1]
    elif field=='jelite':
        mu0 = 4.0E-7*np.pi
        #R0 = sim.fc.r0
        #R0 is taken as the center of the vacuum vessel, as described in Tom Osborne's notes
        deviceR0 = {'nstx': 0.85, 'diiid': 1.6955}
        R0=deviceR0[device.lower()]
        s = np.sign(fc.current[-1])
        #print(s)
        #print(R0)
        if not np.all(np.sign(fc.current)==s):
            fpyl.printerr('ERROR: Current changes sign!')
            return
        #ToDo: it would be cleaner to do the flux averaging of the result, rather than multiplying so many flux averages?
        psi_pp,pprime = flux_average('p', coord='scalar', deriv=2, sim=sim, fcoords=fcoords, points=points, units='mks')
        psi_f,f = flux_average('f', coord='scalar', sim=sim, fcoords=fcoords, points=points, units='mks')
        psi_fp,fprime = flux_average('f', coord='scalar', deriv=2, sim=sim, fcoords=fcoords, points=points, units='mks')
        psi_B2,B2_inv_avg = flux_average('1/B2', coord='scalar', deriv=0, sim=sim, fcoords=fcoords, points=points, units='mks')
        
        #plt.figure()
        #plt.plot(psi_pp,pprime,lw=1)
        #plt.title('pprime')
        
        #plt.figure()
        #plt.plot(psi_f,f,lw=1)
        #plt.title('f')
        
        #plt.figure()
        #plt.plot(psi_fp,fprime,lw=1)
        #plt.title('fprime')
        
        #plt.figure()
        #plt.plot(psi_fp,-f*fprime/R0,lw=1)
        #plt.title('f*fprime')
        
        #plt.figure()
        #plt.plot(psi_B2,B2_inv_avg,lw=1)
        #plt.title('1/B2')
        
        jelite = s * ( R0 * pprime * (f/R0)**2 * B2_inv_avg + f*fprime/(mu0*R0) )
        #psi_jphi,jphi = flux_average('j', coord='phi', sim=sim, fcoords=fcoords, points=points, units='mks')
        #psi_jav,jav = flux_average('jav', coord='scalar', sim=sim, fcoords=fcoords, points=points, units='mks')
        #print(jav[-1])
        #plt.figure()
        #plt.plot(psi_pp,jelite,lw=1)
        #plt.plot(psi_pp,jelite/(2*jav[-1]),lw=1)
        #plt.plot(psi_jphi,jphi,lw=1)
        #plt.title('jelite')
        
        #jelite1 = s * ( R0 * pprime * (f/R0)**2 * B2_inv_avg)
        #jelite2 = s * ( f*fprime/(mu0*R0) )
        #plt.figure()
        #plt.plot(psi_pp,jelite1,lw=1)
        #plt.plot(psi_pp,jelite2,lw=1)
        #plt.title('jelite terms')
        
        fa = jelite
    elif field in ['fs-area']:
        fa = fc.area
    elif field in ['polarea']:
        fa = fc.polarea
    elif field in ['V','volume']:
        fa = fc.V
    elif field=='alpha':
        torphi = np.zeros_like(fc.rpath)
        p = eval_field('p', fc.rpath, torphi, fc.zpath, coord='scalar', sim=sim)
        pavg = flux_average_field(p,fc.j,fc.n,units,'p',sim)
        
        pp = fpyl.deriv(pavg,fc.psi)
        dV = fc.dV_dchi / fc.dpsi_dchi
        alpha = -dV/(2.0*(np.pi)**2) * np.sqrt(fc.V/(2.0*(np.pi)**2*fc.r0)) * pp
        fa = alpha
        #fig = plt.figure()
        #plt.plot(nflux,alpha)
        #cont = plt.contourf(fc.rpath,fc.zpath,p,50)
        #ax = plt.gca()
        #fig.colorbar(cont,ax=ax)
        #plt.axis('equal')
    elif field=='shear':
        q = np.abs(fc.q)
        dqdV = fpyl.deriv(q, fc.V)
        fa = 2.0*fc.V*dqdV/q
    elif field=='elongation':
        fa = None
    elif field=='dqdrho':
        q = np.abs(fc.q)
        fa = fpyl.deriv(q, fc.rho)
    elif field=='f':
        torphi = np.zeros_like(fc.rpath)
        field_val = fc.rpath*eval_field('B', fc.rpath, torphi, fc.zpath, coord='phi', sim=sim)
        fa = flux_average_field(field_val,fc.j,fc.n,units,field,sim)
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
        ffa = flux_average_field(field_val,fc.j,fc.n,units,field,sim)
        fa = ffa
    #fig = plt.figure()
    #plt.plot(nflux,fa)
    #plt.grid(True)
    
    if deriv==1:
        fa = fpyl.deriv(fa, nflux)
    elif deriv==2:
        fa = fpyl.deriv(fa, fc.psi)
    elif deriv==3:
        fa = fpyl.deriv(fa, fc.flux_pol)
    
    return nflux, np.asarray(fa)



def flux_average_field(field,jac,n,units,fieldname,sim):
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
    
    **sim**
    fpy simulation object
    """
    fa = np.zeros(n)
    if units.lower()=='m3dc1':
        field = fpyl.get_conv_field(units,fieldname,field,sim=sim)
    for i in range(n):
        fa[i] = np.sum(field[:,i]*jac[:,i])/np.sum(jac[:,i])
    return fa
