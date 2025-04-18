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
from matplotlib import path
from m3dc1.read_h5 import readParameter
from m3dc1.eval_field import eval_field
from m3dc1.flux_coordinates import flux_coordinates
from m3dc1.unit_conv import unit_conv
from m3dc1.get_shape import get_shape

#ToDo: Implement unit conversion
#ToDo: Allow for linear option, i.e. flux averaging of a difference between two time slides
def flux_average(field,coord='scalar',sim=None, fcoords=None, linear=False, deriv=0, points=200, phit=0.0, filename='C1.h5', time=0, psin_range=None, device=None, units='m3dc1'):
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

    **points**
    Number of flux surfaces between psin = 0 to 1, where flux average is calculated.

    **phit**
    Toroidal angle where flux average is calculated

    **filename**
    File name which will be read, i.e. "../C1.h5".

    **time**
    The time-slice which will be used for the flux average

    **psin_range**
    Range of normalized flux where a flux surface average will be performed. If None, the flux
    average will be done from the magnetic axis to the last closed flux surface.

    **device**
    Device ('nstx', 'diiid', 'sparc', 'mast' or 'mastu') for which flux coordinates are being calculated.
    This determines the major radius, which is taken as machine specific quantity.

    **units**
    Units in which the result will be calculated
    """
    
    #if deriv==True:
    #    deriv=1
    #    fpyl.printnote('Note: deriv=True is not valid. Using deriv=1 (i.e. dy/dpsin) instead.')
    #elif deriv==False:
    #    deriv=0
    #    fpyl.printnote('Note: deriv=False is not valid. Using deriv=0 instead.')
        
    
    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename=filename,time=time)
    else:
        if fcoords is None and sim.fc is not None:
            fcoords = sim.fc.fcoords
    
    # Calculate flux coodinates if it was not calculated yet or a different flux coordinate system than sim.fc.fcoords is desired
    if (not isinstance(sim.fc,fpy.flux_coordinates)) or (fcoords is not None and (sim.fc.fcoords!=fcoords)) or (sim.fc.points!=points):
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
    elif field in ['current','Ip']:
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
        deviceR0 = {'nstx': 0.85, 'diiid': 1.6955, 'sparc': 1.85, 'mast':0.9,'mastu':0.9}
        if device is None:
            fpyl.printerr('ERROR: You need to specify a device to calculate jelite!')
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
    elif field in ['nueff','collisionality']:
        h5file = sim._all_attrs
        version = readParameter('version',h5file=h5file)
        if version >= 23:
            Zeff = readParameter('z_ion',h5file=h5file)
        else:
            Zeff = readParameter('zeff',h5file=h5file)
        
        #R0 = Zeff = readParameter('rzero',h5file=h5file)
        
        #Calculate minor radius
        a,R0 = get_shape(sim)
        
        epsilon = a/R0
        print('Minor radius: a='+str(a))
        print('Major radius: R0='+str(R0))
        #print(a,R0new,R0,epsilon,0.64787234)
        
        q = flux_average('q', coord='scalar', sim=sim, fcoords=fcoords, points=points, units=units)[1]
        ne = flux_average('ne', coord='scalar', sim=sim, fcoords=fcoords, points=points, units='mks')[1]
        Te = flux_average('te', coord='scalar', sim=sim, fcoords=fcoords, points=points, units='mks')[1]
        coulomb = 31.3 - np.log(np.sqrt(ne)/Te)
        fa = 6.921E-18 * (q*R0*ne*Zeff*coulomb)/(Te**2*epsilon**(3/2))
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
    elif field=='eta_spitzer':
        Te = flux_average('te', coord='scalar', sim=sim, fcoords=fcoords, points=points, units='mks')[1]
        fa = Te**(-3/2)
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
    elif field=='ffprime':
        f = flux_average('f', coord='scalar', sim=sim, fcoords=fcoords, points=points, units='mks')[1]
        fprime = flux_average('f', coord='scalar', sim=sim, deriv=True, fcoords=fcoords, points=points, units='mks')[1]
        fa = f*fprime
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



def flux_average_at_psin(psin,field,coord='scalar',sim=None, fcoords=None, deriv=0, points=200, phit=0.0, filename='C1.h5', time=0, device='nstx', units='m3dc1'):
    """
    Returns flux average at specified value of normalized poloidal flux
    
    Arguments:

    **psin**
    Value of psin where the flux average is desired

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

    **points**
    Number of flux surfaces between psin = 0 to 1, where flux average is calculated.

    **phit**
    Toroidal angle where flux average is calculated

    **filename**
    File name which will be read, i.e. "../C1.h5".

    **time**
    The time-slice which will be used for the flux average

    **device**
    Device (e.g. 'nstx', 'diiid') for which flux coordinates are being calculated.
    This determines the major radius, which is taken as machine specific quantity.
    Default is 'nstx'.

    **units**
    Units in which the result will be calculated
    """
    from scipy.interpolate import interp1d
    
    flux, fa = flux_average(field,coord=coord,sim=sim, deriv=deriv, points=points, phit=phit, filename=filename, time=time, device=device, fcoords=fcoords, units=units)
    
    fs_inter = interp1d(flux, fa, kind='cubic')
    
    return fs_inter(psin)
