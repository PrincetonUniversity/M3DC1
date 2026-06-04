#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on July 07 2021

@author: Andreas Kleiner
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import fpy
import m3dc1.fpylib as fpyl
from m3dc1.read_h5 import readParameter
from m3dc1.flux_average import flux_average
from m3dc1.eval_field import get_shape
#from m3dc1.plot_flux_average    import plot_flux_average

def get_collisionality(ped_top,q95,sim=None,filename='C1.h5',time=0,fcoords=None,points=200):
    """
    Calculate pedestal collisionality defined in https://doi.org/10.1088/0741-3335/48/5A/S16
    The n_e and T_e profiles are evaluated at the pedestal top

    Arguments:

    **ped_top**
    Location of pedestal top in normalized poloidal flux coordinates.

    **q95**
    Value of q95.

    **sim**
    simulation sim_data object. Can also be list of such objects. If None is provided, plot_shape will read a file and create
    an object.

    **filename**
    File name which will be read, i.e. "../C1.h5"
    Can also be a list of multiple file paths

    **time**
    The timeslice which will be used for the shape plot.
    Should be a list, if filename is a list.

    **fcoords**
    Name of desired flux coordinate system : 'pest', 'boozer', 'hamada', canonical, ''
    If not specified (i.e. ''), the geometric poloidal angle will be used.

    **points**
    Number of points in radial and poloidal direction, i.e.
    number of flux surfaces and poloidal points
    """
    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename=filename,time=time)
    else:
        if fcoords is None and sim.fc is not None:
            fcoords = sim.fc.fcoords
    
    h5file = sim._all_attrs
    version = readParameter('version',h5file=h5file)
    if version >= 23:
        Zeff = readParameter('z_ion',h5file=h5file)
    else:
        Zeff = readParameter('zeff',h5file=h5file)
    
    #Calculate minor radius
    shape = get_shape(sim=sim,quiet=True)
    a,R0 = (shape["a"],shape["R0"])
    
    epsilon = a/R0
    print('Minor radius: a='+str(a))
    print('Major radius: R0='+str(R0))
    
    psi_ne, ne = flux_average('ne', coord='scalar', sim=sim, fcoords=fcoords, points=points, units='mks')
    psi_Te, Te = flux_average('te', coord='scalar', sim=sim, fcoords=fcoords, points=points, units='mks')
    ne_interp = interp1d(psi_ne, ne,kind='cubic',fill_value="extrapolate")
    Te_interp = interp1d(psi_Te, Te,kind='cubic',fill_value="extrapolate")
    
    ne_ped = ne_interp(ped_top)
    Te_ped = Te_interp(ped_top)
    
    coulomb = 31.3 - np.log(np.sqrt(ne_ped)/Te_ped)
    nueff_ped = 6.921E-18 * (q95*R0*ne_ped*Zeff*coulomb)/(Te_ped**2*epsilon**(3/2))
    
    return nueff_ped



def get_pedestal_structure(profile,fitrange=[0.9,1.0],guess=None, sim=None,filename='C1.h5',time=0,fcoords=None,points=200,units='mks'):
    """
    Calculate pedestal height and width for a given profile (p,Te,ne,etc.)
    by fitting a modified tanh.
    
    Arguments:

    **profile**
    Profile for which to calculate the pedestal structure, e.g. 'p', 'ne', etc.

    **fitrange**
    psin range to consider for modified tanh fit

    **guess**
    Initial guess for fot parameters.

    **sim**
    simulation sim_data object. Can also be list of such objects. If None is provided, plot_shape will read a file and create
    an object.

    **filename**
    File name which will be read, i.e. "../C1.h5"
    Can also be a list of multiple file paths

    **time**
    The timeslice which will be used for the shape plot.
    Should be a list, if filename is a list.

    **fcoords**
    Name of desired flux coordinate system : 'pest', 'boozer', 'hamada', canonical, ''
    If not specified (i.e. ''), the geometric poloidal angle will be used.

    **points**
    Number of points in radial and poloidal direction, i.e.
    number of flux surfaces and poloidal points

    **units**
    System of units, can be either 'm3dc1' or 'mks'
    """
    psin, prof = flux_average(profile, coord='scalar', sim=sim, fcoords=fcoords, points=points, units=units)
    
    #plot_flux_average(profile, coord='scalar', sim=sim, fcoords=fcoords, points=points, units=units)
    
    #Initial guess for fit
    print(type(guess))
    print(type(guess) == float)
    if type(guess) is not list:
        if guess is None:
            ped_top_guess = 0.8
        elif type(guess) == float:
            ped_top_guess = guess
        else:
            fpyl.printerr('ERROR: Please check initial guess!')
            return
        
        psin_guess = fpyl.get_ind_at_val(psin,fpyl.find_nearest(psin,ped_top_guess))
        A0 = prof[psin_guess]
        B0 = prof[-1]/4
        s0 = 2000
        w0 = (1.0 - ped_top_guess)/2
        p0 = ped_top_guess + w0
        
        guess = [A0,B0,s0,p0,w0]

    
    i = (np.argwhere((psin>=fitrange[0]) & (psin<=fitrange[1]))).flatten()
    
    xf = psin[i]
    yf = prof[i]
    

    
    #s1 = fpyl.deriv(yf,xf)[int(0.4*len(xf))]
    #print(psin[int(0.4*len(xf))],s1)
    
    print('Initial guess:')
    print(guess)
    #yfit,yerr = curve_fit(mtanh, xf, yf, p0=p0, sigma=w, absolute_sigma=True)
    try:
        yfit,yerr = curve_fit(mtanh, xf, yf, p0=guess)
    except:
        fpyl.printerr('Initial guess [' + ', '.join(guess) + '] did not work')
        return
    #print(yfit,yerr)
    
    ne_interp = interp1d(psin, prof,kind='cubic',fill_value="extrapolate")
    
    
    ped_top = yfit[3]-yfit[4]
    ped_height = ne_interp(ped_top)
    
    y_new = mtanh(psin,yfit[0],yfit[1],yfit[2],yfit[3],yfit[4])
    print('Pedestal slope: '+str(yfit[2]))
    print('Pedestal center: '+str(yfit[3]))
    print('Pedestal top: '+str(ped_top))
    print('Pedestal width: '+str(2*yfit[4]))
    print('Pedestal height (A+B): '+str(yfit[0]+yfit[1]))
    print('Pedestal height (interpol): '+str(ped_height))
    
    
    #print(y_new)
    #print(len(psin),len(y_new))
    plt.figure()
    plt.plot(psin,prof)
    plt.plot(psin,y_new)
    
    return



def mtanh(x, A,B,s,p,w):
    z = (p - x)/w
    mt = ((1.+s*z)*np.exp(z) - np.exp(-z))/(np.exp(z) + np.exp(-z))
    f = A*mt + B
    return f
