#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Andreas Kleiner
"""
import fpy
import numpy as np
import matplotlib.pyplot as plt
import m3dc1.fpylib as fpyl
from m3dc1.flux_average import flux_average
from m3dc1.read_h5 import readParameter
from m3dc1.unit_conv import unit_conv


def omegastari(sim=None,filename='C1.h5',time=None,units='mks',points=400,n=1,pion=False,fcoords='pest',
               makeplot=False):
    """
    Calculate ion diamagnetic frequency as a function of psi_n
    
    Arguments:

    **sim**
    fpy simulation object.

    **filename**
    Name or path to C1.h5 file to read

    **time**
    The time-slice that will be used if no sim object is given

    **units**
    System of units, can be either 'm3dc1' or 'mks'

    **points**
    Number of radial grid points

    **n**
    Toroidal mode number

    **pion**
    If True, use ion pressure from peqdsk file.

    **fcoords**
    Name of flux coordinate system : 'pest', 'boozer', 'hamada', canonical, ''
    If not specified (i.e. ''), the geometric poloidal angle will be used.

    **makeplot**
    If True, show plot of omegastari.

    """
    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename=filename)
    psin,psi = flux_average('psi',coord='scalar',sim=sim, fcoords=fcoords, linear=False, deriv=0, points=points, phit=0.0, filename=filename, time=time, psin_range=None, units='mks')
    psi = psi*2*np.pi # Because psi is the poloidal flux per radiant as in B = grad(psi) x grad(phi) + B_phi
    ni = flux_average('ni',coord='scalar',sim=sim, fcoords=fcoords, linear=False, deriv=0, points=points, phit=0.0, filename=filename, time=time, psin_range=None, units='mks')[1]
    if pion:
        pi = flux_average('pi',coord='scalar',sim=sim, fcoords=fcoords, linear=False, deriv=0, points=points, phit=0.0, filename=filename, time=time, psin_range=None, units='mks')[1]
    else:
        pi = 0.5*flux_average('p',coord='scalar',sim=sim, fcoords=fcoords, linear=False, deriv=0, points=points, phit=0.0, filename=filename, time=time, psin_range=None, units='mks')[1]
    
    #print(psi)
    #print(ni)
    #print(pi)
    dpidpsi = fpyl.deriv(pi,psi)
    
    Zeff = readParameter('z_ion',sim=sim,listc=False)
    ei = Zeff*1.602176634E-19
    
    omegasi = (n / (ei * ni))*dpidpsi
    if units=='m3dc1':
        omegasi = unit_conv(omegasi,arr_dim='mks',sim=sim,time=-1)
    if makeplot:
        plt.figure()
        plt.plot(psin,omegasi,lw=2)
        ax = plt.gca()
        ax.grid(True,zorder=10,alpha=0.5) #There seems to be a bug in matplotlib that ignores the zorder of the grid #Uncomment for CLT paper
        plt.xlabel(r'$\psi_N$',fontsize=12)
        plt.ylabel(r'$\omega_{*i}$',fontsize=12)
        plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    
    return psin, omegasi
