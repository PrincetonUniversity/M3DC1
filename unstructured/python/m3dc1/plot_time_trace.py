#!/usr/bin/env python3

import numpy as np
import pint
import os
import matplotlib.pyplot as plt
from matplotlib import rc
import fpy
from m3dc1.unit_conv import unit_conv
rc('text', usetex=True)

def plot_time_trace(trace,sim=None,units='mks',filename='C1.h5',growth=False,renorm=False,yscale='linear'):
    """
    Plots the time trace of some quantity. All available
    time traces can be found in the M3DC1 documentation.
    
    Arguments:

    **trace**
    String containing the trace to be plotted

    **sim**
    simulation sim_data object. If none is provided,
    a new object will be created.

    **units**
    The units in which the trace will be plotted

    **filename**
    Contains the path to the C1.h5 file.

    **growth**
    Determines wether to calculate the derivative.
    True/False

    **renorm**
    Removes spikes that are caused by renormalizations
    in linear stability calculations. Interpolates at
    the locations of the spike. Should only be used if
    growth=True.
    
    **yscale**
    Scale of y axis, e.g. linear, log
    """
    if sim is None:
        sim = fpy.sim_data(filename)
    constants = sim.get_constants()

    plt.ylabel(trace)
    time = sim.get_time_trace('time').values
    plt.xlabel(r'time $\tau_A$')
    if units=='mks':
        time = unit_conv(time,sim=sim,time=1)
        plt.xlabel(r'time $[s]$')

    if trace=='me':
        y_axis = sim.get_time_trace('E_MP').values + \
                 sim.get_time_trace('E_MT').values
        if units=='mks':
            y_axis = unit_conv(y_axis, arr_dim='M3DC1', sim=sim, energy=1)
            plt.ylabel(r'Magnetic energy $[J]$')


    if trace=='ke':
        y_axis = sim.get_time_trace('E_KP').values + \
                 sim.get_time_trace('E_KT').values + \
                 sim.get_time_trace('E_K3').values
        if growth:
            plt.ylabel(r'$\gamma/\omega_A$')
        if units=='mks':
            y_axis = unit_conv(y_axis, arr_dim='M3DC1', sim=sim, energy=1)
            if growth:
                plt.ylabel(r'$\gamma$ $[s^{-1}]$')
            else:
                plt.ylabel(r'Kinetic energy $[J]$')
        
    if trace=='it':
        y_axis = sim.get_time_trace('toroidal_current').values
        if units=='mks':
            y_axis = unit_conv(y_axis, arr_dim='M3DC1', sim=sim, current=1)
            plt.ylabel(r'Current $[A]$')


    if trace=='ip':
        y_axis = sim.get_time_trace('toroidal_current_p').values
        if units=='mks':
            y_axis = unit_conv(y_axis, arr_dim='M3DC1', sim=sim, current=1)
            plt.ylabel(r'Current $[A]$')

    
    if trace=='vl':
        y_axis = sim.get_time_trace('loop_voltage').values
        if units=='mks':
            y_axis = unit_conv(y_axis, arr_dim='M3DC1', sim=sim, voltage=1)
            plt.ylabel(r'Voltage $[V]$')
            
    
    if trace=='iw':
        y_axis = sim.get_time_trace('toroidal_current_w').values
        if units=='mks':
            y_axis = unit_conv(y_axis, arr_dim='M3DC1', sim=sim, current=1)
            plt.ylabel(r'Current $[A]$')


    if trace=='itot':
        y_axis = sim.get_time_trace('toroidal_current_w').values + \
                 sim.get_time_trace('toroidal_current').values
        if units=='mks':
            y_axis = unit_conv(y_axis, arr_dim='M3DC1', sim=sim, current=1)
            plt.ylabel(r'Current $[A]$')

    ## Formally not right, we would like BT as a function of time
    if trace=='bt':
        Ave_P  = sim.get_time_trace('Ave_P').values
        B0     = constants.B0
        y_axis = 2*Ave_P/B0**2
        
        
    if trace=='bp':
        if constants.numvar < 3:
            raise Exception('Must have numvar=3 for a poloidal beta calculation!')
        gamma = constants.gamma
        E_P = sim.get_time_trace('E_P').values
        tor_cur = sim.get_time_trace('toroidal_current').values
        y_axis  = 2*(gamma-1)*E_P/(tor_cur**2)
        
        
    if trace=='bn':
        raise Exception('Not coded yet')

    
    if trace=='beta':
        gamma = constants.gamma
        E_P   = sim.get_time_trace('E_P').values
        E_MP  = sim.get_time_trace('E_MP').values
        E_MT  = sim.get_time_trace('E_MT').values
        y_axis= (gamma-1)*E_P/(E_MP + E_MT)
        
        
    if trace in ['psibound', 'psilim']:
        y_axis = sim.get_time_trace('psi_lcfs').values
        if units=='mks':
            y_axis = unit_conv(y_axis, arr_dim='M3DC1', sim=sim, magnetic_field=1, length=2)
            plt.ylabel(r'Current $[T][m]^2$')

    
    if trace=='psimin':
        y_axis = sim.get_time_trace('psimin').values
        if units=='mks':
            y_axis = unit_conv(y_axis, arr_dim='M3DC1', sim=sim, magnetic_field=1, length=2)
            plt.ylabel(r'Current $[T][m]^2$')
    

    if trace=='dt':
        y_axis = sim.get_time_trace('dt').values
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,time=1)
            plt.ylabel(r'Time $[s]$')
        
    
    if trace=='volume':
        y_axis = sim.get_time_trace('volume').values
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,length=3)
            plt.ylabel(r'Volume $[m]^3$')


    if trace=='toroidal flux':
        y_axis = sim.get_time_trace('toroidal_flux_p').values
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,magnetic_field=1,length=2)
            plt.ylabel(r'Current $[T][m]^2$')


    if trace=='reconnected flux':
        y_axis = sim.get_time_trace('reconnected_flux').values
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,magnetic_field=1,length=2)
            plt.ylabel(r'Current $[T][m]^2$')
    
    
    if trace=='temax':
        y_axis = sim.get_time_trace('temax').values
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,temperature=1)
            plt.ylabel(r'Temperature $[eV]$')
    

    if trace=='p':
        y_axis = sim.get_time_trace('E_P').values
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,energy=1)
            plt.ylabel(r'Pressure $[J]$')


    if trace=='pe':
        y_axis = sim.get_time_trace('E_PE').values
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,energy=1)
            plt.ylabel(r'Electron pressure $[J]$')

    ## IDL has L=3, but isn't this already integrated over volume?
    if trace=='n':
        y_axis = sim.get_time_trace('particle_number').values
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,particles=1)
            plt.ylabel(r'Number of particles')

    ## IDL has L=3, but isn't this already integrated over volume?
    if trace=='ne':
        y_axis = sim.get_time_trace('electron_number').values
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,particles=1)
            plt.ylabel(r'Number of particles')

    
    if trace=='angular momentum':
        y_axis = sim.get_time_trace('angular_momentum').values
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,energy=1,time=1)
            plt.ylabel(r'Angular momentum $[J][s]$')


    if trace=='vorticity':
        y_axis = sim.get_time_trace('circulation').values
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,length=2,time=-1)
            plt.ylabel(r'Vorticity $[m]^2/[s]$')

            
    if trace=='bwb2':
        amupar = constants.amupar
        y_axis = 3.0/(4.0*amupar)*sim.get_time_trace('parallel_viscous_heating').values
        if units=='mks':
            # Not too sure about the units here...
            y_axis = unit_conv(y_axis, sim=sim,energy=1,time=-2,length=3)
            plt.ylabel(r'Please verify units in code!')


    if trace=='li':
        R0       = constants.R0
        psi_lcfs = sim.get_time_trace('psi_lcfs').values
        psimin   = sim.get_time_trace('psimin').values
        tor_cur  = sim.get_time_trace('toroidal_current_p').values
        y_axis   = -4.0 * np.pi * (psi_lcfs - psimin) / (R0*tor_cur)
        plt.ylabel(r'Normalized inductance')


    if trace=='li3':
        R0       = constants.R0
        W_M      = sim.get_time_trace('W_M').values
        tor_curp = sim.get_time_trace('toroidal_current_p').values
        y_axis   = 4.0 * W_M / (tor_curp**2 * R0)
        plt.ylabel(r'Normalized inductance')


    if trace=='flux':
        psi_lcfs = sim.get_time_trace('psi_lcfs').values
        psimin   = sim.get_time_trace('psimin').values
        y_axis   = -2*np.pi*(psi_lcfs-psimin)
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,length=2,magnetic_field=-1)
            plt.ylabel(r'Normalized inductance')


    if trace=='xmag':
        y_axis = sim.get_time_trace('xmag').values
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,length=1)
            plt.ylabel(r'R coordinate of magnetic axis $[m]$')


    if trace=='zmag':
        y_axis = sim.get_time_trace('zmag').values
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,length=1)
            plt.ylabel(r'Z coordinate of magnetic axis $[m]$')


    if trace=='runaways':
        y_axis = sim.get_time_trace('runaways').values
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,particles=1)
            plt.ylabel(r'Number of runaways')


    if trace=='radiation':
        y_axis = sim.get_time_trace('radiation').values
        # Not too sure about dimensions
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,energy=1,time=-1)
            plt.ylabel(r'Radation losses $[J]/[s]$')

            
    if trace=='POhm':
        e_mpd = sim.get_time_trace('E_MPD').values
        e_mtd = sim.get_time_trace('E_MTD').values
        y_axis= -(e_mpd+e_mtd)
        # Not too sure about dimensions
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,energy=1,time=-1)
            plt.ylabel(r'Ohmic heating $[J]/[s]$')


    if trace=='pelr':
        y_axis = sim.get_time_trace('pellet_rate').values
        # Not too sure about dimensions
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,particles=1,time=-1)
            plt.ylabel(r'Pellet rate $1/[s]$')


    if trace=='pelvar':
        y_axis = sim.get_time_trace('pellet_var').values
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,length=1)
            plt.ylabel(r'Pellet variance $[m]$')


    if trace=='pelrad':
        if constants.version < 26:
            y_axis = sim.get_time_trace('r_p2').values # RP2 is my fav space
        else:
            y_axis = sim.get_time_trace('r_p').values
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,length=1)
            plt.ylabel(r'Pellet radius $[m]$')


    if trace=='pelrpos':
        if constants.version < 26:
            y_axis = sim.get_time_trace('pellet_x').values
        else:
            y_axis = sim.get_time_trace('pellet_r').values
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,length=1)
            plt.ylabel(r'Pellet R-coordinate $[m]$')


    if trace=='pelzpos':
        y_axis = sim.get_time_trace('pellet_z').values
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,length=1)
            plt.ylabel(r'Pellet Z-coordinate $[m]$')


    if trace=='bharmonics':
        y_axis = sim.get_time_trace('bharmonics').values
        time   = sim.get_time_trace('bharmonics').time
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,energy=1)
            time   = unit_conv(time, sim=sim,time=1)
            plt.ylabel(r'Magnetic energy harmonics $[J]$')
            plt.xlabel(r'Time $[s]$')
        for n in range(0,np.shape(y_axis)[1]):
            plt.plot(time,y_axis[:,n],label='n = {}'.format(n))
            plt.legend(loc='best')
        plt.show()
        return

    if trace=='keharmonics':
        y_axis = sim.get_time_trace('keharmonics').values
        time   = sim.get_time_trace('keharmonics').time
        if units=='mks':
            y_axis = unit_conv(y_axis, sim=sim,energy=1)
            time   = unit_conv(time, sim=sim,time=1)
            plt.ylabel(r'Kinetic energy harmonics $[J]$')
            plt.xlabel(r'Time $[s]$')
        for n in range(0,np.shape(y_axis)[1]):
            plt.plot(time,y_axis[:,n],label='n = {}'.format(n))
            plt.legend(loc='best')
        plt.show()
        return

    #if growth==True:
    #    y_axis = np.gradient(y_axis, time)
    if growth:
        y_axis = 1.0/y_axis[1:] * np.diff(y_axis)/np.diff(time)
    
    if renorm:
        for i in range(len(y_axis)-1):
            if(abs(y_axis[i+1]/y_axis[i]) < 1E-9):
                print('Renormalization found at '+str(time[i]))
                print(y_axis[i],y_axis[i-1]+y_axis[i+1])
                y_axis[i] = (y_axis[i-1] + y_axis[i+1])/2.0

    # If one array has new data but the other one doesn't 
    # plot only previous data
    if y_axis.shape != time.shape:
        ymax   = np.amax(y_axis.shape)
        tmax   = np.amax(time.shape)
        maxidx = np.amin([ymax,tmax])
        time   = time[0:maxidx]
        y_axis = y_axis[0:maxidx]
    
    plt.grid(True)
    plt.plot(time,y_axis)
    plt.yscale(yscale)
    plt.ticklabel_format( axis='y', style='sci',useOffset=False)
    plt.show()
