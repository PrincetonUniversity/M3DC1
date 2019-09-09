import matplotlib.pyplot as plt
import pint
import os
import fpy
from unit_conv import unit_conv
import numpy as np
import sys
from matplotlib import rc
rc('text', usetex=True)

def plot_time_trace(trace,units='mks',file_name='C1.h5'):
    """
    Plots the time trace of some quantity. All available
    time traces can be found in the M3DC1 documentation.
    
    Arguments:

    **trace**
    String containing the trace to be plotted

    **units**
    The units in which the trace will be plotted

    **file_name**
    Contains the path to the C1.h5 file.
    """
    
    sim =  fpy.sim_data(file_name)
    plt.ylabel(trace)
    time = sim.get_time_traces('time').values
    plt.xlabel(r'time $\tau_A$')
    if units=='mks':
        time = unit_conv(time,time=1)
        plt.xlabel(r'time $[s]$')

    if trace=='me':
        y_axis = sim.get_time_traces('E_MP').values + \
                 sim.get_time_traces('E_MT').values
        if units=='mks':
            y_axis = unit_conv(y_axis, arr_dim='M3DC1', energy=1)
            plt.ylabel(r'Magnetic energy $[J]$')


    if trace=='ke':
        y_axis = sim.get_time_traces('E_KP').values + \
                 sim.get_time_traces('E_KT').values + \
                 sim.get_time_traces('E_K3').values
        if units=='mks':
            y_axis = unit_conv(y_axis, arr_dim='M3DC1', energy=1)
            plt.ylabel(r'Kinetic energy $[J]$')
        
    if trace=='it':
        y_axis = sim.get_time_traces('toroidal_current').values
        if units=='mks':
            y_axis = unit_conv(y_axis, arr_dim='M3DC1', current=1)
            plt.ylabel(r'Current $[A]$')


    if trace=='ip':
        y_axis = sim.get_time_traces('toroidal_current_p').values
        if units=='mks':
            y_axis = unit_conv(y_axis, arr_dim='M3DC1', current=1)
            plt.ylabel(r'Current $[A]$')

    
    if trace=='vl':
        y_axis = sim.get_time_traces('loop_voltage').values
        if units=='mks':
            y_axis = unit_conv(y_axis, arr_dim='M3DC1', voltage=1)
            plt.ylabel(r'Voltage $[V]$')
            
    
    if trace=='iw':
        y_axis = sim.get_time_traces('toroidal_current_w').values
        if units=='mks':
            y_axis = unit_conv(y_axis, arr_dim='M3DC1', current=1)
            plt.ylabel(r'Current $[A]$')


    if trace=='itot':
        y_axis = sim.get_time_traces('toroidal_current_w').values + \
                 sim.get_time_traces('toroidal_current').values
        if units=='mks':
            y_axis = unit_conv(y_axis, arr_dim='M3DC1', current=1)
            plt.ylabel(r'Current $[A]$')

    ## Formally not right, we would like BT as a function of time
    if trace=='bt':
        Ave_P  = sim.get_time_traces('Ave_P').values
        B0     = sim.get_constants().B0
        y_axis = 2*Ave_P/B0**2
        
        
    if trace=='bp':
        if sim.get_constants().numvar < 3:
            raise Exception('Must have numvar=3 for a poloidal beta calculation!')
        gamma = sim.get_constants().gamma
        E_P = sim.get_time_traces('E_P').values
        tor_cur = sim.get_time_traces('toroidal_current').values
        y_axis  = 2*(gamma-1)*E_P/(tor_cur**2)
        
        
    if trace=='bn':
        raise Exception('Not coded yet')

    
    if trace=='beta':
        gamma = sim.get_constants().gamma
        E_P   = sim.get_time_traces('E_P').values
        E_MP  = sim.get_time_traces('E_MP').values
        E_MT  = sim.get_time_traces('E_MT').values
        y_axis= (gamma-1)*E_P/(E_MP + E_MT)
        
        
    if trace=='psibound' or trace=='psilim':
        y_axis = sim.get_time_traces('psi_lcfs').values
        if units=='mks':
            y_axis = unit_conv(y_axis, arr_dim='M3DC1', magnetic_field=1, length=2)
            plt.ylabel(r'Current $[T][m]^2$')

    
    if trace=='psimin':
        y_axis = sim.get_time_traces('psimin').values
        if units=='mks':
            y_axis = unit_conv(y_axis, arr_dim='M3DC1', magnetic_field=1, length=2)
            plt.ylabel(r'Current $[T][m]^2$')
    

    if trace=='dt':
        y_axis = sim.get_time_traces('dt').values
        if units=='mks':
            y_axis = unit_conv(y_axis,time=1)
            plt.ylabel(r'Time $[s]$')
        
    
    if trace=='volume':
        y_axis = sim.get_time_traces('volume').values
        if units=='mks':
            y_axis = unit_conv(y_axis,length=3)
            plt.ylabel(r'Volume $[m]^3$')


    if trace=='toroidal flux':
        y_axis = sim.get_time_traces('toroidal_flux_p').values
        if units=='mks':
            y_axis = unit_conv(y_axis,magnetic_field=1,length=2)
            plt.ylabel(r'Current $[T][m]^2$')


    if trace=='reconnected flux':
        y_axis = sim.get_time_traces('reconnected_flux').values
        if units=='mks':
            y_axis = unit_conv(y_axis,magnetic_field=1,length=2)
            plt.ylabel(r'Current $[T][m]^2$')
    
    
    if trace=='temax':
        y_axis = sim.get_time_traces('temax').values
        if units=='mks':
            y_axis = unit_conv(y_axis,temperature=1)
            plt.ylabel(r'Temperature $[eV]$')
    

    if trace=='p':
        y_axis = sim.get_time_traces('E_P').values
        if units=='mks':
            y_axis = unit_conv(y_axis,energy=1)
            plt.ylabel(r'Pressure $[J]$')


    if trace=='pe':
        y_axis = sim.get_time_traces('E_PE').values
        if units=='mks':
            y_axis = unit_conv(y_axis,energy=1)
            plt.ylabel(r'Electron pressure $[J]$')

    ## IDL has L=3, but isn't this already integrated over volume?
    if trace=='n':
        y_axis = sim.get_time_traces('particle_number').values
        if units=='mks':
            y_axis = unit_conv(y_axis,particles=1)
            plt.ylabel(r'Number of particles')

    ## IDL has L=3, but isn't this already integrated over volume?
    if trace=='ne':
        y_axis = sim.get_time_traces('electron_number').values
        if units=='mks':
            y_axis = unit_conv(y_axis,particles=1)
            plt.ylabel(r'Number of particles')

    
    if trace=='angular momentum':
        y_axis = sim.get_time_traces('angular_momentum').values
        if units=='mks':
            y_axis = unit_conv(y_axis,energy=1,time=1)
            plt.ylabel(r'Angular momentum $[J][s]$')


    if trace=='vorticity':
        y_axis = sim.get_time_traces('circulation').values
        if units=='mks':
            y_axis = unit_conv(y_axis,length=2,time=-1)
            plt.ylabel(r'Vorticity $[m]^2/[s]$')

            
    if trace=='bwb2':
        amupar = sim.get_constants().amupar
        y_axis = 3.0/(4.0*amupar)*sim.get_time_traces('parallel_viscous_heating').values
        if units=='mks':
            # Not too sure about the units here...
            y_axis = unit_conv(y_axis,energy=1,time=-2,length=3)
            plt.ylabel(r'Please verify units in code!')


    if trace=='li':
        R0       = sim.get_constants().R0
        psi_lcfs = sim.get_time_traces('psi_lcfs').values
        psimin   = sim.get_time_traces('psimin').values
        tor_cur  = sim.get_time_traces('toroidal_current_p').values
        y_axis   = -4.0 * np.pi * (psi_lcfs - psimin) / (R0*tor_cur)
        plt.ylabel(r'Normalized inductance')


    if trace=='li3':
        R0       = sim.get_constants().R0
        W_M      = sim.get_time_traces('W_M').values
        tor_curp = sim.get_time_traces('toroidal_current_p').values
        y_axis   = 4.0 * W_M / (tor_curp**2 * R0)
        plt.ylabel(r'Normalized inductance')


    if trace=='flux':
        psi_lcfs = sim.get_time_traces('psi_lcfs').values
        psimin   = sim.get_time_traces('psimin').values
        y_axis   = -2*np.pi*(psi_lcfs-psimin)
        if units=='mks':
            y_axis = unit_conv(y_axis,length=2,magnetic_field=-1)
            plt.ylabel(r'Normalized inductance')


    if trace=='xmag':
        y_axis = sim.get_time_traces('xmag').values
        if units=='mks':
            y_axis = unit_conv(y_axis,length=1)
            plt.ylabel(r'R coordinate of magnetic axis $[m]$')


    if trace=='zmag':
        y_axis = sim.get_time_traces('zmag').values
        if units=='mks':
            y_axis = unit_conv(y_axis,length=1)
            plt.ylabel(r'Z coordinate of magnetic axis $[m]$')


    if trace=='runaways':
        y_axis = sim.get_time_traces('runaways').values
        if units=='mks':
            y_axis = unit_conv(y_axis,particles=1)
            plt.ylabel(r'Number of runaways')


    if trace=='radiation':
        y_axis = sim.get_time_traces('radiation').values
        # Not too sure about dimensions
        if units=='mks':
            y_axis = unit_conv(y_axis,energy=1,time=-1)
            plt.ylabel(r'Radation losses $[J]/[s]$')

            
    if trace=='POhm':
        e_mpd = sim.get_time_traces('E_MPD').values
        e_mtd = sim.get_time_traces('E_MTD').values
        y_axis= -(e_mpd+e_mtd)
        # Not too sure about dimensions
        if units=='mks':
            y_axis = unit_conv(y_axis,energy=1,time=-1)
            plt.ylabel(r'Ohmic heating $[J]/[s]$')


    if trace=='pelr':
        y_axis = sim.get_time_traces('pellet_rate').values
        # Not too sure about dimensions
        if units=='mks':
            y_axis = unit_conv(y_axis,particles=1,time=-1)
            plt.ylabel(r'Pellet rate $1/[s]$')


    if trace=='pelvar':
        y_axis = sim.get_time_traces('pellet_var').values
        if units=='mks':
            y_axis = unit_conv(y_axis,length=1)
            plt.ylabel(r'Pellet variance $[m]$')


    if trace=='pelrad':
        if sim.get_constants().version < 26:
            y_axis = sim.get_time_traces('r_p2').values # RP2 is my fav space
        else:
            y_axis = sim.get_time_traces('r_p').values
        if units=='mks':
            y_axis = unit_conv(y_axis,length=1)
            plt.ylabel(r'Pellet radius $[m]$')


    if trace=='pelrpos':
        if sim.get_constants().version < 26:
            y_axis = sim.get_time_traces('pellet_x').values
        else:
            y_axis = sim.get_time_traces('pellet_r').values
        if units=='mks':
            y_axis = unit_conv(y_axis,length=1)
            plt.ylabel(r'Pellet R-coordinate $[m]$')


    if trace=='pelzpos':
        y_axis = sim.get_time_traces('pellet_z').values
        if units=='mks':
            y_axis = unit_conv(y_axis,length=1)
            plt.ylabel(r'Pellet Z-coordinate $[m]$')


    if trace=='bharmonics':
        y_axis = sim.get_time_traces('bharmonics').values
        time   = sim.get_time_traces('bharmonics').time
        if units=='mks':
            y_axis = unit_conv(y_axis,energy=1)
            time   = unit_conv(time,time=1)
            plt.ylabel(r'Magnetic energy harmonics $[J]$')
            plt.xlabel(r'Time $[s]$')
        for n in range(0,np.shape(y_axis)[1]):
            plt.plot(time,y_axis[:,n],label='n = {}'.format(n))
            plt.legend(loc='best')
        plt.show()
        sys.exit(0)

    if trace=='keharmonics':
        y_axis = sim.get_time_traces('keharmonics').values
        time   = sim.get_time_traces('keharmonics').time
        if units=='mks':
            y_axis = unit_conv(y_axis,energy=1)
            time   = unit_conv(time,time=1)
            plt.ylabel(r'Kinetic energy harmonics $[J]$')
            plt.xlabel(r'Time $[s]$')
        for n in range(0,np.shape(y_axis)[1]):
            plt.plot(time,y_axis[:,n],label='n = {}'.format(n))
            plt.legend(loc='best')
        plt.show()
        sys.exit(0)



    

    # If one array has new data but the other one doesn't 
    # plot only previous data
    if y_axis.shape != time.shape:
        ymax   = np.amax(y_axis.shape)
        tmax   = np.amax(time.shape)
        maxidx = np.amin([ymax,tmax])
        time   = time[0:maxidx]
        y_axis = y_axis[0:maxidx]

    plt.plot(time,y_axis)
    plt.ticklabel_format( axis='y', style='sci',useOffset=False)
    plt.show()
