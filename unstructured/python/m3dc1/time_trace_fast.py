#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 14:24:12 2019

@author: Andreas Kleiner
"""
import sys
import numpy as np
import os
import glob
import re
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rc
import matplotlib.gridspec as gridspec
from scipy.fftpack import fft, ifft, fftshift, fftfreq
from scipy import signal

import fpy
import m3dc1.fpylib as fpyl
from m3dc1.unit_conv import unit_conv
from m3dc1.plot_field  import plot_field
from m3dc1.read_h5 import readC1File, readParameter, openH5File
from m3dc1.gamma_file import Gamma_file
from m3dc1.gamma_data import Gamma_data
from m3dc1.flux_average import flux_average
from m3dc1.flux_coordinates import flux_coordinates
rc('text', usetex=True)
plt.rcParams.update({'figure.max_open_warning': 40})


def get_timetrace(trace,file_name='C1.h5',h5file=None,ipellet=0,
                  units='m3dc1',growth=False,renorm=False,quiet=False):
    """
    Read a time trace directly from an hdf5 file. This function does not use fusion-io.
    
    Arguments:

    **trace**
    Name of trace (scalar)

    **file_name**
    Name or path to C1.h5 file to read

    **units**
    The units in which the time trace will be returned

    **growth**
    If True, return growth rate of trace. If false, return trace

    **renorm**
    Remove spikes due to renormalization in linear runs that happens when
    the kinetic energy gets too large.

    **quiet**
    If True, do not print renormalization times to screen.
    """

    if h5file is None:
        h5file = openH5File(file_name)

    itor = readParameter('itor', h5file=h5file)
    version = readParameter('version', h5file=h5file)
    gamma = readParameter('gam', h5file=h5file)

    # Direct transformation of one name to another
    transform = {'toroidal current':'toroidal_current', 'it':'toroidal_current',
                 'plasma current':'toroidal_current_p', 'ip':'toroidal_current_p',
                 'wall current':'toroidal_current_w', 'iw':'toroidal_current_w',
                 'volume':'volume_p', 'plasma volume':'volume_p',
                 'volume_d':'volume', 'domain volume':'volume',
                 'toroidal flux': 'toroidal_flux_p',
                 'time step': 'dt',
                 'psibound':'psi_lcfs', 'psilim':'psi_lcfs',
                 'loop voltage':'loop_voltage', 'vl':'loop_voltage',
                 'poloidal magnetic energy':'E_MP', 'Wm':'E_MP',
                 'thermal energy':'E_P', 'p':'E_P',
                 'electron thermal energy':'E_PE', 'pe':'E_PE',
                 'particles':'particle_number', 'n':'particle_number',
                 'electrons':'electron_number', 'ne':'electron_number',
                 'angular momentum': 'angular_momentum',
                 'vorticity': 'circulation',
                 'parallel viscous heating': 'parallel_viscous_heating',
                 'IZ': 'M_IZ',
                 'ave_p': 'Ave_P',
                 'total current': 'itot',
                 'halo current': 'ih',
                 'kinetic energy': 'ke',
                 'magnetic energy': 'me',
                 'prad': 'radiation',
                 'pline': 'line_rad',
                 'pbrem': 'brem_rad',
                 'pion': 'ion_loss',
                 'preck': 'reck_rad',
                 'precp': 'recp_rad',
                 'POhm': 'pohm',
                 'pelr': 'pellet_rate',
                 'pellet ablation rate': 'pellet_ablrate',
                 'pelablr': 'pellet_ablrate',
                 'pellet var': 'pellet_var',
                 'pelvar': 'pellet_var',
                 'pellet radius': 'r_p',
                 'pelrad': 'r_p',
                 'pellet R position': 'pellet_r',
                 'pelrpos': 'pellet_r',
                 'pellet_x': 'pellet_r',
                 'pellet phi position': 'pellet_phi',
                 'pelphipos': 'pellet_phi',
                 'pellet Z position': 'pellet_z',
                 'pelzpos': 'pellet_z',
                 'poloidal beta': 'betap',
                 'bp': 'betap',
                 }

    # Simple linear combinations
    combos = {'itot':([('toroidal_current_w',1.),('toroidal_current',1.)],
                      {'current':1}),
              'ih':([('toroidal_current',1.),('toroidal_current_p',-1.)],
                    {'current':1}),
              'ke':([('E_KP',1.),('E_KT',1.),('E_K3',1.)], {'energy':1}),
              'me':([('E_MP',1.),('E_MT',1.)], {'energy':1}),
              'energy':([('E_KP',1.),('E_KT',1.),('E_K3',1.),
                         ('E_MP',1.),('E_MT',1.),('E_P',1.)], {'energy':1}),
              'flux':([('psimin',2.*np.pi),('psi_lcfs',-2.*np.pi)],
                      {'magnetic_field':1,'length':2}),
              'radiation':([('radiation',-1.)], None),
              'line_rad':([('line_rad',-1.)], None),
              'brem_rad':([('brem_rad',-1.)], None),
              'ion_loss':([('ion_loss',-1.)], None),
              'reck_rad':([('reck_rad',-1.)], None),
              'recp_rad':([('recp_rad',-1.)], None),
              'rec_rad':([('reck_rad',-1.),('recp_rad',-1.)],
                         {'energy':1,'time':-1}),
              'pohm':([('E_MPD',-1),('E_MTD',-1)], {'energy':1,'time':-1}),
              }

    if trace in transform:
        trace = transform[trace]

    if trace == 'reconnected flux':
        time,scalar = readC1File(scalar='reconnected_flux',h5file=h5file)
        scalar = np.abs(scalar)
        custom = {'magnetic_field':1,'length':1+itor}

    elif trace == 'r_p':
        if (version < 26):
            time,scalar = readC1File(scalar='r_p2',h5file=h5file)
        else:
            time,scalar = readC1File(scalar='r_p',h5file=h5file)
        custom = None

    elif trace == 'pellet_r':
        if (version < 26):
            time,scalar = readC1File(scalar='pellet_x',h5file=h5file)
        else:
            time,scalar = readC1File(scalar='pellet_r',h5file=h5file)
        custom = None

    elif trace == 'beta':
        if (version < 26):
            time,scalar = readC1File(scalar='E_P',h5file=h5file)
        else:
            time,scalar = readC1File(scalar='W_P',h5file=h5file)
        t, E_MP= readC1File(scalar='E_MP',h5file=h5file)
        t, E_MT = readC1File(scalar='E_MT',h5file=h5file)
        scalar *= (gamma-1.)/(E_MP + E_MT)
        custom = None

    elif trace == 'betap':
        if (version < 26):
            time,scalar = readC1File(scalar='E_P',h5file=h5file)
            t,it = readC1File(scalar='toroidal_current',h5file=h5file)
            scalar *= 2.*(gamma-1.)/it**2
        else:
            time,scalar = readC1File(scalar='W_P',h5file=h5file)
            t, W_M = readC1File(scalar='W_M',h5file=h5file)
            scalar *= (gamma-1.)/W_M
        custom = None

    elif trace in ['betan','betat']:
        raise RuntimeError("'%s' not yet implemented; need shape information"%trace)

    elif trace == 'electron_number':
        if version <= 20:
            zeff = readParameter('zeff', h5file=h5file)
            time,scalar = readC1File(scalar='particle_number',h5file=h5file)
            scalar *= zeff
        else:
            time,scalar = readC1File(scalar=trace,h5file=h5file)

        custom = {'particles':1}

    elif trace == 'bwb2':
        amupar = readParameter('amupar', h5file=h5file)
        time,scalar = readC1File(scalar='parallel_viscous_heating',h5file=h5file)
        scalar *= 4./(3.*amupar)
        custom = {'length':3, 'time':-2}

    elif trace == 'li':
        rzero = readParameter('rzero', h5file=h5file)
        time,psi_lcfs = readC1File(scalar='psi_lcfs',h5file=h5file)
        time,psimin = readC1File(scalar='psimin',h5file=h5file)
        time,ip = readC1File(scalar='toroidal_current_p',h5file=h5file)

        scalar = -4.*np.pi*(psi_lcfs - psimin)/(rzero*ip)
        custom = None

    elif trace == 'li3':
        rzero = readParameter('rzero', h5file=h5file)
        time,W_M = readC1File(scalar='W_M',h5file=h5file)
        time,ip = readC1File(scalar='toroidal_current_p',h5file=h5file)

        scalar = 4.*W_M/(rzero*ip**2)
        custom = None

    elif trace in combos:
        # trace is linear combination of native scalars
        combo, custom = combos[trace]
        for i, (name, fac) in enumerate(combo):
            t, y = readC1File(scalar=name,h5file=h5file)
            if i==0:
                time = t
                scalar = fac*y
            else:
                scalar += fac*y

    else:
        # trace is a native scalar
        time,scalar = readC1File(scalar=trace,h5file=h5file)
        custom = None


    if ('pellet_' in trace) or (trace in ['cauchy_fraction','cloud_pel','r_p']):
        # if ipellet is given, get just that pellet's data
        if (ipellet != 'all') and (scalar.ndim==2):
            scalar = scalar[ipellet,:]

    if units=='mks':
        time = unit_conv(time, arr_dim='M3DC1', h5file=h5file, time=1)
        scalar = fpyl.get_conv_trace('mks',trace,scalar,h5file=h5file,itor=itor,custom=custom)
        
    if growth == True:
        scalar = 1.0/scalar[1:] * np.diff(scalar)/np.diff(time)
        time = time[:-1]
    
    if renorm == True:
        renormlist = []
        for i in range(len(scalar)-1):
            if(abs(scalar[i+1]/scalar[i]) < 1E-9):
                renormlist.append(str(time[i]))
                #print(scalar[i],scalar[i-1]+scalar[i+1])
                # Only average value if growth rate is calculated
                if growth == True:
                    scalar[i] = (scalar[i-1] + scalar[i+1])/2.0
        # When growth rate is calculated, check for normalization at last time
        #   step and drop this point, since it carres no information.
        if growth == True:
            if(abs(scalar[-2]/scalar[-1]) < 1E-9):
                renormlist.append(str(time[-1]))
                scalar = scalar[:-1]
                time = time[:-1]
        renormstr = ", ".join(renormlist)
        if quiet==False:
            if len(renormstr) > 0:
                print('Renormalization found at '+renormstr)
    return time, scalar



def avg_time_trace(trace,units='m3dc1',file_name='C1.h5',h5file=None,
                   growth=False,renorm=True,start=None,time_low_lim=500):
    """
    Calculates the mean and standard deviation of a M3DC1 scalar (time trace) starting from a certain point in time
    
    Arguments:

    **trace**
    Name of trace (scalar)

    **file_name**
    Name or path to C1.h5 file to read

    **units**
    The units in which the time trace will be returned

    **growth**
    If True, return growth rate of trace. If false, return trace

    **renorm**
    Remove spikes due to renormalization in linear runs that happens when
    the kinetic energy gets too large.

    **start**
    time at which averaging begins

    **time_low_lim**
    lower limit for starting time in terms of Alfven times. If start < time_low_lim,
    start time is moved to the earliest time larger than time_low_lim
    """

    if h5file is None:
        h5file = openH5File(file_name)
    time,scalar = get_timetrace(trace,h5file=h5file,units=units,growth=growth,
                                renorm=renorm,quiet=False)
    
    if start==None:
        start_ind = int(np.floor(len(scalar)/2))
        start_time = time[start_ind]
    else:
        start_time = start
        start_ind = int(fpyl.find_nearest(time,start_time))
        print(start_time,start_ind)
    
    
    if units.lower() == 'mks':
        time_low_lim = unit_conv(time_low_lim,arr_dim='m3dc1',h5file=h5file,time=1)
    
    if start_time < time_low_lim:
        if time_low_lim < time[-1]:
            start_ind = np.argmax(time>time_low_lim)
            fpyl.printwarn('WARNING: Start of trace averaging has been moved to t='+str(time[start_ind])+'.')
        else:
            fpyl.printwarn('WARNING: time_low_lim > time[-1]. Start of trace averaging has been moved to t='+str(time[start_ind])+'. Please verify validity of results.')
    
    scalar_short = scalar[start_ind:]
    time_short = time[start_ind:]
    
    avg = np.mean(scalar_short)
    std = np.std(scalar_short)
    return avg, std,time_short,scalar_short


def growth_rate(n=None,units='m3dc1',file_name='C1.h5',h5file=None,
                time_low_lim=500,slurm=True,plottrace=False):
    """
    Evaluates kinetic energy growth rate. The growth rate is the mean of the logarithmic derivative of ke.
    The mean is taken over the second half of the simulation time. Returns 
    
    Arguments:

    **units**
    System of units for the plot, can be either 'm3dc1' or 'mks'

    **time_low_lim**
    Minimum time in Alfven times that will be considered for averaging.
    Should be chosen such that the noise in the beginning of the simulation is avoided

    **slurm**
    If True, reads the Slurm log file for the M3DC1 run and checks for GS errors and convergence
    
    **plottrace**
    Show and save plots of growth rate in directory
    """
    if h5file is None:
        h5file = openH5File(file_name)

    if n==None:
        n = readParameter('ntor',h5file=h5file,listc=False)
        fpyl.printnote('Set n=ntor='+"{:d}".format(n)+' as read from C1.h5 file.')
    gamma, dgamma,time,gamma_trace = avg_time_trace('ke',units,h5file=h5file,
                                                    growth=True,renorm=True,
                                                    start=None,
                                                    time_low_lim=time_low_lim)
    print(n,gamma,dgamma)
    
    not_noisy = 1
    gamma_set_manu = 0
    flat = 1
    
    if abs(dgamma/gamma) > 0.1:
        fpyl.printwarn('WARNING: gamma is not constant for n='+str(n)+'!')
        flat = 0
        
        # Identify local maxima in the growth rate and determine the frequency in which these occur
        maxima_ind = signal.argrelmax(gamma_trace)[0]
        frequency = 1.0/np.diff(maxima_ind)
        
        # When frequency is too high, mask it as noise
        noise_mask = np.greater(frequency,0.1)
        
        perform_manual_check = False
        
        # Check whether growth rate is noisy for the whole time, partially,
        #   or not noisy. In case of noise, do not consider the growth rate.
        if np.all(noise_mask) == True:
            fpyl.printwarn('WARNING: gamma is completely noisy for n='+str(n)+'!')
            not_noisy = 0
        elif np.any(noise_mask) == True:

            fpyl.printwarn('WARNING: gamma is partially noisy for n='+str(n)+'!')
            if (gamma < 0 and (np.all(np.sign(gamma_trace)==-1.0))):
                perform_manual_check = True
            not_noisy = 0
        else:
            # If there is no noise, do additional check ups to determine if
            #   the system is stable or unstable
            fpyl.printnote('NOTE: no noise detected for n='+str(n)+'!')
            # Check whether the kinetic energy is overall increasing or
            #   decreasing by looking at the evolution of the maxima.
            ke_time,ke = get_timetrace('ke',h5file=h5file,units=units,
                                       growth=False,renorm=True)
            
            start_ind = int(np.floor(len(ke)/2))
            start_time = ke_time[start_ind]
            if units.lower() == 'mks':
                time_low_lim = unit_conv(time_low_lim,arr_dim='m3dc1',
                                         h5file=h5file,time=1)
            if start_time < time_low_lim:
                start_ind = np.argmax(ke_time>time_low_lim)
            ke_short = ke[start_ind:]
            ke_time_short = ke_time[start_ind:]
            
            # Identify renormalizations in the M3D-C1 run
            renormlist = []
            for i in range(len(ke_short)-1):
                if(abs(ke_short[i+1]/ke_short[i]) < 1E-9):
                    renormlist.append(i)
            if(abs(ke_short[-2]/ke_short[-1]) < 1E-9):
                renormlist.append(len(ke_short)-1)
            ke_linear = []
            time_linear = []
            
            # Get kinetic energy traces in between the renormalizations
            if len(renormlist)>0:
                for i,t in enumerate(renormlist):
                    if i==0:
                        ke_linear.append(ke_short[:t])
                        time_linear.append(ke_time_short[:t])
                    elif i>0:
                        ke_linear.append(ke_short[renormlist[i-1]+1:t])
                        time_linear.append(ke_time_short[renormlist[i-1]+1:t])
                ke_linear.append(ke_short[renormlist[-1]+1:])
                time_linear.append(ke_time_short[renormlist[-1]+1:])

            else:
                ke_linear = [np.asarray(ke_short)]
            
            # Identify maxima in the kinetic energy
            ke_maxima_ind = signal.argrelmax(ke_short)[0]
            
            # Check if ke is monotonic in between all renormalizations
            if ((all(fpyl.strict_monotonic(kel) == True for kel in ke_linear)) or
                (len(ke_maxima_ind)<2)):
                # If it is monotonic (which is a good sign), check the growth rate manually.
                #   This is to ensure that oscillations are not too strong,
                #   or the value for the growth rate is extracted correctly.
                perform_manual_check = True
            else:
                # Look for local maxima and minima among the maxima to identify
                #   a transition from decay to growth or vice versa
                max_ke = ke_short[ke_maxima_ind]
                ke_time_max = ke_time_short[ke_maxima_ind]
            
                max_ke_max = signal.argrelmax(max_ke)[0]
                max_ke_min = signal.argrelmin(max_ke)[0]
                
                if len(max_ke_max) > 0:
                    print(max_ke_max)
                if len(max_ke_min) > 0:
                    print(max_ke_min)
                
                # Calculate time derivative of ke maxima
                dke_dt = np.gradient(max_ke,ke_time_max)
                
                
                # make different decisions depending on the sign of the
                #   growth rate and time trace of ke
                if ((np.sign(gamma) == np.sign(dke_dt[-1])) and
                    (len(max_ke_min)==0 and len(max_ke_max)==0)):
                    # ke and growth rate have same sign, and sign of ke maxima
                    #   does not change
                    # Now check whether growth rate changes sign
                    if (np.all(np.sign(gamma_trace)==1.0) and gamma > 0):
                        fpyl.printnote('NOTE: No change of sign found. Plasma appears linearly UNSTABLE.')
                        perform_manual_check = False
                    elif (np.all(np.sign(gamma_trace)==-1.0) and gamma < 0):
                        fpyl.printnote('NOTE: No change of sign found. Plasma appears linearly STABLE.')
                        perform_manual_check = False
                    else:
                        fpyl.printwarn('WARNING: growth rate crosses zero.')
                        perform_manual_check = True
                else:
                    fpyl.printwarn('WARNING: This case needs some attention.')
                    perform_manual_check = True
        

        # Prompts for manual check of growth rate
        if perform_manual_check==True:
            # Show plots of growth rate and kinetic energy for analysis
            double_plot_time_trace_fast('ke',renorm=True,title='',rescale=True,units=units)
            mono_input = fpyl.prompt('Is the calculated growth rate acceptable? (y/n) : ',['y','n'])
            if mono_input == 'y':
                gamma_set_manu = 1
            else:
                gr_input = fpyl.prompt('Please choose:\n  [1] set a different start time for mean calculation\n  [2] directly enter a value for gamma\n  [3] discard this simulation\n : ',['1','2','3'])
                if gr_input == '1':
                    #ToDo: Test averaging from start_time in various units using
                    #         various values for time and limits
                    start_time = fpyl.prompt('Please enter a start time for mean calculation : ',float)
                    gamma, dgamma,_,_ = avg_time_trace('ke',units,h5file=h5file,
                                                       growth=True,renorm=True,
                                                       start=start_time,
                                                       time_low_lim=time_low_lim)
                    gamma_set_manu = 1
                    print('gamma corrected to new avg. value: '+str(gamma))
                elif gr_input == '2':
                    gamma = fpyl.prompt('Please enter a value for the growth rate : ',float)
                    gamma_set_manu = 1
                    print('gamma corrected to value: '+str(gamma))
                elif gr_input == '3':
                    not_noisy = 0


    #------------------------------------------------
    #Check equilibrium convergence and extract Final GS Error from Slurm log file
    #------------------------------------------------
    if slurm == True:
        try:
            # Read Slurm log files. Should there be multiple log files in the directory,
            #   choose the file with largest Slurm Job ID
            slurmfiles = glob.glob(os.getcwd()+"/slurm*.out")
            if len(slurmfiles) > 1:
                fpyl.printwarn('WARNING: More than 1 Slurm log file found. Using the latest one.')
                slurmfiles.sort()
                #print(slurmfiles)
            slurmfile = slurmfiles[-1]
            
            # Read data from Slurm log file
            maxit = 0
            with open(slurmfile, 'r') as sf:
                for line in sf:
                    if 'igs   ' in line:
                        igsline = line.split()
                        igs = int(igsline[2])
                    elif 'tol_gs' in line:
                        tolgsline = line.split()
                        tolgs = float(tolgsline[2])
                    elif 'GS iteration' in line:
                        gsiterline = line.split()
                        #print(gsiterline)
                        maxit = int(gsiterline[3])
                        gserr = float(gsiterline[4])
                    elif 'Final error in GS solution' in line:
                        fgserrline = line.split()
                        finalerrgs = float(fgserrline[-1])
                        break
            # Check whether the GS solution converged
            if maxit > 0:
                if maxit < igs:
                    gsconvgd = 1
                elif maxit == igs:
                    if gserr / tolgs < 3:
                        fpyl.printwarn('WARNING: Please check GS solution. It might not be converged.')
                        gsconvgd = 1
                    else:
                        gsconvgd = 0
                else:
                    gsconvgd = 0
            else:
                fpyl.printwarn('WARNING: Number of GS iterations could not be determined.')
                gsconvgd = 0
            return gamma, dgamma, n, flat, not_noisy, gamma_set_manu, gsconvgd, finalerrgs
        except:
            fpyl.printerr('ERROR: Could not process Slurm log file!')
            return gamma, dgamma, n, flat, not_noisy, gamma_set_manu, None, None
    else:
        return gamma, dgamma, n, flat, not_noisy, gamma_set_manu, None, None
    



def scan_n(nmin=1,nmax=10,units='m3dc1',file_name='C1.h5',h5file=None,time_low_lim=500,slurm=True,plottrace=False):
    """
    Traverses all subdirectories named nXX in a directory (where XX is the toroidal mode number),
    and reads the growth rate.
    
    Arguments:

    **nmin**
    Lowest mode number to include in the scan

    **nmax**
    Largest mode number to include in the scan

    **units**
    System of units for the plot, can be either 'm3dc1' or 'mks'

    **time_low_lim**
    Minimum time in Alfven times that will be considered for averaging.
    Should be chosen such that the noise in the beginning of the simulation is avoided

    **slurm**
    If True, reads the Slurm log file for the M3DC1 run and checks for GS errors and convergence
    
    **plottrace**
    Show and save plots of growth rate in directory
    """

    if h5file is None:
        h5file = openH5File(file_name)
    
    n_list = []
    gamma_list = []
    dgamma_list = []
    flat_list = []
    not_noisy_list = []
    gamma_manual_list = []
    
    gsconvgd_list = []
    finalerrgs_list = []
    print('Calculating growth rates for n='+str(nmin)+' to n='+str(nmax))
    
    for n in range(nmin,nmax+1):
        if n==nmin:
            path = 'n'+str(n).zfill(2)
        else:
            path = '../n'+str(n).zfill(2)
        os.chdir(path)
        
        print('----------------------------------------------')
        print('n='+str(n)+':')

        gamma, dgamma, n, flat, not_noisy, gamma_set_manu, gsconvgd, finalerrgs = growth_rate(n,units=units,h5file=h5file,time_low_lim=time_low_lim,slurm=slurm,plottrace=True)
        
        gamma_list.append(gamma)
        dgamma_list.append(dgamma)
        n_list.append(n)
        flat_list.append(flat)
        not_noisy_list.append(not_noisy)
        gamma_manual_list.append(gamma_set_manu)
        gsconvgd_list.append(gsconvgd)
        finalerrgs_list.append(finalerrgs)
        
        
        if n==nmax:
            print('----------------------------------------------')
            os.chdir('../')
            
    results = Gamma_data(n_list, gamma_list, dgamma_list, flat_list, not_noisy_list, gamma_manual_list, gsconvgd_list, finalerrgs_list)
    
    return results



def create_plot_gamma_n(n_list, gamma_list,fignum=None,lw=1,c=None,ls=None,marker=None,ms=36,lbl=None,units='m3dc1'):
    plt.figure(num=fignum)
    plt.plot(n_list,gamma_list,lw=lw,c=c,ls=ls,marker=marker,ms=ms,label=lbl)
    plt.grid(True)
    plt.xlabel('n',fontsize=18)
    if units.lower()=='m3dc1':
        plt.ylabel(r'$\gamma/\omega_A$',fontsize=18)
    else:
        plt.ylabel(r'$\gamma$ $[s^{-1}]$',fontsize=18)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.tight_layout()
    if lbl != None:
        plt.legend(loc=0)
    return



def plot_gamma_n(nmin=1,nmax=10,units='m3dc1',fignum=None,c=None,ls=None,lbl=None,slurm=True,plottrace=False):
    files = glob.glob("growth_rates*.dat")
    
    if len(files)>0:
        openf_input = fpyl.prompt('Previous results found. Do you want to open the existing file? (y/n) : ',['y','n'])
        if openf_input == 'y':
            if len(files)>1:
                    files.sort(key=os.path.getmtime)
                    print('More than 1 file found. Opening newest file.')
            f = files[0]
            results = Gamma_file(f)
        else:
            results = scan_n(nmin,nmax,units=units,slurm=True,plottrace=plottrace)
    else:
        results = scan_n(nmin,nmax,units=units,slurm=True,plottrace=plottrace)
    
    # Identify simulations where the growth rate was not calculated reliably. These are highlighted in the plot.
    n_okay = []
    gamma_okay = []
    n_bad = []
    gamma_bad = []
    for n in results.n_list:
        #n_ind = n-1
        n_ind = fpyl.get_ind_at_val(results.n_list,n,unique=True)
        if (results.flat_list[n_ind] != 1 and results.gamma_manual_list[n_ind] != 1 and results.not_noisy_list[n_ind] != 1):
            n_bad.append(results.n_list[n_ind])
            gamma_bad.append(results.gamma_list[n_ind])
        else:
            n_okay.append(results.n_list[n_ind])
            gamma_okay.append(results.gamma_list[n_ind])
    
    create_plot_gamma_n(results.n_list, results.gamma_list, fignum, c=c, ls=ls, marker='.', ms=8, lbl=lbl, units=units)
    if len(n_bad)>0:
        cfn = plt.gcf().number #Current figure number
        create_plot_gamma_n(n_bad, gamma_bad, fignum=cfn, lw=0, marker='.', ms=10, c='r', units=units)

    return



def write_gamma_n(nmin=1,nmax=10,results=None,units='m3dc1'):
    if (nmin >= 0 and nmax >=0 and nmax > nmin):
        #n_list, gamma_list, dgamma_list, not_noisy_list, gsconvgd, finalerrgs = scan_n(nmin,nmax,slurm=True)
        results = scan_n(nmin,nmax,units=units,slurm=True,plottrace=True)
    elif (nmin < 0 and nmax < 0) and (results==None):
        raise Exception('nmin and nmax must be >=0 with nmax > nmin.')
    
    pwd = os.getcwd()
    pathdirs = pwd.split('/')
    vpnum = pathdirs[-2]
    simdir = pathdirs[-1]
    
    eta = re.search('eta_x(\d*[.]?\d*)', simdir).group(1)
    
    flmodel = simdir[0:2]
    if flmodel != '1f' and flmodel != '2f':
        fpyl.printerr('ERROR: Fluid model not recognized!')
    
    if 'eqrotnc' in simdir:
        rotation=-1
    elif 'eqrot_' in simdir:
        rotation=1
    elif 'norot' in simdir:
        rotation=0
    else:
        fpyl.printwarn('WARNING: Equilibrium rotation unknown!')
    
    if 'B' in simdir:
        bscale = re.search('B(.*)_', simdir).group(1)
    else:
        bscale = '1.0'
    
    #try:
    #    jpdata = []
    #    with open('eq_pedestal_parameters.dat', 'r') as data:
    #        for line in data:
    #            jpdata.append(float(line))
    #except:
    #    fpyl.printerr('ERROR: Could not read current and pprime!')
    if nmin<0:
        ndir = 'n01'
    else:
        ndir = 'n' + str(nmin).zfill(2)
    os.chdir(ndir)
    sim0 = fpy.sim_data(time=0)
    sim0 = flux_coordinates(sim=sim0)
    
    alpha_max = np.amax(flux_average('alpha', coord='scalar', sim=sim0, fcoords='', units='m3dc1')[1]) #ToDo: Check if it's really a local maximum
    j_max = np.amax(flux_average('j', coord='phi', sim=sim0, fcoords='', units='m3dc1')[1]) #ToDo: Check if it's really a local maximum
    os.chdir('../')
    
    pbmode = 1
    
    
    #print(results.n_list)
    outfile = 'growth_rates_'+vpnum+'_'+simdir+'.dat'
    with open(outfile, 'w') as f:
        f.write(vpnum+'    '+str("{:f}".format(float(eta)))+'    '+bscale+'    '+str(rotation)+'    '+flmodel+'\n')
        #f.write(str(jpdata[1])+'    '+str(jpdata[0])+'\n') #write j and pprime
        f.write(str(j_max)+'    '+str(alpha_max)+'\n') #write j and pprime
    
        f.write('    n      gamma         sig_gamma     flat     smooth   manu     conv  Fin. GS Err   PB\n')
        for i in range(len(results.gamma_list)):
            wstr = '    ' + str("{:d}".format(results.n_list[i])).ljust(2,' ')
            wstr = wstr + str("{0:.8f}".format(results.gamma_list[i])).rjust(15,' ')
            wstr = wstr + str("{0:.8f}".format(results.dgamma_list[i])).rjust(14,' ')
            wstr = wstr + '    ' + str("{:d}".format(results.flat_list[i])).ljust(5,' ')
            wstr = wstr + '    ' + str("{:d}".format(results.not_noisy_list[i])).ljust(5,' ')
            wstr = wstr + '    ' + str("{:d}".format(results.gamma_manual_list[i])).ljust(5,' ')
            wstr = wstr + '    ' + str("{:d}".format(results.gsconvgd[i])).ljust(2,' ')
            wstr = wstr + '    ' + str("{0:.8f}".format(results.finalerrgs[i]))
            wstr = wstr + '    ' + str(pbmode)+'\n'
            f.write(wstr)
            
    print("Growth rates written to file '"+str(outfile)+"'.")
    
    return


def eval_growth_n(nmin=1,nmax=10,plotef=True,units='m3dc1'):
    #First check if pedestal parameter file is there
    #eq_ped_file = glob.glob("eq_pedestal_parameters.dat")
    #if len(eq_ped_file)!=1:
    #    raise Exception('eq_pedestal_parameters.dat not found!')
    
    
    results = scan_n(nmin,nmax,units=units,slurm=True,plottrace=True)
    #plt.close('all')
    write_gamma_n(-1, -1, results)
    
    
    # Plot gamma as a function of n
    # Identify simulations where the growth rate was not calculated reliably. These are highlighted in the plot.
    n_okay = []
    gamma_okay = []
    n_bad = []
    gamma_bad = []
    for n in results.n_list:
        n_ind = fpyl.get_ind_at_val(results.n_list,n,unique=True)
        if (results.flat_list[n_ind] != 1 and results.gamma_manual_list[n_ind] != 1 and results.not_noisy_list[n_ind] != 1):
            n_bad.append(results.n_list[n_ind])
            gamma_bad.append(results.gamma_list[n_ind])
        else:
            n_okay.append(results.n_list[n_ind])
            gamma_okay.append(results.gamma_list[n_ind])
    
    create_plot_gamma_n(results.n_list, results.gamma_list, fignum=None, marker='.', ms=8, lbl=None, units=units)
    cfn = plt.gcf().number #Current figure number
    create_plot_gamma_n(n_bad, gamma_bad, fignum=cfn, lw=0, marker='.', ms=10, lbl=None, units=units)
    
    
    
    for n in range(nmin,nmax+1):
        if n==nmin:
            path = 'n'+str(n).zfill(2)
        else:
            path = '../n'+str(n).zfill(2)
        os.chdir(path)
        
        if plotef==True:
            print('Plotting eigenfunction for n='+str(n))
            plot_field('p',time='last',linear=True,bound=True,lcfs=True,save=True,savedir='../')
            plt.close()
        if n==nmax:
            os.chdir('../')
    
    return





def create_plot_time_trace_fast(time,scalar,trace,units='mks',file_name='C1.h5',
                                h5file=None,growth=False,yscale='linear',
                                rescale=False,save=False,savedir=None,pub=False):

    if h5file is None:
        h5file = openH5File(file_name)

    # If one array has new data but the other one doesn't 
    # plot only previous data
    if scalar.shape != time.shape:
        ymax   = np.amax(scalar.shape)
        tmax   = np.amax(time.shape)
        maxidx = np.amin([ymax,tmax])
        time   = time[0:maxidx]
        scalar = scalar[0:maxidx]
    
    ntor = readParameter('ntor',h5file=h5file)
    
    plt.figure()
    
    # Set font sizes and plot style parameters
    if pub==False:
        axlblfs = None
        titlefs = None
        ticklblfs = None
        linew = 1
    elif pub==True:
        axlblfs = 20
        titlefs = 18
        ticklblfs = 18
        linew = 2
    
    
    if units=='mks':
        plt.xlabel(r'time $[s]$',fontsize=axlblfs)
    else:
        plt.xlabel(r'time $[\tau_A]$',fontsize=axlblfs)
    
    if trace == 'ke':
        if units=='mks':
            if growth == True:
                plt.ylabel(r'$\gamma$ $[s^{-1}]$')
            else:
                plt.ylabel(r'Kinetic energy $[J]$')
        elif units.lower()=='m3dc1':
            if growth == True:
                plt.ylabel(r'$\gamma/\omega_A$')
            else:
                plt.ylabel(r'Kinetic energy (M3DC1 units)')
    else:
        plt.ylabel(trace,fontsize=axlblfs)
    
    plt.plot(time,scalar,lw=linew)
    plt.grid(True)
    #Determine y-axis limits
    if rescale==True:
        if np.amax(scalar[1:]) < scalar[0]:
            start_time=250
            if units=='mks':
                start_time = unit_conv(start_time,arr_dim='m3dc1',h5file=h5file,time=1)
            start_ind = int(fpyl.find_nearest(time,start_time))
            top_lim=1.1*np.amax(scalar[start_ind:])
            plt.ylim([0,top_lim])
    plt.yscale(yscale)
    ax = plt.gca()
    ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
    plt.ticklabel_format( axis='y', style='sci',useOffset=False)
    plt.title('n='+str(ntor),fontsize=titlefs)
    plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    plt.show()
    
    if save == True:
        tracestr = trace
        if growth==True:
            tracestr = tracestr + '-growth'
        if savedir != None:
            tracestr = savedir + tracestr
        plt.savefig(tracestr+'_n'+"{:d}".format(ntor)+'.pdf', format='pdf',bbox_inches='tight')
    return



def double_plot_time_trace_fast(trace,file_name='C1.h5',h5file=None,
                                renorm=False,rescale=False,units='m3dc1',
                                title=None,pub=False):
    """
    Plots a trace and its growth rate in two subplots side by side
    
    Arguments:

    **trace**
    String containing the trace to be plotted

    **file_name**
    Contains the path to the C1.h5 file.

    **renorm**
    Removes spikes that are caused by renormalizations
    in linear stability calculations. Interpolates at
    the locations of the spike. Should only be used if
    growth=True.

    **rescale**
    Rescale y-axis such that noise in the beginning of the simulation is not considered for axis limits,
    i.e. plot is zoomed in to show relevant values.

    **units**
    The units in which the trace will be plotted

    **title**
    Plot title

    **pub**
    If True, format figure for publication (larger labels and thicker lines)
    """

    if h5file is None:
        h5file = openH5File(file_name)

    time,scalar = get_timetrace(trace,h5file=h5file,units=units,growth=False,renorm=renorm,quiet=True)
    time_growth,scalar_growth = get_timetrace(trace,h5file=h5file,units=units,growth=True,renorm=renorm,quiet=True)
    
    if units=='mks':
        xlbl = r'time $[s]$'
        ylbl_left = r'Kinetic energy $[J]$'
        ylbl_right = r'$\gamma$ $[s^{-1}]$'
    else:
        xlbl = r'time $[\tau_A]$'
        ylbl_left = r'Kinetic energy (M3DC1 units)'
        ylbl_right = r'$\gamma/\omega_A$'
    
    #Determine y-axis limits
    if rescale==True:
        if np.amax(scalar[1:]) < scalar[0]:
            start_time=250
            if units=='mks':
                start_time = unit_conv(start_time,arr_dim='m3dc1',h5file=h5file,time=1)
            start_ind = int(fpyl.find_nearest(time,start_time))
            top_lim=1.1*np.amax(scalar[start_ind:])
        else:
            top_lim=None
    
    fig = plt.figure(constrained_layout=True,figsize=(10,5))
    spec2 = gridspec.GridSpec(ncols=2, nrows=1, figure=fig)
    f2_ax1 = fig.add_subplot(spec2[0, 0])
    f2_ax2 = fig.add_subplot(spec2[0, 1])
    
    f2_ax1.plot(time,scalar)
    f2_ax1.set_xlabel(xlbl)
    f2_ax1.set_ylabel(ylbl_left)
    if rescale==True:
        f2_ax1.set_ylim([0,top_lim])
    f2_ax1.grid()
    f2_ax2.plot(time_growth,scalar_growth)
    f2_ax2.set_xlabel(xlbl)
    f2_ax2.set_ylabel(ylbl_right)
    f2_ax2.grid()
    
    ntor = readParameter('ntor',h5file=h5file)
    fig.suptitle('n='+str(ntor), size=20)
    
    return


def plot_time_trace_fast(trace,units='mks',file_name='C1.h5',h5file=None,
                         growth=False,renorm=False,yscale='linear',
                         rescale=False,save=False,savedir=None,pub=False):
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
    
    **rescale**
    Rescale y-axis such that noise in the beginning of the simulation is not considered for axis limits,
    i.e. plot is zoomed in to show relevant values.
    
    **save**
    If True, save plot to file
    
    **savedir**
    Directory where plot will be saved
    
    **pub**
    If True, format figure for publication (larger labels and thicker lines)
    """
    if h5file is None:
        h5file = openH5File(file_name)
    time,scalar = get_timetrace(trace,h5file=h5file,units=units,
                                growth=growth,renorm=renorm)
    create_plot_time_trace_fast(time,scalar,trace,units=units,h5file=h5file,
                                growth=growth,yscale=yscale,rescale=rescale,
                                save=save,savedir=savedir,pub=pub)
    return
