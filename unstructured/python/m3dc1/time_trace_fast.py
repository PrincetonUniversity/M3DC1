#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 14:24:12 2019

@author: Andreas Kleiner
"""

import numpy as np
import os
import glob
from pathlib import Path
import re
import math
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
#from scipy.fftpack import fft, ifft, fftshift, fftfreq
from scipy import signal
try:
    from scipy.integrate import trapezoid as trapz
    from scipy.integrate import cumulative_trapezoid as cumtrapz
except:
    from scipy.integrate import trapz
    from scipy.integrate import cumtrapz
from termcolor import colored
from labellines import labelLines

import fpy
import m3dc1.fpylib as fpyl
from m3dc1.unit_conv import unit_conv
from m3dc1.plot_field  import plot_field
from m3dc1.read_h5 import readC1File
from m3dc1.read_h5 import readParameter
from m3dc1.gamma_file import Gamma_file
from m3dc1.gamma_data import Gamma_data
from m3dc1.flux_average import flux_average
from m3dc1.get_time_of_slice import get_time_of_slice

from m3dc1.pedestal_finder      import get_ped_structure

from m3dc1.get_timetrace import get_timetrace

#from m3dc1.flux_coordinates import flux_coordinates
from m3dc1.eigenfunction import eigenfunction
from m3dc1.eigenfunction import mode_type

rc('text', usetex=True)
plt.rcParams.update({'figure.max_open_warning': 40})



def avg_time_trace(trace,units='m3dc1',sim=None,filename='C1.h5',
                   growth=False,renorm=True,start=None,time_low_lim=500):
    """
    Calculates the mean and standard deviation of a M3DC1 scalar (time trace) starting from a certain point in time
    
    Arguments:

    **trace**
    Name of trace (scalar)

    **units**
    The units in which the time trace will be returned

    **sim**
    fpy simulation object.

    **filename**
    Name or path to C1.h5 file to read.

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
    if isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename=filename)
    time,values,_,_ = get_timetrace(trace,sim=sim,units=units,growth=growth,
                                renorm=renorm,quiet=False)
    
    if start is None:
        start_ind = int(np.floor(len(values)/2))
        start_time = time[start_ind]
    else:
        start_time = start
        start_ind = int(fpyl.get_ind_near_val(time,start_time))
        print(start_time,start_ind)
    
    
    if units.lower() == 'mks':
        time_low_lim = unit_conv(time_low_lim,arr_dim='m3dc1',sim=sim,time=1)
    
    #print(start_time,time_low_lim)
    
    if start_time < time_low_lim:
        if time_low_lim < time[-1]:
            start_ind = np.argmax(time>time_low_lim)
            #print(start_ind)
            fpyl.printwarn('WARNING: Start of trace averaging has been moved to t='+str(time[start_ind])+'.')
        else:
            fpyl.printwarn('WARNING: time_low_lim > time[-1]. Start of trace averaging has been moved to t='+str(time[start_ind])+'. Please verify validity of results.')
    
    values_short = values[start_ind:]
    time_short = time[start_ind:]
    
    avg = np.mean(values_short)
    std = np.std(values_short)
    return avg, std,time_short,values_short




def growth_rate(n=None,units='m3dc1',sim=None,filename='C1.h5',
                time_low_lim=500,slurm=True,plottrace=False,pub=False):
    """
    Evaluates kinetic energy growth rate. The growth rate is the mean of the logarithmic derivative of ke.
    The mean is taken over the second half of the simulation time. Returns 
    
    Arguments:

    **n**
    Toroidal mode number. Optional, as n will be read from M3D-C1 output. Should be
    specified, if it cannot be read from output.

    **units**
    System of units for the plot, can be either 'm3dc1' or 'mks'

    **sim**
    fpy simulation object.

    **filename**
    Name or path to C1.h5 file to read.

    **time_low_lim**
    Minimum time in Alfven times that will be considered for averaging.
    Should be chosen such that the noise in the beginning of the simulation is avoided

    **slurm**
    If True, reads the Slurm log file for the M3DC1 run and checks for GS errors and convergence
    
    **plottrace**
    Show and save plots of growth rate in directory.

    **pub**
    If True, format plot for publication.
    """
    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename=filename)
    
    if slurm:
        #try:
        # Read Slurm log files. Should there be multiple log files in the directory, choose the file with largest Slurm Job ID
        C1inputfiles = glob.glob(os.getcwd()+"/C1input")
        if len(C1inputfiles) > 1:
            fpyl.printwarn('WARNING: More than 1 C1input file found. Using the latest one.')
            slurmfiles.sort(key=os.path.getmtime,reverse=True)
        C1inputfile = C1inputfiles[0]
        
        # Read n from Slurm log file
        with open(C1inputfile, 'r') as sf:
            for line in sf:
                if 'ntor   ' in line:
                    ntorline = line.split()
                    n = int(ntorline[2])
                    break
        #except:
        #    n = fpyl.prompt('Not able to detemine ntor. Please enter value for n : ',int)
    else:
        n = readParameter('ntor',sim=sim,listc=False)
        fpyl.printnote('Using n=ntor='+"{:d}".format(n)+' as read from C1.h5 file.')
    
    
    
    gamma, dgamma,time,gamma_trace = avg_time_trace('ke',units,sim=sim,growth=True,renorm=True,start=None,time_low_lim=time_low_lim)
    print(n,gamma,dgamma)
    
    not_noisy = 1
    gamma_set_manu = 0
    flat = 1
    perform_manual_check = False
    
    if abs(dgamma/gamma) > 0.5: #Original value was 0.1
        fpyl.printwarn('WARNING: gamma is not constant for n='+str(n)+'!')
        flat = 0
        
        # Identify local maxima in the growth rate and determine the frequency in which these occur
        maxima_ind = signal.argrelmax(gamma_trace)[0]
        frequency = 1.0/np.diff(maxima_ind)
        #print(gamma_trace[maxima_ind])
        #plt.figure()
        #plt.plot(maxima_ind[:-1],frequency,lw=0,marker='.')
        #plt.plot(time,gamma_trace,lw=1)
        #plt.plot(time[maxima_ind],gamma_trace[maxima_ind],lw=0,marker='.')
        #plt.xlabel('time')
        #plt.ylabel('frequency')
        #plt.title('n='+str(n))
        
        # When frequency is too high, mask it as noise
        noise_mask = np.greater(frequency,0.1)
        
        # Check whether growth rate is noisy for the whole time, partially, or not noisy. In case of noise, do not consider the growth rate.
        if np.all(noise_mask):
            fpyl.printwarn('WARNING: gamma is completely noisy for n='+str(n)+'!')
            not_noisy = 0
        elif np.any(noise_mask):
            fpyl.printwarn('WARNING: gamma is partially noisy for n='+str(n)+'!')
            if (gamma < 0 and (np.all(np.sign(gamma_trace)==-1.0))):
                perform_manual_check = True
            not_noisy = 0
        else:
            # If there is no noise, do additional check ups to determine if the system is stable or unstable
            fpyl.printnote('NOTE: no noise detected for n='+str(n)+'!')
            # Check whether the kinetic energy is overall increasing or decreasing by looking at the evolution of the maxima.
            ke_time,ke,_,_ = get_timetrace('ke',sim=sim,units=units,growth=False,renorm=True)
            
            start_ind = int(np.floor(len(ke)/2))
            start_time = ke_time[start_ind]
            if units.lower() == 'mks':
                time_low_lim = unit_conv(time_low_lim,arr_dim='m3dc1',sim=sim,time=1)
            #print(start_time,time_low_lim)
            if start_time < time_low_lim:
                start_ind = np.argmax(ke_time>time_low_lim)
                #print(start_ind)
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
                
                #plt.figure()
                #for i in range(len((ke_linear))):
                #    plt.plot(time_linear[i],ke_linear[i],lw=1)
                #plt.grid(True)
                #plt.xlabel('time')
                #plt.ylabel('ke')
                #plt.title('n='+str(n))
            else:
                ke_linear = [np.asarray(ke_short)]
            
            # Identify maxima in the kinetic energy
            ke_maxima_ind = signal.argrelmax(ke_short)[0]
            
            # Check if ke is monotonic in between all renormalizations
            if ((all(fpyl.strict_monotonic(kel) for kel in ke_linear)) or(len(ke_maxima_ind)<2)):
                # If it is monotonic (which is a good sign), check the growth rate manually. This is to ensure
                # that oscillations are not too strong, or the value for the growth rate is extracted correctly.
                perform_manual_check = True
            else:
                # Look for local maxima and minima among the maxima to identify a transition from decay to growth or vice versa
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
                
                
                # make different decisions depending on the sign of the growth rate and time trace of ke
                if ((np.sign(gamma) == np.sign(dke_dt[-1])) and (len(max_ke_min)==0 and len(max_ke_max)==0)): # ke and growth rate have same sign, and sign of ke maxima does not change
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
                    #plt.figure()
                    #plt.plot(ke_time_max,max_ke,lw=0,marker='.')
                    #plt.plot(ke_time_max,dke_dt,lw=0,marker='.')
                    #plt.grid(True)
                    #plt.xlabel('time')
                    #plt.ylabel('ke and dke_dt')
                    #plt.title('n='+str(n))
        
        # Prompts for manual check of growth rate
        if perform_manual_check:
            # Show plots of growth rate and kinetic energy for analysis
            double_plot_time_trace_fast('ke',renorm=True,title='',rescale=True,units=units,n=n)
            #plot_time_trace_fast('ke',units=units,filename=filename,growth=False,renorm=True,yscale='linear',rescale=True,save=False)
            #plot_time_trace_fast('ke',units=units,filename=filename,growth=True,renorm=True,yscale='linear',save=False)
            stability = colored('UNSTABLE','yellow') if gamma > 0 else colored('STABLE','green')
            print(stability + ' with gamma='+str(gamma))
            mono_input = fpyl.prompt('Is the calculated growth rate correct? (y/n) : ',['y','n'])
            if mono_input == 'y':
                gamma_set_manu = 1
            else:
                gr_input = fpyl.prompt('Please choose:\n  [1] set a different start time for mean calculation\n  [2] directly enter a value for gamma\n  [3] discard this simulation\n : ',['1','2','3'])
                if gr_input == '1':
                    #ToDo: Test averaging from start_time in various units using various values for time and limits
                    start_time = fpyl.prompt('Please enter a start time for mean calculation : ',float)
                    gamma, dgamma,_,_ = avg_time_trace('ke',units,sim=sim,growth=True,renorm=True,start=start_time,time_low_lim=time_low_lim)
                    gamma_set_manu = 1
                    print('gamma corrected to new avg. value: '+str(gamma))
                elif gr_input == '2':
                    gamma = fpyl.prompt('Please enter a value for the growth rate : ',float)
                    gamma_set_manu = 1
                    print('gamma corrected to value: '+str(gamma))
                elif gr_input == '3':
                    not_noisy = 0
            
            #plt.figure()
            #plt.plot(ke_time_max,max_ke,lw=0,marker='.')
            #plt.plot(ke_time_max,dke_dt,lw=0,marker='.')
            #plt.grid(True)
            #plt.xlabel('time')
            #plt.ylabel('ke')
            #plt.title('n='+str(n))
            
        
        #trace = fft(gamma_trace)
        #trace = fftshift(gamma_trace)
        #trace = np.abs(trace)**2
        #timestep = time[1]-time[0]
        #time = fftshift(fftfreq(len(time), d=timestep))
        #plt.figure()
        #plt.plot(time,trace)
        #plt.title('Power spectrum n='+str(n))
    
    # Show plots of growth rate and kinetic energy if requested
    if (plottrace and (not perform_manual_check)):
        double_plot_time_trace_fast('ke',renorm=True,title='',rescale=True,units=units,n=n,pub=pub,showtitle=True)
    
    #------------------------------------------------
    #Check equilibrium convergence and extract Final GS Error from Slurm log file
    #------------------------------------------------
    if slurm:
        try:
            # Read Slurm log files. Should there be multiple log files in the directory, choose the file with largest Slurm Job ID
            if os.path.isfile('gs_slurm.out'):
                slurmfile = glob.glob(os.getcwd()+"/gs_slurm.out")[0]
                gsoutfile=True
            else:
                gsoutfile=False
                slurmfiles = glob.glob(os.getcwd()+"/slurm*.out")
                if len(slurmfiles) > 1:
                    fpyl.printwarn('WARNING: More than 1 Slurm log file found. Using the latest one.')
                    slurmfiles.sort(key=os.path.getmtime,reverse=True)
                    #print(slurmfiles)
                slurmfile = slurmfiles[0]
            
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
                        maxit = int(gsiterline[3])
                        gserr = float(gsiterline[4])
                    elif 'Final error in GS solution' in line:
                        fgserrline = line.split()
                        finalerrgs = float(fgserrline[-1])
                        if not gsoutfile:
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
                fpyl.printwarn('WARNING: Cannot determine number of GS iterations.')
                gsconvgd = 0
                finalerrgs = 0.0
            return gamma, dgamma, n, flat, not_noisy, gamma_set_manu, gsconvgd, finalerrgs
        except:
            fpyl.printerr('ERROR: Cannot process Slurm log file!')
            return gamma, dgamma, n, flat, not_noisy, gamma_set_manu, 0, 0
    else:
        return gamma, dgamma, n, flat, not_noisy, gamma_set_manu, 0, 0



def scan_n(nmin=1,nmax=20,nstep=1,units='m3dc1',filename='C1.h5',time_low_lim=500,slurm=True,plottrace=False):
    """
    Traverses all subdirectories named nXX in a directory (where XX is the toroidal mode number),
    and reads the growth rate.
    
    Arguments:

    **nmin**
    Lowest mode number to include in the scan

    **nmax**
    Largest mode number to include in the scan

    **nstep**
    Toroidal mode number step size.

    **units**
    System of units for the plot, can be either 'm3dc1' or 'mks'

    **filename**
    Name or path to C1.h5 file to read.

    **time_low_lim**
    Minimum time in Alfven times that will be considered for averaging.
    Should be chosen such that the noise in the beginning of the simulation is avoided

    **slurm**
    If True, reads the Slurm log file for the M3DC1 run and checks for GS errors and convergence
    
    **plottrace**
    Show and save plots of growth rate in directory
    """
    
    n_list = []
    gamma_list = []
    dgamma_list = []
    flat_list = []
    not_noisy_list = []
    gamma_manual_list = []
    
    gsconvgd_list = []
    finalerrgs_list = []
    pblist = []
    ped_loc = []
    print('Calculating growth rates for n='+str(nmin)+' to n='+str(nmax))
    
    nmin_not_found = False
    
    for n in np.arange(nmin,nmax+1,nstep):
        if n==nmin or nmin_not_found:
            path = 'n'+str(n).zfill(2)
        else:
            path = '../n'+str(n).zfill(2)
        try:
            os.chdir(path)
            nmin_not_found = False
        except:
            fpyl.printerr('Cannot find directory '+path)
            nmin_not_found = True
            continue
            #raise Exception('Cannot find directory '+path)
        
        print('----------------------------------------------')
        print('Directory '+os.getcwd().split('/')[-1])

        try:
            gamma, dgamma, n, flat, not_noisy, gamma_set_manu, gsconvgd, finalerrgs = growth_rate(n,units=units,sim=None,filename=filename,time_low_lim=time_low_lim,slurm=slurm,plottrace=plottrace)
        except:
            gamma, dgamma, n, flat, not_noisy, gamma_set_manu, gsconvgd, finalerrgs = (0.0,0.0,n,0,0,0,0,0)
            fpyl.printwarn('WARNING: Did not determine growth rate for n='+str(n))
        
        gamma_list.append(gamma)
        dgamma_list.append(dgamma)
        n_list.append(n)
        flat_list.append(flat)
        not_noisy_list.append(not_noisy)
        gamma_manual_list.append(gamma_set_manu)
        gsconvgd_list.append(gsconvgd)
        finalerrgs_list.append(finalerrgs)
        pblist.append(-100)
        ped_loc.append(-100.)
        
        
        if n==nmax:
            print('----------------------------------------------')
            os.chdir('../')
            
    results = Gamma_data(n_list, gamma_list, dgamma_list, flat_list, not_noisy_list, gamma_manual_list, gsconvgd_list, finalerrgs_list, pblist, ped_loc)
    return results



def create_plot_gamma_n(n_list, gamma_list,norm_dia=False,fignum=None,figsize=None,lw=1,c=None,ls=None,marker=None,ms=36,lbl=None,units='m3dc1',xtick_min=None, xtick_step=1,legfs=None,leglblspace=None,leghandlen=None,title=None,export=False,txtname=None,pub=False):
    
    # Set font sizes and plot style parameters
    if pub:
        axlblfs = 20
        titlefs = 18
        ticklblfs = 16
        linew = 2
        inplttxtfs = 20
        if legfs is None:
            legfs = 16
    else:
        axlblfs = None
        titlefs = None
        ticklblfs = None
        linew = 1
        inplttxtfs = 16
    
    plt.figure(num=fignum,figsize=figsize)
    if marker in ['s','d','D']:
        ms = int(ms*0.75)
    elif marker == 'v':
        ms = int(ms*0.7)
        print(ms)
    if ls==':' and lw is not None:
        lw = lw+1
    temp = plt.plot(n_list,gamma_list,lw=lw,c=c,ls=ls,marker=marker,ms=ms,label=lbl)
    plt.grid(True)
    plt.xlabel('n',fontsize=axlblfs)
    if norm_dia:
        ylbl = r'$\gamma/(\omega_{*i}/2)$'
    else:
        if units.lower()=='m3dc1':
            ylbl = r'$\gamma/\omega_A$'
        else:
            ylbl = r'$\gamma$ $[s^{-1}]$'
    plt.ylabel(ylbl,fontsize=axlblfs)
    ax = plt.gca()
    
    ax.set_xlim([np.amin(n_list)-0.5,np.amax(n_list)+0.5])
    
    #if xtick_step == 1:
    nmax = np.amax(n_list)+1
    #else:
        #if np.amin(n_list) + xtick_step <= np.amax(n_list)+1
        #nmax = np.amax(n_list)+1
    if xtick_min is None: xtick_min = np.amin(n_list)
    ax.xaxis.set_ticks(np.arange(xtick_min, nmax, xtick_step))
    if isinstance(xtick_step,(int,float)) and xtick_step>1:
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    #ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))
    plt.tick_params(axis='both', which='major', labelsize=ticklblfs)
    if title is not None:
        plt.title(title,fontsize=titlefs)
    if lbl is not None:
        plt.legend(loc=0,fontsize=legfs,labelspacing=leglblspace,handlelength=leghandlen)
    else:
        
        artistnorot = plt.Line2D((0,1),(0,0), color='k', marker='v', linestyle='-')
        artistrot = plt.Line2D((0,1),(0,0), color='k', marker='s', linestyle=':')
        artisteta01 = plt.Line2D((0,1),(0,0), color='C0',ls='-')
        artisteta1 = plt.Line2D((0,1),(0,0), color='C1',ls='-')
        artisteta2 = plt.Line2D((0,1),(0,0), color='C2',ls='-')
        artisteta10 = plt.Line2D((0,1),(0,0), color='C3',ls='-')

        #Create legend from custom artist/label lists
        #ax.legend([artistnorot,artistrot,artisteta01,artisteta1,artisteta2,artisteta10],[r'$\omega_0 = 0$', r'$\omega_0 =\omega_i$',r'$\eta \times 0.1$',r'$\eta \times 1$',r'$\eta \times 2$',r'$\eta \times 10$'])
    plt.tight_layout()
    
    if export:
        data_points = temp[0].get_data()
        print(np.transpose(data_points))
        np.savetxt(txtname,data_points,delimiter='   ')
        print(txtname)
    return



def plot_gamma_n(nmin=1,nmax=20,nstep=1,norm_dia=False,units='m3dc1',fignum=None,figsize=None,c=None,lw=None,ls=None,mark='.',plot_crosses=True,xtick_min=None,xtick_step=1,lbl=None,slurm=True,plottrace=False,legfs=None,leglblspace=None,leghandlen=None,ylimits=None,title=None,export=False,txtname=None,no_prompt=False,pub=False):
    # Plot gamma as a function of n
    # Identify simulations where the growth rate was not calculated reliably. These are highlighted in the plot.
    #print(os.getcwd())
    files = glob.glob("growth_rates*.dat")
    
    if len(files)>0:
        if len(files)>1:
            files.sort(key=os.path.getmtime,reverse=True)
            print('More than 1 file found. Opening newest file: '+files[0])
        f = files[0]
        results = Gamma_file(f)
    else:
        if not no_prompt:
            eval_input = fpyl.prompt('No results found. Do you want to determine the growth rates? (y/n) : ',['y','n'])
            if eval_input == 'y':
                results = scan_n(nmin,nmax,nstep,units=units,slurm=True,plottrace=plottrace)
                if norm_dia:
                    norm_dia = False
                    fpyl.printwarn('WARNING: Diamagnetic frequency has not been calculated. Setting norm_dia=False.')
            else:
                return
        else:
            return
    
    # Identify simulations where the growth rate was not calculated reliably. These are highlighted in the plot.
    n_all = []
    gamma_all = []
    n_bad = []
    gamma_bad = []
    for n in results.n_list:
        #n_ind = n-1
        n_ind = fpyl.get_ind_at_val(results.n_list,n,unique=True)
        
        n_all.append(results.n_list[n_ind])
        if norm_dia:
            #gamma_all.append(results.gamma_list[n_ind]/(results.n_list[n_ind]*results.omegsti_max/4))
            gamma_all.append(results.gamma_dia_list[n_ind])
        else:
            gamma_all.append(results.gamma_list[n_ind])
        
        if (results.flat_list[n_ind] != 1 and results.gamma_manual_list[n_ind] != 1 and results.not_noisy_list[n_ind] != 1):
            n_bad.append(results.n_list[n_ind])
            if norm_dia:
                #gamma_bad.append(results.gamma_list[n_ind]/(results.n_list[n_ind]*results.omegsti_max/4))
                gamma_bad.append(results.gamma_dia_list[n_ind])
            else:
                gamma_bad.append(results.gamma_list[n_ind])
    
    if len(n_bad)>0 and plot_crosses:
        if fignum is None:
            fignum = plt.gcf().number #Current figure number
        create_plot_gamma_n(n_bad, gamma_bad, norm_dia=norm_dia, fignum=fignum, figsize=figsize, lw=0, marker='x', ms=10, c='r', units=units,xtick_min=xtick_min,xtick_step=xtick_step,legfs=None,leglblspace=None,title=title)
    
    create_plot_gamma_n(n_all, gamma_all, norm_dia=norm_dia, fignum=fignum, figsize=figsize,c=c, lw=lw, ls=ls, marker=mark, ms=10, lbl=lbl, units=units,xtick_min=xtick_min,xtick_step=xtick_step,legfs=legfs,leglblspace=leglblspace,leghandlen=leghandlen,title=title,export=export,txtname=txtname,pub=pub)
    
    
    cfig = plt.gcf()
    ax = cfig.gca()
    if ylimits != None:
        ax.set_ylim(ylimits)
    
    return




def compare_gamma_n(dirs,nmin=1,nmax=20,nstep=1,norm_dia=False,units='m3dc1',labels=None,plot_crosses=True,col=None,lwid=None,lsty=None,markers=None,xtick_min=None,xtick_step=1,fignum=None,figsize=None,legfs=None,leglblspace=None,leghandlen=None,ylimits=None,title=None,export=False,no_prompt=False,quiet=False,pub=False,L0=None,R0=None,B0=None):
    if isinstance(labels, (tuple, list)):
        if len(dirs)!=len(labels):
            fpyl.printerr('ERROR: Number of directories not equal to number of labels.')
            return
    pwd = os.getcwd()
    for i,d in enumerate(dirs):
        if os.path.isdir(d):
            os.chdir(d)
        else:
            if not quiet:
                fpyl.printerr('ERROR: Directory ' + d + ' does not exist!')
            continue
        
        if isinstance(labels, (tuple, list)):
            lbl = labels[i]
        else:
            lbl = d
        if isinstance(col, (tuple, list)):
            c = col[i]
        else:
            c = None
        if isinstance(lwid, (tuple, list)):
            lw = lwid[i]
        else:
            lw = None
        if isinstance(lsty, (tuple, list)):
            ls = lsty[i]
        else:
            ls = None
        if isinstance(markers, (tuple, list)):
            mark = markers[i]
        else:
            mark = '.'
        plot_gamma_n(nmin,nmax,nstep,norm_dia=norm_dia,units=units,fignum=fignum,figsize=figsize,xtick_min=xtick_min,xtick_step=xtick_step,c=c,lw=lw,ls=ls,mark=mark,plot_crosses=plot_crosses,lbl=lbl,slurm=True,plottrace=False,legfs=legfs,leglblspace=leglblspace,leghandlen=leghandlen,ylimits=ylimits,title=title,export=export,txtname='gamma_'+d.replace('/','')+'.txt',no_prompt=no_prompt,pub=pub)
        os.chdir(pwd)
    
    return



def write_gamma_n(results,ped_param, ipres, psin_ped_top,ped_structure=None,units='m3dc1',fix=False):
    #if (nmin >= 0 and nmax >=0 and nmax > nmin):
    #    #n_list, gamma_list, dgamma_list, not_noisy_list, gsconvgd, finalerrgs = scan_n(nmin,nmax,slurm=True)
    #    results = scan_n(nmin,nmax,units=units,slurm=True,plottrace=True)
    #elif (nmin < 0 and nmax < 0) and (results is None):
    #    raise Exception('nmin and nmax must be >=0 with nmax > nmin.')
    
    pwd = os.getcwd()
    pathdirs = pwd.split('/')
    vpnum = pathdirs[-2]
    if vpnum == 'convergence_study':
        vpnum = pathdirs[-3]
    if '_w' in vpnum:
        temp = vpnum.split('_w')
        vpnum = temp[0]
        width = temp[1]
    else:
        width = -1
    
    simdir = pathdirs[-1]
    
    eta = re.search(r'eta_x(\d*[.]?\d*)', simdir).group(1)
    
    flmodel = simdir[0:2]
    if flmodel not in ['1f', '2f']:
        fpyl.printerr('ERROR: Fluid model not recognized!')
    
    if 'eqrotnc' in simdir:
        rotation=-1
    elif 'eqrot_' in simdir:
        rotation=1
    elif 'norot' in simdir:
        rotation=0
    else:
        fpyl.printwarn('WARNING: Equilibrium rotation model unknown!')
    
    if 'B' in simdir:
        bscale = re.search('B(.*)_', simdir).group(1)
    else:
        bscale = '1.0'
    
    
    alpha_max,j_max,jelite_N,omegsti_max=ped_param
    
    #print(results.n_list)
    
    if float(width) > 0:
        outfile = 'growth_rates_'+vpnum+'_' + width + '_' +simdir+'.dat'
    else:
        outfile = 'growth_rates_'+vpnum+'_'+simdir+'.dat'
    if fix:
        outfile = outfile + '_fix'
    if units == 'mks':
        outfile = 'mks_'+outfile
    
    try:
        with open(outfile, 'w') as f:
            f.write(vpnum+'    '+"{:f}".format(float(eta))+'    '+bscale+'    '+str(rotation)+'    '+flmodel+'    '+str(ipres)+'\n')
            #f.write(str(jpdata[1])+'    '+str(jpdata[0])+'\n') #write j and pprime
            f.write(str(alpha_max)+'    '+str(j_max)+'    '+str(jelite_N)+'    '+str(omegsti_max)+'    '+str(psin_ped_top)+'    '+str(width)+'\n') #write alpha, j_max, j_elite and omega_*
            
            if ped_structure is not None:
                f.write(str(ped_structure[0])+'    '+str(ped_structure[1])+'\n') #write pedestal height and width as determined by pedestal finder
                
            f.write('    n      gamma         sig_gamma     flat     smooth   manu     conv  Fin. GS Err   PB\n')
            
            gam_format = "{0:.8f}" if units == 'm3dc1' else "{0:.8e}"
            for i in range(len(results.gamma_list)):
                wstr = '    ' + "{:d}".format(results.n_list[i]).ljust(2,' ')
                wstr = wstr + gam_format.format(results.gamma_list[i]).rjust(15,' ')
                wstr = wstr + "{0:.8f}".format(results.dgamma_list[i]).rjust(14,' ')
                wstr = wstr + '    ' + "{:d}".format(results.flat_list[i]).ljust(5,' ')
                wstr = wstr + '    ' + "{:d}".format(results.not_noisy_list[i]).ljust(5,' ')
                wstr = wstr + '    ' + "{:d}".format(results.gamma_manual_list[i]).ljust(5,' ')
                wstr = wstr + '    ' + "{:d}".format(results.gsconvgd[i]).ljust(2,' ')
                wstr = wstr + '    ' + "{0:.8f}".format(results.finalerrgs[i])
                wstr = wstr + '    ' + "{:d}".format(results.pblist[i])
                wstr = wstr + '    ' + "{0:.8f}".format(results.ped_loc[i])+'\n'
                f.write(wstr)
        print("Growth rate data written to '"+str(outfile)+"'.")
    except Exception as e:
        print(e)
    return



def omegastari(sim=None,filename='C1.h5',time=None,units='mks',points=400,n=1,pion=False,fcoords='pest',makeplot=False):
    #Calculates ion diamagnetic frequency
    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename=filename)
    psin,psi = flux_average('psi',coord='scalar',sim=sim, fcoords=fcoords, linear=False, deriv=0, points=points, phit=0.0, filename=filename, time=time, psin_range=None, units='mks')
    psi = psi*2*math.pi # Because psi is the poloidal flux per radiant as in B = grad(psi) x grad(phi) + B_phi
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



def get_ped_param(sim,filename='C1.h5',time=None,points=400,pion=False,fcoords='pest',psin_ped_top=0.86,psin_var_j=0.85,use_max_j=False,device=None):
    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename=filename)
    #psinedgelim = 0.86 #Min value of psin that is considered edge region
    print('get_ped_param time:'+str(time))
    # Determine pedestal alpha
    psi_a,alpha = flux_average('alpha', coord='scalar', sim=sim, time=time, fcoords=fcoords, points=points, units='m3dc1')
    psinedge = fpyl.find_nearest(psi_a,psin_ped_top)
    psinedge_ind = fpyl.get_ind_at_val(psi_a,psinedge)
    alpha_max = np.amax(alpha[psinedge_ind:])
    #Check if it's really a local maximum inside the pedestal:
    alpha_max_ind = fpyl.get_ind_at_val(alpha,alpha_max)
    #print(alpha_max_ind,psi_a[alpha_max_ind])
    if (alpha_max_ind < len(alpha)-1 and (alpha_max >= alpha[alpha_max_ind-1]) and (alpha_max >= alpha[alpha_max_ind+1])) or (alpha_max_ind==len(alpha)-1): #Check if the maximum is a relative maximum or occurs at the lcfs
        #print('PEDESTAL ALPHA IS ABSOLUTE MAXIMUM')
        print('Pedestal alpha = '+str(alpha_max))
    else:
        alpha_short = alpha[psinedge_ind:]
        psin_a_short = psi_a[psinedge_ind:]
        alpha_rel_max = signal.argrelmax(alpha_short)
        if len(alpha_rel_max)==0:
            fpyl.printwarn('WARNING: maximum of alpha not found!')
            return None,None,None,None
        maxima = np.take(alpha_short,alpha_rel_max)
        maxima_pos = np.take(psin_a_short,alpha_rel_max)
        alpha_max = np.amax(maxima)
        #print(maxima,maxima_pos)
        #print('PEDESTAL ALPHA IS RELATIVE MAXIMUM')
        print('Pedestal alpha = '+str(alpha_max))
    
    # Determine pedestal average toroidal current density
    psi_j,j = flux_average('j', coord='phi', sim=sim, time=time, fcoords=fcoords, points=points, units='m3dc1')
    psinedge = fpyl.find_nearest(psi_j,psin_ped_top)
    psinedge_ind = fpyl.get_ind_at_val(psi_j,psinedge)
    j_max = np.amax(j[psinedge_ind:])
    #Check if it's really a local maximum inside the pedestal:
    j_max_ind = fpyl.get_ind_at_val(j,j_max)
    if j_max_ind < len(j)-1 and (j_max >= j[j_max_ind-1]) and (j_max >= j[j_max_ind+1]):
        print('Pedestal j_max = '+str(j_max))
    else:
        j_max = -1
        fpyl.printwarn('WARNING: j_max has no maximum inside the pedestal!')
        #return None,None,None,None
    
    #Calculate pedestal parallel current density as in ELITE:
    jav = flux_average('jav', coord='scalar', sim=sim, time=time, fcoords=fcoords, points=points, units='mks')[1]
    psi_j,jelite = flux_average('jelite', coord='scalar', sim=sim, time=time, fcoords=fcoords, points=points, units='mks',device=device)
    psinedge = fpyl.find_nearest(psi_j,psin_ped_top)
    psinedge_ind = fpyl.get_ind_at_val(psi_j,psinedge)
    
    if use_max_j:
        #If maximum jelite inside the pedestal region is to be used:
        jelite_max = np.amax(jelite[psinedge_ind:])
    else:
        #Check if jelite peaks inside the pedestal. If it does, take peak value as jelite. If it does
        #not peak, then use jelite(psin=psin_var_j).
        #Calculate absolute maximum in pedestal region:
        jelite_max = np.amax(jelite[psinedge_ind:])
        jelite_max_ind = fpyl.get_ind_at_val(jelite,jelite_max)
        #print(jelite_max_ind, psi_j[jelite_max_ind], jelite_max,jelite[jelite_max_ind-1],jelite[jelite_max_ind+1])
        #If local maximum
        if jelite_max_ind < len(jelite)-1 and (jelite_max >= jelite[jelite_max_ind-1]) and (jelite_max >= jelite[jelite_max_ind+1]):
            jelite_rel_max = jelite_max
            fpyl.printnote('Used MAXIMUM for jelite')
        else:
            jelite_rel_max = 0.0
        #Calculate average value in pedestal region
        jelite_avg = np.average(jelite[psinedge_ind:])
        #print(jelite_rel_max,jelite_avg)
        j_threshold=1.2
        #If jelite_rel_max > j_threshold*jelite_avg then we consider a peak inside the pedestal.
        #Otherwise, use jelite(psin=psin_var_j):
        if not jelite_rel_max > j_threshold*jelite_avg:
            psinj = fpyl.find_nearest(psi_j,psin_var_j)
            psinj_ind = fpyl.get_ind_at_val(psi_j,psinj)
            jelite_max = jelite[psinj_ind]
            fpyl.printwarn('Used psin_var_j for jelite')
    jelite_sep = jelite[-1]
    jelite_N = (jelite_max+jelite_sep)/(2.0*jav[-1])
    print('Current density jelite = '+str(jelite_N))
    
    
    # Determine diamagnetic frequency
    psi_o,omegsti = omegastari(sim=sim,time=time,units='m3dc1',n=1,points=points,pion=pion,fcoords=fcoords,makeplot=False)
    psinedge = fpyl.find_nearest(psi_o,psin_ped_top)
    psinedge_ind = fpyl.get_ind_at_val(psi_o,psinedge)
    omegsti_max = np.amax(omegsti[psinedge_ind:]) #ToDo: Check if it's really a local maximum
    print('Ion diamagnetic frequency = '+str(omegsti_max))
    ped_param = [alpha_max,j_max,jelite_N,omegsti_max]
    return ped_param







def create_plot_time_trace_fast(time,scalar,trace,units='mks',millisec=False,sim=None,filename='C1.h5',
                                growth=False,diff=False,yscale='linear',rescale=False,save=False,in_plot_txt=None,
                                time_marks=[],ts_marks=[],ts_marks_all=False,savedir=None,title=False,pub=False,
                                ylbl=None,show_legend=False,leglbl=None,fignum=None,figsize=None,skip_n0=False,
                                export=False,txtname=None):

    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename)
        
    # If one array has new data but the other one doesn't 
    # plot only previous data
    if scalar.shape != time.shape:
        ymax   = np.amax(scalar.shape)
        tmax   = np.amax(time.shape)
        maxidx = np.amin([ymax,tmax])
        time   = time[0:maxidx]
        scalar = scalar[0:maxidx]
    
    ntor = readParameter('ntor',sim=sim)
    
    plt.figure(num=fignum,figsize=figsize)
    
    # Set font sizes and plot style parameters
    if pub:
        axlblfs = 20
        titlefs = 18
        ticklblfs = 16
        linew = 2
        legfs = 14
        inplttxtfs = 20
    else:
        axlblfs = None
        titlefs = None
        ticklblfs = None
        linew = 1
        legfs = None
        inplttxtfs = 16
    
    if units=='mks':
        #time = unit_conv(time,filename=filename,time=1)
        if millisec:
            time = time*1000
            plt.xlabel(r'time $(ms)$',fontsize=axlblfs)
        else:
            plt.xlabel(r'time $(s)$',fontsize=axlblfs)
    else:
        plt.xlabel(r'time $(\tau_A)$',fontsize=axlblfs)
    
    if diff:
        plt.ylabel('d '+ylbl+' / d t',fontsize=axlblfs)
    else:
        if trace == 'ke':
            if units=='mks':
                #scalar = unit_conv(scalar, arr_dim='M3DC1', filename=filename, energy=1)
                if growth:
                    plt.ylabel(r'$\gamma$ $[s^{-1}]$')
                else:
                    plt.ylabel(r'Kinetic energy $[J]$')
            elif units.lower()=='m3dc1':
                if growth:
                    plt.ylabel(r'$\gamma/\omega_A$')
                else:
                    plt.ylabel(r'Kinetic energy (M3DC1 units)')
        else:
            plt.ylabel(ylbl,fontsize=axlblfs)
    
    if export:
        plot_data = np.column_stack([time])
    
    if len(scalar.shape)>1:
        leglabels = ["n={:2}".format(n) for n in range(scalar.shape[1])]
        for i in range(scalar.shape[1]):
            if not (skip_n0 and i==0):
                plt.plot(time,scalar[:,i],lw=linew,label=leglabels[i])
        
        ncol = 2 if scalar.shape[1] > 5 else 1
        if show_legend:
            plt.legend(loc=0, ncol=ncol,fontsize=legfs)
        else:
            labelLines(plt.gca().get_lines(), zorder=2.5)
    else:
        plt.plot(time,scalar,lw=linew,label=leglbl)
        if show_legend and leglbl is not None:
            plt.legend(loc=0,fontsize=legfs)
            
    if export:
        plot_data = np.column_stack([plot_data, scalar])
    plt.grid(True)
    #Determine y-axis limits
    if rescale:
        if np.amax(scalar[1:]) < scalar[0]:
            start_time=250
            if units=='mks':
                start_time = unit_conv(start_time,arr_dim='m3dc1',sim=sim,time=1)
            start_ind = int(fpyl.find_nearest(time,start_time))
            top_lim=1.1*np.amax(scalar[start_ind:])
            plt.ylim([0,top_lim])
    
    #Plot vertical lines to mark certain points in time specified by time_marks.
    if len(time_marks)>0:
        for t in time_marks:
            plt.axvline(x=t,c='m')
    
    if ts_marks_all:
        slices = glob.glob('time*.h5')
        slices.sort(key=os.path.getmtime)
        temp = [s.replace('time_', '') for s in slices]
        temp = [s.replace('.h5', '') for s in temp]
        ts_marks = np.array(temp,dtype=int)
    
    if len(ts_marks)>0:
        for ts in ts_marks:
            ts_marks_time = get_time_of_slice(ts,filename=filename,units=units,millisec=millisec)
            plt.axvline(x=ts_marks_time,c='m')
    
    ax = plt.gca()
    if in_plot_txt is not None:
        plt.text(0.03, 0.95,in_plot_txt, ha='left', va='top', transform=ax.transAxes,fontsize=inplttxtfs)
    
    plt.yscale(yscale)
    ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
    if pub:
        yaxis_magn = ax.yaxis.get_offset_text()
        yaxis_magn.set_size(ticklblfs)
    if yscale == 'linear':
        plt.ticklabel_format( axis='y', style='sci',useOffset=False)
    if title:
        plt.title('n='+str(ntor),fontsize=titlefs)
    plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    plt.show()
    
    if save:
        tracestr = trace
        if growth:
            tracestr = tracestr + '-growth'
        if savedir is not None:
            tracestr = savedir + tracestr
        plt.savefig(tracestr+'_n'+"{:d}".format(ntor)+'.pdf', format='pdf',bbox_inches='tight')
        
    if export:
        print(plot_data)
        np.savetxt(txtname,plot_data,delimiter='   ')
    return



def double_plot_time_trace_fast(trace,sim=None,filename='C1.h5',renorm=False,rescale=False,
                                units='m3dc1',title=None,n=None,pub=False,showtitle=True):
    """
    Plots a trace and its growth rate in two subplots side by side
    
    Arguments:

    **trace**
    String containing the trace to be plotted

    **filename**
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

    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename)
    
    time,scalar,_,_ = get_timetrace(trace,sim=sim,units=units,growth=False,renorm=renorm,quiet=True)
    time_growth,scalar_growth,_,_ = get_timetrace(trace,sim=sim,units=units,growth=True,renorm=renorm,quiet=True)
    
    # Set font sizes and plot style parameters
    if pub:
        axlblfs = 20
        titlefs = 20
        ticklblfs = 20
        lcfslw = 2
    else:
        axlblfs = None
        titlefs = 20
        ticklblfs = None
        lcfslw = 1
    
    if units=='mks':
        xlbl = r'time $/s$'
        ylbl_left = r'Kinetic energy $/J$'
        ylbl_right = r'$\gamma$ $/s^{-1}$'
    else:
        xlbl = r'time $/\tau_A$'
        ylbl_left = r'Kinetic energy (M3DC1 units)'
        ylbl_right = r'$\gamma/\omega_A$'
    
    #Determine y-axis limits
    if rescale:
        if np.amax(scalar[1:]) < scalar[0]:
            start_time=250
            if units=='mks':
                start_time = unit_conv(start_time,arr_dim='m3dc1',sim=sim,time=1)
            start_ind = int(fpyl.find_nearest(time,start_time))
            top_lim=1.1*np.amax(scalar[start_ind:])
        else:
            top_lim=None
    
    fig = plt.figure(constrained_layout=True,figsize=(10,5))
    spec2 = gridspec.GridSpec(ncols=2, nrows=1, figure=fig)
    f2_ax1 = fig.add_subplot(spec2[0, 0])
    f2_ax2 = fig.add_subplot(spec2[0, 1])
    
    f2_ax1.plot(time,scalar,lw=lcfslw)
    f2_ax1.set_xlabel(xlbl,fontsize=axlblfs)
    f2_ax1.set_ylabel(ylbl_left,fontsize=axlblfs)
    if rescale:
        f2_ax1.set_ylim([0,top_lim])
    f2_ax1.grid()
    f2_ax1.tick_params(axis='both', which='major', labelsize=ticklblfs)
    
    f2_ax2.plot(time_growth,scalar_growth,lw=lcfslw)
    f2_ax2.set_xlabel(xlbl,fontsize=axlblfs)
    f2_ax2.set_ylabel(ylbl_right,fontsize=axlblfs)
    f2_ax2.grid()
    f2_ax2.tick_params(axis='both', which='major', labelsize=ticklblfs)
    if showtitle:
        if n is None:
            ntor = readParameter('ntor',sim=sim)
        else:
            ntor=n
        fig.suptitle('n='+str(ntor), size=titlefs)
    
    return


def plot_time_trace_fast(trace,units='mks',millisec=False,sim=None,filename='C1.h5',diff=False,
                         growth=False,renorm=False,yscale='linear',unitlabel=None,fac=1,
                         show_legend=False,leglbl=None,in_plot_txt=None,time_marks=[],ts_marks=[],
                         ts_marks_all=False,rescale=False,save=False,savedir=None,pub=False,
                         fignum=None,figsize=None,drop_time_steps=None,skip_n0=False,export=False,txtname=None):
    """
    Plots a scalar quantity vs. time. All available
    time traces can be found in the M3D-C1 documentation.
    
    Arguments:

    **trace**
    String containing the trace to be plotted

    **units**
    The units in which the trace will be plotted

    **millisec**
    True/False. If True and units='mks' plot will be in terms of milliseconds, instead of seconds.

    **sim**
    fpy simulation object.

    **filename**
    Contains the path to the C1.h5 file.

    **diff**
    True / False. Plot derivative of scalar.

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

    **unitlabel**
    Specify custom unitlabel. If not specified, default label will be used.

    **fac**
    Factor to apply to field when using mks units. fac=1.0E-3 converts to kilo..., fac=1.0E-6 to Mega...

    **show_legend**
    If True, show plot legend.

    **leglbl**
    Legend label for plot curve.

    **in_plot_txt**
    Overlay text to be shown inside of plot window.

    **time_marks**
    List of times that will be marked in terms of a vertical line in plot.

    **ts_marks**
    List of time slice numbers. A vertical line will be added for the time corresponding
    to each time slice listed.

    **ts_marks_all**
    True / False. If True, indicate times of all time slices.

    **rescale**
    Rescale y-axis such that noise in the beginning of the simulation is not considered for axis limits,
    i.e. plot is zoomed in to show relevant values.

    **save**
    If True, save plot to file

    **savedir**
    Directory where plot will be saved

    **pub**
    If True, format figure for publication (larger labels and thicker lines)

    **fignum**
    Figure number.

    **figsize**
    If True, format figure for publication (larger labels and thicker lines)

    **drop_time_steps**
    Number of time steps at the end of time trace to drop from plot.

    **skip_n0**
    When plotting energy spectrum, do not plot energy for n=0 mode.

    """
    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename)
    time,y_axis,label, unitlabel = get_timetrace(trace,sim=sim,units=units,growth=growth,diff=diff,renorm=renorm,unitlabel=unitlabel,fac=fac)
    if drop_time_steps is not None:
        time = time[:-drop_time_steps]
        y_axis = y_axis[:-drop_time_steps]
    if unitlabel is None or unitlabel == '':
        ylbl = label
    else:
        ylbl = label+' ('+unitlabel+')'
    
    filename = sim.filename
    create_plot_time_trace_fast(time,y_axis,trace,units=units,millisec=millisec,sim=sim,filename=filename,growth=growth,diff=diff,show_legend=show_legend,
                                leglbl=leglbl,yscale=yscale,rescale=rescale,save=save,savedir=savedir,pub=pub,in_plot_txt=in_plot_txt,
                                time_marks=time_marks,ts_marks=ts_marks,ts_marks_all=ts_marks_all,ylbl=ylbl,fignum=fignum,figsize=figsize,skip_n0=skip_n0,
                                export=export,txtname=txtname)
    return



def integrate_time_trace(trace,nts=None,method=None,units='mks',sim=None,
                         filename='C1.h5',growth=False,renorm=False,yscale='linear',
                         rescale=False,save=False,savedir=None,makeplot=True,fignum=None,leglbl=None,show_legend=False,pub=False):
    """
    Integrate OR cumulatively integrate time trace from zero to time step specified by nts.
    
    Arguments:

    **trace**
    String containing the trace to be plotted

    **method**
    method = trapz: integrate using composite trapezoidal rule
    method = cumtrapz: cumulatively integrate using composite trapezoidal rule

    **nts**
    Upper boundary for integration in terms of time step. If None, integration will be over whole simulation.

    **sim**
    simulation sim_data objects. If none is provided, plot_field will read a file and create
    an object.

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

    **rescale**
    Rescale y-axis such that noise in the beginning of the simulation is not considered for axis limits,
    i.e. plot is zoomed in to show relevant values.

    **save**
    If True, save plot to file

    **savedir**
    Directory where plot will be saved

    **makeplot**
    If True, show plot of integrand.

    **pub**
    If True, format figure for publication (larger labels and thicker lines)
    """
    if method is None:
        fpyl.printerr('ERROR: Please specify method! Possible choices are trapz and cumtrapz.')
    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename)
    
    #scalar = sim.get_time_trace(trace)
    scalar,_,_ = get_timetrace(trace,sim=sim,units=units,growth=growth,renorm=renorm,returnas='time_trace')

    if fignum is None:
        fignum = plt.gcf().number + 1
    if makeplot:
        create_plot_time_trace_fast(scalar.time,scalar.values,trace,units=units,sim=sim,growth=growth,
                                    yscale=yscale,rescale=rescale,save=save,savedir=savedir,pub=pub,
                                    leglbl=leglbl,show_legend=show_legend,fignum=fignum+100)
    
    if nts is None:
        temp = -1
    else:
        temp = nts
    print('Integrating from t=0 to t='+str(scalar.time[temp]))
    if method == 'trapz':
        trace_integrated = trapz(scalar.values[:nts],scalar.time[:nts])
        print('Integral of '+trace+' in '+units+' units: '+str(trace_integrated))
    elif method == 'cumtrapz':
        scalar_int = scalar.cum_int(nts)
        if makeplot:
            create_plot_time_trace_fast(scalar_int.time,scalar_int.values,trace,units=units,sim=sim,growth=growth,
                                    yscale=yscale,rescale=rescale,save=save,savedir=savedir,pub=pub,
                                    leglbl=leglbl,show_legend=show_legend,fignum=fignum)
        trace_integrated = scalar_int
    
    return trace_integrated


