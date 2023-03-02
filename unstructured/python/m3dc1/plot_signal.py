#!/usr/bin/env python3
#
# plot_signal: plots diagnostics signal as a function of time or its power spectrum
#
# Coded on 12/18/2019 by:
# Andreas Kleiner:    akleiner@pppl.gov

import fpy
import numpy as np
from scipy.fftpack import fft, fftshift, fftfreq
import matplotlib.pyplot as plt
from matplotlib import rc
import m3dc1.fpylib as fpyl
from m3dc1.eval_field import eval_field
from m3dc1.read_h5 import readC1File
from m3dc1.unit_conv import unit_conv
from m3dc1.compensate_renorm import compensate_renorm
from m3dc1.time_trace_fast import growth_rate
rc('text', usetex=False)



def plot_signal(signal='mag_probes', filename='C1.h5', sim=None, renorm=False,
                scale=False, deriv=False, pspec=False, units='mks', pts_per_probe=1,
                time_low_lim=500, pub=False):
    """
    Plots diagnostics signal
    
    Arguments:

    **signal**
    Type of diagnostics to plot the signal of, e.g. 'mag_probes' or 'flux_loops'

    **filename**
    File name which will be read, i.e. "../C1.h5".

    **renorm**
    Remove spikes due to renormalization in linear runs that happens when
    the kinetic energy gets too large.

    **scale**
    Scale signal by the signal of a perfectly linear growth

    **deriv**
    Calculate time derivative of signal

    **pspec**
    Calculate Fourier spectrum of signal time trace.

    **units**
    Units in which the result will be calculated

    **time_low_lim**
    lower limit for starting time in terms of Alfven times. If start < time_low_lim,
    start time is moved to the earliest time larger than time_low_lim

    **pub**
    If True, format figure for publication (larger labels and thicker lines)
    """
    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename)
    
    time = np.asarray(sim._all_traces['time'])
    trace = sim.get_signal('mag_probes').sigvalues
    #print(trace.shape)
    trace = np.transpose(trace)
    #print(trace.shape)
    nsamples = trace.shape[0]
    #print(nsamples)
    if float.is_integer(nsamples/pts_per_probe):
        nprobes = int(nsamples/pts_per_probe)
    else:
        fpyl.printerr('ERROR: Number of magnetic probes is not an integer. Check value of pts_per_probe!')
        return

    #print(trace.shape)
    
    #print(trace2.shape)


    if units.lower()=='mks':
        time = unit_conv(time, filename=filename, time=1)
    
    
    if pts_per_probe > 1:
        print(trace.shape)
        trace = np.mean(np.asarray(np.split(trace,nprobes,axis=0)),axis=1)
        print(trace.shape)
    
    
    if renorm:
        for i in range(nprobes):
            trace[i,:] = compensate_renorm(trace[i,:])
    
    # ToDo: Add this calculation
    if scale:
        gamma = 0.5*growth_rate(n=None,units=units,sim=sim,time_low_lim=time_low_lim,slurm=False,plottrace=False)[0]
        print(gamma)
        for i in range(nprobes):
            tg = time*gamma
            tg = tg.astype(np.float64)
            trace[i,:] = trace[i,:]/np.exp(tg)/trace[-1,i]
    
    if deriv:
        for i in range(nprobes):
            trace[i,:] = fpyl.deriv(trace[i,:],time)
    
    
    
    # Set font sizes and plot style parameters
    if pub:
        axlblfs = 20
        #titlefs = 18
        ticklblfs = 18
        legfs = 12
        linew = 2
    else:
        axlblfs = None
        #titlefs = None
        ticklblfs = None
        legfs = None
        linew = 1
    
    plt.figure()
    tracefft = np.zeros_like(trace)
    if pspec:
        timestep = time[1]-time[0]
        time = fftshift(fftfreq(len(time), d=timestep))
        for i in range(nprobes):
            trace[i,:] = trace[i,:]/np.amax(np.abs(trace[i,:]))#use max of absolute value, because freq can be negative and thus max would be 0
            temp = 1.0/len(trace[i,:])*fft(trace[i,:])
            if i==1:
                print(trace[i,:])
            tracefft[i,:] = np.abs(fftshift(temp*np.amax(trace[i,:])))**2
            if i==1:
                print(tracefft[i,:])
            peak_ind = np.argmax(tracefft[i,:])
            if units.lower()=='mks':
                print('Probe {i}: Signal peaks at a frequency of {freq:.2f} Hz'.format(i = i, freq = time[peak_ind]))
            elif units.lower()=='m3dc1':
                print('Probe {i}: Signal peaks at a frequency of {freq:.6f} / Alfven time'.format(i = i, freq = time[peak_ind]))
            
            #print(np.any(np.iscomplex(fft(trace[i,:]))),np.iscomplex(np.amax(trace[i,:])))
            #tracefft[i,:] = tracefft[i,:]*np.amax(trace[i,:])
            #tracefft[i,:] = fftshift(tracefft[i,:])
            #tracefft[i,:] = np.abs(tracefft[i,:])**2
        trace = tracefft
        
        
        if units.lower()=='mks':
            plt.xlabel(r'frequency (Hz)',fontsize=axlblfs)
            plt.ylabel(r'spectrum $B \cdot n (T)$',fontsize=axlblfs)
        elif units.lower()=='m3dc1':
            plt.xlabel(r'frequency (m3dc1 units)',fontsize=axlblfs)
            plt.ylabel(r'spectrum $B \cdot n (B0)$',fontsize=axlblfs)
    else:
        plt.xlabel(r'time',fontsize=axlblfs)
        plt.ylabel(signal+' signal',fontsize=axlblfs)
    

    for i in range(nprobes):
        plt.plot(time,trace[i,:],lw=linew,label='probe: '+str(i))
    #plt.plot(time,np.imag(tracefft))
    plt.grid(True)
    ax = plt.gca()
    ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
    plt.ticklabel_format(axis='y', style='sci',useOffset=False)
    if nprobes < 24:
        plt.legend(fontsize=legfs,ncol=nprobes%6)
    else:
        fpyl.printnote('Not showing legend because of too many probes!')
    plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    plt.show()
    
    return
