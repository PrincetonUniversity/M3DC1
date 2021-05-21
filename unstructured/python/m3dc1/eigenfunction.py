#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on February 29 2020

@author: Andreas Kleiner
"""
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
#from scipy.signal import find_peaks
#from scipy.signal import find_peaks_cwt
import fpy
import m3dc1.fpylib as fpyl
from m3dc1.eval_field import eval_field
from m3dc1.flux_coordinates import flux_coordinates



def eigenfunction(sim=None,time=1,phit=0.0,filename='C1.h5',fcoords=None,points=200,fourier=True,units='m3dc1',makeplot=True,norm_to_unity=False,nummodes=10,cmap='jet',pub=False,n=None):
    """
    Calculates the linear eigenfunction ~(p1-p0)

    Arguments:

    **sim**
    Simulation object(s) at equilibrium and/or time where the eigenfunction shall be calculated.
    Can bei either the object itself, or a list of two objects. Order does not matter.

    **time**
    If sim=None, time slice to read

    **phit**
    Toroidal angle where eigenfunction will be calculated

    **filename**
    If sim=None, name of file to read

    **fcoords**
    Name of flux coordinate system: 'pest', 'boozer', 'hamada', or '' (geometric angle)

    **points**
    Number of points in theta and psi_n considered for flux coordinate calculation

    **fourier**
    If True, calculate poloidal Fourier spectrum of the eigenfunction and return spectrum.
    If False, the eigenfunction on the flux-aligned R-Z grid is returned

    **units**
    units in which the result will be calculated

    **makeplot**
    If True, show plot of eigenfunction in R-Z plane and if specified also the Fourier spectrum

    **nummodes**
    Number of Fourier mode to show in the plot legend

    **cmap**
    Name of colormap for countour plot of eigenfunction

    **pub**
    If True, format figure for publication (larger labels and thicker lines)
    """
    
    # make simulation object iterable if it is a single object and not if it is list of objects
    if sim is not None:
        if not isinstance(sim, (tuple, list)):
            if isinstance(sim,fpy.sim_data):
                if sim.timeslice==-1 or sim.timeslice==0:
                    simlin = fpy.sim_data(filename,time=time)
                    sims = [sim,simlin]
                else:
                    simeq = fpy.sim_data(filename,time=-1)
                    sims = [simeq,sim]
            else:
                raise Exception('sim is not a fpy.sim_data object.')
        else:
            if len(sim)>2:
                raise Exception('Please provide not more than 2 simulation objects.')
            else:
                if isinstance(sim[0],fpy.sim_data) and isinstance(sim[1],fpy.sim_data):
                    if sim[0].timeslice==-1 or sim[0].timeslice==0:
                        sims=sim
                    elif sim[1].timeslice==-1 or sim[1].timeslice==0:
                        sims = []
                        sims.append(sim[1])
                        sims.append(sim[0])
                        sims = np.asarray(sims)
                    else:
                        raise Exception('Please provide 1 simulation at time=-1!')
                else:
                    raise Exception('sim is not a list of fpy.sim_data objects.')
    else:
        simeq = fpy.sim_data(filename,time=-1)
        if (isinstance(time, str) and time=='last') or (isinstance(time, int) and time > 0):
            simlin = fpy.sim_data(filename,time=time)
        else:
            raise Exception('Please provide a time slice larger than 0.')
        sims = [simeq,simlin]
    
    # Calculate flux coodinates if it was not calculated yet or a different flux coordinate system than sim.fc.fcoords is desired
    
    if isinstance(sims[0].fc,fpy.flux_coordinates)==False or (fcoords!=None and (sims[0].fc.fcoords!=fcoords)):
            sims[0] = flux_coordinates(sim=sims[0], fcoords=fcoords, phit=phit, points=points)
    else:
        if sims[0].fc.m != points:
            sims[0] = flux_coordinates(sim=sims[0], fcoords=fcoords, phit=phit, points=points)
    fc = sims[0].fc
    
    torphi = np.zeros_like(fc.rpath)
    if phit != 0.0:
        torphi.fill(phit)
    
    # Evaluate fields
    p1 = eval_field('p', fc.rpath, torphi, fc.zpath, coord='scalar', sim=sims[1], filename=filename, time=sims[1].timeslice)

    p0 = eval_field('p', fc.rpath, torphi, fc.zpath, coord='scalar', sim=sims[0],time=sims[0].timeslice)
    
    ef = p1 - p0
    if units.lower()=='m3dc1':
        ef = fpyl.get_conv_field(units,'p',ef,sim=sims[0])
    pathshape = ef.shape
    
    # Plot eigenfunction in R-Z plane
    if makeplot:
        # Set font sizes and plot style parameters
        if pub:
            axlblfs = 20
            titlefs = 16
            cbarlblfs = 14
            cbarticklblfs = 14
            ticklblfs = 18
            legfs = 12
            linew = 2
        else:
            axlblfs = None
            titlefs = None
            cbarlblfs = None
            cbarticklblfs = None
            ticklblfs = None
            legfs = None
            linew = 1
        fig = plt.figure()
        
        ef_field = np.concatenate((ef,np.reshape(ef[0,:],(1,len(ef[0,:])))))
        if norm_to_unity:
            fac = 1.0*10**(-int(math.log10(np.amax(ef_field))))
        else:
            fac = 1.0
        cont = plt.contourf(np.concatenate((fc.rpath,np.reshape(fc.rpath[0,:],(1,len(fc.rpath[0,:]))))),np.concatenate((fc.zpath,np.reshape(fc.zpath[0,:],(1,len(fc.zpath[0,:]))))),fac*ef_field,100,cmap=cmap)
        ax = plt.gca()
        ax.grid(True,zorder=10,alpha=0.5) #There seems to be a bug in matplotlib that ignores the zorder of the grid #Uncomment for CLT paper
        ax.set_xlim([fpyl.get_axlim(np.amin(fc.rpath),'min',0.1),fpyl.get_axlim(np.amax(fc.rpath),'max',0.1)])
        ax.set_ylim([fpyl.get_axlim(np.amin(fc.zpath),'min',0.1),fpyl.get_axlim(np.amax(fc.zpath),'max',0.1)])
        plt.xlabel(r'R (m)',fontsize=axlblfs)
        plt.ylabel(r'Z (m)',fontsize=axlblfs)
        ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
        if n is not None:
            plt.title('n='+str(int(n)),fontsize=titlefs)
        
        sfmt=ticker.ScalarFormatter()
        sfmt.set_powerlimits((-3,4))
        cbar = fig.colorbar(cont,ax=ax,format=sfmt)
        cbar.ax.tick_params(labelsize=cbarticklblfs)
        cbar.ax.yaxis.offsetText.set(size=cbarticklblfs)
        fieldlabel,unitlabel = fpyl.get_fieldlabel(units,'eigenfunction')
        unitlabel = fieldlabel + ' (' + unitlabel + ')'
        cbar.set_label(unitlabel,fontsize=cbarlblfs)
        ax.set_aspect('equal')
        plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    
    # Calculate and plot poloidal Fourier spectrum
    if fourier:
        if points % 2 == 0:
            nspec = int(points/2)+1
        else:
            nspec = int((points+1)/2)
        spec = np.zeros((nspec,pathshape[1]))
        
        fs = list(range(0,pathshape[1]))
        for i in fs:
            spec[:,i] = np.abs(np.fft.rfft(ef[:,i]))/pathshape[1] #Applying normalization by number of elements on forward transformation
        
        #Identify largest modes
        mmax = np.asarray([np.amax(spec[j,:]) for j in range(nspec)])
        mmax_ind = mmax.argsort()
        
        # Plot Fourier spectrum
        if makeplot:
            plt.figure()
            #print(mmax_ind[-nummodes:])
            for j in range(nspec):
                if nummodes>0:
                    if j in mmax_ind[-nummodes:]:
                        plt.plot(fc.psi_norm,fac*spec[j,:],lw=linew,label='m='+str(j))
                    else:
                        plt.plot(fc.psi_norm,fac*spec[j,:],lw=linew)
                else:
                    plt.plot(fc.psi_norm,fac*spec[j,:],lw=linew)
            plt.xlabel(r'$\psi_N$',fontsize=axlblfs)
            plt.ylabel('eigenfunction (a. u.)',fontsize=axlblfs)
            ax = plt.gca()
            ax.grid(True,zorder=10,alpha=0.5) #There seems to be a bug in matplotlib that ignores the zorder of the grid #Uncomment for CLT paper
            ax.yaxis.set_major_formatter(sfmt)
            ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
            
            #if pub!=True:
            if n is not None:
                plt.title('n='+str(int(n)),fontsize=titlefs)
            if nummodes>0:
                plt.legend(ncol=2,fontsize=legfs)
            plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    
    
    
    
    # For test purposes, plot eigenfunction over theta for each radial point 
    #plt.figure()
    #for i in fs:
    #plt.plot(fc.theta,ef[:,-1])
    
    # For test purposes plot flux surfaces
    #plt.figure()
    #plt.plot(fc.rpath,fc.zpath,marker='.')
    #plt.axis('equal')
    if fourier:
        return spec
    else:
        return ef



def mode_type(spec,sim,psin_ped_top=0.86):
    """
    Determines whether a mode is an edge or code mode

    Arguments:

    **spec**
    Fourier spectrum of the mode calculated by eigenfunction().

    **sim**
    Simulation object containing the flux coordinate data.

    **psin_ped_top**
    Value of normalized psi where pedestal top is assumed.
    """
    pathshape = spec.shape
    fc = sim.fc
    nspec = pathshape[0]
    
    # Sum up all Fourier components without considering the phase in order to locate the perturbation.
    efsum = np.zeros(pathshape[1])
    for i in range(pathshape[1]):
        efsum[i] = np.sum(spec[:,i])
    # Plot sum of Fourier components:
    #plt.figure()
    #plt.plot(fc.psi_norm,efsum,lw=2,c='C2')
    #plt.xlabel(r'$\psi_N$',fontsize=20)
    #plt.ylabel(r'$\sum_{m} |\xi_{m}|(r)$ (a.u.)',fontsize=20)
    #ax = plt.gca()
    #ax.grid(True,zorder=10,alpha=0.5) #There seems to be a bug in matplotlib that ignores the zorder of the grid #Uncomment for CLT paper
    #ax.tick_params(axis='both', which='major', labelsize=18)
    #plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    
    # Calculate derivative of sum of modes
    efsumprime = fpyl.deriv(efsum,fc.psi_norm)
    #plt.figure()
    #plt.plot(fc.psi_norm,efsumprime)
    
    
    # Determine location of total maximum
    efmax = np.amax(efsum)
    efmax_ind = np.argmax(efsum)
    psinatmax = fc.psi_norm[efmax_ind]
    #print(efmax)
    #print(psinatmax)
    
    
    #psin_ped_top = 0.86 #Min value of psin that is considered edge region
    
    # Determine maximum value in edge region
    psinedge = fpyl.find_nearest(fc.psi_norm,psin_ped_top)
    psinedge_ind = fpyl.get_ind_at_val(fc.psi_norm,psinedge)
    #print(psinedge_ind)
    efmax_edge = np.amax(efsum[psinedge_ind:])
    efmax_core = np.amax(efsum[:psinedge_ind])
    rel_mode_ampl = efmax_edge/efmax_core
    
    #print(efmax_edge)
    print('Relative amplitude of edge mode: '+str(rel_mode_ampl))
    
    if rel_mode_ampl>8:
        pbmode = 1
    elif rel_mode_ampl<8 and rel_mode_ampl>1:
        pbmode = -1
        fpyl.printwarn('Dominant PB mode. Weak core mode detected.')
    elif rel_mode_ampl<=1 and rel_mode_ampl>0.1:
        pbmode = -2
        fpyl.printwarn('Dominant core mode. PB mode is weak.')
    else:
        pbmode = 0
        fpyl.printwarn('Mode is not a PB mode.')
    
    
    #ToDO: Analyze individual Fourier modes and see where they peak
    
    
    
    
    
    
    #pwid = 0.05*len(efsum)
    #pprom = efmax/100
    #pheight = efmax/100
    #temp = find_peaks(efsum,height=pheight,width=pwid)
    #print(temp[0])
    
    #peaks = find_peaks_cwt(efsum,widths=[pwid])
    #print(peaks)
    
    ## Determine eigenfunction at location of peaks
    #ef_at_peaks = [efsum[p] for p in peaks]
    #psin_at_peaks = [fc.psi_norm[p] for p in peaks]
    #print(psin_at_peaks)
    #sw = 10 # Size of window around peak
    #ef_around_peaks = np.zeros(len(peaks))
    #for i,p in enumerate(peaks):
        #if p<sw:
            #ef_around_peaks[i] = np.amax(efsum[0:p+sw])
        #elif p>len(efsum)-sw+1:
            #ef_around_peaks[i] = np.amax(efsum[p-10:])
        #else:
            #ef_around_peaks[i] = np.amax(efsum[p-sw:p+sw])
    
    #print('Relative amplitude of edge mode: '+str(ef_around_peaks[-1]/ef_around_peaks[0]))
    ## Check each Fourier mode for location of its maxima
    
    #isatedge = np.zeros(nspec, dtype=bool)
    #plt.figure()
    #for j in range(nspec):
        #max_ind = np.argmax(spec[j,:])
        #psinmodeatmax = fc.psi_norm[max_ind]
        ##print(psinmodeatmax)
        #if psinmodeatmax>=psin_ped_top:
            #isatedge[j] = True
            #plt.plot(fc.psi_norm,spec[j,:])
            
    ##print(*isatedge)
    #efnoedgesum = np.zeros(pathshape[1])
    #for i in range(pathshape[1]):
        #for j in range(nspec):
            #if isatedge[j]==False:
                #efnoedgesum[i] = efnoedgesum[i]+spec[j,i]
    
    #plt.figure()
    #plt.plot(fc.psi_norm,efnoedgesum)
    
    #Subtract all modes with peak at edge from total.

    return pbmode
