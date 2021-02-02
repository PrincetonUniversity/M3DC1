#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on February 29 2020

@author: Andreas Kleiner
"""
import math
import numpy as np
import matplotlib.pyplot as plt
import fpy
import m3dc1.fpylib as fpyl
from m3dc1.eval_field import eval_field
from m3dc1.flux_coordinates import flux_coordinates
from m3dc1.unit_conv import unit_conv




def eigenfunction(sim=None,time=1,phit=0.0,filename='C1.h5',fcoords=None,points=200,fourier=True,units='m3dc1',makeplot=True,nummodes=10,cmap='jet',pub=False):
    """
    Calculates the linear eigenfunction ~(p1-p0)

    Arguments:

    **sim**
    Simulation object at time where the eigenfunction shall be calculated.

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
    if sim != None:
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
        if time > 0 or time=='last':
            simlin = fpy.sim_data(filename,time=time)
        else:
            raise Exception('Please provide a time slice larger than 0.')
        sims = [simeq,simlin]
    
    time = sims[1].timeslice
    
    
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
    print('Evaluating first field... ', end=' ', flush=True)
    p1 = eval_field('p', fc.rpath, torphi, fc.zpath, coord='scalar', sim=sims[1], filename=filename, time=time)
    print('[DONE]')

    print('Evaluating second field... ', end=' ', flush=True)
    p0 = eval_field('p', fc.rpath, torphi, fc.zpath, coord='scalar', sim=sims[0])
    print('[DONE]')
    
    ef = p1 - p0
    if units.lower()=='m3dc1':
        ef = fpyl.get_conv_field(units,'p',ef)
    pathshape = ef.shape
    
    # Plot eigenfunction in R-Z plane
    if makeplot==True:
        # Set font sizes and plot style parameters
        if pub==False:
            axlblfs = None
            titlefs = None
            cbarlblfs = None
            cbarticklblfs = None
            ticklblfs = None
            legfs = None
            linew = 1
        elif pub==True:
            axlblfs = 20
            titlefs = 18
            cbarlblfs = 14
            cbarticklblfs = 14
            ticklblfs = 18
            legfs = 12
            linew = 2
        fig = plt.figure()
        cont = plt.contourf(fc.rpath,fc.zpath,ef,100,cmap=cmap)
        ax = plt.gca()
        ax.grid(True,zorder=10,alpha=0.5) #There seems to be a bug in matplotlib that ignores the zorder of the grid #Uncomment for CLT paper
        ax.set_xlim([fpyl.get_axlim(np.amin(fc.rpath),'min',0.1),fpyl.get_axlim(np.amax(fc.rpath),'max',0.1)])
        ax.set_ylim([fpyl.get_axlim(np.amin(fc.zpath),'min',0.1),fpyl.get_axlim(np.amax(fc.zpath),'max',0.1)])
        plt.xlabel(r'R',fontsize=axlblfs)
        plt.ylabel(r'Z',fontsize=axlblfs)
        ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
        cbar = fig.colorbar(cont,ax=ax)
        cbar.ax.tick_params(labelsize=cbarticklblfs)
        cbar.ax.yaxis.offsetText.set(size=cbarticklblfs)
        fieldlabel,unitlabel = fpyl.get_fieldlabel(units,'p')
        unitlabel = fieldlabel + ' (' + unitlabel + ')'
        cbar.set_label(unitlabel,fontsize=cbarlblfs)
        ax.set_aspect('equal')
        plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    
    # Calculate and plot poloidal Fourier spectrum
    if fourier==True:
        if points % 2 == 0:
            nspec = int(points/2)+1
        else:
            nspec = int((points+1)/2)
        spec = np.zeros((nspec,pathshape[1]))
        
        fs = list(range(0,pathshape[1]))
        for i in fs:
            spec[:,i] = np.abs(np.fft.rfft(ef[:,i]))/pathshape[1] #Using normalization by number of elements on forward transformation
        
        #Identify largest modes
        mmax = np.asarray([np.amax(spec[j,:]) for j in range(nspec)])
        mmax_ind = mmax.argsort()
        
        # Plot Fourier spectrum
        if makeplot==True:
            plt.figure()
            print(mmax_ind[-nummodes:])
            for j in range(nspec):
                if j in mmax_ind[-nummodes:]:
                    plt.plot(fc.psi_norm,spec[j,:],lw=linew,label='m='+str(j))
                else:
                    plt.plot(fc.psi_norm,spec[j,:],lw=linew)
            plt.xlabel(r'$\psi_N$',fontsize=axlblfs)
            plt.ylabel('eigenfunction',fontsize=axlblfs)
            ax = plt.gca()
            ax.grid(True,zorder=10,alpha=0.5) #There seems to be a bug in matplotlib that ignores the zorder of the grid #Uncomment for CLT paper
            ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
            plt.title('Eigenfunction - Poloidal spectrum',fontsize=titlefs)
            plt.legend(ncol=2,fontsize=legfs)
            #plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    
    # For test purposes, plot eigenfunction over theta for each radial point 
    #plt.figure()
    #for i in fs:
    #plt.plot(fc.theta,ef[:,-1])
    
    # For test purposes plot flux surfaces
    #plt.figure()
    #plt.plot(fc.rpath,fc.zpath,marker='.')
    #plt.axis('equal')
    
    if fourier==True:
        return spec
    else:
        return ef
