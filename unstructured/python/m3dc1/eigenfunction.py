#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on February 29 2020

@author: Andreas Kleiner
"""
import math
import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.fft import rfft
from scipy.fft import irfft
#from scipy.signal import find_peaks
#from scipy.signal import find_peaks_cwt
from scipy.interpolate import interp1d
from scipy.integrate import trapezoid
import fpy
import m3dc1.fpylib as fpyl
from m3dc1.eval_field import eval_field
from m3dc1.flux_coordinates import flux_coordinates
from m3dc1.flux_average import flux_average
from m3dc1.plot_mesh import plot_mesh
from m3dc1.plot_coils import plot_coils

def eigenfunction(field='p',coord='scalar',sim=None,time=1,phit=0.0,filename='C1.h5',fcoords=None,points=200,fourier=True,units='m3dc1',makeplot=True,show_res=False,
                  device='nstx',norm_to_unity=False,drop_low_m=-1,nummodes=10,cmap='jet',coils=False,mesh=False,bound=False, phys=False,pub=False,n=None,titlestr=None,save=False,savedir=None,xlimits=[None,None],colorbounds=None,extend_cbar='neither',
                  in_plot_txt=None,export=False,figsize=None,quiet=False):
    """
    Calculates the linear eigenfunction ~(p1-p0)

    Arguments:

    **field**
    The field that is to be plotted, i.e. 'B', 'j', etc.

    **coord**
    The vector component of the field to be plotted, options are:
    'phi', 'R', 'Z', 'poloidal', 'radial', 'scalar', 'vector', 'tensor'.
    Scalar is reserved for scalar fields.
    'Poloidal' takes the magnetic axis at time=0 as the origin, and
    defines anti-clockwise as the positive direction.
    Radial also takes the magnetic axis as the origin.

    **sim**
    Simulation object(s) at equilibrium and/or time where the eigenfunction shall be calculated.
    Can bei either the object itself, or a list of two objects. Order does not matter.

    **time**
    If sim=None, time slice to read.

    **phit**
    Toroidal angle where eigenfunction will be calculated

    **filename**
    If sim=None, name of file to read.

    **fcoords**
    Name of flux coordinate system: 'pest', 'boozer', 'hamada', or '' (geometric angle)

    **points**
    Number of points in theta and psi_n considered for flux coordinate calculation

    **fourier**
    If True, calculate poloidal Fourier spectrum of the eigenfunction and return spectrum.
    If False, the eigenfunction on the flux-aligned R-Z grid is returned

    **units**
    units in which the result will be calculated

    **phys**
    Use True for plotting in physical (stellarator) geometry.

    **device**
    Device (e.g. 'nstx', 'diiid') for which flux coordinates are being calculated.
    This determines the major radius, which is taken as machine specific quantity.
    Default is 'nstx'.

    **norm_to_unity**
    Normalize eigenfunction such that the maximum is in the order of unity.

    **drop_low_m**
    Removes poloidal modes from m = 0 to value of drop_low_m from poloidal spectrum and plots
    the eigenfunction with these modes subtracted (by dropping from the inverse Fourier transform).
    Can be useful, when n=0 modes are so strong that higher n/m perturbations are not visible in plot.

    **makeplot**
    If True, show plot of eigenfunction in R-Z plane and if specified also the Fourier spectrum

    **show_res**
    If True, show radial locations of resonant surfaces.

    **mesh**
    True/False
    Overlay mesh on top of plot.

    **bound**
    True/False
    Plot boundary (wall and domain boundary) without the mesh itself.

    **coils**
    True/False. Show contour of magnetic field coils in plot (imported from coil.dat).

    **nummodes**
    Number of Fourier mode to show in the plot legend

    **cmap**
    Name of colormap for countour plot of eigenfunction

    **pub**
    If True, format figure for publication (larger labels and thicker lines)

    **titlestr**
    Plot title. If None, a default title will be generated.

    **xlimits**
    x-axis limits, array of length 2, e.g. [0,1].

    **colorbounds**
    List of two values that define the lower and upper values of the colormap.

    **extend_cbar**
    Extend colorbar, either: 'neither', 'both', 'min', 'max'
    Make colorbar ends for out of range values.

    **save**
    True/False. Save plot as png file.

    **savedir**
    Relative path to directory where plot is saved. If None, use current working
    directory.

    **export**
    Export plot data to text file.

    **quiet**
    If True, suppress most output to terminal.
    """
    
    # make simulation object iterable if it is a single object and not if it is list of objects
    sims = check_sim_object(sim,time,filename)
    
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
    
    
    if show_res:
        psin,q = flux_average('q',sim=sims[0], fcoords=fcoords, linear=False, deriv=0, points=points, phit=0.0, time=sim[1].timeslice, psin_range=None, device=device, units=units)
        q_res_min = np.ceil(np.amin(q))
        q_res_max = np.floor(np.amax(q))
        q_res = np.arange(q_res_min,q_res_max+1)
        q_interp = interp1d(q,psin,kind='cubic',fill_value="extrapolate")
        psin_res = q_interp(q_res)
        print(q_res)
        print(psin_res)
        # Color cycler:
        def col_cycler(cols):
            count = 0
            while True:
                yield cols[count]
                count = (count + 1)%len(cols)
        cols = col_cycler(['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'])
    
    # Evaluate fields
    p1 = eval_field(field, fc.rpath, torphi, fc.zpath, coord=coord, sim=sims[1], filename=filename, time=sims[1].timeslice)

    p0 = eval_field(field, fc.rpath, torphi, fc.zpath, coord=coord, sim=sims[0],time=sims[0].timeslice)
    
    ef = p1 - p0
    if units.lower()=='m3dc1':
        ef = fpyl.get_conv_field(units,'p',ef,sim=sims[0])
    pathshape = ef.shape
    
    
    # Calculate and plot poloidal Fourier spectrum
    if fourier or drop_low_m>-1:
        if points % 2 == 0:
            nspec = int(points/2)+1
        else:
            nspec = int((points+1)/2)
        spec = np.zeros((nspec,pathshape[1]))
        #spec_abs = np.zeros((nspec,pathshape[1]))
        
        fs = list(range(0,pathshape[1]))
        for i in fs:
            spec[:,i] = np.abs(rfft(ef[:,i],norm='forward'))#/pathshape[1] #Applying normalization by number of elements on forward transformation
        
    
    if drop_low_m>-1:
        for i in fs:
            print(spec[:drop_low_m+1,i])
            spec[:drop_low_m+1,i] = 0.0
            ef[:,i] = irfft(spec[:,i],norm='forward')
    
    if fourier:
        for i in fs:
            spec[:,i] = np.abs(spec[:,i])
        spec_abs = np.abs(spec)
    
    
    #Determine toroidal mode number n if not specified via argument
    if n is None:
        path = sims[0].filename[:-6]
        paramfile = fpyl.get_input_parameter_file(directory=path)
        #print(slurmfile)
        n = fpyl.get_parameter_from_ascii('ntor',paramfile,quiet=True)
    
    
    # Plot eigenfunction in R-Z plane
    if makeplot:
        # Set font sizes and plot style parameters
        if pub:
            axlblfs = 20
            titlefs = 20
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
        fig = plt.figure(figsize=figsize)
        
        ef_field = np.concatenate((ef,np.reshape(ef[0,:],(1,len(ef[0,:])))))
        
        if norm_to_unity:
            #fac = 1.0*10**(-int(math.log10(np.amax(ef_field)))) # This normalization was used in early presentations, and is such that the maximum of the perturbation is of order 1e-1
            efmax = np.amax(ef_field)
            efmin = np.amin(ef_field)
            
            fac = 1.0/np.amax([efmax,np.abs(efmin)])
            plot_field = fac*ef_field
        else:
            plot_field = ef_field
        
        R_data = np.concatenate((fc.rpath,np.reshape(fc.rpath[0,:],(1,len(fc.rpath[0,:])))))
        Z_data = np.concatenate((fc.zpath,np.reshape(fc.zpath[0,:],(1,len(fc.zpath[0,:])))))
        
        if export:
            np.savetxt('ef_contours.txt', np.transpose([np.ravel(R_data),np.ravel(Z_data),np.ravel(plot_field)]))
        
        if colorbounds is None:
            cont = plt.contourf(R_data,Z_data,plot_field,100,cmap=cmap)
        else:
            cont = plt.contourf(R_data,Z_data,plot_field,100,cmap=cmap,vmin=colorbounds[0],vmax=colorbounds[1], extend=extend_cbar)
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
        if norm_to_unity:
            ticks = np.linspace(-1, 1, 11, endpoint=True)
        else:
            ticks = None
        cbar = fig.colorbar(cont,ax=ax,format=sfmt,ticks=ticks)
        cbar.ax.tick_params(labelsize=cbarticklblfs)
        cbar.ax.yaxis.offsetText.set(size=cbarticklblfs)
        fieldlabel,unitlabel = fpyl.get_fieldlabel(units,'eigenfunction')
        if norm_to_unity:
            unitlabel = 'a.u.'
        unitlabel = fieldlabel + ' (' + unitlabel + ')'
        cbar.set_label(unitlabel,fontsize=cbarlblfs)
        ax.set_aspect('equal')
        print(cbar.mappable.get_clim())
        
        
        if coils:
            plot_coils(ax=ax)
        if mesh or bound:
            if mesh and bound:#Make sure that whole mesh is plotted. If both are true, plot_mesh only plots boundary.
                bound = False
            mesh_ob      = sims[0].get_mesh(quiet=quiet)
            #mesh_pts     = mesh_ob.elements
            meshplt = plot_mesh(mesh_ob,boundary=bound,ax=ax,meshcol='C1',pub=pub,phys=phys)
        
        plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    
    
    if fourier:
        #Identify largest modes
        mmax = np.asarray([np.amax(spec_abs[j,:]) for j in range(nspec)])
        mmax_ind = mmax.argsort()
        
        # Plot Fourier spectrum
        if makeplot:
            #Plots for early presentations do not use this normalization for the Fourier harmonics
            if norm_to_unity:
                specmax = np.amax(spec_abs)
                fac = 1.0/specmax
            else:
                fac = 1.0
            plt.figure(figsize=figsize)
            #print(mmax_ind[-nummodes:])
            if export:
                np.savetxt('ef_spectrum.txt', np.insert(fac*spec_abs,0,fc.psi_norm,0))
            
            for j in range(nspec):
                if nummodes>0:
                    if j in mmax_ind[-nummodes:]:
                        plt.plot(fc.psi_norm,fac*spec_abs[j,:],lw=linew,label='m='+str(j))
                    else:
                        plt.plot(fc.psi_norm,fac*spec_abs[j,:],lw=linew)
                else:
                    plt.plot(fc.psi_norm,fac*spec_abs[j,:],lw=linew)
            
            ax = plt.gca()
            if show_res:
                for i,pr in enumerate(psin_res):
                    pltcol = next(cols)
                    plt.axvline(x=pr,linestyle='--',color=pltcol)
                    ypos = 0.98*ax.get_ylim()[-1]
                    ax.annotate('q='+str(q_res[i])+', m='+str(q_res[i]*n), xy=(pr, ypos), rotation=90, va='top', ha='left')
            
            plt.xlabel(r'$\psi_N$',fontsize=axlblfs)
            plt.ylabel('eigenfunction (a. u.)',fontsize=axlblfs)
            ax.set_xlim(left=xlimits[0],right=xlimits[1])
            ax.grid(True,zorder=10,alpha=0.5) #There seems to be a bug in matplotlib that ignores the zorder of the grid #Uncomment for CLT paper
            ax.yaxis.set_major_formatter(sfmt)
            ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
            
            #if pub!=True:
            if n is not None and titlestr is None:
                titlestr = 'n='+str(int(n))
            plt.title(titlestr,fontsize=titlefs)
            
            if in_plot_txt is not None:
                plt.text(0.06, 0.9,in_plot_txt, ha='center', va='center', transform=ax.transAxes,fontsize=titlefs+2)
            
            if nummodes>0:
                plt.legend(ncol=2,fontsize=legfs)
            plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
            
            if save:
                timestr = 't'
                timestr = timestr + str(sim[1].timeslice)
                imgname = 'spectrum'
                if savedir is not None:
                    imgname = savedir + imgname
                plt.savefig(imgname + '_' + timestr + '_n'+"{:d}".format(n)+'.png', format='png',dpi=300,bbox_inches='tight')
    
    
    
    
    # For test purposes, plot eigenfunction over theta for each radial point 
    #plt.figure()
    #for i in fs:
    #plt.plot(fc.theta,ef[:,-1])
    
    # For testing purposes plot flux surfaces
    #plt.figure()
    #plt.plot(fc.rpath,fc.zpath,marker='.')
    #plt.axis('equal')
    if fourier:
        return spec_abs
    else:
        return ef



def check_sim_object(sim,time,filename):
    """
    Checks simulation object and makes simulation object iterable
    if it is a single object but not if it is list of objects.

    Arguments:

    **sim**
    Simulation object(s) at equilibrium and/or time where the eigenfunction shall be calculated.
    Can bei either the object itself, or a list of two objects. Order does not matter.

    **time**
    If sim=None, time slice to read.

    **filename**
    If sim=None, name of file to read.
    """
    if sim is not None:
        if not isinstance(sim, (tuple, list,np.ndarray)):
            
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
    return sims


def eigenfunction_vs_psin(sim=None,time=1,filename='C1.h5',fcoords=None,points=801,units='m3dc1',psin_ped_top=0.86,phys=False,
                          norm_to_unity=False,drop_low_m=-1,makeplot=True,quiet=False, n=None,titlestr=None,save=False,
                          savedir=None,in_plot_txt=None,figsize=None):
    """
    Returns and optionally plots the total pressure perturbation
    as a function of psin_n.

    Arguments:

    **sim**
    Simulation object(s) at equilibrium and/or time where the eigenfunction shall be calculated.
    Can bei either the object itself, or a list of two objects. Order does not matter.

    **time**
    If sim=None, time slice to read.

    **filename**
    If sim=None, name of file to read.

    **fcoords**
    Name of flux coordinate system: 'pest', 'boozer', 'hamada', or '' (geometric angle)

    **points**
    Number of points in theta and psi_n considered for flux coordinate calculation

    **units**
    units in which the result will be calculated

    **psin_ped_top**
    Value of normalized psi where pedestal top is assumed.

    **phys**
    Use True for plotting in physical (stellarator) geometry.

    **norm_to_unity**
    Normalize eigenfunction such that the maximum is in the order of unity.

    **drop_low_m**
    Removes poloidal modes from m = 0 to value of drop_low_m from poloidal spectrum and plots
    the eigenfunction with these modes subtracted (by dropping from the inverse Fourier transform).
    Can be useful, when n=0 modes are so strong that higher n/m perturbations are not visible in plot.

    **makeplot**
    If True, show radial plot of total perturbation.

    **quiet**
    If True, suppress most output to terminal.

    **titlestr**
    Plot title. If None, a default title will be generated.
    """
    sims = check_sim_object(sim,time,filename)
    
    spec = eigenfunction(sim=sims,time=time,filename=filename,fcoords=fcoords,points=points,units=units,fourier=True,norm_to_unity=norm_to_unity,drop_low_m=drop_low_m,quiet=quiet,phys=phys,makeplot=False)
    
    _,_,efsum = mode_type(spec,sims[0],psin_ped_top=psin_ped_top,makeplot=makeplot)
    
    fc = sims[0].fc
    
    return fc.psi_norm, efsum


def mode_type(spec,sim,psin_ped_top=0.86,makeplot=False):
    """
    Determines whether a mode is an edge or code mode.
    Classification is done based on radial peak of the
    perturbation and comparing how much of the integral of
    the perturbation lies in the pedestal and the core.

    Arguments:

    **spec**
    Fourier spectrum of the mode calculated by eigenfunction().

    **sim**
    Simulation object containing the flux coordinate data.

    **psin_ped_top**
    Value of normalized psi where pedestal top is assumed.

    **makeplot**
    If True, plot total perturbation.
    """
    
    spec_abs = np.abs(spec)
    
    pathshape = spec_abs.shape
    fc = sim.fc
    nspec = pathshape[0]
    
    # Sum up all Fourier components without considering the phase in order to locate the perturbation.
    efsum = np.zeros(pathshape[1])
    for i in range(pathshape[1]):
        efsum[i] = np.sum(spec_abs[:,i])
    # Plot sum of Fourier components:
    if makeplot:
        plt.figure()
        plt.plot(fc.psi_norm,efsum,lw=2,c='C2')
        plt.xlabel(r'$\psi_N$',fontsize=20)
        plt.ylabel(r'$\sum_{m} |\xi_{m}|(r)$ (a.u.)',fontsize=20)
        ax = plt.gca()
        ax.grid(True,zorder=10,alpha=0.5) #There seems to be a bug in matplotlib that ignores the zorder of the grid #Uncomment for CLT paper
        ax.tick_params(axis='both', which='major', labelsize=18)
        plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    
    # Calculate derivative of sum of modes
    #efsumprime = fpyl.deriv(efsum,fc.psi_norm)
    #plt.figure()
    #plt.plot(fc.psi_norm,efsumprime)
    
    
    # Determine location of total maximum
    efmax = np.amax(efsum)
    efmax_ind = np.argmax(efsum)
    psinatmax = fc.psi_norm[efmax_ind]
    #print(efmax)
    #print(psinatmax)
    
    
    # Determine maximum value in edge region
    
    psinedge = fpyl.find_nearest(fc.psi_norm,psin_ped_top)
    psinedge_ind = fpyl.get_ind_at_val(fc.psi_norm,psinedge)
    print('--------------------------')
    print(psin_ped_top,psinedge_ind)
    print('--------------------------')
    efmax_edge = np.amax(efsum[psinedge_ind:])
    efmax_core = np.amax(efsum[:psinedge_ind])
    rel_mode_ampl = efmax_edge/efmax_core
    
    #print(efmax_edge)
    print('Relative amplitude of edge mode: '+str(rel_mode_ampl))
    
    
    
    #Determine how much of the mode lies inside and outside the pedestal
    #by integrating the perturbation from 0 to psin_ped_top and
    #from psin_ped_top to 1
    pert = np.zeros(pathshape[1])
    for i in range(pathshape[1]):
        pert[i] = np.sum(spec[:,i])
    
    ped_top_ind = fpyl.get_ind_at_val(fc.psi_norm,fpyl.find_nearest(fc.psi_norm,psin_ped_top))
    print(ped_top_ind,fc.psi_norm[ped_top_ind])
    
    temp1 = trapezoid(pert[:ped_top_ind+1],fc.psi_norm[:ped_top_ind+1])
    temp2 = trapezoid(pert[ped_top_ind+1:],fc.psi_norm[ped_top_ind+1:])
    ped_loc = temp2/(temp1+temp2)
    print(temp1,temp2)
    #fc.psi_norm,pert
    
    
    
    if rel_mode_ampl>8:
        mtype = 1
    elif rel_mode_ampl<8 and rel_mode_ampl>1:
        if ped_loc > 0.5:
            mtype = -1
        else:
            mtype = -2
        fpyl.printwarn('Dominant PB mode. Weak core mode detected.')
    elif rel_mode_ampl<=1 and rel_mode_ampl>0.1:
        mtype = -2
        fpyl.printwarn('Dominant core mode. PB mode is weak.')
    else:
        mtype = 0
        fpyl.printwarn('Mode is not a PB mode.')
    
    
    
    
    
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

    return mtype,ped_loc,efsum
