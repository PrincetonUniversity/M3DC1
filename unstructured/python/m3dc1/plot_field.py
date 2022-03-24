#!/usr/bin/env python3
#
# Coded on 11/26/2019 by:
# Andreas Kleiner:    akleiner@pppl.gov
# Ralf Mackenbach:    rmackenb@pppl.gov
# Chris Smiet    :    csmiet@pppl.gov

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
import m3dc1.fpylib as fpyl
from m3dc1.eval_field import eval_field
from m3dc1.plot_mesh import plot_mesh
#rc('text', usetex=True)




def plot_field(field, coord='scalar', row=1, sim=None, filename='C1.h5', time=None, phi=0, linear=False,
               diff=False, tor_av=1, mesh=False, bound=False, lcfs=False, units='mks',res=250, prange=None,
               cmap='viridis', cmap_midpt=None, quiet=False,
               save=False, savedir=None, pub=False, showtitle=True, shortlbl=False, ntor=None,phys=False):
    """
    Plots the field of a file. 
    
    Arguments:

    **field**
    The field that is to be plotted, i.e. 'B', 'j', etc.

    **coord**
    The chosen part of a vector field to be plotted, options are:
    'phi', 'R', 'Z', 'poloidal', 'radial', 'scalar', 'vector', 'tensor'.
    Scalar is reserved for scalar fields.
    Poloidal takes the magnetic axis at time=0 as the origin, and
    defines anti-clockwise as the positive direction.
    Radial also takes the magnetic axis as the origin.

    **row**
    For tensor fields only, row of tensor to plot. Possibilities: 1,2,3

    **sim**
    simulation sim_data object or list of sim_data objects. If none is provided, the object will be created.

    **filename**
    File name which will be read, i.e. "../../C1.h5"
    Can also be a list of two filepaths when used for diff

    **time**
    The time-slice which will be used for the field plot. If time='last', the last time slice will be used.

    **phi**
    The toroidal cross-section coordinate.

    **linear**
    Plot the linear part of the field (so the equilibrium is subtracted).
    True/False

    **diff**
    Plot the difference of two fields. 
    This could be the difference of two files (filename=['a/C1.h5','b/C1.h5']),
    or the difference between two time-slices (time=[t1,t2])
    If list for both time and filename are given file1 will be evaluated at time1,
    and file2 at time2

    **tor_av**
    Calculates the average field over tor_av number of toroidal planes

    **mesh**
    True/False
    Overlay mesh on top of plot.

    **bound**
    True/False
    Only plot boundary. Only works when mesh is true

    **lcfs**
    True/False
    Overlay last closed flux surface on top of the plot.

    **res**
    Resolution in R and Z direction.

    **cmap**
    Color map to use for the contour plot, e.g. viridis, seismic, jet, nipy_spectral

    **cmap_midpt**
    If not None, set midpoint of colormap to cmap_midpt.
    
    **save**
    True/False
    Save plot as png file.

    **savedir**
    Relative path to directory where plot is saved. If None, use current working
    directory.

    **pub**
    If True, format figure for publication (larger labels and thicker lines)

    **phys**
    Use True for plotting in physical (stellarator) geometry
    """
    
    sim, time = fpyl.setup_sims(sim,filename,time,linear,diff)
    field_idx = fpyl.get_field_idx(coord)


    # Make 3D grid based on the mesh points
    mesh_ob      = sim[0].get_mesh(quiet=quiet)
    mesh_pts     = mesh_ob.elements
    R_mesh       = mesh_pts[:,4]
    Z_mesh       = mesh_pts[:,5]
    phi0 = phi*1
    R_range      = [np.amin(R_mesh),np.amax(R_mesh)]
    Z_range      = [np.amin(Z_mesh),np.amax(Z_mesh)]
    R_linspace   = np.linspace(R_range[0], R_range[1], res, endpoint=True)
    phi_linspace = np.linspace(phi,         (360+phi), tor_av, endpoint=False)
    Z_linspace   = np.linspace(Z_range[0], Z_range[1], res, endpoint=True)
    R, phi, Z    = np.meshgrid(R_linspace, phi_linspace,Z_linspace)
    if phys:
        rst = eval_field('rst', R, phi, Z, sim=sim, file_name=file_name, time=time)
        zst = eval_field('zst', R, phi, Z, sim=sim, file_name=file_name, time=time)
        R_mesh = eval_field('rst', mesh_pts[:,4], phi0*np.ones_like(mesh_pts[:,4]), mesh_pts[:,5], sim=sim, file_name=file_name, time=time)
        Z_mesh = eval_field('zst', mesh_pts[:,4], phi0*np.ones_like(mesh_pts[:,4]), mesh_pts[:,5], sim=sim, file_name=file_name, time=time)

    

    # Get the magnetic axis at time zero, which will be used for poloidal coordinates
    if coord in ['poloidal', 'radial']:
        R_mag =  sim[0].get_time_trace('xmag').values[0]
        Z_mag =  sim[0].get_time_trace('zmag').values[0]


    # Evaluate usual vector components
    if coord not in ['poloidal', 'radial', 'vector', 'tensor']:
        # Evaluate field
        field1 = eval_field(field, R, phi, Z, coord=coord, sim=sim[0], time=time[0],quiet=quiet)
        # Evaluate second field and calculate difference between two if linear or diff is True
        if diff or linear:
            field2 = eval_field(field, R, phi, Z, coord=coord, sim=sim[1], time=time[1],quiet=quiet)
            field1 = field1 - field2
            if not quiet:
                print('[DONE]')


    # Evaluate poloidal/radial field components or all field components (coord='vector')
    if coord in ['poloidal', 'radial', 'vector']:
        field1R, field1phi, field1Z  = eval_field(field, R, phi, Z, coord='vector', sim=sim[0], time=time[0],quiet=quiet)
        if diff or linear:
            field2R, field2phi, field2Z  = eval_field(field, R, phi, Z, coord='vector', sim=sim[1], time=time[1],quiet=quiet)
            if not quiet:
                print('[DONE]')
    elif coord =='tensor':
        field1RR, field1phiR, field1ZR, field1Rphi, field1phiphi, field1Zphi, field1RZ, field1phiZ, field1ZZ  = \
            eval_field(field, R, phi, Z, coord='tensor', sim=sim[0], time=time[0],quiet=quiet)
        if row == 1:
            field1R = field1RR
            field1phi = field1phiR
            field1Z = field1ZR
        elif row == 2:
            field1R = field1Rphi
            field1phi = field1phiphi
            field1Z = field1Zphi
        elif row == 3:
            field1R = field1RZ
            field1phi = field1phiZ
            field1Z = field1ZZ
        if diff or linear:
            field2RR, field2phiR, field2ZR, field2Rphi, field2phiphi, field2Zphi, field2RZ, field2phiZ, field2ZZ  = \
                eval_field(field, R, phi, Z, coord='tensor', sim=sim[1], time=time[1],quiet=quiet)
            if row == 1:
                field2R = field2RR
                field2phi = field2phiR
                field2Z = field2ZR
            elif row == 2:
                field2R = field2Rphi
                field2phi = field2phiphi
                field2Z = field2Zphi
            elif row == 3:
                field2R = field2RZ
                field2phi = field2phiZ
                field2Z = field2ZZ

    # Evaluate poloidal component
    if coord == 'poloidal':
        theta  = np.arctan2(Z-Z_mag,R-R_mag)
        field1 = -np.sin(theta)*field1R + np.cos(theta)*field1Z
        
        # Evaluate second field and calculate difference between two if linear or diff is True
        if diff or linear:
            field2 = -np.sin(theta)*field2R + np.cos(theta)*field2Z


    # Evaluate radial component
    if coord == 'radial':
        theta  = np.arctan2(Z-Z_mag,R-R_mag)
        field1 = np.cos(theta)*field1R + np.sin(theta)*field1Z
        
        # Evaluate second field and calculate difference between two if linear or diff is True
        if diff or linear:
            field2 = np.cos(theta)*field2R + np.sin(theta)*field2Z


    # Calculate difference between two if linear or diff is True
    if diff or linear:
        if coord in ['poloidal', 'radial']:
            field1 = field1 - field2
        elif coord in ['vector', 'tensor']:
            field1R   = field1R - field2R
            field1phi = field1phi - field2phi
            field1Z   = field1Z - field2Z

    
    # Calculate average over phi of the field. For a single slice the average is equal to itself
    if coord in ['vector', 'tensor']:
        field1R_ave   = np.average(field1R, 0)
        field1phi_ave = np.average(field1phi, 0)
        field1Z_ave   = np.average(field1Z, 0)
        field1_ave    = [field1R_ave,field1phi_ave,field1Z_ave]
    else:
        field1_ave    = [np.average(field1, 0)]
    R_ave = np.average(R, 0)
    Z_ave = np.average(Z, 0)
    if phys:
        rst_ave    = np.average(rst, 0)
        zst_ave    = np.average(zst, 0)
        R_ave = np.where(np.isnan(rst_ave), R_ave, rst_ave)
        Z_ave = np.where(np.isnan(zst_ave), Z_ave, zst_ave)

    if units.lower()=='m3dc1':
        field1_ave = fpyl.get_conv_field(units,field,field1_ave,sim=sim[0])

    fieldlabel,unitlabel = fpyl.get_fieldlabel(units,field,shortlbl=shortlbl)
    cbarlbl = fieldlabel + ' (' + unitlabel + ')'


    ### Plot the field ###
    # Set font sizes and plot style parameters
    if pub:
        axlblfs = 20
        titlefs = 18
        cbarlblfs = 14
        cbarticklblfs = 14
        ticklblfs = 18
        lcfslw = 2
    else:
        axlblfs = None
        titlefs = None
        cbarlblfs = None
        cbarticklblfs = None
        ticklblfs = None
        lcfslw = 1
    
    if coord not in ['vector', 'tensor']:
        fig, axs = plt.subplots(1, 1)
        fig.set_figheight(7)
        axs = np.asarray([axs])
    else:
        fig, axs = plt.subplots(1, 3, sharey=True,figsize=(14,7))
        comp = ['R','\phi','Z']
    
    if coord in ['vector', 'tensor']:
        titlestr = field + ' at time=' + str(sim[0].timeslice)
    else:
        titlestr = field + ' at time=' + str(sim[0].timeslice)
    if linear:
        titlestr = titlestr + ' linear'
    elif diff:
        titlestr = titlestr + ' - time=' + str(sim[1].timeslice)
    try:
        species = sim[0].available_fields[field][2]
    except:
        species = None
        fpyl.printwarn('WARNING: Field not found in available_fields!')
    if species is not None:
        titlestr = titlestr+' - Species: '+str(species)
    # ToDo: Main title does not show up in vector plot
    if showtitle:
        plt.title(titlestr,fontsize=titlefs)
    
    if mesh or bound:
        meshplt = plot_mesh(mesh_ob,boundary=bound,ax=axs,meshcol='C1',pub=pub,phys=phys)
    
    for i,ax in enumerate(axs):
        if cmap_midpt is not None:
            field1_ave_clean = field1_ave[i][np.logical_not(np.isnan(field1_ave[i]))]
            norm = colors.DivergingNorm(vmin=np.amin(field1_ave_clean), vcenter=cmap_midpt, vmax=np.amax(field1_ave_clean))
            cont = ax.contourf(R_ave, Z_ave, field1_ave[i],100, cmap=cmap,norm=norm)
        else:
            if isinstance(prange,(tuple,list)):
                norm = colors.DivergingNorm(vmin=prange[0], vcenter=(prange[1]+prange[0])/2, vmax=prange[1])
                cont = ax.contourf(R_ave, Z_ave, field1_ave[i],100, cmap=cmap,norm=norm)
            else:
                cont = ax.contourf(R_ave, Z_ave, field1_ave[i],100, cmap=cmap)
        # Set and format axes limits and labels
        if phys:
            ax.set_xlim([fpyl.get_axlim(np.amin(R_mesh),'min',0.1),fpyl.get_axlim(np.amax(R_mesh),'max',0.1)])
            ax.set_ylim([fpyl.get_axlim(np.amin(Z_mesh),'min',0.1),fpyl.get_axlim(np.amax(Z_mesh),'max',0.1)])
        else:
            ax.set_xlim([fpyl.get_axlim(np.amin(R_ave),'min',0.1),fpyl.get_axlim(np.amax(R_ave),'max',0.1)])
            ax.set_ylim([fpyl.get_axlim(np.amin(Z_ave),'min',0.1),fpyl.get_axlim(np.amax(Z_ave),'max',0.1)])
        ax.set_xlabel(r'$R/m$',fontsize=axlblfs)
        ax.set_ylabel(r'$Z/m$',fontsize=axlblfs)
        ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
        if coord in ['vector', 'tensor']:
            compstr = r'$'+field+'_{'+comp[i]+'}$'
            ax.set_title(compstr,fontsize=titlefs)

        #plt.suptitle(titlestr,fontsize=titlefs)
        
        # Style and formatting of colorbar
        #cbar = fig.colorbar(cont,format=ticker.FuncFormatter(fpyl.fmt),ax=ax,cax=cax,orientation=cbarorient)
        sfmt=ticker.ScalarFormatter()
        sfmt.set_powerlimits((-3,4))
        if coord not in ['vector', 'tensor']:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.1)
            cbarorient = "vertical"
            cbarrot = 0
            cbar = fig.colorbar(cont,format=sfmt,ax=ax,cax=cax,orientation=cbarorient)
        else:
            cbarorient = "horizontal"
            cbarrot = 270
            cbar = fig.colorbar(cont,format=sfmt,ax=ax,orientation=cbarorient)
        #cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),rotation=cbarrot,fontsize=cbarticklblfs)
        cbar.ax.tick_params(labelsize=cbarticklblfs)
        cbar.ax.yaxis.offsetText.set(size=cbarticklblfs)
        cbar.set_label(cbarlbl,fontsize=cbarlblfs)
        
        # Fix for white lines in contourf plot when exported as PDF
        for c in cont.collections:
            c.set_edgecolor("face")

        ax.set_aspect('equal')
    plt.rcParams["axes.axisbelow"] = False #Allows mesh to be on top of contour plot. This option conflicts with zorder (matplotlib bug).
    
    
    if lcfs:
        if linear or diff:
            try:
                lcfs_ts_ind = np.argwhere(np.asarray([sim[0].timeslice, sim[1].timeslice]) < 1).flatten()[0]
            except:
                lcfs_ts_ind = 0
                #fpyl.printwarn('WARNING: LCFS plotted at timeslice 0')
        else:
            lcfs_ts_ind = 0
        print('LCFS time: ' + str(time[lcfs_ts_ind]))
        psi_lcfs = sim[lcfs_ts_ind].get_time_trace('psi_lcfs').values[0]
        if not quiet:
            print("Psi at LCFS: "+str(psi_lcfs))
        
        #Aphi = eval_field('A', R, phi, Z, coord='phi', sim=sim[lcfs_ts_ind], time=time[lcfs_ts_ind])
        #psifield = R_ave*Aphi
        psifield = eval_field('psi', R, phi, Z, coord='scalar', sim=sim[lcfs_ts_ind], filename=sim[lcfs_ts_ind].filename, time=time[lcfs_ts_ind])
        for i,ax in enumerate(axs):
            cont = ax.contour(R_ave, Z_ave, np.average(psifield,0),[psi_lcfs],colors='magenta',linewidths=lcfslw,zorder=10)
    
    
    plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    #if coord=='vector':
    #    fig.subplots_adjust(top=0.2)
    plt.show()
    
    if save:
        timestr = 't'
        if diff:
            timestr = timestr + str(time[0]) + '-' + str(time[1])
        else:
            timestr = timestr + str(time[0])
        if linear:
            fieldstr = 'delta_' + str(field)
        else:
            fieldstr = str(field)
        if savedir is not None:
            fieldstr = savedir + fieldstr
        if ntor is None:
            nout=sim[0].ntor
        else:
            nout=ntor
        
        plt.savefig(fieldstr + '_' + timestr + '_n'+"{:d}".format(nout)+'.png', format='png',dpi=900,bbox_inches='tight')
    
    return
