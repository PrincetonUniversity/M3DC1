#!/usr/bin/env python3
#
# Coded on 11/26/2019 by:
# Andreas Kleiner:    akleiner@pppl.gov
# Ralf Mackenbach:    rmackenb@pppl.gov
# Chris Smiet    :    csmiet@pppl.gov

import fpy
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




def plot_field(field, coord='scalar', row=1, sim=None, file_name='C1.h5', time=0, phi=0, linear=False, diff=False, tor_av=1, mesh=False, bound=False, lcfs=False, units='mks',res=250, prange=None, cmap='viridis', cmap_midpt=None, save=False, savedir=None,pub=False,showtitle=True,n=None,phys=False):
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

    **file_name**
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
    This could be the difference of two files (file_name=['a/C1.h5','b/C1.h5']),
    or the difference between two time-slices (time=[t1,t2])
    If list for both time and file_name are given file1 will be evaluated at time1,
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
    # make file name iterable if it is a string and not a list of strings
    file_name = (file_name,) if not isinstance(file_name, (tuple, list)) else file_name
    
    # make time iterable if it is a single int and not if it is list of ints
    time = [time,] if not isinstance(time, (tuple, list)) else time
    
    
    ### Input error handling ###
    #if linear == True and (time[0] == 0 and not isinstance(sim,fpy.sim_data)):
    #    raise Exception('linear=True. Please enter a time slice > 0.')
    
    # make simulation object iterable if it is a single object and not if it is list of objects
    if sim != None:
        if not isinstance(sim, (tuple, list)):
            if isinstance(sim,fpy.sim_data):
                sims = [sim,None]
                time[0] = sims[0].timeslice
            else:
                raise Exception('sim is not a fpy.sim_data object.')
        else:
            if isinstance(sim[0],fpy.sim_data) and isinstance(sim[1],fpy.sim_data):
                time = np.zeros(2)
                # If linear=True, make sure that the equilibrium time slice is the second element in sims
                if linear == True:
                    if sim[0].timeslice==-1 or sim[0].timeslice==0:
                        sims = []
                        sims.append(sim[1])
                        sims.append(sim[0])
                        sims = np.asarray(sims)
                    elif sim[1].timeslice==-1 or sim[1].timeslice==0:
                        sims=sim
                    else:
                        raise Exception('Please provide 1 simulation at time=-1!')
                else:
                    sims=sim
                time[0] = int(sims[0].timeslice)
                time[1] = int(sims[1].timeslice)
            else:
                raise Exception('sim is not a list of fpy.sim_data objects.')
    else:
        sim = fpy.sim_data(file_name[0],time=time[0])
        sims = [sim,None]
        time[0] = int(sims[0].timeslice)
    
    ### More input error handling ###
    if linear == True and time[0] == -1:
        raise Exception('linear=True. Please enter a time slice > 0.')
    if diff==True:
        if len(time)!=2 and (isinstance(sims[0],fpy.sim_data)==False and isinstance(sims[1],fpy.sim_data)==False) :
            raise Exception('Please input two times for differences or specify two sim_data objects.')
    
    if diff==True and linear==True:
        raise Exception('Please choose diff or linear, not both.')

    if diff==False:
        if (len(file_name)>1 or len(time)>1 and (isinstance(sims[0],fpy.sim_data)==False and isinstance(sims[1],fpy.sim_data)==False)):
            raise Exception('Multiple file/time slices detected. Please set diff=True or input single slices')
    
    

    if (coord == 'R' or coord == 'scalar'):
        field_idx = 0
    elif coord == 'phi':
        field_idx = 1
    elif coord == 'Z':
        field_idx = 2
    elif coord=='poloidal':
        field_idx = None
    elif coord=='radial':
        field_idx = None
    elif coord == 'vector':
        field_idx = None
    elif coord == 'tensor':
        field_idx = None
    else:
        raise Exception('Please enter valid coordinate. Accepted: \'R\', \'phi\', \'Z\', \'poloidal\', \'radial\', \'scalar\', \'vector\'')
    ### End input error handling ###
    
    
    # If linear, the difference needs to be taken with time -1
    if linear==True:
        time = [int(time[0]),int(-1)]
    
    
    # If either file_name or time is a list, we will convert both of them to lists of length two.
    if (len(file_name)==2 and len(time)==1):
        time = [time, time]
    if (len(file_name)==1 and len(time)==2):
        file_name = [file_name[0],file_name[0]]
    
    # Make 3D grid based on the mesh points
    mesh_ob      = sims[0].get_mesh(int(time[0]))
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
    if phys==True:
        rst = eval_field('rst', R, phi, Z, sim=sim, file_name=file_name, time=time)
        zst = eval_field('zst', R, phi, Z, sim=sim, file_name=file_name, time=time)
        R_mesh = eval_field('rst', mesh_pts[:,4], phi0*np.ones_like(mesh_pts[:,4]), mesh_pts[:,5], sim=sim, file_name=file_name, time=time)
        Z_mesh = eval_field('zst', mesh_pts[:,4], phi0*np.ones_like(mesh_pts[:,4]), mesh_pts[:,5], sim=sim, file_name=file_name, time=time)

    

    # Get the magnetic axis at time zero, which will be used for poloidal coordinates
    if coord == 'poloidal' or coord == 'radial':
        R_mag =  sims[0].get_time_traces('xmag').values[0]
        Z_mag =  sims[0].get_time_traces('zmag').values[0]


    # Evaluate usual vector components
    if coord != 'poloidal' and coord !='radial' and coord !='vector' and coord !='tensor':
        # Evaluate field
        print('Evaluating field... ', end=' ', flush=True)
        field1 = eval_field(field, R, phi, Z, coord=coord, sim=sims[0], file_name=file_name[0], time=time[0])
        print('[DONE]')
        # Evaluate second field and calculate difference between two if linear or diff is True
        if (diff == True or linear==True):
            print('Evaluating second field... ', end=' ', flush=True)
            field2 = eval_field(field, R, phi, Z, coord=coord, sim=sims[1], file_name=file_name[1], time=time[1])
            field1 = field1 - field2
            print('[DONE]')


    # Evaluate poloidal/radial field components or all field components (coord='vector')
    if coord == 'poloidal' or coord == 'radial' or coord == 'vector':
        print('Evaluating ' + str(coord) + ' field... ')
        field1R, field1phi, field1Z  = eval_field(field, R, phi, Z, coord='vector', sim=sims[0], file_name=file_name[0], time=time[0])
        if (diff == True or linear==True):
            print('Evaluating second field... ', end=' ', flush=True)
            field2R, field2phi, field2Z  = eval_field(field, R, phi, Z, coord='vector', sim=sims[1], file_name=file_name[1], time=time[1])
            print('[DONE]')
    elif coord =='tensor':
        print('Evaluating ' + str(coord) + ' field... ')
        field1RR, field1phiR, field1ZR, field1Rphi, field1phiphi, field1Zphi, field1RZ, field1phiZ, field1ZZ  = eval_field(field, R, phi, Z, coord='tensor', sim=sims[0], file_name=file_name[0], time=time[0])
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
        if (diff == True or linear==True):
            print('Evaluating second field... ', end=' ', flush=True)
            field2RR, field2phiR, field2ZR, field2Rphi, field2phiphi, field2Zphi, field2RZ, field2phiZ, field2ZZ  = eval_field(field, R, phi, Z, coord='tensor', sim=sims[1], file_name=file_name[1], time=time[1])
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
            print('[DONE]')

    # Evaluate poloidal component
    if coord == 'poloidal':
        theta                        = np.arctan2(Z-Z_mag,R-R_mag)
        field1                       = -np.sin(theta)*field1R + np.cos(theta)*field1Z
        
        # Evaluate second field and calculate difference between two if linear or diff is True
        if (diff == True or linear==True):
            field2                       = -np.sin(theta)*field2R + np.cos(theta)*field2Z


    # Evaluate radial component
    if coord == 'radial':
        theta                        = np.arctan2(Z-Z_mag,R-R_mag)
        field1                       = np.cos(theta)*field1R + np.sin(theta)*field1Z
        
        # Evaluate second field and calculate difference between two if linear or diff is True
        if (diff == True or linear==True):
            field2                       = np.cos(theta)*field2R + np.sin(theta)*field2Z


    # Calculate difference between two if linear or diff is True
    if (diff == True or linear==True):
        if coord == 'poloidal' or coord == 'radial':
            field1                       = field1 - field2
        elif coord == 'vector' or coord =='tensor':
            field1R                       = field1R - field2R
            field1phi                     = field1phi - field2phi
            field1Z                       = field1Z - field2Z

    
    # Calculate average over phi of the field. For a single slice the average is equal to itself
    if coord == 'vector' or coord =='tensor':
        field1R_ave = np.average(field1R, 0)
        field1phi_ave = np.average(field1phi, 0)
        field1Z_ave = np.average(field1Z, 0)
        field1_ave = [field1R_ave,field1phi_ave,field1Z_ave]
    else:
        field1_ave     = [np.average(field1, 0)]

    R_ave          = np.average(R, 0)
    Z_ave          = np.average(Z, 0)
    if phys==True:
        rst_ave    = np.average(rst, 0)
        zst_ave    = np.average(zst, 0)
        R_ave = np.where(np.isnan(rst_ave), R_ave, rst_ave)
        Z_ave = np.where(np.isnan(zst_ave), Z_ave, zst_ave)

    if units.lower()=='m3dc1':
        field1_ave = fpyl.get_conv_field(units,field,field1_ave)

    fieldlabel,unitlabel = fpyl.get_fieldlabel(units,field)
    if units.lower()=='m3dc1':
        unitlabel = fieldlabel + ' (' + unitlabel + ')'


    ### Plot the field ###
    # Set font sizes and plot style parameters
    if pub==False:
        axlblfs = None
        titlefs = None
        cbarlblfs = None
        cbarticklblfs = None
        ticklblfs = None
        lcfslw = 1
    elif pub==True:
        axlblfs = 20
        titlefs = 18
        cbarlblfs = 14
        cbarticklblfs = 14
        ticklblfs = 18
        lcfslw = 2
    
    if coord != 'vector' and coord != 'tensor':
        fig, axs = plt.subplots(1, 1)
        fig.set_figheight(7)
        axs = np.asarray([axs])
    else:
        fig, axs = plt.subplots(1, 3, sharey=True,figsize=(14,7))
        comp = ['R','\phi','Z']
    
    if coord == 'vector' or coord == 'tensor':
        titlestr = field + ' at time=' + str(time[0])
    else:
        titlestr = field + ' at time=' + str(time[0])
    if linear==True:
        titlestr = titlestr + ' linear'
    elif diff==True:
        titlestr = titlestr + ' - time=' + str(time[1])
    species = sims[0].available_fields[field][2]
    if species != None:
        titlestr = titlestr+' - Species: '+str(species)
    # ToDo: Main title does not show up in vector plot
    if showtitle==True:
        plt.title(titlestr,fontsize=titlefs)
    
    if mesh == True or bound==True:
        meshplt = plot_mesh(mesh_ob,boundary=bound,ax=axs,pub=pub,phys=phys)
    
    for i,ax in enumerate(axs):
        if cmap_midpt!=None:
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
        ax.set_xlim([fpyl.get_axlim(np.amin(R_mesh),'min',0.1),fpyl.get_axlim(np.amax(R_mesh),'max',0.1)])
        ax.set_ylim([fpyl.get_axlim(np.amin(Z_mesh),'min',0.1),fpyl.get_axlim(np.amax(Z_mesh),'max',0.1)])
        ax.set_xlabel(r'$R/m$',fontsize=axlblfs)
        ax.set_ylabel(r'$Z/m$',fontsize=axlblfs)
        ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
        if coord == 'vector' or coord == 'tensor':
            compstr = r'$'+field+'_{'+comp[i]+'}$'
            ax.set_title(compstr,fontsize=titlefs)

        #plt.suptitle(titlestr,fontsize=titlefs)
        #plt.title('(b) M3D-C1',fontsize=titlefs) #For CLT paper
        
        # Style and formatting of colorbar
        #cbar = fig.colorbar(cont,format=ticker.FuncFormatter(fpyl.fmt),ax=ax,cax=cax,orientation=cbarorient)
        sfmt=ticker.ScalarFormatter()
        #sfmt.set_powerlimits((-6, -6)) # For CLT paper
        sfmt.set_powerlimits((-3,4))
        if coord != 'vector' and coord != 'tensor':
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
        cbar.set_label(unitlabel,fontsize=cbarlblfs) #ToDo: Comment this line for CLT paper
        
        # Fix for white lines in contourf plot when exported as PDF
        for c in cont.collections:
            c.set_edgecolor("face")

        ax.set_aspect('equal')
        ax.grid(True,zorder=10,alpha=0.5) #There seems to be a bug in matplotlib that ignores the zorder of the grid #Uncomment for CLT paper
    plt.rcParams["axes.axisbelow"] = False #Allows mesh to be on top of contour plot. This option conflicts with zorder (matplotlib bug).
    
    
    if lcfs == True:
        psi_lcfs = sims[0].get_time_traces('psi_lcfs').values[0]
        print("Psi at LCFS: "+str(psi_lcfs))
        Aphi = eval_field('A', R, phi, Z, coord='phi', sim=sims[0], file_name=file_name[0], time=sims[0].timeslice)
        
        psifield = R_ave*Aphi
        for i,ax in enumerate(axs):
            cont = ax.contour(R_ave, Z_ave, np.average(psifield,0),[psi_lcfs],colors='magenta',linewidths=lcfslw,zorder=10)
    
    
    plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    #if coord=='vector':
    #    fig.subplots_adjust(top=0.2)
    plt.show()
    
    if save == True:
        timestr = 't'
        if diff == True:
            timestr = timestr + str(time[0]) + '-' + str(time[1])
        else:
            timestr = timestr + str(time[0])
        if linear == True:
            fieldstr = 'delta_' + str(field)
        else:
            fieldstr = str(field)
        if savedir != None:
            fieldstr = savedir + fieldstr
        if n==None:
            nout=sims[0].ntor
        else:
            nout=n
        
        plt.savefig(fieldstr + '_' + timestr + '_n'+"{:d}".format(nout)+'.png', format='png',dpi=900,bbox_inches='tight')
    
    return
