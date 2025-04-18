#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11/23/2022

@author: Andreas Kleiner
"""

import numpy as np
#import sys
import matplotlib.pyplot as plt
from matplotlib import colors
#from matplotlib import path
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
#from distutils.version import StrictVersion
from sys import modules as sysmod
from packaging import version
import m3dc1.fpylib as fpyl
from m3dc1.get_field import get_field_vs_phi
from m3dc1.plot_mesh import plot_mesh
#rc('text', usetex=True)




def plot_field_vs_phi(field, cutr=None, cutz=None, coord='scalar', row=1, sim=None, filename='C1.h5', time=None, linear=False,
               diff=False, phi_res=100, res=250, mesh=False, bound=False, planes=False, units='mks',cont_levels=100,
               prange=None, cmap='viridis', cmap_midpt=None, quiet=False, save=False, savedir=None, pub=False,
               titlestr=None, showtitle=True, shortlbl=False, ntor=None, phys=False, export=False,txtname=None,):
    """
    Plots a field as a function of either R or Z and toroidal angle phi.
    
    Arguments:

    **field**
    The field that is to be plotted, i.e. 'B', 'j', etc.

    **cutr**
    Value of R, where field is being plotted as a function of Z and phi.
    cutz has to be None.

    **cutz**
    Value of Z, where field is being plotted as a function of R and phi.
    cutr has to be None.

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

    **phi_res**
    Number of points in direction of phi.

    **res**
    Resolution in direction of R or Z.

    **linear**
    Plot the linear part of the field (so the equilibrium is subtracted).
    True/False

    **diff**
    Plot the difference of two fields. 
    This could be the difference of two files (filename=['a/C1.h5','b/C1.h5']),
    or the difference between two time-slices (time=[t1,t2])
    If list for both time and filename are given file1 will be evaluated at time1,
    and file2 at time2

    **units**
    Units in which the field will be plotted.

    **mesh**
    True/False
    Overlay mesh on top of plot.

    **bound**
    True/False
    Only plot boundary. Only works when mesh is true

    **planes**
    True/False
    Overlay location of toroidal planes on top of plot.

    **cont_levels**
    Number of contour levels.

    **cmap**
    Color map to use for the contour plot, e.g. viridis, seismic, jet, nipy_spectral

    **cmap_midpt**
    If not None, set midpoint of colormap to cmap_midpt.

    **titlestr**
    Plot title. If None, a default title will be generated.

    **showtitle**
    True/False. Show plot title.

    **shortlbl**
    True/False. If True, axis labels will be shortened.

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

    **quiet**
    If True, suppress output to terminal.
    """
    
    sim, time, mesh_ob, R, phi_list, Z, R_mesh, Z_mesh, R_ave, Z_ave, phi_ave, field1_ave = get_field_vs_phi(field=field, coord=coord, row=row, sim=sim, filename=filename, time=time,
                                                                                                    cutr=cutr, cutz=cutz, linear=linear,  diff=diff, phi_res=phi_res,
                                                                                                    units=units, res=res, quiet=quiet, phys=phys)

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
    else:
        axlblfs = None
        titlefs = None
        cbarlblfs = None
        cbarticklblfs = None
        ticklblfs = None
    
    if coord not in ['vector', 'tensor']:
        fig, axs = plt.subplots(1, 1,figsize=(6.4,4.8))
        #fig.set_figheight(7)
        axs = np.asarray([axs])
    else:
        fig, axs = plt.subplots(1, 3, sharey=True,figsize=(14,7))
        comp = ['R',r'\phi','Z']
    
    if titlestr is None:
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
        if mesh and bound:#Make sure that whole mesh is plotted. If both are true, plot_mesh only plots boundary.
            bound = False
        meshplt = plot_mesh(mesh_ob,boundary=bound,ax=axs,meshcol='C1',pub=pub,phys=phys)
    
    if planes:
        tp_file = open("plane_positions", "r")
        data = tp_file.readlines()
        tp_file.close()
        temp = [d.replace('\n', '') for d in data]
        tor_planes = np.asarray(temp).astype(float)
        for tp in tor_planes:
            if cutr is not None and cutz is None:
                plt.axvline(x=tp,c='C1')
            elif cutz is not None and cutr is None:
                plt.axhline(y=tp,c='C1')
    
    for i,ax in enumerate(axs):
        if cmap_midpt is not None:
            field1_ave_clean = field1_ave[i][np.logical_not(np.isnan(field1_ave[i]))]
            if StrictVersion(sys.modules[plt.__package__].__version__) < StrictVersion('3.2'):#DivergingNorm became TwoSlopeNorm in matplotlib version >=3.2
                norm = colors.DivergingNorm(vmin=np.amin(field1_ave_clean), vcenter=cmap_midpt, vmax=np.amax(field1_ave_clean))
            else:
                norm = colors.TwoSlopeNorm(vmin=np.amin(field1_ave_clean), vcenter=cmap_midpt, vmax=np.amax(field1_ave_clean))
            if cutr is not None and cutz is None:
                cont = ax.contourf(phi_ave, Z_ave, field1_ave[i],cont_levels, cmap=cmap,norm=norm)
            elif cutz is not None and cutr is None:
                cont = ax.contourf(R_ave, phi_ave, field1_ave[i],cont_levels, cmap=cmap,norm=norm)
        else:
            if isinstance(prange,(tuple,list)):
                #if StrictVersion(sys.modules[plt.__package__].__version__) < StrictVersion('3.2'):#DivergingNorm became TwoSlopeNorm in matplotlib version >=3.2
                #    norm = colors.DivergingNorm(vmin=prange[0], vcenter=(prange[1]+prange[0])/2, vmax=prange[1])
                #else:
                #    norm = colors.TwoSlopeNorm(vmin=prange[0], vcenter=(prange[1]+prange[0])/2, vmax=prange[1])
                if cutr is not None and cutz is None:
                    cont = ax.contourf(phi_ave, Z_ave, field1_ave[i],600, cmap=cmap,vmin=prange[0], vmax=prange[1])
                elif cutz is not None and cutr is None:
                    cont = ax.contourf(R_ave, phi_ave, field1_ave[i],600, cmap=cmap,vmin=prange[0], vmax=prange[1])
            else:
                if cutr is not None and cutz is None:
                    cont = ax.contourf(phi_ave, Z_ave, field1_ave[i],cont_levels, cmap=cmap)
                elif cutz is not None and cutr is None:
                    #print(R_ave.shape,phi_ave.shape,field1_ave[i].shape)
                    cont = ax.contourf(R_ave, phi_ave, field1_ave[i],cont_levels, cmap=cmap)
        # Set and format axes limits and labels
        if phys:
            ax.set_xlim([fpyl.get_axlim(np.amin(R_mesh),'min',0.1),fpyl.get_axlim(np.amax(R_mesh),'max',0.1)])
            ax.set_ylim([fpyl.get_axlim(np.amin(Z_mesh),'min',0.1),fpyl.get_axlim(np.amax(Z_mesh),'max',0.1)])
        #else:
        #    if cutr is not None and cutz is None:
        #        ax.set_xlim([fpyl.get_axlim(np.amin(R_ave),'min',0.1),fpyl.get_axlim(np.amax(R_ave),'max',0.1)])
        #        ax.set_ylim([fpyl.get_axlim(np.amin(Z_ave),'min',0.1),fpyl.get_axlim(np.amax(Z_ave),'max',0.1)])
        #    elif cutz is not None and cutr is None:
        #        ax.set_xlim([fpyl.get_axlim(np.amin(R_ave),'min',0.1),fpyl.get_axlim(np.amax(R_ave),'max',0.1)])
        #        ax.set_ylim([fpyl.get_axlim(np.amin(Z_ave),'min',0.1),fpyl.get_axlim(np.amax(Z_ave),'max',0.1)])
        if cutr is not None and cutz is None:
            xlbl = r'Toroidal angle $\phi$ in rad'
            ylbl = r'$Z/m$'
        elif cutz is not None and cutr is None:
            xlbl = r'$R/m$'
            ylbl = r'Toroidal angle $\phi$ in rad'
        ax.set_xlabel(xlbl,fontsize=axlblfs)
        ax.set_ylabel(ylbl,fontsize=axlblfs)
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
        if version.parse(sysmod[plt.__package__].__version__) >= version.parse('3.10.0'):#Attribute collections was removed in matplotlib 3.10
            cont.set_edgecolor("face")
        else:
            for c in cont.collections:
                c.set_edgecolor("face")

    plt.rcParams["axes.axisbelow"] = False #Allows mesh to be on top of contour plot. This option conflicts with zorder (matplotlib bug).
    
    
    
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
    
    if export:
        #print(phi_ave.shape, Z_ave.shape, field1_ave.shape)
        if cutr is not None and cutz is None:
            np.savetxt(txtname+'_phi',np.vstack(phi_ave),delimiter='   ')
            np.savetxt(txtname+'_Z',np.vstack(Z_ave),delimiter='   ')
            np.savetxt(txtname+'_field',np.vstack(field1_ave),delimiter='   ')
        elif cutz is not None and cutr is None:
            np.savetxt(txtname+'_phi',np.vstack(phi_ave),delimiter='   ')
            np.savetxt(txtname+'_Z',np.vstack(Z_ave),delimiter='   ')
            np.savetxt(txtname+'_field',np.vstack(field1_ave),delimiter='   ')
    
    return
