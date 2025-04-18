#!/usr/bin/env python3
#
# Coded on 11/26/2019 by:
# Andreas Kleiner:    akleiner@pppl.gov
# Ralf Mackenbach:    rmackenb@pppl.gov
# Chris Smiet    :    csmiet@pppl.gov

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
#from matplotlib import path
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
#from distutils.version import StrictVersion
from sys import modules as sysmod
from packaging import version
import m3dc1.fpylib as fpyl
from m3dc1.get_field import get_field
from m3dc1.eval_field import eval_field
from m3dc1.plot_mesh import plot_mesh
from m3dc1.plot_coils import plot_coils
#rc('text', usetex=True)




def plot_field(field, coord='scalar', row=1, sim=None, filename='C1.h5', time=None, phi=0, linear=False,
               diff=False, scale=1.0, tor_av=1, mesh=False, bound=False, lcfs=False, coils=False, units='mks',points=250,
               fac=1.,prange=None, cmap='viridis', cmap_midpt=None, clvl=100, quiet=False, save=False, savedir=None,
               pub=False, titlestr=None, showtitle=True, shortlbl=False, ntor=None, phys=False,export=False,txtname=None, R_range=None, Z_range=None):
    """
    Plots a field in the R,Z plane.
    
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

    **units**
    Units in which the field will be plotted.

    **mesh**
    True/False
    Overlay mesh on top of plot.

    **bound**
    True/False
    Plot boundary (wall and domain boundary) without the mesh itself.

    **lcfs**
    True/False
    Overlay last closed flux surface on top of the plot.

    **coils**
    True/False
    Overlay coils from coil.dat file on top of the plot.

    **points**
    Resolution in R and Z direction.

    **cmap**
    Color map to use for the contour plot, e.g. viridis, seismic, jet, nipy_spectral

    **clvl**
    Number of contour levels.

    **prange**
    Plot range to renormalize color map. Should be a list of length 2, defining the
    minimum and maximum values of the field. Values outside of this range will have
    the same color as the endpoints of the colormap.

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
    Use True for plotting in physical (stellarator) geometry.

    **quiet**
    If True, suppress output to terminal.
    """
    
    sim, time, mesh_ob, R, phi_list, Z, R_mesh, Z_mesh, R_ave, Z_ave, field1_ave = get_field(field=field, coord=coord, row=row, sim=sim, filename=filename, time=time, phi=phi, linear=linear,
               diff=diff, scale=scale, tor_av=tor_av, units=units, res=points, quiet=quiet, phys=phys, R_range=R_range, Z_range=Z_range)
    

    fieldlabel,unitlabel = fpyl.get_fieldlabel(units,field,shortlbl=shortlbl,fac=fac)
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
    
    if coils:
        plot_coils(ax=axs)
    if mesh or bound:
        if mesh and bound:#Make sure that whole mesh is plotted. If both are true, plot_mesh only plots boundary.
            bound = False
        meshplt = plot_mesh(mesh_ob,boundary=bound,ax=axs,meshcol='C1',pub=pub,phys=phys,quiet=quiet)
    
    if fac != 1.:
        field1_ave = field1_ave*fac
    
    for i,ax in enumerate(axs):
        if cmap_midpt is not None:
            field1_ave_clean = field1_ave[i][np.logical_not(np.isnan(field1_ave[i]))]
            vmin = np.amin(field1_ave_clean)
            vmid = cmap_midpt
            vmax = np.amax(field1_ave_clean)
            if not np.all(np.diff([vmin,vmid,vmax]) >= 0): #Check if values are in ascending order
                vmid = vmin + (vmax-vmin)/2
                fpyl.printwarn('WARNING: cmap_midpt is either too low or too high! Resetting cmap_midpt to vmin + (vmax-vmin)/2.')
            #if StrictVersion(sys.modules[plt.__package__].__version__) < StrictVersion('3.2'):#DivergingNorm became TwoSlopeNorm in matplotlib version >=3.2
            if version.parse(sysmod[plt.__package__].__version__) >= version.parse('3.2'):#DivergingNorm became TwoSlopeNorm in matplotlib version >=3.2
                norm = colors.DivergingNorm(vmin=vmin, vcenter=vmid, vmax=vmax)
            else:
                norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=vmid, vmax=vmax)
            cont = ax.contourf(R_ave, Z_ave, field1_ave[i],clvl, cmap=cmap,norm=norm)
        else:
            if isinstance(prange,(tuple,list)):
                #version.parse(sysmod[plt.__package__].__version__) >= version.parse('3.10.0'):#DivergingNorm became TwoSlopeNorm in matplotlib version >=3.2
                #    norm = colors.DivergingNorm(vmin=prange[0], vcenter=(prange[1]+prange[0])/2, vmax=prange[1])
                #else:
                #    norm = colors.TwoSlopeNorm(vmin=prange[0], vcenter=(prange[1]+prange[0])/2, vmax=prange[1])
                cont = ax.contourf(R_ave, Z_ave, field1_ave[i],np.linspace(prange[0],prange[1],clvl), cmap=cmap,vmin=prange[0], vmax=prange[1])
            else:
                cont = ax.contourf(R_ave, Z_ave, field1_ave[i],clvl, cmap=cmap)
        # Set and format axes limits and labels
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
        if version.parse(sysmod[plt.__package__].__version__) >= version.parse('3.10.0'):#Attribute collections was removed in matplotlib 3.10
            cont.set_edgecolor("face")
        else:
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
        psi_lcfs = sim[lcfs_ts_ind].get_time_trace('psi_lcfs').values[0] #ToDo: BUG: get values of psi at actual time
        if not quiet:
            print("Psi at LCFS: "+str(psi_lcfs))
        
        #Aphi = eval_field('A', R, phi, Z, coord='phi', sim=sim[lcfs_ts_ind], time=time[lcfs_ts_ind])
        #psifield = R_ave*Aphi
        psifield = eval_field('psi', R, phi_list, Z, coord='scalar', sim=sim[lcfs_ts_ind], filename=sim[lcfs_ts_ind].filename, time=time[lcfs_ts_ind])
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
    
    if export:
        np.savetxt(txtname+'_R',np.vstack(R_ave),delimiter='   ')
        np.savetxt(txtname+'_Z',np.vstack(Z_ave),delimiter='   ')
        np.savetxt(txtname+'_field',np.vstack(field1_ave),delimiter='   ')
    
    return
