#!/usr/bin/env python3
#
# Coded on 08/27/2019 by:
# Andreas Kleiner:    akleiner@pppl.gov
# Ralf Mackenbach:    rmackenb@pppl.gov
# Chris Smiet    :    csmiet@pppl.gov

import math
import numpy as np
import matplotlib.axes as mplax
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import fpy
import m3dc1.fpylib as fpyl
from m3dc1.plot_coils import plot_coils


def plot_mesh(elms=None,time=None,filename='C1.h5',sim=None,boundary=False,coils=False,ax=None,fignum=None,meshcol='C0',zoom=False,pub=False,quiet=False,phys=False,phi=0.,save=False):
    """
    plot_mesh: Creates a plot of the mesh from a M3D-C1 time slice.
    plot_mesh can take the mesh object as input. This is better for large
    meshes since it can take a while to calculate the mesh connectivity.
    Example:
    sim = fpy.sim_data(time=0)
    elms = sim.get_mesh(time=0,filename=filename)
    m.plot_mesh(elms=elms)

    Arguments:

    **elms**
    mesh object

    **time**
    If no mesh object is provided, read mesh at this time slice

    **filename**
    If no mesh object is provided, read file with this name

    **boundary**
    If True, only boundary (resistive wall and vacuum boundary) will be ploted

    **ax**
    Matplotlib axes object. If not None, 

    **fignum**
    Figure number for mesh plot

    **meshcol**
    Color of mesh lines

    **pub**
    If True, plot will be formatted for publication
    """
    if not isinstance(elms,fpy.sim_data.mesh):
        if not isinstance(sim,fpy.sim_data):
            sim = fpy.sim_data(filename)
        elms = sim.get_mesh(time=time,quiet=quiet)
    
    mesh = elms.elements
    version = elms.version
    nplanes = elms.nplanes
    
    meshshape = mesh.shape
    
    nelms = meshshape[0]
    
    threed = 1 if meshshape[1]>8 else 0
    
    boundary = boundary and (version >= 3)
    
    if phys:
        rst = eval_field('rst', mesh[:,4], phi*np.ones_like(mesh[:,4]), mesh[:,5], sim=simplot, file_name=file_name, time=time)
        zst = eval_field('zst', mesh[:,4], phi*np.ones_like(mesh[:,4]), mesh[:,5], sim=simplot, file_name=file_name, time=time)
        minr = np.amin(rst)
        maxr = np.amax(rst)
        minz = np.amin(zst)
        maxz = np.amax(zst)
    else:
        minr = np.amin(mesh[:,4])
        maxr = np.amax(mesh[:,4])
        minz = np.amin(mesh[:,5])
        maxz = np.amax(mesh[:,5])


    # Create figure and plot mesh points only. If you want them to be visible,
    # increase marker size.
    
    # Set font sizes and plot style parameters
    if pub:
        axlblfs = 20
        ticklblfs = 18
        linew = 2
        bdlw = 1
    else:
        axlblfs = 12
        ticklblfs = 12
        linew = 1
        bdlw = 1
    
    if not isinstance(ax, (np.ndarray, mplax._axes.Axes)):
        if (fignum is None) or (fignum is not None and not plt.fignum_exists(fignum)):
            fig, ax = plt.subplots(num=fignum)
            fig.set_figheight(8)
        else:
            ax = plt.gca()
        if phys:
            plt.plot(rst,zst,lw=0,marker='.',ms=0)
        else:
            plt.plot(mesh[:,4],mesh[:,5],lw=0,marker='.',ms=0)
        if coils:
            plot_coils(ax=ax)
        plt.grid(True)
        ax.set_aspect('equal',adjustable='box')
        plt.xlabel(r'$R$',fontsize=axlblfs)
        plt.ylabel(r'$Z$',fontsize=axlblfs)
        ax.tick_params(labelsize=ticklblfs)
        ax.set_xlim(fpyl.get_axlim(minr,'min'),fpyl.get_axlim(maxr,'max'))
        ax.set_ylim(fpyl.get_axlim(minz,'min'),fpyl.get_axlim(maxz,'max'))
        plt.tight_layout()
    
    if version==0:
        fpyl.printwarn("WARNING: xzero and zzero not set!")
    else:
        xzero = 0.
        zzero = 0.
    
    czone = np.asarray([0, 1, 1, 2, 2, 3, 3]) # Color zones: Needed to plot lines of walls and boundary in correct color
    nczone = czone.size
    colors=['C0','C1','C2','C3'] #Colors of the default matplotlib color cycle
    
    # Empty arrays that will be filled with start and end points of the mesh lines
    plot1x = []
    plot1y = []
    plot2x = []
    plot2y = []
    plot3x = []
    plot3y = []
    plot4x = []
    plot4y = []
    
    #The following loop 
    for i in range(nelms):
        if (threed == 1) and (i >= (nelms/nplanes)):
            break
        
        # Mesh information
        a = mesh[i,0]
        b = mesh[i,1]
        c = mesh[i,2]
        t = mesh[i,3]
        x = mesh[i,4]
        y = mesh[i,5]
        bound = math.floor(mesh[i,6])

        # Start and end points of the mesh lines
        p1 = np.asarray([x, y])
        p2 = p1 + np.asarray([(b+a) * math.cos(t), (b+a) * math.sin(t)])
        p3 = p1 + np.asarray([b * math.cos(t) - c * math.sin(t), b * math.sin(t) + c * math.cos(t)])
        delta = 0.0
        if phys: # calculate physical coordinates
            p1 = np.asarray([eval_field('rst', p1[0], phi*np.ones_like(p1[0]), p1[1], sim=simplot, file_name=file_name, time=time),\
                             eval_field('zst', p1[0], phi*np.ones_like(p1[0]), p1[1], sim=simplot, file_name=file_name, time=time)])
            p2 = np.asarray([eval_field('rst', p2[0], phi*np.ones_like(p2[0]), p2[1], sim=simplot, file_name=file_name, time=time),\
                             eval_field('zst', p2[0], phi*np.ones_like(p2[0]), p2[1], sim=simplot, file_name=file_name, time=time)])
            p3 = np.asarray([eval_field('rst', p3[0], phi*np.ones_like(p3[0]), p3[1], sim=simplot, file_name=file_name, time=time),\
                             eval_field('zst', p3[0], phi*np.ones_like(p3[0]), p3[1], sim=simplot, file_name=file_name, time=time)])
        q1 = (1.-2.*delta)*p1 + delta*p2 + delta*p3
        q2 = delta*p1 + (1.-2.*delta)*p2 + delta*p3
        q3 = delta*p1 + delta*p2 + (1.-2.*delta)*p3
        
        # Identify lines that are part of the wall or boundary. These will be plotted in different colors.
        pp=bound if boundary else 7
        
        if((pp & 1) == 1):
            if((bound & 1) == 1):
                izone = (bound & 120)/2**3 + 1
                col = czone[int(izone % nczone)]
                if col == 1:
                    plot2x.append([q1[0],q2[0]])
                    plot2y.append([q1[1],q2[1]])
                elif col == 2:
                    plot3x.append([q1[0],q2[0]])
                    plot3y.append([q1[1],q2[1]])
                elif col ==3:
                    plot4x.append([q1[0],q2[0]])
                    plot4y.append([q1[1],q2[1]])
            else:
                plot1x.append([p1[0],p2[0]])
                plot1y.append([p1[1],p2[1]])
        
        if((pp & 2) == 2):
            if((bound & 2) == 2):
                izone = (bound & 1920)/2**7 + 1
                col = czone[int(izone % nczone)]
                if col == 1:
                    plot2x.append([q2[0],q3[0]])
                    plot2y.append([q2[1],q3[1]])
                elif col == 2:
                    plot3x.append([q2[0],q3[0]])
                    plot3y.append([q2[1],q3[1]])
                elif col ==3:
                    plot4x.append([q2[0],q3[0]])
                    plot4y.append([q2[1],q3[1]])
            else:
                plot1x.append([p2[0],p3[0]])
                plot1y.append([p2[1],p3[1]])
        
        if((pp & 4) == 4):
            if((bound & 4) == 4):
                izone = (bound & 30720)/2**11 + 1
                col = czone[int(izone % nczone)]
                if col == 1:
                    plot2x.append([q3[0],q1[0]])
                    plot2y.append([q3[1],q1[1]])
                elif col == 2:
                    plot3x.append([q3[0],q1[0]])
                    plot3y.append([q3[1],q1[1]])
                elif col ==3:
                    plot4x.append([q3[0],q1[0]])
                    plot4y.append([q3[1],q1[1]])
            else:
                plot1x.append([p3[0],p1[0]])
                plot1y.append([p3[1],p1[1]])
    # Plot mesh lines. Plotting all lines of a certain color at once is the fasted way.
    # Calls to any plotting function are slow and should thus be minimized.
    
    #Check if an axis object or a an array of axis objects was passed to this routine and plot in each axis.
    axarray = np.atleast_1d(ax)
    
    for ax in axarray:
        if not boundary:
            if len(plot1x) > 0:
                pltreg = ax.add_collection(LineCollection(np.stack((plot1x,plot1y), axis=2),linewidths=0.25, colors=meshcol,zorder=4))
            else:
                pltreg = None
        if len(plot2x) > 0:
            pltwin = ax.add_collection(LineCollection(np.stack((plot2x,plot2y), axis=2),linewidths=bdlw, colors='r',zorder=5))
            #pltwin = ax.add_collection(LineCollection(np.stack((plot2x,plot2y), axis=2),linewidths=bdlw, colors=r'#1f77b4', linestyles='--', zorder=5)) #For CLT paper
        else:
            pltwin = None
        if len(plot3x) > 0:
            pltwout = ax.add_collection(LineCollection(np.stack((plot3x,plot3y), axis=2),linewidths=bdlw, colors='g',zorder=5))
        else:
            pltwout = None
        if len(plot4x) > 0:
            pltbd = ax.add_collection(LineCollection(np.stack((plot4x,plot4y), axis=2),linewidths=bdlw, colors='c',zorder=5))
        else:
            pltbd = None
    
    if zoom:
        axins = zoomed_inset_axes(ax, 2, loc=1) # zoom = 6
        # sub region of the original image
        x1, x2, y1, y2 = 1.2, 1.75, -0.4, 0.4
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)

        plt.xticks(visible=False)
        plt.yticks(visible=False)

        # draw a bbox of the region of the inset axes in the parent axes and
        # connecting lines between the bbox and the inset axes area
        mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    
    if save:
        plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
        plt.savefig('mesh.png',bbox_inches='tight',dpi=250)
        plt.close()
    else:
        if not isinstance(ax, (np.ndarray, mplax._axes.Axes)):
            plt.show()
    
    if boundary:
        return pltwin, pltwout, pltbd
    else:
        return pltreg, pltwin, pltwout, pltbd
