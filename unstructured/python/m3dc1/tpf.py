#!/usr/bin/env python3
#
# Coded on 02/24/2023 by:
# Andreas Kleiner:    akleiner@pppl.gov

import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
from scipy.interpolate import RectBivariateSpline
from scipy.integrate import dblquad
import fpy
import m3dc1.fpylib as fpyl
from m3dc1.get_field import get_field
from m3dc1.get_time_of_slice import get_time_of_slice
from m3dc1.get_timetrace import get_timetrace
from m3dc1.read_h5 import readParameter



def tpf(field, method='sum', coord='scalar', row=1, sim=None, filename='C1.h5', time=None, phi_res=24, linear=False,
               diff=False, units='mks',millisec=False,res=250, quiet=False, pub=False, titlestr=None, showtitle=True, phys=False):
    """
    Calculate toroidal peaking factor (TPF) for a given field.
    
    Arguments:

    **field**
    The field that is to be plotted, i.e. 'B', 'j', etc.

    **method**
    Method for integration. Options are 'sum' and 'spline'.
    'sum' integrates by summing up array elements/
    'spline' uses the integral method of RectBivariateSpline.

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

    **res**
    Resolution in R and Z direction.

    **phi_res**
    Number of points in direction of the toroidal angle.

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

    **millisec**
    True/False. If True and units='mks' plot will be in terms of milliseconds, instead of seconds.

    **titlestr**
    Plot title. If None, a default title will be generated.

    **showtitle**
    True/False. Show plot title.

    **pub**
    If True, format figure for publication (larger labels and thicker lines)

    **phys**
    Use True for plotting in physical (stellarator) geometry.

    **quiet**
    If True, suppress output to terminal.
    """
    
    if isinstance(sim,fpy.sim_data):
        time = sim.timeslice
    time_at_ts = get_time_of_slice(time,units=units,millisec=millisec)
    print('Calculating TPF at time: ' + str(time_at_ts))
    
    tor_av=1
    tpf_list = []
    phi_planes = np.linspace(0,2*np.pi,phi_res,endpoint=False)
    
    for phi in phi_planes:
        _, time, mesh_ob, R, phi_list, Z, R_mesh, Z_mesh, R_ave, Z_ave, field1_ave = get_field(field=field, coord=coord, row=row, sim=sim, filename=filename, time=time, phi=phi, linear=linear,
                diff=diff, tor_av=tor_av, units=units, res=res, quiet=quiet, phys=phys)
        
        fieldlabel,unitlabel = fpyl.get_fieldlabel(units,field)
        cbarlbl = fieldlabel + ' (' + unitlabel + ')'
        
        R_1d = np.unique(R)
        Z_1d = np.unique(Z)
        dR = R_1d[1]-R_1d[0]
        dZ = Z_1d[1]-Z_1d[0]
        
        field_in_plane = field1_ave[0,:,:]
        
        if method == 'sum':
            A=dR*dZ
            integral = A*np.nansum(field_in_plane)
        elif method == 'spline':
            spline = RectBivariateSpline(R_1d,Z_1d,np.nan_to_num(field_in_plane),kx=3,ky=3)
            integral = spline.integral(np.amin(R_1d),np.amax(R_1d),np.amin(Z_1d),np.amax(Z_1d))  
        else:
            fpyl.printerr('ERROR: method not recognized!')
            return
        
        print(phi,integral)
        tpf_list.append(integral)
    
    
    pf = np.amax(np.abs(tpf_list))/np.mean(tpf_list)
    print(pf)
    
    
    #def fun(x,y):
    #    val = spline(x,y)
    #    if np.isnan(val):
    #        return 0.0
    #    else:
    #        print(val)
    #        return val
    #integral = dblquad(fun,np.amin(R_1d),np.amax(R_1d),np.amin(Z_1d),np.amax(Z_1d))
    
    
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
    
    plt.figure()
    ax = plt.gca()
    plt.plot(phi_planes,tpf_list)
    ax.set_xlabel(r'$\phi$',fontsize=axlblfs)
    ax.set_ylabel(r'$f(\phi) = \int P_{rad}(R,Z,\phi) dR dZ$',fontsize=axlblfs)
    ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
    plt.grid(True)
    
    return time_at_ts,pf








def tpf_vs_t(field, method='sum', coord='scalar', row=1, filename='C1.h5', phi_res=24, ts_range=[],linear=False,
               diff=False, units='mks',millisec=False,res=250, quiet=False, pub=False, titlestr=None,
               showtitle=True, phys=False,write_to_file=True,filepath=''):
    """
    Calculate toroidal peaking factor (TPF) for a given field.
    
    Arguments:

    **field**
    The field that is to be plotted, i.e. 'B', 'j', etc.

    **method**
    Method for integration. Options are 'sum' and 'spline'.
    'sum' integrates by summing up array elements/
    'spline' uses the integral method of RectBivariateSpline.

    **coord**
    The chosen part of a vector field to be plotted, options are:
    'phi', 'R', 'Z', 'poloidal', 'radial', 'scalar', 'vector', 'tensor'.
    Scalar is reserved for scalar fields.
    Poloidal takes the magnetic axis at time=0 as the origin, and
    defines anti-clockwise as the positive direction.
    Radial also takes the magnetic axis as the origin.

    **row**
    For tensor fields only, row of tensor to plot. Possibilities: 1,2,3

    **filename**
    File name which will be read, i.e. "../../C1.h5"
    Can also be a list of two filepaths when used for diff

    **res**
    Resolution in R and Z direction.

    **phi_res**
    Number of points in direction of the toroidal angle.

    **ts_range**
    Range of time slices. Can be a tupel/list of two or three numbers:
    (start, end) or (start, end, step).

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

    **millisec**
    True/False. If True and units='mks' plot will be in terms of milliseconds, instead of seconds.

    **titlestr**
    Plot title. If None, a default title will be generated.

    **showtitle**
    True/False. Show plot title.

    **pub**
    If True, format figure for publication (larger labels and thicker lines)

    **phys**
    Use True for plotting in physical (stellarator) geometry.

    **write_to_file**
    If True, write TPF as function of time to ASCII file.

    **filepath**
    Directory where TPF output is written. Default is current working directory.

    **quiet**
    If True, suppress output to terminal.
    """
    
    # Get number of time slices
    sim0 = fpy.sim_data(filename=filename,time=-1)
    nts = sim0.ntime
    
    time_at_ts = []
    tpf_list = []
    
    #Set list of time slices at which to calculate TPF
    if len(ts_range)>0:
        try:
            ts_list = range(ts_range[0],ts_range[1]+1,ts_range[2])
        except:
            ts_list = range(ts_range[0],ts_range[1]+1)
    else:
        ts_list = range(nts)
    
    # Open file to write TPF into. Because TPF calculation is time consuming, this is a safe guard to store already
    # calculated TPFs in case calculation fails before completing.
    if write_to_file:
        outfile = filepath + 'TPF_'+field+'_t_'+str(ts_list[0])+'-'+str(ts_list[-1])+'.dat'
        f = open(outfile, 'w')
        f.close()
    
    #Calculate TPF at each chosen time slice
    for ts in ts_list:
        #try:
        sim = fpy.sim_data(filename,time=ts) 
        t,pf = tpf(field=field, method=method, coord=coord, row=row, sim=sim, phi_res=phi_res, linear=linear,diff=diff, units=units,millisec=millisec,res=res, quiet=quiet, pub=pub, titlestr=titlestr, showtitle=showtitle, phys=phys)
        #except:
        #    continue
        
        time_at_ts.append(t)
        tpf_list.append(pf)
        if write_to_file:
            f = open(outfile, 'a')
            f.write(str(ts) + '    ' + str(t) + '    ' + str(pf) + '\n')
            f.close()
    
    #if write_to_file:
    #    f.close()
    
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
    
    plt.figure()
    ax = plt.gca()
    plt.plot(time_at_ts,tpf_list)
    
    if units=='mks':
        if millisec:
            ax.set_xlabel(r'time $(ms)$',fontsize=axlblfs)
        else:
            ax.set_xlabel(r'time $(s)$',fontsize=axlblfs)
    else:
        ax.set_xlabel(r'time $(\tau_A)$',fontsize=axlblfs)
    
    ax.set_ylabel(r'TPF ('+field+')',fontsize=axlblfs)
    ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
    plt.grid(True)
    plt.tight_layout()

    return


def multiply_tpf_time_trace(tpf_file,trace,units='mks',millisec=False,fac=1,growth=False,renorm=False,makeplot=True,
                            fignum=None,figsize=None,leglbl=None,show_legend=False,quiet=False,pub=False,in_plot_txt=None,
                            export=False,txtname=None):

    #Read file with TPF vs time
    f = open(tpf_file, 'r')
    data_str = f.readlines()
    f.close()
    tpfdata = np.asarray([line.split() for line in data_str],dtype=float)
    
    _,time_trace,label,unitlabel = get_timetrace(trace=trace,units=units,growth=growth,renorm=renorm,fac=fac,quiet=quiet)
    
    values = []
    #print(tpfdata)
    if units=='mks':
        time = tpfdata[:,1]
        if millisec:
            time = time*1000
    elif units=='m3dc1':
        time = []
    
    for i in range(len(tpfdata)):
        ts_filename = 'time_'+str(int(tpfdata[i,0])).zfill(3)+'.h5'
        ts = readParameter('ntimestep',fname=ts_filename)
        #print(time_trace)
        #print(ts)
        values.append(tpfdata[i,2] * time_trace[ts])
        if units=='m3dc1':
            time.append(readParameter('time',fname=ts_filename))
    
    if makeplot:
        # Set font sizes and plot style parameters
        if pub:
            axlblfs = 18
            titlefs = 16
            ticklblfs = 16
            linew = 2
            legfs = 12
            inplttxtfs = 20
        else:
            axlblfs = None
            titlefs = None
            ticklblfs = None
            linew = 1
            legfs = None
            inplttxtfs = 16
        
        plt.figure(num=fignum,figsize=figsize)
        
        if units=='mks':
            #time = unit_conv(time,filename=filename,time=1)
            if millisec:
                plt.xlabel(r'time $(ms)$',fontsize=axlblfs)
            else:
                plt.xlabel(r'time $(s)$',fontsize=axlblfs)
        else:
            plt.xlabel(r'time $(\tau_A)$',fontsize=axlblfs)
        
        if trace == 'prad':
            label = r'$P_{rad}$'
        ylbl = r'TPF $\cdot$ ' + label+' ('+unitlabel+')'
        plt.ylabel(ylbl,fontsize=axlblfs)
        
        ax = plt.gca()
        ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
        
        if in_plot_txt is not None:
            plt.text(0.05, 0.95,in_plot_txt, ha='left', va='top', transform=ax.transAxes,fontsize=inplttxtfs)
        
        plt.plot(time,values,lw=linew,label=leglbl)
        plt.grid(True)
        if show_legend:
            ax.legend(loc=0,fontsize=legfs)
        plt.tight_layout()
    
    if export:
        np.savetxt(txtname,np.vstack((time,values)).transpose(),delimiter='   ',fmt='%f')
    return
