#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 14:24:12 2019

@author: Andreas Kleiner
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
#from scipy.fftpack import fft, ifft, fftshift, fftfreq

import fpy
import m3dc1.fpylib as fpyl
from m3dc1.read_h5 import readParameter

rc('text', usetex=True)
plt.rcParams.update({'figure.max_open_warning': 40})




def get_timetrace(trace,sim=None,filename='C1.h5',units='m3dc1',ipellet=0,diff=False,
                  growth=False,renorm=False,quiet=False,returnas='tuple',unitlabel=None,fac=1):
    """
    Read a time trace directly from an hdf5 file. This function does not use fusion-io.
    
    Arguments:

    **trace**
    Name of trace (scalar)

    **sim**
    fpy simulation object.

    **filename**
    Name or path to C1.h5 file to read

    **units**
    The units in which the time trace will be returned

    **growth**
    If True, return growth rate of trace. If false, return trace

    **renorm**
    Remove spikes due to renormalization in linear runs that happens when
    the kinetic energy gets too large.

    **quiet**
    If True, do not print renormalization times to screen.

    **returnas**
    Determines how time trace is being returned.
    'tuple': Returns a tuple of (time, values, label, unitlabel)
    'time_trace': returns time trace as object inside a tuple
                  (fpy.sim_data.time_trace(values,time=time), label, unitlabel)

    **unitlabel**
    Deprecated.

    **fac**
    Scale factor for time trace. Returned time trace values will be multiplied
    by fac. If fac equals 1.0E-3, 1.0E-6 or 1.0E-9, a 'k', 'M' or 'G' will be
    prepended to the unitlabel to reflect the correct order of magnitude.
    """
    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename=filename)
    constants = sim.get_constants()
    itor    = constants.itor
    version = constants.version
    gamma   = constants.gamma

    # Direct transformation of one name to another
    transform = {'toroidal current':'toroidal_current',
                 'it':'toroidal_current',
                 'plasma current':'toroidal_current_p',
                 'ip':'toroidal_current_p',
                 'wall current':'toroidal_current_w',
                 'iw':'toroidal_current_w',
                 'volume':'volume_p',
                 'plasma volume':'volume_p',
                 'volume_d':'volume',
                 'domain volume':'volume',
                 'toroidal flux': 'toroidal_flux_p',
                 'time step': 'dt',
                 'psibound':'psi_lcfs',
                 'psilim':'psi_lcfs',
                 'loop voltage':'loop_voltage',
                 'vl':'loop_voltage',
                 'poloidal magnetic energy':'E_MP',
                 'Wm':'E_MP',
                 'thermal energy':'E_P',
                 'p':'E_P',
                 'electron thermal energy':'E_PE',
                 'pe':'E_PE',
                 'particles':'particle_number',
                 'n':'particle_number',
                 'electrons':'electron_number',
                 'ne':'electron_number',
                 'angular momentum': 'angular_momentum',
                 'vorticity': 'circulation',
                 'parallel viscous heating': 'parallel_viscous_heating',
                 'IZ': 'M_IZ',
                 'ave_p': 'Ave_P',
                 'total current': 'itot',
                 'halo current': 'ih',
                 'kinetic energy': 'ke',
                 'magnetic energy': 'me',
                 'prad': 'radiation',
                 'pline': 'line_rad',
                 'pbrem': 'brem_rad',
                 'pion': 'ion_loss',
                 'preck': 'reck_rad',
                 'precp': 'recp_rad',
                 'POhm': 'pohm',
                 'pelr': 'pellet_rate',
                 'pellet ablation rate': 'pellet_ablrate',
                 'pelablr': 'pellet_ablrate',
                 'pellet var': 'pellet_var',
                 'pelvar': 'pellet_var',
                 'pellet radius': 'r_p',
                 'pelrad': 'r_p',
                 'pellet R position': 'pellet_r',
                 'pelrpos': 'pellet_r',
                 'pellet_x': 'pellet_r',
                 'pellet phi position': 'pellet_phi',
                 'pelphipos': 'pellet_phi',
                 'pellet Z position': 'pellet_z',
                 'pelzpos': 'pellet_z',
                 'poloidal beta': 'betap',
                 'bp': 'betap',
                 }

    # Simple linear combinations
    combos = {'itot':([('toroidal_current_w',1.),('toroidal_current',1.)],
                      {'current':1}, 'Total toroidal current', 'A'),
              'ih':([('toroidal_current',1.),('toroidal_current_p',-1.)],
                    {'current':1}, 'Halo-region toroidal current', 'A'),
              'ke':([('E_KP',1.),('E_KT',1.),('E_K3',1.)], {'energy':1}, 
                    'Kinetic energy', 'J'),
              'me':([('E_MP',1.),('E_MT',1.)], {'energy':1},
                    'Magnetic energy', 'J'),
              'energy':([('E_KP',1.),('E_KT',1.),('E_K3',1.),
                         ('E_MP',1.),('E_MT',1.),('E_P',1.)], {'energy':1},
                        'Total energy', 'J'),
              'flux':([('psimin',2.*np.pi),('psi_lcfs',-2.*np.pi)],
                      {'magnetic_field':1,'length':2}, 'Flux', r'T$\cdot$m$^2$'),
              'radiation':([('radiation',-1.)], None, None, None),
              'line_rad':([('line_rad',-1.)], None, None, None),
              'brem_rad':([('brem_rad',-1.)], None, None, None),
              'ion_loss':([('ion_loss',-1.)], None, None, None),
              'reck_rad':([('reck_rad',-1.)], None, None, None),
              'recp_rad':([('recp_rad',-1.)], None, None, None),
              'rec_rad':([('reck_rad',-1.),('recp_rad',-1.)],
                         {'energy':1,'time':-1}, 'Recombination radiated power', 'W'),
              'pohm':([('E_MPD',-1),('E_MTD',-1)], {'energy':1,'time':-1},
                      'Ohmic heating power', 'W'),
              }

    if trace in transform:
        trace = transform[trace]

    if trace == 'reconnected flux':
        scalar = abs(sim.get_time_trace('reconnected_flux'))
        custom = {'magnetic_field':1,'length':1+itor}
        label = None
        unitlabel = None

    elif trace == 'r_p':
        if (version < 26):
            scalar = sim.get_time_trace('r_p2')
        else:
            scalar = sim.get_time_trace('r_p')
        custom = None
        label = None
        unitlabel = None

    elif trace == 'pellet_r':
        if (version < 26):
            scalar = sim.get_time_trace('pellet_x')
        else:
            scalar = sim.get_time_trace('pellet_r')
        custom = None
        label = None
        unitlabel = None

    elif trace == 'beta':
        if (version < 26):
            scalar = sim.get_time_trace('E_P')
        else:
            scalar = sim.get_time_trace('W_P')
        E_MP = sim.get_time_trace('E_MP')
        E_MT = sim.get_time_trace('E_MT')
        scalar *= (gamma-1.)/(E_MP + E_MT)
        custom = None
        label = r'$\beta$'
        unitlabel = None

    elif trace == 'betap':
        if (version < 26):
            scalar = sim.get_time_trace('E_P')
            it     = sim.get_time_trace('toroidal_current')
            scalar *= 2.*(gamma-1.)/it**2
        else:
            scalar = sim.get_time_trace('W_P')
            W_M    = sim.get_time_trace('W_M')
            scalar *= (gamma-1.)/W_M
        custom = None
        label = r'Poloidal $\beta$'
        unitlabel = None

    elif trace in ['betan','betat']:
        raise RuntimeError("'%s' not yet implemented; need shape information"%trace)

    elif trace == 'electron_number':
        if version <= 20:
            zeff = readParameter('zeff', sim=sim)
            scalar = sim.get_time_trace('particle_number')
            scalar *= zeff
        else:
            scalar = sim.get_time_trace(trace)

        custom = None
        label = None
        unitlabel = None

    elif trace == 'bwb2':
        amupar = constants.amupar
        scalar = sim.get_time_trace('parallel_viscous_heating')
        scalar *= 4./(3.*amupar)
        custom = {'length':3, 'time':-2}
        label = r'$(b\cdot W\cdot b)^2$'
        unitlabel = r'm$^3$/s$^2$'

    elif trace == 'li':
        R0 = constants.R0
        psi_lcfs = sim.get_time_trace('psi_lcfs')
        psimin   = sim.get_time_trace('psimin')
        ip       = sim.get_time_trace('toroidal_current_p')

        scalar = -4.*np.pi*(psi_lcfs - psimin)/(R0*ip)
        custom = None
        label = r'Internal inductance: $l_i$'
        unitlabel = None

    elif trace == 'li3':
        R0 = constants.R0
        W_M = sim.get_time_trace('W_M')
        ip = sim.get_time_trace('toroidal_current_p')

        scalar = 4.*W_M/(R0*ip**2)
        custom = None
        label = r'Internal inductance: $l_i(3)$'
        unitlabel = None

    elif trace in ['kprad_n0','kprad_n']:
        scalar = sim.get_time_trace(trace)*1.0E20
        label = None
        unitlabel = None
        custom = None

    elif trace in ['IZ','M_IZ']:
        scalar = sim.get_time_trace('M_IZ')/sim.get_time_trace('toroidal_current_p')
        label = None
        unitlabel = None
        custom = None



    elif trace == 'Wcond':
        Wcond = integrate_time_trace('Flux_thermal',nts=None,method='cumtrapz',units=units,sim=sim,growth=False,renorm=False,makeplot=False)
        scalar = Wcond
        label = r'$W_{cond}$'
        unitlabel = None
        custom = None

    elif trace == 'frad':
        Wrad = integrate_time_trace('prad',nts=None,method='cumtrapz',units=units,sim=sim,growth=False,renorm=False,makeplot=False)
        #Wth = sim.get_time_trace('E_P')
        Wth,_,_ = get_timetrace('p',sim=sim,units=units,growth=False,renorm=False,returnas='time_trace')
        delta_Wth = Wth.values[0] - Wth
        plt.figure(10)
        plt.plot(delta_Wth.time,delta_Wth.values)
        Wohm = integrate_time_trace('POhm',nts=None,method='cumtrapz',units=units,sim=sim,growth=False,renorm=False,makeplot=False)
        scalar = Wrad/(delta_Wth + Wohm)
        label = r'$W_{rad}/(\Delta W_{th} + W_{ohm})$'
        unitlabel = None
        custom = None

    elif trace == 'energy_balance':
        Wrad = integrate_time_trace('prad',nts=None,method='cumtrapz',units=units,sim=sim,growth=False,renorm=False,makeplot=False)
        #Wth = sim.get_time_trace('E_P')
        Wth,_,_ = get_timetrace('p',sim=sim,units=units,growth=False,renorm=False,returnas='time_trace')
        delta_Wth = Wth.values[0] - Wth
        plt.figure(10)
        plt.plot(delta_Wth.time,delta_Wth.values)
        Wohm = integrate_time_trace('POhm',nts=None,method='cumtrapz',units=units,sim=sim,growth=False,renorm=False,makeplot=False)
        scalar = (Wrad + Wohm)/(delta_Wth)
        label = r'$(W_{rad} + W_{ohm})/\Delta W_{th}$'
        unitlabel = None
        custom = None

    elif trace == 'energy_balance_th':
        Wrad = integrate_time_trace('prad',nts=None,method='cumtrapz',units=units,sim=sim,growth=False,renorm=False,makeplot=False)
        #Wth = sim.get_time_trace('E_P')
        Wth,_,_ = get_timetrace('p',sim=sim,units=units,growth=False,renorm=False,returnas='time_trace')
        delta_Wth = Wth.values[0] - Wth
        Wcond = integrate_time_trace('Flux_thermal',nts=None,method='cumtrapz',units=units,sim=sim,growth=False,renorm=False,makeplot=False)
        
        #plt.figure(11)
        #plt.plot(delta_Wmag.time,delta_Wmag.values)
        #Wohm = integrate_time_trace('POhm',nts=None,method='cumtrapz',units=units,sim=sim,growth=False,renorm=False,makeplot=False)
        scalar = (Wrad + Wcond)/(delta_Wth)
        label = r'$(W_{rad} + W_{cond})/(\Delta W_{th})$'
        unitlabel = None
        custom = None

    elif trace == 'energy_balance_cond':
        Wrad = integrate_time_trace('prad',nts=None,method='cumtrapz',units=units,sim=sim,growth=False,renorm=False,makeplot=False)
        #Wth = sim.get_time_trace('E_P')
        Wth,_,_ = get_timetrace('p',sim=sim,units=units,growth=False,renorm=False,returnas='time_trace')
        delta_Wth = Wth.values[0] - Wth
        Wcond = integrate_time_trace('Flux_thermal',nts=None,method='cumtrapz',units=units,sim=sim,growth=False,renorm=False,makeplot=False)
        
        plt.figure(11)
        plt.plot(Wth.time,Wth.values)
        Wohm = integrate_time_trace('POhm',nts=None,method='cumtrapz',units=units,sim=sim,growth=False,renorm=False,makeplot=False)
        scalar = (Wrad + Wcond - Wohm)/delta_Wth
        label = r'$(W_{rad} + W_{cond} - W_{ohm})/\Delta W_{th}$'
        unitlabel = None
        custom = None


    elif trace == 'energy_sum':
        Wrad = integrate_time_trace('prad',nts=None,method='cumtrapz',units=units,sim=sim,growth=False,renorm=False,makeplot=False)
        #Wth = sim.get_time_trace('E_P')
        Wth,_,_ = get_timetrace('p',sim=sim,units=units,growth=False,renorm=False,returnas='time_trace')
        delta_Wth = Wth - Wth.values[0]
        Wmag,_,_ = get_timetrace('me',sim=sim,units=units,growth=False,renorm=False,returnas='time_trace')
        Wk,_,_ = get_timetrace('ke',sim=sim,units=units,growth=False,renorm=False,returnas='time_trace')
        delta_Wmag = Wmag - Wmag.values[0]
        Wcond = integrate_time_trace('Flux_thermal',nts=None,method='cumtrapz',units=units,sim=sim,growth=False,renorm=False,makeplot=False)
        plt.figure(12)
        plt.plot(Wcond.time,Wcond.values)
        #Wohm = integrate_time_trace('POhm',nts=None,method='cumtrapz',units=units,sim=sim,growth=False,renorm=False,makeplot=False)
        scalar = delta_Wmag    -Wrad + Wk + Wcond + Wmag + Wth
        label = r'$W_{rad} + W_{th} + W_{mag} + W_{k} + W_{cond}$'
        unitlabel = None
        custom = None


    elif trace == 'energy_diff':
        Wrad = integrate_time_trace('prad',nts=None,method='cumtrapz',units=units,sim=sim,growth=False,renorm=False,makeplot=False)
        #Wth = sim.get_time_trace('E_P')
        Wth,_,_ = get_timetrace('p',sim=sim,units=units,growth=False,renorm=False,returnas='time_trace')
        delta_Wth = Wth - Wth.values[0]
        Wmag,_,_ = get_timetrace('me',sim=sim,units=units,growth=False,renorm=False,returnas='time_trace')
        Wk,_,_ = get_timetrace('ke',sim=sim,units=units,growth=False,renorm=False,returnas='time_trace')
        delta_Wmag = Wmag - Wmag.values[0]
        Wcond = integrate_time_trace('Flux_thermal',nts=None,method='cumtrapz',units=units,sim=sim,growth=False,renorm=False,makeplot=False)
        plt.figure(12)
        plt.plot(Wcond.time,Wcond.values)
        Wohm = integrate_time_trace('POhm',nts=None,method='cumtrapz',units=units,sim=sim,growth=False,renorm=False,makeplot=False)
        scalar = Wrad + delta_Wmag + delta_Wth + Wk - Wcond
        label = r'$\Delta W_{th} + \Delta W_{mag} + W_{rad} + W_{k} + W_{cond}$'
        unitlabel = None
        custom = None

    elif trace == 'sideways_force':
        force_x = get_timetrace('Wall_Force_n1_x',sim=sim,units=units,growth=False,renorm=False,returnas='tuple')
        force_y = get_timetrace('Wall_Force_n1_y',sim=sim,units=units,growth=False,renorm=False,returnas='tuple')
        scalar = sim.get_time_trace('E_P')
        scalar.values = np.sqrt(force_x[1]*force_x[1] + force_y[1]*force_y[1])
        label = r'sideways force'
        unitlabel = None
        custom = None

    elif trace in combos:
        # trace is linear combination of native scalars
        combo, custom, label, unitlabel = combos[trace]
        for i, (name, fact) in enumerate(combo):
            y = sim.get_time_trace(name)
            if i==0:
                scalar = fact*y
            else:
                scalar += fact*y

    else:
        # trace is a native scalar
        scalar = sim.get_time_trace(trace)
        custom = None
        label = None
        #unitlabel = None

    if ('pellet_' in trace) or (trace in ['cauchy_fraction','cloud_pel','r_p']):
        # if ipellet is given, get just that pellet's data
        if (ipellet != 'all') and (scalar.values.ndim==2):
            scalar.values = scalar.values[:,ipellet]

    label, unitlabel = fpyl.get_tracelabel(units, trace, label=label, unitlabel=unitlabel,fac=fac)
    if units=='mks' and trace not in ['energy_balance','energy_balance_cond','energy_balance_th','energy_diff','energy_sum','frad']:
        scalar = fpyl.get_conv_trace('mks',trace,scalar,sim=sim,itor=itor,custom=custom)
    
    # now separate time and values arrays
    time = scalar.time
    values = scalar.values
    
    if growth:
        values = 1.0/values[1:] * np.diff(values)/np.diff(time) #Used until 2021-05-12
        #values = 1.0/values[1:] * fpyl.deriv(values,time)
        time = time[:-1]
    
    if diff:
        values = np.diff(values)/np.diff(time)
        time = time[:-1]
    
    if renorm:
        renormlist = []
        for i in range(len(values)-1):
            if(abs(values[i+1]/values[i]) < 1E-9):
                renormlist.append(str(time[i]))
                #print(values[i],values[i-1]+values[i+1])
                # Only average value if growth rate is calculated
                if growth:
                    values[i] = (values[i-1] + values[i+1])/2.0
        # When growth rate is calculated, check for normalization at last time
        # step and drop this point, since it carries no information.
        if growth:
            if(abs(values[-2]/values[-1]) < 1E-9):
                renormlist.append(str(time[-1]))
                values = values[:-1]
                time = time[:-1]
        renormstr = ", ".join(renormlist)
        if not quiet:
            if len(renormstr) > 0:
                print('Renormalization found at '+renormstr)
    
    if returnas=='tuple':
        return time, values*fac, label, unitlabel
    elif returnas=='time_trace':
        return fpy.sim_data.time_trace(values*fac,time=time), label, unitlabel


