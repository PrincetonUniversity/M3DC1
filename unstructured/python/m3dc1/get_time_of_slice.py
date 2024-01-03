#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 14:24:12 2019

@author: Andreas Kleiner
"""

import numpy as np
import fpy
from m3dc1.unit_conv import unit_conv
from m3dc1.read_h5 import readParameter


def get_time_of_slice(time,sim=None,filename='C1.h5',units='mks',millisec=False,quiet=False):
    """
    Return time of time slice in s or M3D-C1 units.
    
    Arguments:

    **sim**
    fpy simulation object.

    **filename**
    Name or path to C1.h5 file to read

    **units**
    The units in which the time trace will be returned

    **millisec**
    True/False. If True and units='mks' plot will be in terms of milliseconds, instead of seconds.

    **quiet**
    If True, do not print renormalization times to screen.
    """
    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename=filename,time=time)
    else:
        time = sim.timeslice
    
    ts = sim.timeslice
    
    h5name = 'time_' + str(ts).zfill(3) + '.h5' if time > -1 else 'equilibrium.h5'
    #print(h5name)
    #Append path to directory to h5name, in case a C1.h5 file in a different directory was specified
    if len(filename)>5:
        h5name = filename[:-5] + h5name
    time = readParameter('time',fname=h5name) #Don't pass sim, because it expects C1.h5 file instead of time slice
    if units == 'mks':
        time = unit_conv(time, sim=sim,time=1)
        
        unit_lbl = 's'
        if millisec:
            time = time*1000
            unit_lbl = 'm'+unit_lbl
    else:
        unit_lbl = 'tau_A'
        
    if not quiet:
        print('Time ('+unit_lbl+') of time slice ' + str(ts) +': ' + str(time))
    return time
