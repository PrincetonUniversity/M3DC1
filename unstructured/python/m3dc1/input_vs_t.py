#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Feb  7 15:17:00 2023

@author: Andreas Kleiner
"""
import os
import glob
import traceback
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator, MultipleLocator, FormatStrFormatter
from m3dc1.get_time_of_slice import get_time_of_slice

def input_vs_t(param,wd=None,units='m3dc1',millisec=False,makeplot=True,fignum=None):
    """
    Reads M3D-C1 input parameter from all Slurm log files within simulation directory
    and plots it vs. time. This is useful to track the value of input parameters
    that are not written to the HDF5 output in restarted simulations.
    
    Arguments:

    **param**
    Parameter name as used in C1input, e.g. denm

    **wd**
    Path to simulation directory. If None, current working directory is used.

    **units**
    Units in which the result will be calculated

    **millisec**
    True/False. If True and units='mks' plot will be in terms of milliseconds, instead of seconds.

    **makeplot**
    If True, show plot, otherwise results will be printed to screen and returned.

    **fignum**
    Figure number.
    """
    if wd is None:
        wd = os.getcwd()
    slurmfiles = glob.glob(wd+"/slurm*.out")
    slurmfiles.sort(key=os.path.getmtime,reverse=False)
    print(slurmfiles)
    
    times = []
    values = []
    
    for slurmfile in slurmfiles:
        with open(slurmfile, 'r') as sf:
            found_restart = False
            found_param = False
            try:
                for line in sf:
                    if 'irestart_slice' in line:
                        restart_line = line.replace("\n", "")
                        found_restart = True
                    if param+' ' in line:
                        line_split = line.split()
                        if param == line_split[0]:
                            param_line = line.replace("\n", "")
                            value = line_split[-1]
                        
                            found_param = True
                    if found_restart and found_param:
                        break
            except Exception as err:
                print(slurmfile)
                traceback.print_exc()
                print(err)
                
        
        values.append(float(value))
        ts = int(restart_line.split()[-1])
        time = get_time_of_slice(ts,units=units,millisec=millisec)
        times.append(time)
        print(restart_line)
        print(param_line)
        print('------------------------------------')
    
    print(times)
    print(values)
    if makeplot:
        plt.figure(num=fignum)
        plt.plot(times,values)
        plt.xlabel('time')
        plt.ylabel(param)
        
    return

