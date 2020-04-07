#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 14:24:12 2019

@author: Andreas Kleiner
"""
import numpy as np
import h5py


def openH5File(fname,listc=False):
    f = h5py.File(fname,'r')
    if listc==True:
        print('=========================================================\nGroups and datasets of file '+str(fname)+'\n---------------------------------------------------------')
        f.visititems(print)
        
        print('\n\n=========================================================\nAttributes in the root group of file '+str(fname)+'\n---------------------------------------------------------')
        for item in f.attrs.keys():
            print(item + ":", f.attrs[item])
    return f



def readParameter(pname,fname='C1.h5',listc=False):
    """
    Read single parameter from M3DC1 C1.h5 output file.
    
    Arguments:

    **pname**
    Name of parameter to read.

    **fname**
    Name of file to read.

    **listc**
    Print content of h5 file to screen.
    """
    if fname=='time_-01.h5':
        fname='equilibrium.h5'
    f = openH5File(fname,listc)
    param = f.attrs[pname]
    return param



def readC1File(scalar=None,signal=None,fname='C1.h5',listc=False):
    """
    Read M3DC1 C1.h5 output file and either return data or show content.
    
    Arguments:

    **scalar**
    Name of scalar to read. Returns two arrays: time and scalar

    **signal**
    Name of signal (diagnostics probe) to read. Returns two arrays: time and signal

    **fname**
    Name of file to read.

    **listc**
    Print content of h5 file to screen.
    
    """
    if ((scalar != None and signal != None) or (scalar == None and signal == None)) and (listc==False):
        raise Exception('Please provide either a scalar or a signal name!')
    f = openH5File(fname,listc)
    #if listc==True:
    #    print('\n\n=========================================================\nAttributes in the root group of file '+str(fname)+'\n---------------------------------------------------------')
    #    for item in f.attrs.keys():
    #        print(item + ":", f.attrs[item])
    
    if scalar != None or signal!=None:
        time = np.asarray(f['scalars/time'])
        if scalar != None and signal==None:
            if scalar == 'ke':
                tracedata1 = list(f['scalars/E_KP'])
                tracedata2 = list(f['scalars/E_KT'])
                tracedata3 = list(f['scalars/E_K3'])
                trace = np.asarray(tracedata1) + np.asarray(tracedata2) + np.asarray(tracedata3)
            elif scalar == 'me':
                tracedata1 = list(f['scalars/E_MP'])
                tracedata2 = list(f['scalars/E_MT'])
                trace = np.asarray(tracedata1) + np.asarray(tracedata2)
            else:
                trace = np.asarray(f['scalars/'+scalar])
            return time, trace
        elif signal != None and scalar==None:
            probes = np.asarray(f[signal+'/value'])
            nprobes = probes.shape[1]
            #print(probes.shape)
            trace = np.zeros((nprobes,len(time)))
            for i in range(nprobes):
                for j in range(len(time)):
                    trace[i,j] = probes[j,i]
            #trace = trace.flatten()
            return time, trace



def readTimeFile(field=None,fname='time_000.h5',listc=False):
    """
    Read M3DC1 time_xxx.h5 output files containing field information.
    
    Arguments:

    **field**
    Name of field to read.

    **fname**
    Name of file to read.

    **listc**
    Print content of h5 file to screen.
    """
    f = openH5File(fname,listc)
    
    if listc==True:
        for item in f.attrs.keys():
            print(item + ":", f.attrs[item])
    
    if type(field) == str:
        field = list(f['fields/'+field])
    #v = list(f['fields/V'])
    return field
