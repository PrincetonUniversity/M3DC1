#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 14:24:12 2019

@author: Andreas Kleiner
"""
import numpy as np
import h5py


def list_contents(fname='C1.h5',h5file=None):
    if h5file is None:
        h5file = openH5File(fname)
    print('=========================================================\nGroups and datasets of file '+str(fname)+'\n---------------------------------------------------------')
    h5file.visititems(print)
    print('\n\n=========================================================\nAttributes in the root group of file '+str(fname)+'\n---------------------------------------------------------')
    for item in h5file.attrs.keys():
        print(item + ":", h5file.attrs[item])
    return

def openH5File(fname):
    if fname=='time_-01.h5':
        fname='equilibrium.h5'
    return h5py.File(fname,'r')

def readParameter(pname,fname='C1.h5',h5file=None,listc=False):
    """
    Read single parameter from M3DC1 C1.h5 output file.
    
    Arguments:

    **pname**
    Name of parameter to read.

    **fname**
    Name of file to read.

    **h5file**
    Provide output of h5py.File if already opened.

    **listc**
    Print content of h5 file to screen.
    """
    if h5file is None:
        h5file = openH5File(fname)
    if listc:
        list_contents(h5file=h5file)
    param = h5file.attrs[pname]
    return param


def readC1File(scalar=None,signal=None,fname='C1.h5',h5file=None,listc=False):
    """
    Read M3DC1 C1.h5 output file and either return data or show content.
    
    Arguments:

    **scalar**
    Name of scalar to read. Returns two arrays: time and scalar

    **signal**
    Name of signal (diagnostics probe) to read. Returns two arrays: time and signal

    **fname**
    Name of file to read.

    **h5file**
    Provide output of h5py.File if already opened.

    **listc**
    Print content of h5 file to screen.
    
    """
    if (scalar is None) == (signal is None):
        raise RuntimeError('Please provide either a scalar or a signal name!')
    if h5file is None:
        h5file = openH5File(fname)
    if listc:
        list_contents(h5file=h5file)

    version = readParameter('version', h5file=h5file)

    time = np.asarray(h5file['scalars/time'])
    Nt = len(time)
    values = None
    trace = None
    if scalar is not None:
        if scalar in h5file['scalars']:
            trace = np.asarray(h5file['scalars/'+scalar])
        elif ('pellet' in h5file) and (scalar in h5file['pellet']):
            values = np.asarray(h5file['pellet/'+scalar])
    elif signal in h5file:
        values = np.asarray(h5file[signal+'/value'])

    if values is not None:
        Np = values.shape[1]
        trace = np.zeros((Np, Nt))
        for i in range(Np):
            for j in range(Nt):
                trace[i,j] = values[j,i]

    if trace is None:
        if scalar is not None:
            raise RuntimeError("Scalar '%s' not found in %s"%(scalar,h5file.filename))
        else:
            raise RuntimeError("Signal '%s' not found in "%(signal,h5file.filename))
    
    return time, trace


def readTimeFile(field=None,fname='time_000.h5',h5file=None,listc=False):
    """
    Read M3DC1 time_xxx.h5 output files containing field information.
    
    Arguments:

    **field**
    Name of field to read.

    **fname**
    Name of file to read.

    **h5file**
    Provide output of h5py.File if already opened.

    **listc**
    Print content of h5 file to screen.
    """
    if h5file is None:
        h5file = openH5File(fname)
    if listc:
        list_contents(h5file=h5file)
    
    if type(field) == str:
        field = list(h5file['fields/'+field])
    return field
