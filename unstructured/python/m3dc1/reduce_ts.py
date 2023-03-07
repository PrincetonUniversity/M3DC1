#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  26 18:31:00 2021

@author: Andreas Kleiner
"""
#import os
#import glob
import numpy as np
import h5py
import os
#import math
import m3dc1.fpylib as fpyl
from m3dc1.read_h5 import readC1File
from m3dc1.read_h5 import readParameter


#def list_data(name, obj,test):
#    return name
#    #print(test)
#    #for key, val in obj.attrs.items():
#    #    print("    %s: %s" % (key, val))


def reduce_ts(fields, time=0, write_mesh=True):
    """
    Reads existing time slice and writes specified fields into a new HDF5 file.

    Arguments:

    **fields**
    List of fields to keep in the new time slice file.

    **time**
    The integer number identifying the time slice for which the fields will be written.

    **write_mesh**
    If True, retain mesh information in new file. If False, only the fields will be written.
    """
    new_dir = 'time_slices_reduced'
    fname = 'time_'+str(time).zfill(3) if time>=0 else 'equilibrium'
    fname_src = fname + '.h5'
    fname_out = new_dir + '/' + fname + '.h5'
    
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)
    
    if os.path.exists(fname_out):
        fpyl.printwarn('WARNING: File "' + fname_out + '" already exists.')
        proceed = fpyl.prompt('Do you want to overwrite it? (y/n) : ',['y','n'])
    else:
        proceed = 'y'
    
    if proceed=='y':
        os.system('cp ' + fname_src + ' ' + fname_out)
        
        #f = h5py.File(fname_src, "r")
        
        #all_fields = list(f['fields'].keys())
        #print(fields)
        #print(all_fields)
        ##f = h5py.File(fname_out, "r")
        
        #with h5py.File(fname_out,  "a") as fout:
            #for fi in all_fields:
                #if fi not in fields:
                    #del fout['fields/'+fi]
        
        
        fs = h5py.File(fname_src, 'r')
        fd = h5py.File(fname_out, 'w')
        for a in fs.attrs:
            fd.attrs[a] = fs.attrs[a]
        
        try:
            if fields[0] is not None:
                fd.create_group('fields')
                for fi in fields:
                    fs.copy(fs['fields/'+fi], fd['fields'],fi)
            
            if write_mesh:
                fs.copy('mesh', fd)
            
            fd.visititems(print)
        except Exception as e:
            print(e)
        
        fs.close()
        fd.close()
    return
