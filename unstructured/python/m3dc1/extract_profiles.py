#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on July 16 2021

@author: Andreas Kleiner
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d
import m3dc1.fpylib as fpyl


def extract_profiles(filename):
    """
    Extracts profiles from p-file-like text file and converts pressure to Pascal.

    Arguments:

    **filename**
    File name which will be read, i.e. "../../C1.h5"
    Can also be a list of two filepaths when used for diff
    """
    os.system('extract_profiles.sh '+filename)
    
    if os.path.isfile('profile_p'):
        convert_p()
        print('profile_p converted to Pascal')
    
    return


def convert_p():
    
    f = open('profile_p', 'r')
    data = f.readlines()
    f.close()
    nlines = len(data)
    
    profile = np.empty((nlines,3), dtype=float)
    
    for i,line in enumerate(data):
        profile[i,:] = list(map(float, line.split()))
    profile[:,1] = profile[:,1]*1000
    profile[:,2] = profile[:,2]*1000
    print(profile)
    
    with open('profile_p', 'w') as pf:
        for ln in profile:
            pf.write(" {0:8.6f}   {1:>9.3f}   {2:>9.3f}\n".format(*ln))
    
    return



def smooth_profile(psin,values,psirange):
    psirange = np.asarray(psirange)
    
    profile = np.column_stack([psin,values])
    profile_old = np.copy(profile)
    
    if len(psirange.shape)==1:
        psirange = np.asarray([psirange])
    
    psi_replace = []
    psi_replace_ind = []
    for i,data in enumerate(profile):
        for interval in psirange:
            #print(interval)
            if data[0]>=interval[0] and data[0]<=interval[1]:
                psi_replace.append(data[0])
                psi_replace_ind.append(i)
                #profile = np.delete(profile,i,0)
    
    #print(psi_replace_ind)
    profile_removed = np.delete(profile,psi_replace_ind,0)
    #print(profile)
    
    f = interp1d(profile_removed[:,0], profile_removed[:,1], kind='cubic')
    new_values = f(psi_replace)
    #print(psi_replace)
    #print(new_values)
    
    for i in range(len(psi_replace_ind)):
        profile[psi_replace_ind[i],1] = new_values[i]
    
    
    
    fig = plt.figure(constrained_layout=True,figsize=(10,5))
    spec2 = gridspec.GridSpec(ncols=2, nrows=1, figure=fig)
    f2_ax1 = fig.add_subplot(spec2[0, 0])
    f2_ax2 = fig.add_subplot(spec2[0, 1])
    
    
    f2_ax1.plot(profile_removed[:,0], profile_removed[:,1],lw=0,marker='.',label='original (points removed)')
    f2_ax1.plot(profile_old[:,0],profile_old[:,1],ls='--',lw=2,label='original')
    f2_ax1.plot(profile[:,0],profile[:,1],lw=2,label='smoothed profile')
    f2_ax1.set_xlabel(r'$\psi_N$')
    f2_ax1.set_ylabel('Rotation profile')
    f2_ax1.grid()
    #plt.legend(loc=0)
    
    f2_ax2.plot(profile[:,0],fpyl.deriv(profile[:,1],profile[:,0]))
    f2_ax2.set_xlabel(r'$\psi_N$')
    f2_ax2.set_ylabel('Rotation profile derivative')
    f2_ax2.grid()
    
    
    
    #plt.plot(profile_old[:,0],test_values)
    
    return profile
