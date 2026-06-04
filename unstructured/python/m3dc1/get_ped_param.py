#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Andreas Kleiner
"""
import fpy
import numpy as np
from scipy import signal
import m3dc1.fpylib as fpyl
from m3dc1.flux_average import flux_average
from m3dc1.omegastari import omegastari


def get_ped_param(sim,filename='C1.h5',time=None,points=400,pion=False,fcoords='pest',psin_ped_top=0.86,psin_var_j=0.85,use_max_j=False,device=None):
    if not isinstance(sim,fpy.sim_data):
        sim = fpy.sim_data(filename=filename)
    #psinedgelim = 0.86 #Min value of psin that is considered edge region
    print('get_ped_param time:'+str(time))
    # Determine pedestal alpha
    psi_a,alpha = flux_average('alpha', coord='scalar', sim=sim, time=time, fcoords=fcoords, points=points, units='m3dc1')
    psinedge = fpyl.find_nearest(psi_a,psin_ped_top)
    psinedge_ind = fpyl.get_ind_at_val(psi_a,psinedge)
    alpha_max = np.amax(alpha[psinedge_ind:])
    #Check if it's really a local maximum inside the pedestal:
    alpha_max_ind = fpyl.get_ind_at_val(alpha,alpha_max)
    #print(alpha_max_ind,psi_a[alpha_max_ind])
    if (alpha_max_ind < len(alpha)-1 and (alpha_max >= alpha[alpha_max_ind-1]) and (alpha_max >= alpha[alpha_max_ind+1])) or (alpha_max_ind==len(alpha)-1): #Check if the maximum is a relative maximum or occurs at the lcfs
        #print('PEDESTAL ALPHA IS ABSOLUTE MAXIMUM')
        print('Pedestal alpha = '+str(alpha_max))
    else:
        alpha_short = alpha[psinedge_ind:]
        psin_a_short = psi_a[psinedge_ind:]
        alpha_rel_max = signal.argrelmax(alpha_short)
        print(len(alpha_rel_max[0]))
        print(alpha_rel_max[0].size)
        if len(alpha_rel_max[0])==0:
            fpyl.printwarn('WARNING: No local maximum of alpha found! Check alpha profile and make sure value for psin_ped_top is correct!')
            return None,None,None,None
        maxima = np.take(alpha_short,alpha_rel_max)
        maxima_pos = np.take(psin_a_short,alpha_rel_max)
        alpha_max = np.amax(maxima)
        #print(maxima,maxima_pos)
        #print('PEDESTAL ALPHA IS RELATIVE MAXIMUM')
        print('Pedestal alpha = '+str(alpha_max))
    
    # Determine pedestal average toroidal current density
    psi_j,j = flux_average('j', coord='phi', sim=sim, time=time, fcoords=fcoords, points=points, units='m3dc1')
    psinedge = fpyl.find_nearest(psi_j,psin_ped_top)
    psinedge_ind = fpyl.get_ind_at_val(psi_j,psinedge)
    j_max = np.amax(j[psinedge_ind:])
    #Check if it's really a local maximum inside the pedestal:
    j_max_ind = fpyl.get_ind_at_val(j,j_max)
    if j_max_ind < len(j)-1 and (j_max >= j[j_max_ind-1]) and (j_max >= j[j_max_ind+1]):
        print('Pedestal j_max = '+str(j_max))
    else:
        j_max = -1
        fpyl.printwarn('WARNING: j_max has no maximum inside the pedestal!')
        #return None,None,None,None
    
    #Calculate pedestal parallel current density as in ELITE:
    jav = flux_average('jav', coord='scalar', sim=sim, time=time, fcoords=fcoords, points=points, units='mks')[1]
    psi_j,jelite = flux_average('jelite', coord='scalar', sim=sim, time=time, fcoords=fcoords, points=points, units='mks',device=device)
    psinedge = fpyl.find_nearest(psi_j,psin_ped_top)
    psinedge_ind = fpyl.get_ind_at_val(psi_j,psinedge)
    
    if use_max_j:
        #If maximum jelite inside the pedestal region is to be used:
        jelite_max = np.amax(jelite[psinedge_ind:])
    else:
        #Check if jelite peaks inside the pedestal. If it does, take peak value as jelite. If it does
        #not peak, then use jelite(psin=psin_var_j).
        #Calculate absolute maximum in pedestal region:
        jelite_max = np.amax(jelite[psinedge_ind:])
        jelite_max_ind = fpyl.get_ind_at_val(jelite,jelite_max)
        #print(jelite_max_ind, psi_j[jelite_max_ind], jelite_max,jelite[jelite_max_ind-1],jelite[jelite_max_ind+1])
        #If local maximum
        if jelite_max_ind < len(jelite)-1 and (jelite_max >= jelite[jelite_max_ind-1]) and (jelite_max >= jelite[jelite_max_ind+1]):
            jelite_rel_max = jelite_max
            fpyl.printnote('Used MAXIMUM for jelite')
        else:
            jelite_rel_max = 0.0
        #Calculate average value in pedestal region
        jelite_avg = np.average(jelite[psinedge_ind:])
        #print(jelite_rel_max,jelite_avg)
        j_threshold=1.2
        #If jelite_rel_max > j_threshold*jelite_avg then we consider a peak inside the pedestal.
        #Otherwise, use jelite(psin=psin_var_j):
        if not jelite_rel_max > j_threshold*jelite_avg:
            psinj = fpyl.find_nearest(psi_j,psin_var_j)
            psinj_ind = fpyl.get_ind_at_val(psi_j,psinj)
            jelite_max = jelite[psinj_ind]
            fpyl.printwarn('Used psin_var_j for jelite')
    jelite_sep = jelite[-1]
    jelite_N = (jelite_max+jelite_sep)/(2.0*jav[-1])
    print('Current density jelite = '+str(jelite_N))
    
    
    # Determine diamagnetic frequency
    psi_o,omegsti = omegastari(sim=sim,time=time,units='m3dc1',n=1,points=points,pion=pion,fcoords=fcoords,makeplot=False)
    psinedge = fpyl.find_nearest(psi_o,psin_ped_top)
    psinedge_ind = fpyl.get_ind_at_val(psi_o,psinedge)
    omegsti_max = np.amax(omegsti[psinedge_ind:]) #ToDo: Check if it's really a local maximum
    print('Ion diamagnetic frequency = '+str(omegsti_max))
    ped_param = [alpha_max,j_max,jelite_N,omegsti_max]
    return ped_param
