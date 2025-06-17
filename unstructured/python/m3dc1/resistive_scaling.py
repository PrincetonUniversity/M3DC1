#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 14:24:12 2019

@author: Andreas Kleiner
"""
import os
import glob
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import m3dc1.fpylib as fpyl
from m3dc1.gamma_file import Gamma_file


def resistive_scaling(dirs,eta_fac_list,Sval,ntor,gamma_list=None,col='C1',mk='.',ls='-',logscale=None,fignum=1,figsize=None,max_fit=-1,plot_analytical=True,connect_fit_data=True,label='',in_plot_txt=None,fit_guess=[0.12,-3/5]):
    
    gammas = []
    
    if dirs is not None:
        pwd = os.getcwd()
        
        for i,d in enumerate(dirs):
            os.chdir(d)
            print(d)
            
            files = glob.glob("growth_rates*.dat")
            
            if len(files)>0:
                if len(files)>1:
                    files.sort(key=os.path.getmtime,reverse=True)
                    print('More than 1 file found. Opening newest file: '+files[0])
                f = files[0]
                results = Gamma_file(f)
            else:
                fpyl.printerr('ERROR: no growth rates found for ' + str(d))
                os.chdir(pwd)
                return
            
            # Identify simulations where the growth rate was not calculated reliably. These are highlighted in the plot.
            #n_okay = []
            #gamma_okay = []
            #n_bad = []
            #gamma_bad = []
            for n in results.n_list:
                #n_ind = n-1
                n_ind = fpyl.get_ind_at_val(results.n_list,n,unique=True)
                if n==ntor:
                    if (results.flat_list[n_ind] != 1 and results.gamma_manual_list[n_ind] != 1 and results.not_noisy_list[n_ind] != 1):
                        fpyl.printwarn('WARNING: Check growth rate!')
                        #n_bad.append(results.n_list[n_ind])
                        #gamma_bad.append(results.gamma_list[n_ind])
                    #else:
                    gammas.append(results.gamma_list[n_ind])
            
            os.chdir(pwd)
    elif gamma_list is not None:
        for i,g in enumerate(gamma_list):
            gammas.append(g[ntor-1])
    
    Slist = Sval/np.asarray(eta_fac_list)
    #print(Slist)
    
    a = fit_guess
    yfit,yerr = curve_fit(poly, Slist[:max_fit], gammas[:max_fit], p0=a)
    #print(yfit)
    
    xlist = np.linspace(Slist[0],Slist[-1],50)
    ylist = poly(xlist,yfit[0],yfit[1])
    if plot_analytical:
        tearing_scaling = (yfit[0]*Slist[0]**(yfit[1])/Slist[0]**(-3/5))*xlist**(-3/5)
        infernal_scaling = (yfit[0]*Slist[0]**(yfit[1])/Slist[0]**(-3/13))*xlist**(-3/13)
    
    if mk in ['s','D']:
        ms=6
    elif mk in ['x','v','P']:
        ms=8
    else:
        ms=14
    plt.figure(num=fignum,figsize=figsize)
    
    if connect_fit_data:
        fitlw = 2
        fitls = '--'
    else:
        fitlw = 0
        fitls = '-'
    
    if label=='':
        label = r'n= '+str(ntor)+' : $\\gamma \\propto S^{'+"{0:.3f}".format(yfit[1])+'}$'
    if logscale==True:
        if plot_analytical:
            plt.loglog(xlist,tearing_scaling,lw=3,c='C0',ls=':',label=r'$\gamma \propto S^{-3/5}$ tearing')
            plt.loglog(xlist,infernal_scaling,lw=2,c='C9',ls='-.',label=r'$\gamma \propto S^{-3/13}$ resistive infernal')
        plt.loglog(Slist,gammas,lw=fitlw,c=col,ls=fitls,marker=mk,ms=ms,label=label)#-3/4
        plt.loglog(xlist,ylist,lw=2,c=col,ls=ls)

    else:
        if plot_analytical:
            plt.plot(xlist,tearing_scaling,lw=3,c='C0',ls=':',label=r'$\gamma \propto S^{-3/5}$ tearing')
            plt.plot(xlist,infernal_scaling,lw=2,c='C9',ls='-.',label=r'$\gamma \propto S^{-3/13}$ resistive infernal')
        plt.plot(Slist,gammas,lw=fitlw,c=col,ls=fitls,marker=mk,ms=ms,label=label)#-3/4
        plt.plot(xlist,ylist,lw=2,c=col,ls=ls)
        
    plt.grid(True)
    ax = plt.gca()
    legend = ax.get_legend()
    if legend is not None:
        leglines = legend.get_lines()
        leglabels = [text.get_text() for text in legend.get_texts()]
    else:
        leglines = []
        leglabels = []
    leglines.append(Line2D([0,1],[0,1],linestyle=ls, linewidth=2, color=col,marker=mk,ms=ms))
    leglabels.append(label)
    #print(leglines)
    #print(leglabels)
    plt.legend(leglines,leglabels,loc=0,fontsize=11,ncol=2)
    ax = plt.gca()
    ax.set_xlim([None,np.amax(Slist)])
    ax.tick_params(axis='both', which='major', labelsize=14)
    plt.xlabel(r'Lundquist number $S$',fontsize=18)
    plt.ylabel(r'$\gamma/\omega_A$',fontsize=18)
    
    if in_plot_txt is not None:
        plt.text(0.06, 0.9,in_plot_txt, ha='center', va='center', transform=ax.transAxes,fontsize=20)
    
    plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    return



# Function to fit when minval is not specified.
def poly(x, a,b):
    f = a*x**b
    return f





def plot_gamma_etafac(dirs,ntor,logscale=True,fignum=1,figsize=None,calculate_scaling=False,max_fit=-1):
    
    gammas = []
    etalist = []
    
    pwd = os.getcwd()
    print(dirs)
    for i,d in enumerate(dirs):
        print(dirs)
        os.chdir(d)
        
        
        files = glob.glob("growth_rates*.dat")
        
        if len(files)>0:
            if len(files)>1:
                files.sort(key=os.path.getmtime,reverse=True)
                print('More than 1 file found. Opening newest file: '+files[0])
            f = files[0]
            results = Gamma_file(f)
        else:
            fpyl.printerr('ERROR: No results found.')
            return
            #eval_input = fpyl.prompt('No results found. Do you want to determine the growth rates? (y/n) : ',['y','n'])
            #if eval_input == 'y':
            #    results = scan_n(nmin,nmax,nstep,units=units,slurm=True,plottrace=plottrace)
            #else:
            #    return
        
        etalist.append(results.eta)
        
        # Identify simulations where the growth rate was not calculated reliably. These are highlighted in the plot.
        n_okay = []
        gamma_okay = []
        n_bad = []
        gamma_bad = []
        for n in results.n_list:
            #n_ind = n-1
            n_ind = fpyl.get_ind_at_val(results.n_list,n,unique=True)
            if n==ntor:
                if (results.flat_list[n_ind] != 1 and results.gamma_manual_list[n_ind] != 1 and results.not_noisy_list[n_ind] != 1):
                    fpyl.printwarn('WARNING: Check growth rate!')
                    #n_bad.append(results.n_list[n_ind])
                    #gamma_bad.append(results.gamma_list[n_ind])
                #else:
                gammas.append(results.gamma_list[n_ind])
        
        os.chdir(pwd)
    
    if calculate_scaling:
        a = [0.12,-3/5]
        yfit,yerr = curve_fit(poly, etalist[:max_fit], gammas[:max_fit], p0=a)
        print(yfit)
    
        xlist = np.linspace(etalist[0],etalist[-1],50)
        ylist = poly(xlist,yfit[0],yfit[1])
        #tearing_scaling = (yfit[0]*etalist[0]**(yfit[1])/etalist[0]**(-3/5))*xlist**(-3/5)
        #infernal_scaling = (yfit[0]*etalist[0]**(yfit[1])/etalist[0]**(-3/13))*xlist**(-3/13)
    
    plt.figure(num=fignum,figsize=figsize)
    plt.plot(etalist,gammas,lw=2,marker='.',ms=14,label=r'$n='+str(ntor)+'$')
    if calculate_scaling:
        plt.plot(xlist,ylist)
    ax = plt.gca()
    if logscale=='x':
        ax.set_xscale('log')
    elif logscale=='y':
        ax.set_yscale('log')
    elif logscale=='both':
        ax.set_xscale('log')
        ax.set_yscale('log')
    plt.grid(True)
    plt.legend(loc=0,fontsize=12)
    
    ax.set_xlim([None,np.amax(etalist)])
    ax.tick_params(axis='both', which='major', labelsize=14)
    plt.xlabel(r'resistivity scale factor',fontsize=18)
    plt.ylabel(r'$\gamma/\omega_A$',fontsize=18)
    plt.tight_layout() #adjusts white spaces around the figure to tightly fit everything in the window
    return

