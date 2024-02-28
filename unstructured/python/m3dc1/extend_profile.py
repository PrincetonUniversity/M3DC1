#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 16 2020

@author: Andreas Kleiner
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
import m3dc1.fpylib as fpyl



def extend_profile(filename,psimax=1.05,psimatch=-1,fitrange=None,minval=None,match=True,smooth=0,weighted=True,suffix='.extended',back_to_maindir=False):
    """
    Extends profile beyond the last closed flux surface (psi_n=1) based on a tanh fit. If desired,
    the fit parameters are adjusted such that the profile and its first derivative are continuous
    across the lcfs.

    Arguments:

    **filename**
    File containing the profile information, e.g. 'profile_ne'

    **psimax**
    Value of normalized psi to which the profile will be extended

    **fitrange**
    Range of normalized psi that is considered for tanh fit

    **minval**
    Asymptotic value that is approached at psi_N -> infinity

    **match**
    If True, the profile and its first derivative are matched at the lcfs.
    When using this option set smooth=0.

    **smooth**
    If smooth>0 a boxcar average of width=smooth is performed. Use only
    if match=False.

    **weighted**
    If True, weight fit by 1/y

    **suffix**
    File extension for output file
    """
    x,y = fpyl.ReadTwoColFile(filename)
    
    if isinstance(fitrange, float) or isinstance(fitrange, int):
            fitrange = [fitrange,max(x)]
    elif not isinstance(fitrange, (tuple, list,float,int)):
        fitrange = [0.95,max(x)]
    
    #Identify index of psi value where extended profile begins.
    if psimatch>0:
        x_match = fpyl.get_ind_near_val(x, psimatch,unique=True)
    else:
        x_match = len(x)-1
    
    print('Fitting points in range '+str(fitrange))
        
    if psimax <= max(x):
        fpyl.printerr('ERROR: profile already extends to '+str(psimax))
        return
    
    i = (np.argwhere((x>=fitrange[0]) & (x<=fitrange[1]))).flatten()
    c = len(i)
    if c == 0 :
        fpyl.printerr('ERROR: no data points in range')
        return
    #print(c)
    
    xf = x[i]
    yf = y[i]
    
    # Fit profile
    w = yf
    try:
        if minval is None:
            a = [0.98, 1./0.01, max(yf), 0., min(yf)]
            if weighted:
                yfit,yerr = curve_fit(tanhfit2, xf, yf, p0=a, sigma=w, absolute_sigma=True)
            else:
                yfit,yerr = curve_fit(tanhfit2, xf, yf, p0=a)
        else:
            a = [0.98, 1./0.01, max(yf), 0.]
            if weighted:
                yfit,yerr = curve_fit(tanhfit, xf, yf-minval, p0=a, sigma=w, absolute_sigma=True)
            else:
                yfit,yerr = curve_fit(tanhfit, xf, yf-minval, p0=a)
    except Exception as e:
        if back_to_maindir:
            dirs = os.getcwd().split('/')
            os.chdir('../../')
            os.system('rm -rf '+'/'.join(dirs[-2:]))
        raise Exception(e)
        return
        
    #print(yfit)
    
    print('---------------------------------------------------------')
    print('Fit center: '+str(yfit[0]))
    print('Fit width: '+str(1.0/yfit[1]))
    print('Fit height: '+str(yfit[2]))
    if minval is None:
        print('Fit floor: '+str(yfit[4]))
        if yfit[4]<=0:
            fpyl.printerr('ERROR: fit approaches negative values for large psi_n')
            fpyl.printerr('Please provide a value for minval!')
            if back_to_maindir:
                dirs = os.getcwd().split('/')
                os.chdir('../../')
                os.system('rm -rf '+'/'.join(dirs[-2:]))
            raise Exception(e)
            return
    print('---------------------------------------------------------')
    
    if minval is None:
        minval = yfit[4]
    
    x_bd = x[x_match]
    yp = fpyl.deriv(y,x)
    yp_bd = yp[x_match]
    
    if match:
        print('\nMatching profile and its derivative with fit at boundary value...')
        
        b_new,c_new = fsolve(equations,(yfit[1],yfit[2]),args=(x_bd,yfit[0],yfit[3],y[x_match],yp_bd,minval),xtol=1.0e-10)
        yfit[1] = b_new
        yfit[2] = c_new
        
        print('...fit parameters corrected to:')
        print('---------------------------------------------------------')
        print('Fit width: '+str(1.0/b_new))
        print('Fit height: '+str(c_new))
        print('---------------------------------------------------------')
    
    
    
    # Extend profile
    print('Extending profile to psi_n='+str(psimax))
    n = len(x)
    dx = x[x_match]-x[x_match-1]
    
    #print(m)
    
    m = int(round((psimax - x[x_match])/dx))
    x_ext = x[x_match] + (np.arange(m)+1)*dx
    
    #Truncate 
    if psimatch>0:
        x_in = x[:x_match+1]
        y_in = y[:x_match+1]
    else:
        x_in = x
        y_in = y
    
    if minval is None:
        y_ext = tanhfit(x_ext, yfit[0],yfit[1],yfit[2],yfit[3],yfit[4])+minval
    else:
        y_ext = tanhfit(x_ext, yfit[0],yfit[1],yfit[2],yfit[3])+minval
    
    xnew = np.concatenate((x_in,x_ext))
    ynew = np.concatenate((y_in,y_ext))
    if smooth > 0:
        if not match:
            ynew = fpyl.smooth(ynew,w=smooth,nan='replace')
        else:
            raise Exception('Sorry, cannot match and smooth at the same time!')
    ynewp = fpyl.deriv(ynew,xnew)
    
    # Check if derivative of profile extension stays negative everywhere
    #print(fpyl.deriv(y_ext,x_ext))
    if np.all(fpyl.deriv(y_ext,x_ext)<0.0):
        print('Derivative negative in extended region.')
    else:
        fpyl.printwarn('WARNING: Derivative NOT negative everywhere in extended region.')
    
    
    
    # Plot profiles
    
    fig = plt.figure(constrained_layout=True,figsize=(10,5))
    spec2 = gridspec.GridSpec(ncols=2, nrows=2, figure=fig)
    f2_ax1 = fig.add_subplot(spec2[0, 0])
    f2_ax2 = fig.add_subplot(spec2[1, 0])
    f2_ax3 = fig.add_subplot(spec2[0, 1])
    f2_ax4 = fig.add_subplot(spec2[1, 1])
    
    
    f2_ax1.plot(x_in,y_in,c='C0')
    f2_ax1.plot(xnew,ynew,c='C1',ls='--')
    if psimatch>0:
        f2_ax1.plot(x[x_match+1:],y[x_match+1:],c='C7',ls='--')
    f2_ax1.plot(xf,yf,c='C1',lw=0,marker='.')
    
    f2_ax3.plot(x_in,y_in,c='C0')
    f2_ax3.plot(xnew,ynew,c='C1',ls='--')
    f2_ax3.plot(xf,yf,c='C1',lw=0,marker='.')
    
    
    if psimatch>0:
        pltbd = x_match+1
    else:
        pltbd = None
    print(x_in)
    print(yp[:pltbd])
    f2_ax2.plot(x_in,yp[:pltbd],c='C0')
    f2_ax2.plot(xnew,ynewp,c='C1',ls='--')
    
    f2_ax4.plot(x_in,yp[:pltbd],c='C0')
    f2_ax4.plot(xnew,ynewp,c='C1',ls='--')
    
    f2_ax1.set_xlabel(r'$\psi_N$')
    f2_ax1.set_ylabel('profile')
    f2_ax1.grid()
    f2_ax3.set_xlabel(r'$\psi_N$')
    f2_ax3.set_ylabel('profile')
    f2_ax3.set_xlim([fitrange[0],psimax])
    f2_ax3.set_ylim([np.amin(ynew)-(np.abs(0.05*np.amax(ynew))),np.amax(yf)*1.05])
    f2_ax3.grid()
    
    f2_ax2.set_xlabel(r'$\psi_N$')
    f2_ax2.set_ylabel('profile derivative')
    f2_ax2.grid()
    f2_ax4.set_xlabel(r'$\psi_N$')
    f2_ax4.set_ylabel('profile derivative')
    f2_ax4.set_xlim([fitrange[0],psimax])
    f2_ax4.set_ylim([np.amin(ynewp)-(np.abs(0.05*np.amin(ynewp))),np.amax(ynewp)+0.05*np.sign(np.amax(ynewp))*np.amin(ynewp)])
    f2_ax4.grid()
    
    fig.suptitle(filename, size=12)
    
    #psi_n_smooth,prof_smooth = fpyl.ReadTwoColFile('profile_te.extended_smooth')
    #deriv_smooth = fpyl.deriv(prof_smooth,psi_n_smooth)
    #f2_ax1.plot(psi_n_smooth,prof_smooth,c='C2')
    #f2_ax2.plot(psi_n_smooth,deriv_smooth,c='C2')
    
    # Plot fitted tanh
    #xx = np.linspace(fitrange[0],psimax,500)
    #yy = tanhfit(xx, yfit[0],yfit[1],yfit[2],yfit[3])+minval
    #yyp = fpyl.deriv(yy,xx)
    #plt.plot(xx,yy)
    #plt.plot(xx,yyp)
    
    
    
    # Write new profile to file
    outfile = filename+suffix
    with open(outfile, 'w') as f:
        for i in range(len(xnew)):
            f.write('{:12.6f}{:12.6f}'.format(xnew[i],ynew[i]) +'\n')
    print('Profile written to file '+outfile)
    return yfit



# Function to solve in fsolve. It is solved for parameters b and c.
def equations(p,x,a,d,y_bd,yp_bd,minval):
    b,c=p
    f1 = 0.5*c*(1.0+d*(1-x))*(1.0-np.tanh(b*(x-a))) - y_bd + minval
    f2 = -0.5*c*(d*(1.0-np.tanh(b*(x-a))) + (1.0+d*(1.0-x))*b/(np.cosh(b*(x-a))**2)) - yp_bd
    return (f1,f2)


# Function to fit when minval is specified.
def tanhfit(x, a,b,c,d):
    t = np.tanh(b*(x-a))
    s = 0.5*c*(1.+d*(1.-x))
    f = s*(1.-t)
    return f

# Function to fit when minval is not specified.
def tanhfit2(x, a,b,c,d,e):
    t = np.tanh(b*(x-a))
    s = 0.5*c*(1.+d*(1.-x))
    f = s*(1.-t) + e
    return f
