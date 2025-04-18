#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on October 04 2023

@author: Andreas Kleiner
"""
import os
import glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from m3dc1.get_time_of_slice import get_time_of_slice


def run_trace(time_slices,dR=0.1, dZ=0.1, pts=51, trans=100):
    """
    Runs trace to create data for Poincare plots for specified time slices.
    
    Arguments:

    **time_slices**
    Array, list or other structure such as range that contains time slices for
    Poincare plot.
    
    **dR**
    R-spacing of seed points
    
    **dZ**
    Z-spacing of seed points
    
    **pts**
    Number of seed points
    
    **trans**
    Toroidal transits per seed point
    """
    
    wd = os.getcwd()
    
    # If the directory does not exist, create it
    if not os.path.exists(wd+'/poincare'):
        os.makedirs('poincare', exist_ok=True)
        os.chdir('poincare')
    
    time_slices = np.atleast_1d(time_slices)
    
    for ts in time_slices:
        tsdir = 'time_'+str(ts).zfill(3)
        os.makedirs(wd+'/poincare/'+tsdir, exist_ok=True)
        os.chdir(wd+'/poincare/'+tsdir)
        cmd_trace = os.environ.get('FIO_INSTALL_DIR')+'/bin/trace -m3dc1 ../../C1.h5 -m3dc1 ../../C1.h5 '+str(ts)+' -dZ '+str(dZ)+' -dR '+str(dR)+' -p '+str(pts)+' -t '+str(trans)
        print(cmd_trace)
        os.system(cmd_trace)
    
    return



def plot_poincare(time=None,fignum=None,pub=False,Rlim=[None,None],Zlim=[None,None],short_title=False,savefig=False,savepath=''):
    """
    Show Poincare plot.

    Arguments:

    **time**
    Time slice

    **fignum**
    Figure number

    **pub**
    True/False. Format figure for publication.

    **Rlim**
    Axis limits in R direction.

    **Zlim**
    Axis limits in Z direction.

    **short_title**
    If True, only time will shown in title. If False, time slice and time will
    be shown.

    **savefig**
    If True, save figure to png and pdf files.

    **savepath**
    Directory path for saved figures.
    """
    
    # Set font sizes and plot style parameters
    if pub:
        axlblfs = 20
        ticklblfs = 18
        titlefs = 18
        linew = 1
    else:
        axlblfs = 12
        ticklblfs = 12
        titlefs = None
        linew = 1
    
    #Change to correct directory if time is specified
    pwd = os.getcwd()
    change_dir=False
    if time is not None:
        if not 'poincare/time_' in pwd:
            change_dir=True
            if 'poincare' in pwd:
                os.chdir('time_'+str(time).zfill(3))
            else:
                os.chdir('poincare/time_'+str(time).zfill(3))
    outfiles = sorted(glob.glob(os.getcwd()+"/out*"))
    
    fig = plt.figure(num=fignum)
    ax = plt.gca()
    
    for of in outfiles:
        if os.path.getsize(of) > 0:
            x,y = np.loadtxt(of, usecols=(1,2), unpack=True)
            plt.plot(x,y,lw=0,marker='.',ms=1)
        #else:
        #    print(of + 'is empty')
    
    plt.grid(True)
    ax.set_aspect('equal',adjustable='box')
    ax.set_xlim(Rlim[0],Rlim[1])
    ax.set_ylim(Zlim[0],Zlim[1])
    plt.xlabel(r'$R$',fontsize=axlblfs)
    plt.ylabel(r'$Z$',fontsize=axlblfs)
    
    if time is None:
        time = pwd.split('/')[-1][-3:]

    unitlabel = 'ms'
    figtitle = 'time='+"{0:.2f}".format(get_time_of_slice(int(time),filename='../../C1.h5',units='mks',millisec=True,quiet=True))+' '+unitlabel
    if not short_title:
        figtitle = 'ts: '+time + ',    ' + figtitle
    plt.title(figtitle,fontsize=titlefs)
    plt.tight_layout()
    

    
    if savefig:
        plt.savefig(savepath+'time_'+str(time).zfill(3)+'.png', format='png',dpi=900,bbox_inches='tight')
        plt.savefig(savepath+'time_'+str(time).zfill(3)+'.pdf', format='pdf',bbox_inches='tight')
    
    if change_dir:
        os.chdir(pwd)
    return fig



def poincare_movie(pub=False,Rlim=[None,None],Zlim=[None,None],short_title=False):
    """
    Save Poincare plot for all time slices after trace was run. A movie can be created
    subsequently by using ffmpeg, e.g.
    ffmpeg -r 2 -pattern_type glob -i '*.png' -vcodec libx264 -crf 25 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" poincare.mp4

    Arguments:

    **pub**
    True/False. Format figure for publication.

    **Rlim**
    Axis limits in R direction.

    **Zlim**
    Axis limits in Z direction.

    **short_title**
    If True, only time will shown in title. If False, time slice and time will
    be shown.
    """
    pwd = os.getcwd()
    if not 'poincare' in pwd:
        os.chdir('poincare')
    
    if not os.path.isdir('movie'):
        os.mkdir('movie')
    
    ts_dirs = sorted(glob.glob(os.getcwd()+"/time_*"))
    
    for ts in ts_dirs:
        os.chdir(ts)
        plot_poincare(time=None,fignum=None,pub=pub,Rlim=Rlim,Zlim=Zlim,short_title=short_title,savefig=True,savepath='../movie/')
    
    #os.chdir('movie')
    #os.system('ffmpeg -r 2 -i time_%03d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p poincare.mp4')
    #os.system('ffmpeg -r 2 -pattern_type glob -i '*.png' -vcodec libx264 -crf 25 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" poincare.mp4')

    return
