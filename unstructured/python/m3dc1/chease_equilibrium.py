7#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 10 2020

@author: Andreas Kleiner
"""
import os
import glob
import math
import numpy as np
import shutil
from termcolor import colored
#import matplotlib.pyplot as plt

from m3dc1.gfile import read_gfile
import m3dc1.fpylib as fpyl


def convert_eq_chease(res=513,tol=1e-2,quiet=True):
    """
    Uses CHEASE to reculate any GEQDSK-based equilibrium.

    Arguments:

    **res**
    res=NW=NH. Number of grid points in output GEQDSK file.

    **tol**
    Convergence tolerance for pprime. Is error for pprime < tol,
    the equilibrium is considered to be converged to the input
    equilibrium. Errors are also calculated for ffprime and the
    boundary shape. The thresholds for these errors are scaled
    with tol.

    **quiet**
    If True, do not print equilibrium details to screen.
    """
    
    #ToDo: CHECK IF CHEASE MODULE IS LAODED
    
    #main_dir = os.getcwd()
    
    tol_ffprime = tol#10.*tol
    tol_lcfs = 2000.*tol
    
    gfiles = sorted(glob.glob('g*.*'))
    ngf = len(gfiles)
    
    chease_dir = 'chease'
    if not os.path.exists(chease_dir):
        os.mkdir(chease_dir)
        
    os.chdir(chease_dir)
    
    chease_namelist = []
    chease_namelist.append("***\n")
    chease_namelist.append("***\n")
    chease_namelist.append("***\n")
    chease_namelist.append("***\n")
    chease_namelist.append(" &EQDATA\n")
    chease_namelist.append("AFBS = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \n")
    chease_namelist.append("AFBS2 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \n")
    chease_namelist.append("AP = 0.300, 0.400, 0.210, 0, 0.210, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \n")
    chease_namelist.append("AP = -0.060, 1.6, 0.75, 0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \n")
    chease_namelist.append("APLACE = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \n")
    chease_namelist.append("AP2 = 0.100, 0.500, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \n")
    chease_namelist.append("ASPCT = 0.1,\n")
    chease_namelist.append("AT = 0.6, 1.6, 0.75, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \n")
    chease_namelist.append("AWIDTH = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \n")
    chease_namelist.append("BEANS = 0,\n")
    chease_namelist.append("BPLACE = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \n")
    chease_namelist.append("BSFRAC = 0.500,\n")
    chease_namelist.append("BWIDTH = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \n")
    chease_namelist.append("CETA = 0,\n")
    chease_namelist.append("CFBAL = 1,\n")
    chease_namelist.append("CPLACE = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \n")
    chease_namelist.append("CPRESS = 1.0,\n")
    chease_namelist.append("CSSPEC = 0,\n")
    chease_namelist.append("CURRT = 0.0,\n")
    chease_namelist.append("CWIDTH = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \n")
    chease_namelist.append("DELTA = 0,\n")
    chease_namelist.append("DPLACE = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \n")
    chease_namelist.append("DWIDTH = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \n")
    chease_namelist.append("ELONG = 1,\n")
    chease_namelist.append("EPLACE = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \n")
    chease_namelist.append("ETAEI = 1.500,\n")
    chease_namelist.append("EWIDTH = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, \n")
    chease_namelist.append("EPSLON = 1.000E-8,\n")
    chease_namelist.append("PANGLE = 0,\n")
    chease_namelist.append("PREDGE = 0,\n")
    chease_namelist.append("PSISCL = 1,\n")
    chease_namelist.append("QSPEC = 11.1,\n")
    chease_namelist.append("RNU = 0,\n")
    chease_namelist.append("RZION = 1,\n")
    chease_namelist.append("R0 = 1.0\n")
    chease_namelist.append("R0EXP = 1,\n")
    chease_namelist.append("B0EXP = 1.0,\n")
    chease_namelist.append("RZ0 = 0.0,\n")
    chease_namelist.append("SGMA = 0,\n")
    chease_namelist.append("SOLPDA = 0.600,\n")
    chease_namelist.append("SOLPDB = 0,\n")
    chease_namelist.append("SOLPDC = 0,\n")
    chease_namelist.append("SOLPDD = 0,\n")
    chease_namelist.append("SOLPDE = 0,\n")
    chease_namelist.append("THETA0 = 0,\n")
    chease_namelist.append("TRIANG = 0,\n")
    chease_namelist.append("TRIPLT = 0,\n")
    chease_namelist.append("XI = 0,\n")
    chease_namelist.append("NBLC0 = 1,\n")
    chease_namelist.append("NBLOPT = 0,\n")
    chease_namelist.append("NBSFUN = 1,\n")
    chease_namelist.append("NBSOPT = 0,\n")
    chease_namelist.append("NBSTRP = 1,\n")
    chease_namelist.append("NCHI = 200,\n")
    chease_namelist.append("NCSCAL = 2,\n")
    chease_namelist.append("NDIFPS = 0,\n")
    chease_namelist.append("NDIFT = 0,\n")
    chease_namelist.append("NEGP = -1,\n")
    chease_namelist.append("NER = 1,\n")
    chease_namelist.append("NEQDSK = 1,\n")
    chease_namelist.append("NFUNC = 4,\n")
    chease_namelist.append("NIDEAL=5,\n")
    chease_namelist.append("NIPR = 1,\n")
    chease_namelist.append("NISO = 120,\n")
    chease_namelist.append("NMESHA = 0,\n")
    chease_namelist.append("NMESHB = 0,\n")
    chease_namelist.append("NMESHC = 0,\n")
    chease_namelist.append("NMESHD = 0,\n")
    chease_namelist.append("NMESHE = 0,\n")
    chease_namelist.append("NMGAUS = 8,\n")
    chease_namelist.append("nplot=1,\n")
    chease_namelist.append("NPOIDA = 0,\n")
    chease_namelist.append("NPOIDB = 0,\n")
    chease_namelist.append("NPOIDC = 0,\n")
    chease_namelist.append("NPOIDD = 0,\n")
    chease_namelist.append("NPOIDE = 0,\n")
    chease_namelist.append("NPOIDQ = 2,\n")
    chease_namelist.append("NPP = 1,\n")
    chease_namelist.append("NPPFUN = 4,\n")
    chease_namelist.append("NPPR = 30,\n")
    chease_namelist.append("NPSI = 120,\n")
    chease_namelist.append("NZBOX = "+"{:d}".format(res)+",\n")
    chease_namelist.append("NRBOX = "+"{:d}".format(res)+",\n")
    chease_namelist.append("RBOXLFT = 0.2,\n")
    chease_namelist.append("RBOXLEN = 1.3,\n")
    chease_namelist.append("ZBOXMID = 0.0,\n")
    chease_namelist.append("ZBOXLEN = 2.6,\n")
    chease_namelist.append("NRSCAL = 0,\n")
    chease_namelist.append("NS = 120,\n")
    chease_namelist.append("NSGAUS = 4,\n")
    chease_namelist.append("NSOUR = 2,\n")
    chease_namelist.append("NSTTP = 2,\n")
    chease_namelist.append("NSURF = 6,\n")
    chease_namelist.append("NT = 120,\n")
    chease_namelist.append("NTGAUS=4,\n")
    chease_namelist.append("NTMF0 = 1,\n")
    chease_namelist.append("NTNOVA = 64,\n")
    chease_namelist.append("NTURN = 20,\n")
    chease_namelist.append(" /\n")
    chease_namelist.append(" &NEWRUN\n")
    chease_namelist.append(" /\n")
    
    f = open('results.dat', 'a',2)
    f.write("   ID      g-file name      pprime error     ffprime error        LCFS error    converged (tol={:.8e}, tol_ffprime={:.8e}, tol_lcfs={:.8e})\n".format(tol,tol_ffprime,tol_lcfs))
    f.close()
    for i,gf in enumerate(gfiles):
        fpyl.printnote("Calculating equilibrium " + gf + " (" + str(i+1) + "/" + str(ngf) + ") ...")
        shutil.copyfile('../'+gf,'EXPEQ')
        gfin = read_gfile('EXPEQ',quiet=quiet)
        currt = 4*math.pi*1e-7*gfin.current
        
        chease_namelist[23] = "CURRT = "+"{0:.20f}".format(currt)+",\n"
        
        with open('chease_namelist', 'w') as namelist:
            for line in chease_namelist:
                namelist.write(line)
        
        os.system('chease > out'+ gf)
        shutil.move('EQDSK_COCOS_02_POS.OUT',gf)
        
        
        #Read output equilibrium for sanity checks
        gfout = read_gfile(gf,quiet=quiet)
        
        #Determine which radial grid to use. It makes sense to use psin of
        #whatever equilibrium has the lowest resolution, as interpolating
        #a low resolution solution and comparing at the functions at these
        #interpolated points artificially inflates the error, even though
        #the functions agree very well on the actual grid points.
        
        
        if gfout.nw != gfin.nw:
            if gfin.nw > gfout.nw:
                psin_grid = gfout.psin
                pprime_in = np.interp(psin_grid,gfin.psin,gfin.pprim)
                ffprime_in = np.interp(psin_grid,gfin.psin,gfin.ffprim)
                pprime_out = gfout.pprim
                ffprime_out = gfout.ffprim
            else:
                psin_grid = gfin.psin
                pprime_in = gfin.pprim
                ffprime_in = gfin.ffprim
                pprime_out = np.interp(psin_grid,gfout.psin,gfout.pprim)
                ffprime_out = np.interp(psin_grid,gfout.psin,gfout.ffprim)
        else:
            psin_grid = gfout.psin
            pprime_in = gfin.pprim
            ffprime_in = gfin.ffprim
            pprime_out = gfout.pprim
            ffprime_out = gfout.ffprim
        
        #Calculate errors and check if output equilibrium is close enough to input
        pprim_err = calculate_error(psin_grid,pprime_in,pprime_out)
        ffprim_err = calculate_error(psin_grid,ffprime_in,ffprime_out)
        lcfs_err = calculate_error_2d(gfin.rbbbs,gfout.rbbbs,gfin.zbbbs,gfout.zbbbs)
        
        converged = True if (pprim_err < tol and ffprim_err < tol_lcfs and lcfs_err < tol_lcfs) else False
        
        errors = [pprim_err, ffprim_err, lcfs_err]
        tols = [tol, tol_ffprime, tol_lcfs]
        
        for j,st in enumerate(['pprime error: ', 'ffprime error: ', 'LCFS error: ']):
            color = 'green' if errors[j] < tols[j] else 'red'
            print(colored(st+str(errors[j]),color))
        
        f = open('results.dat', 'a',2)
        f.write("{0:5d}    {1}    {2:.8e}    {3:.8e}    {4:.8e}    {5:5d}\n".format(i,gf,pprim_err,ffprim_err,lcfs_err,converged))
        print("{0:5d}    {1}    {2:.8e}    {3:.8e}    {4:.8e}    {5:5d}".format(i,gf,pprim_err,ffprim_err,lcfs_err,converged))
        f.close()
        
    #f.close()
        
    os.chdir('../')
    return


def calculate_error(x,y1,y2):
    """
    Calculates the relative error between two functions y1 and y2.
    The relative error is defined as the norm ||y2-y1|| divided
    by the norm ||y1||.

    Arguments:

    **x**
    List of abscissa coordinates. Both functions need to use the
    same values.

    **y1**
    Values of first function.

    **y2**
    Values of second function.
    """
    abs_error = np.sqrt(np.sum((np.asarray(y1) - np.asarray(y2))**2*np.diff(x,append=x[-1])))
    norm = np.sqrt(np.sum(np.asarray(y1)**2*np.diff(x,append=x[-1])))
    return abs_error / norm

def calculate_error_2d(x1,x2,y1,y2):
    """
    Calculates an error estimate for a 2D function, e.g. LCFS by
    adding up the Euclidian distance between each point of the
    input and output curves.

    Arguments:

    **x1**
    List of x values of first data set.

    **x2**
    List of x values of second data set.

    **y1**
    List of y values of first data set.
    
    **y2**
    List of y values of second data set.
    """
    abs_error = 1.0/len(x1)*np.sum(np.sqrt(np.asarray(x2-x1)**2 + np.asarray(y2-y1)**2))
    return abs_error
