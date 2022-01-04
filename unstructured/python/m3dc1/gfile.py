# gfile: read and plot EFIT g-file
#
# Coded on 02/18/2020 by:
# Andreas Kleiner:    akleiner@pppl.gov

import math
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from scipy.io import loadmat
from scipy.interpolate import interp1d
import efit.efitlib as elib
from efit.flux_average import flux_average


#-------------------------------------------
# Class definitions used to process g-file
#-------------------------------------------

class Gfile():
    def __init__(self):
        self.shot = None
        self.time = None
        self.fittype = None
        self.date = None
        self.idum = None
        self.nw = None
        self.nh = None
        self.rdim = None
        self.zdim =  None
        self.rcentr =  None
        self.rleft =  None
        self.zmid =  None
        self.rg = None
        self.zg = None
        self.rmaxis = None
        self.zmaxis = None
        self.simag = None
        self.sibry = None
        self.bcentr = None
        self.current = None
        self.fpol = None
        self.Ipol = None
        self.pres = None
        self.ffprim = None
        self.pprim = None
        self.psirz = None
        self.psirzn = None
        self.qpsi = None
        self.nbbbs = None
        self.limitr = None
        self.rbbbs = None
        self.zbbbs = None
        self.rlim = None
        self.zlim = None
        self.kvtor = None
        self.rvtor = None
        self.nmass = None
        self.psin = None
        self.amin = None


#-------------------------------------------
# Read one or multiple g-files
#-------------------------------------------

def read_gfile(fname):
    """
    Reads information from EFIT g-file
    
    Arguments:

    **fname**
    Name of file that will be read, i.e. "g123456.00789"
    """
    f = open(fname, 'r')
    data = f.readlines()
    f.close()
    nlines = len(data)
    
    gfile_data = Gfile()
    
    temp = data[0].split()
    gfile_data.fittype = temp[0]
    if len(temp)>5:
        gfile_data.date = temp[1]
        gfile_data.shot = int(temp[2][1:])
        gfile_data.time = float(temp[3][:-2]) #in ms
    
        gfile_data.idum = int(temp[4])
        gfile_data.nw = int(temp[5])
        gfile_data.nh = int(temp[6])
    else:
        if gfile_data.fittype == 'EFIT++':
            date_shot_time = temp[1]
            
            gfile_data.date = temp[1][:8]
            
            shotnum = re.search('#(.*)-', temp[1])
            gfile_data.shot = int(shotnum.group(1))
            timestr = re.search('-(.*)ms', temp[1])
            gfile_data.time = float(timestr.group(1)) #in ms
            
            gfile_data.idum = int(temp[2])
            gfile_data.nw = int(temp[3])
            gfile_data.nh = int(temp[4])
            
            print(timestr.group(1), shotnum.group(1), gfile_data.date)
            
        else:
            try:
                gfile_data.idum = int(temp[-3])
                gfile_data.nw = int(temp[-2])
                gfile_data.nh = int(temp[-1])
            except:
                elib.printerr('Error reading g-file in line 1!')
    
    #temp = data[1].split()
    temp = elib.read_floats(data[1])
    gfile_data.rdim = float(temp[0])
    gfile_data.zdim =  float(temp[1])
    gfile_data.rcentr =  float(temp[2])
    gfile_data.rleft =  float(temp[3])
    gfile_data.zmid =  float(temp[4])
    gfile_data.rg = np.linspace(gfile_data.rleft, gfile_data.rleft+gfile_data.rdim, gfile_data.nw)
    gfile_data.zg = np.linspace(-gfile_data.zdim/2.0 - gfile_data.zmid, gfile_data.zdim/2.0 + gfile_data.zmid, gfile_data.nw)
    
    print('rmin: '+str(gfile_data.rleft)+' , rmax: '+str(gfile_data.rleft+gfile_data.rdim))
    print('zmin: '+str(-gfile_data.zdim/2.0 - gfile_data.zmid)+' , zmax: '+str(gfile_data.zdim/2.0 + gfile_data.zmid))
    
    temp=elib.read_floats(data[2])
    
    gfile_data.rmaxis = float(temp[0])
    gfile_data.zmaxis = float(temp[1])
    gfile_data.simag = float(temp[2])
    gfile_data.sibry = float(temp[3])
    gfile_data.bcentr = float(temp[4])

    gfile_data.current = elib.read_floats(data[3])[0]
    
    
    flatten = lambda l: [item for sublist in l for item in sublist]
    nl = int(math.ceil(gfile_data.nw/5)) # number of lines in g-file containing nw values
    def read_list(start,nlins):
        l = []
        for i in range(nlins):
            l.append(elib.read_floats(data[start+i]))
        l = np.asarray(flatten(l))
        return l
    
    
    # Skip line 5, fpol starts in line 6
    gfile_data.fpol = read_list(5,nl)
    gfile_data.Ipol = 2.0*math.pi/(math.pi*4.0e-7)*np.abs(gfile_data.fpol - gfile_data.fpol[0])
    
    gfile_data.pres = read_list(5+nl,nl)
    gfile_data.ffprim = read_list(5+2*nl,nl)
    gfile_data.pprim = read_list(5+3*nl,nl)
    
    nlbig = int(math.ceil(gfile_data.nw*gfile_data.nh/5)) # number of lines in g-file containing nw*nh values
    psirz = read_list(5+4*nl,nlbig)
    psirz = np.reshape(psirz,(gfile_data.nw,gfile_data.nh))
    gfile_data.psirz = psirz#.transpose()
    gfile_data.psirzn = (gfile_data.psirz - gfile_data.simag)/(gfile_data.sibry - gfile_data.simag) # Normalized stream function
    
    gfile_data.qpsi = read_list(5+4*nl+nlbig,nl)
    
    pos = 5+5*nl+nlbig
    temp = data[pos].split()
    gfile_data.nbbbs = int(temp[0])
    gfile_data.limitr = int(temp[1])
    
    nlbdry = int(math.ceil(2*gfile_data.nbbbs/5)) # number of lines containing boundary values
    temp = read_list(pos+1,nlbdry)
    temp = np.reshape(temp,(gfile_data.nbbbs,2)).transpose()
    gfile_data.rbbbs = temp[0,:]
    gfile_data.zbbbs = temp[1,:]
    
    pos = pos+1+nlbdry
    nllimitr = int(math.ceil(2*gfile_data.limitr/5)) # number of lines containing limiter values
    temp = read_list(pos,nllimitr)
    temp = np.reshape(temp,(gfile_data.limitr,2)).transpose()
    gfile_data.rlim = temp[0,:]
    gfile_data.zlim = temp[1,:]
    
    try:
        pos = pos+nllimitr
        #print(temp)
        temp = data[pos].split()
        gfile_data.kvtor = int(temp[0])
        gfile_data.rvtor = float(temp[1])
        gfile_data.nmass = int(temp[2])
    except:
        elib.printwarn('WARNING: Did not read kvtor, rvtor and nmass. This data may be missing in the g-file.')
    
    #ToDo: Implement rotation stuff
    
    gfile_data.psin = np.linspace(0,1,gfile_data.nw)
    gfile_data.amin = (np.max(gfile_data.rbbbs) - np.min(gfile_data.rbbbs))/2.0;
    
    return gfile_data



#-------------------------------------------
# Plot g-file
#-------------------------------------------

def plot_gfile(fname,fignum=None):
    gfdat = read_gfile(fname)
    
    fa = flux_average(gfdat)
    
    #x = loadmat('/u/akleiner/codes/matlab/efit/nstx_obj_12nov_6565.mat')
    
    fig = plt.figure(constrained_layout=True,figsize=(15,7),num=fignum)
    spec2 = gridspec.GridSpec(ncols=3, nrows=2, figure=fig)
    f_ax1 = fig.add_subplot(spec2[:, 0])
    f_ax2 = fig.add_subplot(spec2[0, 1])
    f_ax3 = fig.add_subplot(spec2[0, 2])
    f_ax4 = fig.add_subplot(spec2[1, 1])
    f_ax5 = fig.add_subplot(spec2[1, 2])
    
    
    cont = f_ax1.contour(gfdat.rg,gfdat.zg,gfdat.psirzn,np.linspace(0.1,0.9,9),linewidths=0.7,colors='C0')
    f_ax1.contour(gfdat.rg,gfdat.zg,gfdat.psirzn,np.linspace(1.1,np.amax(gfdat.psirzn),25),linewidths=0.7,colors='C0')
    f_ax1.plot(gfdat.rmaxis,gfdat.zmaxis,lw=0,marker='+',markersize=10,color='C0')
    f_ax1.plot(gfdat.rbbbs,gfdat.zbbbs,color='m')
    f_ax1.plot(gfdat.rlim,gfdat.zlim,color='C1')
    #cbar = fig.colorbar(cont,ax=ax)
    f_ax1.set_aspect('equal')
    f_ax1.set_xlabel('R/m')
    f_ax1.set_ylabel('Z/m')
    f_ax1.grid(True,zorder=10,alpha=0.5)
    
    f_ax2.plot(gfdat.psin,gfdat.pprim/1e6*gfdat.rcentr,lw=2,label=r"$R_0 p'$ (MA/m$^2$)")
    f_ax2.plot(gfdat.psin,gfdat.ffprim/(math.pi*4e-7)/1e6/gfdat.rcentr,lw=2,label=r"$FF'/\mu_0 R_0$ (MA/m$^2$)")
    f_ax2.set_xlabel(r'$\psi_N$')
    f_ax2.grid(True,zorder=10,alpha=0.5)
    f_ax2.legend(loc=0)
    
    f_ax3.plot(gfdat.psin,gfdat.pres/1e3,lw=2)
    f_ax3.set_xlabel(r'$\psi_N$')
    f_ax3.set_ylabel(r'Pressure (kPa)')
    f_ax3.grid(True,zorder=10,alpha=0.5)
    
    f_ax4.plot(gfdat.psin,gfdat.qpsi,lw=2,label=r"$q$")
    #f_ax4.plot(fa.values,fa.qpsi,lw=2,label=r"$q$")
    f_ax4.plot(fa.values,fa.s,lw=2,label=r"$s$")
    f_ax4.plot(fa.values,fa.alpha,lw=2,label=r"$\alpha$")
    f_ax4.set_xlabel(r'$\psi_N$')
    f_ax4.grid(True,zorder=10,alpha=0.5)
    f_ax4.legend(loc=0)
    
    f_ax5.plot(fa.values,fa.jpar/1e6,lw=2,label=r"$j_{\parallel}$ (MA/m$^2$)")
    f_ax5.plot(fa.values,fa.jtor/1e6,lw=2,label=r"$j_{tor}$ (MA/m$^2$)")
    f_ax5.plot(gfdat.psin,gfdat.Ipol/1e6,lw=2,label=r"$I_{pol}$ (MA/10)")
    f_ax5.set_xlabel(r'$\psi_N$')
    f_ax5.grid(True,zorder=10,alpha=0.5)
    f_ax5.legend(loc=0)
    
    fig.suptitle('#'+str(gfdat.shot)+', t='+str(gfdat.time)+'ms')
    
    return


def plot_gfile_p(fname,fignum=None,show_legend=True):
    gfdat = read_gfile(fname)
    
    plt.figure(num=fignum)
    plt.plot(gfdat.psin,gfdat.pres/1e3,lw=2,label=fname)
    ax = plt.gca()
    ax.set_xlabel(r'$\psi_N$')
    ax.set_ylabel(r'Pressure (kPa)')

    if show_legend==True:
        plt.legend(loc=0)
    plt.grid(True)
    #ax.tick_params(axis='both', which='major', labelsize=ticklblfs)
    #plt.title(lbl)
    plt.tight_layout()
    return

def write_wall(fname,outfile):
    gfdat = read_gfile(fname)
    wall = np.column_stack([gfdat.rlim,gfdat.zlim])
    plt.figure()
    plt.plot(gfdat.rlim,gfdat.zlim)
    with open(outfile, 'w') as f:
        f.write('     '+str(len(gfdat.rlim))+'\n')
        for w in wall:
            f.write('     '+"{0:.6f}".format(w[0]) + '     ' + "{0:.6f}".format(w[1])+'\n')
        
    return

def write_p(fname,length=256):
    gfile = read_gfile(fname)
    
    if len(gfile.pres) != length:
        f = interp1d(gfile.psin,gfile.pres)
        psin = np.linspace(0.0, 1.0, num=length, endpoint=True)
        p_new = f(psin)
        temp = np.column_stack((psin,p_new/1000))
    else:
        temp = np.column_stack((gfile.psin,gfile.pres/1000))
    
    for ln in temp: 
        print(" {:8.6f}   {:8.6f}".format(*ln))
    return
