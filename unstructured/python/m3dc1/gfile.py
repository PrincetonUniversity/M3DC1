# gfile: read and plot EFIT g-file
#
# Coded on 02/18/2020 by:
# Andreas Kleiner:    akleiner@pppl.gov

import math
import re
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from matplotlib import path
#from scipy.io import loadmat
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
from scipy.integrate import trapz
from scipy.integrate import cumtrapz
import m3dc1.fpylib as fpyl


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

def read_gfile(fname,quiet=False):
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
        if gfile_data.fittype == 'FREEGS':
            gfile_data.shot = 0#int(temp[3])
            gfile_data.time = float(temp[-4][:-2]) #in ms
            gfile_data.idum = int(temp[-3])
            gfile_data.nw = int(temp[-2])
            gfile_data.nh = int(temp[-1])
        else:
            if len(temp)>=7:
                gfile_data.shot = 'n/a'
                gfile_data.time = 'n/a'
                gfile_data.idum = int(temp[-3])
                gfile_data.nw = int(temp[-2])
                gfile_data.nh = int(temp[-1])
            else:
                gfile_data.shot = int(temp[2][1:])
                gfile_data.time = float(temp[3][:-2]) #in ms
                gfile_data.idum = int(temp[4])
                gfile_data.nw = int(temp[-2])
                gfile_data.nh = int(temp[-1])
    else:
        if gfile_data.fittype == 'EFIT++':
            date_shot_time = temp[1]
            
            gfile_data.date = temp[1][:8]
            
            shotnum = re.search('#(.*)-', temp[1])
            gfile_data.shot = int(shotnum.group(1))
            timestr = re.search('-(.*)ms', temp[1])
            gfile_data.time = float(timestr.group(1)) #in ms
            
            gfile_data.idum = int(temp[2])
            gfile_data.nw = int(temp[-2])
            gfile_data.nh = int(temp[-1])
            if not quiet:
                print(timestr.group(1), shotnum.group(1), gfile_data.date)
            
        else:
            try:
                gfile_data.idum = int(temp[-3])
                gfile_data.nw = int(temp[-2])
                gfile_data.nh = int(temp[-1])
            except:
                fpylib.printerr('Error reading g-file in line 1!')
                return
    
    #temp = data[1].split()
    temp = fpyl.read_floats(data[1])
    gfile_data.rdim = float(temp[0])
    gfile_data.zdim =  float(temp[1])
    gfile_data.rcentr =  float(temp[2])
    gfile_data.rleft =  float(temp[3])
    gfile_data.zmid =  float(temp[4])
    gfile_data.rg = np.linspace(gfile_data.rleft, gfile_data.rleft+gfile_data.rdim, gfile_data.nw)
    gfile_data.zg = np.linspace(-gfile_data.zdim/2.0 - gfile_data.zmid, gfile_data.zdim/2.0 + gfile_data.zmid, gfile_data.nh)
    
    if not quiet:
        print('rmin: '+str(gfile_data.rleft)+' , rmax: '+str(gfile_data.rleft+gfile_data.rdim))
        print('zmin: '+str(-gfile_data.zdim/2.0 - gfile_data.zmid)+' , zmax: '+str(gfile_data.zdim/2.0 + gfile_data.zmid))
    
    temp=fpyl.read_floats(data[2])
    
    gfile_data.rmaxis = float(temp[0])
    gfile_data.zmaxis = float(temp[1])
    gfile_data.simag = float(temp[2])
    gfile_data.sibry = float(temp[3])
    gfile_data.bcentr = float(temp[4])

    gfile_data.current = fpyl.read_floats(data[3])[0]
    
    
    flatten = lambda l: [item for sublist in l for item in sublist]
    nl = int(math.ceil(gfile_data.nw/5)) # number of lines in g-file containing nw values
    def read_list(start,nlins):
        l = []
        for i in range(nlins):
            l.append(fpyl.read_floats(data[start+i]))
        l = np.asarray(flatten(l))
        return l
    
    
    # Skip line 5, fpol starts in line 6
    gfile_data.fpol = read_list(5,nl)
    gfile_data.Ipol = 2.0*math.pi/(math.pi*4.0e-7)*np.abs(gfile_data.fpol - gfile_data.fpol[0])
    
    gfile_data.pres = read_list(5+nl,nl)
    gfile_data.ffprim = read_list(5+2*nl,nl)
    gfile_data.pprim = read_list(5+3*nl,nl)
    
    #read psi on R-Z grid:
    nlbig = int(math.ceil(gfile_data.nw*gfile_data.nh/5)) # number of lines in g-file containing nw*nh values
    psirz = read_list(5+4*nl,nlbig)
    #print(psirz)
    psirz = np.reshape(psirz,(gfile_data.nh,gfile_data.nw))
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
        fpyl.printwarn('WARNING: Did not read kvtor, rvtor and nmass. This data may be missing in the g-file.')
    
    if not quiet:
        rmin_lcfs = np.amin(gfile_data.rbbbs)
        rmax_lcfs = np.amax(gfile_data.rbbbs)
        R_geo=(rmax_lcfs+rmin_lcfs)/2
        a = (rmax_lcfs-rmin_lcfs)/2
        Z_max = np.amax(gfile_data.zbbbs)
        Z_max_ind = fpyl.get_ind_at_val(gfile_data.zbbbs, Z_max, unique=True)
        Z_min = np.amin(gfile_data.zbbbs)
        Z_min_ind = fpyl.get_ind_at_val(gfile_data.zbbbs, Z_min, unique=True)
        R_upper = gfile_data.rbbbs[Z_max_ind]
        R_lower = gfile_data.rbbbs[Z_min_ind]
        print('a         =  '+str(a))
        print('R0        =  '+str(R_geo))
        print('A = R0/a  =  '+str((R_geo)/(a)))
        print('delta_u   =  '+str((R_geo-R_upper)/a))
        print('delta_l   =  '+str((R_geo-R_lower)/a))
        print('delta     =  '+str( ( (R_geo-R_upper)/a + (R_geo-R_lower)/a ) / 2))
        print('kappa     =  '+str((Z_max-Z_min)/(2*a)))
    #ToDo: Implement rotation stuff
    
    gfile_data.psin = np.linspace(0,1,gfile_data.nw)
    gfile_data.amin = (np.max(gfile_data.rbbbs) - np.min(gfile_data.rbbbs))/2.0;
    
    return gfile_data



#-------------------------------------------
# Plot g-file
#-------------------------------------------

def plot_gfile(fname,fignum=None):
    gfdat = read_gfile(fname)
    
    fa = flux_average_gfile(gfdat)
    
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



def plot_flux(fname,vrange=[],fignum=None,cmap='jet'):
    gfdat = read_gfile(fname)
    
    #x = loadmat('/u/akleiner/codes/matlab/efit/nstx_obj_12nov_6565.mat')
    
    fig = plt.figure(num=fignum)
    ax = plt.gca()
    
    if len(vrange)>0:
        cont = ax.contourf(gfdat.rg,gfdat.zg,gfdat.psirz,np.linspace(vrange[0],vrange[1],200),cmap=cmap)
    else:
        cont = ax.contourf(gfdat.rg,gfdat.zg,gfdat.psirz,200,cmap=cmap)
    
    ax.contour(gfdat.rg,gfdat.zg,gfdat.psirzn,np.linspace(1.1,np.amax(gfdat.psirzn),25),linewidths=0.7,colors='C0')
    ax.plot(gfdat.rmaxis,gfdat.zmaxis,lw=0,marker='+',markersize=10,color='C0')
    ax.plot(gfdat.rbbbs,gfdat.zbbbs,color='m')
    ax.plot(gfdat.rlim,gfdat.zlim,color='C1')
    cbar = fig.colorbar(cont,ax=ax)
    #cbar = fig.colorbar(cont,ax=ax)
    ax.set_aspect('equal')
    ax.set_xlabel('R/m')
    ax.set_ylabel('Z/m')
    ax.grid(True,zorder=10,alpha=0.5)
    
    plt.title('#'+str(gfdat.shot)+', t='+str(gfdat.time)+'ms')
    
    return



def plot_jphi(fname,vrange=[],fignum=None,cmap='jet'):
    gfdat = read_gfile(fname)
    
    dpsidr = np.zeros_like(gfdat.psirz)
    dpsi2dr2 = np.zeros_like(gfdat.psirz)
    dpsidz = np.zeros_like(gfdat.psirz)
    dpsi2dz2 = np.zeros_like(gfdat.psirz)
    
    for i in range(gfdat.psirz.shape[0]):
        dpsidr[i,:] = fpyl.deriv(gfdat.psirz[i,:],gfdat.rg)
        dpsi2dr2[i,:] = fpyl.deriv(dpsidr[i,:],gfdat.rg)
    
    for i in range(gfdat.psirz.shape[1]):
        dpsidz[:,i] = fpyl.deriv(gfdat.psirz[:,i],gfdat.zg)
        dpsi2dz2[:,i] = fpyl.deriv(dpsidz[:,i],gfdat.zg)
    
    R, Z    = np.meshgrid(gfdat.rg, gfdat.zg)
    
    
    left = 1.0/R**2 * dpsidr + 1.0/R * dpsi2dr2
    
    jphi = left + 1.0/R * dpsi2dz2
    
    fig = plt.figure(num=fignum)
    ax = plt.gca()
    
    if len(vrange)>0:
        cont = ax.contourf(gfdat.rg,gfdat.zg,jphi,np.linspace(vrange[0],vrange[1],200),cmap=cmap)
    else:
        cont = ax.contourf(gfdat.rg,gfdat.zg,jphi,200,cmap=cmap)
    
    ax.contour(gfdat.rg,gfdat.zg,gfdat.psirzn,np.linspace(1.1,np.amax(gfdat.psirzn),25),linewidths=0.7,colors='C0')
    ax.plot(gfdat.rmaxis,gfdat.zmaxis,lw=0,marker='+',markersize=10,color='C0')
    ax.plot(gfdat.rbbbs,gfdat.zbbbs,color='m')
    ax.plot(gfdat.rlim,gfdat.zlim,color='C1')
    cbar = fig.colorbar(cont,ax=ax)
    #cbar = fig.colorbar(cont,ax=ax)
    ax.set_aspect('equal')
    ax.set_xlabel('R/m')
    ax.set_ylabel('Z/m')
    ax.grid(True,zorder=10,alpha=0.5)
    
    plt.title('#'+str(gfdat.shot)+', t='+str(gfdat.time)+'ms')
    
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





class Fluxavg():
    def __init__(self):
        self.values = None
        self.Bdl = None
        self.Bp = None
        self.Ba = None
        self.Vprime = None
        self.B = None
        self.B_1 = None
        self.B2 = None
        self.Bp2 = None
        self.R02_R2 = None
        self.R2_1 = None
        self.R0_R = None
        self.A = None
        self.qpsi = None
        self.jpar = None
        self.jpar2 = None
        self.jtor = None
        self.jps = None
        self.jefit = None
        self.jpol = None
        self.jelite = None
        self.Vp = None
        self.V = None
        self.alpha = None
        self.qprime = None
        self.s = None
        self.wp = None
        self.betap = None
        self.betat = None
        self.beta = None
        self.betan = None
        
        self.Rmin = None
        self.Rmax = None
        self.q95 = None


#@numba.jit
def PolyArea(x,y):
    return 0.5*np.abs(np.dot(np.ascontiguousarray(x),np.roll(y,1))-np.dot(np.ascontiguousarray(y),np.roll(x,1)))


def flux_average_gfile(gf):
    """
    Calculates the flux average of a quantity
    
    Arguments:

    **field**
    Name of the field to flux average
    """
    muo = math.pi*4e-7
    
    dr = np.mean(np.diff(gf.rg))
    dz = np.mean(np.diff(gf.zg))
    R,Z = np.meshgrid(gf.rg,gf.zg)
    Rflat = R.ravel()
    Zflat = Z.ravel()
    RZpoints = np.vstack((Rflat,Zflat)).T
    
    bbbs = np.vstack((gf.rbbbs,gf.zbbbs)).T
    
    p = path.Path(bbbs)
    ind = np.reshape(p.contains_points(RZpoints),(gf.nh,gf.nh))
    #print(ind)
    dpsidr,dpsidz = np.gradient(gf.psirz,dr,dz)
    
    nans = np.empty((gf.nh,gf.nh))
    nans[:,:] = np.nan
    zeros = np.zeros((gf.nh,gf.nh))
    
    pprime = np.reshape(griddata(gf.psin,gf.pprim,gf.psirzn),(gf.nh,gf.nh))
    pprime = np.where(ind==True,pprime,nans)
    
    ffprime = np.reshape(griddata(gf.psin,gf.ffprim,gf.psirzn),(gf.nh,gf.nh))
    ffprime = np.where(ind==True,ffprime,nans)
    
    fpol = np.reshape(griddata(gf.psin,gf.fpol,gf.psirzn),(gf.nh,gf.nh))
    fpol = np.where(ind==True,fpol,nans)
    
    #plt.figure()
    #plt.contourf(gf.rg,gf.zg,pprime)
    
    dfpoldr,dfpoldz = np.gradient(fpol,dr,dz)
    
    
    jr = -1.0/muo/R*dfpoldz
    jr = np.where(ind==True,jr,zeros)
    jz =  1.0/muo/R*dfpoldr
    jz = np.where(ind==True,jz,zeros)
    jphi = -np.sign(gf.current)*(R*pprime + 1.0/muo/R*ffprime)
    jphi = np.where(ind==True,jphi,zeros)
    
    BR =  np.sign(gf.current)/R*dpsidz
    BZ = -np.sign(gf.current)/R*dpsidr
    Bpol = np.sqrt(BR**2 + BZ**2)
    
    
    
    fig = plt.figure()
    values = gf.psin[1:] #np.linspace(0.0,1.0)
    #values[0] = values[1]/2
    cs = plt.contour(gf.rg,gf.zg,gf.psirzn,levels=values)
    n_points = len(cs.collections)
    
    if n_points != len(values):
        print('ERROR: not enough contour lines')
        return
    
    points=[]
    #Create a nested list of points of each contour line
    # cs.collections[i] is the i-th contour line
    # cs.collections[i].get_paths()[j] is the j-th connected path of the i-th contour line
    for i in range(n_points-1):
        paths = cs.collections[i].get_paths()
        pp = cs.collections[i].get_paths()[0]
        v = pp.vertices
        # Merge disconnected path to one list of points
        if len(paths)>1:
            for j in range(1,len(paths)):
                pp = cs.collections[i].get_paths()[j]
                v = np.concatenate((v,pp.vertices))
        ind = p.contains_points(v)
        temp = []
        # Append points as list, so indexing is more straightforward
        for k in range(len(ind)):
            if ind[k]==True:
                temp.append(list(v[k,:]))
        points.append(np.asarray(temp))
    points.append(bbbs)
    #plt.figure()
    #for i in range(n_points):
    #    plt.plot(points[i][:,0],points[i][:,1])
    plt.close(fig)
    
    jpar = np.zeros(n_points)
    jpar2 = np.zeros(n_points)
    jtor = np.zeros(n_points)
    jpol = np.zeros(n_points)
    jps = np.zeros(n_points)
    jefit = np.zeros(n_points)
    jelite = np.zeros(n_points)
    
    Rgeo = np.zeros(n_points)
    Bdl = np.zeros(n_points)
    Bp = np.zeros(n_points)
    Ba = np.zeros(n_points)
    Vprime = np.zeros(n_points)
    B = np.zeros(n_points)
    B_1 = np.zeros(n_points)
    B_2 = np.zeros(n_points)
    B2 = np.zeros(n_points)
    Bp2 = np.zeros(n_points)
    R02_R2 = np.zeros(n_points)
    R2_1 = np.zeros(n_points)
    R0_R = np.zeros(n_points)
    A = np.zeros(n_points)
    qpsi = np.zeros(n_points)
    
    for ii in range(n_points):
        Rgeo[ii] = (np.amin(points[ii][:,0]) + np.amax(points[ii][:,0]))/2
        dl = [np.sqrt((points[ii][i+1,0]-points[ii][i,0])**2 + (points[ii][i+1,1]-points[ii][i,1])**2) for i in range(len(points[ii])-1)]
        dl.insert(0, 0.0)
        
        #print(dl)
        l = np.cumsum(dl)
        
        # Averages
        #ToDo: griddata takes about 0.47s to execute, and thus slows down the code. Think about alternative
        bpol = griddata(RZpoints,Bpol.flatten(),points[ii])
        #bphi = griddata(gf.psin,gf.fpol,values[ii])/points[ii][:,0]
        #pprime_fs = griddata(gf.psin,gf.pprim,values[ii])
        #ffprime = griddata(gf.psin,gf.ffprim,values[ii])
        #fpol = griddata(gf.psin,gf.fpol,values[ii])
        bphi = gf.fpol[ii+1]/points[ii][:,0]
        pprime_fs = gf.pprim[ii+1]
        ffprime = gf.ffprim[ii+1]
        fpol = gf.fpol[ii+1]
    
        Norm = trapz(1.0/bpol,l)
        #print(ii,Norm)
        Bdl[ii] = trapz(bpol,l)
        Bp[ii] = trapz(bpol/bpol,l)/Norm
        Ba[ii] = 4e-7*math.pi*gf.current/trapz(bpol/bpol,l)
        Vprime[ii] = 2.0*math.pi*Norm
        B[ii] = trapz(np.sqrt(bpol**2 + bphi**2)/bpol,l)/Norm
        B_1[ii] = trapz(1./np.sqrt(bpol**2 + bphi**2)/bpol,l)/Norm
        B2[ii] = trapz((bpol**2 + bphi**2)/bpol,l)/Norm
        B_2[ii] = trapz(1./(bpol**2 + bphi**2)/bpol,l)/Norm
        Bp2[ii] = trapz(bpol**2/bpol,l)/Norm
        R02_R2[ii] = trapz(((gf.rcentr/points[ii][:,0])**2)/bpol,l)/Norm
        R2_1[ii] = trapz(((1./points[ii][:,0])**2)/bpol,l)/Norm
        R0_R[ii] = trapz((gf.rcentr/points[ii][:,0])/bpol,l)/Norm
        A[ii] = PolyArea(points[ii][:,0],points[ii][:,1])
        qpsi[ii] = np.sign(gf.bcentr)*fpol*Vprime[ii]/(2*math.pi)**2*R2_1[ii]
        
        # Current densities
        jpar[ii]   = abs(fpol*pprime_fs*B_1[ii] + ffprime/fpol/4e-7/math.pi*B[ii])
        jpar2[ii]  = abs(gf.rcentr*pprime_fs + ffprime/gf.rcentr/4e-7/math.pi/((fpol/gf.rcentr)**2)*B2[ii])
        jtor[ii]   = abs(gf.rcentr*pprime_fs + ffprime/gf.rcentr/4e-7/math.pi*R02_R2[ii])
        jps[ii]    = abs(fpol*pprime_fs*(B_1[ii] - B[ii]/B2[ii]))
        jefit[ii]  = abs(jtor[ii]/R0_R[ii])
        jpol[ii]   = abs(jpar2[ii] - jtor[ii])
        
    
    
    Vp = Vprime/2/math.pi
    
    V = cumtrapz(Vprime,values,initial=0)*(gf.sibry - gf.simag)
    
    alpha  = 2e-7/math.pi*Vprime*abs(gf.pprim[1:])*np.sqrt(V/2/math.pi/math.pi/gf.rcentr)
    #print('pprim')
    #print(gf.pprim[1:])
    #print('Vprime')
    #plt.figure()
    #plt.plot(values,Vprime)
    #print(Vprime)
    #print('V')
    #plt.figure()
    #plt.plot(values,V)
    #print(V)
    #qprime = elib.deriv(qpsi/(gf.sibry - gf.simag),values)
    startTime = time.time()
    qprime = fpyl.deriv(gf.qpsi[1:],gf.psin[1:])
    executionTime = (time.time() - startTime)
    print('Execution time in seconds: ' + str(executionTime))
    
    #s = 2.0*V/Vprime*qprime/qpsi
    s = qprime*gf.psin[1:]/gf.qpsi[1:]
    wp = 3/2*trapz(gf.pres[1:],V)
    betap = trapz(gf.pres[1:],V)/V[-1]*2*4e-7*math.pi/Ba[-1]**2
    betat = 100*trapz(gf.pres[1:],V)/V[-1]*2*4e-7*math.pi/gf.bcentr**2
    beta  = betap*betat/(betap + betat/100)
    betan = abs(beta*gf.amin*gf.bcentr/(gf.current/1e6))
    
    jelite = -np.sign(gf.current)*(gf.rcentr*gf.pprim[1:]*(fpol/gf.rcentr)**2.*B_2 + ffprime/muo/gf.rcentr)
    jeliten = (np.amax(jelite[np.where(gf.psin>0.8)[0][0]:]) + jelite[-1])/2/(gf.current/A[-1])
    #print(jelite[np.where(gf.psin>0.8)[0][0]:],jelite[-1],gf.current,A[-1])
    #print(jeliten)
    #plt.figure(10)
    #plt.plot(values,jelite)
    
    #plt.figure(11)
    #plt.plot(values,s)
    
    #plt.figure(12)
    #plt.plot(values,qpsi)

    #plt.figure(13)
    #plt.plot(values,jpar)
    
    #plt.figure(14)
    #plt.plot(values,jtor)
    
    fa = Fluxavg()
    fa.values = values
    fa.Bdl = Bdl
    fa.Bp = Bp
    fa.Ba = Ba
    fa.Vprime = Vprime
    fa.B = B
    fa.B_1 = B_1
    fa.B2 = B2
    fa.Bp2 = Bp2
    fa.R02_R2 = R02_R2
    fa.R2_1 = R2_1
    fa.R0_R = R0_R
    fa.A = A
    fa.qpsi = qpsi
    fa.jpar = jpar
    fa.jpar2 = jpar2
    fa.jtor = jtor
    fa.jps = jps
    fa.jefit = jefit
    fa.jpol = jpol
    fa.jelite = jelite
    fa.Vp = Vp
    fa.V = V
    fa.alpha = alpha
    fa.qprime = qprime
    fa.s = s
    fa.wp = wp
    fa.betap = betap
    fa.betat = betat
    fa.beta = beta
    fa.betan = betan
    
    fa.Rmin = np.amin(gf.rbbbs)
    fa.Rmax = np.amax(gf.rbbbs)
    q_interp = interp1d(values, qpsi,kind='cubic',fill_value="extrapolate")
    fa.q95 = q_interp(0.95)
    
    return fa

