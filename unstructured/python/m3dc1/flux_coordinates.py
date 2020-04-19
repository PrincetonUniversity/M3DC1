#!/usr/bin/env python3
# flux_coordinates: calculates flux coordinates
#
# Based on Nate Ferraro's IDL routines
# Coded on 02/19/2020 by:
# Andreas Kleiner:    akleiner@pppl.gov

import fpy
import math
import numpy as np
from scipy.interpolate import interp1d
import h5py
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import m3dc1.fpylib as fpyl
from m3dc1.read_h5 import readParameter

#ToDo: Allow psin_range to start at a value > 0
def flux_coordinates(sim=None, file_name='C1.h5', time=0, fcoords='', phit=0.0, points=200, fbins=None, tbins=None, itor=None, r0=None, psin_range=None, njac=False, makeplot=False,fignum=501):
    """
    Calculates flux coordinates and returns a fpy.sim_data object containing a flux
    coordinate object as property sim.fc. Possibilities are PEST, Boozer, Hamada coordinates
    and flux coordinates using the geometrical poloidal angle.
    
    Arguments:

    **sim**
    simulation data object. If none is provided, file_name and time must be specified.

    **file_name**
    Name of file that will be read, i.e. "../C1.h5". Used when no sim object is given.

    **time**
    The time-slice that will be used if no sim object is given
    
    **fcoords**
    Name of desired flux coordinate system : 'pest', 'boozer', 'hamada', canonical, ''
    If not specified (i.e. ''), the geometric poloidal angle will be used.
    
    **phit**
    Toroidal angle of the plane where flux coordinates shall be calculated.
    
    **points**
    Number of points in radial and poloidal direction, i.e.
    number of flux surfaces and poloidal points
    
    **fbins**
    Number of radial points (flux surfaces)
    
    **tbins**
    Number of poloidal points
    
    **itor**
    Same as itor in M3DC1
    
    **r0**
    Same as rzero in M3DC1
    
    **psin_range**
    Range of normalized psi in which the flux coordinates are calculated.
    Has to start at 0.0
    
    **njac**
    Calculate Jacobian numerically
    
    **makeplot**
    Show plot flux coordinates
    
    **fignum**
    Figure number for flux coordinate plot
    """
    if isinstance(sim,fpy.sim_data)==False:
        sim = fpy.sim_data(file_name,time=time)

    # Read parameters from attributes in HDF5 file
    if itor==None or r0==None:
        f = h5py.File(file_name,'r')
        if itor == None:
            itor = f.attrs["itor"]
        if r0 == None:
            r0 = f.attrs["rzero"]
    
    if(r0 == 0):
        print('Error: r0 = 0. Using r0 = 1')
        r0 = 1.
    
    
    if(itor == 0):
        period = 2.*math.pi*r0
    else:
        period = 2.*math.pi

    if points==None:
        if(x==None or z==None): #ToDo: x and y are currently not input parameters, either remove or add as input parameters
            points=200
        else:
            points=sqrt(len(x)*len(z))
    if fbins==None or fbins==0:
        fbins=points
    if tbins==None or tbins==0:
        tbins=points

    geo = False
    
    if fcoords.lower()=='pest':
        print('Creating flux coordinates using PEST angle')
        fast = False
    elif fcoords.lower()=='boozer':
        print('Creating flux coordinates using BOOZER angle')
        fast = False
    elif fcoords.lower()=='hamada':
        print('Creating flux coordinates using HAMADA angle')
        fast = False
    elif fcoords.lower()=='canonical':
        print('Creating flux coordinates using CANONICAL (equal arc length) angle')
        fast = False
    else:
        print('Creating flux coordinates using GEOMETRIC angle')
        geo = True
        fast=False
    
    print('Fast mode: '+str(fast))
    print('Using FC resolution '+str(tbins)+' , '+str(fbins))

    #ToDo: Reconsider mesh
    elms = sim.get_mesh(time=sim.timeslice)
    mp = elms.elements
    R_linspace_grid = mp[:,4]
    R_linspace = R_linspace_grid


    # Read fields
    print('Reading fields...', end=' ', flush=True)
    A = sim.get_field('A',sim.timeslice)
    gradA = sim.get_field('gradA',sim.timeslice)
    B = sim.get_field('B',sim.timeslice)

    # Get flux on axis and lcfs and magnetic axis position at the given time slice
    timeslice = sim.timeslice
    ntimestep = readParameter('ntimestep',fname='time_'+str(timeslice).zfill(3)+'.h5')
    R_magax = sim.get_time_traces('xmag').values[ntimestep]
    Z_magax = sim.get_time_traces('zmag').values[ntimestep]
    axis = [R_magax, Z_magax]
    psi_0 = R_magax*A.evaluate((R_magax,phit,Z_magax))[1]
    psi_s = sim.get_time_traces('psi_lcfs').values[ntimestep]
    print('[Done]')
    print('Magnetic axis = '+str(axis))
    print('psi_0, psi_s = '+str(psi_0)+' , '+str(psi_s))

    m = tbins # number of poloidal points
    n = fbins # number of radial points
    
    dpsi = psi_s - psi_0
    

    def field_at_point(field,R,phi,Z):
        if field == 'psi0':
            Aphi = A.evaluate((R,phi,Z))[1]
            val = R*Aphi
        elif field == 'psi0_r':
            Aphi = A.evaluate((R,phi,Z))[1]
            dAphidR = gradA.evaluate((R,phi,Z))[1]
            val = Aphi+R*dAphidR
        elif field == 'psi0_z':
            dAphidZ = gradA.evaluate((R,phi,Z))[7]
            val = R*dAphidZ
        elif field == 'psin0':
            psi0 = field_at_point('psi0',R,phi,Z)
            val = (psi0-psi_0)/dpsi
        elif field == 'psin0_r':
            psi0_r = field_at_point('psi0_r',R,phi,Z)
            val = psi0_r/dpsi
        elif field == 'psin0_z':
            psi0_z = field_at_point('psi0_z',R,phi,Z)
            val = psi0_z/dpsi
        elif field == 'gradpsi':
            psi0_r = field_at_point('psi0_r',R,phi,Z)
            psi0_z = field_at_point('psi0_z',R,phi,Z)
            val = np.sqrt(psi0_r**2 + psi0_z**2)
        elif field == 'i0':
            Bphi = B.evaluate((R,phi,Z))[1]
            val = R*Bphi
        elif field == 'b2r2':
            i0 = field_at_point('i0',R,phi,Z)
            psi0_r = field_at_point('psi0_r',R,phi,Z)
            psi0_z = field_at_point('psi0_z',R,phi,Z)
            val = i0**2 + psi0_r**2 + psi0_z**2
        return val
    
    
    if psin_range==None:
        psi = (dpsi)*(np.arange(float(n))+0.5)/n + psi_0
    else:
        p0 = psi_0 + (dpsi)*psin_range[0]
        p1 = psi_0 + (dpsi)*psin_range[1]
        psi = (p1 - p0)*(np.arange(float(n))+0.5)/n + p0
    
    psin = (psi - psi_0)/(dpsi)
    dpsi_dpsin = dpsi
    print('dpsi_dpsin = '+str(dpsi_dpsin))

    theta = 2.*math.pi*np.arange(float(m))/m
    rpath = np.zeros((m,n))
    zpath = np.zeros((m,n))
    jac = np.zeros((m,n))
    omega = np.zeros((m,n))
    q = np.zeros(n)
    dV = np.zeros(n)
    V = np.zeros(n)
    phi = np.zeros(n)
    area = np.zeros(n)
    current = np.zeros(n)

    tol_psin = 6e-5
    tol_r = 1e-4*(np.amax(R_linspace) - np.amin(R_linspace))
    maxits = 20
    rho_old = 0.
    psin_old = 0.

    tot = int(0) #long datatype is deprecated in Python3, but closely resembled by int
    
    print('Finding points on magnetic surfaces...', end=' ', flush=True)
    # find points on magnetic surfaces
    for i in range(m):
        # minus sign is because theta increases clockwise about axis
        co = np.cos(-theta[i])
        sn = np.sin(-theta[i])
        dpsin_drho = 0.
        
        for j in range(n):
            # do newton iterations to find (R,Z) at (psin, theta)
            converged = False
            if(j == 0 or dpsin_drho == 0.):
                rho = 0.01
            else:
                rho = rho + (psin[j]-psin[j-1])/dpsin_drho
            
            for k in range(maxits):
                rpath[i,j] = rho*co + axis[0]
                zpath[i,j] = rho*sn + axis[1]
                psin_x = field_at_point('psin0',rpath[i,j],phit,zpath[i,j])

                if abs(psin_x - psin[j]) < tol_psin:
                    converged = True
                    kk=k
                    #print('converged in iteration'+str(i)+','+str(j)+','+str(k))
                    break

                if fast==True:
                    if((k == 0 and j == 0) or abs(psin_x - psin_old) < tol_psin or abs(rho - rho_old) < tol_r):
                        psin_r = field_at_point('psin0_r',rpath[i,j],phit,zpath[i,j])
                        psin_z = field_at_point('psin0_z',rpath[i,j],phit,zpath[i,j])
                        dpsin_drho = psin_r*co + psin_z*sn
                    else:
                        dpsin_drho = (psin_x - psin_old) / (rho - rho_old)

                    rho_old = rho
                    psin_old = psin_x

                else:
                    psin_r = field_at_point('psin0_r',rpath[i,j],phit,zpath[i,j])
                    psin_z = field_at_point('psin0_z',rpath[i,j],phit,zpath[i,j])
                    dpsin_drho = psin_r*co + psin_z*sn

                drho = (psin[j] - psin_x)/dpsin_drho
                if(rho + drho < 0.):
                    rho = rho / 2.
                else:
                    rho = rho + drho
                #print('rho = '+str(rho))
            tot = tot + kk
            if converged == False:
                print('Error at (psi_n, theta) = '+str(psin[j])+' , '+str(theta[i]))
                print('Did not converge after '+str(maxits)+' iterations. , '+str(rho)+' , '+str(rpath[i,j])+' , '+str(zpath[i,j]))
                rho = 0.01
                return
    print('[Done]')
    print('Total Newton iterations '+str(tot))
    print('Average Newton iterations '+str(float(tot)/float(int(m)*int(n))))

    #ToDo: Remove Plot
    #plt.figure()
    #pathshape = rpath.shape
    #for i in range(pathshape[1]):
    #    plt.plot(rpath[:,i],zpath[:,i])
    #plt.axis('equal')
    
#  find pest angle
    if(psi_s < psi_0):
        print('Ip > 0')
        print('grad(psi) is inward')
        fac = -1.
    else:
        print('Ip < 0')
        print('grad(psi) is outward')
        fac = 1.

    theta_sfl = np.zeros((m,n))
    f = np.zeros(n)
    fjr2 = np.zeros(m)
    dthetadl = np.zeros(m)
    bp = np.zeros(m)

    rp = np.copy(rpath)
    if(itor == 0):
        rp.fill(1.)
    for j in range(n):
        #rx = [rpath[m-1,j],rpath[:,j],rpath[0,j]]
        rx = np.insert(rpath[:,j],0,rpath[m-1,j])
        rx = np.append(rx,rpath[0,j])
        #zx = [zpath[m-1,j],zpath[:,j],zpath[0,j]]
        zx = np.insert(zpath[:,j],0,zpath[m-1,j])
        zx = np.append(zx,zpath[0,j])
        drx = fpyl.deriv(rx)
        dzx = fpyl.deriv(zx)
        dr = drx[1:m+1]
        dz = dzx[1:m+1]
        
        dl = np.sqrt(dr**2 + dz**2)
        br = np.zeros_like(rpath[:,j])
        bz = np.zeros_like(rpath[:,j])
        for kk in range(len(rpath[:,j])):
            br[kk] = -field_at_point('psi0_z',rpath[kk,j],phit,zpath[kk,j])/rp[kk,j]
            bz[kk] = field_at_point('psi0_r',rpath[kk,j],phit,zpath[kk,j])/rp[kk,j]
        
        
        
        if fcoords.lower()=='canonical':
            for i in range(m):
                bp[i] = field_at_point('gradpsi',rpath[i,j],phit,zpath[i,j])/rp[i,j]
                ix = field_at_point('i0',rpath[i,j],phit,zpath[i,j])
                fjr2[i] = -fac*ix/(rp[i,j]**2*bp[i])
        
        
        for i in range(m):
            if not fcoords.lower()=='canonical':
                bp[i] = field_at_point('gradpsi',rpath[i,j],phit,zpath[i,j])/rp[i,j]

            if fast == False:
                if not fcoords.lower()=='canonical':
                    ix = field_at_point('i0',rpath[i,j],phit,zpath[i,j])
                    fjr2[i] = -fac*ix/(rp[i,j]**2*bp[i])
                # dtheta/dl ~ 1./(Bp*Jac)
                if fcoords.lower()=='pest':
                    dthetadl[i] = -fac*ix/(rp[i,j]**2*bp[i])
                elif fcoords.lower()=='boozer':
                    dthetadl[i] = -fac*(bp[i]**2 + (ix/rp[i,j])**2) / bp[i]
                elif fcoords.lower()=='hamada':
                    dthetadl[i] = -fac*1./bp[i]
                elif fcoords.lower()=='canonical':
                    dthetadl[i] = -fac*ix/(np.sum(fjr2*dl)/period)/(rp[i,j]**2*bp[i])
            
                if(i == 0):
                    theta_sfl[i,j] = 0.
                else:
                    theta_sfl[i,j] = theta_sfl[i-1,j] + dl[i]*dthetadl[i]/2. + dl[i-1]*dthetadl[i-1]/2.
        
        current[j] = np.sum(br*dr + bz*dz)
        area[j] = period*np.sum(dl*rp[:,j])
        dV[j] = period*np.sum(fac*dpsi_dpsin*dl/bp)
        if(j == 0):
            V[j] = dV[j]*psin[j]/2.
        else:
            V[j] = V[j-1] + (dV[j]+dV[j-1])*(psin[j]-psin[j-1])/2.

        # calculate q
        if fast == False:
            # normalize theta to 2 pi
            if(geo == False):
                f[j] = np.sum(dthetadl*dl)/(2.*math.pi)
                theta_sfl[:,j] = theta_sfl[:,j]/f[j]
            q[j] = np.sum(fjr2*dl)/period
            
            # calculate toroidal flux
            if(j == 0):
                phi[j] = -period*q[j]*(psi[j]-psi_0)
            else:
                phi[j] = phi[j-1] - period*(q[j]+q[j-1])*(psi[j]-psi[j-1])/2.
        
        # sanity checks
        if(V[j] < 0):
            print('ERROR, volume is negative')
            print(j, V[j], dV[j], psi[j]-psi_0, bp[j])
            return 0
        if fast == False:
            if(theta_sfl[0,j] > theta_sfl[m-1,j]):
                print('ERROR, theta_sfl is clockwise')
                return 0
            if(fac*q[j]*ix > 0.):
                print('ERROR, q has wrong sign')
                return 0
            if(phi[j]*ix < 0.):
                print('ERROR, phi has wrong sign')
                return 0
    
    if(geo == False):
        for j in range(n):
            # interpolate fields to be evenly spaced in PEST angle
            #newm = interpol(findgen(m),theta_sfl[*,j],theta) #IDL
            ##temp = sciint.interp1d(theta_sfl[:,j],np.arange(float(m)))
            ##newm = temp(theta) #CONTINUE HERE
            #if j==0:
            #    plt.figure()
            #    plt.plot(rpath[:,j])
            #rpath[:,j] = interpolate(rpath[*,j],newm)
            #zpath[:,j] = interpolate(zpath[*,j],newm)
            tempr = interp1d(theta_sfl[:,j],rpath[:,j])
            rpath[:,j] = tempr(theta)
            tempz = interp1d(theta_sfl[:,j],zpath[:,j])
            zpath[:,j] = tempz(theta)
            if(itor == 1):
                rp[:,j] = rpath[:,j]
            
            # use analytic expression for Jacobian
            if fcoords.lower()=='pest':
                for i in range(m):
                    jac[i,j] = -dpsi_dpsin*rp[i,j]**2*q[j] / field_at_point('i0',rpath[i,j],phit,zpath[i,j])
            elif fcoords.lower()=='boozer':
                for i in range(m):
                    jac[i,j] = -dpsi_dpsin*f[j]*rp[i,j]**2 / field_at_point('b2r2',rpath[i,j],phit,zpath[i,j])
            elif fcoords.lower()=='hamada':
                jac[:,j] = -dpsi_dpsin*f[j]
            elif fcoords.lower()=='canonical':
                for i in range(m):
                    jac[i,j] = -dpsi_dpsin*f[j]/field_at_point('gradpsi',rpath[i,j],phit,zpath[i,j])

    #  Calculate Jacobian if requested (or if analytic expression isn't available)
    if (geo == True or njac==True):
        print('Calculating Jacobian numerically')
        # calculate jacobian
        dr_dpsi = np.copy(rpath)
        dz_dpsi = np.copy(zpath)
        dr_dtheta = np.copy(rpath)
        dz_dtheta = np.copy(zpath)
        for i in range(m):
            dr_dpsi[i,:] = fpyl.deriv(rpath[i,:],psin)
            dz_dpsi[i,:] = fpyl.deriv(zpath[i,:],psin)
        for j in range(n):
            dr_dtheta[:,j] = fpyl.deriv(rpath[:,j],theta)
            dz_dtheta[:,j] = fpyl.deriv(zpath[:,j],theta)

        jac_test = -(dr_dpsi*dz_dtheta - dr_dtheta*dz_dpsi)
        if (itor == 1):
            jac_test = jac_test*rpath

        if (np.mean(jac_test) < 0):
            print('ERROR: numerical jacobian is negative!')
            return 0

        jac = jac_test

    if (np.mean(jac) < 0):
        print('ERROR: jacobian is negative!')
        return 0


    # calculate deviation of toroidal angle from geometric toroidal angle: omega := zeta - phi
    if (geo == False and fast==False):
        for j in range(n):
            for i in range(m):
                jac_pest = -dpsi_dpsin*rp[i,j]**2*q[j] / field_at_point('i0',rpath[i,j],phit,zpath[i,j])
                if i == 0:
                    omega[i,j] = 0.
                else:
                    dtheta = theta[i] - theta[i-1]
                    omega[i,j] = omega[i-1,j] + dtheta*((1. - jac[i,j]/jac_pest)/2. +(1. - jac[i-1,j]/jac_pest_old)/2.)
                jac_pest_old = jac_pest
            omega[:,j] = omega[:,j]*q[j]

    # Create flux coordinates object
    # psi = poloidal flux
    # phi = toroidal flux
    fc = fpy.flux_coordinates(m,n,rpath,zpath,axis,omega,psi,psin,period,theta,jac,q,area,dV,fcoords,V,phi,itor,r0,current,dpsi_dpsin)
    sim.fc = fc

    if makeplot==True:
        fig = plt.figure(constrained_layout=True,figsize=(15,7),num=fignum)
        spec2 = gridspec.GridSpec(ncols=4, nrows=3, figure=fig)
        f_ax1 = fig.add_subplot(spec2[:, 0])
        f_ax2 = fig.add_subplot(spec2[:2, 1])
        f_ax3 = fig.add_subplot(spec2[0, 2])
        f_ax4 = fig.add_subplot(spec2[0, 3])
        f_ax5 = fig.add_subplot(spec2[1, 2])
        f_ax6 = fig.add_subplot(spec2[1, 3])
        f_ax7 = fig.add_subplot(spec2[2, 2])
        f_ax8 = fig.add_subplot(spec2[2, 3])
        
        pathshape = rpath.shape
        fs = list(range(0,pathshape[1],10))
        if fs[-1] != pathshape[1]-1:
            fs.append(pathshape[1]-1)
        for i in fs:
            rp_plot = np.append(rpath[:,i],rpath[0,i])
            zp_plot = np.append(zpath[:,i],zpath[0,i])
            f_ax1.plot(rp_plot,zp_plot,c='C0')
        for i in range(0,pathshape[0],10):
            f_ax1.plot(rpath[i,:],zpath[i,:],c='C1')
        
        f_ax1.set_xlim([fpyl.get_axlim(np.amin(rpath),'min',0.1),fpyl.get_axlim(np.amax(rpath),'max',0.1)])
        f_ax1.set_ylim([fpyl.get_axlim(np.amin(zpath),'min',0.1),fpyl.get_axlim(np.amax(zpath),'max',0.1)])
        f_ax1.set_aspect('equal')
        f_ax1.set_xlabel(r'R')
        f_ax1.set_ylabel(r'Z')
        f_ax1.grid()
        f_ax1.set_title('Flux surfaces')
        
        cont = f_ax2.contourf(theta,psin,jac.transpose(),200,cmap='jet')
        cbar = fig.colorbar(cont,ax=f_ax2,orientation="horizontal", pad=0.02)
        cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(),rotation=270)
        f_ax2.set_aspect(2*math.pi)
        f_ax2.set_xlabel(r'$\theta$')
        f_ax2.set_ylabel(r'$\psi_N$')
        f_ax2.set_title('Jacobian')
        
        f_ax3.plot(psin,q)
        f_ax3.set_xlabel(r'$\psi_N$')
        f_ax3.set_ylabel(r'$q$')
        f_ax3.grid()
        f_ax4.plot(psin,V)
        f_ax4.set_xlabel(r'$\psi_N$')
        f_ax4.set_ylabel(r'$V$')
        f_ax4.grid()
        f_ax5.plot(psin,area)
        f_ax5.set_xlabel(r'$\psi_N$')
        f_ax5.set_ylabel(r'Area')
        f_ax5.grid()
        f_ax6.plot(psin,current)
        f_ax6.set_xlabel(r'$\psi_N$')
        f_ax6.set_ylabel(r'current')
        f_ax6.grid()
        f_ax7.plot(theta,omega[:,20])
        f_ax7.set_xlabel(r'$\theta$')
        f_ax7.set_ylabel(r'omega')
        f_ax7.grid()
        f_ax8.plot(psin,phi)
        f_ax8.set_xlabel(r'$\psi_N$')
        f_ax8.set_ylabel(r'$\phi$')
        f_ax8.grid()

    return sim
