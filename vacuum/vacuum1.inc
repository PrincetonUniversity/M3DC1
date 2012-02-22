ccc      cliche pestlst
c
c.....    nfs0 = 401, ntsin0=256 for Manickam's latest 3/93, TPX cases. 
c                                and on 04-29-04.        Will
c                                use this consistently.
c.....    nfs0 = 401, ntsin0=128 for Manickam's latest 11/5, FSX cases.
c.....    nfs0 = 513, ntsin0=256 for Manickam-Meynard NTSX cases.5/10/95
c.....    nfs0 = 513, ntsin0=256 for Manickam-Meynard NTSX cases.5/10/95
c.....    nsf0 = 201, ntsin0=402, for PBX cases 1/93.
c.....    nfs0 = 201, ntsin0=128 for Manickam's PBX cases 1/93 & 5/3/93 using
c..........                      PEST1 (GO10YB)
c.....    nsf0 = 201, for pest1, 193 otherwise.
c.....    nsf0 = 201, ntsin0 = 128  for Ming Chu's mapping
c.....    nsf0 = 401, ntsin0 = 128  for the Resistive Wall Stuff.
c.....    nsf0 = 101, ntsin0 = 64 for manickam case, 201,128 otherwise
c.....    nsf0 = 1201 ntsin0 = 512 for STX
c.....    nfm = 41 for pest1, 31 otherwise.?
c.....    nfm, ntot increased to 51 for TPX largfe cases 9/21/93
c.....     ntsin, ntsin0 for input from mapping.
c.....     nths, nths0 for working grid.
c
c$$$      parameter ( nsf0=401, nths0=1001, nfm=401, mtot=401, !  08/14/2009
c$$$      parameter ( nsf0=401, nths0=1001, nfm=401, mtot=401,
c$$$      parameter ( nsf0=401, nths0=2060, nfm=401, mtot=401,  !! temporarily
      parameter ( nsf0=401, nths0=1040, nfm=401, mtot=401,   ! 10/27/2010
     .            ntsin0=256, ntsin=ntsin0+5,
     .            nsf=nsf0+1, nfe=1+nsf/2, nths=nths0+5,
     .            nths2=2*nths, nfm2=2*nfm, nwrkl=nths2**2+3*nths2,
     .            nthsq=nths*nths, nthsq4=4*nthsq, nfnt4=4*nths*nfm,
     .            nfm4=4*nfm, nfmsq=nfm*nfm, nfm21=2*nfm + 1 )
      parameter ( nwm = 201, nwm2 = 2*nwm, nwm4 = 4*nwm,
     $     nwmsq = nwm*nwm )
      PARAMETER ( npfg11=1, npfg12=(nfm/2)+1, npfg13=npfg11+nfm,
     $     npfg14=npfg12+nfm, npfg2x=nfm*(nfm-2*(nfm/2))+1 )
      PARAMETER ( nptg11=1, nptg12=(nths/2)+1, nptg13=nptg11+nths,
     $     nptg14=nptg12+nths, nptg2x=nths*(nths-2*(nths/2))+1 )
c
      parameter ( maxct = 24 )
      real n, no2pi, no2pi2
      logical lsymz, check1, check2, lanal, lpest1, lpless, wall,
     .        lnova, checks
      integer outpest, outmod, outmap1
c
      common / control / lsymz, check1, check2, checks, lpestasym,
     $     lanal, lpest1, lnova, ipshp, ladj, ldcon, lgato, lrnim3d,
     $     literr, lm3dc1,
     $     lrgato, lrnova, lreig7, lrwall7, lgaus, ltnrm, ltfil, lkplt,
     $     lspark, lfele, lfour,
     $     lkdis,ldisc,ldqdtw,ld3dflx,
     $     lhighn,lreshel,lfbcoil,lfbcoila,lfbcoilb,
     $     lfarw,leqarcw,lsymwv,
     $     lplot1, lplot2, lplot3, lplot4, lplot5, irepeat, nrepeat
c
      common / consts / zero, half, one, two, three, four, five, seven,
     .                  pye, twopi, twopi2, bit, no2pi, no2pi2, amu0,
     $     epszer
c
      integer datime(6)
      character*(maxct) time0,date0, timem,datem, timev,datev,
     $       date_array, mp0,mp1, jobid, filein2
      character*60 seps
      common / datim0 / seps, time0,date0, timem,datem, timev,datev,
     $     date_array
      common / datim1 / datime
      common / map111 / mp0, mp1, jobid, filein2, lenif
c
      common / chanel / inpest, outpest, iomode, inmode, outmod,
     .                  outmap1, intty, iotty, nout(3), nout0(3),
     $       ioshel, iovac, iodsk, iosvd, idsk, lzio
c
      common / gridin / dthin, mthin1, mthin2, mdiv, mp
      common / gridm / dth, mth, mth1, mth2, dmfel, mfel, mfel1,
     $     mthout, mthout1, mthout2,
     $     dm0, ndfel, dx0, dx1
c
      common / geom / r2, r4, r6, rgato
c
      common / index0 / lmax(nfe), lmin(nfe), lnint,lxint, lnsav,lxsav
      common / displ / xiinr(nfm),  xiini(nfm), xirc(nfm), xirs(nfm),
     $     xiic(nfm), xiis(nfm), xigr(nths), xigi(nths)
      common / index1 /  n, wall, m, nmin, nmax
      common / stable / eigval(100), nunst,neigvl,neigvc, iloop
c
      common / equil / fa1, ga1, qa1, s_scale, psipr
c
      common / equ5 / ntitle(20),dat,nx,nz,nosurf,mthin,lpless,alx,alz,
     .                xzero, xma, r, p0, gp0, psimin, psilim, psipls,
     .                betag, betap, upsiln, nj,mj,lj, zma, betav,betai,
     $     ctroy, z47,z48,z49,eqdt,eqdr,eqpi,x000,z54,z55,z56,z57,
     $     nzd1,nzd2
c
      dimension grpssq(nths), xsq(nths), gpsdth(nths), xsqdth(nths),
     $     xjacob(nths),enqth(nths),qth0(nths),qth0pp(nths),
     $     delta(nths),xjdtxj(nths),
     $     xsdtxs(nths),gpdtgp(nths)
      common / surf / grpssq, xsq, gpsdth, xsqdth, xjacob,
     $     enqth,qth0,qth0pp,delta,
     .          xjdtxj, xsdtxs, gpdtgp, slngth(nths),ellp(nths),
     $     sqgrps(nths)
c
c      **************************************************************
c
