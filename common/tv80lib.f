c**************************************************************************
c
c  contur.F - contouring routines
c
c  contents	acontr
c		rcontr
c		contur
c		xcontr
c		lcurvs
c		q7plot
c		q8plot
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

c**************************************************************************
c
c  acontr
c
c  calls        xcontr
c
c**************************************************************************

      subroutine acontr (jjj1, bc, jjj2, ba, max, imin, imax, istop
     * , jmin, jmax, jstop)
c=======================================================================
c 12/28/2001 ler pppl - change last two arguments in call to xcontr to
c                        double precision
c=======================================================================
c
cc variable declarations:
c
      dimension bc(*), ba(*)
      real*8 bc, ba
c     dimension ba(1), bc(1)

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c If called by user set routine top routine name

      if (tgname .eq. 'NONE') tgname = 'ACONTR'



c
cc program statements:
c
c% The following call to xcontur used to only have 12 arguments.  The
c% last 2 arguments were set to -1 to make the call consistent with
c% the definition.  The last 2 arguments should not be used by xcontr.
c
      i = 1
      call xcontr (i, jjj1, bc(1), jjj2, ba, max, imin, imax, istop
     * , jmin, jmax, jstop, -1.0d0, -1.0d0)
c
9999  if (tgname .eq. 'ACONTR') tgname = 'NONE'
      return
      END

c**************************************************************************
c
c  rcontr
c
c  calls        xcontr
c
c**************************************************************************

      subroutine rcontr (jjj1, bc, jjj2, ba, max, axarr, imin, imax
     * , istop, ayarr, jmin, jmax, jstop)
c
cc variable declarations:
c
c     dimension ba(1), bc(1), axarr(1), ayarr(1)
      dimension bc(*), ba(*), axarr(*), ayarr(*)
      real*8 bc, ba, axarr, ayarr
      save

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c If called by user set routine top routine name

      if (tgname .eq. 'NONE') tgname = 'RCONTR'



c
cc program statements:
c
      i = 2
      call xcontr (i, jjj1, bc(1), jjj2, ba, max, imin, imax
     * , istop, jmin, jmax, jstop, axarr, ayarr)
c
9999  if (tgname .eq. 'RCONTR') tgname = 'NONE'
      return
      END

c**************************************************************************
c
c  contur
c
c  calls        gqcntn
c               gqnt
c               setcrt
c               q7plot
c
c**************************************************************************

      subroutine contur (jjj1, bc, jjj2, ba, max, axarr, imin, imax
     * , istop, ayarr, jmin, jmax, jstop)
c======================================================================
c 12/28/2001 ler pppl - change arguments to "call setcrt" to double
c                        precision
c=======================================================================
c
cc purpose: draw a contour plot
c           (it is assumed that ba(i,j) = fcn(axarr(i,j),ayarr(i,j)))
c
cc author: f. n. fritsch
c
      common /tvgxx1/ ksent, kclip, krastr, ktyp, cxfact, cxcons
     * , cyfact, cycons, cxmi, cxma, cymi, cyma, klx, kly, kpol
     * , knamsb, knerr, crasmx, karrsz, cxarr(256)
     * , cyarr(256), cxlast, cylast, kspace, kcfx, kcfy
c
      common /q7quad/ kp, ctest, cxc(5), cyc(5), czc(5)
     * , int, k01, k02, c01, cdf, czmin, czmax
c
c      dimension ba(2), bc(2), axarr(2), ayarr(2)
      dimension ba(*), axarr(*), ayarr(*), bc(*)
      real*8 ba, axarr, ayarr, bc
        dimension wd(4),vp(4)
c
      equivalence (ima, nmax)
      equivalence (i, i1m1)

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c If called by user set routine top routine name

      if (tgname .eq. 'NONE') tgname = 'CONTUR'



c
c set min and max NDC window
c
      call gqcntn (ie,nt)
      call gqnt (nt,ie,wd,vp)
      cxmi = vp(1)
      cxma = vp(2)
      cymi = vp(3)
      cyma = vp(4)
c
cc program statements:
c
      call setcrt (0.0d0, 0.0d0)
      kp = 0
      k02 = jjj2
      k01 = jjj1
      if (k01) 1, 50, 3
c
c get contours for negative k01 option...
    1 k01 = -k01
      bc(k01) = bc(2)
      i1m1 = k01-1
      cdf = (bc(k01)-bc(1))/float(i1m1)
      do 2 j = 2, i1m1
         bc(j) = bc(j-1)+cdf
    2 continue
c
c set up addresses for q7plot..
    3 continue
c
c begin logic for sweeping through grid...
    4 continue
      istep = istop
      ima = imax-istep
      jstep = jstop
      jma = jmax-jstep
      nstep = jstep*max
      nmin = (jmin-1)*max
      nmax = nmin+ima
      nmin = nmin+imin
c
      do 40 itemz = nmin, nmax, istep
c
c element (i,jmim)
      cxc(1) = axarr(itemz)
      cyc(1) = ayarr(itemz)
      czc(1) = ba(itemz)
c
c element (i+istep,jmin)
      cxc(4) = axarr(itemz+istep)
      cyc(4) = ayarr(itemz+istep)
      czc(4) = ba(itemz+istep)
      if (czc(4) .ge. czc(1)) go to 22
      czmax = czc(1)
      czmin = czc(4)
      go to 23
c
   22 czmax = czc(4)
      czmin = czc(1)
   23 continue
      karrsz = itemz
c
      do 35 j = jmin, jma, jstep
      karrsz = karrsz+nstep
c
c element (i,j+jstep)
      cxc(2) = axarr(karrsz)
      cyc(2) = ayarr(karrsz)
      czc(2) = ba(karrsz)
c
c element (i+istep,j+jstep)
      cxc(3) = axarr(karrsz+istep)
      cyc(3) = ayarr(karrsz+istep)
      czc(3) = ba(karrsz+istep)
      if (czc(3) .ge. czc(2)) go to 24
      czmax2 = czc(2)
      czmin2 = czc(3)
      go to 25
c
   24 czmax2 = czc(3)
      czmin2 = czc(2)
c
   25 if (czmax .ge. czmax2) go to 26
      czmax = czmax2
c
   26 if (czmin2 .ge. czmin) go to 27
      czmin = czmin2
c
   27 call q7plot (bc(1))
      if (kp .eq. 1) go to 9999
      cxc(1) = cxc(2)
      cyc(1) = cyc(2)
      czc(1) = czc(2)
      cxc(4) = cxc(3)
      cyc(4) = cyc(3)
      czc(4) = czc(3)
      czmax = czmax2
      czmin = czmin2
   35 continue
   40 continue
c
c     plotting completed...
   30 continue
      go to 9999
c
c     set up for k01 = 0 option...
   50 c01 = bc(1)
      cdf = bc(2)
      if (cdf .gt. 0) go to 4
C
c     error exit if  cdf .le. 0 ...
      print *,'error 3: cdf .le. 0'
c
9999  if (tgname .eq. 'CONTUR') tgname = 'NONE'
      return

      END

c**************************************************************************
c
c  xcontr
c
c  calls        gqcntn
c               gqnt
c               setcrt
c               q7plot
c
c  description  This routine was designed to be called by the contour routines
c               only.
c
c**************************************************************************

      subroutine xcontr (ii, jjj1, bc, jjj2, ba, max, imin, imax
     * , istop, jmin, jmax, jstop, axarr, ayarr)
c======================================================================
c 12/28/2001 ler pppl - change arguments to "call setcrt" to double
c                        precision
c=======================================================================
c
cc author: f. n. fritsch
c
cc notes:
c
c 1) it is assumed that a(i,j) = fcn(axarr(i),ayarr(j)).
c
cc variable declarations:
c
      common /tvgxx1/ ksent, kclip, krastr, ktyp, cxfact, cxcons
     * , cyfact, cycons, cxmi, cxma, cymi, cyma, klx, kly, kpol
     * , knamsb, knerr, crasmx, karrsz, cxarr(256)
     * , cyarr(256), cxlast, cylast, kspace, kcfx, kcfy
c
c     dimension ba(2), bc(2), axarr(2), ayarr(2)
      dimension ba(*), bc(*), axarr(*), ayarr(*)
      real*8 ba, bc, axarr, ayarr
       dimension wd(4),vp(4)
c
      common /q7quad/ kp, ctest, cxc(5), cyc(5), czc(5)
     * , int, k01, k02, c01, cdf, czmin, czmax
c
      equivalence (i,i1m1)
      save

c
c set min and max NDC window
c
        call gqcntn (ie,nt)
        call gqnt (nt,ie,wd,vp)
        cxmi = vp(1)
        cxma = vp(2)
        cymi = vp(3)
        cyma = vp(4)
c
cc program statements:
c
      call setcrt (0.0d0, 0.0d0)
      kp = 0
      k02 = jjj2
      k01 = jjj1
      if (k01) 1, 50, 3
c
c  get contours for negative k01 option...
    1 k01 = -k01
      bc(k01) = bc(2)
      i1m1 = k01-1
      cdf = (bc(k01)-bc(1))/float(i1m1)
      do 2 j = 2, i1m1
         bc(j) = bc(j-1)+cdf
    2 continue
c
c  set up addresses for q7plot..
    3 continue
c
c  begin logic for sweeping through grid...
    4 continue
      istep = istop
      ima = imax-istep
      jstep = jstop
      jma = jmax-jstep
      nstep = jstep*max
      nmin = (jmin-1)*max
      do 40  i = imin, ima, istep
      if (ii .eq. 1) go to 10
      cxc(1) = axarr(i)
      cxc(2) = axarr(i)
      cxc(3) = axarr(i+istep)
      cxc(4) = axarr(i+istep)
      go to 15
c
   10 cxc(1) = i
      cxc(2) = i
      cxc(3) = i+istep
      cxc(4) = i+istep
c
   15 continue
c
c element (i,jmin)
      czc(1) = ba(nmin+i)
c
c element (i+istep,jmin)
      czc(4) = ba(nmin+i+istep)
      if (czc(4) .ge. czc(1)) go to 22
      czmax = czc(1)
      czmin = czc(4)
      go to 23
c
   22 czmax = czc(4)
      czmin = czc(1)
   23 continue
      itemz = nmin+i
      do 35 j = jmin, jma, jstep
      if (ii .eq. 1) go to 2350
      cyc(1) = ayarr(j)
      cyc(2) = ayarr(j+jstep)
      cyc(3) = ayarr(j+jstep)
      cyc(4) = ayarr(j)
      go to 2375
c
 2350 cyc(1) = j
      cyc(4) = j
      cyc(2) = j+jstep
      cyc(3) = j+jstep
 2375 continue
      itemz = itemz+nstep
c
c element (i,j+jstep)
      czc(2) = ba(itemz)
c
c element (i+istep,j+jstep)
      czc(3) = ba(itemz+istep)
      if (czc(3) .ge. czc(2)) go to 24
      czmax2 = czc(2)
      czmin2 = czc(3)
      go to 25
c
   24 czmax2 = czc(3)
      czmin2 = czc(2)
c
   25 if (czmax .ge. czmax2) go to 26
      czmax = czmax2
c
   26 if (czmin2 .ge. czmin) go to 27
      czmin = czmin2
c
   27 call q7plot (bc(1))
      if (kp .eq. 1) go to 9999
c
   34 czc(1) = czc(2)
      czc(4) = czc(3)
      czmax = czmax2
      czmin = czmin2
   35 continue
   40 continue
c
c     plotting completed...
   30 continue
      go to 9999
c
c     set up for k01 = 0 option...
   50 c01 = bc(1)
      cdf = bc(2)
      if (cdf .gt. 0) go to 4
C
c     error exit if  cdf .le. 0 ...
      print *,'error 3: cdf .le. 0'
c
9999  return
      END

c**************************************************************************
c
c  lcurvs
c
c  calls        gqcntn
c               gqnt
c               setcrt
c               q7plot
c
c**************************************************************************

      subroutine lcurvs (jjj1, bc, jjj2, fcn, nx, ny)
      real*8 fcn
      external fcn
c=======================================================================
c 12/28/2001 ler pppl - explicitly define passed function 'fcn' as 
c                        double precision
c                     - change first argument in first 'fcn' invocation 
c                        from single to double precision
c=======================================================================
c
cc notes:
c
c 1) fourth argument (fcn) should be typed generic according
c    to the original documentation in ltss chapter 304.
c    however, for tv80cray, it is assumed that it is an ascii string
c
cc variable declarations:
c
C% None in tvgxx1 are used
C
      common /tvgxx1/ ksent, kclip, krastr, ktyp, cxfact, cxcons
     * , cyfact, cycons, cxmi, cxma, cymi, cyma, klx, kly, kpol
     * , knamsb, knerr, crasmx, karrsz, cxarr(256)
     * , cyarr(256), cxlast, cylast, kspace, kcfx, kcfy
c
      common /tvgxx2/ kxarr(511), kyarr(511), kitcnt, klchar
     * , kntenv, kdashv
     * , knumch, kcoff, kcpt, kcprm(7,4), klindx, klindy
     * , kchdx, kchdy, kindyo
     * , klabel(4,7), kltv80(9), krasmx, kact
     * , kdstat(8), kname(8), klenfi(8), klevd(8)
     * , kdbuf(8), kdbufs(8), knfrm(8)
     * , kepflg(8), kmxfil(8), kminbf(8), klerad(8)
     * , kdfact(8), knumfi(8)
c
      common /tvmap1/ kia, knarg, cx1, cx2, cy1, cy2
c
      common /q7quad/ kp, ctest, cxc(5), cyc(5), czc(5)
     * , int, k01, k02, c01, cdf, czmin, czmax
c
c      dimension bc(2)
      dimension bc(*)
      real*8 bc, x, y, xx, yy, dcx1

c        dimension wd(4),vp(4)
c
      equivalence (i, i1m1)

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c If called by user set routine top routine name

      if (tgname .eq. 'NONE') tgname = 'LCURVS'



c
c get min and max NDC window coordinates and user coordinates
c
      call tggetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,it)
c
c The user coordinate should not really be in log form.  But
c if they are, tggetmap will return exponents.  Convert the
c exponents into coordinates.
c
      if (it .eq. 2 .or. it .eq. 4) then
        cy1 = 10.0**cy1
	cy2 = 10.0**cy2
      endif
      if (it .eq. 3 .or. it .eq. 4) then
        cx1 = 10.0**cx1
	cx2 = 10.0**cx2
      endif
c
cc program statements:
c
      kp = 0
      k02 = jjj2
      k01 = jjj1
      if (k01) 1, 50, 3
c
c  get contours for negative k01 option...
    1 k01 = -k01
      bc(k01) = bc(2)
      i1m1 = k01-1
      cdf = (bc(k01)-bc(1))/float(i1m1)
      do 2 j = 2, i1m1
         bc(j) = bc(j-1)+cdf
    2 continue
c
c  set up addresses for q7plot..
    3 continue
c
c  begin logic for sweeping through grid...
    4 continue
      dx = (cx2-cx1)/nx
      dy = (cy2-cy1)/ny
      y = cy1
      do 40 i = 1, ny
      cyc(1) = y
      cyc(2) = y+dy
      cyc(3) = y+dy
      cyc(4) = y
      yy = y+dy
      dcx1 = cx1
      czc(1) = fcn(dcx1,y)
      czc(2) = fcn(dcx1,yy)
      if (czc(2) .ge. czc(1)) go to 22
      czmax = czc(1)
      czmin = czc(2)
      go to 23
c
   22 continue
      czmax = czc(2)
      czmin = czc(1)
   23 continue
      x = cx1
      do 35 j = 1, nx
      cxc(1) = x
      cxc(2) = x
      cxc(4) = x+dx
      cxc(3) = cxc(4)
      xx = cxc(4)
      czc(4) = fcn(xx,y)
      czc(3) = fcn(xx,yy)
      if (czc(3) .ge. czc(4)) go to 24
      czmacx = czc(4)
      czmin2 = czc(3)
      go to 25
c
   24 continue
      czmacx = czc(3)
      czmin2 = czc(4)
   25 continue
      if (czmax .ge. czmacx) go to 26
      czmax = czmacx
   26 if (czmin2 .ge. czmin) go to 27
      czmin = czmin2
   27 call q7plot (bc(1))
      if (kp .eq. 1) go to 9999
      czc(1) = czc(4)
      czc(2) = czc(3)
      czmax = czmacx
      czmin = czmin2
      x = xx
   35 continue
      y = yy
   40 continue
c
c     plotting completed...
   30 continue
      go to 9999
c
c     set up for k01 = 0 option...
   50 c01 = bc(1)
      cdf = bc(2)
      if (cdf .gt. 0) go to 4
C
c     error exit if cdf .le. 0 ...
      print *,'error 3: cdf .le. 0'
c
9999  if (tgname .eq. 'LCURVS') tgname = 'NONE'
      return
      END

c**************************************************************************
c
c  q7plot
c
c  calls        q8plot
c
c  description  This was designed to only be called by the contour routines.
c
c**************************************************************************

      subroutine q7plot (bc)
c
cc author: f. n. fritsch
c
cc revised: dec14.1978 to conform to 7600 version
c
cc variable declarations:
c
      common /tvgxx1/ ksent, kclip, krastr, ktyp, cxfact, cxcons
     * , cyfact, cycons, cxmi, cxma, cymi, cyma, klx, kly, kpol
     * , knamsb, knerr, crasmx, karrsz, cxarr(256)
     * , cyarr(256), cxlast, cylast, kspace, kcfx, kcfy
c
      common /tvgxx2/ kxarr(511), kyarr(511), kitcnt, klchar
     * , kntenv, kdashv
     * , knumch, kcoff, kcpt, kcprm(7,4), klindx, klindy
     * , kchdx, kchdy, kindyo
     * , klabel(4,7), kltv80(9), krasmx, kact
     * , kdstat(8), kname(8), klenfi(8), klevd(8)
     * , kdbuf(8), kdbufs(8), knfrm(8)
     * , kepflg(8), kmxfil(8), kminbf(8), klerad(8)
     * , kdfact(8), knumfi(8)
c
c
      dimension bc(*)
      real*8 bc
c     dimension bc(2)
c
      common /q7quad/ kp, ctest, cxc(5), cyc(5), czc(5)
     * , int, k01, k02, c01, cdf, czmin, czmax
c
c next line revised dec14.78 to conform to 7600 version
      save toolrg
      data toolrg /100.0/
c
cc program statements:
c
      cxc(5) = cxc(1)
      cyc(5) = cyc(1)
      czc(5) = czc(1)
c
      if (k01 .ne. 0) go to 7
c
c  plot all possible contours at spacing fdf from c01.
      fdf = cdf
c
c next 3 lines revised 14-dec-78 to conform to 7600 version
      ncont = (czmax-czmin)/fdf
      if ((czmax-czmin) .gt. (toolrg*fdf)) go to 50
c     if (toolrg .eq. 100.0) go to 200
c
      if ((czmax-czmin)/fdf-2.0) 200, 200, 1000
c
 1000 continue
c
c  increase fdf if 'too many' contours pass through box...  (fdfdet)
      do 201 j1 = 1, 4
      rh = abs(czc(j1+1)-czc(j1))
      if (rh .eq. 0) go to 201
      fnp = fdf/rh
c
c 1/fnp is number of contours passing through.
      if (fnp .ge. 0.5) go to 201
      dx = abs(cxc(j1+1)-cxc(j1))
      dy = abs(cyc(j1+1)-cyc(j1))
      dmx = dx/(cxma-cxmi)
      dmy = dy/(cyma-cymi)
      pp = 0.005/amax1(dmx,dmy,0.707*(dmx+dmy))
      if (fnp .ge. pp) go to 201
      nr = (rh*pp)/fdf
c
c find closest power-of-two approximation to nr...
      jkk = 0
  210 nr = nr/2
      if (nr .eq. 0) go to 220
      jkk = jkk+1
      go to 210
c
  220 fdf = (2**jkk)*fdf
  201 continue
c
  200 continue
      xn = (czmin-c01)/fdf
      karrz = xn
      if (xn) 21, 22, 22
c
   21 karrz = karrz-1
      int = 0
      go to 23
c
   22 int = 1
   23 xn = karrz
      ctest = c01+xn*fdf
c
c the only way to end this loop is for ctest to exceed czmax.
   24 if (ctest-czmax) 25, 25, 9999
c
   25 call q8plot
   28 ctest = ctest+fdf
      karrz = karrz+1
      if (karrz) 24, 29, 24
c
   29 int = 1
      go to 24
c
c plot from an array of ctest values.
    7 if (k02) 8, 8, 9
c
    8 int = 1
      go to 10
c
    9 int = 0
   10 do 18  j1 = 1, k01
      if (j1-k02) 12, 11, 12
c
   11 int = 1
   12 ctest = bc(j1)
      if (ctest .gt. czmax) go to 18
      if (czmin .gt. ctest) go to 18
      call q8plot
c
   18 continue
      go to 9999
c
c next 3 lines revised 14-dec-78 to conform to 7600 version
   50 print *,"q7plot error: too many contours. change bc(2)    "
      kp = 1
c
9999  return

      END

c**************************************************************************
c
c  q8plot
c
c  calls        trace
c               tracep
c
c**************************************************************************

      subroutine q8plot
c
cc author: f. n. fritsch
c
cc variable declarations:
c
      common /q7quad/ kp, ctest, cxc(5), cyc(5), czc(5)
     * , int, k01, k02, c01, cdf, czmin, czmax
C
c coding for quad begins here...
      dimension ba(5)
      dimension bb(5)
      dimension t(5)
      real*8 ba, bb
      save
c
cc program statements:
c
c load t
      do 1 i = 1, 4
    1 t(i) = czc(i)-ctest
      t(5) = t(1)
      j1 = 1
      isplat = 0
c
c begin loop
      do 9 knumch = 1, 4
    2 if (t(knumch)) 3, 4, 5
c
    3 if (t(knumch+1)) 9, 9, 6
c
    4 xi = cxc(knumch)
      yi = cyc(knumch)
      isplat = isplat+1
      go to 7
c
    5 if (t(knumch+1)) 6, 9, 9
c
c interpolate
    6 tt =   t(knumch)/(t(knumch)-t(knumch+1))
      xi = (cxc(knumch+1)-cxc(knumch))*tt+cxc(knumch)
      yi = (cyc(knumch+1)-cyc(knumch))*tt+cyc(knumch)
    7 ba(j1) = xi
      bb(j1) = yi
      j1 = j1+1
    9 continue
C
c switch for closure and correct number of lines..
   12 go to (9999, 9999, 14, 15, 17), j1
c
   14 lpts = 2
      go to 30
c
   15 lpts = j1
      ba(j1) = ba(1)
      bb(j1) = bb(1)
      go to 30
c
c decide which two lines to plot..
c jump if this quadrilateral is a plateau
   17 if (isplat .eq. 4) go to 9999
      ave = 0.25*(t(1)+t(2)+t(3)+t(4))
      if (ave*t(1)) 20, 15, 22
c
   20 ba(4) = ba(2)
      bb(4) = bb(2)
      ba(2) = xi
      bb(2) = yi
   22 lpts = 2
c
c end of subroutine quad
c
c now plot the contours that have been found...
      if (int .eq. 0) go to 23
      call trace (ba(3), bb(3), 2, -1, -1, 0.0, 0.0)
      go to 31
c
   23 call tracep (ba(3), bb(3), 2, 3, -1, -1)
      go to 32
C
   30 if (int .eq. 0) go to 32
   31 call trace (ba, bb, lpts, -1, -1, 0.0, 0.0)
      go to 9999
c
   32 call tracep (ba, bb, lpts, 3, -1, -1)
c
9999  return

      END

      subroutine line4(x14, y14, x24, y24)
      real*4 x14, y14, x24, y24
      real*8 x18, y18, x28, y28
      x18 = x14
      x28 = x24
      y18 = y14
      y28 = y24
      call line(x18, y18, x28, y28)
      return
      end

      subroutine setcrt4(x4, y4)
      real*4 x4, y4
      real*8 x8, y8
      x8 = x4
      y8 = y4
      call setcrt(x8, y8)
      return
      end

      subroutine vector4(x4, y4)
      real*4 x4, y4
      real*8 x8, y8
      x8 = x4
      y8 = y4
      call vector(x8, y8)
      return
      end

      subroutine point4(x4, y4)
      real*4 x4, y4
      real*8 x8, y8
      x8 = x4
      y8 = y4
      call point(x8, y8)
      return
      end

      subroutine points4(xvec4, yvec4, num, incx, incy, delx4, dely4)
      real*4 xvec4(*), yvec4(*), delx4, dely4
      real*8 xvec8(1000), yvec8(1000), delx8, dely8
      do i = 1, num
         xvec8(i) = xvec4(i)
         yvec8(i) = yvec4(i)
      enddo
      delx8 = delx4
      dely8 = dely4
      call points(xvec8, yvec8, num, incx, incy, delx8, dely8)
      return
      end

      
c**************************************************************************
c
c  gaxisa.F - axis labeling routine
c
c  contents	gaxisa - axis parameter calculations
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

      subroutine gaxisa (fxb, fyb, fxe, fye, iori, isize, iside
     *  , form, idiv, cdivin)
     
      real*8 fxb, fyb, fxe, fye, dcxb, dcyb, dcxe, dcye

      character*(*) form
      character*(*) cdivin(*)
      
      integer lenfor
      character*20 frmstr
c
cc author: mike archuleta dec05.1975
c
cc modifier: steven williams sep26.1980
c
c this subroutine enables the sophisticated user to label
c plots or graphs in just about any method desired. this
c routine is intimately tied to tv80cray and baselib,
c since it assumes the existence of a common block from tv80cray
c (tvgxx1), subroutines setch, crtbcd, and line from tv80cray,
c and subroutines zmovechr, zcetoa, zcftoa, zcitoa and zcotoa
c from baselib.
c assumed inline functions are alog10, sign and sqrt
c
c fxb    the beginning x coordinate of the label
c fyb    the beginning y coordinate of the label
c fxe    the end x coordinate of the label
c fye    the end y coordinate of the label
c iori    the orientation of the text string to be plotted
c         a 0 implies horizontal and a 1 is vertical
c isize   the size of the characters (0-3)
c iside   which side of the axis the label are to be drawn.
c         a 0 implies the left and a 1 implies the right.
c iform   the format of the plotted string. typical uses
c         would be 3ha10, 4hf5.2, 5he20.8, 3hi10, 2ho4
c idiv    the number of labeled marks
c idivin  an array containing information to be plotted
c
c idiv is the key which specifies the type of labeling to do.
c if idiv is negative, then the array adivin contains iabs(idiv)
c values. if the iform is 'a', then the adivin array
c contains the alphanumeric text to be plotted. if the
c iform is not 'a', then the adivin array values are the
c numbers to be used as label at their location.
c
c if idiv is zero, then nice numbered labels are generated
c between adivin(1) and adivin(2).
c
c if idiv is positive, then idiv labels will be plotted
c between adivin(1) and adivin(2).
c
cc variable declarations:
c
      logical klr
      common /tvaxis/ cxb, cyb, cxe, cye, ctickx, cticky
     * , coffx, coffy, cac, cbc, kori, ksize, kn1, klr
     * , klxsav, klysav, cpx(30), cpy(30), kdiv, kepclp
c
      real cx1,cx2,cy1,cy2,cxmi,cxma,cymi,cyma
      integer itype
      dimension cliprect(4)
      dimension bsl(4)
c
c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c
      save bsl
      data bsl(1) /128./
      data bsl(2) /85./
      data bsl(3) /64./
      data bsl(4) /42./
c
c Set library entry routine name

      if (tgname .eq. 'NONE') tgname = 'GAXISA'



c
cc program statements:
c
c store the arguments into local variables
      cxb = fxb
      cyb = fyb
      cxe = fxe
      cye = fye
      kori = iori
      ksize = isize
      kdiv = idiv
      side = -1.
      jside = iside
      if (jside .eq. 1) side = 1.
      
      lenfor = len (form)
c
c Get world and NDC coordinates as well as map type, then set the
c map type to be linear.  Set the allignment to left, centered.
c turn the clipping off so text outside the viewport will show up.
c
      call tggetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,itype)
      call tgsetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,1)
      call gqtxal (ierr,ioldhz,ioldvt)
      call gstxal (1,3)
      call gqclip (ierr,kepclp,cliprect)
      call gsclip (0)
c
c If coordinates in x/y were in log, set klxsav/klysav to 1 else 0
c
      klxsav = 0
      klysav = 0
      if (itype .eq. 3 .or. itype .eq. 4) klxsav = 1
      if (itype .eq. 2 .or. itype .eq. 4) klysav = 1
c
c change the argument which specifies the side of the
c axis the label is to be drawn on if the writing is going
c top to bottom or right to left
c
c The following code uses operands that don't exist in fortran.
c They are rewritten to something that should be equivalent.
c
c      if (kori .eq. 0 .and. cye .lt. cyb) jside = .not. jside
c      if (kori .eq. 1 .and. cxe .lt. cxb) jside = .not. jside
c      klr = .not. (kori .xor. jside)
c
      if (kori .eq. 0 .and. cye .lt. cyb) then
        if (jside .eq. 1) then
	  jside = 0
	else
	  jside = 1
	endif
      endif
      if (kori .eq. 1 .and. cxe .lt. cxb) then
        if (jside .eq. 1) then
	  jside = 0
	else
	  jside = 1
	endif
      endif
      klr = (kori .eq. jside)
c calculate the character offset
      coffx = (cx2-cx1)/((cxma-cxmi)*bsl(ksize+1))
      coffy = (cy2-cy1)/((cyma-cymi)*bsl(ksize+1))
c figure out how long the tick marks should be
      ctickx = (cx2 - cx1) * .01
      cticky = (cy2 - cy1) * .01
c calculate the normal to the axis
      cac = cyb-cye
      cbc = cxe-cxb
      ss = sqrt(cac*cac+cbc*cbc)
      if (ss .eq. 0) go to 64
c normalize the normal and point it in the right direction
      cac = side*cac/ss
      cbc = side*cbc/ss
c
c This version of gaxis requires integer input for the labels
c
c Check the input format and get the length of the format.  If 'I10' is
c passed, then the length is 10.
c
      if (form(1:1) .ne. 'a' .and. form(1:1) .ne. 'A') then
        print *,'GAXISA: invalid format'
	goto 999
      endif
      
      kn1 = 0
      do 3 j = 2,lenfor
        if (form(j:j) .lt. '0' .or. form(j:j) .gt. '9') then
	  print *,'GAXISA: invalid format'
	  goto 999
	endif
	kn1 = kn1 * 10 + ichar (form(j:j)) - 48
 3    continue
      frmstr = '(' // form(1:lenfor) // ')'
c
c jump if the format was bad
c
      if (kn1 .eq. 0) go to 60
      if (kn1 .gt. 80) go to 60
c
c use the array bdiv and its information to plot the label
c
   10 kdiv = iabs(kdiv)
      if (kdiv .gt. 30) go to 61
      if (kdiv .lt. 2) go to 62
c space out the ascii information in equal increments
   16 gxs = (cxe-cxb)/(kdiv-1)
      gys = (cye-cyb)/(kdiv-1)
      cpx(1) = cxb
      cpy(1) = cyb
c fill the value array (bpp), position x array (cpx)
c and the position y array (cpy)
      do 17 j = 2, kdiv
         cpx(j) = cpx(j-1)+gxs
         cpy(j) = cpy(j-1)+gys
   17 continue
c
c Print out the labels
c
   40 continue
      
      do 45 j = 1, kdiv
       call gaxdrw (cpx(j), cpy(j), cdivin(j))
   45 continue
c
c draw the axis line
   50 call tgsetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,itype)
      call line (cxb*1.0D0, cyb*1.0D0, cxe*1.0D0, cye*1.0D0)
      goto 9999
c
c these are the error messages
c
   60 call tgerror (2)
      go to 999
   61 call tgerror (3)
      go to 999
   62 call tgerror (4)
      go to 999
   63 call tgerror (5)
      go to 999
   64 call tgerror (6)
c
 999  call tgsetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,itype)
c
9999  call gsclip (kepclp)
      call gstxal (ioldhz,ioldvt)
      if (tgname .eq. 'GAXISA') tgname = 'NONE'
      return

      END


c**************************************************************************
c
c  gaxisf.F - axis labeling routine
c
c  contents	gaxisf - axis parameter calculations
c		gaxdrw - axis label drawing from calculations
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

      subroutine gaxisf (fxb, fyb, fxe, fye, iori, isize, iside
     *  , form, idiv, adivin)

      real*8 fxb, fyb, fxe, fye, dcxb, dcyb, dcxe, dcye

      character*(*) form
      real*8 adivin(*)
      
c      real adivin(*)
      
      integer lenfor
      character*80 outstr
      character*20 frmstr
c
cc author: mike archuleta dec05.1975
c
cc modifier: steven williams sep26.1980
c
c this subroutine enables the sophisticated user to label
c plots or graphs in just about any method desired. this
c routine is intimately tied to tv80cray and baselib,
c since it assumes the existence of a common block from tv80cray
c (tvgxx1), subroutines setch, crtbcd, and line from tv80cray,
c and subroutines zmovechr, zcetoa, zcftoa, zcitoa and zcotoa
c from baselib.
c assumed inline functions are alog10, sign and sqrt
c
c fxb    the beginning x coordinate of the label
c fyb    the beginning y coordinate of the label
c fxe    the end x coordinate of the label
c fye    the end y coordinate of the label
c iori    the orientation of the text string to be plotted
c         a 0 implies horizontal and a 1 is vertical
c isize   the size of the characters (0-3)
c iside   which side of the axis the label are to be drawn.
c         a 0 implies the left and a 1 implies the right.
c iform   the format of the plotted string. typical uses
c         would be 3ha10, 4hf5.2, 5he20.8, 3hi10, 2ho4
c idiv    the number of labeled marks
c idivin  an array containing information to be plotted
c
c idiv is the key which specifies the type of labeling to do.
c if idiv is negative, then the array adivin contains iabs(idiv)
c values. if the iform is 'a', then the adivin array
c contains the alphanumeric text to be plotted. if the
c iform is not 'a', then the adivin array values are the
c numbers to be used as label at their location.
c
c if idiv is zero, then nice numbered labels are generated
c between adivin(1) and adivin(2).
c
c if idiv is positive, then idiv labels will be plotted
c between adivin(1) and adivin(2).
c
cc variable declarations:
c
      logical klr
      common /tvaxis/ cxb, cyb, cxe, cye, ctickx, cticky
     * , coffx, coffy, cac, cbc, kori, ksize, kn1, klr
     * , klxsav, klysav, cpx(30), cpy(30), kdiv, kepclp
c
      real cx1,cx2,cy1,cy2,cxmi,cxma,cymi,cyma
      integer itype
      dimension cliprect(4)
      dimension bdiv(30)
      dimension bpp(3,35)
      dimension bsk(4)
      dimension bsl(4)

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c
      save bsk
      data bsk(1) /17./
      data bsk(2) /12./
      data bsk(3) /9./
      data bsk(4) /7./
c
      save bsl
      data bsl(1) /128./
      data bsl(2) /85./
      data bsl(3) /64./
      data bsl(4) /42./

c Set library entry routine name

      if (tgname .eq. 'NONE') tgname = 'GAXISF'



c
cc program statements:
c
c store the arguments into local variables
      cxb = fxb
      cyb = fyb
      cxe = fxe
      cye = fye
      kori = iori
      ksize = isize
      kdiv = idiv
      side = -1.
      jside = iside
      if (jside .eq. 1) side = 1.
      
      lenfor = len (form)
c
c Get world and NDC coordinates as well as map type, then set the
c map type to be linear.  Set the allignment to left, centered.
c turn the clipping off so text outside the viewport will show up.
c
      call tggetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,itype)
      call tgsetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,1)
      call gqtxal (ierr,ioldhz,ioldvt)
      call gstxal (1,3)
      call gqclip (ierr,kepclp,cliprect)
      call gsclip (0)
c
c If coordinates in x/y were in log, set klxsav/klysav to 1 else 0
c
      klxsav = 0
      klysav = 0
      if (itype .eq. 3 .or. itype .eq. 4) klxsav = 1
      if (itype .eq. 2 .or. itype .eq. 4) klysav = 1
c
c change the argument which specifies the side of the
c axis the label is to be drawn on if the writing is going
c top to bottom or right to left
c
c The following code uses operands that don't exist in fortran.
c They are rewritten to something that should be equivalent.
c
c      if (kori .eq. 0 .and. cye .lt. cyb) jside = .not. jside
c      if (kori .eq. 1 .and. cxe .lt. cxb) jside = .not. jside
c      klr = .not. (kori .xor. jside)
c
      if (kori .eq. 0 .and. cye .lt. cyb) then
        if (jside .eq. 1) then
	  jside = 0
	else
	  jside = 1
	endif
      endif
      if (kori .eq. 1 .and. cxe .lt. cxb) then
        if (jside .eq. 1) then
	  jside = 0
	else
	  jside = 1
	endif
      endif
      klr = (kori .eq. jside)
c calculate the character offset
      coffx = (cx2-cx1)/((cxma-cxmi)*bsl(ksize+1))
      coffy = (cy2-cy1)/((cyma-cymi)*bsl(ksize+1))
c figure out how long the tick marks should be
      ctickx = (cx2 - cx1) * .01
      cticky = (cy2 - cy1) * .01
c calculate the normal to the axis
      cac = cyb-cye
      cbc = cxe-cxb
      ss = sqrt(cac*cac+cbc*cbc)
      if (ss .eq. 0) go to 64
c normalize the normal and point it in the right direction
      cac = side*cac/ss
      cbc = side*cbc/ss
c
c This version of gaxis requires floating point input for the labels
c
c Check the input format and get the length of the format.  If 'F10' is
c passed, then the length is 10.
c
      if (form(1:1) .ne. 'f' .and. form(1:1) .ne. 'F' .and.
     +    form(1:1) .ne. 'e' .and. form(1:1) .ne. 'E') then
        print *,'GAXISF: invalid format'
	goto 999
      endif
      
      kn1 = 0
      do 3 j = 2,lenfor
        if (form(j:j) .eq. '.') goto 4
        if (form(j:j) .lt. '0' .or. form(j:j) .gt. '9') then
	  print *,'GAXISF: invalid format'
	  goto 999
	endif
	kn1 = kn1 * 10 + ichar (form(j:j)) - 48
 3    continue
 4    frmstr = '(' // form(1:lenfor) // ')'
c
c jump if the format was bad
c
      if (kn1 .eq. 0) go to 60
      if (kn1 .gt. 80) go to 60
c
c jump according to the type of labeling to be done
      if (kdiv) 10, 20, 30
c
c use the array bdiv and its information to plot the label
c
   10 kdiv = iabs(kdiv)
c
      do 11 j = 1, kdiv
        bdiv(j) = adivin(j)
   11 continue
c
c fill the value array (bpp), position x array (cpx)
c and the position y array (cpy)
      divs = bdiv(kdiv)-bdiv(1)
      if (divs .le. 0) go to 63
      do 15 j = 1, kdiv
         bpp(1,j) = bdiv(j)
         gq = (bdiv(j)-bdiv(1))/divs
         cpx(j) = (cxe-cxb)*gq+cxb
         cpy(j) = (cye-cyb)*gq+cyb
   15 continue
      go to 40
c
c generate nice numbered labels between bdiv(1) and bdiv(2)
   20 kdiv = amax1(cxma-cxmi,cyma-cymi)*bsk(ksize+1)
      bdiv(1) = adivin(1)
      bdiv(2) = adivin(2)
   21 gdx = (bdiv(2)-bdiv(1))/kdiv
c jump if bdiv(2) comes before bdiv(1)
      if (gdx .le. 0) go to 63
      gux = alog10(gdx)
      if (gux .lt. 0.) gux = gux-1.
c get the floor of the step size
      j20 = gux
c only use steps of .1, .2, or .5
      g1 = gdx*(10.**(-j20))
      j10 = g1+0.5
      if (j10-5) 22, 25, 23
   22    if (j10-2) 25, 25, 24
   23 j10 = 5
      go to 25
   24 j10 = 2
   25 gdx = j10*(10.**j20)
c change gdx if it would create more than kdiv
c labels.
      if (gdx*kdiv .le. bdiv(2)-bdiv(1)) gdx = 2*gdx
c find the starting point. must be the first nice number
c either on or after the specified starting point
      start = int(bdiv(1)/gdx)*gdx
c     write(*,*) start
      if (start .lt. bdiv(1)) start = start + gdx
c see if we need more than kdiv labels
      kdiv = kdiv/2
   26 end = kdiv*gdx+start
      kdiv = kdiv+1
c the end point has to be the first nice number on or before
c the specified end point
      if (end - bdiv(2)) 26, 27, 265
  265 end = end - gdx
      kdiv = kdiv - 1
c fill the value array (bpp), position x array (cpx)
c and the position y array (cpy)
   27 do 28 j = 1, kdiv
         bpp(1,j) = start+(j-1)*gdx
         gq = (bpp(1,j)-bdiv(1))/(bdiv(2)-bdiv(1))
         cpx(j) = (cxe-cxb)*gq+cxb
         cpy(j) = (cye-cyb)*gq+cyb
   28 continue
      go to 40
c
c generate kdiv labels between bdiv(1) and bdiv(2)
c jump if bad number of divisions
   30 if (kdiv .gt. 30) go to 61
      if (kdiv .lt. 2) go to 62
c calcualte the x and y step sizes
      gxs = (cxe-cxb)/(kdiv-1)
      gys = (cye-cyb)/(kdiv-1)
      bdiv(1) = adivin(1)
      bdiv(2) = adivin(2)
   31 gds = (bdiv(2)-bdiv(1))/(kdiv-1)
c jump if bdiv(2) comes before bdiv(1)
      if (gds .le. 0) go to 63
      bpp(1,1) = bdiv(1)
      cpx(1) = cxb
      cpy(1) = cyb
c fill the value array (bpp), position x array (cpx)
c and the position y array (cpy)
      do 32 j = 2, kdiv
         bpp(1,j) = bpp(1,j-1)+gds
         cpx(j) = cpx(j-1)+gxs
         cpy(j) = cpy(j-1)+gys
   32 continue
c
c Print out the labels
c
   40 continue

c this is the floating point format
      do 45 j = 1, kdiv
	 write (outstr,frmstr) bpp(1,j)
         call gaxdrw (cpx(j), cpy(j), outstr)
   45 continue
c
c draw the axis line
   50 call tgsetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,itype)
      call line (cxb*1.0D0, cyb*1.0D0, cxe*1.0D0, cye*1.0D0)
      goto 9999
c
c these are the error messages
c
   60 call tgerror (2)
      go to 999
   61 call tgerror (3)
      go to 999
   62 call tgerror (4)
      go to 999
   63 call tgerror (5)
      go to 999
   64 call tgerror (6)
c
 999  call tgsetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,itype)
c
9999  call gsclip (kepclp)
      call gstxal (ioldhz,ioldvt)
      if (tgname .eq. 'GAXISF') tgname = 'NONE'
      return

      END


# 1 "gaxisi.F" 
c**************************************************************************
c
c  gaxisi.F - axis labeling routine
c
c  contents	gaxisi - axis parameter calculations
c		gaxdrw - axis label drawing from calculations
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

# 428 "gaxisi.F" 

      subroutine gaxisi (fxb, fyb, fxe, fye, iori, isize, iside
     *  , form, idiv, idivin)
     
      real*8 fxb, fyb, fxe, fye, dcxb, dcyb, dcxe, dcye

      character*(*) form
      integer idivin(*)
      
      integer lenfor
      character*80 outstr
      character*20 frmstr
c
cc author: mike archuleta dec05.1975
c
cc modifier: steven williams sep26.1980
c
c this subroutine enables the sophisticated user to label
c plots or graphs in just about any method desired. this
c routine is intimately tied to tv80cray and baselib,
c since it assumes the existence of a common block from tv80cray
c (tvgxx1), subroutines setch, crtbcd, and line from tv80cray,
c and subroutines zmovechr, zcetoa, zcftoa, zcitoa and zcotoa
c from baselib.
c assumed inline functions are alog10, sign and sqrt
c
c fxb    the beginning x coordinate of the label
c fyb    the beginning y coordinate of the label
c fxe    the end x coordinate of the label
c fye    the end y coordinate of the label
c iori    the orientation of the text string to be plotted
c         a 0 implies horizontal and a 1 is vertical
c isize   the size of the characters (0-3)
c iside   which side of the axis the label are to be drawn.
c         a 0 implies the left and a 1 implies the right.
c iform   the format of the plotted string. typical uses
c         would be 3ha10, 4hf5.2, 5he20.8, 3hi10, 2ho4
c idiv    the number of labeled marks
c idivin  an array containing information to be plotted
c
c idiv is the key which specifies the type of labeling to do.
c if idiv is negative, then the array adivin contains iabs(idiv)
c values. if the iform is 'a', then the adivin array
c contains the alphanumeric text to be plotted. if the
c iform is not 'a', then the adivin array values are the
c numbers to be used as label at their location.
c
c if idiv is zero, then nice numbered labels are generated
c between adivin(1) and adivin(2).
c
c if idiv is positive, then idiv labels will be plotted
c between adivin(1) and adivin(2).
c
cc variable declarations:
c
      logical klr
      common /tvaxis/ cxb, cyb, cxe, cye, ctickx, cticky
     * , coffx, coffy, cac, cbc, kori, ksize, kn1, klr
     * , klxsav, klysav, cpx(30), cpy(30), kdiv, kepclp
c
      real cx1,cx2,cy1,cy2,cxmi,cxma,cymi,cyma
      integer itype
      dimension cliprect(4)
      dimension bdiv(30)
      dimension bpp(3,35)
      dimension bsk(4)
      dimension bsl(4)
c      dimension bvstr(3)
c      dimension ndiv(30)
c      dimension npp(3,35)
c      dimension nstr(3)
c
c      equivalence (bpp, npp)
c
c      equivalence (bdiv(1), ndiv(1))
c
c      equivalence (bvstr(1), nstr(1))
c
c      equivalence (gq, jq)
c
c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,
     +    scaley,rotat,centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
      logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c
      save bsk
      data bsk(1) /17./
      data bsk(2) /12./
      data bsk(3) /9./
      data bsk(4) /7./
c
      save bsl
      data bsl(1) /128./
      data bsl(2) /85./
      data bsl(3) /64./
      data bsl(4) /42./
c
c      save jbpc
c      data jbpc /8/
c      save jcbla
c      data jcbla /#20/
c      data jcbla /32/
c      save jcnul
c      data jcnul /#00/
c      data jcnul /0/
c      save jcpw
c      data jcpw /8/

c Set library entry routine name

      if (tgname .eq. 'NONE') tgname = 'gaxisi'

c
cc program statements:
c
c store the arguments into local variables
      cxb = fxb
      cyb = fyb
      cxe = fxe
      cye = fye
      kori = iori
      ksize = isize
      kdiv = idiv
      side = -1.
      jside = iside
      if (jside .eq. 1) side = 1.
      
      lenfor = len (form)
c
c Get world and NDC coordinates as well as map type, then set the
c map type to be linear.  Set the allignment to left, centered.
c turn the clipping off so text outside the viewport will show up.
c
      call tggetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,itype)
      call tgsetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,1)
      call gqtxal (ierr,ioldhz,ioldvt)
      call gstxal (1,3)
      call gqclip (ierr,kepclp,cliprect)
      call gsclip (0)
c
c If coordinates in x/y were in log, set klxsav/klysav to 1 else 0
c
      klxsav = 0
      klysav = 0
      if (itype .eq. 3 .or. itype .eq. 4) klxsav = 1
      if (itype .eq. 2 .or. itype .eq. 4) klysav = 1
c
c change the argument which specifies the side of the
c axis the label is to be drawn on if the writing is going
c top to bottom or right to left
c
c The following code uses operands that don't exist in fortran.
c They are rewritten to something that should be equivalent.
c
c      if (kori .eq. 0 .and. cye .lt. cyb) jside = .not. jside
c      if (kori .eq. 1 .and. cxe .lt. cxb) jside = .not. jside
c      klr = .not. (kori .xor. jside)
c
      if (kori .eq. 0 .and. cye .lt. cyb) then
        if (jside .eq. 1) then
	  jside = 0
	else
	  jside = 1
	endif
      endif
      if (kori .eq. 1 .and. cxe .lt. cxb) then
        if (jside .eq. 1) then
	  jside = 0
	else
	  jside = 1
	endif
      endif
      klr = (kori .eq. jside)
c calculate the character offset
      coffx = (cx2-cx1)/((cxma-cxmi)*bsl(ksize+1))
      coffy = (cy2-cy1)/((cyma-cymi)*bsl(ksize+1))
c figure out how long the tick marks should be
      ctickx = (cx2 - cx1) * .01
      cticky = (cy2 - cy1) * .01
c calculate the normal to the axis
      cac = cyb-cye
      cbc = cxe-cxb
      ss = sqrt(cac*cac+cbc*cbc)
      if (ss .eq. 0) go to 64
c normalize the normal and point it in the right direction
      cac = side*cac/ss
      cbc = side*cbc/ss
c
c This version of gaxis requires integer input for the labels
c
c Check the input format and get the length of the format.  If 'I10' is
c passed, then the length is 10.
c
      if (form(1:1) .ne. 'i' .and. form(1:1) .ne. 'I' .and.
     +    form(1:1) .ne. 'o' .and. form(1:1) .ne. 'O') then
        print *,'GAXISI: invalid format'
	return
      endif
      
      kn1 = 0
      do 3 j = 2,lenfor
        if (form(j:j) .lt. '0' .or. form(j:j) .gt. '9') then
	  print *,'GAXISI: invalid format'
	  return
	endif
	kn1 = kn1 * 10 + ichar (form(j:j)) - 48
 3    continue
      frmstr = '(I' // form(2:lenfor) // ')'
c
c jump if the format was bad
c
      if (kn1 .eq. 0) go to 60
      if (kn1 .gt. 80) go to 60
c
c jump according to the type of labeling to be done
      if (kdiv) 10, 20, 30
c
c The idivin array contains iabs(kdiv) elements which are to be plotted
c at even intervals along the axis.
c
   10 kdiv = iabs(kdiv)
c
c jump if bad number of divisions
c
      if (kdiv .gt. 30) go to 61
      if (kdiv .lt. 2) go to 62
c
        if (form(1:1) .eq. 'o' .or. form(1:1) .eq. 'O') then
      do 11 j = 1, kdiv
c * * * fix for compiler bug
C * * *          jj = (j-1)*2 + 1
c * * *	  jnum = idivin(jj)
          jnum = idivin(j)
          jmod = 8
          jfact = 1
          jnew = 0
  112     continue
          jdiv = jnum/jmod+.5
          jdif = jnum-jdiv*jmod
          jnew = jnew+jdif*jfact
          if (jdiv .gt. 0) then
            jnum = jdiv
            jfact = jfact*10
            go to 112
          endif
          bdiv(j) = jnew
   11 continue
        else
      do 12 j = 1, kdiv
c * * * fix for compiler bug
c * * *       jj = (j-1)*2 + 1
c * * *	  bdiv(j) = idivin(jj)
          bdiv(j) = idivin(j)
   12 continue
        endif
c
c fill the value array (bpp), position x array (cpx)
c and the position y array (cpy)
      divs = bdiv(kdiv)-bdiv(1)
      if (divs .le. 0) go to 63
      do 15 j = 1, kdiv
         bpp(1,j) = bdiv(j)
         gq = (bdiv(j)-bdiv(1))/divs
         cpx(j) = (cxe-cxb)*gq+cxb
         cpy(j) = (cye-cyb)*gq+cyb
   15 continue
      go to 40




c
c generate nice numbered labels between bdiv(1) and bdiv(2)
   20 kdiv = amax1(cxma-cxmi,cyma-cymi)*bsk(ksize+1)
      bdiv(1) = idivin(1)
c * * * fix for compiler bug
c * * *      bdiv(2) = idivin(3)
      bdiv(2) = idivin(2)
      gdx = (bdiv(2)-bdiv(1))/kdiv
c jump if bdiv(2) comes before bdiv(1)
      if (gdx .le. 0) go to 63
      gux = alog10(gdx)
      if (gux .lt. 0.) gux = gux-1.
c get the floor of the step size
      j20 = gux
c only use steps of .1, .2, or .5
      g1 = gdx*(10.**(-j20))
      j10 = g1+0.5
      if (j10-5) 22, 25, 23
   22    if (j10-2) 25, 25, 24
   23 j10 = 5
      go to 25
   24 j10 = 2
   25 gdx = j10*(10.**j20)
c change gdx if it would create more than kdiv
c labels.
      if (gdx*kdiv .le. bdiv(2)-bdiv(1)) gdx = 2*gdx
c find the starting point. must be the first nice number
c either on or after the specified starting point
      start = int(bdiv(1)/gdx)*gdx
      if (start .lt. bdiv(1)) start = start + gdx
c see if we need more than kdiv labels
      kdiv = kdiv/2
   26 end = kdiv*gdx+start
      kdiv = kdiv+1
c the end point has to be the first nice number on or before
c the specified end point
      if (end - bdiv(2)) 26, 27, 265
  265 end = end - gdx
      kdiv = kdiv - 1
c fill the value array (bpp), position x array (cpx)
c and the position y array (cpy)
   27 do 28 j = 1, kdiv
         bpp(1,j) = start+(j-1)*gdx
         gq = (bpp(1,j)-bdiv(1))/(bdiv(2)-bdiv(1))
         cpx(j) = (cxe-cxb)*gq+cxb
         cpy(j) = (cye-cyb)*gq+cyb
   28 continue
      go to 40
c
c generate kdiv labels between bdiv(1) and bdiv(2)
c jump if bad number of divisions
   30 if (kdiv .gt. 30) go to 61
      if (kdiv .lt. 2) go to 62
c calcualte the x and y step sizes
      gxs = (cxe-cxb)/(kdiv-1)
      gys = (cye-cyb)/(kdiv-1)
      bdiv(1) = idivin(1)
c * * * fix for compiler bug
c * * *      bdiv(2) = idivin(3)
      bdiv(2) = idivin(2)
   31 gds = (bdiv(2)-bdiv(1))/(kdiv-1)
c jump if bdiv(2) comes before bdiv(1)
      if (gds .le. 0) go to 63
      bpp(1,1) = bdiv(1)
      cpx(1) = cxb
      cpy(1) = cyb
c fill the value array (bpp), position x array (cpx)
c and the position y array (cpy)
      do 32 j = 2, kdiv
         bpp(1,j) = bpp(1,j-1)+gds
         cpx(j) = cpx(j-1)+gxs
         cpy(j) = cpy(j-1)+gys
   32 continue
c
c Print out the labels
c
   40 continue

c this is the integer format
      do 45 j = 1, kdiv
         j1 = bpp(1,j) + .05
	 write (outstr,frmstr) j1
         call gaxdrw (cpx(j), cpy(j), outstr)
   45 continue
c
c draw the axis line
   50 call tgsetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,itype)
      dcxb = cxb
      dcyb = cyb
      dcxe = cxe
      dcye = cye
      call line (dcxb, dcyb, dcxe, dcye)
      goto 9999
c
c these are the error messages
c
   60 call tgerror (2)
      go to 999
   61 call tgerror (3)
      go to 999
   62 call tgerror (4)
      go to 999
   63 call tgerror (5)
      go to 999
   64 call tgerror (6)
c
 999  call tgsetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,itype)
c
9999  call gsclip (kepclp)
      call gstxal (ioldhz,ioldvt)
      if (tgname .eq. 'gaxisi') tgname = 'NONE'
      return
      END


# 929 "gaxisi.F" 




      subroutine gaxdrw (gxp, gyp, txt)
      
      real   gxp, gyp
      real*8 dgxp,dgyp,dgxs,dgys
      character*(*) txt
c
c this subroutine will draw the text and the tick mark for gaxis
c
      logical klr
      common /tvaxis/ cxb, cyb, cxe, cye, ctickx, cticky
     * , coffx, coffy, cac, cbc, kori, ksize, kn1, klr
     * , klxsav, klysav, cpx(30), cpy(30), kdiv, kepclp
c
cc program statements:
c
      if (klxsav .ne. 0) gxp = alog10(gxp)
      if (klysav .ne. 0) gyp = alog10(gyp)
c draw the tick mark
      call line (gxp*1.0D0, gyp*1.0D0, 
     *  (gxp+cac*ctickx)*1.0D0, (gyp+cbc*cticky)*1.0D0)
c jump if horizontal text on a horizontal axis
      if (kori .eq. 0 .and. abs(cac) .le. .01) go to 1
c jump if vertical text on a vertical axis
      if (kori .eq. 1 .and. abs(cbc) .le. .01) go to 2
c get the starting position of the text
      gxs = gxp+cac*coffx
      gys = gyp+cbc*coffy
c change the starting position if the orientation requires it
      if (klr) go to 3
      if (kori .eq. 0) gxs = gxs-kn1*coffx
      if (kori .eq. 1) gys = gys-kn1*coffy
         go to 3
c
    1    gxs = gxp-kn1*coffx/2.
         gys = gyp+2.*sign(coffy,cbc)
         go to 3
c
    2    gxs = gxp+2.*sign(coffx,cac)
         gys = gyp-kn1*coffy/2.
c position the beam and draw the text
    3 dgxs = gxs
      dgys = gys
      call setlch (dgxs, dgys, 0, ksize, kori, -1)
      call gtext (txt, kn1, 0)
c
 9999 return
      END

c      block data
c      common / hatchcm0 / ndimen
c      data ndimen /101/
c      end

c**************************************************************************
c
c  hatch.F - area shading subroutine
c
c  contents	hatch - draws lines in specified area
c		hatchmap
c		hatchunm
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

      subroutine hatch( xvert, yvert, npoints, phi, ngrid, mode )
C
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     H A T C H
c
c     Provide shading for a general polygonal region.  There is absolutely no
c     assumption made about convexity.  A polygon is specified by its vertices,
c     given in either a clockwise or counter-clockwise order.  The density of
c     the shading lines (or points) and the angle for the shading lines are
c     both determined by the parameters passed to the subroutine.
c
c     The calling sequence is:
c
c        call hatch( xvert, yvert, npoints, phi, ngrid, mode )
c
c     The input parameters are interpreted as follows:
c
c        xvert    -  An array or x coordinates for the polygon vertices
c
c        xvert    -  An array or x coordinates for the polygon vertices
c
c        npoints  -  The number of vertices in the polygon
c
c        phi      -  The angle for the shading, measured counter-clockwise
c                    in radians from the positive x-axis
c
c        ngrid    -  A parameter  determining the shading density:
c                      -n => every n-th shading line
c                        ....
c                      -2 => every other shading line
c                      -1 => every shading line
c                       0 => every shading line
c                      +1 => every other shading line
c                      +2 => every fourth shading line
c                        ....
c                      +n => every 2**n-th shading line
c
c        mode     -  A parameter determining the shading mode:
c                      -1 =>    no shading - boundary drawn
c                       0 =>  line-shading - boundary not drawn
c                       1 => point-shading - boundary not drawn
c                       2 =>  line-shading - boundary drawn
c                       3 => point-shading - boundary drawn
c
c     All coordinates are assumed to be in the user coordinate system as
c     specified in the most recent call to the tv80lib mapping subroutines.
c     Either Cartesian or polar coordinates are acceptable, although weird
c     (but aesthetic) results are produced using log-log or semi-log plots.
c     Direct tv80lib raster units (integer) may also be used.
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      integer npoints, ngrid
c     real xvert(npoints), yvert(npoints), phi
      logical points,bound

      common /hacomm/ sinphi,cosphi

c
c     This subroutine has to maintain an internal array of the transformed
c     coordinates.  This establishes a storage limitation.  The parameter
c     "nsize" is the maximum number of vertices allowed.  If the user wants
c     to increase the size of the arrays, make sure that the three common
c     blocks - hatchcom1, hatchcom2, hatchcom3 - are declared before this
c     routine and that the variable "ndimen" which is in common block
c     - hatchcom0 - is initialized at run time (not data-loaded) to the
c     correct dimension.  The maximum number of vertices is one less than
c     the dimension of the work arrays.
C
c      parameter ( nsize = 101 )
C
      common / hatchcm0 / ndimen
c      data ndimen /101/
      common / hatchcm1 / cxvert(101)
      common / hatchcm2 / cyvert(101)
      common / hatchcm3 / xintercept(101)

c include the standard common block

      dimension xvert(*), yvert(*)
      real*8 xvert, yvert, phi, x1, x2, y1, y2

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



C

c Set library entry routine name
        ndimen = 101
      if (tgname .eq. 'NONE') tgname = 'HATCH'
	


c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
c     Check for valid number of vertices.
c
      if (npoints .lt. 2 .or. npoints .gt. ndimen-1) then
        goto 9999
      endif
C
c     See if point plotting was specified and whether to draw the
c     outline of the polygons.  Set flags for later.
C
      points = (mode .eq. 1) .or. (mode .eq. 3)
      bound  = (mode .eq. 2) .or. (mode .eq. 3) .or. (mode .eq. -1)
c
c If the last vertex is not the same as the first, set nvert to npoints+1
c to signify that we need to add a vertex.
c
      if (xvert(npoints) .ne. xvert(1) .or.
     +    yvert(npoints) .ne. yvert(1)) then
        nvert = npoints + 1
      else
        nvert = npoints
      endif
c
c Draw the border if required.  Return if no other processing is required.
c
      if (bound) then
        call trace (xvert,yvert,npoints,-1,-1,0.0,0.0)
        if (nvert .ne. npoints)
     +    call line (xvert(1),yvert(1),xvert(npoints),yvert(npoints))
      endif
      if (mode .eq. -1) goto 9999
c
c Save the old mapping and set values for rotation.
c
      call tggetmap (vl,vr,vb,vt,wl,wr,wb,wt,it)
      sinphi = sin(phi)
      cosphi = cos(phi)
c
c Convert the users coordinate first into linear and then into
c normalized coordinates.  Various problems arise when the users
c input values have a very high or very low magnitude.  The easiest
c way to solve these problems is to convert the coordinates into
c normalized coordinats.
c

      do 5 i = 1,npoints
        call culxy (xvert(i),yvert(i),cxvert(i),cyvert(i))
        call clnxy (cxvert(i),cyvert(i),cxvert(i),cyvert(i))
5     continue

c
c Now that our input coordinate are all normalized, change the mapping
c into a normalized map.
c

      call tgsetmap (vl,vr,vb,vt,vl,vr,vb,vt,1)

c
c Compute the spacing from the parameters.  The finest spacing is 1024
c points across the screen.  This corresponds to either an input of
c ngrid = -1 or ngrid = 0.
c
      if (ngrid .ge. 0) then
        step = real(2 ** (min0(ngrid,10))) / 1023.
      else
        step = real(abs (ngrid)) / 1023.
      endif
c
c Find the maximum and minimum x and y coordinates after rotation.  The
c user specifies an angle of rotation for the hatch lines.  It is easier
c to rotate the vertices and then compute horizontal lines in some sort
c of scanline algorithm.  The resulting 'scan' lines are then unrotated
c and drawn to produce the desired shading lines.  hatchmap will rotate
c the vertex. hatchunm will unrotate.
c
      call hatchmap (cxvert(1),cyvert(1))
      xmin = cxvert(1)
      xmax = xmin
      ymin = cyvert(1)
      ymax = ymin
      do 10 i = 2,npoints
        call hatchmap (cxvert(i),cyvert(i))
        xmax = max (cxvert(i),xmax)
        xmin = min (cxvert(i),xmin)
        ymax = max (cyvert(i),ymax)
        ymin = min (cyvert(i),ymin)
10    continue
      if (nvert .ne. npoints) then
        cxvert(nvert) = cxvert(1)
        cyvert(nvert) = cyvert(1)
      endif
c
c Do the scan line algorithm on the rotated vertices stored in cxvert
c and cyvert.
c
      do 20 y = ymin,ymax,step
        ivert = 0
        icount = 0
        do 30 i = 1,nvert-1
          yhead = y - cyvert(i+1)
          ytail = y - cyvert(i)
          if (sign(1.,yhead) .ne. sign(1.,ytail)) then
            icount = icount + 1
            delx = cxvert(i+1) - cxvert(i)
            dely = cyvert(i+1) - cyvert(i)
            delta = delx/dely * yhead
            xintercept (icount) = delta + cxvert(i+1)
          endif
30      continue

c
c       Sort the x intercept values.  Use a bubblesort because there aren't
c       very many of them (usually only two).
c

        do 40 i = 1,icount
          xkey = xintercept(i)
          do 50 j = 1,i-1
            if (xintercept(j) .gt. xkey) then
              xtemp = xkey
              xkey = xintercept(j)
              xintercept(j) = xtemp
            endif
50        continue
          xintercept(i) = xkey
40      continue

C
c All of the x coordinates for the shading segments along the current
c shading line are now known and are in sorted order.  All that remains
c is to draw them.  Process the x coordinates two at a time.
C
        do 60 i = 1, icount, 2
C
          x1 = xintercept(i)
          x2 = xintercept(i+1)
          y1 = y
          y2 = y
c
c Rotate back to original direction.
c
          call hatchunm (x1, y1)
          call hatchunm (x2, y2)
C
c See if plotting lines or points.
c
          if (points) then
            call linep (x1,y1, x2,y2, ngrid)
          else
            call line (x1, y1, x2, y2)
          endif
60      continue
20    continue

c
c Reset the map to what it was before hatch screwed it up. All of
c the 'goto 9999' statements in this subroutine are performed before
c the change to the mapping. So 9999 is places after the map reset.
c

      call tgsetmap (vl,vr,vb,vt,wl,wr,wb,wt,it)

9999  if (tgname .eq. 'HATCH') tgname = 'NONE'
      return

      end
 
C
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     H A T C H U N M A P
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      subroutine hatchunm (xmap, ymap)
      real*8 xmap, ymap, xunmap, yunmap
      common /hacomm/ sinphi,cosphi
      
c Rotate through an angle of phi to original orientation.

      xunmap = cosphi*xmap - sinphi*ymap
      yunmap = sinphi*xmap + cosphi*ymap
      xmap = xunmap
      ymap = yunmap
      return

      END
C
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     H A T C H M A P
c
c     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      subroutine hatchmap (xuser, yuser)
      real xuser,yuser,xmap,ymap
      common /hacomm/ sinphi,cosphi
c
c Rotate through an angle of -phi.
c
      xmap = cosphi*xuser + sinphi*yuser
      ymap = - sinphi*xuser + cosphi*yuser
      xuser = xmap
      yuser = ymap
      return

      END
c**************************************************************************
c
c  mapdrw.f - drawing routine for gaxis
c
c  contents	mapdrw - draws text and tick marks for gaxis
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************


c**************************************************************************

	subroutine mapdrw (string,numchr,ioffset,axis,chrpos,isx,itic)
	character*(*) string
	integer numchr,ioffset
	real axis,chrpos
	integer isx,itic
	real*8 x1, y1, x2, y2
c
c numchr - number of characters in the string
c
c ioffset - offset in the string to the first character
c
c axis - value on the axis to print on
c
c chrpos - The character position at axis
c
c isx - Axis number
c       = 0 (y axis)
c       = 1 (x axis)
c       = 2 (polar axis)
c
c itic - Print tic or grid
c       = 0 (print grid line)
c       = 1 (print tic mark)
c       = 2 (omit marks)
c

c Get the current mapping, the type (it) has to be linear

	call tggetmap (xvmin,xvmax,yvmin,yvmax,
     +    xwmin,xwmax,ywmin,ywmax,it)

c Check to see if we are plotting on the x axis

	if (isx .eq. 1) then
	  sdata = (axis-xwmin) * (xvmax-xvmin)/(xwmax-xwmin) + xvmin
	  sdata = (sdata * 42.5) + 1.5
	  x1 = sdata
	  y1 = chrpos
	  call setch (x1,y1,0,1,1,-1)
cccccc          call setch (sdata,chrpos,0,1,1,-1)
	  call gtext (string,numchr,ioffset)
	  if (itic .eq. 0) then
	    x1 = axis
	    x2 = axis
	    y1 = ywmin
	    y2 = ywmax
	    call line (x1, y1, x2, y2)
cccccc            call line (axis,ywmin,axis,ywmax)
	  else if (itic .eq. 1) then
	    ticlen = .01 * (ywmax-ywmin)
	    x1 = axis
	    x2 = axis
	    y1 = ywmin
	    y2 = ywmin+ticlen
	    call line (x1, y1, x2, y2)
cccccc            call line (axis,ywmin,axis,ywmin+ticlen)
	    y1 = ywmax
	    y2 = ywmax-ticlen
	    call line (x1, y1, x2, y2)
cccccc            call line (axis,ywmax,axis,ywmax-ticlen)
	  endif

c Check to see if we are plotting on the y axis

	else if (isx .eq. 0) then
	  sdata = (axis-ywmin) * (yvmax-yvmin)/(ywmax-ywmin) + yvmin
	  sdata = (sdata * 42.5) + .5
	  x1 = chrpos
	  y1 = sdata
	  call setch (x1,y1,0,1,0,-1)
cccccc          call setch (chrpos,sdata,0,1,0,-1)
	  call gtext (string,numchr,ioffset)
	  if (itic .eq. 0) then
	    x1 = xwmin
	    y1 = axis
	    x2 = xwmax
	    y2 = axis
	    call line (x1, y1, x2, y2)
cccccc            call line (xwmin,axis,xwmax,axis)
	  else if (itic .eq. 1) then
	    ticlen = .01 * (xwmax-xwmin)
	    x1 = xwmin
	    y1 = axis
	    x2 = xwmin+ticlen
	    y2 = axis
	    call line (x1, y1, x2, y2)
cccccc            call line (xwmin,axis,xwmin+ticlen,axis)
	    x1 = xwmax
	    x2 = xwmax-ticlen
	    call line (x1, y1, x2, y2)
cccccc            call line (xwmax,axis,xwmax-ticlen,axis)
	  endif

c We are not on the x or y axis, we are plotting polar, put ticks in x direction

	else
	  sdata = (axis-xwmin) * (xvmax-xvmin)/(xwmax-xwmin) + xvmin
	  sdata = (sdata * 42.5) + 1.5
	  x1 = sdata
	  y1 = chrpos
	  call setch (x1,y1,0,1,1,-1)
cccccc          call setch (sdata,chrpos,0,1,1,-1)
	  call gtext (string,numchr,ioffset)
	  if (itic .eq. 1) then
	    ticlen = .01 * (ywmax-ywmin)
	    x1 = axis
	    y1 = 0.
	    x2 = axis
	    y2 = ticlen
	    call line (x1, y1, x2, y2)
cccccc            call line (axis,0.,axis,ticlen)
	  endif
	endif

	return
	end




      function log10a(f1)
c
cc purpose: get an integer approximation of an antilogarithm
c
cc program statements:
c
      g1 = f1
      if (f1 .lt. 0.0) g1 = -f1
      if (f1 .eq. 0.0) then
        print *,'log10awarning: trying to take antilog'
        print *,'of 0.0 - will return 0'
        log10a= 0
        return
      endif
c
c get antilog
c
      g1 = alog10(g1)
      log10a= 0
      if (g1 .lt. 0.0) log10a= g1-0.9999999999
      if (g1 .ge. 0.0) log10a= g1+0.0000000001
c
      return
      END


      function inta (f1)
c
cc purpose: to get an integer approximation of a floater
c
cc program statements:
c
      inta = int(f1+0.9999999999)
      if (f1 .lt. 0.0) inta = int(f1-0.9999999999)
      return

      END

      function icomp (f1, f2)
c
cc purpose: to determine the numeric relationship
c           between two floaters
c
cc program statements:
c
      icomp = 0
      if (f1 .eq. f2) return
      icomp = 1
      if (f1 .lt. f2) icomp = -1
c
c see if floaters are equivalent within a reasonable value
      g1 = 1.00001
      if (icomp .eq. 1 .and. 0.0 .le. f1
     * .or. icomp .eq. -1 .and. f1 .lt. 0.0)
     * g1 = 0.99999
      g2 = g1*f1
      if (icomp .eq. 1 .and. g2 .lt. f2
     * .or. icomp .eq. -1 .and. f2 .le. g2)
     * icomp = 0
c
      return
      END


      subroutine maplab (fmin, fmax, islog, isx, itic)
c=======================================================================
c 12/19/2001 ler pppl - fix rounding when shifting digits by dividing by
c                        10.0d0 instead of 10.0
c=======================================================================

	real*8 x1, y1, x2, y2
c
cc purpose: set number of divisions and get ascii value
c           that should be displayed at tick marks
c
cc revised: oct20.1980 by steven williams
c
cc arguments:
c
c fmin = minimum x or y axis value
c
c fmax = maximum x or y axis value
c
c islog = axis log indicator
c      = 1 (do log scaling)
c      = 0 (do linear scaling)
c
c isx = x or y axis
c      = 2 (polar x axis)
c      = 1 (x axis)
c      = 0 (y axis)
c
c itic = print tic marks
c      = 1 (print tics)
c      = 0 (omit tics)
c
cc enhancements:
c
c gdx = distance between successive tick marks
c gend = axis end point
c gsavsp = saved starting point coordinate
c gstart = axis start point
c jblaw = word filled with ascii blank characters
c jdshl = number of character positions to shift decimal point left
c jdxe = integer exponent of gdx
c jende = integer exponent of gend
c jlimch = "limit characters"
c        = upper limit on desirable number of characters to print
c          in tick mark labels
c jmaxch = "maximum characters"
c        = maximum number of characters that can be printed
c          in a tick mark label
c jmaxe = maximum of jstare and jende
c jnice = 0 (=> mean data - need an exponent factor at bottom of axis)
c       = 1 (=> nice data - dont need an exponent factor at bottom of axis)
c jprinn = number of trailing characters not to print in tick labels
c jprint = base for number of characters to print in tick labels
c jsavst = saved string
c jstare = integer exponent of gstart
c jzeroe = "zero end"
c        = -1 (=> gend .lt. 0.0)
c        = 0  (=> gend .eq. 0.0)
c        = 1  (=> gend .gt. 0.0)
c jzeros = "zero start"
c        = -1 (=> gstart .lt. 0.0)
c        = 0  (=> gstart .eq. 0.0)
c        = 1  (=> gstart .gt. 0.0)
c ndotz = "dot zeros"
c       = a 2-word character string starting with a dot
c         and followed by zeros
c
cc variable declarations:
c
c        integer kascii(3),kprint,kskip,knp
        integer kprint,kskip,knp
	character*20 kstring
        real cdata,cspoin
c
c      dimension ndotz(2)
	character*20 ndotz
c      dimension ns(10)
	character*20 ns
      dimension cliprect(4)
c
      save
c      data jblaw /"        "/
      data jlimch /6/
c can only have 12 digits to right of decimal point with baselib
c subroutine zcetoa - need 7 more for sign, one digit to left of
c decimal point, decimal point and e+## trailer phrase => 19 characters
      data jmaxch /18/
      data ndotz /'.0000000000000000000'/
      data ccpl1 /42.5/
c
c Turn off clipping so the labels will appear
c
        call gqclip (ierr,kepclp,cliprect)
        call gsclip (0)
c
c Do the initial setup and draw a border around the window
c
        call tggetmap (vl,vr,vb,vt,wl,wr,wb,wt,it)
        call tgsetmap (vl,vr,vb,vt,wl,wr,wb,wt,1)

        if (isx .eq. 2) then
          knp = (vr-vl) * 11.
          cspoin = ccpl1*(vt+vb)-5.0
          fscale = .5 * (vt-vb)/(wr+wr)
          call line4 (wl,(wt+wb)/2.,wr,(wt+wb)/2.)
          call line4 ((wr+wl)/2.,wb,(wr+wl)/2.,wt)
        else
          if (isx .eq. 1) then
            knp = (vr-vl) * 13.
            cspoin = 2*ccpl1*vb-4.0
            fscale = (vr-vl)/(wr-wl)
          else
            knp = (vt-vb) * 13.
            cspoin = 2*ccpl1*vl-4.0
            fscale = (vb-vt)/(wb-wt)
          endif
          wl8 = wl
          wr8 = wr
          wb8 = wb
          wt8 = wt
          call setcrt4 (wl,wb)
          call vector4 (wr,wb)
          call vector4 (wr,wt)
          call vector4 (wl,wt)
          call vector4 (wl,wb)
        endif

        if (cspoin .lt. 1) cspoin = 1.0

c jump to 5000 if log axis scaling is to be performed

      if (islog .eq. 1) go to 5000
c
      if (knp .le. 0) knp = 1
      gdx = (fmax-fmin)/knp
c
c get the floor of the step size
      jdxe = log10a(gdx)
c
c we will only use step sizes of 0.1, 0.2 and 0.5
      i1 = gdx*(10.0**(-jdxe))+0.000001
      if (i1-5) 1100, 1400, 1200
 1100 if (i1-2) 1400, 1400, 1300
 1200 i1 = 5
      go to 1400
c
 1300 i1 = 2
 1400 gdx = i1*(10.0**jdxe)
c
c get the starting position of the text
      gstart = inta(fmin/gdx)*gdx
      if (icomp(gstart,fmin) .eq. -1) gstart = gstart+gdx
c
c get the new number of divisions
      knp = inta((fmax-gstart)/gdx)
      go to 1550
 1500 knp = knp+1
 1550 gend = knp*gdx+gstart
      if (knp .gt. 34) go to 1600
      if (icomp(gend,fmax)) 1500, 1700, 1600
 1600 knp = knp-1
c
c get exponent of start point
 1700 jzeros = icomp(gstart,0.0)
      if (jzeros .eq. 0) jstare = 0
      if (jzeros .ne. 0) jstare = log10a(gstart)
c
c get exponent of end point
      jzeroe = icomp(gend,0.0)
      if (jzeroe .eq. 0) jende = 0
      if (jzeroe .ne. 0) jende = log10a(gend)
c
c specify mean or nice numbers
      jnice = 1
      if ((jdxe .le. -(jlimch))
     * .or. (jdxe .eq. -(jlimch-1) .and. gstart .lt. 0.0)
     * .or. (jdxe .eq. -(jlimch-1) .and. gend .lt. 0.0)
     * .or. (jstare .eq. (jlimch) .and. gstart .lt. 0.0)
     * .or. (jende .eq. (jlimch) .and. gend .lt. 0.0)
     * .or. ((jlimch) .lt. jstare)
     * .or. ((jlimch) .lt. jende))
     * jnice = 0
c
c get maximum number of characters to print
      jmaxe = max0(jstare,jende)
      if (jnice .eq. 0 .and. jmaxe .eq. 0 .and. jzeros .eq. 0)
     * jmaxe = jende
      if (jnice .eq. 0 .and. jmaxe .eq. 0 .and. jzeroe .eq. 0)
     * jmaxe = jstare
      if (jnice .eq. 0) then
         jprint = jmaxe-jdxe+1
      elseif (jnice .gt. 0) then
         if (jmaxe .lt. 0) then
            jprint = 1-jdxe
         elseif (jdxe .gt. 0) then
            jprint = jmaxe+1
         else
            jprint = jmaxe-jdxe+2
         endif
      endif
c     add 1 character for leading minus sign
      if (gstart .lt. 0. .or. gend .lt. 0.) jprint = jprint+1
c
c output warning message
      if (jprint .gt. jmaxch) then
         jprint = jmaxch
         print *,'warning: can only print ',jmaxch,
     +           ' characters in tick mark labels'
      endif
c
      if (jnice .ge. 1) go to 2000
c
c get step exponent value
c
c divide the original calculation '10.0**float(jdxe)' by 10 in order
c for the internal write to write out the same exponent as zcetoa.
c
      diff = 1.0/(ccpl1*fscale)
      cdata = gstart-diff
      g1 = (10.0**float(jdxe)) / 10.0d0
c     g1 = 10.0d0**(jdxe-1)
      gexp1 = 1.00001*g1
c
c convert floater [gexp1] to ascii e<jmaxch>.<jmaxch-7> format string
c      kascii(1) = jblaw
c      kascii(2) = jblaw
c      kascii(3) = jblaw
c      call zcetoa (kascii(1), 0, gexp1, jmaxch, jmaxch-7)
      kstring = ' '
      write (kstring,'(E18.11)') gexp1
c
c print step exponent value with no tic or grid marks
      kprint = 4
      kskip = jmaxch-4
      kexp = 1
      cspoin = cspoin-2.0
      call mapdrw (kstring,kprint,kskip,cdata,cspoin,isx,2)
      cspoin = cspoin+2.0
c
c now go through the number of divisions loop
 2000 n = knp+1
      do 4099 i = 1, n
c
c get the ascii value to be plotted
      g1 = -(i-1)*gdx
      cdata =(gstart-g1)*1.000001
c
      if (icomp(gstart,g1) .eq. 0) then
c
c        get (dummy) exponent for zero
         jdatae = min0(jdxe,0)
c
c        get character code for zero
c
c I think that the original code actually generated one more character
c than necessary.
c
c         kascii(1) = 8h 0.00000
c         kascii(2) = 8h0000000e
c         kascii(3) = 8h+00
         kstring = ' 0.00000000000e+00'
      else
c
c        get exponent for data
         jdatae = log10a(cdata)
c
c        get character code for data
c        convert floater [cdata] to ascii e<jmaxch>.<jmaxch-7> format string
c         kascii(1) = jblaw
c         kascii(2) = jblaw
c         kascii(3) = jblaw
c         call zcetoa (kascii(1), 0, cdata, jmaxch, jmaxch-7)
c
c This code isn't clear in some parts what it is doing.  Convert the
c number to the same format that zcetoa uses so that code changes
c aren't necessary.
c
         kstring = ' '
         write (kstring,'(E18.11)') cdata / 10.0d0
         kstring(2:jmaxch-5) = kstring(3:jmaxch-4)
         kstring(2:2) = kstring(3:3)
         kstring(3:3) = '.'
      endif
c
c right justify character code for data
      if (jnice .ge. 1) go to 2600
c
c have mean data
c      call zmovechr (kascii, 2, kascii, 3, jmaxch-7)
c
c This appears to move the fractional portion over the decimal.
c so 1.23456789012e+12 becomes 1234567890122e+12
c
      kstring(3:jmaxch-5) = kstring(4:jmaxch-4)
      if (icomp(cdata,0.) .eq. 0) then
         jprint = 1
      else
         jprint = jdatae-jdxe+1
      endif
      go to 3000
c
c have nice data
 2600 if (jdatae .le. -1) go to 2700
c
c nice data has non-negative exponent
c
c This appears to do the same as the above zmovechr except with 'nice'
c data.  So  1.23456789012e+12 may become 1234567890122e+12
c
c      call zmovechr (kascii, 2, kascii, 3, jdatae)
      kstring(3:2+jdatae) = kstring(4:3+jdatae)
c      call zmovechr (kascii, jdatae+2, ndotz, 0, 1)
      kstring(jdatae+3:jdatae+3) = ndotz(1:1)
      go to 2800
c
c nice data has negative exponent
 2700 jdshl = -jdatae-1
c      call zmovechr (ns, 0, kascii, 0, jmaxch)
      ns(1:jmaxch) = kstring(1:jmaxch)
c      call zmovechr (kascii, jdshl+3, ns, 3, jmaxch-7-jdshl)
      kstring(jdshl+4:jmaxch-7) = ns(4:jmaxch-jdshl-4)
c      call zmovechr (kascii, jdshl+2, ns, 1, 1)
      kstring(jdshl+3:jdshl+3) = ns(2:2)
c      call zmovechr (kascii, 1, ndotz, 0, jdshl+1)
      kstring(2:jdshl+2) = ndotz(1:jdshl+1)
c
c shift data right
 2800 jprinn = jprint-(max0(jdatae,-1)+1)+min0(jdxe,0)
      if (jdxe .le. -1) jprinn = jprinn-1
c      call zmovechr (ns, 0, kascii, 0, jmaxch)
      ns(1:jmaxch) = kstring(1:jmaxch)
c      call zmovechr (kascii, jprinn, ns, 0, jprint)
      kstring(jprinn+1:jprinn+jprint) = ns(1:jprint)
c      call zmovechr (kascii, 0, jblaw, 0, jprinn)
      kstring(1:jprinn) = ' '
c
c get number of significant digits to print for data
 3000 kprint = jprint
      if (cdata .lt. 0.0) kprint = kprint+1
c
c get number of characters to skip for data
      kskip = 1
      if (cdata .lt. 0.0) kskip = 0
c
c indicate no exponent
 3100 kexp = 0
c
c print data
      gsavsp = cspoin
      cspoin = cspoin-float(kprint)+4.
      if (cspoin .lt. 1.0) cspoin = 1.0
      call mapdrw (kstring,kprint,kskip,cdata,cspoin,isx,itic)
c
c restore cspoin value
      cspoin = gsavsp
c
 4099 continue
c
      goto 9999
c
c this is log scaling
 5000 jnp = (2.0*knp)/(fmax-fmin)
      if (jnp .gt. 8) jnp = 8
      kexp = 0
c
c now go through the number of divisions loop
      i1 = fmax-fmin+1
      do 6099 i = 1, i1
c
c write out exponent with line
      expsta = fmin+i-1.0
      cdata = expsta
      data1 = 10.0**cdata
c
c     cft-dependent encode statement has been converted to baselib calls
c     encode (12, 111, kascii) data1
c
c     convert floater [data1] to ascii e12.5 format string
c      call zmovechr (kascii(1), 0, "                ", 0, 16)
c      call zcetoa (kascii(1), 0, data1, 12, 5)
c
c Only the exponent seems to be of interest here.  Divide by 10
c to get the same exponent returned from zcetoa.
c
      kstring = ' '
      write (kstring,'(E12.5)') data1 / 10.0d0
c
      kprint = 4
      kskip = 8
      call mapdrw (kstring,kprint,kskip,cdata,cspoin,isx,itic)
c
      if (jnp .le. 0) go to 6099
c
c specify to print out rightmost character in kstring array with line
      kprint = 1
      kskip = 7
      cspoin = cspoin+3.0
c
      if (i .eq. i1) goto 9999
c
c now go through parts loop
c
c The following call to mapdrw seems to only print one character at
c a time.  The previous internal write does not seem to influence
c this loop.
c
      jnp1 = jnp+1
      do 7099 j = 2, jnp1
      cdata = expsta+alog10(float(j))
c      kascii(1) = 48+j
      kstring(1:8) = '        '
      kstring(8:8) = char (48+j)
      call mapdrw  (kstring,kprint,kskip,cdata,cspoin,isx,itic)
 7099 continue
c
      cspoin = cspoin-3.0
c
 6099 continue
c
c Reset the map to what it was before entering maplab
c
 9999   call tgsetmap (vl,vr,vb,vt,wl,wr,wb,wt,it)
        call gsclip (kepclp)
        return
        END
      subroutine setusr(xc1,xcm,yc1,ycn)
c
c * * set viewport coordinates
c
      call cpsetr('XC1', xc1)
      call cpsetr('XCM', xcm)
      call cpsetr('YC1', yc1)
      call cpsetr('YCN', ycn)
c
      return
      end
      subroutine setvpt(vpl,vpr,vpb,vpt)
c
c  * * set viewport limits
c
      call cpsetr('VPL', vpl)
      call cpsetr('VPR', vpr)
      call cpsetr('VPB', vpb)
      call cpsetr('VPT', vpt)
c
      return
      end
c**************************************************************************
c
c  taxis.F - axis drawing companion for gaxis
c
c  contents	taxis - draws axis
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

      subroutine taxis (fxb, fyb, fww)
c=======================================================================
c 12/28/2001 ler pppl - dimension cliprect
c=======================================================================
c
      real*8 fxb, fyb, fww, dcxb, dcyb, gxp, gxe, gyp, gye
cc created: may17.1976 by mick archuleta
c
cc revised: sep26.1980 by steven williams
c
c this routine can be called after a call to gaxis has been made.
c its purpose is to draw an axis line without tick marks. the
c information for the axis line is still available from the last call
c to gaxis. thus, fxb and fyb are used to specify the begin point
c of this tick line. the end point is calculated by taxis and it
c is in the same relative direction as the line drawn by gaxis. fww is
c a variable which tells taxis how long to make the tick
c marks and in which direction. if the same direction and length, then
c it should be 1, if the opposite direction and same length, then
c it should be -1.
c
cc variable declarations:
c
      logical klr
      common /tvaxis/ cxb, cyb, cxe, cye, ctickx, cticky
     * , coffx, coffy, cac, cbc, kori, ksize, kn1, klr
     * , klxsav, klysav, cpx(30), cpy(30), kdiv, kepclp
c
c xdiff,ydiff   - difference between this and last point
c fxbl,fybl     - linear version of fxb,fyb
c cxbl,cybl     - linear version of cxb,cyb
c
      real xdiff,ydiff
      real fxbl,fybl,cxbl,cybl

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace

      real cliprect(4)



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'taxis'

c
c Convert the starting point from taxis and gaxis to linear and find the
c difference in the x and y directions.
c
      call culxy (fxb,fyb,fxbl,fybl)
      dcxb = cxb
      dcyb = cyb
      call culxy (dcxb,dcyb,cxbl,cybl)
      xdiff = fxbl - cxbl
      ydiff = fybl - cybl
c
c Get world and NDC coordinates as well as map type, then set the
c map type to be linear.
c
      call tggetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,itype)
      call tgsetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,1)
c
c Turn clipping off so tic marks outside of clipping rect will be drawn
c
      call gqclip (ierr,kepclp,cliprect)
      call gsclip (0)
c
c Loop through the array cpx and cpy to place tic marks at each location.
c
      do 80 j = 1, kdiv
         gxp = cpx(j) + xdiff
         gyp = cpy(j) + ydiff
         call line (gxp, gyp, gxp+cac*ctickx*fww, gyp+cbc*cticky*fww)
   80 continue
c
c Reset the map and draw the axis
c
      call tgsetmap (cxmi,cxma,cymi,cyma,cx1,cx2,cy1,cy2,itype)
      call line (fxb, fyb, fxb+cxe-cxb, fyb+cye-cyb)
c
c Return Clipping back to its original setting
c
      call gsclip (kepclp)
c
      if (tgname .eq. 'taxis') tgname = 'NONE'
      return
      END


c**************************************************************************
c
c  tggks.f  -  Routines to initialize the gks drivers
c
c  contents	The routines in this file open GKS if it is not already
c		opened and then open a specific workstation given a
c		workstation identifier and a connection id.  The names
c		and the workstation that is opened will change depending
c		on what workstations are available with your implementation
c		of GKS.  After initializing a workstation, tginit should
c		be called to put the workstation in a state that is
c		consistent with what the library thinks it should be.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

c**************************************************************************
c
c  cgmbin -  CGM binary
c
c  description	Open GKS if necessary, get the identity of the workstation
c		set the output file name, open and then activate the
c		workstation.
c
c**************************************************************************

        subroutine tgchset (ifont,index,iorient)

        real size,upx,upy
        parameter (PI = 3.1415926535897932)

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



        ichindx = index
        ichfont = ifont

c font with character precision

        call gstxfp (ichfont,2)

c Set the rotation based on the angle in iorient

        if (iorient .eq. 0) then
          chrot = 0.
          chupx = cos (chrot + PI/2.)
          chupy = sin (chrot + PI/2.)
	  if (dotrans) then
	    upx = cos (chrot + rotat + PI/2.)
	    upy = sin (chrot + rotat + PI/2.)
	  else
	    upx = chupx
	    upy = chupy
	  endif
          call gschup (upx,upy)
          size = chparm(1,ichindx+1)
          if (dotrans) size = size * scaley
          call gschh (size)
          chaddx = - chparm(2,ichindx+1) * chupx
          chaddy = - chparm(2,ichindx+1) * chupy
          size = chparm(3,ichindx+1)
          if (dotrans) size = size * scalex/scaley
          call gschxp (size)
        else if (iorient .eq. 1) then
          chrot = PI/2.
          chupx = cos (chrot + PI/2.)
          chupy = sin (chrot + PI/2.)
	  if (dotrans) then
	    upx = cos (chrot + rotat + PI/2.)
	    upy = sin (chrot + rotat + PI/2.)
	  else
	    upx = chupx
	    upy = chupy
	  endif
          call gschup (upx,upy)
          size = chparm(1,ichindx+1)
          if (dotrans) size = size * scalex
          call gschh (size)
          chaddx = - chparm(2,ichindx+1) * chupx
          chaddy = - chparm(2,ichindx+1) * chupy
          size = chparm(3,ichindx+1)
          if (dotrans) size = size * scaley/scalex
          call gschxp (size)
        endif

        return
        end

c*************************************************************************
c
c  tgchclip  -  Turn character clipping on or off
c
c  synopsis	call tgchclip (iclip)
c		integer iclip		Clipping flag
c
c  description	All character output is funneled through GTEXT.  When a
c		user calls GTEXT the characters are never clipped, when
c		internal routines need character output, it may be
c		necessary to turn the clipping on temporarily.  TGCHCLIP
c		will turn the clipping on or off depending on whether
c		iclip is 1 or 0 respectively.
c
c*************************************************************************

        subroutine tgchclip (iclflag)

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



        ichclip = iclflag

        return
        end

c*************************************************************************
c
c  setch  - set character attributes for absolute grid
c
c  synopsis     call setch (x,y)
c               call setch (x,y,icase)
c               call setch (x,y,icase,isize)
c               call setch (x,y,icase,isize,iorient)
c               call setch (x,y,icase,isize,iorient,ifont)
c
c               real x,y        Position
c               integer icase   Case for future calls
c               integer isize   Size of characters
c               integer iorient Orientation of strings
c               integer itype   Font to use
c
c  description  x,y will be the position in absolute coordinates.
c
c*************************************************************************

        subroutine setch (x,y,icase,isize,iorient,ifont)

	real*8 x, y

        integer itsize,itorient,itfont
        real x2,x3,x4,y2,y3,y4

        logical tgisfont

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'SETCH'



c Set temporary variables to the old values of the character attributes

        itsize = ichindx
        itorient = ichangle
        itfont = ichfont

c If a font was passed, check to see if the font number passed to this
c routine is valid.  If it is not, then itfont will not be changed.

        if (ifont .ne. -1) then
          if (tgisfont (ifont)) itfont = ifont
        endif

c angle to rotate

        if (iorient .ne. -1) itorient = iorient

c specify sizes from 0 to 3

        if (isize .ne. -1) itsize = isize

c specify which case to use

        if (icase .ne. -1) ichcase = icase

        call tgchset (itfont,itsize,itorient)

        if (itorient .eq. 0) then
          chx = (x-1.0) * chparm(4,ichindx+1)
          chy = (y-1.0) * chparm(2,ichindx+1)
        else
          chx = (x-1.0) * chparm(2,ichindx+1)
          chy = (y-1.0) * chparm(4,ichindx+1)
        endif

c If autofeed is not on, then don't reposition character. Go to end.

        if (.not. autofeed) goto 9999

c Calculate the approximate corners of the rectangle that bounds the initial
c character of the string in clockwise order.  tgchset will have set chupx to
c the cosine of the angle of the up vector and chupy to the sine of that angle.
c The up vector is at an angle of chrot+PI/2.

        xtemp = cos(chrot)
        ytemp = sin(chrot)
        x2 = (chparm(2,itsize+1)*chupx) + chx
        y2 = (chparm(2,itsize+1)*chupy) + chy
        x3 = (chparm(4,itsize+1)*xtemp) + x2
        y3 = (chparm(4,itsize+1)*ytemp) + y2
        x4 = (chparm(4,itsize+1)*xtemp) + chx
        y4 = (chparm(4,itsize+1)*ytemp) + chy

c If the rectangle is outside the NDC area, move it into this area.

        xmax = max (chx,max(x2,max(x3,x4)))
        xmin = min (chx,min(x2,min(x3,x4)))
        ymax = max (chy,max(y2,max(y3,y4)))
        ymin = min (chy,min(y2,min(y3,y4)))

        if (xmin .lt. 0.) chx = (chx - xmin)
        if (xmax .gt. 1.) chx = 1. - (xmax - chx)
        if (ymin .lt. 0.) chy = (chy - ymin)
        if (ymax .gt. 1.) chy = 1. - (ymax - chy)

9999    if (tgname .eq. 'SETCH') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  setlch  - set position, intensity and size of characters
c
c  synopsis     call setlch (x,y)
c               call setlch (x,y,icase)
c               call setlch (x,y,icase,isize)
c               call setlch (x,y,icase,isize,iorient)
c               call setlch (x,y,icase,isize,iorient,ifont)
c
c               real x,y        Position
c               integer icase   Case for future calls
c               integer isize   Size of characters
c               integer iorient Orientation of strings
c               integer itype   Font to use
c
c  description  x,y will be the position in user coordinates for the
c               next string to be plotted.
c
c               intens is the intensity of the characters.
c                 0 - (default) low
c                 1 - high
c
c*************************************************************************

        subroutine setlch (x,y,icase,isize,iorient,ifont)

	real*8 x, y

        integer itsize,itorient,itfont
        real x2,x3,x4,y2,y3,y4

        logical tgisfont

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'SETLCH'



c Set temporary variables to the old values of the character attributes

        itsize = ichindx
        itorient = ichangle
        itfont = ichfont

c If a font was passed, check to see if the font number passed to this
c routine is valid.  If it is not, then itfont will not be changed.

        if (ifont .ne. -1) then
          if (tgisfont (ifont)) itfont = ifont
        endif

c angle to rotate

        if (iorient .ne. -1) itorient = iorient

c specify sizes from 0 to 3

        if (isize .ne. -1) itsize = isize

c specify which case to use

        if (icase .ne. -1) ichcase = icase

c convert x,y to window linear coordinates and then to NDC coordinates

10      call culxy (x,y,chx,chy)
        call clnxy (chx,chy,chx,chy)
        call tgchset (itfont,itsize,itorient)

c If autofeed is not on, then don't reposition character. Go to end.

        if (.not. autofeed) goto 9999

c Calculate the approximate corners of the rectangle that bounds the initial
c character of the string in clockwise order.  tgchset will have set chupx to
c the cosine of the angle of the up vector and chupy to the sine of that angle.
c The up vector is at an angle of chrot+PI/2.

        xtemp = cos(chrot)
        ytemp = sin(chrot)
        x2 = (chparm(2,itsize+1)*chupx) + chx
        y2 = (chparm(2,itsize+1)*chupy) + chy
        x3 = (chparm(4,itsize+1)*xtemp) + x2
        y3 = (chparm(4,itsize+1)*ytemp) + y2
        x4 = (chparm(4,itsize+1)*xtemp) + chx
        y4 = (chparm(4,itsize+1)*ytemp) + chy

c If the rectangle is outside the NDC area, move it into this area.

        xmax = max (chx,max(x2,max(x3,x4)))
        xmin = min (chx,min(x2,min(x3,x4)))
        ymax = max (chy,max(y2,max(y3,y4)))
        ymin = min (chy,min(y2,min(y3,y4)))

        if (xmin .lt. 0.) chx = 0.
        if (xmax .gt. 1.) chx = 1. - (xmax-xmin)
        if (ymin .lt. 0.) chy = 0.
        if (ymax .gt. 1.) chy = 1. - (ymax-ymin)

9999    if (tgname .eq. 'SETLCH') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  gtext  -  plot a string
c
c  synopsis     call gtext (string)
c               call gtext (string,num)
c               call gtext (string,num,ioffset)
c
c  description  Plots string of num characters starting at ioffset.
c
c*************************************************************************

        subroutine gtext (string,num,ioffset)
        character*(*) string
        integer num,ioffset

        integer inum,ioff
        character*256 strbuf
        dimension cliprect(4)
	real x,y

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'GTEXT'



        if (num .gt. 256) then
          call tgerror (7)
          goto 9999
        endif

c Get the parameters into local variables.  The defaults are the length
c of the string and an offset of 0

        if (num .eq. -1) then
	  inum = len (string)
	else
	  inum = num
	endif
	
	if (ioffset .eq. -1) then
	  ioff = 1
	else
	  ioff = ioffset + 1
	endif

c Convert the string to the appropriate case and save in strbuf

        if (ichcase .eq. 0) then
          call tgtoup (string(ioff:ioff+inum-1),strbuf,inum)
        elseif (ichcase .eq. 1) then
          call tgtolow (string(ioff:ioff+inum-1),strbuf,inum)
        else
          strbuf = string(ioff:ioff+inum-1)
        endif

c Reposition the currrent x,y in normalized coordinates if needed

        if (autofeed .and. (chx .gt. 1.0 .or. chx .lt. 0.0 .or.
     +      chy .gt. 1.0 .or. chy .lt. 0.0)) then
          if (chx .gt. 1.0) then
            chx = abs (chaddx)
          else if (chx .lt. 0.0) then
            chx = 1.0 - abs (chaddx)
          endif
          if (chy .gt. 1.0) then
            chy = abs (chaddy)
          else if (chy .lt. 0.0) then
            chy = 1.0 - abs (chaddy)
          endif
          call frame (-1)
        endif

c Turn off clipping so the labels will appear if character clipping is off

        if (ichclip .eq. 0) then
          call gqclip (ierr,kepclp,cliprect)
          call gsclip (0)
        endif

c Write string and update the current x,y in normalized coordinates

        if (dotrans) then
          call cntxy (chx,chy,x,y)
          call gtx (x,y,strbuf(1:inum))
        else
          call gtx (chx,chy,strbuf(1:inum))
        endif

c Return clipping to its original state if character clipping is off

        if (ichclip .eq. 0) call gsclip (kepclp)

        chx = chx + chaddx
        chy = chy + chaddy

9999    if (tgname .eq. 'GTEXT') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  wrtstr  -  write an array of strings
c
c  synopsis     call wrtstr (string,num)
c               integer num                      Number of strings
c               character*(*) string(num)        Array of strings
c
c*************************************************************************

        subroutine wrtstr (string,num)
        character*(*) string(*)
        integer num
        integer length

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'WRTSTR'



        length = len(string(1))

        do 10 i = 1,num
          call gtext (string(i),length,0)
10      continue

9999    if (tgname .eq. 'WRTSTR') tgname = 'NONE'
        return

        end
c*************************************************************************
c
c  tginit.F  -  initialization and necessary routines in TV80GKS
c
c  contents     plote -  end plotting
c		plotea -  flush buffers
c               endpl -  flush buffer and reset mapping
c               frame -  Flush buffers and frame advance
c               colori -  Set the current color by index
c               colora -  Set the current color by name
c               tginit -  Initialize the interface
c		tgreset -  Reinitialize the interface
c               tgdefs - block data for tv80 to gks shell
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c*************************************************************************

c*************************************************************************
c
c  plote  -  end plotting
c
c  synopsis     call plote
c
c  description  Flush the line buffer, deactivate all workstations that
c               are currently active, close all that are open and then
c               close gks.
c
c*************************************************************************

        subroutine plote

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'PLOTE'



        call tgflush

c While GKS is in the state of workstation active (3), get the workstation
c identifier and then deactivate the workstation.  When all workstations
c are deactivated, the state will be workstation open (2).

10	call gqops (istate)
	if (istate .eq. 3) then
	  call gqacwk (1,ierr,numws,iws)
	  call gdawk (iws)
	  goto 10
	endif

c While GKS is in the state of workstation open (2), get the workstation
c identifier and then close the workstation.  When all workstations
c are closed, the state will be workstation gks open (1).

20	call gqops (istate)
	if (istate .eq. 2) then
	  call gqopwk (1,ierr,numws,iws)
	  call gclwk (iws)
	  goto 20
	endif

        call gclks

        doinit = .TRUE.

        if (tgname .eq. 'PLOTE') tgname = 'NONE'
        return

c Some sort of disaster happened. Do an emergency close and return.

999     call geclks

c Reset init variable so if the user wants to use the library without
c restarting his program, the library will initialize again.

        doinit = .TRUE.

        if (tgname .eq. 'PLOTE') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  plotea  -  flush buffers
c
c  synopsis     call plotea ()
c
c  description  Flushes TV80GKS buffers.  Under TV80LIB, this routine had
c		a single integer argument which was called 'iempty'. The
c		argument was not used so it was removed.
c
c*************************************************************************

        subroutine plotea ()

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'PLOTEA'



        call tgflush

	call gqacwk (0,ierr,numws,iws)
        if (ierr .eq. 0) then
	  do 10 i=1,numws
	    call gqacwk (i,ierr,ierr,iws)
	    call guwk (iws,1)
10	  continue
	endif

        if (tgname .eq. 'PLOTEA') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  endpl  -  flush buffer and reset mapping
c
c  synopsis     call endpl
c
c  description  Flush the spps buffers and reset mapping.
c
c*************************************************************************

        subroutine endpl

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'ENDPL'



        call tgflush
        call map (0.d0,1.d0,0.d0,1.d0, 0.d0,1.d0,0.d0,1.d0)

        if (tgname .eq. 'ENDPL') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  frame  -  Flush buffers and frame advance
c
c  synopsis     call frame ()
c               call frame (nframe)
c               integer nframe          Number of empty frames
c
c  description  Buffers are flushed and empty frames are inserted if
c		requested.  max(nframe,1) frame advances will be done.
c
c*************************************************************************

        subroutine frame (nframe)

        integer nframe
        integer iempty

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'FRAME'



        call tgflush

c if nframe was passed to this routine and it was positive, then the
c number of frames will be inserted.

        if (nframe .ge. 1) then
          iempty = nframe
        else
          iempty = 1
        endif

c Create the empty frames

        call gqacwk (0,ierr,numwrk,iwkid)
        do 10 i = 1,numwrk
          call gqacwk (i,ierr,idummy,iwkid)
          do 20 j = 1,iempty
20          call gclrwk (iwkid,1)
10      continue

        if (tgname .eq. 'FRAME') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  colori  -  Set the current color by index
c
c
c*************************************************************************

        subroutine colori (index)

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'COLORI'



c search for the index that represents the color.

        if (index .ge. 0) then
	  icurclr = index
        else
          call tgerror (1)
          if (tgname .eq. 'COLORI') tgname = 'NONE'
          return
        endif

        call tgflush

        call gsplci (icurclr)
        call gspmci (icurclr)
        call gstxci (icurclr)
        call gsfaci (icurclr)

        if (tgname .eq. 'COLORI') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  colora  -  Set the current color by name
c
c
c*************************************************************************

        subroutine colora (colorstr)

	character*(*) colorstr

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'COLORA'



c search for the string that represents the color.

        if (colorstr(1:5) .eq. 'black' .or.
     +      colorstr(1:5) .eq. 'BLACK') then
          icurclr = 0
        else if (colorstr(1:5) .eq. 'white' .or.
     +      colorstr(1:5) .eq. 'WHITE') then
          icurclr = 1
        else if (colorstr(1:3) .eq. 'red' .or.
     +      colorstr(1:3) .eq. 'RED') then
          icurclr = 2
        else if (colorstr(1:5) .eq. 'green' .or.
     +      colorstr(1:5) .eq. 'GREEN') then
          icurclr = 3
        else if (colorstr(1:4) .eq. 'blue' .or.
     +      colorstr(1:4) .eq. 'BLUE') then
          icurclr = 4
        else if (colorstr(1:4) .eq. 'cyan' .or.
     +      colorstr(1:4) .eq. 'CYAN') then
          icurclr = 5
        else if (colorstr(1:7) .eq. 'magenta' .or.
     +      colorstr(1:7) .eq. 'MAGENTA') then
          icurclr = 6
        else if (colorstr(1:6) .eq. 'yellow' .or.
     +      colorstr(1:6) .eq. 'YELLOW') then
          icurclr = 7
        else
          call tgerror (1)
          if (tgname .eq. 'COLORA') tgname = 'NONE'
          return
        endif

        call tgflush

        call gsplci (icurclr)
        call gspmci (icurclr)
        call gstxci (icurclr)
        call gsfaci (icurclr)

        if (tgname .eq. 'COLORA') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  tginit  -  Initialize the interface
c
c  synopsis     call tginit (iwkid)
c		integer iwkid		workstation identifier
c
c  description  Initialize the tv80 to gks shell so that gks will be in
c               a state similar to the default state of tv80lib.  This
c               must be called after a every workstation is opened.
c		iwkid is the workstation identifier associated with the
c		workstation.
c
c*************************************************************************

        subroutine tginit (iwkid)

c=======================================================================
c 12/28/2001 pppl ler - change 1 to variable i_one in call to gqdsp 
c                       ... this change done earlier than 12/28/2001
c=======================================================================
c The tgdefs BLOCK DATA subprogram must be mentioned in a called routine
c in order for it to be included into the program.  tginit must be
c called for the shell to work properly, so it is placed here.

        external tgdefs

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c include the site dependent info

c***************************************************************************
c
c  tgsite.inc - site dependent information
c
c***************************************************************************

c The minimum and maximum font numbers

	integer minfnt,maxfnt
	parameter (minfnt = 1, maxfnt = 20)

c The default font

	integer deffnt
	parameter (deffnt = 1)

c The character height and expansion factor for fonts 1, 2, 3 and 4

	parameter (chrht1 = .009, chexp1 = 0.87)
	parameter (chrht2 = .013, chexp2 = 0.90)
	parameter (chrht3 = .016, chexp3 = 0.98)
	parameter (chrht4 = .018, chexp4 = 1.28)


c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'TGINIT'



c If the this is the first time the initialization routine was called,
c set up the defaults of the library.

        if (doinit) then

          doinit = .FALSE.
          dotrans = .FALSE.

c set global variables

          ichfont = deffnt
          numpnt = 0
          autofeed = .TRUE.

c set the default mapping

          call map (0.d0,1.d0,0.d0,1.d0,0.d0,1.d0,0.d0,1.d0)

c tv80lib defaults to no clipping, turn it off

          call dders (1)

c set line attributes

          call setpch (0,1,ichfont,100)

c set the colormap and current color for the device

          call colora ('WHITE')

c set default font, tv80 font size index 0, rotation 0, clipping off

          call setch (0.0D0,0.0D0,0,0,0,ichfont)
          call tgchclip (0)

c set text path to left to right and allignment to bottom left

          call gstxp (0)
          call gstxal (1,5)

c set transformation to none

	  call init2d ()

        endif

c Alter the workstation to conform to tv80 by centering the window

        call gqwkc (iwkid,ierr,dummy,itype)
        i_one = 1
        call gqdsp (itype,ierr,i_one,rx,ry,lx,ly)

        if (rx .gt. ry) then
          diff = real(rx - ry)/2.
          call gswkvp (iwkid,diff,real(ry)+diff,0.,real(ry))
        else if (rx .lt. ry) then
          diff = real(ry - rx)/2.
          call gswkvp (iwkid,0.,real(rx),diff,real(rx)+diff)
        endif

c set up to 16 colors

        call gqeci (iwkid,0,ierr,numclr,dummy)
        if (numclr .eq. 0) numclr = 16
        do 10 i = 1,min(numclr,16)
          call gscr (iwkid,i-1,red(i),green(i),blue(i))
10      continue

        if (tgname .eq. 'TGINIT') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  tgreset  -  Reinitialize the interface
c
c  synopsis     call tgreset ()
c
c  description  This routine allows the user to muck with gks and then
c		call this routine to reset gks to what the library thinks
c		it should be.
c
c*************************************************************************

        subroutine tgreset ()

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'TGRESET'



c Make sure that gks was opened (operation state .gt. 0)

        call gqops (iopstate)
        if (iopstate .eq. 0) then
          call tgerror (20)
          call exit (1)
        endif

c Make sure the user has called tginit to initialize the library

        if (doinit) then
          call tgerror (21)
          call exit (1)
        endif

c Set the viewport, window and current transform

        call gsvp (1,xvmin,xvmax,yvmin,yvmax)
        call gswn (1,xvmin,xvmax,yvmin,yvmax)
        call gselnt (1)

c Set the clipping

        call gsclip (iclip)

c Set font number, character precision to stroke, up vector, character height,
c and character expansion factor.

        call gstxfp (ichfont,2)
        call gschup (chupx,chupy)
        call gschh (chparm(1,ichindx+1))
        call gschxp (chparm(3,ichindx+1))

c Set text path to right and allignment to left bottom

        call gstxp (0)
        call gstxal (1,5)

c Set color

        call gsplci (icurclr)
        call gspmci (icurclr)
        call gstxci (icurclr)
        call gsfaci (icurclr)

        if (tgname .eq. 'TGRESET') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  tgdefs - block data for tv80 to gks shell
c
c  description  Sets the common blocks to something meaningful
c
c*************************************************************************

        block data tgdefs

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c include the site dependent info

c***************************************************************************
c
c  tgsite.inc - site dependent information
c
c***************************************************************************

c The minimum and maximum font numbers

	integer minfnt,maxfnt
	parameter (minfnt = 1, maxfnt = 20)

c The default font

	integer deffnt
	parameter (deffnt = 1)

c The character height and expansion factor for fonts 1, 2, 3 and 4

	parameter (chrht1 = .009, chexp1 = 0.87)
	parameter (chrht2 = .013, chexp2 = 0.90)
	parameter (chrht3 = .016, chexp3 = 0.98)
	parameter (chrht4 = .018, chexp4 = 1.28)


        data doinit /.TRUE./

        data tgname/"NONE    "/
        data minfont/minfnt/
        data maxfont/maxfnt/

        data chparm/chrht1,.0156250,chexp1,.00781250,
     +              chrht2,.0235294,chexp2,.0117647,
     +              chrht3,.0312500,chexp3,.0156250,
     +              chrht4,.0476190,chexp4,.0238095/
	
        data red   /0.,1.,1.,0.,0.,0.,1.,1.,.0,.5,.5,.0,.0,.0,.5,.5/
        data green /0.,1.,0.,1.,0.,1.,0.,1.,.0,.5,.0,.5,.0,.5,.0,.5/
        data blue  /0.,1.,0.,0.,1.,1.,1.,0.,.0,.5,.0,.0,.5,.5,.5,.0/

        end

c*************************************************************************
c
c  tgline.F  -  Line generation routines
c
c  contents     setcrt  -  Set the current point
c		vector  -  Draw from current point
c		tgflush  -  Flush the line buffer
c		line  -  Draw a line
c		linep  -  Draw a line with points
c		point  -  Draw a point
c		points  -  Draw a series of points
c		pointc  -  Draw a series of lines
c		trace  -  Draw a series of lines
c		tracep  -  Draw a series of dotted lines.
c		tracec  -  Draw a series of lines
c		setpch  -  Set parameters for pointc and tracec.
c		plotv  -  Draw an arrow
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c*************************************************************************

c*************************************************************************
c
c  setcrt  -  Set the current point
c
c  synopsis     call setcrt (x,y)
c
c               real x,y                Point
c
c  description  Sets current point.  If the user has specified a
c               transformation then transform the point before moving.
c
c*************************************************************************

        subroutine setcrt (x,y)

	real*8 x,y

cccc        real x,y

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'SETCRT'



        if (numpnt .gt. 0) call tgflush
        call cufxy (x,y,xpnt(1),ypnt(1))
        numpnt = 1

        if (tgname .eq. 'SETCRT') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  vector  -  Draw from current point
c
c  synopsis     call vector (x,y)
c
c               real x,y                Point
c
c  description  Draw vector.  If the user has specified a transformation,
c               then convert the point before drawing a vector to it.
c
c*************************************************************************

        subroutine vector (x,y)

	real*8 x,y

c        real x,y

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'VECTOR'



        numpnt = numpnt + 1
        call cufxy (x,y,xpnt(numpnt),ypnt(numpnt))
        if (numpnt .eq. MAXPNT) then
          call tgflush
          numpnt = 1
          xpnt(1) = xpnt(MAXPNT)
          ypnt(1) = ypnt(MAXPNT)
        endif

        if (tgname .eq. 'VECTOR') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  tgflush  -  Flush the line buffer
c
c  synopsis     call tgflush
c
c  description  Creates a polyline from the line buffer.  This routine
c               was designed to be called only from tv80lib.  To do a
c               flush, the user should call plotea.
c
c*************************************************************************

        subroutine tgflush ()

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



        if (numpnt .gt. 1) call gpl (numpnt,xpnt,ypnt)
        numpnt = 0

        return
        end

c*************************************************************************
c
c  line  -  Draw a line
c
c  synopsis     call line (x1,y1, x2,y2)
c
c               real x1,y1              First point
c               real x2,y2              Second
c
c  description  Draws a line.
c
c*************************************************************************

        subroutine line (x1,y1, x2,y2)
	real*8 x1,y1,x2,y2


c        real x1,y1,x2,y2

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'LINE'



        call setcrt (x1,y1)
        call vector (x2,y2)

        if (tgname .eq. 'LINE') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  linep  -  Draw a line with points
c
c  synopsis     call linep (x1,y1, x2,y2)
c               call linep (x1,y1, x2,y2, iexponent)
c
c               real x1,y1              First point
c               real x2,y2              Second
c               integer iexponent       Exponent to use for screen res
c
c  description  Draws a line from x1,y1 to x2,y2 using a series of points.
c               The points are assumed to be a distance of 2**iexponent
c               apart.  Where the screen is assumed to be 2**10 x 2**10
c               (1024 x 1024).  Initially iexponent is assumed to be 2.
c
c*************************************************************************

        subroutine linep (x1,y1, x2,y2, iexponent)

 	real*8 x1, y1, x2, y2
        real x1p,y1p,x2p,y2p
	real*8 dx1p, dy1p, dx2p, dy2p, ddeltax, ddeltay
        real vl,vr,vb,vt,wl,wr,wb,wt
        real distx,disty,dist
        integer it

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'LINEP'



c Convert to plotter coordinates

        call culxy (x1,y1,x1p,y1p)
        call clnxy (x1p,y1p,x1p,y1p)
        call culxy (x2,y2,x2p,y2p)
        call clnxy (x2p,y2p,x2p,y2p)

c Save the old mapping and reset the mapping to normalized coordinates

        call tggetmap (vl,vr,vb,vt,wl,wr,wb,wt,it)
        call tgsetmap (vl,vr,vb,vt,vl,vr,vb,vt,1)

c find the distance between successive points by calculating distance
c between endpoints and dividing by 2**k.

        distx = x2p - x1p
        disty = y2p - y1p
        dist = sqrt (distx*distx + disty*disty)
	if (iexponent .gt. 0) iptspace = iexponent
        k = iptspace
        numpts = nint (dist / (2.**k / 1024.))
        if (numpts .gt. 0) then
          deltax = distx / real(numpts)
          deltay = disty / real(numpts)
        endif

c Draw the line of points.  To avoid any round off error, the last point
c is explicitly plotted rather than calculated from a series of additions
c by the points routine.

        call point4 (x2p,y2p)
        if (numpts .gt. 0)
     +    call points4 (x1p,y1p,numpts,0,0,deltax,deltay)

c Reset the mapping to the users coordinates

        call tgsetmap (vl,vr,vb,vt,wl,wr,wb,wt,it)

        if (tgname .eq. 'LINEP') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  point  -  Draw a point
c
c  synopsis     call point (x,y)
c
c               real x,y                Point
c
c  description  Draws a point in the user coordinate system.
c
c*************************************************************************

        subroutine point (x,y)
	real*8 x, y


c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'POINT'



        call setcrt (x,y)
        call vector (x,y)

        if (tgname .eq. 'POINT') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  points  -  Draw a series of points
c
c  synopsis     call points (x,y,num)
c               call points (x,y,num,incx)
c               call points (x,y,num,incx,incy)
c               call points (x,y,num,incx,incy,  0.0)
c               call points (x,y,num,incx,incy,  0.0, 0.0)
c               call points (x,y,num,0   ,incy, delx)
c               call points (x,y,num,0   ,0   , delx,dely)
c               call points (x,y,num,incx,0   ,  0.0,dely)
c
c               real x,y                Points
c               integer num             Number of points
c               integer incx,incy       Increment between successive x
c                                       and y storage locations
c               real delx,dely          The amount to be added to the
c                                       first element in the array when
c                                       the increment is zero, in order
c                                       to generate succesive locations
c
c  description  Draws a series of points.  The default value for incx and
c               incy is 1.  If incx or incy is 0, then we generate the
c               succesive points by adding delx or dely to the previous
c               points.
c
c*************************************************************************

        subroutine points (x,y,num,incx,incy,delx,dely)

c        real xpos,ypos
	real*8 xpos,ypos
        integer index,iaddx,iaddy
c        real deltx,delty,x(1),y(1)
	real*8 delx,dely,x(1),y(1)

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'POINTS'



c Set the value of iaddx and deltx from arguments

        if (incx .eq. -1) then
	  iaddx = 1
	  deltx = 0.0
	else if (incx .eq. 0) then
	  iaddx = 0
	  deltx = delx
	else
	  iaddx = incx
	  deltx = 0.0
	endif
	
c Set the value of iaddy and delty from arguments

        if (incy .eq. -1) then
	  iaddy = 1
	  delty = 0.0
	else if (incy .eq. 0) then
	  iaddy = 0
	  delty = dely
	else
	  iaddy = incy
	  delty = 0.0
	endif
	
c If 3 argument version used, we loop through the points in the array

        if (iaddx .eq. 1 .and. iaddy .eq. 1) then

          do 10 index = 1,num
            call point (x(index),y(index))
10        continue

        else

c If more than 3 arguments were passed, life gets a bit more complicated, the
c user has control of the amount to add to the index and the difference to
c add to the points to generate the new points.


          do 20 index = 0,num-1
            xpos = x(1+index*iaddx)+deltx*real(index)
            ypos = y(1+index*iaddy)+delty*real(index)
            call point (xpos,ypos)
20        continue

        endif

9999    if (tgname .eq. 'POINTS') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  pointc  -  Draw a series of lines
c
c  synopsis     call pointc (ch,x,y,num)
c               call pointc (ch,x,y,num,incx)
c               call pointc (ch,x,y,num,incx,incy)
c               call pointc (ch,x,y,num,incx,incy,  0.0)
c               call pointc (ch,x,y,num,incx,incy,  0.0, 0.0)
c               call pointc (ch,x,y,num,0   ,incy, delx)
c               call pointc (ch,x,y,num,0   ,0   , delx,dely)
c               call pointc (ch,x,y,num,incx,0   ,  0.0,dely)
c
c               character ch            Character to plot at each point
c               real x,y                Points
c               integer num             Number of points
c               integer incx,incy       Increment between successive x
c                                       and y storage locations
c               real delx,dely          The amount to be added to the
c                                       first element in the array when
c                                       the increment is zero, in order
c                                       to generate succesive locations
c
c  description  Draws a series of points.  The default value for incx and
c               incy is 1.  If incx or incy is 0, then we generate the
c               succesive endpoint by adding delx or dely to the previous
c               endpoint.  The character ch will be place at each point
c               if it is more than ilnspace points farther away from the
c               last point where the character was plotted.
c
c*************************************************************************

        subroutine pointc (ch,x,y,num,incx,incy,delx,dely)

        character*1 ch
        integer index,iaddx,iaddy
c       real xpos,ypos,xlin,ylin
        real xlin,ylin
        real xlastc,ylastc,xnextc,ynextc,xdiffc,ydiffc
	real*8 delx,dely,x(1),y(1)
	real*8 xpos,ypos

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'POINTC'



c allignment is centered in horizontal and vertical directions, character
c clipping is turned on if clipping was turned on

        call gqtxal (ierr,ioldhz,ioldvt)
        call gstxal (2,3)
        if (iclip .eq. 1) call tgchclip (1)
        autofeed = .FALSE.

c calculate constants for conversion of linear user coordinates to 1024 x 1024

        pxmin = xvmin * 1024.
        pxmax = xvmax * 1024.
        pymin = yvmin * 1024.
        pymax = yvmax * 1024.
        xmult = (pxmax-pxmin)/(xwmax-xwmin)
        xadd = pxmin - xwmin * xmult
        ymult = (pymax-pymin)/(ywmax-ywmin)
        yadd = pymin - ywmin * ymult
        isqrsp = ilnspace*ilnspace

c Set the value of iaddx and deltx from arguments

        if (incx .eq. -1) then
	  iaddx = 1
	  deltx = 0.0
	else if (incx .eq. 0) then
	  iaddx = 0
	  deltx = delx
	else
	  iaddx = incx
	  deltx = 0.0
	endif
	
c Set the value of iaddy and delty from arguments

        if (incy .eq. -1) then
	  iaddy = 1
	  delty = 0.0
	else if (incy .eq. 0) then
	  iaddy = 0
	  delty = dely
	else
	  iaddy = incy
	  delty = 0.0
	endif
	
        call point (x(1),y(1))
        call culxy (x(1),y(1),xlin,ylin)
        xlastc = xlin*xmult + xadd
        ylastc = ylin*ymult + yadd
        call setlch (x(1),y(1),ilncase,ilnindx,0,ilnfont)
        call gtext (ch,1,0)

        do 10 index = 1,num-1
          xpos = x(1+index*iaddx)+deltx*real(index)
          ypos = y(1+index*iaddy)+delty*real(index)
          call point (xpos,ypos)
          call culxy (xpos,ypos,xlin,ylin)
          xnextc = xlin*xmult + xadd
          ynextc = ylin*ymult + yadd
          xdiffc = xnextc - xlastc
          ydiffc = ynextc - ylastc
          if ((xdiffc*xdiffc + ydiffc*ydiffc) .ge. isqrsp) then
            call setlch (xpos,ypos,ilncase,ilnindx,0,ilnfont)
            call gtext (ch,1,0)
            xlastc = xnextc
            ylastc = ynextc
          endif
10      continue

c Reset the text horizontal, vertical alignment and character clipping

        call gstxal (ioldhz,ioldvt)
        call tgchclip (0)
        autofeed = .TRUE.

        if (tgname .eq. 'POINTC') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  trace  -  Draw a series of lines
c
c  synopsis     call trace (x,y,num)
c               call trace (x,y,num,incx)
c               call trace (x,y,num,incx,incy)
c               call trace (x,y,num,incx,incy,  0.0)
c               call trace (x,y,num,incx,incy,  0.0, 0.0)
c               call trace (x,y,num,0   ,incy, delx)
c               call trace (x,y,num,0   ,0   , delx,dely)
c               call trace (x,y,num,incx,0   ,  0.0,dely)
c
c               real x,y                Endpoint of the lines
c               integer num             Number of endpoints
c               integer incx,incy       Increment between successive x
c                                       and y storage locations
c               real delx,dely          The amount to be added to the
c                                       first element in the array when
c                                       the increment is zero, in order
c                                       to generate succesive locations
c
c  description  Draws a series of vectors.  The default value for incx and
c               incy is 1.  If incx or incy is 0, then we generate the
c               succesive endpoint by adding delx or dely to the previous
c               endpoint.
c
c*************************************************************************

        subroutine trace (x,y,num,incx,incy,delx,dely)

        integer index,iaddx,iaddy
c        real deltx,delty,x(1),y(1)
c        real xpos,ypos
	real*8 deltx,delty,x(1),y(1)
	real*8 xpos,ypos

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'TRACE'



c Set the value of iaddx and deltx from arguments

        if (incx .eq. -1) then
	  iaddx = 1
	  deltx = 0.0
	else if (incx .eq. 0) then
	  iaddx = 0
	  deltx = delx
	else
	  iaddx = incx
	  deltx = 0.0
	endif
	
c Set the value of iaddy and delty from arguments

        if (incy .eq. -1) then
	  iaddy = 1
	  delty = 0.0
	else if (incy .eq. 0) then
	  iaddy = 0
	  delty = dely
	else
	  iaddy = incy
	  delty = 0.0
	endif

c If 3 argument version used, we  loop through the points in the array

        if (iaddx .eq. 1 .and. iaddy .eq. 1) then

          call setcrt (x(1),y(1))
          do 10 index = 2,num
            call vector (x(index),y(index))
10        continue

        else

c If more than 3 arguments were passed, life gets a bit more complicated, the
c user has control of the amount to add to the index and the difference to
c add to the points to generate the new points.


          call setcrt (x(1),y(1))
          do 20 index = 1,num-1
            xpos = x(1+index*iaddx)+deltx*real(index)
            ypos = y(1+index*iaddy)+delty*real(index)
            call vector (xpos,ypos)
20        continue

        endif

9999    if (tgname .eq. 'TRACE') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  tracep  -  Draw a series of dotted lines.
c
c  synopsis     call tracep (xarray,yarray,num)
c               call tracep (xarray,yarray,num,k)
c               call tracep (xarray,yarray,num,k,incx)
c               call tracep (xarray,yarray,num,k,incx,incy)
c
c               real xarray(num)        X coordinates
c               real yarray(num)        Y coordinates
c               integer num             Number of coordinates
c               integer k               Exponent to pass to linep
c               integer incx            Increment between x coordinates
c               integer incy            Increment between y coordinates
c
c  description  Plots a series of points between the points given in
c               xarray and yarray.  The points are drawn by making a
c               call to linep.
c
c*************************************************************************

        subroutine tracep (xarray,yarray,num,k,incx,incy)
c       real xarray(1)
c        real yarray(1)
	real*8 xarray(1)
	real*8 yarray(1)
        integer num
        integer k
        integer incx
        integer incy
c
c The value of k needs to be stored between successive calls to tracep.
c The initial value of k is 2.
c
        integer ix,iy

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'TRACEP'



        ix = 1
        iy = 1

        if (incx .ne. -1) ix = incx
	if (incy .ne. -1) iy = incy
	if (k .ne. -1) iptspace = k

        ixindx = 1
        iyindx = 1
        do 10 i = 1,num-1
          call linep (xarray(ixindx),yarray(iyindx),
     +      xarray(ixindx+ix),yarray(iyindx+iy),iptspace)
          ixindx = ixindx + ix
          iyindx = iyindx + iy
10      continue

        if (tgname .eq. 'TRACEP') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  tracec  -  Draw a series of lines
c
c  synopsis     call tracec (ch,x,y,num)
c               call tracec (ch,x,y,num,incx)
c               call tracec (ch,x,y,num,incx,incy)
c               call tracec (ch,x,y,num,incx,incy,  0.0)
c               call tracec (ch,x,y,num,incx,incy,  0.0, 0.0)
c               call tracec (ch,x,y,num,0   ,incy, delx)
c               call tracec (ch,x,y,num,0   ,0   , delx,dely)
c               call tracec (ch,x,y,num,incx,0   ,  0.0,dely)
c
c               character ch            Character to plot at each point
c               real x,y                Endpoint of the lines
c               integer num             Number of endpoints
c               integer incx,incy       Increment between successive x
c                                       and y storage locations
c               real delx,dely          The amount to be added to the
c                                       first element in the array when
c                                       the increment is zero, in order
c                                       to generate succesive locations
c
c  description  Draws a series of vectors.  The default value for incx and
c               incy is 1.  If incx or incy is 0, then we generate the
c               succesive endpoint by adding delx or dely to the previous
c               endpoint.  The character ch will be place at each point
c               if it is more than ilnspace points farther away from the
c               last point where the character was plotted.
c
c*************************************************************************

        subroutine tracec (ch,x,y,num,incx,incy,delx,dely)

        character*1 ch
        integer index,iaddx,iaddy
c        real deltx,delty,x(1),y(1)
c        real xpos,ypos,xlin,ylin
	real*8 delx,dely,x(1),y(1)
	real*8 xpos,ypos
	real xlin,ylin
        real xlastc,ylastc,xnextc,ynextc,xdiffc,ydiffc

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'TRACEC'



c allignment is centered in horizontal and vertical directions, character
c clipping is turned on if clipping was turned on.

        call gqtxal (ierr,ioldhz,ioldvt)
        call gstxal (2,3)
        if (iclip .eq. 1) call tgchclip (1)
        autofeed = .FALSE.

c calculate constants for conversion of linear user coordinates to 1024 x 1024

        pxmin = xvmin * 1024.
        pxmax = xvmax * 1024.
        pymin = yvmin * 1024.
        pymax = yvmax * 1024.
        xmult = (pxmax-pxmin)/(xwmax-xwmin)
        xadd = pxmin - xwmin * xmult
        ymult = (pymax-pymin)/(ywmax-ywmin)
        yadd = pymin - ywmin * ymult
        isqrsp = ilnspace*ilnspace

c Set the value of iaddx and deltx from arguments

        if (incx .eq. -1) then
	  iaddx = 1
	  deltx = 0.0
	else if (incx .eq. 0) then
	  iaddx = 0
	  deltx = delx
	else
	  iaddx = incx
	  deltx = 0.0
	endif
	
c Set the value of iaddy and delty from arguments

        if (incy .eq. -1) then
	  iaddy = 1
	  delty = 0.0
	else if (incy .eq. 0) then
	  iaddy = 0
	  delty = dely
	else
	  iaddy = incy
	  delty = 0.0
	endif

c plot the first point.  A character is plotted at this point if possible.

        call setcrt (x(1),y(1))
        call culxy (x(1),y(1),xlin,ylin)
        xlastc = xlin*xmult + xadd
        ylastc = ylin*ymult + yadd
        if (xlastc .ge. 0. .and. xlastc .lt. 1024. .and.
     +      ylastc .ge. 0. .and. ylastc .lt. 1024.) then
          call setlch (x(1),y(1),ilncase,ilnindx,0,ilnfont)
          call gtext (ch,1,0)
        endif

        do 10 index = 1,num-1
          xpos = x(1+index*iaddx)+deltx*real(index)
          ypos = y(1+index*iaddy)+delty*real(index)
          call vector (xpos,ypos)
          call culxy (xpos,ypos,xlin,ylin)
          xnextc = xlin*xmult + xadd
          ynextc = ylin*ymult + yadd
          xdiffc = xnextc - xlastc
          ydiffc = ynextc - ylastc
          if ((xdiffc*xdiffc + ydiffc*ydiffc) .ge. isqrsp) then
            if (xnextc .ge. 0. .and. xnextc .lt. 1024. .and.
     +          ynextc .ge. 0. .and. ynextc .lt. 1024.) then
              call setlch (xpos,ypos,ilncase,ilnindx,0,ilnfont)
              call gtext (ch,1,0)
            endif
            xlastc = xnextc
            ylastc = ynextc
          endif
10      continue

c Reset the text horizontal, vertical alignment and character clipping.

        call gstxal (ioldhz,ioldvt)
        call tgchclip (0)
        autofeed = .TRUE.

        if (tgname .eq. 'TRACEC') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  setpch  -  Set parameters for pointc and tracec.
c
c  synopsis     call setpch (icase)
c               call setpch (icase,isize)
c               call setpch (icase,isize,itype)
c               call setpch (icase,isize,itype,kspace)
c
c
c  description  Sets the parameters depended upon by pointc and tracec.
c
c*************************************************************************

        subroutine setpch (icase,isize,itype,kspace)

        logical tgisfont

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'SETPCH'



        if (kspace .ne. -1) ilnspace = kspace
	if (itype .ne. -1) then
          if (tgisfont (itype)) ilnfont = itype
        endif
        if (isize .ne. -1) ilnindx = isize
	if (icase .ne. -1) ilncase = icase

        if (tgname .eq. 'SETPCH') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  plotv  -  Draw an arrow
c
c  synopsis     call plotv (x1,y1,x2,y2)
c               call plotv (x1,y1,x2,y2,fact)
c               real x1,y1      Start of arrow
c               real x2,y2      Head of arrow
c               real fact       Some sort of factor to indicate head size
c
c  description  Creates an arrow pointing from x1,y1 to x2,y2.  fact
c               specifies how large the arrowhead will be.  fact should
c               be between 0. and 1.  This routine was extracted from
c               tv80lib and rewritten.
c
c*************************************************************************

	subroutine plotv (xx1, yy1, xx2, yy2, fact)
c=======================================================================
c 12/28/2001 ler pppl - change 'call line' and 'call vector' arguments
c                        from single precision to double precision
c=======================================================================
c
	real*8 xx1, yy1, xx2, yy2, fact, dgx1, dgx2, dgy1, dgy2
c
cc purpose: plot a vector (pat crowleys method - may, 1968)
c           (argument h, if given, is fixed length of arrowhead -
c            otherwise range)
c           (length of arrowhead proportional to length of v)
c
c  Default legth of arrowhead

        parameter (con=0.12)

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'PLOTV'



c Convert the coordinates into normalized coordinates.  Then convert the map
c into a normalized linear map.

        call culxy (xx1,yy1,gx1,gy1)
        call clnxy (gx1,gy1,gx1,gy1)
        call culxy (xx2,yy2,gx2,gy2)
        call clnxy (gx2,gy2,gx2,gy2)
        call tggetmap (vl,vr,vb,vt,wl,wr,wb,wt,it)
        call tgsetmap (vl,vr,vb,vt,vl,vr,vb,vt,1)
c
        dx = gx2-gx1
        dy = gy2-gy1
        pro = con
        if (fact .ne. -1.0) then
          pro = dx*dx+dy*dy
          if (pro .gt. 0) pro = fact/sqrt(pro+pro)
        endif
c
        dgx1 = gx1
        dgy1 = gy1
        dgx2 = gx2
        dgy2 = gy2
        call line (dgx1, dgy1, dgx2, dgy2)
        ax = (dy-dx)*pro
        by = -(dy-dx)*pro
        ay = (-dy-dx)*pro
        bx = (-dy-dx)*pro
        call line (dgx2+ax, dgy2+ay, dgx2, dgy2)
        call vector (dgx2+bx, dgy2+by)

c Reset the map to the original map

        call tgsetmap (vl,vr,vb,vt,wl,wr,wb,wt,it)

        if (tgname .eq. 'PLOTV') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  tgmap.f  -  tv80gks mapping routines
c
c  contents	tggetmap  -  Get the linear map used by the plotting library
c		tgsetmap  -  Set the linear map used by the plotting library
c		tgchkmap - check to see if the specified mapping is legal
c		map  -  Define an x and y axes in linear coordinates
c		mapll  -  Define an x and y axes in log-log coordinates
c		mapsl  -  Define an x and y in linear-log coordinates
c		mapls  -  Define an x and y axes in log-linear coordinates
c		maps  -  Define an x and y axes in linear coordinates and 
c			draw ticks
c		mapsll  -  Define an x and y axes in log coordinates and draw
c			ticks
c		mapssl  -  Define an x,y axes in linear-log coordinates and
c			draw ticks
c		mapsls  -  Define an x,y axes in log-linear coordinates and
c			draw ticks
c		mapg  -  Define an x and y axes in linear coordinates and draw
c			grid
c		mapgll  -  Define an x and y axes in log coordinates and
c			draw grid
c		mapgsl  -  Define an x,y axes in linear-log coordinates and
c			draw grid
c		mapgls  -  Define an x,y axes in log-linear coordinates and
c			draw grid
c		mapx  -  Select mapping routine at run time
c		mapp  -  Define an x,y axes as polar
c		maplog  -  Convert log coordinates of a map
c		dders  -  Set clipping on or off
c
c  description  The mapping routines let you define x and y axes and the
c               window in the picture plane into which you wish to draw.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c*************************************************************************

c*************************************************************************
c
c  tggetmap  -  Get the linear map used by the plotting library
c
c  synopsis     call tggetmap (vl,vr,vb,vt,wl,wr,wb,wt,itype)
c               real vl,vr      Viewport x range
c               real vb,vt      Viewport y range
c               real wl,wr      Window x range
c               real wb,wt      Window y range
c               integer itype   Type of map
c
c  description  Allows a single routine to get the map that will be used
c               convert the users coordinates into window the workstation
c               window linear units.
c
c*************************************************************************

        subroutine tggetmap (vl,vr,vb,vt,wl,wr,wb,wt,itype)

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



        vl = xvmin
        vr = xvmax
        vb = yvmin
        vt = yvmax
        wl = xwmin
        wr = xwmax
        wb = ywmin
        wt = ywmax
        itype = maptyp

        return
        end

c*************************************************************************
c
c  tgsetmap  -  Set the linear map used by the plotting library
c
c  synopsis     call tgsetmap (vl,vr,vb,vt,wl,wr,wb,wt,itype)
c               real vl,vr      Viewport x range
c               real vb,vt      Viewport y range
c               real wl,wr      Window x range
c               real wb,wt      Window y range
c               integer itype   Type of map
c
c  description  Allows a single routine to set the map that will be used
c               convert the users coordinates into window the workstation
c               window linear units.
c
c*************************************************************************

        subroutine tgsetmap (vl,vr,vb,vt,wl,wr,wb,wt,itype)
c=======================================================================
c 12/28/2001 ler pppl - eliminate unneeded calls to setusr and setvpt
c=======================================================================

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



        call tgflush

        xvmin = vl
        xvmax = vr
        yvmin = vb
        yvmax = vt
        xwmin = wl
        xwmax = wr
        ywmin = wb
        ywmax = wt
        maptyp = itype

c       call setusr(vl,vr,vb,vt)
c       call setvpt(wl,wr,wb,wt)
c     previous two lines may have to be commented

        call gsvp (1,xvmin,xvmax,yvmin,yvmax)
        call gswn (1,xvmin,xvmax,yvmin,yvmax)
        call gselnt (1)

        return
        end

c*************************************************************************
c
c  tgchkmap - check to see if the specified mapping is legal
c
c  synopsis     call tgchkmap (xleft,xright,ybot,ytop,xmin,xmax,
c		  ymin,ymax,itype)
c
c		real xleft,xright	Users x range
c		real ybot,ytop		Users y range
c		real xmin,xmax		Virtual x range
c		real ymin,ymax		Virtual y range
c		integer itype		The type of mapping
c
c  description  Checks to see if the desired mapping is legal.  Sets the
c		mapping to one if it is not.  The type of maps specified
c		by itype are 1 (linear-linear), 2 (linear-log), 3 (log-
c		linear) and 4 (log-log)
c
c*************************************************************************

        integer function tgchkmap (xl,xr,yb,yt,vl,vr,vb,vt,itype)
        real xl,xr,yb,yt
        real vl,vr,vb,vt
        integer itype
        real epsilon


        epsilon = 1.0e-6


        tgchkmap = 1

c x range should not be equivalent

        if (xl .eq. xr) then
          print *,'left: ',xl,'  right: ',xr
          call tgerror (12)
          return
        endif
        if (vl .eq. vr) then
          print *,'normalized left: ',vl,'  normalized right: ',vr
          call tgerror (12)
          return
        endif

c y range should not be equivalent

        if (yb .eq. yt) then
          print *,'bottom: ',yb,'  top: ',yt
          call tgerror (13)
          return
        endif
        if (vb .eq. vt) then
          print *,'normalized bottom: ',vb,'  normalized top: ',vt
          call tgerror (13)
          return
        endif

c minimum x should be less than maximum	

        if (xl .gt. xr .or. vl .gt. vr) then
          call tgerror (14)
          return
        endif

c minimum y should be less than maximun

        if (yb .gt. yt .or. vb .eq. vt) then
          call tgerror (15)
          return
        endif

c log x should be greater than zero

        if (itype .eq. 3 .or. itype .eq. 4) then
          if (xl .le. 0 .or. xr .le. 0) then
            call tgerror (16)
            return
          endif
        endif

c log y should be greater than zero

        if (itype .eq. 2 .or. itype .eq. 4) then
          if (yb .le. 0 .or. yt .le. 0) then
            call tgerror (17)
            return
          endif
        endif

c virtual x should be in the range of 0 to 1

        if (vl .lt. 0.0 .or. vr .gt. 1.0) then
          call tgerror (18)
          return
        endif

c virtual y should be in the range of 0 to 1

        if (vb .lt. 0.0 .or. vt .gt. 1.0) then
          call tgerror (19)
          return
        endif

c make sure the range has enough significant bits

        denom = 0.5 * (abs(xl) + abs(xr))
        if (abs(xr-xl)/denom .lt. epsilon) then
          print *,'left: ',xl,'  right: ',xr
          call tgerror (22)
          return
        endif

        denom = 0.5 * (abs(yb) + abs(yt))
        if (abs(yt-yb)/denom .lt. epsilon) then
          print *,'bottom: ',yb,'  top: ',yt
          call tgerror (23)
          return
        endif

c mapping has my blessing

        tgchkmap = 0
        return
        end

c*************************************************************************
c
c  map  -  Define an x and y axes in linear coordinates
c
c  synopsis     call map (xleft,xright, ybot,ytop)
c               call map (xleft,xright, ybot,ytop, xmin,xmax, ymin,ymax)
c
c               float xleft,xright,ybot,ytop    World coordinates
c               float xmin,xmax,ymin,ymax       Window coordinates
c
c  description  The world coordinates are mapped into the range specified
c               by the Window.  The window x and y range should be between
c               0.0 and 1.0.  0.0 and 1.0 is the default value.  This
c               routine is very similar to SPPS SET.  The difference is
c               the last argument to set.  A 1 will specify linear x and
c               y mapping.
c
c*************************************************************************

        subroutine map (xleft,xright, ybot,ytop, xmin,xmax, ymin,ymax)

	real*8 xleft, xright, ybot, ytop, xmin, xmax, ymin, ymax
        real xvl,xvr,yvb,yvt
        real xwl,xwr,ywb,ywt
        integer tgchkmap

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'MAP'



	xvl = 0.0
	yvb = 0.0
	xvr = 1.0
	yvt = 1.0
	if (xmin .ne. -1.0) xvl = xmin
	if (ymin .ne. -1.0) yvb = ymin
	if (xmax .ne. -1.0) xvr = xmax
	if (ymax .ne. -1.0) yvt = ymax

        xwl = xleft
        ywb = ybot
        xwr = xright
        ywt = ytop

        if (tgchkmap (xwl,xwr,ywb,ywt,xvl,xvr,yvb,yvt,1) .ne. 0) then
          goto 9999
        endif

        call tgsetmap (xvl,xvr,yvb,yvt,xwl,xwr,ywb,ywt,1)

9999    if (tgname .eq. 'MAP') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  mapll  -  Define an x and y axes in log-log coordinates
c
c  synopsis     call mapll (xleft,xright, ybot,ytop)
c               call mapll (xleft,xright, ybot,ytop, xmin,xmax, ymin,ymax)
c
c               float xleft,xright,ybot,ytop    World coordinates
c               float xmin,xmax,ymin,ymax       Window coordinates
c
c  description  The world coordinates are mapped into the range specified
c               by the Window.  The window x and y range should be between
c               0.0 and 1.0.  0.0 and 1.0 is the default value.
c
c*************************************************************************

        subroutine mapll (xleft,xright, ybot,ytop, xmin,xmax, ymin,ymax)

        real*8 xleft,xright,ybot,ytop
        real*8 xmin,xmax,ymin,ymax
        integer tgchkmap

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'MAPLL'



	xvl = 0.0
	yvb = 0.0
	xvr = 1.0
	yvt = 1.0
	if (xmin .ne. -1.0) xvl = xmin
	if (ymin .ne. -1.0) yvb = ymin
	if (xmax .ne. -1.0) xvr = xmax
	if (ymax .ne. -1.0) yvt = ymax

        xwl = xleft
        ywb = ybot
        xwr = xright
        ywt = ytop

        if (tgchkmap (xwl,xwr,ywb,ywt,xvl,xvr,yvb,yvt,4) .ne. 0) then
          goto 9999
        endif

        call maplog (xwl,xwr)
        call maplog (ywb,ywt)

        call tgsetmap (xvl,xvr,yvb,yvt,xwl,xwr,ywb,ywt,4)

9999    if (tgname .eq. 'MAPLL') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  mapsl  -  Define an x and y in linear-log coordinates
c
c  synopsis     call mapsl (xleft,xright, ybot,ytop)
c               call mapsl (xleft,xright, ybot,ytop, xmin,xmax, ymin,ymax)
c
c               float xleft,xright,ybot,ytop    World coordinates
c               float xmin,xmax,ymin,ymax       Window coordinates
c
c  description  The world coordinates are mapped into the range specified
c               by the Window.  The window x and y range should be between
c               0.0 and 1.0.  0.0 and 1.0 is the default value.
c
c*************************************************************************

        subroutine mapsl (xleft,xright, ybot,ytop, xmin,xmax, ymin,ymax)

        real*8 xleft,xright,ybot,ytop
        real*8 xmin,xmax,ymin,ymax
        integer tgchkmap

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'MAPSL'



	xvl = 0.0
	yvb = 0.0
	xvr = 1.0
	yvt = 1.0
	if (xmin .ne. -1.0) xvl = xmin
	if (ymin .ne. -1.0) yvb = ymin
	if (xmax .ne. -1.0) xvr = xmax
	if (ymax .ne. -1.0) yvt = ymax

        xwl = xleft
        ywb = ybot
        xwr = xright
        ywt = ytop

        if (tgchkmap (xwl,xwr,ywb,ywt,xvl,xvr,yvb,yvt,2) .ne. 0) then
          goto 9999
        endif

        call maplog (ywb,ywt)

        call tgsetmap (xvl,xvr,yvb,yvt,xwl,xwr,ywb,ywt,2)

9999    if (tgname .eq. 'MAPSL') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  mapls  -  Define an x and y axes in log-linear coordinates
c
c  synopsis     call mapls (xleft,xright, ybot,ytop)
c               call mapls (xleft,xright, ybot,ytop, xmin,xmax, ymin,ymax)
c
c               float xleft,xright,ybot,ytop    World coordinates
c
c  description  The world coordinates are mapped into the range specified
c               by the Window.  The window x and y range should be between
c               0.0 and 1.0.  0.0 and 1.0 is the default value.  This
c               routine is very similar to SPPS SET.  The difference is
c               the last argument to set.  A 3 will specify logorithmic x
c               and linear y mapping.
c
c*************************************************************************

        subroutine mapls (xleft,xright, ybot,ytop, xmin,xmax, ymin,ymax)
c=======================================================================
c 12/28/2001 ler pppl - change subroutine parameters from single  
c                        precision to double precision
c=======================================================================

        real*8 xleft,xright,ybot,ytop
        real*8 xmin,xmax,ymin,ymax
        integer tgchkmap

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'MAPLS'



	xvl = 0.0
	yvb = 0.0
	xvr = 1.0
	yvt = 1.0
	if (xmin .ne. -1.0) xvl = xmin
	if (ymin .ne. -1.0) yvb = ymin
	if (xmax .ne. -1.0) xvr = xmax
	if (ymax .ne. -1.0) yvt = ymax

        xwl = xleft
        ywb = ybot
        xwr = xright
        ywt = ytop

        if (tgchkmap (xwl,xwr,ywb,ywt,xvl,xvr,yvb,yvt,3) .ne. 0) then
          goto 9999
        endif

        call maplog (xwl,xwr)

        call tgsetmap (xvl,xvr,yvb,yvt,xwl,xwr,ywb,ywt,3)

9999    if (tgname .eq. 'MAPLS') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  maps  -  Define an x and y axes in linear coordinates and draw ticks
c
c  synopsis     call maps (xleft,xright, ybot,ytop)
c               call maps (xleft,xright, ybot,ytop, xmin,xmax, ymin,ymax)
c
c               float xleft,xright,ybot,ytop    World coordinates
c               float xmin,xmax,ymin,ymax       Window coordinates
c
c  description  The world coordinates are mapped into the range specified
c               by the Window.  The window x and y range should be between
c               0.0 and 1.0.  0.0 and 1.0 is the default value.  This
c               routine is very similar to SPPS SET.  The difference is
c               the last argument to set.  A 1 will specify linear x and
c               y mapping.
c
c*************************************************************************

        subroutine maps (xleft,xright,ybot,ytop,xmin,xmax,ymin,ymax)

        real*8 xleft,xright,ybot,ytop
        real*8 xmin,xmax,ymin,ymax
        integer tgchkmap

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'MAPS'



	xvl = 0.11328
	yvb = 0.11328
	xvr = 1.0
	yvt = 1.0
	if (xmin .ne. -1.0) xvl = xmin
	if (ymin .ne. -1.0) yvb = ymin
	if (xmax .ne. -1.0) xvr = xmax
	if (ymax .ne. -1.0) yvt = ymax

        xwl = xleft
        ywb = ybot
        xwr = xright
        ywt = ytop

        if (tgchkmap (xwl,xwr,ywb,ywt,xvl,xvr,yvb,yvt,1) .ne. 0) then
          goto 9999
        endif

        call tgsetmap (xvl,xvr,yvb,yvt,xwl,xwr,ywb,ywt,1)

        call maplab (xwl,xwr,0,1,1)
        call maplab (ywb,ywt,0,0,1)

9999    if (tgname .eq. 'MAPS') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  mapsll  -  Define an x and y axes in log coordinates and draw ticks
c
c  synopsis     call mapsll (xleft,xright, ybot,ytop)
c               call mapsll (xleft,xright,ybot,ytop,xmin,xmax,ymin,ymax)
c
c               float xleft,xright,ybot,ytop    World coordinates
c               float xmin,xmax,ymin,ymax       Window coordinates
c
c  description  The world coordinates are mapped into the range specified
c               by the Window.  The window x and y range should be between
c               0.0 and 1.0.  0.0 and 1.0 is the default value.  This
c               routine is very similar to SPPS SET.  The difference is
c               the last argument to set.  A 1 will specify linear x and
c               y mapping.
c
c*************************************************************************

        subroutine mapsll (xleft,xright,ybot,ytop,xmin,xmax,ymin,ymax)

        real*8 xleft,xright,ybot,ytop
        real*8 xmin,xmax,ymin,ymax
        integer tgchkmap

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'MAPSLL'



	xvl = 0.11328
	yvb = 0.11328
	xvr = 1.0
	yvt = 1.0
	if (xmin .ne. -1.0) xvl = xmin
	if (ymin .ne. -1.0) yvb = ymin
	if (xmax .ne. -1.0) xvr = xmax
	if (ymax .ne. -1.0) yvt = ymax

        xwl = xleft
        ywb = ybot
        xwr = xright
        ywt = ytop

        if (tgchkmap (xwl,xwr,ywb,ywt,xvl,xvr,yvb,yvt,4) .ne. 0) then
          goto 9999
        endif

        call maplog (xwl,xwr)
        call maplog (ywb,ywt)

        call tgsetmap (xvl,xvr,yvb,yvt,xwl,xwr,ywb,ywt,4)

        call maplab (xwl,xwr,1,1,1)
        call maplab (ywb,ywt,1,0,1)

9999    if (tgname .eq. 'MAPSLL') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  mapssl  -  Define an x,y axes in linear-log coordinates and draw ticks
c
c  synopsis     call mapssl (xleft,xright, ybot,ytop)
c               call mapssl (xleft,xright,ybot,ytop,xmin,xmax,ymin,ymax)
c
c               float xleft,xright,ybot,ytop    World coordinates
c               float xmin,xmax,ymin,ymax       Window coordinates
c
c  description  The world coordinates are mapped into the range specified
c               by the Window.  The window x and y range should be between
c               0.0 and 1.0.  0.0 and 1.0 is the default value.  This
c               routine is very similar to SPPS SET.  The difference is
c               the last argument to set.  A 1 will specify linear x and
c               y mapping.
c
c*************************************************************************

        subroutine mapssl (xleft,xright,ybot,ytop,xmin,xmax,ymin,ymax)

        real*8 xleft,xright,ybot,ytop
        real*8 xmin,xmax,ymin,ymax
        integer tgchkmap

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'MAPSSL'



	xvl = 0.11328
	yvb = 0.11328
	xvr = 1.0
	yvt = 1.0
	if (xmin .ne. -1.0) xvl = xmin
	if (ymin .ne. -1.0) yvb = ymin
	if (xmax .ne. -1.0) xvr = xmax
	if (ymax .ne. -1.0) yvt = ymax

        xwl = xleft
        ywb = ybot
        xwr = xright
        ywt = ytop

        if (tgchkmap (xwl,xwr,ywb,ywt,xvl,xvr,yvb,yvt,2) .ne. 0) then
          goto 9999
        endif

        call maplog (ywb,ywt)

        call tgsetmap (xvl,xvr,yvb,yvt,xwl,xwr,ywb,ywt,2)

        call maplab (xwl,xwr,0,1,1)
        call maplab (ywb,ywt,1,0,1)

9999    if (tgname .eq. 'mapssl') tgname = 'NONE'
        return
        end

c*************************************************************************
c
c  mapsls  -  Define an x,y axes in log-linear coordinates and draw ticks
c
c  synopsis     call mapsls (xleft,xright, ybot,ytop)
c               call mapsls (xleft,xright,ybot,ytop,xmin,xmax,ymin,ymax)
c
c               float xleft,xright,ybot,ytop    World coordinates
c               float xmin,xmax,ymin,ymax       Window coordinates
c
c  description  The world coordinates are mapped into the range specified
c               by the Window.  The window x and y range should be between
c               0.0 and 1.0.  0.0 and 1.0 is the default value.  This
c               routine is very similar to SPPS SET.  The difference is
c               the last argument to set.  A 1 will specify linear x and
c               y mapping.
c
c*************************************************************************

        subroutine mapsls (xleft,xright,ybot,ytop,xmin,xmax,ymin,ymax)

        real*8 xleft,xright,ybot,ytop
        real*8 xmin,xmax,ymin,ymax
        integer tgchkmap

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'MAPSLS'



	xvl = 0.11328
	yvb = 0.11328
	xvr = 1.0
	yvt = 1.0
	if (xmin .ne. -1.0) xvl = xmin
	if (ymin .ne. -1.0) yvb = ymin
	if (xmax .ne. -1.0) xvr = xmax
	if (ymax .ne. -1.0) yvt = ymax

        xwl = xleft
        ywb = ybot
        xwr = xright
        ywt = ytop

        if (tgchkmap (xwl,xwr,ywb,ywt,xvl,xvr,yvb,yvt,3) .ne. 0) then
          goto 9999
        endif

        call maplog (xwl,xwr)

        call tgsetmap (xvl,xvr,yvb,yvt,xwl,xwr,ywb,ywt,3)

        call maplab (xwl,xwr,1,1,1)
        call maplab (ywb,ywt,0,0,1)

9999    if (tgname .eq. 'MAPSLS') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  mapg  -  Define an x and y axes in linear coordinates and draw grid
c
c  synopsis     call mapg (xleft,xright, ybot,ytop)
c               call mapg (xleft,xright, ybot,ytop, xmin,xmax, ymin,ymax)
c
c               float xleft,xright,ybot,ytop    World coordinates
c               float xmin,xmax,ymin,ymax       Window coordinates
c
c  description  The world coordinates are mapped into the range specified
c               by the Window.  The window x and y range should be between
c               0.0 and 1.0.  0.0 and 1.0 is the default value.  This
c               routine is very similar to SPPS SET.  The difference is
c               the last argument to set.  A 1 will specify linear x and
c               y mapping.
c
c*************************************************************************

        subroutine mapg (xleft,xright,ybot,ytop,xmin,xmax,ymin,ymax)

        real*8 xleft,xright,ybot,ytop
        real*8 xmin,xmax,ymin,ymax
        integer tgchkmap

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'MAPG'



	xvl = 0.11328
	yvb = 0.11328
	xvr = 1.0
	yvt = 1.0
	if (xmin .ne. -1.0) xvl = xmin
	if (ymin .ne. -1.0) yvb = ymin
	if (xmax .ne. -1.0) xvr = xmax
	if (ymax .ne. -1.0) yvt = ymax

        xwl = xleft
        ywb = ybot
        xwr = xright
        ywt = ytop

        if (tgchkmap (xwl,xwr,ywb,ywt,xvl,xvr,yvb,yvt,1) .ne. 0) then
          goto 9999
        endif

        call tgsetmap (xvl,xvr,yvb,yvt,xwl,xwr,ywb,ywt,1)

        call maplab (xwl,xwr,0,1,0)
        call maplab (ywb,ywt,0,0,0)

9999    if (tgname .eq. 'MAPG') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  mapgll  -  Define an x and y axes in log coordinates and draw grid
c
c  synopsis     call mapgll (xleft,xright, ybot,ytop)
c               call mapgll (xleft,xright,ybot,ytop,xmin,xmax,ymin,ymax)
c
c               float xleft,xright,ybot,ytop    World coordinates
c               float xmin,xmax,ymin,ymax       Window coordinates
c
c  description  The world coordinates are mapped into the range specified
c               by the Window.  The window x and y range should be between
c               0.0 and 1.0.  0.0 and 1.0 is the default value.  This
c               routine is very similar to SPPS SET.  The difference is
c               the last argument to set.  A 1 will specify linear x and
c               y mapping.
c
c*************************************************************************

        subroutine mapgll (xleft,xright,ybot,ytop,xmin,xmax,ymin,ymax)

        real*8 xleft,xright,ybot,ytop
        real*8 xmin,xmax,ymin,ymax
        integer tgchkmap

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'MAPGLL'


	xvl = 0.11328
	yvb = 0.11328
	xvr = 1.0
	yvt = 1.0
	if (xmin .ne. -1.0) xvl = xmin
	if (ymin .ne. -1.0) yvb = ymin
	if (xmax .ne. -1.0) xvr = xmax
	if (ymax .ne. -1.0) yvt = ymax

        xwl = xleft
        ywb = ybot
        xwr = xright
        ywt = ytop

        if (tgchkmap (xwl,xwr,ywb,ywt,xvl,xvr,yvb,yvt,4) .ne. 0) then
          goto 9999
        endif

        call maplog (xwl,xwr)
        call maplog (ywb,ywt)

        call tgsetmap (xvl,xvr,yvb,yvt,xwl,xwr,ywb,ywt,4)

        call maplab (xwl,xwr,1,1,0)
        call maplab (ywb,ywt,1,0,0)

9999    if (tgname .eq. 'MAPGLL') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  mapgsl  -  Define an x,y axes in linear-log coordinates and draw grid
c
c  synopsis     call mapgsl (xleft,xright, ybot,ytop)
c               call mapgsl (xleft,xright,ybot,ytop,xmin,xmax,ymin,ymax)
c
c               float xleft,xright,ybot,ytop    World coordinates
c               float xmin,xmax,ymin,ymax       Window coordinates
c
c  description  The world coordinates are mapped into the range specified
c               by the Window.  The window x and y range should be between
c               0.0 and 1.0.  0.0 and 1.0 is the default value.  This
c               routine is very similar to SPPS SET.  The difference is
c               the last argument to set.  A 1 will specify linear x and
c               y mapping.
c
c*************************************************************************

        subroutine mapgsl (xleft,xright,ybot,ytop,xmin,xmax,ymin,ymax)

        real*8 xleft,xright,ybot,ytop
        real*8 xmin,xmax,ymin,ymax
        integer tgchkmap

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'MAPGSL'



	xvl = 0.11328
	yvb = 0.11328
	xvr = 1.0
	yvt = 1.0
	if (xmin .ne. -1.0) xvl = xmin
	if (ymin .ne. -1.0) yvb = ymin
	if (xmax .ne. -1.0) xvr = xmax
	if (ymax .ne. -1.0) yvt = ymax

        xwl = xleft
        ywb = ybot
        xwr = xright
        ywt = ytop

        if (tgchkmap (xwl,xwr,ywb,ywt,xvl,xvr,yvb,yvt,2) .ne. 0) then
          goto 9999
        endif

        call maplog (ywb,ywt)

        call tgsetmap (xvl,xvr,yvb,yvt,xwl,xwr,ywb,ywt,2)

        call maplab (xwl,xwr,0,1,0)
        call maplab (ywb,ywt,1,0,0)

9999    if (tgname .eq. 'MAPGSL') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  mapgls  -  Define an x,y axes in log-linear coordinates and draw grid
c
c  synopsis     call mapgls (xleft,xright, ybot,ytop)
c               call mapgls (xleft,xright,ybot,ytop,xmin,xmax,ymin,ymax)
c
c               float xleft,xright,ybot,ytop    World coordinates
c               float xmin,xmax,ymin,ymax       Window coordinates
c
c  description  The world coordinates are mapped into the range specified
c               by the Window.  The window x and y range should be between
c               0.0 and 1.0.  0.0 and 1.0 is the default value.  This
c               routine is very similar to SPPS SET.  The difference is
c               the last argument to set.  A 1 will specify linear x and
c               y mapping.
c
c*************************************************************************

        subroutine mapgls (xleft,xright,ybot,ytop,xmin,xmax,ymin,ymax)

        real*8 xleft,xright,ybot,ytop
        real*8 xmin,xmax,ymin,ymax
        integer tgchkmap

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'MAPGLS'



	xvl = 0.11328
	yvb = 0.11328
	xvr = 1.0
	yvt = 1.0
	if (xmin .ne. -1.0) xvl = xmin
	if (ymin .ne. -1.0) yvb = ymin
	if (xmax .ne. -1.0) xvr = xmax
	if (ymax .ne. -1.0) yvt = ymax

        xwl = xleft
        ywb = ybot
        xwr = xright
        ywt = ytop

        if (tgchkmap (xwl,xwr,ywb,ywt,xvl,xvr,yvb,yvt,3) .ne. 0) then
          goto 9999
        endif

        call maplog (xwl,xwr)

        call tgsetmap (xvl,xvr,yvb,yvt,xwl,xwr,ywb,ywt,3)

        call maplab (xwl,xwr,1,1,0)
        call maplab (ywb,ywt,0,0,0)

9999    if (tgname .eq. 'MAPGLS') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  mapx  -  Select mapping routine at run time
c
c  synopsis     call mapx (imap,xleft,xright, ybot,ytop)
c               call mapx (imap,xleft,xright,ybot,ytop,xmin,xmax,ymin,ymax)
c
c               integer imap                    Mappimg routine
c               float xleft,xright,ybot,ytop    World coordinates
c               float xmin,xmax,ymin,ymax       Window coordinates
c
c  description  The world coordinates are mapped into the range specified
c               by the Window.  The window x and y range should be between
c               0.0 and 1.0.  0.0 and 1.0 is the default value.  This
c               routine is very similar to SPPS SET.  The difference is
c               the last argument to set.  A 1 will specify linear x and
c               y mapping.
c
c               imap Routine  imap Routine  imap Routine
c               ------------  ------------  ------------
c               1    map      5    mapg     9    maps
c               2    mapll    6    mapgll   10   mapsll
c               3    mapsl    7    mapgsl   11   mapssl
c               4    mapls    8    mapgls   12   mapsls
c
c*************************************************************************

        subroutine mapx (imap,xleft,xright,ybot,ytop,
     +                   xmin,xmax,ymin,ymax)

        integer imap
        real*8 xleft,xright,ybot,ytop
        real*8 xmin,xmax,ymin,ymax
        real*8 wl,wr,wb,wt,vl,vr,vb,vt

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'MAPX'



        if (imap .lt. 1 .or. imap .gt. 12) then
          call tgerror (11)
          goto 9999
        endif

        if (imap .gt. 4) then
          vl = 0.11328
          vb = 0.11328
          vr = 1.0
          vt = 1.0
        else
          vl = 0.0
          vb = 0.0
          vr = 1.0
          vt = 1.0
        endif
	if (xmin .ne. -1.0) vl = xmin
	if (ymin .ne. -1.0) vb = ymin
	if (xmax .ne. -1.0) vr = xmax
	if (ymax .ne. -1.0) vt = ymax

        wl = xleft
        wb = ybot
        wr = xright
        wt = ytop

        goto (1,2,3,4,5,6,7,8,9,10,11,12),imap

1       call map (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999

2       call mapll (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999

3       call mapsl (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999

4       call mapls (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999

5       call mapg (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999

6       call mapgll (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999

7       call mapgsl (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999

8       call mapgls (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999

9       call maps (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999

10      call mapsll (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999

11      call mapssl (wl,wr,wb,wt,vl,vr,vb,vt)
        goto 9999

12      call mapsls (wl,wr,wb,wt,vl,vr,vb,vt)

9999    if (tgname .eq. 'MAPX') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  mapp  -  Define an x,y axes as polar
c
c  synopsis     call mapp (rmax)
c               call mapp (rmax,xmin)
c               call mapp (rmax,xmin,xmax)
c               call mapp (rmax,xmin,xmax,ymin)
c
c               float rmax                      Max radius
c               float xmin,xmax,ymin            Window coordinates
c
c  description  rmax is the maximum radial distance.  If positive, a grid
c               is drawn.  xmin and xmax define the range in NDC where the
c               plot will be placed.  ymin is the minimum y in NDC.  ymax
c               is ymin+(xmax-xmin) since the plot must be square.  The
c               default size for xmin,xmax,ymin,ymax is .001 to .999
c
c*************************************************************************

        subroutine mapp (rmax,xmin,xmax,ymin)

        real*8 rmax
        real*8 xmin,xmax,ymin
        integer tgchkmap

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'MAPP'



        xvl = 0.001
        yvb = 0.001
        xvr = 0.999

        xwr = abs(rmax)
        ywt = xwr
        xwl = -xwr
        ywb = -xwr

        if (ymin .ne. -1.0) yvb = ymin
        if (xmax .ne. -1.0) xvr = xmax
	if (xmin .ne. -1.0) xvl = xmin

        yvt = yvb+(xvr-xvl)

c Although this is a polar mapping, check it as if it were linear

        if (tgchkmap (xwl,xwr,ywb,ywt,xvl,xvr,yvb,yvt,1) .ne. 0) then
          goto 9999
        endif

        call tgsetmap (xvl,xvr,yvb,yvt,xwl,xwr,ywb,ywt,5)

        if (rmax .gt. 0) call maplab (0.,xwr,0,2,1)

9999    if (tgname .eq. 'MAPP') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  maplog  -  Convert log coordinates of a map
c
c  synopsis     call maplog (xlow,xhi)
c
c               float xlow                      Low end of range
c               float xhi                       Hi end of range
c
c  description  Converts the low and high end of a range into a log
c               range that will have integer exponents.  This algorithm
c               was taken directly from tv80lib.
c
c*************************************************************************

        subroutine maplog (xlow, xhi)
c
cc purpose: make sure that xlow is less than xhi
c           and generate the log of the values if necessary
c
c if xlow eq xhi, make xhi a smidge bigger
c
        if (xlow .gt. xhi) then
          t = xlow
          xlow = xhi
          xhi = t
        else if (xlow .eq. xhi) then
          xhi = xlow + 1.0
        endif
c
c reset xlow and xhi to non-negative values
c
        if (xlow .lt. 0.0) xlow = 0.0
        if (xhi .lt. 0.0) xhi = 0.0
c
c reset xlow and xhi to default values if xhi has an bad value
c
        if (xhi .eq. 0) then
          xlow = .00001
          xhi = 100000.
          call tgerror (8)
          go to 110
        endif

c reset xlow with respect to xhi

        if (xlow .eq. 0.0) then
          xlow = xhi*0.0000000001
          call tgerror (9)
        endif

c get log value for xlow

110     xlow = log10a(xlow)

c get log value for xhi

        g1 = alog10(xhi)
        if (g1 .gt. 0.0) g1 = g1+0.9999999999
        xhi = float(int(g1))

c reset xlow

        if (xlow .eq. xhi) xlow = xhi-1.0
c
        return
        end

c*************************************************************************
c
c  dders  -  Set clipping on or off
c
c  synopsis     call dders (iflag)
c
c               integer iflag                   Clip indicator
c
c  description  If iflag is -1, the clipping to the view port is turned on.
c               If iflag is 0 (default), the clipping is turned off.  The
c               original version of tv80lib had +1 which indicated no
c               clipping at all.
c
c*************************************************************************

        subroutine dders (iflag)

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'DDERS'



c Flush buffers that are assuming current clipping mode

        call tgflush

c If -1 or 0 then turn clipping on, if 1 then turn clipping off.

        if (iflag .eq. -1 .or. iflag .eq. 0) then
          iclip = 1
          call gsclip (iclip)
        else if (iflag .eq. 1) then
          iclip = 0
          call gsclip (iclip)
        else
          call tgerror (10)
        endif

9999    if (tgname .eq. 'DDERS') tgname = 'NONE'
        return

        end

c**************************************************************************
c
c  tggks.f  -  Routines to initialize the gks drivers
c
c  contents	The routines in this file open GKS if it is not already
c		opened and then open a specific workstation given a
c		workstation identifier and a connection id.  The names
c		and the workstation that is opened will change depending
c		on what workstations are available with your implementation
c		of GKS.  After initializing a workstation, tginit should
c		be called to put the workstation in a state that is
c		consistent with what the library thinks it should be.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************


c**************************************************************************
c
c  ncarcgm -  NCAR CGM binary
c
c  description	NCARGKS will generate a GKS binary metafile.  It does not
c		set the filename in the standar gks fashion.  Version 3.1
c		uses an environment variable called NCARG_GKS_OUTPUT.  By
c		setting this variable, I can name the output file.  The
c		standard GKS method of calling GSCNID does not work.  It
c		does not exist in the GKS written by NCAR.  The connection
c		identifier in NCAR is the unit number that will be used
c		when the file is created.
c
c**************************************************************************

        subroutine ncarcgm (iwk,connam)
        character*(*) connam

	character*80 filnam

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'NCARCGM'



c Open gks if it was not done already

        call gqops (iopstate)
        if (iopstate .eq. 0) call gopks (6,0)

c Open the workstation.  The second parameter is the connection identifier.
c In NCAR GKS, this is actually the number that will be used as the unit
c number for the output file when the OPEN statement is performed.  Use
c the same number for the connection identifier as for the workstation
c identifier.  Note that on many systems, unit 5 and 6 are reserved for
c output and input.

        filnam = connam
        call gesc (-1391,1,filnam,0,idum,filnam)
        call gopwk (iwk,iwk,1)
        call gacwk (iwk)

        call tginit (iwk)

        if (tgname .eq. 'NCARCGM') tgname = 'NONE'
        return

        end

c*************************************************************************
c
c  tgstat.F  -  Status changing routines
c
c  contents     gstat -  Change status of device
c
c  description  These routines change the status of some parameters inside
c               of NCAR.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c*************************************************************************

c*************************************************************************
c
c  gstat  -  Change status of device
c
c  synopsis     call gstat (iwkid,istat)
c
c               integer iwkid           Workstation id
c               integer istat           Status
c
c  description  If istat is negative, then close the workstation.  If it
c               is 0, then activate the workstation.  If it is positive,
c               then activate it.
c
c*************************************************************************

        subroutine gstat (iwkid, istat)

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'GSTAT'



        call tgflush

c gdacwk requires a workstation be in an active state.  Check first before
c deactivation.

        if (istat .lt. 0) then
          call gqwks (iwkid,ierr,iwkstate)
          if (iwkstate .eq. 1) call gdawk (iwkid)
          call gclwk (iwkid)
        elseif (istat .eq. 0) then
          call gqwks (iwkid,ierr,iwkstate)
          if (iwkstate .eq. 1) call gdawk (iwkid)
        else
          call gacwk (iwkid)
        endif

9999    if (tgname .eq. 'GSTAT') tgname = 'NONE'
        return

        end

c**********************************************************************
c
c  tgtran.F  -  transformation routines
c
c  contents	cntxy  -  convert normalized point to a transformed point
c		init2d  -  Initialize transformation matrix
c		tran2d  -  Set translation
c		scal2d  -  Set scale
c		rot2d  -  Set rotation
c		center  -  Set the center of scalings and rotations
c		catmat  -  Concatenate matrices
c		makidmat  -  Make an identity matrix
c		multmat  -  Multiply the transformation matrix by another
c			matrix
c
c  description  All routines in tglib which do transformations will
c               need to be fed through the transformation routines
c               to obtain the correct coordinates.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**********************************************************************

c**********************************************************************
c
c  cntxy  -  convert normalized point to a transformed point
c
c  synopsis     call cntxy (x,y,xnew,ynew)
c               real xnew,ynew       transformed x,y
c               real x,y             Point to transform
c
c**********************************************************************

        subroutine cntxy (x,y,xnew,ynew)

c The temporary variables have two purposes.  Obviously they save the
c converted value of x and y so that the conversion routines don't
c have to be called twice.  By saving the values in temporary variables
c the input x,y can be the same as the output x,y.  For example:
c   call cntxy (x,y,x,y)

        real xtemp,ytemp

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



        if (.not. dotrans) then
          xnew = x
          ynew = y
          return
        endif

        if (.not. matmade) call catmat

        xtemp = x
        ytemp = y
        xnew = xtemp*trnmat(1,1) + ytemp*trnmat(2,1) + trnmat(3,1)
        ynew = xtemp*trnmat(1,2) + ytemp*trnmat(2,2) + trnmat(3,2)

        return
        end

c**********************************************************************
c
c  init2d  -  Initialize transformation matrix
c
c  synopsis     call init2d
c
c**********************************************************************

        subroutine init2d ()

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'INIT2D'



        dotrans = .FALSE.
        transx = 0.0
        transy = 0.0
        scalex = 1.0
        scaley = 1.0
        rotat = 0.0
        centrx = 0.0
        centry = 0.0

c The angle of text changes with the rotation.

        call tgchset (ichfont,ichindx,ichangle)

9999    if (tgname .eq. 'INIT2D') tgname = 'NONE'
        return

        end

c**********************************************************************
c
c  tran2d  -  Set translation
c
c  synopsis     call tran2d (tx,ty)
c               real tx,ty      Translation in x,y
c
c  description  Sets variable to let the library know that it must
c               do the transformation and sets the variables.  Since
c               the amount of translation is simply an offset to add
c               to all coordinates, it really didn't make sense for
c               it to be a log value when the users coordinate system
c               is logorithmic.  It is therefore linear.  The amount
c               to translate in a log system can be found by taking
c               the log of the position that you would like to move
c               to and the position to move from and then subtract
c               the two values.
c
c**********************************************************************

        subroutine tran2d (tx,ty)

	real*8 tx, ty


c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'TRAN2D'



        dotrans = .TRUE.
        matmade = .FALSE.

        transx = tx
        transy = ty

9999    if (tgname .eq. 'TRAN2D') tgname = 'NONE'
        return

        end

c**********************************************************************
c
c  scal2d  -  Set scale
c
c  synopsis     call scal2d (sx,sy)
c               real sx,sy      Scale factor in x,y
c
c  description  Sets variable to let the library know that it must
c               do the scaling and sets the variables
c
c**********************************************************************

        subroutine scal2d (sx,sy)

	real*8 sx, sy


c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'SCAL2D'



        dotrans = .TRUE.
        matmade = .FALSE.

        scalex = sx
        scaley = sy

9999    if (tgname .eq. 'SCAL2D') tgname = 'NONE'
        return

        end

c**********************************************************************
c
c  rot2d  -  Set rotation
c
c  synopsis     call rot2d (angle)
c               real angle      angle of rotation about the center
c
c  description  Sets variable to let the library know that it must
c               do the rotation and sets the angle in radians.  The
c               user specifies the angle in degrees.
c
c**********************************************************************

        subroutine rot2d (angle)

	real*8 angle


c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'ROT2D'



        dotrans = .TRUE.
        matmade = .FALSE.

        rotat = angle * 3.14159265358979323846/180.0

c The angle of text changes with the rotation.

        call tgchset (ichfont,ichindx,ichangle)

9999    if (tgname .eq. 'ROT2D') tgname = 'NONE'
        return

        end

c**********************************************************************
c
c  center  -  Set the center of scalings and rotations
c
c  synopsis     call center (cx,cy)
c               real cx,cy      x,y of the center
c
c  description  Sets variable to let the library know that it must
c               do the transformation and sets the center.  cx and cy
c               are in the users coordinate system.  They are converted
c               to normalized coordinated when the transformation
c               matrix is created in catmat.
c
c**********************************************************************

        subroutine center (cx,cy)

	real*8 cx, cy


c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'CENTER'



        dotrans = .TRUE.
        matmade = .FALSE.

        centrx = cx
        centry = cy

9999    if (tgname .eq. 'CENTER') tgname = 'NONE'
        return

        end

c**********************************************************************
c
c  catmat  -  Concatenate matrices
c
c  synopsis     call catmat ()
c
c  description  Concatenates the matrices to create one transformation
c               matrix.  The effect will be to move the center to the
c               origin, scale, rotate, move back out to the original
c               position, and then move the amount specified in the
c               translation.
c
c**********************************************************************

        subroutine catmat
c=======================================================================
c 12/28/2001 ler pppl - Change first two argments in 'call culxy' from
c                        single precision to double precision
c=======================================================================

        real trxnew,trynew
        real ctxnew,ctynew

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



        real matrix (3,3)
	real*8 dcentrx, dcentry


c convert translation and center into the proper units

        call clnxy (xwmin+transx,ywmin+transy,trxnew,trynew)
        trxnew = trxnew - xvmin
        trynew = trynew - yvmin

        dcentrx = centrx
        dcentry = centry
        call culxy (dcentrx,dcentry,ctxnew,ctynew)
        call clnxy (ctxnew,ctynew,ctxnew,ctynew)

c erase old trasformation matrix

        call makidmat (trnmat)

c move to the origin

        call makidmat (matrix)
        matrix (3,1) = -ctxnew
        matrix (3,2) = -ctynew
        call multmat (matrix)

c scale

        call makidmat (matrix)
        matrix (1,1) = scalex
        matrix (2,2) = scaley
        call multmat (matrix)

c rotate

        call makidmat (matrix)
        matrix (1,1) = cos (rotat)
        matrix (1,2) = sin (rotat)
        matrix (2,1) = - matrix(1,2)
        matrix (2,2) = matrix(1,1)
        call multmat (matrix)

c move center back where it was, add in the translation

        call makidmat (matrix)
        matrix (3,1) = ctxnew + trxnew
        matrix (3,2) = ctynew + trynew
        call multmat (matrix)

        matmade = .TRUE.

        return
        end

c**********************************************************************
c
c  makidmat  -  Make an identity matrix
c
c  synopsis     call makidmat (matrix)
c               real matrix(3,3)        Matrix to change into identity
c
c**********************************************************************

        subroutine makidmat (matrix)

        real matrix(3,3)

        do 10 i = 1,3
          do 20 j = 1,3
            if (i .eq. j) then
              matrix(i,j) = 1.0
            else
              matrix(i,j) = 0.0
            endif
20        continue
10      continue

        return
        end

c**********************************************************************
c
c  multmat  -  Multiply the transformation matrix by another matrix
c
c  synopsis     call multmat (matrix)
c               real matrix(3,3)        Matrix to multiply by
c
c  description  Combines the passed matrix into the transformation
c               matrix by multiplying the two.
c
c**********************************************************************

        subroutine multmat (matrix)

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



        real matrix(3,3), temp(3,3)

        do 10 i = 1,3
          do 20 j = 1,3
            temp(i,j) = trnmat(i,1)*matrix(1,j) + trnmat(i,2)
     +        *matrix(2,j) + trnmat(i,3)*matrix(3,j)
20        continue
10      continue

        do 30 i = 1,3
          do 40 j = 1,3
            trnmat(i,j) = temp(i,j)
40        continue
30      continue

        return
        end

c**************************************************************************
c
c  tgtv80.f - Routines used in the tv80 lib that shouldn't be there at all
c
c  contents	cartmm - calculate min and max of an array
c
c  description  The routines in this file are here for backwards support only.
c               These routines should not be used by any TV80 program.
c		However, some programmers used them, so they are here.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        subroutine cartmm (n, xmin, xmax, x, inc)
c
cc variable declarations:
c
        real x(*)

        tmin = x(1)
        tmax = x(1)
        do 1 i = 1, n, inc
          if (x(i) .lt. tmin) tmin = x(i)
          if (x(i) .gt. tmax) tmax = x(i)
    1   continue
        xmin = tmin
        xmax = tmax
9999    return
        END


c*************************************************************************
c
c  tgutil.F  -  Tv80 to gks shell utility routines
c
c  contents     tgtoup          Convert string to uppercase
c               tgtolow         Convert string to lowercase
c               clnxy           Convert linear to normalized
c               cufxy           Convert user to transformed normalized
c               culxy           Convert user to linear
c		tgerror		Print an error message
c		tgisfont	Test for validity of a font
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c*************************************************************************

c*************************************************************************
c
c  tgtoup  -  convert a string to upper case
c
c  synopsis     call tgtoup (string)
c               character*(*) string    String to convert
c
c*************************************************************************

        subroutine tgtoup (string,outstr,num)
        character*(*) string
        character*(*) outstr

        do 10 i=1,num
          if (ichar(string(i:i)) .ge. ichar('a') .and.
     +        ichar(string(i:i)) .le. ichar('z')) then
            outstr(i:i) = char (ichar(string(i:i))-32)
          else
            outstr(i:i) = string(i:i)
          endif
10      continue

        return
        end

c*************************************************************************
c
c  tgtolow  -  convert a string to lower case
c
c  synopsis     call tgtolow (string)
c               character*(*) string    String to convert
c
c*************************************************************************

        subroutine tgtolow (string,outstr,num)
        character*(*) string
        character*(*) outstr

        do 10 i=1,num
          if (ichar(string(i:i)) .ge. ichar('A') .and.
     +        ichar(string(i:i)) .le. ichar('Z')) then
            outstr(i:i) = char (ichar(string(i:i))+32)
          else
            outstr(i:i) = string(i:i)
          endif
10      continue

        return
        end

c*************************************************************************
c
c  clnxy  -  Convert linear coordinates to normalized
c
c  synopsis     call clnxy (x,y,xnew,ynew)
c               real x,y        Linear x,y
c               real xnew,ynew  Normalized point
c
c  description  someday
c
c*************************************************************************

        subroutine clnxy (x,y,xnew,ynew)
        real x,y,xnew,ynew

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



        xnew = (x - xwmin) / (xwmax - xwmin) * (xvmax - xvmin) + xvmin
        ynew = (y - ywmin) / (ywmax - ywmin) * (yvmax - yvmin) + yvmin

        return
        end

c*************************************************************************
c
c  cufxy  -  Convert user coordinates to something resonable
c
c  synopsis     call cufxy (x,y,xnew,ynew)
c               real x,y        Users x,y
c               real xnew,ynew  Transformed point
c
c  description  someday
c
c*************************************************************************

        subroutine cufxy (x,y,xnew,ynew)
c        real x,y,xnew,ynew
	real*8 x,y
	real xnew,ynew



c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



        call culxy (x,y,xnew,ynew)
        call clnxy (xnew,ynew,xnew,ynew)

        if (dotrans) call cntxy (xnew,ynew,xnew,ynew)

        return
        end

c*************************************************************************
c
c  culxy  -  Convert user coordinates to linear-linear
c
c  synopsis     call culxy (x,y,xnew,ynew)
c               real x,y        Users x,y
c               real xnew,ynew  Transformed point
c
c  description  someday
c
c*************************************************************************

        subroutine culxy (x,y,xnew,ynew)
c       real x,y,xnew,ynew
	real*8 x,y
        real xtemp,ytemp

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



        if (maptyp .eq. 1) then

          xtemp = x
          ytemp = y

        else if (maptyp .eq. 2) then

          xtemp = x
          if (y .le. 0.0) then
            print *,'Error in CULXY: Invalid log coordinate'
            ytemp = ywmin
          else
c           ytemp = alog10(y)
           ytemp = dlog10(y)
          endif

        else if (maptyp .eq. 3) then

          if (x .le. 0.0) then
            print *,'Error in CULXY: Invalid log coordinate'
            xtemp = xwmin
          else
c           xtemp = alog10(x)
           xtemp = dlog10(x)
          endif
          ytemp = y

        else if (maptyp .eq. 4) then

          if (x .le. 0.0 .or. y .le. 0.0) then
            print *,'Error in CULXY: Invalid log coordinate'
            xtemp = xwmin
            ytemp = ywmin
          else
c           xtemp = alog10(x)
           xtemp = dlog10(x)
c            ytemp = alog10(y)
            ytemp = dlog10(y)
          endif

        else if (maptyp .eq. 5) then

          xtemp = x * cos(y)
          ytemp = x * sin(y)

        endif

        xnew = xtemp
        ynew = ytemp

        return
        end

c*************************************************************************
c
c  tgerror  -  Print an error message
c
c  synopsis     call tgerror (ierror)
c               integer ierror	Error number
c
c  description  Prints the name stored in tgname, followed by the
c		error message.  Negative error numbers indicate error
c		messages that are somehow dependent on LRLTRAN.  Positive
c		numbers indicate more realistic errors.
c
c*************************************************************************

        subroutine tgerror (ierror)
        integer ierror

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



        if (ierror .eq. -1) then
          print *,tgname,': Invalid number of arguments'

        else if (ierror .eq. -2) then
          print *,tgname,': Invalid camera id'

        else if (ierror .eq. -3) then
          print *,tgname,': Unused workstation identifier not found'

        else if (ierror .eq. 1) then
          print *,tgname,': Invalid color'

        else if (ierror .eq. 2) then
          print *,tgname,': Bad format'

        else if (ierror .eq. 3) then
          print *,tgname,': no more than 30 labels'

        else if (ierror .eq. 4) then
          print *,tgname,': no less than 2 labels'

        else if (ierror .eq. 5) then
          print *,tgname,': bdiv(2) .le. bdiv(1)'

        else if (ierror .eq. 6) then
          print *,tgname,': axis line is a point'

        else if (ierror .eq. 7) then
          print *,tgname,': string exceeded maximum of 256'

        else if (ierror .eq. 8) then
          print *,tgname,
     +      ': Negative coordinates reset to 1e-5 to 1e+5'

        else if (ierror .eq. 9) then
          print *,tgname,': Minimum log map value reset'

        else if (ierror .eq. 10) then
          print *,tgname,': Invalid value for argument'

        else if (ierror .eq. 11) then
          print *,tgname,': Invalid map number'

        else if (ierror .eq. 12) then
          print *,tgname,': Bad mapping. left .eq. right'

        else if (ierror .eq. 13) then
          print *,tgname,': Bad mapping. bottom .eq. top'

        else if (ierror .eq. 14) then
          print *,tgname,': Bad mapping. left .gt. right'

        else if (ierror .eq. 15) then
          print *,tgname,': Bad mapping. bottom .gt. top'

        else if (ierror .eq. 16) then
          print *,tgname,
     +      ': Bad mapping. Log of negative x is undefined'

        else if (ierror .eq. 17) then
          print *,tgname,
     +      ': Bad mapping. Log of negative y is undefined'

        else if (ierror .eq. 18) then
          print *,tgname,
     +      ': Bad mapping. Virtual x should be 0 <= x <= 1'

        else if (ierror .eq. 19) then
          print *,tgname,
     +      ': Bad mapping. Virtual y should be 0 <= y <= 1'

        else if (ierror .eq. 20) then
          print *,tgname,': GKS not opened'

        else if (ierror .eq. 21) then
          print *,tgname,': Library not initialized'

        else if (ierror .eq. 22) then
          print *,tgname,': Insufficient significant bits in x range'

        else if (ierror .eq. 23) then
          print *,tgname,': Insufficient significant bits in y range'

        endif

        return
        end

c*************************************************************************
c
c  tgisfont  -  test if a font is valid
c
c  synopsis	logical tgisfont
c		if (tgisfont (ifont)) newfont = ifont
c               integer ifont		Font number to test
c
c  description  Checks to see if a font is valid.  Returns .FALSE. if it
c		isn't and .TRUE. if it is.  This routine is dependent
c		upon a specific implementation of GKS.
c
c*************************************************************************

        logical function tgisfont (ifont)
        integer ifont

c include the site dependent info

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



        tgisfont = .TRUE.
        if (ifont .lt. minfont .or. ifont .gt. maxfont) then
	  tgisfont = .FALSE.
	endif

        return
        end


c**************************************************************************
c
c  tggks.f  -  Routines to initialize the gks drivers
c
c  contents	The routines in this file open GKS if it is not already
c		opened and then open a specific workstation given a
c		workstation identifier and a connection id.  The names
c		and the workstation that is opened will change depending
c		on what workstations are available with your implementation
c		of GKS.  After initializing a workstation, tginit should
c		be called to put the workstation in a state that is
c		consistent with what the library thinks it should be.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

c**************************************************************************
c
c  cgmbin -  CGM binary
c
c  description	Open GKS if necessary, get the identity of the workstation
c		set the output file name, open and then activate the
c		workstation.
c
c**************************************************************************

        subroutine cgmbin (iwk,connam)
        character*(*) connam

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'CGMBIN'



c Open gks if it was not done already

        call gqops (iopstate)
        if (iopstate .eq. 0) call gopks (6,0)

        call ncarcgm (iwk,connam)
        call gopwk (iwk,iwk,10100)

        call gacwk (iwk)

        call tginit (iwk)

        if (tgname .eq. 'CGMBIN') tgname = 'NONE'
        return

        end

c**************************************************************************
c
c  tgopen -  open some GKS workstation
c
c  description	Open GKS if necessary, open and then activate the
c		workstation.
c
c**************************************************************************

        subroutine tgopen (iwk,icnid,itype,connam)
        character*(*) connam

c include the standard common block

c**************************************************************************
c
c  tgcommon  -  TV80 to GKS common blocks
c
c  description  This is used by routines in tv80gks to make a common
c               place for accessing commons.
c
c    National Energy Research Supercomputer Center
c    Lawrence Livermore National Laboratory
c    University of California
c    Livermore, California 94550, USA
c
c    (C) Copyright 1991 The Regents of the University of California.
c    All Rights Reserved.
c
c**************************************************************************

        parameter (MAXPNT=200)
        integer WRDSIZ
        parameter (WRDSIZ=8)

        logical dotrans
        real transx,transy
        real scalex,scaley
        real rotat
        real centrx,centry
        real trnmat(3,3)
        logical matmade
        character*8 tgname
        logical doinit

        common /tgcblk/ dotrans,transx,transy,scalex,scaley,rotat,
     +    centrx,centry,trnmat,matmade,tgname,doinit

        real xvmin,xvmax,yvmin,yvmax
        real xwmin,xwmax,ywmin,ywmax
        integer maptyp,iclip

        common /tgmap/ xvmin,xvmax,yvmin,yvmax,xwmin,xwmax,ywmin,ywmax,
     +    maptyp,iclip

        real chx,chy
        integer ichcase
        integer ichindx
        real chrot
        integer ichangle
        real chparm(4,4)
        real chupx,chupy
        integer ichfont
        real chaddx,chaddy
        logical autofeed

        common /tgchar/ chx,chy,ichcase,ichindx,chrot,chparm,
     +                  ichfont,chupx,chupy,ichangle,chaddx,chaddy,
     +                  ichclip,minfont,maxfont,autofeed

        real red(16),green(16),blue(16)
        integer icurclr

        common /tgcolr/ icurclr,red,green,blue

        integer numpnt
        real xpnt(MAXPNT),ypnt(MAXPNT)
        integer ilncase,ilnfont,ilnindx,ilnspace,iptspace
        common /tgpnts/ numpnt,xpnt,ypnt,ilncase,ilnfont,
     +                  ilnindx,ilnspace,iptspace



c Set library entry routine name

        if (tgname .eq. 'NONE') tgname = 'TGOPEN'



c Open gks if it was not done already

        call gqops (iopstate)
        if (iopstate .eq. 0) call gopks (6,0)

        call ncarcgm (icnid,connam)
        call gopwk (iwk,icnid,itype)

        call gacwk (iwk)

        call tginit (iwk)

        if (tgname .eq. 'TGOPEN') tgname = 'NONE'
        return

        end

c
c....added autofeed to tgchar common in gaxisi
c....added block data statement to define ndimen
c...changed inta to int on lines 3168 and 3172
