c     
c..........................................................
      subroutine testvec 
c..........................................................
c     
      include 'vacuum1.inc'
      include 'vacuum2.inc'
      include 'vacuum5.inc'
c     
      dimension xpla(nths), zpla(nths),
     $     xob(nths),zob(nths), the(nths)
      dimension axrl(nths,nfm), axil(nths,nfm),
     $     azrl(nths,nfm), azil(nths,nfm), bdotnr(nfm), bdotni(nfm)
c     
      equivalence (xpla,xinf), (zpla,zinf)
c     
      jmax1 = lmax(1) - lmin(1) + 1
      mth12 = 2*mth
c
c.....Put these in input and defolt later.
c
      mdel = 16
      jdel = 8
c     
      do 10 i = 1, mth1
         the(i) = (i-1) * dth
 10   continue
c     
c     if ( farwal ) mth12 = mth
c     
c.....set up window.
c$$$c     
c$$$      write ( outmod, '("I, xp, zp, xw, zw, xpp, zpp, xwp, zwp =")' )
c$$$      write ( iotty,  '("I, xp, zp, xw, zw, xpp, zpp, xwp, zwp =")' )
c$$$      do i = 1, mth1, 8
c$$$      write ( outmod, '(i3, 1x, 8f8.2)' )
c$$$     $        i, xpla(i),  zpla(i),  xwal(i),  zwal(i),
c$$$     $           xplap(i), zplap(i), xwalp(i), zwalp(i)
c$$$      write ( iotty,  '(i3, 1x, 8f8.2)' )
c$$$     $        i, xpla(i),  zpla(i),  xwal(i),  zwal(i),
c$$$     $           xplap(i), zplap(i), xwalp(i), zwalp(i)
c$$$      end do
c$$$c
      call bounds(xpla,zpla,1,mth,xmnp,xmxp,zmnp,zmxp)
      xmin = xmnp
      xmax = xmxp
      zmin = zmnp
      zmax = zmxp
c     
      plrad = 0.5 * ( xmxp - xmnp )
      xmaj  = 0.5 * ( xmxp + xmnp )
c     
c     if ( farwal ) go to 20
c     
      call bounds(xwal,zwal,1,mw,xmnw,xmxw,zmnw,zmxw)
c     
      xmin = amin1(xmnp,xmnw)
      xmax = amax1(xmxp,xmxw)
      zmin = amin1(zmnp,zmnw)
      zmax = amax1(zmxp,zmxw)
c     
 20   continue
c     
      dx = xmax - xmin
c     dz = zmax - zmin
      dz = 2.0 * zmax
      dxz = amax1(dx,dz)
      xmaxw = 0.5 * ( xmax + xmin ) + 0.50*dxz
      xminw = 0.5 * ( xmax + xmin ) - 0.50*dxz
      zmaxw = 0.50 * dxz
      zminw = - zmaxw
c     
      xming = xminw - 0.05*(xmaxw-xminw)
      xmaxg = xmaxw + 0.05*(xmaxw-xminw)
      zming = zminw - 0.05*(zmaxw-zminw)
      zmaxg = zmaxw + 0.05*(zmaxw-zminw)
c     
      call maps(xming,xmaxg,zming,zmaxg,0.1,0.815,0.285,1.0)
c     
c.... plot conductors and plasma shape..
c     
      call points ( xpla,zpla, mth, -1,-1,0.,0. )
      call points ( xwal,zwal, mth, -1,-1,0.,0. )
c     
      call frame(1)
c
c      phi = pye / 4.0
      phi = 30.0 * pye / 180.0
c     
      nobs = 0
      mthdel = mth/mdel
      jmxdel = jmax1/jdel
      mdeld = max0(1,mthdel)
      jdeld = max0(1,jmxdel)
c 
      write ( outmod,'(a)' ) char(12)
c
      write ( outmod, 8 ) n, phi, mdeld, jdeld
      write ( iotty, 8 ) n, phi, mdeld, jdeld
 8    format ( "n = ",f8.2, " phi = ",f8.3, " mdeld = ", i3,
     $     " jdeld = ", i3, / )
c
      do 100 iobs = 1, mth, mdeld
         nobs = nobs + 1
         xob(nobs) = xwal(iobs)
         zob(nobs) = zwal(iobs)
 100  continue
c
      if ( .not. checks ) go to 300
c     
      isg = -1
      call vecpot ( xwal,zwal,mth, xpla,zpla,xplap,zplap,
     $     mth, isg, phi,axrl,axil, azrl,azil )
c
      do 23 jj = 1, jmax1
         axrl(mth1,jj) = axrl(1,jj)
         axil(mth1,jj) = axil(1,jj)
         azrl(mth1,jj) = azrl(1,jj)
         azil(mth1,jj) = azil(1,jj)
 23   continue
c
      do 200 iobs = 1, mth, mdeld
         write ( outmod, 11 ) iobs, (axrl(iobs,ll), ll = 1,jmax1,jdeld)
 11      format ( i4,("   axrl ",1p10e11.4,/) )
         write ( outmod, 12 ) iobs, (axil(iobs,ll), ll = 1,jmax1,jdeld)
 12      format ( i4,("   axil ",1p10e11.4,/) )
         write ( outmod, 15 ) iobs, (azrl(iobs,ll), ll = 1,jmax1,jdeld)
 15      format ( i4,("   azrl ",1p10e11.4,/) )
         write ( outmod, 16 ) iobs, (azil(iobs,ll), ll = 1,jmax1,jdeld)
 16      format ( i4,("   azil ",1p10e11.4,/) )
c     
c......Calculate normal derivative...
c     
         do 80 ll = 1, jmax1
            bdotnr(ll) = - n * ( axil(iobs,ll)*xwalp(iobs)
     $           + azil(iobs,ll)*zwalp(iobs) )
            bdotni(ll) =   n * ( axrl(iobs,ll)*xwalp(iobs)
     $           + azrl(iobs,ll)*zwalp(iobs) )
 80      continue
         write ( outmod, 13 ) iobs, (bdotnr(ll), ll = 1,jmax1,jdeld)
 13      format ( i4,(" bdotnr ",1p10e11.4,/) )
         write ( outmod, 33 ) iobs, (bdotni(ll), ll = 1,jmax1,jdeld)
 33      format ( i4,(" bdotni ",1p10e11.4,/) )
c     
 200  continue
c
 300  continue
c     
      if ( (lspark .eq. 1) .or. (lspark .eq. 2) ) call spark2d

c  Turn off SPARK3D for now to save space
c      Will take ot out of the MAKEFILES    March 09, 2000

c$$$      if ( (lspark .eq. 1) .or. (lspark .eq. 3) ) call spark3d
c
      return
      end
c     
c.........................................................
      subroutine vecpot ( xobs,zobs,nobs, xsce,zsce,xscp,zscp,nsce,
     $     isg, phi, axrl,axil, azrl,azil )
c..........................................................
c     
      include 'vacuum1.inc'
      include 'vacuum3.inc'
      include 'vacuum5.inc'
c     
      dimension xobs(*),zobs(*),xsce(*),zsce(*),xscp(*),zscp(*),
     $     axrl(nths,nfm),axil(nths,nfm),
     $     azrl(nths,nfm),azil(nths,nfm),
     $     dchxr(nths,nfm),dchxi(nths,nfm),
     $     dchzr(nths,nfm),dchzi(nths,nfm)
c     
      jmax1 = lmax(1) - lmin(1) + 1
c     
      sinnph = sin(n*phi)
      cosnph = cos(n*phi)
c     
      call grchi(xobs,zobs,nobs, xsce,zsce,xscp,zscp,nsce,
     $     isg,cplar,cplai,cslth,snlth, 1,1,
     $     cpwr,cpwi, dchxr,dchxi,dchzr,dchzi)
c     
      if ( abs(n) .lt. 1.e-10 ) then
         write ( outmod, 11 ) 
         write ( iotty,  11 ) 
 11      format ( /,1x, "!!! N= 0 in subroutine VECPOT !!! " )
         return
      end if
c     
      do 200 iobs = 1, nobs
         do 100 ll = 1, jmax1
            axrl(iobs,ll) = - xobs(iobs) * ( dchzr(iobs,ll)*sinnph
     $           - dchzi(iobs,ll)*cosnph ) / n
            axil(iobs,ll) = - xobs(iobs) * ( dchzr(iobs,ll)*cosnph
     $           + dchzi(iobs,ll)*sinnph ) / n
c     
            azrl(iobs,ll) =   xobs(iobs) * ( dchxr(iobs,ll)*sinnph
     $           - dchxi(iobs,ll)*cosnph ) / n
            azil(iobs,ll) =   xobs(iobs) * ( dchxr(iobs,ll)*cosnph
     $           + dchxi(iobs,ll)*sinnph ) / n
 100     continue
 200  continue
c     
      return
      end
c     
c........................................................
      subroutine spark2d
c........................................................
c     
      include 'vacuum1.inc'
      include 'vacuum2.inc'
      include 'vacuum3.inc'
      include 'vacuum5.inc'
      include 'vacuum6.inc'
c
      dimension grdgre(nths2,nths2), grwp(nths,nths), grri(nths2,nfm2)
      common / bigv1 / grdgre, grwp, grri
c     
      dimension xpla(nths),zpla(nths),bxr(nths,nfm),bxi(nths,nfm),
     $     bzr(nths,nfm),bzi(nths,nfm),
     $     bwr(nths,nfm), bwi(nths,nfm),
     $     dcc(nfm,nfm),dss(nfm,nfm),sdsc(nfm,nfm),sdsci(nfm,nfm),
     $     unit(nfm,nfm),uvpr(nfm,nfm),uvpi(nfm,nfm),
     $     uvwwr(nfm,nfm),uvwwi(nfm,nfm),
     $     uvwr(nfm,nfm),uvwi(nfm,nfm),
     $     uvpwr(nfm,nfm),uvpwi(nfm,nfm),
     $     uvp0(nfm,nfm),uvw0(nfm,nfm),uvpw0(nfm,nfm),
     $     vcmtp0(nfm,nfm),vcmtw0(nfm,nfm),
     $     vacmts(nfm,nfm),dvacmat(nfm,nfm),
     $     wrki(nfm),waa(nfm,nfm),wbb(nfm,nfm),copmat(nfm,nfm),
     $     waa1(nfmsq),wbb1(nfmsq),wrk0(nfm,nfm)
      dimension vacpstr(nfmsq),vacpsti(nfmsq)
      dimension wkxx(nths), wkyy(nths)
c     
      equivalence (xpla,xinf), (zpla,zinf)
      equivalence (waa,waa1), (wbb,wbb1)
      equivalence (vacpstr,wrk0), (vacpsti,copmat)
c     
      ntso = 9
      jmax1 = lmax(1) - lmin(1) + 1
      lrnge = jmax1
      ln = lmin(1)
      lx = lmax(1)
c
c************************   STEP I   ***********************************
c
c               The response matrices CPLAR, CPLAI, from the plasma
c               surface with no wall. Already done in the main part of 
c               the vacuum calculation.
c
      write ( outmod,'(a)' ) char(12)
c     
      write ( outmod, 131 ) 
 131  format (/,1x,"STEP I. The response matrices CPLAR, CPLAI,",/,
     $     "       from the plasma surface with no wall.",/,
     $     "       Already done in the main part of the vacuum",/,
     $     "       calculation." )
c     
      if ( checks ) then
      call matwrtn ( cplar,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"cplar",
     $     outmod, 0 )
      call matwrtn ( cplai,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"cplai",
     $     outmod,0 )
c
c.. Fourier in combination:
      call fanal2 ( cplar,cplai,nths,nfm,0, waa,wbb,nfm,nfm )
         call matwrtn(waa,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Real Cpla(l,l)", outmod,0)
         call matwrtn(wbb,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Imag Cpla(l,l)", outmod,0 )
      end if
c
c.....Save p-p contrib. in cpp. Total will be in cpwt w-p in STEP V
c
      do i = 1, nths
         do m =1, nfm
            cppr(i,m) = cplar(i,m)
            cppi(i,m) = cplai(i,m)
         end do
      end do
c
c************************   STEP II  **********************************     
c
c.......STEP II. Chi and grad-chi at the virtual wall position using the 
c                plasma sources.  Application of Green's Second 
c                Identity again.
c
c      if ( checks ) write ( outmod,'(a)' ) char(12)
c     
      write ( outmod, 132 ) 
 132  format (/,1x,"STEP II. Chi and Grchi at virtual wall" )
c     
      isg = -1
      call grchi ( xwal,zwal,mth, xpla,zpla,xplap,zplap,mth,
     $     isg, cplar,cplai,
     $     cslth,snlth, 1,1, cpwr,cpwi, bxr,bxi,bzr,bzi )
c
      if ( lff .eq. 1 ) then
   
c... Multiply by FEEDBACK FACTOR, ff.
      call scalmul ( bxr, nths,nfm, ff )
      call scalmul ( bxi, nths,nfm, ff )
      call scalmul ( bzr, nths,nfm, ff )
      call scalmul ( bzi, nths,nfm, ff )
      call scalmul ( cpwr, nths,nfm, ff )
      call scalmul ( cpwi, nths,nfm, ff )

      end if
c
      if ( lff .eq. 2 ) then
c
c... Multiply by FEEDBACK FACTOR, ff. Except for l=0 rows.
      call scalmul0 ( bxr, nths,nfm, ln,lx, ff )
      call scalmul0 ( bxi, nths,nfm, ln,lx, ff )
      call scalmul0 ( bzr, nths,nfm, ln,lx, ff )
      call scalmul0 ( bzi, nths,nfm, ln,lx, ff )
      call scalmul0 ( cpwr, nths,nfm, ln,lx, ff )
      call scalmul0 ( cpwi, nths,nfm, ln,lx, ff )
c
      end if
c
      if ( lff .eq. 3 ) then
c
c... Multiply by FEEDBACK VECTOR FV.
      call vecmula ( bxr, nths,nfm, ln,lx, fv )
      call vecmula ( bxi, nths,nfm, ln,lx, fv )
      call vecmula ( bzr, nths,nfm, ln,lx, fv )
      call vecmula ( bzi, nths,nfm, ln,lx, fv )
      call vecmula ( cpwr, nths,nfm, ln,lx, fv )
      call vecmula ( cpwi, nths,nfm, ln,lx, fv )
c
      end if
c
c....... Store the p-w contrib. of the scalar potential in cwwt 
c        for diagnostics.
c..      The w-w contrib. will be added in STEP IV. The signs in cwwt 
c        make this correspond to the Direct method.
c
c......  Store also the field from the plasma at the wall position. Used for 
c        feedback studies.
c
      do i = 1, nths
         do m = 1, nfm
            cwwtr(i,m) = - cpwr(i,m)
            cwwti(i,m) = - cpwi(i,m)
            bxpwr(i,m) = bxr(i,m)
            bxpwi(i,m) = bxi(i,m)
            bzpwr(i,m) = bzr(i,m)
            bzpwi(i,m) = bzi(i,m)
         end do
      end do
c     
      if ( checks ) then
      call matwrtn ( cpwr,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"cpwr",
     $     outmod,0 )
      call matwrtn ( cpwi,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"cpwi",
     $     outmod,0 )
c
c.. Fourier in combination:
      call fanal2 ( cpwr,cpwi,nths,nfm,0, waa,wbb,nfm,nfm )
         call matwrtn(waa,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Real Cpw(l,l)",outmod,0 )
         call matwrtn(wbb,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Imag Cpw(l,l)",outmod,0 )
c
      call matwrtn ( bxr,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"bxr",
     $     outmod,0 )
      call matwrtn ( bxi,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"bxi",
     $     outmod,0 )
      call matwrtn ( bzr,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"bzr",
     $     outmod,0 )
      call matwrtn ( bzi,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"bzi",
     $     outmod,0 )
      end if
c     
c*************************   STEP III       **************************
c
c                 The magnetic field normal components along the surface of
c                 the virtual wall position.
c
c      if ( checks ) write ( outmod,'(a)' ) char(12)
c     
      write ( outmod, 133 ) 
 133  format (/,1x,"STEP III. Normal B field at wall." )
c     
      call bwri ( xwal,zwal,mth, xwalp,zwalp, bxr,bxi,bzr,bzi,
     $     jmax1, bwr,bwi )
c     
c.....  Store normal field at virtual wall due the plasma. 
c       Used for feedback studies.
c
      do i = 1, nths
         do m = 1, nfm
            bnpwr(i,m) = bwr(i,m)
            bnpwi(i,m) = bwi(i,m)
         end do
      end do
c     
      if ( checks ) then
      call matwrtn ( bwr,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"bwr",
     $     outmod,0 )
      call matwrtn ( bwi,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"bwi",
     $     outmod,0 )
c
c.. Fourier in combination:
      call fanal2 ( bwr,bwi,nths,nfm,0, waa,wbb,nfm,nfm )
         call matwrtn(waa,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Real Bw(l,l)", outmod,0 )
         call matwrtn(wbb,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Imag Bw(l,l)", outmod,0 )
      end if
c
c     
c*************************   STEP IV      ***************************
c
c                The chi wall-wall response at the wall using BWR, BWI
c                from step III as boundary conditions, and assuming 
c                no plasma.
c
c      if ( checks ) write ( outmod,'(a)' ) char(12)
c     
      write ( outmod, 134 ) 
 134  format (/,1x,"STEP IV.  Chi wall-wall response at wall" )
c     
      call chiwal ( xwal,zwal, xwal,zwal, bwr,bwi,2,2,1,1,
     $     cwalr,cwali )
c     
c.......Add w-w contrib. to p-w contrib. from STEP II to the total chiw:-
c
      do i = 1, nths
         do m = 1, nfm
            cwwtr(i,m) = cwwtr(i,m) + cwalr(i,m)
            cwwti(i,m) = cwwti(i,m) + cwali(i,m)
         end do
      end do
c           
      if ( checks ) then
      call matwrtn ( cwalr,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"cwalr",
     $     outmod,0 )
      call matwrtn ( cwali,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"cwali",
     $     outmod,0 )
c Fourier analyse wall response at setp IV.
c Cos and SIn parts of both cwalr and cwali:
      call fanal0 ( cwalr,nths,nfm,0, waa,wbb,nfm,nfm )
         call matwrtn(waa,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Cos-Cwalr(l,l)", outmod,0 )
         call matwrtn(wbb,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Sin-Cwalr(l,l)", outmod,0 )
      call fanal0 ( cwali,nths,nfm,0, waa,wbb,nfm,nfm )
         call matwrtn(waa,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Cos-Cwali(l,l)", outmod,0 )
         call matwrtn(wbb,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Sin-Cwali(l,l)", outmod,0 )
c
c.. Fourier in combination:
      call fanal2 ( cwalr,cwali,nths,nfm,0, waa,wbb,nfm,nfm )
         call matwrtn(waa,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Real Cwal(l,l)", outmod,0 )
         call matwrtn(wbb,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Imag Cwal(l,l)", outmod,0 )
c
      end if
c
c..... The wall self-energy.
c
      call uvacww ( bwr,bwi, cwalr,cwali, jmax1, uvwwr,uvwwi )
c     
c$$$      call matwrtn ( uvwwr,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
c$$$     $     "uvwwr at VI", outmod,0 )
c$$$c     
c$$$      call vacasym ( uvwwr, nfm,jmax1,"Uvwwr", outmod,iotty )
c$$$      call masym0 ( uvwwr, nfm,jmax1,ln,ln,0,0,
c$$$     $     "Vac_WWr", outmod,iotty )
c$$$c
c$$$      call matwrtn ( uvwwi,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
c$$$     $     "uvwwi at VI", outmod,0 )
c     
c************************   STEP V    ******************************
c
c               Chi and Grad-chi at the virtual plasma position using 
c               the wall sources calculated from  step IV. Green's 
c               second identity again.
c
c      if ( checks ) write ( outmod,'(a)' ) char(12)
c     
      write ( outmod, 135 ) 
 135  format (/,1x,"STEP V. Chi and Grchi at plasma position." )
c     
      isg = -1
      ipin = 1
      ipm = -1
      write ( outmod,101 ) isg, ipin,ipm
      write ( iotty, 101 ) isg, ipin, ipm
 101  format ( /,1x, "isg, ipin, ipm = ", 3i3 )
      call grchi ( xpla,zpla,mth, xwal,zwal,xwalp,zwalp,mth,
     $     isg, cwalr,cwali,
     $     bwr,bwi, ipin,ipm, cpwr,cpwi, bxr,bxi,bzr,bzi )
c     
c...  Add in w-p contrib. to the p-p contrib. from STEP I to the total chip:-
c     Save the chi w-p contribution for feedback, 10-20-97
c
      do i = 1, nths
         do m = 1, nfm
            cpwtr(i,m) = cppr(i,m) + cpwr(i,m)
            cpwti(i,m) = cppi(i,m) + cpwi(i,m)
c
c$$$            cwptr(i,m) = cpwtr(i,m) ! Save the total chi in the plasma.
c$$$            cwpti(i,m) = cpwti(i,m) ! Save the total chi in the plasma.
            cwptr(i,m) = cpwr(i,m)  ! Save the p-w contribution.
            cwpti(i,m) = cpwi(i,m)  ! Save the p-w contribution.
         end do
      end do

      if ( checks ) then
      call matwrtn ( cpwr,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"cpwr at V",
     $     outmod,0 )
      call matwrtn ( cpwi,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"cpwi at V",
     $     outmod,0 )
      call matwrtn ( bxr,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"bxr at V",
     $     outmod,0 )
      call matwrtn ( bxi,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"bxi at V",
     $     outmod,0 )
      call matwrtn ( bzr,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"bzr at V",
     $     outmod,0 )
      call matwrtn ( bzi,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"bzi at V",
     $     outmod,0 )
      end if
c     
c************************   STEP VI   ************************
c
c              The unnormalized vacuum delta-W.(wrt. b.c.'s)

c      if ( checks ) write ( outmod,'(a)' ) char(12)
c     
      write ( outmod, 136 ) 
 136  format (/,1x,"STEP VI. Vacuum Delta-w w/o b.c.'s " )
c     
      call uvacw ( jmax1, uvpr,uvpi,uvwr,uvwi,uvpwr,uvpwi )
c     
      if ( checks ) then
      call matwrtn ( uvpr,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "uvpr at VI", outmod,0 )
      call vacasym ( uvpr, nfm,jmax1,"Uvpr", outmod,iotty )
      call matwrtn ( uvpi,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "uvpi at VI", outmod,0 )
c
      call matwrtn ( uvwr,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $     "uvwr at VI", outmod,0 )
      call vacasym ( uvwr, nfm,jmax1,"Uvwr", outmod,iotty )
      call matwrtn ( uvwi,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $     "uvwi at VI", outmod,0 )
c
c      write ( outmod,'(a)' ) char(12)
      call matwrtn ( uvpwr,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $     "Unnormalized uVPWr at VI", outmod,0 )
      call vacasym ( uvpwr, nfm,jmax1,"uVPWr", outmod,iotty )
      call matwrtn ( uvpwi,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $     "uVPWi at VI", outmod,0 )
      end if
c     
c************************   STEP VII   ***************************
c
c               The magnetic field components at the virtual surface of
c               the plasma.
c     
      write ( outmod, 137 ) 
 137  format (/,1x,"STEP VII. B components at plasma surface." )
c     
      call bwri ( xpla,zpla,mth, xplap,zplap, bxr,bxi,bzr,bzi,
     $     jmax1, bwr,bwi )
c     
c.....  Store normal field at virtual plasma due the wall. 
c       Used for feedback studies.
c
      do i = 1, nths
         do m = 1, nfm
            bnwpr(i,m) = bwr(i,m)
            bnwpi(i,m) = bwi(i,m)
         end do
      end do
c     
      if ( checks ) then
      call matwrtn ( bwr,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"bwr",
     $     outmod,0 )
      call matwrtn ( bwi,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"bwi",
     $     outmod,0 )
      end if
c     
c************************   STEP VIII    **************************
c
c                  Fourier analysis for the normalizing D-matrices.
c     
      write ( outmod, 138 ) 
 138  format (/,1x,"STEP VIII. Fourier analysis for D-matrices." )
c     
      call dmats ( bwr,bwi, dcc,dss, sdsc )
c     
      if ( checks ) then
      call matwrtn ( dcc,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,"dcc",
     $     outmod,0 )
      call matwrtn ( dss,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,"dss",
     $     outmod,0 )
      end if
c      write ( outmod,'(a)' ) char(12)
      call matwrtn ( sdsc,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $     "Normalizing D matrix, SDSC",
     $     outmod,0 )
c     
c***********************   STEP IX   ****************************
c
c                  Vacuum Delta-W. First, without (l-nq)'s.
c
c      if ( checks ) write ( outmod,'(a)' ) char(12)
c     
      write ( outmod, 139 ) 
 139  format (/,1x,"STEP IX. Vacuum Delta-W. Without (l-nq)." )
c
c......Normalize the vacuum matrices for the correct boundary conditions.
c
      ifail = 0
      call macopy ( sdsc,nfm, copmat,nfm )
c      call gelima ( copmat,nfm,uvpr,nfm,jmax1,jmax1,uvp0,nfm,wrki,
c     $     waa,nfm,wbb,nfm,ifail )
      call f04aef ( copmat,nfm,uvpr,nfm,jmax1,jmax1,uvp0,nfm,wrki,
     $     waa,nfm,wbb,nfm,ifail )
c     
      write ( outmod, 109 ) ifail
      write ( 6, 109 ) ifail
 109  format ( /,1x, "ifail for uvp0 = ", i4 )
c     
      ifail = 0
      call macopy ( sdsc,nfm, copmat,nfm )
c      call gelima ( copmat,nfm,uvwr,nfm,jmax1,jmax1,uvw0,nfm,wrki,
c     $     waa,nfm,wbb,nfm,ifail )
      call f04aef ( copmat,nfm,uvwr,nfm,jmax1,jmax1,uvw0,nfm,wrki,
     $     waa,nfm,wbb,nfm,ifail )
c     
      write ( outmod, 111 ) ifail
      write ( 6, 111 ) ifail
 111  format ( 1x, "ifail for uvw0 = ", i4 )
c     
c.....use transpose of uvpwr, or sdsc
c
      call mtrans ( uvpwr, nfm, jmax1 )
      write ( outmod, '(/,"UVPWR transposed",/)' )
      write ( iotty,  '(/,"UVPWR transposed",/)' )
      call mtrans ( sdsc,  nfm, jmax1 )
      write ( outmod, '(/,"SDSC transposed",/)' )
      write ( iotty,  '(/,"SDSC transposed",/)' )
      ifail = 0
      call macopy ( sdsc,nfm, copmat,nfm )
c      call gelimb ( copmat,nfm,uvpwr,nfm,jmax1,jmax1,uvpw0,nfm,wrki,
c     $     ifail )
      call f04aaf ( copmat,nfm,uvpwr,nfm,jmax1,jmax1,uvpw0,nfm,wrki,
     $     ifail )
c     $     waa,nfm,wbb,nfm,ifail )
c     
      write ( outmod, 112 ) ifail
      write ( 6, 112 ) ifail
 112  format ( 1x, "ifail for uvpw0 = ", i4 )
c     
      if ( checks ) then
      call matwrtn ( uvp0,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $     "Vacuum P-Delta_w w/o l-nq", outmod,0 )
      call vacasym ( uvp0, nfm,jmax1,"Uvp0", outmod,iotty )
c     
      call matwrtn ( uvw0,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $     "Vacuum W-Delta_w w/o l-nq", outmod,0 )
      call vacasym ( uvw0, nfm,jmax1,"Uvw0", outmod,iotty )
c
c      write ( outmod,'(a)' ) char(12)
      call matwrtn ( uvpw0,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $     "Vacuum PW-Delta_w w/o l-nq", outmod,0 )
      call vacasym ( uvpw0, nfm,jmax1,"Vacuum PW-Delta_w w/o l-nq",
     $     outmod,iotty )
      call masym0 ( uvpw0, nfm,jmax1,ln,ln,0,0, 1.0,
     $     "Vacuum PW-Delta_w w/o l-nq", outmod,iotty )
      call vacasymi ( uvpw0, nfm,jmax1,"Vacuum PW-Delta_w w/o l-nq",
     $     outmod,iotty )
      end if
c     
c.....Normalize the total chi on the wall:-
c........The real part and imaginary part
c        SDSC already transposed.
c
      call chinrm ( cwwtr, nths,nfm, mth,jmax1, sdsc,copmat,
     $     wrktr1,wrktr2,wrki,outmod )
      call chinrm ( cwwti, nths,nfm, mth,jmax1, sdsc,copmat,
     $     wrktr1,wrktr2,wrki,outmod )
c
c.....Now, the chi on the plasma form the plasma:-
c.... cplar will be overwritten in Main. FIx later.
c
      call chinrm ( cplar, nths,nfm, mth,jmax1, sdsc,copmat,
     $     wrktr1,wrktr2,wrki,outmod )
      call chinrm ( cplai, nths,nfm, mth,jmax1, sdsc,copmat,
     $     wrktr1,wrktr2,wrki,outmod )
      call chinrm ( cppr, nths,nfm, mth,jmax1, sdsc,copmat,
     $     wrktr1,wrktr2,wrki,outmod )
      call chinrm ( cppi, nths,nfm, mth,jmax1, sdsc,copmat,
     $     wrktr1,wrktr2,wrki,outmod )
c
c.....Now, the total chi on the plasma:-
c
      call chinrm ( cpwtr, nths,nfm, mth,jmax1, sdsc,copmat,
     $     wrktr1,wrktr2,wrki,outmod )
      call chinrm ( cpwti, nths,nfm, mth,jmax1, sdsc,copmat,
     $     wrktr1,wrktr2,wrki,outmod )
c
c.....Now, the fields on the wall from the plasma:-
c
      call chinrm ( bxpwr, nths,nfm, mth,jmax1, sdsc,copmat,
     $     wrktr1,wrktr2,wrki,outmod )
      call chinrm ( bxpwi, nths,nfm, mth,jmax1, sdsc,copmat,
     $     wrktr1,wrktr2,wrki,outmod )
      call chinrm ( bzpwr, nths,nfm, mth,jmax1, sdsc,copmat,
     $     wrktr1,wrktr2,wrki,outmod )
      call chinrm ( bzpwi, nths,nfm, mth,jmax1, sdsc,copmat,
     $     wrktr1,wrktr2,wrki,outmod )
      call chinrm ( bnpwr, nths,nfm, mth,jmax1, sdsc,copmat,
     $     wrktr1,wrktr2,wrki,outmod )
      call chinrm ( bnpwi, nths,nfm, mth,jmax1, sdsc,copmat,
     $     wrktr1,wrktr2,wrki,outmod )
c
c.....Now, the fields on the plasma from the wall:-
c
      call chinrm ( cwptr, nths,nfm, mth,jmax1, sdsc,copmat,
     $     wrktr1,wrktr2,wrki,outmod )
      call chinrm ( cwpti, nths,nfm, mth,jmax1, sdsc,copmat,
     $     wrktr1,wrktr2,wrki,outmod )
      call chinrm ( bnwpr, nths,nfm, mth,jmax1, sdsc,copmat,
     $     wrktr1,wrktr2,wrki,outmod )
      call chinrm ( bnwpi, nths,nfm, mth,jmax1, sdsc,copmat,
     $     wrktr1,wrktr2,wrki,outmod )
c
c***********************  STEP X    *************************
c
c                    Vacuum Delta-W with the (l-nq)'s
c     
      write ( outmod, 1395 ) 
 1395 format (/,1x,"STEP X. Vacuum Delta-W with (l-nq)." )
c     
      call vacdw ( uvp0, jmax1, vcmtp0 )
      call vacdw ( uvw0, jmax1, vcmtw0 )
      call vacdw ( uvpw0, jmax1, vacmts )
c     
c***********************   STEP X   ******************************
c
c                        Compare Vacuum matrices.
c     
c$$$      write ( outmod, '(a)' ) char(12)
c$$$      write ( outmod, 140 ) 
c$$$ 140  format (/,1x,"STEP X. COMPARE MATRICES." )
c$$$c
c$$$      if ( checks ) then
c$$$      call matwrtn ( vacmatu,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
c$$$     $        "Vacmat_u", outmod,iotty )
c$$$c     
c$$$      call compv ( vacmatu,vacmts, nfm,jmax1, dvacmat )
c$$$c     
c$$$      call vacasym ( vacmatu, nfm,jmax1,"Vacmatu", outmod,iotty )
c$$$
c$$$      call matwrtn ( dvacmat,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
c$$$     $     "norm. difference", outmod,iotty )
c$$$      end if
c
c..........Now revert to original lmin lmax
c
      ldu = lnsav - lmin(1)
      ldl = lnsav - lmin(1)
      mdu = 0
c   
      jxsav = lxsav - lnsav + 1
c
      call shiftmat ( cplar, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( cplai, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( cppr, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( cppi, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( cwalr, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( cwali, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( cwallr, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( cwalli, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( cwwtr, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( cwwti, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( cpwtr, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( cpwti, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( cwptr, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( cwpti, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( bxpwr, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( bxpwi, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( bzpwr, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( bzpwi, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( bnpwr, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( bnpwi, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( bnwpr, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( bnwpi, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( sinlt, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( coslt, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( snlth, nths,nfm, nths,jxsav, ldl,mdu )
      call shiftmat ( cslth, nths,nfm, nths,jxsav, ldl,mdu )
c
      call shiftmat ( vcmtp0, nfm,nfm, jxsav,jxsav, ldl,ldu )
      call shiftmat ( vcmtw0, nfm,nfm, jxsav,jxsav, ldl,ldu )
      call shiftmat ( vacmts, nfm,nfm, jxsav,jxsav, ldl,ldu )
c
      lmax(1) = lxsav
      lmin(1) = lnsav
      jmax1 = lmax(1) - lmin(1) + 1
      lrnge = jmax1
      ln = lmin(1)
      lx = lmax(1)
c     
c...  Abscissa for l values:
      do j1 = 1, jmax1
         lfm(j1) = lmin(1) - 1 + j1
         xlfm(j1) = lfm(j1)
      end do
c..................................................................
c     
      if ( checks ) then
      call matwrtn ( vcmtp0,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $     "Vac-SparkP", outmod,0 )
      call vacasym ( vcmtp0, nfm,jmax1,"Vac-SparkP", outmod,iotty )
      call matwrtn ( vcmtw0,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $     "Vac-SparkW",outmod,0 )
      call vacasym ( vcmtw0, nfm,jmax1,"Vac-SparkW", outmod,iotty )      
      end if
c
c      write ( outmod,'(a)' ) char(12)
c  
      write ( outmod, '(/ 20x, "LFF,FEEDBACK FACTOR = ", i4,
     $     1pe13.5, / )') lff, ff
      write ( iotty,  '(/ 20x, "LFF, FEEDBACK FACTOR = ", i4,
     $     1pe13.5, / )') lff, ff
c
      call matwrtn ( vacmts,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $     "Vac-Spark", outmod,iotty )
      call matwrt9 ( vacmts, nfm,nfm,ln,lx, "Vac-Spark", outmod,iotty )
      call vacasym ( vacmts, nfm,jmax1,"Vac-Spark", outmod,iotty )
      call masym0 ( vacmts, nfm,jmax1,ln,ln,0,0, 1.0,
     $     "Vac-Spark", outmod,iotty )
c      call vacasymi ( vacmts, nfm,jmax1,"Vac-Spark", outmod,iotty )
c
      if ( (lchecke .gt. 0 ) .or. (lff .ne. 0) ) then
         write ( outmod, '(/,20x, "Vac-Spark With l-nq :")' )
            call mateig2a ( vacmts, nfm, jmax1, wrk0, wrki,
     $        ln,lx,ff, 2,lchecke, jobid, outmod )
c$$$            call mateig2 ( vacmts, nfm, jmax1,wrk0, wrki,waa1,wbb1,
c$$$     $        ln,lx,ff, 2,lchecke, wkxx,wkyy, jobid, outmod )
      end if
c
c..... Store the chi in array GRRI to correspond to the direct method.
c
      do i = 1, mth
         do m = 1, jmax1
            grri(i,m)           = cpwtr(i,m)
            grri(i,jmax1+m)     = cpwti(i,m)
            grri(mth+i,m)       = cwwtr(i,m)
            grri(mth+i,jmax1+m) = cwwti(i,m)
         end do
      end do
c     
c$$$      call matwrtn ( cpwtr,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"cpwtr",
c$$$     $     outmod, 0 )
c$$$      call matwrtn ( cpwti,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"cpwti",
c$$$     $     outmod, 0 )
c$$$      call matwrtn ( cwwtr,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"cwwtr",
c$$$     $     outmod, 0 )
c$$$      call matwrtn ( cpwti,nths,nfm,1,ln,mth,jmax1,mdel,jdel,"cpwti",
c$$$     $     outmod, 0 )
c     
      if ( lsymz ) then
c     
c.....symmetrize.
c     
         do  l1 = 1, jmax1
            do  l2 = l1, jmax1
               vacmts(l1,l2) = 0.5 * ( vacmts(l1,l2)+vacmts(l2,l1) )
            end do
         end do
         do  l1 = 1, jmax1
            do  l2 = l1, jmax1
               vacmts(l2,l1) = vacmts(l1,l2)
            end do
         end do
c     
      end if
c
c.....write vacmat to disk for the pest code.
c     
      j12 = 1
c     
      do  j2 = 1, jmax1
         do  j1 = 1, jmax1
            vacpstr(j12) = vacmts(j1,j2)
            vacpsti(j12) = 0.0
            j12 = j12 + 1
         end do
      end do
c     
      call wtopest ( vacpstr, vacpsti, xwal, zwal )
c
      return
      end
c
c......................................................................
      subroutine chinrm ( aio, nths,nfm, mth,jmax1, sdsc,copmat,
     $     wrktr1,wrktr2,wrki,outmod )
c.......................................................................
c
c... Normalizes the matrix aio with the D matrix and outputs to aio.
c
      integer outmod
c
      dimension aio(nths,nfm), sdsc(nfm,nfm), copmat(nfm,nfm),
     $     wrktr1(nfm,nths), wrktr2(nfm,nths), wrki(*)
c
c      First, transpose the chi to wrk.. then Gauss then back to chi..
c
      call mtransr ( aio, nths,nfm, mth,jmax1, wrktr1 )
      call macopy ( sdsc,nfm, copmat,nfm )
      ifail = 0
c      call gelimb ( copmat,nfm,wrktr1,nfm,jmax1,mth,wrktr2,nfm,wrki,
c     $     ifail )      
      call f04aaf ( copmat,nfm,wrktr1,nfm,jmax1,mth,wrktr2,nfm,wrki,
     $     ifail )      
      write ( outmod, '( 1x, "ifail in chinrm  = ", i4 )' ) ifail
      write ( 6, '( 1x, "ifail in chinrm = ", i4 )' ) ifail
      call mtransr ( wrktr2, nfm,nths, jmax1,mth, aio )
c
      return
      end

c..........................................................
      subroutine chiwal ( xobs,zobs,xsce,zsce, bwr,bwi, j1,j2,iopw,iops,
     $     cresr,cresi ) 
c...........................................................
c     
      include 'vacuum1.inc'
      include 'vacuum3.inc'
      include 'vacuum2.inc'
      include 'vacuum5.inc'
c     
c      parameter ( nfmsq = nfm*nfm )
      dimension work1(nths2)
      dimension bwr(nths,nfm),bwi(nths,nfm)
      dimension xobs(*),zobs(*),xsce(*),zsce(*),
     $     cresr(nths,nfm),cresi(nths,nfm),grro(nths2,nfm2)
c
      common / bigv1 / grdgre(nths2,nths2), grwp(nths,nths),
     $     grri(nths2,nfm2)
c
      ier = 0
c     
      jmax = 2*lmax(1) + 1
      jmax1 = lmax(1) - lmin(1) + 1
      q = qa1
      nq = n*q
c     
c.....j1v and j2v are input to vacout and are the size of the matrices.
c     
      j1v = nfm
      j2v = nfm
c
      do 5 i = 1, nths2
         do 4 j = 1, nfm2
            grri(i,j) = zero
 4       continue
 5    continue
c     
      do 10 i = 1, nths2
         do 9 j = 1, nths2
            grdgre(i,j) = zero
 9       continue
 10   continue
c
      ksgn = 2*j2 - 3
c      ksgn = -1

c...Changed 04-04-2004
c$$$      call kernel(xwal,zwal,xwal,zwal,grdgre,grwp,j1,j2,ksgn,iopw,iops)
      call kernel(xobs,zobs,xsce,zsce,grdgre,grwp,j1,j2,ksgn,iopw,iops)

c
c.......Repack GRDGRE from the 2-2 block to the 1-1 block. 
c
      call repack ( grdgre )
c
      if ( checks ) then
      call matwrtn ( grwp,nths,nths,1,1,mth,mth,16,8,
     $        "grwp at IV", outmod, iotty )         
      call matwrtn ( grdgre,nths2,nths2,1,1,mth,mth,16,8,
     $        "grdgre at IV", outmod, iotty )   
      end if
c
      jdc = mth / 8
      write ( iotty,  '(/,"Sum of COLUMNS in GRDGRE:",/)')
      write ( outmod, '(/, "Sum of COLUMNS in GRDGRE:",/)')
      do j1c = 1, mth,jdc 
         sumg = 0.0
         do i = 1, mth
            sumg = sumg + grdgre( mth+i,mth+j1c )
         end do
         write( iotty,  '("jrow, sumg= ",i3,1pe11.3)' ) j1c,sumg
         write( outmod, '("jrow, sumg= ",i3,1pe11.3)' ) j1c,sumg
      end do
c
      write ( iotty,  '(/, "Sum of ROWS in GRDGRE:",/)')
      write ( outmod, '(/, "Sum of ROWS in GRDGRE:",/)')
      do j1c = 1, mth,jdc
         sumg = 0.0
         do i = 1, mth
            sumg = sumg + grdgre( mth+j1c,mth+i )
         end do
         write( iotty,  '("jcol, sumg= ",i3,1pe11.3)' ) j1c,sumg
         write( outmod, '("jcol, sumg= ",i3,1pe11.3)' ) j1c,sumg
      end do
c     
c.......add arbitrary constant to grdgre to make matrix nonsingular.
c     
c      civ = 1.0
c$$$c     
c$$$      write (outmod,8180) civ
c$$$ 8180 format ( /, 5x, "Constant,civ, added to K(i,i) = ", f7.3 )
c$$$c     
c$$$      do  i = 1, mth
c$$$         grdgre(i,i) = grdgre(i,i) + civ
c$$$      end do
c     
      write (outmod,8180) civ
 8180 format ( /, 5x, "Constant,civ, added to K(i,j) = ", f7.3 )
c     
      do 181 i = 1, mth
         do 180 j = 1, mth
            grdgre(i,j) = grdgre(i,j) + civ
 180     continue
 181  continue
c     
c.....fourier analyse source points. real and imag. parts.
c     pack into one-d array. ( ... not for leqt1f or leqt2f... )
c     
c.....the entire cos(lth) matrix is 1-d stored before sin(lth) matrix.
c     
      do 140 l1 = 1, jmax1
         do 135 j = 1, mth
            do 130 i = 1, mth
               grri(i,l1) = grri(i,l1) + bwr(j,l1)*grwp(i,j)
               grro(i,l1) = grro(i,l1) + bwi(j,l1)*grwp(i,j)
               grri(i,jmax1+l1) = grri(i,jmax1+l1) + bwi(j,l1)*grwp(i,j)
 130        continue
 135     continue
 140  continue
c
c$$$      CALL atpoint ( "BEFORE CHECKS ", "j1", j1, "sumg", sumg,
c$$$     $     iotty, outmod )

      if ( checks ) then
      call matwrtn ( grri,nths2,nfm2,1,lmin(1),mth,jmax1,mdel,jdel,
     $     "grrir", outmod,0 )
      call matwrtn ( grro,nths2,nfm2,1,lmin(1),mth,jmax1,mdel,jdel,
     $     "grrii", outmod,0 )
c
      write ( iotty,  '(/,
     $     "Sum of Columns GRRI at Step IV:",/)')
      write ( outmod, '(/,
     $     "Sum of Columns GRRI at Step IV:",/)')

      do jz1 = 1, 2*jmax1
         sumg = 0.0
         do i = 1, mth
            sumg = sumg + grri(i,jz1)
         end do
         write( iotty,  '("jz1, sumg= ",i3,1pe11.3)' ) jz1, sumg
         write( outmod, '("jz1, sumg= ",i3,1pe11.3)' ) jz1, sumg
      end do

      end if

c$$$      CALL atpoint ( "AFTER CHECKS ", "j2", j2, "sumg", sumg,
c$$$     $     iotty, outmod )

c
c.....gaussian elimination.
      
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      call timer ( outmod, iotty, "before leqt1f" )
c     
      lmax2 = 2*jmax1
c     
      ier = 0
c      call gelimb ( grdgre,nths2,grri,nths2,mth,lmax2,grro,nths2,
c     $     work1,ier )
      call f04aaf ( grdgre,nths2,grri,nths2,mth,lmax2,grro,nths2,
     $     work1,ier )
c
      call timer ( outmod, iotty, "after leqt1f" )
c
      write ( outmod,8050 ) ier
      write ( iotty, 8050 ) ier
 8050 format (/, 1x, " ier in gelimb = ", i5 )
c     
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     
c
      if ( checks ) then
      write ( iotty,  '(/,
     $     "Sum of Columns GRRO at Step IV:",/)')
      write ( outmod, '(/,
     $     "Sum of Columns GRRO at Step IV:",/)')
      do jz1 = 1, lmax2
         sumg = 0.0
         do i = 1, mth
            sumg = sumg + grro(i,jz1)
         end do
         write( iotty,  '("jz1, sumg= ",i3,1pe11.3)' ) jz1, sumg
         write( outmod, '("jz1, sumg= ",i3,1pe11.3)' ) jz1, sumg
      end do
      end if
c
c.....Now store wall response matrices.
c
      do 160 l1 = 1, jmax1
         do 150 i = 1, mth
            cresr(i,l1) = grro(i,l1) 
            cresi(i,l1) = grro(i,jmax1+l1)
 150     continue
 160  continue
c
      return
      end
c
c........................................................
      subroutine repack ( grdgre )
c........................................................
c
      include 'vacuum1.inc'
c
      dimension grdgre(nths2,nths2)
c
      do 200 i2 = 1, mth
         do 100 i1 = 1, mth
            grdgre(i1,i2) =  grdgre(mth+i1,mth+i2)
 100     continue
 200  continue
c
      return
      end
c
c.........................................................
      subroutine bwri (xobs,zobs,nobs, xobp,zobp, bxr,bxi,bzr,bzi,
     $     jmax1, bwr,bwi )
c..........................................................
c
c.......Calculates the magnetic field components along some curve.
c.......From both cosine and sine sources. Provides input to GRCHI.
c
      include 'vacuum1.inc'
c     
      dimension xobs(*),zobs(*), xobp(*),zobp(*),
     $     bxr(nths,nfm),bxi(nths,nfm), bzr(nths,nfm),bzi(nths,nfm),
     $     bwr(nths,nfm), bwi(nths,nfm)
c
      do 200 io = 1, nobs+1
         do 100 ll = 1, jmax1
            bwr(io,ll) = xobs(io)*( xobp(io)*bzr(io,ll)
     $           - zobp(io)*bxr(io,ll) )
            bwi(io,ll) = xobs(io)*( xobp(io)*bzi(io,ll)
     $           - zobp(io)*bxi(io,ll) )
 100     continue
 200  continue
c     
      return
      end
c     
c.............................................................
      subroutine grchi(xobs,zobs,nobs, xsce,zsce,xscp,zscp,ns,
     $     isg,creal,cimag,
     $     cdriv,sdriv, ip,ipm, cpwr,cpwi, dchxr,dchxi,dchzr,dchzi)
c.............................................................
c     
c......2nd application of Green's identity. Returns chi, Bx and Bz.
c     
      include 'vacuum1.inc'
      include 'vacuum3.inc'
      include 'vacuum2.inc'
c     
      dimension xobs(*),zobs(*), xsce(*),zsce(*),xscp(*),zscp(*)
      dimension creal(nths,nfm), cimag(nths,nfm)
      dimension cdriv(nths,nfm), sdriv(nths,nfm)
      dimension chrl(5,nfm), chil(5,nfm), wrkr(5),wrki(5)
      dimension cpwr(nths,nfm),cpwi(nths,nfm),
     $     dchxr(nths,nfm),dchxi(nths,nfm),
     $     dchzr(nths,nfm),dchzi(nths,nfm)
c     
      real nq
c
c........   Factpi for dividing into BVAL to make consistent with AVAL
c
      factpi = twopi
c     
      jmax1 = lmax(1) - lmin(1) + 1
      q = qa1
c     
c.......overwrite q on surface
c     
c     q = 1.50
c     
c...................................................
c     
      nq = n * q
      dtpw = twopi / ns
c     
      do 500 io = 1, nobs
         do 450 ixz = 1, 2
            yes1 = 2 - ixz
            yes2 = ixz - 1
c     
c.....Initialize chrl and chil
c     
            do 5 io5 = 1, 5
c.....io5 = 3 gives Chi at xs,zs...
               do 4 ll = 1, jmax1
                  chrl(io5,ll) = 0.0
                  chil(io5,ll) = 0.0
 4             continue
 5          continue
c     
            do 400 io5 = 1, 5
               xs = yes1*(xobs(io) + (io5-3)*delx) + yes2*xobs(io)
               zs = yes2*(zobs(io) + (io5-3)*delz) + yes1*zobs(io)        
c     
               if ( (io5 .eq. 3) .and. (ixz .eq. 2) ) go to 400
c     
               do 200 is = 1, ns
c     
                  xt = xsce(is)
                  zt = zsce(is)
                  xtp = xscp(is)
                  ztp = zscp(is)
c     
                  call green
c
                  bval = bval / factpi
c     
                  do 80 l1 = 1, jmax1
c     
                     ll = lmin(1) - 1 + l1
c
                     if ( ip .eq. 3 ) go to 50
c
                     chrl(io5,l1) = chrl(io5,l1)  + creal(is,l1)*aval
                     chil(io5,l1) = chil(io5,l1)  + cimag(is,l1)*aval
c
 50                  continue
c     
                     if ( ip .eq. 0 ) go to 60
c
c##########change plus to minus bval  #########
                     chrl(io5,l1)  = chrl(io5,l1)  +
     $                    ipm * bval * cdriv(is,l1)
                     chil(io5,l1)  = chil(io5,l1)  +
     $                    ipm * bval * sdriv(is,l1)
c     
 60                  continue
 80               continue
 200           continue
c     
               do 220 l1 = 1, jmax1
                  chrl(io5,l1)  = 0.5 * isg*dtpw * chrl(io5,l1) 
                  chil(io5,l1)  = 0.5 * isg*dtpw * chil(io5,l1) 
 220           continue
c     
 400        continue
            do 430 l1 = 1, jmax1
               do 420 id = 1, 5
                  wrkr(id) = chrl(id,l1)
                  wrki(id) = chil(id,l1)
 420           continue
               if ( yes1 .eq. 1 ) then
                  cpwr(io,l1) = chrl(3,l1)
                  cpwi(io,l1) = chil(3,l1)
                  call diff5 ( wrkr, delx, dchxr(io,l1) )
                  call diff5 ( wrki, delx, dchxi(io,l1) )
               end if
               if ( yes2 .eq. 1 ) then
                  call diff5 ( wrkr, delz, dchzr(io,l1) )
                  call diff5 ( wrki, delz, dchzi(io,l1) )
               end if
 430        continue
 450     continue
 500  continue
c     
      return
      end
c     
c............................................................
      subroutine uvacww ( bwr,bwi, cwalr,cwali, jmax1, uvwwr,uvwwi )
c...........................................................
c     
      include 'vacuum1.inc'
c     
      dimension bwr(nths,nfm), bwi(nths,nfm),
     $     cwalr(nths,nfm), cwali(nths,nfm),
     $     uvwwr(nfm,nfm), uvwwi(nfm,nfm),
     $     uvwwr1(nfm,nfm), uvwwr2(nfm,nfm)
c     
      tpdth = twopi*dth
c
      do 20 j2 = 1, jmax1
         do 10 j1 = 1, jmax1
            uvwwr(j1,j2) = 0.0
            uvwwi(j1,j2) = 0.0
            uvwwr1(j1,j2) = 0.0
            uvwwr2(j1,j2) = 0.0
 10      continue
 20   continue
c     
      do 200 j2 = 1, jmax1
         do 100 j1 = 1, jmax1
            do 50 i = 1,mth
               uvwwr1(j1,j2) = uvwwr1(j1,j2) + bwr(i,j1)*cwalr(i,j2)
               uvwwr2(j1,j2) = uvwwr2(j1,j2) + bwi(i,j1)*cwali(i,j2)
               uvwwi(j1,j2) = uvwwi(j1,j2) + bwr(i,j1)*cwali(i,j2)
     $              - bwi(i,j1)*cwalr(i,j2)
 50         continue
            uvwwr1(j1,j2) = tpdth * uvwwr1(j1,j2)
            uvwwr2(j1,j2) = tpdth * uvwwr2(j1,j2)
            uvwwr(j1,j2) = tpdth * uvwwr(j1,j2)
            uvwwi(j1,j2) = tpdth * uvwwi(j1,j2)
            uvwwr(j1,j2) = uvwwr1(j1,j2) + uvwwr2(j1,j2)
 100     continue
 200  continue
c
      mdel = 16
      jdel = 8
      ln = lmin(1)
      lx = lmax(1)
      if ( checks ) then
      call matwrtn ( uvwwr1,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $     "uvwwr1 at IV", outmod,0 )
      call matwrtn ( uvwwr2,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $     "uvwwr2 at IV", outmod,0 )
      end if
c
      return
      end
c......................................................
      subroutine uvacw ( jm, uvpr,uvpi, uvwr,uvwi, uvpwr,uvpwi )
c......................................................
c
c*****   Separate out Plasma and Wall contributions for now!!!!! *******
c     
      include 'vacuum1.inc'
      include 'vacuum5.inc'
c     
      dimension uvpr(nfm,nfm),uvpi(nfm,nfm),
     $     uvwr(nfm,nfm),uvwi(nfm,nfm),
     $     uvpwr(nfm,nfm),uvpwi(nfm,nfm)
c
      tpdth = twopi*dth
c     
      do 20 jb = 1, jm
         do 10 ja = 1, jm
            uvpr(ja,jb) = 0.0
            uvpi(ja,jb) = 0.0
            uvwr(ja,jb) = 0.0
            uvwi(ja,jb) = 0.0
            uvpwr(ja,jb) = 0.0
            uvpwi(ja,jb) = 0.0
 10      continue
 20   continue
c     
      do 420 jb = 1, jm
         do 400 ja = 1, jm
            do 200 i = 1, mth
               uvpr(jb,ja) = uvpr(jb,ja) +
     $              cplar(i,ja) * cslth(i,jb) + 
     $              cplai(i,ja) * snlth(i,jb)
               uvpi(jb,ja) = uvpi(jb,ja) -
     $              cplai(i,ja) * cslth(i,jb) +
     $              cplar(i,ja) * snlth(i,jb)
               uvwr(jb,ja) = uvwr(jb,ja) +
     $              cpwr(i,ja) * cslth(i,jb) + 
     $              cpwi(i,ja) * snlth(i,jb)
               uvwi(jb,ja) = uvwi(jb,ja) -
     $              cpwi(i,ja) * cslth(i,jb) +
     $              cpwr(i,ja) * snlth(i,jb)
               uvpwr(jb,ja) = uvpwr(jb,ja) +
     $              ( cplar(i,ja) + cpwr(i,ja) ) * cslth(i,jb) + 
     $              ( cplai(i,ja) + cpwi(i,ja) ) * snlth(i,jb)
               uvpwi(jb,ja) = uvpwi(jb,ja) -
     $              ( cplai(i,ja) + cpwi(i,ja) ) * cslth(i,jb) +
     $              ( cplar(i,ja) + cpwr(i,ja) ) * snlth(i,jb)
 200        continue
            uvpr(jb,ja) = tpdth * uvpr(jb,ja)
            uvpi(jb,ja) = tpdth * uvpi(jb,ja)
            uvwr(jb,ja) = tpdth * uvwr(jb,ja)
            uvwi(jb,ja) = tpdth * uvwi(jb,ja)
            uvpwr(jb,ja) = tpdth * uvpwr(jb,ja)
            uvpwi(jb,ja) = tpdth * uvpwi(jb,ja)
 400     continue
 420  continue
c     
      return
      end
c
c......................................................
      subroutine vacdw ( vac0, jm, vac )
c......................................................
c
      include 'vacuum1.inc'
      dimension vac0(nfm,nfm), vac(nfm,nfm)
c
      anq = n * qa1
c
      do 20 jb = 1, jm
         lb = lmin(1) + jb - 1
         blnq = lb - anq
         do 10 ja = 1, jm
            la = lmin(1) + ja - 1
            alnq = la - anq
            vac(ja,jb) = alnq*blnq* vac0(jb,ja)
 10      continue
 20   continue
c
      return
      end
c     
c.....................................................
      subroutine compv ( vacina,vacinb, nfm,nj, dvac )
c.....................................................
c     
      dimension vacina(nfm,nfm), vacinb(nfm,nfm),dvac(nfm,nfm)
c     
      do 200 jb = 1, nj
         do 100 ja = 1, nj
            dnom =  vacinb(ja,jb) + vacina(ja,jb)
            dvac(ja,jb) = 1.0
            if ( dnom .eq. 0 ) go to 100
            dvac(ja,jb) = ( vacinb(ja,jb) - vacina(ja,jb) ) / dnom
 100     continue
 200  continue
c     
      return
      end
c
c.......................................................
      subroutine dmats ( bc,bs, dcc,dss, sdsc )
c.......................................................
c     
      include  'vacuum1.inc'
      include  'vacuum5.inc'
c     
      dimension bc(nths,nfm), bs(nths,nfm),
     $     dcc(nfm,nfm),dss(nfm,nfm)
      dimension dcs(nfm,nfm),dsc(nfm,nfm),ddsc(nfm,nfm),sdsc(nfm,nfm)
c     
      jmax1 = lmax(1) - lmin(1) + 1
      ln = lmin(1)
      dsign = -1.0
c     
      do 20 la = 1, jmax1
         do 18 lb = 1, jmax1
            dcc(la,lb) = 0.0
            dss(la,lb) = 0.0
            dcs(la,lb) = 0.0
            dsc(la,lb) = 0.0
 18      continue
 20   continue
c     
      do 200 la = 1, jmax1
         do 180 lb = 1, jmax1
            do 150 i = 1, mth
               dcc(lb,la) = dcc(lb,la) + bc(i,la) * cslth(i,lb)
               dss(lb,la) = dss(lb,la) + bs(i,la) * snlth(i,lb)
               dcs(lb,la) = dcs(lb,la) + bc(i,la) * snlth(i,lb)
               dsc(lb,la) = dsc(lb,la) + bs(i,la) * cslth(i,lb)
 150        continue
            dcc(lb,la) = dcc(lb,la)/mth
            dss(lb,la) = dss(lb,la)/mth
            dcs(lb,la) = dcs(lb,la)/mth
            dsc(lb,la) = dsc(lb,la)/mth
            ddsc(lb,la) = dsc(lb,la) - dcs(lb,la)
            sdsc(lb,la) = dsign * ( dcc(lb,la) + dss(lb,la) )
            if ( la .eq. lb ) sdsc(lb,la) = 1.0 + sdsc(lb,la)
 180     continue
 200  continue
c     
      if ( checks ) then
      call matwrtn ( dcs,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,"DCS",
     $     outmod,0 )
      call matwrtn ( dsc,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,"DSC",
     $     outmod,0 )
      call matwrtn ( ddsc,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,"DDSC",
     $     outmod,0 )
      end if
c     
      return
      end
c
c.........................................................
      subroutine bcomps ( chir, nd, dchl )
c.........................................................
c
      include 'vacuum1.inc'
      include 'vacuum3.inc'
      include 'vacuum5.inc'
c
      dimension chir(*), wrk(5)
c     
         do 418 id = 1, nd
            wrk(id) = chir(id)
 418     continue
c
         call diff5 ( wrk, delx, df5 )
         dchl = df5
c
      return
      end
c
