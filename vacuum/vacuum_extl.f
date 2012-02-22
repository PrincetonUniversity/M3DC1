c...................................................
c
c     3.1   conducting wall specification.
c.......................................................................
      subroutine wwall(nqnqnq,xwal,zwal)
c.......................................................................
c      conducting wall specification
c
c        uses a, b, aw, bw from namelist in modin.
c
c        if a less than -100.0 the wall is ellipsoidal with
c        horizontal radius hrad equal to xplamax(1.0+aw), and
c        vertical radius vrad equal to zplamax(1.0+bw).
c
c        if a greater than equal 10, then cond. wall at infinity
c        if a less than equal -10. then wall at a fixed distance
c        given by b.
c        if not then set wall by the expn.
c         xw = a + xma + aw*cos(theta + b*sin(theta))
c         zw = zma + bw*sin(theta)
c
c        all lengths are scaled to the larger of the max. vertical
c        distance ,max. horizontal distance.
c
c      in the third option 'a' defines the shift of the wall shape 'axis',
c        'b' the triangularity, 'aw' the minor radius and ' bw'the major radius
c
c...........................................................
c
      include 'vacuum1.inc'
      include 'vacuum3.inc'
      include 'vacuum2.inc'
      logical infwal, lfix, insect
      dimension xwal(*), zwal(*)
      dimension iop(2),xpp(nths),zpp(nths),ww1(nths),ww2(nths),
     .    ww3(nths),thet(nths),tabx(3),tabz(3)
c
      REAL, DIMENSION(:), ALLOCATABLE :: xwal0, zwal0
      INTEGER :: mthin
      CHARACTER*(80) tit_walin

      REAL, DIMENSION(:), ALLOCATABLE :: vecinw_x, vecinw_z

      character*(80) stringv(5)
c
      data iplt/0/
c
      call atpoint ( "WWALL","ishape", ishape,"xinf(i)", xinf(1),
     $     iotty, outmod )
c
      xpp(1) = 1.
      zpp(1) = 1.
      ww3(1) = 1.
      tabx(1) = 1.
      tabz(1) = 1.
      awsave = aw
      bwsave = bw
      insect = .false.
c      infwal = .false.
c
c.....define wall for spheromak type plasma.
c
      isph = 0
      if ( a .ge. -100.0 ) go to 9
      isph = 1
      ishape = -10
c
c.....find maximum x and z values for the plasma surface.
c.....and compute horizontal and vertical radius of wall.
c
      call bounds(xinf,zinf,1,mth,xmnp,xmxp,zmnp,zmxp)
      xmin = xmnp
      xmax = xmxp
      zmin = zmnp
      zmax = zmxp
c
      plrad = 0.5 * ( xmxp - xmnp )
      xmaj = 0.5 * ( xmxp + xmnp )
      zmid = 0.5 * ( zmnp + zmxp )
c
      hrad = xmax + aw*(xmax-xmaj)
      vrad = zmax + bw*(zmax-zmid)
c
      do 8 i = 1, mth1
c
      xi = xinf(i) - xmaj
      zeta = zinf(i) - zmid
      bbb = (xi*vrad)**2 + (zeta*hrad)**2
      ccc = -xmaj*vrad*xi + hrad*sqrt( bbb - (zeta*xmaj)**2 )
c
      xwal(i) = xmaj + xi*vrad*ccc/bbb
      zwal(i) = zmid + zeta*vrad*ccc/bbb
c
    8 continue
c
      go to 145
c
c.....end of ellipsoidal wall.
c
    9 continue
c
c     dw = b
      if(a .gt. -10.)lfix=.true.
      if( a .lt. 10. ) go to 10
      infwal = .true.
      go to 2000
   10 continue
      xshift = a
c
c$$$      go to 20
c$$$c
c$$$      if( .not. lfix) go to 15
c$$$c
c$$$c      check bw
c$$$c
c$$$      if( bw .gt. 1.0 ) go to 20
c$$$      write(iotty,1000)bw
c$$$      write(outpest,1000)bw
c$$$      write(outmod,1000)bw
c$$$ 1000 format(" error bw is less than 1 =",e12.5," wall will intersect
c$$$     .  the plasma"/)
c$$$      call errmes(outpest,'vacdat')
c$$$c
c$$$c      check b. should be greater than 0.0
c$$$c
c$$$   15 continue
c$$$      if( b .ge. 0.0) go to 20
c$$$      write(iotty,1100)b
c$$$      write(outpest,1100)b
c$$$      write(outmod,1100)b
c$$$      call errmes(outpest,'vacdat')
c$$$ 1100 format(" **error** b is less than zero=",e12.5)
c
   20 continue
c
      mthalf = mth2 / 2
c
      xmax = MAXVAL ( xinf(1:mth) )
      zmax = MAXVAL ( zinf(1:mth) )
      xmin = MINVAL ( xinf(1:mth) )
      zmin = MINVAL ( zinf(1:mth) )

c$$$      zmin = huge(0.0D0)
c$$$      zmax = -huge(0.0D0)
c$$$      xmin = huge(0.0D0)
c$$$      xmax = -huge(0.0D0)
c$$$      do 40 i = 1, mth
c$$$      if(xmax .lt. xinf(i))xmax = xinf(i)
c$$$      if(xmin .gt. xinf(i))xmin = xinf(i)
c$$$      if(zmax .lt. zinf(i)) zmax = zinf(i)
c$$$      if(zmin .gt. zinf(i))zmin = zinf(i)
c$$$   40 continue
c
      plrad = 0.5 * ( xmax-xmin )
      xmaj = 0.5 * ( xmax + xmin )
      zmid = 0.5 * ( zmax + zmin )
      zrad = 0.5 * ( zmax - zmin )
c
      scale =  ( zmax - zmin )
      if(((xmax-xmin)/2.0).gt.scale) scale = (xmax-xmin)/2.0
c
c......   set sacle = 1.0 temporarily.....
c
      scale = 1.0
c
      aw = aw * scale
      bw = bw * scale
      delta1 = dw * (xinf(1) - xma)
c
      if ( ishape .ne. 2 ) go to 295
c
c  Wall is Confocal ellipse at horizontal distance of a.
c
      zh = sqrt ( abs(zrad**2 - plrad**2) )   ! h metric
      zah = a / zh
      zph = plrad / zh
      zmup = 0.5*alog ((zrad+plrad)/(zrad-plrad))  ! mu-plas
      zmuw = alog ( zah + sqrt(zah**2 + 1) )  ! mu-wall
      zxmup = exp(zmup)
      zxmuw = exp(zmuw)
      zbwal = zh * cosh ( zmuw )              ! Major radius of wall
      bw = zbwal / a                          ! Elongation of wall
      do i = 1, mth2
         the = (i-1) * dth
         xwal(i) =   xmaj + a*cos( the )
         zwal(i) = - bw * a*sin( the )
      end do
c
      write ( outmod, '(/,"Confocal Ellipse:"/,
     $     "mup, expmup = ", 1p2e13.5,/,
     $     "muw, expmuw = ", 1p2e13.5,/)') zmup, zxmup, zmuw, zxmuw
      write ( iotty,  '(/,"Confocal Ellipse:"/,
     $     "mup, expmup = ", 1p2e13.5,/,
     $     "muw, expmuw = ", 1p2e13.5,/)') zmup, zxmup, zmuw, zxmuw
c
 295  continue
c
      IF ( ishape .NE. 3 ) GO TO 300

      call atpoint ( "nstxhel-bef","ishape", ishape,"xwal(1)", xwal(1),
     $     iotty, outmod )

      CALL nstx_shell ( xwal, zwal, mth, lsymwv, ioshel,
     $     outmod, iotty )

      call atpoint ( "nstxhel-aft","ishape", ishape,"zwal(1)", zwal(1),
     $     iotty, outmod )

      WRITE ( outmod, '(/,10x,
     $     "Xwal, Zwal, after call to nstx_shell" )' )
      DO i = 1, mth2, 4
         WRITE ( outmod, '(i4, 2es11.3)' ) i, xwal(i), zwal(i)
      END DO

  300 continue
c
      if ( ishape .ne. 4 ) go to 320
c
c.....  dee-shaped wall independent of plasma position.
c.....Added parameters for NSTX added 5/10/95
c   +tw squares up the outer corners. +aw rounds the inner corners.
c
      wcentr = cw
c
      do i = 1, mth2
         the0 = (i-1) * dth
c         the = the0 + aw*sin(2.0*the0)
         the = the0
         sn2th = sin(2.0*the)
         xwal(i) =   cw + rw*cos( the + dw*sin(the) )
         zwal(i) = - bw * rw*sin( the + tw*sn2th ) - aw*sn2th
      end do
c
 320  continue
c
      IF ( ishape /= 45 ) GO TO 345

c...  Rectangular shell with smoothed corners.
c..   NOTE: users must ensure that there is a grid point at the 
c     corners of the reference rectangle (i.e., before the smoothing).
c     Should initialize the following four input variables in NAMELIST later.

c$$$      xcen = 3.0 
c$$$      zcen = 0.0 
c$$$      xrad = 1.5 
c$$$      zrad = 1.5 

      xcen = cw
      zcen = dw
      xrad = aw
      zrad = bw

c...  Define edges 
      xtl = xcen - xrad 
      xtr = xcen + xrad 
      xbr = xcen + xrad 
      xbl = xcen - xrad 

      ztl = zcen + zrad 
      ztr = zcen + zrad 
      zbr = zcen - zrad 
      zbl = zcen - zrad 

      acir = 4.0 * ( xrad + zrad ) 
      dl = acir / mth 

      ALLOCATE ( xwal0(mth+5), zwal0(mth+5) )
      
      nrb = (zcen-zbr) / dl + 0.1
      xwal0(1:nrb+1) = (/ (xbr, i=1,nrb+1) /) 
      zwal0(1:nrb+1) = 
     $     (/ (zcen-(i-1)*dl, i = 1,nrb+1) /)

      nb = (xbr-xbl) / dl + 0.1
      xwal0(nrb+1:nrb+nb+1) = 
     $     (/ (xbr-(i-1)*dl, i = 1, nb+1) /)
      zwal0(nrb+1:nrb+nb+1) = (/ (zbr, i = nrb+1,nrb+nb+1) /)

      nl = (ztl-zbl) / dl + 0.1 
      xwal0(nrb+nb+1:nrb+nb+nl+1) =
     $     (/ (xbl, i = nrb+nb+1, nrb+nb+nl+1) /)
      zwal0(nrb+nb+1:nrb+nb+nl+1) = 
     $     (/ (zbl+(i-1)*dl, i = 1, nl+1) /)

      nt = (xtr-xtl) / dl + 0.1
      xwal0(nrb+nb+nl+1:nrb+nb+nl+nt+1) = 
     $     (/ (xtl+(i-1)*dl, i = 1, nt+1) /)
      zwal0(nrb+nb+nl+1:nrb+nb+nl+nt+1) =
     $     (/ (ztl, i = nrb+nb+nl+1, nrb+nb+nl+nt+1) /)

      nrt = (xtr-xcen) / dl + 0.1
      xwal0(nrb+nb+nl+nt+1:nrb+nb+nl+nt+nrt+1) =
     $     (/ (xtr, i = nrb+nb+nl+nt+1, nrb+nb+nl+nt+nrt+1) /)
      zwal0(nrb+nb+nl+nt+1:nrb+nb+nl+nt+nrt+1) = 
     $     (/ (ztr-(i-1)*dl, i = 1, nrt+1) /)

c.. Check number of points.
      nrect = nrb+nb+nl+nt+nrt 
      WRITE ( iotty, '(/, 10x, "nrect, mth = ", 2i5 )' )
     $     nrect, mth
      WRITE ( outmod, '(/, 10x, "nrect, mth = ", 2i5 )' )
     $     nrect, mth

      CALL  msc_smooths1y ( xwal0, xwal, mth, 8 )
      CALL  msc_smooths1y ( zwal0, zwal, mth, 8 )
 
      xwal(mth1) = xwal(1)
      zwal(mth1) = zwal(1)
      xwal(mth2) = xwal(2)
      zwal(mth2) = zwal(2)

      DEALLOCATE ( xwal0, zwal0 )

 345  CONTINUE

      if ( ishape .ne. 5 ) go to 340
c
c.....   dee-shaped wall.
c
      wcentr = xmaj + cw*plrad
c
      do 330 i = 1, mth2
c
         the0 = (i-1) * dth
c         the = the0 + aw*sin(2.0*the0)
          the = the0
        sn2th = sin(2.0*the)
         xwal(i) = xmaj + cw*plrad +
     $        plrad*(1.0+a-cw)*cos(the+dw*sin(the))
         zwal(i) = - bw*plrad*(1.0+a-cw) * sin( the + tw*sn2th ) -
     $        aw*plrad*sn2th
  330 continue
c
  340 continue
c
      if ( ishape .ne. 6 ) go to 338
c
c.....  wall at constant distance from plasma.
c
      wcentr = xmaj
c
      do i = 2, mth1
         alph = atan2m ( xinf(i+1)-xinf(i-1), zinf(i-1)-zinf(i+1) )
         xwal(i) = xinf(i) + a*plrad * cos(alph)
         zwal(i) = zinf(i) + a*plrad * sin(alph)
      end do
c
      xwal(1) = xwal(mth1)
      zwal(1) = zwal(mth1)
      xwal(mth2) = xwal(2)
      zwal(mth2) = zwal(2)
c
  338 continue
c
      if ( ishape .ne. 7 ) go to 348
c
c..... enclosing bean shaped walls.
c
      cwr = cw * pye / 180.0
c
      do i = 1, mth2
         the0 = (i-1)*dth
c     the  = the0 - bw*sin(the0)
         the = the0
         rho = aw * ( 1.0 + bw*cos(the) )
         the2 = cwr * sin(the)
         xofsw = xmax + a*plrad - aw*(1.0+bw)
         xwal(i) = xofsw + rho*cos(the2)
         zwal(i) = - b*rho*sin(the2)
      end do
c
  348 continue
c
      if ( ishape .ne. 8 ) go to 349
c
      call d3dwall ( xwal, zwal, mth, outmod, iotty )
c
 349  continue

      IF ( ishape /= 81 ) go to 3491

      CALL d3d_shell ( xwal, zwal, mth, lsymwv, ioshel, outmod, iotty )

      WRITE ( outmod, '(/,10x, "Xwal, Zwal, after call to d3dwall" )' )
      DO i = 1, mth2, 4
         WRITE ( outmod, '(i4, 2es11.3)' ) i, xwal(i), zwal(i)
      END DO

 3491 CONTINUE
c
      if ( ishape .ne. 9 ) go to 350
c
c... Wall shape read in as xwin, zwin:
      do i = 1, nths
         xwal(i) = xwin(i)
         zwal(i) = zwin(i)
      end do

 350  continue

c... ITER wall
      IF ( ishape .NE. 11 ) GO TO 355

      call atpoint ( "ITERshel-bef","ishape", ishape,"xwal(1)", xwal(1),
     $     iotty, outmod )

      CALL iter_sym_shell ( xwal, zwal, mth, ioshel, outmod, iotty )

      call atpoint ( "ITERshel-aft","ishape", ishape,"zwal(1)", zwal(1),
     $     iotty, outmod )

      WRITE ( outmod, '(/,10x,
     $     "Xwal, Zwal, after call to iter_sym_shell" )' )
      DO i = 1, mth2, 4
         WRITE ( outmod, '(i4, 2es11.3)' ) i, xwal(i), zwal(i)
      END DO

c... Expand this wall by awx to simulate blanket. And shift the wall
c    vertically by awv.

      do i = 2, mth1
         alph = atan2m ( xwal(i+1)-xwal(i-1), zwal(i-1)-zwal(i+1) )
         ww1(i) = xwal(i) + awx*plrad * cos(alph)
         ww2(i) = zwal(i) + awx*plrad * sin(alph)
      end do
c
      xwal(2:mth1) = ww1(2:mth1)
      zwal(2:mth1) = ww2(2:mth1)
      xwal(1) = xwal(mth1)
      zwal(1) = zwal(mth1)
      xwal(mth2) = xwal(2)
      zwal(mth2) = zwal(2)
c
 355  continue

c.. READS wall data from disk, WALL_IN.

      IF ( ishape /= 15 ) GO TO 365

      OPEN ( UNIT = 365, FILE='wall_in' )
      
      WRITE ( OUTMOD, '(/"***** reading XWIN, ZWIN from WALL_IN",/)' )
      WRITE ( IOTTY,  '(/"***** reading XWIN, ZWIN from WALL_IN",/)' )

      READ(365,*) mthin

c-----------------------------------------------------------------------
c     read X,Z arrays.
c-----------------------------------------------------------------------

      ALLOCATE ( vecinw_x(mthin+5), vecinw_z(mthin+5) )

      READ( 365,'(A)') tit_walin

      WRITE ( OUTMOD, '("mthin = ", I5, A ,/)' ) mthin, tit_walin
      WRITE ( IOTTY,  '("mthin = ", I5, A ,/)' ) mthin, tit_walin

      DO i = 1, mthin
         READ(365, * ) vecinw_x(i), vecinw_z(i)
      END DO

      dx0 = 0.0
      dx1 = 0.0
      CALL transdx ( vecinw_x,mthin, xwal,mth, dx0,dx1 )
      CALL transdx ( vecinw_z,mthin, zwal,mth, dx0,dx1 )
      xwal(mth1) = xwal(1)
      xwal(mth2) = xwal(2)
      zwal(mth1) = zwal(1)
      zwal(mth2) = zwal(2)

      DEALLOCATE ( vecinw_x, vecinw_z )

      CLOSE ( UNIT = 365 )

 365  CONTINUE
c
      if ( ishape .ne. 110 ) go to 400
c
c.... wall is finite d-shaped coil.
c
      do i = 1, mth2
c
         the = (i-1) * dth
         plrad = 0.5 * ( xmax - xmin )
         xwal(i) = xmax + plrad * ( a + aw - aw*cos(the + dw*sin(the)) )
         zwal(i) = - plrad * bw * sin(the)
      end do
c
  400 continue
c
      if ( ishape .ne. 120 ) go to 500
c
c.....  wall is finite bean shaped conductor on right.
c.....      offset by xmaj+xofs+wd*plrad.   a*plrad away from plasma edge.
c.....      half width of coil is aw*plrad.
c.....     subtending angle is b.  elongation factor is bw.
c.....      offset from xmaj when a=0 is cw
c
      plrad = 0.5 * ( xmax-xmin )
      xmaj = 0.5 * ( xmax + xmin )
      a0 = plrad * ( 1.0 + aw - cw + a )
      brad = b * pye / 180.0
c
      do i = 1, mth2
         the0 = (i-1) * dth
c     the = the0 - aw*plrad * sin(2.0*the0) / a0
         the = the0
         rho = a0 - aw*plrad*cos(the)
         the2 = brad * sin(the)
         xwal(i) = xmaj + cw*plrad + rho * cos(the2)
         zwal(i) = - bw * rho * sin(the2)
      end do
c
  500 continue
c
      if ( ishape .ne. 130 ) go to 600
c
c.....  wall is finite bean shaped conductor on left.
c.....      offset by xmaj+xofs+wd*plrad.   a*plrad away from plasma edge.
c.....      half width of coil is aw*plrad.
c.....     subtending angle is b.  elongation factor is bw.
c.....      offset from xmaj when a=0 is cw
c
      plrad = 0.5 * ( xmax-xmin )
      xmaj = 0.5 * ( xmax + xmin )
      a0 = plrad * ( 1.0 + aw - cw + a )
      brad = b * pye / 180.0
c
      do 550 i = 1, mth2
c
      the0 = (i-1) * dth
c      the = the0 - aw*plrad * sin(2.0*the0) / a0
          the = the0
      rho = a0 + aw*plrad*cos(the)
      the2 = brad * sin(the)
      xwal(i) = xmaj + cw*plrad - rho * cos(the2)
      zwal(i) = - bw * rho * sin(the2)
c
  550 continue
c
  600 continue
c
      if ( ishape .ne. 210 ) go to 700
c
c.....  wall is solid plate on right normalized to plasma radius.
c.....      offset by xmaj+xofs+wd*plrad.   a*plrad away from plasma edge.
c.....      half width of coil is aw*plrad.
c.....      subtending angle is b.  elongation factor is bw.
c.....      offset from xmaj when a=0 is cw
c.....      triangularity is dw.   sharpness of tips is tw
c.....      size of bulginess along X-axis is abulg*plrad.
c.....      half-angle of bulge in degrees is bbulg.
c.....      sharpness of bulge corners is tbulg.
c
      plrad = 0.5 * ( xmax-xmin )
      xmaj = 0.5 * ( xmax + xmin )
      a0 = plrad * ( 1.0 + aw - cw + a )
      a0b = (a0 + plrad*aw)*bw
      brad0 = b * pye / 180.0
c      call adjustb ( brad0, brad, a, bw, cw, dw, xmaj, plrad,
c     $     ishape )
      brad = brad0
      blgrad0 = bbulg * pye / 180.0
      wcentr = xmaj + cw*plrad
c
      call adjustb ( blgrad0, blgrado, a, bw, cw, dw, xmaj, plrad,
     $     ishape )
      dthb = ( 2.0*aw*plrad / a0b )
     $        * ( 1.0 - sin(blgrado) ) / cos(blgrado)
      blgrad0 = blgrad0 - dthb
      call adjustb ( blgrad0, blgradi, a, bw, cw, dw, xmaj, plrad,
     $     ishape )
c
      do 650 i = 1, mth2
c
      the0 = (i-1) * dth
      if ( the0 .gt. 0.5*pye .and. the0 .lt. 1.5*pye ) then
         thbulg = blgrado
      else
         thbulg = blgradi
      end if
c
      cost2b = cos(2*thbulg)
c      the = the0 - aw*plrad * sin(2.0*the0) / a0
      the = the0
      cost = cos(the)
      ferm = +1.0 - 2.0 / ( exp(cost/tw) + 1.0 )
      rho = a0 - aw*plrad*ferm
      the2 = brad * sin(the)
      cost2 = cos(2.0*the2)
      fermb = 1.0 / ( exp( (cost2b - cost2)/tbulg ) + 1.0 )
      bulge = abulg*plrad*fermb
      xwal(i) = xmaj + cw*plrad + rho * cos(the2+dw*sin(the2))
     $     + bulge
      zwal(i) = - bw * rho * sin(the2)
c
  650 continue
c
  700 continue
c
      if ( ishape .ne. 240 ) go to 800
c
c.....  wall is solid plate  on left.
c.....      offset by xmaj+xofs+wd*plrad.   a*plrad away from plasma edge.
c.....      half width of coil is aw*plrad.
c.....     subtending angle is b.  elongation factor is bw.
c.....      offset from xmaj when a=0 is cw
c.....      triang. is dw.  sharpness of tips is tw.
c
      plrad = 0.5 * ( xmax-xmin )
      xmaj = 0.5 * ( xmax + xmin )
      a0 = plrad * ( 1.0 + aw - cw + a )
      brad = b * pye / 180.0
c
      wcentr = xmaj + cw*plrad
c
      do 750 i = 1, mth2
c
      the0 = (i-1) * dth
c      the = the0 - aw*plrad * sin(2.0*the0) / a0
          the = the0
      cost = cos(the)
      ferm = +1.0 - 2.0 / ( exp(cost/tw) + 1.0 )
      rho = a0 + aw*plrad*ferm
      the2 = brad * sin(the)
      xwal(i) = xmaj + cw*plrad - rho * cos(the2-dw*sin(the2))
      zwal(i) = - bw * rho * sin(the2)
c
  750 continue
c
  800 continue
c
      if ( ishape .ne. 310 ) go to 1700
c
c.....  wall is solid plate  on right independent of plasma position.
c.....      center of wall is at cw.
c.....      half width of coil is aw
c.....      radius of inner side of the wall is a.
c.....      subtending angle is b.  elongation factor is bw.
c.....      triangularity is dw.   sharpness of tips is tw.
c.....      size of bulginess along the X-axis is abulg.
c.....      half-angle of bulge in degrees is bbulg.
c.....      sharpness of bulge corners is tbulg.
c
      a0 = a + aw
      a0b = (a0 + aw)*bw
      brad0 = b * pye / 180.0
c      call adjustb ( brad0, brad, a, bw, cw, dw, xmaj, plrad,
c     $     ishape )
      brad = brad0
      blgrad0 = bbulg * pye / 180.0
      wcentr = cw
c
      call adjustb ( blgrad0, blgrado, a, bw, cw, dw, xmaj, plrad,
     $     ishape )
      dthb = ( 2.0*aw / a0b )
     $        * ( 1.0 - sin(blgrado) ) / cos(blgrado)
      blgrad0 = blgrad0 - dthb
      call adjustb ( blgrad0, blgradi, a, bw, cw, dw, xmaj, plrad,
     $     ishape )
c
      do 1650 i = 1, mth2
c
      the0 = (i-1) * dth
      if ( the0 .gt. 0.5*pye .and. the0 .lt. 1.5*pye ) then
         thbulg = blgrado
      else
         thbulg = blgradi
      end if
c
      cost2b = cos(2.0*thbulg)
c      the = the0 - aw * sin(2.0*the0) / a0
         the = the0
      cost = cos(the)
      ferm = +1.0 - 2.0 / ( exp(cost/tw) + 1.0 )
      rho = a0 - aw*ferm
      the2 = brad * sin(the)
      cost2 = cos(2.0*the2)
      fermb = 1.0 / ( exp( (cost2b - cost2)/tbulg ) + 1.0 )
      bulge = abulg*fermb
      xwal(i) = cw + rho * cos(the2+dw*sin(the2)) + bulge
      zwal(i) = - bw * rho * sin(the2)
c
 1650 continue
c
 1700 continue
c
      if ( ishape .ne. 340 ) go to 1800
c
c.....  wall is solid plate  on left independent of plasma position.
c.....      center of wall is at cw.
c.....      half width of coil is aw
c....       radius of inner side of the wall is a.
c.....      subtending angle is b.  elongation factor is bw.
c.....      triangularity is dw.   sharpness of tips is tw
c
      a0 = a + aw
      brad = b * pye / 180.0
c
      wcentr = cw
c
      do 1750 i = 1, mth2
c
      the0 = (i-1) * dth
c      the = the0 - aw * sin(2.0*the0) / a0
         the = the0
      cost = cos(the)
      ferm = +1.0 - 2.0 / ( exp(cost/tw) + 1.0 )
      rho = a0 + aw*ferm
      the2 = brad * sin(the)
      xwal(i) = cw - rho * cos(the2-dw*sin(the2))
      zwal(i) = - bw * rho * sin(the2)
c
 1750 continue
c
 1800 continue
c
c     write(outmod,1460)(i,xinf(i),zinf(i),xwal(i),zwal(i),i=1,mthalf)
c     if(xwal(1) .lt. xinf(1))insect=.true.
c     if(xwal(mthalf).gt.xinf(mthalf))insect=.true.
      xmx = xma + xshift
c
c.......... bypass this section.  not suitable for beans.
c
c$$$      go to 145
c$$$c
c$$$c&&&&& test change  &&&&&
c$$$c    fit spline to get xw,zw at equal intervals of theta
c$$$c     periodic boundary conditions
c$$$      iop(1) = 4
c$$$      iop(2) = 4
c$$$      do 110 il=mthalf+1,mth1
c$$$      ilm = mth1 - il + 1
c$$$      thet(il) = twopi - thet(ilm)
c$$$  110 continue
c$$$c
c$$$      call spl1d1(mth1,thet,xwal,xpp,iop,1,ww1,ww2,ww3)
c$$$      call spl1d1(mth1,thet,zwal,zpp,iop,1,ww1,ww2,ww3)
c$$$c
c$$$      do 125 i=2,mthalf-1
c$$$      xs = xinf(i)
c$$$      zs = zinf(i)
c$$$      xt = xwal(i)
c$$$      zt = zwal(i)
c$$$      tt = thet(i)
c$$$c      now use a newton-rafson iteration scheme to minimize d(xd)/d(theta)
c$$$c      where xd = (xs-xt)**2 + (zs-zt)**2
c$$$c
c$$$      do 120 k=1,20
c$$$      call spl1d2(mth1,thet,xwal,xpp,1,tt,tabx)
c$$$      call spl1d2(mth1,thet,zwal,zpp,1,tt,tabz)
c$$$      xt = tabx(1)
c$$$      zt = tabz(1)
c$$$      xmx1 = xs - xt
c$$$      zmz1 = zs - zt
c$$$      f = xmx1 * tabx(2) + zmz1 * tabz(2)
c$$$      fp = xmx1*tabx(3) + zmz1*tabz(3) - (tabx(2))**2 - (tabz(2))**2
c$$$      delt = f / fp
c$$$      if(abs(delt).lt.1.0e-4 .and. abs(f) .lt. 1.0e-4) go to 124
c$$$      tt = tt - delt
c$$$  120 continue
c$$$      write(iotty,1430)i,xs,zs,xt,zt,tt,f,fp
c$$$      write(outpest,1430)i,xs,zs,xt,zt,tt,f,fp
c$$$      call errmes(outpest,'vacdat')
c$$$ 1430 format(" newton scheme for finding the nearest wall point",/,
c$$$     .  "  did not converge ",/,
c$$$     .  " i=",i4,"x/zplasma=",2e15.8,"x/zwall=",2e15.8,
c$$$     . " theta=",1e10.3," f=",1e15.8,"fp="1e15.8)
c$$$  124 continue
c$$$      ww1(i) = xt
c$$$      ww2(i) = zt
c$$$  125 continue
c$$$c
c$$$      do 140 i = 2, mthalf - 1
c$$$      iq1 = mth1 - i + 1
c$$$      xwal(i) = ww1(i)
c$$$      xwal(iq1) = ww1(i)
c$$$      zwal(i) = ww2(i)
c$$$      zwal(iq1) = - ww2(i)
c$$$  140 continue
c$$$c
  145 continue
c
      iplt = 1
c$$$  146 continue
c$$$      if(.not. insect) go to 2000
c$$$      write(iotty,1450)inside
c$$$      write(outpest,1450)inside
c$$$      write(outmod,1450)inside
c$$$      call errmes(outpest,'vacdat')
c$$$c
c$$$ 1460 format(i3," xinf,zinf=",2e10.3," xwal,zwal=",2e10.3)
c$$$ 1450 format(" there are at least ",i3," wall points in the plasma")

 2000 continue

      xwal(mth2) = xwal(2)
      zwal(mth2) = zwal(2)
c
c.... Equal arcs on wall. Use xpp,zpp,ww1,ww2,ww3 as storage.
c

      if ( ( .not. infwal) .and. (leqarcw .eq. 1) ) then
         call atpoint ( "Equal Wall Arcs","Ishape", ishape,"a", a,
     $        iotty, outmod )
         call eqarcw ( xwal,zwal, xpp,zpp, ww1,ww2,ww3, mth1 )
         do i = 1, mth1
            xwal(i) = xpp(i)
            zwal(i) = zpp(i)
         end do
      end if

      xwal(mth2) = xwal(2)
      zwal(mth2) = zwal(2)
c
c...   Shift wall vertically by awv.          ! MSC May 05, 2007
      zwal(1:mth2) = zwal(1:mth2) + awv
c
c......  reset aw, bw
c
      aw = awsave
      bw = bwsave
c
      call atpoint ( "END OF WWALL","mth", mth,"a", a,
     $     iotty, outmod )

      RETURN
      END
c
c..................................................................
      SUBROUTINE iter_sym_shell ( xsh_out, zsh_out, nshl, ioshel,
     $     nout1, nout2 )
c...................................................................
c...ITER shell smoothed and symmeterized from Y.Gribov's data,
c   obtained with Matlab scripts.
c.. Uses 513  points. Smoothed symmetric shell. 
c.. Interpolates and outputs nshl points on 0 =< theta < 2pi in
c     xsh_out, zsh_out.

      REAL, DIMENSION(*) :: xsh_out, zsh_out

      REAL, DIMENSION(:), ALLOCATABLE :: xshelld0, zshelld0, angwal
  
      ALLOCATE ( xshelld0(520), zshelld0(520) )
 
      xshelld0 = (/ 
     $  8.9733266945,   8.9734601541,   8.9738531521, 
     $  8.9744877692,   8.9753353916,   8.9763582680, 
     $  8.9775104401,   8.9787403810,   8.9799939070, 
     $  8.9812157952,   8.9823513473,   8.9833490243, 
     $  8.9841614422,   8.9847460980,   8.9850658019, 
     $  8.9850891138,   8.9847909784,   8.9841511996, 
     $  8.9831543718,   8.9817891375,   8.9800482602, 
     $  8.9779273467,   8.9754243148,   8.9725389917, 
     $  8.9692726815,   8.9656279461,   8.9616080129, 
     $  8.9572166631,   8.9524580304,   8.9473364600, 
     $  8.9418563194,   8.9360218999,   8.9298373292, 
     $  8.9233065246,   8.9164330028,   8.9092199048, 
     $  8.9016699185,   8.8937852546,   8.8855675611, 
     $  8.8770178552,   8.8681365878,   8.8589236053, 
     $  8.8493781955,   8.8394989477,   8.8292839145, 
     $  8.8187306213,   8.8078361073,   8.7965969526, 
     $  8.7850092906,   8.7730689409,   8.7607714288, 
     $  8.7481120430,   8.7350858463,   8.7216877546, 
     $  8.7079125673,   8.6937549896,   8.6792096548, 
     $  8.6642711303,   8.6489339397,   8.6331925602, 
     $  8.6170414246,   8.6004749116,   8.5834873388, 
     $  8.5660729594,   8.5482259529,   8.5299404258, 
     $  8.5112103756,   8.4920297120,   8.4723922460, 
     $  8.4522916931,   8.4317216653,   8.4106756724, 
     $  8.3891471439,   8.3671294332,   8.3446158322, 
     $  8.3215995750,   8.2980738695,   8.2740319106, 
     $  8.2494668946,   8.2243720349,   8.1987405850, 
     $  8.1725658429,   8.1458411600,   8.1185599421, 
     $  8.0907156626,   8.0623018521,   8.0333120918, 
     $  8.0037400102,   7.9735792771,   7.9428236090, 
     $  7.9114667611,   7.8795025398,   7.8469248094, 
     $  7.8137275226,   7.7799047557,   7.7454507337, 
     $  7.7103598698,   7.6746267752,   7.6382463868, 
     $  7.6012139193,   7.5635249026,   7.5251751733, 
     $  7.4861608949,   7.4464785205,   7.4061246861, 
     $  7.3650961431,   7.3233896788,   7.2810019447, 
     $  7.2379292980,   7.1941676623,   7.1497123991, 
     $  7.1045581550,   7.0586985686,   7.0121263521, 
     $  6.9648331981,   6.9168098363,   6.8680458989, 
     $  6.8185302613,   6.7682513110,   6.7171972458, 
     $  6.6653563772,   6.6127179455,   6.5592726008, 
     $  6.5050130511,   6.4499344909,   6.3940357625, 
     $  6.3373201152,   6.2797953954,   6.2214746875, 
     $  6.1623761493,   6.1025252164,   6.0419528844, 
     $  5.9806962056,   5.9187979254,   5.8563067285, 
     $  5.7932769214,   5.7297671936,   5.6658404004, 
     $  5.6015628160,   5.5370042379,   5.4722366168, 
     $  5.4073337365,   5.3423707800,   5.2774239774, 
     $  5.2125702520,   5.1478868482,   5.0834512228, 
     $  5.0193409197,   4.9556334914,   4.8924064876, 
     $  4.8297374992,   4.7677042191,   4.7063845033, 
     $  4.6458565977,   4.5861992220,   4.5274917440, 
     $  4.4698142815,   4.4132480164,   4.3578753791, 
     $  4.3037800586,   4.2510471040,   4.1997628335, 
     $  4.1500150767,   4.1018923877,   4.0554836074, 
     $  4.0108772058,   3.9681602757,   3.9274164890, 
     $  3.8887247528,   3.8521571931,   3.8177782046, 
     $  3.7856389925,   3.7557771451,   3.7282144196, 
     $  3.7029552766,   3.6799848657,   3.6592648312, 
     $  3.6407372103,   3.6243235027,   3.6099266969, 
     $  3.5974294596,   3.5866999021,   3.5775948139, 
     $  3.5699625682,   3.5636462720,   3.5584884230, 
     $  3.5543351138,   3.5510387878,   3.5484602968, 
     $  3.5464729846,   3.5449643489,   3.5438355153, 
     $  3.5430019657,   3.5423921929,   3.5419508531, 
     $  3.5416331091,   3.5414045111,   3.5412392082, 
     $  3.5411190021,   3.5410314382,   3.5409677338, 
     $  3.5409220364,   3.5408904175,   3.5408703174, 
     $  3.5408596849,   3.5408568126,   3.5408601610, 
     $  3.5408682645,   3.5408797381,   3.5408933044, 
     $  3.5409078356,   3.5409223728,   3.5409361966, 
     $  3.5409488030,   3.5409598798,   3.5409692870, 
     $  3.5409770264,   3.5409832162,   3.5409880260, 
     $  3.5409916594,   3.5409943303,   3.5409962421, 
     $  3.5409975722,   3.5409984720,   3.5409990636, 
     $  3.5409994435,   3.5409996790,   3.5409998203, 
     $  3.5409999025,   3.5409999490,   3.5409999745, 
     $  3.5409999877,   3.5409999943,   3.5409999974, 
     $  3.5409999989,   3.5409999996,   3.5409999998, 
     $  3.5409999999,   3.5410000000,   3.5410000000, 
     $  3.5410000000,   3.5410000000,   3.5410000000, 
     $  3.5410000000,   3.5410000000,   3.5410000000, 
     $  3.5410000000,   3.5410000000,   3.5410000000, 
     $  3.5410000000,   3.5410000000,   3.5410000000, 
     $  3.5410000000,   3.5410000000,   3.5410000000, 
     $  3.5410000000,   3.5410000000,   3.5410000000, 
     $  3.5410000000,   3.5410000000,   3.5410000000, 
     $  3.5410000000,   3.5410000000,   3.5410000000, 
     $  3.5410000000,   3.5410000000,   3.5410000000, 
     $  3.5410000000,   3.5410000000,   3.5410000000, 
     $  3.5410000000,   3.5410000000,   3.5410000000, 
     $  3.5410000000,   3.5410000000,   3.5409999999, 
     $  3.5409999998,   3.5409999996,   3.5409999989, 
     $  3.5409999974,   3.5409999943,   3.5409999877, 
     $  3.5409999745,   3.5409999490,   3.5409999025, 
     $  3.5409998203,   3.5409996790,   3.5409994435, 
     $  3.5409990636,   3.5409984720,   3.5409975722, 
     $  3.5409962421,   3.5409943303,   3.5409916594, 
     $  3.5409880260,   3.5409832162,   3.5409770264, 
     $  3.5409692870,   3.5409598798,   3.5409488030, 
     $  3.5409361966,   3.5409223728,   3.5409078356, 
     $  3.5408933044,   3.5408797381,   3.5408682645, 
     $  3.5408601610,   3.5408568126,   3.5408596849, 
     $  3.5408703174,   3.5408904175,   3.5409220364, 
     $  3.5409677338,   3.5410314382,   3.5411190021, 
     $  3.5412392082,   3.5414045111,   3.5416331091, 
     $  3.5419508531,   3.5423921929,   3.5430019657, 
     $  3.5438355153,   3.5449643489,   3.5464729846, 
     $  3.5484602968,   3.5510387878,   3.5543351138, 
     $  3.5584884230,   3.5636462720,   3.5699625682, 
     $  3.5775948139,   3.5866999021,   3.5974294596, 
     $  3.6099266969,   3.6243235027,   3.6407372103, 
     $  3.6592648312,   3.6799848657,   3.7029552766, 
     $  3.7282144196,   3.7557771451,   3.7856389925, 
     $  3.8177782046,   3.8521571931,   3.8887247528, 
     $  3.9274164890,   3.9681602757,   4.0108772058, 
     $  4.0554836074,   4.1018923877,   4.1500150767, 
     $  4.1997628335,   4.2510471040,   4.3037800586, 
     $  4.3578753791,   4.4132480164,   4.4698142815, 
     $  4.5274917440,   4.5861992220,   4.6458565977, 
     $  4.7063845033,   4.7677042191,   4.8297374992, 
     $  4.8924064876,   4.9556334914,   5.0193409197, 
     $  5.0834512228,   5.1478868482,   5.2125702520, 
     $  5.2774239774,   5.3423707800,   5.4073337365, 
     $  5.4722366168,   5.5370042379,   5.6015628160, 
     $  5.6658404004,   5.7297671936,   5.7932769214, 
     $  5.8563067285,   5.9187979254,   5.9806962056, 
     $  6.0419528844,   6.1025252164,   6.1623761493, 
     $  6.2214746875,   6.2797953954,   6.3373201152, 
     $  6.3940357625,   6.4499344909,   6.5050130511, 
     $  6.5592726008,   6.6127179455,   6.6653563772, 
     $  6.7171972458,   6.7682513110,   6.8185302613, 
     $  6.8680458989,   6.9168098363,   6.9648331981, 
     $  7.0121263521,   7.0586985686,   7.1045581550, 
     $  7.1497123991,   7.1941676623,   7.2379292980, 
     $  7.2810019447,   7.3233896788,   7.3650961431, 
     $  7.4061246861,   7.4464785205,   7.4861608949, 
     $  7.5251751733,   7.5635249026,   7.6012139193, 
     $  7.6382463868,   7.6746267752,   7.7103598698, 
     $  7.7454507337,   7.7799047557,   7.8137275226, 
     $  7.8469248094,   7.8795025398,   7.9114667611, 
     $  7.9428236090,   7.9735792771,   8.0037400102, 
     $  8.0333120918,   8.0623018521,   8.0907156626, 
     $  8.1185599421,   8.1458411600,   8.1725658429, 
     $  8.1987405850,   8.2243720349,   8.2494668946, 
     $  8.2740319106,   8.2980738695,   8.3215995750, 
     $  8.3446158322,   8.3671294332,   8.3891471439, 
     $  8.4106756724,   8.4317216653,   8.4522916931, 
     $  8.4723922460,   8.4920297120,   8.5112103756, 
     $  8.5299404258,   8.5482259529,   8.5660729594, 
     $  8.5834873388,   8.6004749116,   8.6170414246, 
     $  8.6331925602,   8.6489339397,   8.6642711303, 
     $  8.6792096548,   8.6937549896,   8.7079125673, 
     $  8.7216877546,   8.7350858463,   8.7481120430, 
     $  8.7607714288,   8.7730689409,   8.7850092906, 
     $  8.7965969526,   8.8078361073,   8.8187306213, 
     $  8.8292839145,   8.8394989477,   8.8493781955, 
     $  8.8589236053,   8.8681365878,   8.8770178552, 
     $  8.8855675611,   8.8937852546,   8.9016699185, 
     $  8.9092199048,   8.9164330028,   8.9233065246, 
     $  8.9298373292,   8.9360218999,   8.9418563194, 
     $  8.9473364600,   8.9524580304,   8.9572166631, 
     $  8.9616080129,   8.9656279461,   8.9692726815, 
     $  8.9725389917,   8.9754243148,   8.9779273467, 
     $  8.9800482602,   8.9817891375,   8.9831543718, 
     $  8.9841511996,   8.9847909784,   8.9850891138, 
     $  8.9850658019,   8.9847460980,   8.9841614422, 
     $  8.9833490243,   8.9823513473,   8.9812157952, 
     $  8.9799939070,   8.9787403810,   8.9775104401, 
     $  8.9763582680,   8.9753353916,   8.9744877692, 
     $  8.9738531521,   8.9734601541,   8.9733266945
     $      /)                             
      
      zshelld0 = (/
     $  0.0000000000,  -0.0339545925,  -0.0679205991, 
     $ -0.1019091241,  -0.1359307565,  -0.1699954015, 
     $ -0.2041120210,  -0.2382883975,  -0.2725312178, 
     $ -0.3068459894,  -0.3412371287,  -0.3757077251, 
     $ -0.4102599062,  -0.4448949017,  -0.4796131649, 
     $ -0.5144144733,  -0.5492980418,  -0.5842627580, 
     $ -0.6193072542,  -0.6544300116,  -0.6896294559, 
     $ -0.7249040702,  -0.7602524481,  -0.7956733439, 
     $ -0.8311657073,  -0.8667287819,  -0.9023620665, 
     $ -0.9380653461,  -0.9738386816,  -1.0096824516, 
     $ -1.0455973222,  -1.0815842039,  -1.1176442372, 
     $ -1.1537787515,  -1.1899892694,  -1.2262774081, 
     $ -1.2626448554,  -1.2990933293,  -1.3356245368, 
     $ -1.3722401146,  -1.4089416044,  -1.4457304213, 
     $ -1.4826078442,  -1.5195749489,  -1.5566326318, 
     $ -1.5937816002,  -1.6310223746,  -1.6683552848, 
     $ -1.7057804564,  -1.7432978582,  -1.7809073007, 
     $ -1.8186084566,  -1.8564008479,  -1.8942838777, 
     $ -1.9322568408,  -1.9703189272,  -2.0084692298, 
     $ -2.0467067236,  -2.0850302873,  -2.1234386892, 
     $ -2.1619305868,  -2.2005044981,  -2.2391587890, 
     $ -2.2778916696,  -2.3167011736,  -2.3555851633, 
     $ -2.3945412444,  -2.4335668139,  -2.4726590306, 
     $ -2.5118148161,  -2.5510308277,  -2.5903034446, 
     $ -2.6296288101,  -2.6690028254,  -2.7084211776, 
     $ -2.7478793026,  -2.7873724560,  -2.8268957204, 
     $ -2.8664440217,  -2.9060121381,  -2.9455946903, 
     $ -2.9851861643,  -3.0247808909,  -3.0643730433, 
     $ -3.1039565713,  -3.1435251834,  -3.1830723208, 
     $ -3.2225911249,  -3.2620744222,  -3.3015146221, 
     $ -3.3409038111,  -3.3802337532,  -3.4194959423, 
     $ -3.4586816433,  -3.4977820657,  -3.5367885304, 
     $ -3.5756926638,  -3.6144865152,  -3.6531631584, 
     $ -3.6917168382,  -3.7301433915,  -3.7684405277, 
     $ -3.8066084246,  -3.8446504681,  -3.8825732440, 
     $ -3.9203870132,  -3.9581056453,  -3.9957478642, 
     $ -4.0333362477,  -4.0708970309,  -4.1084595756, 
     $ -4.1460558153,  -4.1837188620,  -4.2214811052, 
     $ -4.2593725846,  -4.2974196375,  -4.3356412369, 
     $ -4.3740467400,  -4.4126340844,  -4.4513876874, 
     $ -4.4902772458,  -4.5292517122,  -4.5682426118, 
     $ -4.6071622885,  -4.6459048755,  -4.6843441682, 
     $ -4.7223364135,  -4.7597238415,  -4.7963364617, 
     $ -4.8319951237,  -4.8665138535,  -4.8997059739, 
     $ -4.9313867741,  -4.9613764944,  -4.9895042625, 
     $ -5.0156126736,  -5.0395582498,  -5.0612136316, 
     $ -5.0804673194,  -5.0972289331,  -5.1114249476, 
     $ -5.1229985119,  -5.1319082340,  -5.1381273053, 
     $ -5.1416427956,  -5.1424511600,  -5.1405573696, 
     $ -5.1359728875,  -5.1287148000,  -5.1188032087, 
     $ -5.1062601658,  -5.0911088118,  -5.0733726315, 
     $ -5.0530746152,  -5.0302369454,  -5.0048807592, 
     $ -4.9770260497,  -4.9466913812,  -4.9138939421, 
     $ -4.8786497425,  -4.8409737209,  -4.8008799659, 
     $ -4.7583820373,  -4.7134938413,  -4.6662304092, 
     $ -4.6166087605,  -4.5646493706,  -4.5103787435, 
     $ -4.4538306754,  -4.3950485693,  -4.3340862725, 
     $ -4.2710142931,  -4.2059192960,  -4.1389059775, 
     $ -4.0700979433,  -3.9996390697,  -3.9276961493, 
     $ -3.8544535954,  -3.7801129210,  -3.7048897143, 
     $ -3.6290122982,  -3.5527146544,  -3.4762321035, 
     $ -3.3997970782,  -3.3236353938,  -3.2479575640, 
     $ -3.1729576675,  -3.0988095762,  -3.0256658549, 
     $ -2.9536524267,  -2.8828688619,  -2.8133910944, 
     $ -2.7452715090,  -2.6785421324,  -2.6132114177, 
     $ -2.5492735566,  -2.4867092007,  -2.4254884823, 
     $ -2.3655729839,  -2.3069189247,  -2.2494798451, 
     $ -2.1932078347,  -2.1380547472,  -2.0839735842, 
     $ -2.0309189961,  -1.9788474990,  -1.9277175804, 
     $ -1.8774897675,  -1.8281267178,  -1.7795926584, 
     $ -1.7318533640,  -1.6848759246,  -1.6386288461, 
     $ -1.5930817577,  -1.5482053393,  -1.5039712982, 
     $ -1.4603523071,  -1.4173222039,  -1.3748557231, 
     $ -1.3329285588,  -1.2915172946,  -1.2505995040, 
     $ -1.2101536599,  -1.1701590244,  -1.1305956431, 
     $ -1.0914442421,  -1.0526863869,  -1.0143042168, 
     $ -0.9762804739,  -0.9385984443,  -0.9012419659, 
     $ -0.8641954038,  -0.8274435515,  -0.7909716313, 
     $ -0.7547652390,  -0.7188104104,  -0.6830935082, 
     $ -0.6476012162,  -0.6123205164,  -0.5772386771, 
     $ -0.5423432641,  -0.5076220676,  -0.4730631067, 
     $ -0.4386546027,  -0.4043849881,  -0.3702428955, 
     $ -0.3362170830,  -0.3022964586,  -0.2684702226, 
     $ -0.2347272215,  -0.2010573496,  -0.1674496894, 
     $ -0.1338920252,  -0.1003807205,  -0.0668908398, 
     $ -0.0334267623,   0.0000000000,   0.0334267623, 
     $  0.0668908398,   0.1003807205,   0.1338920252, 
     $  0.1674496894,   0.2010573496,   0.2347272215, 
     $  0.2684702226,   0.3022964586,   0.3362170830, 
     $  0.3702428955,   0.4043849881,   0.4386546027, 
     $  0.4730631067,   0.5076220676,   0.5423432641, 
     $  0.5772386771,   0.6123205164,   0.6476012162, 
     $  0.6830935082,   0.7188104104,   0.7547652390, 
     $  0.7909716313,   0.8274435515,   0.8641954038, 
     $  0.9012419659,   0.9385984443,   0.9762804739, 
     $  1.0143042168,   1.0526863869,   1.0914442421, 
     $  1.1305956431,   1.1701590244,   1.2101536599, 
     $  1.2505995040,   1.2915172946,   1.3329285588, 
     $  1.3748557231,   1.4173222039,   1.4603523071, 
     $  1.5039712982,   1.5482053393,   1.5930817577, 
     $  1.6386288461,   1.6848759246,   1.7318533640, 
     $  1.7795926584,   1.8281267178,   1.8774897675, 
     $  1.9277175804,   1.9788474990,   2.0309189961, 
     $  2.0839735842,   2.1380547472,   2.1932078347, 
     $  2.2494798451,   2.3069189247,   2.3655729839, 
     $  2.4254884823,   2.4867092007,   2.5492735566, 
     $  2.6132114177,   2.6785421324,   2.7452715090, 
     $  2.8133910944,   2.8828688619,   2.9536524267, 
     $  3.0256658549,   3.0988095762,   3.1729576675, 
     $  3.2479575640,   3.3236353938,   3.3997970782, 
     $  3.4762321035,   3.5527146544,   3.6290122982, 
     $  3.7048897143,   3.7801129210,   3.8544535954, 
     $  3.9276961493,   3.9996390697,   4.0700979433, 
     $  4.1389059775,   4.2059192960,   4.2710142931, 
     $  4.3340862725,   4.3950485693,   4.4538306754, 
     $  4.5103787435,   4.5646493706,   4.6166087605, 
     $  4.6662304092,   4.7134938413,   4.7583820373, 
     $  4.8008799659,   4.8409737209,   4.8786497425, 
     $  4.9138939421,   4.9466913812,   4.9770260497, 
     $  5.0048807592,   5.0302369454,   5.0530746152, 
     $  5.0733726315,   5.0911088118,   5.1062601658, 
     $  5.1188032087,   5.1287148000,   5.1359728875, 
     $  5.1405573696,   5.1424511600,   5.1416427956, 
     $  5.1381273053,   5.1319082340,   5.1229985119, 
     $  5.1114249476,   5.0972289331,   5.0804673194, 
     $  5.0612136316,   5.0395582498,   5.0156126736, 
     $  4.9895042625,   4.9613764944,   4.9313867741, 
     $  4.8997059739,   4.8665138535,   4.8319951237, 
     $  4.7963364617,   4.7597238415,   4.7223364135, 
     $  4.6843441682,   4.6459048755,   4.6071622885, 
     $  4.5682426118,   4.5292517122,   4.4902772458, 
     $  4.4513876874,   4.4126340844,   4.3740467400, 
     $  4.3356412369,   4.2974196375,   4.2593725846, 
     $  4.2214811052,   4.1837188620,   4.1460558153, 
     $  4.1084595756,   4.0708970309,   4.0333362477, 
     $  3.9957478642,   3.9581056453,   3.9203870132, 
     $  3.8825732440,   3.8446504681,   3.8066084246, 
     $  3.7684405277,   3.7301433915,   3.6917168382, 
     $  3.6531631584,   3.6144865152,   3.5756926638, 
     $  3.5367885304,   3.4977820657,   3.4586816433, 
     $  3.4194959423,   3.3802337532,   3.3409038111, 
     $  3.3015146221,   3.2620744222,   3.2225911249, 
     $  3.1830723208,   3.1435251834,   3.1039565713, 
     $  3.0643730433,   3.0247808909,   2.9851861643, 
     $  2.9455946903,   2.9060121381,   2.8664440217, 
     $  2.8268957204,   2.7873724560,   2.7478793026, 
     $  2.7084211776,   2.6690028254,   2.6296288101, 
     $  2.5903034446,   2.5510308277,   2.5118148161, 
     $  2.4726590306,   2.4335668139,   2.3945412444, 
     $  2.3555851633,   2.3167011736,   2.2778916696, 
     $  2.2391587890,   2.2005044981,   2.1619305868, 
     $  2.1234386892,   2.0850302873,   2.0467067236, 
     $  2.0084692298,   1.9703189272,   1.9322568408, 
     $  1.8942838777,   1.8564008479,   1.8186084566, 
     $  1.7809073007,   1.7432978582,   1.7057804564, 
     $  1.6683552848,   1.6310223746,   1.5937816002, 
     $  1.5566326318,   1.5195749489,   1.4826078442, 
     $  1.4457304213,   1.4089416044,   1.3722401146, 
     $  1.3356245368,   1.2990933293,   1.2626448554, 
     $  1.2262774081,   1.1899892694,   1.1537787515, 
     $  1.1176442372,   1.0815842039,   1.0455973222, 
     $  1.0096824516,   0.9738386816,   0.9380653461, 
     $  0.9023620665,   0.8667287819,   0.8311657073, 
     $  0.7956733439,   0.7602524481,   0.7249040702, 
     $  0.6896294559,   0.6544300116,   0.6193072542, 
     $  0.5842627580,   0.5492980418,   0.5144144733, 
     $  0.4796131649,   0.4448949017,   0.4102599062, 
     $  0.3757077251,   0.3412371287,   0.3068459894, 
     $  0.2725312178,   0.2382883975,   0.2041120210, 
     $  0.1699954015,   0.1359307565,   0.1019091241, 
     $  0.0679205991,   0.0339545925,   0.0000000000
     $       /)      

      npts = 512
      npts1 = npts + 1
      mth = nshl
      mth1 = mth + 1
      mth2 = mth1 + 1
      pi = acos(-1.0)

      WRITE ( ioshel,
     $     '(/,20x,"******  ITER SYMMETRIC SHELL *******",/)' )  
      WRITE ( ioshel, '("**** npts1, mth1 = ", 2i5,/)' ) npts1, mth1

c.. Now interpolate to MTH points:

      WRITE ( ioshel,
     $     '(//,"#",10x, "Interpolating to MTH points.",/,"#")' )

      CALL transdx ( xshelld0,npts, xsh_out,mth, 0.0, 0.0 )
      CALL transdx ( zshelld0,npts, zsh_out,mth, 0.0, 0.0 )

      xsh_out(mth2) = xsh_out(2)
      zsh_out(mth2) = zsh_out(2)

      xcensh = ( xsh_out(1) + xsh_out(mth/2+1) ) / 2.0
      zcensh = 0.0
      WRITE ( ioshel, '("#",10x, "Amtn", 7x, "Xwal-mth",5x,"Zwal-mth",
     $     5x,
     $     "Angwal",/,"#")' )

c Scale mth grid to the npts abscissae. Get the angle that the wall
c   points subtends form the center of the vessel.

      ALLOCATE ( angwal(mth+5) )

      DO i = 1, mth1
         amtnp = (i-1)*REAL(npts)/mth + 1
         angwal(i) =  - ATAN2 ( zsh_out(i)-zcensh,xsh_out(i)-xcensh )
         if ( angwal(i) < 0.0 ) angwal(i) = angwal(i) + 2.0*pi
         angwal(i) = 180.0 * angwal(i) / pi
         WRITE ( ioshel, '(2x, i4, 4es13.5 )' ) i, amtnp, xsh_out(i),
     $        zsh_out(i), angwal(i)
      END DO

      DEALLOCATE ( xshelld0, zshelld0, angwal )

c      CLOSE ( UNIT = ioshel )
c      CLOSE ( UNIT = 21 )

      END SUBROUTINE iter_sym_shell


c..................................................................
      SUBROUTINE nstx_shell ( xsh_out, zsh_out, nshl, lsymwv, ioshel,
     $     nout1, nout2 )
c...................................................................

c.. Uses the 1026  points provided by Jon Menard. Thers is an option to 
c   symmetrize the shell using the upper points (lsymwv = 1) or lower
c   points (lsymwv = 2). Using 513 points in both cases.
c   Smoothes the data if lsmth_nstx is == 1. This is currently set below
c      in this subroutine. Can be input data if necessary.
c.. Interpolates and outputs nshl points on 0 =< theta < 2pi in
c     xsh_out, zsh_out.
c.. nshl = no of output points on the open domain.
c.   Points are clockwise.


      REAL, DIMENSION(*) :: xsh_out, zsh_out

      REAL, DIMENSION(:), ALLOCATABLE :: xwal, zwal, angwal, thgrd,
     $           xwalq, zwalq, angwalq, ellwq, thgrq

      REAL, DIMENSION(:), ALLOCATABLE :: anglp, angx, angz, angxz

      REAL, DIMENSION(:), ALLOCATABLE :: xshell0, zshell0,
     $     xshell1, zshell1, xshell2, zshell2, xshell3, zshell3
 
      REAL, DIMENSION(1026) :: xshelld0 = (/ 
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6957191598e+00,
     $     1.6929338401e+00, 1.6884469961e+00, 1.6823776986e+00,
     $     1.6753957019e+00, 1.6678306860e+00, 1.6599620571e+00,
     $     1.6519457820e+00, 1.6438539382e+00, 1.6357214862e+00,
     $     1.6275625467e+00, 1.6194044961e+00, 1.6112431777e+00,
     $     1.6030735634e+00, 1.5949149108e+00, 1.5867512589e+00,
     $     1.5785978174e+00, 1.5704303798e+00, 1.5622669997e+00,
     $     1.5541094735e+00, 1.5459684486e+00, 1.5378598381e+00,
     $     1.5298203745e+00, 1.5219524103e+00, 1.5144476372e+00,
     $     1.5075461273e+00, 1.5017790945e+00, 1.4974955271e+00,
     $     1.4946714820e+00, 1.4920011114e+00, 1.4894410660e+00,
     $     1.4868861685e+00, 1.4843294522e+00, 1.4817718729e+00,
     $     1.4792164840e+00, 1.4766589899e+00, 1.4741027022e+00,
     $     1.4715461164e+00, 1.4689899215e+00, 1.4664321545e+00,
     $     1.4638773323e+00, 1.4613192607e+00, 1.4587630304e+00,
     $     1.4562072045e+00, 1.4536492758e+00, 1.4510942956e+00,
     $     1.4485361649e+00, 1.4459807317e+00, 1.4434232560e+00,
     $     1.4408677907e+00, 1.4383098016e+00, 1.4357542593e+00,
     $     1.4331977754e+00, 1.4306400642e+00, 1.4280845795e+00,
     $     1.4255276560e+00, 1.4229712449e+00, 1.4204138371e+00,
     $     1.4178596879e+00, 1.4153282011e+00, 1.4128311997e+00,
     $     1.4103534293e+00, 1.4078866859e+00, 1.4054213030e+00,
     $     1.4029511631e+00, 1.4004669583e+00, 1.3979615364e+00,
     $     1.3954287221e+00, 1.3928598807e+00, 1.3902569438e+00,
     $     1.3876319107e+00, 1.3849899057e+00, 1.3823423292e+00,
     $     1.3796918956e+00, 1.3770454676e+00, 1.3743871014e+00,
     $     1.3717176921e+00, 1.3690500145e+00, 1.3663979603e+00,
     $     1.3637339462e+00, 1.3604641515e+00, 1.3571491286e+00,
     $     1.3536937708e+00, 1.3501214467e+00, 1.3464317289e+00,
     $     1.3426237695e+00, 1.3386915365e+00, 1.3346186990e+00,
     $     1.3303971366e+00, 1.3260179101e+00, 1.3214776169e+00,
     $     1.3167649366e+00, 1.3118853502e+00, 1.3068312272e+00,
     $     1.3016248588e+00, 1.2969063344e+00, 1.2921678715e+00,
     $     1.2874692485e+00, 1.2827549443e+00, 1.2780437300e+00,
     $     1.2733318253e+00, 1.2686202134e+00, 1.2639084566e+00,
     $     1.2591961092e+00, 1.2544848748e+00, 1.2497748079e+00,
     $     1.2450592282e+00, 1.2403513529e+00, 1.2356377945e+00,
     $     1.2309260726e+00, 1.2262146538e+00, 1.2215035491e+00,
     $     1.2167908421e+00, 1.2120792704e+00, 1.2073680780e+00,
     $     1.2026562775e+00, 1.1979444802e+00, 1.1932320365e+00,
     $     1.1885216619e+00, 1.1838085862e+00, 1.1790987223e+00,
     $     1.1743874805e+00, 1.1696794803e+00, 1.1650373358e+00,
     $     1.1605620489e+00, 1.1562177900e+00, 1.1519789187e+00,
     $     1.1478751437e+00, 1.1438654088e+00, 1.1399317760e+00,
     $     1.1361475101e+00, 1.1324603924e+00, 1.1288665387e+00,
     $     1.1254068660e+00, 1.1219287475e+00, 1.1183317682e+00,
     $     1.1146388179e+00, 1.1108688895e+00, 1.1069501669e+00,
     $     1.1029649802e+00, 1.0988516601e+00, 1.0946183307e+00,
     $     1.0903034927e+00, 1.0858267145e+00, 1.0812483953e+00,
     $     1.0765874428e+00, 1.0716803945e+00, 1.0665795200e+00,
     $     1.0612487703e+00, 1.0555631411e+00, 1.0495697871e+00,
     $     1.0432932465e+00, 1.0367134929e+00, 1.0298003246e+00,
     $     1.0225112579e+00, 1.0149669640e+00, 1.0073698477e+00,
     $     9.9977373912e-01, 9.9217717902e-01, 9.8458394292e-01,
     $     9.7698487044e-01, 9.6939099849e-01, 9.6179471878e-01,
     $     9.5419833030e-01, 9.4660178828e-01, 9.3900710872e-01,
     $     9.3141009855e-01, 9.2381501340e-01, 9.1621805287e-01,
     $     9.0862230551e-01, 9.0102714762e-01, 8.9342937497e-01,
     $     8.8583608744e-01, 8.7823777600e-01, 8.7064249627e-01,
     $     8.6304680514e-01, 8.5545060349e-01, 8.4785477879e-01,
     $     8.4025849935e-01, 8.3266329374e-01, 8.2506801888e-01,
     $     8.1747023993e-01, 8.0987707221e-01, 8.0227969425e-01,
     $     7.9468419446e-01, 7.8708820134e-01, 7.7949335377e-01,
     $     7.7189624086e-01, 7.6430183397e-01, 7.5670737072e-01,
     $     7.4910274250e-01, 7.4148096002e-01, 7.3383080829e-01,
     $     7.2615746038e-01, 7.1845860783e-01, 7.1073569976e-01,
     $     7.0297634977e-01, 6.9519296487e-01, 6.8738553162e-01,
     $     6.7954986600e-01, 6.7168955022e-01, 6.6380721740e-01,
     $     6.5590737568e-01, 6.4799444946e-01, 6.4005389230e-01,
     $     6.3210347752e-01, 6.2412808032e-01, 6.1612423195e-01,
     $     6.0807863536e-01, 6.0001701497e-01, 5.9193566433e-01,
     $     5.8384145240e-01, 5.7574116212e-01, 5.6762716643e-01,
     $     5.5950715829e-01, 5.5138083360e-01, 5.4324560235e-01,
     $     5.3510587969e-01, 5.2695938429e-01, 5.1880879528e-01,
     $     5.1065473172e-01, 5.0249648014e-01, 4.9433444035e-01,
     $     4.8617288352e-01, 4.7801149792e-01, 4.6984204671e-01,
     $     4.6168281772e-01, 4.5351791657e-01, 4.4535427006e-01,
     $     4.3718994937e-01, 4.2902709085e-01, 4.2086277885e-01,
     $     4.1269915178e-01, 4.0453598119e-01, 3.9637282147e-01,
     $     3.8820466126e-01, 3.8003875705e-01, 3.7190666167e-01,
     $     3.6381138237e-01, 3.5607714140e-01, 3.4840228124e-01,
     $     3.4073623447e-01, 3.3312675600e-01, 3.2559823706e-01,
     $     3.1863037365e-01, 3.1214555079e-01, 3.0631624791e-01,
     $     3.0114097070e-01, 2.9631584550e-01, 2.9216104786e-01,
     $     2.8829216069e-01, 2.8485563911e-01, 2.8171504009e-01,
     $     2.7889805829e-01, 2.7630110166e-01, 2.7401515771e-01,
     $     2.7181958159e-01, 2.6999464757e-01, 2.6819033661e-01,
     $     2.6666754532e-01, 2.6525138936e-01, 2.6390392038e-01,
     $     2.6283591707e-01, 2.6177341417e-01, 2.6083645834e-01,
     $     2.6007241579e-01, 2.5930434827e-01, 2.5853380765e-01,
     $     2.5775967049e-01, 2.5699576052e-01, 2.5633329103e-01,
     $     2.5567703255e-01, 2.5507368990e-01, 2.5468955699e-01,
     $     2.5430887134e-01, 2.5400000811e-01, 2.5400000811e-01,
     $     2.5400000811e-01, 2.5400000811e-01, 2.5400000811e-01,
     $     2.5400000522e-01, 2.5399919287e-01, 2.5399462206e-01,
     $     2.5398988167e-01, 2.5393469132e-01, 2.5377790590e-01,
     $     2.5361926661e-01, 2.5345740714e-01, 2.5329381140e-01,
     $     2.5312693613e-01, 2.5286838584e-01, 2.5260798722e-01,
     $     2.5226494879e-01, 2.5185479992e-01, 2.5140576393e-01,
     $     2.5082927439e-01, 2.5024356876e-01, 2.4948285310e-01,
     $     2.4871855073e-01, 2.4775531301e-01, 2.4677361475e-01,
     $     2.4557459398e-01, 2.4431450647e-01, 2.4285927545e-01,
     $     2.4125386822e-01, 2.3951389587e-01, 2.3747882466e-01,
     $     2.3527032922e-01, 2.3285557922e-01, 2.3007640741e-01,
     $     2.2700366494e-01, 2.2364450350e-01, 2.2015329642e-01,
     $     2.1679569138e-01, 2.1342700571e-01, 2.1025701761e-01,
     $     2.0736374654e-01, 2.0473369610e-01, 2.0234270468e-01,
     $     2.0017501223e-01, 1.9826434143e-01, 1.9658811229e-01,
     $     1.9509272713e-01, 1.9376029523e-01, 1.9263267660e-01,
     $     1.9168300838e-01, 1.9086633545e-01, 1.9019405252e-01,
     $     1.8971303170e-01, 1.8933882300e-01, 1.8908084444e-01,
     $     1.8899382198e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8895000219e-01, 1.8895000219e-01,
     $     1.8895000219e-01, 1.8891701714e-01, 1.8888427811e-01,
     $     1.8885176971e-01, 1.8881965641e-01, 1.8878800848e-01,
     $     1.8875690186e-01, 1.8874134650e-01, 1.8873010388e-01,
     $     1.8872052709e-01, 1.8871107386e-01, 1.8869813789e-01,
     $     1.8868743185e-01, 1.8867828350e-01, 1.8866930843e-01,
     $     1.8866704922e-01, 1.8866652040e-01, 1.8866651556e-01,
     $     1.8866651555e-01, 1.8866651555e-01, 1.8866651555e-01,
     $     1.8869053424e-01, 1.8874520811e-01, 1.8884246245e-01,
     $     1.8910357413e-01, 1.8946399532e-01, 1.8995144202e-01,
     $     1.9063004584e-01, 1.9142597737e-01, 1.9234402289e-01,
     $     1.9348186079e-01, 1.9478057072e-01, 1.9623839788e-01,
     $     1.9787609394e-01, 1.9979120754e-01, 2.0191475987e-01,
     $     2.0425253600e-01, 2.0681741153e-01, 2.0963998023e-01,
     $     2.1279822505e-01, 2.1613877162e-01, 2.1949832134e-01,
     $     2.2295746085e-01, 2.2639322695e-01, 2.2954461030e-01,
     $     2.3232515370e-01, 2.3481365342e-01, 2.3708903788e-01,
     $     2.3914525847e-01, 2.4092673200e-01, 2.4258706904e-01,
     $     2.4404299836e-01, 2.4535238867e-01, 2.4655211170e-01,
     $     2.4757831296e-01, 2.4854701043e-01, 2.4934642880e-01,
     $     2.5010785794e-01, 2.5072906411e-01, 2.5130623776e-01,
     $     2.5178737399e-01, 2.5219819862e-01, 2.5257020957e-01,
     $     2.5283064474e-01, 2.5308945588e-01, 2.5327491846e-01,
     $     2.5343849131e-01, 2.5360142736e-01, 2.5375995150e-01,
     $     2.5391742017e-01, 2.5399987151e-01, 2.5400000811e-01,
     $     2.5400000811e-01, 2.5400000811e-01, 2.5400000811e-01,
     $     2.5400000811e-01, 2.5400000811e-01, 2.5400000811e-01,
     $     2.5400000811e-01, 2.5423683052e-01, 2.5461451192e-01,
     $     2.5499850683e-01, 2.5555155028e-01, 2.5620433217e-01,
     $     2.5686568068e-01, 2.5760982784e-01, 2.5838274498e-01,
     $     2.5915404984e-01, 2.5992246230e-01, 2.6068727999e-01,
     $     2.6156457234e-01, 2.6262805689e-01, 2.6369473197e-01,
     $     2.6497944786e-01, 2.6638965754e-01, 2.6783853514e-01,
     $     2.6963760989e-01, 2.7146249919e-01, 2.7357332714e-01,
     $     2.7585307528e-01, 2.7835155130e-01, 2.8116115846e-01,
     $     2.8418800014e-01, 2.8761603591e-01, 2.9135440794e-01,
     $     2.9550057044e-01, 3.0017555122e-01, 3.0517932319e-01,
     $     3.1099674354e-01, 3.1731755928e-01, 3.2412662765e-01,
     $     3.3164642020e-01, 3.3923665829e-01, 3.4689721694e-01,
     $     3.5455510115e-01, 3.6227369666e-01, 3.7031540781e-01,
     $     3.7844758099e-01, 3.8660551835e-01, 3.9477668266e-01,
     $     4.0293920918e-01, 4.1110276015e-01, 4.1926604450e-01,
     $     4.2743048063e-01, 4.3559350498e-01, 4.4375749157e-01,
     $     4.5192134520e-01, 4.6008884196e-01, 4.6823875868e-01,
     $     4.7641693140e-01, 4.8458441016e-01, 4.9273108814e-01,
     $     5.0090559340e-01, 5.0906447677e-01, 5.1723613654e-01,
     $     5.2539647459e-01, 5.3355457263e-01, 5.4171909004e-01,
     $     5.4988899838e-01, 5.5804668048e-01, 5.6619418821e-01,
     $     5.7433095096e-01, 5.8245483532e-01, 5.9055939503e-01,
     $     5.9863762375e-01, 6.0669509133e-01, 6.1472571186e-01,
     $     6.2274149218e-01, 6.3068489505e-01, 6.3863136015e-01,
     $     6.4656098437e-01, 6.5443790094e-01, 6.6228871236e-01,
     $     6.7011043016e-01, 6.7788495268e-01, 6.8564646868e-01,
     $     6.9333318531e-01, 7.0097015521e-01, 7.0854413460e-01,
     $     7.1612555124e-01, 7.2372226872e-01, 7.3132306313e-01,
     $     7.3891470339e-01, 7.4651166794e-01, 7.5410733420e-01,
     $     7.6170354699e-01, 7.6929863271e-01, 7.7689497231e-01,
     $     7.8449048336e-01, 7.9208610553e-01, 7.9968135069e-01,
     $     8.0727892015e-01, 8.1487268237e-01, 8.2246858429e-01,
     $     8.3006764792e-01, 8.3765804804e-01, 8.4525763246e-01,
     $     8.5285244085e-01, 8.6044837910e-01, 8.6804487512e-01,
     $     8.7563908491e-01, 8.8323777419e-01, 8.9083190985e-01,
     $     8.9842816384e-01, 9.0602499495e-01, 9.1361919080e-01,
     $     9.2121709353e-01, 9.2881202776e-01, 9.3640850700e-01,
     $     9.4400440319e-01, 9.5159925514e-01, 9.5919781250e-01,
     $     9.6679085534e-01, 9.7438863793e-01, 9.8198456387e-01,
     $     9.8958003064e-01, 9.9717506554e-01, 1.0047725990e+00,
     $     1.0123661619e+00, 1.0199630177e+00, 1.0275621653e+00,
     $     1.0349445250e+00, 1.0420226728e+00, 1.0487810068e+00,
     $     1.0552206106e+00, 1.0612456073e+00, 1.0669501556e+00,
     $     1.0723737043e+00, 1.0774781229e+00, 1.0822097436e+00,
     $     1.0868049898e+00, 1.0912782304e+00, 1.0956104030e+00,
     $     1.0998613041e+00, 1.1039749604e+00, 1.1079756397e+00,
     $     1.1119106543e+00, 1.1156808279e+00, 1.1193865332e+00,
     $     1.1229900781e+00, 1.1264831938e+00, 1.1299128565e+00,
     $     1.1332024939e+00, 1.1365264388e+00, 1.1399856762e+00,
     $     1.1435461116e+00, 1.1471721089e+00, 1.1509506251e+00,
     $     1.1548091896e+00, 1.1587586599e+00, 1.1629207197e+00,
     $     1.1672384736e+00, 1.1717332203e+00, 1.1764228946e+00,
     $     1.1811414006e+00, 1.1858521359e+00, 1.1905638859e+00,
     $     1.1952744686e+00, 1.1999873878e+00, 1.2046988331e+00,
     $     1.2094102915e+00, 1.2141219019e+00, 1.2188334505e+00,
     $     1.2235468341e+00, 1.2282552750e+00, 1.2329700303e+00,
     $     1.2376382795e+00, 1.2422615936e+00, 1.2468439732e+00,
     $     1.2513894002e+00, 1.2559009182e+00, 1.2603862360e+00,
     $     1.2648559787e+00, 1.2693096469e+00, 1.2737685547e+00,
     $     1.2782294400e+00, 1.2827152241e+00, 1.2872202738e+00,
     $     1.2917374691e+00, 1.2962588771e+00, 1.3007758826e+00,
     $     1.3052770191e+00, 1.3097823738e+00, 1.3143137358e+00,
     $     1.3188549317e+00, 1.3233927813e+00, 1.3279137255e+00,
     $     1.3323886647e+00, 1.3367532251e+00, 1.3409912478e+00,
     $     1.3450885865e+00, 1.3490396854e+00, 1.3528312450e+00,
     $     1.3564612043e+00, 1.3599374450e+00, 1.3632725901e+00,
     $     1.3664737176e+00, 1.3695523911e+00, 1.3725191447e+00,
     $     1.3753838894e+00, 1.3781579873e+00, 1.3808482142e+00,
     $     1.3834662928e+00, 1.3860258765e+00, 1.3885829720e+00,
     $     1.3911387217e+00, 1.3936945788e+00, 1.3962536731e+00,
     $     1.3988078454e+00, 1.4013654150e+00, 1.4039220592e+00,
     $     1.4064779214e+00, 1.4090358069e+00, 1.4115910172e+00,
     $     1.4141490036e+00, 1.4167032843e+00, 1.4192632757e+00,
     $     1.4218161772e+00, 1.4243757060e+00, 1.4269307152e+00,
     $     1.4294877763e+00, 1.4320439161e+00, 1.4346015152e+00,
     $     1.4371568522e+00, 1.4397141355e+00, 1.4422710996e+00,
     $     1.4448264704e+00, 1.4473842375e+00, 1.4499402720e+00,
     $     1.4524965816e+00, 1.4550540953e+00, 1.4576097440e+00,
     $     1.4601665482e+00, 1.4627234610e+00, 1.4652799964e+00,
     $     1.4678358406e+00, 1.4703934682e+00, 1.4729496995e+00,
     $     1.4755058147e+00, 1.4780625889e+00, 1.4806198773e+00,
     $     1.4831755356e+00, 1.4857325920e+00, 1.4882902942e+00,
     $     1.4908422079e+00, 1.4934058202e+00, 1.4960822625e+00,
     $     1.4992746094e+00, 1.5039866401e+00, 1.5099885370e+00,
     $     1.5169851753e+00, 1.5245385223e+00, 1.5323955163e+00,
     $     1.5404131492e+00, 1.5485050947e+00, 1.5566368890e+00,
     $     1.5647893528e+00, 1.5729510897e+00, 1.5811204714e+00,
     $     1.5892832902e+00, 1.5974121464e+00, 1.6056319198e+00,
     $     1.6137666750e+00, 1.6219288066e+00, 1.6300933364e+00,
     $     1.6382339675e+00, 1.6463390079e+00, 1.6543825114e+00,
     $     1.6623168804e+00, 1.6700473528e+00, 1.6773929516e+00,
     $     1.6840885962e+00, 1.6896816032e+00, 1.6936677033e+00,
     $     1.6960800757e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00,
     $     1.6969500780e+00, 1.6969500780e+00, 1.6969500780e+00
     $      /)                             
      
      REAL, DIMENSION(1026) :: zshelld0 = (/
     $      0.0000000000e+00, -8.1637245383e-03, -1.6314626970e-02,
     $     -2.4506010571e-02, -3.2622390605e-02, -4.0840317837e-02,
     $     -4.8968589922e-02, -5.7143555003e-02, -6.5307127055e-02,
     $     -7.3468464527e-02, -8.1635827225e-02, -8.9796123401e-02,
     $     -9.7960690907e-02, -1.0612683622e-01, -1.1428531150e-01,
     $     -1.2245634992e-01, -1.3061381826e-01, -1.3877904443e-01,
     $     -1.4694488922e-01, -1.5510765199e-01, -1.6326908878e-01,
     $     -1.7143656533e-01, -1.7959703175e-01, -1.8776319234e-01,
     $     -1.9592523880e-01, -2.0408950648e-01, -2.1225326309e-01,
     $     -2.2041793004e-01, -2.2857980748e-01, -2.3674505563e-01,
     $     -2.4490707231e-01, -2.5307363034e-01, -2.6123451677e-01,
     $     -2.6939966194e-01, -2.7756330123e-01, -2.8572577205e-01,
     $     -2.9389255368e-01, -3.0205269726e-01, -3.1021839328e-01,
     $     -3.1838166510e-01, -3.2654538785e-01, -3.3470950994e-01,
     $     -3.4287334990e-01, -3.5103476815e-01, -3.5920255249e-01,
     $     -3.6736250460e-01, -3.7552881735e-01, -3.8369119473e-01,
     $     -3.9185429673e-01, -4.0002067182e-01, -4.0818116567e-01,
     $     -4.1634767006e-01, -4.2450933991e-01, -4.3267427150e-01,
     $     -4.4083690411e-01, -4.4900212180e-01, -4.5716519597e-01,
     $     -4.6532716124e-01, -4.7349468927e-01, -4.8165463221e-01,
     $     -4.8982220599e-01, -4.9798473290e-01, -5.0613524220e-01,
     $     -5.1431990897e-01, -5.2248044743e-01, -5.3058509055e-01,
     $     -5.3830934835e-01, -5.4520384824e-01, -5.5074652466e-01,
     $     -5.5503116814e-01, -5.5815249354e-01, -5.6036259448e-01,
     $     -5.6193126280e-01, -5.6303016346e-01, -5.6376198112e-01,
     $     -5.6416585406e-01, -5.6422254012e-01, -5.6413073011e-01,
     $     -5.6391696620e-01, -5.6370356199e-01, -5.6348493943e-01,
     $     -5.6326157111e-01, -5.6319771664e-01, -5.6331472416e-01,
     $     -5.6365426789e-01, -5.6427204077e-01, -5.6524722850e-01,
     $     -5.6669500609e-01, -5.6893192964e-01, -5.7218542848e-01,
     $     -5.7663090836e-01, -5.8247583632e-01, -5.8949597931e-01,
     $     -5.9715315472e-01, -6.0488500371e-01, -6.1263428154e-01,
     $     -6.2038225619e-01, -6.2813576563e-01, -6.3589185755e-01,
     $     -6.4364131278e-01, -6.5139715930e-01, -6.5914934186e-01,
     $     -6.6690243544e-01, -6.7465433104e-01, -6.8241098459e-01,
     $     -6.9015871402e-01, -6.9791630290e-01, -7.0566831703e-01,
     $     -7.1341910217e-01, -7.2117624883e-01, -7.2892447712e-01,
     $     -7.3668224027e-01, -7.4443183409e-01, -7.5218762644e-01,
     $     -7.5993731393e-01, -7.6769466396e-01, -7.7544459310e-01,
     $     -7.8319739594e-01, -7.9095389885e-01, -7.9870365901e-01,
     $     -8.0645777779e-01, -8.1421033717e-01, -8.2196593274e-01,
     $     -8.2971598362e-01, -8.3747941673e-01, -8.4525156005e-01,
     $     -8.5303041868e-01, -8.6081142699e-01, -8.6859585103e-01,
     $     -8.7637547157e-01, -8.8415237737e-01, -8.9192299491e-01,
     $     -8.9968187031e-01, -9.0743319003e-01, -9.1516991675e-01,
     $     -9.2289849775e-01, -9.3062731531e-01, -9.3834473272e-01,
     $     -9.4606849951e-01, -9.5379180356e-01, -9.6151030995e-01,
     $     -9.6922557542e-01, -9.7694307010e-01, -9.8465146578e-01,
     $     -9.9238159842e-01, -9.9988587714e-01, -1.0073094949e+00,
     $     -1.0147248298e+00, -1.0220659135e+00, -1.0293450053e+00,
     $     -1.0365693287e+00, -1.0437219652e+00, -1.0507982406e+00,
     $     -1.0577861425e+00, -1.0646780795e+00, -1.0714597973e+00,
     $     -1.0781314006e+00, -1.0846724027e+00, -1.0910843442e+00,
     $     -1.0973996837e+00, -1.1040105833e+00, -1.1107152246e+00,
     $     -1.1173635040e+00, -1.1240339658e+00, -1.1307000596e+00,
     $     -1.1373671330e+00, -1.1440337891e+00, -1.1507006467e+00,
     $     -1.1573683430e+00, -1.1640344688e+00, -1.1706989350e+00,
     $     -1.1773711950e+00, -1.1840325694e+00, -1.1907019793e+00,
     $     -1.1973687915e+00, -1.2040351721e+00, -1.2107011044e+00,
     $     -1.2173693060e+00, -1.2240359074e+00, -1.2307019648e+00,
     $     -1.2373688829e+00, -1.2440357977e+00, -1.2507036280e+00,
     $     -1.2573685313e+00, -1.2640372552e+00, -1.2707014388e+00,
     $     -1.2773685443e+00, -1.2840547645e+00, -1.2907422502e+00,
     $     -1.2975882217e+00, -1.3044969894e+00, -1.3114753222e+00,
     $     -1.3185309696e+00, -1.3256431001e+00, -1.3327971598e+00,
     $     -1.3400317714e+00, -1.3473206373e+00, -1.3546313992e+00,
     $     -1.3620430130e+00, -1.3694351607e+00, -1.3767494661e+00,
     $     -1.3840358009e+00, -1.3912756503e+00, -1.3984391104e+00,
     $     -1.4055648128e+00, -1.4126152180e+00, -1.4195961955e+00,
     $     -1.4265278929e+00, -1.4333552456e+00, -1.4401192362e+00,
     $     -1.4468012259e+00, -1.4533541886e+00, -1.4597200460e+00,
     $     -1.4659004648e+00, -1.4717676762e+00, -1.4773158604e+00,
     $     -1.4825414190e+00, -1.4873753849e+00, -1.4917351233e+00,
     $     -1.4954365249e+00, -1.4985510210e+00, -1.5015433435e+00,
     $     -1.5045352723e+00, -1.5075273741e+00, -1.5105180597e+00,
     $     -1.5135104814e+00, -1.5165008585e+00, -1.5194921890e+00,
     $     -1.5224835599e+00, -1.5254749894e+00, -1.5284656887e+00,
     $     -1.5314573115e+00, -1.5344481765e+00, -1.5374397839e+00,
     $     -1.5404309131e+00, -1.5434218104e+00, -1.5464137410e+00,
     $     -1.5494039056e+00, -1.5523960563e+00, -1.5553870140e+00,
     $     -1.5583781355e+00, -1.5613694606e+00, -1.5643606586e+00,
     $     -1.5673525964e+00, -1.5703441649e+00, -1.5733357596e+00,
     $     -1.5763283410e+00, -1.5793191092e+00, -1.5823115354e+00,
     $     -1.5853032210e+00, -1.5882951013e+00, -1.5912865318e+00,
     $     -1.5942788496e+00, -1.5972701045e+00, -1.6002610566e+00,
     $     -1.6032448428e+00, -1.6061596809e+00, -1.6090095689e+00,
     $     -1.6117968405e+00, -1.6145178703e+00, -1.6171573603e+00,
     $     -1.6197000982e+00, -1.6221619674e+00, -1.6245439331e+00,
     $     -1.6268452434e+00, -1.6290407355e+00, -1.6311697952e+00,
     $     -1.6332330325e+00, -1.6352297195e+00, -1.6371401061e+00,
     $     -1.6389881925e+00, -1.6407317748e+00, -1.6423328189e+00,
     $     -1.6437512017e+00, -1.6450073773e+00, -1.6461669232e+00,
     $     -1.6472510033e+00, -1.6482510590e+00, -1.6491611618e+00,
     $     -1.6500023832e+00, -1.6507751562e+00, -1.6514787664e+00,
     $     -1.6520874479e+00, -1.6526253762e+00, -1.6530908130e+00,
     $     -1.6534817699e+00, -1.6537878398e+00, -1.6539999865e+00,
     $     -1.6541368578e+00, -1.6542000771e+00, -1.6542000771e+00,
     $     -1.6542000771e+00, -1.6542000771e+00, -1.6542000771e+00,
     $     -1.6542000771e+00, -1.6542000771e+00, -1.6542000771e+00,
     $     -1.6542000771e+00, -1.6542000771e+00, -1.6542000771e+00,
     $     -1.6542000578e+00, -1.6541863230e+00, -1.6537561914e+00,
     $     -1.6523136026e+00, -1.6496944576e+00, -1.6469191153e+00,
     $     -1.6440921009e+00, -1.6412516407e+00, -1.6379244587e+00,
     $     -1.6336908346e+00, -1.6287146613e+00, -1.6229661417e+00,
     $     -1.6166498280e+00, -1.6100578354e+00, -1.6030183997e+00,
     $     -1.5958324532e+00, -1.5884204770e+00, -1.5808867578e+00,
     $     -1.5732198192e+00, -1.5654814490e+00, -1.5576417409e+00,
     $     -1.5497784190e+00, -1.5418196644e+00, -1.5338573364e+00,
     $     -1.5258346361e+00, -1.5177909340e+00, -1.5097600992e+00,
     $     -1.5016456941e+00, -1.4935513166e+00, -1.4854364494e+00,
     $     -1.4773293741e+00, -1.4691806412e+00, -1.4610629731e+00,
     $     -1.4529290663e+00, -1.4448250773e+00, -1.4366605499e+00,
     $     -1.4285312253e+00, -1.4203970530e+00, -1.4122325551e+00,
     $     -1.4040803922e+00, -1.3959257809e+00, -1.3877571406e+00,
     $     -1.3795977622e+00, -1.3714218441e+00, -1.3632716698e+00,
     $     -1.3551028404e+00, -1.3469412256e+00, -1.3387769577e+00,
     $     -1.3306120711e+00, -1.3224427577e+00, -1.3142859771e+00,
     $     -1.3061296563e+00, -1.2979627329e+00, -1.2897988789e+00,
     $     -1.2816347503e+00, -1.2734802521e+00, -1.2653175150e+00,
     $     -1.2571614518e+00, -1.2490088825e+00, -1.2408568746e+00,
     $     -1.2327132756e+00, -1.2245707470e+00, -1.2164421870e+00,
     $     -1.2083132688e+00, -1.2002068345e+00, -1.1921026580e+00,
     $     -1.1840259507e+00, -1.1759604985e+00, -1.1679262912e+00,
     $     -1.1599224533e+00, -1.1519450305e+00, -1.1440378002e+00,
     $     -1.1361784864e+00, -1.1283807655e+00, -1.1206956274e+00,
     $     -1.1131364532e+00, -1.1056947407e+00, -1.0983087730e+00,
     $     -1.0908853073e+00, -1.0834154587e+00, -1.0759288951e+00,
     $     -1.0682691021e+00, -1.0605446246e+00, -1.0527472212e+00,
     $     -1.0448557981e+00, -1.0369412343e+00, -1.0289363145e+00,
     $     -1.0209116416e+00, -1.0128605523e+00, -1.0047719310e+00,
     $     -9.9666415173e-01, -9.8854163828e-01, -9.8040352511e-01,
     $     -9.7225504731e-01, -9.6410161021e-01, -9.5593304310e-01,
     $     -9.4777837812e-01, -9.3962557707e-01, -9.3143035477e-01,
     $     -9.2328919393e-01, -9.1512501232e-01, -9.0695346682e-01,
     $     -8.9879463295e-01, -8.9063388666e-01, -8.8246065484e-01,
     $     -8.7430495332e-01, -8.6614017864e-01, -8.5796981340e-01,
     $     -8.4981610099e-01, -8.4164637771e-01, -8.3348299401e-01,
     $     -8.2532042327e-01, -8.1715645540e-01, -8.0899275075e-01,
     $     -8.0082858514e-01, -7.9266396735e-01, -7.8450188694e-01,
     $     -7.7634537872e-01, -7.6815749693e-01, -7.6002349284e-01,
     $     -7.5184681702e-01, -7.4367373208e-01, -7.3552861446e-01,
     $     -7.2734685540e-01, -7.1919846870e-01, -7.1102430302e-01,
     $     -7.0286455194e-01, -6.9470484155e-01, -6.8652820850e-01,
     $     -6.7837721759e-01, -6.7021350901e-01, -6.6204285976e-01,
     $     -6.5387599158e-01, -6.4572737551e-01, -6.3754712472e-01,
     $     -6.2939368281e-01, -6.2123293614e-01, -6.1304589111e-01,
     $     -6.0491865479e-01, -5.9672820152e-01, -5.8857153101e-01,
     $     -5.8040011698e-01, -5.7226195562e-01, -5.6407251825e-01,
     $     -5.5591013735e-01, -5.4776375634e-01, -5.3958348176e-01,
     $     -5.3143334383e-01, -5.2325050676e-01, -5.1510934872e-01,
     $     -5.0692941283e-01, -4.9877157472e-01, -4.9060641818e-01,
     $     -4.8244970297e-01, -4.7427596676e-01, -4.6611046514e-01,
     $     -4.5795985998e-01, -4.4978492001e-01, -4.4163348563e-01,
     $     -4.3344698849e-01, -4.2530599169e-01, -4.1713700456e-01,
     $     -4.0896031156e-01, -4.0082147049e-01, -3.9261860441e-01,
     $     -3.8450637217e-01, -3.7629574100e-01, -3.6816344312e-01,
     $     -3.5997387852e-01, -3.5183953559e-01, -3.4365141141e-01,
     $     -3.3550046654e-01, -3.2732537315e-01, -3.1918293850e-01,
     $     -3.1098963175e-01, -3.0285252577e-01, -2.9467513469e-01,
     $     -2.8651163615e-01, -2.7835371502e-01, -2.7018483855e-01,
     $     -2.6202744999e-01, -2.5385457156e-01, -2.4569135214e-01,
     $     -2.3754610660e-01, -2.2935416671e-01, -2.2121086216e-01,
     $     -2.1303745448e-01, -2.0488481645e-01, -1.9670293152e-01,
     $     -1.8855217444e-01, -1.8039419101e-01, -1.7221249465e-01,
     $     -1.6406123143e-01, -1.5589671602e-01, -1.4772482001e-01,
     $     -1.3957973023e-01, -1.3139024449e-01, -1.2324050049e-01,
     $     -1.1509249633e-01, -1.0688793963e-01, -9.8764447449e-02,
     $     -9.0588814368e-02, -8.2402540016e-02, -7.4276606073e-02,
     $     -6.6078644124e-02, -5.7941118202e-02, -4.9764344351e-02,
     $     -4.1594708020e-02, -3.3444989763e-02, -2.5280553853e-02,
     $     -1.7091591288e-02, -8.9743469226e-03, -7.6561473687e-04,
     $      7.3793603255e-03,  1.5526813413e-02,  2.3734770343e-02,
     $      3.1849329119e-02,  4.0040019886e-02,  4.8200277332e-02,
     $      5.6360591640e-02,  6.4528966206e-02,  7.2692131658e-02,
     $      8.0837514137e-02,  8.9038005449e-02,  9.7169358985e-02,
     $      1.0534265049e-01,  1.1351733361e-01,  1.2166566668e-01,
     $      1.2983109317e-01,  1.3801352399e-01,  1.4614476166e-01,
     $      1.5435112905e-01,  1.6246233652e-01,  1.7067154649e-01,
     $      1.7882431563e-01,  1.8695794553e-01,  1.9517421876e-01,
     $      2.0327732702e-01,  2.1151405528e-01,  2.1958864284e-01,
     $      2.2784488519e-01,  2.3593058638e-01,  2.4414975619e-01,
     $      2.5227686815e-01,  2.6046675041e-01,  2.6861119633e-01,
     $      2.7678286167e-01,  2.8496508010e-01,  2.9306448207e-01,
     $      3.0133474246e-01,  3.0938899084e-01,  3.1762179407e-01,
     $      3.2575971425e-01,  3.3393639935e-01,  3.4209063368e-01,
     $      3.5024129103e-01,  3.5843844259e-01,  3.6658536387e-01,
     $      3.7472210216e-01,  3.8293786813e-01,  3.9106575222e-01,
     $      3.9924052386e-01,  4.0739065326e-01,  4.1557771460e-01,
     $      4.2373624451e-01,  4.3187646254e-01,  4.4006970352e-01,
     $      4.4820686704e-01,  4.5640186926e-01,  4.6455047594e-01,
     $      4.7268538160e-01,  4.8090274197e-01,  4.8902430443e-01,
     $      4.9721801434e-01,  5.0534675574e-01,  5.1355834837e-01,
     $      5.2165639052e-01,  5.2989652601e-01,  5.3799631271e-01,
     $      5.4620726500e-01,  5.5432990283e-01,  5.6250531549e-01,
     $      5.7070907917e-01,  5.7882029077e-01,  5.8700862131e-01,
     $      5.9515765299e-01,  6.0334757123e-01,  6.1148767158e-01,
     $      6.1965473588e-01,  6.2783586551e-01,  6.3597239245e-01,
     $      6.4415753917e-01,  6.5232471951e-01,  6.6044077373e-01,
     $      6.6869052325e-01,  6.7677400281e-01,  6.8497537640e-01,
     $      6.9313360893e-01,  7.0130106661e-01,  7.0945002266e-01,
     $      7.1763305870e-01,  7.2578392895e-01,  7.3395313468e-01,
     $      7.4211581081e-01,  7.5025837638e-01,  7.5848899170e-01,
     $      7.6655915764e-01,  7.7479284756e-01,  7.8292958278e-01,
     $      7.9109653282e-01,  7.9926343180e-01,  8.0741839055e-01,
     $      8.1559510182e-01,  8.2374886275e-01,  8.3191514783e-01,
     $      8.4007826870e-01,  8.4824600987e-01,  8.5640606014e-01,
     $      8.6457096001e-01,  8.7273278428e-01,  8.8089654505e-01,
     $      8.8906275062e-01,  8.9722926197e-01,  9.0538452487e-01,
     $      9.1355424686e-01,  9.2172503795e-01,  9.2986141102e-01,
     $      9.3805272553e-01,  9.4621521703e-01,  9.5436304428e-01,
     $      9.6253204101e-01,  9.7068879795e-01,  9.7883428050e-01,
     $      9.8697267568e-01,  9.9509759981e-01,  1.0032085474e+00,
     $      1.0112958816e+00,  1.0193546380e+00,  1.0273794506e+00,
     $      1.0353997602e+00,  1.0433131779e+00,  1.0512083700e+00,
     $      1.0590335056e+00,  1.0667657209e+00,  1.0744657089e+00,
     $      1.0819526535e+00,  1.0894293174e+00,  1.0968610237e+00,
     $      1.1042545379e+00,  1.1116667158e+00,  1.1191894141e+00,
     $      1.1268767771e+00,  1.1346498053e+00,  1.1424904617e+00,
     $      1.1503925500e+00,  1.1583602751e+00,  1.1663531698e+00,
     $      1.1743874207e+00,  1.1824447933e+00,  1.1905212980e+00,
     $      1.1986199410e+00,  1.2067258187e+00,  1.2148508977e+00,
     $      1.2229795690e+00,  1.2311193612e+00,  1.2392628852e+00,
     $      1.2474127274e+00,  1.2555663011e+00,  1.2637196039e+00,
     $      1.2718831119e+00,  1.2800385466e+00,  1.2882016001e+00,
     $      1.2963617619e+00,  1.3045358023e+00,  1.3126860769e+00,
     $      1.3208484904e+00,  1.3290108089e+00,  1.3371804840e+00,
     $      1.3453427642e+00,  1.3535056565e+00,  1.3616725360e+00,
     $      1.3698273311e+00,  1.3779966567e+00,  1.3861607766e+00,
     $      1.3943266751e+00,  1.4024853672e+00,  1.4106374024e+00,
     $      1.4187991290e+00,  1.4269409490e+00,  1.4350696211e+00,
     $      1.4432280557e+00,  1.4513425693e+00,  1.4594757866e+00,
     $      1.4675861816e+00,  1.4757382527e+00,  1.4838535528e+00,
     $      1.4919603934e+00,  1.5000664543e+00,  1.5081707579e+00,
     $      1.5162229737e+00,  1.5242601449e+00,  1.5322990389e+00,
     $      1.5402626269e+00,  1.5482212071e+00,  1.5561069309e+00,
     $      1.5639483381e+00,  1.5717194454e+00,  1.5793880966e+00,
     $      1.5869695982e+00,  1.5943843301e+00,  1.6016394348e+00,
     $      1.6086843844e+00,  1.6153741307e+00,  1.6218401962e+00,
     $      1.6275961256e+00,  1.6327696386e+00,  1.6372705177e+00,
     $      1.6406076807e+00,  1.6435357271e+00,  1.6463639980e+00,
     $      1.6491610344e+00,  1.6518653950e+00,  1.6534730014e+00,
     $      1.6541770536e+00,  1.6542000560e+00,  1.6542000771e+00,
     $      1.6542000771e+00,  1.6542000771e+00,  1.6542000771e+00,
     $      1.6542000771e+00,  1.6542000771e+00,  1.6542000771e+00,
     $      1.6542000771e+00,  1.6542000771e+00,  1.6542000771e+00,
     $      1.6542000771e+00,  1.6542000771e+00,  1.6542000771e+00,
     $      1.6542000771e+00,  1.6542000771e+00,  1.6542000771e+00,
     $      1.6542000771e+00,  1.6542000771e+00,  1.6542000771e+00,
     $      1.6542000771e+00,  1.6538750710e+00,  1.6533765230e+00,
     $      1.6527220265e+00,  1.6519119325e+00,  1.6509074902e+00,
     $      1.6497408565e+00,  1.6484267055e+00,  1.6469682387e+00,
     $      1.6453469936e+00,  1.6435520315e+00,  1.6416554798e+00,
     $      1.6396681245e+00,  1.6375618288e+00,  1.6353097380e+00,
     $      1.6329581220e+00,  1.6305179146e+00,  1.6279270485e+00,
     $      1.6251966864e+00,  1.6223208974e+00,  1.6192487640e+00,
     $      1.6162461051e+00,  1.6132534338e+00,  1.6102596617e+00,
     $      1.6072694934e+00,  1.6042772309e+00,  1.6012854784e+00,
     $      1.5982935093e+00,  1.5953019845e+00,  1.5923099709e+00,
     $      1.5893182809e+00,  1.5863265459e+00,  1.5833349609e+00,
     $      1.5803424595e+00,  1.5773514571e+00,  1.5743596138e+00,
     $      1.5713665258e+00,  1.5683768467e+00,  1.5653835518e+00,
     $      1.5623925986e+00,  1.5594013716e+00,  1.5564099315e+00,
     $      1.5534193953e+00,  1.5504270948e+00,  1.5474365942e+00,
     $      1.5444452616e+00,  1.5414537036e+00,  1.5384631864e+00,
     $      1.5354712068e+00,  1.5324804007e+00,  1.5294889870e+00,
     $      1.5264978062e+00,  1.5235070411e+00,  1.5205148186e+00,
     $      1.5175247611e+00,  1.5145328459e+00,  1.5115416633e+00,
     $      1.5085503190e+00,  1.5055588289e+00,  1.5025663387e+00,
     $      1.4995754127e+00,  1.4965831988e+00,  1.4935900734e+00,
     $      1.4900971607e+00,  1.4860162958e+00,  1.4814257839e+00,
     $      1.4764043003e+00,  1.4708802717e+00,  1.4650412505e+00,
     $      1.4589339903e+00,  1.4525685682e+00,  1.4459034012e+00,
     $      1.4391532577e+00,  1.4323320160e+00,  1.4254089596e+00,
     $      1.4184378991e+00,  1.4113877419e+00,  1.4042712711e+00,
     $      1.3971161562e+00,  1.3898765063e+00,  1.3826023690e+00,
     $      1.3752757206e+00,  1.3678907761e+00,  1.3605031963e+00,
     $      1.3530126041e+00,  1.3455529192e+00,  1.3381691452e+00,
     $      1.3308192109e+00,  1.3235047983e+00,  1.3162672107e+00,
     $      1.3090730987e+00,  1.3019271313e+00,  1.2949051640e+00,
     $      1.2879764843e+00,  1.2811473366e+00,  1.2744875687e+00,
     $      1.2678111095e+00,  1.2611456910e+00,  1.2544788370e+00,
     $      1.2478136392e+00,  1.2411451364e+00,  1.2344787121e+00,
     $      1.2278122692e+00,  1.2211456122e+00,  1.2144790455e+00,
     $      1.2078098865e+00,  1.2011477146e+00,  1.1944766175e+00,
     $      1.1877810300e+00,  1.1810540184e+00,  1.1742972149e+00,
     $      1.1675142374e+00,  1.1607116677e+00,  1.1538922019e+00,
     $      1.1470561526e+00,  1.1402205111e+00,  1.1333745854e+00,
     $      1.1265446288e+00,  1.1197200303e+00,  1.1129115349e+00,
     $      1.1061127913e+00,  1.0993157679e+00,  1.0925135831e+00,
     $      1.0857039121e+00,  1.0788974129e+00,  1.0721048359e+00,
     $      1.0653215550e+00,  1.0585354403e+00,  1.0517347332e+00,
     $      1.0449116520e+00,  1.0380104001e+00,  1.0310308470e+00,
     $      1.0239690833e+00,  1.0168277622e+00,  1.0095957035e+00,
     $      1.0022807519e+00,  9.9489730513e-01,  9.8744292648e-01,
     $      9.7993379392e-01,  9.7237281806e-01,  9.6476623442e-01,
     $      9.5712227070e-01,  9.4944213582e-01,  9.4173752882e-01,
     $      9.3400370218e-01,  9.2625093075e-01,  9.1849633320e-01,
     $      9.1074577273e-01,  9.0299488500e-01,  8.9523418298e-01,
     $      8.8748840602e-01,  8.7973232541e-01,  8.7197905038e-01,
     $      8.6422814236e-01,  8.5647109710e-01,  8.4872217390e-01,
     $      8.4096483249e-01,  8.3321872800e-01,  8.2545529731e-01,
     $      8.1771337583e-01,  8.0995134766e-01,  8.0220303976e-01,
     $      7.9444851814e-01,  7.8669677639e-01,  7.7894061463e-01,
     $      7.7119131228e-01,  7.6343609749e-01,  7.5568184907e-01,
     $      7.4793243363e-01,  7.4017575584e-01,  7.3242431635e-01,
     $      7.2467206389e-01,  7.1691616324e-01,  7.0916590907e-01,
     $      7.0141213836e-01,  6.9365805455e-01,  6.8590511183e-01,
     $      6.7815426190e-01,  6.7039799687e-01,  6.6264598374e-01,
     $      6.5489432894e-01,  6.4714066692e-01,  6.3938545153e-01,
     $      6.3163517221e-01,  6.2388064976e-01,  6.1612416539e-01,
     $      6.0838524054e-01,  6.0061804943e-01,  5.9290504541e-01,
     $      5.8537655950e-01,  5.7864883420e-01,  5.7303903644e-01,
     $      5.6878207174e-01,  5.6562482639e-01,  5.6338642544e-01,
     $      5.6181149367e-01,  5.6071177389e-01,  5.5997906432e-01,
     $      5.5953636870e-01,  5.5932839026e-01,  5.5939528540e-01,
     $      5.5948506653e-01,  5.5958293131e-01,  5.5967765969e-01,
     $      5.5977514531e-01,  5.5976893225e-01,  5.5950286076e-01,
     $      5.5885734149e-01,  5.5786624144e-01,  5.5645001292e-01,
     $      5.5450211960e-01,  5.5184711072e-01,  5.4823514058e-01,
     $      5.4350531343e-01,  5.3748571500e-01,  5.3030003362e-01,
     $      5.2245022219e-01,  5.1431952627e-01,  5.0615066354e-01,
     $      4.9797077543e-01,  4.8982567638e-01,  4.8165485473e-01,
     $      4.7349463019e-01,  4.6532716415e-01,  4.5716518882e-01,
     $      4.4900212531e-01,  4.4083689947e-01,  4.3267427305e-01,
     $      4.2450933507e-01,  4.1634755733e-01,  4.0818174213e-01,
     $      4.0001947272e-01,  3.9185553937e-01,  3.8369054685e-01,
     $      3.7552894938e-01,  3.6736250359e-01,  3.5920255067e-01,
     $      3.5103476675e-01,  3.4287334821e-01,  3.3470950834e-01,
     $      3.2654538224e-01,  3.1838168484e-01,  3.1021834796e-01,
     $      3.0205274029e-01,  2.9389252945e-01,  2.8572577505e-01,
     $      2.7756329966e-01,  2.6939966037e-01,  2.6123451519e-01,
     $      2.5307362877e-01,  2.4490707073e-01,  2.3674505405e-01,
     $      2.2857980591e-01,  2.2041792847e-01,  2.1225326152e-01,
     $      2.0408950491e-01,  1.9592523723e-01,  1.8776319077e-01,
     $      1.7959703018e-01,  1.7143656376e-01,  1.6326908720e-01,
     $      1.5510765041e-01,  1.4694488764e-01,  1.3877904285e-01,
     $      1.3061381668e-01,  1.2245634835e-01,  1.1428530992e-01,
     $      1.0612683465e-01,  9.7960689333e-02,  8.9796121826e-02,
     $      8.1635825650e-02,  7.3468462952e-02,  6.5307183628e-02,
     $      5.7142984514e-02,  4.8970616326e-02,  4.0836737071e-02,
     $      3.2625770224e-02,  2.4504361469e-02,  1.6314950993e-02,
     $      8.1637245383e-03,  0.0000000000e+00, -8.1637245383e-03
     $       /)      

	
c$$$      xcensh = xshelld0(513)
c$$$      zcensh = zshelld0(513)
      xcensh = ( MAXVAL(xshelld0(1:1026)) +
     $     MINVAL(xshelld0(1:1026)) ) / 2.0
      xcensh = ( MAXVAL(zshelld0(1:1026)) +
     $     MINVAL(zshelld0(1:1026)) ) / 2.0

      pi = acos (-1.0 )

      nptsd = 1024
      nptsd1 = nptsd + 1
      mth = nshl
      mth1 = mth + 1
      mth2 = mth1 + 1

c..  Use [1:513] = 513 points on the upper half, Z >= 0.
c    and [513:1025] = 513 points on the lower, Z <= 0.

      npu = 513
      npl = 513

c.. The following are for vertical symmetrizing using the upper data.
      nptusy = 2*(npu-1)
      nptusy1 = nptusy + 1
      nptlsy = 2*(npl-1)
      nptlsy1 = nptlsy + 1
      npdim = max(nptusy,nptlsy,nptsd)

      ALLOCATE ( xshell0(npdim+5), zshell0(npdim+5),
     $     xshell1(npdim+5), zshell1(npdim+5),
     $     xshell2(npdim+5), zshell2(npdim+5),
     $     xshell3(npdim+5), zshell3(npdim+5) )

      xshell1 = 0.0
      xshell2 = 0.0
      xshell3 = 0.0
      zshell1 = 0.0
      zshell2 = 0.0
      zshell3 = 0.0

c..  Don't overwrite data arrays. So put X,Z elsewhere.

      xshell0(1:nptsd) = xshelld0(1:nptsd)
      zshell0(1:nptsd) = zshelld0(1:nptsd)
      xshell0(nptsd1) = xshelld0(1)
      zshell0(nptsd1) = zshelld0(1)

c... Now vertically symmetrize if LSYMWV is 1(upper), 2(lower):

      npts = nptsd
      npts1 = npts + 1
      
      IF ( lsymwv == 1 ) THEN

c.. Symmetry using upper points

         zshell0(1:nptusy) =
     $        (/ (-zshell0(i),i=nptsd1,npu,-1), zshell0(npu+1:nptsd) /)
         xshell0(1:nptusy) =
     $        (/  (xshell0(i),i=nptsd1,npu,-1), xshell0(npu+1:nptsd) /)
         xshell0(nptusy1) = xshell0(1)
         zshell0(nptusy1) = zshell0(1)
         
         npts = nptusy
         npts1 = npts + 1

      END IF

      IF ( lsymwv == 2 ) THEN

c.. Symmetry using lower points

         zshell0(1:nptlsy) =
     $        (/ zshell0(1:npl), (-zshell0(i),i=npl-1,2,-1) /)
         xshell0(1:nptlsy) =
     $        (/ xshell0(1:npl),  (xshell0(i),i=npl-1,2,-1) /)
         xshell0(nptlsy1) = xshell0(1)
         zshell0(nptlsy1) = zshell0(1)
         
         npts = nptlsy
         npts1 = npts + 1

      END IF

c.. For various degrees of smoothing, remove the appropriate comments,
c    and replace the appropriate inputs in the CALLs to transdx.

      lsmth_nstx = 0

      IF ( lsmth_nstx == 1 ) THEN
         CALL smooth1 ( xshell0, xshell1, npts1 )
         CALL smooth1 ( zshell0, zshell1, npts1 )
         
         CALL smooth1 ( xshell1, xshell2, npts1 )
         CALL smooth1 ( zshell1, zshell2, npts1 )
         
         CALL smooth1 ( xshell2, xshell3, npts1 )
         CALL smooth1 ( zshell2, zshell3, npts1 )
      END IF 

      WRITE ( ioshel, '(/,20x,"******  NSTX SHELL *******",/)' ) 
      WRITE ( ioshel, '("**** lsymwv = ", i5)' ) lsymwv
      WRITE ( ioshel, '("**** lsmth_nstx = ", i5)' ) lsmth_nstx
      WRITE ( ioshel, '("**** npts1, mth1 = ", 2i5,/)' ) npts1, mth1

      WRITE ( ioshel, '("#",9x, "Xshell_0",5x,"Zshell_0" ,5x,
     $     "Xshell_1",5x,"Zshell_1",5x,"Xshell_2",5x,"Zshell_2", 5x,
     $     "Xshell_3",5x,"Zshell_3")' )

      DO i = 1, npts1
         WRITE ( ioshel, '(2x, i4, 8es13.5 )' )
     $        i, xshell0(i), zshell0(i), xshell1(i), zshell1(i),
     $        xshell2(i), zshell2(i), xshell3(i), zshell3(i)
      END DO

c.. Now interpolate thrice smoothed data to MTH points:

      ALLOCATE ( xwal(mth+5), zwal(mth+5), angwal(mth+5),
     $     thgrd(mth+5) )
      ALLOCATE ( xwalq(mth+5), zwalq(mth+5), angwalq(mth+5),
     $     ellwq(mth+5),  thgrq(mth+5) )

      IF ( lsmth_nstx == 0 ) THEN
         xshell3(1:npts1) = xshell0(1:npts1)
         zshell3(1:npts1) = zshell0(1:npts1)
      END IF
      CALL transdx ( xshell3,npts, xwal,mth, 0.0, 0.0 )
      CALL transdx ( zshell3,npts, zwal,mth, 0.0, 0.0 )

      WRITE ( ioshel,
     $     '(//,"#",10x, "Interpolating to MTH points.",/,"#")' )
      WRITE ( ioshel, '("#",10x, "Amtn", 7x, "Xwal-mth",5x,"Zwal-mth",
     $     5x,
     $     "Angwal",/,"#")' )

c Scale mth grid to the npts abscissae. Get the angle that the wall
c   points subtends form the center of the vessel.

      DO i = 1, mth1
         amtnp = (i-1)*REAL(npts)/mth + 1
         angwal(i) =  - ATAN2 ( zwal(i) - zcensh,xwal(i)-xcensh )
         if ( angwal(i) < 0.0 ) angwal(i) = angwal(i) + 2.0*pi
         angwal(i) = 180.0 * angwal(i) / pi
         WRITE ( ioshel, '(2x, i4, 4es13.5 )' ) i,amtnp, xwal(i),
     $        zwal(i),
     $        angwal(i)
      END DO
c.. Ensure monoticity of the angles
      DO i = 1, mth/8
         j = mth1-i+1
         IF ( angwal(i) > 90.0 ) angwal(i) = angwal(i) - 360.0
         IF ( angwal(j) < 90.0 ) angwal(j) = 360.0 - angwal(j)
      END DO

c...Output thrice smoothed data expanded by the expansion factor, 
c     shfx = 1.05

      shfx = 1.00

      xsh_out(1:mth1) = xcensh + shfx*( xwal(1:mth1)-xcensh)
      zsh_out(1:mth1) = zcensh + shfx*( zwal(1:mth1)-zcensh)
      xsh_out(mth2) = xsh_out(2)
      zsh_out(mth2) = zsh_out(2)

      DEALLOCATE ( xshell0, zshell0, xshell1, zshell1, xshell2,
     $     zshell2, xshell3, zshell3, xwal, zwal, angwal )

      DEALLOCATE ( xwalq, zwalq, angwalq, ellwq, thgrq )

c      CLOSE ( UNIT = ioshel )
c      CLOSE ( UNIT = 21 )

      END SUBROUTINE nstx_shell

c..................................................................
      SUBROUTINE d3d_shell ( xsh_out, zsh_out, nshl, lsymwv, ioshel,
     $     nout1, nout2 )
c...................................................................
 
      REAL, DIMENSION(*) :: xsh_out, zsh_out

      REAL, DIMENSION(1001) :: thdum

      REAL, DIMENSION(:), ALLOCATABLE :: xshelld, zshelld
      REAL, DIMENSION(:), ALLOCATABLE :: xwal, zwal, angwal, thgrd,
     $           xwalq, zwalq, angwalq, ellwq, thgrq

      REAL, DIMENSION(:), ALLOCATABLE :: anglp, angx, angz, angxz

      REAL, DIMENSION(:), ALLOCATABLE :: xshell0, zshell0,
     $     xshell1, zshell1, xshell2, zshell2, xshell3, zshell3
 
      REAL, DIMENSION(105) :: xshelld0 = (/ 
     $     1.01600,
     $     1.016, 1.016, 1.016,                ! Extra points
     $     1.01600, 1.01800, 1.01800, 1.01800,
     $     1.01800, 1.01800, 1.01800, 1.01800, 1.01800,
     $     1.01600, 1.01200, 1.00100, 1.02900, 1.04200,
     $     1.04600, 1.05600, 1.09700, 1.10800, 1.11600,
     $     1.13400, 1.14800, 1.16200, 1.18100, 1.18200,
     $     1.18500, 1.19000, 1.19500, 1.20100, 1.20900,
     $     1.21500, 1.22200, 1.22800, 1.23400, 1.23900,
     $     1.24200, 1.24800, 1.25800, 1.26300, 1.28000,
     $     1.28000, 1.28000, 1.31000, 1.32800, 1.36100,
     $     1.38000, 1.41900, 1.41900, 1.37200, 1.37200,
     $     1.431, 1.490, 1.549,                 ! Extra points
     $     1.60800, 1.64700, 1.78500,
     $     1.9275,                              ! Extra point
     $     2.07000, 2.12800,
     $     2.1865,                              ! Extra point
     $     2.24500,
     $     2.284,                               ! Extra point
     $     2.32300, 2.37700, 2.35200, 2.35100,
     $     2.35100, 2.35400, 2.35400, 2.35400, 2.35300,
     $     2.35100, 2.35300, 2.37700,
     $     2.316, 2.2555, 2.19475,              ! Extra points
     $     2.13400,
     $     1.96,                                ! Extra point
     $     1.78600,
     $     1.76800, 1.76800, 1.68200, 1.67400, 1.67400,
     $     1.68100, 1.75400, 1.75400, 1.68900, 1.68900,
     $     1.55500, 1.41200, 1.27300, 1.15300,
     $     1.0845,                              ! Extra point
     $     1.01600,
     $     1.016,                               ! Extra point
     $     1.01600,
     $     1.016,                               ! Extra point
     $     1.01600, 1.01600,
     $     1.016 /)                             ! Extra point
      
      REAL, DIMENSION(105) :: zshelld0 = (/
     $     0.00000,
     $     0.242, 0.484, 0.726,                ! Extra points
     $     0.96400,  0.96800,  1.00100,  1.01900,
     $     1.07700,  1.07000,  1.09600,  1.11300,  1.13800, 
     $     1.14700,  1.16500,  1.21700,  1.21700,  1.16240, 
     $     1.16238,  1.16260,  1.16450,  1.16594,  1.16591, 
     $     1.16896,  1.17175,  1.17556,  1.18300,  1.18350, 
     $     1.18500,  1.18800,  1.19100,  1.19600,  1.20200, 
     $     1.20800,  1.21400,  1.22100,  1.23100,  1.23800, 
     $     1.24400,  1.25400,  1.27800,  1.29000,  1.33100, 
     $     1.34700,  1.34800,  1.34800,  1.34800,  1.34800, 
     $     1.34800,  1.34800,  1.31000,  1.31000,  1.29200,
     $     1.24275, 1.1935, 1.14425,            ! Extra points
     $     1.09500,  1.07700,  1.07700,
     $     1.0585,                              ! Exra point
     $     1.04000,  0.993000,
     $     0.851,                               ! Extra point
     $     0.709000,
     $     0.614,                               ! Extra point
     $     0.519000, 0.389000, 0.40000,  0.337000, 
     $     0.205000, 0.068000, 0.00000, -0.06800, -0.205000, 
     $     -0.33700, -0.40000, -0.38900,
     $     -0.535, -0.681, -0.827,              ! Extra points
     $     -0.97300,
     $     -1.0735,                             ! Extra point
     $     -1.17400,
     $     -1.21100, -1.25400, -1.25400, -1.26200, -1.32900, 
     $     -1.33900, -1.33900, -1.36800, -1.36800, -1.36600, 
     $     -1.36600, -1.36600, -1.36600, -1.36600,
     $     -1.2975,                             ! Extra point
     $     -1.22900,
     $     -1.0145,                             ! Extra point
     $     -0.80000,
     $     -0.6075,                             ! Extra point
     $     -0.41500, -0.400000,
     $     -0.20 /)                             ! Extra point

c$$$      OPEN ( UNIT = 21, FILE = "read_lao", STATUS = "OLD",
c$$$     $     FORM = "FORMATTED" )

c      lsymwv = 2

      xcensh = 1.6955
      zcensh = 0.0

      pi = acos (-1.0 )

      nptsd = 105
      nptsd1 = nptsd + 1
      mth = nshl
      mth1 = mth + 1
      mth2 = mth1 + 1

c..  There are [1:71] = 71 points on the upper half, Z >= 0.
c    and [71:105] = 35 points on the lower, Z <= 0.

      npu = 71
      npl = 35

c.. The following are for vertical symmetrizing using the upper data.
      nptusy = 2*(npu-1)
      nptusy1 = nptusy + 1
      nptlsy = 2*(npl-1)
      nptlsy1 = nptlsy + 1
      npdim = max(nptusy,nptlsy,nptsd)

      ALLOCATE ( xshelld(nptsd+5), zshelld(nptsd+5),
     $     xshell0(npdim+5), zshell0(npdim+5),
     $     xshell1(npdim+5), zshell1(npdim+5),
     $     xshell2(npdim+5), zshell2(npdim+5),
     $     xshell3(npdim+5), zshell3(npdim+5) )

c..  Change origin to the outer midplane. Don't overwrite data arrays 
c     since this subroutine may be called twice.

      xshelld(1:nptsd) = (/ xshelld0(npu:nptsd), xshelld0(1:npu-1) /)
      zshelld(1:nptsd) = (/ zshelld0(npu:nptsd), zshelld0(1:npu-1) /)

      xshell0(1:nptsd) = xshelld(1:nptsd)
      zshell0(1:nptsd) = zshelld(1:nptsd)
      xshell0(nptsd1) = xshell0(1)
      zshell0(nptsd1) = zshell0(1)

      DEALLOCATE ( xshelld, zshelld )

c... Now vertically symmetrize if LSYMWV is 1(upper), 2(lower):

      npts = nptsd
      npts1 = npts + 1
c$$$ 
c$$$      xshell1 = xshell0
c$$$      zshell1 = zshell0
      
      IF ( lsymwv == 1 ) THEN

c.. Symmetry using upper points

c$$$         DO i = 1, npu
c$$$         zshell0(i) =
c$$$     $        -zshell1(nptsd1+1-i)
c$$$         xshell0(i) =
c$$$     $         xshell1(nptsd1+1-i)
c$$$         END DO
c$$$
c$$$         zshell0(npu+1:nptusy) = zshell1(npl+1:nptsd)
c$$$         xshell0(npu+1:nptusy) = xshell1(npl+1:nptsd)
c$$$
c$$$         xshell0(nptusy1) = xshell0(1)
c$$$         zshell0(nptusy1) = zshell0(1)

c.. Why doesn't the following work????
c$$$         
c$$$         zshell0(1:nptusy) =
c$$$     $        (/ -zshell0(nptsd1:npl:-1), zshell0(npl+1:nptsd) /)
c$$$         xshell0(1:nptusy) =
c$$$     $        (/  xshell0(nptsd1:npl:-1), xshell0(npl+1:nptsd) /)

         zshell0(1:nptusy) =
     $        (/ (-zshell0(i),i=nptsd1,npl,-1), zshell0(npl+1:nptsd) /)
         xshell0(1:nptusy) =
     $        (/  (xshell0(i),i=nptsd1,npl,-1), xshell0(npl+1:nptsd) /)
         xshell0(nptusy1) = xshell0(1)
         zshell0(nptusy1) = zshell0(1)
         
         npts = nptusy
         npts1 = npts + 1

      END IF

      IF ( lsymwv == 2 ) THEN

c.. Symmetry using lower points

         zshell0(1:nptlsy) =
     $        (/ zshell0(1:npl), (-zshell0(i),i=npl-1,2,-1) /)
         xshell0(1:nptlsy) =
     $        (/ xshell0(1:npl),  (xshell0(i),i=npl-1,2,-1) /)
         xshell0(nptlsy1) = xshell0(1)
         zshell0(nptlsy1) = zshell0(1)
         
         npts = nptlsy
         npts1 = npts + 1

      END IF

      CALL smooth1 ( xshell0, xshell1, npts1 )
      CALL smooth1 ( zshell0, zshell1, npts1 )

      CALL smooth1 ( xshell1, xshell2, npts1 )
      CALL smooth1 ( zshell1, zshell2, npts1 )

      CALL smooth1 ( xshell2, xshell3, npts1 )
      CALL smooth1 ( zshell2, zshell3, npts1 )

c$$$      WRITE ( 23, '("#",9x, "Xshell_0",5x,"Zshell_0" ,5x,
c$$$     $     "Xshell_1",5x,"Zshell_1",5x,"Xshell_2",5x,"Zshell_2", 5x,
c$$$     $     "Xshell_3",5x,"Zshell_3")' )
c$$$
c$$$      DO i = 1, npts1
c$$$         WRITE ( 23, '(2x, i4, 8es13.5 )' )
c$$$     $        i, xshell0(i), zshell0(i), xshell1(i), zshell1(i),
c$$$     $        xshell2(i), zshell2(i), xshell3(i), zshell3(i)
c$$$      END DO

      WRITE ( ioshel, '("#",9x, "Xshell_0",5x,"Zshell_0" ,5x,
     $     "Xshell_1",5x,"Zshell_1",5x,"Xshell_2",5x,"Zshell_2", 5x,
     $     "Xshell_3",5x,"Zshell_3")' )

      DO i = 1, npts1
         WRITE ( ioshel, '(2x, i4, 8es13.5 )' )
     $        i, xshell0(i), zshell0(i), xshell1(i), zshell1(i),
     $        xshell2(i), zshell2(i), xshell3(i), zshell3(i)
      END DO

c.. Now interpolate thrice smoothed data to MTH points:

      ALLOCATE ( xwal(mth+5), zwal(mth+5), angwal(mth+5),
     $     thgrd(mth+5) )
      ALLOCATE ( xwalq(mth+5), zwalq(mth+5), angwalq(mth+5),
     $     ellwq(mth+5),  thgrq(mth+5) )

      CALL transdx ( xshell3,npts, xwal,mth, 0.0, 0.0 )
      CALL transdx ( zshell3,npts, zwal,mth, 0.0, 0.0 )

      WRITE ( ioshel,
     $     '(//,"#",10x, "Interpolating to MTH points.",/,"#")' )
      WRITE ( ioshel, '("#",10x, "Amtn", 7x, "Xwal-mth",5x,"Zwal-mth",
     $     5x,
     $     "Angwal",/,"#")' )

c Scale mth grid to the npts abscissae. Get the angle that the wall
c   points subtends form the center of the vessel.

      DO i = 1, mth1
         amtnp = (i-1)*REAL(npts)/mth + 1
         angwal(i) =  - ATAN2 ( zwal(i) - zcensh,xwal(i)-xcensh )
         if ( angwal(i) < 0.0 ) angwal(i) = angwal(i) + 2.0*pi
         angwal(i) = 180.0 * angwal(i) / pi
         WRITE ( ioshel, '(2x, i4, 4es13.5 )' ) i,amtnp, xwal(i),
     $        zwal(i),
     $        angwal(i)
      END DO
c.. Ensure monoticity of the angles
      DO i = 1, mth/8
         j = mth1-i+1
         IF ( angwal(i) > 90.0 ) angwal(i) = angwal(i) - 360.0
         IF ( angwal(j) < 90.0 ) angwal(j) = 360.0 - angwal(j)
      END DO


c Now read in Mirnov loops data. Comment out for now.

      CALL lao ( ioshel )

c. Calling EQARCW here gives a "floating invalid" error. 
c  Need to find out why. Fortunately, don't need it.
c$$$      GO TO 1234
c$$$
c$$$c.. Interpolate to equal-arcs grid. Use thdum as dummy
c$$$
c$$$      CALL eqarcw ( xwal,zwal, xwalq,zwalq, ellwq, thgrq, thdum, mth1 )
c$$$
c$$$      WRITE ( 6, '(//,"#",10x, "Interpolating to Equal arcs.",/,"#")' )
c$$$      WRITE ( ioshel,
c$$$     $     '(//,"#",10x, "Interpolating to Equal arcs.",/,"#")' )
c$$$      WRITE ( ioshel, '("#",10x, "Amtn", 7x, "Xwalq-mth",4x,"Zwalq-mth",
c$$$     $     5x,
c$$$     $     "Ellwq", 8x, "Angwal",/,"#")' )
c$$$
c$$$c Scale mth grid to the npts abscissae. Get the angle that the wall
c$$$c   points subtends form the center of the vessel.
c$$$
c$$$      DO i = 1, mth1
c$$$         amtnp = (i-1)*REAL(npts)/mth + 1
c$$$         angwalq(i) =  - ATAN2 ( zwalq(i) - zcensh,xwalq(i)-xcensh )
c$$$         if ( angwalq(i) < 0.0 ) angwalq(i) = angwalq(i) + 2.0*pi
c$$$         angwalq(i) = 180.0 * angwalq(i) / pi
c$$$         WRITE ( ioshel, '(2x, i4, 5es13.5 )' ) i,amtnp, xwalq(i),
c$$$     $        zwalq(i),
c$$$     $        ellwq(i), angwalq(i)
c$$$      END DO
c$$$
c$$$c.. Ensure monoticity of the angles
c$$$      DO i = 1, mth/8
c$$$         j = mth1-i+1
c$$$         IF ( angwalq(i) > 90.0 ) angwalq(i) = angwalq(i) - 360.0
c$$$         IF ( angwalq(j) < 90.0 ) angwalq(j) = 360.0 - angwalq(j)
c$$$      END DO
c$$$
c$$$c.. Interpolate to marc mang1 loops at equal angle points.
c$$$
c$$$      nangl = 32
c$$$      nangl1 = nangl + 1
c$$$      
c$$$      ALLOCATE ( anglp(nangl+5), angx(nangl+5), angz(nangl+5),
c$$$     $     angxz(nangl+5) )
c$$$
c$$$      thgrd(1:mth1) = (/ (REAL(i-1)/mth,i=1,mth1) /)
c$$$c      PRINT '((10es11.3),/))', (thgrd(i), i = 1, mth1)
c$$$
c$$$      angwal(mth1) = 1.0
c$$$      DO i = 1, nangl
c$$$         agrdi = REAL(i-1)*360.0/nangl
c$$$         CALL lag ( angwal,thgrd,mth1, 4, agrdi, anglp(i), da, 0 )
c$$$      END DO
c$$$      anglp(nangl1) = 1.0
c$$$      
c$$$      DO i = 1, nangl
c$$$         CALL lagpe4 ( xwal,mth, anglp(i), angx(i), da, 0 )
c$$$         CALL lagpe4 ( zwal,mth, anglp(i), angz(i), da, 0 )
c$$$         angxz(i) =  - ATAN2 ( angz(i) - zcensh,angx(i)-xcensh )
c$$$         if ( angxz(i) < 0.0 ) angxz(i) = angxz(i) + 2.0*pi
c$$$         angxz(i) = 180.0 * angxz(i) / pi
c$$$      END DO
c$$$
c$$$      angx(nangl1) = angx(1)
c$$$      angz(nangl1) = angz(1)
c$$$
c$$$      WRITE ( ioshel, '(//,"#",16x,"Loops at equal angles",/,"#")' )
c$$$
c$$$      DO i = 1, nangl1
c$$$         WRITE ( ioshel, '(2x, i4, 4es13.5 )' ) i, angx(i), angz(i),
c$$$     $        angxz(i), anglp(i)
c$$$      END DO
c$$$
c$$$      DEALLOCATE ( anglp, angx, angz, angxz )
c$$$
 1234 continue

c...Output thrice smoothed data expanded by the expansion factor, 
c     shfx = 1.05

      shfx = 1.05

      xsh_out(1:mth1) = xcensh + shfx*( xwal(1:mth1)-xcensh)
      zsh_out(1:mth1) = zcensh + shfx*( zwal(1:mth1)-zcensh)
      xsh_out(mth2) = xsh_out(2)
      zsh_out(mth2) = zsh_out(2)

      DEALLOCATE ( xshell0, zshell0, xshell1, zshell1, xshell2,
     $     zshell2, xshell3, zshell3, xwal, zwal, angwal )

      DEALLOCATE ( xwalq, zwalq, angwalq, ellwq, thgrq )

c      CLOSE ( UNIT = ioshel )
c      CLOSE ( UNIT = 21 )

      END SUBROUTINE d3d_shell

c.....................................................
      SUBROUTINE lao ( ioshel )
c.........................................................

      CHARACTER(10), DIMENSION(50) :: namp067, namp322,  namp157
      CHARACTER(10), DIMENSION(10) :: nads180, nads157

      REAL, DIMENSION(50) :: xmp067, xmp322, xmp157
      REAL, DIMENSION(50) :: zmp067, zmp322, zmp157
      REAL, DIMENSION(50) :: amp067, amp322, amp157
      REAL, DIMENSION(50) :: smp067, smp322, smp157
      REAL, DIMENSION(10) :: xds180, xds157
      REAL, DIMENSION(10) :: zds180, zds157
      REAL, DIMENSION(10) :: ads180, ads157
      REAL, DIMENSION(10) :: sds180, sds157

      REAL, DIMENSION(:), ALLOCATABLE :: tmpall        
      REAL, DIMENSION(:), ALLOCATABLE :: tmp067, tmp322, tmp157
      REAL, DIMENSION(:), ALLOCATABLE :: tds180, tds157
      
      REAL, DIMENSION(24) :: rsisvs = (/ 
     $     .0019209,.0019209,.0019209,.0036784,.0037424,.0046093,
     $     .0048872,.0055595,.0058032,.0061324,.0055482,.00503212,
     $     .0019209,.0019209,.0019209,.0036784,.0030229,.0038912,
     $     .0063351,.0068921,.0058032,.0061324,.0055482,.00503212 /)

      REAL, DIMENSION(18) ::  turnfc = (/
     $     (58.0, i=1,5), (55.0, i=1,2), 58.0, 55.0,
     $     (58.0, i=1,5), (55.0, i=1,2), 58.0, 55.0 /)

      CHARACTER(6),  DIMENSION(24)::  vsname = (/
     $     'V-1A ','V-2A ','V-3A ','V-4A ','V-5A ','V-6A ',
     $     'V-7A ','V-8A ','V-9A ','V-10A','V-11A','V-12A',
     $     'V-1B ','V-2B ','V-3B ','V-4B ','V-5B ','V-6B ',
     $     'V-7B ','V-8B ','V-9B ','V-10B','V-11B','V-12B' /)

      CHARACTER(10), DIMENSION(41) :: lpname = (/
     $     'PSF1A  ','PSF2A  ','PSF3A  ','PSF4A  ','PSF5A  ',
     $     'PSF6NA ','PSF7NA ','PSF8A  ','PSF9A  ','PSF1B  ',
     $     'PSF2B  ','PSF3B  ','PSF4B  ','PSF5B  ','PSF6NB ',
     $     'PSF7NB ','PSF8B  ','PSF9B  ','PSI11M ','PSI12A ',
     $     'PSI23A ','PSI34A ','PSI45A ','PSI58A ','PSI9A  ',
     $     'PSF7FA ','PSI7A  ','PSF6FA ','PSI6A  ','PSI12B ',
     $     'PSI23B ','PSI34B ','PSI45B ','PSI58B ','PSI9B  ',
     $     'PSF7FB ','PSI7B  ','PSF6FB ','PSI6B  ','PSI89FB',
     $     'PSI89NB' /)

      CHARACTER(12), DIMENSION(68) :: mpnam2 = (/
     $     'MPI11M067 ','MPI1A067  ','MPI2A067  ','MPI3A067  ',
     $     'MPI4A067  ',
     $     'MPI5A067  ','MPI8A067  ','MPI9A067  ','MPI79A067 ',
     $     'MPI7FA067 ',
     $     'MPI7NA067 ','MPI67A067 ','MPI6FA067 ','MPI6NA067 ',
     $     'MPI66M067 ',
     $     'MPI1B067  ','MPI2B067  ','MPI3B067  ','MPI4B067  ',
     $     'MPI5B067  ',
     $     'MPI8B067  ','MPI89B067 ','MPI9B067  ','MPI79B067 ',
     $     'MPI7FB067 ',
     $     'MPI7NB067 ','MPI67B067 ','MPI6FB067 ','MPI6NB067 ',
     $     'MPI8A322  ','MPI89A322 ','MPI9A322  ','MPI79FA322',
     $     'MPI79NA322',
     $     'MPI7FA322 ','MPI7NA322 ','MPI67A322 ','MPI6FA322 ',
     $     'MPI6NA322 ',
     $     'MPI66M322 ','MPI6NB322 ','MPI6FB322 ','MPI67B322 ',
     $     'MPI7NB322 ',
     $     'MPI7FB322 ','MPI79B322 ','MPI9B322  ','MPI89B322 ',
     $     'MPI8B322  ',
     $     'MPI5B322  ','MPI4B322  ','MPI3B322  ','MPI2B322  ',
     $     'MPI1B322  ',
     $     'MPI11M322 ','MPI1A322  ','MPI2A322  ','MPI3A322  ',
     $     'MPI4A322  ',
     $     'MPI5A322  ',
     $     'MPI1U157  ','MPI2U157  ','MPI3U157  ','MPI4U157  ',
     $     'DSL1U180  ','DSL2U180  ','DSL3U180  ','DSL4U157  ' /)

      REAL, DIMENSION(41) :: rsi = (/
     $     0.8929,  0.8936,  0.8950,  0.8932,  1.0106,
     $     2.5110,  2.2850,  1.2517,  1.6885,  0.8929,
     $     0.8928,  0.8933,  0.8952,  1.0152,  2.5090,
     $     2.2970,  1.2491,  1.6882,  0.9247,  0.9247,
     $     0.9247,  0.9247,  0.9621,  1.1212,  1.6978,
     $     2.2090,  2.1933,  2.5010,  2.4347,  0.9247,
     $     0.9247,  0.9247,  0.9611,  1.1199,  1.8559,
     $     2.2150,  2.1933,  2.5010,  2.4328,  1.3828,
     $     1.5771 /)

      REAL, DIMENSION(41) :: zsi = (/
     $     0.1683,  0.5092,  0.8500,  1.1909,  1.4475,
     $     0.3320,  1.0180,  1.5517,  1.5169, -0.1726,
     $     -0.5135, -0.8543, -1.1953, -1.4536, -0.3350,
     $     -1.0100, -1.5527, -1.5075,  0.0000,  0.3429,
     $     0.6858,  1.0287,  1.3208,  1.4646,  1.4646,
     $     1.2130,  1.0629,  0.5270,  0.4902, -0.3429,
     $     -0.6858, -1.0287, -1.3198, -1.4625, -1.4625,
     $     -1.2140, -1.0602, -0.5240, -0.4909, -1.4625,
     $     -1.4625 /)

      REAL, DIMENSION(68) :: xmp2 = (/
     $     0.9720,  0.9717,  0.9712,  0.9717,  0.9714,
     $     1.0432,  1.2615,  1.5340,  1.7861,  2.0514,
     $     2.2055,  2.2855,  2.3685,  2.4157,  2.4168,
     $     0.9718,  0.9711,  0.9712,  0.9708,  1.0490,
     $     1.2530,  1.4660,  1.6987,  1.8978,  2.0887,
     $     2.2057,  2.2914,  2.3697,  2.4164,
     $     1.2194,  1.4017,  1.5841,  1.7825,  1.9243,
     $     2.0673,  2.2066,  2.2873,  2.3681,  2.4159,
     $     2.4182,  2.4162,  2.3695,  2.2897,  2.2029,
     $     2.0848,  1.8940,  1.6989,  1.4769,  1.2541,
     $     1.0479,  0.9724,  0.9740,  0.9748,  0.9737,
     $     0.9732,  0.9741,  0.9742,  0.9745,  0.9722,  
     $     1.0511,
     $     1.489,  1.560,  1.723,  1.873,
     $     1.506,  1.590,  1.730,  1.897 /)

      REAL, DIMENSION(68) :: ymp2 = (/
     $     -0.0043,  0.1813,  0.5090,  0.8705,  1.1437,
     $     1.3223,  1.4072,  1.4101,  1.3233,  1.1068,
     $     0.9041,  0.7070,  0.5107,  0.2521,  0.0008,
     $     -0.1878, -0.5216, -0.8601, -1.1524, -1.3362,
     $     -1.4088, -1.4086, -1.4099, -1.3303, -1.0969,
     $     -0.8977, -0.6911, -0.5037, -0.2435,
     $     1.4055,  1.4070,  1.4084,  1.3230,  1.2064,
     $     1.0904,  0.9015,  0.7036,  0.5079,  0.2492,
     $     -0.0009, -0.2438, -0.5041, -0.6907, -0.8979,
     $     -1.1020, -1.3327, -1.4063, -1.4055, -1.4047,
     $     -1.3302, -1.1589, -0.8542, -0.5122, -0.1870,
     $     -0.0023,  0.1817,  0.5116,  0.8503,  1.1609,
     $     1.3304,
     $     1.250,  1.187,  1.112,  1.116,
     $     1.270,  1.198,  1.142,  1.133 /)

      REAL, DIMENSION(68) :: amp2 = (/
     $    -269.8424, -269.8633, -269.9790,  89.8003, -269.7477,
     $     45.8000,    0.3500,    0.6667, -39.4000,  -39.2167,
     $     -68.1500,  -67.6500,  -67.3667, -89.6333,  -90.0167,
     $     89.9264, -269.9369, -269.9369,  89.6320, -223.5167,
     $     -180.0500, -180.3167, -180.2500, -129.800, -128.9833,
     $     -111.9667, -112.7833, -112.9167, -90.1000,
     $     0.3000,    0.6670,    0.6170, -39.2830,  -39.1670,
     $     -39.4170,  -67.9000,  -67.8000, -67.3670,  -89.3330,
     $     -89.9000,  269.2500,  246.1830, 247.2830,  247.5670,
     $     230.8170,  230.4670,  179.8670, 180.0830,  180.0330,
     $     136.3830,   89.5670,   89.7500,  90.2670,   90.1170,
     $     89.9330,   90.0170,   89.7000,  90.1500,   90.4830,
     $     44.6170,
     $     -40.55,  -40.77,  -3.30,   0.08,
     $     49.53,   49.53,  92.30,  90.50 /)

      REAL, DIMENSION(68) :: smp2 = (/
     $     0.1016,  0.1384,  0.1384,  0.1384,  0.1384,
     $     0.1016,  0.1321,  0.1321,  0.1107,  0.1107,
     $     0.1107,  0.1107,  0.1107,  0.1321,  0.1321,
     $     0.1384,  0.1384,  0.1384,  0.1384,  0.1016,
     $     0.1384,  0.1384,  0.1384,  0.1107,  0.1107,
     $     0.1107,  0.1107,  0.1107,  0.1321,
     $     0.1394,  0.1145,  0.1403,  0.1409,  0.1141,
     $     0.1398,  0.1403,  0.1407,  0.1402,  0.1407,
     $     0.1399,  0.1399,  0.1403,  0.1406,  0.1407,
     $     0.1403,  0.1406,  0.1415,  0.1403,  0.1403,
     $     0.1403,  0.1405,  0.1402,  0.1404,  0.1400,
     $     0.1153,  0.1400,  0.1403,  0.1404,  0.1404,
     $     0.1407,
     $     0.054,   0.027,   0.027,   0.054,
     $     -0.094,  -0.088,  -0.106,  -0.165 /)

      REAL, DIMENSION(60) :: patmp2 = (/
     $     0.1845500,  0.2566500,  0.3446000,  0.3173500,
     $     0.2413664,  0.2202499,  0.2728091,  0.2806501,  0.3106111,
     $     0.3109208,  0.2295793,  0.2129159,  0.2292131,  0.2710449,
     $     0.2478000,  0.2586500,  0.3361500,  0.3154000,  0.2556976,
     $     0.2131110,  0.2341576,  0.2228506,  0.2486520,  0.2551194,
     $     0.2669137,  0.2292701,  0.2133847,  0.2214800,  0.2718855,
     $     0.1838351,  0.1823557,  0.1839847,  0.2282504,  0.1838572,
     $     0.2136775,  0.2273982,  0.2127223,  0.2300735,  0.2691039,
     $     0.2465000,  0.2759520,  0.1962927,  0.2346333,  0.2421203,
     $     0.2582166,  0.2455210,  0.2451018,  0.2224015,  0.2393365,
     $     0.2123904,  0.2517095,  0.3233501,  0.3336000,  0.2549500,
     $     0.1843500,  0.2569500,  0.3343000,  0.3246499,  0.2469479,
     $     0.2170768  /)

      xcensh = 1.6955
      zcensh = 0.0

      pi = acos (-1.0 )

      nmp067 = 29
      nmp322 = 31 
      nmp157 = 4
      nds180 = 3      
      nds157 = 1

      nm0 = 0
      nm1 = nmp067
      nm2 = nm1 + nmp322 
      nm3 = nm2 + nmp157
      nm4 = nm3 + nds180
      nm5 = nm4 + nds157
      nmp = nm5

      namp067(1:nmp067) = mpnam2(nm0+1:nm1)
      namp322(1:nmp322) = mpnam2(nm1+1:nm2)
      namp157(1:nmp157) = mpnam2(nm2+1:nm3)
      nads180(1:nds180) = mpnam2(nm3+1:nm4)
      nads157(1:nds157) = mpnam2(nm4+1:nm5)

      xmp067(1:nmp067) = xmp2(nm0+1:nm1)
      xmp322(1:nmp322) = xmp2(nm1+1:nm2)
      xmp157(1:nmp157) = xmp2(nm2+1:nm3)
      xds180(1:nds180) = xmp2(nm3+1:nm4)
      xds157(1:nds157) = xmp2(nm4+1:nm5)

      zmp067(1:nmp067) = ymp2(nm0+1:nm1)
      zmp322(1:nmp322) = ymp2(nm1+1:nm2)
      zmp157(1:nmp157) = ymp2(nm2+1:nm3)
      zds180(1:nds180) = ymp2(nm3+1:nm4)
      zds157(1:nds157) = ymp2(nm4+1:nm5)

      amp067(1:nmp067) = amp2(nm0+1:nm1)
      amp322(1:nmp322) = amp2(nm1+1:nm2)
      amp157(1:nmp157) = amp2(nm2+1:nm3)
      ads180(1:nds180) = amp2(nm3+1:nm4)
      ads157(1:nds157) = amp2(nm4+1:nm5)

      smp067(1:nmp067) = smp2(nm0+1:nm1)
      smp322(1:nmp322) = smp2(nm1+1:nm2)
      smp157(1:nmp157) = smp2(nm2+1:nm3)
      sds180(1:nds180) = smp2(nm3+1:nm4)
      sds157(1:nds157) = smp2(nm4+1:nm5)

      ALLOCATE ( tmpall(nm5) )
      ALLOCATE ( tmp067(nmp067), tmp322(nmp322), tmp157(nmp157) )
      ALLOCATE ( tds180(nds180), tds157(nds157) )
      

c Get ray angle of loops subtended from the center of the shell:      

      tmpall(1:nm5) = atan2 ( ymp2(1:nm5)-zcensh,xmp2(1:nm5)-xcensh ) 
      tmpall = 180.0 * tmpall / pi

c.. Now separate into categories:

      tmp067(1:nmp067) = tmpall(nm0+1:nm1)
      tmp322(1:nmp322) = tmpall(nm1+1:nm2)
      tmp157(1:nmp157) = tmpall(nm2+1:nm3)
      tds180(1:nds180) = tmpall(nm3+1:nm4)
      tds157(1:nds157) = tmpall(nm4+1:nm5)
         
      write ( ioshel, '(//,"#",30x,"Mirnov Loops" )' )

      write ( ioshel, '("#",/,"#",8x,"Name",11x, "X", 12x, "Z", 9x,
     $     "Orient.", 9x, "L", 7x, "Ray Angle" )' )
      
      DO i = 1, nmp067
         write ( ioshel, '(1x, i4, 2x, a, 4es13.5, es11.3)' )
     $        i, namp067(i),
     $        xmp067(i), zmp067(i), amp067(i), smp067(i), tmp067(i)
      END DO
      
      WRITE ( ioshel, '(/)' )

      DO i = 1, nmp322
         write ( ioshel, '(1x, i4, 2x, a, 4es13.5, es11.3)' )
     $        i, namp322(i),
     $        xmp322(i), zmp322(i), amp322(i), smp322(i), tmp322(i)
      END DO

      WRITE ( ioshel, '(/)' )

      DO i = 1, nmp157
         write ( ioshel, '(1x, i4, 2x, a, 4es13.5, es11.3)' )
     $        i, namp157(i),
     $        xmp157(i), zmp157(i), amp157(i), smp157(i), tmp157(i)
      END DO

      WRITE ( ioshel, '(/)' )

      DO i = 1, nds180
         write ( ioshel, '(1x, i4, 2x, a, 4es13.5, es11.3)' )
     $        i, nads180(i),
     $        xds180(i), zds180(i), ads180(i), sds180(i), tds180(i)
      END DO

      WRITE ( ioshel, '(/)' )

      DO i = 1, nds157
         write ( ioshel, '(1x, i4, 2x, a, 4es13.5, es11.3)' )
     $        i, nads157(i),
     $        xds157(i), zds157(i), ads157(i), sds157(i), tds157(i)
      END DO

      call atpoint ( "End of Lao","nmp157", nmp157,
     $     "tmpall(1)",tmpall(1), 6, 23 )

      DEALLOCATE ( tmpall, tmp067, tmp322, tmp157, tds180, tds157 )

      END SUBROUTINE lao

c.........................................................................
      subroutine d3dwall ( xwall, ywall, mth, iomod, iotty )
c........................................................................
c
c...This from Alan Turnbull for parameterizing the DIID vessel.
c
c               For use in VACUUM:
c    Replace xwall,ywall by xwal,zwal.
c            mth = nwalp
c            ncdf = 31
c       Try  zlim = 0.0
c
      parameter ( ncdf = 26 )
      dimension xwall(*), ywall(*), rwi(ncdf), zwi(ncdf)
c    End of for use in VACUUM
c
c --------------------------------------------
c initialize DIII-D wall fourier coefficients
c --------------------------------------------
      data nwcoef/ncdf/
      data rwall0/ 1.6400000/, zwall0/ 0.0000000/, awall0/ 0.8839410/
     &    ,ewall0/ 1.4037020/
      data (rwi(k),k=1,ncdf)
     &    / 0.1000000e+01, 0.5526794e-01,-0.1738114e+00, 0.1850757e-01
     &    , 0.3714965e-01,-0.2882647e-01,-0.2357329e-02, 0.9548103e-02
     &    ,-0.1214923e-01,-0.1853416e-02, 0.6837493e-02,-0.1711245e-02
     &    , 0.2270762e-02, 0.3689963e-02,-0.3959393e-02,-0.1098017e-02
     &    , 0.3745465e-02,-0.2157904e-03,-0.3977743e-03,-0.2725623e-03
     &    ,-0.1005857e-02,-0.4579016e-05, 0.2396789e-02,-0.7057043e-03
     &    , 0.1158347e-02, 0.3552319e-03/
c
      data (zwi(k),k=1,ncdf)
     &    / 0.1000000e+01,-0.3236632e-01,-0.1629422e+00, 0.6013983e-01
     &    , 0.1167756e-01,-0.2579542e-01, 0.1626464e-01,-0.2085857e-02
     &    ,-0.9098639e-02, 0.1022163e-01,-0.4388253e-02,-0.9367258e-02
     &    , 0.8308497e-02, 0.4765150e-02,-0.4611675e-02,-0.1121423e-02
     &    ,-0.2501100e-03, 0.4282634e-03, 0.2669702e-02,-0.1073800e-02
     &    ,-0.2191338e-02, 0.1328267e-02, 0.5050959e-03,-0.5758863e-03
     &    , 0.9348883e-03, 0.7094351e-03/
c
c 4.1 Construct the wall position and place in (xwall,ywall)
c
c 4.1.1 Initialize the wall parameters
c       The starting point is at the point corresponding to
c       (xlim,zlim)
c
c$$        nwalp     = min0(ncv,nwp-1)
        rwll      = rwall0
        zwll      = zwall0
        awll      = awall0*rext
        ewll      = ewall0
        zstart    = zlim
c
c Lines for VACUUM to overwrite some variables
      zlim = 0.0
      zstart = 0.0
      nwalp = mth
      rext = 1.0
      awll = awall0*rext
c End of lines for VACUUM
c
c 4.1.2 Construct the wall position
c
        call d3dvesl(rwll,zwll,awll,ewll,rwi,zwi,nwcoef,zstart
     &              ,xwall,ywall,nwalp,ier)
c
        write ( iomod, '("ier in d3dwall = ", i3)' ) ier
        write ( iotty,  '("ier in d3dwall = ", i3)' ) ier
        xwall(mth+1) = xwall(1)
        ywall(mth+1) = ywall(1)
        xwall(mth+2) = xwall(2)
        ywall(mth+2) = ywall(2)
c        
        return
        end
c
c...................................................................
      subroutine d3dvesl(r0,z0,a0,e0,ar,az,nval,zst,r,z,npts,ier)
c...................................................................
c
c ------------------------------------------------------------------
c
c     Compute the DIII-D vessel position and store in (r(j),z(j))
c
c ------------------------------------------------------------------
c
      dimension ar(nval),az(nval)
      dimension r (npts),z (npts)
c
c
c -----------------------------------------------------------
c  Define the starting angle to coincide with z = zst on the
c  outboard side
c -----------------------------------------------------------
c
      pii    =  3.1415926535897932385
      ier    = 0
c      if(z0 .eq. zst) then
      if ( abs(z0-zst) .le. 1.0e-6 ) then
c        arcst    = 0.0
      else
        if(z0 .lt. zst) isgn  = +1
        if(z0 .gt. zst) isgn  = -1
        zpcmp    = (zst-z0)/(a0*e0)
c
        arci     = 0.0
        arcf     = 0.5*pii
        dfi      = (arcf-arci)/npts
c
c
        arca     = arci
        zza      = 0.0
        do 20g j  = 2,npts
        arcb     = arca + isgn*dfi
c
        zzb      = 0.0
        do  5 k  = 1,nval
        ackb     = k*arcb
        zzb      = zzb + az(k)*sin(ackb)
  5     continue
c
        if((zza-zpcmp)*(zzb-zpcmp) .le. 0.0) then
          arcc     = arcb + isgn*dfi
          zzc      = 0.0
          do 10 k  = 1,nval
          ackc     = k*arcc
          zzc      = zzc + az(k)*sin(ackc)
  10      continue
          go to 25
        else
          arca     = arcb
          zza      = zzb
        end if
  20    continue
        ier    = 1
        return
c
  25    continue
        dzp    = zzc - zzb
        dzm    = zzb - zza
        dzt    = zzc - zza
        dcf1   = dfi*(dzm/dzp + dzp/dzm)/dzt
        dcf2   = dfi*(1.0/dzp - 1.0/dzm)/dzt
        zdf    = zpcmp - zzb
        arcst  = arcb + dcf1*zdf + dcf2*zdf*zdf
      end if
c
c
c
c --------------------------------------------------------------
c  Loop through npts equal angles and compute the corresponding
c  wall position
c  The loop is anticlockwise from the point where z = zst
c  and the final periodic point is omitted
c  The angle arc increases with index j:  darc is positive
c --------------------------------------------------------------
c
      arc0     =  arcst
      arc1     =  arcst + 2.0*pii
      darc     = (arc1-arc0)/npts
c
      do 100 j = 1,npts
      arc      = arc0 + (j-1.)*darc
c
c  Sum the fourier contributions to the wall shape
      sumr     = 0.0
      sumz     = 0.0
      do 50 k  = 1,nval
      arck     = k*arc
      sumr     = sumr + ar(k)*cos(arck)
      sumz     = sumz + az(k)*sin(arck)
  50  continue
c
c  Define the wall position
      rpval    = r0  +    a0*sumr
      zpval    = z0  - e0*a0*sumz
      r(j)     = rpval
      z(j)     = zpval
 100  continue
c
      return
      end
c......................................................................
      SUBROUTINE waleig
c     ......................................................................
c     
      include 'vacuum1.inc'
      include 'vacuum5.inc'
      include 'vacuum7.inc'
      include 'vacuum8.inc'
      dimension wkxx(nths), wkyy(nths)
c     
c-----------------------------------------------------------------------
c     reads eigenarrays from wal-eig..dat
c-----------------------------------------------------------------------
c     
      CHARACTER(5) :: za
      CHARACTER(7) :: za7
      CHARACTER(12) :: ardrz, ardre, ardnp, ardne
      CHARACTER(17) :: ardtr

      nout1 = iotty
      nout2 = 26
c     
c     nvl = 100
c     nvl = 42
c     nvc = 50
c     nvc = 100
c     nvc = 300
c     nvc1 = nvc + 1
c     
c     lwmn = lmin(1)
c     lwmx = lmax(1)
      jwal1 = lwmx - lwmn + 1
      jwal2 = 2*jwal1
      jwal4 = 2*jwal2

c...  Abscissa for l values on wall
      do j1 = 1, jwal2
         xlwm(j1) = j1 
      end do

c-----------------------------------------------------------------------
c     read data
c-----------------------------------------------------------------------

      OPEN ( UNIT = 25, FILE = 'wal_eig.dat')
      OPEN ( UNIT = 26, FILE = 'wal_eig_o.dat', STATUS = "REPLACE")

      WRITE ( NOUT1, '(/"READING WAL_EIG.DAT",/)' )
      WRITE ( NOUT2, '(/"reading wal_eig.dat",/)' )
      WRITE ( OUTMOD, '(/"reading wal_eig.dat",/)' )
c      read(25,'(//)')

      READ (25,'(a12, 1p2e13.5)') ardrz, rcen, zcen
      READ (25,'(a12, 1p2e13.5)') ardre, radius,elong
      READ (25,'(a17, 1p2e13.5)') ardtr, trian, achance

      WRITE (NOUT1,'(a, 1p2e13.5)') ardrz, rcen, zcen
      WRITE (NOUT1,'(a, 1p2e13.5)') ardre, radius, elong
      WRITE (NOUT1,'(a, 1p2e13.5)') ardtr, trian, achance

      WRITE (NOUT2,'(a, 1p2e13.5)') ardrz, rcen, zcen
      WRITE (NOUT2,'(a, 1p2e13.5)') ardre, radius, elong
      WRITE (NOUT2,'(a, 1p2e13.5)') ardtr, trian, achance

      WRITE (OUTMOD,'(a, 1p2e13.5)') ardrz, rcen, zcen
      WRITE (OUTMOD,'(a, 1p2e13.5)') ardre, radius, elong
      WRITE (OUTMOD,'(a, 1p2e13.5)') ardtr, trian, achance

      READ (25,'(a12, i10)') ardnp, nvc1
      READ (25,'(a12, i10)') ardne, nvlh

      WRITE (NOUT1,'(a, i10)') ardnp, nvc1
      WRITE (NOUT2,'(a, i10)') ardnp, nvc1
      WRITE (OUTMOD,'(a, i10)') ardnp, nvc1

      WRITE (NOUT1,'(a, i10)') ardne, nvlh
      WRITE (NOUT2,'(a, i10)') ardne, nvlh
      WRITE (OUTMOD,'(a, i10)') ardne, nvlh

      nvc = nvc1 - 1
      nvl = 2 * nvlh

      DO ird = 1, nvc1
         READ ( 25,'(a7,i3,1p3e13.5)') za7,i,zxs,zzs,zts
      END DO

      DO ivl = 1, nvl
c     read(25,*) za, neig(ivl), egval(ivl)
         READ(25,'(a5,i3, 1pe13.5)') za, neig(ivl), rwval(ivl)
         READ(25,*) ( rwvec(ivc,ivl),ivc=1,nvc1 )
         rwvec(ivc,ivl) = rwvec(ivc,ivl)
         WRITE(26,*) za, neig(ivl), rwval(ivl)
         WRITE(26,'(1p6e13.5)') ( rwvec(ivc,ivl),ivc=1,nvc1 )
      END DO
      
c     Rearrange to separate even and odd solutions.
c     
c$$$  c     Leave DC eigenmode in first address.
c$$$  nwk(1) = neig(1)
c$$$  rwvl(1) = rwval(1)
c$$$  do ivc = 1, nvc
c$$$  rwvc(ivc,1) = rwvec(ivc,1)
c$$$  end do
c     
c..   Even modes first
      icnte = 0
      do ivl = 1, nvl
         if ( rwvec(2,ivl)*rwvec(nvc,ivl) .gt. 0.0 ) then
            icnte = icnte + 1
            nwke(icnte) = neig(ivl)
            rwvle0(icnte) = rwval(ivl)
            do ivc = 1, nvc
               rwvce0(ivc,icnte) = rwvec(ivc,ivl)
            end do
            rwvce0(nvc+1,icnte) = rwvce0(1,icnte)
            rwvce0(nvc+2,icnte) = rwvce0(2,icnte)
         end if
      end do
c..   Now odd modes
      icnto = 0
      do ivl = 1, nvl
         if ( rwvec(2,ivl)*rwvec(nvc,ivl) .le. 0.0 ) then
            icnto = icnto + 1
            nwko(icnto) = neig(ivl)
            rwvlo0(icnto) = rwval(ivl)
            do ivc = 1, nvc
               rwvco0(ivc,icnto) = rwvec(ivc,ivl)
            end do
            rwvco0(nvc+1,icnto) = rwvco0(1,icnto)
            rwvco0(nvc+2,icnto) = rwvco0(2,icnto)
         end if
      end do
c     
c$$$  do ivl = 1, nvl
c$$$  write(26,*) za, neig(ivl), rwval(ivl)
c$$$  write(26,'(1p6e13.5)') ( rwvec(ivc,ivl),ivc=1,nvc )
c$$$  end do
c     
      write ( 26,'(/,20x,"Even-Odd Reordering",/)' )
      write ( 26,'(/,20x,"Even Modes",/)' )
      do ivl = 1, icnte
         write(26,*) ivl, "  ", za, nwke(ivl), rwvle0(ivl)
         write(26,'(1p6e13.5)') ( rwvce0(ivc,ivl),ivc=1,nvc1 )
      end do
      write ( 26,'(/,20x,"Odd Modes",/)' )
      do ivl = 1, icnto
         write(26,*) ivl, "  ", za, nwko(ivl), rwvlo0(ivl)
         write(26,'(1p6e13.5)') ( rwvco0(ivc,ivl),ivc=1,nvc1 )
      end do
c     
      write ( outmod,'(/,20x,"Eigenvalues: Even-Odd Reordering",/)' )
      write ( outmod,'(/,20x,"Even Modes",/)' )
      do ivl = 1, icnte
         write(outmod,*) ivl, "  ", za, nwke(ivl), rwvle0(ivl)
      end do
      write ( outmod,'(/,20x,"Odd Modes",/)' )
      do ivl = 1, icnto
         write(outmod,*) ivl, "  ", za, nwko(ivl), rwvlo0(ivl)
      end do

c     
c..   Orthogonality checks:
c     . Use rwvec for storage of result:
c$$$  call orthchk( rwvce0,nths,nwm, rwvce0,nths,nwm, rwvec,nths,nwm2,
c$$$  $     21,nvc,21, 21,21, "rwvce0-or", 0, outmod )
c$$$  call orthchk( rwvco0,nths,nwm, rwvco0,nths,nwm, rwvec,nths,nwm2,
c$$$  $     21,nvc,21, 21,21, "rwvco0-or", 0, outmod )
c$$$  call orthchk( rwvce0,nths,nwm, rwvco0,nths,nwm, rwvec,nths,nwm2,
c$$$  $     21,nvc,21, 21,21, "rwvce0-rwvco0-or", 0, outmod )
c$$$  c
c$$$  call evplts ( rwvce0, nths,nwm, rwvle0,ff, 1,nvc1, 1,icnte,
c$$$  $     jobid, "Wall Modes" )
c$$$  call evplts ( rwvco0, nths,nwm, rwvlo0,ff, 1,nvc1, 1,icnto,
c$$$  $     jobid, "Wall Modes" )
c     
c...  Interpolate to mth points. Multiply by sqrt(J grad-psi)
c     Use workg1, workg2 for workspace.
c     
c     Even modes..
      do ivl = 1, icnte
         do ivc = 1, nvc
            workl(ivc) = rwvce0(ivc,ivl)
            call trans ( workl,nvc, workl2,mth )
         end do
         do i = 1, mth
            rwvce0(i,ivl) = workl2(i)
c     if ( ivl .eq. 1) rwvce0(i,ivl) = - workl2(i)
         end do
         rwvce0(mth+1,ivl) = rwvce0(1,ivl)
         rwvce0(mth+2,ivl) = rwvce0(2,ivl)
      end do
c     Odd modes..
      do ivl = 1, icnto
         do ivc = 1, nvc
            workl(ivc) = rwvco0(ivc,ivl)
            call trans ( workl,nvc, workl2,mth )
         end do
         do i = 1, mth
            rwvco0(i,ivl) = workl2(i)
         end do
         rwvco0(mth+1,ivl) = rwvco0(1,ivl)
         rwvco0(mth+2,ivl) = rwvco0(2,ivl)
      end do

c$$$c..Experiment:  Multiply by XWAL to test for Ming's norm. 
c$$$      DO ivl = 1, icnte
c$$$         rwvce0(1:mth+2,ivl) = rwvce0(1:mth+2,ivl)*sqrt(xwal(1:mth+2))
c$$$      END DO
c$$$      DO ivl = 1, icnto
c$$$         rwvco0(1:mth+2,ivl) = rwvco0(1:mth+2,ivl)*sqrt(xwal(1:mth+2))
c$$$      END DO
     
c...  Fourier Analyse the Wall Eigenfunctions.
      write ( 26,'(/,20x,"Fourier Analyze the Wall Eigenfuncts.",/)' )
c     
      pii = 1.0 / pye
      jdel1 = 21
      jdel2 = 21
      jmax1 = lmax(1) - lmin(1) + 1
      jwal10 = 21
      rnge = twopi
      zzfac = pii * rnge / mth
c...  Replace by Fortran 90:
c$$$  call foranv2 ( rwvce0,nths,nwm,jwal10, wrkvr,nfm,nfm,
c$$$  $     coslt,nths,nfm,jmax1, mth, 0,0, rnge,pii )
c$$$  call foranv2 ( rwvce0,nths,nwm,jwal10, wrkvi,nfm,nfm,
c$$$  $     sinlt,nths,nfm,jmax1, mth, 0,0, rnge,pii )
      wrkvr(1:jwal10,1:jmax1) = zzfac *
     $     matmul ( transpose(rwvce0(1:mth,1:jwal10)),
     $     coslt(1:mth,1:jmax1) )
      wrkvi(1:jwal10,1:jmax1) = zzfac *
     $     matmul ( transpose(rwvce0(1:mth,1:jwal10)),
     $     sinlt(1:mth,1:jmax1) )
c$$$  call matwrtn ( wrkvr, nfm,nfm, ln,ln,jwal10,jwal10,jdel1,jdel2,
c$$$  $     "Cosine of Evens", outmod,0 )
c$$$  call matwrtn ( wrkvi, nfm,nfm, ln,ln,jwal10,jwal10,jdel1,jdel2,
c$$$  $     "Sine of Evens", outmod,0 )
c     
c$$$  call foranv2 ( rwvco0,nths,nwm,jwal10, wrkvr,nfm,nfm,
c$$$  $     coslt,nths,nfm,jmax1, mth, 0,0, rnge,pii )
c$$$  call foranv2 ( rwvco0,nths,nwm,jwal10, wrkvi,nfm,nfm,
c$$$  $     sinlt,nths,nfm,jmax1, mth, 0,0, rnge,pii  )
      wrkvr(1:jwal10,1:jmax1) = zzfac *
     $     matmul ( transpose(rwvco0(1:mth,1:jwal10)),
     $     coslt(1:mth,1:jmax1) )
      wrkvi(1:jwal10,1:jmax1) = zzfac *
     $     matmul ( transpose(rwvco0(1:mth,1:jwal10)),
     $     sinlt(1:mth,1:jmax1) )

c$$$  call matwrtn ( wrkvr, nfm,nfm, ln,ln,jwal10,jwal10,jdel1,jdel2,
c$$$  $     "Cosine of Odds", outmod,0 )
c$$$  call matwrtn ( wrkvi, nfm,nfm, ln,ln,jwal10,jwal10,jdel1,jdel2,
c$$$  $     "Sine of Odds", outmod,0 )
c     
c..   Orthogonality checks:
c     . Use rwvec for storage of result:
      call orthchk( rwvce0,nths,nwm, rwvce0,nths,nwm, xwal,
     $     rwvec,nths,nwm2,
     $     21,mth,21, 21,21, "rwvce0-or", 0, outmod )
      call orthchk( rwvco0,nths,nwm, rwvco0,nths,nwm, xwal,
     $     rwvec,nths,nwm2,
     $     21,mth,21, 21,21, "rwvco0-or", 0, outmod )
      call orthchk( rwvce0,nths,nwm, rwvco0,nths,nwm, xwal,
     $     rwvec,nths,nwm2,
     $     21,mth,21, 21,21, "rwvce0-rwvco0-or", 0, outmod )
c     
      call evplts ( rwvce0, nths,nwm, rwvle0,ff, 1,mth1, 1,icnte,
     $     jobid, "Wall Modes-E" )
      call evplts ( rwvco0, nths,nwm, rwvlo0,ff, 1,mth1, 1,icnto,
     $     jobid, "Wall Modes-O" )
c     
c..   Construct Functions for negative l_w values.
c.... Bypass for Now
c     
c$$$  do i = 1, mth1
c$$$  do jw = 1, jwal1
c$$$  jj = lwmn - 1 + jw
c$$$  ijj = iabs(jj)
c$$$  rwvce(i,jw) = rwvce0(i,ijj+1)
c$$$  rwvle(jw) = rwvle0(ijj+1)
c$$$  if ( jj .ne. 0 ) then
c$$$  rwvco(i,jw) = isign(1,jj)* rwvco0(i,ijj)
c$$$  rwvlo(jw) = rwvlo0(ijj)
c$$$  else
c$$$  rwvco(i,jw) = 0.0
c$$$  rwvlo(jw) = 1.0
c$$$  end if
c$$$  end do
c$$$  end do
c..   Construct Functions for negative l_w values.
c..   Use same indexing for positive non zero values.
c     
      do i = 1, mth1
         do jw = 1, jwal1
            jj = lwmn - 1 + jw
            ijj = iabs(jj)
            rwvce(i,jw) = rwvce0(i,ijj)
            rwvle(jw) = rwvle0(ijj)
            if ( jj .ne. 0 ) then
               rwvco(i,jw) = isign(1,jj)* rwvco0(i,ijj)
               rwvlo(jw) = rwvlo0(ijj)
            else
               rwvco(i,jw) = 0.0
               rwvlo(jw) = 1.0
            end if
         end do
      end do
c     
      CLOSE ( UNIT = 25 )
      CLOSE ( UNIT = 26 )
c     
c-----------------------------------------------------------------------
c     terminate routine.
c-----------------------------------------------------------------------
      RETURN
      END

c.....................................................................
      SUBROUTINE fbcoil2 (ctip, ntip, xcw,zcw, xcwp, zcwp,
     $     ccoil, cgrd, zlc, zlen, chgt, cname ) 
c......................................................................
c     
c... Heaviside distribution for the coils in CCOIL.  Parameterize it with 
c    distance along the coil.  (x,z) coordnates of the coil tips input in 
c    ctip(x1,z1,x2,z2). 
c    ZANGLE is the angle of the normal 
c    wrt to the x-axis. ZLC is length along the coil. ZLEN is the 
c    length of the coil.

c     Outputs: xcw, zcw, chgt, ccoil
c     
c$$$      INCLUDE 'vacuum1.inc'
c$$$      INCLUDE 'vacuum2.inc'
c$$$      INCLUDE 'vacuum3.inc'
c$$$      INCLUDE 'vacuum5.inc'
c$$$      INCLUDE 'vacuum8.inc'
c     
      DIMENSION ctip(*), xcw(*), zcw(*),xcwp(*), zcwp(*),
     $     ccoil(*), cgrd(*), zlc(*), chgt(*)
      CHARACTER*(*) cname

      pi = ACOS ( -1.0 )
      twopi = 2.0 * pi

      zd = 1.0 / ( ntip - 1 )
      zdx = ( ctip(3) - ctip(1) ) * zd
      zdz = ( ctip(4) - ctip(2) ) * zd
      zangle = ATAN2 ( zdz, zdx ) + pi / 2.0
      IF ( zangle > twopi ) zangle = zangle - twopi
      zsang = SIN(zangle)
      zcang = COS(zangle)
      
      DO i = 1, ntip
         xcw(i) = ctip(1) + (i-1)*zdx
         zcw(i) = ctip(2) + (i-1)*zdz
         cgrd(i) = (i-1) * zd
         ccoil(i) = 1.0
      END DO

c...  Multiply coils by mu_0 here. Careful that we don't need the current
c       though, since only the fields should be multiplied by mu_0.

c$$$      zmu0 =  4.0 * pi * 1.0e-7
c$$$ 
c$$$      ccoil(1:ntip) = zmu0 * ccoil(1:ntip)

      zlc(1:ntip) = ( xcw(1:ntip)-xcw(1) ) * zsang -
     $     ( zcw(1:ntip)-zcw(1) ) * zcang
      zlen = zlc(ntip) - zlc(1)
      chgt(1:4) = ctip(1:4)
      
      xcwp(1:ntip) =   zsang
      zcwp(1:ntip) = - zcang
      
      WRITE ( 28, '(/, 10x, a, "... in subroutine fbcoil2" )') cname
      WRITE ( 28, '( 1x, "Npoints =", i5, 2x, "ZANGLE = ", es12.4 )' )
     $      ntip, zangle*180.0/pi
      WRITE ( 28, '( 1X, "ZLC = ", /, (8ES12.4) )') ( zlc(i),i=1,ntip ) 

      END SUBROUTINE fbcoil2
c.....................................................................
      SUBROUTINE fbcoil ( angfc, taufc, centfc, mthfc, ncgr, iopfc,
     $     icshp, afc, bfc, dfc, drfc, dzfc, xcw,zcw, ccoil,
     $     cgrd, chgt ) 
c......................................................................
c     
c     Outputs: xcw, zcw, chgt, ccoil
c     
      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum5.inc'
      INCLUDE 'vacuum8.inc'
c     
      REAL nq
      DIMENSION xcw(*), zcw(*), ccoil(*), cgrd(*), chgt(*)
      DIMENSION xpla(nths),zpla(nths)

      REAL, DIMENSION(:), ALLOCATABLE :: ellc, ellqc
c      REAL, DIMENSION(:), ALLOCATABLE :: ellqc, qtgrc, ellc, cgrd
c     
      COMMON / bigv1 / grdgre(nths2,nths2), grwp(nths,nths),
     $     grri(nths2,nfm2)
c     
      EQUIVALENCE (xpla,xinf), (zpla,zinf)


c..   Put coil parametrization in wxgrd, wzgrd. Transform to equal arc
c     in xcw, zcw.
c     
c     The input angles are in units of pi for some reason.

      mthfc1 = mthfc + 1
      cenpi = centfc * pye
      angpi = angfc * pye
      cang = angpi
      ccang = cos(cang)
      acpl = afc * plrad
      rcoil = xmaj + drfc * plrad
      zdtc = twopi/mthfc

      ALLOCATE  ( ellc(nths), ellqc(nths)  )
c      ALLOCATE  ( ellqc(nths), qtgrc(nths), ellc(nths), cgrd(nths) )

c$$$      DO i = 1, mthfc1
c$$$         cgrd(i) = ( 1.0/mthfc ) * (i-1)
c$$$      END DO

c...  Parameterize the containing coil shell. Several options.
c     DEE shaped, sized relative to plrad.
c     Equidistance from the resistive shell.

      IF ( iopfc .EQ. 4 ) THEN

c... Independent of plasma or wall

         DO i = 1, mthfc1
            zth = (i-1)*zdtc
            zsth = sin(zth)
            wxgrd(i) = drfc + afc * cos( zth + dfc*zsth )
            wzgrd(i) = dzfc - bfc * afc * zsth
         END DO
         xcw(1:nths) = wxgrd(1:nths)
         zcw(1:nths) = wzgrd(1:nths)
c     Get the coil tips.
         zcgrd1 = cenpi + angpi
         chgt(1) =  drfc + afc * cos ( zcgrd1 + dfc * sin(zcgrd1) )
         chgt(2) = dzfc - bfc * afc * sin(zcgrd1)
         
         zcgrd1 = cenpi - angpi
         chgt(3) = drfc + afc * cos ( zcgrd1 + dfc * sin(zcgrd1) )
         chgt(4) = dzfc - bfc * afc * sin(zcgrd1)
         
      END IF


      IF ( iopfc .EQ. 5 ) THEN

         DO i = 1, mthfc1
            zth = (i-1)*zdtc
            zsth = sin(zth)
            wxgrd(i) = rcoil + acpl * cos( zth + dfc*zsth )
            wzgrd(i) = - bfc * acpl * zsth
         END DO
         xcw(1:nths) = wxgrd(1:nths)
         zcw(1:nths) = wzgrd(1:nths)
c     Get the coil tips.
         zcgrd1 = cenpi + angpi
         chgt(1) = rcoil + acpl * cos ( zcgrd1 + dfc * sin(zcgrd1) )
         chgt(2) = - bfc * acpl * sin(zcgrd1)
         
         zcgrd1 = cenpi - angpi
         chgt(3) = rcoil + acpl * cos ( zcgrd1 + dfc * sin(zcgrd1) )
         chgt(4) = - bfc * acpl * sin(zcgrd1)
         
      END IF

      IF ( iopfc .EQ. 6 ) THEN   
c     !  Coil at constant distance from wall.
         wcentr = xmaj
         DO iz = 2, mth1
            alph = atan2m ( xwal(iz+1)-xwal(iz-1),
     $           zwal(iz-1)-zwal(iz+1) )
            xcw(iz) = xwal(iz) + afc*plrad * cos(alph)
            zcw(iz) = zwal(iz) + afc*plrad * sin(alph)
         END DO
c     
         xcw(1) = xcw(mth1)
         zcw(1) = zcw(mth1)
         xcw(mth2) = xcw(2)
         zcw(mth2) = zcw(2)
         
c     . .Interpolate to mthfc grid:
         CALL trans ( xcw, mth, wxgrd, mthfc )
         CALL trans ( zcw, mth, wzgrd, mthfc )
c     Coil tips:
         zcgrd1 = ( cenpi+angpi ) / twopi
         IF ( zcgrd1 .le. 0 ) zcgrd1 = 1.0 + zcgrd1
         CALL lagp ( cgrd,wxgrd, mthfc1,4, zcgrd1,chgt(1),dum, 0,1 )
         CALL lagp ( cgrd,wzgrd, mthfc1,4, zcgrd1,chgt(2),dum, 0,1 )
         zcgrd1 = ( cenpi-angpi ) / twopi
         IF ( zcgrd1 .le. 0 ) zcgrd1 = 1.0 + zcgrd1
         CALL lagp ( cgrd,wxgrd, mthfc1,4, zcgrd1,chgt(3),dum, 0,1 )
         CALL lagp ( cgrd,wzgrd, mthfc1,4, zcgrd1,chgt(4),dum, 0,1 )

      END IF

c..   Set up the coil structure

      IF ( icshp .EQ. 1 ) THEN
c..   Normalize to value at theta = 0
         coil0 = 1.0 / ( 1.0/ exp( (1.0-ccang)/taufc ) + 1.0 )
c     
         DO i = 1, mthfc1
            zth = (i-1)*zdtc
c            zthp = zth + cenpi
c ... Change to minus
            zthp = zth - cenpi
            zcth = cos(zthp)
            zexp = exp( ( zcth-ccang)/taufc )
            ccoil(i) = zexp / ( zexp + 1.0 )
            ccoil(i) = ccoil(i) / coil0
         END DO

      END IF

      IF ( icshp .EQ. 2 ) THEN
c..   Normalize to value at theta = 0
         coil0 = 1.0 / ( 1.0/ exp( (1.0-ccang)/taufc ) + 1.0 )
c..   new Normalization
         coil0 = 1.0 / ( 1.0/ exp( cang / taufc ) + 1.0 )
c     
         DO i = 1, mthfc1
            zth = (i-1)*zdtc
            zcth = cos(zth)
            zsth = sin(zth)
c     zexp = exp( ( zcth-ccang)/taufc )
c     try this msc 6/3/01
            zthp = zth 
            if (abs(zthp).le. pye) then
               zexp = exp( ( cang-abs(zthp))/taufc )
            else if (abs(zthp).ge. pye) then
               zexp = exp( ( cang-(2*pye-abs(zthp)))/taufc)
            end if
c     try the above and still use the same normalization
            ccoil(i) = zexp / ( zexp + 1.0 )
            ccoil(i) = ccoil(i) / coil0
c     try this msc 6/3/01
c     if  (zsth .lt. 0.) then
c     ccoil(i)=-ccoil(i)
c     else if ( zsth.eq.0.) then
c     ccoil(i)=0.
c     end if
c     the above added to see if the e matrix goes to zero
c     when the coil angle goes to zero
c     ccoilt(i) = zsth * ( zexp / ( taufc * (zexp+1.0)**2 )
         end do

      END IF

      If ( icshp .EQ. 3 ) THEN

c..   Ming's Modification

c..   Put coil parametrization in wxgrd, wzgrd. Transform to equal arc
c     in xcw, zcw.

c     Replace: taucoil -> taufc. coilang -> angfc. coilcen -> centfc.
c     mtcoil -> mthfc. acoil -> afc. bcoil -> bfc. 
c     dcoil -> dfc. drcoil -> drfc.

c...  Replace statements with newer variables:
c     cang = coilang * pye
c     cangcen = coilcen * pye
         cangcen = cenpi
         ccang = cos(cang)
c     acpl = acoil * plrad
         acpl = afc * plrad
c     rcoil = xmaj + drcoil*(plrad)
         rcoil = xmaj + drfc * plrad
c     zdtc = twopi/mtcoil
         zdtc = twopi/mthfc

c..   Normalize to value at theta = 0
         coil0 = 1.0 / ( 1.0/ exp( (1.0-ccang)/taufc ) + 1.0 )
c..   new Normalization
         coil0 = 1.0 / ( 1.0/ exp( cang / taufc ) + 1.0 )
c     
         do i = 1, mthfc1
            zth = (i-1)*zdtc
            zcth = cos(zth)
            zsth = sin(zth)
c     wxgrd(i) = rcoil + acpl * cos( zth+dcoil*zsth )
c     wzgrd(i) = - bcoil * acpl * zsth
c     zexp = exp( ( zcth-ccang)/taufc )
c     try this msc 6/3/01
            if (cangcen.eq.0.) then
               if (zth.le. pye) then
                  zexp = exp( ( cang-abs(zth))/taufc )
               else if (zth.ge. pye) then
                  zexp = exp( ( cang-(2*pye-abs(zth)))/taufc)
               end if
               ccoil(i) = zexp / ( zexp + 1.0 )
               ccoil(i) = ccoil(i) / coil0
            else if (cangcen.gt.0.) then
               cangbig=cangcen+cang
               cangsml=cangcen-cang
               if (zth.le. cangcen) then
                  zexp = exp( ( zth-cangsml)/taufc )
               else if (zth.ge. cangcen) then
                  zexp = exp( ( cangbig-zth)/taufc)
               end if
               ccoil(i) = zexp / ( zexp + 1.0 )
               ccoil(i) = ccoil(i) / coil0
            else if (cangcen.lt.0.) then
               cangbig=2*pye+cangcen+cang
               cangsml=2*pye+cangcen-cang
               if (zth.le. 2*pye+cangcen) then
                  zexp = exp( ( zth-cangsml)/taufc )
               else if (zth.ge. 2*pye+cangcen) then
                  zexp = exp( ( cangbig-zth)/taufc)
               end if
               ccoil(i) = zexp / ( zexp + 1.0 )
               ccoil(i) = ccoil(i) / coil0
            end if
c     try this msc 6/3/01
c     if  (zsth .lt. 0.) then
c     ccoil(i)=-ccoil(i)
c     else if ( zsth.eq.0.) then
c     ccoil(i)=0.
c     end if 
c     the above added to see if the e matrix goes to zero 
c     when the coil angle goes to zero
            cgrd(i) = ( 1.0/mthfc ) * (i-1)
c     ccoilt(i) = zsth * ( zexp / ( taufc * (zexp+1.0)**2 )
         end do

      END IF

c...  Multiply coils by mu_0 here. Careful that we don't need the current
c       though, since only the fields should be multiplied by mu_0.

c$$$      zmu0 =  4.0 * pye * 1.0e-7
c$$$ 
c$$$      ccoil(1:mthfc1) = zmu0 * ccoil(1:mthfc1)

c$$$c... No equal arc for now.
c$$$      CAll eqarcw ( wxgrd,wzgrd, xcw,zcw,ellqc,wkgrd,zork1, mthfc1 )

c$$$      CALL drawc3(xwal,zwal,xinf,zinf,xcw,zcw,1,mth1,1,mthfc1,
c$$$     $     "z","x",xmaj,xma,zma,
c$$$     $     plrad,1,1,idot, jobid, ishape,a,b,aw,bw,cw,dw,tw,
c$$$     $     abulg, bbulg, tbulg, wcentr, pcircm,wcircm,
c$$$     $     pelong, pdness, wmaj,wlrad, welong, wdness, chgt )

      DEALLOCATE ( ellc,ellqc )

      RETURN
      END SUBROUTINE fbcoil

c............................................................
      SUBROUTINE shfcoilx ( zxin, zzin, ccoil, cgrd, mthfc,
     $     zxout, zzout, ellc, ccout)
c...........................................................

      INCLUDE 'vacuum1.inc'

      DIMENSION zxin(*), zzin(*), ccoil(*), cgrd(*)
      DIMENSION zxout(*), zzout(*), ellc(*), ccout(*)

      REAL, DIMENSION(:), ALLOCATABLE :: ellqc, qtgrc
      REAL, DIMENSION(:), ALLOCATABLE :: wxgrd, wzgrd, wkgrd
      REAL, DIMENSION(:), ALLOCATABLE :: zork1

      mthfc1 = mthfc + 1
      mthfc5 = mthfc + 5

      ALLOCATE (  wxgrd(mthfc5), wzgrd(mthfc5), wkgrd(mthfc5) )
      ALLOCATE  ( ellqc(mthfc5), qtgrc(mthfc5) )
      ALLOCATE  ( zork1(mthfc5) )

      CAll eqarcw ( zxin,zzin, wxgrd,wzgrd, ellqc,wkgrd,zork1, mthfc1 )
      CALL atpoint ( "EQARC at COIL in SHFCOIL","mthfc1", mthfc1,
     $     "ellqc(mthfc1)",ellqc(mthfc1), 6, 23 )
c     
c...  Get equal arcs stuff for the coil function
      call qgrid ( wxgrd, wzgrd, ellqc, qtgrc, mthfc1, wkgrd )
c     
c     Shift needed variables by twopi/2 so that theta=0 is on the inside.
      CALL shfgrb ( ccoil, ccout, mthfc )
c     
c     Get arc length on coilwall, ellc, shifted coordinates...
      zxout(1:mthfc1) = (/ wxgrd(mthfc/2+1:mthfc),wxgrd(1:mthfc/2+1) /)
      zzout(1:mthfc1) = (/ wzgrd(mthfc/2+1:mthfc),wzgrd(1:mthfc/2+1) /)
c     
      ellc(1) = 0.0
      DO ic = 1, mthfc1
         IF ( ic .eq. 1 ) GO TO 9
         thet = ( cgrd(ic)+cgrd(ic-1) ) / 2.0
         CALL lag ( cgrd,zxout,mthfc1,3,thet,f,df,1 )
         xtzt = df
         CALL lag ( cgrd,zzout,mthfc1,3,thet,f,df,1 )
         xtzt = sqrt ( xtzt**2 + df**2 )
         ellc(ic) = ellc(ic-1) + xtzt/(mthfc1-1)
 9       CONTINUE
      END DO
c..   Plot the coil current:
c     
c..   Imaginary part of coil current function is zero:
      zork4 = 0.0

c$$$  c.. Transform coil to equal-arcs grid.
c$$$  c
c$$$  call eqarcw ( zork3,zork4, wxgrd,wzgrd, wkgrd, wkdxt,wkdzt,
c$$$  $        mthfc1 )

      CALL vecwrt ( mthfc, ccout, "CCOIL in shifted coords", 1,mthfc,
     $     outmod, iotty )

      WRITE ( outmod,'( /,
     $     10x,"Theta, phi coil current -- ")')

      DEALLOCATE ( ellqc, qtgrc )
      DEALLOCATE (  wxgrd, wzgrd, wkgrd )
      DEALLOCATE  (  zork1 )

      RETURN
      END SUBROUTINE shfcoilx

