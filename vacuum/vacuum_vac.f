c..................................................................
      subroutine kernel0(xobs,zobs,xsce,zsce,grdgre,gren,j1,j2,isgn,
     .                 iopw,iops,lfele ,nths2,nths)
c......................................................................
c
      dimension grdgre(nths2,*), gren(nths,*)
      dimension xobs(*), zobs(*), xsce(*), zsce(*)
c
      if ( lfele .eq. 3 )  then
         call kernel3(xobs,zobs,xsce,zsce,grdgre,gren,j1,j2,isgn,
     $        iopw,iops)
      else 
         call kernel(xobs,zobs,xsce,zsce,grdgre,gren,j1,j2,isgn,
     .        iopw,iops) 
      end if
c
      return
      end
c
c..................................................................
      subroutine kernel(xobs,zobs,xsce,zsce,grdgre,gren,j1,j2,isgn,
     .                 iopw,iops)
c......................................................................
c
c..... calculates the kernels of the integral equation for laplace's
c     equation for a torus.
c     j1, j2 correspond to the four blocks in grdgre(j1,j2).  they
c     respectively indicate observer and source.
c     1 and 2 are plasma and wall blocks, respectively.
c     iopw = 1 calculates bval terms on rhs.
c     iops = 1 subtract and add log behavior in bval.  log
c              contribution is done analytically.
c
      include 'vacuum1.inc'
      include 'vacuum3.inc'
      include 'vacuum2.inc'
      include 'vacuum5.inc'

      COMMON / grvac / aval0, phex0, zphfc0, iscs, isco,
     $     j1com,j2com

      dimension grdgre(nths2,nths2), gren(nths,nths)
      dimension xobs(*), zobs(*), xsce(*), zsce(*)
      dimension tgaus(8), wgaus(8), tlog(4), wlog(4)
      dimension iop(2), the(nths), xpp(nths), zpp(nths), tab(3)
      dimension work(nths), ww1(nths), ww2(nths), ww3(nths)
      dimension xpr(nths), zpr(nths)
c
      data ischk /0/, isph/0/       ! 08/26/2009, since not SAVEd
c
      call atpoint ( "KERNEL","j1", j1,"xobs(1)", xobs(1),
     $     iotty, outmod )
c      
      factpi = twopi
c      isph = 0
c
      xpp(1) = 0.
      zpp(1) = 0.
      tab(1) = 0.
      ww1(1) = 0.
      ww2(1) = 0.
      ww3(1) = 0.
      the = 0.0                   ! 08/21/2009
c
c AK0I will contain integral of sign*K0 at j = jres.
c.This is only for verifying that the integrals over K_0 are correct.
c ak0i is not used in the algorithm.
c
      write ( outmod,
     $     '(10x,"*************************************************")' )
      if (j1 .eq. 1 ) write ( outmod, '(20x, "j1 = ",i2,
     $     " Observer = Plasma")' ) j1
      if ( j1 .eq. 2 ) write ( outmod, '(20x, "j1 = ",i2,
     $     " Observer = Wall")' ) j1
      if ( j2 .eq. 1 ) write ( outmod, '(20x, "j2 = ",i2,
     $     " Source   = Plasma")' ) j2
      if ( j2 .eq. 2 ) write ( outmod, '(20x, "j2 = ",i2,
     $     " Source   = Wall")' ) j2
      write ( outmod,
     $     '(10x,"*************************************************")' )
c
      ak0i = 0.0
      jres = 1
c
      mthm = mth - 1
c
c.....isgn is positive over wall and negative over plasma...
c
c      isgn = 2*j2 - 3
c
      do 7 j = 1, mth1
      do 7 i = 1, mth1
      gren(i,j) = zero
    7 continue
c
c.....Weights for quadratures.
c
c  Set weights equal to dth for finite elements in GATO
c
      if ( lfele .eq. 1 ) then
         wsimpb1 = dth / 2.0
         wsimpb2 = dth
         wsimpb4 = dth 
         write ( outmod, '(/,20x,"Weights_f = dth")')
         write ( iotty , '(/,20x,"Weights_f = dth")')
c
c$$$         wsimpb1 = dth 
c$$$         wsimpb2 = dth
c$$$         wsimpb4 = dth 
c$$$         write ( outmod, '(/,20x,"Weights_f = dth,dth")')
c$$$         write ( iotty , '(/,20x,"Weights_f = dth,dth")')
c
c$$$         wsimpb1 = dth / three
c$$$         wsimpb2 = two * dth / three
c$$$         wsimpb4 = four * dth / three
c$$$         write ( outmod, '(/,20x,"Weights_f = Simpson")')
c$$$         write ( iotty , '(/,20x,"Weights_f = Simpson")')
      else
c     
c$$$      wsimpb1 = dth/2
c$$$      wsimpb2 = dth
c$$$      wsimpb4 = dth
c$$$      write ( outmod, '(/,20x,"Weights_b = Trapezoidal")')
c$$$      write ( iotty , '(/,20x,"Weights_b = Trapezoidal")')
c
c$$$         wsimpb1 = one
c$$$         wsimpb2 = one
c$$$         wsimpb4 = one
c$$$         write ( outmod, '(/,20x,"Weights_b = Unity")')
c$$$         write ( iotty , '(/,20x,"Weights_b = Unity")')
c
         wsimpb1 = dth / three
         wsimpb2 = two * dth / three
         wsimpb4 = four * dth / three
         write ( outmod, '(/,20x,"Weights_b = Simpson")')
         write ( iotty , '(/,20x,"Weights_b = Simpson")')
              end if
c     
c$$$      wsimpa1 = dth/2
c$$$      wsimpa2 = dth
c$$$      wsimpa4 = dth
c$$$      write ( outmod, '(/,20x,"Weights_a = Trapezoidal")')
c$$$      write ( iotty , '(/,20x,"Weights_a = Trapezoidal")')
c
c$$$         wsimpa1 = one
c$$$         wsimpa2 = one
c$$$         wsimpa4 = one
c$$$         write ( outmod, '(/,20x,"Weights_a = unity")')
c$$$         write ( iotty , '(/,20x,"Weights_a = unity")')
c
         wsimpa1 = dth / three
         wsimpa2 = two * dth / three
         wsimpa4 = four * dth / three
         write ( outmod, '(/,20x,"Weights_a = Simpson")')
         write ( iotty , '(/,20x,"Weights_a = Simpson")')
c
c.....weights for eight point gaussian quadrature.
c
      wgaus(1) = 0.101228536290376
      wgaus(2) = 0.222381034453374
      wgaus(3) = 0.313706645877887
      Wgaus(4) = 0.362683783378362
      wgaus(5) = wgaus(4)
      wgaus(6) = wgaus(3)
      wgaus(7) = wgaus(2)
      wgaus(8) = wgaus(1)
c
c.....calculate some quantities for the singular region
c     for the analytic integral over the log behavior of bval.
c
      third = one/three
      algdth = alog(dth)
      slog1m = third * dth * ( algdth - third )
      slog0 = four * third * dth * ( algdth - four*third )
      slog1p = slog1m
      tdth = two * dth
      alg = alog(tdth)
      alg0 = 16.0*dth*(alg - 68.0/15.0) / 15.0
      alg1 = 128.0*dth*(alg - 8.0/15.0) / 45.0
      alg2 = 4.0*dth*(7.0*alg - 11.0/15.0) / 45.0
c
c.....get x coordinates for wall to check for spher case. put them in
c     temporary vector ww1.  then get jtop and jbot.  these
c      correspond to the grid points just to the right
c     of the major axis.

c..!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c   UP-DOWN SYMMETRY ASSUMPTION HERE
C..!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      DO i = 1, mth1
         the(i) = (i-1) * dth
      END DO

      if ( ischk .gt. 0 ) go to 9
      jbot = mth/2 + 1
      jtop = mth/2 + 1
      isph = 0

      IF ( farwal ) GO TO 1003   ! Don't check for spheromak.

      WRITE ( outmod, '("**** Checking  spher case at ischk. " )')
c.. Instead of calling SUB wwall again just use xwal, zwal:
      ww1(1:mth1) = xwal(1:mth1)
      ww2(1:mth1) = zwal(1:mth1)
c$$$      call wwall ( mth1, ww1, ww2 )
c$$$      if ( check2 )
c$$$          write(outmod,2000) (i,xinf(i),ww1(i),zinf(i),ww2(i),i=1,mth1)
c$$$ 2000 format (//,5x,"ith, xpla, xwal, zpla, zwal = ",
c$$$     .          /, 200(2x,i4, 4e12.4,/) )
c
      do 10 i = 1, mth1
c$$$      the(i) = (i-1) * dth    ! 08/21/2009
c
c......take care of isph case when points on wall are near to the
c     centre line.  that is, lagrange interpolation points are always to
c     right of the center line.
c
      if ( i .eq. mth1 ) go to 10
      if ( ww1(i)*ww1(i+1) .gt. zero ) go to 10
      if ( ww1(i) .gt. zero ) jbot = i
      if ( ww1(i) .lt. zero ) jtop = i+1
      isph = 1
   10 continue
c
    9 continue
      ischk = 1
c
      if ( check2 ) then
      write ( 23,1001 ) j1, j2, jtop, jbot
 1001 format ( /,1x, "j1, j2, jtop, jbot = ", 4i4,/ )
      write ( 23,1002 ) ww1(jtop), ww1(jbot), ww2(jtop), ww2(jbot)
 1002 format ( 1x,"xtop, xbot = ",2e11.4,2x,"ztop, zbot = ",2e11.4,/)
      end if

 1003 CONTINUE                   ! Spheromak check end
c
c     periodic boundary conditions.
c
      iop(1)  = 4
      iop(2)  = 4
      call spl1d1( mth1, the, xsce, xpp, iop, 1, ww1, ww2, ww3 )
c     ( note.. ww1, ww2, ww3 just used as temporary working storage here )
      call spl1d1( mth1, the, zsce, zpp, iop, 1, ww1, ww2, ww3 )
c
      do 11 i = 1, mth1
c
      theta = (i-1)*dth
      call spl1d2( mth1, the, xsce, xpp, 1, theta, tab )
      xpr(i) = tab(2)
      call spl1d2( mth1, the, zsce, zpp, 1, theta, tab )
      zpr(i) = tab(2)
c
   11 continue
c
c.....index j runs over the observer points.
c
      do 200 j = 1, mth
c
      xs = xobs(j)
      zs = zobs(j)
      thes = the(j)
c.....calculate non-singular part first.
c
c.....  array work will contain grdgre kernel.
c
      do 12 iw = 1, mth1
      work(iw) = zero
   12 continue
c
c.....there is no observer for xobs negative.  put zero in these
c     positions in grdgre.
c
      if ( xs .lt. zero ) go to 175
c
      aval1 = zero
c
c$$$
c$$$ In the next do loop, ic indicates the theta values of the
c$$$ source points in the non singular region.   i is the just the
c$$$ integration index.  ic = j would be the singular point. mths is 
c$$$ the number of integration points.

c$$$
c$$$  iend ensures that the nonsingular region starts two grid points from the
c$$$  singular point. Thus iend is 2 for the plasma (j1 = 1) and for a wall
c$$$  which does not intersect the major axis. If the wall does, like the
c$$$  spherical (isph = 1 ) wall then iend is calculated for each case.
c
      iend = 2
      if ( .not. (isph .eq. 1 .and. j2 .eq. 2) ) go to 15
      if ( jbot-j .eq. 1 ) iend = 3
      if ( jbot-j .eq. 0 ) iend = 4
      if ( j-jtop .eq. 0 ) iend = 0
      if ( j-jtop .eq. 1 ) iend = 1
   15 continue
      istart = 4 - iend
c
      mths = mth - (istart+iend-1)
c
      do 25 i = 1, mths
c
      ic = j+i+istart-1
      if ( ic .ge. mth1 ) ic = ic - mth
c     ic = mod ( j+i-1, mth ) +1
      theta = (ic-1) * dth
c
      xt = xsce(ic)
      zt = zsce(ic)
c
c.....there are no sources from negative x. for isph=1, the virtual
c     sources are taken care of through iend above.
c
      if ( xt .lt. zero ) go to 25
c
      xtp  = xpr(ic)
      ztp  = zpr(ic)
c
c.....ic should never be equal to j unless ispher=1, in which
c     case we are on the major axis. negligible contribution there.
c
      if ( ic .eq. j ) go to 25
c
      call green
c
c.....weights for simpson integration of grdgre and grwp.
c..... note that mths must be odd.
c
c
      wsimpb = wsimpb2
      if ( (i/2)*2 .eq. i ) wsimpb = wsimpb4
      if ( i .eq. 1 .or. i .eq. mths ) wsimpb = wsimpb1
c
      wsimpa = wsimpa2
      if ( (i/2)*2 .eq. i ) wsimpa = wsimpa4
      if ( i .eq. 1 .or. i .eq. mths ) wsimpa = wsimpa1
c
      work(ic) = work(ic) +  isgn * aval * wsimpa
      gren(j,ic) = gren(j,ic) + bval * wsimpb
c
ccc      work(ic) = isgn * aval * wsimpa
ccc      gren(j,ic) = bval * wsimpb
c
      aval1 = aval1 + aval0 * wsimpa
c
   25 continue
c
c.....end of non-singular part.
c
c.....dont compute singular contributions in off-diagonal blocks
c     when xt less than zero.
c
      j1j2 = j1 + j2
      if (j1j2 .eq. 2 ) go to 28
      if ( (isph.eq.1) .and. (j.gt.jbot) .and. (j.lt.jtop) ) go to 175
   28 continue
c
c.....An eight point  gaussian quadrature is used.
c
      thes = the(j)
      js1 = mod ( j-iend+mth-1, mth ) + 1
      js2 = mod ( j-iend+mth, mth ) + 1
      js3 = mod ( j-iend+mth+1, mth ) + 1
      js4 = mod ( j-iend+mth+2, mth ) + 1
      js5 = mod ( j-iend+mth+3, mth ) + 1
c
      do 165 ilr = 1, 2
c
c.....  xl, xu is upper and lower limits of the gaussian integration
c.....  in the singular region. its done in two parts.
c
      xl = thes + (2*ilr - iend - 2)*dth
      xu = xl + tdth
c
      agaus = half*(xu+xl)
      bgaus = half*(xu-xl)
c
      c1gaus = 0.960289856497536 * bgaus
      c2gaus = 0.796666477413627 * bgaus
      c3gaus = 0.525532409916329 * bgaus
      c4gaus = 0.183434642495650 * bgaus
c
      tgaus(1) = agaus - c1gaus
      tgaus(2) = agaus - c2gaus
      tgaus(3) = agaus - c3gaus
      tgaus(4) = agaus - c4gaus
      tgaus(5) = agaus + c4gaus
      tgaus(6) = agaus + c3gaus
      tgaus(7) = agaus + c2gaus
      tgaus(8) = agaus + c1gaus
c
      do 160 ig = 1, 8
c
      tgaus0 = tgaus(ig)
      if ( tgaus0 .lt. zero ) tgaus0 = twopi + tgaus0
      if ( tgaus0 .ge. twopi ) tgaus0 = tgaus0 - twopi
      call spl1d2 ( mth1, the, xsce, xpp, 1, tgaus0, tab )
      xt = tab(1)
      xtp = tab(2)
      call spl1d2 ( mth1, the, zsce, zpp, 1, tgaus0, tab )
      zt = tab(1)
      ztp = tab(2)
c
      call green
c
c..... subtract out log behavior if iops is one.
c
c... Don't use modulo 2pi for tgaus here!
      bval = bval + iops * alog( ( thes-tgaus(ig) )**2) / xs
c
c.....chi is approximated in the wall plasma region by
c     a five point lagrange interpolation formula.
c
      pgaus = ( tgaus(ig) - thes-(2-iend)*dth ) / dth
      pgaus2 = pgaus*pgaus
      wgbg = wgaus(ig) * bgaus
c
      amm = (pgaus2-one) * pgaus * ( pgaus - two ) / 24.0
      amm = amm * wgbg
      work(js1) = work(js1) + isgn*aval*amm
c
      am = -(pgaus-one) * pgaus * (pgaus2-four) / 6.0
      am = am * wgbg
      work(js2) = work(js2) + isgn*aval*am
c
      a0 = (pgaus2-one) * (pgaus2-four) / four
      work(js3) = work(js3) + isgn*a0*aval*wgbg
c
      ap = -(pgaus+one) * pgaus * (pgaus2-four) / 6.0
      ap = ap * wgbg
      work(js4) = work(js4) + isgn*aval*ap
c
      app = (pgaus2-one) * pgaus * ( pgaus+two) / 24.0
      app = app * wgbg
      work(js5) = work(js5) + isgn*aval*app
c
c.....add in integral over k0 in singular region.
c
c... comment out this statement if NOT adding K_0 contribution....***1 of 3***
      work(j) = work(j) - isgn*aval0*wgbg
c
      if ( j .eq. jres ) ak0i = ak0i - isgn*aval0*wgbg
c
c..... dont need gren if iopw is zero.
c
      if ( iopw .eq. 0 ) go to 160
c
      gren(j,js1) = gren(j,js1) + bval*amm
      gren(j,js2)  = gren(j,js2)  + bval*am
      gren(j,js3)   = gren(j,js3)   + bval*a0*wgbg
      gren(j,js4)  = gren(j,js4)  + bval*ap
      gren(j,js5) = gren(j,js5) + bval*app
c
  160 continue
  165 continue
c
c..... add in integral over k0 in nonsingular region.
c

c... Do only the K_0 residues here. The del**2 G will be taken care 
c    of later.    06/16/00 -> 6/21/2011
      
c... For shell not enclosing the plasma..     !! June/21/2011
      residu = 0.0
      if ( j1 .eq. j2 ) residu = 1.0
c$$$  if ( j1 .eq. j2 ) residu = 2.0

c Add in residues:
c
c$$$      if ( ishape .lt. 100 ) then           !! June 21/2011
      if ( ishape .lt. 200 ) then               !! June 21/2011
c.. Shell enclosing the plasma - usual cae      !! June 21/2011
c. Residue from del**2 G minus half contours::
         resdg = (2-j1)*(2-j2) + (j1-1)*(j2-1)
c. Residue from K^0 integrals:
         resk0 = (2-j1)*(2-j2) + (j1-3)*(j2-1)
c... If need to take out K_0 contribution: Set resk0 = 0.0....***2 of 3***
c         resk0 = 0.0
c$$$         residu = resdg + resk0
         residu = resk0
c  Note: residu can be written as 2*(j1-2)*(2*j2-3) = 2*isgn*(j1-2) as above.
      end if
c
c...If need to take out k_0 contribution: Set aval1 = 0.0....***3 of 3***
c      aval1 = 0.0
      work(j) = work(j) - isgn*aval1 + residu
c
      if ( j .eq. jres ) ak0i = ak0i - isgn*aval1
c
      if ( iops .ne. 1 .or. iopw .eq. 0 ) go to 170
c
c.......add in analytical integrals over log behavior.
c
      gren(j,js1) = gren(j,js1) - alg2/xs
      gren(j,js2)  = gren(j,js2) - alg1/xs
      gren(j,js3)   = gren(j,js3) - alg0/xs
      gren(j,js4)  = gren(j,js4) - alg1/xs
      gren(j,js5) = gren(j,js5) - alg2/xs
c
  170 continue
c
  175 continue
c
      if ( (xs .lt. zero) .and. (j2 .eq. 2) ) work(j) = 1.0
c
c$$$      if ( .not. checkd ) go to 176
c$$$c
c$$$      if( (j.eq.1) .or. (j.eq.jbot-1) .or. (j.eq.jbot)
c$$$     .  .or. (j.eq.50) .or. (j.eq.mth) )
c$$$     .write(23,300) j, (work(iwr),iwr=1,mth1)
c$$$  300 format(1x,"j = ", i4,/, 20(1x,10e11.4,/) )
c
  176 continue
c
c.....now put work vector in appropriate place in grdgre .....
c
c                   6/27/94.....
c       Divide gren by twopi to be consistent with definition of grdgre. 
c
      do 180 ic = 1, mth
c
c     j1j2 = mth * ( 2*(j2-1)*mth + j1-1 + 2*(ic-1) ) + j
c     grdgre(j1j2) = work(ic)
      grdgre ( (j1-1)*mth+j, (j2-1)*mth+ic ) = work(ic)
c
      gren(j,ic) = gren(j,ic) / factpi
c
  180 continue
c$$$            write ( outmod, '("work: j =" i3, / (1p10e12.4) )' )
c$$$     $           j, (work(iw), iw=1, 20)
c$$$            write ( outmod, '("gren: j =" i3, / (1p10e12.4) )' )
c$$$     $           j, (gren(j,iw), iw=1, 20)
c
  200 continue
c End of observer do-loop

c      Add in del**2 residue here  ! 06/16/00
 
         IF ( j1 == j2 ) THEN
            indx0 = (j1-1) * mth 
            DO iker = 1, mth
               grdgre(indx0+iker,indx0+iker) =
     $              grdgre(indx0+iker,indx0+iker) + 1.0
            END DO
         END IF

       write ( 23,330 ) jres, j1,j2, ak0i,
     $     resdg, resk0, residu, j1,j2
       write ( iotty,330 ) jres, j1,j2, ak0i,
     $     resdg, resk0, residu, j1,j2
  330 format ( /,1x,"jres, j1,j2, isgn*AK0I,  = ", 3i3,1pe14.5,/
     $     " Resdg, resk0, Residue added = ",
     $     1p3e11.2, " in ",2i2, "  block", / )
c
      return
      end
c     
c..................................................................
      subroutine kernel3(xobs,zobs,xsce,zsce,grdgre,gren,j1,j2,isgn,
     .     iopw,iops)
c......................................................................
c     
c.....calculates the kernels of the integral equation for laplace's
c     equation for a torus for finite elements.
c     j1, j2 correspond to the four blocks in grdgre(j1,j2).  they
c     respectively indicate observer and source.
c     1 and 2 are plasma and wall blocks, respectively.
c     iopw = 1 calculates bval terms on rhs.
c     iops = 1 subtract and add log behavior in bval.  log
c     contribution is done analytically.
c     
c     The integral over the source is local on each element but the 
c     integral of the singularity treatment of K is done over 2pi
c     to compensate for the full residue. 
  
      include 'vacuum1.inc'
      include 'vacuum3.inc'
      include 'vacuum2.inc'
      include 'vacuum5.inc'

      COMMON / grvac / aval0, phex0, zphfc0, iscs, isco,
     $     j1com,j2com

      dimension grdgre(nths2,nths2), gren(nths,nths)
      dimension xobs(*), zobs(*), xsce(*), zsce(*)
      dimension wgaus3(3), cgaus3(3), algaus(10)
      dimension wgaus5(5), cgaus5(5)
      dimension wgausn(5), cgausn(5)
      dimension iop(2), the(nths), xpp(nths), zpp(nths), tab(3)
      dimension xppo(nths), zppo(nths)
      dimension work(nths), worg(nths), ww1(nths), ww2(nths), ww3(nths)
      dimension xpr(nths), zpr(nths)
c     
      data ischk /0/
      j1com = j1
      j2com = j2
c     
      call atpoint ( "KERNEL3","j1", j1,"xobs(1)", xobs(1),
     $     iotty, outmod )
c     
c$$$  write ( outmod, '("iops, iopw, isgn = ", 3i4 )')
c$$$  $     iops, iopw, isgn
c     
      factpi = twopi
c     
c...  Choice of either 3 or 5 point Gaussian for observer point in 
c     each element. Choose ngaus = 3, or 5 here. Put in input maybe.
c     Doesn't seem to make any difference!
c     
c...  ngaus is the first digit of lgaus:
      ngaus = lgaus / 10
      ngm1 = (ngaus-1)/2 + 1              ! Gaussian mid point 
c     
      tdth = two * dth

      xpp(1) = 0.
      zpp(1) = 0.
      tab(1) = 0.
      ww1(1) = 0.
      ww2(1) = 0.
      ww3(1) = 0.
c     
c     AK0I will contain integral of sign*K0 at j = jres, jres2
c     .This is only for verifying that the integrals over K_0 are correct.
c     ak0i is not used in the algorithm.
c     
      write ( outmod,
     $     '(10x,"*************************************************")' )
      if (j1 .eq. 1 ) write ( outmod, '(20x, "j1 = ",i2,
     $     " Observer = Plasma")' ) j1
      if ( j1 .eq. 2 ) write ( outmod, '(20x, "j1 = ",i2,
     $     " Observer = Wall")' ) j1
      if ( j2 .eq. 1 ) write ( outmod, '(20x, "j2 = ",i2,
     $     " Source   = Plasma")' ) j2
      if ( j2 .eq. 2 ) write ( outmod, '(20x, "j2 = ",i2,
     $     " Source   = Wall")' ) j2
      write ( outmod,
     $     '(10x,"*************************************************")' )
c     
      ak0i = 0.0
      jres = 1
      jres2 = 5
c     
      mthm = mth - 1
c     
c.....isgn is positive over wall and negative over plasma...
c     
c     isgn = 2*j2 - 3
c     
      DO i = 1, mth
         DO j = 1, mth
            indx1 = (j1-1) * mth + i
            indx2 = (j2-1) * mth + j
            IF ( lhighn == 2 ) THEN
               indx1 = ( 2*(j1-1)+(isco-1) ) * mfel + i
               indx2 = ( 2*(j2-1)+(iscs-1) ) * mfel + j
            END IF
            gren(i,j) = zero
            grdgre ( indx1, indx2 ) = 0.0
         END DO
      END DO

c.....Weights for quadratures.
c     
c.....weights for three point gaussian quadrature. Used for observer points
c     
         wgaus3(1) = 0.555555555555556 / 2.0
         wgaus3(2) = 0.888888888888889 / 2.0
         wgaus3(3) = wgaus3(1)
c     
         cgaus3(1) = 0.774596669241483 
         cgaus3(2) = 0.0 
         cgaus3(3) = - cgaus3(1)
c     
c.....weights for five point gaussian quadrature.
c     
         wgaus5(1) = 0.236926885056189 / 2.0
         wgaus5(2) = 0.478628670499366 / 2.0
         wgaus5(3) = 0.568888888888889 / 2.0
         wgaus5(4) = wgaus5(2)
         wgaus5(5) = wgaus5(1)
c     
         cgaus5(1) = 0.906179845938664
         cgaus5(2) = 0.538469310105683
         cgaus5(3) = 0.0
         cgaus5(4) = - cgaus5(2)
         cgaus5(5) = - cgaus5(1)
c     
c.... Put Gaussian weights etc. in working arrays:
c     
         do i = 1, ngaus
            if ( ngaus .eq. 3 ) then
               wgausn(i) = wgaus3(i)
               cgausn(i) = cgaus3(i)
            end if
            if ( ngaus .eq. 5 ) then
               wgausn(i) = wgaus5(i)
               cgausn(i) = cgaus5(i)
            end if
         end do
c     
c.....calculate some quantities for the singular region
c     for the analytic integral over the log behavior of bval.
c     
c...........For finite elements:
c     
         hdel = pye / mfel                              ! half grid
         thdel = 2.0 * hdel
         alghdl = alog ( hdel )
         flgii = thdel * ( alghdl - 1.0 )               ! log integral over
                                                        ! center grid
         flgmp = thdel * ( alghdl + 0.6479184335 )      ! adjacent grids
c     
         write ( outmod, '(/,"flgii = ", 1pe12.4 )' ) flgii
c     
c...  Log integrals for arbitrary x0, x1,x2: int_x1^x2 ln(x-x0) dx:
c      x0(i) are Gaussian observer points. Will use ngaus of them in each 
c      obserever finite element.
c....  Change variable st. x0 = 0.0, x1 = -(x0-x1), x2 = +(x2-x0):
c     
         do i = 1, ngaus
            x0 = 0.0
            x1 = - ( 1.0 + cgausn(i) ) * hdel
            x2 =   ( 1.0 - cgausn(i) ) * hdel
            call alg012 ( x0,x1,x2, algaus(i) )
         end do
c     
         write ( outmod, '(/,"algaus = ", 1p3e12.4 )' )
     $        (algaus(i), i=1,ngaus )
c     
c.....get x coordinates for wall to check for spher case. put them in
c     temporary vector ww1.  then get jtop and jbot.  these
c     correspond to the grid points just to the right
c     of the major axis.
c     
         if ( ischk .gt. 0 ) go to 9
         jbot = mth/2 + 1
         jtop = mth/2 + 1
         isph = 0

      WRITE ( outmod, '("**** Checking  spher case at ischk. " )')
c.. Instead of calling SUB wwall again just use xwal, zwal:
      ww1(1:mth1) = xwal(1:mth1)
      ww2(1:mth1) = zwal(1:mth1)

c$$$         call wwall ( mth1, ww1, ww2 )
c$$$  if ( check2 )
c$$$  $     write(outmod,2000) (i,xinf(i),ww1(i),zinf(i),ww2(i),i=1,mth1)
c$$$  2000 format (//,5x,"ith, xpla, xwal, zpla, zwal = ",
c$$$  .          /, 200(2x,i4, 4e12.4,/) )
c     
         do 10 i = 1, mth1
            the(i) = (i-1) * dth
c     
c......take care of isph case when points on wall are near to the
c     centre line.  that is, lagrange interpolation points are always to
c     right of the center line.
c     
            if ( i .eq. mth1 ) go to 10
            if ( ww1(i)*ww1(i+1) .gt. zero ) go to 10
            if ( ww1(i) .gt. zero ) jbot = i
            if ( ww1(i) .lt. zero ) jtop = i+1
            isph = 1
 10      continue
c     
 9       continue
         ischk = 1
c     
         if ( check2 ) then
            write ( 23,1001 ) j1, j2, jtop, jbot
 1001       format ( /,1x, "j1, j2, jtop, jbot = ", 4i4,/ )
            write ( 23,1002 ) ww1(jtop), ww1(jbot), ww2(jtop),
     $           ww2(jbot)
 1002       format ( 1x,"xtop, xbot = ",2e11.4,2x,
     $           "ztop, zbot = ",2e11.4,/)
         end if
c     
c     periodic boundary conditions.
c     
         iop(1)  = 4
         iop(2)  = 4
         call spl1d1( mth1, the, xsce, xpp, iop, 1, ww1, ww2, ww3 )
c     ( note.. ww1, ww2, ww3 just used as temporary working storage here )
         call spl1d1( mth1, the, zsce, zpp, iop, 1, ww1, ww2, ww3 )
c     
c...  Need to interopolate on observer points because of Gauss.
         call spl1d1( mth1, the, xobs, xppo, iop, 1, ww1, ww2, ww3 )
c     ( note.. ww1, ww2, ww3 just used as temporary working storage here )
         call spl1d1( mth1, the, zobs, zppo, iop, 1, ww1, ww2, ww3 )
c     
         do 11 i = 1, mth1
c     
            theta = (i-1)*dth
            call spl1d2( mth1, the, xsce, xpp, 1, theta, tab )
            xpr(i) = tab(2)
            call spl1d2( mth1, the, zsce, zpp, 1, theta, tab )
            zpr(i) = tab(2)
c     
 11      continue

c$$$
c$$$      WRITE ( OUTMOD, '("qth0 = ",/, (10es12.4))' )
c$$$     $     (qth0(i), i = 1, mth2)  
c$$$      WRITE ( OUTMOD, '("qth0pp = ",/, (10es12.4))' )
c$$$     $     (qth0pp(i), i = 1, mth2)  

c.....index j runs over the observer points.

         do 200 j = 1, mth

c.... j3 do loop is for gaussian points on the observer.

            do 180 j3 = 1, ngaus
c$$$  c
c$$$  write ( outmod, '(/,20x,"DOING J3 = ", i3,/)' ) j3
c$$$  write ( iotty,  '(/,20x,"DOING J3 = ", i3,/)' ) j3
c     
c.....The array work will contain rows for the grdgre kernel.
c               wrog for the gren kernel
c    These will be filled for each Gaussian observation index j3.
c     
               do 12 iw = 1, mth1
                  work(iw) = zero
                  worg(iw) = zero
 12            continue
c     
c... Get quantities for the Gaussian points about observer at J
c    phex0 is the secular phase exponent calculated at the 
c    observer point (including the Gaussian ones) for the high-n option. 
c    This will be passed to the kernel integration.

               xs = xobs(j)
               zs = zobs(j)
               thes = the(j)
               thesg0 = thes
               qthj3 = qth0(j)
               phex0 = 0.0
               phex00 = qa1 * ( qthj3 + thes )
               IF ( (lhighn == 2) .AND. (j1com == 1) )
     $              phex0 = phex00
               IF ( j3 .NE. ngm1 ) THEN
                  thes = the(j) + cgausn(j3) * hdel
                  thesg0 = thes
                  IF ( thesg0 .lt. zero ) thesg0 = twopi + thesg0
                  IF ( thesg0 .ge. twopi ) thesg0 = thesg0 - twopi
                  call spl1d2( mth1, the, xobs, xppo, 1, thesg0, tab )
                  xs = tab(1)
                  call spl1d2( mth1, the, zobs, zppo, 1, thesg0, tab )
                  zs = tab(1)
                  phex0 = 0.0
                  IF ( lhighn == 2 ) THEN
                     call spl1d2( mth1,the, qth0,qth0pp, 1,thesg0,tab )
                     qthj3 = tab(1)
                     phex00 = qa1 * ( qthj3 + thes )
                  END IF
                  IF ( (lhighn == 2) .AND. (j1com == 1) ) THEN
                     phex0 = phex00
                  END IF
               END IF

c...  Now setup source highn phase for singularity treatment
c     The phase is evaluated at the observer theta.
c     Phext0 is the phase due to the SOURCE but at the 
c     singular angle, i.e., the obserever theta.

               phext0 = 0.0
               zphfc0 = 1.0
               IF ( lhighn == 2 ) THEN
                  IF (j2com == 1) phext0 = phex00
                  dphex0 = phext0 - phex0
c                  IF ( (j1com == 2) .AND. (j2com == 2) ) dphex0 = 0.0
                  IF ( iscs == 1 ) zphfc0 = cos (n*dphex0)
                  IF ( iscs == 2 ) zphfc0 = sin (n*dphex0)
               END IF

c$$$               write ( outmod, '("j, j3, qthj3, thes, phex0 = ",
c$$$     $              2i4,1x, 3es12.4)' ) j, j3,qthj3, thes, phex0 

c.....There is no observer for xobs negative.  Leave zero in these
c     positions in grdgre.
c     
               if ( xs .lt. zero ) go to 175
c     
               aval1 = zero
               ak00 = zero
               ak0i = zero
c     
c$$$  
c$$$  In the next do loop, ic indicates the theta values of the
c$$$  source points in the non singular region.   i is the just the
c$$$  integration index.  ic = j would be the singular point. mths is 
c$$$  the number of integration points.

c$$$  
c$$$  iend ensures that the nonsingular region starts two grid points from the
c$$$  singular point. Thus iend is 2 for the plasma (j1 = 1) and for a wall
c$$$  which does not intersect the major axis. If the wall does, like the
c$$$  spherical (isph = 1 ) wall then iend is calculated for each case.
c     
               iend = 2
               if ( .not. (isph .eq. 1 .and. j2 .eq. 2) ) go to 15
               if ( jbot-j .eq. 1 ) iend = 3
               if ( jbot-j .eq. 0 ) iend = 4
               if ( j-jtop .eq. 0 ) iend = 0
               if ( j-jtop .eq. 1 ) iend = 1
 15            continue
c     
               mths = mth
c     
               do 25 i = 1, mths
c     
                  ic = j + i - 1
                  if ( ic .ge. mth1 ) ic = ic - mth
                  theta = (ic-1) * dth
c     
c.....xl, xu is lower and upper limits of the gaussian integration
c     
                  xl = theta - dth/2.0
                  xu = xl + dth
                  ising = 0
                  if ( ic .eq. j ) ising = 1
c     
c..   Choose gausi4 or gausi6 here for 4 of 5 points Gaussian integration 
c     for source points. Not much difference!
c     
c..   Switch lg46 is the second digit of lgaus:
                  lg46 = lgaus - (lgaus/10)*10
                  if (lg46 .eq. 4 )
     $                 call gausi4 ( the, thes, xl,xu, xsce,zsce,
     $                 xpp,zpp, isgn,iopw,iops,ising,aker,bker,aker0 )
                  if (lg46 .eq. 6 )
     $                 call gausi6 ( the, thes, xl,xu, xsce,zsce,
     $                 xpp,zpp, isgn,iopw,iops,ising,aker,bker,aker0 )
c     
                  work(ic) = isgn*aker   ! localized on finite element
                  worg(ic) = bker        ! localized on finite element
c...  aval1 contains unsigned integral of K^0.
                  aval1 = aval1 + aker0  ! summed over whole domain
                  if ( (j .eq. jres) .or. (j .eq. jres2) ) then
                     ak00 = ak00 - ising*isgn*aker0
                     ak0i = ak0i - (1-ising)*isgn*aker0
                  end if
c     
c$$$  write (outmod, '(
c$$$  $"jres,j,i,ic,ising,xl,theta,xu,aker,bker,aker0, aval1,ak0i=",
c$$$  $        5i4,/, 1p8e12.4)' ) jres,j,i,ic,ising,xl,theta,xu,
c$$$  $              aker,bker,aker0, aval1,ak0i
 25            continue
c     
c.....don't compute singular contributions in off-diagonal blocks
c     when xt less than zero.
c     
               j1j2 = j1 + j2
               if (j1j2 .eq. 2 ) go to 28
               if ( (isph.eq.1) .and. (j.gt.jbot) .and. (j.lt.jtop) )
     $              go to 175
 28            continue
c     
c.....singular region is done only if iops = 1. 
c     
c$$$  js1 = mod ( j-iend+mth-1, mth ) + 1
c$$$  js2 = mod ( j-iend+mth, mth ) + 1
c$$$  js3 = mod ( j-iend+mth+1, mth ) + 1
c$$$  js4 = mod ( j-iend+mth+2, mth ) + 1
c$$$  js5 = mod ( j-iend+mth+3, mth ) + 1
c     
c.....Subtract error: Integral k0 - resk^0

c... Do only the K_0 residues here. The del**2 G will be taken care 
c    of in calling routine.    06/16/00
c     
               residu = 0.0
c$$$         if ( j1 .eq. j2 ) residu = 2.0
               if ( j1 .eq. j2 ) residu = 1.0

c     Add in residues:

               if ( ishape .lt. 100 ) then
c     . Residue from del**2 G minus half contours::
                  resdg = (2-j1)*(2-j2) + (j1-1)*(j2-1)
c     . Residue from K^0 integrals:
                  resk0 = (2-j1)*(2-j2) + (j1-3)*(j2-1)
c$$$                  residu = resdg + ( zphfc0 * resk0 )
                  residu = resk0 
c     Note: residu can be written as 
c           2*(j1-2)*(2*j2-3) = 2*isgn*(j1-2) as above.
               end if

c.. Now, compensate for the singularities at the source integrals. 
c    At the index corresponding to the observer at index J.

               work(j) = work(j) - zphfc0 * ( isgn*aval1 - residu )
c     
               if ( iops .ne. 1 .or. iopw .eq. 0 ) go to 170
c     
c.......add in analytical integrals over log behavior, modified by high-n
               
c     worg(j-1)   = worg(j) - 2.0*flgmp/xs
               worg(j)     = worg(j) - 2.0*zphfc0*algaus(j3)/xs
c     worg(j+1)   = worg(j) - 2.0*flgmp/xs

 170           continue
               if ( (j .eq. jres) .or. (j .eq. jres2) ) then
                  write ( iotty,  '("J, J3 = ", 2i3)' ) j, j3
                  write ( outmod, '("J, J3 = ", 2i3)' ) j, j3
                  write(iotty,'(1x, "ak00, ak0i, aval1 = ",
     $                 1p3e12.4)') ak00, ak0i, aval1
                  write(outmod,'(1x, "ak00, ak0i, aval1 = ",
     $                 1p3e12.4)') ak00, ak0i, aval1
               end if
               if ( (j .eq. jres) .or. (j .eq. jres2) ) then 
                  write(iotty,'(1x,
     $                 "ISCS= ", i3, " PHEXT0, PHEX0, ZPHFC0 = ",
     $                 1p3e12.4)') iscs, phext0, phex0, zphfc0
                  write(outmod,'(1x,
     $                 "ISCS= ", i3, " PHEXT0, PHEX0, ZPHFC0 = ",
     $                 1p3e12.4)') iscs, phext0, phex0, zphfc0
               end if

 175           continue
c     
               if ( (xs .lt. zero) .and. (j2 .eq. 2) ) work(j) = 1.0
c     
c$$$  if ( .not. checkd ) go to 176
c$$$  c
c$$$  if( (j.eq.1) .or. (j.eq.jbot-1) .or. (j.eq.jbot)
c$$$  .  .or. (j.eq.50) .or. (j.eq.mth) )
c$$$  .write(23,300) j, (work(iwr),iwr=1,mth1)
c$$$  300 format(1x,"j = ", i4,/, 20(1x,10e11.4,/) )
c     
 176           continue
c     
c.....now put work vector in appropriate place in grdgre .....

c     6/27/94.....
c     Divide gren by twopi to be consistent with definition of grdgre. 
c     
               do ic = 1, mth
c     
c     The next grdgre and gren sum will include the Gaussian observer points 
c        of the do index j3. Put in high-n phase factor for the 
c        observer point here. New indices
c        of GRDGRE(indx1,indx2) to incorporate real and inaginary parts
c        coming from observer and source. ISCO, ISCS = 1 for old style.
c        These are passed from VACCAL through the common block GRVAC
c        Indices of GREN is taken care of in calls

c     j1j2 = mth * ( 2*(j2-1)*mth + j1-1 + 2*(ic-1) ) + j
c     grdgre(j1j2) = work(ic)
c                  grdgre ( (j1-1)*mth+j, (j2-1)*mth+ic ) = 
c     $                 grdgre ( (j1-1)*mth+j, (j2-1)*mth+ic ) +

                  indx1 = (j1-1) * mth + j
                  indx2 = (j2-1) * mth + ic
                  IF ( lhighn == 2 ) THEN
                     indx1 = ( 2*(j1-1)+(isco-1) ) * mth + j
                     indx2 = ( 2*(j2-1)+(iscs-1) ) * mth + ic
                  END IF

                  grdgre (indx1,indx2) = 
     $                 grdgre ( indx1,indx2 ) +
     $                 work(ic) * wgausn(j3)
c     
                  gren(j,ic) = gren(j,ic)
     $                 + worg(ic) * wgausn(j3) / factpi
c     
               end do
c     
c$$$  write ( outmod, '(2i3," >>> work: j, j3, thesg0, xs, wgaus",
c$$$  $           1p3e12.4, /, (1p10e12.4) )' )
c$$$  $           j,j3, thesg0, xs, wgausn(j3), (work(iw), iw=1, 20)
c$$$  write ( outmod, '(11x,"worg:", /, (1p10e12.4) )' )
c$$$  $           (worg(iw), iw=1, 20)
c     
 180        continue
c     
c$$$  write ( outmod, '(/,20x,"grdgre: j =" i3, / (1p10e12.4) )' )
c$$$  $           j, (grdgre((j1-1)*mth+j,(j2-1)*mth+iw), iw=1, 20)
c$$$  write ( outmod, '(20x,"gren: j =" i3, / (1p10e12.4) )' )
c$$$  $           j, (gren(j,iw), iw=1, 20)
c$$$  write ( outmod,  '(/)' )
c     
 200     continue
c     End of observer do-loop

c      Add in del**2 residue here  ! 06/16/00
 
         IF ( (j1 == j2) .AND. (iscs == isco) ) THEN
            indx0 = (j1-1) * mth
            IF ( lhighn == 2 ) indx0 = ( 2*(j1-1)+(isco-1) ) * mth 
            DO iker = 1, mth
               grdgre(indx0+iker,indx0+iker) =
     $              grdgre(indx0+iker,indx0+iker) + 1.0
            END DO
         END IF

         if ( check2 ) write ( 23,330 ) jres, jres2, j1,j2, ak0i,
     $        resdg, resk0, residu, j1,j2
 330     format ( /,1x,"jres, jres2, j1,j2, isgn*AK0I,  = ", 4i3,1pe14.5,/
     $        " Resdg, resk0, Residue added = ",
     $        1p3e11.2, " in ",2i2, "  block", / )
c     
         call atpoint ( "KERNEL3-END","j1", j1,"xobs(1)", xobs(1),
     $        iotty, outmod )
c     
         return
         end
c     
c.......................................................
      subroutine gausi4 ( the, thobs, xl,xu, xsce,zsce,
     $     xpp,zpp, isgn,iopw,iops,ising, aker,bker,aker0 )
c.......................................................
c     
c.... The is independent array of theta. 
c.....xl, xu is lower and upper limits of the gaussian integration
c     kernel K, G, stored in aker and bker.
c     aker0 is the integral over K^0.
c     
      include 'vacuum1.inc'
      include 'vacuum3.inc'
      include 'vacuum2.inc'

      COMMON / grvac / aval0, phex0, zphfc0, iscs, isco,
     $     j1com,j2com

      dimension  the(*), xsce(*), zsce(*)
      dimension tgaus4(4), wgaus4(4)
      dimension xpp(*), zpp(*), tab(3)
c     
c.....Weights for four points Gaussian quadrature         
c     
      wgaus4(1) = 0.347854845137454
      wgaus4(2) = 0.652145154862546
      wgaus4(3) = wgaus4(2)
      wgaus4(4) = wgaus4(1)
c     
      agaus4 = half*(xu+xl)
      bgaus4 = half*(xu-xl)
c     
      c1gaus4 = 0.861136311594053 * bgaus4
      c2gaus4 = 0.339981043584856 * bgaus4
c     
      tgaus4(1) = agaus4 - c1gaus4
      tgaus4(2) = agaus4 - c2gaus4
      tgaus4(3) = agaus4 + c2gaus4
      tgaus4(4) = agaus4 + c1gaus4
c     
c     thobs = the (iobs)
      aker = 0.0
      aker0 = 0.0
      bker = 0.0
c$$$  write ( outmod, '("iops, iopw, isgn = ", 3i4 )')
c$$$  $     iops, iopw, isgn

c  ZPHFAC is the phase stuff for the high-n calculation.

      do 160 ig = 1, 4
c     
         tgaus40 = tgaus4(ig)
         if ( tgaus40 .lt. zero ) tgaus40 = twopi + tgaus40
         if ( tgaus40 .ge. twopi ) tgaus40 = tgaus40 - twopi
         call spl1d2 ( mth1, the, xsce, xpp, 1, tgaus40, tab )
         xt = tab(1)
         xtp = tab(2)
         call spl1d2 ( mth1, the, zsce, zpp, 1, tgaus40, tab )
         zt = tab(1)
         ztp = tab(2)

         phext = 0.0
         zphfac = 1.0

         IF ( (lhighn == 2) .AND. (j2com == 1) ) THEN
            call spl1d2 ( mth1, the, qth0, qth0pp, 1, tgaus40, tab )
            zphex = tab(1)
            phext =  qa1 * ( zphex + tgaus4(ig) )
         END IF
         IF ( lhighn == 2 ) THEN
            dphex = phext - phex0
            IF ( iscs == 1 ) zphfac = cos (n*dphex)
            IF ( iscs == 2 ) zphfac = sin (n*dphex)
         END IF

c$$$         write ( outmod,'("ig, j2com, qa1,tgaus4(ig),phext = ",
c$$$     $        2i4,1x, 3es12.4)' ) ig, j2com, qa1, tgaus4(ig), phext 

c 157     CONTINUE
c     
         call green

         wgbg = wgaus4(ig) *  bgaus4
c     
         aker = aker + aval * zphfac * wgbg
c     
c.... If ising eq 1 then treat singularity
c     
c.... aker0 will contribute to the integral of K^0.
c     
         aker0 = aker0 + aval0 * wgbg
c     
c.....don't need gren if iopw is zero,
c     and subtract out log behavior if iops is one.
c     
         if ( iopw .eq. 0 ) go to 160
c     
c...  Don't use modulo 2pi on tgaus4 here!
c    The high-n factor ZPHFC0 was calculated at observer point.

         bval = bval * zphfac +  iops*ising
     $        * zphfc0*alog((thobs-tgaus4(ig))**2) / xs
c     
         bker = bker + bval * wgbg
c     
 160  continue
c     
      return
      end
c     c     
c.......................................................
      subroutine gausi6 ( the, thobs, xl,xu, xsce,zsce,
     $     xpp,zpp, isgn,iopw,iops,ising, aker,bker,aker0 )
c.......................................................
c     
c.... The is independent array of theta.
c.....xl, xu is lower and upper limits of the gaussian integration
c     kernel K, G, stored in aker and bker.
c     aker0 is the integral over K^0.
c     
      include 'vacuum1.inc'
      include 'vacuum3.inc'
      include 'vacuum2.inc'

      COMMON / grvac / aval0, phex0, zphfc0, iscs, isco,
     $     j1com,j2com

      dimension  the(*), xsce(*), zsce(*)
      dimension tgaus6(6), wgaus6(6)
      dimension xpp(*), zpp(*), tab(3)
c     
c.....Weights for six points Gaussian quadrature         
c     
      wgaus6(1) = 0.171324492379170
      wgaus6(2) = 0.360761573048139
      wgaus6(3) = 0.467913934572691
      wgaus6(4) = wgaus6(3)
      wgaus6(5) = wgaus6(2)
      wgaus6(6) = wgaus6(1)
c     
      agaus6 = half*(xu+xl)
      bgaus6 = half*(xu-xl)
c     
      c1gaus6 = 0.932469514203152 * bgaus6
      c2gaus6 = 0.661209386466265 * bgaus6
      c3gaus6 = 0.238619186083197 * bgaus6
c     
      tgaus6(1) = agaus6 - c1gaus6
      tgaus6(2) = agaus6 - c2gaus6
      tgaus6(3) = agaus6 - c3gaus6
      tgaus6(4) = agaus6 + c3gaus6
      tgaus6(5) = agaus6 + c2gaus6
      tgaus6(6) = agaus6 + c1gaus6

c     thobs = the (iobs)
      aker = 0.0
      aker0 = 0.0
      bker = 0.0
c$$$  write ( outmod, '("iops, iopw, isgn = ", 3i4 )')
c$$$  $     iops, iopw, isgn
c     
      DO 160 ig = 1, 6
c     
         tgaus60 = tgaus6(ig)
         IF ( tgaus60 .lt. zero ) tgaus60 = twopi + tgaus60
         IF ( tgaus60 .ge. twopi ) tgaus60 = tgaus60 - twopi
         call spl1d2 ( mth1, the, xsce, xpp, 1, tgaus60, tab )
         xt = tab(1)
         xtp = tab(2)
         call spl1d2 ( mth1, the, zsce, zpp, 1, tgaus60, tab )
         zt = tab(1)
         ztp = tab(2)

         phext = 0.0
         zphfac = 1.0

         IF ( (lhighn == 2) .AND. (j2com == 1) ) THEN
            call spl1d2 ( mth1, the, qth0, qth0pp, 1, tgaus60, tab )
            zphex = tab(1)
            phext =  qa1 * ( zphex + tgaus6(ig) )
         END IF
         IF ( lhighn == 2 ) THEN
            dphex = phext - phex0
            IF ( iscs == 1 ) zphfac = cos (n*dphex)
            IF ( iscs == 2 ) zphfac = sin (n*dphex)
         END IF

         call green

         wgbg = wgaus6(ig) * bgaus6
c     
         aker = aker + aval * zphfac * wgbg
c     
c.... If ising eq 1 then treat singularity
c     
c.... aker0 will contribute to the integral of K^0.
c     
         aker0 = aker0 + aval0 * wgbg
c     
c.....don't need gren if iopw is zero,
c     and subtract out log behavior if iops is one.
c     
         if ( iopw .eq. 0 ) go to 160
c     
c...  Don't use modulo 2pi on tgaus6 here!
c    The high-n factor ZPHFC0 was calculated at observer point.

         bval = bval * zphfac + iops*ising
     $        * zphfc0*alog((thobs-tgaus6(ig))**2) / xs
c     
         bker = bker + bval * wgbg
c     
 160  continue
c     
      return
      end
c     
c...............................................
      subroutine alg012 ( t0, t1, t2, alg )
c...............................................
c     
c.....Calculates the integral:- 
c     
c     int_{t1}^{t2} \log |t-t0| dt.
c     
      a20 = abs ( t2-t0 )
      a10 = abs ( t1-t0 )
c     
      alg = (t2-t0) * alog(a20) + (t0-t1) * alog(a10) - (t2-t1)
c     
      return
      end
c
c........................................................
      SUBROUTINE framep ( jobid, ff )
c.......................................................
c
      CHARACTER*(10) logo(2), ntim, mach, nsfx
      CHARACTER*(8) ndate
      CHARACTER*(*) jobid
c
      DATA icall /0/
c
c      placate compiler
c
      DATA mach,nsfx/" "," "/
c
c
c$$$      GO TO 100
c
      CALL map( 0.0, 1.0, 0.0, 1.0, 0.86, 1.0, 0.85, 1.0 )
cc    now plot the box
      call setcrt(0.,0.)
      call vector(0.,1.)
      call vector(1.,1.)
      call vector(1.,0.)
      call vector(0.,0.)
c     now enter the logo
      call setlch( 0.1, 0.75, 0, 1, 0, -1 )
      CALL date_and_time(ndate,ntim)
c$$$      call clock(ntim)
c$$$      logo(1) = ntim
c$$$      logo(2) = ndate
      CALL gtext ( ntim, -1, -1 )
      CALL setlch( 0.1, 0.55, 0, 1, -1, -1 )
      CALL gtext( ndate, -1, -1 )
      CALL setlch( 0.1, 0.35, 0, 1, -1, -1 )
      CALL gtext ( "vacuum", -1, -1 )
      CALL setlch( 0.02, 0.15, 2,0, -1, -1 )
      CALL gtext ( jobid, -1, -1 )
c     , 1 )
      CALL dders(1)
c
  100 CONTINUE
c
      CALL FRAME(1)
c
      RETURN
      END
c
c..................................................................
      subroutine coneig ( eigmat, ndim, eigval, l1,l2 )
c.................................................................
c     
      dimension eigmat(ndim,ndim), eigval(*), xcon(101), ycon(101)
      dimension c(101)
      character*80 strg
c     
      ncon = 21
      jmax1 = l2 - l1 + 1
c     
      zemin =  eigmat(1,1)
      zemax =  eigmat(1,1)
      zxmin = l1
      zxmax = l2
      zymin = 1.0
      zymax = jmax1
c     
      do i = 1, jmax1
         do j = 1, jmax1
            if ( zemin .gt. eigmat(i,j) ) zemin = eigmat(i,j)
            if ( zemax .lt. eigmat(i,j) ) zemax = eigmat(i,j)
         end do
         ll = l1 - 1 + i
         xcon(i) = ll
         ycon(i) = i
      end do
c     
      delx = zxmax - zxmin
      dely = zymax - zymin
      lsqre = 0
      if ( lsqre .ne. 1 ) go to 5
      zxmax = zxmin + amax1(delx,dely)
      zymax = zymin + amax1(delx,dely)
    5 continue
c     
      dzx = (zxmax-zxmin)/20.0
      dzy = (zymax-zymin)/20.0
      zxmin = zxmin - dzx
      zymin = zymin - dzy
      zxmax = zxmax + dzx
      zymax = zymax + dzy
c     
      call maps(zxmin,zxmax,zymin,zymax,0.142,0.858,0.285,1.0)
c     
c     do 15 i = 1, mth1
c     call pointc ( 1h+, xedg(i), zedg(i),1,-1,-1,0.0,0.0)
c     15 continue
c     
      k1 = 11
      k2 = 0.9
      delev = ( zemax-zemin )/(ncon-1)
      c(2) = delev
      do ik = 1, k1
         c(ik) = 0.5*zemax + (ik-1)*0.5*zemax / (k1-1)
      end do
c
      call rcontr(k1,c,k2,eigmat,ndim,xcon,1,jmax1,1,ycon,1,jmax1,1)
c     
      xwrt = zxmax + (zxmax-zxmin)/40.0
      ywrt = zymax - (zymax-zymin)/6.0
      call setlch ( xwrt,ywrt,0,0,0,-1)
      write (strg, 101 )
      call gtext(strg,-1,-1)
 101  format ( "E-Vectors" )
      do 17 i = 1, ncon
         c(i) = zemin + (i-1)*delev
         write(strg,16) i, c(i)
         call gtext(strg,-1,-1)
 16      format ( i3, 1pe12.4 )
 17   continue
      call frame(1)
c
      return
      end

c ...............................................................
      subroutine gatonorm ( vacin,nd1, gatovac, nd2, rgato,mfel,mth,
     $     qa1,twopi,taugam, zfac, zlabel )
c ...............................................................
c
c... Normalize vacuum matrix to GATO's norm:
c
      CHARACTER(*) :: zlabel
      CHARACTER(3) blanks
      dimension  vacin(nd1,nd1), gatovac(nd2,nd2)
c
      zfac =  rgato / (twopi*mfel)
c.. futher normalized with PI/(2*\mu_0 * R)
      zfac = zfac * 10**7 / (8.0*rgato)
      do m1 = 1, mfel
         do m2 = 1, mfel
            gatovac(m1,m2) = zfac * vacin(m1,m2) 
         end do
      end do
 
      blanks(1:3) = "   "
      lenzl = INDEX(zlabel, "   ", .false. ) - 1
c
      OPEN ( UNIT=36, FILE = 'vacgato', STATUS = "REPLACE" )

      WRITE ( 36, '(3x, "VACGATO: ",/,3x, a )' ) zlabel(1:lenzl)
      write ( 36, '("rgato, mfel,mth, qa1, taugam = ",/,
     $     1pe13.5, 2i4, 1p2e13.5 / )' ) rgato, mfel,mth, qa1, taugam
      do m1 = 1, mfel
         write ( 36,'(/, i5 )') m1
         write ( 36,'(10e15.7)') (gatovac(m1,m2), m2 = 1, mfel)
      end do
      CLOSE ( UNIT = 36 )
c
      return
      end

c..............................................................
      subroutine adjustb ( betin, betout, a, bw, cw, dw, xmaj, plrad,
     $     ishape )
c..............................................................
c     
c....... Adjusts for the skewing of the subtending bulge angle due to
c         elongation and triangularity. Only good for small angles.
c
      if ( ishape .eq. 310 ) then
         r0 = cw
         r  = a
      end if
c
      if ( ishape .eq. 210 ) then
         r0 = xmaj + cw*plrad
         r  = plrad * ( 1.0 + a - cw )
      end if
c
c...bypass triangularity...not very accurate???
c
c.......Adjust for triangularity...
c
c      x2 =  r * cos( betin )
c      x3 =  r * cos( betin + dw * sin(betin) )
c
c      bet2 = atan( x3*tan(betin)/x2 )
c
      bet2 = betin
c
c.......Adjust for elongation...
c
      betout = abs ( atan ( tan(bet2) / bw ) )
c
      return
      end
c
c...............................................................
      SUBROUTINE mafill ( vin, nd1, min, vout, nd2, mout )
c...............................................................

c..Fills in the excess diagonal elements of V_ll' with the circular 
c    Cylindrical values. If ltfil = 2 then the filled values starts 
c    smoothly from the endpoints of V_ll'
c  This gives V_ll' a skeleton so that the transformation to V_kk' 
c    has support. The higher |l|  values needs to be large. This 
c    avoids zero eigenvalues in the excess values of V_kk' Could 
c    perhaps just use the unit dyadic for the excess. 
c  Note that we assume that mfel > jmax1. 

      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum5.inc'
      
      REAL, DIMENSION(nd1,nd1) :: vin
      REAL, DIMENSION(nd2,nd2) :: vout
      
      vout = 0.0
      
      ln = lmin(1)
      lx = lmax(1)
      jmax1 = lx - ln + 1

c..Copy vin into vout in the appropriate position:      

      lnf = - mfel/2
      lxf = + mfel/2
      jmax1f = lxf - lnf + 1
      
      vout( ln-lnf+1:lx-lnf+1, ln-lnf+1:lx-lnf+1 )
     $     = vin(1:jmax1,1:jmax1)

c.. Now fill in excess diagonal elements with cylindrical values
c..  twopi2 is 2pi*2pi

      fpsqr = twopi2 / xmaj
      
      IF ( ltfil == 1 ) facm = 1.0

      IF ( ltfil == 2 ) THEN
c        Calibrate to fit continuously from VIN to VOUT at lmin
         zabl = (1.0+a)**(2.0*abs(ln))
         vinm = fpsqr * (zabl+1.0)/(zabl-1.0) / abs(ln)
         facm = vin(1,1) / vinm
         write ( outmod,
     $        '( /, 5x, "Calib. factor, facm = ", 1pe13.5 )' ) facm
         write ( iotty,
     $        '( /, 5x, "Calib. factor, facm = ", 1pe13.5 )' ) facm
      END IF

      DO la = 1, ln - lnf 
         ll = lnf + la - 1
         zabl = (1.0 + a) ** (2.0*abs(ll))
         vout(la,la) =
     $        facm * fpsqr * (zabl + 1.0) / (zabl - 1.0) / abs(ll)
      END DO

      IF ( ltfil == 1 ) facx = 1.0

      IF ( ltfil == 2 ) THEN
c        Calibrate to fit continuously from VIN to VOUT at lmax
         zabl = (1.0+a)**(2.0*abs(lx))
         vinx = fpsqr * (zabl+1.0)/(zabl-1.0) / abs(lx)
         facx = vin(jmax1,jmax1) / vinx
         write ( outmod,
     $        '( /, 5x, "Calib. factor, facx = ", 1pe13.5 )' ) facx
         write ( iotty,
     $        '( /, 5x, "Calib. factor, facx = ", 1pe13.5 )' ) facx
      END IF

      DO la = 1, lxf - lx 
         ll = lx + la
         zabl = ( 1.0 + a ) ** (2.0*abs(ll))
         vout(lx-lnf+1+la,lx-lnf+1+la) =
     $        facx * fpsqr * (zabl + 1.0) / (zabl - 1.0) / abs(ll)
      END DO
      
      END SUBROUTINE mafill

c..................................................................
      SUBROUTINE new_basis ( vin, nd1, mfelz,jmax1z, ltfil, zslt,zclt,
     $     nd01,nd02, nout1,nout2)
c...................................................................

      REAL, DIMENSION(nd1,nd1) :: vin, zrk1
      REAL, DIMENSION(nd01,nd02) :: zslt, zclt
      REAL, DIMENSION(:,:), ALLOCATABLE :: vout

      jmax1f = jmax1z

      IF ( ltfil /= 0 ) THEN

         lnf = - mfelz/2
         lxf = + mfelz/2
         jmax1f = lxf - lnf + 1

      call atpoint ( "Newbasis-1","jmax1f", jmax1f,
     $     "vin(1,1)",vin(1,1), 6, 23 )

         ALLOCATE ( vout(jmax1f,jmax1f) )

      call atpoint ( "Newbasis-2","mfelz", mfelz,
     $     "vin(1,1)",vin(1,1), 6, 23 )

         CALL mafill ( vin, nd1, min, vout, jmax1f, mout )

      call atpoint ( "Newbasis-3","nd1", nd1,
     $     "vin(1,1)",vin(1,1), 6, 23 )
         vin(1:jmax1f,1:jmax1f) = vout (1:jmax1f,1:jmax1f)

      call atpoint ( "Newbasis-4","nd1", nd1,
     $     "vin(1,1)",vin(1,1), 6, 23 )

         DEALLOCATE ( vout )
         
      call atpoint ( "Newbasis-5","nd1", nd1,
     $     "vin(1,1)",vin(1,1), 6, 23 )

      CALL matwrtn ( vin,nd1,nd1,1,1,jmax1f,jmax1f,8,8,"VFILLed",
     $     nout1,nout2 )

      END IF

      iopc = 1
      iops = 2
      CALL fotofi ( vin, zrk1, nd1, mfelz,jmax1f, zclt,nd01,nd02, iopc )
      CALL matwrtn ( zrk1,nd1,nd1,1,1,mfelz,mfelz,8,8,"V_NEW_C",
     $     nout1,nout2 )
      CALL fotofi ( vin, vin, nd1, mfelz,jmax1f,zslt,nd01,nd02, iops )
      CALL matwrtn ( vin,nd1,nd1,1,1,mfelz,mfelz,8,8,"V_NEW_S",
     $     nout1,nout2 )
c     
      vin(1:mfelz,1:mfelz) =
     $     vin(1:mfelz,1:mfelz) + zrk1(1:mfelz,1:mfelz)

      END SUBROUTINE new_basis
c     
c.........................................................
      SUBROUTINE fotofi ( vin,vout,nd1, mfelz,jmax1z,
     $     scnlth,nd01,nd02, iopsc )
c.........................................................
c     
c     Transforms the Vacuum matrix form Fourier-l space to finite
c     elements-i space.
c     
c     Use automatic arrays for wrk1, wrk2, wrk3
c     
      include 'vacuum1.inc'
c     
      dimension vin(nd1,nd1), vout(nd1,nd1), scnlth(nd01,nd02)
      REAL, DIMENSION(mfelz,jmax1z) :: wrk1, wrk2
      REAL, DIMENSION(jmax1z,mfelz) :: wrk3
c     
c      jmax1 = lmax(1) - lmin(1) + 1
      jmax1f = jmax1z
c     
c     Get left side transformation matrix, T_il:
      call tmat ( scnlth,nd01,nd02, wrk1,mfel,jmax1f, iopsc )
c     
c     Multiply V from the left by T:
      wrk2(1:mfel,1:jmax1f) =
     $     matmul ( wrk1(1:mfel,1:jmax1f),vin(1:jmax1f,1:jmax1f) )
c$$$  call matmul3 ( wrk1, vin, nd1,nd1, mfel,jmax1,jmax1, wrk2,nd1 ) 
c     
c     Get right side transformation rectangular matrix, T_li.
c     Transform whole square matrix including the zeros:
      wrk3 = transpose ( wrk1 )
c$$$  call mtrans ( wrk1, nd1, nd1 )
c     
c     Now get V in finite element basis:
      vout(1:mfel,1:mfel) =
     $     matmul ( wrk2(1:mfel,1:jmax1f),wrk3(1:jmax1f,1:mfel) )
c$$$  call matmul3 ( wrk2, wrk1, nd1,nd1, mfel,jmax1,mfel, vout,nd1 )
c     
      return
      end
c
c.........................................................
c      subroutine orchek ( wrkr, wrki, wrkrt, wrkit )
      subroutine orchek ( mfelz, jmax1z )
c.........................................................
c
c Checks for the orthogonality of the transformation matrix.
c Real and imaginary parts and their transpose are needed.
c
      include 'vacuum1.inc'
      include 'vacuum5.inc'
c
      REAL, DIMENSION(mfelz,jmax1z)  :: wrkr, wrki, wrkrt, wrkit
      REAL, DIMENSION(mfelz,mfelz)   :: wrko1, wrko2
      REAL, DIMENSION(jmax1z,jmax1z) :: wrkl1, wrkl2
      REAL, DIMENSION(jmax1z,mfelz)  :: wrkrtt, wrkitt
c$$$      dimension  wrkr(nfm,nfm), wrki(nfm,nfm),
c$$$     $     wrkrt(nfm,nfm), wrkit(nfm,nfm),
c$$$     $     wrko1(nfm,nfm), wrko2(nfm,nfm)
c     $     wrka(*),wrkb(*)
c      equivalence ( wrka(1),wrko1(1,1) ), ( wrkb(1),wrko2(1,1) )
c
      lnf = lmin(1)
      lxf = lmax(1)

      IF ( ltfil > 0 ) THEN
         lnf = - mfel / 2
         lxf = + mfel / 2
      END IF
      
      jmax1f = lxf - lnf + 1
c
c Get Real and Imaginary  matrix, T_im
      iopc = 1
      iops = 2
      nd1 = nths
      nd2 = nfm
      call tmat ( coslt,nd1,nd2, wrkrt,mfel,jmax1f, iopc )
      call tmat ( sinlt,nd1,nd2, wrkit,mfel,jmax1f, iops )
c
c....... Checks on the T matrix:
c
c      call matwrtn ( wrkrt,nfm,nfm,1,lmin(1),mfel,jmax1f,mfel,jdel,
c     $     "TMAT-COS", outmod,iotty )
c
      go to 8881
c$$$      write ( outmod, '(/,5x,"Sum over the finite els. for each  L:")')
c$$$c...... Sum the columns of TMAT for each row:
c$$$      do irow = 1, jmax1
c$$$         ll = lmin(1) - 1 + irow
c$$$         sumcol = 0.0
c$$$         do icol = 1, mfel
c$$$            sumcol = sumcol + wrkrt(icol,irow)
c$$$         end do
c$$$         write ( outmod, '(2x,"L  = ", i4, " Sum = ", e11.4)' )
c$$$     $        ll, sumcol
c$$$      end do
 8881 continue
c
c      call matwrtn ( wrkit,mfelz,jmax1z,1,lmin(1),mfel,jmax1f,mfel,jdel,
c     $     "TMAT-SIN", outmod,iotty )
c
      go to 8882
c$$$      write ( outmod, '(/,5x,"Sum over the finite els. for each  L:")')
c$$$c...... Sum the columns of TMAT for each row:
c$$$      do irow = 1, jmax1
c$$$         ll = lmin(1) - 1 + irow
c$$$         sumcol = 0.0
c$$$         do icol = 1, mfel
c$$$            sumcol = sumcol + wrkit(icol,irow)
c$$$         end do
c$$$         write ( outmod, '(2x,"L  = ", i4, " Sum = ", e11.4)' )
c$$$     $        ll, sumcol
c$$$      end do
 8882 continue
c
c Save Real and Imag. T_im
      wrkr(1:mfel,1:jmax1f) = wrkrt(1:mfel,1:jmax1f)
      wrki(1:mfel,1:jmax1f) = wrkit(1:mfel,1:jmax1f)
c
c$$$      do i = 1, mfel
c$$$         do m = 1, jmax1
c$$$            wrkr(i,m) = wrkrt(i,m)
c$$$            wrki(i,m) = wrkit(i,m)
c$$$         end do
c$$$      end do
c
c Get transpose of T:
      wrkrtt = transpose ( wrkrt )
      wrkitt = transpose ( wrkit )
c
c$$$      call mtrans ( wrkrt, nfm, nfm )
c$$$      call mtrans ( wrkit, nfm, nfm )
c
c Test for orthogonality in del-kk'. Real part first.
c
      wrko1(1:mfel,1:mfel) =
     $     matmul ( wrkr(1:mfel,1:jmax1f),wrkrtt(1:jmax1f,1:mfel) )
      wrko2(1:mfel,1:mfel) =
     $     matmul ( wrki(1:mfel,1:jmax1f),wrkitt(1:jmax1f,1:mfel) )
      wrko1(1:mfel,1:mfel) =
     $     wrko1(1:mfel,1:mfel) + wrko2(1:mfel,1:mfel)
c
c$$$      call matmul3 ( wrkr,wrkrt, nfm,nfm, mfel,jmax1,mfel, wrko1,nfm )
c$$$      call matmul3 ( wrki,wrkit, nfm,nfm, mfel,jmax1,mfel, wrko2,nfm )
c
c$$$      do i = 1, mfel
c$$$         do j = 1, mfel
c$$$            wrko1(i,j) = wrko1(i,j) + wrko2(i,j)
c$$$         end do
c$$$      end do
c
      call matwrtn ( wrko1,mfelz,mfelz,1,1,mfel,mfel,8,8,
     $     "Real TMAT*trans-TMAT: kk", outmod,iotty )
c
c Imaginary Part:
c$$$      call matmul3 ( wrkr,wrkit, nfm,nfm, mfel,jmax1,mfel, wrko1,nfm )
c$$$      call matmul3 ( wrki,wrkrt, nfm,nfm, mfel,jmax1,mfel, wrko2,nfm )
c$$$      do i = 1, mfel
c$$$         do j = 1, mfel
c$$$            wrko1(i,j) = wrko1(i,j) + wrko2(i,j)
c$$$         end do
c$$$      end do
c$$$      call matwrtn ( wrko1,nfm,nfm,1,1,mfel,mfel,8,8,
c$$$     $     "Imag TMAT*trans-TMAT: kk", outmod,iotty )
c
c
c Test for orthogonality in del-ll'. Real part first.
c
      wrkl1(1:jmax1f,1:jmax1f) =
     $     matmul ( wrkrtt(1:jmax1f,1:mfel),wrkr(1:mfel,1:jmax1f) )
      wrkl2(1:jmax1f,1:jmax1f) =
     $     matmul ( wrkitt(1:jmax1f,1:mfel),wrki(1:mfel,1:jmax1f) )
      wrkl1(1:jmax1f,1:jmax1f) =
     $     wrkl1(1:jmax1f,1:jmax1f) + wrkl2(1:jmax1f,1:jmax1f)
c
c$$$      call matmul3 ( wrkrt,wrkr, nfm,nfm, jmax1,mfel,jmax1, wrko1,nfm )
c$$$      call matmul3 ( wrkit,wrki, nfm,nfm, jmax1,mfel,jmax1, wrko2,nfm )
c$$$      do i = 1, jmax1
c$$$         do j = 1, jmax1
c$$$            wrko1(i,j) = wrko1(i,j) + wrko2(i,j)
c$$$         end do
c$$$      end do
c
      call matwrtn ( wrkl1,jmax1z,jmax1z,1,1,jmax1f,jmax1f,8,8,
     $     "Real TMAT*trans-TMAT: ll", outmod,iotty )
c
c Imaginary Part:
c$$$      call matmul3 ( wrkrt,wrki, nfm,nfm, jmax1,mfel,jmax1, wrko1,nfm )
c$$$      call matmul3 ( wrkit,wrkr, nfm,nfm, jmax1,mfel,jmax1, wrko2,nfm )
c$$$      do i = 1, jmax1
c$$$         do j = 1, jmax1
c$$$            wrko1(i,j) = wrko1(i,j) + wrko2(i,j)
c$$$         end do
c$$$      end do
c$$$      call matwrtn ( wrko1,nfm,nfm,1,1,jmax1,jmax1,8,8,
c$$$     $     "Imag TMAT*trans-TMAT: ll", outmod,iotty )
c
      return
      end
c     
c................................................................
      subroutine tmat ( sil,nd1,nd2, tll,mfelz,jmax1z, iop )
c................................................................
c     
c     Calculates the transformation matrix for converting the vacuum 
c     matrix from Fourier to finite element basis. 
c      Finite element is Sqrt N if t_k - pi/N <= t <= t_k + pi/N
c                                with    t_k = 2(k-1+dm0) * pi/N,  k = 1, N
c     
c
      include 'vacuum1.inc'
c     
      real nq
      dimension sil(nd1,nd2)
      REAL, DIMENSION (mfelz,jmax1z) :: tll
c     
c     Normalize finite elements to Sqrt(MFEL). Put in integration weights
c     for convenience.
c     
      znorm = sqrt( float(mfel) )
      jmax1 = lmax(1) - lmin(1) + 1
      pi = pye
c     
c      dm0 = 0.0
c
      tll = 0.0
c$$$      do i = 1, mfelz
c$$$         do j = 1, jmax1z
c$$$            tll(i,j) = 0.0 
c$$$         end do
c$$$      end do
c     
      if ( iop .ne. 0 ) go to 1000
c     
c..   Numerically:-
c     
      zwt1 = 0.5 * dth * znorm
      zwt2 = dth * znorm
      zws1 = dth * znorm /3.0
      zws2 = 2.0 * zws1
      zws4 = 2.0 * zws2
c     
      nzdel = ndfel
      nzdel1 = nzdel + 1
c     
c..   Loop to integrate over finite elements.
      do 135 l1 = 1, mfel
c     
         mzl = (l1-1) * nzdel + 1
         mzr = mzl + nzdel
c     
c..   Loop over obvserver integral.
         do 130 j = 1, nzdel1
c     
            jth0 = mzl + j - 2
            jth = mod(jth0,mth)
            if ( jth .lt. 0 ) jth = mth + jth
            jth1 = jth + 1
c     
c..   Loop over original fourier elements.
            do 140 l2 = 1, jmax1
               if ( nzdel .eq. 1 ) then
                  tll(l1,l2) = tll(l1,l2) + zwt1 * sil(jth1,l2)
               else
                  zwt = zws2
                  if ( (j/2)*2 .eq. j ) zwt = zws4
                  if ( (j .eq. 1) .or. (j .eq. nzdel1) ) zwt = zws1
                  tll(l1,l2) = tll(l1,l2) + zwt * sil(jth1,l2)
               end if
 140        continue
 130     continue
 135  continue
c     
      do mm = 1, mfel
         do ll = 1, jmax1
            tll(mm,ll) = tll(mm,ll) / twopi
         end do
      end do
c     
 1000 continue
c     
c..   Analytical expressions for T matrix:
c     
      znorm0 = 1.0 / znorm

c.... Try a different norm for discrete Fourier case:
c....   This is equivalent to a point transformation.

      IF ( ltnrm == 2 ) THEN
         znorm1 = znorm0
      END IF

c.. Fill in T matrix if ltfil > 0. 

      lnf = lmin(1)
      lxf = lmax(1)
      jmax1f = jmax1

      IF ( ltfil > 0 ) THEN
         lnf = - mfel/2
         lxf = + mfel/2
         jmax1f = lxf - lnf + 1
      END IF

      DO ll = 1, jmax1f
         la = lnf - 1 + ll
         IF ( ltnrm == 1 ) THEN
            IF ( la .ne. 0 )
     $           zslpn = znorm * sin(la*pi/mfel) / ( la * pi )
            znorm1 = zslpn
         END IF
         DO mm = 1, mfel
c..   Real part:-
            IF ( (iop .eq. 1) .and. (la .ne. 0) )
     $           tll(mm,ll) = cos(la*pi*2*(mm-1+dm0)/mfel) * znorm1
            IF ( (iop .eq. 1) .and. (la .eq. 0) )
     $           tll(mm,ll) = 1.0 * znorm0
c..   Imaginary part:-
            IF ( (iop .eq. 2) .and. (la .ne. 0 ) )
     $           tll(mm,ll) = sin(la*pi*2*(mm-1+dm0)/mfel) * znorm1
            IF ( (iop .eq. 2) .and. (la .eq. 0) )
     $           tll(mm,ll) = 0.0
         END DO
      END DO
c     
      return
      end
c
c........................................................
      subroutine wtopest ( vacpstr, vacpsti, xwal, zwal )
c........................................................
c
      include 'vacuum1.inc'
      
      dimension vacpstr(*), vacpsti(*), xwal(*), zwal(*)
c
      jmax1 = lmax(1) - lmin(1) + 1
      len = 2*nfm*nfm + 2*mth2 + 10
      ndsk = 1
      write (  iotty, '("Writing to unit iovac: iovac, jmax1, len = ",
     $     /, 3i10 )' ) iovac, jmax1, len
      write ( outmod, '("writing to unit iovac: iovac, jmax1, len = ",
     $     /, 3i10 )' ) iovac, jmax1, len
      call zop(iovac,"vacmat",len,ndsk,iiff,999)
      lmn = lmin(1)
      lmx = lmax(1)
      lgivup = 1
      nadres = 1
      call zwr(iovac,lmn,1,nadres,lgivup,999)
      nadres = nadres + 1
      call zwr(iovac,lmx,1,nadres,lgivup,999)
      nadres = nadres + 1
      call trans ( xwal,mth, xjdtxj,mthin )
      call zwr(iovac,xjdtxj(1),mthin2,nadres,lgivup,999)
      nadres = nadres + mthin2
      call trans ( zwal,mth, xjdtxj,mthin )
      call zwr(iovac,xjdtxj(1),mthin2,nadres,lgivup,999)
      nadres = nadres + mthin2
      length = jmax1**2
      call zwr(iovac,vacpstr(1),length,nadres,lgivup,999)
      nadres = nadres + jmax1**2
      call zwr(iovac,vacpsti(1),length,nadres,lgivup,999)
      write (  iotty, '("Wrote to iovac, lmn, lmx, nadres = ",
     $     /, 3i10 )' ) lmn, lmx, nadres
      write ( outmod, '("Wrote to iovac, lmn, lmx, nadres = ",
     $     /, 3i10 )' ) lmn, lmx, nadres
      call zcl( iovac, 999 )
c     
      return
 999  call errmes ( outpest, 'wtopest' )
      end
c
c*****************************************************************


      
