c......................................................................
      subroutine inglo
c......................................................................
c     this subroutine reads control logicals for the package
c     and opens necessary disk files..
c     
c......................................................................
c     
      include 'vacuum1.inc'
c     
      character*(80) string
      character*20 fort18, fort19
      nmap1 = 0
      nmpdsk = 0
      ss = 0.
c     
      if ( ladj .eq. 1 ) then
         open ( 50, file='inadjv', status='old', form='formatted' )
         return
      end if
c

      IF ( lrnim3d .eq. 1 ) return
      IF ( lrnim3d .eq. 2 ) return

      if ( ldcon .eq. 1 ) return
c
      if ( lrgato .gt. 0 ) return
c     
      if ( lzio .eq. 1 ) then
c     
c.......For Cray change file names to fort.18, fort.19 for wopen in zio...
c ....  For Hydra no change necessary. Subroutine SHELLB in Cray or Hydra
c       module.
c     
         call shellb
c
c$$$         fort18 = "fort.18"
c$$$         fort19 = "fort.19"
c$$$         write ( string, 115 ) mp0, fort19
c$$$         istat = ishell ( string )
c$$$         write ( string, 115 ) mp1, fort18
c$$$         istat = ishell ( string )
c$$$c     
c$$$ 115     format ( 'cp ',a, a )
c     
c         write ( iotty,  '( a "copied to ", a)' ) mp0, fort19
c         write ( outmod, '( a "copied to ", a)' ) mp0, fort19
c     
         call zop ( outmap1, mp1, nmap1, ndsk, ss, 100 )
         call zop ( iomode, mp0, nmpdsk, ndsk, ss, 100 )
c
      end if
c     
      return
 100  call errmes ( outmod, 'inglo' )
      end
c     
c     1.2   read data from cards.
c.......................................................................
c
      subroutine cardmo
c.......................................................................
      include 'vacuum1.inc'
      include 'vacuum2.inc'
      include 'vacuum3.inc'
      include 'vacuum7.inc'
c
c      real under
      CHARACTER(LEN=8) :: under
      data under / "--------" /
c
      namelist / modes  / lmax,lmin,lnint,lxint,mfel,m,mth,n,mdiv,
     $     mthout, xiinr,xiini, lsymz,lfunin,leqarcw, lsymwv,
     $     lpest1, lnova, ladj, ldcon, lrnim3d, lm3dc1, literr,
     $     lgato, lrgato, lrnova, lgaus,
     $     ltnrm, ltfil, lhighn, ldqdtw, lspark, ismth, lzio, mp0,mp1,
     $     dm0, dx0, dx1, s_scale
      namelist / rshel / lreshel, lwmn, lwmx, mshel, nwcur, taugam,
     $     gmat, irotplasma
      namelist / fcoil / drcoil, acoil, bcoil, dcoil, coilang, dzcoil,
     $     taucoil, coilcen, gainc, iopfc, icshp, mtcoil, ncgr,
     $     lfbcoil
      namelist / fcoila / drcoila, acoila, bcoila, dcoila, coilanga,
     $     dzcoila, taucoila, coilcena, gainca, iopfca, icshpa,
     $     mtcoila, lfbcoila
      namelist / fcoilb / drcoilb, acoilb, bcoilb, dcoilb, coilangb,
     $     dzcoilb, taucoilb, coilcenb, gaincb, iopfcb, icshpb,
     $     mtcoilb, lfbcoilb 
      namelist / scoil / isloop, asloop, sang1, sang2, nslpt,
     $     lsensplot, lrswcomp, lpolsensor
      namelist / debugs / checkd, lchecke, check1, check2, checks,
     $     wall, lkplt, lplot1, lplot2, lplot3, lplot4, lplot5
      namelist / vacdat / ishape,aw,awx,awv,bw,cw,dw,rw,tw,nsing,epsq,
     $     noutv,delg, idgt, idot, delfac, idsk, cn0
      namelist / shape  / ipshp, xpl, apl,bpl, dpl, epl, fpl, gpl, hpl,
     $     zpoff, a, b, r, abulg, bbulg, tbulg, qain,
     $     ldelta, adel
      namelist / diagns / lpestasym,
     $     lkdis, ieig, iloop, xloopnl, zloopnl,
     $     nloop,nloopr, ld3dflx, nxlpin, nzlpin, epslp,
     $     xlpmin, xlpmax, zlpmin, zlpmax,
     $     linterior, lchkinextp, lclosepts,
     .     lpsub, nphil, nphse, mx, mz, nph, xofsl,
     $     aloop, bloop, dloop, rloop, ntloop, deloop, ldisc
      namelist / sprk / nminus, nplus, mphi, lwrt11,civ,
     $     sp2sgn1, sp2sgn2, sp2sgn3, sp2sgn4, sp2sgn5,
     $     sp3sgn1, sp3sgn2, sp3sgn3, sp3sgn4, sp3sgn5,
     $     lff, ff, fv
c
c          checkd.... dumps vacuum matrices
c          lchecke....calculates eigenvalues and eigenvalues
c          check1....writes out pot,kin in delpla and minim.
c          symvac....forces vacuum matrix to be symmetric
c
c     1.2.1   job title
c             ... .....
c
      rewind inmode
      read ( inmode, 9100 )   (ntitle(i),i=1,20)
      write ( outmod, 9000 )   ntitle, ( under,i=1,2 )
c
c     1.2.2   read data from card.
c     ............................
c
      write ( outmod, 9001 )
c
 1221 continue
      rsave  = r
      write ( outmod,9002 )
c
      read(inmode,modes)
      write (*,'(" inmode after modes")')
      read(inmode,rshel)
      write (*,'(" inmode after rshel")')
      read(inmode,fcoil)
      write (*,'(" inmode after fcoil")')
      read(inmode,fcoila)
      write (*,'(" inmode after fcoila")')
      read(inmode,fcoilb)
      write (*,'(" inmode after fcoilb")')
      read(inmode,scoil)
      write (*,'(" inmode after scoil")')
      read(inmode,debugs)
      write (*,'(" inmode after debugs")')
      read(inmode,vacdat)
      write (*,'(" inmode after vacdat")')
      read(inmode,shape)
      write (*,'(" inmode after shape")')
      read(inmode,diagns)
      write (*,'(" inmode after diagns")')
      read(inmode,sprk)
      write (*,'(" inmode after sprk")')
c
      write(outmod,modes)
      write(outmod,rshel)
      write(outmod,fcoil)
      write(outmod,fcoila)
      write(outmod,fcoilb)
      write(outmod,scoil)
      write(outmod,debugs)
      write(outmod,vacdat)
      write(outmod,shape)
      write(outmod,diagns)
      write(outmod,sprk)
c
c..Flags On Input Data moved to VACUUM_MA.F
c
      r      = rsave
 1200 continue
c
c     output data just read in
c     ...... .... .... .... ...
c
c     write ( outmod, 9001 )
c     write ( outmod, modes )
      write ( outpest, 9003 )lmax(1),lmin(1),m,mdiv,n
c     write ( outmod, vacdat )
c     write ( outmod, shape )
c     write ( outmod, profil )
c     write ( outmod, debugs )
c
c     1.2.3   set any additional constants.
c     .....................................
c
c     f0     = gp
c     alphap = two * sqrt(p0)
      mp     = m + 1
      nosurf = ( mp - 1 ) * mdiv + 1
      mth1   = mth + 1
      mth2   = mth1 + 1
      mthout1 = mthout + 1
      mthout2 = mthout + 2
      no2pi  = n / twopi
      no2pi2 = no2pi * no2pi
      mw = mth1
      mtcoil1 = mtcoil + 1
      mtcoila1 = mtcoila + 1
      mtcoilb1 = mtcoilb + 1
c
      dth    = twopi / mth
      r2     = r * r
      r4     = r2 * r2
      r6     = r4 * r2
c
      open (  iodsk, file='vdata', status='new', form='formatted' )     
c
      lfour = 1
      lfele = 0
      if ( lgato .eq. 1 ) then
         lfour = 0
         lfele = 1
      end if
      if ( lgato .eq. 3 ) then
         lfour = 0
         lfele = 3
         mth = mfel
         mth1 = mth + 1 
         mth2 = mth + 2
         mw = mth1
         mfel1 = mfel + 1
         dmfel = twopi / mfel
         dth = twopi / mth
      end if
c
      lfarw = 0
      if ( a .gt. 10.0 ) lfarw = 1
c
c.. Save the values of lmin, lmax.
         lnsav = lmin(1)
         lxsav = lmax(1)
         lrnge = lmax(1) - lmin(1) + 1
c
      if ( lspark .ne. 0 ) then
c
         write ( iodsk,8001) mp1
         write ( iodsk,8002) 
         write ( iodsk,8002)
 8001    format ( a20 )
 8002    format ( a60 )
         write ( iodsk, 600 ) lmin(1), lmax(1), n
 600     format ( 2i5, 1pe12.5 )
c
c....Extend range of lmin, lmax for internal calculations.
c         lmin(1) = - lmax(1)
         lmin(1) = lnint
         lmax(1) = lxint
         lrnge = lmax(1) - lmin(1) + 1
c
      end if
      write ( iotty, '(/,"lreshel = ", i4 )' ) lreshel
      write ( outmod, '(/,"lreshel = ", i4 )' ) lreshel
c
      return
 9000 format ( 1x, 20a4, / 1x, 2a10 / )
 9001 format ( 1x, " form of data input" )
 9002 format ( 1x, " card data input" )
 9003 format(1x," lmax,lmin=",2i4," m,mdiv=",2i4,3x," n=",e12.4,/)
 9100 format ( 20a4 )
      end

c......................................................
c     
c     1.3   disk data input.
c.......................................................................
c     
      subroutine dskmd1
c.......................................................................
      include 'vacuum1.inc'
      include 'vacuum3.inc'
      include 'vacuum2.inc'
      integer t, surf
      real dat
c     
      dimension vecin(ntsin)

      REAL, DIMENSION(nths) :: zth
c     
c     1.3.1   check if disk data required.
c     
c     
c     1.3.2   read equilibrium quantities from disk storage.
c     first enough information to do the vacuum calculation
c     
c.......default delta and xjacob to zero.
c     
c Put lcdf in common and namelist later
c 
      lcdf = 0
c
      do 5 i = 1, mth2
         delta(i) = 0.0
         xjacob(i) = 0.0
 5    continue
c     
c......gato inputs...
c     
      if ( lrgato .eq. 0 ) go to 200
c
         lzio = 0
c
c.. GATO data is shifted by half gird point so:
c         dx0 = 0.5
c
c         if ( (lgato .eq. 1) .or. (lgato. eq. 2) ) then
         if ( lrgato .eq. 1 ) then
            CALL readvacin ( mfel,rgato,ndum2,ngato,qa1,xinf,zinf,
     $           delta, xigr,xigi, mth,mth1,mth2, ndfel,dx0,dx1,
     $           lgato,ieig, outmod, iotty )
            mfel1 = mfel + 1
            dth = twopi / mth
            dmfel = twopi / mfel
            n = ngato
         end if
c
c$$$         if ( lgato .eq. 3 ) then
c$$$            call readvacin ( mthin,rgato,ndum2,ngato,qa1,xinf,zinf,
c$$$     $           delta, xigr,xigi, mth,mth1,mth2, ndfel,dx0,dx1,
c$$$     $           lgato,ieig, outmod, iotty )
c$$$            mfel1 = mfel + 1
c$$$            dth = twopi / mth
c$$$            dmfel = twopi / mfel
c$$$            n = ngato
c$$$         end if
c
c... Set up arrays for plotting.
c
         jmax1 = lrnge
c
         do ii = 1, mth+1
            delta(ii) = 0.0
         end do
c
         write ( outmod, '(/,"******* DELTA set to 0.0 ******" )' )
         write ( iotty,  '(/,"******* DELTA set to 0.0 ******" )' )
c     
         write ( iotty, 701 ) mfel,rgato,ndum2,ngato, ga1,fa1,qa1
         write ( outmod,701 ) mfel,rgato,ndum2,ngato, ga1,fa1,qa1
 701     format (/, 'mfel, rgato,ndum2, ngato, ga1, fa1, qa1 = ',/,
     $        i5,1pe13.5,2i5,1p3e13.5, / )
c
      if ( (lrgato .eq. 2) .or. (lrgato .eq. 3)
     $        .OR. (lrgato .eq. 4) ) then
c
         lzio = 0
c
c.. GATO data is shifted by half gird point so:
c
         WRITE ( OUTMOD,
     $        '(/,5x,"DSKMD1: Shifts, DX0, DX1 = ", 2es11.3)' )
     $        dx0, dx1
         WRITE ( IOTTY,
     $        '(/,5x,"DSKMD1: Shifts, DX0, DX1 = ", 2es11.3)' )
     $        dx0, dx1
c
c$$$         if ( (lgato .eq. 1) .or. (lgato. eq. 2) ) then
c$$$            call readvacin ( mfel,rgato,ndum2,ngato,qa1,xinf,zinf,
c$$$     $           delta, xigr,xigi, mth,mth1,mth2, ndfel,dx0,dx1,
c$$$     $           lgato,ieig, outmod, iotty )
c$$$            mfel1 = mfel + 1
c$$$            dth = twopi / mth
c$$$            dmfel = twopi / mfel
c$$$            n = ngato
c$$$         end if
c
         IF ( lrgato == 2 ) THEN
            CALL readvacin2 ( mthin,rgato,ndum2,ngato,qa1,xinf,zinf,
     $           xwin, zwin, enqth,delta, xigr,xigi, bnkr,bnki,
     $           mth,mth1,mth2, ndfel,dx0,dx1,
     $           lgato,ieig, outmod, iotty )
            CALL atpoint ( "after readvacin2","mfel",mfel,"qa1",qa1,
     $           iotty, outmod )

c  Can calculate a q*delta from enqth. Call it qth0
            IF ( lhighn == 2 ) THEN
               zth(1:mth+2) = (twopi/mth) * (/ (i-1, i = 1, mth+2) /)
               qth0(1:mth1) =
     $              twopi - enqth(1:mth1) - zth(1:mth1)
               qth0(mth1) = qth0(1)
               qth0(mth2) = qth0(2)
            END IF
         END IF

         IF ( lrgato == 3 ) THEN
            CALL readvacin3 ( mthin,rgato,ndum2,ngato,qa1,xinf,zinf,
     $           xwin, zwin, enqth,delta, xigr,xigi, bnkr,bnki,
     $           mth,mth1,mth2, ndfel,dx0,dx1,
     $           lgato,ieig, outmod, iotty )
            CALL atpoint ( "after readvacin3","mfel",mfel,"qa1",qa1,
     $           iotty, outmod )
c  Can calculate a q*delta from enqth. Call it qth0
            IF ( lhighn == 2 ) THEN
               zth(1:mth+2) = (twopi/mth) * (/ (i-1, i = 1, mth+2) /) 
               qth0(1:mth1) =
     $              twopi - enqth(1:mth1) - zth(1:mth1)
               qth0(mth1) = qth0(1)
               qth0(mth2) = qth0(2)
            END IF
         END IF
 
         IF ( lrgato == 4 ) THEN
            CALL readvacin4 ( mthin,rgato,ndum2,ngato,qa1,xinf,zinf,
     $           xwin, zwin, enqth,delta, xigr,xigi, bnkr,bnki,
     $           mth,mth1,mth2, ndfel,dx0,dx1,
     $           lgato,ieig, outmod, iotty )
            CALL atpoint ( "after readvacin4","mfel",mfel,"qa1",qa1,
     $           iotty, outmod )
c  Can calculate a q*delta from enqth. Call it qth0
            IF ( lhighn == 2 ) THEN
               zth(1:mth+2) = (twopi/mth) * (/ (i-1, i = 1, mth+2) /) 
               qth0(1:mth1) =
     $              twopi - enqth(1:mth1) - zth(1:mth1)
               qth0(mth1) = qth0(1)
               qth0(mth2) = qth0(2)
            END IF
         END IF

         mfel1 = mfel + 1
         dth = twopi / mth
         dmfel = twopi / mfel
         n = ngato
c
c... Set up arrays for plotting.
c
         jmax1 = lrnge
c
         do ii = 1, mth+1
            delta(ii) = 0.0
         end do
c
         write ( outmod, '(/,"******* DELTA set to 0.0 ******" )' )
         write ( iotty,  '(/,"******* DELTA set to 0.0 ******" )' )
c
         write ( iotty, 701 ) mfel,rgato,ndum2,ngato, ga1,fa1,qa1
         write ( outmod,701 ) mfel,rgato,ndum2,ngato, ga1,fa1,qa1
c 701     format (/, 'mfel, rgato,ndum2, ngato, ga1, fa1, qa1 = ',/,
c     $        i5,1pe13.5,2i5,1p3e13.5, / )
c     
      end if
c
 200  continue  ! End of GATO inputs

c  Now input from NIMROD:

      IF ( lrnim3d == 1 ) THEN

         CALL readnim ( mthin, n_nim, xinf,zinf, mth, dx0,dx1 )
         
         n = n_nim
         mthin1 = mthin + 1
         mthin2 = mthin1 + 1
         mfel1 = mfel + 1
         dmfel = twopi / mfel
         dth = twopi / mth

         WRITE ( outmod, '(/,10x, "Nim3D inputs: n, mthin = ",
     $        f10.3, i5)' ) n, mthin
         WRITE ( iotty,  '(/,10x, "Nim3D inputs: n, mthin = ",
     $        f10.3, i5)' ) n, mthin

         jmax1 = lrnge
c
         do ii = 1, mth+1
            delta(ii) = 0.0
         end do
c
         write ( outmod, '(/,"******* DELTA set to 0.0 ******" )' )
         write ( iotty,  '(/,"******* DELTA set to 0.0 ******" )' )

      END IF

      IF ( lrnim3d == 2 ) THEN

c... Input DIII-D wall as the inner boundary for Nimrod. That is,
c... treat it as the plasma boundary.

         CALL d3dwall( xinf, zinf, mth, outmod, iotty )
         CALL eqarc ( xinf, zinf, mth1 )
         xinf(mth2) = xinf(2)
         zinf(mth2) = zinf(2)
         
         WRITE ( outmod, '(/,10x, "lrnim3d = ", i5)' ) lrnim3d
         WRITE (  iotty, '(/,10x, "lrnim3d = ", i5)' ) lrnim3d
c$$$         WRITE ( outmod, '(/, 1x,  "xinf = ", /(8es12.4) )' )
c$$$     $        (xinf(i), i = 1, mw, 10)
c$$$         WRITE ( outmod, '(/, 1x, "zinf = ", /(8es12.4) )' )
c$$$     $        (zinf(i), i = 1, mw, 10)
c$$$         WRITE ( iotty, '(/, 1x,  "xinf = ", /(8es12.4) )' )
c$$$     $        (xinf(i), i = 1, mw, 10)
c$$$         WRITE ( iotty, '(/, 1x,  "zinf = ", /(8es12.4) )' )
c$$$     $        (zinf(i), i = 1, mw, 10)

         DO ii = 1, mth+1
            delta(ii) = 0.0
         END DO
c
         WRITE ( OUTMOD, '(/,"******* DELTA set to 0.0 ******" )' )
         WRITE ( IOTTY,  '(/,"******* DELTA set to 0.0 ******" )' )

      END IF

      if ( lgato .ne. 0 ) then
c
         if ( (lgato .eq. 1) .or. (lgato .eq. 3) ) then
c..         Use same indexing for finite elements as Fourier.
c already set            lfele = 1
            lmin(1) = 1
            lmax(1) = mfel
         end if
c
c..  MTH already adjusted in READVACIN. If not then adjust it.
c
         if (  (lrgato .eq. 0) .and. (lgato .eq. 1) )
     $        call adjustm ( mth, mfel, mth1,mth2, ndfel, iotty,outmod )
         dth = twopi/mth
c
      write ( outmod, '(/,5x, "lgato, mth, mfel, ndfel, dmfel, dth = ",
     $     4i5, 2es13.5 )' ) lgato, mth, mfel, ndfel, dmfel, dth
      write ( iotty,  '(/,5x, "lgato, mth, mfel, ndfel, dmfel, dth = ",
     $     4i5, 2es13.5 )' ) lgato, mth, mfel, ndfel, dmfel, dth
c
      end if
c
c......adj inputs...
c     
      if ( ladj .eq. 1 ) then
c     
         read ( 50, 601 ) mthin1,lmin(1),lmax(1),nadj, qa1
 601     format ( 4i5, e13.5 )
c     
         mthin = mthin1 - 1
         mthin2 = mthin1 + 1
         n = nadj
c         fa1 = r*ga1 / qa1
c     
         write ( iotty, 602 ) mthin,lmin(1),lmax(1),nadj, ga1,fa1,qa1
         write ( outmod,602 ) mthin,lmin(1),lmax(1),nadj, ga1,fa1,qa1
 602     format ( /,'mthin, lmin,lmax, nadj, ga1, fa1, qa1 = ',/,
     $        4i5, 1p3e13.5, / )
c     
         read ( 50,605 ) ( vecin(i), i=1,mthin1 )
         call trans ( vecin,mthin, xinf,mth )
         read ( 50,605 ) ( vecin(i), i = 1,mthin1 )
         call trans ( vecin,mthin, zinf,mth )
         read ( 50,605 ) ( vecin(i), i = 1,mthin1 )
         call trans ( vecin,mthin, delta,mth )
c     
 605     format ( 10e13.5 )
c     
         close (50)
c     
         go to 1111
c     
      end if
c     
 1001 continue
c     
c......dcon inputs...
c     
      if ( ldcon .eq. 1 ) then
c
         lzio = 1
c
         CALL readahg ( mthin,lmin(1),lmax(1),ndcon,qa1,xinf,zinf,
     $     delta, mth )
c     
         mthin1 = mthin + 1
         mthin2 = mthin1 + 1
         n = ndcon
c         fa1 = r*ga1 / qa1
c
         write ( iotty, 702 ) mthin,lmin(1),lmax(1),ndcon, ga1,fa1,qa1
         write ( outmod,702 ) mthin,lmin(1),lmax(1),ndcon, ga1,fa1,qa1
 702     format (/, 'mthin, lmin,lmax, ndcon, ga1, fa1, qa1 = ',/,
     $        4i5, 1p3e13.5, / )
c     
c$$$c... Set DELTA = 0 for now:
c$$$
c$$$         delta(1:nths) = 2.0*delta(1:nths)
c$$$         WRITE ( iotty, '(/,"****** DELTA =2*DELTA *******",/)' )
c$$$         WRITE ( outmod, '(/,"****** DELTA = 2*DELTA *******",/)' )

         go to 1111
c     
      end if
c     
c......dcon-jkp inputs...
c     
      IF ( ldcon == 5 ) THEN
c
         lzio = 1
c
         dx0 = 0.0
         dx1 = 0.0
         CALL readvacin5 ( nxlpin,nzlpin, xlpmin,xlpmax, zlpmin,zlpmax,
     $        mthin, rscale, lmin(1),lmax(1),ntor, qa1,
     $     xinf,zinf, xwin,zwin, delta, bnlr,bnli,
     $     mth,mth1,mth2, dx0,dx1, lread, ieig, iotty, outmod )

         mthin1 = mthin + 1
         mthin2 = mthin1 + 1
         n = ntor

         WRITE ( IOTTY, 5702 ) mthin,lmin(1),lmax(1),n, qa1
         WRITE ( OUTMOD,5702 ) mthin,lmin(1),lmax(1),n, qa1
 5702    FORMAT (/, 'mthin, lmin,lmax, n, qa1 = ',/,
     $        3I5, 2ES13.5, / )

         GO TO 1111
C     
      END IF

c... M3dC1 Inputs

      IF ( lm3dc1 == 1 ) THEN
         IF ( ipshp == 1 ) GO TO 1111
         CALL readvacin_c1 ( mthin,ntor, xinf,zinf,mth,iotty,outmod )
         mthin1 = mthin + 1
         mthin2 = mthin + 2
         n = ntor
         delta = 0.0
         
         GO TO 1111

      END IF

c..... ITER error field input
c     
      IF ( literr == 1 ) THEN
c
c         lzio = 1
c
         dx0 = 0.0
         dx1 = 0.0
         CALL readvacin6c ( mthin, ntot, nmin, nmax, ntor, xinf,zinf,
     $        xwin,zwin,
     $        bnkr,bnki, mth,mth1,mth2, ieig, irepeat, iotty, outmod )

         mthin1 = mthin + 1
         mthin2 = mthin1 + 1
         n = ntor
         nrepeat = ntot

         WRITE ( IOTTY, 6702 ) mthin, ntot, n
         WRITE ( OUTMOD,6702 ) mthin, ntot, n
 6702    FORMAT (5x, 'mthin,  ntot, n, = ',/,5x 2I5, ES12.2 / )

c$$$         WRITE ( IOTTY, '("X, Z = ", /,(2ES12.4))' )
c$$$     $        (xinf(i), zinf(i), i = 1, mth2) 
c$$$         WRITE ( IOTTY, '("Br, Bi = ", /,(2ES12.4))' )
c$$$     $        (bnkr(i), bnki(i), i = 1, mth2) 

         GO TO 1111
C     
      END IF

      IF ( lrnova == 1 ) THEN

         dx0 = 0.0
         dx1 = 0.0
         CALL readvacin7 ( mthin, rdum, lmin(1),lmax(1),ntor, qa1,
     $     xinf,zinf, xwin,zwin, delta, bnlr,bnli,
     $     mth,mth1,mth2, dx0,dx1, lread, lreig7, lrwall7,
     $     iotty,outmod )

         mthin1 = mthin + 1
         mthin2 = mthin1 + 1
         n = ntor
         IF ( lreig7 == 1 ) ieig = 7
         WRITE ( IOTTY, 6703 ) lreig7, mthin, lmin(1), lmax(1), n, qa1
         WRITE ( OUTMOD,6703 ) lreig7, mthin, lmin(1), lmax(1), n, qa1
 6703    FORMAT (5x, 'lreig7, mthin,  lmin, lmax , n, qa1 = ',/,
     $        5x 4I5, 2ES12.2 / )

         GO TO 1111

      END IF

c     first, some constants.
c
      if ( lzio .eq. 1 ) then
         lgivup=1
         nadres = 1
         call zrd(iomode,ntitle(1),43,nadres,lgivup,999)
         write ( outmod, 9000 )   ntitle,dat
c
c+++++++ ++++++++ overwrite lj and zma to be consistent with old maps...
         lj = 0 
         zma = 0.0
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
         write ( iodsk, 8011 ) nosurf,mthin, lj,mj,nj, xzero, r,
     $        upsiln, xma, zma
 8011    format ( 5i4, 1p5e14.6 )
c
         r2     = r * r
         r4     = r2 * r2
         r6     = r4 * r2
c     
         mthin1 = mthin + 1
         mthin2 = mthin + 2
c     
      else
c         if ( lcdf .eq. 1 ) call mdskrd0
         continue
      end if
c     
c     then, all arrays.
c     read interface grid points...
c     
      if ( lzio .eq. 1 ) then
c     
         nadres = 50
         call zrd(iomode,vecin(1),mthin2,nadres,lgivup,999)
         write ( iodsk, 8021 ) ( vecin(i), i = 1, mthin2 )
 8021    format ( 1p10e14.6 )
         call trans ( vecin,mthin, xinf,mth )
         nadres = nadres + mthin2
         call zrd(iomode,vecin(1),mthin2,nadres,lgivup,999)
         write ( iodsk, 8021 ) ( vecin(i), i = 1, mthin2 )
         call trans ( vecin,mthin, zinf,mth )
c     
         length = ntsin * nsf
         ladres = 50 + 2*mthin2 + nosurf
         nadres = ladres + (nosurf-1)*ntsin
         lgivup = 1
         call zrd(iomode,vecin(1),mthin2,nadres,lgivup,999)
         write ( iodsk, 8021 ) ( vecin(i), i = 1, mthin2 )
         call trans ( vecin,mthin, grpssq,mth )
c     
         if ( .not. lpest1 ) then
c     
            nadres = nadres + 8*length
            call zrd(iomode,vecin(1),mthin2,nadres,lgivup,999)
         write ( iodsk, 8021 ) ( vecin(i), i = 1, mthin2 )
            call trans ( vecin,mthin, xjacob,mth )
            length = ntsin * nsf
            ladres = 50 + 2*mthin2 + nosurf
            nadres = ladres + (nosurf-1)*ntsin + 10*length
            call zrd(iomode,vecin(1),mthin2,nadres,lgivup,999)
         write ( iodsk, 8021 ) ( vecin(i), i = 1, mthin2 )
            call trans ( vecin,mthin, delta,mth )
c     
         end if
c     
         if ( lpest1 ) then
            do 100 i = 1, mth2
               xjacob(i) = upsiln * xinf(i)**2 / ( twopi*r )
 100        continue
         end if
c     
c     finally get poloidal flux values...
c     
c     
c     
c     read the equilibrium functions if mapping done before.
c     first ntitle, dat, nx, nz, alx,... etc.
c     nadres = 1
c     call zrd(outmap1,ntitle(1),43,nadres,lgivup,999)
         nadres = 50
         nadres = nadres + nosurf*3 - 1
         call zrd(outmap1,qa1,1,nadres,lgivup,999)
         nadres = nadres + nosurf*2
         call zrd(outmap1,ga1,1,nadres,lgivup,999)
         nadres = nadres + nosurf*2
         call zrd(outmap1,fa1,1,nadres,lgivup,999)
c     
         write ( outmod,9600 ) mp1, mp0, qa1, fa1, ga1
         write ( iotty, 9600 ) mp1, mp0, qa1, fa1, ga1
 9600    format ( /, 1x, " mp1, mp0 = " 2a, /
     $         1x, " qa1, fa1, ga1 = ", 1p3e12.5,/ )
c     
         write ( iodsk, 8031 ) qa1, ga1, fa1
 8031    format ( 1p3e14.6 )
c     
ccc   read array of theta values on the last magnetic surface
c     
         if ( lnova ) then
c     
            WRITE ( outmod, '("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")' )
            WRITE ( outmod, '("NOVA assumes up-down symmetry here!!")' )
            WRITE ( outmod, '("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")' )
            WRITE ( iotty, '("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")' )
            WRITE ( iotty, '("NOVA assumes up-down symmetry here!!")' )
            WRITE ( iotty, '("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")' )

            length = ntsin * nsf
            ladres = 50 + 2*mthin2 + nosurf
            nadres = ladres + (nosurf-1)*ntsin
            lgivup = 1
            call zrd(iomode,vecin(1),mthin2,nadres,lgivup,999)
            call trans ( vecin,mthin, grpssq,mth )
            nadres=nadres+length
            call zrd(iomode,vecin(1),mthin2,nadres,lgivup,999)
            call trans ( vecin,mthin, xsq,mth )
            nadres=nadres+5*length
            call zrd(iomode,vecin(1),mthin2,nadres,lgivup,999)
            call trans ( vecin,mthin, gpsdth,mth )
            nadres=nadres+length
            call zrd(iomode,vecin(1),mthin2,nadres,lgivup,999)
            call trans ( vecin,mthin, xsqdth,mth )
            nadres=nadres+length
            call zrd(iomode,vecin(1),mthin2,nadres,lgivup,999)
            call trans ( vecin,mthin, xjacob,mth )
c     
            mthd2p1=mth/2+1
            do 10 i=1,mthd2p1
               xsdtxs(i)=xsqdth(i)/xsq(i)
               gpdtgp(i)=gpsdth(i)/grpssq(i)
               xjdtxj(i)=0.5*(mj*xsdtxs(i)-nj*gpdtgp(i))
ccc   even parity
               xjacob(mth2-i)=xjacob(i)
ccc   odd parity
               xjdtxj(mth2-i)=-xjdtxj(i)
               delta(mth2-i)=-delta(i)
 10         continue
            xjdtxj(1)=0.
            xjdtxj(mthd2p1)=0.
            delta(1)=0.
            delta(mthd2p1)=0.
         end if
c     
c     close files
c     
c         call zcl ( outmap1, 999 )
c         call zcl ( iomode, 999 )
c     
 9000    format (1x,20a4,/,
     $        1x, " equilibrium from disk, calculated on date",a)
 9100    format(20a4,a10)
 9200    format(10i5)
 9300    format(4e20.13)
c     
c         jobid = mp1(3:)
c     
      else
c         if ( lcdf .eq. 1 ) call mdskrd1
         continue
      end if
c
 1111 continue
c
cc.....arbitrary plasma shape..
c     
      if ( ipshp .eq. 1 ) then
c     
         qa1  = qain
         rgato = xpl
c
         write ( outmod, 111 ) ipshp, qa1
         write ( iotty,  111 ) ipshp, qa1
 111     format ( /,"<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>",/,
     $        1x, "ipshp = ", i3, " qa1 = ", e13.5,/,
     $              "<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>",/ )
c.. +ve gpl values squares up outside corners
c   +ve hpl values rounds off inside corners.
c       epl redistributes points.
         DO  i = 1, mth2
            theta0 = (i-1) * dth
            theta = theta0 + epl*SIN(theta0)
            snthp =  SIN(theta)
            sn2thp = SIN(2*theta)
            xinf(i) = xpl + apl * ( COS(theta+dpl*snthp)
     $                    + fpl * sn2thp )
            zinf(i) =     - bpl* apl * SIN(theta+gpl*sn2thp)
     $                    - hpl * sn2thp + zpoff * apl
            delta(i) = 0.0
            IF ( ldelta == 1) then
               delta(i) = - adel * SIN (theta0)
            END IF
            xjacob(i) = 0.0
         END DO

      end if
c
c Plot some quantities.
c     
      write ( outmod,9500 ) nsf0, ntsin0, nths0, nfm, mtot,
     $     r, upsiln, mthin, nosurf
      write ( iotty,9500 ) nsf0, ntsin0, nths0, nfm, mtot,
     $     r, upsiln, mthin, nosurf
 9500 format ( /, 1x, "Mapping parameters, nsf0, ntsin0 = ", 2i5,/,
     $     1x, "Working parameter, nths0  = ", i5,/,
     $     1x, "nfm, mtot = ", 2i5,/,
     $     1x, "r, upsiln, mthin, nosurf = ", 1p2e12.5, 2i5,/ )
c
c      if ( lplot1 .eq. 0 ) go to 133
c
c..... plot x,z, for plasma versus theta.
c
      do i = 1, mth1
         pigrd(i) = (i-1)*dth
      end do
      call pospl1 ( pigrd,xinf, 1,mth1, 1,"xpla",mth1 )
      call pospl1 ( pigrd,zinf, 1,mth1, 2,"zpla",mth1 )
      call pospl1 ( pigrd,grpssq, 1,mth1, 3,"grpssq",mth1 )
      call pospl1 ( pigrd,delta, 1,mth1, 4,"delta",mth1 )
      call framep(jobid, ff)
c
 133  continue

      RETURN
 999  call errmes(outpest,'dskmd1')
      END
c     
c......................................................................
      SUBROUTINE efun ( zilr, zili, omsq )
c......................................................................
c     Reads in input from PEST-1. Comples xi.

      include 'vacuum1.inc'
      include 'vacuum3.inc'
      include 'vacuum6.inc'
      logical swall
      common /saved/saved,mdivs,msave,mat,swall
      dimension zilr(nfm),zili(nfm), lfm(nfm), zdel(nfm21)
cc%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NNN
c........  Hardwire lpestasym here.   NOW AN INPUT IN NAMELIST, April 14, 2011
c      lpestasym = 1 for up down asymmetry
c      lpestasym = 0 for symmetry
c      lpestasym = 0                    
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      do i = 1, nfm
         zilr(i) = 0.0
         zili(i) = 0.0
      end do
c
      call atpoint ( "Top of Efun", "mw", mw, "omsq",omsq,
     $     iotty, outmod )
c     
      length = nsf * ntsin
c      find the number of cases stored
      nadres = 50 + 8*nosurf + 6*length
      call zrd(outmap1,istore,1,nadres,1,999)
      if(istore.lt.1)then
      write(iotty,1000)
      write(outmod,1000)
 1000 format(" there are no eigenvectors stored in mapout1 ")
c      call exit
      stop
      end if
      do 100 is=1, istore
      iunst = is
      nadres = nadres + 1
      call zrd(outmap1,saved,5,nadres,1,999)
      mp = msave + 1
      if(swall)then
      write(iotty,1100)iunst
      write(outmod,1100)iunst
 1100 format(" case number ",i2," had a wall on the plasma edge",/)
      end if
      if(saved.ne.n)then
      write(iotty,1200)n,saved
      write(outmod,1200)n,saved
 1200 format(" mismatch in n value: input= ",f7.1," disk= ",f7.1)
c      call exit
      stop
      end if
      nadres = nadres + 5
      call zrd(outmap1,lminsv,1,nadres,1,999)
      if(lminsv.ne.lmin(1))then
      write(iotty,1300)lmin(1),lminsv
      write(outmod,1300)lmin(1),lminsv
 1300 format(" mismatch in lmin: input= ",i3," disk= ",i3)
c      call exit
      stop
      end if
      nadres = nadres + mp
      call zrd(outmap1,lmaxsv,1,nadres,1,999)
      if(lmaxsv.ne.lmax(1))then
      write(iotty,1400)lmax(1),lmaxsv
      write(outmod,1400)lmax(1),lmaxsv
 1400 format(" mismatch in lmax: input= ",i3," disk= ",i3)
c      call exit
      stop
      end if
      nadres = nadres + 2*mp
      call atpoint ( "A in EFUN", "mp", mp, "omsq-1",eigval(1),
     $     iotty, outmod )

      call zrd(outmap1,nunst,3,nadres,1,999)
      nadres = nadres + 3
      if(neigvl.eq.0) go to 100
      call zrd(outmap1,eigval(1),neigvl,nadres,1,999)
      call atpoint ( "B of Efun", "neigvc", neigvc, "eigval(1)",
     $     eigval(1), iotty, outmod )

      WRITE ( OUTMOD, '(/, "lminsav, lmaxsav, mat = ", 2i4, i10 )' )
     $     lminsv, lmaxsv, mat
      WRITE (  IOTTY, '(/, "lminsav, lmaxsav, mat = ", 2i4, i10 )' )
     $     lminsv, lmaxsv, mat

      write(outmod,1450)neigvl, nunst
      write( iotty,1450)neigvl, nunst
 1450 format(" there are (neigvl)",i2," eigenvalues stored",/,
     $     " nunst = ", i3,/ )
      write( iotty,1500)(eigval(i),i=1,neigvl)
      write(outmod,1500)(eigval(i),i=1,neigvl)
 1500 format(" the eigenvalues are : ",/,10(e10.3))
      nadres = nadres + neigvl
      if(nunst.eq.0)go to 100
      do 90 i=1, neigvc
      ik = i
      if((neigvc-ik).ge.nunst)go to 90
      jsub0 =  lmax(1)-lmin(1)+1
      jsub1 = ( lpestasym+1 ) * ( lmax(1)-lmin(1)+1 )  ! 09/17/2009
      jsub2 = 2*jsub1
      jsub4 = 4*jsub1           ! 08/18/2009 to acommodate complex xi
      nadres = nadres + mat - jsub2 - jsub1
c... Step back another 3*jusb1 for complex xi.
c$$$      nadres = nadres - 3*jsub1                  ! 8/31/2009 ! 09/17/2009
      call zrd(outmap1,zdel(1),jsub2,nadres,1,999) 
      nadres = nadres + jsub2   ! 08/18/2009 to acommodate complex xi
      il = 0
      do 10 ll=lmin(1),lmax(1)
      il=il+1
      ild=il+jsub1
      zilr(il) = zdel(ild) + twopi * ll * zdel(il)
      zili = 0.0
      IF ( lpestasym == 1 )
     $     zili(il) = zdel(ild+jsub0) + twopi * ll * zdel(il+jsub0)
      lfm(il)= ll
c
      write ( iotty,1600) ll, zilr(il), zili(il)
      write (outmod,1600) ll, zilr(il), zili(il)
 1600 format ( 1x, i4, 1p2e12.4 )
c
   10 continue
c
      omsq = eigval(i)
c
c$$$      if ( .not. swall ) then
c$$$c     
c$$$         psipr = 1.0
c$$$      do  l = 1, jsub1
c$$$         bnli(l) =   psipr * ( lfm(l)-n*qa1 ) * xilr(l)
c$$$         bnlr(l) = - psipr * ( lfm(l)-n*qa1 ) * xili(l)
c$$$      end do
c$$$c
c$$$      if ( lkdis ) call kdis ( omsq, bnlr,bnli, lfm )
c$$$c
c$$$c      if ( (iloop .ne. 0) .and. (lspark .eq. 0) )
c$$$      if ( iloop .ne. 0 )
c$$$     $     call pickup ( xilr, xili )
c$$$c
c$$$      end if
c
   90 continue
c
  100 continue
      return
  999 call errmes( outmod,'efun' )
      end


c......................................................................
      subroutine efunr ( zilr, zili, omsq )
c......................................................................
c...  Computes eddy currents for given perturbations
c       This is now renamed from efun. Replaced by efun, which is 
c        modified for complex xi.  08/18/2009
c.............
      include 'vacuum1.inc'
      include 'vacuum3.inc'
      include 'vacuum6.inc'
      logical swall
      common /saved/saved,mdivs,msave,mat,swall
      dimension zilr(nfm),zili(nfm), lfm(nfm), zdel(nfm21)
c
      do i = 1, nfm
         zilr(i) = 0.0
         zili(i) = 0.0
      end do
c
      call atpoint ( "To of Efun", "mw", mw, "omsq",omsq,
     $     iotty, outmod )
c     
      length = nsf * ntsin
c      find the number of cases stored
      nadres = 50 + 8*nosurf + 6*length
      call zrd(outmap1,istore,1,nadres,1,999)
      if(istore.lt.1)then
      write(iotty,1000)
      write(outmod,1000)
 1000 format(" there are no eigenvectors stored in mapout1 ")
c      call exit
      stop
      end if
      do 100 is=1, istore
      iunst = is
      nadres = nadres + 1
      call zrd(outmap1,saved,5,nadres,1,999)
      mp = msave + 1
      if(swall)then
      write(iotty,1100)iunst
      write(outmod,1100)iunst
 1100 format(" case number ",i2," had a wall on the plasma edge",/)
      end if
      if(saved.ne.n)then
      write(iotty,1200)n,saved
      write(outmod,1200)n,saved
 1200 format(" mismatch in n value: input= ",f7.1," disk= ",f7.1)
c      call exit
      stop
      end if
      nadres = nadres + 5
      call zrd(outmap1,lminsv,1,nadres,1,999)
      if(lminsv.ne.lmin(1))then
      write(iotty,1300)lmin(1),lminsv
      write(outmod,1300)lmin(1),lminsv
 1300 format(" mismatch in lmin: input= ",i3," disk= ",i3)
c      call exit
      stop
      end if
      nadres = nadres + mp
      call zrd(outmap1,lmaxsv,1,nadres,1,999)
      if(lmaxsv.ne.lmax(1))then
      write(iotty,1400)lmax(1),lmaxsv
      write(outmod,1400)lmax(1),lmaxsv
 1400 format(" mismatch in lmax: input= ",i3," disk= ",i3)
c      call exit
      stop
      end if
      nadres = nadres + 2*mp
      call atpoint ( "A in EFUN", "mp", mp, "omsq-1",eigval(1),
     $     iotty, outmod )

      call zrd(outmap1,nunst,3,nadres,1,999)
      nadres = nadres + 3
      if(neigvl.eq.0) go to 100
      call zrd(outmap1,eigval(1),neigvl,nadres,1,999)
      call atpoint ( "B of Efun", "neigvl", neigvl, "omsq-1",eigval(1),
     $     iotty, outmod )

      write(outmod,1450)neigvl, nunst
      write( iotty,1450)neigvl, nunst
 1450 format(" there are ",i2," eigenvalues stored",/,
     $     " nunst = ", i3,/ )
      write( iotty,1500)(eigval(i),i=1,neigvl)
      write(outmod,1500)(eigval(i),i=1,neigvl)
 1500 format(" the eigenvalues are : ",/,10(e10.3))
      nadres = nadres + neigvl
      if(nunst.eq.0)go to 100
      do 90 i=1, neigvc
      ik = i
      if((neigvc-ik).ge.nunst)go to 90
      jsub1 = lmax(1)-lmin(1)+1
      jsub2 = 2*jsub1
      nadres = nadres + mat - jsub2 - jsub1
      call zrd(outmap1,zdel(1),jsub2,nadres,1,999)
      nadres = nadres + jsub2
      il = 0
      do 10 ll=lmin(1),lmax(1)
      il=il+1
      ild=il+jsub1
      zilr(il) = zdel(ild) + twopi * ll * zdel(il)
      lfm(il)= ll
c
      write ( iotty,1600) ll, zilr(il), zili(il)
      write (outmod,1600) ll, zilr(il), zili(il)
 1600 format ( 1x, i4, 1p2e12.4 )
c
   10 continue
c
      omsq = eigval(i)
c
c$$$      if ( .not. swall ) then
c$$$c     
c$$$         psipr = 1.0
c$$$      do  l = 1, jsub1
c$$$         bnli(l) =   psipr * ( lfm(l)-n*qa1 ) * xilr(l)
c$$$         bnlr(l) = - psipr * ( lfm(l)-n*qa1 ) * xili(l)
c$$$      end do
c$$$c
c$$$      if ( lkdis ) call kdis ( omsq, bnlr,bnli, lfm )
c$$$c
c$$$c      if ( (iloop .ne. 0) .and. (lspark .eq. 0) )
c$$$      if ( iloop .ne. 0 )
c$$$     $     call pickup ( xilr, xili )
c$$$c
c$$$      end if
c
   90 continue
c
  100 continue
      return
  999 call errmes( outmod,'efun' )
      end

c.................................................................     
      SUBROUTINE readvacin5 ( nlx,nlz, xll,xlr, zlb,zlt,
     $     mthin, rscale, lmin,lmax,ntor, qa1,
     $     xinf,zinf, xwin,zwin, delta, bnlr,bnli,
     $     mth,mth1,mth2, dx0,dx1, lread, ireig5, nout1, nout2 )
c................................................................

c.... Reading DCON type input.   October 2, 2007

      CHARACTER*80 titvacin
      INTEGER mthin, ntor, lmin,lmax, ith,jl, jmax1, lread
      DIMENSION xinf(*), zinf(*), xwin(*),zwin(*),bnlr(*),bnli(*),
     $      delta(*)

      REAL, DIMENSION(:), ALLOCATABLE :: vecin

c-----------------------------------------------------------------------
c     read scalars.
c-----------------------------------------------------------------------
      OPEN ( UNIT=55, FILE='vacin5' )

      WRITE ( NOUT1, '(/"reading VACIN5 in READVACIN5",/)' )
      WRITE ( NOUT2, '(/"reading VACIN5 in READVACIN5",/)' )

      READ(55,'(a,/)') titvacin
      WRITE(NOUT1,'(a,/)') titvacin
      READ(55,*) nlx
      READ(55,*) nlz
      READ(55,*) xll
      READ(55,*) xlr
      READ(55,*) zlb
      READ(55,*) zlt
      READ(55,*) mthin
c$$$      READ(55,*) rscale
      READ(55,*) lmin
      READ(55,*) lmax
      READ(55,*) ntor
      READ(55,*) qa1

      mthin1 = mthin + 1
      jmax1 = lmax - lmin + 1

      ALLOCATE ( vecin(mthin+5) )

c-----------------------------------------------------------------------
c     read arrays.
c-----------------------------------------------------------------------

c.... Poloidal Coordinate Theta
      READ(55,'(a,/)') titvacin
      READ(55,*)(vecin(ith),ith=1,mthin1)
      WRITE(NOUT1,'(a,/)') titvacin

c.... Polatr Angle Eta
      READ(55,'(a,/)') titvacin
      READ(55,*)(vecin(ith),ith=1,mthin1)
      WRITE(NOUT1,'(a,/)') titvacin

c.... Xplasma
      READ(55,'(a,/)') titvacin
      READ(55,*)(vecin(ith),ith=1,mthin1)
      CALL transdx ( vecin,mthin, xinf,mth, dx0,dx1 )
      WRITE(NOUT1,'(a,/)') titvacin

c.... Zplasma
      READ(55,'(a,/)') titvacin
      READ(55,*)(vecin(ith),ith=1,mthin1)
      CALL transdx ( vecin,mthin, zinf,mth, dx0,dx1 )
      WRITE(NOUT1,'(a,/)') titvacin

c$$$c.... Xwall
c$$$      READ(55,'(a,/)') titvacin
c$$$      READ(55,*)(vecin(ith),ith=1,mthin1)
c$$$      CALL transdx ( vecin,mthin, xwin,mth, dx0,dx1 )
c$$$      WRITE(NOUT1,'(a,/)') titvacin

c$$$c.... Zwall
c$$$      READ(55,'(a,/)') titvacin
c$$$      READ(55,*)(vecin(ith),ith=1,mthin1)
c$$$      CALL transdx ( vecin,mthin, zwin,mth, dx0,dx1 )
c$$$      WRITE(NOUT1,'(a,/)') titvacin

c.... Delta
      READ(55,'(a,/)') titvacin
      READ(55,*)(vecin(ith),ith=1,mthin1)
      CALL transdx ( vecin,mthin, delta,mth, dx0,dx1 )
      WRITE(NOUT1,'(a,/)') titvacin

c$$$c.. Skip the J-B_pol and  Jacobian for now.      

c$$$c... J-B_pol
c$$$      READ(55,'(a,/)') titvacin
c$$$      READ(55,*)(vecin(ith),ith=1,mthin1)
c$$$      CALL transdx ( vecin,mthin, zinf,mth, dx0,dx1 )
c$$$      WRITE(NOUT1,'(a,/)') titvacin

c$$$c... Jacobian
c$$$      READ(55,'(a,/)') titvacin
c$$$      READ(55,*)(vecin(ith),ith=1,mthin1)
c$$$      CALL transdx ( vecin,mthin, zinf,mth, dx0,dx1 )
c$$$      WRITE(NOUT1,'(a,/)') titvacin

      IF ( ireig5 == 0 ) GO TO 200

c.... Real B_plasma_l
      READ(55,'(a,/)') titvacin
      READ(55,*)(bnlr(jl),jl=1,jmax1)
      WRITE(NOUT1,'(a,/)') titvacin

c.... Imag B_plasma_l
      READ(55,'(a,/)') titvacin
      READ(55,*)(bnli(jl),jl=1,jmax1)
      WRITE(NOUT1,'(a,/)') titvacin

 200  CONTINUE

      DEALLOCATE ( vecin )

      CLOSE ( UNIT=55 )
c-----------------------------------------------------------------------
c     terminate routine.
c-----------------------------------------------------------------------
      RETURN
      END

c.................................................................     
      SUBROUTINE readvacin7 ( mthin, rdum, lmin,lmax,ntor, qa1,
     $     xinf,zinf, xwin,zwin, delta, bnlr,bnli,
     $     mth,mth1,mth2, dx0,dx1, lread, lreig7, lrwall7,
     $     nout1, nout2 )
c................................................................

c.... Reading NOVA type input.  August 3, 2011

      CHARACTER*80 titvacin
      INTEGER mthin, ntor, lmin,lmax, ith,jl, jmax1, lread
      DIMENSION xinf(*), zinf(*), xwin(*),zwin(*),bnlr(*),bnli(*),
     $      delta(*)

      REAL, DIMENSION(:), ALLOCATABLE :: vecin

c-----------------------------------------------------------------------
c     read scalars.
c-----------------------------------------------------------------------
      OPEN ( UNIT=57, FILE='vacin7' )

      WRITE ( NOUT1, '(/"reading VACIN7 in READVACIN7",/)' )
      WRITE ( NOUT2, '(/"reading VACIN7 in READVACIN7",/)' )

      READ(57,'(a,/)') titvacin      ! Header
      WRITE(NOUT1,'(a,/)') titvacin
      READ(57,*) lrwall7             ! Equal 1 if wall 
      READ(57,*) lreig               ! Equal 1 if B_n 
      READ(57,*) mthin               ! Poloidal grid, open domain
c$$$      READ(57,*) rdum
      READ(57,*) lmin                ! Min Poloidal mode number
      READ(57,*) lmax                ! Max poloidal mode number
      READ(57,*) ntor                ! Toroidal mose number
      READ(57,*) qa1                 ! Safety factor at surface.

      mthin1 = mthin + 1
      jmax1 = lmax - lmin + 1

      ALLOCATE ( vecin(mthin+5) )

c-----------------------------------------------------------------------
c     read arrays.
c-----------------------------------------------------------------------

c$$$c.... Poloidal Coordinate Theta
c$$$      READ(57,'(a,/)') titvacin
c$$$      READ(57,*)(vecin(ith),ith=1,mthin1)
c$$$      WRITE(NOUT1,'(a,/)') titvacin
c$$$
c$$$c.... Polatr Angle Eta
c$$$      READ(57,'(a,/)') titvacin
c$$$      READ(57,*)(vecin(ith),ith=1,mthin1)
c$$$      WRITE(NOUT1,'(a,/)') titvacin

c.... Xplasma
      READ(57,'(a,/)') titvacin
      READ(57,*)(vecin(ith),ith=1,mthin1)
      CALL transdx ( vecin,mthin, xinf,mth, dx0,dx1 )
      WRITE(NOUT1,'(a,/)') titvacin

c.... Zplasma
      READ(57,'(a,/)') titvacin
      READ(57,*)(vecin(ith),ith=1,mthin1)
      CALL transdx ( vecin,mthin, zinf,mth, dx0,dx1 )
      WRITE(NOUT1,'(a,/)') titvacin

      IF ( lrwall7 == 0 ) GO TO 100
c.... Xwall
      READ(57,'(a,/)') titvacin
      READ(57,*)(vecin(ith),ith=1,mthin1)
      CALL transdx ( vecin,mthin, xwin,mth, dx0,dx1 )
      WRITE(NOUT1,'(a,/)') titvacin

c.... Zwall
      READ(57,'(a,/)') titvacin
      READ(57,*)(vecin(ith),ith=1,mthin1)
      CALL transdx ( vecin,mthin, zwin,mth, dx0,dx1 )
      WRITE(NOUT1,'(a,/)') titvacin

c.... Delta
      READ(57,'(a,/)') titvacin
      READ(57,*)(vecin(ith),ith=1,mthin1)
      CALL transdx ( vecin,mthin, delta,mth, dx0,dx1 )
      WRITE(NOUT1,'(a,/)') titvacin
 100  CONTINUE

      IF ( lreig7 == 0 ) GO TO 200

c.... Real B_plasma_l
      READ(57,'(a,/)') titvacin
      READ(57,*)(bnlr(jl),jl=1,jmax1)
      WRITE(NOUT1,'(a,/)') titvacin

c.... Imag B_plasma_l
      READ(57,'(a,/)') titvacin
      READ(57,*)(bnli(jl),jl=1,jmax1)
      WRITE(NOUT1,'(a,/)') titvacin

 200  CONTINUE

      DEALLOCATE ( vecin )

      CLOSE ( UNIT=57 )
c-----------------------------------------------------------------------
c     terminate routine.
c-----------------------------------------------------------------------
      RETURN
      END


c
c-----------------------------------------------------------------------
c     module 2. readahg.
c     reads boundary data for DCON.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE readahg ( mthin,lmin,lmax,ndcon,qa1,xinf,zinf,
     $     delta, mth )
c      implicit none
c
      integer mthin,nthin,lmin,lmax,ndcon,ith
      dimension xinf(*), zinf(*), delta(*)

      REAL, DIMENSION(:), ALLOCATABLE :: vecin

c-----------------------------------------------------------------------
c     read scalars.
c-----------------------------------------------------------------------
      open(unit=3,file='invacuum')
      read(3,*)mthin
      read(3,*)lmin
      read(3,*)lmax
      read(3,*)ndcon
      read(3,*)qa1
      mthin1 = mthin + 1

      ALLOCATE ( vecin(mthin+5) )

c-----------------------------------------------------------------------
c     read arrays.
c-----------------------------------------------------------------------
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call trans ( vecin,mthin, xinf,mth )
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call trans ( vecin,mthin, zinf,mth )
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call trans ( vecin,mthin, delta,mth )
c Try changing sign of delta
c      do 100 i = 1, mth+2
c         delta(i) = -delta(i)
c 100  continue
      close(unit=3)
c-----------------------------------------------------------------------
c     terminate routine.
c-----------------------------------------------------------------------

      DEALLOCATE ( vecin )

      return
      end
c     
c-----------------------------------------------------------------------
      SUBROUTINE readvacin_c1 ( mthin,nc1, xinf,zinf,mth,nout1,nout2 )
c......................................................................
c     implicit none

      INTEGER mthin, nc1
      DIMENSION xinf(*), zinf(*)

      REAL, DIMENSION(:), ALLOCATABLE :: vecinx, vecinz, id

c-----------------------------------------------------------------------
c     read scalars.
c-----------------------------------------------------------------------
      OPEN ( UNIT=3, FILE='vacin_c1' )

      READ( 3, '(i8)' ) nc1
      
      READ( 3, '(I8)') mthin

c-----------------------------------------------------------------------
c     read arrays.
c-----------------------------------------------------------------------

      ALLOCATE ( vecinx(mthin+5) )
      ALLOCATE ( vecinz(mthin+5) )
      ALLOCATE ( id(mthin+5) )

      DO i = 1, mthin
         READ ( 3, '(I8,2f12.8)' ) id(i), vecinx(i), vecinz(i)
      END DO

      WRITE(NOUT1, '("N_TOR, M3DC1= ", I8)') nc1
      WRITE(NOUT2, '("N_TOR, M3DC1= ", I8)') nc1

      WRITE(NOUT1, '("NODES, M3DC1= ", I8)') mthin
      WRITE(NOUT2, '("NODES, M3DC1= ", I8)') mthin

      DO i=1, mthin
         WRITE(NOUT1, '(I8,2f12.8)') id(i), vecinx(i), vecinz(i)
         WRITE(NOUT2, '(I8,2f12.8)') id(i), vecinx(i), vecinz(i)
      END DO

      mthin1 = mthin + 1
      mthin2 = mthin + 2

      vecinx(mthin1) = vecinx(1)
      vecinz(mthin1) = vecinz(1)

c..   Check for clockwise:
      IF ( vecinz(2) < vecinz(1) ) THEN
         WRITE ( NOUT1,
     $        '("****** Z1,Z2 from VACIN_C1 is CLOCKWISE")' )
         WRITE ( NOUT2,
     $        '("****** Z1,Z2 from VACIN_C1 is CLOCKWISE")' )
      ELSE
         WRITE ( NOUT1,
     $        '("****** Z1,Z2 from VACIN_C1 is COUNTER CLOCKWISE")' )
         WRITE ( NOUT2,
     $        '("****** Z1,Z2 from VACIN_C1 is COUNTER CLOCKWISE")' )
      END IF

c...  Now make theta go clockwise if not

      IF ( vecinz(2) > vecinz(1) ) THEN
         WRITE ( NOUT1,
     $        '("****** MAKING IT CLOCKWISE *****")' )
         WRITE ( NOUT2,
     $        '("****** MAKING IT CLOCKWISE *****")' )

         vecinx(1:mthin1) = (/ (vecinx(mthin1-i), i=0,mthin) /)
         vecinz(1:mthin1) = (/ (vecinz(mthin1-i), i=0,mthin) /)

         WRITE(NOUT1, '(/,"Clockwise:")') 
         WRITE(NOUT2, '(/,"Clockwise:")') 
         DO i=1, mthin1
            WRITE(NOUT1, '(I8,2f12.8)') i, vecinx(i), vecinz(i)
            WRITE(NOUT2, '(I8,2f12.8)') i, vecinx(i), vecinz(i)
         END DO

      END IF

c$$$  read(3,'(//)')

      CALL trans ( vecinx,mthin, xinf,mth )
      CALL trans ( vecinz,mthin, zinf,mth )

c     Try changing sign of delta
c     do 100 i = 1, mth+2
c     delta(i) = -delta(i)
c     100  continue

      CLOSE ( UNIT = 3 )
c-----------------------------------------------------------------------
c     terminate routine.
c-----------------------------------------------------------------------

      DEALLOCATE ( vecinx, vecinz, id )

      RETURN
      END

c-----------------------------------------------------------------------
c     reads boundary data for GATO or DCON. New Format.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE readvacin2 ( mthin,rgato,ndum2,ngato,qa1,xinf,zinf,
     $     xwin,zwin, enqth, delta, xigr,xigi, bnkr,bnki,
     $     mth,mth1,mth2, ndfel, dx0,dx1,
     $     lgato, ireig, nout1, nout2 )
c
c... Ming Chu's's format: Readvacin2 and Readvacin3 has same data
c
c      implicit none
c
      character*80 titvacin
      integer mthin,nthin,ndum1,ndum2,ngato,ith
      dimension xinf(*), zinf(*), xwin(*),zwin(*),bnkr(*),bnki(*),
     $     enqth(*), delta(*), xigr(*),xigi(*)

      REAL, DIMENSION(:), ALLOCATABLE :: vecin

c-----------------------------------------------------------------------
c     read scalars.
c-----------------------------------------------------------------------
      open(unit=3,file='vacin',form='formatted',status='old')
      write ( nout1, '(/"reading VACIN in VACIN2",/)' )
      write ( nout2, '(/"reading VACIN in VACIN2",/)' )
      read(3,'(a,/)')titvacin
      read(3,'(i10)')mthin
      read(3,'(1pe15.8)')rgato
      read(3,'(i10)')ndum2
      read(3,'(i10)')ngato
      read(3,'(1pe15.8)')qa1
      write(nout1,'(a,/)')titvacin
      write(nout1,'(i10)')mthin
      write(nout1,'(1pe15.8)')rgato
      write(nout1,'(i10)')ndum2
      write(nout1,'(i10)')ngato
      write(nout1,'(1pe15.8)')qa1
      mthin1 = mthin + 1
c
c..   Adjust mth for finite elements, if lgato .eq. 1
      if ( lgato .eq. 1)
     $     call adjustm ( mth, mthin, mth1,mth2,  ndfel, nout1,nout2 )

      ALLOCATE ( vecin(mthin+5) )

c-----------------------------------------------------------------------
c     read arrays.
c-----------------------------------------------------------------------
      read(3,'(a,/)')titvacin
      read(3,'(1p4e15.8)')(vecin(ith),ith=1,mthin1)
      write(nout1,'(a,/)')titvacin
      write(nout1,'(1p4e15.8)')(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, xinf,mth, dx0,dx1 )
      read(3,'(a,/)')titvacin
      read(3,'(1p4e15.8)')(vecin(ith),ith=1,mthin1)
      write(nout1,'(a,/)')titvacin
      write(nout1,'(1p4e15.8)')(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, zinf,mth, dx0,dx1 )
      read(3,'(a,/)')titvacin
      read(3,'(1p4e15.8)')(vecin(ith),ith=1,mthin1)
      write(nout1,'(a,/)')titvacin
      write(nout1,'(1p4e15.8)')(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, enqth,mth, dx0,dx1 )
      read(3,'(a,/)')titvacin
      read(3,'(1p4e15.8)')(vecin(ith),ith=1,mthin1)
      write(nout1,'(a,/)')titvacin
      write(nout1,'(1p4e15.8)')(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, xwin,mth, dx0,dx1 )
      read(3,'(a,/)')titvacin
      read(3,'(1p4e15.8)')(vecin(ith),ith=1,mthin1)
      write(nout1,'(a,/)')titvacin
      write(nout1,'(1p4e15.8)')(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, zwin,mth, dx0,dx1 )
c
c.. Skip the PEST angle stuff for now.      
c      read(3,'(//)')
c      read(3,*)(vecin(ith),ith=1,mthin1)
c
      if ( ireig .eq. 0 ) go to 200
c
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, xigr,mth, dx0,dx1 )
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, xigi,mth, dx0,dx1 )
      read(3,'(//)')
c..Reading Magnetic Perturbation. Don't translate and interpolate
c   until after vacuum energy is calculated.
      read(3,*)(bnkr(ith),ith=1,mthin1)
c      call transdx ( vecin,mthin, bnkr,mth, dx0,dx1 )
      read(3,'(//)')
      read(3,*)(bnki(ith),ith=1,mthin1)
c      call transdx ( vecin,mthin, bnki,mth, dx0,dx1 )
c
 200  continue
c Try changing sign of delta
c      do 100 i = 1, mth+2
c         delta(i) = -delta(i)
c 100  continue
      close(unit=3)

      DEALLOCATE ( vecin )

c-----------------------------------------------------------------------
c     terminate routine.
c-----------------------------------------------------------------------
      RETURN
      END
c
c-----------------------------------------------------------------------
c     reads boundary data for GATO or DCON. New Format.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE readvacin3 ( mthin,rgato,ndum2,ngato,qa1,xinf,zinf,
     $     xwin,zwin, enqth, delta, xigr,xigi, bnkr,bnki,
     $     mth,mth1,mth2, ndfel, dx0,dx1,
     $     lgato, ireig, nout1, nout2 )
c
c... Alan Turnbull's format: Readvacin2 and Readvacin3 has same data
c
c      implicit none
c
      integer mthin,nthin,ndum1,ndum2,ngato,ith
      dimension xinf(*), zinf(*), xwin(*),zwin(*),bnkr(*),bnki(*),
     $     enqth(*), delta(*), xigr(*),xigi(*)

      REAL, DIMENSION(:), ALLOCATABLE :: vecin

c-----------------------------------------------------------------------
c     read scalars.
c-----------------------------------------------------------------------
      open(unit=3,file='vacin')
      write ( nout1, '(/"reading VACIN in VACIN3",/)' )
      write ( nout2, '(/"reading VACIN in VACIN3",/)' )
      read(3,'(//)')
      read(3,*)mthin
      read(3,*)rgato
      read(3,*)ndum2
      read(3,*)ngato
      read(3,*)qa1
      mthin1 = mthin + 1
c
c..   Adjust mth for finite elements, if lgato .eq. 1
      if ( lgato .eq. 1)
     $     call adjustm ( mth, mthin, mth1,mth2,  ndfel, nout1,nout2 )

      ALLOCATE ( vecin(mthin+5) )

c-----------------------------------------------------------------------
c     read arrays.
c-----------------------------------------------------------------------
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, xinf,mth, dx0,dx1 )
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, zinf,mth, dx0,dx1 )
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, xwin,mth, dx0,dx1 )
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, zwin,mth, dx0,dx1 )
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, enqth,mth, dx0,dx1 )
c
c.. Skip the PEST angle stuff for now.      
c      read(3,'(//)')
c      read(3,*)(vecin(ith),ith=1,mthin1)
c
      if ( ireig .eq. 0 ) go to 200
c
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, xigr,mth, dx0,dx1 )
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, xigi,mth, dx0,dx1 )
      read(3,'(//)')
c..Reading Magnetic Perturbation. Don't translate and interpolate
c   until after vacuum energy is calculated.
      read(3,*)(bnkr(ith),ith=1,mthin1)
c      call transdx ( vecin,mthin, bnkr,mth, dx0,dx1 )
      read(3,'(//)')
      read(3,*)(bnki(ith),ith=1,mthin1)
c      call transdx ( vecin,mthin, bnki,mth, dx0,dx1 )
c
 200  continue
c Try changing sign of delta
c      do 100 i = 1, mth+2
c         delta(i) = -delta(i)
c 100  continue
      close(unit=3)

      DEALLOCATE ( vecin )

c-----------------------------------------------------------------------
c     terminate routine.
c-----------------------------------------------------------------------
      return
      end

      SUBROUTINE readvacin4 ( mthin,rgato,ndum2,ngato,qa1,xinf,zinf,
     $     xwin,zwin, enqth, delta, xigr,xigi, bnkr,bnki,
     $     mth,mth1,mth2, ndfel, dx0,dx1,
     $     lgato, ireig, nout1, nout2 )
c
c... Alan Turnbull's new format June 27, 2000
c    Now contains J*B_pol**2, Jacobian and delta-psi
c      implicit none
c
      integer mthin,nthin,ndum1,ndum2,ngato,ith
      dimension xinf(*), zinf(*), xwin(*),zwin(*),bnkr(*),bnki(*),
     $     enqth(*), delta(*), xigr(*),xigi(*)

      REAL, DIMENSION(:), ALLOCATABLE :: vecin

c-----------------------------------------------------------------------
c     read scalars.
c-----------------------------------------------------------------------
      open(unit=3,file='vacin')
      write ( nout1, '(/"reading VACIN in VACIN3",/)' )
      write ( nout2, '(/"reading VACIN in VACIN3",/)' )
      read(3,'(//)')
      read(3,*)mthin
      read(3,*)rgato
      read(3,*)ndum2
      read(3,*)ngato
      read(3,*)qa1
      mthin1 = mthin + 1
c
c..   Adjust mth for finite elements, if lgato .eq. 1
      if ( lgato .eq. 1)
     $     call adjustm ( mth, mthin, mth1,mth2,  ndfel, nout1,nout2 )

      ALLOCATE ( vecin(mthin+5) )

c-----------------------------------------------------------------------
c     read arrays.
c-----------------------------------------------------------------------
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, xinf,mth, dx0,dx1 )
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, zinf,mth, dx0,dx1 )
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, xwin,mth, dx0,dx1 )
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, zwin,mth, dx0,dx1 )
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, enqth, mth, dx0,dx1 )
c
c.. Skip the J-B_pol and  Jacobian for now.      
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)

      if ( ireig .eq. 0 ) go to 200
c
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, xigr,mth, dx0,dx1 )
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, xigi,mth, dx0,dx1 )

c   Skip the perturbed flux for now
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)

      read(3,'(//)')
c..Reading Magnetic Perturbation. Don't translate and interpolate
c   until after vacuum energy is calculated.
      read(3,*)(bnkr(ith),ith=1,mthin1)
c      call transdx ( vecin,mthin, bnkr,mth, dx0,dx1 )
      read(3,'(//)')
      read(3,*)(bnki(ith),ith=1,mthin1)
c      call transdx ( vecin,mthin, bnki,mth, dx0,dx1 )
c
 200  continue
c Try changing sign of delta
c      do 100 i = 1, mth+2
c         delta(i) = -delta(i)
c 100  continue
      close(unit=3)

      DEALLOCATE ( vecin )

c-----------------------------------------------------------------------
c     terminate routine.
c-----------------------------------------------------------------------
      return
      end

c.........................................................................
      SUBROUTINE readvacin6c ( mthin, ntot, nmin, nmax, ntor, xinf,zinf,
     $     xwin,zwin,
     $     bnkr,bnki, mth,mth1,mth2, ireig, izrepeat, nout1, nout2 )
c........................................................................

c Reads ITER saddle loops and magnetic field data.

      INTEGER :: mthin, nthin, ndum1, ndum2, ngato, ith
      CHARACTER*(80) titvacin
      REAL, DIMENSION(*) :: xinf, zinf, bnki, bnkr, xwin, zwin

      REAL, DIMENSION(:), ALLOCATABLE :: vecin6r, vecin6i,
     $     vecout6r, vecout6i


c-----------------------------------------------------------------------
c     read scalars.
c-----------------------------------------------------------------------

c...Only open vacin6_1 on the frist pass...
      If ( izrepeat == 0 )  OPEN( UNIT = 661, FILE = 'vacin6_1')

      WRITE ( nout1, '(/"***** reading vacin6_1 in READVACIN6C",/)' )
      WRITE ( nout2, '(/"***** reading vacin6_1 in READVACIN6C",/)' )

      IF ( izrepeat == 0 ) READ(661,*) mthin

c-----------------------------------------------------------------------
c     read X,Z arrays.
c-----------------------------------------------------------------------

      ALLOCATE ( vecin6r(mthin+5), vecout6r(mth+5),
     $     vecin6i(mthin+5), vecout6i(mth+5) )

      IF ( izrepeat > 0 ) GO TO 7771 ! need only n and fields

      READ( 661,'(A)') titvacin

      WRITE ( nout1, '("mthin = ", I5, A ,/)' ) mthin, titvacin
      WRITE ( nout2, '("mthin = ", I5, A ,/)' ) mthin, titvacin

      DO i = 1, mthin
         READ(661, * ) vecin6r(i), vecin6i(i)
      END DO

      dx0 = 0.0
      dx1 = 0.0
      CALL transdx ( vecin6r,mthin, vecout6r,mth, dx0,dx1 )
      CALL transdx ( vecin6i,mthin, vecout6i,mth, dx0,dx1 )
      vecout6r(mth1) = vecout6r(1)
      vecout6r(mth2) = vecout6r(2)
      vecout6i(mth1) = vecout6i(1)
      vecout6i(mth2) = vecout6i(2)

      xinf(1:mth2) = vecout6r(1:mth2)
      zinf(1:mth2) = vecout6i(1:mth2)

      READ(661, '(2I5)' ) nmin, nmax
      ntot = nmax - nmin + 1
 
 7771 CONTINUE

      READ(661,*) ntor
      n = ndum

      IF ( ireig == 0 ) GO TO 200

c-----------------------------------------------------------------------
c     read arrays of b_normal values.
c-----------------------------------------------------------------------
      READ( 661,'(A)') titvacin

      WRITE ( nout1, '("mthin, nmin, nmax, ntor = ", 4I5, A /)' )
     $     mthin, nmin, nmax, ntor, titvacin
      WRITE ( nout2, '("mthin, nmin, nmax, ntor = ", 3I5, A /)' )
     $     mthin, nmin, nmax, ntor, titvacin

      DO i = 1, mthin
         READ(661, * ) vecin6r(i), vecin6i(i)
      END DO

      CALL transdx ( vecin6r,mthin, vecout6r,mth, dx0,dx1 )
      CALL transdx ( vecin6i,mthin, vecout6i,mth, dx0,dx1 )
      vecout6r(mth1) = vecout6r(1)
      vecout6r(mth2) = vecout6r(2)
      vecout6i(mth1) = vecout6i(1)
      vecout6i(mth2) = vecout6i(2)

      bnkr(1:mth2) = vecout6r(1:mth2)
      bnki(1:mth2) = vecout6i(1:mth2)

c$$$      CALL transdxc ( vecin6c,mthin, vecout6c,mth, dx0,dx1 )
c$$$      vecout6c(mth1) = vecout6c(1)
c$$$      vecout6c(mth2) = vecout6c(2)
c$$$      bnkr(1:mth2) = REAL( vecout6c(1:mth2) )
c$$$      bnki(1:mth2) = AIMAG( vecout6c(1:mth2) )

      DEALLOCATE ( vecin6r, vecin6i, vecout6r, vecout6i )

c$$$
c$$$      if ( ireig .eq. 0 ) go to 200
c$$$
c$$$      READ(3,'(//)')
c$$$c..Reading Magnetic Perturbation. Don't translate and interpolate
c$$$c   until after vacuum energy is calculated.
c$$$      READ(3,*)(bnkr(ith),ith=1,mthin1)
c$$$c      CALL transdx ( vecin,mthin, bnkr,mth, dx0,dx1 )
c$$$      READ(3,'(//)')
c$$$      READ(3,*)(bnki(ith),ith=1,mthin1)
c$$$c      CALL transdx ( vecin,mthin, bnki,mth, dx0,dx1 )
c
 200  CONTINUE

c$$$      CLOSE ( UNIT = 661 )

c-----------------------------------------------------------------------
c     terminate routine.
C-----------------------------------------------------------------------
      RETURN
      END


c     reads boundary data for GATO or DCON.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE readvacin ( mthin,rgato,ndum2,ngato,qa1,xinf,zinf,
     $     delta, xigr,xigi, mth,mth1,mth2, ndfel, dx0,dx1,
     $     lgato, ireig, nout1, nout2 )
c      implicit none
c
      integer mthin,nthin,ndum1,ndum2,ngato,ith
      dimension xinf(*), zinf(*), delta(*), xigr(*),xigi(*)

      REAL, DIMENSION(:), ALLOCATABLE :: vecin

c-----------------------------------------------------------------------
c     read scalars.
c-----------------------------------------------------------------------
      OPEN(unit=3,file='vacin')
      WRITE ( nout1, '(/"reading VACIN",/)' )
      WRITE ( nout2, '(/"reading VACIN",/)' )
      read(3,'(//)')
      read(3,*)mthin
      read(3,*)rgato
      read(3,*)ndum2
      read(3,*)ngato
      read(3,*)qa1
      mthin1 = mthin + 1
c..   Adjust mth for finite elements, if lgato .eq. 1
      if ( lgato .eq. 1)
     $     call adjustm ( mth, mthin, mth1,mth2,  ndfel, nout1,nout2 )

      ALLOCATE ( vecin(mthin+5) )

c-----------------------------------------------------------------------
c     read arrays.
c-----------------------------------------------------------------------
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, xinf,mth, dx0,dx1 )
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, zinf,mth, dx0,dx1 )
      read(3,'(//)')
      read(3,*)(vecin(ith),ith=1,mthin1)
      call transdx ( vecin,mthin, delta,mth, dx,dx1 )
      if ( ireig .eq. 8 ) then
         read(3,'(//)')
         read(3,*)(xigr(ith),ith=1,mthin1)
         read(3,'(//)')
         read(3,*)(xigi(ith),ith=1,mthin1)
      end if
c Try changing sign of delta
c      do 100 i = 1, mth+2
c         delta(i) = -delta(i)
c 100  continue
      close(unit=3)

      DEALLOCATE ( vecin )

c-----------------------------------------------------------------------
c     terminate routine.
c-----------------------------------------------------------------------
      return
      end

c........................................................
      SUBROUTINE readnim ( mthin, n_nim, xinf,zinf, mth, dx0,dx1 )
c........................................................

c  Read input from NIMROD.
c         Toroidal mode number and [X,Z] for the resistive shell.
c         The shell will be treated as a plasma surface in "xinf,zinf".

c      implicit none
c
      INTEGER ::  mthin, n_nim
      REAL, DIMENSION(*) :: xinf, zinf

      REAL, DIMENSION(:), ALLOCATABLE :: zvec

c-----------------------------------------------------------------------
c     read scalars.
c-----------------------------------------------------------------------
      OPEN( UNIT=3, FILE = 'nim2msc.out')
      READ(3,*) mthin
      READ(3,*) n_nim

      mthin1 = mthin + 1
      mthdim = mthin + 5       ! need mthin+3 in lage4(transdx). But
                                ! use +10 to prevent overflow???

      ALLOCATE ( zvec(mthdim) )

c-----------------------------------------------------------------------
c     read arrays.
c-----------------------------------------------------------------------
      READ(3,'(//)')
      READ(3,*)( zvec(ith),ith=1,mthin1 )
      CALL transdx ( zvec,mthin, xinf,mth, dx0,dx1 )

      READ(3,'(//)')
      READ(3,*)( zvec(ith),ith=1,mthin1 )
      CALL transdx ( zvec,mthin, zinf,mth, dx0,dx1 )

      DEALLOCATE ( zvec )
     
      CLOSE( UNIT = 3)
c-----------------------------------------------------------------------
c     terminate routine.
c-----------------------------------------------------------------------

      END SUBROUTINE readnim


c.........................................................
c     
c      subroutine mdskrd0
cccc   
cccc   .............................................................
cccc   
c      include 'netcdf.inc'
c      include 'vacuum1.inc'
c      include 'vacuum2.inc'
c      include 'vacuum3.inc'
c      include 'vacuum4.inc'
cccc   
c      integer iargc, lena, counta
c      external iargc, getarg
c      character*72 namea, datype(6)
cccc   
cccc   ...... Store Data Types:
cccc   
c      datype(1) = "byte:      Eight-bit data. For saving space."
c      datype(2) = "character: Synonymous with byte. ASCII characters."
c      datype(3) = "short:     16-bit integers."
c      datype(4) = "long:      32-bit integers."
c      datype(5) = "float:     32-bit IEEE floating-point."
c      datype(6) = "double:    64-bit IEEE floating-point."
cccc   
cccc   ...... Open a NetCDF File for Access. 5.3
cccc   
c      counta = iargc()
c      if ( counta .ge. 1 ) then
c         call getarg ( 1, namea )
c         lena = index ( namea, ' ' ) - 1
c         cdfin = namea(1:lena)
c         cdfout = 'vac'//namea(4:lena)
cc         jobid = namea(4:lena)
c         jobid = namea(1:lena)
c         if ( nout(3) .gt. 0 ) then
c            open ( nout(3), file='vac'//namea(4:lena)//'.out',
c     $           status='unknown', form='formatted' )
c         end if
c      else
c         cdfin = "mapdsk.cdf"
c         cdfout = "vacdsk.cdf"
c         jobid = "mapdsk.cdf"
c         if ( nout(3) .gt. 0 ) then
c            open ( nout(3), file="vacdsk.cdf.out",
c     $           status='unknown', form='formatted' )
c         end if
c      end if
cccc   
c      call writg1 ( "Opening a NetCDF File (5.3).", seps, 3, nout ) 
cccc   
c      cdfid = ncopn ( cdfin, ncnowrit, rcode )
c      call nwopn ( cdfid, cdfin, "ncnowrit", rcode, nout )
cccc   
c      call writg1 ( "Inquiring about a NetCDF File (5.7).",
c     $     seps, 3, nout )
cccc   
cccc   ...... Inquire about an Open NetCDF file. 5.7
cccc   
c      call ncinq ( cdfid, ndims, nvars, ngatts, recid, rcode )
c      call nwinq ( cdfid, ndims, nvars, ngatts, recid, rcode, nout )
cccc   
c      call writg1 ( "Dimensions' Information (6.3).", seps, 1, nout )
cccc   
cccc   .....  Inquire about a Dimension. 6.3.  All ndims of them.
cccc   
c      do idims = 1, ndims
cccc   ..... get name and size, (even though we may already know name)
c         call ncdinq ( cdfid, idims, vname, dimsiz(idims), rcode )
c         call nwdinq ( cdfid, idims, vname, dimsiz(idims), rcode,
c     $        nout )
c      end do
cccc   
c      call writg1 ( "The Record Dimension info., if any.. ",
c     $     seps, 3, nout )
cccc   
cccc   .........    Already have the record dimension  ................
cccc   ..... get ID of record dimension (among other things)
cccc   call ncinq( cdfid, ndims, nvars, ngatts, recid, rcode )
cccc   ................................................................
cccc   ..... get record dimension name and current size. Should be none.
cccc   
c      if ( recid .ne. -1 ) then
c         call ncdinq ( cdfid, recid, recnam, nrecs, rcode )
c         call nwdinq ( cdfid, recid, recnam, nrecs, rcode, nout )
c      end if
cc     
cccc   ...............................CHECK1 . >>>>>>>>>>>>>>>>>>>>>>>>>
cc     
c      if ( check1 ) then
cccc   
c         call writg1 ( "The Global Atrtributes (8.3): ", seps, 3, nout )
cccc   
cccc   .....  Get the Global Attributes' values. 8.3
cccc   
c         do igatts = 1, ngatts
c            call ncanam ( cdfid, ncglobal, igatts, attnam, rcode )
cccc   call nwanam ( cdfid, ncglobal, igatts, attnam, rcode, nout )
c            call ncainq ( cdfid, ncglobal, attnam, attype, attlen,
c     $           rcode )
cccc   call nwainq ( cdfid, ncglobal, attnam, attype, attlen,
cccc   $        rcode, nout )
c            if ( attype .eq. 2 ) then
c               call ncagtc ( cdfid, ncglobal, attnam, astrng, maxa1,
c     $              rcode )
c               call nwagtc ( cdfid, ncglobal, attnam, astrng, attlen,
c     $              maxa1, rcode, nout )
c            else
c               call ncagt ( cdfid, ncglobal, attnam, attval, rcode )
c               call nwagt ( cdfid, ncglobal, attnam, attval, attlen,
c     $              rcode, nout )
c            end if
c            call writg1 ( ".............................", seps, 0,
c     $           nout )
c         end do
cccc   
c         call writg1 ( "Variable's Information (7.3): ", seps, 1, nout )
cccc   
c         call writg1 ( "        Data Types:", seps, 1, nout )
cccc   
c         WRITE ( 6, '(/,3x,"No.",15x,"Type",/, (i4, 2x, a72) )' )
c     $        ( idat, datype(idat), idat = 1, 6 )
c         WRITE ( outmod, '(/,3x,"No.",15x,"Type",/, (i4, 2x, a72) )' )
c     $        ( idat, datype(idat), idat = 1, 6 )
cccc   
c      end if
cc     
cc...  <<<<<<<<<<<<<<<<<<<<<<<...CHECK1.................................
cccc   
cccc   ... Get Information about a Variable from Its ID. 7.3. All of them.
cccc   
c      do ivars  = 1, nvars
c         call ncvinq ( cdfid, ivars, vname, vtype, vn(ivars),
c     $        vdims, vnatt, rcode )
cccc   call nwvinq ( cdfid, ivars, vname, vtype, vn(ivars),
cccc   $        vdims, vnatt, rcode, nout )
c         vrname(ivars) = vname
c         vrtype(ivars) = vtype
c         vrnat(ivars) = vnatt
c         if ( vn(ivars) .ne. 0 ) then
c            do ivn = 1, vn(ivars)
c               vvdims(ivars,ivn) = vdims(ivn)
c            end do
c         end if
cccc   
c      end do
cccc   
cc...........................CHECK1 >>>>>>>>>>>>>>>>>>>>>>>>>>>.
cc     
c      if ( check1 ) then
cccc   
c         WRITE ( 6, 1001 )
c         WRITE ( outmod, 1001 )
c 1001    format ( /,3x, "Id #", 2x, "Variable", 3x, "V-Type", 1x,
c     $        "Attributes" 1x, "Dimensions" )
c         do ivars = 1, nvars
c            write ( 6, 1002 ) ivars, vrname(ivars), vrtype(ivars),
c     $           vrnat(ivars),
c     $           ( dimsiz(vvdims(ivars,ivn)), ivn = 1, vn(ivars) )
c            write ( outmod, 1002 ) ivars, vrname(ivars), vrtype(ivars),
c     $           vrnat(ivars),
c     $           ( dimsiz(vvdims(ivars,ivn)), ivn = 1, vn(ivars) )
c         end do
c 1002    format ( 2x, i5, 2x, a10, i5,4x, i5,3x, 5i5 ) 
cccc   
cccc   
c      end if
cc     
cc     .<<<<<<<<<<<<<<<<<<<<<<<<<<<..CHECK1.............................
cccc   
c      call writg1 ( "Reading the Variables.", seps, 3, nout )
cccc   
cccc   
c      vid = ncvid ( cdfid, 'ctitle', rcode )
c      call nwvid ( cdfid, 'ctitle', vid, rcode, nout )
cccc   
cccc   ..... Read The Character String 'ctitle'. 7.7.
cccc   
c      start(1) = 1
c      count(1) = dimsiz( vvdims(vid,1) )
c      nd = 1
c      length = maxc1
c      call ncvgtc ( cdfid, vid, start,count, ctitle, length,
c     $     rcode )
cc     call nwvgtc ( cdfid, vid, start,count,nd, ctitle, length,
cc     $     rcode, nout )
cccc   
c      vid = ncvid ( cdfid, 'date0', rcode )
c      start(1) = 1
c      count(1) = dimsiz( vvdims(vid,1) )
c      nd = 1
c      length = maxct
c      call ncvgtc ( cdfid, vid, start,count, date0, length,
c     $     rcode )
cccc   
c      vid = ncvid ( cdfid, 'time0', rcode )
c      start(1) = 1
c      count(1) = dimsiz( vvdims(vid,1) )
c      nd = 1
c      length = maxct
c      call ncvgtc ( cdfid, vid, start,count, time0, length,
c     $     rcode )
cccc   
c      vid = ncvid ( cdfid, 'datem', rcode )
c      start(1) = 1
c      count(1) = dimsiz( vvdims(vid,1) )
c      nd = 1
c      length = maxct
c      call ncvgtc ( cdfid, vid, start,count, datem, length,
c     $     rcode )
cccc   
c      vid = ncvid ( cdfid, 'timem', rcode )
c      start(1) = 1
c      count(1) = dimsiz( vvdims(vid,1) )
c      nd = 1
c      length = maxct
c      call ncvgtc ( cdfid, vid, start,count, timem, length,
c     $     rcode )
cccc   
c      vid = ncvid ( cdfid, 'dskout', rcode )
c      start(1) = 1
c      count(1) = dimsiz( vvdims(vid,1) )
c      nd = 1
c      length = nccl3
c      call ncvgtc ( cdfid, vid, start,count, dskout, length,
c     $     rcode )
cccc   
cccc   
c      vid = ncvid ( cdfid, 'nx', rcode )
c      vindx(1) = 1
c      nd = vn(vid)
c      call ncvgt1 ( cdfid, vid, vindx, nx, rcode )
cccc   
c      vid = ncvid ( cdfid, 'nz', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, nz, rcode )
cccc   
c      vid = ncvid ( cdfid, 'nosurf', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, nosurf, rcode )
cccc   
c      vid = ncvid ( cdfid, 'mth', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, mthin, rcode )
cccc   
c      vid = ncvid ( cdfid, 'lpless', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, lpless, rcode )
cccc   
c      vid = ncvid ( cdfid, 'alx', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, alx, rcode )
cccc   
c      vid = ncvid ( cdfid, 'alz', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, alz, rcode )
cccc   
c      vid = ncvid ( cdfid, 'xzpst', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, xzero, rcode )
cccc   
c      vid = ncvid ( cdfid, 'xmag', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, xma, rcode )
cccc   
c      vid = ncvid ( cdfid, 'rpst', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, r, rcode )
cccc   
c      vid = ncvid ( cdfid, 'p0pst', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, p0, rcode )
cccc   
c      vid = ncvid ( cdfid, 'gppst', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, gp0, rcode )
cccc   
c      vid = ncvid ( cdfid, 'pminpst', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, psimin, rcode )
cccc   
c      vid = ncvid ( cdfid, 'plimpst', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, psilim, rcode )
cccc   
c      vid = ncvid ( cdfid, 'plspst', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, psipls, rcode )
cccc   
c      vid = ncvid ( cdfid, 'betapst', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, betag, rcode )
cccc   
c      vid = ncvid ( cdfid, 'betag', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, betap, rcode )
cccc   
c      vid = ncvid ( cdfid, 'rjacps', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, rjacps, rcode )
cccc   
c      if ( lpest1 ) then
c         vid = ncvid ( cdfid, 'upsiln', rcode )
c         call ncvgt1 ( cdfid, vid, vindx, upsiln, rcode )
c      end if
cccc   
c      vid = ncvid ( cdfid, 'mjacx', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, mj, rcode )
cccc   
c      vid = ncvid ( cdfid, 'njacg', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, nj, rcode )
cccc   
c      vid = ncvid ( cdfid, 'ljacb', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, lj, rcode )
cccc   
cc     vid = ncvid ( cdfid, '   ', rcode )
cc     call ncvgt1 ( cdfid, vid, vindx, nxy(4), rcode )
cccc   
cc     vid = ncvid ( cdfid, '   ', rcode )
cc     call ncvgt1 ( cdfid, vid, vindx, nxy(5), rcode )
cccc   
cc     vid = ncvid ( cdfid, '   ', rcode )
cc     call ncvgt1 ( cdfid, vid, vindx, nxy(6), rcode )
cccc   
cc     vid = ncvid ( cdfid, '   ', rcode )
cc     call ncvgt1 ( cdfid, vid, vindx, nxy(7), rcode )
cccc   
cc     vid = ncvid ( cdfid, '   ', rcode )
cc     call ncvgt1 ( cdfid, vid, vindx, nxy(8), rcode )
cccc   
c      vid = ncvid ( cdfid, 'nzd1', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, nzd1map, rcode )
cccc   
c      vid = ncvid ( cdfid, 'nzd2', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, nzd2map, rcode )
cccc   
c      vid = ncvid ( cdfid, 'dth', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, dtmap, rcode )
cccc   
c      vid = ncvid ( cdfid, 'dr', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, drmap, rcode )
cccc   
c      vid = ncvid ( cdfid, 'pi', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, pimap, rcode )
cccc   
c      vid = ncvid ( cdfid, 'xzero', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, xzeromap, rcode )
cccc   
c      vid = ncvid ( cdfid, 'zmag', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, zma, rcode )
cccc   
cc     vid = ncvid ( cdfid, 'qstar', rcode )
cc     call ncvgt1 ( cdfid, vid, vindx, qstar, rcode )
cccc   
cc     vid = ncvid ( cdfid, 'betav', rcode )
cc     call ncvgt1 ( cdfid, vid, vindx, betav, rcode )
cccc   
cc     vid = ncvid ( cdfid, 'betai', rcode )
cc     call ncvgt1 ( cdfid, vid, vindx, betai, rcode )
cccc   
cc     vid = ncvid ( cdfid, 'ctroy', rcode )
cc     call ncvgt1 ( cdfid, vid, vindx, ctroy, rcode )
cccc   
c      vid = ncvid ( cdfid, 'upsiln', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, upsiln, rcode )
cccc   
cc     vid = ncvid ( cdfid, '   ', rcode )
cc     call ncvgt1 ( cdfid, vid, vindx, axy(10), rcode )
cccc   
c      vid = ncvid ( cdfid, 'betat', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, betat, rcode )
cccc   
c      vid = ncvid ( cdfid, 'betap', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, betap, rcode )
cccc   
c      vid = ncvid ( cdfid, 'ctroy', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, ctroy, rcode )
cccc   
c      vid = ncvid ( cdfid, 'betatot', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, betatot, rcode )
cccc   
c      vid = ncvid ( cdfid, 'betapo', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, betapo, rcode )
cccc   
c      vid = ncvid ( cdfid, 'betato', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, betato, rcode )
cccc   
c      vid = ncvid ( cdfid, 'betats', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, betats, rcode )
cccc   
c      vid = ncvid ( cdfid, 'btor', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, btor, rcode )
cccc   
c      vid = ncvid ( cdfid, 'aima', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, aima, rcode )
cccc   
c      vid = ncvid ( cdfid, 'aioab', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, aioab, rcode )
cccc   
c      vid = ncvid ( cdfid, 'bt2dv', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, bt2dv, rcode )
cccc   
c      vid = ncvid ( cdfid, 'pdvamu0', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, pdvamu0, rcode )
cccc   
c      vid = ncvid ( cdfid, 'dv', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, dv, rcode )
cccc   
c      vid = ncvid ( cdfid, 'xliga', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, xliga, rcode )
cccc   
c      vid = ncvid ( cdfid, 'dlp', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, dlp, rcode )
cccc   
c      vid = ncvid ( cdfid, 'qstar', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, qstar, rcode )
cccc   
c      vid = ncvid ( cdfid, 'qratio', rcode )
c      call ncvgt1 ( cdfid, vid, vindx, qratio, rcode )
cccc   
c      vid = ncvid ( cdfid, 'q', rcode )
c      call ncvgt1 ( cdfid, vid, 1, q00, rcode )
cccc   
c      vid = ncvid ( cdfid, 'p', rcode )
c      call ncvgt1 ( cdfid, vid, 1, p01, rcode )
cccc   
c      vid = ncvid ( cdfid, 'fb', rcode )
c      call ncvgt1 ( cdfid, vid, 1, f0, rcode )
cccc   
c      vid = ncvid ( cdfid, 'fb', rcode )
c      call ncvgt1 ( cdfid, vid, nosurf, fa1, rcode )
cccc   
c      vid = ncvid ( cdfid, 'g', rcode )
c      call ncvgt1 ( cdfid, vid, 1, g0, rcode )
cccc   
c      vid = ncvid ( cdfid, 'g', rcode )
c      call ncvgt1 ( cdfid, vid, nosurf, ga1, rcode )
cccc   
c      vid = ncvid ( cdfid, 'q', rcode )
c      call ncvgt1 ( cdfid, vid, nosurf, qa1, rcode )
cccc   
cccc   
c      q0 = q00
c      q1 = qa1
c      ljacb = lj
c      mjacx = mj
c      njacg = nj
cc     
c      write ( outmod, 9000 ) ctitle, jobid
c      write ( iotty,  9000 ) ctitle, jobid
cc     
c 9000 format (/, 1x, a,/, 1x, " jobid = ", a )
cccc   
c      write(outmod,2002) datev,timev, datem,timem, date0,time0,
c     $     dskout, nosurf, mthin,
c     $     mjacx, njacg, ljacb,
c     $     mjacx, njacg, ljacb
cc$$$  write(96,2002) datev,timev, datem,timem, date0,time0,
cc$$$  $        dskout, nosurf, mthin,
cc$$$  $        mjacx, njacg, ljacb,
cc$$$  $        mjacx, njacg, ljacb
c      write( 6,2002) datev,timev, datem,timem, date0,time0,
c     $     dskout, nosurf, mthin,
c     $     mjacx, njacg, ljacb,
c     $     mjacx, njacg, ljacb
c 2002 format(/,'           VACUUM on ', A10, ' at ', a10,/,
c     $     '          MAPPING on ', a10, ' at ', a10,/,
c     $     ' From EQUILIBRIUM on ', a10, ' at ', a10,
c     $     ' to dsk: ', a ,/,
c     $     ' with NOSURF = ',i3,' MTHIN = ',i4,/,
c     $     '      MJACX, NJACG, LJACB =',3i2,/,
c     $     ' Jacobian =  X ** ',i1,
c     $     '/ ( Grad-PSI ** ',i1,' * B ** ',i1,')',/ )
cccc   
c      write(nout(2),3066) betap,betat,betatot
c      write(6,3066) betap,betat,betatot
c 3066 format(/," betap,betat,beta=",f6.2,f8.2,"%   ",f8.2,"% ")
c      write(outmod,3067) betapo,betato,betatot
c      write(6,3067) betapo,betato,betatot
c 3067 format(" betap,betat,beta(old def.)=",f6.2,f8.2,"%   ",f8.2,"% ")
c      write(outmod,3068) betats
c      write(6,3068) betats
c 3068 format(" beta star =",f8.2,"%")
cc     
c      write(outmod,3777) btor,aima,aioab,ctroy,q00,q1
c      write(6,3777) btor,aima,aioab,ctroy,q00,q1
c      write(outmod,3778) qstar,qratio
c      write(6,3778) qstar,qratio
c      write(outmod,3779) bt2dv,pdvamu0,dv
c      write(6,3779) bt2dv,pdvamu0,dv
c      write(outmod,3780) xliga,dlp
c      write(6,3780) xliga,dlp
c 3777 format(" tor field=  ",f6.2,"(T)     IP=",f7.3,"(MA) ",/,
c     1     " I(MA)/A(m)B(T)= ",f6.2,"   troyon factor",f6.2,/,
c     2     " q(axis)  =  ",f6.2,"        q(edge)",f6.2)
c 3778 format(" qstar    =  ",f6.2,"  qstar/q(1) = ",f6.2)
c 3779 format(" bt2dv= ",e12.4," pdv= ",e12.4," dv= ",e12.4)
c 3780 format(" li(GA)= ",e12.4," dlp= ",e12.4)
cccc   
cc     
c      write ( outmod,9600 ) qa1, fa1, ga1
c      write ( iotty, 9600 ) qa1, fa1, ga1
c 9600 format ( /, 1x, " qa1, fa1, ga1 = ", 1p3e12.5,/ )
c      write ( outmod,9500 ) r, upsiln, xma
c      write ( iotty, 9500 ) r, upsiln, xma
c 9500 format ( /, 1x, "r, upsiln, xma  = ", 1p3e12.5,/ )
cccc   
c      mthin1 = mthin  + 1
c      mthin2 = mthin1 + 1
cccc   
c      r2     = r * r
c      r4     = r2 * r2
c      r6     = r4 * r2
cccc   
c      return
c      end
cc     
cc................................................................
c      subroutine mdskrd1
cc... ............................................................
cc
c      include 'netcdf.inc'
c      include 'vacuum1.inc'
c      include 'vacuum2.inc'
c      include 'vacuum3.inc'
c      include 'vacuum4.inc'
cc     
c      dimension vecin(ntsin)
cccc   
cc     ..... iref is index to set reference values of g and gp in
cc     .....  subroutine varymo  each time the mapping is read.
cc     
c      iref = 0
cc     
c      mthin1 = mthin + 1
c      mthin2 = mthin1 + 1
cc     
cc     Start reading from netCDF file.
cc     
cccc   ....... 2D  Arrays......
cc     
c      vid = ncvid ( cdfid, 'x', rcode )
cccc   
c      nd = 2
c      start(1) = 1
c      start(2) = nosurf
c      count(1) = dimsiz ( vvdims(vid,1) )
c      count(2) = 1
cccc   
c      call ncvgt ( cdfid, vid, start, count, vecin, rcode )
c      call nwvgt ( cdfid, vid, start,count,nd, vecin,vals, nd1,nd2,
c     $     rcode, nout0 )
c      call trans ( vecin,mthin, xinf,mth )
cccc   
c      vid = ncvid ( cdfid, 'z', rcode )
c      call ncvgt ( cdfid, vid, start, count, vecin, rcode )
c      call nwvgt ( cdfid, vid, start,count,nd, vecin,vals, nd1,nd2,
c     $     rcode, nout0 )
c      call trans ( vecin,mthin, zinf,mth )
cccc   
c      vid = ncvid ( cdfid, 'grpssq', rcode )
c      call ncvgt ( cdfid, vid, start, count, vecin, rcode )
c      call nwvgt ( cdfid, vid, start,count,nd, vecin,vals, nd1,nd2,
c     $     rcode, nout0 )
c      call trans ( vecin,mthin, grpssq,mth )
cccc   
c      if ( .not. lpest1 ) then
c         vid = ncvid ( cdfid, 'xjacob', rcode )
c         call ncvgt ( cdfid, vid, start, count, vecin, rcode )
c         call nwvgt ( cdfid, vid, start,count,nd, vecin,vals, nd1,nd2,
c     $        rcode, nout0 )
c         call trans ( vecin,mthin, xjacob,mth )
cccc   
c         vid = ncvid ( cdfid, 'delta', rcode )
c         call ncvgt ( cdfid, vid, start, count, vecin, rcode )
c         call nwvgt ( cdfid, vid, start,count,nd, vecin,vals, nd1,nd2,
c     $        rcode, nout0 )
c         call trans ( vecin,mthin, delta,mth )
c      end if
cc     
c      if ( lpest1 ) then
c         do i = 1, mth2
c            xjacob(i) = upsiln * xinf(i)**2 / ( twopi*r )
c         end do
c      end if
cc     
cccc   
c      if ( lnova ) then
c         vid = ncvid ( cdfid, 'xsq', rcode )
c         call ncvgt ( cdfid, vid, start, count, vecin, rcode )
c         call nwvgt ( cdfid, vid, start,count,nd, vecin,vals, nd1,nd2,
c     $        rcode, nout0 )
c         call trans ( vecin,mthin, xsq,mth )
cccc   
c         vid = ncvid ( cdfid, 'gpsdth', rcode )
c         call ncvgt ( cdfid, vid, start, count, vecin, rcode )
c         call nwvgt ( cdfid, vid, start,count,nd, vecin,vals, nd1,nd2,
c     $        rcode, nout0 )
c         call trans ( vecin,mthin, gpsdth,mth )
cccc   
c         vid = ncvid ( cdfid, 'xsqdth', rcode )
c         call ncvgt ( cdfid, vid, start, count, vecin, rcode )
c         call nwvgt ( cdfid, vid, start,count,nd, vecin,vals, nd1,nd2,
c     $        rcode, nout0 )
c         call trans ( vecin,mthin, xsqdth,mth )
cccc   
c         vid = ncvid ( cdfid, 'xjacob', rcode )
c         call ncvgt ( cdfid, vid, start, count, vecin, rcode )
c         call nwvgt ( cdfid, vid, start,count,nd, vecin,vals, nd1,nd2,
c     $        rcode, nout0 )
c         call trans ( vecin,mthin, xjacob,mth )
cc     
c         mthd2p1=mth/2+1
cc     
c         do i = 1, mthd2p1
c            zbsq = ( grpssq(i) + r2*ga1 ) / xsq(i)
c            zdbsqb = ( gpsdth(i) - zbsq*xsqdth(i) ) / (xsq(i)*zbsq )
c            xsdtxs(i) = xsqdth(i)/xsq(i)
c            gpdtgp(i) = gpsdth(i)/grpssq(i)
c            xjdtxj(i) = 0.5 * ( mj*xsdtxs(i)-nj*gpdtgp(i) - lj*zdbsqb )
cccc   even parity
c            xjacob(mth2-i)=xjacob(i)
cccc   odd parity
c            xjdtxj(mth2-i)=-xjdtxj(i)
c            delta(mth2-i)=-delta(i)
c         end do
cc     
c         xjdtxj(1)=0.
c         xjdtxj(mthd2p1)=0.
c         delta(1)=0.
c         delta(mthd2p1)=0.
cc     
c      end if
cccc   
c      return
c      end
c     
c     ..................................................................
      subroutine adjustm ( mth, mfel, mth1,mth2, ndfel, nout1,nout2 )
c     ..................................................................
c..   Make sure that ndfel, the number of grid points spanning a 
c     finite elment is even, unless its one.
c 
      if ( mth .eq. mfel ) then
         ndfel = 1
         return
      end if
c
      mth00 = mth
      ndfel = mth / mfel
      if ( ndfel .lt. 1 ) ndfel = 1
      ndfel = ( (ndfel+1)/2 ) * 2
      mth = ndfel * mfel
c     
      if ( mth00 .ne. mth ) then
         write ( nout1,
     $        '(/,1x,"****** WARNING: mth00 .ne. mth ******",/ )' )
         write ( nout2,
     $        '(/,1x,"****** WARNING: mth00 .ne. mth ******",/ )' )
         mth1 = mth + 1
         mth2 = mth1 + 1
      end if
c     
      write ( nout1, '(/,5x, "mth00, mth, mfel, ndfel = ", 4i5 )' )
     $     mth00, mth, mfel, ndfel
      write ( nout2, '(/,5x, "mth00, mth, mfel, ndfel = ", 4i5 )' )
     $     mth00, mth, mfel, ndfel
c     
      return
      end
c.....................................................................














