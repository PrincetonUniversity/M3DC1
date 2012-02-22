c....................................................
      program vacuum
c....................................................
c
      include 'vacuum1.inc'
      include 'vacuum2.inc'
      include 'vacuum3.inc'
      include 'vacuum5.inc'
      include 'vacuum7.inc'
      include 'vacuum8.inc'
c
      REAL, DIMENSION(:), ALLOCATABLE :: cgrd, cgrda, cgrdb
      REAL, DIMENSION(:), ALLOCATABLE :: shfcoil, shfcoila, shfcoilb

      integer iargc, lena, counta
      external iargc, getarg
      character*72 namea, name2
c      character*20 filein2
c      CHARACTER*24 date_array
      character*8 cgmname
      LOGICAL backl

      CHARACTER*(80) stringv(5), string1
c
      cgmname = 'vac.cgm'
 
      irepeat = -1 ! initialize index for repeating the code at statement 7777

      IF ( lreshel /= 0 ) THEN
         iosen = 28
         OPEN (  iosen, FILE='SENSOR-MEAS', STATUS='REPLACE',
     $        FORM='FORMATTED' )
         WRITE ( iosen, '(/, 10x, "IOSEN = ", I5 )' ) iosen
      END IF
c
      call defglo
c
      call shella
c
c......initialize graphics, etc
c
      counta = iargc()
      backl = .false.

      outmod = 23
      WRITE ( IOTTY,  '(/, 10x, "OUTMOD = ", I5 )' ) outmod

      OPEN (  UNIT=outmod, FILE='modovmc', STATUS='new',
     $     FORM='formatted' )

      WRITE ( OUTMOD, '(/, 10x, "OUTMOD = ", I5 )' ) outmod

      WRITE (  IOTTY, '(/, 10x, "No. OF ARGUMENTS = ", I20 )' )
     $     counta
      WRITE (  OUTMOD, '(/, 10x, "NO. OF ARGUMENTS = ", I20 )' )
     $     counta

      IF ( counta .GE. 1 ) THEN
         CALL GETARG ( 1, namea )
         lena = LEN_TRIM ( namea )
c         lena = index ( namea, '  ', backl ) - 1
         jobid = namea(1:lena)
      ELSE
         jobid = "No Jobid"
      END IF

      WRITE (  IOTTY, '(/, 10x, "JOBID = ", a, )' )
     $     jobid
      WRITE (  OUTMOD, '(/, 10x, "JOBID = ", a, )' )
     $     jobid
      IF ( lreshel /= 0 ) WRITE (  IOSEN, '(/, 10x, "JOBID = ", a, )' )
     $     jobid
      
      IF ( counta .ge. 2 ) THEN
         CALL GETARG ( 2, name2 )
c         lena = index ( name2, '  ', backl ) - 1 
         lena = LEN_TRIM ( name2 )
         lenif = lena
         filein2 = name2(1:lena)
      ELSE
         lenif = 7
         filein2  = "modivmc"
      END IF

      WRITE (  IOTTY, '(/, 10x, "INPUT FILE = ", a, / )' )
     $     filein2
      WRITE (  OUTMOD, '(/, 10x, "INPUT FILE = ", a, / )' )
     $     filein2
      IF ( lreshel /= 0 )
     $     WRITE (  IOSEN, '(/, 10x, "INPUT FILE = ", a, / )' )
     $     filein2

      open ( outpest, file='pestotv', status='new', form='formatted' )
      open (  inmode, file=filein2, status='old', form='formatted' )

      OPEN ( UNIT = ioshel, FILE = "SHELL_DATA", STATUS = "REPLACE" )

      call date_time ( date_array, datev, timev, iotty,outmod )

c$$$      write (  iotty, '("date and time = ", a10,1x,a10)' )
c$$$     $     datev, timev
c$$$      write ( outmod, '("date and time = ", a10,1x,a10)' )
c$$$     $     datev, timev

      write ( iotty,
     $   '(/,"Parameters: nsf0, nths0, ntsin0, nths,   nfm  = ",/,
     $     9x,5I7 )' ) nsf0, nths0, ntsin0, nths, nfm

      write ( outmod,
     $   '(/,"Parameters: nsf0, nths0, ntsin0, nths,   nfm  = ",/,
     $     9x,5I7 )' ) nsf0, nths0, ntsin0, nths, nfm

      write ( iotty,
     $   '("Parameters: npfg11, npfg12, npfg13, npfg14, npfg2x  = "/,
     $     8x,5I8 )' ) npfg11, npfg12, npfg13, npfg14, npfg2x

      write ( outmod,
     $   '("Parameters: npfg11, npfg12, npfg13, npfg14, npfg2x  = ",/,
     $     8x,5I8 )' ) npfg11, npfg12, npfg13, npfg14, npfg2x

      write ( iotty,
     $   '("Parameters: nptg11, nptg12, nptg13, nptg14, nptg2x  = ",/,
     $     8x,5I8 )' ) nptg11, nptg12, nptg13, nptg14, nptg2x

      write ( outmod,
     $   '("Parameters: nptg11, nptg12, nptg13, nptg14, nptg2x  = ",/,
     $     8x,5I8 )' ) nptg11, nptg12, nptg13, nptg14, nptg2x
c
      i_one = 1
      call ncarcgm ( i_one, cgmname)
c         CALL GSCR(1, 0, 1.0, 1.0, 1.0)
c         CALL GSCR(1, 61, 0.0, 0.0, 0.0)
c
      CALL setlch ( 0.2, 0.8, 0, 2, 0, -1 )

      WRITE ( stringv, '(1x, "date and time = ", a )' ) date_array
      CALL wrtstr ( stringv, 1 ) 

      WRITE ( stringv, '(/, 1x, "JOBID = ", a, / )' )  jobid
      CALL wrtstr ( stringv, 2 )

      WRITE ( stringv, '(/, 1x, "INPUT FILE = ", a, / )' ) filein2
      CALL wrtstr ( stringv,2 )

      CALL framep(jobid, ff)

      call timer ( outmod, iotty, "top of main" )
c
c**********************************************************************
c         entry point for segment 33.
c**********************************************************************
c
c      input data from cards.
c
      call cardmo
c
      call inglo
c
c      first read from disk.

c...&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c.. This is the return point for repeating VACUUM for ITER Diagnostics

 7777 CONTINUE
ccc&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      WRITE ( IOTTY,  '(/, 5x, "irepeat = ", I3 )' ) irepeat
      WRITE ( OUTMOD, '(/, 5x, "irepeat = ", I3 )' ) irepeat

      irepeat = irepeat + 1

      call dskmd1
c
c      set up equilibrium functions
c
      IF ( irepeat == 0 ) then
         CALL funint
         CALL varymo
      END IF

c... Flags on input data:

      WRITE ( iotty,
     $     '(/,15x,"@@@@@@ Flags On Input Data @@@@@@",/)' )
      WRITE ( outmod,
     $     '(/,15x,"@@@@@@ Flags On Input Data @@@@@@",/)' )
c
      IF ( lgato == 1 .OR. lgato == 3 ) THEN
         IF ( (lkdis .gt. 0) .OR. (iloop .gt. 0) ) THEN
            PRINT '(/, a,i2,a,i2,a)', "  lkdis = ",lkdis,
     $           " iloop = ",iloop, " <-- ONLY ALLOWED WITH LGATO = 2"
            WRITE (OUTMOD, '(/, a,i2,a,i2,a)') "  lkdis = ",lkdis,
     $           " iloop = ",iloop, " <-- ONLY ALLOWED WITH LGATO = 2"
         END IF
         PRINT '(/,a)', "   SETTING lkdis = 0,  iloop = 0"
         WRITE ( OUTMOD, '(/,a)') "   SETTING lkdis = 0,  iloop = 0"
         lkdis = 0
         iloop = 0
      END IF
c
c... Check on LTFIL option:
      IF ( (lmin(1) < - mfel/2) .OR. (lmax(1) > mfel/2 )
     $     .AND. (ltfil > 0) ) THEN
         ltfil = 0
         WRITE ( iotty, '
     $        ("Inconsistency with ltfil, lmin, lmax and mfel",/,
     $        "Setting ltfil = 0" )' )
         WRITE ( outmod, '
     $        ("Inconsistency with ltfil, lmin, lmax and mfel",/,
     $        "Setting ltfil = 0" )')
      END IF

      WRITE ( iotty,
     $     '(/,15x,"^^^^^ Flags On Input Data ^^^^^^",/)' )
      WRITE ( outmod,
     $     '(/,15x,"^^^^^ Flags On Input Data ^^^^^^",/)' )

c      if ( wall ) go to 10

c...Main Vacuum Calculations here.

      WRITE ( outmod, '("**** Calling SUBROUTINE wwall in MAIN here" )')

      if ( lfarw .eq. 0 ) call wwall ( mth1, xwal, zwal )
c
c....Set up arrays. First call.
c
      call arrays ( 1 )

c....Plot Field Line at Surface.

      nblines = 51
      CALL b_lines ( delta, qa1, mth1, nblines, jobid )

      if ( lfarw .eq. 0 ) call plgeom
      if ( lfarw .eq. 1 ) call plpgeom
c
      ALLOCATE ( cgrd(nths), cgrda(nths), cgrdb(nths) )

      DO i = 1, mtcoil1
         cgrd(i) = ( 1.0/mtcoil ) * (i-1)
      END DO
      DO i = 1, mtcoila1
         cgrda(i) = ( 1.0/mtcoila ) * (i-1)
      END DO
      DO i = 1, mtcoilb1
         cgrdb(i) = ( 1.0/mtcoilb ) * (i-1)
      END DO

c... Shfcoil is the coil function in shifted coordinates

      ALLOCATE ( shfcoil(nths), shfcoila(nths), shfcoilb(nths) )


      IF ( lfbcoila .EQ. 1 ) THEN
         CALL fbcoil ( coilanga,taucoila,coilcena, mtcoila, ncgr,
     $        iopfca, icshpa, acoila,bcoila,dcoila,drcoila,dzcoila,
     $        xcwa,zcwa, ccoila, cgrda, chgta )

         wcentr=0.
         CALL drawc3(xwal,zwal,xinf,zinf,xcwa,zcwa,1,mth1,1,mtcoila1,
     $        "z","x",xmaj,xma,zma,
     $        plrad,1,1,idot, jobid, ishape,a,b,aw,bw,cw,dw,tw,
     $        abulg, bbulg, tbulg, wcentr, pcircm,wcircm,
     $        pelong, pdness, wmaj,wlrad, welong, wdness, chgta )

         CALL shfcoilx ( xcwa, zcwa, ccoila, cgrda, mtcoila,
     $        zork1, zork2, ellc, shfcoila )

         zork4 = 0.0
         CALL kcur0 ( zork1,zork2,ellc,cgrda,mtcoila1,shfcoila,zork4,n,
     $        ncgr,wkgrd,wxgrd,wzgrd,wkdxt,wkdzt,
     $        wkt2,wkp2,nths,0,"Coil Current",ncgr, outmod )

c......Get derivative on coil points.
          CALL difsplt ( mtcoila, cgrda, xcwa, xcwap )
         CALL difsplt ( mtcoila, cgrda, zcwa, zcwap )
         xcwap(mtcoila+1) = xcwap(1)
         zcwap(mtcoila+1) = zcwap(1)
         xcwap(mtcoila+2) = xcwap(2)
         zcwap(mtcoila+2) = zcwap(2)
         xcwap = xcwap / twopi
         zcwap = zcwap / twopi
      END IF

      IF ( lfbcoila .EQ. 2 ) THEN
         CALL fbcoil2 ( ctip_l, ntip_l, xcwa,zcwa, xcwap, zcwap,
     $        ccoila, cgrda, clena_l, clen_l, chgta, "I-Coil-L" )
c$$$         xcwap = xcwap / twopi
c$$$         zcwap = zcwap / twopi

         wcentr=0.
c$$$         CALL drawc3(xwal,zwal,xinf,zinf,xcwa,zcwa,1,mth1,1,mtcoila1,
c$$$     $        "z","x",xmaj,xma,zma,
c$$$     $        plrad,1,1,idot, jobid, ishape,a,b,aw,bw,cw,dw,tw,
c$$$     $        abulg, bbulg, tbulg, wcentr, pcircm,wcircm,
c$$$     $        pelong, pdness, wmaj,wlrad, welong, wdness, chgta )
c$$$
c$$$         CALL shfcoilx ( xcwa, zcwa, ccoila, cgrda, mtcoila,
c$$$     $        zork1, zork2, ellc, shfcoila )
c$$$
c$$$         zork4 = 0.0
c$$$         CALL kcur0 ( zork1,zork2,ellc,cgrda,mtcoila1,shfcoila,zork4,n,
c$$$     $        ncgr,wkgrd,wxgrd,wzgrd,wkdxt,wkdzt,
c$$$     $        wkt2,wkp2,nths,0,"Coil Current",ncgr, outmod )
      END IF

      IF ( lfbcoilb .EQ. 1 ) THEN
         CALL fbcoil ( coilangb,taucoilb,coilcenb, mtcoilb, ncgr,
     $        iopfcb, icshpb, acoilb,bcoilb,dcoilb,drcoilb,dzcoilb,
     $        xcwb,zcwb, ccoilb, cgrdb, chgtb )

         CALL drawc3(xwal,zwal,xinf,zinf,xcwb,zcwb,1,mth1,1,mtcoilb1,
     $        "z","x",xmaj,xma,zma,
     $        plrad,1,1,idot, jobid, ishape,a,b,aw,bw,cw,dw,tw,
     $        abulg, bbulg, tbulg, wcentr, pcircm,wcircm,
     $        pelong, pdness, wmaj,wlrad, welong, wdness, chgtb )

         CALL shfcoilx ( xcwb, zcwb, ccoilb, cgrdb, mtcoilb,
     $        zork1, zork2, ellc, shfcoilb)

         zork4 = 0.0
         CALL kcur0 ( zork1,zork2,ellc,cgrdb,mtcoilb1,shfcoilb,zork4,n,
     $        ncgr,wkgrd,wxgrd,wzgrd,wkdxt,wkdzt,
     $        wkt2,wkp2,nths,0,"Coil Current",ncgr, outmod )

c......Get derivative on coil points.
         CALL difsplt ( mtcoilb, cgrdb, xcwb, xcwbp )
         CALL difsplt ( mtcoilb, cgrdb, zcwb, zcwbp )
         xcwbp(mtcoilb+1) = xcwbp(1)
         zcwbp(mtcoilb+1) = zcwbp(1)
         xcwbp(mtcoilb+2) = xcwbp(2)
         zcwbp(mtcoilb+2) = zcwbp(2)
         xcwbp = xcwbp / twopi
         zcwbp = zcwbp / twopi
      END IF

      IF ( lfbcoilb .EQ. 2 ) THEN
         CALL fbcoil2 ( ctip_u, ntip_u, xcwb,zcwb, xcwbp, zcwbp,
     $        ccoilb, cgrdb, clena_u, clen_u, chgtb, "I-coil-U" )
c$$$         xcwbp = xcwbp / twopi
c$$$         zcwbp = zcwbp / twopi

         wcentr=0.
c$$$         CALL drawc3(xwal,zwal,xinf,zinf,xcwb,zcwb,1,mth1,1,mtcoilb1,
c$$$     $        "z","x",xmaj,xma,zma,
c$$$     $        plrad,1,1,idot, jobid, ishape,a,b,aw,bw,cw,dw,tw,
c$$$     $        abulg, bbulg, tbulg, wcentr, pcircm,wcircm,
c$$$     $        pelong, pdness, wmaj,wlrad, welong, wdness, chgta )
c$$$
c$$$         CALL shfcoilx ( xcwa, zcwa, ccoila, cgrda, mtcoila,
c$$$     $        zork1, zork2, ellc, shfcoila )
c$$$
c$$$         zork4 = 0.0
c$$$         CALL kcur0 ( zork1,zork2,ellc,cgrda,mtcoila1,shfcoila,zork4,n,
c$$$     $        ncgr,wkgrd,wxgrd,wzgrd,wkdxt,wkdzt,
c$$$     $        wkt2,wkp2,nths,0,"Coil Current",ncgr, outmod )
      END IF

c.  Plot both coila and coilb together with a dclphse phase shift:

      IF ( lfbcoila == 1 .AND. lfbcoilb == 1 ) THEN

         clphsea =  2.0 * pye / 3.0
         clphseb = -2.0 * pye / 3.0

         zork3(1:mtcoila1) = shfcoila(1:mtcoila1) * cos(clphsea)  
     $        + shfcoilb(1:mtcoila1) * cos(clphseb)
         zork4(1:mtcoila1) = - shfcoila(1:mtcoila1) * sin(clphsea)
     $        - shfcoilb(1:mtcoila1) * sin(clphseb)

         CALL kcur0 ( zork1,zork2,ellc,cgrdb,mtcoilb1,zork3,zork4,n,
     $        ncgr,wkgrd,wxgrd,wzgrd,wkdxt,wkdzt,
     $        wkt2,wkp2,nths,0,"Coil Current",ncgr, outmod )

      END IF
      
      IF ( lfbcoil .EQ. 1 ) THEN

         CALL fbcoil ( coilang, taucoil, coilcen, mtcoil, ncgr,
     $        iopfc, icshp, acoil, bcoil, dcoil, drcoil,dzcoil,
     $        xcw,zcw, ccoil, cgrd,chgt )

c......Get derivative on coil points.
         CALL difsplt ( mtcoil, cgrd, xcw, xcwp )
         CALL difsplt ( mtcoil, cgrd, zcw, zcwp )
         xcwp(mtcoil+1) = xcwp(1)
         zcwp(mtcoil+1) = zcwp(1)
         xcwp(mtcoil+2) = xcwp(2)
         zcwp(mtcoil+2) = zcwp(2)
         xcwp = xcwp / twopi
         zcwp = zcwp / twopi

         CALL drawc3(xwal,zwal,xinf,zinf,xcw,zcw,1,mth1,1,mtcoil1,
     $        "z","x",xmaj,xma,zma,
     $        plrad,1,1,idot, jobid, ishape,a,b,aw,bw,cw,dw,tw,
     $        abulg, bbulg, tbulg, wcentr, pcircm,wcircm,
     $        pelong, pdness, wmaj,wlrad, welong, wdness, chgt )

         CALL shfcoilx ( xcw, zcw, ccoil, cgrd, mtcoil,
     $        zork1, zork2, ellc, shfcoil)

         zork4 = 0.0
         CALL kcur0 ( zork1,zork2,ellc,cgrd,mtcoil1,shfcoil,zork4,n,
     $        ncgr,wkgrd,wxgrd,wzgrd,wkdxt,wkdzt,
     $        wkt2,wkp2,nths,0,"Coil Current",ncgr, outmod )

      END IF
      IF ( lfbcoil .EQ. 2 ) THEN
         CALL fbcoil2 ( ctip_c, ntip_c, xcw,zcw, xcwp, zcwp,
     $        ccoil, cgrd, clena_c, clen_c, chgt, "C-Coil" )
c$$$         xcwp = xcwp / twopi
c$$$         zcwp = zcwp / twopi

         wcentr=0.
c$$$         CALL drawc3(xwal,zwal,xinf,zinf,xcwb,zcwb,1,mth1,1,mtcoilb1,
c$$$     $        "z","x",xmaj,xma,zma,
c$$$     $        plrad,1,1,idot, jobid, ishape,a,b,aw,bw,cw,dw,tw,
c$$$     $        abulg, bbulg, tbulg, wcentr, pcircm,wcircm,
c$$$     $        pelong, pdness, wmaj,wlrad, welong, wdness, chgta )
c$$$
c$$$         CALL shfcoilx ( xcwa, zcwa, ccoila, cgrda, mtcoila,
c$$$     $        zork1, zork2, ellc, shfcoila )
c$$$
c$$$         zork4 = 0.0
c$$$         CALL kcur0 ( zork1,zork2,ellc,cgrda,mtcoila1,shfcoila,zork4,n,
c$$$     $        ncgr,wkgrd,wxgrd,wzgrd,wkdxt,wkdzt,
c$$$     $        wkt2,wkp2,nths,0,"Coil Current",ncgr, outmod )
      END IF

      IF ( lfbcoil > 0 )
     $     CALL pospl1 ( pigrd,ccoil,  1,mtcoil, 1,"CCOIL", mtcoil1 )
      IF ( lfbcoila > 0 )
     $     CALL pospl1 ( pigrd,ccoila,  1,mtcoila, 2,"CCOIL-A",
     $     mtcoila1 )
      IF ( lfbcoilb > 0 )
     $     CALL pospl1 ( pigrd,ccoilb,  1,mtcoilb, 3,"CCOIL-B",
     $     mtcoilb1 )

      CALL framep(jobid,ff)

 21   CONTINUE
 
      DEALLOCATE ( cgrd, cgrda, cgrdb )
      DEALLOCATE ( shfcoil, shfcoila, shfcoilb )

c
      call vaccal
c
c....Update arrays. Second call.
c
      call arrays (2)
c
c   10 continue
c
      If ( lspark .ne. 0 ) call testvec
c
c--------------------------------------------------------------
c     Diagnostic Section
c--------------------------------------------------------------
c
      if ( ieig .eq. 0 ) go to 100
c
c....Update arrays again if LSPARK is on. Third call.
c
      if ( lspark .ne. 0 ) call arrays ( 0 )
c
      jmax1 = lrnge
c
      go to ( 41,42,43,44,45,46,47,48,49,50 ), ieig
c
 41   continue
c
c... PEST type input
c
      call efun ( xilr, xili, omsq )
c
      call makebp ( "PESTs EFUN" )
c
c..Calculate vacuum energy using PEST's perturbations
      call venergy ( "PEST-41" )
c
      GO TO 500
C
 42   CONTINUE

c...  Read B grid data.

      CALL makebl0 ( "READVACIN6C" )
c
      WRITE ( OUTMOD, '(//, "Had read B_i from VACIN61" )' )
      WRITE ( IOTTY,  '(//, "Had read B_i from VACIN61" )' )

      GO TO 500
C
 43   CONTINUE
      GO TO 500
C
 44   CONTINUE
C
c... Read input from file OUTIDST
c
      open ( 60, file='outidst', status='old', form='formatted' )
c
 25   read ( 60, '(i5)' ) mflag
      if ( mflag .ne. 0 ) go to 25
      read ( 60, '( e12.5,/,(8e12.5) )' ) omsq, ( xilr(l),l = 1,jmax1 )
c
      GO TO 500
C
 45   CONTINUE

c... Reads input from VACIN5, DCON like inputs. Reads Q in Fourier space.
c
      WRITE ( OUTMOD, '(//, "Had read B_l from VACIN5" )' )
      WRITE ( IOTTY,  '(//, "Had read B_l from VACIN5" )' )
      
      CALL make_bltobp ( "VACIN5 - DOCN-like" )

      GO TO 500

 46   CONTINUE
c

c  Reads DCON output

      CALL dcnvac

c$$$      IF (irotplasma.eq.0) THEN
c$$$         CALL resshel_2
c$$$      ELSE  IF (irotplasma.eq.1) THEN
c$$$c         CALL resshel_3
c$$$      END IF

      IF (irotplasma==0 .AND. lreshel==1) THEN
         CALL resshel_2
      ELSE   IF (irotplasma==1 .AND. lreshel==1) THEN
c         CALL resshel_3
      END IF

      CALL makebp ( "DCONs Output" )
c
      GO TO 500
C
 47   CONTINUE

c... Reads input from VACIN7, NOVA inputs. 
c    ieig = 7 if lreig7 == 1
c    Read wall parameters if lrwall == 1
c
      WRITE ( OUTMOD, '(//, "Had read B_l from VACIN7" )' )
      WRITE ( IOTTY,  '(//, "Had read B_l from VACIN7" )' )
      
      CALL make_bltobp ( "VACIN7 - DOCN-like" )

      GO TO 500
C
 48   CONTINUE
c
c...Reads from GATO's file VACIN.
c
      WRITE ( OUTMOD, '(//, "Had read xi from GATOs vacin" )' )
      WRITE ( IOTTY,  '(//, "Had read xi from GATOs vacin" )' )
c
      call makebg ( "GATO-48" )
 49   continue
c
c...Reads from GATO's file VACIN in READVACIN2, gets B directly.
c
      write ( outmod, '(//, "Had read B from GATOs vacin" )' )
c
c.. Calculate the vacuum energy using GATO's perturbations.
      if ( lgato .gt. 0 ) call venergy ( "GATO-49-I" )
c
c.. Now interpolate and shift the Magnetic Perturbations from READVACIN2
c
      dx00 = 0.5
      call transdxx ( bnkr, mthin, zork1,mth, dx00 )
      call transdxx ( bnki, mthin, zork1,mth, dx00 )
      write ( iotty,  '(//, "Had read B from GATOs vacin" )' )
c
      call makebl ( "GATO-49" )
c
c..Calculate vacuum energy using PEST's perturbations
      if ( lgato .eq. 0 ) call venergy ( "GATO-49_II" )
c
      go to 500

 50   CONTINUE
c
c... Reads input from namelist variables xiinr(l), xiini(l).
c
      DO j1 = 1, jmax1
         xilr(j1) = xiinr(j1)
         xili(j1) = xiini(j1)
      END DO
c
      IF (irotplasma==0 .AND. lreshel==1) THEN
         CALL timer ( outmod, iotty, "Before resshel_2" )
         CALL resshel_2
         CALL timer ( outmod, iotty, "After resshel_2" )
      ELSE   IF (irotplasma==1 .AND. lreshel==1) THEN
c         CALL resshel_3
      END IF
c
      CALL makebp ( "XIINr,i-50" )
c
c..Calculate vacuum energy using PEST's perturbations
      call venergy ( "PEST-45" )
c
      GO TO 500
c
 500   CONTINUE

      CALL diaplt

      IF ( lfarw == 0 ) CALL svdinpts
c
      if ( lreshel .gt. 0 ) call cplt ( bnlr, bnli )
      if ( lkdis .ne. 0 ) call kdis ( omsq, bnlr, bnli )
      if ( iloop .ne. 0 ) call pickup ( bnlr, bnli )
c

 60   CONTINUE
C
 100  CONTINUE

c... &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c.. Repeat VACUUM from here.
      IF ( irepeat < nrepeat - 1 ) GO TO 7777

c.. UNIT 661 is VACIN6_1 read in READVACIN6C. 
      CLOSE ( UNIT = 661 )
cc &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

c...CLOSE FILES

      CALL plote( )
c
      if ( lzio .eq. 1 ) then
         call zcl ( outmap1, 999 )
         call zcl ( iomode, 999 )
      end if
c
      call timer ( outmod, iotty, "end of main" )
c
      CLOSE ( UNIT = ioshel )
      CLOSE ( UNIT = outmod )
      CLOSE ( UNIT = iodsk )

      CALL EXIT(0)
c      stop
      END
c
c..............................................................
      SUBROUTINE svdinpts
c...............................................................

c... Writes out quantities needed for SVD solution for the inverse
c    Response matrix.
      
      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum5.inc'

      CHARACTER(1), PARAMETER :: tab = CHAR(9)

c      iosvd = 86

      nzn = n

      jmax1 = lmax(1) - lmin(1) + 1

      OPEN ( UNIT=86, FILE='svd-in_out',
     $     STATUS='REPLACE', FORM='FORMATTED' )

      WRITE ( iosvd, '(/,20x,"SVD INPUTS",/)' )

      WRITE ( iosvd, '(5x, i5, a )' ) mth,
     $     tab//tab//"mth"//tab//"[No. of plasma points]"
      WRITE ( iosvd, '(5x, i5, a )' ) mth,
     $     tab//tab//"mth"//tab//"[No. of observer points]"

      WRITE ( iosvd, '(5x, i5, a )' ) nzn,
     $     tab//tab//"n"//tab//"[Toroidal mode number]"
      WRITE ( iosvd, '(5x, 2i5, a )' ) lmin(1), lmax(1),
     $     tab//tab//"lmin, lmax"//tab//"[Fourier Range]"
      WRITE ( iosvd, '(5x, es15.7, a )' ) qa1,
     $     tab//tab//"qa1"//tab//"[Safety Factor]"

      WRITE ( iosvd, '(5x, /a/ )' )  "Plasma X coordinate:"
      WRITE ( iosvd, '(10es15.7)' ) (xinf(i), i = 1, mth1)

      WRITE ( iosvd, '(5x, /a/ )' )  "Plasma Z coordinate:"
      WRITE ( iosvd, '(10es15.7)' ) (zinf(i), i = 1, mth1)

      WRITE ( iosvd, '(5x, /a/ )' )  "Shell X coordinate:"
      WRITE ( iosvd, '(10es15.7)' ) (xwal(i), i = 1, mth1)

      WRITE ( iosvd, '(5x, /a/ )' )  "Shell Z coordinate:"
      WRITE ( iosvd, '(10es15.7)' ) (zwal(i), i = 1, mth1)

      WRITE ( iosvd, '(5x, /a/ )' )  "Delta: coordinate distortion:"
      WRITE ( iosvd, '(10es15.7)' ) (delta(i), i = 1, mth1)

      WRITE ( iosvd, '(5x, /a/ )' )  "GPSJP: Jacobian * Grad-Psi_plas"
      WRITE ( iosvd, '(10es15.7)' ) (gpsjp(i), i = 1, mth1)

      WRITE ( iosvd, '(5x, /a/ )' )  "GPSJW: Jacobian * Grad-Psi_wall"
      WRITE ( iosvd, '(10es15.7)' ) (gpsjw(i), i = 1, mth1)

      WRITE ( iosvd, '(5x, /a/ )' )  "Real Response(obs,srce in l):"
      DO i = 1, mth1
         WRITE ( iosvd, '(10es15.7)' ) ( cwallr(i,l1), l1=1,jmax1 )
      END DO 

      WRITE ( iosvd, '(5x, /a/ )' )  "Imag. Reponse(obs,srce in l):"
      DO i = 1, mth1
         WRITE ( iosvd, '(10es15.7)' ) ( cwalli(i,l1), l1=1,jmax1 )
      END DO 

c.... Plasma magnetic perturbations:

      WRITE ( iosvd, '(5x, /a/ )' )  "BNLR, real plasma harmonic:"
      WRITE ( iosvd, '(10es15.7)' ) (bnlr(i),i = 1, jmax1)

      WRITE ( iosvd, '(5x, /a/ )' )  "BNLI, imag plasma harmonic:"
      WRITE ( iosvd, '(10es15.7)' ) (bnli(i), i = 1, jmax1)

      WRITE ( iosvd, '(5x, /a/ )' )  "BNKR, real plasma perturbation:"
      WRITE ( iosvd, '(10es15.7)' ) (bnkr(i),i = 1, mth1)

      WRITE ( iosvd, '(5x, /a/ )' )  "BNKI, imag plasma perturbation:"
      WRITE ( iosvd, '(10es15.7)' ) (bnki(i), i = 1, mth1)

c... Need magnetic field components at wall. Get them from Subroutine
c    PICKUP after call to BTANG.

c      CLOSE ( UNIT = 86 )

      RETURN
      END SUBROUTINE svdinpts

c...............................................................
      SUBROUTINE dcnvac
c..............................................................

c Reads the output from DCON. Use it's coding, mostly.

      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum3.inc'

      INTEGER out_unit, ipert1, isol1

      CHARACTER*80 ztitle
      CHARACTER*9 ztit9
c      DIMENSION xilr(*), xili(*)
      INTEGER, DIMENSION(:), ALLOCATABLE :: izmax, mzz, mzzz
      REAL, DIMENSION(:), ALLOCATABLE :: ep, ev, et, evac, singfac
      REAL, DIMENSION(:), ALLOCATABLE :: ep0, ev0, et0
      REAL, DIMENSION(:,:), ALLOCATABLE :: wtr, wti, wtri, awt,
     $     zwrk, wpabs
      CHARACTER(1),DIMENSION(:,:),ALLOCATABLE :: star
      COMPLEX, DIMENSION(:,:), ALLOCATABLE :: wppc

      out_unit = 165
      OPEN ( out_unit, FILE='dcon_ev', STATUS='OLD', FORM='FORMATTED' )

c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(1x,"Energies: plasma = ",1pe10.3,", vacuum = ",e10.3,
     $     ", total = ",e10.3)
 20   FORMAT(/3x,"isol",3x,"plasma",5x,"vacuum",5x,"total"/)
 30   FORMAT(i6,1p3e11.3)
 40   FORMAT(/3x,"isol",2x,"imax",3x,"plasma",5x,"vacuum",5x,"total"/)
 50   FORMAT(2i6,1p3e11.3)
 60   FORMAT(/2x,"ipert",4x,"m",4x,"re wt",6x,"im wt",6x,"abs wt"/)
 70   FORMAT(2i6,1p3e11.3,2x,a)
 80   FORMAT(/3x,"isol",3x,"plasma",5x,"vacuum"/)
 90   FORMAT(i6,1p2e11.3)
 81   FORMAT(/3x,"isol",3x,"plasma",5x,"vacdcon", 4x,"vacvac",
     $     5x,"total"/)
 91   FORMAT(i6,1p2e11.3,e12.4,1e11.3)

c-----------------------------------------------------------------------

c  Define the "undefined" variables necessary to read file 'dcon_ev'.  Chance

      mlow = lmin(1)
      mpert = lmax(1) - lmin(1) + 1
      mpert2 = 2*mpert

      ALLOCATE ( izmax(mpert), mzz(mpert), mzzz(mpert) )
      ALLOCATE ( ep0(mpert), ev0(mpert), et0(mpert) )
      ALLOCATE ( ep(mpert), ev(mpert), et(mpert), singfac(mpert) )
      ALLOCATE ( evac(mpert), zwrk(mpert2,mpert2) )
      ALLOCATE ( wppc(mpert,mpert),
     $     wtr(mpert,mpert), wti(mpert,mpert),wtri(mpert2,mpert2),
     $     wpabs(mpert,mpert))
      ALLOCATE ( awt(mpert,mpert), star(mpert,mpert) )

c-----------------------------------------------------------------------
c     write to screen and copy to output.
c-----------------------------------------------------------------------
c$$$      WRITE(*,10)ep(1),ev(1),et(1)
c$$$      plasma1=ep(1)
c$$$      vacuum1=ev(1)
c$$$      total1=et(1)
c-----------------------------------------------------------------------
c     Read eigenvalues from file.
c-----------------------------------------------------------------------
c      READ(out_unit,*)"Total Energy Eigenvalues:"
      READ(out_unit,'(/,a,/)') ztitle
      WRITE(outmod,'(/,a,/)') ztitle
      WRITE(iotty,'(/,a,/)') ztitle
      READ(out_unit,'(a)')ztitle
      READ(out_unit,'(/,a,/)')ztitle
      WRITE(iotty,'(/,a,/)') ztitle

      WRITE(iotty,'("il=",i5," ih=",i5," mpert=",i5)')
     &     lmin(1),lmax(1),mpert

      DO isol=1,mpert
         READ(out_unit,30)iso,ep0(isol),ev0(isol),et0(isol)
         WRITE(iotty,30)iso,ep0(isol),ev0(isol),et0(isol)
      END DO

c$$$      READ(out_unit,30)
c$$$     $     (isol,ep0(isol),ev0(isol),et0(isol),iso=1,mpert)

      READ(out_unit,'(/,a,/)')ztitle
      WRITE(iotty,'(/,a,/)') ztitle
c-----------------------------------------------------------------------
c     read eigenvectors to file.
c-----------------------------------------------------------------------
c      READ(out_unit,*)"Total Energy Eigenvectors:"
      READ(out_unit,'(a)')ztitle
      mzz=mlow+(/(isol,isol=0,mpert-1)/)
      DO isol=1,mpert
         READ(out_unit,'(/,a,/)')ztitle
         READ(out_unit,50)isolz,izmax(1),ep(isol),ev(isol),et(isol)
         WRITE(iotty,50)isolz,izmax(1),ep(isol),ev(isol),et(isol)
         READ(out_unit,'(/,a,/)')ztitle
c$$$         WRITE(iotty,'(/,a,/)')ztitle
         DO iper=1,mpert
         READ(out_unit,70)ipert,mzz(iper),wtr(iper,isol),
     $        wti(iper,isol),
     $        awt(iper,isol),star(iper,isol)
c$$$         WRITE(iotty,70)ipert,mzz(iper),wtr(iper,isol),
c$$$     $        wti(iper,isol),
c$$$     $        awt(iper,isol),star(iper,isol)
         END DO
         READ(out_unit,'(/,a,/)')ztitle
      END DO

c Store the first eigenvectors in xilr, xili

      xilr(1:mpert) = wtr(1:mpert,1)
      xili(1:mpert) = wti(1:mpert,1)

c..  Construct array of real and imaginary parts in xilri:
      xilri(1:mpert) = xilr(1:mpert)
      xilri(mpert+1:mpert2) = xili(1:mpert)

c.. Now put the real an imaginary parts of wt in wtri:
      DO isol = 1, mpert
         wtri(1:mpert,isol) = wtr(1:mpert,isol)
         wtri(mpert+1:mpert2,isol) = wti(1:mpert,isol)
      END DO

c Read DCON's plasma matrix:

c$$$c.. Mine:
c$$$      READ (out_unit,'(a)') ztit9
c$$$      WRITE(iotty,'(/,a)') ztit9
c$$$c..End if mine

c..  Ming's:
      READ(out_unit,'(/a/)') ztitle
      WRITE(iotty,'(/,a)') ztitle
c..  End of Mings

c      mzz=mlow+(/(isol,isol=0,mpert-1)/)
      DO isol=1,mpert
c........Comment out to use ming's

c$$$         READ(out_unit,'(a9, i5)') ztit9, mzzz(isol)
c$$$         READ(out_unit,'(1p5e12.5)')
c$$$     $        (wppc(ipert,isol),ipert=1,mpert)
c         WRITE(iotty,'(/,"at A",a,2i5)') ztit9, isol, mzzz(isol)

c....   End comment

c.. Ming's modification: for new DCON read:
         READ(out_unit,'(1x,a7,i3,a6,i3)')isol1,mzzz(isol)
         READ(out_unit,'(/a/)')ztitle
         READ(out_unit,'(i3,1p,3e13.5)')
     $        (ipert1,wppc(ipert,isol),wpabs(ipert,isol),ipert=1,mpert)
         READ(out_unit,'(/a/)')ztitle
         WRITE(iotty,'(1x,2(a,i3))')"isol = ",isol,", m = ",mzzz(isol)
         WRITE(iotty,'(/2x,"i",5x,"re wp",8x,"im wp",8x,"abs wp"/)')
         WRITE(iotty,'(i3,1p,3e13.5)')
     $      (ipert,wppc(ipert,isol),wpabs(ipert,isol),ipert=1,mpert)
         WRITE(iotty,'(/2x,"i",5x,"re wp",8x,"im wp",8x,"abs wp"/)')
c... End of ming's modification

      END DO

      wpp(1:mpert,1:mpert) = REAL ( wppc(1:mpert,1:mpert) )
      wppi(1:mpert,1:mpert) = IMAG ( wppc(1:mpert,1:mpert) )

      WRITE(outmod,*)"Read Plasma Energy Matrix:"

      CALL matwrtn ( wpp, nfm,nfm, mlow,mlow,mpert,mpert,8,8,
     $     "DCONs plasma matrix", outmod,0 )
      call matwrts ( wpp, nfm,nfm, 10,10, 18, 18, 9,9,
     $     10,10,"DCONs plasma matrix", outmod,0 )
      call matwrts ( wpp, nfm,nfm, 19, 19, 27, 27, 9,9,
     $     19,19,"DCONs plasma matrix", outmod,0 )

c. Use \calB as basis for wpp0:

      WRITE ( iotty,  '(/, "mlow, n, qa1 = ", i5, 2es12.4)' )
     $     mlow, n, qa1
      WRITE ( outmod, '(/, "mlow, n, qa1 = ", i5, 2es12.4)' )
     $     mlow, n, qa1

      singfac =
     $     1.0 / ( mlow - n*qa1 + (/ (ipert,ipert=0,mpert-1) /) )
      WRITE ( iotty, '(/, "l-nq = ",/,(10es11.4))' )
     $     ( singfac(i),i=1,mpert )
      WRITE ( outmod, '(/, "l-nq = ",/,(10es11.4))' )
     $     ( singfac(i),i=1,mpert )

      DO ipert = 1, mpert
         wpp0(ipert,1:mpert) = wpp(ipert,1:mpert)*singfac(1:mpert)
         wpp0i(ipert,1:mpert) = wppi(ipert,1:mpert)*singfac(1:mpert)
      END DO
      DO ipert = 1, mpert
         wpp0(1:mpert,ipert) = singfac(1:mpert)*wpp0(1:mpert,ipert)
         wpp0i(1:mpert,ipert) = singfac(1:mpert)*wpp0i(1:mpert,ipert)
      END DO

      CALL matwrtn ( wpp0, nfm,nfm, mlow,mlow,mpert,mpert,8,8,
     $     "DCONs plasma matrix w/o l-nq", outmod,0 )
      call matwrts ( wpp0, nfm,nfm, 10, 10, 18, 18, 9,9,
     $     10,10,"DCONs plasma matrix w/o l-nq", outmod,0 )
      call matwrts ( wpp0, nfm,nfm, 19, 19, 27, 27, 9,9,
     $     19,19,"DCONs plasma matrix w/o l-nq", outmod,0 )

c$$$      DO isol=1,mpert
c$$$         WRITE(outmod,'(" m column", i5)')mzzz(isol)
c$$$         WRITE(outmod,'(1p5e12.5)')(wppc(ipert,isol),ipert=1,mpert)
c$$$      END DO

c-----------------------------------------------------------------------
c     compute and print separate plasma and vacuum eigenvalues.
c-----------------------------------------------------------------------

c.. Reconstruct the vacuum energy with the "total" eigenvectors.
c    Use msc's vacmat and ahg's "wt" eignevectors.
c Assume only real XI (aka wtr) for now.

c... Comment out up-down symmetric stuff
c$$$      DO isol = 1, mpert
c$$$         zwrk(1:mpert,isol) =
c$$$     $        MATMUL ( vacmat(1:mpert,1:mpert), wtr(1:mpert,isol) )
c$$$         evac(isol) =
c$$$     $        DOT_PRODUCT ( wtr(1:mpert,isol),zwrk(1:mpert,isol) )
c$$$      END DO

c... Use this for up-down asymmetry:
      DO isol = 1, mpert
         zwrk(1:mpert2,isol) =
     $        MATMUL ( vcmattu(1:mpert2,1:mpert2),wtri(1:mpert2,isol) )
         evac(isol) =
     $        DOT_PRODUCT ( wtri(1:mpert2,isol),zwrk(1:mpert2,isol) )
      END DO

c$$$      CALL zheev('V','U',mpert,wp,mpert,ep,work,lwork,rwork,info)
c$$$      CALL zheev('V','U',mpert,wv,mpert,ev,work,lwork,rwork,info)
c      READ(out_unit,*)"Separate Energy Eigenvalues:"
      READ(out_unit,'(a)')
      READ(out_unit,'(/,a,/)')

      DO isol = 1, mpert
         READ(out_unit,90)iso,ep(isol),ev(isol)
      END DO

c$$$      READ(out_unit,90)(isol,ep(isol),ev(isol),iso=1,mpert)

      READ(out_unit,'(/,a,/)')

      WRITE(outmod,*)"Separate Energy Eigenvalues: DCON"
      WRITE(outmod,80)
      WRITE(outmod,90)(isol,ep(isol),ev(isol),isol=1,mpert)
      WRITE(outmod,80)

      WRITE(outmod,*)"Vacuum Energy: VACUUM, Comparing Top Values"
      WRITE(outmod,81)
      WRITE(outmod,91)(isol,ep0(isol),ev0(isol),evac(isol),et0(isol),
     $     isol=1,mpert)
      WRITE(outmod,81)

      WRITE(iotty,*)"Vacuum Energy: VACUUM, Comparing Top Values"
      WRITE(iotty,81)
      WRITE(iotty,91)(isol,ep0(isol),ev0(isol),evac(isol),et0(isol),
     $     isol=1,mpert)
      WRITE(iotty,81)

      DEALLOCATE ( izmax, mzz,mzzz, ep0, ev0, et0, ep, ev, et )
      DEALLOCATE ( wppc, wtr, wti, wtri, awt, star, singfac )
      DEALLOCATE ( evac, zwrk )
      WRITE(outmod,*)"END OF dcnvac"

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dcnvac
c-----------------------------------------------------------------------

c................................................................
      subroutine defglo
c......................................................................
c      this subroutine defaults logical unit numbers
c      and logical control switches for package.
c......................................................................
c
      include 'vacuum1.inc'
      include 'vacuum2.inc'
      include 'vacuum3.inc'
      include 'vacuum5.inc'
      include 'vacuum7.inc'
c
c     1.1.1   constants
c
      zero   = 0.0e0
      pt1    = 1.e-1
      half   = 0.5e0
      one    = 1.0e0
      two    = 2.0e0
      three  = 3.0e0
      four   = 4.0e0
      five   = 5.0e0
      seven  = 7.0e0
      epsq = 1.0e-5
      epszer = 1.0e-30
c
      pye    = 3.141592653589e0
      twopi  = two * pye
      twopi2 = twopi * twopi
      tpisqi = 1.0 / (twopi*pye)
      alx     = 1
      alz     = 1
      n       = 1
      m       = 2
      mp     = 3
      minc   = 10
      mth    = 128
      mth1   = mth + 1
      mth2   = mth1 + 1
      mthin  = 128
      mthin1 = mthin + 1
      mthin2 = mthin + 2
      mthout = 100
      mw = 129
      mdiv   = 2
      dm0 = 0.0
      dx0 = 0.0
      dx1 = 0.0
      idgt = 0
      nosurf = ( mp - 1 ) * mdiv + 1
      n1surf = 4
      npsurf = 6
      dpsi   = one / mdiv / m
      bit    = 1.0e-8
      jwal1 = 1
      jwal2 = 2*jwal1
c
      awx = 0.0
      awv = 0.0
      nwcur = 1
      amu0   = four * pye * 1.e-7
c..   Resistivity and thickness of DIII-D vacuum vessel in MKS
      eincnl = 1.03e-6
      thdiiid = 1.0e-2
      eta_vess = eincnl
      th_vess = thdiiid
      gmat=10000.
      irotplasma=0
      taugam = 10.0
c
      dth    = twopi / mth
      dthinc = dth  / minc
      gamma  = five / three
      p0     = zero
      upsiln = one
      lj = 0
      mj = 0
      nj = 0
      r      = 20.0e0
      r2     = r * r
      r4     = r2 * r2
      r6     = r4 *r2
      rgato = r
      fa1 = 1.0
      ga1 = 1.0
      qa1 = 1.0
      s_scale = 1.0
      zpoff = 0.0         ! Vertical offset of plasma when ipshp=1
      psipr = 1.0
      isymz  = 2
      xiinr = 0.0
      xiini = 0.0
      xma = 1.0
      zma = 0.0
      xzero = 1.0
      ntloop = 8
      deloop = .001
      adel = 0.0
      nxlpin = 21
      nzlpin = 21
      xlpmin = 0.5
      xlpmax = 2.5
      zlpmin = -1.5
      zlpmax =  1.5
      wcentr = 0.0
      wmaj = 0.0
      wlrad = 0.0
      welong = 0.0
      wdness = 0.0
      epslp = 0.02
      linterior = 0
      lchkinextp = 0
      lclosepts = 0
      fpl = 0.0
      epl = 0.0
c
      mtcoil = 512
      mtcoil1 = mtcoil + 1
      coilang = 0.25
      coilcen= 0.0
      taucoil = 0.02
      acoil = 2.0
      bcoil = 1.5
      dcoil = 0.4
      drcoil = 1.0
      dzcoil = 0.0
      dzcoila = 0.0
      dzcoilb = 0.0
      
c.... Default feedback coil values. B is upper, A is lower
      
      lfbcoil = 2
      ctip_c(1) =  3.23
      ctip_c(2) =  0.8
      ctip_c(3) =  3.23
      ctip_c(4) = -0.8
      ntip_c = 51
      
      lfbcoilb = 2
      ctip_u(1) =  2.184
      ctip_u(2) =  1.012
      ctip_u(3) =  2.394
      ctip_u(4) =  0.504
      ntip_u = 51

      lfbcoila = 2
      ctip_l(1) =  2.394
      ctip_l(2) = -0.504
      ctip_l(3) =  2.184
      ctip_l(4) = -1.012
      ntip_l = 51

c      default to force i/o into separate files for equilibrium,
c      mapping and normal modes.
c
      idsk = 1
      intty = 5
      iotty = 6
      inpest = 2
      outpest = 3
      iomode = 19
      inmode = 22
      outmod = 23
      iodsk = 16
      ioshel = 20
      outmap1 = 18
      iovac = 36
      iosvd = 86
      do i = 1, 3
         nout(i) = 0
         nout0(1) = 0
      end do
      nout(1) = 6
      nout(2) = outmod
      mp0 = 'mapdsk'
      mp1 = 'mpout1'
c
c      logicals
c
      lpestasym = 0
      lzio = 1
      lsymz  = .false.
      check1 = .false.
      check2 = .false.
      lanal   = .false.
      lkdis = 1
      ldisc = 0
      lpest1 = .false.
      lnova = .false.
      lrnova = 0
      lrwall7 = 0
      lreig7 = 0
      lspark= 0
      lhighn = 0
      lreshel = 0
      leqarcw = 0
      lsymwv = 0
      ladj = 0
      ldcon = 0
      lrnim3d = 0
      ldqdtw = 0
      ltnrm = 1
      ltfil = 0
      lsensplot = 1
      lrswcomp = 1
      lpolsensor = 0
      wall   = .false.
      lkplt = 0
      ldelta = 0
c
      do ich = 1, 60
         seps(ich:ich) = '.'
      end do
c
      oplp = 0.0

      lmin = 0.0
      lmax = 0.0

      xplap = 0.0
      zplap = 0.0
      xinf = 0.0
      zinf = 0.0
      xwal = 1.0
      zwal = 1.0
      xwalp = 0.0
      zwalp = 0.0

      do i = 1, nths
         grpssq(i) = 2.5
         delta(i) = 0.0
         xjacob(i) = 1.0
      end do
c
      do j1 = 1, nfm
         xilr(j1)  = 0.0
         xili(j1) = 0.0
         bnlr(j1) = 0.0
         bnli(j1) = 0.0
      end do
      do j1 = 1, nfm2
         xilri(j1)  = 0.0
         bnlri(j1) = 0.0
      end do

      xloop = 0.0
      zloop = 0.0
      zloopnl = 0.0
      xloopnl = 0.0

      fv = 0.0

c
c      the following are pointers to tell how many times the code has
c      gone through equilibrium, mapping and normal modes calculations.
c
c     intialize the plotting....
c     pname = 8rgraphics
c     call paper ( pname )
c
      return
      end
c
c....................................................................
c     1.4   function initialization.
c.......................................................................
c
c     this subroutine will set up various arrays of quantities needed
c     in setting up the matrix elements.  if lfunin is set to true,
c     certain functions which otherwise would have been taken from
c     disk storage will be overwritten here.
c.......................................................................
c***********************************************************
c       version for numerically calculated equilibrium
c***********************************************************
c
c.......................................................................
      SUBROUTINE funint
c.......................................................................
      INCLUDE 'vacuum1.inc'

      ga1 = r*ga1
c$$$      ga1p = r*ga1p

      RETURN
      END
c     2.1   variation of parameters.
c.......................................................................
c
      SUBROUTINE varymo
c.......................................................................
      INCLUDE 'vacuum1.inc'
      
      INTEGER :: iref = 0
c
      IF ( iref > 0 ) GO TO 10
      
      IF ( lpest1 ) THEN
      WRITE ( OUTMOD, '(/, " UPSILN = ", ES12.5, " FA1 = ", ES12.5,
     $     " GA1 = ", ES12.5, " QA1(upsg/2pif?) = ", ES12.5 )' )
     $     upsiln, fa1, ga1, qa1
      WRITE ( IOTTY, '(/, " UPSILN = ", ES12.5, " FA1 = ", ES12.5,
     $     " GA1 = ", ES12.5, " QA1(upsg/2pif?) = ", ES12.5 )' )
     $     upsiln, fa1, ga1, qa1

c..   Redefine fa1. fa1 not used in rest of VACUUM anyway. !! MSC 6/17/2011
c$$$         fa1 = twopi*r*fa1 / upsiln
c$$$         qa1 = ga1/fa1
c... QA1 is off a bit so use:
         fa1 = ga1 / qa1
      END IF
      WRITE ( OUTMOD, '(/," FA1 = ", ES12.5,
     $     " GA1 = ", ES12.5, " QA1(g/f) = ", ES12.5 )' )
     $     fa1, ga1, qa1
      WRITE ( IOTTY, '(/, " FA1 = ", ES12.5,
     $     " GA1 = ", ES12.5, " QA1(g/f) = ", ES12.5 )' )
     $     fa1, ga1, qa1

      IF ( ipshp == 1 ) qa1 = qain
c         f02 is f squared at interface..
      f02 = fa1**2

      ga1ref = ga1
c$$$      ga1pref = ga1p

      iref = iref + 1

   10 CONTINUE

      IF ( ABS(s_scale-one) <  1.0e-6 ) GO TO 50

c ..... now use the s scaling parameter...

      s2 = s_scale * s_scale
      ga1 = SQRT( ga1ref*ga1ref - (one - s2 ) *r2)

c$$$      ga1p = ga1ref * ga1pref / ga1
c$$$      z1 = ga1/ga1ref
c$$$      z2 = (ga1ref*ga1p - ga1pref*ga1) / (fa1*ga1ref)
c$$$      DO 11 i=1,nths
c$$$   11 qdelp(i) = z1*qdelp(i) + z2*delta(i)

   50 CONTINUE

      qa1 = ga1 / fa1

c$$$      qp = ( fa1*ga1p - ga1*fa1p ) / ( fa1*fa1 )
c$$$      fa12 = fa1*fa1
c$$$      ga12 = ga1*ga1
c$$$      q2 = q*q
c$$$      IF ( qpmax < qp ) qpmax = qp

       WRITE ( OUTMOD, 15 ) s_scale, ga1, fa1, qa1
 15    FORMAT ( /,1x, "s_scale, ga1, fa1, qa1, = ",/,
     $      (5ES12.5) )

      RETURN
      END

c$$$      subroutine funint
c$$$c.......................................................................
c$$$      include 'vacuum1.inc'
c$$$      include 'vacuum3.inc'
c$$$      include 'vacuum2.inc'
c$$$      integer t, surf
c$$$c
c$$$c     1.4.1   functions not fixed by the specific equilibrium.
c$$$c
c$$$c          upsiln is scaled out of the expressions for deltaw
c$$$c         by working with f/upsiln and with metric elements
c$$$c          multiplied by upsiln**2.
c$$$c
c$$$      upsil2 = upsiln**2
c$$$c
c$$$c      qa1 = r * ga1 / fa1
c$$$
c$$$c...Comment out next statement since r and upsiln seems to be undetermined
c$$$c             .....04-30-2004
c$$$
c$$$c      if ( lpest1 ) qa1 = upsiln*ga1/(twopi*fa1)
c$$$      if ( ipshp .eq. 1 ) qa1 = qain
c$$$c
c$$$c         f02 is f squared at interface..
c$$$      f02 = fa1**2
c$$$c
c$$$      return
c$$$      end
c
c.......................................................
      SUBROUTINE arrays ( icall )
c........................................................
c
      include 'vacuum1.inc'
      include 'vacuum2.inc'
      include 'vacuum3.inc'
      include 'vacuum5.inc'
      include 'vacuum7.inc'
      include 'vacuum8.inc'
c
      data izplot /1/
c
      real nq
      REAL, SAVE :: xmin, xmax, zmin, zmax    !  08/26/2009
      dimension xpla(nths),zpla(nths)
      dimension zorkx(nths2), zorkz(nths2),
     $     dlenth(nths)
c      dimension ellqc(nths), qtgrc(nths)
c
      common / bigv1 / grdgre(nths2,nths2), grwp(nths,nths),
     $     grri(nths2,nfm2)
c
      equivalence (xpla,xinf), (zpla,zinf)
c
c..   Arrays is called several times. Before VACCAL, after VACCAL ,
c     then before the diagnostic section if LSPARK is turned on.
c
      call atpoint ( "ARRAYS", "ICALL", icall, "XINF(i)",xinf(1),
     $     iotty, outmod )
c$$$c
c$$$      if ( icall .eq. 2 ) go to 100
c
c..   New variable, lfarw for SPARK2 to reset farwal.
c
c$$$      lfarw = 0
c$$$      if ( a .gt. 10.0 ) lfarw = 1
c$$$      write ( iotty, '("a = ", 1pe12.3, " lfarw in ARRAYS = ", i4 )' )
c$$$     $     a, lfarw
c$$$      write ( outmod, '("a = ", 1pe12.3, " lfarw in ARRAYS = ", i4 )' )
c$$$     $     a, lfarw
c
      mth1 = mth + 1
      mth2 = mth + 2
      mth3 = mth + 3
      mth4 = mth + 4
      mth5 = mth + 5
      mth12 = 2*mth
c
c      jmax1 = lxsav - lnsav + 1
      jmax1 = lmax(1) - lmin(1) + 1
      lrnge = jmax1
      q = qa1
      nq = n*q
c
      if ( icall .eq. 2 ) go to 100   ! moved here on 09-30-2008
c
c...  Abscissa for l values:
      do j1 = 1, jmax1
         lfm(j1) = lmin(1) - 1 + j1
         xlfm(j1) = lfm(j1)
      end do

c
      do i = 1, nths
         zerov(i) = 0.0
      end do
c
      do i = 1, mth2
         pgrd(i) = ( 1.0/mth ) * (i-1)
         pigrd(i) = (i-1) * dth
      end do
c
      call boundsi(xpla,zpla,1,mth,xmnp,xmxp,zmnp,zmxp,
     $     ixn,ixx,izn,izx)
      xmin = xmnp
      xmax = xmxp
      zmin = zmnp
      zmax = zmxp
c
      plrad = 0.5 * ( xmxp - xmnp )
      xmaj = 0.5 * ( xmxp + xmnp )
      delx = plrad * delfac
      delz = plrad * delfac
      pelong = (zmxp-zmnp) / (xmxp-xmnp)
      pdness = asin ( (xmaj - xpla(izx))/plrad )
      write ( iotty, 3 ) delx, delz
      write ( outmod, 3 ) delx, delz
 3    format ( 1x, "  delx, dely = ", 1p2e13.5 )
c
c.....calculate arc length on surface..
c
c......Get derivative on plasma points.
c......Calculate Jacobian*Grad-Psi. Store in gpsjp and gpsjw
c
      CALL difspl ( mth, pigrd, xpla, xplap )
      CALL difspl ( mth, pigrd, zpla, zplap )

c... Strore away dlplasma/dtheta:

      dlpdth = SQRT ( xplap**2 + zplap**2 )

c.. Get table of second derivatives of qth0 for the high n stuff:

c      qth0(1:mth2) = 0.5 * sin ( pigrd(1:mth2) )
      qth0 = 0.0
      qth0pp = 0.0
      enqth(1:mth2) = qa1 * ( qth0(1:mth2) + pigrd(1:mth2) )
c$$$      WRITE ( OUTMOD, '("qth0 = ",/, (10es12.4))' )
c$$$     $     (qth0(i), i = 1, mth2)
c$$$      WRITE ( OUTMOD, '("ENQTH = ",/, (10es12.4))' )
c$$$     $     (enqth(i), i = 1, mth2)

      IF ( lhighn /= 0 )
     $     CALL difspl2 ( mth, pigrd, qth0, zorkx, qth0pp )

c$$$      WRITE ( OUTMOD, '("qth0pp = ",/, (10es12.4))' )
c$$$     $     (qth0pp(i), i = 1, mth2)

c      CALL pospl1 ( pigrd,zth,    1,mth1, 1,"zth",   mth1 )
      CALL pospl1 ( pigrd,enqth,  1,mth1, 2,"enqth", mth1 )
      CALL pospl1 ( pigrd,qth0,   1,mth1, 3,"qth0",  mth1 )
      CALL pospl1 ( pigrd,qth0pp, 1,mth1, 4,"qth0pp",mth1 )
      CALL framep(jobid,ff)

c
c$$$      call intset ( xpla, zorkx, mth, twopi, 4 )
c$$$      call intset ( zpla, zorkz, mth, twopi, 4 )
      do i = 1, mth1
c$$$         call lagpe5 ( xpla, zorkx, mth, twopi, pigrd(i), dum,
c$$$     $        xplap(i), 1, 4 )
c$$$         call lagpe5 ( zpla, zorkz, mth, twopi, pigrd(i), dum,
c$$$     $        zplap(i), 1, 4 )
c$$$         call lagp ( pigrd,xpla,mth1,2,pigrd(i),dum,xplap(i),1, 1 )
c$$$         call lagp ( pigrd,zpla,mth1,2,pigrd(i),dum,zplap(i),1, 1 )
         dlenth(i) = sqrt ( xplap(i)**2 + zplap(i)**2 )
         gpsjp(i) = xpla(i) * dlenth(i)
      end do

c.. Set up extra points for the indefinite integral:

      do i = 1, mth1
         zork1(i+2) = dlenth(i)
      end do
c
      zork1(1) = dlenth(mth-1)
      zork1(2) = dlenth(mth)
      zork1(mth3) = dlenth(1)
      zork1(mth4) = dlenth(2)
      zork1(mth5) = dlenth(3)

c.. Integrate from 3 to mth3 and store properly:

      call indef4 ( zork1, zork2, dth, 3,mth3, alen, 0 )
c
      do i = 1, mth1
         slngth(i) = zork2(i+2)
      end do
c
c... Ellp will contain the arc length of the plasma surface starting
c    from the inside:

      mthh1 = mth/2 + 1
      CALL shflen ( slngth, ellp, mth1, mthh1 )

      pcircm = slngth(mth1)
      write( iotty,'("Circumference of Plasma = ",1pe11.3)')
     $     pcircm
      write(outmod,'("Circumference of Plasma = ",1pe11.3)')
     $     pcircm
c
c.......store away sin and cos. only need for source points.
c
      do is = 1, mth1
         theta = (is-1) * dth
         znqd = nq*delta(is)
         cnqd(is) = cos(znqd)
         snqd(is) = sin(znqd)
         amxdelta = MAXVAL(ABS(delta))
         if ( grpssq(is) .lt. 1.e-20 ) then
            grpssq(is) = 1.0e-20
            write ( outmod, '(1x,"grpssq(i) .lt. 1.e-20 at ", i3)' ) is
            write ( iotty,  '(1x,"grpssq(i) .lt. 1.e-20 at ", i3)' ) is
         end if
         sqgrps(is) = sqrt ( grpssq(is) )
         do l1 =  1, jmax1
            ll = lmin(1) - 1 + l1
            elth = ll * theta
            elthnq = ll * theta + znqd
            sinlt(is,l1) = sin(elth)
            coslt(is,l1) = cos(elth)
            snlth(is,l1) = sin(elthnq)
            cslth(is,l1) = cos(elthnq)
         end do
      end do

      WRITE ( outmod, '(/,"Max Val of Delta = ", es11.4,/)' ) amxdelta
      WRITE ( iotty , '(/,"Max Val of Delta = ", es11.4,/)' ) amxdelta
c
      if ( lplot1 .eq. 0 ) go to 50
c
c... Plot only once.
c
      if ( izplot .gt. 1 ) go to 50
c
      izplot = izplot + 1
c
c... Plot some quantities versus grid. Shift so that mth1/2 is outside.
c
      call shftpi ( pigrd,xjacob, zork1,zork2, mth1 )
      call pospl2 ( zork1,zork2,zerov, 1,mth1, 1,'xjacob',mth1 )
      call shftpi ( pigrd,delta, zork1,zork2, mth1 )
      call pospl2 ( zork1,zork2, zerov, 1,mth1, 2,'delta',mth1 )
      call shftpi ( pigrd,sqgrps, zork1,zork2, mth1 )
      call pospl2 ( zork1,zork2,zerov, 1,mth1, 3,'grad-psi',mth1 )
      call shftpi ( pigrd,gpsjp, zork1,zork2, mth1 )
      call pospl2 ( zork1,zork2, zerov, 1,mth1, 4,'j*gpsi-p',mth1 )
      call framep(jobid,ff)
c
 50   continue
c
      if ( lfarw .eq. 0 ) then
c
         dthw = dth
c
         do i = 1, mw
            wgrd(i) = ( 1.0/(mw-1) ) * (i-1)
            piwgrd(i) = (i-1) * dthw
         end do

c......Get derivative on wall points. Calculate J-GradPsi.
c
         call difspl ( mth, pigrd, xwal, xwalp )
         call difspl ( mth, pigrd, zwal, zwalp )
         xwalp(mth+1) = xwalp(1)
         zwalp(mth+1) = zwalp(1)
c
         do i = 1, mth1
            gpsjw(i) = xwal(i) * sqrt ( xwalp(i)**2 + zwalp(i)**2 )
         end do
c
c     Get arc length on wall, ell, shifted coordinates...
         call shfgrb ( xwal, zork1, mth )
         call shfgrb ( zwal, zork2, mth )
c
         ell(1) = 0.0
         do iw = 1, mw
            if ( iw .eq. 1 ) go to 8
            thet = ( wgrd(iw)+wgrd(iw-1) ) / 2.0
            call lag ( wgrd,zork1,mw,3,thet,f,df,1 )
            xtzt = df
            call lag ( wgrd,zork2,mw,3,thet,f,df,1 )
            xtzt = sqrt ( xtzt**2 + df**2 )
            ell(iw) = ell(iw-1) + xtzt/(mw-1)
 8          continue
         end do
         write(outmod,19)
     $        (i,wgrd(i),zork1(i),zork2(i),ell(i),i=1,mw,4)
 19      format(//,2x,"wgrd , zork1, zork2, ell  of wall =",/,
     $        200(2x,i4,1p4e12.4,/))
c
c. Need Equal arcs for wall even if option is not turned on
c. Get equal wall properties if wall is not yet so parametrized.
c    Store in xwalq, zwalq. Derivatives in xwalpq, zwalpq

         xwalq = xwal
         zwalq = zwal
         xwalpq = xwalp
         zwalpq = zwalp

         IF ( leqarcw == 1 ) GO TO 40

         CALL eqarcw ( xwal,zwal, xwalq,zwalq, elwal,zork2,wkgrd, mth1 )

c......Get derivative on equal arc wall points.
c
         CALL difspl ( mth, pigrd, xwalq, xwalpq )
         CALL difspl ( mth, pigrd, zwalq, zwalpq )
         xwalpq(mth+1) = xwalpq(1)
         zwalpq(mth+1) = zwalpq(1)

 40      CONTINUE

         wcircm = ell(mw)
         write( iotty,'("Circumference of Wall = ",1pe11.3)')
     $        wcircm
         write(outmod,'("Circumference of Wall = ",1pe11.3)')
     $        wcircm

      end if

      if ( lreshel .eq. 0 ) go to 60
c
c...Resistive shell Eigenmodes, even and odd pairs
c      : Read in from disc: Dummy stuff for now
c
      call waleig
c
c$$$      jwal1 = jmax1
c$$$      do j = 1, jwal1
c$$$         rsval(j) = 1.0
c$$$         do i = 1, mth1
c$$$            rwvece(i,j) = coslt(i,j)
c$$$            rwveco(i,j) = sinlt(i,j)
c$$$         end do
c$$$      end do
c
 60   continue
c
      if ( icall .eq. 1 ) return
c
 100  continue
c
      jmax1 = lmax(1) - lmin(1) + 1   ! 08/26/2009, since maybe not SAVED
      if ( lfarw .eq. 0 ) then
c
c     Get arc length on wall, ell, in passed shifted coordinates...
c.      This will overwrite previous evaluation!!
c
          ell(1) = 0.0
         do iw = 1, mw
            if ( zpass(iw) .gt. zmax ) zmax = zpass(iw)
            if ( iw .eq. 1 ) go to 15
            thet = ( wgrd(iw)+wgrd(iw-1) ) / 2.0
            call lag ( wgrd,xpass,mw,3,thet,f,df,1 )
            xtzt = df
            call lag ( wgrd,zpass,mw,3,thet,f,df,1 )
            xtzt = sqrt ( xtzt**2 + df**2 )
            ell(iw) = ell(iw-1) + xtzt/(mw-1)
 15         continue
         end do
c
c         if ( check2 )
             write(outmod,16)
     $        (i,wgrd(i),xpass(i),zpass(i),ell(i),i=1,mw,4)
 16      format(//,2x,"wgrd , xpass, zpass, ell =",/,
     $        200(2x,i4,1p4e12.4,/))
c
         end if
c
c     Store C response matrices from plasma sources. 
c
      do l1 = 1, jmax1
         do i = 1, mth
            cplar(i,l1) = grri(i,l1)
            cplai(i,l1) = grri(i,jmax1+l1)
            if ( lfarw .eq. 0 ) then
               cwallr(i,l1) = grri(mth+i,l1)
               cwalli(i,l1) = grri(mth+i,jmax1+l1)
            end if
         end do
         cplar(mth1,l1) = cplar(1,l1)
         cplai(mth1,l1) = cplai(1,l1)
         if ( lfarw .eq. 0 ) then
            cwallr(mth1,l1) = cwallr(1,l1)
            cwalli(mth1,l1) = cwalli(1,l1)
         end if
      end do

c... If resisitve wall then B_n is nonzero so store the response from
c    the even and odd wall eigenfunctions

      IF ( lreshel == 0 ) GO TO 23

      DO l1 = 1, jwal1
         DO i = 1, mth
            cplae(i,l1) = grri(i,2*jmax1+l1)
            cplao(i,l1) = grri(i,2*jmax1+jwal1+l1)
            IF ( lfarw == 0 ) THEN
               cwwale(i,l1) = grri(mth+i,2*jmax1+l1)
               cwwalo(i,l1) = grri(mth+i,2*jmax1+jwal1+l1)
            END IF
         END DO
         cplae(mth1,l1) = cplae(1,l1)
         cplao(mth1,l1) = cplao(1,l1)
         IF ( lfarw == 0 ) THEN
            cwwale(mth1,l1) = cwwale(1,l1)
            cwwalo(mth1,l1) = cwwalo(1,l1)
         END IF
      END DO


 23   CONTINUE

      call atpoint ( "End of Arrays", "mw", mw, "qa1", qa1,
     $     iotty, outmod )
c
      return
      end
c     
c...................................................................
      SUBROUTINE b_lines ( zdelta, zqa, nzpts1, nzlin, jobid )
c..................................................................
c
c.. Plots equilibrium magnetic field lines on a rectilinear phi-theta
c   surface in configuration space. 
c.. zqa: safety factor value.
c.. zdelta: deviation from the PEST theta. I.e., the delta from PEST II
c.  zdelta is shifted here to make the outer major radius centered in the 
c   plots.
c   nzpts: number of points in zdelta on closed region. mth1, say
c   nzlin: number of field lines.

      DIMENSION zdelta(*)
      REAL, DIMENSION(:), ALLOCATABLE :: ztheta, zphi, zdeltas, zzero
      CHARACTER*80 strg
      CHARACTER*(*) jobid
      CHARACTER*(80) stringv(5)
      CHARACTER*(20) pllabel0, pllabel
      LOGICAL backl
c
c.....plot each eigenvector as function of Fourier index
c
      zxmin = -0.1
      zymin = -0.1
      zxmax =  1.1
      zymax =  1.1

      nzpts = nzpts1 - 1
      ztpi = 2.0 * acos (-1.0)
      backl = .false.
c
      ALLOCATE ( ztheta(nzpts+5), zphi(nzpts+5), zdeltas(nzpts+5),
     $     zzero(nzpts+5) )

c...Abscissa

      ztheta = (/ (REAL (i), i = 0, nzpts+4) /) / nzpts

c... Shift delta by half grid grid, ie. maybe pi

      zdeltas(1:nzpts1) =
     $     (/ zdelta(nzpts/2+1:nzpts),zdelta(1:nzpts/2+1) /)
      zzero = 0.0
      
      pllabel0 = "Delta: "//jobid
c      lena = index ( pllabel0, '  ', backl ) - 1
      lena = LEN_TRIM( pllabel0 )
      pllabel = pllabel0(1:lena)
      CALL pospl2s ( ztheta, zdeltas, zzero, 1,nzpts, 1,pllabel,nzpts )

      CALL maps(zxmin,zxmax, zymin,zymax,.114,.714, 0.12,0.720 )

      DO  i = 1, nzlin
         zphi0 = -zqa + (1.0+zqa)*(i-1)/(nzlin-1)  ! Origin for zphi.
         DO j = 1, nzpts1
            zphi(j) = zphi0 + zqa * ( ztheta(j) + zdeltas(j)/ztpi )
            IF ( zphi(j) < 0.0 ) zphi(j) = 0.0
            IF ( zphi(j) > 1.0 ) zphi(j) = 1.0
         END DO
         CALL trace ( ztheta, zphi, nzpts1, 1,1,0.,0. )
      END DO

      DEALLOCATE ( ztheta, zphi, zdeltas, zzero )

      xwrt = zxmax + (zxmax-zxmin)/40.0
      ywrt = zymax - (zymax-zymin)/6.0
      CALL setlch ( xwrt,ywrt,0,0,0,-1)
      WRITE (strg, 101 )
      CALL gtext(strg,-1,-1)
 101  FORMAT ( "Field lines" )

      xwrt = zxmin + (zxmax-zxmin)/40.0
      ywrt = zymin - (zymax-zymin)/6.0

      CALL setlch ( xwrt,ywrt, 0,1,0,-1 )

      WRITE ( stringv, 3001 ) zqa
 3001 FORMAT ( 1x, "Magnetic field lines:", 3x,
     $     "q at surface = ", f7.3 )
      CALL wrtstr ( stringv, 1 )

      CALL framep( jobid, ff )
c
      RETURN
      END SUBROUTINE B_LINES

c..................................................................
      subroutine orthchk ( a,nda1,nda2, b,ndb1,ndb2, zwght,
     $     c,ndc1,ndc2,
     $     la,lb,lc, ld1,ld2, label, nout1,nout2 )
c.................................................................
c
c... Checks orthomormality of matrices A and B.
c    Just ordinary sum: correspond to Trapezoidal rule.
c... zwght is the weight function for the orthonormality. 

      dimension a(nda1,nda2), b(ndb1,ndb2), c(ndc1,ndc2)
      REAL, DIMENSION(*) :: zwght
      character*(*) label
c
      do i = 1, la
         do j = 1, lc
            sum = 0.0
            do k = 1, lb
c$$$               sum = sum + a(k,i)* b(k,j) * zwght(k)
               sum = sum + a(k,i)* b(k,j)
            end do
            c(i,j) = sum / lb
         end do
      end do
c
       call matwrtn ( c,ndc1,ndc2, 1,1, la,lc, ld1,ld2,
     $     label, nout1,nout2 )
c
       if ( nout1 .ne. 0 )
     $ write ( nout1, '(10x,"Diagonal Elements = ",/,(1x,1p10e11.4))' )
     $      (c(i,i), i = 1, lc)
       if ( nout2 .ne. 0 )
     $ write ( nout2, '(10x,"Diagonal Elements = ",/,(1x,1p10e11.4))' )
     $      (c(i,i), i = 1, lc)
c
      return
      end

c..............................................
      subroutine plpgeom
c...................................................................
c
      include 'vacuum1.inc'
      include 'vacuum2.inc'
      include 'vacuum3.inc'
      include 'vacuum5.inc'
c
c$$$      call boundsi(xwal,zwal,1,mth,xmnp,xmxp,zmnp,zmxp,
c$$$     $     ixn,ixx,izn,izx)
c$$$c
c$$$      wlrad = 0.5 * ( xmxp - xmnp )
c$$$      wmaj = 0.5 * ( xmxp + xmnp )
c$$$      welong = (zmxp-zmnp) / (xmxp-xmnp)
c$$$      wdarg = ( wmaj - xwal(izx) ) / wlrad
c$$$c
c$$$      write ( iotty, '(/,"xmnp,xmxp,zmnp,zmxp,xwal(izx),wdarg= ",
c$$$     $     1p6e12.4,/)') xmnp,xmxp,zmnp,zmxp,xwal(izx),wdarg
c$$$c
c$$$      if ( abs(wdarg) .ge. 1.0000000 ) then
c$$$         write ( iotty, '(/,"WDARG = ", 1pe12.4)' ) wdarg
c$$$         wdarg = 0.9999999999
c$$$         go to 132
c$$$      end if
c$$$c
c$$$      wdness = asin ( wdarg )
c$$$c
c$$$ 132  continue
c$$$c
c$$$      write ( iotty, '(/,"wlrad,wmaj,welong,wdness= ",1p4e12.4/)')
c$$$     $     wlrad,wmaj,welong,wdness
c
      xmx = xmaj
ccc      xma = xmaj
      zma = 0.0
         write( 6,'("Circumference of Plasma,vac = ",1pe11.3)')
     $        pcircm
         write(23,'("Circumference of Plasma,vac = ",1pe11.3)')
     $        pcircm
      call drawc1(xinf,zinf,1,mth1,"z","x",xmx,xma,zma,
     .     plrad,1,1,idot, jobid, ishape,a,b,aw,bw,cw,dw,tw,
     $     abulg, bbulg, tbulg, wcentr, pcircm,
     $     pelong, pdness, wmaj,wlrad, welong, wdness )
c
c..... plot x,z, for plasma and wall versus theta.
c
      call atpoint ( "PLPGEOM-1", "mth", mth, "pcircm",pcircm,
     $     iotty, outmod )

      call pospl1 ( pigrd,xinf, 1,mth1, 1,"xpla",mth1 )
      call pospl1 ( pigrd,zinf, 1,mth1, 2,"zpla",mth1 )
      call framep(jobid,ff)

      call atpoint ( "PLPGEOM-2", "mth", mth, "plrad",plrad,
     $     iotty, outmod )

 50   continue
c
      return
      end
c     c
c...................................................................
      subroutine plgeom
c...................................................................
c
      include 'vacuum1.inc'
      include 'vacuum2.inc'
      include 'vacuum3.inc'
      include 'vacuum5.inc'
      include 'vacuum7.inc'
c
c     dimension xpp(nths),zpp(nths),ww1(nths),ww2(nths),
c     $     ww3(nths)
c
      call boundsi(xwal,zwal,1,mth,xmnp,xmxp,zmnp,zmxp,
     $     ixn,ixx,izn,izx)
c
      wlrad = 0.5 * ( xmxp - xmnp )
      wmaj =  0.5 * ( xmxp + xmnp )
      wzma =  0.5 * ( zmxp + zmnp )
      welong = (zmxp-zmnp) / (xmxp-xmnp)
      wdarg = ( wmaj - xwal(izx) ) / wlrad

      angwal = 0.0
      angwalq = 0.0

      write ( iotty, '(/,"xmnp,xmxp,zmnp,zmxp,xwal(izx),wdarg= ",
     $     1p6e12.4,/)') xmnp,xmxp,zmnp,zmxp,xwal(izx),wdarg
c
      if ( abs(wdarg) .ge. 1.0000000 ) then
         write ( iotty, '(/,"WDARG = ", 1pe12.4)' ) wdarg
         wdarg = 0.9999999999
         go to 132
      end if
c
      wdness = asin ( wdarg )
c
 132  continue
c
      write ( iotty, '(/,"wlrad, wmaj, wzma, welong, wdness= ",/,
     $     1p5e12.4/)')
     $     wlrad,wmaj,wzma,welong,wdness
      write ( outmod, '(/,"wlrad, wmaj, wzma, welong, wdness= ",/,
     $     1p5e12.4/)')
     $     wlrad,wmaj,wzma,welong,wdness

c
      xmx = xmaj
ccc   xma = xmaj
      zma = 0.0
      write( 6,'("Circumference of Plasma,vac = ",1pe11.3)')
     $     pcircm
      write(23,'("Circumference of Plasma,vac = ",1pe11.3)')
     $     pcircm

      DO i = 1, mth1
         angwal(i) =  - ATAN2 ( zwal(i) - wzma,xwal(i)-wmaj )
         IF ( angwal(i) < 0.0 ) angwal(i) = angwal(i) + twopi
         angwal(i) = 180.0 * angwal(i) / pye
      END DO

c.. Ensure monoticity of the angles
      DO i = 1, mth/8
         j = mth1-i+1
         IF ( angwal(i) > 90.0 ) angwal(i) = angwal(i) - 360.0
         IF ( angwal(j) < 90.0 ) angwal(j) = 360.0 - angwal(j)
      END DO

c$$$      WRITE ( outmod, '(/,1x,"Angwal = "/, (10es11.3) )' )
c$$$     $        (angwal(i), i = 1, mth1)

      angwalq = angwal

      IF ( leqarcw == 1 ) go to 43

      DO i = 1, mth1
         angwalq(i) =  - ATAN2 ( zwalq(i) - wzma,xwalq(i)-wmaj )
         IF ( angwalq(i) < 0.0 ) angwalq(i) = angwalq(i) + twopi
         angwalq(i) = 180.0 * angwalq(i) / pye
      END DO

c.. Ensure monoticity of the angles. Check a few points at end points.
      DO i = 1, mth/8
         j = mth1-i+1
         IF ( angwalq(i) > 90.0 ) angwalq(i) = angwalq(i) - 360.0
         IF ( angwalq(j) < 90.0 ) angwalq(j) = 360.0 - angwalq(j)
      END DO

      WRITE ( outmod, '(/,1x,"Angwalq = "/, (10es11.3) )' )
     $     (angwalq(i), i = 1, mth1)

 43   CONTINUE

      call drawc2(xwal,zwal,xinf,zinf,1,mth1,"z","x",xmx,xma,zma,
     $     plrad,1,1,idot, jobid, ishape,a,b,aw,bw,cw,dw,tw,
     $     abulg, bbulg, tbulg, wcentr, pcircm, wcircm,
     $     pelong, pdness, wmaj,wlrad, welong, wdness )
c
c$$$      if ( lfbcoil .eq. 1 ) then
c$$$         call drawc3(xwal,zwal,xinf,zinf,xcw,zcw,1,mth1,1,mtcoil1,
c$$$     $        "z","x",xmx,xma,zma,
c$$$     $        plrad,1,1,idot, jobid, ishape,a,b,aw,bw,cw,dw,tw,
c$$$     $        abulg, bbulg, tbulg, wcentr, pcircm,
c$$$     $        pelong, pdness, wmaj,wlrad, welong, wdness, chgt )
c$$$      end if
c
      if ( lplot1 .eq. 0 ) go to 50
c
c.....plot x,z, for plasma and wall versus theta.
c
      call pospl1 ( pigrd,xinf, 1,mth1, 1,"xpla",mth1 )
      call pospl1 ( pigrd,zinf, 1,mth1, 2,"zpla",mth1 )
      call pospl1 ( pigrd,xwal, 1,mth1, 3,"xwal",mth1 )
      call pospl1 ( pigrd,zwal, 1,mth1, 4,"zwal",mth1 )
      call framep(jobid,ff)
c
 50   continue
c
      return
      end
c
c...................................................................
      subroutine qgrid ( xin, zin, ell, thgr, mw, thlag )
c................................................................
c
c...  Redistribute points for equal arcs.
c    Inputs:   xin, zin, mw(no. of points on closed domain)
c..  Outputs:  Arc length on raw grid in ell.
c              Theta for equal arcs on raw grid in thgr
c    Dummy:    thlag.
c
      dimension xin(*), zin(*), ell(*), thgr(*), thlag(*)
c
c.....set up dummy array
c
      do  iw = 1, mw
         thlag(iw) = (1.0/(mw-1))*(iw-1)
      end do
c
c     get arc length on wall, ell, on raw grid.
c
      ell(1) = 0.0
      do  iw = 1, mw
         if ( iw .eq. 1 ) go to 15
         thet = ( thlag(iw)+thlag(iw-1) ) / 2.0
         call lag ( thlag,xin,mw,3,thet,f,df,1 )
         xtzt = df
         call lag ( thlag,zin,mw,3,thet,f,df,1 )
         xtzt = sqrt ( xtzt**2 + df**2 )
         ell(iw) = ell(iw-1) + xtzt/(mw-1)
 15      continue
      end do
c
c.... Get grid for equal arcs: Set thlag as a function of ell,
c..   then interpolate for equal arcs.
c
      do i = 1, mw
         elgr = ( ell(mw)/(mw-1) ) * (i-1)
         call lag ( ell,thlag,mw,3,elgr,f,df,0 )
         thgr(i) = f
      end do
c
      return
      end
c
c...................................................................
      subroutine qarc ( fin, qgrd, fout, mw, thlag )
c................................................................
c
c...  Redistribute points for equal arcs.
c
      dimension fin(*), fout(*), qgrd(*), thlag(*)
c
c.....set up dummy array
c
      do  iw = 1, mw
         thlag(iw) = (1.0/(mw-1))*(iw-1)
      end do
c
c......Now get fout for the equal arcs grid from fin.
c
      do  i = 1, mw
         ttt = qgrd(i)
         call lag ( thlag,fin,mw,3,ttt,f,df,0 )
         fout(i) = f
      end do
c
      return
      end
c
c...................................................................
      SUBROUTINE eqarcwx ( xin, zin, xout, zout, ell, thgr, mw )
c................................................................

c...  Redistribute points for equal arcs.

      REAL, DIMENSION(:), ALLOCATABLE :: thlag

      dimension xin(*), xout(*), zin(*), zout(*), ell(*), thgr(*)

      ALLOCATE ( thlag(mw+5) )

c.....set up dummy array

      do  iw = 1, mw
         thlag(iw) = (1.0/(mw-1))*(iw-1)
      end do

c     get arc length on wall, ell, on raw grid.

      ell(1) = 0.0
      do  iw = 1, mw
         if ( iw .eq. 1 ) go to 15
         thet = ( thlag(iw)+thlag(iw-1) ) / 2.0
         call lag ( thlag,xin,mw,3,thet,f,df,1 )
         xtzt = df
         call lag ( thlag,zin,mw,3,thet,f,df,1 )
         xtzt = sqrt ( xtzt**2 + df**2 )
         ell(iw) = ell(iw-1) + xtzt/(mw-1)
 15      continue
      end do

c.... Get grid for equal arcs: Set thlag as a function of ell,
c..   then interpolate for equal arcs.

      do i = 1, mw
         elgr = ( ell(mw)/(mw-1) ) * (i-1)
         call lag ( ell,thlag,mw,3,elgr,f,df,0 )
         thgr(i) = f
      end do

c......Now get xout for the equal arcs grid from xin. Same for zin.

      do  i = 1, mw
         ttt = thgr(i)
         call lag ( thlag,xin,mw,3,ttt,f,df,0 )
         xout(i) = f
         call lag ( thlag,zin,mw,3,ttt,f,df,0 )
         zout(i) = f
      end do

      DEALLOCATE ( thlag )

      call atpoint ( "END of EQARCW","mw", mw,"ell(mw)",ell(mw), 6, 23 )

      END SUBROUTINE eqarcwx

c...................................................................
      SUBROUTINE eqarc ( xin, zin, mw )
c................................................................
c
c...  Redistribute points for equal arcs. 
c
      REAL, DIMENSION(*) :: xin, zin
      REAL, DIMENSION(:), ALLOCATABLE :: xout, zout, ell, thgr, thlag

      CALL atpoint ( "Top of EQARCW","mw", mw,"xin(1)",xin(1), 6, 23 )

      mw5 = mw + 5
      ALLOCATE ( xout(mw5),zout(mw5), ell(mw5), thgr(mw5),thlag(mw5) )

c.....set up dummy array

      DO  iw = 1, mw
         thlag(iw) = (1.0/(mw-1))*(iw-1)
      END DO

c     get arc length on surface , ell, on raw grid.

      ell(1) = 0.0
      DO  iw = 1, mw
         IF ( iw == 1 ) GO TO 15
         thet = ( thlag(iw)+thlag(iw-1) ) / 2.0
         CALL lag ( thlag,xin,mw,3,thet,f,df,1 )
         xtzt = df
         CALL lag ( thlag,zin,mw,3,thet,f,df,1 )
         xtzt = sqrt ( xtzt**2 + df**2 )
         ell(iw) = ell(iw-1) + xtzt/(mw-1)
 15      CONTINUE
      END DO
      CALL atpoint ( "Mid of EQARCW","mw", mw,"ell(mw)",ell(mw), 6, 23 )

c.... Get grid for equal arcs: Set thlag as a function of ell,
c..   then interpolate for equal arcs.

      DO i = 1, mw
         elgr = ( ell(mw)/(mw-1) ) * (i-1)
         CALL lag ( ell,thlag,mw,3,elgr,f,df,0 )
         thgr(i) = f
      END DO
      CALL atpoint ( "Mid2 of EQARCW","mw", mw,"ell(mw)",ell(mw), 6,23 )

c......Now get xout for the equal arcs grid from xin. Same for zin.

      DO  i = 1, mw
         ttt = thgr(i)
         CALL lag ( thlag,xin,mw,3,ttt,f,df,0 )
         xout(i) = f
         CALL lag ( thlag,zin,mw,3,ttt,f,df,0 )
         zout(i) = f
      END DO

      xin(1:mw) = xout(1:mw)
      zin(1:mw) = zout(1:mw)

      CALL atpoint ( "END of EQARCW","mw", mw,"ell(mw)",ell(mw), 6, 23 )

      DEALLOCATE ( xout, zout, ell, thgr, thlag )

      END SUBROUTINE eqarc

c...................................................................
      subroutine eqarcw ( xin, zin, xout, zout, ell, thgr, thlag, mw )
c................................................................
c
c...  Redistribute points for equal arcs.
c
      dimension xin(*), xout(*), zin(*), zout(*),
     $     ell(*), thgr(*), thlag(*)
c
      call atpoint ( "Top of EQARCW","mw", mw,"xin(1)",xin(1), 6, 23 )
c
c.....set up dummy array
c
      do  iw = 1, mw
         thlag(iw) = (1.0/(mw-1))*(iw-1)
      end do
c
c     get arc length on wall, ell, on raw grid.
c
      ell(1) = 0.0
      do  iw = 1, mw
         if ( iw .eq. 1 ) go to 15
         thet = ( thlag(iw)+thlag(iw-1) ) / 2.0
         call lag ( thlag,xin,mw,3,thet,f,df,1 )
         xtzt = df
         call lag ( thlag,zin,mw,3,thet,f,df,1 )
         xtzt = sqrt ( xtzt**2 + df**2 )
         ell(iw) = ell(iw-1) + xtzt/(mw-1)
 15      continue
      end do
      call atpoint ( "Mid of EQARCW","mw", mw,"ell(mw)",ell(mw), 6, 23 )
c
c
c.... Get grid for equal arcs: Set thlag as a function of ell,
c..   then interpolate for equal arcs.
c
      do i = 1, mw
         elgr = ( ell(mw)/(mw-1) ) * (i-1)
         call lag ( ell,thlag,mw,3,elgr,f,df,0 )
         thgr(i) = f
      end do
      call atpoint ( "Mid2 of EQARCW","mw", mw,"ell(mw)",ell(mw), 6,23 )
c
c......Now get xout for the equal arcs grid from xin. Same for zin.
c
      do  i = 1, mw
         ttt = thgr(i)
         call lag ( thlag,xin,mw,3,ttt,f,df,0 )
         xout(i) = f
         call lag ( thlag,zin,mw,3,ttt,f,df,0 )
         zout(i) = f
      end do
c
      call atpoint ( "END of EQARCW","mw", mw,"ell(mw)",ell(mw), 6, 23 )
c
      return
      end

c...................................................................
      subroutine eqarcw1 ( xin, zin, mwin, xout, zout, mwout,
     $     ell, thgr, thlag )
c................................................................
c
c...  Redistribute points for equal arcs.
c     Input and output grids may differ: mwin and mwout on closed domain.
c
      dimension xin(*), xout(*), zin(*), zout(*),
     $     ell(*), thgr(*), thlag(*)
      call atpoint ( "Top of EQARCW","mwin", mwin,
     $     "ell(mwin)",ell(mwin), 6, 23 )
c
c.....set up dummy array
c
      do  iw = 1, mwin
         thlag(iw) = (1.0/(mwin-1))*(iw-1)
      end do
c
c     get arc length on wall, ell, on raw grid.
c
      ell(1) = 0.0
      do  iw = 1, mwin
         if ( iw .eq. 1 ) go to 15
         thet = ( thlag(iw)+thlag(iw-1) ) / 2.0
         call lag ( thlag,xin,mwin,3,thet,f,df,1 )
         xtzt = df
         call lag ( thlag,zin,mwin,3,thet,f,df,1 )
         xtzt = sqrt ( xtzt**2 + df**2 )
         ell(iw) = ell(iw-1) + xtzt/(mwin-1)
 15      continue
      end do
      call atpoint ( "Mid of EQARCW","mwin", mwin,
     $     "ell(mwin)",ell(mwin), 6, 23 )
c
c
c.... Get grid for equal arcs: Set thlag as a function of ell,
c..   then interpolate for equal arcs on a mwout grid.
c
      do i = 1, mwout
         elgr = ( ell(mwin)/(mwout-1) ) * (i-1)
         call lag ( ell,thlag,mwin,3,elgr,f,df,0 )
         thgr(i) = f
      end do
      call atpoint ( "Mid2 of EQARCW","mwout", mwout,
     $     "ell(mwin)",ell(mwin), 6,23 )
c
c......Now get xout for the equal arcs grid from xin. Same for zin.
c
      do  i = 1, mwout
         ttt = thgr(i)
         call lag ( thlag,xin,mwin,3,ttt,f,df,0 )
         xout(i) = f
         call lag ( thlag,zin,mwin,3,ttt,f,df,0 )
         zout(i) = f
      end do
c
      call atpoint ( "END of EQARCW","mwout", mwout,
     $     "ell(mwin)",ell(mwin), 6, 23 )
c
      return
      end


c......................................................
      subroutine atpoint ( point, iname, ivar, aname, avar,
     $     nout1,nout2 )
c.....................................................
c
      character*(*) point, iname, aname
c
      if ( nout1 .ne. 0 )
     $     write ( nout1, '(/,
     $     "ooooooooooooooooooooooooooooooooooooooooooooooooo",
     $     /,5x, "Reached: ", a, ", "a, " = ", i6,
     $     ", ", a," = ", 1pe12.4, /
     $     "ooooooooooooooooooooooooooooooooooooooooooooooooo",
     $     /)' ) point, iname, ivar,
     $     aname, avar
      if ( nout2 .ne. 0 )
     $     write ( nout2, '(/,
     $     "ooooooooooooooooooooooooooooooooooooooooooooooooo",
     $     /,5x, "Reached: ", a, ", "a, " = ", i6,
     $     ", ", a," = ", 1pe12.4, /
     $     "ooooooooooooooooooooooooooooooooooooooooooooooooo",
     $     /)' ) point, iname, ivar,
     $     aname, avar
c
      return
      end
c.....................................................
