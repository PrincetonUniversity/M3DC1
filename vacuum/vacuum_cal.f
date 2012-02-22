c
c     3.2          solution of the vacuum integral equations.
c
      SUBROUTINE vaccal
c
      include 'vacuum1.inc'
      include 'vacuum2.inc'
      include 'vacuum3.inc'
      include 'vacuum5.inc'
      include 'vacuum7.inc'
      include 'vacuum8.inc'

      COMMON / grvac / aval0, phex0, zphfc0, iscs, isco,
     $     j1com,j2com

      CHARACTER(132) :: zlabel
      CHARACTER(24) :: datimv
      LOGICAL backl
      COMPLEX ziii

      dimension xpla(nths), zpla(nths)
c      dimension gwpi(nths,nths)
c
      REAL, DIMENSION(:), ALLOCATABLE :: zwk, zwkpl
      REAL, DIMENSION(:,:), ALLOCATABLE :: gwpi, v_cpwwp, grrin
      REAL, DIMENSION(:,:), ALLOCATABLE :: wrkvr2
      REAL, DIMENSION(:,:), ALLOCATABLE :: ajll, rmatr, rmati, umtrx
      REAL, DIMENSION(:,:), ALLOCATABLE :: workg1, workg2
c      REAL, DIMENSION(:,:), ALLOCATABLE :: akker
 
      REAL, DIMENSION(:), ALLOCATABLE :: zgri, zgro, zgrop, zdlent
      COMPLEX, DIMENSION(:,:), ALLOCATABLE :: zgrij0, zgrij, zgrijb,
     $     zgrchi, zgrbth, zijbp, zgrbph, cgrri, expmilt
      REAL, DIMENSION(:,:), ALLOCATABLE :: zgric, zgris, zgrcc, zgrss,
     $     zgrll 
c
c     dimension grgr(nths,nfm), grgi(nths,nfm)
      dimension arr(nfm,nfm), aii(nfm,nfm), ari(nfm,nfm), air(nfm,nfm)
      dimension vacmti(nfm,nfm), summ(2)
c      dimension gatovac(nfm,nfm)
c      dimension umtrx(nfm,nfm)
      dimension wkxx(nths), wkyy(nths)
c
      real nq
      integer tmth
c
      common / bigv1 / grdgre(nths2,nths2), grwp(nths,nths),
     $     grri(nths2,nfm2)
c
c     parameter(nfmsq=nfm*nfm,neqv1=1,neqv2=neqv1+nfmsq,
c     .  neqv3=neqv2+nfmsq,neqv4=neqv3+nfmsq,neqv5=neqv4+nfmsq)
      equivalence (xpla,xinf), (zpla,zinf)
c     equivalence ( grdgre(1,1), vacmti(1,1) )
c     equivalence ( grdgre(neqv2), arr(1,1) )
c     equivalence ( grdgre(neqv3), aii(1,1) )
c     equivalence ( grdgre(neqv4), air(1,1) )
c     equivalence ( grdgre(neqv5), ari(1,1) )
c
      parameter ( neqv1=1, neqv2=(nfm-1)/2+1,neqv3=nfm+1,
     .     neqv4=3*(nfm+1)/2,neqv5=2*nfm+1 )
c
c      equivalence ( grdgre(1,1), vacmti(1,1) )
      equivalence ( grdgre(1,1), arr(1,1) )
      equivalence ( grdgre(npfg2x,npfg12), aii(1,1) )
      equivalence ( grdgre(1,npfg13), air(1,1) )
      equivalence ( grdgre(npfg2x,npfg14), ari(1,1) )
c
      dimension vacpstr(nfmsq),vacpsti(nfmsq)
c
c     equivalence ( vacmti(1,1),vacpst(1) )

      ALLOCATE ( gwpi(mth2,mth2), rmatr(nfm,nfm), rmati(nfm,nfm) )
c
      call atpoint ( "VACCAL", "lrnge",lrnge,"a", a, iotty, outmod )
c
c........Vacmat-like quantities will be multiplied by twopi to compensate
c     for the division of bval(=gren)  by twopi in kernel..
c
      write ( outmod, '(/,10x, "Matrix Storage: K(obs_ji,sou_ji):",//,
     $     10x, "j = observer points.  :i = source points.",/,
     $     10x, "  ie. K operates on chi from the left."//,
     $     10x, "Observer Source  Block",/,
     $     10x, " plasma  plasma   1  1",/,
     $     10x, " plasma  wall     1  2",/,
     $     10x, " wall    plasma   2  1",/,
     $     10x, " wall    wall     2  2",/ )' )
c
      backl = .false.

      datimv(1:24) = date_array(1:24)
c      datimv = datev(1:10)//","//timev(1:10)
      lenj = index(jobid, "  ", backl ) - 1
      ier = 0
      factpi = twopi

      jmax = 2*lmax(1) + 1
      jmax1 = lmax(1) - lmin(1) + 1
      jmax2 = 2*jmax1
      lmax1 = lmax(1) + 1
c
      lmax2 = 2*jmax1
c      if ( lfele .ne. 0 ) lmax2 = jmax1

      ln = lmin(1)
      lx = lmax(1)
      mfel2 = 2*mfel
      mfjm2 = max ( mfel2,jmax2 )

      jdel = 8
      q = qa1
      nq = n*q
c
      tmth = 2*mth
      mthsq = tmth * tmth
      lmth = tmth * 2*jmax1
c
c...  Use Fourier on the shell for now, with same as plasma:
c      jwal1 = jmax1
c      jwal2 = 2 * jwal1
c
      write ( iotty, '(/,"lreshel, jwal1 = ", 2i4 )' ) lreshel, jwal1
      write ( outmod, '(/,"lreshel, jwal1 = ", 2i4 )' ) lreshel, jwal1

      farwal = .false.
      if ( (a .ge. 10.0) .or. (lspark .ne. 0) ) farwal = .true.

      iwp = 2
      IF ( farwal ) iwp = 1

      IF ( lhighn == 2 ) THEN
         iwp = 4
         IF ( farwal ) iwp = 2
      END IF

      mth12 = iwp*mth
      nouk = mth12   ! Number of unknowns in the Gaussian Elimination.

c.....j1v and j2v are input to vacout and are the size of the matrices.
c
      j1v = nfm
      j2v = nfm

c... Allocate interim storage for the K and G kernels to accomodate complex
c    elements for the high-n calculaions.

c      ALLOCATE ( akker(iwp*(mfel+2),iwp*(mfel+2)) )
c     $     agker(4*(mfel+2),2*(mfel+2)) )

c.. Use grrin for the G kernel:

      malgr = iwp*mth+10
      jalgr = lmax2 + jwal2 + 10
      IF ( lhighn == 2 ) THEN
         malgr = iwp*mfel + 10
         jalgr = 2*mfel
      END IF
      IF ( ldqdtw == 2 ) THEN
         malgr = iwp * mth
         jalgr = iwp * mth
      END IF

      ALLOCATE ( grrin(malgr,jalgr) )

c... Initialize Matrices

      grri = zero
      gwin = zero
      gwot = zero
      grdgre = zero
      akker = 0.0
      grrin = 0.0

      rmatr = 0.0
      rmati = 0.0
c

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      if ( check1 )
     .     call timer ( outmod,iotty, "before kernels" )

c... Feedback coils, A and B:
c
      if ( lfbcoila .eq. 0 ) go to 8001
c
c..------------------------------------------------------------
c                          Plasma-Coila
c..------------------------------------------------------------
      call atpoint ( "Plasma-Coila","mtcoila", mtcoila,"xcwa(1)",
     $     xcwa(1), iotty, outmod )
c
c...       Use off appropriate diagonal block of GRDGRE as storage.
c            Don't need GREEN (a.k.a grwp) here.
c
      j1 = 1
      j2 = 2
      ksgn = 2*j2 - 3
      iopw = 0
      iops = 0

      j1a = (j1-1)*mth + 1
      j1b = j1*mth
      j2a = (j2-1)*mth + 1
      j2b = (j2-1)*mth + mtcoila

      call kernel0(xpla,zpla,xcwa,zcwa,grdgre,grwp,j1,j2,ksgn,
     $     iopw,iops,lfele, nths2,nths )
c
c..       Can integrate over the coil, used as source later:
c
      skcoilpa = 0.0
         skcoilpa(1:mth) = gainca * matmul (
     $        grdgre(j1a:j1b,j2a:j2b),ccoila(1:mtcoila) )
c
      call vecwrt ( mth, skcoilpa, "skcintpa", 1, mth, outmod, iotty )
c
c..------------------------------------------------------------
c                          Wall-Coila
c..------------------------------------------------------------
      call atpoint ( "Wall-Coila","mtcoila", mtcoila,"xcwa(1)",
     $     xcwa(1), iotty, outmod )
c
c...       Use off appropriate diagonal block of GRDGRE as storage.
c            Don't need GREEN (a.k.a grwp) here.
c
      j1 = 2
      j2 = 1
      ksgn = 2*j2 - 3
      iopw = 0
      iops = 0

      j1a = (j1-1)*mth + 1
      j1b = j1*mth
      j2a = (j2-1)*mth + 1
      j2b = (j2-1)*mth + mtcoila

      call kernel0(xwal,zwal,xcwa,zcwa,grdgre,grwp,j1,j2,ksgn,
     $     iopw,iops,lfele,nths2,nths)
c
c..        Integrate over the coil, used as source later:
c
      skcoilwa = 0.0
         skcoilwa(1:mth) = gainca * matmul (
     $        grdgre(j1a:j1b,j2a:j2b),ccoila(1:mtcoila) )
c
      call vecwrt ( mth, skcoilwa, "skcintwa", 1, mth, outmod, iotty )
c

 8001 continue

c
      if ( lfbcoilb .eq. 0 ) go to 8002
c
c..------------------------------------------------------------
c                          Plasma-Coilb
c..------------------------------------------------------------
      call atpoint ( "Plasma-Coilb","mtcoilb", mtcoilb,"xcwb(1)",
     $     xcwb(1), iotty, outmod )
c
c...       Use off appropriate diagonal block of GRDGRE as storage.
c            Don't need GREEN (a.k.a grwp) here.
c
      j1 = 1
      j2 = 2
      ksgn = 2*j2 - 3
      iopw = 0
      iops = 0

      j1a = (j1-1)*mth + 1
      j1b = j1*mth
      j2a = (j2-1)*mth + 1
      j2b = (j2-1)*mth + mtcoilb

      call kernel0(xpla,zpla,xcwb,zcwb,grdgre,grwp,j1,j2,ksgn,
     $     iopw,iops,lfele,nths2,nths)
c
c..        Integrate over the coil, used as source later:
c
      skcoilpb = 0.0
         skcoilpb(1:mth) = gaincb * matmul (
     $        grdgre(j1a:j1b,j2a:j2b),ccoilb(1:mtcoilb) )
c
      call vecwrt ( mth, skcoilpb, "skcintpb", 1, mth, outmod, iotty )

c..------------------------------------------------------------
c                          Wall-Coilb
c..------------------------------------------------------------
      call atpoint ( "Wall-Coilb","mtcoilb", mtcoilb,"xcwb(1)",
     $     xcwb(1), iotty, outmod )
c
c...       Use off appropriate diagonal block of GRDGRE as storage.
c            Don't need GREEN (a.k.a grwp) here.
c
      j1 = 2
      j2 = 1
      ksgn = 2*j2 - 3
      iopw = 0
      iops = 0

      j1a = (j1-1)*mth + 1
      j1b = j1*mth
      j2a = (j2-1)*mth + 1
      j2b = (j2-1)*mth + mtcoilb

      call kernel0(xwal,zwal,xcwb,zcwb,grdgre,grwp,j1,j2,ksgn,
     $     iopw,iops,lfele,nths2,nths)
c
c..        Integrate over the coil, used as source later:
c
      skcoilwb = 0.0
         skcoilwb(1:mth) = gaincb * matmul (
     $        grdgre(j1a:j1b,j2a:j2b),ccoilb(1:mtcoilb) )
c
      call vecwrt ( mth, skcoilwb, "skcintwb", 1, mth, outmod, iotty )
c
 8002 continue

c.. ------------------------------------------------------------
c                           Plasma-Plasma
c.. ------------------------------------------------------------
      call atpoint ( "Plasma-Plasma","jwal1", jwal1,"xinf(1)", xinf(1),
     $     iotty, outmod )

      j1 = 1
      j2 = 1

c... Pretend that the plasma is the wall for the interior problem

      IF ( linterior == 1 ) THEN
         j1 = 2
         j2 = 2
      END IF

c  The following indices will shift the elements of the kernels the
c  correct amount for storage.

      j1f = 2*(j1-1) * mfel
      j2f = 2*(j2-1) * mfel
      j1f1 = j1f + 1
      j2f1 = j2f + 1
      j1ff = j1f + mfel
      j2ff = j2f + mfel
      j1ff1 = j1ff + 1
      j2ff1 = j2ff + 1
      j1fff = j1ff + mfel
      j2fff = j2ff + mfel

      ksgn = 2*j2 - 3
      iopw = 1
      iops = 1

c..  These DO LOOPS will cover the real and imaginary parts of
c    observer and sources. ISCO, ISCS are passed in COMMON GRVAC
c..  GRRIN used for temporary storage here.

      nzobsr = 1
      nzsrce = 1

      IF ( lhighn == 2 ) THEN
         nzobsr = 1
         nzsrce = 2
      END IF

      DO izobsr = 1, nzobsr
         isco = izobsr
         DO izsrce = 1, nzsrce
            iscs = izsrce
            CALL kernel0(xpla,zpla,xpla,zpla,grdgre,grwp,j1,j2,ksgn,
     $           iopw,iops,lfele,nths2,nths)
            IF ( lhighn == 2 ) THEN
               DO i = 1, mfel
                  DO j = 1, mfel
                     indx1 = ( 2*(j1-1)+(isco-1) ) * mfel + i
                     indx2 = ( 2*(j2-1)+(iscs-1) ) * mfel + j
                     grrin(indx1,indx2) = grwp(i,j)
                  END DO
               END DO
            END IF
         END DO
      END DO

c$$$      call matwrtn ( grdgre,nths2,nths2,1,1,mth12,mth12,16,8,
c$$$     $        "grdgre before repack", outmod, iotty )
c$$$
c$$$      CALL matwrts ( grdgre,nths2,nths2,1,1,tmth,tmth,16,8,
c$$$     $        mth1,mth1,"grdgre at 2,2", outmod, iotty )
c$$$
c$$$      write ( iotty,  '(/,
c$$$     $     "Sum of first COLUMN in each block of GRDGRE:",/)')
c$$$      write ( outmod, '(/,
c$$$     $     "Sum of first COLUMN in each block of GRDGRE:",/)')
c$$$      do j1 = 1,2
c$$$         do j2 = 1,2
c$$$            sumg = 0.0
c$$$            do i = 1, mth
c$$$               sumg = sumg + grdgre( (j1-1)*mth+i,(j2-1)*mth+1 )
c$$$            end do
c$$$            write( iotty,  '("j1,j2, sumg= ",2i3,1pe11.3)' ) j1,j2,sumg
c$$$            write( outmod, '("j1,j2, sumg= ",2i3,1pe11.3)' ) j1,j2,sumg
c$$$         end do
c$$$      end do
c$$$c
c$$$c
c$$$      write ( iotty,  '(/,
c$$$     $     "Sum of first ROW in each block of GRDGRE:",/)')
c$$$      write ( outmod, '(/,
c$$$     $     "Sum of first ROW in each block of GRDGRE:",/)')
c$$$      do j1 = 1,2
c$$$         do j2 = 1,2
c$$$            sumg = 0.0
c$$$            do i = 1, mth
c$$$               sumg = sumg + grdgre( (j1-1)*mth+1,(j2-1)*mth+i )
c$$$            end do
c$$$            write( iotty,  '("j1,j2, sumg= ",2i3,1pe11.3)' ) j1,j2,sumg
c$$$            write( outmod, '("j1,j2, sumg= ",2i3,1pe11.3)' ) j1,j2,sumg
c$$$         end do
c$$$      end do
c$$$c


c... Repack the K kernel from the (2,2) block to the (1,1) block for the 
c      interior problem.
      IF ( linterior == 1 ) THEN
         grdgre ( 1:mth, 1:mth ) = grdgre ( mth+1:2*mth, mth+1:2*mth )
      END IF


c$$$      call matwrtn ( grdgre,nths2,nths2,1,1,mth12,mth12,16,8,
c$$$     $        "grdgre after repack", outmod, iotty )
c$$$
c$$$      IF ( .NOT. farwal )
c$$$     $     CALL matwrts ( grdgre,nths2,nths2,1,1,tmth,tmth,16,8,
c$$$     $        mth1,mth1,"grdgre at 2,2", outmod, iotty )
c$$$
c$$$      write ( iotty,  '(/,
c$$$     $     "Sum of first COLUMN in each block of GRDGRE:",/)')
c$$$      write ( outmod, '(/,
c$$$     $     "Sum of first COLUMN in each block of GRDGRE:",/)')
c$$$      do j1 = 1,2
c$$$         do j2 = 1,2
c$$$            sumg = 0.0
c$$$            do i = 1, mth
c$$$               sumg = sumg + grdgre( (j1-1)*mth+i,(j2-1)*mth+1 )
c$$$            end do
c$$$            write( iotty,  '("j1,j2, sumg= ",2i3,1pe11.3)' ) j1,j2,sumg
c$$$            write( outmod, '("j1,j2, sumg= ",2i3,1pe11.3)' ) j1,j2,sumg
c$$$         end do
c$$$      end do
c$$$c
c$$$c
c$$$      write ( iotty,  '(/,
c$$$     $     "Sum of first ROW in each block of GRDGRE:",/)')
c$$$      write ( outmod, '(/,
c$$$     $     "Sum of first ROW in each block of GRDGRE:",/)')
c$$$      do j1 = 1,2
c$$$         do j2 = 1,2
c$$$            sumg = 0.0
c$$$            do i = 1, mth
c$$$               sumg = sumg + grdgre( (j1-1)*mth+1,(j2-1)*mth+i )
c$$$            end do
c$$$            write( iotty,  '("j1,j2, sumg= ",2i3,1pe11.3)' ) j1,j2,sumg
c$$$            write( outmod, '("j1,j2, sumg= ",2i3,1pe11.3)' ) j1,j2,sumg
c$$$         end do
c$$$      end do
c

      IF ( lhighn /= 2) GO TO 31
c$$$c.. Now put the K kernel in the correct position. Store CS submatrix
c$$$c   in AKKER
c$$$
c$$$      akker(j1f1:j1ff,j2f1:j2ff) = grdgre(j1f1:j1ff,j2f1:j2ff) +
c$$$     $     grdgre(j1ff1:j1fff,j2ff1:j2fff)
c$$$      akker(j1ff1:j1fff,j2ff1:j2fff) =  akker(j1f1:j1ff,j2f1:j2ff)
c$$$      akker(j1f1:j1ff,j2ff1:j2fff) =
c$$$     $     grdgre(j1ff1:j1fff,j2f1:j2ff) - grdgre(j1f1:j1ff,j2ff1:j2fff)
c$$$      akker(j1ff1:j1fff,j2f1:j2ff) = - akker(j1f1:j1ff,j2ff1:j2fff)
c$$$
c. Now do the same for the G kernel. Use GRRI directly. Will re-store
c  in GRRIN for the Gaussian later.
c  Normalize to the SQRT(mfel) here also. So no need to call FELANV3

c$$$      grri(j1f1:j1ff,j2f1:j2ff) = grrin(j1f1:j1ff,j2f1:j2ff) +
c$$$     $     grrin(j1ff1:j1fff,j2ff1:j2fff)
      grri(j1f1:j1ff,j2f1:j2ff) =     grrin(j1f1:j1ff,j2f1:j2ff)
      grri(j1ff1:j1fff,j2ff1:j2fff) = grrin(j1f1:j1ff,j2f1:j2ff)
c$$$      grri(j1f1:j1ff,j2ff1:j2fff) =
c$$$     $     grrin(j1ff1:j1fff,j2f1:j2ff) - grrin(j1f1:j1ff,j2ff1:j2fff)
      grri(j1f1:j1ff,j2ff1:j2fff) =   grrin(j1f1:j1ff,j2ff1:j2fff)
      grri(j1ff1:j1fff,j2f1:j2ff) = - grrin(j1f1:j1ff,j2ff1:j2fff)

      znorm = SQRT( FLOAT(mfel) )
      grri(j1f1:j1fff,j2f1:j2fff) = znorm * grri(j1f1:j1fff,j2f1:j2fff)

 31   CONTINUE

c$$$      call matwrts ( grdgre,nths2,nths2,1,1,mth,mth,8,8,
c$$$     $        1,1,"grdgre at 1,1", outmod, iotty )
c$$$      call matwrts ( grdgre,nths2,nths2,j1f1,j2f1,j1fff,j2fff,8,8,
c$$$     $        j1f1,j2f1,"grdgre at 1,1,ff", outmod, iotty )
c$$$      call matwrtn ( grrin,malgr,jalgr,1,1,nouk,lmax2,16,8,
c$$$     $        "GRRIN AT 1 1", outmod, iotty )
c$$$      call matwrtn ( grri,nths2,nfm2,1,1,nouk,lmax2,16,8,
c$$$     $        "GRRI AT 1 1", outmod, iotty )

      if ( checkd ) then
      call matwrtn ( grwp,nths,nths,1,1,mth,mth,16,8,
     $        "grwp at p-p", outmod, iotty )
      call matwrtn ( grdgre,nths2,nths2,1,1,mth,mth,mth,mth,
     $        "grdgre at 1,1", outmod, iotty )
      end if
c
c Store G(mth1) = G(1)

      grwp(mth1:mth2,1:mth) = grwp(1:2,1:mth)
      grwp(1:mth2,mth1:mth2) = grwp(1:mth2,1:2)

c... Set Imaginary grwp = 0 until we put in exp(-in\nu) in gren
      gwpi = 0.0
c
c.....fourier analyse source points. real and imag. parts.
c     pack into one-d array. ( ... not for leqt1f or leqt2f... )
c
      IF ( lfele .eq. 1 ) THEN
         CALL felang ( grwp, grri, cnqd, 0,0 )
         CALL felang ( grwp, grri, snqd, 0,jmax1 )
      END IF

      IF ( lfele .eq. 3 ) THEN
         CALL atpoint ( "felang3a","j1", j1,"xinf(1)", xinf(1),
     $        iotty, outmod )
         IF ( lhighn /= 2 ) THEN    ! already normalized for lhighn=2
            CALL felang3 ( grwp,nths,nths, grri, 0,0 )
            CALL felang3 ( gwpi,mth2,mth2, grri, 0,jmax1 )
         END IF
      END IF

      IF ( ldqdtw == 2 ) THEN
         grri(1:mth2,1:mth2) = grwp(1:mth2,1:mth2)
         nrhs = mth
         GO TO 32
      END IF
      IF ( lfour .eq. 1 ) THEN
         CALL fouran ( grwp, grri, cslth, 0,0 )
         CALL fouran ( grwp, grri, snlth, 0,jmax1 )
      END IF

 32   CONTINUE

      call atpoint ( "After Fouran","jwal1", jwal1,"xinf(1)", xinf(1),
     $     iotty, outmod )
c
c Get Fourier representation of K(1,1). Use grwp as storage:
c
      if ( checkd .and. (lfele .eq. 0) ) then
         call foura2 ( grdgre,0,0, grwp, 0 )
         call fanal1d ( grwp,nths,nths,0, wrkvr,wrkvi,nfm,nfm )
c
         call matwrtn(wrkvr,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Real Kpp(l,l)", outmod,0 )
         call matwrtn(wrkvi,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Imag Kpp(l,l)", outmod,0 )
c
c......fourier analyse gpp matrix....
c
         call fanal1d ( grri,nths2,nfm2,0, wrkvr,wrkvi,nfm,nfm )
c
         call matwrtn(wrkvr,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Real GPP_ll", outmod,0 )
         call matwrtn(wrkvi,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Imag GPP_ll", outmod,0 )
c
      end if

      IF ( farwal ) GO TO 34   !!! About line 871
c
c$$$c.....get wall coordinates....
c$$$c
c$$$      call wwall ( mth1, xwal, zwal )
c
c.. ------------------------------------------------------------
c                           Plasma-Wall
c.. ------------------------------------------------------------
      call atpoint ( "Plasma-Wall","jwal1", jwal1,"xinf(1)", xinf(1),
     $     iotty, outmod )
c
      j1 = 1
      j2 = 2

c  The following indices will shift the elements of the kernels the
c  correct amount for storage.

c$$$      j1f = 2*(j1-1) * mfel
c$$$      j2f = 2*(j2-1) * mfel
c$$$      j1f1 = j1f + 1
c$$$      j2f1 = j2f + 1
c$$$      j1ff = j1f + mfel
c$$$      j2ff = j2f + mfel
c$$$      j1ff1 = j1ff + 1
c$$$      j2ff1 = j2ff + 1
c$$$      j1fff = j1ff + mfel
c$$$      j2fff = j2ff + mfel

      ksgn = 2*j2 - 3
      iopw = 0
      if ( lreshel .eq. 1 ) iopw = 1
      iops = 0

      DO izobsr = 1, nzobsr
         isco = izobsr
         DO izsrce = 1, nzsrce
            iscs = izsrce
            CALL kernel0(xpla,zpla,xwal,zwal,grdgre,grwp,j1,j2,ksgn,
     $           iopw,iops,lfele,nths2,nths)
         END DO
      END DO
c$$$      call matwrts ( grdgre,nths2,nths2,1,mth1,mth,tmth,8,8,
c$$$     $        1,mth1,"grdgre at 1,2", outmod, iotty )
c$$$      call matwrts ( grdgre,nths2,nths2,j1f1,j2f1,j1fff,j2fff,8,8,
c$$$     $        j1f1,j2f1,"grdgre at 1,2,,ff", outmod, iotty )
c$$$      call matwrtn ( grri,nths2,nfm2,1,1,nouk,lmax2,16,8,
c$$$     $        "GRRI AT 1 2", outmod, iotty )

c$$$      IF ( lhighn /= 2) GO TO 41
c$$$c.. Now put the K kernel in the correct position. Store CS submatrix
c$$$c   in AKKER
c$$$
c$$$      akker(j1f1:j1ff,j2f1:j2ff) = grdgre(j1f1:j1ff,j2f1:j2ff) +
c$$$     $     grdgre(j1ff1:j1fff,j2ff1:j2fff)
c$$$      akker(j1ff1:j1fff,j2ff1:j2fff) =  akker(j1f1:j1ff,j2f1:j2ff)
c$$$      akker(j1f1:j1ff,j2ff1:j2fff) =
c$$$     $     grdgre(j1ff1:j1fff,j2f1:j2ff) - grdgre(j1f1:j1ff,j2ff1:j2fff)
c$$$      akker(j1ff1:j1fff,j2f1:j2ff) = - akker(j1f1:j1ff,j2ff1:j2fff)
c$$$
c$$$ 41   CONTINUE

c$$$            call kernel0(xpla,zpla,xwal,zwal,grdgre,grwp,j1,j2,ksgn,
c$$$     $           iopw,iops,lfele,nths2,nths)

c$$$      call matwrtn ( grwp,nths,nths,1,1,mth,mth,16,8,
c$$$     $        "grwp at p-w", outmod, iotty )

      IF ( lreshel .eq. 1 ) THEN
c
c Store G(mth1) = G(1)
         grwp(1:mth2,mth1:mth2) = grwp(1:mth2,1:2)

c Fourier analyze source index: no dth here. It is consistent with the LHS
c  of the Gaussian elimination.

         IF ( lfour .eq. 1 ) THEN
            gwin(1:mth,1:jwal1) =
     $           matmul( grwp(1:mth,1:mth),rwvce(1:mth,1:jwal1) )
            gwin(1:mth,jwal1+1:jwal2) =
     $           matmul( grwp(1:mth,1:mth),rwvco(1:mth,1:jwal1) )
c$$$            call fouran ( grwp, gwin, rwvce, 0,0 )
c$$$            call fouran ( grwp, gwin, rwvco, 0,jwal1 )
         END IF

      END IF
c
c.. ------------------------------------------------------------
c                           Wall-Wall -  minus
c.. ------------------------------------------------------------
      call atpoint ( "Wall-Wall","jwal1", jwal1,"xinf(1)", xinf(1),
     $     iotty, outmod )
      j1 = 2
      j2 = 2

c$$$c  The following indices will shift the elements of the kernels the
c$$$c  correct amount for storage.
c$$$
c$$$      j1f = 2*(j1-1) * mfel
c$$$      j2f = 2*(j2-1) * mfel
c$$$      j1f1 = j1f + 1
c$$$      j2f1 = j2f + 1
c$$$      j1ff = j1f + mfel
c$$$      j2ff = j2f + mfel
c$$$      j1ff1 = j1ff + 1
c$$$      j2ff1 = j2ff + 1
c$$$      j1fff = j1ff + mfel
c$$$      j2fff = j2ff + mfel

      ksgn = 2*j2 - 3
      iopw = 0
      IF ( lreshel .eq. 1 ) iopw = 1
      iops = 0
      IF ( lreshel .eq. 1 ) iops = 1

      DO izobsr = 1, nzobsr
         isco = izobsr
         DO izsrce = 1, nzsrce
            iscs = izsrce
            CALL kernel0(xwal,zwal,xwal,zwal,grdgre,grwp,j1,j2,ksgn,
     $           iopw,iops,lfele,nths2,nths)
         END DO
      END DO
c$$$      call matwrts ( grdgre,nths2,nths2,mth1,mth1,tmth,tmth,8,8,
c$$$     $        mth1,mth1,"grdgre at 2,2", outmod, iotty )
c$$$      call matwrts ( grdgre,nths2,nths2,j1f1,j2f1,j1fff,j2fff,8,8,
c$$$     $        j1f1,j2f1,"grdgre at 2,2,ff", outmod, iotty )
c$$$      call matwrtn ( grri,nths2,nfm2,1,1,nouk,lmax2,16,8,
c$$$     $        "GRRI AT 2 2", outmod, iotty )

c$$$      IF ( lhighn /= 2) GO TO 51
c$$$c.. Now put the K kernel in the correct position. Store CS submatrix
c$$$c   in AKKER
c$$$
c$$$      akker(j1f1:j1ff,j2f1:j2ff) = grdgre(j1f1:j1ff,j2f1:j2ff) +
c$$$     $     grdgre(j1ff1:j1fff,j2ff1:j2fff)
c$$$      akker(j1ff1:j1fff,j2ff1:j2fff) =  akker(j1f1:j1ff,j2f1:j2ff)
c$$$      akker(j1f1:j1ff,j2ff1:j2fff) =
c$$$     $     grdgre(j1ff1:j1fff,j2f1:j2ff) - grdgre(j1f1:j1ff,j2ff1:j2fff)
c$$$      akker(j1ff1:j1fff,j2f1:j2ff) = - akker(j1f1:j1ff,j2ff1:j2fff)
c$$$
c$$$ 51   CONTINUE

c$$$      call kernel0(xwal,zwal,xwal,zwal,grdgre,grwp,j1,j2,ksgn,
c$$$     $     iopw,iops,lfele,nths2,nths)
c$$$      call matwrtn ( grwp,nths,nths,1,1,mth,mth,16,8,
c$$$     $        "grwp at w-w", outmod, iotty )

      IF ( lreshel .eq. 1 ) THEN

c Store G(mth1) = G(1)
         grwp(1:mth2,mth1:mth2) = grwp(1:mth2,1:2)

         IF ( lfour .eq. 1 ) THEN
            gwin(mth+1:2*mth,1:jwal1) =
     $           matmul( grwp(1:mth,1:mth),rwvce(1:mth,1:jwal1) )
            gwin(mth+1:2*mth,jwal1+1:jwal2) =
     $           matmul( grwp(1:mth,1:mth),rwvco(1:mth,1:jwal1) )
c$$$            call fouran ( grwp, gwin, rwvce, mth,0 )
c$$$            call fouran ( grwp, gwin, rwvco, mth,jwal1 )
         END IF

      END IF
c
c Get Fourier representation of K(2,2). Use grwp as storage:
c  Note: Check nths is larger than 2*nfm for grwp
c
      IF ( checkd .and. (lfele .eq. 0) ) THEN
c
c......fourier analyse Kww matrix....
c
         call foura2 ( grdgre,mth,mth, grwp, 0 )
         call fanal1d ( grwp,nths,nths,0, wrkvr,wrkvi,nfm,nfm )
c
         call matwrtn(wrkvr,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Real Kww(l,l)", outmod,0 )
         call matwrtn(wrkvi,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Imag Kww(l,l)", outmod,0 )
c
      END IF
c
c.. ------------------------------------------------------------
c                           Wall-Plasma
c.. ------------------------------------------------------------
      j1 = 2
      j2 = 1

c  The following indices will shift the elements of the kernels the
c  correct amount for storage.

      j1f = 2*(j1-1) * mfel
      j2f = 2*(j2-1) * mfel
      j1f1 = j1f + 1
      j2f1 = j2f + 1
      j1ff = j1f + mfel
      j2ff = j2f + mfel
      j1ff1 = j1ff + 1
      j2ff1 = j2ff + 1
      j1fff = j1ff + mfel
      j2fff = j2ff + mfel

      ksgn = 2*j2 - 3
      iopw = 1
      iops = 0

      DO izobsr = 1, nzobsr
         isco = izobsr
         DO izsrce = 1, nzsrce
            iscs = izsrce
            CALL kernel0(xwal,zwal,xpla,zpla,grdgre,grwp,j1,j2,ksgn,
     $           iopw,iops,lfele,nths2,nths)
            IF ( lhighn == 2) THEN
               DO i = 1, mfel
                  DO j = 1, mfel
                     indx1 = ( 2*(j1-1)+(isco-1) ) * mfel + i
                     indx2 = ( 2*(j2-1)+(iscs-1) ) * mfel + j
                     grrin(indx1,indx2) = grwp(i,j)
                  END DO
               END DO
            END IF
         END DO
      END DO
c$$$      call matwrts ( grdgre,nths2,nths2,mth1,1,tmth,mth,8,8,
c$$$     $        mth1,1,"grdgre at 2,1", outmod, iotty )
c$$$      call matwrts ( grdgre,nths2,nths2,j1f1,j2f1,j1fff,j2fff,8,8,
c$$$     $        j1f1,j2f1,"grdgre at 2,1,ff", outmod, iotty )

      IF ( lhighn /= 2) GO TO 61
c$$$c.. Now put the K kernel in the correct position. Store CS submatrix
c$$$c   in AKKER
c$$$
c$$$      akker(j1f1:j1ff,j2f1:j2ff) = grdgre(j1f1:j1ff,j2f1:j2ff) +
c$$$     $     grdgre(j1ff1:j1fff,j2ff1:j2fff)
c$$$      akker(j1ff1:j1fff,j2ff1:j2fff) =  akker(j1f1:j1ff,j2f1:j2ff)
c$$$      akker(j1f1:j1ff,j2ff1:j2fff) =
c$$$     $     grdgre(j1ff1:j1fff,j2f1:j2ff) - grdgre(j1f1:j1ff,j2ff1:j2fff)
c$$$      akker(j1ff1:j1fff,j2f1:j2ff) = - akker(j1f1:j1ff,j2ff1:j2fff)

c. Now do the same for the G kernel. Use GRRI directly

c$$$      grri(j1f1:j1ff,j2f1:j2ff) = grrin(j1f1:j1ff,j2f1:j2ff) +
c$$$     $     grrin(j1ff1:j1fff,j2ff1:j2fff)
      grri(j1f1:j1ff,j2f1:j2ff) =     grrin(j1f1:j1ff,j2f1:j2ff)
      grri(j1ff1:j1fff,j2ff1:j2fff) = grrin(j1f1:j1ff,j2f1:j2ff)
c$$$      grri(j1f1:j1ff,j2ff1:j2fff) =
c$$$     $     grrin(j1ff1:j1fff,j2f1:j2ff) - grrin(j1f1:j1ff,j2ff1:j2fff)
      grri(j1f1:j1ff,j2ff1:j2fff) =   grrin(j1f1:j1ff,j2ff1:j2fff)
      grri(j1ff1:j1fff,j2f1:j2ff) = - grrin(j1f1:j1ff,j2ff1:j2fff)

      grri(j1f1:j1fff,j2f1:j2fff) = znorm * grri(j1f1:j1fff,j2f1:j2fff)

      call matwrtn ( grri,nths2,nfm2,1,1,nouk,lmax2,16,8,
     $        "GRRI AT 2 1", outmod, iotty )

 61   CONTINUE

c$$$      call kernel0(xwal,zwal,xpla,zpla,grdgre,grwp,j1,j2,ksgn,
c$$$     $     iopw,iops,lfele,nths2,nths)
c$$$      call matwrtn ( grwp,nths,nths,1,1,mth,mth,16,8,
c$$$     $        "grwp at w-p", outmod, iotty )
c      call matwrtn(grwp,nths,nths,9,9,"grwp", outmod,0 )

c Store G(mth1) = G(1), etc.

      grwp(1:mth2,mth1:mth2) = grwp(1:mth2,1:2)

c.....fourier analyse source points. real and imag. parts.
c     pack into one-d array. ( ... not for leqt1f or leqt2f... )

      if ( lfele .eq. 1 ) then
         call felang ( grwp, grri, cnqd,  mth,0 )
         call felang ( grwp, grri, snqd,  mth,jmax1 )
      end if
c            ! ! already normalized for lhighn=2
      IF ( (lfele == 3) .AND. ( lhighn /= 2 ) ) THEN
      CALL atpoint ( "felang3b","j1", j1,"zinf(1)", zinf(1),
     $     IOTTY, OUTMOD )
         CALL felang3 ( grwp,nths,nths, grri, mth,0 )
         CALL felang3 ( gwpi,mth2,mth2, grri, mth,jmax1 )
      END IF
      if ( lfour .eq. 1 ) then
         call fouran ( grwp, grri, cslth, mth,0 )
         call fouran ( grwp, grri, snlth, mth,jmax1 )
      end if
c
      if ( checkd .and. (lfele .eq. 0) ) then
c
c......fourier analyse gwp matrix....
c
         call fanal1d ( grri,nths2,nfm2,mth, wrkvr,wrkvi,nfm,nfm )
c
         call matwrtn(wrkvr,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Real GWP_ll", outmod,0 )
         call matwrtn(wrkvi,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Imag GWP_ll", outmod,0 )
c
      end if

 34   CONTINUE           ! if Farwal switch ends here.

      DEALLOCATE ( gwpi )
c
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      if ( check1 )
     .     call timer ( outmod, iotty, "aft kern and fourier" )

      if ( n .gt. 0.1 ) go to 46
c
c.....examine possibility of sum rules......
c
c
      do 45 j = 1, mth12
c
         do 44 ii = 1, 2
c
            summ(ii) = zero
            i1 = (ii-1) * mth + 1
            i2 = ii * mth
c
            do 40 i = i1, i2
c
c     is = 2*(i-1)*mth + j
c     summ(ii) = summ(ii) + grdgre(is)
               summ(ii) = summ(ii) + grdgre ( j, i )
c
 40         continue
c
 44      continue
c
         sum11 = summ(1)
         sum12 = summ(2)
         if ( check2 )
     $        write ( outmod, 9002 ) j, sum11, sum12
 9002    format (  1x, i4, 3x, 1p5e14.5 )
c
 45   continue
c
 46   continue
c
c.....repack grdgre if wall is far away.
c
c     if ( .not. farwal ) go to 49
c
c     do 48 j = 1, mth-1
c     do 48 i = 1, mth
c
c     ii = mth*j + i
c     jj = 2*mth*j + i
c     grdgre(ii) = grdgre(jj)
c
c     48 continue
c     49 continue
      if ( check2 ) then
c$$$         write(outmod,8010) (grdgre(i,1),i=1,9)
c$$$ 8010 format ( 1x, "grdgre = ",/, 10e11.4 )
c$$$          call matwrtn ( grdgre,nths2,nths2,1,1,mth12,mth12,16,8,
c$$$     $        "grdgre", outmod, 0 )
      end if
c
      if ( checkd ) then
      call matwrtn ( grwp,nths,nths,1,1,mth,mth,16,8,
     $        "grwp at end", outmod, iotty )
      end if

      call matwrtn ( grdgre,nths2,nths2,1,1,mth12,mth12,16,8,
     $        "grdgre at end", outmod, iotty )

      IF ( .NOT. farwal )
     $     CALL matwrts ( grdgre,nths2,nths2,1,1,tmth,tmth,16,8,
     $        mth1,mth1,"grdgre at 2,2", outmod, iotty )

      write ( iotty,  '(/,
     $     "Sum of first COLUMN in each block of GRDGRE:",/)')
      write ( outmod, '(/,
     $     "Sum of first COLUMN in each block of GRDGRE:",/)')
      do j1 = 1,2
         do j2 = 1,2
            sumg = 0.0
            do i = 1, mth
               sumg = sumg + grdgre( (j1-1)*mth+i,(j2-1)*mth+1 )
            end do
            write( iotty,  '("j1,j2, sumg= ",2i3,1pe11.3)' ) j1,j2,sumg
            write( outmod, '("j1,j2, sumg= ",2i3,1pe11.3)' ) j1,j2,sumg
         end do
      end do
c
c
      write ( iotty,  '(/,
     $     "Sum of first ROW in each block of GRDGRE:",/)')
      write ( outmod, '(/,
     $     "Sum of first ROW in each block of GRDGRE:",/)')
      do j1 = 1,2
         do j2 = 1,2
            sumg = 0.0
            do i = 1, mth
               sumg = sumg + grdgre( (j1-1)*mth+1,(j2-1)*mth+i )
            end do
            write( iotty,  '("j1,j2, sumg= ",2i3,1pe11.3)' ) j1,j2,sumg
            write( outmod, '("j1,j2, sumg= ",2i3,1pe11.3)' ) j1,j2,sumg
         end do
      end do
c
c.....gaussian elimination for inside shell.

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      if ( check1 )
     .     call timer ( outmod, iotty, "before leqt1f" )

c     if ( abs(n) .gt. 1.0e-5 ) go to 190
      if ( (abs(n) .gt. 1.e-5) .or. farwal .or. (ishape .gt. 100) )
     .     go to 190
c
c.......add arbitrary constant to grdgre to make matrix nonsingular.
c
c      cn0 = 1.0
c
      write (outmod,8180) cn0
 8180 format ( /, 5x, "constant,cn0, added to grdgre = ", f7.3 )
c
      grdgre(1:mth12,1:mth12) = grdgre(1:mth12,1:mth12) + cn0
c
 190  continue
c
c.... Store grri in grrin for Gaussian Elimination..
c$$$
c$$$      malgr = 2*mth+10
c$$$      jalgr = lmax2 + jwal2 + 10
c$$$
c$$$      ALLOCATE ( grrin(malgr,jalgr) )

      call atpoint ( "nouk=","mth12", mth12,"zinf(1)", zinf(1),
     $     iotty, outmod )

c... Number of right-hand sides, nrhs:

      nrhs = lmax2
c 88888888888888888888888888888888888888888888888888888
c ...Check next statement. Should be lmax2 instead of jmax1???
c8888888888888888888888888888888888888888888888888888
      IF ( farwal .and. lreshel .eq. 0 ) nrhs = lmax2
      IF ( ldqdtw == 2 ) nrhs = mth

      grrin(1:nouk,1:nrhs) = grri(1:nouk,1:nrhs)

c.. If lhighn == 2, regroup GRDGRE:
c.. And subtract unity because of doubling the residues when
c    comnbining CosCos and SinSin contributions

      call atpoint ( "after nouk=","nouk", nouk,"zinf(1)", zinf(1),
     $     iotty, outmod )

      IF ( lhighn == 2 ) THEN
         CALL regrupk ( grdgre, nths2,nths2, mfel )
      END IF

      CALL matwrtn ( grdgre,nths2,nths2,1,1,mth12,mth12,16,8,
     $        "GRDGRE AFTER REGRUPK", outmod, iotty )
      CALL matwrtn ( grrin,malgr,jalgr,1,1,nouk,lmax2,16,8,
     $        "GRRIN AT END", outmod, iotty )

c      DEALLOCATE ( akker )
      CALL atpoint ( "deall akker","mth12", mth12,"zinf(1)", zinf(1),
     $     iotty, outmod )

c...... If resistive shell, append shell sources:
c  .... Also, if internal coils append coil sources.
c.. Update nhrs also
c
      IF ( lreshel .eq. 1 ) THEN
         nrhs = lmax2 + jwal2
         grrin(1:mth12,lmax2+1:nrhs) = gwin(1:mth12,1:jwal2)
      END IF

      IF ( lfbcoila .eq. 1 ) THEN
         nrhs = nrhs + 1
         grrin(1:mth,nrhs) = skcoilpa(1:mth)
         IF ( .NOT. farwal ) grrin(mth+1:2*mth,nrhs) = skcoilwa(1:mth)
      END IF

      IF ( lfbcoilb .eq. 1 ) THEN
         nrhs = nrhs + 1
         grrin(1:mth,nrhs) = skcoilpb(1:mth)
         IF ( .NOT. farwal ) grrin(mth+1:2*mth,nrhs) = skcoilwb(1:mth)
      END IF
c
c     call gelg ( grri, grdgre, mth12, lmax2, epsq, ier )
c     call leqt1f(grdgre,lmax2,mth12,nths2,grri,idgt,workl,ier)

      ier = 0

c$$$      call gelima ( grdgre,nths2,grrin,nths2,mth12,lmax2,
c$$$     $     grri,nths2,work1,workg1,nths2,workg2,nths2,ier )
c$$$      write ( outmod,8050 ) ier
c$$$      write ( iotty, 8050 ) ier
c$$$ 8050 format (/, 1x, " ier in f04aee = ", i5 )

c**********************************************************************
c*********************************************************************
c
c  IV-138(6): Solution in region I for \chi
c
c**********************************************************************
c***********************************************************************
      CALL atpoint ( "before glimb","nouk", nouk,"qa1", qa1,
     $     iotty, outmod )

c      CALL gelimb ( grdgre,nths2,grrin,malgr,nouk,nrhs,
c     $     grri,nths2,work1,ier )
      CALL f04aaf( grdgre,nths2,grrin,malgr,nouk,nrhs,
     $     grri,nths2,work1,ier )

      WRITE ( outmod,8050 ) ier
      WRITE ( iotty, 8050 ) ier
 8050 FORMAT (/, 1x, " ier in f04aae = ", i5 )

c-----------------------------------------------------------------------
c.... Now: grri(:,lmax2+jwal2+1) contains the pa and wa internal coils.
c          grri(:,lmax2+jwal2+2) contains the pb and wb internal coils.
c-----------------------------------------------------------------------

      IF ( .NOT. checkd ) GO TO 57   !!! About line 1161

c.. Write out the response matrix in Region 1 to disk for diagnostics.
      
      ioresp_1 = 176
      OPEN (  ioresp_1, FILE='RESPONSE_REG_1', STATUS='REPLACE',
     $     FORM='FORMATTED' )

      WRITE ( ioresp_1, '("%", 5x, "Response Matrix in Region 1" )' )

      WRITE ( ioresp_1, '("%",1x, "Toroidal mode number = ", /,
     $     f6.2 )' ) n

      WRITE ( ioresp_1, '("%",1x,
     $     "Number of plasma points on the open domain, mth = ",
     $     /, i5 )' )
     $     mth
      WRITE ( ioresp_1, '("%",1x,
     $     "Number of wall points on the open domain, mthw = ",
     $     /, i5 )' )
     $     mth
      WRITE ( ioresp_1, '("%",1x, "Plasma modes: lmin = ", /,
     $     i5 )' ) lmin(1)
      WRITE ( ioresp_1, '("%",1x, "Plasma modes: lmax = ", /,
     $     i5 )' ) lmax(1)

      WRITE ( ioresp_1, '("%",1x, "Wall modes: evens  = ", /,
     $     i5 )' ) jwal1
      WRITE ( ioresp_1, '("%",1x, "Wall modes: odds = ", /,
     $     i5 )' ) jwal1

      nnzc = nrhs - (jmax2+jwal2)
      WRITE ( ioresp_1, '("%",1x,
     $     "Number of coils in region 1, nc1:", /, i4 )' ) nnzc

      WRITE ( ioresp_1, '("%",1x, "Number of RHS = ", /,
     $     I5 )' ) nrhs

      WRITE ( ioresp_1, '("%",1x,
     $     "Plasma parameterization, X(mth):", i5 )' ) mth
      WRITE ( ioresp_1, '( (1x, 8es14.6) )' ) (xinf(i), i=1, mth)

      WRITE ( ioresp_1, '("%",1x,
     $     "Plasma parameterization, Z(mth):", i5 )' ) mth
      WRITE ( ioresp_1, '( (1x, 8es14.6) )' ) (zinf(i), i=1, mth)

      WRITE ( ioresp_1, '("%",1x,
     $     "Wall parameterization, X(mthw):", i5 )' ) mth
      WRITE ( ioresp_1, '( (1x, 8es14.6) )' ) (xwal(i), i=1, mth)

      WRITE ( ioresp_1, '("%",1x,
     $     "Wall parameterization, Z(mthw):", i5 )' ) mth
      WRITE ( ioresp_1, '( (1x, 8es14.6) )' ) (zwal(i), i=1, mth)

      WRITE ( ioresp_1, '("%",1x,
     $     "Response Matrix(theta_p + theta_w, modes_p + modes_w + nc)"
     $     " has nc coils responses appended.",/, "%", 1x,
     $     "Inner do loop is over first index:",/, "%", 1x,
     $     "resp_1[mth+mthw,2*(lmax-lmin+1+jwal1)+nc1]: ",
     $     "(2*",i4,",",i4,"+",i4,"+",i4")" )' ) mth, jmax2, jwal2, nnzc

      DO izl = 1, nrhs
         WRITE ( ioresp_1, '( (1x, 8es14.6) )' )
     $        (grri(i,izl), i=1, mth12)
      END DO

      CLOSE ( UNIT = ioresp_1 )
 
 57   CONTINUE    !!! IF .NOT. checkd 

c Store the Chi on the surfaces of the plasma and the wall due to
c     the surface-coil in skcoil arrays:

      indxrhs = lmax2
      IF ( lreshel == 1 ) indxrhs = lmax2 + jwal2

      IF (lfbcoila == 1 ) THEN
         skcoilpa(1:mth) = grri(1:mth,indxrhs+1)
         CALL vecwrt ( mth, skcoilpa, "skcoilpa", 1, mth, outmod,iotty )
      END IF

      IF (lfbcoilb == 1 ) skcoilpb(1:mth) = grri(1:mth,indxrhs+2)
      IF ( .NOT. farwal ) THEN
         IF ( lfbcoila == 1)
     $        skcoilwa(1:mth) = grri(mth+1:2*mth,indxrhs+1)
         IF ( lfbcoilb ==1 )
     $        skcoilwb(1:mth) = grri(mth+1:2*mth,indxrhs+2)
      END IF

      CALL matwrtn ( grri,nths2,nfm2,1,ln,mth,nrhs,16,8,
     $        "Chi-response matrix", outmod, 0 )
      CALL matwrts ( grri,nths2,nfm2,1,ln, 25,9, 25,9,
     $     1,1,"Chi-response matrix", outmod,0 )

      ALLOCATE ( zwkpl(mth2) )

      DO izpl = 1, 4
         izplz = (izpl-1)*(nrhs/4) + 1
         zwkpl(1:mth) = grri(1:mth,izplz)
         zwkpl(mth1) = zwkpl(1)
         CALL pospl1 ( pigrd,zwkpl, 1,mth1, izpl,"grri_resp", izplz )
      END DO

      CALL framep(jobid, ff)

      DEALLOCATE ( zwkpl )

      IF ( lm3dc1 == 0 ) GO TO 685   !!! ABOUT LINE 1426

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c... Construct Response Matrix for M3D_C1 from GRRI
c   Use complex algebra
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      IF ( (lm3dc1 == 1) .AND. (ipshp == 0) ) mthout = mthin
      mthout1 = mthout+1
      mthout2 = mthout+2

      mthmx = MAX(mth,mthout)
      mtxdim = mthmx + 5

      WRITE ( OUTMOD, '(/, 5x, "MTHIN, MTH, MTHOUT, MTHMX = ", 4I5 )' )
     $     mthin, mth, mthout, mthmx
      WRITE ( IOTTY,  '(/, 5x, "MTHIN, MTH, MTHOUT, MTHMX = ", 4I5 )' )
     $     mthin, mth, mthout, mthmx

      ALLOCATE ( zgrij0(mtxdim,mtxdim), zgrij(mtxdim,mtxdim),
     $     zgrijb(mtxdim,mtxdim) )
      ALLOCATE ( zdlent(mtxdim) )
      ALLOCATE ( zgrchi(mtxdim,mtxdim), zgrbth(mtxdim,mtxdim) )
      ALLOCATE ( zgri(mtxdim), zgro(mtxdim), zgrop(mtxdim) )
      ALLOCATE ( zijbp(mtxdim,mtxdim), zgrbph(mtxdim,mtxdim) ) 
      ALLOCATE ( expmilt(mtxdim,jmax1+5), cgrri(mtxdim,jmax1+5) )

c... Construct exp(-il\theta):
      expmilt(1:mth,1:jmax1) =
     $     CMPLX( coslt(1:mth,1:jmax1),-sinlt(1:mth,1:jmax1) )

c... Make the chi response complex:
      cgrri(1:mth,1:jmax1) =
     $     CMPLX( grri(1:mth,1:jmax1),grri(1:mth,jmax1+1:jmax2) )

c... Now transform to real space. Sum over exp(-il\theta):
      zgrij0(1:mth,1:mth) = MATMUL ( cgrri(1:mth,1:jmax1),
     $     TRANSPOSE ( expmilt(1:mth,1:jmax1) ) )

c.. Complete for closed domain:
      DO i = 1, mth
         zgrij0(i,mth1) = zgrij0(i,1)
      END DO

      DO j = 1, mth1
         zgrij0(mth1,j) = zgrij0(1,j)
      END DO

c$$$      CALL matwrtn(zgrij0,mtxdim,mtxdim,1,1,mth1,mth1,8,8,
c$$$     $     "ZGRIJ(i,j)", outmod,0 )
c$$$      CALL vacasyma ( zgrij0,mtxdim,mth1, 1.0, "REAL zgrij",
c$$$     $     outmod,iotty )

c... Multiply source by J*grad Z for response to delta B rather than CalB

c... Ensure periodicity:
      xplap(mth1) = xplap(1)
      zplap(mth1) = zplap(1)
      xplap(mth+2) = xplap(2)
      zplap(mth+2) = zplap(2)

      zdlent(1:mth2) = SQRT ( xplap(1:mth2)**2 + zplap(1:mth2)**2 )

      DO i = 1, mth1
         zgrij(i,1:mth1) =
     $        zgrij0(i,1:mth1) * xpla(1:mth1) * zdlent(1:mth1)
      END DO

c$$$      CALL matwrtn(zgrij,mtxdim,mtxdim,1,1,mth1,mth1,8,8,
c$$$     $     "ZGRIJ(i,j)*zdlent", outmod,0 )
c$$$      CALL vacasyma ( zgrij,mtxdim,mth1, 1.0, "REAL zgrij*zdlent",
c$$$     $     outmod,iotty )

c... To check, transform back to Fourier Space:

      jzzdim = jmax1 + 5
c$$$
c$$$      ALLOCATE ( zgric(mtxdim,jzzdim), zgris(mtxdim,jzzdim),
c$$$     $     zgrcc(jzzdim,jzzdim), zgrss(jzzdim,jzzdim),
c$$$     $     zgrll(jzzdim,jzzdim) )
c$$$
c$$$      zgric(1:mth2,1:jmax1) = MATMUL ( zgrij0(1:mth2,1:mth),
c$$$     $     coslt(1:mth,1:jmax1) )
c$$$      zgrcc(1:jmax1,1:jmax1) =
c$$$     $     MATMUL ( TRANSPOSE ( coslt(1:mth,1:jmax1) ),
c$$$     $     zgric(1:mth,1:jmax1) )
c$$$
c$$$      zgris(1:mth2,1:jmax1) = MATMUL ( zgrij0(1:mth2,1:mth),
c$$$     $     sinlt(1:mth,1:jmax1) )
c$$$      zgrss(1:jmax1,1:jmax1) =
c$$$     $     MATMUL ( TRANSPOSE ( sinlt(1:mth,1:jmax1) ),
c$$$     $     zgris(1:mth,1:jmax1) )
c$$$
c$$$      zgrll = ( zgrcc + zgrss ) / (mth*mth)
c$$$
c$$$      CALL matwrtn(zgrll,jzzdim,jzzdim,ln,ln,jmax1,jmax1,jdel,jdel,
c$$$     $     "ZGRLL", outmod, 0 )
c$$$      CALL vacasyma ( zgrll,jzzdim,jmax1, 1.0, "REAL zgrll",
c$$$     $     outmod,iotty )
c$$$
c$$$      DEALLOCATE ( zgric, zgris, zgrcc, zgrss, zgrll )


c.. Now calculate B_theta response: Differentiate wrt to observer index,
c    i.e., column.
c.. Put it in ZGRIJB. Use ZGRO and ZGROP as temporary storage. 
c    Do for each source point.
c.. Need a grid in dl/dtheta first:

      DO i = 1, mth2
c.. Real Part:
         zgro(1:mth2) = REAL( zgrij(1:mth2, i) )
         CALL difspl ( mth, pigrd, zgro, zgrop )
         zgrop(mth1) = zgrop(1)
         zgrop(mth2) = zgrop(2)
         zgrop(1:mth2) = zgrop(1:mth2) / zdlent(1:mth2)
c... Save zgrop:
         zgrijb(1:mth2, i) = zgrop(1:mth2)
c.. Imaginary Part:
         zgro(1:mth2) = AIMAG( zgrij(1:mth2, i) )
         CALL difspl ( mth, pigrd, zgro, zgrop )
c... Store Real zgrop back into zgro:
         zgro(1:mth2) = zgrijb(1:mth2,i)
         zgrop(mth1) = zgrop(1)
         zgrop(mth2) = zgrop(2)
         zgrop(1:mth2) = zgrop(1:mth2) / zdlent(1:mth2)
c.. Now construct Complex zgrijb:
         zgrijb(1:mth2,i) = CMPLX( zgro(1:mth2),zgrop(1:mth2) )
c         gpsjp(i) = xpla(i) * dlenth(i)
      END DO

c... Since B_theta points the other way, this needs to change sign:
      zgrijb = - zgrijb
     

c$$$      CALL matwrtn(zgrijb,mtxdim,mtxdim,1,1,mth,mth,8,8,
c$$$     $     "Bth(i,j), mth pts", outmod,0 )

c....   Calculate B_phi  response. This is -i*n*Chi/X_observer.

      ziii = CMPLX( 0.0,1.0 )
      DO j = 1, mth2
         zijbp(1:mth2,j) = - ziii*n*zgrij(1:mth2,j) / xpla(1:mth2)
      END DO

c$$$
c$$$      CALL matwrtn(zijbp,mtxdim,mtxdim,1,1,mth,mth,8,8,
c$$$     $     "Bph(i,j), mth pts", outmod,0 )

c... Now interpolate chi to mthout points

      CALL transc( zgrij, zgrchi, mtxdim,mtxdim, mth,mthout )

c$$$      CALL matwrtn(zgrchi,mtxdim,mtxdim,1,1,mthout1,mthout1,8,8,
c$$$     $     "ZGRCHI(i,j), mthout pts", outmod,0 )
c$$$      CALL vacasyma ( zgrchi,mtxdim,mthout1, 1.0, "REAL zgrchi",
c$$$     $     outmod,iotty )

c.... Now Bth

      CALL transc( zgrijb, zgrbth, mtxdim,mtxdim, mth,mthout )

c$$$      CALL matwrtn(zgrbth,mtxdim,mtxdim,1,1,mthout1,mthout1,8,8,
c$$$     $     "ZGRBTH(i,j), mthout pts", outmod,0 )

c.... Now Bph

      CALL transc( zijbp, zgrbph, mtxdim,mtxdim, mth,mthout )

c$$$      CALL matwrtn(zgrbph,mtxdim,mtxdim,1,1,mthout1,mthout1,8,8,
c$$$     $     "ZGRBPH(i,j)/in, mthout pts", outmod,0 )

c$$$      CALL matwrtn(zgrijb,mtxdim,mtxdim,1,1,mth1,mth1,8,8,
c$$$     $     "IMAG ZGRIJ(i,j)", outmod,0 )
c$$$      CALL vacasyma ( zgrijb,mtxdim,mth1, 1.0, "IMAG zgriji",
c$$$     $     outmod,iotty )

c.. Revert back to Counterclockwise.

      CALL  reversematc (  zgrchi, mtxdim,mtxdim, mthout )
      CALL  reversematc (  zgrbth, mtxdim,mtxdim, mthout )
      CALL  reversematc (  zgrbph, mtxdim,mtxdim, mthout )
      
c.. Write out quantities for M3D-C1
      
      iomc1 = 376
      OPEN (  iomc1, FILE='RESPONSE-M3DC1', STATUS='REPLACE',
     $     FORM='FORMATTED' )

      WRITE ( iomc1,
     $     '(5x, "Response Matrix in theta space for M3D_C1 on:",
     $     2x, A )' ) date_array

      WRITE ( iomc1, '(/,1x, "No of internal points used = ", i5 )' )
     $     mth1
      WRITE ( iomc1, '(/,1x, "Min Fourier mode used = ", i5 )' ) lmin(1)
      WRITE ( iomc1, '(  1x, "Max Fourier mode used = ", i5 )' ) lmax(1)

      nzzz = n+0.1
      WRITE ( iomc1, '(/,1x, "Toroidal mode no. = ", i5 )' ) nzzz

      WRITE ( iomc1, '(/,1x,
     $     "Number of surface points on the closed domain = ", i5 )' )
     $     mthout1

      WRITE ( iomc1,
     $     '(/,1x, "B_theta response Matrix, Rth(obs,srce):" )' )

      DO izi = 1, mthout1
         WRITE ( iomc1, '(/, 1x, "i_obs = ", i5 )' ) izi
         WRITE ( iomc1, '( (1x, 8es14.6) )' )
     $        (zgrbth(izi,j), j=1, mthout1)
      END DO

      WRITE ( iomc1,
     $     '(/,1x, "B_phi response Matrix, Rph(obs,src):" )' )

      DO izi = 1, mthout1
         WRITE ( iomc1, '(/, 1x, "i_obs = ", i5 )' ) izi
         WRITE ( iomc1, '( (1x, 8es14.6) )' )
     $        (zgrbph(izi,j), j=1, mthout1)
      END DO

      CLOSE ( UNIT = 376 )

      DEALLOCATE ( zgri, zgro, zgrop, cgrri, expmilt )
      DEALLOCATE ( zgrij, zgrij0, zgrijb, zgrchi, zgrbth, zdlent )
      DEALLOCATE ( zijbp, zgrbph ) 

 685  CONTINUE          !  IF ( lm3dc1 == 0 )

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      IF ( lrnim3d == 0 ) GO TO 695

c.. Write out quantities for Nimrod
      
      ionim = 276
      OPEN (  ionim, FILE='RESPONSE-NIM', STATUS='REPLACE',
     $     FORM='FORMATTED' )

      WRITE ( ionim, '(5x, "Response Matrix in External Region" )' )

      WRITE ( ionim, '(/,1x,
     $     "Number of surface points on the open domain = ", i5 )' )
     $     mth
      WRITE ( ionim, '(/,1x, "Fourier index: lmin = ", i5 )' ) lmin(1)
      WRITE ( ionim, '(/,1x, "Fourier index: lmax = ", i5 )' ) lmax(1)

      WRITE ( ionim, '(/,1x, "Surface parameterization, X:" )' )
      WRITE ( ionim, '( (1x, 8es14.6) )' ) (xinf(i), i=1, mth)

      WRITE ( ionim, '(/,1x, "Surface parameterization, Z:" )' )
      WRITE ( ionim, '( (1x, 8es14.6) )' ) (zinf(i), i=1, mth)

      WRITE ( ionim, '(/,1x, "Cosine response Matrix, CR(i,l):" )' )
      DO izl = 1, jmax1
         lzz = lmin(1) - 1 +  izl
         WRITE ( ionim, '(/, 1x, "l = ", i5 )' )  lzz
         WRITE ( ionim, '( (1x, 8es14.6) )' ) (grri(i,izl), i=1, mth)
      END DO

      WRITE ( ionim, '(/,1x, "Sine response Matrix, SR(i,l):" )' )
      DO izl = jmax1 + 1, jmax2
         lzz = lmin(1) - 1 +  izl - jmax1
         WRITE ( ionim, '(/, 1x, "l = ", i5 )' )  lzz
         WRITE ( ionim, '( (1x, 8es14.6) )' ) (grri(i,izl), i=1, mth)
      END DO

      CLOSE ( UNIT = 276 )

c..Wall response to dQ/dt for Nimrod or M3D
      lmnf = -10
      lmxf = 20

      IF ( ldqdtw == 1 ) THEN
         CALL dqdtwl ( grri,nths2,nfm2, mth,nrhs, n,
     $        pigrd, xpla, xplap,zplap, lmnf,lmxf, outmod )
      END IF
      IF ( ldqdtw == 2 .OR. ldqdtw == 3 ) THEN
         CALL dqdtw ( grri,nths2,nfm2, mth,nrhs, n,
     $        pigrd, xpla, xplap,zplap, ldqdtw, outmod )
      END IF
c$$$      IF ( ldqdtw == 3 ) THEN
c$$$         CALL dqdtw ( grri,nths2,nfm2, mth,nrhs, n,
c$$$     $        pigrd, xpla, xplap,zplap, ldqdtw, outmod )
c$$$      END IF

 695  CONTINUE

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      if ( check1 )
     .     call timer ( outmod, iotty, "after gaus-elim" )
c
      if ( lreshel .eq. 0 )  go to 33  !!! About line 1780
c
c.. Extract wall source contribution for both plasma and wall observer pts
c
      gwin(1:mth12, 1:jwal2) = grri(1:mth12, lmax2+1:lmax2+jwal2)

c... ....   Region II: Wall to Coil region. Treat wall as if its
c           the interior boundary:

c.. ------------------------------------------------------------
c                           Wall-Wall - plus
c.. ------------------------------------------------------------
      call atpoint ( "Wall-Wall + ","jwal1", jwal1,"xinf(1)", xinf(1),
     $     iotty, outmod )
      j1 = 1
      j2 = 1
      ksgn = 2*j2 - 3
      iopw = 1
      iops = 1
      CALL kernel0(xwal,zwal,xwal,zwal,grdgre,grwp,j1,j2,ksgn,
     $     iopw,iops,lfele,nths2,nths)
c$$$      call matwrtn ( grwp,nths,nths,1,1,mth,mth,16,8,
c$$$     $        "grwp at w-w + ", outmod, iotty )
c
c Store G(mth1) = G(1)
      grwp(1:mth2,mth1:mth2) = grwp(1:mth2,1:2)
c
      IF ( lfour .eq. 1 ) THEN
         gwot(1:mth,1:jwal1) =
     $        matmul( grwp(1:mth,1:mth),rwvce(1:mth,1:jwal1) )
         gwot(1:mth,jwal1+1:jwal2) =
     $        matmul( grwp(1:mth,1:mth),rwvco(1:mth,1:jwal1) )
      END IF
c
      IF ( lfbcoil .eq. 0 ) GO TO 8010
c
c..------------------------------------------------------------
c                          Wall-Coil
c..------------------------------------------------------------
      CALL atpoint ( "Wall-Coil","mtcoil", mtcoil,"xcw(1)", xcw(1),
     $     iotty, outmod )
c
c...       Use off diagonal block of GRDGRE as storage.
c            Don't need GREEN (a.k.a grwp) here.
c
      j1 = 1
      j2 = 2
      ksgn = 2*j2 - 3
      iopw = 0
      iops = 0
      CALL kernel0(xwal,zwal,xcw,zcw,grdgre,grwp,j1,j2,ksgn,
     $     iopw,iops,lfele,nths2,nths)
c$$$      call matwrtn ( grwp,nths,nths,1,1,mth,mth,16,8,
c$$$     $        "grwp at w-c ", outmod, iotty )
c
c..        Integrate over the coil:
c
      skcoil = 0.0
      IF ( lfbcoil .eq. 1 ) THEN
         skcoil(1:mth) = gainc * matmul (
     $        grdgre(1:mth,mth+1:mth+mtcoil),ccoil(1:mtcoil) )
      END IF
c
      CALL vecwrt ( mth, skcoil, "skcint", 1, mth, outmod, iotty )
c
 8010 CONTINUE
c
c     .... Region II Gaussian Elimination for Wall Source.
c
c.... Store gwot in grrin for RHS of  Gaussian Elimination..
c...   Append coil source if needed. Treat the coil as a  separate
c      source for now. Energize it with different sensors later.
c
      grrin(1:mth,1:jwal2) = gwot(1:mth,1:jwal2)
      jwal2c = jwal2
      IF ( lfbcoil .eq. 1 ) THEN
         jwal2c = jwal2 + 1
         grrin(1:mth,jwal2c) = skcoil(1:mth)
      END IF
c
         ier = 0
c         CALL gelimb ( grdgre,nths2,grrin,malgr,mth,jwal2c,
c     $        gwot,nths2,work1,ier )
         CALL f04aaf ( grdgre,nths2,grrin,malgr,mth,jwal2c,
     $        gwot,nths2,work1,ier )
         WRITE ( outmod, '(/, 1x,"ier in gauss-walin = ", i5)' ) ier
         WRITE ( iotty,  '(/, 1x,"ier in gauss-walin = ", i5)' ) ier
         WRITE ( OUTMOD, '(/, 1X,"COIL SOURCE IN GWOT(:,JWAL2+1)")' )
         WRITE ( iotty,  '(/, 1x,"Coil source in GWOT(:,jwal2+1)")' )
c
c..   Remember GWOT(:,jwal2+1) now contains the coil source:
c    ----------------------------------------------------

c   Write out the response of the plasma and shell to the coils

c   First, plasma to coila in file R-PCA, then plasma to coilb:

         IF ( lfbcoila == 1 ) THEN

            clrsp(1:mth) = skcoilpa(1:mth)

            OPEN ( UNIT=53, FILE='R-PCA' )

            WRITE ( 53, '("Response of the Plasma to Coil-a")' )
            WRITE ( 53, '("IN PLASMAS GRID, FIRST")' )
            WRITE ( 53, '("mth = ", I4,/)' ) mth
            WRITE ( 53, '(2x, "I", 7x,"XPLA",11x,"ZPLA",10x,"XPLAP",10x,
     $           "ZPLAP",11x,"R-PCA",/)' )
            WRITE ( 53, '(i4, 5E15.7)' )
     $           ( i, xpla(i),zpla(i), xplap(i),zplap(i),
     $           clrsp(i), i = 1, mth )

            CLOSE ( UNIT=53 )

c.. Plot the current response on the plasma due to coil-a.
c     Shift the points first

            wkxx(1:mth1) = (/ clrsp(mth/2+1:mth),clrsp(1:mth/2+1) /)
            wkyy(1:nths) = 0.0
            workl(1:mth1)  = (/ xpla(mth/2+1:mth),xpla(1:mth/2+1) /)
            workl2(1:mth1) = (/ zpla(mth/2+1:mth),zpla(1:mth/2+1) /)

            nbngr = ncgr
            iwrto = outmod

c.. Can plot the surface component of the magnetic field if 
c    K_phi -> Q_theta, and K_theta -> - Q_phi. 
c  ikb = 1 for this  option.

            call kcur0 ( workl,workl2,ellp,pgrd,mth1,wkxx,wkyy,n,
     $           nbngr,wkgrd,wxgrd,wzgrd,wkdxt,wkdzt,
     $           wkt2,wkp2,nths,1,"P_Resp._to_C-a",lfbcoila, iwrto )

         END IF

         IF ( lfbcoilb == 1 ) THEN

            clrsp(1:mth) = skcoilpb(1:mth)

            OPEN ( UNIT=53, FILE='R-PCB' )

            WRITE ( 53, '("Response of the Plasma to Coil-b")' )
            WRITE ( 53, '("IN PLASMAS GRID, FIRST")' )
            WRITE ( 53, '("mth = ", I4,/)' ) mth
            WRITE ( 53, '(2x, "I", 7x,"XPLA",11x,"ZPLA",10x,"XPLAP",10x,
     $           "ZPLAP",11x,"R-PCB",/)' )
            WRITE ( 53, '(i4, 5E15.7)' )
     $           ( i, xpla(i),zpla(i), xplap(i),zplap(i),
     $           clrsp(i), i = 1, mth )

            CLOSE ( UNIT=53 )

c.. Plot the current response on the plasma due to coil-b.
c     Shift the points first

            wkxx(1:mth1) = (/ clrsp(mth/2+1:mth),clrsp(1:mth/2+1) /)
            wkyy(1:nths) = 0.0
            workl(1:mth1)  = (/ xpla(mth/2+1:mth),xpla(1:mth/2+1) /)
            workl2(1:mth1) = (/ zpla(mth/2+1:mth),zpla(1:mth/2+1) /)

            nbngr = ncgr
            iwrto = outmod
            call kcur0 ( workl,workl2,ellp,pgrd,mth1,wkxx,wkyy,n,
     $           nbngr,wkgrd,wxgrd,wzgrd,wkdxt,wkdzt,
     $           wkt2,wkp2,nths,1,"P_Resp._to_C-b",lfbcoilb, iwrto )

         END IF

         IF ( farwal ) GO TO 651

c.   Now, wall to coila in file R-SCA, then wall to coilb:

         IF ( lfbcoila ==1 ) THEN

            clrsp(1:mth) = skcoilwa(1:mth)
            
            OPEN ( UNIT=53, FILE='R-SCA' )

            WRITE ( 53, '("Response of the Shell to Coil-a")' )
            WRITE ( 53, '("IN SHELLS GRID, FIRST")' )
            WRITE ( 53, '("mth = ", I4,/)' ) mth
            WRITE ( 53, '(2x, "I", 7x,"XWAL",11x,"ZWAL",10x,"XWALP",10x,
     $           "ZWALP",11x,"R-SCA",/)' )
            WRITE ( 53, '(i4, 5E15.7)' )
     $           ( i, xwal(i),zwal(i), xwalp(i),zwalp(i),
     $           clrsp(i), i = 1, mth )

            CLOSE ( UNIT=53 )
         
c.. Plot the current response on the shell due to coil-a.
c     Shift the points first

            wkxx(1:mth1) = (/ clrsp(mth/2+1:mth),clrsp(1:mth/2+1) /)
            wkyy(1:nths) = 0.0
            workl(1:mth1)  = (/ xwal(mth/2+1:mth),xwal(1:mth/2+1) /)
            workl2(1:mth1) = (/ zwal(mth/2+1:mth),zwal(1:mth/2+1) /)

            nbngr = ncgr
            iwrto = 0
            call kcur0 ( workl,workl2,ell,wgrd,mth1,wkxx,wkyy,n,
     $           nbngr,wkgrd,wxgrd,wzgrd,wkdxt,wkdzt,
     $           wkt2,wkp2,nths,0,"S_Resp._to_C-a",lfbcoila, iwrto )

         END IF

         IF ( lfbcoilb ==1 ) THEN

            clrsp(1:mth) = skcoilwb(1:mth)
            
            OPEN ( UNIT=53, FILE='R-SCB' )
            
            WRITE ( 53, '("Response of the Shell to Coil-b")' )
            WRITE ( 53, '("IN SHELLS GRID, FIRST")' )
            WRITE ( 53, '("mth = ", I4,/)' ) mth
            WRITE ( 53, '(2x, "I", 7x,"XWAL",11x,"ZWAL",10x,"XWALP",10x,
     $           "ZWALP",11x,"R-SCB",/)' )
            WRITE ( 53, '(i4, 5E15.7)' )
     $           ( i, xwal(i),zwal(i), xwalp(i),zwalp(i),
     $           clrsp(i), i = 1, mth )
            
            CLOSE ( UNIT=53 )

c.. Plot the current response on the shell due to coil-b.
c     Shift the points first

            wkxx(1:mth1) = (/ clrsp(mth/2+1:mth),clrsp(1:mth/2+1) /)
            wkyy(1:nths) = 0.0
            workl(1:mth1)  = (/ xwal(mth/2+1:mth),xwal(1:mth/2+1) /)
            workl2(1:mth1) = (/ zwal(mth/2+1:mth),zwal(1:mth/2+1) /)

            nbngr = ncgr
            iwrto = 0
            call kcur0 ( workl,workl2,ell,wgrd,mth1,wkxx,wkyy,n,
     $           nbngr,wkgrd,wxgrd,wzgrd,wkdxt,wkdzt,
     $           wkt2,wkp2,nths,0,"S_Resp._to_C-b",lfbcoilb, iwrto )

         END IF
         
 651     CONTINUE

c   Now the external coil:

         IF ( lfbcoil ==1 ) THEN

            clrsp(1:mth) = gwot(1:mth,jwal2+1)
            
            OPEN ( UNIT=53, FILE='R-SC' )

            WRITE ( 53, '("Response of the Shell to the Coil")' )
            WRITE ( 53, '("IN SHELLS GRID, FIRST")' )
            WRITE ( 53, '("mth = ", I4,/)' ) mth
            WRITE ( 53, '(2x, "I", 7x,"XWAL",11x,"ZWAL",10x,"XWALP",10x,
     $           "ZWALP",11x,"R-SC",/)' )
            WRITE ( 53, '(i4, 5E15.7)' )
     $           ( i, xwal(i),zwal(i), xwalp(i),zwalp(i),
     $           gwot(i,jwal2+1), i = 1, mth )

            CLOSE ( UNIT=53 )

c.. Plot the current response on the shell due to the coil.
c     Shift the points first

            wkxx(1:mth1) = (/ gwot(mth/2+1:mth,jwal2+1),
     $           gwot(1:mth/2+1,jwal2+1) /)
            wkyy(1:nths) = 0.0
            workl(1:mth1)  = (/ xwal(mth/2+1:mth),xwal(1:mth/2+1) /)
            workl2(1:mth1) = (/ zwal(mth/2+1:mth),zwal(1:mth/2+1) /)

            nbngr = ncgr
            iwrto = 0
            call kcur0 ( workl,workl2,ell,wgrd,mth1,wkxx,wkyy,n,
     $           nbngr,wkgrd,wxgrd,wzgrd,wkdxt,wkdzt,
     $           wkt2,wkp2,nths,0,"S_Resp._to_C",lfbcoil, iwrto )

         END IF

      call atpoint ( "bef-RESSHEL", "jwal2c",jwal2c,"a", a,
     $        iotty,outmod )
c
      ALLOCATE ( v_cpwwp(mfjm2,mfjm2) )

      CALL timer ( outmod, iotty, "Before resshel" )
      CALL resshel ( v_cpwwp, mfjm2 )
      CALL timer ( outmod, iotty, "After resshel" )
c
      call atpoint ( "aft-RESSHEL", "jwal1",jwal1,"a", a,
     $        iotty,outmod )
c
 33   continue     !! lrshel from about line 1496

         DEALLOCATE ( grrin )

c.... Write Chi to disk, Vdata (iodsk)
c
      if ( lspark .ne. 0 ) then
         write ( iodsk, 8312 )
 8312    format ( " Writing Chi-i-l " )
         write ( iodsk, 8311 ) jmax1, mth12
 8311    format ( 2i5 )
         do 8320 jwdsk = 1, lmax2
            write ( iodsk, 8313 ) ( grri(iwdsk,jwdsk),iwdsk=1,mth12 )
 8313       format ( 1p10e14.6 )
 8315       continue
 8320    continue
      end if
c
c.....calculate sum of chi over surfaces for sin and cosine
c
      if ( check2 ) write (outmod,215)
 215  format ( /, 3x,"l",5x,"sumpc",5x,"sumps",5x,
     $     "sumwc",5x,"sumws" )
c
      do 220 l1 = 1, jmax1
c
         sumpc = 0.0
         sumps = 0.0
         sumwc = 0.0
         sumws = 0.0
c
         do 216 i = 1, mth
            sumpc = sumpc + grri(i,l1)
            sumps = sumps + grri(i,jmax1+l1)
            sumwc = sumwc + grri(mth+i,l1)
            sumws = sumws + grri(mth+i,jmax1+l1)
 216     continue
c
         ll = l1 - 1 + lmin(1)
         if ( check2 ) write ( outmod,217 ) ll, sumpc,sumps,
     $        sumwc,sumws
 217     format ( 1x, i4, 1p4e12.4 )
c
 220  continue
c
c......store away chi on wall for input to kdis.  both the cos(lth)
c     and the sin(lth) contribution.
c
c      call wwall ( mth1, xwal, zwal )
      mw = mth - (jtop-jbot-1)
c
      do 210 l1 = 1, jmax1
c
         iw = mth - jtop + 2
         icnt = 0
c
         do 200 i = 1, mth
c
            il11 = (l1-1)*mth12 + mth + i
            il12 = jmax1*mth12 + il11
c
            if (  (i.gt.jbot) .and. (i.lt.jtop) ) go to 200
c
            if ( i .eq. jtop ) iw = 1
c
            chiwc(iw,l1) = grri(mth+i,l1)
            chiws(iw,l1) = grri(mth+i,jmax1+l1)
            xpass(iw) = xwal(i)
            zpass(iw) = zwal(i)
c
            iw = iw + 1
            icnt = icnt + 1
c
 200     continue

c.....fill in last point if torus..
c
         if ( jtop .ne. jbot ) go to 210
         chiwc(mw,l1) = chiwc(1,l1)
         chiws(mw,l1) = chiws(1,l1)
         xpass(mw) = xpass(1)
         zpass(mw) = zpass(1)
c
 210  continue
c
      if ( jtop .ne. jbot ) go to 260
c
c      if ( checkd .and. (lfele .eq. 0) ) then
c
c......fourier analyse chi response function
c     for plasma.  use wrkvr and wrkvi as
c     temporary storage.
c
         call fanal1d ( grri,nths2,nfm2,0, wrkvr,wrkvi,nfm,nfm )
c
         call matwrtn(wrkvr,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Chiplar(l_srce,l_obs)", outmod,0 )
         call matwrtn(wrkvi,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Chiplai(l_srce,l_obs)", outmod,0 )
c
c......fourier analyse chi response function for wall.
c
         call fanal1d ( grri,nths2,nfm2,mth, wrkvr,wrkvi,nfm,nfm )
c
         call matwrtn(wrkvr,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Chiwal(l_srce,l_obs)", outmod,0 )
         call matwrtn(wrkvi,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "Chiwali(l_srce,l_obs)", outmod,0 )
c
c      end if
c
 260  continue
c
      if ( check2 )
     $     write ( outmod,8020 ) jtop, jbot, mw, icnt
 8020 format ( 1x,/, "jtop, jbot, mw, icount = ", 4i5, / )
c
      IF ( lfele .ne. 0 ) THEN
         WRITE ( outmod, '(/,20x,"Finite Elements: ARR, AII, etc.",/)' )
c..... Integrate over finite elements for GATO, etc.  .....
c      FACTPI factor included in FELANV
c
      CALL atpoint ( "FELANV","j1", j1,"xinf(1)", xinf(1),
     $     iotty, outmod )
         IF ( lfele .eq. 1 ) THEN
            CALL felanv ( grri,arr, cnqd, 0,0 )
            CALL felanv ( grri,aii, snqd, 0,jmax1 )
            CALL felanv ( grri,ari, snqd, 0,0 )
            CALL felanv ( grri,air, cnqd, 0,jmax1 )
         ENDIF
         IF ( (lfele == 3) .and. (lhighn /= 2) ) THEN
            CALL felanv3 ( grri,arr, cnqd, 0,0 )
            CALL felanv3 ( grri,aii, snqd, 0,jmax1 )
            CALL felanv3 ( grri,ari, snqd, 0,0 )
            CALL felanv3 ( grri,air, cnqd, 0,jmax1 )
         END IF
         IF ( lhighn /= 2 ) THEN
            DO j1 = 1, jmax1
               DO j2 = 1, jmax1
                  vacmat0(j1,j2) = ( arr(j1,j2) + aii(j1,j2) )
                  vacmti(j1,j2) = ( air(j1,j2) - ari(j1,j2) )
               END DO
            END DO
         END IF
         IF ( lhighn == 2 ) THEN
            grri(1:2*mfel,1:2*mfel) =
     $           znorm*twopi*dth*grri(1:2*mfel,1:2*mfel)
            vacmat0(1:mfel,1:mfel) = grri(1:mfel,1:mfel)
            vacmti(1:mfel,1:mfel) = grri(mfel+1:2*mfel,1:mfel)
         END IF
      CALL vacasym ( vacmat0, nfm,jmax1,"Vfinel-13 without l-nq",
     $     outmod,iotty )
         call matwrtn ( vacmti, nfm,nfm, ln,ln,jmax1,jmax1,jdel,jdel,
     $        "vacmti", outmod,0 )
c$$$      CALL vacasym ( vacmti, nfm,jmax1,"Vfineli-13 without l-nq",
c$$$     $     outmod,iotty )
         GO TO 320      !!! About line 2146
      END IF
c
c.....Fourier Analyze over observer points with trapezoidal rule.
c     FACTPI factor included in FORANV
c
      zzfac = dth * twopi   ! This will introduce a factor of 4*pi*pi.
      arr(1:jmax1,1:jmax1) = zzfac *
     $     matmul ( transpose(cslth(1:mth,1:jmax1)),
     $     grri(1:mth,1:jmax1) )
      aii(1:jmax1,1:jmax1) = zzfac *
     $     matmul ( transpose(snlth(1:mth,1:jmax1)),
     $     grri(1:mth,jmax1+1:2*jmax1) )
      ari(1:jmax1,1:jmax1) = zzfac *
     $     matmul ( transpose(snlth(1:mth,1:jmax1)),
     $     grri(1:mth,1:jmax1) )
      air(1:jmax1,1:jmax1) = zzfac *
     $     matmul ( transpose(cslth(1:mth,1:jmax1)),
     $     grri(1:mth,jmax1+1:2*jmax1) )

c      if ( checkd ) then
         call matwrtn ( arr, j1v,j2v, ln,ln,jmax1,jmax1,jdel,jdel,
     $        "arr", outmod,0 )
         call matwrtn ( aii, j1v,j2v, ln,ln,jmax1,jmax1,jdel,jdel,
     $        "aii", outmod,0 )
         call matwrtn ( ari, j1v,j2v, ln,ln,jmax1,jmax1,jdel,jdel,
     $        "ari", outmod,0 )
         call matwrtn ( air, j1v,j2v, ln,ln,jmax1,jmax1,jdel,jdel,
     $        "air", outmod,0 )
c      endif
c
c################################################################
c
      if ( ldisc .eq. 0 ) go to 8909  !!!! About line 2061
c
c...        Calculate Chi discontinuity at "rational" surface.
c           Outer solution already calculated with wall boundary
c           condition. Now calculate interior region by using
c           the interior settings of the wall for the "rational"
c           surface.
c           Can use subroutine chiwal from sprk2. Use cnqd and snqd
c           sources. Store wall response in cwalr, cwali.
c           Use plasma for both observer and source.
c
      call chiwal ( xpla,zpla,xpla,zpla, cslth,snlth,2,2,1,1,
     $     cwalr,cwali )
c... Set j's as if for a wall
c      j1 = 2
c      j2 = 2
c   force ksgn here to simulate the interior problem.
c      ksgn = 1
c      iopw = 1
c      iops = 1
c
c     As above:
c.....Fourier Analyze over observer points with trapezoidal rule.
c     FACTPI factor included in FORANV.
c     Use vacmat, vacmti, rmatr, rmati for storage:
c
      call foran1 ( cwalr,vacmat, cslth, 0,0 )
      call foran1 ( cwali,rmatr,  snlth, 0,0 )
      call foran1 ( cwalr,vacmti, snlth, 0,0 )
      call foran1 ( cwali,rmati,  cslth, 0,0 )
c
      if ( checkd ) then
         call matwrtn ( vacmat, j1v,j2v, ln,ln,jmax1,jmax1,jdel,jdel,
     $        "vacmat", outmod,0 )
         call matwrtn ( rmatr, nfm,nfm, ln,ln,jmax1,jmax1,jdel,jdel,
     $        "rmatr", outmod,0 )
         call matwrtn ( vacmti, j1v,j2v, ln,ln,jmax1,jmax1,jdel,jdel,
     $        "vacmti", outmod,0 )
         call matwrtn ( rmati, nfm,nfm, ln,ln,jmax1,jmax1,jdel,jdel,
     $        "rmati", outmod,0 )
      endif
c
c...Store difference in arr and aii
c   Assume for now that cross terms vanish:
      do j1 = 1, jmax1
         do j2 = 1, jmax1
            arr(j1,j2) = arr(j1,j2) + aii(j1,j2)
     $           - vacmat(j1,j2) - rmatr(j1,j2)
         end do
      end do
c
      call matwrtn ( arr, j1v,j2v, ln,ln,jmax1,jmax1,jdel,jdel,
     $        "d-chi", outmod,0 )
c
c.. Now, invert arr into aii. Use rmatr as workspace:
c
c   Define a unit matrix

      ALLOCATE ( umtrx(jmax1,jmax1) )

      umtrx = 0.0
      do j1 = 1, nfm
            umtrx(j1,j1) = 1.0
      end do

c      call gelimb ( arr,nfm, umtrx,jmax1, jmax1,jmax1, aii,nfm, rmatr,
c     $     ier )
      call f04aaf ( arr,nfm, umtrx,jmax1, jmax1,jmax1, aii,nfm, rmatr,
     $     ier )

      DEALLOCATE ( umtrx )

      write ( outmod, '("ier in del-chi = ", i5)' ) ier
c
      call matwrtn ( aii, j1v,j2v, ln,ln,jmax1,jmax1,jdel,jdel,
     $        "inv-d-chi", outmod,iotty )
      call vacasym ( aii, nfm,jmax1,"Inverse-D-Chi", outmod,iotty )
c
c.. Normalize each row to unity on the diagonal:
c
      do j1 = 1, jmax1
         aj1j1 = aii(j1,j1)
         do j2 = 1, jmax1
            aii(j1,j2) = aii(j1,j2) / aj1j1
         end do
      end do
c
      call matwrtn ( aii, j1v,j2v, ln,ln,jmax1,jmax1,jdel,jdel,
     $        "norm-inv-chi", outmod,iotty )
c
      return
c
 8909 continue
c
c #########################  END OF LDISC ############################
c
      if ( .not. lnova ) go to 1133

c  Allocate some storage for NOVA:

      ALLOCATE ( ajll(jmax1,jmax1) )
c
      do 27 i = 1, jmax1
         do 26 j = 1, jmax1
            rmatr(i,j) = zero
            rmati(i,j) = zero
            ajll(i,j) = zero
 26      continue
 27   continue
c
      do 101 l1 = 1, jmax1
         ll1 = l1 - 1 + lmin(1)
         al1nq = ll1 - nq
         do 100 l2 = 1, jmax1
            ll = l2 - 1 + lmin(1)
            al2nq = ll - nq
            do 80 i = 1, mth
               theta = (i-1)*dth
c     il11 = (l1-1) * mth12 + i
c     il12 = jmax1 * mth12 + il11
c
               ar = ( xjdtxj(i)*grri(i,l1)-
     $              al2nq*grri(i,jmax1+l1) ) * factpi
               ai = ( xjdtxj(i)*grri(i,jmax1+l1)+
     $              al2nq*grri(i,l1) ) * factpi
               ar1 = (ai*cslth(i,l2)-ar*snlth(i,l2))
     $              /xjacob(i)
               ai1 = -(ar*cslth(i,l2)+ai*snlth(i,l2))
     $              /xjacob(i)
c
               rmatr(l2,l1) = rmatr(l2,l1)
     $              + ar1
               rmati(l2,l1) = rmati(l2,l1)
     $              + ai1
c
c......calculate inverse jacobian matrix.
c
               el1l2t = ( ll1-ll ) * theta
               ajll(l2,l1) = ajll(l2,l1)
     $              + cos(el1l2t)
     $              / (twopi*xjacob(i))
c
 80         CONTINUE
c
            ajll(l2,l1) = dth * ajll(l2,l1)
            rmatr(l2,l1) = rmatr(l2,l1)
     $           * dth*al1nq/twopi2
            rmati(l2,l1) = rmati(l2,l1)
     $           * dth*al1nq/twopi2
c
 100     CONTINUE
 101  CONTINUE

 1133 CONTINUE   ! IF ( NOT lnova ) switch
c.........................................................................

      DO j1 = 1, jmax1
         DO j2 = 1, jmax1
            vacmat0(j1,j2) = arr(j1,j2) + aii(j1,j2)
            vacmat00(j1,j2) = arr(j1,j2) + aii(j1,j2)
            vacmt00i(j1,j2) = air(j1,j2) - ari(j1,j2)
            vacmti(j1,j2) = air(j1,j2) - ari(j1,j2)
         END DO
      END DO

c.. Construct the total vacuum matrix for up-down asymmetry in vcmatu0,
c    Unsymmetrized and no l-nq

      vcmattu0(1:jmax1,1:jmax1) = vacmat0(1:jmax1,1:jmax1)
      vcmattu0(jmax1+1:jmax2,jmax1+1:jmax2) = vacmat0(1:jmax1,1:jmax1)
      vcmattu0(1:jmax1,jmax1+1:jmax2) = - vacmti(1:jmax1,1:jmax1)
      vcmattu0(jmax1+1:jmax2,1:jmax1) = + vacmti(1:jmax1,1:jmax1)
c
 320  CONTINUE   ! Comes from IF ( LFELE .NE. 0 ) around  line 1940
c
c$$$      write ( outmod,'(a)' ) char(12)
c$$$      call matwrtn ( vacmat0,nfm,nfm,ln,ln,jmax1,jmax1,8,8,
c$$$     $     "Vacmat0 without l-nq", outmod,iotty )
c$$$      call vacasym ( vacmat0, nfm,jmax1,"Vacmat without l-nq",
c$$$     $     outmod,iotty )
c$$$      call masym0 ( vacmat0, nfm,jmax1,ln,ln, 0,0, 1.0,
c$$$     $     "Vacmat0 without l-nq", outmod,iotty )
c$$$c      call vacasymi ( vacmat0, nfm,jmax1,"Vacmat0 without l-nq",
c$$$c     $     outmod,iotty )

      write ( outmod, '(a)' ) char(12)

      call matwrtn ( vacmat0,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $     "Vacmat0 without l-nq", outmod,iotty )
      IF ( jmax1 > 8)
     $     CALL matwrts ( vacmat0,nfm,nfm,ln,ln, 25,9, 25,9,
     $     1,1,"Vacmat0 without l-nq", outmod,iotty )

      CALL matwrtn ( vcmattu0,nfm2,nfm2,ln,ln,jmax1,jmax1,jdel,jdel,
     $     "Vcmattu0 without l-nq", outmod,iotty )
      IF ( jmax1 > 8 )
     $     CALL matwrts ( vcmattu0,nfm2,nfm2,ln,ln, 25,9, 25,9,
     $     1,1,"Vcmattu0 without l-nq", outmod,iotty )
      CALL vacasyma ( vcmattu0, nfm2,jmax2,1.0, "Vcmattu0",
     $     outmod,iotty )
      CALL masym0 ( vcmattu0, nfm2,jmax2,ln,ln,0,0, 1.0,
     $     "Vcmattu0", outmod,iotty )

c$$$c.. Compare response with LAR cylinder...
c$$$
c$$$      DO jzj = 1, jmax1
c$$$         lzl = lmin(1) - 1 + jzj
c$$$         chlarc = 1.0 / ABS(lzl) / xpl
c$$$         chlarc = twopi2 * chlarc
c$$$         IF ( chlarc < 1.0E-60 ) chlarc = 1.0E-20
c$$$         dchlarc = ( chlarc - vcmattu0(jzj,jzj) ) / chlarc
c$$$         WRITE ( outmod, '( 5x, "D - LAR Cylinder = ", ES14.7 )' )
c$$$     $        dchlarc
c$$$         WRITE ( iotty,  '( 5x, "D - LAR Cylinder = ", ES14.7 )' )
c$$$     $        dchlarc
c$$$      END DO

      IF ( checkd ) THEN
         ALLOCATE ( wrkvr2(nfm2,nfm2) )
         CALL mateig2a ( vcmattu0, nfm2,jmax2, wrkvr2,work0,
     $        ln,lx,ff, 0, 2, jobid, "Vcmattu0", outmod )
         DEALLOCATE ( wrkvr2 )
      END IF


      ALLOCATE ( zwk(jmax1+5) )
      DO izpl = 1, 4
         izplz = (izpl-1)*(jmax1/4) + 1
         zwk(1:jmax1) = vacmat0(1:jmax1,izplz)
         CALL pospl1 ( pigrd,zwk, 1,jmax1, izpl,"Vacmat0", izplz )
      END DO
      CALL framep(jobid, ff)
      DEALLOCATE ( zwk )
         
      call matwrt9 ( vacmat0, nfm,nfm,ln,lx, "Vacmat0 without  l-nq",
     $     outmod,iotty )
c
      if ( .not. checkd ) go to 105
      call matwrtn ( arr, j1v,j2v, ln,ln,jmax1,jmax1,jdel,jdel,
     $     "arr", outmod,0 )
      call matwrtn ( aii, j1v,j2v, ln,ln,jmax1,jmax1,jdel,jdel,
     $     "aii", outmod,0 )
      call matwrtn ( ari, j1v,j2v, ln,ln,jmax1,jmax1,jdel,jdel,
     $     "ari", outmod,0 )
      call matwrtn ( air, j1v,j2v, ln,ln,jmax1,jmax1,jdel,jdel,
     $     "air", outmod,0 )
 105  continue
c
      write ( iotty,  '(/,1x," n, mth, lmin, lmax = ",f4.1, 3i5)')
     $     n, mth, lmin(1), lmax(1)
      write ( outmod, '(/,1x," n, mth, lmin, lmax = ",f4.1, 3i5)')
     $     n, mth, lmin(1), lmax(1)
c
      IF ( (lgato .eq. 1) .OR. (lgato .eq. 3) ) THEN

         zlabel =
     $     datimv//": GATOVAC-13, ID= "//jobid(1:lenj)
         call gatonorm ( vacmat0,nfm, gatovac,nfm, rgato,mfel,mth,
     $        qa1,twopi,taugam, vtog, zlabel )
         write ( outmod, '("Conversion factor = ", 1pe12.5)') vtog
         write ( iotty,  '("Conversion factor = ", 1pe12.5)') vtog
         call matwrtn ( gatovac,nfm,nfm,ln,ln,jmax1,jmax1,8,8,
     $        "GATOVAC-13", outmod,iotty )
         call matwrtx ( gatovac,nfm,nfm,1,8,1,8,
     $        "GATOVAC-13", outmod,iotty )

         call vacasym ( gatovac, nfm,jmax1,"GATOVAC-13 without l-nq",
     $        outmod,iotty )

c  Allocate storage
c
c   Assume for now for lhighn=1 that the imaginary part of W_vac is zero:

         IF ( lhighn > 0 ) THEN

            ALLOCATE ( workg1(mfel2,mfel2), workg2(mfel2,mfel2) )

            IF ( lhighn == 1 ) THEN
               DO i = 1, mfel
                  DO j = 1, mfel
                     cnqth =  cos ( n*(enqth(i) - enqth(j)) )
                     snqth =  sin ( n*(enqth(i) - enqth(j)) )
                     workg1(i,j) =           gatovac(i,j) * cnqth
                     workg1(i,j+mfel) =      gatovac(i,j) * snqth
                     workg1(i+mfel,j) =    - gatovac(i,j) * snqth
                     workg1(i+mfel,j+mfel) = gatovac(i,j) * cnqth
                  END DO
               END DO
            END IF

c... Use GRRI here to calculate the vacuum matrix. And multiply
c     by the conversion factor, VTOG, which was calculated above.

            IF ( lhighn == 2 ) THEN
               workg1(1:2*mfel,1:2*mfel) =
     $              vtog * grri(1:2*mfel,1:2*mfel)
            END IF

            CALL matwrtn ( workg1,mfel2,mfel2,1,1,mfel2,mfel2,8,8,
     $           "V-HIN", outmod,iotty )
            CALL matwrtx ( workg1,mfel2,mfel2,1,8,1,8,
     $           "V-HIN", outmod,iotty )
            CALL vacasym ( workg1, mfel2,mfel2,"V-HIN",
     $           outmod,iotty )

            WRITE (  IOTTY,
     $      '("Writing to 46, VACGAHIN:  mfel,mth, rgato, n,qa1 = ",/,
     $           2i4, 1p3e13.5 /  )' ) mfel,mth, rgato, n,qa1
            WRITE (  OUTMOD,
     $      '("Writing to 46, VACGAHIN:  mfel,mth, rgato, n,qa1 = ",/,
     $           2i4, 1p3e13.5 /  )' ) mfel,mth, rgato, n,qa1

            OPEN ( UNIT=46,FILE='vacgahin' )

            WRITE ( 46, '("VACGAHIN: rgato, mfel,mth, qa1 = ",/,
     $           2i4, 1p3e13.5 / )' )  mfel,mth, rgato, n,qa1

            DO m1 = 1, mfel2
               WRITE ( 46,'(/, i5 )') m1
               WRITE ( 46,'(10e15.8)') (workg1(m1,m2), m2 = 1, mfel2)
            END DO

            CLOSE ( UNIT=46 )
c
            write ( outmod, '(/,20x,"lgato = 3, GATOVAC with nqth:")' )
            call mateig2a ( workg1,mfel2,mfel2, workg2,work0,
     $           1,mfel2,ff, 0, lchecke, jobid, outmod )
c$$$  call mateig2 ( workg1, mfel2,mfel2, workg2,work0,work,work1,
c$$$     $           1,mfel2,ff, 0, lchecke, wkxx,wkyy, jobid, outmod )

            DEALLOCATE ( workg1, workg2 )

         END IF  ! if ( lhighn > 0 )
      END IF  ! if (lgato .eq. 1 .or. lgato .eq. 3 )
c
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      if ( check1 )
     .     call timer ( outmod, iotty, "end of vacmat" )
c
c Transform the Fourier basis to finite element basis.
c Use AIR, ARI as work space. Reconstruct VACMAT w/o l-nq in RMATR
c Careful that RMATR, RMATI are  being overwritten.
c
      if ( lgato .eq. 2 )  then
c
c Check orthogonality of the T matrices.

         lnf = lmin(1)
         lxf = lmax(1)

         IF ( ltfil > 0 ) THEN
            lnf = - mfel/2
            lxf = + mfel/2
         END IF

         jmax1f = lxf - lnf + 1

         call orchek ( mfel, jmax1f)
c
c Get vacuum matrix without l-nq.
         rmatr(1:jmax1,1:jmax1) =
     $        arr(1:jmax1,1:jmax1) + aii(1:jmax1,1:jmax1)
         call matwrtn ( rmatr,nfm,nfm,1,1,jmax1,jmax1,8,8,
     $        "RMATR", outmod,iotty )
         call matwrtx ( rmatr,nfm,nfm,1,8,1,8,
     $        "RMATR", outmod,iotty )
c
c.... Symmetrize:
         if (lsymz ) then
            rmatr(1:jmax1,1:jmax1) =
     $           0.5* ( rmatr(1:jmax1,1:jmax1) +
     $           transpose(rmatr(1:jmax1,1:jmax1)))
         end if
         write ( outmod, '(/,20x, "lgato = 2, Fourier Matrix:")' )
         call mateig2a ( rmatr,nfm,jmax1, wrkvr,work0,
     $        ln,lx,ff, 0, lchecke, jobid, outmod )
c$$$         call mateig2 ( rmatr, nfm,jmax1, wrkvr,work0,work,work1,
c$$$     $        ln,lx,ff, 0, lchecke, wkxx,wkyy, jobid, outmod )
c
c$$$         iopc = 1
c$$$         iops = 2
c$$$         call fotofi ( rmatr, rmati,mfel,mfel,jmax1, sinlt,nths, iops )
c$$$         call matwrtn ( rmati,nfm,nfm,1,1,mfel,mfel,8,8,"VFinels",
c$$$     $        outmod,iotty )
c$$$         call fotofi ( rmatr, rmatr,mfel,mfel,jmax1, coslt,nths, iopc )
c$$$         call matwrtn ( rmatr,nfm,nfm,1,1,mfel,mfel,8,8,"VFinelc",
c$$$     $        outmod,iotty )
c$$$c
c$$$         rmatr(1:mfel,1:mfel) =
c$$$     $        rmatr(1:mfel,1:mfel) + rmati(1:mfel,1:mfel)
c
         call new_basis ( rmatr, nfm, mfel,jmax1,ltfil,
     $        sinlt,coslt,nths,nfm,  outmod,iotty )

         call matwrtn ( rmatr,nfm,nfm,1,1,mfel,mfel,8,8,"VFinel",
     $        outmod,iotty )
         call matwrtx ( rmatr,nfm,nfm,1,8,1,8, "VFinel",
     $        outmod,iotty )
c
         call vacasym ( rmatr, nfm,mfel,"VFinel-2", outmod,iotty )
c
         zlabel =
     $     datimv//": GATOVAC, ID= "//jobid(1:lenj)
         call gatonorm ( rmatr,nfm, gatovac,nfm, rgato,mfel,mth,
     $        qa1,twopi,taugam, vtog, zlabel )
         write ( outmod, '("Conversion factor = ", 1pe12.5)') vtog
         write ( iotty,  '("Conversion factor = ", 1pe12.5)') vtog
         call matwrtn ( gatovac,nfm,nfm,1,1,mfel,mfel,8,8,
     $        "GATOVAC", outmod,iotty )
         call matwrtx ( gatovac,nfm,nfm,1,8,1,8,
     $        "GATOVAC", outmod,iotty )
      end if
c
      write ( outmod, 500 ) n,q, nj,mj,lj
 500  format (//,1x,'n,q, nj,mj,lj = ',1p2e12.4,3i5,/ )
      if ( (.not. lpest1) .and. check2 )
     $     write ( outmod, 501 )  ( delta(i),i=1,mth1 )
 501  format ( 1x, " delta =",/, (1x,10e11.4) )
c
 630  continue
c
      if ( lchecke .gt. 0 ) then
         if ( lgato .eq. 1 ) then
            write ( outmod, '(/,20x, "lgato = 1:")' )
            call mateig2a ( vacmat0, nfm,mfel, wrkvr,work0,
     $           ln,lx,ff, 0, lchecke, jobid, outmod )
c$$$            call mateig2 ( vacmat0, nfm,mfel, wrkvr,work0,work,work1,
c$$$     $           ln,lx,ff, 0, lchecke, wkxx,wkyy, jobid, outmod )
         end if
         if ( lgato .eq. 2 ) then
            write ( outmod, '(/,20x, "lgato = 2: RMATR:")' )
            call mateig2a ( rmatr,nfm,mfel, wrkvr,work0,
     $           ln,lx,ff, 0, lchecke, jobid, outmod )
c$$$            call mateig2 ( rmatr, nfm,mfel, wrkvr,work0,work,work1,
c$$$     $            ln,lx,ff, 0, lchecke, wkxx,wkyy, jobid, outmod )
         end if
         if ( lgato .eq. 2 ) then
            write ( outmod, '(/,20x, "lgato = 2: GATOVAC:")' )
            call mateig2a ( gatovac,nfm,mfel, wrkvr,work0,
     $           ln,lx,ff, 0, lchecke, jobid, outmod )
c$$$            call mateig2 ( gatovac, nfm,mfel, wrkvr,work0,work,work1,
c$$$     $           ln,lx,ff, 0, lchecke, wkxx,wkyy, jobid, outmod )

         end if
         IF ( (lgato == 3) .AND. (lhighn==0) ) THEN
            WRITE ( OUTMOD, '(/,20x, "lgato = 3: GATOVAC")' )
            CALL mateig2a ( gatovac,nfm,mfel, wrkvr,work0,
     $           ln,lx,ff, 0, lchecke, jobid, outmod )
c$$$            call mateig2 ( gatovac, nfm,mfel, wrkvr,work0,work,work1,
c$$$     $           ln,lx,ff, 0, lchecke, wkxx,wkyy, jobid, outmod )
         END IF
c$$$         if ( lgato .eq. 0 ) then
c$$$            write ( outmod, '(/,20x, "Fourier Without l-nq :")' )
c$$$            call mateig2 ( vacmat0, nfm,jmax1, wrkvr,work0,work,work1,
c$$$     $           ln,lx,ff, 0, lchecke, wkxx,wkyy, jobid, outmod )
c$$$         end if
      end if
c
      if ( lgato .ne. 0 ) go to 4111
c
c$$$c.........find the difference between the two methods.
c$$$c
c$$$      do 551 l1 = 1, jmax1
c$$$         do 550 l2 = 1, jmax1
c$$$            ajll(l1,l2) = ( rmati(l1,l2)-rmatr(l1,l2) )
c$$$c......switch symmetry in rmatr.
c$$$c     rmati(l1,l2) = rmatr(l2,l1)
c$$$ 550     continue
c$$$ 551  continue
c
c.....calculate the vacuum matrix.
c
      do 3012 j1 = 1, jmax1
c
         l1 = j1 + lmin(1) - 1
         alnq1 = l1 - nq
c
         do 3011 j2 = 1, jmax1
c
            l2 = j2 + lmin(1) - 1
            alnq2 = l2 - nq
c
            vacmat(j1,j2) = alnq1*alnq2 * vacmat0(j1,j2)
            vacmti(j1,j2) = alnq1*alnq2 * vacmti(j1,j2)
c
c......Store away unsymmetric form of Vacmat.
            vacmatu(j1,j2) = vacmat(j1,j2)
            vacmtiu(j1,j2) = vacmti(j1,j2)
c
 3011    continue
 3012 continue

c.. Construct the total vacuum matrix for up-down asymmetry in vcmattu,
c    With  l-nq

      vcmattu(1:jmax1,1:jmax1) = vacmatu(1:jmax1,1:jmax1)
      vcmattu(jmax1+1:jmax2,jmax1+1:jmax2) = vacmatu(1:jmax1,1:jmax1)
      vcmattu(1:jmax1,jmax1+1:jmax2) = -vacmtiu(1:jmax1,1:jmax1)
      vcmattu(jmax1+1:jmax2,1:jmax1) = + vacmtiu(1:jmax1,1:jmax1)
c
      IF ( LNOVA ) THEN
c
         call matwrtn(rmati,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "rmati", outmod,0 )
c
c.......calculate p-xi matrix.
c
         call matmul0 ( ajll,vacmat, jmax1,nfm,jmax1, rmati,nfm )
c
c.......normalize.
c
         do 401 l1 = 1, jmax1
            do 400 l2 = 1, jmax1
               rmati(l1,l2) = rmati(l1,l2) / twopi2
 400        continue
 401     continue
c
         call matwrtn ( ajll,jmax1,jmax1,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "inverse j matrix", outmod,0 )
         call matwrtn ( rmati,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "p-xi matrix ", outmod,0 )
c
      end if
c
      write ( outmod, '(a)' ) char(12)
      call matwrtn ( vacmat,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $     "Vacmat with l-nq", outmod,iotty )
      call matwrt9 ( vacmat, nfm,nfm,ln,lx, "Vacmat with l-nq",
     $     outmod,iotty )
      call matwrtn(vacmti,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $     "vacmti", outmod,0 )
c
c.......calculate asymmetries...
c
      call vacasyma ( vacmat, nfm,jmax1,1.0, "Vacmat", outmod,iotty )
      call masym0 ( vacmat, nfm,jmax1,ln,ln,0,0, 1.0,
     $     "Vacmat", outmod,iotty )
      call vacasyma ( vacmti, nfm,jmax1,-1.0, "Vacmti", outmod,iotty )
      call masym0 ( vacmti, nfm,jmax1,ln,ln,0,0, -1.0,
     $     "Vacmti", outmod,iotty )

c      call vacasymi ( vacmat, nfm,jmax1,"Vacmat", outmod,iotty )
c
      if ( lchecke .gt. 0  ) then
         write ( outmod, '(/,20x, "Fourier With l-nq :")' )
c$$$            call mateig2 ( vacmat,nfm,jmax1, wrkvr,work0,work,work1,
c$$$     $        ln,lx, ff, 2, lchecke, wkxx,wkyy, jobid, outmod )

            call mateig2a ( vacmat,nfm,jmax1, wrkvr,work0,
     $        ln,lx, ff, 2, lchecke, jobid, outmod )
      end if
c
 4111 continue
c
c  Write vacuum matrix for stability codes. Symmetrize if needed.
c
      IF ( lsymz ) THEN
c
c.....symmetrize.
c
         do 601 l1 = 1, jmax1
            do 600 l2 = l1, jmax1
               vacmat(l1,l2) = 0.5 * ( vacmat(l1,l2)+vacmat(l2,l1) )
               vacmti(l1,l2) = 0.5 * ( vacmti(l1,l2)-vacmti(l2,l1) )
               rmatr(l1,l2) = 0.5 * ( rmatr(l1,l2)+rmatr(l2,l1) )
 600        continue
 601     continue
         do 621 l1 = 1, jmax1
            do 620 l2 = l1, jmax1
               vacmat(l2,l1) = vacmat(l1,l2)
               IF ( l1 /= l2 ) vacmti(l2,l1) = - vacmti(l1,l2)
               IF ( L1 == L2 ) vacmti(l2,l1) = 0.0
               rmatr(l2,l1) = rmatr(l1,l2)
 620        continue
 621     continue
c
      write ( outmod, '(a)' ) char(12)
      call matwrtn ( vacmat,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $     "Vacmat with l-nq, sym", outmod,iotty )
      call matwrt9 ( vacmat, nfm,nfm,ln,lx, "Vacmat with l-nq, sym",
     $     outmod,iotty )
      call matwrtn(vacmti,j1v,j2v,ln,ln,jmax1,jmax1,jdel,jdel,
     $     "vacmti, sym", outmod,0 )

      END IF
c
      if ( ladj .eq. 1 ) then
c
         open ( 60, file='vacadj', status='new', form='formatted' )
         nadj = n + 0.001
         write ( 60, 554 ) mthin1,lmin(1),lmax(1),nadj,qa1
 554     format ( 4i5, e13.5 )
         do 556 j1 = 1, jmax1
            write ( 60, 555 )  ( vacmat(j1,j2), j2=1,jmax1 )
 555        format ( 10e13.5 )
 556     continue
c
         close ( 60 )
         go to 710
c
      end if
c
      if ( ldcon .eq. 1 ) then
c
         open ( 62, file='vacdcon', status='new', form='formatted' )
         ndcon = n + 0.001
         write ( 62,'("mthin1, lmin, lmax, n, q = ")' )
         write ( 62, '(4i5, e15.7)' ) mthin1,lmin(1),lmax(1),ndcon,qa1
         write ( 62, '("Real Vacmat:")' )
         do 557 j1 = 1, jmax1
            write ( 62,'(10e15.7)' )  ( vacmat(j1,j2), j2=1,jmax1 )
 557     continue
         write ( 62, '("Imag. Vacmat:")' )
         do 558 j1 = 1, jmax1
            write ( 62,'(10e15.7)' )  ( vacmti(j1,j2), j2=1,jmax1 )
 558     continue
c
         close ( 62 )
         go to 710
c
      end if
c
      if ( lnova ) then
c
         call matwrtn(rmatr,nfm,nfm,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "rmatr", outmod,0 )
         call matwrtn ( ajll,jmax1,jmax1,ln,ln,jmax1,jmax1,jdel,jdel,
     $        "pll-mc minus pll-fc ", outmod,0 )
c
         mtots = mtot**2
         len=mtots+5
         ndsk=1
         call zop (iovac,"vacout",len,ndsk,iiff,999)
         lgivup=1
         call zwr(iovac,rmatr,mtots,1,lgivup,999)
         call zcl ( iovac, 999 )

         DEALLOCATE ( ajll )

         GO TO 710
c
      end if
c
c.....write vacmat to disk for the pest code.
c
      call atpoint ( "WTOPEST1","ln", ln,"xwal(1)", xwal(1),
     $     iotty, outmod )
      j12 = 1
c
      do 701 j2 = 1, jmax1
         do 700 j1 = 1, jmax1
            vacpstr(j12) = vacmat(j1,j2)
            vacpsti(j12) = vacmti(j1,j2)
            j12 = j12 + 1
 700     continue
 701  continue
c
      IF ( lzio==1 ) THEN
         CALL atpoint ( "WTOPEST2","ln", ln,"zwal(1)", zwal(1),
     $        iotty, outmod )
         CALL wtopest ( vacpstr, vacpsti, xwal, zwal )
      END IF
c
 710  continue
c
      if ( check2 ) then
      write ( outmod, '(2x,"I", 6x,"xp",8x,"zp",8x,"xw",8x,"zw",8x,
     $        "xpp",7x,"zpp",7x,"xwp",7x,"zwp =")' )
      write ( iotty,  '("I, xp, zp, xw, zw, xpp, zpp, xwp, zwp =")' )
      do i = 1, mth1, 8
      write ( outmod, '(i4, 1x, 8f10.3)' )
     $        i, xpla(i),  zpla(i),  xwal(i),  zwal(i),
     $           xplap(i), zplap(i), xwalp(i), zwalp(i)
      write ( iotty,  '(i3, 1x, 8f8.2)' )
     $        i, xpla(i),  zpla(i),  xwal(i),  zwal(i),
     $           xplap(i), zplap(i), xwalp(i), zwalp(i)
      end do
      end if

      IF ( lreshel == 1 ) THEN
c     dont forget to turn it back on if turned off 5/25/01
         CALL timer( outmod, iotty, "Before resshel_1" )
         CALL resshel_1 ( v_cpwwp, mfjm2 )
         CALL timer ( outmod, iotty, "After resshel_1" )
      END IF
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c      DEALLOCATE ( gwpi )
      DEALLOCATE ( rmatr, rmati )
      IF ( lreshel == 1 ) DEALLOCATE ( v_cpwwp )
c
c     if ( checkd ) call exit
c
      return
 999  call errmes ( outpest,'vacuum' )
      end

c...................................................................
      SUBROUTINE regrupk ( grdgre, nd1,nd2, mfel )
c....................................................................

c..   Put the K kernel subgroups in the correct position for the
c     lhighn = 2 option

      REAL, DIMENSION(nd1,nd2) :: grdgre
c      REAL, DIMENSION(:,:), ALLOCATABLE :: zwork

c     The following indices will shift the elements of the kernels the
c     correct amount for storage.

      DO j1 = 1, 2
         DO j2 = 1, 2

            j1f = 2*(j1-1) * mfel
            j2f = 2*(j2-1) * mfel
            j1f1 = j1f + 1
            j2f1 = j2f + 1
            j1ff = j1f + mfel
            j2ff = j2f + mfel
            j1ff1 = j1ff + 1
            j2ff1 = j2ff + 1
            j1fff = j1ff + mfel
            j2fff = j2ff + mfel

c$$$            grdgre(j1f1:j1ff,j2f1:j2ff) =
c$$$     $           grdgre(j1f1:j1ff,j2f1:j2ff) +
c$$$     $           grdgre(j1ff1:j1fff,j2ff1:j2fff)
            grdgre(j1ff1:j1fff,j2ff1:j2fff) =
     $           grdgre(j1f1:j1ff,j2f1:j2ff)
c$$$            grdgre(j1f1:j1ff,j2ff1:j2fff) =
c$$$     $           grdgre(j1ff1:j1fff,j2f1:j2ff) -
c$$$     $           grdgre(j1f1:j1ff,j2ff1:j2fff)
            grdgre(j1ff1:j1fff,j2f1:j2ff) = -
     $           grdgre(j1f1:j1ff,j2ff1:j2fff)

         END DO
      END DO

      END SUBROUTINE regrupk

c.......................................................................
      SUBROUTINE dqdtw ( crsp,nd1,nd2, mobs,msrc, zn, zind,
     $     zx,zxp,zzp,lopt, nout1 )
c.......................................................................

c.. Calculates the response of dQ/dt at a resistive shell. Basically, the
c    Surface Laplacian of the scalar potential.
c   lopt == 3 is for tophat finite element representation
c   lopt == 2 is for the collacation response.

      REAL, DIMENSION(*) :: zind, zx, zxp, zzp
      REAL, DIMENSION(nd1,nd2) :: crsp
      REAL, DIMENSION(:,:), ALLOCATABLE :: zrsp
      REAL, DIMENSION(:), ALLOCATABLE :: zxz, zxxz1, zxxz2, zwk, zdwk,
     $     zrsp1
      INTEGER :: nzn, mobs, msrc, m1tst, m2tst, lopt
      REAL bw_test, znrm, zr, za, zarg, zsarg, zchim, zpi
      CHARACTER(1), PARAMETER :: tab = CHAR(9)

      mobs1 = mobs + 1
      msrc1 = msrc + 1
      mobs4 = mobs + 4
      msrc4 = msrc + 4
      zpi = ACOS(-1.0)

c.. For finite elements. NEGATE THE SQRT(MFEL) FACTOR IN FELANG3. 
c   THEN THE SOURCE AND OBSERVER CAN BE IN k OR theta_k space.      

      IF ( lopt==3 )crsp = crsp / SQRT(FLOAT(msrc))

c.. m1tst and m2tstr are fourier arguments to the test function
      m1tst = 0
      m2tst = 20
      WRITE ( 6,     '(/,5x,"Testing with m1ts, m2tst = ", 2i5)' )
     $     m1tst, m2tst
      WRITE ( nout1, '(/,5x,"Testing with m1ts, m2tst = ", 2i5)' )
     $     m1tst, m2tst

      nzn = zn + 0.01

      ALLOCATE ( zxz(mobs1+5), zxxz1(mobs1+5), zxxz2(mobs1+5),
     $     zwk(mobs+5), zdwk(mobs+5) )
      ALLOCATE ( zrsp(mobs4,msrc4 ), zrsp1(mobs4) )

      zxz(1:mobs1) = SQRT ( zxp(1:mobs1)**2 + zzp(1:mobs1)**2 )
      zxxz1(1:mobs1) = zx(1:mobs1) / zxz(1:mobs1)
      zxxz2(1:mobs1) = 1.0 / ( zx(1:mobs1)*zxz(1:mobs1) )
      
      CALL pospl1 ( zind,zzp,  1,mobs1, 1,"zzp",   mobs1 )
      CALL pospl1 ( zind,zxz,  1,mobs1, 2,"zxz",   mobs1 )
      CALL pospl1 ( zind,zxxz1,1,mobs1, 3,"zxxz1", mobs1 )
      CALL pospl1 ( zind,zxxz2,1,mobs1, 4,"zxxz2", mobs1 )
      CALL framep(jobid, ff)

c... Do a test case on an example input field, bw_test:
c       Integrate over the source first:

c   Normalize the test case to:
c.. The LAR factor of the chi response, zchim, to a single harmonic:
c     Get R and a first

      zr = ( zx(1)+zx(mobs/2+1) ) / 2.0
      za = ( zx(1)-zx(mobs/2+1) ) / 2.0

      zarg = zpi * m2tst / msrc
      zsarg = ( SIN(zarg) )**2 / (msrc*zarg**2)
c.. Assume B(t_k) const over \tau(t): zsarg = 1
      zsarg = 1.0
      IF ( lopt == 3 ) zchim = zsarg / (zr*m2tst) 
      IF ( lopt == 2 ) zchim = 1.0 / (zr*m2tst) 
      
c... For a  LAR circle, dqdt should scale like 1/(m*Ra^2)[m**2 + (na/R)**2]
c      for single harmonic and no wall.:

      znrm = - ( m2tst**2 + (nzn*za/zr)**2 ) / ( m2tst*zr*za**2 )
      IF ( lopt == 3 ) znrm = zsarg * znrm

      WRITE ( nout1, '(20x, "zsarg, znrm, zchim = ", 3es12.4)' )
     $     zsarg, znrm, zchim
      WRITE ( 6, '(20x, "zsarg, znrm, zchim = ", 3es12.4)' )
     $     zsarg, znrm, zchim

      zwk(1:msrc1) = (/ (bw_test(zind(i),m1tst,m2tst),i=1,msrc1) /)

      zrsp1(1:mobs) = MATMUL ( crsp(1:mobs,1:msrc),zwk(1:msrc) )
      zrsp1(mobs1) = zrsp1(1)

c.  Normalize zrsp for testing. Unnormalize later for further calculations.
      zrsp1 = zrsp1 / zchim

      CALL vecwrt ( msrc1, zwk, "bw_test", 1, msrc1, nout1, 6 )
      CALL vecwrt ( mobs1, zrsp1, "chi", 1, mobs1, nout1, 6 )
      CALL pospl1 ( zind,zwk,  1,msrc1, 1,"Bw_test", msrc1 )
      CALL pospl1 ( zind,zrsp1, 1,mobs1, 2,"Chi",     mobs1 )

c.. Subtract LAR analytic solution:
      zwk(1:mobs1) = (/ (bw_test(zind(jj),m1tst,m2tst),jj=1,mobs1) /)
      zdwk(1:mobs1) = zrsp1(1:mobs1) - zwk(1:mobs1)
      CALL vecwrt ( mobs1, zdwk, "chi - chi_an", 1, mobs1, nout1, 6 )
      CALL pospl1 ( zind,zdwk,  1,mobs1, 3,"Chi - Chi_an", mobs )

c.. Unnormalize zrsp1 before calculating dqdtw
      zrsp1 = zrsp1 * zchim

      zwk(1:mobs1) = zrsp1(1:mobs1)
      CALL difspl ( mobs, zind, zwk, zdwk )
      zdwk(mobs1) = zdwk(1)
      zwk(1:mobs1) = zxxz1(1:mobs1) * zdwk(1:mobs1)
      CALL difspl ( mobs, zind, zwk, zdwk )
      zwk(1:mobs) = zxxz2(1:mobs) * zdwk(1:mobs)
     $        - (zn/zx(1:mobs))**2 * zrsp1(1:mobs)
      zwk(mobs1) = zwk(1)
      zwk = zwk / znrm
      CALL vecwrt ( mobs1, zwk, "dqdtw_a", 1, mobs1, nout1, 6 )
      CALL pospl1 ( zind,zwk, 1,mobs1, 4,"Dqdtw_a", mobs )
      CALL framep(jobid, ff)
      
c... End of test case.

      DO i = 1, msrc
         zwk(1:mobs) = crsp(1:mobs,i)
         zwk(mobs1) = zwk(1)
         CALL difspl ( mobs, zind, zwk, zdwk )
         zdwk(mobs1) = zdwk(1)
         IF ( i == 1 .or. i == mobs/2 + 1 ) THEN
            CALL pospl1 ( zind,zwk,  1,mobs1, 1,"crsp",   i )
            CALL pospl1 ( zind,zdwk,  1,mobs1, 2,"dcrsp", i )
         END IF
         zwk(1:mobs1) = zxxz1(1:mobs1) * zdwk(1:mobs1)
         CALL difspl ( mobs, zind, zwk, zdwk )
         zrsp(1:mobs,i) = zxxz2(1:mobs) * zdwk(1:mobs)
     $        - (zn/zx(1:mobs))**2 * crsp(1:mobs,i)
         zrsp(mobs1,i) = zrsp(1,i)
         IF ( i == 1 .or. i == msrc/2 + 1 ) THEN
            zdwk(mobs1) = zdwk(1)
            zrsp1(1:mobs1) = zrsp(1:mobs1,i)
            CALL pospl1 ( zind,zdwk, 1,mobs1, 3,"2nd-der", i )
            CALL pospl1 ( zind,zrsp1, 1,mobs1, 4,"dqdt-resp",    i )
            CALL framep(jobid, ff)
         END IF
      END DO

      DO i = 1, 4
         zrsp1(1:mobs1) = zrsp(1:mobs1,i)
         CALL pospl1 ( zind,zrsp1, 1,mobs1, i,"dqdt-resp",    i )
      END DO
         CALL framep(jobid, ff)

c.. Now Integrate over bw_test to compare with above:

      zwk(1:msrc1) = (/ (bw_test(zind(jj),m1tst,m2tst),jj=1,msrc1) /)
      zdwk(1:mobs1) =
     $     MATMUL ( zrsp(1:mobs1,1:msrc),zwk(1:msrc) )
      zdwk = zdwk / znrm
      CALL vecwrt ( msrc1, zwk, "bw_test", 1, msrc1, nout1, 6 )
      CALL vecwrt ( mobs1, zdwk, "dqdtw_b", 1, mobs1, nout1, 6 )
      CALL pospl1 ( zind,zdwk, 1,mobs1, 1,"Dqdtw_b", mobs )

c.. Subtract LAR analytic solution:
      zwk(1:mobs1) = (/ (bw_test(zind(jj),m1tst,m2tst),jj=1,mobs1) /)
      zdwk(1:mobs1) = zdwk(1:mobs1) - zwk(1:mobs1)
      CALL vecwrt ( mobs1, zdwk, "Dq_b - Dq_b_an", 1, mobs1, nout1, 6 )
      CALL pospl1 ( zind,zdwk, 1,mobs1, 2,"Dq_b - Dq_b_an", mobs )
      CALL framep(jobid, ff)
            
c.. Since NIMROD uses Q dot n rather than the "natural norm", multiply
c   the source points by J.grad r:

      DO i = 1, msrc1
         zrsp(1:mobs1,i) = zrsp(1:mobs1,i) / zxxz2(i)
      END DO

      DEALLOCATE ( zxz, zxxz1, zxxz2, zwk, zdwk, zrsp1 )

      OPEN ( UNIT=76, FILE='dqdt_out',
     $     STATUS='REPLACE', FORM='FORMATTED' )
      
      WRITE ( 76, '(/,20x,"DQDT response",/)' )
      WRITE ( 76, '(5x, i5, a )' ) mobs1,
     $     tab//tab//"mobs1"//tab//"[No. of observer points]"
      WRITE ( 76, '(5x, i5, a )' ) msrc1,
     $     tab//tab//"msrc1"//tab//"[No. of source points]"

      WRITE ( 76, '(5x, i5, a )' ) nzn,
     $     tab//tab//"n"//tab//"Toroidal mode number"
      WRITE ( 76, '(5x, /a/ )' )  "Reponse Matrix(obs,srce):"

c.. Periodicity:
      zrsp(1:mobs1,msrc1) = zrsp(1:mobs1,1)

c.. Reverse the direction of source and observer points:
      zrsp(1:mobs1,1:msrc1) = zrsp(mobs1:1:-1,msrc1:1:-1)

c.. For NIMROD, pass on only half of the boundary source points.

      zrsp(:,1) = 0.5 * zrsp(:,1)
      zrsp(:,msrc1) = 0.5 * zrsp(:,msrc1)

      DO i = 1, mobs1
         WRITE ( 76, '(10es15.7)' ) ( zrsp(i,j), j=1,msrc1 )
      END DO 

      CALL matwrtn ( zrsp,mobs4,msrc4,1,1,mobs,msrc,16,8,
     $        "DQDTW-response matrix", nout1, 0 )
      CALL matwrts ( zrsp,mobs4,msrc4, 1,1, 25,9, 25,9,
     $     1,1,"DQDT-response matrix", nout1,0 )

      CLOSE ( UNIT = 76 )
      DEALLOCATE ( zrsp )
      
      END SUBROUTINE dqdtw
c......................................................................
c.......................................................................
      SUBROUTINE dqdtwl ( crsp,nd1,nd2, mobs,msrc, zn, zind,
     $     zx,zxp,zzp, lmnf, lmxf, nout1 )
c.......................................................................
      
      REAL, DIMENSION(*) :: zind, zx, zxp, zzp
      REAL, DIMENSION(nd1,nd2) :: crsp
      REAL, DIMENSION(:), ALLOCATABLE :: zxz, zxxz1, zxxz2
      REAL, DIMENSION(:,:), ALLOCATABLE :: zclft, zslft, zwklc, zwkls,
     $     zrspc, zrsps
      REAL, DIMENSION(:), ALLOCATABLE :: zwkc, zwks, zdwkc, zdwks
      INTEGER :: nzn, mobs, msrc
      CHARACTER(1), PARAMETER :: tab = CHAR(9)

      zpi = ACOS(-1.0)
      ztpii = 1.0 / ( 2.0*zpi )
      mobs1 = mobs + 1
      mobs2 = mobs1 + 1
      msrc1 = msrc + 1
      msrc2 = msrc1 + 1
      jmxf = lmxf - lmnf + 1

      ALLOCATE ( zxz(mobs1+5), zxxz1(mobs1+5), zxxz2(mobs1+5) )
      ALLOCATE ( zclft(mobs+5,jmxf), zslft(mobs+5,jmxf) )
      ALLOCATE ( zwkc(mobs+5), zwks(mobs+5),
     $     zdwkc(mobs+5),zdwks(mobs+5) )
      ALLOCATE ( zwklc(mobs+5,jmxf), zwkls(mobs+5,jmxf),
     $     zrspc(mobs+5,jmxf), zrsps(mobs+5,jmxf) )

      zxz(1:mobs1) = SQRT ( zxp(1:mobs1)**2 + zzp(1:mobs1)**2 )
      zxxz1(1:mobs1) = zx(1:mobs1) / zxz(1:mobs1)
      zxxz2(1:mobs1) = 1.0 / ( zx(1:mobs1)*zxz(1:mobs1) )

      CALL pospl1 ( zind,zzp,  1,mobs1, 1,"zzp",   mobs1 )
      CALL pospl1 ( zind,zxz,  1,mobs1, 2,"zxz",   mobs1 )
      CALL pospl1 ( zind,zxxz1,1,mobs1, 3,"zxxz1", mobs1 )
      CALL pospl1 ( zind,zxxz2,1,mobs1, 4,"zxxz2", mobs1 )
      CALL framep(jobid, ff)
      
c..   Store Sine and Cosine of source variable:
      DO l1 = 1, jmxf
         zclft(1:msrc2,l1) = COS ( (lmnf-1+l1)*zind(1:msrc2) )
         zslft(1:msrc2,l1) = SIN ( (lmnf-1+l1)*zind(1:msrc2) )
      END DO

c..   Fourier Analyze the response matrix  source points:
      zwklc(1:mobs,1:jmxf) = ztpii *
     $     MATMUL ( crsp(1:mobs,1:msrc),zclft(1:msrc,1:jmxf) )
      zwkls(1:mobs,1:jmxf) = ztpii * 
     $     MATMUL ( crsp(1:mobs,1:msrc),zslft(1:msrc,1:jmxf) )

c..   Perform Surface Laplacian on observer for each
c     Cos and Sine Fourier component.      
      DO i = 1, jmxf
         zwkc(1:mobs) = zwklc(1:mobs,i)
         zwkc(mobs1) = zwkc(1)
         CALL difspl ( mobs, zind, zwkc, zdwkc )
         zdwkc(mobs1) = zdwkc(1)
         zwks(1:mobs) = zwkls(1:mobs,i)
         zwks(mobs1) = zwks(1)
         CALL difspl ( mobs, zind, zwks, zdwks )
         zdwks(mobs1) = zdwks(1)
         IF ( i == 1 .or. i == jmxf/2 + 1 ) THEN
            CALL pospl2 ( zind,zwkc,zwks,  1,mobs1, 1,"crsp",   i )
            CALL pospl2 ( zind,zdwkc,zdwks,  1,mobs1, 2,"dcrsp", i )
         END IF
         zwkc(1:mobs1) = zxxz1(1:mobs1) * zdwkc(1:mobs1)
         zwks(1:mobs1) = zxxz1(1:mobs1) * zdwks(1:mobs1)
         CALL difspl ( mobs, zind, zwkc, zdwkc )
         CALL difspl ( mobs, zind, zwks, zdwks )
         zdwkc(mobs1) = zdwkc(1)
         zdwks(mobs1) = zdwks(1)
         IF ( i == 1 .or. i == jmxf/2 + 1 ) THEN
            CALL pospl2 ( zind,zdwkc,zdwks, 1,mobs1, 3,"2nd-der", i )
         END IF
         zrspc(1:mobs,i) = zxxz2(1:mobs) * zdwkc(1:mobs)
     $        - (zn/zx(1:mobs))**2 * zwklc(1:mobs,i)
         zrspc(mobs1,i) = zrspc(1,i)
         zrsps(1:mobs,i) = zxxz2(1:mobs) * zdwks(1:mobs)
     $        - (zn/zx(1:mobs))**2 * zwkls(1:mobs,i)
         zrsps(mobs1,i) = zrsps(1,i)
         IF ( i == 1 .or. i == jmxf/2 + 1 ) THEN
            zdwkc(1:mobs1) = zrspc(1:mobs1,i)
            zdwks(1:mobs1) = zrsps(1:mobs1,i)
            CALL pospl2 ( zind,zdwkc,zdwks, 1,mobs1, 4,"dqdt",    i )
            CALL framep(jobid, ff)
         END IF
      END DO

      DEALLOCATE ( zdwkc,zdwks )
      DEALLOCATE ( zxz, zxxz1, zxxz2, zwklc,zwkls )
      ALLOCATE ( zwklc(mobs+3,msrc+3),zwkls(mobs+3,msrc+3) )
      
c.. Reconstruct response in source theta space:
      zwklc(1:mobs,1:msrc) = MATMUL ( zrspc(1:mobs,1:jmxf),
     $     TRANSPOSE ( zclft(1:msrc,1:jmxf) ) ) + 
     $     MATMUL ( zrsps(1:mobs,1:jmxf),
     $     TRANSPOSE ( zslft(1:msrc,1:jmxf) ) )
      zwkls(1:mobs,1:msrc) = MATMUL ( zrspc(1:mobs,1:jmxf),
     $     TRANSPOSE ( zslft(1:msrc,1:jmxf) ) ) -
     $     MATMUL ( zrsps(1:mobs,1:jmxf),
     $     TRANSPOSE ( zclft(1:msrc,1:jmxf) ) )

      DEALLOCATE ( zrspc,zrsps, zclft,zslft )

      OPEN ( UNIT=76, FILE='dqdt_out',
     $     STATUS='REPLACE', FORM='FORMATTED' )
      
      WRITE ( 76, '(/,20x,"DQDT response",/)' )
      WRITE ( 76, '(5x, i5, a )' ) mobs,
     $     tab//tab//"mobs"//tab//"No. of observer points"
      WRITE ( 76, '(5x, i5, a )' ) msrc,
     $     tab//tab//"msrc"//tab//"No. of source points"
      nzn = zn + 0.01
      WRITE ( 76, '(5x, i5, a )' ) nzn,
     $     tab//tab//"n"//tab//"Toroidal mode number"
      WRITE ( 76, '(5x, i5, a )' ) lmnf,
     $     tab//tab//"n"//tab//"lmnf - Resolution"
      WRITE ( 76, '(5x, i5, a )' ) lmxf,
     $     tab//tab//"n"//tab//"lmxf - Resolution"

      WRITE ( 76, '(5x, /a/ )' )  "REAL_Reponse Matrix(obs,srce):"
      DO i = 1, mobs
         WRITE ( 76, '(10es15.7)' ) ( zwklc(i,j), j=1,msrc )
      END DO 

      WRITE ( 76, '(5x, /a/ )' )  "IMAG_Reponse Matrix(obs,srce):"
      DO i = 1, mobs
         WRITE ( 76, '(10es15.7)' ) ( zwkls(i,j), j=1,msrc )
      END DO 

      CALL matwrtn ( zwklc,mobs,msrc,1,1,mobs,msrc,16,8,
     $     "REAL_DQDTW-response matrix", nout1, 0 )
      CALL matwrts ( zwklc,mobs,msrc, 1,1, 25,9, 25,9,
     $     1,1,"REAL_DQDT-response matrix", nout1,0 )
      CALL matwrtn ( zwkls,mobs,msrc,1,1,mobs,msrc,16,8,
     $     "IMAG_DQDTW-response matrix", nout1, 0 )
      CALL matwrts ( zwkls,mobs,msrc, 1,1, 25,9, 25,9,
     $     1,1,"IMAG_DQDT-response matrix", nout1,0 )

      DO izpl = 1, 4
         izplz = (izpl-1)*(msrc/4) + 1
         zwkc(1:mobs) = zwklc(1:mobs,izplz)
         zwks(1:mobs) = zwkls(1:mobs,izplz)
         zwkc(mobs1) = zwkc(1)
         zwks(mobs1) = zwks(1)
         CALL pospl2 ( zind,zwkc,zwks, 1,mobs1, izpl,
     $        "Resp_matrix", izplz )
      END DO
      CALL framep(jobid, ff)
      CLOSE ( UNIT = 76 )
      DEALLOCATE ( zwkc,zwks )
      DEALLOCATE ( zwklc,zwkls )
      
      END SUBROUTINE dqdtwl

c......................................................................
      REAL FUNCTION bw_test ( t, m1, m2 )
c......................................................................

      IMPLICIT NONE
      REAL a, b, c, t
      INTEGER m1, m2
      
      a = 0.0
      b = 0.0
      c = 1.0
      bw_test = a + b * SIN(m1*t) + c * COS(m2*t)

      END FUNCTION bw_test


