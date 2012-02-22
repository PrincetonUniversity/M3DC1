c............................................................
      SUBROUTINE sensorp ( slx, slxt, slz, slzt, nzspt, lpols,
     $     blr, bli, flxr,flxi, pflxr,pflxi )
c............................................................

c.... The magnetic perturbations and the fluxes at the sensor
c      due to the plasma contrbution to the perturbaton.

      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum5.inc'
      INCLUDE 'vacuum7.inc'
c     
      DIMENSION xpla(nths), zpla(nths)

      DIMENSION grdgre(nths2,nths2),grwp(nths,nths),
     .     grri(nths2,nfm2)

      REAL, DIMENSION(*) :: slx, slxt, slz, slzt
      REAL, DIMENSION(*) :: blr, bli
c     
c      REAL, DIMENSION(:), ALLOCATABLE :: slx, slxt, slz ,slzt
      REAL, DIMENSION(:), ALLOCATABLE :: bxr, bxi, bzr, bzi
      REAL, DIMENSION(:,:), ALLOCATABLE :: chir, chii, cwrkr,cwrki
      REAL, DIMENSION(:,:), ALLOCATABLE :: zchc, zchs

      COMMON / bigv1 / grdgre, grwp, grri
c     
      COMMON / gridc / xobs(ndimlpc),zobs(ndimlpc)
c     
      EQUIVALENCE (xpla,xinf), (zpla,zinf)

      CHARACTER(130), DIMENSION(10) :: string

c     nzspt is the number of observer points.

      jmax1 = lmax(1) - lmin(1) + 1
      mth12 = 2*mth
      IF ( lfarw .gt. 0 ) mth12 = mth
c     
c$$$      ALLOCATE ( slx(nzspt), slxt(nzspt), slz(nzspt), slzt(nzspt) )
      ALLOCATE ( chir(4,nzspt),chii(4,nzspt),
     $     cwrkr(4,nzspt),cwrki(4,nzspt) )
c$$$
c$$$      CALL sloops ( slx, slxt, slz ,slzt, nzspt)
c     
      ns = mth
      delx = plrad * deloop
      delz = plrad * deloop
c     
c.....get array of observer points, and initialize chi matrix.
c     
      ALLOCATE ( zchc(mth+1,jmax1), zchs(mth+1,jmax1) )
c     
      DO nsew = 1, 4
c     
         DO i = 1, nzspt
c     
            chir(nsew,i) = 0.0
            chii(nsew,i) = 0.0
            cwrkr(nsew,i) = 0.0
            cwrki(nsew,i) = 0.0
c     
            GO TO ( 51, 52, 53, 54), nsew
c     
 51         CONTINUE
c......north
c     
            xobs(i) = slx(i)
            zobs(i) = slz(i) + delz
            GO TO 90
c     
 52         CONTINUE
c.....south
c     
            xobs(i) = slx(i)
            zobs(i) = slz(i) - delz
            GO TO 90
c     
 53         CONTINUE
c.....east
c     
            xobs(i) = slx(i) + delx
            zobs(i) = slz(i)
            GO TO 90
c     
 54         CONTINUE
c.....west
c     
            xobs(i) = slx(i) - delx
            zobs(i) = slz(i)
c     
 90         CONTINUE

         END DO                 ! Loop on I, nzspt
c     
c.....pass chi on source surface.
c     
c     Integrate over plasma first.
c     

         DO l1 = 1, jmax1
            DO i = 1, mth
               zchc(i,l1) = grri(i,l1)
               zchs(i,l1) = grri(i,jmax1+l1)
            END DO
         END DO
c     
         isg = 1
         WRITE ( outmod,'("isgn for plasma source contrib. = ", i3)')
     $        isg
         WRITE ( iotty,'("isgn for plasma source contrib. = ", i3)')
     $        isg
c     
         CALL chip1 ( xpla,zpla,xplap,zplap,isg, zchc,zchs,jmax1,nzspt,
     $        ns,1,cwrkr,cwrki,nsew, blr,bli )
         CALL atpoint("After Chip1","lfarw",lfarw,"delx",delx,
     $        iotty, outmod )

c..   Store away the plasma result.
c     
         DO  i = 1, nzspt
            chir(nsew,i) = cwrkr(nsew,i)
            chii(nsew,i) = cwrki(nsew,i)
         END DO
c     
c     
         write (outmod,'("Plasma contribution to probes",i4)')nsew
         DO  i = nzspt/2, nzspt/2+2
            write(outmod,'("i,chir,chii,cwrkr,cwrki",i4,1p4e12.3)')
     $           i,chir(nsew,i),chii(nsew,i),cwrkr(nsew,i),cwrki(nsew,i)
         END DO
c     
c     
         IF ( lfarw .gt. 0 ) GO TO 200
c     
c.... Integrate over Shell:
c     
         DO l1 = 1, jmax1
            DO i = 1, mw
               zchc(i,l1) = grri(mth+i,l1)
               zchs(i,l1) = grri(mth+i,jmax1+l1)
            END DO
         END DO
C     
         DO i = 1, nzspt
            cwrkr(nsew,i) = 0.0
            cwrki(nsew,i) = 0.0
         END DO
c     
         isg = -1
         WRITE ( outmod,
     $        '("isgn for shell integ. contrib. = ", i3)') isg
c     
         CALL chip1 ( xwal,zwal,xwalp,zwalp,isg,zchc,zchs,jmax1,nzspt,
     $        ns,0,cwrkr,cwrki,nsew, blr,bli )
c     
         DO  i = 1, nzspt
            chir(nsew,i) = chir(nsew,i) + cwrkr(nsew,i)
            chii(nsew,i) = chii(nsew,i) + cwrki(nsew,i)
         END DO
C     
c     
         write (outmod,'("Plasma contribution to probes",i4)') nsew
         DO  i = nzspt/2, nzspt/2+2
            write(outmod,'("i,chir,chii,cwrkr,cwrki",i4,1p4e12.3)')
     $           i,chir(nsew,i),chii(nsew,i),cwrkr(nsew,i),cwrki(nsew,i)
         END DO
c     
 200     CONTINUE

      END DO                    ! Loop over NSEW
c     
c$$$  WRITE ( OUTMOD, '(/,36x,"****   Mirnov Loops   ****")' )
c$$$  WRITE ( OUTMOD,
c$$$  $     '(1x,"iobs",18x,"chir(nsew)",35x,"chii(nsew)" )' )
c$$$  DO i = 1, nzspt
c$$$  write (outmod,'(1x, i4, 1p4e11.3, 2x, 1p4e11.3)') i,
c$$$  $        (chir(ii,i),ii=1,4), (chii(ii,i), ii= 1, 4)
c$$$  END DO
c     
c.....calculate bx and bz
c     
      ALLOCATE ( bxr(nzspt), bxi(nzspt), bzr(nzspt), bzi(nzspt) )
c     
      DO i = 1, nzspt
         bxr(i) = 0.0
         bxi(i) = 0.0
         bzr(i) = 0.0
         bzi(i) = 0.0
      END DO

      DO i = 1, nzspt
c     
         bxr(i) = ( chir(3,i) - chir(4,i) ) / (2.0*delx)
         bxi(i) = ( chii(3,i) - chii(4,i) ) / (2.0*delx)
         bzr(i) = ( chir(1,i) - chir(2,i) ) / (2.0*delz)
         bzi(i) = ( chii(1,i) - chii(2,i) ) / (2.0*delz)
         WRITE ( outmod,'("bxr,ir,zr,zi = ", i3, 1p4e12.3)')
     $        i, bxr(i),bxi(i),bzr(i),bzi(i)
      END DO

         CALL rflux ( slx,slxt, slzt, bxr,bzr, nzspt, dtsl, flxr )
         CALL rflux ( slx,slxt, slzt, bxi,bzi, nzspt, dtsl, flxi )
         CALL rfluxp ( slx,slxt, slzt, bxr,bzr, nzspt, dtsl, pflxr )
         CALL rfluxp ( slx,slxt, slzt, bxi,bzi, nzspt, dtsl, pflxi )

c$$$      DEALLOCATE ( slx, slxt, slz, slzt )
      DEALLOCATE ( bxr, bxi, bzr, bzi )
      DEALLOCATE ( chir,chii, cwrkr,cwrki )
      DEALLOCATE ( zchc, zchs )

      WRITE ( 6, '( "XXXXX  SENSORP nzspt = ", i5 )' ) nzspt
      WRITE ( 23, '( "XXXXX SENSORP nzspt = ", i5 )' ) nzspt

 500  CONTINUE

      END SUBROUTINE sensorp
c     
c............................................................
      SUBROUTINE sensorw ( slx, slxt, slz, slzt, nzspt, lpols,
     $     ble, blo, flxr,flxi, pflxr,pflxi )
c............................................................

c.... The magnetic perturbations and the fluxes at the sensor
c      due to the wall contrbution to the perturbaton.

c... BLE(jwal2) and BLO(jwal2) contains the (EC,OC) and (ES,OS) 
c     perturbations, respecitvely

      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum5.inc'
      INCLUDE 'vacuum7.inc'
c     
      DIMENSION xpla(nths), zpla(nths)

      REAL, DIMENSION(nths2,nths2):: grdgre
      REAL, DIMENSION(nths ,nths ):: grwp
      REAL, DIMENSION(nths2,nfm2 ):: grri

      REAL, DIMENSION(*) :: slx, slxt, slz, slzt
      REAL, DIMENSION(*) :: ble, blo
c     
c$$$      REAL, DIMENSION(:), ALLOCATABLE :: slx, slxt, slz ,slzt
      REAL, DIMENSION(:), ALLOCATABLE :: bxr, bxi, bzr, bzi
      REAL, DIMENSION(:,:), ALLOCATABLE :: chir, chii, cwrkr,cwrki
      REAL, DIMENSION(:,:), ALLOCATABLE :: zchc, zchs

      COMMON / bigv1 / grdgre, grwp, grri
c     
      COMMON / gridc / xobs(ndimlpc),zobs(ndimlpc)
c     
      EQUIVALENCE (xpla,xinf), (zpla,zinf)

      CHARACTER(130), DIMENSION(10) :: string

c     nzspt is the number of observer points.

      jmax1 = lmax(1) - lmin(1) + 1
      mth12 = 2*mth
      IF ( lfarw .gt. 0 ) mth12 = mth
c     
c$$$      ALLOCATE ( slx(nzspt), slxt(nzspt), slz(nzspt) ,slzt(nzspt))
      ALLOCATE ( chir(4,nzspt),chii(4,nzspt),
     $     cwrkr(4,nzspt),cwrki(4,nzspt) )
c$$$
c$$$      CALL sloops ( slx, slxt, slz ,slzt, nzspt)
c     
      ns = mth
      delx = plrad * deloop
      delz = plrad * deloop
c     
c.....get array of observer points, and initialize chi matrix.
c     
      ALLOCATE ( zchc(mth+1,jwal1), zchs(mth+1,jwal1) )
c     
      DO nsew = 1, 4
c     
         DO i = 1, nzspt
c     
            chir(nsew,i) = 0.0
            chii(nsew,i) = 0.0
            cwrkr(nsew,i) = 0.0
            cwrki(nsew,i) = 0.0
c     
            GO TO ( 51, 52, 53, 54), nsew
c     
 51         CONTINUE
c......north
c     
            xobs(i) = slx(i)
            zobs(i) = slz(i) + delz
            GO TO 90
c     
 52         CONTINUE
c.....south
c     
            xobs(i) = slx(i)
            zobs(i) = slz(i) - delz
            GO TO 90
c     
 53         CONTINUE
c.....east
c     
            xobs(i) = slx(i) + delx
            zobs(i) = slz(i)
            GO TO 90
c     
 54         CONTINUE
c.....west
c     
            xobs(i) = slx(i) - delx
            zobs(i) = slz(i)
c     
 90         CONTINUE

         END DO                 ! loop on i, nzspt
c     
c.....pass chi on source surface.
c     
c     Integrate Over Shell First.
c     
c     DO i=1,mw
c     DO l1 = 1, jwal1
c     zchc(1:mw,l1) = gwin(mth+1:mth+mw,l1)
c     zchs(1:mw,l1) = gwin(mth+1:mth+mw,jwal1+l1)
         zchc(1:mw,1:jwal1) =
     $        grri(mth+1:mth+mw,2*jmax1+1:2*jmax1+jwal1)
         zchs(1:mw,1:jwal1) =
     $        grri(mth+1:mth+mw,2*jmax1+jwal1+1:2*jmax1+jwal2)
c     END DO
c     END DO
c     
         isg = 1
         WRITE ( outmod,'("isgn for shell source contrib. = ", i3)')
     $        isg
         WRITE ( iotty,'("isgn for shell source contrib. = ", i3)')
     $        isg
c     
         CALL chis1 ( xwal,zwal,xwalp,zwalp,isg, zchc,zchs,jwal1,nzspt,
     $        ns,1,cwrkr,cwrki,nsew, ble,blo )
c     $        ns,0,cwrkr,cwrki,nsew, ble,blo )
         CALL atpoint("After Chis1","lfarw",lfarw,"delx",delx,
     $        iotty, outmod )

c..   Store Away the Shell Result.
c     
         DO  i = 1, nzspt
            chir(nsew,i) = cwrkr(nsew,i)
            chii(nsew,i) = cwrki(nsew,i)
         END DO
c     
c     
         write (outmod,'("Wall contribution to probes, nsew",i4)')nsew
         DO  i = nzspt/2, nzspt/2+2
            write(outmod,'("i,chir,chii,cwrkr,cwrki",i4,1p4e12.3)')
     $           i,chir(nsew,i),chii(nsew,i),cwrkr(nsew,i),cwrki(nsew,i)
         END DO
c     
c     
c.... Integrate over Plasma:
c     
c     DO i=1,mth
c     DO l1 = 1, jwal1
c     zchc(1:mth,l1) = gwin(1:mth,l1)
c     zchs(1:mth,l1) = gwin(1:mth,jwal1+l1)
         zchc(1:mth,1:jwal1) =
     $        grri(1:mth,2*jmax1+1:2*jmax1+jwal1)
         zchs(1:mth,1:jwal1) =
     $        grri(1:mth,2*jmax1+jwal1+1:2*jmax1+jwal2)
c     END DO
c     END DO
C     
         DO i = 1, nzspt
            cwrkr(nsew,i) = 0.0
            cwrki(nsew,i) = 0.0
         END DO
c     
         isg = -1
         WRITE ( outmod,
     $        '("isgn for shell integ. contrib. = ", i3)') isg
c     
         CALL chis1 ( xpla,zpla,xplap,zplap,isg,zchc,zchs,jwal1,nzspt,
     $        ns,0,cwrkr,cwrki,nsew, ble,blo )
c     
         DO  i = 1, nzspt
            chir(nsew,i) = chir(nsew,i) + cwrkr(nsew,i)
            chii(nsew,i) = chii(nsew,i) + cwrki(nsew,i)
         END DO
c     
         write (outmod,'("Wall contribution to probes, nsew",i4)')nsew
         DO  i = nzspt/2, nzspt/2+2
            write(outmod,'("i,chir,chii,cwrkr,cwrki",i4,1p4e12.3)')
     $           i,chir(nsew,i),chii(nsew,i),cwrkr(nsew,i),cwrki(nsew,i)
         END DO
c     
      END DO                    ! Loop over NSEW, 4
c     
c$$$  WRITE ( OUTMOD, '(/,36x,"****   Mirnov Loops   ****")' )
c$$$  WRITE ( OUTMOD,
c$$$  $     '(1x,"iobs",18x,"chir(nsew)",35x,"chii(nsew)" )' )
c$$$  DO i = 1, nzspt
c$$$  write (outmod,'(1x, i4, 1p4e11.3, 2x, 1p4e11.3)') i,
c$$$  $        (chir(ii,i),ii=1,4), (chii(ii,i), ii= 1, 4)
c$$$  END DO
c     
c.....calculate bx and bz
c     
      ALLOCATE ( bxr(nzspt), bxi(nzspt), bzr(nzspt), bzi(nzspt) )
c     
      DO i = 1, nzspt
         bxr(i) = 0.0
         bxi(i) = 0.0
         bzr(i) = 0.0
         bzi(i) = 0.0
      END DO

      DO i = 1, nzspt
c     
         bxr(i) = ( chir(3,i) - chir(4,i) ) / (2.0*delx)
         bxi(i) = ( chii(3,i) - chii(4,i) ) / (2.0*delx)
         bzr(i) = ( chir(1,i) - chir(2,i) ) / (2.0*delz)
         bzi(i) = ( chii(1,i) - chii(2,i) ) / (2.0*delz)
         WRITE ( outmod,'("bxr,ir,zr,zi = ", i3, 1p4e12.3)')
     $        i, bxr(i),bxi(i),bzr(i),bzi(i)
      END DO

         CALL rflux ( slx,slxt, slzt, bxr,bzr, nzspt, dtsl, flxr )
         CALL rflux ( slx,slxt, slzt, bxi,bzi, nzspt, dtsl, flxi )
         CALL rfluxp ( slx,slxt, slzt, bxr,bzr, nzspt, dtsl, pflxr )
         CALL rfluxp ( slx,slxt, slzt, bxi,bzi, nzspt, dtsl, pflxi )

c$$$      DEALLOCATE ( slx, slxt, slz, slzt )
      DEALLOCATE ( bxr, bxi, bzr, bzi )
      DEALLOCATE ( chir,chii, cwrkr,cwrki )
      DEALLOCATE ( zchc, zchs )

      WRITE ( 6, '( "XXXXX  SENSORW nzspt = ", i5 )' ) nzspt
      WRITE ( 23, '( "XXXXX SENSORW nzspt = ", i5 )' ) nzspt

 500  CONTINUE
      RETURN
      END SUBROUTINE sensorw
c     
c......................................................................
      SUBROUTINE chip1(xsce,zsce,xscp,zscp,isg,creal,cimag,ndimc,nobs,
     $     ns,ip,chir,chii,nsew,blr,bli)
c......................................................................
c     
      include 'vacuum1.inc'
      include 'vacuum3.inc'
      include 'vacuum2.inc'
      include 'vacuum5.inc'
c     
      dimension iop(2), ww1(nths),ww2(nths),ww3(nths),tab(3)
      dimension blr(*),bli(*), xsce(*),zsce(*),xscp(*),zscp(*)
      dimension creal(ns+1,*), cimag(ns+1,*)
      dimension chir(4,*), chii(4,*)
c     
      real nq
c     
      common / gridc / xobs(ndimlpc),zobs(ndimlpc)
c     
c     Factpi for dividing BVAL to be consistent with AVAL.
c     
      factpi = twopi
      jmax1 = lmax(1) - lmin(1) + 1
      q = qa1
c     
      nq = n * q
      dtpw = twopi / ns
c     
      ns1 = ns + 1
c     nobs = nloop + 3*nloopr
c     
      do io = 1, nobs
c     
         xs = xobs(io)
         zs = zobs(io)
c     
         do is = 1, ns

            xt = xsce(is)
            zt = zsce(is)
            xtp = xscp(is)
            ztp = zscp(is)
c     
            call green

c     Divide BVAL by twopi to account for definition of script G to be
c     onsistent with K.
c     
            bval = bval / factpi
c     
            do l1 = 1, jmax1
c     
               zbr = blr(l1)
               zbi = bli(l1)
               chir(nsew,io) = chir(nsew,io) +
     $              aval * ( creal(is,l1)*zbr - cimag(is,l1)*zbi )
               chii(nsew,io) = chii(nsew,io) +
     $              aval * ( cimag(is,l1)*zbr + creal(is,l1)*zbi )
c     
               if ( ip .eq. 0 ) go to 60
c     
c##########change plus to minus bval  #########
               ipm = + 1.0
               chir(nsew,io) = chir(nsew,io) + ipm * bval
     $              * ( cslth(is,l1)*zbr - snlth(is,l1)*zbi )
               chii(nsew,io) = chii(nsew,io) + ipm * bval
     $              * ( snlth(is,l1)*zbr + cslth(is,l1)*zbi )
c     
 60            continue
            end do              ! Loop over l1: 1, jmax1
         end do                 ! Loop over IS: 1, ns, sources
c     
         chir(nsew,io) = 0.5 * isg*dtpw * chir(nsew,io)
         chii(nsew,io) = 0.5 * isg*dtpw * chii(nsew,io)
c     
      end do                    ! Loop over Io: 1, nobs, observers
c     
      RETURN
      END SUBROUTINE chip1
c     
c......................................................................
      SUBROUTINE chis1(xsce,zsce,xscp,zscp,isg,creal,cimag,ndimc,nobs,
     $     ns,ip,chir,chii,nsew,blr,bli)
c......................................................................

c... BLR(jwal2) and BLI(jwal2) contains the (EC,OC) and (ES,OS) 
c     perturbations, respecitvely

c     
      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum5.inc'
      INCLUDE 'vacuum7.inc'
c     
      dimension iop(2), ww1(nths),ww2(nths),ww3(nths),tab(3)
      dimension blr(*),bli(*), xsce(*),zsce(*),xscp(*),zscp(*)
      dimension creal(ns+1,*), cimag(ns+1,*)
      dimension chir(4,*), chii(4,*)
c     
      real nq
c     
      common / gridc / xobs(ndimlpc),zobs(ndimlpc)
c     
c     Factpi for dividing BVAL to be consistent with AVAL.
c     
      factpi = twopi
      q = qa1
c     
      nq = n * q
      dtpw = twopi / ns
c     
      ns1 = ns + 1
c     nobs = nloop + 3*nloopr
c     
      do io = 1, nobs
c     
         xs = xobs(io)
         zs = zobs(io)
c     
         do is = 1, ns

            xt = xsce(is)
            zt = zsce(is)
            xtp = xscp(is)
            ztp = zscp(is)
c     
            call green

c     Divide BVAL by twopi to account for definition of script G to be
c     onsistent with K.
c     
            bval = bval / factpi
c     
            do l1 = 1, jwal1
c     
               zbr = blr(l1)
               zbi = bli(l1)
               zboc = blr(l1+jwal1)
               zbos = bli(l1+jwal1)

c     chir(nsew,io) = chir(nsew,io) +
               chir(nsew,io) = chir(nsew,io) -
c     $              aval * ( creal(is,l1)*zbr + cimag(is,l1)*zbi )/2
c     $              aval * ( creal(is,l1)*zbr - cimag(is,l1)*zbi )/2
c     $              aval * ( creal(is,l1) * zbr ) / 2.0
     $              aval * ( creal(is,l1)*zbr + cimag(is,l1)*zboc )/2.0
c     chii(nsew,io) = chii(nsew,io) +
               chii(nsew,io) = chii(nsew,io) -
c     $              aval * ( cimag(is,l1)*zbr - creal(is,l1)*zbi )/2
c     $              aval * ( cimag(is,l1)*zbr + creal(is,l1)*zbi )/2
c     $             aval * ( cimag(is,l1) * zbi ) / 2.0
     $              aval * ( cimag(is,l1)*zbos + creal(is,l1)*zbi )/2.0
c     
               if ( ip .eq. 0. ) go to 60
c     
c##########change plus to minus bval  #########
               ipm = + 1.0
c     it turned out that a should change sign on the wall
c     ipm= - 1.0
c$$$               chir(nsew,io) = chir(nsew,io) + ipm * bval
c$$$     $              * rwvce(is,l1)*zbr /2.
c$$$               chii(nsew,io) = chii(nsew,io) + ipm * bval
c$$$     $              * rwvco(is,l1)*zbi /2.
               chir(nsew,io) = chir(nsew,io) + ipm * bval
     $              * ( rwvce(is,l1)*zbr + rwvco(is,l1)*zboc )/2.0
               chii(nsew,io) = chii(nsew,io) + ipm * bval
     $              * ( rwvco(is,l1)*zbos + rwvce(is,l1)*zbi )/2.0

c$$$  chir(nsew,io) = chir(nsew,io) + ipm * bval
c$$$  $              * ( cslth(is,l1)*zbr - snlth(is,l1)*zbi )
c$$$  chii(nsew,io) = chii(nsew,io) + ipm * bval
c$$$  $              * ( snlth(is,l1)*zbr + cslth(is,l1)*zbi )
c     
 60            continue
            end do              ! Loop over l1: 1, jwal1
         end do                 ! Loop over IS: 1, ns, sources
c     
         chir(nsew,io) = 0.5 * isg*dtpw * chir(nsew,io)
         chii(nsew,io) = 0.5 * isg*dtpw * chii(nsew,io)
c     
      end do                    ! Loop over Io: 1, nobs, observers
c     
      RETURN
      END SUBROUTINE chis1

c............................................
      SUBROUTINE diii_loop_meas ( nzspt, zdtsl, lpols )
c............................................

c.. DIII-D Diagnostic loops and Sensor loops. Only radial loops, so far.

c.. Storage:
c...  (x-center, z-center, phi location, axis angle wrt x-axis, 
c      length of coil, width in phi, NA:number_of_loops area_cm^2)

c... Actual DIII-D loops:
      REAL, DIMENSION(:,:), ALLOCATABLE :: r_udb_i, r_r0_i, r_rp1_i,
     $     r_rm1_i, r_67_o, r_157_o, r_rp1_o, r_r0_o, r_rm1_o,
     $     R_r0_e, r_rp1_e, r_rm1_e
      REAL, DIMENSION(:,:), ALLOCATABLE:: p_r0_ta

c... Others:
      REAL, DIMENSION(:,:), ALLOCATABLE ::  horz_sens, vert_sens_01

      REAL, DIMENSION(:), ALLOCATABLE :: slx, slxt, slz, slzt

      ALLOCATE ( r_udb_i(6,7), r_r0_i(6,7), r_rp1_i(6,7),
     $     r_rm1_i(6,7), r_67_o(11,7), r_157_o(11,7), r_rp1_o(4,7),
     $     r_r0_o(6,7), r_rm1_o(4,7), r_r0_e(6,7), r_rp1_e(12,7),
     $     r_rm1_e(12,7) )

      ALLOCATE ( p_r0_ta(11,7) )

      ALLOCATE ( horz_sens(5,7), vert_sens_01(3,7) )

C... Radial sensors:

      r_udb_i = TRANSPOSE ( RESHAPE ( (/   
     $     1.506, 1.270,  185.0,  49.5, 0.094,  60.1,   1391.2,
     $     1.590, 1.198,  185.0,  49.5, 0.088,  60.2,   1484.3,
     $     1.730, 1.142,  185.2,  92.3, 0.106,  59.2,   1890.9,
     $     1.897, 1.133,  157.6,  90.5, 0.165,   6.2,    642.2, 
     $     1.185, 1.286,  168.7, 156.6, 0.091, 108.0,   2046.6,
     $     1.086, 1.228,  168.7,  90.4, 0.107, 107.6,   2178.2 /), 
     $     (/ 7, 6 /) ) )

      r_r0_i =  TRANSPOSE ( RESHAPE ( (/   
     $     2.431, 0.000,  17.2,  0.0, 0.800,  60.3,   20466.3,
     $     2.431, 0.000,  72.4,  0.0, 0.800,  50.1,   17006.6,
     $     2.431, 0.000, 132.6,  0.0, 0.800,  70.3,   23849.9,
     $     2.431, 0.000, 197.6,  0.0, 0.800,  59.6,   20237.8,
     $     2.431, 0.000, 252.4,  0.0, 0.800,  50.1,   17006.8,
     $     2.431, 0.000, 312.3,  0.0, 0.800,  69.6,   23632.8  /),
     $     (/ 7, 6 /) ) )

      r_rp1_i =  TRANSPOSE ( RESHAPE ( (/    
     $     2.300,  0.714,  17.8,  22.6, 0.680,  61.5, 16288.6,
     $     2.300,  0.714,  73.0,  22.6, 0.680,  48.9, 13344.2,
     $     2.300,  0.714, 133.0,  22.6, 0.680,  71.1, 19396.3,
     $     2.300,  0.714, 198.6,  22.6, 0.680,  60.0, 16370.3,
     $     2.300,  0.714, 253.0,  22.6, 0.680,  48.9, 13344.4,
     $     2.300,  0.714, 312.3,  22.6, 0.680,  69.6, 19485.8  /), 
     $     (/ 7, 6 /) ) )

      r_rm1_i =  TRANSPOSE ( RESHAPE ( (/    
     $     2.300, -0.714,  18.7, -22.6,  0.680, 59.7, 16778.5,
     $     2.300, -0.714,  73.0, -22.6,  0.680, 48.9, 13344.2,
     $     2.300, -0.714, 133.0, -22.6,  0.680, 71.1, 19396.3,
     $     2.300, -0.714, 198.6, -22.6,  0.680, 60.0, 16370.3,
     $     2.300, -0.714, 253.0, -22.6,  0.680, 48.9, 13344.4,
     $     2.300, -0.714, 313.2, -22.6,  0.680, 71.4, 18997.3  /), 
     $     (/ 7, 6 /) ) )

      r_67_o = TRANSPOSE ( RESHAPE ( (/    
     $     0.925,  0.343, 67.0,  180.0, 0.682, 30.0,  3301.1,
     $     0.942,  1.003, 67.0,  176.9, 0.631, 30.0,  3113.0,
     $     1.320,  1.393, 67.0,  101.4, 0.729, 29.7,  4990.6,
     $     1.946,  1.263, 67.0,   51.0, 0.639, 29.6,  6426.2,
     $     2.315,  0.780, 67.0,   22.4, 0.557, 29.7,  6687.5,
     $     2.433,  0.000, 52.0,   -0.1, 1.036, 24.8,  10912.,
     $     2.313, -0.781, 67.0,  -22.4, 0.530, 29.9,  6405.5,
     $     2.025, -1.262, 67.0,  -39.9, 0.519, 29.9,  5488.3,
     $     1.405, -1.392, 67.0,  -99.1, 0.896, 29.9,  6574.8,
     $     0.942, -1.003, 67.0, -176.9, 0.631, 30.0,  3113.0,
     $     0.925, -0.343, 67.0,  180.0, 0.682, 30.0,  3301.1  /), 
     $     (/ 7, 11 /) ) )
      
      r_157_o = TRANSPOSE ( RESHAPE ( (/    
     $     0.925,  0.343, 157.0,  180.0,  0.682, 30.0, 3301.1,
     $     0.942,  1.003, 157.0,  176.9,  0.631, 30.0, 3113.0,
     $     1.320,  1.393, 157.0,  101.4,  0.729, 30.0, 5035.9,
     $     1.944,  1.265, 157.0,   51.0,  0.635, 29.9, 6430.7,
     $     2.314,  0.783, 157.0,   22.4,  0.557, 30.1, 6771.8,
     $     2.433,  0.001, 152.0,    0.0,  1.037, 29.1, 12821.,
     $     2.313, -0.781, 157.0,  -22.4,  0.549, 30.0, 6653.3,
     $     2.025, -1.262, 157.0,  -39.9,  0.523, 30.0, 5550.0,
     $     1.406, -1.392, 157.0,  -99.1,  0.897, 30.0, 6599.5,
     $     0.942, -1.003, 157.0, -176.9,  0.631, 30.0, 3113.0,
     $     0.925, -0.343, 157.0,  180.0,  0.682, 30.0, 3301.1  /), 
     $     (/ 7, 11 /) ) )
      
      r_rp1_o = TRANSPOSE ( RESHAPE ( (/    
     $     2.315,  0.780,  67.0,  22.4,  0.557,  29.7, 6687.5,
     $     2.314,  0.783, 157.0,  22.4,  0.557,  30.1, 6771.8,
     $     2.289,  0.845, 345.0,  22.4,  0.432,  21.0, 3613.3,
     $     2.395,  0.587, 345.0,  22.4,  0.124,  21.0, 1085.8  /), 
     $     (/ 7, 4 /) ) )
      
      r_r0_o = TRANSPOSE ( RESHAPE ( (/    
     $     2.433,  0.000,  52.0, -0.1,  1.036,  24.8,  10912.,
     $     2.449,  0.261, 132.0,  3.5,  0.519,  13.1,  2895.8,
     $     2.447, -0.262, 132.0, -3.8,  0.522,  13.1,  2914.0,
     $     2.433,  0.001, 152.0,  0.0,  1.037,  29.1,  12821.,
     $     2.449,  0.260, 312.0,  3.4,  0.519,  12.8,  2848.1,
     $     2.447, -0.261, 312.0, -3.7,  0.522,  12.8,  2862.7 /) ,
     $     (/ 7, 6 /) ) )

      r_rm1_o =  TRANSPOSE ( RESHAPE ( (/    
     $     2.393, -0.588,  15.0, -22.4, 0.124, 15.0,  776.5, 
     $     2.287, -0.845,  15.0, -22.4, 0.430, 15.0, 2575.1,
     $     2.313, -0.781,  67.0, -22.4, 0.530, 29.9, 6405.5,
     $     2.313, -0.781, 157.0, -22.4, 0.549, 30.0, 6653.3  /), 
     $     (/ 7, 4 /) ) )
      
      r_r0_e = TRANSPOSE ( RESHAPE ( (/    
     $     2.477, 0.000,  19.0, 0.0, 1.194, 60.0, 30959.8,
     $     2.477, 0.000,  79.0, 0.0, 1.194, 60.0, 30959.8,
     $     2.477, 0.000, 139.0, 0.0, 1.194, 60.0, 30959.8,
     $     2.477, 0.000, 199.0, 0.0, 1.194, 60.0, 30959.8,
     $     2.477, 0.000, 259.0, 0.0, 1.194, 60.0, 30959.8,
     $     2.477, 0.000, 319.0, 0.0, 1.194, 60.0, 30959.8  /), 
     $     (/ 7, 6 /) ) )

      r_rp1_e = TRANSPOSE ( RESHAPE ( (/    
     $     2.389, 0.806,   4.0, 24.1, 0.454, 30.0, 22696.9,
     $     2.389, 0.806,  34.0, 24.1, 0.454, 30.0, 22696.9,
     $     2.389, 0.806,  64.0, 24.1, 0.454, 30.0, 22696.9,
     $     2.389, 0.806,  94.0, 24.1, 0.454, 30.0, 22696.9,
     $     2.389, 0.806, 124.0, 24.1, 0.454, 30.0, 22696.9,
     $     2.389, 0.806, 154.0, 24.1, 0.454, 30.0, 22696.9,
     $     2.389, 0.806, 184.0, 24.1, 0.454, 30.0, 22696.9,
     $     2.389, 0.806, 214.0, 24.1, 0.454, 30.0, 22696.9,
     $     2.389, 0.806, 244.0, 24.1, 0.454, 30.0, 22696.9,
     $     2.389, 0.806, 274.0, 24.1, 0.454, 30.0, 22696.9,
     $     2.389, 0.806, 304.0, 24.1, 0.454, 30.0, 22696.9,
     $     2.389, 0.806, 334.0, 24.1, 0.454, 30.0, 22696.9  /), 
     $     (/ 7, 12 /) ) )

      r_rm1_e = TRANSPOSE ( RESHAPE ( (/    
     $     2.389, -0.806,   4.0, -24.1, 0.454, 30.0, 22696.9,
     $     2.389, -0.806,  34.0, -24.1, 0.454, 30.0, 22696.9,
     $     2.389, -0.806,  64.0, -24.1, 0.454, 30.0, 22696.9,
     $     2.389, -0.806,  94.0, -24.1, 0.454, 30.0, 22696.9,
     $     2.389, -0.806, 124.0, -24.1, 0.454, 30.0, 22696.9,
     $     2.389, -0.806, 154.0, -24.1, 0.454, 30.0, 22696.9,
     $     2.389, -0.806, 184.0, -24.1, 0.454, 30.0, 22696.9,
     $     2.389, -0.806, 214.0, -24.1, 0.454, 30.0, 22696.9,
     $     2.389, -0.806, 244.0, -24.1, 0.454, 30.0, 22696.9,
     $     2.389, -0.806, 274.0, -24.1, 0.454, 30.0, 22696.9,
     $     2.389, -0.806, 304.0, -24.1, 0.454, 30.0, 22696.9,
     $     2.389, -0.806, 334.0, -24.1, 0.454, 30.0, 22696.9  /), 
     $     (/ 7, 12 /) ) )

c... Poloidal Sensors. \phi width set to 0.0 for now:
			
      p_r0_ta = TRANSPOSE ( RESHAPE ( (/  
     $	2.413,	 0.003,	 67.5,	-89.9,	0.141,	0.0,	621.8,
     $	2.413,	-0.005,	 97.4,	-89.8,	0.137,	0.0,	608.7,
     $	2.413,	 0.002,	127.9,	-90.0,	0.140,	0.0,	598.4,
     $	2.416,	 0.006,	137.4,	-90.4,	0.054,	0.0,	231.0,
     $	2.413,	-0.001,	157.6,	-89.8,	0.137,	0.0,	610.9,
     $	2.413,	-0.003,	246.4,	-90.5,	0.140,	0.0,	615.0,
     $	2.413,  -0.009,	277.5,	-89.7,	0.137,	0.0,	619.4,
     $	2.413,	 0.001,	307.0,	-90.2,	0.140,	0.0,	598.4,
     $	2.411,	-0.001,	312.4,	-90.0,	0.054,	0.0,	229.7,
     $	2.418,	-0.001,	317.4,	-89.9,	0.140,	0.0,	604.0,
     $	2.413,	-0.002,	339.7,	-89.8,	0.140,	0.0,	607.2  /),
     $  (/ 7, 11 /) ) )
						
c... Ad Hoc sensors for checking:

      horz_sens = TRANSPOSE ( RESHAPE ( (/
     $     1.750, -0.50, 0.0, 90.0, 2.5, 30.0, 100.0,
     $     1.750, -0.25, 0.0, 90.0, 2.5, 30.0, 100.0,
     $     1.750, -0.00, 0.0, 90.0, 2.5, 30.0, 100.0,
     $     1.750,  0.25, 0.0, 90.0, 2.5, 30.0, 100.0,
     $     1.750,  0.50, 0.0, 90.0, 2.5, 30.0, 100.0  /),
     $     (/ 7, 5 /) ) )

      vert_sens_01 = TRANSPOSE ( RESHAPE ( (/
     $     2.376,  0.00, 0.0,  00.0, 0.5240, 30.0, 100.0,
     $     2.388,  0.00, 0.0,  00.0, 0.5240, 30.0, 100.0,
     $     2.400,  0.00, 0.0,  00.0, 0.5240, 30.0, 100.0  /),
     $     (/ 7, 3 /) ) )

      nzspt = 51

      ALLOCATE ( slx(nzspt), slxt(nzspt), slz(nzspt), slzt(nzspt) )

      WRITE ( 28, '( /, 1x,"---------------------------------" )' ) 
      WRITE ( 28, '( 1x,"..... R_R0_I.... nzsp = ", i5 )' ) nzspt

c$$$      CALL r_loops_row ( r_r0_i, 6,7, 1, nzspt,
c$$$     $     slx,slxt, slz,slzt, zphil, zzpol, zdtsl, "r_r0_i" )

c..   zzpol aligns probe's axis so that measurements are along the 
c           correct direction. 
c           zzpol = 0.0 for radial probes. 1.0 for poloidal probes.

c..  The fields and flux in Bpol - R0 toroidal array.

      zzpol = 1.0
      DO izcase = 1, 11
         CALL r_loops_row ( p_r0_ta, 11,7, izcase, nzspt,
     $        slx,slxt, slz,slzt, zdtsl, zphil, zzpol, "Bp_r0_ta" )

         WRITE ( 28, '( /, "slx, slxt, slz, slzt(1), zdtsl = ",
     $        /, 5es13.5 )' ) slx(1),slxt(1),slz(1),slzt(1), zdtsl
         
         CALL sensflxi (slx, slxt, slz, slzt, nzspt, zphil, lpols,
     $        "Bp_r0_ta" )
      END DO

c..  The fields and flux in the vertical sensor 01

      zzpol = 0.0
      DO izcase = 1, 1
         CALL r_loops_row ( vert_sens_01, 3,7, izcase, nzspt,
     $        slx,slxt, slz,slzt, zdtsl, zphil, zzpol, "vert_sens_01" )

         WRITE ( 28, '( /, "slx, slxt, slz, slzt(1), zdtsl = ",
     $        /, 5es13.5 )' ) slx(1),slxt(1),slz(1),slzt(1), zdtsl
         
         CALL sensflxi (slx, slxt, slz, slzt, nzspt, zphil, lpols,
     $        "vert_sens_01" )
      END DO

c..  Check the fields along x at constant z:

      zzpol = 1.0
      DO izcase = 1, 1
         CALL r_loops_row ( horz_sens, 5,7, izcase, nzspt,
     $        slx,slxt, slz,slzt, zdtsl, zphil, zzpol, "horz_sens" )

         WRITE ( 28, '( /, "slx, slxt, slz, slzt(1), zdtsl = ",
     $        /, 5es13.5 )' ) slx(1),slxt(1),slz(1),slzt(1), zdtsl
         
         CALL sensflxi (slx, slxt, slz, slzt, nzspt, zphil, lpols,
     $        "horz_sens" )
         
         WRITE ( 6, '( /, "slx, slxt, slz, slzt(1), zdtsl = ",
     $        /, 5es13.5 )' ) slx(1),slxt(1),slz(1),slzt(1), zdtsl
         WRITE ( 23, '( /, "slx, slxt, slz, slzt(1), zdtsl = ",
     $        /, 5es13.5 )' ) slx(1),slxt(1),slz(1),slzt(1), zdtsl

      END DO

      WRITE ( 6, '( "XXXXX  DIII nzspt = ", i5 )' ) nzspt
      WRITE ( 23, '( "XXXXX  DIII nzspt = ", i5 )' ) nzspt

      DEALLOCATE ( slx, slxt, slz, slzt )

      DEALLOCATE ( r_udb_i, r_r0_i, r_rp1_i,
     $     r_rm1_i, r_67_o, r_157_o, r_rp1_o, r_r0_o, r_rm1_o,
     $     R_r0_e, r_rp1_e, r_rm1_e )

      DEALLOCATE ( p_r0_ta )

      DEALLOCATE ( horz_sens, vert_sens_01 )

      END SUBROUTINE diii_loop_meas

c.......................................................
      SUBROUTINE r_loops_row ( r_in, nrow,ncol, irow, nzspt,
     $     slx,slxt, slz,slzt, zdtsl, zphil, zpol, s_name )
c......................................................
c.
c... x, z, and their derivatives on a grid over the loop face
c          at any angle, zalfa.
c    nzspt grid points in theta (dimension length here).
c.. lzpol = 1 aligns the poloidal sensors so that dl is along the 
c             axis of the probe. lzpol = 0 for radial sensors.

      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum5.inc'
      INCLUDE 'vacuum8.inc'
      INCLUDE 'vacuum7.inc'
c     

      REAL, DIMENSION(*) :: slx, slxt, slz, slzt
      CHARACTER*(*) s_name
      DIMENSION r_in(nrow,ncol)

      zpi = ACOS ( -1.0 )

      zx0 =   r_in ( irow, 1 ) * 1.00
      zz0 =   r_in ( irow, 2 )
      zphil = r_in (irow, 3 )
      zphilr = zphil * zpi / 180.0
      zalfa = r_in ( irow, 4 ) + zpol*90.0
      zalfa = zalfa * zpi / 180.0
      zlen =  r_in ( irow, 5 ) * 1.00
      
      zsalfa = SIN ( zalfa )
      zcalfa = COS ( zalfa )
      
      zdpts = 1.0/(nzspt-1)
      zdtsl = zlen * zdpts
c...Put zdtsl in vacuum7.inc now
      dtsl = zdtsl
      DO i = 1, nzspt
         zth = - zlen * ( 1.0/2.0 - (i-1)*zdpts )
         slx(i) = zx0 + zth * zsalfa
         slz(i) = zz0 - zth * zcalfa
         slxt(i) = + zsalfa
         slzt(i) = - zcalfa
      END DO

      WRITE ( 6, '( /, "zx0, zz0, zalfa, zlen, slzt(1), zphil = ",
     $     /, 6es13.5 )' ) zx0, zz0, zalfa, zlen, slzt(1), zphil
      WRITE ( 23, '( /, "zx0, zz0, zalfa, zlen, slzt(1), zphil = ",
     $     /, 6es13.5 )' ) zx0, zz0, zalfa, zlen, slzt(1), zphil
      WRITE ( 28, '( /, "zx0, zz0, zalfa, zlen, slzt(1), zphil = ",
     $     /, 6es13.5 )' ) zx0, zz0, zalfa, zlen, slzt(1), zphil

      WRITE ( 6, '( "XXXXX  R_LOOPS nzspt = ", i5 )' ) nzspt
      WRITE ( 23, '( "XXXXX  RLOOPS nzspt = ", i5 )' ) nzspt

         CALL drawc6(xwal,zwal,xinf,zinf,xcw,zcw,xcwa,zcwa,xcwb,zcwb,
     $     slx,slz,
     $        1,mth1,1,ntip_c,1,ntip_u,1,ntip_l,1,nzspt,
     $        "z","x",xmx,xma,zma,
     $        plrad,1,1,idot, s_name, ishape,a,b,aw,bw,cw,dw,tw,
     $        abulg, bbulg, tbulg, wcentr, pcircm, wcircm,
     $        pelong, pdness, wmaj,wlrad, welong, wdness,
     $        chgt,chgta,chgtb )

      END SUBROUTINE r_loops_row

c........................................................
      SUBROUTINE sensflxi (slx, slxt, slz, slzt, nzspt, zphil, lpols,
     $     s_name )
c........................................................

c... The radial or poloidal flux in the sensor loops. There are 
c    contributions from the feedback coil, plasma, and wall.

      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum7.inc'

      REAL, DIMENSION(:), ALLOCATABLE :: gccoila, gccoilb, gccoil
      REAL, DIMENSION(*) :: slx, slxt, slz, slzt
      CHARACTER*(*) s_name

      ALLOCATE ( gccoila(mtcoila + 5), gccoilb(mtcoilb + 5),
     $     gccoil(mtcoil + 5) )

c..   null vectors for testing:
c$$$  ccoila = 0.0
c$$$  ccoilb = 0.0
c$$$  skcoilpa = 0.0
c$$$  skcoilpb = 0.0
c$$$  skcoilwa = 0.0
c$$$  skcoilwb = 0.0

c...  Multiply the ccoil by the gain to be consistent with skcoil:

      gccoila(1:mtcoila+1) = gainca * ccoila(1:mtcoila+1)

c... Separate contributions from the the coil, plasma and wall.

      zlens = clen_u
      IF ( lfbcoila == 1 ) THEN
         mtcoil0 = mtcoila
         CALL sensorcc ( slx, slxt, slz, slzt, nzspt, zphil, lpols,
     $        gccoila,
     $        xcwa, zcwa, xcwap, zcwap, zlens, mtcoil0, rsensflxaa,
     $        psensflxaa, s_name, "Sensflx-L-L" )
         IF ( lpols == 0 ) THEN
            sensflxaa = rsensflxaa
         ELSE IF ( lpols ==1 ) THEN
            sensflxaa = psensflxaa
         END IF
      END IF
      IF ( lfbcoila == 2 ) THEN
         mtcoil0 = ntip_u
         CALL sensorcc2 ( slx, slxt, slz, slzt, nzspt, zphil, lpols,
     $        gccoila,
     $        xcwa, zcwa, xcwap, zcwap, zlens, mtcoil0, rsensflxaa,
     $        psensflxaa, s_name, "Sensflx-L-L" )
         IF ( lpols == 0 ) THEN
            sensflxaa = rsensflxaa
         ELSE IF ( lpols ==1 ) THEN
            sensflxaa = psensflxaa
         END IF
      END IF
      CALL sensorcp ( slx, slxt, slz, slzt, nzspt, zphil, lpols,
     $     skcoilpa,
     $     xcwa, zcwa, xcwap, zcwap, mtcoila, rsensflxap, psensflxap,
     $     s_name, "Sensflx-U-P" )
         IF ( lpols == 0 ) THEN
            sensflxap = rsensflxap
         ELSE IF ( lpols ==1 ) THEN
            sensflxap = psensflxap
         END IF
      CALL sensorcw ( slx, slxt, slz, slzt, nzspt, zphil, lpols,
     $        skcoilwa,
     $     xcwa, zcwa, xcwap, zcwap, mtcoila, rsensflxaw, psensflxaw,
     $     s_name, "Sensflx-L-W" )
      IF ( lpols == 0 ) THEN
         sensflxaw = rsensflxaw
      ELSE IF ( lpols ==1 ) THEN
         sensflxaw = psensflxaw
      END IF

      gccoilb(1:mtcoilb+1) = gaincb * ccoilb(1:mtcoilb+1)

      zlens = clen_l
      IF ( lfbcoilb == 1 ) THEN
         mtcoil0 = mtcoilb
         CALL sensorcc ( slx, slxt, slz, slzt, nzspt, zphil, lpols,
     $        gccoilb,
     $        xcwb, zcwb, xcwbp, zcwbp, zlens, mtcoil0, rsensflxbb,
     $        psensflxbb, s_name, "Sensflx-U-U" )
         IF ( lpols == 0 ) THEN
            sensflxbb = rsensflxbb
         ELSE IF ( lpols ==1 ) THEN
            sensflxbb = psensflxbb
         END IF
      END IF
      IF ( lfbcoilb == 2 ) THEN
         mtcoil0 = ntip_l
         CALL sensorcc2 ( slx, slxt, slz, slzt, nzspt, zphil, lpols,
     $        gccoilb,
     $        xcwb, zcwb, xcwbp, zcwbp, zlens, mtcoil0, rsensflxbb,
     $        psensflxbb, s_name, "Sensflx-U-U" )
         IF ( lpols == 0 ) THEN
            sensflxbb = rsensflxbb
         ELSE IF ( lpols ==1 ) THEN
            sensflxbb = psensflxbb
         END IF
      END IF
      CALL sensorcp ( slx, slxt, slz, slzt, nzspt, zphil, lpols,
     $     skcoilpb,
     $     xcwb, zcwb, xcwbp, zcwbp, mtcoilb, rsensflxbp, psensflxbp,
     $     s_name, "Sensflx-L-P" )
      IF ( lpols == 0 ) THEN
         sensflxbp = rsensflxbp
      ELSE IF ( lpols ==1 ) THEN
         sensflxbp = psensflxbp
      END IF
      CALL sensorcw ( slx, slxt, slz, slzt, nzspt, zphil, lpols,
     $     skcoilwb,
     $     xcwb, zcwb, xcwbp, zcwbp, mtcoilb, rsensflxbw, psensflxbw,
     $     s_name, "Sensflx-U-W" )
      IF ( lpols == 0 ) THEN
         sensflxbw = rsensflxbw
      ELSE IF ( lpols ==1 ) THEN
         sensflxbw = psensflxbw
      END IF


c... C-coil:

      gccoil(1:mtcoil+1) = gainc * ccoil(1:mtcoil+1)

      zlens = clen_c
      IF ( lfbcoil == 1 ) THEN
         mtcoil0 = mtcoil
         CALL sensorcc ( slx, slxt, slz, slzt, nzspt, zphil, lpols,
     $        gccoil,
     $        xcw, zcw, xcwp, zcwp, zlens, mtcoil0, rsensflxcc,
     $        psensflxcc, s_name, "Sensflx-C-C" )
         IF ( lpols == 0 ) THEN
            sensflxcc = rsensflxcc
         ELSE IF ( lpols ==1 ) THEN
            sensflxcc = psensflxcc
         END IF
      END IF         
      IF ( lfbcoil == 2 ) THEN
         mtcoil0 = ntip_c
         CALL sensorcc2 ( slx, slxt, slz, slzt, nzspt, zphil, lpols,
     $        gccoil,
     $        xcw, zcw, xcwp, zcwp, zlens, ntip_c, rsensflxcc,
     $        psensflxcc, s_name, "Sensflx-C-C" )
         IF ( lpols == 0 ) THEN
            sensflxcc = rsensflxcc
         ELSE IF ( lpols ==1 ) THEN
            sensflxcc = psensflxcc
         END IF
      END IF

      sensflxcp = 0.0
      sensflxcw = 0.0

      DEALLOCATE ( gccoila, gccoilb, gccoil )

c   Add contributions from A and B coils.

      sensflxa = sensflxaa + sensflxap + sensflxaw
      sensflxb = sensflxbb + sensflxbp + sensflxbw
      sensflxc = sensflxcc + sensflxcp + sensflxcw

      WRITE ( outmod, '(/,5x,"RADIAL SENSOR FLUX From LL Coil = ",
     $     es11.4 )' ) rsensflxaa
      WRITE ( 28, '(/,5x,"RADIAL SENSOR FLUX From LL Coil = ",
     $     es11.4 )' ) rsensflxaa

      WRITE ( outmod, '(5x,"RADIAL SENSOR FLUX From UU Coil = ",
     $     es11.4 )' ) rsensflxbb
      WRITE ( 28, '(5x,"RADIAL SENSOR FLUX From UU Coil = ",
     $     es11.4 )' ) rsensflxbb

      WRITE ( outmod, '(5x,"RADIAL SENSOR FLUX From CC Coil = ",
     $     es11.4, / )' ) rsensflxcc
      WRITE ( 28, '(5x,"RADIAL SENSOR FLUX From CC Coil = ",
     $     es11.4, /  )' ) rsensflxcc
       
      WRITE ( outmod, '( 5x,/, "LPOLSENSOR = ", i5 )' ) lpols
      WRITE (  iotty, '( 5x,/, "LPOLSENSOR = ", i5 )' ) lpols
      WRITE (     28, '( 5x,/, "LPOLSENSOR = ", i5 )' ) lpols
      
      WRITE ( outmod, '( /, 7x, "flx-LL", 7x, "flx-LP", 7x, "flx-LW",
     $     7x, "flx-A" )' )
      WRITE ( outmod, '( 1x, 4es13.5 )' )
     $     sensflxaa, sensflxap, sensflxaw, sensflxa

      WRITE ( outmod, '( /, 7x, "flx-UU", 7x, "flx-UP", 7x, "flx-UW",
     $     7x, "flx-B" )' )
      WRITE ( outmod, '( 1x, 4es13.5 )' )
     $     sensflxbb, sensflxbp, sensflxbw, sensflxb
      
      WRITE ( outmod, '( /, 7x, "flx-CC", 7x, "flx-CP", 7x, "flx-CW",
     $     7x, "flx-C" )' )
      WRITE ( outmod, '( 1x, 4es13.5 )' )
     $     sensflxcc, sensflxcp, sensflxcw, sensflxc
      
      WRITE ( iotty, '( /, 7x, "flx-LL", 7x, "flx-LP", 7x, "flx-LW",
     $     7x, "flx-A" )' )
      WRITE ( iotty, '( 1x, 4es13.5 )' )
     $     sensflxaa, sensflxap, sensflxaw, sensflxa

      WRITE ( iotty, '( /, 7x, "flx-UU", 7x, "flx-UP", 7x, "flx-UW",
     $     7x, "flx-B" )' )
      WRITE ( iotty, '( 1x, 4es13.5 )' )
     $     sensflxbb, sensflxbp, sensflxbw, sensflxb

      WRITE ( iotty, '( /, 7x, "flx-CC", 7x, "flx-CP", 7x, "flx-CW",
     $     7x, "flx-C" )' )
      WRITE ( iotty, '( 1x, 4es13.5 )' )
     $     sensflxcc, sensflxcp, sensflxcw, sensflxc

      WRITE ( 28, '( /, 7x, "flx-LL", 7x, "flx-LP", 7x, "flx-LW",
     $     7x, "flx-A" )' )
      WRITE ( 28, '( 1x, 4es13.5 )' )
     $     sensflxaa, sensflxap, sensflxaw, sensflxa

      WRITE ( 28, '( /, 7x, "flx-UU", 7x, "flx-UP", 7x, "flx-UW",
     $     7x, "flx-B" )' )
      WRITE ( 28, '( 1x, 4es13.5 )' )
     $     sensflxbb, sensflxbp, sensflxbw, sensflxb

      WRITE ( 28, '( /, 7x, "flx-CC", 7x, "flx-CP", 7x, "flx-CW",
     $     7x, "flx-C" )' )
      WRITE ( 28, '( 1x, 4es13.5 )' )
     $     sensflxcc, sensflxcp, sensflxcw, sensflxc

      WRITE ( outmod, '()' )
      WRITE ( outmod, '( "Sensor Flux from Coil A: ", es13.5)' )
     $     sensflxa
      WRITE ( outmod, '( "Sensor Flux from Coil B: ", es13.5)' )
     $     sensflxb

      WRITE ( iotty, '()' )
      WRITE ( iotty, '( "Sensor Flux from Coil A: ", es13.5)' )
     $     sensflxa
      WRITE ( iotty, '( "Sensor Flux from Coil B: ", es13.5)' )
     $     sensflxb

      WRITE ( 6, '( "XXXXX  SENSFLXI nzspt = ", i5 )' ) nzspt
      WRITE ( 23, '( "XXXXX  SENSFLXI nzspt = ", i5 )' ) nzspt

      END SUBROUTINE sensflxi
      
c............................................................
      SUBROUTINE sensorcc ( slx, slxt, slz, slzt, nzspt, zphil, lpols,
     $     hcr, zxcl,zzcl,zxclp,zzclp,
     $     zlens, mzcoil, flxr, pflxr, s_name, cname )
c............................................................
c     
      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum5.inc'
      INCLUDE 'vacuum7.inc'
c     
      CHARACTER*(*) s_name, cname
      CHARACTER(8), DIMENSION(nths) :: zclab
      CHARACTER(8), DIMENSION(1) :: zhlab

      DIMENSION xpla(nths), zpla(nths)

      DIMENSION grdgre(nths2,nths2),grwp(nths,nths),
     .     grri(nths2,nfm2)

      REAL, DIMENSION(*) :: hcr, zxcl,zzcl, zxclp,zzclp
      REAL, DIMENSION(*) :: slx, slxt, slz, slzt

c$$$      REAL, DIMENSION(:), ALLOCATABLE :: slx, slxt, slz ,slzt
      REAL, DIMENSION(:), ALLOCATABLE :: bxr, bzr, bpi, chir0,
     $     zxgrd, zbamp, zbphs
      REAL, DIMENSION(:,:), ALLOCATABLE :: chir, cwrkr
c      REAL, DIMENSION(1) :: zhmax

      COMMON / bigv1 / grdgre, grwp, grri
c     
      COMMON / gridc / xobs(ndimlpc),zobs(ndimlpc)
c     
      EQUIVALENCE (xpla,xinf), (zpla,zinf)

      CHARACTER(130), DIMENSION(10) :: string

c     nzspt is the number of observer points.

      WRITE ( 28, '()' )
      WRITE ( 28, '("****** Calculating: ", a, 2x, a)' )
     $     s_name, cname

      WRITE ( outmod, '()' )
      WRITE ( outmod, '("****** Calculating: ", a, 2x, a)' )
     $     s_name, cname
      WRITE ( iotty, '()' )
      WRITE ( iotty, '("****** Calculating: ", a, 2x, a)' )
     $     s_name, cname

      CALL vecwrt( mzcoil, hcr, "hcr", 1, mzcoil, 28, 0 )

      zhmax = MAXVAL ( hcr(1:mzcoil) )
      WRITE ( 28, '( "CURRENT PEAK = ", 8ES12.4, "amperes" )' ) zhmax

      zclab = " "

      mth12 = 2*mth
      IF ( lfarw .gt. 0 ) mth12 = mth
c     
c$$$      ALLOCATE ( slx(nzspt), slxt(nzspt), slz(nzspt), slzt(nzspt) )
      ALLOCATE ( chir(4,nzspt), cwrkr(4,nzspt), chir0(nzspt) )

c$$$      CALL sloops ( slx, slxt, slz ,slzt, nzspt)
c     
      delx = plrad * deloop
      delz = plrad * deloop
c     
c.....get array of observer points, and initialize chi matrix.
c     
c     
      DO nsew = 1, 4
c     
         DO i = 1, nzspt
c     
            chir(nsew,i) = 0.0
            cwrkr(nsew,i) = 0.0
c     
            GO TO ( 51, 52, 53, 54), nsew
c     
 51         CONTINUE
c......north
c     
            xobs(i) = slx(i)
            zobs(i) = slz(i) + delz
            GO TO 90
c     
 52         CONTINUE
c.....south
c     
            xobs(i) = slx(i)
            zobs(i) = slz(i) - delz
            GO TO 90
c     
 53         CONTINUE
c.....east
c     
            xobs(i) = slx(i) + delx
            zobs(i) = slz(i)
            GO TO 90
c     
 54         CONTINUE
c.....west
c     
            xobs(i) = slx(i) - delx
            zobs(i) = slz(i)
c     
 90         CONTINUE

         END DO                 ! Loop on I, nzspt

c     Integrate over coil.

         isg = 1
         WRITE ( outmod,'("isgn for coil source contrib. = ", i3)')
     $        isg
         WRITE ( iotty,'("isgn for coil source contrib. = ", i3)')
     $        isg
c     
         DO i = 1, nzspt
            cwrkr(nsew,i) = 0.0
         END DO

         ns = mzcoil
         CALL chic1 ( zxcl,zzcl, zxclp,zzclp, isg,hcr,nzspt,
     $        ns,1,cwrkr,nsew )
         CALL atpoint("After Chic1","lfarw",lfarw,"delx",delx,
     $        iotty, outmod )

c..   Store away the coil result.
c     
         DO  i = 1, nzspt
            chir(nsew,i) = cwrkr(nsew,i)
         END DO
c     
         WRITE (outmod,'("Coil Contrib.to Probes",i4)') nsew
         DO  i = nzspt/2, nzspt/2+2
            WRITE(outmod,'("i,chir, cwrkr,",i4,1p2e12.3)')
     $           i,chir(nsew,i),cwrkr(nsew,i)
         END DO

c     
 200     CONTINUE

      END DO                    ! Loop over NSEW

c... Calculate Chi at observer points.
     
      DO i = 1, nzspt
         xobs(i) = slx(i)
         zobs(i) = slz(i)
      END DO
      DO i = 1, nzspt
         cwrkr(1,i) = 0.0
      END DO
      ns = mzcoil
      CALL chic1 ( zxcl,zzcl, zxclp,zzclp, isg,hcr,nzspt,
     $     ns,1,cwrkr,1 )
      DO  i = 1, nzspt
         chir0(i) = cwrkr(1,i)
      END DO

c.....calculate bx and bz
c     
      ALLOCATE ( bxr(nzspt), bzr(nzspt), bpi(nzspt) )
      ALLOCATE ( zxgrd(nzspt), zbamp(nzspt), zbphs(nzspt) )
c     
      DO i = 1, nzspt
         bxr(i) = 0.0
         bzr(i) = 0.0
         bpi(i) = 0.0
      END DO

      DO i = 1, nzspt
c     
         bxr(i) = ( chir(3,i) - chir(4,i) ) / (2.0*delx)
         bzr(i) = ( chir(1,i) - chir(2,i) ) / (2.0*delz)
         bpi(i) = - n * chir0(i) / slx(i)
         WRITE ( outmod,'("bxr,bzr,bpi = ", i3, 1p4e12.3)')
     $        i, bxr(i),bzr(i),bpi(i)
      END DO

      zdeg = 180.0/pye
      zmu0 =  4.0 * pye * 1.0e-7
      zmu0f = zmu0

      bxr = zmu0f * bxr
      bzr = zmu0f * bzr
      bpi = zmu0f * bpi

         CALL rflux ( slx,slxt, slzt, bxr,bzr, nzspt, dtsl, flxr )
         CALL rfluxp ( slx,slxt, slzt, bxr,bzr, nzspt, dtsl, pflxr )

      WRITE ( 6, '( "XXXXX  SENSORCC nzspt = ", i5 )' ) nzspt
      WRITE ( 23, '( "XXXXX  SENSORCC nzspt = ", i5 )' ) nzspt


c... Write out about 50 or so values to SENSOR-MEAS, and flag mid value

      zclab( (nzspt-1)/2+1 ) = " xx"//cname(8:11)
      zhlab(1) = "  h"//cname(8:11)

      WRITE ( 28, '(20x, ".......",  A, 2x, A, "......." )' )
     $     s_name, cname
      WRITE ( 28, '("%", 4x, "i", 6x, "X", 11x, "Z", 11x, "Phi", 11x,
     $     "B_x", 10x, "B_z", 8x, "B_p", a )' ) zhlab(1)
      DO i = 1, nzspt
         zbamp(i) = SQRT ( bxr(i)**2 + bzr(i)**2 + bpi(i)**2 )
         zbphs(i) = ATAN2 ( bzr(i),bxr(i) ) * zdeg
         WRITE ( 28, '( 1x, i5, 6ES12.3, a )' ) i, slx(i), slz(i),
     $        zphil,  bxr(i), bzr(i), bpi(i), zclab(i)
      END DO

      zangl = ATAN2 ( slz(nzspt)-slz(1),slx(nzspt)-slx(1) )
      zxgrd(1:nzspt) =
     $     slx(1:nzspt)*COS(zangl) + slz(1:nzspt)*SIN(zangl)

      CALL pospl1 ( zxgrd, bxr, 1,nzspt, 1,"B_x", nzspt )
      CALL pospl1 ( zxgrd, bzr, 1,nzspt, 2,"B_z", nzspt )
      CALL pospl1 ( zxgrd, bpi, 1,nzspt, 3,"I-B_P", nzspt )
      CALL pospl1 ( zxgrd, zbphs, 1,nzspt, 4,"PHASE-B", nzspt )
      
      call frame( 1 )

c$$$      DEALLOCATE ( slx, slxt, slz, slzt )
      DEALLOCATE ( bxr, bzr, bpi, zxgrd )
      DEALLOCATE ( chir, chir0, cwrkr )

 500  CONTINUE

      END SUBROUTINE sensorcc

c............................................................
      SUBROUTINE sensorcc2 ( slx, slxt, slz, slzt, nzspt, zphil, lpols,
     $     hcr, zxcl,zzcl,zxclp,zzclp,
     $     zlens, mzcoil, flxr, pflxr, s_name, cname )
c............................................................
c     
      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum5.inc'
      INCLUDE 'vacuum7.inc'
c     
      CHARACTER*(*)  s_name, cname
      CHARACTER(8), DIMENSION(nths) :: zclab
      CHARACTER(8), DIMENSION(1) :: zhlab

      DIMENSION xpla(nths), zpla(nths)

      DIMENSION grdgre(nths2,nths2),grwp(nths,nths),
     .     grri(nths2,nfm2)

      REAL, DIMENSION(*) :: hcr, zxcl,zzcl, zxclp,zzclp
      REAL, DIMENSION(*) :: slx, slxt, slz, slzt

c$$$      REAL, DIMENSION(:), ALLOCATABLE :: slx, slxt, slz ,slzt
      REAL, DIMENSION(:), ALLOCATABLE :: bxr, bzr, bpi, chir0,
     $     zxgrd, zbamp, zbphs
      REAL, DIMENSION(:,:), ALLOCATABLE :: chir, cwrkr

      COMMON / bigv1 / grdgre, grwp, grri
c     
      COMMON / gridc / xobs(ndimlpc),zobs(ndimlpc)
c     
      EQUIVALENCE (xpla,xinf), (zpla,zinf)

      CHARACTER(130), DIMENSION(10) :: string

c     nzspt is the number of observer points.

      WRITE ( 28, '()' )
      WRITE ( 28, '("****** Calculating: ", a, 2x, a)' )
     $     s_name, cname

      WRITE ( outmod, '()' )
      WRITE ( outmod, '("****** Calculating: ", a, 2x, a)' )
     $     s_name, cname
      WRITE ( iotty, '()' )
      WRITE ( iotty, '("****** Calculating: ", a, 2x, a)' )
     $     s_name, cname

      WRITE ( 28, '( "CURRENT = ", 8ES12.4, "amperes" )' ) hcr(1)

      zclab = " "

      mth12 = 2*mth
      IF ( lfarw .gt. 0 ) mth12 = mth
c     
c$$$      ALLOCATE ( slx(nzspt), slxt(nzspt), slz(nzspt), slzt(nzspt) )
      ALLOCATE ( chir(4,nzspt), cwrkr(4,nzspt), chir0(nzspt) )

c$$$      CALL sloops ( slx, slxt, slz ,slzt, nzspt)
c     
      delx = plrad * deloop
      delz = plrad * deloop
c     
c.....get array of observer points, and initialize chi matrix.
c     
c     
      DO nsew = 1, 4
c     
         DO i = 1, nzspt
c     
            chir(nsew,i) = 0.0
            cwrkr(nsew,i) = 0.0
c     
            GO TO ( 51, 52, 53, 54), nsew
c     
 51         CONTINUE
c......north
c     
            xobs(i) = slx(i)
            zobs(i) = slz(i) + delz
            GO TO 90
c     
 52         CONTINUE
c.....south
c     
            xobs(i) = slx(i)
            zobs(i) = slz(i) - delz
            GO TO 90
c     
 53         CONTINUE
c.....east
c     
            xobs(i) = slx(i) + delx
            zobs(i) = slz(i)
            GO TO 90
c     
 54         CONTINUE
c.....west
c     
            xobs(i) = slx(i) - delx
            zobs(i) = slz(i)
c     
 90         CONTINUE

         END DO                 ! Loop on I, nzspt

c     Integrate over coil.

         isg = 1
         WRITE ( outmod,'("isgn for coil source contrib. = ", i3)')
     $        isg
         WRITE ( iotty,'("isgn for coil source contrib. = ", i3)')
     $        isg
c     
         DO i = 1, nzspt
            cwrkr(nsew,i) = 0.0
         END DO

         ns = mzcoil-1
         CALL chic2 ( zxcl,zzcl, zxclp,zzclp, isg,hcr,nzspt,
     $        zlens, ns,1,cwrkr,nsew )
         CALL atpoint("After Chic2","lfarw",lfarw,"delx",delx,
     $        iotty, outmod )

c..   Store away the coil result.
c     
         DO  i = 1, nzspt
            chir(nsew,i) = cwrkr(nsew,i)
         END DO
c     
         WRITE (outmod,'("Coil Contrib.to Probes",i4)') nsew
         DO  i = nzspt/2, nzspt/2+2
            WRITE(outmod,'("i,chir, cwrkr,",i4,1p2e12.3)')
     $           i,chir(nsew,i),cwrkr(nsew,i)
         END DO

c     
 200     CONTINUE

      END DO                    ! Loop over NSEW

c... Calculate Chi at observer points.
     
      DO i = 1, nzspt
         xobs(i) = slx(i)
         zobs(i) = slz(i)
      END DO
      DO i = 1, nzspt
         cwrkr(1,i) = 0.0
      END DO
      ns = mzcoil-1
      CALL chic2 ( zxcl,zzcl, zxclp,zzclp, isg,hcr,nzspt,
     $     zlens,ns,1,cwrkr,1 )
      DO  i = 1, nzspt
         chir0(i) = cwrkr(1,i)
      END DO

c.....calculate bx and bz
c     
      ALLOCATE ( bxr(nzspt), bzr(nzspt), bpi(nzspt) )
      ALLOCATE ( zxgrd(nzspt), zbamp(nzspt), zbphs(nzspt) )
c     
      DO i = 1, nzspt
         bxr(i) = 0.0
         bzr(i) = 0.0
         bpi(i) = 0.0
      END DO

      DO i = 1, nzspt
c     
         bxr(i) = ( chir(3,i) - chir(4,i) ) / (2.0*delx)
         bzr(i) = ( chir(1,i) - chir(2,i) ) / (2.0*delz)
         bpi(i) = - n * chir0(i) / slx(i)
         WRITE ( outmod,'("bxr,bzr,bpi = ", i3, 1p4e12.3)')
     $        i, bxr(i),bzr(i),bpi(i)
      END DO

      zdeg = 180.0/pye
      zmu0 =  4.0 * pye * 1.0e-7
      zmu0f = zmu0 

      bxr = zmu0f * bxr
      bzr = zmu0f * bzr
      bpi = zmu0f * bpi

      WRITE ( outmod, '( 5x, "LPOLSENSOR = ", i5 )' ) lpols
      WRITE ( iotty, '( 5x, "LPOLSENSOR = ", i5 )' ) lpols
      WRITE (    28, '( 5x, "LPOLSENSOR = ", i5 )' ) lpols

      CALL rflux ( slx,slxt, slzt, bxr,bzr, nzspt, dtsl, flxr )
      CALL rfluxp ( slx,slxt, slzt, bxr,bzr, nzspt, dtsl, pflxr )

      WRITE ( 6, '( "XXXXX  SENSORCC nzspt = ", i5 )' ) nzspt
      WRITE ( 23, '( "XXXXX  SENSORCC nzspt = ", i5 )' ) nzspt

c... Write out about 50 or so values to SENSOR-MEAS, and flag mid value.

      zclab( (nzspt-1)/2+1 ) = " xx"//cname(8:11)
      zhlab(1) = "  h"//cname(8:11)

      WRITE ( 28, '(20x, ".......", A, 2x, A, "......." )' )
     $     s_name, cname
      WRITE ( 28, '("%", 4x, "i", 6x, "X", 11x, "Z", 11x, "Phi", 9x,
     $     "B_x", 10x, "B_z", 8x, "B_p", a  )' ) zhlab(1)
      DO i = 1, nzspt
         zbamp(i) = SQRT ( bxr(i)**2 + bzr(i)**2 + bpi(i)**2 )
         zbphs(i) = ATAN2 ( bzr(i),bxr(i) ) * zdeg
         WRITE ( 28, '( 1x, i5, 6ES12.3, a )' ) i, slx(i), slz(i),
     $        zphil, bxr(i), bzr(i), bpi(i), zclab(i)
      END DO

      zangl = ATAN2 ( slz(nzspt)-slz(1),slx(nzspt)-slx(1) )
      zxgrd(1:nzspt) =
     $     slx(1:nzspt)*COS(zangl) + slz(1:nzspt)*SIN(zangl)

      CALL pospl1 ( zxgrd, bxr, 1,nzspt, 1,"B_x", nzspt )
      CALL pospl1 ( zxgrd, bzr, 1,nzspt, 2,"B_z", nzspt )
      CALL pospl1 ( zxgrd, bpi, 1,nzspt, 3,"I-B_P", nzspt )
      CALL pospl1 ( zxgrd, zbphs, 1,nzspt, 4,"PHASE-B", nzspt )
      
      call frame( 1 )

c$$$      DEALLOCATE ( slx, slxt, slz, slzt )
      DEALLOCATE ( bxr, bzr, bpi, zxgrd )
      DEALLOCATE ( chir, chir0, cwrkr )

 500  CONTINUE

      END SUBROUTINE sensorcc2

c............................................................
      SUBROUTINE sensorcp ( slx, slxt, slz, slzt, nzspt, zphil, lpols,
     $     hpr, zxcl,zzcl,zxclp,zzclp,
     $     mzcoil, flxr, pflxr, s_name, cname )
c............................................................
c     
      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum5.inc'
      INCLUDE 'vacuum7.inc'
c     
      CHARACTER*(*) s_name, cname

      DIMENSION xpla(nths), zpla(nths)

      DIMENSION grdgre(nths2,nths2),grwp(nths,nths),
     .     grri(nths2,nfm2)

      REAL, DIMENSION(*) :: hpr, zxcl,zzcl, zxclp,zzclp
      REAL, DIMENSION(*) :: slx, slxt, slz, slzt

c$$$      REAL, DIMENSION(:), ALLOCATABLE :: slx, slxt, slz ,slzt
      REAL, DIMENSION(:), ALLOCATABLE :: bxr, bzr
      REAL, DIMENSION(:,:), ALLOCATABLE :: chir, cwrkr

      COMMON / bigv1 / grdgre, grwp, grri
c     
      COMMON / gridc / xobs(ndimlpc),zobs(ndimlpc)
c     
      EQUIVALENCE (xpla,xinf), (zpla,zinf)

      CHARACTER(130), DIMENSION(10) :: string

c     nzspt is the number of observer points.

      WRITE ( outmod, '()' )
      WRITE ( outmod, '("****** Calculating: ", a, 2x, a)' )
     $     s_name, cname
      WRITE ( iotty, '()' )
      WRITE ( iotty, '("****** Calculating: ", a, 2x, a)' )
     $     s_name, cname

      mth12 = 2*mth
      IF ( lfarw .gt. 0 ) mth12 = mth
c     
c$$$      ALLOCATE ( slx(nzspt), slxt(nzspt), slz(nzspt), slzt(nzspt) )
      ALLOCATE ( chir(4,nzspt), cwrkr(4,nzspt) )

c$$$      CALL sloops ( slx, slxt, slz ,slzt, nzspt)
c     
      delx = plrad * deloop
      delz = plrad * deloop
c     
c.....get array of observer points, and initialize chi matrix.
c     
c     
      DO nsew = 1, 4
c     
         DO i = 1, nzspt
c     
            chir(nsew,i) = 0.0
            cwrkr(nsew,i) = 0.0
c     
            GO TO ( 51, 52, 53, 54), nsew
c     
 51         CONTINUE
c......north
c     
            xobs(i) = slx(i)
            zobs(i) = slz(i) + delz
            GO TO 90
c     
 52         CONTINUE
c.....south
c     
            xobs(i) = slx(i)
            zobs(i) = slz(i) - delz
            GO TO 90
c     
 53         CONTINUE
c.....east
c     
            xobs(i) = slx(i) + delx
            zobs(i) = slz(i)
            GO TO 90
c     
 54         CONTINUE
c.....west
c     
            xobs(i) = slx(i) - delx
            zobs(i) = slz(i)
c     
 90         CONTINUE

         END DO                 ! Loop on I, nzspt

c     Integrate over plasma.

         isg = 1
         WRITE ( outmod,'("isgn for plasma source contrib. = ", i3)')
     $        isg
         WRITE ( iotty,'("isgn for plasma source contrib. = ", i3)')
     $        isg
c     
         DO i = 1, nzspt
            cwrkr(nsew,i) = 0.0
         END DO

         ns = mth
         CALL chic1 ( xpla,zpla,xplap,zplap,isg,hpr,nzspt,
     $        ns,1,cwrkr,nsew )
         CALL atpoint("After Chic1","lfarw",lfarw,"delx",delx,
     $        iotty, outmod )

c..   Add the plasma result.

         DO  i = 1, nzspt
            chir(nsew,i) = cwrkr(nsew,i)
         END DO
c     
         WRITE (outmod,'("Coil-Plasma Contrib.to Probes",i4)')nsew
         DO  i = nzspt/2, nzspt/2+2
            WRITE(outmod,'("i,chir, cwrkr,",i4,1p2e12.3)')
     $           i,chir(nsew,i),cwrkr(nsew,i)
         END DO
c     
 200     CONTINUE

      END DO                    ! Loop over NSEW

c.....calculate bx and bz
c     
      ALLOCATE ( bxr(nzspt), bzr(nzspt) )
c     
      DO i = 1, nzspt
         bxr(i) = 0.0
         bzr(i) = 0.0
      END DO

      DO i = 1, nzspt
c     
         bxr(i) = ( chir(3,i) - chir(4,i) ) / (2.0*delx)
         bzr(i) = ( chir(1,i) - chir(2,i) ) / (2.0*delz)
         WRITE ( outmod,'("bxr,zr = ", i3, 1p4e12.3)')
     $        i, bxr(i),bzr(i)
      END DO

      zmu0 =  4.0 * pye * 1.0e-7
      zmu0f = zmu0

      bxr = zmu0f * bxr
      bzr = zmu0f * bzr
      bpi = zmu0f * bpi

      CALL rflux ( slx,slxt, slzt, bxr,bzr, nzspt, dtsl, flxr )
      CALL rfluxp ( slx,slxt, slzt, bxr,bzr, nzspt, dtsl, pflxr )

      WRITE ( 6, '( "XXXXX  SENSORCP nzspt = ", i5 )' ) nzspt
      WRITE ( 23, '( "XXXXX  SENSORCP nzspt = ", i5 )' ) nzspt

c$$$      DEALLOCATE ( slx, slxt, slz, slzt )
      DEALLOCATE ( bxr, bzr )
      DEALLOCATE ( chir, cwrkr )

 500  CONTINUE

      END SUBROUTINE sensorcp

c............................................................
      SUBROUTINE sensorcw ( slx, slxt, slz, slzt, nzspt, zphil, lpols,
     $     hwr, zxcl,zzcl,zxclp,zzclp,
     $     mzcoil, flxr, pflxr, s_name, cname )
c............................................................
c     
      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum5.inc'
      INCLUDE 'vacuum7.inc'
c     
      CHARACTER*(*) s_name, cname

      DIMENSION xpla(nths), zpla(nths)

      DIMENSION grdgre(nths2,nths2),grwp(nths,nths),
     .     grri(nths2,nfm2)

      REAL, DIMENSION(*) :: hwr, zxcl,zzcl, zxclp,zzclp
      REAL, DIMENSION(*) :: slx, slxt, slz, slzt

c$$$      REAL, DIMENSION(:), ALLOCATABLE :: slx, slxt, slz ,slzt
      REAL, DIMENSION(:), ALLOCATABLE :: bxr, bzr
      REAL, DIMENSION(:,:), ALLOCATABLE :: chir, cwrkr

      COMMON / bigv1 / grdgre, grwp, grri
c     
      COMMON / gridc / xobs(ndimlpc),zobs(ndimlpc)
c     
      EQUIVALENCE (xpla,xinf), (zpla,zinf)

      CHARACTER(130), DIMENSION(10) :: string

c     nzspt is the number of observer points.

      WRITE ( outmod, '()' )
      WRITE ( outmod, '("****** Calculating: ", a, 2x, a)' )
     $     s_name, cname
      WRITE ( iotty, '()' )
      WRITE ( iotty, '("****** Calculating: ", a, 2x, a)' )
     $     s_name, cname

      mth12 = 2*mth
      IF ( lfarw .gt. 0 ) mth12 = mth
c     
c$$$      ALLOCATE ( slx(nzspt), slxt(nzspt), slz(nzspt), slzt(nzspt) )
      ALLOCATE ( chir(4,nzspt), cwrkr(4,nzspt) )

c$$$      CALL sloops ( slx, slxt, slz ,slzt, nzspt)
c     
      delx = plrad * deloop
      delz = plrad * deloop
c     
c.....get array of observer points, and initialize chi matrix.
c     
c     
      DO nsew = 1, 4
c     
         DO i = 1, nzspt
c     
            chir(nsew,i) = 0.0
            cwrkr(nsew,i) = 0.0
c     
            GO TO ( 51, 52, 53, 54), nsew
c     
 51         CONTINUE
c......north
c     
            xobs(i) = slx(i)
            zobs(i) = slz(i) + delz
            GO TO 90
c     
 52         CONTINUE
c.....south
c     
            xobs(i) = slx(i)
            zobs(i) = slz(i) - delz
            GO TO 90
c     
 53         CONTINUE
c.....east
c     
            xobs(i) = slx(i) + delx
            zobs(i) = slz(i)
            GO TO 90
c     
 54         CONTINUE
c.....west
c     
            xobs(i) = slx(i) - delx
            zobs(i) = slz(i)
c     
 90         CONTINUE

         END DO                 ! Loop on I, nzspt

         IF ( lfarw .gt. 0 ) GO TO 200
c     
c.... Integrate over Shell:

         isg = 1
         WRITE ( outmod,
     $        '("isgn for shell integ. contrib. = ", i3)') isg
c     
         DO i = 1, nzspt
            cwrkr(nsew,i) = 0.0
         END DO

         ns = mth
         CALL chic1 ( xwal,zwal,xwalp,zwalp,isg,hwr,nzspt,
     $        ns,0,cwrkr,nsew )

c..   Add the shell result.

         DO  i = 1, nzspt
            chir(nsew,i) = cwrkr(nsew,i)
         END DO
c     
         WRITE (outmod,'("Coil-Wall contrib.to probes",i4)') nsew
         DO  i = nzspt/2, nzspt/2+2
            WRITE(outmod,'("i,chir,cwrkr",i4,1p4e12.3)')
     $           i,chir(nsew,i),cwrkr(nsew,i)
         END DO
c     
 200     CONTINUE

      END DO                    ! Loop over NSEW
c     
c.....calculate bx and bz
c     
      ALLOCATE ( bxr(nzspt), bzr(nzspt) )
c     
      DO i = 1, nzspt
         bxr(i) = 0.0
         bzr(i) = 0.0
      END DO

      DO i = 1, nzspt
c     
         bxr(i) = ( chir(3,i) - chir(4,i) ) / (2.0*delx)
         bzr(i) = ( chir(1,i) - chir(2,i) ) / (2.0*delz)
         WRITE ( outmod,'("bxr,zr = ", i3, 1p4e12.3)')
     $        i, bxr(i),bzr(i)
      END DO

      zmu0 =  4.0 * pye * 1.0e-7
      zmu0f = zmu0

      bxr = zmu0f * bxr
      bzr = zmu0f * bzr
      bpi = zmu0f * bpi
      
      CALL rflux ( slx,slxt, slzt, bxr,bzr, nzspt, dtsl, flxr )
      CALL rfluxp ( slx,slxt, slzt, bxr,bzr, nzspt, dtsl, pflxr )

      WRITE ( 6, '( "XXXXX  SENSORCW nzspt = ", i5 )' ) nzspt
      WRITE ( 23, '( "XXXXX  SENSORCW nzspt = ", i5 )' ) nzspt

c$$$      DEALLOCATE ( slx, slxt, slz, slzt )
      DEALLOCATE ( bxr, bzr )
      DEALLOCATE ( chir, cwrkr )

 500  CONTINUE

      END SUBROUTINE sensorcw

c......................................................................
      SUBROUTINE chic1(xsce,zsce,xscp,zscp,isg,creal,nobs,
     $     ns,ip,chir,nsew )
c......................................................................
c     
      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum5.inc'
c     
      DIMENSION iop(2), ww1(nths),ww2(nths),ww3(nths),tab(3)
      DIMENSION xsce(*),zsce(*),xscp(*),zscp(*)
      DIMENSION creal(*)
      DIMENSION chir(4,*)
c     
      COMMON / gridc / xobs(ndimlpc),zobs(ndimlpc)
c     
      dtpw = twopi / ns
c     
      ns1 = ns + 1
c     nobs = nloop + 3*nloopr

      DO io = 1, nobs
c     
         xs = xobs(io)
         zs = zobs(io)
c     
         DO is = 1, ns

            xt = xsce(is)
            zt = zsce(is)
            xtp = xscp(is)
            ztp = zscp(is)
c     
            CALL green

            chir(nsew,io) = chir(nsew,io) + aval * creal(is) 
c     
         END DO                 ! Loop over sources, IS: 1, ns
c     
         chir(nsew,io) = 0.5 * isg*dtpw * chir(nsew,io)
c     
      END DO                    ! Loop over observer Is: 1, nobs

      RETURN
      END SUBROUTINE chic1
c......................................................................
      SUBROUTINE chic2(xsce,zsce,xscp,zscp,isg,creal,nobs,
     $     zlens,ns,ip,chir,nsew )
c......................................................................
c.. zlens is the length of the source. ns is the open number of 
c    source points.     

      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum5.inc'
c     
      DIMENSION iop(2), ww1(nths),ww2(nths),ww3(nths),tab(3)
      DIMENSION xsce(*),zsce(*),xscp(*),zscp(*)
      DIMENSION creal(*)
      DIMENSION chir(4,*)
c     
      COMMON / gridc / xobs(ndimlpc),zobs(ndimlpc)
c     
      dtpw = zlens / ns
c     
      ns1 = ns + 1
c     nobs = nloop + 3*nloopr

      DO io = 1, nobs
c     
         xs = xobs(io)
         zs = zobs(io)
c     
         DO is = 1, ns
            is1 = is + 1
            xt = 0.5 * ( xsce(is) + xsce(is1) )
            zt = 0.5 * ( zsce(is) + zsce(is1) ) 
            xtp = 0.5 * ( xscp(is) + xscp(is1) ) 
            ztp = 0.5 * ( zscp(is) + zscp(is1) ) 
c     
            CALL green

            chir(nsew,io) = chir(nsew,io) + aval * creal(is) 
c     
         END DO                 ! Loop over source, IS: 1, ns
c     
         chir(nsew,io) = 0.5 * isg*dtpw * chir(nsew,io)
c     
      END DO                    ! Loop over observer, Io: 1, nobs

      RETURN
      END SUBROUTINE chic2

c........................................................
      SUBROUTINE sloops ( slx, slxt, slz, slzt, nzspt)
c........................................................

c... x, z, and their derivatives on a grid over the sensor loop face.
c    nzspt grid points in theta.

      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum5.inc'
      INCLUDE 'vacuum8.inc'
      INCLUDE 'vacuum7.inc'
c     
      REAL, DIMENSION(*) :: slx, slxt, slz, slzt
c     
      DIMENSION xpla(nths), zpla(nths)
c     
      EQUIVALENCE (xpla,xinf), (zpla,zinf)

      call boundsi ( xpla,zpla, 1,mth, xmnp,xmxp, zmnp,zmxp,
     $     ixn,ixx,izn,izx )
c     
      plrad = 0.5 * ( xmxp-xmnp )
      xmaj = 0.5 * ( xmxp+xmnp )
c     
      WRITE (OUTMOD,'("isloop",i5)') isloop
      GO TO ( 1,2,3,4 ), isloop
c     
    1 CONTINUE
c     
c.....nzspt points at constant distance from plasma.
c     use zork1, zork2 as storage
c     
c     .. Equal arcs on plasma first:
c     
      CALL eqarcw (xinf,zinf, wxgrd,wzgrd, wkgrd,wkth,workl,mth1 )
      wxgrd(mth2) = wxgrd(2)
      wzgrd(mth2) = wzgrd(2)
      DO i = 2, mth1
         alph = atan2m ( wxgrd(i+1)-wxgrd(i-1), wzgrd(i-1)-wzgrd(i+1) )
         zork1(i) = wxgrd(i) + asloop*plrad * cos(alph)
         zork2(i) = wzgrd(i) + asloop*plrad * sin(alph)
         WRITE (OUTMOD,'("zork1,zork2",i4,1p2e12.3)')
      END DO
c     
      zork1(1) = zork1(mth1)
      zork2(1) = zork2(mth1)
      zork1(mth2) = zork1(2)
      zork2(mth2) = zork2(2)
c     
      GO TO 150
c     
 2    CONTINUE

 3    CONTINUE

    4 CONTINUE

      IF ( lfarw > 0 ) THEN
         WRITE ( IOTTY, '(5x, "isloop = 4 inconsistent with no wall.
     $        WILL ABORT!!")' )
         WRITE ( OUTMOD, '(5x, "isloop = 4 inconsistent with no wall.
     $        WILL ABORT!!")' )
         CALL EXIT (0)
      END IF
c     
c.....nzspt points at constant distance from wall
c     use zork1, zork2 as storage
c     
      xwal(mth2) = xwal(2)
      zwal(mth2) = zwal(2)
      DO i = 2, mth1
         alph = atan2m ( xwal(i+1)-xwal(i-1), zwal(i-1)-zwal(i+1) )
         zork1(i) = xwal(i) - asloop*plrad * cos(alph)
         zork2(i) = zwal(i) - asloop*plrad * sin(alph)
c$$$  WRITE (OUTMOD,'("zork1,zork2",i4,1p2e12.3)')
c$$$  $      i,zork1(i),zork2(i)
      END DO
c     
      zork1(1) = zork1(mth1)
      zork2(1) = zork2(mth1)
      zork1(mth2) = zork1(2)
      zork2(mth2) = zork2(2)
      go to 150
c     
 150  CONTINUE
c     
c..   Now interpolate for NZSPT points in the range [SANG1,SANG2] in pi's
c     Need X, X_theta, Z_theta. PIWGRD is independent variable.
c     Store in SLX(i), SLXT(i), SLZ(i), SLZT(i).
c     
      zsang1 = sang1 * pye
      zsang2 = sang2 * pye
      dzsang = (zsang2 - zsang1) / (nzspt-1)
      dtsl = dzsang

      WRITE ( outmod, '()' )
      DO i = 1, nzspt
         zsang = zsang1 + (i-1)*dzsang
         if (zsang.lt.0.) zsang=zsang+2.*pye
         CALL lag (piwgrd,zork1,mth1, 5, zsang, slx(i), slxt(i), 2 )
         CALL lag (piwgrd,zork2,mth1, 5, zsang, slz(i), slzt(i), 2 )
         WRITE (OUTMOD,'("slx,slz",i4,1p2e12.3)')i,slx(i),slz(i)
      END DO
      call flush(OUTMOD)

      if ( lsensplot .eq. 1 ) then
         lsensplot=0
         call drawc4(xwal,zwal,xinf,zinf,xcw,zcw,slx,slz,
     $        1,mth1,1,mtcoil1,1,nzspt,
     $        "z","x",xmx,xma,zma,
     $        plrad,1,1,idot, jobid, ishape,a,b,aw,bw,cw,dw,tw,
     $        abulg, bbulg, tbulg, wcentr, pcircm, wcircm,
     $        pelong, pdness, wmaj,wlrad, welong, wdness, chgt )

         CALL drawc5(xwal,zwal,xinf,zinf,xcw,zcw,xcwa,zcwa,slx,slz,
     $        1,mth1,1,mtcoil1,1,mtcoila1,1,nzspt,
     $        "z","x",xmx,xma,zma,
     $        plrad,1,1,idot, jobid, ishape,a,b,aw,bw,cw,dw,tw,
     $        abulg, bbulg, tbulg, wcentr, pcircm, wcircm,
     $        pelong, pdness, wmaj,wlrad, welong, wdness,
     $        chgt,chgta,chgtb )

      end if


      END SUBROUTINE sloops

c.............................................................
      SUBROUTINE rflux ( zx,zxp, zzp, zbx,zbz, npts,zdt, rflx )
c.............................................................

      REAL, DIMENSION(*) :: zx,zxp, zzp, zbx,zbz
      REAL, DIMENSION(:), ALLOCATABLE:: zintl,zintd

      ALLOCATE ( zintl(npts), zintd(npts) )

      DO i = 1, npts
         zintd(i) = zx(i) * ( zxp(i)*zbz(i) - zzp(i)*zbx(i) )
      END DO

      CALL indef4 ( zintd, zintl, zdt, 1, npts, rflx, 1 )

      WRITE ( 6, '( "XXXXX  RFLUX npts = ", i5 )' ) npts
      WRITE ( 23, '( "XXXXX  RFLUX npts = ", i5 )' ) npts

      DEALLOCATE ( zintl ,zintd)

      END SUBROUTINE rflux

c.........................................................
      SUBROUTINE rfluxp ( zx,zxp, zzp, zbx,zbz, npts,zdt, rflx )
c.............................................................

      REAL, DIMENSION(*) :: zx,zxp, zzp, zbx,zbz
      REAL, DIMENSION(:), ALLOCATABLE:: zintl,zintd

      ALLOCATE ( zintl(npts), zintd(npts) )

      DO i = 1, npts
         zintd(i) = zx(i) * ( - zzp(i)*zbz(i) - zxp(i)*zbx(i) )
      END DO

      CALL indef4 ( zintd, zintl, zdt, 1, npts, rflx, 1 )

      DEALLOCATE ( zintl ,zintd)

      END SUBROUTINE rfluxp

