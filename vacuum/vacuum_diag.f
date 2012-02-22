c
c........................................................
      SUBROUTINE makebp ( inlabel )
c........................................................
c     
      include 'vacuum1.inc'
      include 'vacuum2.inc'
      include 'vacuum3.inc'
      include 'vacuum5.inc'
c
      character*(*) inlabel
      character*(40) pllabel
c 
      dimension chlagdy(nths,nfm)
      dimension thmgr(nths), z1tmp(nths), z2tmp(nths)

      COMPLEX, DIMENSION(:), ALLOCATABLE :: vecin6
c     
      write ( outmod, '(//, "Reads xil from: ", a )' ) inlabel
      write ( iotty , '(//, "Reads xil from: ", a )' ) inlabel
c     
      psipr = 1.0
      do  l = 1, lrnge
         xirlnq = ( lfm(l)-n*qa1 ) * xilr(l)
         xiilnq = ( lfm(l)-n*qa1 ) * xili(l)
         bnlr(l) = - psipr * xiilnq
         bnli(l) = + psipr * xirlnq
      end do

      zomsq = 0.0
      write ( outmod,222 ) zomsq, ( l,lfm(l),xilr(l),bnli(l),
     $     xili(l),bnlr(l), l = 1,lrnge )
 222  format ( //,4x, "zomsq=",e12.5,/,
     .     4x, "i",3x,"l",3x,"xilr-edge",7x,"bnli",4x,
     $     "xili-edge",8x,"-bnlr",
     .     /, 101(1x,2i4,1p4e12.4,/),/ )
c     
c... Construct xi(theta):
c
      do i = 1, mth1
         zgr = 0.0
         zgi = 0.0
         do l = 1, lrnge
            zgr = zgr + xilr(l)*cslth(i,l) - xili(l)*snlth(i,l)
            zgi = zgi + xilr(l)*snlth(i,l) + xili(l)*cslth(i,l)
         end do
         xigr(i) = zgr
         xigi(i) = zgi
      end do
c     
c... Construct bnk(theta):
c
      do i = 1, mth1
         zgr = 0.0
         zgi = 0.0
         do l = 1, lrnge
            zgr = zgr + bnlr(l)*cslth(i,l) - bnli(l)*snlth(i,l)
            zgi = zgi + bnlr(l)*snlth(i,l) + bnli(l)*cslth(i,l)
         end do
         bnkr(i) = zgr
         bnki(i) = zgi
      end do
c
      if ( lplot2 .eq. 0 ) go to 50
c
c....Plot bnkr  and bnki:
c
      call shftpi ( pigrd,xigr, z1tmp,z2tmp, mth1 )
      pllabel = inlabel//": xigr"
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 1,
     $     pllabel,mth )
      call shftpi ( pigrd,xigi, z1tmp,z2tmp, mth1 )
      pllabel = inlabel//": xigi"
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 2,
     $     pllabel,mth )
c
      call shftpi ( pigrd,bnkr, z1tmp,z2tmp, mth1 )
      pllabel = inlabel//": bnkr"
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 3,
     $     pllabel,mth )
      call shftpi ( pigrd,bnki, z1tmp,z2tmp, mth1 )
      pllabel = inlabel//": bnki"
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 4,
     $     pllabel,mth )
      CALL framep (jobid, ff)
c
 50   continue
c
c...Fourier analyse Xigr, xigi:
c
      l11 = lmin(1)
      l22 = lmax(1)
      call fanal ( xigr, mth, xirc, xirs, l11,l22, pye,0.0 )
      call fanal ( xigi, mth, xiic, xiis, l11,l22, pye,0.0 )
c     
      return
      end
c
c........................................................
      SUBROUTINE make_bltobp ( inlabel )
c........................................................
c     
      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum5.inc'
c
      CHARACTER*(*) inlabel
      CHARACTER*(40) pllabel
c 
      DIMENSION chlagdy(nths,nfm)
      DIMENSION thmgr(nths), z1tmp(nths), z2tmp(nths)
c     
      WRITE ( OUTMOD, '(//, "Reads bnlr, bnli from: ", a )' ) inlabel
      WRITE ( IOTTY , '(//, "Reads bnlr, bnli from: ", a )' ) inlabel

      psipr = 1.0
      DO  l = 1, lrnge
         psilnq = psipr * ( lfm(l)-n*qa1 )
         xilr(l) =  bnli(l) / psilnq
         xili(l) = -bnlr(l) / psilnq
      END DO
c     
      zomsq = 0.0
      WRITE ( OUTMOD,222 ) zomsq, ( l,lfm(l),xilr(l),bnli(l),
     $     xili(l),bnlr(l), l = 1,lrnge )
 222  FORMAT ( //,4x, "zomsq=",e12.5,/,
     .     4x, "i",3x,"l",3x,"xilr-edge",7x,"bnli",4x,
     $     "xili-edge",8x,"-bnlr",
     .     /, 101(1x,2i4,1p4e12.4,/),/ )
c     
c... Construct xi(theta): Neglect nqdelta for test
c$$$c
c$$$      DO I = 1, mth1
c$$$         zgr = 0.0
c$$$         zgi = 0.0
c$$$         DO l = 1, lrnge
c$$$            zgr = zgr + xilr(l)*coslt(i,l) - xili(l)*sinlt(i,l)
c$$$            zgi = zgi + xilr(l)*sinlt(i,l) + xili(l)*coslt(i,l)
c$$$         END DO
c$$$         xigr(i) = zgr
c$$$         xigi(i) = zgi
c$$$      END DO
c$$$c     
c$$$c... Construct bnk(theta):
c$$$c
c$$$      DO i = 1, mth1
c$$$         zgr = 0.0
c$$$         zgi = 0.0
c$$$         DO l = 1, lrnge
c$$$            zgr = zgr + bnlr(l)*coslt(i,l) - bnli(l)*sinlt(i,l)
c$$$            zgi = zgi + bnlr(l)*sinlt(i,l) + bnli(l)*coslt(i,l)
c$$$         END DO
c$$$         bnkr(i) = zgr
c$$$         bnki(i) = zgi
c$$$      END DO
c... Construct xi(theta):
c
      DO I = 1, mth1
         zgr = 0.0
         zgi = 0.0
         DO l = 1, lrnge
            zgr = zgr + xilr(l)*cslth(i,l) - xili(l)*snlth(i,l)
            zgi = zgi + xilr(l)*snlth(i,l) + xili(l)*cslth(i,l)
         END DO
         xigr(i) = zgr
         xigi(i) = zgi
      END DO
c     
c... Construct bnk(theta):
c
      DO i = 1, mth1
         zgr = 0.0
         zgi = 0.0
         DO l = 1, lrnge
            zgr = zgr + bnlr(l)*cslth(i,l) - bnli(l)*snlth(i,l)
            zgi = zgi + bnlr(l)*snlth(i,l) + bnli(l)*cslth(i,l)
         END DO
         bnkr(i) = zgr
         bnki(i) = zgi
      END DO

c
      IF ( lplot2 == 0 ) GO TO 50
c
c....Plot xir, xii, bnkr and bnki as function of \theta
c
      CALL shftpi ( pigrd,xigr, z1tmp,z2tmp, mth1 )
      pllabel = inlabel//": xigr"
      CALL pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 1,
     $     pllabel,mth )
      CALL shftpi ( pigrd,xigi, z1tmp,z2tmp, mth1 )
      pllabel = inlabel//": xigi"
      CALL pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 2,
     $     pllabel,mth )
c
      CALL shftpi ( pigrd,bnkr, z1tmp,z2tmp, mth1 )
      pllabel = inlabel//": bnkr"
      CALL pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 3,
     $     pllabel,mth )
      CALL shftpi ( pigrd,bnki, z1tmp,z2tmp, mth1 )
      pllabel = inlabel//": bnki"
      CALL pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 4,
     $     pllabel,mth )

      CALL framep(jobid, ff)
c
 50   CONTINUE
c
c...Fourier analyse Xigr, xigi:
c
      l11 = lmin(1)
      l22 = lmax(1)
      CALL fanal ( xigr, mth, xirc, xirs, l11,l22, pye,0.0 )
      CALL fanal ( xigi, mth, xiic, xiis, l11,l22, pye,0.0 )
c     
      return
      end

c........................................................
      SUBROUTINE makebl0 ( inlabel )
c........................................................
c     

      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum5.inc'
c
      CHARACTER*(*) inlabel
      CHARACTER*(40) pllabel
c
      DIMENSION thmgr(nths), z1tmp(nths), z2tmp(nths)
c
      WRITE ( OUTMOD, '(//, "Reads B from: ", a )' ) inlabel
      WRITE ( IOTTY,  '(//, "Reads B from: ", a )' ) inlabel
c
c..Multiply B by Jacobian*abs(grad-psi):
      DO i = 1, mth1
         bnkr(i) = gpsjp(i) * bnkr(i)
         bnki(i) = gpsjp(i) * bnki(i)
      END DO

c...  Plot real and imaginary parts of B(theta). Outer side is centered.
c...  Shift by mfel/2:

      IF ( lplot2 == 0 ) GO TO 25
c
c....Plot bnkr  and bnki:
c
      CALL shftpi ( pigrd,bnkr, z1tmp,z2tmp, mth1 )
      pllabel = inlabel//": bnkr"
      CALL pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 1,
     $     pllabel,mth )
      CALL shftpi ( pigrd,bnki, z1tmp,z2tmp, mth1 )
      pllabel = inlabel//": bnki"
      CALL pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 2,
     $     pllabel,mth )

 25   CONTINUE

      CALL vecwrt ( mth, bnkr, "bnkr", 1, mth1, outmod,0 )
      CALL vecwrt ( mth, bnki, "bnki", 1, mth1, outmod,0 )

c...  Fourier analyse the magnetic field  Put in xirc and xiic, etc.

      l11 = lmin(1)
      l22 = lmax(1)
      CALL fanal ( bnkr, mth, xirc, xirs, l11,l22, pye,-0.0 )
      CALL fanal ( bnki, mth, xiic, xiis, l11,l22, pye,-0.0 )

      DO j1 = 1, lrnge
         bnlr(j1) =   xirc(j1) + xiis(j1)
         bnli(j1) =   xiic(j1) - xirs(j1)
      END DO

      zomsq = 0.0
      WRITE ( outmod,222 ) zomsq, ( l,lfm(l),bnlr(l),bnli(l),
     $     xirc(l),xirs(l), xiic(l),xiis(l),
     $     l = 1,lrnge )
 222  FORMAT ( //,4x, "zomsq=",e12.5,/,
     .     4x, "i",3x,"l",3x,"blnr-edge",1x,"bnli-edge", 2x,
     $     "xirc,    xirs,     xiic,    xiis",
     .     /, 101(1x,2i4,1p6e12.4,/),/ )

      IF ( lplot2 == 0 ) GO TO 50

c$$$      CALL shftpi ( thmgr,bnkr, z1tmp,z2tmp, mfel1 )
c$$$      pllabel = inlabel//": xigr on mfel"
c$$$      CALL pospl2 ( z1tmp,z2tmp,zerov, 1,mfel1, 1,
c$$$     $     pllabel,mfel )
c$$$      CALL shftpi ( thmgr,bnki, z1tmp,z2tmp, mfel1 )
c$$$      pllabel = inlabel//": xigi on mfel"
c$$$      CALL pospl2 ( z1tmp,z2tmp,zerov, 1,mfel1, 2,
c$$$     $     pllabel,mfel )
c
      CALL framep(jobid, ff)

 50   CONTINUE

      RETURN
      END

c........................................................
      SUBROUTINE makebl ( inlabel )
c........................................................
c     

      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum5.inc'
c
      character*(*) inlabel
      character*(40) pllabel
c
      dimension thmgr(nths), z1tmp(nths), z2tmp(nths)
c
      WRITE ( outmod, '(//, "Reads B from: ", a )' ) inlabel
      WRITE ( iotty,  '(//, "Reads B from: ", a )' ) inlabel
c
c..Multiply Alan's B by Jacobian*abs(grad-psi):
      do i = 1, mth1
         bnkr(i) = gpsjp(i) * bnkr(i)
         bnki(i) = gpsjp(i) * bnki(i)
      end do

c     
c...  Plot real and imaginary parts of B(theta). Outer side is centered.
c...  Shift by mfel/2:
c     
      WRITE ( 6, '(/,/,"dmfel, mfel1 = ", 1pe12.4, i4,/,/)' )
     $     dmfel, mfel1
      do i = 1, mfel1
         thmgr(i) = (i-1) * dmfel
      end do
c
      if ( lplot2 .eq. 0 ) go to 25
c
c....Plot bnkr  and bnki:
c
      call shftpi ( pigrd,bnkr, z1tmp,z2tmp, mth1 )
      pllabel = inlabel//": bnkr"
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 1,
     $     pllabel,mth )
      call shftpi ( pigrd,bnki, z1tmp,z2tmp, mth1 )
      pllabel = inlabel//": bnki"
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 2,
     $     pllabel,mth )
c
 25   continue
c     
      call vecwrt ( mth, bnkr, "bnkr", 1, mth1, outmod,16 )
      call vecwrt ( mth, bnki, "bnki", 1, mth1, outmod,16 )
c     
c...  Fourier analyse Gato's magnetic field  Put in xirc and xiic, etc.
c
      l11 = lmin(1)
      l22 = lmax(1)
      call fanal ( bnkr, mth, xirc, xirs, l11,l22, pye,-0.0 )
      call fanal ( bnki, mth, xiic, xiis, l11,l22, pye,-0.0 )
c     
c... Since ist should be exp - il\theta, need to change sign on sines:
c      Changed 10/2004

c      xirs = - xirs
c      xiis = - xiis 

      do j1 = 1, lrnge
c$$$  xi(j1)  = xirc(j1) + xiis(j1)
c$$$  xii(j1) = xiic(j1) - xirs(j1)
c.... Since GATO uses exp(in\phi), take comp. conj of B:
c     so xii -> - xii, xir -> + xir
c$$$         bnlr(j1) =   xirc(j1) - xiis(j1)
c$$$         bnli(j1) = - xiic(j1) - xirs(j1)
c... Test other combination for SVD:
         bnlr(j1) =   xirc(j1) + xiis(j1)
         bnli(j1) =   xiic(j1) - xirs(j1)
      end do
c     
      zomsq = 0.0
      WRITE ( outmod,222 ) zomsq, ( l,lfm(l),bnlr(l),bnli(l),
     $     xirc(l),xirs(l), xiic(l),xiis(l),
     $     l = 1,lrnge )
 222  format ( //,4x, "zomsq=",e12.5,/,
     .     4x, "i",3x,"l",3x,"blnr-edge",1x,"bnli-edge", 2x,
     $     "xirc,    xirs,     xiic,    xiis",
     .     /, 101(1x,2i4,1p6e12.4,/),/ )
c     
c
      if ( lplot2 .eq. 0 ) go to 50
c
      call shftpi ( pigrd,xigr, z1tmp,z2tmp, mth1 )
      pllabel = inlabel//": xigr"
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 3,
     $     pllabel,mth )
      call shftpi ( pigrd,xigi, z1tmp,z2tmp, mth1 )
      pllabel = inlabel//": xigi"
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 4,
     $     pllabel,mth )
c
c$$$      call shftpi ( thmgr,bnkr, z1tmp,z2tmp, mfel1 )
c$$$      pllabel = inlabel//": xigr on mfel"
c$$$      call pospl2 ( z1tmp,z2tmp,zerov, 1,mfel1, 1,
c$$$     $     pllabel,mfel )
c$$$      call shftpi ( thmgr,bnki, z1tmp,z2tmp, mfel1 )
c$$$      pllabel = inlabel//": xigi on mfel"
c$$$      call pospl2 ( z1tmp,z2tmp,zerov, 1,mfel1, 2,
c$$$     $     pllabel,mfel )
c
      CALL framep(jobid, ff)
c
 50   continue
c     
      return
      end
c     
c........................................................
      subroutine makebg ( inlabel )
c........................................................
c     

      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum5.inc'
c
      character*(*) inlabel
      character*(40) pllabel
c
      dimension thmgr(nths), z1tmp(nths), z2tmp(nths)
c
      WRITE ( outmod, '(//, "Reads xi from: ", a )' ) inlabel
      WRITE ( iotty,  '(//, "Reads xi from: ", a )' ) inlabel
c     
c...  Plot real and imaginary parts of Xi(theta). Outer side is centered.
c...  Shift by mfel/2:
c     
      WRITE ( 6, '(/,/,"dmfel, mfel1 = ", 1pe12.4, i4,/,/)' )
     $     dmfel, mfel1
      do i = 1, mfel1
         thmgr(i) = (i-1) * dmfel
      end do
c
      if ( lplot2 .eq. 0 ) go to 25
c
      call shftpi ( thmgr,xigr, z1tmp,z2tmp, mfel1 )
      pllabel = inlabel//": xigr on mfel"
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mfel1, 1,
     $     pllabel,mfel )
      call shftpi ( thmgr,xigi, z1tmp,z2tmp, mfel1 )
      pllabel = inlabel//": xigi on mfel"
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mfel1, 2,
     $     pllabel,mfel )
c
 25   continue
c     
      call vecwrt ( mfel, xigr, "xigr", 1, mfel1, outmod,16 )
      call vecwrt ( mfel, xigi, "xigi", 1, mfel1, outmod,16 )
c     
c...  Fourier analyse Gato's displacement. Put in xirc and xiic, etc.
c
      l11 = lmin(1)
      l22 = lmax(1)
      call fanal ( xigr, mfel, xirc, xirs, l11,l22, pye,-0.5 )
      call fanal ( xigi, mfel, xiic, xiis, l11,l22, pye,-0.5 )
c     
      do j1 = 1, lrnge
c$$$  xi(j1)  = xirc(j1) + xiis(j1)
c$$$  xii(j1) = xiic(j1) - xirs(j1)
c.... Since GATO uses exp(in\phi), take comp. conj of xi:
c     so xii -> - xii, xir -> + xir
         xilr(j1) =   xirc(j1) - xiis(j1)
         xili(j1) = - xiic(j1) - xirs(j1)
      end do
c     
      psipr = 1.0
      do  l = 1, lrnge
         xirlnq = ( lfm(l)-n*qa1 ) * xilr(l)
         xiilnq = ( lfm(l)-n*qa1 ) * xili(l)
         bnlr(l) = - psipr * xiilnq
         bnli(l) = + psipr * xirlnq
      end do
c     
      zomsq = 0.0
      WRITE ( outmod,222 ) zomsq, ( l,lfm(l),xilr(l),bnli(l),
     $     xili(l),bnlr(l), l = 1,lrnge )
 222  format ( //,4x, "zomsq=",e12.5,/,
     .     4x, "i",3x,"l",3x,"xilr-edge",7x,"bnli",4x,
     $     "xili-edge",8x,"-bnlr",
     .     /, 101(1x,2i4,1p4e12.4,/),/ )
c     
c...  Construct xi(theta) on mth grid.
c     
      do i = 1, mth1
         zgr = 0.0
         zgi = 0.0
         do l = 1, lrnge
            zgr = zgr + xilr(l)*cslth(i,l) - xili(l)*snlth(i,l)
            zgi = zgi + xilr(l)*snlth(i,l) + xili(l)*cslth(i,l)
         end do
         xigr(i) = zgr
         xigi(i) = zgi
      end do
c     
c...  Construct bnk(theta):
c     
      do i = 1, mth1
         zgr = 0.0
         zgi = 0.0
         do l = 1, lrnge
            zgr = zgr + bnlr(l)*cslth(i,l) - bnli(l)*snlth(i,l)
            zgi = zgi + bnlr(l)*snlth(i,l) + bnli(l)*cslth(i,l)
         end do
         bnkr(i) = zgr
         bnki(i) = zgi
      end do
c
      if ( lplot2 .eq. 0 ) go to 50
c
c....Plot bnkr  and bnki:
c
      call shftpi ( pigrd,bnkr, z1tmp,z2tmp, mth1 )
      pllabel = inlabel//": bnkr"
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 3,
     $     pllabel,mth )
      call shftpi ( pigrd,bnki, z1tmp,z2tmp, mth1 )
      pllabel = inlabel//": bnki"
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 4,
     $     pllabel,mth )
      CALL framep(jobid, ff)
c
 50   continue
c     
      return
      end
c     
c..........................................................
      subroutine diaplt
c..........................................................
c     
      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum5.inc'
      INCLUDE 'vacuum6.inc'
      INCLUDE 'vacuum8.inc'
c     
      dimension z1tmp(nths), z2tmp(nths), zorkr(nths),zorki(nths),
     $     zorkpr(nths), zorkpi(nths), chlagdy(nths,nfm),
     $     thph(nths), cppgr(nths),cppgi(nths), 
     $     cplgr(nths), cplgi(nths), cplgtr(nths), cplgti(nths),
     $     chwr1(nths),chwi1(nths),
     $     dxdt(nths), dzdt(nths), zkt(nths,2), zkp(nths,2)
c     
      call atpoint("DIAPLT","lrnge",lrnge,"bnli(1)",bnli(1),
     $     iotty, outmod )

      jmax1 = lrnge
      laxis = lrnge

      zorkpr = 0.0
      zorkpi = 0.0

c$$$      IF ( laxis == 1 ) laxis =2   !  COMMENTED OUT 09-30-2008
c
c.. New variable, lfarw for SPARK2 to reset farwal.
c
c$$$      lfarw = 0
c$$$      if ( a .gt. 10.0 ) lfarw = 1
c$$$      WRITE ( iotty, '("a = ", 1pe12.3, " lfarw in DIAPLT = ", i4 )' )
c$$$     $     a, lfarw
c$$$      WRITE ( outmod, '("a = ", 1pe12.3, " lfarw in DIAPLT = ", i4 )' )
c$$$     $     a, lfarw

      zomsq = 0.0
      WRITE ( outmod,222 ) zomsq, ( l,lfm(l),xilr(l),bnli(l),
     $     xili(l),bnlr(l), l = 1,lrnge )
 222  format ( //,4x, "zomsq=",e12.5,/,
     .     4x, "i",3x,"l",3x,"xilr-edge",7x,"bnli",4x,
     $     "xili-edge",8x,"-bnlr",
     .     /, 101(1x,2i4,1p4e12.4,/),/ )
c     
      if ( lplot2 .eq. 0 ) go to 100
c
c..   Plot xilr, xili, bnlr, and bnli versus l:
c     
      call pospl2 ( xlfm,xilr, zerov, 1, laxis, 1,'Xilr vs. l', laxis )
      call pospl2 ( xlfm,xili, zerov, 1, laxis, 2,'Xili vs. l', laxis )
      call pospl2 ( xlfm,bnlr, zerov, 1, laxis, 3,'bnlr vs. l', laxis )
      call pospl2 ( xlfm,bnli, zerov, 1, laxis, 4,'bnli vs. l', laxis )
      call framep(jobid,ff)
c
 100  continue
c     
      l11 = lmin(1)
c     
      WRITE ( outmod,'(/,5x,"FOURIER COEFF. OF XIR AND XII",/,
     $     3x,"L", 5x, "XIRC",7x, "XIRS",7x, "XIIC",7x, "XIIS",/,
     $     (i4, 1p4e11.3))' )
     $     (l11-1+i, xirc(i),xirs(i),xiic(i),xiis(i),i=1, lrnge)
c     
c...  Shift by mth/2:
c     
c
      if ( lplot2 .eq. 0 ) go to 200
c
      call shftpi ( pigrd,xigr, z1tmp,z2tmp, mth1 )
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 1,"xigr:mth",mth )
      call shftpi ( pigrd,xigi, z1tmp,z2tmp, mth1 )
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 2,"xigi:mth",mth )
      call framep(jobid,ff)
c     
      call pospl2 ( xlfm,xirc,zerov, 1,laxis, 1,'Xirc vs. m', laxis )
      call pospl2 ( xlfm,xirs,zerov, 1,laxis, 2,'Xirs vs. m', laxis )
      call pospl2 ( xlfm,xiic,zerov, 1,laxis, 3,'Xiic vs. m', laxis )
      call pospl2 ( xlfm,xiis,zerov, 1,laxis, 4,'Xiis vs. m', laxis )
      call framep(jobid,ff)
c     
c...  Plot Xi/grad-psi, and Bn/Jacob-grad-psi
c     
      do i = 1, mth1
         zorkr(i) = xigr(i) / sqgrps(i)
         zorki(i) = xigi(i) / sqgrps(i)
      end do
c     
      call shftpi ( pigrd,zorkr, z1tmp,z2tmp, mth1 )
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 1,"Xir/sqrt",mth )
      call shftpi ( pigrd,zorki, z1tmp,z2tmp, mth1 )
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 2,"Xii/sqrt",mth )
c     
      do i = 1, mth1
         zorkr(i) = bnkr(i) / gpsjp(i)
         zorki(i) = bnki(i) / gpsjp(i)
      end do
c     
      call shftpi ( pigrd,zorkr, z1tmp,z2tmp, mth1 )
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 3,"Bnkr/gpsjp",mth )
      call shftpi ( pigrd,zorki, z1tmp,z2tmp, mth1 )
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 4,"Bnki/gpsjp",mth )
      call framep(jobid,ff)
c
 200  continue
c     
c     ... Calculate the normal magnetic field at the surface. 
c     Note that a factor of at least f/R is dropped from PEST-I stuff.
c     
      do  i = 1, mth1
         bnpptr(i) = 0.0
         bnppti(i) = 0.0
         do j1 = 1, lrnge
            bnpptr(i) = bnpptr(i) + ( bnlr(j1)*cslth(i,j1)
     $           - bnli(j1)*snlth(i,j1) ) / gpsjp(i)
            bnppti(i) = bnppti(i) + ( bnlr(j1)*snlth(i,j1)
     $           + bnli(j1)*cslth(i,j1) )/ gpsjp(i)
         end do
      end do
c     
c------------------------ Commented out 030109
c$$$      WRITE ( OUTMOD,
c$$$     $     '(/, "X, Z, Bf_real,imag ,gpsjp, Bf_imag, B_r, B_i:-", /, 
c$$$     $     (i4, 7es13.5) )' )
c$$$     $     (i, xinf(i),zinf(i), bnkr(i),bnki(i),gpsjp(i),
c$$$     $     bnpptr(i), bnppti(i), i = 1, mth1)

      if ( lplot2 .eq. 0 ) go to 300
c
      call shftpi ( pigrd,bnpptr, z1tmp,z2tmp, mth1 )
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 1,"Bnpptr",mth )
      call shftpi ( pigrd,bnppti, z1tmp,z2tmp, mth1 )
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 2,"Bnppti",mth )
      call shftpi ( pigrd,gpsjp, z1tmp,z2tmp, mth1 )
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 3,"GPSJP",mth )
      IF ( lfarw /= 1  ) THEN
         CALL shftpi ( pigrd,gpsjw, z1tmp,z2tmp, mth1 )
         CALL pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 4,"GPSJW",mth )
      END IF
      CALL FRAMEP(JOBID,FF)
c     
c...plot versus length.
c
      call shftpi ( slngth,bnpptr, z1tmp,z2tmp, mth1 )
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 1,"Bnpptr vs. lth",mth )
      call shftpi ( slngth,bnppti, z1tmp,z2tmp, mth1 )
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 2,"Bnppti vs. lth",mth )
      call shftpi ( slngth,gpsjp, z1tmp,z2tmp, mth1 )
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 3,"GPSJP vs. lth",mth )
      IF ( lfarw /= 1  ) THEN
         CALL shftpi ( slngth,gpsjw, z1tmp,z2tmp, mth1 )
         CALL pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 4,
     $        "GPSJW vs. lth",mth )
      END IF
      CALL FRAMEP(JOBID,FF)
c
 300  continue
c     
c.....plot chi on plasma and on conductor.
c     
c     
      do i = 1, mth
         chipr(i) = 0.0
         chipi(i) = 0.0
         do  j1 = 1, jmax1
            zrrc = cplar(i,j1)
            zrrs = cplai(i,j1)
            chipr(i) = chipr(i) + zrrc*bnlr(j1) - zrrs*bnli(j1)
            chipi(i) = chipi(i) + zrrs*bnlr(j1) + zrrc*bnli(j1)
         end do
      end do
      chipr(mth1) = chipr(1)
      chipi(mth1) = chipi(1)
      chipr(mth2) = chipr(2)
      chipi(mth2) = chipi(2)

c
c.....Calculate B theta on the plasma ...
c     
      DO i = 1, mth1
         CALL lagp ( slngth,chipr, mth1,5, slngth(i), f,zorkpr(i), 1,1 )
         CALL lagp ( slngth,chipi, mth1,5, slngth(i), f,zorkpi(i), 1,1 )
         zork3(i) = sqrt ( zorkpr(i)**2 + zorkpi(i)**2 )
      END DO
      bthpr = zorkpr
      bthpi = zorkpi

c calculate b_phi. Use z1tmp, z2tmp as storage.
      z1tmp(1:mth2) =  n*chipi(1:mth2)/xinf(1:mth2)
      z2tmp(1:mth2) = -n*chipr(1:mth2)/xinf(1:mth2)

c---------------Commented out 030109
c$$$      WRITE ( OUTMOD, '(/,40x,"Plasma Surface Fields")' )
c$$$      WRITE ( OUTMOD, '(1x, A4, 9A13)' )
c$$$     $     "i", "X", "Z" ,"gpsjp", "Bn_r","Bn_i", "Bth_r","Bht_i",
c$$$     $     "Bph_r","Bph_i"
c$$$      WRITE ( OUTMOD, '(1x, I4, 9ES13.5)' )
c$$$     $     (i, xinf(i),zinf(i),gpsjp(i), bnpptr(i),bnppti(i),
c$$$     $     bthpr(i),bthpi(i),  z1tmp(i),z2tmp(i), i = 1, mth1)
c$$$      WRITE ( OUTMOD, '(1x, A4, 9A13)' )
c$$$     $     "i", "X", "Z" ,"gpsjp", "Bn_r","Bn_i", "Bth_r","Bht_i",
c$$$     $     "Bph_r","Bph_i"

      IF ( lplot3 == 0 ) GO TO 400

c.... shift so that outside is centered.
c     
      call shftpi ( slngth,zorkpr, z1tmp,z2tmp, mth1 )
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 1,
     $     'B-th cos vs. lth', mth1 )
      call shftpi ( slngth,zorkpi, z1tmp,z2tmp, mth1 )
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 2,
     $     'B-th sin vs. lth', mth1 )
      call shftpi ( slngth,zork3, z1tmp,z2tmp, mth1 )
      call pospl2 ( z1tmp,z2tmp,zerov, 1,mth1, 3,
     $     'B-th ampl vs. lth', mth1 )
c     
c.....plot arc lengths vs.theta.
c     
      call pospl2 (pigrd,slngth,zerov,1,mth1,4,
     $     'lp vs.th,out=0',mth1)
c     
      call framep(jobid,ff)
c
c      if ( ipshp .ne. 1 ) then
         call shftpi ( slngth,chipr, z1tmp,zorkr, mth1 )
         call shftpi ( slngth,chipi, z1tmp,zorki, mth1 )
         call pospl2 ( z1tmp,zorkr,zorki,1,mth1, 1,
     $        'chi-plasma vs. l',jmax1)
c      end if
      call pospl1(chipr,chipi, 1,mth1, 2,'chpr-chpi',jmax1)
c     
      call framep(jobid,ff)
c
 400  continue
c     
      pi = pye
      tpn = 2.0 * pi * n
c
      if ( lfarw .eq. 1 ) go to 5000
c
      ngr = noutv
c     
c     ****************************** for outboard gap *****dylee,03/18/97
c     
      if ((ishape .eq. 340) .or. (ishape .eq. 240)) ngr=2*noutv-1
c     
c     ***********************************************************
c     
c...  if ( ishape .gt. 100 ) ngr = 2*noutv - 1
      ngr2 = (ngr-1)/2 + 1
      ngr4 = (ngr-1)/4 + 1
      ngr34 = 3*(ngr-1)/4 + 1
      f = 0.0
      df = 0.0
c     
      isph = 1
      if ( jtop .eq. jbot ) isph = 0
c     
      do i = 1, ngr
         thph(i) = ( 1.0/(ngr-1) ) * (i-1)
      end do
c     
c......Extrtact chi-l on the wall -- both the cos(lth) and
c     the sin(lth) contribution. Note that the index is shifted
c     here.
c     
      mw = mth - (jtop-jbot-1)
c     
      do l1 = 1, jmax1
         iw = mth - jtop + 2
         icnt = 0
         do  i = 1, mth
            if (  (i.gt.jbot) .and. (i.lt.jtop) ) go to 1200
            if ( i .eq. jtop ) iw = 1
            chiwc(iw,l1) = cwallr(i,l1)
            chiws(iw,l1) = cwalli(i,l1)
            iw = iw + 1
            icnt = icnt + 1
 1200       continue
         end do
c.....fill in last point if torus..
         if ( jtop .ne. jbot ) go to 1210
         chiwc(mw,l1) = chiwc(1,l1)
         chiws(mw,l1) = chiws(1,l1)
 1210    continue
      end do
c     
      if ( checkd ) then
c     
c.......plot chiwc and chiws for each l
c     
         do j1 = 1, jmax1
            ipos = mod(j1,4)
            if (ipos .eq. 0 ) ipos = 4
c     
            do iw = 1, mw
               zorkr(iw) = chiwc(iw,j1)*bnlr(j1)-chiws(iw,j1)*bnli(j1)
               zorki(iw) = chiws(iw,j1)*bnlr(j1)+chiwc(iw,j1)*bnli(j1)
c     chlagdy for later use to write out
               chlagdy(iw,j1) = zorkr(iw)
            end do
c     
c$$$  c     ******* chiwc,chiws multiplied by xilmnq,by dylee, 4/17/97 ******
c$$$  do 11 iw = 1, mw
c$$$  chlagc(iw) = chiwc(iw,j1)*xilmnq
c$$$  chlags(iw) = chiws(iw,j1)*xilmnq
c$$$  c     
c$$$  c     chlagdy for later use to write out
c$$$  chlagdy(iw,j1) = chiwc(iw,j1)*xilmnq      
c$$$  11         continue
c     *****************************************************************
c     
            call pospl2 (wgrd,zorkr,zorki,1,mw,ipos,'chi-l-wall',jmax1)
c     
c     ****  plot only one graph, real or imaginary part, dylee,4/17/97 *****
c     call pospl1(wgrd,zorkr, 1,mw, ipos,'chi-l-wall',jmax1)
c     
            if ( ipos .eq. 4 ) call framep(jobid,ff)
c     
         end do
c     
         if ( ipos .ne. 4 ) call framep(jobid,ff)
c     
c     ******* write out chlagc for several l values, for dylee, 4/18/97 ****
c     ** this should be commented out for usual ordinary runs **
c     do iw = 1, mw
c     write (outmod,2123) wgrd(iw), chlagdy(iw,6),chlagdy(iw,7),
c     $         chlagdy(iw,8),chlagdy(iw,9),chlagdy(iw,10),
c     $         chlagdy(iw,11),chlagdy(iw,12),chlagdy(iw,13)
c     end do
c     write (outmod, '( / )')
c     do iw = 1, mw
c     write (outmod,2123) wgrd(iw), chlagdy(iw,14),chlagdy(iw,15),
c     $         chlagdy(iw,16),chlagdy(iw,17),chlagdy(iw,18),
c     $         chlagdy(iw,19),chlagdy(iw,20),chlagdy(iw,21)
c     end do
c     2123 format (1x, f7.3, 1p8e11.3)
c     **********************************************************************
      end if
c     
c.....plot arc length vs.theta
c     
c     get min and max of wall coords.
c     
      zmax = -1.0
      xmax = -1.0
      do iw = 1, mw
c     ************************************************dylee,03/18/97
         if ( xpass(iw) .gt. xmax) then
            xmax = xpass(iw)
            if (zpass(iw) .ge. 0.) imwxmax = iw
         end if
c     ************************************************
c     if ( xpass(iw) .gt. xmax ) xmax = xpass(iw)
         if ( zpass(iw) .gt. zmax ) zmax = zpass(iw)
      end do
c     
      rmax = xmax
      if ( zmax .gt. xmax ) rmax = zmax
c     
c
      if ( lplot3 .eq. 0 ) go to 500
c
      call pospl2 ( wgrd,ell,zerov, 1,mw, 1, "lw vs. wallgr", mw )
      IF ( lfarw /= 1  ) THEN
         CALL shftpi ( pigrd,gpsjw, zork1,zork2, mth1 )
         CALL pospl2 ( zork1,zork2, zerov, 1,mth1, 2,'j*gpsi-w',mth1 )
      END IF
c
 500  continue
c     
      do  i = 1, mw
         chwr1(i) = 0.0
         chwi1(i) = 0.0
         do  j1 = 1, jmax1
            chwr1(i) = chwr1(i) + chiwc(i,j1)*bnlr(j1)
     $           - chiws(i,j1)*bnli(j1)
            chwi1(i) = chwi1(i) + chiws(i,j1)*bnlr(j1)
     $           + chiwc(i,j1)*bnli(j1)
         end do
      end do
c     
      if ( lplot3 .eq. 0 ) go to 600
c
      call pospl2(ell,chwr1,chwi1, 1,mw,3,'chi-wall',jmax1)
      call pospl1(chwr1,chwi1, 1,mw,4,'chwc-chws', jmax1)
c     
      call framep(jobid,ff)
c
 600  continue
c     
c     ******* write out chlagc for several l values, for dylee, 4/21/97 ****
c     ** this should be commented out for usual ordinary runs **
c     
      if (checkd) then
c     
         write (outmod,
     $        '(/,5x,"ell, chiwl, j1=6 to 14,l=j1+lmin-1",/)')
         do iw = 1, mw, 4
            write (outmod,2123) ell(iw),
     $           chlagdy(iw,6),chlagdy(iw,7),
     $           chlagdy(iw,8),chlagdy(iw,9),chlagdy(iw,10),
     $           chlagdy(iw,11),chlagdy(iw,12),
     $           chlagdy(iw,13),chlagdy(iw,14)
         end do
c     
         write (outmod, '(/,5x,"ell, chiwl, j1=15 to 21,l=j1+lmin-1,
     $        chiwc,chiws",/)')
         do iw = 1, mw, 4
            write (outmod,2123) ell(iw),
     $           chlagdy(iw,15),chlagdy(iw,16),
     $           chlagdy(iw,17),chlagdy(iw,18),chlagdy(iw,19),
     $           chlagdy(iw,20),chlagdy(iw,21),
     $           chwr1(iw),chwi1(iw)
         end do
 2123    format (1x, f7.3, 1p9e11.3)
c     
      end if
c     ***********************************************************************
c     
 5000 continue  ! .............  FARWAL ....................
c
      if ( lspark .eq. 0 ) then
c
c... Get plasma response
c
c     Sum over l with bnlr and bnli here plasma response..
c     Use cplgt for storage as in the spark2 case.
c     
      do i = 1, mth1
         cplgtr(i) = 0.0
         cplgti(i) = 0.0 
         do  j1 = 1, jmax1
            cplgtr(i) = cplgtr(i) + cplar(i,j1)*bnlr(j1)
     $           - cplai(i,j1)*bnli(j1)
            cplgti(i) = cplgti(i) + cplai(i,j1)*bnlr(j1)
     $           + cplar(i,j1)*bnli(j1)
         end do
      end do
c   
      nbngr = 129
c
c...  Shift length coordinate on plasma.
c     
         mthh1 = mth/2 + 1
         call shflen ( slngth, zork1, mth1, mthh1 )
c     
         if ( lplot4 .eq. 0 ) go to 700
c
         IF ( lfarw .NE. 1 ) THEN
c     Write out information for the theta and phi components 
c     of the skin current at the shell.
c     
            write ( outmod,'( /,
     $           10x,"Theta, phi skin current components at shell-- ")')
c     
c...  chwr1 and chwi1 already shifted.
c     
            call kcur0 ( xpass,zpass,ell,wgrd,mw,chwr1,chwi1,n,
     $           nbngr,wkgrd,zorkr,zorki,dxdt,dzdt,
     $           zkt,zkp,nths,0,"Tot.Shell",lrnge, outmod )
            
         END IF
c     
c...  Now the skin current on the plasma.
c     First, shift variables so that theta=0 is on the inside.
c     
         call shftgr ( cplgtr, mth )
         call shftgr ( cplgti, mth )
         call shfgrb ( xinf, zork3, mth )
         call shfgrb ( zinf, zork4, mth )
c     
         write ( outmod,'( /,
     $        10x,"Theta, phi skin current components at plasma-- ")')
c     
         call kcur0 ( zork3,zork4,zork1,pgrd,mth1,cplgtr,cplgti,n,
     $        nbngr,wkgrd,zorkr,zorki,dxdt,dzdt,
     $        zkt,zkp,nths,0,"Tot. Plas",lrnge, outmod )
c
 700     continue
c     
      end if
c.................... lspark = 0 ......................
c
      if ( lspark .eq. 2 ) then
c
c     Sum over l with bnlr and bnli  here for pp response..
c     
      do i = 1, mth1
         cppgr(i)  = 0.0
         cppgi(i)  = 0.0 
         do  j1 = 1, jmax1
            cppgr(i) = cppgr(i) + cppr(i,j1)*bnlr(j1)
     $           - cppi(i,j1)*bnli(j1)
            cppgi(i) = cppgi(i) + cppi(i,j1)*bnlr(j1)
     $           + cppr(i,j1)*bnli(j1)
         end do
      end do

c     Sum over l with bnlr and bnli here for pw and pw+pp response..
c     
      do i = 1, mth1
         cplgr(i) = 0.0
         cplgi(i) = 0.0 
         cplgtr(i) = 0.0
         cplgti(i) = 0.0 
         do  j1 = 1, jmax1
            cplgr(i) = cplgr(i) + cwptr(i,j1)*bnlr(j1)
     $           - cwpti(i,j1)*bnli(j1)
            cplgi(i) = cplgi(i) + cwpti(i,j1)*bnlr(j1)
     $           + cwptr(i,j1)*bnli(j1)
            cplgtr(i) = cplgtr(i) + cpwtr(i,j1)*bnlr(j1)
     $           - cpwti(i,j1)*bnli(j1)
            cplgti(i) = cplgti(i) + cpwti(i,j1)*bnlr(j1)
     $           + cpwtr(i,j1)*bnli(j1)
         end do
      end do
c     
c.....Calculate magnetic field components at virtual wall position.
c     Also at plasma position due to wall. Bug fix on 10/17/97
c     
         do  i = 1, mth1
            bnpwtr(i) = 0.0
            bnpwti(i) = 0.0
            bnwptr(i) = 0.0
            bnwpti(i) = 0.0
c     
c...  Get rid of jacob*grpsi in bperp
c     
            do  j1 = 1, jmax1
               bnpwtr(i) = bnpwtr(i) + bnpwr(i,j1)*bnlr(j1)
     $              - bnpwi(i,j1)*bnli(j1)
               bnpwti(i) = bnpwti(i) + bnpwi(i,j1)*bnlr(j1)
     $              + bnpwr(i,j1)*bnli(j1)
               bnwptr(i) = bnwptr(i) + bnwpr(i,j1)*bnlr(j1)
     $              - bnwpi(i,j1)*bnli(j1)
               bnwpti(i) = bnwpti(i) + bnwpi(i,j1)*bnlr(j1)
     $              + bnwpr(i,j1)*bnli(j1)
            end do
            bnpwtr(i) = bnpwtr(i) / gpsjw(i)
            bnpwti(i) = bnpwti(i) / gpsjw(i)
            bnwptr(i) = bnwptr(i) / gpsjp(i)
            bnwpti(i) = bnwpti(i) / gpsjp(i)
         end do
c     
         bnpwtr(mth1) = bnpwtr(1) 
         bnpwti(mth1) = bnpwti(1) 
         bnwptr(mth1) = bnwptr(1) 
         bnwpti(mth1) = bnwpti(1) 
c     
c.... Write out Bnorm on plasma, on wall due to plasma, 
c     and vice versa on equal arcs grid
c     
c.....Shift the grid first
c     
         call shftgr ( bnpptr, mth )
         call shftgr ( bnppti, mth )
         call shftgr ( bnpwtr, mth )
         call shftgr ( bnpwti, mth )
         call shftgr ( bnwptr, mth )
         call shftgr ( bnwpti, mth )
c     
         nbngr = 129
c     
         write ( outmod,'( /,
     $        10x,"Normal Magnetic field at virtual wall -- ")' )
c     
         call bnsrf ( xpass,zpass,ell,wgrd,mw,bnpwtr,bnpwti,
     $        nbngr,piwgrd,zorkr,zorki,
     $        z1tmp,z2tmp,zorkpr,zorkpi,outmod )
c     
c......Now get bnr, bni on plasma due to plasma and due to wall 
c     for the equal arcs grid.
c     Shift the grid so that outer side is centered.
c     Use zorks as storage.
c     
         call shftpi ( zork1, xinf, zorkr, zork3, mth1 )
         call shftpi ( zork1, zinf, zorkr, zork4, mth1 )
c     
c...  Shift length coordinate on plasma.
c     
         mthh1 = mth/2 + 1
         call shflen ( slngth, zork1, mth1, mthh1 )
c     
         write ( outmod,'(/,
     $        10x,"Normal B-field at plasma due to plasma -- ")')
c     
         call bnsrf ( zork3,zork4,zork1,pgrd,mth1,bnpptr,bnppti,
     $        nbngr,pigrd,zorkr,zorki,z1tmp,z2tmp,zorkpr,zorkpi,outmod )
c     
         write ( outmod,'(/,
     $        10x,"Normal B-field at plasma due to wall -- ")')
c     
         call bnsrf ( zork3,zork4,zork1,pgrd,mth1,bnwptr,bnwpti,
     $        nbngr,pigrd,zorkr,zorki,z1tmp,z2tmp,zorkpr,zorkpi,outmod )
c     
c
         if ( lplot4 .eq. 0 ) go to 800
c
c     Write out information for the theta and phi components 
c     of the skin current at the shell.
c     
         write ( outmod,'( /,
     $        10x,"Theta, phi skin current components at shell-- ")')
c
c... chwr1 and chwi1 already shifted.
c   
         call kcur0 ( xpass,zpass,ell,wgrd,mw,chwr1,chwi1,n,
     $        nbngr,wkgrd,zorkr,zorki,dxdt,dzdt,
     $        zkt,zkp,nths,0,"Shell",lrnge, outmod )
c     
c...  Now the effective skin current on the plasma due to the wall current.
c     First, shift variables so that theta=0 is on the inside.
c     
         call shftgr ( cppgr, mth )
         call shftgr ( cppgi, mth )
         call shftgr ( cplgr, mth )
         call shftgr ( cplgi, mth )
         call shftgr ( cplgtr, mth )
         call shftgr ( cplgti, mth )
         call shfgrb ( xinf, zork3, mth )
         call shfgrb ( zinf, zork4, mth )
c     
         write ( outmod,'( /,
     $        10x," plasma-plasma theta, phi skin current  -- ")')
c     
         call kcur0 ( zork3,zork4,zork1,pgrd,mth1,cppgr,cppgi,n,
     $        nbngr,wkgrd,zorkr,zorki,dxdt,dzdt,
     $        zkt,zkp,nths,0,"pla-pla",lrnge, outmod )
c     
         write ( outmod,'( /,
     $        10x,"plasma from wall theta, phi skin current   -- ")')
c     
         call kcur0 ( zork3,zork4,zork1,pgrd,mth1,cplgr,cplgi,n,
     $        nbngr,wkgrd,zorkr,zorki,dxdt,dzdt,
     $        zkt,zkp,nths,0,"p_from_w",lrnge, outmod )
c     
         write ( outmod,'( /,
     $        10x,"total plasma theta, phi skin current   -- ")')
c     
         call kcur0 ( zork3,zork4,zork1,pgrd,mth1,cplgtr,cplgti,n,
     $        nbngr,wkgrd,zorkr,zorki,dxdt,dzdt,
     $        zkt,zkp,nths,0,"Tot. Plas",lrnge, outmod )
c
 800      continue
c     
      end if
c     
c***************end of "if ( lspark .eq. 2 ) then" ************
c     
      return
      end
c
c............................................................
      subroutine kdis ( omsq, blr,bli )
c............................................................
c
      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum5.inc'
      INCLUDE 'vacuum6.inc'
c
      dimension blr(*),bli(*), thph(nths), thgr(nths),
     $     chwr1(nths), chwi1(nths),
     $     dxdt(nths),dzdt(nths),xgrd(nths),zgrd(nths)
      dimension skinx(nths,nths),skinz(nths,nths),
     $     skinp(nths,nths),skint(nths,nths)
c
      dimension chlagdy(nths,nfm)
c     
      common / bigv1 / grdgre(nths2,nths2), grwp(nths,nths),
     $     grri(nths2,nfm2)
c
      real nq
c
c$$$      parameter ( neqv1=(nths-1)/2+1,neqv2=nths+1,neqv3=3*(nths+1)/2 )
c
c$$$      equivalence ( grdgre(1,1), skinx(1,1) )
c$$$      equivalence ( grdgre(neqv2,neqv1), skinp(1,1) )
c$$$      equivalence ( grdgre(1,neqv2), skinz(1,1) )
c$$$      equivalence ( grdgre(neqv2,neqv3), skint(1,1) )
c
      equivalence ( grdgre(1,1), skinx(1,1) )
      equivalence ( grdgre(nptg2x,nptg12), skinp(1,1) )
      equivalence ( grdgre(1,nptg13), skinz(1,1) )
      equivalence ( grdgre(nptg2x,nptg14), skint(1,1) )

      character*(80) string(10)
c     
c     
      call atpoint("KDIS","lrnge",lrnge,"bli(1)",bli(1),
     $     iotty, outmod )
c
c.. New variable, lfarw for SPARK2 to reset farwal.
c
c$$$      lfarw = 0

c$$$      if ( a .gt. 10.0 ) lfarw = 1
c$$$      write ( iotty, '("a = ", 1pe12.3, " lfarw in KDIS = ", i4 )' )
c$$$     $     a, lfarw
c$$$      write ( outmod, '("a = ", 1pe12.3, " lfarw in KDIS = ", i4 )' )
c$$$     $     a, lfarw
c
      if ( lfarw .eq. 0 ) then
c
c     get min and max of wall coords.
c
         zmax = -1.0
         xmax = -1.0
         do  iw = 1, mw
c************************************************dylee,03/18/97
            if ( xpass(iw) .gt. xmax) then
               xmax = xpass(iw)
               if (zpass(iw) .ge. 0.) imwxmax = iw
            end if
c************************************************
c     if ( xpass(iw) .gt. xmax ) xmax = xpass(iw)
            if ( zpass(iw) .gt. zmax ) zmax = zpass(iw)
         end do
         rmax = xmax
         if ( zmax .gt. xmax ) rmax = zmax
c
      end if
c     
      pi = pye
      tpn = twopi * n
      ngr = noutv
      jmax1 = lrnge
c 
c ****************************** for outboard gap *****dylee,03/18/97
c
      if ((ishape .eq. 340) .or. (ishape .eq. 240)) ngr=2*noutv-1
c
c ***********************************************************
c
c...      if ( ishape .gt. 100 ) ngr = 2*noutv - 1
      ngr2 = (ngr-1)/2 + 1
      ngr4 = (ngr-1)/4 + 1
      ngr34 = 3*(ngr-1)/4 + 1
      f = 0.0
      df = 0.0
c
      isph = 1
      if ( jtop .eq. jbot ) isph = 0
c
      do  i = 1, ngr
         thph(i) = ( 1.0/(ngr-1) ) * (i-1)
      end do
c
c......   Extrtact chi-l on the wall -- both the cos(lth) and
c         the sin(lth) contribution. Note that the index was shifted
c         here.
c
      mth12 = 2.0*mth
      if ( lfarw .gt. 0 ) mth12 = mth
c
      mw = mth - (jtop-jbot-1)
c
      do  l1 = 1, jmax1
         iw = mth - jtop + 2
         icnt = 0
         do i = 1, mth
            if (  (i.gt.jbot) .and. (i.lt.jtop) ) go to 1200
            if ( i .eq. jtop ) iw = 1
            chiwc(iw,l1) = grri(mth+i,l1)
            chiws(iw,l1) = grri(mth+i,jmax1+l1)
            iw = iw + 1
            icnt = icnt + 1
 1200       continue
         end do
c.....fill in last point if torus..
         if ( jtop .ne. jbot ) go to 1210
         chiwc(mw,l1) = chiwc(1,l1)
         chiws(mw,l1) = chiws(1,l1)
 1210    continue
      end do
c
c************************************************************************
      if ( lfarw .gt. 0 ) call framep(jobid,ff)
      if ( lfarw .gt. 0 ) return
c************************************************************************
c
c   Sum over l with blr and bli  here.
c
      do  i = 1, mw
         chwr1(i) = 0.0
         chwi1(i) = 0.0
         do  j1 = 1, jmax1
            chwr1(i) = chwr1(i) + chiwc(i,j1)*blr(j1)
     $           - chiws(i,j1)*bli(j1)
            chwi1(i) = chwi1(i) + chiws(i,j1)*blr(j1)
     $           + chiwc(i,j1)*bli(j1)
         end do
      end do
c
      if ( ismth .lt. 1 ) go to 2060
c
c.....    smooth chi on wall.
c
      call smooth ( chwr1, mth1 )
      call smooth ( chwi1, mth1 )
      call pospl2(wgrd,chwr1,chwi1, 1,mth1,1,'chw-smooth',ismth)
      call pospl1(chwr1,chwi1, 1,mth1,2,'chwc-chws', ismth)
c
      if ( ismth .lt. 2) go to 2060
      call smooth ( chwr1, mth1 )
      call smooth ( chwi1, mth1 )
      call pospl2(wgrd,chwr1,chwi1, 1,mth1,3,'chw-smooth',ismth)
      call pospl1(chwr1,chwi1, 1,mth1,4,'chwc-chws', ismth)
      call framep(jobid,ff)
 2060 continue
c
c........    plot current map on raw and equal arc grids...
c... If lkplt = 0 only equal arc. If 1 both raw and eqjual arc.
c... The theta grid is in thgr(i) ...                 !! MSC June 7, 2011
c
      ipl12 = 1
      if ( lkplt .eq. 0 ) ipl12 = 2
c
      do 1900 imap = ipl12, 2
c
      do 40 i = 1, ngr
c
      if ( imap .eq. 2 ) go to 30
c
c.....   raw grid......
c
      thgr(i) = ( 1.0/(ngr-1) ) * (i-1)
c
      go to 35
c
c....    equal arcs grid
c
   30 continue
c
      elgr = ( ell(mw)/(ngr-1) ) * (i-1)
      call lag ( ell,wgrd,mw,3,elgr,f,df,0 )
      thgr(i) = f
c
   35 continue
   40 continue
c
c......   now get x, z, dx, dz, for the appropriate grid.
c
      do 60 i = 1, ngr
      ttt = thgr(i)
      call lag ( wgrd,xpass,mw,3,ttt,f,df,2 )
      xgrd(i) = f
c      if ( abs(xgrd(i)) .lt. 1.0e-8 ) xgrd(i) = 1.0e-8
      dxdt(i) = df
c      if ( abs(dxdt(i)) .le. 1.e-20 ) dxdt(i) = 1.e-20
      call lag ( wgrd,zpass,mw,3,ttt,f,df,2 )
      zgrd(i) = f
      dzdt(i) = df
   60 continue
c
      if ( check2 ) write(outmod,61) (i,thgr(i),xgrd(i),zgrd(i),i=1,ngr)
   61 format(//,2x,"thgr, xgrd, zgrd = ",/,200(2x,i4,1p3e12.4,/))
c
c.....  find grid points where z is maximum and minimum...
c
      zmax = -100.0
      zmin = 100.0
c
      do 70 i = 1, ngr
      if ( zgrd(i) .lt. zmax ) go to 65
      zmax = zgrd(i)
      izmax = i
   65 if ( zgrd(i) .gt. zmin ) go to 70
      zmin = zgrd(i)
      izmin = i
   70 continue
c
      if ( check2 ) write ( outmod, 75 ) izmin, izmax, zmin, zmax
   75 format ( 1x, "izmin, izmax, zmin, zmax = ", 2i5, 1p2e12.4,/ )
c
      do 900 iphse = 1, nphse
c
      do 80 i = 1, ngr
      do 80 j = 1, ngr
      skinx(i,j) = 0.0
      skinp(i,j) = 0.0
      skinz(i,j) = 0.0
      skint(i,j) = 0.0
   80 continue
c
      zan = n
      if ( abs(n) .lt. 1.e-5 ) zan = 1.0
      do  iph = 1, ngr
c
         phi = thph(iph) - (iphse-1)/(8.0*zan)
         phisav = (iphse-1) / (4.0*zan)
c
         tpnph = tpn*phi
         ctpnp = cos(tpn*phi)
         stpnp = sin(tpn*phi)
c
c... Calculate skin currents.
c
         do ith = 1, ngr
c
            thet = thgr(ith)
            xtzt = sqrt ( dxdt(ith)**2 + dzdt(ith)**2 )
c     
c.....chwr1 part.
c     
            call lag(wgrd,chwr1,mw,3,thet,f,df,2)
            skinx(ith,iph) = skinx(ith,iph) + n*f*stpnp*dxdt(ith)
     .           / ( xgrd(ith)*xtzt )
            skinp(ith,iph) = skinp(ith,iph) + df*ctpnp
     .           / xtzt
            skinz(ith,iph) = skinz(ith,iph) + n*f*stpnp*dzdt(ith)
     .           / ( xgrd(ith)*xtzt )
            skint(ith,iph) = skint(ith,iph) + n*f*stpnp / xgrd(ith)
c     
c.....chwi1 part.
c     
            call lag(wgrd,chwi1,mw,3,thet,f,df,2)
            skinx(ith,iph) = skinx(ith,iph) - n*f*ctpnp*dxdt(ith)
     .           / ( xgrd(ith)*xtzt )
            skinp(ith,iph) = skinp(ith,iph) + df*stpnp
     .           / xtzt
            skinz(ith,iph) = skinz(ith,iph) - n*f*ctpnp*dzdt(ith)
     .           / ( xgrd(ith)*xtzt )
            skint(ith,iph) = skint(ith,iph) - n*f*ctpnp / xgrd(ith)
c     
         end do
      end do
c
c......  find largest current.
c
      cmax = -1.0
      do i = 1, ngr
         do  j = 1, ngr
            cur = skinx(i,j)**2 + skinp(i,j)**2 + skinz(i,j)**2
            if ( cur .lt. cmax ) go to 250
            cmax = cur
            imx = i
            jmx = j
 250        continue
         end do
      end do
c
      cmax = sqrt(cmax)
c
      if ( iphse .gt. 1 ) go to 500
c
c....... normalize such that largest current
c        is proportional to grid spacing.
c        plot the current vectors.
c
      idelg = delg
      dartl = float(idelg)/10.0
      darth0 = delg - 10.0*dartl
c
      cnorm = dartl / ( 2.0*cmax*(ngr-1) )
      darth = darth0 / sqrt(4.0*cmax*cnorm)
      elt = ell(mw)
c     ccc = 0.5 * sqrt ( elt**2 + (2.0*pi*r)**2 )
c
      write ( outmod,1250 ) imx,jmx,cmax,cnorm,dartl,darth
 1250 format ( 1x, "imx,jmx,cmax,cnorm,dartl,darth in KDIS= ",/,
     $     2i4, 1p4e13.5,// )
c
      amax = 1.1
      amin = -.1
      bmax = 1.1
      bmin = -.1
c
c      if ( ishape .gt. 100 ) then
c         amax = 1.1/2.0
c         amin = -.1/2.0
c         bmax = 1.1/2.0
c         bmin = -.1/2.0
c      end if
c
      xw = amax + (amax-amin)/100.0
      yw = bmax - (bmax-bmin)/20.0
c
      call maps ( amin,amax, bmin,bmax, .114,.828, .285,1. )
c     
c *****************************************************************
c ******* added for outboard gap **** dylee 03/18/97 *******
c 
      if ((ishape .eq. 340) .or. (ishape .eq. 240)) then
c
c ....  resolution already doubled above : ngr=2*noutv-1
c
c ..... find xmaxzpl (xmax,z>0) and corresponding i, ixmaxzpl
c
         xmaxzpl = -1.0
c
         do 270 i = 1, ngr
          if ((xgrd(i) .gt. xmaxzpl) .and. (zgrd(i) .ge. 0.)) then
           xmaxzpl=xgrd(i)
           ixmaxzpl=i
          end if
 270    continue
c
c.... plot over region from i=1 to (xmax,z>0) point
c
        do 278 j = 1, ngr, 4
        do 278 i = 1, ixmaxzpl
c          x0 = (1.0/(ngr-1))*(i-1)
          x0 = (1.0/(ixmaxzpl-1))*(i-1)
          y0 = (1.0/(ngr-1))*(j-1)
c
c          idy = (ngr-1)/2 + 1 - i
          idy = i
          xtzt = sqrt ( dxdt(idy)**2+dzdt(idy)**2 )
c
c....ct and cp are the theta and phi components, respectively.
c
       ct = skint(idy,j)
       cp = skinp(idy,j)
c
c..... take care of distortion of coordinate grid and ensure that
c      length of vectors is proportional to the magnitude of the
c      original current.
c
       elp = 2.0 * pi * xgrd(i)
       edy = ell(imwxmax)
       cc = sqrt ( (cp*cp+ct*ct)/(edy*edy*cp*cp+elp*elp*ct*ct) )
c......   set cc = constant. for now
c      cc = 1.0
       cndl = 6.0*cnorm
       vp = cndl * cc * edy * cp
       vt = cndl * cc * elp * ct
       call dart(x0,y0, vt,vp, darth )
 278   continue
c
      go to 996
c
      end if
c *******************************************************************
c
      do 280 j = 1, ngr
      do 280 i = 1, ngr
      x0 = thph(i)
      y0 = thph(j)
      xtzt = sqrt ( dxdt(i)**2+dzdt(i)**2 )
c
c....ct and cp are the theta and phi components, respectively.
c
      ct = skint(i,j)
      cp = skinp(i,j)
c
c..... take care of distortion of coordinate grid and ensure that
c      length of vectors is proportional to the magnitude of the
c      original current.
c
      elp = 2.0 * pi * xgrd(i)
      cc = sqrt ( (cp*cp+ct*ct)/(elt*elt*cp*cp+elp*elp*ct*ct) )
c......   set cc = constant. for now
c     cc = 1.0
      vp = cnorm * cc * elt * cp
      vt = cnorm * cc * elp * ct
      call dart(x0,y0, vt,vp, darth )
  280 continue
c     
c
  996 continue
c
      if ( imap .eq. 1 ) 
     .call wrtgr1 ( jobid, xw,yw, cmax, n, q, omsq, lfm,
     $     xili,"xili(l)", jmax1 )
      if ( imap .eq. 2) 
     .call wrtgr1 ( jobid, xw,yw, cmax, n, q, omsq, lfm,
     $     xilr,"xilr(l)", jmax1 )
c
      xw = amin + (amax-amin)/30.0
      yw = bmax - (bmax-bmin)/25.0
      call setlch ( xw, yw, 1,1,0,-1 )
      write ( string, '("Phi vs. Theta Projection, Dart-l,h = ",
     $     1p2e10.2)' ) dartl, darth
      call wrtstr ( string, 1 )
c
      xw = amin - (amax-amin)/10.0
      yw = bmin - (bmax-bmin)/10.0
      call setlch ( xw, yw, 1,1,0,-1 )
c
      if ( imap .eq. 1 ) then
         write ( string,1004 ) ishape
         call wrtstr ( string, 1 )
      end if
 1004 format ( 1x, "ishape = ", i3, " Computation Grid" )
c
c
c ************************************************************************
c ************** modified for outboard gap *** dylee 03/17/97 *****
c ... ell(imwxmax) = wall length from i=1 to imwxmax (= ixmaxzpl roughly)
c
      if ((ishape .eq. 340) .or. (ishape .eq. 240)) then
         elll = ell(imwxmax)
      else
         elll = ell(mw)
      end if
c
      if ( imap .eq. 2 ) then
         write ( string, 1005 ) ishape, elll
         call wrtstr ( string, 1 )
      end if
 1005 format ( 1x, "ishape = ", i3,
     $     " Equal Arcs, Wall Length = ", f8.3 )
c
c ************************************************************************
c
      ngr1 = 1
      write ( string,1003 ) thph(ngr1),xgrd(ngr1),zgrd(ngr1),
     $     thph(ngr4),xgrd(ngr4),zgrd(ngr4),
     $     thph(ngr2),xgrd(ngr2),zgrd(ngr2)
      call wrtstr ( string, 3 )
 1003 format ( 1x," x,z at",f7.3, " = ", 2f7.3,/,
     $     1x," x,z at",f7.3, " = ", 2f7.3,/,
     $     1x," x,z at",f7.3, " = ", 2f7.3 )
c
      call framep(jobid,ff)
c
      if ( (imap .eq. 1) .or. (lkplt .eq. 0) ) go to 900
c
c..... now look at it from directly above the pole.
c
      amin = -1.1 * rmax
      amax =  1.1 * rmax
      bmin = -1.1 * rmax
      bmax =  1.1 * rmax
c
      xw = amax + (amax-amin)/100.0
      yw = bmax - (bmax-bmin)/30.0
c
      cnorm = sqrt(dartl)*rmax / ( cmax*(ngr-1) )
      darth = darth0 / sqrt(4.0*cmax*cnorm)
      call maps ( amin, amax, bmin, bmax, .114,.828, .285,1. )
c
      ngrw = ngr2
      if ( ishape .gt. 100 ) ngrw = ngr4
      do 380 j = 1, ngr
      do 380 i = 1, ngrw
      phi = 2.0*pi*thph(j)
      x0 = xgrd(i) * cos(phi)
      y0 = xgrd(i) * sin(phi)
      vx = cnorm * ( skinx(i,j)*cos(phi)
     .               - skinp(i,j)*sin(phi) )
      vy = cnorm * ( skinx(i,j)*sin(phi)
     .               + skinp(i,j)*cos(phi) )
      call dart ( x0, y0, vx, vy, darth )
  380 continue
c     
      call atpoint("380","ngrw",ngrw,"cmax",cmax,
     $     iotty, outmod )
c
      xw1 = amin + (amax-amin) / 5.0
      yw1 = bmin - (bmax-bmin) / 8.0
      call setlch ( xw1,yw1, 0,1,0,-1 )
      call gtext ("View of the top from above", -1, 0 )
c
c      call setlch ( xw, yw, 1,0,0,0 )
c      write(100,1001) q, n, omsq, ( lfm(i),xilr(i),i=1,jmax1 )
c
      call wrtgr1 ( jobid, xw,yw, cmax, n, q, omsq, lfm,
     $     xili,"xili(l)", jmax1 )
c
      call framep(jobid,ff)
c
c..... now look at the bottom, but still from the top.
c
      call maps ( amin, amax, bmin, bmax, .114,.828, .285,1. )
c
      ngrw = ngr2
      if ( ishape .gt. 100 ) ngrw = ngr34
      do 480 j = 1, ngr
      do 480 i = ngrw, ngr
      phi = 2.0*pi*thph(j)
      x0 = xgrd(i) * cos(phi)
      y0 = xgrd(i) * sin(phi)
      vx = cnorm * ( skinx(i,j)*cos(phi)
     .               - skinp(i,j)*sin(phi) )
      vy = cnorm * ( skinx(i,j)*sin(phi)
     .               + skinp(i,j)*cos(phi) )
      call dart ( x0, y0, vx, vy, darth )
  480 continue
c
      xw1 = amin + (amax-amin) / 5.0
      yw1 = bmin - (bmax-bmin) / 8.0
      call setlch ( xw1,yw1, 0,1,0,-1 )
      call gtext ("View of the bottom from above", -1, 0 )
c
c      call setlch ( xw, yw, 1,0,0,0 )
c      write(100,1001) q, n, omsq, ( lfm(i),xilr(i),i=1,jmax1 )
c
      call wrtgr1 ( jobid, xw,yw, cmax, n, q, omsq, lfm,
     $     xilr,"xilr(l)", jmax1 )
c
      call framep(jobid,ff)
c
  500 continue
c
      if ( imap .eq. 1 ) go to 900
c
      call maps ( amin, amax, bmin, bmax, .114,.828, .285,1. )
c
      do 580 j = 1, ngr2
      do 580 i = izmax, izmin
c
      phi = 2.0*pi*thph(j)
      x0 = -xgrd(i) * cos(phi)
      z0 = zgrd(i)
      vx = -cnorm * ( skinx(i,j)*cos(phi)
     .               -skinp(i,j)*sin(phi) )
      vz = cnorm * skinz(i,j)
      call dart ( x0, z0, vx, vz, darth )
  580 continue
c
      xw1 = amin + (amax-amin) / 10.0
      yw1 = bmax - (bmax-bmin) / 10.0
      call setlch ( xw1, yw1, 0,1,0,-1 )
      write ( string, '("Phi = ", f7.2," Pi")' ) phisav
      call wrtstr ( string, 1 )
c
c      call setlch ( xw, yw, 1,0,0,0 )
c      write(100,1001) q, n, omsq, ( lfm(i),xilr(i),i=1,jmax1 )
c
      call wrtgr1 ( jobid, xw,yw, cmax, n, q, omsq, lfm,
     $     xilr,"xilr(l)", jmax1 )
c
      call framep(jobid,ff)
c
c...... look at small major radius side of torus...
c
      if ( isph .eq. 1 ) go to 900
c
      call maps ( amin, amax, bmin, bmax, .114,.828, .285,1. )
c
      do 680 j = ngr2, ngr
      phi = 2.0*pi*thph(j)
      do 660 i = izmin, ngr
c
      x0 = -xgrd(i) * cos(phi)
      z0 = zgrd(i)
      vx = -cnorm * ( skinx(i,j)*cos(phi)
     .              -skinp(i,j)*sin(phi) )
      vz = cnorm * skinz(i,j)
      call dart ( x0, z0, vx, vz, darth )
  660 continue
      do 680 i = 2, izmax
      x0 = -xgrd(i) * cos(phi)
      z0 = zgrd(i)
      vx = -cnorm * ( skinx(i,j)*cos(phi)
     .              -skinp(i,j)*sin(phi) )
      vz = cnorm * skinz(i,j)
      call dart ( x0, z0, vx, vz, darth )
c
c      xw1 = amin + (amax-amin) / 10.0
c      yw1 = bmax - (bmax-bmin) / 10.0
      call setlch ( xw1, yw1, 0,1,0,-1 )
      write ( string, '("twopin-phi = ", f7.2)' ) tpnph
      call wrtstr ( string, 1 )
c
c.....outline outer edge of torus.
c
      if ( .not. ( (j .eq. ngr) .or. (j .eq. ngr2) ) ) go to 680
      do 670 ii = izmax, izmin
      x0 = xgrd(ii) * cos(phi)
      z0 = zgrd(ii)
      call pointc ( '+', x0, z0, 1, -1, -1, 0., 0. )
  670 continue
  680 continue
c
c      call setlch ( xw, yw, 1,0,0,0 )
c      write(100,1001) q, n, omsq, ( lfm(i),xilr(i),i=1,jmax1 )
c
      call wrtgr1 ( jobid, xw,yw, cmax, n, q, omsq, lfm,
     $     xilr,"xilr(l)", jmax1 )
c
 1001 format(1x,"qs =",f5.2,/,2x,"n =",f4.1,/
     .       1x,"o=",e12.5,/,3x,"l",3x,"xilr(l)",/,
     .       100(i4,1pe11.3,/) )
c
      call framep(jobid,ff)
c
  900 continue
c
 1900 continue
c
c     
      call atpoint("end KDIS","lrnge",lrnge,"bnli(1)",bnli(1),
     $     iotty, outmod )
      return
      end
c     
c......................................................................
      subroutine bnsrf ( xpass,zpass,ell,thlag,mw,bnpwtr,bnpwti,
     $     nbngr,thgr,xgrd,zgrd,bngrdr,bngrdi,tbibr,abibr,nout1 )
c.....................................................................
c     
c...thgr,xgrd,zgrd,bngrdr,bngrdi,tbibr,abibr are dummy variables.
c
      dimension xpass(*), zpass(*), ell(*), thlag(*),
     $     bnpwtr(*), bnpwti(*), thgr(*), xgrd(*), zgrd(*),
     $     bngrdr(*), bngrdi(*), tbibr(*), abibr(*)
c     
c.....Get equal arcs grid.
c     
      do i = 1, nbngr
         elgr = ( ell(mw)/(nbngr-1) ) * (i-1)
         call lag ( ell,thlag,mw,3,elgr,f,df,0 )
         thgr(i) = f
      end do
c     
c......now get x, z, bnr, bni on surface for the equal arcs grid.
c     
      do  i = 1, nbngr
         ttt = thgr(i)
         call lag ( thlag,xpass,mw,3,ttt,f,df,0 )
         xgrd(i) = f
c         if ( abs(xgrd(i)) .lt. 1.0e-8 ) xgrd(i) = 1.0e-8
         call lag ( thlag,zpass,mw,3,ttt,f,df,0 )
         zgrd(i) = f
         call lag ( thlag,bnpwtr,mw,3,ttt,f,df,0 )
         bngrdr(i) = f
         call lag ( thlag,bnpwti,mw,3,ttt,f,df,0 )
         bngrdi(i) = f
         tbibr(i) = atan2m ( bngrdi(i),bngrdr(i) )
         abibr(i) = sqrt( bngrdi(i)**2 + bngrdr(i)**2 )
      end do
c     
c     **** find max. of abibr, then normalize, dylee, 4/17/97 ******
      abibrmax=-10.0
      do i=1,nbngr
         if (abibr(i) .gt. abibrmax) then
            abibrmax = abibr(i)
         end if
      end do
c     **************************************************************
c     
      write ( nout1,
     $     '( 10x,"EQUAL ARCS GRID -- THETA = 0 ON INSIDE", /,
     $     "b_norm(i,phi) = bnr(i)*cos(n*phi) + bni(i)*sin(n*phi)",/,
     $     4x,"i",5x,"x",7x,"z",8x,"bnr",8x,"bni",5x,"absbibr",5x,
     $     "atanbibr" )' )
c     
      do i = 1, nbngr
         write ( nout1,2246 ) i, xgrd(i),zgrd(i),
     $        bngrdr(i),bngrdi(i),abibr(i)/abibrmax,tbibr(i)
 2246    format ( 1x, i4, 2f8.3, 1p4e11.3 )
      end do
c     
      return
      end

c...................................................................
      SUBROUTINE MPHASE ( xpass,zpass,ell,thlag,mw,chreal,chimag,n,
     $     nbngr,zccen,zccena,zccenb,
     $     thgr,xgrd,zgrd,dxdt,dzdt,zkt,zkp,nths,ikb,
     $     inlabel,iznum,nout1,nout2 )
c..................................................................
c     
c  Called in RESSHEL_2. Plots the relative phases of the wall mode pattern
c   at the three positions where the coils are centered.

c     Write out information for the theta and phi components 
c     of the skin current at a surface  on an equal arc grid.
c     Only the theta dependence: the phi dependence is analytical.
c     Calls PLOTK to plot currents. 
c     Chi on surface = [Chreal, Chimag].
c     Thlag(mw): Array = [0,1] in theta.
c     ELL(mw): Arc length of shell.
c$$$  
c$$$  To reconstruct: for i along the surface in theta,
c$$$  K_theta(i,phi) = Ktr(i)*sin(n*phi) - Kti(i)*cos(n*phi)
c$$$  K_phi(i,phi)   = Kpr(i)*cos(n*phi) + Kpi(i)*sin(n*phi)
c$$$  Ktr, Kti stored in zkt(i,1) zkt(i,2)
c$$$  Kpr, Kpi stored in zkp(i,1) zkp(i,2)
c
c     Input: xpass, zpass, ell, thlag, mw, chreal, chimag, n, nbngr,
c            nout1, nths
c     Dummy: thgr, xgrd, zgrd, dxdt, dzdt, zkt, zkp.
c  
      REAL n
c
      DIMENSION xpass(*), zpass(*), ell(*), thlag(*), chreal(*),
     $     chimag(*), thgr(*), xgrd(*), zgrd(*), dxdt(*), dzdt(*),
     $     zkt(nths,2), zkp(nths,2)
c
      REAL, DIMENSION(:), ALLOCATABLE :: zkt1, zkt2, zkp1, zkp2
      REAL, DIMENSION(:), ALLOCATABLE :: zktit1, zktit2, zkpit1, zkpit2
      REAL, DIMENSION(:), ALLOCATABLE ::  ktphse, kpphse, zccent
      REAL, DIMENSION(:), ALLOCATABLE ::  zmagt, zmagp
      REAL, DIMENSION(:), ALLOCATABLE ::  xobs, zobs   ! local here
      REAL, DIMENSION(:,:), ALLOCATABLE :: zkt3, zkp3
      INTEGER, DIMENSION(:), ALLOCATABLE :: nzt

      CHARACTER*(*) inlabel
      CHARACTER*(20) pllabel0, pllabel
      CHARACTER(130), DIMENSION(10) :: string
      LOGICAL backl
c     
      zpi = acos(-1.0)
      ztwopi = 2.0 * zpi
      backl = .false.

c... Interpolate to nbngr points.:

      DO  i = 1, nbngr
c     
c.... equal arcs grid
c     
         elgr = ( ell(mw)/(nbngr-1) ) * (i-1)
         call lag ( ell,thlag,mw,5,elgr,f,df,0 )
         thgr(i) = f
c     
      END DO
c     
c......get x, z, dx, dz, for equal arcs grid.
c     
      DO  i = 1, nbngr
         ttt = thgr(i)
         call lag ( thlag,xpass,mw,5,ttt,f,df,2 )
         xgrd(i) = f
c         if ( abs(xgrd(i)) .lt. 1.0e-8 ) xgrd(i) = 1.0e-8
         dxdt(i) = df
c         if ( abs(dxdt(i)) .le. 1.e-20 ) dxdt(i) = 1.e-20
         call lag ( thlag,zpass,mw,5,ttt,f,df,2 )
         zgrd(i) = f
         dzdt(i) = df
      END DO
c     
c     Contribution of theta and phi components from real part of Xil 
c     Store in zkt(1,2), and zkp(1,2). 
c     Already summed over B(l)
c     
      DO  ith = 1, nbngr
c     
         thet = thgr(ith)
         xtzt = sqrt ( dxdt(ith)**2 + dzdt(ith)**2 )
c     
c...Contribution from Real part
         CALL lag(thlag,chreal,mw,5,thet,f,df,2)
         zkt(ith,1) = n*f / xgrd(ith)
         zkp(ith,1) = df / xtzt
c     
c...Contribution from Imaginary part.
         CALL lag(thlag,chimag,mw,5,thet,f,df,2)
         zkt(ith,2) = n*f / xgrd(ith)
         zkp(ith,2) = df / xtzt
c     
      END DO
c
      IF ( nout1 .ne. 0 ) THEN
         WRITE ( nout1, '( 
     $        "To reconstruct: for i along the surface in theta,",/,
     $    "K_theta(i,phi) = Ktr(i)*sin(n*phi) - Kti(i)*cos(n*phi)",/,
     $    "K_phi(i,p`hi)  = Kpr(i)*cos(n*phi) + Kpi(i)*sin(n*phi)")' )
c     
         WRITE ( nout1,'( 15x, "Skin Current on ", a , " No. = ", i5,/,
     $        10x,"EQUAL ARCS GRID -- THETA = 0 ON INSIDE", /,
     $        4x,"i",5x,"x",7x,"z",8x,"Ktr",8x,"Kti",8x,"Kpr",
     $        8x,"Kpi" )' ) inlabel, iznum
c     
         DO i = 1, nbngr
            write ( nout1, '( 1x, i4, 2f8.3, 1p4e11.3 )' )
     $           i, xgrd(i),zgrd(i),
     $           zkt(i,1),zkt(i,2), zkp(i,1),zkp(i,2)
         END DO
      END IF
c
c.... Plot Kt and Kp. Use xgrd,zkt1,zkt2, zkp1,zkp2, as storage.
c
      ALLOCATE ( zkt1(nths),zkt2(nths), zkp1(nths),zkp2(nths) )

      DO i = 1, nbngr
         thgr(i) = (i-1) * ztwopi / (nbngr -1) - zpi
      END DO
c
      ALLOCATE ( zkt3(nths,11), zkp3(nths,11), zmagt(11), zmagp(11),
     $     xobs(nths), zobs(nths), nzt(11) )
      ALLOCATE ( zktit1(11), zktit2(11), zkpit1(11), zkpit2(11) )
      ALLOCATE ( ktphse(11), kpphse(11), zccent(11) )

c      Get abscissa starting from -pi. Remember the coil centers
c        are in original coordinates with origin at outer side.


c.. Calculate phase of the theta and phi components, according to:
c    Kt = - A * \cos ( n\phi + \delta ) 
c    Kp =   A * \sin ( n\phi - \delta ) 

      DO i = 1, nbngr
c         thgr(i) = i
         zkt1(i) = ATAN2 ( zkt(i,1),zkt(i,2) )
         zkt2(i) = SQRT ( zkt(i,1)**2 + zkt(i,2)**2 )
         zkp1(i) = ATAN2 ( zkp(i,2),zkp(i,1) )
         zkp2(i) = SQRT ( zkp(i,1)**2 + zkp(i,2)**2 )
      END DO

c.... Get values at the coilcenters.

      zccent(1) = zccena * zpi
      zccent(2) = zccen * zpi
      zccent(3) = zccenb * zpi

      nzcoil = 3

      DO it = 1, nzcoil
         CALL lag ( thgr,zkt1, nbngr, 3, zccent(it), f,df, 0 )
         ktphse(it) = f
         CALL lag ( thgr,zkt2, nbngr, 3, zccent(it), f,df, 0 )
         zmagt(it) = f
         CALL lag ( thgr,zkp1, nbngr, 3, zccent(it), f,df, 0 )
         kpphse(it) = f
         CALL lag ( thgr,zkp2, nbngr, 3, zccent(it), f,df, 0 )
         zmagp(it) = f
      END DO

      nzcoilh = (nzcoil-1) / 2 + 1

      zmagth = zmagt(nzcoilh)
      zmagph = zmagp(nzcoilh)

      DO it = 1, nzcoil
         zmagt(it) = zmagt(it) / zmagth
         zmagp(it) = zmagp(it) / zmagph
      END DO

      WRITE ( nout1, '()' )
      WRITE ( nout1, '( 3x, "K-theta coil mag & phase", 6x,
     $     "K-phi coil mag & phase")' )
      DO it = 1, nzcoil
         WRITE ( nout1, '( 1x, es13.5, 1x, es13.5, 2x,
     $        es13.5, 2x, es13.5 )' )
     $        zmagt(it), ktphse(it),  zmagp(it), kpphse(it)
      END DO

      WRITE ( nout2, '()' )
      WRITE ( nout2, '( 15x, "K-theta coil phase,  K-phi coil phase")' )

      DO it = 1, nzcoil
         WRITE ( nout2, '(17x, es13.5, 5x, es13.5 )' )
     $        ktphse(it), kpphse(it)
      END DO

      zkt2 = 0.0
      zkp2 = 0.0

      pllabel0 = "Tht: "//inlabel
      lena = index ( pllabel0, '  ', backl ) - 1
      pllabel = pllabel0(1:lena)
      CALL pospl2sc ( thgr, zkt1,zkt2, 1,nbngr, 1,pllabel,
     $     zccent, ktphse,nzcoil, iznum )

      pllabel0 = "Phi: "//inlabel
      lena = index ( pllabel0, '  ', backl ) - 1
      pllabel = pllabel0(1:lena)
      CALL pospl2sc ( thgr, zkp1,zkp2, 1,nbngr, 2,pllabel,
     $      zccent, kpphse,nzcoil, iznum )

c... Normalize to pi and the center coil:

      zktphc = ktphse(nzcoilh)
      zkpphc = kpphse(nzcoilh)

      DO it = 1,  nzcoil
         ktphse(it) = ( ktphse(it) - zktphc ) / zpi
         kpphse(it) = ( kpphse(it) - zkpphc ) / zpi
      END DO

      WRITE ( nout1, '()' )
      WRITE ( nout1, '( 15x, "Rel.K-t phase/pi,  Rel. K-p phase/pi")' )

      DO it = 1, nzcoil
         WRITE ( nout1, '( 17x, es13.5, 5x, es13.5 )' )
     $        ktphse(it), kpphse(it)
      END DO

      WRITE ( nout2, '()' )
      WRITE ( nout2, '( 15x, "Rel.K-t phase/pi,  Rel. K-p phase/pi")' )

      DO it = 1, nzcoil
         WRITE ( nout2, '( 17x, es13.5, 5x, es13.5 )' )
     $        ktphse(it), kpphse(it)
      END DO

c. Plot Kt(phi), Kp(phi) at the coil center theta values.
c      WITH AMPLITUDES NORMALIZED TO THE CENTER COIL.
      
      DO i = 1, nbngr
         zkt1(i) = zkt(i,1)
         zkt2(i) = zkt(i,2)
         zkp1(i) = zkp(i,1)
         zkp2(i) = zkp(i,2)
      END DO

      DO it = 1, nzcoil
         CALL lag ( thgr,zkt1, nbngr, 3, zccent(it), f,df, 0 )
         zktit1(it) = f
         CALL lag ( thgr,zkt2, nbngr, 3, zccent(it), f,df, 0 )
         zktit2(it) = f
         CALL lag ( thgr,zkp1, nbngr, 3, zccent(it), f,df, 0 )
         zkpit1(it) = f
         CALL lag ( thgr,zkp2, nbngr, 3, zccent(it), f,df, 0 )
         zkpit2(it) = f
      END DO
      
      azbt = sqrt ( zktit1(nzcoilh)**2 + zktit2(nzcoilh)**2 )
      azbp = sqrt ( zkpit1(nzcoilh)**2 + zkpit2(nzcoilh)**2 )

      DO i = 1, nbngr
         zphi = ztwopi * (i-1) / (nbngr - 1)
         xobs(i) = zphi
         zphi = n * zphi
         sinnph = sin(zphi)
         cosnph = cos(zphi)
         DO it = 1, nzcoil
            zkt3(i,it) =
     $           ( zktit1(it)*sinnph - zktit2(it)*cosnph ) / (2.0*azbt)
            zkp3(i,it) =
     $           ( zkpit1(it)*cosnph + zkpit2(it)*sinnph ) / (2.0*azbp)
         END DO
      END DO


      nobs = nzcoil
      nfobs = 1
      nfobs1 = 1
      cycl = ztwopi
      ncyc = nbngr

c     
      zz1 = nobs + 1.0
      zz0 = nfobs - 1.0
      CALL mapg(0.0,cycl,zz0,zz1,0.1,.4,.30,0.60)
c     
      DO i = nfobs1, nobs
         DO iph = 1, ncyc
            zobs(iph) = zkt3(iph,i) + i 
         END DO
         CALL trace ( xobs, zobs, ncyc, 1,1,0.,0. )
      END DO

      xw = 0.0
      yw = zz0 - 0.30*(nobs-nfobs)

      CALL setlch ( xw, yw, 1,1,0,-1 )

      WRITE ( string,
     $     '(10x,"K_theta and K_phi on shell vs. phi. Rel. amp.")' )

      CALL wrtstr ( string, 2 )

      xw = 0.0
      yw = zz0 - 0.4*(nobs-nfobs)

      CALL setlch ( xw, yw, 1,1,0,-1 )

      WRITE ( string, '()' )
      CALL wrtstr ( string, 1 )

      WRITE ( string, '( 5x, "Rel.K-t phase/pi", 17x,
     $     "Rel. K-p phase/pi")' )
      CALL wrtstr ( string, 1 )

      DO it = 1, nzcoil
         WRITE ( string, '( 5x, es13.5, 20x, es13.5 )' )
     $        ktphse(it), kpphse(it)
      CALL wrtstr ( string, 1 )
      END DO

      CALL mapg(0.0,cycl,zz0,zz1,.5,.8,.30,0.60)
c     
      DO i = nfobs1, nobs
         DO iph = 1, ncyc
            zobs(iph) = zkp3(iph,i) + i 
         END DO
         CALL trace ( xobs, zobs, ncyc, 1,1,0.,0. )
      END DO

      CALL framep(jobid,ff)

      DEALLOCATE ( zkt1, zkt2, zkp1, zkp2, zmagt,zmagp )
      DEALLOCATE ( zkt3, zkp3,xobs,zobs, nzt )
      DEALLOCATE ( zktit1, zktit2, zkpit1, zkpit2 )
      DEALLOCATE ( ktphse, kpphse, zccent )

      RETURN
      END
c     

c...................................................................
      subroutine kcur0 ( xpass,zpass,ell,thlag,mw,chreal,chimag,n,
     $     nbngr,thgr,xgrd,zgrd,dxdt,dzdt,zkt,zkp,nths,ikb,
     $     inlabel,iznum,nout1 )
c..................................................................
c     
c     Write out information for the theta and phi components 
c     of the skin current at a surface  on an equal arc grid.
c     Only the theta dependence: the phi dependence is analytical.
c     Calls PLOTK to plot currents. 
c     Chi on surface = [Chreal, Chimag].
c     Thlag(mw): Array = [0,1] in theta.
c     ELL(mw): Arc length of surface.

c.. Can plot the surface component of the magnetic field if 
c    K_phi -> Q_theta, and K_theta -> - Q_phi. 
c  ikb = 1 for this  option. That is:

c$$$      IF ( ikb == 1 ) THEN
c$$$         DO i = 1, nbngr
c$$$            zkp(i,1) =   zkt2(i)
c$$$            zkp(i,2) = - zkt1(i)
c$$$            zkt(i,1) =   zkp2(i)
c$$$            zkt(i,2) = - zkp1(i)
c$$$         END DO

c$$$  To reconstruct: for i along the surface in theta,
c$$$  K_theta(i,phi) = Ktr(i)*sin(n*phi) - Kti(i)*cos(n*phi)
c$$$  K_phi(i,phi)   = Kpr(i)*cos(n*phi) + Kpi(i)*sin(n*phi)
c$$$  Ktr, Kti stored in zkt(i,1) zkt(i,2)
c$$$  Kpr, Kpi stored in zkp(i,1) zkp(i,2)
c
c     Input: xpass, zpass, ell, thlag, mw, chreal, chimag, n, nbngr,
c            nout1, nths
c     Dummy: thgr, xgrd, zgrd, dxdt, dzdt, zkt, zkp.
c  
      real n
c
      dimension xpass(*), zpass(*), ell(*), thlag(*), chreal(*),
     $     chimag(*), thgr(*), xgrd(*), zgrd(*), dxdt(*), dzdt(*),
     $     zkt(nths,2), zkp(nths,2)
c
      REAL, DIMENSION(:), ALLOCATABLE :: zkt1, zkt2, zkp1, zkp2
      REAL, DIMENSION(:), ALLOCATABLE :: zgrdpi

      character*(*) inlabel
      character*(25) pllabel0, pllabel
      LOGICAL backl

      backl = .false.
c     
      do  i = 1, nbngr
c     
c.... equal arcs grid
c     
         elgr = ( ell(mw)/(nbngr-1) ) * (i-1)
         call lag ( ell,thlag,mw,3,elgr,f,df,0 )
         thgr(i) = f
c     
      end do
c     
c......get x, z, dx, dz, for equal arcs grid.
c     
      do  i = 1, nbngr
         ttt = thgr(i)
         call lag ( thlag,xpass,mw,3,ttt,f,df,2 )
         xgrd(i) = f
c         if ( abs(xgrd(i)) .lt. 1.0e-8 ) xgrd(i) = 1.0e-8
         dxdt(i) = df
c         if ( abs(dxdt(i)) .le. 1.e-20 ) dxdt(i) = 1.e-20
         call lag ( thlag,zpass,mw,3,ttt,f,df,2 )
         zgrd(i) = f
         dzdt(i) = df
      end do
c     
c     Contribution of theta and phi components from real part of Xil 
c     Store in zkt(1,2), and zkp(1,2). 
c     Already summed over B(l)
c     
      do  ith = 1, nbngr
c     
         thet = thgr(ith)
         xtzt = sqrt ( dxdt(ith)**2 + dzdt(ith)**2 )
c     
c...Contribution from Real part
         call lag(thlag,chreal,mw,3,thet,f,df,2)
         zkt(ith,1) = n*f / xgrd(ith)
         zkp(ith,1) = df / xtzt
c     
c...Contribution from Imaginary part.
         call lag(thlag,chimag,mw,3,thet,f,df,2)
         zkt(ith,2) = n*f / xgrd(ith)
         zkp(ith,2) = df / xtzt
c     
      end do
c
      if ( nout1 .ne. 0 ) then
         Write ( nout1, '( 
     $        "To reconstruct: for i along the surface in theta,",/,
     $    "K_theta(i,phi) = Ktr(i)*sin(n*phi) - Kti(i)*cos(n*phi)",/,
     $    "K_phi(i,p`hi)  = Kpr(i)*cos(n*phi) + Kpi(i)*sin(n*phi)")' )
c     
         write ( nout1,'( 15x, "Skin Current on ", a , " No. = ", i5,/,
     $        10x,"EQUAL ARCS GRID -- THETA = 0 ON INSIDE", /,
     $        4x,"i",5x,"x",7x,"z",8x,"Ktr",8x,"Kti",8x,"Kpr",
     $        8x,"Kpi" )' ) inlabel, iznum
c     
         do i = 1, nbngr
            write ( nout1, '( 1x, i4, 2f8.3, 1p4e11.3 )' )
     $           i, xgrd(i),zgrd(i),
     $           zkt(i,1),zkt(i,2), zkp(i,1),zkp(i,2)
         end do
      end if
c
c.... Plot Kt and Kp. Use xgrd,zkt1,zkt2, zkp1,zkp2, as storage.
c
      ALLOCATE ( zkt1(nths),zkt2(nths), zkp1(nths),zkp2(nths) )
      ALLOCATE ( zgrdpi(nths) )

      zpi = acos(-1.0)
      ztwopi = 2.0 * zpi

      DO i = 1, nbngr
         zgrdpi(i) = (i-1) * ztwopi / ( nbngr-1 ) - zpi
         zkt1(i) = zkt(i,1)
         zkt2(i) = zkt(i,2)
         zkp1(i) = zkp(i,1)
         zkp2(i) = zkp(i,2)
      END DO

      pllabel0 = "K-t: "//inlabel
      lena = LEN_TRIM( pllabel0 )
c      lena = index ( pllabel0, '  ', backl ) - 1
      pllabel = pllabel0(1:lena)
      CALL pospl2s ( zgrdpi, zkt1,zkt2, 1,nbngr, 1,pllabel,iznum )
      pllabel0 = "K-p: "//inlabel
c      lena = index ( pllabel0, '  ', backl ) - 1
      lena = LEN_TRIM( pllabel0 )
      pllabel = pllabel0(1:lena)
      CALL pospl2s ( zgrdpi, zkp1,zkp2, 1,nbngr, 2,pllabel,iznum )

c
c.....Plot current vectors on phi-theta plane. Use dxdt, dzdt as dummy.
c
      elt = ell(mw)
c$$$      call plotk ( zkt, zkp, xgrd,zgrd, elt, nbngr, 2, inlabel,
c$$$     $     dxdt,dzdt ) 

c.. Can plot the surface component of the magnetic field if 
c    K_phi -> Q_theta, and K_theta -> - Q_phi


      IF ( ikb == 0 ) THEN
         call plotk ( zkt, zkp, xgrd,zgrd, elt, nbngr, 4, inlabel,
     $        iznum,dxdt,dzdt )
      END IF

      IF ( ikb == 1 ) THEN
         DO i = 1, nbngr
            zkp(i,1) =   zkt2(i)
            zkp(i,2) = - zkt1(i)
            zkt(i,1) =   zkp2(i)
            zkt(i,2) = - zkp1(i)
         END DO
         call plotk ( zkt, zkp, xgrd,zgrd, elt, nbngr, 4, inlabel,
     $        iznum,dxdt,dzdt )
      END IF

      DEALLOCATE ( zkt1, zkt2, zkp1, zkp2, zgrdpi )

      RETURN
      END
c     
c............................................................
      subroutine plotk (zkt, zkp, xgrd,zgrd,elt, ngrd,ndg,inlabel,
     $     iznum,thgr,thph )
c............................................................
c     
c..  Plots the skin current from zkt, zkp (the theta and phi components)
c     data  given on a ngrd grid. in xgrd(theta), zgrd(theta).
c    It plots every ndg points. 

c.. Can plot the surface component of the magnetic field if 
c    it is called with:
c    K_phi -> Q_theta, and K_theta -> - Q_phi

c
      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum7.inc'
c     
      common / bigv1 / grdgre(nths2,nths2), grwp(nths,nths),
     $     grri(nths2,nfm2)

      dimension zkt(nths,2), zkp(nths,2), thph(*), thgr(*),
     $     xgrd(*),zgrd(*)
      dimension skint(nths,nths),skinp(nths,nths)
c     
c$$$      parameter ( neqv1=(nths-1)/2+1,neqv2=nths+1,neqv3=3*(nths+1)/2 )
c$$$      equivalence ( grdgre(1,1), skint(1,1) )
c$$$      equivalence ( grdgre(neqv2,neqv1), skinp(1,1) )
c      equivalence ( grdgre(1,neqv2), skinz(1,1) )

c      equivalence ( grdgre(1,1), skinx(1,1) )
      equivalence ( grdgre(nptg2x,nptg12), skinp(1,1) )
c      equivalence ( grdgre(1,nptg13), skinz(1,1) )
      equivalence ( grdgre(nptg2x,nptg14), skint(1,1) )
c     
      character*(*) inlabel
      character*(80) string(10)
c     
      pi = pye
      tpn = twopi * n
c      ngr = noutv
      ngr = (ngrd-1) / ndg + 1
      jmax1 = lrnge
c     
      ngr2 = (ngr-1)/2 + 1
      ngr4 = (ngr-1)/4 + 1
      ngr34 = 3*(ngr-1)/4 + 1
c     
      do  i = 1, ngr
         thgr(i) = ( 1.0/(ngr-1) ) * (i-1)
         thph(i) = ( 1.0/(ngr-1) ) * (i-1)
         do  j = 1, ngr
            skint(i,j) = 0.0
            skinp(i,j) = 0.0
         end do
      end do
c     
      do iph = 1, ngr
         phi = thph(iph)
         ctpnp = cos(tpn*phi)
         stpnp = sin(tpn*phi)
c...  Calculate skin currents.
         do ith = 1, ngr
            it0 = (ith-1)*ndg + 1
            skint(ith,iph) = zkt(it0,1)*stpnp - zkt(it0,2)*ctpnp
            skinp(ith,iph) = zkp(it0,1)*ctpnp + zkp(it0,2)*stpnp
         end do
      end do
c     
c......find largest current.
c     
      cmax = -1.0
      do i = 1, ngr
         do  j = 1, ngr
            cur = skint(i,j)**2 + skinp(i,j)**2
            if ( cur .lt. cmax ) go to 250
            cmax = cur
            imx = i
            jmx = j
 250        continue
         end do
      end do
c     
      if ( cmax .lt. 1.0e-20) cmax = 1.0e-20
      cmax = sqrt(cmax)
c     
c.......normalize such that largest current
c     is proportional to grid spacing.
c     plot the current vectors.
c     
      idelg = delg
      dartl = float(idelg)/10.0
      darth0 = delg - 10.0*dartl
c
      call atpoint("kplot 1","ngr",ngr,"cmax",cmax,
     $     iotty, outmod )
c     
      cnorm = dartl / ( 2.0*cmax*(ngr-1) )
      darth = darth0 / sqrt(4.0*cmax*cnorm)
c
c     ccc = 0.5 * sqrt ( elt**2 + (2.0*pi*r)**2 )
c     
      write ( outmod,1250 ) inlabel, iznum, imx,jmx,cmax,cnorm,
     $     dartl,darth
 1250 format ( 1x, "imx,jmx,cmax,cnorm,dartl,darth in: ", a, " No. ",
     $     i4, " = " ,/, 2i4, 1p4e13.5,// ) 
c     
      amax = 1.1
      amin = -.1
      bmax = 1.1
      bmin = -.1
c     
c     if ( ishape .gt. 100 ) then
c     amax = 1.1/2.0
c     amin = -.1/2.0
c     bmax = 1.1/2.0
c     bmin = -.1/2.0
c     end if
c     
      xw = amax + (amax-amin)/100.0
      yw = bmax - (bmax-bmin)/20.0
c     
      call maps ( amin,amax, bmin,bmax, .114,.714, 0.12,0.72 )
c     
      do  j = 1, ngr
         do  i = 1, ngr
            x0 = thph(i)
            y0 = thph(j)
c     
c.... ct and cp are the theta and phi components, respectively.
c     
            ct = skint(i,j)
            cp = skinp(i,j)
c     
c.....take care of distortion of coordinate grid and ensure that
c     length of vectors is proportional to the magnitude of the
c     original current.
c
c$$$      call atpoint("kplot 2","i",i,"ct",ct,
c$$$     $     iotty, outmod )
c$$$      call atpoint("kplot 2","j",j,"cp",cp,
c$$$     $     iotty, outmod )
c     
            i0 = (i-1)*ndg + 1
            elp = 2.0 * pi * xgrd(i0)
            cc0 = cp*cp + ct*ct
            cc1 = elt*elt*cp*cp + elp*elp*ct*ct
            if ( cc1 .le. epszer ) then
               cc1 = epszer
               cc = 1.0
c$$$      call atpoint("kplot 2","j",j,"cc",cc,
c$$$     $     iotty, outmod )
            else
               cc = sqrt ( cc0/cc1 )
            end if
c
c......set cc = constant. for now
c     cc = 1.0
            vp = cnorm * cc * elt * cp
            vt = cnorm * cc * elp * ct
c$$$      call atpoint("kplot 3","j",j,"vp",vp,
c$$$     $     iotty, outmod )            
c$$$      call atpoint("kplot 3","j",j,"vt",vt,
c$$$     $     iotty, outmod )            
            call dart(x0,y0, vt,vp, darth )
c$$$c           call atpoint("kplot 4","j",j,"vp",vp,
c$$$     $     iotty, outmod )            
         end do
      end do
c     
      call wrtgr1 ( jobid, xw,yw, cmax, n, q, omsq, lfm,
     $     bnli,"bnli(l)", jmax1 )
c     
      xw = amin + (amax-amin)/30.0
      yw = bmax - (bmax-bmin)/25.0
      call setlch ( xw, yw, 1,1,0,-1 )
c$$$      write ( string, '("Phi vs. Theta Projection, Dart-l,h = ",
c$$$     $     1p2e10.2)' ) dartl, darth
      write ( string,'( 1x, " Skin Current for: ", a, " No. ", i4 )' )
     $     inlabel, iznum
      call wrtstr ( string, 1 )
c     
      xw = amin - (amax-amin)/10.0
      yw = bmin - (bmax-bmin)/9.0
      call setlch ( xw, yw, 1,1,0,-1 )
c     
      ngr1 = 1
      ng20 = (ngr2-1)*ndg + 1
      ng40 = (ngr4-1)*ndg + 1
      write ( string,1003 ) thph(ngr1),xgrd(ngr1),zgrd(ngr1),
     $     thph(ngr4),xgrd(ng40),zgrd(ng40),
     $     thph(ngr2),xgrd(ng20),zgrd(ng20)
      call wrtstr ( string, 1 )
 1003 format ( 1x,"(L x z) = ", "(",3f6.3,"),",
     $     1x, "(",3f6.3,")," 1x,"(",3f6.3,")" )
c     
      if ( lfbcoil .eq. 1 ) then
         write ( string,
     $        '(1x,"taucoil = ",f7.3,2x,"coilang = ", f7.3)' )
     $        taucoil, coilang
         call wrtstr ( string, 1 )
      end if

 1001 format(1x,"qs =",f5.2,/,2x,"n =",f4.1,/
     .     1x,"o=",e12.5,/,3x,"l",3x,"bnlr(l)",/,
     .     100(i4,1pe11.3,/) )
c     
      call framep(jobid,ff)
c     
      return
      end
c     
c............................................................
      SUBROUTINE PICKUP ( blr, bli )
c............................................................

c.... Calculates the fields at the pickup loops and at the shell.
c     Called from VACUUM_MA.F with argulents ( bnlr, bnli )

      USE vacuum11_mod
      USE vacuum13_mod

      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum5.inc'
     
      REAL, DIMENSION(:,:), ALLOCATABLE ::
     $     chir, chii, cwrkr, cwrki
      REAL, DIMENSION(:), ALLOCATABLE ::
     $     zchipr, zchipi, chirr,
     $     bxr, bxi, bzr, bzi,
     $     btr, bti, bpr, bpi,
     $     bphir, bphii, bphi,
     $     phbt, phbp, phbx, phbz

      DIMENSION xpla(nths), zpla(nths),
     $     btxr(nths),btxi(nths),btzr(nths),btzi(nths),
     $     bthr(nths),bthi(nths),btpr(nths),btpi(nths),
     $     wrk1(nths), wrk2(nths)

      DIMENSION grdgre(nths2,nths2),grwp(nths,nths),
     .     grri(nths2,nfm2)
      DIMENSION blr(*), bli(*)
c     
      INTEGER, DIMENSION(:), ALLOCATABLE :: igdl

      common / bigv1 / grdgre, grwp, grri
c     
c$$$      common / gridc / xobs(ndimlp),zobs(ndimlp)
c     
      equivalence (xpla,xinf), (zpla,zinf)
c     
c$$$      parameter(neqv1=(nths-1)/2+1, neqv2=nths+1, neqv3=3*(nths+1)/2  )
c$$$c     
c$$$c$$$      equivalence ( grdgre(1,1), chii(1,1) )
c$$$c$$$      equivalence ( grdgre(neqv2,neqv1), chir(1,1) )

c$$$      equivalence ( grdgre(1,1), chii(1,1) )
c$$$      equivalence ( grdgre(nptg2x,nptg12), chir(1,1) )

      CHARACTER(130), DIMENSION(10) :: string

c     
c     ntloop is the number of loop positions  on thw shell.
c     
c     
c..   New variable, lfarw for SPARK2 to reset farwal.
c     
c$$$  lfarw = 0
c$$$  if ( a .gt. 10.0 ) lfarw = 1
c$$$  write ( iotty, '("a = ", 1pe12.3, " lfarw in PICKUP = ", i4 )' )
c$$$  $     a, lfarw
c$$$  write ( outmod, '("a = ", 1pe12.3, " lfarw in PICKUP = ", i4 )' )
c$$$  $     a, lfarw
c     
      ndlp = mth / ntloop
c     
      jmax1 = lmax(1) - lmin(1) + 1
      mth12 = 2*mth
      if ( lfarw .gt. 0 ) mth12 = mth
c     
c.....set up window.
c     
      CALL bounds(xpla,zpla,1,mth,xmnp,xmxp,zmnp,zmxp)
      xmin = xmnp
      xmax = xmxp
      zmin = zmnp
      zmax = zmxp
c     
      plrad = 0.5 * ( xmxp - xmnp )
      xmaj = 0.5 * ( xmxp + xmnp )

c     
      if ( lfarw .gt. 0 ) go to 20
c     
      call bounds(xwal,zwal,1,mw,xmnw,xmxw,zmnw,zmxw)
c     
      xmin = amin1(xmnp,xmnw)
      xmax = amax1(xmxp,xmxw)
      zmin = amin1(zmnp,zmnw)
      zmax = amax1(zmxp,zmxw)
c     
 20   CONTINUE

c$$$      ALLOCATE ( xloop(ndimlp1), zloop(ndimlp1) ) 
c$$$
c$$$      xloop = 0.0
c$$$      zloop = 0.0

      CALL loops

      nobs = nloop + 3*nloopr
c       mdimlp1 must accommodate the tangential field on the shell.
      ndimlp1 = nobs + ntloop + 5  ! 08/25/2009

      ALLOCATE (
     $     chir(5,ndimlp1),chii(5,ndimlp1),
     $     zchipr(ndimlp1),zchipi(ndimlp1),chirr(ndimlp1),
     $     cwrkr(5,ndimlp1),cwrki(5,ndimlp1),
     $     bxr(ndimlp1),bxi(ndimlp1),bzr(ndimlp1),bzi(ndimlp1),
     $     btr(ndimlp1),bti(ndimlp1), bpr(ndimlp1),bpi(ndimlp1),
     $     bphir(ndimlp1),bphii(ndimlp1),bphi(ndimlp1)
     $     )

c... The following needs larger space since they are also needed to store 
c     the fields at mth1 surface points. Should be ndimlp1 + mth1.
      
      ndimlp2 = ndimlp1 + mth1
      ALLOCATE (
     $     phbt(ndimlp2), phbp(ndimlp2), phbx(ndimlp2), phbz(ndimlp2)
     $     )
      bxr = 0.0
      bxi = 0.0
      bzr = 0.0
      bzi = 0.0
      bphir = 0.0
      bphii = 0.0
      bphi = 0.0
      chirr = 0.0
      btr = 0.0
      bti = 0.0
      bpr = 0.0
      bpi = 0.0

      phbt = 0.0
      phbp = 0.0
      phbx = 0.0
      phbz = 0.0
      btxr = 0.0
      btxi = 0.0
      btzr = 0.0
      btzi = 0.0
      bthr = 0.0
      bthi = 0.0
      btpr = 0.0
      btpi = 0.0
      zchipr = 0.0
      zchipi = 0.0

      CALL vacuum11_alloc ( ndimlp1 )  ! xobs, zobs

      call bounds ( xloop,zloop,1, nobs, xmn,xmx, zmn,zmx )
      xmin = amin1 ( xmin,xmn )
      xmax = amax1 ( xmax,xmx )
      zmin = amin1 ( zmin,zmn )
      zmax = amax1 ( zmax,zmx )
c     
      dtpw = dth
      ns = mth
      delx = plrad * deloop
      delz = plrad * deloop
c     
c.....get array of observer points, and initialize chi matrix.
c     
      chir = 0.0
      chii = 0.0
c     
      DO NSEW = 1, 5                !!! To about line 2901

         DO i = 1, nobs
c     
            chir(nsew,i) = 0.0
            chii(nsew,i) = 0.0
            cwrkr(nsew,i) = 0.0
            cwrki(nsew,i) = 0.0
c     
            go to ( 51, 52, 53, 54, 55), nsew
c     
 51         continue
c......north
c     
            xobs(i) = xloop(i)
            zobs(i) = zloop(i) + delz
            go to 90
c     
 52         continue
c.....south
c     
            xobs(i) = xloop(i)
            zobs(i) = zloop(i) - delz
            go to 90
c     
 53         continue
c.....east
c     
            xobs(i) = xloop(i) + delx
            zobs(i) = zloop(i)
            go to 90
c     
 54         continue
c.....west
c     
            xobs(i) = xloop(i) - delx
            zobs(i) = zloop(i)
            GO TO 90

 55         CONTINUE
c.....center
            xobs(i) = xloop(i) 
            zobs(i) = zloop(i)
            
 90         CONTINUE
         END DO                                   ! nobs do loop
c     
c.....pass chi on source surface.
c     
c.....plasma contribution.
c     
         do l1 = 1, jmax1
            do i = 1, mth
               chiwc(i,l1) = grri(i,l1)
               chiws(i,l1) = grri(i,jmax1+l1)
            end do
         end do
c     
c$$$         isg = 1
c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c... change sign to -1 since it seems to match the cylinder.
c    change the sign for the wall contribution also
c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
         isg = -1
         write ( outmod,'("isgn for plasma contribution = ", i3)')
     $        isg
         write ( iotty,'("isgn for plasma contribution = ", i3)')
     $        isg
c     
c$$$         call chi ( xpla,zpla,xplap,zplap,isg, chiwc,chiws, nobs,ns,1,
c$$$     $        cwrkr,cwrki,nsew, blr,bli )     

         call chi ( xpla,zpla,xplap,zplap,isg, chiwc,chiws, ns,1,
     $        cwrkr,cwrki,ndimlp1, nsew, blr,bli )     
      call atpoint("After Chi","lfarw",lfarw,"delx",delx,
     $     iotty, outmod )

c     
c..   Store away the plasma result.
c     
         do  i = 1, nobs
            chir(nsew,i) = cwrkr(nsew,i)
            chii(nsew,i) = cwrki(nsew,i)
         end do
c     
         if ( lfarw .gt. 0 ) go to 200
c     
c.....conductor contribution.
c     
         do l1 = 1, jmax1
            do i = 1, mw
               chiwc(i,l1) = grri(mth+i,l1)
               chiws(i,l1) = grri(mth+i,jmax1+l1)
            end do
         end do
c     
         do i = 1, nobs
            cwrkr(nsew,i) = 0.0
            cwrki(nsew,i) = 0.0
         end do
c     
c$$$         isg = -1 
c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c... change sign to +1 since it seems to match the cylinder.
c    change the sign for the plasma contribution also
c   $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

         isg = 1
c
         write ( outmod,
     $        '("isgn for conductor contribution = ", i3)') isg
c     
c$$$         call chi ( xwal,zwal,xwalp,zwalp,isg,chiwc,chiws, nobs,ns,0,
c$$$     $        cwrkr,cwrki,nsew, blr,bli )

         call chi ( xwal,zwal,xwalp,zwalp,isg,chiwc,chiws, ns,0,
     $        cwrkr,cwrki,ndimlp1, nsew, blr,bli )
c     
         do  i = 1, nobs
            chir(nsew,i) = chir(nsew,i) + cwrkr(nsew,i)
            chii(nsew,i) = chii(nsew,i) + cwrki(nsew,i)
         end do
c     
 200     CONTINUE
      END DO                                ! NSEW do loop
c     
      WRITE ( OUTMOD, '(/,36x,"****   Mirnov Loops   ****")' )
c$$$      WRITE ( OUTMOD,
c$$$     $     '(1x,"iobs",24x,"chir(nsew)",41x,"chii(nsew)" )' )
c$$$      DO i = 1, nobs
c$$$         WRITE (OUTMOD,'(1x, i4, 1p5e11.3, 2x, 1p5e11.3)') i,
c$$$     $        (chir(ii,i),ii=1,4), (chii(ii,i), ii= 1, 5)
c$$$      END DO
c     
c.....calculate b_x, b_z and b_\phi
c   
c..  At this point, the chi's in the whole region are calculated, properly
c    approximating the vanishing ones in the appropriate interior or
c    exterior region. The following coding attempts to actually not 
c    calculate the B fields in the latter, unwanted region.
c    It also Avoid points near the plasma surface. But calculate these fields
c     at those points on the surface.
c.... isgchi is sign for quantities involving Chi. This may be -1 for all 
c     calculations involving Chi. i.e., B_theta and B_phi. 

      ALLOCATE ( igdl(nobs) )

      igdl = 8
      isgchi = -1

c.. DO loop over all observers:
      DO i = 1, nobs               !! To about line 3044

         IF ( lchkinextp == 0 ) GO TO 238

c.. Check whether observers are interior or interior:
         fintjj = 0.0
         DO jj  = 1, mth
            dxjj = xloop(i) - xinf(jj)
            dzjj = zloop(i) - zinf(jj)
            rhoj2 = dxjj**2 + dzjj**2
c.. If rhoj2 too small, skip test.
            IF ( rhoj2 < 1.0e-16 ) GO TO 238
            fintjj = fintjj +
     $           ( zplap(jj)*dxjj - xplap(jj)*dzjj ) / rhoj2
         END DO
         fintjj = fintjj / mth
c$$$         WRITE ( OUTMOD, '("I, fintjj = ", I4, ES15.5)' )
c$$$     $        i, fintjj

c..   Neglect points interior to the plasma, if exterior solution required:
         IF ( fintjj > 0.1 ) THEN ! interior
            igdl(i) = 1
            IF ( linterior == 0 ) GO TO 239 
         END IF

c..   Neglect points exterior to the plasma, if interior solution required:
         IF ( fintjj < 0.1 ) THEN ! exterior
            igdl(i) = 0
            IF ( linterior == 1 ) GO TO 239
         END IF
            
 238     CONTINUE

c... Only the closed wanted region remains here. Treat points close to 
c    the boundary differently.

         IF ( lclosepts == 0 ) GO TO 2385

         DO j = 1, mth
            IF ( (xinf(j)-xloop(i))**2 + (zinf(j)-zloop(i))**2 
     $           < (epslp*plrad)**2 )  THEN ! points close to surface.
               igdl(i) = -1
               IF ( nloop < 1001 ) WRITE ( IOTTY,
     $              '("CLOSE: iobs,jpla, xinf, xloop, zinf, zloop = ",
     $              2I5, 4ES12.4)' )
     $              i,j, xinf(j),xloop(i) ,zinf(j),zloop(i)
               IF ( nloop < 1001 ) WRITE ( OUTMOD,
     $           '("CLOSE: iobs,jpla, xinf, xloop, zinf, zloop = ",
     $           2I5, 4ES12.4)' )
     $           i,j, xinf(j),xloop(i) ,zinf(j),zloop(i)
c..  Observer for these point is on the plasma surface
            xloop(i) = xinf(j)
            zloop(i) = zinf(j)
c.. Now calculate the fields at i for these points
c.. bnkr, bnki were angular fluxes. Divide by gpsjp to get actual field.
c.. bthpr, bthpi calculated at line ~625 in SUBROUTINE diaplt.
    
c... ireg takes care of interior and exterior cases:
c... don't overwrite isgchi by itself
            ireg = 1 - 2*linterior
            isgchix = ireg*isgchi
            alph = atan2m ( xinf(j+1)-xinf(j-1), zinf(j-1)-zinf(j+1) )
            cosalph = COS(alph)
            sinalph = SIN(alph) 
            bxr(i) = bnkr(j)*cosalph/gpsjp(j) + isgchix*bthpr(j)*sinalph
            bxi(i) = bnki(j)*cosalph/gpsjp(j) + isgchix*bthpi(j)*sinalph
            bzr(i) = bnkr(j)*sinalph/gpsjp(j) - isgchix*bthpr(j)*cosalph
            bzi(i) = bnki(j)*sinalph/gpsjp(j) - isgchix*bthpi(j)*cosalph
            zchipr(i) = isgchix*chipr(j)
            zchipi(i) = isgchix*chipi(j)
            bphir(i) =   n * zchipi(i) / xloop(i)
            bphii(i) = - n * zchipr(i) / xloop(i)

             IF ( nloop < 1001 ) WRITE ( IOTTY,
     $      '("i,chipr,chipi,bnkr,bnki,bthpr,bthpi,bphir,bphii,xloop=",/
     $           i4, 9ES13.5)' ) i, chipr(j),chipi(j),
     $           bnkr(j)/gpsjp(j),bnki(j)/gpsjp(j),
     $           bthpr(j),bthpi(j),bphir(i), bphii(i),xloop(i)
             IF ( nloop < 1001 ) WRITE ( OUTMOD,
     $      '("i,chipr,chipi,bnkr,bnki,bthpr,bthpi,bphir,bphii,xloop=",/
     $           i4, 9ES13.5)' ) i, chipr(j),chipi(j),
     $           bnkr(j)/gpsjp(j),bnki(j)/gpsjp(j),
     $           bthpr(i),bthpi(i),bphir(i), bphii(i),xloop(i)
            GO TO 239      ! Other j  points may be closer but refine later.
         END IF  ! Points very close to surface
      END DO     ! mth, over boundary points.

 2385 CONTINUE      !!! lclosepts

c.. The following are the wanted points not on the surface:


         bxr(i) = ( chir(3,i) - chir(4,i) ) / (2.0*delx)
         bxi(i) = ( chii(3,i) - chii(4,i) ) / (2.0*delx)
         bzr(i) = ( chir(1,i) - chir(2,i) ) / (2.0*delz)
         bzi(i) = ( chii(1,i) - chii(2,i) ) / (2.0*delz)
         zchipr(i) = chir(5,i)
         zchipi(i) = chii(5,i)
         bphir(i) =   n * chii(5,i) / xloop(i)
         bphii(i) = - n * chir(5,i) / xloop(i)
c     
         phbx(i) = 0.0
         phbz(i) = 0.0
c     
c$$$  if ( i .eq. nloop+nloopr+1 ) go to 234
c$$$  if ( i .eq. nloop+2*nloopr+1 ) go to 234
c$$$  if ( i .gt. nloop+1 ) go to 235
c$$$  c     
c$$$  234     continue
         phbx(i) = ATAN2M ( bxi(i),bxr(i) )
         phbz(i) = ATAN2M ( bzi(i),bzr(i) )
         IF ( phbx(i) < 0.0 ) phbx(i) = phbx(i) + twopi
         IF ( phbz(i) < 0.0 ) phbz(i) = phbz(i) + twopi

 239     CONTINUE
 
      END DO     ! nobs over all observers.

      WRITE ( OUTMOD, '(/, 2x, 12(a10) )' ) "xloop", "zloop",
     $     "bxr", "bxi", "bzr", "bzi",
     $     "zchipr", "zchipi", "bphir", "bphii", "phbx", "phbz"
      DO i = 1, nobs
         WRITE ( OUTMOD, '(1x, I3, 12ES10.2)' ) i, xloop(i), zloop(i),
     $        bxr(i), bxi(i),
     $        bzr(i), bzi(i), zchipr(i), zchipi(i),
     $        bphir(i), bphii(i), phbx(i), phbz(i)
      END DO

c...  Write results to file:

      IF ( iloop /= 5 ) GO TO 240     !  08/24/2009

      iobv = 76
      OPEN (  iobv, FILE='B-VACUUM', STATUS='REPLACE',
     $     FORM='FORMATTED' )

      WRITE(iobv,'(1x,2(a4,I4))') "nr:",nxlpin, "nz:",nzlpin
      WRITE(iobv,'(1x,a2,8(a16))')"l","r","z","re(b_r)","im(b_r)",
     $     "re(b_z)","im(b_z)","re(b_phi)","im(b_phi)"

      DO i = 1, nxlpin
         DO j = 1, nzlpin
            indxl = (i-1)*nzlpin + j
            WRITE ( iobv, '(1x, I2, 8(E16.8) )' )
     $           igdl(indxl),
     $           xloop(indxl), zloop(indxl),
     $           bxr(indxl), bxi(indxl),
     $           bzr(indxl), bzi(indxl),
     $           bphir(indxl), bphii(indxl)
         END DO
      END DO

      CLOSE ( UNIT = iobv )

 240  CONTINUE

      DEALLOCATE ( igdl )

      IF ( iloop /= 6 ) GO TO 260     !  08/31/2010

      IF ( irepeat == 0 ) THEN

         iobv = 76
         OPEN ( iobv, FILE='B-VAC-ITER', STATUS='REPLACE',
     $     FORM='FORMATTED' )
         
         WRITE ( IOBV,
     $        '(5x, "Magnetic Fields at Saddle Loops Calculated on:",
     $        2x, A )' ) date_array
         WRITE ( IOBV,
     $        '(/, 5x, "Job id:", 2x, A,
     $        "Input Control File:", 2x, A)' ) jobid, filein2
         
         WRITE ( IOBV, '(/,1x, I5,  "  No of internal points." )' )
     $        mth1
         WRITE ( IOBV, '(/,1x, I5, "  Min Fourier mode." )' )
     $        lmin(1)
         WRITE ( IOBV, '(  1x, I5, "  Max Fourier mode.", / )' )
     $        lmax(1)
         WRITE ( IOBV, '(  1x, E16.8, "  Plasma ~Radius." )' )
     $        plrad
         WRITE ( IOBV, '(  1x, E16.8,
     $        "  Fraction of Plasma Radius used for derivatives.", / )')
     $        deloop
         WRITE ( IOBV, '(  1x, I5,  "  No of toroidal modes.", / )' )
     $        nrepeat
         WRITE ( IOBV, '(  1x, I5,  "  nmin.", / )' ) nmin
         WRITE ( IOBV, '(  1x, I5,  "  nmax.", / )' ) nmax

      END IF


      nztor = n
      WRITE ( IOBV, '(1x, I4, A10)' )
     $     nztor, "   n_tor"
      WRITE ( IOBV, '(1x, I5, A10 )' )
     $     nloop, "   n_loops"

      WRITE ( IOBV, '(4x, A2, A12, 7(A16))' )
     $     "i", "r  ", "z  ", "re(b_r)", "im(b_r)",
     $     "re(b_z)", "im(b_z)", "re(b_phi)", "im(b_phi)"

      DO indxl = 1, nloop
         WRITE ( iobv, '(1x, I5, 8(E16.8) )' )
     $        indxl,
     $        xloop(indxl), zloop(indxl),
     $        bxr(indxl), bxi(indxl),
     $        bzr(indxl), bzi(indxl),
     $        bphir(indxl), bphii(indxl)
      END DO

      IF ( irepeat == nrepeat - 1 ) CLOSE ( UNIT = iobv )

 260  CONTINUE

      nobs0 = nobs          ! = nloop + 3* nloopr

      IF ( lfarw > 0 ) GO TO 241       

c     Include the tangential field on the wall. nobs will be incremented as 
c     needed.
c     
c  Subroutine BTANG will calculate on all points on the shell.


      CALL btang ( blr,bli, bthr,bthi, btxr,btxi, btzr,btzi,
     $     btpr,btpi, wrk1, wrk2 )

c... Multiply by isghchi to be consistent wiith the free space loops.
c       see above where isgchi is defined.  08/26/2009

      bthr = isgchi * bthr
      bthi = isgchi * bthi
      btxr = isgchi * btxr
      btxi = isgchi * btxi
      btzr = isgchi * btzr
      btzi = isgchi * btzi
      btpr = isgchi * btpr
      btpi = isgchi * btpi

c... Write out these components to SVD-IN_OUT:


      WRITE ( iosvd, '(5x, /a/ )' )  "BTHR, real theta component:"
      WRITE ( iosvd, '(10es15.7)' ) (bthr(i), i = 1, mth1)

      WRITE ( iosvd, '(5x, /a/ )' )  "BTHI, imag theta component:"
      WRITE ( iosvd, '(10es15.7)' ) (bthi(i), i = 1, mth1)

      WRITE ( iosvd, '(5x, /a/ )' )  "BTPR, real phi component:"
      WRITE ( iosvd, '(10es15.7)' ) (btpr(i), i = 1, mth1)

      WRITE ( iosvd, '(5x, /a/ )' )  "BTPI, imag phi component:"
      WRITE ( iosvd, '(10es15.7)' ) (btpi(i), i = 1, mth1)

      CLOSE ( UNIT = 86 )


c...   Write out B_tangential  for all  points

      WRITE ( OUTMOD, '(15x, "Tangential Field on Shell",/,
     $     3x, "I", 5x, "xwall", 8x, "zwall", 8x, "Bthr", 8x,
     $     "Bthi" )' )
      DO i = 1, mth
         WRITE( OUTMOD, '(i4, 4es12.4)' ) i, xwal(i), zwal(i),
     $        bthr(i), bthi(i)
      END DO
c     
c...  Plot phases of the surface loops: 
c     X and Z components
c     Theta and phi components 
c     Magnitude of B_th
c     Complex B_th
c     
c... The storage for the following were increased to accommodate mth fields.
c      See the allocate statement above.     ! 08/25/2009
      DO i = 1, mth1
         phbt(i) = atan2m ( bthi(i),bthr(i) )
         phbp(i) = atan2m ( btpi(i),btpr(i) )
         phbx(i) = atan2m ( btxi(i),btxr(i) )
         phbz(i) = atan2m ( btzi(i),btzr(i) )
c$$$         if ( phbt(i) .lt. 0.0 ) phbt(i) = phbt(i) + twopi
c$$$         if ( phbp(i) .lt. 0.0 ) phbp(i) = phbp(i) + twopi
c$$$         if ( phbx(i) .lt. 0.0 ) phbx(i) = phbx(i) + twopi
c$$$         if ( phbz(i) .lt. 0.0 ) phbz(i) = phbz(i) + twopi
      END DO
c     
      call pospl2 ( piwgrd,bthr,bthi, 1,mth1,1,
     $     "Bthr and Bthi", mth1)
      call pospl1(bthr,bthi, 1,mth1,2,'bthi vs. bthr', mth1)
      call posplv2 ( piwgrd,phbt,phbp, wrk1,wrk2, 1,mth1,3,
     $     "Phase: Bth & Bphi", mth1)
      call posplv2 ( angwalq,phbt,phbp, wrk1,wrk2, 1,mth1,4,
     $     "Phase vs Angle", mth1)
c      call posplv2 ( piwgrd,phbx,phbz, wrk1,wrk2, 1,mth1,1,
c     $     "Phase: Bt_x,Bt_z", mth1)

      CALL MAP ( 0.0,1.0, 0.0,1.0, 0.0,1.0, 0.0,1.0 )
      xw = 0.4
      yw = 0.2
      call setlch ( xw, yw, 0,1,0,-1 )
      write ( string,'(a)' ) jobid
      call wrtstr ( string, 1 )
      
      CALL FRAMEP(JOBID,FF)

c...  Now do ONLY NTLOOPS for plotting on shell and writing.

      DO i = 1, mth, ndlp
         nobs = nobs + 1                  ! This needed ntloop more storage
         xloop(nobs) = xwalq(i)
         zloop(nobs) = zwalq(i)
         bxr(nobs) = btxr(i)
         bxi(nobs) = btxi(i)
         bzr(nobs) = btzr(i)
         bzi(nobs) = btzi(i)
         btr(nobs) = bthr(i)
         bti(nobs) = bthi(i)
         bpr(nobs) = btpr(i)
         bpi(nobs) = btpi(i)
      END DO

 241  CONTINUE                   ! IF farwal

c  Get number of free space loops, exact number of ntloops, 
c   and starting index for ntloops:

      nfobs = nobs0
      ntobs = nobs - nobs0
      nfobs1 = nobs0 + 1
c     
      nobs = nobs0

      IF ( lfarw > 0 ) GO TO 243
c     
      do i = 1, mth, ndlp
         nobs = nobs + 1
         phbt(nobs) = atan2m ( bti(nobs),btr(nobs) )
         phbp(nobs) = atan2m ( bpi(nobs),bpr(nobs) )
         phbx(nobs) = atan2m ( bxi(nobs),bxr(nobs) )
         phbz(nobs) = atan2m ( bzi(nobs),bzr(nobs) )
c$$$         if ( phbt(nobs) .lt. 0.0 ) phbt(nobs) = phbt(nobs) + twopi
c$$$         if ( phbp(nobs) .lt. 0.0 ) phbp(nobs) = phbp(nobs) + twopi
c$$$         if ( phbx(nobs) .lt. 0.0 ) phbx(nobs) = phbx(nobs) + twopi
c$$$         if ( phbz(nobs) .lt. 0.0 ) phbz(nobs) = phbz(nobs) + twopi
      end do

 243  CONTINUE              !  IF farwal
c     
      WRITE ( OUTMOD, 210 )
 210  FORMAT ( /,4x,"i",6x,"x",7x,"z",7x,"bxr",8x,"bxi",8x,
     $     "bzr",8x,"bzi",8x,"bpr",8x,"bpi",4x,"phbt",4x,
     $     "phbx,mu",2x,"phbz,mu" )
c     
      DO i = 1, nobs0
         WRITE ( OUTMOD, 245 ) i, xloop(i),zloop(i), bxr(i),bxi(i),
     .        bzr(i),bzi(i), bpr(i),bpi(i), phbt(i),phbx(i),phbz(i)
 245     FORMAT ( 1x, i4, 1x, 2f8.3, 1p6e11.3, 0p3f8.3 )
      END DO

      IF (  lfarw > 0 ) GO TO 246

      WRITE ( OUTMOD, '(20x,"field on shell")' )
      DO i = nobs0+1, nobs
         WRITE ( OUTMOD, 245 ) i, xloop(i),zloop(i), bxr(i),bxi(i),
     .        bzr(i),bzi(i), bpr(i),bpi(i), phbt(i),phbx(i),phbz(i)
      END DO

c     
 246  CONTINUE

      if ( nloopr == 0 ) go to 2391
c     
c........calculate power of algebraic fall-off.  store in phbx and phbz.
c........Take care of situation where xloop < x_major-axis.
c     
      WRITE ( OUTMOD, '(/, 4x, "I", 2x,"X-Falloff",4x,"Z-Falloff" )' )
      DO i = nloop+2, nloop + 3*nloopr
         phbx(i) = 0.0
         phbz(i) = 0.0
         xldif = abs ( xloop(i) - xloop(i-1) )
         if ( (i .eq. nloop + nloopr+1) .or.
     $        (i .eq. nloop+2*nloopr+1) ) go to 2390
         if ( i .gt. (nloop+2*nloopr) ) go to 2380
         if ( (i .gt. nloop+nloopr) .and. (i .le. nloop+2*nloopr)
     $        .and. (xldif .le. 1.e-6) ) then
            phbx(i) = -1.111
            phbz(i) = -1.111
            go to 2390
         end if
         phbx(i) = alog ( abs(bxr(i)/bxr(i-1)) )
     .        / alog ( (xloop(i)-xmaj)/(xloop(i-1)-xmaj) )
         GO TO 2390
 2380    CONTINUE
         phbz(i) = ALOG( ABS(bzi(i)/bzi(i-1)) )
     $        / alog( zloop(i)/zloop(i-1) )
c     
 2390    continue
         write ( outmod, '(1x, i4, 1p2e12.4)' ) i, phbx(i), phbz(i)
      END DO

 2391 continue
c     
c......plot b vector for different phi's.  use xobs and zobs for storage.
c....  Increment phi by dphi/n 
c     
      dphi = pye / nphil

      IF ( abs(n) > 0.001 ) THEN
         dphi = dphi / n
      END IF

      bmax = -1.0e10
c     
      DO iph = 1, nphil
c     
         phi = (iph-1) * dphi
         phi_deg = 180.0 * phi / pye
         sinnph = sin ( n*phi )
         cosnph = cos( n*phi )
c     
         DO i = 1, nobs
            chirr(i) = zchipr(i)*cosnph + zchipi(i)*sinnph
            xobs(i) = bxr(i)*cosnph + bxi(i)*sinnph
            zobs(i) = bzr(i)*cosnph + bzi(i)*sinnph
            bphi(i) = bphir(i)*cosnph + bphii(i)*sinnph
            bmod = SQRT ( xobs(i)**2 + zobs(i)**2 )
            IF ( iph > 1 ) GO TO 250
            IF ( bmax < bmod ) bmax = bmod
 250        CONTINUE
         END DO
c     
         dx = xmax - xmin
         dz = zmax - zmin
c         dz = 2.0 * zmax
         dxz = amax1(dx,dz)
         xmaxw = 0.5 * ( xmax + xmin ) + 0.50*dxz
         xminw = 0.5 * ( xmax + xmin ) - 0.50*dxz
         zmaxw = 0.5 * ( zmax + zmin ) + 0.50*dxz
         zminw = 0.5 * ( zmax + zmin ) - 0.50*dxz
c$$$         zmaxw = 0.50 * dxz
c$$$         zminw = - zmaxw
c     
         xming = xminw - 0.05*(xmaxw-xminw)
         xmaxg = xmaxw + 0.05*(xmaxw-xminw)
         zming = zminw - 0.05*(zmaxw-zminw)
         zmaxg = zmaxw + 0.05*(zmaxw-zminw)
c     
         CALL MAPS(xming,xmaxg,zming,zmaxg,0.1,0.815,0.285,1.0)
c     
         CALL SETPCH ( 0, 0, -1, 1 )
         CALL POINTC ('+', xloop,zloop,nobs, -1,-1, 0.0, 0.0 )
c     
         idelg = delg
         dartl = FLOAT(idelg)/10.0
         darth0 = delg - 10.0*dartl
         WRITE ( OUTMOD, '("dxz, bmax, delg, darth0 = ",
     $        1p4e13.5)' ) dxz, bmax, delg, darth0
         WRITE ( IOTTY , '("dxz, bmax, delg, darth0 = ",
     $        1p4e13.5)' ) dxz, bmax, delg, darth0

         n_int = n
         WRITE ( OUTMOD, '(/, "n_int = ",I2, 1x, "iph = ",I2, 1x,
     $        "phi_deg = ", ES11.4,/ )' ) n_int, iph, phi_deg

         WRITE ( IOTTY, '(/, "n_int = ",I2, 1x, "iph = ",I2, 1x,
     $        "phi_deg = ", ES11.4,/ )' ) n_int, iph, phi_deg

c.... Calculate B vectors for plotting. Store in VX and VZ,
         WRITE ( IOTTY,
     $        '(1x," Plot vectors are in Vx and Vz: ")' )
        WRITE ( OUTMOD,
     $        '(1x," Plot vectors are in Vx and Vz: ")' )
        WRITE ( IOTTY, '(1x, "Iob", 6x, "X", 12x, "Z", 11x, "Vx", 12x,
     $       "Vz", 11x, "Chi", 10x, "B_phi", 9x, "Bx", 11x, "B_z" )' )
        WRITE ( OUTMOD, '(1x, "Iob", 6x, "X", 12x, "Z", 11x, "Vx", 12x,
     $       "Vz", 11x, "Chi", 10x, "B_phi", 9x, "Bx", 11x, "B_z" )' )

c$$$        rdfac = 1.0
        rdfac = 0.15
         DO i = 1, nobs
            x0 = xloop(i)
            z0 = zloop(i)
            vx = rdfac * dxz * xobs(i) / (2.0*bmax)
            vz = rdfac * dxz * zobs(i) / (2.0*bmax)
c     darth = .01
c$$$            IF (iph == 1 ) WRITE ( IOTTY,  '(i4,1p8e13.5 )' )
c$$$            WRITE ( IOTTY,  '(i4,1p8e13.5 )' )
c$$$     $           i, x0, z0, vx, vz, chirr(i),bphi(i),xobs(i), zobs(i)
c$$$            IF (iph == 1 ) WRITE ( OUTMOD, '(i4,1p8e13.5 )' )

c....   Only write this if nloop is less than 1001:
            IF ( nloop < 1001) 
     $           WRITE ( OUTMOD, '(i4,1p8e13.5 )' )
     $           i, x0, z0, vx, vz, chirr(i),bphi(i),xobs(i), zobs(i)
            CALL dart ( x0,z0, vx,vz, darth0 )
         END DO
        WRITE ( IOTTY, '(1x, "Iob", 6x, "X", 12x, "Z", 11x, "Vx", 12x,
     $       "Vz", 11x, "Chi", 10x, "B_phi", 9x, "Bx", 11x, "B_z" )' )
        WRITE ( OUTMOD, '(1x, "Iob", 6x, "X", 12x, "Z", 11x, "Vx", 12x,
     $       "Vz", 11x, "Chi", 10x, "B_phi", 9x, "Bx", 11x, "B_z" )' )

c.... plot conductors and plasma shape..
c     
         CALL points ( xpla,zpla, mth, -1,-1,0.,0. )
         IF ( lfarw == 0 )
     .        CALL points ( xwal,zwal, mth, -1,-1,0.,0. )

      CALL MAP ( 0.0,1.0, 0.0,1.0, 0.0,1.0, 0.0,1.0 )
      xw = 0.1
      yw = 0.2
      n_int = n
      CALL SETLCH ( xw, yw, 0,1,0,-1 )
      WRITE ( STRING,'(a, 2x, "n_tor = ",I2,/, "nloop = ",I3, 1x,
     $     "aloop = ",ES11.4, 1x, "phi_deg = ", ES11.4)' )
     $     jobid, n_int, nloop, aloop, phi_deg
      CALL WRTSTR ( STRING, 2 )

      CALL FRAMEP(JOBID,FF)
c     
      END DO

c... Don't plot \phi variation if too many loops. This is needed since 
c    chiwc and chiws, being used as dummies here are only dimensioned to nths.

      IF ( nloop >= nths-100 ) GO TO 300

c.....Plot each channel on shell as function of phi.   
c     Use chiwc, chiws as storage.

c   Plot X and Z components first.

      ncyc = 51
      IF ( ndimlp1 < 51 ) ncyc = ndimlp1
      cycl = 2.0
      dphi = cycl*twopi / (ncyc-1)
      bmax = -1.0e10
c     
      DO iph = 1, ncyc
         phi = (iph-1) * dphi
         sinnph = SIN(n*phi)
         cosnph = COS(n*phi)
         xobs(iph) = phi/twopi
         DO i = 1, nloop
            chiwc(iph,i) = bxr(i)*cosnph + bxi(i)*sinnph
            chiws(iph,i) = bzr(i)*cosnph + bzi(i)*sinnph
            bmod = SQRT ( chiwc(iph,i)**2+chiws(iph,i)**2 )
            IF ( bmax > bmod ) bmax = bmod
         END DO
      END DO
c     
      zz1 = nloop + 0.5
      call maps(0.0,cycl,0.0,zz1,0.1,.4,.285,1.0)
c     
      do i = 1, nloop
         do iph = 1, ncyc
            zobs(iph) = chiwc(iph,i)/ (2.0*bmax) + i
         end do
         call trace ( xobs, zobs, ncyc, 1,1,0.,0. )
      end do
     
      xw = 0.0
      yw = - 0.15*nloop
      call setlch ( xw, yw, 1,1,0,-1 )
      write ( string,
     $     '(a,/,"B_x and B_z of Free loops around Plasma")' ) jobid
      call wrtstr ( string, 2 )

      call maps(0.0,cycl,0.0,zz1,.5,.8,.285,1.0)
c     
      do i = 1, nloop
         do iph = 1, ncyc
            zobs(iph) = chiws(iph,i)/ (2.0*bmax) + i
         end do
         call trace ( xobs, zobs, ncyc, 1,1,0.,0. )
      end do

      call framep(jobid,ff)

c.. Now the free loops WITH CONSTANT AMPLITUDES

      do iph = 1, ncyc
         phi = (iph-1) * dphi
         sinnph = sin(n*phi)
         cosnph = cos(n*phi)
         xobs(iph) = phi/twopi
         do i = 1, nloop
            abx = sqrt ( bxr(i)**2 + bxi(i)**2 )
            abz = sqrt ( bzr(i)**2 + bzi(i)**2 )
            chiwc(iph,i) = (bxr(i)*cosnph + bxi(i)*sinnph)/(2.0*abx)
            chiws(iph,i) = (bzr(i)*cosnph + bzi(i)*sinnph)/(2.0*abz)
         end do
      end do
c     
      zz1 = nloop + 0.5
      call maps(0.0,cycl,0.0,zz1,0.1,.4,.285,1.0)
c     
      do i = 1, nloop
         do iph = 1, ncyc
            zobs(iph) = chiwc(iph,i) + i
         end do
         call trace ( xobs, zobs, ncyc, 1,1,0.,0. )
      end do
     
      xw = 0.0
      yw = - 0.15*nloop
      call setlch ( xw, yw, 1,1,0,-1 )
      write ( string,
     $     '(a,/, "B_x and B_z of Free loops around Plasma: "
     $     "constant amplitude")' ) jobid
      call wrtstr ( string, 2 )

      call maps(0.0,cycl,0.0,zz1,.5,.8,.285,1.0)
c     
      do i = 1, nloop
         do iph = 1, ncyc
            zobs(iph) = chiws(iph,i) + i
         end do
         call trace ( xobs, zobs, ncyc, 1,1,0.,0. )
      end do

      call framep(jobid,ff)

 300  CONTINUE     ! End of skipping plots of \phi variation.

      if ( nloopr == 0 ) go to 2302

c  Now the radial loops:

c     ncyc = 51
c     cycl = 2.0 * twopi
c     dphi = cycl / (ncyc-1)
      bmax = -1.0e10
c     
      do iph = 1, ncyc
         phi = (iph-1) * dphi
         sinnph = sin(n*phi)
         cosnph = cos(n*phi)
c     xobs(iph) = phi
         do i = nloop+1, nfobs
            chiwc(iph,i) = bxr(i)*cosnph + bxi(i)*sinnph
            chiws(iph,i) = bzr(i)*cosnph + bzi(i)*sinnph
            bmod = sqrt ( chiwc(iph,i)**2+chiws(iph,i)**2 )
            if ( bmax .lt. bmod ) bmax = bmod
         end do
      end do
c     
      zz1 = nloop + 3*nloopr + 0.5
      zz0 = nloop - 0.5
      call maps(0.0,cycl,zz0,zz1,0.1,.4,.285,1.0)
c     
      do i = nloop+1, nfobs
         do iph = 1, ncyc
            zobs(iph) = chiwc(iph,i)/ (2.0*bmax) + i
         end do
         call trace ( xobs, zobs, ncyc, 1,1,0.,0. )
      end do

      xw = 0.0
      yw = zz0 - 0.15*(nfobs-nloop)
      call setlch ( xw, yw, 1,1,0,-1 )
      write ( string,
     $     '(a,/,"B_x and B_z of Radial loops around Plasma")' ) jobid
      call wrtstr ( string, 2 )

      call maps(0.0,cycl,zz0,zz1,.5,.8,.285,1.0)
c     
      do i = nloop+1, nfobs
         do iph = 1, ncyc
            zobs(iph) = chiws(iph,i)/ (2.0*bmax) + i
         end do
         call trace ( xobs, zobs, ncyc, 1,1,0.,0. )
      end do
     
      call framep(jobid,ff)

 2302 continue

      IF ( lfarw > 0 ) GO TO 500

c  Now Plot the Theta and Phi components of shell fields.
c 
      bmax = -1.0e10

      do iph = 1, ncyc
         phi = (iph-1) * dphi
         sinnph = sin(n*phi)
         cosnph = cos(n*phi)
c         xobs(iph) = phi/twopi
         do i = nfobs1, nobs
            chiwc(iph,i) = btr(i)*cosnph + bti(i)*sinnph
            chiws(iph,i) = bpr(i)*cosnph + bpi(i)*sinnph
            bmod = sqrt ( chiwc(iph,i)**2+chiws(iph,i)**2 )
            if ( bmax .lt. bmod ) bmax = bmod
         end do
      end do
c     
      zz1 = nobs + 0.5
      zz0 = nfobs - 0.5
      call maps(0.0,cycl,zz0,zz1,0.1,.4,.285,1.0)
c     
      do i = nfobs1, nobs
         do iph = 1, ncyc
            zobs(iph) = chiwc(iph,i)/ (2.0*bmax) + i
         end do
         call trace ( xobs, zobs, ncyc, 1,1,0.,0. )
      end do

      xw = 0.0
      yw = zz0 - 0.15*(nobs-nfobs)
      call setlch ( xw, yw, 1,1,0,-1 )
      write ( string, '(a,/,"B_theta and B_phi on shell")' ) jobid
      call wrtstr ( string, 2 )

      call maps(0.0,cycl,zz0,zz1,.5,.8,.285,1.0)
c     
      do i = nfobs1, nobs
         do iph = 1, ncyc
            zobs(iph) = chiws(iph,i)/ (2.0*bmax) + i
         end do
         call trace ( xobs, zobs, ncyc, 1,1,0.,0. )
      end do

      call framep(jobid,ff)

c     c  Now Plot the Theta and Phi components, WITH EQUAL AMPLITUDES

      do iph = 1, ncyc
         phi = (iph-1) * dphi
         sinnph = sin(n*phi)
         cosnph = cos(n*phi)
c         xobs(iph) = phi/twopi
         do i = nfobs1, nobs
            abt = sqrt ( btr(i)**2 + bti(i)**2 )
            abp = sqrt ( bpr(i)**2 + bpi(i)**2 )
            chiwc(iph,i) = (btr(i)*cosnph + bti(i)*sinnph)/(2.0*abt)
            chiws(iph,i) = (bpr(i)*cosnph + bpi(i)*sinnph)/(2.0*abp)
         end do
      end do
c     
      zz1 = nobs + 0.5
      zz0 = nfobs - 0.5
      call maps(0.0,cycl,zz0,zz1,0.1,.4,.285,1.0)
c     
      do i = nfobs1, nobs
         do iph = 1, ncyc
            zobs(iph) = chiwc(iph,i) + i 
         end do
         call trace ( xobs, zobs, ncyc, 1,1,0.,0. )
      end do

      xw = 0.0
      yw = zz0 - 0.15*(nobs-nfobs)
      call setlch ( xw, yw, 1,1,0,-1 )
      write ( string,
     $     '(a,/,"B_theta and B_phi on shell: const. amplitude")' )
     $     jobid
      call wrtstr ( string, 2 )

      call maps(0.0,cycl,zz0,zz1,.5,.8,.285,1.0)
c     
      do i = nfobs1, nobs
         do iph = 1, ncyc
            zobs(iph) = chiws(iph,i) + i 
         end do
         call trace ( xobs, zobs, ncyc, 1,1,0.,0. )
      end do

      call framep(jobid,ff)

 500  CONTINUE
      
      DEALLOCATE(
     $     chir,chii,
     $     zchipr,zchipi,chirr,
     $     cwrkr,cwrki,
     $     bxr,bxi,bzr,bzi,
     $     btr,bti, bpr,bpi,
     $     bphir,bphii,bphi,
     $     phbt,phbp, phbx, phbz
     $     )

      CALL vacuum11_dealloc

      CALL vacuum13_dealloc

      RETURN
      END
c
c.............................................................
c
      subroutine btang ( blr,bli, bthr,bthi,btxr,btxi,btzr,btzi,
     $     btpr,btpi, wrk1,wrk2 )
c.............................................................

c Calculates the tangential magnetic field along the shell.
c            Do it in equal arcs.

c  B_theta in (BTHR,BTHI), B_phi in (BTPR,BTPI), 
c  B_x     in (BTXR,BTXI), B_z   in (BTZR,BTZI)

      include 'vacuum1.inc'
      include 'vacuum2.inc'
      include 'vacuum3.inc'
      include 'vacuum5.inc'
c
      dimension blr(*),bli(*), btxr(*),btxi(*), btzr(*),btzi(*),
     $     bthr(*),bthi(*), btpr(*),btpi(*), wrk1(*),wrk2(*)

      REAL, DIMENSION(:), ALLOCATABLE :: zthe, zchiwr, zchiwi,
     $     zchiwrq, zchiwiq

      ALLOCATE ( zthe(mth+5), zchiwr(mth+5), zchiwi(mth+5),
     $     zchiwrq(mth+5), zchiwiq(mth+5) )
c
      anq = n * qa1
      jmax1 = lmax(1) - lmin(1) + 1
c
      do i = 1, mth1
         wrk1(i) =  xwalpq(i)**2 + zwalpq(i)**2
         zthe(i) = (i-1) * dth
      end do
c      write ( iotty, '("wrk1 = ", 1p10e11.3)' ) (wrk1(i), i = 1, mth1)
c
      do i = 1, mth1
         zchiwr(i) = 0.0
         zchiwi(i) = 0.0
      end do
      do l1 = 1, jmax1
         do i = 1, mth1
            zbr = blr(l1)
            zbi = bli(l1)
            zchiwr(i) =  zchiwr(i) + ( cwallr(i,l1)*zbr
     $           - cwalli(i,l1)*zbi )
            zchiwi(i) =  zchiwi(i) + ( cwallr(i,l1)*zbi
     $           + cwalli(i,l1)*zbr )
         end do
      end do

c..Need to transform chiwall to equal arcs if not done already.

      zchiwrq(1:mth1) = zchiwr(1:mth1)
      zchiwiq(1:mth1) = zchiwi(1:mth1)

      IF ( leqarcw == 1 ) GO TO 120

      DO i = 1, mth1
         ttt = elwal(mth1)*(i-1)/mth
         CALL lag ( elwal, zchiwr, mth1,3,ttt,zchiwrq(i),dum, 0 )
         CALL lag ( elwal, zchiwi, mth1,3,ttt,zchiwiq(i),dum, 0 )
      END DO

 120  CONTINUE

      btpr(1:mth1) = + n * zchiwiq(1:mth1) / xwalq(1:mth1)
      btpi(1:mth1) = - n * zchiwrq(1:mth1) / xwalq(1:mth1)

      call difspl ( mth, zthe, zchiwrq, wrk2 )
c     
      do i = 1, mth1
         z12 = wrk2(i) / wrk1(i)
         bthr(i) = wrk2(i) / sqrt(wrk1(i))
         btxr(i) = xwalpq(i) * z12
         btzr(i) = zwalpq(i) * z12
      end do
c     
      call difspl ( mth, zthe, zchiwiq, wrk2 )
c     pq: Command not found.

      do i = 1, mth1
         z12 = wrk2(i) / wrk1(i)
         bthi(i) = wrk2(i) / sqrt(wrk1(i))
         btxi(i) = xwalpq(i) * z12
         btzi(i) = zwalpq(i) * z12
      end do

       DEALLOCATE ( zthe, zchiwr, zchiwi, zchiwrq, zchiwiq )
 
      return
      end
c     
c......................................................................
      subroutine chi(xsce,zsce,xscp,zscp,isg,creal,cimag,ns,ip,
     $     chir,chii,ndimlp1, nsew,blr,bli)
c......................................................................

      USE vacuum11_mod     ! xobs, zobs

      include 'vacuum1.inc'
      include 'vacuum3.inc'
      include 'vacuum2.inc'
      include 'vacuum5.inc'
c     
      dimension iop(2), ww1(nths),ww2(nths),ww3(nths),tab(3)
      dimension blr(*),bli(*), xsce(*),zsce(*),xscp(*),zscp(*)
      dimension creal(nths,nfm), cimag(nths,nfm)
      dimension chir(5,ndimlp1), chii(5,ndimlp1)
c     
      real nq
c     
c$$$      common / gridc / xobs(ndimlp),zobs(ndimlp)
c     
c     Factpi for dividing BVAL to be consistent with AVAL.
c     
      factpi = twopi
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
      ns1 = ns + 1
      nobs = nloop + 3*nloopr
c     
      do io = 1, nobs
c     
         xs = xobs(io)
         zs = zobs(io)
c     
c$$$      call atpoint("IN ioloop","io",io,"qa1",qa1,
c$$$     $     iotty, outmod )
c$$$      call atpoint("IN ioloop","is",is,"factpi",factpi,
c$$$     $     iotty, outmod )

         do is = 1, ns
            
            xt = xsce(is)
            zt = zsce(is)
            xtp = xscp(is)
            ztp = zscp(is)
c     
c$$$            IF ( (io == 17) .AND. (nsew == 1) ) THEN
c$$$               call atpoint("IN ioloop","io",io,"xs",xs,
c$$$     $              iotty, outmod )
c$$$               call atpoint("IN ioloop","is",is,"zt",zt,
c$$$     $              iotty, outmod )
c$$$            END IF

            call green

c$$$            IF ( (io == 17) .AND. (nsew == 1) ) THEN
c$$$               call atpoint("AFTER GREEN","is",is,"zt",zt,
c$$$     $              iotty, outmod )
c$$$            END IF
            
c     
c     Divide BVAL by twopi to account for definition of script G to be 
c     consistent with K.
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
               ipm =  + 1.0
c.... Change sign t0 - 1.0 for interior problem:
               IF ( linterior == 1 ) ipm = - 1.0
               chir(nsew,io) = chir(nsew,io) + ipm * bval
     $              * ( cslth(is,l1)*zbr - snlth(is,l1)*zbi )
               chii(nsew,io) = chii(nsew,io) + ipm * bval
     $              * ( snlth(is,l1)*zbr + cslth(is,l1)*zbi )
c     
 60            continue
            end do
         end do
c     
         chir(nsew,io) = 0.5 * isg*dtpw * chir(nsew,io)
         chii(nsew,io) = 0.5 * isg*dtpw * chii(nsew,io)
c     
      end do
c     
      return
      end
c     
c........................................................
      SUBROUTINE loops
c........................................................
c     
      USE vacuum13_mod

      include 'vacuum1.inc'
      include 'vacuum2.inc'
      include 'vacuum3.inc'
      include 'vacuum5.inc'
      include 'vacuum8.inc'

c$$$      REAL, DIMENSION(*) :: xloop, zloop ! This moved to VACUUM13_ALLOC
      CHARACTER*(130) titvacin

      DIMENSION xpla(nths), zpla(nths)

      REAL, DIMENSION(:,:), ALLOCATABLE :: xlooptem, zlooptem
      REAL, DIMENSION(:), ALLOCATABLE :: sourcemat

      EQUIVALENCE (xpla,xinf), (zpla,zinf)

      call boundsi ( xpla,zpla, 1,mth, xmnp,xmxp, zmnp,zmxp,
     $     ixn,ixx,izn,izx )

      plrad = 0.5 * ( xmxp-xmnp )
      xmaj = 0.5 * ( xmxp+xmnp )

c$$$      if ( lpsub .ne. 1 ) go to 150 ! 10-04-2008

      IF ( lpsub == 0 ) THEN
         nzdimlp = nloop + 3*nloopr + ntloop + 5
         CALL vacuum13_alloc ( nzdimlp )   ! xloop, zloop
c.. Use namelist inputs here:
         xloop(1:nloop) = xloopnl(1:nloop)
         zloop(1:nloop) = zloopnl(1:nloop)
         GO TO 150
      END IF

      GO TO ( 1,2,3,4,5,6 ), iloop

    1 continue

      nzdimlp = nloop + 3*nloopr + ntloop + 5

      CALL vacuum13_alloc ( nzdimlp )   ! xloop, zloop

      dtl = twopi / nloop

      do i = 1, nloop
         the = (i-1)*dtl
         xloop(i) = xmaj + xofsl*plrad
     $        + plrad*(1.0+aloop-xofsl) * cos(the+dloop*sin(the))
         zloop(i) = bloop*plrad*(1.0+aloop-xofsl) * sin(the)
      end do

      go to 150

    2 continue

c......PBX loops.

      nloop = 16

      nzdimlp = nloop + 3*nloopr + ntloop + 5
      CALL vacuum13_alloc ( nzdimlp )     ! xloop, zloop

      xloop(1)  = 2.17
      xloop(2)  = 2.17
      xloop(3)  = 2.17
      xloop(4)  = 1.45
      xloop(5)  = 1.25
      xloop(6)  = 0.96
      xloop(7)  = 0.82
      xloop(8)  = 1.05
      xloop(9)  = 0.82
      xloop(10) = 0.96
      xloop(11) = 1.25
      xloop(12) = 1.45
      xloop(13) = 1.05
      xloop(14) = 2.17
      xloop(15) = 1.05
      xloop(16) = 2.17
c     
      zloop(1)  = -0.26
      zloop(2)  =  0.00
      zloop(3)  =  0.26
      zloop(4)  =  0.67
      zloop(5)  =  0.67
      zloop(6)  =  0.67
      zloop(7)  =  0.45
      zloop(8)  = -0.05
      zloop(9)  = -0.45
      zloop(10) = -0.67
      zloop(11) = -0.67
      zloop(12) = -0.67
      zloop(13) = -0.05
      zloop(14) =  0.00
      zloop(15) = -0.05
      zloop(16) =  0.0
c     
      go to 150
c     
    3 continue
c
c.....  nloop loops at constant distance from plasma.
c       use zork1, zork2 as storage
c
c .. Equal arcs on plasma first:

      nzdimlp = nloop + 3*nloopr + ntloop + 5
      CALL vacuum13_alloc ( nzdimlp )   ! xloop, zloop

      call eqarcw (xinf,zinf, wxgrd,wzgrd, wkgrd,wkth,workl,mth1 )
      wxgrd(mth2) = wxgrd(2)
      wzgrd(mth2) = wzgrd(2)
      do i = 2, mth1
         alph = atan2m ( wxgrd(i+1)-wxgrd(i-1), wzgrd(i-1)-wzgrd(i+1) )
         zork1(i) = wxgrd(i) + aloop*plrad * cos(alph)
         zork2(i) = wzgrd(i) + aloop*plrad * sin(alph)
      end do
c
      zork1(1) = zork1(mth1)
      zork2(1) = zork2(mth1)
      zork1(mth2) = zork1(2)
      zork2(mth2) = zork2(2)
c
      call trans ( zork1,mth, xloop,nloop )
      call trans ( zork2,mth, zloop,nloop )
c
      go to 150
c     
    4 continue

      IF ( lfarw > 0 ) THEN
         WRITE ( IOTTY, '(5x, "iloop = 4 inconsistent with no wall.
     $        WILL ABORT!!")' )
         WRITE ( OUTMOD, '(5x, "iloop = 4 inconsistent with no wall.
     $        WILL ABORT!!")' )
         CALL EXIT (0)
      END IF
c
c.....  nloop loops at constant distance from wall
c       use zork1, zork2 as storage

      nzdimlp = nloop + 3*nloopr + ntloop + 5
      CALL vacuum13_alloc ( nzdimlp )      ! xloop, zloop

      xwal(mth2) = xwal(2)
      zwal(mth2) = zwal(2)
      do i = 2, mth1
         alph = atan2m ( xwal(i+1)-xwal(i-1), zwal(i-1)-zwal(i+1) )
         zork1(i) = xwal(i) - aloop*plrad * cos(alph)
         zork2(i) = zwal(i) - aloop*plrad * sin(alph)
      end do
c
      zork1(1) = zork1(mth1)
      zork2(1) = zork2(mth1)
      zork1(mth2) = zork1(2)
      zork2(mth2) = zork2(2)
c
      call trans ( zork1,mth, xloop,nloop )
      call trans ( zork2,mth, zloop,nloop )

      GO TO 150

    5 CONTINUE

cc..   Loop positions on a rectangular (X-Z) grid, read in input from
c      VACIN5

c$$$      nxlin = 21
c$$$      nzlin = 21
      dxlin = 1.0 / (nxlpin-1)
      dzlin = 1.0 / (nzlpin-1)
      nxzlin = nxlpin * nzlpin
      nloop = nxzlin       ! Note: nloop not from namelist now.

      nzdimlp = nloop +3* nloopr + ntloop + 5
      CALL vacuum13_alloc ( nzdimlp )   ! xloop.loop

      ALLOCATE ( xlooptem(nxlpin,nzlpin), zlooptem(nxlpin,nzlpin) )
      ALLOCATE ( sourcemat(nxzlin) )

c... Dummy Data:

      sourcemat = (/ (((i-1)*dxlin, i=1,nxlpin),j=1,nzlpin) /)
      sourcemat = sourcemat * (xlpmax-xlpmin) + xlpmin

      xlooptem(1:nxlpin,1:nzlpin) =
     $     RESHAPE ( sourcemat(1:nxzlin), SHAPE = (/ nxlpin,nzlpin /) )

      sourcemat = (/ (((j-1)*dzlin, i=1,nxlpin),j=1,nzlpin) /)
      sourcemat = sourcemat * (zlpmax-zlpmin) + zlpmin

      zlooptem(1:nxlpin,1:nzlpin) =
     $     RESHAPE ( sourcemat(1:nxzlin), SHAPE = (/ nxlpin,nzlpin /) )

      xloop(1:nxzlin) =
     $     (/ ( (xlooptem(i,j), j=1,nzlpin), i=1,nxlpin ) /)
      zloop(1:nxzlin) =
     $     (/ ( (zlooptem(i,j), j=1,nzlpin), i=1,nxlpin ) /)
      
c$$$      CALL matwrt ( xlooptem, nxlpin,nzlpin, nxlpin,nzlpin, "xlooptem" )
c$$$      CALL matwrt ( zlooptem, nxlpin,nzlpin, nxlpin,nzlpin, "zlooptem" )

      DEALLOCATE ( xlooptem, zlooptem, sourcemat )

      GO TO 150

 6    CONTINUE

c..   Loop positions read in from file VACIN6_2

      OPEN ( UNIT = 62, FILE = 'vacin6_2')
      WRITE ( IOTTY,  '(/,5x,"reading vacin6_2 in LOOPS",/)' )
      write ( OUTMOD, '(/,5x,"reading vacin6_2 in LOOPS",/)' )

      READ(62,*) nloop
      READ(62,'(a)') titvacin

      WRITE ( IOTTY, '(/"nloop = ", I5, A ,/)' ), nloop, titvacin
      WRITE ( OUTMOD, '(/"nloop = ", I5, A ,/)' ), nloop, titvacin

c-----------------------------------------------------------------------
c     read arrays. MUST ALLOCATE STORAGE FOR LOOPS SO CALL VACUUM13_ALLOC
c-----------------------------------------------------------------------

      nzdimlp = nloop + 3*nloopr + ntloop + 5
      CALL vacuum13_alloc ( nzdimlp )     ! xloop, zloop

      DO i = 1, nloop
         READ ( 62, * )  xloop(i), zloop(i)
      END DO

c$$$      WRITE ( IOTTY, '("LOOPS: ",/, (2ES12.4))' )
c$$$     $     (xloop(i),zloop(i), i = 1, nloop)

      GO TO 150

 150  CONTINUE

      IF ( nloopr == 0 ) go to 200

      IF ( a < 10.0 ) THEN
         drl = plrad * a / (nloopr)
      ELSE 
         drl = 2.0 * plrad / nloopr
      ENDIF

      DO i = 1, nloopr
         xloop(nloop+i) = xmaj + plrad + i * drl - drl/2.0
         zloop(nloop+i) = zpla(ixx)
      END DO
c     
      DO i = 1, nloopr
         il = nloop + nloopr + i
         xloop(il) = xmaj - plrad - i * drl + drl/2.0
         IF ( xloop(il) <= .01*xmaj ) xloop(il) = .01*xmaj
         zloop(il) = zpla(ixn)
      END DO
c     
      DO i = 1, nloopr
         il = nloop + 2*nloopr + i
         xloop(il) = xpla(izx)
         zloop(il) = zmxp + i * drl - drl/2.0
      END DO

 200  CONTINUE

      CLOSE ( UNIT = 62 )

      RETURN
      END
C     
C......................................
      SUBROUTINE chicon ( xil )
c.....................................

      USE vacuum12_mod

      include 'vacuum1.inc'
      include 'vacuum3.inc'
      include 'vacuum2.inc'
c     
      dimension xwal(nths), zwal(nths), xpla(nths), zpla(nths),
     .     chia(101,101), chiar(101,101), chiai(101,101),
     $     cwrkr(101,101), cwrki(101,101)
      dimension grdgre(nths2,nths2),grwp(nths,nths),
     .     grri(nths2,nfm2)
      dimension ch(101)
      dimension xil(*)
c     
      common / bigv1 / grdgre, grwp, grri
c
c$$$      common / gridc / xobs(ndimlp),zobs(ndimlp)
c     
      equivalence (xpla,xinf), (zpla,zinf)
c     
c$$$      parameter(neqv1=(nths-1)/2+1, neqv2=nths+1, neqv3=3*(nths+1)/2  )
c$$$c     parameter(neqv1=nthsq+1, neqv2=neqv1+nthsq, neqv3=neqv2+nthsq,
c$$$c     .          neqv4=neqv3+nths*nfm)
c$$$c     equivalence ( grdgre(1,1), chiai(1,1) )
c$$$c     equivalence ( grdgre(neqv1), chiar(1,1) )
c$$$c     equivalence ( grdgre(neqv2), chia(1,1) )
c$$$c     equivalence ( grdgre(neqv3), chiwc(1,1) )
c$$$c     equivalence ( grdgre(neqv4), chiws(1,1) )
c$$$c     
c$$$      equivalence ( grdgre(1,1), chiai(1,1) )
c$$$      equivalence ( grdgre(neqv2,neqv1), chiar(1,1) )
c$$$      equivalence ( grdgre(1,neqv2), chia(1,1) )
c     
      equivalence ( grdgre(1,1), chiai(1,1) )
      equivalence ( grdgre(nptg2x,nptg12), chiar(1,1) )
      equivalence ( grdgre(1,nptg13), chia(1,1) )

      if ( nph .eq. 0 ) return
c     c
c.. New variable, lfarw for SPARK2 to reset farwal.
c
c$$$      lfarw = 0
c$$$      if ( a .gt. 10.0 ) lfarw = 1
c$$$      write ( iotty, '("a = ", 1pe12.3, " lfarw in CHICON = ", i4 )' )
c$$$     $     a, lfarw
c$$$      write ( outmod, '("a = ", 1pe12.3, " lfarw in CHICON = ", i4 )' )
c$$$     $     a, lfarw
      jmax1 = lmax(1) - lmin(1) + 1
      mth12 = 2*mth
      if ( lfarw .gt. 0 ) mth12 = mth
c     
c.....set up window.
c     
      call bounds(xpla,zpla,1,mth,xmnp,xmxp,zmnp,zmxp)
      xmin = xmnp
      xmax = xmxp
      zmin = zmnp
      zmax = zmxp
c     
      plrad = 0.5 * ( xmxp - xmnp )
      xmaj = 0.5 * ( xmxp + xmnp )
c     
      if ( lfarw .gt. 0 ) go to 20
c     
      call wwall ( mth1, xwal,zwal )
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
      dz = zmax - zmin
      dxz = amax1(dx,dz)
      xmaxw = 0.5 * ( xmax + xmin ) + 0.75*dxz
      xminw = 0.5 * ( xmax + xmin ) - 0.75*dxz
      zmaxw = 0.75 * dxz
      zminw = - zmaxw
c     
      xming = xminw - 0.05*(xmaxw-xminw)
      xmaxg = xmaxw + 0.05*(xmaxw-xminw)
      zming = zminw - 0.05*(zmaxw-zminw)
      zmaxg = zmaxw + 0.05*(zmaxw-zminw)
c     
c.....get array of observer points, and initialize chi matrix.
c     
      CALL vacuum12_alloc ( mx, mz )

      do i = 1, mx
         xobs(i) = xminw + (i-1) * (xmaxw-xminw) / (mx-1)
         do j = 1, mz
            chiar(i,j) = 0.0
            chiai(i,j) = 0.0
            cwrkr(i,j) = 0.0
            cwrki(i,j) = 0.0
         end do
      end do 
c     
      do 210 j = 1, mz
         zobs(j) = zminw + (j-1) * (zmaxw-zminw) / (mz-1)
 210  continue
c     
c.....pass chi on source surface.
c     
c.....plasma contribution.
c     
      dtpw = dth
      ns = mth
c     
      do l1 = 1, jmax1
         do i = 1, mth
c     il11 = (l1-1)*mth12 + i
c     il12 = jmax1*mth12 + il11
            chiwc(i,l1) = grri(i,l1)
            chiws(i,l1) = grri(i,jmax1+l1)
         end do
      end do
c     
      isg = 1
      call chimat ( xpla,zpla,isg, chiwc,chiws, ns,1, cwrkr,cwrki, xil
     $     )
      do i = 1, mx
         do j = 1, mz
            chiar(i,j) = cwrkr(i,j)
            chiai(i,j) = cwrki(i,j)
            cwrkr(i,j) = 0.0
            cwrki(i,j) = 0.0
         end do
      end do
c     
c.....conductor contribution.
c     
      do l1 = 1, jmax1
         do  i = 1, mw
c     il11 = (l1-1)*mth12 + mth + i
c     il12 = jmax1 * mth12 + il11
            chiwc(i,l1) = grri(mth+i,l1)
            chiws(i,l1) = grri(mth+i,jmax1+l1)
         end do
      end do

c     
      isg = -1
      call chimat ( xwal,zwal,isg,chiwc,chiws, ns,0, cwrkr,cwrki
     $     ,xil)
c     
      do i = 1, mx
         do j = 1, mz
            chiar(i,j) = chiar(i,j) + cwrkr(i,j)
            chiai(i,j) = chiai(i,j) + cwrki(i,j)
         end do
      end do 
c     
c......plot contours for different phi's
c     
      dphi = 0.5 * pye / nph
c     
      do 300 iph = 1, nph
c     
         phi = (iph-1) * dphi
         sinnph = sin ( n*phi )
         cosnph = cos( n*phi )
         cmin = 1.0e10
         cmax = -1.0e10
c     
         do i = 1, mx
            do j = 1, mz
               chia(i,j) = chiar(i,j)*sinnph - chiai(i,j)*cosnph
               if ( cmax .lt. chia(i,j) ) cmax = chia(i,j)
               if ( cmin .gt. chia(i,j) ) cmin = chia(i,j)
            end do
         end do
c     
         call maps(xming,xmaxg,zming,zmaxg,0.1,0.815,0.285,1.0)
c     
         k1 = -21
         k2 = 0
         cfac = 1.0
         ch(1) = cfac*cmin
         ch(2) = cfac*cmax
c     
         write ( iotty, 1010 ) cmin, cmax
         write ( outmod, 1010 ) cmin, cmax
 1010    format ( 1x, ' cmin = ',1pe12.4, ' cmax = ',1pe12.4,/
     $        )
c     
c.... plot conductors and plasma shape..
c     
         call points ( xpla,zpla, mth, -1,-1,0.,0. )
         call points ( xwal,zwal, mw, -1,-1,0.,0. )
c     
         call rcontr(k1,ch,k2, chia,nths, xobs,1,mx,1, zobs,1
     $        ,mz,1)
c     
         call framep(jobid,ff)  
c     
 300  continue

      CALL vacuum12_dealloc
c     
      RETURN
      END
C     
c......................................................
c     
      SUBROUTINE chimat(xsce,zsce,isg,creal,cimag,ns,ip,chiar,chiai,xil)
c......................................................

c.. Only Needed by CHICON

      USE vacuum12_mod

      include 'vacuum1.inc'
      include 'vacuum3.inc'
      include 'vacuum2.inc'
c     
      dimension the(nths),xpp(nths),zpp(nths),xscp(nths),zscp(nths)
      dimension iop(2), ww1(nths),ww2(nths),ww3(nths),tab(3)
      dimension xil(*), xsce(*),zsce(*)
      dimension creal(nths,nfm), cimag(nths,nfm)
      dimension chiar(101,101), chiai(101,101)
c     
      real nq
c     
c$$$      common / gridc / xobs(ndimlp),zobs(ndimlp)

c     
c     Factpi for dividing BVAL to be consistent with AVAL.
c     
      factpi = twopi
      jmax1 = lmax(1) - lmin(1) + 1
      q = qa1
      nq = n * q
      dtpw = twopi / ns
c     
      ns1 = ns + 1
c     
      do 10 i = 1, ns1
         the(i) = (i-1) * dth
 10   continue
c     
c......get derivatives. replace splines later.
c     
      iop(1) = 4
      iop(2) = 4
      call spl1d1(ns1,the,xsce,xpp,iop,1,ww1,ww2,ww3)
      call spl1d1(ns1,the,zsce,zpp,iop,1,ww1,ww2,ww3)
c
      do 20 i = 1, ns1
         theta = (i-1) * dth
         call spl1d2 ( ns1,the, xsce,xpp, 1, theta, tab )
         xscp(i) = tab(2)
         call spl1d2 ( ns1,the, zsce,zpp, 1, theta, tab )
         zscp(i) = tab(2)
 20   continue
c     
      do 210 i = 1, mx
c     
         xs = xobs(i)
c     
         do 200 j = 1, mz
c     
            zs = zobs(j)
c     
            do 100 is = 1, ns
c     
               xt = xsce(is)
               zt = zsce(is)
               xtp = xscp(is)
               ztp = zscp(is)
               theta = (is-1) * dtpw
c     
               call green
c     
c     Divide BVAL by twopi to account for definition of script G to be
c     consistent
c     with K.
c     
               bval = bval / factpi
c     
               do 80 l1 = 1, jmax1
c     
                  ll = lmin(1) - 1 + l1
                  alnq = ll - nq
                  xilnq = xil(l1) * alnq
                  elthnq = ll * theta + nq*delta(is)
                  sinlth = sin(elthnq)
                  coslth = cos(elthnq)
c     
                  chiar(i,j) = chiar(i,j) + creal(is,l1) * aval*xilnq
                  chiai(i,j) = chiai(i,j) + cimag(is,l1) * aval*xilnq
c     
                  if ( ip .eq. 0 ) go to 60
c     
                  chiar(i,j) = chiar(i,j) + bval * coslth * xilnq
                  chiai(i,j) = chiai(i,j) + bval * sinlth * xilnq
c     
 60               continue
c     
 80            continue
 100        continue
c     
            chiar(i,j) = 0.5 * isg * dtpw * chiar(i,j)
            chiai(i,j) = 0.5 * isg * dtpw * chiai(i,j)
c     
 200     continue
 210  continue
c     
      call matwrt(chiar,101,101,mx,mz,"chiar               " )
      call matwrt(chiai,101,101,mx,mz,"chiai               " )
c     
      return
      end
c
c.............................................................
      function atan2m ( z,x )
c.............................................................
c
      include 'vacuum1.inc'
c
      if ( abs(x) .le. epszer ) then
         atan2m = sign ( pye/2.0, z )
         return
      end if
c
      atan2m = atan2 (z,x)
c
      return
      end
c.............................................................
      
