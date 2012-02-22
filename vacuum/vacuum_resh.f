c.................................................................
      SUBROUTINE resshel ( v_cpwwp, nd1 )
c.................................................................

      include 'vacuum1.inc'
      include 'vacuum2.inc'
      include 'vacuum3.inc'
      include 'vacuum5.inc'
      include 'vacuum6.inc'
      include 'vacuum7.inc'
      include 'vacuum8.inc'
c     
      REAL, DIMENSION(:,:), ALLOCATABLE :: wprr, wpii, wpri, wpir,
     $     zrr, zii, zri, zir, wbig, wbig2, wrkb3, wrkb4, wrkb5
c     
      REAL, DIMENSION(nd1,nd1) :: v_cpwwp

c...  Allocate the matrix, oplp, for storing the open loop stuff.
c     REAL, DIMENSION(:,:), ALLOCATABLE :: oplp
c     
      dimension xpla(nths), zpla(nths)
      dimension wrr(nfm,nfm), wii(nfm,nfm), wri(nfm,nfm), wir(nfm,nfm)
      dimension vacmti(nfm,nfm)
      dimension wkxx(nths), wkyy(nths)
c     
      integer tmth
c     
      common / bigv1 / grdgre(nths2,nths2), grwp(nths,nths),
     $     grri(nths2,nfm2)
c     
      parameter ( neqv1=1, neqv2=(nfm-1)/2+1,neqv3=nfm+1,
     .     neqv4=3*(nfm+1)/2,neqv5=2*nfm+1 )
c$$$  c
c$$$  equivalence ( grdgre(1,1), vacmti(1,1) )
c$$$  equivalence ( grdgre(1,neqv2), wrr(1,1) )
c$$$  equivalence ( grdgre(1,neqv3), wii(1,1) )
c$$$  equivalence ( grdgre(1,neqv4), wir(1,1) )
c$$$  equivalence ( grdgre(1,neqv5), wri(1,1) )

      equivalence ( grdgre(1,1), wrr(1,1) )
      equivalence ( grdgre(npfg2x,npfg12), wii(1,1) )
      equivalence ( grdgre(1,npfg13), wir(1,1) )
      equivalence ( grdgre(npfg2x,npfg14), wri(1,1) )
c     
      dimension vacpstr(nfmsq),vacpsti(nfmsq)
c     
      call atpoint ( "IN RESHEL","jwal1", jwal1,"xinf(1)", xinf(1),
     $     iotty, outmod )
c     
      ln = lmin(1)
      lx = lmax(1)
      jmax1 = lx - ln + 1
      jmax2 = 2*jmax1
      j1v = nfm
      j2v = nfm
c     
      ALLOCATE ( wprr(jmax1,jmax1), wpii(jmax1,jmax1),
     $     wpri(jmax1,jmax1), wpir(jmax1,jmax1),
     $     zrr(jmax1,jmax1), zii(jmax1,jmax1), zri(jmax1,jmax1),
     $     zir(jmax1,jmax1), wbig(jmax2,jmax2), wbig2(jmax2,jmax2),
     $     wrkb3(jmax2,jmax2), wrkb4(jmax2,jmax2), wrkb5(jmax2,jmax2) )

      wrkb3 = 0.0

c     ALLOCATE ( oplp(jmax2+jwal2,jmax2+jwal2) )
c     
c..   tpisqi is optional factor to take care norm and to
c     negate the TWOPI factor was  in SUBROUTINE FORANV.
c     tpisqi = 1.0 / (twopi*pye)
c...  Set tpisqi = 1.0 so that wall observer is consistent with other matrices.
      tpisqi = 1.0
c     
c     Get delta Cww:
c     
      write ( outmod,
     $     '(/,"Getting  delta C_(lw,lw) in shell functions.")' )
      write ( iotty,
     $     '(/,"Getting  delta C_(lw,lw) in shell functions.")' )
c$$$  c
c$$$  c..Uncertainty in sign so experiment:
c$$$  c
c$$$  csign = +1.0
c$$$  write ( iotty,  '(/,20x,"CSIGN for W-OUT + CSIGN * W-IN = ",
c$$$  $     1pe12.3)' ) csign
c$$$  write ( outmod, '(/,20x,"CSIGN for W-OUT + CSIGN * W-IN = ",
c$$$  $     1pe12.3)' ) csign
c$$$  do i = 1, mth
c$$$  do j = 1, jwal2
c$$$  gwot(i,j) = gwot(i,j) + csign * gwin(mth+i,j)
c$$$  end do
c$$$  end do
c     
c     . Expand Delta C_+ww  obs. pts in shell functions.
c     Use arr, aii, ari, ari, etc
c     as storage. Twopi factor was included in FORANV.
c     Remove Twopi factor and divide by another pi to
c     account for the orthonormality.
c     
c..   Remember GWOT(:,jwal2+1) contains the coil source:
c     ----------------------------------------------------
c     
c     Append coil source in lw space to WII(:,jwal1+1), WIR(:,jwal1+1):
c     
      zzfac = dth * twopi * tpisqi
      wrr(1:jwal1,1:jwal1) = zzfac *
     $     matmul ( transpose(rwvce(1:mth,1:jwal1)),
     $     gwot(1:mth,1:jwal1) )
      wii(1:jwal1,1:jwal1+1) = zzfac *
     $     matmul ( transpose(rwvco(1:mth,1:jwal1)),
     $     gwot(1:mth,jwal1+1:jwal2+1) )
      wri(1:jwal1,1:jwal1) = zzfac *
     $     matmul ( transpose(rwvco(1:mth,1:jwal1)),
     $     gwot(1:mth,1:jwal1) )
      wir(1:jwal1,1:jwal1+1) = zzfac *
     $     matmul ( transpose(rwvce(1:mth,1:jwal1)),
     $     gwot(1:mth,jwal1+1:jwal2+1) )
c     
      call atpoint ( "At:C_+(lw,lw)-rr","jwal1",jwal1,
     $     "wrr(1,1)", wrr(1,1), iotty, outmod )
      call matwrtn ( wrr, j1v,j2v, 1,1,jwal1,jwal1,jdel,jdel,
     $     "C_+(lw,lw)-rr", outmod,0 )
      call atpoint ( "At C_+(lw,lw)-ii","jmax1",jmax1,
     $     "wii(1,1)", wii(1,1), iotty, outmod )
      call matwrtn ( wii, j1v,j2v, 1,1,jwal1,jwal1+1,jdel,jdel,
     $     "C_+(lw,lw)-ii + COIL", outmod,0 )
      call matwrtn ( wri, j1v,j2v, 1,1,jwal1,jwal1,jdel,jdel,
     $     "C_+(lw,lw)-ri", outmod,0 )
      call matwrtn ( wir, j1v,j2v, 1,1,jwal1,jwal1+1,jdel,jdel,
     $     "C_+(lw,lw)-ir + COIL", outmod,0 )
c     
      call vacasym ( wrr, nfm,jwal1,"wrr", outmod,iotty )
      call vacasym ( wii, nfm,jwal1,"wii", outmod,iotty )
      call vacasym ( wri, nfm,jwal1,"wri", outmod,iotty )
      call vacasym ( wir, nfm,jwal1,"wir", outmod,iotty )
c     
c..   Now C_-ww:
c..   Store gwin in gwot:
c     
      gwot(1:mth,1:jwal2) = gwin(mth+1:2*mth,1:jwal2)
c     
      zzfac = dth * twopi * tpisqi
      zrr(1:jwal1,1:jwal1) = zzfac *
     $     matmul ( transpose(rwvce(1:mth,1:jwal1)),
     $     gwot(1:mth,1:jwal1) )
      zii(1:jwal1,1:jwal1) = zzfac *
     $     matmul ( transpose(rwvco(1:mth,1:jwal1)),
     $     gwot(1:mth,jwal1+1:jwal2) )
      zri(1:jwal1,1:jwal1) = zzfac *
     $     matmul ( transpose(rwvco(1:mth,1:jwal1)),
     $     gwot(1:mth,1:jwal1) )
      zir(1:jwal1,1:jwal1) = zzfac *
     $     matmul ( transpose(rwvce(1:mth,1:jwal1)),
     $     gwot(1:mth,jwal1+1:jwal2) )
c     
      call matwrtn ( zrr, jmax1,jmax1, 1,1,jwal1,jwal1,jdel,jdel,
     $     "C_-(lw,lw)-rr", outmod,0 )
      call matwrtn ( zii, jmax1,jmax1, 1,1,jwal1,jwal1,jdel,jdel,
     $     "C_-(lw,lw)-ii", outmod,0 )
      call matwrtn ( zri, jmax1,jmax1, 1,1,jwal1,jwal1,jdel,jdel,
     $     "C_-(lw,lw)-ri", outmod,0 )
      call matwrtn ( zir, jmax1,jmax1, 1,1,jwal1,jwal1,jdel,jdel,
     $     "C_-(lw,lw)-ir", outmod,0 )
c     
      call vacasym ( zrr, jmax1,jwal1,"zrr", outmod,iotty )
      call vacasym ( zii, jmax1,jwal1,"zii", outmod,iotty )
      call vacasym ( zri, jmax1,jwal1,"zri", outmod,iotty )
      call vacasym ( zir, jmax1,jwal1,"zir", outmod,iotty )
c     
c     
c..   Uncertainty in sign in differencing C, so experiment:
c     
c..   Extract the feedback coil contrib. before proceeding. RWVCE first.
c     
      coil00(1:jwal2) = (/ wir(1:jwal1,jwal1+1),wii(1:jwal1,jwal1+1) /)
      call vecwrt ( jwal2, coil00, "coil00 in lw basis", 1,jwal2,
     $     outmod, iotty )
c     
      csign = +1.0
      write ( iotty,  '(/,20x,"CSIGN for W-OUT + CSIGN * W-IN = ",
     $     1pe12.3)' ) csign
      write ( outmod, '(/,20x,"CSIGN for W-OUT + CSIGN * W-IN = ",
     $     1pe12.3)' ) csign

      wrr(1:jwal1,1:jwal1) = wrr(1:jwal1,1:jwal1)
     $     + csign * zrr(1:jwal1,1:jwal1)
      wii(1:jwal1,1:jwal1) = wii(1:jwal1,1:jwal1)
     $     + csign * zii(1:jwal1,1:jwal1)
      wri(1:jwal1,1:jwal1) = wri(1:jwal1,1:jwal1)
     $     + csign * zri(1:jwal1,1:jwal1)
      wir(1:jwal1,1:jwal1) = wir(1:jwal1,1:jwal1)
     $     + csign * zir(1:jwal1,1:jwal1)
c     
c...  Now at end of Equation IV-141-6,7
c     
c...  Store even-odd in  4x4 matrix, wbig. Normalise out the eigenvalues:
c     
      write ( outmod,
     $     '(/,"Normalize D-C_(lw,lw) with sqrt(eigenvalues):")' )
      write ( iotty,
     $     '(/,"Normalize D-C_(lw,lw) with sqrt(eigenvalues):")' )
c     
      wbig(1:jwal1,1:jwal1) = reshape
     $     ((/ ((sqrt(rwvle(i)),i=1,jwal1),j=1,jwal1)/),(/jwal1,jwal1/))
     $     * wrr(1:jwal1,1:jwal1) * reshape
     $     ((/ ((sqrt(rwvle(j)),i=1,jwal1),j=1,jwal1)/),(/jwal1,jwal1/))

      wbig(jwal1+1:2*jwal1,jwal1+1:2*jwal1) = reshape
     $     ((/ ((sqrt(rwvlo(i)),i=1,jwal1),j=1,jwal1)/),(/jwal1,jwal1/))
     $     * wii(1:jwal1,1:jwal1) * reshape
     $     ((/ ((sqrt(rwvlo(j)),i=1,jwal1),j=1,jwal1)/),(/jwal1,jwal1/))

      wbig(1:jwal1,jwal1+1:2*jwal1) = reshape
     $     ((/ ((sqrt(rwvle(i)),i=1,jwal1),j=1,jwal1)/),(/jwal1,jwal1/))
     $     * wri(1:jwal1,1:jwal1) * reshape
     $     ((/ ((sqrt(rwvlo(j)),i=1,jwal1),j=1,jwal1)/),(/jwal1,jwal1/))

      wbig(jwal1+1:2*jwal1,1:jwal1) = reshape
     $     ((/ ((sqrt(rwvlo(i)),i=1,jwal1),j=1,jwal1)/),(/jwal1,jwal1/))
     $     * wir(1:jwal1,1:jwal1) * reshape
     $     ((/ ((sqrt(rwvle(j)),i=1,jwal1),j=1,jwal1)/),(/jwal1,jwal1/))
c     
      call matwrtn ( wbig,jmax2,jmax2, 1,1,jwal2,jwal2,jdel,jdel,
     $     "Norm Big D-C_(lw,lw)", outmod,0 )
c     
      call vacasym ( wbig, jmax2,jwal2,"Norm D-C_(lw,lw)",
     $     outmod,iotty )

c..   Store this in oppl(l_w,l_w) for the open loop matrix:
      oplp(jmax2+1:jmax2+jwal2,jmax2+1:jmax2+jwal2) =
     $     wbig(1:jwal2,1:jwal2)
c     
c..   Now the coil:
      coil00(1:jwal2) = coil00(1:jwal2) *
     $     (/ sqrt(rwvle(1:jwal1)),sqrt(rwvlo(1:jwal1)) /)
      call vecwrt ( jwal2, coil00, "coil00_lw * sqrt(lmda)", 1,jwal2,
     $     outmod, iotty )
c     
c...  IV-148(1) to IV-148(2):  C(p,w)  oooooooooooooooooooooooooooooooo
c     
c...  Transform C(p,lw) to C(lp,lw). Store in wrkb5.
c     
c...  zzfac to correspond to that in FORANV.
c     
c     Storage:  |   Cos-Even      Cos-Odd |
c               | - Sin-Even    - Sin-Odd |
c     
      zzfac = dth * twopi
      wrkb5(1:jmax1,1:jwal1) = zzfac *
     $     matmul ( transpose(cslth(1:mth,1:jmax1)),
     $     gwin(1:mth,1:jwal1) )
      wrkb5(jmax1+1:jmax2,1:jwal1) = - zzfac *
     $     matmul ( transpose(snlth(1:mth,1:jmax1)),
     $     gwin(1:mth,1:jwal1) )
      wrkb5(1:jmax1,jwal1+1:jwal2) = zzfac *
     $     matmul ( transpose(cslth(1:mth,1:jmax1)),
     $     gwin(1:mth,jwal1+1:jwal2) )
      wrkb5(jmax1+1:jmax2,jwal1+1:jwal2) = - zzfac *
     $     matmul ( transpose(snlth(1:mth,1:jmax1)),
     $     gwin(1:mth,jwal1+1:jwal2) )
c     
      call matwrtn ( wrkb5,jmax2,jmax2, 1,1,jmax2,jwal2,jdel,jdel,
     $     "C(lp,lw)", outmod,0 )
c     
c     . IV-149(7): Right multiply by sqrt(eigenvalues), to get \bar C(lp,lw)
c     
      do i = 1, jmax1
         wrkb5(i,1:jwal1) =
     $        wrkb5(i,1:jwal1) * sqrt(rwvle(1:jwal1))
         wrkb5(i,jwal1+1:jwal2) =
     $        wrkb5(i,jwal1+1:jwal2) * sqrt(rwvlo(1:jwal1))
         wrkb5(jmax1+i,1:jwal1) =
     $        wrkb5(jmax1+i,1:jwal1) * sqrt(rwvle(1:jwal1))
         wrkb5(jmax1+i,jwal1+1:jwal2) =
     $        wrkb5(jmax1+i,jwal1+1:jwal2) * sqrt(rwvlo(1:jwal1))
      end do
c     
      call matwrtn ( wrkb5,jmax2,jmax2, 1,1,jmax2,jwal2,jdel,jdel,
     $     "BarC(lp,lw)", outmod,0 )
      call matwrts ( wrkb5,jmax2,jmax2, 1,1, 9,9, 9,9,
     $     1,1,"BarC(lp,lw)", outmod,0 )
      call matwrts ( wrkb5,jmax2,jmax2, 10,10, 18,18, 9,9,
     $     10,10,"BarC(lp,lw)", outmod,0 )
      call matwrts ( wrkb5,jmax2,jmax2, 1,10, 9,18, 9,9,
     $     1,10,"BarC(lp,lw)", outmod,0 )
      call matwrts ( wrkb5,jmax2,jmax2, 10,1, 18,9, 9,9,
     $     10,1,"BarC(lp,lw)", outmod,0 )

      call matwrts ( wrkb5,jmax2,jmax2, 19,1, 27,9, 9,9,
     $     19,1,"BarC(lp,lw)", outmod,0 )
      call matwrts ( wrkb5,jmax2,jmax2, 28,10, 36,18, 9,9,
     $     28,10,"BarC(lp,lw)", outmod,0 )
      call matwrts ( wrkb5,jmax2,jmax2, 19,10, 27,18, 9,9,
     $     19,10,"BarC(lp,lw)", outmod,0 )
      call matwrts ( wrkb5,jmax2,jmax2, 28,1, 36,9, 9,9,
     $     28,1,"BarC(lp,lw)", outmod,0 )


c     . Store this in oplp(l_p,l_w) for the open loop matrix:
      oplp(1:jmax2,jmax2+1:jmax2+jwal2) = wrkb5(1:jmax2,1:jwal2)
c$$$      oplp(1:jmax1,jmax1+1:jmax1+jwal1) = wrkb5(1:jmax1,1:jwal1)
c$$$      oplp(1:jmax1,jmax1+jwal1+1:jmax1+jwal2) =
c$$$     $     wrkb5(jmax1+1:jmax2,jwal1+1:jwal2)

c     
c...  Will transform this later to the diagonal space (lp,ls).
c     
c...  ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

c     . Get Eigenmodes of Norm-Delta-C_(lw,lw) -> to ls basis.
c     
      write ( outmod, '(/,20x, "Eigen. of Norm-C(l_w,l_w.",/
     $     "To ls basis:")' )
      write ( iotty,  '(/,20x, "Eigen. of Norm-C(l_w,l_w.",/
     $     "To ls basis:")' )
c     
c..   Even Modes first
      wrkb4(1:jwal1,1:jwal1) = wbig(1:jwal1,1:jwal1)

c..   Transformation matrix (a.k.a column eigenvectors) in WRKB3.
c     Eigenvalues in WORK0:

      call mateig2a ( wrkb4, jmax2,jwal1, wrkb3,work0, 1,jwal1,
     $     ff, 2, 2, jobid, "Norm_Del-C-e", outmod )
      call matwrtn ( wrkb3, jmax2,jmax2, 1,1,jwal1,jwal1,jdel,jdel,
     $     "T_(lw,lw) - Even modes.", outmod,0 )
c     
c...  Calculate Time Constants for EVEN Eigenmodes of the Shell:
      write ( iotty, '(/,10x,
     $     "Calculate Time Consts for EVEN Eigenmodes of the Shell:")' )
      write ( outmod, '(/,10x,
     $     "Calculate Time Consts for EVEN Eigenmodes of the Shell:")' )
c     
      tauhat = amu0 * th_vess / eta_vess
      write ( iotty, '(/, 5x,
     $     "Tau_hat = Mu*DeltaZ/Eta = ", 1pe13.4 )' ) tauhat
      write ( outmod, '(/, 5x,
     $     "Tau_hat = Mu*DeltaZ/Eta = ", 1pe13.4 )' ) tauhat
c     
      write ( iotty, '(5x,"eta = ", 1pe13.4," Ohm-m", 3x,
     $     "Thickness = ",1pe13.4," m",//,
     $     5x, "Eig. No", 5x, "Eigenvalue", 3x, "Time, msec." )' )
     $     eta_vess, th_vess
      write ( outmod,'(5x,"eta = ", 1pe13.4," Ohm-m", 3x,
     $     "Thickness = ",1pe13.4," m",//,
     $     5x, "Eig. No", 5x, "Eigenvalue", 3x, "Time, msec." )' )
     $     eta_vess, th_vess
c     
      do i = 1, jwal1
         if ( abs(work0(i)) .le. epszer )
     $        work0(i) = sign(epszer,work0(i))
         tc_vess(i) = 1000.0 * tauhat / work0(i)
      end do
      do i = 1, jwal1
         write ( iotty,  '(5x, i4, 5x, 1p2e13.4)' )
     $        i, work0(i), tc_vess(i)*twopi**2
         write ( outmod, '(5x, i4, 5x, 1p2e13.4)' )
     $        i, work0(i), tc_vess(i)*twopi**2
      end do
c     
c..   Check diagonalization:
c     
      call atpoint ( "WRTN1", "jdel",jdel,"epszer",
     $     epszer, iotty, outmod )
c     call mtransr ( wrkb3,nwm2,nwm2, nwm2,nwm2, wrkb4 )
      wrkb4 = transpose (wrkb3)
c$$$  call atpoint ( "WRTN2", "jdel",jdel,"epszer",
c$$$  $     epszer, iotty, outmod )
c     call dvmat ( work0,1.0, wrkb4,nwm2,nwm2, jwal1,jwal1, 1 )
      wrkb4(1:jwal1,1:jwal1) =
     $     reshape((/ ((work0(i),i=1,jwal1),j=1,jwal1) /),
     $     (/ jwal1,jwal1 /)) * wrkb4(1:jwal1,1:jwal1)
c$$$  call atpoint ( "WRTN3", "jdel",jdel,"epszer",
c$$$  $     epszer, iotty, outmod )
c     call matmul0 ( wrkb3,wrkb4,nwm2,nwm2, jwal1, wbig2,nwm2 )
      wbig2(1:jwal1,1:jwal1) =
     $     matmul(wrkb3(1:jwal1,1:jwal1),wrkb4(1:jwal1,1:jwal1))
c$$$  call atpoint ( "WRTN4", "jdel",jdel,"epszer",
c$$$  $     epszer, iotty, outmod )
      call matwrtn ( wbig2, jmax2,jmax2, 1,1,jwal1,jwal1,jdel,jdel,
     $     "D-C(lw,lw)???", outmod,0 )
c     
c..   Calculate the Current Eigenmodes. Store in wbig2:
c     
      write ( iotty,  '(/,5x,"The EVEN Current Eigenmodes.")' )
      write ( outmod, '(/,5x,"The EVEN Current Eigenmodes.")' )

      wbig2(1:jwal1,1:jwal1) =
     $     transpose ( wrkb3(1:jwal1,1:jwal1) / reshape
     $     ((/ ((sqrt(rwvle(i)),i=1,jwal1),j=1,jwal1) /),
     $     (/ jwal1,jwal1 /)) )

      call matwrtn ( wbig2, jmax2,jmax2, 1,1,jwal1,jwal1,jdel,jdel,
     $     "EVEN Current, J_(ls,lw)", outmod,0 )
c     
c..   Plot the First Few Eigencurrents
c     Assume Imaginary part of current coeff. is zero:
c     .  Use wkxx, wkyy as storage.
c     
c     First, shift variables so that theta=0 is on the inside.
c     Use workl, workl2 as storage.
c     
      workl(1:mth1)  = (/ xwal(mth/2+1:mth),xwal(1:mth/2+1) /)
      workl2(1:mth1) = (/ zwal(mth/2+1:mth),zwal(1:mth/2+1) /)
c     
      nbngr = 129
      do js = jwal1-nwcur+1, jwal1
         write ( outmod, '("workl,workl2,ell,wgrd,wkxx,wkyy=")')
         wkxx = 0.0
         wkyy = 0.0
         wkxx(1:mth1) = matmul ( wbig2(js,1:jwal1),
     $        transpose(rwvce(1:mth1,1:jwal1)) )
         wkxx(1:mth1) = (/ wkxx(mth/2+1:mth),wkxx(1:mth/2+1) /)
         wkyy(1:mth1) = (/ wkyy(mth/2+1:mth),wkyy(1:mth/2+1) /)

         iwrto = 0
c     iwrto = outmod
         call kcur0 ( workl,workl2,ell,wgrd,mth1,wkxx,wkyy,n,
     $        nbngr,wkgrd,wxgrd,wzgrd,wkdxt,wkdzt,
     $        wkt2,wkp2,nths,0,"EVEN Eig-Cur",js, iwrto )
      end do
c..................End of Evens  ...................
c     Now Odd Modes  ...............
c     
      wrkb4(1:jwal1,1:jwal1) = wbig(jwal1+1:jwal2,jwal1+1:jwal2)

      call mateig2a ( wrkb4, jmax2,jwal1, wrkb3,work0,
     $     1,jwal1, ff, 2, 2, jobid, "Norm-Del-C-o", outmod )

      call matwrtn ( wrkb3, jmax2,jmax2, 1,1,jwal1,jwal1,jdel,jdel,
     $     "T_(lw,lw) - Odd modes.", outmod,0 )
c     
c...  Calculate Time Constants for ODD Eigenmodes of the Shell:
      write ( iotty, '(/,10x,
     $     "Time Consts for ODD Eigenmodes of the Shell:")' )
      write ( outmod, '(/,10x,
     $     "Time Consts. For ODD Eigenmodes of the Shell:")' )
c     
      tauhat = amu0 * th_vess / eta_vess
      write ( iotty, '(/, 5x,
     $     "Tau_hat = Mu*DeltaZ/Eta = ", 1pe13.4 )' ) tauhat
      write ( outmod, '(/, 5x,
     $     "Tau_hat = Mu*DeltaZ/Eta = ", 1pe13.4 )' ) tauhat
c     
      write ( iotty, '(5x,"eta = ", 1pe13.4," Ohm-m", 3x,
     $     "Thickness = ",1pe13.4," m",//,
     $     5x, "Eig. No", 5x, "Eigenvalue", 3x, "Time, msec." )' )
     $     eta_vess, th_vess
      write ( outmod,'(5x,"eta = ", 1pe13.4," Ohm-m", 3x,
     $     "Thickness = ",1pe13.4," m",//,
     $     5x, "Eig. No", 5x, "Eigenvalue", 3x, "Time, msec." )' )
     $     eta_vess, th_vess
c     
      do i = 1, jwal1
         if ( abs(work0(i)) .le. epszer )
     $        work0(i) = sign(epszer,work0(i))
         tc_vess(jwal1+i) = 1000.0 * tauhat / work0(i)
      end do
      do i = 1, jwal1
         write ( iotty,  '(5x, i4, 5x, 1p2e13.4)' )
     $        i, work0(i), tc_vess(jwal1+i)*twopi**2
         write ( outmod, '(5x, i4, 5x, 1p2e13.4)' )
     $        i, work0(i), tc_vess(jwal1+i)*twopi**2
      end do
c     
c..   Check diagonalization:
c     
      wrkb4 = transpose (wrkb3)
      wrkb4(1:jwal1,1:jwal1) =
     $     reshape((/ ((work0(i),i=1,jwal1),j=1,jwal1) /),
     $     (/ jwal1,jwal1 /)) * wrkb4(1:jwal1,1:jwal1)
      wbig2(1:jwal1,1:jwal1) =
     $     matmul(wrkb3(1:jwal1,1:jwal1),wrkb4(1:jwal1,1:jwal1))
      call matwrtn ( wbig2, jmax2,jmax2, 1,1,jwal1,jwal1,jdel,jdel,
     $     "D-C(lw,lw)???", outmod,0 )
c     
c..   Calculate the Current Eigenmodes. Store in wbig2:
c     
      write ( iotty,  '(/,5x,"The ODD Current Eigenmodes.")' )
      write ( outmod, '(/,5x,"The ODD Current Eigenmodes.")' )
      wbig2(1:jwal1,1:jwal1) =
     $     transpose ( wrkb3(1:jwal1,1:jwal1) / reshape
     $     ((/ ((sqrt(rwvlo(i)),i=1,jwal1),j=1,jwal1) /),
     $     (/ jwal1,jwal1 /)) )

      call matwrtn ( wbig2, jmax2,jmax2, 1,1,jwal1,jwal1,jdel,jdel,
     $     "ODD Current, J_(ls,lw)", outmod,0 )
c     
c..   Plot the First Few Eigencurrents
c     Assume Imaginary part of current coeff. is zero:
c     .  Use wkxx, wkyy as storage.
c     
c     First, shift variables so that theta=0 is on the inside.
c     Use workl, workl2 as storage.
c     
      workl(1:mth1)  = (/ xwal(mth/2+1:mth),xwal(1:mth/2+1) /)
      workl2(1:mth1) = (/ zwal(mth/2+1:mth),zwal(1:mth/2+1) /)
c     
      nbngr = 129
      do js = jwal1-nwcur+1, jwal1
         write ( outmod, '("workl,workl2,ell,wgrd,wkxx,wkyy=")')
         wkxx = 0.0
         wkyy = 0.0
         wkxx(1:mth1) = matmul ( wbig2(js,1:jwal1),
     $        transpose(rwvco(1:mth1,1:jwal1)) )

         wkxx(1:mth1) = (/ wkxx(mth/2+1:mth),wkxx(1:mth/2+1) /)
         wkyy(1:mth1) = (/ wkyy(mth/2+1:mth),wkyy(1:mth/2+1) /)

         iwrto = 0
         call kcur0 ( workl,workl2,ell,wgrd,mth1,wkxx,wkyy,n,
     $        nbngr,wkgrd,wxgrd,wzgrd,wkdxt,wkdzt,
     $        wkt2,wkp2,nths,0,"ODD_Eig-Cur",js, iwrto )
      end do
c     
c....................End of Odd Modes
c     
c     . Get Eigenmodes of Norm-Delta-C_(lw,lw) -> to ls basis.
c     .    Both Even and Odd together
c     
      WRITE ( OUTMOD, '(/,20x, "Eigen. of Norm-C(l_w,l_w.",/
     $     "To ls basis:")' )
      WRITE ( IOTTY,  '(/,20x, "Eigen. of Norm-C(l_w,l_w.",/
     $     "To ls basis:")' )
c     
      CALL mateig2a ( wbig, jmax2,jwal2, wrkb3,work0,
     $     1,jwal2, ff, 2, 2, jobid, "Norm-Del-C-eo",  outmod )
      CALL matwrtn ( wrkb3, jmax2,jmax2, 1,1,jwal2,jwal2,jdel,jdel,
     $     "T_(lw,lw) - evens and odds", outmod,0 )

c...  Calculate Time Constants for Eigenmodes of the Shell:
      WRITE ( IOTTY, '(/,10x,
     $     "Time Consts for Eigenmodes of the Shell:")' )
      WRITE ( OUTMOD, '(/,10x,
     $     "Time Consts. For Eigenmodes of the Shell:")' )
c     
      tauhat = amu0 * th_vess / eta_vess
      WRITE ( IOTTY, '(/, 5x,
     $     "Tau_hat = Mu*DeltaZ/Eta = ", 1pe13.4 )' ) tauhat
      WRITE ( OUTMOD, '(/, 5x,
     $     "Tau_hat = Mu*DeltaZ/Eta = ", 1pe13.4 )' ) tauhat
c     
      WRITE ( IOTTY, '(5x,"eta = ", 1pe13.4," Ohm-m", 3x,
     $     "Thickness = ",1pe13.4," m",//,
     $     5x, "Eig. No", 5x, "Eigenvalue", 3x, "Time, msec." )' )
     $     eta_vess, th_vess
      WRITE ( OUTMOD,'(5x,"eta = ", 1pe13.4," Ohm-m", 3x,
     $     "Thickness = ",1pe13.4," m",//,
     $     5x, "Eig. No", 5x, "Eigenvalue", 3x, "Time, msec." )' )
     $     eta_vess, th_vess
c     
      DO i = 1, jwal2
         IF ( abs(work0(i)) .le. epszer )
     $        work0(i) = sign(epszer,work0(i))
         tc_vess(i) = 1000.0 * tauhat / work0(i)
      END DO
      DO i = 1, jwal2
         WRITE ( IOTTY,  '(5x, i4, 5x, 1p2e13.4)' )
     $        i, work0(i), tc_vess(i)
      END DO

c     . Also expand the C(w_i,p) terms of the shell functions to C(l_w,l_p).
c     They are stored in GRRI after MTH rows.

      WRITE ( IOTTY,  '(/,"Expand C_lp(w_i) terms of shell functions",/,
     $     "to C_(l_w,l_p).")' )
      WRITE ( OUTMOD, '(/,"Expand C_lp(w_i) terms of shell functions",/,
     $     "to C_(l_w,l_p).")' )
c     
c...  Replace with Fortran 90:
c$$$  c88888888888888888888888888888888888888888888888888888888
c$$$  c... Careful here! Using jwal1 wall modes:
c$$$  c8888888888888888888888888888888888888888888888888888888
c     
      zzfac = dth * twopi * tpisqi
      wprr(1:jwal1,1:jmax1) = zzfac *
     $     matmul ( transpose(rwvce(1:mth,1:jwal1)),
     $     grri(mth+1:2*mth,1:jmax1) )                  ! even-cos
      wpii(1:jwal1,1:jmax1) = zzfac *
     $     matmul ( transpose(rwvco(1:mth,1:jwal1)),
     $     grri(mth+1:2*mth,jmax1+1:2*jmax1) )          ! odd-sin
      wpir(1:jwal1,1:jmax1) = zzfac *
     $     matmul ( transpose(rwvco(1:mth,1:jwal1)),
     $     grri(mth+1:2*mth,1:jmax1) )                  ! odd-cos
      wpri(1:jwal1,1:jmax1) = zzfac *
     $     matmul ( transpose(rwvce(1:mth,1:jwal1)),
     $     grri(mth+1:2*mth,jmax1+1:2*jmax1) )          ! even-sin
c     
      call matwrtn ( wprr, jmax1,jmax1, 1,ln,jwal1,17,8,8,
     $     "C_(lw,lp)-rr", outmod,0 )
      call matwrtn ( wpii, jmax1,jmax1, 1,ln,jwal1,17,8,8,
     $     "C_(lw,lp)-ii", outmod,0 )
c$$$  call matwrtn ( wprr, jmax1,jmax1, 1,ln,jwal1,jmax1,jdel,jdel,
c$$$  $     "C_(lw,lp)-rr", outmod,0 )
c$$$  call matwrtn ( wpii, jmax1,jmax1, 1,ln,jwal1,jmax1,jdel,jdel,
c$$$  $     "C_(lw,lp)-ii", outmod,0 )

      call matwrtn ( wpri, jmax1,jmax1, 1,ln,jwal1,jmax1,jdel,jdel,
     $     "C_(lw,lp)-ri", outmod,0 )
      call matwrtn ( wpir, jmax1,jmax1, 1,ln,jwal1,jmax1,jdel,jdel,
     $     "C_(lw,lp)-ir", outmod,0 )
c     
c     . Add rr and ii parts:
c$$$  do j1 = 1, jwal1
c$$$  do j2 = 1, jwal1
c$$$  wprr(j1,j2) = wprr(j1,j2) + wpii(j1,j2)
c$$$  end do
c$$$  end do
c     
      write ( outmod, '(/,"Store C_(lw,lp) in 4x4 matrix and left",/,
     $     "multiply by sqrt(evalues)" )')
c     
c     Storage:  | Even-Cos    - Even-Sin |
c               | Odd-Cos     - Odd-Sin |
c     
      do i = 1, jmax1
         wbig(1:jwal1,i) =
     $        sqrt(rwvle(1:jwal1)) * wprr(1:jwal1,i)    !   even-cos
         wbig(jwal1+1:jwal2,jmax1+i) =
     $        - sqrt(rwvlo(1:jwal1)) * wpii(1:jwal1,i)  ! - odd-sin
         wbig(1:jwal1,jmax1+i) =
     $        - sqrt(rwvle(1:jwal1)) * wpri(1:jwal1,i)  ! - even-sin
         wbig(jwal1+1:jwal2,i) =
     $        sqrt(rwvlo(1:jwal1)) * wpir(1:jwal1,i)    !   odd-cos
      end do
c     
      call matwrtn ( wbig,jmax2,jmax2, 1,1,jwal2,2*jmax1,jdel,jdel,
     $     "Norm Big C_(lw,lp)", outmod,0 )
      call matwrts ( wbig,jmax2,jmax2, 1,1, 9,9, 9,9,
     $     1,1, "Norm Big C_(lw,lp)", outmod,0 )
      call matwrts ( wbig,jmax2,jmax2, 10,10, 18,18, 9,9,
     $     10,10, "Norm Big C_(lw,lp)", outmod,0 )
      call matwrts ( wbig,jmax2,jmax2, 1,10, 9,18, 9,9,
     $     1,10, "Norm Big C_(lw,lp)", outmod,0 )
      call matwrts ( wbig,jmax2,jmax2, 10,1, 18,9, 9,9,
     $     10,1, "Norm Big C_(lw,lp)", outmod,0 )

c     . Store this in oplp(l_w,l_p) for the open loop matrix:
      oplp(jmax2+1:jmax2+jwal2,1:jmax2) = wbig(1:jwal2,1:jmax2)
c$$$      oplp(jmax1+1:jmax1+jwal1,1:jmax1) = wbig(1:jwal1,1:jmax1)
c$$$      oplp(jmax1+jwal1+1:jmax1+jwal2,1:jmax1) =
c$$$     $     wbig(jwal1+1:jwal2,jmax1+1:jmax2)
c     
c     . Transform C_(lw,lp) to C_(ls,lp). Use WPII for storage
c     Then transform WPRR to WPII
c     
      write ( outmod,
     $     '(/,20x,"Get transpose of the eigenfunction matrix.",/,
     $     "Transform C_(lw,lp) to C_(ls,lp)" )' )
      write ( iotty,
     $     '(/,20x,"Get transpose of the eigenfunction matrix.",/,
     $     "Transform C_(lw,lp) to C_(ls,lp)" )' )
c     
c     . Get transpose of the eigenfunction matrix, WRKB3. Put into WRKB4.
c     Then transform WBIG to WBIG2

      wrkb4 = transpose ( wrkb3 )
      wbig2(1:jwal2,1:2*jmax1) =
     $     matmul ( wrkb4(1:jwal2,1:jwal2),wbig(1:jwal2,1:2*jmax1) )
      CALL matwrtn ( wbig2, jmax2,jmax2, 1,1,jwal2,2*jmax1,jdel,jdel,
     $     "C_(ls,lp)", outmod,0 )
c     
c..   oooooo   IV-149(6):   Get C(lp,ls) ooooooooooooooooooooooooooooooooo
c     
c     Right multiply \bar C(lp,lw) by T:
c     
      wrkb5(1:jmax2,1:jwal2) =
     $     matmul ( wrkb5(1:jmax2,1:jwal2),wrkb3(1:jwal2,1:jwal2) )
      CALL matwrtn ( wrkb5, jmax2,jmax2, 1,1,jmax2,jwal2,jdel,jdel,
     $     "C_(lp,ls)", outmod,0 )
c..   ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
c     
c     IV-150(1):  Get contribution to Chi-plasma from the wall:  ooooooo
c     Store in wrkb5:
c     
c     .  Normalize such that taugamma is dimensionless, i.e.
c     to the lowest eigenvalue of M_ls, which is:

      amlsmin = minval ( work0(1:jwal2) )
      tgmlsmn = taugam * amlsmin
      WRITE ( iotty,   '(/,5x, "amlsmin = ", es13.5,
     $     " tgmlsmn = ", es13.5 )' ) amlsmin, tgmlsmn
      WRITE ( outmod,  '(/,5x, "amlsmin = ", es13.5,
     $     " tgmlsmn = ", es13.5 )' ) amlsmin, tgmlsmn
c     
c..   First, divide C(ls,lp) by tau*gamma + M:
c     Use wbig for storage.
c     
c     taugam = 10.0

      DO jp = 1, jmax2
         wbig(1:jwal2,jp) =
     $        wbig2(1:jwal2,jp) / ( tgmlsmn + work0(1:jwal2) )
      END DO
c     
      wrkb5(1:jmax2,1:jmax2) =
     $     matmul ( wrkb5(1:jmax2,1:jwal2),wbig(1:jwal2,1:jmax2) )
      CALL matwrtn ( wrkb5, jmax2,jmax2, 1,1,jmax2,jmax2,jdel,jdel,
     $     "C_(lp,ls)C_(ls,lp)", outmod,0 )
      CALL matwrtx ( wrkb5,jmax2,jmax2,1,8,1,8,
     $     "C_(lp,ls)C_(ls,lp)", outmod,iotty )
c     
      CALL vacasym ( wrkb5, jmax2,jmax2,"C_(lp,ls)C_(ls,lp)",
     $     outmod,iotty )
c     
c     .  Store in v_cpwwp
      v_cpwwp(1:jmax2,1:jmax2) = wrkb5(1:jmax2,1:jmax2)
c     
c     
c     oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
c     
c..   Now transform the coil into ls space:
c     
      coil00(1:jwal2) =
     $     matmul ( wrkb4(1:jwal2,1:jwal2),coil00(1:jwal2) )
      CALL vecwrt ( jwal2, coil00, "coil00 in ls basis", 1,jwal2,
     $     outmod, iotty )
c     
      WRITE ( outmod,
     $     '(/,"Solve for B_ls. Easy, since coeff. is diagonal.")' )
      WRITE ( iotty,
     $     '(/,"Solve for B_ls. Easy, since coeff. is diagona.l")' )
c     
c.....Put in WBIG. Take care of vanishing eigenvalues:
c     
      WHERE ( abs (work0(1:jwal2)) <= epszer )
         work0(1:jwal2) = sign ( epszer, work0(1:jwal2) )
      END WHERE

      DO jp = 1, 2*jmax1
         wbig(1:jwal2,jp) = wbig2(1:jwal2,jp) / work0(1:jwal2)
      END DO

c..   Store B_ls in blsp for use in subroutine blsp:
      blsp(1:jwal2,1:2*jmax1) = wbig(1:jwal2,1:2*jmax1)
c     
c$$$  do js = 1, jwal2
c$$$  do jp = 1, 2*jmax1
c$$$  if ( abs(work0(js)) .le. epszer )
c$$$  $        work0(js) = sign(epszer,work0(js))
c$$$  wbig(js,jp) = wbig2(js,jp)/work0(js)
c$$$  end do
c$$$  end do
c     
      CALL matwrtn ( wbig, jmax2,jmax2, 1,1,jwal2,2*jmax1,jdel,jdel,
     $     "B_(ls,lp)", outmod,0 )
c     
      WRITE ( outmod, '("Take B_ls.back to lw representation.")' )
      WRITE ( iotty,  '("Take B_ls.back to lw representation.")' )
c     
      wrkb4(1:jwal2,1:2*jmax1) =
     $     matmul (wrkb3(1:jwal2,1:jwal2),wbig(1:jwal2,1:2*jmax1) )
      call matwrtn ( wrkb4, jmax2,jmax2, 1,1,jwal2,2*jmax1,jdel,jdel,
     $     "bar-B_(lw,lp)", outmod,0 )
c     
c..   Unnormalize from sqrt. eigenvalues:

      DO jp = 1, 2*jmax1
         wbig(1:jwal1,jp) = wrkb4(1:jwal1,jp) * sqrt(rwvle(1:jwal1))
         wbig(jwal1+1:jwal2,jp) = wrkb4(jwal1+1:jwal2,jp)
     $        * sqrt(rwvlo(1:jwal1))
      END DO

      CALL matwrtn ( wbig, jmax2,jmax2, 1,1,jwal2,2*jmax1,jdel,jdel,
     $     "B_(lw,lp)", outmod,0 )
c     
      WRITE ( outmod, '("Get Current Response at Large Time.")' )
      WRITE ( iotty, '("Get Current Response at Large Time.")' )
c     
      DO jp = 1, 2*jmax1
         wbig2(1:jwal1,jp) = wbig(1:jwal1,jp) / rwvle(1:jwal1)
         wbig2(jwal1+1:2*jwal1,jp) = wbig(jwal1+1:jwal2,jp)
     $        / rwvlo(1:jwal1)
      END DO

      CALL matwrtn ( wbig2, jmax2,jmax2, 1,1,jwal2,2*jmax1,jdel,jdel,
     $     "I_(lw,lp)", outmod,0 )
c     
c..   Construct Wall Response Matrix to Plasma Modes.
c     .   Add Even and Odd Contributions:
c     
c...  Store in GWIN, GWOT
c     
      gwin(1:mth1,1:jmax1) =
     $     matmul ( rwvce(1:mth1,1:jwal1),
     $     wbig2(1:jwal1,1:jmax1) ) +
     $     matmul ( rwvco(1:mth1,1:jwal1),
     $     wbig2(jwal1+1:2*jwal1,1:jmax1) )
      gwot(1:mth1,1:jmax1) =
     $     matmul ( rwvce(1:mth1,1:jwal1),
     $     wbig2(1:jwal1,jmax1+1:2*jmax1) ) +
     $     matmul ( rwvco(1:mth1,1:jwal1),
     $     wbig2(jwal1+1:2*jwal1,jmax1+1:2*jmax1) )
c     
      CALL matwrtn ( gwin,nths2,nwm2, 1,1,mth1,jmax1,jdel,jdel,
     $     "GWIN_(i,lp)", outmod,0 )
      CALL matwrtn ( gwot,nths2,nwm2, 1,1,mth1,jmax1,jdel,jdel,
     $     "GWOT_(i,lp)", outmod,0 )
      WRITE ( outmod, '(/, "Get Response matrix for B on Wall",/,
     $     "as B_lp(theta_w)" )')
      WRITE ( iotty,  '(/, "Get Response matrix for B on Wall",/,
     $     "as B_lp(theta_w)" )')
c     
      DEALLOCATE ( wprr, wpii, wpri, wpir,
     $     zrr, zii, zri, zir, wbig, wbig2, wrkb3, wrkb4, wrkb5 )

c..   Model for coil current:
c     
      nbngr = 129
c     
      RETURN
      END SUBROUTINE resshel

c.....................................................................
      SUBROUTINE resshel_1 ( v_cpwwp, nd1 )
c.................................................................
c     
      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum5.inc'
      INCLUDE 'vacuum7.inc'
c     
      REAL, DIMENSION(:,:), ALLOCATABLE :: zwk1
      REAL, DIMENSION(:), ALLOCATABLE :: zwk2
      REAL, DIMENSION(nd1,nd1) :: v_cpwwp
      CHARACTER(132) :: zlabel
      CHARACTER(24) :: datimv
      LOGICAL backl
c     
c      datimv = datev(1:10)//","//timev(1:10)
      datimv(1:24) = date_array(1:24)

      ln = lmin(1)
      lx = lmax(1)
      jmax1 = lx - ln + 1
      jmax2 = 2*jmax1
      backl = .false.
c     
c..   Transform cpp to finite elements basis.
c     
      lenj = INDEX(jobid, '  ', backl ) - 1
      zlabel =
     $     datimv//": C_PP-fel, ID= "//jobid(1:lenj)//
     $     ", "//filein2(1:lenif)

      CALL new_basis ( vacmat0, nfm, mfel,jmax1, ltfil,
     $     sinlt,coslt, nths,nfm, outmod,iotty )
      CALL matwrtn ( vacmat0,nfm,nfm, 1,1,mfel,mfel,jdel,jdel,
     $     "C_pp-fel", outmod,0 )

      CALL gatonorm ( vacmat0, nfm, gatovac,nfm, rgato,mfel,mth,
     $     qa1,twopi,taugam, vtog, zlabel )

      WRITE ( outmod, '("Conversion factor = ", 1pe12.5)') vtog
      WRITE ( iotty,  '("Conversion factor = ", 1pe12.5)') vtog

      CALL matwrtn ( gatovac,nfm,nfm,1,1,mfel,mfel,8,8,
     $     "Norm VACMAT0-fel in RESSHEL_1", outmod,iotty )

      ndz = mfel
      leigen = 1

      ALLOCATE ( zwk1(nfm,nfm), zwk2(nfm) )

      call mateig2a ( gatovac,nfm, ndz, zwk1,zwk2,
     $     ln,lx, ff, 0, leigen, jobid, "GatoVac-fel", outmod )

c     Now the V_CPWWP Matrix with Taugam dependence:
c     Combine RR and II parts of V_CPWWP:

      v_cpwwp(1:jmax1,1:jmax1) = v_cpwwp(1:jmax1,1:jmax1) +
     $     v_cpwwp(jmax1+1:jmax2,jmax1+1:jmax2)

c..   Transform v_cpwwp to finite elements basis.

      zlabel =
     $     datimv//": C_PWWPP-fel, ID= "//jobid(1:lenj)//
     $     ", "//filein2(1:lenif)
      lenzl = INDEX(zlabel, "   " ) - 1
      write ( outmod,'(/,a)' ) zlabel(1:lenzl)
      write ( iotty,'(/,a)' ) zlabel(1:lenzl)

      call atpoint ( "v_cpwwp", "ndz",ndz,"vtog",
     $     vtog, iotty, outmod )

      CALL new_basis ( v_cpwwp, nd1, mfel,jmax1, ltfil,
     $     sinlt,coslt, nths,nfm, outmod,iotty )

      CALL matwrtn ( v_cpwwp,nd1,nd1, 1,1,mfel,mfel,jdel,jdel,
     $     "C_PWWP-fel", outmod,0 )

      CALL gatonorm ( v_cpwwp, nd1, gatovac,nfm, rgato,mfel,mth,
     $     qa1,twopi,taugam, vtog, zlabel )

      WRITE ( outmod, '("Conversion factor = ", 1pe12.5)') vtog
      WRITE ( iotty,  '("Conversion factor = ", 1pe12.5)') vtog

      CALL matwrtn ( gatovac,nfm,nfm,1,1,mfel,mfel,8,8,
     $     "Norm C_PWWP-fel in RESSHEL_1", outmod,iotty )

      call mateig2a ( gatovac,nfm, ndz, zwk1,zwk2,
     $     ln,lx, ff, 0, leigen, jobid, "Norm C_pwwp-fel", outmod )

c     Now add CPP to CPWWP:

      zlabel =
     $     datimv//": (C_pp + C_pwwp)-fel, ID= "//jobid(1:lenj)//
     $     ", "//filein2(1:lenif)

      v_cpwwp(1:mfel,1:mfel) = - v_cpwwp(1:mfel,1:mfel) +
     $     vacmat0(1:mfel,1:mfel)

      CALL matwrtn ( v_cpwwp,nd1,nd1, 1,1,mfel,mfel,jdel,jdel,
     $     "(C_PP + C_PWWP)-fel", outmod,0 )

      CALL gatonorm ( v_cpwwp, nd1, gatovac,nfm, rgato,mfel,mth,
     $     qa1,twopi,taugam, vtog, zlabel )

      WRITE ( outmod, '("Conversion factor = ", 1pe12.5)') vtog
      WRITE ( iotty,  '("Conversion factor = ", 1pe12.5)') vtog

      CALL matwrtn ( gatovac,nfm,nfm,1,1,mfel,mfel,8,8,
     $     "Norm (C_PP+C_PWWP)-fel in RESSHEL_1", outmod,iotty )
      CALL matwrtx ( gatovac,nfm,nfm,1,8,1,8,
     $     "Norm (C_PP+C_PWWP)-fel in RESSHEL_1", outmod,iotty )
      CALL vacasym ( gatovac, nfm,mfel,
     $     "Norm (C_PP+C_PWWP)-fel in RESSHEL_1",
     $     outmod,iotty )

      CALL mateig2a ( gatovac,nfm, ndz, zwk1,zwk2,
     $     ln,lx, ff, 0, leigen, jobid, "C_pp+C_pwwp", outmod )


      DEALLOCATE (zwk1, zwk2 )

      END SUBROUTINE resshel_1

c......................................................
      SUBROUTINE resshel_2
c......................................................
c     
      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum5.inc'
      INCLUDE 'vacuum7.inc'
      INCLUDE 'vacuum8.inc'
c     
      REAL, DIMENSION(:,:), ALLOCATABLE :: zwk1
      REAL, DIMENSION(:,:), ALLOCATABLE :: zbwk,zwk4
      REAL, DIMENSION(:), ALLOCATABLE :: zwk2, zwk3, bwze, bwzo
      REAL, DIMENSION(:), ALLOCATABLE :: bpzc, bpzs, bpzoc,bpzos
      REAL, DIMENSION(:), ALLOCATABLE :: kwze, kwzo, xibpe, xize, xizo
      REAL, DIMENSION(:), ALLOCATABLE :: z1tmp, z2tmp , xibpo
      REAL, DIMENSION(:), ALLOCATABLE :: xibwe, xibwo
      REAL, DIMENSION(:), ALLOCATABLE :: chiplasr, chiplasi
      REAL, DIMENSION(:), ALLOCATABLE :: chiwallr, chiwalli
      REAL, DIMENSION(:), ALLOCATABLE :: b_thet_w_r, b_thet_w_i,
     $     b_phi_w_r, b_phi_w_i
      REAL, DIMENSION(:), ALLOCATABLE :: flbpe, flbpo,flbwe,flbwo
      REAL, DIMENSION(:), ALLOCATABLE :: fltte, fltto,fltamp,fltang
      REAL, DIMENSION(:), ALLOCATABLE :: exo, exe, examp, exang
      REAL, DIMENSION(:), ALLOCATABLE :: exawo, exawe, exawamp, exawang
      REAL, DIMENSION(:), ALLOCATABLE :: exbwo, exbwe, exbwamp, exbwang
      REAL, DIMENSION(:), ALLOCATABLE :: exaps, exapc, exapamp, exapang
      REAL, DIMENSION(:), ALLOCATABLE :: exbps, exbpc, exbpamp, exbpang
      REAL, DIMENSION(:), ALLOCATABLE :: dwp,dwe,dwo,dwpe,dwpo
      DIMENSION slx(501), slxt(501), slz(501), slzt(501)
c$$$      REAL, DIMENSION(:), ALLOCATABLE :: slx, slxt, slz, slzt

      CHARACTER(132) :: zlabel
      CHARACTER(24) :: datimv
c     

      ioop = 66
      OPEN (  ioop, FILE='OPLP-eigen', STATUS='REPLACE',
     $     FORM='FORMATTED' )
c$$$
      iosen = 28
c$$$      OPEN (  iosen, FILE='SENSOR-MEAS', STATUS='REPLACE',
c$$$     $     FORM='FORMATTED' )
      
      datimv(1:24) = date_array(1:24)
c     datimv = datev(1:10)//","//timev(1:10)

      ln = lmin(1)
      lx = lmax(1)
      jmax1 = lx - ln + 1
      jmax2 = 2*jmax1
      jwal3 = 3*jwal1           ! put in common later.

c..   Fill in plasma and vacuum matrices in oplp(l_p,l_p).
c     Use half of each for RR and II for now. Zfac is for whether
c     we need the full complement or only half in each block.
c...  12-14-05: No uncertainty now for up-down asymmetry. Fill in 
c     all jmax2 by jmax2 elements for the plasma-plasma block

c..   .The ww. pw, wp parts come from Subroutine RESSHEL. Approximate lines:
c     ww at line 229, pw at line 299, wp at line 638.

      zfac = 1.0
      oplp(1:jmax1,1:jmax1) = zfac * ( wpp0(1:jmax1,1:jmax1) +
     $     vacmat00(1:jmax1,1:jmax1) )

      oplp(jmax1+1:jmax2,jmax1+1:jmax2) = zfac *
     $     ( wpp0(1:jmax1,1:jmax1) + vacmat00(1:jmax1,1:jmax1) )

c..   Experiment with sign of vacmt00i: 
      oplp(1:jmax1,jmax1+1:jmax2) = 
     $     - ( wpp0i(1:jmax1,1:jmax1) + vacmt00i(1:jmax1,1:jmax1) )
     $     * zfac

      oplp(jmax1+1:jmax2,1:jmax1) =
     $     ( wpp0i(1:jmax1,1:jmax1) + vacmt00i(1:jmax1,1:jmax1) )
     $     * zfac

      ln = 1
      lx = jmax2 + jwal4        ! For up-down-asymmetry.
      ndpw2 = nfm2 + nwm2
      ndpw4 = nfm2 + nwm4       ! For up-down-asymmetry.
      nspw2 = jmax2 + jwal2     ! For up-down asymmetry.
      nspw4 = jmax2 + jwal4     ! For up-down asymmetry.

c$$$  call matwrts ( oplp,ndpw4,ndpw4, 19,1, 27,9, 9,9,
c$$$  $     19,1,"OPLP(lplw,lplw)", outmod,0 )
c$$$  call matwrts ( oplp,ndpw4,ndpw4, 28,10, 36,18, 9,9,
c$$$  $     28,10,"OPLP(lplw,lplw)", outmod,0 )
c$$$  call matwrts ( oplp,ndpw4,ndpw4, 19,10, 27,18, 9,9,
c$$$  $     19,10,"OPLP(lplw,lplw)", outmod,0 )
c$$$  call matwrts ( oplp,ndpw4,ndpw4, 28,1, 36,9, 9,9,
c$$$  $     28,1,"OPLP(lplw,lplw)", outmod,0 )

      ALLOCATE ( zwk1(ndpw4,ndpw4), zwk2(nspw4), zwk3(nspw4) )
      ALLOCATE ( zbwk(nspw4,nspw4) )

c..   The RHS norm is 0.0 for the jmax2 plasma matrix and jwal4 for the wall.

      zbwk = 0.0
      epszb = 1.0E-10
      DO i = 1, jmax2
         zbwk(i,i) = epszb
      END DO

c$$$  DO i = 1, jwal1
c$$$  c     zbwk(jmax2+i,jmax2+i) = 1.0/rwvle(i)
c$$$  c     zbwk(jmax2+jwal1+i,jmax2+jwal1+i) = 1.0/rwvlo(i)
c$$$  zbwk(jmax2+i,jmax2+i) = 1.0
c$$$  zbwk(jmax2+jwal1+i,jmax2+jwal1+i) = 1.0
c$$$  END DO

c...  Now the norm is 4*jwal1 = jwal4 for the up-down case.
      DO i = 1, jwal4
         zbwk(jmax2+i,jmax2+i) = 1.0
      END DO

c...  wp and pw blocks:
      DO i = 1, jmax2
         oplp(jmax2+1:jmax2+jwal1,i) = 1.0 *
     $        oplp(jmax2+1:jmax2+jwal1,i)
         oplp(i,jmax2+1:jmax2+jwal1) = 1.0 *
     $        oplp(i,jmax2+1:jmax2+jwal1)
         oplp(jmax2+jwal1+1:jmax2+jwal2,i) = 1.0 *
     $        oplp(jmax2+jwal1+1:jmax2+jwal2,i)
         oplp(i,jmax2+jwal1+1:jmax2+jwal2) = 1.0 *
     $        oplp(i,jmax2+jwal1+1:jmax2+jwal2)
      END DO

c..   ww blocks;
      DO i = jmax2+1, jmax2+jwal2
         oplp(i,jmax2+1:jmax2+jwal2) = 1.0 *
     $        oplp(i,jmax2+1:jmax2+jwal2)
      END DO
c$$$  DO i = 1, jmax1
c$$$  oplp(jmax1+1:jmax1+jwal1,i) = 1./sqrt(rwvle(1:jwal1)) *
c$$$  $        oplp(jmax1+1:jmax1+jwal1,i)
c$$$  oplp(i,jmax1+1:jmax1+jwal1) = 1./sqrt(rwvle(1:jwal1)) *
c$$$  $        oplp(i,jmax1+1:jmax1+jwal1)
c$$$  oplp(jmax1+jwal1+1:jmax1+jwal2,i) = 1./sqrt(rwvlo(1:jwal1)) *
c$$$  $        oplp(jmax1+jwal1+1:jmax1+jwal2,i)
c$$$  oplp(i,jmax1+jwal1+1:jmax1+jwal2) = 1./sqrt(rwvlo(1:jwal1)) *
c$$$  $        oplp(i,jmax1+jwal1+1:jmax1+jwal2)
c$$$  END DO
c$$$  DO i = 1, jwal2
c$$$  oplp(jmax1+1:jmax1+jwal2,jmax1+i) = 1./
c$$$  $        (/ sqrt(rwvle(1:jwal1)),sqrt(rwvlo(1:jwal1)) /) *
c$$$  $        oplp(jmax1+1:jmax1+jwal2,jmax1+i)
c$$$  oplp(jmax1+i,jmax1+1:jmax1+jwal2) = 1./
c$$$  $        (/ sqrt(rwvle(1:jwal1)),sqrt(rwvlo(1:jwal1)) /) *
c$$$  $        oplp(jmax1+i,jmax1+1:jmax1+jwal2)
c$$$  c$$$
c$$$  c$$$         oplp(jmax1+jwal1+1:jmax1+jwal2,jmax1+i) =
c$$$  c$$$     $        sqrt(rwvlo(1:jwal1)) *
c$$$  c$$$     $        oplp(jmax1+jwal1+1:jmax1+jwal2,jmax1+i)
c$$$  c$$$         oplp(jmax1+i,jmax1+jwal1+1:jmax1+jwal2) =
c$$$  c$$$     $        sqrt(rwvlo(1:jwal1)) *
c$$$  c$$$     $        oplp(jmax1+i,jmax1+jwal1+1:jmax1+jwal2)
c$$$  END DO

c...  Now include the extra plasma-wall and wall terms for up-down asymmetry:
c     OPLP now initialized to 0.0 in Main Program

      DO i = 1, jmax1
         oplp (jmax2+jwal2+1:jmax2+jwal4,i) =
     $        - oplp(jmax2+1:jmax2+jwal2,jmax1+i)
         oplp (jmax2+jwal2+1:jmax2+jwal4,jmax1+i) =
     $        oplp(jmax2+1:jmax2+jwal2,i)
         oplp (i,jmax2+jwal2+1:jmax2+jwal4) =
     $        - oplp(jmax1+i,jmax2+1:jmax2+jwal2)
         oplp (jmax1+i,jmax2+jwal2+1:jmax2+jwal4) =
     $        oplp(i,jmax2+1:jmax2+jwal2)
      END DO
c..   Fiallly, the ww diagonal terms.      
      DO i = 1, jwal2
         oplp(jmax2+jwal2+1:jmax2+jwal4, jmax2+jwal2+i) =
     $        oplp(jmax2+1:jmax2+jwal2, jmax2+i)
      END DO

      call matwrtn ( oplp,ndpw4,ndpw4, 1,1,nspw4,nspw4,9,9,
     $     "OPLP", outmod,0 )
      call vacasym ( oplp, ndpw4,nspw4,"OPLP", outmod,iotty )

      call matwrtn ( oplp,ndpw4,ndpw4, 1,1,nspw4,nspw4,9,9,
     $     "OPLP(lplw,lplw)", outmod,0 )

c     call matwrts ( oplp,ndpw4,ndpw4, 1,1, 9,9, 9,9,
c     $     1,1,"OPLP(lplw,lplw)", outmod,0 )
c     call matwrts ( oplp,ndpw4,ndpw4, 37,1, 45,10, 9,9,
c     $     37,1,"OPLP(lplw,lplw)", outmod,0 )
c     call matwrts ( oplp,ndpw4,ndpw4, 1,37, 9,45, 9,9,
c     $     1,37,"OPLP(lplw,lplw)", outmod,0 )
c     call matwrts ( oplp,ndpw4,ndpw4, 37,37, 45,45, 9,9,
c     $     37,37,"OPLP(lplw,lplw)", outmod,0 )
      call matwrts ( oplp,ndpw4,ndpw4, 10,10, 18,18, 9,9,
     $     10,10,"OPLP(lplw,lplw)", outmod,0 )
      call matwrts ( oplp,ndpw4,ndpw4, 19,19, 27,27, 9,9,
     $     19,19,"OPLP(lplw,lplw)", outmod,0 )
c     call matwrts ( oplp,ndpw4,ndpw4, 35,35, 43,43, 9,9,
c     $     35,35,"OPLP(lplw,lplw)", outmod,0 )
      call matwrts ( oplp,ndpw4,ndpw4, 28,28, 36,36, 9,9,
     $     28,28,"OPLP(lplw,lplw)", outmod,0 )
c     call matwrts ( oplp,ndpw4,ndpw4, 55,55, 63,63, 9,9,
c     $     55,55,"OPLP(lplw,lplw)", outmod,0 )
      call matwrts ( oplp,ndpw4,ndpw4, 48,48, 56,56, 9,9,
     $     48,48,"OPLP(lplw,lplw)", outmod,0 )
c     call matwrts ( oplp,ndpw4,ndpw4,  9,35, 17,43, 9,9,
c     $     9,35,"OPLP(lplw,lplw)", outmod,0 )
      call matwrts ( oplp,ndpw4,ndpw4,  9,28, 17,36, 9,9,
     $     9,28,"OPLP(lplw,lplw)", outmod,0 )
c     call matwrts ( oplp,ndpw4,ndpw4, 17,35, 25,43, 9,9,
c     $     17,35,"OPLP(lplw,lplw)", outmod,0 )
      call matwrts ( oplp,ndpw4,ndpw4, 17,28, 25,36, 9,9,
     $     17,28,"OPLP(lplw,lplw)", outmod,0 )
c     call matwrts ( oplp,ndpw4,ndpw4,  9,55, 17,63, 9,9,
c     $      9,55,"OPLP(lplw,lplw)", outmod,0 )
      call matwrts ( oplp,ndpw4,ndpw4,  9,48, 17,56, 9,9,
     $     9,48,"OPLP(lplw,lplw)", outmod,0 )
c     call matwrts ( oplp,ndpw4,ndpw4, 17,55, 25,63, 9,9,
c     $     17,55,"OPLP(lplw,lplw)", outmod,0 )
      call matwrts ( oplp,ndpw4,ndpw4, 17,48, 25,56, 9,9,
     $     17,48,"OPLP(lplw,lplw)", outmod,0 )


      leval = 1
c     
      ALLOCATE (dwp(1:nspw4))
      ALLOCATE (dwe(1:nspw4))
      ALLOCATE (dwo(1:nspw4))
      ALLOCATE (dwpe(1:nspw4))
      ALLOCATE (dwpo(1:nspw4))
c     
c...  Force Symmetry on OPLP:

      WRITE ( outmod, '( "OPLP is symmetrized " )' )
      WRITE ( iotty,  '( "OPLP is symmetrized " )' )

      DO i = 2, jmax2 + jwal4
         DO j = 1, i-1
            oplpav = (oplp(i,j) + oplp(j,i))/2.0
            oplp(i,j) = oplpav
            oplp(j,i) = oplpav
         END DO
      END DO

      CALL mateig2b ( oplp,ndpw4, zbwk,nspw4, nspw4, epszb, zwk1,
     $     zwk2,zwk3, ln,lx,ff, 2, leval, jobid, "OPLP Eigs",
     $     outmod ,ell(mth1), jwal2, dwp, dwe, dwo, dwpe , dwpo )

      ALLOCATE ( zwk4(ndpw4,ndpw4) )
      ALLOCATE ( bwze(mth1),bwzo(mth1) )
      ALLOCATE ( bpzc(mth1),bpzs(mth1),bpzoc(mth1),bpzos(mth1) )
      ALLOCATE ( kwze(mth1),kwzo(mth1) )
      ALLOCATE ( chiplasr(mth1),chiplasi(mth1) )
      ALLOCATE ( z1tmp(mth1),z2tmp(mth1) )
      ALLOCATE ( xibpe(jmax1),xibpo(jmax1),xize(jmax1),xizo(jmax1) )
      ALLOCATE ( xibwe(jwal2),xibwo(jwal2) )
      ALLOCATE ( flbpe(lrswcomp),flbpo(lrswcomp),flbwe(lrswcomp))
      ALLOCATE ( flbwo(lrswcomp),fltte(lrswcomp),fltto(lrswcomp))
      ALLOCATE ( exo(lrswcomp)  ,exe(lrswcomp)  )
      ALLOCATE ( fltamp(lrswcomp)  ,fltang(lrswcomp)  )
      ALLOCATE ( examp(lrswcomp),exang(lrswcomp))
      ALLOCATE ( exawo(lrswcomp), exawe(lrswcomp),
     $     exawamp(lrswcomp), exawang(lrswcomp) )
      ALLOCATE ( exbwo(lrswcomp), exbwe(lrswcomp),
     $     exbwamp(lrswcomp), exbwang(lrswcomp) )
      ALLOCATE ( exaps(lrswcomp), exapc(lrswcomp),
     $     exapamp(lrswcomp), exapang(lrswcomp) )
      ALLOCATE ( exbps(lrswcomp), exbpc(lrswcomp),
     $     exbpamp(lrswcomp), exbpang(lrswcomp) )

      WRITE ( outmod, '( "Coil Excitation: " )' )
      WRITE ( iotty,  '( "Coil Excitation: " )' )
c     
c     note that at this point the zwk1 eigenfunctions are the coefficients
c     for the modified wall eigen functions
c     (rwvce, rwvco) sqrt(rwvle, rwvlo) / sqrt(2.)
c     therefore to make matters simple, we sould like to transform
c     to the base (rwvce, rwvco)
c     then the magnetic field amplitudes are stored in zwk4 and the
c     current potentials are stored in zwk1
c     
      DO i = 1, nspw4
         zwk4(jmax2+1:jmax2+jwal1,i) = zwk1(jmax2+1:jmax2+jwal1,i)
     $        * sqrt(rwvle(1:jwal1))
         zwk4(jmax2+jwal1+1:jmax2+jwal2,i) =
     $        zwk1(jmax2+jwal1+1:jmax2+jwal2,i)
     $        * sqrt(rwvlo(1:jwal1))
         zwk1(jmax2+1:jmax2+jwal1,i) = zwk1(jmax2+1:jmax2+jwal1,i)
     $        / sqrt(rwvle(1:jwal1))
         zwk1(jmax2+jwal1+1:jmax2+jwal2,i) =
     $        zwk1(jmax2+jwal1+1:jmax2+jwal2,i)
     $        / sqrt(rwvlo(1:jwal1))
c...  Extra up-down terms..
         zwk4(nspw2+1:nspw2+jwal1,i) = zwk1(nspw2+1:nspw2+jwal1,i)
     $        * sqrt(rwvle(1:jwal1))
         zwk4(nspw2+jwal1+1:nspw2+jwal2,i) =
     $        zwk1(nspw2+jwal1+1:nspw2+jwal2,i)
     $        * sqrt(rwvlo(1:jwal1))
         zwk1(nspw2+1:nspw2+jwal1,i) = zwk1(nspw2+1:nspw2+jwal1,i)
     $        / sqrt(rwvle(1:jwal1))
         zwk1(nspw2+jwal1+1:nspw2+jwal2,i) =
     $        zwk1(nspw2+jwal1+1:nspw2+jwal2,i)
     $        / sqrt(rwvlo(1:jwal1))

      END DO

      izcnt = 0
c ===================  SENSOR FLUXES FROM THE COILS  ================

c... Need the sensor loops in the next do loop. Put it here outside the loop.
c      And also calculate the fluxes in the sensor from the feedback
c      coils here.

         CALL sloops ( slx, slxt, slz, slzt, nslpt )

c..   Calculate the flux through the sensor loops due to the feedback coil
c     interacting with the plasma and the wall. This is independent OPLP

c.. Write out the coordinates of the coil tips:

      WRITE ( 28, '(20x, "X", 12x, "Z", 12x, "X", 12x, "Z" )' )
      WRITE ( 28, '(1x, "C-COIL tips: ", 4ES13.5 )' )
     $     (chgt(i) , i = 1,4)
      WRITE ( 28, '(1x, "L-COIL tips: ", 4ES13.5 )' )
     $     (chgta(i) , i = 1,4)
      WRITE ( 28, '(1x, "U-COIL tips: ", 4ES13.5 )' )
     $     (chgtb(i) , i = 1,4)

      WRITE ( 28, '( /, 1x,"---------------------------------" )' ) 
      WRITE ( 28, '( 1x,".....SENSOR.... nslpt = ", i5 )' ) nslpt

      zphil = 0.0           ! This is the \phi location of the loop.
      lpols = lpolsensor
      CALL sensflxi ( slx, slxt, slz, slzt, nslpt, zphil, lpols,
     $     "vert" )

      WRITE ( 6, '( "XXXXX  RESSHEL-A nslpt = ", i5 )' ) nslpt
      WRITE ( 23, '( "XXXXX RESSHEL-A nslpt = ", i5 )' ) nslpt
      WRITE ( 28, '( "XXXXX RESSHEL-A nslpt = ", i5 )' ) nslpt

      IF ( ld3dflx == 1 ) THEN
c.. The fluxes in DIII-D radial loops:

         WRITE ( 28, '( 20x,
     $        "#################################################" )' ) 
      WRITE ( outmod, '(/, 10x, "Entering DIII-D loop measurements" )' )
      WRITE ( iotty,  '(/, 10x, "Entering DIII-D loop measurements" )' )
      WRITE ( 28,  '(/, 10x, "Entering DIII-D loop measurements" )' )

c... Comment this out for checking with res3
      CALL diii_loop_meas (nzpts, dtsl, lpols )

      WRITE ( 6, '( "XXXXX  RESSHEL-B nzspt = ", i5 )' ) nzpts
      WRITE ( 23, '( "XXXXX RESSHEL-B nzspt = ", i5 )' ) nzpts
      WRITE ( 28, '( "XXXXX RESSHEL-B nzspt = ", i5 )' ) nzpts
      END IF

      DO i = 1, lrswcomp
         
         WRITE ( outmod, '( 20x,
     $        "#################################################" )' ) 
         WRITE ( iotty, '( 20x,
     $        "#################################################" )' ) 
         
         WRITE ( outmod,
     $        '( 25x, "Eigenmode of OPLP No. ", I4, ": " )' ) i
         WRITE ( iotty,
     $        '( 25x, "Eigenmode of OPLP No. ", I4, ": " )' ) i
         
         WRITE ( outmod, '( 20x,
     $        "-------------------------------------------------",/ )' ) 
         WRITE ( iotty, '( 20x,
     $        "-------------------------------------------------",/ )' ) 

         izcnt = izcnt + 1
c     izm = nspw2 - i + 1
c     izm = jmax1 + 4 + i
         izm = i
c     
c     kwze and kwzo are the the value of current potentials
c     in along the wall. Used for plotting Skin current.
c     
         kwze(1:mth1) =
     $        MATMUL ( rwvce(1:mth1,1:jwal1),
     $        zwk1(jmax2+1:jmax2+jwal1,izm) )
         kwze(1:mth1) = kwze(1:mth1) + ! Odd-Cosine term
     $        MATMUL ( rwvco(1:mth1,1:jwal1),
     $        zwk1(jmax2+jwal1+1:jmax2+jwal2,izm) )
         

c     xibpe and xibpo are the amplitudes of magnetic perturbation
c     on the plasma. Divide by the resonance operator to get the
c     perturbations xize and xizo for the even and odd displacements

         xibpe(1:jmax1)=zwk1(1:jmax1,izm)

         DO ij=1,jmax1
            xize(ij)=xibpe(ij)/(lmin(1)+ij-1-n*qa1)
         END DO

         iksgn=0

c     if (kwze(1).lt.0.) then
c     iksgn=1
c     kwze(1:mth1)=-kwze(1:mth1)
c     end if
c     the current potential should be proportional to the growth
c     rate

         kwze(1:mth1) = zwk2(izm)*kwze(1:mth1)
c     CALL pospl1 ( pigrd, kwze, 1,mth1, izcnt,"kwz-e", i )
c     
c     bwze and bwzo are the actual values of the magnetice perturbation
c     on the wall in real space
c     xibwe and xibwo are the amplitudes wrt to the input wall
c     eigenfunctions
c     
c$$$  bwze(1:mth1) =
c$$$  $        MATMUL ( rwvce(1:mth1,1:jwal1),
c$$$  $        zwk4(jmax2+1:jmax2+jwal1,izm) )*2.0
         bwze(1:mth1) =
     $        MATMUL ( rwvce(1:mth1,1:jwal1),
     $        zwk4(jmax2+1:jmax2+jwal1,izm) )*2.0
         bwze(1:mth1) = bwze(1:mth1) + ! Odd-Cosine terms
     $        MATMUL ( rwvco(1:mth1,1:jwal1),
     $        zwk4(jmax2+jwal1+1:jmax2+jwal2,izm) )*2.0

         xibwe(1:jwal2)=zwk4(jmax2+1:jmax2+jwal2,izm)*2.0
         ibsgn=0
c     if (bwze(1).lt.0.) then
c     ibsgn=1
c     bwze(1:mth1)=-bwze(1:mth1)
c     end if
c     
c     we dont need to divide bwze by major radius because
c     it has already been taken care of in the jacobian
c     
c     bwze(1:mth1) = bwze(1:mth1)/xwal(1:mth1)

c...  bpzc and bpzs are the actual values of the magnetic perturbation 
c     on the plasmsa, cosine and sine terms, in theta space.

         bpzc(1:mth1) = MATMUL ( cslth(1:mth1,1:jmax1),xibpe(1:jmax1) )
         bpzs(1:mth1) = MATMUL ( snlth(1:mth1,1:jmax1),xibpe(1:jmax1) )

         excw  = dth * DOT_PRODUCT ( clrsp(1:mth),bwze(1:mth) )

         CALL shftpi (pigrd , bwze, z1tmp, z2tmp , mth1)
         CALL pospl1 ( z1tmp, z2tmp, 1,mth1, izcnt,"bwz-e", i )

         exe(i)=excw

         WRITE ( outmod, '( "izm, xibpe = " , 2i5, es13.5 )' )
     $        (izm, ij, xibpe(ij), ij= 1, jmax1)
         WRITE ( iotty,  '( "izm, xibpe = " , 2i5, es13.5 )' )
     $        (izm, ij, xibpe(ij), ij= 1, jmax1)

c...  Even Wall Excitations for the internal coils:

         excw  = dth * DOT_PRODUCT ( skcoilwa(1:mth),bwze(1:mth) )
         exawe(i) = excw

         excw  = dth * DOT_PRODUCT ( skcoilwb(1:mth),bwze(1:mth) )
         exbwe(i) = excw

         CALL pospl1 ( xlfm, xibpe, 1,jmax1, izcnt+2,"bp_l-e", i )
c     WRITE ( outmod, '( "xize = " , i5, es13.5 )' )
c     $        (ij, xize(ij), ij= 1, jmax1)
c     WRITE ( iotty,  '( "xize = " , i5, es13.5 )' )
c     $        (ij, xize(ij), ij= 1, jmax1)
         CALL pospl1 ( xlfm, xize, 1,jmax1, izcnt+3,"xi_l-e", i )
         izcnt = izcnt + 1
c     
c     the odd part of the eigenfunction, for the moment,
c     this odd part comes from the representation of complex
c     exponential function for the magnetic perturbation on the plasma
c     surface but real functions on the wall
c     
c...  Kwzo used for plotting skin current
         kwzo(1:mth1) =
     $        MATMUL ( rwvco(1:mth1,1:jwal1),
     $        zwk1(jmax2+jwal3+1:jmax2+jwal4,izm) )
         kwzo(1:mth1) =  kwzo(1:mth1) + ! Even-Sine terms
     $        MATMUL ( rwvce(1:mth1,1:jwal1),
     $        zwk1(jmax2+jwal2+1:jmax2+jwal3,izm) )

         if (iksgn.eq.1)kwzo(1:mth1)=-kwzo(1:mth1)
         kwzo(1:mth1) = zwk2(izm)* kwzo(1:mth1)
c     CALL pospl1 ( pigrd, kwzo, 1,mth1, izcnt,"kwz-o", i )

         bwzo(1:mth1) =
     $        MATMUL ( rwvco(1:mth1,1:jwal1),
     $        zwk4(jmax2+jwal3+1:jmax2+jwal4,izm) )*2.0
         bwzo(1:mth1) = bwzo(1:mth1) + !  Even-Sine terms.
     $        MATMUL ( rwvce(1:mth1,1:jwal1),
     $        zwk4(jmax2+jwal2+1:jmax2+jwal3,izm) )*2.0

         IF ( i == 1 ) THEN

            ALLOCATE ( chiwallr(mth2), chiwalli(mth2) )
            ALLOCATE (b_thet_w_r(mth2), b_thet_w_i(mth2),
     $           b_phi_w_r(mth2), b_phi_w_i(mth2) )

            iobwal = 33
            OPEN (  iobwal, FILE='B-SHELL', STATUS='REPLACE',
     $           FORM='FORMATTED' )
            WRITE ( iobwal, '("No. of points = ", i5 )' ) mth
            WRITE ( iobwal, '("REAL-B_normal = ", / (6ES13.5) )' )
     $           ( bwze(iiz), iiz = 1, mth )
            WRITE ( iobwal, '("IMAG-B_normal = ", / (6ES13.5) )' )
     $           ( bwzo(iiz), iiz = 1, mth )
c.. Now, calculate chi(wall) from the response matrices and the OPLP
c        eigenfunction.
c.. Plasma sources:
            chiwallr(1:mth1) =
     $           MATMUL(cwallr(1:mth1,1:jmax1),zwk1(1:jmax1,i)) -
     $           MATMUL(cwalli(1:mth1,1:jmax1),zwk1(jmax1+1:jmax2,i))
            chiwalli(1:mth1) =
     $           MATMUL(cwallr(1:mth1,1:jmax1),zwk1(jmax1+1:jmax2,i)) +
     $           MATMUL(cwalli(1:mth1,1:jmax1),zwk1(1:jmax1,i))
c.. Wall sources:
            chiwallr(1:mth1) = chiwallr(1:mth1) + 
     $           MATMUL(cwwale(1:mth1,1:jwal1),
     $           zwk1(jmax2+1:jmax2+jwal1,i)) +
     $           MATMUL(cwwalo(1:mth1,1:jwal1),
     $           zwk1(jmax2+jwal1+1:jmax2+jwal2,i))
            chiwalli(1:mth1) = chiwalli(1:mth1) + 
     $           MATMUL(cwwale(1:mth1,1:jwal1),
     $           zwk1(jmax2+jwal2+1:jmax2+jwal3,i)) +
     $           MATMUL(cwwalo(1:mth1,1:jwal1),
     $           zwk1(jmax2+jwal3+1:jmax2+jwal4,i))
c.. Calculate B_theta on shell:
            call difspl ( mth, pigrd, chiwallr, b_thet_w_r )
            call difspl ( mth, pigrd, chiwalli, b_thet_w_i )
            b_thet_w_r(mth1) = b_thet_w_r(1)
            b_thet_w_i(mth1) = b_thet_w_i(1)
            b_thet_w_r (1:mth1) = xwal(1:mth1) * b_thet_w_r(1:mth1) /
     $           gpsjw(1:mth1)
            b_thet_w_i (1:mth1) = xwal(1:mth1) * b_thet_w_i(1:mth1) /
     $           gpsjw(1:mth1)
c.. Calculate B_phi on shell:
            b_phi_w_r(1:mth1) =   n * chiwalli(1:mth1) / xwal(1:mth1)
            b_phi_w_i(1:mth1) = - n * chiwallr(1:mth1) / xwal(1:mth1)
c.. Write to disk:
            WRITE ( iobwal, '("REAL-chi = ", / (6ES13.5) )' )
     $           ( chiwallr(iiz), iiz = 1, mth )
            WRITE ( iobwal, '("IMAG-chi = ", / (6ES13.5) )' )
     $           ( chiwalli(iiz), iiz = 1, mth )
            WRITE ( iobwal, '("REAL-B_theta = ", / (6ES13.5) )' )
     $           ( b_thet_w_r(iiz), iiz = 1, mth )
            WRITE ( iobwal, '("IMAG-B_theta = ", / (6ES13.5) )' )
     $           ( b_thet_w_i(iiz), iiz = 1, mth )
            WRITE ( iobwal, '("REAL-B_phi = ", / (6ES13.5) )' )
     $           ( b_phi_w_r(iiz), iiz = 1, mth )
            WRITE ( iobwal, '("IMAG-B_phi = ", / (6ES13.5) )' )
     $           ( b_phi_w_i(iiz), iiz = 1, mth )

c....Plots:

            CALL framep(jobid, ff)

            CALL shftpi ( pigrd, chiwallr, zork1, zork2, mth1 )
            CALL shftpi ( pigrd, chiwalli, zork1, zork3, mth1 )
            CALL pospl2 ( zork1, zork2, zork3, 1, mth1, 1,
     $           "chiwall", mth )
            CALL shftpi ( pigrd, bwze, zork1, zork2, mth1 )
            CALL shftpi ( pigrd, bwzo, zork1, zork3, mth1 )
            CALL pospl2 ( zork1, zork2, zork3, 1, mth1, 2,
     $           "B-wall_n", mth )
            CALL shftpi ( pigrd, b_thet_w_r, zork1, zork2, mth1 )
            CALL shftpi ( pigrd, b_thet_w_i, zork1, zork3, mth1 )
            CALL pospl2 ( zork1, zork2, zork3, 1, mth1, 3,
     $           "B-wall_theta", mth )
            CALL shftpi ( pigrd, b_phi_w_r, zork1, zork2, mth1 )
            CALL shftpi ( pigrd, b_phi_w_i, zork1, zork3, mth1 )
            CALL pospl2 ( zork1, zork2, zork3, 1, mth1, 4,
     $           "B-wall_phi", mth )
            
            CALL framep(jobid, ff)
            
            CLOSE ( UNIT = 33 )
            DEALLOCATE ( chiwallr, chiwalli, b_thet_w_r, b_thet_w_i,
     $           b_phi_w_r, b_phi_w_i )

      END IF

         xibwo(1:jwal2)=zwk4(jmax2+jwal2+1:jmax2+jwal4,izm)*2.0

         xibpo(1:jmax1)=zwk1(jmax1+1:jmax2,izm)

         DO ij=1,jmax1
            xizo(ij)=xibpo(ij)/(lmin(1)+ij-1-n*qa1)
         END DO

         WRITE ( outmod, '( "izm, xibpo = " , 2i5, es13.5 )' )
     $        (izm, ij, xibpo(ij), ij= 1, jmax1)
         WRITE ( iotty,  '( "izm, xibpo = " , 2i5, es13.5 )' )
     $        (izm, ij, xibpo(ij), ij= 1, jmax1)

c..   Now the imaginary contribs. from up-down asymmetry
c...  bpzoc and bpzos are the actual values of the magnetic perturbation 
c     on the plasmsa, cosine and sine terms, in theta space.

         bpzoc(1:mth1) = MATMUL ( cslth(1:mth1,1:jmax1),xibpo(1:jmax1) )
         bpzos(1:mth1) = MATMUL ( snlth(1:mth1,1:jmax1),xibpo(1:jmax1) )


         if (ibsgn.eq.1)bwzo(1:mth1)=-bwzo(1:mth1)
c     do not need to divide by xwal
c     bwzo(1:mth1) = bwzo(1:mth1)/xwal(1:mth1)
c     
         call shftpi (pigrd , bwzo, z1tmp , z2tmp , mth1)
         CALL pospl1 ( z1tmp, z2tmp, 1,mth1, izcnt,"bwz-o", i )
         excw  = dth * DOT_PRODUCT ( clrsp(1:mth),bwzo(1:mth) )
         exo(i)=excw

         WRITE ( outmod, '()' )
         WRITE ( outmod, '( "excw-e = " , i5, es13.5 )' ) i, exe(i)
         WRITE ( iotty,  '( "excw-e = " , i5, es13.5 )' ) i, exe(i)

         WRITE ( outmod, '( "excw-o = " , i5, es13.5 )' ) i, exo(i)
         WRITE ( iotty,  '( "excw-o = " , i5, es13.5 )' ) i, exo(i)

         examp(i)=sqrt(exe(i)**2+exo(i)**2)
         exang(i)=atan2(exo(i),exe(i))

c...  Odd  Wall Excitations for the internal coils:

         excw  = dth * DOT_PRODUCT ( skcoilwa(1:mth),bwzo(1:mth) )
         exawo(i) = excw

         excw  = dth * DOT_PRODUCT ( skcoilwb(1:mth),bwzo(1:mth) )
         exbwo(i) = excw

         WRITE ( outmod, '()' )
         WRITE ( outmod, '( "excw-e-wa = " , i5, es13.5 )' ) i, exawe(i)
         WRITE ( iotty,  '( "excw-e-wa = " , i5, es13.5 )' ) i, exawe(i)
         WRITE ( outmod, '( "excw-e-wb = " , i5, es13.5 )' ) i, exbwe(i)
         WRITE ( iotty,  '( "excw-e-wb = " , i5, es13.5 )' ) i, exbwe(i)

         WRITE ( outmod, '( "excw-o-wa = " , i5, es13.5 )' ) i, exawo(i)
         WRITE ( iotty,  '( "excw-o-wa = " , i5, es13.5 )' ) i, exawo(i)
         WRITE ( outmod, '( "excw-o-wb = " , i5, es13.5 )' ) i, exbwo(i)
         WRITE ( iotty,  '( "excw-o-wb = " , i5, es13.5 )' ) i, exbwo(i)

c..   Amplitudes and phase, wall:

         exawamp(i)=sqrt(exawe(i)**2+exawo(i)**2)
         exawang(i)=atan2(exawo(i),exawe(i))
         exbwamp(i)=sqrt(exbwe(i)**2+exbwo(i)**2)
         exbwang(i)=atan2(exbwo(i),exbwe(i))

c...  Plasma Excitations for the internal a nd b coils, both cosine and sine.
c...  Add in the contribs. from up-down asymmetry.

         exapc(i) = dth * DOT_PRODUCT ( skcoilpa(1:mth),
     $        ( bpzc(1:mth)-bpzos(1:mth) ) )
         exbpc(i) = dth * DOT_PRODUCT ( skcoilpb(1:mth), 
     $        ( bpzc(1:mth)-bpzos(1:mth) ) )

         exaps(i) = dth * DOT_PRODUCT ( skcoilpa(1:mth), 
     $        ( bpzs(1:mth)+bpzoc(1:mth) ) )
         exbps(i) = dth * DOT_PRODUCT ( skcoilpb(1:mth), 
     $        ( bpzs(1:mth)+bpzoc(1:mth) ) )

         WRITE ( outmod, '( )' )
         WRITE ( outmod, '( "excp-c-pa = " , i5, es13.5 )' ) i, exapc(i)
         WRITE ( iotty,  '( "excp-c-pa = " , i5, es13.5 )' ) i, exapc(i)
         WRITE ( outmod, '( "excp-c-pb = " , i5, es13.5 )' ) i, exbpc(i)
         WRITE ( iotty,  '( "excp-c-pb = " , i5, es13.5 )' ) i, exbpc(i)

         WRITE ( outmod, '( "excp-s-pa = " , i5, es13.5 )' ) i, exaps(i)
         WRITE ( iotty,  '( "excp-s-pa = " , i5, es13.5 )' ) i, exaps(i)
         WRITE ( outmod, '( "excp-s-pb = " , i5, es13.5 )' ) i, exbps(i)
         WRITE ( iotty,  '( "excp-s-pb = " , i5, es13.5 )' ) i, exbps(i)
         WRITE ( outmod, '()' )

c..   Amplitudes and phase, plasma:

         exapamp(i)=sqrt(exapc(i)**2+exaps(i)**2)
         exapang(i)=atan2(exaps(i),exapc(i))
         exbpamp(i)=sqrt(exbpc(i)**2+exbps(i)**2)
         exbpang(i)=atan2(exbps(i),exbpc(i))

         CALL framep(jobid, ff)

c...  Plot the plasma magnetic field as a function of theta.

c..   Combine BPZC = BPZC - BPZOS and BPZS = BPZS + BPOC

         bpzc = bpzc - bpzos
         bpzs = bpzs + bpzoc

         call shftpi (pigrd , bpzc, z1tmp, z2tmp , mth1)
         CALL pospl1 ( z1tmp, z2tmp, 1,mth1, 1,"bpz-c", i )

         call shftpi (pigrd , bpzs, z1tmp, z2tmp , mth1)
         CALL pospl1 ( z1tmp, z2tmp, 1,mth1, 2,"bpz-s", i )

c..   Plot Wall Magnetic field as a function of jwal.

         CALL pospl1 ( xlwm, xibwe, 1,jwal2, 3,"bw-e_j", i )
         CALL pospl1 ( xlwm, xibwo, 1,jwal2, 4,"bw-o_j", i )

         CALL framep(jobid, ff)

c =============================================
c..   The sensor fluxes from the modes:
c =============================================

c$$$         ALLOCATE ( slx(nslpt), slxt(nslpt), slz(nslpt), slzt(nslpt) )

         lpols = lpolsensor
         CALL sensorp(slx, slxt, slz, slzt, nslpt, lpols,
     $        xibpe,xibpo, rflbpe, rflbpo, pflbpe, pflbpo)
         CALL sensorw(slx, slxt, slz, slzt, nslpt, lpols,
     $        xibwe,xibwo, rflbwe, rflbwo, pflbwe, pflbwo)
         CALL chiplasma(xibpe,xibpo,xibwe,xibwo,chiplasr,
     $        chiplasi)
      IF ( lpols == 0 ) THEN
         flbpe(izm) = rflbpe
         flbpo(izm) = rflbpo
         flbwe(izm) = rflbwe
         flbwo(izm) = rflbwo
      ELSE IF ( lpols ==1 ) THEN
         flbpe(izm) = pflbpe
         flbpo(izm) = pflbpo
         flbwe(izm) = pflbwe
         flbwo(izm) = pflbwo
      END IF

c....   Save the radial sensor flux for the first mode:
      IF ( izm ==1 ) THEN
         fltamp1 = SQRT ( (rflbpe+rflbwe)**2 + (rflbpo+rflbwo)**2 )
         WRITE ( outmod, '(/,5x,"RADIAL SENSOR FLUX From First Mode = ",
     $        es11.4, / )' ) fltamp1
         WRITE ( 28, '(/,5x,"RADIAL SENSOR FLUX From First Mode = ",
     $        es11.4, /  )' ) fltamp1
      END IF
         fltte(izm)=flbpe(izm)+flbwe(izm)
         fltto(izm)=flbpo(izm)+flbwo(izm)
         fltamp(i)=sqrt(fltte(i)**2+fltto(i)**2)
         fltang(i)=atan2(fltto(i),fltte(i))
         WRITE ( outmod, '( "flpwte-o=",i2,6es11.3 )' ) 
     &        i, flbpe(i),flbwe(i),fltte(i),flbpo(i),flbwo(i),fltto(i)
         WRITE ( iotty,  '( "flpwte-o=",i2,6es11.3 )' )
     &        i, flbpe(i),flbwe(i),fltte(i),flbpo(i),flbwo(i),fltto(i)
c..   Plot displacement on the plasma surface
         xili(1:jmax1)=-xize(1:jmax1)
         xilr(1:jmax1)= xizo(1:jmax1)
         call makebp ( "RWM Eigenfunction" )

c..   Plot skin current on shell.
c     
c     First, shift variables so that theta=0 is on the inside.
c     
         workl(1:mth1)  = (/ xwal(mth/2+1:mth),xwal(1:mth/2+1) /)
         workl2(1:mth1) = (/ zwal(mth/2+1:mth),zwal(1:mth/2+1) /)
c     
         nbngr = 129

         kwze(1:mth1) = (/ kwze(mth/2+1:mth),kwze(1:mth/2+1) /)
         kwzo(1:mth1) = (/ kwzo(mth/2+1:mth),kwzo(1:mth/2+1) /)

c     iwrto = outmod
         iwrto = 0

         nbngr = 129
         call kcur0 ( workl,workl2,ell,wgrd,mth1,kwze,kwzo,n,
     $        nbngr,wkgrd,wxgrd,wzgrd,wkdxt,wkdzt,
     $        wkt2,wkp2,nths,0,"Shell Eig-Cur",i, iwrto )

c...  Plot the relative phases of the mode pattern on the shell at the
c     coil centers.

         iwrto1 = outmod
         iwrto2 = iotty
         call mphase ( workl,workl2,ell,wgrd,mth1,kwze,kwzo,n,
     $        nbngr,coilcen,coilcena,coilcenb,
     $        wkgrd,wxgrd,wzgrd,wkdxt,wkdzt,
     $        wkt2,wkp2,nths,0,"Eig-Phase",i, iwrto1, iwrto2 )

         IF ( izcnt == 2 ) THEN
c     CALL framep(jobid, ff)
            izcnt = 0
         END IF

c...  Write out the OPLP perturbations in real space:
         
         WRITE ( ioop, '(30x, "OPLP Eigenfunction(theta) No. ", i4 )' )
     $        izm
         WRITE ( ioop, '(3x, "I", 4x, "B-plas-Cos", 5x, "B-plas-Sin",
     $        5x, "B-wall-Cos",5x,"B-wall-Sin")' ) 
         
         DO iw = 1, mth1
            WRITE ( ioop, '(i4, 4ES15.7)' ) iw, bpzc(iw), bpzs(iw),
     $           bwze(iw), bwzo(iw)
         END DO

c..   Plot virtual skin current on plasma using OPLP's B_n as source

c     First, shift variables so that theta=0 is on the inside.

         workl(1:mth1)  = (/ xinf(mth/2+1:mth),xinf(1:mth/2+1) /)
         workl2(1:mth1) = (/ zinf(mth/2+1:mth),zinf(1:mth/2+1) /)
c     
         nbngr = 129

c...  Get b_normal only, i.e., divide by j*grpssq 
         bpzc(1:mth1) = bpzc(1:mth1) !  / gpsjp(1:mth1)
         bpzs(1:mth1) = bpzs(1:mth1) !  / gpsjp(1:mth1)

         bpzc(1:mth1) = (/ bpzc(mth/2+1:mth),bpzc(1:mth/2+1) /)
         bpzs(1:mth1) = (/ bpzs(mth/2+1:mth),bpzs(1:mth/2+1) /)

c     iwrto = outmod
         iwrto = 0

         nbngr = 129
         call kcur0 ( workl,workl2,ellp,wgrd,mth1,bpzc,bpzs,n,
     $        nbngr,wkgrd,wxgrd,wzgrd,wkdxt,wkdzt,
     $        wkt2,wkp2,nths,0,"Fake OPLP K_p",i, iwrto )

      END DO                    ! End Loop i = lrswcomp.

      IF ( lreshel /= 0 ) CLOSE ( UNIT = iosen )

      CLOSE ( UNIT = ioop )

      WRITE ( outmod, '( 20x,
     $     "#################################################" )' ) 
      WRITE ( iotty, '( 20x,
     $     "#################################################" )' ) 

      WRITE ( outmod, '( 30x, "End of OPLP Eigenmode Loop" )' )
      WRITE ( iotty,  '( 30x, "End of OPLP Eigenmode Loop" )' )

      WRITE ( outmod, '( 20x,
     $     "-------------------------------------------------" )' ) 
      WRITE ( iotty, '( 20x,
     $     "-------------------------------------------------" )' ) 


c$$$  DO i = 5, 8
c$$$  izm = nspw2 - i + 1
c$$$  bwz(1:mth) =
c$$$  $        MATMUL ( rwvce(1:mth,1:jwal1),
c$$$  $        zwk1(jmax1+1:jmax1+jwal1,izm) )
c$$$  bwz(1:mth) = bwz(1:mth) +
c$$$  $        MATMUL ( rwvco(1:mth,1:jwal1),
c$$$  $        zwk1(jmax1+jwal1+1:jmax1+jwal2,izm) )
c$$$  
c$$$  excw  = dth * DOT_PRODUCT ( clrsp(1:mth),bwz(1:mth) )
c$$$  
c$$$  CALL pospl1 ( pigrd, bwz, 1,mth, i-4,"bwz", i )
c$$$  
c$$$  WRITE ( outmod, '( "excw = " , i5, es13.5 )' ) i, excw
c$$$  WRITE ( iotty,  '( "excw = " , i5, es13.5 )' ) i, excw
c$$$  END DO
c$$$  CALL framep(jobid, ff)

c..   Write out excitations.

c..   External coil:

      WRITE ( outmod, '()' )

      DO i=1,lrswcomp,1
         WRITE( outmod, '( "exe,exo,examp,exang",i5,4es13.5)')
     &        i,exe(i),exo(i),examp(i),exang(i)
      END DO

c..   Internal Coils:

      WRITE ( outmod, '()' )
      DO i=1,lrswcomp,1
         WRITE( outmod, '( "exawe,exawo,exawamp,exawang",i5,4es13.5)')
     &        i,exawe(i),exawo(i),exawamp(i),exawang(i)
      END DO
      
      WRITE ( outmod, '()' )
      DO i=1,lrswcomp,1
         WRITE( outmod, '( "exbwe,exbwo,exbwamp,exbwang",i5,4es13.5)')
     &        i,exbwe(i),exbwo(i),exbwamp(i),exbwang(i)
      END DO

      WRITE ( outmod, '()' )
      DO i=1,lrswcomp,1
         WRITE( outmod, '( "exapc,exaps,exapamp,exapang",i5,4es13.5)')
     &        i,exapc(i),exaps(i),exapamp(i),exapang(i)
      END DO
      
      WRITE ( outmod, '()' )
      DO i=1,lrswcomp,1
         WRITE( outmod, '( "exbpc,exbps,exbpamp,exbpang",i5,4es13.5)')
     &        i,exbpc(i),exbps(i),exbpamp(i),exbpang(i)
      END DO

c..   Sensor loop stuff:
      WRITE ( outmod, '()' )
      DO i=1,lrswcomp,1
         WRITE( outmod, '( "fle,flo,flamp,flang",i5,4es13.5)')
     &        i,fltte(i),fltto(i),fltamp(i),fltang(i)
      END DO

c$$$      DEALLOCATE ( slx, slxt, slz, slzt )
      DEALLOCATE (zwk1, zwk2, zwk3, zwk4, bwze, bwzo ,kwze,kwzo )
      DEALLOCATE ( bpzc, bpzs,  bpzoc, bpzos )
      DEALLOCATE (zbwk, xize, xizo, xibpe ,z1tmp, z2tmp, xibpo, xibwe )
      DEALLOCATE (xibwo, flbpe, flbpo, flbwe, flbwo )
      DEALLOCATE (exe, exo, examp, exang )
      DEALLOCATE ( exawo, exawe, exawamp, exawang )
      DEALLOCATE ( exbwo, exbwe, exbwamp, exbwang )
      DEALLOCATE ( exaps, exapc, exapamp, exapang )
      DEALLOCATE ( exbps, exbpc, exbpamp, exbpang )
      DEALLOCATE ( dwp, dwe, dwo, dwpe, dwpo )

      END SUBROUTINE resshel_2

c.......................................................
      function cth(t)
c.......................................................
c     
      a = 0.0
      b = 1.0
      tau = 0.02
      f = a + b*cos(t)
      exf = exp(f/tau)
      ff = exf + 1.0
      cth = 1.0 - 1.0 / ff
c     
      entry cph(t)
c     
      fp = - b*sin(t)
      cph = fp * exf / (tau * ff**2)
c     
      return
      end
c
c................................................................
      subroutine venergy ( zlabel )
c................................................................
c
c..   Calculates the vacuum energy using the surface perturbation
c     from PEST or GATO
c
      include 'vacuum1.inc'
      include 'vacuum2.inc'
      include 'vacuum3.inc'
      include 'vacuum8.inc'
      character*(*) zlabel
c
      venerg = 0.0
c
      if ( lgato .gt. 0 ) then
         work(1:mfel)  = matmul ( gatovac(1:mfel,1:mfel),bnkr(1:mfel) )
         work1(1:mfel) = matmul ( gatovac(1:mfel,1:mfel),bnki(1:mfel) )
c
         venerg = dot_product ( bnkr(1:mfel),work(1:mfel) )
     $        + dot_product ( bnki(1:mfel),work1(1:mfel) )

         venerg = pye * venerg

         print "(/,2A,1pe13.5)", zlabel, ": Vacuum Energy = ", venerg
         write (outmod, "(/,2A,1pe13.5)")
     $        zlabel, ": Vacuum Energy = ", venerg
      end if
c
      if ( lgato .eq. 0 ) then
         work(1:lrnge)  =
     $        matmul ( vacmat0(1:lrnge,1:lrnge),bnlr(1:lrnge) )
         work1(1:lrnge) =
     $        matmul ( vacmat0(1:lrnge,1:lrnge),bnli(1:lrnge) )
c
         venerg = dot_product ( bnlr(1:lrnge),work(1:lrnge) )
     $        + dot_product ( bnli(1:lrnge),work1(1:lrnge) )

         venerg = pye * venerg
         if ( lrgato .eq. 2 ) then
c            zfac =  rgato / (twopi*mfel)
c..   futher normalized with PI/(2*\mu_0 * R)
            zfac = 10**7 / (8.0*rgato)
            venerg = zfac * venerg
         end if

         print "(/,2A,1pe13.5)", zlabel,
     $        ": Vacuum Energy = ", venerg
         write (outmod, "(/,2A,1pe13.5)")
     $        zlabel, ": Vacuum Energy = ", venerg
      end if
c
      return
      end
c
c................................................................
      subroutine cplt ( blr, bli )
c................................................................
c
      include 'vacuum1.inc'
      include 'vacuum2.inc'
      include 'vacuum3.inc'
      include 'vacuum5.inc'
      include 'vacuum6.inc'
      include 'vacuum7.inc'
      include 'vacuum8.inc'
c
      dimension blr(*), bli(*)
      dimension wkxx(nths), wkyy(nths)
c
      jmax1 = lrnge
c
c.. Calculate B-normal for each diagonalized wall function on the shell
c    driven by the plasma response.

c     Sum over l with blr and bli  here.
c
      wkxx(1:jwal2) = matmul ( blsp(1:jwal2,1:jmax1),blr(1:jmax1) ) -
     $     matmul ( blsp(1:jwal2,jmax1+1:2*jmax1),bli(1:jmax1) )
      wkyy(1:jwal2) = matmul ( blsp(1:jwal2,1:jmax1),bli(1:jmax1) ) +
     $     matmul ( blsp(1:jwal2,jmax1+1:2*jmax1),blr(1:jmax1) )

      WRITE ( OUTMOD,
     $     '(/,5x, "Wall Response to Plasma, and Time Constants",/,
     $     2x,"i",5x,"Real", 10x, "Imag", 9x, "T-msec")' )

      WRITE ( OUTMOD, '(i4, 1p3e13.5)' )
     $     (i, wkxx(i), wkyy(i), tc_vess(i), i = 1, jwal2 )

c.. Calculate quantities for plotting the time asymptotic shell current:

c     Sum over l with blr and bli  here.
c
      wkxx(1:mw) = matmul ( gwin(1:mw,1:jmax1),blr(1:jmax1) ) -
     $     matmul ( gwot(1:mw,1:jmax1),bli(1:jmax1) )
      wkyy(1:mw) = matmul ( gwot(1:mw,1:jmax1),blr(1:jmax1) ) +
     $     matmul ( gwin(1:mw,1:jmax1),bli(1:jmax1) )
c
      write ( outmod, '(/,5x, "Long time Response Current",/,
     $     2x, "i",5x, "Real", 10x, "Imag")' )
c
      write ( outmod, '(i4, 1p2e13.5)' )
     $     (i, wkxx(i), wkyy(i), i = 1, mw, 5 )
c
      write ( outmod, '(/, 5x, "Magnetic Source Components",/,
     $     2x, "i",5x, "Real", 10x, "Imag")' )
c
      write ( outmod, '(i4, 1p2e13.5)' )
     $     (i, blr(i), bli(i), i = 1, jmax1 )
c
c..   Plot the First Few Eigencurrents
c     Assume Imaginary part of current coeff. is zero:
c     .  Use wkxx, wkyy as storage.
c
c     First, shift variables so that theta=0 is on the inside.
c     Use workl, workl2 as storage.
c
      workl(1:mth1)  = (/ xwal(mth/2+1:mth),xwal(1:mth/2+1) /)
      workl2(1:mth1) = (/ zwal(mth/2+1:mth),zwal(1:mth/2+1) /)
c
      nbngr = 129

      wkxx(1:mth1) = (/ wkxx(mth/2+1:mth),wkxx(1:mth/2+1) /)
      wkyy(1:mth1) = (/ wkyy(mth/2+1:mth),wkyy(1:mth/2+1) /)
c
c     iwrto = outmod
      iwrto = 0
      call kcur0 ( workl,workl2,ell,wgrd,mth1,wkxx,wkyy,n,
     $     nbngr,wkgrd,wxgrd,wzgrd,wkdxt,wkdzt,
     $     wkt2,wkp2,nths,0,"Asym-Current",lrnge, iwrto )
c
      return
      end

c$$$c.................................................................
c$$$      subroutine mateig ( zvec, work, work1 )
c$$$c.................................................................
c$$$c
c$$$      include 'vacuum1.inc'
c$$$c
c$$$      dimension zvec(nfm,nfm), work(*), work1(*)
c$$$c
c$$$c     find eigenvalues of vecin....
c$$$c
c$$$      jmax1 = lmax(1) - lmin(1) + 1
c$$$c
c$$$c$$$      do 9195 i = 1, nfmsq
c$$$c$$$         work1(i) = zero
c$$$c$$$ 9195 continue
c$$$c
c$$$      work1(1:jmax1*(jmax1+1)/2) = (/ ((zvec(i,j),j=1,i),i=1,jmax1) /)
c$$$c
c$$$c$$$      kk = 1
c$$$c$$$      do 9201 i = 1, jmax1
c$$$c$$$         do 9200 j = 1, i
c$$$c$$$            work1(kk) = zvec(i,j)
c$$$c$$$            kk = kk + 1
c$$$c$$$ 9200    continue
c$$$c$$$ 9201 continue
c$$$c
c$$$      call eigen ( work1,work, jmax1, 0 )
c$$$      kk = 1
c$$$      do 9210 i = 1, jmax1
c$$$         ii = ( i*  (i+1) ) / 2
c$$$         j = (i-1)*jmax1 + 1
c$$$         jj = j + jmax1 - 1
c$$$c
c$$$         write(outmod,9003) i, work1(ii), ( work(j1), j1=j,jj )
c$$$c
c$$$ 9003    format(1x,/,1x, "eigenvector no.",i4,
c$$$     $        "  eigenvalue = ",e12.5,/, 10(1x,1p10e11.4,/ ) )
c$$$c
c$$$ 9210 continue
c$$$c
c$$$      return
c$$$      end
c$$$c.................................................................
c$$$      subroutine mateig2 ( zvec, nd,msiz, zwk,work0,work,work1, l1,l2,
c$$$     $     ff,lcone, lopt, wkxx,wkyy,jobid, nout )
c$$$c.................................................................
c$$$c
c$$$c... Calls EIGEN to find eigensolutions of ZVEC. Puts eigenvalues in
c$$$c    WORK0, and eigenvectors in zwk.
c$$$c    LCONE = 1: Contour plots of eigenvectors(not working?)
c$$$c    LCONE = 2  Plots eigenvectors versus expansion index.
c$$$c
c$$$      character*(*) jobid
c$$$      dimension zvec(nd,nd), zwk(nd,nd), work(*), work1(*), work0(*)
c$$$      dimension wkxx(*), wkyy(*)
c$$$c
c$$$c     find eigenvalues of zvec...
c$$$c
c$$$      mfsq = msiz*msiz
c$$$c
c$$$c$$$      do 9195 i = 1, mfsq
c$$$c$$$         work1(i) = 0.0
c$$$c$$$ 9195 continue
c$$$c
c$$$      work1(1:msiz*(msiz+1)/2) = (/ ((zvec(i,j),j=1,i),i=1,msiz) /)
c$$$c
c$$$c$$$      kk = 1
c$$$c$$$      do 9201 i = 1, msiz
c$$$c$$$         do 9200 j = 1, i
c$$$c$$$            work1(kk) = zvec(i,j)
c$$$c$$$            kk = kk + 1
c$$$c$$$ 9200    continue
c$$$c$$$ 9201 continue
c$$$c
c$$$      call eigen ( work1,work, msiz, 0 )
c$$$      kk = 1
c$$$      do 9210 i = 1, msiz
c$$$         ii = ( i*  (i+1) ) / 2
c$$$         j = (i-1)*msiz + 1
c$$$         jj = j + msiz - 1
c$$$         work0(i) = work1(ii)
c$$$c
c$$$         if ( lopt .eq. 1 ) then
c$$$            write(nout,9002) i, work1(ii)
c$$$ 9002       format(1x, "eigenvector no.",i4,
c$$$     $           "  eigenvalue = ",e12.5 )
c$$$         end if
c$$$c
c$$$         if ( lopt .eq. 2 ) then
c$$$            write(nout,9003) i, work1(ii), ( work(j1), j1=j,jj )
c$$$ 9003       format(1x,/,1x, "eigenvector no.",i4,
c$$$     $           "  eigenvalue = ",e12.5,/, 10(1x,1p10e11.4,/ ) )
c$$$         end if
c$$$c
c$$$ 9210 continue
c$$$c
c$$$      zwk = reshape ( work(1:msiz*msiz), (/ msiz, msiz /) )
c$$$c
c$$$c$$$      do i = 1, msiz
c$$$c$$$         ii = ( i*(i+1) ) / 2
c$$$c$$$         j = (i-1)*msiz + 1
c$$$c$$$         jj = j + msiz - 1
c$$$c$$$         kk = 1
c$$$c$$$         do j1 = j, jj
c$$$c$$$            zwk(kk,i) = work(j1)
c$$$c$$$            kk = kk + 1
c$$$c$$$         end do
c$$$c$$$      end do
c$$$c
c$$$      if ( lcone .eq. 1 ) then
c$$$         call coneig ( zwk, nd, work0, l1,l2 )
c$$$      end if
c$$$c
c$$$      if ( lcone .eq. 2 ) then
c$$$         call evplts ( zwk, nd,nd, work0,ff, l1,l2, l1,l2,
c$$$     $        wkxx, wkyy, jobid )
c$$$      end if
c$$$c
c$$$      return
c$$$      end
c
c...................................................................
      SUBROUTINE mateig2a ( zvec,nd, msiz, zwk, work0, l1,l2, ff,
     $     lcone, lopt, jobid, locid, nout )
c...................................................................

c...  Calls EIGEN to find eigensolutions of ZVEC. Puts eigenvalues in
c     WORK0, and eigenvectors in zwk.
c     LCONE = 1: Contour plots of eigenvectors(not working?)
c     LCONE = 2  Plots eigenvectors versus expansion index.

      CHARACTER*(*) jobid, locid

      REAL, DIMENSION(*) :: work0
      REAL, DIMENSION(nd,nd) :: zvec, zwk

      REAL, DIMENSION(:), ALLOCATABLE :: work, work1

c     find eigenvalues of zvec...

      mfsq = msiz*msiz

      ALLOCATE ( work(mfsq),work1(mfsq) )

      work1(1:( msiz*(msiz+1) )/2) = (/ ((zvec(i,j),j=1,i),i=1,msiz) /)
c
      CALL eigen ( work1,work, msiz, 0 )

      zwk(1:msiz,1:msiz) =
     $     RESHAPE ( work(1:msiz*msiz), (/ msiz, msiz /) )
      work0(1:msiz) = (/ ( work1( (i*(i+1))/2), i = 1,msiz ) /)

         IF ( lopt == 2 ) THEN
            DO i = 1, msiz
               WRITE( NOUT,9003 ) i, work0(i), ( zwk(j1,i), j1=1,msiz )
 9003          FORMAT( 1x,/,1x, "eigenvector no.",i4,
     $              "  eigenvalue = ",e12.5,/, 10(1x,1p10e11.4,/ ) )
            END DO
         END IF

         IF ( lopt == 1 ) THEN
            WRITE ( NOUT,9002 ) ( work0(i),i = 1, msiz)
 9002       FORMAT ( 1x, "  EIGENVALUES = ",/, (10ES12.5) )
         END IF

      DEALLOCATE ( work1 )

      IF ( lcone .eq. 1 ) THEN
         CALL coneig ( zwk, nd, work0, l1,l2 )
      END IF

      IF ( lcone .eq. 2 ) THEN
         CALL evplts ( zwk, nd,nd, work0, ff, l1,l2, l1,l2, jobid,
     $        locid )
      END IF

      DEALLOCATE ( work )

      END SUBROUTINE mateig2a

c...................................................................
      SUBROUTINE mateig2b ( zvec,nd, zbvc,ndb,msiz, epszb, zwk,
     $     workr,worki, l1,l2, ff, lcone, lopt, jobid, locid, nout,
     $     arclenw, jwal2, dwp, dwe, dwo, dwpe , dwpo  )
c...................................................................

c...  Calls F02BJF to find eigensolutions of ZVEC*X = LAMBDA*ZBVC*X
c.     Puts eigenvalues in
c     WORK0, and eigenvectors in zwk.
c     LCONE = 1: Contour plots of eigenvectors(not working?)
c     LCONE = 2  Plots eigenvectors versus expansion index.

      CHARACTER*(*) jobid, locid

      REAL, DIMENSION(*) :: workr, worki, dwp, dwo, dwe, dwpe, dwpo
      REAL, DIMENSION(nd,nd) :: zvec, zwk
      REAL, DIMENSION(ndb,ndb) :: zbvc
      REAL, DIMENSION(:), ALLOCATABLE :: zalfar, zalfai, zbeta, tmp
      REAL, DIMENSION(:,:), ALLOCATABLE :: zvecsv
      INTEGER, DIMENSION(:), ALLOCATABLE :: izter

      REAL, DIMENSION(:), ALLOCATABLE :: zvr, zvi, zzz1, zzz2
      COMPLEX, DIMENSION(:), ALLOCATABLE :: zvc, zvc2
      COMPLEX :: zphasei
      INTEGER, DIMENSION(1) :: imax

      REAL, DIMENSION(:), ALLOCATABLE :: work, ziter
      LOGICAL :: lzmatv

c     find eigenvalues of zvec...

      mfsq = msiz*msiz
c$$$      jmax1= msiz-jwal2
      jwal4 = 2*jwal2
      jmax2 = msiz - jwal4
      jmax1 = jmax2/2
      jwal1= jwal2/2

      WRITE ( nout, '(/, 4x, "msiz, jmax1, jwal1, jwal2, in mateig2b = ",
     $      4I5,/ )' ) msiz, jmax1, jwal1, jwal2

      ALLOCATE ( work(mfsq) )
      ALLOCATE ( izter(msiz) )
      ALLOCATE ( tmp(msiz), ziter(msiz) )

c      work1(1:( msiz*(msiz+1) )/2) = (/ ((zvec(i,j),j=1,i),i=1,msiz) /)

      ALLOCATE ( zalfar(msiz), zalfai(msiz),zbeta(msiz) )
      ALLOCATE ( zvecsv(msiz,msiz) )

c     save the zvec array

      DO i=1,msiz
         zvecsv(i,1:msiz)=zvec(i,1:msiz)
      END DO

      zeps1 = -1.0
      izfail = 0
      lzmatv = .TRUE.
      izter = 0

      ndzz = nd
c      CALL evalef ( msiz,zvec,nd,zbvc,ndb,zeps1,zalfar,zalfai,
c     $     zbeta,lzmatv,zwk,ndzz,izter,izfail )
      CALL f02bjf ( msiz,zvec,nd,zbvc,ndb,zeps1,zalfar,zalfai,
     $     zbeta,lzmatv,zwk,ndzz,izter,izfail )
      DO i = 1, msiz
         ziter(i) = REAL(izter(i))
      END DO

c     Put back the zvec array

      DO i=1,msiz
         zvec(i,1:msiz)=zvecsv(i,1:msiz)
      END DO

      WRITE ( 6, '(/, "izfail in mateig2b = ", i4,/ )' ) izfail
      WRITE ( NOUT, '(/, "izfail = in mateig2b ", i4,/ )' ) izfail

c... Calculate eigenvalues from (zalfar + i*zalfai )/ zbeta

      zepsb = 1.0e-30
      DO i = 1, msiz
         IF ( ABS(zbeta(i)) <= zepsb ) THEN
            WRITE ( 6,
     $           '(5x, "ZBETA = ", es12.4, " at i = ",i4)' )
     $           zbeta(i), i
            zbeta(i) = sign (zepsb,zbeta(i))
            WRITE ( nout,
     $           '(5x, "ZBETA = ", es12.4, " at i = ",i4)' )
     $           zbeta(i), i
            zbeta(i) = sign (zepsb,zbeta(i))
         END IF
         workr(i) = zalfar(i) / zbeta(i)
         worki(i) = zalfai(i) / zbeta(i)
      END DO

      WRITE ( NOUT, '(/, "***** Before Reorder *****",/ )' )

         ztwopi = 2.0 * ACOS(-1.0)
         IF ( lopt == 1 ) THEN
            WRITE ( NOUT, '(/,3x, "i",2x, "izter", 7x, "ziter", 5x,
     $           "alfar",7x,
     $           "alfai",7x,"beta",8x,"evalr",7x,"evali",7x,"gtauw",/,
     $           (i4,i8, 7es12.4) )' )
     $           ( i, izter(i), ziter(i),zalfar(i),zalfai(i),zbeta(i),
     $           workr(i), worki(i), workr(i)/ztwopi**3*arclenw,
     $           i = 1, msiz )
         END IF

      DO i= 1,jmax2+jwal4
         anorm=DOT_PRODUCT(zwk(1:msiz,i),zwk(1:msiz,i)) - 1.0
c         zwk(1:msiz,i)=zwk(1:msiz,i)/sqrt(anorm)
c         anorm1=DOT_PRODUCT(zwk(jmax2+1:msiz,i),zwk(jmax2+1:msiz,i))
         write(nout,'("#i, zwk*zwk-1 ",i4,es12.3)')i,anorm
         write(   6,'("#i, zwk*zwk-1 ",i4,es12.3)')i,anorm
      END DO

      WRITE ( NOUT, '(/, "Check Orthogonality before reorder",/ )' )

      WRITE(NOUT,'("#i,sump0,sumw0,sumpw0-1, sump,sumw,sumpw = ")')
      DO i= 1, jmax2+jwal4
         ajsump = 0.0
         ajsumw = 0.0
         sump0 = DOT_PRODUCT(zwk(1:jmax2,i),zwk(1:jmax2,i))
         sumw0 = DOT_PRODUCT(zwk(jmax2+1:msiz,i),zwk(jmax2+1:msiz,i))
         sumpw0 = sump0 + sumw0 - 1.0
         DO j = 1, msiz
            IF ( I == J ) GO TO 455
            anorm=DOT_PRODUCT(zwk(1:jmax2,i),zwk(1:jmax2,j))
             ajsump = ajsump + anorm
            anorm=DOT_PRODUCT(zwk(jmax2+1:msiz,i),zwk(jmax2+1:msiz,j))
            ajsumw = ajsumw + anorm
 455        CONTINUE
         END DO
         ajsumpw = epszb*ajsump + ajsumw
         WRITE(NOUT,'("#i, ",i4, 3es12.3, 1x, 3es12.3, es20.11)')
     $        i,sump0,sumw0,sumpw0, ajsump,ajsumw,ajsumpw, workr(i)
      END DO
c
      CALL reorder(zwk,workr,worki,zalfar,zbeta,izter,
     &     nd,msiz,nout)

      WRITE ( NOUT, '(/, "******** After Reorder ********", )' )

         ztwopi = 2.0 * acos(-1.0)
         IF ( lopt == 1 ) THEN
            WRITE ( NOUT, '(/,3x, "i",2x, "izter", 7x, "ziter", 5x,
     $           "alfar",7x,
     $           "alfai",7x,"beta",8x,"evalr",7x,"evali",7x,"gtauw",/,
     $           (i4,i8, 7es12.4) )' )
     $           ( i, izter(i), ziter(i),zalfar(i),zalfai(i),zbeta(i),
     $           workr(i), worki(i), workr(i)/ztwopi**3*arclenw,
     $           i = 1, msiz )
         END IF

      DO i= 1,jmax2+jwal4
         anorm=DOT_PRODUCT(zwk(1:msiz,i),zwk(1:msiz,i)) - 1.0
c         zwk(1:msiz,i)=zwk(1:msiz,i)/sqrt(anorm)
c         anorm1=DOT_PRODUCT(zwk(jmax2+1:msiz,i),zwk(jmax2+1:msiz,i))
         write(nout,'("#i, zwk*zwk-1 ",i4,es12.3)')i,anorm
         write(   6,'("#i, zwk*zwk-1",i4,es12.3)')i,anorm
      END DO

      WRITE ( NOUT, '(/, "Check Orthonormality after reorder",/ )' )

      WRITE(NOUT,'("#i,sump0,sumw0,sumpw0-1, sump,sumw,sumpw = ")')
      DO i= 1, jmax2+jwal4
         ajsump = 0.0
         ajsumw = 0.0
         sump0 = DOT_PRODUCT(zwk(1:jmax2,i),zwk(1:jmax2,i))
         sumw0 = DOT_PRODUCT(zwk(jmax2+1:msiz,i),zwk(jmax2+1:msiz,i))
         sumpw0 = sump0 + sumw0- 1.0
         DO j = 1, msiz
            IF ( I == J ) GO TO 465
            anorm=DOT_PRODUCT(zwk(1:jmax2,i),zwk(1:jmax2,j))
             ajsump = ajsump + anorm
            anorm=DOT_PRODUCT(zwk(jmax2+1:msiz,i),zwk(jmax2+1:msiz,j))
            ajsumw = ajsumw + anorm
 465        CONTINUE
         END DO
         ajsumpw = epszb*ajsump + ajsumw
         WRITE(NOUT,'("#i, ",i4, 3es12.3, 1x, 3es12.3, es20.11)')
     $        i,sump0,sumw0,sumpw0, ajsump,ajsumw,ajsumpw, workr(i)
      END DO

c      check normalization

c$$$      DO i=1,jwal2
c$$$         anorm=DOT_PRODUCT(zwk(jmax1+1:msiz,i),zwk(jmax1+1:msiz,i))
c$$$         zwk(1:msiz,i)=zwk(1:msiz,i)/sqrt(anorm)
c$$$         anorm1=DOT_PRODUCT(zwk(jmax1+1:msiz,i),zwk(jmax1+1:msiz,i))
c$$$         write(nout,'("#i,norm,norm1 ",i4,1p2e12.3)')i,anorm,anorm1
c$$$         write(6 ,'("#i,norm,norm1 ",i4,1p2e12.3)')i,anorm,anorm1
c$$$      END DO
c$$$      DO i=1,jwal4
c$$$         anorm=DOT_PRODUCT(zwk(jmax2+1:msiz,i),zwk(jmax2+1:msiz,i))
c$$$         zwk(1:msiz,i)=zwk(1:msiz,i)/sqrt(anorm)
c$$$         anorm1=DOT_PRODUCT(zwk(jmax2+1:msiz,i),zwk(jmax2+1:msiz,i))
c$$$         write(nout,'("#i,norm,norm1 ",i4,1p2e12.3)')i,anorm,anorm1
c$$$         write(6 ,'("#i,norm,norm1 ",i4,1p2e12.3)')i,anorm,anorm1
c$$$      END DO

c...             Schmidt's  Orthonormalization:

c... Assume that the degeneracies are grouped nicely in ordered pairs,
c    so that we can orthonormalize pair by pair.

c$$$      DO i = 1, msiz-1, 2
c$$$         ainorm = DOT_PRODUCT(zwk(1:msiz,i),zwk(1:msiz,i))
c$$$         work(1:msiz) = zwk(1:msiz,i) / SQRT(ainorm)
c$$$         zwk(1:msiz,i) = work(1:msiz)
c$$$         work(1:msiz) = zwk(1:msiz,i+1) - zwk(1:msiz,i) *
c$$$     $        DOT_PRODUCT(work(1:msiz),zwk(1:msiz,i+1))
c$$$         ainorm = DOT_PRODUCT(work(1:msiz),work(1:msiz))
c$$$         zwk(1:msiz,i+1) = work(1:msiz) / SQRT(ainorm)
c$$$      END DO
      
c     Take the RHS norm into account.

      DO i = 1, msiz-1, 2
         sump0 = DOT_PRODUCT(zwk(1:jmax2,i),zwk(1:jmax2,i+1))
         sumw0 = DOT_PRODUCT(zwk(jmax2+1:msiz,i),
     $        zwk(jmax2+1:msiz,i+1))
         coef12 = epszb*sump0 + sumw0
         sump0 = DOT_PRODUCT(zwk(1:jmax2,i),zwk(1:jmax2,i))
         sumw0 = DOT_PRODUCT(zwk(jmax2+1:msiz,i),zwk(jmax2+1:msiz,i))
         coef12 = coef12 / (epszb*sump0 + sumw0)
         zwk(1:msiz,i+1) = zwk(1:msiz,i+1) - coef12*zwk(1:msiz,i)
         zwk(1:msiz,i+1) = zwk(1:msiz,i+1) /
     $        SQRT( DOT_PRODUCT(zwk(1:msiz,i+1),zwk(1:msiz,i+1)) )
      END DO

      WRITE ( NOUT, '(/, "After Orthonormalization",/ )' )
         
      DO i= 1,jmax2+jwal4
         anorm=DOT_PRODUCT(zwk(1:msiz,i),zwk(1:msiz,i)) - 1.0
c         zwk(1:msiz,i)=zwk(1:msiz,i)/sqrt(anorm)
c         anorm1=DOT_PRODUCT(zwk(jmax2+1:msiz,i),zwk(jmax2+1:msiz,i))
         write(nout,'("#i, zwk*zwk-1 ",i4,es12.3)')i,anorm
         write(   6,'("#i, zwk*zwk-1 ",i4,es12.3)')i,anorm
      END DO

      WRITE ( NOUT, '(/,
     $     "Check Orthonormality after Schmidt, epszb = ", e12.3/ )' )
     $     epszb

      WRITE(NOUT,'("#i,sump0,sumw0,sumpw0-1, sump,sumw,sumpw = ")')
      DO i= 1, jmax2+jwal4
         ajsump = 0.0
         ajsumw = 0.0
         sump0 = DOT_PRODUCT(zwk(1:jmax2,i),zwk(1:jmax2,i))
         sumw0 = DOT_PRODUCT(zwk(jmax2+1:msiz,i),zwk(jmax2+1:msiz,i))
         sumpw0 = sump0 + sumw0 - 1.0
         DO j = 1, msiz
            IF ( I == J ) GO TO 475
            anorm=DOT_PRODUCT(zwk(1:jmax2,i),zwk(1:jmax2,j))
             ajsump = ajsump + anorm
            anorm=DOT_PRODUCT(zwk(jmax2+1:msiz,i),zwk(jmax2+1:msiz,j))
            ajsumw = ajsumw + anorm
 475        CONTINUE
         END DO
         ajsumpw = epszb*ajsump + ajsumw
         WRITE(NOUT,'("#i, ",i4, 3es12.3, 1x, 3es12.3, es20.11)')
     $        i,sump0,sumw0,sumpw0, ajsump,ajsumw,ajsumpw, workr(i)
      END DO

c
c$$$      zwk(1:msiz,1:msiz) =
c$$$     $     RESHAPE ( work(1:msiz*msiz), (/ msiz, msiz /) )
c$$$      work0(1:msiz) = work1( (i*(i+1))/2, i = 1,msiz )

c$$$      kk = 1
c$$$      DO 9210 i = 1, msiz
c$$$         ii = ( i*(i+1) ) / 2
c$$$         j = (i-1)*msiz + 1
c$$$         jj = j + msiz - 1
c$$$         work0(i) = work1(ii)
c$$$
c$$$         IF ( lopt == 2 ) THEN
c$$$            WRITE( NOUT,9003 ) i, work1(ii), ( work(j1), j1=j,jj )
c$$$ 9003       FORMAT( 1x,/,1x, "eigenvector no.",i4,
c$$$     $           "  eigenvalue = ",e12.5,/, 10(1x,1p10e11.4,/ ) )
c$$$         END IF
c$$$c
c$$$ 9210 CONTINUE

c...  Rotate eigenvectors by using a linear combination of degenerate pairs.
c      Use complex arithmetic first to find the component with the 
c      largest magnitude.
c      Look for the largest component.

      msizh = msiz/2
      ALLOCATE ( zvr(msizh), zvi(msizh), zvc(msizh), zvc2(msizh) )
      ALLOCATE ( zzz1(msiz), zzz2(msiz) )
      
      zvr = 0.0
      zvi = 0.0
      zvc = 0.0
      zvc2 = 0.0

      WRITE ( NOUT, '(/, "After Rotation",/ )' )
      
      WRITE(NOUT, '( /,3x, "i", 1x, "imax", 6x, "Z-before", 11x,
     $  "phase", 11x, "Z-after", 13x,
     $     "Tot Re,Im P", 11x, "Tot Re,Im W", 6x, "Total-1.0" )' )

      DO i = 1, msiz, 2

c... Get the phase at i, phi1
c... Restrict to plasma components only.
c$$$         zvr(1:msizh) = (/ zwk(1:jmax1,i),zwk(jmax2+1:jmax2+jwal2,i) /)
c$$$         zvi(1:msizh) = (/ zwk(jmax1+1:jmax2,i),
c$$$     $        zwk(jmax2+jwal2+1:msiz,i) /)
         zvr(1:jmax1) = zwk(1:jmax1,i)
         zvi(1:jmax1) = zwk(jmax1+1:jmax2,i)
         zvc = CMPLX ( zvr, zvi )
c... Plasma components only:
c$$$         imax = MAXLOC ( ABS(zvc(:) ) )
         imax = MAXLOC ( ABS(zvc(1:jmax1) ) )
         ix1 = imax(1)
         zvcor1 = REAL(zvc(imax(1)))      ! First of Pair: Real part.
         zvcoi1 = AIMAG(zvc(imax(1)))    ! Imag part
         zphi1 = ATAN2 ( zvi(ix1),zvr(ix1) )  ! Phase of first of pair

c... now get the phase at i+1, phi2
c  No need to change anything here to deal with plasma components only.
         zvr(1:msizh) =
     $        (/ zwk(1:jmax1,i+1),zwk(jmax2+1:jmax2+jwal2,i+1) /)
         zvi(1:msizh) = (/ zwk(jmax1+1:jmax2,i+1),
     $        zwk(jmax2+jwal2+1:msiz,i+1) /)
         zvc = CMPLX ( zvr, zvi )
         zvcor2 = REAL(zvc(imax(1)))     ! Second of Pair: Real part
         zvcoi2 = AIMAG(zvc(imax(1)))    ! Imag part
         zphi2 = ATAN2 ( zvi(ix1),zvr(ix1) )  ! Phase

         sphi1 = SIN(zphi1)
         sphi2 = SIN(zphi2)
         zzz1(1:msiz) = zwk(1:msiz,i)*sphi2 - zwk(1:msiz,i+1)*sphi1
         zzz2(1:msiz) = zwk(1:msiz,i)*sphi1 + zwk(1:msiz,i+1)*sphi2
         zwk(1:msiz,i) = zzz1(1:msiz)
         zwk(1:msiz,i+1) = zzz2(1:msiz)
 
         IF ( ix1 <= jmax1 ) THEN
            zvcfr = zwk(ix1,i)         ! Real part after rotation at i
            zvcfi = zwk(jmax1+ix1,i)   ! Imag part after rotation at i
         ELSE
            zvcfr = zwk(jmax1+ix1,i)
            zvcfi = zwk(jmax1+jwal2+ix1,i)
         END IF
c... Ensure that the maximum component is positive. Then change all signs.
         IF ( zvcfr < 0.0 ) THEN
            zvcfr = - zvcfr
            zvcfi = - zvcfi
            zwk(1:msiz,i) = - zwk(1:msiz,i)
c$$$            zwk(1:msiz,i+1) = - zwk(1:msiz,i+1) ! Only the 1st mode. 080807
         END IF
c.. Sum of the real and imaginary elements:
         zrsump = SUM(ABS(zwk(1:jmax1,i))**2)
         zrsumw = SUM(ABS(zwk(jmax2+1:jmax2+jwal2,i))**2)
         zisump = SUM(ABS(zwk(jmax1+1:jmax2,i))**2) 
         zisumw = SUM(ABS(zwk(jmax2+jwal2+1:jmax2+jwal4,i))**2)
         zsumt = zrsump + zrsumw + zisump + zisumw - 1.0
         WRITE ( nout, '(2i4, 10es11.3 )' )
     $        i, imax, zvcor1, zvcoi1, zphi1, zvcfr, zvcfi,
     $        zrsump, zisump, zrsumw, zisumw, zsumt
c  now, i+1:
         IF ( ix1 <= jmax1 ) THEN
            zvcfr = zwk(ix1,i+1)
            zvcfi = zwk(jmax1+ix1,i+1)
         ELSE
            zvcfr = zwk(jmax1+ix1,i+1)
            zvcfi = zwk(jmax1+jwal2+ix1,i+1)
         END IF
c.   Ensure Imaginary component for 2nd mode is positive .....Aug. 08, 2007
c... Ensure that the maximum component is positive. Then change all signs.
         IF ( zvcfi < 0.0 ) THEN
            zvcfr = - zvcfr
            zvcfi = - zvcfi
c$$$            zwk(1:msiz,i) = - zwk(1:msiz,i)    ! only the 2nd mode
            zwk(1:msiz,i+1) = - zwk(1:msiz,i+1)
         END IF

         zrsump = SUM(ABS(zwk(1:jmax1,i+1))**2) 
         zrsumw = SUM(ABS(zwk(jmax2+1:jmax2+jwal2,i+1))**2)
         zisump = SUM(ABS(zwk(jmax1+1:jmax2,i+1))**2 )
         zisumw = SUM(ABS(zwk(jmax2+jwal2+1:jmax2+jwal4,i+1))**2)
         zsumt = zrsump + zrsumw + zisump + zisumw - 1.0
         WRITE ( nout, '(2i4, 10es11.3 )' )
     $        i, ix1, zvcor2, zvcoi2, zphi2, zvcfr, zvcfi,
     $        zrsump, zisump, zrsumw, zisumw, zsumt

      END DO
      
      DEALLOCATE ( zvr,zvi,zvc,zvc2, zzz1, zzz2 )

      WRITE ( NOUT, '(/, "After Orthonormalization and Rotation",/ )' )
         
      DO i= 1,jmax2+jwal4
         anorm=DOT_PRODUCT(zwk(1:msiz,i),zwk(1:msiz,i)) - 1.0
c         zwk(1:msiz,i)=zwk(1:msiz,i)/sqrt(anorm)
c         anorm1=DOT_PRODUCT(zwk(jmax2+1:msiz,i),zwk(jmax2+1:msiz,i))
         write(nout,'("#i, zwk*zwk-1 ",i4,es12.3)')i,anorm
         write(   6,'("#i, zwk*zwk-1 ",i4,es12.3)')i,anorm
      END DO

      WRITE ( NOUT, '(/,
     $     "Orthonorm.: Schmidt and rotated, epszb = ", e12.3/ )')
     $     epszb

      WRITE(NOUT,'("#i,sump0,sumw0,sumpw0-1, sump,sumw,sumpw = ")')
      DO i= 1, jmax2+jwal4
         ajsump = 0.0
         ajsumw = 0.0
         sump0 = DOT_PRODUCT(zwk(1:jmax2,i),zwk(1:jmax2,i))
         sumw0 = DOT_PRODUCT(zwk(jmax2+1:msiz,i),zwk(jmax2+1:msiz,i))
         sumpw0 = sump0 + sumw0 - 1.0
         DO j = 1, msiz
            IF ( I == J ) GO TO 485
            anorm=DOT_PRODUCT(zwk(1:jmax2,i),zwk(1:jmax2,j))
             ajsump = ajsump + anorm
            anorm=DOT_PRODUCT(zwk(jmax2+1:msiz,i),zwk(jmax2+1:msiz,j))
            ajsumw = ajsumw + anorm
 485        CONTINUE
         END DO
         ajsumpw = epszb*ajsump + ajsumw
         WRITE(NOUT,'("#i, ",i4, 3es12.3, 1x, 3es12.3, es20.11)')
     $        i,sump0,sumw0,sumpw0, ajsump,ajsumw,ajsumpw, workr(i)
      END DO

c     the energies of normalized displacements. i, i+1 separately.

c.. First degenerate pair:

      WRITE ( NOUT, '(/,"......First of the degenerate pair...." )' )
      WRITE (NOUT, '("  i",6x,"dwp",9x,"dwpi",8x,"dwe",9x,"dwo",9x,
     $     "dwpe",8x,"dwpo", 9x,"dwtot",5x,"RHS-norm",5x,"dwtotn")')
      DO i=1,msiz, 2
         izwrt = (i+1)/2
         dwp(i)=0.
         tmp(1:msiz)=0.
         DO j=1,jmax2
            tmp(j)=DOT_PRODUCT(zvec(j,1:jmax2),zwk(1:jmax2,i))
         END DO
         dwp(i)=DOT_PRODUCT(tmp(1:jmax2),zwk(1:jmax2,i))
c
         dwpi=0.
         tmp(1:msiz)=0.
         DO j=jmax1+1,jmax2
            tmp(j)=DOT_PRODUCT(zvec(j,jmax1+1:jmax2),
     $           zwk(jmax1+1:jmax2,i))
         END DO
         dwpi=DOT_PRODUCT(tmp(jmax1+1:jmax2),zwk(jmax1+1:jmax2,i))
c
         dwe(i)=0.
         tmp(1:msiz)=0.
         DO j=jmax2+1,jmax2+jwal2
            tmp(j)=DOT_PRODUCT(zvec(j,jmax2+1:jmax2+jwal2)
     $        ,zwk(jmax2+1:jmax2+jwal2,i))
         END DO
         dwe(i)=DOT_PRODUCT(tmp(jmax2+1:jmax2+jwal2)
     $        ,zwk(jmax2+1:jmax2+jwal2,i))
c
         dwo(i)=0.
         tmp(1:msiz)=0.
         DO j=jmax2+jwal2+1,jmax2+jwal4
            tmp(j)=DOT_PRODUCT(zvec(j,jmax2+jwal2+1:jmax2+jwal4)
     $        ,zwk(jmax2+jwal2+1:jmax2+jwal4,i))
         END DO
         dwo(i)=DOT_PRODUCT(tmp(jmax2+jwal2+1:jmax2+jwal4)
     $        ,zwk(jmax2+jwal2+1:jmax2+jwal4,i))
c
         dwpe(i)=0.
         tmp(1:msiz)=0.
         DO j=1,jmax2
            tmp(j)=DOT_PRODUCT(zvec(j,jmax2+1:jmax2+jwal2)
     $        ,zwk(jmax2+1:jmax2+jwal2,i))
         END DO
         dwpe(i)=DOT_PRODUCT(tmp(1:jmax2)
     $        ,zwk(1:jmax2,i))
c
         dwpo(i)=0.
         tmp(1:msiz)=0.
         DO j=1,jmax2
            tmp(j)=DOT_PRODUCT(zvec(j,jmax2+jwal2+1:jmax2+jwal4)
     $        ,zwk(jmax2+jwal2+1:jmax2+jwal4,i))
         END DO
         dwpo(i)=DOT_PRODUCT(tmp(1:jmax2)
     $        ,zwk(1:jmax2,i))

         dwtot = dwp(i) + dwe(i) + dwo(i) + 2.0*( dwpe(i) + dwpo(i) )

c.... Get RHS norm:
         ajsumpw=epszb * DOT_PRODUCT(zwk(1:jmax2,i),zwk(1:jmax2,i))
     $        + DOT_PRODUCT(zwk(jmax2+1:msiz,i),zwk(jmax2+1:msiz,i))
         dwtotn = dwtot / ajsumpw

         WRITE (NOUT,'(i5, 9es12.5)')izwrt,dwp(i),dwpi,dwe(i)
     $        ,dwo(i),dwpe(i),dwpo(i), dwtot, ajsumpw, dwtotn
      END DO

c.. Next degenerate pair:

      WRITE ( NOUT, '(/,".....Second of the degenerate pair......." )' )
      WRITE (NOUT, '("  i",6x,"dwp",9x,"dwpi",8x,"dwe",9x,"dwo",9x,
     $     "dwpe",8x,"dwpo",9x,"dwtot",5x,"RHS-norm",5x,"dwtotn")')
      DO i=2,msiz, 2
         izwrt = i/2
         dwp(i)=0.
         tmp(1:msiz)=0.
         DO j=1,jmax2
            tmp(j)=DOT_PRODUCT(zvec(j,1:jmax2),zwk(1:jmax2,i))
         END DO
         dwp(i)=DOT_PRODUCT(tmp(1:jmax2),zwk(1:jmax2,i))
c
         dwpi=0.
         tmp(1:msiz)=0.
         DO j=jmax1+1,jmax2
            tmp(j)=DOT_PRODUCT(zvec(j,jmax1+1:jmax2),
     $           zwk(jmax1+1:jmax2,i))
         END DO
         dwpi=DOT_PRODUCT(tmp(jmax1+1:jmax2),zwk(jmax1+1:jmax2,i))
c
         dwe(i)=0.
         tmp(1:msiz)=0.
         DO j=jmax2+1,jmax2+jwal2
            tmp(j)=DOT_PRODUCT(zvec(j,jmax2+1:jmax2+jwal2)
     $        ,zwk(jmax2+1:jmax2+jwal2,i))
         END DO
         dwe(i)=DOT_PRODUCT(tmp(jmax2+1:jmax2+jwal2)
     $        ,zwk(jmax2+1:jmax2+jwal2,i))
c
         dwo(i)=0.
         tmp(1:msiz)=0.
         DO j=jmax2+jwal2+1,jmax2+jwal4
            tmp(j)=DOT_PRODUCT(zvec(j,jmax2+jwal2+1:jmax2+jwal4)
     $        ,zwk(jmax2+jwal2+1:jmax2+jwal4,i))
         END DO
         dwo(i)=DOT_PRODUCT(tmp(jmax2+jwal2+1:jmax2+jwal4)
     $        ,zwk(jmax2+jwal2+1:jmax2+jwal4,i))
c
         dwpe(i)=0.
         tmp(1:msiz)=0.
         DO j=1,jmax2
            tmp(j)=DOT_PRODUCT(zvec(j,jmax2+1:jmax2+jwal2)
     $        ,zwk(jmax2+1:jmax2+jwal2,i))
         END DO
         dwpe(i)=DOT_PRODUCT(tmp(1:jmax2)
     $        ,zwk(1:jmax2,i))
c
         dwpo(i)=0.
         tmp(1:msiz)=0.
         DO j=1,jmax2
            tmp(j)=DOT_PRODUCT(zvec(j,jmax2+jwal2+1:jmax2+jwal4)
     $        ,zwk(jmax2+jwal2+1:jmax2+jwal4,i))
         END DO
         dwpo(i)=DOT_PRODUCT(tmp(1:jmax2)
     $        ,zwk(1:jmax2,i))

         dwtot = dwp(i) + dwe(i) + dwo(i) + 2.0*( dwpe(i) + dwpo(i) )

c.... Get RHS norm:
         ajsumpw=epszb * DOT_PRODUCT(zwk(1:jmax2,i),zwk(1:jmax2,i))
     $        + DOT_PRODUCT(zwk(jmax2+1:msiz,i),zwk(jmax2+1:msiz,i))
         dwtotn = dwtot / ajsumpw
         
         WRITE (NOUT,'(i5, 9es12.5)')izwrt,dwp(i),dwpi,dwe(i)
     $        ,dwo(i),dwpe(i),dwpo(i),dwtot, ajsumpw, dwtotn
      END DO

c..  Now renormalize the eigenvectors for the feedback:
      DO i=1,jwal4
         anorm=DOT_PRODUCT(zwk(jmax2+1:msiz,i),zwk(jmax2+1:msiz,i))
         zwk(1:msiz,i)=zwk(1:msiz,i)/sqrt(anorm)
         anorm1=DOT_PRODUCT(zwk(jmax2+1:msiz,i),zwk(jmax2+1:msiz,i))
         write(nout,'("#i,norm,norm1 ",i4,1p2e12.5)')i,anorm,anorm1
         write(6 ,'("#i,norm,norm1 ",i4,1p2e12.5)')i,anorm,anorm1
      END DO

         IF ( lopt == 2 ) THEN
            DO i = 1, msiz
               WRITE( NOUT,9003 )
     $              i, workr(i), worki(i), ( zwk(j1,i), j1=1,msiz )
 9003          FORMAT( 1x,/,1x, "eigenvector no.",i4,
     $              "  eigenvalue = ",2es12.5,/, 10(1x,1p10e11.4,/ ) )
            END DO
         END IF

c$$$         IF ( lopt == 1 ) THEN
c$$$            WRITE ( NOUT,9002 ) ( work0(i),i = 1, msiz)
c$$$ 9002       FORMAT ( 1x, "  EIGENVALUES = ",/, (10ES12.5) )
c$$$         END IF

c$$$      zwk(1:msiz,1:msiz) =
c$$$     $     RESHAPE ( work(1:msiz*msiz), (/ msiz, msiz /) )
c
c$$$  do i = 1, msiz
c$$$  ii = ( i*(i+1) ) / 2
c$$$  j = (i-1)*msiz + 1
c$$$  jj = j + msiz - 1
c$$$  kk = 1
c$$$  do j1 = j, jj
c$$$  zwk(kk,i) = work(j1)
c$$$  kk = kk + 1
c$$$  end do
c$$$  end do

      IF ( lcone .eq. 1 ) THEN
         CALL coneig ( zwk, nd, workr, l1,l2 )
      END IF

      IF ( lcone .eq. 2 ) THEN
         ly1 = 1      ! plot plasma portion first
         ly2 = jmax1
         CALL evplts ( zwk, nd,nd, workr, ff, l1,l2, ly1,ly2, jobid,
     $        locid )
         ly1 = jmax1 + 1      !  now plot wall portion
         ly2 = msiz
         CALL evplts ( zwk, nd,nd, workr, ff, l1,l2, ly1,ly2, jobid,
     $        locid )

      END IF

      DEALLOCATE ( zalfar, zalfai, zbeta )
      DEALLOCATE ( work, izter, ziter )
      DEALLOCATE ( zvecsv, tmp )

      END SUBROUTINE mateig2b

c...................................................................
      subroutine evplts ( eigmat, ndm1,ndm2, eigval, ff,
     $     lx1,lx2, lz1,lz2, jobid, locid )
c..................................................................
c
c..   ndm1,ndm2: dimensions of eigenvector (stored as columns) matrix, eigmat
c     l1, l2: start and end of eigenvector abscissa.
c
      DIMENSION eigmat(ndm1,ndm2), eigval(*)
      REAL, DIMENSION(:), ALLOCATABLE :: xxxx, yeig
      DIMENSION c(201)
      CHARACTER*80 strg
      CHARACTER*(*) jobid, locid
c
c.....plot each eigenvector as function of Fourier index
c
      jmax1 = lz2 - lz1 + 1
      msiz = jmax1
      nx = lx2 - lx1 + 1
c
      zxmin = lx1 - 0.1
      zymin = -0.1
      zxmax = lx2 + 0.1
      zymax = msiz + 1.1
c
      call mapg(zxmin,zxmax, zymin,zymax ,0.1,0.858,0.10,1.0)
c
      ALLOCATE ( xxxx(nx), yeig(nx) )

      do ix = 1, nx
         xxxx(ix) = lx1 - 1 + ix
      end do
c
      do  i = 1, msiz
         do j1 = 1, nx
            yeig(j1) = eigmat(j1,i)/2.0 + i
         end do
         call trace ( xxxx, yeig, nx, 1,1,0.,0. )
      end do

      DEALLOCATE ( xxxx, yeig )
c
c     if ( l1*l2 .lt. 0 ) then
c     call linep ( 0.0, zymin, 0.0, zymax, 2 )
c     end if
c
      xwrt = zxmax + (zxmax-zxmin)/40.0
      ywrt = zymax - (zymax-zymin)/6.0
      call setlch ( xwrt,ywrt,0,0,0,-1)
      write (strg, 101 ) locid
      call gtext(strg,-1,-1)
 101  format ( a )
      do 17 i = 1, jmax1
         iw = jmax1 + 1 - i
         write(strg,16) iw, eigval(iw)
         call gtext(strg,-1,-1)
 16      format ( i3, 1pe12.4 )
 17   continue
c
      call framep( jobid, ff )
c
      return
      end
c
c...........................................................
      SUBROUTINE reorder(zwk,workr,worki,alfar,beta,izter,
     &     nd,msiz,outmod)
c..............................................................

c     reorder the elements in workr, worki and zwk, so that the
c     smallest of the value of workr is the first one

      IMPLICIT NONE
      INTEGER               ::nd,msiz,i,j1,k,outmod,iztertmp
      INTEGER,DIMENSION(1)  ::j
      INTEGER,DIMENSION(*)  ::izter
      REAL, DIMENSION(nd,nd)::zwk
      REAL, DIMENSION(*)    ::workr,worki,alfar,beta
      REAL                  ::tmpr,tmpi,alfartmp,betatmp
      REAL,DIMENSION(:),ALLOCATABLE::atmp
      LOGICAL,DIMENSION(:),ALLOCATABLE::mask
c
      ALLOCATE(atmp(1:msiz))
      ALLOCATE(mask(1:msiz))
c
c     write the input before being ordered
c
c$$$      DO i=1,msiz
c$$$         WRITE (outmod, '(/,"input to reorder",/)')
c$$$         WRITE (outmod, '("i=",i5," iz=", i5," rv=",1pe12.3,
c$$$     &        " iv=",1pe12.3)')i,izter(i),workr(i),worki(i)
c$$$         WRITE (outmod, '(1p6e12.3)')(zwk(j1,i),j1=1,msiz)
c$$$      END DO
c
      tmpr=1.e30
c      write (*,'("reorder i workr",i5,1pe12.3)')(i,workr(i),i=1,msiz)
      DO i=1,msiz
c         if ( workr(i).lt.tmpr) then
c            j=MAXLOC(workr(1:msiz))
         j=MAXLOC(workr(1:msiz-i+1),mask=(workr(1:msiz-i+1).le.tmpr))
c         end if
c         write (*,'("reorder i,j1, tmpr",2i5,1pe12.3)')i,j(1),tmpr
         j1=j(1)
         tmpr=workr(j1)
         tmpi=worki(j1)
         alfartmp=alfar(j1)
         betatmp=beta(j1)
         iztertmp=izter(j1)
         atmp(1:msiz)=zwk(1:msiz,j1)
         DO k=j1,msiz-i
            workr(k)=workr(k+1)
            worki(k)=worki(k+1)
            alfar(k)=alfar(k+1)
            beta (k)=beta (k+1)
            izter(k)=izter(k+1)
            zwk(1:msiz,k)=zwk(1:msiz,k+1)
         END DO
         workr(msiz-i+1)=tmpr
         worki(msiz-i+1)=tmpi
         alfar(msiz-i+1)=alfartmp
         beta (msiz-i+1)=betatmp
         izter(msiz-i+1)=iztertmp
         zwk(1:msiz,msiz-i+1)=atmp(1:msiz)
      END DO
c
c     write the output after being ordered
c
c$$$      DO i=1,msiz
c$$$         WRITE (outmod, '(/,"output from reorder",/)')
c$$$         WRITE (outmod, '("i=",i5," iz=", i5," rv=",1pe12.3,
c$$$     &        " iv=",1pe12.3)')i,izter(i),workr(i),worki(i)
c$$$         WRITE (outmod, '(1p6e12.3)')(zwk(j1,i),j1=1,msiz)
c$$$      END DO
c
      DEALLOCATE(atmp,mask)
c
      RETURN
      END SUBROUTINE reorder

c....................................................................
      SUBROUTINE mateig2c ( zvec,nd, zbvc,ndb,msiz, zwk, workr,worki,
     $     l1,l2, ff, lcone, lopt, jobid, locid, nout, arclenw,
     $     jwal2 ,gmat)
c...................................................................

c...  Calls F02BJF to find eigensolutions of ZVEC*X = LAMBDA*ZBVC*X
c.     Puts eigenvalues in
c     WORK0, and eigenvectors in zwk.
c     LCONE = 1: Contour plots of eigenvectors(not working?)
c     LCONE = 2  Plots eigenvectors versus expansion index.

      CHARACTER*(*) jobid, locid

      REAL, DIMENSION(*) :: workr, worki
      REAL, DIMENSION(nd,nd) :: zvec, zwk
      REAL, DIMENSION(ndb,ndb) :: zbvc
      REAL, DIMENSION(:), ALLOCATABLE :: zalfar, zalfai, zbeta
      INTEGER, DIMENSION(:), ALLOCATABLE :: izter
      INTEGER,DIMENSION(1)     ::k 

      REAL, DIMENSION(:), ALLOCATABLE :: work
      LOGICAL :: lzmatv
      REAL, DIMENSION(:,:),ALLOCATABLE::oplw,opln,orthm,zbvcn
      REAL, DIMENSION(:,:),ALLOCATABLE::zwkn
      REAL, DIMENSION(:),  ALLOCATABLE::eigwal
c 
c     find eigenvalues of zvec...

      jmax1=msiz-jwal2

      mfsq = msiz*msiz

      ALLOCATE ( work(mfsq) )
      ALLOCATE ( izter(msiz) )
      ALLOCATE ( eigwal(msiz) )
c
      ALLOCATE ( oplw(jwal2,jwal2) )
      ALLOCATE ( opln(jwal2,jwal2) )
      ALLOCATE ( orthm(jwal2,jwal2) )
      ALLOCATE ( zbvcn(jwal2,jwal2) )
      ALLOCATE ( zwkn(jwal2,jwal2) )

c      work1(1:( msiz*(msiz+1) )/2) = (/ ((zvec(i,j),j=1,i),i=1,msiz) /)

      ALLOCATE ( zalfar(msiz), zalfai(msiz),zbeta(msiz) )

c     store away the matrix so that we may use it later
      do i=1,jwal2
         oplw(1:jwal2,i)=zvec(jmax1+1:msiz,jmax1+i)
      end do
c
      zeps1 = -1.0
      izfail = 0
      lzmatv = .TRUE.
      izter = 0

c      CALL evalef ( msiz,zvec,nd,zbvc,ndb,zeps1,zalfar,zalfai,
c     $     zbeta,lzmatv,zwk,nd,izter,izfail )
      CALL f02bjf ( msiz,zvec,nd,zbvc,ndb,zeps1,zalfar,zalfai,
     $     zbeta,lzmatv,zwk,nd,izter,izfail )

      WRITE ( 6, '(/, "izfail in mateig2c = ", i4,/ )' ) izfail
      WRITE ( NOUT, '(/, "izfail = in mateig2c ", i4,/ )' ) izfail

c... Calculate eigaenvalues from (zalfar + i*zalfai )/ zbeta

      zepsb = 1.0e-30
      DO i = 1, msiz
         IF ( ABS(zbeta(i)) <= zepsb ) THEN
            WRITE ( 6,
     $           '(5x, "ZBETA = ", es12.4, " at i = ",i4)' )
     $           zbeta(i), i
            zbeta(i) = sign (zepsb,zbeta(i))
         END IF
         workr(i) = zalfar(i) / zbeta(i)
         worki(i) = zalfai(i) / zbeta(i)
         write (nout,'("i,eigoriginal",i4,1p2e12.3)')i,workr(i),worki(i)
      END DO
c
c     use the unstable wall mode eigenfunction to create the orthogonality
c     matrix
c
c      gmat=1000000.
      tmpr=1.e30
      k=MINLOC(workr(1:msiz),mask=workr(1:msiz).lt.tmpr)
      jmin=k(1)
      eigwal(1:jwal2)=zwk(jmax1+1:msiz,jmin)
      workrs=workr(jmin)
      workis=worki(jmin)
      zalfars=zalfar(jmin)
      zalfais=zalfai(jmin)
      zbetas=zbeta(jmin)
      izters=izter(jmin)
c
      k=MINLOC(eigwal(1:jwal2),mask=eigwal(1:jwal2).lt.tmpr)
      kmin=k(1)
      k=MAXLOC(eigwal(1:jwal2),mask=eigwal(1:jwal2).lt.tmpr)
      kmax=k(1)
      eigmax=eigwal(kmax)
      if (abs(eigwal(kmin)).gt.abs(eigmax))then
         kmax=kmin
         eigmax=eigwal(kmax)
      end if
c
      do i=1,jwal2
         opln(1:jwal2,i)=oplw(1:jwal2,i)
         zbvcn(1:jwal2,i)=0.
         zbvcn(i,1:jwal2)=0.
         zbvcn(i,i)=1.
      end do
c
      do i=1,jwal2
         do j=1,jwal2
            opln(i,j)=opln(i,j)-oplw(i,kmax)*eigwal(j)/eigmax
     $           -oplw(kmax,j)*eigwal(i)/eigmax
     $           +oplw(kmax,kmax)*eigwal(i)*eigwal(j)/eigmax**2
            zbvcn(i,j)=zbvcn(i,j)+eigwal(i)*eigwal(j)/eigmax**2
         end do
      end do
      do i=1,jwal2
         opln(i,kmax)=0.
         opln(kmax,i)=0.
         zbvcn(i,kmax)=0.
         zbvcn(kmax,i)=0.
      end do
         opln(kmax,kmax)=1.
         zbvcn(kmax,kmax)=1.e-10
c         
      call vacasym(opln,jwal2,jwal2,"OPLN",nout,6)
      call vacasym(zbvcn,jwal2,jwal2,"ZBVCN",nout,6)
c
      zeps1 = -1.0
      izfail = 0
      lzmatv = .TRUE.
      izter(1:msiz) = 0
c
c      CALL evalef ( jwal2,opln,jwal2,zbvcn,jwal2,zeps1,zalfar,zalfai,
c     $     zbeta,lzmatv,zwkn,jwal2,izter,izfail )
      CALL f02bjf ( jwal2,opln,jwal2,zbvcn,jwal2,zeps1,zalfar,zalfai,
     $     zbeta,lzmatv,zwkn,jwal2,izter,izfail )
c
c
      DO i = 1, jwal2
         IF ( ABS(zbeta(i)) <= zepsb ) THEN
            WRITE ( 6,
     $           '(5x, "ZBETA = ", es12.4, " at i = ",i4)' )
     $           zbeta(i), i
            zbeta(i) = sign (zepsb,zbeta(i))
         END IF
         sumchk=DOT_PRODUCT(eigwal(1:jwal2),zwkn(1:jwal2,i))
         zwkn(kmax,i)=-sumchk/eigmax
         workr(i) = zalfar(i) / zbeta(i)
         worki(i) = zalfai(i) / zbeta(i)
c         workr1=workr(i)-sumchk**2*gmat
         write (nout,'("i,eigc,ortho,amp",i4,1p3e12.3)')i,workr(i),
     $        sumchk,zwkn(kmax,i)
c         workr(i)=workr1
      END DO
c
      call reorder(zwkn,workr,worki,zalfar,zbeta,izter,
     &     jwal2,jwal2,nout)
c      transfer and shift the eigenfunctions anc eigenvalues
      zwk(1:msiz,1)=zwk(1:msiz,jmin)
c
      do i=1,jwal2-1
         i1=i+1
         zwk(1:jmax1,i1)=0.
         zwk(jmax1+1:msiz,i1)=zwkn(1:jwal2,i)
      end do
      do i=jwal2,2,-1
         workr(i)=workr(i-1)
         zalfar(i)=zalfar(i-1)
         zbeta(i)=zbeta(i-1)
         worki(i)=worki(i-1)
         izter(i)=izter(i-1)
      end do
c      check normalization
      do i=1,jwal2
         anorm=DOT_PRODUCT(zwk(jmax1+1:msiz,i),zwk(jmax1+1:msiz,i))
         zwk(1:msiz,i)=zwk(1:msiz,i)/sqrt(anorm)
         anorm1=DOT_PRODUCT(zwk(jmax1+1:msiz,i),zwk(jmax1+1:msiz,i))
         write(nout,'("#i,norm,norm1 ",i4,1p2e12.3)')i,anorm,anorm1
         write(6 ,'("#i,norm,norm1 ",i4,1p2e12.3)')i,anorm,anorm1
      end do
      workr(1)=workrs
      worki(1)=workis
      zalfar(1)=zalfars
      zalfai(1)=zalfais
      zbeta(1)=zbetas
      izter(1)=izters
c
c$$$      zwk(1:msiz,1:msiz) =
c$$$     $     RESHAPE ( work(1:msiz*msiz), (/ msiz, msiz /) )
c$$$      work0(1:msiz) = work1( (i*(i+1))/2, i = 1,msiz )

c$$$      kk = 1
c$$$      DO 9210 i = 1, msiz
c$$$         ii = ( i*(i+1) ) / 2
c$$$         j = (i-1)*msiz + 1
c$$$         jj = j + msiz - 1
c$$$         work0(i) = work1(ii)
c$$$
c$$$         IF ( lopt == 2 ) THEN
c$$$            WRITE( NOUT,9003 ) i, work1(ii), ( work(j1), j1=j,jj )
c$$$ 9003       FORMAT( 1x,/,1x, "eigenvector no.",i4,
c$$$     $           "  eigenvalue = ",e12.5,/, 10(1x,1p10e11.4,/ ) )
c$$$         END IF
c$$$c
c$$$ 9210 CONTINUE

         IF ( lopt == 2 ) THEN
            DO i = 1, jwal2
               WRITE( NOUT,9003 )
     $              i, workr(i), worki(i), ( zwkn(j1,i), j1=1,jwal2 )
 9003          FORMAT( 1x,/,1x, "eigenvector no.",i4,
     $              "  eigenvalue = ",2es12.5,/, 10(1x,1p10e11.4,/ ) )
            END DO
         END IF
         twopi=2*acos(-1.)
         IF ( lopt == 1 ) THEN
            WRITE ( NOUT, '(/,3x, "i",2x, "izter", 5x,"alfar",7x,
     $           "alfai",7x,"beta",8x,"evalr",7x,"evali",7x,"gtauw"/,
     $           (i4,i8, 6es12.4) )' )
     $           ( i, izter(i), zalfar(i),zalfai(i),zbeta(i),
     $         workr(i), worki(i),workr(i)/twopi**3*arclenw,i=1,jwal2)
            WRITE ( NOUT,'(/,"first eigenfunction",/)')
            WRITE ( NOUT,'(1p5e12.3)')(zwk(i,1),i=1,msiz)
         END IF

c$$$         IF ( lopt == 1 ) THEN
c$$$            WRITE ( NOUT,9002 ) ( work0(i),i = 1, msiz)
c$$$ 9002       FORMAT ( 1x, "  EIGENVALUES = ",/, (10ES12.5) )
c$$$         END IF

c$$$      zwk(1:msiz,1:msiz) =
c$$$     $     RESHAPE ( work(1:msiz*msiz), (/ msiz, msiz /) )
c
c$$$  do i = 1, msiz
c$$$  ii = ( i*(i+1) ) / 2
c$$$  j = (i-1)*msiz + 1
c$$$  jj = j + msiz - 1
c$$$  kk = 1
c$$$  do j1 = j, jj
c$$$  zwk(kk,i) = work(j1)
c$$$  kk = kk + 1
c$$$  end do
c$$$  end do

      IF ( lcone .eq. 1 ) THEN
         CALL coneig ( zwk, nd, workr, l1,l2 )
      END IF

      IF ( lcone .eq. 2 ) THEN
         CALL evplts ( zwk, nd,nd, workr, ff, l1,l2, l1,l2, jobid,
     $        locid )
      END IF

      DEALLOCATE ( zalfar, zalfai, zbeta )
      DEALLOCATE ( work, izter )

      END SUBROUTINE mateig2c

c.............................................................
      SUBROUTINE chiplasma ( blrp, blip, blrw,bliw,chir,chii )
c............................................................
c     
      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum5.inc'
      INCLUDE 'vacuum7.inc'
c     
      DIMENSION xpla(nths), zpla(nths)

      DIMENSION grdgre(nths2,nths2),grwp(nths,nths),
     .     grri(nths2,nfm2)

      REAL, DIMENSION(*) :: blrp, blip, blrw, bliw, chir, chii
      REAL, DIMENSION(:), ALLOCATABLE:: xpb,zpb,cwrkr,cwrki
      REAL, DIMENSION(:,:), ALLOCATABLE:: zchc, zchs
c     
      COMMON / bigv1 / grdgre, grwp, grri
c     
c     
      EQUIVALENCE (xpla,xinf), (zpla,zinf)

      CHARACTER(130), DIMENSION(10) :: string

c     mth is the number of observer points.

      jmax1 = lmax(1) - lmin(1) + 1
      mth12 = 2*mth

      IF ( lfarw .GT. 0 ) mth12 = mth
c     
      ALLOCATE ( xpb(mth),zpb(mth),cwrkr(mth),cwrki(mth) )
c     
      ns = mth
c     
c.....get array of observer points, and initialize chi matrix.
c     
      ALLOCATE ( zchc(mth+1,jmax1), zchs(mth+1,jmax1) )
c
c     
      DO i = 1, mth
c     
         chir(i) = 0.0
         chii(i) = 0.0
         cwrkr(i) = 0.0
         cwrki(i) = 0.0
         IF (i.LT.mth) THEN
            xpb(i) = (xpla(i)+xpla(i+1))/2.
            zpb(i) = (zpla(i)+zpla(i+1))/2.
         ELSE IF (i.EQ.mth) THEN
            xpb(mth)=(xpla(mth)+xpla(1))/2.
            zpb(mth)=(zpla(mth)+zpla(1))/2.
         END IF
      END DO                    ! Loop on I, mth
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
     $     isg
      WRITE ( iotty,'("isgn for plasma source contrib. = ", i3)')
     $     isg
c     
      CALL chiplp1(xpla,zpla,xplap,zplap,isg, zchc,zchs,jmax1,mth,
     $     ns,1,cwrkr,cwrki,xpb,zpb, blrp,blip )     
      CALL atpoint("After Chiplp1","mth",mth,"xobs",xpb(1),
     $     iotty, outmod )
      
c..   Store away the plasma result.
c     
      DO  i = 1, mth
         chir(i) = cwrkr(i)
         chii(i) = cwrki(i)
      END DO
c
c
      write (outmod,'("Plasma contribution to plasmachi")')
      DO  i = mth/2, mth/2+2
         write(outmod,'("i,chir,chii,cwrkr,cwrki",i4,1p4e12.3)')
     $        i,chir(i),chii(i),cwrkr(i),cwrki(i)
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
      DO i = 1, mth
         cwrkr(i) = 0.0
         cwrki(i) = 0.0
      END DO
c     
      isg = -1 
      WRITE ( outmod,
     $     '("isgn for shell integ. contrib. = ", i3)') isg
c     
      CALL chiplp1(xwal,zwal,xwalp,zwalp,isg,zchc,zchs,jmax1,mth,
     $     ns,0,cwrkr,cwrki,xpb,zpb,blrp,blip )
c     
      DO  i = 1, mth
         chir(i) = chir(i) + cwrkr(i)
         chii(i) = chii(i) + cwrki(i)
      END DO
C     
c
      write (outmod,'("Plasma contribution to plasmachi")') 
      DO  i = mth/2, mth/2+2
         write(outmod,'("i,chir,chii,cwrkr,cwrki",i4,1p4e12.3)')
     $        i,chir(i),chii(i),cwrkr(i),cwrki(i)
      END DO
c
 200  CONTINUE
c
c     Integrate Over Shell First.
c     
      zchc(1:mw,1:jwal1) =
     $     grri(mth+1:mth+mw,2*jmax1+1:2*jmax1+jwal1)    
      zchs(1:mw,1:jwal1) =
     $     grri(mth+1:mth+mw,2*jmax1+jwal1+1:2*jmax1+jwal2)
c     
      isg = 1
      WRITE ( outmod,'("isgn for shell source contrib. = ", i3)')
     $     isg
      WRITE ( iotty,'("isgn for shell source contrib. = ", i3)')
     $     isg
c     
      DO i = 1, mth
         cwrkr(i) = 0.0
         cwrki(i) = 0.0
      END DO
c
      CALL chiplw1(xwal,zwal,xwalp,zwalp,isg, zchc,zchs,jwal1,mth,
     $     ns,1,cwrkr,cwrki,xpb,zpb,blrw,bliw )     
c     
      CALL atpoint("After Chiplw1","ns",ns,"zobs",zpb(1),
     $     iotty, outmod )
      
c..   Store Away the Shell Result.
c     
      DO  i = 1, mth
         chir(i) = chir(i)+cwrkr(i)
         chii(i) = chii(i)+cwrki(i)
      END DO
c
c
      write (outmod,'("Wall contribution to plasmachi")')
      DO  i = mth/2, mth/2+2
         write(outmod,'("i,chir,chii,cwrkr,cwrki",i4,1p4e12.3)')
     $        i,chir(i),chii(i),cwrkr(i),cwrki(i)
      END DO
c
c     
c.... Integrate over Plasma:
c     
      zchc(1:mth,1:jwal1) =
     $     grri(1:mth,2*jmax1+1:2*jmax1+jwal1)
      zchs(1:mth,1:jwal1) =
     $     grri(1:mth,2*jmax1+jwal1+1:2*jmax1+jwal2)
C     
      DO i = 1, mth
         cwrkr(i) = 0.0
         cwrki(i) = 0.0
      END DO
c     
      isg = -1 
      WRITE ( outmod,
     $     '("isgn for shell integ. contrib. = ", i3)') isg
c     
      CALL chiplw1(xpla,zpla,xplap,zplap,isg,zchc,zchs,jwal1,mth,
     $     ns,0,cwrkr,cwrki,xpb,zpb,blrw,bliw )
c     
      DO  i = 1, mth
         chir(i) = chir(i) + cwrkr(i)
         chii(i) = chii(i) + cwrki(i)
      END DO
c
      write (outmod,'("Wall contribution to plasmachi")')
      DO  i = mth/2, mth/2+2
         write(outmod,'("i,chir,chii,cwrkr,cwrki",i4,1p4e12.3)')
     $        i,chir(i),chii(i),cwrkr(i),cwrki(i)
      END DO
c     
      chir(mth+1)=chir(1)
      chii(mth+1)=chii(1)
      DEALLOCATE ( xpb,zpb, cwrkr,cwrki,zchc,zchs )
c
      RETURN
c
      END SUBROUTINE chiplasma

c.............................................................
      SUBROUTINE chiplp1(xsce,zsce,xscp,zscp,isg,creal,cimag,ndimc,
     $     nobs,ns,ip,chir,chii,xobs,zobs,blr,bli)
c......................................................................
c     
      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum5.inc'
c     
      DIMENSION blr(*),bli(*), xsce(*),zsce(*),xscp(*),zscp(*)
      DIMENSION xobs(*),zobs(*)
      DIMENSION creal(ns+1,*), cimag(ns+1,*)
      DIMENSION chir(*), chii(*)
c     
      REAL nq
c     
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
c     
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

c     Divide BVAL by twopi to account for definition of script G to be 
c     onsistent with K.
c     
            bval = bval / factpi
c     
            DO l1 = 1, jmax1
c     
               zbr = blr(l1)
               zbi = bli(l1)
               chir(io) = chir(io) + 
     $              aval * ( creal(is,l1)*zbr - cimag(is,l1)*zbi )
               chii(io) = chii(io) +
     $              aval * ( cimag(is,l1)*zbr + creal(is,l1)*zbi )
c     
               IF ( ip .EQ. 0 ) GO TO 60
c     
c##########change plus to minus bval  #########
               ipm = + 1.0
               chir(io) = chir(io) + ipm * bval
     $              * ( cslth(is,l1)*zbr - snlth(is,l1)*zbi )
               chii(io) = chii(io) + ipm * bval
     $              * ( snlth(is,l1)*zbr + cslth(is,l1)*zbi )
c     
 60            CONTINUE
            END DO              ! Loop over l1: 1, jmax1
         END DO                 ! Loop over IS: 1, nobs
c     
         chir(io) = 0.5 * isg*dtpw * chir(io)
         chii(io) = 0.5 * isg*dtpw * chii(io)
c     
      END DO                    ! Loop over Io: 1, nobs
c     
      RETURN
      END SUBROUTINE chiplp1

c.....................................................................
      SUBROUTINE chiplw1(xsce,zsce,xscp,zscp,isg,creal,cimag,ndimc,
     $     nobs,ns,ip,chir,chii,xobs,zobs,blr,bli)
c......................................................................
c     
      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum3.inc'
      INCLUDE 'vacuum2.inc'
      INCLUDE 'vacuum5.inc'
      INCLUDE 'vacuum7.inc'
c     
      DIMENSION blr(*),bli(*), xsce(*),zsce(*),xscp(*),zscp(*)
      DIMENSION xobs(*),zobs(*)
      DIMENSION creal(ns+1,*), cimag(ns+1,*)
      DIMENSION chir(*), chii(*)
c     
      REAL nq
c     
c     
c     Factpi for dividing BVAL to be consistent with AVAL.
c     
      factpi = twopi
      q = qa1
c     
      nq = n * q
      dtpw = twopi / ns
c     
c     
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

c     Divide BVAL by twopi to account for definition of script G to be 
c     onsistent with K.
c     
            bval = bval / factpi
c     
            DO l1 = 1, jwal1
c     
               zbr = blr(l1)
               zbi = bli(l1)
c               chir(io) = chir(io) + 
               chir(io) = chir(io) - 
c     $              aval * ( creal(is,l1)*zbr + cimag(is,l1)*zbi )/2
c     $              aval * ( creal(is,l1)*zbr - cimag(is,l1)*zbi )/2
     $              aval * ( creal(is,l1)*zbr ) /2.
c               chii(io) = chii(io) +
               chii(io) = chii(io) -
c     $              aval * ( cimag(is,l1)*zbr - creal(is,l1)*zbi )/2
c     $              aval * ( cimag(is,l1)*zbr + creal(is,l1)*zbi )/2
     $              aval * ( cimag(is,l1)*zbi ) /2.
c     
               IF ( ip .EQ. 0 ) GO TO 60
c     
c##########change plus to minus bval  #########
               ipm = + 1.0
c     it turned out that a should change sign on the wall
c               ipm= - 1.0 
               chir(io) = chir(io) + ipm * bval
     $              * rwvce(is,l1)*zbr /2.
               chii(io) = chii(io) + ipm * bval
     $              * rwvco(is,l1)*zbi /2.

c$$$               chir(io) = chir(io) + ipm * bval
c$$$     $              * ( cslth(is,l1)*zbr - snlth(is,l1)*zbi )
c$$$               chii(io) = chii(io) + ipm * bval
c$$$     $              * ( snlth(is,l1)*zbr + cslth(is,l1)*zbi )
c     
 60            CONTINUE
            END DO              ! Loop over l1: 1, jwal1
         END DO                 ! Loop over IS: 1, nobs
c     
         chir(io) = 0.5 * isg*dtpw * chir(io)
         chii(io) = 0.5 * isg*dtpw * chii(io)
c     
      END DO                    ! Loop over Io: 1, nobs
c     
      RETURN
      END SUBROUTINE chiplw1
