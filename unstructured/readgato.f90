      Program readgato
!=======================================================================
!     read in equilibrium quantities from file dskgato
!     write files needed by M3D-C1
!     scj    10/25/08
!=======================================================================
!
!--> Part 0.0  declare variables
      IMPLICIT  none
      REAL*8, ALLOCATABLE  :: xnorm(:,:), znorm(:,:)
      REAL*8, ALLOCATABLE  :: psiflux(:), fnorm(:), ffpnorm(:), ponly(:)
      REAL*8, ALLOCATABLE  :: pponly(:), qsf(:), dnorm(:), dpdz(:), dpdr(:)
      REAL*8, ALLOCATABLE  :: sn(:),xn(:,:),zn(:,:),snorm(:)
      REAL*8, ALLOCATABLE  :: snorm2(:), xnew(:,:), znew(:,:)
      REAL*8 :: d0, xmag, zmag,btorgato,rcgato, eaxe, curtot, ds, fac, dela, delb, denom, delpsi, amu0, pi
      REAL*8, ALLOCATABLE :: psinormt(:), g4big0t(:), g4bigt(:), g4bigp(:),g4bigpp(:)
      REAL*8, ALLOCATABLE :: fbig0t(:), fbigt(:), fbigp(:), fbigpp(:)
      REAL*8, ALLOCATABLE :: psinorm(:), g4big0(:), g4big(:), fbig0(:), fbig(:)
      INTEGER :: i,j,ngato, nthet, npsi, nt, np, jjold, jstart, jj
      INTEGER :: ii, im, iz, iiold, istart, nthe, iprint, isym
      REAL*8, allocatable :: fspl(:,:)
      REAL*8 :: bcxmin,bcxmax,fval(3),gzero,g1,g2,psi30,psi40,g30,g40,g4bigim1,g4save
      integer :: ict(3),ibcxmin,ibcxmax,ilinx,ier,i40,i50

      real, allocatable :: norm(:,:), curv(:)
!
      pi = acos(-1.)
      amu0 = pi*4.e-7
!     define number of output points
      np = 51
      nt = 50
      iprint = 0
     
!--> Part 1.0 read data from file dskgato
!     open output file
      ngato=73
      open(unit=ngato,file="dskgato",status="old")
!   
!.....symmetry of data file
      isym=1
     
      read(ngato,100) npsi,nthet  ! number of mesh points
      read(ngato,200) rcgato, xmag, zmag, btorgato
      read(ngato,200) curtot,eaxe,d0
!
!     NOTE:
!     for (isym.eq.0) ...no symmetry:  points 1 and nthe+1 are the same
!         (isym.eq.1) ......symmetry:  points 1 and nthe/2 + 1 are on midplane
!
      if(isym.eq.0) then
        nthe = nthet
      else
        nthe = 2*(nthet-1)
      endif
      write(*,50) isym,npsi, nthet, nthe
      write(*,51) rcgato, btorgato, curtot
!
      allocate(psiflux(npsi))
      allocate(fnorm(npsi))
      allocate(ffpnorm(npsi))
      allocate(ponly(npsi))
      allocate(pponly(npsi))
      allocate(qsf(npsi))
      allocate(dnorm(npsi))
      allocate(dpdz(nthe))
      allocate(dpdr(nthe))
      allocate(xnorm(npsi,nthe))
      allocate(znorm(npsi,nthe))
      allocate(snorm(npsi),sn(np))
      allocate(xn(np,nthe),zn(np,nthe),snorm2(nthe+1),xnew(np,nt),znew(np,nt))
      allocate(norm(2,nthe))
      allocate(curv(nthe))
!
      read(ngato,200) (psiflux(j),j=1,npsi) ! poloidal flux
      read(ngato,200) (fnorm(j),j=1,npsi) ! f
      read(ngato,200) (ffpnorm(j),j=1,npsi) ! ff'
      read(ngato,200) (ponly(j),j=1,npsi) ! mu0 p
      read(ngato,200) (pponly(j),j=1,npsi) ! mu0 p'
      read(ngato,200) (qsf(j),j=1,npsi) ! safety factor
!
      write(*,300) qsf(1),qsf(npsi)
!
      read(ngato,200) (dnorm(j),j=1,npsi) ! density set = dnorm(i)
      write(*,301) dnorm(1),dnorm(npsi)

      write(*,303) amu0*curtot,  amu0*ponly(1)
!
      read(ngato,200) (dpdz(i),i=1,nthet) ! dpsi/dz on boundary
      read(ngato,200) (dpdr(i),i=1,nthet) ! dpsi/dr on boundary
!
      if(iprint.eq.1) write(*,400)
      read(ngato,200) ((xnorm(j,i),j=1,npsi),i=1,nthet) ! R
!
      write(*,500) xnorm(1,1), xmag
      read(ngato,200) ((znorm(j,i),j=1,npsi),i=1,nthet) ! Z
!
      write(*,600) znorm(1,1), zmag
!
!     check on the periodicity conditions
      if(iprint.eq.1) then 
        write(*,602)
        do j=1,npsi,10
          write(*,601) j,xnorm(j,nthet-1),xnorm(j,nthet),xnorm(j,1),xnorm(j,2)
          write(*,601) j,znorm(j,nthet-1),znorm(j,nthet),znorm(j,1),znorm(j,2)
        enddo
      endif
!
      close(ngato)
!
!.....reflect points if isym.eq.1
      if(isym.eq.1) then
        do j=1,npsi
          do i=nthet+1,nthe
            xnorm(j,i) = xnorm(j,2*nthet-i)
            znorm(j,i) =-znorm(j,2*nthet-i)
          enddo
        enddo
      endif
!
!--> Part 2.0 Interpolate onto finer grid
!
!     2.1 for each i, interpolate in j as an intermediate step
      do i=1,nthe
        snorm(1) = 0.
        do j=2,npsi
          snorm(j) = snorm(j-1) + sqrt((xnorm(j,i)-xnorm(j-1,i))**2    &
                                      +(znorm(j,i)-znorm(j-1,i))**2)
        enddo
        ds = snorm(npsi)/(np-1)
        sn(1) = 0.
        xn(1,i) = xnorm(1,i)
        zn(1,i) = znorm(1,i)
        jjold = 2
        do j=2,np
          sn(j) = (j-1)*ds
          jstart = jjold
          do jj = jstart,npsi
            jjold = jj
            if(snorm(jj).ge.sn(j)) exit
          enddo
!
!       new point is between old points jj-1 and jj
          fac = (sn(j)-snorm(jj-1))/(snorm(jj)-snorm(jj-1))
          xn(j,i) = xnorm(jj-1,i) + fac*(xnorm(jj,i)-xnorm(jj-1,i))
          zn(j,i) = znorm(jj-1,i) + fac*(znorm(jj,i)-znorm(jj-1,i))
        enddo
      enddo
      if(iprint.eq.1) then
        write(*,402) 
        do j=1,np
          do i=1,nthe,10
             write(*,401) i,j,xn(j,i),zn(j,i)
          enddo
        enddo
      endif
!
!     2.2 for each new j, interpolate in i
      do i=1,nt
        xnew(1,i) = xn(1,1)
        znew(1,i) = zn(1,1)
      enddo
      do j=2,np
        snorm2(1) = 0.
        do ii=2,nthe+1
          i = ii
          im = ii-1
          if (ii.eq.nthe+1) i = 1
          snorm2(ii) = snorm2(im) + sqrt((xn(j,i)-xn(j,im))**2    &
                                      +(zn(j,i)-zn(j,im))**2)
        enddo
        ds = snorm2(nthe+1)/(nt)
        sn(1) = 0.
        xnew(j,1) = xn(j,1)
        znew(j,1) = zn(j,1)
        iiold = 2
        do i=2,nt
          sn(i) = (i-1)*ds
          istart = iiold
          do ii = istart,nthe+1

            iiold = ii
            if(snorm2(ii).ge.sn(i)) exit
          enddo
          im = ii-1
          iz = ii
          if(ii.eq.nthe+1) iz = 1
!
!         new point is between old points ii-1 and ii
          fac = (sn(i)-snorm2(im))/(snorm2(ii)-snorm2(im))
          xnew(j,i) = xn(j,im) + fac*(xn(j,iz)-xn(j,im))
          znew(j,i) = zn(j,im) + fac*(zn(j,iz)-zn(j,im))
        enddo
      enddo
      if(iprint.eq.1) then
        do j=1,np
          do i=1,nt
             write(*,401) i,j,xnew(j,i),znew(j,i)
          enddo
        enddo
      endif
!
      ! write POLAR file

      open(unit=74,file="POLAR",status="unknown")
      write(74,701) np,nt
      do j=1,np-1
        do i=1,nt
          write(74,702) xnew(j,i),znew(j,i)
        enddo
      enddo
      call calc_norms(xnew(np,:),znew(np,:),nt,norm,curv)
      do i=1,nt
         write(74,703) xnew(np,i),znew(np,i),norm(1,i),norm(2,i),curv(i)
      end do

      close(74)
!
      open(unit=75,file="profiles",status="unknown")
      write(75,801) npsi
      do j=1,npsi
        write(75,802) j,psiflux(j),fnorm(j),ffpnorm(j),ponly(j),pponly(j), qsf(j)
      enddo
      close(75)
!
!....prepare arrays needed to compute p and g functions in M3D-C1
      allocate(psinormt(npsi),g4big0t(npsi),g4bigt(npsi),g4bigp(npsi),g4bigpp(npsi))
      allocate(fbig0t(npsi),fbigt(npsi),fbigp(npsi),fbigpp(npsi))
      allocate(psinorm(npsi),g4big0(npsi),g4big(npsi),fbig0(npsi),fbig(npsi))
      allocate(fspl(0:1,npsi))
!
      delpsi = psiflux(npsi) - psiflux(1)
      do j=1,npsi
        psinormt(j) = (psiflux(j)-psiflux(1))/(delpsi)
        g4big0t(j) = 0.5*(fnorm(j)**2 - (rcgato*btorgato)**2)
        g4bigt(j)= ffpnorm(j)
        fbig0t(j) = ponly(j)/ponly(1)
        fbigt(j) = pponly(j)*amu0
      enddo
!
      do j=1,npsi
        fspl(0,j) = g4bigt(j)
      enddo
      ibcxmin = 4
      bcxmin = 0.
      ibcxmax = 0
      bcxmax = 0.
      call R8mkspline(psinormt,npsi,fspl,ibcxmin,bcxmin,ibcxmax,bcxmax,ilinx,ier)
      if(ier.ne.0) write(*,5001) ier
 5001 format("ier= ",i3,"  after call to R8mkspline")
     ict(1) = 1
     ict(2) = 1
     ict(3) = 1
     do j=1,npsi
     call r8evspline(psinormt(j),psinormt,npsi,ilinx,fspl,ict,fval,ier)
     if(ier.ne.0) write(*,5002) ier
 5002 format("ier=",i3,"  after call to r8evspline")
     g4big(j) = fval(1)
     g4bigp(j) = fval(2)
     g4bigpp(j) = fval(3)
     enddo
!
!.....smooth derivatives near origin
      i40 = 40
      i50 = 50
      denom = psinormt(i50)-psinormt(i40)
      gzero = (g4bigp(i40)*psinormt(i50)-g4bigp(i50)*psinormt(i40))/denom
      g1 = (g4bigp(i50)-g4bigp(i40))/denom
      do i=1,i40
      g4bigp(i) = gzero + g1*psinormt(i)
      g4bigpp(i) = psinormt(i)*g4bigpp(i40)/psinormt(i40)
      enddo
!
!.....additional smoothing for second derivative
      g4bigim1 = g4bigpp(1)
      do i=2,npsi-1
      g4save = g4bigpp(i)
      g4bigpp(i) = .5*g4bigpp(i) + .25*(g4bigpp(i+1)+g4bigim1)
      g4bigim1 = g4save
      enddo
!
      do j=1,npsi
        fspl(0,j) = fbigt(j)
      enddo
      ibcxmin = 0.
      bcxmin = 0.
      ibcxmax = 0
      bcxmax = 0.
      call R8mkspline(psinormt,npsi,fspl,ibcxmin,bcxmin,ibcxmax,bcxmax,ilinx,ier)
      if(ier.ne.0) write(*,5001) ier
     ict(1) = 1
     ict(2) = 1
     ict(3) = 1
     do j=1,npsi
     call r8evspline(psinormt(j),psinormt,npsi,ilinx,fspl,ict,fval,ier)
     if(ier.ne.0) write(*,5002) ier
     fbig(j) = fval(1)
     fbigp(j) = fval(2)
     fbigpp(j) = fval(3)
     enddo
!
!.....smooth derivatives near origin
      denom = psinormt(i50)-psinormt(i40)
      gzero = (fbigp(i40)*psinormt(i50)-fbigp(i50)*psinormt(i40))/denom
      g1 = (fbigp(i50)-fbigp(i40))/denom
      do i=1,i40
      fbigp(i) = gzero + g1*psinormt(i)
      fbigpp(i) = psinormt(i)*fbigpp(i40)/psinormt(i40)
      enddo
!
!.....additional smoothing for second derivative
      g4bigim1 = fbigpp(1)
      do i=2,npsi-1
      g4save = fbigpp(i)
      fbigpp(i) = .5*fbigpp(i) + .25*(fbigpp(i+1)+g4bigim1)
      g4bigim1 = g4save
      enddo
!
      open(unit=76,file="profiles-p",status="unknown")
      write(76,803) npsi
      do j=1,npsi
        write(76,802) j,psinormt(j),fbig0t(j),fbig(j),fbigp(j),fbigpp(j)
      enddo
      close(76)
!
      open(unit=77,file="profiles-g",status="unknown")
      write(77,804) npsi
      do j=1,npsi
        write(77,802) j,psinormt(j),g4big0t(j),g4big(j),g4bigp(j),g4bigpp(j)
      enddo
      close(77)
!

      deallocate(psiflux)
      deallocate(fnorm)
      deallocate(ffpnorm)
      deallocate(ponly)
      deallocate(pponly)
      deallocate(qsf)
      deallocate(dnorm)
      deallocate(dpdz)
      deallocate(dpdr)
      deallocate(xnorm)
      deallocate(znorm)
      deallocate(snorm,sn)
      deallocate(xn,zn,snorm2,xnew,znew)
      deallocate(norm,curv)



      stop 0
   50 format(" isym, npsi, nthet, nthe  = ",4i5)
   51 format(" Rzero, Btor, I_P = ", 1p3e16.8)
  100 format(2i5)
  200 format(4e19.12)
  400 format(" read xnorm and znorm")
  301 format(" dnorm(1) and dnorm(npsi)",1p2e12.4)
  300 format(" q(0) and q(edge) ",1p2e12.4)
  303 format("  normalized current and central pressure", 1p2e12.4)
  401        format(2i3,1p2e12.4)
  402   format(" intermediate step")
  500 format(" x(1,1), xmag  = ", 1p2e12.4)
  600 format(" z(1,1), zmag  = ", 1p2e12.4)
  601   format(i4,1p4e12.4)
  602   format(" periodicty check")
  701 format(2i5)
  702 format(1p2e18.10)
  703 format(1p5e18.10)
  801 format(i5,"         psi                f                ff'",     &
                 "                p                 p'")
  803 format(i5,"      psinorm            fbig0             fbig",   &
                "            fbigp             fbigpp")
  804 format(i5,"      psinorm           g4big0            g4big",   &
                "           g4bigp            g4bigpp")
  802 format(i5,1p6e18.10)
      end


subroutine calc_norms(x, z, n, norm, curv)

  implicit none

  integer, intent(in) :: n
  real*8, intent(in), dimension(n) :: x, z
  real, intent(out), dimension(2,n) :: norm
  real, intent(out), dimension(n) :: curv

  integer :: i, i1, i2
  real :: dx1, dx2, dz1, dz2, l, t1, t2, t3
  real :: l1(n), l2(n), norm1(2), norm2(2)

  ! calculate normals
  do i=1,n
     i1 = mod(i+n-2,n)+1
     i2 = mod(i,n)+1
     dx1 = x(i) - x(i1)
     dx2 = x(i2) - x(i)
     dz1 = z(i) - z(i1)
     dz2 = z(i2) - z(i)

     l1(i) = 1./sqrt(dx1**2 + dz1**2)
     l2(i) = 1./sqrt(dx2**2 + dz2**2)
     norm1(1) =  dz1*l1(i)
     norm1(2) = -dx1*l1(i)
     norm1(1) =  dz2*l2(i)
     norm1(2) = -dx2*l2(i)

     ! perform weigted average of adjacent edge normals
     norm(:,i) = (l1(i)*norm1 + l2(i)*norm2)/(l1(i)+l2(i))

     ! normalize normal
     l = sqrt(norm(1,i)**2 + norm(2,i)**2)
     norm(:,i) = norm(:,i)/l
  enddo

  ! calculate curvature
  do i=1,n
     i1 = mod(i+n-2,n)+1
     i2 = mod(i,n)+1

     t1 = atan2(norm(2,i1),norm(1,i1))
     t2 = atan2(norm(2,i ),norm(1,i ))
     t3 = atan2(norm(2,i2),norm(1,i2))
     
     ! perform weigted average of dtheta/dl
     curv(i) = ((t2-t1)*l1(i)**2 + (t3-t2)*l2(i)**2)/(l1(i) + l2(i))
  end do
  

end subroutine calc_norms
