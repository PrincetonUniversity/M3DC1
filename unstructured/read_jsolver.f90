module jsolver
  
  implicit none
  
  real, allocatable :: p_jsv(:), ppxx_jsv(:), q_jsv(:), qp_jsv(:)
  real, allocatable :: gxx_jsv(:), gpx_jsv(:), fb_jsv(:), fbp_jsv(:)
  real, allocatable :: f_jsv(:), fp_jsv(:), psival_jsv(:)
  real, allocatable :: x_jsv(:,:), z_jsv(:,:), aj3_jsv(:,:), aj_jsv(:,:)
  real, allocatable :: te_jsv(:), ane_jsv(:), ti_jsv(:)
  real, allocatable :: avez_jsv(:), zeffa_jsv(:)
 
  integer :: nz1, nz2, nthe, npsi_jsv, isym, ibcxmin, ibcxmax
  real :: dr, dt, xzero_jsv, p0_jsv, plimst, pminst, bcxmin, bcxmax
  
contains
 
  subroutine load_jsolver
 
    implicit none
 
    integer, parameter :: ifile = 33
    integer :: i,j
 
    open(unit=ifile,file='fixed',status='old')
 
    read(ifile,1098) nz1,nz2,nthe,npsi_jsv,isym,dr,dt
    read(ifile,1099) xzero_jsv,p0_jsv,plimst,pminst
 
  1098 format(5i5,1p2e12.4)
  1099 format(1p5e16.8)
 
    allocate(p_jsv(nz2), ppxx_jsv(nz2), q_jsv(nz2), qp_jsv(nz2), &
         gxx_jsv(nz2), gpx_jsv(nz2), fb_jsv(nz2), fbp_jsv(nz2), &
         f_jsv(nz2), fp_jsv(nz2), psival_jsv(nz2), &
         x_jsv(nz1,nz2), z_jsv(nz1,nz2), aj3_jsv(nz1,nz2), aj_jsv(nz1,nz2), &
         te_jsv(nz2), ane_jsv(nz2), ti_jsv(nz2), &
         avez_jsv(nz2), zeffa_jsv(nz2))
 
    read(ifile,1099) (p_jsv(j),j=1,nz2)
    read(ifile,1099) (ppxx_jsv(j),j=1,nz2)
    read(ifile,1099) (q_jsv(j),j=1,nz2)
    read(ifile,1099) (qp_jsv(j),j=1,nz2)
    read(ifile,1099) (gxx_jsv(j),j=1,nz2)
    read(ifile,1099) (gpx_jsv(j),j=1,nz2)
    read(ifile,1099) (fb_jsv(j),j=1,nz2)
    read(ifile,1099) (fbp_jsv(j),j=1,nz2)
    read(ifile,1099) (f_jsv(j),j=1,nz2)
    read(ifile,1099) (fp_jsv(j),j=1,nz2)
    read(ifile,1099) (psival_jsv(j),j=1,nz2)
    read(ifile,1099) ((x_jsv(i,j),i=1,nz1),j=1,nz2)
    read(ifile,1099) ((z_jsv(i,j),i=1,nz1),j=1,nz2)
    read(ifile,1099) ((aj3_jsv(i,j),i=1,nz1),j=1,nz2)
    read(ifile,1099) ((aj_jsv(i,j),i=1,nz1),j=1,nz2)
    read(ifile,1099) (te_jsv(j),j=1,nz2)
    read(ifile,1099) (ane_jsv(j),j=1,nz2)
    read(ifile,1099) (ti_jsv(j),j=1,nz2)
    read(ifile,1099) (avez_jsv(j),j=1,nz2)
    read(ifile,1099) (zeffa_jsv(j),j=1,nz2)
  !
  !.....note on the centering of the variables:
  !  defined at flux value psival_jsv(j):
  !  ppxx_jsv(j), gpx_jsv(j)
  !
  !  defined at flux value .5*(psival_jsv(j)+psival_jsv(j+1)):
  !  p_jsv(j), gxx_jsv(j), f(j)
 
    close(ifile)
 
  end subroutine load_jsolver
 
 
  subroutine unload_jsolver
    implicit none
 
    deallocate(p_jsv, ppxx_jsv, q_jsv, qp_jsv, &
         gxx_jsv, gpx_jsv, fb_jsv, fbp_jsv, &
         f_jsv, fp_jsv, psival_jsv, &
         x_jsv, z_jsv, aj3_jsv, aj_jsv, &
         te_jsv, ane_jsv, ti_jsv, &
         avez_jsv, zeffa_jsv)
 
  end subroutine unload_jsolver
 
 
end module jsolver
 
 
#ifdef READ_JSOLVER
  
Program readjsolver
  use jsolver
  use polar

  implicit none

  REAL*8, ALLOCATABLE :: psinormt(:), g4big0t(:), g4bigt(:), g4bigp(:),g4bigpp(:)
  REAL*8, ALLOCATABLE :: fbig0t(:), fbigt(:), fbigp(:), fbigpp(:)
  REAL*8, ALLOCATABLE :: psinorm(:), g4big0(:), g4big(:), fbig0(:), fbig(:), ffpnorm(:)
  REAL*8, allocatable :: fspl(:,:)

  real :: delpsi, xzero, bzero, denom, gzero, g1, g4bigim1, g4save, fval(3)
  real :: pzero
  integer :: ilinx, ier, ict(3), i40, i50, i, j, l

  call load_jsolver
  
  call write_polar(x_jsv(3:nthe+2,1:npsi_jsv), z_jsv(3:nthe+2,1:npsi_jsv), &
       nthe, npsi_jsv, 0)
  
  allocate(ffpnorm(npsi_jsv))
!
!
!
  xzero = xzero_jsv
  bzero = gxx_jsv(npsi_jsv)

  ! remove factor of 1/xzero_jsv in definition of g
  gxx_jsv = gxx_jsv*xzero_jsv
  gpx_jsv = gpx_jsv*xzero_jsv
  pzero = 1.5*p_jsv(1)   - .5*p_jsv(2)
  gzero = 1.5*gxx_jsv(1) - .5*gxx_jsv(2)
! shift gxx and p to be defined at psi_j locations
! instead of at psi_{j+1/2} locations
  do l=npsi_jsv,2,-1
     p_jsv(l) = (p_jsv(l-1) + p_jsv(l))/2.
     gxx_jsv(l) = (gxx_jsv(l-1) + gxx_jsv(l))/2.
  end do
  p_jsv(1)   = pzero
  gxx_jsv(1) = gzero
!......define ffpnorm    ff'
  ffpnorm = gpx_jsv*gxx_jsv
!
!
!
  write(*,2000) xzero, bzero, gzero, pzero,x_jsv(3,1),z_jsv(3,1)
2000 format( "xzero = ",1pe12.4,"   bzero = ",1pe12.4,"   gzero = ",1pe12.4,"   pzero ",e12.4,/," xmag,zmag = ",1p2e12.4)
  ! write profiles
  open(unit=75,file="profiles",status="unknown")
  write(75,801) npsi_jsv
  do j=1,npsi_jsv
     write(75,812) j,psival_jsv(j),gxx_jsv(j),ffpnorm(j),p_jsv(j),ppxx_jsv(j)
  enddo
  close(75)
!
!....prepare arrays needed to compute p and g functions in M3D-C1
  allocate(psinormt(npsi_jsv),g4big0t(npsi_jsv),g4bigt(npsi_jsv),g4bigp(npsi_jsv),g4bigpp(npsi_jsv))
  allocate(fbig0t(npsi_jsv),fbigt(npsi_jsv),fbigp(npsi_jsv),fbigpp(npsi_jsv))
  allocate(psinorm(npsi_jsv),g4big0(npsi_jsv),g4big(npsi_jsv),fbig0(npsi_jsv),fbig(npsi_jsv))
  allocate(fspl(0:1,npsi_jsv))
!
  delpsi = psival_jsv(npsi_jsv) - psival_jsv(1)
  do j=1,npsi_jsv
     psinormt(j) = (psival_jsv(j)-psival_jsv(1))/(delpsi)
     g4big0t(j) = 0.5*(gxx_jsv(j)**2 - (gxx_jsv(npsi_jsv))**2)
     g4bigt(j)= ffpnorm(j)
     fbig0t(j) = p_jsv(j)     !NOTE:  now fbig0t is total pressure
     fbigt(j) = ppxx_jsv(j)
  enddo
  !
  do j=1,npsi_jsv
     fspl(0,j) = g4bigt(j)
  enddo
  ibcxmin = 4
  bcxmin = 0.
  ibcxmax = 0
  bcxmax = 0.
  call R8mkspline(psinormt,npsi_jsv,fspl,ibcxmin,bcxmin,ibcxmax,bcxmax,ilinx,ier)
  if(ier.ne.0) write(*,5001) ier
5001 format("ier= ",i3,"  after call to R8mkspline")
  ict(1) = 1
  ict(2) = 1
  ict(3) = 1
  do j=1,npsi_jsv
     call r8evspline(psinormt(j),psinormt,npsi_jsv,ilinx,fspl,ict,fval,ier)
     if(ier.ne.0) write(*,5002) ier
5002 format("ier=",i3,"  after call to r8evspline")
     g4big(j) = gxx_jsv(j)*gpx_jsv(j)
     g4bigp(j) = fval(1)
     g4bigpp(j) = fval(2)
  enddo
  !
  !.....smooth derivatives near origin
  i40 = npsi_jsv/6
  i50 = npsi_jsv/5
  denom = psinormt(i50)-psinormt(i40)
  gzero = (g4bigp(i40)*psinormt(i50)-g4bigp(i50)*psinormt(i40))/denom
  g1 = (g4bigp(i50)-g4bigp(i40))/denom
  do i=1,i40
     ! g4bigp(i) = gzero + g1*psinormt(i)
     ! g4bigpp(i) = psinormt(i)*g4bigpp(i40)/psinormt(i40)
  enddo
  !
  !.....additional smoothing for second derivative
  g4bigim1 = g4bigpp(1)
  do i=2,npsi_jsv-1
     g4save = g4bigpp(i)
    ! g4bigpp(i) = .5*g4bigpp(i) + .25*(g4bigpp(i+1)+g4bigim1)
     g4bigim1 = g4save
  enddo
  !
  do j=1,npsi_jsv
     fspl(0,j) = fbigt(j)
  enddo
  ibcxmin = 0.
  bcxmin = 0.
  ibcxmax = 0
  bcxmax = 0.
  call R8mkspline(psinormt,npsi_jsv,fspl,ibcxmin,bcxmin,ibcxmax,bcxmax,ilinx,ier)
  if(ier.ne.0) write(*,5001) ier
  ict(1) = 1
  ict(2) = 1
  ict(3) = 1
  do j=1,npsi_jsv
     call r8evspline(psinormt(j),psinormt,npsi_jsv,ilinx,fspl,ict,fval,ier)
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
! !.....additional smoothing for second derivative
! g4bigim1 = fbigpp(1)
! do i=2,npsi_jsv-1
!    g4save = fbigpp(i)
!    fbigpp(i) = .5*fbigpp(i) + .25*(fbigpp(i+1)+g4bigim1)
!    g4bigim1 = g4save
! enddo
  !
  open(unit=76,file="profiles-p",status="unknown")
  write(76,803) npsi_jsv
  do j=1,npsi_jsv
     write(76,802)  psinormt(j),fbig0t(j),fbigt(j),fbigp(j),fbigpp(j)
  enddo
  close(76)
  !
  open(unit=77,file="profiles-g",status="unknown")
  write(77,804) npsi_jsv
  do j=1,npsi_jsv
     write(77,802)  psinormt(j),g4big0t(j),g4big(j),gxx_jsv(j), gpx_jsv(j)
  enddo
802 format(5x,6e18.9)
812 format(i5,1p6e18.10)
  close(77)
  !
  
  deallocate(ffpnorm)
  call unload_jsolver
  stop 0
  
51 format(" Rzero, Btor, I_P = ", 1p3e16.8)
100 format(2i5)
200 format(4e19.12)
801 format(i5,"         psi                f                ff'",     &
         "                p                 p'")
803 format(i5,"      psinorm            fbig0             fbig",   &
         "            fbigp             fbigpp")
804 format(i5,"      psinorm           g4big0            g4big")


  
end Program readjsolver

#endif

