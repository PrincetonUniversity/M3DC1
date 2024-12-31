! This module reads in and processes data from VMEC files
module read_vmec 
  implicit none
  character(len=256) :: vmec_filename
  real :: bloat_factor     ! factor to expand VMEC domain 
  real :: bloat_distance   ! distance to expand VMEC domain 
  integer :: nzer_factor      ! Zernike resolution parameter
  integer :: nzer_manual      ! Zernike resolution parameter

#ifdef USEST
  integer :: nfp, lasym
  real, allocatable :: rbc(:), zbs(:)
  real, allocatable :: rmnc(:,:), zmns(:,:)!, lmns(:,:), gmnc(:,:)
  real, allocatable :: rmncz(:,:), zmnsz(:,:)!, lmnsz(:,:), gmncz(:,:)
  real, allocatable :: rbs(:), zbc(:)
  real, allocatable :: rmns(:,:), zmnc(:,:)!, lmnc(:,:)
  real, allocatable :: rmnsz(:,:), zmncz(:,:)!, lmncz(:,:)
  real, allocatable :: presf(:)!, phiv(:), chiv(:)
  real, allocatable :: bsupumnc(:,:), bsupvmnc(:,:)
  real, allocatable :: bsupumns(:,:), bsupvmns(:,:)
!  real, allocatable :: bsupumncz(:,:), bsupvmncz(:,:)
!  real, allocatable :: bsupumnsz(:,:), bsupvmnsz(:,:)
!  real, allocatable :: raxiscc(:), zaxiscs(:)
  integer, allocatable :: mb(:), nb(:), mb_nyq(:), nb_nyq(:)
  real, allocatable :: xmv(:), xnv(:), xnv_nyq(:), xmv_nyq(:)
  integer :: mn_mode, mn_mode_nyq, ns, n_tor, n_zer, m_pol, n_quad
  integer :: m_nyq, n_nyq, n_zer_nyq  
  real, allocatable :: s_vmec(:), quad(:,:) 

contains
  ! read VMEC data and put on Zernike or spline basis
  subroutine process_vmec(myrank)
    use math
    implicit none
    
    integer, intent(in) :: myrank
    integer :: i 

    call read_vmec_nc(myrank)
    ! real to integer
    mb = nint(xmv)
    nb = nint(xnv)
    mb_nyq = nint(xmv_nyq)
    nb_nyq = nint(xnv_nyq)
    ! boundary coefficients
    rbc = rmnc(:,ns)
    zbs = zmns(:,ns)
    if(lasym.eq.1) then
      rbs = rmns(:,ns)
      zbc = zmnc(:,ns)
    endif
    ! Change from half to full mesh
    if(myrank.eq.0) print *, 'half mesh to full'
    call half2full(mn_mode_nyq,mb_nyq,bsupumnc)
    call half2full(mn_mode_nyq,mb_nyq,bsupvmnc)
    if(lasym.eq.1) then
      call half2full(mn_mode_nyq,mb_nyq,bsupumns)
      call half2full(mn_mode_nyq,mb_nyq,bsupvmns)
    end if

    ! put VMEC data on Zernike basis
    ! radial grid
    do i = 1, ns
      s_vmec(i) = 1.*(i-1)/(ns-1)
    end do

    ! bloat VMEC domain
    if(bloat_factor.ne.0 .or. bloat_distance.ne.0) then
      call bloat_domain(myrank)
    end if

    ! calculate Gauss quadradure
    call gaussquad(n_quad,quad)

    ! perform Zernike transform
    call zernike_transform(mn_mode,mb,rmnc,rmncz)
    if(myrank.eq.0) print *, 'rmnc transformed'
    call zernike_transform(mn_mode,mb,zmns,zmnsz)
    if(myrank.eq.0) print *, 'zmns transformed'
!    call zernike_transform(mn_mode_nyq,mb_nyq,bsupumnc,bsupumncz)
!    if(myrank.eq.0) print *, 'bsupumnc transformed'
!    call zernike_transform(mn_mode_nyq,mb_nyq,bsupvmnc,bsupvmncz)
!    if(myrank.eq.0) print *, 'bsupvmnc transformed'
    if(lasym.eq.1) then
      call zernike_transform(mn_mode,mb,rmns,rmnsz)
      if(myrank.eq.0) print *, 'rmns transformed'
      call zernike_transform(mn_mode,mb,zmnc,zmncz)
      if(myrank.eq.0) print *, 'zmnc transformed'
!      call zernike_transform(mn_mode_nyq,mb_nyq,bsupumns,bsupumnsz)
!      if(myrank.eq.0) print *, 'bsupumns transformed'
!      call zernike_transform(mn_mode_nyq,mb_nyq,bsupvmns,bsupvmnsz)
!      if(myrank.eq.0) print *, 'bsupvmns transformed'
    endif
    ! Make sure pressure is positive
    if (presf(ns).lt.0) presf = presf - presf(ns)
  end subroutine process_vmec

  subroutine allocate_vmec()
    implicit none

    allocate(xmv(mn_mode))
    allocate(xnv(mn_mode))
    allocate(mb(mn_mode))
    allocate(nb(mn_mode))
    allocate(xmv_nyq(mn_mode_nyq))
    allocate(xnv_nyq(mn_mode_nyq))
    allocate(mb_nyq(mn_mode_nyq))
    allocate(nb_nyq(mn_mode_nyq))
    allocate(presf(ns))
    allocate(rmnc(mn_mode,ns))
    allocate(zmns(mn_mode,ns))
    allocate(bsupumnc(mn_mode_nyq,ns))
    allocate(bsupvmnc(mn_mode_nyq,ns))
    allocate(rbc(mn_mode))
    allocate(zbs(mn_mode))
! Set Zernike polynomial resolution
! Set default
    if(bloat_factor.ne.0 .or. bloat_distance.ne.0) then
      n_zer = m_pol*1 ! for free-boundary
    else 
      n_zer = m_pol*2 ! for fixed-boundary
    end if
! Manual resolution
    if(nzer_manual.ge.0 .and. nzer_manual.ge.n_zer) then
      n_zer = nzer_manual
    end if
! Use scale factor only if manual resolution not set
    if(nzer_factor.ge.0 .and. nzer_manual.lt.0) then
      n_zer = m_pol*nzer_factor
    end if
    allocate(rmncz(mn_mode,n_zer+1))
    allocate(zmnsz(mn_mode,n_zer+1))
    !allocate(bsupumncz(mn_mode_nyq,n_zer+1))
    !allocate(bsupvmncz(mn_mode_nyq,n_zer+1))
    if(lasym.eq.1) then
      allocate(rmns(mn_mode,ns))
      allocate(zmnc(mn_mode,ns))
      allocate(rbs(mn_mode))
      allocate(zbc(mn_mode))
      allocate(rmnsz(mn_mode,n_zer+1))
      allocate(zmncz(mn_mode,n_zer+1))
      allocate(bsupumns(mn_mode_nyq,ns))
      allocate(bsupvmns(mn_mode_nyq,ns))
      !allocate(bsupumnsz(mn_mode_nyq,n_zer+1))
      !allocate(bsupvmnsz(mn_mode_nyq,n_zer+1))
    endif
    allocate(s_vmec(ns))
    n_quad = 1*n_zer + 1
    allocate(quad(2,n_quad))
  end subroutine allocate_vmec

  subroutine read_vmec_nc(myrank)
    use netcdf 

    implicit none
    
    integer, intent(in) :: myrank 
    integer :: ierr , ncid, id 

    ! Open NetCDF file
    if(myrank.eq.0) print *, 'Opening VMEC file'
    ierr = nf90_open(trim(vmec_filename), nf90_nowrite, ncid)
    if(ierr.ne.0) then 
      if(myrank.eq.0) print *, 'Failed to open VMEC file'
      call safestop(5)  
    end if

    ! Get dimension data
    ierr = nf90_inq_dimid(ncid, "mn_mode", id)
    ierr = ierr + nf90_inquire_dimension(ncid, id, len=mn_mode)
    if(ierr.ne.0) call safestop(5) 
    if(myrank.eq.0) print *, 'mn_mode = ', mn_mode
     
    ierr = nf90_inq_dimid(ncid, "mn_mode_nyq", id)
    ierr = ierr + nf90_inquire_dimension(ncid, id, len=mn_mode_nyq)
    if(ierr.ne.0) call safestop(5) 
    if(myrank.eq.0) print *, 'mn_mode_nyq = ', mn_mode_nyq

    ! Get constants
    ierr = nf90_inq_varid(ncid, "ns", id)
    ierr = ierr + nf90_get_var(ncid, id, ns)
    if(ierr.ne.0) call safestop(5) 
    if(myrank.eq.0) print *, 'ns = ', ns
    
    ierr = nf90_inq_varid(ncid, "ntor", id)
    ierr = ierr + nf90_get_var(ncid, id, n_tor)
    if(ierr.ne.0) call safestop(5) 
    if(myrank.eq.0) print *, 'n_tor = ', n_tor

    ierr = nf90_inq_varid(ncid, "mpol", id)
    ierr = ierr + nf90_get_var(ncid, id, m_pol)
    if(ierr.ne.0) call safestop(5) 
    if(myrank.eq.0) print *, 'm_pol = ', m_pol

    ierr = nf90_inq_varid(ncid, "nfp", id)
    ierr = ierr + nf90_get_var(ncid, id, nfp)
    if(ierr.ne.0) call safestop(5) 
    if(myrank.eq.0) print *, 'nfp = ', nfp 

    ierr = nf90_inq_varid(ncid, "lasym__logical__", id)
    ierr = ierr + nf90_get_var(ncid, id, lasym)
    if(ierr.ne.0) call safestop(5) 
    if(myrank.eq.0) print *, 'lasym = ', lasym 

    call allocate_vmec()
    if(myrank.eq.0) print *, 'n_zer = ', n_zer
    if(myrank.eq.0) print *, 'n_quad = ', n_quad
    ! Get 1D arrays xmv & xnv
    ierr = nf90_inq_varid(ncid, "xm", id)
    ierr = ierr + nf90_get_var(ncid, id, xmv)
    if(ierr.ne.0) call safestop(5) 
    if(myrank.eq.0) print *, 'xmv read'
    ierr = nf90_inq_varid(ncid, "xn", id)
    ierr = ierr + nf90_get_var(ncid, id, xnv)
    if(ierr.ne.0) call safestop(5) 
    xnv = -xnv ! change sign for consistency
    if(myrank.eq.0) print *, 'xnv read'
    ierr = nf90_inq_varid(ncid, "xm_nyq", id)
    ierr = ierr + nf90_get_var(ncid, id, xmv_nyq)
    if(ierr.ne.0) call safestop(5) 
    if(myrank.eq.0) print *, 'xmv_nyq read'
    ierr = nf90_inq_varid(ncid, "xn_nyq", id)
    ierr = ierr + nf90_get_var(ncid, id, xnv_nyq)
    if(ierr.ne.0) call safestop(5) 
    xnv_nyq = -xnv_nyq ! change sign for consistency
    if(myrank.eq.0) print *, 'xnv_nyq read'

    ! Get 1D array presf, phiv, chiv
    ierr = nf90_inq_varid(ncid, "presf", id)
    ierr = ierr + nf90_get_var(ncid, id, presf)
    if(ierr.ne.0) call safestop(5) 
    if(myrank.eq.0) print *, 'presf read'

    ! Get 2D array rmnc, zmns, lmns 
    ierr = nf90_inq_varid(ncid, "rmnc", id)
    ierr = ierr + nf90_get_var(ncid, id, rmnc)
    if(ierr.ne.0) call safestop(5) 
    if(myrank.eq.0) print *, 'rmnc read'
    ierr = nf90_inq_varid(ncid, "zmns", id)
    ierr = ierr + nf90_get_var(ncid, id, zmns)
    if(ierr.ne.0) call safestop(5) 
    if(myrank.eq.0) print *, 'zmns read'
    if(lasym.eq.1) then
      ierr = nf90_inq_varid(ncid, "rmns", id)
      ierr = ierr + nf90_get_var(ncid, id, rmns)
      if(ierr.ne.0) call safestop(5) 
      if(myrank.eq.0) print *, 'rmns read'
      ierr = nf90_inq_varid(ncid, "zmnc", id)
      ierr = ierr + nf90_get_var(ncid, id, zmnc)
      if(ierr.ne.0) call safestop(5) 
      if(myrank.eq.0) print *, 'zmnc read'
    endif

    ! Get 2D array gmnc, bsupumnc, bsupsmnv 
    ierr = nf90_inq_varid(ncid, "bsupumnc", id)
    ierr = ierr + nf90_get_var(ncid, id, bsupumnc)
    if(ierr.ne.0) call safestop(5) 
    if(myrank.eq.0) print *, 'bsupumnc read'
    ierr = nf90_inq_varid(ncid, "bsupvmnc", id)
    ierr = ierr + nf90_get_var(ncid, id, bsupvmnc)
    if(ierr.ne.0) call safestop(5) 
    if(myrank.eq.0) print *, 'bsupvmnc read'
    if(lasym.eq.1) then
      ierr = nf90_inq_varid(ncid, "bsupumns", id)
      ierr = ierr + nf90_get_var(ncid, id, bsupumns)
      if(ierr.ne.0) call safestop(5) 
      if(myrank.eq.0) print *, 'bsupumns read'
      ierr = nf90_inq_varid(ncid, "bsupvmns", id)
      ierr = ierr + nf90_get_var(ncid, id, bsupvmns)
      if(ierr.ne.0) call safestop(5) 
      if(myrank.eq.0) print *, 'bsupvmns read'
    endif
  end subroutine read_vmec_nc

  ! read boundary geometry from file
  subroutine read_boundary_geometry(myrank)
    use read_ascii 
    implicit none

    integer, intent(in) :: myrank 
    character(len=256) :: boundary_fname = "boundary"

    mn_mode = 0
    call read_ascii_column(boundary_fname, nb, mn_mode, icol=1)
    call read_ascii_column(boundary_fname, mb, mn_mode, icol=2)
    call read_ascii_column(boundary_fname, rbc, mn_mode, icol=3)
    call read_ascii_column(boundary_fname, zbs, mn_mode, icol=6)
    if(mn_mode.eq.0) then 
      if(myrank.eq.0) print *, 'Error: could not find geometry'
      call safestop(5)
    else
      if(myrank.eq.0) print *, 'Boundary geometry read'
    end if
  endsubroutine read_boundary_geometry

  subroutine destroy_vmec
    implicit none

    deallocate(xmv)
    deallocate(xnv)
    deallocate(mb)
    deallocate(nb)
    deallocate(xmv_nyq)
    deallocate(xnv_nyq)
    deallocate(mb_nyq)
    deallocate(nb_nyq)
    deallocate(rmnc)
    deallocate(zmns)
    deallocate(rmncz)
    deallocate(zmnsz)
    deallocate(bsupumnc)
    deallocate(bsupvmnc)
    deallocate(rbc)
    deallocate(zbs)
!    deallocate(bsupumncz)
!    deallocate(bsupvmncz)
    if(lasym.eq.1) then
      deallocate(rmns)
      deallocate(zmnc)
      deallocate(rmnsz)
      deallocate(zmncz)
      deallocate(rbs)
      deallocate(zbc)
      deallocate(bsupumns)
      deallocate(bsupvmns)
!      deallocate(bsupumnsz)
!      deallocate(bsupvmnsz)
    endif
    deallocate(s_vmec)
    deallocate(quad)
    deallocate(presf)
  end subroutine destroy_vmec

  ! evaluate zernike-representation at r
  pure subroutine zernike_evaluate(r,mn,m,fmn,fout)
    implicit none

    real, intent(in) :: r 
    integer, intent(in) :: mn 
    integer, dimension(mn), intent(in) :: m 
    real, dimension(mn,n_zer+1), intent(in) :: fmn
    real, dimension(mn), intent(out) :: fout
    integer :: i, j
    real :: zmn 

    fout = 0.
    do i = 1, mn
       do j = m(i), n_zer, 2
          call zernike_polynomial(r,j,m(i),zmn) 
          fout(i) = fout(i) + fmn(i,j+1)*zmn
       end do
    end do
  end subroutine zernike_evaluate


  ! radial zernike transform on VMEC data
  pure subroutine zernike_transform(mn,m,fmn,fout)
    implicit none

    integer, intent(in) :: mn 
    integer, dimension(mn), intent(in) :: m 
    real, dimension(mn,ns), intent(in) :: fmn
    real, dimension(mn) :: ftemp
    real, dimension(mn,n_zer+1), intent(out) :: fout
    integer :: i, j, k, nt
    real :: rho, zmn 

    fout = 0.

    do k = 1, n_quad
       rho = .5*(1.+quad(1,k))
       call vmec_interpl(rho,mn,m,fmn,ftemp)
       do i = 1, mn
          do j = m(i), n_zer-2, 2
             call zernike_polynomial(rho,j,m(i),zmn) 
             fout(i,j+1) = fout(i,j+1) + ftemp(i)*zmn*(j+1)*rho*quad(2,k)
          end do
       end do
    end do

    ! make sure the coefficients sum up to the boundary value 
    do i = 1, mn
       zmn = fmn(i,ns) 
       do j = m(i), n_zer, 2
          if (j.lt.(n_zer-1)) then
             zmn = zmn - fout(i,j+1)
          else
             fout(i,j+1) = zmn
          end if 
       end do
    end do
       
  end subroutine zernike_transform

  pure recursive subroutine zernike_polynomial(r, n, m, z)
    implicit none

    real, intent(in) :: r 
    integer, intent(in) :: n, m 
    real, intent(out) :: z
    real :: z1, z2 

    if (n.eq.m) then
       z = r**m
    else if (n.eq.(m+2)) then
       z = (m+2)*r**(m+2) - (m+1)*r**m
    else if (mod(n-m, 2) .eq. 0) then
       call zernike_polynomial(r,n-2,m,z1)
       call zernike_polynomial(r,n-4,m,z2)
       z = ((4*(n-1)*r**2-(n-2.+m)**2/(n-2)-(n-m+0.)**2/n)*z1 &
            -((n-2.)**2-m**2)/(n-2)*z2)*n/(n**2-m**2)
    else 
       z = 0.
    end if    

  end subroutine zernike_polynomial

  pure subroutine vmec_interpl(r,mn,m,fmn,fout)
    implicit none

    real, intent(in) :: r
    integer, intent(in) :: mn 
    integer, dimension(mn), intent(in) :: m 
    real :: r2n, ds 
    integer :: i, js
    real, dimension(mn,ns), intent(in) :: fmn
    real, dimension(mn,ns) :: f
    real, dimension(mn), intent(out) :: fout
    real, dimension(mn) :: df, df0, df1 
    real, dimension(mn,4) :: a 

    if (r.eq.0.) then
       fout = fmn(:,1)
    else
       f(:,1) = 0.
       do i = 1, ns-1
          f(:,i+1) = fmn(:,i+1)*(1.*i/(ns-1))**(-.5*m+1.)
       end do
       do i = 1, mn
          if(m(i).ne.0.) f(i,1) = 0.
       end do
       r2n = r**2*(ns-1)
       js = ceiling(r2n)
       if(js.gt.ns-1) js = ns-1
       ds = r2n - (js-1)
       df = (f(:,js+1) - f(:,js))
       if(js.eq.1) then
          df0 = df
          df1 = (f(:,js+2) - f(:,js  ))/2.
       else if(js.eq.ns-1) then
          df0 = (f(:,js+1) - f(:,js-1))/2.
          df1 = df
       else
          df0 = (f(:,js+1) - f(:,js-1))/2.
          df1 = (f(:,js+2) - f(:,js  ))/2.
       end if
       a(:,1) = f(:,js)
       a(:,2) = df0
       a(:,3) = (3.*df - (2.*df0+df1))
       a(:,4) = (-2.*df + (df0+df1))
       fout = a(:,1) + a(:,2)*ds + a(:,3)*ds**2 + a(:,4)*ds**3
       fout = fout*r**(m-2)
    end if

  end subroutine vmec_interpl

  ! half mesh VMEC data to full mesh
  subroutine half2full(mn,m,fmn)
    implicit none

    integer, intent(in) :: mn 
    integer, dimension(mn), intent(in) :: m 
    real, dimension(mn,ns), intent(inout) :: fmn
    real, dimension(mn,ns) :: ftemp
    integer :: i 

    ftemp(:,1) = fmn(:,2)*1.5 - fmn(:,3)*.5
    ftemp(:,ns) = fmn(:,ns)*1.5 - fmn(:,ns-1)*.5
    do i = 2, ns-1
       ftemp(:,i) = .5*(fmn(:,i+1) + fmn(:,i)) 
    end do
    fmn = ftemp
    do i = 1, mn
       if(m(i).ne.0) fmn(i,1) = 0
    end do 
  end subroutine half2full

  ! Adapted from https://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature#Fortran
  subroutine gaussquad(n,quad)
    use math
    implicit none
    integer, intent(in) :: n
    real, intent(inout) :: quad(2, n)
    integer :: i, k, iter
    real :: x, f, df, dx
    real :: p0(n+1), p1(n+1), tmp(n+1)
  
    p0 = 0.
    p1 = 0.
    tmp = 0.
   
    p0(1) = 1.
    p1(1:2) = [1., 0.]
   
    do k = 2, n
      tmp(1:k+1) = ((2*k-1)*p1(1:k+1)-(k-1)*[0., 0.,p0(1:k-1)])/k
      p0 = p1; p1 = tmp
    end do
    do i = 1, n
      x = cos(pi*(i-0.25)/(n+0.5))
      do iter = 1, 10
        f = p1(1); df = 0.
        do k = 2, size(p1)
          df = f + x*df
          f  = p1(k) + x * f
        end do
        dx =  f / df
        x = x - dx
        if (abs(dx) .lt. 10.*epsilon(dx)) exit
      end do
      quad(1,i) = x
      quad(2,i) = 2./((1.-x**2)*df**2)
    end do

  end subroutine gaussquad 

  subroutine bloat_domain(myrank)
    use math
    use spline 
    implicit none

    real, allocatable :: rreal(:,:,:), zreal(:,:,:)
    real, dimension(mn_mode) :: co, sn 
    real, dimension(ns) :: s_bloat 
    real, allocatable :: theta(:), zeta(:)
    real, allocatable :: rmnc_temp(:,:), zmns_temp(:,:)
    real, allocatable :: rmns_temp(:,:), zmnc_temp(:,:)
    integer :: i, j, k, l, nt, nz 
    integer, intent(in) :: myrank
    type(spline1d) :: rreal_spline, zreal_spline 
    real :: dr, dz

    ! grid in theta and zeta
    nt = 2*m_pol-1
    nz = 2*n_tor+1
    allocate(theta(nt))
    allocate(zeta(nz))
    do i = 1, nt 
      theta(i) = twopi*(i-1)/nt
    end do
    do i = 1, nz 
      zeta(i) = twopi*(i-1)/(nfp*nz)
    end do

    ! calculate R and Z in real space
    allocate(rreal(ns,nt,nz))
    allocate(zreal(ns,nt,nz))
    rreal = 0.
    zreal = 0.
    do i = 1, ns 
      do j = 1, nt 
        do k = 1, nz 
          co = cos(xmv*theta(j)+xnv*zeta(k))
          sn = sin(xmv*theta(j)+xnv*zeta(k))
          dr = 0.
          dz = 0.
          do l = 1, mn_mode  
            rreal(i,j,k) = rreal(i,j,k) + rmnc(l,i)*co(l)
            zreal(i,j,k) = zreal(i,j,k) + zmns(l,i)*sn(l)
            dr = dr - rmnc(l,i)*sn(l)*xmv(l)
            dz = dz + zmns(l,i)*co(l)*xmv(l)
            if(lasym.eq.1) then
              rreal(i,j,k) = rreal(i,j,k) + rmns(l,i)*sn(l)
              zreal(i,j,k) = zreal(i,j,k) + zmnc(l,i)*co(l)
              dr = dr + rmns(l,i)*co(l)*xmv(l)
              dz = dz - zmnc(l,i)*sn(l)*xmv(l)
            end if
          end do
          if(bloat_distance.ne.0) then
            rreal(i,j,k) = rreal(i,j,k) &
              + bloat_distance*sqrt(s_vmec(i))*dz/sqrt(dz**2+dr**2)
            zreal(i,j,k) = zreal(i,j,k) &
              - bloat_distance*sqrt(s_vmec(i))*dr/sqrt(dz**2+dr**2)
          end if 
        end do
      end do
    end do

    if(bloat_distance.ne.0) then
      bloat_factor = 0  ! bloat_distance overrides bloat_factor
      if(myrank.eq.0) print *, 'expanding domain by distance'
    else if(bloat_factor.ne.0) then
      if(myrank.eq.0) print *, 'expanding domain by ratio'
      ! extrapolate R and Z in real space
      s_bloat = s_vmec*bloat_factor**2
      do j = 1, nt 
        do k = 1, nz 
          call create_spline(rreal_spline, ns, s_vmec, rreal(:,j,k))
          call create_spline(zreal_spline, ns, s_vmec, zreal(:,j,k))
          do i = 2, ns 
            call evaluate_spline(rreal_spline, s_bloat(i), rreal(i,j,k),extrapolate=1)
            call evaluate_spline(zreal_spline, s_bloat(i), zreal(i,j,k),extrapolate=1)
          end do
          call destroy_spline(rreal_spline)
          call destroy_spline(zreal_spline)
        end do
      end do
    end if

    ! transform back into Fourier space
    allocate(rmnc_temp(mn_mode,ns))
    allocate(zmns_temp(mn_mode,ns)) 
    rmnc_temp = 0.
    zmns_temp = 0.
    if(lasym.eq.1) then
      allocate(rmns_temp(mn_mode,ns))
      allocate(zmnc_temp(mn_mode,ns)) 
      rmns_temp = 0.
      zmnc_temp = 0.
    end if
    do j = 1, nt 
      do k = 1, nz
        co = cos(xmv*theta(j)+xnv*zeta(k))
        sn = sin(xmv*theta(j)+xnv*zeta(k))
        do i = 2, ns 
          do l = 1, mn_mode  
            rmnc_temp(l,i) = rmnc_temp(l,i) + rreal(i,j,k)*co(l)
            zmns_temp(l,i) = zmns_temp(l,i) + zreal(i,j,k)*sn(l)
            if(lasym.eq.1) then
              rmns_temp(l,i) = rmns_temp(l,i) + rreal(i,j,k)*sn(l)
              zmnc_temp(l,i) = zmnc_temp(l,i) + zreal(i,j,k)*co(l)
            end if
          end do
        end do
      end do
    end do
    rmnc_temp = rmnc_temp*2./(nt*nz)
    zmns_temp = zmns_temp*2./(nt*nz)
    rmnc_temp(1,:) = rmnc_temp(1,:)*.5
    zmns_temp(1,:) = zmns_temp(1,:)*.5
    rmnc(:,2:) = rmnc_temp(:,2:)
    zmns(:,2:) = zmns_temp(:,2:)
    deallocate(rmnc_temp,zmns_temp)
    if(lasym.eq.1) then
      rmns_temp = rmns_temp*2./(nt*nz)
      zmnc_temp = zmnc_temp*2./(nt*nz)
      rmns_temp(1,:) = rmns_temp(1,:)*.5
      zmnc_temp(1,:) = zmnc_temp(1,:)*.5
      rmns(:,2:) = rmns_temp(:,2:)
      zmnc(:,2:) = zmnc_temp(:,2:)
      deallocate(rmns_temp,zmnc_temp)
    end if
    deallocate(theta,zeta)
    deallocate(rreal,zreal)

  end subroutine bloat_domain

#endif
end module read_vmec 
