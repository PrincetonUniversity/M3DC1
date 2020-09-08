! This module reads     presf = presf*pi*in data from VMEC files
module read_vmec 
  use spline
  implicit none
  character(len=256) :: vmec_filename

#ifdef USEST
  integer :: nfp
  real, allocatable :: rbc(:), zbs(:)
  real, allocatable :: rmnc(:,:), zmns(:,:), lmns(:,:), gmnc(:,:)
  real, allocatable :: rmncz(:,:), zmnsz(:,:), lmnsz(:,:), gmncz(:,:)
  real, allocatable :: bsupumnc(:,:), bsupvmnc(:,:)
  real, allocatable :: bsupumncz(:,:), bsupvmncz(:,:)
  real, allocatable :: presf(:), phiv(:), chiv(:)
!  real, allocatable :: raxiscc(:), zaxiscs(:)
  integer, allocatable :: mb(:), nb(:), mb_nyq(:), nb_nyq(:)
  real, allocatable :: xmv(:), xnv(:), xnv_nyq(:), xmv_nyq(:)
  integer :: mn_mode, mn_mode_nyq, ns, n_tor, n_zer, m_pol 
  real, allocatable :: s_vmec(:) 
  type(spline1d) :: presf_spline      ! total pressure
  type(spline1d) :: phiv_spline       ! toroidal flux
  type(spline1d) :: chiv_spline       ! poloidal flux
   

contains
  ! read VMEC data and put on Zernike or spline basis
  subroutine process_vmec(myrank)
    use math
    implicit none
    
    integer, intent(in) :: myrank
    integer :: i 

    !call read_vmec_h5(myrank)
    call read_vmec_nc(myrank)
    ! real to integer
    mb = nint(xmv)
    nb = nint(xnv)
!    mb_nyq = nint(xmv_nyq)
!    nb_nyq = nint(xnv_nyq)
    ! normalize pressure
    presf = presf*pi*4e-7
    ! boundary coefficients
    rbc = rmnc(:,ns)
    zbs = zmns(:,ns)

    ! Change from half to full mesh
    call half2full(mn_mode,mb,lmns)
!    call half2full(mn_mode_nyq,mb_nyq,gmnc)
!    call half2full(mn_mode_nyq,mb_nyq,bsupumnc)
!    call half2full(mn_mode_nyq,mb_nyq,bsupvmnc)
    if(myrank.eq.0) print *, 'half mesh to full'

    ! put VMEC data on Zernike basis
    ! radial grid
    do i = 1, ns
       s_vmec(i) = 1.*(i-1)/(ns-1)
    end do
    ! perform Zernike transform
    call zernike_transform(rmnc,rmncz)
    if(myrank.eq.0) print *, 'rmnc transformed'
    call zernike_transform(zmns,zmnsz)
    if(myrank.eq.0) print *, 'zmns transformed'
    call zernike_transform(lmns,lmnsz)
    if(myrank.eq.0) print *, 'lmns transformed'
!    call zernike_transform(mn_mode_nyq,mb_nyq,gmnc,gmncz)
!    if(myrank.eq.0) print *, 'gmnc transformed'
!    call zernike_transform(mn_mode_nyq,mb_nyq,bsupumnc,bsupumncz)
!    if(myrank.eq.0) print *, 'bsupumnc transformed'
!    call zernike_transform(mn_mode_nyq,mb_nyq,bsupvmnc,bsupvmncz)
!    if(myrank.eq.0) print *, 'bsupvmnc transformed'
    ! Make sure pressure is positive
    if (presf(ns).lt.0) presf = presf - presf(ns)
    ! 1D spline for pressure
    call create_spline(presf_spline, ns, s_vmec, presf)
    call create_spline(phiv_spline, ns, s_vmec, phiv)
    call create_spline(chiv_spline, ns, s_vmec, chiv)

!    do i = 1, mn_mode
!        print *, mb(i), nb(i)
!        print *, sum(rmncz(i,:)) 
!        print *, rmnc(i,ns) 
!       do j = mb(i), n_zer, 2
!          !if (abs(rmncz(i,j+1).gt.0.1)) then
!            print *, j, mb(i), nb(i)
!            print *, rmncz(i,j+1)
!          !end if 
!          do k = ns, ns
!             rho = sqrt((k-1.)/(ns-1))
!             call zernike_polynomial(rho,j,mb(i),zmn) 
!             !print *, rho, j, mb(i), zmn
!          end do
!       end do
!    end do


  end subroutine process_vmec

  subroutine allocate_vmec()
    implicit none

    allocate(xmv(mn_mode))
    allocate(xnv(mn_mode))
    allocate(mb(mn_mode))
    allocate(nb(mn_mode))
!    allocate(xmv_nyq(mn_mode_nyq))
!    allocate(xnv_nyq(mn_mode_nyq))
!    allocate(mb_nyq(mn_mode_nyq))
!    allocate(nb_nyq(mn_mode_nyq))
    allocate(presf(ns))
    allocate(phiv(ns))
    allocate(chiv(ns))
    allocate(rmnc(mn_mode,ns))
    allocate(zmns(mn_mode,ns))
    allocate(lmns(mn_mode,ns))
!    allocate(gmnc(mn_mode_nyq,ns))
!    allocate(bsupumnc(mn_mode_nyq,ns))
!    allocate(bsupvmnc(mn_mode_nyq,ns))
    allocate(rbc(mn_mode))
    allocate(zbs(mn_mode))
    n_zer = m_pol*3
    allocate(rmncz(mn_mode,n_zer+1))
    allocate(zmnsz(mn_mode,n_zer+1))
    allocate(lmnsz(mn_mode,n_zer+1))
!    allocate(gmncz(mn_mode_nyq,n_zer+1))
!    allocate(bsupumncz(mn_mode_nyq,n_zer+1))
!    allocate(bsupvmncz(mn_mode_nyq,n_zer+1))
    allocate(s_vmec(ns))
  end subroutine allocate_vmec

  subroutine read_vmec_nc(myrank)
    use netcdf 

    implicit none
    
    integer, intent(in) :: myrank 
    integer :: ierr , ncid, id 

    ! Open NetCDF file
    ierr = nf90_open(trim(vmec_filename), nf90_nowrite, ncid)
    if(ierr.ne.0) call safestop(5) 
    if(myrank.eq.0) print *, 'Opening VMEC file'
    
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

    call allocate_vmec()
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
!    ierr = nf90_inq_varid(ncid, "xm_nyq", id)
!    ierr = ierr + nf90_get_var(ncid, id, xmv_nyq)
!    if(ierr.ne.0) call safestop(5) 
!    if(myrank.eq.0) print *, 'xmv_nyq read'
!    ierr = nf90_inq_varid(ncid, "xn_nyq", id)
!    ierr = ierr + nf90_get_var(ncid, id, xnv_nyq)
!    if(ierr.ne.0) call safestop(5) 
!    xnv_nyq = -xnv_nyq ! change sign for consistency
!    if(myrank.eq.0) print *, 'xnv_nyq read'


    ! Get 1D array presf, phiv, chiv
    ierr = nf90_inq_varid(ncid, "presf", id)
    ierr = ierr + nf90_get_var(ncid, id, presf)
    if(ierr.ne.0) call safestop(5) 
    if(myrank.eq.0) print *, 'presf read'
    ierr = nf90_inq_varid(ncid, "phi", id)
    ierr = ierr + nf90_get_var(ncid, id, phiv)
    if(ierr.ne.0) call safestop(5) 
    if(myrank.eq.0) print *, 'phiv read'
    ierr = nf90_inq_varid(ncid, "chi", id)
    ierr = ierr + nf90_get_var(ncid, id, chiv)
    if(ierr.ne.0) call safestop(5) 
    if(myrank.eq.0) print *, 'chiv read'

    ! Get 2D array rmnc, zmns, lmns 
    ierr = nf90_inq_varid(ncid, "rmnc", id)
    ierr = ierr + nf90_get_var(ncid, id, rmnc)
    if(ierr.ne.0) call safestop(5) 
    if(myrank.eq.0) print *, 'rmnc read'
    ierr = nf90_inq_varid(ncid, "zmns", id)
    ierr = ierr + nf90_get_var(ncid, id, zmns)
    if(ierr.ne.0) call safestop(5) 
    if(myrank.eq.0) print *, 'zmns read'
    ierr = nf90_inq_varid(ncid, "lmns", id)
    ierr = ierr + nf90_get_var(ncid, id, lmns)
    if(ierr.ne.0) call safestop(5) 
    if(myrank.eq.0) print *, 'lmns read'

    ! Get 2D array gmnc, bsupumnc, bsupsmnv 
!    ierr = nf90_inq_varid(ncid, "gmnc", id)
!    ierr = ierr + nf90_get_var(ncid, id, gmnc)
!    if(ierr.ne.0) call safestop(5) 
!    if(myrank.eq.0) print *, 'gmnc read'
!    ierr = nf90_inq_varid(ncid, "bsupumnc", id)
!    ierr = ierr + nf90_get_var(ncid, id, bsupumnc)
!    if(ierr.ne.0) call safestop(5) 
!    if(myrank.eq.0) print *, 'bsupumnc read'
!    ierr = nf90_inq_varid(ncid, "bsupvmnc", id)
!    ierr = ierr + nf90_get_var(ncid, id, bsupvmnc)
!    if(ierr.ne.0) call safestop(5) 
!    if(myrank.eq.0) print *, 'bsupvmnc read'

  end subroutine read_vmec_nc

  subroutine read_vmec_h5(myrank)
    use hdf5
    
    implicit none
    

    integer, intent(in) :: myrank 
    integer :: error 

    integer(HID_T) :: file_id, dset_id, attr_id
    integer(HSIZE_T), dimension(1) :: dim0 = 1
    integer(HSIZE_T), dimension(1) :: dim1 
    integer(HSIZE_T), dimension(2) :: dim2

    vmec_filename = "geometry.h5"
    call h5open_f(error)
    call h5fopen_f(vmec_filename, H5F_ACC_RDONLY_F, file_id, error)

    ! read attributes ns, n_tor, m_pol, nfp, and mn_mode
    call h5gopen_f(file_id, '/', dset_id, error)
    call h5aopen_name_f(dset_id, 'ns', attr_id, error)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, ns, dim0, error)
    call h5aclose_f(attr_id, error)
    call h5aopen_name_f(dset_id, 'ntor', attr_id, error)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, n_tor, dim0, error)
    call h5aclose_f(attr_id, error)
    call h5aopen_name_f(dset_id, 'mpol', attr_id, error)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, m_pol, dim0, error)
    call h5aclose_f(attr_id, error)
    call h5aopen_name_f(dset_id, 'nfp', attr_id, error)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, nfp, dim0, error)
    call h5aclose_f(attr_id, error)
    call h5aopen_name_f(dset_id, 'mnmode', attr_id, error)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, mn_mode, dim0, error)
    call h5aclose_f(attr_id, error)
    call h5gclose_f(dset_id, error)
    mn_mode_nyq = mn_mode
    if(myrank.eq.0) print *, 'attributes read'
    if(myrank.eq.0) print *, mn_mode, mn_mode_nyq, ns, nfp

    call allocate_vmec()

    dim1(1) = mn_mode
    ! read 1d arrays xmv and xnv 
    call h5dopen_f(file_id, 'xm', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, xmv, dim1, error)
    call h5dclose_f(dset_id, error)
    call h5dopen_f(file_id, 'xn', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, xnv, dim1, error)
    call h5dclose_f(dset_id, error)
    if(myrank.eq.0) print *, 'xmv, xnv read'

    dim1(1) = ns
    ! read 1d array presf
    call h5dopen_f(file_id, 'presf', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, presf, dim1, error)
    call h5dclose_f(dset_id, error)
    if(myrank.eq.0) print *, 'presf read'

!    dim1(1) = n_tor+1
!    allocate(raxiscc(n_tor+1))
!    allocate(zaxiscs(n_tor+1))
!    ! read 1d arrays raxisicc and zaxiscs 
!    call h5dopen_f(file_id, 'raxiscc', dset_id, error)
!    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, raxiscc, dim1, error)
!    call h5dclose_f(dset_id, error)
!    call h5dopen_f(file_id, 'zaxiscs', dset_id, error)
!    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, zaxiscs, dim1, error)
!    call h5dclose_f(dset_id, error)
!    if(myrank.eq.0) print *, 'raxiscc, zaxiscs read'
!    if(myrank.eq.0) print *, raxiscc(:) 
!    if(myrank.eq.0) print *, zaxiscs(:) 

    dim2(1) = mn_mode
    dim2(2) = ns
    ! read 2d arrays rmnc and zmns
    call h5dopen_f(file_id, 'rmnc', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, rmnc, dim2, error)
    call h5dclose_f(dset_id, error)
    call h5dopen_f(file_id, 'zmns', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, zmns, dim2, error)
    call h5dclose_f(dset_id, error)
    if(myrank.eq.0) print *, 'rmnc, zmns read'
    call h5dopen_f(file_id, 'bsupumnc', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, bsupumnc, dim2, error)
    call h5dclose_f(dset_id, error)
    call h5dopen_f(file_id, 'bsupvmnc', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, bsupvmnc, dim2, error)
    call h5dclose_f(dset_id, error)
    if(myrank.eq.0) print *, 'bsupumnc, bsupvmnc read'

    call h5fclose_f(file_id, error)
    call h5close_f(error)

    if(mn_mode.eq.0) then 
      if(myrank.eq.0) print *, 'Error: could not find geometry'
      call safestop(5)
    else
      if(myrank.eq.0) print *, 'VMEC geometry read'
    end if

  endsubroutine read_vmec_h5

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
!    deallocate(xmv_nyq)
!    deallocate(xnv_nyq)
!    deallocate(mb_nyq)
!    deallocate(nb_nyq)
    deallocate(rmnc)
    deallocate(zmns)
    deallocate(lmns)
    deallocate(rmncz)
    deallocate(zmnsz)
    deallocate(lmnsz)
!    deallocate(gmnc)
!    deallocate(gmncz)
!    deallocate(bsupumnc)
!    deallocate(bsupvmnc)
!    deallocate(bsupumncz)
!    deallocate(bsupvmncz)
    deallocate(rbc)
    deallocate(zbs)
    deallocate(s_vmec)
    deallocate(presf)
    deallocate(phiv)
    deallocate(chiv)
    call destroy_spline(presf_spline)
    call destroy_spline(phiv_spline)
    call destroy_spline(chiv_spline)
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
  !pure subroutine zernike_transform(mn,m,fmn,fout)
  pure subroutine zernike_transform(fmn,fout)
    implicit none

    !integer, intent(in) :: mn 
    !integer, dimension(mn), intent(in) :: m 
    real, dimension(mn_mode,ns), intent(in) :: fmn
    real, dimension(mn_mode) :: ftemp
    real, dimension(mn_mode,n_zer+1), intent(out) :: fout
    integer :: i, j, k, nt
    real :: rho, zmn 

    fout = 0.

    nt = 1024*n_zer 

    do k = 1, nt
       rho = sqrt((k-1.)/(nt-1))
       call vmec_interpl(rho,fmn,ftemp)
       do i = 1, mn_mode
          do j = mb(i), n_zer-2, 2
             call zernike_polynomial(rho,j,mb(i),zmn) 
             if (k.eq.1 .or. k.eq.nt) then
                fout(i,j+1) = fout(i,j+1) + ftemp(i)*zmn*(j+1)/(2*nt-2) 
             else
                fout(i,j+1) = fout(i,j+1) + ftemp(i)*zmn*(j+1)/(nt-1) 
             end if 
          end do
          ! calculate integral in s, substract end points
          !fout(i,j+1) = (sum(fmn(i,:)*zmn)-fmn(i,1)*zmn(1)/2&
          !            -fmn(i,ns)*zmn(ns)/2)/(ns-1)*(j+1)
       end do
    end do
    ! make sure the coefficients sum up to the boundary value 
    do i = 1, mn_mode
       zmn = fmn(i,ns) 
       do j = mb(i), n_zer, 2
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
       z = ((4*(n-1)*r**2-(n-2.+m)**2/(n-2)-(n-m+0.)**2/n)*z1&
           -((n-2.)**2-m**2)/(n-2)*z2)*n/(n**2-m**2)
    else 
       z = 0.
    end if    

!    z = 0.
!    if (mod(n-m, 2) .eq. 0 .and. r.le.1) then
!    !if (mod(n-m, 2) .eq. 0) then
!       do k = 0, (n-m)/2
!          z = z + r**(n-2*k)*(-1)**k*&
!              factorial(n-k)/(factorial(k)*&
!              factorial((n+m)/2-k)*factorial((n-m)/2-k))
!              binomial(n-k,k)*binomial(n-2*k,(n-m)/2-k)
!       end do
!    end if 

  end subroutine zernike_polynomial

  pure subroutine vmec_interpl(r,fmn,fout)
    implicit none

    real, intent(in) :: r
    !integer, intent(in) :: mn 
    !integer, dimension(mn), intent(in) :: m 
    real :: r2n, ds 
    integer :: i, js
    real, dimension(mn_mode,ns), intent(in) :: fmn
    real, dimension(mn_mode,ns) :: f
    real, dimension(mn_mode), intent(out) :: fout
    real, dimension(mn_mode) :: df, df0, df1 
    real, dimension(mn_mode,4) :: a 

    if (r.eq.0) then
       fout = fmn(:,1)
    else
       do i = 0, ns-1
          f(:,i+1) = fmn(:,i+1)*(1.*i/(ns-1))**(-.5*mb+1) 
       end do
       !f(:,1)=fmn(:,1)
       do i = 1, mn_mode
          if(mb(i).ne.0.) f(i,1) = 0.
       end do
       r2n = r**2*(ns-1)
       js = ceiling(r2n)
       ds = r2n - (js-1)
       df = (f(:,js+1) - f(:,js))
       if(js.eq.1) then
          df0 = df
          df1 = (f(:,js+2) - f(:,js  ))/2
       else if(js.ge.ns-1) then
          df0 = (f(:,js+1) - f(:,js-1))/2
          df1 = df
       else
          df0 = (f(:,js+1) - f(:,js-1))/2
          df1 = (f(:,js+2) - f(:,js  ))/2
       end if
       a(:,1) = f(:,js)
       a(:,2) = df0
       a(:,3) = (3.*df - (2.*df0+df1))
       a(:,4) = (-2.*df + (df0+df1))
       fout = a(:,1) + a(:,2)*ds + a(:,3)*ds**2 + a(:,4)*ds**3
       fout = fout*r**(mb-2)
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

!    fmn(:,1) = fmn(:,2)*1.5 - fmn(:,3)*.5
!    do i = 2, ns-1
!       fmn(:,i) = .5*(fmn(:,i+1) + fmn(:,i)) 
!    end do
!    fmn(:,ns) = fmn(:,ns-1)*2 - fmn(:,ns-2)
  
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


#endif
end module read_vmec 
