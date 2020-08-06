! This module reads     presf = presf*pi*in data from VMEC files
module read_vmec 
  use math
  use spline
  implicit none

#ifdef USEST
  integer :: nfp
  real, allocatable :: rbc(:), zbs(:)
  real, allocatable :: rmnc(:,:), zmns(:,:)
  real, allocatable :: rmncz(:,:), zmnsz(:,:)
  real, allocatable :: bsupumnc(:,:), bsupvmnc(:,:)
  real, allocatable :: bsupumncz(:,:), bsupvmncz(:,:)
  real, allocatable :: presf(:)
!  real, allocatable :: raxiscc(:), zaxiscs(:)
  integer, allocatable :: mb(:), nb(:)
  real, allocatable :: xmv(:), xnv(:)
  integer :: mn_mode, ns, n_tor, n_zer, m_pol 
  real, allocatable :: s_vmec(:), rho_vmec(:) 
  type(spline1d) :: presf_spline      ! Total pressure
   

contains

  subroutine read_vmec_h5(myrank)
    use hdf5
    
    implicit none
    
    character(len=256) :: vmec_filename = "geometry.h5"

    integer, intent(in) :: myrank 
    integer :: error , i, j, k

!    real :: rho, zmn
    
    integer(HID_T) :: file_id, dset_id, attr_id
    integer(HSIZE_T), dimension(1) :: dim0 = 1
    integer(HSIZE_T), dimension(1) :: dim1 
    integer(HSIZE_T), dimension(2) :: dim2

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
    if(myrank.eq.0) print *, 'attributes read'
    if(myrank.eq.0) print *, mn_mode, ns, nfp

    dim1(1) = mn_mode
    allocate(xmv(mn_mode))
    allocate(xnv(mn_mode))
    ! read 1d arrays xmv and xnv 
    call h5dopen_f(file_id, 'xm', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, xmv, dim1, error)
    call h5dclose_f(dset_id, error)
    call h5dopen_f(file_id, 'xn', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, xnv, dim1, error)
    call h5dclose_f(dset_id, error)
    if(myrank.eq.0) print *, 'xmv, xnv read'
    allocate(mb(mn_mode))
    allocate(nb(mn_mode))
    mb = nint(xmv)
    nb = nint(xnv)

    dim1(1) = ns
    allocate(presf(ns))
    ! read 1d array presf
    call h5dopen_f(file_id, 'presf', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, presf, dim1, error)
    call h5dclose_f(dset_id, error)
    presf = presf*pi*4e-7
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
    allocate(rmnc(mn_mode,ns))
    allocate(zmns(mn_mode,ns))
    allocate(bsupumnc(mn_mode,ns))
    allocate(bsupvmnc(mn_mode,ns))
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
    allocate(rbc(mn_mode))
    allocate(zbs(mn_mode))
    rbc = rmnc(:,ns)
    zbs = zmns(:,ns)

    if(mn_mode.eq.0) then 
      if(myrank.eq.0) print *, 'Error: could not find geometry'
      call safestop(5)
    else
      if(myrank.eq.0) print *, 'VMEC geometry read'
    end if

    ! radial grid
    allocate(s_vmec(ns))
!    allocate(rho_vmec(ns))
    do i = 1, ns
       s_vmec(i) = 1.*(i-1)/(ns-1)
!    s_vmec(ns) = 1
!    rho_vmec = sqrt(s_vmec)
    end do
    ! Perform Zernike transform
    n_zer = m_pol
    allocate(rmncz(mn_mode,n_zer+1))
    allocate(zmnsz(mn_mode,n_zer+1))
    call zernike_transform(rmnc,rmncz)
    call zernike_transform(zmns,zmnsz)
    allocate(bsupumncz(mn_mode,n_zer+1))
    allocate(bsupvmncz(mn_mode,n_zer+1))
    call zernike_transform(bsupumnc,bsupumncz)
    call zernike_transform(bsupvmnc,bsupvmncz)
    call create_spline(presf_spline, ns, s_vmec, presf)

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
    deallocate(rmnc)
    deallocate(zmns)
    deallocate(rmncz)
    deallocate(zmnsz)
    deallocate(bsupumnc)
    deallocate(bsupvmnc)
    deallocate(bsupumncz)
    deallocate(bsupvmncz)
    deallocate(rbc)
    deallocate(zbs)
    call destroy_spline(presf_spline)
    deallocate(s_vmec)
!    deallocate(rho_vmec)
  end subroutine destroy_vmec

  ! evaluate zernike-representation at r
  pure subroutine zernike_evaluate(r,fmn,fout)
    implicit none

    real, intent(in) :: r 
    real, dimension(mn_mode,n_zer+1), intent(in) :: fmn
    real, dimension(mn_mode), intent(out) :: fout
    integer :: i, j
    real :: zmn 

    fout = 0.
    do i = 1, mn_mode
       do j = mb(i), n_zer, 2
          call zernike_polynomial(r,j,mb(i),zmn) 
          fout(i) = fout(i) + fmn(i,j+1)*zmn
       end do
    end do
  end subroutine zernike_evaluate


  ! radial zernike transform on VMEC data
  pure subroutine zernike_transform(fmn,fout)
    implicit none

    real, dimension(mn_mode,ns), intent(in) :: fmn
    real, dimension(mn_mode) :: ftemp
    real, dimension(mn_mode,n_zer+1), intent(out) :: fout
    integer :: i, j, k, nt
    real :: rho, zmn 

    fout = 0.

    nt = 64*n_zer 

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
          !fout(i,j) = (sum(fmn(i,:)*zmn)-fmn(i,1)*zmn(1)/2&
          !            -fmn(i,ns)*zmn(ns)/2)/(ns-1)*(j+1)
       end do
    end do
    ! make sure the coefficients sum to the boundary value 
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
!    integer :: k 

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
          f(:,i+1) = fmn(:,i+1)*(1.*i/(ns-1))**(-xmv/2+1) 
       end do
   !    f(:,1)=fmn(:,1)
       do i = 1, mn_mode
          if(xmv(i).ne.0.) f(i,1) = 0.
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
       fout = fout*r**(xmv-2)
   end if 
!        s1 = 1.*js/(ns-1) 
!        s0 = 1.*(js-1)/(ns-1) 
!        rstc = (rmnc(:,js+1)*(1-ds)*s1**(-xmv/2+1)& 
!              + rmnc(:,js)*ds*s0**(-xmv/2+1))*r**(xmv-2) 
  end subroutine vmec_interpl
#endif
end module read_vmec 
