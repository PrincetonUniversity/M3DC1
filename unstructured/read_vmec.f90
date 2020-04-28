! This module reads in data from VMEC files
module read_vmec 
  implicit none

#ifdef USEST
  integer :: nfp
  real, allocatable :: rbc(:), zbs(:)
  real, allocatable :: rstc(:), zsts(:)
  real, allocatable :: rmnc(:,:), zmns(:,:)
  real, allocatable :: bsupumnc(:,:), bsupvmnc(:,:)
  real, allocatable :: presf(:)
!  real, allocatable :: raxiscc(:), zaxiscs(:)
  integer, allocatable :: mb(:), nb(:)
  real, allocatable :: xmv(:), xnv(:)
  integer :: mn_mode, ns, n_tor 

contains

  subroutine read_vmec_h5(myrank)
    use hdf5
    
    implicit none
    
    character(len=256) :: vmec_filename = "geometry.h5"

    integer, intent(in) :: myrank 
    integer :: error
    
    integer(HID_T) :: file_id, dset_id, attr_id
    integer(HSIZE_T), dimension(1) :: dim0 = 1
    integer(HSIZE_T), dimension(1) :: dim1 
    integer(HSIZE_T), dimension(2) :: dim2

    call h5open_f(error)
    call h5fopen_f(vmec_filename, H5F_ACC_RDONLY_F, file_id, error)

    ! read attributes ns, n_tor, nfp, and mn_mode
    call h5gopen_f(file_id, '/', dset_id, error)
    call h5aopen_name_f(dset_id, 'ns', attr_id, error)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, ns, dim0, error)
    call h5aclose_f(attr_id, error)
    call h5aopen_name_f(dset_id, 'ntor', attr_id, error)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, n_tor, dim0, error)
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

    dim1(1) = ns
    allocate(presf(ns))
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
    allocate(rstc(mn_mode))
    allocate(zsts(mn_mode))
    rbc = rmnc(:,ns)
    zbs = zmns(:,ns)

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

#endif
end module read_vmec 
