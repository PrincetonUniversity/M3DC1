module hdf5_output
  
  use hdf5

  implicit none

  integer(HID_T) :: file_id
  integer :: offset, global_elms
  integer :: times_output
  logical, private :: initialized = .false.
  character(LEN=7), parameter, private :: hdf5_filename = "C1.h5"

  integer :: idouble_out

contains

  ! hdf5_flush
  ! ==========
  subroutine hdf5_flush(error)
    use hdf5
  
    implicit none

    integer, intent(out) :: error

    ! Flush the data to disk
    call h5fflush_f(file_id, H5F_SCOPE_GLOBAL_F, error)
  end subroutine hdf5_flush


  ! hdf5_initialize
  ! ===============
  subroutine hdf5_initialize(restart, error)
    use hdf5

    implicit none

    include 'mpif.h'

    logical, intent(in) :: restart     ! if true, do not overwrite file
    integer, intent(out) :: error

    integer(HID_T) :: root_id, plist_id
    integer :: info

    call h5open_f(error)
    if(error.lt.0) then
       print *, "Error: could not initialize HDF5 library: ", error
       return
    endif

    ! Set up the file access property list with parallel I/O
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    info = MPI_INFO_NULL
    call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, info, error)

    if(.not.restart) then
       ! create hdf5 file
       call h5fcreate_f(hdf5_filename, H5F_ACC_TRUNC_F, file_id, error, &
            access_prp = plist_id)
       if(error.lt.0) then
          print *, "Error: could not open ", hdf5_filename, &
               " for HDF5 output.  error = ", error
          return
       endif

       call h5gopen_f(file_id, "/", root_id, error)
       call write_int_attr(root_id, "ntime", 0, error)
       call h5gclose_f(root_id, error)

       times_output = 0
    else
       ! open hdf5 file
       call h5fopen_f(hdf5_filename, H5F_ACC_RDWR_F, file_id, error, &
            access_prp = plist_id)
       if(error.lt.0) then
          print *, "Error: could not open ", &
               hdf5_filename, " for HDF5 output: ", error
       endif

       call h5gopen_f(file_id, "/", root_id, error)
       call read_int_attr(root_id, "ntime", times_output, error)
       call h5gclose_f(root_id, error)

       ! overwrite the last time slice
       times_output = times_output - 1
    endif

    call h5pclose_f(plist_id, error)

    initialized = .true.
  end subroutine hdf5_initialize


  ! hdf5_finalize
  ! =============
  subroutine hdf5_finalize(error)
    use hdf5
    
    implicit none
    
    integer, intent(out) :: error
    
    error = 0

    if(.not. initialized) return

    ! Close the file.
    call h5fclose_f(file_id, error)
    if(error .lt. 0) print *, "Error closing hdf5 file"

    call h5close_f(error)
    if(error .lt. 0) print *, "Error closing hdf5 library"
    
  end subroutine hdf5_finalize

  ! read_int_attr
  ! =============
  subroutine read_int_attr(parent_id, name, value, error)
    use hdf5
    
    implicit none
    
    integer(HID_T), intent(in) :: parent_id
    character(LEN=*), intent(in) :: name
    integer, intent(out)  :: value
    integer, intent(out) :: error
    
    integer(HID_T) :: attr_id
    integer(HSIZE_T), dimension(1) :: dims = 1

    call h5aopen_name_f(parent_id, name, attr_id, error)
    call h5aread_f(attr_id, H5T_NATIVE_INTEGER, value, dims, error)
    call h5aclose_f(attr_id, error)
    
  end subroutine read_int_attr

  ! write_int_attr
  ! ==============
  subroutine write_int_attr(parent_id, name, value, error)
    use hdf5
    
    implicit none
    
    integer(HID_T), intent(in) :: parent_id
    character(LEN=*), intent(in) :: name
    integer, intent(in)  :: value
    integer, intent(out) :: error
    
    integer(HID_T) :: dspace_id, attr_id
    integer(HSIZE_T), dimension(1) :: dims = 1

    call h5screate_f(H5S_SCALAR_F, dspace_id, error)
    call h5acreate_f(parent_id, name, H5T_NATIVE_INTEGER, dspace_id, attr_id, error)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, value, dims, error)
    call h5aclose_f(attr_id, error)
    call h5sclose_f(dspace_id, error)
    
  end subroutine write_int_attr

  ! update_int_attr
  ! ===============
  subroutine update_int_attr(parent_id, name, value, error)
    use hdf5
    
    implicit none
    
    integer(HID_T), intent(in) :: parent_id
    character(LEN=*), intent(in) :: name
    integer, intent(in)  :: value
    integer, intent(out) :: error
    
    integer(HID_T) :: dspace_id, attr_id
    integer(HSIZE_T), dimension(1) :: dims = 1

    call h5screate_f(H5S_SCALAR_F, dspace_id, error)
    
    call h5aopen_name_f(parent_id, name, attr_id, error)
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, value, dims, error)
    call h5aclose_f(attr_id, error)
    call h5sclose_f(dspace_id, error)
    
  end subroutine update_int_attr

  ! write_real_attr
  ! ===============
  subroutine write_real_attr(parent_id, name, value, error)
    use hdf5
    
    implicit none
    
    integer(HID_T), intent(in) :: parent_id
    character(LEN=*), intent(in) :: name
    real, intent(in)  :: value
    integer, intent(out) :: error
    
    integer(HID_T) :: dspace_id, attr_id
    integer(HSIZE_T), dimension(1) :: dims = 1

    call h5screate_f(H5S_SCALAR_F, dspace_id, error)
    
    call h5acreate_f(parent_id, name, H5T_NATIVE_DOUBLE, dspace_id, attr_id, &
         error)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, value, dims, error)
    call h5aclose_f(attr_id, error)
    call h5sclose_f(dspace_id, error)
    
  end subroutine write_real_attr


  ! write_vec_attr
  ! ==============
  subroutine write_vec_attr(parent_id, name, values, len, error)
    use hdf5
    
    implicit none
    
    integer(HID_T), intent(in) :: parent_id
    character(LEN=*), intent(in) :: name
    integer, intent(in) :: len
    real, dimension(len), intent(in)  :: values
    integer, intent(out) :: error
    
    integer(HID_T) :: dspace_id, attr_id
    integer(HSIZE_T), dimension(1) :: dims
#ifdef USETAU
    integer :: dummy     ! this is necessary to prevent TAU from
    dummy = 0            ! breaking formatting requirements
#endif

    dims(1) = len

    call h5screate_simple_f(1, dims, dspace_id, error)
    call h5acreate_f(parent_id, name, H5T_NATIVE_DOUBLE, dspace_id, attr_id, error)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, values, dims, error)
    call h5aclose_f(attr_id, error)
    call h5sclose_f(dspace_id, error)
    
  end subroutine write_vec_attr


  ! output_field
  ! ============
  subroutine output_field(parent_id, name, values, ndofs, nelms, error)
    use hdf5
    
    implicit none
    
    integer(HID_T), intent(in) :: parent_id
    character(LEN=*), intent(in) :: name
    integer, intent(in) :: ndofs, nelms
    real, dimension(ndofs, nelms), intent(in) :: values
    integer, intent(out) :: error
    
    integer, parameter ::  rank = 2
    integer(HID_T) :: filespace, memspace, dset_id, plist_id
    integer(HSIZE_T), dimension(rank) :: local_dims, global_dims
    integer(HSSIZE_T), dimension(rank) :: off

#ifdef USETAU
    integer :: dummy     ! this is necessary to prevent TAU from
    dummy = 0            ! breaking formatting requirements
#endif
    
    local_dims(1) = ndofs
    local_dims(2) = nelms
    global_dims(1) = ndofs
    global_dims(2) = global_elms
    off(1) = 0
    off(2) = offset

    ! Create global dataset
    call h5screate_simple_f(rank, global_dims, filespace, error)
         if(error.ne.0) then
           write(*,*) error,rank," after h5screate_simple_f"
           call safestop(101)
         endif
    if(idouble_out.eq.1) then
       call h5dcreate_f(parent_id, name, H5T_NATIVE_DOUBLE, filespace, &
            dset_id, error)
    else
       call h5dcreate_f(parent_id, name, H5T_NATIVE_REAL, filespace, &
            dset_id, error)
    end if
         if(error.ne.0) then
           write(*,*) error,rank," after h5dcreate_f"
           call safestop(101)
         endif
    call h5sclose_f(filespace, error)
         if(error.ne.0) then
           write(*,*) error,rank," after hsclose_f"
           call safestop(101)
         endif
    
    ! Select local hyperslab within dataset
    call h5screate_simple_f(rank, local_dims, memspace, error)
         if(error.ne.0) then
           write(*,*) error,rank," after h5screate_simple_f"
           call safestop(102)
         endif
    call h5dget_space_f(dset_id, filespace, error)
         if(error.ne.0) then
           write(*,*) error,rank," after h5dget_space_f"
           call safestop(102)
         endif
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, off, local_dims, &
         error)
         if(error.ne.0) then
           write(*,*) error,rank," after h5sselect_hyperslab_f"
           call safestop(102)
         endif
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
         if(error.ne.0) then
           write(*,*) error,rank," after h5pcreate_f"
           call safestop(102)
         endif
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
         if(error.ne.0) then
           write(*,*) error,rank," after h5pset_dxpl_mpio_f"
           call safestop(102)
         endif
  
    ! Write the dataset
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, values, global_dims, error, &
         file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
         if(error.ne.0) then
           write(*,*) error,rank," h5dwrite_f"
           call safestop(103)
         endif

    ! Close HDF5 handles
    call h5sclose_f(filespace, error)
         if(error.ne.0) then
           write(*,*) error,rank," h5sclose_f"
           call safestop(104)
         endif
    call h5sclose_f(memspace, error)
         if(error.ne.0) then
           write(*,*) error,rank," h5sclose_f"
           call safestop(105)
         endif
    
    call h5dclose_f(dset_id, error)
         if(error.ne.0) then
           write(*,*) error,rank," h5dclose_f"
           call safestop(104)
         endif
    call h5pclose_f(plist_id, error)
         if(error.ne.0) then
           write(*,*) error,rank," h5dclose_f"
           call safestop(105)
         endif

  end subroutine output_field

  ! output_scalar
  ! =============
  subroutine output_scalar(parent_id, name, value, t, error)
    use basic
    use hdf5

    implicit none
    
    integer(HID_T), intent(in) :: parent_id
    character(LEN=*), intent(in) :: name
    real, intent(in) :: value
    integer, intent(in) :: t
    integer, intent(out) :: error

    integer(HSIZE_T) :: chunk_size(1) = (/ 100 /)
    integer(HSIZE_T) :: dims(1)
    integer(HSIZE_T) :: maxdims(1)
    integer(HSIZE_T), parameter :: local_dims(1) = (/ 1 /)
    integer(HSIZE_T), dimension(1,1) :: coord
    integer(SIZE_T), parameter :: num_elements = 1
    integer(HID_T) :: memspace, filespace, dset_id, p_id, plist_id
    real :: values(1)

#ifdef USETAU
    integer :: dummy     ! this is necessary to prevent TAU from
    dummy = 0            ! breaking formatting requirements
#endif

    dims(1) = t+1
    maxdims(1) = H5S_UNLIMITED_F
    values(1) = value
    coord(1,1) = t + 1
    
    if(t.eq.0) then
       call h5screate_simple_f(1, dims, filespace, error, maxdims)
       call h5pcreate_f(H5P_DATASET_CREATE_F, p_id, error)
       call h5pset_chunk_f(p_id, 1, chunk_size, error)
       if(idouble_out.eq.1) then
          call h5dcreate_f(parent_id, name, H5T_NATIVE_DOUBLE, &
               filespace, dset_id, error, p_id)
       else
          call h5dcreate_f(parent_id, name, H5T_NATIVE_REAL, &
               filespace, dset_id, error, p_id)
       end if
       call h5pclose_f(p_id, error)
       call h5sclose_f(filespace, error)
    else
       call h5dopen_f(parent_id, name, dset_id, error)
       call h5dextend_f(dset_id, dims, error)
    endif

    if(myrank.eq.0) then
       call h5screate_simple_f(1, local_dims, memspace, error)
       call h5dget_space_f(dset_id, filespace, error)
       call h5sselect_elements_f(filespace, H5S_SELECT_SET_F, 1, &
            num_elements, coord, error)
       
       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)

       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, values, local_dims, error, &
            file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)

       ! Close HDF5 handles
       call h5pclose_f(plist_id, error)
       call h5sclose_f(filespace, error)
       call h5sclose_f(memspace, error)
    endif
    call h5dclose_f(dset_id, error)

  end subroutine output_scalar


  ! output_1dvector
  ! ============
  subroutine output_1darr(parent_id, name, values, ndim, t, error)
    use hdf5

    implicit none
    
    integer(HID_T), intent(in) :: parent_id
    character(LEN=*), intent(in) :: name
    integer, intent(in) :: ndim, t
    real, dimension(ndim), intent(in) :: values
    integer, intent(out) :: error
    
    integer, parameter ::  rank = 1
    integer(HID_T) :: filespace, memspace, dset_id, plist_id
    integer(HSIZE_T), dimension(rank) :: local_dims, global_dims
    integer(HSSIZE_T), dimension(rank) :: off

#ifdef USETAU
    integer :: dummy     ! this is necessary to prevent TAU from
    dummy = 0            ! breaking formatting requirements
#endif
    
    local_dims(1) = ndim
    global_dims(1) = ndim
    off(1) = 0
    
    ! Create global dataset
    call h5screate_simple_f(rank, global_dims, filespace, error)
         if(error.ne.0) then
           write(*,*) error,rank," after h5sceate_simple_f"
           call safestop(101)
         endif
    if(idouble_out.eq.1) then 
       call h5dcreate_f(parent_id, name, H5T_NATIVE_DOUBLE, filespace, &
            dset_id, error)
    else
       call h5dcreate_f(parent_id, name, H5T_NATIVE_REAL, filespace, &
            dset_id, error)
    end if
         if(error.ne.0) then
           write(*,*) error,rank," after h5dcreate_f"
           call safestop(101)
         endif
    call h5sclose_f(filespace, error)
         if(error.ne.0) then
           write(*,*) error,rank," after hsclose_f"
           call safestop(101)
         endif
    
    ! Select local hyperslab within dataset
    call h5screate_simple_f(rank, local_dims, memspace, error)
         if(error.ne.0) then
           write(*,*) error,rank," after h5screate_simple_f"
           call safestop(102)
         endif
    call h5dget_space_f(dset_id, filespace, error)
         if(error.ne.0) then
           write(*,*) error,rank," after h5dget_space_f"
           call safestop(102)
         endif
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, off, local_dims, &
         error)
         if(error.ne.0) then
           write(*,*) error,rank," after h5sselect_hyperslab_f"
           call safestop(102)
         endif
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
         if(error.ne.0) then
           write(*,*) error,rank," after h5pcreate_f"
           call safestop(102)
         endif
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
         if(error.ne.0) then
           write(*,*) error,rank," after h5pset_dxpl_mpio_f"
           call safestop(102)
         endif
  
    ! Write the dataset
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, values, global_dims, error, &
         file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
         if(error.ne.0) then
           write(*,*) error,rank," h5dwrite_f"
           call safestop(103)
         endif

    ! Close HDF5 handles
    call h5sclose_f(filespace, error)
         if(error.ne.0) then
           write(*,*) error,rank," h5sclose_f"
           call safestop(104)
         endif
    call h5sclose_f(memspace, error)
         if(error.ne.0) then
           write(*,*) error,rank," h5sclose_f"
           call safestop(105)
         endif
    
    call h5dclose_f(dset_id, error)
         if(error.ne.0) then
           write(*,*) error,rank," h5dclose_f"
           call safestop(104)
         endif
    call h5pclose_f(plist_id, error)
         if(error.ne.0) then
           write(*,*) error,rank," h5dclose_f"
           call safestop(105)
         endif

  end subroutine output_1darr


  ! output_1dextendarr
  ! =============
  subroutine output_1dextendarr(parent_id, name, value, NMAX, t, error)
    use basic
    use hdf5

    implicit none
    
    integer(HID_T), intent(in) :: parent_id
    character(LEN=*), intent(in) :: name
    integer, intent(in) :: NMAX, t
    real, intent(in) :: value(NMAX)
    integer, intent(out) :: error

    integer, parameter ::  rank = 2
    integer(HSIZE_T) :: chunk_size(2)
    integer(HSIZE_T) :: dims(2), maxdims(2), local_dims(2), off(2)
    integer(SIZE_T) :: num_elements
    integer(HID_T) :: memspace, filespace, dset_id, p_id, plist_id

    integer :: j

#ifdef USETAU
    integer :: dummy     ! this is necessary to prevent TAU from
    dummy = 0            ! breaking formatting requirements
#endif

    dims(1) = NMAX
    dims(2) = t+1
    maxdims(1) = NMAX
    maxdims(2) = H5S_UNLIMITED_F
    local_dims(1) = NMAX
    local_dims(2) = 1
    chunk_size(1) = NMAX
    chunk_size(2) = 1
    off(1) = 0
    off(2) = t
    num_elements = NMAX
    
    if(t.eq.0) then
       call h5screate_simple_f(rank, dims, filespace, error, maxdims)
       call h5pcreate_f(H5P_DATASET_CREATE_F, p_id, error)
       call h5pset_chunk_f(p_id, rank, chunk_size, error)
       if(idouble_out.eq.1) then
          call h5dcreate_f(parent_id, name, H5T_NATIVE_DOUBLE, &
               filespace, dset_id, error, p_id)
       else
          call h5dcreate_f(parent_id, name, H5T_NATIVE_REAL, &
               filespace, dset_id, error, p_id)
       end if
       call h5pclose_f(p_id, error)
       call h5sclose_f(filespace, error)
    else
       call h5dopen_f(parent_id, name, dset_id, error)
       call h5dextend_f(dset_id, dims, error)
    endif

    if(myrank.eq.0) then
       call h5screate_simple_f(rank, local_dims, memspace, error)
       call h5dget_space_f(dset_id, filespace, error)
       call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, off, local_dims, &
            error)
       
       call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
       call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)

       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, value, local_dims, error, &
            file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)

       ! Close HDF5 handles
       call h5pclose_f(plist_id, error)
       call h5sclose_f(filespace, error)
       call h5sclose_f(memspace, error)
    endif
    call h5dclose_f(dset_id, error)

  end subroutine output_1dextendarr


end module hdf5_output
