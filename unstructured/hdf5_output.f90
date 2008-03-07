module hdf5_output
  
  use hdf5

  implicit none

  integer(HID_T) :: file_id, offset, global_elms
  integer :: times_output
  logical :: initialized = .false.
  character(LEN=7), parameter :: hdf5_filename = "C1.h5"

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
  subroutine hdf5_initialize(error)
    use basic
    use hdf5

    implicit none

    include 'mpif.h'
   
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

    ! if irestart.eq.0 then create hdf5 file
    if(irestart.eq.0) then
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
    ! else open hdf5 file
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
    
    local_dims(1) = ndofs
    local_dims(2) = nelms
    global_dims(1) = ndofs
    global_dims(2) = global_elms
    off(1) = 0
    off(2) = offset
    
    ! Create global dataset
    call h5screate_simple_f(rank, global_dims, filespace, error)
    call h5dcreate_f(parent_id, name, H5T_NATIVE_REAL, filespace, &
         dset_id, error)
    call h5sclose_f(filespace, error)
    
    ! Select local hyperslab within dataset
    call h5screate_simple_f(rank, local_dims, memspace, error)
    call h5dget_space_f(dset_id, filespace, error)
    call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, off, local_dims, error)
    call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
  
    ! Write the dataset
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, values, global_dims, error, &
         file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

    ! Close HDF5 handles
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    
    call h5dclose_f(dset_id, error)
    call h5pclose_f(plist_id, error)

  end subroutine output_field

  ! output_scalar
  ! =============
  subroutine output_scalar(parent_id, name, value, time, error)
    use hdf5

    implicit none
    
    integer(HID_T), intent(in) :: parent_id
    character(LEN=*), intent(in) :: name
    real, intent(in) :: value
    integer, intent(in) :: time
    integer, intent(out) :: error

    integer(HSIZE_T) :: chunk_size(1) = (/ 100 /)
    integer(HSIZE_T) :: dims(1)
    integer(HSIZE_T) :: maxdims(1)
    integer(HSIZE_T) :: local_dims(1) = (/ 1 /)
    integer(HSIZE_T), dimension(1,1) :: coord
    real :: values(1)
    integer(HID_T) :: memspace, filespace, dset_id, p_id

    dims(1) = time+1
    maxdims(1) = H5S_UNLIMITED_F
    values(1) = value
    coord(1,1) = time + 1
    
    if(time.eq.0) then
       call h5screate_simple_f(1, dims, filespace, error, maxdims)
       call h5pcreate_f(H5P_DATASET_CREATE_F, p_id, error)
       call h5pset_chunk_f(p_id, 1, chunk_size, error)
       call h5dcreate_f(parent_id, name, H5T_NATIVE_REAL, &
            filespace, dset_id, error, p_id)
       call h5pclose_f(p_id, error)
       call h5sclose_f(filespace, error)
    else
       call h5dopen_f(parent_id, name, dset_id, error)
       call h5dextend_f(dset_id, dims, error)
    endif

    call h5screate_simple_f(1, local_dims, memspace, error)
    call h5dget_space_f(dset_id, filespace, error)
#ifdef _AIX
    call h5sselect_elements_f(filespace, H5S_SELECT_SET_F, 1, &
         1, coord, error)
#else
    call h5sselect_elements_f(filespace, H5S_SELECT_SET_F, 1, &
         local_dims(1), coord, error)
#endif
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, values, local_dims, error, &
         file_space_id = filespace, mem_space_id = memspace)

    ! Close HDF5 handles
    call h5sclose_f(filespace, error)
    call h5sclose_f(memspace, error)
    call h5dclose_f(dset_id, error)

  end subroutine output_scalar

end module hdf5_output


! hdf5_write_parameters
! =====================
subroutine hdf5_write_parameters(error)
  use hdf5
  use hdf5_output
  use basic

  implicit none

  integer, intent(out) :: error

  integer(HID_T) :: root_id

  call h5gopen_f(file_id, "/", root_id, error)

  call write_int_attr (root_id, "numvar"     , numvar,     error)
  call write_int_attr (root_id, "idens"      , idens,      error)
  call write_int_attr (root_id, "ipres"      , ipres,      error)
  call write_int_attr (root_id, "itor"       , itor,       error)
  call write_int_attr (root_id, "gyro"       , gyro,       error)
  call write_int_attr (root_id, "linear"     , linear,     error)
  call write_int_attr (root_id, "eqsubtract" , eqsubtract, error)
  call write_int_attr (root_id, "iper"       , iper,       error)
  call write_int_attr (root_id, "jper"       , jper,       error)
  call write_int_attr (root_id, "imask"      , imask,      error)
  call write_int_attr (root_id, "integrator" , integrator, error)
  call write_int_attr (root_id, "ipellet"    , ipellet,    error)
  call write_int_attr (root_id, "ivform"     , ivform,     error)
  call write_real_attr(root_id, "db"         , db,         error)
  call write_real_attr(root_id, "xzero"      , xzero,      error)
  call write_real_attr(root_id, "zzero"      , zzero,      error)
  call write_real_attr(root_id, "xlim"       , xlim,       error)
  call write_real_attr(root_id, "zlim"       , zlim,       error)
  call write_real_attr(root_id, "xmag"       , xmag,       error)
  call write_real_attr(root_id, "zmag"       , zmag,       error)
  call write_real_attr(root_id, "vloop"      , vloop,      error)
  call write_real_attr(root_id, "gam"        , gam,        error)
  call write_real_attr(root_id, "thimp"      , thimp,      error)
  call write_real_attr(root_id, "bzero"      , bzero,      error)
  call write_real_attr(root_id, "gravr"      , gravr,      error)
  call write_real_attr(root_id, "gravz"      , gravz,      error)
  call write_real_attr(root_id, "amu"        , amu,        error)
  call write_real_attr(root_id, "amuc"       , amuc,       error)
  call write_real_attr(root_id, "etar"       , etar,       error)
  call write_real_attr(root_id, "eta0"       , eta0,       error)
  call write_real_attr(root_id, "kappar"     , kappar,     error)
  call write_real_attr(root_id, "kappa0"     , kappa0,     error)
  call write_real_attr(root_id, "kappat"     , kappat,     error)
  call write_real_attr(root_id, "denm"       , denm,       error)
  call write_real_attr(root_id, "pellet_rate", pellet_rate,error)
  call write_real_attr(root_id, "pellet_var" , pellet_var, error)
  call write_real_attr(root_id, "pellet_x"   , pellet_x,   error)
  call write_real_attr(root_id, "pellet_z"   , pellet_z,   error)
  call write_real_attr(root_id, "ln"         , ln,         error)

  call h5gclose_f(root_id, error)

end subroutine hdf5_write_parameters

! hdf5_write_scalars
! ==================
subroutine hdf5_write_scalars(error)
  use basic
  use diagnostics
  use hdf5_output

  implicit none

  integer, intent(out) :: error

  integer(HID_T) :: root_id, scalar_group_id

  real :: temp

  if(myrank.eq.1) print *, 'writing scalars for ntime = ', ntime

  call h5gopen_f(file_id, "/", root_id, error)

  if(ntime.eq.0 .and. irestart.eq.0) then
     call h5gcreate_f(root_id, "scalars", scalar_group_id, error)
  else
     call h5gopen_f(root_id, "scalars", scalar_group_id, error)
  endif

  call output_scalar(scalar_group_id, "time" , time  , ntime, error)

  call output_scalar(scalar_group_id, "loop_voltage"    , vloop , ntime, error)
  call output_scalar(scalar_group_id, "psi_lcfs"        , psilim, ntime, error)

  call output_scalar(scalar_group_id, "area"            , area  , ntime, error)
  call output_scalar(scalar_group_id, "toroidal_flux"   , tflux , ntime, error)
  call output_scalar(scalar_group_id, "toroidal_current", totcur, ntime, error)
  call output_scalar(scalar_group_id, "particle_number" , totden, ntime, error)
  call output_scalar(scalar_group_id, "angular_momentum", tmom  , ntime, error)
  call output_scalar(scalar_group_id, "circulation"     , tvor  , ntime, error)

  call output_scalar(scalar_group_id, "area_p"            , parea, ntime, error)
  call output_scalar(scalar_group_id, "toroidal_flux_p"   , pflux, ntime, error)
  call output_scalar(scalar_group_id, "toroidal_current_p", pcur , ntime, error)
  call output_scalar(scalar_group_id, "particle_number_p" , pden , ntime, error)
  call output_scalar(scalar_group_id, "angular_momentum_p", pmom , ntime, error)
  call output_scalar(scalar_group_id, "circulation_p"     , pvor , ntime, error)
  
  call output_scalar(scalar_group_id, "E_MP" , emagp , ntime, error)
  call output_scalar(scalar_group_id, "E_KP" , ekinp , ntime, error)
  call output_scalar(scalar_group_id, "E_MPD", emagpd, ntime, error)
  call output_scalar(scalar_group_id, "E_KPD", ekinpd, ntime, error)
  call output_scalar(scalar_group_id, "E_MPH", emagph, ntime, error)
  call output_scalar(scalar_group_id, "E_KPH", ekinph, ntime, error)

  if(numvar.ge.2) then
     call output_scalar(scalar_group_id, "E_MT" , emagt , ntime, error)
     call output_scalar(scalar_group_id, "E_KT" , ekint , ntime, error)
     call output_scalar(scalar_group_id, "E_MTD", emagtd, ntime, error)
     call output_scalar(scalar_group_id, "E_KTD", ekintd, ntime, error)
     call output_scalar(scalar_group_id, "E_MTH", emagth, ntime, error)
     call output_scalar(scalar_group_id, "E_KTH", ekinth, ntime, error)
  endif

  if(numvar.ge.3) then
     call output_scalar(scalar_group_id, "E_P" , emag3, ntime, error)
     call output_scalar(scalar_group_id, "E_K3", ekin3, ntime, error)
     call output_scalar(scalar_group_id, "E_PD", emag3d, ntime, error)
     call output_scalar(scalar_group_id, "E_K3D", ekin3d, ntime, error)
     call output_scalar(scalar_group_id, "E_PH", emag3h, ntime, error)
     call output_scalar(scalar_group_id, "E_K3H", ekin3h, ntime, error)
  endif

  call output_scalar(scalar_group_id, "Flux_pressure ", efluxp, ntime, error)
  call output_scalar(scalar_group_id, "Flux_kinetic  ", efluxk, ntime, error)
  call output_scalar(scalar_group_id, "Flux_poynting ", efluxs, ntime, error)
  call output_scalar(scalar_group_id, "Flux_thermal  ", efluxt, ntime, error)
  call output_scalar(scalar_group_id, "E_grav        ", epotg,  ntime, error)

  call output_scalar(scalar_group_id, "Particle_Flux_diffusive", &
       nfluxd, ntime, error)
  call output_scalar(scalar_group_id, "Particle_Flux_convective", &
       nfluxv, ntime, error)
  call output_scalar(scalar_group_id, "Particle_source", &
       nsource, ntime, error)

  call output_scalar(scalar_group_id, "Torque_em",   tau_em,   ntime, error)
  call output_scalar(scalar_group_id, "Torque_sol",  tau_sol,  ntime, error)
  call output_scalar(scalar_group_id, "Torque_com",  tau_com,  ntime, error)
  call output_scalar(scalar_group_id, "Torque_visc", tau_visc, ntime, error)
  call output_scalar(scalar_group_id, "Torque_gyro", tau_gyro, ntime, error)
  call output_scalar(scalar_group_id, "Torque_denm", tau_denm, ntime, error)


  if(itaylor.eq.3) then
     temp = reconnected_flux()
     call output_scalar(scalar_group_id, "Reconnected_Flux", temp, ntime, error)
  endif

  call h5gclose_f(scalar_group_id, error)
  call h5gclose_f(root_id, error)

end subroutine hdf5_write_scalars


! hdf5_write_timings
! ==================
subroutine hdf5_write_timings(error)
  use basic
  use diagnostics
  use hdf5_output

  implicit none

  integer, intent(out) :: error

  integer(HID_T) :: root_id, timing_group_id

  if(maxrank.gt.1) call distribute_timings

  call h5gopen_f(file_id, "/", root_id, error)

  if(ntime.eq.0 .and. irestart.eq.0) then
     call h5gcreate_f(root_id, "timings", timing_group_id, error)

     ! for grad-shafranov equilibrium, output gs times
     if(itor.eq.1 .and. itaylor.eq.1) then
        call write_real_attr(timing_group_id, "t_gs", t_gs, error) 
        call write_real_attr(timing_group_id, "t_gs_magaxis", t_gs_magaxis, error) 
        call write_real_attr(timing_group_id, "t_gs_fundef", t_gs_fundef, error)
        call write_real_attr(timing_group_id, "t_gs_solve", t_gs_solve, error) 
        call write_real_attr(timing_group_id, "t_gs_init", t_gs_init, error) 
     endif
  else
     call h5gopen_f(root_id, "timings", timing_group_id, error)
  endif

  call output_scalar(timing_group_id, "t_ludefall"    , t_ludefall    , ntime, error)
  call output_scalar(timing_group_id, "t_sources"     , t_sources     , ntime, error)
  call output_scalar(timing_group_id, "t_smoother"    , t_smoother    , ntime, error)
  call output_scalar(timing_group_id, "t_aux"         , t_aux         , ntime, error)
  call output_scalar(timing_group_id, "t_solve_v"     , t_solve_v     , ntime, error)
  call output_scalar(timing_group_id, "t_solve_b"     , t_solve_b     , ntime, error)
  call output_scalar(timing_group_id, "t_solve_n"     , t_solve_n     , ntime, error)
  call output_scalar(timing_group_id, "t_solve_p"     , t_solve_p     , ntime, error)
  call output_scalar(timing_group_id, "t_output_cgm"  , t_output_cgm  , ntime, error)
  call output_scalar(timing_group_id, "t_output_hdf5" , t_output_hdf5 , ntime, error)
  call output_scalar(timing_group_id, "t_output_reset", t_output_reset, ntime, error)
  call output_scalar(timing_group_id, "t_mvm"         , t_mvm         , ntime, error)
  call output_scalar(timing_group_id, "t_onestep"     , t_onestep     , ntime, error)

  call h5gclose_f(timing_group_id, error)
  call h5gclose_f(root_id, error)

end subroutine hdf5_write_timings


! hdf5_write_time_slice
! =====================
subroutine hdf5_write_time_slice(equilibrium, error)
  use hdf5
  use hdf5_output
  use basic

  implicit none
  
  include 'mpif.h'

  integer, intent(out) :: error
  integer, intent(in) :: equilibrium

  character(LEN=19) :: time_group_name
  integer(HID_T) :: time_group_id, root_id
  integer :: nelms

!  integer :: global_nodes, global_edges, global_regions
   
  call numfac(nelms)

  ! Calculate offset of current process
  call mpi_scan(nelms, offset, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, error)
  offset = offset - nelms
!  print *, "Offset of ", myrank, " = ", offset
!  call numglobalents(global_nodes, gobal_edges, global_elms, global_regions)
  call mpi_allreduce(nelms, global_elms, 1, MPI_INTEGER, &
       MPI_SUM, MPI_COMM_WORLD, error)


  ! Create the time group
  ! ~~~~~~~~~~~~~~~~~~~~~
  
  ! create the name of the group
  if(equilibrium.eq.1) then
     time_group_name = "equilibrium"
  else
     write(time_group_name, '("time_",I3.3)') times_output
     ! remove the time group if it already exists
     ! (from before a restart, for example)
     if(irestart.eq.1) then
        call h5gunlink_f(file_id, time_group_name, error)
     endif
  endif

  if(myrank.eq.1 .and. iprint.eq.1) &
       print *, 'Writing time slice ', time_group_name

  ! create the group
  call h5gcreate_f(file_id, time_group_name, time_group_id, error)

  ! Write attributes
  ! ~~~~~~~~~~~~~~~~
  call write_real_attr(time_group_id, "time", time, error)
  call write_int_attr(time_group_id, "nspace", 2, error)

  ! Output the mesh data
  call output_mesh(time_group_id, nelms, error)

  ! Output the field data 
  call output_fields(time_group_id, equilibrium, error)


  ! Close the time group
  ! ~~~~~~~~~~~~~~~~~~~~
  call h5gclose_f(time_group_id, error)  

  if(equilibrium.eq.0) times_output = times_output + 1
  call h5gopen_f(file_id, "/", root_id, error)
  call update_int_attr(root_id, "ntime", times_output, error)
  call h5gclose_f(root_id, error)

end subroutine hdf5_write_time_slice


! output_mesh
! ===========
subroutine output_mesh(time_group_id, nelms, error)
  use hdf5
  use hdf5_output
  use t_data

  implicit none

  integer(HID_T), intent(in) :: time_group_id
  integer, intent(in) :: nelms
  integer, intent(out) :: error

  integer(HID_T) :: mesh_group_id
  integer :: i
  real, dimension(6,nelms) :: elm_data
  double precision, dimension(3) :: coords
  integer, dimension(4) :: nodeids
  real :: alx, alz

  ! Create the group
  call h5gcreate_f(time_group_id, "mesh", mesh_group_id, error) 

  ! Write attributes
  call write_int_attr(mesh_group_id, "nelms", global_elms, error)
  call getboundingboxsize(alx, alz)
  call write_real_attr(mesh_group_id, "width", alx, error)
  call write_real_attr(mesh_group_id, "height", alz, error)

  ! Output the mesh data
  do i=1, nelms
     call nodfac(i,nodeids)
     call xyznod(nodeids(1), coords)

     elm_data(1,i) = atri(i)
     elm_data(2,i) = btri(i)
     elm_data(3,i) = ctri(i)
     elm_data(4,i) = ttri(i)
     elm_data(5,i) = coords(1)
     elm_data(6,i) = coords(2)
  end do
  call output_field(mesh_group_id, "elements", elm_data, 6, nelms, error)

  ! Close the group
  call h5gclose_f(mesh_group_id, error)
end subroutine output_mesh


! output_fields
! =============
subroutine output_fields(time_group_id, equilibrium, error)
  use hdf5
  use hdf5_output
  use basic
  use arrays
  
  implicit none

  integer(HID_T), intent(in) :: time_group_id
  integer, intent(out) :: error
  integer, intent(in) :: equilibrium
  
  integer(HID_T) :: group_id
  integer :: i, nelms, nfields
  vectype, allocatable :: dum(:,:)
  vectype, pointer :: f_ptr(:)

  if(equilibrium.eq.1) then
     f_ptr => field0
  else
     f_ptr => field
  endif

  nfields = 0
  call numfac(nelms)
  
  allocate(dum(20,nelms))

  ! Create the fields group
  call h5gcreate_f(time_group_id, "fields", group_id, error)

  ! Output the fields
  ! ~~~~~~~~~~~~~~~~~
  
  ! psi
  do i=1, nelms
     call calcavector(i, f_ptr, psi_g, num_fields, dum(:,i))
  end do
  call output_field(group_id, "psi", real(dum), 20, nelms, error)
  nfields = nfields + 1

  ! u
  do i=1, nelms
     call calcavector(i, f_ptr, u_g, num_fields, dum(:,i))
  end do
  call output_field(group_id, "phi", real(dum), 20, nelms, error)
  nfields = nfields + 1

  ! jphi
  do i=1, nelms
     call calcavector(i, jphi, 1, 1, dum(:,i))
  end do
  call output_field(group_id, "jphi", real(dum), 20, nelms, error)
  nfields = nfields + 1

  ! vor
  do i=1, nelms
     call calcavector(i, vor, 1, 1, dum(:,i))
  end do
  call output_field(group_id, "vor", real(dum), 20, nelms, error)
  nfields = nfields + 1

  ! eta
  do i=1, nelms
     call calcavector(i, resistivity, 1, 1, dum(:,i))
  end do
  call output_field(group_id, "eta", real(dum), 20, nelms, error)
  nfields = nfields + 1

  ! visc
  do i=1, nelms
     call calcavector(i, visc, 1, 1, dum(:,i))
  end do
  call output_field(group_id, "visc", real(dum), 20, nelms, error)
  nfields = nfields + 1

  ! tempvar
  do i=1, nelms
     call calcavector(i, tempvar, 1, 1, dum(:,i))
  end do
  call output_field(group_id, "tempvar", real(dum), 20, nelms, error)
  nfields = nfields + 1

  ! I
  do i=1, nelms
     call calcavector(i, f_ptr, bz_g, num_fields, dum(:,i))
  end do
  call output_field(group_id, "I", real(dum), 20, nelms, error)
  nfields = nfields + 1
  
  ! BF
  if(i3d.eq.1) then
     do i=1, nelms
        call calcavector(i, bf, 1, 1, dum(:,i))
     end do
     call output_field(group_id, "f", real(dum), 20, nelms, error)
     nfields = nfields + 1
  endif

  ! V
  do i=1, nelms
     call calcavector(i, f_ptr, vz_g, num_fields, dum(:,i))
  end do
  call output_field(group_id, "V", real(dum), 20, nelms, error)
  nfields = nfields + 1

  ! P and Pe
  if(ipres.eq.1) then
     do i=1, nelms
        call calcavector(i, f_ptr, pe_g, num_fields, dum(:,i))
     end do
     call output_field(group_id, "Pe", real(dum), 20, nelms, error)
     nfields = nfields + 1
     do i=1, nelms
        call calcavector(i, f_ptr, p_g, num_fields, dum(:,i))
     end do
     call output_field(group_id, "P", real(dum), 20, nelms, error)
     nfields = nfields + 1
  else
     do i=1, nelms
        call calcavector(i, f_ptr, pe_g, num_fields, dum(:,i))
     end do
     call output_field(group_id, "Pe", pefac*real(dum), 20, nelms, error)
     nfields = nfields + 1
     call output_field(group_id, "P", real(dum), 20, nelms, error)
     nfields = nfields + 1
  endif
     
  ! chi
  do i=1, nelms
     call calcavector(i, f_ptr, chi_g, num_fields, dum(:,i))
  end do
  call output_field(group_id, "chi", real(dum), 20, nelms, error)
  nfields = nfields + 1

  ! com
  do i=1, nelms
     call calcavector(i, com, 1, 1, dum(:,i))
  end do
  call output_field(group_id, "com", real(dum), 20, nelms, error)
  nfields = nfields + 1

  ! kappa
  do i=1, nelms
     call calcavector(i, kappa, 1, 1, dum(:,i))
  end do
  call output_field(group_id, "kappa", real(dum), 20, nelms, error)
  nfields = nfields + 1

  ! visc_c
  do i=1, nelms
     call calcavector(i, visc_c, 1, 1, dum(:,i))
  end do
  call output_field(group_id, "visc_c", real(dum), 20, nelms, error)
  nfields = nfields + 1

  ! den
  do i=1, nelms
     call calcavector(i, f_ptr, den_g, num_fields, dum(:,i))
  end do
  call output_field(group_id, "den", real(dum), 20, nelms, error)
  nfields = nfields + 1
  
  if(ipellet.eq.1 .or. ionization.eq.1) then
     do i=1, nelms
        call calcavector(i, sigma, 1, 1, dum(:,i))
     end do
     call output_field(group_id, "sigma", real(dum), 20, nelms, error)
     nfields = nfields + 1
  endif

  if(isources.eq.1) then
     do i=1, nelms
        call calcavector(i, sb1, 1, 1, dum(:,i))
     end do
     call output_field(group_id, "sb1", real(dum), 20, nelms, error)
     nfields = nfields + 1
     if(numvar.ge.2) then
        do i=1, nelms
           call calcavector(i, sb2, 1, 1, dum(:,i))
        end do
        call output_field(group_id, "sb2", real(dum), 20, nelms, error)
        nfields = nfields + 1
     endif
     if(numvar.ge.3) then
        do i=1, nelms
           call calcavector(i, sp1, 1, 1, dum(:,i))
        end do
        call output_field(group_id, "sp1", real(dum), 20, nelms, error)
        nfields = nfields + 1
     endif
  endif

  if(gyro.eq.1) then
     ! gyro_tau
     do i=1, nelms
        call calcavector(i, gyro_tau, 1, 1, dum(:,i))
     end do
     call output_field(group_id, "gyro_tau", real(dum), 20, nelms, error)
     nfields = nfields + 1
  end if


  call write_int_attr(group_id, "nfields", nfields, error)

  ! Close the mesh group
  call h5gclose_f(group_id, error)

  deallocate(dum)
end subroutine output_fields
