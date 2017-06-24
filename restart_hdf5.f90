module restart_hdf5
  implicit none

contains
  
  subroutine rdrestart_hdf5()
    use basic
    use hdf5_output
    use hdf5
    use pellet

    implicit none

    integer :: error
    integer(HID_T) :: root_id, scalar_group_id, time_id, eq_time_id
    character(LEN=19) :: time_group_name

    integer :: times_output_in, icomplex_in

    if(myrank.eq.0) print *, 'Reading HDF5 file for restart.'

    call h5gopen_f(file_id, "/", root_id, error)

    call read_int_attr(root_id, "icomplex", icomplex_in, error)
    if(icomplex.ne.icomplex_in) then
       if(myrank.eq.0) print *, 'Error: restart conflict with icomplex'
       call safestop(2)
    end if

    call read_int_attr(root_id, "version", version_in, error)
    if(version_in.lt.16) then
       if(myrank.eq.0) print *, 'Error: HDF5 file is from too old a version to use iread_hdf5=1.'
       call safestop(1)
    end if

    ! Read Time Slice
    if(myrank.eq.0) print *, 'Reading data from time slice', times_output
    call read_int_attr(root_id, "ntime", times_output_in, error)
    if(times_output_in .le. times_output) then
       print *, 'Error: Requested restart time slice exceeds number of time slices in file.'
       call safestop(2)
    end if

    write(time_group_name, '("time_",I3.3)') times_output
    call h5gopen_f(root_id, time_group_name, time_id, error)
    call read_int_attr(time_id, "ntimestep", ntime, error)

    ! Read Scalars
    call h5gopen_f(root_id, "scalars", scalar_group_id, error)

    if(myrank.eq.0 .and. iprint.ge.1) print *, 'Reading scalars from time step', ntime

    ! Time step
    call read_scalar(scalar_group_id, "time"        , time    , ntime, error)

    ! Magnetic Geometry
    call read_scalar(scalar_group_id, "xnull"       , xnull   , ntime, error)
    call read_scalar(scalar_group_id, "znull"       , znull   , ntime, error)
    call read_scalar(scalar_group_id, "xnull2"      , xnull2  , ntime, error)
    call read_scalar(scalar_group_id, "znull2"      , znull2  , ntime, error)
    call read_scalar(scalar_group_id, "xmag"        , xmag    , ntime, error)
    call read_scalar(scalar_group_id, "zmag"        , zmag    , ntime, error)

    ! Pellet stuff
    call read_scalar(scalar_group_id, "pellet_x",       pellet_x,      ntime, error)
    call read_scalar(scalar_group_id, "pellet_phi",     pellet_phi,    ntime, error)
    call read_scalar(scalar_group_id, "pellet_z",       pellet_z,      ntime, error)
    call read_scalar(scalar_group_id, "pellet_velx",    pellet_velx,   ntime, error)
    call read_scalar(scalar_group_id, "pellet_velphi",  pellet_velphi, ntime, error)
    call read_scalar(scalar_group_id, "pellet_velz",    pellet_velz,   ntime, error)
    call read_scalar(scalar_group_id, "pellet_var",     pellet_var,    ntime, error)
    call read_scalar(scalar_group_id, "r_p",            r_p,           ntime, error)
    call read_scalar(scalar_group_id, "r_p2",           r_p2,          ntime, error)
    call read_scalar(scalar_group_id, "pellet_rate",    pellet_rate,   ntime, error)
    call read_scalar(scalar_group_id, "pellet_rate1",   pellet_rate1,  ntime, error)
    call read_scalar(scalar_group_id, "pellet_rate2",   pellet_rate2,  ntime, error)
    call read_scalar(scalar_group_id, "pellet_ablrate", pellet_ablrate,ntime, error)
    
    ! Controllers
    call read_scalar(scalar_group_id, "loop_voltage",                  vloop,     ntime, error)
    call read_scalar(scalar_group_id, "i_control%err_i",     i_control%err_i,     ntime, error)
    call read_scalar(scalar_group_id, "i_control%err_p_old", i_control%err_p_old, ntime, error)
    call read_scalar(scalar_group_id, "n_control%err_i",     n_control%err_i,     ntime, error)
    call read_scalar(scalar_group_id, "n_control%err_p_old", n_control%err_p_old, ntime, error)

    call h5gclose_f(scalar_group_id, error)

    ! Read Fields
    call read_fields(time_id, 0, error)

    call h5gclose_f(time_id, error)

    ! Read Equilibrium Fields
    if(eqsubtract.eq.1) then
       call h5gopen_f(root_id, "equilibrium", eq_time_id, error)
       call read_fields(eq_time_id, 1, error)
       call h5gclose_f(eq_time_id, error)
    end if

    call h5gclose_f(root_id, error)

  end subroutine rdrestart_hdf5


  subroutine read_fields(time_group_id, equilibrium, error)
    use basic
    use hdf5
    use mesh_mod
    use field
    use arrays
    use hdf5_output

    implicit none

    integer(HID_T), intent(in) :: time_group_id
    integer, intent(out) :: error
    integer, intent(in) :: equilibrium

    integer(HID_T) :: group_id
    integer :: nelms, ilin

    ilin = 1 - equilibrium
    error = 0
    call hdf5_get_local_elms(nelms, error)

    call h5gopen_f(time_group_id, "fields", group_id, error)


    call h5r_read_field(group_id, "I",    bz_field(ilin), nelms, error)

    call h5r_read_field(group_id, "phi",   u_field(ilin), nelms, error)
    call h5r_read_field(group_id, "V",    vz_field(ilin), nelms, error)
    call h5r_read_field(group_id, "chi", chi_field(ilin), nelms, error)

    if(icsubtract.eq.1) then
       call h5r_read_field(group_id, "psi_coil", psi_coil_field, nelms, error)
    end if
    
    if(icsubtract.eq.1 .or. &
         (extsubtract.eq.1 .and. (ilin.eq.1 .or. eqsubtract.eq.0))) then
       call h5r_read_field(group_id, "psi_plasma", psi_field(ilin), nelms, error)
    else
       call h5r_read_field(group_id, "psi", psi_field(ilin), nelms, error)
    end if

    if(extsubtract.eq.1 .and. (ilin.eq.1 .or. eqsubtract.eq.0)) then
       call h5r_read_field(group_id, "I_plasma", bz_field(ilin), nelms, error)
       if(ifout.eq.1) then
          call h5r_read_field(group_id, "f_plasma", bf_field(ilin), nelms, error)
       end if
    else
       call h5r_read_field(group_id, "I", bz_field(ilin), nelms, error)
       if(ifout.eq.1) then
          call h5r_read_field(group_id, "f", bf_field(ilin), nelms, error)
       end if
    end if

    if(use_external_fields) then
       call h5r_read_field(group_id, "psi_ext", psi_ext, nelms, error)
       call h5r_read_field(group_id,   "I_ext",  bz_ext, nelms, error)
       call h5r_read_field(group_id,   "f_ext",  bf_ext, nelms, error)       
    end if

    if(jadv.eq.0) then
       call h5r_read_field(group_id, "potential", e_field(ilin), nelms, error)
    endif

    call h5r_read_field(group_id, "P",   p_field(ilin),   nelms, error)
    call h5r_read_field(group_id, "Pe",  pe_field(ilin),  nelms, error)
    call h5r_read_field(group_id, "den", den_field(ilin), nelms, error)
    call h5r_read_field(group_id, "te",  te_field(ilin),  nelms, error)
    call h5r_read_field(group_id, "ti",  ti_field(ilin),  nelms, error)
    
    call h5gclose_f(group_id, error)

  end subroutine read_fields

  subroutine h5r_read_field(group_id, name, f, nelms, error)
    use hdf5
    use field
    use hdf5_output

    implicit none

    integer(HID_T) :: group_id
    character(LEN=*) :: name
    type(field_type), intent(inout) :: f
    integer, intent(in) :: nelms
    integer, intent(out) :: error

    real, dimension(coeffs_per_element,nelms) :: dum
    vectype, dimension(coeffs_per_element,nelms) :: zdum
    integer :: i

    error = 0

    call read_field(group_id, name, dum, coeffs_per_element, &
         nelms, error)
    zdum = dum
#ifdef USECOMPLEX
    call read_field(group_id,name//"_i", dum, coeffs_per_element, &
         nelms,error)
    zdum = zdum + (0.,1.)*dum
#endif
    f = 0.
    do i=1, nelms
       call setavector(i, f, zdum(:,i))
    end do   
    
  end subroutine h5r_read_field

end module restart_hdf5
