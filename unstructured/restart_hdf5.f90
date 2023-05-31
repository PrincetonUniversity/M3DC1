module restart_hdf5
  implicit none

  integer, private :: icomplex_in, eqsubtract_in, extsubtract_in, ifin, nplanes_in
  integer, private :: ikprad_in, kprad_z_in, ipellet_in
  real, private, allocatable :: phi_in(:)
  real, private :: toroidal_period_in

contains
  
  subroutine rdrestart_hdf5()
    use basic
    use hdf5_output
    use hdf5
    use pellet
    use arrays
    use kprad_m3dc1
    use init_common
    use rmp

    implicit none

    include 'mpif.h'

    integer :: error
    integer(HID_T) :: root_id, scalar_group_id, time_id, eq_time_id, pel_group_id, mesh_id
    character(LEN=19) :: time_group_name

    integer :: times_output_in, i3d_in, istartnew, i
    real :: xnullt,znullt,xnull2t,znull2t

    call h5gopen_f(file_id, "/", root_id, error)

    call read_int_attr(root_id, "version", version_in, error)
    
    if (myrank.eq.0) print *, 'Reading HDF5 file: version=', version_in
  
    if(version_in.lt.16) then
       if(myrank.eq.0) print *, 'Error: HDF5 file is from too old a version to use iread_hdf5=1.'
       call h5gclose_f(root_id, error)
       call safestop(1)
    end if


    ! Read Time Slice
    ! ~~~~~~~~~~~~~~~
    if(irestart_slice.ge.0 .and. irestart_slice.lt.times_output) then
       ! unlink output files beyond irestart_slice
       do i = irestart_slice+1, times_output
          write(time_group_name, '("time_",I3.3)') i
          call h5gunlink_f(file_id, time_group_name, error)
       end do
       times_output = irestart_slice
    end if
    if(myrank.eq.0) print *, 'Reading data from time slice', times_output
    call read_int_attr(root_id, "ntime", times_output_in, error)
    if(times_output_in .le. times_output) then
       print *, 'Error: Requested restart time slice exceeds number of time slices in file.'
       call safestop(2)
    end if

    write(time_group_name, '("time_",I3.3)') times_output
    call h5gopen_f(root_id, time_group_name, time_id, error)
    call read_int_attr(time_id, "ntimestep", ntime, error)


    ! Read Attributes
    ! ~~~~~~~~~~~~~~~

    ! For version >= 17, the version number can change with timeslice
    if(version_in.ge.17) then
       call read_int_attr(time_id, "version", version_in, error)
    end if
    
    call read_int_attr(root_id, "eqsubtract", eqsubtract_in, error)
    call read_int_attr(root_id, "extsubtract", extsubtract_in, error)
    call read_int_attr(root_id, "icomplex", icomplex_in, error)
    call read_int_attr(root_id, "3d", i3d_in, error)

    call h5gopen_f(time_id, "mesh", mesh_id, error)
    call read_int_attr(mesh_id, "nplanes", nplanes_in, error)
    call read_real_attr(mesh_id, "period", toroidal_period_in, error)

    if(version_in.ge.34 .and. i3d_in.eq.1) then
       allocate(phi_in(nplanes_in))
       call read_vec_attr(mesh_id, "phi", phi_in, nplanes_in, error)
       if(myrank.eq.0) then
          print *, 'Read phi_in = ', phi_in
       end if
    end if

    call h5gclose_f(mesh_id, error)

    if(i3d_in.eq.1 .or. icomplex_in.eq.1) then
       ifin = 1
       if(i3d.eq.0) then
          if(myrank.eq.0) then
             print *, 'Error: cannot start an axisymmetric calculation'
             print *, '       from non-axisymmetric data'
          end if
          call h5gclose_f(root_id, error)
          call safestop(1)
       end if
       if(i3d_in.eq.1 .and. icomplex.eq.1) then
          if(myrank.eq.0) then
             print *, 'Error: cannot start a complex calculation'
             print *, '       from a non-axisymmetric equilibrium'
          end if
          call h5gclose_f(root_id, error)
          call safestop(1)
       end if

       ! check if fp is present
       if(version_in.ge.38 .or. (version_in.ge.35 .and. numvar.gt.1)) then
          irestart_fp = 1
          if(myrank.eq.0 .and. iprint.ge.2) print *, " fp is present at restart"
       else
          irestart_fp = 0
          if(myrank.eq.0 .and. iprint.ge.2) print *, " fp is absent at restart"
       endif
    else
       ifin = 0
    end if

    if(version_in.ge.18) then
       call read_real_attr(root_id, "xlim" , xlim,  error)
       call read_real_attr(root_id, "zlim" , zlim,  error)
       call read_real_attr(root_id, "xlim2", xlim2, error)
       call read_real_attr(root_id, "zlim2", zlim2, error)
    end if

    if(version_in.ge.19) then
       call read_int_attr(root_id, "ikprad", ikprad_in,  error)
       call read_int_attr(root_id, "kprad_z", kprad_z_in,  error)
    else
       ikprad_in = 0
    end if

    ! Read Scalars
    call h5gopen_f(root_id, "scalars", scalar_group_id, error)

    if(myrank.eq.0 .and. iprint.ge.1) print *, 'Reading scalars from time step', ntime

    ! Time step
    call read_scalar(scalar_group_id, "time"        , time    , ntime, error)

    ! Magnetic Geometry
    call read_scalar(scalar_group_id, "xnull"       , xnullt   , ntime, error)
    call read_scalar(scalar_group_id, "znull"       , znullt   , ntime, error)
    call read_scalar(scalar_group_id, "xnull2"      , xnull2t  , ntime, error)
    call read_scalar(scalar_group_id, "znull2"      , znull2t  , ntime, error)
    call read_scalar(scalar_group_id, "xmag"        , xmag     , ntime, error)
    call read_scalar(scalar_group_id, "zmag"        , zmag     , ntime, error)
    call read_scalar(scalar_group_id, "psi_lcfs"    , psibound , ntime, error)
    call read_scalar(scalar_group_id, "psimin"      , psimin   , ntime, error)

    if(mod_null_rs .eq.0) then
       xnull = xnullt
       znull = znullt
    endif
    if(mod_null_rs2 .eq.0) then
       xnull2 = xnull2t
       znull2 = znull2t
    endif
   
    if(ipellet.ne.0) then
       ! Pellet stuff
       if(version_in.le.30) then
          ! pellets are scalars
          npellets = 1
          if(version_in.le.25) then
             call read_scalar(scalar_group_id, "pellet_x",       pellet_r(1),      ntime, error)
          else
             call read_scalar(scalar_group_id, "pellet_r",       pellet_r(1),      ntime, error)
          end if
          call read_scalar(scalar_group_id, "pellet_phi",     pellet_phi(1),    ntime, error)
          call read_scalar(scalar_group_id, "pellet_z",       pellet_z(1),      ntime, error)
          if(version_in.le.25) then
             call read_scalar(scalar_group_id, "pellet_velx",    pellet_velr(1),   ntime, error)
          else
             call read_scalar(scalar_group_id, "pellet_velr",    pellet_velr(1),   ntime, error)
          end if
          call read_scalar(scalar_group_id, "pellet_velphi",  pellet_velphi(1), ntime, error)
          call read_scalar(scalar_group_id, "pellet_velz",    pellet_velz(1),   ntime, error)
          if(version_in.ge.26) then
             call read_scalar(scalar_group_id, "pellet_vx",    pellet_vx(1),   ntime, error)
             call read_scalar(scalar_group_id, "pellet_vy",    pellet_vy(1),   ntime, error)
          end if
          call read_scalar(scalar_group_id, "pellet_var",     pellet_var(1),    ntime, error)
          call read_scalar(scalar_group_id, "r_p",            r_p(1),           ntime, error)
          call read_scalar(scalar_group_id, "pellet_rate",    pellet_rate(1),   ntime, error)
       else
          ! pellets are arrays
          call read_int_attr(root_id, "ipellet", ipellet_in,  error)
          if(ipellet_in.ne.0) then
             call read_int_attr(root_id, "npellets", npellets,  error)
             call h5gopen_f(root_id, "pellet", pel_group_id, error)
             
             if(irestart_pellet.eq.0) then     
                 ! Read all pellet parameters from restart .h5 file
                 call read_1dextendarr(pel_group_id, "pellet_r",       pellet_r,       npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "pellet_phi",     pellet_phi,     npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "pellet_z",       pellet_z,       npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "pellet_rate",    pellet_rate,    npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "pellet_rate_D2", pellet_rate_D2, npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "pellet_velr",    pellet_velr,    npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "pellet_velphi",  pellet_velphi,  npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "pellet_velz",    pellet_velz,    npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "pellet_vx",      pellet_vx,      npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "pellet_vy",      pellet_vy,      npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "r_p",            r_p,            npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "pellet_mix",     pellet_mix,     npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "pellet_var_tor", pellet_var_tor, npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "cloud_pel",      cloud_pel,      npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "pellet_var",     pellet_var,     npellets, ntime, error)
                 if(version_in.ge.33) then
                    call read_1dextendarr(pel_group_id, "cauchy_fraction",cauchy_fraction, npellets, ntime, error)
                 end if
                 if(version_in.ge.37) then
                    call read_1dextendarr(pel_group_id, "pellet_te", temp_pel, npellets, ntime, error)
                    temp_pel = temp_pel*p0_norm/(n0_norm*1.6022e-12)
                    call read_1dextendarr(pel_group_id, "pellet_ne", nsource_pel, npellets, ntime, error)
                 end if

             else if(irestart_pellet.eq.1) then
                 ! Control ablated cloud size and others parameters from C1input
                 call read_1dextendarr(pel_group_id, "pellet_r",       pellet_r,       npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "pellet_phi",     pellet_phi,     npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "pellet_z",       pellet_z,       npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "pellet_velr",    pellet_velr,    npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "pellet_velphi",  pellet_velphi,  npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "pellet_velz",    pellet_velz,    npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "pellet_vx",      pellet_vx,      npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "pellet_vy",      pellet_vy,      npellets, ntime, error)
                 call read_1dextendarr(pel_group_id, "r_p",            r_p,            npellets, ntime, error)
                 if(version_in.ge.37) then
                    call read_1dextendarr(pel_group_id, "pellet_te", temp_pel, npellets, ntime, error)
                    call read_1dextendarr(pel_group_id, "pellet_ne", nsource_pel, npellets, ntime, error)
                 end if
             else
                 if(myrank.eq.0) then
                    print *, 'Error: irestart_pellet not available for this value: ',irestart_pellet
                    call safestop(30)
                end if         
             end if
             call h5gclose_f(pel_group_id, error)
          end if
       end if
    endif

    ! Only read vloop if Ip is under current control
    if(control_type.ne.-1) then
       call read_scalar(scalar_group_id, "loop_voltage", vloop, ntime, error)
       if(version_in.ge.17 .and. version_in.lt.39 .and. itor.eq.0) then
          if(myrank.eq.0) then
             print *, "Multiplying VLOOP by RZERO to account for new definition."
          end if
          vloop = vloop*rzero
       end if
    end if
    call read_scalar(scalar_group_id, "i_control%err_i",     i_control%err_i,     ntime, error)
    call read_scalar(scalar_group_id, "i_control%err_p_old", i_control%err_p_old, ntime, error)
    call read_scalar(scalar_group_id, "n_control%err_i",     n_control%err_i,     ntime, error)
    call read_scalar(scalar_group_id, "n_control%err_p_old", n_control%err_p_old, ntime, error)

    call h5gclose_f(scalar_group_id, error)

    ! Read Fields
    if(myrank.eq.0 .and. iprint.ge.1) print *, 'Reading fields from time step', ntime

    call read_fields(time_id, 0, error)

    call h5gclose_f(time_id, error)

    ! Read Equilibrium Fields
    if(eqsubtract_in.eq.1) then
       if(myrank.eq.0 .and. iprint.ge.1) print *, 'Reading equilibrium fields'

       call h5gopen_f(root_id, "equilibrium", eq_time_id, error)
       call read_fields(eq_time_id, 1, error)
       call h5gclose_f(eq_time_id, error)
    end if

    call h5gclose_f(root_id, error)


    ! If type of calculation has changed (i.e. real to complex)
    ! then overwrite output
    istartnew = 0

    ! If eqsubtract = 1 but eqsubtract_in = 0, then
    ! old fields become new equilibrium fields
    if(eqsubtract_in.eq.0 .and. eqsubtract.eq.1) then
       if(myrank.eq.0) then
          print *, 'Restarting while changing eqsubtract from 0 to 1'
          print *, 'Previous data will be overwritten.'
       end if
       field0_vec = field_vec
       field_vec = 0.
       bfp_field(0) = bfp_field(1)
       bfp_field(1) = 0. 
       istartnew = 1
       time = 0
    end if
    if(icomplex.eq.1 .and. icomplex_in.eq.0) then
       if(myrank.eq.0) then
          print *, 'Starting complex calculation from 2D real calculation.'
          print *, 'Previous data will be overwritten.'
       end if
       istartnew = 1
    end if
    if(nplanes_in.eq.1 .and. nplanes.gt.1) then
       if(myrank.eq.0) then
          print *, 'Starting 3D calculation from 2D calculation.'
          print *, 'Previous data will be overwritten.'
       end if
       istartnew = 1
    end if
    if(toroidal_period_in.ne.toroidal_period) then
       if(myrank.eq.0) then
          print *, 'Starting full torus calculation from one field period.'
          print *, 'Previous data will be overwritten.'
       end if
       istartnew = 1
       time = 0
    end if

    if (istartnew.eq.1) then
       ntime = 0
       irestart = 0
       call hdf5_finalize(error)
       call hdf5_initialize(.false., error)

       if(eqsubtract.eq.0) then
         ! move data to field0 so that equilibrium data is written correctly
         field0_vec = field_vec
         field_vec = 0.
         nre_field(0) = nre_field(1)
         nre_field(1) = 0.
       endif
       call init_perturbations
       if(eqsubtract.eq.0) then
         call add_field_to_field(nre_field(1),nre_field(0))
         nre_field(0) = 0.
      endif

      ! For RMP and error fields
      if(irmp.ge.1 .or. iread_ext_field.ge.1 .or. &
           tf_tilt.ne.0. .or. tf_shift.ne.0. .or. &
           any(pf_tilt.ne.0.) .or. any(pf_shift.ne.0.)) then
         ! External fields already loaded for itaylor = 41
         if(itaylor.eq.41) then
            if(myrank.eq.0 .and. iprint.ge.2) print *, &
                 "Skipping: RMP specification not currently implemented for ST."
         else
            call rmp_per
         end if
      end if

    end if

    if(allocated(phi_in)) deallocate(phi_in)
  end subroutine rdrestart_hdf5


  subroutine read_fields(time_group_id, equilibrium, error)
    use basic
    use hdf5
    use h5lt
    use mesh_mod
    use field
    use arrays
    use hdf5_output
    use kprad_m3dc1

    implicit none

    integer(HID_T), intent(in) :: time_group_id
    integer, intent(out) :: error
    integer, intent(in) :: equilibrium

    integer(HID_T) :: group_id
    integer :: nelms, ilin, i
    character(len=64) :: field_name

    ilin = 1 - equilibrium
    error = 0
    call hdf5_get_local_elms(nelms, error)

!    if(nplanes.ne.nplanes_in) then
!       if(myrank.eq.0) print *, '3D gloabl_nelms: ', global_elms
!       global_elms = global_elms * nplanes_in / nplanes
!       if(myrank.eq.0) print *, '2D global_nelms: ', global_elms
!       offset = modulo(offset, global_elms)
!    end if

    call h5gopen_f(time_group_id, "fields", group_id, error)

    call h5r_read_field(group_id, "I",    bz_field(ilin), nelms, error)

    call h5r_read_field(group_id, "phi",   u_field(ilin), nelms, error)
!    call m3dc1_field_write(u_field(ilin)%vec%id, "phi-r", 0)

    call h5r_read_field(group_id, "V",    vz_field(ilin), nelms, error)
    call h5r_read_field(group_id, "chi", chi_field(ilin), nelms, error)
    if(irunaway.gt.0) then
      call h5r_read_field(group_id, "nre", nre_field(ilin), nelms, error)
    endif

    if(icsubtract.eq.1) then
       call h5r_read_field(group_id, "psi_coil", psi_coil_field, nelms, error, .true.)
    end if
    
    if(icsubtract.eq.1 .or. &
         (extsubtract_in.eq.1 .and. (ilin.eq.1 .or. eqsubtract_in.eq.0))) then
       call h5r_read_field(group_id, "psi_plasma", psi_field(ilin), nelms, error)
    else
       call h5r_read_field(group_id, "psi", psi_field(ilin), nelms, error)
    end if

    if(extsubtract_in.eq.1 .and. (ilin.eq.1 .or. eqsubtract_in.eq.0)) then
       call h5r_read_field(group_id, "I_plasma", bz_field(ilin), nelms, error)
       if(ifin.eq.1) then
          call h5r_read_field(group_id, "f_plasma", bf_field(ilin), nelms, error)
       end if
       if(irestart_fp.eq.1) then
          call h5r_read_field(group_id, "fp_plasma", bfp_field(ilin), nelms, error)
       end if
    else
       call h5r_read_field(group_id, "I", bz_field(ilin), nelms, error)
       if(ifin.eq.1) then
          call h5r_read_field(group_id, "f", bf_field(ilin), nelms, error)
       end if
       if(irestart_fp.eq.1) then
          call h5r_read_field(group_id, "fp", bfp_field(ilin), nelms, error)
       end if
    end if

    if (extsubtract_in.eq.1) then
       call h5r_read_field(group_id, "psi_ext", psi_ext, nelms, error)
       call h5r_read_field(group_id,   "I_ext",  bz_ext, nelms, error)
       call h5r_read_field(group_id,   "f_ext",  bf_ext, nelms, error)       
       if(irestart_fp.eq.1) then
          call h5r_read_field(group_id, "fp_ext", bfp_ext, nelms, error)
       end if
    end if

    if(jadv.eq.0) then
       call h5r_read_field(group_id, "potential", e_field(ilin), nelms, error)
    endif

    call h5r_read_field(group_id, "P",   p_field(ilin),   nelms, error)
    call h5r_read_field(group_id, "Pe",  pe_field(ilin),  nelms, error)
    call h5r_read_field(group_id, "den", den_field(ilin), nelms, error)
    call h5r_read_field(group_id, "ne",  ne_field(ilin),  nelms, error)
    call h5r_read_field(group_id, "te",  te_field(ilin),  nelms, error)
    call h5r_read_field(group_id, "ti",  ti_field(ilin),  nelms, error)

    if(ikprad.eq.1 .and. ikprad_in.eq.1) then
       do i=0, kprad_z
          write(field_name, '(A,I2.2)') "kprad_n_", i
          call h5r_read_field(group_id,trim(field_name),kprad_n(i),nelms,error)
       end do
       call h5r_read_field(group_id,"kprad_sigma_e",kprad_sigma_e,nelms,error)
       call h5r_read_field(group_id,"kprad_sigma_i",kprad_sigma_i,nelms,error)
       call h5r_read_field(group_id,"kprad_rad",    kprad_rad,    nelms,error)
       if(version.ge.22) then 
          call h5r_read_field(group_id,"kprad_brem",kprad_brem,   nelms,error)
       end if
       if(version.eq.23) then
          call h5r_read_field(group_id,"kprad_ion",kprad_ion,   nelms,error)
          call h5r_read_field(group_id,"kprad_rec",kprad_reck,   nelms,error)
       else if(version.ge.24) then
          call h5r_read_field(group_id,"kprad_ion",kprad_ion,   nelms,error)
          call h5r_read_field(group_id,"kprad_reck",kprad_reck,   nelms,error)
          call h5r_read_field(group_id,"kprad_recp",kprad_recp,   nelms,error)
       end if
    end if
    
    call h5gclose_f(group_id, error)

  end subroutine read_fields

  subroutine h5r_read_field(group_id, name, f, nelms, error, isreal)
    use basic
    use hdf5
    use field
    use hdf5_output

    implicit none

    integer(HID_T) :: group_id
    character(LEN=*) :: name
    type(field_type), intent(inout) :: f
    integer, intent(in) :: nelms
    integer, intent(out) :: error
    logical, intent(in), optional :: isreal

    real, dimension(coeffs_per_element,nelms) :: dum
    vectype, dimension(coeffs_per_element,nelms) :: zdum
    vectype, dimension(coeffs_per_element) :: kdum
    integer :: i, coefs
    logical :: ir
    integer :: elms_per_plane, new_plane, old_plane, plane_fac, k
    integer :: offset_in, global_elms_in
    real :: dphi, shift, phi_new
    logical :: transform
!    vectype, dimension(dofs_per_element, dofs_per_element) :: trans_mat

    if(myrank.eq.0 .and. iprint.ge.2) print *, 'Reading ', name

    if(present(isreal)) then
       ir = isreal
    else
       ir = .false.
    end if

    error = 0

    if(nplanes_in.eq.1) then
       coefs = coeffs_per_tri
       global_elms_in = global_elms * nplanes_in / nplanes
       offset_in = modulo(offset, global_elms_in)
       dum = 0.
       transform = .false.
    else if(nplanes_in .ne. nplanes) then 
       if((mod(nplanes,nplanes_in).ne.0 .or. nplanes .lt. nplanes_in) &
            .and. version_in.lt.34) then 
          if(myrank.eq.0) then
             print *, 'Error: new nplanes must be integer multiple of existing nplanes.'
             call safestop(42)
          end if
       else
          if(myrank.eq.0) print *, "Restarting ", nplanes_in, " planes case with ", nplanes, " planes"
       end if
       coefs = coeffs_per_element
       plane_fac = nplanes / nplanes_in
       elms_per_plane = global_elms / nplanes
       new_plane = offset / elms_per_plane
       global_elms_in = global_elms / plane_fac

       if(version_in.lt.34) then
          old_plane = new_plane / plane_fac
          dphi = toroidal_period / nplanes_in
          k = new_plane - old_plane * plane_fac
          shift = k*dphi / plane_fac
       else
          call m3dc1_plane_getphi(new_plane, phi_new)
          if(toroidal_period_in.ne.toroidal_period) then
             phi_new = modulo(phi_new, toroidal_period_in)
          end if
          call find_plane_and_shift(nplanes_in, phi_in, phi_new, old_plane, shift)
          old_plane = old_plane - 1    ! switch to 0-based indexing
       end if
       offset_in = offset - elms_per_plane*(new_plane - old_plane)
       transform = .true.
    else
       coefs = coeffs_per_element
       offset_in = offset
       global_elms_in = global_elms
       transform = .false.
    end if

    call read_field(group_id, name, dum(1:coefs,:), coefs, nelms, &
         offset_in, global_elms_in, error)
    zdum = dum
    if(icomplex_in.eq.1 .and. .not.ir) then
       call read_field(group_id,name//"_i", dum(1:coefs,:), coefs, nelms, &
            offset_in, global_elms_in, error)
#ifdef USECOMPLEX
       zdum = zdum + (0.,1.)*dum
#endif
    end if
    f = 0.

    do i=1, nelms
       if(transform) then
          call transform_coeffs_nplanes(zdum(:,i),shift,kdum)
          call setavector(i, f, kdum)
       else
          call setavector(i, f, zdum(:,i))
       end if
    end do
    
  end subroutine h5r_read_field

  subroutine find_plane_and_shift(nplanes_old, phi_old, phi_new, iplane, shift)
    use mesh_mod
    implicit none

    integer, intent(in) :: nplanes_old
    real, intent(in), dimension(nplanes_old) :: phi_old
    real, intent(in) :: phi_new
    integer, intent(out) :: iplane
    real, intent(out) :: shift
    real :: tol

    integer :: i

    tol = toroidal_period*1e-5

    iplane = -1
    do i=1, nplanes_old-1
       if(phi_new.ge.phi_old(i)-tol .and. phi_new.lt.phi_old(i+1)-tol) then
          iplane = i
          exit
       end if
    end do

    if(iplane.lt.0) iplane=nplanes_old
    shift = phi_new - phi_old(iplane)
  end subroutine find_plane_and_shift

end module restart_hdf5
