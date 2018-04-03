module m3dc1_output

  integer, parameter, private :: ke_file = 11
  character*(*), parameter, private :: ke_filename = 'C1ke'

  integer :: iwrite_transport_coeffs
  integer :: iwrite_aux_vars

contains

  subroutine initialize_output
    use basic
    use hdf5_output
    implicit none

    integer :: ier
    
    call hdf5_initialize(irestart.ne.0,ier)
    if(ier.lt.0) then 
       print *, "Error initializing HDF5"
       call safestop(5)
    end if

    if(myrank.eq.0) then
       open(unit=ke_file, file=ke_filename, status='unknown')
    endif
  end subroutine initialize_output

  subroutine finalize_output
    use basic
    use hdf5_output
    implicit none

    integer :: ier

    if(myrank.eq.0) then
       close(ke_file)
    end if

    call hdf5_finalize(ier)
    if(ier.ne.0) print *, 'Error finalizing HDF5:',ier
  end subroutine finalize_output

  ! ======================================================================
  ! output
  ! ~~~~~~
  !
  ! writes output and restart files
  ! ======================================================================
  subroutine output
    use basic
    use hdf5_output
    use diagnostics
    use auxiliary_fields

    implicit none

    include 'mpif.h'

    integer :: ier,i
    real :: tstart, tend, diff

    if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
    call hdf5_write_scalars(ier)

#ifdef USE3D
    if(ike_harmonics .gt. 0) call hdf5_write_keharmonics(ier)
    if(ibh_harmonics .gt. 0) call hdf5_write_bharmonics(ier)
    call hdf5_write_kspits(ier)
#endif

    if(myrank.eq.0 .and. itimer.eq.1) then
      call second(tend)
      diff = tend - tstart
      t_output_hdf5 = t_output_hdf5 + diff
      if(iprint.ge.1) write(*,1003) ntime,diff,t_output_hdf5
1003  format("OUTPUT: hdf5_write_scalars   ", I5, 1p2e16.8)
    endif

    ! only write field data evey ntimepr timesteps
    if(mod(ntime-ntime0,ntimepr).eq.0) then
       if(iwrite_aux_vars.eq.1) then
          if(myrank.eq.0 .and. iprint.ge.2) print *, "  calculating aux fields"
          call calculate_auxiliary_fields(eqsubtract)
       end if

       if(myrank.eq.0 .and. iprint.ge.2) print *, "  writing timeslice"
       if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

       call hdf5_write_time_slice(0,ier)
       if(myrank.eq.0 .and. itimer.eq.1) then
          call second(tend)
          diff = tend - tstart
          t_output_hdf5 = t_output_hdf5 + diff
          if(iprint.ge.1) write(*,1004) ntime,diff,t_output_hdf5
1004      format("OUTPUT: hdf5_write_time_slice", I5, 1p2e16.8)
       endif
    endif

    if(itimer.eq.1) then
       if(myrank.eq.0) call second(tstart)
       if(myrank.eq.0 .and. iprint.ge.2) print *, "  writing timings"
       call hdf5_write_timings(ier)
       if(myrank.eq.0) then
         call second(tend)
         diff = tend - tstart
         t_output_hdf5 = t_output_hdf5 + diff
         if(iprint.ge.1) write(*,1005) ntime,diff,t_output_hdf5,ier
  1005 format("OUTPUT: hdf5_write_timings   ", I5, 1p2e16.8,i5)
       endif
       call reset_timings
    end if
    
    ! flush hdf5 data to disk
    if(myrank.eq.0 .and. iprint.ge.2) print *, "  flushing data to file"
    if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
    call hdf5_flush(ier)
    if(myrank.eq.0 .and. itimer.eq.1) then
       call second(tend)
       diff = tend - tstart
       t_output_hdf5 = t_output_hdf5 + diff
       if(iprint.ge.1) write(*,1002) ntime,diff,t_output_hdf5,ier
  1002 format("OUTPUT: hdf5_flush           ", I5, 1p2e16.8,i5)
    end if    

    ! only write restart file evey ntimers timesteps
    if(iwrite_restart.eq.1 .and. mod(ntime-ntime0,ntimers).eq.0 .and. ntime.ne.ntime0) then
       if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
       if(myrank.eq.0 .and. iprint.ge.1) print *, "  writing restart files"
       if(iglobalout.eq.1) then
          call wrrestartglobal
       else
          if(iwrite_adios.eq.1) then
             call wrrestart_adios
          else
!...........sequential restart writing
             do i=0,maxrank-1
                if(myrank.eq.i) call wrrestart
                call MPI_Barrier(MPI_COMM_WORLD,ier)
             enddo
          end if
       endif
       if(myrank.eq.0 .and. itimer.eq.1) then
          call second(tend)
          diff = tend - tstart
          if(iprint.ge.1) write(*,1006) ntime,diff,t_output_hdf5
1006      format("OUTPUT: wrrestart            ", I5, 1p2e16.8)
       endif
    endif
 

    ! Write C1ke data
    if(myrank.eq.0) then
       if((ekin+ekino)*dtold.eq.0. .or. ekin.eq.0.) then
          gamma_gr = 0.
       else
          gamma_gr = (ekin - ekino)/((ekin+ekino)*dtold)
       endif
       write(ke_file, '(I8, 1p3e12.4,2x,1p3e12.4,2x,1p3e12.4,2x,1pe13.5)') &
            ntime, time, ekin, gamma_gr, &
            ekinp,ekint,ekin3, emagp, emagt, emag3, etot
    endif
  end subroutine output


! hdf5_write_parameters
! =====================
subroutine hdf5_write_parameters(error)
  use hdf5
  use hdf5_output
  use basic
  use pellet
  use bootstrap
  use diagnostics
  use resistive_wall
  use kprad_m3dc1

  implicit none

  integer, intent(out) :: error

  integer(HID_T) :: root_id

  if(myrank.eq.0 .and. iprint.ge.2) print *, 'Writing attributes'

  call h5gopen_f(file_id, "/", root_id, error)

#ifdef USE3D
  call write_int_attr (root_id, "3d"         , 1,        error)
#else
  call write_int_attr (root_id, "3d"         , 0,        error)
#endif

  call write_int_attr (root_id, "icomplex"   , icomplex,   error)

  call write_int_attr (root_id, "nplanes"    , nplanes,    error)

  call write_int_attr (root_id, "version"    , version,    error)
  call write_int_attr (root_id, "numvar"     , numvar,     error)
  call write_int_attr (root_id, "idens"      , idens,      error)
  call write_int_attr (root_id, "ipres"      , ipres,      error)
  call write_int_attr (root_id, "itemp"      , itemp,      error)
  call write_int_attr (root_id, "ipressplit" , ipressplit, error)
  call write_int_attr (root_id, "itime_independent", itime_independent, error)
  call write_int_attr (root_id, "itor"       , itor,       error)
  call write_int_attr (root_id, "gyro"       , gyro,       error)
  call write_int_attr (root_id, "linear"     , linear,     error)
  call write_int_attr (root_id, "kinetic"    , kinetic,    error)
  call write_int_attr (root_id, "eqsubtract" , eqsubtract, error)
  call write_int_attr (root_id, "extsubtract", extsubtract,error)
  call write_int_attr (root_id, "icsubract"  , icsubtract, error)
  call write_int_attr (root_id, "iper"       , iper,       error)
  call write_int_attr (root_id, "jper"       , jper,       error)
  call write_int_attr (root_id, "integrator" , integrator, error)
  call write_int_attr (root_id, "ipellet"    , ipellet,    error)
  call write_int_attr (root_id, "ipellet_abl", ipellet_abl,error)
  call write_int_attr (root_id, "ivform"     , ivform,     error)
  call write_int_attr (root_id, "ntor"       , ntor,       error)
  call write_int_attr (root_id, "nonrect"    , nonrect,    error)
  call write_int_attr (root_id, "ifixedb"    , ifixedb,    error)
  call write_int_attr (root_id, "imulti_region", imulti_region, error)
  call write_real_attr(root_id, "db"         , db,         error)
  call write_real_attr(root_id, "rzero"      , rzero,      error)
  call write_real_attr(root_id, "gam"        , gam,        error)
  call write_real_attr(root_id, "thimp"      , thimp,      error)
  call write_real_attr(root_id, "bzero"      , bzero,      error)
  call write_real_attr(root_id, "gravr"      , gravr,      error)
  call write_real_attr(root_id, "gravz"      , gravz,      error)
  call write_real_attr(root_id, "amu"        , amu,        error)
  call write_real_attr(root_id, "amuc"       , amuc,       error)
  call write_real_attr(root_id, "amupar"     , amupar,     error)
  call write_real_attr(root_id, "etar"       , etar,       error)
  call write_real_attr(root_id, "eta0"       , eta0,       error)
  call write_real_attr(root_id, "kappar"     , kappar,     error)
  call write_real_attr(root_id, "kappa0"     , kappa0,     error)
  call write_real_attr(root_id, "kappat"     , kappat,     error)
  call write_real_attr(root_id, "denm"       , denm,       error)
  call write_real_attr(root_id, "ln"         , ln,         error)
  call write_real_attr(root_id, "hyper"      , hyper,      error)
  call write_real_attr(root_id, "hyperi"     , hyperi,     error)
  call write_real_attr(root_id, "hyperv"     , hyperv,     error)
  call write_real_attr(root_id, "hyperc"     , hyperc,     error)
  call write_real_attr(root_id, "hyperp"     , hyperp,     error)
  call write_real_attr(root_id, "b0_norm"    , b0_norm,    error)
  call write_real_attr(root_id, "n0_norm"    , n0_norm,    error)
  call write_real_attr(root_id, "l0_norm"    , l0_norm,    error)  
  call write_real_attr(root_id, "eta_wall"   , eta_wall,   error)
  call write_real_attr(root_id, "zeff"       , zeff,       error)
  call write_real_attr(root_id, "ion_mass"   , ion_mass,   error)
  call write_real_attr(root_id, "frequency"  , frequency,  error)
  call write_int_attr (root_id, "ibootstrap_model", ibootstrap_model, error)
  call write_real_attr(root_id, "bootstrap_alpha", bootstrap_alpha, error)
  call write_real_attr(root_id, "eta_te_offset", eta_te_offset, error)
  call write_int_attr (root_id, "imag_probes", imag_probes, error)
  call write_int_attr (root_id, "iflux_loops", iflux_loops, error)
  call write_real_attr(root_id, "tflux0"     , tflux0,      error)
  call write_real_attr(root_id, "xlim"       , xlim,        error)
  call write_real_attr(root_id, "zlim"       , zlim,        error)
  call write_real_attr(root_id, "xlim2"      , xlim2,       error)
  call write_real_attr(root_id, "zlim2"      , zlim2,       error)
  call write_int_attr (root_id, "ikprad"     , ikprad,      error)
  call write_int_attr (root_id, "kprad_z"    , kprad_z,     error)

  call h5gclose_f(root_id, error)

end subroutine hdf5_write_parameters

subroutine hdf5_reconcile_version(ver, error)
  use basic
  use hdf5
  use hdf5_output
  use kprad_m3dc1

  implicit none

  integer, intent(in) :: ver
  integer, intent(out) :: error
  integer(HID_T) :: root_id, scalar_group_id

  if(ver.ge.version) return

  if(myrank.eq.0 .and. iprint.ge.1) then
     print *, 'Reconciling version data from ', version_in, ' to ', version
  end if

  call h5gopen_f(file_id, "/", root_id, error)

  if(ver.lt.17) then
     call h5gopen_f(root_id, "scalars", scalar_group_id, error)
     call write_int_attr(scalar_group_id, "ntimestep", ntime, error)
     call h5gclose_f(scalar_group_id, error)

     ! This is necessary because if version < 17 in the restart file,
     ! we can't find the version of the latest restart.  This will cause
     ! subsequent restarts will fail due to the above
     ! write_int_attr statement.  However, in general we will leave the 
     ! "version" variable in root_id alone when restarting with a new version
     ! (instead, the version change is reflected in the time slice "version").
     call update_int_attr(root_id, "version", 17, error)
  end if

  if(ver.lt.18) then
     call write_real_attr(root_id, "xlim", xlim, error)
     call write_real_attr(root_id, "zlim", zlim, error)
     call write_real_attr(root_id, "xlim2", xlim2, error)
     call write_real_attr(root_id, "zlim2", zlim2, error)
  end if

  if(ver.lt.19) then
     call write_int_attr(root_id, "ikprad", ikprad, error)
     call write_int_attr(root_id, "kprad_z", kprad_z, error)
  end if

  call h5gclose_f(root_id, error)

end subroutine hdf5_reconcile_version

! hdf5_write_scalars
! ==================
subroutine hdf5_write_scalars(error)
  use basic
  use diagnostics
  use hdf5_output
  use pellet

  implicit none

  integer, intent(out) :: error

  integer(HID_T) :: root_id, scalar_group_id, fl_group_id, mp_group_id

  real :: temp

  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, " Writing scalars"

  call h5gopen_f(file_id, "/", root_id, error)

  if(ntime.eq.0 .and. irestart.eq.0) then
     call h5gcreate_f(root_id, "scalars", scalar_group_id, error)
     call write_int_attr(scalar_group_id, "ntimestep", ntime, error)
     if(imag_probes.ne.0) call h5gcreate_f(root_id, "mag_probes", mp_group_id, error)
     if(iflux_loops.ne.0) call h5gcreate_f(root_id, "flux_loops", fl_group_id, error)
  else
     call h5gopen_f(root_id, "scalars", scalar_group_id, error)
     call update_int_attr(scalar_group_id, "ntimestep", ntime, error)
     if(imag_probes.ne.0) call h5gopen_f(root_id, "mag_probes", mp_group_id, error)
     if(iflux_loops.ne.0) call h5gopen_f(root_id, "flux_loops", fl_group_id, error)
  endif

  ! State Variables (needed for restart)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! Time step
  call output_scalar(scalar_group_id, "time" , time, ntime, error)
  call output_scalar(scalar_group_id, "dt" ,     dt, ntime, error)

  ! Magnetic geometry
  call output_scalar(scalar_group_id, "xnull"   ,xnull   ,ntime,error)
  call output_scalar(scalar_group_id, "znull"   ,znull   ,ntime,error)
  call output_scalar(scalar_group_id, "xnull2"  ,xnull2  ,ntime,error)
  call output_scalar(scalar_group_id, "znull2"  ,znull2  ,ntime,error)
  call output_scalar(scalar_group_id, "xmag"    ,xmag    ,ntime,error)
  call output_scalar(scalar_group_id, "zmag"    ,zmag    ,ntime,error)

  ! Pellet stuff
  call output_scalar(scalar_group_id, "pellet_x",   pellet_x,   ntime, error)
  call output_scalar(scalar_group_id, "pellet_phi", pellet_phi, ntime, error)
  call output_scalar(scalar_group_id, "pellet_z",   pellet_z,   ntime, error)
  call output_scalar(scalar_group_id, "pellet_velx", pellet_velx, ntime, error)
  call output_scalar(scalar_group_id, "pellet_velphi", pellet_velphi, ntime, error)
  call output_scalar(scalar_group_id, "pellet_velz", pellet_velz, ntime, error)
  call output_scalar(scalar_group_id, "pellet_var", pellet_var, ntime, error)
  call output_scalar(scalar_group_id, "r_p",        r_p,        ntime, error)
  call output_scalar(scalar_group_id, "r_p2",       r_p2,       ntime, error)
  call output_scalar(scalar_group_id, "pellet_rate", pellet_rate, ntime, error)
  call output_scalar(scalar_group_id, "pellet_rate1",  pellet_rate1, ntime, error)
  call output_scalar(scalar_group_id, "pellet_rate2", pellet_rate2, ntime, error)
  call output_scalar(scalar_group_id, "pellet_ablrate", pellet_ablrate, ntime, error)

  ! Controllers
  call output_scalar(scalar_group_id, "loop_voltage",        vloop,               ntime, error)
  call output_scalar(scalar_group_id, "i_control%err_i",     i_control%err_i,     ntime, error)
  call output_scalar(scalar_group_id, "i_control%err_p_old", i_control%err_p_old, ntime, error)
  call output_scalar(scalar_group_id, "n_control%err_i",     n_control%err_i,     ntime, error)
  call output_scalar(scalar_group_id, "n_control%err_p_old", n_control%err_p_old, ntime, error)


  ! Diagnostics (Not needed for restart)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  call output_scalar(scalar_group_id, "psi0", psi0, ntime, error)
  call output_scalar(scalar_group_id, "psimin"  ,psimin  ,ntime,error)
  call output_scalar(scalar_group_id, "temax"  ,temax  ,ntime,error)
  call output_scalar(scalar_group_id, "runaways", totre, ntime, error)
  call output_scalar(scalar_group_id, "psi_lcfs"        , psibound,ntime,error)

  call output_scalar(scalar_group_id, "area"            , area  , ntime, error)
  call output_scalar(scalar_group_id, "toroidal_flux"   , tflux , ntime, error)
  call output_scalar(scalar_group_id, "toroidal_current", totcur, ntime, error)
  call output_scalar(scalar_group_id, "particle_number" , totden, ntime, error)
  call output_scalar(scalar_group_id, "electron_number" , totne , ntime, error)
  call output_scalar(scalar_group_id, "angular_momentum", tmom  , ntime, error)
  call output_scalar(scalar_group_id, "circulation"     , tvor  , ntime, error)
  call output_scalar(scalar_group_id, "volume"          , volume, ntime, error)

  call output_scalar(scalar_group_id, "area_p"            , parea,ntime, error)
  call output_scalar(scalar_group_id, "toroidal_flux_p"   , pflux,ntime, error)
  call output_scalar(scalar_group_id, "toroidal_current_p", pcur ,ntime, error)
  call output_scalar(scalar_group_id, "particle_number_p" , pden ,ntime, error)
  call output_scalar(scalar_group_id, "angular_momentum_p", pmom ,ntime, error)
  call output_scalar(scalar_group_id, "volume_p"          , pvol ,ntime, error)
  call output_scalar(scalar_group_id, "M_IZ"              , m_iz ,ntime, error)
  call output_scalar(scalar_group_id, "M_IZ_co"           , m_iz_co ,ntime, error)
  call output_scalar(scalar_group_id, "M_IZ_sn"           , m_iz_sn ,ntime, error)
  call output_scalar(scalar_group_id, "IP_co"             , pcur_co ,ntime, error)
  call output_scalar(scalar_group_id, "IP_sn"             , pcur_sn ,ntime, error)

  call output_scalar(scalar_group_id, "toroidal_current_w",wallcur,ntime,error)
  
  call output_scalar(scalar_group_id, "E_MP" , emagp , ntime, error)
  call output_scalar(scalar_group_id, "E_KP" , ekinp , ntime, error)
  call output_scalar(scalar_group_id, "E_MPD", emagpd, ntime, error)
  call output_scalar(scalar_group_id, "E_KPD", ekinpd, ntime, error)
  call output_scalar(scalar_group_id, "E_MPH", emagph, ntime, error)
  call output_scalar(scalar_group_id, "E_KPH", ekinph, ntime, error)

  call output_scalar(scalar_group_id, "E_MT" , emagt , ntime, error)
  call output_scalar(scalar_group_id, "E_KT" , ekint , ntime, error)
  call output_scalar(scalar_group_id, "E_MTD", emagtd, ntime, error)
  call output_scalar(scalar_group_id, "E_KTD", ekintd, ntime, error)
  call output_scalar(scalar_group_id, "E_MTH", emagth, ntime, error)
  call output_scalar(scalar_group_id, "E_KTH", ekinth, ntime, error)

  call output_scalar(scalar_group_id, "E_P" , emag3, ntime, error)
  call output_scalar(scalar_group_id, "Ave_P" , avep, ntime, error)
  call output_scalar(scalar_group_id, "E_K3", ekin3, ntime, error)
  call output_scalar(scalar_group_id, "E_PD", emag3d, ntime, error)
  call output_scalar(scalar_group_id, "E_K3D", ekin3d, ntime, error)
  call output_scalar(scalar_group_id, "E_PH", emag3h, ntime, error)
  call output_scalar(scalar_group_id, "E_K3H", ekin3h, ntime, error)
  call output_scalar(scalar_group_id, "E_PE", w_pe ,ntime, error)

  call output_scalar(scalar_group_id, "Flux_pressure ", efluxp, ntime, error)
  call output_scalar(scalar_group_id, "Flux_kinetic  ", efluxk, ntime, error)
  call output_scalar(scalar_group_id, "Flux_poynting ", efluxs, ntime, error)
  call output_scalar(scalar_group_id, "Flux_thermal  ", efluxt, ntime, error)
  call output_scalar(scalar_group_id, "E_grav        ", epotg,  ntime, error)

  if(rad_source) then
     call output_scalar(scalar_group_id, "radiation"       , totrad, ntime, error)
  endif

  if(xray_detector_enabled.eq.1) then
     call output_scalar(scalar_group_id,"xray_signal",xray_signal,ntime,error)
  end if

  if(itaylor.eq.3) then
     temp = reconnected_flux()
     call output_scalar(scalar_group_id, "Reconnected_Flux", temp, ntime, error)
  endif

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
  call output_scalar(scalar_group_id, "Torque_parvisc",tau_parvisc,ntime,error)

  call output_scalar(scalar_group_id, "Wall_Force_n0_x", wall_force_n0_x,ntime,error)
  call output_scalar(scalar_group_id, "Wall_Force_n0_y", wall_force_n0_y,ntime,error)
  call output_scalar(scalar_group_id, "Wall_Force_n0_z", wall_force_n0_z,ntime,error)
  call output_scalar(scalar_group_id, "Wall_Force_n1_x", wall_force_n1_x,ntime,error)
  call output_scalar(scalar_group_id, "Wall_Force_n1_y", wall_force_n1_y,ntime,error)

  call output_scalar(scalar_group_id, "Parallel_viscous_heating",bwb2,ntime,error)


  ! Probes
  if(imag_probes.ne.0) then
     call output_1dextendarr(mp_group_id, "value", mag_probe_val, imag_probes, &
          ntime, error)
  end if

  if(iflux_loops.ne.0) then
     call output_1dextendarr(fl_group_id, "value", flux_loop_val, iflux_loops, &
          ntime, error)
  end if

  call h5gclose_f(scalar_group_id, error)
  if(imag_probes.ne.0) call h5gclose_f(mp_group_id, error)
  if(iflux_loops.ne.0) call h5gclose_f(fl_group_id, error)

  call h5gclose_f(root_id, error)

  if(myrank.eq.0 .and. iprint.ge.1) print *, 'Done writing scalars'

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
  call output_scalar(timing_group_id, "t_kprad"       , t_kprad       , ntime, error)

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
  integer(HID_T) :: root_id
  integer :: nelms

  character(LEN=19) :: time_file_name
  integer(HID_T) :: time_file_id, time_root_id, plist_id
  integer :: info

  call hdf5_get_local_elms(nelms, error)

  ! create the name of the group
  if(equilibrium.eq.1) then
     time_group_name = "equilibrium"
     time_file_name = "equilibrium.h5"
  else
     write(time_group_name, '("time_",I3.3)') times_output
     write(time_file_name, '("time_",I3.3,".h5")') times_output
  endif

  ! remove the time group link if it already exists
  ! from before a restart
  if(irestart.ne.0 .and. ntime.eq.ntime0) then
     call h5gunlink_f(file_id, time_group_name, error)
  endif


  ! Create new file for timeslice
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! Set up the file access property list with parallel I/O
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  info = MPI_INFO_NULL
  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, info, error)
     
  ! Open the new file
  call h5fcreate_f(time_file_name, H5F_ACC_TRUNC_F, time_file_id, error, &
       access_prp = plist_id)
  if(error.lt.0) then
     print *, "Error: could not open ", time_file_name, &
          " for HDF5 output.  error = ", error
     return
  endif
     
  ! open the root group
  call h5gopen_f(time_file_id, "/", time_root_id, error)
  
  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, ' Writing time slice file ', time_file_name
  
  ! Write attributes
  if(myrank.eq.0 .and. iprint.ge.1) print *, '  Writing attr '
  call write_real_attr(time_root_id, "time", time, error)
#ifdef USE3D
  call write_int_attr(time_root_id, "nspace", 3, error)
#else
  call write_int_attr(time_root_id, "nspace", 2, error)
#endif
  ! Time step associated with this time slice
  if(equilibrium.eq.1) then
     call write_int_attr(time_root_id, "ntimestep", 0, error)
  else
     call write_int_attr(time_root_id, "ntimestep", ntime, error)
  end if
  
  ! Write version number
  call write_int_attr(time_root_id, "version", version, error)
  
  ! Output the mesh data
  if(myrank.eq.0 .and. iprint.ge.1) print *, '  Writing mesh '
  call output_mesh(time_root_id, nelms, error)
  
  ! Output the field data 
  if(myrank.eq.0 .and. iprint.ge.1) print *, '  Writing fields '
  call output_fields(time_root_id, equilibrium, error)
  if(myrank.eq.0 .and. iprint.ge.1) print *, '  Done writing fields ', error
  
  ! Close the file
  call h5gclose_f(time_root_id, error)
  call h5fclose_f(time_file_id, error)  
  call h5pclose_f(plist_id, error)


  ! Add timeslice link in main file
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! create a link to the file
  if(myrank.eq.0) print *, 'linking ', time_file_name
  call h5gopen_f(file_id, "/", root_id, error)
  call h5lcreate_external_f(time_file_name, "/", root_id, time_group_name, &
       error)
  
  ! update number of time slices
  if(equilibrium.eq.0) times_output = times_output + 1
  call update_int_attr(root_id, "ntime", times_output, error)

  ! close root group
  call h5gclose_f(root_id, error)
  
  if(myrank.eq.0 .and. iprint.ge.1) print *, '  End of hdf5_write_time_slice '

end subroutine hdf5_write_time_slice


! output_mesh
! ===========
subroutine output_mesh(time_group_id, nelms, error)
  use hdf5
  use hdf5_output
  use mesh_mod
  use basic
  use boundary_conditions

  implicit none

  integer(HID_T), intent(in) :: time_group_id
  integer, intent(in) :: nelms
  integer, intent(out) :: error

  type(element_data) :: d
  integer(HID_T) :: mesh_group_id
  integer :: i
#ifdef USE3D
  integer, parameter :: vals_per_elm = 10
#else
  integer, parameter :: vals_per_elm = 8
#endif
  real, dimension(vals_per_elm,nelms) :: elm_data
  integer, dimension(nodes_per_element) :: nodeids
  real :: alx, alz

  integer :: is_edge(3)
  real :: normal(2,3)
  integer :: idim(3)
  real :: bound
  integer :: izone

  ! Create the group
  call h5gcreate_f(time_group_id, "mesh", mesh_group_id, error) 

  ! Write attributes
  call write_int_attr(mesh_group_id, "nelms", global_elms, error)
  call get_bounding_box_size(alx, alz)
  call write_real_attr(mesh_group_id, "width", alx, error)
  call write_real_attr(mesh_group_id, "height", alz, error)
#ifdef USE3D
  call write_int_attr(mesh_group_id, "3D", 1, error)
#else
  call write_int_attr(mesh_group_id, "3D", 0, error)
#endif
  call write_int_attr(mesh_group_id, "nplanes", nplanes, error)
  call write_real_attr(mesh_group_id, "period", toroidal_period, error)

  ! Output the mesh data
  do i=1, nelms
     call get_element_nodes(i,nodeids)

     ! don't call boundary_edge if iadapt != 0
     ! because bug in scorec software causes crash when querying 
     ! normal/curvature at newly created boundary nodes
     if(iadapt.eq.0) call boundary_edge(i, is_edge, normal, idim)

     bound = 0.
     if(is_edge(1).ne.0) bound = bound + 1. + (is_edge(1)-1)*2**3
     if(is_edge(2).ne.0) bound = bound + 2. + (is_edge(2)-1)*2**7
     if(is_edge(3).ne.0) bound = bound + 4. + (is_edge(3)-1)*2**11

     call get_element_data(i, d)

     call get_zone(i, izone)

     elm_data( 1,i) = d%a
     elm_data( 2,i) = d%b
     elm_data( 3,i) = d%c
     elm_data( 4,i) = atan2(d%sn,d%co)
     elm_data( 5,i) = d%R
     elm_data( 6,i) = d%Z
     elm_data( 7,i) = bound
     elm_data( 8,i) = izone
#ifdef USE3D
     elm_data( 9,i) = d%d
     elm_data(10,i) = d%Phi
#endif
  end do
  call output_field(mesh_group_id, "elements", elm_data, vals_per_elm, &
       nelms, error)

  ! Close the group
  call h5gclose_f(mesh_group_id, error)
end subroutine output_mesh


subroutine write_field(group_id, name, f, nelms, error, isreal)
  use hdf5
  use basic
  use element
  use field
  use hdf5_output

  implicit none
  integer(HID_T) :: group_id
  character(LEN=*) :: name
  type(field_type), intent(in) :: f
  integer, intent(in) :: nelms
  integer, intent(out) :: error
  logical, intent(in), optional :: isreal

  vectype, dimension(coeffs_per_element,nelms) :: dum
  integer :: i
  logical :: ir

  if(myrank.eq.0 .and. iprint.ge.2) print *, 'Writing ', name

  if(present(isreal)) then
     ir = isreal
  else
     ir = .false.
  end if

  error = 0

  do i=1, nelms
     call calcavector(i, f, dum(:,i))
  end do
  call output_field(group_id, name, real(dum), coeffs_per_element, &
       nelms, error)
#ifdef USECOMPLEX
  if(.not.ir) then
     call output_field(group_id,name//"_i",aimag(dum),coeffs_per_element, &
          nelms,error)
  end if
#endif
end subroutine write_field


! output_fields
! =============
subroutine output_fields(time_group_id, equilibrium, error)
  use hdf5
  use hdf5_output
  use basic
  use arrays
  use time_step
  use auxiliary_fields
  use transport_coefficients
  use kprad_m3dc1

  implicit none

  integer(HID_T), intent(in) :: time_group_id
  integer, intent(out) :: error
  integer, intent(in) :: equilibrium
  
  integer(HID_T) :: group_id
  integer :: i, nelms, ilin
  vectype, allocatable :: dum(:,:), dum2(:,:)
  character(len=64) :: field_name

  ilin = 1 - equilibrium

  if(myrank.eq.0 .and. iprint.ge.1) print *, 'In output_fields.  ilin=', ilin

  nelms = local_elements()
  error = 0
  
  allocate(dum(coeffs_per_element,nelms))

  ! Create the fields group
  call h5gcreate_f(time_group_id, "fields", group_id, error)

  ! Output the fields
  ! ~~~~~~~~~~~~~~~~~

  ! psi_plasma
  if(icsubtract.eq.1 .or. &
       (extsubtract.eq.1 .and. (ilin.eq.1 .or. eqsubtract.eq.0))) then
     call write_field(group_id, "psi_plasma", psi_field(ilin), nelms, error)
  end if

  ! psi
  do i=1, nelms
     call calcavector(i, psi_field(ilin), dum(:,i))
  end do
  if(icsubtract.eq.1 .and. (ilin.eq.0 .or. eqsubtract.eq.0)) then
     allocate(dum2(coeffs_per_element,nelms))
     do i=1, nelms
        call calcavector(i, psi_coil_field, dum2(:,i))
     end do
     dum = dum + dum2
     deallocate(dum2)
  endif
  if(extsubtract.eq.1 .and. (ilin.eq.1 .or. eqsubtract.eq.0)) then 
     allocate(dum2(coeffs_per_element,nelms))
     do i=1, nelms
        call calcavector(i, psi_ext, dum2(:,i))
     end do
     dum = dum + dum2
     deallocate(dum2)
  end if
  call output_field(group_id, "psi", real(dum), coeffs_per_element, &
       nelms, error)
#ifdef USECOMPLEX
  call output_field(group_id,"psi_i",aimag(dum),coeffs_per_element, &
       nelms,error)
#endif

  ! u
  call write_field(group_id, "phi", u_field(ilin), nelms, error)

  ! electrostatic potential
  if(jadv.eq.0) then
     call write_field(group_id, "potential", e_field(ilin), nelms, error)
  endif

  ! I
  do i=1, nelms
     call calcavector(i, bz_field(ilin), dum(:,i))
  end do
  if(extsubtract.eq.1 .and. (ilin.eq.1 .or. eqsubtract.eq.0)) then 
     call output_field(group_id, "I_plasma", real(dum), coeffs_per_element, &
          nelms, error)
#ifdef USECOMPLEX
     call output_field(group_id,"I_plasma_i",aimag(dum),coeffs_per_element, &
          nelms,error)
#endif

     allocate(dum2(coeffs_per_element,nelms))
     do i=1, nelms
        call calcavector(i, bz_ext, dum2(:,i))
     end do
     dum = dum + dum2
     deallocate(dum2)
  end if
  call output_field(group_id, "I", real(dum), coeffs_per_element, &
       nelms, error)
#ifdef USECOMPLEX
  call output_field(group_id,"I_i",aimag(dum),coeffs_per_element, &
       nelms,error)
#endif
    
    ! BF
  if(ifout.eq.1) then
     do i=1, nelms
        call calcavector(i, bf_field(ilin), dum(:,i))
     end do
     if(extsubtract.eq.1 .and. (ilin.eq.1 .or. eqsubtract.eq.0)) then 
        call output_field(group_id, "f_plasma", real(dum), coeffs_per_element, &
             nelms, error)
#ifdef USECOMPLEX
        call output_field(group_id,"f_plasma_i",aimag(dum),coeffs_per_element, &
             nelms,error)
#endif

        allocate(dum2(coeffs_per_element,nelms))
        do i=1, nelms
           call calcavector(i, bf_ext, dum2(:,i))
        end do
        dum = dum + dum2
        deallocate(dum2)
     end if
     call output_field(group_id, "f", real(dum), coeffs_per_element, &
          nelms, error)
#ifdef USECOMPLEX
     call output_field(group_id,"f_i",aimag(dum),coeffs_per_element, &
          nelms,error)
#endif
  endif
    
  call write_field(group_id, "V",    vz_field(ilin), nelms, error)
  call write_field(group_id, "Pe",   pe_field(ilin), nelms, error)
  call write_field(group_id, "P",     p_field(ilin), nelms, error)
  call write_field(group_id, "chi", chi_field(ilin), nelms, error)
  call write_field(group_id, "den", den_field(ilin), nelms, error)
  call write_field(group_id, "te",   te_field(ilin), nelms, error)
  call write_field(group_id, "ti",   ti_field(ilin), nelms, error)
  call write_field(group_id, "ne",   ne_field(ilin), nelms, error)
  
  if(icsubtract.eq.1) then
     call write_field(group_id, "psi_coil", psi_coil_field, nelms, error, .true.)
  end if

#ifdef USEPARTICLES
  if (kinetic.eq.1) then
     if (associated(p_i_perp%vec)) then
        !Perpendicular component of hot ion pressure tensor
        call write_field(group_id, "p_i_perp", p_i_perp, nelms, error)
     endif

     if (associated(p_i_par%vec)) then
        !Parallel component of hot ion pressure tensor
        call write_field(group_id, "p_i_par", p_i_par, nelms, error)
     endif
  endif
#endif

  if(use_external_fields) then 
     call write_field(group_id, "psi_ext", psi_ext, nelms, error)
     call write_field(group_id, "I_ext", bz_ext, nelms, error)    
     call write_field(group_id, "f_ext", bf_ext, nelms, error)
  endif

  if(ikprad.eq.1) then
     do i=0, kprad_z
        write(field_name, '(A,I2.2)') "kprad_n_", i
        call write_field(group_id, trim(field_name), kprad_n(i), nelms, error)
     end do
     call write_field(group_id, "kprad_sigma_e", kprad_sigma_e, nelms, error)
     call write_field(group_id, "kprad_sigma_i", kprad_sigma_i, nelms, error)
     call write_field(group_id, "kprad_rad", kprad_rad, nelms, error)
  end if

  ! transport coefficients do not change with time in linear calculations
  ! so don't output linear perturbations
  if(iwrite_transport_coeffs.eq.1 .and. &
     (linear.eq.0 .or. equilibrium.eq.1)) then

     call write_field(group_id, "eta", resistivity_field, nelms, error, .true.)
     call write_field(group_id, "visc", visc_field, nelms, error, .true.)
     call write_field(group_id, "visc_c", visc_c_field, nelms, error)
     call write_field(group_id, "kappa", kappa_field, nelms, error)
     
     ! poloidal force and mach number
     if(ipforce.gt.0) then
        call write_field(group_id, "pforce", pforce_field, nelms, error,.true.)
        call write_field(group_id, "pmach",  pmach_field,  nelms, error,.true.)
     endif
     
     if(ibootstrap.gt.0) then
        call write_field(group_id, "visc_e", visc_e_field, nelms, error,.true.)
     endif
  end if !(iwrite_transport_coeffs.eq.1)

  if(irunaway.ne.0) then
     call write_field(group_id, "n_re", nre_field, nelms, error)
  end if

  if(iwrite_aux_vars.eq.1) then 
    call write_field(group_id, "wall_dist", wall_dist, nelms, error, .true.)
    call write_field(group_id, "jphi", jphi_field, nelms, error)
    call write_field(group_id, "vor", vor_field, nelms, error)
    call write_field(group_id, "com", com_field, nelms, error)
    call write_field(group_id, "torque_em", torque_density_em, nelms, error)
    call write_field(group_id,"torque_ntv", torque_density_ntv, nelms, error)
    call write_field(group_id, "bdotgradp", bdotgradp, nelms, error)
    call write_field(group_id, "bdotgradt", bdotgradt, nelms, error)
    call write_field(group_id, "zeff", z_effective, nelms, error)
    if(ikprad.eq.1) then
       call write_field(group_id, "kprad_totden", kprad_totden, nelms, error)
    end if
    if(itemp_plot .eq. 1) then
       call write_field(group_id, "vdotgradt", vdotgradt, nelms, error)
       call write_field(group_id, "adv1", adv1, nelms, error)
       call write_field(group_id, "adv2", adv2, nelms, error)
       call write_field(group_id, "adv3", adv3, nelms, error)
       call write_field(group_id, "deldotq_perp", deldotq_perp, nelms, error)
       call write_field(group_id, "deldotq_par", deldotq_par, nelms, error)
       call write_field(group_id, "eta_jsq", eta_jsq, nelms, error, .true.)
       call write_field(group_id, "vpar", vpar_field, nelms, error)
       call write_field(group_id, "f1vplot", f1vplot, nelms, error)
       call write_field(group_id, "f1eplot", f1eplot, nelms, error)
       call write_field(group_id, "f2vplot", f2vplot, nelms, error)
       call write_field(group_id, "f2eplot", f2eplot, nelms, error)
       call write_field(group_id, "f3vplot", f3vplot, nelms, error)
       call write_field(group_id, "f3eplot", f3eplot, nelms, error)
       call write_field(group_id, "jdbobs", jdbobs, nelms, error)

       if(jadv.eq.0) then
          call write_field(group_id, "psidot", psidot, nelms, error, .true.)
          call write_field(group_id, "veldif", veldif, nelms, error, .true.)
          call write_field(group_id, "eta_jdb", eta_jdb, nelms, error, .true.)
          call write_field(group_id, "bdgp", bdgp, nelms, error, .true.)
          call write_field(group_id, "vlbdgp", vlbdgp, nelms, error, .true.)
       endif

    endif
    
    ! sigma
    if(density_source) then
       call write_field(group_id, "sigma", sigma_field, nelms, error, .true.)
    endif
    
    ! momentum source
    if(momentum_source) then
       call write_field(group_id, "force_phi", Fphi_field, nelms, error,.true.)
    endif
    
    ! heat source
    if(heat_source) then
       call write_field(group_id, "heat_source", Q_field, nelms, error, .true.)
    endif
    
    ! radiation source
    if(rad_source) then
       call write_field(group_id, "rad_source", rad_field, nelms, error,.true.)
    endif
    
    ! current drive source
    if(icd_source.gt.0) then
       call write_field(group_id, "cd_source", cd_field, nelms, error, .true.)
    endif
    
    if(xray_detector_enabled.eq.1) then 
       ! chord_mask
       do i=1, nelms
          call calcavector(i, chord_mask, dum(:,i))
       end do
       call output_field(group_id,"chord_mask",real(dum),coeffs_per_element,&
            nelms, error)
    end if

    ! magnetic region
    call write_field(group_id, "magnetic_region", mag_reg, nelms, error,.true.)

    ! mesh zone
    call write_field(group_id, "mesh_zone", mesh_zone, nelms, error, .true.)

    ! electric_field
    call write_field(group_id, "E_R", ef_r, nelms, error)
    call write_field(group_id, "E_PHI", ef_phi, nelms, error)
    call write_field(group_id, "E_Z", ef_z, nelms, error)
    call write_field(group_id, "E_par", ef_par, nelms, error)

    call write_field(group_id, "eta_J", eta_j, nelms, error, .true.)

 endif !(iwrite_aux_vars.eq.1)

!!$  if(equilibrium.eq.1) then
!!$     ! partition
!!$     dum = 0
!!$     dum(1,:) = myrank
!!$     call output_field(group_id, "part", real(dum), coeffs_per_element, &
!!$          nelms, error)
!!$  end if
     
  ! Close the mesh group
  call h5gclose_f(group_id, error)

  deallocate(dum)

end subroutine output_fields


! hdf5_write_keharmonics
! ==================
subroutine hdf5_write_keharmonics(error)
  use basic
  use diagnostics
  use hdf5_output

  implicit none

  integer, intent(out) :: error

  integer(HID_T) :: root_id, keharmonics_group_id
  real, allocatable :: dum(:)
  integer :: i

  if(.not.allocated(keharmonic)) return

  allocate(dum(NMAX+1))
  call h5gopen_f(file_id, "/", root_id, error)

  if(ntime.eq.0 .and. irestart.eq.0) then
     call h5gcreate_f(root_id, "keharmonics", keharmonics_group_id, error)
  else
     call h5gopen_f(root_id, "keharmonics", keharmonics_group_id, error)
  endif

  ! keharmonic
  do i = 0, NMAX
     dum(i+1) = keharmonic(i)
  enddo
  
  call output_1dextendarr(keharmonics_group_id, "keharmonics" , dum, NMAX+1, ntime, error)
  
  call h5gclose_f(keharmonics_group_id, error)
  call h5gclose_f(root_id, error)

  deallocate(dum)

end subroutine hdf5_write_keharmonics


! hdf5_write_bharmonics
! ==================
subroutine hdf5_write_bharmonics(error)
  use basic
  use diagnostics
  use hdf5_output

  implicit none

  integer, intent(out) :: error

  integer(HID_T) :: root_id, bharmonics_group_id
  real, allocatable :: dum(:)
  integer :: i

  if(.not.allocated(bharmonic)) return

  allocate(dum(BNMAX+1))
  call h5gopen_f(file_id, "/", root_id, error)

  if(ntime.eq.0 .and. irestart.eq.0) then
     call h5gcreate_f(root_id, "bharmonics", bharmonics_group_id, error)
  else
     call h5gopen_f(root_id, "bharmonics", bharmonics_group_id, error)
  endif
  if(myrank.eq.0 .and. iprint.ge.1) print *, error, 'before output_1dextendarr'

  ! bharmonic
  do i = 0, BNMAX
     dum(i+1) = bharmonic(i)
  enddo
  
  call output_1dextendarr(bharmonics_group_id, "bharmonics" , dum, BNMAX+1, ntime, error)
  if(myrank.eq.0 .and. iprint.ge.1) print *, error, 'after output_1dextendarr', error
  
  call h5gclose_f(bharmonics_group_id, error)
  call h5gclose_f(root_id, error)

  deallocate(dum)

end subroutine hdf5_write_bharmonics


! hdf5_write_kspits
! ==================
subroutine hdf5_write_kspits(error)
  use basic
  use diagnostics
  use hdf5_output
  use matrix_mod, ONLY : maxnumofsolves, kspits

  implicit none

!  integer :: i, maxnumofsolves
!  real, allocatable:: kspits(:)
  integer, intent(out) :: error
  integer(HID_T) :: root_id, kspits_group_id


!  maxnumofsolves=3
!  if(.not.allocated(kspits)) allocate(kspits(1:maxnumofsolves))
!  do i=1,maxnumofsolves
!  kspits(i)=i
!  enddo

  call h5gopen_f(file_id, "/", root_id, error)
  if(ntime.eq.0 .and. irestart.eq.0) then
     call h5gcreate_f(root_id, "kspits", kspits_group_id, error)
  else
     call h5gopen_f(root_id, "kspits", kspits_group_id, error)
  endif
  if(myrank.eq.0 .and. iprint.ge.1) print *, error, 'before output_1dextendarr'

  ! ksp iteration number for solve #5(velocity) #17(pressure) #6(field) stored in array kspits(3)
  ! kspits(1) : #5(velocity)
  ! kspits(2) : #1(pressure)
  ! kspits(3) : #17(pressure)
  ! kspits(4) : #6(field)
  call output_1dextendarr(kspits_group_id, "kspits" , kspits, maxnumofsolves, ntime, error)
  if(myrank.eq.0 .and. iprint.ge.1) print *, error, 'after output_1dextendarr', error

  call h5gclose_f(kspits_group_id, error)
  call h5gclose_f(root_id, error)
end subroutine hdf5_write_kspits



end module m3dc1_output
