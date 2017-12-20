Program Reducedquintic

!   Ref:  [1] Strang and Fix, An Analysis of the Finite Element Method, page 83
!         [2] G.R. Cowper, et al, AIAA Journal, Vol 7, no. 10 page 19
!         [3] JCP Vol 147 p 318-336 (1998)

  use basic
  use arrays
  use newvar_mod
  use sparse
  use hdf5_output
  use diagnostics
  use boundary_conditions
  use time_step
  use m3dc1_output
  use auxiliary_fields
  use pellet
  use scorec_mesh_mod
  use adapt
  use particles
  use math

  implicit none

#include "finclude/petsc.h"

  integer :: ier, i, adapt_flag
  real :: tstart, tend, dtsave, period
  character*10 :: datec, timec
  character*256 :: arg

  ! Initialize MPI
  call MPI_Init(ier)
  if (ier /= 0) then
     print *, 'Error in MPI_Init', ier
     call safestop(1)
  endif
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ier)
  if (ier /= 0) then
     print *,'Error in MPI_Comm_rank:',ier
     call safestop(1)
  endif
  call MPI_Comm_size(MPI_COMM_WORLD,maxrank,ier)
  if (ier /= 0) then
     print *,'Error in MPI_Comm_size:',ier
     call safestop(1)
  endif

  print_help = .false.
  do i=1, command_argument_count()
     call get_command_argument(i, arg)
     if(trim(arg) == '--help') print_help = .true.
  end do

  ! Write version information
  if(myrank.eq.0) then
     print *, '=============================================================='
     print *, 'RELEASE VERSION: ', RELEASE_VERSION
     print *, 'BUILD DATE: ', DATE_BUILT
!     print *, BUILD_INFO
     call date_and_time( datec, timec)
     write(*,1001) datec(1:4),datec(5:6),datec(7:8), &
          timec(1:2),timec(3:4),timec(5:8)
1001 format("START DATE: ", a4,1x,a2,1x,a2,3x,"TIME: ",a2,":",a2,":",a4,/)
#ifdef USECOMPLEX
     print *, 'COMPLEX VERSION'
#else
     print *, 'REAL VERSION'
#endif
#ifdef USE3D
     print *, '3D VERSION'
#else
     print *, '2D VERSION'
#endif
  endif

#ifdef USESCOREC
  call m3dc1_domain_init()
#endif

  ! Initialize PETSc
  if(myrank.eq.0) print *, 'Initializing PETSc...'
  call PetscInitialize(PETSC_NULL_CHARACTER, ier)
  if (ier /= 0) then
     print *,'Error in PetscInitialize:',ier
     call safestop(1)
  endif

  ! read input file
  if(myrank.eq.0) print *, ' Reading input'
  call input

  ! load mesh
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Loading mesh'
  call m3dc1_matrix_setassembleoption(imatassemble)
  if(itor.eq.0) then 
     period = twopi*rzero
  else
     period = twopi
  end if
  call load_mesh(period)

!  call print_node_data
!  call safestop(1)

  ! allocate arrays
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Allocating arrays'
  call space(1)

  sparse_initialized = .true.

  ! initialize variables
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Initializing variables'
  call init

  ! output info about simulation to be run
  call print_info

  ! create the newvar matrices
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Generating newvar matrices'
  call create_newvar_matrices


  ! Set initial conditions either from restart file
  ! or from initialization routine
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  select case (irestart)
  case(0)
     ! Initialize from routine

     ptot = 0.
     ntime = 0
     time = 0.

     if(myrank.eq.0 .and. iprint.ge.1) print *, ' Defining initial conditions'
     call initial_conditions

     ! correct for left-handed coordinates
     if(iflip.eq.1) then
        if(myrank.eq.0 .and. iprint.ge.1) &
             print *, ' Flipping coordinate system handedness'
        call flip_handedness
     endif

     ! combine the equilibrium and perturbed fields of linear=0
     ! unless eqsubtract = 1
     if(eqsubtract.eq.0) then
        call add(field_vec, field0_vec)
        field0_vec = 0.
     endif

     ! initialize feedback systems
     i_control%err_i = 0.
     i_control%err_p_old = 0.
     n_control%err_i = 0.
     n_control%err_p_old = 0.

  case(1)
!
!....save timestep from input file (needed if not a variable timestep run)
     dtsave = dt
     ! Read restart file(s)

     if(myrank.eq.0 .and. iprint.ge.1) print *, ' Reading restart file(s)'
     if(iglobalin.eq.1) then
        call rdrestartglobal
     else
#ifdef USEADIOS
        call rdrestart_adios
#else
        call rdrestart
#endif
     endif
!
!....use timestep from input file if not a variable timestep run
    if(dtkecrit.eq.0) dt = dtsave

  case(3)
!....read 2D real RL=1 restart file to start 2D complex COM=1 run
!
!....save timestep from input file
     dtsave = dt
        call rdrestart_cplx
     dt = dtsave

  end select                     !  end of the branch on restart/no restart

  ntime0 = ntime

  ! initialize output
  call initialize_output

  ! zero-out scalar data
  call reset_scalars


  if(itimer.eq.1) call reset_timings

  ! output simulation parameters and equilibrium
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(ntime.eq.0) then
     if(myrank.eq.0 .and. iprint.ge.1) &
          print *, " Writing simulation parameters"
     call hdf5_write_parameters(ier)

     if(eqsubtract.eq.1) then
        call derived_quantities(0)
     end if

     if(iwrite_aux_vars.eq.1) then
        call calculate_auxiliary_fields(0)
     end if

     ! Output the equilibrium
     call hdf5_write_time_slice(1,ier)
  end if

  ! Calculate all quantities derived from basic fields
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  call derived_quantities(1)


  ! Adapt the mesh
  ! ~~~~~~~~~~~~~~
#ifdef USESCOREC
  if (iadapt .eq. 1) then
   if(iprint.ge.1 .and. myrank.eq.0) write(*,*) "before adapt_by_psi call:  psibound, psimin", psibound, psimin
    call adapt_by_psi
  end if
#endif

  if(irestart.eq.0  .or. iadapt.gt.0) then
     tflux0 = tflux
     totcur0 = totcur
  endif

  ! output initial conditions
  call output

  ! if there are no timesteps to calculate, then skip time loop
  if(ntimemax.le.ntime) call safestop(0)

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Initializing timestep'
  call initialize_timestep


  if (0.eq.1) then
     call particle_test
     call safestop(0)
  endif

  ! main time loop
  ! ~~~~~~~~~~~~~~
  do ntime=ntime+1,ntimemax

     if(myrank.eq.0) print *, 'TIME STEP: ', ntime

     ! check for error
     if(ekin.ne.ekin .or. emag.ne.emag) then
        print *, "Error: energy is NaN"
        exit
     endif

     ! re-scale solution if energy is too large
     if(linear.eq.1) call scaleback


     ! take time step
     if(myrank.eq.0 .and. iprint.ge.1) print *, " Calling onestep"
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call onestep
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_onestep = t_onestep + tend - tstart
      write(*,1002) ntime, t_onestep, tend, tstart
 1002 format(" LOOP TIME", i5,1p3e16.8)
     endif

     if(linear.eq.0 .and. eqsubtract.eq.0 .and. i_control%icontrol_type .ge. 0) then
     ! feedback control on toroidal current
          if(myrank.eq.0 .and. iprint.ge.1) &
             print *, " Applying current feedback", &
             vloop, totcur, i_control%p, &
             i_control%target_val, i_control%err_p_old, i_control%err_i

          call control(totcur, vloop,       i_control, dt)

          if(myrank.eq.0 .and. iprint.ge.1) &
             print *, " After current feedback", &
             vloop, totcur, i_control%p, &
             i_control%target_val, i_control%err_p_old, i_control%err_i
     endif

     if(linear.eq.0 .and. eqsubtract.eq.0 .and. n_control%icontrol_type .ge. 0) then
     ! feedback control on density source
          if(myrank.eq.0 .and. iprint.ge.1) &
             print *, " Applying density feedback", &
             pellet_rate, totden, n_control%p, &
             n_control%target_val, n_control%err_p_old, n_control%err_i

          call control(totden, pellet_rate, n_control, dt)

          if(myrank.eq.0 .and. iprint.ge.1) &
             print *, " After density feedback", &
             pellet_rate, totden, n_control%p, &
             n_control%target_val, n_control%err_p_old, n_control%err_i
     endif

     ! Write output
     if(myrank.eq.0 .and. iprint.ge.1) print *, " Writing output."
     call output
     call run_adapt(adapt_flag)
     if(adapt_flag .eq. 1 .and. iadapt .gt. 1) call adapt_by_error
  enddo ! ntime

  if(myrank.eq.0 .and. iprint.ge.1) print *, "Done time loop."

  call safestop(0)

end Program Reducedquintic


!============================================================
! init
! ~~~~
! initialize variables
!============================================================
subroutine init
  use basic
  use mesh_mod
  
  implicit none
 
  rfac = (0,1)*ntor
  if(itor.eq.0) rfac = rfac/rzero

  ! define properties of triangles
  call tridef

  gbound = 0.

  ! Set up PID controllers
  i_control%p = control_p
  i_control%i = control_i
  i_control%d = control_d
  i_control%target_val = tcur
  i_control%icontrol_type = control_type

  n_control%p = n_control_p
  n_control%i = n_control_i
  n_control%d = n_control_d
  n_control%target_val = n_target
  n_control%icontrol_type = n_control_type

   if(itor.eq.0 .and. itaylor.ge.21) call init_qp
end subroutine init


!============================================
! print_info
! ~~~~~~~~~~
! print basic info about simulation options
!============================================
subroutine print_info
  use basic

  implicit none

  ! velocity form
  if(myrank.eq.0) then
     select case(ivform)
     case(0)
        print*, "V = grad(U)xgrad(phi) + V grad(phi) + grad(chi)"
     case(1)
        print*, "V = R^2 grad(U)xgrad(phi) + R^2 V grad(phi) + grad(chi)/R^2"
     end select
  endif

  ! check time-integration options
  select case(integrator)
  case(1)
     ! For BDF2 integration, first timestep is Crank-Nicholson with thimp=1,
     ! and subsequent timesteps are BDF2.
     if(myrank.eq.0) print *, "Time integration: BDF2."
  case default
     if(myrank.eq.0) print *, "Time integration: Crank-Nicholson."
  end select
end subroutine print_info

!=========================================
! safestop
! ~~~~~~~~
! stop program
!=========================================
subroutine safestop(iarg)

  use basic
  use sparse
  use m3dc1_output
  use time_step
  use auxiliary_fields

  implicit none
#include "finclude/petsc.h"
      
  integer, intent(in) :: iarg
  integer :: ier
  character*10 :: datec, timec

  call destroy_auxiliary_fields

  call finalize_timestep

  ! close hdf5 file
  if(myrank.eq.0 .and. iprint.ge.2) print *, "  finalizing output..."
  call finalize_output

  if(myrank.eq.0 .and. iprint.ge.2) print *, "  deleting matrices..."
  call delete_matrices

  if(myrank.eq.0 .and. iprint.ge.2) print *, "  unloading mesh..."
  call unload_mesh

  if(myrank.eq.0 .and. iprint.ge.2) print *, "  finalizing PETSC..."
  call PetscFinalize(ier)
  if(ier.ne.0) print *, 'Error finalizing PETSC:', ier

 ! Write time information
  if(myrank.eq.0) then
     print *, '=============================================================='
     call date_and_time( datec, timec)
     write(*,1002) datec(1:4),datec(5:6),datec(7:8), &
          timec(1:2),timec(3:4),timec(5:8)
1002 format(" END DATE: ", a4,1x,a2,1x,a2,3x,"TIME: ",a2,":",a2,":",a4,/)
  endif

  if(myrank.eq.0 .and. iprint.ge.2) print *, "  finalizing MPI..."
  call MPI_Finalize(ier)
  if (ier.ne.0) print *, 'Error finalizing MPI:', ier
  
  if(myrank.eq.0) print *, "Stopped at", iarg
  stop
end subroutine safestop


! ======================================================================
! smooth_velocity
! ~~~~~~~~~~~~~~~
!
! applies smoothing operators to velocity
! ======================================================================
subroutine smooth_velocity(uin, chiin)
  use basic
  use arrays
  use newvar_mod
  use diagnostics
  use sparse

  implicit none

  type(field_type), intent(inout) :: uin, chiin
  real :: tstart, tend
  type(vector_type) :: temp_vec
  type(field_type) :: vor_new

  if(hyperc.eq.0) return

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' smoothing velocity...'
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

  call create_vector(temp_vec,2)
  call associate_field(vor_new, temp_vec, 2)

  ! smooth vorticity
  call solve_newvar1(mass_mat_lhs_dc, vor_field, gs_mat_rhs_dc, uin)
  vor_new = vor_field
  call solve_newvar_axby(s5_mat,temp_vec,d5_mat,temp_vec)
  uin = vor_new
  
  ! smooth compression
  if(numvar.ge.3) then
     if(com_bc.eq.0) then
        call solve_newvar1(mass_mat_lhs,    com_field, lp_mat_rhs,    chiin)
     else
        call solve_newvar1(mass_mat_lhs_dc, com_field, lp_mat_rhs_dc, chiin)
     endif
     temp_vec = 0.
     vor_new = com_field
     call solve_newvar_axby(s7_mat,temp_vec,d7_mat,temp_vec)
     chiin = vor_new
  endif

  call destroy_vector(temp_vec)
  
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_smoother = t_smoother + tend - tstart
  endif
  
end subroutine smooth_velocity

! ======================================================================
! smooth_fields
! ~~~~~~~~~~~~~
!
! applies smoothing operators to fields
! ======================================================================
subroutine smooth_fields(psiin)
  use basic
  use arrays
  use newvar_mod
  use diagnostics
  use sparse

  implicit none

  type(field_type), intent(inout) :: psiin

  type(vector_type) :: temp_vec
  type(field_type) :: j_new, psi_new
  real :: tstart, tend

  if(jadv.eq.0 .or. hyper.eq.0. .or. (jadv.eq.1 .and. imp_hyper.ge.1)) return

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' smoothing fields...'
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

  ! smooth current density
  call create_vector(temp_vec,2)
  call associate_field(j_new, temp_vec, 1)
  call associate_field(psi_new, temp_vec, 2)

  j_new = 0.
  psi_new = psiin
  call solve_newvar_axby(s10_mat,temp_vec,d10_mat,temp_vec)
  psiin = psi_new

  call destroy_vector(temp_vec)

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_smoother = t_smoother + tend - tstart
  endif
     
end subroutine smooth_fields


! ======================================================================
! derived_quantities
! ~~~~~~~~~~~~~~~~~~
! calculates all derived quantities, including auxiliary fields
! and scalars
! ======================================================================
subroutine derived_quantities(ilin)
  use basic
  use arrays
  use newvar_mod
  use diagnostics
  use sparse
  use transport_coefficients
  use auxiliary_fields

  implicit none

  type(field_type) :: psi_temp, te_temp
  integer, intent(in) :: ilin    ! 0 for equilibrium fields, 1 for perturbed
  integer :: ier
  real :: tstart, tend

  vectype :: temp

  ! Find lcfs
  ! ~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.2) print *, "  finding lcfs"
  if(eqsubtract.eq.1) then
     if(linear.eq.1) then 
        if(ntime.eq.ntime0) call lcfs(psi_field(0))
     else
        call create_field(psi_temp)
        psi_temp = psi_field(0)
        call add_field_to_field(psi_temp, psi_field(1))
        call lcfs(psi_temp)
        call destroy_field(psi_temp)
     endif
  else
     call lcfs(psi_field(1))
  endif


  ! Find maximum temperature:  te_max
  ! ~~~~~~~~~
  ier = 0
  if(myrank.eq.0 .and. iprint.ge.2) print *, "  finding temax"
  if(eqsubtract.eq.1) then
     if(linear.eq.1) then 
        if(ntime.eq.ntime0) call te_max(xmag,zmag,te_field(0),temax,0,ier)
     else
        call create_field(te_temp)
        te_temp = te_field(0)
        call add_field_to_field(te_temp, te_field(1))
        call te_max(xmag,zmag,te_temp,temax,0,ier)
        call destroy_field(te_temp)
     endif
  else
     call te_max(xmag,zmag,te_field(1),temax,0,ier)
  endif
  if(ier.eq.0) then
    if(myrank.eq.0 .and. iprint.ge.1) then
      write(*,'(A, E12.4)') 'max te', temax
    endif
  else
    if(myrank.eq.0 .and. iprint.ge.1) then
      write(*,'(A,2e12.4)') ' no temperatue maximum found near ',xmag,zmag
    endif
  endif


  ! Define auxiliary fields
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

  if(myrank.eq.0 .and. iprint.ge.2) print *, "  transport coefficients"
  call define_transport_coefficients

  if(itemp.eq.0 .and. (numvar.eq.3 .or. ipres.gt.0)) then
     if(myrank.eq.0 .and. iprint.ge.2) print *, "  temperatures"
     call calculate_temperatures(ilin, te_field(ilin), ti_field(ilin), &
          eqsubtract)
  end if

  !   toroidal current
  if(myrank.eq.0 .and. iprint.ge.2) print *, "  toroidal current"
  if(inocurrent_tor.eq.1) then
     call solve_newvar1(mass_mat_lhs_dc,jphi_field,gs_mat_rhs_dc, &
          psi_field(ilin))
  else
     call solve_newvar1(mass_mat_lhs,jphi_field,gs_mat_rhs, &
          psi_field(ilin))
  endif

  if(hyperc.ne.0.) then
     !   vorticity
     if(myrank.eq.0 .and. iprint.ge.2) print *, "  vorticity"
     call solve_newvar1(mass_mat_lhs_dc,vor_field,gs_mat_rhs_dc,u_field(ilin))

     !   compression
     if(numvar.ge.3) then
        if(myrank.eq.0 .and. iprint.ge.2) print *, "  compression"
        if(com_bc.eq.1) then
           call solve_newvar1(mass_mat_lhs_dc,com_field,lp_mat_rhs_dc,& 
                chi_field(ilin))
        else
           call solve_newvar1(mass_mat_lhs,   com_field,lp_mat_rhs, &
                chi_field(ilin))
        endif
     else
        com_field = 0.
     endif
  endif

  ! vector potential stream function
  if(imp_bf.eq.0 .or. ilin.eq.0 .or. ntime.eq.0) then
     if((i3d.eq.1 .or. ifout.eq.1) .and. numvar.ge.2) then
        if(myrank.eq.0 .and. iprint.ge.2) print *, "  f"
        if(ilin.eq.0 .and. eqsubtract.eq.1) then
           if(itor.eq.0) then
              temp = bzero
           else
              temp = bzero*rzero
           end if
           call add(bz_field(ilin),-temp)
        endif
        call solve_newvar1(bf_mat_lhs,bf_field(ilin),mass_mat_rhs_bf, &
             bz_field(ilin), bf_field(ilin))
        if(ilin.eq.0 .and. eqsubtract.eq.1) call add(bz_field(ilin), temp)
     endif
  end if

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_aux = t_aux + tend - tstart
  endif

  ! calculate scalars
  ! ~~~~~~~~~~~~~~~~~
  if(icalc_scalars.eq.1 .and. ilin.eq.1) then
     if(myrank.eq.0 .and. iprint.ge.2) print *, "  scalars"
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call calculate_scalars
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_sources = t_sources + tend - tstart
     endif
  end if

end subroutine derived_quantities

!======================================================================
! conserve_flux
! ~~~~~~~~~~~~~
! adjusts the boundary conditions to conserve toroidal flux
!======================================================================
subroutine conserve_flux

  use basic
  use arrays
  use diagnostics

  implicit none
  
  ! adjust toroidal field and boundary condition to conserve toroidal flux
  if(numvar.ge.2 .and. iconstflux.eq.1) then
     gbound = (tflux0-tflux)/area
     call add(bz_field(1), gbound)
     if(myrank.eq.0) then
        print *, "Correction to toroidal flux: ", gbound*area
     end if
  endif

  return
end subroutine conserve_flux


!======================================================================
! flip_handedness
! ~~~~~~~~~~~~~~~
! Flips coordinate system handedness by flipping sign of
! psi, u, vz, and bz
!======================================================================
subroutine flip_handedness

  use basic
  use arrays
  use diagnostics

  implicit none

  integer :: ilin

  vectype, parameter :: temp = -1

  do ilin=0,1
     call mult(psi_field(ilin), temp)
     call mult(bz_field(ilin),  temp)
     call mult(u_field(ilin),   temp)
     call mult(vz_field(ilin),  temp)
     if(icsubtract.eq.1) call mult(psi_coil_field, temp)
  end do

  psimin = -psimin
  psilim = -psilim
  psibound = -psibound

end subroutine flip_handedness


!======================================================================
! magnetic_region
! ~~~~~~~~~~~~~~~
! determines what magnetic region the point x, z is in
! 0: inside plasma
! 1: scrape-off layer
! 2: private flux
!======================================================================
integer function magnetic_region(psi, x, z)
  use basic

  implicit none

  vectype, intent(in), dimension(dofs_per_node) :: psi
  real, intent(in) :: x, z 
  real :: psii, dpsii, pl, rl, al

  magnetic_region = 0

  dpsii = psibound - psimin

  psii = (real(psi(1)) - psimin)/dpsii
  if(psii .gt. 1.) then
     ! if Psi > 1, we are in scrape-off layer
     magnetic_region = 1
  else if(psii .lt. 0.) then
     magnetic_region = 0
  else
     ! if Psi < 1, but flux is increasing, we are in private flux region
     pl = sqrt(real(psi(2))**2 + real(psi(3))**2)
     rl = sqrt((x-xmag)**2 + (z-zmag)**2)
     if(pl.eq.0. .or. rl.eq.0.) return
     
     al = (real(psi(2))*(x-xmag) + real(psi(3))*(z-zmag))/(pl*rl)
     if(al*dpsii/abs(dpsii) .lt. 0.3) then
        magnetic_region = 2
     end if

     ! if z is far above or below x-point, we're in private flux region
     if(znull.ne.0. .and. xnull.gt.0.) then
        if((z-zmag)/(znull-zmag).gt.1.03) then
           magnetic_region = 2
        endif
     end if
     if(znull2.ne.0. .and. xnull2.gt.0.) then
        if((z-zmag)/(znull2-zmag).gt.1.03) then
           magnetic_region = 2
        endif
     end if
  end if
end function magnetic_region


!============================================================
! rotation
! ~~~~~~~~
! calculates the rotation matrix rot given angle theta
!============================================================
subroutine rotation(rot,ndim,theta)
  implicit none

  integer, intent(in) :: ndim
  real, intent(in) :: theta
  real, intent(out) :: rot(ndim,*)

  integer :: i, j
  real :: r1(6,6), co, sn

  co = cos(theta)
  sn = sin(theta)
  do i=1,6
     do j=1,6
        r1(i,j) = 0.
     enddo
  enddo

  r1(1,1) = 1.

  r1(2,2) = co
  r1(2,3) = sn

  r1(3,2) = -sn
  r1(3,3) = co

  r1(4,4) = co**2
  r1(4,5) = 2.*sn*co
  r1(4,6) = sn**2

  r1(5,4) = -sn*co
  r1(5,5) = co**2-sn**2
  r1(5,6) = sn*co

  r1(6,4) = sn**2
  r1(6,5) = -2.*sn*co
  r1(6,6) = co**2

  do i=1,18
     do j=1,18
        rot(i,j) = 0.
     enddo
  enddo

  do i=1,6
     do j=1,6
        rot(i,j)       = r1(i,j)
        rot(i+6,j+6)   = r1(i,j)
        rot(i+12,j+12) = r1(i,j)
     enddo
  enddo

  return
end subroutine rotation


  !============================================================
  ! tridef
  ! ~~~~~~
  ! populates the *tri arrays
  !============================================================
  subroutine tridef
    use basic
    use math
    use mesh_mod

    implicit none

    include 'mpif.h'
  
    type(element_data) :: d
    integer :: itri, i, j, k, ii, jj, numelms, numnodes, ndofs, ierr, m, n
    real, dimension(coeffs_per_tri,coeffs_per_tri) :: ti 
    real, dimension(dofs_per_tri, dofs_per_tri) :: rot, newrot
    real :: sum, theta, mean_area, tot_area, mean_len, temp
    real :: norm(2), curv, x, z
    integer :: inode(nodes_per_element)
    logical :: is_boundary
    integer :: izone, izonedim
    integer :: tot_elms
    real, dimension(dofs_per_tri) :: temp_vec, temp_vec2

    real, dimension(nodes_per_element) :: node_sz


    numelms = local_elements()
    numnodes = local_nodes()
    ndofs = numnodes*dofs_per_node

    ! start the loop over triangles within a rectangular region
    do itri=1,numelms

       ! define a,b,c and theta
       call get_element_data(itri,d)
       
       ! define the Inverse Transformation Matrix that enforces the 
       ! condition that the normal slope between triangles has only 
       ! cubic variation
       call tmatrix(ti,coeffs_per_tri,d%a,d%b,d%c)
       
       ! calculate the rotation matrix rot
       theta = atan2(d%sn,d%co)
       call rotation(rot,dofs_per_tri,theta)
       
       newrot = 0.
       call get_element_nodes(itri, inode)

       do i=1, 3
          call boundary_node(inode(i), &
               is_boundary, izone, izonedim, norm, curv, x, z, &
               all_boundaries)

          k = (i-1)*6 + 1
          if(is_boundary) then
             newrot(k  ,k  ) = 1.
             newrot(k+1,k+1) =  norm(1)
             newrot(k+1,k+2) =  norm(2)
             newrot(k+1,k+3) =  curv*norm(2)**2
             newrot(k+1,k+4) = -curv*norm(1)*norm(2)
             newrot(k+1,k+5) =  curv*norm(1)**2
             newrot(k+2,k+1) = -norm(2)
             newrot(k+2,k+2) =  norm(1)
             newrot(k+2,k+3) =  2.*curv*norm(1)*norm(2)
             newrot(k+2,k+4) = -curv*(norm(1)**2 - norm(2)**2) 
             newrot(k+2,k+5) = -2.*curv*norm(1)*norm(2)
             newrot(k+3,k+3) =  norm(1)**2 
             newrot(k+3,k+4) =  norm(1)*norm(2)
             newrot(k+3,k+5) =  norm(2)**2
             newrot(k+4,k+3) = -2.*norm(1)*norm(2)
             newrot(k+4,k+4) =  norm(1)**2 - norm(2)**2
             newrot(k+4,k+5) =  2.*norm(1)*norm(2)
             newrot(k+5,k+3) =  norm(2)**2
             newrot(k+5,k+4) = -norm(1)*norm(2)
             newrot(k+5,k+5) =  norm(1)**2
          else
             do j=1, 6
                newrot(k+j-1,k+j-1) = 1.
             end do
          end if
       end do     

       ! form the matrix g by multiplying ti and rot
       do k=1, coeffs_per_tri
          do j=1, dofs_per_tri
             sum = 0.
             do ii = 1, dofs_per_tri
                do jj=1, dofs_per_tri
                   sum = sum + newrot(j,jj)*ti(k,ii)*rot(ii,jj)
                end do
             enddo
             gtri(k,j,itri) = sum
          enddo
       enddo

       htri(1,1,itri) = 1.
#ifdef USE3D
       htri(2,1,itri) = 0.
       htri(3,1,itri) =-3./d%d**2
       htri(4,1,itri) = 2./d%d**3

       htri(1,2,itri) = 0.
       htri(2,2,itri) = 1.
       htri(3,2,itri) =-2./d%d
       htri(4,2,itri) = 1./d%d**2

       htri(1,3,itri) = 0.
       htri(2,3,itri) = 0.
       htri(3,3,itri) = 3./d%d**2
       htri(4,3,itri) =-2./d%d**3

       htri(1,4,itri) = 0.
       htri(2,4,itri) = 0.
       htri(3,4,itri) =-1./d%d
       htri(4,4,itri) = 1./d%d**2
#endif
       if(iprecompute_metric.eq.1) then
          call local_coeff_vector(itri,ctri(:,:,itri))
       end if

    end do

    select case(equilibrate)
    case(1)
       tot_area = 0.
       do itri=1, numelms
          call get_element_data(itri,d)
          tot_area = tot_area + (d%c)*(d%a + d%b)/2.
       end do

       call mpi_allreduce(numelms, tot_elms, 1, MPI_INTEGER, &
            MPI_SUM, MPI_COMM_WORLD, ierr)      
       call mpi_allreduce(tot_area, mean_area, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, MPI_COMM_WORLD, ierr)
       if(nplanes.le.1) then 
          mean_len = 1.
       else
          mean_len = 2.*pi/(nplanes-1)
       endif
       
       if(myrank.eq.0 .and. iprint.ge.1) then 
          print *, ' Total mesh area: ', mean_area
          print *, ' Total elements: ', tot_elms
       endif
       
       mean_area = mean_area/tot_elms
    
       if(myrank.eq.0 .and. iprint.ge.1) then 
          print *, ' Average area: ', mean_area
          print *, ' 1/mean_area = ', 1./mean_area
          print *, ' 1/mean_area**2 = ', 1./mean_area**2
       end if
       
       do i=1, pol_nodes_per_element
          equil_fac(1+dofs_per_node*(i-1),:) = 1./mean_area
          equil_fac(2+dofs_per_node*(i-1),:) = 1./sqrt(mean_area**3)
          equil_fac(3+dofs_per_node*(i-1),:) = 1./sqrt(mean_area**3)
          equil_fac(4+dofs_per_node*(i-1),:) = 1./mean_area**2
          equil_fac(5+dofs_per_node*(i-1),:) = 1./mean_area**2
          equil_fac(6+dofs_per_node*(i-1),:) = 1./mean_area**2
#ifdef USE3D
          do j=1, 6
             equil_fac(j+6+dofs_per_node*(i-1),:) = &
                  equil_fac(j+dofs_per_node*(i-1),:)/mean_len**2
          end do
#endif
       end do
       
    case(2)

       temp_vec = 0.
       do itri=1, numelms
#ifdef USESCOREC
          !call getelmsizes(itri, node_sz)
          node_sz = 1.
#else
          node_sz = 1.
#endif
          if(myrank.eq.0 .and. itri.eq.1) print *, 'node_sz = ', node_sz

          do i=1, nodes_per_element
             equil_fac(1+dofs_per_node*(i-1),:) = 1./node_sz(i)**2
             equil_fac(2+dofs_per_node*(i-1),:) = 1./node_sz(i)**3
             equil_fac(3+dofs_per_node*(i-1),:) = 1./node_sz(i)**3
             equil_fac(4+dofs_per_node*(i-1),:) = 1./node_sz(i)**4
             equil_fac(5+dofs_per_node*(i-1),:) = 1./node_sz(i)**4
             equil_fac(6+dofs_per_node*(i-1),:) = 1./node_sz(i)**4
          end do
          if(myrank.eq.0 .and. itri.eq.1) print *, equil_fac(:,itri)
       end do
    end select
  end subroutine tridef

!============================================================
! space
! ~~~~~
! allocates space for big arrays
!
! ifirstcall = 1 if this is the first call and 0 otherwise
!============================================================
subroutine space(ifirstcall)

  use element
  use mesh_mod
  use basic
  use arrays
  use sparse
  use time_step
  use auxiliary_fields
  use transport_coefficients

  implicit none

  integer, intent(in) :: ifirstcall

  integer :: numelms

#ifdef USESCOREC
  integer :: i, maxdofs
  character(len=32) :: field_name
#endif

  if(myrank.eq.0 .and. iprint.ge.1) print *, " Entering space..."

!.....create numberings
#ifdef USESCOREC
  if(ifirstcall .eq. 1) then
     do i=1, num_fields
       write(field_name,"(I2,A)")  i,0
#ifdef USECOMPLEX
       call m3dc1_field_create (i, trim(field_name), i, 1, dofs_per_node)
#else
       call m3dc1_field_create (i, trim(field_name), i, 0, dofs_per_node)
#endif
     end do
  endif ! on firstcall
#endif
  
  numelms = local_elements()

! arrays defined at all vertices
! createvec will delete the arrays if they have already been allocated
  if(ifirstcall.eq.1) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, 'Allocating...'

     ! Physical Variables
     call create_vector(field_vec , num_fields)
     call create_vector(field0_vec, num_fields)
     !if(iadapt .ne. 0) then
        call create_vector(field_vec_pre, 2)
     !end if

     ! Auxiliary Variables
     call create_field(jphi_field)
     call create_field(vor_field)
     call create_field(com_field)
     call create_field(resistivity_field)
     call create_field(kappa_field)
     call create_field(visc_field)
     call create_field(visc_c_field)
     if(ipforce.gt.0) call create_field(pforce_field)
     if(ipforce.gt.0) call create_field(pmach_field)
     if(density_source) call create_field(sigma_field)
     if(momentum_source) call create_field(Fphi_field)
     if(heat_source) call create_field(Q_field)
     if(icd_source.gt.0) call create_field(cd_field)
     call create_field(bf_field(0))
     call create_field(bf_field(1))
     if(ibootstrap.gt.0) call create_field(visc_e_field)

     call create_field(psi_coil_field)

     call create_auxiliary_fields
  endif

  ! arrays associated with the triangles
  if(ifirstcall.eq.0) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, ' deallocating...'
     deallocate(gtri,htri)
     if(iprecompute_metric.eq.1) deallocate(ctri)
     if(equilibrate.ne.0) deallocate(equil_fac)
  endif
  
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Allocating tri...'
  allocate(gtri(coeffs_per_tri,dofs_per_tri,numelms))
  allocate(htri(coeffs_per_dphi,dofs_per_dphi,numelms))
  if(iprecompute_metric.eq.1) &
       allocate(ctri(dofs_per_element,coeffs_per_element,numelms))
  if(equilibrate.ne.0) allocate(equil_fac(dofs_per_element,numelms))

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' associating...'
  call associate_field(u_field(1),   field_vec, u_g)
  call associate_field(vz_field(1),  field_vec, vz_g)
  call associate_field(chi_field(1), field_vec, chi_g)
  call associate_field(psi_field(1), field_vec, psi_g)
  call associate_field(bz_field(1),  field_vec, bz_g)
  call associate_field(pe_field(1),  field_vec, pe_g)
  call associate_field(den_field(1), field_vec, den_g)
  call associate_field(p_field(1),   field_vec, p_g)
  call associate_field(te_field(1),  field_vec, te_g)
  call associate_field(ti_field(1),  field_vec, ti_g)
  call associate_field(e_field(1),  field_vec, e_g)

  call associate_field(u_field(0),   field0_vec, u_g)
  call associate_field(vz_field(0),  field0_vec, vz_g)
  call associate_field(chi_field(0), field0_vec, chi_g)
  call associate_field(psi_field(0), field0_vec, psi_g)
  call associate_field(bz_field(0),  field0_vec, bz_g)
  call associate_field(pe_field(0),  field0_vec, pe_g)
  call associate_field(den_field(0), field0_vec, den_g)
  call associate_field(p_field(0),   field0_vec, p_g)
  call associate_field(te_field(0),  field0_vec, te_g)
  call associate_field(ti_field(0),  field0_vec, ti_g)

  call allocate_kspits

  !if (iadapt .ne. 0)  then
     call associate_field(u_field_pre,   field_vec_pre, u_g)
     call associate_field(psi_field_pre, field_vec_pre, psi_g)
  !end if


  if(myrank.eq.0 .and. iprint.ge.1) print *, " Exiting space."

  return
end subroutine space

! subroutine adapt_mesh
!   use basic
!   use arrays
!   use mesh_mod

!   implicit none

!   integer :: izone, izonedim, inode(nodes_per_element), numelms, itri, i
!   integer :: numnodes
!   real :: x, phi, z
!   vectype, dimension(dofs_per_node) :: dat
!   integer :: magnetic_region

! #ifdef USESCOREC
!   select case(iadapt)
!   case(1)
!      print *, 'adapting mesh...'
! !     call hessianadapt(resistivity_field%vec%data,1, 0, ntime, &
! !          adapt_factor, adapt_hmin, adapt_hmax)
!      print *, 'done adapting.'
!      call safestop(0)
     
!   case(2)
!     call create_field(temporary_field)
!     if(eqsubtract.eq.1) then
!        temporary_field = psi_field(0)
!     else
!        temporary_field = psi_field(1)
!     end if
!     if(icsubtract.eq.1) call add(temporary_field, psi_coil_field)
!     call straighten_field(temporary_field)

!     numnodes = local_nodes()
!     do i=1, numnodes
!        call get_node_pos(i, x, phi, z)
!        call get_node_data(temporary_field, i, dat)

!        if(magnetic_region(dat,x,z).eq.2) then
!           dat = 2. - (dat - psimin) / (psibound - psimin)
!           dat = (psibound - psimin)*dat + psimin          
!        end if
       
!        call set_node_data(temporary_field,i,dat)
!     end do
!     call sum_shared(temporary_field%vec)

!     if((adapt_psin_vacuum.ne.0 .or. adapt_psin_wall.ne.0) &
!          .and. imulti_region.eq.1) then
!        dat = 0.
!        numelms = local_elements()
!        do itri=1, numelms
!           call zonfac(itri,izone,izonedim)
!           if(izone.ge.3 .and. adapt_psin_vacuum.ne.0) then
!              call nodfac(itri,inode)
!              do i=1,3
!                 dat(1) = (psibound - psimin)*adapt_psin_vacuum + psimin
!                 call set_node_data(temporary_field,inode(i),dat)
!              end do
!           end if
!           if(izone.eq.2 .and. adapt_psin_wall.ne.0) then
!              call nodfac(itri,inode)
!              do i=1,3
!                 dat(1) = (psibound - psimin)*adapt_psin_wall + psimin
!                 call set_node_data(temporary_field,inode(i),dat)
!              end do
!           end if
!        end do
!        call sum_shared(temporary_field%vec)
!     end if

! !    if (iprint.ge.2) then
! !    write(25,1003) temporary_field%vec%data
! !1003 format(1p10E12.4)
! !    end if

!     print *, 'initializing solution transfer...'
!     call initsolutiontransfer(0)
!     print *, 'setting smoothing factor...', adapt_smooth
!     call setsmoothfact(adapt_smooth)
!     print *, 'adapting mesh...', psimin, psibound
! !    call adapt(temporary_field%vec%data,psimin,psibound)
!     print *, 'done adapting.'
!     call destroy_field(temporary_field)
!     call safestop(0)

!   end select
! #endif
! end subroutine adapt_mesh
