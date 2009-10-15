Program Reducedquintic

!   Ref:  [1] Strang and Fix, An Analysis of the Finite Element Method, page 83
!         [2] G.R. Cowper, et al, AIAA Journal, Vol 7, no. 10 page 19
!         [3] JCP Vol 147 p 318-336 (1998)

  use p_data
  use t_data
  use basic
  use arrays
  use newvar_mod
  use sparse
  use hdf5_output
  use diagnostics

  implicit none
#ifdef _AIX
  include 'mpif.h'
#endif
#include "finclude/petsc.h"

  integer :: ier
  real :: tstart, tend
  vectype :: temp

  PetscTruth :: flg
  PetscInt :: mpetscint,npetscint

  ! Start up message passing, SuperLU process grid
  call MPI_Init(ier)
  if (ier /= 0) then
     print *,'Error initializing MPI:',ier
     call safestop(1)
  endif
  call PetscInitialize(PETSC_NULL_CHARACTER, ier)
  call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-m',mpetscint,flg,ier)
  call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-n',npetscint,flg,ier)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ier)
  if (ier /= 0) then
     print *,'Error in MPI_Comm_rank:',ier
     call safestop(1)
  endif

  if(myrank.eq.0) then
     call date_and_time( datec, timec)
     write(*,1001) datec(1:4),datec(5:6),datec(7:8), &
          timec(1:2),timec(3:4),timec(5:8)
1001 format("M3D-C1 DATE: ", a4,1x,a2,1x,a2,3x,"TIME: ",a2,":",a2,":",a4,/)
#ifdef USECOMPLEX
     write(*,*) 'COMPLEX VERSION'
#else
     write(*,*) 'REAL VERSION'
#endif
  endif

  ! initialize the SUPERLU process grid
  if(myrank.eq.0) print *, 'Setting up SuperLU process grid'
  call initsolvers
  call MPI_Comm_size(MPI_COMM_WORLD,maxrank,ier)
  if (ier /= 0) then
     print *,'Error in MPI_Comm_size:',ier
     call safestop(1)
  endif
  ! initialize autopack
  if(myrank.eq.0) print *, 'Initializing Autopack'
  call AP_INIT()

  if(myrank.eq.0) print *, 'Loading Mesh'
  call loadmesh("struct.dmg", "struct-dmg.sms")
! call loaddebugger()

  ! read input file
  call input

  ! allocate arrays
  call space(1)

  ! initialize variables
  call init

  ! output info about simulation to be run
  call print_info

  ! create the newvar matrices
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.1) print *, "Generating newvar matrices..."
  call create_matrix(mass_matrix_lhs_dc,    NV_DCBOUND, NV_I_MATRIX,  NV_LHS)
  call create_matrix(mass_matrix_lhs,       NV_NOBOUND, NV_I_MATRIX,  NV_LHS)
  call create_matrix(gs_matrix_rhs_dc,      NV_DCBOUND, NV_GS_MATRIX, NV_RHS)
  if(hyperc.ne.0) then
     call create_matrix(s5matrix_sm, NV_SVBOUND, NV_SV_MATRIX, NV_LHS)
     call create_matrix(d5matrix_sm, NV_SVBOUND, NV_SV_MATRIX, NV_RHS)
     if(numvar.ge.3) then
        call create_matrix(s7matrix_sm, NV_SCBOUND, NV_SC_MATRIX, NV_LHS)
        call create_matrix(d7matrix_sm, NV_SCBOUND, NV_SC_MATRIX, NV_RHS)
        if(com_bc.eq.0) then
           call create_matrix(lp_matrix_rhs,   NV_NOBOUND,NV_LP_MATRIX,NV_RHS)
        else
           call create_matrix(lp_matrix_rhs_dc,NV_DCBOUND,NV_LP_MATRIX,NV_RHS)
        endif
     end if
  endif
  if((i3d.eq.1 .or. ifout.eq.1) .and. numvar.ge.2) then
     call create_matrix(mass_matrix_rhs_nm, NV_NMBOUND, NV_I_MATRIX,  NV_RHS)
     call create_matrix(bf_matrix_lhs_nm,   NV_NMBOUND, NV_BF_MATRIX, NV_LHS)
     call create_matrix(bf_matrix_rhs,      NV_NOBOUND, NV_BF_MATRIX, NV_RHS)
  endif
  if(jadv.eq.1 .and. hyper.ne.0.) then
     call create_matrix(s10matrix_sm, NV_SJBOUND, NV_SJ_MATRIX, NV_LHS)
     call create_matrix(d10matrix_sm, NV_SJBOUND, NV_SJ_MATRIX, NV_RHS)
  endif
  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, "Done generating newvar matrices."


  ! Set initial conditions either from restart file
  ! or from initialization routine
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  select case (irestart)
  case(1)
     ! Read restart file(s)

     if(myrank.eq.1 .and. iprint.ge.1) print *, 'Reading restart file(s)...'
     if(iglobalin.eq.1) then
        call rdrestartglobal
     else
        call rdrestart
     endif

  case(0)
     ! Initialize from routine

     ptot = 0.
     ntime = 0
     time = 0.
     if(myrank.eq.0 .and. iprint.ge.1) &
          print *, 'defining initial conditions...'
     call initial_conditions
     if(myrank.eq.0 .and. iprint.ge.1) print *, 'done initial conditions'

     if(eqsubtract.eq.1) then
        fieldi = field
     else
        fieldi = field0
     endif

     ! correct for left-handed coordinates
     if(iflip.eq.1) call flip_handedness

     ! combine the equilibrium and perturbed fields of linear=0
     ! unless eqsubtract = 1
     if(eqsubtract.eq.0) then
        field = field + field0
        field0 = 0.
     endif

     ! initialize t(n-1) values    
     phiold = phi

     if(isplitstep.eq.1) then
        velold = vel
        if(idens.eq.1) denold = den
        if(ipres.eq.1) presold = pres
     endif
  case(2)
     ! Read restart file(s) and use these to initialize grad-shafranov solve

     if(myrank.eq.1 .and. iprint.ge.1) print *, 'Reading restart file(s)...'
     if(iglobalin.eq.1) then
        call rdrestartglobal
     else
        call rdrestart
     endif
     ptot = 0.
     ntime = 0
     time = 0.
     if(myrank.eq.0 .and. iprint.ge.1) &
          print *, 'defining initial conditions...'
     call initial_conditions
     if(myrank.eq.0 .and. iprint.ge.1) print *, 'done initial conditions'

     if(eqsubtract.eq.1) then
        fieldi = field
     else
        fieldi = field0
     endif

     if(iflip.eq.1) call flip_handedness

     ! combine the equilibrium and perturbed fields of linear=0
     ! unless eqsubtract = 1
     if(eqsubtract.eq.0) then
        field = field + field0
        field0 = 0.
     endif

     ! initialize t(n-1) values
     phiold = phi

     if(isplitstep.eq.1) then
        velold = vel
        if(idens.eq.1) denold = den
        if(ipres.eq.1) presold = pres
     endif  
  end select                     !  end of the branch on restart/no restart

  ntime0 = ntime

  ! initialize hdf5
  if(myrank.eq.0 .and. iprint.ge.1) print *, "Initializing HDF5."
  call hdf5_initialize(ier)
  if(ier.lt.0) then 
     print *, "Error initializing HDF5"
     call safestop(5)
  end if
  if(myrank.eq.0 .and. iprint.ge.1) print *, "Done initializing HDF5."

  if(itimer.eq.1) call reset_timings


  ! output simulation parameters and equilibrium
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(irestart.eq.0  .or. iadapt.gt.0) then
     if(myrank.eq.0 .and. iprint.ge.1) &
          print *, "Writing simulation parameters."
     call hdf5_write_parameters(ier)

     if(eqsubtract.eq.1) then
        temp = bzero*rzero
        call scalar_operation(field0,bz_g,num_fields,OP_PLUS, -temp)
        call derived_quantities(field0, bf0)
        call scalar_operation(field0,bz_g,num_fields,OP_PLUS,  temp)
     end if

     ! Output the equilibrium
     if(eqsubtract.eq.1) call hdf5_write_time_slice(1,ier)
  end if

  ! Calculate all quantities derived from basic fields
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  call derived_quantities(field, bf)


  ! Adapt the mesh
  ! ~~~~~~~~~~~~~~
  select case(iadapt)
  case(1)
     print *, 'adapting mesh...'
     call hessianadapt(resistivity,1, 0, ntime, &
          adapt_factor, adapt_hmin, adapt_hmax)
     print *, 'done adapting.'
     call space(0)
     call tridef
     
!     call updatenormalcurvature
     call write_normlcurv
     
  case(2)
    call createvec(temporary_vector,1)
    call copyvec(field0, psi_g, num_fields, temporary_vector, 1, 1)
    print *, 'adapting mesh...'
    call adapt(temporary_vector,psimin,psibound)
    print *, 'done adapting.'
    call space(0)
    call tridef

!    call updatenormalcurvature
    call write_normlcurv
  end select

  if(irestart.eq.0  .or. iadapt.gt.0) then
     tflux0 = tflux
     totcur0 = totcur
  endif

  ! output initial conditions
  call output

  ! if there are no timesteps to calculate, then skip time loop
  if(ntimemax.le.ntime) go to 101


  ! main time loop
  ! ~~~~~~~~~~~~~~
  do ntime=ntime+1,ntimemax

     if(myrank.eq.0) print *, 'TIME STEP: ',ntime

     ! check for error
     if(ekin.ne.ekin .or. emag.ne.emag) then
        print *, "Error: energy is NaN"
        goto 101
     endif

     ! re-scale solution if energy is too large
     if(linear.eq.1) call scaleback

     if(myrank.eq.0 .and. iprint.ge.1) print *, "ntime = ", ntime

     if(myrank.eq.0 .and. iprint.ge.1) print *, "Before onestep"
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call onestep
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_onestep = t_onestep + tend - tstart
     endif
     if(myrank.eq.0 .and. iprint.ge.1) print *, "After onestep"


     ! feedback control on toroidal current
!!$     if(itor.eq.1 .and. itaylor.eq.1) call control_pid

     if(myrank.eq.0 .and. iprint.ge.1) print *, "Applying feedback.."
     call control(totcur, vloop,       i_control, dt)
     call control(totden, pellet_rate, n_control, dt)
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Done applying feedback."

     ! Write output
     call output
  enddo ! ntime


101 continue

  call safestop(2)

end Program Reducedquintic


!============================================================
! init
! ~~~~
! initialize variables
!============================================================
subroutine init
  use basic
  use t_data
  
  implicit none

  if(myrank.eq.0 .and. iprint.ge.1) print *, " Entering init..."
  
  ! define properties of triangles
  call tridef

  ! define 1/r vector
  if(itor.eq.1) call rinvdef(rinv)


  fbound = 0.
  gbound = 0.


  ! Set up PID controllers
  i_control%p = control_p
  i_control%i = control_i
  i_control%d = control_d
  i_control%target_val = tcur

  n_control%p = n_control_p
  n_control%i = n_control_i
  n_control%d = n_control_d
  n_control%target_val = n_target


  if(myrank.eq.0 .and. iprint.ge.1) print *, " Exiting init."

  return
end subroutine init


subroutine print_info
  use basic

  implicit none

  integer :: ndofs, numelms, numnodes, j


  if(myrank.eq.0) then
     select case(ivform)
     case(0)
        print*, "V = grad(U)xgrad(phi) + V grad(phi) + grad(chi)"
     case(1)
        print*, "V = R^2 grad(U)xgrad(phi) + R^2 V grad(phi) + grad(chi)/R^2"
     end select
  endif

  ! Output information about local dofs, nodes, etc.
  if(iprint.ge.1) then
     if(myrank.eq.0) then
        call numglobaldofs(1, ndofs)
        print *, 'global dofs = ', ndofs
     end if

     call numfac(numelms)
     call numnod(numnodes)
     call numprocdofs(1, j)
     call numdofs(1, ndofs)
     print *, 'proc, owned dofs, needed dofs',myrank, j, ndofs
     print *, 'proc, numnodes, numfaces', myrank, numnodes,numelms
  endif

!!$  if(nonrect.eq.1) call updatenormalcurvature
!!$  call write_normlcurv

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


! Stop program
subroutine safestop(iarg)

  use basic
  use sparse
  use hdf5_output

  implicit none
      
  integer, intent(in) :: iarg

  include 'mpif.h'

  integer :: ier
      
  ! close hdf5 file
  print *, "finalizing hdf5..."
  call hdf5_finalize(ier)
  print *, "done."
  
  print *, "finalizing SCOREC software..."
  call delete_matrices()
  call deletesearchstructure()
  call clearscorecdata()
  print *, "done."
  
  call finalizesolvers
  call PetscFinalize(ier)
  call MPI_Finalize(ier)
  if (ier.ne.0) print *,'Error terminating MPI:',ier
  
  write(*,*) "stopped at", iarg
  stop
end subroutine safestop


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

  implicit none

  integer :: ier
  real :: tstart, tend

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  if(myrank.eq.0 .and. iprint.ge.1) print *, "Writing output"
  if(myrank.eq.0 .and. iprint.ge.1) print *, "-Scalars (HDF5)."
  call hdf5_write_scalars(ier)
  if(mod(ntime-ntime0,ntimepr).eq.0) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, "-Time slice (HDF5)."
     call hdf5_write_time_slice(0,ier)
     if(myrank.eq.0 .and. iprint.ge.1) print *, "-Restart file(s)."
     if(iglobalout.eq.1) then
        call wrrestartglobal
     else
        call wrrestart
     endif
  endif
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_output_hdf5 = t_output_hdf5 + tend - tstart
  end if
  
  if(itimer.eq.1) then
     if(myrank.eq.0) call second(tstart)
     if(myrank.eq.0 .and. iprint.ge.1) print *, "-Timings (HDF5)."
     call hdf5_write_timings(ier)
     call reset_timings
  end if
  
  ! flush hdf5 data to disk
  if(myrank.eq.0 .and. iprint.ge.1) print *, "Flushing data to HDF5 file."
  call hdf5_flush(ier)
  if(myrank.eq.0 .and. iprint.ge.1) print *, "Done output."
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_output_hdf5 = t_output_hdf5 + tend - tstart
  end if

end subroutine output


! ======================================================================
! smooth_velocity
! ~~~~~~~~~~~~~~~
!
! applies smoothing operators to velocity
! ======================================================================
subroutine smooth_velocity(vec)
  use basic
  use arrays
  use newvar_mod
  use diagnostics
  use sparse

  implicit none

  vectype, dimension(*), intent(inout) :: vec
  real :: tstart, tend
  vectype, allocatable :: temp(:)

  if(hyperc.eq.0) return

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

  call createvec(temp,2)

  ! smooth vorticity
  call newvar1(mass_matrix_lhs_dc,vor,vec,u_i,vecsize_vel, &
       gs_matrix_rhs_dc,NV_DCBOUND,vor)

  temp = 0.
  call copyvec(vor,1,1,temp,2,2)
  call newvar(s5matrix_sm,temp,d5matrix_sm,temp,NV_SVBOUND,temp)
  call copyvec(temp,2,2,vec,u_i,vecsize_vel)
  
  ! smooth compression
  if(numvar.ge.3) then
     if(com_bc.eq.0) then
        call newvar1(mass_matrix_lhs   ,com,vec,chi_i,vecsize_vel, &
             lp_matrix_rhs,   NV_NOBOUND,com)
     else
        call newvar1(mass_matrix_lhs_dc,com,vec,chi_i,vecsize_vel, &
             lp_matrix_rhs_dc,NV_DCBOUND,com)
     endif

     temp = 0.
     call copyvec(com,1,1,temp,2,2)
     call newvar(s7matrix_sm,temp,d7matrix_sm,temp,NV_SCBOUND,temp)
     call copyvec(temp,2,2,vec,chi_i,vecsize_vel)
  endif

  call deletevec(temp)
  
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
subroutine smooth_fields(vec)
  use basic
  use arrays
  use newvar_mod
  use diagnostics
  use sparse

  implicit none

  vectype, dimension(*), intent(inout) :: vec
  vectype, allocatable :: temp(:)
  real :: tstart, tend

  if(jadv.eq.0 .or. hyper.eq.0.) return

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

  ! smooth current density
  call newvar1(mass_matrix_lhs_dc,jphi,vec,psi_i,vecsize_phi, &
       gs_matrix_rhs_dc,NV_DCBOUND,jphi)

  call createvec(temp,2)
  temp = 0.
  call copyvec(jphi,1,1,temp,2,2)
  call newvar(s10matrix_sm,temp,d10matrix_sm,temp,NV_SJBOUND,temp)
  call copyvec(temp,2,2,vec,psi_i,vecsize_phi)
  call deletevec(temp)

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
subroutine derived_quantities(vec, bfv)
  use basic
  use arrays
  use newvar_mod
  use diagnostics
  use sparse

  implicit none

  real :: tstart, tend
  integer :: ndofs
  vectype, dimension(*), intent(in) :: vec
  vectype, dimension(*), intent(inout) :: bfv

  if(myrank.eq.0 .and. iprint.ge.1) print *, "Calculating derived fields..."

  ! Find lcfs
  ! ~~~~~~~~~
  if(eqsubtract.eq.1) then
     if(ntime.eq.ntime0) &
          call lcfs(field0,psi_g,num_fields)
  else
     call lcfs(field,psi_g,num_fields)
  endif


  ! Define auxiliary fields
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.1) print *, "Defining auxiliary fields."
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

  if(myrank.eq.0 .and. iprint.ge.1) print *, "-Transport coefficients"
  call define_transport_coefficients

  !   toroidal current
  if(myrank.eq.0 .and. iprint.ge.1) print *, "-Toroidal current"
  call newvar1(mass_matrix_lhs_dc,jphi,vec,psi_g,num_fields, &
       gs_matrix_rhs_dc,NV_DCBOUND,jphi)

  if(hyperc.ne.0.) then
     !   vorticity
     if(myrank.eq.0 .and. iprint.ge.1) print *, "-Vorticity"
     call newvar1(mass_matrix_lhs_dc,vor,vec,u_g,num_fields, &
          gs_matrix_rhs_dc,NV_DCBOUND,vor)

     !   compression
     if(numvar.ge.3) then
        if(myrank.eq.0 .and. iprint.ge.1) print *, "-Compression"
        if(com_bc.eq.1) then
           call newvar1(mass_matrix_lhs_dc,com,vec,chi_g,num_fields, &
                lp_matrix_rhs_dc,NV_DCBOUND,com)
        else
           call newvar1(mass_matrix_lhs,   com,vec,chi_g,num_fields, &
                lp_matrix_rhs,   NV_NOBOUND,com)
        endif
     else
        com = 0.
     endif
  endif

  ! vector potential stream function
  if((i3d.eq.1 .or. ifout.eq.1) .and. numvar.ge.2) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, "-f"
     call numdofs(1, ndofs)
     bfi = bfv(1:ndofs)
     call newvar1(bf_matrix_lhs_nm,bfv,vec,bz_g,num_fields, &
          mass_matrix_rhs_nm,NV_NMBOUND,bfi)
  endif

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_aux = t_aux + tend - tstart
  endif

  ! calculate scalars
  ! ~~~~~~~~~~~~~~~~~
  if(icalc_scalars.eq.1) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Calculating scalars"
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
  
  integer :: numnodes, j, ibegin, iendplusone, ndofs
  integer, allocatable :: itemp(:)

  call numdofs(num_fields,ndofs)
  allocate(itemp(ndofs))
  itemp = 1

  ! adjust toroidal field and boundary condition to conserve toroidal flux
  if(numvar.ge.2 .and. iconstflux.eq.1) then
     gbound = (tflux0-tflux)/area
     call numnod(numnodes)
     do j=1,numnodes
        call entdofs(num_fields, j, 0, ibegin, iendplusone)
        if(itemp(ibegin+(bz_g-1)*6) .eq. 0) cycle
        field(ibegin+(bz_g-1)*6) = field(ibegin+(bz_g-1)*6) + gbound
        itemp(ibegin+(bz_g-1)*6) = 0
     enddo
     if(myrank.eq.0) then
        print *, "Correction to toroidal flux: ", gbound*area
     end if
  endif

  deallocate(itemp)

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

  integer :: i, ndofs, numnodes, ibegin, iendplusone
  integer, allocatable ::  itemp(:)

  call numnod(numnodes)
  call numdofs(num_fields, ndofs)
  allocate(itemp(ndofs))
  itemp = 1
  do i=1,numnodes
     call entdofs(num_fields, i, 0, ibegin, iendplusone)
     if(itemp(ibegin) .eq. 1) then
        call assign_local_pointers(i)
        
        psi0_l = -psi0_l
        psi1_l = -psi1_l
        psis_l = -psis_l
          u0_l = -  u0_l
          u1_l = -  u1_l
          us_l = -  us_l
         bz0_l = - bz0_l
         bz1_l = - bz1_l
         bzs_l = - bzs_l
         vz0_l = - vz0_l
         vz1_l = - vz1_l
         vzs_l = - vzs_l
        
        itemp(ibegin) = 0
     endif
  enddo
  deallocate(itemp)

  psimin = -psimin
  psilim = -psilim
  psibound = -psibound

end subroutine flip_handedness


logical function inside_lcfs(psi, x, z, exclude_pf)
  use basic

  vectype, intent(in), dimension(6) :: psi
  real, intent(in) :: x, z
  logical :: exclude_pf
  real :: dpsii

  dpsii = psibound - psimin

  if((real(psi(1)) - psimin)/dpsii .gt. 1.) then
     inside_lcfs = .false.
     return
  endif

  if(exclude_pf) then
     if((real(psi(2))*(x-xmag) + real(psi(3))*(z-zmag))*dpsii .lt. 0.) then
        inside_lcfs = .false.
        return
     endif
  endif

  inside_lcfs = .true.
  return
end function inside_lcfs
