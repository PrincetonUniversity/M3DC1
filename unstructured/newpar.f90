Program Reducedquintic

!   Ref:  [1] Strang and Fix, An Analysis of the Finite Element Method, page 83
!         [2] G.R. Cowper, et al, AIAA Journal, Vol 7, no. 10 page 19
!         [3] JCP Vol 147 p 318-336 (1998)

  use p_data
  use basic
  use arrays
  use newvar_mod
  use sparse
  use hdf5_output
  use diagnostics
  use vacuum_interface
  use boundary_conditions
  use time_step

  implicit none

#include "finclude/petsc.h"

  integer :: ier
  real :: tstart, tend

  ! Start up message passing, SuperLU process grid
  call MPI_Init(ier)
  if (ier /= 0) then
     print *,'Error initializing MPI:',ier
     call safestop(1)
  endif
  call PetscInitialize(PETSC_NULL_CHARACTER, ier)
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

  ! load mesh
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Loading mesh'
  call load_mesh

  ! read input file
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Reading input'
  call input

  ! allocate arrays
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Allocating arrays'
  call space(1)

  ! initialize variables
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Initializing variables'
  call init

  ! load resistive wall response matrix
  if(eta_wall.ne.0. .and. itaylor.ne.10 .and. itaylor.ne.12) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, ' Loading VACUUM data'
     call load_vacuum_data(ntor, ier)
     if(ier.ne.0) call safestop(7)
  end if

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

     ! initialize t(n-1) values
     phiold_vec = phi_vec

     if(isplitstep.eq.1) then
        velold_vec = vel_vec
        if(idens.eq.1) denold_vec = den_vec
        if(ipres.eq.1) presold_vec = pres_vec
     endif

  case(1)
     ! Read restart file(s)

     if(myrank.eq.1 .and. iprint.ge.1) print *, ' Reading restart file(s)'
     if(iglobalin.eq.1) then
        call rdrestartglobal
     else
        call rdrestart
     endif
  end select                     !  end of the branch on restart/no restart

  ntime0 = ntime

  ! initialize hdf5
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Initializing HDF5'
  call hdf5_initialize(ntime.gt.0,ier)
  if(ier.lt.0) then 
     print *, "Error initializing HDF5"
     call safestop(5)
  end if

  if(itimer.eq.1) call reset_timings

  ! output simulation parameters and equilibrium
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(irestart.eq.0  .or. iadapt.gt.0) then
     if(myrank.eq.0 .and. iprint.ge.1) &
          print *, " Writing simulation parameters"
     call hdf5_write_parameters(ier)

     if(eqsubtract.eq.1) then
        call derived_quantities(0)
     end if

     ! Output the equilibrium
     if(myrank.eq.0 .and. iprint.ge.2) print *, ' Writing equilibrium'
     call hdf5_write_time_slice(1,ier)
  end if


  ! Calculate all quantities derived from basic fields
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  call derived_quantities(1)


  ! Adapt the mesh
  ! ~~~~~~~~~~~~~~
#ifdef USESCOREC
  select case(iadapt)
  case(1)
     print *, 'adapting mesh...'
     call hessianadapt(resistivity_field%vec%data,1, 0, ntime, &
          adapt_factor, adapt_hmin, adapt_hmax)
     print *, 'done adapting.'
     call space(0)
     call tridef
     
!     call updatenormalcurvature
     call write_normlcurv
     
  case(2)
    call create_field(temporary_field)
    temporary_field = psi_field(0)

    print *, 'adapting mesh...'
    call adapt(temporary_field%vec%data,psimin,psibound)
    print *, 'done adapting.'

    call destroy_field(temporary_field)

    call space(0)
    call tridef

!    call updatenormalcurvature
    call write_normlcurv
  end select
#endif

  if(irestart.eq.0  .or. iadapt.gt.0) then
     tflux0 = tflux
     totcur0 = totcur

     if(itaylor.ne.10 .and. itaylor.ne.11 .and. &
          itaylor.ne.12 .and. itaylor.ne.13) then
        external_psi_field = psi_field(1)
        external_bz_field = bz_field(1)
        external_bf_field = bf_field(1)
     endif
  endif

  ! output initial conditions
  call output

  ! if there are no timesteps to calculate, then skip time loop
  if(ntimemax.le.ntime) go to 101

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Initializing timestep'
  call initialize_timestep

  ! main time loop
  ! ~~~~~~~~~~~~~~
  do ntime=ntime+1,ntimemax

     if(myrank.eq.0) print *, 'TIME STEP: ',ntime

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
     endif

     ! feedback control on toroidal current
     if(myrank.eq.0 .and. iprint.ge.1) print *, " Applying feedback"
     call control(totcur, vloop,       i_control, dt)
     call control(totden, pellet_rate, n_control, dt)

     ! Write output
     if(myrank.eq.0 .and. iprint.ge.1) print *, " Writing output."
     call output
  enddo ! ntime

  if(myrank.eq.0 .and. iprint.ge.1) print *, "Done time loop."

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

  n_control%p = n_control_p
  n_control%i = n_control_i
  n_control%d = n_control_d
  n_control%target_val = n_target
end subroutine init


!============================================
! print_info
! ~~~~~~~~~~
! print basic info about simulation options
!============================================
subroutine print_info
  use basic

  implicit none

  integer :: ndofs, numelms, numnodes, j

  ! velocity form
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
#ifdef USESCOREC
     if(myrank.eq.0) then
        call numglobaldofs(1, ndofs)
        print *, 'global dofs = ', ndofs
     end if

     numelms = local_elements()
     numnodes = owned_nodes()
     call numprocdofs(1, j)
     call numdofs(1, ndofs)
     print *, 'proc, owned dofs, needed dofs',myrank, j, ndofs
     print *, 'proc, numnodes, numfaces', myrank, numnodes,numelms
#endif
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
  use hdf5_output
  use vacuum_interface
  use time_step

  implicit none
#include "finclude/petsc.h"
      
  integer, intent(in) :: iarg
  integer :: ier
      
  call finalize_timestep

  ! unload resistive wall response matrix
  call unload_vacuum_data

  ! close hdf5 file
  if(myrank.eq.0 .and. iprint.ge.2) print *, "  finalizing hdf5..."
  call hdf5_finalize(ier)
  if(ier.ne.0) print *, 'Error finalizing HDF5:',ier

  if(myrank.eq.0 .and. iprint.ge.2) print *, "  deleting matrices..."
  call delete_matrices

  if(myrank.eq.0 .and. iprint.ge.2) print *, "  unloading mesh..."
  call unload_mesh

  if(myrank.eq.0 .and. iprint.ge.2) print *, "  finalizing PETSC..."
  call PetscFinalize(ier)
  if(ier.ne.0) print *, 'Error finalizing PETSC:', ier

  if(myrank.eq.0 .and. iprint.ge.2) print *, "  finalizing MPI..."
  call MPI_Finalize(ier)
  if (ier.ne.0) print *, 'Error finalizing MPI:', ier
  
  print *, "Stopped at", iarg
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
  if(myrank.eq.0 .and. iprint.ge.2) print *, "  writing scalars"
  call hdf5_write_scalars(ier)
  if(mod(ntime-ntime0,ntimepr).eq.0) then
     if(myrank.eq.0 .and. iprint.ge.2) print *, "  writing timeslice"
     call hdf5_write_time_slice(0,ier)

     if(iwrite_restart.eq.1) then
        if(myrank.eq.0 .and. iprint.ge.2) print *, "  writing restart files"
        if(iglobalout.eq.1) then
           call wrrestartglobal
        else
           call wrrestart
        endif
     endif
  endif
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_output_hdf5 = t_output_hdf5 + tend - tstart
  end if
  
  if(itimer.eq.1) then
     if(myrank.eq.0) call second(tstart)
     if(myrank.eq.0 .and. iprint.ge.2) print *, "  writing timings"
     call hdf5_write_timings(ier)
     call reset_timings
  end if
  
  ! flush hdf5 data to disk
  if(myrank.eq.0 .and. iprint.ge.2) print *, "  flushing data to file"
  call hdf5_flush(ier)
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
  type(field_type) :: j_new
  real :: tstart, tend

  if(jadv.eq.0 .or. hyper.eq.0.) return

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

  ! smooth current density
  call solve_newvar1(mass_mat_lhs_dc, jphi_field, gs_mat_rhs_dc, psiin)

  call create_vector(temp_vec,2)
  call associate_field(j_new, temp_vec, 2)

  j_new = jphi_field
  call solve_newvar_axby(s10_mat,temp_vec,d10_mat,temp_vec)
  psiin = j_new

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

  implicit none

  integer, intent(in) :: ilin    ! 0 for equilibrium fields, 1 for perturbed

  real :: tstart, tend

  vectype :: temp

  ! Find lcfs
  ! ~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.2) print *, "  finding lcfs"
  if(eqsubtract.eq.1) then
     if(ntime.eq.ntime0) call lcfs(psi_field(0))
  else
     call lcfs(psi_field(1))
  endif

  ! Define auxiliary fields
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

  if(myrank.eq.0 .and. iprint.ge.2) &
       print *, "  transport coefficients"
  call define_transport_coefficients

  !   toroidal current
  if(myrank.eq.0 .and. iprint.ge.2) print *, "  toroidal current"
  call solve_newvar1(mass_mat_lhs_dc,jphi_field,gs_mat_rhs_dc,psi_field(ilin))

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
  if((i3d.eq.1 .or. ifout.eq.1) .and. numvar.ge.2) then
     if(myrank.eq.0 .and. iprint.ge.2) print *, "  f"
     if(ilin.eq.0) then 
        temp = bzero*rzero
        call add(bz_field(ilin),-temp)
     endif
     call solve_newvar1(bf_mat_lhs,bf_field(ilin),mass_mat_rhs_bf, &
          bz_field(ilin), bf_field(ilin))
     if(ilin.eq.0) call add(bz_field(ilin), temp)
  endif

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_aux = t_aux + tend - tstart
  endif

  ! calculate scalars
  ! ~~~~~~~~~~~~~~~~~
  if(icalc_scalars.eq.1) then
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
  end do

  psimin = -psimin
  psilim = -psilim
  psibound = -psibound

end subroutine flip_handedness


!======================================================================
! inside_lcfs
! ~~~~~~~~~~~
! Flips coordinate system handedness by flipping sign of
! psi, u, vz, and bz
!======================================================================
logical function inside_lcfs(psi, x, z, exclude_pf)
  use basic

  vectype, intent(in), dimension(dofs_per_node) :: psi
  real, intent(in) :: x, z 
  logical :: exclude_pf    ! if true, count private flux region as outside
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
! evaluate
! ~~~~~~~~
! calculates the value ans of field dum at global coordinates
! (x,z).  itri is the element containing (x,z).  (If this
! element does not reside on this process, itri=-1).
!============================================================
subroutine evaluate(x,z,ans,ans2,fin,itri)
  
  use p_data
  use mesh_mod
  use basic
  use nintegrate_mod
  use field

  implicit none

  include 'mpif.h'

  integer, intent(inout) :: itri
  real, intent(in) :: x, z
  type(field_type), intent(in) :: fin

  real, intent(out) :: ans, ans2

  type(element_data) :: d
  integer :: p, nodeids(nodes_per_element), ier
  real :: x1, z1
  vectype, dimension(20) :: avector
  real :: ri, si, eta
  real :: term1, term2
  real, dimension(2) :: temp1, temp2
  integer :: hasval, tothasval

  ! evaluate the solution to get the value [ans] at one point (x,z)

  ! first find out what triangle x,z is in.  whattri
  ! returns itri, x1, and z1 with x1 and z1 being
  ! the coordinates of the first node/vertex

  if(itri.eq.0) then
     call whattri(x,z,itri,x1,z1)
  else
     call get_element_nodes(itri,nodeids)
     call get_node_pos(nodeids(1), x1, z1)
  endif

  ans = 0.
  ans2 = 0.


  ! if this process contains the point, evaluate the field at that point.
  if(itri.gt.0) then

     call get_element_data(itri, d)

     ! calculate local coordinates
     call global_to_local(d, x, z, si, eta)

     ! calculate the inverse radius
     if(itor.eq.1) then
        ri = 1./x
     else
        ri = 1.
     endif

     ! calculate the value of the function
     call calcavector(itri, fin, avector)
     
     do p=1,20
     
        term1 = si**mi(p)*eta**ni(p)
        term2 = 0.
        
        if(mi(p).ge.1) then
           if(itor.eq.1) then
              term2 = term2 - 2.*d%co*(mi(p)*si**(mi(p)-1) * eta**ni(p))*ri
           endif
           
           if(mi(p).ge.2) then
              term2 = term2 + si**(mi(p)-2)*(mi(p)-1)*mi(p) * eta**ni(p)
           endif
        endif
     
        if(ni(p).ge.1) then
           if(itor.eq.1) then
              term2 = term2 + 2.*d%sn*(si**mi(p) * eta**(ni(p)-1)*ni(p))*ri
           endif
           
           if(ni(p).ge.2) then
              term2 = term2 + si**mi(p) * eta**(ni(p)-2)*(ni(p)-1)*ni(p)
           endif
        endif
     
        ans = ans + avector(p)*term1
        ans2 = ans2 + avector(p)*term2
        hasval = 1
     enddo
  else
     hasval = 0
  endif
     

  ! Distribute the result if necessary
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(maxrank.gt.1) then
     ! Determine number of processes whose domains contain this point
     call mpi_allreduce(hasval, tothasval, 1, MPI_INTEGER, MPI_SUM, &
          MPI_COMM_WORLD, ier)

     ! Find the average value at this point over all processes containing
     ! the point.  (Each value should be identical.)
     temp1(1) = ans
     temp1(2) = ans2
     call mpi_allreduce(temp1, temp2, 2, MPI_DOUBLE_PRECISION, MPI_SUM, &
          MPI_COMM_WORLD, ier)
     ans = temp2(1)/tothasval
     ans2 = temp2(2)/tothasval
  endif

end subroutine evaluate

  !============================================================
  ! tridef
  ! ~~~~~~
  ! populates the *tri arrays
  !============================================================
  subroutine tridef

    use math
    use arrays
    use mesh_mod
    use time_step

    implicit none
  
    type(element_data) :: d
    integer :: itri, i, j, k, ii, jj, numelms, numnodes, ndofs
    integer, dimension(num_fields, dofs_per_node) :: ind
    real :: ti(20,20),rot(18,18),newrot(18,18), sum, theta
    logical, allocatable :: is_set(:)
    real :: norm(2), curv, x, z
    integer :: inode(nodes_per_element)
    logical :: is_boundary
    integer :: izone, izonedim
  
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
       call tmatrix(ti,20,d%a,d%b,d%c)
       
       ! calculate the rotation matrix rot
       theta = atan2(d%sn,d%co)
       call rotation(rot,18,theta)
       
       newrot = 0.
       call get_element_nodes(itri, inode)
       do i=1, nodes_per_element
          call boundary_node(inode(i), &
               is_boundary, izone, izonedim, norm, curv, x, z)
          k = (i-1)*dofs_per_node + 1
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

             if(abs(norm(1)**2 + norm(2)**2 - 1.) > 0.01) &
                  print *, 'warning: bad normals', norm(1), norm(2)
          else
             do j=1, dofs_per_node
                newrot(k+j-1,k+j-1) = 1.
             end do
          end if
       end do

       ! form the matrix g by multiplying ti and rot
       do k=1,20
          do j=1,18
             sum = 0.
             do ii = 1,18
                sum = sum + ti(k,ii)*rot(ii,j)
             enddo
             gtri_old(k,j,itri) = sum

             sum = 0.
             do ii = 1,18
                do jj=1, 18
                   sum = sum + newrot(j,jj)*ti(k,ii)*rot(ii,jj)
                end do
             enddo

             gtri(k,j,itri) = sum
          enddo
       enddo
    end do
  end subroutine tridef

! write_normlcurv
subroutine write_normlcurv

  use basic
  use mesh_mod
  
  implicit none
  
  integer :: numnodes, i, j, inode(4), nbound, numelms, itri
  integer :: i1, i2
  real :: dx1, dx2, dz1, dz2, norm1(2), norm2(2), l, dl1, dl2, dl
  real :: vx1, vx2, vz1, vz2, ax, az
  integer, allocatable :: id(:), adjacent(:,:), nn(:)
  real, allocatable :: x(:), z(:), norm(:,:), curv(:), curv_new(:)

  if(maxrank.gt.1) then
     if(myrank.eq.0) print *, 'write_normlcurv can only be called in serial.'
     return
  end if
     
  numnodes = local_nodes()
  numelms = local_elements()

  ! calculate number of boundary nodes
  nbound = 0
  do i=1, numnodes
     if(is_boundary_node(i)) nbound = nbound + 1
  end do
  
  ! allocate memory
  allocate(x(nbound), z(nbound), adjacent(2,nbound), norm(2,nbound), &
       id(nbound), nn(numnodes), curv(nbound), curv_new(nbound))

  nbound = 0
  nn = 0
  do i=1, numnodes
     if(.not.is_boundary_node(i)) cycle
     nbound = nbound + 1
     
     call get_node_pos(i,x(nbound),z(nbound))
     id(nbound) = i
     nn(i) = nbound
  end do
  
  ! determine adjacent nodes
  adjacent = 0
  do itri=1, numelms

     call get_element_nodes(itri,inode)

     do i=1,nodes_per_element
        j = mod(i,nodes_per_element) + 1
     
        if((nn(inode(i)).eq.0).or.(nn(inode(j)).eq.0)) cycle
        
        if(adjacent(1,nn(inode(i))).eq.0) then
           adjacent(1,nn(inode(i))) = nn(inode(j))
        else if(adjacent(2,nn(inode(i))).eq.0) then
           adjacent(2,nn(inode(i))) = nn(inode(j))
        else
           print *, "Error in write_normlcurv", 1
        endif
        if(adjacent(1,nn(inode(j))).eq.0) then
           adjacent(1,nn(inode(j))) = nn(inode(i))
        else if(adjacent(2,nn(inode(j))).eq.0) then
           adjacent(2,nn(inode(j))) = nn(inode(i))
        else
           print *, "Error in write_normlcurv", 1
        endif
     end do
  end do

  ! calculate normals
  do i=1, nbound
     i1 = adjacent(1,i)
     i2 = adjacent(2,i)

     if(i1.eq.0 .or. i2.eq.0) then
        print *, 'Error in write_normlcurv', 3
        cycle
     endif

     dx1 = x(i) - x(i1)
     dx2 = x(i2) - x(i)
     dz1 = z(i) - z(i1)
     dz2 = z(i2) - z(i)

     dl1 = sqrt(dx1**2 + dz1**2)
     dl2 = sqrt(dx2**2 + dz2**2)
     norm1(1) =  dz1/dl1
     norm1(2) = -dx1/dl1
     norm2(1) =  dz2/dl2
     norm2(2) = -dx2/dl2

     ! perform weigted average of adjacent edge normals
     norm(:,i) = (norm1/dl1 + norm2/dl2) / (1./dl1 + 1./dl2)

     ! normalize normal
     l = sqrt(norm(1,i)**2 + norm(2,i)**2)
     norm(:,i) = norm(:,i)/l

     ! calculate curvature
     vx1 = dx1/dl1
     vx2 = dx2/dl2
     vz1 = dz1/dl1
     vz2 = dz2/dl2

     dl = .5*(dl1 + dl2)

     ax = (vx2 - vx1) / dl
     az = (vz2 - vz1) / dl

     curv(i) = sqrt(ax**2 + az**2)

     ! make sure normal is pointing outward
     if(norm(1,i)*(x(i)-xmag) + norm(2,i)*(z(i)-zmag) .lt. 0) then
        norm(:,i) = -norm(:,i)
     endif
  end do

  ! write normlcurv
  open(unit=23, file='normlcurv_new', status='unknown')

  do j=1,numnodes
     i = nn(j)
     if(i.eq.0) cycle
     write(23,'(I5,5F10.6)') j, x(i), z(i), norm(1,i), norm(2,i), curv(i)
  end do

  close(23)
  
  ! free memory
  deallocate(x,z,adjacent,id,nn,norm,curv,curv_new)

end subroutine write_normlcurv


!============================================================
! space
! ~~~~~
! allocates space for big arrays
!
! ifirstcall = 1 if this is the first call and 0 otherwise
!============================================================
subroutine space(ifirstcall)

  use p_data
  use mesh_mod
  use basic
  use arrays
  use sparse
  use time_step

  implicit none

  integer, intent(in) :: ifirstcall

  integer :: maxdofs1, maxdofs2, maxdofs3, maxdofsn
  integer :: numelms

  if(myrank.eq.1 .and. iprint.ge.1) print *, " Entering space..."

  if(isplitstep.eq.1) then
     vecsize_phi = numvar
     vecsize_vel = numvar
     vecsize_n = 1
     vecsize_p = 1
  else
     vecsize_phi  = numvar*2 + idens + ipres
  endif

! add electrostatic potential equation
#ifdef USECOMPLEX
  if(jadv.eq.0) vecsize_phi = vecsize_phi + 1
#endif

  if(isplitstep.eq.0) then
     vecsize_vel  = vecsize_phi
     vecsize_n    = vecsize_phi
     vecsize_p    = vecsize_phi
  endif

!
!.....create numberings
#ifdef USESCOREC
  if(ifirstcall .eq. 1) then
     call createdofnumbering(numvar1_numbering, iper, jper, &
          6, 0, 0, 0, maxdofs1)
     call createdofnumbering(numvar2_numbering, iper, jper, &
          12, 0, 0, 0, maxdofs2)
     call createdofnumbering(numvar3_numbering, iper, jper, &
          18, 0, 0, 0, maxdofs3)
     if(num_fields.gt.3) then
        call createdofnumbering(num_fields, iper, jper, &
             num_fields*6, 0, 0, 0, maxdofsn)
     endif
     if(vecsize_phi.gt.3 .and. vecsize_phi.ne.num_fields) then
        call createdofnumbering(vecsize_phi, iper, jper, &
             vecsize_phi*6, 0, 0, 0, maxdofsn)
     endif
     if(vecsize_vel.gt.3 .and. vecsize_vel.ne.vecsize_phi) then
        call createdofnumbering(vecsize_vel, iper, jper, &
             vecsize_vel*6, 0, 0, 0, maxdofsn)
     endif

#ifdef USERW
     ! set resistive wall bcs for relevant numberings
     if(eta_wall.ne.0) then
        call setresistivewallbcstate(1,3)
        if(vecsize_phi .ne. 3) then
           call setresistivewallbcstate(1,vecsize_phi)
        endif
     endif
#endif

  endif ! on firstcall
#endif
  
  numelms = local_elements()

! arrays defined at all vertices
! createvec will delete the arrays if they have already been allocated
  if(ifirstcall.eq.1) then
     if(myrank.eq.0 .and. iprint.eq.1) print *, 'Allocating...'

     ! Physical Variables
     call create_vector(field_vec , num_fields)
     call create_vector(field0_vec, num_fields)

     ! Auxiliary Variables
     call create_field(jphi_field)
     call create_field(vor_field)
     call create_field(com_field)
     call create_field(resistivity_field)
     call create_field(kappa_field)
     call create_field(visc_field)
     call create_field(visc_c_field)
     call create_field(sigma_field)
     call create_field(tempvar_field)
     call create_field(bf_field(0))
     call create_field(bf_field(1))
     if(ibootstrap.gt.0) call create_field(visc_e_field)

     ! External fields
     call create_vector(external_field, 3)
     call associate_field(external_psi_field, external_field, 1)
     call associate_field(external_bz_field,  external_field, 2)
     call associate_field(external_bf_field,  external_field, 3)

     ! Arrays for implicit time advance
     call create_vector(phi_vec,      vecsize_phi)
     call create_vector(phiold_vec,   vecsize_phi)
     call create_vector(phip_vec,     vecsize_phi)
     call create_vector(q4_vec,       vecsize_phi)
     
     call create_vector(b1_phi, vecsize_phi)
     call create_vector(b2_phi, vecsize_phi)

     if(isplitstep.eq.1) then
        call create_vector(vel_vec,      vecsize_vel)
        call create_vector(velold_vec,   vecsize_vel)
        call create_vector(veln_vec,     vecsize_vel)
        call create_vector(veloldn_vec,  vecsize_vel)
        call create_vector(r4_vec,       vecsize_vel)

        if(ipres.eq.1) then
           call create_vector(pres_vec,    vecsize_p)
           call create_vector(presold_vec, vecsize_p)
           call create_vector(qp4_vec,     vecsize_p)
        endif
        
        call create_vector(den_vec,    vecsize_n)        
        call create_vector(denold_vec, vecsize_n)
        call create_vector(qn4_vec,    vecsize_n)
        
        call create_vector(b1_vel, vecsize_vel)
        call create_vector(b2_vel, vecsize_vel)
     endif
  endif

  ! arrays associated with the triangles
  if(ifirstcall.eq.0) then
     if(myrank.eq.0 .and. iprint.eq.1) print *, ' deallocating...'
     deallocate(gtri,gtri_old)
  endif
  
  if(myrank.eq.0 .and. iprint.eq.1) print *, ' Allocating tri...'
  allocate(gtri(20,18,numelms),gtri_old(20,18,numelms)) 

  if(myrank.eq.0 .and. iprint.eq.1) print *, ' associating...'
  call associate_field(u_field(1),   field_vec, u_g)
  call associate_field(vz_field(1),  field_vec, vz_g)
  call associate_field(chi_field(1), field_vec, chi_g)
  call associate_field(psi_field(1), field_vec, psi_g)
  call associate_field(bz_field(1),  field_vec, bz_g)
  call associate_field(pe_field(1),  field_vec, pe_g)
  call associate_field(den_field(1), field_vec, den_g)
  call associate_field(p_field(1),   field_vec, p_g)
  call associate_field(u_field(0),   field0_vec, u_g)
  call associate_field(vz_field(0),  field0_vec, vz_g)
  call associate_field(chi_field(0), field0_vec, chi_g)
  call associate_field(psi_field(0), field0_vec, psi_g)
  call associate_field(bz_field(0),  field0_vec, bz_g)
  call associate_field(pe_field(0),  field0_vec, pe_g)
  call associate_field(den_field(0), field0_vec, den_g)
  call associate_field(p_field(0),   field0_vec, p_g)


  ! assign pointers to proper vectors
  if(myrank.eq.0 .and. iprint.eq.1) print *, ' assinging...'
  call assign_variables

  if(myrank.eq.0 .and. iprint.ge.1) print *, " Exiting space."

  return
end subroutine space


#ifdef USESCOREC
subroutine resizevec(vec, ivecsize)
  use arrays
  implicit none

  integer :: ivecsize
  double precision :: vec

  call arrayresizevec(vec, ivecsize)

  return
end subroutine resizevec

subroutine arrayresizevec(vec, ivecsize)
  use arrays
  use time_step

  implicit none
  integer :: ivecsize, i
  double precision :: vec

  print *, "In arrayresizevec!", ivecsize

  call checkppplveccreated(vec, i)
  if(i .eq. 0) then
     call printfpointer(vec)
     write(*,*) 'trying to resize a vector that has not been created'
     call safestop(8844)
  endif


  call checksameppplvec(field_vec%data, vec, i)
  if(i .eq. 1) then
     print *, "field"
     if(allocated(field_vec%data)) deallocate(field_vec%data, STAT=i)
     allocate(field_vec%data(ivecsize))
     field_vec%data = 0.
     call updateids(vec, field_vec%data)
     print *, 'done'
     return
  endif

  call checksameppplvec(field0_vec%data, vec, i)
  if(i .eq. 1) then
     print *, "field0"
     if(allocated(field0_vec%data)) deallocate(field0_vec%data, STAT=i)
     allocate(field0_vec%data(ivecsize))
     field0_vec%data = 0.
     call updateids(vec, field0_vec%data)
     return
  endif

  call checksameppplvec(jphi_field%vec%data, vec, i)
  if(i .eq. 1) then
     print *, "jphi"
     if(allocated(jphi_field%vec%data)) deallocate(jphi_field%vec%data, STAT=i)
     allocate(jphi_field%vec%data(ivecsize))
     jphi_field%vec%data = 0.
     call updateids(vec, jphi_field%vec%data)
     return
  endif
  
  call checksameppplvec(vor_field%vec%data, vec, i)
  if(i .eq. 1) then
     print *, "vor"
     if(allocated(vor_field%vec%data)) deallocate(vor_field%vec%data, STAT=i)
     allocate(vor_field%vec%data(ivecsize))
     vor_field%vec%data = 0.
     call updateids(vec, vor_field%vec%data)
     return
  endif
  
  call checksameppplvec(com_field%vec%data, vec, i)
  if(i .eq. 1) then
     print *, "com"
     if(allocated(com_field%vec%data)) deallocate(com_field%vec%data, STAT=i)
     allocate(com_field%vec%data(ivecsize))
     com_field%vec%data = 0.
     call updateids(vec, com_field%vec%data)
     return
  endif
  
  call checksameppplvec(resistivity_field%vec%data, vec, i)
  if(i .eq. 1) then
     print *, "resistivity_field%vec%data"
     if(allocated(resistivity_field%vec%data)) deallocate(resistivity_field%vec%data, STAT=i)
     allocate(resistivity_field%vec%data(ivecsize))
     resistivity_field%vec%data = 0.
     call updateids(vec, resistivity_field%vec%data)
     return
  endif
  
  call checksameppplvec(tempvar_field%vec%data, vec, i)
  if(i .eq. 1) then
     print *, "tempvar"
     if(allocated(tempvar_field%vec%data)) deallocate(tempvar_field%vec%data, STAT=i)
     allocate(tempvar_field%vec%data(ivecsize))
     tempvar_field%vec%data = 0.
     call updateids(vec, tempvar_field%vec%data)
     return
  endif
  
  call checksameppplvec(kappa_field%vec%data, vec, i)
  if(i .eq. 1) then
     print *, "kappa_field%vec%data"
     if(allocated(kappa_field%vec%data)) deallocate(kappa_field%vec%data, STAT=i)
     allocate(kappa_field%vec%data(ivecsize))
     kappa_field%vec%data = 0.
     call updateids(vec, kappa_field%vec%data)
     return
  endif

  call checksameppplvec(sigma_field%vec%data, vec, i)
  if(i .eq. 1) then
     print *, "sigma_field%vec%data"
     if(allocated(sigma_field%vec%data)) deallocate(sigma_field%vec%data, STAT=i)
     allocate(sigma_field%vec%data(ivecsize))
     sigma_field%vec%data = 0.
     call updateids(vec, sigma_field%vec%data)
     return
  endif
        
  call checksameppplvec(visc_field%vec%data, vec, i)
  if(i .eq. 1) then
     print *, "visc_field%vec%data"
     if(allocated(visc_field%vec%data)) deallocate(visc_field%vec%data, STAT=i)
     allocate(visc_field%vec%data(ivecsize))
     visc_field%vec%data = 0.
     call updateids(vec, visc_field%vec%data)
     return
  endif
  
  call checksameppplvec(visc_c_field%vec%data, vec, i)
  if(i .eq. 1) then
     print *, "visc_c"
     if(allocated(visc_c_field%vec%data)) deallocate(visc_c_field%vec%data, STAT=i)
     allocate(visc_c_field%vec%data(ivecsize))
     visc_c_field%vec%data = 0.
     call updateids(vec, visc_c_field%vec%data)
     return
  endif

  call checksameppplvec(bf_field(1)%vec%data, vec, i)
  if(i .eq. 1) then
     print *, "bf"
     if(allocated(bf_field(1)%vec%data)) deallocate(bf_field(1)%vec%data, STAT=i)
     allocate(bf_field(1)%vec%data(ivecsize))
     bf_field(1)%vec%data = 0.
     call updateids(vec, bf_field(1)%vec%data)
     return
  endif

  call checksameppplvec(bf_field(0)%vec%data, vec, i)
  if(i .eq. 1) then
     print *, "bf"
     if(allocated(bf_field(0)%vec%data)) deallocate(bf_field(0)%vec%data, STAT=i)
     allocate(bf_field(0)%vec%data(ivecsize))
     bf_field(0)%vec%data = 0.
     call updateids(vec, bf_field(0)%vec%data)
     return
  endif


  call checksameppplvec(phi_vec%data, vec, i)
  if(i .eq. 1) then
     print *, "phi"
     if(allocated(phi_vec%data)) deallocate(phi_vec%data, STAT=i)
     allocate(phi_vec%data(ivecsize))
     phi_vec%data = 0.
     call updateids(vec, phi_vec%data)
     return
  endif
     
  call checksameppplvec(phiold_vec%data, vec, i)
  if(i .eq. 1) then
     print *, "phiold"
     if(allocated(phiold_vec%data)) deallocate(phiold_vec%data, STAT=i)
     allocate(phiold_vec%data(ivecsize))
     phiold_vec%data = 0.
     call updateids(vec, phiold_vec%data)
     return
  endif

  call checksameppplvec(vel_vec%data, vec, i)
  if(i .eq. 1) then
     print *, "vel"
     if(allocated(vel_vec%data)) deallocate(vel_vec%data, STAT=i)
     allocate(vel_vec%data(ivecsize))
     vel_vec%data = 0.
     call updateids(vec, vel_vec%data)
     return
  endif
      
  call checksameppplvec(velold_vec%data, vec, i)
  if(i .eq. 1) then
     print *, "velold"
     if(allocated(velold_vec%data)) deallocate(velold_vec%data, STAT=i)
     allocate(velold_vec%data(ivecsize))
     velold_vec%data = 0.
     call updateids(vec, velold_vec%data)
     return
  endif
   
     
  call checksameppplvec(den_vec%data, vec, i)
  if(i .eq. 1) then
     print *, "den"
     if(allocated(den_vec%data)) deallocate(den_vec%data, STAT=i)
     allocate(den_vec%data(ivecsize))
     den_vec%data = 0.
     call updateids(vec, den_vec%data)
     return
  endif
  
  call checksameppplvec(denold_vec%data, vec, i)
  if(i .eq. 1) then
     print *, "denold"
     if(allocated(denold_vec%data)) deallocate(denold_vec%data, STAT=i)
     allocate(denold_vec%data(ivecsize))
     denold_vec%data = 0.
     call updateids(vec, denold_vec%data)
     return
  endif
  
  call checksameppplvec(pres_vec%data, vec, i)
  if(i .eq. 1) then
     print *, "pres"
     if(allocated(pres_vec%data)) deallocate(pres_vec%data, STAT=i)
     allocate(pres_vec%data(ivecsize))
     pres_vec%data = 0.
     call updateids(vec, pres_vec%data)
     return
  endif
  
  call checksameppplvec(presold_vec%data, vec, i)
  if(i .eq. 1) then
     print *, "presold"
     if(allocated(presold_vec%data)) deallocate(presold_vec%data, STAT=i)
     allocate(presold_vec%data(ivecsize))
     presold_vec%data = 0.
     call updateids(vec, presold_vec%data)
     return
  endif
     
  call checksameppplvec(q4_vec%data, vec, i)
  if(i .eq. 1) then
     print *, "q4"
     if(allocated(q4_vec%data)) deallocate(q4_vec%data, STAT=i)
     allocate(q4_vec%data(ivecsize))
     q4_vec%data = 0.
     call updateids(vec, q4_vec%data)
     return
  endif

  call checksameppplvec(r4_vec%data, vec, i)
  if(i .eq. 1) then
     print *, "r4"
     if(allocated(r4_vec%data)) deallocate(r4_vec%data, STAT=i)
     allocate(r4_vec%data(ivecsize))
     r4_vec%data = 0.
     call updateids(vec, r4_vec%data)
     return
  endif

  call checksameppplvec(qn4_vec%data, vec, i)
  if(i .eq. 1) then
     print *, "qn4"
     if(allocated(qn4_vec%data)) deallocate(qn4_vec%data, STAT=i)
     allocate(qn4_vec%data(ivecsize))
     qn4_vec%data = 0.
     call updateids(vec, qn4_vec%data)
     return
  endif
  
  call checksameppplvec(qp4_vec%data, vec, i)
  if(i .eq. 1) then
     print *, "qp4"
     if(allocated(qp4_vec%data)) deallocate(qp4_vec%data, STAT=i)
     allocate(qp4_vec%data(ivecsize))
     qp4_vec%data = 0.
     call updateids(vec, qp4_vec%data)
     return
  endif

  call checksameppplvec(veln_vec%data, vec, i)
  if(i .eq. 1) then
     print *, "veln"
     if(allocated(veln_vec%data)) deallocate(veln_vec%data, STAT=i)
     allocate(veln_vec%data(ivecsize))
     veln_vec%data = 0.
     call updateids(vec, veln_vec%data)
     return
  endif
      
  call checksameppplvec(veloldn_vec%data, vec, i)
  if(i .eq. 1) then
     print *, "veloldn"
     if(allocated(veloldn_vec%data)) deallocate(veloldn_vec%data, STAT=i)
     allocate(veloldn_vec%data(ivecsize))
     veloldn_vec%data = 0.
     call updateids(vec, veloldn_vec%data)
     return
  endif

  call checksameppplvec(phip_vec%data, vec, i)
  if(i .eq. 1) then
     print *, "phip"
     if(allocated(phip_vec%data)) deallocate(phip_vec%data, STAT=i)
     allocate(phip_vec%data(ivecsize))
     phip_vec%data = 0.
     call updateids(vec, phip_vec%data)
     return
  endif
  
  call checksameppplvec(b1_phi%data, vec, i)
  if(i .eq. 1) then
     print *, "b1_phi"
     if(allocated(b1_phi%data)) deallocate(b1_phi%data, STAT=i)
     allocate(b1_phi%data(ivecsize))
     b1_phi%data = 0.
     call updateids(vec, b1_phi%data)
     return
  endif
  
  call checksameppplvec(b2_phi%data, vec, i)
  if(i .eq. 1) then
     print *, "b2_phi"
     if(allocated(b2_phi%data)) deallocate(b2_phi%data, STAT=i)
     allocate(b2_phi%data(ivecsize))
     b2_phi%data = 0.
     call updateids(vec, b2_phi%data)
     return
  endif
  
  call checksameppplvec(b1_vel%data, vec, i)
  if(i .eq. 1) then
     print *, "b1_vel"
     if(allocated(b1_vel%data)) deallocate(b1_vel%data, STAT=i)
     allocate(b1_vel%data(ivecsize))
     b1_vel%data = 0.
     call updateids(vec, b1_vel%data)
     return
  endif
  
  call checksameppplvec(b2_vel%data, vec, i)
  if(i .eq. 1) then
     print *, "b2_vel"
     if(allocated(b2_vel%data)) deallocate(b2_vel%data, STAT=i)
     allocate(b2_vel%data(ivecsize))
     b2_vel%data = 0.
     call updateids(vec, b2_vel%data)
     return
  endif
    
  call checksameppplvec(temporary_field%vec%data, vec, i)
  if(i .eq. 1) then
     print *, "temporary_vector"
     if(allocated(temporary_field%vec%data)) deallocate(temporary_field%vec%data, STAT=i)
     allocate(temporary_field%vec%data(ivecsize))
     temporary_field%vec%data = 0.
     call updateids(vec, temporary_field%vec%data)
     return
  endif

  print *, "Error: unknown vector"

end subroutine arrayresizevec
#endif
