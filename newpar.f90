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

  integer :: j, i, ier, numelms, numnodes
  integer :: ndofs, ibegin, iendplusone

  real :: tstart, tend
  real :: factor, hmin, hmax  

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
  ! initialize the SUPERLU process grid
  call initsolvers
  call MPI_Comm_size(MPI_COMM_WORLD,maxrank,ier)
  if (ier /= 0) then
     print *,'Error in MPI_Comm_size:',ier
     call safestop(1)
  endif
  ! initialize autopack
  call AP_INIT()

  if(myrank.eq.0 .and. itimer.ge.1) call second(tstart)
  call loadmesh("struct.dmg", "struct-dmg.sms")
  if(myrank.eq.0 .and. itimer.ge.1) then
     call second(tend)
     print *, 'Time spent in loadmesh: ', tend-tstart
  endif
  if(myrank.eq.0) then
     call date_and_time( datec, timec)
     write(*,1001) datec(1:4),datec(5:6),datec(7:8),                        &
          timec(1:2),timec(3:4),timec(5:8)
1001 format("M3D-C1 VERSION 1.0    DATE: ", a4,1x,a2,1x,a2,3x,               &
          "TIME: ",a2,":",a2,":",a4,/)
  endif

  if(myrank.eq.0) then 
#ifdef NEW_VELOCITY
     print*, "V = grad(U)xgrad(phi) + V grad(phi) + grad(chi)"
#else
     print*, "V = r^2 grad(U)xgrad(phi) + r V grad(phi) + grad(chi)"
#endif
  endif


  ! Output information about local dofs, nodes, etc.
  if(iprint.ge.1) then
     call numfac(numelms)
     call numnod(numnodes)
     call numprocdofs(1, j)
     call numdofs(1, ndofs)
     print *, 'proc, owned dofs, needed dofs',myrank, j, ndofs
     print *, 'proc, numnodes, numfaces', myrank, numnodes,numelms
  endif

  pi = acos(-1.)

  ! initialize needed variables and define geometry and triangles
  call init

  if(linear.eq.1) eqsubtract = 1

  ! calculate pfac (pe*pfac = electron pressure)
  if(ipres.eq.1) then
     pefac = 1.
  else
     if(p0.gt.0.) then
        pefac = (p0-pi0)/p0
     else
        pefac = 0.
     endif
     if(myrank.eq.0 .and. iprint.ge.1) print *, "pefac = ", pefac
  endif


  ! check time-integration options
  select case(integrator)
  case(1)
     ! For BDF2 integration, first timestep is Crank-Nicholson with thimp=1,
     ! and subsequent timesteps are BDF2.
     if(myrank.eq.0) print *, "Time integration: BDF2."
     thimp = 1.
     thimp_ohm = 1.
  case default
     if(myrank.eq.0) print *, "Time integration: Crank-Nicholson."
  end select


  if(irestart.eq.1) then
     ! Read restart file(s)

     if(myrank.eq.1 .and. iprint.ge.1) print *, 'Reading restart file(s)...'
     if(iglobalin.eq.1) then
        call rdrestartglobal
     else
        call rdrestart
     endif

  else
     ! Initialize with initial conditions

     ptot = 0.
     ntime = 0
     time = 0.
     if(myrank.eq.0 .and. iprint.ge.1) print *, 'defining initial conditions...'
     call initial_conditions
     if(myrank.eq.0 .and. iprint.ge.1) print *, 'done initial conditions'

     if(eqsubtract.eq.1) then
        fieldi = field
     else
        fieldi = field0
     endif

     ! correct for left-handed coordinates
     if(myrank.eq.0 .and. iprint.ge.1) &
          print *, "adjusting fields for left-handed coordinates"
     call flip_handedness

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
     
  endif                     !  end of the branch on restart/no restart

  ! initialize hdf5
  if(myrank.eq.0 .and. iprint.ge.1) print *, "Initializing HDF5."
  call hdf5_initialize(ier)
  if(ier.lt.0) then 
     print *, "Error initializing HDF5"
     call safestop(5)
  end if
  if(myrank.eq.0 .and. iprint.ge.1) print *, "Done initializing HDF5."
   
  ! create the newvar matrices
  if(myrank.eq.0 .and. iprint.ge.1) print *, "Generating newvar matrices..."
  call create_matrix(mass_matrix_lhs_dc,    NV_DCBOUND, NV_I_MATRIX,  NV_LHS)
  call create_matrix(mass_matrix_lhs,       NV_NOBOUND, NV_I_MATRIX,  NV_LHS)
  call create_matrix(gs_matrix_rhs_dc,      NV_DCBOUND, NV_GS_MATRIX, NV_RHS)
  if(numvar.ge.3 .and. hyperc.ne.0. .and. com_bc.eq.0) then
     call create_matrix(lp_matrix_rhs,      NV_NOBOUND, NV_LP_MATRIX, NV_RHS)
  endif
  if(numvar.ge.3 .and. hyperc.ne.0. .and. com_bc.eq.1) then
     call create_matrix(lp_matrix_rhs_dc,   NV_DCBOUND, NV_LP_MATRIX, NV_RHS)
  end if
  if(i3d.eq.1) then
     call create_matrix(poisson_matrix_lhs, NV_DCBOUND, NV_LP_MATRIX, NV_LHS)
     call create_matrix(bf_matrix_rhs_dc,   NV_DCBOUND, NV_BF_MATRIX, NV_RHS)
  endif
  if(gyro.eq.1 .and. numvar.ge.2) then
     call zeromultiplymatrix(gyro_torque_sm,icomplex,vecsize_vel)
     call finalizematrix(gyro_torque_sm)
  endif
  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, "Done generating newvar matrices."


  if(itimer.eq.1) call reset_timings

  ! Calculate all quantities derived from basic fields
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  call derived_quantities

  ! Adapt the mesh
  ! ~~~~~~~~~~~~~~
  if(iadapt.eq.1) then
!     call outputfield(phi, numvar, 0, ntime, 123) 
     if(maxrank .eq. 1) then
!        call outputfield(phi, numvar, 0, ntime, 123) 
!        call writefieldatnodes(resistivity, 1, 1) 
        factor = 0.2
        hmin = .005
        hmax = 0.05

        print *, 'adapting mesh...'
!!$#ifdef USECOMPLEX
        call hessianadapt(resistivity,1, 0, ntime, factor, hmin, hmax) 
!!$#else
!!$        call hessianadapt(resistivity,1, ntime, factor, hmin, hmax)
!!$#endif
        print *, 'done adapting.'
        call space(0)
        call tridef
     endif
  endif


  if(irestart.eq.0) then
     tflux0 = tflux
     totcur0 = totcur
  endif

  ! output simulation parameters
  if(irestart.eq.0) then
     if(myrank.eq.0 .and. iprint.ge.1) &
          print *, "Writing simulation parameters."
     call hdf5_write_parameters(ier)

     ! Output the equilibrium
     if(eqsubtract.eq.1) call hdf5_write_time_slice(1,ier)
  end if

  ! output initial conditions
  call output

  ! if there are no timesteps to calculate, then skip time loop
  if(ntimemax.le.ntime) go to 101


  ! main time loop
  ! ~~~~~~~~~~~~~~
  do ntime=ntime+1,ntimemax

     ! check for error
     if(ekin.ne.ekin .or. emag.ne.emag) then
        print *, "Error: energy is NaN"
        call safestop(3)
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
     if(itor.eq.1 .and. itaylor.eq.1) call control_pid


     ! Write output
     call output
  enddo ! ntime


101 continue


!     free memory from sparse matrices
  call deletematrix(mass_matrix_lhs)
  call deletematrix(mass_matrix_lhs_dc)
  call deletematrix(gs_matrix_rhs_dc)
  call deletematrix(lp_matrix_rhs)
  call deletematrix(lp_matrix_rhs_dc)
  call deletematrix(gsmatrix_sm)
  call deletematrix(s7matrix_sm)
  call deletematrix(s4matrix_sm)
  call deletematrix(s5matrix_sm)
  call deletematrix(s1matrix_sm)
  call deletematrix(s2matrix_sm)
  call deletematrix(d1matrix_sm)
  call deletematrix(d2matrix_sm)
  call deletematrix(d4matrix_sm)
  call deletematrix(q1matrix_sm)
  call deletematrix(r2matrix_sm)
  call deletematrix(q2matrix_sm)
  call deletematrix(r14matrix_sm)
  if(idens.eq.1) then
     call deletematrix(s8matrix_sm)
     call deletematrix(d8matrix_sm)
     call deletematrix(q8matrix_sm) 
     call deletematrix(r8matrix_sm)
  end if
  if(ipres.eq.1) then 
     call deletematrix(s9matrix_sm)
     call deletematrix(d9matrix_sm)
     call deletematrix(r9matrix_sm)
     call deletematrix(q9matrix_sm)
  endif
  if(i3d.eq.1) then
     call deletematrix(poisson_matrix_lhs)
     call deletematrix(bf_matrix_rhs_dc)
  end if
  if(gyro.eq.1) then
     call deletematrix(gyro_torque_sm)
  end if
  if(idens.eq.1 .and. linear.eq.1) then
     call deletematrix(q42matrix_sm)
  end if
#ifdef USECOMPLEX
  call deletematrix(o1matrix_sm)
  call deletematrix(o2matrix_sm)
#endif
  call deletesearchstructure()
!  free memory for numberings
  call deletedofnumbering(1)
  call deletedofnumbering(2)
  if(num_fields.gt.2) call deletedofnumbering(num_fields)
  if(vecsize_phi.gt.2 .and. vecsize_phi.ne.num_fields) then
     call deletedofnumbering(vecsize_phi)
  end if
  if(vecsize_vel.gt.2 .and. vecsize_vel.ne.num_fields &
       .and. vecsize_vel.ne.vecsize_phi) then
     call deletedofnumbering(vecsize_vel)
  end if

  call safestop(2)

end Program Reducedquintic

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
  if(mod(ntime,ntimepr).eq.0) then
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
! smooth
! ~~~~~~
!
! applies smoothing operators
!
! ======================================================================
subroutine smooth
  use basic
  use arrays
  use newvar_mod
  use diagnostics
  use sparse

  implicit none

  real :: tstart, tend
  integer :: numnodes

  if(hyperc.eq.0.) return

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

  call numnod(numnodes)
     
  ! smooth vorticity
  call newvar(mass_matrix_lhs_dc,vor,vel,1,vecsize_vel, &
       gs_matrix_rhs_dc,NV_DCBOUND)
  call smoother1(vor,vel,numnodes,vecsize_vel,1)

  ! smooth compression
  if(numvar.ge.3) then
     if(com_bc.eq.1) then
        call newvar(mass_matrix_lhs_dc,com,vel,3,vecsize_vel, &
             lp_matrix_rhs_dc,NV_DCBOUND)
     else
        call newvar(mass_matrix_lhs   ,com,vel,3,vecsize_vel, &
             lp_matrix_rhs,   NV_NOBOUND)
     endif
             
     call smoother3(com,vel,numnodes,vecsize_vel,3)
  endif

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_smoother = t_smoother + tend - tstart
  endif

end subroutine smooth


!======================================================================
! copyvec
! ~~~~~~~
! copies a field from inarr to outarr
!======================================================================
subroutine copyvec(inarr, inpos, insize, outarr, outpos, outsize)

  implicit none

  vectype, intent(in) :: inarr(*)
  integer, intent(in) :: insize, inpos
  vectype, intent(out) :: outarr(*)
  integer, intent(in) :: outsize, outpos

  integer :: ibegini, iendplusonei
  integer :: ibegino, iendplusoneo
  integer :: l, numnodes, in_i, out_i

  in_i = (inpos-1)*6
  out_i = (outpos-1)*6

  call numnod(numnodes)

  do l=1,numnodes
     call entdofs(insize, l, 0, ibegini, iendplusonei)
     call entdofs(outsize, l, 0, ibegino, iendplusoneo)
    
     outarr(ibegino+out_i:ibegino+out_i+5) = &
          inarr(ibegini+in_i:ibegini+in_i+5)
  enddo
end subroutine copyvec

! ======================================================================
! derived_quantities
! ~~~~~~~~~~~~~~~~~~
! calculates all derived quantities, including auxiliary fields
! and scalars
! ======================================================================
subroutine derived_quantities
  use basic
  use arrays
  use newvar_mod
  use diagnostics
  use sparse

  implicit none

  real :: tstart, tend

  ! Define auxiliary fields
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.1) print *, "Defining auxiliary fields."
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

  if(myrank.eq.0 .and. iprint.ge.1) print *, "-Transport coefficients"
  call define_transport_coefficients

  !   toroidal current
  if(myrank.eq.0 .and. iprint.ge.1) print *, "-Toroidal current"
  call newvar(mass_matrix_lhs_dc,jphi,field,psi_g,num_fields, &
       gs_matrix_rhs_dc,NV_DCBOUND)

  if(hyperc.ne.0.) then
     !   vorticity
     if(myrank.eq.0 .and. iprint.ge.1) print *, "-Vorticity"
     call newvar(mass_matrix_lhs_dc,vor,field,u_g,num_fields, &
          gs_matrix_rhs_dc,NV_DCBOUND)

     !   compression
     if(numvar.ge.3) then
        if(myrank.eq.0 .and. iprint.ge.1) print *, "-Compression"
        if(com_bc.eq.1) then
           call newvar(mass_matrix_lhs_dc,com,field,chi_g,num_fields, &
                lp_matrix_rhs_dc,NV_DCBOUND)
        else
           call newvar(mass_matrix_lhs,   com,field,chi_g,num_fields, &
                lp_matrix_rhs,   NV_NOBOUND)
        endif
     else
        com = 0.
     endif
  endif

  ! vector potential stream function
  if(i3d.eq.1) then
     call newvar(poisson_matrix_lhs,bf,field,bz_g,num_fields, &
          bf_matrix_rhs_dc,NV_DCBOUND)
  endif

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_aux = t_aux + tend - tstart
  endif
  

!!$  ! find lcfs
!!$  ! ~~~~~~~~~
!!$  call lcfs(psi+psi0,1)

  ! calculate scalars
  ! ~~~~~~~~~~~~~~~~~
  if(myrank.eq.0 .and. iprint.ge.1) print *, "Calculating scalars"
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call calculate_scalars
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_sources = t_sources + tend - tstart
  endif

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
  
  integer :: ndofs, j

  ! adjust toroidal field and boundary condition to conserve toroidal flux
  if(numvar.ge.2 .and. iconstflux.eq.1) then
     gbound = (tflux0-tflux)/area
     call numdofs(1, ndofs)
     do j=1+(bz_g-1)*6,ndofs,6*num_fields
        field(j) = field(j) + gbound
     enddo
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
end subroutine flip_handedness


