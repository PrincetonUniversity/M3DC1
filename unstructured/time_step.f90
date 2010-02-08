subroutine matvecmult(imatrix,vin,vout)

  use basic
  use diagnostics

  implicit none

  integer, intent(in) :: imatrix
  vectype, dimension(*), intent(in) :: vin
  vectype, dimension(*), intent(out) :: vout

  real :: tstart, tend

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call matrixvectormult(imatrix,vin,vout)
  if(myrank.eq.0 .and. itimer.eq.1) then 
     call second(tend)
     t_mvm = t_mvm + tend - tstart
  end if
     
end subroutine matvecmult


!============================================================
! ONESTEP
! ~~~~~~~
! advance solution to time ntime+1
!============================================================
subroutine onestep

  use basic
  use diagnostics
  use arrays

  implicit none

  integer :: calc_matrices
  logical, save :: first_time = .true.

  real :: tstart, tend

  ! apply loop voltage
  fbound = fbound + dt*vloop/(2.*pi)


  ! Determine whether matrices should be re-calculated
  if(first_time &
       .or. (linear.eq.0 .and. mod(ntime,nskip).eq.0) &
       .or. (integrator.eq.1 .and. ntime.eq.1)) then
     calc_matrices = 1
  else
     calc_matrices = 0
  endif

  ! calculate matrices for time advance
  if(calc_matrices.eq.1) then
     if(myrank.eq.0 .and. iprint.eq.1) print *, "Defining matrices"
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

     ! in linear case, eliminate second-order terms from matrix
     call ludefall(1-istatic, idens, ipres, 1-iestatic)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_ludefall = t_ludefall + tend - tstart
     endif
     if(myrank.eq.0 .and. iprint.eq.1) print *, "Done defining matrices."
  endif


  ! copy field data to time-advance vectors
  if(myrank.eq.0 .and. iprint.eq.2) print *, "Importing time advance vectors.."
  call import_time_advance_vectors


  ! advance time
  if(myrank.eq.0 .and. iprint.eq.1) print *, "Advancing times..."
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  if(isplitstep.eq.1) then
     call split_step(calc_matrices)
  else
     call unsplit_step(calc_matrices)
  end if
  if(myrank.eq.0 .and. itimer.eq.1) then 
     call second(tend)
     if(iprint.ge.1) print *, "Time spent in *_step: ", tend-tstart
  end if

  time = time + dt
  dt = dt + ddt


  ! copy time advance vectors to field data
  if(myrank.eq.0 .and. iprint.eq.2) print *, "Exporting time advance vectors.."
  call export_time_advance_vectors


  ! Calculate all quantities derived from basic fields
  call derived_quantities(field, bf, 1)


  ! Conserve toroidal flux
  if(iconstflux.eq.1 .and. numvar.ge.2) then
     call conserve_flux
     tflux = tflux + gbound*area
  endif


  first_time = .false.

end subroutine onestep


!======================================================================
! import_time_advance_vectors
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
! construct vectors for time-advance
!======================================================================
subroutine import_time_advance_vectors
  use basic
  use arrays

  implicit none

  integer :: ibeginnv, iendplusonenv, ibegin, iendplusone
  integer :: numnodes, l

  call numnod(numnodes)

  do l=1,numnodes
     call assign_local_pointers(l)

     call entdofs(vecsize_phi, l, 0, ibeginnv, iendplusonenv)
     psi_v(ibeginnv+psi_off:ibeginnv+psi_off+5) = psi1_l
     if(numvar.ge.2) then
        bz_v(ibeginnv+bz_off:ibeginnv+bz_off+5) =  bz1_l
     endif
     if(numvar.ge.3) then
        pe_v(ibeginnv+ pe_off:ibeginnv+ pe_off+5) =  pe1_l
     endif

     call entdofs(vecsize_vel, l, 0, ibeginnv, iendplusonenv)   
     u_v(ibeginnv+  u_off:ibeginnv+  u_off+5) =   u1_l
     if(numvar.ge.2) then
        vz_v(ibeginnv+vz_off:ibeginnv+vz_off+5) =  vz1_l
     endif
     if(numvar.ge.3) then
        chi_v(ibeginnv+chi_off:ibeginnv+chi_off+5) = chi1_l
     endif

     if(ipres.eq.1) then
        call entdofs(vecsize_p, l, 0, ibegin, iendplusone)
        p_v(ibegin+p_off:ibegin+p_off+5) =   p1_l
     end if
     if(idens.eq.1) then
        call entdofs(vecsize_n, l, 0, ibegin, iendplusone)
        den_v(ibegin+den_off:ibegin+den_off+5) = den1_l
     endif
  enddo
end subroutine import_time_advance_vectors


!======================================================================
! export_time_advance_vectors
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
! extract values from time-advance vectors
!======================================================================
subroutine export_time_advance_vectors
  use basic
  use arrays

  implicit none

  integer :: ibeginnv, iendplusonenv, ibegin, iendplusone
  integer :: numnodes, l

  call numnod(numnodes)

  do l=1,numnodes
     call assign_local_pointers(l)

     call entdofs(vecsize_phi, l, 0, ibeginnv, iendplusonenv)
     psi1_l = psi_v(ibeginnv+psi_off:ibeginnv+psi_off+5)
     if(numvar.ge.2) then
        bz1_l = bz_v(ibeginnv+bz_off:ibeginnv+bz_off+5)
     endif
     if(numvar.ge.3) then
         pe1_l =  pe_v(ibeginnv+ pe_off:ibeginnv+ pe_off+5)
     endif

     call entdofs(vecsize_vel, l, 0, ibeginnv, iendplusonenv)
     u1_l = u_v(ibeginnv+  u_off:ibeginnv+  u_off+5)   
     if(numvar.ge.2) then
        vz1_l = vz_v(ibeginnv+vz_off:ibeginnv+vz_off+5)  
     endif
     if(numvar.ge.3) then
        chi1_l = chi_v(ibeginnv+chi_off:ibeginnv+chi_off+5)
     endif


     if(ipres.eq.1) then
        call entdofs(vecsize_p, l, 0, ibegin, iendplusone)
        p1_l =   p_v(ibegin  +  p_off:ibegin  +  p_off+5)
     end if
     if(idens.eq.1) then
        call entdofs(vecsize_n, l, 0, ibegin, iendplusone)
        den1_l = den_v(ibegin  +den_off:ibegin  +den_off+5)
     endif
  enddo
end subroutine export_time_advance_vectors

!============================================================
! SPILT_STEP
! ~~~~~~~~~~
! advance fields using split time step
!============================================================
subroutine split_step(calc_matrices)

  use p_data
  use t_data
  use basic
  use arrays
  use sparse
  use diagnostics

  implicit none
#include "finclude/petsc.h"

  integer, intent(in) :: calc_matrices
  real :: tstart, tend, t_bound
  integer :: i, l, numnodes, ndofs
  integer :: ibegin, iendplusone, ibeginnv, iendplusonenv
  integer :: jer, ier
  integer, allocatable:: itemp(:)
  vectype, allocatable :: temp(:), vel_temp(:), veln_temp(:), veloldn_temp(:)
  real :: x, z

  PetscTruth :: flg_petsc, flg_solve2
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-ipetsc', flg_petsc,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve2', flg_solve2,ier)

  t_bound = 0.

  call numnod(numnodes)

  ! Store current-time velocity matrices for use in field advance
  veln = vel
  veloldn = velold


  if(istatic.eq.0) then

     ! Advance Velocity
     ! ================
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Velocity"
  
     ! d1matrix_sm * vel(n)
     call matvecmult(d1matrix_sm,vel,b1_vel)
  
     ! q1matrix_sm * phi(n)
     if(ipres.eq.1 .and. numvar.ge.3) then
        ! replace electron pressure with total pressure
        do l=1,numnodes
           call entdofs(vecsize_p,   l, 0, ibegin,   iendplusone)
           call entdofs(vecsize_phi, l, 0, ibeginnv, iendplusonenv)
           
           phip(ibeginnv   :ibeginnv+11) = phi(ibeginnv:ibeginnv+11)
           phip(ibeginnv+12:ibeginnv+17) = pres(ibegin:ibegin+5)
        enddo
        call matvecmult(q1matrix_sm, phip, b1_phi)
     else
        call matvecmult(q1matrix_sm, phi , b1_phi)
     endif
     
     ! r14matrix_sm * den(n)
     if(idens.eq.1 .and. (gravr.ne.0 .or. gravz.ne.0)) then
        ! make a larger vector that can be multiplied by a numvar=3 matrix
        phip = 0.
        do l=1,numnodes
           call assign_local_pointers(l)
           call entdofs(vecsize_phi, l, 0, ibeginnv, iendplusonenv)       
           phip(ibeginnv:ibeginnv+5) = den1_l
        enddo

        call matvecmult(r14matrix_sm,phip,b2_phi)
        b1_phi = b1_phi + b2_phi
     endif

     ! o1matrix_sm * bf(n)
     if(numvar.ge.2 .and. i3d.eq.1) then
        ! make a larger vector that can be multiplied by a numvar=3 matrix
        phip = 0.
        call copyvec(bf,1,1,phip,1,vecsize_phi)
        call matvecmult(o1matrix_sm,phip,b2_phi)
        b1_phi = b1_phi + b2_phi
     endif

     ! Construct right-hand side
     call numdofs(vecsize_vel,ndofs)
     allocate(itemp(ndofs)) ! this is used to make sure that we 
     ! don't double count the sum for periodic dofs
     itemp = 1
     do l=1,numnodes
        call entdofs(vecsize_vel, l, 0, ibegin, iendplusone)
        call entdofs(vecsize_phi, l, 0, ibeginnv, iendplusonenv)
        do i=0,iendplusone-ibegin-1
           b1_vel(ibegin+i) = b1_vel(ibegin+i) + itemp(ibegin+i) * &
                (b1_phi(ibeginnv+i) + r4(ibegin+i))

           itemp(ibegin+i) = 0
        enddo
     enddo
     deallocate(itemp)

  
     ! apply boundary conditions
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     if(calc_matrices.eq.1) then
        call boundary_vel(s1matrix_sm, b1_vel)
        call finalizematrix(s1matrix_sm)
     else
        call boundary_vel(0, b1_vel)
     endif
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_bound = t_bound + tend - tstart
     end if

  
     ! solve linear system with rhs in vtemp (note LU-decomp done first time)
     if(myrank.eq.0 .and. iprint.ge.1) print *, "solving velocity advance..."
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     if(flg_petsc.eq.PETSC_TRUE .and. flg_solve2.eq.PETSC_TRUE) then
        call solve2(s1matrix_sm, b1_vel, jer)
     else
        call solve(s1matrix_sm, b1_vel, jer)
     endif
     if(myrank.eq.0 .and. iprint.ge.1) print *, "done solve."
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_solve_v = t_solve_v + tend - tstart
     endif
     if(jer.ne.0) then
        write(*,*) 'Error in velocity solve', jer
        call safestop(42)
     endif
  
     if(integrator.eq.1 .and. ntime.gt.1) then
        b1_vel = (2.*b1_vel - velold)/3.
     endif
  
     !.....new velocity solution at time n+1 (or n* for second order advance)
     velold = vel
     vel = b1_vel
   
     ! apply smoothing operators
     ! ~~~~~~~~~~~~~~~~~~~~~~~~~
     call smooth_velocity(vel)
  else
     velold = vel
  end if

  
  ! Advance Density
  ! ===============
  if(idens.eq.1) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Density"
     
     ! r8matrix_sm * vel(n+1)
     call matvecmult(r8matrix_sm,vel,b1_vel)

     ! r8matrix_sm * vel(n-1)
     if(integrator.eq.1 .and. ntime.gt.1) then
        call matvecmult(r8matrix_sm,veloldn,b2_vel)
        b1_vel = 1.5*b1_vel + 0.5*b2_vel
     endif

     ! q8matrix_sm * vel(n)
     call matvecmult(q8matrix_sm,veln,b2_vel)
     b1_vel = b1_vel + b2_vel
          
     ! temp = d8matrix_sm * phi(n)
     call createvec(temp, numvar1_numbering)
     temp = 0.
     
     call matvecmult(d8matrix_sm,den,temp)
     
     call numdofs(vecsize_vel,ndofs)
     allocate(itemp(ndofs)) ! this is used to make sure that we 
     ! don't double count the sum for periodic dofs
     itemp = 1
     
     do l=1,numnodes
        call entdofs(vecsize_n,   l, 0, ibegin, iendplusone)
        call entdofs(vecsize_vel, l, 0, ibeginnv, iendplusonenv)
        do i=0,iendplusone-ibegin-1
           temp(ibegin+i) = temp(ibegin+i) + itemp(ibegin+i) * &
                (b1_vel(ibeginnv+i) + qn4(ibegin+i))
           
           itemp(ibegin+i) = 0
        enddo
     enddo
     deallocate(itemp)
     
     ! Insert boundary conditions
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     if(calc_matrices.eq.1) then
        call boundary_den(s8matrix_sm, temp)
        call finalizematrix(s8matrix_sm)
     else
        call boundary_den(0, temp)
     endif
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_bound = t_bound + tend - tstart
     end if

     
     if(myrank.eq.0 .and. iprint.ge.1) print *, " solving..."
     
     
     ! solve linear system...LU decomposition done first time
     ! -- okay to here      call printarray(temp, 150, 0, 'vtemp on')
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     if(flg_petsc.eq.PETSC_TRUE .and. flg_solve2.eq.PETSC_TRUE) then
        call solve2(s8matrix_sm, temp, jer)
     else
        call solve(s8matrix_sm, temp, jer)
     endif
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_solve_n = t_solve_n + tend - tstart
     endif
     if(jer.ne.0) then
        write(*,*) 'Error in density solve', jer
        call safestop(29)
     endif
     
     ! new field solution at time n+1 (or n* for second order advance)
     if(integrator.eq.1 .and. ntime.gt.1) then
        temp = (2.*temp - denold)/3.
     endif
     denold = den
     den = temp
     call deletevec(temp)

     if(irecalc_eta.eq.1) then
        call export_time_advance_vectors
        call define_transport_coefficients
     end if
  endif
     
  !
  ! Advance Pressure
  ! ================
  if(ipres.eq.1) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Pressure"
     
     ! r9matrix_sm * vel(n+1)
     call matvecmult(r9matrix_sm,vel,b1_vel)
        
     ! r9matrix_sm * vel(n-1)
     if(integrator.eq.1 .and. ntime.gt.1) then
        call matvecmult(r9matrix_sm,veloldn,b2_vel)
        b1_vel = 1.5*b1_vel + 0.5*b2_vel
     endif
     
     ! q9matrix_sm * vel(n)
     call matvecmult(q9matrix_sm,veln,b2_vel)
     b1_vel = b1_vel + b2_vel
     
     ! temp = d8matrix_sm * pres(n)
     call createvec(temp, numvar1_numbering)
     temp = 0.
     call matvecmult(d9matrix_sm,pres,temp)
     

     ! Construct right-hand side
     call numdofs(vecsize_vel,ndofs)
     allocate(itemp(ndofs)) ! this is used to make sure that we 
     ! don't double count the sum for periodic dofs
     itemp = 1
     do l=1,numnodes
        call entdofs(vecsize_p, l, 0, ibegin, iendplusone)
        call entdofs(vecsize_vel, l, 0, ibeginnv, iendplusonenv)
        do i=0,iendplusone-ibegin-1
           temp(ibegin+i) = temp(ibegin+i) + itemp(ibegin+i) * &
                (b1_vel(ibeginnv+i) + qp4(ibegin+i))

           itemp(ibegin+i) = 0
        enddo
     enddo
     deallocate(itemp)
     
     ! Insert boundary conditions
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     if(calc_matrices.eq.1) then
        call boundary_pres(s9matrix_sm, temp)
        call finalizematrix(s9matrix_sm)
     else
        call boundary_pres(0, temp)
     endif
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_bound = t_bound + tend - tstart
     end if

     
     ! solve linear system...LU decomposition done first time
     ! -- okay to here      call printarray(temp, 150, 0, 'vtemp on')
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call solve(s9matrix_sm, temp, jer)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_solve_p = t_solve_p + tend - tstart
     endif
     if(jer.ne.0) then
        write(*,*) 'Error in pressure solve', jer
        call safestop(29)
     endif
     
     ! new field solution at time n+1 (or n* for second order advance)
     if(integrator.eq.1 .and. ntime.gt.1) then
        temp = (2.*temp - presold)/3.
     endif
     presold = pres
     pres = temp
     call deletevec(temp)
  endif

     
  !
  ! Advance Fields
  ! ==============

  if(iestatic.eq.0) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Fields"
   
     ! make new velocity vectors of the correct length
     call createvec(vel_temp, vecsize_phi)
     call createvec(veln_temp, vecsize_phi)
     call createvec(veloldn_temp, vecsize_phi)

     call copyfullvec(vel, vecsize_vel, vel_temp, vecsize_phi)
     call copyfullvec(veln, vecsize_vel, veln_temp, vecsize_phi)
     call copyfullvec(veloldn, vecsize_vel, veloldn_temp, vecsize_phi)

     ! r2matrix_sm * vel(n+1)
     call matvecmult(r2matrix_sm,vel_temp,b1_phi)
   
     ! r2matrix_sm * vel(n-1)
     if(integrator.eq.1 .and. ntime.gt.1) then
        call matvecmult(r2matrix_sm,veloldn_temp,b2_phi)
        b1_phi = 1.5*b1_phi + 0.5*b2_phi
     endif

     ! q2matrix_sm * vel(n)
     call matvecmult(q2matrix_sm,veln_temp,b2_phi)
     b1_phi = -b1_phi + b2_phi

     ! d2matrix_sm * phi(n)
     call matvecmult(d2matrix_sm,phi,b2_phi)
     b1_phi = b1_phi + b2_phi

     ! Include linear n^-1 terms and B^-2 terms
     if(eqsubtract.eq.1) then
        
        ! make a larger vector that can be multiplied by a numvar=3 matrix
        phip = 0.
        do l=1,numnodes
           call entdofs(num_fields, l, 0, ibegin, iendplusone)
           call entdofs(vecsize_phi, l, 0, ibeginnv, iendplusonenv)
           
           i = ibegin+6*(den_g-1)
           call nodcoord(l,x,z)

           call calc_ni(phip(ibeginnv:ibeginnv+5),field0(i:i+5),field(i:i+5))
           call calc_b2i(phip(ibeginnv+6:ibeginnv+11),x, &
                field0(ibegin+6*(psi_g-1):ibegin+6*(psi_g-1)+5), &
                field (ibegin+6*(psi_g-1):ibegin+6*(psi_g-1)+5), &
                field0(ibegin+6*( bz_g-1):ibegin+6*( bz_g-1)+5), &
                field (ibegin+6*( bz_g-1):ibegin+6*( bz_g-1)+5))
        enddo
        call matvecmult(q42matrix_sm,phip,b2_phi)
        b1_phi = b1_phi + b2_phi
     endif

  
     ! Include linear f terms
     if(numvar.ge.2 .and. i3d.eq.1) then
        ! b2vector = r15 * bf(n)

        ! make a larger vector that can be multiplied by a numvar=3 matrix
        phip = 0.        
        call copyvec(bf,1,1,phip,1,vecsize_phi)
        call matvecmult(o2matrix_sm,phip,b2_phi)
        b1_phi = b1_phi + b2_phi
     endif

     b1_phi = b1_phi + q4
       
     ! Insert boundary conditions
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     if(calc_matrices.eq.1) then
        call boundary_mag(s2matrix_sm, b1_phi)
        call finalizematrix(s2matrix_sm)
     else 
        call boundary_mag(0, b1_phi)
     endif
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_bound = t_bound + tend - tstart
     end if


!     call writematrixtofile(s2matrix_sm, 2)

     ! solve linear system...LU decomposition done first time
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

     if(flg_petsc.eq.PETSC_TRUE .and. flg_solve2.eq.PETSC_TRUE) then
        call solve2(s2matrix_sm, b1_phi, jer)
     else
        call solve(s2matrix_sm, b1_phi, jer)
     endif
     
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_solve_b = t_solve_b + tend - tstart
     endif
     if(jer.ne.0) then
        write(*,*) 'Error in field solve', jer
        call safestop(29)
     endif
  
     ! new field solution at time n+1 (or n* for second order advance)
     if(integrator.eq.1 .and. ntime.gt.1) then
        b1_phi = (2.*b1_phi - phiold)/3.
     endif

     ! Iterate field solve using re-defined transport coefficients
     if(iteratephi.eq.1) then
      
        ! temporarily advance fields to new values
        b2_phi = phi
        phi = b1_phi
        call export_time_advance_vectors
        ! redefine transport coefficients with new den/pe values
        call lcfs(field,psi_g,num_fields)
        call define_transport_coefficients
        ! revert fields to old values
        phi = b2_phi
        
        ! recalculate field advance matrix
        ! (advanced velocity variables will be used in defining matrix)
        if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
        call ludefall(0, 0, 0, 1)
        if(myrank.eq.0 .and. itimer.eq.1) then
           call second(tend)
           t_ludefall = t_ludefall + tend - tstart
        endif

        if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Fields Again"

        ! r2matrix_sm * vel(n+1)
        call matvecmult(r2matrix_sm,vel_temp,b1_phi)
   
        ! r2matrix_sm * vel(n-1)
        if(integrator.eq.1 .and. ntime.gt.1) then
           call matvecmult(r2matrix_sm,veloldn_temp,b2_phi)
           b1_phi = 1.5*b1_phi + 0.5*b2_phi
        endif

        ! q2matrix_sm * vel(n)
        call matvecmult(q2matrix_sm,veln_temp,b2_phi)
        b1_phi = -b1_phi + b2_phi

        ! d2matrix_sm * phi(n)
        call matvecmult(d2matrix_sm,phi,b2_phi)
        b1_phi = b1_phi + b2_phi
      
        ! Include linear f terms
        if(numvar.ge.2 .and. i3d.eq.1) then
           ! b2vector = r15 * bf(n)
           
           ! make a larger vector that can be multiplied by a numvar=3 matrix
           phip = 0.        
           call copyvec(bf,1,1,phip,1,vecsize_phi)
           call matvecmult(o2matrix_sm,phip,b2_phi)
           b1_phi = b1_phi + b1_phi
        endif

        b1_phi = b1_phi + q4
       
        ! Insert boundary conditions
        if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
        if(calc_matrices.eq.1) then
           call boundary_mag(s2matrix_sm, b1_phi)
           call finalizematrix(s2matrix_sm)
        else 
           call boundary_mag(0, b1_phi)
        endif
        if(myrank.eq.0 .and. itimer.eq.1) then
           call second(tend)
           t_bound = t_bound + tend - tstart
        end if

  
        ! solve linear system...LU decomposition done first time
        if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
        
        if(flg_petsc.eq.PETSC_TRUE .and. flg_solve2.eq.PETSC_TRUE) then
           call solve2(s2matrix_sm, b1_phi, jer)
        else
           call solve(s2matrix_sm, b1_phi, jer)
        endif
        
        if(myrank.eq.0 .and. itimer.eq.1) then
           call second(tend)
           t_solve_b = t_solve_b + tend - tstart
        endif
        if(jer.ne.0) then
           write(*,*) 'Error in field solve', jer
           call safestop(29)
        endif
  
        ! new field solution at time n+1 (or n* for second order advance)
        if(integrator.eq.1 .and. ntime.gt.1) then
           b1_phi = (2.*b1_phi - phiold)/3.
        endif
     end if

     phiold = phi
     phi = b1_phi

     call deletevec(vel_temp)
     call deletevec(veln_temp)
     call deletevec(veloldn_temp)

     ! apply smoothing operators
     ! ~~~~~~~~~~~~~~~~~~~~~~~~~
     call smooth_fields(phi) 
  end if


  if(myrank.eq.0 .and. iprint.eq.1 .and. itimer.eq.1) then
     print *, " split_step: Time solve: ", &
          t_solve_v + t_solve_b + t_solve_n + t_solve_p
     print *, " split_step: Time bcs: ", t_bound

  end if
  
end subroutine split_step


!============================================================
! UNSPILT_STEP
! ~~~~~~~~~~~~
! advance fields using unsplit time step
!============================================================
subroutine unsplit_step(calc_matrices)

  use p_data
  use t_data
  use basic
  use arrays
  use sparse
  use diagnostics

  implicit none
#include "finclude/petsc.h" 

  integer, intent(in) :: calc_matrices
  integer :: jer, ier
  
  real :: tstart, tend
  PetscTruth :: flg_petsc, flg_solve2
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-ipetsc', flg_petsc,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve2', flg_solve2,ier)

  if(myrank.eq.0 .and. iprint.ge.1) print *, "Solving matrix equation..."
  
  ! vtemp = d1matrix_sm * phi(n)
  call matvecmult(d1matrix_sm,phi,b1_phi) 
  b1_phi = b1_phi + q4

  ! Include linear f terms
  if(numvar.ge.2 .and. i3d.eq.1) then
     ! b2vector = r15 * bf(n)
     
     ! make a larger vector that can be multiplied by a vecsize matrix
     phip = 0.
     call copyvec(bf,1,1,phip,1,vecsize_phi)
     call matvecmult(o1matrix_sm,phip,b2_phi)
     b1_phi = b1_phi + b2_phi
  endif
    
  ! Insert boundary conditions
  if(calc_matrices.eq.1) then
     call boundary_mag(s1matrix_sm, b1_phi)
     call boundary_vel(s1matrix_sm, b1_phi)
     if(idens.eq.1) call boundary_den(s1matrix_sm, b1_phi)
     if(ipres.eq.1) call boundary_pres(s1matrix_sm, b1_phi)
     call finalizematrix(s1matrix_sm)
  else 
     call boundary_mag(0, b1_phi)
     call boundary_vel(0, b1_phi)
     if(idens.eq.1) call boundary_den(0, b1_phi)
     if(ipres.eq.1) call boundary_pres(0, b1_phi)
  endif
  
  ! solve linear system...LU decomposition done first time
  if(myrank.eq.0 .and. iprint.eq.1) print *, "solving.."
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  if(flg_petsc.eq.PETSC_TRUE .and. flg_solve2.eq.PETSC_TRUE) then
     call solve2(s1matrix_sm, b1_phi, jer)
  else
     call solve(s1matrix_sm, b1_phi, jer)
  endif
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_solve_b = t_solve_b + tend - tstart
  endif
  if(jer.ne.0) then
     write(*,*) 'Error in field solve', jer
     call safestop(29)
  endif
  if(myrank.eq.0 .and. iprint.eq.1) print *, "done solve."
  
  ! new field solution at time n+1
  if(integrator.eq.1 .and. ntime.gt.1) then
     b1_phi = (2.*b1_phi - phiold)/3.
  endif

  ! apply smoothing operators
  call smooth_velocity(phi)
  call smooth_fields(phi)

  if(iteratephi.eq.1 .and. iestatic.eq.0) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, "secondary advance..."

     ! temporarily advance fields to new values
     b2_phi = phi
     phi = b1_phi
     call export_time_advance_vectors
     ! redefine transport coefficients with new den/pe values
     call lcfs(field,psi_g,num_fields)
     call define_transport_coefficients
     ! revert fields to old values
     phi = b2_phi
        
     ! recalculate field advance matrix
     ! (advanced velocity variables will be used in defining matrix)
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call ludefall(1-istatic, idens, ipres, 1-iestatic)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_ludefall = t_ludefall + tend - tstart
     endif

     ! vtemp = d1matrix_sm * phi(n)
     call matvecmult(d1matrix_sm,phi,b1_phi) 
     b1_phi = b1_phi + q4

     ! Include linear f terms
     if(numvar.ge.2 .and. i3d.eq.1) then
        ! b2vector = r15 * bf(n)
     
        ! make a larger vector that can be multiplied by a vecsize matrix
        phip = 0.
        call copyvec(bf,1,1,phip,bf_i,vecsize_phi)
        call matvecmult(o1matrix_sm,phip,b2_phi)
        b1_phi = b1_phi + b2_phi
     endif
    
     ! Insert boundary conditions
     if(calc_matrices.eq.1) then
        call boundary_mag(s1matrix_sm, b1_phi)
        call boundary_vel(s1matrix_sm, b1_phi)
        if(idens.eq.1) call boundary_den(s1matrix_sm, b1_phi)
        if(ipres.eq.1) call boundary_pres(s1matrix_sm, b1_phi)
        call finalizematrix(s1matrix_sm)
     else 
        call boundary_mag(0, b1_phi)
        call boundary_vel(0, b1_phi)
        if(idens.eq.1) call boundary_den(0, b1_phi)
        if(ipres.eq.1) call boundary_pres(0, b1_phi)
     endif
     
     ! solve linear system...LU decomposition done first time
     if(myrank.eq.0 .and. iprint.eq.1) print *, "solving.."
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     if(flg_petsc.eq.PETSC_TRUE .and. flg_solve2.eq.PETSC_TRUE) then
        call solve2(s1matrix_sm, b1_phi, jer)
     else
        call solve(s1matrix_sm, b1_phi, jer)
     endif
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_solve_b = t_solve_b + tend - tstart
     endif
     if(jer.ne.0) then
        write(*,*) 'Error in field solve', jer
        call safestop(29)
     endif
     
     ! new field solution at time n+1
     if(integrator.eq.1 .and. ntime.gt.1) then
        b1_phi = (2.*b1_phi - phiold)/3.
     endif
     
     ! apply smoothing operators
     call smooth_velocity(phi)
     call smooth_fields(phi)
  endif
     
  phiold = phi
  phi = b1_phi

  if(myrank.eq.0 .and. iprint.ge.1) print *, "Done solving matrix equation."
end subroutine unsplit_step


subroutine calc_ni(vec, n0, n1)
  use arrays

  implicit none

  vectype, dimension(6), intent(out) :: vec
  vectype, dimension(6), intent(in) :: n0, n1

  vec(1) = -n1(1)/n0(1)**2
  vec(2) = &
       -(n1(2)*n0(1) - 2.*n1(1)*n0(2))/n0(1)**3
  vec(3) = &
       -(n1(3)*n0(1) - 2.*n1(1)*n0(3))/n0(1)**3
  vec(4) = &
       -n1(4)/n0(1)**2 &
       +2.*(2.*n1(2)*n0(2) &
       +   n1(1)*n0(4))/n0(1)**3 &
       -6.*n1(1)*n0(2)**2/n0(1)**4
  vec(5) = &
       -n1(5)/n0(1)**2 &
       +2.*(n1(2)*n0(3) + n1(3)*n0(2) &
       +n1(1)*n0(5))/n0(1)**3 &
       -6.*n1(1)*n0(2)*n0(3)/n0(1)**4
  vec(6) = &
       -n1(6)/n0(1)**2 &
       +2.*(2.*n1(3)*n0(3) &
       +   n1(1)*n0(6))/n0(1)**3 &
       -6.*n1(1)*n0(3)**2/n0(1)**4
end subroutine calc_ni

subroutine calc_b2i(vec, x, psi0, psi1, b0, b1)
  use arrays

  implicit none

  real, intent(in) :: x
  vectype, dimension(6), intent(out) :: vec
  vectype, dimension(6), intent(in) :: psi0, psi1, b0, b1

  vectype, dimension(6) :: b2i
  vectype :: b20

  b20 = psi0(2)**2 + psi0(3)**2 + b0(1)**2

  b2i(1) = x**2/b20
  b2i(2) = 2.*x/b20 &
       - (2.*x**2/b20**2)*(psi0(2)*psi0(4) + psi0(3)*psi0(5) + b0(1)*b0(2))
  b2i(3) = &
       - (2.*x**2/b20**2)*(psi0(2)*psi0(5) + psi0(3)*psi0(6) + b0(1)*b0(3))
  b2i(4:6) = 0.

  vec(1) = -2.*b2i(1)**2/x**2 * &
       (b0(1)*b1(1) + psi0(2)*psi1(2) + psi0(3)*psi1(3))
  vec(2) = -2.*b2i(2)**2/x**2 * &
       (b0(1)*b1(1) + psi0(2)*psi1(2) + psi0(3)*psi1(3)) &
           -2.*b2i(1)**2/x**2 * &
       (b0(2)*b1(1) + psi0(4)*psi1(2) + psi0(5)*psi1(3)  &
       +b0(1)*b1(2) + psi0(2)*psi1(4) + psi0(3)*psi1(5)) &
           +4.*b2i(1)**2/x**3 * &
       (b0(1)*b1(1) + psi0(2)*psi1(2) + psi0(3)*psi1(3))
  vec(3) = -2.*b2i(3)**2/x**2 * &
       (b0(1)*b1(1) + psi0(2)*psi1(2) + psi0(3)*psi1(3)) &
           -2.*b2i(1)**2/x**2 * &
       (b0(3)*b1(1) + psi0(5)*psi1(2) + psi0(6)*psi1(3)  &
       +b0(1)*b1(3) + psi0(2)*psi1(5) + psi0(3)*psi1(6))
  vec(4:6) = 0.

end subroutine calc_b2i
