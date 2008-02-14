!============================================================
! ONESTEP
! ~~~~~~~
! advance solution to time ntime+1
!============================================================
subroutine onestep

  use p_data
  use t_data
  use basic
  use arrays
  use sparse
  use diagnostics
  use gradshafranov

  implicit none

  integer :: l, jer, i
  integer :: ibegin, iendplusone, ibeginnv, iendplusonenv
  integer, allocatable:: itemp(:)
  integer :: ndofs, numnodes

  integer :: calc_matrices
  logical :: first_time = .true.

  real :: tstart, tend
  vectype, allocatable :: temp(:), temp2(:)
  
  call numnod(numnodes)

  ! apply loop voltage
  ! ~~~~~~~~~~~~~~~~~~
  fbound = fbound + dt*vloop/(2.*pi)


  ! Determine whether matrices should be re-calculated
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(first_time &
       .or. (linear.eq.0 .and. mod(ntime,nskip).eq.0) &
       .or. (integrator.eq.1 .and. ntime.eq.1)) then
     calc_matrices = 1
  else
     calc_matrices = 0
  endif

  ! calculate matrices for time advance
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(calc_matrices.eq.1) then 
     if(myrank.eq.0 .and. iprint.eq.1) print *, "Defining matrices"
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call ludefall
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_ludefall = t_ludefall + tend - tstart
     endif
     if(myrank.eq.0 .and. iprint.eq.1) print *, "Done defining matrices."
  endif


  ! construct vectors for time-advance
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do l=1,numnodes
     call entdofs(vecsize1, l, 0, ibegin, iendplusone)
     call entdofs(vecsize , l, 0, ibeginnv, iendplusonenv)
     
     call assign_local_pointers(l)

       u_v(ibeginnv+  u_off:ibeginnv+  u_off+5) =   u1_l
     psi_v(ibeginnv+psi_off:ibeginnv+psi_off+5) = psi1_l
     if(numvar.ge.2) then
        vz_v(ibeginnv+vz_off:ibeginnv+vz_off+5) =  vz1_l
        bz_v(ibeginnv+bz_off:ibeginnv+bz_off+5) =  bz1_l
     endif
     if(numvar.ge.3) then
        chi_v(ibeginnv+chi_off:ibeginnv+chi_off+5) = chi1_l
         pe_v(ibeginnv+ pe_off:ibeginnv+ pe_off+5) =  pe1_l
         if(ipres.eq.1) p_v(ibegin  +  p_off:ibegin  +  p_off+5) =   p1_l
     endif
     if(idens.eq.1) then
        den_v(ibegin+den_off:ibegin+den_off+5) = den1_l
     endif
  enddo

  ! advance time
  ! ~~~~~~~~~~~~
  if(isplitstep.eq.1) then
     call split_step(calc_matrices)
  else
     call unsplit_step(calc_matrices)
  end if

  ! extract values from time-advance vectors
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  do l=1,numnodes
     call entdofs(vecsize1, l, 0, ibegin, iendplusone)
     call entdofs(vecsize , l, 0, ibeginnv, iendplusonenv)
     
     call assign_local_pointers(l)

       u1_l = u_v(ibeginnv+  u_off:ibeginnv+  u_off+5)   
     psi1_l = psi_v(ibeginnv+psi_off:ibeginnv+psi_off+5)
     if(numvar.ge.2) then
        vz1_l = vz_v(ibeginnv+vz_off:ibeginnv+vz_off+5)  
        bz1_l = bz_v(ibeginnv+bz_off:ibeginnv+bz_off+5)
     endif
     if(numvar.ge.3) then
        chi1_l = chi_v(ibeginnv+chi_off:ibeginnv+chi_off+5)
         pe1_l =  pe_v(ibeginnv+ pe_off:ibeginnv+ pe_off+5)
         if(ipres.eq.1) p1_l =   p_v(ibegin  +  p_off:ibegin  +  p_off+5)
     endif
     if(idens.eq.1) then
        den1_l = den_v(ibegin  +den_off:ibegin  +den_off+5)
     endif
  enddo


  ! Calculate all quantities derived from basic fields
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  call derived_quantities


  ! Conserve toroidal flux
  ! ~~~~~~~~~~~~~~~~~~~~~~
  if(iconstflux.eq.1 .and. numvar.ge.2) then
     call conserve_flux
     tflux = tflux + gbound*area
  endif


  first_time = .false.

end subroutine onestep



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

  integer, intent(in) :: calc_matrices
  real :: tstart, tend
  integer :: i, l, numnodes, ndofs
  integer :: ibegin, iendplusone, ibeginnv, iendplusonenv
  integer :: jer
  integer, allocatable:: itemp(:)
  vectype, allocatable :: temp(:), temp2(:)

  call numnod(numnodes)

  phip = phi

  ! Store current-time velocity matrices for use in field advance
  veln = vel
  veloldn = velold


  if(istatic.eq.0) then
     ! Advance Velocity
     ! ================
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Velocity"
  
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  
     ! b1vector = q1matrix_sm * phi(n)
     if(ipres.eq.1 .and. numvar.ge.3) then
        ! replace electron pressure with total pressure
        do l=1,numnodes
           call entdofs(1, l, 0, ibegin, iendplusone)
           call entdofs(vecsize, l, 0, ibeginnv, iendplusonenv)
           
           phip(ibeginnv   :ibeginnv+11) = phi(ibeginnv:ibeginnv+11)
           phip(ibeginnv+12:ibeginnv+17) = pres(ibegin:ibegin+5)
        enddo
        call matrixvectormult(q1matrix_sm, phip, b1vector)
     else
        call matrixvectormult(q1matrix_sm, phi , b1vector)
     endif
  
     ! vtemp = d1matrix_sm * vel(n)
     vtemp = 0.
     call matrixvectormult(d1matrix_sm,vel,vtemp)
  
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_mvm = t_mvm + tend - tstart
     endif
  
     vtemp = vtemp + b1vector + r4
  
     ! Include linear density terms
     if(idens.eq.1 .and. (gravr.ne.0 .or. gravz.ne.0)) then
        ! b2vector = r14 * den(n)
     
        ! make a larger vector that can be multiplied by a numvar=3 matrix
        phip = 0.
        do l=1,numnodes
           call assign_local_pointers(l)
           call entdofs(vecsize, l, 0, ibeginnv, iendplusonenv)       
           phip(ibeginnv:ibeginnv+5) = den1_l
        enddo

        call matrixvectormult(r14matrix_sm,phip,b2vector)
        vtemp = vtemp + b2vector
     endif

     ! Include linear f terms
     if(numvar.ge.2 .and. i3d.eq.1) then
        ! b2vector = r15 * bf(n)
     
        ! make a larger vector that can be multiplied by a numvar=3 matrix
        phip = 0.
        do l=1,numnodes
           call entdofs(1, l, 0, ibegin, iendplusone)
           call entdofs(vecsize, l, 0, ibeginnv, iendplusonenv)
           
           phip(ibeginnv  :ibeginnv+5) = bf(ibegin:ibegin+5)
        enddo
        call matrixvectormult(o1matrix_sm,phip,b2vector)
        vtemp = vtemp + b2vector
     endif

  
     ! apply boundary conditions
     if(calc_matrices.eq.1) then
        call boundary_vel(s1matrix_sm, vtemp)
        call finalizematrix(s1matrix_sm)
     else
        call boundary_vel(0, vtemp)
     endif
  
     ! solve linear system with rhs in vtemp (note LU-decomp done first time)
     if(myrank.eq.0) print *, "solving velocity advance..."
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call solve(s1matrix_sm, vtemp, jer)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_solve_v = t_solve_v + tend - tstart
     endif
     if(jer.ne.0) then
        write(*,*) 'Error in velocity solve', jer
        call safestop(42)
     endif
  
     if(integrator.eq.1 .and. ntime.gt.1) then
        vtemp = (2.*vtemp - velold)/3.
     endif
  
     !.....new velocity solution at time n+1 (or n* for second order advance)
     velold = vel
     vel = vtemp
    
     ! apply smoothing operators
     ! ~~~~~~~~~~~~~~~~~~~~~~~~~
     call smooth
  else
     velold = vel
  end if

  
  ! Advance Density
  ! ===============
  if(idens.eq.1) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Density"
     
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     
     ! b2vector = r8matrix_sm * vel(n+1)
     call matrixvectormult(r8matrix_sm,vel,b2vector)
     
     ! b3vector = q8matrix_sm * vel(n)
     call matrixvectormult(q8matrix_sm,veln,b3vector)
     
     ! b1vector = r8matrix_sm * vel(n-1)
     if(integrator.eq.1 .and. ntime.gt.1) then
        b2vector = 1.5*b2vector
        call matrixvectormult(r8matrix_sm,veloldn,b1vector)
        b1vector = 0.5*b1vector
     else
        b1vector = 0.
     endif
     
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_mvm = t_mvm + tend - tstart
     endif
     
     ! temp = d8matrix_sm * phi(n)
     call createvec(temp, numvar1_numbering)
     temp = 0.
     
     call matrixvectormult(d8matrix_sm,den,temp)
     
     call numdofs(numvar,ndofs)
     allocate(itemp(ndofs)) ! this is used to make sure that we 
     ! don't double count the sum for periodic dofs
     itemp = 1
     
     do l=1,numnodes
        call entdofs(1, l, 0, ibegin, iendplusone)
        call entdofs(vecsize, l, 0, ibeginnv, iendplusonenv)
        do i=0,iendplusone-ibegin-1
           temp(ibegin+i) = temp(ibegin+i) + itemp(ibegin+i) * &
                (b2vector(ibeginnv+i) + b3vector(ibeginnv+i) + qn4(ibegin+i)) &
                +b1vector(ibeginnv+i)
           
           itemp(ibegin+i) = 0
        enddo
     enddo
     deallocate(itemp)
     
     ! Insert boundary conditions
     if(calc_matrices.eq.1) then
        call boundary_den(s8matrix_sm, temp)
        call finalizematrix(s8matrix_sm)
     else
        call boundary_den(0, temp)
     endif
     
     if(myrank.eq.0 .and. iprint.ge.1) print *, " solving..."
     
     
     ! solve linear system...LU decomposition done first time
     ! -- okay to here      call printarray(temp, 150, 0, 'vtemp on')
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call solve(s8matrix_sm, temp, jer)
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
  endif
     
  !
  ! Advance Pressure
  ! ================
  if(ipres.eq.1) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Pressure"
     
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     
     ! b2vector = r9matrix_sm * vel(n+1)
     call matrixvectormult(r9matrix_sm,vel,b2vector)
     
     ! b3vector = q9matrix_sm * vel(n)
     call matrixvectormult(q9matrix_sm,veln,b3vector)
     
     ! b1vector = r9matrix_sm * vel(n-1)
     if(integrator.eq.1 .and. ntime.gt.1) then
        b2vector = 1.5*b2vector
        call matrixvectormult(r9matrix_sm,veloldn,b1vector)
        b1vector = 0.5*b1vector
     else
        b1vector = 0.
     endif
     
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_mvm = t_mvm + tend - tstart
     endif
     
     ! temp = d8matrix_sm * pres(n)
     call createvec(temp, numvar1_numbering)
     temp = 0.
     call matrixvectormult(d9matrix_sm,pres,temp)
     
     call numdofs(vecsize,ndofs)
     allocate(itemp(ndofs)) ! this is used to make sure that we 
     ! don't double count the sum for periodic dofs
     itemp = 1
     
     do l=1,numnodes
        call entdofs(1, l, 0, ibegin, iendplusone)
        call entdofs(vecsize, l, 0, ibeginnv, iendplusonenv)
        do i=0,iendplusone-ibegin-1
           temp(ibegin+i) = temp(ibegin+i) + itemp(ibegin+i) * &
                (b2vector(ibeginnv+i) + b3vector(ibeginnv+i) + qp4(ibegin+i)) &
                +b1vector(ibeginnv+i)
           itemp(ibegin+i) = 0
        enddo
     enddo
     deallocate(itemp)
     
     ! Insert boundary conditions
     if(calc_matrices.eq.1) then
        call boundary_pres(s9matrix_sm, temp)
        call finalizematrix(s9matrix_sm)
     else
        call boundary_pres(0, temp)
     endif
     
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
  
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  
     ! b2vector = r2matrix_sm * vel(n+1)
     call matrixvectormult(r2matrix_sm,vel,b2vector)
  
     ! b3vector = q2matrix_sm * vel(n)
     call matrixvectormult(q2matrix_sm,veln,b3vector)
  
     ! b1vector = r2matrix_sm * vel(n-1)
     if(integrator.eq.1 .and. ntime.gt.1) then
        b2vector = 1.5*b2vector
        call matrixvectormult(r2matrix_sm,veloldn,b1vector)
        b1vector = 0.5*b1vector
     else
        b1vector = 0.
     endif
  
     ! vtemp = d2matrix_sm * phi(n)
     vtemp = 0.
     call matrixvectormult(d2matrix_sm,phi,vtemp)
     
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_mvm = t_mvm + tend - tstart
     endif
  
     vtemp = vtemp + b2vector + b3vector + q4  + b1vector   

     ! Include linear f terms
     if(numvar.ge.2 .and. i3d.eq.1) then
        ! b2vector = r15 * bf(n)
        
        ! make a larger vector that can be multiplied by a numvar=3 matrix
        phip = 0.
        do l=1,numnodes
           call entdofs(1, l, 0, ibegin, iendplusone)
           call entdofs(vecsize, l, 0, ibeginnv, iendplusonenv)
           
           phip(ibeginnv  :ibeginnv+5) = bf(ibegin:ibegin+5)
        enddo
        call matrixvectormult(o2matrix_sm,phip,b2vector)
        vtemp = vtemp + b2vector
     endif
  
     ! Insert boundary conditions
     if(calc_matrices.eq.1) then
        call boundary_mag(s2matrix_sm, vtemp)
        call finalizematrix(s2matrix_sm)
     else 
        call boundary_mag(0, vtemp)
     endif
  
     ! solve linear system...LU decomposition done first time
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  
     call solve(s2matrix_sm, vtemp, jer)
     
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
        vtemp = (2.*vtemp - phiold)/3.
     endif
     phiold = phi
     phi = vtemp
  else
     phiold = phi
  endif

  
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

  integer, intent(in) :: calc_matrices
  integer :: l, numnodes, jer
  integer :: ibegin, iendplusone, ibeginnv, iendplusonenv
  
  real :: tstart, tend

  if(myrank.eq.0 .and. iprint.ge.1) print *, "Solving matrix equation."
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  
  ! vtemp = d1matrix_sm * phi(n)
  call matrixvectormult(d1matrix_sm,phi,vtemp)
  
  vtemp = vtemp + q4

  ! Include linear f terms
  if(numvar.ge.2 .and. i3d.eq.1) then
     ! b2vector = r15 * bf(n)
     
     ! make a larger vector that can be multiplied by a vecsize matrix
     phip = 0.
     do l=1,numnodes
        call entdofs(1, l, 0, ibegin, iendplusone)
        call entdofs(vecsize, l, 0, ibeginnv, iendplusonenv)
        
        phip(ibeginnv  :ibeginnv+5) = bf(ibegin:ibegin+5)
     enddo
     call matrixvectormult(o1matrix_sm,phip,b2vector)
     vtemp = vtemp + b2vector
  endif
  
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_mvm = t_mvm + tend - tstart
  endif
  
  ! Insert boundary conditions
  if(calc_matrices.eq.1) then
     call boundary_mag(s1matrix_sm, vtemp)
     call boundary_vel(s1matrix_sm, vtemp)
     if(idens.eq.1) call boundary_den(s1matrix_sm, vtemp)
     if(ipres.eq.1) call boundary_pres(s1matrix_sm, vtemp)
     call finalizematrix(s1matrix_sm)
  else 
     call boundary_mag(0, vtemp)
     call boundary_vel(0, vtemp)
     if(idens.eq.1) call boundary_den(0, vtemp)
     if(ipres.eq.1) call boundary_pres(0, vtemp)
  endif
  
  ! solve linear system...LU decomposition done first time
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call solve(s1matrix_sm, vtemp, jer)
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
     vtemp = (2.*vtemp - phiold)/3.
  endif
  phiold = phi
  phi = vtemp
end subroutine unsplit_step
