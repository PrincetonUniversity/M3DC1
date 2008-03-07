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


  call import_time_advance_vectors

  ! advance time
  ! ~~~~~~~~~~~~
  if(isplitstep.eq.1) then
     call split_step(calc_matrices)
  else
     call unsplit_step(calc_matrices)
  end if

  call export_time_advance_vectors


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
     if(implicit_eta.eq.1) then
        call entdofs(1, l, 0, ibegin, iendplusone)
        eta_v(ibeginnv+eta_off:ibeginnv+eta_off+5) = &
             resistivity(ibegin:ibegin+5)
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
        p_v(ibegin  +  p_off:ibegin  +  p_off+5) =   p1_l
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
     if(implicit_eta.eq.1) then
        call entdofs(1, l, 0, ibegin, iendplusone)
        resistivity(ibegin:ibegin+5) = &
             eta_v(ibeginnv+eta_off:ibeginnv+eta_off+5)
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

  integer, intent(in) :: calc_matrices
  real :: tstart, tend
  integer :: i, l, numnodes, ndofs
  integer :: ibegin, iendplusone, ibeginnv, iendplusonenv
  integer :: jer
  integer, allocatable:: itemp(:)
  vectype, allocatable :: temp(:)

  integer :: istaticold, idensold, ipresold

  call numnod(numnodes)

  ! Store current-time velocity matrices for use in field advance
  veln = vel
  veloldn = velold


  if(istatic.eq.0) then

     ! Advance Velocity
     ! ================
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Velocity"
  
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

     ! d1matrix_sm * vel(n)
     call matrixvectormult(d1matrix_sm,vel,b1_vel)
  
     ! q1matrix_sm * phi(n)
     if(ipres.eq.1 .and. numvar.ge.3) then
        ! replace electron pressure with total pressure
        do l=1,numnodes
           call entdofs(vecsize_p,   l, 0, ibegin,   iendplusone)
           call entdofs(vecsize_phi, l, 0, ibeginnv, iendplusonenv)
           
           phip(ibeginnv   :ibeginnv+11) = phi(ibeginnv:ibeginnv+11)
           phip(ibeginnv+12:ibeginnv+17) = pres(ibegin:ibegin+5)
           if(implicit_eta.eq.1) then
              phip(ibeginnv+18:ibeginnv+23) = phi(ibeginnv+18:ibeginnv+23)
           end if
        enddo
        call matrixvectormult(q1matrix_sm, phip, b1_phi)
     else
        call matrixvectormult(q1matrix_sm, phi , b1_phi)
     endif
   
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_mvm = t_mvm + tend - tstart
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

        call matrixvectormult(r14matrix_sm,phip,b2_phi)
        b1_phi = b1_phi + b2_phi
     endif

     ! o1matrix_sm * bf(n)
     if(numvar.ge.2 .and. i3d.eq.1) then
        ! make a larger vector that can be multiplied by a numvar=3 matrix
        phip = 0.
        do l=1,numnodes
           call entdofs(1, l, 0, ibegin, iendplusone)
           call entdofs(vecsize_phi, l, 0, ibeginnv, iendplusonenv)
           
           phip(ibeginnv  :ibeginnv+5) = bf(ibegin:ibegin+5)
        enddo
        call matrixvectormult(o1matrix_sm,phip,b2_phi)
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
     if(calc_matrices.eq.1) then
        call boundary_vel(s1matrix_sm, b1_vel)
        call finalizematrix(s1matrix_sm)
     else
        call boundary_vel(0, b1_vel)
     endif
  
     ! solve linear system with rhs in vtemp (note LU-decomp done first time)
     if(myrank.eq.0) print *, "solving velocity advance..."
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call solve(s1matrix_sm, b1_vel, jer)
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
     call smooth
  else
     velold = vel
  end if

  
  ! Advance Density
  ! ===============
  if(idens.eq.1) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Density"
     
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

     ! r8matrix_sm * vel(n+1)
     call matrixvectormult(r8matrix_sm,vel,b1_vel)
         
     ! r8matrix_sm * vel(n-1)
     if(integrator.eq.1 .and. ntime.gt.1) then
        call matrixvectormult(r8matrix_sm,veloldn,b2_vel)
        b1_vel = 1.5*b1_vel + 0.5*b2_vel
     endif

     ! q8matrix_sm * vel(n)
     call matrixvectormult(q8matrix_sm,veln,b2_vel)
     b1_vel = b1_vel + b2_vel
     
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_mvm = t_mvm + tend - tstart
     endif
     
     ! temp = d8matrix_sm * phi(n)
     call createvec(temp, numvar1_numbering)
     temp = 0.
     
     call matrixvectormult(d8matrix_sm,den,temp)
     
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
     
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     
     ! r9matrix_sm * vel(n+1)
     call matrixvectormult(r9matrix_sm,vel,b1_vel)
        
     ! r9matrix_sm * vel(n-1)
     if(integrator.eq.1 .and. ntime.gt.1) then
        call matrixvectormult(r9matrix_sm,veloldn,b2_vel)
        b1_vel = 1.5*b1_vel + 0.5*b2_vel
     endif
     
     ! q9matrix_sm * vel(n)
     call matrixvectormult(q9matrix_sm,veln,b2_vel)
     b1_vel = b1_vel + b2_vel

     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_mvm = t_mvm + tend - tstart
     endif
     
     ! temp = d8matrix_sm * pres(n)
     call createvec(temp, numvar1_numbering)
     temp = 0.
     call matrixvectormult(d9matrix_sm,pres,temp)
     

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
  
     ! r2matrix_sm * vel(n+1)
     call matrixvectormult(r2matrix_sm,vel,b1_vel)
   
     ! r2matrix_sm * vel(n-1)
     if(integrator.eq.1 .and. ntime.gt.1) then
        call matrixvectormult(r2matrix_sm,veloldn,b2_vel)
        b1_vel = 1.5*b1_vel + 0.5*b2_vel
     endif

     ! q2matrix_sm * vel(n)
     call matrixvectormult(q2matrix_sm,veln,b2_vel)
     b1_vel = b1_vel + b2_vel


     ! d2matrix_sm * phi(n)
     call matrixvectormult(d2matrix_sm,phi,b1_phi)

    
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_mvm = t_mvm + tend - tstart
     endif
  
     ! Include linear f terms
     if(numvar.ge.2 .and. i3d.eq.1) then
        ! b2vector = r15 * bf(n)
        
        ! make a larger vector that can be multiplied by a numvar=3 matrix
        phip = 0.
        do l=1,numnodes
           call entdofs(1, l, 0, ibegin, iendplusone)
           call entdofs(vecsize_phi, l, 0, ibeginnv, iendplusonenv)
           
           phip(ibeginnv  :ibeginnv+5) = bf(ibegin:ibegin+5)
        enddo
        call matrixvectormult(o2matrix_sm,phip,b2_phi)
        b1_phi = b1_phi + b1_phi
     endif

     ! Construct right-hand side
     call numdofs(vecsize_phi,ndofs)
     allocate(itemp(ndofs)) ! this is used to make sure that we 
     ! don't double count the sum for periodic dofs
     itemp = 1
     do l=1,numnodes
        call entdofs(vecsize_phi, l, 0, ibegin, iendplusone)
        call entdofs(vecsize_vel, l, 0, ibeginnv, iendplusonenv)
        do i=0,iendplusone-ibegin-1
           b1_phi(ibegin+i) = b1_phi(ibegin+i) + itemp(ibegin+i) * &
                (b1_vel(ibeginnv+i) + q4(ibegin+i))

           itemp(ibegin+i) = 0
        enddo
     enddo
     deallocate(itemp)

  
     ! Insert boundary conditions
     if(calc_matrices.eq.1) then
        call boundary_mag(s2matrix_sm, b1_phi)
        call finalizematrix(s2matrix_sm)
     else 
        call boundary_mag(0, b1_phi)
     endif
  
     ! solve linear system...LU decomposition done first time
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  
     call solve(s2matrix_sm, b1_phi, jer)
     
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
        call define_transport_coefficients
        ! revert fields to old values
        phi = b2_phi
        
        ! recalculate field advance matrix
        ! (advanced velocity variables will be used in defining matrix)
        istaticold=istatic
        idensold=idens
        ipresold=ipres
        istatic=1
        idens=0
        ipres=0
        if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
        call ludefall
        if(myrank.eq.0 .and. itimer.eq.1) then
           call second(tend)
           t_ludefall = t_ludefall + tend - tstart
        endif
        istatic=istaticold
        idens=idensold
        ipres=ipresold

        if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Fields Again"

        if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  
        ! r2matrix_sm * vel(n+1)
        call matrixvectormult(r2matrix_sm,vel,b1_vel)
   
        ! r2matrix_sm * vel(n-1)
        if(integrator.eq.1 .and. ntime.gt.1) then
           call matrixvectormult(r2matrix_sm,veloldn,b2_vel)
           b1_vel = 1.5*b1_vel + 0.5*b2_vel
        endif

        ! q2matrix_sm * vel(n)
        call matrixvectormult(q2matrix_sm,veln,b2_vel)
        b1_vel = b1_vel + b2_vel

        ! d2matrix_sm * phi(n)
        call matrixvectormult(d2matrix_sm,phi,b1_phi)
    
        if(myrank.eq.0 .and. itimer.eq.1) then
           call second(tend)
           t_mvm = t_mvm + tend - tstart
        endif
  
        ! Include linear f terms
        if(numvar.ge.2 .and. i3d.eq.1) then
           ! b2vector = r15 * bf(n)
           
           ! make a larger vector that can be multiplied by a numvar=3 matrix
           phip = 0.
           do l=1,numnodes
              call entdofs(1, l, 0, ibegin, iendplusone)
              call entdofs(vecsize_phi, l, 0, ibeginnv, iendplusonenv)
              
              phip(ibeginnv  :ibeginnv+5) = bf(ibegin:ibegin+5)
           enddo
           call matrixvectormult(o2matrix_sm,phip,b2_phi)
           b1_phi = b1_phi + b1_phi
        endif

        ! Construct right-hand side
        call numdofs(vecsize_phi,ndofs)
        allocate(itemp(ndofs)) ! this is used to make sure that we 
        ! don't double count the sum for periodic dofs
        itemp = 1
        do l=1,numnodes
           call entdofs(vecsize_phi, l, 0, ibegin, iendplusone)
           call entdofs(vecsize_vel, l, 0, ibeginnv, iendplusonenv)
           do i=0,iendplusone-ibegin-1
              b1_phi(ibegin+i) = b1_phi(ibegin+i) + itemp(ibegin+i) * &
                   (b1_vel(ibeginnv+i) + q4(ibegin+i))
              
              itemp(ibegin+i) = 0
           enddo
        enddo
        deallocate(itemp)
        
        ! Insert boundary conditions
        if(calc_matrices.eq.1) then
           call boundary_mag(s2matrix_sm, b1_phi)
           call finalizematrix(s2matrix_sm)
        else 
           call boundary_mag(0, b1_phi)
        endif
  
        ! solve linear system...LU decomposition done first time
        if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
        
        call solve(s2matrix_sm, b1_phi, jer)
        
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

  integer, intent(in) :: calc_matrices
  integer :: l, numnodes, jer
  integer :: ibegin, iendplusone, ibeginnv, iendplusonenv
  
  real :: tstart, tend

  if(myrank.eq.0 .and. iprint.ge.1) print *, "Solving matrix equation."
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  
  ! vtemp = d1matrix_sm * phi(n)
  call matrixvectormult(d1matrix_sm,phi,b1_phi)
  
  b1_phi = b1_phi + q4

  ! Include linear f terms
  if(numvar.ge.2 .and. i3d.eq.1) then
     ! b2vector = r15 * bf(n)
     
     ! make a larger vector that can be multiplied by a vecsize matrix
     phip = 0.
     do l=1,numnodes
        call entdofs(1, l, 0, ibegin, iendplusone)
        call entdofs(vecsize_phi, l, 0, ibeginnv, iendplusonenv)
        
        phip(ibeginnv  :ibeginnv+5) = bf(ibegin:ibegin+5)
     enddo
     call matrixvectormult(o1matrix_sm,phip,b2_phi)
     b1_phi = b1_phi + b2_phi
  endif
  
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_mvm = t_mvm + tend - tstart
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
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  call solve(s1matrix_sm, b1_phi, jer)
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
  phiold = phi
  phi = b1_phi
end subroutine unsplit_step
