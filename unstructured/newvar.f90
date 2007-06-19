module newvar_mod

  implicit none

  integer, parameter :: NV_NOBOUND = 0
  integer, parameter :: NV_DCBOUND = 1

  integer, parameter :: NV_LP = 0
  integer, parameter :: NV_GS = 1


contains


! create_newvar_matrix
! ====================
subroutine create_newvar_matrix(matrix, ibound)

  use p_data
  use t_data
  use arrays
  use sparse
  use nintegrate_mod

  implicit none

  integer, intent(in) :: matrix, ibound

  integer :: numelms, itri, i, j, ione, jone, izone, izonedim
  real :: temp
  real, allocatable :: rhs(:)

  call numfac(numelms)

  ! populate matrix
  call zeroarray4solve(matrix, numvar1_numbering)
  do itri=1,numelms

     call area_to_local(79,                                            &
          alpha_79,beta_79,gamma_79,area_weight_79,                    &
          atri(itri), btri(itri), ctri(itri),                          &
          si_79, eta_79, weight_79)

     call calcr(itri, si_79, eta_79, 79, r_79)
     ri_79 = 1./r_79

     do i=1,18
        call eval_ops(gtri(:,i,itri), si_79, eta_79, ttri(itri), ri_79, 79, g79(:,:,i))
     end do

     if(ijacobian.eq.1) weight_79 = weight_79 * r_79

     do j=1,18
        jone = isval1(itri,j)
        do i=1,18
           ione = isval1(itri,i)
           temp = int2(g79(:,OP_1,i),g79(:,OP_1,j),weight_79,79)
           call insertval(matrix, temp, ione, jone, 1)
        enddo
     enddo
  enddo

  ! apply boundary conditions
  if(ibound.eq.NV_DCBOUND) then
     call createvec(rhs, 1)
     call boundary_dc(matrix, rhs)
     call deletevec(rhs)
  end if

  call finalizearray4solve(matrix)

end subroutine create_newvar_matrix


! newvar_d2
! =========
subroutine newvar_d2(inarray,outarray,itype,ibound,gs)

  use basic
  use t_data
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itype, ibound
  real, intent(in) :: inarray(*) ! length using numvard ordering
  real, intent(out) :: outarray(*) ! length using numvar=1 ordering
  integer, intent(in) :: gs ! NV_GS for grad-shafranov operator, NV_LP for laplacian

  integer :: ndof, numelms, itri, i, j, ione, j1

  real :: sum

  call numdofs(1, ndof)
  do i=1, ndof
     outarray(i) = 0.
  end do

  ! Calculate RHS
  call numfac(numelms)
  do itri=1,numelms

     ! calculate the local sampling points and weights for numerical integration
     call area_to_local(79,                                            &
          alpha_79,beta_79,gamma_79,area_weight_79,                    &
          atri(itri), btri(itri), ctri(itri),                          &
          si_79, eta_79, weight_79)

     call calcr(itri, si_79, eta_79, 79, r_79)
     ri_79 = 1./r_79

     if(ijacobian.eq.1) weight_79 = weight_79*r_79

     do i=1,18
        call eval_ops(gtri(:,i,itri), si_79, eta_79, ttri(itri), ri_79, 79, g79(:,:,i))
     end do

     do i=1,18
        ione = isval1(itri,i)
        sum = 0.
        do j=1,18
           j1 = isvaln(itri,j)

           sum = sum - inarray(j1 + 6*(itype-1)) * &
                (int2(g79(:,OP_DR,i),g79(:,OP_DR,j),weight_79,79) &
                +int2(g79(:,OP_DZ,i),g79(:,OP_DZ,j),weight_79,79))
           if(itor.eq.1 .and. gs.eq.NV_GS) then
              sum = sum - inarray(j1 + 6*(itype-1)) * &
                   2.*int3(ri_79,g79(:,OP_1,i),g79(:,OP_DR,j),weight_79,79)
           endif
        end do
        outarray(ione) = outarray(ione) + sum
     end do
  enddo

  ! solve linear equation
  call solve_newvar(outarray, ibound)

end subroutine newvar_d2


! define_transport_coefficients
! =============================
subroutine define_transport_coefficients()

  use basic
  use t_data
  use arrays
  use nintegrate_mod

  integer :: i, itri
  integer :: ione, numelms, def_fields

  resistivity = 0.
  kappa = 0.

  def_fields = 0.

  if(idens.eq.1)  def_fields = def_fields + FIELD_N
  if(numvar.ge.3) def_fields = def_fields + FIELD_PE + FIELD_P

  ! Calculate RHS
  call numfac(numelms)
  do itri=1,numelms

     call define_fields_79(itri, def_fields)
        
     ! resistivity
     if(eta0.eq.0.) then
        temp79a = 0.
     else       
        ! resistivity = 1/T**(3/2) = sqrt((n/p)**3)
        if(numvar.ge.3) then         
           if(idens.eq.0) then
              temp79a = sqrt((1./(pefac*pet79(:,OP_1)))**3)
           else
              temp79a = sqrt((nt79(:,OP_1)/(pefac*pet79(:,OP_1)))**3)
           endif
        else
           if(idens.eq.0) then
              temp79a = sqrt((1./(p0-pi0))**3)
           else
              temp79a = sqrt((nt79(:,OP_1)/(p0-pi0))**3)
           endif
        endif
     endif

     ! thermal conductivity
     if(kappa0.eq.0. .or. numvar.lt.3) then
        temp79b = 0.
     else
        ! kappa = p/T**(3/2) = sqrt(n**3/p)
        if(idens.eq.0) then
           temp79b = sqrt(1./pt79(:,OP_1))
        else
           temp79b = sqrt(nt79(:,OP_1)**3/pt79(:,OP_1))
        endif        
     endif
        
     do i=1,18
        ione = isval1(itri,i)       
        
        resistivity(ione) = resistivity(ione) &
             + etar*int1(g79(:,OP_1,i),weight_79,79) &
             + eta0*int2(g79(:,OP_1,i),temp79a, weight_79,79)
        kappa(ione) = kappa(ione) &
             + kappat*int1(g79(:,OP_1,i),weight_79,79) &
             + kappa0*int2(g79(:,OP_1,i),temp79b, weight_79,79)
     end do
  end do

  call solve_newvar(resistivity, NV_NOBOUND)
  call solve_newvar(kappa, NV_NOBOUND)

end subroutine define_transport_coefficients



!!$! newvar_eta
!!$! ==========
!!$subroutine newvar_eta()
!!$
!!$  use basic
!!$  use t_data
!!$  use arrays
!!$  use nintegrate_mod
!!$
!!$  use gradshafranov
!!$
!!$  implicit none
!!$
!!$  real :: minpe, ajlim, xguess, zguess
!!$  integer :: i, itri
!!$  integer :: ione, numelms, def_fields
!!$!  integer :: numnodes, ibegin, iendplusone, ibegin1, iendplusone1
!!$
!!$  resistivity = 0.
!!$
!!$  def_fields = 0.
!!$  if(idens.eq.1)  def_fields = def_fields + FIELD_N
!!$  if(numvar.ge.3) def_fields = def_fields + FIELD_PE
!!$
!!$  if(itor.eq.1 .and. itaylor.eq.1 .and. numvar.lt.3 .and. eta0.ne.0) then  
!!$     itri = 0.
!!$     call evaluate(xlim-xzero,zlim-zzero,psilim,ajlim,phi,1,numvar,itri)
!!$
!!$     xguess = xmag - xzero
!!$     zguess = zmag - zzero
!!$     if(linear.eq.1 .or. eqsubtract.eq.1) then
!!$        call magaxis(xguess,zguess,phi+phi0,numvar)
!!$     else
!!$        call magaxis(xguess,zguess,phi,numvar)
!!$     endif
!!$     xmag = xguess + xzero
!!$     zmag = zguess + zzero
!!$
!!$     def_fields = def_fields + FIELD_PSI
!!$  endif
!!$
!!$  ! Calculate RHS
!!$  call numfac(numelms)
!!$  do itri=1,numelms
!!$
!!$     if(eta0.eq.0) then
!!$        temp79a = 0.
!!$     else 
!!$
!!$        call define_fields_79(itri, def_fields)
!!$
!!$        ! for the grad-shafranov simulation with numvar < 3,
!!$        ! calculate the pressure assuming that p(psi) = p0(psi)
!!$        if(numvar.lt.3 .and. itor.eq.1 .and. itaylor.eq.1) then
!!$           temp79c = (pst79(:,OP_1) - psimin)/(psilim - psimin)
!!$           
!!$           do i=1,79
!!$              if(temp79c(i).lt.0) then
!!$                 pet79(i,OP_1) = p0-pi0*ipres
!!$              else if(temp79c(i).gt.1) then
!!$                 pet79(i,OP_1) = pedge*(p0-pi0*ipres)/p0
!!$              else
!!$                 pet79(i,OP_1) = pedge*(p0-pi0*ipres)/p0 + &
!!$                      (p0-pi0*ipres)* &
!!$                      (1.+p1*temp79c(i)+p2*temp79c(i)**2 &
!!$                      -(20. + 10.*p1 + 4.*p2)*temp79c(i)**3 &
!!$                      +(45. + 20.*p1 + 6.*p2)*temp79c(i)**4 &
!!$                      -(36. + 15.*p1 + 4.*p2)*temp79c(i)**5 &
!!$                      +(10. + 4.*p1 + p2)*temp79c(i)**6)
!!$              endif
!!$           end do
!!$        endif
!!$        
!!$        if(itor.eq.1) then
!!$           if(idens.eq.0) then
!!$              temp79a = sqrt((1./(pefac*pet79(:,OP_1)))**3)
!!$           else
!!$              temp79a = sqrt((nt79(:,OP_1)/(pefac*pet79(:,OP_1)))**3)
!!$           endif
!!$        else
!!$           if(numvar.ge.3) then         
!!$              if(idens.eq.0) then
!!$                 temp79a = sqrt((1./(pefac*pet79(:,OP_1)))**3)
!!$              else
!!$                 temp79a = sqrt((nt79(:,OP_1)/(pefac*pet79(:,OP_1)))**3)
!!$              endif
!!$           else
!!$              if(idens.eq.0) then
!!$                 temp79a = sqrt((1./(p0-pi0))**3)
!!$              else
!!$                 temp79a = sqrt((nt79(:,OP_1)/(p0-pi0))**3)
!!$              endif
!!$           endif
!!$        endif
!!$     endif
!!$        
!!$     do i=1,18
!!$        ione = isval1(itri,i)       
!!$        
!!$        resistivity(ione) = resistivity(ione) &
!!$             + etar*int1(g79(:,OP_1,i),weight_79,79) &
!!$             + eta0*int2(g79(:,OP_1,i),temp79a, weight_79,79)
!!$     end do
!!$  enddo
!!$
!!$  ! solve linear equation
!!$  call solve_newvar(resistivity, NV_NOBOUND)
!!$end subroutine newvar_eta


! solve_newvar
! ============
subroutine solve_newvar(rhs, ibound)

  use sparse

  implicit none

  integer, intent(in) :: ibound
  real, dimension(*), intent(inout) :: rhs

  integer :: i, ier

  call sumshareddofs(rhs)

  if(ibound.eq.NV_DCBOUND) then
     call boundary_dc(0,rhs)
     call solve(s6matrix_sm,rhs,ier)
  else
     call solve(s3matrix_sm,rhs,ier)
  endif

end subroutine solve_newvar

end module newvar_mod
