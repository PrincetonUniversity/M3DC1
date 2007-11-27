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

  use basic
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
  call zerosuperluarray(matrix, numvar1_numbering)
  do itri=1,numelms

     call area_to_local(79,                                            &
          alpha_79,beta_79,gamma_79,area_weight_79,                    &
          atri(itri), btri(itri), ctri(itri),                          &
          si_79, eta_79, weight_79)

     call calcpos(itri, si_79, eta_79, 79, x_79, z_79)
     if(itor.eq.1) then
        r_79 = x_79
     else
        r_79 = 1.
     endif
     ri_79 = 1./r_79

     do i=1,18
        call eval_ops(gtri(:,i,itri), si_79, eta_79, ttri(itri), &
             ri_79, 79, g79(:,:,i))
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

  call finalizearray(matrix)

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

     call calcpos(itri, si_79, eta_79, 79, x_79, z_79)
     if(itor.eq.1) then
        r_79 = x_79
     else
        r_79 = 1.
     endif
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
  double precision :: coords(3)

  real :: factor

  resistivity = 0.
  kappa = 0.
  sigma = 0.
  visc = 0.
  tempvar = 0.

  temp79c = 0.

  def_fields = 0.

  if(idens.eq.1)  def_fields = def_fields + FIELD_N
  if(numvar.ge.3) def_fields = def_fields + FIELD_PE + FIELD_P

  ! Calculate RHS
  call numfac(numelms)
  do itri=1,numelms

     call define_fields_79(itri, def_fields)   
        
     ! resistivity
     ! ~~~~~~~~~~~
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
     ! ~~~~~~~~~~~~~~~~~~~~
     if((kappa0.eq.0. .and. kappah.eq.0) .or. numvar.lt.3) then
        temp79b = 0.
     else
        ! kappa = p/T**(3/2) = sqrt(n**3/p)
        if(idens.eq.0) then
!!$           temp79c = (eta0*temp79a/2.)**2 * &
!!$                (9.*(pet79(:,OP_DZ)**2+pet79(:,OP_DR)**2)/pet79(:,OP_1)**2)

!!$           temp79b = kappa0/sqrt(pt79(:,OP_1)) &
!!$                + kappah/(1.+sqrt(temp79c))

           temp79b = kappa0/sqrt(pt79(:,OP_1)) &
                + kappah/(1.+pt79(:,OP_LP)**2)
        else
           temp79c = (eta0*temp79a/2.)**2 * &
                    ((nt79(:,OP_DZ)**2 + nt79(:,OP_DR)**2)/ nt79(:,OP_1)**2 &
                +9.*(pet79(:,OP_DZ)**2 +pet79(:,OP_DR)**2)/pet79(:,OP_1)**2 &
                - 6.*(nt79(:,OP_DZ)*pet79(:,OP_DZ) &
                     +nt79(:,OP_DR)*pet79(:,OP_DR)) &
                    /(nt79(:,OP_1 )*pet79(:,OP_1 )))

           temp79b = kappa0*sqrt(nt79(:,OP_1)**3/pt79(:,OP_1)) &
                + kappah/(1.+sqrt(temp79c))
        endif        
     endif
 
     ! density source
     ! ~~~~~~~~~~~~~~
     if(idens.eq.1) then
        temp79c = 0.
        if(ipellet.eq.1) then
           temp79c = temp79c + ri_79*pellet_rate/(2.*pi*pellet_var**2) & 
                *exp(-((x_79 - pellet_x)**2 + (z_79 - pellet_z)**2) &
                      /(2.*pellet_var**2))
        endif
        if(ionization.eq.1) then
           if(numvar.ge.3) then
              temp79d = pt79(:,OP_1)
           else
              temp79d = p0
           endif
           if(idens.eq.1) then
              temp79d = temp79d / nt79(:,OP_1)
           endif

           do i=1, 79
              if(temp79d(i) .gt. ionization_temp) then
                 temp79e(i) = exp(-(temp79d(i) - ionization_temp) &
                                   / ionization_depth)
              else
                 temp79e(i) = 1.
              endif
           enddo

           temp79c = temp79c + ionization_rate * temp79e * &
                exp(-ionization_temp / temp79d)

        endif
     endif

     ! visc
     ! ~~~~
     do i=1,79
        call mask(x_79(i)-xzero, z_79(i)-zzero, factor)
        temp79d(i) = amu_edge*(1.-factor) + 1.
     end do

     ! tempvar
     ! ~~~~~~~
     call interpolate_size_field(itri)
     temp79e = sz79(:,OP_1)

     do i=1,18
        ione = isval1(itri,i)       
        
        resistivity(ione) = resistivity(ione) &
             + etar*int1(g79(:,OP_1,i),weight_79,79) &
             + eta0*int2(g79(:,OP_1,i),temp79a, weight_79,79)
        kappa(ione) = kappa(ione) &
             + kappat*int1(g79(:,OP_1,i),weight_79,79) &
             +        int2(g79(:,OP_1,i),temp79b, weight_79,79)
        sigma(ione) = sigma(ione) &
             + int2(g79(:,OP_1,i),temp79c, weight_79,79)
        visc(ione) = visc(ione) &
             + int2(g79(:,OP_1,i),temp79d, weight_79,79)
        tempvar(ione) = tempvar(ione) &
             + int2(g79(:,OP_1,i),temp79e, weight_79,79)
     end do
  end do

  call solve_newvar(resistivity, NV_NOBOUND)
  call solve_newvar(kappa, NV_NOBOUND)
  call solve_newvar(sigma, NV_NOBOUND)
  call solve_newvar(visc, NV_NOBOUND)
  call solve_newvar(tempvar, NV_NOBOUND)

  visc_c = amuc*visc
  visc = amu*visc

end subroutine define_transport_coefficients


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
