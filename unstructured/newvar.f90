module newvar_mod

  implicit none

  integer, parameter :: NV_NOBOUND = 0
  integer, parameter :: NV_DCBOUND = 1

  integer, parameter :: NV_I_MATRIX = 0
  integer, parameter :: NV_LP_MATRIX = 1
  integer, parameter :: NV_GS_MATRIX = 2
  integer, parameter :: NV_BF_MATRIX = 3

  integer, parameter :: NV_RHS = 0
  integer, parameter :: NV_LHS = 1

contains

!============================================
! create_matrix
! ~~~~~~~~~~~~~
! creates a matrix.
! ibound:
!   NV_NOBOUND: no boundary conditions
!   NV_DCBOUND: dirichlet boundary conditions
! itype: operator (NV_I_MATRIX, etc..)
!============================================
subroutine create_matrix(matrix, ibound, itype, isolve)

  use basic
  use p_data
  use t_data
  use arrays
  use sparse
  use nintegrate_mod

  implicit none

  integer, intent(in) :: matrix, ibound, itype, isolve

  integer :: numelms, itri, i, j, ione, jone, izone, izonedim
  vectype :: temp
  vectype, allocatable :: rhs2(:)

  call numfac(numelms)

  ! populate matrix
  if(isolve.eq.NV_LHS) then
     call zerosuperlumatrix(matrix, icomplex, numvar1_numbering)
  else
     call zeromultiplymatrix(matrix, icomplex, numvar1_numbering)
  end if

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
           selectcase(itype)

           case(NV_I_MATRIX)
              temp = int2(g79(:,OP_1,i),g79(:,OP_1,j),weight_79,79)

           case(NV_LP_MATRIX)
              temp = - &
                   (int2(g79(:,OP_DR,i),g79(:,OP_DR,j),weight_79,79) &
                   +int2(g79(:,OP_DZ,i),g79(:,OP_DZ,j),weight_79,79))

           case(NV_GS_MATRIX)
              temp = - &
                   (int2(g79(:,OP_DR,i),g79(:,OP_DR,j),weight_79,79) &
                   +int2(g79(:,OP_DZ,i),g79(:,OP_DZ,j),weight_79,79))
              if(itor.eq.1) then
                 temp = temp - &
                      2.*int3(ri_79,g79(:,OP_1,i),g79(:,OP_DR,j),weight_79,79)
              endif

           case(NV_BF_MATRIX)
              temp = int3(ri2_79,g79(:,OP_1,i),g79(:,OP_1,j),weight_79,79)
           end select
           call insertval2(matrix, temp, icomplex, ione, jone, 1)
        enddo
     enddo
  enddo

  ! apply boundary conditions
  if(isolve.eq.NV_LHS) then
     if(ibound.eq.NV_DCBOUND) then
        call createvec(rhs2, numvar1_numbering)
        call boundary_dc(matrix, rhs2)
        call deletevec(rhs2)
     end if
  end if

  call finalizematrix(matrix)


end subroutine create_matrix


!=====================================================================
! newvar
! ~~~~~~
! creates rhs and solves matrix equation:
!  A x = B y
!
! ilhs: lhs matrix (A)
! outarray: array to contain the result (x)
! inarray: array containing rhs vector (y)
! iplace: index of rhs vector field in inarray
! numvari: number of fields in inarray
! irhs: rhs matrix(B)
! ibound: boundary conditions to apply
!   NV_NOBOUND: no boundary conditions
!   NV_DCBOUND: dirichlet boundary conditions
!======================================================================
subroutine newvar(ilhsmat,outarray,inarray,iplace,numvari,irhsmat,ibound)

  use basic
  use t_data
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: iplace, ibound, numvari, ilhsmat, irhsmat
  vectype, intent(in) :: inarray(*)   ! size=numvari
  vectype, intent(out) :: outarray(*) ! size=1

  integer :: numnodes, ier, i
  integer :: ibegin, iendplusone, ibegin1, iendplusone1

  vectype, allocatable :: temp(:)

  ! if inarray is bigger than vecsize=1, then 
  ! create vecsize=1 for matrix multiplication
  if(numvari.gt.1) then
     call createvec(temp,1)

     call numnod(numnodes)     
     do i=1,numnodes
        call entdofs(1, i, 0, ibegin1, iendplusone1)
        call entdofs(numvari, i, 0, ibegin, iendplusone)
        ibegin = ibegin + (iplace-1)*6
        temp(ibegin1:ibegin1+5) = inarray(ibegin:ibegin+5)
     enddo
  end if

  call matrixvectormult(irhsmat, temp, outarray)

  if(ibound.eq.NV_DCBOUND) then
     call boundary_dc(0, outarray)
  end if

  call solve(ilhsmat,outarray,ier)

  if(numvari.gt.1) call deletevec(temp)


end subroutine newvar


!=====================================================
! solve_newvar
! ~~~~~~~~~~~~
! Solves equation imatrix*x = rhs
! with ibound boundary conditions applied to rhs.
! rhs is overwritten with result (x) on return.
! Be sure imatrix was also generated using ibound.
!=====================================================
subroutine solve_newvar(rhs, ibound, imatrix)

  use sparse

  implicit none

  integer, intent(in) :: ibound, imatrix
  vectype, dimension(*), intent(inout) :: rhs

  integer :: ier

  call sumsharedppplvecvals(rhs)

  if(ibound.eq.NV_DCBOUND) call boundary_dc(0,rhs)
  call solve(imatrix,rhs,ier)

end subroutine solve_newvar

end module newvar_mod


! define_transport_coefficients
! =============================
subroutine define_transport_coefficients()

  use basic
  use t_data
  use arrays
  use nintegrate_mod
  use newvar_mod
  use sparse

  implicit none

  integer :: i, itri, ibegin, iendplusone, numnodes
  integer :: ione, numelms, def_fields
  double precision :: coords(3)

  real :: factor

  logical, save :: first_time = .true.
  logical :: solve_sigma, solve_kappa, solve_visc, solve_resistivity

  if((linear.eq.1).and.(.not.first_time)) return
  first_time = .false.

  solve_resistivity = eta0.ne.0
  solve_visc = amu_edge.ne.0
  solve_kappa = numvar.ge.3 .and. (kappa0.ne.0 .or. kappah.ne.0)
  solve_sigma = idens.eq.1 .and. (ipellet.eq.1 .or. ionization.eq.1)

  if(implicit_eta.eq.0 .or. ntime.eq.0 .or. iconst_eta.eq.0) then
     resistivity = 0.
  else
     solve_resistivity = .false.
  end if
  kappa = 0.
  sigma = 0.
  visc = 0.
  tempvar = 0.

  temp79c = 0.

  def_fields = FIELD_N + FIELD_PE + FIELD_P

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' defining...'

  ! Calculate RHS
  call numfac(numelms)
  do itri=1,numelms

     call define_fields_79(itri, def_fields)   
        
     ! resistivity
     ! ~~~~~~~~~~~
     if(solve_resistivity) then
        ! resistivity = 1/T**(3/2) = sqrt((n/p)**3)
#ifdef USECOMPLEX
        temp79a = sqrt((nt79(:,OP_1)/(pefac*pet79(:,OP_1)))**3)
#else
        temp79a = sqrt(max((nt79(:,OP_1)/(pefac*pet79(:,OP_1)))**3,0.))
#endif
     endif

     ! thermal conductivity
     ! ~~~~~~~~~~~~~~~~~~~~
     if(solve_kappa) then
        ! kappa = p/T**(3/2) = sqrt(n**3/p)
        temp79b = kappa0*sqrt(nt79(:,OP_1)**3/pt79(:,OP_1))

        if(kappah.ne.0) then
           temp79c = (eta0*temp79a/2.)**2 * &
                ((nt79(:,OP_DZ)**2 + nt79(:,OP_DR)**2)/ nt79(:,OP_1)**2 &
                +9.*(pet79(:,OP_DZ)**2 +pet79(:,OP_DR)**2)/pet79(:,OP_1)**2 &
                - 6.*(nt79(:,OP_DZ)*pet79(:,OP_DZ) &
                     +nt79(:,OP_DR)*pet79(:,OP_DR)) &
                /(nt79(:,OP_1 )*pet79(:,OP_1 )))
           temp79b = temp79b + kappah/(1.+sqrt(temp79c))
        end if
     endif
 
     ! density source
     ! ~~~~~~~~~~~~~~
     if(solve_sigma) then
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
              if(real(temp79d(i)) .gt. ionization_temp) then
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
     if(solve_visc) then
        do i=1,79
           call mask(x_79(i)-xzero, z_79(i)-zzero, factor)
           temp79d(i) = amu_edge*(1.-factor)
        end do
     endif

     ! tempvar
     ! ~~~~~~~
!!$     call interpolate_size_field(itri)
     temp79e = sz79(:,OP_1)

     do i=1,18
        ione = isval1(itri,i)

        if(solve_resistivity) then
           resistivity(ione) = resistivity(ione) &
                + eta0*int2(g79(:,OP_1,i),temp79a, weight_79,79)
        endif

        if(solve_kappa) then
           kappa(ione) = kappa(ione) &
                + int2(g79(:,OP_1,i),temp79b, weight_79,79)           
        endif

        if(solve_sigma) then
           sigma(ione) = sigma(ione) &
                + int2(g79(:,OP_1,i),temp79c, weight_79,79)
        endif

        if(solve_visc) then
           visc(ione) = visc(ione) &
                + int2(g79(:,OP_1,i),temp79d, weight_79,79)
        end if

        tempvar(ione) = tempvar(ione) &
             + int2(g79(:,OP_1,i),temp79e, weight_79,79)
     end do
  end do

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' solving...'

  if(solve_resistivity) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  resistivity'
     call solve_newvar(resistivity, NV_NOBOUND, mass_matrix_lhs)
  end if

  if(solve_kappa) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  kappa'
     call solve_newvar(kappa, NV_NOBOUND, mass_matrix_lhs)
  endif

  if(solve_sigma) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  sigma'
     call solve_newvar(sigma, NV_NOBOUND, mass_matrix_lhs)
  endif

  if(solve_visc) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  viscosity'
     call solve_newvar(visc, NV_NOBOUND, mass_matrix_lhs)
  endif

  if(myrank.eq.0 .and. iprint.ge.1) print *, '  size field'
  call solve_newvar(tempvar, NV_NOBOUND, mass_matrix_lhs)

  ! add in constant components
  call numnod(numnodes)
  do i=1,numnodes
     call entdofs(1,i,0,ibegin,iendplusone)
     if(implicit_eta.eq.0 .or. ntime.eq.0 .or. iconst_eta.eq.0) then
        resistivity(ibegin) = resistivity(ibegin) + etar
     end if
     visc(ibegin) = visc(ibegin) + 1.
     kappa(ibegin) = kappa(ibegin) + kappat
  enddo

  visc_c = amuc*visc
  visc = amu*visc

end subroutine define_transport_coefficients
