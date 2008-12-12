module newvar_mod

  implicit none

  integer, parameter :: NV_NOBOUND = 0
  integer, parameter :: NV_DCBOUND = 1
  integer, parameter :: NV_NMBOUND = 2

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
!   NV_DCBOUND: Dirichlet boundary conditions
!   NV_NMBOUND: Neumann boundary conditions
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

#include "finclude/petsc.h"

  integer, intent(in) :: matrix, ibound, itype, isolve

  integer :: numelms, itri, i, j, ione, jone
  vectype :: temp
  vectype, allocatable :: rhs2(:)

  integer :: ier
  PetscTruth :: flg_petsc, flg_solve2, flg_solve1

  call numfac(numelms)

  ! populate matrix default linear solver superlu cj-april-09-2008
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-ipetsc', flg_petsc,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve2', flg_solve2,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve1', flg_solve1,ier) 

  if(isolve.eq.NV_LHS) then
     if(flg_petsc.eq.PETSC_TRUE) then
        call zeropetscmatrix(matrix, icomplex, numvar1_numbering)
        if(iprint.ge.1) &
        print *, "	newvar_create_matrix zeropetscmatrix", matrix
     else
        call zerosuperlumatrix(matrix, icomplex, numvar1_numbering)
        if(iprint.ge.1) &
        print *, "	newvar_create_matrix zerosuperlumatrix", matrix
     endif
  else
     call zeromultiplymatrix(matrix, icomplex, numvar1_numbering)
  end if

  do itri=1,numelms

     call define_fields(itri,0,25,1)

     do j=1,18
        jone = isval1(itri,j)
        do i=1,18
           ione = isval1(itri,i)
           selectcase(itype)

           case(NV_I_MATRIX)
              temp = int2(g79(:,OP_1,i),g79(:,OP_1,j))

           case(NV_LP_MATRIX)
              temp = - &
                   (int2(g79(:,OP_DR,i),g79(:,OP_DR,j)) &
                   +int2(g79(:,OP_DZ,i),g79(:,OP_DZ,j)))
              if(ibound.eq.NV_NMBOUND) then
                temp = temp - regular*int2(g79(:,OP_1,i),g79(:,OP_1,j))
              endif

           case(NV_GS_MATRIX)
              temp = - &
                   (int2(g79(:,OP_DR,i),g79(:,OP_DR,j)) &
                   +int2(g79(:,OP_DZ,i),g79(:,OP_DZ,j)))
              if(itor.eq.1) then
                 temp = temp - &
                      2.*int3(ri_79,g79(:,OP_1,i),g79(:,OP_DR,j))
              endif

           case(NV_BF_MATRIX)
              temp = int3(ri2_79,g79(:,OP_1,i),g79(:,OP_1,j))
           end select
           call insval(matrix, temp, icomplex, ione, jone, 1)
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
     if(ibound.eq.NV_NMBOUND) then
        call createvec(rhs2, numvar1_numbering)
        call boundary_nm(matrix, rhs2)
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
!   NV_DCBOUND: Dirichlet boundary conditions
!======================================================================
subroutine newvar(ilhsmat,outarray,inarray,iplace,numvari,irhsmat,ibound)

  use basic
  use t_data
  use arrays
  use nintegrate_mod

  implicit none
#include "finclude/petsc.h" 

  integer, intent(in) :: iplace, ibound, numvari, ilhsmat, irhsmat
  vectype, intent(in) :: inarray(*)   ! size=numvari
  vectype, intent(out) :: outarray(*) ! size=1

  vectype, allocatable :: temp(:)

  integer :: ier
  PetscTruth :: flg_petsc, flg_solve2, flg_solve1
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-ipetsc', flg_petsc,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve2', flg_solve2,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve1', flg_solve1,ier) 

  ! if inarray is bigger than vecsize=1, then 
  ! create vecsize=1 for matrix multiplication
  if(numvari.gt.1) then
     call createvec(temp,1)
     call copyvec(inarray,iplace,numvari,temp,1,1)
     call matrixvectormult(irhsmat, temp, outarray)
     call deletevec(temp)
  else
     call matrixvectormult(irhsmat, inarray, outarray)
  end if

  if(ibound.eq.NV_DCBOUND) call boundary_dc(0, outarray)
  if(ibound.eq.NV_NMBOUND) call boundary_nm(0, outarray)

  if(flg_petsc.eq.PETSC_TRUE .and. flg_solve1.eq.PETSC_TRUE) then 
  call solve1(ilhsmat,outarray,ier)
  else
  call solve(ilhsmat,outarray,ier)
  endif
  
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
#include "finclude/petsc.h" 

  integer, intent(in) :: ibound, imatrix
  vectype, dimension(*), intent(inout) :: rhs

  integer :: ier
  PetscTruth :: flg_petsc, flg_solve2, flg_solve1
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-ipetsc', flg_petsc,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve2', flg_solve2,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve1', flg_solve1,ier) 

  call sumsharedppplvecvals(rhs)

  if(ibound.eq.NV_DCBOUND) call boundary_dc(0,rhs)
  if(ibound.eq.NV_NMBOUND) call boundary_nm(0,rhs)
  if(flg_petsc.eq.PETSC_TRUE .and. flg_solve1.eq.PETSC_TRUE) then 
  call solve1(imatrix,rhs,ier)
  else
  call solve(imatrix,rhs,ier)
  endif

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

  real :: factor

  logical, save :: first_time = .true.
  logical :: solve_sigma, solve_kappa, solve_visc, solve_resistivity, &
       solve_visc_e

  if((linear.eq.1).and.(.not.first_time)) return
  first_time = .false.

  solve_resistivity = eta0.ne.0.
  solve_visc = amu_edge.ne.0.
  solve_kappa = numvar.ge.3 .and. (kappa0.ne.0. .or. kappah.ne.0.)
  solve_sigma = idens.eq.1 .and. &
       (ipellet.eq.1 .or. ionization.eq.1 .or. isink.gt.0)
  solve_visc_e = ibootstrap.gt.0

  resistivity = 0.
  kappa = 0.
  sigma = 0.
  visc = 0.
  visc_e = 0.
  tempvar = 0.

  temp79c = 0.

  def_fields = FIELD_N + FIELD_PE + FIELD_P + FIELD_PSI + FIELD_I


  if(myrank.eq.0 .and. iprint.ge.1) print *, ' defining...'

  ! Calculate RHS
  call numfac(numelms)
  do itri=1,numelms

     call define_fields(itri, def_fields, int_pts_aux, 1)
        
     ! resistivity
     ! ~~~~~~~~~~~
     if(solve_resistivity) then
     select case (iresfunc)
     case(0)
        ! resistivity = 1/Te**(3/2) = sqrt((n/pe)**3)
        if(linear.eq.1) then
          temp79a = eta0*sqrt((n079(:,OP_1)/(pefac*pe079(:,OP_1)))**3)
        else
          temp79a = eta0*sqrt((nt79(:,OP_1)/(pefac*pet79(:,OP_1)))**3)
        endif
!
!     added 08/05/08 for stability benchmarking
      case(1)
#ifdef USECOMPLEX
          temp79a = eta0*.5*(1. + tanh((real(ps079(:,OP_1))-(psilim + etaoff*(psilim-psimin)))      &
                                 /(etadelt*(psilim-psimin))))
#else
          temp79a = eta0*.5*(1. + tanh((ps079(:,OP_1)-(psilim + etaoff*(psilim-psimin)))      &
                                 /(etadelt*(psilim-psimin))))
#endif
      end select

     endif

     ! thermal conductivity
     ! ~~~~~~~~~~~~~~~~~~~~
     if(solve_kappa) then
        ! kappa = p/T**(3/2) = sqrt(n**3/p)
        temp79b = kappa0*sqrt(nt79(:,OP_1)**3/pt79(:,OP_1))

        if(kappah.ne.0.) then
           temp79c = ( temp79a/2.)**2 * &
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
     if(isink.ge.1) then
        temp79c = temp79c &
             - nt79(:,OP_1)*ri_79*sink1_rate/(2.*pi*sink1_var**2) & 
             *exp(-((x_79 - sink1_x)**2 + (z_79 - sink1_z)**2) &
             /(2.*sink1_var**2))
     endif
     if(isink.ge.2) then
        temp79c = temp79c &
             - nt79(:,OP_1)*ri_79*sink2_rate/(2.*pi*sink2_var**2) & 
             *exp(-((x_79 - sink2_x)**2 + (z_79 - sink2_z)**2) &
             /(2.*sink2_var**2))
     endif

     ! visc
     ! ~~~~
     if(solve_visc) then
!
!....added 10/18/08  to make viscosity function more like resistivity for iresfunc=1  scj
     select case (iresfunc)
     case(0)
        do i=1,npoints
           call mask(x_79(i)-xzero, z_79(i)-zzero, factor)
           temp79d(i) = amu*amu_edge*(1.-factor)
        end do
     case(1)
#ifdef USECOMPLEX
          temp79d = amu_edge*.5*(1. + tanh((real(ps079(:,OP_1))-(psilim + etaoff*(psilim-psimin)))      &
                                 /(etadelt*(psilim-psimin))))
#else
          temp79d = amu_edge*.5*(1. + tanh((ps079(:,OP_1)-(psilim + etaoff*(psilim-psimin)))      &
                                 /(etadelt*(psilim-psimin))))
#endif
      endselect



     endif

     ! tempvar
     ! ~~~~~~~
     temp79e = sz79(:,OP_1)
!!$     temp79e = ri_79* &
!!$          (bzt79(:,OP_DZ)*pst79(:,OP_DR) - bzt79(:,OP_DR)*pst79(:,OP_DZ))


     ! electron viscosity
     ! ~~~~~~~~~~~~~~~~~~
     temp79f = -amue * r2_79 * &
          (bzt79(:,OP_DZ)*pst79(:,OP_DZ) + bzt79(:,OP_DR)*pst79(:,OP_DR)) &
          / (nt79(:,OP_1)*(pst79(:,OP_DZ)**2 + pst79(:,OP_DR)**2 + 1e-1)**2)

     do i=1,18
        ione = isval1(itri,i)

        if(solve_resistivity) then
           resistivity(ione) = resistivity(ione) &
                +  int2(g79(:,OP_1,i),temp79a)
        endif

        if(solve_kappa) then
           kappa(ione) = kappa(ione) &
                + int2(g79(:,OP_1,i),temp79b)           
        endif

        if(solve_sigma) then
           sigma(ione) = sigma(ione) &
                + int2(g79(:,OP_1,i),temp79c)
        endif

        if(solve_visc) then
           visc(ione) = visc(ione) &
                + int2(g79(:,OP_1,i),temp79d)
        endif

        if(solve_visc_e) then
           visc_e(ione) = visc_e(ione) &
                + int2(g79(:,OP_1,i),temp79f)
        endif

        tempvar(ione) = tempvar(ione) &
             + int2(g79(:,OP_1,i),temp79e)
     end do
  end do

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' solving...'

  if(solve_resistivity) then
     if(myrank.eq.0 .and. iprint.ge.1) then
        print *, '  resistivity'
        print *, "psilim,psimin,etaoff,etadelt,eta0,etar", psilim,psimin,etaoff,etadelt,eta0,etar
     endif
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

  if(solve_visc_e) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, '  electron viscosity'
     call solve_newvar(visc_e, NV_NOBOUND, mass_matrix_lhs)
  endif

  if(myrank.eq.0 .and. iprint.ge.1) print *, '  size field'
  call solve_newvar(tempvar, NV_NOBOUND, mass_matrix_lhs)

  visc_c = visc

  ! add in constant components
  call numnod(numnodes)
  do i=1,numnodes
     call entdofs(1,i,0,ibegin,iendplusone)
     resistivity(ibegin) = resistivity(ibegin) + etar
     visc(ibegin) = visc(ibegin) + amu
     visc_c(ibegin) = visc_c(ibegin) + amuc
     kappa(ibegin) = kappa(ibegin) + kappat
  enddo

end subroutine define_transport_coefficients
