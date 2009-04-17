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
     if(flg_petsc) then
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

     call define_triangle_quadrature(itri, 25)
     call define_fields(itri,0,1)

     do j=1,18
        jone = isval1(itri,j)
        do i=1,18
           ione = isval1(itri,i)
           selectcase(itype)

           case(NV_I_MATRIX)
              temp = int2(g79(:,OP_1,i),g79(:,OP_1,j))

           case(NV_LP_MATRIX)
              temp = int2(g79(:,OP_1,i),g79(:,OP_LP,j))
              if(ibound.eq.NV_NMBOUND) then
                 temp = temp - regular*int2(g79(:,OP_1,i),g79(:,OP_1,j))
              endif

           case(NV_GS_MATRIX)              
              temp = int2(g79(:,OP_1,i),g79(:,OP_GS,j))

           case(NV_BF_MATRIX)
              temp = int3(r2_79,g79(:,OP_1,i),g79(:,OP_LP,j)) &
                   - regular*int3(r2_79,g79(:,OP_1,i),g79(:,OP_1,j))
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

  if(flg_petsc .and. flg_solve1) then 
     call solve1(ilhsmat,outarray,ier)
  else
     call solve(ilhsmat,outarray,ier)
  endif
  if(ier.ne.0) then
     if(myrank.eq.0) print *, 'Error in newvar solve'
     call safestop(10)
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

  use basic
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
  if(flg_petsc .and. flg_solve1) then 
     call solve1(imatrix,rhs,ier)
  else
     call solve(imatrix,rhs,ier)
  endif

  if(ier .ne. 0) then
     if(myrank.eq.0) print *, "Error in solve_newvar solve."
     call safestop(4)
  endif

end subroutine solve_newvar

end module newvar_mod
