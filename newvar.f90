module newvar_mod

  implicit none

  integer, parameter :: NV_NOBOUND = 0
  integer, parameter :: NV_DCBOUND = 1
  integer, parameter :: NV_NMBOUND = 2
  integer, parameter :: NV_SJBOUND = 3
  integer, parameter :: NV_SVBOUND = 4
  integer, parameter :: NV_SCBOUND = 5

  integer, parameter :: NV_I_MATRIX = 0
  integer, parameter :: NV_LP_MATRIX = 1
  integer, parameter :: NV_GS_MATRIX = 2
  integer, parameter :: NV_BF_MATRIX = 3
  integer, parameter :: NV_SJ_MATRIX = 4  ! Current density smoother
  integer, parameter :: NV_SV_MATRIX = 5  ! Vorticity smoother
  integer, parameter :: NV_SC_MATRIX = 6  ! Compression smoother

  integer, parameter :: NV_RHS = 0
  integer, parameter :: NV_LHS = 1

contains


subroutine apply_bc(imatrix, rhs, ibound)
  implicit none
  integer, intent(in) :: imatrix, ibound
  vectype, dimension(*), intent(inout) :: rhs

  select case(ibound)
  case(NV_DCBOUND)
     call boundary_dc(imatrix, rhs)
  case(NV_NMBOUND)
     call boundary_nm(imatrix, rhs)
  case(NV_SJBOUND)
     call boundary_jphi(imatrix, rhs)
  case(NV_SVBOUND)
     call boundary_vor(imatrix, rhs)
  case(NV_SCBOUND)
     call boundary_com(imatrix, rhs)
  end select

end subroutine apply_bc

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

  integer :: numelms, itri, i, ii, iii, j, jj, jjj, in, jn, m, n
  integer :: isize, ibegin, iendplusone, jbegin, jendplusone, ier
  vectype, allocatable :: temp(:,:)
  vectype, allocatable :: rhs2(:)
  real :: hyp

  PetscTruth :: flg_petsc, flg_solve2, flg_solve1

  ! determine the size of the matrix
  select case(itype)
  case(NV_SJ_MATRIX,NV_SV_MATRIX,NV_SC_MATRIX)
     isize = 2
  case default
     isize = 1
  end select

  hyp = hypc
  if(itype.eq.NV_SV_MATRIX .and. ihypamu.eq.1) hyp = hypc*amu
  if(itype.eq.NV_SC_MATRIX .and. ihypamu.eq.1) hyp = hypc*amuc

  ! populate matrix default linear solver superlu cj-april-09-2008
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-ipetsc', flg_petsc,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve2', flg_solve2,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve1', flg_solve1,ier) 

  if(isolve.eq.NV_LHS) then
     if(flg_petsc.eq.PETSC_TRUE) then
        call zeropetscmatrix(matrix, icomplex, isize)
        if(iprint.ge.1) &
        print *, "	newvar_create_matrix zeropetscmatrix", matrix
     else
        call zerosuperlumatrix(matrix, icomplex, isize)
        if(iprint.ge.1) &
        print *, "	newvar_create_matrix zerosuperlumatrix", matrix
     endif
  else
     call zeromultiplymatrix(matrix, icomplex, isize)
  end if

  allocate(temp(isize, isize))
  call numfac(numelms)

  do itri=1,numelms

     call define_triangle_quadrature(itri, 25)
     call define_fields(itri,0,1,0)

     do iii=1,3
     call entdofs(isize,  ist(itri,iii)+1, 0, ibegin, iendplusone)
     do ii=1,6
        i = (iii-1)*6 + ii
        in = ibegin + ii - 1

        do jjj=1,3
        call entdofs(isize,  ist(itri,jjj)+1, 0, jbegin, jendplusone)
        do jj=1,6
           j = (jjj-1)*6 + jj
           jn = jbegin + jj - 1
                 
           temp = 0.

           selectcase(itype)
              
           case(NV_I_MATRIX)
              temp(1,1) = int2(g79(:,OP_1,i),g79(:,OP_1,j))
              
           case(NV_LP_MATRIX)
              temp(1,1) = int2(g79(:,OP_1,i),g79(:,OP_LP,j))
              if(ibound.eq.NV_NMBOUND) then
                 temp(1,1) = temp(1,1) &
                      - regular*int2(g79(:,OP_1,i),g79(:,OP_1,j))
              endif
              
           case(NV_GS_MATRIX)              
              temp(1,1) = int2(g79(:,OP_1,i),g79(:,OP_GS,j))
              
           case(NV_BF_MATRIX)
              temp(1,1) = int3(r2_79,g79(:,OP_1,i),g79(:,OP_LP,j))  &
                   - regular*int3(r2_79,g79(:,OP_1,i),g79(:,OP_1,j))

           case(NV_SJ_MATRIX)
              if(isolve.eq.NV_LHS) then
                 temp(1,1) =  int2(g79(:,OP_1,i),g79(:,OP_1,j))
                 temp(1,2) = -int2(g79(:,OP_1,i),g79(:,OP_GS,j))
                 temp(2,1) = dt*hypf*thimpsm* &
                      int2(g79(:,OP_GS,i),g79(:,OP_GS,j))
                 temp(2,2) = -temp(1,2)
              else
                 temp(2,2) = int2(g79(:,OP_1,i),g79(:,OP_1,j)) &
                      -dt*hypf*(1.-thimpsm)* &
                      int2(g79(:,OP_GS,i),g79(:,OP_GS,j))
              end if

           case(NV_SV_MATRIX)
              if(isolve.eq.NV_LHS) then
                 temp(1,1) =  int2(g79(:,OP_1,i),g79(:,OP_1,j))
                 temp(1,2) = -int2(g79(:,OP_1,i),g79(:,OP_GS,j))
                 temp(2,1) = dt*hyp*thimpsm* &
                      int2(g79(:,OP_GS,i),g79(:,OP_GS,j))
                 temp(2,2) = -temp(1,2)
                 if(inoslip_pol.eq.0) temp(2,2) = temp(2,2) - regular*temp(1,1)
              else
                 temp(2,2) = int2(g79(:,OP_1,i),g79(:,OP_1,j)) &
                      -dt*hypf*(1.-thimpsm)* &
                      int2(g79(:,OP_GS,i),g79(:,OP_GS,j))
              end if

           case(NV_SC_MATRIX)
              if(isolve.eq.NV_LHS) then
                 temp(1,1) =  int2(g79(:,OP_1,i),g79(:,OP_1,j))
                 temp(1,2) = -int2(g79(:,OP_1,i),g79(:,OP_LP,j))
                 temp(2,1) = dt*hyp*thimpsm* &
                      int2(g79(:,OP_LP,i),g79(:,OP_LP,j))
                 temp(2,2) = -temp(1,2)
                 if(inoslip_pol.eq.0) temp(2,2) = temp(2,2) - regular*temp(1,1)
              else
                 temp(2,2) = int2(g79(:,OP_1,i),g79(:,OP_1,j)) &
                      -dt*hypf*(1.-thimpsm)* &
                      int2(g79(:,OP_LP,i),g79(:,OP_LP,j))
              end if


           end select
           do m=1,isize
              do n=1,isize
                 call insval(matrix, temp(m,n), icomplex, &
                      in+6*(m-1), jn+6*(n-1), 1)
              end do
           end do
        end do
        end do
     end do
     end do
  end do

  deallocate(temp)

  ! apply boundary conditions
  if(isolve.eq.NV_LHS) then
     call createvec(rhs2, isize)
     call apply_bc(matrix, rhs2, ibound)
     call deletevec(rhs2)
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
! irhs: rhs matrix(B)
! isize: size of matrix equation
! ibound: boundary conditions to apply
!   NV_NOBOUND: no boundary conditions
!   NV_DCBOUND: Dirichlet boundary conditions
!======================================================================
subroutine newvar(ilhsmat,outarray,irhsmat,inarray,ibound)

  use basic
  use arrays

  implicit none
#include "finclude/petsc.h" 

  integer, intent(in) :: ibound, ilhsmat, irhsmat
  vectype, intent(in) :: inarray(*)
  vectype, intent(out) :: outarray(*)

  integer :: ier
  PetscTruth :: flg_petsc, flg_solve2, flg_solve1
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-ipetsc', flg_petsc,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve2', flg_solve2,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve1', flg_solve1,ier) 

  call matrixvectormult(irhsmat, inarray, outarray)

  call apply_bc(0,outarray,ibound)

  if(flg_petsc.eq.PETSC_TRUE .and. flg_solve1.eq.PETSC_TRUE) then 
     call solve1(ilhsmat,outarray,ier)
  else
     call solve(ilhsmat,outarray,ier)
  endif
  if(ier.ne.0) then
     if(myrank.eq.0) print *, 'Error in newvar solve'
     call safestop(10)
  endif
  
end subroutine newvar


subroutine newvar1(ilhsmat,outarray,inarray,iplace,numvari,irhsmat,ibound)

  use arrays

  implicit none

  integer, intent(in) :: iplace, ibound, numvari, ilhsmat, irhsmat
  vectype, intent(in) :: inarray(*)   ! size=numvari
  vectype, intent(out) :: outarray(*) ! size=1

  vectype, allocatable :: temp(:)

  ! if inarray is bigger than vecsize=1, then 
  ! create vecsize=1 for matrix multiplication
  if(numvari.gt.1) then
     call createvec(temp,1)
     call copyvec(inarray,iplace,numvari,temp,1,1)
     call newvar(ilhsmat,outarray,irhsmat,temp,ibound)
     call deletevec(temp)
  else
     call newvar(ilhsmat,outarray,irhsmat,inarray,ibound)
  end if
  
end subroutine newvar1


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

  call apply_bc(0,rhs,ibound)

  if(flg_petsc.eq.PETSC_TRUE .and. flg_solve1.eq.PETSC_TRUE) then 
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
