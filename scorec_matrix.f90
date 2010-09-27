module scorec_matrix_mod

  use element

  implicit none

  type scorec_matrix
     integer :: imatrix
     integer :: icomplex
     integer :: isize, m, n
     logical :: lhs
  end type scorec_matrix

  interface create_mat
     module procedure create_scorec_matrix
  end interface

  interface clear_mat
     module procedure clear_scorec_matrix
  end interface

  interface destroy_mat
     module procedure destroy_scorec_matrix
  end interface

  interface set_matrix_index
     module procedure set_scorec_index
  end interface

  interface matvecmult
     module procedure matvecmult_scorec
  end interface

  interface newsolve
     module procedure solve_scorec
  end interface

  interface finalize
     module procedure finalize_scorec
  end interface

  interface flush
     module procedure flush_scorec
  end interface

  interface get_dof_indices
     module procedure scorec_matrix_get_dof_indices
  end interface

  interface get_global_dof_indices
     module procedure scorec_matrix_get_global_dof_indices
  end interface

  interface insert
     module procedure insert_real_scorec
#ifdef USECOMPLEX
     module procedure insert_complex_scorec
#endif
  end interface

  interface insert_global
     module procedure insert_global_real_scorec
#ifdef USECOMPLEX
     module procedure insert_global_complex_scorec
#endif
  end interface


  interface insert_block
     module procedure matrix_insert_block
  end interface

  integer, parameter :: MAT_SET = 0
  integer, parameter :: MAT_ADD = 1

contains

  !====================================================================
  ! set_scorec_index
  ! ~~~~~~~~~~~~~~~~
  ! sets the index of a scorec matrix
  ! this must be called before the matrix is created
  !====================================================================
  subroutine set_scorec_index(mat, imat)
    implicit none
    type(scorec_matrix) :: mat
    integer, intent(in) :: imat

    mat%imatrix = imat
  end subroutine set_scorec_index

  !====================================================================
  ! create_scorec_matrix
  ! ~~~~~~~~~~~~~~~~~~~~
  ! creates a scorec solve matrix with index imat and size isize
  !====================================================================
  subroutine create_scorec_matrix(mat, m, n, icomplex, lhs)
    implicit none

    type(scorec_matrix) :: mat
    integer, intent(in) :: m, n, icomplex
    logical, intent(in) :: lhs

#include "finclude/petsc.h"

    PetscTruth :: flg_petsc
    integer :: ierr

    if(mat%imatrix .le. 0) then
       print*, 'Error: scorec matrix index not set!'
       return
    endif

    mat%m = m
    mat%n = n
    mat%isize = max(m, n)
    mat%icomplex = icomplex
    mat%lhs = lhs

    if(lhs) then
       call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-ipetsc',flg_petsc,ierr)

       if(flg_petsc.eq.PETSC_TRUE) then
          call zeropetscmatrix(mat%imatrix, mat%icomplex, mat%isize)
       else
          call zerosuperlumatrix(mat%imatrix, mat%icomplex, mat%isize)
       endif
    else
       call zeromultiplymatrix(mat%imatrix, mat%icomplex, mat%isize)
    endif
  end subroutine create_scorec_matrix


  !====================================================================
  ! clear_scorec_matrix
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~
  ! zeroes out a scorec solve matrix
  !====================================================================
  subroutine clear_scorec_matrix(mat)
    implicit none

    type(scorec_matrix), intent(inout) :: mat

    call create_scorec_matrix(mat,mat%m,mat%n,mat%icomplex,mat%lhs)
  end subroutine clear_scorec_matrix


  !====================================================================
  ! destroy_scorec_matrix
  ! ~~~~~~~~~~~~~~~~~~~~~
  ! destroys and scorec solve matrix
  !====================================================================
  subroutine destroy_scorec_matrix(mat)
    implicit none
    
    type(scorec_matrix) :: mat
    call deletematrix(mat%imatrix)
  end subroutine destroy_scorec_matrix


  !====================================================================
  ! matvecmult_scorec
  ! ~~~~~~~~~~~~~~~~~
  ! matrix vector multiply with scorec data structures
  !====================================================================
  subroutine matvecmult_scorec(mat,vin,vout)

    use vector_mod

    implicit none

    type(scorec_matrix), intent(in) :: mat
    type(vector_type), intent(in), target :: vin
    type(vector_type), intent(inout), target :: vout

    type(vector_type), pointer :: vin_ptr, vout_ptr
    type(vector_type), target :: temp_in, temp_out

    if(mat%m .ne. vout%isize) then
       print *, 'Error: sizes for matvacmult do not conform'
       return
    endif
    if(mat%n .ne. vin%isize) then
       print *, 'Error: sizes for matvacmult do not conform'
       return
    endif

    ! scorec matvecmult only works for square matrices
    if(vin%isize .eq. mat%isize) then
       vin_ptr => vin
    else
       call create_vector(temp_in, mat%isize)
       temp_in = vin
       vin_ptr => temp_in
    endif
    if(vout%isize .eq. mat%isize) then
       vout_ptr => vout
    else
       call create_vector(temp_out, mat%isize)
       vout_ptr => temp_out
    endif

    call matrixvectormult(mat%imatrix,vin_ptr%data,vout_ptr%data)

    nullify(vin_ptr, vout_ptr)

    if(vin%isize .ne. mat%isize) call destroy_vector(temp_in)
    if(vout%isize .ne. mat%isize) then
       vout = temp_out
       call destroy_vector(temp_out)
    endif
     
  end subroutine matvecmult_scorec


  !====================================================================
  ! insert_real_scorec
  ! ~~~~~~~~~~~~~~~~~~
  ! inserts (or increments) a matrix element with real value
  !====================================================================
  subroutine insert_real_scorec(mat,val,i,j,iop)
    implicit none

    type(scorec_matrix), intent(in) :: mat
    real, intent(in) :: val
    integer, intent(in) :: i, j, iop

    call insertval(mat%imatrix, val, 0, i, j, iop)
  end subroutine insert_real_scorec


#ifdef USECOMPLEX
  !====================================================================
  ! insert_complex_scorec
  ! ~~~~~~~~~~~~~~~~~~~~~
  ! inserts (or increments) a matrix element with complex value
  !====================================================================
  subroutine insert_complex_scorec(mat,val,i,j,iop)
    implicit none

    type(scorec_matrix), intent(in) :: mat
    complex, intent(in) :: val
    integer, intent(in) :: i, j, iop

    call insertval(mat%imatrix, val, 1, i, j, iop)
  end subroutine insert_complex_scorec
#endif

  !====================================================================
  ! insert_global_real_scorec
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~
  ! inserts (or increments) a matrix element with real value
  ! using global indices
  !====================================================================
  subroutine insert_global_real_scorec(mat,val,i,j,iop)

    implicit none

    type(scorec_matrix), intent(in) :: mat
    real, intent(in) :: val
    integer, intent(in) :: i, j, iop

    call globalinsertval(mat%imatrix, val, 0, i, j, iop)
  end subroutine insert_global_real_scorec


#ifdef USECOMPLEX
  !====================================================================
  ! insert_global_complex_scorec
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! inserts (or increments) a matrix element with complex value
  ! using local indices
  !====================================================================
  subroutine insert_global_complex_scorec(mat,val,i,j,iop)

    implicit none

    type(scorec_matrix), intent(in) :: mat
    complex, intent(in) :: val
    integer, intent(in) :: i, j, iop

    call globalinsertval(mat%imatrix, val, 1, i, j, iop)
  end subroutine insert_global_complex_scorec
#endif


  !====================================================================
  ! solve_scorec
  ! ~~~~~~~~~~~~
  ! linear matrix solve with scorec data structures
  !====================================================================
  subroutine solve_scorec(mat, v, ierr)
    use vector_mod
    
    implicit none

#include "finclude/petsc.h"
    
    type(scorec_matrix), intent(in) :: mat
    type(vector_type), intent(inout) :: v
    integer, intent(out) :: ierr

    PetscTruth :: flg_petsc, flg_solve1

    call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-ipetsc', flg_petsc ,ierr)
    call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve1', flg_solve1,ierr)

    if(flg_petsc.eq.PETSC_TRUE .and. flg_solve1.eq.PETSC_TRUE) then
       call solve1(mat%imatrix,v%data,ierr)
    else
       call solve(mat%imatrix,v%data,ierr)
    endif
  end subroutine solve_scorec

  !====================================================================
  ! finalize_scorec
  ! ~~~~~~~~~~~~~~~
  ! finalizes matrix 
  !====================================================================
  subroutine finalize_scorec(mat)
    implicit none
    type(scorec_matrix) :: mat   

    call finalizematrix(mat%imatrix)
  end subroutine finalize_scorec


  !====================================================================
  ! flush_scorec
  ! ~~~~~~~~~~~~
  ! flushes matrix 
  !====================================================================
  subroutine flush_scorec(mat)
    implicit none
    type(scorec_matrix) :: mat
  end subroutine flush_scorec

  subroutine scorec_matrix_get_dof_indices(mat, inode, irow, icol)
    use vector_mod
    implicit none
    type(scorec_matrix), intent(in) :: mat
    integer, intent(in) :: inode
    integer, intent(out), dimension(mat%m,dofs_per_node) :: irow
    integer, intent(out), dimension(mat%n,dofs_per_node) :: icol

    integer, dimension(mat%isize,dofs_per_node) :: ind
    integer :: i

    call get_dof_indices(mat%isize,inode,ind)
    do i=1,mat%m
       irow(i,:) = ind(i,:)
    end do
    do i=1,mat%n
       icol(i,:) = ind(i,:)
    end do    
  end subroutine scorec_matrix_get_dof_indices

  subroutine scorec_matrix_get_global_dof_indices(mat, inode, irow, icol)
    use vector_mod
    implicit none
    type(scorec_matrix), intent(in) :: mat
    integer, intent(in) :: inode
    integer, intent(out), dimension(mat%m,dofs_per_node) :: irow
    integer, intent(out), dimension(mat%n,dofs_per_node) :: icol

    integer, dimension(mat%isize,dofs_per_node) :: ind
    integer :: i

    call get_global_dof_indices(mat%isize,inode,ind)
    do i=1,mat%m
       irow(i,:) = ind(i,:)
    end do
    do i=1,mat%n
       icol(i,:) = ind(i,:)
    end do
  end subroutine scorec_matrix_get_global_dof_indices


  subroutine get_indices(mat, itri, irow, icol)
    use vector_mod
    implicit none
    type(scorec_matrix), intent(in) :: mat
    integer, intent(in) :: itri
    integer, intent(out), dimension(mat%m,dofs_per_element) :: irow
    integer, intent(out), dimension(mat%n,dofs_per_element) :: icol

    integer, dimension(mat%isize,dofs_per_element) :: ind
    integer :: i
    call get_basis_indices(mat%isize,itri,ind)
    do i=1,mat%m
       irow(i,:) = ind(i,:)
    end do
    do i=1,mat%n
       icol(i,:) = ind(i,:)
    end do
  end subroutine get_indices

  subroutine matrix_insert_block(mat, itri, m, n, val, iop)
    implicit none
    type(matrix_type) :: mat
    integer, intent(in) :: itri, m, n, iop
    vectype, intent(in), dimension(dofs_per_element,dofs_per_element) :: val

    integer, dimension(mat%m,dofs_per_element) :: irow
    integer, dimension(mat%n,dofs_per_element) :: icol
    integer :: i, j

    call get_indices(mat, itri, irow, icol)

    do i=1, dofs_per_element
       do j=1, dofs_per_element
          call insert(mat,val(i,j),irow(m,i),icol(n,j),iop)
       end do
    end do
  end subroutine matrix_insert_block

  subroutine identity_row(mat, irow)
    implicit none
    type(scorec_matrix) :: mat
    integer, intent(in) :: irow

    call setdiribc(mat%imatrix, irow)
    
  end subroutine identity_row

  subroutine set_row_vals(mat, irow, ncols, icols, vals)
    use vector_mod

    type(scorec_matrix) :: mat
    integer, intent(in) :: irow, ncols
    integer, intent(in), dimension(ncols) :: icols
    vectype, intent(in), dimension(ncols) :: vals

#ifdef USECOMPLEX
    call setgeneralbc(mat%imatrix, irow, ncols, icols, vals, 1)
#else
    call setgeneralbc(mat%imatrix, irow, ncols, icols, vals, 0)
#endif

!!$    integer :: i
!!$
!!$    do i=1, ncols
!!$       call insert(mat,vals(i),irow,icols(i),MAT_ADD)
!!$    end do
  end subroutine set_row_vals

  subroutine write_matrix(mat)
    type(scorec_matrix), intent(in) :: mat
    
    call writematrixtofile(mat%imatrix, mat%imatrix)
  end subroutine write_matrix


end module scorec_matrix_mod
