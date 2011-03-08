module scorec_matrix_mod

  use element

  implicit none

  type scorec_matrix
     integer :: imatrix
     integer :: icomplex
     integer :: isize, m, n
     logical :: lhs
  end type scorec_matrix

  integer, parameter :: MAT_SET = 0
  integer, parameter :: MAT_ADD = 1

  interface clear_mat
     module procedure scorec_matrix_clear
  end interface

  interface create_mat
     module procedure scorec_matrix_create
  end interface

  interface destroy_mat
     module procedure scorec_matrix_destroy
  end interface

  interface finalize
     module procedure scorec_matrix_finalize
  end interface

  interface flush
     module procedure scorec_matrix_flush
  end interface

  interface get_global_node_indices
     module procedure scorec_matrix_get_global_node_indices
  end interface

  interface get_element_indices
     module procedure scorec_matrix_get_element_indices
  end interface

  interface get_node_indices
     module procedure scorec_matrix_get_node_indices
  end interface

  interface insert
     module procedure scorec_matrix_insert_real
#ifdef USECOMPLEX
     module procedure scorec_matrix_insert_complex
#endif
  end interface

  interface insert_block
     module procedure scorec_matrix_insert_block
  end interface

  interface insert_global
     module procedure scorec_matrix_insert_global_real
#ifdef USECOMPLEX
     module procedure scorec_matrix_insert_global_complex
#endif
  end interface

  interface matvecmult
     module procedure scorec_matrix_matvecmult
  end interface

  interface newsolve
     module procedure scorec_matrix_solve
  end interface

  interface set_matrix_index
     module procedure scorec_matrix_set_index
  end interface

  interface write_matrix
     module procedure scorec_matrix_write
  end interface

contains

  !====================================================================
  ! set_matrix_index
  ! ~~~~~~~~~~~~~~~~
  ! sets the index of a scorec matrix
  ! this must be called before the matrix is created
  !====================================================================
  subroutine scorec_matrix_set_index(mat, imat)
    implicit none
    type(scorec_matrix) :: mat
    integer, intent(in) :: imat

    mat%imatrix = imat
  end subroutine scorec_matrix_set_index

  !====================================================================
  ! create
  ! ~~~~~~
  ! creates a scorec solve matrix with index imat and size isize
  !====================================================================
  subroutine scorec_matrix_create(mat, m, n, icomplex, lhs)
    implicit none

    type(scorec_matrix) :: mat
    integer, intent(in) :: m, n, icomplex
    logical, intent(in) :: lhs

#include "finclude/petsc.h"

#ifdef PetscDEV
    PetscBool :: flg_petsc
#else
    PetscTruth :: flg_petsc
#endif
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
  end subroutine scorec_matrix_create


  !====================================================================
  ! clear
  ! ~~~~~
  ! zeroes out a scorec solve matrix
  !====================================================================
  subroutine scorec_matrix_clear(mat)
    implicit none

    type(scorec_matrix), intent(inout) :: mat

    call create_mat(mat,mat%m,mat%n,mat%icomplex,mat%lhs)
  end subroutine scorec_matrix_clear


  !====================================================================
  ! destroy
  ! ~~~~~~~
  ! destroys and scorec solve matrix
  !====================================================================
  subroutine scorec_matrix_destroy(mat)
    implicit none
    
    type(scorec_matrix) :: mat
    call clear_mat(mat)
    call deletematrix(mat%imatrix)
  end subroutine scorec_matrix_destroy


  !====================================================================
  ! matvecmult
  ! ~~~~~~~~~~
  ! matrix vector multiply with scorec data structures
  !====================================================================
  subroutine scorec_matrix_matvecmult(mat,vin,vout)

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
     
  end subroutine scorec_matrix_matvecmult


  !====================================================================
  ! insert_real
  ! ~~~~~~~~~~~
  ! inserts (or increments) a matrix element with real value
  !====================================================================
  subroutine scorec_matrix_insert_real(mat,val,i,j,iop)
    implicit none

    type(scorec_matrix), intent(in) :: mat
    real, intent(in) :: val
    integer, intent(in) :: i, j, iop

    if(val.eq.0.) then
       if(i.ne.j) return
    endif
    call insertval(mat%imatrix, val, 0, i, j, iop)
  end subroutine scorec_matrix_insert_real


#ifdef USECOMPLEX
  !====================================================================
  ! insert_complex
  ! ~~~~~~~~~~~~~~
  ! inserts (or increments) a matrix element with complex value
  !====================================================================
  subroutine scorec_matrix_insert_complex(mat,val,i,j,iop)
    implicit none

    type(scorec_matrix), intent(in) :: mat
    complex, intent(in) :: val
    integer, intent(in) :: i, j, iop

    if(val.eq.0.) then
       if(i.ne.j) return
    endif
    call insertval(mat%imatrix, val, 1, i, j, iop)
  end subroutine scorec_matrix_insert_complex
#endif

  !====================================================================
  ! insert_global_real
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~
  ! inserts (or increments) a matrix element with real value
  ! using global indices
  !====================================================================
  subroutine scorec_matrix_insert_global_real(mat,val,i,j,iop)

    implicit none

    type(scorec_matrix), intent(in) :: mat
    real, intent(in) :: val
    integer, intent(in) :: i, j, iop

    call globalinsertval(mat%imatrix, val, 0, i, j, iop)
  end subroutine scorec_matrix_insert_global_real


#ifdef USECOMPLEX
  !====================================================================
  ! insert_global_complex
  ! ~~~~~~~~~~~~~~~~~~~~~
  ! inserts (or increments) a matrix element with complex value
  ! using local indices
  !====================================================================
  subroutine scorec_matrix_insert_global_complex(mat,val,i,j,iop)

    implicit none

    type(scorec_matrix), intent(in) :: mat
    complex, intent(in) :: val
    integer, intent(in) :: i, j, iop

    call globalinsertval(mat%imatrix, val, 1, i, j, iop)
  end subroutine scorec_matrix_insert_global_complex
#endif


  !====================================================================
  ! solve
  ! ~~~~~
  ! linear matrix solve with scorec data structures
  !====================================================================
  subroutine scorec_matrix_solve(mat, v, ierr)
    use vector_mod
    
    implicit none

#include "finclude/petsc.h"
    
    type(scorec_matrix), intent(in) :: mat
    type(vector_type), intent(inout) :: v
    integer, intent(out) :: ierr
#ifdef CJ_MATRIX_DUMP
    integer :: ndof, gndof, i
    real rms, grms
#endif 

#ifdef PetscDEV
    PetscBool :: flg_petsc, flg_solve2, flg_pdslin
#else
    PetscTruth :: flg_petsc, flg_solve2, flg_pdslin
#endif

    call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-ipetsc', flg_petsc ,ierr)
    call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve2', flg_solve2,ierr)
    call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-pdslin', flg_pdslin,ierr)

    if(flg_solve2.eq.PETSC_TRUE) then  ! use pppl petsc
       call solve2(mat%imatrix,v%data,ierr)

    else if(flg_pdslin.eq.PETSC_TRUE) then  ! use pdslin
#ifdef USEHYBRID
       call hybridsolve(mat%imatrix,v%data,ierr)
#else
       print *, 'Error: not compiled with PDSLIN.  Using default solve.'
       call solve(mat%imatrix,v%data,ierr)
#endif
    else  ! use scorec superlu or petsc (-ipetsc)
       call solve(mat%imatrix,v%data,ierr)
    endif

#ifdef CJ_MATRIX_DUMP
    call numdofs(v%isize, ndof) 
    call mpi_allreduce(ndof, gndof, 1, MPI_INTEGER, &
         MPI_SUM, MPI_COMM_WORLD, ierr) 
    rms=0.
    do i=1, ndof
       rms = rms + v%data(i) * v%data(i)
       !print *, "scorec_matrix_solve sol=", mat%imatrix, i, v%data(i)
    enddo
    call mpi_allreduce(rms, grms, 1, MPI_DOUBLE_PRECISION, &
         MPI_SUM, MPI_COMM_WORLD, ierr) 
    write(*,'(a,2i7,2e20.10)') "scorec_matrix_solve sol_rms=", &
                          mat%imatrix, ndof, grms, sqrt(grms)/real(gndof) 
#endif 
  end subroutine scorec_matrix_solve

  !====================================================================
  ! finalize
  ! ~~~~~~~~
  ! finalizes matrix 
  !====================================================================
  subroutine scorec_matrix_finalize(mat)
    implicit none
    type(scorec_matrix) :: mat   

    call finalizematrix(mat%imatrix)
  end subroutine scorec_matrix_finalize


  !====================================================================
  ! flush
  ! ~~~~~
  ! flushes matrix 
  !====================================================================
  subroutine scorec_matrix_flush(mat)
    implicit none
    type(scorec_matrix) :: mat
  end subroutine scorec_matrix_flush

  !======================================================================
  ! get_node_indices
  ! ~~~~~~~~~~~~~~~~
  ! given a matrix mat and node inode, returns:
  ! irow(i,j): the local row index associated with dof j associated 
  !           with field i
  ! icol(i,j): the local column index associated with dof j associated 
  !           with field i
  !======================================================================
  subroutine scorec_matrix_get_node_indices(mat, inode, irow, icol)
    use vector_mod
    implicit none
    type(scorec_matrix), intent(in) :: mat
    integer, intent(in) :: inode
    integer, intent(out), dimension(mat%m,dofs_per_node) :: irow
    integer, intent(out), dimension(mat%n,dofs_per_node) :: icol

    integer, dimension(mat%isize,dofs_per_node) :: ind
    integer :: i

    call get_node_indices(mat%isize,inode,ind)
    do i=1,mat%m
       irow(i,:) = ind(i,:)
    end do
    do i=1,mat%n
       icol(i,:) = ind(i,:)
    end do    
  end subroutine scorec_matrix_get_node_indices

  !======================================================================
  ! get_global_node_indices
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  ! given a matrix mat and node inode, returns:
  ! irow(i,j): the global row index associated with dof j associated 
  !           with field i
  ! icol(i,j): the global column index associated with dof j associated 
  !           with field i
  !======================================================================
  subroutine scorec_matrix_get_global_node_indices(mat, inode, irow, icol)
    use vector_mod
    implicit none
    type(scorec_matrix), intent(in) :: mat
    integer, intent(in) :: inode
    integer, intent(out), dimension(mat%m,dofs_per_node) :: irow
    integer, intent(out), dimension(mat%n,dofs_per_node) :: icol

    integer, dimension(mat%isize,dofs_per_node) :: ind
    integer :: i

    call get_global_node_indices(mat%isize,inode,ind)
    do i=1,mat%m
       irow(i,:) = ind(i,:)
    end do
    do i=1,mat%n
       icol(i,:) = ind(i,:)
    end do
  end subroutine scorec_matrix_get_global_node_indices

  !======================================================================
  ! get_element_indices
  ! ~~~~~~~~~~~~~~~~~~~
  ! given a matrix mat and element itri, returns:
  ! irow(i,j): the local row index associated with dof j associated 
  !           with field i
  ! icol(i,j): the local column index associated with dof j associated 
  !           with field i
  !======================================================================
  subroutine scorec_matrix_get_element_indices(mat, itri, irow, icol)
    use vector_mod
    implicit none
    type(scorec_matrix), intent(in) :: mat
    integer, intent(in) :: itri
    integer, intent(out), dimension(mat%m,dofs_per_element) :: irow
    integer, intent(out), dimension(mat%n,dofs_per_element) :: icol

    integer, dimension(mat%isize,dofs_per_element) :: ind
    integer :: i
    call get_element_indices(mat%isize,itri,ind)
    do i=1,mat%m
       irow(i,:) = ind(i,:)
    end do
    do i=1,mat%n
       icol(i,:) = ind(i,:)
    end do
  end subroutine scorec_matrix_get_element_indices

  subroutine scorec_matrix_insert_block(mat, itri, m, n, val, iop)
    implicit none
    type(matrix_type) :: mat
    integer, intent(in) :: itri, m, n, iop
    vectype, intent(in), dimension(dofs_per_element,dofs_per_element) :: val

    integer, dimension(mat%m,dofs_per_element) :: irow
    integer, dimension(mat%n,dofs_per_element) :: icol
    integer :: i, j

    call get_element_indices(mat, itri, irow, icol)

    do i=1, dofs_per_element
       do j=1, dofs_per_element
          call insert(mat,val(i,j),irow(m,i),icol(n,j),iop)
       end do
    end do
  end subroutine scorec_matrix_insert_block

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
  end subroutine set_row_vals

  subroutine scorec_matrix_write(mat, file)
    type(scorec_matrix), intent(in) :: mat
    character(len=*) :: file

    call writematrixtofile(mat%imatrix, mat%imatrix)
  end subroutine scorec_matrix_write


end module scorec_matrix_mod
