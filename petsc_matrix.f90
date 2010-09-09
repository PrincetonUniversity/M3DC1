module petsc_matrix_mod

  use element

  implicit none

#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"

  type :: petsc_matrix
     Mat :: data
     KSP :: context
     integer :: m, n
     logical :: lhs
  end type petsc_matrix

  interface create_mat
     module procedure create_petsc_matrix
  end interface

  interface clear_mat
     module procedure clear_petsc_matrix
  end interface

  interface destroy_mat
     module procedure destroy_petsc_matrix
  end interface

  interface set_matrix_index
     module procedure set_petsc_index
  end interface

  interface matvecmult
     module procedure matvecmult_petsc
  end interface

  interface newsolve
     module procedure solve_petsc
  end interface

  interface flush
     module procedure flush_petsc
  end interface

  interface finalize
     module procedure finalize_petsc
  end interface

  interface get_dof_indices
     module procedure petsc_matrix_get_dof_indices
  end interface

  interface get_global_dof_indices
     module procedure petsc_matrix_get_dof_indices
  end interface

  interface insert
     module procedure insert_real_petsc
#ifdef USECOMPLEX
     module procedure insert_complex_petsc
#endif
  end interface

  interface insert_global
     module procedure insert_real_petsc
#ifdef USECOMPLEX
     module procedure insert_complex_petsc
#endif
  end interface

  interface insert_block
     module procedure petsc_matrix_insert_block
  end interface

  integer, parameter :: MAT_SET = 0
  integer, parameter :: MAT_ADD = 1

contains

  !====================================================================
  ! set_petsc_index
  ! ~~~~~~~~~~~~~~~~
  ! sets the index of a petsc matrix
  ! this must be called before the matrix is created
  !====================================================================
  subroutine set_petsc_index(mat, imat)
    implicit none
    type(petsc_matrix) :: mat
    integer, intent(in) :: imat

  end subroutine set_petsc_index

  !====================================================================
  ! create_petsc_matrix
  ! ~~~~~~~~~~~~~~~~~~~
  ! creates a petsc solve matrix with index imat and size isize
  !====================================================================
  subroutine create_petsc_matrix(mat, m, n, icomplex, lhs)
    use mesh_mod
    use vector_mod
    implicit none
#include "finclude/petsc.h"
#include "finclude/petscmat.h"

    type(petsc_matrix) :: mat
    integer, intent(in) :: m, n, icomplex
    logical, intent(in) :: lhs

    integer :: ierr, local_m, local_n, global_m, global_n, i
    integer, allocatable :: gid(:)

    mat%m = m
    mat%n = n
    mat%lhs = lhs

    global_m = m*global_dofs()
    global_n = n*global_dofs()
    local_m = m*owned_dofs()
    local_n = n*owned_dofs()

    call MatCreate(PETSC_COMM_WORLD, mat%data, ierr)
    call MatSetSizes(mat%data, local_m, local_n, global_m, global_n, ierr)
    call MatSetFromOptions(mat%data,ierr)
    call MatMPIAIJSetPreallocation(mat%data, &
         7*dofs_per_node*mat%n, PETSC_NULL_INTEGER, &
         0, PETSC_NULL_INTEGER, ierr)
    call MatSeqAIJSetPreallocation(mat%data, &
         7*dofs_per_node*mat%n, PETSC_NULL_INTEGER, ierr)
    call MatMPIBAIJSetPreallocation(mat%data,dofs_per_node, &
         7*mat%n,PETSC_NULL_INTEGER,&
         0,PETSC_NULL_INTEGER,ierr)
    call MatSeqBAIJSetPreallocation(mat%data,dofs_per_node, &
         7*mat%n,PETSC_NULL_INTEGER,ierr)


!    call MatSetOption(mat%data,MAT_KEEP_ZEROED_ROWS,PETSC_TRUE,ierr)

    if(lhs) then
       call KSPCreate(PETSC_COMM_WORLD,mat%context,ierr)
       call KSPSetOperators(mat%context,mat%data,mat%data, &
            SAME_NONZERO_PATTERN,ierr)
       call KSPSetFromOptions(mat%context,ierr)
    end if

  end subroutine create_petsc_matrix


  !====================================================================
  ! clear_petsc_matrix
  ! ~~~~~~~~~~~~~~~~~~
  ! zeroes out a petsc matrix
  !====================================================================
  subroutine clear_petsc_matrix(mat)
    implicit none
    type(petsc_matrix), intent(inout) :: mat
    integer :: ierr

    call MatZeroEntries(mat%data, ierr)
  end subroutine clear_petsc_matrix


  !====================================================================
  ! destroy_petsc_matrix
  ! ~~~~~~~~~~~~~~~~~~~~
  ! destroys a petsc matrix
  !====================================================================
  subroutine destroy_petsc_matrix(mat)
    implicit none   
    type(petsc_matrix) :: mat
    integer :: ierr

    if(mat%lhs) call KSPDestroy(mat%context,ierr)
    call MatDestroy(mat%data, ierr)
  end subroutine destroy_petsc_matrix


  !====================================================================
  ! matvecmult_petsc
  ! ~~~~~~~~~~~~~~~~
  ! matrix vector multiply with petsc data structures
  !====================================================================
  subroutine matvecmult_petsc(mat,vin,vout)

    use vector_mod

    implicit none

    type(petsc_matrix), intent(in) :: mat
    type(vector_type), intent(in), target :: vin
    type(vector_type), intent(inout), target :: vout

    integer :: ierr

    if(mat%m .ne. vout%isize) then
       print *, 'Error: sizes for matvacmult do not conform'
       return
    endif
    if(mat%n .ne. vin%isize) then
       print *, 'Error: sizes for matvacmult do not conform'
       return
    endif

    call VecAssemblyBegin(vout%vec, ierr)
    call VecAssemblyEnd(vout%vec, ierr)
    call MatMult(mat%data,vin%vec,vout%vec,ierr)
     
  end subroutine matvecmult_petsc


  !====================================================================
  ! insert_real_petsc
  ! ~~~~~~~~~~~~~~~~~
  ! inserts (or increments) a matrix element
  !====================================================================
  subroutine insert_real_petsc(mat,val,i,j,iop)
    use vector_mod
    use mesh_mod
    implicit none
#include "finclude/petscvec.h"

    type(petsc_matrix), intent(in) :: mat
    real, intent(in) :: val
    integer, intent(in) :: i, j, iop

    PetscScalar :: data(1)
    integer :: im(1), in(1), ierr
    data(1) = val
    im(1) = i - 1
    in(1) = j- 1
    
    select case(iop)
    case(MAT_ADD)
       call MatSetValues(mat%data,1,im,1,in,data,ADD_VALUES,ierr)
    case(MAT_SET)
       call MatSetValues(mat%data,1,im,1,in,data,INSERT_VALUES,ierr)
    end select
  end subroutine insert_real_petsc


#ifdef USECOMPLEX
  !====================================================================
  ! insert_real_petsc
  ! ~~~~~~~~~~~~~~~~~
  ! inserts (or increments) a matrix element
  !====================================================================
  subroutine insert_complex_petsc(mat,val,i,j,iop)
    use vector_mod
    implicit none
#include "finclude/petscvec.h"

    type(petsc_matrix), intent(in) :: mat
    complex, intent(in) :: val
    integer, intent(in) :: i, j, iop

    PetscScalar :: data(1)
    integer :: im(1), in(1), ierr
    data(1) = val
    im(1) = i - 1
    in(1) = j - 1

    select case(iop)
    case(MAT_ADD)
       call MatSetValues(mat%data,1,im,1,in,data,ADD_VALUES,ierr)
    case(MAT_SET)
       call MatSetValues(mat%data,1,im,1,in,data,INSERT_VALUES,ierr)
    end select     
  end subroutine insert_complex_petsc
#endif

  !====================================================================
  ! solve_petsc
  ! ~~~~~~~~~~~
  ! linear matrix solve with petsc data structures
  !====================================================================
  subroutine solve_petsc(mat, v, ierr)
    use vector_mod
    
    implicit none
    type(petsc_matrix), intent(in) :: mat
    type(petsc_vector), intent(inout) :: v
    integer, intent(out) :: ierr

!    call finalize(v)
    call KSPSolve(mat%context,v%vec,v%vec,ierr)
!    call finalize(v)

  end subroutine solve_petsc

  !====================================================================
  ! finalize_petsc
  ! ~~~~~~~~~~~~~~
  ! finalizes matrix 
  !====================================================================
  subroutine finalize_petsc(mat)
    implicit none
#include "finclude/petscmat.h"
    type(petsc_matrix) :: mat
    integer :: ierr

    call MatAssemblyBegin(mat%data, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(mat%data, MAT_FINAL_ASSEMBLY, ierr)
  end subroutine finalize_petsc

  !====================================================================
  ! flush_petsc
  ! ~~~~~~~~~~~
  ! flush matrix 
  !====================================================================
  subroutine flush_petsc(mat)
    implicit none
#include "finclude/petscmat.h"
    type(petsc_matrix) :: mat
    integer :: ierr

    call MatAssemblyBegin(mat%data, MAT_FLUSH_ASSEMBLY, ierr)
    call MatAssemblyEnd(mat%data, MAT_FLUSH_ASSEMBLY, ierr)
  end subroutine flush_petsc


  subroutine get_indices(mat, itri, irow, icol)
    use vector_mod
    implicit none
    type(petsc_matrix), intent(in) :: mat
    integer, intent(in) :: itri
    integer, intent(out), dimension(mat%m,dofs_per_element) :: irow
    integer, intent(out), dimension(mat%n,dofs_per_element) :: icol

    call get_basis_indices(mat%m,itri,irow)
    call get_basis_indices(mat%n,itri,icol)
    
  end subroutine get_indices

  subroutine petsc_matrix_get_dof_indices(mat, inode, irow, icol)
    use vector_mod
    implicit none
    type(petsc_matrix), intent(in) :: mat
    integer, intent(in) :: inode
    integer, intent(out), dimension(mat%m,dofs_per_node) :: irow
    integer, intent(out), dimension(mat%n,dofs_per_node) :: icol

    call get_dof_indices(mat%m,inode,irow)
    call get_dof_indices(mat%n,inode,icol)
    
  end subroutine petsc_matrix_get_dof_indices


  subroutine petsc_matrix_insert_block(mat, itri, m, n, val, iop)
    use mesh_mod
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
  end subroutine petsc_matrix_insert_block

  subroutine identity_row(mat, irow)
    use mesh_mod
    use vector_mod
    implicit none
#include "finclude/petscvec.h"

    type(petsc_matrix) :: mat
    integer, intent(in) :: irow

    integer :: irow_global, ierr
    PetscScalar, parameter :: diag = 1.

    irow_global = irow - 1
    call MatSetValue(mat%data, irow_global, irow_global, diag, &
         INSERT_VALUES, ierr)

  end subroutine identity_row

  subroutine set_row_vals(mat, irow, ncols, icols, vals)
    use mesh_mod
    use vector_mod
#include "finclude/petscvec.h"

    type(petsc_matrix) :: mat
    integer, intent(in) :: irow, ncols
    integer, intent(in), dimension(ncols) :: icols
    vectype, intent(in), dimension(ncols) :: vals

    integer :: irow_local(1), irow_global(1), ierr, i
    integer, dimension(ncols) :: icols_global

    irow_global(1) = irow - 1
    do i=1, ncols
       icols_global(i) = icols(i) - 1
    end do
    call MatSetValues(mat%data, 1, irow_global, ncols, icols_global, &
         vals, INSERT_VALUES, ierr)
    
  end subroutine set_row_vals

  subroutine write_matrix(m, file)
    implicit none
#include "finclude/petsc.h"    
    type(matrix_type), intent(in) :: m
    character(len=*) :: file

    Vec :: vl
    PetscViewer :: pv
    integer :: ierr
    
    
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD, file, pv, ierr)
    call MatView(m%data, pv, ierr)
    call PetscViewerDestroy(pv, ierr)

  end subroutine write_matrix


end module petsc_matrix_mod
