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

  integer, parameter :: MAT_SET = 0
  integer, parameter :: MAT_ADD = 1

  interface clear_mat
     module procedure petsc_matrix_clear
  end interface

  interface create_mat
     module procedure petsc_matrix_create
  end interface

  interface destroy_mat
     module procedure petsc_matrix_destroy
  end interface

  interface flush
     module procedure petsc_matrix_flush
  end interface

  interface finalize
     module procedure petsc_matrix_finalize
  end interface

  interface get_element_indices
     module procedure petsc_matrix_get_element_indices
  end interface

  interface get_global_node_indices
     module procedure petsc_matrix_get_node_indices
  end interface

  interface get_node_indices
     module procedure petsc_matrix_get_node_indices
  end interface

  interface insert
     module procedure petsc_matrix_insert_real
#ifdef USECOMPLEX
     module procedure petsc_matrix_insert_complex
#endif
  end interface

  interface insert_block
     module procedure petsc_matrix_insert_block
  end interface

  interface insert_global
     module procedure petsc_matrix_insert_real
#ifdef USECOMPLEX
     module procedure petsc_matrix_insert_complex
#endif
  end interface

  interface matvecmult
     module procedure petsc_matrix_matvecmult
  end interface

  interface newsolve
     module procedure petsc_matrix_solve
  end interface

  interface set_matrix_index
     module procedure petsc_matrix_set_matrix_index
  end interface

  interface write_matrix
     module procedure petsc_matrix_write
  end interface

contains

  !====================================================================
  ! set_matrix_index
  ! ~~~~~~~~~~~~~~~~
  ! sets the index of a petsc matrix
  ! this must be called before the matrix is created
  !====================================================================
  subroutine petsc_matrix_set_matrix_index(mat, imat)
    implicit none
    type(petsc_matrix) :: mat
    integer, intent(in) :: imat

  end subroutine petsc_matrix_set_matrix_index

  !====================================================================
  ! petsc_matrix_create
  ! ~~~~~~~~~~~~~~~~~~~
  ! creates a petsc matrix with index imat and size isize
  !====================================================================
  subroutine petsc_matrix_create(mat, m, n, icomplex, lhs)
    use mesh_mod
    use vector_mod
    implicit none
#include "finclude/petsc.h"
#ifndef PETSC_31
#include "finclude/petscmat.h"
#endif

    type(petsc_matrix) :: mat
    integer, intent(in) :: m, n, icomplex
    logical, intent(in) :: lhs

    integer :: ierr, local_m, local_n, global_m, global_n

    ! average number of nodes coupled to any given node (including itself)
#ifdef USE3D
    integer, parameter :: neighbors = 7*3
#else
    integer, parameter :: neighbors = 7
#endif

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
         neighbors*dofs_per_node*mat%n, PETSC_NULL_INTEGER, &
         0, PETSC_NULL_INTEGER, ierr)
    call MatSeqAIJSetPreallocation(mat%data, &
         neighbors*dofs_per_node*mat%n, PETSC_NULL_INTEGER, ierr)
    call MatMPIBAIJSetPreallocation(mat%data,dofs_per_node, &
         neighbors*mat%n,PETSC_NULL_INTEGER,&
         0,PETSC_NULL_INTEGER,ierr)
    call MatSeqBAIJSetPreallocation(mat%data,dofs_per_node, &
         neighbors*mat%n,PETSC_NULL_INTEGER,ierr)

!    call MatSetOption(mat%data,MAT_KEEP_ZEROED_ROWS,PETSC_TRUE,ierr)

    if(lhs) then
       call KSPCreate(PETSC_COMM_WORLD,mat%context,ierr)
       call KSPSetOperators(mat%context,mat%data,mat%data, &
            SAME_NONZERO_PATTERN,ierr)
       call KSPSetFromOptions(mat%context,ierr)
    end if

  end subroutine petsc_matrix_create


  !====================================================================
  ! clear
  ! ~~~~~
  ! zeroes out a petsc matrix
  !====================================================================
  subroutine petsc_matrix_clear(mat)
    implicit none
    type(petsc_matrix), intent(inout) :: mat
    integer :: ierr

    call MatZeroEntries(mat%data, ierr)
  end subroutine petsc_matrix_clear


  !====================================================================
  ! destroy
  ! ~~~~~~~
  ! destroys a petsc matrix
  !====================================================================
  subroutine petsc_matrix_destroy(mat)
    implicit none   
    type(petsc_matrix) :: mat
    integer :: ierr

    if(mat%lhs) call KSPDestroy(mat%context,ierr)
    call MatDestroy(mat%data, ierr)
  end subroutine petsc_matrix_destroy


  !====================================================================
  ! matvecmult
  ! ~~~~~~~~~~
  ! matrix vector multiply with petsc data structures
  !====================================================================
  subroutine petsc_matrix_matvecmult(mat,vin,vout)

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
     
  end subroutine petsc_matrix_matvecmult


  !====================================================================
  ! insert_real
  ! ~~~~~~~~~~~
  ! inserts (or increments) a matrix element
  !====================================================================
  subroutine petsc_matrix_insert_real(mat,val,i,j,iop)
    use vector_mod
    use mesh_mod
    implicit none
#include "finclude/petscvec.h"

    type(petsc_matrix), intent(in) :: mat
    real, intent(in) :: val
    integer, intent(in) :: i, j, iop

    PetscScalar :: data(1)
    integer :: im(1), in(1), ierr

    if(val.eq.0) then
       if(i.ne.j) return
    end if

    data(1) = val
    im(1) = i - 1
    in(1) = j- 1
    
    select case(iop)
    case(MAT_ADD)
       call MatSetValues(mat%data,1,im,1,in,data,ADD_VALUES,ierr)
    case(MAT_SET)
       call MatSetValues(mat%data,1,im,1,in,data,INSERT_VALUES,ierr)
    end select
  end subroutine petsc_matrix_insert_real


#ifdef USECOMPLEX
  !====================================================================
  ! insert_complex
  ! ~~~~~~~~~~~~~~
  ! inserts (or increments) a matrix element
  !====================================================================
  subroutine petsc_matrix_insert_complex(mat,val,i,j,iop)
    use vector_mod
    implicit none
#include "finclude/petscvec.h"

    type(petsc_matrix), intent(in) :: mat
    complex, intent(in) :: val
    integer, intent(in) :: i, j, iop

    PetscScalar :: data(1)
    integer :: im(1), in(1), ierr

    if(val.eq.0.) then
       if(i.ne.j) return
    endif

    data(1) = val
    im(1) = i - 1
    in(1) = j - 1

    select case(iop)
    case(MAT_ADD)
       call MatSetValues(mat%data,1,im,1,in,data,ADD_VALUES,ierr)
    case(MAT_SET)
       call MatSetValues(mat%data,1,im,1,in,data,INSERT_VALUES,ierr)
    end select     
  end subroutine petsc_matrix_insert_complex
#endif

  !====================================================================
  ! solve
  ! ~~~~~
  ! linear matrix solve with petsc data structures
  !====================================================================
  subroutine petsc_matrix_solve(mat, v, ierr)
    use vector_mod
    
    implicit none
    type(petsc_matrix), intent(in) :: mat
    type(petsc_vector), intent(inout) :: v
    integer, intent(out) :: ierr

    call finalize(v)
    call KSPSolve(mat%context,v%vec,v%vec,ierr)
!    call finalize(v)

  end subroutine petsc_matrix_solve

  !====================================================================
  ! finalize
  ! ~~~~~~~~
  ! finalizes matrix 
  !====================================================================
  subroutine petsc_matrix_finalize(mat)
    implicit none
#include "finclude/petscmat.h"
    type(petsc_matrix) :: mat
    integer :: ierr

    call MatAssemblyBegin(mat%data, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(mat%data, MAT_FINAL_ASSEMBLY, ierr)
  end subroutine petsc_matrix_finalize

  !====================================================================
  ! flush
  ! ~~~~~
  ! flush matrix operations 
  !====================================================================
  subroutine petsc_matrix_flush(mat)
    implicit none
#include "finclude/petscmat.h"
    type(petsc_matrix) :: mat
    integer :: ierr

    call MatAssemblyBegin(mat%data, MAT_FLUSH_ASSEMBLY, ierr)
    call MatAssemblyEnd(mat%data, MAT_FLUSH_ASSEMBLY, ierr)
  end subroutine petsc_matrix_flush

  !======================================================================
  ! get_element_indices
  ! ~~~~~~~~~~~~~~~~~~~
  ! given a matrix mat and element itri, returns:
  ! irow(i,j): the local row index associated with dof j associated 
  !           with field i
  ! icol(i,j): the local column index associated with dof j associated 
  !           with field i
  !======================================================================
  subroutine petsc_matrix_get_element_indices(mat, itri, irow, icol)
    use vector_mod
    implicit none
    type(petsc_matrix), intent(in) :: mat
    integer, intent(in) :: itri
    integer, intent(out), dimension(mat%m,dofs_per_element) :: irow
    integer, intent(out), dimension(mat%n,dofs_per_element) :: icol

    call get_element_indices(mat%m,itri,irow)
    call get_element_indices(mat%n,itri,icol)
    
  end subroutine petsc_matrix_get_element_indices

  !======================================================================
  ! get_node_indices
  ! ~~~~~~~~~~~~~~~~
  ! given a matrix mat and node inode, returns:
  ! irow(i,j): the local row index associated with dof j associated 
  !           with field i
  ! icol(i,j): the local column index associated with dof j associated 
  !           with field i
  !======================================================================
  subroutine petsc_matrix_get_node_indices(mat, inode, irow, icol)
    use vector_mod
    implicit none
    type(petsc_matrix), intent(in) :: mat
    integer, intent(in) :: inode
    integer, intent(out), dimension(mat%m,dofs_per_node) :: irow
    integer, intent(out), dimension(mat%n,dofs_per_node) :: icol

    call get_node_indices(mat%m,inode,irow)
    call get_node_indices(mat%n,inode,icol)
    
  end subroutine petsc_matrix_get_node_indices

  !======================================================================
  ! insert_block
  ! ~~~~~~~~~~~~
  ! inserts values val into matrix mat, into the
  ! rows associated with field m dofs of itri
  ! cols associated with field n dofs of itri
  !======================================================================
  subroutine petsc_matrix_insert_block(mat, itri, m, n, val, iop)
    use mesh_mod
    implicit none
!#include "finclude/petscvec.h"

    type(matrix_type) :: mat
    integer, intent(in) :: itri, m, n, iop
    vectype, intent(in), dimension(dofs_per_element,dofs_per_element) :: val

    integer, dimension(mat%m,dofs_per_element) :: irow
    integer, dimension(mat%n,dofs_per_element) :: icol
    integer :: i, j

!    integer, dimension(1) :: im, in
!    integer :: ierr

    call get_element_indices(mat, itri, irow, icol)

    do i=1, dofs_per_element
       do j=1, dofs_per_element
          call insert(mat,val(i,j),irow(m,i),icol(n,j),iop)
       end do
    end do
    
!    im(1) = m-1
!    in(1) = n-1
!    select case(iop)
!    case(MAT_ADD)
!       call MatSetValuesBlocked(mat%data,1,im,1,in,val,ADD_VALUES,ierr)
!    case(MAT_SET)
!       call MatSetValuesBlocked(mat%data,1,im,1,in,val,INSERT_VALUES,ierr)
!    end select

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

    integer :: irow_global(1), ierr, i
    integer, dimension(ncols) :: icols_global

    irow_global(1) = irow - 1
    do i=1, ncols
       icols_global(i) = icols(i) - 1
    end do
    call MatSetValues(mat%data, 1, irow_global, ncols, icols_global, &
         vals, INSERT_VALUES, ierr)
    
  end subroutine set_row_vals

  subroutine petsc_matrix_write(m, file)
    implicit none
#include "finclude/petsc.h"    
    type(matrix_type), intent(in) :: m
    character(len=*) :: file

    PetscViewer :: pv
    integer :: ierr
    
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD, file, pv, ierr)
    call MatView(m%data, pv, ierr)
    call PetscViewerDestroy(pv, ierr)

  end subroutine petsc_matrix_write


end module petsc_matrix_mod
