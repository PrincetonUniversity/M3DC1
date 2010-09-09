module petsc_vector_mod

  implicit none

#include "finclude/petscvecdef.h"

  integer, parameter :: VEC_SET = 0
  integer,  parameter :: VEC_ADD = 1

  type :: petsc_vector
     Vec :: vec
     integer :: isize
     integer :: insert_mode
     Vec :: vec_local
     PetscScalar, pointer :: data(:)
     logical :: is_finalized
  end type petsc_vector

  
  interface create_vector
     module procedure create_petsc_vec
  end interface

  interface destroy_vector
     module procedure destroy_petsc_vec
  end interface

  interface assignment(=)
     module procedure copy_petsc_vec 
     module procedure const_petsc_vec_real
  end interface

  interface add
     module procedure add_petsc_vectors
  end interface

  interface mult
     module procedure multiply_petsc_vector_by_real
#ifdef USECOMPLEX
     module procedure multiply_petsc_vector_by_complex
#endif
  end interface

  interface insert
     module procedure insert_real_petsc_vec
#ifdef USECOMPLEX
     module procedure insert_complex_petsc_vec
#endif
  end interface

  interface finalize
     module procedure finalize_petsc_vector
  end interface

  interface sum_shared
     module procedure sum_shared_petsc_vector
  end interface

  interface node_index
     module procedure node_index_vector
  end interface

  interface global_node_index
     module procedure node_index_vector
  end interface

  interface get_dof_indices
     module procedure get_dof_indices_petsc
  end interface

  interface get_global_dof_indices
     module procedure get_dof_indices_petsc
  end interface

  interface get_node_index
     module procedure get_node_index_petsc
  end interface

  interface vector_insert_block
     module procedure petsc_vector_insert_block
  end interface

  interface is_nan
     module procedure petsc_is_nan
  end interface

  interface set_vector_node_data
     module procedure set_petsc_vector_node_data_real
#ifdef USECOMPLEX
     module procedure set_petsc_vector_node_data_complex
#endif
  end interface
  interface get_vector_node_data
     module procedure get_petsc_vector_node_data_real
#ifdef USECOMPLEX
     module procedure get_petsc_vector_node_data_complex
#endif
  end interface

contains

  !========================================================
  ! dof_index
  ! ~~~~~~~~~
  ! Get the dof index where ifirst is the first index associated with the node
  ! iplace is the field place, and idof is the local dof number
  !========================================================
  integer function dof_index(ifirst, iplace, idof)
    use element
    integer, intent(in) :: ifirst, iplace, idof
    dof_index = ifirst + (iplace-1)*dofs_per_node + idof-1
  end function dof_index


  !======================================================================
  ! create_petsc_vec
  ! ~~~~~~~~~~~~~~~~
  ! creates a vector of size n
  !======================================================================
  subroutine create_petsc_vec(v,n)
    use mesh_mod
    implicit none
#include "finclude/petsc.h"

    type(petsc_vector), intent(inout) :: v
    integer, intent(in) :: n
    integer :: i, j, k, l, dofs, ierr, nghost
    integer, pointer :: pghost(:)
    integer, allocatable :: vghost(:)

    call get_ghost_nodes(pghost, nghost)

    allocate(vghost(nghost*dofs_per_node*n))
    l = 1
    do i=1, nghost
       do j=1, n
          do k=1, dofs_per_node
             vghost(l) = (pghost(i) - 1)*dofs_per_node*n &
                  + (j - 1)*dofs_per_node + k                  
             l = l + 1
          end do
       end do
    end do
    vghost = vghost - 1     ! PETSc uses 0-based global index

    nghost = nghost*dofs_per_node*n
    dofs = local_nodes()*dofs_per_node*n - nghost

    call VecCreateGhost(PETSC_COMM_WORLD,dofs,PETSC_DECIDE, &
         nghost,vghost,v%vec,ierr)
    deallocate(vghost)

    v%isize = n
    v%is_finalized = .true.
  end subroutine create_petsc_vec

  subroutine destroy_petsc_vec(v)
    implicit none
    type(petsc_vector), intent(inout) :: v
    integer :: ierr

    call VecDestroy(v%vec, ierr)
  end subroutine destroy_petsc_vec


  !======================================================================
  ! copy_petsc_vec
  ! ~~~~~~~~~~~~~~
  ! copy data from vin to vout
  !======================================================================
  subroutine copy_petsc_vec(vout,vin)
    use mesh_mod
    implicit none
#include "finclude/petscvec.h"

    type(petsc_vector), intent(inout) :: vout    
    type(petsc_vector), intent(in) :: vin
    
    integer :: numnodes, dofs, i, j, ierr
    integer, allocatable :: iin(:), iout(:)
    PetscScalar, allocatable :: vals(:)

    if(vin%isize .ne. vout%isize) then
       numnodes = owned_nodes()
       dofs = min(vin%isize, vout%isize)*dofs_per_node
       allocate(iin(dofs), iout(dofs), vals(dofs))
       
       do i=1, numnodes
          iin(1) = node_index(vin, i) - 1
          iout(1) = node_index(vout, i) - 1
          do j=2, dofs
             iin(j) = iin(1) + j - 1
             iout(j) = iout(1) + j - 1
          end do
          call VecGetValues(vin%vec, dofs, iin, vals, ierr)
          call VecSetValues(vout%vec, dofs, iout, vals, INSERT_VALUES, ierr)
       end do
       deallocate(iin, iout, vals)
    else
       call VecCopy(vin%vec,vout%vec,ierr)
    endif

    call finalize(vout)
  end subroutine copy_petsc_vec


  !======================================================================
  ! const_petsc_vec_real
  ! ~~~~~~~~~~~~~~~~~~~~
  ! set all elements of v to s
  !======================================================================
  subroutine const_petsc_vec_real(v,s)
    implicit none
#include "finclude/petscvec.h"

    type(petsc_vector), intent(inout) :: v
    real, intent(in) :: s
    PetscScalar :: val
    integer :: ierr

    val = s
    call VecSet(v%vec, val, ierr)
    call finalize(v)
  end subroutine const_petsc_vec_real



  subroutine add_petsc_vectors(vout,vin)
    use mesh_mod
    implicit none
#include "finclude/petscvec.h"
    type(petsc_vector), intent(inout) :: vout    
    type(petsc_vector), intent(in) :: vin
    
    integer :: numnodes, dofs, i, j, ierr
    integer, allocatable :: iin(:), iout(:)
    PetscScalar :: alpha
    PetscScalar, allocatable :: vals(:)

    if(vin%isize .ne. vout%isize) then
       dofs = min(vin%isize, vout%isize)*dofs_per_node
       allocate(iin(dofs), iout(dofs), vals(dofs))

       numnodes = owned_nodes()
       do i=1, numnodes
          iin(1) = node_index(vin, i) - 1
          iout(1) = node_index(vout, i) - 1
          do j=2, dofs
             iin(j) = iin(1) + j - 1
             iout(j) = iout(1) + j - 1
          end do
          call VecGetValues(vin%vec, dofs, iin, vals, ierr)
          call VecSetValues(vout%vec, dofs, iout, vals, ADD_VALUES, ierr)
       end do
       deallocate(iin, iout, vals)
    else
       alpha = 1.
       call VecAXPY(vout%vec, alpha, vin%vec, ierr)
    endif
    call finalize(vout)
  end subroutine add_petsc_vectors


  subroutine multiply_petsc_vector_by_real(v, s)
    implicit none
#include "finclude/petscvec.h"
    type(petsc_vector), intent(inout) :: v
    real, intent(in) :: s
    integer :: ierr
    PetscScalar :: val
    val = s
    call VecScale(v%vec, s, ierr)
    call finalize(v)
end subroutine multiply_petsc_vector_by_real

#ifdef USECOMPLEX
  subroutine multiply_petsc_vector_by_complex(v, s)
    implicit none
    type(petsc_vector), intent(inout) :: v
    complex, intent(in) :: s
    PetscScalar :: val
    val = s
    call VecScale(v%vec, s, ierr)
    call finalize(v)
  end subroutine multiply_petsc_vector_by_complex
#endif



  !======================================================================
  ! insert_petsc_vec
  ! ~~~~~~~~~~~~~~~~
  ! set element i of v to s
  !======================================================================
  subroutine insert_real_petsc_vec(v, i, s, iop)
    implicit none
#include "finclude/petscvec.h"
    
    type(petsc_vector), intent(inout) :: v
    real, intent(in) :: s
    integer, intent(in) :: i, iop

    integer :: i1, ierr, ind(1)
    PetscScalar :: val(1)

    i1 = 1
    ind(1) = i - 1       ! PETSc uses zero-based arrays
    val(1) = s

    select case(iop)
    case(VEC_SET)
       call VecSetValues(v%vec,i1,ind,val,INSERT_VALUES,ierr)
    case(VEC_ADD)
       call VecSetValues(v%vec,i1,ind,val,ADD_VALUES,ierr)
    case default
       print *, 'Error: invalid vector operation'
    end select
    v%is_finalized = .false.
  end subroutine insert_real_petsc_vec

#ifdef USECOMPLEX
  subroutine insert_complex_petsc_vec(v, i, s, iop)
    implicit none
#include "finclude/petscvec.h"
    
    type(petsc_vector), intent(inout) :: v
    complex, intent(in) :: s
    integer, intent(in) :: i, iop

    integer :: i1, ierr, ind(1)
    PetscScalar :: val(1)

    i1 = 1
    ind(1) = i - 1       ! PETSc uses zero-based arrays
    val(1) = s

    select case(iop)
    case(VEC_SET)
       call VecSetValues(v%vec,i1,ind,val,INSERT_VALUES,ierr)
    case(VEC_ADD)
       call VecSetValues(v%vec,i1,ind,val,ADD_VALUES,ierr)
    case default
       print *, 'Error: invalid vector operation'
    end select
    v%is_finalized = .false.
  end subroutine insert_complex_petsc_vec
#endif


  subroutine sum_shared_petsc_vector(v)
    implicit none
#include "finclude/petscvec.h"

    type(petsc_vector), intent(inout) :: v
    integer :: ierr

!    if(v%is_finalized) return

    call VecAssemblyBegin(v%vec,ierr)
    call VecAssemblyEnd(v%vec,ierr)

    call VecGhostUpdateBegin(v%vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecGhostUpdateEnd(v%vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
    v%is_finalized = .true.
  end subroutine sum_shared_petsc_vector

  subroutine finalize_petsc_vector(v)
    implicit none
#include "finclude/petscvec.h"

    type(petsc_vector), intent(inout) :: v
    integer :: ierr

    if(v%is_finalized) print *, 'multiple finalizes'

    call VecAssemblyBegin(v%vec,ierr)
    call VecAssemblyEnd(v%vec,ierr)

    call VecGhostUpdateBegin(v%vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
    call VecGhostUpdateEnd(v%vec, INSERT_VALUES, SCATTER_FORWARD, ierr)
    v%is_finalized = .true.
  end subroutine finalize_petsc_vector


  !========================================================
  ! get_dof_indices
  ! ~~~~~~~~~~~~~~~~
  ! returns the indices of node inode for each field
  ! of a vector of size isize
  !========================================================
  subroutine get_dof_indices_petsc(isize, inode, ind)
    use mesh_mod
    implicit none
    integer, intent(in) :: isize, inode
    integer, intent(out), dimension(isize,dofs_per_node) :: ind

    integer :: ibegin, i, j
    ibegin = (global_node_id(inode) - 1)*isize*dofs_per_node + 1

    do i=1, isize
       do j=1, dofs_per_node
          ind(i,j) = ibegin + (i-1)*dofs_per_node + j - 1
       end do
    end do

  end subroutine get_dof_indices_petsc

  !========================================================
  ! get_node_index
  ! ~~~~~~~~~~~~~~
  ! returns the index of node inode in field iplace
  ! of a vector of size isize
  !========================================================
  subroutine get_node_index_petsc(inode, iplace, isize, ind)
    use element
    use mesh_mod
    implicit none
    integer, intent(in) :: iplace, isize, inode
    integer, intent(out) :: ind

    ind = ((global_node_id(inode) - 1)*isize + (iplace - 1))*dofs_per_node + 1
  end subroutine get_node_index_petsc


  subroutine petsc_vector_insert_block(v, itri, m, val, iop)
    use element

    implicit none
    type(petsc_vector) :: v
    integer, intent(in) :: itri, m, iop
    vectype, intent(in), dimension(dofs_per_element) :: val

    integer, dimension(v%isize, dofs_per_element) :: irow
    integer :: i

    call get_basis_indices(v%isize, itri, irow)

    do i=1,dofs_per_element
       call insert(v, irow(m,i), val(i), iop)
    end do
  end subroutine petsc_vector_insert_block

  logical function petsc_is_nan(v)
    implicit none
#include "finclude/petsc.h"
    type(petsc_vector), intent(in) :: v
    integer :: ierr
    PetscScalar :: y(1)
    PetscInt :: ni, ix(1)

    ni = 1
    ix(1) = global_dof_id(v%isize, 1) - 1

    call VecGetValues(v%vec, ni, ix, y, ierr)

    petsc_is_nan = y(1).ne.y(1)
  end function petsc_is_nan

  subroutine get_petsc_vector_node_data_real(v, iplace, inode, data)
    use element
    use mesh_mod
    implicit none
#include "finclude/petsc.h"

    type(petsc_vector), intent(in) :: v
    integer, intent(in) :: inode, iplace
    real, intent(out), dimension(dofs_per_node) :: data
    
    PetscScalar, dimension(dofs_per_node) :: vals
    integer :: ind(dofs_per_node), i, ierr
    Vec :: vl

    if(.not.v%is_finalized) then
       if(inode.gt.owned_nodes()) then
          print *, 'WARNING: vec not finalized!', inode, owned_nodes()
       endif
    endif

    ind(1) = (inode-1)*v%isize*dofs_per_node + (iplace-1)*dofs_per_node + 1
    do i=2, dofs_per_node
       ind(i) = ind(1) + i - 1
    end do
    ind = ind - 1
    
    call VecGhostGetLocalForm(v%vec, vl, ierr)
    call VecGetValues(vl, dofs_per_node, ind, vals, ierr)
    data = vals
    call VecGhostRestoreLocalForm(v%vec, vl, ierr)

  end subroutine get_petsc_vector_node_data_real


#ifdef USECOMPLEX
  subroutine get_petsc_vector_node_data_complex(v, iplace, inode, data)
    use element
    implicit none
#include "finclude/petsc.h"
    type(petsc_vector), intent(in) :: v
    integer, intent(in) :: inode, iplace
    complex, intent(out), dimension(dofs_per_node) :: data
    
    PetscScalar, dimension(dofs_per_node) :: vals
    integer :: ind(dofs_per_node), i, ierr

    Vec :: vl

    if(.not.v%is_finalized) then
       if(inode.gt.owned_nodes()) then
          print *, 'Warning: vec not finalized!', inode, owned_nodes()
       endif
    endif

    ind(1) = (inode-1)*v%isize*dofs_per_node + (iplace-1)*dofs_per_node + 1
    do i=2, dofs_per_node
       ind(i) = ind(1) + i - 1
    end do
    ind = ind - 1
    
    call VecGhostGetLocalForm(v%vec, vl, ierr)
    call VecGetValues(vl, dofs_per_node, ind, vals, ierr)
    data = vals
    call VecGhostRestoreLocalForm(v%vec, vl, ierr)
  end subroutine get_petsc_vector_node_data_complex
#endif

  subroutine set_petsc_vector_node_data_real(v, iplace, inode, data)
    use element
    implicit none
#include "finclude/petscvec.h"
    type(petsc_vector), intent(inout) :: v
    integer, intent(in) :: inode, iplace
    real, intent(in), dimension(dofs_per_node) :: data

    integer :: ind(dofs_per_node), i, ierr
    PetscScalar, dimension(dofs_per_node) :: vals

    vals = data
    call get_node_index(inode, iplace, v%isize, ind(1))    
    do i=2, dofs_per_node
       ind(i) = ind(1) + i - 1
    end do
    ind = ind - 1

    call VecSetValues(v%vec, dofs_per_node, ind, vals, INSERT_VALUES, ierr)
    v%is_finalized = .false.
  end subroutine set_petsc_vector_node_data_real

#ifdef USECOMPLEX
  subroutine set_petsc_vector_node_data_complex(v, iplace, inode, data)
    use element
    implicit none
#include "finclude/petscvec.h"
    type(petsc_vector), intent(inout) :: v
    integer, intent(in) :: inode, iplace
    complex, intent(in), dimension(dofs_per_node) :: data

    integer :: ind(dofs_per_node), i, ierr
    PetscScalar, dimension(dofs_per_node) :: vals

    vals = data
    call get_node_index(inode, iplace, v%isize, ind(1))
    do i=2, dofs_per_node
       ind(i) = ind(1) + i - 1
    end do
    ind = ind - 1

    call VecSetValues(v%vec, dofs_per_node, ind, vals, INSERT_VALUES, ierr)
    v%is_finalized = .false.
  end subroutine set_petsc_vector_node_data_complex
#endif

!============================================================================

  !======================================================================
  ! node_index
  ! ~~~~~~~~~~
  !======================================================================
  integer function node_index_vector(v, inode, ifield)
    implicit none
    
    type(vector_type), intent(in) :: v
    integer, intent(in) :: inode
    integer, intent(in), optional :: ifield
    integer :: iplace

    if(present(ifield)) then
       iplace = ifield
    else
       iplace = 1
    endif

    call get_node_index(inode, iplace, v%isize, node_index_vector)
  end function node_index_vector

  integer function global_dof_id(isize, idof_local)
    use element
    use mesh_mod
    implicit none
    integer, intent(in) :: isize, idof_local
    integer :: itemp, inode

    itemp = isize*dofs_per_node*owned_nodes()

    if(idof_local.le.itemp) then
       global_dof_id = (global_node_id(1)-1)*isize*dofs_per_node + idof_local
    else
       inode = (idof_local - itemp - 1) / (isize*dofs_per_node) + 1
       itemp = idof_local - itemp - (inode-1)*(isize*dofs_per_node)
       print *, 'inode, itemp', owned_nodes(), idof_local, inode, itemp
       global_dof_id = (global_node_id(inode)-1)*isize*dofs_per_node + itemp
    endif
  end function global_dof_id

  subroutine write_vector(v, file)
    implicit none
#include "finclude/petsc.h"    
    type(vector_type), intent(in) :: v
    character(len=*) :: file

    PetscViewer :: pv
    integer :: ierr
    
    call PetscViewerASCIIOpen(PETSC_COMM_WORLD, file, pv, ierr)
    call VecView(v%vec, pv, ierr)
    call PetscViewerDestroy(pv, ierr)
  end subroutine write_vector


!========================================================================
  !===================================================================
  ! get_basis_indices
  ! ~~~~~~~~~~~~~~~~~
  ! Return a list of local indices ind
  ! for each dof associated with element itri
  ! in a vector of size isize
  !===================================================================
  subroutine get_basis_indices(isize, itri, ind)
    use mesh_mod

    implicit none
    integer, intent(in) :: isize, itri
    integer, intent(out), dimension(isize,dofs_per_element) :: ind
    
    integer, dimension(isize, dofs_per_node) :: temp
    integer, dimension(nodes_per_element) :: inode
    integer :: i, iii

    call get_element_nodes(itri, inode)

    i = 1
    do iii=1,nodes_per_element
       call get_dof_indices(isize, inode(iii), temp)
       ind(:,i:i+dofs_per_node-1) = temp
       i = i + dofs_per_node
    end do
       
  end subroutine get_basis_indices
    
end module petsc_vector_mod
