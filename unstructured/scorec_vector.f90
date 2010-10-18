module scorec_vector_mod
  use element

  type scorec_vector
     vectype, allocatable :: data(:)
     integer :: isize
  end type scorec_vector

  integer, parameter :: VEC_SET = 0
  integer, parameter :: VEC_ADD = 1

  interface assignment(=)
     module procedure scorec_vector_copy
     module procedure scorec_vector_const_real
  end interface

  interface add
     module procedure scorec_vector_add
  end interface

  interface create_vector
     module procedure scorec_vector_create
  end interface

  interface destroy_vector
     module procedure scorec_vector_destroy
  end interface

  interface finalize
     module procedure scorec_vector_finalize
  end interface

  interface get_element_indices
     module procedure scorec_vector_get_element_indices
  end interface

  interface get_global_node_index
     module procedure scorec_vector_get_global_node_index
  end interface

  interface get_global_node_indices
     module procedure scorec_vector_get_global_node_indices
  end interface

  interface get_node_data
     module procedure scorec_vector_get_node_data_real
#ifdef USECOMPLEX
     module procedure scorec_vector_get_node_data_complex
#endif
  end interface

  interface get_node_index
     module procedure scorec_vector_get_node_index
  end interface

  interface get_node_indices
     module procedure scorec_vector_get_node_indices
  end interface

  interface global_node_index
     module procedure scorec_vector_global_node_index
  end interface

  interface insert
     module procedure scorec_vector_insert_real
#ifdef USECOMPLEX
     module procedure scorec_vector_insert_complex
#endif
  end interface

  interface is_nan
     module procedure scorec_vector_is_nan
  end interface

  interface mult
     module procedure scorec_vector_multiply_real
#ifdef USECOMPLEX
     module procedure scorec_vector_multiply_complex
#endif
  end interface

  interface node_index
     module procedure node_index_vector
  end interface

  interface set_node_data
     module procedure scorec_vector_set_node_data_real
#ifdef USECOMPLEX
     module procedure scorec_vector_set_node_data_complex
#endif
  end interface

  interface sum_shared
     module procedure scorec_vector_sum_shared
  end interface

  interface vector_insert_block
     module procedure scorec_vector_insert_block
  end interface

  interface write_vector
     module procedure scorec_vector_write
  end interface

contains

  !========================================================
  ! dof_index
  ! ~~~~~~~~~
  ! Get the dof index where 
  ! ifirst is the first index associated with the node 
  ! iplace is the field place, and 
  ! idof is the local dof number
  !========================================================
  integer function dof_index(ifirst, iplace, idof)
    integer, intent(in) :: ifirst, iplace, idof
    dof_index = ifirst + (iplace-1)*dofs_per_node + idof-1
  end function dof_index
  
  !========================================================
  ! get_node_index
  ! ~~~~~~~~~~~~~~
  ! returns the index of node inode in field iplace
  ! of a vector of size isize
  !========================================================
  subroutine scorec_vector_get_node_index(inode, iplace, isize, ind)

    implicit none
    integer, intent(in) :: iplace, isize, inode
    integer, intent(out) :: ind

    integer :: ibegin, iendplusone
    call entdofs(isize, inode, 0, ibegin, iendplusone)
    ind = ibegin + (iplace-1)*dofs_per_node
  end subroutine scorec_vector_get_node_index

  !========================================================
  ! get_global_node_index_scorec
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! returns the index of node inode in field iplace
  ! of a vector of size isize
  !========================================================
  subroutine scorec_vector_get_global_node_index(inode, iplace, isize, ind)

    implicit none
    integer, intent(in) :: iplace, isize, inode
    integer, intent(out) :: ind

    integer :: ibegin, iendplusone
    call globalentdofs(isize, inode, 0, ibegin, iendplusone)
    ind = ibegin + (iplace-1)*dofs_per_node
  end subroutine scorec_vector_get_global_node_index


  !========================================================
  ! get_node_indices
  ! ~~~~~~~~~~~~~~~~
  ! returns the global indices of node inode for each field
  ! of a vector of size isize
  !========================================================
  subroutine scorec_vector_get_node_indices(isize, inode, ind)

    implicit none
    integer, intent(in) :: isize, inode
    integer, intent(out), dimension(isize,dofs_per_node) :: ind

    integer :: ibegin, iendplusone, i, j
    call entdofs(isize, inode, 0, ibegin, iendplusone)
    do i=1, isize
       do j=1, dofs_per_node
          ind(i,j) = ibegin + (i-1)*dofs_per_node + j - 1
       end do
    end do

  end subroutine scorec_vector_get_node_indices

  !========================================================
  ! get_global_nodes_indices
  ! ~~~~~~~~~~~~~~~~~~~~~~~~
  ! returns the indices of dofs associated with node inode
  ! for each field of a vector of size isize
  !========================================================
  subroutine scorec_vector_get_global_node_indices(isize, inode, ind)

    implicit none
    integer, intent(in) :: isize, inode
    integer, intent(out), dimension(isize,dofs_per_node) :: ind

    integer :: ibegin, iendplusone, i, j
    call globalentdofs(isize, inode, 0, ibegin, iendplusone)
    do i=1, isize
       do j=1, dofs_per_node
          ind(i,j) = ibegin + (i-1)*dofs_per_node + j - 1
       end do
    end do
  end subroutine scorec_vector_get_global_node_indices


  !======================================================================
  ! add_scorec_vectors
  ! ~~~~~~~~~~~~~~~~~~
  ! adds v2 to v1.  result is stored in v1
  !======================================================================
  subroutine scorec_vector_add(v1, v2)
    use mesh_mod

    implicit none

    type(scorec_vector), intent(inout) :: v1
    type(scorec_vector), intent(in) :: v2

    integer :: numnodes, i1, i2, dofs, i

    if(.not.allocated(v1%data)) print *, 'Error: vector not allocated'
    if(.not.allocated(v2%data)) print *, 'Error: vector not allocated'

    if(v1%isize.eq.v2%isize) then
       v1%data = v1%data + v2%data
    else
       numnodes = local_nodes()
       dofs = min(v1%isize, v2%isize)*dofs_per_node

       do i=1, numnodes
          i1 = node_index(v2, i, 1)
          i2 = node_index(v2, i, 1)
          v1%data(i1:i1+dofs-1) = v1%data(i1:i1+dofs-1) + v2%data(i2:i2+dofs-1)
       end do
    endif
  end subroutine scorec_vector_add

  subroutine scorec_vector_multiply_real(v, s)
    implicit none

    type(scorec_vector), intent(inout) :: v
    real :: s

    if(.not.allocated(v%data)) print *, 'Error: vector not allocated'
    v%data = v%data * s
  end subroutine scorec_vector_multiply_real

#ifdef USECOMPLEX
  subroutine scorec_vector_multiply_complex(v, s)
    implicit none

    type(scorec_vector), intent(inout) :: v
    complex :: s

    if(.not.allocated(v%data)) print *, 'Error: vector not allocated'
    v%data = v%data * s
  end subroutine scorec_vector_multiply_complex
#endif


  !======================================================================
  ! const_scorec_vec_real
  ! ~~~~~~~~~~~~~~~~~~~~~
  ! set all elements of v to s
  !======================================================================
  subroutine scorec_vector_const_real(v,s)
    implicit none

    type(scorec_vector), intent(inout) :: v
    real, intent(in) :: s

    if(.not.allocated(v%data)) print *, 'Error: vector not allocated'
    v%data = s
  end subroutine scorec_vector_const_real


  subroutine scorec_vector_get_node_data_real(v, iplace, inode, data)
    implicit none
    type(scorec_vector), intent(in) :: v
    integer, intent(in) :: inode, iplace
    real, intent(out), dimension(dofs_per_node) :: data
    
    integer :: index
    call get_node_index(inode, iplace, v%isize, index)
    data = v%data(index:index+dofs_per_node-1)
  end subroutine scorec_vector_get_node_data_real


#ifdef USECOMPLEX
  subroutine scorec_vector_get_node_data_complex(v, iplace, inode, data)
    implicit none
    type(scorec_vector), intent(in) :: v
    integer, intent(in) :: inode, iplace
    vectype, intent(out), dimension(dofs_per_node) :: data
    
    integer :: index
    call get_node_index(inode, iplace, v%isize, index)
    data = v%data(index:index+dofs_per_node-1)
  end subroutine scorec_vector_get_node_data_complex
#endif

  subroutine scorec_vector_set_node_data_real(v, iplace, inode, data)
    implicit none
    type(scorec_vector), intent(inout) :: v
    integer, intent(in) :: inode, iplace
    real, intent(in), dimension(dofs_per_node) :: data
    
    integer :: index
    call get_node_index(inode, iplace, v%isize, index)
    v%data(index:index+dofs_per_node-1) = data
  end subroutine scorec_vector_set_node_data_real

#ifdef USECOMPLEX
  subroutine scorec_vector_set_node_data_complex(v, iplace, inode, data)
    implicit none
    type(scorec_vector), intent(inout) :: v
    integer, intent(in) :: inode, iplace
    complex, intent(in), dimension(dofs_per_node) :: data
    
    integer :: index
    call get_node_index(inode, iplace, v%isize, index)
    v%data(index:index+dofs_per_node-1) = data
  end subroutine scorec_vector_set_node_data_complex
#endif



  !======================================================================
  ! copy
  ! ~~~~
  ! copy data from vin to vout
  !======================================================================
  subroutine scorec_vector_copy(vout,vin)
    use mesh_mod

    implicit none

    type(scorec_vector), intent(inout) :: vout    
    type(scorec_vector), intent(in) :: vin
    
    integer :: numnodes, dofs, iin, iout, i

    if(vin%isize .ne. vout%isize) then
       numnodes = local_nodes()
       dofs = min(vin%isize, vout%isize)*dofs_per_node

       do i=1, numnodes
          iin = node_index(vin, i, 1)
          iout = node_index(vout, i, 1)
          vout%data(iout:iout+dofs-1) = vin%data(iin:iin+dofs-1)
       end do
    else
       vout%data = vin%data
    endif

  end subroutine scorec_vector_copy

  !======================================================================
  ! insert
  ! ~~~~~~
  ! set element i of v to s
  !======================================================================
  subroutine scorec_vector_insert_real(v, i, s, iop)
    implicit none
    
    type(scorec_vector), intent(inout) :: v
    real, intent(in) :: s
    integer, intent(in) :: i, iop

    if(.not.allocated(v%data)) print *, 'Error: vector not allocated'
    select case(iop)
    case(VEC_SET)
       v%data(i) = s
    case(VEC_ADD)
       v%data(i) = v%data(i) + s
    case default
       print *, 'Error: invalid vector operation'
    end select
  end subroutine scorec_vector_insert_real

#ifdef USECOMPLEX
  subroutine scorec_vector_insert_complex(v, i, s, iop)
    implicit none
    
    type(scorec_vector), intent(inout) :: v
    complex, intent(in) :: s
    integer, intent(in) :: i, iop

    if(.not.allocated(v%data)) print *, 'Error: vector not allocated'
    select case(iop)
    case(VEC_SET)
       v%data(i) = s
    case(VEC_ADD)
       v%data(i) = v%data(i) + s
    case default
       print *, 'Error: invalid vector operation'
    end select
  end subroutine scorec_vector_insert_complex

#endif 

  !======================================================================
  ! create
  ! ~~~~~~
  ! creates a vector of size n
  !======================================================================
  subroutine scorec_vector_create(f,n)
    implicit none

    type(scorec_vector), intent(inout) :: f
    integer, intent(in) :: n
    integer :: ndof

    call numdofs(n, ndof)
    allocate(f%data(ndof))

#ifdef USECOMPLEX
    call createppplvec(f%data, n, 1)
#else
    call createppplvec(f%data, n, 0)
#endif
    f%data = 0.
    f%isize = n
  end subroutine scorec_vector_create

  !======================================================================
  ! destroy
  ! ~~~~~~~
  ! frees the memory allocated with create_field
  !======================================================================
  subroutine scorec_vector_destroy(f)
    implicit none

    type(scorec_vector), intent(inout) :: f
    integer :: i

    call checkppplveccreated(f%data, i)
    if(i .ne. 0) then
       call deleteppplvec(f%data)
       deallocate(f%data)
    endif
  end subroutine scorec_vector_destroy

  subroutine scorec_vector_insert_block(v, itri, m, val, iop)
    implicit none
    type(scorec_vector) :: v
    integer, intent(in) :: itri, m, iop
    vectype, intent(in), dimension(dofs_per_element) :: val

    integer, dimension(v%isize, dofs_per_element) :: irow
    integer :: i

    call get_element_indices(v%isize, itri, irow)

    do i=1,dofs_per_element
       call insert(v, irow(m,i), val(i), iop)
    end do
  end subroutine scorec_vector_insert_block

  subroutine scorec_vector_get_element_indices(isize, itri, ind)
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
       call get_node_indices(isize, inode(iii), temp)
       ind(:,i:i+dofs_per_node-1) = temp
       i = i + dofs_per_node
    end do
       
  end subroutine scorec_vector_get_element_indices

  logical function scorec_vector_is_nan(v)
    implicit none
    type(scorec_vector) :: v

    scorec_vector_is_nan = v%data(1).ne.v%data(1)
  end function scorec_vector_is_nan

  subroutine scorec_vector_finalize(v)
    implicit none
    type(scorec_vector) :: v
  end subroutine scorec_vector_finalize


  subroutine scorec_vector_sum_shared(v)
    implicit none
    type(scorec_vector) :: v

    call sumsharedppplvecvals(v%data)
  end subroutine scorec_vector_sum_shared


  subroutine scorec_vector_write(v, file)
    implicit none
    type(scorec_vector), intent(in) :: v
    character(len=*) :: file    

    open(unit=29, file=file, status='unknown')
    write(29, '(1e20.10)') v%data
    close(29)

  end subroutine scorec_vector_write


!==========================================================================


  !======================================================================
  ! node_index_vector
  ! ~~~~~~~~~~~~~~~~~
  !======================================================================
  integer function node_index_vector(v, inode, iplace)
    implicit none
    
    type(vector_type), intent(in) :: v
    integer, intent(in) :: inode
    integer, intent(in) :: iplace

    call get_node_index(inode, iplace, v%isize, node_index_vector)
  end function node_index_vector

  !======================================================================
  ! global_node_index
  ! ~~~~~~~~~~~~~~~~~
  !======================================================================
  integer function scorec_vector_global_node_index(v, inode, iplace)
    implicit none
    
    type(vector_type), intent(in) :: v
    integer, intent(in) :: inode
    integer, intent(in) :: iplace

    call get_global_node_index(inode, iplace, v%isize, &
         scorec_vector_global_node_index)
  end function scorec_vector_global_node_index

end module scorec_vector_mod
