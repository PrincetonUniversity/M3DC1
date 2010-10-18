module scorec_vector_mod
  use element

  type scorec_vector
     vectype, allocatable :: data(:)
     integer :: isize
  end type scorec_vector

  integer, parameter :: VEC_SET = 0
  integer, parameter :: VEC_ADD = 1

  interface assignment(=)
     module procedure copy_scorec_vec 
     module procedure const_scorec_vec_real
  end interface

  interface add
     module procedure add_scorec_vectors
  end interface

  interface mult
     module procedure multiply_scorec_vector_by_real
#ifdef USECOMPLEX
     module procedure multiply_scorec_vector_by_complex
#endif
  end interface

  interface create_vector
     module procedure create_scorec_vec
  end interface

  interface destroy_vector
     module procedure destroy_scorec_vec
  end interface

  interface insert
     module procedure insert_real_scorec_vec
#ifdef USECOMPLEX
     module procedure insert_complex_scorec_vec
#endif
  end interface

  interface vector_insert_block
     module procedure scorec_vector_insert_block
  end interface

  interface get_basis_indices
     module procedure get_basis_indices_scorec
  end interface

  interface node_index
     module procedure node_index_vector
  end interface

  interface global_node_index
     module procedure global_node_index_vector
  end interface

  interface get_dof_indices
     module procedure get_dof_indices_scorec
  end interface

  interface get_global_dof_indices
     module procedure get_global_dof_indices_scorec
  end interface

  interface get_node_index
     module procedure get_node_index_scorec
  end interface

  interface get_global_node_index
     module procedure get_global_node_index_scorec
  end interface

  interface set_vector_node_data
     module procedure set_vector_node_data_real
#ifdef USECOMPLEX
     module procedure set_vector_node_data_complex
#endif
  end interface
  interface get_vector_node_data
     module procedure get_vector_node_data_real
#ifdef USECOMPLEX
     module procedure get_vector_node_data_complex
#endif
  end interface

  interface is_nan
     module procedure scorec_is_nan
  end interface

  interface finalize
     module procedure finalize_scorec_vector
  end interface

  interface flush
     module procedure flush_scorec_vector
  end interface

  interface sum_shared
     module procedure sum_shared_scorec_vector
  end interface

contains

  !========================================================
  ! dof_index
  ! ~~~~~~~~~
  ! Get the dof index where ifirst is the first index associated with the node
  ! iplace is the field place, and idof is the local dof number
  !========================================================
  integer function dof_index(ifirst, iplace, idof)
    integer, intent(in) :: ifirst, iplace, idof
    dof_index = ifirst + (iplace-1)*dofs_per_node + idof-1
  end function dof_index
  

  !========================================================
  ! get_node_index_scorec
  ! ~~~~~~~~~~~~~~~~~~~~~
  ! returns the index of node inode in field iplace
  ! of a vector of size isize
  !========================================================
  subroutine get_node_index_scorec(inode, iplace, isize, ind)

    implicit none
    integer, intent(in) :: iplace, isize, inode
    integer, intent(out) :: ind

    integer :: ibegin, iendplusone
    call entdofs(isize, inode, 0, ibegin, iendplusone)
    ind = ibegin + (iplace-1)*dofs_per_node
  end subroutine get_node_index_scorec

  !========================================================
  ! get_global_node_index_scorec
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! returns the index of node inode in field iplace
  ! of a vector of size isize
  !========================================================
  subroutine get_global_node_index_scorec(inode, iplace, isize, ind)

    implicit none
    integer, intent(in) :: iplace, isize, inode
    integer, intent(out) :: ind

    integer :: ibegin, iendplusone
    call globalentdofs(isize, inode, 0, ibegin, iendplusone)
    ind = ibegin + (iplace-1)*dofs_per_node
  end subroutine get_global_node_index_scorec


  !========================================================
  ! get_dof_indices
  ! ~~~~~~~~~~~~~~~
  ! returns the indices of node inode for each field
  ! of a vector of size isize
  !========================================================
  subroutine get_dof_indices_scorec(isize, inode, ind)

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

  end subroutine get_dof_indices_scorec

  !========================================================
  ! get_global_dof_indices
  ! ~~~~~~~~~~~~~~~~~~~~~~
  ! returns the indices of node inode for each field
  ! of a vector of size isize
  !========================================================
  subroutine get_global_dof_indices_scorec(isize, inode, ind)

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
  end subroutine get_global_dof_indices_scorec


  !======================================================================
  ! add_scorec_vectors
  ! ~~~~~~~~~~~~~~~~~~
  ! adds v2 to v1.  result is stored in v1
  !======================================================================
  subroutine add_scorec_vectors(v1, v2)
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
  end subroutine add_scorec_vectors

  subroutine multiply_scorec_vector_by_real(v, s)
    implicit none

    type(scorec_vector), intent(inout) :: v
    real :: s

    if(.not.allocated(v%data)) print *, 'Error: vector not allocated'
    v%data = v%data * s
  end subroutine multiply_scorec_vector_by_real

#ifdef USECOMPLEX
  subroutine multiply_scorec_vector_by_complex(v, s)
    implicit none

    type(scorec_vector), intent(inout) :: v
    complex :: s

    if(.not.allocated(v%data)) print *, 'Error: vector not allocated'
    v%data = v%data * s
  end subroutine multiply_scorec_vector_by_complex
#endif


  !======================================================================
  ! const_scorec_vec_real
  ! ~~~~~~~~~~~~~~~~~~~~~
  ! set all elements of v to s
  !======================================================================
  subroutine const_scorec_vec_real(v,s)
    implicit none

    type(scorec_vector), intent(inout) :: v
    real, intent(in) :: s

    if(.not.allocated(v%data)) print *, 'Error: vector not allocated'
    v%data = s
  end subroutine const_scorec_vec_real



  subroutine get_vector_node_data_real(v, iplace, inode, data)
    implicit none
    type(scorec_vector), intent(in) :: v
    integer, intent(in) :: inode, iplace
    real, intent(out), dimension(dofs_per_node) :: data
    
    integer :: index
    call get_node_index(inode, iplace, v%isize, index)
    data = v%data(index:index+dofs_per_node-1)
  end subroutine get_vector_node_data_real


#ifdef USECOMPLEX
  subroutine get_vector_node_data_complex(v, iplace, inode, data)
    implicit none
    type(scorec_vector), intent(in) :: v
    integer, intent(in) :: inode, iplace
    vectype, intent(out), dimension(dofs_per_node) :: data
    
    integer :: index
    call get_node_index(inode, iplace, v%isize, index)
    data = v%data(index:index+dofs_per_node-1)
  end subroutine get_vector_node_data_complex
#endif

  subroutine set_vector_node_data_real(v, iplace, inode, data)
    implicit none
    type(scorec_vector), intent(inout) :: v
    integer, intent(in) :: inode, iplace
    real, intent(in), dimension(dofs_per_node) :: data
    
    integer :: index
    call get_node_index(inode, iplace, v%isize, index)
    v%data(index:index+dofs_per_node-1) = data
  end subroutine set_vector_node_data_real

#ifdef USECOMPLEX
  subroutine set_vector_node_data_complex(v, iplace, inode, data)
    implicit none
    type(scorec_vector), intent(inout) :: v
    integer, intent(in) :: inode, iplace
    complex, intent(in), dimension(dofs_per_node) :: data
    
    integer :: index
    call get_node_index(inode, iplace, v%isize, index)
    v%data(index:index+dofs_per_node-1) = data
  end subroutine set_vector_node_data_complex
#endif



  !======================================================================
  ! copy_scorec_vec
  ! ~~~~~~~~~~~~~~~
  ! copy data from vin to vout
  !======================================================================
  subroutine copy_scorec_vec(vout,vin)
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

  end subroutine copy_scorec_vec

  !======================================================================
  ! insert_scorec_vec
  ! ~~~~~~~~~~~~~~~~~
  ! set element i of v to s
  !======================================================================
  subroutine insert_real_scorec_vec(v, i, s, iop)
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
  end subroutine insert_real_scorec_vec

#ifdef USECOMPLEX
  subroutine insert_complex_scorec_vec(v, i, s, iop)
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
  end subroutine insert_complex_scorec_vec

#endif 

  !======================================================================
  ! create_scorec_vec
  ! ~~~~~~~~~~~~~~~~~
  ! creates a vector of size n
  !======================================================================
  subroutine create_scorec_vec(f,n)
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
  end subroutine create_scorec_vec

  !======================================================================
  ! destroy_scorec_vec
  ! ~~~~~~~~~~~~~~~~~~
  ! frees the memory allocated with create_field
  !======================================================================
  subroutine destroy_scorec_vec(f)
    implicit none

    type(scorec_vector), intent(inout) :: f
    integer :: i

    call checkppplveccreated(f%data, i)
    if(i .ne. 0) then
       call deleteppplvec(f%data)
       deallocate(f%data)
    endif
  end subroutine destroy_scorec_vec

  subroutine scorec_vector_insert_block(v, itri, m, val, iop)
    implicit none
    type(scorec_vector) :: v
    integer, intent(in) :: itri, m, iop
    vectype, intent(in), dimension(dofs_per_element) :: val

    integer, dimension(v%isize, dofs_per_element) :: irow
    integer :: i

    call get_basis_indices(v%isize, itri, irow)

    do i=1,dofs_per_element
       call insert(v, irow(m,i), val(i), iop)
    end do
  end subroutine scorec_vector_insert_block

  subroutine get_basis_indices_scorec(isize, itri, ind)
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
       
  end subroutine get_basis_indices_scorec

  logical function scorec_is_nan(v)
    implicit none
    type(scorec_vector) :: v

    scorec_is_nan = v%data(1).ne.v%data(1)
  end function scorec_is_nan

  subroutine finalize_scorec_vector(v)
    implicit none
    type(scorec_vector) :: v
  end subroutine finalize_scorec_vector

  subroutine flush_scorec_vector(v)
    implicit none
    type(scorec_vector) :: v
  end subroutine flush_scorec_vector


  subroutine sum_shared_scorec_vector(v)
    implicit none
    type(scorec_vector) :: v

    call sumsharedppplvecvals(v%data)
  end subroutine sum_shared_scorec_vector


  subroutine write_vector(v, file)
    implicit none
    type(scorec_vector), intent(in) :: v
    character(len=*) :: file    

    open(unit=29, file=file, status='unknown')
    write(29, '(1e20.10)') v%data
    close(29)

  end subroutine write_vector


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
  ! global_node_index_vector
  ! ~~~~~~~~~~~~~~~~~~~~~~~~
  !======================================================================
  integer function global_node_index_vector(v, inode, iplace)
    implicit none
    
    type(vector_type), intent(in) :: v
    integer, intent(in) :: inode
    integer, intent(in) :: iplace

    call get_global_node_index(inode, iplace, v%isize, &
         global_node_index_vector)
  end function global_node_index_vector

end module scorec_vector_mod
