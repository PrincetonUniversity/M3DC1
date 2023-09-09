module scorec_vector_mod
  use element

  type scorec_vector
     integer :: isize
     integer :: id
  end type scorec_vector

!!$  type trilinos_options
!!$     character(len=:), allocatable :: str
!!$  end type trilinos_options

  integer, parameter :: VEC_SET = 0
  integer, parameter :: VEC_ADD = 1
  real :: solver_tol
  character(len=50) :: krylov_solver
  character(len=50) :: preconditioner
  character(len=50) :: sub_dom_solver
  integer :: solver_type, num_iter
  integer :: graph_fill, subdomain_overlap, poly_ord
  real :: ilu_drop_tol, ilu_fill, ilu_omega

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

  interface mark_vector_for_solutiontransfer
     module procedure scorec_vector_mark_for_solutiontransfer
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
    call m3dc1_ent_getlocaldofid(0, inode-1, isize, ibegin,iendplusone)
    ibegin = ibegin + 1
    ind = ibegin + (iplace-1)*dofs_per_node
  end subroutine scorec_vector_get_node_index

  !========================================================
  ! get_global_node_index_scorec
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! returns the index of node inode in field iplace
  ! of a vector of size isize
  !========================================================
  subroutine scorec_vector_get_global_node_index(inode, iplace, isize, ind)

    implicit none
    integer, intent(in) :: iplace, isize, inode
    integer, intent(out) :: ind

    integer :: ibegin, iendplusone
    call m3dc1_ent_getglobaldofid(0, inode-1, isize, ibegin,iendplusone)
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
    call m3dc1_ent_getlocaldofid(0, inode-1, isize, ibegin,iendplusone)
    ibegin = ibegin +1
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
    call m3dc1_ent_getglobaldofid(0, inode-1, isize, ibegin,iendplusone)
    ibegin = ibegin + 1
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
    call m3dc1_field_add (v1%id, v2%id)
  end subroutine scorec_vector_add

  subroutine scorec_vector_multiply_real(v, s)
    implicit none

    type(scorec_vector), intent(inout) :: v
    real :: s
    call m3dc1_field_mult(v%id, s, 0)
  end subroutine scorec_vector_multiply_real

#ifdef USECOMPLEX
  subroutine scorec_vector_multiply_complex(v, s)
    implicit none

    type(scorec_vector), intent(inout) :: v
    complex :: s
    call m3dc1_field_mult(v%id, s, 1)
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
    call m3dc1_field_assign(v%id, s, 0)
  end subroutine scorec_vector_const_real


  subroutine scorec_vector_get_node_data_real(v, iplace, inode, data, rotate)
    use mesh_mod
    implicit none
    type(scorec_vector), intent(in) :: v
    integer, intent(in) :: inode, iplace
    real, intent(out), dimension(dofs_per_node) :: data
    logical, intent(in), optional :: rotate
    
    logical :: is_boundary
    integer :: izone, izonedim
    real :: normal(2), curv(3), x, phi, z
    vectype, dimension(dofs_per_node) :: temp1, temp2
    vectype, dimension(v%isize*dofs_per_node) :: dofs
    integer :: index, numDofs
    logical :: r

    if(present(rotate)) then 
       r = rotate
    else 
       r = .true.
    end if
    call m3dc1_ent_getdofdata ( 0, inode-1, v%id, numDofs, dofs)
    if(numDofs .ne. v%isize*dofs_per_node) then
       print *, "Error scorec_vector_get_node_data_real"
       call safestop(1)
    end if
    index = 1+(iplace-1)*dofs_per_node

    if(r) then
       temp1 = dofs(index:index+dofs_per_node-1)

       ! if node is on boundary, rotate data from n,t to R,Z
       call boundary_node(inode, is_boundary, izone, izonedim, &
            normal, curv, x, phi, z, BOUND_ANY)
       if(is_boundary) then
          call rotate_dofs(temp1, temp2, normal, curv, -1)
          data = temp2
       else
          data = temp1
       end if
    else
       data = dofs(index:index+dofs_per_node-1)
    endif

  end subroutine scorec_vector_get_node_data_real


#ifdef USECOMPLEX
  subroutine scorec_vector_get_node_data_complex(v, iplace, inode, data, &
       rotate)
    use mesh_mod
    implicit none
    type(scorec_vector), intent(in) :: v
    integer, intent(in) :: inode, iplace
    complex, intent(out), dimension(dofs_per_node) :: data
    logical, intent(in), optional :: rotate

    logical :: is_boundary
    integer :: izone, izonedim
    real :: normal(2), curv(3), x, phi, z
    vectype, dimension(dofs_per_node) :: temp1, temp2   
    vectype, dimension(v%isize*dofs_per_node) :: dofs
    integer :: index, numDofs
    logical :: r

    if(present(rotate)) then 
       r = rotate
    else 
       r = .true.
    end if

    call m3dc1_ent_getdofdata ( 0, inode-1, v%id, numDofs, dofs)
    index = 1+(iplace-1)*dofs_per_node

    if(r) then
       temp1 = dofs(index:index+dofs_per_node-1)

       ! if node is on boundary, rotate data from n,t to R,Z
       call boundary_node(inode, is_boundary, izone, izonedim, &
            normal, curv, x, phi, z, BOUND_ANY)
       if(is_boundary) then
          call rotate_dofs(temp1, temp2, normal, curv, -1)
          data = temp2
       else
          data = temp1
       end if
    else
       data = dofs(index:index+dofs_per_node-1) 
    endif
  end subroutine scorec_vector_get_node_data_complex
#endif

  subroutine scorec_vector_set_node_data_real(v, iplace, inode, data, rotate)
    use mesh_mod
    implicit none
    type(scorec_vector), intent(inout) :: v
    integer, intent(in) :: inode, iplace
    real, intent(in), dimension(dofs_per_node) :: data
    logical, intent(in), optional :: rotate

    logical :: is_boundary
    integer :: izone, izonedim
    real :: normal(2), curv(3), x, phi, z
    vectype, dimension(dofs_per_node) :: temp1, temp2

    vectype, dimension(v%isize*dofs_per_node) :: dofs
    integer :: index, numDofs
    logical :: r

    numDofs = v%isize*dofs_per_node
    if(present(rotate)) then
       r = rotate
    else
       r = .true.
    end if
    index = 1+(iplace-1)*dofs_per_node
    call m3dc1_ent_getdofdata ( 0, inode-1, v%id, numDofs, dofs)
    ! if node is on boundary, rotate data from R,Z to n,t
    if(r) then
       call boundary_node(inode, is_boundary, izone, izonedim, &
            normal, curv, x, phi, z, BOUND_ANY)
       if(is_boundary) then
          temp1 = data
          call rotate_dofs(temp1, temp2, normal, curv, 1)
          dofs(index:index+dofs_per_node-1) = temp2
          !call m3dc1_ent_setdofdata ( 0, inode-1, v%id, numDofs, dofs)
       else
          dofs(index:index+dofs_per_node-1) = data
       endif
    else
       dofs(index:index+dofs_per_node-1) = data
    endif
    call m3dc1_ent_setdofdata ( 0, inode-1, v%id, numDofs, dofs)
  end subroutine scorec_vector_set_node_data_real

#ifdef USECOMPLEX
  subroutine scorec_vector_set_node_data_complex(v, iplace, inode, data, rotate)
    use mesh_mod
    implicit none
    type(scorec_vector), intent(inout) :: v
    integer, intent(in) :: inode, iplace
    complex, intent(in), dimension(dofs_per_node) :: data
    logical, intent(in), optional :: rotate

    logical :: is_boundary
    integer :: izone, izonedim
    real :: normal(2), curv(3), x, phi, z
    vectype, dimension(dofs_per_node) :: temp1, temp2    
    vectype, dimension(v%isize*dofs_per_node) :: dofs
    integer :: index, numDofs
    logical :: r

    if(present(rotate)) then 
       r = rotate
    else 
       r = .true.
    end if
    numDofs = v%isize*dofs_per_node
    index = 1+(iplace-1)*dofs_per_node
    call m3dc1_ent_getdofdata ( 0, inode-1, v%id, numDofs, dofs)
    if(r) then
       ! if node is on boundary, rotate data from R,Z to n,t
       call boundary_node(inode, is_boundary, izone, izonedim, &
            normal, curv, x, phi, z, BOUND_ANY)
       if(is_boundary) then
          temp1 = data
          call rotate_dofs(temp1, temp2, normal, curv, 1)
          dofs(index:index+dofs_per_node-1) = temp2
       else
          dofs(index:index+dofs_per_node-1) = data
       endif
    else
       dofs(index:index+dofs_per_node-1) = data
    end if
    call m3dc1_ent_setdofdata ( 0, inode-1, v%id, numDofs, dofs)
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
    call m3dc1_field_copy(vout%id, vin%id)  
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
    call m3dc1_field_insert (v%id, i-1, 1, s, 0, iop) 
  end subroutine scorec_vector_insert_real

#ifdef USECOMPLEX
  subroutine scorec_vector_insert_complex(v, i, s, iop)
    implicit none
    
    type(scorec_vector), intent(inout) :: v
    complex, intent(in) :: s
    integer, intent(in) :: i, iop
    call m3dc1_field_insert (v%id, i-1, 1, s, 1,iop)
  end subroutine scorec_vector_insert_complex

#endif 

  !======================================================================
  ! create
  ! ~~~~~~
  ! creates a vector of size n with name prefix (optional)
  !======================================================================
  subroutine scorec_vector_create(f,n,prefix)
    implicit none

    type(scorec_vector), intent(inout) :: f
    integer, intent(in) :: n
    character(len=*), intent(in), optional :: prefix
    integer :: dataType
    character(len=32) :: field_name
    call m3dc1_field_genid (f%id)

    if(present(prefix)) then
      write(field_name,"(A16,A)") prefix, 0
      field_name = adjustl(field_name)
    else
      write(field_name,"(I0,A)") f%id, 0
    end if
    dataType = 0
#ifdef USECOMPLEX
    dataType = 1
#endif
    call m3dc1_field_create (f%id, trim(field_name), n, dataType, dofs_per_node);
    f%isize = n
  end subroutine scorec_vector_create

  subroutine scorec_vector_mark_for_solutiontransfer(f)
    implicit none

    type(scorec_vector), intent(inout) :: f
#ifdef ADAPT
    call m3dc1_field_mark4tx(f%id)
#endif
  end subroutine scorec_vector_mark_for_solutiontransfer

  !======================================================================
  ! destroy
  ! ~~~~~~~
  ! frees the memory allocated with create_field
  !======================================================================
  subroutine scorec_vector_destroy(f)
    implicit none

    type(scorec_vector), intent(inout) :: f
    integer :: i

    call m3dc1_field_exist(f%id, i)
    if(i .ne. 0) then
       call m3dc1_field_delete(f%id)
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
    integer :: isnan
    scorec_vector_is_nan = .false.
    call m3dc1_field_isnan(v%id, isnan)
    !scorec_vector_is_nan = isnan 
    if(isnan .ne. 0) scorec_vector_is_nan = .true. 
  end function scorec_vector_is_nan

  subroutine scorec_vector_finalize(v)
    implicit none
    type(scorec_vector) :: v
    call m3dc1_field_sync (v%id)
  end subroutine scorec_vector_finalize


  subroutine scorec_vector_sum_shared(v)
    implicit none
    type(scorec_vector) :: v

    call m3dc1_field_sum(v%id)
  end subroutine scorec_vector_sum_shared


  subroutine scorec_vector_write(v, file)
    implicit none
    type(scorec_vector), intent(in) :: v
    character(len=*) :: file    
    call m3dc1_field_write(v%id, file, 0);
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
