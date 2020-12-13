module field
  implicit none

  type field_type
     real, allocatable :: data(:)
  end type field_type

  interface assignment (=)
     module procedure copy_field
     module procedure const_field
  end interface

  interface add
     module procedure add_field_to_field
     module procedure add_real_to_node
     module procedure add_real_to_field
  end interface

  interface mult
     module procedure multiply_field_by_real
  end interface

  interface get_node_data
     module procedure field_get_node_data
  end interface

  interface set_node_data
     module procedure field_set_node_data
  end interface

  interface node_index
     module procedure node_index_field
  end interface

contains 


  !======================================================================
  ! create_field
  ! ~~~~~~~~~~~~
  ! creates a field and allocates its associated size-1 vector 
  !====================================================================== 
   subroutine create_field(f)
     use element
     use mesh
     implicit none

     type(field_type) :: f
     
     allocate(f%data(local_nodes()*dofs_per_node))
  end subroutine create_field

  !======================================================================
  ! destroy_field
  ! ~~~~~~~~~~~~~
  ! destroys a field created with create_field
  !====================================================================== 
  subroutine destroy_field(f)
    implicit none
     type(field_type) :: f

    deallocate(f%data)
  end subroutine destroy_field

  !======================================================================
  ! copy_field
  ! ~~~~~~~~~~
  ! destroys a field created with create_field
  !====================================================================== 
  subroutine copy_field(fout, fin)
    implicit none
     type(field_type), intent(out) :: fout
     type(field_type), intent(in) :: fin

     fout%data = fin%data
  end subroutine copy_field


  !======================================================================
  ! node_index
  ! ~~~~~~~~~~
  ! returns index of first dof of f associated with node inode
  !======================================================================
  integer function node_index_field(f, inode)
    use element
    implicit none
    
    type(field_type), intent(in) :: f
    integer, intent(in) :: inode

    node_index_field = (inode-1)*dofs_per_node + 1
  end function node_index_field

  !======================================================================
  ! set_node_data
  ! ~~~~~~~~~~~~~
  ! copies array data into dofs of f associated with node inode
  !======================================================================
  subroutine field_set_node_data(f, inode, data)
    use element
    implicit none
    
    type(field_type), intent(inout) :: f
    integer, intent(in) :: inode
    real, dimension(dofs_per_node), intent(in) :: data

    integer :: index
    index = node_index_field(f, inode)
    f%data(index:index+dofs_per_node-1) = data
  end subroutine field_set_node_data

  !======================================================================
  ! get_node_data
  ! ~~~~~~~~~~~~~
  ! populates data with dofs of f associated with node inode
  !======================================================================
  subroutine field_get_node_data(f, inode, data)
    use element
    implicit none
    
    type(field_type), intent(in) :: f
    integer, intent(in) :: inode
    real, dimension(dofs_per_node), intent(out) :: data

    integer :: index
    index = node_index_field(f, inode)
    data = f%data(index:index+dofs_per_node-1)
  end subroutine field_get_node_data

  !======================================================================
  ! add_to_node
  ! ~~~~~~~~~~~
  ! adds data to inode of field f
  !======================================================================
  subroutine add_real_to_node(f, inode, data)
    use element
    implicit none

    type(field_type), intent(inout) :: f
    integer, intent(in) :: inode
    real, dimension(dofs_per_node), intent(in) :: data
    real, dimension(dofs_per_node) :: d

    call get_node_data(f, inode, d)
    d = d + data
    call set_node_data(f, inode, d)
  end subroutine add_real_to_node

  !======================================================================
  ! add_scalar_to_field
  ! ~~~~~~~~~~~~~~~~~~~
  ! adds scalar val to f
  !======================================================================
  subroutine add_real_to_field(f, val)
    use element
    use mesh

    implicit none

    type(field_type), intent(inout) :: f
    real, intent(in) :: val

    integer :: inode, numnodes
    real, dimension(dofs_per_node) :: d

    d = 0.
    d(1) = val
    
    numnodes = owned_nodes()
    do inode=1,numnodes
       call add(f, inode, d)
    enddo
  end subroutine add_real_to_field


  !======================================================================
  ! add_field_to_field
  ! ~~~~~~~~~~~~~~~~~~
  ! adds field fin to fout
  !======================================================================
  subroutine add_field_to_field(fout, fin, factor)
    use mesh

    implicit none

    type(field_type), intent(inout) :: fout
    type(field_type), intent(in) :: fin
    real, optional :: factor

    integer :: inode, numnodes
    real, dimension(dofs_per_node) :: datain, dataout


    numnodes = owned_nodes()
    do inode=1,numnodes
       call get_node_data(fout, inode, dataout)
       call get_node_data(fin, inode, datain)
       if(present(factor)) datain = datain*factor
       dataout = dataout + datain
       call set_node_data(fout, inode, dataout)
    enddo
  end subroutine add_field_to_field


  !======================================================================
  ! multiply_field_by_real
  ! ~~~~~~~~~~~~~~~~~~~~~~
  ! multiplies fout by real val
  !======================================================================
  subroutine multiply_field_by_real(f, val)
    use mesh
    implicit none

    type(field_type), intent(inout) :: f
    real, intent(in) :: val

    integer :: inode, numnodes
    real, dimension(dofs_per_node) :: data

    numnodes = owned_nodes()
    do inode=1,numnodes
       call get_node_data(f, inode, data)
       data = data*val
       call set_node_data(f, inode, data)
    enddo
  end subroutine multiply_field_by_real

  !======================================================================
  ! const_field
  ! ~~~~~~~~~~~
  ! sets field fout to constant value val
  !======================================================================
  subroutine const_field(fout, val)
    use element
    use mesh

    implicit none

    type(field_type), intent(inout) :: fout
    real, intent(in) :: val

    integer :: inode, numnodes
    real, dimension(dofs_per_node) :: data

    data(1) = val
    data(2:dofs_per_node) = 0.

    numnodes = owned_nodes()
    do inode=1,numnodes
       call set_node_data(fout, inode, data)
    enddo
  end subroutine const_field


  !===========================================================
  ! get_element_dofs
  ! ~~~~~~~~~~~~~~~~
  ! get dofs associated with element itri
  !===========================================================
  subroutine get_element_dofs(fin, itri, dofs)
    use element
    use mesh

    implicit none

    type(field_type), intent(in) :: fin
    integer, intent(in) :: itri
    real, dimension(dofs_per_element), intent(out) :: dofs

    integer :: i, iii
    integer, dimension(nodes_per_element) :: inode
    
    logical :: is_boundary
    integer :: izone, izonedim
    real :: normal(2), curv, x, z
    real, dimension(dofs_per_element) :: temp

    call get_element_nodes(itri, inode)
    i = 1
    do iii=1, nodes_per_element
       call get_node_data(fin, inode(iii), dofs(i:i+dofs_per_node-1))
       
       call boundary_node(inode(iii), is_boundary, izone, izonedim, &
            normal, curv, x, z)
       if(is_boundary) then
          temp = dofs(i:i+dofs_per_node-1)
          call rotate_dofs(temp, dofs(i:i+dofs_per_node-1), normal, curv, 1)
       endif

       i = i + dofs_per_node
    enddo
    
  end subroutine get_element_dofs

  !==========================================================
  ! calcavector
  ! ~~~~~~~~~~~
  !
  ! calculates the polynomial coefficients avector
  ! of field fin in element itri.
  !==========================================================
  subroutine calcavector(itri, fin, avector)
    use element
    implicit none
    
    integer, intent(in) :: itri
    type(field_type), intent(in) :: fin

    real, dimension(coeffs_per_element), intent(out) :: avector
    real, dimension(dofs_per_element) :: dofs

    call get_element_dofs(fin, itri, dofs)
    call local_coeffs(itri, dofs, avector)
  end subroutine calcavector


end module field
