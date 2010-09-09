module field
  use vector_mod

  implicit none

  type field_type
     type(vector_type), pointer :: vec
     integer :: index
  end type field_type

  interface assignment (=)
     module procedure copy_field
     module procedure const_field_real
#ifdef USECOMPLEX
     module procedure const_field_complex
#endif
  end interface

  interface add
     module procedure add_field_to_field
     module procedure add_real_to_node
     module procedure add_real_to_field
#ifdef USECOMPLEX
     module procedure add_complex_to_node
     module procedure add_complex_to_field
#endif
  end interface

  interface mult
     module procedure multiply_field_by_real
#ifdef USECOMPLEX
     module procedure multiply_field_by_complex
#endif
  end interface

  interface pow
     module procedure raise_field_to_real_power
  end interface

  interface set_node_data
     module procedure set_node_data_real
#ifdef USECOMPLEX
     module procedure set_node_data_complex
#endif
  end interface

  interface node_index
     module procedure node_index_field
  end interface

  interface matvecmult
     module procedure matvecmult_field_vec
     module procedure matvecmult_vec_field
  end interface

  interface is_nan
     module procedure field_is_nan
  end interface
contains 

!!$  subroutine field_get_local_pointer(f, inode, ptr)
!!$    implicit none
!!$
!!$    type(field_type) :: f
!!$    integer, intent(in) :: inode
!!$    vectype, pointer, intent(out) :: ptr(:)
!!$
!!$    call vector_get_local_pointer(f%vec, f%index, inode, ptr)
!!$  end subroutine field_get_local_pointer


  !======================================================================
  ! create_field
  ! ~~~~~~~~~~~~
  ! creates a field and allocates its associated size-1 vector 
  !====================================================================== 
   subroutine create_field(f)
    implicit none

    type(field_type), intent(inout) :: f

    allocate(f%vec)
    call create_vector(f%vec, 1)
    f%index = 1
  end subroutine create_field

  !======================================================================
  ! destroy_field
  ! ~~~~~~~~~~~~~
  ! destroys a field created with create_field
  !====================================================================== 
  subroutine destroy_field(f)
    implicit none

    type(field_type), intent(inout) :: f

    call destroy_vector(f%vec)
    if(associated(f%vec)) deallocate(f%vec)
    nullify(f%vec)
  end subroutine destroy_field

  !======================================================================
  ! associate_field
  ! ~~~~~~~~~~~~~~~
  ! associates field f with the iplaceth field in vector vec, 
  ! where vec contains isize total fields
  !====================================================================== 
  subroutine associate_field(f, vec, iplace)
    implicit none
    
    type(field_type), intent(inout) :: f
    type(vector_type), target, intent(in) :: vec
    integer, intent(in) :: iplace
   
    if(iplace.gt.vec%isize) then
       print *, 'Error: field index out of range.'
       return
    endif

    f%vec => vec
    f%index = iplace
  end subroutine associate_field

  !======================================================================
  ! node_index
  ! ~~~~~~~~~~
  ! returns index of first dof of f associated with node inode
  !======================================================================
  integer function node_index_field(f, inode)
    implicit none
    
    type(field_type), intent(in) :: f
    integer, intent(in) :: inode

    node_index_field = node_index(f%vec, inode, f%index)
  end function node_index_field

  !======================================================================
  ! set_node_data
  ! ~~~~~~~~~~~~~
  ! copies array data into dofs of f associated with node inode
  !======================================================================
  subroutine set_node_data_real(f, inode, data)
    use element
    implicit none
    
    type(field_type), intent(inout) :: f
    integer, intent(in) :: inode
    real, dimension(dofs_per_node), intent(in) :: data

    call set_vector_node_data(f%vec,f%index,inode,data)
  end subroutine set_node_data_real

#ifdef USECOMPLEX
  subroutine set_node_data_complex(f, inode, data)
    use element
    implicit none
    
    type(field_type), intent(inout) :: f
    integer, intent(in) :: inode
    complex, dimension(dofs_per_node), intent(in) :: data
    
    call set_vector_node_data(f%vec,f%index,inode,data)
  end subroutine set_node_data_complex
#endif


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

    call get_vector_node_data(f%vec, f%index, inode, d)
    d = d + data
    call set_vector_node_data(f%vec, f%index, inode, d)
  end subroutine add_real_to_node

#ifdef USECOMPLEX
  subroutine add_complex_to_node(f, inode, data)
    use element
    implicit none

    type(field_type), intent(inout) :: f
    integer, intent(in) :: inode
    complex, dimension(dofs_per_node), intent(in) :: data
    complex, dimension(dofs_per_node) :: d

    call get_vector_node_data(f%vec, f%index, inode, d)
    d = d + data
    call set_vector_node_data(f%vec, f%index, inode, d)
  end subroutine add_complex_to_node
#endif


  !======================================================================
  ! get_node_data
  ! ~~~~~~~~~~~~~
  ! populates data with dofs of f associated with node inode
  !======================================================================
  subroutine get_node_data(f, inode, data)
    use element
    implicit none
    
    type(field_type), intent(in) :: f
    integer, intent(in) :: inode
    vectype, dimension(dofs_per_node), intent(out) :: data

    call get_vector_node_data(f%vec, f%index, inode, data)
  end subroutine get_node_data


  !======================================================================
  ! add_scalar_to_field
  ! ~~~~~~~~~~~~~~~~~~~
  ! adds scalar val to f
  !======================================================================
  subroutine add_real_to_field(f, val)
    use element
    use mesh_mod

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
    call finalize(f%vec)
  end subroutine add_real_to_field

#ifdef USECOMPLEX
  subroutine add_complex_to_field(f, val)
    use element
    use mesh_mod

    implicit none

    type(field_type), intent(inout) :: f
    complex, intent(in) :: val

    integer :: inode, numnodes
    complex, dimension(dofs_per_node) :: d

    d = 0.
    d(1) = val

    numnodes = owned_nodes()
    do inode=1,numnodes
       call add(f, inode, d)
    enddo
    call finalize(f%vec)
  end subroutine add_complex_to_field
#endif

  !======================================================================
  ! add_field_to_field
  ! ~~~~~~~~~~~~~~~~~~
  ! adds field fin to fout
  !======================================================================
  subroutine add_field_to_field(fout, fin, factor)
    use mesh_mod

    implicit none

    type(field_type), intent(inout) :: fout
    type(field_type), intent(in) :: fin
    real, optional :: factor

    integer :: inode, numnodes
    vectype, dimension(dofs_per_node) :: datain, dataout


    numnodes = owned_nodes()
    do inode=1,numnodes
       call get_node_data(fout, inode, dataout)
       call get_node_data(fin, inode, datain)
       if(present(factor)) datain = datain*factor
       dataout = dataout + datain
       call set_node_data(fout, inode, dataout)
    enddo
    call finalize(fout%vec)
  end subroutine add_field_to_field


  !======================================================================
  ! multiply_field_by_real
  ! ~~~~~~~~~~~~~~~~~~~~~~
  ! multiplies fout by real val
  !======================================================================
  subroutine multiply_field_by_real(f, val)
    use mesh_mod
    implicit none

    type(field_type), intent(inout) :: f
    real, intent(in) :: val

    integer :: inode, numnodes
    vectype, dimension(dofs_per_node) :: data

    numnodes = owned_nodes()
    do inode=1,numnodes
       call get_node_data(f, inode, data)
       data = data*val
       call set_node_data(f, inode, data)
    enddo
    call finalize(f%vec)
  end subroutine multiply_field_by_real

#ifdef USECOMPLEX
  subroutine multiply_field_by_complex(f, val)
    use mesh_mod

    implicit none

    type(field_type), intent(inout) :: f
    complex, intent(in) :: val

    integer :: inode, numnodes
    vectype, dimension(dofs_per_node) :: data

    numnodes = owned_nodes()
    do inode=1,numnodes
       call get_node_data(f, inode, data)
       data = data*val
       call set_node_data(f, inode, data)
    enddo
    call finalize(f%vec)
  end subroutine multiply_field_by_complex
#endif

  subroutine raise_field_to_real_power(f, val)
    use element
    use mesh_mod
    implicit none

    type(field_type), intent(inout) :: f
    real, intent(in) :: val
    vectype, dimension(dofs_per_node) :: data, new_data

    integer :: inode, numnodes

    numnodes = owned_nodes()
    do inode=1,numnodes
       call get_node_data(f, inode, data)
       
       if(data(1).eq.0.) cycle
          
       if(val.eq.0.) then
          new_data(1) = 1.
          new_data(2:6) = 0.
       else
          new_data(1) =     data(1)**val
          new_data(2) = val*new_data(1)/data(1) * data(2)
          new_data(3) = val*new_data(1)/data(1) * data(3)
          new_data(4) = val*(val-1.)*new_data(1)/data(1)**2 * data(2)**2 &
               + val*new_data(1)/data(1) * data(4)
          new_data(5) = val*(val-1.)*new_data(1)/data(1)**2 * data(2)*data(3) &
               + val*new_data(1)/data(1) * data(5)
          new_data(6) = val*(val-1.)*new_data(1)/data(1)**2 * data(3)**2 &
               + val*new_data(1)/data(1) * data(6)
       endif
       call set_node_data(f, inode, new_data)
    enddo
    call finalize(f%vec)
  end subroutine raise_field_to_real_power



  !======================================================================
  ! const_field
  ! ~~~~~~~~~~~
  ! sets field fout to constant value val
  !======================================================================
  subroutine const_field_real(fout, val)
    use element
    use mesh_mod

    implicit none

    type(field_type), intent(inout) :: fout
    real, intent(in) :: val

    integer :: inode, numnodes
    vectype, dimension(dofs_per_node) :: data

    data(1) = val
    data(2:dofs_per_node) = 0.

    numnodes = owned_nodes()
    do inode=1,numnodes
       call set_node_data(fout, inode, data)
    enddo
    call finalize(fout%vec)
  end subroutine const_field_real

#ifdef USECOMPLEX
  subroutine const_field_complex(fout, val)
    use element
    use mesh_mod

    implicit none

    type(field_type), intent(inout) :: fout
    complex, intent(in) :: val

    integer :: inode, numnodes
    vectype, dimension(dofs_per_node) :: data

    data(1) = val
    data(2:dofs_per_node-1) = 0.

    numnodes = owned_nodes()
    
    do inode=1,numnodes
       call set_node_data(fout, inode, data)
    enddo
    call finalize(fout%vec)
  end subroutine const_field_complex
#endif

  !======================================================================
  ! copy_field
  ! ~~~~~~~~~~
  ! copies data from fin to fout
  !======================================================================
  subroutine copy_field(fout, fin)
    use element
    use vector_mod
    use mesh_mod

    implicit none

    type(field_type), intent(inout) :: fout
    type(field_type), intent(in) :: fin

    integer :: inode, numnodes
    vectype, dimension(dofs_per_node) :: data

    if(fin%vec%isize.eq.1 .and. fout%vec%isize.eq.1) then
       fout%vec = fin%vec
    else 
       numnodes = owned_nodes()
    
       do inode=1,numnodes
          call get_node_data(fin,inode,data)
          call set_node_data(fout,inode,data)
       enddo

       call finalize(fout%vec)
    end if
  end subroutine copy_field


  logical function field_is_nan(f)
    implicit none
    type(field_type), intent(in) :: f

    field_is_nan = is_nan(f%vec)
  end function field_is_nan

  !==========================================================
  ! calcavector
  ! ~~~~~~~~~~~
  !
  ! calculates the 20 polynomial coefficients avector
  ! of field inarr in element itri.
  !==========================================================
  subroutine calcavector(itri, fin, avector)
    use element
    use mesh_mod
    
    implicit none
    
    integer, intent(in) :: itri
    type(field_type), intent(in) :: fin
    vectype, dimension(20), intent(out) :: avector
    
    integer :: i, iii, k
    integer, dimension(nodes_per_element) :: inode
    vectype, dimension(nodes_per_element*dofs_per_node) :: wlocal

    call get_element_nodes(itri, inode)
    i = 1
    do iii=1, nodes_per_element
       call get_node_data(fin, inode(iii), wlocal(i:i+dofs_per_node-1))
       i = i + dofs_per_node
    enddo
    
    ! calculate the function value corresponding to this point
    do i=1,20
       avector(i) = 0.
       do k=1,nodes_per_element*dofs_per_node
          avector(i) = avector(i) + gtri_old(i,k,itri)*wlocal(k)
       enddo
    enddo
  end subroutine calcavector

  subroutine matvecmult_field_vec(mat,fin,vout)
    use vector_mod
    use matrix_mod

    implicit none

    type(matrix_type), intent(in) :: mat
    type(field_type), intent(in) :: fin
    type(vector_type), intent(inout) :: vout
    type(field_type) :: fin_new

    if(fin%vec%isize .eq. 1) then
       call matvecmult(mat, fin%vec, vout)
    else
       call create_field(fin_new)
       fin_new = fin
       call matvecmult(mat, fin_new%vec, vout)
       call destroy_field(fin_new)
    endif
  end subroutine matvecmult_field_vec

  subroutine matvecmult_vec_field(mat,vin,fout)
    use vector_mod
    use matrix_mod

    implicit none

    type(matrix_type), intent(in) :: mat
    type(vector_type), intent(in) :: vin
    type(field_type), intent(inout) :: fout
    type(field_type) :: fout_new

    if(fout%vec%isize .eq. 1) then
       call matvecmult(mat, vin, fout%vec)
    else
       call create_field(fout_new)
       call matvecmult(mat, vin, fout_new%vec)
       fout = fout_new
       call destroy_field(fout_new)
    endif
  end subroutine matvecmult_vec_field


end module field
