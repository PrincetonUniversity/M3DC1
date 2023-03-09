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

  interface get_node_data
     module procedure field_get_node_data
  end interface

  interface set_node_data
     module procedure field_set_node_data_real
#ifdef USECOMPLEX
     module procedure field_set_node_data_complex
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
  subroutine create_field(f,prefix)
    implicit none

    type(field_type), intent(inout) :: f
    character(len=*), intent(in), optional :: prefix

    allocate(f%vec)
    if (present(prefix)) then
      call create_vector(f%vec, 1, prefix)
    else
      call create_vector(f%vec, 1)
    end if
    f%index = 1
  end subroutine create_field

  !======================================================================
  ! mark_field_for_solutiontransfer
  ! ~~~~~~~~~~~~
  ! creates a field and allocates its associated size-1 vector
  !======================================================================
   subroutine mark_field_for_solutiontransfer(f)
    implicit none

    type(field_type), intent(inout) :: f
    call mark_vector_for_solutiontransfer(f%vec)
  end subroutine mark_field_for_solutiontransfer

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
  subroutine field_set_node_data_real(f, inode, data, rotate)
    use element
    implicit none
    
    type(field_type), intent(inout) :: f
    integer, intent(in) :: inode
    real, dimension(dofs_per_node), intent(in) :: data
    logical, intent(in), optional :: rotate

    call set_node_data(f%vec,f%index,inode,data,rotate)
  end subroutine field_set_node_data_real

#ifdef USECOMPLEX
  subroutine field_set_node_data_complex(f, inode, data, rotate)
    use element
    implicit none
    
    type(field_type), intent(inout) :: f
    integer, intent(in) :: inode
    complex, dimension(dofs_per_node), intent(in) :: data
    logical, intent(in), optional :: rotate
    
    call set_node_data(f%vec,f%index,inode,data,rotate)
  end subroutine field_set_node_data_complex
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

    call get_node_data(f%vec, f%index, inode, d)
    d = d + data
    call set_node_data(f%vec, f%index, inode, d)
  end subroutine add_real_to_node

#ifdef USECOMPLEX
  subroutine add_complex_to_node(f, inode, data)
    use element
    implicit none

    type(field_type), intent(inout) :: f
    integer, intent(in) :: inode
    complex, dimension(dofs_per_node), intent(in) :: data
    complex, dimension(dofs_per_node) :: d

    call get_node_data(f%vec, f%index, inode, d)
    d = d + data
    call set_node_data(f%vec, f%index, inode, d)
  end subroutine add_complex_to_node
#endif


  !======================================================================
  ! get_node_data
  ! ~~~~~~~~~~~~~
  ! populates data with dofs of f associated with node inode
  !======================================================================
  subroutine field_get_node_data(f, inode, data, rotate)
    use element
    implicit none
    
    type(field_type), intent(in) :: f
    integer, intent(in) :: inode
    vectype, dimension(dofs_per_node), intent(out) :: data
    logical, intent(in), optional :: rotate

    call get_node_data(f%vec, f%index, inode, data, rotate)
  end subroutine field_get_node_data


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

    integer :: inode, numnodes, icounter_t 
    real, dimension(dofs_per_node) :: d

    d = 0.
    d(1) = val
    
    numnodes = owned_nodes()
    do icounter_t=1,numnodes
       inode = nodes_owned(icounter_t)
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

    integer :: inode, numnodes, icounter_t
    complex, dimension(dofs_per_node) :: d

    d = 0.
    d(1) = val

    numnodes = owned_nodes()
    do icounter_t=1,numnodes
       inode = nodes_owned(icounter_t)
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

    integer :: inode, numnodes, icounter_t
    vectype, dimension(dofs_per_node) :: datain, dataout


    numnodes = owned_nodes()
    do icounter_t=1,numnodes
       inode = nodes_owned(icounter_t)
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

    integer :: inode, numnodes, icounter_t
    vectype, dimension(dofs_per_node) :: data

    numnodes = owned_nodes()
    do icounter_t=1,numnodes
       inode = nodes_owned(icounter_t)
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

    integer :: inode, numnodes, icounter_t
    vectype, dimension(dofs_per_node) :: data

    numnodes = owned_nodes()
    do icounter_t=1,numnodes
       inode = nodes_owned(icounter_t)
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

    integer :: inode, numnodes, icounter_t

    numnodes = owned_nodes()
    do icounter_t=1,numnodes
       inode = nodes_owned(icounter_t)
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

    integer :: inode, numnodes, icounter_t
    vectype, dimension(dofs_per_node) :: data

    data(1) = val
    data(2:dofs_per_node) = 0.

    numnodes = owned_nodes()
    do icounter_t=1,numnodes
       inode = nodes_owned(icounter_t)
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

    integer :: inode, numnodes,icounter_t
    vectype, dimension(dofs_per_node) :: data

    data(1) = val
    data(2:dofs_per_node-1) = 0.

    numnodes = owned_nodes()
    do icounter_t=1,numnodes
       inode = nodes_owned(icounter_t) 
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

    integer :: inode, numnodes, icounter_t
    vectype, dimension(dofs_per_node) :: data

    if(fin%vec%isize.eq.1 .and. fout%vec%isize.eq.1) then
       fout%vec = fin%vec
    else 
       numnodes = owned_nodes()
       do icounter_t=1,numnodes
          inode = nodes_owned(icounter_t) 
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

  !===========================================================
  ! get_element_dofs
  ! ~~~~~~~~~~~~~~~~
  ! get dofs associated with element itri
  !===========================================================
  subroutine get_element_dofs(fin, itri, dofs)
    use element
    use mesh_mod

    implicit none

    type(field_type), intent(in) :: fin
    integer, intent(in) :: itri
    vectype, dimension(dofs_per_element), intent(out) :: dofs

    integer :: i, iii
    integer, dimension(nodes_per_element) :: inode

    call get_element_nodes(itri, inode)
    i = 1
    do iii=1, nodes_per_element
       call get_node_data(fin, inode(iii), dofs(i:i+dofs_per_node-1), .false.)

       i = i + dofs_per_node
    enddo
    
  end subroutine get_element_dofs

  !===========================================================
  ! set_element_dofs
  ! ~~~~~~~~~~~~~~~~
  ! set dofs associated with element itri
  !===========================================================
  subroutine set_element_dofs(fout, itri, dofs)
    use element
    use mesh_mod

    implicit none

    type(field_type), intent(inout) :: fout
    integer, intent(in) :: itri
    vectype, dimension(dofs_per_element), intent(in) :: dofs

    integer :: i, iii
    integer, dimension(nodes_per_element) :: inode

    call get_element_nodes(itri, inode)
    i = 1
    do iii=1, nodes_per_element
       call set_node_data(fout, inode(iii), dofs(i:i+dofs_per_node-1), .true.)

       i = i + dofs_per_node
    enddo    
  end subroutine set_element_dofs


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

    vectype, dimension(coeffs_per_element), intent(out) :: avector
    vectype, dimension(dofs_per_element) :: dofs

    call get_element_dofs(fin, itri, dofs)
    call local_coeffs(itri, dofs, avector)
  end subroutine calcavector

  !==========================================================
  ! setavector
  ! ~~~~~~~~~~
  !
  ! sets dofs of field of field fout in element itri 
  ! given the polynomial coefficients avector
  !==========================================================
  subroutine setavector(itri, fout, avector)
    use mesh_mod
    implicit none
    
    integer, intent(in) :: itri
    type(field_type), intent(inout) :: fout

    vectype, dimension(coeffs_per_element), intent(in) :: avector
    vectype, dimension(dofs_per_element) :: dofs

    call local_dofs(itri, dofs, avector)
    call set_element_dofs(fout, itri, dofs)
  end subroutine setavector


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



  subroutine check_axisymmetry(fin, name)
    use element
    use mesh_mod
    implicit none

    type(field_type), intent(in) :: fin
    character(len=*), intent(in) :: name

    integer :: inode(nodes_per_element), nelm, itri, i
    vectype, dimension(dofs_per_node) :: data1, data2
    logical :: is_axisymmetric

    is_axisymmetric = .true.
#ifdef USE3D
    nelm = local_elements()
    do itri=1, nelm
       call get_element_nodes(itri, inode)

       do i=1, pol_nodes_per_element
          call get_node_data(fin, inode(i), data1)
          call get_node_data(fin, inode(i+pol_nodes_per_element), data2)
          if(maxval(abs(data1(1:6) - data2(1:6))) .gt. 1e-8) then
             print *, name, ' NOT AXISYMMETRIC AT ', itri, i, inode(i)
             write(*, '(6e12.4)') data1(1:6), data2(1:6)
             is_axisymmetric = .false.
          endif
       end do
    end do
#endif

    if(is_axisymmetric) then
       print *, name, ' IS AXISYMMETRIC'
    endif

  end subroutine check_axisymmetry


  !======================================================================
  ! straighten_field
  ! ~~~~~~~~~~~~~~~~
  ! rotates dofs at boundary nodes to (R,Z) coordinates
  !======================================================================
  subroutine straighten_field(f)
    use mesh_mod
    implicit none

    type(field_type), intent(inout) :: f

    integer :: inode, numnodes, icounter_t
    vectype, dimension(dofs_per_node) :: data

    numnodes = owned_nodes()
    do icounter_t=1,numnodes
       inode = nodes_owned(icounter_t)

       call get_node_data(f, inode, data, .true.)   ! rotate to R,Z
       call set_node_data(f, inode, data, .false.)  ! don't rotate back to n,t
    enddo
    call finalize(f%vec)
  end subroutine straighten_field

  subroutine unstraighten_field(f)
   use mesh_mod
    implicit none

    type(field_type), intent(inout) :: f

    integer :: inode, numnodes, icounter_t
    vectype, dimension(dofs_per_node) :: data

    numnodes = owned_nodes()
    do icounter_t=1,numnodes
       inode = nodes_owned(icounter_t)

       call get_node_data(f, inode, data, .false.)   ! rotate to R,Z
       call set_node_data(f, inode, data, .true.)  ! don't rotate back to n,t
    enddo
    call finalize(f%vec)
  end subroutine unstraighten_field
end module field
