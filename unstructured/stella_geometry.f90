! This module provides routines for defining the stellarator geometry,
! specifically, the logical to physical coordinate transformation.
module stella_gemoetry
  implicit none

contains

  subroutine calc_stella_geometry
    use basic
    use arrays
    use sparse
    use mesh_mod
    use m3dc1_nint
    use matrix_mod
    use field
    use boundary_conditions

    implicit none
    type(matrix_type) :: lp_matrix
    integer :: itri, numelms, ibound, ier
    vectype, dimension(dofs_per_element, dofs_per_element) :: mat_dofs

    ! Create a matrix for solving -del*(A)/r = B
    call set_matrix_index(lp_matrix, lp_mat_index)
    call create_mat(lp_matrix, 1, 1, icomplex, 1)

    ! Define coordinate mappings to be solved for 
    call create_field(Rs)
    Rs = 0.
    call create_field(Zs)
    Zs = 0.

    ibound = BOUNDARY_DIRICHLET

    numelms = local_elements()
    do itri=1,numelms

       call define_element_quadrature(itri,int_pts_main,int_pts_tor)
       call define_fields(itri,0,1,0)

       mat_dofs = intxx3(mu79(:,:,OP_1),nu79(:,:,OP_GS),ri_79)
       
       call apply_boundary_mask(itri, ibound, mat_dofs, tags=domain_boundary)

       call insert_block(lp_matrix, itri, 1, 1, mat_dofs, MAT_ADD)
    enddo

    call boundary_stella_geomatry(Rs%vec, Zs%vec, lp_matrix)
    call finalize(lp_matrix)
    
    call newsolve(lp_matrix,Rs%vec,ier)
    call newsolve(lp_matrix,Zs%vec,ier)

    call destroy_mat(lp_matrix)

  end subroutine calc_stella_geomatry

  subroutine destroy_stella_geomatry
    use arrays

    implicit none

    call destroy_field(Rs)
    call destroy_field(Zs)
  end subroutine destroy_stella_geomatry


  subroutine boundary_stella_geomatry(rhs, mat)
    use basic
    use field
    use arrays
    use matrix_mod
    use boundary_conditions
    implicit none
    
    type(vector_type) :: rhs
    type(matrix_type), optional :: mat
    
    integer :: i, izone, izonedim, numnodes, icounter_t
    real :: normal(2), curv, x, z, phi
    logical :: is_boundary
    vectype, dimension(dofs_per_node) :: temp
    
    integer :: index
    
    if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_stella_geomatry called"
    
    numnodes = owned_nodes()
    do icounter_t=1,numnodes
       i = nodes_owned(icounter_t)
       index = node_index(rhs, i, 1)
       
       call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,phi,z, &
            inner_wall)
       if(.not.is_boundary) cycle
       
       temp = 0.
       temp(2) = -normal(1)
       temp(3) = -normal(2)
       call set_dirichlet_bc(index,rhs,temp,normal,curv,izonedim,mat)
       call set_normal_bc(index,rhs,temp,normal,curv,izonedim,mat)
    end do
  end subroutine boundary_stella_geomatry
  
end module wall
