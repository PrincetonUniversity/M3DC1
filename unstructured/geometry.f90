! This module provides routines for defining the stellarator geometry,
! specifically, the logical to physical coordinate transformation.
module geometry
  implicit none

contains

  subroutine calc_geometry
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

    ! Create a matrix for solving Laplace Eq. 
    call set_matrix_index(lp_matrix, lp_mat_index)
    call create_mat(lp_matrix, 1, 1, icomplex, 1)

    ! Define coordinate mappings to be solved for 
    call create_field(Rst)
    Rst = 0.
    call create_field(Zst)
    Zst = 0.

    ibound = BOUNDARY_DIRICHLET

    numelms = local_elements()
    do itri=1,numelms

       call define_element_quadrature(itri,int_pts_main,int_pts_tor)
       call define_fields(itri,0,1,0)

       mat_dofs = intxx3(mu79(:,:,OP_1),nu79(:,:,OP_LP),ri_79)\
                 +intxx3(mu79(:,:,OP_1),nu79(:,:,OP_DPP),ri3_79)
       
       call apply_boundary_mask(itri, ibound, mat_dofs, tags=domain_boundary)

       call insert_block(lp_matrix, itri, 1, 1, mat_dofs, MAT_ADD)
    enddo

    call boundary_geometry(Rst%vec, Zst%vec, lp_matrix)
    call finalize(lp_matrix)
    
    call newsolve(lp_matrix,Rst%vec,ier)
    call newsolve(lp_matrix,Zst%vec,ier)

    call destroy_mat(lp_matrix)

  end subroutine calc_geometry

  subroutine destroy_geometry
    use arrays

    implicit none

    call destroy_field(Rst)
    call destroy_field(Zst)
  end subroutine destroy_geometry

  ! set Dirichlet BC for solving Laplace Eq.
  subroutine boundary_geometry(Rst, Zst, mat)
    use basic
    use field
    use arrays
    use matrix_mod
    use boundary_conditions
    implicit none
    
    type(vector_type) :: Rst, Zst
    type(matrix_type), optional :: mat
    
    integer :: i, izone, izonedim, numnodes, icounter_t
    real :: normal(2), curv, x, z, phi
    logical :: is_boundary
    vectype, dimension(dofs_per_node) :: Rtemp, Ztemp
    
    integer :: index
    
    if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_geometry called"
    
    numnodes = owned_nodes()
    do icounter_t=1,numnodes
       i = nodes_owned(icounter_t)
       index = node_index(Rst, i, 1)
       
       call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,phi,z, &
            inner_wall)
       if(.not.is_boundary) cycle
       
       Rtemp = 0.
       Ztemp = 0.
       call get_geometry(Rtemp(1),Ztemp(1),x,phi,z) 
       call set_dirichlet_bc(index,Rst,Rtemp,normal,curv,izonedim,mat)
       call set_dirichlet_bc(index,Zst,Ztemp,normal,curv,izonedim,mat)
    end do
  end subroutine boundary_geometry

  ! Calculate Rst and Zst given x, phi, z 
  subroutine get_geometry(Rout, Zout, x, phi, z)
    use basic
    use field
    use arrays
    use matrix_mod
    use boundary_conditions
    implicit none

    real, intent(in) :: x, phi, z
    real, intent(out) :: Rout, Zout 
    real :: r, theta

    r = sqrt((x - xzero)**2 + (z - zzero)**2 + regular**2)
    theta = -atan2(z - zzero, x - xzero)
    Rout = 3.0 + 0.3*cos(theta) - 0.06*cos(theta-2*phi)
    Zout = 0.3*sin(theta) +0.06*sin(2*phi) - 0.06*sin(theta-2*phi)

  end subroutine get_geometry
  
end module geometry 
