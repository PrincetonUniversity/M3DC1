! This module provides routines for defining the stellarator geometry,
! specifically, the logical to physical coordinate transformation.
module geometry
  implicit none

#ifdef USEST

contains

  subroutine calc_geometry
    use basic
    use arrays
    use sparse
    use mesh_mod
    use m3dc1_nint
    use matrix_mod
    use newvar_mod
    use field
    use boundary_conditions

    implicit none
    type(matrix_type) :: lp_matrix
    integer :: itri, numelms, ibound, ier
    real, dimension(dofs_per_element) :: dofs
    real, dimension(MAX_PTS) :: rst_79, zst_79 
    vectype, dimension(dofs_per_element, dofs_per_element) :: mat_dofs

    ! Define coordinate mappings to be solved for 
    call create_field(rst)
    rst = 0.
    call create_field(zst)
    zst = 0.

    numelms = local_elements()

    select case(igeometry)

    case(1) ! prescribed geometry
    do itri=1,numelms

       call define_element_quadrature(itri,int_pts_main,int_pts_tor)
       call define_fields(itri,0,1,0)

       call prescribe_geometry(rst_79, zst_79, x_79, phi_79, z_79)

       dofs = intx2(mu79(:,:,OP_1),rst_79)
       call vector_insert_block(rst%vec, itri, 1, dofs, VEC_ADD)

       dofs = intx2(mu79(:,:,OP_1),zst_79)
       call vector_insert_block(zst%vec, itri, 1, dofs, VEC_ADD)

    enddo

    call newvar_solve(rst%vec, mass_mat_lhs)
    call newvar_solve(zst%vec, mass_mat_lhs)

    case(2) ! solve geometry from Laplace's equation

    ! Create a matrix for solving Laplace Eq. 
    call set_matrix_index(lp_matrix, lp_mat_index)
    call create_mat(lp_matrix, 1, 1, icomplex, 1)

    ibound = BOUNDARY_DIRICHLET

    do itri=1,numelms

       call define_element_quadrature(itri,int_pts_main,int_pts_tor)
       call define_fields(itri,0,1,0)

       mat_dofs = -intxx3(mu79(:,:,OP_1),nu79(:,:,OP_LP),ri_79) 
                 !+intxx3(mu79(:,:,OP_1),nu79(:,:,OP_DPP),ri3_79)
       
       call apply_boundary_mask(itri, ibound, mat_dofs, tags=inner_wall)

       call insert_block(lp_matrix, itri, 1, 1, mat_dofs, MAT_ADD)
    enddo

    call boundary_geometry(rst%vec, zst%vec, lp_matrix)
    call finalize(lp_matrix)
   
    call newsolve(lp_matrix,rst%vec,ier)
    call newsolve(lp_matrix,zst%vec,ier)

    call destroy_mat(lp_matrix)

    case default ! identity geometry
    do itri=1,numelms

       call define_element_quadrature(itri,int_pts_main,int_pts_tor)
       call define_fields(itri,0,1,0)

       rst_79 = x_79
       zst_79 = z_79

       dofs = intx2(mu79(:,:,OP_1),rst_79)
       call vector_insert_block(rst%vec, itri, 1, dofs, VEC_ADD)

       dofs = intx2(mu79(:,:,OP_1),zst_79)
       call vector_insert_block(zst%vec, itri, 1, dofs, VEC_ADD)

    enddo

    call newvar_solve(rst%vec, mass_mat_lhs)
    call newvar_solve(zst%vec, mass_mat_lhs)

    end select

  end subroutine calc_geometry

  subroutine destroy_geometry
    use arrays

    implicit none

    call destroy_field(rst)
    call destroy_field(zst)
  end subroutine destroy_geometry

  ! set Dirichlet BC for solving Laplace equation 
  subroutine boundary_geometry(rst, zst, mat)
    use basic
    use field
    use arrays
    use matrix_mod
    use boundary_conditions
    implicit none
    
    type(vector_type) :: rst, zst
    type(matrix_type), optional :: mat
    
    integer :: i, izone, izonedim, numnodes, icounter_t
    real :: normal(2), curv, x, z, phi
    logical :: is_boundary
    vectype, dimension(dofs_per_node) :: rtemp, ztemp
    
    integer :: ibegin
    
    if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_geometry called"
    
    numnodes = owned_nodes()
    do icounter_t=1,numnodes
       i = nodes_owned(icounter_t)
       
       call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,phi,z, &
            inner_wall)
       if(.not.is_boundary) cycle
       
       ibegin = node_index(rst, i, 1)
       rtemp = 0.
       ztemp = 0.
       call get_boundary_geometry(rtemp(1),ztemp(1),x,phi,z) 
       call set_dirichlet_bc(ibegin,rst,rtemp,normal,curv,izonedim,mat)
       call set_dirichlet_bc(ibegin,zst,ztemp,normal,curv,izonedim)
    end do
  end subroutine boundary_geometry

  ! specify geometry of boundary
  subroutine get_boundary_geometry(rout, zout, x, phi, z)
    use basic
    use field
    use arrays
    use matrix_mod
    use boundary_conditions
    implicit none

    real, intent(in) :: x, phi, z
    real, intent(out) :: rout, zout 
    real :: theta

    theta = -atan2(z - zzero, x - xzero)
    rout = 1.0 + 0.4*cos(theta) 
    zout = 0.8*sin(theta)

  end subroutine get_boundary_geometry

  ! Calculate rst and zst given x, phi, z 
  subroutine prescribe_geometry(rout, zout, x, phi, z)
    use basic
    use field
    use arrays
    use matrix_mod
    use boundary_conditions
    implicit none

    real, intent(in), dimension(MAX_PTS) :: x, phi, z
    real, intent(out), dimension(MAX_PTS) :: rout, zout 
    real, dimension(MAX_PTS) :: r, theta

    r = sqrt((x - xzero)**2 + (z - zzero)**2 + regular**2)
    theta = -atan2(z - zzero, x - xzero)
    rout = 1.0 + 0.6*r*cos(theta) 
    zout = 0.4*r*sin(theta) 

  end subroutine prescribe_geometry
  
#endif
end module geometry 
