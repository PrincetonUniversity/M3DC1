! This module provides routines for defining the stellarator geometry,
! specifically, the logical to physical coordinate transformation.
module geometry
  use mesh_mod
  use basic
  use matrix_mod
  use field
  use arrays 
  implicit none

#ifdef USEST
 
contains

  subroutine calc_geometry
    use sparse
    use m3dc1_nint
    use boundary_conditions

    implicit none
    type(matrix_type) :: st_matrix 
    integer :: itri, numelms, ibound, ier
    integer :: inode, numnodes, ii
    vectype, dimension(dofs_per_element) :: dofs
    vectype, dimension(dofs_per_element, dofs_per_element) :: mat_dofs

    ! Define coordinate mappings to be solved for 
    call create_field(rst)
    rst = 0.
    call create_field(zst)
    zst = 0.

    numelms = local_elements()

    ! Create a matrix for solving geometry 
    call set_matrix_index(st_matrix, st_mat_index)
    call create_mat(st_matrix, 1, 1, icomplex, 1)

    select case(igeometry)

    case(1) ! prescribed geometry when igeometry==1

    do itri=1,numelms

      call define_element_quadrature(itri,int_pts_main,int_pts_tor)
      call define_fields(itri,0,1,0,ilog=1)

      call prescribe_geometry(rst79(:,OP_1), zst79(:,OP_1), x_79, phi_79, z_79)

      dofs = intx2(must79(:,:,OP_1),rst79(:,OP_1))
      call vector_insert_block(rst%vec, itri, 1, dofs, VEC_ADD)

      dofs = intx2(must79(:,:,OP_1),zst79(:,OP_1))
      call vector_insert_block(zst%vec, itri, 1, dofs, VEC_ADD)

      mat_dofs = intxx2(mu79(:,:,OP_1),nu79(:,:,OP_1))
      call insert_block(st_matrix, itri, 1, 1, mat_dofs, MAT_ADD)

    enddo

!    case(2) ! solve geometry from Laplace's equation

!    ibound = BOUNDARY_DIRICHLET

!    do itri=1,numelms

!      call define_element_quadrature(itri,int_pts_main,int_pts_tor)
!      call define_fields(itri,0,1,0,ilog=1)

!      mat_dofs = -intxx3(must79(:,:,OP_1),nust79(:,:,OP_LP),ri_79) 
!                !+intxx3(must79(:,:,OP_1),nust79(:,:,OP_DPP),ri3_79)
       
!      call apply_boundary_mask(itri, ibound, mat_dofs, tags=inner_wall)

!      call insert_block(st_matrix, itri, 1, 1, mat_dofs, MAT_ADD)
!    enddo

!    call boundary_geometry(rst%vec, zst%vec, st_matrix)

    case default ! identity geometry

    do itri=1,numelms

      call define_element_quadrature(itri,int_pts_main,int_pts_tor)
      call define_fields(itri,0,1,0,ilog=1)

      rst79(:,OP_1) = x_79
      zst79(:,OP_1) = z_79

      dofs = intx2(must79(:,:,OP_1),rst79(:,OP_1))
      call vector_insert_block(rst%vec, itri, 1, dofs, VEC_ADD)

      dofs = intx2(must79(:,:,OP_1),zst79(:,OP_1))
      call vector_insert_block(zst%vec, itri, 1, dofs, VEC_ADD)

      mat_dofs = intxx2(mu79(:,:,OP_1),nu79(:,:,OP_1))
      call insert_block(st_matrix, itri, 1, 1, mat_dofs, MAT_ADD)

    enddo
    end select

    call finalize(st_matrix)

    call sum_shared(rst%vec)
    call sum_shared(zst%vec)
    if(myrank.eq.0 .and. iprint.ge.2) print *, "calculating geometry"
    call newsolve(st_matrix,rst%vec,ier)
    call newsolve(st_matrix,zst%vec,ier)
    if(myrank.eq.0 .and. iprint.ge.2) print *, "geometry calculated"

    call destroy_mat(st_matrix)

    if (igeometry.eq.1) then
      numnodes = owned_nodes()
      allocate(rnode(dofs_per_node, numnodes))
      allocate(znode(dofs_per_node, numnodes))
      do inode=1,numnodes
        ii=nodes_owned(inode) 
        call get_node_data(rst, ii, rnode(:,inode))
        call get_node_data(zst, ii, znode(:,inode))
      end do
    end if

  end subroutine calc_geometry

  subroutine destroy_geometry
    implicit none

    call destroy_field(rst)
    call destroy_field(zst)
    if (igeometry.eq.1) then
      deallocate(rnode)
      deallocate(znode)
    end if 
  end subroutine destroy_geometry

  ! set Dirichlet BC for solving Laplace equation 
  subroutine boundary_geometry(rst, zst, mat)
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
       call get_boundary_geometry(rtemp,ztemp,x,phi,z) 
       call set_dirichlet_bc(ibegin,rst,rtemp,normal,curv,izonedim,mat)
       ! do not pass mat on the second call 
       call set_dirichlet_bc(ibegin,zst,ztemp,normal,curv,izonedim) 
    end do
  end subroutine boundary_geometry

  ! specify geometry of boundary
  subroutine get_boundary_geometry(rout, zout, x, phi, z)
    implicit none

    real, intent(in) :: x, phi, z
    vectype, intent(out), dimension(dofs_per_node) :: rout, zout 
    real :: theta

    theta = atan2(z - zzero, x - xzero)
    rout(1) = x 
    rout(2) = 1. 
    zout(1) = z 
    zout(3) = 1.

  end subroutine get_boundary_geometry

  ! Calculate rst and zst given x, phi, z 
  subroutine prescribe_geometry(rout, zout, x, phi, z)
    implicit none

    real, intent(in), dimension(MAX_PTS) :: x, phi, z
    vectype, intent(out), dimension(MAX_PTS) :: rout, zout 
    !real, dimension(MAX_PTS) :: r, theta
    integer :: i
    
    do i=1, MAX_PTS
      call physical_geometry(rout(i), zout(i), x(i), phi(i), z(i))
    enddo
  end subroutine prescribe_geometry
  
#endif
end module geometry 
