! This module sets up initial conditions using VMEC data
module init_vmec 
  use mesh_mod
  use basic
  use matrix_mod
  use field
  use arrays 
  use read_vmec 
  implicit none

#ifdef USEST
 
contains

  subroutine vmec_init()
    use sparse
    use newvar_mod 
    use m3dc1_nint
    use boundary_conditions

    implicit none
    type(field_type) :: p_vec 
    integer :: itri, numelms, ibound, ier
    vectype, dimension(dofs_per_element) :: dofs
!    vectype, dimension(dofs_per_element, dofs_per_element) :: mat_dofs

    ! Create fields 
    call create_field(p_vec)
    p_vec = 0.

    numelms = local_elements()

    ! Create a matrix for solving geometry 
!    call set_matrix_index(st_matrix, st_mat_index)
!    call create_mat(st_matrix, 1, 1, icomplex, 1)

    if(myrank.eq.0 .and. iprint.ge.1) print *, 'Defining VMEC Equilibrium'

    select case(igeometry)

    case(1) ! prescribed geometry when igeometry==1

    do itri=1,numelms

      call define_element_quadrature(itri,int_pts_main,int_pts_tor)
      call define_fields(itri,0,1,0)

      call vmec_pressure(temp79a, xl_79, phi_79, zl_79)

      dofs = intx2(mu79(:,:,OP_1),temp79a)
      call vector_insert_block(p_vec%vec, itri, 1, dofs, VEC_ADD)

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

!    case default ! identity geometry

!    do itri=1,numelms

!      call define_element_quadrature(itri,int_pts_main,int_pts_tor)
!      call define_fields(itri,0,1,0,ilog=1)

!      rst79(:,OP_1) = x_79
!      zst79(:,OP_1) = z_79

!      dofs = intx2(must79(:,:,OP_1),rst79(:,OP_1))
!      call vector_insert_block(rst%vec, itri, 1, dofs, VEC_ADD)

!      dofs = intx2(must79(:,:,OP_1),zst79(:,OP_1))
!      call vector_insert_block(zst%vec, itri, 1, dofs, VEC_ADD)

!      mat_dofs = intxx2(mu79(:,:,OP_1),nu79(:,:,OP_1))
!      call insert_block(st_matrix, itri, 1, 1, mat_dofs, MAT_ADD)

!    enddo
    end select

    call newvar_solve(p_vec%vec,mass_mat_lhs)
    p_field(0) = p_vec
    call destroy_field(p_vec)

!    call finalize(st_matrix)

!    call sum_shared(rst%vec)
!    call sum_shared(zst%vec)
!    if(myrank.eq.0 .and. iprint.ge.2) print *, "calculating geometry"
!    call newsolve(st_matrix,rst%vec,ier)
!    call newsolve(st_matrix,zst%vec,ier)
!    if(myrank.eq.0 .and. iprint.ge.2) print *, "geometry calculated"

!    call destroy_mat(st_matrix)

  end subroutine vmec_init

  ! Calculate pressure given x, phi, z 
  elemental subroutine vmec_pressure(pout, x, phi, z)
    implicit none

    real, intent(in) :: x, phi, z
    vectype, intent(out) :: pout
    real :: r, r2n, ds 
    integer :: js 
    
    r = sqrt((x - xzero)**2 + (z - zzero)**2 + 0e-6)
!    theta = atan2(z - zzero, x - xzero)
    pout = 0
    r2n = r**2*(ns-1)
    js = ceiling(r2n)
    if (js>(ns-1)) js = ns-1 
    ds = js - r2n 
    pout = presf(js+1)*(1-ds) + presf(js)*ds
  end subroutine vmec_pressure
  
#endif
end module init_vmec 
