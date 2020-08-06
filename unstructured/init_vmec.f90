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

    case(1) 

    do itri=1,numelms

      call define_element_quadrature(itri,int_pts_main,int_pts_tor)
      call define_fields(itri,0,1,0)

      call vmec_pressure(temp79a, xl_79, phi_79, zl_79)

      dofs = intx2(mu79(:,:,OP_1),temp79a)
      call vector_insert_block(p_vec%vec, itri, 1, dofs, VEC_ADD)

    enddo

!    case(2) 

    case default 

    end select

    call newvar_solve(p_vec%vec,mass_mat_lhs)
    p_field(0) = p_vec
    call destroy_field(p_vec)

  end subroutine vmec_init

  ! Calculate pressure given x, phi, z 
  elemental subroutine vmec_pressure(pout, x, phi, z)
    implicit none

    real, intent(in) :: x, phi, z
    vectype, intent(out) :: pout
    real :: r, r2n, ds 
    integer :: js 
    
    r = sqrt((x - xcenter)**2 + (z - zcenter)**2 + 0e-6)
    pout = 0
    r2n = r**2*(ns-1)
    js = ceiling(r2n)
    if (js>(ns-1)) js = ns-1 
    ds = js - r2n 
    pout = presf(js+1)*(1-ds) + presf(js)*ds
  end subroutine vmec_pressure

  ! Calculate VMEC fields given x, phi, z 
  elemental subroutine vmec_fields(x, phi, z, br, bphi, bz, p)
    implicit none

    real, intent(in) :: x, phi, z
    real, intent(out) :: p, br, bphi, bz
    real :: r, r2n, ds, rout, bu, bv, theta 
    integer :: js, i 
    real, dimension(mn_mode) :: rstc, zsts, co, sn, buc, bvc 
    real :: dr, dz, dr1, dz1, phis

    phis = phi*mf+mesh_phase
    
    r = sqrt((x - xcenter)**2 + (z - zcenter)**2 + 0e-6)
    theta = atan2(z - zcenter, x - xcenter)
    p = 0
    r2n = r**2*(ns-1)
    js = ceiling(r2n)
    if (js>(ns-1)) js = ns-1 
    ds = js - r2n 
    co = cos(xmv*theta+xnv*phis)
    sn = sin(xmv*theta+xnv*phis)
!    p = presf(js+1)*(1-ds) + presf(js)*ds
    call evaluate_spline(presf_spline, r**2, p)
    call zernike_evaluate(r,rmncz,rstc)
    call zernike_evaluate(r,zmnsz,zsts)
    call zernike_evaluate(r,bsupumncz,buc)
    call zernike_evaluate(r,bsupvmncz,bvc)
!    rstc = rmnc(:,js+1)*(1-ds) + rmnc(:,js)*ds
!    zsts = zmns(:,js+1)*(1-ds) + zmns(:,js)*ds
!    buc = bsupumnc(:,js+1)*(1-ds) + bsupumnc(:,js)*ds
!    bvc = bsupvmnc(:,js+1)*(1-ds) + bsupvmnc(:,js)*ds
    rout = 0.
    dr = 0.
    dz = 0.
    dr1 = 0.
    dz1 = 0.
    bu = 0.
    bv = 0.
    do i = 1, mn_mode 
      rout = rout + rstc(i)*co(i)
      dr = dr - rstc(i)*sn(i)*xmv(i)
      dz = dz + zsts(i)*co(i)*xmv(i)
      dr1 = dr - rstc(i)*sn(i)*xnv(i)*mf
      dz1 = dz + zsts(i)*co(i)*xnv(i)*mf
      bu = bu + buc(i)*co(i) 
      bv = bv + bvc(i)*co(i) 
    end do
    br = bu*dr + bv*dr1    
    bphi = rout*bv
    bz = bu*dz + bv*dz1
  end subroutine vmec_fields
  
#endif
end module init_vmec 
