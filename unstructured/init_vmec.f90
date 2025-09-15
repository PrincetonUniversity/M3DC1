! This module sets up initial conditions using VMEC data
module init_vmec 
  use mesh_mod
  use basic
  use matrix_mod
  use field
  use arrays 
  use math 
  use read_vmec 
  use spline 
  implicit none

  type(spline1d), private :: den_spline        ! density 
  type(spline1d), private :: temper_spline        ! temperature
  type(spline1d), private :: press_ext_spline        ! pressure
  type(spline1d), private :: presf_spline      ! total pressure
#ifdef USEST
 
contains

  subroutine vmec_init()
    use sparse
    use newvar_mod 
    use m3dc1_nint
    use boundary_conditions
    use init_common
    use read_ascii

    implicit none


    type(matrix_type) :: br_mat
    type(vector_type) :: fppsi_vec
    type(field_type) :: psi_f, bf_f, bfp_f, bz_f, den_f, te_f, pexternal_f
    type(field_type) :: p_f 
    integer :: itri, numelms, ifpbound, ier, ipsibound, ipsifpbound, i, k, k1
    integer :: inode(nodes_per_element)
    vectype, dimension(dofs_per_element) :: dofs
    vectype, dimension(dofs_per_element,dofs_per_element,2,2) :: temp
    vectype, dimension(dofs_per_element,2) :: temp2
    real :: fzero
    real, allocatable :: xvals(:), yvals(:)
    integer :: nvals

    ! normalize pressure
    presf = presf*pi*40/b0_norm**2
    ! 1D spline for pressure
    call create_spline(presf_spline, ns, s_vmec, presf)

    if(itor.eq.0) then 
      fzero = bzero
    else 
      fzero = bzero*rzero
    end if

    ! Create fields 
    call create_field(p_f)
    call create_field(bz_f)
    call create_field(bf_f)
    p_f = 0.
    bz_f = 0.
    bf_f = 0.
    if(iread_ne.eq.21) then 
      call create_field(den_f)
      den_f = 0.
      nvals = 0
      call read_ascii_column('n_profile', xvals, nvals, icol=1)
      call read_ascii_column('n_profile', yvals, nvals, icol=2)
      if(nvals.eq.0) call safestop(5)
      yvals = yvals / 1e6 / n0_norm / z_ion
      call create_spline(den_spline, nvals, xvals, yvals)
      deallocate(xvals, yvals)
    endif

    if(iread_te.eq.21) then 
      call create_field(te_f)
      te_f = 0.
      nvals = 0
      call read_ascii_column('te_profile', xvals, nvals, icol=1)
      call read_ascii_column('te_profile', yvals, nvals, icol=2)
      if(nvals.eq.0) call safestop(5)
      yvals = yvals * 1.6022e-9 / (b0_norm**2/(4.*pi*n0_norm))
      call create_spline(temper_spline , nvals, xvals, yvals)
      deallocate(xvals, yvals)
    endif

    if(iread_p.eq.21) then 
      call create_field(pexternal_f)
      pexternal_f = 0.
      nvals = 0
      call read_ascii_column('p_profile', xvals, nvals, icol=1)
      call read_ascii_column('p_profile', yvals, nvals, icol=2)
      if(nvals.eq.0) call safestop(5)
      yvals = yvals * 10. / (b0_norm**2/(4.*pi))
      call create_spline(press_ext_spline , nvals, xvals, yvals)
      deallocate(xvals, yvals)
    endif

    call create_vector(fppsi_vec,2)
    call associate_field(bfp_f,fppsi_vec,1)
    call associate_field(psi_f,fppsi_vec,2)

    call set_matrix_index(br_mat, br_mat_index)
    call create_mat(br_mat, 2, 2, icomplex, 1)

    ! boundary condition on psi, g, and f
    !ipsifpbound = BOUNDARY_NONE
    !ipsibound = BOUNDARY_DIRICHLET
    ipsibound = BOUNDARY_NONE
    !ifpbound = BOUNDARY_NONE
    ifpbound = BOUNDARY_DIRICHLET
    !ifpbound = BOUNDARY_NEUMANN

    numelms = local_elements()

    do itri=1,numelms

      call define_element_quadrature(itri,int_pts_main,int_pts_tor)
      call define_fields(itri,0,1,0)
      !vmec_fields(x, phi, z, br, bphi, bz, p, den,temper)
      do k=1,MAX_PTS
       call vmec_fields(xl_79(k), phi_79(k), zl_79(k), temp79a(k), &
            temp79b(k), temp79c(k), temp79d(k), temp79e(k), temp79f(k))
      enddo
      ! fp equation
      temp(:,:,1,1) = &
          -intxx2(mu79(:,:,OP_DR),nu79(:,:,OP_DR)) &
          -intxx2(mu79(:,:,OP_DZ),nu79(:,:,OP_DZ)) &
          + regular*intxx3(mu79(:,:,OP_1),nu79(:,:,OP_1),ri2_79)
      temp(:,:,1,2) = intxx3(mu79(:,:,OP_DZ),nu79(:,:,OP_DR),ri_79) &
                    - intxx3(mu79(:,:,OP_DR),nu79(:,:,OP_DZ),ri_79)
      temp2(:,1) = intx2(mu79(:,:,OP_DR),temp79a) &
                 + intx2(mu79(:,:,OP_DZ),temp79c) 

      ! psi equation
      temp(:,:,2,1) = -intxx3(mu79(:,:,OP_DZ),nu79(:,:,OP_DR),ri_79) &
                      +intxx3(mu79(:,:,OP_DR),nu79(:,:,OP_DZ),ri_79) 
      temp(:,:,2,2) = -intxx3(mu79(:,:,OP_DZ),nu79(:,:,OP_DZ),ri2_79) &
                      -intxx3(mu79(:,:,OP_DR),nu79(:,:,OP_DR),ri2_79) &
                      -regular*intxx3(mu79(:,:,OP_1),nu79(:,:,OP_1),ri4_79) 
      temp2(:,2) = intx3(mu79(:,:,OP_DZ),temp79a,ri_79) &
                 - intx3(mu79(:,:,OP_DR),temp79c,ri_79) 

      call apply_boundary_mask(itri, ifpbound, temp(:,:,1,1), &
          tags=BOUND_DOMAIN)
      call apply_boundary_mask(itri, ifpbound, temp(:,:,1,2), &
          tags=BOUND_DOMAIN)
      call apply_boundary_mask(itri, ipsibound, temp(:,:,2,1), &
          tags=BOUND_DOMAIN)
      call apply_boundary_mask(itri, ipsibound, temp(:,:,2,2), &
          tags=BOUND_DOMAIN)

      call insert_block(br_mat, itri, 1, 1, temp(:,:,1,1), MAT_ADD)
      call insert_block(br_mat, itri, 1, 2, temp(:,:,1,2), MAT_ADD)
      call insert_block(br_mat, itri, 2, 1, temp(:,:,2,1), MAT_ADD)
      call insert_block(br_mat, itri, 2, 2, temp(:,:,2,2), MAT_ADD)

      call vector_insert_block(fppsi_vec, itri, 1, temp2(:,1), MAT_ADD)
      call vector_insert_block(fppsi_vec, itri, 2, temp2(:,2), MAT_ADD)

      ! pressure p
      if(iread_p.eq.21) then 
        dofs = intx2(mu79(:,:,OP_1),temp79d)
        call vector_insert_block(pexternal_f%vec, itri, 1, dofs, VEC_ADD)
      else
        dofs = intx2(mu79(:,:,OP_1),temp79d)
        call vector_insert_block(p_f%vec, itri, 1, dofs, VEC_ADD)
      end if

      ! density n 
      if(iread_ne.eq.21) then 
        dofs = intx2(mu79(:,:,OP_1),temp79e)
        call vector_insert_block(den_f%vec, itri, 1, dofs, VEC_ADD)
      end if

      ! temperature te ti 
      if(iread_te.eq.21) then 
        dofs = intx2(mu79(:,:,OP_1),temp79f)
        call vector_insert_block(te_f%vec, itri, 1, dofs, VEC_ADD)
      end if


      ! F = R*B_phi 
      dofs = intx3(mu79(:,:,OP_1),temp79b,r_79) 
      call vector_insert_block(bz_f%vec, itri, 1, dofs, VEC_ADD)
      call vector_insert_block(bf_f%vec, itri, 1, dofs, VEC_ADD)
    end do

    if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving fp & psi..."
    call sum_shared(fppsi_vec)
    if(iread_p.eq.21) then 
        call boundary_vmec(fppsi_vec,br_mat,pexternal_f)
    else
        call boundary_vmec(fppsi_vec,br_mat,p_f)
    end if
    call finalize(br_mat)
    call newsolve(br_mat,fppsi_vec,ier)
    if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving psi: ier = ", ier

    bfp_field(0) = bfp_f 
    psi_field(0) = psi_f

    if(iread_ne.eq.21) then
      if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving den..."
      call newvar_solve(den_f%vec,mass_mat_lhs)
      den_field(0) = den_f 
    end if
    if(iread_te.eq.21) then
      if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving temp..."
      call newvar_solve(te_f%vec,mass_mat_lhs)
      te_field(0) = te_f 
      ti_field(0) = te_f 
    end if

    
    if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving p & per..."
    if(iread_p.eq.21) then
     call newvar_solve(pexternal_f%vec,mass_mat_lhs)
     p_field(0)=pexternal_f
    else
     call newvar_solve(p_f%vec,mass_mat_lhs)

     p_field(0) = p_f 
    endif

    pe_field(0) = p_field(0)
    call mult(pe_field(0),pefac)
    call evaluate_spline(presf_spline, 0., p0)
    p0 = p0 + pedge 
    pi0 = p0*(1-pefac)

    if(myrank.eq.0 .and. iprint.ge.2) print *, "Solving bz & bf..."
    call newvar_solve(bz_f%vec,mass_mat_lhs)
    call newvar_solve(bf_f%vec,bf_mat_lhs)

    bz_field(0) = bz_f
    bf_field(0) = bf_f

    call destroy_field(p_f)
    call destroy_field(bf_f)
    call destroy_field(bz_f)
    if(iread_ne.eq.21) then
      call destroy_field(den_f)
      call destroy_spline(den_spline)
    end if
    if(iread_te.eq.21) then
      call destroy_field(te_f)
      call destroy_spline(temper_spline)
    end if
    if(iread_p.eq.21) then
      call destroy_field(pexternal_f)
      call destroy_spline(press_ext_spline)
    end if
    call destroy_vector(fppsi_vec)
    call destroy_mat(br_mat)

    call init_perturbations

  end subroutine vmec_init

  subroutine boundary_vmec(rhs, mat, vec)
    use vector_mod
    use matrix_mod
    use boundary_conditions
  
    implicit none
  
    type(field_type) :: vec 
    type(vector_type) :: rhs
    type(matrix_type), optional :: mat
    
    integer :: i, izone, izonedim, i_fp, i_psi, numnodes, icounter_t
    real :: normal(2), curv(3)
    real :: x, z, phi
    logical :: is_boundary
    vectype, dimension(dofs_per_node) :: temp
  
    if(iper.eq.1 .and. jper.eq.1) return
    if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_vmec called"
  
    temp = 0.
  
    numnodes = owned_nodes()
    do icounter_t=1,numnodes
       i = nodes_owned(icounter_t)
       call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,phi,z)
       if(.not.is_boundary) cycle
  
       i_fp = node_index(rhs, i, 1)
       i_psi = node_index(rhs, i, 2)
  
       !call get_node_data(vec, i, temp)
       !call set_normal_bc(i_fp, rhs,temp, normal,curv,izonedim,mat)
       call set_dirichlet_bc(i_fp, rhs,temp, normal,curv,izonedim,mat)
       !call set_dirichlet_bc(i_psi, rhs,temp, normal,curv,izonedim,mat)
    end do
  end subroutine boundary_vmec

  ! Calculate VMEC fields given x, phi, z 
  elemental subroutine vmec_fields(x, phi, z, br, bphi, bz, p, den,temper)
    implicit none

    real, intent(in) :: x, phi, z
    real, intent(out) :: p, br, bphi, bz, den,temper
    real :: r, r2n, ds, rout, bu, bv, theta 
    integer :: js, i 
    real, dimension(mn_mode) :: rstc, zsts, co, sn, ls, lc, rsts, zstc 
    real, dimension(mn_mode_nyq) :: co_nyq, sn_nyq, buc, bvc, bus, bvs 
    real :: dr, dz, dr1, dz1, phis

    phis = phi*mf+mesh_phase
    
    r = sqrt((x - xcenter)**2 + (z - zcenter)**2 + regular**2)
    theta = atan2(z - zcenter, x - xcenter)
    co = cos(xmv*theta+xnv*phis)
    sn = sin(xmv*theta+xnv*phis)
    co_nyq = cos(xmv_nyq*theta+xnv_nyq*phis)
    sn_nyq = sin(xmv_nyq*theta+xnv_nyq*phis)

    if(iread_p.eq.21) then
      call evaluate_spline(press_ext_spline, r**2, p)
    else
      call evaluate_spline(presf_spline, r**2, p)
    endif
   ! if(iread_ne.eq.21) call evaluate_spline(den_spline, r**2, den)
   ! if(iread_te.eq.21) call evaluate_spline(temper_spline, r**2, temper)

    if(iread_ne.eq.21) then
      call evaluate_spline(den_spline, r**2, den)
    else if(iread_te.eq.21) then
      call evaluate_spline(temper_spline, r**2, temper)
      if(iread_p.eq.21) then
        call evaluate_spline(press_ext_spline, r**2, p)
      else
        call evaluate_spline(presf_spline, r**2, p)
      end if
      den = temper/p
    else
      den=den0
    endif

    call zernike_evaluate(r,mn_mode,mb,rmncz,rstc)
    call zernike_evaluate(r,mn_mode,mb,zmnsz,zsts)
    !call zernike_evaluate(r,mn_mode_nyq,mb_nyq,bsupumncz,buc)
    !call zernike_evaluate(r,mn_mode_nyq,mb_nyq,bsupvmncz,bvc)
    !call vmec_interpl(r,mn_mode,rmnc,rstc)
    !call vmec_interpl(r,mn_mode,zmns,zsts)
    call vmec_interpl(r,mn_mode_nyq,mb_nyq,bsupumnc,buc)
    call vmec_interpl(r,mn_mode_nyq,mb_nyq,bsupvmnc,bvc)
    if(lasym.eq.1) then 
      call zernike_evaluate(r,mn_mode,mb,rmnsz,rsts)
      call zernike_evaluate(r,mn_mode,mb,zmncz,zstc)
      call vmec_interpl(r,mn_mode_nyq,mb_nyq,bsupumns,bus)
      call vmec_interpl(r,mn_mode_nyq,mb_nyq,bsupvmns,bvs)
    end if
    rout = 0.
    dr = 0.
    dz = 0.
    dr1 = 0.
    dz1 = 0.
    bu = 0.
    bv = 0.
    do i = 1, mn_mode 
      if (xmv(i)<m_max .and. abs(xnv(i))<n_max) then
        rout = rout + rstc(i)*co(i)
        dr = dr - rstc(i)*sn(i)*xmv(i)
        dz = dz + zsts(i)*co(i)*xmv(i)
        dr1 = dr1 - rstc(i)*sn(i)*xnv(i)*mf
        dz1 = dz1 + zsts(i)*co(i)*xnv(i)*mf
        if(lasym.eq.1) then 
          rout = rout + rsts(i)*sn(i)
          dr = dr + rsts(i)*co(i)*xmv(i)
          dz = dz - zstc(i)*sn(i)*xmv(i)
          dr1 = dr1 + rsts(i)*co(i)*xnv(i)*mf
          dz1 = dz1 - zstc(i)*sn(i)*xnv(i)*mf
        end if 
      end if 
    end do
    do i = 1, mn_mode_nyq 
      if (xmv_nyq(i)<m_max .and. abs(xnv_nyq(i))<n_max) then
        bu = bu + buc(i)*co_nyq(i) 
        bv = bv + bvc(i)*co_nyq(i) 
        if(lasym.eq.1) then 
          bu = bu + bus(i)*sn_nyq(i) 
          bv = bv + bvs(i)*sn_nyq(i) 
        end if 
      end if 
    end do
    br = bu*dr + bv*dr1    
    bphi = rout*bv 
    bz = bu*dz + bv*dz1
    p = p + pedge  
  end subroutine vmec_fields
  
#endif
end module init_vmec 
