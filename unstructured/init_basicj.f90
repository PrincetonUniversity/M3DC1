!===============================
! Basicj Equilibrium
! ~~~~~~~~~~~~~~~~~~
!===============================

module basicj

  implicit none

  real :: basicj_nu
  real :: basicj_j0
  real :: basicj_voff
  real :: basicj_vdelt
  real :: basicj_dexp
  real :: basicj_dvac

contains

  subroutine basicj_current(x, z, jphi)
    use basic

    implicit none

    real, intent(in), dimension(MAX_PTS) :: x, z
    vectype, intent(out), dimension(MAX_PTS) :: jphi
    
    real, dimension(MAX_PTS) :: r2

    r2 = (x - xmag)**2 + (z - zmag)**2
    where(r2.lt.ln**2)
       jphi = basicj_j0*(1.-(r2/ln**2))**basicj_nu
    elsewhere
       jphi = 0.
    end where
  end subroutine basicj_current

  subroutine basicj_vz(x, z, vz)
    use basic

    implicit none

    real, intent(in), dimension(MAX_PTS) :: x, z
    vectype, intent(out), dimension(MAX_PTS) :: vz
    
    real, dimension(MAX_PTS) :: r

    r = sqrt((x - xmag)**2 + (z - zmag)**2)
    where(r.lt.basicj_voff)
       vz = vzero
    elsewhere
       vz = vzero*exp(-(r - basicj_voff)**2/(basicj_vdelt*ln)**2) &
            + alpha1*exp(-(r**2 - basicj_voff**2)**2/(basicj_vdelt*ln)**4)
    end where
    where(r.lt.ln)
       vz = vz + alpha0*(1.-(r/ln)**2)
    end where
  end subroutine basicj_vz


  subroutine basicj_init()
    use basic
    use arrays
    use m3dc1_nint
    use mesh_mod
    use newvar_mod
    use field
    use diagnostics
    use sparse
    use boundary_conditions
    use matrix_mod
    use init_common

    implicit none

    type(field_type) :: jphi_vec, f_vec, vz_vec
    integer :: itri, numelms, i, j, ibound, ierr
    integer, dimension(dofs_per_element) :: imask
    vectype, dimension(dofs_per_element) :: dofs, dofs_vz
    vectype, dimension(dofs_per_element,dofs_per_element) :: temp
    type(matrix_type) :: dr_matrix, lp_matrix

    if(myrank.eq.0 .and. iprint.ge.1) print *, 'Defining BASICJ Equilibrium'

    call create_field(jphi_vec)

    jphi_vec = 0.

    ! Create a matrix for solving -del*(A)/r = B
    call set_matrix_index(lp_matrix, lp_mat_index)
    call create_mat(lp_matrix, 1, 1, icomplex, 1)

    ibound = BOUNDARY_DIRICHLET

    ! Define Jphi
    numelms = local_elements()
    do itri=1,numelms

       call define_element_quadrature(itri,int_pts_main,int_pts_tor)
       call define_fields(itri,0,1,0)

       call basicj_current(x_79, z_79, temp79a)
              
       call get_boundary_mask(itri, ibound, imask, domain_boundary)
       do i=1,dofs_per_element
          if(imask(i).eq.0) then
             temp(i,:) = 0.
             dofs(i) = 0.
          else
             do j=1,dofs_per_element
                temp(i,j) = -int3(ri_79,mu79(:,OP_1,i),nu79(:,OP_GS,j))
             enddo
             
             dofs(i) = int2(mu79(:,OP_1,i),temp79a)
          endif
       enddo
       
       call insert_block(lp_matrix, itri, 1, 1, temp, MAT_ADD)
       call vector_insert_block(jphi_vec%vec, itri, 1, dofs, MAT_ADD)
    enddo
    
    call sum_shared(jphi_vec%vec)
    call boundary_dc(jphi_vec%vec, jphi_vec%vec, lp_matrix)
    call finalize(lp_matrix)

    if(myrank.eq.0 .and. iprint.ge.2) print *, 'Solving PSI'
    call newsolve(lp_matrix,jphi_vec%vec,ierr)
    psi_field(0) = jphi_vec

    call destroy_mat(lp_matrix)


    ! Create a matrix for solving (dA/dpsi)/(R) = Jphi
    call set_matrix_index(dr_matrix, dr_mat_index)
    call create_mat(dr_matrix, 1, 1, icomplex, 1)
    
    jphi_vec = 0.
    do itri=1,numelms

       call define_element_quadrature(itri,int_pts_main,int_pts_tor)
       call define_fields(itri,0,1,0)

       call eval_ops(itri, psi_field(0), ps079, rfac)

       call basicj_current(x_79, z_79, temp79a)
       temp79b = ps079(:,OP_DR)**2 + ps079(:,OP_DZ)**2
       
       call get_boundary_mask(itri, ibound, imask, domain_boundary)
       do i=1,dofs_per_element
          if(imask(i).eq.0) then
             temp(i,:) = 0.
             dofs(i) = 0.
          else
             do j=1,dofs_per_element
                temp(i,j) = int4(ri_79,mu79(:,OP_1,i),nu79(:,OP_DR,j),ps079(:,OP_DR)) &
                     +      int4(ri_79,mu79(:,OP_1,i),nu79(:,OP_DZ,j),ps079(:,OP_DZ))
             enddo
             dofs(i) = int3(mu79(:,OP_1,i),temp79a,temp79b)
          endif
       enddo
       
       call insert_block(dr_matrix, itri, 1, 1, temp, MAT_ADD)
       call vector_insert_block(jphi_vec%vec, itri, 1, dofs, MAT_ADD)
    enddo
    
    call sum_shared(jphi_vec%vec)
    call boundary_dc(jphi_vec%vec, jphi_vec%vec, dr_matrix)
    call finalize(dr_matrix)
 
    if(myrank.eq.0 .and. iprint.ge.2) print *, 'Solving G'
    call newsolve(dr_matrix,jphi_vec%vec,ierr)
    ! jphi_vec now contains G = (F^2 - F0^2)/2


    if(myrank.eq.0 .and. iprint.ge.2) print *, 'Calculating BZ, VZ'
    ! Calculate BZ from G = (F^2 - F0^2)/2
    call create_field(f_vec)    
    call create_field(vz_vec)
    f_vec = 0.
    vz_vec = 0.

    do itri=1,numelms

       call define_element_quadrature(itri,int_pts_main,int_pts_tor)
       call define_fields(itri,0,1,0)

       call eval_ops(itri, jphi_vec, bf079, rfac)

       call basicj_vz(x_79, z_79, temp79b)
       if(itor.eq.1) then 
          temp79a = sqrt(2.*bf079(:,OP_1) + (bzero*rzero)**2)
       else
          temp79a = sqrt(2.*bf079(:,OP_1) + bzero**2)
       end if

       do i=1,dofs_per_element
          dofs(i) = int2(mu79(:,OP_1,i),temp79a)
          dofs_vz(i) = int2(mu79(:,OP_1,i),temp79b)
       end do
       call vector_insert_block(f_vec%vec, itri, 1, dofs, MAT_ADD)
       call vector_insert_block(vz_vec%vec, itri, 1, dofs_vz, MAT_ADD)
    enddo
    
!    call sum_shared(f_vec%vec)

    if(myrank.eq.0 .and. iprint.ge.2) print *, 'Solving BZ'
    call newvar_solve(f_vec%vec, mass_mat_lhs)
    bz_field(0) = f_vec

    if(myrank.eq.0 .and. iprint.ge.2) print *, 'Solving VZ'
    call newvar_solve(vz_vec%vec, mass_mat_lhs)
    vz_field(0) = vz_vec

    call destroy_field(f_vec)
    call destroy_field(jphi_vec)
    call destroy_field(vz_vec)

    if(bzero.lt.0.) call mult(bz_field(0),-1.)

    p_field(0) = p0
    pe_field(0) = p0 - pi0
    den_field(0) = den0
    
    xlim = xmag + ln
    zlim = zmag

    call lcfs(psi_field(0))


    if(irmp.ne.2) call init_perturbations
end subroutine basicj_init

elemental real function basicj_dscale(rsq)
  use basic

  implicit none

  real, intent(in) :: rsq

  basicj_dscale = (1. + (sqrt(basicj_dvac) - 1.)*(sqrt(rsq)/ln)**basicj_dexp)**2
end function basicj_dscale

end module basicj
