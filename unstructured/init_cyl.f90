!===============================
! Cylinder Equilibrium
! ~~~~~~~~~~~~~~~~~~
!===============================
module cylinder

  use spline

  implicit none

  type(spline1d), private :: p_spline   ! pressure as a function of r
  type(spline1d), private :: j_spline   ! J_phi as a function of r
  type(spline1d), private :: f_spline   ! B_phi as a function of r
  type(spline1d), private :: v_spline   ! V_phi as a function of r
  type(spline1d), private :: den_spline ! Main ion density as a function of r
  type(spline1d), private :: te_spline  ! Te as a function of r

contains

  subroutine cyl_init()
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
    use read_ascii
    use math

    implicit none

    type(field_type) :: jphi_vec, f_vec, vz_vec, p_vec, den_vec
    integer :: itri, numelms, ibound, ierr
    integer, dimension(dofs_per_element) :: imask
    vectype, dimension(dofs_per_element) :: dofs_j, dofs_bz, dofs_vz, dofs_p, &
         dofs_den
    vectype, dimension(dofs_per_element,dofs_per_element) :: temp
    type(matrix_type) :: lp_matrix

    real, allocatable :: xvals(:), yvals(:)
    integer :: nvals
    integer :: i, iout
    real, dimension(MAX_PTS) :: r
    real :: sp, sq

    if(myrank.eq.0 .and. iprint.ge.1) print *, 'Defining BASICJ Equilibrium'

    ! If iread_j == 1, then read the current density from the file
    if(iread_j.eq.1) then
       
       ! Read in A/m^2
       nvals = 0
       call read_ascii_column('profile_j', xvals, nvals, icol=1)
       call read_ascii_column('profile_j', yvals, nvals, icol=2)
       if(nvals.eq.0) then 
          if(myrank.eq.0) print *, 'Error: could not find profile_j'
          call safestop(5)
       end if
       yvals = yvals / j0_norm * 3.e5

       call create_spline(j_spline, nvals, xvals, yvals)
       deallocate(xvals, yvals)
    end if


    ! If iread_p == 1, then read the pressure from the file
    if(iread_p.eq.1) then
       
       ! Read in A/m^2
       nvals = 0
       call read_ascii_column('profile_p', xvals, nvals, icol=1)
       call read_ascii_column('profile_p', yvals, nvals, icol=2)
       if(nvals.eq.0) then 
          if(myrank.eq.0) print *, 'Error: could not find profile_p'
          call safestop(5)
       end if
       yvals = yvals / p0_norm * 10.
       yvals = yvals + pedge

       call create_spline(p_spline, nvals, xvals, yvals)
       deallocate(xvals, yvals)
    end if     


    ! If iread_f == 1, then read the toroidal field from file
    if(iread_f.eq.1) then
       
       ! Read in A/m^2
       nvals = 0
       call read_ascii_column('profile_f', xvals, nvals, icol=1)
       call read_ascii_column('profile_f', yvals, nvals, icol=2)
       if(nvals.eq.0) then 
          if(myrank.eq.0) print *, 'Error: could not find profile_f'
          call safestop(5)
       end if
       yvals = yvals / b0_norm * 1e4

       bzero = yvals(nvals)

       call create_spline(f_spline, nvals, xvals, yvals)
       deallocate(xvals, yvals)
    end if     

    ! If iread_omega == 1, then read the rotation from file
    if(iread_omega.eq.1) then
       
       ! Read in krad/s
       nvals = 0
       call read_ascii_column('profile_omega', xvals, nvals, icol=1)
       call read_ascii_column('profile_omega', yvals, nvals, icol=2)
       if(nvals.eq.0) then 
          if(myrank.eq.0) print *, 'Error: could not find profile_omega'
          call safestop(5)
       end if
       yvals = yvals / t0_norm * 1e3

       ! convert from angular frequency to velocity
       yvals = yvals * toroidal_period / (2.*pi)

       call create_spline(v_spline, nvals, xvals, yvals)
       deallocate(xvals, yvals)
    end if     

    ! If iread_ne == 1, then read the electron density from file
    if(iread_ne.eq.1) then
       
       ! Read in /m^3
       nvals = 0
       call read_ascii_column('profile_ne', xvals, nvals, icol=1)
       call read_ascii_column('profile_ne', yvals, nvals, icol=2)
       if(nvals.eq.0) then 
          if(myrank.eq.0) print *, 'Error: could not find profile_ne'
          call safestop(5)
       end if
       yvals = yvals / n0_norm / 1.e6
       yvals = yvals / z_ion      ! convert from ne to ni

       call create_spline(den_spline, nvals, xvals, yvals)
       deallocate(xvals, yvals)
    end if     

    ! If iread_te == 1, then read the electron density from file
    if(iread_te.eq.1) then
       
       ! Read in /m^3
       nvals = 0
       call read_ascii_column('profile_te', xvals, nvals, icol=1)
       call read_ascii_column('profile_te', yvals, nvals, icol=2)
       if(nvals.eq.0) then 
          if(myrank.eq.0) print *, 'Error: could not find profile_te'
          call safestop(5)
       end if
       yvals = yvals * 1.6022e-9 / (b0_norm**2/(4.*pi*n0_norm))

       call create_spline(te_spline, nvals, xvals, yvals)
       deallocate(xvals, yvals)
    end if     


    call create_field(jphi_vec)
    call create_field(p_vec)
    call create_field(f_vec)
    call create_field(vz_vec)
    call create_field(den_vec)
    jphi_vec = 0.
    p_vec = 0.
    f_vec = 0.
    vz_vec = 0.
    den_vec = 0.

    ! Create a matrix for solving -del*(A)/r = B
    call set_matrix_index(lp_matrix, lp_mat_index)
    call create_mat(lp_matrix, 1, 1, icomplex, 1)

    ibound = BOUNDARY_DIRICHLET

    ! Define Jphi
    numelms = local_elements()
    do itri=1,numelms

       call define_element_quadrature(itri,int_pts_main,int_pts_tor)
       call define_fields(itri,0,1,0)

       call get_boundary_mask(itri, ibound, imask, domain_boundary)
       
       r = sqrt((x_79 - xmag)**2 + (z_79 - zmag)**2)

       ! Define lapacian matrix for calculating psi
       temp = -intxx3(mu79(:,:,OP_1),nu79(:,:,OP_GS),ri_79)

       
       ! Define right-hand sides
       if(iread_j.eq.0) then
          temp79a = 0.
       else 
          do i=1, MAX_PTS
             call evaluate_spline(j_spline, r(i), sp, iout=iout)
             if(iout.eq.1) sp = 0.
             temp79a(i) = sp
          end do
       end if
       dofs_j = intx2(mu79(:,:,OP_1),temp79a)

       if(iread_p.eq.0) then 
          temp79a = p0
       else 
          do i=1, MAX_PTS
             call evaluate_spline(p_spline, r(i), sp)
             temp79a(i) = sp
          end do
       end if
       dofs_p = intx2(mu79(:,:,OP_1),temp79a)

       if(iread_f.eq.0) then 
          temp79a = bzero
       else 
          do i=1, MAX_PTS
             call evaluate_spline(f_spline, r(i), sp)
             temp79a(i) = sp
          end do
       end if
       dofs_bz = intx2(mu79(:,:,OP_1),temp79a)

       if(iread_omega.eq.0) then 
          temp79a = vzero
       else 
          do i=1, MAX_PTS
             call evaluate_spline(v_spline, r(i), sp)
             temp79a(i) = sp
          end do
       end if
       dofs_vz = intx2(mu79(:,:,OP_1),temp79a)

       if(iread_ne.eq.1) then
          do i=1, MAX_PTS
             call evaluate_spline(den_spline, r(i), sp)
             temp79a(i) = sp
          end do
       else if(iread_te.eq.1) then
          do i=1, MAX_PTS
             call evaluate_spline(te_spline, r(i), sp)
             if(iread_p.eq.1) then
                call evaluate_spline(p_spline, r(i), sq)
             else
                sq = p0
             end if
             temp79a(i) = sp/sq
          end do
       else
          temp79a = den0
       end if
       dofs_den = intx2(mu79(:,:,OP_1),temp79a)


       call apply_boundary_mask(itri, ibound, temp, imask)
       call apply_boundary_mask_vec(itri, ibound, dofs_j, imask)
       
       call insert_block(lp_matrix, itri, 1, 1, temp, MAT_ADD)

       call vector_insert_block(jphi_vec%vec, itri, 1, dofs_j,  VEC_ADD)
       call vector_insert_block(p_vec%vec,    itri, 1, dofs_p,  VEC_ADD)
       call vector_insert_block(f_vec%vec,    itri, 1, dofs_bz, VEC_ADD)
       call vector_insert_block(vz_vec%vec,   itri, 1, dofs_vz, VEC_ADD)
       call vector_insert_block(den_vec%vec,  itri, 1, dofs_den,VEC_ADD)
    enddo
    
    call sum_shared(jphi_vec%vec)
    call boundary_dc(jphi_vec%vec, jphi_vec%vec, lp_matrix)
    call finalize(lp_matrix)

    if(myrank.eq.0 .and. iprint.ge.2) print *, 'Solving PSI'
    call newsolve(lp_matrix,jphi_vec%vec,ierr)
    psi_field(0) = jphi_vec

    call destroy_mat(lp_matrix)

    ! Evaluate pressure
    call newvar_solve(p_vec%vec, mass_mat_lhs)
    p_field(0) = p_vec
    call destroy_field(p_vec)

    ! Evaluate toroidal field
    call newvar_solve(f_vec%vec, mass_mat_lhs)
    bz_field(0) = f_vec
    call destroy_field(f_vec)

    ! Evaluate toroidal rotation
    call newvar_solve(vz_vec%vec, mass_mat_lhs)
    vz_field(0) = vz_vec
    call destroy_field(vz_vec)

    ! Evaluate density
    call newvar_solve(den_vec%vec, mass_mat_lhs)
    den_field(0) = den_vec
    call destroy_field(den_vec)


    pe_field(0) = p_field(0)
    call mult(pe_field(0),pefac)

    
    if(xlim.eq.0) then
       xlim = xmag + ln
       zlim = zmag
    end if

    call destroy_spline(j_spline)
    call destroy_spline(p_spline)
    call destroy_spline(f_spline)
    call destroy_spline(v_spline)
    call destroy_spline(den_spline)
    call destroy_spline(te_spline)

    call lcfs(psi_field(0))

    if(irmp.ne.2) call init_perturbations
end subroutine cyl_init

end module cylinder
