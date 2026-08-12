module hypervisc
  use matrix_mod

  implicit none

  type(matrix_type) :: hyperv_lhs
  type(matrix_type) :: hyperv_rhs

  logical :: hyperv_finalized
  
contains

  ! Create the LHS and RHS hyperviscosity matrices, associated with the
  ! equations
  !
  ! dU/dt = hypv*del^4(U)
  ! dV/dt = hypv*del^4(V)
  ! dchi/dt = hypv*del^4(chi)
  !
  ! Note that while these are independent equations, they may be coupled
  ! through mixed boundary conditions, and are therefore solved together
  subroutine init_hyperv_mat()
    use basic
    use sparse
    use m3dc1_nint
    use boundary_conditions
    use vector_mod
    use matrix_mod
    use mesh_mod
    use model
    
    implicit none

    real :: thimp_sm
    integer :: itri, numelms
    integer, dimension(dofs_per_element) :: imask_vor, imask_vz, imask_chi

    vectype, dimension(dofs_per_element,dofs_per_element,numvar,numvar) :: lhs, rhs
    vectype, dimension(dofs_per_element,dofs_per_element) :: temp 

    if(hyperc.eq.0. .and. hyperv.eq.0.) return

    if(myrank.eq.0. .and. iprint.ge.1) &
         print *, 'Initializing hypervsicosity matrices'
    
    ! Create the matrices
    call set_matrix_index(hyperv_lhs, hypv_lhs_index)
    call set_matrix_index(hyperv_rhs, hypv_rhs_index)
    call create_mat(hyperv_lhs, numvar, numvar, icomplex, 1)
    call create_mat(hyperv_rhs, numvar, numvar, icomplex, 0)
    call clear_mat(hyperv_lhs)
    call clear_mat(hyperv_rhs)

    thimp_sm = thimp

    ! Populate the matrices
    numelms = local_elements()
    do itri=1,numelms

       call get_vor_mask(itri, imask_vor)
       call get_vz_mask (itri, imask_vz)
       call get_chi_mask(itri, imask_chi)

       call define_element_quadrature(itri,int_pts_main,5)
       call define_fields(itri,0,1,0)

       temp = intxx2(mu79(:,:,OP_1),nu79(:,:,OP_1))
       lhs(:,:,1,1) = temp
       rhs(:,:,1,1) = temp
       if(numvar.ge.2) then
          lhs(:,:,2,2) = temp
          rhs(:,:,2,2) = temp
       end if
       if(numvar.ge.3) then
          lhs(:,:,3,3) = temp
          rhs(:,:,3,3) = temp
       end if

       temp = intxx2(mu79(:,:,OP_LP),nu79(:,:,OP_LP))
#ifdef USE3D
       temp = temp + &
            intxx3(mu79(:,:,OP_DPP),nu79(:,:,OP_DPP),ri4_79)
#endif
       lhs(:,:,1,1) = lhs(:,:,1,1) +  thimp_sm    *hypv*dt*temp
       rhs(:,:,1,1) = rhs(:,:,1,1) + (thimp_sm-1.)*hypv*dt*temp
       if(numvar.ge.2) then 
          lhs(:,:,2,2) = lhs(:,:,2,2) +  thimp_sm    *hypv*dt*temp
          rhs(:,:,2,2) = rhs(:,:,2,2) + (thimp_sm-1.)*hypv*dt*temp
       end if
       if(numvar.ge.3) then
          lhs(:,:,3,3) = lhs(:,:,3,3) +  thimp_sm    *hypc*dt*temp
          rhs(:,:,3,3) = rhs(:,:,3,3) + (thimp_sm-1.)*hypc*dt*temp
       end if

       call apply_boundary_mask(itri,0,lhs(:,:,1,1),imask_vor,tags=BOUND_ANY)
       call apply_boundary_mask(itri,0,rhs(:,:,1,1),imask_vor,tags=BOUND_ANY)
       call insert_block(hyperv_lhs, itri, 1, 1, lhs(:,:,1,1), MAT_ADD)
       call insert_block(hyperv_rhs, itri, 1, 1, rhs(:,:,1,1), MAT_ADD)
       
       if(numvar.ge.2) then
          call apply_boundary_mask(itri,0,lhs(:,:,2,2),imask_vz,tags=BOUND_ANY)
          call apply_boundary_mask(itri,0,rhs(:,:,2,2),imask_vz,tags=BOUND_ANY)
          call insert_block(hyperv_lhs, itri, 2, 2, lhs(:,:,2,2), MAT_ADD)
          call insert_block(hyperv_rhs, itri, 2, 2, rhs(:,:,2,2), MAT_ADD)
       end if
       if(numvar.ge.3) then
          call apply_boundary_mask(itri,0,lhs(:,:,3,3),imask_chi,tags=BOUND_ANY)
          call apply_boundary_mask(itri,0,rhs(:,:,3,3),imask_chi,tags=BOUND_ANY)
          call insert_block(hyperv_lhs, itri, 3, 3, lhs(:,:,3,3), MAT_ADD)
          call insert_block(hyperv_rhs, itri, 3, 3, rhs(:,:,3,3), MAT_ADD)
       end if
    enddo

    ! insert boundary conditions
    call flush(hyperv_lhs)
    call flush(hyperv_rhs)
    call finalize(hyperv_rhs)

    hyperv_finalized = .false.
   
  end subroutine init_hyperv_mat


  
  ! Apply the hyperviscosity step
  subroutine apply_hyperv(u_f, vz_f, chi_f)
    use basic
    use vector_mod
    use field
    use model
    
    implicit none

    type(field_type), intent(inout) :: u_f, vz_f, chi_f

    integer :: ierr
    type(vector_type) :: vel, rhs_vec
    type(field_type) :: u_tmp, vz_tmp, chi_tmp

    if(hypv.eq.0. .and. hypc.eq.0.) return
    
    call create_vector(vel, numvar)
    call create_vector(rhs_vec, numvar)

    call associate_field(u_tmp,vel,1)
    if(numvar.ge.2) call associate_field(vz_tmp,  vel,2)
    if(numvar.ge.3) call associate_field(chi_tmp, vel,3)

    ! Copy data into temporary velocity vector
    u_tmp = u_f
    if(numvar.ge.2) vz_tmp = vz_f
    if(numvar.ge.3) chi_tmp = chi_f

    ! Calculate RHS vector
    call matvecmult(hyperv_rhs,vel,rhs_vec)

    if(hyperv_finalized) then
       call boundary_vel(rhs_vec, u_tmp, vz_tmp, chi_tmp)
    else
       call boundary_vel(rhs_vec, u_tmp, vz_tmp, chi_tmp, hyperv_lhs)
       call finalize(hyperv_lhs)
       hyperv_finalized = .true.
    end if

    call newsolve(hyperv_lhs, rhs_vec, ierr)
    vel = rhs_vec
    
    u_f = u_tmp
    if(numvar.ge.2) vz_f = vz_tmp
    if(numvar.ge.3) chi_f = chi_tmp
    
    call destroy_vector(vel)
    call destroy_vector(rhs_vec)
    
  end subroutine apply_hyperv
  
end module hypervisc
