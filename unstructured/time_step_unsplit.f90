module time_step_unsplit
  use field
  use matrix_mod
  use model

  type(vector_type), private :: phi_vec, phip_vec
  type(vector_type), private :: vel_vec, veln_vec
  type(vector_type), private :: pres_vec
  type(vector_type), private :: den_vec
  type(vector_type), private :: pret_vec

  ! the following pointers point to the vector containing the named field.
  ! set by assign_variables()
  type(field_type), private ::   u_v
  type(field_type), private ::  vz_v
  type(field_type), private :: chi_v
  type(field_type), private :: psi_v
  type(field_type), private ::  bz_v
  type(field_type), private ::  pe_v
  type(field_type), private :: den_v
  type(field_type), private ::   p_v
  type(field_type), private ::   e_v
  type(field_type), private ::  bf_v

  integer, private :: vecsize_phi

  ! temporary vectors
  type(vector_type), private :: b1_vel, b2_vel, b1_phi, b2_phi

  logical, private :: initialized

contains

  subroutine initialize_timestep_unsplit
    use sparse
    use basic
    implicit none

    vecsize_phi  = numvar*2 + idens + ipres + imp_bf
    
    ! add electrostatic potential equation
    if((jadv.eq.0 .and. i3d.eq.1).or.(jadv.eq.1 .and. imp_hyper.ge.1)) &
         vecsize_phi = vecsize_phi + 1
    
    ! Vectors
    call create_vector(phi_vec,      vecsize_phi)
    call create_vector(phip_vec,     vecsize_phi)
    call create_vector(q4_vec,       vecsize_phi)
     
    call create_vector(b1_phi, vecsize_phi)
    call create_vector(b2_phi, vecsize_phi)

    ! Matrices
    call set_matrix_index(s1_mat, s1_mat_index)
    call set_matrix_index(d1_mat, d1_mat_index)
    call create_mat(s1_mat, vecsize_phi, vecsize_phi, icomplex, 1)
    call create_mat(d1_mat, vecsize_phi, vecsize_phi, icomplex, 0)
#ifdef CJ_MATRIX_DUMP
    print *, "create_mat time_step s1_mat", s1_mat%imatrix     
    print *, "create_mat time_step d1_mat", d1_mat%imatrix     
#endif 
    if(i3d.eq.1) then
       call set_matrix_index(o1_mat, o1_mat_index)
       call create_mat(o1_mat, vecsize_phi, 1, icomplex, 0)
#ifdef CJ_MATRIX_DUMP
       print *, "create_mat time_step o1_mat", o1_mat%imatrix     
#endif 
    endif

    initialized = .true.
  end subroutine initialize_timestep_unsplit

  subroutine finalize_timestep_unsplit
    use basic
    use sparse
    implicit none

    if(.not.initialized) return

    call destroy_mat(s1_mat)
    call destroy_mat(d1_mat)
    if(i3d.eq.1) call destroy_mat(o1_mat)

  end subroutine finalize_timestep_unsplit


  !==========================================================
  ! assign_variables
  ! ~~~~~~~~~~~~~~~
  ! Assigns variables to appropriate vectors for time advance
  !==========================================================
  subroutine assign_variables_unsplit()
    use basic

    implicit none
        
    u_i = 1
    psi_i = 2
    vz_i = 3
    bz_i = 4
    chi_i = 5
    den_i = 2*numvar+1
    if(numvar.ge.3) then
       p_i = 6
       if(ipres.eq.1) pe_i = 2*numvar+idens+1
    else
       if(ipres.eq.1) p_i = 2*numvar+idens+1
    end if
    if(imp_bf.eq.1) then
       bf_i = 2*numvar+idens+ipres+1
       e_i = 2*numvar+idens+ipres+2
    else
       bf_i = 1
       e_i = 2*numvar+idens+ipres+1
    end if

    call associate_field(u_v,    phi_vec,      u_i)
    call associate_field(psi_v,  phi_vec,    psi_i)
    
    if(numvar.ge.2) then
       call associate_field(vz_v,  phi_vec,    vz_i)
       call associate_field(bz_v,  phi_vec,    bz_i)
    endif
    
    if(numvar.ge.3) then
       call associate_field(chi_v,  phi_vec,    chi_i)
    endif
    
    if(ipres.eq.1 .and. numvar.ge.3) then
       call associate_field(pe_v,  phi_vec,    pe_i)
    end if
    if(ipres.eq.1 .or. numvar.ge.3) then 
       call associate_field(p_v,  phi_vec,    p_i)
    end if
    
    if(idens.eq.1) then
       call associate_field(den_v,  phi_vec,    den_i)
    end if

    if(imp_bf.eq.1) then
       call associate_field(bf_v, phi_vec, bf_i)
    end if
    
    if((jadv.eq.0 .and. i3d.eq.1).or.(jadv.eq.1 .and. imp_hyper.ge.1)) then
       call associate_field(e_v, phi_vec, e_i)
    end if       
  end subroutine assign_variables_unsplit

  subroutine clear_matrices_unsplit
    use basic
    implicit none

    call clear_mat(s1_mat)
    call clear_mat(d1_mat)
    if(i3d.eq.1) call clear_mat(o1_mat)
    q4_vec = 0.
  end subroutine clear_matrices_unsplit

  subroutine finalize_matrices_unsplit
    use basic
    implicit none

    ! Finalize matrices for multiplication
    call flush(s1_mat)
    call finalize(d1_mat)
    if(i3d.eq.1) call finalize(o1_mat)
    call sum_shared(q4_vec)
  end subroutine finalize_matrices_unsplit


!======================================================================
! import_time_advance_vectors
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
! construct vectors for time-advance
!======================================================================
subroutine import_time_advance_vectors_unsplit
  use basic
  use arrays

  implicit none

  u_v = u_field(1)
  psi_v = psi_field(1)
  if(numvar.ge.2) then
     vz_v = vz_field(1)
     bz_v = bz_field(1)
  end if
  if(numvar.ge.3) then 
     chi_v = chi_field(1)
  endif

  if(ipres.eq.1 .or. numvar.ge.3) then
     p_v = p_field(1)
  end if

  if(ipres.eq.1 .and. numvar.ge.3) then
     pe_v = pe_field(1)
  end if

  if((jadv.eq.0 .and. i3d.eq.1).or.(jadv.eq.1 .and. imp_hyper.ge.1)) then
     e_v = e_field(1)
  end if

  if(idens.eq.1) den_v = den_field(1)
  if(imp_bf.eq.1) bf_v = bf_field(1)
end subroutine import_time_advance_vectors_unsplit


!======================================================================
! export_time_advance_vectors
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
! extract values from time-advance vectors
!======================================================================
subroutine export_time_advance_vectors_unsplit
  use basic
  use arrays

  implicit none

  u_field(1) = u_v
  psi_field(1) = psi_v
  if(numvar.ge.2) then
     vz_field(1) = vz_v
     bz_field(1) = bz_v
  end if
  if(numvar.ge.3) then
     chi_field(1) = chi_v
  endif

  if(ipres.eq.1 .or. numvar.ge.3) then
     p_field(1) = p_v
  end if

  if(ipres.eq.1 .and. numvar.ge.3) then
     pe_field(1) = pe_v
  else
     if(ipres.eq.1 .or. numvar.ge.3) then
        pe_field(1) = p_v
        call mult(pe_field(1), pefac)
     end if
  end if

  if((jadv.eq.0 .and. i3d.eq.1).or.(jadv.eq.1 .and. imp_hyper.ge.1)) then
     e_field(1) = e_v
  end if

  if(idens.eq.1) den_field(1) = den_v
  if(imp_bf.eq.1) bf_field(1) = bf_v
end subroutine export_time_advance_vectors_unsplit


!============================================================
! STEP_UNSPLIT
! ~~~~~~~~~~~~
! advance fields using unsplit time step
!============================================================
subroutine step_unsplit(calc_matrices)

  use boundary_conditions
  use basic
  use arrays
  use sparse
  use diagnostics
  use matrix_mod
  use model
  use auxiliary_fields

  implicit none

#ifdef CJ_MATRIX_DUMP
  integer :: counter
#endif 


  integer, intent(in) :: calc_matrices
  integer :: jer
  
  real :: tstart, tend

  if(myrank.eq.0 .and. iprint.ge.1) print *, "Solving matrix equation..."

  b1_phi = q4_vec

  if(itime_independent.eq.0) then
     ! vtemp = d1matrix_sm * phi(n)
     call matvecmult(d1_mat,phi_vec,b2_phi)
     call add(b1_phi, b2_phi)

     ! Include linear f terms
     if(numvar.ge.2 .and. i3d.eq.1 .and. imp_bf.eq.0) then
        ! b2vector = r15 * bf(n)
        call matvecmult(o1_mat,bf_field(1)%vec,b2_phi)
        call add(b1_phi, b2_phi)
     endif
  end if
   
  ! Insert boundary conditions
  if(calc_matrices.eq.1) then
     call boundary_mag(b1_phi, psi_v, bz_v, bf_v, e_v, s1_mat)
     call boundary_vel(b1_phi, u_v, vz_v, chi_v, s1_mat)
     if(idens.eq.1) call boundary_den(b1_phi, den_v, s1_mat)
     if(ipres.eq.1 .and. numvar.ge.3) call boundary_pe(b1_phi, pe_v, s1_mat)
     if(ipres.eq.1 .or. numvar.ge.3) call boundary_p(b1_phi, p_v, s1_mat)
     call finalize(s1_mat)
  else
     call boundary_mag(b1_phi, psi_v, bz_v, bf_v, e_v)
     call boundary_vel(b1_phi, u_v, vz_v, chi_v)
     if(idens.eq.1) call boundary_den(b1_phi, den_v)
     if(ipres.eq.1 .and. numvar.ge.3) call boundary_pe(b1_phi, pe_v)
     if(ipres.eq.1 .or. numvar.ge.3)  call boundary_p(b1_phi, p_v)
  endif

  ! solve linear system...LU decomposition done first time
  if(myrank.eq.0 .and. iprint.ge.1) print *, "solving.."
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
#ifdef CJ_MATRIX_DUMP
  call get_counter( s1_mat%imatrix, counter) 
  if(counter.le.0) then 
     call write_matrix(s1_mat,'s1_mat')
     call write_vector(b1_phi, 's1_mat_rhs.out')
  endif
#endif 
  call newsolve(s1_mat, b1_phi, jer)
  if(linear.eq.0 .and. iteratephi.eq.0) call clear_mat(s1_mat)
#ifdef CJ_MATRIX_DUMP
  if(counter.le.0) then 
     call write_vector(b1_phi, 's1_mat_sol.out')
  endif
#endif 
  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     t_solve_b = t_solve_b + tend - tstart
  endif
  if(jer.ne.0) then
     write(*,*) 'Error in field solve', jer
     call safestop(29)
  endif
  if(myrank.eq.0 .and. iprint.ge.1) print *, "done solve."
  
  ! new field solution at time n+1     
  phi_vec = b1_phi

  ! apply smoothing operators
  call smooth_velocity(u_v, chi_v)
  call smooth_fields(psi_v)

  ! populate temperature fields
!!$  if(ipres.eq.1 .or. numvar.ge.3) then
!!$     call get_temperatures(den_v, p_v, p_v, te_field(1), ti_field(1))
!!$  end if

  if(myrank.eq.0 .and. iprint.ge.1) print *, "Done solving matrix equation."
end subroutine step_unsplit


end module time_step_unsplit
