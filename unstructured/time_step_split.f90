module time_step_split
  use field
  use matrix_mod
  use model

  type(vector_type), private :: phi_vec, phip_vec
  type(vector_type), private :: vel_vec, veln_vec
  type(vector_type), private :: pres_vec
  type(vector_type), private :: den_vec, denold_vec
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
  type(field_type), private ::  te_v
  type(field_type), private ::  ti_v

  integer :: vecsize_vel, vecsize_phi, vecsize_n, vecsize_p, vecsize_t

  ! temporary vectors
  type(vector_type), private :: b1_vel, b2_vel, b1_phi, b2_phi

  logical, private :: initialized

contains

  subroutine initialize_timestep_split
    use sparse
    use basic
    implicit none

    vecsize_phi = numvar
    vecsize_vel = numvar
    vecsize_n = 1
    vecsize_p = 1
!   if(ipressplit.eq.1 .and. numvar.eq.3 .and. linear.eq.0 .and. eqsubtract.eq.0) then
    if(ipressplit.eq.1 .and. numvar.eq.3) then
       !split pressure solve from field solve
       if(ipres.eq.0 .and. itemp.eq.0) then  
          imode = 1
          vecsize_phi = 2
          vecsize_p   = 1
       endif
       !Solve for Temperature instead of Pressure
       if(ipres.eq.0 .and. itemp.eq.1) then  
          imode = 2
          vecsize_phi = 2
          vecsize_p   = 1
       endif
       !electron and total pressures solved together
       if(ipres.eq.1 .and. itemp.eq.0) then  
          imode = 3
          vecsize_phi = 2
          vecsize_p   = 2
       endif
       !electron and ion Temperatures solved together
       if(ipres.eq.1 .and. itemp.eq.1) then  
          imode = 4
          vecsize_phi = 2
          vecsize_p   = 2
       endif
    endif
    vecsize_t = vecsize_p

    ! add electrostatic potential equation
    if(jadv.eq.0 .and. i3d.eq.1) vecsize_phi = vecsize_phi + 1
  
    ! add bf equation
    if(imp_bf.eq.1) vecsize_phi = vecsize_phi + 1

    ! Vectors
    call create_vector(phi_vec,      vecsize_phi)
    call create_vector(phip_vec,     vecsize_phi)
    call create_vector(q4_vec,       vecsize_phi)
     
    call create_vector(b1_phi, vecsize_phi)
    call create_vector(b2_phi, vecsize_phi)

    call create_vector(vel_vec,      vecsize_vel)
    call create_vector(veln_vec,     vecsize_vel)
    call create_vector(r4_vec,       vecsize_vel)

    if(ipres.eq.1 .or. ipressplit.eq.1) then
       call create_vector(pres_vec,    vecsize_p)
       call create_vector(qp4_vec,     vecsize_p)
    endif

    call create_vector(den_vec,    vecsize_n)
    call create_vector(denold_vec, vecsize_n)
    call create_vector(qn4_vec,    vecsize_n)

    call create_vector(b1_vel, vecsize_vel)
    call create_vector(b2_vel, vecsize_vel)

    call create_vector(pret_vec,    vecsize_t)

    ! Matrices
    call set_matrix_index(s1_mat, s1_mat_index)
    call set_matrix_index(d1_mat, d1_mat_index)
    call set_matrix_index(q1_mat, q1_mat_index)
    call set_matrix_index(r14_mat, r14_mat_index)
       
    call create_mat(s1_mat, vecsize_vel, vecsize_vel, icomplex, .true.)
    call create_mat(d1_mat, vecsize_vel, vecsize_vel, icomplex, .false.)
    call create_mat(q1_mat, vecsize_vel, vecsize_phi, icomplex, .false.)
    call create_mat(r14_mat, vecsize_vel, vecsize_n, icomplex, .false.)

#ifdef CJ_MATRIX_DUMP
    print *, "create_mat time_step s1_mat", s1_mat%imatrix     
    print *, "create_mat time_step d1_mat", d1_mat%imatrix     
    print *, "create_mat time_step q1_mat", q1_mat%imatrix     
    print *, "create_mat time_step r14_mat", r14_mat%imatrix     
#endif 
    if(i3d.eq.1) then
       call set_matrix_index(o1_mat, o1_mat_index)
       call create_mat(o1_mat, vecsize_vel, 1, icomplex, .false.)
#ifdef CJ_MATRIX_DUMP
       print *, "create_mat time_step o1_mat", o1_mat%imatrix     
#endif 
    endif
    if((ipres.eq.1 .and. numvar.lt.3) .or. ipressplit.gt.0) then
       call set_matrix_index(p1_mat, p1_mat_index)
       call create_mat(p1_mat, vecsize_vel, vecsize_p, icomplex, .false.)
    end if

    call set_matrix_index(s2_mat, s2_mat_index)
    call set_matrix_index(d2_mat, d2_mat_index)
    call set_matrix_index(r2_mat, r2_mat_index)
    call set_matrix_index(q2_mat, q2_mat_index)
    call create_mat(s2_mat, vecsize_phi, vecsize_phi, icomplex, .true.)
    call create_mat(d2_mat, vecsize_phi, vecsize_phi, icomplex, .false.)
    call create_mat(r2_mat, vecsize_phi, vecsize_vel, icomplex, .false.)
    call create_mat(q2_mat, vecsize_phi, vecsize_vel, icomplex, .false.)

#ifdef CJ_MATRIX_DUMP
    print *, "create_mat time_step s2_mat", s2_mat%imatrix     
    print *, "create_mat time_step d2_mat", d2_mat%imatrix     
    print *, "create_mat time_step r2_mat", r2_mat%imatrix     
    print *, "create_mat time_step q2_mat", q2_mat%imatrix     
#endif 

    if(idens.eq.1) then
       call set_matrix_index(r42_mat, r42_mat_index)
       call set_matrix_index(q42_mat, q42_mat_index)
       call create_mat(r42_mat, vecsize_phi, 1, icomplex, .false.)
       call create_mat(q42_mat, vecsize_phi, 1, icomplex, .false.)
#ifdef CJ_MATRIX_DUMP
       print *, "create_mat time_step r42_mat", r42_mat%imatrix
       print *, "create_mat time_step q42_mat", q42_mat%imatrix
#endif 
    endif
    
    if(i3d.eq.1) then
       call set_matrix_index(o2_mat, o2_mat_index)
       call create_mat(o2_mat, vecsize_phi, 1, icomplex, .false.)
       if(ipressplit.eq.1 .or. ipres.eq.1) then
          call set_matrix_index(o3_mat, o3_mat_index)
          call create_mat(o3_mat, vecsize_p, 1, icomplex, .false.)
       endif
         
#ifdef CJ_MATRIX_DUMP
       print *, "create_mat time_step o2_mat", o2_mat%imatrix     
#endif 
    endif
       
    if(idens.eq.1) then
       call set_matrix_index(s8_mat, s8_mat_index)
       call set_matrix_index(d8_mat, d8_mat_index)
       call set_matrix_index(r8_mat, r8_mat_index)
       call set_matrix_index(q8_mat, q8_mat_index)
       call create_mat(s8_mat, vecsize_n, vecsize_n, icomplex, .true.)
       call create_mat(d8_mat, vecsize_n, vecsize_n, icomplex, .false.)
       call create_mat(r8_mat, vecsize_n, vecsize_vel, icomplex, .false.)
       call create_mat(q8_mat, vecsize_n, vecsize_vel, icomplex, .false.)
#ifdef CJ_MATRIX_DUMP
       print *, "create_mat time_step s8_mat", s8_mat%imatrix     
       print *, "create_mat time_step d8_mat", d8_mat%imatrix     
       print *, "create_mat time_step r8_mat", r8_mat%imatrix     
       print *, "create_mat time_step q8_mat", q8_mat%imatrix     
#endif 
    endif

    if(ipres.eq.1 .or. ipressplit.eq.1) then
       call set_matrix_index(s9_mat, s9_mat_index)
       call set_matrix_index(d9_mat, d9_mat_index)
       call set_matrix_index(r9_mat, r9_mat_index)
       call set_matrix_index(q9_mat, q9_mat_index)
       call set_matrix_index(o9_mat, o9_mat_index)
       call create_mat(s9_mat, vecsize_p, vecsize_p, icomplex, .true.)
       call create_mat(d9_mat, vecsize_p, vecsize_p, icomplex, .false.)
       call create_mat(r9_mat, vecsize_p, vecsize_vel, icomplex, .false.)
       call create_mat(q9_mat, vecsize_p, vecsize_vel, icomplex, .false.)
       call create_mat(o9_mat, vecsize_p, vecsize_phi, icomplex, .false.)
#ifdef CJ_MATRIX_DUMP
       print *, "create_mat time_step s9_mat", s9_mat%imatrix
       print *, "create_mat time_step d9_mat", d9_mat%imatrix
       print *, "create_mat time_step r9_mat", r9_mat%imatrix
       print *, "create_mat time_step q9_mat", q9_mat%imatrix
       print *, "create_mat time_step o9_mat", o9_mat%imatrix
#endif 
    endif
       
!      if(ipressplit.eq.1) then
!         call set_matrix_index(s11_mat, s11_mat_index)
!         call set_matrix_index(d11_mat, d11_mat_index)
!         call create_mat(s11_mat, vecsize_n, vecsize_n, icomplex, .true.)
!         call create_mat(d11_mat, vecsize_n, vecsize_n, icomplex, .false.)
#ifdef CJ_MATRIX_DUMP
!         print *, "create_mat time_step s11_mat", s11_mat%imatrix
!         print *, "create_mat time_step d11_mat", d11_mat%imatrix
#endif 
!         call set_matrix_index(s12_mat, s12_mat_index)
!         call set_matrix_index(d12_mat, d12_mat_index)
!         call create_mat(s12_mat, vecsize_n, vecsize_n, icomplex, .true.)
!         call create_mat(d12_mat, vecsize_n, vecsize_n, icomplex, .false.)
#ifdef CJ_MATRIX_DUMP
!         print *, "create_mat time_step s12_mat", s12_mat%imatrix
!         print *, "create_mat time_step d12_mat", d12_mat%imatrix
#endif 
!      endif

       initialized = .true.
  end subroutine initialize_timestep_split

  subroutine finalize_timestep_split
    use basic
    implicit none

    if(.not.initialized) return
      
    call destroy_mat(s1_mat)
    call destroy_mat(d1_mat)
    call destroy_mat(q1_mat)
    call destroy_mat(r14_mat)
    if(i3d.eq.1) call destroy_mat(o1_mat)

    call destroy_mat(s2_mat)
    call destroy_mat(d2_mat)
    call destroy_mat(r2_mat)
    call destroy_mat(q2_mat)
    if(idens.eq.1) then
       call destroy_mat(r42_mat)
       call destroy_mat(q42_mat)
    end if
    if(i3d.eq.1) then
      call destroy_mat(o2_mat)
      if(ipres.eq.1 .or. ipressplit.gt.0) call destroy_mat(o3_mat)
    endif
    if((ipres.eq.1 .and. numvar.lt.3) .or. ipressplit.gt.0) call destroy_mat(p1_mat)
       
    if(idens.eq.1) then
       call destroy_mat(s8_mat)
       call destroy_mat(d8_mat)
       call destroy_mat(r8_mat)
       call destroy_mat(q8_mat)
    endif

    if(ipres.eq.1 .or. ipressplit.eq.1) then
       call destroy_mat(s9_mat)
       call destroy_mat(d9_mat)
       call destroy_mat(r9_mat)
       call destroy_mat(q9_mat)
       call destroy_mat(o9_mat)
    endif
       
!      if(ipressplit.eq.1) then
!         call destroy_mat(s11_mat)
!         call destroy_mat(d11_mat)
!         call destroy_mat(s12_mat)
!         call destroy_mat(d12_mat)
!      endif
  end subroutine finalize_timestep_split


  subroutine assign_variables_split()
    use basic

    implicit none
        
    u_i = 1
    psi_i = 1
    vz_i = 2
    bz_i = 2
    chi_i = 3
    den_i = 1
    if(imp_bf.eq.1) then
       if(jadv.eq.0 .and. i3d.eq.1) then
         bf_i = vecsize_phi - 1
         e_i =  vecsize_phi
       else
         bf_i = vecsize_phi
         e_i = 1
       endif
    else
       if(jadv.eq.0 .and. i3d.eq.1) then
         bf_i = 1
         e_i = vecsize_phi
       else
         bf_i = 1
         e_i = 1
       endif
    end if
    if(ipressplit.eq.1) then
       p_i  = 1
       te_i = 1
       if(ipres.eq.1) then
          pe_i = 2
          ti_i = 2
       endif
    else  ! ipressplit.eq.0
       pe_i = 3
       if(ipres.eq.1) then
          p_i = 1
       end if
       te_i = 1
    endif

    call associate_field(u_v,    vel_vec,      u_i)
    call associate_field(psi_v,  phi_vec,    psi_i)
    
    if(numvar.ge.2) then
       call associate_field(vz_v,  vel_vec,    vz_i)
       call associate_field(bz_v,  phi_vec,    bz_i)
    endif
    
    if(numvar.ge.3) then
       call associate_field(chi_v,  vel_vec,    chi_i)
    endif
    
    if(numvar.ge.3 .or. ipres.eq.1) then
       call associate_field(te_v, pret_vec, te_i)
    endif
    
    if(ipressplit.eq.1) then
       call associate_field(p_v, pres_vec, p_i)
       if(ipres.eq.1) then
          call associate_field(pe_v, pres_vec, pe_i)
          call associate_field(ti_v, pret_vec, ti_i)
       endif
    else   !  on ipressplit.eq.1
       if(numvar.ge.3) then
          call associate_field(pe_v,   phi_vec,     pe_i)
       endif
       if(ipres.eq.1) then
          call associate_field(p_v,  pres_vec,    p_i)
       end if
    endif  !  on ipressplit.eq.1
    
    if(idens.eq.1) then
       call associate_field(den_v,  den_vec,    den_i)
    end if
    
    if(jadv.eq.0 .and. i3d.eq.1) then
       call associate_field(e_v, phi_vec, e_i)
    end if
    
    if(imp_bf.eq.1) then
       call associate_field(bf_v, phi_vec, bf_i)
    end if
         
    u_off = (u_i-1)*dofs_per_node
    psi_off = (psi_i-1)*dofs_per_node
    vz_off = (vz_i-1)*dofs_per_node
    bz_off = (bz_i-1)*dofs_per_node
    chi_off = (chi_i-1)*dofs_per_node
    pe_off = (pe_i-1)*dofs_per_node
    den_off = (den_i-1)*dofs_per_node
    p_off = (p_i-1)*dofs_per_node
    bf_off = (bf_i-1)*dofs_per_node
    e_off = (e_i-1)*dofs_per_node

  end subroutine assign_variables_split

  subroutine clear_matrices_split
    use basic
    implicit none

    call clear_mat(s1_mat)
    call clear_mat(d1_mat)
    call clear_mat(q1_mat)
    call clear_mat(r14_mat)
    if(i3d.eq.1) call clear_mat(o1_mat)
    if((ipres.eq.1 .and. numvar.lt.3) .or. ipressplit.eq.1) &
         call clear_mat(p1_mat)
    r4_vec = 0.

    
    call clear_mat(s2_mat)
    call clear_mat(d2_mat)
    call clear_mat(r2_mat)
    call clear_mat(q2_mat)
    if(idens.eq.1) then
       call clear_mat(r42_mat)
       call clear_mat(q42_mat)
    end if
    if(i3d.eq.1) then
       call clear_mat(o2_mat)
       if(ipres.eq.1 .or. ipressplit.eq.1) call clear_mat(o3_mat)
    endif
    q4_vec = 0.
    
    if(idens.eq.1) then
       call clear_mat(s8_mat)
       call clear_mat(d8_mat)
       call clear_mat(r8_mat)
       call clear_mat(q8_mat)
       qn4_vec = 0.
    end if
    
    if(ipres.eq.1 .or. ipressplit.eq.1) then
       call clear_mat(s9_mat)
       call clear_mat(d9_mat)
       call clear_mat(r9_mat)
       call clear_mat(q9_mat)
       call clear_mat(o9_mat)
       qp4_vec = 0.
    end if
  end subroutine clear_matrices_split

  subroutine finalize_matrices_split
    use basic
    implicit none

    call flush(s1_mat)
    call finalize(d1_mat)
    call finalize(q1_mat)
    call finalize(r14_mat)
    if(i3d.eq.1) call finalize(o1_mat)
    if((ipres.eq.1 .and. numvar.lt.3) .or. ipressplit.eq.1) &
         call finalize(p1_mat)
    call sum_shared(r4_vec)

    if(myrank.eq.0 .and. iprint.ge.1) &
         print *, " before field finalize..."
    
    call flush(s2_mat)
    call finalize(d2_mat)
    call finalize(r2_mat)
    call finalize(q2_mat)
    if(i3d.eq.1) then
      call finalize(o2_mat)
      if(ipres.eq.1 .or. ipressplit.eq.1) call finalize(o3_mat)
    endif
    if(idens.eq.1) then
       call finalize(r42_mat)
       call finalize(q42_mat)
    end if
    call sum_shared(q4_vec)
    
    if(idens.eq.1) then
       call flush(s8_mat)
       call finalize(d8_mat)
       call finalize(q8_mat)
       call finalize(r8_mat)
       call sum_shared(qn4_vec)
       
       if(myrank.eq.0 .and. iprint.ge.1) &
            print *, " before pressure finalize...", vecsize_vel, vecsize_p, vecsize_phi   
    end if
       
    if(ipres.eq.1 .or. ipressplit.eq.1) then
       call flush(s9_mat)
       call finalize(d9_mat)
       call finalize(q9_mat)
       call finalize(r9_mat)
       if(myrank.eq.0 .and. iprint.ge.1) print *, 'before call to finalize(o9_mat)'
       call finalize(o9_mat)
       if(myrank.eq.0 .and. iprint.ge.1) print *, 'before call to sum_shared'
       call sum_shared(qp4_vec)
    end if
end subroutine finalize_matrices_split


!======================================================================
! import_time_advance_vectors
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
! construct vectors for time-advance
!======================================================================
subroutine import_time_advance_vectors_split
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

  if(ipressplit.eq.0) then
     if(numvar.ge.3) then
        if(ipres.eq.1) then
        ! store electron pressure in pe_v
           pe_v = pe_field(1)
           p_v  = p_field(1)
        else
        ! store total pressure in pe_v
           pe_v = p_field(1)
        end if
     else 
        if(ipres.eq.1) then
           p_v = p_field(1)
        endif
     end if
  else    ! on ipressplit.eq.0
     p_v  = p_field(1)
     te_v = te_field(1)
     if(ipres.eq.1) then
        pe_v = pe_field(1)
        ti_v = ti_field(1)
     endif
  endif   ! on ipressplit.eq.0

  if(idens.eq.1) den_v = den_field(1)
  if(imp_bf.eq.1) bf_v = bf_field(1)

end subroutine import_time_advance_vectors_split


!======================================================================
! export_time_advance_vectors
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
! extract values from time-advance vectors
!======================================================================
subroutine export_time_advance_vectors_split
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
 
  if(ipressplit.eq.0) then
     if((numvar.ge.3 .or. ipres.eq.1) .and. isplitstep.ge.1) then
        te_field(1) = te_v
     endif
     if(numvar.ge.3) then

        if(ipres.eq.1) then
        ! electron pressure is stored in pe_v
           pe_field(1) = pe_v
           p_field(1)  = p_v
        else
        ! total pressure is stored in pe_v
           pe_field(1) = pe_v
           call mult(pe_field(1), pefac)
           p_field(1)  = pe_v
        end if
     else 
        if(ipres.eq.1) then
           pe_field(1) = p_v
           call mult(pe_field(1), pefac)
           p_field(1) = p_v
        endif
     endif
  else    ! on ipressplit.eq.0
      p_field(1) = p_v 
      te_field(1) = te_v
      if(ipres.eq.0) then
         pe_field(1) = p_v
         call mult(pe_field(1), pefac)
         ti_field(1) = te_v
         call mult(ti_field(1), (1.-pefac)/pefac)
      else
         pe_field(1) = pe_v
         ti_field(1) = ti_v 
      endif

  endif   ! on ipressplit.eq.0

  if(idens.eq.1) den_field(1) = den_v
  if(imp_bf.eq.1) bf_field(1) = bf_v

end subroutine export_time_advance_vectors_split


!============================================================
! STEP_SPLIT
! ~~~~~~~~~~
! advance fields using split time step
!============================================================
subroutine step_split(calc_matrices)

  use basic
  use arrays
  use sparse
  use diagnostics
  use matrix_mod
  use boundary_conditions
  use mesh_mod 
  use model
  use transport_coefficients

  implicit none

#ifdef CJ_MATRIX_DUMP
  integer :: counter
#endif


  integer, intent(in) :: calc_matrices
  real :: tstart, tend, t_bound
  real :: ans
  integer :: jer, i, numnodes
  type(vector_type) :: temp, temp2

  type(field_type) :: phip_1, phip_2, phip_3, temp_field_1, temp_field_2
  call associate_field(phip_1, phip_vec, 1)
  if(numvar.ge.2) call associate_field(phip_2, phip_vec, 2)
  if(numvar.ge.3 .and. ipressplit.eq.0) call associate_field(phip_3, phip_vec, 3)

  t_bound = 0.

  ! Store current-time velocity matrices for use in field advance
  veln_vec = vel_vec

  if(istatic.ne.1) then

     ! Advance Velocity
     ! ================
     if(myrank.eq.0 .and. iprint.ge.1) print *, " Advancing velocity"
  
     ! d1matrix_sm * vel(n)
     call matvecmult(d1_mat,vel_vec,b1_vel)
  
     ! q1matrix_sm * phi(n)
     if(ipres.eq.1 .and. numvar.ge.3 .and. ipressplit.eq.0) then
        ! replace electron pressure with total pressure

        phip_vec = phi_vec
        phip_3 = p_v
        call matvecmult(q1_mat, phip_vec, b2_vel)
     else
        call matvecmult(q1_mat, phi_vec , b2_vel)
     endif
     call add(b1_vel, b2_vel)

     ! If ipres==1 and numvar<3 or ipressplit==1, need to add in pressure contribution
     ! separately
     if((ipres.eq.1 .and. numvar.lt.3) .or. ipressplit.gt.0) then
        call matvecmult(p1_mat,pres_vec,b2_vel)
        call add(b1_vel, b2_vel)
     endif
  
     ! r14matrix_sm * den(n)
     if(idens.eq.1 .and. (gravr.ne.0 .or. gravz.ne.0)) then
        call matvecmult(r14_mat,den_v%vec,b2_vel)
        call add(b1_vel, b2_vel)
     endif

     ! o1matrix_sm * bf(n)
     if(numvar.ge.2 .and. i3d.eq.1 .and. imp_bf .eq. 0) then
        call matvecmult(o1_mat,bf_field(1)%vec,b2_vel)
        call add(b1_vel, b2_vel)
     endif

     call add(b1_vel, r4_vec)
  
     ! apply boundary conditions
     if(myrank.eq.0 .and. iprint.ge.2) print *, "  inserting bcs"
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     if(calc_matrices.eq.1) then
        call boundary_vel(b1_vel, u_v, vz_v, chi_v, s1_mat)
        call finalize(s1_mat)
     else
        call boundary_vel(b1_vel, u_v, vz_v, chi_v)
     endif

     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_bound = t_bound + tend - tstart
     end if

     ! solve linear system with rhs in vtemp (note LU-decomp done first time)
     if(myrank.eq.0 .and. iprint.ge.2) print *, "  solving"
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

#ifdef CJ_MATRIX_DUMP
     call get_counter( s1_mat%imatrix, counter) 
     if(counter.le.0) then 
        call write_matrix(s1_mat,'s1_mat')
        call write_vector(b1_vel, 's1_mat_rhs.out')
     endif
#endif

     call newsolve(s1_mat, b1_vel, jer)
     if(linear.eq.0) call clear_mat(s1_mat)

     if(idifv .gt.0) then
        call add(b1_vel,vel_vec)
     endif

#ifdef CJ_MATRIX_DUMP
     if(counter.le.0) then 
        call write_vector(b1_vel, 's1_mat_sol.out')
     endif
#endif 

     if(myrank.eq.0 .and. iprint.ge.2) print *, "  done solve"
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_solve_v = t_solve_v + tend - tstart
     endif
     if(jer.ne.0) then
        write(*,*) 'Error in velocity solve', jer
        call safestop(42)
     endif
  
     !.....new velocity solution at time n+1 (or n* for second order advance)
     vel_vec = b1_vel
   
     ! apply smoothing operators
     ! ~~~~~~~~~~~~~~~~~~~~~~~~~
     call smooth_velocity(u_v, chi_v)
  end if
  
  ! Advance Density
  ! ===============
  if(idens.eq.1) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, " Advancing density"
     
     call create_vector(temp, vecsize_n)
     call create_vector(temp2, vecsize_n)

     ! r8matrix_sm * vel(n+1)
     call matvecmult(r8_mat,vel_vec,temp)

     ! q8matrix_sm * vel(n)
     call matvecmult(q8_mat,veln_vec,temp2)
     call add(temp, temp2)
          
     ! temp = d8matrix_sm * phi(n)    
     call matvecmult(d8_mat,den_vec,temp2)
     call add(temp, temp2)
     
     call add(temp, qn4_vec)
     
     ! Insert boundary conditions
     if(myrank.eq.0 .and. iprint.ge.2) print *, "  inserting bcs"
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     if(calc_matrices.eq.1) then
        call boundary_den(temp, den_v, s8_mat)
        call finalize(s8_mat)
     else
        call boundary_den(temp, den_v)
     endif
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_bound = t_bound + tend - tstart
     end if

     if(myrank.eq.0 .and. iprint.ge.2) print *, "  solving"
     ! solve linear system...LU decomposition done first time
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

#ifdef CJ_MATRIX_DUMP
  if(ntime.eq.2) then 
     call write_matrix(s8_mat,'s8_mat')
     call write_vector(temp, 's8_mat_rhs.out')
  endif
#endif 

     call newsolve(s8_mat, temp, jer)

     if(idiff .gt. 0) then
         call add(temp,den_vec)  ! add time n density to increment to get time n+1 value
     endif

     if(linear.eq.0) call clear_mat(s8_mat)

#ifdef CJ_MATRIX_DUMP
  if(ntime.eq.2) then
     call write_vector(temp, 's8_mat_sol.out')
  endif
#endif 

     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_solve_n = t_solve_n + tend - tstart
     endif
     if(jer.ne.0) then
        write(*,*) 'Error in density solve', jer
        call safestop(29)
     endif
     if(myrank.eq.0 .and. iprint.ge.2) print *, "  done solve"
     
     ! new field solution at time n+1 (or n* for second order advance)
     denold_vec = den_vec
     den_vec = temp
     call destroy_vector(temp)
     call destroy_vector(temp2)

     if(irecalc_eta.eq.1) then
        call export_time_advance_vectors_split
        call define_transport_coefficients
     end if
  endif    ! on idens=1
     
  !
  ! Advance Pressure
  ! ================
  if((ipressplit.eq.0 .and. ipres.eq.1) .or. (ipressplit.eq.1 .and. itemp.eq.0)) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Pressure"

     call create_vector(temp, vecsize_p)
     call create_vector(temp2, vecsize_p)
     
     ! r9matrix_sm * vel(n+1)
     call matvecmult(r9_mat,vel_vec,temp)
     call mult(temp, -1.)     ! added 1/1/2013   SCJ
            
     ! q9matrix_sm * vel(n)
     call matvecmult(q9_mat,veln_vec,temp2)
     call add(temp, temp2)
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Pressure -- before o9matrix"

     ! o9matrix_sm * phi(n)
     call matvecmult(o9_mat,phi_vec,temp2)
     call add(temp, temp2)

     ! temp = d9matrix_sm * pres(n)
     call matvecmult(d9_mat,pres_vec,temp2)
     call add(temp, temp2)

     call add(temp, qp4_vec)
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Pressure -- before boundary_pres"

     
     ! Insert boundary conditions
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     if(calc_matrices.eq.1) then
        call boundary_p(temp, p_v, s9_mat)
        if(ipressplit.eq.1 .and. ipres.eq.1) &
             call boundary_pe(temp, pe_v, s9_mat)
        call finalize(s9_mat)
     else
        call boundary_p(temp, p_v)
        if(ipressplit.eq.1 .and. ipres.eq.1) &
             call boundary_pe(temp, pe_v)
     endif
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_bound = t_bound + tend - tstart
     end if

     
     ! solve linear system...LU decomposition done first time
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

#ifdef CJ_MATRIX_DUMP
  if(ntime.eq.2) then 
     call write_matrix(s9_mat,'s9_mat')
     call write_vector(temp, 's9_mat_rhs.out')
  endif
#endif 

     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Pressure--before newsolve"
     call newsolve(s9_mat, temp, jer)
     if(linear.eq.0) call clear_mat(s9_mat)

    if(idiff .gt. 0) then
         call add(temp,pres_vec)
     endif


#ifdef CJ_MATRIX_DUMP
  if(ntime.eq.2) then
     call write_vector(temp, 's9_mat_sol.out')
  endif
#endif 

     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_solve_p = t_solve_p + tend - tstart
     endif
     if(jer.ne.0) then
        write(*,*) 'Error in pressure solve', jer
        call safestop(29)
     endif
     
     ! new field solution at time n+1 (or n* for second order advance)
     pres_vec = temp
     call destroy_vector(temp)
     call destroy_vector(temp2)
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Pressure--before get_temperatures"
     call get_temperatures(den_v, p_v, pe_v, te_v, ti_v)
  endif

  ! Advance Temperature
  ! ================
  if(ipressplit.eq.1 .and. itemp.eq.1) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Temperature"

     call create_vector(temp, vecsize_p)
     call create_vector(temp2, vecsize_p)
     
     ! r9matrix_sm * vel(n+1)
     call matvecmult(r9_mat,vel_vec,temp)
     call mult(temp, -1.)    ! added 1/1/2013   SCJ

     ! q9matrix_sm * vel(n)
     call matvecmult(q9_mat,veln_vec,temp2)
     call add(temp, temp2)
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Temperature -- before o9matrix"

     ! o9matrix_sm * phi(n)
     call matvecmult(o9_mat,phi_vec,temp2)
     call add(temp, temp2)

     
     ! temp = d9matrix_sm * temp(n)
     call matvecmult(d9_mat,pret_vec,temp2)
     call add(temp, temp2)

 
     ! Include linear f terms
     if(i3d.eq.1 .and. (ipres.eq.1 .or. ipressplit.eq.1)) then
        call matvecmult(o3_mat,bf_field(1)%vec,temp2)
        call add(temp, temp2)
     endif

    
     call add(temp, qp4_vec)
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Temperature -- before boundary_temp"
     
     ! Insert boundary conditions
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     if(calc_matrices.eq.1) then
        call boundary_te(temp, te_v, s9_mat)
        if(ipres.eq.1) call boundary_ti(temp, ti_v, s9_mat)
        call finalize(s9_mat)
     else
        call boundary_te(temp, te_v)
        if(ipres.eq.1) call boundary_ti(temp, ti_v)
     endif
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_bound = t_bound + tend - tstart
     end if

     
     ! solve linear system...LU decomposition done first time
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

#ifdef CJ_MATRIX_DUMP
  if(ntime.eq.2) then 
     call write_matrix(s9_mat,'s9_mat')
     call write_vector(temp, 's9_mat_rhs.out')
  endif
#endif 

     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Temperature--before newsolve"
     call newsolve(s9_mat, temp, jer)
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Temperature--after newsolve"
     if(linear.eq.0) call clear_mat(s9_mat)

#ifdef CJ_MATRIX_DUMP
  if(ntime.eq.2) then
     call write_vector(temp, 's9_mat_sol.out')
  endif
#endif 

     if(idiff .gt. 0) then
         call add(temp,pret_vec)
     endif

     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_solve_p = t_solve_p + tend - tstart
     endif
     if(jer.ne.0) then
        write(*,*) 'Error in pressure solve', jer
        call safestop(29)
     endif
     
     ! new field solution at time n+1 (or n* for second order advance)
     pret_vec = temp
     call destroy_vector(temp)
     call destroy_vector(temp2)
     call get_pressures(den_v, te_v, ti_v, p_v, pe_v)
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Temperature--end"
  endif

     
  ! Advance Fields
  ! ==============
  if(iestatic.eq.0) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, " Advancing Fields"

     ! r2matrix_sm * vel(n+1)
     call matvecmult(r2_mat,vel_vec,b1_phi)
     call mult(b1_phi, -1.)

     ! q2matrix_sm * vel(n)
     call matvecmult(q2_mat,veln_vec,b2_phi)
     call add(b1_phi, b2_phi)

     ! d2matrix_sm * phi(n)

     call matvecmult(d2_mat,phi_vec,b2_phi)
     call add(b1_phi, b2_phi)

     ! Inculde density terms
     if(idens.eq.1) then
        call matvecmult(r42_mat,den_vec,b2_phi)
        call add(b1_phi, b2_phi)
        call matvecmult(q42_mat,denold_vec,b2_phi)
        call add(b1_phi, b2_phi)
     end if

     ! Include linear f terms
     if(numvar.ge.2 .and. i3d.eq.1 .and. imp_bf.eq.0) then
        ! b2vector = r15 * bf(n)
        call matvecmult(o2_mat,bf_field(1)%vec,b2_phi)
        call add(b1_phi, b2_phi)
     endif

     call add(b1_phi, q4_vec)
 
     ! Insert boundary conditions
     if(myrank.eq.0 .and. iprint.ge.2) print *, "  inserting bcs"
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     if(calc_matrices.eq.1) then
        call boundary_mag(b1_phi, psi_v, bz_v, bf_v, e_v, s2_mat)
        if(ipressplit.eq.0 .and. numvar.ge.3) then
           if(ipres.eq.1) then
              call boundary_pe(b1_phi, pe_v, s2_mat)
           else
              call boundary_p(b1_phi, pe_v, s2_mat)
           endif
        endif
        call finalize(s2_mat)
     else 
        call boundary_mag(b1_phi, psi_v, bz_v, bf_v, e_v)
        if(ipressplit.eq.0 .and. numvar.ge.3) then
           if(ipres.eq.1) then
              call boundary_pe(b1_phi, pe_v)
           else
              call boundary_p(b1_phi, pe_v)
           endif
        endif
     endif
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_bound = t_bound + tend - tstart
     end if

     if(myrank.eq.0 .and. iprint.ge.2) print *, "  solving"
     ! solve linear system...LU decomposition done first time
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

#ifdef CJ_MATRIX_DUMP
  if(ntime.eq.2) then 
     call write_matrix(s2_mat,'s2_mat')
     call write_vector(b1_phi, 's2_mat_rhs.out')
  endif
#endif 

     call newsolve(s2_mat, b1_phi, jer)
     if(linear.eq.0 .and. iteratephi.eq.0) call clear_mat(s2_mat)


   if(idiff .gt. 0) then
         call add(b1_phi,phi_vec)
   endif

#ifdef CJ_MATRIX_DUMP
  if(ntime.eq.2) then
     call write_vector(b1_phi, 's2_mat_sol.out')
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

     if(myrank.eq.0 .and. iprint.ge.2) print *, "  done solve"

     ! Iterate field solve using re-defined transport coefficients
     if(iteratephi.eq.1) then
      
        if(myrank.eq.0 .and. iprint.ge.1) &
             print *, " Advancing fields... corrector step"
        ! temporarily advance fields to new values
        b2_phi = phi_vec
        phi_vec = b1_phi
!
        call export_time_advance_vectors_split
        ! redefine transport coefficients with new den/pe values
        call lcfs(psi_field(1))
        call define_transport_coefficients
        ! revert fields to old values
        phi_vec = b2_phi
        
        ! recalculate field advance matrix
        ! (advanced velocity variables will be used in defining matrix)
        if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
        call ludefall(0, 0, 0, 0, 1)
        if(myrank.eq.0 .and. itimer.eq.1) then
           call second(tend)
           t_ludefall = t_ludefall + tend - tstart
        endif

        ! r2matrix_sm * vel(n+1)
        call matvecmult(r2_mat,vel_vec,b1_phi)
        call mult(b1_phi, -1.)
   
        ! q2matrix_sm * vel(n)
        call matvecmult(q2_mat,veln_vec,b2_phi)
        call add(b1_phi, b2_phi)

        ! d2matrix_sm * phi(n)
        call matvecmult(d2_mat,phi_vec,b2_phi)
        call add(b1_phi, b2_phi)
      
        ! Include density terms
        if(idens.eq.1) then
           call matvecmult(r42_mat,den_vec,b2_phi)
           call add(b1_phi, b2_phi)
           call matvecmult(q42_mat,denold_vec,b2_phi)
           call add(b1_phi, b2_phi)
        end if

        ! Include linear f terms
        if(numvar.ge.2 .and. i3d.eq.1 .and. imp_bf.eq.0) then
           ! b2vector = r15 * bf(n)
           call matvecmult(o2_mat,bf_field(1)%vec,b2_phi)
           call add(b1_phi, b2_phi)
        endif

        call add(b1_phi, q4_vec)
       
        ! Insert boundary conditions
        if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
        if(calc_matrices.eq.1) then
           call boundary_mag(b1_phi, psi_v, bz_v, bf_v, e_v, s2_mat)
           if(ipressplit.eq.0 .and. numvar.ge.3) then
              if(ipres.eq.1) then
                 call boundary_pe(b1_phi, pe_v, s2_mat)
              else
                 call boundary_p(b1_phi, pe_v, s2_mat)
              endif
           endif

           call finalize(s2_mat)
        else 
           call boundary_mag(b1_phi, psi_v, bz_v, bf_v, e_v)
           if(ipressplit.eq.0 .and. numvar.ge.3) then
              if(ipres.eq.1) then
                 call boundary_pe(b1_phi, pe_v)
              else
                 call boundary_p(b1_phi, pe_v)
              endif
           endif
        endif
        if(myrank.eq.0 .and. itimer.eq.1) then
           call second(tend)
           t_bound = t_bound + tend - tstart
        end if

  
        ! solve linear system...LU decomposition done first time
        if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
        
#ifdef CJ_MATRIX_DUMP
  if(ntime.eq.2) then 
     call write_matrix(s2_mat,'s2_mat')
     call write_vector(b1_phi, 's2_mat_rhs.out')
  endif
#endif 

        call newsolve(s2_mat, b1_phi, jer)
        if(linear.eq.0) call clear_mat(s2_mat)
        
#ifdef CJ_MATRIX_DUMP
  if(ntime.eq.2) then
     call write_vector(b1_phi, 's2_mat_sol.out')
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
     end if     !...on iteratephi

     phi_vec = b1_phi

     ! apply smoothing operators
     ! ~~~~~~~~~~~~~~~~~~~~~~~~~
     call smooth_fields(psi_v) 

!!$     if(ipressplit.eq.0 .and. (numvar.ge.3 .or. ipres.eq.1)) &
!!$          call get_temperatures(den_v, p_v, pe_v, te_v, ti_v)
  end if       !...on iestatic


  if(myrank.eq.0 .and. iprint.ge.1 .and. itimer.eq.1) then
     print *, " split_step: Time solve: ", &
          t_solve_v + t_solve_b + t_solve_n + t_solve_p
     print *, " split_step: Time bcs: ", t_bound

  end if
end subroutine step_split


subroutine subtract_axi

    use basic
    use mesh_mod
    use arrays
    implicit none
    integer:: l,numnodes
    vectype, dimension (dofs_per_node) :: vec_l, axi_l
    type(field_type) :: axi

    call create_field(axi)
    numnodes = owned_nodes()

    do l=1,numnodes
        call get_node_data(u_v,l,vec_l)
        call set_node_data(axi,l,vec_l)
    enddo

    call finalize(axi%vec)
    call sum_vec_planes(axi%vec%data)

    do l=1,numnodes
        call get_node_data(axi,l,axi_l)
        call get_node_data(u_v,l,vec_l)
        vec_l = vec_l - axi_l/nplanes
        call set_node_data(u_v,l,vec_l)
    enddo
    call finalize(u_v%vec)
    if(numvar.lt.2) return
    
    do l=1,numnodes
        call get_node_data(vz_v,l,vec_l)
        call set_node_data(axi,l,vec_l)
    enddo

    call finalize(axi%vec)
    call sum_vec_planes(axi%vec%data)

    do l=1,numnodes
        call get_node_data(axi,l,axi_l)
        call get_node_data(vz_v,l,vec_l)
        vec_l = vec_l - axi_l/nplanes
        call set_node_data(vz_v,l,vec_l)
    enddo
    call finalize(vz_v%vec)
    if(numvar.lt.3) return

    do l=1,numnodes
        call get_node_data(chi_v,l,vec_l)
        call set_node_data(axi,l,vec_l)
    enddo

    call finalize(axi%vec)
    call sum_vec_planes(axi%vec%data)

    do l=1,numnodes
        call get_node_data(axi,l,axi_l)
        call get_node_data(chi_v,l,vec_l)
        vec_l = vec_l - axi_l/nplanes
        call set_node_data(chi_v,l,vec_l)
    enddo
    call finalize(chi_v%vec)

end subroutine subtract_axi
end module time_step_split
