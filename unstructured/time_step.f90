module time_step

  use field

  type(vector_type) :: phi_vec, phiold_vec, phip_vec
  type(vector_type) :: vel_vec, velold_vec, veln_vec, veloldn_vec
  type(vector_type) :: pres_vec, presold_vec
  type(vector_type) :: den_vec, denold_vec
  type(vector_type) :: pret_vec, pretold_vec
  type(vector_type), target :: q4_vec, r4_vec, qp4_vec, qn4_vec

  ! the following pointers point to the vector containing the named field.
  ! set by assign_variables()
  type(field_type) ::   u_v,   uo_v
  type(field_type) ::  vz_v,  vzo_v
  type(field_type) :: chi_v, chio_v
  type(field_type) :: psi_v, psio_v
  type(field_type) ::  bz_v,  bzo_v
  type(field_type) ::  pe_v,  peo_v
  type(field_type) :: den_v, deno_v
  type(field_type) ::   p_v,   po_v
  type(field_type) ::   e_v,   eo_v
  type(field_type) ::  bf_v,  bfo_v
  type(field_type) :: te_v, teo_v
  type(field_type) :: ti_v, tio_v

  ! the offset (relative to the node offset) of the named field within
  ! their respective vectors
  integer :: u_off, vz_off, chi_off
  integer :: psi_off, bz_off, pe_off
  integer :: den_off, p_off
  integer :: bf_off, e_off
  integer :: vecsize_vel, vecsize_phi, vecsize_n, vecsize_p, vecsize_t

  ! temporary vectors
  type(vector_type) :: b1_vel, b2_vel, b1_phi, b2_phi

  integer :: u_i, vz_i, chi_i
  integer :: psi_i, bz_i, pe_i
  integer :: den_i, p_i
  integer :: bf_i, e_i
  integer :: te_i, ti_i

  logical, private :: initialized = .false.

  integer, allocatable :: global_dof_ids_1(:)
  integer, allocatable :: global_dof_ids_row(:), global_dof_ids_col(:)

contains

  subroutine initialize_timestep
    use sparse
    use matrix_mod
    use basic
    implicit none

    select case(isplitstep)
    case(0)
       call set_matrix_index(s1_mat, s1_mat_index)
       call set_matrix_index(d1_mat, d1_mat_index)
       call create_mat(s1_mat, vecsize_vel, vecsize_vel, icomplex, .true.)
       call create_mat(d1_mat, vecsize_vel, vecsize_vel, icomplex, .false.)
#ifdef USERW
       if(eta_wall.ne.0.) then
          call setmatrixrwb(d1_mat_index, 1)
       endif
#endif
#ifdef CJ_MATRIX_DUMP
       print *, "create_mat time_step s1_mat", s1_mat%imatrix     
       print *, "create_mat time_step d1_mat", d1_mat%imatrix     
#endif 
       if(i3d.eq.1) then
          call set_matrix_index(o1_mat, o1_mat_index)
          call create_mat(o1_mat, vecsize_vel, 1, icomplex, .false.)
#ifdef CJ_MATRIX_DUMP
          print *, "create_mat time_step o1_mat", o1_mat%imatrix     
#endif 
       endif

    case(1)
       
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

#ifdef USERW
       if(eta_wall.ne.0.) then
          call setmatrixrwb(d2_mat_index, 1)
       endif
#endif

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
    end select

    initialized = .true.
  end subroutine initialize_timestep

  subroutine finalize_timestep
    use matrix_mod
    use basic
    use sparse
    implicit none

    if(.not.initialized) return

    select case(isplitstep)
    case(0)
       call destroy_mat(s1_mat)
       call destroy_mat(d1_mat)
       if(i3d.eq.1) call destroy_mat(o1_mat)

    case(1)
       
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
       if(i3d.eq.1) call destroy_mat(o2_mat)
       if((ipres.eq.1 .and. numvar.lt.3) .or. ipressplit.gt.0) call destroy_mat(p1_mat)
       
       if(idens.eq.1) then
          call destroy_mat(s8_mat)
          call destroy_mat(d8_mat)
          call destroy_mat(r8_mat)
          call destroy_mat(q8_mat)
       endif

       if(ipres.eq.1) then
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
    end select
  end subroutine finalize_timestep

  

!==========================================================
! assign_variables
! ~~~~~~~~~~~~~~~
! Assigns variables to appropriate vectors for time advance
!==========================================================
  subroutine assign_variables()

    use basic

    implicit none
        
    if(isplitstep.eq.1) then
       
       u_i = 1
       psi_i = 1
       vz_i = 2
       bz_i = 2
       chi_i = 3
       den_i = 1
       if(imp_bf.eq.1) then
          bf_i = numvar - ipressplit + 1
          e_i = numvar - ipressplit + 2
       else
          bf_i = 1
          e_i = numvar - ipressplit + 1
       end if
       if(ipressplit.eq.1) then
            p_i  = 1
            te_i = 1
            if(ipres.eq.1) then
               pe_i = 2
               ti_i = 2
            endif
       else  ! ipresssplit.eq.0
          pe_i = 3
          if(ipres.eq.1) then
             p_i = 1
          end if
          te_i = 1
       endif


       call associate_field(u_v,    vel_vec,      u_i)
       call associate_field(uo_v,   velold_vec,   u_i)
       call associate_field(psi_v,  phi_vec,    psi_i)
       call associate_field(psio_v, phiold_vec, psi_i)
       
       if(numvar.ge.2) then
          call associate_field(vz_v,  vel_vec,    vz_i)
          call associate_field(vzo_v, velold_vec, vz_i)
          call associate_field(bz_v,  phi_vec,    bz_i)
          call associate_field(bzo_v, phiold_vec, bz_i)
       endif

       if(numvar.ge.3) then
          call associate_field(chi_v,  vel_vec,    chi_i)
          call associate_field(chio_v, velold_vec, chi_i)
       endif
       
       if(numvar.ge.3 .or. ipres.eq.1) then
          call associate_field(te_v, pret_vec, te_i)
          call associate_field(teo_v, pretold_vec, te_i)
       endif

       if(ipressplit.eq.1) then
          call associate_field(p_v, pres_vec, p_i)
          call associate_field(po_v, presold_vec, p_i)
          if(ipres.eq.1) then
             call associate_field(pe_v, pres_vec, pe_i)
             call associate_field(peo_v, presold_vec, pe_i)
             call associate_field(ti_v, pret_vec, ti_i)
             call associate_field(tio_v, pretold_vec, ti_i)
          endif
       else   !  on ipressplit.eq.1
          if(numvar.ge.3) then
             call associate_field(pe_v,   phi_vec,     pe_i)
             call associate_field(peo_v,  phiold_vec,  pe_i)
          endif
          if(ipres.eq.1) then
             call associate_field(p_v,  pres_vec,    p_i)
             call associate_field(po_v, presold_vec, p_i)
          end if
       endif  !  on ipressplit.eq.1

       if(idens.eq.1) then
          call associate_field(den_v,  den_vec,    den_i)
          call associate_field(deno_v, denold_vec, den_i)
       end if

       if(jadv.eq.0 .and. i3d.eq.1) then
          call associate_field(e_v, phi_vec, e_i)
       end if

       if(imp_bf.eq.1) then
          call associate_field(bf_v, phi_vec, bf_i)
          call associate_field(bfo_v, phiold_vec, bf_i)
       end if

    else  !  ... on isplitstep.eq.1
       u_i = 1
       psi_i = 2
       vz_i = 3
       bz_i = 4
       chi_i = 5
       pe_i = 6    
       den_i = 2*numvar+1
       if(ipres.eq.1) then
          p_i = 2*numvar+idens+1
       else
          p_i = pe_i
       endif
       if(imp_bf.eq.1) then
          bf_i = 2*numvar+idens+ipres+1
          e_i = 2*numvar+idens+ipres+2
       else
          bf_i = 1
          e_i = 2*numvar+idens+ipres+1
       end if
       te_i = 1
       ti_i = 2

       call associate_field(u_v,    phi_vec,      u_i)
       call associate_field(uo_v,   phiold_vec,   u_i)
       call associate_field(psi_v,  phi_vec,    psi_i)
       call associate_field(psio_v, phiold_vec, psi_i)
       
       if(numvar.ge.2) then
          call associate_field(vz_v,  phi_vec,    vz_i)
          call associate_field(vzo_v, phiold_vec, vz_i)
          call associate_field(bz_v,  phi_vec,    bz_i)
          call associate_field(bzo_v, phiold_vec, bz_i)
       endif

       if(numvar.ge.3) then
          call associate_field(chi_v,  phi_vec,    chi_i)
          call associate_field(chio_v, phiold_vec, chi_i)
          call associate_field(pe_v,   phi_vec,     pe_i)
          call associate_field(peo_v,  phiold_vec,  pe_i)
       endif

       if(ipres.eq.1) then
          call associate_field(p_v,  phi_vec,    p_i)
          call associate_field(po_v, phiold_vec, p_i)
       end if

       if(idens.eq.1) then
          call associate_field(den_v,  phi_vec,    den_i)
          call associate_field(deno_v, phiold_vec, den_i)
       end if

       if(imp_bf.eq.1) then
          call associate_field(bf_v, phi_vec, bf_i)
          call associate_field(bfo_v, phiold_vec, bf_i)
       end if

       if(jadv.eq.0 .and. i3d.eq.1) then
          call associate_field(e_v, phi_vec, e_i)
       end if

    endif !  ... on isplitstep.eq.1
      
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

  end subroutine assign_variables


!============================================================
! ONESTEP
! ~~~~~~~
! advance solution to time ntime+1
!============================================================
subroutine onestep

  use basic
  use diagnostics
  use arrays

  implicit none

  integer :: calc_matrices
  logical, save :: first_time = .true.

  real :: tstart, tend

  ! Determine whether matrices should be re-calculated
  if(first_time &
       .or. (linear.eq.0 .and. mod(ntime,nskip).eq.0) &
       .or. (integrator.eq.1 .and. ntime.eq.1)) then
     calc_matrices = 1
  else
     calc_matrices = 0
  endif

  ! calculate matrices for time advance
  if(calc_matrices.eq.1) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Defining matrices"
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

     ! in linear case, eliminate second-order terms from matrix
     call ludefall(1-istatic, idens, ipres, ipressplit, 1-iestatic)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_ludefall = t_ludefall + tend - tstart
     endif
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Done defining matrices."
  endif


  ! copy field data to time-advance vectors
  if(myrank.eq.0 .and. iprint.ge.2) print *, "Importing time advance vectors.."
  call import_time_advance_vectors

  ! advance time
  if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing time..."
  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
  if(isplitstep.eq.1) then
     call split_step(calc_matrices)
  else
     call unsplit_step(calc_matrices)
  end if
  if(myrank.eq.0 .and. itimer.eq.1) then 
     call second(tend)
     if(iprint.ge.1) print *, "Time spent in *_step: ", tend-tstart
  end if

  time = time + dt
  if(ntime.gt.1 .and. linear.eq.0) call variable_timestep

  ! copy time advance vectors to field data
  if(myrank.eq.0 .and. iprint.ge.2) print *, "Exporting time advance vectors.."
  call export_time_advance_vectors
!

  ! Calculate all quantities derived from basic fields
  call derived_quantities(1)


  ! Conserve toroidal flux
  if(iconstflux.eq.1 .and. numvar.ge.2) then
     call conserve_flux
     tflux = tflux + gbound*area
  endif


  first_time = .false.

end subroutine onestep


!======================================================================
! import_time_advance_vectors
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
! construct vectors for time-advance
!======================================================================
subroutine import_time_advance_vectors
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


end subroutine import_time_advance_vectors


!======================================================================
! export_time_advance_vectors
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
! extract values from time-advance vectors
!======================================================================
subroutine export_time_advance_vectors
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
     te_field(1) = te_v
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

end subroutine export_time_advance_vectors

!============================================================
! SPILT_STEP
! ~~~~~~~~~~
! advance fields using split time step
!============================================================
subroutine split_step(calc_matrices)

  use p_data
  use basic
  use arrays
  use sparse
  use diagnostics
  use matrix_mod
  use boundary_conditions
  use mesh_mod 

  implicit none

#ifdef CJ_MATRIX_DUMP
  integer :: counter
#endif


  integer, intent(in) :: calc_matrices
  real :: tstart, tend, t_bound
  integer :: jer, i, numnodes
  type(vector_type) :: temp, temp2

  type(field_type) :: phip_1, phip_2, phip_3, temp_field_1, temp_field_2
  call associate_field(phip_1, phip_vec, 1)
  if(numvar.ge.2) call associate_field(phip_2, phip_vec, 2)
  if(numvar.ge.3 .and. ipressplit.eq.0) call associate_field(phip_3, phip_vec, 3)

  t_bound = 0.

  ! Store current-time velocity matrices for use in field advance
  veln_vec = vel_vec
  veloldn_vec = velold_vec


  if(istatic.eq.0) then

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
     if(numvar.ge.2 .and. i3d.eq.1) then
        call matvecmult(o1_mat,bf_field(1)%vec,b2_vel)
        call add(b1_vel, b2_vel)
     endif

     call add(b1_vel, r4_vec)
  
     ! apply boundary conditions
     if(myrank.eq.0 .and. iprint.ge.2) print *, "  inserting bcs"
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     if(calc_matrices.eq.1) then
        call boundary_vel(b1_vel, s1_mat)
        call finalize(s1_mat)
     else
        call boundary_vel(b1_vel)
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
  
#ifdef USESCOREC
     if(integrator.eq.1 .and. ntime.gt.1) then
        b1_vel%data = (2.*b1_vel%data - velold_vec%data)/3.
     endif
#endif
  
     !.....new velocity solution at time n+1 (or n* for second order advance)
     velold_vec = vel_vec
     vel_vec = b1_vel
   
     ! apply smoothing operators
     ! ~~~~~~~~~~~~~~~~~~~~~~~~~
     call smooth_velocity(u_v, chi_v)
  else
     velold_vec = vel_vec
  end if
  
  ! Advance Density
  ! ===============
  if(idens.eq.1) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, " Advancing density"
     
     call create_vector(temp, vecsize_n)
     call create_vector(temp2, vecsize_n)

     ! r8matrix_sm * vel(n+1)
     call matvecmult(r8_mat,vel_vec,temp)

     ! r8matrix_sm * vel(n-1)
     if(integrator.eq.1 .and. ntime.gt.1) then
        call matvecmult(r8_mat,veloldn_vec,temp2)
        call mult(temp, 1.5)
        call mult(temp2, 0.5)
        call add(temp, temp2)
     endif

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
        call boundary_den(temp, s8_mat)
        call finalize(s8_mat)
     else
        call boundary_den(temp)
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
#ifdef USESCOREC
     if(integrator.eq.1 .and. ntime.gt.1) then
        temp%data = (2.*temp%data - denold_vec%data)/3.
     endif
#endif
     denold_vec = den_vec
     den_vec = temp
     call destroy_vector(temp)
     call destroy_vector(temp2)

     if(irecalc_eta.eq.1) then
        call export_time_advance_vectors
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
        
     ! r9matrix_sm * vel(n-1)
     if(integrator.eq.1 .and. ntime.gt.1) then
        call matvecmult(r9_mat,veloldn_vec,temp2)
        call mult(temp, 1.5)
        call mult(temp2, 0.5)
        call add(temp, temp2)
     endif
     
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
        call boundary_pres(temp, s9_mat)
        call finalize(s9_mat)
     else
        call boundary_pres(temp)
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
#ifdef USESCOREC
     if(integrator.eq.1 .and. ntime.gt.1) then
        temp%data = (2.*temp%data - presold_vec%data)/3.
     endif
#endif
     presold_vec = pres_vec
     pres_vec = temp
     call destroy_vector(temp)
     call destroy_vector(temp2)
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Pressure--before get_temperatures"
     call get_temperatures
  endif

  ! Advance Temperature
  ! ================
  if(ipressplit.eq.1 .and. itemp.eq.1) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Temperature"

     call create_vector(temp, vecsize_p)
     call create_vector(temp2, vecsize_p)
     
     ! r9matrix_sm * vel(n+1)
     call matvecmult(r9_mat,vel_vec,temp)
        
     ! r9matrix_sm * vel(n-1)
     if(integrator.eq.1 .and. ntime.gt.1) then
        call matvecmult(r9_mat,veloldn_vec,temp2)
        call mult(temp, 1.5)
        call mult(temp2, 0.5)
        call add(temp, temp2)
     endif
     
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
     
     call add(temp, qp4_vec)
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Temperature -- before boundary_temp"
     
     ! Insert boundary conditions
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     if(calc_matrices.eq.1) then
        call boundary_temp(temp, s9_mat)
        call finalize(s9_mat)
     else
        call boundary_temp(temp)
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

     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_solve_p = t_solve_p + tend - tstart
     endif
     if(jer.ne.0) then
        write(*,*) 'Error in pressure solve', jer
        call safestop(29)
     endif
     
     ! new field solution at time n+1 (or n* for second order advance)
#ifdef USESCOREC
     if(integrator.eq.1 .and. ntime.gt.1) then
        temp%data = (2.*temp%data - pretold_vec%data)/3.
     endif
#endif
     pretold_vec = pret_vec
     pret_vec = temp
     call destroy_vector(temp)
     call destroy_vector(temp2)
     call get_pressures
     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Temperature--end"
  endif

     
  ! Advance Fields
  ! ==============
  if(iestatic.eq.0) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, " Advancing Fields"

     ! r2matrix_sm * vel(n+1)
     call matvecmult(r2_mat,vel_vec,b1_phi)

     ! r2matrix_sm * vel(n-1)
     if(integrator.eq.1 .and. ntime.gt.1) then
        call matvecmult(r2_mat,veloldn_vec,b2_phi)
        call mult(b1_phi, 1.5)
        call mult(b2_phi, 0.5)
        call add(b1_phi, b2_phi)
     endif

     ! q2matrix_sm * vel(n)
     call mult(b1_phi, -1.)
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
        call boundary_mag(b1_phi, s2_mat)
        call finalize(s2_mat)
     else 
        call boundary_mag(b1_phi)
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
  
     ! new field solution at time n+1 (or n* for second order advance
#ifdef USESCOREC
     if(integrator.eq.1 .and. ntime.gt.1) then
        b1_phi%data = (2.*b1_phi%data - phiold_vec%data)/3.
     endif
#endif

     ! Iterate field solve using re-defined transport coefficients
     if(iteratephi.eq.1) then
      
        if(myrank.eq.0 .and. iprint.ge.1) &
             print *, " Advancing fields... corrector step"
        ! temporarily advance fields to new values
        b2_phi = phi_vec
        phi_vec = b1_phi
!
        call export_time_advance_vectors
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
   
        ! r2matrix_sm * vel(n-1)
        if(integrator.eq.1 .and. ntime.gt.1) then
           call matvecmult(r2_mat,veloldn_vec,b2_phi)
           call mult(b1_phi, 1.5)
           call mult(b2_phi, 0.5)
           call add(b1_phi, b2_phi)
        endif

        ! q2matrix_sm * vel(n)
        call mult(b1_phi, -1.)
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
           call boundary_mag(b1_phi, s2_mat)
           call finalize(s2_mat)
        else 
           call boundary_mag(b1_phi)
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
  
        ! new field solution at time n+1 (or n* for second order advance)
#ifdef USESCOREC
        if(integrator.eq.1 .and. ntime.gt.1) then
           b1_phi%data = (2.*b1_phi%data - phiold_vec%data)/3.
        endif
#endif

     end if     !...on iteratephi

     phiold_vec = phi_vec
     phi_vec = b1_phi

     ! apply smoothing operators
     ! ~~~~~~~~~~~~~~~~~~~~~~~~~
     call smooth_fields(psi_v) 

     if(ipressplit.eq.0 .and. numvar.ge.3) call get_temperatures

  end if       !...on iestatic


  if(myrank.eq.0 .and. iprint.ge.1 .and. itimer.eq.1) then
     print *, " split_step: Time solve: ", &
          t_solve_v + t_solve_b + t_solve_n + t_solve_p
     print *, " split_step: Time bcs: ", t_bound

  end if
  
end subroutine split_step


!============================================================
! UNSPILT_STEP
! ~~~~~~~~~~~~
! advance fields using unsplit time step
!============================================================
subroutine unsplit_step(calc_matrices)

  use boundary_conditions
  use p_data
  use basic
  use arrays
  use sparse
  use diagnostics
  use matrix_mod

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
     call boundary_mag(b1_phi, s1_mat)
     call boundary_vel(b1_phi, s1_mat)
     if(idens.eq.1) call boundary_den(b1_phi, s1_mat)
     if(ipres.eq.1) call boundary_pres(b1_phi, s1_mat)
     call finalize(s1_mat)
  else
     call boundary_mag(b1_phi)
     call boundary_vel(b1_phi)
     if(idens.eq.1) call boundary_den(b1_phi)
     if(ipres.eq.1) call boundary_pres(b1_phi)
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
#ifdef USESCOREC
  if(integrator.eq.1 .and. ntime.gt.1) then
     b1_phi%data = (2.*b1_phi%data - phiold_vec%data)/3.
  endif
#endif

  ! apply smoothing operators
  call smooth_velocity(u_v, chi_v)
  call smooth_fields(psi_v)

  if(iteratephi.eq.1 .and. iestatic.eq.0) then
     if(myrank.eq.0 .and. iprint.ge.1) print *, "secondary advance..."

     ! temporarily advance fields to new values
     b2_phi = phi_vec
     phi_vec = b1_phi
     call export_time_advance_vectors
     ! redefine transport coefficients with new den/pe values
     call lcfs(psi_field(1))
     call define_transport_coefficients
     ! revert fields to old values
     phi_vec = b2_phi
        
     ! recalculate field advance matrix
     ! (advanced velocity variables will be used in defining matrix)
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call ludefall(1-istatic, idens, ipres, ipressplit, 1-iestatic)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_ludefall = t_ludefall + tend - tstart
     endif

     ! vtemp = d1matrix_sm * phi(n)
     call matvecmult(d1_mat,phi_vec,b1_phi) 
     call add(b1_phi, q4_vec)

     ! Include linear f terms
     if(numvar.ge.2 .and. i3d.eq.1) then
        ! b2vector = r15 * bf(n)
     
        ! make a larger vector that can be multiplied by a vecsize matrix
        call matvecmult(o1_mat,bf_field(1)%vec,b2_phi)
        call add(b1_phi, b2_phi)
     endif
    
     ! Insert boundary conditions
     if(calc_matrices.eq.1) then
        call boundary_mag(b1_phi, s1_mat)
        call boundary_vel(b1_phi, s1_mat)
        if(idens.eq.1) call boundary_den(b1_phi, s1_mat)
        if(ipres.eq.1) call boundary_pres(b1_phi, s1_mat)
        call finalize(s1_mat)
     else 
        call boundary_mag(b1_phi)
        call boundary_vel(b1_phi)
        if(idens.eq.1) call boundary_den(b1_phi)
        if(ipres.eq.1) call boundary_pres(b1_phi)
     endif
     
     ! solve linear system...LU decomposition done first time
     if(myrank.eq.0 .and. iprint.ge.1) print *, "solving.."
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

#ifdef CJ_MATRIX_DUMP
     if(counter.le.0) then 
        call write_matrix(s1_mat,'s1_mat')
        call write_vector(b1_phi, 's1_mat_rhs.out')
     endif
#endif 

     call newsolve(s1_mat, b1_phi, jer)
     if(linear.eq.0) call clear_mat(s1_mat)

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
     
     ! new field solution at time n+1
#ifdef USESCOREC
     if(integrator.eq.1 .and. ntime.gt.1) then
        b1_phi%data = (2.*b1_phi%data - phiold_vec%data)/3.
     endif
#endif
     
     ! apply smoothing operators
     call smooth_velocity(u_v, chi_v)
     call smooth_fields(psi_v)
  endif
     
  phiold_vec = phi_vec
  phi_vec = b1_phi


  if(myrank.eq.0 .and. iprint.ge.1) print *, "Done solving matrix equation."
end subroutine unsplit_step


subroutine calc_ni(ni_field, n0_field, n1_field)

  use field
  use mesh_mod

  implicit none

  type(field_type), intent(inout) :: ni_field
  type(field_type), intent(in) :: n0_field, n1_field

  integer :: numnodes, inode
  vectype, dimension(dofs_per_node) :: ninv, n0, n1

  numnodes = owned_nodes()

  do inode=1,numnodes
     call get_node_data(n0_field, inode, n0)
     call get_node_data(n1_field, inode, n1)

     ninv(1) = -n1(1)/n0(1)**2
     ninv(2) = &
          -(n1(2)*n0(1) - 2.*n1(1)*n0(2))/n0(1)**3
     ninv(3) = &
          -(n1(3)*n0(1) - 2.*n1(1)*n0(3))/n0(1)**3
     ninv(4) = &
          -n1(4)/n0(1)**2 &
          +2.*(2.*n1(2)*n0(2) &
          +   n1(1)*n0(4))/n0(1)**3 &
          -6.*n1(1)*n0(2)**2/n0(1)**4
     ninv(5) = &
          -n1(5)/n0(1)**2 &
          +2.*(n1(2)*n0(3) + n1(3)*n0(2) &
          +n1(1)*n0(5))/n0(1)**3 &
          -6.*n1(1)*n0(2)*n0(3)/n0(1)**4
     ninv(6) = &
          -n1(6)/n0(1)**2 &
          +2.*(2.*n1(3)*n0(3) &
          +   n1(1)*n0(6))/n0(1)**3 &
          -6.*n1(1)*n0(3)**2/n0(1)**4

     call set_node_data(ni_field, inode, ninv)
  end do
end subroutine calc_ni

subroutine calc_b2i(b2i_field, psi0_field, psi1_field, b0_field, b1_field)
  use basic
  use field
  implicit none

  type(field_type), intent(inout) :: b2i_field
  type(field_type), intent(in) :: psi0_field, psi1_field, b0_field, b1_field

  vectype, dimension(dofs_per_node) :: vec, psi0, psi1, b0, b1
  vectype, dimension(dofs_per_node) :: b2i
  vectype :: b20
  real :: x, phi, z
  integer :: inode, numnodes

  numnodes = owned_nodes()

  do inode=1,numnodes
     call get_node_data(psi0_field, inode, psi0)
     call get_node_data(psi1_field, inode, psi1)
     call get_node_data(b0_field, inode, b0)
     call get_node_data(b1_field, inode, b1)

     if(itor.eq.0) then
        x = 1.
     else
        call get_node_pos(inode, x, phi, z)
     endif

     b20 = psi0(2)**2 + psi0(3)**2 + b0(1)**2

     b2i(1) = x**2/b20
     b2i(2) = 2.*x/b20 &
          - (2.*x**2/b20**2)*(psi0(2)*psi0(4) + psi0(3)*psi0(5) + b0(1)*b0(2))
     b2i(3) = &
          - (2.*x**2/b20**2)*(psi0(2)*psi0(5) + psi0(3)*psi0(6) + b0(1)*b0(3))
     b2i(4:6) = 0.
     
     vec(1) = -2.*b2i(1)**2/x**2 * &
          (b0(1)*b1(1) + psi0(2)*psi1(2) + psi0(3)*psi1(3))
     vec(2) = -2.*b2i(2)**2/x**2 * &
          (b0(1)*b1(1) + psi0(2)*psi1(2) + psi0(3)*psi1(3)) &
          -2.*b2i(1)**2/x**2 * &
          (b0(2)*b1(1) + psi0(4)*psi1(2) + psi0(5)*psi1(3)  &
          +b0(1)*b1(2) + psi0(2)*psi1(4) + psi0(3)*psi1(5)) &
          +4.*b2i(1)**2/x**3 * &
          (b0(1)*b1(1) + psi0(2)*psi1(2) + psi0(3)*psi1(3))
     vec(3) = -2.*b2i(3)**2/x**2 * &
          (b0(1)*b1(1) + psi0(2)*psi1(2) + psi0(3)*psi1(3)) &
          -2.*b2i(1)**2/x**2 * &
          (b0(3)*b1(1) + psi0(5)*psi1(2) + psi0(6)*psi1(3)  &
          +b0(1)*b1(3) + psi0(2)*psi1(5) + psi0(3)*psi1(6))
     vec(4:6) = 0.

     call set_node_data(b2i_field, inode, vec)
  end do
end subroutine calc_b2i


subroutine get_den_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = 0
  if(inograd_n.eq.1) ibound = ior(ibound, BOUNDARY_NEUMANN)
  if(iconst_n.eq.1) ibound = ior(ibound, BOUNDARY_DIRICHLET)
  
  call get_boundary_mask(itri, ibound, imask)
end subroutine get_den_mask

subroutine get_temp_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = 0
  if(inograd_t.eq.1) ibound = ior(ibound, BOUNDARY_NEUMANN)
  if(iconst_t.eq.1) ibound = ior(ibound, BOUNDARY_DIRICHLET)
  
  call get_boundary_mask(itri, ibound, imask)
end subroutine get_temp_mask

subroutine get_pres_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = 0
  if(inograd_p.eq.1) ibound = ior(ibound, BOUNDARY_NEUMANN)
  if(iconst_p.eq.1) ibound = ior(ibound, BOUNDARY_DIRICHLET)
  
  call get_boundary_mask(itri, ibound, imask)
end subroutine get_pres_mask


subroutine get_flux_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  if(eta_wall.eq.0.) then 
     ibound = BOUNDARY_DIRICHLET
  else
     ibound = BOUNDARY_RESISTIVE_WALL
  end if
  if(inocurrent_tor.eq.1) ibound = ior(ibound, BOUNDARY_LAPLACIAN)
  if(inocurrent_norm.eq.1) then
     if(i3d.eq.1) then
        ibound = ior(ibound, BOUNDARY_NEUMANNP)
     endif
  endif

  call get_boundary_mask(itri, ibound, imask)
end subroutine get_flux_mask

subroutine get_bz_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = 0
!  if(eta_wall.ne.0.) ibound = ior(ibound, BOUNDARY_RESISTIVE_WALL)
  if(inocurrent_pol.eq.1) ibound = ior(ibound, BOUNDARY_NEUMANN)
  if(inocurrent_norm.eq.1 .or. iconst_bz.eq.1) &
       ibound = ior(ibound, BOUNDARY_DIRICHLET)
  call get_boundary_mask(itri, ibound, imask)
end subroutine get_bz_mask

subroutine get_q_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = BOUNDARY_DIRICHLET
  call get_boundary_mask(itri, ibound, imask)
end subroutine get_q_mask

subroutine get_bf_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = BOUNDARY_DIRICHLET
  call get_boundary_mask(itri, ibound, imask)
end subroutine get_bf_mask

subroutine get_vor_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = 0
  if(inonormalflow.eq.1) ibound = ior(ibound, BOUNDARY_DIRICHLET)
  if(inoslip_pol.eq.1)   ibound = ior(ibound, BOUNDARY_NEUMANN)
  if(vor_bc.eq.1)        ibound = ior(ibound, BOUNDARY_LAPLACIAN)
  call get_boundary_mask(itri, ibound, imask)
end subroutine get_vor_mask

subroutine get_vz_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = 0
  if(inoslip_tor.eq.1)   ibound = ior(ibound, BOUNDARY_DIRICHLET)
  if(inostress_tor.eq.1) ibound = ior(ibound, BOUNDARY_NEUMANN)
  call get_boundary_mask(itri, ibound, imask)
end subroutine get_vz_mask

subroutine get_chi_mask(itri, imask)
  use element
  use basic
  use boundary_conditions
  implicit none
  integer, intent(in) :: itri
  integer, intent(out), dimension(dofs_per_element) :: imask
  integer :: ibound

  ibound = 0
  if(inonormalflow.eq.1) ibound = ior(ibound, BOUNDARY_NEUMANN)
  if(inoslip_pol.eq.1)   ibound = ior(ibound, BOUNDARY_DIRICHLET)
  if(com_bc.eq.1)        ibound = ior(ibound, BOUNDARY_LAPLACIAN)
  call get_boundary_mask(itri, ibound, imask)
end subroutine get_chi_mask




!=======================================================
! boundary_vel
! ~~~~~~~~~~~~
!
! sets boundary conditions for velocity fields
!=======================================================
subroutine boundary_vel(rhs, mat)
  use basic
  use arrays
  use matrix_mod
  use boundary_conditions

  implicit none
  
  type(vector_type), intent(inout) :: rhs
  type(matrix_type), optional :: mat

  vectype, dimension(dofs_per_node) :: temp 
  real :: normal(2), curv, x, z
  integer :: i, izone, izonedim, numnodes
  integer :: i_u, i_vz, i_chi
  logical :: is_boundary

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_vel called"

  numnodes = owned_nodes()
  do i=1, numnodes

     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     i_u = node_index(u_v, i)
     if(numvar.ge.2) i_vz = node_index(vz_v, i)
     if(numvar.ge.3) i_chi = node_index(chi_v, i)
       
     ! no normal flow
     if(inonormalflow.eq.1) then
        temp = 0.
        call set_dirichlet_bc(i_u,rhs,temp,normal,curv,izonedim,mat)
        if(numvar.ge.3) then
           call set_normal_bc(i_chi,rhs,temp,normal,curv,izonedim,mat)
        endif
     end if
     
     ! no poloidal slip
     if(inoslip_pol.eq.1) then
        temp = 0.
        call set_normal_bc(i_u,rhs,temp,normal,curv,izonedim,mat)
        if(numvar.ge.3) then
           call set_dirichlet_bc(i_chi,rhs,temp,normal,curv,izonedim,mat)
        endif
     end if

     ! toroidal velocity
     if(numvar.ge.2) then
        ! no slip
        if(inoslip_tor.eq. 1) then
           call get_node_data(vz_field(1), i, temp)
           call set_dirichlet_bc(i_vz,rhs,temp,normal,curv,izonedim,mat)
        end if
        
        ! no toroidal stress
        if(inostress_tor.eq.1) then
           temp = 0.
           call set_normal_bc(i_vz,rhs,temp,normal,curv,izonedim,mat)
        end if
     endif
       
     ! no vorticity
     if(vor_bc.eq.1) then
        temp = 0.
        call set_laplacian_bc(i_u,rhs,temp,normal,curv,izonedim,-x,mat)
     endif

     ! no compression
     if(com_bc.eq.1 .and. numvar.ge.3) then
        select case(ivform)
        case(0)
           call set_laplacian_bc(i_chi,rhs,temp,normal,curv,izonedim,x,mat)
        case(1)
           call set_laplacian_bc(i_chi,rhs,temp,normal,curv,izonedim,-x,mat)
        end select
     endif
  end do
end subroutine boundary_vel


!=======================================================
! boundary_mag
! ~~~~~~~~~~~~
!
! sets boundary conditions for magnetic fields
! and electron pressure 
!=======================================================
subroutine boundary_mag(rhs, mat)
  use math
  use basic
  use arrays
  use matrix_mod
  use boundary_conditions

  implicit none
  
  type(vector_type), intent(inout) :: rhs
  type(matrix_type), optional :: mat

  vectype, dimension(dofs_per_node) :: temp, temp2, temp3
  real :: normal(2), curv, x, z
  integer :: i, izone, izonedim,  numnodes
  integer :: i_psi, i_bz, i_pe, i_e, i_bf
  logical :: is_boundary

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_mag called"

  numnodes = owned_nodes()
  do i=1, numnodes

     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     i_psi = node_index(psi_v, i)
     if(numvar.ge.2) i_bz = node_index(bz_v, i)
     if(numvar.ge.3 .and. ipressplit.eq.0) i_pe = node_index(pe_v, i)
     if(jadv.eq.0 .and. i3d.eq.1) i_e = node_index(e_v, i)
     if(imp_bf.eq.1) i_bf = node_index(bf_v, i)

     ! constant normal field = -n.grad(psi)/R - n.grad(f')
     if(iconst_bn.eq.1) then
        call get_node_data(psi_field(1), i, temp)
        ! add loop voltage
        if(igauge.eq.0) temp(1) = temp(1) + dt*vloop/twopi
        call set_dirichlet_bc(i_psi,rhs,temp,normal,curv,izonedim,mat)
     end if

     ! no toroidal current = -Delta*(psi)/R
     if(inocurrent_tor.eq.1) then       
        temp = 0.
        call set_laplacian_bc(i_psi,rhs,temp,normal,curv,izonedim,-x,mat)
     end if

     ! no tangential current = n.Grad(F)/R + t.Grad(psi')/R^2
     if(inocurrent_pol.eq.1 .and. numvar.ge.2) then
        call get_node_data(bz_field(1), i, temp)
        call set_normal_bc(i_bz,rhs,temp,normal,curv,izonedim,mat)
     end if

     ! no normal current = n.Grad(psi')/R^2 - t.Grad(F)/R
     if(inocurrent_norm.eq.1) then
        if(i3d.eq.1) then
!        if(numvar.ge.2) then
!           call get_node_data(bz_field(1), i, temp2)
!        else
!           temp2 = 0.
!        endif
           call get_node_data(psi_field(1), i, temp)
           call set_normalp_bc(i_psi,rhs,temp,normal,curv,izonedim,mat)
        endif

        if(numvar.ge.2) then
           call get_node_data(bz_field(1), i, temp)
           call set_dirichlet_bc(i_bz,rhs,temp,normal,curv,izonedim,mat)
        endif
     else if(iconst_bz.eq.1 .and. numvar.ge.2) then
        call get_node_data(bz_field(1), i, temp)
        call set_dirichlet_bc(i_bz,rhs,temp,normal,curv,izonedim,mat)
     endif

     if(numvar.ge.3 .and. ipressplit.eq.0) then
        if(inograd_p.eq.1) then
           temp = 0.
           call set_normal_bc(i_pe,rhs,temp,normal,curv,izonedim,mat)
        end if
        if(ipres.eq.1) then
           call get_node_data(pe_field(1), i, temp)
        else
           call get_node_data(p_field(1), i, temp)
        end if
        if(iconst_p.eq.1) then
           call set_dirichlet_bc(i_pe,rhs,temp,normal,curv,izonedim,mat)
        else if(iconst_t.eq.1) then
           call get_node_data(den_v, i, temp2)
           call get_node_data(den_field(1), i, temp3)
           temp = temp*temp2(1)/temp3(1)
           call set_dirichlet_bc(i_pe,rhs,temp,normal,curv,izonedim,mat)
        end if
     endif

     if(jadv.eq.0 .and. i3d.eq.1) then
        ! electrostatic potential
        temp = 0.
        call set_dirichlet_bc(i_e,rhs,temp,normal,curv,izonedim,mat)
     endif

     if(imp_bf.eq.1) then
        temp = 0.
        call set_dirichlet_bc(i_bf,rhs,temp,normal,curv,izonedim,mat)
     end if
  end do


  if(eta_wall.ne.0.) call insert_resistive_wall(rhs,mat)
end subroutine boundary_mag


!=======================================================
! boundary_den
! ~~~~~~~~~~~~
!
! sets boundary conditions for density
!=======================================================
subroutine boundary_den(rhs, mat)
  use basic
  use arrays
  use matrix_mod
  use boundary_conditions
  implicit none
  
  type(vector_type) :: rhs
  type(matrix_type), optional :: mat
  
  integer :: i, izone, izonedim, numnodes
  real :: normal(2), curv, x,z
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

  integer :: i_n

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_den called"

  numnodes = owned_nodes()
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     i_n = node_index(den_v, i)

     if(inograd_n.eq.1) then
        temp = 0.
        call set_normal_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
     end if
     if(iconst_n.eq.1) then
        call get_node_data(den_field(1), i, temp)
        call set_dirichlet_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
     end if

  end do

end subroutine boundary_den

!=======================================================
! boundary_te
! ~~~~~~~~~~~~
!
! sets boundary conditions for electron temperature
!=======================================================
subroutine boundary_te(rhs, mat)
  use basic
  use arrays
  use matrix_mod
  use boundary_conditions
  implicit none
  
  type(vector_type) :: rhs
  type(matrix_type), optional :: mat
  
  integer :: i, izone, izonedim, numnodes
  real :: normal(2), curv, x,z
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

  integer :: i_n

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_te called"

  numnodes = owned_nodes()
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     i_n = node_index(te_v, i)

     if(inograd_t.eq.1) then
        temp = 0.
        call set_normal_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
     end if
     if(iconst_t.eq.1) then
        call get_node_data(te_field(1), i, temp)
        call set_dirichlet_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
     end if

  end do

end subroutine boundary_te
subroutine boundary_ti(rhs, mat)
  use basic
  use arrays
  use matrix_mod
  use boundary_conditions
  implicit none
  
  type(vector_type) :: rhs
  type(matrix_type), optional :: mat
  
  integer :: i, izone, izonedim, numnodes
  real :: normal(2), curv, x,z
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

  integer :: i_n

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_ti called"

  numnodes = owned_nodes()
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     i_n = node_index(ti_v, i)

     if(inograd_t.eq.1) then
        temp = 0.
        call set_normal_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
     end if
     if(iconst_t.eq.1) then
        call get_node_data(ti_field(1), i, temp)
        call set_dirichlet_bc(i_n,rhs,temp,normal,curv,izonedim,mat)
     end if

  end do

end subroutine boundary_ti


!=======================================================
! boundary_pres
! ~~~~~~~~~~~~~
!
! sets boundary conditions for pressure
!=======================================================
subroutine boundary_pres(rhs, mat)
  use basic
  use arrays
  use matrix_mod
  use boundary_conditions
  implicit none
  
  type(vector_type) :: rhs
  type(matrix_type), optional :: mat
  
  integer :: i, izone, izonedim, numnodes
  real :: normal(2), curv, x, z
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

  integer :: i_p, i_pe

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.2) print *, "boundary_pres called"

  numnodes = owned_nodes()
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     i_p = node_index(p_v, i)

     if(inograd_p.eq.1) then
        temp = 0.
        call set_normal_bc(i_p,rhs,temp,normal,curv,izonedim,mat)
     end if
     if(iconst_p.eq.1) then
        call get_node_data(p_field(1), i, temp)
        call set_dirichlet_bc(i_p,rhs,temp,normal,curv,izonedim,mat)
     end if
 
     if(ipressplit.eq.1 .and. ipres.eq.1) then
        i_pe = node_index(pe_v, i)
        if(inograd_p.eq.1) then
           temp = 0.
           call set_normal_bc(i_pe,rhs,temp,normal,curv,izonedim,mat)
        end if
        if(iconst_p.eq.1) then
           call get_node_data(pe_field(1), i, temp)
           call set_dirichlet_bc(i_pe,rhs,temp,normal,curv,izonedim,mat)
        end if
     endif

  end do

end subroutine boundary_pres
subroutine boundary_temp(rhs, mat)
  use basic
  use arrays
  use matrix_mod
  use boundary_conditions
  implicit none
  
  type(vector_type) :: rhs
  type(matrix_type), optional :: mat
  
  integer :: i, izone, izonedim, numnodes
  real :: normal(2), curv, x, z
  logical :: is_boundary
  vectype, dimension(dofs_per_node) :: temp

  integer :: i_te, i_ti

  if(iper.eq.1 .and. jper.eq.1) return
  if(myrank.eq.0 .and. iprint.ge.1) print *, "boundary_temp called"

  numnodes = owned_nodes()
  do i=1, numnodes
     call boundary_node(i,is_boundary,izone,izonedim,normal,curv,x,z)
     if(.not.is_boundary) cycle

     i_te = node_index(te_v, i)

     if(inograd_p.eq.1) then
        temp = 0.
        call set_normal_bc(i_te,rhs,temp,normal,curv,izonedim,mat)
     end if
     if(iconst_p.eq.1) then
        call get_node_data(te_field(1), i, temp)
        call set_dirichlet_bc(i_te,rhs,temp,normal,curv,izonedim,mat)
     end if
 
     if(ipressplit.eq.1 .and. ipres.eq.1) then
        i_ti = node_index(ti_v, i)
        if(inograd_p.eq.1) then
           temp = 0.
           call set_normal_bc(i_ti,rhs,temp,normal,curv,izonedim,mat)
        end if
        if(iconst_p.eq.1) then
           call get_node_data(ti_field(1), i, temp)
           call set_dirichlet_bc(i_ti,rhs,temp,normal,curv,izonedim,mat)
        end if
     endif

  end do

end subroutine boundary_temp

  subroutine insert_resistive_wall(rhs,mat)
    use basic
    use arrays
    use sparse
    use vacuum_interface
    use matrix_mod
    use boundary_conditions

    implicit none

    include 'mpif.h'
    
    type(vector_type), intent(inout) :: rhs
    type(matrix_type), optional :: mat
    
    type(vector_type) :: tempout
    type(vector_type), save :: tempout0
    integer :: i, j, ii, jj, ir, m, n, ip, jp
    integer :: irow(2,dofs_per_node), icol(3,dofs_per_node)
    integer :: jrow(2,dofs_per_node), jcol(3,dofs_per_node)
    integer :: ibegin, jbegin, ierr
    real :: fac, thimprw
    real :: xii, zii, xjj, normi(2), normj(2), curvi, curvj
    integer, parameter :: num_rows = 1
    integer :: rows(3)
    
    logical :: is_boundary
    integer :: izone, izonedim

    logical, save :: first_time = .true.
    logical :: use_resistive_bz
    integer, allocatable :: buf(:)

    vectype, dimension(dofs_per_node,dofs_per_node,2,3) :: ssi, ddi
    vectype, dimension(dofs_per_node,dofs_per_node,2,3) :: ssj, ddj, rrj
    vectype, dimension(dofs_per_node,dofs_per_node,2,3) :: ss, dd
    vectype, dimension(dofs_per_node) :: temp, temp0

    type(vector_type) :: mag
    type(field_type) :: mag_psi, mag_bz, mag_bf

    real :: mp, a, k
    mp = 2.
    a = 1.
    k = ntor/100.

    use_resistive_bz = numvar.ge.2 .and. iconst_bz.ne.1
!    use_resistive_bz = .false.

    rows(1) = 1
    rows(2) = 3

    if(myrank.eq.0 .and. iprint.ge.1) print *, 'Inserting resistive wall'

    ! create magnetic field vector
    call create_vector(mag, 3)
    call associate_field(mag_psi, mag, 1)
    mag_psi = psi_v
    call add(mag_psi, external_psi_field, -1.)
    if(numvar.ge.2) then
       call associate_field(mag_bz, mag, 2)
       mag_bz = bz_v
       call add(mag_bz, external_bz_field, -1.)
    endif
    if(i3d.eq.1 .and. numvar.ge.2) then
       call associate_field(mag_bf, mag, 3)
       mag_bf = bf_field(1)
       call add(mag_bf, external_bf_field, -1.)
    endif

    thimprw = thimp
    
    fac = dt*eta_wall/delta_wall
    
    if(first_time) then

       call set_matrix_index(rw_rhs_mat, rw_rhs_mat_index)
       call create_mat(rw_rhs_mat, 2, 3, icomplex, .false.)
       call set_matrix_index(rw_lhs_mat, rw_lhs_mat_index)
       call create_mat(rw_lhs_mat, 2, 3, icomplex, .false.)
#ifdef USERW
       call setmatrixrwb(rw_rhs_mat_index, 1)
       call setmatrixrwb(rw_lhs_mat_index, 1)
#endif
#ifdef CJ_MATRIX_DUMP
       print *, "create_mat time_step rw_rhs_mat", rw_rhs_mat%imatrix     
       print *, "create_mat time_step rw_lhs_mat", rw_lhs_mat%imatrix     
#endif 

       ! calculate global dof ids
       allocate(global_dof_ids_1(nodes))
       allocate(global_dof_ids_row(nodes))
       allocate(global_dof_ids_col(nodes))
       allocate(buf(nodes))

       global_dof_ids_1 = 0
       global_dof_ids_row = 0
       global_dof_ids_col = 0

       do i=1, nodes
          if(local_id(i).le.0) cycle
          global_dof_ids_1(i) = global_node_index(rhs, global_id(i), 1)
          call get_global_node_indices(rw_rhs_mat, global_id(i), irow, icol)
          global_dof_ids_row(i) = irow(1,1)
          global_dof_ids_col(i) = icol(1,1)
       end do
       call mpi_allreduce(global_dof_ids_1, buf, nodes, MPI_INTEGER, &
            MPI_MAX, MPI_COMM_WORLD, ierr)
       global_dof_ids_1 = buf
       call mpi_allreduce(global_dof_ids_row, buf, nodes, MPI_INTEGER, &
            MPI_MAX, MPI_COMM_WORLD, ierr)
       global_dof_ids_row = buf
       call mpi_allreduce(global_dof_ids_col, buf, nodes, MPI_INTEGER, &
            MPI_MAX, MPI_COMM_WORLD, ierr)
       global_dof_ids_col = buf

       ! populate matrices
       ssi = 0.
       ddi = 0.
       ssj = 0.
       ddj = 0.
       rrj = 0.
       ss = 0.
       dd = 0.

       do i=1, nodes
          if(local_id(i).le.0) cycle

          call boundary_node(local_id(i), is_boundary, izone, izonedim, normi,&
               curvi, xii, zii)
          if(itor.eq.0) xii = 1.

          ibegin = global_dof_ids_1(i)

          if(ibegin.le.0) then
             print *, myrank, ': Error: ibegin <= 0', i, local_id(i), global_dof_ids_1(i)
             call safestop(0)
          endif

          do ii=1, 2
             do jj=1, dofs_per_node
                irow(ii,jj) = global_dof_ids_row(i) &
                     + (ii-1)*dofs_per_node + jj-1
             end do
          end do
          do ii=1, 3
             do jj=1, dofs_per_node
                icol(ii,jj) = global_dof_ids_col(i) &
                     + (ii-1)*dofs_per_node + jj-1
             end do
          end do


          do j=1, nodes
             normj(1) = nxnode(j)
             normj(2) = nznode(j)
             xjj = xnode(j)
             if(itor.eq.0) xjj = 1.
          
             jbegin = global_dof_ids_1(j)
             do ii=1, 2
                do jj=1, dofs_per_node
                   jrow(ii,jj) = global_dof_ids_row(j) &
                        + (ii-1)*dofs_per_node + jj-1
                end do
             end do
             do ii=1, 3
                do jj=1, dofs_per_node
                   jcol(ii,jj) = global_dof_ids_col(j) &
                        + (ii-1)*dofs_per_node + jj-1
                end do
             end do

             ! Four indices are:
             ! (dof, operator, equation, field) 

             ! dpsi/dt equation
             ! ~~~~~~~~~~~~~~~~
             ! coefficients of t_j.grad(psi_j)
             rrj(1,3,1,1) = -xii*zgrbth (i,j)/xjj
             rrj(3,3,1,1) = -xii*zgrbthp(i,j)/xjj
             if(itor.eq.1) then
                rrj(3,3,1,1) = rrj(3,3,1,1) + normi(2)*zgrbth(i,j)/xjj
             endif
             
             ! coefficients of n_j.grad(f'_j)
             rrj(1,2,1,3) = -xii*zgrbth (i,j)
             rrj(3,2,1,3) = -xii*zgrbthp(i,j)
             if(itor.eq.1) then
                rrj(3,2,1,3) = rrj(3,2,1,3) + normi(2)*zgrbth(i,j)
             endif

             rrj = rrj*fac
             ssj(:,:,:,1:2) =    -thimprw *rrj(:,:,:,1:2)
             ddj(:,:,:,1:2) = (1.-thimprw)*rrj(:,:,:,1:2)
             ddj(:,:,:,3)   =         rfac*rrj(:,:,:,3)

             if(i.eq.j) then              
                ! dpsi/dt
                ssi(1,1,1,1) = 1.
                ssi(3,3,1,1) = 1.
                ddi(1,1,1,1) = 1.
                ddi(3,3,1,1) = 1.

                ssi(1,2,1,1) =  fac*thimprw
                ssi(3,5,1,1) =  fac*thimprw
                ddi(1,2,1,1) = -fac*(1.-thimprw)
                ddi(3,5,1,1) = -fac*(1.-thimprw)
             endif

!!$             ! cylindrical approximation
!!$             ! ===========================
!!$             ssi = 0.
!!$             ssj = 0.
!!$             ddi = 0.
!!$             ddj = 0.
!!$
!!$             if(i.eq.j) then
!!$                ssi(1,1,1,1) = 1. - fac*(mp/a)**2/k*(-5.e-3)
!!$                
!!$                ssi(1,2,1,1) = fac
!!$                ddi(1,1,1,1) = 1.
!!$             endif
!!$             ! ===========================

             do ir=1, num_rows
                ii = rows(ir)
                ip = ibegin + ii - 1
                
                ! transform from boundary coordinates to global coordinates
                ss(ii,:,:,:) = ssj(ii,:,:,:)
                dd(ii,:,:,:) = ddj(ii,:,:,:)
                if(i.eq.j) then
                   ss(ii,:,:,:) = ss(ii,:,:,:) + ssi(ii,:,:,:)
                   dd(ii,:,:,:) = dd(ii,:,:,:) + ddi(ii,:,:,:)
                end if

                ! insert values into the appropriate matrices
                do jj=1, dofs_per_node
                   jp = jbegin + jj - 1
                   
!                   write(90+myrank,'(2I6,2f12.5)') i, j, ss(ii,jj,1,1)

                   if(present(mat)) then
                      call insert_global(mat,ss(ii,jj,1,1), &
                           ip+psi_off,jp+psi_off,MAT_SET)
                   endif
                   call insert_global(rw_lhs_mat,ss(ii,jj,1,1), &
                        irow(1,ii),jcol(1,jj),MAT_SET)
                   call insert_global(rw_rhs_mat,dd(ii,jj,1,1), &
                        irow(1,ii),jcol(1,jj),MAT_SET)
                   if(numvar.ge.2) then
                      call insert_global(rw_rhs_mat,dd(ii,jj,1,2), &
                        irow(1,ii),jcol(2,jj),MAT_SET)
                   endif
                   if(i3d.eq.1 .and. numvar.ge.2) then
                      call insert_global(rw_rhs_mat,dd(ii,jj,1,3), &
                        irow(1,ii),jcol(3,jj),MAT_SET)
                   end if
                end do
             end do
          end do
       end do

       if(use_resistive_bz) then
          call insert_resistive_wall_bz(rw_lhs_mat, rw_rhs_mat, mat)
       endif
       call finalize(rw_rhs_mat)
       call finalize(rw_lhs_mat)

!       call write_matrix(rw_rhs_mat, 'rw_rhs_mat')
!       call write_matrix(rw_lhs_mat, 'rw_lhs_mat')

       if(present(mat)) first_time = .false.

       deallocate(global_dof_ids_1, global_dof_ids_row, global_dof_ids_col)

       call create_vector(tempout0, 2)
       call matvecmult(rw_lhs_mat, external_field, tempout0)
    endif

    call create_vector(tempout, 2)
    call matvecmult(rw_rhs_mat, mag, tempout)

    ! add contributions to rhs vector
    do i=1, nodes
       if(local_id(i).le.0) cycle

       ibegin = node_index(rhs, local_id(i), psi_i)
       call get_node_data(tempout, 1, local_id(i), temp, rotate=.false.)
       call get_node_data(tempout0, 1, local_id(i), temp0, rotate=.false.)
       temp = temp + temp0
       do j=1, num_rows
          call insert(rhs, ibegin+rows(j)-1, temp(rows(j)), VEC_SET)
       end do

       if(use_resistive_bz) then
          ibegin = node_index(rhs, local_id(i), bz_i)
          call get_node_data(tempout, 2, local_id(i), temp, rotate=.false.)
          call get_node_data(tempout0, 2, local_id(i), temp0, rotate=.false.)
          temp = temp + temp0

          do j=1, dofs_per_node
             call insert(rhs, ibegin+j-1, temp(j), VEC_ADD)
          end do
       endif
    end do

    call destroy_vector(tempout)
    call destroy_vector(mag)
  end subroutine insert_resistive_wall

  subroutine insert_resistive_wall_bz(lhs_mat,rhs_mat,mat)
    use basic
    use arrays
    use sparse
    use vacuum_interface
    use matrix_mod
    use boundary_conditions
    use m3dc1_nint

    implicit none

    type(matrix_type), intent(inout) :: rhs_mat, lhs_mat
    type(matrix_type) :: mat

    integer :: i, j, k, itri, iedge, inode, jnode, idof, jdof
    vectype :: val, temp

    logical :: is_edge(3), is_boundary
    real :: norm(2,3)
    real :: curvi, xii, zii, xjj
    real, dimension(2) :: normi, normj
    integer :: idim(3), izone, izonedim
    integer :: ii, jj
    real :: thimprw

    integer, dimension(nodes_per_element) :: inodes
    integer :: iii

    integer :: ir_lh, ic_lh, ir_rh, ic_rh
    integer :: jr_lh, jc_lh, jr_rh, jc_rh

    thimprw = thimp

    ! cycle through each boundary edge
    do i=1, boundary_edges
       itri = edge_elm(i)
       iedge = edge_edge(i)
       call boundary_edge(itri, is_edge, norm, idim)
       
       call get_element_nodes(itri,inodes)

       call define_boundary_quadrature(itri, iedge, 5, 5, norm, idim)
       call define_fields(itri, 0, 1, linear)

       ! cycle through each node of the element containing the present edge
       do iii=1,nodes_per_element
          inode = global_node_id(inodes(iii))
          call boundary_node(inodes(iii), is_boundary, izone, izonedim, &
               normi, curvi, xii, zii)
          if(itor.eq.0) xjj = 1.

          ! determine the boundary-index of inode
          do j=1, nodes
             if(inodes(iii).eq.local_id(j)) exit
          end do
          if(j.gt.nodes) then
             if(is_boundary_node(inodes(iii))) then
                print *, 'error! boundary node not found'
                call safestop(9)
             endif
             cycle
          endif

          ! these are the global row and column id's of the first dof
          ! associated with boundary-node j for the two matrices
          ir_lh = global_dof_ids_1(j)
          ic_lh = global_dof_ids_1(j)
          ir_rh = global_dof_ids_row(j)
          ic_rh = global_dof_ids_col(j)

          ! cycle through each trial function associated with inode
          do ii=1,dofs_per_node

          ! idof is the element-based dof index
          ! of global node inode
          idof = (iii-1)*dofs_per_node + ii

          ! for each boundary edge, there are 2 associated nodes
          do j=1, 2
             ! we are only calculating the value (not derivatives) of
             ! the vacuum toroidal field; choose the basis function
             ! associated with the value at inode.

             ! jdof is the element-based dof index
             ! of the first dof associated with the current boundary node
             jdof = mod(iedge+j-2,nodes_per_element)*dofs_per_node+1
             
             ! plasma toroidal field term
             do jj=1, dofs_per_node
                temp = dt*eta_wall/delta_wall * &
                     int3(ri2_79,mu79(:,OP_1,idof),nu79(:,OP_1,jdof+jj-1))

                val = thimprw*temp
                call insert_global(mat,val, &
                     dof_index(ir_lh,bz_i,ii), &
                     dof_index(ic_lh,bz_i,jj),MAT_ADD)
                call insert_global(lhs_mat,val, &
                     dof_index(ir_rh,2,ii),dof_index(ic_rh,2,jj),MAT_ADD)
                
                val = (thimprw-1.)*temp
                call insert_global(rhs_mat,val, &
                     dof_index(ir_rh,2,ii),dof_index(ic_rh,2,jj),MAT_ADD)
             end do
            
             ! The coefficient of the chosen basis function is simply
             ! the value of the vacuum toroidal field at inode,
             ! which is a linear combination of the normal field
             ! at every boundary node.
             do k=1, nodes
                ! int[ mu V(x) dS]
                temp = dt*eta_wall/delta_wall * &
                     (zgrbph(edge_nodes(j,i),k) &
                     *int3(ri_79,mu79(:,OP_1,idof),nu79(:,OP_1,jdof)) )! &
!!$                     +zgrbphp(edge_nodes(j,i),k) &
!!$                     *(normi(1) &
!!$                     *int3(ri_79,mu79(:,OP_1,idof),nu79(:,OP_1,jdof+2)) &
!!$                     -normi(2) &
!!$                     *int3(ri_79,mu79(:,OP_1,idof),nu79(:,OP_1,jdof+1))))

                jnode = global_id(k)
                
                normj(1) = nxnode(k)
                normj(2) = nznode(k)
                xjj = xnode(k)
                if(itor.eq.0) xjj = 1.

                ! these are the global row and column id's of the first dof
                ! associated with boundary-node k for the two matrices
                jr_lh = global_dof_ids_1(k)
                jc_lh = global_dof_ids_1(k)
                jr_rh = global_dof_ids_row(k)
                jc_rh = global_dof_ids_col(k)
                
                ! -(1/R_j) * (dpsi/dt)
                val = thimprw*temp/xjj
                call insert_global(mat, val, &
                     dof_index(ir_lh,bz_i,ii),dof_index(jc_lh,psi_i,3),MAT_ADD)
                call insert_global(lhs_mat, val, &
                     dof_index(ir_rh,2,ii),dof_index(jc_rh,1,3),MAT_ADD)

                val = (thimprw-1.)*temp/xjj
                call insert_global(rhs_mat, val, &
                     dof_index(ir_rh,2,ii),dof_index(jc_rh,1,3),MAT_ADD)

                ! -df'/dn
                val = -temp*rfac
                call insert_global(rhs_mat, val, &
                     dof_index(ir_rh,2,ii),dof_index(jc_rh,3,2),MAT_ADD)
             end do
          end do ! on j
       end do    ! on ii
       end do    ! on iii
    end do
  end subroutine insert_resistive_wall_bz



!=============================
! scaleback
! ~~~~~~~~~
! rescale eigenfunction
!=============================
subroutine scaleback

  use basic
  use arrays
  use diagnostics

  implicit none

  vectype, parameter :: scalefac = 1.e-10

  if(ekin.lt.max_ke .or. max_ke.eq.0) return
  if(myrank.eq.0) write(*,*) " =>solution scaled back at time", time

  call mult(field_vec, scalefac)
  call mult(phiold_vec, scalefac)
  if(isplitstep.eq.1) then
     call mult(velold_vec, scalefac)
     call mult(veloldn_vec, scalefac)
     if(idens.eq.1) call mult(denold_vec, scalefac)
     if(ipres.eq.1) call mult(presold_vec, scalefac)
  endif
  if(i3d.eq.1) call mult(bf_field(1), scalefac)
  
end subroutine scaleback


subroutine variable_timestep

  use basic
  use arrays
  use diagnostics

  implicit none
  include 'mpif.h'
  integer :: ierr
!
! increase or decrease timestep based on kinetic energy and gamma_gr,
! but limit change to fraction dtfrac and bound by dtmin and demax
!
  dtold = dt
  if(dtkecrit.eq.0 .or. dtgamma.eq.0) return
  if(myrank.eq.0) then
!
!
    if(ekin.lt.dtkecrit) then
       if(gamma_gr.gt.0) then
          if(dt .gt.dtgamma/gamma_gr) then
            dt = dtold/(1. + dtfrac)
          endif
       else
            dt = dtold*(1. + dtfrac)
       endif
    else
            dt = dtold/(1. + dtfrac)
    endif
    dt = max(dt,dtmin)
    dt = min(dt,dtmax)
  endif
  if(iprint.ge.1) print *,"dtold,dt,gamma_gr,ekin",dtold,dt,gamma_gr,ekin
  call MPI_bcast(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  
end subroutine variable_timestep
 subroutine get_temperatures

  use basic
  use arrays
  use field
  use mesh_mod
  implicit none
  integer :: i, numnodes
   if(numvar.lt.3 .and. ipres.eq.0) return

   numnodes = owned_nodes()
   do i=1,numnodes

     if(idens.eq.1) then
        call get_node_data(den_v,i,den1_l)
        call get_node_data(den_field(0),i,den0_l)
     else
        if(eqsubtract.eq.1) then
           den1_l = 0.
           den0_l(1) = 1.
           den0_l(2:dofs_per_node) = 0.
        else
           den0_l = 0.
           den1_l(1) = 1.
           den1_l(2:dofs_per_node) = 0.
        endif
     endif

     if(numvar.eq.3) then
       if(ipres.eq.1) then
         call get_node_data(p_v,i,p1_l)
         call get_node_data(pe_v,i,pe1_l)
         call get_node_data(p_field(0),i,p0_l)
         call get_node_data(pe_field(0),i,pe0_l)
       else
         if(ipressplit.eq.0) then
            call get_node_data(pe_v,i,p1_l)
            call get_node_data(pe_field(0),i,p0_l)
         else
            call get_node_data(p_v,i,p1_l)
            call get_node_data(p_field(0),i,p0_l)
         endif
         pe1_l = pefac*p1_l
         pe0_l = pefac*p0_l
       endif
     else
       if(ipres.eq.1) then
         call get_node_data(p_v,i,p1_l)
         pe1_l = pefac*p1_l
         call get_node_data(p_field(0),i,p0_l)
         pe0_l = pefac*p0_l
       endif
     endif

     if(linear.eq.1 .or. eqsubtract.eq.1) then
        call calc_lin_electron_temperature(te1_l, pe0_l, den0_l, pe1_l, den1_l)
        call calc_lin_ion_temperature(ti1_l, p0_l, pe0_l, den0_l, p1_l, pe1_l, den1_l)
     else
        call calc_electron_temperature(te1_l,pe1_l,den1_l)
        call calc_ion_temperature(ti1_l, p1_l, pe1_l, den1_l)
     endif

     call set_node_data(te_v,i,te1_l)
     if(ipres.eq.1) then
        if(ipressplit.eq.1) then
           call set_node_data(ti_v,i,ti1_l)
        else
           call set_node_data(ti_field(1),i,ti1_l)
        endif
     endif

   enddo


   return
 end subroutine get_temperatures
 subroutine get_pressures

  use basic
  use arrays
  use field
  use mesh_mod
  implicit none
  integer :: i, numnodes
   if(numvar.lt.3 .and. ipres.eq.0) return

   numnodes = owned_nodes()
   do i=1,numnodes

     if(idens.eq.1) then
        call get_node_data(den_v,i,den1_l)
        call get_node_data(den_field(0),i,den0_l)
     else
        if(eqsubtract.eq.1) then
           den1_l = 0.
           den0_l(1) = 1.
           den0_l(2:dofs_per_node) = 0.
        else
           den0_l = 0.
           den1_l(1) = 1.
           den1_l(2:dofs_per_node) = 0.
        endif
     endif

     call get_node_data(te_v,i,te1_l)
     call get_node_data(te_field(0),i,te0_l)
     call get_node_data(ti_field(0),i,ti0_l)
     if(ipressplit.eq.1 .and. ipres.eq.1) then
        call get_node_data(ti_v,i,ti1_l)
     else
        ti1_l = te1_l*(1.-pefac)/pefac
     endif

     if(linear.eq.1 .or. eqsubtract.eq.1) then
        call calc_lin_electron_pressure(pe1_l, te0_l, den0_l, te1_l, den1_l)
        call calc_lin_pressure(p1_l, te0_l, ti0_l, den0_l, te1_l, ti1_l, den1_l)
     else
        call calc_tot_electron_pressure(pe1_l,te1_l,den1_l)
        call calc_tot_pressure(p1_l, te1_l, ti1_l, den1_l)
     endif

     call set_node_data(p_v,i,p1_l)
     if(ipres.ge.1) then
        call set_node_data(pe_v,i,pe1_l)
     endif


   enddo


   return
 end subroutine get_pressures
end module time_step
! calc_electron_temperature
! ~~~~~~~~~~~~~~~~~~~~~~
!
! calculates the electron temperature
!======================================================================
subroutine calc_electron_temperature(te, pe0, n0)
  use basic

  implicit none

  vectype, intent(out), dimension(dofs_per_node) :: te
  vectype, intent(in), dimension(dofs_per_node) :: pe0, n0

     te(1) = pe0(1)/n0(1)
     te(2) = pe0(2)/n0(1) - pe0(1)*n0(2)/n0(1)**2
     te(3) = pe0(3)/n0(1) - pe0(1)*n0(3)/n0(1)**2
     te(4) = pe0(4)/n0(1) - 2.*pe0(2)*n0(2)/n0(1)**2                       &
           - pe0(1)*n0(4)/n0(1)**2 + 2*pe0(1)*n0(2)*n0(2)/n0(1)**3
     te(5) = pe0(5)/n0(1) - pe0(3)*n0(2)/n0(1)**2  - pe0(2)*n0(3)/n0(1)**2  &
           - pe0(1)*n0(5)/n0(1)**2 + 2*pe0(1)*n0(3)*n0(2)/n0(1)**3
     te(6) = pe0(6)/n0(1) - 2.*pe0(3)*n0(3)/n0(1)**2                       &
           - pe0(1)*n0(6)/n0(1)**2 + 2*pe0(1)*n0(3)*n0(3)/n0(1)**3
     return

end subroutine calc_electron_temperature
! calc_ion_temperature
! ~~~~~~~~~~~~~~~~~~~~~~
!
! calculates the ion temperature
!======================================================================
subroutine calc_ion_temperature(ti, pres0, pe0, n0)
  use basic

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: pres0, pe0, n0
  vectype, intent(out), dimension(dofs_per_node) :: ti
  vectype, dimension(dofs_per_node) :: pion0

     pion0 = pres0 - pe0

     ti(1) = pion0(1)/n0(1)
     ti(2) = pion0(2)/n0(1) - pion0(1)*n0(2)/n0(1)**2
     ti(3) = pion0(3)/n0(1) - pion0(1)*n0(3)/n0(1)**2
     ti(4) = pion0(4)/n0(1) - 2.*pion0(2)*n0(2)/n0(1)**2                       &
           - pion0(1)*n0(4)/n0(1)**2 + 2*pion0(1)*n0(2)*n0(2)/n0(1)**3
     ti(5) = pion0(5)/n0(1) - pion0(3)*n0(2)/n0(1)**2  - pion0(2)*n0(3)/n0(1)**2  &
           - pion0(1)*n0(5)/n0(1)**2 + 2*pion0(1)*n0(3)*n0(2)/n0(1)**3
     ti(6) = pion0(6)/n0(1) - 2.*pion0(3)*n0(3)/n0(1)**2                       &
           - pion0(1)*n0(6)/n0(1)**2 + 2*pion0(1)*n0(3)*n0(3)/n0(1)**3
     return

end subroutine calc_ion_temperature
! calc_lin_electron_temperature
! ~~~~~~~~~~~~~~~~~~~~~~
!
! calculates the linearized electron temperature
!======================================================================
subroutine calc_lin_electron_temperature(te, pe0, n0, pe1, n1)
  use basic

  implicit none

  vectype, intent(out), dimension(dofs_per_node) :: te
  vectype, intent(in), dimension(dofs_per_node) :: pe0, n0, pe1, n1

     te(1) = (pe0(1)+pe1(1))/(n0(1)+n1(1))  &
             -pe0(1)/n0(1)
     te(2) = (pe0(2)+pe1(2))/(n0(1)+n1(1)) &
           - (pe0(1)+pe1(1))*(n0(2)+n1(2))/(n0(1)+n1(1))**2  &
             -pe0(2)/n0(1) + pe0(1)*n0(2)/n0(1)**2
     te(3) = (pe0(3)+pe1(3))/(n0(1)+n1(1)) &
           - (pe0(1)+pe1(1))*(n0(3)+n1(3))/(n0(1)+n1(1))**2  &
             -pe0(3)/n0(1) + pe0(1)*n0(3)/n0(1)**2
     te(4) = (pe0(4)+pe1(4))/(n0(1)+n1(1)) &
        - 2.*(pe0(2)+pe1(2))*(n0(2)+n1(2))/(n0(1)+n1(1))**2    &
           - (pe0(1)+pe1(1))*(n0(4)+n1(4))/(n0(1)+n1(1))**2  &
           + 2*(pe0(1)+pe1(1))*(n0(2)+n1(2))*(n0(2)+n1(2))/(n0(1)+n1(1))**3   &
            - pe0(4)/n0(1) + 2.*pe0(2)*n0(2)/n0(1)**2                         &
           +  pe0(1)*n0(4)/n0(1)**2 - 2*pe0(1)*n0(2)*n0(2)/n0(1)**3
     te(5) = (pe0(5)+pe1(5))/(n0(1)+n1(1)) &
           - (pe0(3)+pe1(3))*(n0(2)+n1(2))/(n0(1)+n1(1))**2  &
           - (pe0(2)+pe1(2))*(n0(3)+n1(3))/(n0(1)+n1(1))**2  &
           - (pe0(1)+pe1(1))*(n0(5)+n1(5))/(n0(1)+n1(1))**2 &
           + 2*(pe0(1)+pe1(1))*(n0(3)+n1(3))*(n0(2)+n1(2))/(n0(1)+n1(1))**3  &
            - pe0(5)/n0(1) + pe0(3)*n0(2)/n0(1)**2  + pe0(2)*n0(3)/n0(1)**2  &
           +  pe0(1)*n0(5)/n0(1)**2 - 2*pe0(1)*n0(3)*n0(2)/n0(1)**3
     te(6) = (pe0(6)+pe1(6))/(n0(1)+n1(1)) &
        - 2.*(pe0(3)+pe1(3))*(n0(3)+n1(3))/(n0(1)+n1(1))**2  &
           - (pe0(1)+pe1(1))*(n0(6)+n1(6))/(n0(1)+n1(1))**2   &
           + 2*(pe0(1)+pe1(1))*(n0(3)+n1(3))*(n0(3)+n1(3))/(n0(1)+n1(1))**3   &
           -  pe0(6)/n0(1) + 2.*pe0(3)*n0(3)/n0(1)**2                         &
           +  pe0(1)*n0(6)/n0(1)**2 - 2*pe0(1)*n0(3)*n0(3)/n0(1)**3
     return

end subroutine calc_lin_electron_temperature
! calc_lin_ion_temperature
! ~~~~~~~~~~~~~~~~~~~~~~
!
! calculates the linearized ion temperature
!======================================================================
subroutine calc_lin_ion_temperature(ti, pres0, pe0, n0, pres1, pe1, n1)
  use basic

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: pres0, pe0, n0, pres1, pe1, n1
  vectype, intent(out), dimension(dofs_per_node) :: ti
  vectype, dimension(dofs_per_node) :: pion0, pion1

     pion0 = pres0 - pe0
     pion1 = pres1 - pe1

     ti(1) = (pion0(1)+pion1(1))/(n0(1)+n1(1))  &
           -pion0(1)/n0(1)
     ti(2) = (pion0(2)+pion1(2))/(n0(1)+n1(1)) &
           - (pion0(1)+pion1(1))*(n0(2)+n1(2))/(n0(1)+n1(1))**2   &
            -pion0(2)/n0(1) + pion0(1)*n0(2)/n0(1)**2
     ti(3) = (pion0(3)+pion1(3))/(n0(1)+n1(1)) &
           - (pion0(1)+pion1(1))*(n0(3)+n1(3))/(n0(1)+n1(1))**2   &
             -pion0(3)/n0(1) + pion0(1)*n0(3)/n0(1)**2
     ti(4) = (pion0(4)+pion1(4))/(n0(1)+n1(1)) &
           - 2.*(pion0(2)+pion1(2))*(n0(2)+n1(2))/(n0(1)+n1(1))**2   &
           - (pion0(1)+pion1(1))*(n0(4)+n1(4))/(n0(1)+n1(1))**2 &
           + 2*(pion0(1)+pion1(1))*(n0(2)+n1(2))*(n0(2)+n1(2))/(n0(1)+n1(1))**3   &
           - pion0(4)/n0(1) + 2.*pion0(2)*n0(2)/n0(1)**2       &              
           + pion0(1)*n0(4)/n0(1)**2 - 2*pion0(1)*n0(2)*n0(2)/n0(1)**3

     ti(5) = (pion0(5)+pion1(5))/(n0(1)+n1(1)) &
           - (pion0(3)+pion1(3))*(n0(2)+n1(2))/(n0(1)+n1(1))**2  &
           - (pion0(2)+pion1(2))*(n0(3)+n1(3))/(n0(1)+n1(1))**2  &
           - (pion0(1)+pion1(1))*(n0(5)+n1(5))/(n0(1)+n1(1))**2  &
         + 2*(pion0(1)+pion1(1))*(n0(3)+n1(3))*(n0(2)+n1(2))/(n0(1)+n1(1))**3   &
            - pion0(5)/n0(1) + pion0(3)*n0(2)/n0(1)**2  &
            + pion0(2)*n0(3)/n0(1)**2  &
            + pion0(1)*n0(5)/n0(1)**2 - 2*pion0(1)*n0(3)*n0(2)/n0(1)**3

     ti(6) = (pion0(6)+pion1(6))/(n0(1)+n1(1)) &
           - 2.*(pion0(3)+pion1(3))*(n0(3)+n1(3))/(n0(1)+n1(1))**2   &
           - (pion0(1)+pion1(1))*(n0(6)+n1(6))/(n0(1)+n1(1))**2 &
           + 2*(pion0(1)+pion1(1))*(n0(3)+n1(3))*(n0(3)+n1(3))/(n0(1)+n1(1))**3  &
           -  pion0(6)/n0(1) + 2.*pion0(3)*n0(3)/n0(1)**2                       &
           +  pion0(1)*n0(6)/n0(1)**2 - 2*pion0(1)*n0(3)*n0(3)/n0(1)**3
    
     return

end subroutine calc_lin_ion_temperature
! ~~~~~~~~~~~~~~~~~~~~~~
!
!======================================================================
subroutine calc_tot_electron_pressure(pe, te0, n0)
  use basic

  implicit none

  vectype, intent(out), dimension(dofs_per_node) :: pe
  vectype, intent(in), dimension(dofs_per_node) :: te0, n0

     pe(1) = te0(1)*n0(1)
     pe(2) = te0(2)*n0(1) + te0(1)*n0(2)
     pe(3) = te0(3)*n0(1) + te0(1)*n0(3)
     pe(4) = te0(4)*n0(1) + 2.*te0(2)*n0(2) + te0(1)*n0(4) 
     pe(5) = te0(5)*n0(1) + te0(3)*n0(2)  &
           + te0(2)*n0(3) + te0(1)*n0(5)
     pe(6) = te0(6)*n0(1) + 2.*te0(3)*n0(3) + te0(1)*n0(6)

     return

end subroutine calc_tot_electron_pressure
! calc_tot_pressure
! ~~~~~~~~~~~~~~~~~~~~~~
!
! calculates the total pressure
!======================================================================
subroutine calc_tot_pressure(pres0, te0, ti0, n0)
  use basic

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: te0, ti0, n0
  vectype, intent(out), dimension(dofs_per_node) :: pres0
  vectype, dimension(dofs_per_node) :: tepti

     tepti = te0 + ti0

     pres0(1) = tepti(1)*n0(1)
     pres0(2) = tepti(2)*n0(1) + tepti(1)*n0(2)
     pres0(3) = tepti(3)*n0(1) + tepti(1)*n0(3)
     pres0(4) = tepti(4)*n0(1) + 2.*tepti(2)*n0(2) + tepti(1)*n0(4) 
     pres0(5) = tepti(5)*n0(1) + tepti(3)*n0(2)  &
              + tepti(2)*n0(3) + tepti(1)*n0(5)
     pres0(6) = tepti(6)*n0(1) + 2.*tepti(3)*n0(3) + tepti(1)*n0(6)

     return

end subroutine calc_tot_pressure
subroutine calc_lin_electron_pressure(pe, te0, n0, te1, n1)
  use basic

  implicit none

  vectype, intent(out), dimension(dofs_per_node) :: pe
  vectype, intent(in), dimension(dofs_per_node) :: te0, n0, te1, n1

     pe(1) = (te0(1)+te1(1))*(n0(1)+n1(1)) &
           -  te0(1)*n0(1)

     pe(2) = (te0(2)+te1(2))*(n0(1)+n1(1)) &
           + (te0(1)+te1(1))*(n0(2)+n1(2)) &
           -  te0(2)*n0(1) - te0(1)*n0(2)
 

     pe(3) = (te0(3)+te1(3))*(n0(1)+n1(1)) &
           + (te0(1)+te1(1))*(n0(3)+n1(3)) &
           -  te0(3)*n0(1) - te0(1)*n0(3)

     pe(4) = (te0(4)+te1(4))*(n0(1)+n1(1)) &
           + 2.*(te0(2)+te1(2))*(n0(2)+n1(2)) &
           + (te0(1)+te1(1))*(n0(4)+n1(4))    &
           - te0(4)*n0(1) - 2.*te0(2)*n0(2) - te0(1)*n0(4)

     pe(5) = (te0(5)+te1(5))*(n0(1)+n1(1)) &
           + (te0(3)+te1(3))*(n0(2)+n1(2)) &
           + (te0(2)+te1(2))*(n0(3)+n1(3)) &
           + (te0(1)+te1(1))*(n0(5)+n1(5)) &
           - te0(5)*n0(1) - te0(3)*n0(2) -  te0(2)*n0(3) - te0(1)*n0(5)

     pe(6) = (te0(6)+te1(6))*(n0(1)+n1(1)) &
           + 2.*(te0(3)+te1(3))*(n0(3)+n1(3)) &
           + (te0(1)+te1(1))*(n0(6)+n1(6)) &
           - te0(6)*n0(1) - 2.*te0(3)*n0(3) - te0(1)*n0(6)

     return

end subroutine calc_lin_electron_pressure
! calc_linearized pressure
! ~~~~~~~~~~~~~~~~~~~~~~
!
! calculates the linearized pressure
!======================================================================
subroutine calc_lin_pressure(pres1, te0, ti0, n0, te1, ti1, n1)
  use basic

  implicit none

  vectype, intent(in), dimension(dofs_per_node)  :: te0, ti0, n0, te1, ti1, n1
  vectype, intent(out), dimension(dofs_per_node) :: pres1
  vectype, dimension(dofs_per_node) :: tepti, tepti0, nt

     tepti0 = te0 + ti0
     tepti  = te0 + te1 + ti0 + ti1
     nt = n0 + n1

     pres1(1) = tepti(1)*nt(1) - tepti0(1)*n0(1)
     pres1(2) = tepti(2)*nt(1) + tepti(1)*nt(2) &
             -  tepti0(2)*n0(1)- tepti0(1)*n0(2)
     pres1(3) = tepti(3)*nt(1) + tepti(1)*nt(3) &
              - tepti0(3)*n0(1)- tepti0(1)*n0(3)
     pres1(4) = tepti(4)*nt(1) + 2.*tepti(2)*nt(2) + tepti(1)*nt(4) &
              - tepti0(4)*n0(1)- 2.*tepti0(2)*n0(2)- tepti0(1)*n0(4)
     pres1(5) = tepti(5)*nt(1) + tepti(3)*nt(2)  &
              + tepti(2)*nt(3) + tepti(1)*nt(5)  &
              - tepti0(5)*n0(1)- tepti0(3)*n0(2) &
              - tepti0(2)*n0(3)- tepti0(1)*n0(5)
     pres1(6) = tepti(6)*nt(1) + 2.*tepti(3)*nt(3) + tepti(1)*nt(6) &
              - tepti0(6)*n0(1)- 2.*tepti0(3)*n0(3)- tepti0(1)*n0(6)




     return

end subroutine calc_lin_pressure
