module time_step_split

  use field
  use matrix_mod
  use model

  type(vector_type), private :: phi_vec, phip_vec
  type(vector_type), private :: vel_vec, veln_vec
  type(vector_type), private :: pres_vec
  type(vector_type), private :: den_vec, denold_vec, ne_vec, neold_vec
  type(vector_type), private :: pret_vec
  type(vector_type), private :: nre_vec, nreold_vec
  type(vector_type), private :: densmooth_vec

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
  type(field_type), private ::  bfp_v
  type(field_type), private ::  te_v
  type(field_type), private ::  ti_v
  type(field_type), private ::  ne_v
  type(field_type), private :: nre_v
  type(field_type), private :: densmooth_v

  integer :: vecsize_vel, vecsize_phi, vecsize_n, vecsize_p, vecsize_t

  ! temporary vectors
  type(vector_type), private :: b1_vel, b2_vel, b1_phi, b2_phi

  logical, private :: initialized

contains

  subroutine initialize_timestep_split
    use sparse
    use basic
    use m3dc1_vel_prof
    implicit none
    integer :: jer

    vecsize_phi = numvar
    vecsize_vel = numvar
    vecsize_n = 1
    vecsize_p = 1
    imode = 1
    if(itemp.eq.1) imode = 2
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

    ! add bf equation
    if(imp_bf.eq.1) vecsize_phi = vecsize_phi + 1

    ! add electrostatic potential equation OR implicit hyper-resistivity
    if((jadv.eq.0) .or. (jadv.eq.1 .and. imp_hyper.ge.1)) &
                    vecsize_phi = vecsize_phi + 1
  
if (ispradapt .eq. 1) then
    call create_vector(phi_vec,      vecsize_phi, "phi_vec")
    call create_vector(phip_vec,     vecsize_phi, "phip_vec")
    call create_vector(q4_vec,       vecsize_phi, "q4_vec")
    call mark_vector_for_solutiontransfer(phi_vec)
    call mark_vector_for_solutiontransfer(phip_vec)
    call mark_vector_for_solutiontransfer(q4_vec)

    call create_vector(b1_phi, vecsize_phi, "b1_phi")
    call create_vector(b2_phi, vecsize_phi, "b2_phi")
    call mark_vector_for_solutiontransfer(b1_phi)
    call mark_vector_for_solutiontransfer(b2_phi)

    call create_vector(vel_vec,      vecsize_vel, "vel_vec")
    call create_vector(veln_vec,     vecsize_vel, "veln_vec")
    call create_vector(r4_vec,       vecsize_vel, "r4_vec")
    call mark_vector_for_solutiontransfer(vel_vec)
    call mark_vector_for_solutiontransfer(veln_vec)
    call mark_vector_for_solutiontransfer(r4_vec)

    if(ipres.eq.1 .or. ipressplit.eq.1) then
       call create_vector(pres_vec,    vecsize_p, "pres_vec")
       call create_vector(qp4_vec,     vecsize_p, "qp4_vec")
       call mark_vector_for_solutiontransfer(pres_vec)
       call mark_vector_for_solutiontransfer(qp4_vec)
    endif

    call create_vector(den_vec,    vecsize_n, "den_vec")
    call create_vector(denold_vec, vecsize_n, "denold_vec")
    call create_vector(ne_vec,     1, "ne_vec")
    call create_vector(neold_vec,  1, "neold_vec")
    call create_vector(qn4_vec,    vecsize_n, "qn4_vec")
    call create_vector(densmooth_vec,    vecsize_n, "densmooth_vec")
    call mark_vector_for_solutiontransfer(den_vec)
    call mark_vector_for_solutiontransfer(denold_vec)
    call mark_vector_for_solutiontransfer(ne_vec)
    call mark_vector_for_solutiontransfer(neold_vec)
    call mark_vector_for_solutiontransfer(qn4_vec)
    call mark_vector_for_solutiontransfer(densmooth_vec)

    if(irunaway .gt. 0) then
      call create_vector(nre_vec,    vecsize_n, "nre_vec")
      call create_vector(nreold_vec, vecsize_n, "nreold_vec")
      call create_vector(qn5_vec,    vecsize_n, "qn5_vec")
      call mark_vector_for_solutiontransfer(nre_vec)
      call mark_vector_for_solutiontransfer(nreold_vec)
      call mark_vector_for_solutiontransfer(qn5_vec)
    endif

    call create_vector(b1_vel, vecsize_vel, "b1_vel")
    call create_vector(b2_vel, vecsize_vel, "b2_vel")
    call mark_vector_for_solutiontransfer(b1_vel)
    call mark_vector_for_solutiontransfer(b2_vel)

    call create_vector(pret_vec,    vecsize_t, "pret_vec")
    call mark_vector_for_solutiontransfer(pret_vec)
else

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
    call create_vector(ne_vec,     1)
    call create_vector(neold_vec,  1)
    call create_vector(qn4_vec,    vecsize_n)
    call create_vector(densmooth_vec,    vecsize_n)

    if(irunaway .gt. 0) then
      call create_vector(nre_vec,    vecsize_n)
      call create_vector(nreold_vec, vecsize_n)
      call create_vector(qn5_vec,    vecsize_n)
    endif

    call create_vector(b1_vel, vecsize_vel)
    call create_vector(b2_vel, vecsize_vel)

    call create_vector(pret_vec,    vecsize_t)
endif

call PetscLogStagePush(stageA,jer)
    ! Matrices associated with velocity advance
    call set_matrix_index(s1_mat, s1_mat_index)
    call set_matrix_index(d1_mat, d1_mat_index)
    call set_matrix_index(q1_mat, q1_mat_index)
    call set_matrix_index(r14_mat, r14_mat_index)
#ifdef REORDERED       
    call create_mat(s1_mat, vecsize_vel, vecsize_vel, icomplex, 1, is1_agg_blk_cnt, is1_agg_scp)
#else
    call create_mat(s1_mat, vecsize_vel, vecsize_vel, icomplex, 1)
#endif
    call create_mat(d1_mat, vecsize_vel, vecsize_vel, icomplex, 0)
    call create_mat(q1_mat, vecsize_vel, vecsize_phi, icomplex, 0)
    call create_mat(r14_mat, vecsize_vel, vecsize_n, icomplex, 0)

    if(i3d.eq.1) then
       call set_matrix_index(o1_mat, o1_mat_index)
       call create_mat(o1_mat, vecsize_vel, 1, icomplex, 0)
    endif

    if((ipres.eq.1 .and. numvar.lt.3) .or. ipressplit.gt.0) then
       call set_matrix_index(p1_mat, p1_mat_index)
       call create_mat(p1_mat, vecsize_vel, vecsize_p, icomplex, 0)
    end if
call PetscLogStagePop(jer)

    ! Matrices associated with magnetic field advance
    call set_matrix_index(s2_mat, s2_mat_index)
    call set_matrix_index(d2_mat, d2_mat_index)
    call set_matrix_index(r2_mat, r2_mat_index)
    call set_matrix_index(q2_mat, q2_mat_index)
    call create_mat(s2_mat, vecsize_phi, vecsize_phi, icomplex, 1)
    call create_mat(d2_mat, vecsize_phi, vecsize_phi, icomplex, 0)
    call create_mat(r2_mat, vecsize_phi, vecsize_vel, icomplex, 0)
    call create_mat(q2_mat, vecsize_phi, vecsize_vel, icomplex, 0)

    if(idens.eq.1) then
       call set_matrix_index(r42_mat, r42_mat_index)
       call set_matrix_index(q42_mat, q42_mat_index)
       call create_mat(r42_mat, vecsize_phi, 1, icomplex, 0)
       call create_mat(q42_mat, vecsize_phi, 1, icomplex, 0)
    endif
    
    if(irunaway .gt. 0) then
       call set_matrix_index(r43_mat, r43_mat_index)
       call set_matrix_index(q43_mat, q43_mat_index)
       call create_mat(r43_mat, vecsize_phi, 1, icomplex, 0)
       call create_mat(q43_mat, vecsize_phi, 1, icomplex, 0)
    endif

    if(i3d.eq.1) then
       call set_matrix_index(o2_mat, o2_mat_index)
       call create_mat(o2_mat, vecsize_phi, 1, icomplex, 0)
       if(ipressplit.eq.1 .or. ipres.eq.1) then
          call set_matrix_index(o3_mat, o3_mat_index)
          call create_mat(o3_mat, vecsize_p, 1, icomplex, 0)
       endif
    endif

    ! Matrices associated with density advance
    if(idens.eq.1) then
       call set_matrix_index(s8_mat, s8_mat_index)
       call set_matrix_index(d8_mat, d8_mat_index)
       call set_matrix_index(r8_mat, r8_mat_index)
       call set_matrix_index(q8_mat, q8_mat_index)
       call create_mat(s8_mat, vecsize_n, vecsize_n, icomplex, 1)
       call create_mat(d8_mat, vecsize_n, vecsize_n, icomplex, 0)
       call create_mat(r8_mat, vecsize_n, vecsize_vel, icomplex, 0)
       call create_mat(q8_mat, vecsize_n, vecsize_vel, icomplex, 0)
    endif

   ! Matrices associated with RE density advance
   if(irunaway .gt. 0) then
       call set_matrix_index(s15_mat, s15_mat_index)
       call set_matrix_index(d15_mat, d15_mat_index)
       call set_matrix_index(r15_mat, r15_mat_index)
       call set_matrix_index(q15_mat, q15_mat_index)
       call set_matrix_index(k15_mat, k15_mat_index)
       call create_mat(s15_mat, vecsize_n, vecsize_n, icomplex, 1)
       call create_mat(d15_mat, vecsize_n, vecsize_n, icomplex, 0)
       call create_mat(r15_mat, vecsize_n, vecsize_phi, icomplex, 0)
       call create_mat(q15_mat, vecsize_n, vecsize_phi, icomplex, 0)
       call create_mat(k15_mat, vecsize_n, vecsize_vel, icomplex, 0)
    endif

    ! Matrices associated with separate pressure advance
    if(ipres.eq.1 .or. ipressplit.eq.1) then
       call set_matrix_index(s9_mat, s9_mat_index)
       call set_matrix_index(d9_mat, d9_mat_index)
       call set_matrix_index(r9_mat, r9_mat_index)
       call set_matrix_index(q9_mat, q9_mat_index)
       call set_matrix_index(o9_mat, o9_mat_index)
       call create_mat(s9_mat, vecsize_p, vecsize_p, icomplex, 1)
       call create_mat(d9_mat, vecsize_p, vecsize_p, icomplex, 0)
       call create_mat(r9_mat, vecsize_p, vecsize_vel, icomplex, 0)
       call create_mat(q9_mat, vecsize_p, vecsize_vel, icomplex, 0)
       call create_mat(o9_mat, vecsize_p, vecsize_phi, icomplex, 0)
       if(idens.eq.1) then
          call set_matrix_index(rp42_mat, rp42_mat_index)
          call set_matrix_index(qp42_mat, qp42_mat_index)
          call create_mat(rp42_mat, vecsize_p, 1, icomplex, 0)
          call create_mat(qp42_mat, vecsize_p, 1, icomplex, 0)
       endif
    endif
       
!      if(ipressplit.eq.1) then
!         call set_matrix_index(s11_mat, s11_mat_index)
!         call set_matrix_index(d11_mat, d11_mat_index)
!         call create_mat(s11_mat, vecsize_n, vecsize_n, icomplex, .true.)
!         call create_mat(d11_mat, vecsize_n, vecsize_n, icomplex, .false.)
!         call set_matrix_index(s12_mat, s12_mat_index)
!         call set_matrix_index(d12_mat, d12_mat_index)
!         call create_mat(s12_mat, vecsize_n, vecsize_n, icomplex, .true.)
!         call create_mat(d12_mat, vecsize_n, vecsize_n, icomplex, .false.)
!      endif

       initialized = .true.
  end subroutine initialize_timestep_split

  subroutine finalize_timestep_split
    use basic
    use m3dc1_vel_prof
    implicit none
    integer :: jer

    if(.not.initialized) return
      
call PetscLogStagePush(stageA,jer)
    call destroy_mat(s1_mat)
    call destroy_mat(d1_mat)
    call destroy_mat(q1_mat)
    call destroy_mat(r14_mat)
    if(i3d.eq.1) call destroy_mat(o1_mat)
call PetscLogStagePop(jer)

    call destroy_mat(s2_mat)
    call destroy_mat(d2_mat)
    call destroy_mat(r2_mat)
    call destroy_mat(q2_mat)
    if(idens.eq.1) then
       call destroy_mat(r42_mat)
       call destroy_mat(q42_mat)
    end if
    if(irunaway .gt. 0) then
       call destroy_mat(r43_mat)
       call destroy_mat(q43_mat)
    end if
    if(i3d.eq.1) then
      call destroy_mat(o2_mat)
      if(ipres.eq.1 .or. ipressplit.gt.0) call destroy_mat(o3_mat)
    endif
    if((ipres.eq.1 .and. numvar.lt.3) .or. ipressplit.gt.0) then        
       call destroy_mat(p1_mat)
    end if
       
    if(idens.eq.1) then
       call destroy_mat(s8_mat)
       call destroy_mat(d8_mat)
       call destroy_mat(r8_mat)
       call destroy_mat(q8_mat)
    endif

    if(irunaway .gt. 0) then
       call destroy_mat(s15_mat)
       call destroy_mat(d15_mat)
       call destroy_mat(r15_mat)
       call destroy_mat(q15_mat)
       call destroy_mat(k15_mat)
    endif

    if(ipres.eq.1 .or. ipressplit.eq.1) then
       call destroy_mat(s9_mat)
       call destroy_mat(d9_mat)
       call destroy_mat(r9_mat)
       call destroy_mat(q9_mat)
       call destroy_mat(o9_mat)
       if(idens.eq.1) then
          call destroy_mat(rp42_mat)
          call destroy_mat(qp42_mat)
       end if
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
    if(irunaway.gt.0) nre_i = 1
    if(imp_bf.eq.1) then
       if((jadv.eq.0).or.(jadv.eq.1 .and. imp_hyper.ge.1)) then
         bf_i = vecsize_phi - 1
         e_i =  vecsize_phi
       else
         bf_i = vecsize_phi
         e_i = vecsize_phi + 1
       endif
    else
       if((jadv.eq.0).or.(jadv.eq.1 .and. imp_hyper.ge.1)) then
         ! bf_i = vecsize_phi + 1
         bf_i = 1     ! changed 11/25/20  SCJ
         e_i = vecsize_phi
       else
         bf_i = 1
         e_i = vecsize_phi + 1
       endif
    end if
    if(ipressplit.eq.1) then
       p_i  = 1
       te_i = 1
       if(ipres.eq.1 .and. numvar.eq.3) then
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
       if(ipres.eq.1 .and. numvar.eq.3) then
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
    
    call associate_field(den_v,  den_vec,    den_i)
    call associate_field(ne_v,    ne_vec,        1)
    if(irunaway.gt.0) call associate_field(nre_v,  nre_vec,    nre_i)
#ifdef USEPARTICLES
    if(itwofluid.gt.0) call associate_field(densmooth_v,  densmooth_vec,    den_i)
#endif
    
    if((jadv.eq.0).or.(jadv.eq.1 .and. imp_hyper.ge.1)) then
       call associate_field(e_v, phi_vec, e_i)
    end if
    
    if(imp_bf.eq.1) then
       call associate_field(bfp_v, phi_vec, bf_i)
    end if
  end subroutine assign_variables_split

  subroutine clear_matrices_split
    use basic
    use m3dc1_vel_prof
    implicit none
    integer :: jer

call PetscLogStagePush(stageA,jer)
#ifdef VERSION1
    if(mod(ntime,pskip)==0) then
        if(myrank.eq.0) print *, " clear_mat s1_mat",ntime, s1_mat%imatrix
    call clear_mat(s1_mat)
    else
        if(myrank.eq.0) print *, " zero_mat  s1_mat",ntime, s1_mat%imatrix
        call zero_mat(s1_mat)
    endif
#else
    if(pskip==0) then
            !orgin solve option: kspSet=1, with KSP Object completely destroyed after each solve
            if(myrank.eq.0) print *, " clear_mat s1_mat",ntime, s1_mat%imatrix
            call clear_mat(s1_mat)
    else if(pskip==1) then
            !default
            !third solve option: kspSet=3, similar to pskip>1, 
            !update both A&P contents at every timestep, but no KSP Object destroyed after each solve
            !to save memory usage and solve setup time.
            if(myrank.eq.0) print *, " update_mat s1_mat",ntime, s1_mat%imatrix
            call update_mat(s1_mat)

    else !pskip>1
            !second solve option: kspSet=2, 
            !update A content only at every timestep, update P every 'pskip' number of timesteps. 
            if(mod(ntime,pskip)==0) then
                    if(myrank.eq.0) print *, " update_mat s1_mat",ntime, s1_mat%imatrix 
                    call update_mat(s1_mat)
            else
                    if(myrank.eq.0) print *, " zero_mat s1_mat",ntime, s1_mat%imatrix 
                    call zero_mat(s1_mat)
            end if
    endif
#endif
    call clear_mat(d1_mat)
    call clear_mat(q1_mat)
    call clear_mat(r14_mat)
    if(i3d.eq.1) call clear_mat(o1_mat)
    if((ipres.eq.1 .and. numvar.lt.3) .or. ipressplit.eq.1) then 
       call clear_mat(p1_mat)
    end if
    r4_vec = 0.
call PetscLogStagePop(jer)

    
    call clear_mat(s2_mat)
    call clear_mat(d2_mat)
    call clear_mat(r2_mat)
    call clear_mat(q2_mat)
    if(idens.eq.1) then
       call clear_mat(r42_mat)
       call clear_mat(q42_mat)
    end if
    if(irunaway .gt. 0) then
       call clear_mat(r43_mat)
       call clear_mat(q43_mat)
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
    if(irunaway .gt. 0) then
       call clear_mat(s15_mat)
       call clear_mat(d15_mat)
       call clear_mat(r15_mat)
       call clear_mat(q15_mat)
       call clear_mat(k15_mat)
       qn5_vec = 0.
    end if    

    if(ipres.eq.1 .or. ipressplit.eq.1) then
       call clear_mat(s9_mat)
       call clear_mat(d9_mat)
       call clear_mat(r9_mat)
       call clear_mat(q9_mat)
       call clear_mat(o9_mat)
       if(idens.eq.1) then
          call clear_mat(rp42_mat)
          call clear_mat(qp42_mat)
       end if
       qp4_vec = 0.
    end if
  end subroutine clear_matrices_split

  subroutine finalize_matrices_split
    use basic
    use m3dc1_vel_prof
    implicit none
    integer :: jer

call PetscLogStagePush(stageA,jer)
    call finalize(d1_mat)
    call finalize(q1_mat)
    call finalize(r14_mat)
    call flush(s1_mat)

    if(i3d.eq.1) call finalize(o1_mat)
    if((ipres.eq.1 .and. numvar.lt.3) .or. ipressplit.eq.1) then
       call finalize(p1_mat)
    end if
    call sum_shared(r4_vec)
call PetscLogStagePop(jer)

    if(myrank.eq.0 .and. iprint.ge.1) &
         print *, " before field finalize..."
    
    call finalize(d2_mat)
    call finalize(r2_mat)
    call finalize(q2_mat)
    call flush(s2_mat)
    if(i3d.eq.1) then
      call finalize(o2_mat)
      if(ipres.eq.1 .or. ipressplit.eq.1) call finalize(o3_mat)
    endif
    if(idens.eq.1) then
       call finalize(r42_mat)
       call finalize(q42_mat)
    end if
    if(irunaway .gt. 0) then
       call finalize(r43_mat)
       call finalize(q43_mat)
    end if
    call sum_shared(q4_vec)
    
    if(idens.eq.1) then
       call finalize(d8_mat)
       call finalize(q8_mat)
       call finalize(r8_mat)
       call flush(s8_mat)
       call sum_shared(qn4_vec)
       
       if(myrank.eq.0 .and. iprint.ge.1) &
            print *, " before pressure finalize...", vecsize_vel, vecsize_p, vecsize_phi   
    end if
       
    if(irunaway .gt. 0) then
       call finalize(d15_mat)
       call finalize(q15_mat)
       call finalize(r15_mat)
       call finalize(k15_mat)
       call flush(s15_mat)
       call sum_shared(qn5_vec)
    endif
    if(ipres.eq.1 .or. ipressplit.eq.1) then
       call finalize(d9_mat)
       call finalize(q9_mat)
       call finalize(r9_mat)
       call flush(s9_mat)
       if(myrank.eq.0 .and. iprint.ge.1) print *, 'before call to finalize(o9_mat)'
       call finalize(o9_mat)
       if(idens.eq.1) then
          call finalize(rp42_mat)
          call finalize(qp42_mat)
       end if

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
  use runaway_mod
  use runaway_advection

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
     if(ipres.eq.1 .and. numvar.ge.3) then
        pe_v = pe_field(1)
        ti_v = ti_field(1)
     endif
  endif   ! on ipressplit.eq.0

  den_v = den_field(1)
  ne_v = ne_field(1)
  if(irunaway .eq. 2) call runaway_advance
  if((irunaway .ge. 1).and.(ra_characteristics.eq.1)) call runaway_advection_step(dt*t0_norm)
  if(irunaway .gt. 0) nre_v = nre_field(1)
#ifdef USEPARTICLES
  if(itwofluid .gt. 0) densmooth_v = densmooth_field
#endif
  if(imp_bf.eq.1) bfp_v = bfp_field(1)
  if((jadv.eq.0) .or. (jadv.eq.1 .and. imp_hyper.ge.1)) e_v = e_field(1)

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

  if(iadapt .ne. 0) then
    u_field_pre = u_field(1)
    psi_field_pre = psi_field(1)
  end if
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
      if(ipres.eq.0 .or. (ipres.eq.1 .and. numvar.lt.3)) then
         if(itemp.eq.1) then
            ti_field(1) = te_v
            call mult(ti_field(1), (1.-pefac)/pefac)
         else
            pe_field(1) = p_v
            call mult(pe_field(1), pefac)
         end if
      else
         pe_field(1) = pe_v
         ti_field(1) = ti_v 
      endif

  endif   ! on ipressplit.eq.0

  if(idens.eq.1) then
     ne_field(1) = ne_v
     den_field(1) = den_v
  end if
  if(irunaway .gt. 0) then
     nre_field(1) = nre_v
  endif
  if(imp_bf.eq.1) bfp_field(1) = bfp_v
  if((jadv.eq.0) .or. (jadv.eq.1 .and. imp_hyper.ge.1)) e_field(1) = e_v

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
  use auxiliary_fields
  use m3dc1_vel_prof
  use hypervisc

  implicit none


  integer, intent(in) :: calc_matrices
  real :: tstart, tend, t_bound
  integer :: jer, i 
  type(vector_type) :: temp, temp2
  type(vector_type) :: xVec_guess

  type(field_type) :: phip_1, phip_2, phip_3

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
  
call PetscLogStagePush(stageA,jer)
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
        call matvecmult(o1_mat,bfp_field(1)%vec,b2_vel)
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
call PetscLogStagePop(jer)

     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_bound = t_bound + tend - tstart
     end if

     ! solve linear system with rhs in vtemp (note LU-decomp done first time)
     if(myrank.eq.0 .and. iprint.ge.2) print *, "  solving"
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

call PetscLogStagePush(stageS,jer)
     if(isolve_with_guess==1) then
        call create_vector(xVec_guess, vecsize_vel)
        xVec_guess = vel_vec ! Set the value from the previous step as the initial guess.
        call newsolve_with_guess(s1_mat, b1_vel, xVec_guess, jer)
        call destroy_vector(xVec_guess)
     else
        call newsolve(s1_mat, b1_vel, jer)
     endif 
call PetscLogStagePop(jer)
     !if(linear.eq.0) call clear_mat(s1_mat)

     if(idifv .gt.0) then
        call add(b1_vel,vel_vec)
     endif

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
     call apply_hyperv(u_v, vz_v, chi_v)
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

     if(isolve_with_guess==1) then
        call create_vector(xVec_guess, vecsize_n)
        xVec_guess = den_vec ! Set the value from the previous step as the initial guess.
        call newsolve_with_guess(s8_mat, temp, xVec_guess, jer)
        call destroy_vector(xVec_guess)
     else
        call newsolve(s8_mat, temp, jer)
     endif

     if(idiff .gt. 0) then
         call add(temp,den_vec)  ! add time n density to increment to get time n+1 value
     endif

     !if(linear.eq.0) call clear_mat(s8_mat)

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

     ! Electron temperature
     neold_vec = ne_vec
     call calculate_ne(1, den_v, ne_v, eqsubtract)
     ne_field(1) = ne_v    ! This is needed so that boundary_te works properly
     den_field(1) = den_v
  endif    ! on idens=1

  ! Advance RE Density
  ! ===============
  if(irunaway .gt. 0) then

     ! Here
     call create_vector(temp, vecsize_n)
     call create_vector(temp2, vecsize_n)

     nreold_vec = nre_vec
     do i=1,ra_cyc
        call matvecmult(q15_mat,phi_vec,temp)
        call matvecmult(k15_mat,veln_vec,temp2)
        call add(temp, temp2)

        call matvecmult(d15_mat,nre_vec,temp2)
        call add(temp, temp2)

        call add(temp, qn5_vec)

        ! Insert boundary conditions
        if(myrank.eq.0 .and. iprint.ge.2) print *, "  inserting bcs"
        if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
        if((calc_matrices.eq.1).and.(i==1)) then
           call boundary_nre(temp, nre_v, s15_mat)
           call finalize(s15_mat)
        else
           call boundary_nre(temp, nre_v)
        endif
        if(myrank.eq.0 .and. itimer.eq.1) then
           call second(tend)
           t_bound = t_bound + tend - tstart
        end if

        if(myrank.eq.0 .and. iprint.ge.2) print *, "  solving"
        ! solve linear system...LU decomposition done first time
        if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

        if(isolve_with_guess==1) then
           call create_vector(xVec_guess, vecsize_n)
           xVec_guess = nreold_vec ! Set the value from the previous step as the initial guess.
           call newsolve_with_guess(s15_mat, temp, xVec_guess, jer)
           call destroy_vector(xVec_guess)
        else
           call newsolve(s15_mat, temp, jer)
        endif    

        if(myrank.eq.0 .and. itimer.eq.1) then
           call second(tend)
           write(0,*) "Time spent", tend-tstart
           t_solve_n = t_solve_n + tend - tstart
        endif

        if(jer.ne.0) then
           write(*,*) 'Error in RE solve', jer
           call safestop(29)
        endif
        if(idiff .gt. 0) then
           call add(temp,nre_vec)  ! add time n RE density to increment to get time n+1 value
        endif
        ! new field solution at time n+1 (or n* for second order advance)
        nre_vec = temp
     enddo
     if(linear.eq.0) call clear_mat(s15_mat)

     if(myrank.eq.0 .and. iprint.ge.2) print *, "  done solve"

     call destroy_vector(temp)
     call destroy_vector(temp2)

     if(irecalc_eta.eq.1) then
        call export_time_advance_vectors_split
        call define_transport_coefficients
     end if

     nre_field(1) = nre_v     
  endif    ! on irunaway>0
  

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

     ! Inculde density terms
     if(idens.eq.1) then
        call matvecmult(rp42_mat,den_vec,temp2)
        call add(temp, temp2)
        call matvecmult(qp42_mat,denold_vec,temp2)
        call add(temp, temp2)
     end if

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
        if(ipressplit.eq.1 .and. ipres.eq.1 .and. numvar.ge.3) &
             call boundary_pe(temp, pe_v, s9_mat)
        call finalize(s9_mat)
     else
        call boundary_p(temp, p_v)
        if(ipressplit.eq.1 .and. ipres.eq.1 .and. numvar.ge.3) &
             call boundary_pe(temp, pe_v)
     endif
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_bound = t_bound + tend - tstart
     end if

     
     ! solve linear system...LU decomposition done first time
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Pressure--before newsolve"
     if(isolve_with_guess==1) then
        call create_vector(xVec_guess, vecsize_p)
        xVec_guess = pres_vec ! Set the value from the previous step as the initial guess.
        call newsolve_with_guess(s9_mat, temp, xVec_guess, jer)
        call destroy_vector(xVec_guess)
     else
        call newsolve(s9_mat, temp, jer)
     endif    

     !if(linear.eq.0) call clear_mat(s9_mat)

    if(idiff .gt. 0) then
         call add(temp,pres_vec)
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
     pres_vec = temp
     call destroy_vector(temp)
     call destroy_vector(temp2)
!     if(myrank.eq.0 .and. iprint.ge.1) print *, "Advancing Pressure--before get_temperatures"
     call export_time_advance_vectors_split
     call calculate_temperatures(1, te_field(1), ti_field(1), pe_field(1), p_field(1), ne_field(1), den_field(1), eqsubtract)
     if(numvar.ge.3.    .or.  ipres.eq.1) te_v = te_field(1)
     if(ipressplit.eq.1 .and. ipres.eq.1 .and. numvar.ge.3) ti_v = ti_field(1)
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

     ! o9matrix_sm * phi(n)
     call matvecmult(o9_mat,phi_vec,temp2)
     call add(temp, temp2)

     
     ! temp = d9matrix_sm * temp(n)
     call matvecmult(d9_mat,pret_vec,temp2)
     call add(temp, temp2)

 
     ! Include linear f terms
     if(i3d.eq.1 .and. numvar.ge.3 .and. (ipres.eq.1 .or. ipressplit.eq.1)) then
        call matvecmult(o3_mat,bfp_field(1)%vec,temp2)
        call add(temp, temp2)
     endif

    
     call add(temp, qp4_vec)
     
     ! Insert boundary conditions
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     if(calc_matrices.eq.1) then
        call boundary_te(temp, te_v, s9_mat)
        if(ipres.eq.1 .and. numvar.ge.3) call boundary_ti(temp, ti_v, s9_mat)
        call finalize(s9_mat)
     else
        call boundary_te(temp, te_v)
        if(ipres.eq.1 .and. numvar.ge.3) call boundary_ti(temp, ti_v)
     endif
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        t_bound = t_bound + tend - tstart
     end if

     
     ! solve linear system...LU decomposition done first time
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

     if(isolve_with_guess==1) then
       call create_vector(xVec_guess, vecsize_p)
       xVec_guess = pret_vec ! Set the value from the previous step as the initial guess.
       call newsolve_with_guess(s9_mat, temp, xVec_guess, jer)
       call destroy_vector(xVec_guess)
     else
       call newsolve(s9_mat, temp, jer)
     endif   
     !if(linear.eq.0) call clear_mat(s9_mat)

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
     call export_time_advance_vectors_split
     if(ipres.eq.1 .and. numvar.ge.3) then
        call calculate_pressures(1, pe_v, p_v, ne_field(1), &
             den_field(1), te_v, ti_v, eqsubtract)
     else
        call calculate_pressures(1, pe_field(1), p_v, ne_field(1), &
             den_field(1), te_v, ti_field(1), eqsubtract)
     end if

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
#ifdef USEPARTICLES
        call matvecmult(r42_mat,densmooth_vec,b2_phi)
        call mult(b2_phi,-1.)
        call add(b1_phi, b2_phi)
        call matvecmult(q42_mat,densmooth_vec,b2_phi)
        call add(b1_phi, b2_phi)
#else
        call matvecmult(r42_mat,ne_vec,b2_phi)
        call add(b1_phi, b2_phi)
        call matvecmult(q42_mat,neold_vec,b2_phi)
        call add(b1_phi, b2_phi)
!!$        call matvecmult(r42_mat,den_vec,b2_phi)
!!$        call add(b1_phi, b2_phi)
!!$        call matvecmult(q42_mat,denold_vec,b2_phi)
!!$        call add(b1_phi, b2_phi)
#endif
     end if

     ! Include RE density terms
     if(irunaway .gt. 0) then
        call matvecmult(r43_mat,nre_vec,b2_phi)
        call mult(b2_phi,-1.)
        call add(b1_phi, b2_phi)
        call matvecmult(q43_mat,nreold_vec,b2_phi)
        call add(b1_phi, b2_phi)
     end if

     ! Include linear f terms
     if(numvar.ge.2 .and. i3d.eq.1 .and. imp_bf.eq.0) then
        ! b2vector = r15 * bf(n)
        call matvecmult(o2_mat,bfp_field(1)%vec,b2_phi)
        call add(b1_phi, b2_phi)
     endif

     call add(b1_phi, q4_vec)
 
     ! Insert boundary conditions
     if(myrank.eq.0 .and. iprint.ge.2) print *, "  inserting bcs"
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     if(calc_matrices.eq.1) then
        call boundary_mag(b1_phi, psi_v, bz_v, bfp_v, e_v, s2_mat)
        if(ipressplit.eq.0 .and. numvar.ge.3) then
           if(ipres.eq.1) then
              call boundary_pe(b1_phi, pe_v, s2_mat)
           else
              call boundary_p(b1_phi, pe_v, s2_mat)
           endif
        endif
        call finalize(s2_mat)
     else 
        call boundary_mag(b1_phi, psi_v, bz_v, bfp_v, e_v)
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

     if(isolve_with_guess==1) then
        call create_vector(xVec_guess, vecsize_phi)
        xVec_guess = phi_vec ! Set the value from the previous step as the initial guess.
        call newsolve_with_guess(s2_mat, b1_phi, xVec_guess, jer)
        call destroy_vector(xVec_guess)
     else
        call newsolve(s2_mat, b1_phi, jer)
     endif    
     !if(linear.eq.0 .and. iteratephi.eq.0) call clear_mat(s2_mat)


   if(idiff .gt. 0) then
         call add(b1_phi,phi_vec)
   endif

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
#ifdef USEPARTICLES
           call matvecmult(r42_mat,densmooth_vec,b2_phi)
           call mult(b1_phi,-1.)
           call add(b1_phi, b2_phi)
           call matvecmult(q42_mat,densmooth_vec,b2_phi)
           call add(b1_phi, b2_phi)
#else
           call matvecmult(r42_mat,den_vec,b2_phi)
           call add(b1_phi, b2_phi)
           call matvecmult(q42_mat,denold_vec,b2_phi)
           call add(b1_phi, b2_phi)
#endif
        end if

        ! Include linear f terms
        if(numvar.ge.2 .and. i3d.eq.1 .and. imp_bf.eq.0) then
           ! b2vector = r15 * bf(n)
           call matvecmult(o2_mat,bfp_field(1)%vec,b2_phi)
           call add(b1_phi, b2_phi)
        endif

        ! Include RE density terms
        if(irunaway .gt. 0) then
           call matvecmult(r43_mat,nre_vec,b2_phi)
           call mult(b2_phi,-1.)
           call add(b1_phi, b2_phi)
           call matvecmult(q43_mat,nreold_vec,b2_phi)
           call add(b1_phi, b2_phi)
        end if

        call add(b1_phi, q4_vec)
       
        ! Insert boundary conditions
        if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
        if(calc_matrices.eq.1) then
           call boundary_mag(b1_phi, psi_v, bz_v, bfp_v, e_v, s2_mat)
           if(ipressplit.eq.0 .and. numvar.ge.3) then
              if(ipres.eq.1) then
                 call boundary_pe(b1_phi, pe_v, s2_mat)
              else
                 call boundary_p(b1_phi, pe_v, s2_mat)
              endif
           endif

           call finalize(s2_mat)
        else 
           call boundary_mag(b1_phi, psi_v, bz_v, bfp_v, e_v)
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
        
        !if(linear.eq.0) call clear_mat(s2_mat)
        if(isolve_with_guess==1) then
            call create_vector(xVec_guess, vecsize_phi)
            xVec_guess = phi_vec ! Set the value from the previous step as the initial guess.
            call newsolve_with_guess(s2_mat, b1_phi, xVec_guess, jer)
            call destroy_vector(xVec_guess)
        else
            call newsolve(s2_mat, b1_phi, jer)
        endif    
        
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
    integer:: l,numnodes, icounter_t
    vectype, dimension (dofs_per_node) :: vec_l, axi_l
    type(field_type) :: axi

    call create_field(axi)
    numnodes = owned_nodes()

    do  icounter_t=1,numnodes
        l = nodes_owned(icounter_t)
        call get_node_data(u_v,l,vec_l)
        call set_node_data(axi,l,vec_l)
    enddo

    call finalize(axi%vec)
    call m3dc1_field_sum_plane(axi%vec%id)

    do  icounter_t=1,numnodes
        l = nodes_owned(icounter_t)
        call get_node_data(axi,l,axi_l)
        call get_node_data(u_v,l,vec_l)
        vec_l = vec_l - axi_l/nplanes
        call set_node_data(u_v,l,vec_l)
    enddo
    call finalize(u_v%vec)
    if(numvar.lt.2) return
    
    do  icounter_t=1,numnodes
        l = nodes_owned(icounter_t)
        call get_node_data(vz_v,l,vec_l)
        call set_node_data(axi,l,vec_l)
    enddo

    call finalize(axi%vec)
    call m3dc1_field_sum_plane(axi%vec%id)

    do  icounter_t=1,numnodes
        l = nodes_owned(icounter_t)
        call get_node_data(axi,l,axi_l)
        call get_node_data(vz_v,l,vec_l)
        vec_l = vec_l - axi_l/nplanes
        call set_node_data(vz_v,l,vec_l)
    enddo
    call finalize(vz_v%vec)
    if(numvar.lt.3) return

    do  icounter_t=1,numnodes
        l = nodes_owned(icounter_t)
        call get_node_data(chi_v,l,vec_l)
        call set_node_data(axi,l,vec_l)
    enddo

    call finalize(axi%vec)
    call m3dc1_field_sum_plane(axi%vec%id)

    do  icounter_t=1,numnodes
        l = nodes_owned(icounter_t)
        call get_node_data(axi,l,axi_l)
        call get_node_data(chi_v,l,vec_l)
        vec_l = vec_l - axi_l/nplanes
        call set_node_data(chi_v,l,vec_l)
    enddo
    call finalize(chi_v%vec)

end subroutine subtract_axi
end module time_step_split
