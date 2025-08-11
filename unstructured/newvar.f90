module newvar_mod

  use matrix_mod

  implicit none

  integer, parameter :: NV_NOBOUND = 0
  integer, parameter :: NV_DCBOUND = 1  ! Dirichlet boundary conditions
  integer, parameter :: NV_NMBOUND = 2  ! Neumann boundary conditions
  integer, parameter :: NV_SJBOUND = 3
  integer, parameter :: NV_CYBOUND = 6  ! Cauchy boundary conditions

  integer, parameter :: NV_I_MATRIX = 0
  integer, parameter :: NV_LP_MATRIX = 1
  integer, parameter :: NV_GS_MATRIX = 2
  integer, parameter :: NV_BF_MATRIX = 3
  integer, parameter :: NV_SJ_MATRIX = 4  ! Current density smoother
  integer, parameter :: NV_DP_MATRIX = 7
#if defined(USECOMPLEX) || defined(USE3D)
  integer, parameter :: NV_IP_MATRIX = 8
#endif

  integer, parameter :: NV_RHS = 0
  integer, parameter :: NV_LHS = 1

  type newvar_matrix
     type(matrix_type) :: mat
     integer :: ibound
  end type newvar_matrix

  type(newvar_matrix) :: mass_mat_lhs
  type(newvar_matrix) :: mass_mat_lhs_dc
  type(newvar_matrix) :: mass_mat_rhs
!  type(newvar_matrix) :: lp_mat_rhs
!  type(newvar_matrix) :: lp_mat_rhs_dc
  type(newvar_matrix) :: gs_mat_rhs_dc
  type(newvar_matrix) :: gs_mat_rhs
  type(newvar_matrix) :: bf_mat_lhs
  type(newvar_matrix) :: pot2_mat_lhs
  type(newvar_matrix) :: mass_mat_rhs_bf
  type(newvar_matrix) :: dp_mat_rhs_bfp
  type(newvar_matrix) :: s10_mat, d10_mat

contains

  !================================================================
  ! set_newvar_indices
  ! ~~~~~~~~~~~~~~~~~~
  ! Sets the scorec indices of the newvar matrices
  ! This must be called before matrices are created or used
  !================================================================
  subroutine set_newvar_indices
    use sparse

    implicit none

    call set_matrix_index(mass_mat_lhs%mat,    mass_mat_lhs_index)
    call set_matrix_index(mass_mat_rhs%mat,    mass_mat_rhs_index)
    call set_matrix_index(mass_mat_lhs_dc%mat, mass_mat_lhs_dc_index)
!!$    call set_matrix_index(lp_mat_rhs%mat,      lp_mat_rhs_index)
!!$    call set_matrix_index(lp_mat_rhs_dc%mat,   lp_mat_rhs_dc_index)
    call set_matrix_index(gs_mat_rhs_dc%mat,   gs_mat_rhs_dc_index)
    call set_matrix_index(gs_mat_rhs%mat,      gs_mat_rhs_index)
    call set_matrix_index(bf_mat_lhs%mat,      bf_mat_lhs_dc_index)
    call set_matrix_index(mass_mat_rhs_bf%mat, mass_mat_rhs_dc_index)
    call set_matrix_index(dp_mat_rhs_bfp%mat,  bfp_mat_rhs_index)
    call set_matrix_index(s10_mat%mat,         s10_mat_index)
    call set_matrix_index(d10_mat%mat,         d10_mat_index)
    call set_matrix_index(pot2_mat_lhs%mat,    pot2_mat_lhs_index)

  end subroutine set_newvar_indices


  !==========================================================
  ! create_newvar_matrices
  ! ~~~~~~~~~~~~~~~~~~~~~~
  ! populate all of the newvar matrices that will be used
  !==========================================================
  subroutine create_newvar_matrices
    use sparse
    use basic

    ! assign the proper reference indicies to each matrix
    call set_newvar_indices

    call create_newvar_matrix(mass_mat_lhs_dc, NV_DCBOUND,NV_I_MATRIX, 1)
    call create_newvar_matrix(mass_mat_lhs,    NV_NOBOUND,NV_I_MATRIX, 1)
    call create_newvar_matrix(gs_mat_rhs_dc,   NV_DCBOUND,NV_GS_MATRIX,0)

#if defined(USE3D) || defined(USECOMPLEX)
    call create_newvar_matrix(dp_mat_rhs_bfp,  NV_NOBOUND,NV_IP_MATRIX,0)
#endif

    if(inocurrent_tor.eq.0) then 
       call create_newvar_matrix(gs_mat_rhs,  NV_NOBOUND,NV_GS_MATRIX, 0)
    endif


    if((i3d.eq.1 .or. ifout.eq.1) .and. numvar.ge.2) then
       if(ifbound.eq.1) then 
          call create_newvar_matrix(bf_mat_lhs, &
               NV_DCBOUND, NV_BF_MATRIX, 1)
          call create_newvar_matrix(mass_mat_rhs_bf, NV_DCBOUND, &
               NV_I_MATRIX,  0)
       else if(ifbound.eq.2) then 
          call create_newvar_matrix(bf_mat_lhs, &
               NV_NMBOUND, NV_BF_MATRIX, 1)
          call create_newvar_matrix(mass_mat_rhs_bf, NV_NMBOUND, &
               NV_I_MATRIX,  0)
       end if
    endif

    if(jadv.eq.1 .and. hyper.ne.0. .and. imp_hyper.eq.0) then
       call create_newvar_matrix(s10_mat, NV_SJBOUND, NV_SJ_MATRIX, 1)
       call create_newvar_matrix(d10_mat, NV_SJBOUND, NV_SJ_MATRIX, 0)
    endif
        call create_newvar_matrix(pot2_mat_lhs,NV_DCBOUND, NV_DP_MATRIX,1)

    if(igs.ne.0) then
       call create_newvar_matrix(mass_mat_rhs,NV_NOBOUND,NV_I_MATRIX, 0)
    endif
  end subroutine create_newvar_matrices

subroutine apply_bc(rhs, ibound, bvec, mat)
  use vector_mod
  use matrix_mod
  use boundary_conditions

  implicit none

  type(matrix_type), optional :: mat
  integer, intent(in) :: ibound
  type(vector_type) :: rhs
  type(vector_type) :: bvec

  if(ibound.eq.NV_NOBOUND) return

  select case(ibound)
  case(NV_DCBOUND)
     call boundary_dc(rhs, bvec, mat)
  case(NV_NMBOUND)
     call boundary_nm(rhs, bvec, mat)
  case(NV_SJBOUND)
     call boundary_jphi(rhs, mat)
  case(NV_CYBOUND)
     call boundary_cy(rhs, mat)
  end select
end subroutine apply_bc


!============================================
! create_newvar_matrix
! ~~~~~~~~~~~~~~~~~~~~
! creates a matrix.
! ibound:
!   NV_NOBOUND: no boundary conditions
!   NV_DCBOUND: Dirichlet boundary conditions
!   NV_NMBOUND: Neumann boundary conditions
! itype: operator (NV_I_MATRIX, etc..)
!============================================
#ifdef REORDERED
  subroutine create_newvar_matrix(mat, ibound, itype, is_lhs, tags, agg_blk_cnt, agg_scp)
#else
  subroutine create_newvar_matrix(mat, ibound, itype, is_lhs, tags)
#endif
  use vector_mod
  use basic
  use m3dc1_nint
  use boundary_conditions
  use m3dc1_omp

  implicit none

  type(newvar_matrix), intent(out) :: mat
  integer, intent(in) :: ibound
  integer, intent(in) :: itype
  integer, intent(in) :: is_lhs
#ifdef REORDERED
  integer, intent(in), optional :: agg_blk_cnt
  integer, intent(in), optional :: agg_scp
#endif
  integer, intent(in), optional :: tags

  integer :: numelms, itri, m, n, isize
  vectype, allocatable :: temp(:,:,:,:)
  type(vector_type) :: rhs2

  ! determine the size of the matrix
  select case(itype)
  case(NV_SJ_MATRIX)
     isize = 2
  case default
     isize = 1
  end select

#ifdef REORDERED
  call create_mat(mat%mat, isize, isize, icomplex, is_lhs, agg_blk_cnt, agg_scp)
#else
  call create_mat(mat%mat, isize, isize, icomplex, is_lhs)
#endif
  mat%ibound = ibound

  allocate(temp(dofs_per_element, dofs_per_element, isize, isize))

  if(myrank.eq.0 .and. iprint.ge.2) print *, ' populating matrix...', &
       mat%mat%m, mat%mat%n
  numelms = local_elements()

!!$OMP PARALLEL DO &
!!$OMP& PRIVATE(temp)
  do itri=1,numelms

     call define_element_quadrature(itri, int_pts_main, 5)
     call define_fields(itri,0,1,0)

     temp = 0.

     selectcase(itype)
        
     case(NV_I_MATRIX)
        temp(:,:,1,1) = intxx2(mu79(:,:,OP_1),nu79(:,:,OP_1))
        
#if defined(USE3D)
     case(NV_IP_MATRIX)
        temp(:,:,1,1) = intxx2(mu79(:,:,OP_1),nu79(:,:,OP_DP))
#elif defined(USECOMPLEX)
     case(NV_IP_MATRIX)
        temp(:,:,1,1) = intxx2(mu79(:,:,OP_1),nu79(:,:,OP_1 ))*rfac
#endif

     case(NV_LP_MATRIX)
        temp(:,:,1,1) = intxx2(mu79(:,:,OP_1),nu79(:,:,OP_LP))
        if(ibound.eq.NV_NMBOUND) then
           temp(:,:,1,1) = temp(:,:,1,1) &
                - regular*intxx2(mu79(:,:,OP_1),nu79(:,:,OP_1))
        endif
              
     case(NV_GS_MATRIX)              
        temp(:,:,1,1) = intxx2(mu79(:,:,OP_1),nu79(:,:,OP_GS))
              
     case(NV_BF_MATRIX)
        temp(:,:,1,1) = intxx3(mu79(:,:,OP_1),nu79(:,:,OP_LP),r2_79)
        if(ifbound.eq.2) then
           temp(:,:,1,1) = temp(:,:,1,1) + &
                regular*intxx2(mu79(:,:,OP_1),nu79(:,:,OP_1))
        end if
              
     case(NV_SJ_MATRIX)
        if(is_lhs .eq. 1) then
           temp(:,:,1,1) = intxx2(mu79(:,:,OP_1),nu79(:,:,OP_1))
           temp(:,:,1,2) = -thimpsm*intxx2(mu79(:,:,OP_1),nu79(:,:,OP_GS))
           temp(:,:,2,1) = dt*hypf*intxx2(mu79(:,:,OP_1),nu79(:,:,OP_GS))
           temp(:,:,2,2) = temp(:,:,1,1)
        else
           temp(:,:,1,2) = (1.-thimpsm)*&
                intxx2(mu79(:,:,OP_1),nu79(:,:,OP_GS))
           temp(:,:,2,2) = intxx2(mu79(:,:,OP_1),nu79(:,:,OP_1))
        end if
        
     case(NV_DP_MATRIX)              
        temp(:,:,1,1) = intxx2(mu79(:,:,OP_DZ),nu79(:,:,OP_DZ)) &
             +          intxx2(mu79(:,:,OP_DR),nu79(:,:,OP_DR))
        
     end select

     select case(ibound)
     case(NV_DCBOUND)
        call apply_boundary_mask(itri, BOUNDARY_DIRICHLET, temp(:,:,1,1), &
             tags=tags)
     case(NV_NMBOUND)
        call apply_boundary_mask(itri, BOUNDARY_NEUMANN, temp(:,:,1,1), &
             tags=tags)
     case(NV_SJBOUND)
        call apply_boundary_mask(itri, BOUNDARY_DIRICHLET, temp(:,:,1,1), &
             tags=tags)
        call apply_boundary_mask(itri, BOUNDARY_DIRICHLET, temp(:,:,1,2), &
             tags=tags)
        call apply_boundary_mask(itri, BOUNDARY_DIRICHLET, temp(:,:,2,1), &
             tags=tags)
        call apply_boundary_mask(itri, BOUNDARY_DIRICHLET, temp(:,:,2,2), &
             tags=tags)

     case(NV_CYBOUND)
        call apply_boundary_mask(itri, & 
             BOUNDARY_DIRICHLET + BOUNDARY_NEUMANN, temp(:,:,1,1), tags=tags)

     end select

!!$OMP CRITICAL
     do m=1,isize
        do n=1,isize
           call insert_block(mat%mat,itri,m,n,temp(:,:,m,n),MAT_ADD)
        end do
     end do
!!$OMP END CRITICAL
  end do
!!$OMP END PARALLEL DO

  deallocate(temp)

  ! apply boundary conditions
  if(is_lhs .eq. 1 .and. ibound.ne.NV_NOBOUND) then
     if(myrank.eq.0 .and. iprint.ge.2) print *, ' applying bcs...'
     call flush(mat%mat)
     call create_vector(rhs2, isize)
     call apply_bc(rhs2, ibound, rhs2, mat%mat)
     call destroy_vector(rhs2)
  end if

  call finalize(mat%mat)

  if(myrank.eq.0 .and. iprint.ge.2) print *, ' done.'

end subroutine create_newvar_matrix


!=====================================================================
! solve_newvar_axby
! ~~~~~~~~~~~~~~~~~
! creates rhs and solves matrix equation:
!  A x = B y
!
! mata: lhs matrix (A)
! vout: vector to contain the result (x)
! vin:  vector y
! matb: rhs matrix(B)
!======================================================================
subroutine solve_newvar_axby(mata,vout,matb,vin,bvec)
  use vector_mod

  implicit none

  type(newvar_matrix), intent(in) :: mata, matb
  type(vector_type), intent(inout), target :: vout
  type(vector_type), intent(in) :: vin
  type(vector_type), target, optional :: bvec
  type(vector_type) :: temp

  type(vector_type), pointer :: bptr
  integer :: ier

  if(present(bvec)) then
     bptr => bvec
  else
     bptr => vout
  end if

  call create_vector(temp, matb%mat%m)
  call matvecmult(matb%mat, vin, temp)
  call apply_bc(temp,mata%ibound,bptr)

  call newsolve(mata%mat,temp,ier)

  if(ier.ne.0) then
     print *, 'Error in newvar solve'
     call safestop(10)
  endif

  vout = temp
  call destroy_vector(temp)
  
end subroutine solve_newvar_axby


!=====================================================================
! solve_newvar1
! ~~~~~~~~~~~~~
! creates rhs and solves matrix equation:
!  A x = B y
! for two fields x and y (i.e. matrix eqn is rank 1)
!
! mata: lhs matrix (A)
! fout: field x
! fin: field y
! matb: rhs matrix(B)
!======================================================================
subroutine solve_newvar1(lhsmat, fout, rhsmat, fin, bvec)

  use field

  implicit none

  type(newvar_matrix), intent(in) :: lhsmat, rhsmat
  type(field_type), intent(in), target :: fin
  type(field_type), intent(inout), target :: fout
  type(field_type), intent(in), target, optional :: bvec

  type(field_type), target :: temp_in, temp_out
  type(vector_type), pointer :: bptr, ptrin, ptrout

  ! if inarray is bigger than vecsize=1, then 
  ! create vecsize=1 for matrix multiplication
  if(fin%vec%isize.gt.1) then
     call create_field(temp_in)
     temp_in = fin
     ptrin => temp_in%vec
  else
     ptrin => fin%vec
  end if

  if(fout%vec%isize.gt.1) then
     call create_field(temp_out)
     temp_out = fout
     ptrout => temp_out%vec
  else
     ptrout => fout%vec
  end if

  if(present(bvec)) then
     bptr => bvec%vec
  else
     bptr => ptrout
  end if 

  call solve_newvar_axby(lhsmat,ptrout,rhsmat,ptrin,bptr)

  if(fin%vec%isize.gt.1)  call destroy_field(temp_in)
  if(fout%vec%isize.gt.1) then
     fout = temp_out
     call destroy_field(temp_out)
  endif
end subroutine solve_newvar1


  !=====================================================
  ! newvar_solve
  ! ~~~~~~~~~~~~
  ! Solves equation mat*x = rhs
  ! with ibound boundary conditions applied to rhs.
  ! rhs is overwritten with result (x) on return.
  ! Be sure imatrix was also generated using ibound.
  !=====================================================
  subroutine newvar_solve(rhs, mat, bvec)
    use vector_mod

    implicit none

    type(vector_type), intent(inout), target :: rhs
    type(newvar_matrix), intent(in) :: mat
    type(vector_type), intent(inout), target, optional :: bvec

    integer :: ier
    type(vector_type), pointer :: bptr

    if(.not.present(bvec)) then
       bptr => rhs
    else
       bptr => bvec
    end if
    call sum_shared(rhs)

    call apply_bc(rhs,mat%ibound,bptr)

    call newsolve(mat%mat,rhs,ier)

    call finalize(rhs)
  end subroutine newvar_solve

  subroutine newvar_solve_with_guess(rhs, xVec_guess, mat, bvec)
    use vector_mod

    implicit none

    type(vector_type), intent(inout), target :: rhs, xVec_guess
    type(newvar_matrix), intent(in) :: mat
    type(vector_type), intent(inout), target, optional :: bvec

    integer :: ier
    type(vector_type), pointer :: bptr

    if(.not.present(bvec)) then
       bptr => rhs
    else
       bptr => bvec
    end if
    call sum_shared(rhs)

    call apply_bc(rhs,mat%ibound,bptr)

    call newsolve_with_guess(mat%mat,rhs,xVec_guess,ier)

    call finalize(rhs)
  end subroutine newvar_solve_with_guess

end module newvar_mod
