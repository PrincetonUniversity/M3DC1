module adapt
  use vector_mod
  use scorec_adapt
  implicit none
  real :: adapt_ke
  integer :: iadapt_ntime
  real :: adapt_target_error
  integer :: iadapt_max_node
  integer :: adapt_control
  real :: iadapt_order_p
  data iadapt_order_p /3/
  data iadapt_max_node /100/
  real :: error_tol = 1.
  real , dimension(2) :: abs_size, rel_size
  real :: adapt_hmin_rel, adapt_hmax_rel
  data rel_size /0.5, 2.0/

  real :: adapt_coil_delta
  real :: adapt_pellet_length, adapt_pellet_delta
  
  real :: adapt_zlow, adapt_zup
  logical :: do_z_coarsen
  integer :: iadaptFaceNumber
  integer, parameter :: maxqs = 32
  real, dimension(maxqs) :: adapt_qs
 
  contains

  subroutine adapt_by_psi
    use basic
    use mesh_mod
    use arrays
    use newvar_mod
    use sparse
    use hdf5_output
    use diagnostics
    use boundary_conditions
    use time_step
    use m3dc1_output
    use auxiliary_fields
    use pellet
    use scorec_mesh_mod
    use m3dc1_nint
    use coils
    use transport_coefficients

    integer :: izone, izonedim, i, j
    integer :: numelms, itri
    character(len=32) :: mesh_file_name
    real :: q, tmp
    integer, dimension(MAX_PTS) :: mr
    vectype, dimension(dofs_per_element) :: dofs

    real, dimension(maxfilaments) :: xc_adapt, zc_adapt
    complex, dimension(maxfilaments) :: ic_adapt
    integer :: numcoils_adapt

    real, allocatable :: xp_adapt(:,:), zp_adapt(:,:)
    integer :: p_steps
    real :: p_dt, p_v
    real :: x0, y0, x, y
    integer :: ip

    real :: psib
    call create_field(temporary_field)

    if (ispradapt .eq. 1) then
     call m3dc1_field_mark4tx(temporary_field)
    endif

    temporary_field = 0.

    if(adapt_pellet_delta.gt.0) then

       allocate(xp_adapt(maxfilaments,npellets))
       allocate(zp_adapt(maxfilaments,npellets))
       
       ! determine pellet path to adapt along
       do ip=1, npellets
          x0 = pellet_r(ip)*cos(pellet_phi(ip))
          y0 = pellet_r(ip)*sin(pellet_phi(ip))
          p_steps = ceiling(adapt_pellet_length/(2.*adapt_pellet_delta))+1
          p_steps = min(p_steps,maxfilaments)
          p_v = sqrt(pellet_vx(ip)**2 + pellet_vy(ip)**2 + pellet_velz(ip)**2)
          if(p_v.gt.0.) then
             p_dt = 2.*adapt_pellet_delta/p_v
          else
             ! if not moving, just adapt along stationary pellet position
             p_dt = 0.
             p_steps = 1
          end if

          do j=1, p_steps
             x = x0 + pellet_vx(ip)*(j-1)*p_dt
             y = y0 + pellet_vy(ip)*(j-1)*p_dt
             xp_adapt(j,ip) = sqrt(x**2 + y**2)
             zp_adapt(j,ip) = pellet_z(ip) + pellet_velz(ip)*(j-1)*p_dt
          end do
       end do
    end if

    if(adapt_coil_delta.gt.0) then
       call load_coils(xc_adapt,zc_adapt,ic_adapt,numcoils_adapt, &
            'adapt_coil.dat','adapt_current.dat')
    end if

    numelms = local_elements()
    do itri=1,numelms
!       call zonfac(itri,izone,izonedim)
       call m3dc1_ent_getgeomclass(2, itri-1, izonedim, izone)

       call define_element_quadrature(itri,int_pts_main,int_pts_tor)
       call define_fields(itri,0,1,0)

       if(eqsubtract.eq.1) then 
          call eval_ops(itri, psi_field(0), ps079)
       else
          call eval_ops(itri, psi_field(1), ps079)
       end if
       if(icsubtract.eq.1) then 
          call eval_ops(itri, psi_coil_field, ps179)
          ps079 = ps079 + ps179
       end if

       ! store psi_N in temp79a
       temp79a = (ps079(:,OP_1) - psimin) / (psibound - psimin)

       ! let temp79b store the "modified" psi_N
       temp79b = temp79a
       
       ! determine magnetic region of each point
       do i=1, npoints
          call magnetic_region(ps079(i,OP_1),ps079(i,OP_DR),ps079(i,OP_DZ), &
               x_79(i),z_79(i),mr(i),psib)

          do_z_coarsen = ((adapt_zup.ne.0.).and.(z_79(i).gt.adapt_zup)).or. &
                         ((adapt_zlow.ne.0.).and.(z_79(i).lt.adapt_zlow))
         
          ! if point is in private flux region, set psi_N -> 2 - psi_N
          if(mr(i).eq.REGION_PF) then
             temp79b(i) = 2.*psib - temp79a(i)
             if(do_z_coarsen) temp79b(i) = 1. + (temp79b(i) - 1.)**0.5
          else if((mr(i).eq.REGION_SOL).and.(do_z_coarsen)) then
             temp79b(i) = 1. + (temp79a(i) - 1.)**0.75
          else if(do_z_coarsen) then
             temp79b(i) = 1. - (1. - temp79a(i))**0.5
          end if

       end do

       ! if adapt_psin_wall or adapt_psin_vacuum is set in multi-region mesh,
       ! set psin in the wall or vacuum region to the appropriate value
       if(imulti_region.eq.1) then
          if(izone.eq.2 .and. adapt_psin_wall.ge.0) then
             temp79b = adapt_psin_wall
          end if
          if(izone.ge.3 .and. adapt_psin_vacuum.ge.0) then
             temp79b = adapt_psin_vacuum
          end if
       end if

       ! change psi_N inside plasma to be sum of Lorentzians around
       ! rational surfaces
       if(iadapt_pack_rationals.gt.0 .and. ntor.ne.0) then
          if(.not.allocated(q_spline%y)) then
             print *, 'Error, iadapt_pack_rationals > 0, but q profile not set'
             call safestop(6)
          end if
          do i=1, npoints
             
             if(mr(i).ne.0) cycle

             call evaluate_spline(q_spline, real(temp79a(i)), q)

             temp79b(i) = 0.
             do j=1, iadapt_pack_rationals
                tmp = (q*ntor - j) / (q*ntor*adapt_pack_factor)
                temp79b(i) = temp79b(i) + 1./(tmp**2 + 1.)
                if(real(temp79b(i)).gt.1.) temp79b(i) = 1.
             end do
          end do
       end if

       ! do adaptation around specified safety factor values
       if(any(adapt_qs.ne.0.)) then
          if(.not.allocated(q_spline%y)) then
             print *, 'Error, iadapt_pack_rationals > 0, but q profile not set'
             call safestop(6)
          end if

          temp79c = 0.
          do i=1, npoints

             if(mr(i).ne.0) cycle

             call evaluate_spline(q_spline, real(temp79a(i)), q)

             do j=1, maxqs
                if(adapt_qs(j).ne.0.) then
                   temp79c(i) = temp79c(i) + &
                        exp(-(q - adapt_qs(j))**2 / (2.*adapt_pack_factor**2))
                end if
             end do

          end do
          where(real(temp79c).gt.1.) temp79c = 1.
          temp79b = temp79b*(1.-temp79c) + temp79c
       end if


       ! do adaptation around coils
       if(adapt_coil_delta.gt.0) then
          temp79c = 0.
          do j=1, numcoils_adapt
             temp79c = temp79c + &
                  exp(-((x_79 - xc_adapt(j))**2 + (z_79 - zc_adapt(j))**2) / &
                       (2.*adapt_coil_delta**2))
          end do
          where(temp79c.ne.temp79c) temp79c = 0.
          where(real(temp79c).gt.1.) temp79c = 1.
          where(real(temp79c).lt.0.) temp79c = 0.
          temp79b = temp79b*(1.-temp79c) + temp79c
       end if

       ! do adaptation along pellet path
       if(adapt_pellet_delta.gt.0) then
          temp79c = 0.
          do ip = 1, npellets
             do j = 1, p_steps
                temp79c = temp79c + &
                     exp(-((x_79 - xp_adapt(j,ip))**2 + (z_79 - zp_adapt(j,ip))**2) / &
                         (2.*adapt_pellet_delta**2))
             end do
          end do
          !BCL: this normalizes the overlapping Gaussians
          if(p_steps.eq.2) temp79c = temp79c/exp(-2.)
          if(p_steps.gt.2) temp79c = temp79c/(2.*exp(-2.))
          where(real(temp79c).gt.1.) temp79c = 1.
          temp79b = temp79b*(1.-temp79c) + temp79c
       end if

       ! convert back to un-normalized psi
       temp79b = (psibound - psimin)*temp79b + psimin

       do i=1, dofs_per_element
          dofs(i) = int2(mu79(i,:,OP_1),temp79b)
       end do
       call vector_insert_block(temporary_field%vec,itri,1,dofs,VEC_ADD)
    end do

    if(adapt_pellet_delta.gt.0) then
       deallocate(xp_adapt)
       deallocate(zp_adapt)
    end if

    call newvar_solve(temporary_field%vec,mass_mat_lhs)

    call straighten_fields()

    if (iadaptFaceNumber.gt.0) then
        call adapt_model_face(temporary_field%vec%id,psimin,psibound,iadaptFaceNumber)
    else
        call adapt_by_field(temporary_field%vec%id,psimin,psibound)
    endif

    write(mesh_file_name,"(A7,A)") 'adapted', 0
    if(iadapt_writevtk .eq. 1) call m3dc1_mesh_write (mesh_file_name,0,ntime)
    if(iadapt_writesmb .eq. 1) call m3dc1_mesh_write (mesh_file_name,1,ntime)

    call destroy_field(temporary_field)
    call space(0)
    call update_nodes_owned()
    call tridef
    call unstraighten_fields()

    call create_newvar_matrices
    if(irestart .ne. 0) return
    field_vec = 0.
    field0_vec = 0.
    if (myrank .eq. 0) print *, "re-calculate equlibrium after adapt .."
    call initial_conditions
    ! combine the equilibrium and perturbed fields of linear=0
    ! unless eqsubtract = 1
    if(eqsubtract.eq.0) then
       call add(field_vec, field0_vec)
       field0_vec = 0.
    endif
    i_control%err_i = 0.
    i_control%err_p_old = 0.
    n_control%err_i = 0.
    n_control%err_p_old = 0.
    i_control%err_i = 0.
    i_control%err_p_old = 0.
    n_control%err_i = 0.
    n_control%err_p_old = 0.
    if(myrank.eq.0 .and. iprint.ge.2) print *, "  transport coefficients"
    call define_transport_coefficients
    call derived_quantities(1)
    !ke_previous = ekin
  end subroutine adapt_by_psi

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine adapt_by_spr(fid,idx,t,ar,maxsize,refinelevel,coarsenlevel,update)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    use diagnostics
    use basic
    use error_estimate
    use scorec_mesh_mod
    use basic
    use mesh_mod
    use arrays
    use newvar_mod
    use sparse
    use time_step
    use auxiliary_fields
    use scorec_mesh_mod
    use transport_coefficients

    character(len=32) :: mesh_file_name
    integer, intent(in) :: fid
    integer, intent(in) :: idx
    real, intent(in) :: ar
    real, intent(in) :: maxsize
    integer, intent(in) :: t
    integer, intent(in) :: refinelevel
    integer, intent(in) :: coarsenlevel
    logical, intent(in) :: update


  !  if (myrank .eq. 0) print *, " error exceeds tolerance, start adapting mesh"
    call straighten_fields()
    call m3dc1_spr_adapt(fid,idx,t,ar,maxsize,refinelevel,coarsenlevel,update)
    call space(0)
    call update_nodes_owned()
    call reset_itris()
    call tridef
    call unstraighten_fields()

    call create_newvar_matrices
    if(irestart .ne. 0 .or. linear .eq. 0) return
    if (myrank .eq. 0) print *, "reset simulation after adapt .."
    field_vec = 0.
    field0_vec = 0.
    jphi_field = 0.
    !vor_field = 0.
    !com_field = 0.
    !resistivity_field = 0.
    !kappa_field = 0.
    !visc_field = 0.
    !visc_c_field = 0.
    if(ipforce.gt.0) pforce_field = 0.
    if(ipforce.gt.0) pmach_field = 0.
    if(momentum_source) Fphi_field = 0.
    if(heat_source) Q_field = 0.
    if(rad_source) then
      Totrad_field = 0.
      Linerad_field = 0.
      Bremrad_field = 0.
      Ionrad_field = 0.
      Reckrad_field = 0.
      Recprad_field = 0.
    endif
    if(icd_source.gt.0) cd_field = 0.
    bf_field(0) = 0.
    bf_field(1) = 0.
    bfp_field(0) = 0.
    bfp_field(1) = 0.
    if(ibootstrap.gt.0) visc_e_field = 0.
    if(ibootstrap.gt.0) Jbs_L31_field = 0.
    if(ibootstrap.gt.0) Jbs_L32_field = 0.
    if(ibootstrap.gt.0) Jbs_L34_field = 0.
    if(ibootstrap.gt.0) Jbs_alpha_field = 0.
    if(ibootstrap.gt.0) Jbs_fluxavg_iBsq_field = 0.
    if(ibootstrap.eq.2 .or. ibootstrap.eq.3 ) Jbs_dtedpsit_field = 0.
    if(ibootstrap.eq.3) Jbs_ftrap_field = 0.
    if(ibootstrap.eq.3) Jbs_qR_field = 0.
    if(ibootstrap.eq.3) Jbs_invAspectRatio_field = 0.
    psi_coil_field = 0.
    !call destroy_auxiliary_fields
    !call create_auxiliary_fields

    call initial_conditions
    ! combine the equilibrium and perturbed fields of linear=0
    ! unless eqsubtract = 1
    if(eqsubtract.eq.0) then
      call add(field_vec, field0_vec)
      field0_vec = 0.
    endif
    i_control%err_i = 0.
    i_control%err_p_old = 0.
    n_control%err_i = 0.
    n_control%err_p_old = 0.
    call reset_scalars
    if(myrank.eq.0 .and. iprint.ge.2) print *, "  transport coefficients"
    call define_transport_coefficients
    if(eqsubtract.eq.1) then
      call derived_quantities(0)
    end if
    call derived_quantities(1)
    meshAdapted =1
  end subroutine adapt_by_spr

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine adapt_by_error
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    use diagnostics
    use basic
    use error_estimate
    use scorec_mesh_mod
    use basic
    use mesh_mod
    use arrays
    use newvar_mod
    use sparse
    use time_step
    use auxiliary_fields
    use scorec_mesh_mod
    use transport_coefficients

    vectype, allocatable :: edge_error(:,:)
    vectype, allocatable :: elm_error(:,:), elm_error_res(:,:), elm_error_sum(:,:)
    real, allocatable :: node_error(:,:)
    integer :: num_adj, elem_dim, num_edge, ii, jj, num_get, num_elm, num_node, ier
    integer, dimension(256) :: elms 

    integer, parameter :: vecsize=1 
    integer :: num_node_total
    real, dimension(2):: max_error, buff
    character(len=32) :: mesh_file_name, file_name1, file_name2, file_name3

    vectype, dimension(dofs_per_node*num_fields) :: max_val, min_val
    vectype :: maxPhi, maxPs
    vectype, dimension(NUMTERM) :: jump_sum

    write(mesh_file_name,"(A7,A)") 'adapted', 0
    write(file_name1, "(A9,I0,A)") 'errorJump', ntime,0
    write(file_name2, "(A8,I0,A)") 'errorElm', ntime,0
    write(file_name3,"(A8,I0,A)") 'errorSum', ntime,0

    call m3dc1_mesh_getnumglobalent (0, num_node_total)

    !if (myrank .eq. 0) print*, "time", ntime, "current max", max_val(1), min_val(1), "mesh size before adapt", num_node_total
    call m3dc1_field_max(field_vec%id, max_val, min_val)
    maxPhi = max(abs(max_val(1+(u_g-1)*dofs_per_node)),abs(min_val(1+(u_g-1)*dofs_per_node)))
    maxPs =  max(abs(max_val((psi_g-1)*dofs_per_node)+1),abs(min_val((psi_g-1)*dofs_per_node)+1))

    call m3dc1_mesh_getnument(1, num_edge)
    allocate(edge_error(num_edge, NUMTERM))
    elem_dim = 2
    num_adj=2
    if (nplanes .gt. 1) then
      elem_dim = 3
      num_adj = 8
    end if
 
    call m3dc1_mesh_getnument(elem_dim, num_elm)

    allocate(elm_error(num_elm, NUMTERM))
    allocate(elm_error_sum(2,num_elm))
    allocate(elm_error_res(2,num_elm))
    elm_error_res = 0
    call  m3dc1_mesh_getnument(0, num_node)
    allocate(node_error(2,num_node))
    node_error=0.
    edge_error=0.
    elm_error=0.
    call jump_discontinuity(edge_error)
    do ii=1, NUMTERM
       jump_sum(ii) =sum(edge_error(:,ii))
    end do

    do ii=1, num_edge
       call m3dc1_ent_getadj (1, ii-1, elem_dim, elms, num_adj, num_get)
       do jj=1, num_get 
          elm_error(elms(jj)+1,:)=elm_error(elms(jj)+1,:)+0.5*edge_error(ii,:);
       end do
    end do

    ! implementated for numvar .eq. 1  
    if (numvar .eq. 1 ) call elem_residule (elm_error_res(1,:), elm_error_res(2,:))

    elm_error_sum (1,:) = elm_error(:,JUMPU) + elm_error_res (1,:) 
    elm_error_sum (2,:) = elm_error(:,JUMPPSI) + elm_error_res(2,:)
    !if(isplitstep .eq. 1) elm_error_sum (:) = elm_error_sum (:) + elm_error(:,TSQ) 
    !call output_face_data (NUMTERM, sqrt(real(elm_error)), file_name1);
    !call output_face_data (2, TRANSPOSE(sqrt(real(elm_error_res))), file_name2);
    !call output_face_data (2, TRANSPOSE(sqrt(real(elm_error_sum))), file_name3);
    call get_node_error_from_elm (real(elm_error_sum(1,:)), 1, node_error(1,:));
    call get_node_error_from_elm (real(elm_error_sum(2,:)), 1, node_error(2,:));

    node_error = sqrt(node_error)
    !print *, edge_error
    deallocate(edge_error)
    deallocate(elm_error)
    deallocate(elm_error_sum)
    deallocate(elm_error_res)
    if (adapt_control .eq. 0) then
      max_error(1)= maxval(node_error(1,:))
      max_error(2)= maxval(node_error(2,:))
      buff = max_error
      call mpi_allreduce (buff, max_error, 2, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier )
    else
      max_error(1) = sum (node_error(1,:))
      max_error(2) = sum (node_error(2,:))
      !max_error= maxval(node_error(:,:))
      buff = max_error
      call mpi_allreduce (buff, max_error, 2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ier )
    end if

    if (max_error(1) .gt. error_tol * adapt_target_error .or. max_error(2) &
         .gt. error_tol * adapt_target_error) then

      if (myrank .eq. 0) then
        print *, "jump_sum:", jump_sum
        print *, "estimated error in engergy norm", max_error
        print *, "estimated error in solution in energy norm", solutionH2Norm
      end if

      !  if (myrank .eq. 0) print *, " error exceeds tolerance, start adapting mesh"
       call straighten_fields()
       abs_size(1) = adapt_hmin
       abs_size(2) = adapt_hmax
       rel_size(1) = adapt_hmin_rel 
       rel_size(2) = adapt_hmax_rel

       call set_mesh_size_bound (abs_size, rel_size)
       call set_adapt_p(iadapt_order_p)
       call adapt_by_error_field(sqrt(node_error(1,:)**2+node_error(2,:)**2), adapt_target_error, &
iadapt_max_node, adapt_control);
       if(iadapt_writevtk .eq. 1) call m3dc1_mesh_write (mesh_file_name,0,ntime)
       if(iadapt_writesmb .eq. 1) call m3dc1_mesh_write (mesh_file_name,1,ntime)
       call space(0)
       call update_nodes_owned()
       call tridef
       call unstraighten_fields()

       call create_newvar_matrices
       if(irestart .ne. 0 .or. linear .eq. 0) return
       if (myrank .eq. 0) print *, "reset simulation after adapt .."
       field_vec = 0.
       field0_vec = 0.
       jphi_field = 0.
       !resistivity_field = 0.
       !kappa_field = 0.
       !visc_field = 0.
       !visc_c_field = 0.
       if(ipforce.gt.0) pforce_field = 0.
       if(ipforce.gt.0) pmach_field = 0.
       if(momentum_source) Fphi_field = 0.
       if(heat_source) Q_field = 0.
       if(rad_source) then 
          Totrad_field = 0.
          Linerad_field = 0.
          Bremrad_field = 0.
          Ionrad_field = 0.
          Reckrad_field = 0.
          Recprad_field = 0.
       endif
       if(icd_source.gt.0) cd_field = 0.
       bf_field(0) = 0.
       bf_field(1) = 0.
       bfp_field(0) = 0.
       bfp_field(1) = 0.
       if(ibootstrap.gt.0) visc_e_field = 0.
       if(ibootstrap.gt.0) Jbs_L31_field = 0.
       if(ibootstrap.gt.0) Jbs_L32_field = 0.
       if(ibootstrap.gt.0) Jbs_L34_field = 0.
       if(ibootstrap.gt.0) Jbs_alpha_field = 0.
       if(ibootstrap.gt.0) Jbs_fluxavg_iBsq_field = 0.
       if(ibootstrap.eq.2 .or. ibootstrap.eq.3 ) Jbs_dtedpsit_field = 0.
       if(ibootstrap.eq.3) Jbs_ftrap_field = 0.
       if(ibootstrap.eq.3) Jbs_qR_field = 0.
       if(ibootstrap.eq.3) Jbs_invAspectRatio_field = 0.
       psi_coil_field = 0.
       !call destroy_auxiliary_fields
       !call create_auxiliary_fields

       call initial_conditions
       ! combine the equilibrium and perturbed fields of linear=0
       ! unless eqsubtract = 1
       if(eqsubtract.eq.0) then
         call add(field_vec, field0_vec)
         field0_vec = 0.
       endif
       i_control%err_i = 0.
       i_control%err_p_old = 0.
       n_control%err_i = 0.
       n_control%err_p_old = 0.
       call reset_scalars
       if(myrank.eq.0 .and. iprint.ge.2) print *, "  transport coefficients"
       call define_transport_coefficients
       if(eqsubtract.eq.1) then
         call derived_quantities(0)
       end if
       call derived_quantities(1)
       meshAdapted =1
    else
      if (myrank .eq. 0) print *, "error did not exceed tolerance -- adaptation not needed"
    end if
    deallocate(node_error)
  end subroutine adapt_by_error
 
  subroutine diagnose_adapt (flag)
    use diagnostics
    use basic
    integer :: flag
    !print *, "check diagnose_adapt: ke_previous, ekin, ntime",ke_previous,ekin,ntime
    flag = 0

    !if iadapt_ntime>0, run adapt_by_error at the end of every N time steps
    if (iadapt_ntime .gt. 0 .and. mod(ntime, iadapt_ntime) .eq. 0) then
      flag=1
      if (myrank .eq. 0) print *, "diagnose_adapt: run adaptation every", iadapt_ntime, "time step(s)"
    endif
    !if non-linear & iadapt_ntime=0, run adapt_by_error at the end of every time step
    if (linear .eq. 0 .and. iadapt_ntime .eq. 0) then
      flag=1
      if (myrank .eq. 0) print *, "diagnose_adapt: run adaptation every time step (non-linear)"
    endif
    !if linear, adapt_ke>0 & ekin>adapt_ke, run adapt_by_error at the end of time step
    if(linear .eq. 1 .and. adapt_ke .gt. 0 .and. ekin .gt. adapt_ke) then
      flag=1
      if (myrank .eq. 0) print *, "diagnose_adapt: run adaptation as ekin is greater than adapt_ke", adapt_ke, "(linear)"
    endif
  end subroutine diagnose_adapt
end module adapt
