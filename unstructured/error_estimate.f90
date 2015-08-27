module error_estimate
  use m3dc1_nint
#include "mpif.h"
  integer, parameter :: EOP_1 = 1
  integer, parameter :: EOP_DX = 2
  integer, parameter :: EOP_DY = 3
  integer, parameter :: EOP_DXX = 4
  integer, parameter :: EOP_DXY = 5
  integer, parameter :: EOP_DYY = 6
  integer, parameter :: EOP_DXXX = 7
  integer, parameter :: EOP_DXXY = 8
  integer, parameter :: EOP_DXYY = 9
  integer, parameter :: EOP_DYYY = 10
  integer, parameter :: EOP_DXXXX = 11
  integer, parameter :: EOP_DXXXY = 12
  integer, parameter :: EOP_DXXYY = 13
  integer, parameter :: EOP_DXYYY = 14
  integer, parameter :: EOP_DYYYY = 15
  integer, parameter :: EOP_NUM = 16
  real, private :: efterm(MAX_PTS, coeffs_per_element, EOP_NUM)
  vectype, dimension(MAX_PTS, EOP_NUM, dofs_per_element) :: emu79, enu79
  vectype, dimension(MAX_PTS, EOP_NUM) :: eph179, eph179_pre, eps179, eps179_pre, en179, evis79, eetar79 
  vectype, dimension(MAX_PTS, EOP_NUM) :: eph079, eps079, en079
  vectype,  dimension(MAX_PTS, EOP_NUM) :: U_bar, U_2bar, U_tild, psi_bar, psi_2bar, psi_tmp,psi_tmp2, psi_tmp3, U_tmp, n_tmp
  vectype,  dimension(MAX_PTS, EOP_NUM) :: U_delta, psi_delta
  vectype, dimension(MAX_PTS) :: er_79, er2_79, er4_79 
  real, private :: max_current, min_current
  integer, parameter :: LPUN = 1 !  mu * partial (laplace U)/ parital n
  integer, parameter :: LPPSIN =2 ! eta * partial (laplace Psi)/ parital n
  integer, parameter :: LPUT = 3 ! ro * laplace U * partial U / partial t
  integer, parameter :: LPU = 4 !  mu * laplace U; 
  integer, parameter :: LPPSI = 5 ! eta * laplace psi;
  integer, parameter :: LPPSPST = 6 ! laplace psi * partial psi  / partial t
  integer, parameter :: CPSIU = 7 ! partial n <psi, U>
  integer, parameter :: R2U = 8
  integer, parameter :: R2PSI = 9
  integer, parameter :: TEST = 10 
  integer, parameter :: JUMPU = 11
  integer, parameter :: JUMPPSI = 12
  integer, parameter :: NUMTERM = 12
  integer, parameter :: npoint_int =3
  integer, parameter :: npoint_int_elm = 12
  vectype, allocatable :: edge_jump (:,:,:)
  vectype, dimension (MAX_PTS) :: fn_eval, fn_eval_tmp, fn_eval_tmp2, fn_eval_tmp3, fn_eval_tmp4, fn_eval_tmp5, fn_eval_tmp6
  vectype, dimension (MAX_PTS) :: fn_eval_elm
  vectype, dimension (MAX_PTS) :: fn_eval_elm_tmp

  real, dimension(2) :: solutionH2Norm
  integer :: iadapt_useH1, iadapt_removeEquiv
  contains
  subroutine getElmL2 (vec, value)
    use m3dc1_nint
    implicit none
    vectype, dimension (npoint_int_elm) :: vec
    real :: value
#ifdef USECOMPLEX
    value = int1(vec*CONJG(vec))
#else
    value = int1(vec**2)
#endif 
  end subroutine getElmL2 
  subroutine getElmH1 (vec, value)
    use m3dc1_nint
    implicit none
    vectype, dimension (npoint_int_elm,EOP_NUM) :: vec
    real :: value, buff

    call getElmL2(vec(:,EOP_DX), buff)
    value = buff
    call getElmL2(vec(:,EOP_DY), buff)
    value = value + buff 
  end subroutine getElmH1

  subroutine getElmH2 (vec, value)
    use m3dc1_nint
    implicit none
    vectype, dimension (npoint_int_elm, EOP_NUM) :: vec
    real :: value, buff
    call getElmL2(vec(:,EOP_DXX), buff)
    value = buff
    call getElmL2(vec(:,EOP_DYY), buff)
    value = value + buff
    call getElmL2(vec(:,EOP_DXY), buff)
  end subroutine getElmH2

  subroutine jump_discontinuity (edge_error)
    use basic
    use arrays
    use scorec_mesh_mod
    use vector_mod
    use m3dc1_nint
    implicit none
    integer :: itri, numelms,numedgs, iedge, ii, node_next, num_get
    vectype, dimension(:,:), intent(out) :: edge_error
    integer, dimension(3) :: edges, nodes, nodes_edge, edge_dir, idimgeo, idimgeo_t, is_bdy
    real, dimension(3) :: edge_len
    real, dimension(2,3) :: normal, normal_t
    real, dimension(3) :: coords
    integer, allocatable :: edge_tag(:)
    type(element_data) :: d
    integer :: ifaczonedim, ifaczone, icylinder
    call m3dc1_mesh_getnument(1, numedgs)

    icylinder = itor
    allocate(edge_jump(npoint_int, numedgs, NUMTERM))
    edge_jump = 0
    !print*, "shape edge_jump", shape(edge_jump)
    idimgeo=2;
    numelms = local_elements()
    do itri=1,numelms 
       call m3dc1_ent_getgeomclass(2, itri-1, ifaczonedim, ifaczone)
       if(ifaczone .ne. 1) continue;
       call get_element_data(itri, d)
       call m3dc1_ent_getadj (2, itri-1, 1, edges, 3, num_get)
       call m3dc1_ent_getadj (2, itri-1, 0, nodes, 3, num_get)
       do iedge=1, 3
          call get_edge_data (edges(iedge), normal(:,iedge), edge_len(iedge))
          call m3dc1_ent_getadj (1, edges(iedge), 0, nodes_edge, 2, num_get)
          !print *, "edge node", edges(iedge),nodes_edge
          edge_dir(iedge)=0
          do ii=1,3
             if (nodes(ii) .eq. nodes_edge(1)) then
                node_next=ii+1;
                if (node_next .gt. 3) node_next=1
                if (nodes(node_next) .eq. nodes_edge(2)) edge_dir(iedge)=1
             end if
          end do
          if (edge_dir(iedge) .eq. 0) then
             normal(:,iedge)=-1.*normal(:,iedge)
          end if
       end do
       call boundary_edge(itri, is_bdy, normal_t, idimgeo_t)
       !print *,myrank, 'edges',edges,'nodes',nodes
       !do ii=1,3
          !call m3dc1_node_getcoord(nodes(ii),coords)
          !print*, myrank, coords
       !end do
       !print *, myrank, 'edge_dir',edge_dir
       !print *, myrank, 'edge_le', edge_len
       !print *, myrank, 'normal'
       !print *, myrank, normal(:,1)
       !print *, myrank, normal(:,2)
       !print *, myrank, normal(:,3)

       do iedge=1,3
          if(is_bdy(iedge).ne.0) cycle
          do ii=1,3
             normal_t(:,ii)=normal(:,iedge)
          end do
          call define_boundary_quadrature(itri, iedge, npoint_int, npoint_int, normal_t, idimgeo)
          call define_fields_error(itri)

          !call rotate_vec(normal(:,iedge), d%sn, d%co )
          !print *, 'ele dir',d%co, d%sn
          !print *, 'rotate', normal(:,iedge)
          ! caculate jump discontinuity
          !if (isplitstep .eq. 1) then
          !   psi_tmp = eps179_pre
          !else
             psi_tmp = psi_2bar
          !end if
          if (linear .eq. 1) then
             psi_tmp2 = eps079
             U_tmp = eph079
          else
             psi_tmp2 = eps179_pre
             U_tmp = eph179_pre
          end if

          !lp_U_n
          call laplace_gs_r (U_bar, fn_eval_tmp, icylinder, 0)
          call laplace_gs_z (U_bar, fn_eval_tmp2, icylinder, 0)

          fn_eval=evis79(:,EOP_1)*(fn_eval_tmp*normal(1,iedge) + fn_eval_tmp2*normal(2,iedge))

          if(edge_dir(iedge) .eq. 0) call reverse_fn(fn_eval(1:npoint_int),npoint_int)
          edge_jump(1:npoint_int,edges(iedge)+1, LPUN)= edge_jump(1:npoint_int,edges(iedge)+1, LPUN)+ fn_eval(1:npoint_int)
          !LPPSIN
          if(jadv .eq. 1 ) then
             call laplace_gs_r (psi_bar, fn_eval_tmp, icylinder, 0)
             call laplace_gs_z (psi_bar, fn_eval_tmp2, icylinder, 0)
             fn_eval=eetar79(:,EOP_1)*(fn_eval_tmp*normal(1,iedge)+ fn_eval_tmp2*normal(2,iedge))
            if(edge_dir(iedge) .eq. 0) call reverse_fn(fn_eval(1:npoint_int),npoint_int)
            edge_jump(1:npoint_int,edges(iedge)+1, LPPSIN)= edge_jump(1:npoint_int,edges(iedge)+1, LPPSIN)+ fn_eval(1:npoint_int)
          end if
          !print *, edges(iedge), 'dir',edge_dir(iedge)
          !print *, fn_eval(1:5)

          !rho lp_U * U_t
          call laplace_gs (U_2bar, fn_eval_tmp, icylinder, 0)
          fn_eval=en179(:,EOP_1)*fn_eval_tmp*(-U_tmp(:,EOP_DX)*normal(2,iedge) + U_tmp(:,EOP_DY)*normal(1,iedge))
          call laplace_gs (U_tmp, fn_eval_tmp, icylinder, 0)
          fn_eval=fn_eval+en179(:,EOP_1)*fn_eval_tmp*(-U_2bar(:,EOP_DX)*normal(2,iedge)+ U_2bar(:,EOP_DY)*normal(1,iedge))

          if(edge_dir(iedge) .eq. 0) call reverse_fn(fn_eval(1:npoint_int),npoint_int)
          edge_jump(1:npoint_int,edges(iedge)+1, LPUT)=edge_jump(1:npoint_int,edges(iedge)+1,LPUT)+0.5*fn_eval(1:npoint_int)

          !lp_psi * Psi_t
          call laplace_gs (psi_tmp, fn_eval_tmp, icylinder, 0)
          call laplace_gs (psi_tmp2, fn_eval_tmp2, icylinder, 0)
          fn_eval = fn_eval_tmp*(-psi_tmp2(:,EOP_DX)*normal(2,iedge)+psi_tmp2(:,EOP_DY)*normal(1,iedge))
          fn_eval = fn_eval + fn_eval_tmp2*(-psi_tmp(:,EOP_DX)*normal(2,iedge)+psi_tmp(:,EOP_DY)*normal(1,iedge))
          if(edge_dir(iedge) .eq. 0) call reverse_fn(fn_eval(1:npoint_int),npoint_int)
          edge_jump(1:npoint_int,edges(iedge)+1, LPPSPST)=edge_jump(1:npoint_int,edges(iedge)+1,LPPSPST)+0.5*fn_eval(1:npoint_int)

          !lp_U
          call laplace_gs (U_bar, fn_eval_tmp, icylinder, 0)
          fn_eval = evis79(:, EOP_1)*fn_eval_tmp
          if (edge_dir(iedge) .eq. 0) then
             call reverse_fn(fn_eval(1:npoint_int),npoint_int)
             fn_eval=-fn_eval;
          end if
          edge_jump(1:npoint_int,edges(iedge)+1, LPU)=edge_jump(1:npoint_int,edges(iedge)+1, LPU)+fn_eval(1:npoint_int)
          !lp_psi
          if(jadv .eq. 1 ) then
            call laplace_gs (psi_bar, fn_eval_tmp, icylinder, 0)
            fn_eval = eetar79(:,EOP_1)*fn_eval_tmp
            if (edge_dir(iedge) .eq. 0) then
               call reverse_fn(fn_eval(1:npoint_int),npoint_int)
               fn_eval=-fn_eval;
            end if
            edge_jump(1:npoint_int,edges(iedge)+1, LPPSI)=edge_jump(1:npoint_int,edges(iedge)+1, LPPSI)+fn_eval(1:npoint_int)
          end if

          !test jump field should be zero
          ! test field dR*dZ
          fn_eval = eps179(:,EOP_DX) * eps179(:,EOP_DY) 
          if (edge_dir(iedge) .eq. 0) then
             call reverse_fn(fn_eval(1:npoint_int),npoint_int)
             fn_eval=-fn_eval;
          end if
          edge_jump(1:npoint_int,edges(iedge)+1, TEST)=edge_jump(1:npoint_int,edges(iedge)+1, TEST)+fn_eval(1:npoint_int)
          !edge_jump(1:npoint_int,edges(iedge)+1, TEST)=eps179(1:npoint_int,EOP_DXX)
          !print *, 'test', edges(iedge)+1, fn_eval(1:npoint_int)
          ! include delta t terms for isplitstep =1
          !if (isplitstep .eq. 1) then
          !   fn_eval =(psi_tmp2(:,EOP_DX)*normal(2,iedge) &
           !           - psi_tmp2(:,EOP_DY)*normal(1,iedge));
           !  psi_tmp3 = psi_tmp2
          !else 
             !if(jadv .eq. 1) then
               !fn_eval(1:npoint_int)=normal(2,iedge)-normal(1,iedge);
               !psi_tmp3 = psi_bar
             !end if
          !end if
          !call laplace_gs_r (U_tild, fn_eval_tmp, icylinder, 0)
          !call laplace_gs_z (U_tild, fn_eval_tmp2, icylinder, 0)
          !call laplace_gs_r (psi_tmp3, fn_eval_tmp3, icylinder, icylinder)
          !call laplace_gs_z (psi_tmp3, fn_eval_tmp4, icylinder, icylinder)
          !fn_eval = fn_eval &
          !            *( fn_eval_tmp * psi_tmp3(:,EOP_DY) &
          !               -fn_eval_tmp2* psi_tmp3(:,EOP_DX) &
          !               + fn_eval_tmp3 * U_tild(:,EOP_DY) &
          !               - fn_eval_tmp4 * U_tild(:,EOP_DX) & 
          !               + 2* U_tild(:,EOP_DXX)*psi_tmp3(:,EOP_DXY) &
          !               - 2 * U_tild(:,EOP_DXY)*psi_tmp3(:,EOP_DYY) &
          !               + 2* U_tild(:,EOP_DXY)*psi_tmp3(:,EOP_DYY) &
          !               - 2 * U_tild(:,EOP_DYY)*psi_tmp3(:,EOP_DXY) )

          fn_eval = (psi_2bar(:, EOP_DXX)*U_tmp(:, EOP_DY)+ psi_2bar(:, EOP_DX)*U_tmp(:, EOP_DXY) &
                     -psi_2bar(:, EOP_DXY)*U_tmp(:, EOP_DX)-psi_2bar(:, EOP_DY)*U_tmp(:, EOP_DXX))*normal(1,iedge) &
                   +(psi_2bar(:, EOP_DXY)*U_tmp(:, EOP_DY)+psi_2bar(:, EOP_DX)*U_tmp(:, EOP_DYY) &
                     -psi_2bar(:, EOP_DYY)*U_tmp(:, EOP_DX)-psi_2bar(:, EOP_DY)*U_tmp(:, EOP_DXY))*normal(2,iedge)
          fn_eval = fn_eval + (psi_tmp2(:, EOP_DXX)*U_2bar(:, EOP_DY)+ psi_tmp(:, EOP_DX)*U_2bar(:, EOP_DXY) &
                     -psi_tmp(:, EOP_DXY)*U_2bar(:, EOP_DX)-psi_tmp(:, EOP_DY)*U_2bar(:, EOP_DXX))*normal(1,iedge) &
                   +(psi_tmp(:, EOP_DXY)*U_2bar(:, EOP_DY)+psi_tmp(:, EOP_DX)*U_2bar(:, EOP_DYY) &
                     -psi_tmp(:, EOP_DYY)*U_2bar(:, EOP_DX)-psi_tmp(:, EOP_DY)*U_2bar(:, EOP_DXY))*normal(2,iedge)
          if (edge_dir(iedge) .eq. 0) then
             call reverse_fn(fn_eval(1:npoint_int),npoint_int)
          end if
          !if(isplitstep .eq. 1) then
            edge_jump(1:npoint_int,edges(iedge)+1, CPSIU) = edge_jump(1:npoint_int,edges(iedge)+1, CPSIU) + 0.5*fn_eval(1:npoint_int) 
          !else
            !edge_jump(1:npoint_int,edges(iedge)+1, TSQ) = edge_jump(1:npoint_int,edges(iedge)+1, TSQ) + dt*fn_eval(1:npoint_int)

          !end if
          edge_jump(1:npoint_int,edges(iedge)+1,R2U) =  edge_jump(1:npoint_int,edges(iedge)+1, LPUN) + edge_jump(1:npoint_int,edges(iedge)+1, LPUT) - edge_jump(1:npoint_int,edges(iedge)+1, LPPSPST) 
          edge_jump(1:npoint_int,edges(iedge)+1,R2PSI) =  edge_jump(1:npoint_int,edges(iedge)+1, CPSIU) + edge_jump(1:npoint_int,edges(iedge)+1, LPPSIN)

          !if(isplitstep .eq. 1) then
            !edge_jump(1:npoint_int,edges(iedge)+1,R1) = edge_jump(1:npoint_int,edges(iedge)+1,R1) + edge_jump(1:npoint_int,edges(iedge)+1, TSQ)
          !end if 
          ! R3
          !if(isplitstep .eq. 1 .or. isplitstep .eq. 0 .and. jadv .eq. 1) then
           ! call laplace_gs (U_tild, fn_eval_tmp, icylinder, 0)
            !call laplace_gs (psi_tmp2, fn_eval_tmp2, icylinder, icylinder)
            !fn_eval = fn_eval_tmp * (psi_tmp2(:,EOP_DX)*normal(2,iedge)-psi_tmp2(:,EOP_DY)*normal(1,iedge)) + fn_eval_tmp2 * (U_tild(:,EOP_DX)*normal(2,iedge)-U_tild(:,EOP_DY)*normal(1,iedge)) 
            !if (edge_dir(iedge) .eq. 0) then
             !  call reverse_fn(fn_eval(1:npoint_int),npoint_int)
            !end if
            !if(isplitstep .eq. 1) then
             ! edge_jump(1:npoint_int,edges(iedge)+1,R3) = edge_jump(1:npoint_int,edges(iedge)+1,R3) + dt*dt*fn_eval(1:npoint_int)
            !else
              !edge_jump(1:npoint_int,edges(iedge)+1,R3) = edge_jump(1:npoint_int,edges(iedge)+1,R3) + dt*fn_eval(1:npoint_int)
            !end if
          !end if
       end do
    end do
    do ii=1,NUMTERM
       !if(ii .eq. 5) print*,myrank, edge_jump(1:npoint_int,:,ii)
#ifdef USECOMPLEX
       call sum_edge_data(edge_jump(1:npoint_int,:,ii),npoint_int*2)
#else
       call sum_edge_data(edge_jump(1:npoint_int,:,ii),npoint_int)
#endif
    end do
    ! integrate over the edge
    allocate(edge_tag(numedgs));
    edge_tag=0
    do itri=1, numelms
       call m3dc1_ent_getadj (2, itri-1, 1, edges, 3, num_get)
       do iedge=1,3
          if(edge_tag(edges(iedge)+1).ne.0) cycle
          edge_tag(edges(iedge)+1)=1
          call define_boundary_quadrature(itri, iedge, npoint_int, npoint_int, normal_t, idimgeo)
          call get_edge_data (edges(iedge), normal(:,1), edge_len(1))
          do ii=1, NUMTERM
#ifdef USECOMPLEX
             edge_error(edges(iedge)+1,ii)=int1(edge_jump(1:npoint_int,edges(iedge)+1,ii)*CONJG(edge_jump(1:npoint_int,edges(iedge)+1,ii)))
             !edge_error(edges(iedge)+1,ii)=real(int1(edge_jump(1:npoint_int,edges(iedge)+1,ii)))**2
#else
             edge_error(edges(iedge)+1,ii)=int1(edge_jump(1:npoint_int,edges(iedge)+1,ii)**2)
#endif
          end do
          edge_error(edges(iedge)+1,LPUN)=edge_error(edges(iedge)+1,LPUN)* edge_len(1)**3 
          edge_error(edges(iedge)+1,LPPSIN)=edge_error(edges(iedge)+1,LPPSIN)* edge_len(1)**3
          edge_error(edges(iedge)+1,LPUT)=edge_error(edges(iedge)+1,LPUT)* edge_len(1)**3 
          edge_error(edges(iedge)+1,LPPSPST)=edge_error(edges(iedge)+1,LPPSPST)* edge_len(1)**3 
          edge_error(edges(iedge)+1,LPU)=edge_error(edges(iedge)+1,LPU)*edge_len(1)
          edge_error(edges(iedge)+1,LPPSI)=edge_error(edges(iedge)+1,LPPSI)*edge_len(1)
          edge_error(edges(iedge)+1,CPSIU)=edge_error(edges(iedge)+1,CPSIU)* edge_len(1)**3 
          edge_error(edges(iedge)+1,R2U)=edge_error(edges(iedge)+1,R2U)*edge_len(1)**3 
          edge_error(edges(iedge)+1,R2PSI)=edge_error(edges(iedge)+1,R2PSI)*edge_len(1)**3
          edge_error(edges(iedge)+1,JUMPU)=edge_error(edges(iedge)+1,R2U)+edge_error(edges(iedge)+1,LPU)
          edge_error(edges(iedge)+1,JUMPPSI)=edge_error(edges(iedge)+1,R2PSI)+edge_error(edges(iedge)+1,LPPSI)
          if(iadapt_useH1 .eq. 1) edge_error(edges(iedge)+1,:) = edge_error(edges(iedge)+1,:) /edge_len(1)**2;
       end do
    end do
    deallocate (edge_tag)
    deallocate(edge_jump)
  end subroutine jump_discontinuity

  subroutine reverse_fn( f, n)
    implicit none
    integer :: n
    vectype, dimension(n) :: f
    vectype :: buff
    integer :: ii
    do ii=1, n/2
       buff=f(ii)
       f(ii)=f(n+1-ii)
       f(n+1-ii)=buff
    end do
  end subroutine reverse_fn

  subroutine get_edge_data (edge, normal, edge_len)
     use scorec_mesh_mod
     implicit none
     integer :: edge
     real, dimension(2) :: normal
     integer, dimension(2) :: nodes
     real :: edge_len
     integer :: num_get, inode
     real, dimension(3,2) :: node_pos 
     call m3dc1_ent_getadj (1, edge, 0, nodes, 2, num_get)
     do inode=1,2
        call get_node_pos(nodes(inode)+1,node_pos(1,inode),node_pos(3,inode),node_pos(2,inode))
     end do
     normal(1)=node_pos(2,2)-node_pos(2,1)
     normal(2)=node_pos(1,1)-node_pos(1,2)
     edge_len=sqrt( normal(1)**2+normal(2)**2)
     normal=normal/edge_len
   end subroutine

  subroutine define_fields_error (itri)
    use basic
    use mesh_mod
    use arrays

    implicit none

    integer, intent(in) :: itri

    integer :: i
    type(element_data) :: d
    vectype, dimension(dofs_per_element,coeffs_per_element) :: cl

    ! calculate the major radius, and useful powers
    call get_element_data(itri, d)
    call local_to_global(d, xi_79, zi_79, eta_79, x_79, phi_79, z_79)
    if(itor.eq.1) then
       er_79 = x_79
    else
       er_79 = 1.
    endif
    er2_79 = er_79**2
    er4_79 = er2_79**2

    !if(ijacobian.eq.1) weight_79 = weight_79 * r_79

    call precalculate_terms_error(xi_79,eta_79,npoints)
    call define_basis_error(itri, d%co, d%sn)

    call eval_ops_error(itri, u_field(1), eph179)
    call eval_ops_error(itri, psi_field(1), eps179)
    call eval_ops_error(itri, den_field(1), en179)

    call eval_ops_error(itri, u_field_pre, eph179_pre)
    call eval_ops_error(itri, psi_field_pre, eps179_pre)

    call eval_ops_error(itri, u_field(0), eph079)
    call eval_ops_error(itri, psi_field(0), eps079)
    call eval_ops_error(itri, den_field(0), en079)

    call eval_ops_error(itri, visc_field, evis79)
    call eval_ops_error(itri, resistivity_field, eetar79)


    !if(eqsubtract.eq.1) then
    !   ph179 = ph079 + ph179
    !else
    !   ph079 = 0.
    !endif

    !if(eqsubtract.eq.1) then
    !   ps179 = ps179 + ps079
    !else
    !   ps079 = 0.
    !endif

    if(eqsubtract.eq.1) then
      en179 = en079;
    end if
    if(iadapt_removeEquiv .eq. 1) then
      eph079 =0.
      eps079 =0.
    end if
    if(ivform .eq. 1 .and. itor .eq. 1) then
      do i=1, EOP_NUM
         eph179(:,i)=eph179(:,i)*er_79**2
         eph179_pre(:,i)=eph179_pre(:,i)*er_79**2
         eph079(:,i)=eph079(:,i)*er_79**2
         !eps179(:,i)=eps179(:,i)/er_79
         !eps179_pre(:,i)=eps179(:,i)/er_79
      end do
    end if   
    !U_bar = thimp * eph179 +(1. -thimp) * eph179_pre
    !U_tild = thimp*thimp*eph179 + thimp*(1-thimp)*eph179_pre
    !U_2bar = 2*thimp* eph179 +(1. -2*thimp) * eph179_pre
    !psi_bar = thimp * eph179 +(1. -thimp) * eph179_pre
    !psi_2bar = 2*thimp* eps179 +(1. -2*thimp) * eps179_pre
    U_bar = eph179
    U_tild = eph179
    U_2bar = eph179
    psi_bar = eps179
    psi_2bar = eps179
    U_delta = eph179 - eph179_pre
    psi_delta = eps179 - eps179_pre
  end subroutine define_fields_error

  subroutine eval_ops_error(itri,fin,outarr)
    use field
    implicit none

    integer, intent(in) :: itri
    type(field_type), intent(in) :: fin
    vectype, dimension(MAX_PTS, EOP_NUM), intent(out) :: outarr

    integer :: i, op
    vectype, dimension(dofs_per_element) :: dofs

    call get_element_dofs(fin, itri, dofs)

    outarr = 0.

    do op=1, EOP_NUM
       do i=1, dofs_per_element
          outarr(:,op) = outarr(:,op) + dofs(i)*enu79(:,op,i)
       end do
    end do
  end subroutine eval_ops_error

  subroutine precalculate_terms_error(xi,eta,npoints)
    use basic

    implicit none
      
    integer, intent(in) :: npoints
    real, dimension(MAX_PTS), intent(in) :: xi, eta 
    real, dimension(MAX_PTS) :: temp

    integer :: p
    real :: xpow(MAX_PTS,-3:5), ypow(MAX_PTS,-3:5)

    xpow(:,-3:-1) = 0.
    ypow(:,-3:-1) = 0.
    xpow(:,0) = 1.
    ypow(:,0) = 1.

    do p=1, 5
       xpow(:,p) = xpow(:,p-1)*xi(:)
       ypow(:,p) = ypow(:,p-1)*eta(:)
    end do
    
    efterm = 0.
    do p=1, coeffs_per_tri
       efterm(:,p,EOP_1) = xpow(:,mi(p))*ypow(:,ni(p))
       ! first order    
       if(mi(p).ge.1) then
          efterm(:,p,EOP_DX)  = mi(p)*xpow(:,mi(p)-1) * ypow(:,ni(p))
       end if
       if(ni(p).ge.1) then
          efterm(:,p,EOP_DY)  = ni(p)*ypow(:,ni(p)-1) * xpow(:,mi(p))
       end if
       ! second order
       if(mi(p).ge.2) then
          efterm(:,p,EOP_DXX) = (mi(p)-1)*mi(p)*xpow(:,mi(p)-2) * ypow(:,ni(p))
       end if
       if(ni(p).ge.2) then
          efterm(:,p,EOP_DYY) = (ni(p)-1)*ni(p)*ypow(:,ni(p)-2) * xpow(:,mi(p))
       end if
       if(mi(p).ge.1 .and. ni(p).ge.1) then
          efterm(:,p,EOP_DXY) = mi(p)*ni(p)*xpow(:,mi(p)-1) * ypow(:,ni(p)-1)
       end if
       ! third order
       if(mi(p).ge.3) then
          efterm(:,p,EOP_DXXX) = (mi(p)-2)*(mi(p)-1)*mi(p)*xpow(:,mi(p)-3) * ypow(:,ni(p))
       end if
       if(mi(p).ge.2 .and. ni(p) .ge. 1) then
          efterm(:,p,EOP_DXXY) = (mi(p)-1)*mi(p)*ni(p)*xpow(:,mi(p)-2) * ypow(:,ni(p)-1)
       end if
       if(mi(p).ge.1 .and. ni(p) .ge. 2) then
          efterm(:,p,EOP_DXYY) = mi(p)*(ni(p)-1)*ni(p)*xpow(:,mi(p)-1) * ypow(:,ni(p)-2)
       end if
       if(ni(p) .ge. 3) then
          efterm(:,p,EOP_DYYY) = (ni(p)-2)*(ni(p)-1)*ni(p)*xpow(:,mi(p)) * ypow(:,ni(p)-3)
       end if
       ! forth order
       if(mi(p).ge.4) then
          efterm(:,p,EOP_DXXXX) =(mi(p)-3)*(mi(p)-2)*(mi(p)-1)*mi(p)*xpow(:,mi(p)-4) * ypow(:,ni(p))
       end if
       if(mi(p).ge.3 .and. ni(p) .ge. 1) then
          efterm(:,p,EOP_DXXXY) =(mi(p)-2)*(mi(p)-1)*mi(p)*ni(p)*xpow(:,mi(p)-3) * ypow(:,ni(p)-1)
       end if
       if(mi(p).ge.2 .and. ni(p) .ge. 2) then
          efterm(:,p,EOP_DXXYY) =(mi(p)-1)*mi(p)*(ni(p)-1)*ni(p)*xpow(:,mi(p)-2) * ypow(:,ni(p)-2)
       end if
       if(mi(p).ge.1 .and. ni(p) .ge. 3) then
          efterm(:,p,EOP_DXYYY) =mi(p)*(ni(p)-2)*(ni(p)-1)*ni(p)*xpow(:,mi(p)-1) * ypow(:,ni(p)-3)
       end if
       if(ni(p) .ge. 4) then
          efterm(:,p,EOP_DYYYY) =(ni(p)-3)*(ni(p)-2)*(ni(p)-1)*ni(p)*xpow(:,mi(p)) * ypow(:,ni(p)-4)
       end if
    end do
  end subroutine precalculate_terms_error

  subroutine define_basis_error(itri, c, s)
    use basic
    implicit none

    integer, intent(in) :: itri
    real, dimension(dofs_per_element,coeffs_per_element) :: cl
    real :: c,s

    integer :: i, p, op

    emu79 = 0.
    call local_coeff_vector(itri, cl)
    do op=1, EOP_NUM
       do i=1, dofs_per_element
          do p=1, coeffs_per_element
             emu79(:,op,i) = emu79(:,op,i) + efterm(:,p,op)*cl(i,p)
          end do
       end do
    end do
    do i=1,dofs_per_element
       call rotate_field(emu79(:,:,i),c,s)
    end do
    enu79 = emu79
  end subroutine define_basis_error

  subroutine rotate_vec ( vec, sn, co)
    implicit none
    real, dimension(2) :: vec, vec_buff
    real :: sn, co
    vec_buff(1) =  vec(1)*co + vec(2)*sn
    vec_buff(2) =  -vec(1)*sn + vec(2)*co
    vec = vec_buff
  end subroutine rotate_vec

  subroutine elem_residule (elm_res_U, elm_res_psi)
    use basic
    use arrays
    use scorec_mesh_mod
    use vector_mod
    use m3dc1_nint
    implicit none

    integer :: itri, numelms, ier
    vectype, dimension(:), intent(out) :: elm_res_U, elm_res_psi
    type(element_data) :: d
    real :: buff, norm_buff
    real, dimension(2) :: buff2Norm
    real, dimension(npoint_int_elm) :: buff2
    integer :: ifaczonedim, ifaczone, icylinder

    icylinder = itor

    numelms = local_elements()
    max_current = 0
    min_current = 0
    do itri=1,numelms
       call m3dc1_ent_getgeomclass(2, itri-1, ifaczonedim, ifaczone)
       if(ifaczone .ne. 1) continue;
       call get_element_data(itri, d)
       call define_element_quadrature(itri, npoint_int_elm, 0)
       call define_fields_error(itri)
       call getElmH1 (eps179,norm_buff)
       solutionH2Norm(2) = solutionH2Norm(2)+norm_buff/dt
       call getElmH2 (eps179,norm_buff)
       solutionH2Norm(2) = solutionH2Norm(2)+etar*norm_buff
       call getElmH1 (ph179,norm_buff)
       solutionH2Norm(1) = solutionH2Norm(1)+norm_buff/dt
       call getElmH2 (ph179,norm_buff)
       solutionH2Norm(1) = solutionH2Norm(1)+amu*norm_buff

       buff2=eps179(1:npoint_int_elm, EOP_DXX)+eps179(1:npoint_int_elm, EOP_DYY)
       buff = maxval(buff2)
       if(buff>max_current) max_current = buff 
       buff = minval(buff2)
       if(buff<min_current) min_current = buff
       fn_eval_elm = 0
       fn_eval_elm_tmp = 0
       call biharmic_gs_op (U_bar, fn_eval_tmp, icylinder, 0)
       fn_eval_elm_tmp  = fn_eval_tmp
       fn_eval_elm_tmp (1:npoint_int_elm) = evis79(1:npoint_int_elm, EOP_1)*fn_eval_elm_tmp (1:npoint_int_elm)
       fn_eval_elm = fn_eval_elm - fn_eval_elm_tmp

       if (isplitstep .eq. 1 .and. .false.) then ! skip
         !<psi, <psi, laplace U > >
         if (linear .eq. 0) then
            psi_tmp = eps179_pre
         else
            psi_tmp = eps079
         end if
         call laplace_gs_r (U_tild, fn_eval_tmp, icylinder, 0)
         call laplace_gs_z (U_tild, fn_eval_tmp2, icylinder, 0)
         call laplace_gs_rr (U_tild, fn_eval_tmp3, icylinder, 0)
         call laplace_gs_zz (U_tild, fn_eval_tmp4, icylinder, 0)
         call laplace_gs_rz (U_tild, fn_eval_tmp5, icylinder, 0)
         fn_eval_elm_tmp = psi_tmp(:, EOP_DX) &
                           *( psi_tmp(:, EOP_DXY) * fn_eval_tmp2  &
                           + psi_tmp(:, EOP_DX) *  fn_eval_tmp4 &
                           -psi_tmp(:, EOP_DYY) * fn_eval_tmp & 
                           +psi_tmp(:, EOP_DY) * fn_eval_tmp5 ) &
                           -psi_tmp(:, EOP_DY) & 
                           *(psi_tmp(:, EOP_DXX) * fn_eval_tmp2 & 
                           +eps179(:, EOP_DX) * fn_eval_tmp5 &
                           -psi_tmp(:, EOP_DXY) * fn_eval_tmp & 
                           +psi_tmp(:, EOP_DY) * fn_eval_tmp3 )
         fn_eval_elm = fn_eval_elm + dt*dt * fn_eval_elm_tmp
         ! <psi, <laplace psi, U> >
         call laplace_gs_r (psi_tmp, fn_eval_tmp, icylinder, icylinder)
         call laplace_gs_z (psi_tmp, fn_eval_tmp2, icylinder, icylinder)
         call laplace_gs_rr (psi_tmp, fn_eval_tmp3, icylinder, icylinder)
         call laplace_gs_zz (psi_tmp, fn_eval_tmp4, icylinder, icylinder)
         call laplace_gs_rz (psi_tmp, fn_eval_tmp5, icylinder, icylinder)
         fn_eval_elm_tmp  = -psi_tmp(:, EOP_DX) &
                            *( ph179(:, EOP_DXY) * fn_eval_tmp2 &
                            +U_tild(:, EOP_DX) * fn_eval_tmp4  &
                            -U_tild(:, EOP_DYY) * fn_eval_tmp &
                            +U_tild(:, EOP_DY) * fn_eval_tmp5 ) &
                            +psi_tmp(:, EOP_DY) &
                            *( U_tild(:, EOP_DXX) * fn_eval_tmp2 &
                            +U_tild(:, EOP_DX) * fn_eval_tmp5 &
                            -U_tild(:, EOP_DXY) * fn_eval_tmp &
                            +U_tild(:, EOP_DY) * fn_eval_tmp3 )
         fn_eval_elm = fn_eval_elm + dt*dt * fn_eval_elm_tmp
         ! <psi, <psi_R, u_R> >
         fn_eval_elm_tmp (:) = psi_tmp(:, EOP_DX) &
                               *( psi_tmp(:, EOP_DXXY) * U_tild(:, EOP_DXY) &
                               +psi_tmp(:, EOP_DXX) * U_tild(:, EOP_DXYY) &
                               +psi_tmp(:, EOP_DXYY) * U_tild(:, EOP_DXX) &
                               +psi_tmp(:, EOP_DXY) * U_tild(:, EOP_DXXY) ) &
                               -psi_tmp(:, EOP_DY) &
                               *( psi_tmp(:, EOP_DXXX) * U_tild(:, EOP_DXY) &
                               +psi_tmp(:, EOP_DXX) * U_tild(:, EOP_DXXY) &
                               +psi_tmp(:, EOP_DXXY) * U_tild(:, EOP_DXX) &
                               +psi_tmp(:, EOP_DXY) * U_tild(:, EOP_DXXX) )
         fn_eval_elm = fn_eval_elm + 2 * dt*dt * fn_eval_elm_tmp

         ! <psi, <psi_Z, u_Z> >
         fn_eval_elm_tmp (:) = psi_tmp(:, EOP_DX) &
                               *( psi_tmp(:, EOP_DXYY) * U_tild(:, EOP_DXY) &
                               + psi_tmp(:, EOP_DXY) * U_tild(:, EOP_DXYY) &
                               + psi_tmp(:, EOP_DYYY) * U_tild(:, EOP_DXX) &
                               + psi_tmp(:, EOP_DYY) * U_tild(:, EOP_DXXY) ) &
                               - psi_tmp(:, EOP_DY) &
                               *( psi_tmp(:, EOP_DXXY) * U_tild(:, EOP_DYY) &
                               +psi_tmp(:, EOP_DXY) * U_tild(:, EOP_DXYY) &
                               +psi_tmp(:, EOP_DXYY) * U_tild(:, EOP_DXY) &
                               +psi_tmp(:, EOP_DYY) * U_tild(:, EOP_DXXY) )
         fn_eval_elm = fn_eval_elm + 2 * dt*dt * fn_eval_elm_tmp
       end if
       ! <laplace U, U>
       if (linear .eq. 0) then
          !U_tmp = eph179_pre
          U_tmp =eph179
       else
          U_tmp = eph079
       end if
       call laplace_gs_r (U_tmp, fn_eval_tmp, icylinder, 0)
       call laplace_gs_z (U_tmp, fn_eval_tmp2, icylinder, 0)
       call laplace_gs_r (U_2bar, fn_eval_tmp3, icylinder, 0)
       call laplace_gs_z (U_2bar, fn_eval_tmp4, icylinder, 0)

       fn_eval_elm_tmp (:) = fn_eval_tmp*U_2bar(:, EOP_DY) &
                             -fn_eval_tmp2*U_2bar(:, EOP_DX) &
                             +fn_eval_tmp3*U_tmp(:, EOP_DY) &
                             -fn_eval_tmp4*U_tmp(:, EOP_DX)
       fn_eval_elm = fn_eval_elm - 0.5* en179 (:, EOP_1) * fn_eval_elm_tmp
       
       ! -<psi, laplace psi>
       !if (isplitstep .eq. 1 ) then
        !  psi_tmp = eps179_pre
       !else
          psi_tmp = psi_2bar
       !end if 

      if (linear .eq. 0) then
          !psi_tmp2 = eps179_pre
          psi_tmp2 = eps179
       else
          psi_tmp2 = eps079
       end if
       call laplace_gs_r (psi_tmp2, fn_eval_tmp, icylinder, 0)
       call laplace_gs_z (psi_tmp2, fn_eval_tmp2, icylinder, 0)
       call laplace_gs_r (psi_tmp, fn_eval_tmp3, icylinder, 0)
       call laplace_gs_z (psi_tmp, fn_eval_tmp4, icylinder,0)
       fn_eval_elm_tmp (:) = fn_eval_tmp*psi_tmp(:, EOP_DY) &
                             -fn_eval_tmp2*psi_tmp(:, EOP_DX) &
                             +fn_eval_tmp3*psi_tmp2(:, EOP_DY) &
                             -fn_eval_tmp4*psi_tmp2(:, EOP_DX)
       fn_eval_elm = fn_eval_elm + 0.5* fn_eval_elm_tmp

       ! time derivative
       call laplace_gs (U_delta, fn_eval_tmp, icylinder, 0)
       fn_eval_elm_tmp (:) = en179(:, EOP_1)*fn_eval_tmp/dt
       fn_eval_elm = fn_eval_elm + fn_eval_elm_tmp

#ifdef USECOMPLEX
       elm_res_U (itri) = int1 (fn_eval_elm) *CONJG (int1 (fn_eval_elm)) * ((d%a+d%b)*d%c)**(2-iadapt_useH1)
       !elm_res_U (itri) = real(int1 (fn_eval_elm)) **2* ((d%a+d%b)*d%c)**2
#else
       elm_res_U (itri) = int1 (fn_eval_elm) **2 * ((d%a+d%b)*d%c)**(2-iadapt_useH1)
#endif
       if (isplitstep .eq. 0 .or. .true.) then
          fn_eval_elm = 0
          ! eta psi
          if(jadv .eq. 0) then
            call laplace_gs (psi_bar, fn_eval_tmp, icylinder, icylinder)
          else
            call biharmic_gs_op(psi_bar, fn_eval_tmp, icylinder, icylinder)
          end if
          fn_eval_elm = fn_eval_elm - eetar79(:,EOP_1)*fn_eval_tmp
          ![psi, U]
          if(jadv .eq. 0) then
            fn_eval_elm_tmp (:) = psi_2bar(:, EOP_DX)*eph179_pre(:, EOP_DY) &
              -psi_2bar(:, EOP_DY)*eph179_pre(:, EOP_DX) &
              +psi_tmp2(:, EOP_DX)*U_2bar(:, EOP_DY) &
              -psi_tmp2(:, EOP_DY)* U_2bar(:, EOP_DX)
              fn_eval_elm = fn_eval_elm - fn_eval_elm_tmp 
          else
         !  <psi_R, u_R> 
           fn_eval_elm_tmp (:) = psi_2bar(:, EOP_DXX) * eph179(:, EOP_DXY) &
                                 - psi_2bar(:, EOP_DXY) * eph179(:, EOP_DXX)
           fn_eval_elm_tmp =   2*fn_eval_elm_tmp
           fn_eval_elm = fn_eval_elm - fn_eval_elm_tmp

           !  <psi_Z, u_Z>
           fn_eval_elm_tmp (:) = psi_2bar(:, EOP_DXY) * eph179(:, EOP_DYY) &
                                 - psi_2bar(:, EOP_DYY) * eph179(:, EOP_DXY)
           fn_eval_elm_tmp =   2*fn_eval_elm_tmp
           fn_eval_elm = fn_eval_elm - fn_eval_elm_tmp

           ! <nalba psi, u>
           call laplace_gs_r (psi_2bar, fn_eval_tmp, icylinder, 0)
           call laplace_gs_z (psi_2bar, fn_eval_tmp2, icylinder, 0)
           fn_eval_elm_tmp (:) = fn_eval_tmp*eph179(:, EOP_DY)-fn_eval_tmp2*eph179(:, EOP_DX)  
           fn_eval_elm = fn_eval_elm - fn_eval_elm_tmp
           ! <psi, nabla u>
           call laplace_gs_r (eph179, fn_eval_tmp, icylinder, 0)
           call laplace_gs_z (eph179, fn_eval_tmp2, icylinder, 0)
           fn_eval_elm_tmp (:) = fn_eval_tmp2*psi_2bar(:, EOP_DX)-fn_eval_tmp*psi_2bar(:, EOP_DY) 
           fn_eval_elm = fn_eval_elm - fn_eval_elm_tmp
          end if
          ! time derivative
          if(jadv .eq. 0) then
            fn_eval_elm = fn_eval_elm +psi_delta(:, EOP_1)/dt
          else
            call laplace_gs (psi_delta, fn_eval_tmp, icylinder, icylinder)
            fn_eval_elm = fn_eval_elm + fn_eval_tmp/dt
          end if
#ifdef USECOMPLEX
          elm_res_psi (itri) = int1 (fn_eval_elm) *CONJG(int1 (fn_eval_elm)) * ((d%a+d%b)*d%c)**(2-iadapt_useH1)
          !elm_res_psi (itri) = real(int1 (fn_eval_elm))**2 * ((d%a+d%b)*d%c)**2
#else
          elm_res_psi (itri) = int1 (fn_eval_elm) **2 * ((d%a+d%b)*d%c)**(2-iadapt_useH1)
#endif
       end if
    end do
    buff = max_current
    call mpi_allreduce (buff, max_current, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier )
    buff = min_current
    call mpi_allreduce (buff, min_current, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier )
    if(myrank .eq. 0) print *, "max current", max_current,min_current
    buff2Norm =  solutionH2Norm
    call mpi_allreduce (buff2Norm, solutionH2Norm, 2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ier )

  end subroutine elem_residule

  subroutine rotate_field (val, c, s)
    implicit none
    real :: c, s
    vectype, dimension(MAX_PTS, EOP_NUM) :: val, valCopy
    valCopy = val

    val(:, EOP_DX) = valCopy(:, EOP_DX)*c - valCopy(:, EOP_DY)*s
    val(:, EOP_DY) = valCopy(:, EOP_DX)*s + valCopy(:, EOP_DY)*c

    val(:, EOP_DXX) = valCopy(:, EOP_DXX)*c**2 - 2.*valCopy(:, EOP_DXY)*s*c+valCopy(:, EOP_DYY)*s**2
    val(:, EOP_DXY) = valCopy(:, EOP_DXX)*c*s + valCopy(:, EOP_DXY)*(c**2-s**2) - valCopy(:, EOP_DYY)*s*c
    val(:, EOP_DYY) = valCopy(:, EOP_DXX)*s**2 + 2.*valCopy(:, EOP_DXY)*s*c+valCopy(:, EOP_DYY)*c**2

    val(:, EOP_DXXX) = valCopy(:, EOP_DXXX)*c**3 - 3.*valCopy(:, EOP_DXXY)*s*c*c+3.*valCopy(:, EOP_DXYY)*s*s*c - valCopy(:, EOP_DYYY)*s**3
    val(:, EOP_DXXY) = valCopy(:, EOP_DXXX)*c*c*s + valCopy(:, EOP_DXXY)*(c**3-2*s*s*c) + valCopy(:, EOP_DXYY)*(s**3-2*s*c*c) + valCopy(:, EOP_DYYY)*s*s*c
    val(:, EOP_DXYY) = valCopy(:, EOP_DXXX)*c*s*s + valCopy(:, EOP_DXXY)*(-s**3+2*s*c*c) + valCopy(:, EOP_DXYY)*(c**3-2*s*s*c) - valCopy(:, EOP_DYYY)*s*c*c
    val(:, EOP_DYYY) = valCopy(:, EOP_DXXX)*s**3 + 3.*valCopy(:, EOP_DXXY)*s*s*c+3.*valCopy(:, EOP_DXYY)*s*c*c + valCopy(:, EOP_DYYY)*c**3

    val(:, EOP_DXXXX) = valCopy(:, EOP_DXXXX)*c**4 - 4.*valCopy(:, EOP_DXXXY)*s*c*c*c+6.*valCopy(:, EOP_DXXYY)*s*s*c*c - 4.*valCopy(:, EOP_DXYYY)*c*s**3+valCopy(:, EOP_DYYYY)*s**4
    val(:, EOP_DXXXY) = valCopy(:, EOP_DXXXX)*c**3*s + valCopy(:, EOP_DXXXY)*(c**4-3.*s*s*c*c) + 3.*valCopy(:, EOP_DXXYY)*(-c**3*s+c*s**3) + valCopy(:, EOP_DXYYY)*(3.*c*c*s*s-s**4) -valCopy(:, EOP_DYYYY)*s**3*c
    val(:, EOP_DXXYY) = valCopy(:, EOP_DXXXX)*c*c*s*s + valCopy(:, EOP_DXXXY)*(2.*s*c**3-2.*s**3*c) + valCopy(:, EOP_DXXYY)*(c**4+s**4-4.*s*s*c*c) + valCopy(:, EOP_DXYYY)*(2.*c*s*s*s-2*s*c**3) +valCopy(:, EOP_DYYYY)*s*s*c*c
    val(:, EOP_DXYYY) = valCopy(:, EOP_DXXXX)*s**3*c + valCopy(:, EOP_DXXXY)*(-s**4+3.*s*s*c*c) + 3.*valCopy(:, EOP_DXXYY)*(-s**3*c+s*c**3) + valCopy(:, EOP_DXYYY)*(-3.*c*c*s*s+c**4) -valCopy(:, EOP_DYYYY)*c**3*s
    val(:, EOP_DYYYY) = valCopy(:, EOP_DXXXX)*s**4 + 4.*valCopy(:, EOP_DXXXY)*c*s*s*s+6.*valCopy(:, EOP_DXXYY)*s*s*c*c + 4.*valCopy(:, EOP_DXYYY)*s*c**3+valCopy(:, EOP_DYYYY)*c**4
  end subroutine rotate_field 
  
  subroutine laplace_gs (val, ans, cylinder, gs)
    implicit none
    integer :: cylinder, gs
    vectype, dimension(MAX_PTS, EOP_NUM) :: val
    vectype, dimension(MAX_PTS) :: ans
    ans(:)=val(:, EOP_DXX)+val(:, EOP_DYY)
    if(cylinder .eq. 1) then
      if(gs .eq. 0) then
        ans(:)=ans(:)+1./er_79*val(:,EOP_DX)
      else
        ans(:)=ans(:)-(1./er_79)*val(:,EOP_DX)
      end if
    end if
  end subroutine laplace_gs 
  
  subroutine laplace_gs_r (val, ans, cylinder, gs)
    implicit none
    integer :: cylinder, gs
    vectype, dimension(MAX_PTS, EOP_NUM) :: val
    vectype, dimension(MAX_PTS) :: ans
    ans(:)=val(:, EOP_DXXX)+val(:, EOP_DXYY)
    if(cylinder .eq. 1) then
      if(gs .eq. 0) then
        ans(:)=ans(:)+1./er_79*val(:,EOP_DXX)-(1./er_79)**2*val(:,EOP_DX)
      else
        ans(:)=ans(:)-1./er_79*val(:,EOP_DXX)+(1./er_79)**2*val(:,EOP_DX)
      end if
    end if
  end subroutine laplace_gs_r 

  subroutine laplace_gs_rr (val, ans, cylinder, gs)
    implicit none
    integer :: cylinder, gs
    vectype, dimension(MAX_PTS, EOP_NUM) :: val
    vectype, dimension(MAX_PTS) :: ans
    ans(:)=val(:, EOP_DXXXX)+val(:, EOP_DXXYY)
    if(cylinder .eq. 1) then
      if(gs .eq. 0) then
        ans(:)=ans(:)+1./er_79*val(:,EOP_DXXX)-(2./er_79)**2*val(:,EOP_DXX)+ (2.0/er_79)**3*val(:,EOP_DX)
      else
        ans(:)=ans(:)-1./er_79*val(:,EOP_DXXX)+(2./er_79)**2*val(:,EOP_DXX)-(2.0/er_79)**3*val(:,EOP_DX)
      end if
    end if
  end subroutine laplace_gs_rr

  subroutine laplace_gs_z (val, ans, cylinder, gs)
    implicit none
    integer :: cylinder, gs
    vectype, dimension(MAX_PTS, EOP_NUM) :: val
    vectype, dimension(MAX_PTS) :: ans
    ans(:)=val(:, EOP_DXXY)+val(:, EOP_DYYY)
    if(cylinder .eq. 1) then
      if(gs .eq. 0) then
        ans(:)=ans(:)+1./er_79*val(:,EOP_DXY)
      else
        ans(:)=ans(:)-1./er_79*val(:,EOP_DXY)
      end if
    end if
  end subroutine laplace_gs_z

 subroutine laplace_gs_zz (val, ans, cylinder, gs)
    implicit none
    integer :: cylinder, gs
    vectype, dimension(MAX_PTS, EOP_NUM) :: val
    vectype, dimension(MAX_PTS) :: ans
    ans(:)=val(:, EOP_DXXYY)+val(:, EOP_DYYYY)
    if(cylinder .eq. 1) then
      if(gs .eq. 0) then
        ans(:)=ans(:)+1./er_79*val(:,EOP_DXYY)
      else
        ans(:)=ans(:)-1./er_79*val(:,EOP_DXYY)
      end if
    end if
  end subroutine laplace_gs_zz 

 subroutine laplace_gs_rz (val, ans, cylinder, gs)
    implicit none
    integer :: cylinder, gs
    vectype, dimension(MAX_PTS, EOP_NUM) :: val
    vectype, dimension(MAX_PTS) :: ans
    ans(:)=val(:, EOP_DXXXY)+val(:, EOP_DXYYY)
    if(cylinder .eq. 1) then
      if(gs .eq. 0) then
        ans(:)=ans(:)+1./er_79*val(:,EOP_DXXY)-(1./er_79)**2*val(:,EOP_DXY)
      else
        ans(:)=ans(:)-1./er_79*val(:,EOP_DXXY)+(1./er_79)**2*val(:,EOP_DXY)
      end if
    end if
  end subroutine laplace_gs_rz

  subroutine biharmic_gs_op (val, ans, cylinder, gs)
    implicit none
    integer :: cylinder, gs
    vectype, dimension(MAX_PTS, EOP_NUM) :: val
    vectype, dimension(MAX_PTS) :: ans
    ans(:)=val(:, EOP_DXXXX)+val(:, EOP_DYYYY)+2*val(:, EOP_DXXYY)
    if(cylinder .eq. 1) then
      if(gs .eq. 0) then
        ans(:)=ans(:)+2./er_79*val(:,EOP_DXXX)-val(:, EOP_DXX)+1./er_79*val(:, EOP_DX)+ 1./er_79 * (val(:, EOP_DXYY)+val(:, EOP_DXYY))
     else
        ans(:)=ans(:)+val(:, EOP_DXX)-1./er_79*val(:, EOP_DX)+ 1./er_79 * (val(:, EOP_DXYY)-val(:, EOP_DXYY))
     end if
    end if
  end subroutine biharmic_gs_op
end module error_estimate 
