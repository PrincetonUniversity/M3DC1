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
  vectype, dimension(MAX_PTS, EOP_NUM) :: eph179, eph179_pre, eps179, eps179_pre, en179, evis79 
  vectype,  dimension(MAX_PTS, EOP_NUM) :: U_bar, U_2bar, U_tild, psi_bar, psi_2bar, psi_tmp 
  vectype,  dimension(MAX_PTS, EOP_NUM) :: U_delta, psi_delta
  real, private :: max_current, min_current
  integer, parameter :: LPUN = 1 ! dt * mu * partial (laplace U)/ parital n
  integer, parameter :: LPUT = 2 ! dt * ro * laplace U * partial U / partial t
  integer, parameter :: LPPSPST = 3 ! dt * laplace psi * partial psi  / partial t
  integer, parameter :: LPU = 4 ! dt * mu  laplace U; R2 term
  integer, parameter :: TSQ = 5  ! dt*dt term in R1
  integer, parameter :: R1 = 6
  integer, parameter :: R3 = 7
  integer, parameter :: TEST = 8 ! psi_dx + psi_dy
  integer, parameter :: NUMTERM = 9
  integer, parameter :: npoint_int =3
  integer, parameter :: npoint_int_elm = 12
  vectype, allocatable :: edge_jump (:,:,:)
  vectype, dimension (npoint_int) :: fn_eval
  vectype, dimension (npoint_int_elm) :: fn_eval_elm
  vectype, dimension (npoint_int_elm) :: fn_eval_elm_tmp
  contains
  subroutine jump_discontinuity (edge_error)
    use basic
    use arrays
    use scorec_mesh_mod
    use vector_mod
    use m3dc1_nint
    integer :: itri, numelms,numedgs, iedge, ii, node_next, num_get 
    vectype, dimension(:,:), intent(out) :: edge_error
    integer, dimension(3) :: edges, nodes, nodes_edge, edge_dir, idimgeo, idimgeo_t, is_bdy
    real, dimension(3) :: edge_len
    real, dimension(2,3) :: normal, normal_t
    real, dimension(3) :: coords
    integer, allocatable :: edge_tag(:)
    type(element_data) :: d

    call m3dc1_mesh_getnument(1, numedgs)

    allocate(edge_jump(npoint_int, numedgs, NUMTERM))
    edge_jump = 0
    !print*, "shape edge_jump", shape(edge_jump)
    idimgeo=2;
    numelms = local_elements()
    do itri=1,numelms 
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

          call rotate_vec(normal(:,iedge), d%sn, d%co )
          !print *, 'ele dir',d%co, d%sn
          !print *, 'rotate', normal(:,iedge)
          ! caculate jump discontinuity

          !lp_U_n
          fn_eval(1:npoint_int)=evis79(1:npoint_int,EOP_1)*((U_bar(1:npoint_int,EOP_DXXX)+ U_bar(1:npoint_int,EOP_DXYY))*normal(1,iedge)+ (U_bar(1:npoint_int,EOP_DXXY) + U_bar(1:npoint_int,EOP_DYYY))*normal(2,iedge))

          if(edge_dir(iedge) .eq. 0) call reverse_fn(fn_eval(1:npoint_int),npoint_int)
          edge_jump(1:npoint_int,edges(iedge)+1, LPUN)= edge_jump(1:npoint_int,edges(iedge)+1, LPUN)+ dt*fn_eval(1:npoint_int)
          !print *, edges(iedge), 'dir',edge_dir(iedge)
          !print *, fn_eval(1:5)

          !rho lp_U * U_t
          fn_eval(1:npoint_int)=en179(1:npoint_int,EOP_1)*(U_2bar(1:npoint_int,EOP_DXX)+U_2bar(1:npoint_int,EOP_DYY))*(eph179_pre(1:npoint_int,EOP_DX)*normal(2,iedge)- eph179_pre(1:npoint_int,EOP_DY)*normal(1,iedge))
          fn_eval(1:npoint_int)=fn_eval(1:npoint_int)+en179(1:npoint_int,EOP_1)*(eph179_pre(1:npoint_int,EOP_DXX)+eph179_pre(1:npoint_int,EOP_DYY))*(U_2bar(1:npoint_int,EOP_DX)*normal(2,iedge)- U_2bar(1:npoint_int,EOP_DY)*normal(1,iedge))

          if(edge_dir(iedge) .eq. 0) call reverse_fn(fn_eval(1:npoint_int),npoint_int)
          edge_jump(1:npoint_int,edges(iedge)+1, LPUT)=edge_jump(1:npoint_int,edges(iedge)+1,LPUT)+0.5*dt*fn_eval(1:npoint_int)

          !lp_psi * Psi_t
          if (isplitstep .eq. 1) then
             psi_tmp = eps179_pre
          else
             psi_tmp = psi_2bar 
          end if
          fn_eval(1:npoint_int) = (psi_tmp(1:npoint_int,EOP_DXX)+ psi_tmp(1:npoint_int,EOP_DYY))*(eps179_pre(1:npoint_int,EOP_DX)*normal(2,iedge)-eps179_pre(1:npoint_int,EOP_DY)*normal(1,iedge))
          fn_eval(1:npoint_int) = fn_eval(1:npoint_int) + (eps179_pre(1:npoint_int,EOP_DXX)+ eps179_pre(1:npoint_int,EOP_DYY))*(psi_tmp(1:npoint_int,EOP_DX)*normal(2,iedge)-psi_tmp(1:npoint_int,EOP_DY)*normal(1,iedge))
          if(edge_dir(iedge) .eq. 0) call reverse_fn(fn_eval(1:npoint_int),npoint_int)
          edge_jump(1:npoint_int,edges(iedge)+1, LPPSPST)=edge_jump(1:npoint_int,edges(iedge)+1,LPPSPST)+0.5*dt*fn_eval(1:npoint_int)

          !lp_U
          fn_eval(1:npoint_int) = evis79(1:npoint_int, EOP_1)*(U_bar(1:npoint_int,EOP_DXX)+U_bar(1:npoint_int,EOP_DYY))
          if (edge_dir(iedge) .eq. 0) then
             call reverse_fn(fn_eval(1:npoint_int),npoint_int)
             fn_eval=-fn_eval;
          end if
          edge_jump(1:npoint_int,edges(iedge)+1, LPU)=edge_jump(1:npoint_int,edges(iedge)+1, LPU)+dt*fn_eval(1:npoint_int)

          !test jump field should be zero
          ! test field dR+dZ
          fn_eval(1:npoint_int) = eps179(1:npoint_int,EOP_DX)*d%co-eps179(1:npoint_int,EOP_DY)*d%sn 
          fn_eval(1:npoint_int) = fn_eval(1:npoint_int) + eps179(1:npoint_int,EOP_DX)*d%sn+eps179(1:npoint_int,EOP_DY)*d%co
          if (edge_dir(iedge) .eq. 0) then
             call reverse_fn(fn_eval(1:npoint_int),npoint_int)
             fn_eval=-fn_eval;
          end if
          edge_jump(1:npoint_int,edges(iedge)+1, TEST)=edge_jump(1:npoint_int,edges(iedge)+1, TEST)+dt*fn_eval(1:npoint_int)
          !print *, 'test', edges(iedge)+1, fn_eval(1:npoint_int)
          ! include delta t terms for isplitstep =1
          if (isplitstep .eq. 1) then
             fn_eval(1:npoint_int) = (eps179_pre(1:npoint_int,EOP_DX)*normal(2,iedge) &
                      - eps179_pre(1:npoint_int,EOP_DY)*normal(1,iedge)) &
                      *( (U_tild(1:npoint_int,EOP_DXXX)+ U_tild(1:npoint_int,EOP_DXYY)) * eps179_pre(1:npoint_int,EOP_DY) &
                         -(U_tild(1:npoint_int,EOP_DXXY)+U_tild(1:npoint_int,EOP_DYYY))* eps179_pre(1:npoint_int,EOP_DX) &
                         + (eps179_pre(1:npoint_int,EOP_DXXX)+eps179_pre(1:npoint_int,EOP_DXYY)) * U_tild(1:npoint_int,EOP_DY) &
                         - (eps179_pre(1:npoint_int,EOP_DXXY)+ eps179_pre(1:npoint_int,EOP_DYYY)) * U_tild(1:npoint_int,EOP_DX) & 
                         + 2* U_tild(1:npoint_int,EOP_DXX)*eps179_pre(1:npoint_int,EOP_DXY) &
                         - 2 * U_tild(1:npoint_int,EOP_DXY)*eps179_pre(1:npoint_int,EOP_DYY) &
                         + 2* U_tild(1:npoint_int,EOP_DXY)*eps179_pre(1:npoint_int,EOP_DYY) &
                         - 2 * U_tild(1:npoint_int,EOP_DYY)*eps179_pre(1:npoint_int,EOP_DXY) )

             if (edge_dir(iedge) .eq. 0) then
                call reverse_fn(fn_eval(1:npoint_int),npoint_int)
             end if
             edge_jump(1:npoint_int,edges(iedge)+1, TSQ) = edge_jump(1:npoint_int,edges(iedge)+1, TSQ) + dt*dt*fn_eval(1:npoint_int) 
          end if
          edge_jump(1:npoint_int,edges(iedge)+1,R1) =  edge_jump(1:npoint_int,edges(iedge)+1, LPUN) &
                                                       + edge_jump(1:npoint_int,edges(iedge)+1, TSQ) &
                                                       - edge_jump(1:npoint_int,edges(iedge)+1, LPUT) &
                                                       + edge_jump(1:npoint_int,edges(iedge)+1, LPPSPST)
          ! R3
          if(isplitstep .eq. 1) then
            fn_eval(1:npoint_int) = (U_tild(1:npoint_int,EOP_DXX) + U_tild(1:npoint_int,EOP_DYY)) * (eps179_pre(1:npoint_int,EOP_DX)*normal(2,iedge)-eps179_pre(1:npoint_int,EOP_DY)*normal(1,iedge)) + (eps179_pre(1:npoint_int,EOP_DXX) + eps179_pre(1:npoint_int,EOP_DYY)) * (U_tild(1:npoint_int,EOP_DX)*normal(2,iedge)-U_tild(1:npoint_int,EOP_DY)*normal(1,iedge)) 
            if (edge_dir(iedge) .eq. 0) then
               call reverse_fn(fn_eval(1:npoint_int),npoint_int)
            end if
            edge_jump(1:npoint_int,edges(iedge)+1,R3) = edge_jump(1:npoint_int,edges(iedge)+1,R3) + dt*dt*fn_eval
          end if
       end do
    end do
    do ii=1,NUMTERM
       !if(ii .eq. 5) print*,myrank, edge_jump(1:npoint_int,:,ii)
       call sum_edge_data(edge_jump(1:npoint_int,:,ii),npoint_int)
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
#else
             edge_error(edges(iedge)+1,ii)=int1(edge_jump(1:npoint_int,edges(iedge)+1,ii)**2)
#endif
          end do
          edge_error(edges(iedge)+1,LPUN)=edge_error(edges(iedge)+1,LPUN)* edge_len(1)**3 
          edge_error(edges(iedge)+1,LPUT)=edge_error(edges(iedge)+1,LPUT)* edge_len(1)**3 
          edge_error(edges(iedge)+1,LPPSPST)=edge_error(edges(iedge)+1,LPPSPST)* edge_len(1)**3 
          edge_error(edges(iedge)+1,LPU)=edge_error(edges(iedge)+1,LPU)* edge_len(1)
          edge_error(edges(iedge)+1,TSQ)=edge_error(edges(iedge)+1,TSQ)* edge_len(1)**3 
          edge_error(edges(iedge)+1,R1)=edge_error(edges(iedge)+1,R1)* edge_len(1)**3 
          edge_error(edges(iedge)+1,R3)=edge_error(edges(iedge)+1,R3)* edge_len(1) 
          edge_error(edges(iedge)+1,NUMTERM)=edge_error(edges(iedge)+1,R1)+edge_error(edges(iedge)+1,LPU)+edge_error(edges(iedge)+1,LPPSPST)+edge_error(edges(iedge)+1,R3)
       end do
    end do
    deallocate (edge_tag)
    deallocate(edge_jump)
  end subroutine jump_discontinuity

  subroutine reverse_fn( f, n)
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
       r_79 = x_79
    else
       r_79 = 1.
    endif
    ri_79 = 1./r_79

    if(ijacobian.eq.1) weight_79 = weight_79 * r_79

    call precalculate_terms_error(xi_79,eta_79,npoints)
    call define_basis_error(itri)

    call eval_ops_error(itri, u_field(1), eph179)
    call eval_ops_error(itri, psi_field(1), eps179)
    call eval_ops_error(itri, den_field(1), en179)

    call eval_ops_error(itri, u_field_pre, eph179_pre)
    call eval_ops_error(itri, psi_field_pre, eps179_pre)
    U_bar = thimp * eph179 +(1. -thimp) * eph179_pre
    U_tild = thimp*thimp*eph179 + thimp*(1-thimp)*eph179_pre
    U_2bar = 2*thimp* eph179 +(1. -2*thimp) * eph179_pre
    psi_bar = thimp * eph179 +(1. -thimp) * eph179_pre
    psi_2bar = 2*thimp* eps179 +(1. -2*thimp) * eps179_pre
    U_delta = eph179 - eph179_pre
    psi_delta = eps179 - eps179_pre
    !if(eqsubtract.eq.1) then
      !n179 = 1 
    !end if
    call eval_ops_error(itri, visc_field, evis79)
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

  subroutine define_basis_error(itri)
    use basic
    implicit none

    integer, intent(in) :: itri
    real, dimension(dofs_per_element,coeffs_per_element) :: cl

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
    integer :: itri, numelms, ier
    vectype, dimension(:), intent(out) :: elm_res_U, elm_res_psi
    type(element_data) :: d
    real :: buff
    real, dimension(npoint_int_elm) :: buff2
    numelms = local_elements()
    max_current = 0
    min_current = 0
    do itri=1,numelms
       call get_element_data(itri, d)
       call define_element_quadrature(itri, npoint_int_elm, 0)
       call define_fields_error(itri)
       buff2=eps179(1:npoint_int_elm, EOP_DXX)+eps179(1:npoint_int_elm, EOP_DYY)
       buff = maxval(buff2)
       if(buff>max_current) max_current = buff 
       buff = minval(buff2)
       if(buff<min_current) min_current = buff
       fn_eval_elm = 0
       fn_eval_elm_tmp = 0
       fn_eval_elm_tmp (1:npoint_int_elm) = U_bar(1:npoint_int_elm, EOP_DXXXX) + U_bar(1:npoint_int_elm, EOP_DYYYY) + 2*U_bar(1:npoint_int_elm, EOP_DXXYY)
       fn_eval_elm_tmp (1:npoint_int_elm) = dt*evis79(1:npoint_int_elm, EOP_1)*fn_eval_elm_tmp (1:npoint_int_elm)
       fn_eval_elm = fn_eval_elm + fn_eval_elm_tmp

       if (isplitstep .eq. 1) then
         !<psi, <psi, laplace U > >
         fn_eval_elm_tmp (1:npoint_int_elm) = eps179_pre(1:npoint_int_elm, EOP_DX) &
                                            *( eps179_pre(1:npoint_int_elm, EOP_DXY) * (U_tild(1:npoint_int_elm, EOP_DXXY)+U_tild(1:npoint_int_elm, EOP_DYYY)) &
                                              +eps179_pre(1:npoint_int_elm, EOP_DX) * (U_tild(1:npoint_int_elm, EOP_DXXYY)+U_tild(1:npoint_int_elm, EOP_DYYYY)) &
                                              -eps179_pre(1:npoint_int_elm, EOP_DYY) * (U_tild(1:npoint_int_elm, EOP_DXXX)+U_tild(1:npoint_int_elm, EOP_DXYY)) & 
                                              +eps179_pre(1:npoint_int_elm, EOP_DY) * (U_tild(1:npoint_int_elm, EOP_DXXXY)+U_tild(1:npoint_int_elm, EOP_DXYYY)) ) &
                                           -eps179_pre(1:npoint_int_elm, EOP_DY) & 
                                            *( eps179_pre(1:npoint_int_elm, EOP_DXX) * (U_tild(1:npoint_int_elm, EOP_DXXY)+U_tild(1:npoint_int_elm, EOP_DYYY)) & 
                                              +eps179(1:npoint_int_elm, EOP_DX) * (U_tild(1:npoint_int_elm, EOP_DXXXY)+U_tild(1:npoint_int_elm, EOP_DXYYY)) &
                                              -eps179_pre(1:npoint_int_elm, EOP_DXY) * (U_tild(1:npoint_int_elm, EOP_DXXX)+U_tild(1:npoint_int_elm, EOP_DXYY)) & 
                                              +eps179_pre(1:npoint_int_elm, EOP_DY) * (U_tild(1:npoint_int_elm, EOP_DXXXX)+U_tild(1:npoint_int_elm, EOP_DXXYY)) )
         fn_eval_elm = fn_eval_elm + dt*dt * fn_eval_elm_tmp
         ! <psi, <laplace psi, U> >
         fn_eval_elm_tmp (1:npoint_int_elm) = -eps179_pre(1:npoint_int_elm, EOP_DX) &
                                            *( ph179(1:npoint_int_elm, EOP_DXY) * (eps179_pre(1:npoint_int_elm, EOP_DXXY)+eps179_pre(1:npoint_int_elm, EOP_DYYY)) &
                                              +U_tild(1:npoint_int_elm, EOP_DX) * (eps179_pre(1:npoint_int_elm, EOP_DXXYY)+eps179_pre(1:npoint_int_elm, EOP_DYYYY)) &
                                              -U_tild(1:npoint_int_elm, EOP_DYY) * (eps179_pre(1:npoint_int_elm, EOP_DXXX)+eps179_pre(1:npoint_int_elm, EOP_DXYY)) &
                                              +U_tild(1:npoint_int_elm, EOP_DY) * (eps179_pre(1:npoint_int_elm, EOP_DXXXY)+eps179_pre(1:npoint_int_elm, EOP_DXYYY)) ) &
                                           +eps179_pre(1:npoint_int_elm, EOP_DY) &
                                            *( U_tild(1:npoint_int_elm, EOP_DXX) * (eps179_pre(1:npoint_int_elm, EOP_DXXY) + eps179_pre(1:npoint_int_elm, EOP_DYYY)) &
                                              +U_tild(1:npoint_int_elm, EOP_DX) * (eps179_pre(1:npoint_int_elm, EOP_DXXXY)+eps179_pre(1:npoint_int_elm, EOP_DXYYY)) &
                                              -U_tild(1:npoint_int_elm, EOP_DXY) * (eps179_pre(1:npoint_int_elm, EOP_DXXX)+eps179_pre(1:npoint_int_elm, EOP_DXYY)) &
                                              +U_tild(1:npoint_int_elm, EOP_DY) * (eps179_pre(1:npoint_int_elm, EOP_DXXXX)+eps179_pre(1:npoint_int_elm, EOP_DXXYY)) )
         fn_eval_elm = fn_eval_elm + dt*dt * fn_eval_elm_tmp
         ! <psi, <psi_R, u_R> >
         fn_eval_elm_tmp (1:npoint_int_elm) = eps179_pre(1:npoint_int_elm, EOP_DX) &
                                            *( eps179_pre(1:npoint_int_elm, EOP_DXXY) * U_tild(1:npoint_int_elm, EOP_DXY) &
                                              +eps179_pre(1:npoint_int_elm, EOP_DXX) * U_tild(1:npoint_int_elm, EOP_DXYY) &
                                              +eps179_pre(1:npoint_int_elm, EOP_DXYY) * U_tild(1:npoint_int_elm, EOP_DXX) &
                                              +eps179_pre(1:npoint_int_elm, EOP_DXY) * U_tild(1:npoint_int_elm, EOP_DXXY) ) &
                                            -eps179_pre(1:npoint_int_elm, EOP_DY) &
                                            *( eps179_pre(1:npoint_int_elm, EOP_DXXX) * U_tild(1:npoint_int_elm, EOP_DXY) &
                                              +eps179_pre(1:npoint_int_elm, EOP_DXX) * U_tild(1:npoint_int_elm, EOP_DXXY) &
                                              +eps179_pre(1:npoint_int_elm, EOP_DXXY) * U_tild(1:npoint_int_elm, EOP_DXX) &
                                              +eps179_pre(1:npoint_int_elm, EOP_DXY) * U_tild(1:npoint_int_elm, EOP_DXXX) )
         fn_eval_elm = fn_eval_elm + 2 * dt*dt * fn_eval_elm_tmp

         ! <psi, <psi_Z, u_Z> >
         fn_eval_elm_tmp (1:npoint_int_elm) = eps179_pre(1:npoint_int_elm, EOP_DX) &
                                            *( eps179_pre(1:npoint_int_elm, EOP_DXYY) * U_tild(1:npoint_int_elm, EOP_DXY) &
                                              + eps179_pre(1:npoint_int_elm, EOP_DXY) * U_tild(1:npoint_int_elm, EOP_DXYY) &
                                              +eps179_pre(1:npoint_int_elm, EOP_DYYY) * U_tild(1:npoint_int_elm, EOP_DXX) &
                                              +eps179_pre(1:npoint_int_elm, EOP_DYY) * U_tild(1:npoint_int_elm, EOP_DXXY) ) &
                                            -eps179_pre(1:npoint_int_elm, EOP_DY) &
                                            *( eps179_pre(1:npoint_int_elm, EOP_DXXY) * U_tild(1:npoint_int_elm, EOP_DYY) &
                                              +eps179_pre(1:npoint_int_elm, EOP_DXY) * U_tild(1:npoint_int_elm, EOP_DXYY) &
                                              +eps179_pre(1:npoint_int_elm, EOP_DXYY) * U_tild(1:npoint_int_elm, EOP_DXY) &
                                              +eps179_pre(1:npoint_int_elm, EOP_DYY) * U_tild(1:npoint_int_elm, EOP_DXXY) )
         fn_eval_elm = fn_eval_elm + 2 * dt*dt * fn_eval_elm_tmp
       end if
       ! <laplace U, U>
       fn_eval_elm_tmp (1:npoint_int_elm) = (eph179_pre(1:npoint_int_elm, EOP_DXXX)+eph179_pre(1:npoint_int_elm, EOP_DXYY))*U_2bar(1:npoint_int_elm, EOP_DY) &
                                            -(eph179_pre(1:npoint_int_elm, EOP_DXXY)+ eph179_pre(1:npoint_int_elm, EOP_DYYY))*U_2bar(1:npoint_int_elm, EOP_DX) &
                                            +(U_2bar(1:npoint_int_elm, EOP_DXXX)+U_2bar(1:npoint_int_elm, EOP_DXYY))*eph179_pre(1:npoint_int_elm, EOP_DY) &
                                            -(U_2bar(1:npoint_int_elm, EOP_DXXY)+U_2bar(1:npoint_int_elm, EOP_DYYY))*eph179_pre(1:npoint_int_elm, EOP_DX)
       fn_eval_elm = fn_eval_elm - 0.5*dt* en179 (1:npoint_int_elm, EOP_1) * fn_eval_elm_tmp
       
       ! -<psi, laplace psi>
       if (isplitstep .eq. 1) then
          psi_tmp = eps179_pre
       else
          psi_tmp = psi_2bar
       end if 

       fn_eval_elm_tmp (1:npoint_int_elm) = (eps179_pre(1:npoint_int_elm, EOP_DXXX)+eps179_pre(1:npoint_int_elm, EOP_DXYY))*psi_tmp(1:npoint_int_elm, EOP_DY) &
                                            -(eps179_pre(1:npoint_int_elm, EOP_DXXY)+eps179_pre(1:npoint_int_elm, EOP_DYYY))*psi_tmp(1:npoint_int_elm, EOP_DX) &
                                            +(psi_tmp(1:npoint_int_elm, EOP_DXXX)+psi_tmp(1:npoint_int_elm, EOP_DXYY))*eps179_pre(1:npoint_int_elm, EOP_DY) &
                                            -(psi_tmp(1:npoint_int_elm, EOP_DXXY)+psi_tmp(1:npoint_int_elm, EOP_DYYY))*eps179_pre(1:npoint_int_elm, EOP_DX)
       fn_eval_elm = fn_eval_elm + 0.5*dt* fn_eval_elm_tmp

       ! time derivative
       fn_eval_elm_tmp (1:npoint_int_elm) = en179(1:npoint_int_elm, EOP_1)*U_delta(1:npoint_int_elm, EOP_DXX)+U_delta(1:npoint_int_elm, EOP_DYY)
       fn_eval_elm = fn_eval_elm - fn_eval_elm_tmp

       elm_res_U (itri) = int1 (fn_eval_elm) **2 * ((d%a+d%b)*d%c)**2
      
       if (isplitstep .eq. 0) then
          fn_eval_elm = 0
          ! eta psi
          fn_eval_elm_tmp (1:npoint_int_elm) = psi_bar(1:npoint_int_elm, EOP_DXX)+psi_bar(1:npoint_int_elm, EOP_DYY)
          fn_eval_elm = fn_eval_elm + eta*dt*fn_eval_elm_tmp
          ![psi, U]
          fn_eval_elm_tmp (1:npoint_int_elm) = psi_2bar(1:npoint_int_elm, EOP_DX)*eph179_pre(1:npoint_int_elm, EOP_DY) &
              -psi_2bar(1:npoint_int_elm, EOP_DY)*eph179_pre(1:npoint_int_elm, EOP_DX) &
              +eps179_pre(1:npoint_int_elm, EOP_DX)*U_2bar(1:npoint_int_elm, EOP_DY) &
              -eps179_pre(1:npoint_int_elm, EOP_DY)* U_2bar(1:npoint_int_elm, EOP_DX) 
          fn_eval_elm = fn_eval_elm - 0.5*dt*fn_eval_elm_tmp
          ! time derivative
          fn_eval_elm = fn_eval_elm -psi_delta(1:npoint_int_elm, EOP_1)
          elm_res_psi (itri) = int1 (fn_eval_elm) **2 * ((d%a+d%b)*d%c)**2
       end if
    end do
    buff = max_current
    call mpi_allreduce (buff, max_current, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier )
    buff = min_current
    call mpi_allreduce (buff, min_current, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier )
    if(myrank .eq. 0) print *, "max current", max_current,min_current
  end subroutine elem_residule
end module error_estimate 
