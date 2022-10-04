      module c1_12_element
        implicit none
        integer, parameter :: DOF_1   = 1
        integer, parameter :: DOF_DR  = 2
        integer, parameter :: DOF_DZ  = 3
        integer, parameter :: DOF_NORAML  = 1

#ifdef USE3D
#endif

        integer, parameter :: maxpol = 3
        integer, parameter :: pol_dofs_per_node = 3
        integer, parameter :: pol_dofs_per_edge = 1 
        integer, parameter :: pol_nodes_per_element = 3
#ifdef USE3D
#else
        integer, parameter :: tor_dofs_per_node = 1
        integer, parameter :: tor_nodes_per_element = 1
        integer, parameter :: maxtor = 1
        integer, parameter :: edges_per_element = 3
        integer, parameter :: coeffs_per_dphi = 1
        integer, parameter :: dofs_per_dphi = 1
#endif
        integer, parameter :: dofs_per_node = tor_dofs_per_node*pol_dofs_per_node
        integer, parameter :: nodes_per_element = &
          pol_nodes_per_element*tor_nodes_per_element
        integer, parameter :: dofs_per_tri = 12
        integer, parameter :: nodes_per_edge = 2
        integer, parameter :: dofs_per_element = nodes_per_element*(dofs_per_node + pol_dofs_per_edge)
        ! for compile
        integer, parameter :: coeffs_per_element=20,coeffs_per_tri=20
        integer :: ni(coeffs_per_tri),mi(coeffs_per_tri)


        real, dimension(3,3) :: coord_glb
        integer :: itri_tag
        real, dimension(3,3) :: dx_dedge, dy_dedge
        real :: jacobi 
        ! norm direction of each edge
        ! consistent for adjencent elements
        ! edge order: 12 23 31
        real, dimension(3,2) :: normal_edges
        ! length of edges
        real, dimension(3):: len_edges
        ! index of shape functions
        integer, dimension(12), private :: idx_node, idx_loc
        ! derivative value of N_i in the normal direction of edge mid
        real, dimension(9,3), private :: normal_mid_edge
        real, dimension(3,3), private :: coord_mid_edge
        ! the sig sign to define edge normal
        ! if outnormal == global normal -1
        ! else 1
        real, dimension(3)  :: sig_edges

        ! for first derivative from (d/db12)-->(d/dxy)
        ! db1/dx db2/dx; db1/dy db2/dy
        real, dimension(2,2) :: db12_dxy_1
        ! for second derivative from (d^2/db12)-->(d^2/dxy)
        real, dimension(3,3) :: db12_dxy_2
        ! c matrix to convert from local to global c1 (b123->xy)
        real, dimension(2,2,3) :: loc2glb_basis
        data idx_node /1,1,1,2,2,2,3,3,3,1,2,3 /
        data idx_loc /1,2,3,1,2,3,1,2,3,4,4,4/
        data coord_mid_edge /0.5,0.5,0.0,0.0,0.5,0.5,0.5,0.0,0.5/
        integer, dimension(3) :: is_bounary_edges
        real, allocatable :: normal_all_edges(:,:)

        contains
        subroutine set_element_coord ( itri_p, coord, is_bounary_edges_p)
          real, dimension(3,3), intent(in) :: coord ! three nodes, xyz
          integer, dimension(3), intent(in) :: is_bounary_edges_p
          integer , intent(in) :: itri_p

          coord_glb=coord
          itri_tag=itri_p
          is_bounary_edges=is_bounary_edges_p

          dx_dedge(1,2) = coord(1,1)-coord(2,1)
          dx_dedge(1,3) = coord(1,1)-coord(3,1)
          dx_dedge(2,1) = -dx_dedge(1,2)
          dx_dedge(2,3) = coord(2,1)-coord(3,1)
          dx_dedge(3,1) = -dx_dedge(1,3)
          dx_dedge(3,2) = -dx_dedge(2,3)

          dy_dedge(1,2) = coord(1,2)-coord(2,2)
          dy_dedge(1,3) = coord(1,2)-coord(3,2)
          dy_dedge(2,1) = -dy_dedge(1,2)
          dy_dedge(2,3) = coord(2,2)- coord(3,2)
          dy_dedge(3,1) = -dy_dedge(1,3)
          dy_dedge(3,2) = -dy_dedge(2,3)
          jacobi=dx_dedge(1,3)*dy_dedge(2,3)+dx_dedge(3,2)*dy_dedge(1,3)
          if( jacobi .le. 0) then
            write(*,*) 'warning! clock-wise node order'
            !call abort
          end if

          len_edges(1)=sqrt(dx_dedge(1,2)**2+dy_dedge(1,2)**2)
          len_edges(2)=sqrt(dx_dedge(2,3)**2+dy_dedge(2,3)**2)
          len_edges(3)=sqrt(dx_dedge(1,3)**2+dy_dedge(1,3)**2)

          call calcu_db12_dxy ()
          call calcu_loc2glb_basis ()
          call calcu_normal_edges ()
          call calcu_normal_mid_edge ()
          !write (*,*), 'elem', itri_p, ' sig_edges', sig_edges
          !write (*,*), 'normal', normal_edges(1,:),normal_edges(2,:),normal_edges(3,:)
        end subroutine set_element_coord

        subroutine calcu_db12_dxy ()
          db12_dxy_1(1,1)=dy_dedge(2,3)
          db12_dxy_1(1,2)=dy_dedge(3,1)
          db12_dxy_1(2,1)=dx_dedge(3,2)
          db12_dxy_1(2,2)=dx_dedge(1,3)
          db12_dxy_1=db12_dxy_1/jacobi

          db12_dxy_2(1,1)=db12_dxy_1(1,1)*db12_dxy_1(1,1)
          db12_dxy_2(1,2)=2*db12_dxy_1(1,1)*db12_dxy_1(1,2)
          db12_dxy_2(1,3)=db12_dxy_1(1,2)*db12_dxy_1(1,2)
 
          db12_dxy_2(2,1)=db12_dxy_1(1,1)*db12_dxy_1(2,1)
          db12_dxy_2(2,2)=db12_dxy_1(1,1)*db12_dxy_1(2,2)+db12_dxy_1(1,2)*db12_dxy_1(2,1)
          db12_dxy_2(2,3)=db12_dxy_1(1,2)*db12_dxy_1(2,2)

          db12_dxy_2(3,1)=db12_dxy_1(2,1)*db12_dxy_1(2,1)
          db12_dxy_2(3,2)=2*db12_dxy_1(2,1)*db12_dxy_1(2,2)
          db12_dxy_2(3,3)=db12_dxy_1(2,2)*db12_dxy_1(2,2)
        end subroutine calcu_db12_dxy
 
        subroutine calcu_normal_edges ()
          integer :: iedge, inode1, inode2, inode3
          real :: edge32_dot_normal_edge
          do iedge=1,3
            inode1=iedge
            inode2=inode1+1
            if (inode2 .gt. 3) inode2=1
            inode3=inode1-1
            if (inode3 .lt. 1) inode3=3
            if (dy_dedge(inode1,inode2) .gt. 0) then
              normal_edges(iedge,1)=dy_dedge(inode1,inode2)/len_edges(iedge)
              normal_edges(iedge,2)=-dx_dedge(inode1,inode2)/len_edges(iedge)
            else
              if (dy_dedge(inode1,inode2) .lt. 0) then
                normal_edges(iedge,1)=-dy_dedge(inode1,inode2)/len_edges(iedge)
                normal_edges(iedge,2)=dx_dedge(inode1,inode2)/len_edges(iedge)
              else
                normal_edges(iedge,1)=0
                normal_edges(iedge,2)=1
              end if
            end if

            ! define sig, compare with out normal
            edge32_dot_normal_edge=dx_dedge(inode2,inode3)*normal_edges(iedge,1)+dy_dedge(inode2,inode3)*normal_edges(iedge,2)
            ! if on the bounary, use outer normal
            if (is_bounary_edges(iedge) .eq. 1) then
              if (edge32_dot_normal_edge .lt. 0) then
                normal_edges(iedge,:)=-1*normal_edges(iedge,:)
              end if
              sig_edges(iedge) = -1
            else
              if (edge32_dot_normal_edge .gt. 0) then 
                sig_edges(iedge) = -1 
              else
                sig_edges(iedge) = 1
              end if
           end if
          end do
        end subroutine calcu_normal_edges

        subroutine calcu_normal_mid_edge ()
          integer :: counter, iedge
          real, dimension(3) ::res
          real, dimension(2) :: gradxy
          do counter = 1,9
            do iedge = 1,3
              call eval_n_i ( counter, coord_mid_edge(:,iedge), 1, res)
              ! first get gradient with respect to xy
              gradxy(1)=db12_dxy_1(1,1)*res(1)+db12_dxy_1(1,2)*res(2)
              gradxy(2)=db12_dxy_1(2,1)*res(1)+db12_dxy_1(2,2)*res(2)
              ! get deriv in edge normal n dot gradxy
              normal_mid_edge(counter,iedge)=normal_edges(iedge,1)*gradxy(1)+normal_edges(iedge,2)*gradxy(2)
              if (sig_edges(iedge) .lt. 0) normal_mid_edge(counter,iedge)=-1*normal_mid_edge(counter,iedge)
            end do
            !write (*,*) "shape fn ", counter, "normal_mid_edge",normal_mid_edge(counter,:)
          end do  
        end subroutine calcu_normal_mid_edge

        subroutine calcu_loc2glb_basis ()
          integer :: inode, edge1, edge2
          do inode=1,3
            edge1=inode+1
            if (edge1 .gt. 3) then
              edge1=1
            end if
            edge2=inode-1
            if (edge2 .lt. 1) then
              edge2=3
            end if
            loc2glb_basis(1,1,inode)=dx_dedge(edge1,inode)
            loc2glb_basis(2,1,inode)=dy_dedge(edge1,inode)
            loc2glb_basis(1,2,inode)=dx_dedge(edge2,inode)
            loc2glb_basis(2,2,inode)=dy_dedge(edge2,inode)
          end do
        end subroutine calcu_loc2glb_basis

        ! xyz_loc: b1, b2
        subroutine glb_xyz_2_loc ( xyz_glb, xyz_loc)
          real, dimension(3), intent(in) :: xyz_glb
          real, dimension(3), intent(out) :: xyz_loc
          xyz_loc(1)=dy_dedge(2,3)*(xyz_glb(1)-coord_glb(3,1))+dx_dedge(3,2)*(xyz_glb(2)-coord_glb(3,2))
          xyz_loc(2)=dy_dedge(3,1)*(xyz_glb(1)-coord_glb(3,1))+dx_dedge(1,3)*(xyz_glb(2)-coord_glb(3,2))
          xyz_loc(1)=xyz_loc(1)/jacobi
          xyz_loc(2)=xyz_loc(2)/jacobi
          !xyz_loc(3)=1.0-xyz_loc(1)-xyz_loc(2)
        end subroutine glb_xyz_2_loc

        subroutine loc_xyz_2_glb ( xyz_loc, xyz_glb)
          real, dimension(3), intent(in) :: xyz_loc
          real, dimension(3), intent(out) :: xyz_glb
          real :: b3
          b3=1.0-xyz_loc(1)-xyz_loc(2)
          xyz_glb(1)=xyz_loc(1)*coord_glb(1,1)+xyz_loc(2)*coord_glb(2,1)+b3*coord_glb(3,1)
          xyz_glb(2)=xyz_loc(1)*coord_glb(1,2)+xyz_loc(2)*coord_glb(2,2)+b3*coord_glb(3,2)
        end subroutine loc_xyz_2_glb

        ! RM Koch 1996
        ! order: 0, value; 1, first deriv, b1, b2
        ! 2, second deriv, b1b1, b1b2, b2b2
        ! 3,  higher order for surface term, to do
        subroutine eval_n_i ( ith, xyz_loc, order, res)
          real, dimension(3), intent(in) :: xyz_loc
          integer, intent(in) :: ith, order
          real, dimension(3), intent(out) :: res
          real, dimension(3) :: res_tmp
          real, dimension(3) :: xyz_loc_tmp
          integer :: inode, ith_loc, r_th, s_th
          real :: r, s

          xyz_loc_tmp=xyz_loc
          xyz_loc_tmp(3)=1.-xyz_loc_tmp(1)-xyz_loc_tmp(2)
          inode=idx_node(ith)
          ith_loc=idx_loc(ith)
          r_th=inode
          s_th=inode+1
          if(s_th .gt. 3) s_th=1

          r=xyz_loc_tmp(r_th)
          s=xyz_loc_tmp(s_th)

          select case (ith_loc)
            case (1)
              call eval_n_1(r,s, order, res)
            case (2)
              call eval_n_2(r,s, order, res)
            case (3)
              call eval_n_3(r,s, order, res)
            case (4)
              call eval_n_10(r,s, order, res)
            case default
              write(*,*) "invalid input ith for function eval_n_i"
              call abort
          end select
             
          select case (order)
            case (0)
              return
            case (1)
              select case (inode)
                case (1)
                  return
                case (2)
                  res_tmp=res
                  res(1)=-res_tmp(2)
                  res(2)=res_tmp(1)-res_tmp(2)
                case (3)
                  res_tmp=res
                  res(1)=-res_tmp(1)+res_tmp(2)
                  res(2)=-res_tmp(1)
                case default
                  write(*,*) "invalid input for function eval_n_i"
                  call abort
              end select  
            case (2)
              select case (inode)
                case (1)
                  return
                case (2)
                  res_tmp=res
                  res(1)=res_tmp(3)
                  res(2)=-res_tmp(2)+res_tmp(3)
                  res(3)=res_tmp(1)-2*res_tmp(2)+res_tmp(3)
                case (3)
                  res_tmp=res
                  res(1)=res_tmp(1)-2*res_tmp(2)+res_tmp(3)
                  res(2)=res_tmp(1)-res_tmp(2)
                  res(3)=res_tmp(1)
                case default
                  write(*,*) "invalid input for function eval_n_i"
                  call abort
              end select
            case default
              write(*,*) "invalid input for function eval_n_i"
              call abort
          end select
        end subroutine eval_n_i

        ! RH Koch 1996
        ! order: 0, value; 1, first deriv, tangent, normal
        ! default: normal= x
        ! if v=x, v=[1,0]; u=y, u=[0,1]
        ! 2, second deriv, nn, nv, vv
        ! 3,  higher order for surface term, to do
        ! r,s,t !=1., or singular at vertex for ith= 10,11,12
        subroutine eval_basis_function( ith, xyz_loc, order, res, isnode, normal)
          real, dimension(3), intent(in) :: xyz_loc
          integer, intent(in) :: ith, order, isnode
          real, dimension(3), intent(out) :: res
          real, dimension(2), optional :: normal
          real, dimension(3) :: res_tmp, res_tmp2, res_tmp3
          ! values of n_10,n_11,n_12 at the r,s
          ! (res, ith)
          real, dimension(3,3) :: singular_n_values
          real, dimension(2,2) :: b2xy_1st
          real :: sig_rs
          integer :: inode, ith_loc,edge1,edge2, counter
          real :: dx_dedge1, dy_dedge1, dx_dedge2, dy_dedge2

          inode=idx_node(ith)
          ith_loc=idx_loc(ith)

          edge1=inode+1
          if (edge1 .gt. 3) then
            edge1=1
          end if
          edge2=inode-1
          if (edge2 .lt. 1) then
            edge2=3
          end if
          dx_dedge1=dx_dedge(edge1,inode)
          dy_dedge1=dy_dedge(edge1,inode)
          dx_dedge2=dx_dedge(edge2,inode)
          dy_dedge2=dy_dedge(edge2,inode)
           
          ! applt new rot if on the bounbdary
          if (present(normal)) then
            dx_dedge1=normal(1)*dx_dedge(edge1,inode)+normal(2)*dy_dedge(edge1,inode)
            dx_dedge2=normal(1)*dx_dedge(edge2,inode)+normal(2)*dy_dedge(edge2,inode)
            dy_dedge1=-normal(2)*dx_dedge(edge1,inode)+normal(1)*dy_dedge(edge1,inode)
            dy_dedge2=-normal(2)*dx_dedge(edge2,inode)+normal(1)*dy_dedge(edge2,inode)
          end if

          ! get 4*jacobi/len_edges(rs)*N_i, i=10,11,12
          if ( ith .lt. 10 .and. isnode .ne. 1) then
           do counter=1,3
              call eval_n_i (9+counter, xyz_loc, order, singular_n_values(:,counter))
              singular_n_values(:,counter)=4.*jacobi/len_edges(counter)*singular_n_values(:,counter)
           end do
          end if

          select case (ith_loc)
            case (1)
              call eval_n_i (ith, xyz_loc, order, res_tmp)
              if ( ith .lt. 10 .and. isnode .ne. 1) then
                do counter=1,3
                  res_tmp=res_tmp-normal_mid_edge(ith,counter)*singular_n_values(:,counter)
                end do
              end if
            case (2)
              call eval_n_i (ith, xyz_loc, order, res_tmp)
              call eval_n_i (ith+1, xyz_loc, order, res_tmp2)
              if ( ith .lt. 10 .and. isnode .ne. 1) then
                do counter=1,3
                  res_tmp=res_tmp-normal_mid_edge(ith,counter)*singular_n_values(:,counter)
                  res_tmp2=res_tmp2-normal_mid_edge(ith+1,counter)*singular_n_values(:,counter)
                end do
              end if
              res_tmp=dx_dedge1*res_tmp+dx_dedge2*res_tmp2
            case (3)
              call eval_n_i (ith-1, xyz_loc, order, res_tmp)
              call eval_n_i (ith, xyz_loc, order, res_tmp2)
              if ( ith .lt. 10 .and. isnode .ne. 1) then
                do counter=1,3
                  res_tmp=res_tmp-normal_mid_edge(ith-1,counter)*singular_n_values(:,counter)
                  res_tmp2=res_tmp2-normal_mid_edge(ith,counter)*singular_n_values(:,counter)
                end do
              end if
              res_tmp=dy_dedge1*res_tmp+dy_dedge2*res_tmp2
            case (4) 
              call eval_n_i (ith, xyz_loc, order, res_tmp)
              res_tmp=res_tmp*4.*jacobi/len_edges(inode)*sig_edges(inode)
          end select
          
          ! transform to dirivatives with respect to xy
          select case (order)
            case (0)
              res(1)=res_tmp(1)
            case (1)
              res(1)=db12_dxy_1(1,1)* res_tmp(1)+db12_dxy_1(1,2)* res_tmp(2)
              res(2)=db12_dxy_1(2,1)* res_tmp(1)+db12_dxy_1(2,2)* res_tmp(2)
            case (2)
              res(1)=db12_dxy_2(1,1)* res_tmp(1)+db12_dxy_2(1,2)* res_tmp(2)+ db12_dxy_2(1,3)* res_tmp(3)
              res(2)=db12_dxy_2(2,1)* res_tmp(1)+db12_dxy_2(2,2)* res_tmp(2)+ db12_dxy_2(2,3)* res_tmp(3)
              res(3)=db12_dxy_2(3,1)* res_tmp(1)+db12_dxy_2(3,2)* res_tmp(2)+ db12_dxy_2(3,3)* res_tmp(3)
          end select
        end subroutine eval_basis_function

        subroutine eval_n_1 ( b1, b2, order, res)
          real, intent(in) :: b1,b2
          integer, intent(in) :: order
          real, dimension(3), intent(out) :: res
          real :: b3

          b3=1.-b1-b2
          select case (order)
            case (0)
              res(1)=b1+b1*b1*b2+b1*b1*b3-b1*b2*b2-b1*b3*b3
            case (1)
              res(1)=-6*b1*b1 - 4*b1*b2 + 6*b1 - 2*b2*b2 + 2*b2
              res(2)=(-2)*b1*(b1 + 2*b2 - 1)
            case (2)
              res(1)=6 - 4*b2 - 12*b1
              res(2)=2 - 4*b2 - 4*b1
              res(3)=(-4)*b1
            case default
              write(*,*) "invalid input order for function eval_n_1"
              call abort
          end select
        end subroutine eval_n_1
        
        subroutine eval_n_2 ( b1, b2, order, res)
          real, intent(in) :: b1,b2
          integer, intent(in) :: order
          real, dimension(3), intent(out) :: res
          real :: b3

          b3=1.-b1-b2
          select case (order)
            case (0)
              res(1)=b1*b1*b2+1.0/2.0*b1*b2*b3
            case (1)
              res(1)=(b2*(2*b1 - b2 + 1))/2.
              res(2)=(b1*(b1 - 2*b2 + 1))/2.
            case (2)
              res(1)=b2
              res(2)=b1 - b2 + 1./2.
              res(3)=-b1
            case default
              write(*,*) "invalid input order for function eval_n_2"
              call abort
            end select
        end subroutine eval_n_2
 
       subroutine eval_n_3 ( b1, b2, order, res)
          real, intent(in) :: b1,b2
          integer, intent(in) :: order
          real, dimension(3), intent(out) :: res
          real :: b3

          b3=1.-b1-b2
          select case (order)
            case (0)
              res(1)=b1*b1*b3+1.0/2.0*b1*b2*b3
            case (1)
              res(1)=-3*b1*b1 - 3*b1*b2 + 2*b1 - b2*b2/2. + b2/2.
              res(2)=-(b1*(3*b1 + 2*b2 - 1))/2.
            case (2)
              res(1)=2 - 3*b2 - 6*b1
              res(2)=1./2. - b2 - 3*b1
              res(3)=-b1
            case default
              write(*,*) "invalid input order for function eval_n_3"
              call abort
            end select
        end subroutine eval_n_3

       subroutine eval_n_10 ( b1, b2, order, res)
          real, intent(in) :: b1,b2
          integer, intent(in) :: order
          real, dimension(3), intent(out) :: res
          real :: b3

          b3=1.-b1-b2
          select case (order)
            case (0)
              res(1)=b1*b1*b2*b2*b3*(1+b3)/(b1+b3)/(b2+b3)
            case (1)
              res(1)=(b1*b2**2*(3*b1**3 + 4*b1**2*b2 - 10*b1**2 + b1*b2**2 & 
                     - 9*b1*b2 + 11*b1 - 2*b2**2 + 6*b2 - 4))/((b1 - 1)**2*(b2 - 1))
              res(2)=(b1**2*b2*(b1**2*b2 - 2*b1**2 + 4*b1*b2**2 - 9*b1*b2 + 6*b1 &
 		     + 3*b2**3 - 10*b2**2 + 11*b2 - 4))/((b1 - 1)*(b2 - 1)**2)
            case (2)
              res(1)=(2*b2**3)/(b1 - 1)**3 + (4*b2**3 - 22*b2**2 &
 		    + (6*b2**2*(3*b2 - 3))/(b2 - 1))/(b2 - 1) + (6*b1*b2**2)/(b2 - 1)
              res(2)=(b1*b2*(3*b1**3*b2 - 6*b1**3 + 8*b1**2*b2**2 - 22*b1**2*b2 &
                     + 20*b1**2 + 3*b1*b2**3 - 22*b1*b2**2 + 38*b1*b2 - 22*b1 &
                     - 6*b2**3 + 20*b2**2 - 22*b2 + 8))/((b1 - 1)**2*(b2 - 1)**2)
              res(3)=(2*b1**3)/(b2 - 1)**3 + (4*b1**3 - 22*b1**2 + &
	             (6*b1**2*(3*b1 - 3))/(b1 - 1))/(b1 - 1) + (6*b1**2*b2)/(b1 - 1)
            case default
              write(*,*) "invalid input order for function eval_n_3"
              call abort
            end select
        end subroutine eval_n_10

        subroutine xy_glb_mid_edge (iedge, xy_mid)
          integer, intent(in):: iedge
          real, dimension(3), intent(out) ::xy_mid
          integer :: node1, node2
          node1=iedge
          node2=iedge+1
          if (node2 .gt. 3) node2=1
          xy_mid(:)=0.5*(coord_glb(node1,:)+coord_glb(node2,:))
        end subroutine xy_glb_mid_edge 
   
        subroutine dof_gradxy2midedge_normal ( iedge, gradxy, midedge_normal, xy_mid)
          integer, intent(in):: iedge
          real, intent(out) :: midedge_normal
          real, dimension(2),intent(in) :: gradxy
          real, dimension(2), intent(out) ::xy_mid
          midedge_normal=normal_edges(iedge,1)*gradxy(1)+normal_edges(iedge,2)*gradxy(2)
        end subroutine dof_gradxy2midedge_normal

        ! evalutate field at local coord (r,s)
        ! dofs: local
        ! order: global (x,y)  0 value; 1 partial x; 2 partial y; 11, xx; 12, xy, 22, yy 
        subroutine eval_field( r, s, dofs, values, inode)
          real, intent(in) :: r,s
          real, dimension(12), intent(in) :: dofs
          real, dimension(6), intent(out) :: values
          integer, optional, intent(in) :: inode
          real :: values_tmp,values_tmp2,values_tmp3
          integer :: counter, num_basis, isnode
          real, dimension(3) :: xyz_loc, res
      
          isnode=0
          xyz_loc(1)=r
          xyz_loc(2)=s
          xyz_loc(3)=1-r-s
          num_basis=12
          if( present(inode) ) then
            num_basis=9
            isnode=1
          end if
          values=0.
          do counter=1,num_basis
            call eval_basis_function( counter, xyz_loc, 0, res, isnode)
            values(1)=values(1)+dofs(counter)*res(1)
            call eval_basis_function( counter, xyz_loc, 1, res, isnode)
            values(2:3)=values(2:3)+dofs(counter)*res(1:2)
            call eval_basis_function( counter, xyz_loc, 2, res, isnode)
            values(4:6)=values(4:6)+dofs(counter)*res(1:3)
          end do 
        end subroutine eval_field
      end module c1_12_element

      subroutine getedgenormal (iedge, normal)
        use c1_12_element
        integer, intent(in) :: iedge
        real, dimension(2), intent(out):: normal
        normal(:)=normal_all_edges(:,iedge)
      end subroutine getedgenormal
