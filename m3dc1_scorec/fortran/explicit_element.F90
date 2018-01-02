      module explicit_element
        implicit none
        real, dimension(3,3) :: coord_glb
        real :: lambda_1, lambda_2, lambda_3
        real, dimension(3,3) :: lambda 
        real, dimension(3,3) :: dx_dedge, dy_dedge, dedge_dx, dedge_dy ! b1b2b3,three node
        real :: jacobi ! b1 b2--> x1 x2, directed area *2
        integer, parameter :: num_terms=3 ! max coffs for each shape fn
        integer :: itri_tag
        ! used for pre caculation
        ! transpose of glb2loc  matrix
        ! glb dofs: interior: xy
        ! boundary: normal, tangent (really tricky here)
        real, allocatable :: glb2loc1st_t_elem(:,:,:,:)
        real, allocatable :: glb2loc2nd_t_elem(:,:,:,:)
        real, allocatable :: glb2loc3rd_t_elem(:,:,:,:)
        integer, allocatable :: boundary_info(:,:)
        contains
        subroutine set_element_coord ( itri_p, coord)
          real, dimension(3,3), intent(in) :: coord ! three nodes, xyz
          integer , intent(in) :: itri_p
          real, dimension(2) :: e1,e2,e3
          coord_glb=coord
          itri_tag=itri_p
          e1(1)=coord(3,1) - coord(2,1)
          e1(2)=coord(3,2) - coord(2,2)
          e2(1)=coord(1,1) - coord(3,1)
          e2(2)=coord(1,2) - coord(3,2)
          e3(1)=coord(2,1) - coord(1,1)
          e3(2)=coord(2,2) - coord(1,2)
          !-dot(e3,e1)/norm(e1)^2
          lambda_1 = -(e3(1)*e1(1)+e3(2)*e1(2))/(e1(1)*e1(1)+e1(2)*e1(2)) 
          !-dot(e1,e2)/norm(e2)^2
          lambda_2 = -(e1(1)*e2(1)+e1(2)*e2(2))/(e2(1)*e2(1)+e2(2)*e2(2))
          !-dot(e2,e3)/norm(e3)^2
          lambda_3 = -(e2(1)*e3(1)+e2(2)*e3(2))/(e3(1)*e3(1)+e3(2)*e3(2))

          lambda(1,2) = 1 - lambda_3
          lambda(2,1) = lambda_3
          lambda(1,3) = lambda_2
          lambda(3,1) = 1 - lambda_2
          lambda(2,3) = 1 - lambda_1
          lambda(3,2) = lambda_1

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

          dedge_dx(1,2) = dy_dedge(2,3)/jacobi
          dedge_dx(1,3) = dy_dedge(2,3)/jacobi
          dedge_dx(2,1) = dy_dedge(3,1)/jacobi 
          dedge_dx(2,3) = dy_dedge(3,1)/jacobi
          dedge_dx(3,1) = dy_dedge(1,2)/jacobi
          dedge_dx(3,2) = dy_dedge(1,2)/jacobi

          dedge_dy(1,2) = dx_dedge(3,2)/jacobi 
          dedge_dy(1,3) = dx_dedge(3,2)/jacobi 
          dedge_dy(2,1) = dx_dedge(1,3)/jacobi
          dedge_dy(2,3) = dx_dedge(1,3)/jacobi
          dedge_dy(3,1) = dx_dedge(2,1)/jacobi
          dedge_dy(3,2) = dx_dedge(2,1)/jacobi

        end subroutine set_element_coord

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

        subroutine glb_dof_2_loc ( dofs_glb, dofs_loc )
          real, dimension(18), intent(in):: dofs_glb
          real, dimension(18), intent(out):: dofs_loc
          integer :: inode,start_idx, edge1, edge2
          real :: dx_dedge1, dx_dedge2,dy_dedge1, dy_dedge2
          do inode=1,3
            start_idx=(inode-1)*6+1
            edge1=inode+1
            if (edge1 .gt. 3) then 
              edge1=edge1-3
            end if
            edge2=inode-1
            if (edge2 .lt. 1) then
              edge2=edge2+3
            end if 
            dx_dedge1=dx_dedge(edge1,inode)
            dy_dedge1=dy_dedge(edge1,inode)
            dx_dedge2=dx_dedge(edge2,inode)
            dy_dedge2=dy_dedge(edge2,inode) 
            ! field value
            dofs_loc(start_idx)=dofs_glb(start_idx)
            ! first derivatives
            dofs_loc(start_idx+1)=dx_dedge1*dofs_glb(start_idx+1)+dy_dedge1*dofs_glb(start_idx+2)
            dofs_loc(start_idx+2)=dx_dedge2*dofs_glb(start_idx+1)+dy_dedge2*dofs_glb(start_idx+2)
            ! second derivatives
            dofs_loc(start_idx+3)=dx_dedge1*dx_dedge1*dofs_glb(start_idx+3) &
                                  +2*dx_dedge1*dy_dedge1*dofs_glb(start_idx+4) &
                                  +dy_dedge1*dy_dedge1*dofs_glb(start_idx+5)
            dofs_loc(start_idx+4)=dx_dedge1*dx_dedge2*dofs_glb(start_idx+3) &
                                  +(dx_dedge1*dy_dedge2+dx_dedge2*dy_dedge1)*dofs_glb(start_idx+4) &
                                  +dy_dedge1*dy_dedge2*dofs_glb(start_idx+5)
            dofs_loc(start_idx+5)=dx_dedge2*dx_dedge2*dofs_glb(start_idx+3) &
                                  +2*dx_dedge2*dy_dedge2*dofs_glb(start_idx+4) &
                                  +dy_dedge2*dy_dedge2*dofs_glb(start_idx+5)
          end do
        end subroutine glb_dof_2_loc

        subroutine loc_dof_2_glb ( dofs_loc, dofs_glb )
          real, dimension(18), intent(in):: dofs_loc
          real, dimension(18), intent(out):: dofs_glb
          integer :: inode,start_idx, edge1, edge2
          real :: dedge1_dx, dedge2_dx,dedge1_dy, dedge2_dy
          do inode=1,3
            start_idx=(inode-1)*6+1
            edge1=inode+1
            if (edge1 .gt. 3) then
              edge1=edge1-3
            end if
            edge2=inode-1
            if (edge2 .lt. 1) then
              edge2=edge2+3
            end if

            dedge1_dx=dedge_dx(edge1,inode)
            dedge1_dy=dedge_dy(edge1,inode)
            dedge2_dx=dedge_dx(edge2,inode)
            dedge2_dy=dedge_dy(edge2,inode)
            ! field value
            dofs_glb(start_idx)=dofs_loc(start_idx)
            ! first derivatives
            dofs_glb(start_idx+1)=dedge1_dx*dofs_loc(start_idx+1)+dedge2_dx*dofs_loc(start_idx+2)
            dofs_glb(start_idx+2)=dedge1_dy*dofs_loc(start_idx+1)+dedge2_dy*dofs_loc(start_idx+2)
            ! second derivatives
            dofs_glb(start_idx+3)=dedge1_dx*dedge1_dx*dofs_loc(start_idx+3) &
                                  +2*dedge1_dx*dedge2_dx*dofs_loc(start_idx+4) &
                                  +dedge2_dx*dedge2_dx*dofs_loc(start_idx+5)
            dofs_glb(start_idx+4)=dedge1_dx*dedge1_dy*dofs_loc(start_idx+3) &
                                  +(dedge1_dx*dedge2_dy+dedge1_dy*dedge2_dx)*dofs_loc(start_idx+4) &
                                  +dedge2_dx*dedge2_dy*dofs_loc(start_idx+5)
            dofs_glb(start_idx+5)=dedge1_dy*dedge1_dy*dofs_loc(start_idx+3) &
                                  +2*dedge1_dy*dedge2_dy*dofs_loc(start_idx+4) &
                                  +dedge2_dy*dedge2_dy*dofs_loc(start_idx+5)
          end do
        end subroutine loc_dof_2_glb
       
       
        ! position funtion at vertex 1
        ! order: 0 value; 1 partial r; 2 partial s; 11, r^2; 12, rs, 22, ss
        ! output terms(1)+lambda21*(terms(2))+lambda31*(terms(3))
        subroutine eval_p1(r,s,order, terms ) 
          real, intent(in) :: r,s
          integer, intent(in) :: order
          real, dimension(:), intent(out) :: terms
          real :: r_sq, r_cb !r^2, r^3
          real :: s_sq, s_cb !s^2, s^3,
          real :: b3, b3_sq, b3_cb
          b3=1-r-s
          b3_sq=b3*b3
          b3_cb=b3_sq*b3
          r_sq=r*r
          r_cb=r_sq*r
          s_sq=s*s
          s_cb=s_sq*s

          select case (order)
            case (0)
              terms(1)=r+r_sq*s+r_sq*b3-r*s_sq-r*b3_sq+3*r_cb*s_sq &
                       +4*r_cb*s*b3+3*r_cb*b3_sq-3*r_sq*s_cb &
                       -2*r*s_cb*b3-3*r_sq*b3_cb-2*r*s*b3_cb &
                       -4*r*s_sq*b3_sq+17*r_sq*s*b3_sq+17*r_sq*s_sq*b3
              terms(2)=-30*r_sq*s_sq*b3
              terms(3)=-30*r_sq*s*b3_sq
            case (1)
              terms(1)=1-r_sq+2*r*s-s_sq+4*r*b3-b3_sq-4*r_cb*s-6*r_cb*b3 &
                       -8*r_sq*s_sq-4*r*s_cb+18*r_sq*b3_sq+40*r*s*b3_sq &
                       +42*r*s_sq*b3-22*r_sq*s*b3-2*s_cb*b3-6*r*b3_cb &
                       -2*s*b3_cb-4*s_sq*b3_sq
              terms(2)=30*r_sq*s_sq-60*r*s_sq*b3
              terms(3)=60*r_sq*s*b3-60*r*s*b3_sq
            case (2)
              terms(1)=-2*r*s+2*r*b3+2*r_cb*s-2*r_cb*b3-26*r_sq*s_sq &
                       +2*r*s_sq*b3+2*r*s_cb+26*r_sq*b3_sq-2*r*b3_cb &
                       -2*r*s*b3_sq
              terms(2)=-60*r_sq*s*b3+30*r_sq*s_sq
              terms(3)=-30*r_sq*b3_sq+60*r_sq*s*b3
            case (11)
              terms(1)=6+10*r_sq*s-54*r_sq*b3+46*s*b3_sq+54*r*b3_sq &
                       -58*r*s_sq-4*s-12*r+6*r_cb-124*r*s*b3-2*s_cb &
                       -6*b3_cb+50*s_sq*b3
              terms(2)=120*r*s_sq-60*s_sq*b3
              terms(3)=-60*s*b3_sq+240*r*s*b3-60*r_sq*s
            case (12)
              terms(1)=-4*s+2-4*r+6*r_sq*s-58*r_sq*b3+2*r_cb-54*r*s_sq &
                       +2*s_sq*b3+2*s_cb+58*r*b3_sq-2*b3_cb-2*s*b3_sq &
                       +4*r*s*b3
              terms(2)=-120*r*s*b3+60*r_sq*s+60*r*s_sq
              terms(3)=-60*r*b3_sq+60*r_sq*b3+120*r*s*b3-60*r_sq*s
            case (22)
              terms(1)=-4*r+4*r_cb-52*r_sq*s+8*r*s*b3+4*r*s_sq-52*r_sq*b3+4*r*b3_sq
              terms(2)=-60*r_sq*b3+120*r_sq*s
              terms(3)=120*r_sq*b3-60*r_sq*s
            case default
              write(*,*) "invalid input order for function eval_p1"
              call abort
            end select
        end subroutine eval_p1

        subroutine eval_p2 (r,s,order, terms )
          real, intent(in) :: r,s
          integer, intent(in) :: order
          real, dimension(:), intent(out) :: terms
          real :: r_loc, s_loc
          r_loc=s
          s_loc=r

          select case (order)
            case (0)
              call eval_p1(r_loc, s_loc,0, terms)
            case (1)
              call eval_p1(r_loc, s_loc,2, terms)
            case (2)
              call eval_p1(r_loc, s_loc,1, terms)
            case (11)
              call eval_p1(r_loc, s_loc,22, terms)
            case (12)
              call eval_p1(r_loc, s_loc,12, terms)
            case (22)
              call eval_p1(r_loc, s_loc,11, terms)
            case default
              write(*,*) "invalid input order for function eval_p2"
              call abort
            end select
        end subroutine eval_p2

        subroutine eval_p3 (r,s,order, terms )
          real, intent(in) :: r,s
          integer, intent(in) :: order
          real, dimension(:), intent(out) :: terms
          real :: r_loc, s_loc
          real, dimension(3) :: terms_tmp,terms_tmp2
          r_loc=1-r-s
          s_loc=s

          select case (order)
            case (0)
              call eval_p1(r_loc, s_loc,0, terms)
            case (1)
              call eval_p1(r_loc, s_loc,1, terms)
              terms(1)=-terms(1)
              terms(2)=-terms(2)
              terms(3)=-terms(3)
            case (2)
              call eval_p1(r_loc, s_loc,1, terms)
              call eval_p1(r_loc, s_loc,2, terms_tmp)
              terms(1)=-terms(1)+terms_tmp(1)
              terms(2)=-terms(2)+terms_tmp(2)
              terms(3)=-terms(3)+terms_tmp(3)
            case (11)
              call eval_p1(r_loc, s_loc,11, terms)
            case (12)
              call eval_p1(r_loc, s_loc,11, terms)
              call eval_p1(r_loc, s_loc,12, terms_tmp)
              terms(1)=terms(1)-terms_tmp(1)
              terms(2)=terms(2)-terms_tmp(2)
              terms(3)=terms(3)-terms_tmp(3)
            case (22)
              call eval_p1(r_loc, s_loc,11, terms)
              call eval_p1(r_loc, s_loc,12, terms_tmp)
              call eval_p1(r_loc, s_loc,22, terms_tmp2)
              terms(1)=terms(1)-2*terms_tmp(1)+terms_tmp2(1)
              terms(2)=terms(2)-2*terms_tmp(2)+terms_tmp2(2)
              terms(3)=terms(3)-2*terms_tmp(3)+terms_tmp2(3)
            case default
              write(*,*) "invalid input order for function eval_p3"
              call abort
            end select
        end subroutine eval_p3
        ! first derivative funtion at vertex 1
        ! order: 0 value; 1 partial r; 2 partial s; 11, r^2; 12, rs, 22, ss
        ! output terms(1)+lambda21*(terms(2))
        subroutine eval_t1_2(r,s,order,terms)
        real, intent(in) :: r,s
          integer, intent(in) :: order
          real, dimension(:), intent(out) :: terms
          real :: r_sq, r_cb !r^2, r^3
          real :: s_sq, s_cb !s^2, s^3,
          real :: b3, b3_sq, b3_cb
          b3=1-r-s
          b3_sq=b3*b3
          b3_cb=b3_sq*b3
          r_sq=r*r
          r_cb=r_sq*r
          s_sq=s*s
          s_cb=s_sq*s
          select case (order)
            case (0)
              terms(1)=r_sq*s+2*r_cb*s_sq+2*r_cb*s*b3-r_sq*s_cb+10*r_sq*s_sq*b3+2*r_sq*s*b3_sq
              terms(2)=-15*r_sq*s_sq*b3
            case (1)
              terms(1)=2*r*s-4*r_sq*s_sq+2*r_sq*s*b3-2*r_cb*s-2*r*s_cb+20*r*s_sq*b3+4*r*s*b3_sq
              terms(2)=-30*r*s_sq*b3+15*r_sq*s_sq
            case (2)
              terms(1)=r_sq+2*r_cb*s+2*r_cb*b3-13*r_sq*s_sq+16*r_sq*s*b3+2*r_sq*b3_sq
              terms(2)=-30*r_sq*s*b3+15*r_sq*s_sq
            case (11)
              terms(1)=2*s-28*r*s_sq-4*r*s*b3-8*r_sq*s-2*s_cb+20*s_sq*b3+4*s*b3_sq
              terms(2)=-30*s_sq*b3+60*r*s_sq
            case (12)
              terms(1)=2*r-10*r_sq*s+2*r_sq*b3-2*r_cb-26*r*s_sq+32*r*s*b3+4*r*b3_sq
              terms(2)=-60*r*s*b3+30*r_sq*s+30*r*s_sq
            case (22)
              terms(1)=-42*r_sq*s+12*r_sq*b3
              terms(2)=-30*r_sq*b3+60*r_sq*s
            case default
              write(*,*) "invalid input order for function eval_t1_2"
              call abort
            end select
        end subroutine eval_t1_2

        subroutine eval_t2_1 (r,s,order, terms )
          real, intent(in) :: r,s
          integer, intent(in) :: order
          real, dimension(:), intent(out) :: terms
          real :: r_loc, s_loc
          r_loc=s
          s_loc=r

          select case (order)
            case (0)
              call eval_t1_2(r_loc, s_loc,0, terms)
            case (1)
              call eval_t1_2(r_loc, s_loc,2, terms)
            case (2)
              call eval_t1_2(r_loc, s_loc,1, terms)
            case (11)
              call eval_t1_2(r_loc, s_loc,22, terms)
            case (12)
              call eval_t1_2(r_loc, s_loc,12, terms)
            case (22)
              call eval_t1_2(r_loc, s_loc,11, terms)
            case default
              write(*,*) "invalid input order for function eval_t2_1"
              call abort
            end select
        end subroutine eval_t2_1

        subroutine eval_t1_3 (r,s,order, terms )
          real, intent(in) :: r,s
          integer, intent(in) :: order
          real, dimension(:), intent(out) :: terms
          real :: r_loc, s_loc
          real, dimension(3) :: terms_tmp,terms_tmp2
          r_loc=r
          s_loc=1-r-s

          select case (order)
            case (0)
              call eval_t1_2(r_loc, s_loc,0, terms)
            case (1)
              call eval_t1_2(r_loc, s_loc,1, terms)
              call eval_t1_2(r_loc, s_loc,2, terms_tmp)
              terms(1)=terms(1)-terms_tmp(1)
              terms(2)=terms(2)-terms_tmp(2)
            case (2)
              call eval_t1_2(r_loc, s_loc,2, terms_tmp)
              terms(1)=-terms_tmp(1)
              terms(2)=-terms_tmp(2)
            case (11)
              call eval_t1_2(r_loc, s_loc,11, terms)
              call eval_t1_2(r_loc, s_loc,12, terms_tmp)
              call eval_t1_2(r_loc, s_loc,22, terms_tmp2)
              terms(1)=terms(1)-2*terms_tmp(1)+terms_tmp2(1)
              terms(2)=terms(2)-2*terms_tmp(2)+terms_tmp2(2)
            case (12)
              call eval_t1_2(r_loc, s_loc,12, terms_tmp)
              call eval_t1_2(r_loc, s_loc,22, terms_tmp2)
              terms(1)=-terms_tmp(1)+terms_tmp2(1)
              terms(2)=-terms_tmp(2)+terms_tmp2(2)
            case (22)
              call eval_t1_2(r_loc, s_loc,22, terms)
            case default
              write(*,*) "invalid input order for function eval_t1_3"
              call abort
            end select
        end subroutine eval_t1_3


        subroutine eval_t2_3 (r,s,order, terms )
          real, intent(in) :: r,s
          integer, intent(in) :: order
          real, dimension(:), intent(out) :: terms
          real :: r_loc, s_loc
          real, dimension(3) :: terms_tmp,terms_tmp2
          r_loc=s
          s_loc=1-r-s

          select case (order)
            case (0)
              call eval_t1_2(r_loc, s_loc,0, terms)
            case (1)
              call eval_t1_2(r_loc, s_loc,2, terms_tmp)
              terms(1)=-terms_tmp(1)
              terms(2)=-terms_tmp(2)
            case (2)
              call eval_t1_2(r_loc, s_loc,1, terms)
              call eval_t1_2(r_loc, s_loc,2, terms_tmp)
              terms(1)=terms(1)-terms_tmp(1)
              terms(2)=terms(2)-terms_tmp(2)
            case (11)
              call eval_t1_2(r_loc, s_loc,22, terms)
            case (12)
              call eval_t1_2(r_loc, s_loc,12, terms_tmp)
              call eval_t1_2(r_loc, s_loc,22, terms_tmp2)
              terms(1)=-terms_tmp(1)+terms_tmp2(1)
              terms(2)=-terms_tmp(2)+terms_tmp2(2)
            case (22)
              call eval_t1_2(r_loc, s_loc,11, terms)
              call eval_t1_2(r_loc, s_loc,12, terms_tmp)
              call eval_t1_2(r_loc, s_loc,22, terms_tmp2)
              terms(1)=terms(1)-2*terms_tmp(1)+terms_tmp2(1)
              terms(2)=terms(2)-2*terms_tmp(2)+terms_tmp2(2)
            case default
              write(*,*) "invalid input order for function eval_t2_3"
              call abort
            end select
        end subroutine eval_t2_3

        subroutine eval_t3_1 (r,s,order, terms )
          real, intent(in) :: r,s
          integer, intent(in) :: order
          real, dimension(:), intent(out) :: terms
          real :: r_loc, s_loc
          real, dimension(3) :: terms_tmp,terms_tmp2
          r_loc=1-r-s
          s_loc=r

          select case (order)
            case (0)
              call eval_t1_2(r_loc, s_loc,0, terms)
            case (1)
              call eval_t1_2(r_loc, s_loc,1, terms)
              call eval_t1_2(r_loc, s_loc,2, terms_tmp)
              terms(1)=-terms(1)+terms_tmp(1)
              terms(2)=-terms(2)+terms_tmp(2)
            case (2)
              call eval_t1_2(r_loc, s_loc,1, terms)
              terms(1)=-terms(1)
              terms(2)=-terms(2)
            case (11)
              call eval_t1_2(r_loc, s_loc,11, terms)
              call eval_t1_2(r_loc, s_loc,12, terms_tmp)
              call eval_t1_2(r_loc, s_loc,22, terms_tmp2)
              terms(1)=terms(1)-2*terms_tmp(1)+terms_tmp2(1)
              terms(2)=terms(2)-2*terms_tmp(2)+terms_tmp2(2)
            case (12)
              call eval_t1_2(r_loc, s_loc,11, terms)
              call eval_t1_2(r_loc, s_loc,12, terms_tmp)
              !call eval_t1_2(r_loc, s_loc,22, terms_tmp2)
              terms(1)=terms(1)-terms_tmp(1)
              terms(2)=terms(2)-terms_tmp(2)
            case (22)
              call eval_t1_2(r_loc, s_loc,11, terms)
            case default
              write(*,*) "invalid input order for function eval_t3_1"
              call abort
            end select
        end subroutine eval_t3_1

        subroutine eval_t3_2 (r,s,order, terms )
          real, intent(in) :: r,s
          integer, intent(in) :: order
          real, dimension(:), intent(out) :: terms
          real :: r_loc, s_loc
          real, dimension(3) :: terms_tmp,terms_tmp2
          r_loc=1-r-s
          s_loc=s

          select case (order)
            case (0)
              call eval_t1_2(r_loc, s_loc,0, terms)
            case (1)
              call eval_t1_2(r_loc, s_loc,1, terms)
              terms(1)=-terms(1)
              terms(2)=-terms(2)
            case (2)
              call eval_t1_2(r_loc, s_loc,1, terms)
              call eval_t1_2(r_loc, s_loc,2, terms_tmp)
              terms(1)=-terms(1)+terms_tmp(1)
              terms(2)=-terms(2)+terms_tmp(2)
            case (11)
              call eval_t1_2(r_loc, s_loc,11, terms)
            case (12)
              call eval_t1_2(r_loc, s_loc,11, terms)
              call eval_t1_2(r_loc, s_loc,12, terms_tmp)
              terms(1)=terms(1)-terms_tmp(1)
              terms(2)=terms(2)-terms_tmp(2)
            case (22)
              call eval_t1_2(r_loc, s_loc,11, terms)
              call eval_t1_2(r_loc, s_loc,12, terms_tmp)
              call eval_t1_2(r_loc, s_loc,22, terms_tmp2)
              terms(1)=terms(1)-2*terms_tmp(1)+terms_tmp2(1)
              terms(2)=terms(2)-2*terms_tmp(2)+terms_tmp2(2)
            case default
              write(*,*) "invalid input order for function eval_t3_2"
              call abort
            end select
        end subroutine eval_t3_2

        ! second derivative funtion at vertex 1
        ! order: 0 value; 1 partial r; 2 partial s; 11, r^2; 12, rs, 22, ss
        ! output terms(1)+lambda21*(terms(2))
        subroutine eval_t1_22(r,s,order,terms)
        real, intent(in) :: r,s
          integer, intent(in) :: order
          real, dimension(:), intent(out) :: terms
          real :: r_sq, r_cb !r^2, r^3
          real :: s_sq, s_cb !s^2, s^3,
          real :: b3, b3_sq, b3_cb
          b3=1-r-s
          b3_sq=b3*b3
          b3_cb=b3_sq*b3
          r_sq=r*r
          r_cb=r_sq*r
          s_sq=s*s
          s_cb=s_sq*s

          terms=0.
          select case (order)
            case (0)
              terms(1)=1./2.*r_cb*s_sq+3./2.*r_sq*s_sq*b3
              terms(2)=-5./2.*r_sq*s_sq*b3
            case (1)
              terms(1)=3*r*s_sq*b3
              terms(2)=-5*r*s_sq*b3+5./2*r_sq*s_sq
            case (2)
              terms(1)=r_cb*s+3*r_sq*s*b3-3./2.*r_sq*s_sq
              terms(2)=-5*r_sq*s*b3+5./2.*r_sq*s_sq
            case (11)
              terms(1)=3*s_sq*b3-3*r*s_sq
              terms(2)=-5*s_sq*b3+10*r*s_sq
            case (12)
              terms(1)=6*r*s*b3-3*r*s_sq
              terms(2)=-10*r*s*b3+5*r*s_sq+5*r_sq*s
            case (22)
              terms(1)=r_cb+3*r_sq*b3-6*r_sq*s
              terms(2)=-5*r_sq*b3+10*r_sq*s
            case default
              write(*,*) "invalid input order for function eval_t1_22"
              call abort
            end select
        end subroutine eval_t1_22


        subroutine eval_t1_33 (r,s,order, terms )
          real, intent(in) :: r,s
          integer, intent(in) :: order
          real, dimension(:), intent(out) :: terms
          real :: r_loc, s_loc
          real, dimension(3) :: terms_tmp,terms_tmp2
          r_loc=r
          s_loc=1-r-s

          select case (order)
            case (0)
              call eval_t1_22(r_loc, s_loc,0, terms)
            case (1)
              call eval_t1_22(r_loc, s_loc,1, terms)
              call eval_t1_22(r_loc, s_loc,2, terms_tmp)
              terms(1)=terms(1)-terms_tmp(1)
              terms(2)=terms(2)-terms_tmp(2)
            case (2)
              call eval_t1_22(r_loc, s_loc,2, terms_tmp)
              terms(1)=-terms_tmp(1)
              terms(2)=-terms_tmp(2)
            case (11)
              call eval_t1_22(r_loc, s_loc,11, terms)
              call eval_t1_22(r_loc, s_loc,12, terms_tmp)
              call eval_t1_22(r_loc, s_loc,22, terms_tmp2)
              terms(1)=terms(1)-2*terms_tmp(1)+terms_tmp2(1)
              terms(2)=terms(2)-2*terms_tmp(2)+terms_tmp2(2)
            case (12)
              call eval_t1_22(r_loc, s_loc,12, terms_tmp)
              call eval_t1_22(r_loc, s_loc,22, terms_tmp2)
              terms(1)=-terms_tmp(1)+terms_tmp2(1)
              terms(2)=-terms_tmp(2)+terms_tmp2(2)
            case (22)
              call eval_t1_22(r_loc, s_loc,22, terms)
            case default
              write(*,*) "invalid input order for function eval_t1_33"
              call abort
            end select
        end subroutine eval_t1_33

        subroutine eval_t2_11 (r,s,order, terms )
          real, intent(in) :: r,s
          integer, intent(in) :: order
          real, dimension(:), intent(out) :: terms
          real :: r_loc, s_loc
          real, dimension(3) :: terms_tmp,terms_tmp2
          r_loc=s
          s_loc=r

          select case (order)
            case (0)
              call eval_t1_22(r_loc, s_loc,0, terms)
            case (1)
              call eval_t1_22(r_loc, s_loc,2, terms)
            case (2)
              call eval_t1_22(r_loc, s_loc,1, terms)
            case (11)
              call eval_t1_22(r_loc, s_loc,22, terms)
            case (12)
              call eval_t1_22(r_loc, s_loc,12, terms)
            case (22)
              call eval_t1_22(r_loc, s_loc,11, terms)
            case default
              write(*,*) "invalid input order for function eval_t2_11"
              call abort
            end select
        end subroutine eval_t2_11

        subroutine eval_t2_33 (r,s,order, terms )
          real, intent(in) :: r,s
          integer, intent(in) :: order
          real, dimension(:), intent(out) :: terms
          real :: r_loc, s_loc
          real, dimension(3) :: terms_tmp,terms_tmp2
          r_loc=s
          s_loc=1-r-s

          select case (order)
            case (0)
              call eval_t1_22(r_loc, s_loc,0, terms)
            case (1)
              call eval_t1_22(r_loc, s_loc,2, terms_tmp)
              terms(1)=-terms_tmp(1)
              terms(2)=-terms_tmp(2)
            case (2)
              call eval_t1_22(r_loc, s_loc,1, terms)
              call eval_t1_22(r_loc, s_loc,2, terms_tmp)
              terms(1)=terms(1)-terms_tmp(1)
              terms(2)=terms(2)-terms_tmp(2)
            case (11)
              call eval_t1_22(r_loc, s_loc,22, terms)
            case (12)
              call eval_t1_22(r_loc, s_loc,12, terms_tmp)
              call eval_t1_22(r_loc, s_loc,22, terms_tmp2)
              terms(1)=-terms_tmp(1)+terms_tmp2(1)
              terms(2)=-terms_tmp(2)+terms_tmp2(2)
            case (22)
              call eval_t1_22(r_loc, s_loc,11, terms)
              call eval_t1_22(r_loc, s_loc,12, terms_tmp)
              call eval_t1_22(r_loc, s_loc,22, terms_tmp2)
              terms(1)=terms(1)-2*terms_tmp(1)+terms_tmp2(1)
              terms(2)=terms(2)-2*terms_tmp(2)+terms_tmp2(2)
            case default
              write(*,*) "invalid input order for function eval_t2_33"
              call abort
            end select
        end subroutine eval_t2_33

        subroutine eval_t3_11 (r,s,order, terms )
          real, intent(in) :: r,s
          integer, intent(in) :: order
          real, dimension(:), intent(out) :: terms
          real :: r_loc, s_loc
          real, dimension(3) :: terms_tmp,terms_tmp2
          r_loc=1-r-s
          s_loc=r

          select case (order)
            case (0)
              call eval_t1_22(r_loc, s_loc,0, terms)
            case (1)
              call eval_t1_22(r_loc, s_loc,1, terms)
              call eval_t1_22(r_loc, s_loc,2, terms_tmp)
              terms(1)=-terms(1)+terms_tmp(1)
              terms(2)=-terms(2)+terms_tmp(2)
            case (2)
              call eval_t1_22(r_loc, s_loc,1, terms)
              terms(1)=-terms(1)
              terms(2)=-terms(2)
            case (11)
              call eval_t1_22(r_loc, s_loc,11, terms)
              call eval_t1_22(r_loc, s_loc,12, terms_tmp)
              call eval_t1_22(r_loc, s_loc,22, terms_tmp2)
              terms(1)=terms(1)-2*terms_tmp(1)+terms_tmp2(1)
              terms(2)=terms(2)-2*terms_tmp(2)+terms_tmp2(2)
            case (12)
              call eval_t1_22(r_loc, s_loc,11, terms)
              call eval_t1_22(r_loc, s_loc,12, terms_tmp)
              terms(1)=terms(1)-terms_tmp(1)
              terms(2)=terms(2)-terms_tmp(2)
            case (22)
              call eval_t1_22(r_loc, s_loc,11, terms)
            case default
              write(*,*) "invalid input order for function eval_t3_11"
              call abort
            end select
        end subroutine eval_t3_11

        subroutine eval_t3_22 (r,s,order, terms )
          real, intent(in) :: r,s
          integer, intent(in) :: order
          real, dimension(:), intent(out) :: terms
          real :: r_loc, s_loc
          real, dimension(3) :: terms_tmp,terms_tmp2
          r_loc=1-r-s
          s_loc=s

          select case (order)
            case (0)
              call eval_t1_22(r_loc, s_loc,0, terms)
            case (1)
              call eval_t1_22(r_loc, s_loc,1, terms)
              terms(1)=-terms(1)
              terms(2)=-terms(2)
            case (2)
              call eval_t1_22(r_loc, s_loc,1, terms)
              call eval_t1_22(r_loc, s_loc,2, terms_tmp)
              terms(1)=-terms(1)+terms_tmp(1)
              terms(2)=-terms(2)+terms_tmp(2)
            case (11)
              call eval_t1_22(r_loc, s_loc,11, terms)
            case (12)
              call eval_t1_22(r_loc, s_loc,11, terms)
              call eval_t1_22(r_loc, s_loc,12, terms_tmp)
              terms(1)=terms(1)-terms_tmp(1)
              terms(2)=terms(2)-terms_tmp(2)
            case (22)
              call eval_t1_22(r_loc, s_loc,11, terms)
              call eval_t1_22(r_loc, s_loc,12, terms_tmp)
              call eval_t1_22(r_loc, s_loc,22, terms_tmp2)
              terms(1)=terms(1)-2*terms_tmp(1)+terms_tmp2(1)
              terms(2)=terms(2)-2*terms_tmp(2)+terms_tmp2(2)
            case default
              write(*,*) "invalid input order for function eval_t3_22"
              call abort
            end select
        end subroutine eval_t3_22

        ! twisted derivative funtion at vertex 1
        ! order: 0 value; 1 partial r; 2 partial s; 11, r^2; 12, rs, 22, ss
        subroutine eval_t1_23(r,s,order,terms)
        real, intent(in) :: r,s
          integer, intent(in) :: order
          real, dimension(:),intent(out) :: terms
          real :: r_sq, r_cb !r^2, r^3
          real :: s_sq, s_cb !s^2, s^3,
          real :: b3, b3_sq, b3_cb
          b3=1-r-s
          b3_sq=b3*b3
          b3_cb=b3_sq*b3
          r_sq=r*r
          r_cb=r_sq*r
          s_sq=s*s
          s_cb=s_sq*s

          terms=0.
          select case (order)
            case (0)
              terms(1)=r_sq*s*b3
            case (1)
              terms(1)=2*r*s*b3-r_sq*s
            case (2)
              terms(1)=r_sq*b3-r_sq*s
            case (11)
              terms(1)=2*b3*s-4*r*s
            case (12)
              terms(1)=2*r*b3-2*r*s-r_sq
            case (22)
              terms(1)=-2*r_sq
            case default
              write(*,*) "invalid input order for function eval_t1_23"
              call abort
          end select
        end subroutine eval_t1_23

        subroutine eval_t2_13 (r,s,order, terms )
          real, intent(in) :: r,s
          integer, intent(in) :: order
          real, dimension(:), intent(out) :: terms
          real :: r_loc, s_loc
          real, dimension(3) :: terms_tmp,terms_tmp2
          r_loc=s
          s_loc=r

          select case (order)
            case (0)
              call eval_t1_23(r_loc, s_loc,0, terms)
            case (1)
              call eval_t1_23(r_loc, s_loc,2, terms)
            case (2)
              call eval_t1_23(r_loc, s_loc,1, terms)
            case (11)
              call eval_t1_23(r_loc, s_loc,22, terms)
            case (12)
              call eval_t1_23(r_loc, s_loc,12, terms)
            case (22)
              call eval_t1_23(r_loc, s_loc,11, terms)
            case default
              write(*,*) "invalid input order for function eval_t2_13"
              call abort
            end select
        end subroutine eval_t2_13

        subroutine eval_t3_12 (r,s,order, terms )
          real, intent(in) :: r,s
          integer, intent(in) :: order
          real, dimension(:), intent(out) :: terms
          real :: r_loc, s_loc
          real, dimension(3) :: terms_tmp,terms_tmp2
          r_loc=1-r-s
          s_loc=r

          select case (order)
            case (0)
              call eval_t1_23(r_loc, s_loc,0, terms)
            case (1)
              call eval_t1_23(r_loc, s_loc,1, terms)
              call eval_t1_23(r_loc, s_loc,2, terms_tmp)
              terms(1)=-terms(1)+terms_tmp(1)
            case (2)
              call eval_t1_23(r_loc, s_loc,1, terms)
              terms(1)=-terms(1)
            case (11)
              call eval_t1_23(r_loc, s_loc,11, terms)
              call eval_t1_23(r_loc, s_loc,12, terms_tmp)
              call eval_t1_23(r_loc, s_loc,22, terms_tmp2)
              terms(1)=terms(1)-2*terms_tmp(1)+terms_tmp2(1)
            case (12)
              call eval_t1_23(r_loc, s_loc,11, terms)
              call eval_t1_23(r_loc, s_loc,12, terms_tmp)
              terms(1)=terms(1)-terms_tmp(1)
            case (22)
              call eval_t1_23(r_loc, s_loc,11, terms)
            case default
              write(*,*) "invalid input order for function eval_t3_12"
              call abort
            end select
        end subroutine eval_t3_12


        ! evaluate nth shape funtion terms at (r,s)
        ! order: 0 value; 1 partial r; 2 partial s; 11, r^2; 12, rs, 22, ss
        ! r=b1, s=b2 
        subroutine eval_basis_function_terms( r, s, nth, order, terms)
          integer, intent(in) :: nth, order
          real, intent(in) :: r, s
          real, dimension(:),intent(out) :: terms

          select case (nth)
            case (1)
              call eval_p1(r,s,order, terms)  
            case (2) 
              call eval_t1_2(r,s,order, terms)
            case (3) 
              call eval_t1_3(r,s,order, terms)
            case (4)
              call eval_t1_22(r,s,order, terms)
            case (5)
              call eval_t1_23(r,s,order, terms)
            case (6) 
              call eval_t1_33(r,s,order, terms)
            case (7)
              call eval_p2(r,s,order, terms)
            case (8)
              call eval_t2_3(r,s,order, terms)
            case (9)
              call eval_t2_1(r,s,order, terms)
            case (10)
              call eval_t2_33(r,s,order, terms)
            case (11)
              call eval_t2_13(r,s,order, terms)
            case (12)
              call eval_t2_11(r,s,order, terms)
            case (13)
              call eval_p3(r,s,order, terms)
            case (14)
              call eval_t3_1(r,s,order, terms)
            case (15)
              call eval_t3_2(r,s,order, terms)
            case (16)
              call eval_t3_11(r,s,order, terms)
            case (17)
              call eval_t3_12(r,s,order, terms)
            case (18)
              call eval_t3_22(r,s,order, terms)
            case default
              write(*,*) "invalid input order for eval_basis_function"
              call abort
          end select
        end subroutine eval_basis_function_terms

        ! evaluate nth shape funtion at (r,s)
        ! not recommended for pre caculation
        ! order: 0 value; 1 partial r; 2 partial s; 11, r^2; 12, rs, 22, ss
        ! r=b1, s=b2 
        subroutine eval_basis_function( r, s, nth, order, value)
          integer, intent(in) :: nth, order
          real, intent(in) :: r, s
          real, intent(out) :: value
          real, dimension(3) :: terms

          select case (nth)
            case (1)
              call eval_p1(r,s,order, terms)  
              value=terms(1)+lambda(2,1)*terms(2)+lambda(3,1)*terms(3)
            case (2) 
              call eval_t1_2(r,s,order, terms)
              value=terms(1)+lambda(2,1)*terms(2)
            case (3) 
              call eval_t1_3(r,s,order, terms)
              value=terms(1)+lambda(3,1)*terms(2)
            case (4)
              call eval_t1_22(r,s,order, terms)
              value=terms(1)+lambda(2,1)*terms(2)
            case (5)
              call eval_t1_23(r,s,order, terms)
              value=terms(1)
            case (6) 
              call eval_t1_33(r,s,order, terms)
              value=terms(1)+lambda(3,1)*terms(2)
            case (7)
              call eval_p2(r,s,order, terms)
              value=terms(1)+lambda(1,2)*terms(2)+lambda(3,2)*terms(3)
            case (8)
              call eval_t2_3(r,s,order, terms)
              value=terms(1)+lambda(3,2)*terms(2)
            case (9)
              call eval_t2_1(r,s,order, terms)
              value=terms(1)+lambda(1,2)*terms(2)
            case (10)
              call eval_t2_33(r,s,order, terms)
              value=terms(1)+lambda(3,2)*terms(2)
            case (11)
              call eval_t2_13(r,s,order, terms)
              value=terms(1)
            case (12)
              call eval_t2_11(r,s,order, terms)
              value=terms(1)+lambda(1,2)*terms(2)
            case (13)
              call eval_p3(r,s,order, terms)
              value=terms(1)+lambda(2,3)*terms(2)+lambda(1,3)*terms(3)
            case (14)
              call eval_t3_1(r,s,order, terms)
              value=terms(1)+lambda(1,3)*terms(2)
            case (15)
              call eval_t3_2(r,s,order, terms)
              value=terms(1)+lambda(2,3)*terms(2)
            case (16)
              call eval_t3_11(r,s,order, terms)
              value=terms(1)+lambda(1,3)*terms(2)
            case (17)
              call eval_t3_12(r,s,order, terms)
              value=terms(1)
            case (18)
              call eval_t3_22(r,s,order, terms)
              value=terms(1)+lambda(2,3)*terms(2)
            case default
              write(*,*) "invalid input order for eval_basis_function"
              call abort
          end select
        end subroutine eval_basis_function


        ! evalutate field at local coord (r,s)
        ! dofs: local
        ! order: global (x,y)  0 value; 1 partial x; 2 partial y; 11, xx; 12, xy, 22, yy
        ! not recommended for precaculation
        subroutine eval_field( r, s, dofs, values)
          real, intent(in) :: r,s
          real, dimension(18), intent(in) :: dofs
          real, dimension(6), intent(out) :: values
          real :: values_tmp,values_tmp2,values_tmp3
          integer :: counter
     
          values=0;
          do counter=1,18
             call eval_basis_function(r,s,counter,0,values_tmp)
             !write(*,*) counter,values_tmp
             values(1)=values(1)+dofs(counter)*values_tmp
          end do
          do counter=1,18
             call eval_basis_function(r,s,counter,1,values_tmp)
             call eval_basis_function(r,s,counter,2,values_tmp2)
             values(2)=values(2)+dofs(counter)*(values_tmp*dedge_dx(1,3)+values_tmp2*dedge_dx(2,3))
             values(3)=values(3)+dofs(counter)*(values_tmp*dedge_dy(1,3)+values_tmp2*dedge_dy(2,3))
          end do
             !write(*,*) dedge_dx(1,3),dedge_dx(2,3),dedge_dy(1,3),dedge_dy(2,3)
          do counter=1,18
             call eval_basis_function(r,s,counter,11,values_tmp)
             call eval_basis_function(r,s,counter,12,values_tmp2)
             call eval_basis_function(r,s,counter,22,values_tmp3)
             values(4)=values(4)+dofs(counter)*(values_tmp*dedge_dx(1,3)*dedge_dx(1,3) &
                       +2*values_tmp2*dedge_dx(2,3)*dedge_dx(1,3)  &
                       +values_tmp3*dedge_dx(2,3)*dedge_dx(2,3))
             values(5)=values(5)+dofs(counter)*(values_tmp*dedge_dx(1,3)*dedge_dy(1,3) &
                       +values_tmp2*(dedge_dx(1,3)*dedge_dy(2,3)+dedge_dx(2,3)*dedge_dy(1,3))  &
                       +values_tmp3*dedge_dx(2,3)*dedge_dy(2,3)) 
             values(6)=values(6)+dofs(counter)*(values_tmp*dedge_dy(1,3)*dedge_dy(1,3) &
                       +2*values_tmp2*dedge_dy(2,3)*dedge_dy(1,3)  &
                       +values_tmp3*dedge_dy(2,3)*dedge_dy(2,3))
          end do
        end subroutine eval_field
      end module explicit_element 
