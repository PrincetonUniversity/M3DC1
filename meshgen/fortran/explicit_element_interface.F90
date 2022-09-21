      subroutine set_element_coord_i ( itri_p, coord)
        use explicit_element
        real, dimension(3,3), intent(in) :: coord ! three nodes, xyz
        integer , intent(in) :: itri_p
        call set_element_coord(itri_p, coord)
      end subroutine set_element_coord_i

      subroutine glb_xyz_2_loc_i ( xyz_glb, xyz_loc)
        use explicit_element
        real, dimension(3), intent(in) :: xyz_glb
        real, dimension(3), intent(out) :: xyz_loc
        call glb_xyz_2_loc( xyz_glb, xyz_loc)
      end subroutine glb_xyz_2_loc_i

      subroutine loc_xyz_2_glb_i ( xyz_loc, xyz_glb)
        use explicit_element
        real, dimension(3), intent(in) :: xyz_loc
        real, dimension(3), intent(out) :: xyz_glb
        call loc_xyz_2_glb ( xyz_loc, xyz_glb )
      end subroutine loc_xyz_2_glb_i

      subroutine glb_dof_2_loc_i ( dofs_glb, dofs_loc )
        use explicit_element
        real, dimension(18), intent(in):: dofs_glb
        real, dimension(18), intent(out):: dofs_loc
        call glb_dof_2_loc ( dofs_glb, dofs_loc )
      end subroutine  glb_dof_2_loc_i

      subroutine loc_dof_2_glb_i ( dofs_loc, dofs_glb )
        use explicit_element
        real, dimension(18), intent(in):: dofs_loc
        real, dimension(18), intent(out):: dofs_glb
        call loc_dof_2_glb( dofs_loc, dofs_glb )
      end subroutine loc_dof_2_glb_i 

      subroutine eval_field_i( r, s, dofs, values)
        use explicit_element
        real, intent(in) :: r,s
        real, dimension(18), intent(in) :: dofs
        real, dimension(6), intent(out) :: values
        call eval_field( r, s, dofs, values)
      end subroutine eval_field_i
      subroutine get_jacobi_i( jac)
        use explicit_element
        real, intent(out) :: jac
        jac=jacobi
      end subroutine get_jacobi_i

