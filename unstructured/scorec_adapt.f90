module scorec_adapt
  use vector_mod
  implicit none
  !adaptation control parameters
  integer :: iadapt_writevtk, iadapt_writesmb

  contains

!obsolete
#ifdef ADAPT_
subroutine adapt_mesh
    use basic
    use arrays
    use newvar_mod
    use transport_coefficients
    use diagnostics

    integer :: i, idx, node_dim=0
    integer :: numnodes, inode
    character(len=32) :: mesh_file_name
    integer, dimension(MAX_PTS) :: mr
    vectype, dimension(dofs_per_node) :: dofs1, dofs2
    type(field_type) :: size1, size2
    real, allocatable :: unit1(:)
    integer, allocatable :: node_ids(:)
    real :: x, phi, z, r, w
    real, parameter :: mins = 0.01, maxs = 0.2, aniso = 2.
    real, parameter :: x0 = 0.9, z0 =  0., r0=0.3, w0 = 0.2
    integer, parameter :: adapt_iters=10

    do i=1,adapt_iters
       if(myrank.eq.0) print *, "ADAPT TEST", i

       call create_field(size1)
       call create_field(size2)
       numnodes = local_nodes()
       allocate(unit1(3*numnodes))
       allocate(node_ids(numnodes))

       size1 = 0.
       size2 = 0.
       dofs1 = 0.
       dofs2 = 0.
       unit1 = 0.

       call m3dc1_ent_getlocalid(node_dim, node_ids,numnodes, numnodes)
 
       do idx=1, numnodes
          inode = node_ids(idx)+1  ! nodes_ids in c++ starts from 0
          call get_node_pos(inode, x, phi, z)

          r = sqrt((x - x0)**2 + (z - z0)**2)
          w = abs(r-r0)
          if(w > w0) then
             dofs1(1) = maxs
             dofs2(1) = maxs
             unit1((inode-1)*3+1) = 1.0
             unit1((inode-1)*3+2) = 0.0
             unit1((inode-1)*3+3) = 0.0
          else
             dofs1(1) = mins + (maxs - mins)*w/w0
             dofs2(1) = mins + aniso*(maxs - mins)*w/w0
             unit1((inode-1)*3+1) = (x-x0)/r
             unit1((inode-1)*3+2) = (z-z0)/r
             unit1((inode-1)*3+3) = 0.0
          end if
          
          call set_node_data(size1,inode,dofs1)
          call set_node_data(size2,inode,dofs2)

       end do
       call finalize(size1%vec)
       call finalize(size2%vec)

       call straighten_fields()
       call m3dc1_mesh_adapt(size1%vec%id, size2%vec%id, unit1)
       write(mesh_file_name,"(A7,A)") 'adapted', 0
       if(iadapt_writesmb .eq. 1) call m3dc1_mesh_write (mesh_file_name,1,i)
       if(iadapt_writevtk .eq. 1) call m3dc1_mesh_write (mesh_file_name,0,i)

       deallocate(unit1)
       deallocate(node_ids)

       call destroy_field(size1)
       call destroy_field(size2)
       call space(0)
       call update_nodes_owned()
       call reset_itris()
       call tridef
       call unstraighten_fields()

       call create_newvar_matrices
       field_vec = 0.
       field0_vec = 0.
       if (myrank .eq. 0) print *, "re-calculate equlibrium after adapt .."
       call initial_conditions
    end do

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

  end subroutine adapt_mesh
#endif

  subroutine straighten_fields ()
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

    integer :: ifield

    do ifield=0,1
       call straighten_field(u_field(ifield))
       call straighten_field(vz_field(ifield))
       call straighten_field(chi_field(ifield))

       call straighten_field(psi_field(ifield))
       call straighten_field(bz_field(ifield))
       call straighten_field(pe_field(ifield))

       call straighten_field(den_field(ifield))
       call straighten_field(p_field(ifield))
       call straighten_field(te_field(ifield))
       call straighten_field(ti_field(ifield))
    end do
  end subroutine straighten_fields

  subroutine unstraighten_fields ()
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

    integer :: ifield

    do ifield=0,1
       call unstraighten_field(u_field(ifield))
       call unstraighten_field(vz_field(ifield))
       call unstraighten_field(chi_field(ifield))

       call unstraighten_field(psi_field(ifield))
       call unstraighten_field(bz_field(ifield))
       call unstraighten_field(pe_field(ifield))

       call unstraighten_field(den_field(ifield))
       call unstraighten_field(p_field(ifield))
       call unstraighten_field(te_field(ifield))
       call unstraighten_field(ti_field(ifield))
    end do
  end subroutine unstraighten_fields
end module scorec_adapt
