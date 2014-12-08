module adapt
  real :: ke_previous
  contains
  subroutine adapt_by_psi
    use basic
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

    vectype, dimension(dofs_per_node) :: dat
    integer :: izone, izonedim, inode(nodes_per_element)
    integer :: numelms, itri, num_adj_ent

    call create_field(temporary_field)
    if(eqsubtract.eq.1) then
       temporary_field = psi_field(0)
       call add(temporary_field, psi_field(1))
    else
       temporary_field = psi_field(1)
    end if
    if(icsubtract.eq.1) call add(temporary_field, psi_coil_field)
    call straighten_field(temporary_field)
    if(imulti_region.eq.1 .and. adapt_psin_vacuum.ne.0.) then
       dat = 0.
       dat(1) = (psibound - psimin)*adapt_psin_vacuum + psimin
       numelms = local_elements()
       do itri=1, numelms
          !call zonfac(itri,izone,izonedim)
          call m3dc1_ent_getgeomclass(2, itri-1, izonedim, izone)
          if(izone.ge.3) then
             !call nodfac(itri,inode)
             call m3dc1_ent_getadj (2, itri-1, 0, inode, 3, num_adj_ent)
             inode = inode+1
             do i=1,3
                call set_node_data(temporary_field,inode(i),dat)
             end do
         end if
       end do
       call sum_shared(temporary_field%vec)
    end if
    call adapt_by_field(temporary_field%vec%id,psimin,psibound)
    call destroy_field(temporary_field)
    call space(0)
    call update_nodes_owned()
    call tridef
    call create_newvar_matrices
    call derived_quantities(1)
    ke_previous = ekin
  end subroutine adapt_by_psi

  integer function run_adapt ()
    use diagnostics
    use basic
    !print *, "check run_adaptA: ke_previous, ekin, ntime",ke_previous,ekin,ntime
    run_adapt = 0
    !if(mod(ntime,1) .eq. 0) then
     ! run_adapt=1
      !ke_previous = ekin;
    !end if
  end function
end module adapt

