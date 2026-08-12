module scorec_adapt
  use vector_mod
  implicit none
  !adaptation control parameters
  integer :: iadapt_writevtk, iadapt_writesmb

  contains

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
