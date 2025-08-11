#ifdef USECOMPLEX
#define CONJUGATE(x) (conjg(x))
#else
#define CONJUGATE(x) (x)
#endif

module auxiliary_fields
  use field

  implicit none

  type(field_type) :: bdotgradp
  type(field_type) :: bdotgradt
  type(field_type) :: torque_density_em
  type(field_type) :: torque_density_ntv
  type(field_type) :: chord_mask
  type(field_type) :: mag_reg
  type(field_type) :: ef_r, ef_phi, ef_z, ef_par, eta_j, psidot, veldif
  type(field_type) :: eta_jdb, bdgp, vlbdgp
  type(field_type) :: vdotgradt, adv1, adv2, adv3, vpar_field
  type(field_type) :: f1vplot,f1eplot,f2vplot,f2eplot,f3vplot,f3eplot,jdbobs
  type(field_type) :: deldotq_perp
  type(field_type) :: deldotq_par
  type(field_type) :: eta_jsq
  type(field_type) :: pot2_field
  type(field_type) :: mesh_zone
  type(field_type) :: z_effective
  type(field_type) :: kprad_totden
  type(field_type) :: Jp_BS_r,Jp_BS_z,Jp_BS_phi,JpdotB,JpdotB_dndpsi,JpdotB_dtedpsi,JpdotB_dtidpsi

  logical, private :: initialized = .false.

contains

subroutine create_auxiliary_fields
  use basic
  use kprad_m3dc1
  implicit none

if (ispradapt .eq. 1) then
  call create_field(bdotgradp, "bdotgradp")
  call create_field(bdotgradt, "bdotgradt")
  call create_field(torque_density_em, "torque_density_em")
  call create_field(torque_density_ntv, "torque_density_ntv")
  call create_field(chord_mask, "chord_mask")
  call create_field(mag_reg, "mag_reg")
  call create_field(ef_r, "ef_r")
  call create_field(ef_phi, "ef_phi")
  call create_field(ef_z, "ef_z")
  call create_field(ef_par, "ef_par")
  call create_field(eta_j, "eta_j")
  call create_field(mesh_zone, "mesh_zone")
  call create_field(z_effective, "z_effective")

  if(ikprad.ne.0) call create_field(kprad_totden, "kprad_totden")
  if(jadv.eq.0) then
     call create_field(psidot, "psidot")
     call create_field(veldif, "veldif")
     call create_field(eta_jdb, "eta_jdb")
     call create_field(bdgp, "bdgp")
     call create_field(vlbdgp, "vlbdgp")
  endif
  if(itemp_plot.eq.1) then
     call create_field(vdotgradt, "vdotgradt")
     call create_field(adv1, "adv1")
     call create_field(adv2, "adv2")
     call create_field(adv3, "adv3")
     call create_field(deldotq_perp, "deldotq_perp")
     call create_field(deldotq_par, "deldotq_par")
     call create_field(eta_jsq, "eta_jsq")
     call create_field(pot2_field, "pot2_field")
     call create_field(vpar_field, "vpar_field")
     call create_field(f1vplot, "f1vplot")
     call create_field(f1eplot, "f1eplot")
     call create_field(f2vplot, "f2vplot")
     call create_field(f2eplot, "f2eplot")
     call create_field(f3vplot, "f3vplot")
     call create_field(f3eplot, "f3eplot")
     call create_field(jdbobs, "jdbobs")
   endif
   if(ibootstrap.ne.0) then
    call create_field(Jp_BS_r, "Jp_BS_r")
    call create_field(Jp_BS_z, "Jp_BS_z")
    call create_field(Jp_BS_phi, "Jp_BS_phi")
    call create_field(JpdotB, "JpdotB")
    call create_field(JpdotB, "JpdotB_dndpsi")
    call create_field(JpdotB, "JpdotB_dtedpsi")
    call create_field(JpdotB, "JpdotB_dtidpsi")
   endif
else
  call create_field(bdotgradp)
  call create_field(bdotgradt)
  call create_field(torque_density_em)
  call create_field(torque_density_ntv)
  call create_field(chord_mask)
  call create_field(mag_reg)
  call create_field(ef_r)
  call create_field(ef_phi)
  call create_field(ef_z)
  call create_field(ef_par)
  call create_field(eta_j)
  call create_field(mesh_zone)
  call create_field(z_effective)
  if(ikprad.ne.0) call create_field(kprad_totden)
  if(jadv.eq.0) then
     call create_field(psidot)
     call create_field(veldif)
     call create_field(eta_jdb)
     call create_field(bdgp)
     call create_field(vlbdgp)
  endif
  if(itemp_plot.eq.1) then
     call create_field(vdotgradt)
     call create_field(adv1)
     call create_field(adv2)
     call create_field(adv3)
     call create_field(deldotq_perp)
     call create_field(deldotq_par)
     call create_field(eta_jsq)
  call create_field(pot2_field)
     call create_field(vpar_field)
     call create_field(f1vplot)
     call create_field(f1eplot)
     call create_field(f2vplot)
     call create_field(f2eplot)
     call create_field(f3vplot)
     call create_field(f3eplot)
     call create_field(jdbobs)
  endif
  if(ibootstrap.ne.0) then
    call create_field(Jp_BS_r)
    call create_field(Jp_BS_z)
    call create_field(Jp_BS_phi)
    call create_field(JpdotB)
    call create_field(JpdotB_dndpsi)
    call create_field(JpdotB_dtedpsi)
    call create_field(JpdotB_dtidpsi)
  endif
endif
  initialized = .true.
end subroutine create_auxiliary_fields

subroutine destroy_auxiliary_fields
  use basic
  use kprad_m3dc1
  implicit none

  if(.not.initialized) return
  call destroy_field(bdotgradp)
  call destroy_field(bdotgradt)
  call destroy_field(torque_density_em)
  call destroy_field(torque_density_ntv)
  call destroy_field(chord_mask)
  call destroy_field(mag_reg)
  call destroy_field(ef_r)
  call destroy_field(ef_phi)
  call destroy_field(ef_z)
  call destroy_field(ef_par)
  call destroy_field(eta_j)
  call destroy_field(mesh_zone)
  call destroy_field(z_effective)
  if(ikprad.ne.0) call destroy_field(kprad_totden)
  if(jadv.eq.0) then
     call destroy_field(psidot)
     call destroy_field(veldif)
     call destroy_field(eta_jdb)
     call destroy_field(bdgp)
     call destroy_field(vlbdgp)
  endif
  if(itemp_plot.eq.1) then
     call destroy_field(vdotgradt)
     call destroy_field(adv1)
     call destroy_field(adv2)
     call destroy_field(adv3)
     call destroy_field(deldotq_perp)
     call destroy_field(deldotq_par)
     call destroy_field(eta_jsq)
     call destroy_field(pot2_field)
     call destroy_field(vpar_field)
     call destroy_field(f1vplot)
     call destroy_field(f1eplot)
     call destroy_field(f2vplot)
     call destroy_field(f2eplot)
     call destroy_field(f3vplot)
     call destroy_field(f3eplot)
     call destroy_field(jdbobs)
  endif
  if(ibootstrap.ne.0) then
    call destroy_field(Jp_BS_r)
    call destroy_field(Jp_BS_z)
    call destroy_field(Jp_BS_phi)
    call destroy_field(JpdotB)
    call destroy_field(JpdotB_dndpsi)
    call destroy_field(JpdotB_dtidpsi)
    call destroy_field(JpdotB_dtedpsi)
  endif
end subroutine destroy_auxiliary_fields
  
subroutine calculate_pressures(ilin, pe, p, ne, nion, te, ti, ieqsub)
  use math
  use basic
  use m3dc1_nint
  use newvar_mod
  use diagnostics
  use metricterms_new
  use field
  use arrays
  use kprad_m3dc1

  implicit none

  type(field_type) :: pe, p, te, ti, ne, nion
  integer, intent(in) :: ilin    ! 0 = calculate only equilibrium part
  integer, intent(in) :: ieqsub  ! 1 = supplied fields are only perturbed parts

  integer :: numelms, itri, i

  type(field_type) :: pe_f, p_f

  vectype, dimension(dofs_per_element) :: dofs_p, dofs_pe

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Calculating pressures'

if (ispradapt .eq. 1) then 
  call create_field(pe_f, " pe_f")
  call create_field(p_f, "p_f")
else
  call create_field(pe_f)
  call create_field(p_f)
endif

  pe_f = 0.
  p_f = 0.

  numelms = local_elements()
  do itri=1,numelms
     call define_element_quadrature(itri, int_pts_aux, 5)
     call define_fields(itri, 0, 1, 0, ieqsub)

     call eval_ops(itri, ne, ne179, rfac)
     call eval_ops(itri, nion, n179, rfac)
     call eval_ops(itri, te, te179, rfac)
     call eval_ops(itri, ti, ti179, rfac)
     if(ieqsub.eq.1 .and. ilin.eq.1) then
        call eval_ops(itri, ne_field(0), ne079, rfac)
        call eval_ops(itri, den_field(0), n079, rfac)
        call eval_ops(itri, te_field(0), te079, rfac)
        call eval_ops(itri, ti_field(0), ti079, rfac)
        nt79  = n179  + n079
        net79 = ne179 + ne079
        tit79 = ti179 + ti079
        tet79 = te179 + te079
     else
        nt79  = n179
        net79 = ne179
        tit79 = ti179
        tet79 = te179
     end if

     ! ion pressure
     if(linear.eq.1 .and. ilin.eq.1) then 
        temp79a = n079(:,OP_1)*ti179(:,OP_1) + n179(:,OP_1)*ti079(:,OP_1)
     else
        ! calculate the full pressure even when eqsubtract=1,
        ! since we don't store the equilibrium impurity density
        temp79a = nt79(:,OP_1)*tit79(:,OP_1)

        ! Add impurity pressure
        if(ikprad.ne.0) then
           do i=1, kprad_z
              call eval_ops(itri, kprad_n(i), n079, rfac)
              temp79a = temp79a + n079(:,OP_1)*tit79(:,OP_1)
           end do
        end if
     end if

     ! electron pressure
     if(linear.eq.1 .and. ilin.eq.1) then 
        temp79b = ne079(:,OP_1)*te179(:,OP_1) + ne179(:,OP_1)*te079(:,OP_1)
     else
        ! calculate the full pressure even when eqsubtract=1,
        ! since we don't store the equilibrium impurity density
        temp79b = net79(:,OP_1)*tet79(:,OP_1)
     end if


     ! Total pressure is ion pressure plus electron pressure
     temp79a = temp79a + temp79b


     ! If eqsubtract = 1, subtract the equilibrium pressure
     if(ieqsub.eq.1 .and. ilin.eq.1 .and. linear.eq.0) then
        call eval_ops(itri, p_field(0), p079, rfac)
        call eval_ops(itri, pe_field(0), pe079, rfac)

        temp79a = temp79a - p079(:,OP_1)
        temp79b = temp79b - pe079(:,OP_1)
     end if

     dofs_p  = intx2(mu79(:,:,OP_1),temp79a)
     dofs_pe = intx2(mu79(:,:,OP_1),temp79b)

     call vector_insert_block(pe_f%vec,itri,1,dofs_pe,VEC_ADD)
     call vector_insert_block(p_f%vec,itri,1,dofs_p,VEC_ADD)
  end do

  call newvar_solve(pe_f%vec, mass_mat_lhs)
  call newvar_solve(p_f%vec, mass_mat_lhs)

  pe = pe_f
  p = p_f

  call destroy_field(pe_f)
  call destroy_field(p_f)

  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, ' Done calculating pressures'
  
end subroutine calculate_pressures


subroutine calculate_temperatures(ilin, te, ti, pe, p, ne, nion, ieqsub)
  use math
  use basic
  use m3dc1_nint
  use newvar_mod
  use diagnostics
  use metricterms_new
  use field
  use arrays
  use kprad_m3dc1

  implicit none

  type(field_type) :: pe, p, te, ti, ne, nion
  integer, intent(in) :: ilin    ! 0 = calculate only equilibrium part
  integer, intent(in) :: ieqsub  ! 1 = supplied fields are only perturbed parts

  integer :: numelms, itri, i

  type(field_type) :: te_f, ti_f

  vectype, dimension(dofs_per_element) :: dofs_ti, dofs_te

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Calculating Temperatures'

if (ispradapt .eq. 1) then
  call create_field(te_f, "te_f")
  call create_field(ti_f, "ti_f")
else
  call create_field(te_f)
  call create_field(ti_f)
endif
  te_f = 0.
  ti_f = 0.

  numelms = local_elements()
  do itri=1,numelms
     call define_element_quadrature(itri, int_pts_aux, 5)
     call define_fields(itri, 0, 1, 0, ieqsub)

     call eval_ops(itri, ne, ne179, rfac)
     call eval_ops(itri, nion, n179,  rfac)
     call eval_ops(itri, p, p179, rfac)
     call eval_ops(itri, pe, pe179, rfac)
     if(ieqsub.eq.1 .and. ilin.eq.1) then
        call eval_ops(itri, ne_field(0), ne079, rfac)
        call eval_ops(itri, den_field(0), n079, rfac)
        call eval_ops(itri, pe_field(0), pe079, rfac)
        call eval_ops(itri, p_field(0),   p079, rfac)
        nt79  = n179  +  n079
        net79 = ne179 + ne079
        pt79  = p179  +  p079
        pet79 = pe179 + pe079
     else
        nt79  =  n179
        net79 = ne179
        pt79  =  p179
        pet79 = pe179
     end if

     ! Electron temperature
     if(linear.eq.1 .and. ilin.eq.1) then 
        where(real(ne079(:,OP_1)).gt.0.)
           temp79a = pe179(:,OP_1)/ne079(:,OP_1) &
                - ne179(:,OP_1)*pe079(:,OP_1)/ne079(:,OP_1)**2
        elsewhere
           temp79a = 0.
        end where
     else
        ! calculate the full temperature even when eqsubtract=1,
        ! since we don't store the equilibrium impurity density
        where(real(net79(:,OP_1)).gt.0.)
           temp79a = pet79(:,OP_1)/net79(:,OP_1)
        elsewhere
           temp79a = 0.
        end where
     end if

     ! Ion temperature
     ! Ti = (p - pe) / [n_i + sum(n_z)]
     if(linear.eq.1 .and. ilin.eq.1) then
        where(real(n079(:,OP_1)).gt.0.)
           temp79b = (p179(:,OP_1) - pe179(:,OP_1))/n079(:,OP_1) &
                - n179(:,OP_1)*(p079(:,OP_1) - pe079(:,OP_1))/n079(:,OP_1)**2
        elsewhere
           temp79b = 0.
        end where
     else
           
        ! ion density
        temp79c = nt79(:,OP_1)
        ! Add impurity density
        if(ikprad.ne.0) then
           do i=1, kprad_z
              call eval_ops(itri, kprad_n(i), n079, rfac)
              temp79c = temp79c + n079(:,OP_1)
           end do
        end if
        
        ! calculate the full temperature even when eqsubtract=1,
        ! since we don't store the equilibrium impurity density
        where(real(temp79c) .gt. 0.)
           temp79b = (pt79(:,OP_1) - pet79(:,OP_1))/temp79c
        elsewhere
           temp79b = 0.
        end where
     end if

     ! If eqsubtract = 1, subtract the equilibrium pressure
     if(ieqsub.eq.1 .and. ilin.eq.1 .and. linear.eq.0) then
        call eval_ops(itri, te_field(0), te079, rfac)
        temp79a = temp79a - te079(:,OP_1)

        if(ipressplit.eq.1 .and. ipres.eq.1) then
           call eval_ops(itri, ti_field(0), ti079, rfac)
           temp79b = temp79b - ti079(:,OP_1)
        end if
     end if

     dofs_te = intx2(mu79(:,:,OP_1),temp79a)
     dofs_ti = intx2(mu79(:,:,OP_1),temp79b)
        
     call vector_insert_block(te_f%vec,itri,1,dofs_te,VEC_ADD)
     call vector_insert_block(ti_f%vec,itri,1,dofs_ti,VEC_ADD)
  end do

  call newvar_solve(te_f%vec, mass_mat_lhs)
  te = te_f
  call destroy_field(te_f)

  call newvar_solve(ti_f%vec, mass_mat_lhs)
  ti = ti_f
  call destroy_field(ti_f)

  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, ' Done calculating temperatures'
  
end subroutine calculate_temperatures



subroutine calculate_ne(ilin, nion, ne, ieqsub)
  use math
  use basic
  use m3dc1_nint
  use newvar_mod
  use diagnostics
  use metricterms_new
  use field
  use kprad_m3dc1
  use arrays

  implicit none

  type(field_type) :: ne, nion
  integer, intent(in) :: ilin, ieqsub

  integer :: numelms
  integer :: itri, i

  type(field_type) :: ne_f

  vectype, dimension(dofs_per_element) :: dofs

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Calculating electron density'

  if(ikprad.eq.0 .or. linear.eq.1) then
     ne = nion
     call mult(ne, z_ion)
     return
  end if
 
if (ispradapt .eq. 1) then
  call create_field(ne_f, "ne_f")
else
  call create_field(ne_f)
endif
  ne_f = 0.

  numelms = local_elements()
  do itri=1,numelms
     call define_element_quadrature(itri, int_pts_aux, 5)
     call define_fields(itri, 0, 1, 0, ieqsub)

     ! always calculate the full electron density (regardless of eqsubtract)
     ! since we don't store the equilibrium impurity density
     call eval_ops(itri, nion, n179, rfac)
     if(ieqsub.eq.1 .and. ilin.eq.1) then
        call eval_ops(itri, den_field(0), n079, rfac)
        call eval_ops(itri, ne_field(0), ne079, rfac)
        nt79 = n179 + n079
     else
        nt79 = n179
     end if

     temp79a = nt79(:,OP_1)*z_ion

     if(ikprad.ne.0) then
        do i=1, kprad_z
           call eval_ops(itri, kprad_n(i), n079, rfac)
           temp79a = temp79a + i*n079(:,OP_1)
        end do
     end if

     ! If eqsubtract = 1, subtract the equilibrium electron density.
     if(ieqsub.eq.1 .and. ilin.eq.1) then
        temp79a = temp79a - ne079(:,OP_1)
     end if

     dofs = intx2(mu79(:,:,OP_1),temp79a)

     call vector_insert_block(ne_f%vec,itri,1,dofs,VEC_ADD)
  end do

  call newvar_solve(ne_f%vec, mass_mat_lhs)

  ne = ne_f

  call destroy_field(ne_f)

  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, ' Done calculating electron density'
  
end subroutine calculate_ne


subroutine calculate_sigma_e(itri)
  use basic
  use kprad
  use kprad_m3dc1
  use m3dc1_nint

  implicit none

  integer, intent(in) :: itri

  sie79 = sig79

  if(ikprad.ne.0) then 
     call eval_ops(itri, kprad_sigma_e, tm79, rfac)

     sie79 = sie79 + tm79
  end if
  
end subroutine calculate_sigma_e


subroutine calculate_sigma_i(itri)
  use basic
  use kprad
  use kprad_m3dc1
  use m3dc1_nint

  implicit none

  integer, intent(in) :: itri

  sii79 = sig79

  if(ikprad.ne.0) then 
     call eval_ops(itri, kprad_sigma_i, tm79, rfac)

     sii79 = sii79 + tm79
  end if
  
end subroutine calculate_sigma_i

subroutine calculate_sigma_rho(itri)
  use basic
  use kprad
  use kprad_m3dc1
  use m3dc1_nint

  implicit none

  integer, intent(in) :: itri

  sir79 = sig79

  if(ikprad.ne.0) then
     call eval_ops(itri, kprad_sigma_i, tm79, rfac)

     sir79 = sir79 + tm79*kprad_mz/ion_mass
  end if

end subroutine calculate_sigma_rho

subroutine calculate_weighted_density(itri)
  use basic
  use m3dc1_nint
  use kprad_m3dc1

  implicit none

  integer, intent(in) :: itri

  real :: ti_over_te
  integer :: i

  ! Ion density
  nw79 = nt79

  if(ikprad.ne.0) then
     do i=1, kprad_z
        call eval_ops(itri, kprad_n(i), tm79, rfac)
        nw79 = nw79 + tm79
     end do
  end if

  if(ipres.eq.0) then
     ! Total weighted density
     ti_over_te = (1. - pefac) / pefac

     nw79 = ti_over_te*nw79 + net79
  end if
end subroutine calculate_weighted_density

! Total impurity ion density
subroutine calculate_kprad_totden(itri, z)
  use basic
  use kprad
  use kprad_m3dc1
  use m3dc1_nint

  implicit none

  integer, intent(in) :: itri

  integer :: i
  vectype, dimension(MAX_PTS), intent(out) :: z

  z = 0.
  if(ikprad.ne.0) then 
     do i=0, kprad_z
        call eval_ops(itri, kprad_n(i), tm79, rfac)
        z = z + tm79(:,OP_1)
     end do
  end if
end subroutine calculate_kprad_totden


subroutine calculate_auxiliary_fields(ilin)
  use math
  use basic
  use m3dc1_nint
  use newvar_mod
  use diagnostics
  use metricterms_new
  use electric_field
  use temperature_plots
  use kprad_m3dc1
  use arrays
  use diagnostics
  use bootstrap

  implicit none

  integer, intent(in) :: ilin

  integer :: def_fields
  integer :: numelms
  integer :: i, itri, izone
  integer, dimension(MAX_PTS) :: mr

  vectype, dimension(dofs_per_element) :: dofs

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Calculating diagnostic fields'

  bdotgradp = 0.
  bdotgradt = 0.
  torque_density_em = 0.
  torque_density_ntv = 0.
  chord_mask = 0.
  mag_reg = 0.
  ef_r = 0.
  ef_phi = 0.
  ef_z = 0.
  ef_par = 0.
  eta_j = 0.
  mesh_zone = 0.
  z_effective = 0.
  if(ikprad.ne.0) kprad_totden = 0.
  if(jadv.eq.0) then
     psidot = 0.
     veldif = 0.
     eta_jdb = 0.
     bdgp = 0.
     vlbdgp = 0.
  endif
  if(itemp_plot.eq.1) then
     vdotgradt = 0.
     adv1 = 0.
     adv2 = 0.
     adv3 = 0.
     deldotq_perp = 0.
     deldotq_par = 0.
     eta_jsq = 0.
     vpar_field = 0.
     f1vplot = 0.
     f1eplot = 0.
     f2vplot = 0.
     f2eplot = 0.
     f3vplot = 0.
     f3eplot = 0.
     jdbobs  = 0.
     pot2_field = 0.
  endif
  if(ibootstrap.ne.0) then
    Jp_BS_r = 0.
    Jp_BS_z = 0.
    Jp_BS_phi = 0.
    JpdotB = 0.
    JpdotB_dndpsi = 0.
    JpdotB_dtedpsi = 0.
    JpdotB_dtidpsi = 0.
  endif

  ! specify which fields are to be evalulated
  def_fields = FIELD_N + FIELD_NI + FIELD_P + FIELD_PSI + FIELD_I
  def_fields = def_fields + FIELD_PHI + FIELD_V + FIELD_CHI
  def_fields = def_fields + FIELD_ETA + FIELD_TE + FIELD_TI + FIELD_KAP
  def_fields = def_fields + FIELD_MU + FIELD_B2I
  if(jadv.eq.0) def_fields = def_fields + FIELD_ES
  if(heat_source .and. itemp_plot.eq.1) def_fields = def_fields + FIELD_Q
  if(rad_source .and. itemp_plot.eq.1) def_fields = def_fields + FIELD_RAD
  if(ibootstrap.ne.0) def_fields = def_fields + FIELD_JBS
  if (irunaway.ge.1) def_fields = def_fields + FIELD_RE

  numelms = local_elements()

  ! EM Torque density
  do itri=1,numelms
     call define_element_quadrature(itri, int_pts_aux, 5)
     call define_fields(itri, def_fields, 1, 0)

     call get_zone(itri, izone)

     temp79a = izone
     dofs = intx2(mu79(:,:,OP_1),temp79a)
     call vector_insert_block(mesh_zone%vec,itri,1,dofs,VEC_ADD)

     ! magnetic torque_density (ignoring toroidal magnetic pressure gradient)
#ifdef USECOMPLEX
     if(numvar.gt.1) then
        dofs = &
             + intx4(mu79(:,:,OP_1),ps179(:,OP_DR),conjg(bz179(:,OP_DZ)),ri_79)&
             - intx4(mu79(:,:,OP_1),ps179(:,OP_DZ),conjg(bz179(:,OP_DR)),ri_79)&
             - intx3(mu79(:,:,OP_1),bfp179(:,OP_DR),conjg(bz179(:,OP_DR))) &
             - intx3(mu79(:,:,OP_1),bfp179(:,OP_DZ),conjg(bz179(:,OP_DZ))) &
             + intx4(mu79(:,:,OP_1),conjg(ps179(:,OP_DR)),bz179(:,OP_DZ),ri_79)&
             - intx4(mu79(:,:,OP_1),conjg(ps179(:,OP_DZ)),bz179(:,OP_DR),ri_79)&
             - intx3(mu79(:,:,OP_1),conjg(bfp179(:,OP_DR)),bz179(:,OP_DR)) &
             - intx3(mu79(:,:,OP_1),conjg(bfp179(:,OP_DZ)),bz179(:,OP_DZ))
     endif
     dofs = dofs / 2.
#else
     dofs = &
          + intx4(mu79(:,:,OP_1),pst79(:,OP_DR),bzt79(:,OP_DZ),ri_79) &
          - intx4(mu79(:,:,OP_1),pst79(:,OP_DZ),bzt79(:,OP_DR),ri_79)
#ifdef USE3D
     if(numvar.gt.1) then
        dofs = dofs &
             - intx3(mu79(:,:,OP_1),bfpt79(:,OP_DR),bzt79(:,OP_DR)) &
             - intx3(mu79(:,:,OP_1),bfpt79(:,OP_DZ),bzt79(:,OP_DZ))
     endif
#endif
     
#endif
     call vector_insert_block(torque_density_em%vec,itri,1,dofs,VEC_ADD)

     ! NTV torque_density
     if(amupar.ne.0.) then
        do i=1, dofs_per_element
           call PVV2(mu79(i,:,:),temp79f)
!           dofs(i) = int1(temp79f)
           dofs(i) = 0.

           if(ilin.eq.1) then
              call PVS1psipsi(CONJUGATE(ph179),ps179,ps079,temp79c)
              call PVS1psib  (CONJUGATE(ph179),ps179,bz079,temp79d)
              call PVS1bb    (CONJUGATE(ph179),bz179,bz079,temp79e)
              temp79a = temp79b + temp79c + temp79d + temp79e
              dofs(i) = int3(vip79(:,OP_1),temp79a,temp79f)
              call PVS1psipsi(CONJUGATE(ph179),ps079,ps179,temp79c)
              call PVS1psib  (CONJUGATE(ph179),ps079,bz179,temp79d)
              call PVS1bb    (CONJUGATE(ph179),bz079,bz179,temp79e)
              temp79a = temp79b + temp79c + temp79d + temp79e
              dofs(i) = int3(vip79(:,OP_1),temp79a,temp79f)
              call PVS1psipsi(ph079,CONJUGATE(ps179),ps179,temp79c)
              call PVS1psib  (ph079,CONJUGATE(ps179),bz179,temp79d)
              call PVS1bb    (ph079,CONJUGATE(bz179),bz179,temp79e)
              temp79a = temp79b + temp79c + temp79d + temp79e
              dofs(i) = dofs(i) + int3(vip79(:,OP_1),temp79a,temp79f)

              if(numvar.ge.2) then
                 call PVS2psipsi(CONJUGATE(vz179),ps179,ps079,temp79c)
                 call PVS2psib  (CONJUGATE(vz179),ps179,bz079,temp79d)
                 call PVS2bb    (CONJUGATE(vz179),bz179,bz079,temp79e)
                 temp79a = temp79b + temp79c + temp79d + temp79e
                 dofs(i) = dofs(i) + int3(vip79(:,OP_1),temp79a,temp79f)
                 call PVS2psipsi(CONJUGATE(vz179),ps079,ps179,temp79c)
                 call PVS2psib  (CONJUGATE(vz179),ps079,bz179,temp79d)
                 call PVS2bb    (CONJUGATE(vz179),bz079,bz179,temp79e)
                 temp79a = temp79b + temp79c + temp79d + temp79e
                 dofs(i) = dofs(i) + int3(vip79(:,OP_1),temp79a,temp79f)
                 call PVS2psipsi(vz079,CONJUGATE(ps179),ps179,temp79c)
                 call PVS2psib  (vz079,CONJUGATE(ps179),bz179,temp79d)
                 call PVS2bb    (vz079,CONJUGATE(bz179),bz179,temp79e)
                 temp79a = temp79b + temp79c + temp79d + temp79e
                 dofs(i) = dofs(i) + int3(vip79(:,OP_1),temp79a,temp79f)
              endif
              
              if(numvar.ge.3) then
                 call PVS3psipsi(CONJUGATE(ch179),ps179,ps079,temp79c)
                 call PVS3psib  (CONJUGATE(ch179),ps179,bz079,temp79d)
                 call PVS3bb    (CONJUGATE(ch179),bz179,bz079,temp79e)
                 temp79a = temp79b + temp79c + temp79d + temp79e
                 dofs(i) = dofs(i) + int3(vip79(:,OP_1),temp79a,temp79f)
                 call PVS3psipsi(CONJUGATE(ch179),ps079,ps179,temp79c)
                 call PVS3psib  (CONJUGATE(ch179),ps079,bz179,temp79d)
                 call PVS3bb    (CONJUGATE(ch179),bz079,bz179,temp79e)
                 temp79a = temp79b + temp79c + temp79d + temp79e
                 dofs(i) = dofs(i) + int3(vip79(:,OP_1),temp79a,temp79f)
                 call PVS3psipsi(ch079,CONJUGATE(ps179),ps179,temp79c)
                 call PVS3psib  (ch079,CONJUGATE(ps179),bz179,temp79d)
                 call PVS3bb    (ch079,CONJUGATE(bz179),bz179,temp79e)
                 temp79a = temp79b + temp79c + temp79d + temp79e
                 dofs(i) = dofs(i) + int3(vip79(:,OP_1),temp79a,temp79f)
              endif
           else
              call PVS1      (pht79,temp79b)
              call PVS1psipsi(pht79,pst79,pst79,temp79c)
              call PVS1psib  (pht79,pst79,bzt79,temp79d)
              call PVS1bb    (pht79,bzt79,bzt79,temp79e)
              temp79a = temp79b + temp79c + temp79d + temp79e
              dofs(i) = int3(vip79(:,OP_1),temp79a,temp79f)

              if(numvar.ge.2) then
                 call PVS2      (vzt79,temp79b)
                 call PVS2psipsi(vzt79,pst79,pst79,temp79c)
                 call PVS2psib  (vzt79,pst79,bzt79,temp79d)
                 call PVS2bb    (vzt79,bzt79,bzt79,temp79e)
                 temp79a = temp79b + temp79c + temp79d + temp79e
                 dofs(i) = dofs(i) + int3(vip79(:,OP_1),temp79a,temp79f)
              endif
              
              if(numvar.ge.3) then
                 call PVS3      (cht79,temp79b)
                 call PVS3psipsi(cht79,pst79,pst79,temp79c)
                 call PVS3psib  (cht79,pst79,bzt79,temp79d)
                 call PVS3bb    (cht79,bzt79,bzt79,temp79e)
                 temp79a = temp79b + temp79c + temp79d + temp79e
                 dofs(i) = dofs(i) + int3(vip79(:,OP_1),temp79a,temp79f)
              endif
           endif
        end do
     else
        dofs = 0.
     endif
     call vector_insert_block(torque_density_ntv%vec,itri,1,dofs,VEC_ADD)

     ! b dot grad p
     if(ilin.eq.1) then 
        dofs = intx4(mu79(:,:,OP_1),p179(:,OP_DZ),ps079(:,OP_DR),ri_79) &
             - intx4(mu79(:,:,OP_1),p179(:,OP_DR),ps079(:,OP_DZ),ri_79) &
             + intx4(mu79(:,:,OP_1),p079(:,OP_DZ),ps179(:,OP_DR),ri_79) &
             - intx4(mu79(:,:,OP_1),p079(:,OP_DR),ps179(:,OP_DZ),ri_79)
#if defined(USECOMPLEX) || defined(USE3D)
        if(numvar.gt.1) then
           dofs = dofs &
                + intx4(mu79(:,:,OP_1),p179(:,OP_DP),bz079(:,OP_1),ri2_79) &
                - intx3(mu79(:,:,OP_1),p079(:,OP_DZ),bfp179(:,OP_DZ)) &
                - intx3(mu79(:,:,OP_1),p079(:,OP_DR),bfp179(:,OP_DR))
        endif
#endif
     else
        dofs = intx4(mu79(:,:,OP_1),pt79(:,OP_DZ),pst79(:,OP_DR),ri_79) &
             - intx4(mu79(:,:,OP_1),pt79(:,OP_DR),pst79(:,OP_DZ),ri_79)
#if defined(USECOMPLEX) || defined(USE3D)
        if(numvar.gt.1) then
           dofs = dofs &
                + intx4(mu79(:,:,OP_1),pt79(:,OP_DP),bzt79(:,OP_1),ri2_79) &
                - intx3(mu79(:,:,OP_1),pt79(:,OP_DZ),bfpt79(:,OP_DZ)) &
                - intx3(mu79(:,:,OP_1),pt79(:,OP_DR),bfpt79(:,OP_DR))
        endif
#endif
     end if

     call vector_insert_block(bdotgradp%vec,itri,1,dofs,VEC_ADD)


     ! b dot grad T        
     if(ilin.eq.1) then 
        dofs = &
             +intx4(mu79(:,:,OP_1),ri_79,te179(:,OP_DZ),ps079(:,OP_DR)) &
             -intx4(mu79(:,:,OP_1),ri_79,te179(:,OP_DR),ps079(:,OP_DZ)) &
             +intx4(mu79(:,:,OP_1),ri_79,te079(:,OP_DZ),ps179(:,OP_DR)) &
             -intx4(mu79(:,:,OP_1),ri_79,te079(:,OP_DR),ps179(:,OP_DZ))

#if defined(USECOMPLEX) || defined(USE3D)
        dofs = dofs + intx4(mu79(:,:,OP_1),ri2_79,te179(:,OP_DP),bz079(:,OP_1))
        if(numvar.gt.1) then
           dofs = dofs &
                -intx3(mu79(:,:,OP_1),te079(:,OP_DZ),bfp179(:,OP_DZ)) &
                -intx3(mu79(:,:,OP_1),te079(:,OP_DR),bfp179(:,OP_DR))
        endif
#endif
     else
        dofs = &
             +intx4(mu79(:,:,OP_1),ri_79,tet79(:,OP_DZ),pst79(:,OP_DR)) &
             -intx4(mu79(:,:,OP_1),ri_79,tet79(:,OP_DR),pst79(:,OP_DZ))
#if defined(USECOMPLEX) || defined(USE3D)
        dofs = dofs &
             +intx4(mu79(:,:,OP_1),ri2_79,tet79(:,OP_DP),bzt79(:,OP_1))

        if(numvar.gt.1) then
           dofs = dofs &
                -intx3(mu79(:,:,OP_1),tet79(:,OP_DZ),bfpt79(:,OP_DZ)) &
                -intx3(mu79(:,:,OP_1),tet79(:,OP_DR),bfpt79(:,OP_DR))
        endif
#endif
     end if
     call vector_insert_block(bdotgradt%vec,itri,1,dofs,VEC_ADD)


     ! x-ray detector signal
     if(xray_detector_enabled.eq.1) then
        call get_chord_mask(xray_r0, xray_phi0*pi/180., xray_z0, &
             x_79, phi_79, z_79, npoints, &
             xray_theta*pi/180., xray_sigma*pi/180., temp79a)
        dofs = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(chord_mask%vec,itri,1,dofs,VEC_ADD)
     end if

     ! magnetic_region
     call magnetic_region(pst79(:,OP_1),pst79(:,OP_DR),pst79(:,OP_DZ),&
          x_79,z_79,mr)
     temp79a = mr
     dofs = intx2(mu79(:,:,OP_1),temp79a)
     call vector_insert_block(mag_reg%vec,itri,1,dofs,VEC_ADD)

     ! electric_field
     call electric_field_r(ilin,temp79a,izone)
     dofs = intx2(mu79(:,:,OP_1),temp79a)
     call vector_insert_block(ef_r%vec,itri,1,dofs,VEC_ADD)

     call electric_field_phi(ilin,temp79a,izone)
     if(itor.eq.0) temp79a = temp79a/rzero
     dofs = intx2(mu79(:,:,OP_1),temp79a)
     call vector_insert_block(ef_phi%vec,itri,1,dofs,VEC_ADD)

     call electric_field_z(ilin,temp79a,izone)
     dofs = intx2(mu79(:,:,OP_1),temp79a)
     call vector_insert_block(ef_z%vec,itri,1,dofs,VEC_ADD)

     call electric_field_par(ilin,temp79a,izone)
     dofs = intx2(mu79(:,:,OP_1),temp79a)
     call vector_insert_block(ef_par%vec,itri,1,dofs,VEC_ADD)

     call electric_field_eta_j(ilin,temp79a)
     dofs = intx2(mu79(:,:,OP_1),temp79a)
     call vector_insert_block(eta_j%vec,itri,1,dofs,VEC_ADD)

     call calculate_zeff(itri, temp79a)
     dofs = intx2(mu79(:,:,OP_1),temp79a)
     call vector_insert_block(z_effective%vec,itri,1,dofs,VEC_ADD)

     if(ikprad.ne.0) then
        call calculate_kprad_totden(itri, temp79a)
        dofs = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(kprad_totden%vec,itri,1,dofs,VEC_ADD)
     end if
     
     if(jadv.eq.0) then
        call electric_field_psidot(ilin,temp79a)
        dofs = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(psidot%vec,itri,1,dofs,VEC_ADD)

        call electric_field_veldif(ilin,temp79a)
        dofs = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(veldif%vec,itri,1,dofs,VEC_ADD)

        call ef_eta_jdb(ilin,temp79a)
        dofs = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(eta_jdb%vec,itri,1,dofs,VEC_ADD)

        call ef_bdgp(ilin,temp79a)
        dofs = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(bdgp%vec,itri,1,dofs,VEC_ADD)

        call ef_vlbdgp(ilin,temp79a)
        dofs = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(vlbdgp%vec,itri,1,dofs,VEC_ADD)
     endif

     if(itemp_plot.eq.1) then
        call advection(temp79a)
        dofs = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(vdotgradt%vec,itri,1,dofs,VEC_ADD)

        call advection1(temp79a)
        dofs = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(adv1%vec,itri,1,dofs,VEC_ADD)

        call advection2(temp79a)
        dofs = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(adv2%vec,itri,1,dofs,VEC_ADD)

        call advection3(temp79a)
        dofs = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(adv3%vec,itri,1,dofs,VEC_ADD)
  
        call hf_perp(dofs)
        call vector_insert_block(deldotq_perp%vec,itri,1,dofs,VEC_ADD)
     
        call hf_par(dofs)
        call vector_insert_block(deldotq_par%vec,itri,1,dofs,VEC_ADD)

        call ohmic(temp79a)
        dofs = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(eta_jsq%vec,itri,1,dofs,VEC_ADD)

        call vpar_get(temp79a)
        dofs = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(vpar_field%vec,itri,1,dofs,VEC_ADD)

        call f1vplot_sub(dofs)
        call vector_insert_block(f1vplot%vec,itri,1,dofs,VEC_ADD)

        call f1eplot_sub(dofs)
        call vector_insert_block(f1eplot%vec,itri,1,dofs,VEC_ADD)

        call f2vplot_sub(dofs)
        call vector_insert_block(f2vplot%vec,itri,1,dofs,VEC_ADD)
        sig79 = 0
        if(density_source) call eval_ops(itri, sigma_field, sig79)
        call calculate_sigma_e(itri)
        call f2eplot_sub(dofs)
        call vector_insert_block(f2eplot%vec,itri,1,dofs,VEC_ADD)
        call calculate_qdfac(itri,qd79)
        call f3vplot_sub(dofs)
        call vector_insert_block(f3vplot%vec,itri,1,dofs,VEC_ADD)

        call f3eplot_sub(dofs)
        call vector_insert_block(f3eplot%vec,itri,1,dofs,VEC_ADD)

        call jdbobs_sub(temp79a)
        dofs = intx2(mu79(:,:,OP_1),temp79a)
        call vector_insert_block(jdbobs%vec,itri,1,dofs,VEC_ADD)

     call potential2(dofs)
     call vector_insert_block(pot2_field%vec, itri,1,dofs,VEC_ADD)

        
     end if  ! on itemp_plot.eq.1
     
     if(ibootstrap.ne.0) then
      
         !temp79a =  <J.B>/<B^2> = jbscommon * bootsrap_alpha =  tempD * 1/<B^2> * bootsrap_alpha
         !temp79b =  <J.B>                                    =  tempD           * bootsrap_alpha    
         !temp79c =  dndpsi_term  or dlnpdpsi_term                
         !temp79d =  dtedpsi_term                       
         !temp79e =  dtidpsi_term   

      if (ibootstrap.eq.1)then
         !using da/dpsi:   da/dpsi  = 1/|Bp|^2 1/R [(1/R psi_z + f'_r) a_z     + (1/R psi_r - f'_z) a_r]
         call calculate_CommonTerm_Lambda(temp79a,temp79b,temp79c,temp79d,temp79e)
      elseif (ibootstrap.eq.2)then
         !using dte/dpsit: da/dpsit = da/dTe dTe/dpsit = (del a.del Te)/(|del Te|^2 + chi^2) dTe/dpsit              
         call calculate_CommonTerm_Lambda_fordtedpsit(temp79a,temp79b,temp79c,temp79d,temp79e)
         
       endif
      

      
      dofs = -intx4(mu79(:,:,OP_1),ri_79,pst79(:,OP_DZ),temp79a) &
             -intx3(mu79(:,:,OP_1),bfpt79(:,OP_DR),temp79a)
      call vector_insert_block(Jp_BS_r%vec,itri,1,dofs,VEC_ADD)
    
      
      dofs = intx4(mu79(:,:,OP_1),ri_79,bzt79(:,OP_1),temp79a) 
      call vector_insert_block(Jp_BS_phi%vec,itri,1,dofs,VEC_ADD)
      
      
      dofs =  intx4(mu79(:,:,OP_1),ri_79,pst79(:,OP_DR),temp79a) &
             -intx3(mu79(:,:,OP_1),bfpt79(:,OP_DZ),temp79a)
      call vector_insert_block(Jp_BS_z%vec,itri,1,dofs,VEC_ADD)

      
      dofs = intx2(mu79(:,:,OP_1),temp79b) 
      call vector_insert_block(JpdotB%vec,itri,1,dofs,VEC_ADD)

      dofs = intx2(mu79(:,:,OP_1),temp79c) 
      call vector_insert_block(JpdotB_dndpsi%vec,itri,1,dofs,VEC_ADD)

      dofs = intx2(mu79(:,:,OP_1),temp79d) 
      call vector_insert_block(JpdotB_dtedpsi%vec,itri,1,dofs,VEC_ADD)

      dofs = intx2(mu79(:,:,OP_1),temp79e) 
      call vector_insert_block(JpdotB_dtidpsi%vec,itri,1,dofs,VEC_ADD)

     end if
  end do

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' before bdotgradp solve'
  call newvar_solve(bdotgradp%vec, mass_mat_lhs)
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' before bdotgradt solve'
  call newvar_solve(bdotgradt%vec, mass_mat_lhs)
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' before torque_density solve'
  call newvar_solve(torque_density_em%vec, mass_mat_lhs)
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' before torque_density_ntv solve'
  call newvar_solve(torque_density_ntv%vec, mass_mat_lhs)

  if(xray_detector_enabled.eq.1) then
     call newvar_solve(chord_mask%vec, mass_mat_lhs)
  end if

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' before mag_reg solve'
  call newvar_solve(mag_reg%vec, mass_mat_lhs)

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' before ef solve'
  call newvar_solve(ef_r%vec, mass_mat_lhs)
  call newvar_solve(ef_phi%vec, mass_mat_lhs)
  call newvar_solve(ef_z%vec, mass_mat_lhs)
  call newvar_solve(ef_par%vec, mass_mat_lhs)
  call newvar_solve(eta_j%vec, mass_mat_lhs)
  call newvar_solve(mesh_zone%vec, mass_mat_lhs)
  call newvar_solve(z_effective%vec, mass_mat_lhs)
  if(ikprad.ne.0) call newvar_solve(kprad_totden%vec, mass_mat_lhs)
  if(jadv.eq.0) then
     call newvar_solve(psidot%vec, mass_mat_lhs)
     call newvar_solve(veldif%vec, mass_mat_lhs)
     call newvar_solve(eta_jdb%vec, mass_mat_lhs)
     call newvar_solve(bdgp%vec, mass_mat_lhs)
     call newvar_solve(vlbdgp%vec, mass_mat_lhs)
  endif
  if(itemp_plot.eq.1) then
      if(myrank.eq.0 .and. iprint.ge.1) print *,' before vdotgradt solve'
     call newvar_solve(vdotgradt%vec, mass_mat_lhs)
     call newvar_solve(adv1%vec, mass_mat_lhs)
     call newvar_solve(adv2%vec, mass_mat_lhs)
     call newvar_solve(adv3%vec, mass_mat_lhs)
     call newvar_solve(deldotq_perp%vec, mass_mat_lhs)
     call newvar_solve(deldotq_par%vec, mass_mat_lhs)
     call newvar_solve(eta_jsq%vec, mass_mat_lhs)
     call newvar_solve(vpar_field%vec, mass_mat_lhs)
     call newvar_solve(f1vplot%vec, mass_mat_lhs)
     call newvar_solve(f1eplot%vec, mass_mat_lhs)
     call newvar_solve(f2vplot%vec, mass_mat_lhs)
     call newvar_solve(f2eplot%vec, mass_mat_lhs)
     call newvar_solve(f3vplot%vec, mass_mat_lhs)
     call newvar_solve(f3eplot%vec, mass_mat_lhs)
     call newvar_solve(jdbobs%vec , mass_mat_lhs)
  if(myrank.eq.0 .and. iprint.ge.1) print *, 'before pot2 solve'
     call newvar_solve(pot2_field%vec, pot2_mat_lhs)

  endif
  if(myrank.eq.0 .and. iprint.ge.1) print *, 'before J parallel bootstrap solve'
   if(ibootstrap.ne.0) then
    call newvar_solve(Jp_BS_r%vec, mass_mat_lhs)
    call newvar_solve(Jp_BS_z%vec, mass_mat_lhs)
    call newvar_solve(Jp_BS_phi%vec, mass_mat_lhs)
    call newvar_solve(JpdotB%vec, mass_mat_lhs)
    call newvar_solve(JpdotB_dndpsi%vec, mass_mat_lhs)
    call newvar_solve(JpdotB_dtedpsi%vec, mass_mat_lhs)
    call newvar_solve(JpdotB_dtidpsi%vec, mass_mat_lhs)
   endif
  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Done calculating diagnostic fields'
  
  end subroutine calculate_auxiliary_fields

  subroutine calculate_electric_fields(ilin)
  use math
  use basic
  use m3dc1_nint
  use newvar_mod
  use diagnostics
  use metricterms_new
  use electric_field
  use temperature_plots
  use kprad_m3dc1

  implicit none

  integer, intent(in) :: ilin

  integer :: def_fields
  integer :: numelms
  integer :: i, itri, izone
  vectype, dimension(dofs_per_element) :: dofs

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Calculating electric fields'
  
  ef_r = 0.
  ef_phi = 0.
  ef_z = 0.
  !if (ilin==0) b0f = 0.
  !b1f = 0.

  ! specify which fields are to be evalulated
  def_fields = FIELD_N + FIELD_NI + FIELD_P + FIELD_PSI + FIELD_I
  def_fields = def_fields + FIELD_PHI + FIELD_V + FIELD_CHI
  def_fields = def_fields + FIELD_ETA + FIELD_TE + FIELD_KAP
  def_fields = def_fields + FIELD_MU + FIELD_B2I
  if (irunaway.ge.1) def_fields = def_fields + FIELD_RE

  numelms = local_elements()

  do itri=1,numelms
     call define_element_quadrature(itri, int_pts_aux, 5)
     call define_fields(itri, def_fields, 1, 0)

     ! electric_field
     call electric_field_r(ilin,temp79a,izone)
     dofs = intx2(mu79(:,:,OP_1),temp79a)
     call vector_insert_block(ef_r%vec,itri,1,dofs,VEC_ADD)

     call electric_field_phi(ilin,temp79a,izone)
     dofs = intx2(mu79(:,:,OP_1),temp79a)
     call vector_insert_block(ef_phi%vec,itri,1,dofs,VEC_ADD)

     call electric_field_z(ilin,temp79a,izone)
     dofs = intx2(mu79(:,:,OP_1),temp79a)
     call vector_insert_block(ef_z%vec,itri,1,dofs,VEC_ADD)
    
  end do

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' before ef solve'
  call newvar_solve(ef_r%vec, mass_mat_lhs)
  call newvar_solve(ef_phi%vec, mass_mat_lhs)
  call newvar_solve(ef_z%vec, mass_mat_lhs)

  if(myrank.eq.0 .and. iprint.ge.1) print *, ' Done calculating electric fields'
  
  end subroutine calculate_electric_fields

end module auxiliary_fields
