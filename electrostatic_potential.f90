module electrostatic_potential

contains

! B4e
! ===
vectype function b4e(e,f)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f
  vectype :: temp

  if(surface_int) then 
     temp = &
          - int4(ri2_79,e(:,OP_1),norm79(:,1),f(:,OP_DR)) &
          - int4(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_DZ))
  else
     temp = &
          - int3(ri2_79,e(:,OP_DZ),f(:,OP_DZ)) &
          - int3(ri2_79,e(:,OP_DR),f(:,OP_DR))
  endif
  
  b4e = temp
  return
end function b4e


! B4psieta
! ========
vectype function b4psieta(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

#ifdef USECOMPLEX
  select case(ivform)
  case(0)
     ! not yet implemented
     temp = 0.
  case(1)
     if(surface_int) then 
        temp = &
             - int5(ri4_79,e(:,OP_1),norm79(:,1),f(:,OP_DRP),g(:,OP_1)) &
             - int5(ri4_79,e(:,OP_1),norm79(:,2),f(:,OP_DZP),g(:,OP_1))
     else
        temp = int4(ri4_79,e(:,OP_DZ),f(:,OP_DZP),g(:,OP_1)) &
             + int4(ri4_79,e(:,OP_DR),f(:,OP_DRP),g(:,OP_1))
     endif
  end select
#else
  temp = 0.
#endif
  
  b4psieta = temp
  return
end function b4psieta


! B4beta
! ======
vectype function b4beta(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     ! not yet implemented
     temp = 0.
  case(1)
     if(surface_int) then 
        temp = int5(ri3_79,e(:,OP_1),norm79(:,1),f(:,OP_DZ),g(:,OP_1)) &
             - int5(ri3_79,e(:,OP_1),norm79(:,2),f(:,OP_DR),g(:,OP_1))
     else
        temp = int4(ri3_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_1)) &
             - int4(ri3_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_1))
     endif
  end select
  
  b4beta = temp
  return
end function b4beta


! B4feta
! ======
vectype function b4feta(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

#ifdef USECOMPLEX
  select case(ivform)
  case(0)
     ! not yet implemented
     temp = 0.
  case(1)
     if(surface_int) then 
        temp = int5(ri3_79,e(:,OP_1),norm79(:,1),f(:,OP_DZPP),g(:,OP_1)) &
             - int5(ri3_79,e(:,OP_1),norm79(:,2),f(:,OP_DRPP),g(:,OP_1))
     else
        temp = int4(ri3_79,e(:,OP_DZ),f(:,OP_DRPP),g(:,OP_1)) &
             - int4(ri3_79,e(:,OP_DR),f(:,OP_DZPP),g(:,OP_1))
     endif
  end select
#else
  temp = 0.
#endif
  
  b4feta = temp
  return
end function b4feta



! B4psiv
! ======
vectype function b4psiv(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     ! not yet implemented
     temp = 0.
  case(1)
     if(surface_int) then 
        temp = int5(ri2_79,e(:,OP_1),norm79(:,1),f(:,OP_DR),g(:,OP_1)) &
             + int5(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_DZ),g(:,OP_1))
     else
        temp = &
             - int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1)) &
             - int4(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_1))
     endif
  end select
  
  b4psiv = temp
  return
end function b4psiv

! B4bu
! ====
vectype function b4bu(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     ! not yet implemented
     temp = 0.
  case(1)
     if(surface_int) then 
        temp = &
             - int5(ri2_79,e(:,OP_1),norm79(:,1),g(:,OP_DR),f(:,OP_1)) &
             - int5(ri2_79,e(:,OP_1),norm79(:,2),g(:,OP_DZ),f(:,OP_1))
     else
        temp = int4(ri2_79,e(:,OP_DZ),g(:,OP_DZ),f(:,OP_1)) &
             + int4(ri2_79,e(:,OP_DR),g(:,OP_DR),f(:,OP_1))
     endif
  end select
  
  b4bu = temp
  return
end function b4bu

! B4bchi
! ======
vectype function b4bchi(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     ! not yet implemented
     temp = 0.
  case(1)
     if(surface_int) then 
        temp = int5(ri5_79,e(:,OP_1),norm79(:,2),g(:,OP_DR),f(:,OP_1)) &
             - int5(ri5_79,e(:,OP_1),norm79(:,1),g(:,OP_DZ),f(:,OP_1))
     else
        temp = int4(ri5_79,e(:,OP_DR),g(:,OP_DZ),f(:,OP_1)) &
             - int4(ri5_79,e(:,OP_DZ),g(:,OP_DR),f(:,OP_1))
     endif
  end select
  
  b4bchi = temp
  return
end function b4bchi

! B4fv
! ====
vectype function b4fv(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

#ifdef USECOMPLEX
  select case(ivform)
  case(0)
     ! not yet implemented
     temp = 0.
  case(1)
     if(surface_int) then 
        temp = int5(ri_79,e(:,OP_1),norm79(:,2),f(:,OP_DRP),g(:,OP_1)) &
             - int5(ri_79,e(:,OP_1),norm79(:,1),f(:,OP_DZP),g(:,OP_1))
     else
        temp = int4(ri_79,e(:,OP_DR),f(:,OP_DZP),g(:,OP_1)) &
             - int4(ri_79,e(:,OP_DZ),f(:,OP_DRP),g(:,OP_1))
     endif
  end select
#else
  temp = 0.
#endif
  
  b4fv = temp
  return
end function b4fv



!======================================================================
! Electrostatic Potential Equation
!======================================================================
subroutine potential_lin(trial, lin, ssterm, ddterm, q_ni, q_bf, r_e)

  use basic
  use arrays
  use nintegrate_mod
  use metricterms_new

  implicit none

  vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: trial, lin 
  vectype, dimension(num_fields), intent(out) :: ssterm, ddterm
  vectype, intent(out) :: r_e
  vectype, intent(out) :: q_ni, q_bf
  vectype :: temp
  real :: thimp_e, thimpb_e, thimpf_e

  vectype, dimension(MAX_PTS, OP_NUM) :: hp
  hp = hypp*sz79

  thimp_e = thimp

  if(imp_mod.eq.0) then
     thimpb_e = thimp_e
  else
     thimpb_e = 1.
  endif
  
  thimpf_e = thimp_e

  ssterm = 0.
  ddterm = 0.
  q_ni = 0.
  q_bf = 0.

  r_e = b4e(trial,lin)

  temp = b4psieta(trial,lin,eta79)
  ssterm(psi_g) = ssterm(psi_g)       - thimp_e     *temp
  ddterm(psi_g) = ddterm(psi_g) + (1. - thimp_e*bdf)*temp

  temp = b4psiv(trial,lin,vzt79)
  ssterm(psi_g) = ssterm(psi_g) -     thimp_e     *temp
  ddterm(psi_g) = ddterm(psi_g) + (1.-thimp_e*bdf)*temp

  if(linear.eq.0) then
     temp = b4bu(trial,bz179,lin)
     ssterm(u_g) = ssterm(u_g) - thimpb_e*temp
     ddterm(u_g) = ddterm(u_g) - thimpb_e*temp*bdf
  endif

  if(eqsubtract.eq.1) then
     temp = b4bu(trial,bz079,lin)
     ssterm(u_g) = ssterm(u_g) -     thimpb_e     *temp
     ddterm(u_g) = ddterm(u_g) + (1.-thimpb_e*bdf)*temp
  endif

  if(numvar.ge.2) then 
     temp = b4beta(trial,lin,eta79)
     ssterm(bz_g) = ssterm(bz_g)       - thimp_e     *temp
     ddterm(bz_g) = ddterm(bz_g) + (1. - thimp_e*bdf)*temp

     temp = b4bu(trial,lin,pht79)
     ssterm(bz_g) = ssterm(bz_g) -     thimpb_e     *temp
     ddterm(bz_g) = ddterm(bz_g) + (1.-thimpb_e*bdf)*temp

     temp = b4bchi(trial,lin,cht79)
     ssterm(bz_g) = ssterm(bz_g) -     thimpb_e     *temp
     ddterm(bz_g) = ddterm(bz_g) + (1.-thimpb_e*bdf)*temp

     if(linear.eq.0) then 
        temp = b4psiv(trial,ps179,lin)
        ssterm(vz_g) = ssterm(vz_g) - thimpb_e*temp
        ddterm(vz_g) = ddterm(vz_g) - thimpb_e*temp*bdf
     endif

     if(eqsubtract.eq.1) then
        temp = b4psiv(trial,ps079,lin)
        ssterm(vz_g) = ssterm(vz_g) -     thimpb_e     *temp
        ddterm(vz_g) = ddterm(vz_g) + (1.-thimpb_e*bdf)*temp
     endif
  endif

  if(numvar.ge.3) then 
     if(linear.eq.0) then
        temp = b4bchi(trial,bz179,lin)
        ssterm(chi_g) = ssterm(chi_g) - thimpb_e*temp
        ddterm(chi_g) = ddterm(chi_g) - thimpb_e*temp*bdf
     endif

     if(eqsubtract.eq.1) then
        temp = b4bchi(trial,bz079,lin)
        ssterm(chi_g) = ssterm(chi_g) -     thimpb_e     *temp
        ddterm(chi_g) = ddterm(chi_g) + (1.-thimpb_e*bdf)*temp
     end if
  endif

  if(i3d.eq.1 .and. numvar.ge.2) then 
     temp = b4fv(trial,bft79,lin)
     ssterm(vz_g) = ssterm(vz_g) - thimpb_e*temp
     ddterm(vz_g) = ddterm(vz_g) - thimpb_e*temp*bdf

     q_bf = q_bf + &
          (b4feta(trial,lin,eta79) &
          +b4fv  (trial,lin,vz079))
  endif

end subroutine potential_lin

end module electrostatic_potential
