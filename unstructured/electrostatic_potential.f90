module electrostatic_potential

contains

! B4e
! ===
vectype function b4e(e,f)

  use basic
  use m3dc1_nint

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
end function b4e


! B4psieta
! ========
vectype function b4psieta(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

#ifdef USECOMPLEX
  if(surface_int) then 
     temp = &
          - int5(ri4_79,e(:,OP_1),norm79(:,1),f(:,OP_DRP),g(:,OP_1)) &
          - int5(ri4_79,e(:,OP_1),norm79(:,2),f(:,OP_DZP),g(:,OP_1))
  else
     temp = int4(ri4_79,e(:,OP_DZ),f(:,OP_DZP),g(:,OP_1)) &
          + int4(ri4_79,e(:,OP_DR),f(:,OP_DRP),g(:,OP_1))
  endif
#else
  temp = 0.
#endif
  
  b4psieta = temp
end function b4psieta


! B4psietahyp
! ===========
vectype function b4psietahyp(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

#ifdef USECOMPLEX
  if(surface_int) then 
     temp = 0.
  else
     if(ihypeta.eq.0) then
        temp79a = e(:,OP_DZZ)
        if(itor.eq.1) temp79a = temp79a -    ri_79*e(:,OP_DR)
        temp79b = e(:,OP_DRR)
        if(itor.eq.1) temp79b = temp79b - 3.*ri_79*e(:,OP_DR)
        temp79c = e(:,OP_DRZ)
        if(itor.eq.1) temp79c = temp79c -    ri_79*e(:,OP_DZ)             
        
        temp = 2.*int4(ri4_79,temp79a,f(:,OP_DZZP),h(:,OP_1)) &
             + 2.*int4(ri4_79,temp79b,f(:,OP_DRRP),h(:,OP_1)) &
             + 4.*int4(ri4_79,temp79c,f(:,OP_DRZP),h(:,OP_1)) &
             -    int4(ri4_79,e(:,OP_GS),f(:,OP_GSP),h(:,OP_1))
        
        if(itor.eq.1) then
           temp = temp &
                - 4.*int4(ri5_79,temp79b,f(:,OP_DRP),h(:,OP_1)) &
                - 4.*int4(ri5_79,temp79c,f(:,OP_DZP),h(:,OP_1))
        endif
     else
        temp = 0.
     endif
  endif
#else
  temp = 0.
#endif
  
  b4psietahyp = temp
end function b4psietahyp



! B4beta
! ======
vectype function b4beta(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(surface_int) then 
     temp = int5(ri3_79,e(:,OP_1),norm79(:,1),f(:,OP_DZ),g(:,OP_1)) &
          - int5(ri3_79,e(:,OP_1),norm79(:,2),f(:,OP_DR),g(:,OP_1))
  else
     temp = int4(ri3_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_1)) &
          - int4(ri3_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_1))
  endif
  
  b4beta = temp
end function b4beta


! B4betahyp
! =========
vectype function b4betahyp(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(surface_int) then 
     temp = 0.
  else
     if(ihypeta.eq.0) then
        temp79a = e(:,OP_DZZ)
        if(itor.eq.1) temp79a = temp79a -    ri_79*e(:,OP_DR)
        temp79b = e(:,OP_DRR)
        if(itor.eq.1) temp79b = temp79b - 3.*ri_79*e(:,OP_DR)
        temp79c = e(:,OP_DRZ)
        if(itor.eq.1) temp79c = temp79c -    ri_79*e(:,OP_DZ)             
        
        temp = 2.*int4(ri3_79,temp79a,f(:,OP_DRZ),h(:,OP_1)) &
             - 2.*int4(ri3_79,temp79b,f(:,OP_DRZ),h(:,OP_1)) &
             - 2.*int4(ri3_79,temp79c,f(:,OP_DZZ),h(:,OP_1)) &
             + 2.*int4(ri3_79,temp79c,f(:,OP_DRR),h(:,OP_1)) 
        
        if(itor.eq.1) then 
           temp = temp &
                + 2.*int4(ri4_79,temp79b,f(:,OP_DZ),h(:,OP_1)) &
                - 2.*int4(ri4_79,temp79c,f(:,OP_DR),h(:,OP_1))
        endif
        
#ifdef USECOMPLEX
        temp = temp &
             + 4.*int4(ri5_79,e(:,OP_DR),f(:,OP_DZPP),h(:,OP_1)) &
             - 4.*int4(ri5_79,e(:,OP_DZ),f(:,OP_DRPP),h(:,OP_1))
#endif
     else
        temp = 0.
     endif
  endif
  
  b4betahyp = temp
end function b4betahyp


! B4feta
! ======
vectype function b4feta(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

#ifdef USECOMPLEX
  if(surface_int) then 
     temp = int5(ri3_79,e(:,OP_1),norm79(:,1),f(:,OP_DZPP),g(:,OP_1)) &
          - int5(ri3_79,e(:,OP_1),norm79(:,2),f(:,OP_DRPP),g(:,OP_1))
  else
     temp = int4(ri3_79,e(:,OP_DZ),f(:,OP_DRPP),g(:,OP_1)) &
          - int4(ri3_79,e(:,OP_DR),f(:,OP_DZPP),g(:,OP_1))
  endif
#else
  temp = 0.
#endif
  
  b4feta = temp
end function b4feta


! B4fetahyp
! =========
vectype function b4fetahyp(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

#ifdef USECOMPLEX
  if(surface_int) then 
     temp = 0.
  else
     if(ihypeta.eq.0) then
        temp79a = e(:,OP_DZZ)
        if(itor.eq.1) temp79a = temp79a -    ri_79*e(:,OP_DR)
        temp79b = e(:,OP_DRR)
        if(itor.eq.1) temp79b = temp79b - 3.*ri_79*e(:,OP_DR)
        temp79c = e(:,OP_DRZ)
        if(itor.eq.1) temp79c = temp79c -    ri_79*e(:,OP_DZ)             
        
        temp = 2.*int4(ri3_79,temp79a,f(:,OP_DRZPP),h(:,OP_1)) &
             - 2.*int4(ri3_79,temp79b,f(:,OP_DRZPP),h(:,OP_1)) &
             - 2.*int4(ri3_79,temp79c,f(:,OP_DZZPP),h(:,OP_1)) &
             + 2.*int4(ri3_79,temp79c,f(:,OP_DRRPP),h(:,OP_1)) 
        
        if(itor.eq.1) then 
           temp = temp &
                + 2.*int4(ri4_79,temp79b,f(:,OP_DZPP),h(:,OP_1)) &
                - 2.*int4(ri4_79,temp79c,f(:,OP_DRPP),h(:,OP_1))
        endif
     else
        temp = 0.
     endif
  endif
#else
  temp = 0.
#endif
  
  b4fetahyp = temp
end function b4fetahyp



! B4psiv
! ======
vectype function b4psiv(e,f,g)

  use basic
  use m3dc1_nint

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
end function b4psiv

! B4bu
! ====
vectype function b4bu(e,f,g)

  use basic
  use m3dc1_nint

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
end function b4bu

! B4bchi
! ======
vectype function b4bchi(e,f,g)

  use basic
  use m3dc1_nint

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
end function b4bchi

! B4fv
! ====
vectype function b4fv(e,f,g)

  use basic
  use m3dc1_nint

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
end function b4fv


! B4psipsid
! ~~~~~~~~~
vectype function b4psipsid(e,f,g,h)

  use basic
  use m3dc1_nint
  
  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itwofluid.eq.0 .or. db.eq.0.) then
     b4psipsid = 0.
     return
  endif

  if(surface_int) then
     temp79a = ri4_79*e(:,OP_1)*g(:,OP_GS)*h(:,OP_1)
     temp = int3(temp79a,norm79(:,1),f(:,OP_DR)) &
          + int3(temp79a,norm79(:,2),f(:,OP_DZ))
  else
     temp = &
          - int5(ri4_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_GS),h(:,OP_1)) &
          - int5(ri4_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_GS),h(:,OP_1))
  endif

  b4psipsid = temp
end function b4psipsid


! B4psibd
! ~~~~~~~
vectype function b4psibd(e,f,g,h)

  use basic
  use m3dc1_nint
  
  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itwofluid.eq.0 .or. db.eq.0.) then
     b4psibd = 0.
     return
  endif

#ifdef USECOMPLEX
  if(surface_int) then
     temp79a = ri5_79*e(:,OP_1)*g(:,OP_1)*h(:,OP_1)
     temp = int3(temp79a,norm79(:,1),f(:,OP_DZP)) &
          - int3(temp79a,norm79(:,2),f(:,OP_DRP))
  else
     temp = int5(ri5_79,e(:,OP_DZ),f(:,OP_DRP),g(:,OP_1),h(:,OP_1)) &
          - int5(ri5_79,e(:,OP_DR),f(:,OP_DZP),g(:,OP_1),h(:,OP_1))
  endif
#else
  temp = 0.
#endif

  b4psibd = temp
end function b4psibd


! B4psifd
! ~~~~~~~
vectype function b4psifd(e,f,g,h)

  use basic
  use m3dc1_nint
  
  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itwofluid.eq.0 .or. db.eq.0.) then
     b4psifd = 0.
     return
  endif

#ifdef USECOMPLEX
  if(surface_int) then
     temp79a = ri3_79*e(:,OP_1)*f(:,OP_GS)*h(:,OP_1)
     temp = int3(temp79a,norm79(:,2),g(:,OP_DRP)) &
          + int3(temp79a,norm79(:,1),g(:,OP_DZP))
  else
     temp = int5(ri3_79,e(:,OP_DR),g(:,OP_DZP),f(:,OP_GS),h(:,OP_1)) &
          - int5(ri3_79,e(:,OP_DZ),g(:,OP_DRP),f(:,OP_GS),h(:,OP_1))
  endif
#else
  temp = 0.
#endif

  b4psifd = temp
end function b4psifd


! B4bbd
! ~~~~~
vectype function b4bbd(e,f,g,h)

  use basic
  use m3dc1_nint
  
  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itwofluid.eq.0 .or. db.eq.0.) then
     b4bbd = 0.
     return
  endif

  if(surface_int) then
     temp79a = ri4_79*e(:,OP_1)*g(:,OP_1)*h(:,OP_1)
     temp = int3(temp79a,norm79(:,1),f(:,OP_DR)) &
          + int3(temp79a,norm79(:,2),f(:,OP_DZ))
  else
     temp = &
          - int5(ri4_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1),h(:,OP_1)) &
          - int5(ri4_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_1),h(:,OP_1))
  endif

  b4bbd = temp
end function b4bbd


! B4bfd
! ~~~~~
vectype function b4bfd(e,f,g,h)

  use basic
  use m3dc1_nint
  
  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itwofluid.eq.0 .or. db.eq.0.) then
     b4bfd = 0.
     return
  endif

#ifdef USECOMPLEX
  if(surface_int) then
     temp79a = ri4_79*e(:,OP_1)*f(:,OP_1)*h(:,OP_1)
     temp = int3(temp79a,norm79(:,1),g(:,OP_DRPP)) &
          + int3(temp79a,norm79(:,2),g(:,OP_DZPP))
  else
     temp = &
          - int5(ri4_79,e(:,OP_DZ),g(:,OP_DZPP),f(:,OP_1),h(:,OP_1)) &
          - int5(ri4_79,e(:,OP_DR),g(:,OP_DRPP),f(:,OP_1),h(:,OP_1))
  endif
#else
  temp = 0.
#endif

  b4bfd = temp
end function b4bfd


! B4ped
! ~~~~~
vectype function b4ped(e,f,g)

  use basic
  use m3dc1_nint
  
  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(itwofluid.eq.0 .or. db.eq.0.) then
     b4ped = 0.
     return
  endif

  if(surface_int) then
     temp = int5(ri2_79,e(:,OP_1),norm79(:,1),f(:,OP_DR),g(:,OP_1)) &
          + int5(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_DZ),g(:,OP_1))
  else
     temp = &
          - int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1)) &
          - int4(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_1))
  endif

  b4ped = temp
end function b4ped



!======================================================================
! Electrostatic Potential Equation
!======================================================================
subroutine potential_lin(trial, lin, ssterm, ddterm, q_ni, r_bf, q_bf, r_e)

  use basic
  use arrays
  use m3dc1_nint
  use metricterms_new

  implicit none

  vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: trial, lin 
  vectype, dimension(num_fields), intent(out) :: ssterm, ddterm
  vectype, intent(out) :: r_e
  vectype, intent(out) :: q_ni, q_bf, r_bf
  vectype :: temp
  real :: thimp_e, thimpb_e, thimpf_e

  vectype, dimension(MAX_PTS, OP_NUM) :: hf
  hf = hypf*sz79

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
  r_bf = 0.
  q_bf = 0.

  r_e = b4e(trial,lin)

  temp = b4psieta   (trial,lin,eta79) &
       + b4psietahyp(trial,lin,eta79,hf)
  ssterm(psi_g) = ssterm(psi_g)      -thimp_e     *temp
  ddterm(psi_g) = ddterm(psi_g) + (1.-thimp_e*bdf)*temp

  temp = b4psiv(trial,lin,vzt79)
  ssterm(psi_g) = ssterm(psi_g) -     thimp_e     *temp
  ddterm(psi_g) = ddterm(psi_g) + (1.-thimp_e*bdf)*temp

  if(linear.eq.0) then
     temp = b4bu(trial,bz179,lin)
     ssterm(u_g) = ssterm(u_g) - thimpb_e*temp
     ddterm(u_g) = ddterm(u_g) - thimpb_e*temp*bdf

     temp = (b4psipsid(trial,lin,ps179,ni79) &
            +b4psipsid(trial,ps179,lin,ni79))*dbf
     ssterm(psi_g) = ssterm(psi_g) -     thimpf_e     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (.5-thimpf_e*bdf)*dt*temp

     temp = b4psibd(trial,lin,bz179,ni79)*dbf
     ssterm(psi_g) = ssterm(psi_g) -     thimpf_e     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (.5-thimpf_e*bdf)*dt*temp
  endif

  if(eqsubtract.eq.1) then
     temp = b4bu(trial,bz079,lin)
     ssterm(u_g) = ssterm(u_g) -     thimpb_e     *temp
     ddterm(u_g) = ddterm(u_g) + (1.-thimpb_e*bdf)*temp

     temp = (b4psipsid(trial,lin,ps079,ni79) &
            +b4psipsid(trial,ps079,lin,ni79))*dbf
     ssterm(psi_g) = ssterm(psi_g) -     thimpf_e     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (1.-thimpf_e*bdf)*dt*temp

     temp = b4psibd(trial,lin,bz079,ni79)*dbf
     ssterm(psi_g) = ssterm(psi_g) -     thimpf_e     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (1.-thimpf_e*bdf)*dt*temp
  endif


  if(numvar.ge.2) then 
     temp = b4beta   (trial,lin,eta79)    &
          + b4betahyp(trial,lin,eta79,hf)
     ssterm(bz_g) = ssterm(bz_g)      -thimp_e     *temp
     ddterm(bz_g) = ddterm(bz_g) + (1.-thimp_e*bdf)*temp

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

        temp = b4psibd(trial,ps179,lin,ni79)*dbf &
             + b4bbd  (trial,bz179,lin,ni79)*dbf &
             + b4bbd  (trial,lin,bz179,ni79)*dbf
        ssterm(bz_g) = ssterm(bz_g) -     thimpf_e     *dt*temp
        ddterm(bz_g) = ddterm(bz_g) + (.5-thimpf_e*bdf)*dt*temp
     endif

     if(eqsubtract.eq.1) then
        temp = b4psiv(trial,ps079,lin)
        ssterm(vz_g) = ssterm(vz_g) -     thimpb_e     *temp
        ddterm(vz_g) = ddterm(vz_g) + (1.-thimpb_e*bdf)*temp

        temp = b4psibd(trial,ps079,lin,ni79)*dbf &
             + b4bbd  (trial,bz079,lin,ni79)*dbf &
             + b4bbd  (trial,lin,bz079,ni79)*dbf
        ssterm(bz_g) = ssterm(bz_g) -     thimpf_e     *dt*temp
        ddterm(bz_g) = ddterm(bz_g) + (1.-thimpf_e*bdf)*dt*temp
     endif
  endif


  if(numvar.ge.3) then 
     temp = b4ped(trial,lin,ni79)*dbf*pefac
     ssterm(pe_g) = ssterm(pe_g) -     thimpf_e     *dt*temp
     ddterm(pe_g) = ddterm(pe_g) + (1.-thimpf_e*bdf)*dt*temp

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


  if(idens.eq.1 .and. eqsubtract.eq.1) then
     q_ni = q_ni + dt* &
          (b4psipsid(trial,ps079,ps079,lin)*dbf &
          +b4psibd  (trial,ps079,bz079,lin)*dbf &
          +b4bbd    (trial,bz079,bz079,lin)*dbf &
          +b4ped    (trial,pe079,lin)*dbf)
  endif


  if(i3d.eq.1 .and. numvar.ge.2) then 
     temp = b4fv(trial,bft79,lin)
     ssterm(vz_g) = ssterm(vz_g) -     thimpb_e *temp
     ddterm(vz_g) = ddterm(vz_g) + (1.-thimpb_e)*temp*bdf

     temp = b4psifd(trial,lin,bft79,ni79)*dbf
     ssterm(psi_g) = ssterm(psi_g) -     thimpf_e     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (1.-thimpf_e*bdf)*dt*temp

     temp = b4bfd(trial,lin,bft79,ni79)*dbf
     ssterm(bz_g) = ssterm(bz_g) -     thimpf_e     *dt*temp
     ddterm(bz_g) = ddterm(bz_g) + (1.-thimpf_e*bdf)*dt*temp

     if(eqsubtract.eq.1) then
        q_bf = q_bf + dt* &
             (b4fv   (trial,lin,vz079)          &
             +b4psifd(trial,ps079,lin,ni79)*dbf &
             +b4bfd  (trial,bz079,lin,ni79)*dbf)
     endif

     q_bf = q_bf + &
          (b4feta   (trial,lin,eta79) &
          +b4fetahyp(trial,lin,eta79,hf))
  endif

end subroutine potential_lin

end module electrostatic_potential
