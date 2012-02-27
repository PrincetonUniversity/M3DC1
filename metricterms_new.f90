module metricterms_new

implicit none

contains      

!============================================================================
! V1 TERMS
!============================================================================
  

! V1umu 
! =====
vectype function v1umu(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  select case(ivform)
  case(0)

     if(surface_int) then
        temp = 0.
     else
        temp79a = e(:,OP_GS)*g(:,OP_1) &
             + e(:,OP_DR)*g(:,OP_DR) + e(:,OP_DZ)*g(:,OP_DZ)
        if(itor.eq.1) then 
           temp79a = temp79a  &
                + 4.*ri_79*e(:,OP_DR)*g(:,OP_1 ) &
                + 2.*ri_79*e(:,OP_1) *g(:,OP_DR)
        endif
        
        temp79b = e(:,OP_DR)*f(:,OP_DR) + e(:,OP_DZ)*f(:,OP_DZ)
        temp79c = g(:,OP_DR)*f(:,OP_DR) + g(:,OP_DZ)*f(:,OP_DZ)

        temp = &
             - int2(temp79a,f(:,OP_GS)) &
             - int2(temp79b,g(:,OP_GS)) &
             - int2(temp79c,e(:,OP_GS))
        
        if(itor.eq.1) then
           temp = temp &
                - 2.*int3(ri_79,temp79b,g(:,OP_DR)) &
                - 4.*int3(ri_79,temp79c,e(:,OP_DR)) &
                - 2.*int4(ri_79,e(:,OP_1),g(:,OP_LP),f(:,OP_DR))
        endif
        
#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp + &
             (int4(ri2_79,e(:,OP_DZ),f(:,OP_DZPP),g(:,OP_1)) &
             +int4(ri2_79,e(:,OP_DR),f(:,OP_DRPP),g(:,OP_1)))
        if(itor.eq.1) then
           temp = temp + &
                2.*int4(ri3_79,e(:,OP_1),f(:,OP_DRPP),g(:,OP_1))
        endif
#endif
     endif

  case(1)

     temp79b = f(:,OP_DZZ) - f(:,OP_DRR)
     if(itor.eq.1) temp79b = temp79b - ri_79*f(:,OP_DR)

     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp = 0.
        else
           temp = int5(r2_79,g(:,OP_1),norm79(:,2),e(:,OP_DZ),temp79b) &
                - int5(r2_79,g(:,OP_1),norm79(:,1),e(:,OP_DR),temp79b) &
                + 2.* &
                (int5(r2_79,g(:,OP_1),norm79(:,1),e(:,OP_DZ),f(:,OP_DRZ)) &
                +int5(r2_79,g(:,OP_1),norm79(:,2),e(:,OP_DR),f(:,OP_DRZ)))
           
           if(itor.eq.1) then
              temp79a = h(:,OP_1) - g(:,OP_1)
              temp = temp &
                   + 2.*int5(r_79,g(:,OP_1),norm79(:,1),e(:,OP_DZ),f(:,OP_DZ))&
                   + 4.* &
                   (int5(r_79,norm79(:,1),e(:,OP_DZ),f(:,OP_DZ),temp79a) &
                   -int5(r_79,norm79(:,2),e(:,OP_DR),f(:,OP_DZ),temp79a)) &
                   + 4.* &
                   (int5(r_79,e(:,OP_1),norm79(:,1),f(:,OP_DZZ),h(:,OP_1)) &
                   -int5(r_79,e(:,OP_1),norm79(:,2),f(:,OP_DRZ),h(:,OP_1)))
           endif

#if defined(USE3D) || defined(USECOMPLEX)
           temp = temp &
                - int4(e(:,OP_1),norm79(:,1),f(:,OP_DRPP),g(:,OP_1)) &
                - int4(e(:,OP_1),norm79(:,2),f(:,OP_DZPP),g(:,OP_1))
#endif
        endif
     else 
        temp79a = e(:,OP_DZZ) - e(:,OP_DRR)
        temp79c = 2.*e(:,OP_DRZ)
        temp79d = 2.*f(:,OP_DRZ)
        if(itor.eq.1) then
           temp79a = temp79a - ri_79*e(:,OP_DR)
           temp79c = temp79c + ri_79*e(:,OP_DZ)
           temp79d = temp79d + ri_79*f(:,OP_DZ)
        endif

        temp = &
             - int4(r2_79,temp79a,temp79b,g(:,OP_1)) &
             - int4(r2_79,temp79c,temp79d,g(:,OP_1))
        
        if(itor.eq.1) then
           temp = temp &
                + 5.*int3(e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1)) &
                - 8.*int3(e(:,OP_DZ),f(:,OP_DZ),h(:,OP_1))
        endif
     
#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp &
             + int3(e(:,OP_DZ),f(:,OP_DZPP),g(:,OP_1)) &
             + int3(e(:,OP_DR),f(:,OP_DRPP),g(:,OP_1))
#endif
     end if

  end select

  v1umu = temp
  return
end function v1umu


! V1vmu
! =====
vectype function v1vmu(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  temp = 0.

#if defined(USE3D) || defined(USECOMPLEX)
  select case(ivform)
  case(0)
     ! not implemented
     temp = 0.

  case(1)
     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp = 0.
        else
           temp = int5(r_79,e(:,OP_1),norm79(:,2),f(:,OP_DRP),g(:,OP_1)) &
                - int5(r_79,e(:,OP_1),norm79(:,1),f(:,OP_DZP),g(:,OP_1))
           if(itor.eq.1) then
              temp = temp &
                   + 2.*int4(e(:,OP_1),norm79(:,2),f(:,OP_DP),g(:,OP_1)) &
                   - 4.*int4(e(:,OP_1),norm79(:,2),f(:,OP_DP),h(:,OP_1))
           endif
        end if
     else 
        temp = int4(r_79,e(:,OP_DR),f(:,OP_DZP),g(:,OP_1)) &
             - int4(r_79,e(:,OP_DZ),f(:,OP_DRP),g(:,OP_1))
        
        if(itor.eq.1) then
           temp = temp &
                + 4.*int3(e(:,OP_DZ),f(:,OP_DP),h(:,OP_1)) &
                - 2.*int3(e(:,OP_DZ),f(:,OP_DP),g(:,OP_1))
        endif
     endif
  end select
#endif
  v1vmu = temp
  return
end function v1vmu


! V1chimu
! =======
vectype function v1chimu(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp79a = g(:,OP_DZ)*f(:,OP_DR) - g(:,OP_DR)*f(:,OP_DZ)
        temp79b = r_79*(e(:,OP_DZ)*g(:,OP_DR) - e(:,OP_DR)*g(:,OP_DZ))
        if(itor.eq.1) temp79b = temp79b - 2.*e(:,OP_1)*g(:,OP_DZ)
        
        temp = int3(r_79,e(:,OP_LP),temp79a) &
             + int4(r_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_GS)) &
             - int4(r_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_GS)) &
             + int2(temp79b,f(:,OP_GS))
        
        if(itor.eq.1) then
           temp = temp &
                + 4.*int2(      e(:,OP_DR),temp79a) &
                + 4.*int3(ri_79,e(:,OP_1 ),temp79a) &
                + 4.*int3(ri_79,f(:,OP_DR),temp79b) &
                - 2.*int3(e(:,OP_1),f(:,OP_DZ),g(:,OP_GS))
        endif
     end if

  case(1)
     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp = 0.
        else
           temp = -2.* &
                (int5(ri_79,e(:,OP_DZ),norm79(:,1),f(:,OP_DRR),g(:,OP_1)) &
                +int5(ri_79,e(:,OP_DZ),norm79(:,2),f(:,OP_DRZ),g(:,OP_1)) &
                -int5(ri_79,e(:,OP_DR),norm79(:,1),f(:,OP_DRZ),g(:,OP_1)) &
                -int5(ri_79,e(:,OP_DR),norm79(:,2),f(:,OP_DZZ),g(:,OP_1))) &
                + 2.* &
                (int5(ri_79,norm79(:,1),e(:,OP_DZ),f(:,OP_GS),g(:,OP_1)) &
                -int5(ri_79,norm79(:,2),e(:,OP_DR),f(:,OP_GS),g(:,OP_1)))
           
           if(itor.eq.1) then
              temp = temp + 2.* &
                   (int5(ri2_79,e(:,OP_DZ),norm79(:,1),f(:,OP_DR),g(:,OP_1)) &
                   +int5(ri2_79,e(:,OP_DZ),norm79(:,2),f(:,OP_DZ),g(:,OP_1)) &
                   +int5(ri2_79,e(:,OP_1),norm79(:,1),f(:,OP_DRZ),g(:,OP_1)) &
                   +int5(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_DZZ),g(:,OP_1)) &
                   +int5(ri2_79,e(:,OP_DZ),norm79(:,1),f(:,OP_DR),g(:,OP_1)) &
                   -int5(ri2_79,e(:,OP_DR),norm79(:,1),f(:,OP_DZ),g(:,OP_1)) &
                   -2.*int5(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_GS),h(:,OP_1)) &
                   -2.*int5(ri3_79,e(:,OP_1),norm79(:,1),f(:,OP_DZ),g(:,OP_1)))
           end if
           
#if defined(USE3D) || defined(USECOMPLEX)
           temp = temp &
                + int5(ri3_79,e(:,OP_1),g(:,OP_1),norm79(:,2),f(:,OP_DRPP)) &
                - int5(ri3_79,e(:,OP_1),g(:,OP_1),norm79(:,1),f(:,OP_DZPP))
#endif
        endif
     else 
        temp79a = e(:,OP_DRZ)
        temp79b = f(:,OP_DZZ) - f(:,OP_DRR)
        temp79c = f(:,OP_DRZ)
        temp79d = e(:,OP_DZZ) - e(:,OP_DRR)
        if(itor.eq.1) then
           temp79a = temp79a + .5*ri_79*e(:,OP_DZ)
           temp79b = temp79b + 2.*ri_79*f(:,OP_DR)
           temp79c = temp79c -    ri_79*f(:,OP_DZ)
           temp79d = temp79d -    ri_79*e(:,OP_DR)
        endif
     
        temp = -2.* &
             (int4(ri_79,temp79a,temp79b,g(:,OP_1)) &
             -int4(ri_79,temp79c,temp79d,g(:,OP_1)))
        
        if(itor.eq.1) then
           temp = temp &
                + 4.*int4(ri2_79,e(:,OP_DZ),f(:,OP_GS),h(:,OP_1)) &
                - 3.*int4(ri2_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_1))
        endif
     
#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp &
             + int4(ri3_79,e(:,OP_DR),f(:,OP_DZPP),g(:,OP_1)) &
             - int4(ri3_79,e(:,OP_DZ),f(:,OP_DRPP),g(:,OP_1))
#endif
     endif
  end select

  v1chimu = temp
  return
end function v1chimu




! V1un
! ====
vectype function v1un(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int3(e(:,OP_DR),f(:,OP_DR),g(:,OP_1)) &
             + int3(e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1))
        if(itor.eq.1) then
           temp = temp + 2.*int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_1))
        endif
     endif

  case(1)
     if(surface_int) then 
        if(inoslip_pol.eq.1) then
           temp = 0.
        else
           temp = &
                - int5(r2_79,e(:,OP_1),norm79(:,1),f(:,OP_DR),g(:,OP_1)) &
                - int5(r2_79,e(:,OP_1),norm79(:,2),f(:,OP_DZ),g(:,OP_1))
        endif
     else
        temp = int4(r2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_1)) &
             + int4(r2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1))
     end if
  end select

  v1un = temp
  return
end function v1un


! V1chin
! ======
vectype function v1chin(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int4(r_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ)) &
             - int4(r_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR))
     endif
  case(1)
     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp = 0.
        else
           temp = int5(ri_79,e(:,OP_1),norm79(:,2),f(:,OP_DR),g(:,OP_1)) &
                - int5(ri_79,e(:,OP_1),norm79(:,1),f(:,OP_DZ),g(:,OP_1))
        end if
     else        
        temp = int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_1)) &
             - int4(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_1))
     endif
  end select

  v1chin = temp
  return
end function v1chin


! V1psipsi
! ========
vectype function v1psipsi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int4(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_GS)) &
             - int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_GS))
        
        if(itor.eq.1) then
           temp = temp - 2.*int4(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_GS))
        endif
     end if
  case(1)
     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp = 0.
        else
           temp = int5(ri_79,e(:,OP_1),norm79(:,1),f(:,OP_DZ),g(:,OP_GS)) &
                - int5(ri_79,e(:,OP_1),norm79(:,2),f(:,OP_DR),g(:,OP_GS))
        endif
     else
        temp = int4(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_GS)) &
             - int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_GS))
     endif
  end select

  v1psipsi = temp
  return
end function v1psipsi


! V1psib
! ======
vectype function v1psib(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_DZP),g(:,OP_1)) &
             + int4(ri2_79,e(:,OP_DR),f(:,OP_DRP),g(:,OP_1))

        if(itor.eq.1) then
           temp = temp &
                + 2.*int4(ri3_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_1))
        endif
     end if

  case(1)
     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp = 0.
        else
           temp = &
                - int5(ri2_79,e(:,OP_1),norm79(:,1),f(:,OP_DRP),g(:,OP_1)) &
                - int5(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_DZP),g(:,OP_1))
        end if
     else
        temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_DZP),g(:,OP_1)) &
             + int4(ri2_79,e(:,OP_DR),f(:,OP_DRP),g(:,OP_1))
     endif
  end select

  v1psib = temp
#else
  v1psib = 0.
#endif
end function v1psib


! V1bb
! ====
vectype function v1bb(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        if(itor.eq.0) then
           temp = 0.
        else
           temp = -2.*int4(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_1)) 
        endif
     endif

  case(1)
     if(surface_int) then
        temp = -0.25* &
             (int5(ri_79,norm79(:,1),e(:,OP_DZ),f(:,OP_1),g(:,OP_1)) &
             -int5(ri_79,norm79(:,2),e(:,OP_DR),f(:,OP_1),g(:,OP_1))) &
             - 0.5 * &
             (int5(ri_79,e(:,OP_1),norm79(:,1),f(:,OP_DZ),g(:,OP_1)) &
             -int5(ri_79,e(:,OP_1),norm79(:,2),f(:,OP_DR),g(:,OP_1)))
     else
        temp = 0.
     endif
  end select

  v1bb = temp
  return
end function v1bb


! V1uun 
! =====
vectype function v1uun(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(inertia.eq.0) then
     v1uun = 0.
     return
  end if

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
  
        temp = int5(ri_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DZ),h(:,OP_1)) &
             - int5(ri_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DR),h(:,OP_1)) &
             + 0.5*(int4(ri_79,temp79a,e(:,OP_DR),h(:,OP_DZ)) &
                   -int4(ri_79,temp79a,e(:,OP_DZ),h(:,OP_DR)))

        if(itor.eq.1) then
           temp = temp + &
                2.*int5(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_DZ),h(:,OP_1)) &
                +  int4(ri2_79,temp79a,e(:,OP_1),h(:,OP_DZ))
        endif
     end if

  case(1)
     if(surface_int) then
        temp = 0.
     else
        if(inoslip_pol.eq.1) then
           temp = 0.
        else
           temp79a = &
                e(:,OP_DZ)*(f(:,OP_DZZ)*g(:,OP_DR) - f(:,OP_DRZ)*g(:,OP_DZ)) &
                + e(:,OP_DR)*(f(:,OP_DRZ)*g(:,OP_DR) - f(:,OP_DRR)*g(:,OP_DZ))
           temp = -int3(r3_79,temp79a,h(:,OP_1))
        endif
     
        if(itor.eq.1) then
           temp = temp &
                + int5(r2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
                + int5(r2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1))
        end if
     end if
  end select

  v1uun = temp
  return
end function v1uun


! V1uvn
! =====
vectype function v1uvn(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0) then
     v1uvn = 0.
     return
  end if
  select case(ivform)
  case(0)
     ! not implemented
     if(surface_int) then
        temp = 0.
     else
        temp = 0.
     end if
  case(1)
     if(surface_int) then
        temp = 0.
     else
        temp = &
             - int5(r2_79,e(:,OP_DZ),f(:,OP_DZP),g(:,OP_1),h(:,OP_1)) &
             - int5(r2_79,e(:,OP_DR),f(:,OP_DRP),g(:,OP_1),h(:,OP_1))
     end if
  end select
#else
  temp = 0.
#endif
  v1uvn = temp
  return
end function v1uvn


! v1uchin
! =======
vectype function v1uchin(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(inertia.eq.0) then
     v1uchin = 0.
     return
  end if


  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)

        temp =-int4(e(:,OP_DZ),f(:,OP_GS),g(:,OP_DZ),h(:,OP_1)) &
             - int4(e(:,OP_DR),f(:,OP_GS),g(:,OP_DR),h(:,OP_1)) &
             + int3(temp79a,e(:,OP_DZ),h(:,OP_DR)) &
             - int3(temp79a,e(:,OP_DR),h(:,OP_DZ))

        if(itor.eq.1) then
           temp = temp - 2.*&
                (int5(ri_79,e(:,OP_1),f(:,OP_GS),g(:,OP_DR),h(:,OP_1)) &
                +int4(ri_79,temp79a,e(:,OP_1),h(:,OP_DZ)))
        endif
     end if

  case(1)
     if(surface_int) then
        temp = 0.
     else
        temp79a = e(:,OP_DR)*(f(:,OP_DR)*g(:,OP_DZZ) + f(:,OP_DRZ)*g(:,OP_DZ) &
                             -f(:,OP_DZ)*g(:,OP_DRZ) + f(:,OP_DRR)*g(:,OP_DR))&
                + e(:,OP_DZ)*(f(:,OP_DZ)*g(:,OP_DRR) + f(:,OP_DRZ)*g(:,OP_DR) &
                             -f(:,OP_DR)*g(:,OP_DRZ) + f(:,OP_DZZ)*g(:,OP_DZ))
        temp = -int2(temp79a,h(:,OP_1))

        if(itor.eq.1) then
           temp79a = &
                e(:,OP_DR)*(f(:,OP_DR)*g(:,OP_DR) + f(:,OP_DZ)*g(:,OP_DZ)) &
                + f(:,OP_DZ)*(e(:,OP_DR)*g(:,OP_DZ) - e(:,OP_DZ)*g(:,OP_DR))
           temp = temp - int3(ri_79,temp79a,h(:,OP_1))
        end if
     end if
  end select

  v1uchin = temp
end function v1uchin


! V1vvn
! =====
vectype function v1vvn(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(inertia.eq.0) then
     v1vvn = 0.
     return
  end if

  temp = 0.
  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        if(itor.eq.1) then
           temp = -int5(ri2_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_1),h(:,OP_1))
        endif
     end if
     
  case(1)
     if(surface_int) then
        temp = 0.
     else
        if(itor.eq.1) then
           temp = -int5(r2_79, e(:,OP_DZ),f(:,OP_1),g(:,OP_1),h(:,OP_1))
        endif
     end if
  end select


  v1vvn = temp
  return
end function v1vvn


! V1vchin
! =======
vectype function v1vchin(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0) then
     v1vchin = 0.
     return
  end if
  select case(ivform)
  case(0)
     ! not implemented
     if(surface_int) then
        temp = 0.
     else
        temp = 0.
     endif
  case(1)
     if(surface_int) then
        temp = 0.
     else
        temp = int5(ri_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DRP),h(:,OP_1)) &
             - int5(ri_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DZP),h(:,OP_1))
     endif
  end select
#else
  temp = 0.
#endif
  v1vchin = temp
  return
end function v1vchin


! v1chichin
! =========
vectype function v1chichin(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(inertia.eq.0) then
     v1chichin = 0.
     return
  end if


  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
        
        temp = 0.5* &
             (int4(r_79,e(:,OP_DR),temp79a,h(:,OP_DZ)) &
             -int4(r_79,e(:,OP_DZ),temp79a,h(:,OP_DR)))
        
        if(itor.eq.1) then
           temp = temp + int3(e(:,OP_1),temp79a,h(:,OP_DZ))
        endif
     end if

  case(1)
     if(surface_int) then
        temp = 0.
     else
        temp79a = g(:,OP_DZ)*(e(:,OP_DR)*f(:,OP_DZZ)-e(:,OP_DZ)*f(:,OP_DRZ)) &
             +    g(:,OP_DR)*(e(:,OP_DR)*f(:,OP_DRZ)-e(:,OP_DZ)*f(:,OP_DRR))

        temp = -int3(ri3_79,temp79a,h(:,OP_1))

        if(itor.eq.1) then
           temp79a = g(:,OP_DR)*(e(:,OP_DZ)*f(:,OP_DR) - e(:,OP_DR)*f(:,OP_DZ))
           temp = temp - 2.*int3(ri4_79,temp79a,h(:,OP_1))
        end if
     endif
  end select
  
  v1chichin = temp

  return
end function v1chichin



! V1upsipsi
! =========
vectype function v1upsipsi(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  ! |u, psi(1)|,r
  temp79b = f(:,OP_DRZ)*g(:,OP_DR ) - f(:,OP_DRR)*g(:,OP_DZ ) &
       +    f(:,OP_DZ )*g(:,OP_DRR) - f(:,OP_DR )*g(:,OP_DRZ)
  ! |u, psi(1)|,z
  temp79c = f(:,OP_DZZ)*g(:,OP_DR ) - f(:,OP_DRZ)*g(:,OP_DZ ) &
       +    f(:,OP_DZ )*g(:,OP_DRZ) - f(:,OP_DR )*g(:,OP_DZZ)

  ! |nu, psi(2)|,r
  temp79e = e(:,OP_DRZ)*h(:,OP_DR ) - e(:,OP_DRR)*h(:,OP_DZ ) &
       +    e(:,OP_DZ )*h(:,OP_DRR) - e(:,OP_DR )*h(:,OP_DRZ)
  ! |nu, psi(2)|,z
  temp79f = e(:,OP_DZZ)*h(:,OP_DR ) - e(:,OP_DRZ)*h(:,OP_DZ ) &
       +    e(:,OP_DZ )*h(:,OP_DRZ) - e(:,OP_DR )*h(:,OP_DZZ)

  select case(ivform)
  case(0)  

     if(surface_int) then
        temp = 0.
     else
        
        temp =-int3(ri2_79,temp79b,temp79e) &
             - int3(ri2_79,temp79c,temp79f) &
             + int4(ri2_79,e(:,OP_DZ),temp79b,h(:,OP_GS)) &
             - int4(ri2_79,e(:,OP_DR),temp79c,h(:,OP_GS))

        if(itor.eq.1) then
           ! |u, psi(1)|
           temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
           
           ! |nu,psi(2)|
           temp79d = e(:,OP_DZ)*h(:,OP_DR) - e(:,OP_DR)*h(:,OP_DZ)
           
           temp = temp &
                -    int4(ri3_79,e(:,OP_DZ),temp79a,h(:,OP_GS)) &
                - 2.*int4(ri3_79,e(:,OP_1 ),temp79c,h(:,OP_GS)) &
                +    int3(ri3_79,temp79e,temp79a) &
                -    int3(ri3_79,temp79d,temp79b) &
                + 2.*int4(ri3_79,e(:,OP_1 ),temp79b,h(:,OP_DRZ)) &
                + 2.*int4(ri3_79,e(:,OP_1 ),temp79c,h(:,OP_DZZ)) &
                + 2.*int4(ri3_79,e(:,OP_DR),temp79b,h(:,OP_DZ )) &
                + 2.*int4(ri3_79,e(:,OP_DZ),temp79c,h(:,OP_DZ )) &
                +    int3(ri4_79,temp79a,temp79d) &
                - 2.*int4(ri4_79,e(:,OP_DR),temp79a,h(:,OP_DZ )) &
                - 2.*int4(ri4_79,e(:,OP_1 ),temp79a,h(:,OP_DRZ))
        endif
     end if

  case(1)
     if(surface_int) then
        temp79a = e(:,OP_DZ)*h(:,OP_DR) - e(:,OP_DR)*h(:,OP_DZ)
        temp = int3(temp79a,norm79(:,1),temp79b) &
             + int3(temp79a,norm79(:,2),temp79c)
        if(itor.eq.1) then
           temp79d = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
           temp = temp + &
                int4(ri_79,temp79a,norm79(:,1),temp79d)
        endif
     else
        temp = -int2(temp79b,temp79e)  &
             - int2(temp79c,temp79f)  &
             + int3(e(:,OP_DZ),temp79b,h(:,OP_GS))  &
             - int3(e(:,OP_DR),temp79c,h(:,OP_GS))

        if(itor.eq.1) then
           ! |u, psi(1)|
           temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
           
           ! |nu,psi(2)|
           temp79d = e(:,OP_DZ)*h(:,OP_DR) - e(:,OP_DR)*h(:,OP_DZ)
           temp = temp             &
                - int3(ri_79,temp79d,temp79b)   &
                - int3(ri_79,temp79a,temp79e)   &
                - int3(ri2_79,temp79a,temp79d)   &
                + int4(ri_79,h(:,OP_GS),temp79a,e(:,OP_DZ))
        endif
     endif
  end select


  v1upsipsi = temp
  return
end function v1upsipsi

! V1upsib
! =======
vectype function v1upsib(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp
#if defined(USE3D) || defined(USECOMPLEX)

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp79a = g(:,OP_DR)*(e(:,OP_DZZ)*f(:,OP_DZP)-e(:,OP_DZ)*f(:,OP_DZZP) &
                             +e(:,OP_DRZ)*f(:,OP_DRP)-e(:,OP_DR)*f(:,OP_DRZP))&
                 -g(:,OP_DZ)*(e(:,OP_DRR)*f(:,OP_DRP)-e(:,OP_DR)*f(:,OP_DRRP) &
                             +e(:,OP_DRZ)*f(:,OP_DZP)-e(:,OP_DZ)*f(:,OP_DRZP))


        temp = -int3(ri3_79,temp79a,h(:,OP_1))

        if(itor.eq.1) then
           temp79a = g(:,OP_DZ)* &
                (e(:,OP_1)*f(:,OP_DRRP) - e(:,OP_DZ)*f(:,OP_DZP) &
                - 2.*e(:,OP_DR)*f(:,OP_DRP)) &
                +    g(:,OP_DR)* &
                (e(:,OP_DZ)*f(:,OP_DRP) - e(:,OP_1 )*f(:,OP_DRZP))
           temp = temp &
                - 2.*int3(ri4_79,temp79a,h(:,OP_1)) &
                + 2.*int5(ri5_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_DZ),h(:,OP_1))
        endif
     end if

  case(1)

     if(surface_int) then
        temp79a = f(:,OP_DZP)*g(:,OP_DR ) - f(:,OP_DRP)*g(:,OP_DZ ) &
             +    f(:,OP_DZ )*g(:,OP_DRP) - f(:,OP_DR )*g(:,OP_DZP)
        temp79b = e(:,OP_1)*h(:,OP_DP)
        temp = int5(ri_79,norm79(:,1),e(:,OP_DR),temp79a,h(:,OP_1)) &
             + int5(ri_79,norm79(:,2),e(:,OP_DZ),temp79a,h(:,OP_1)) &
             + int5(ri_79,temp79b,norm79(:,1),f(:,OP_DRZ),g(:,OP_DR )) &
             - int5(ri_79,temp79b,norm79(:,1),f(:,OP_DRR),g(:,OP_DZ )) &
             + int5(ri_79,temp79b,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRR)) &
             - int5(ri_79,temp79b,norm79(:,1),f(:,OP_DR ),g(:,OP_DRZ)) &
             + int5(ri_79,temp79b,norm79(:,2),f(:,OP_DZZ),g(:,OP_DR )) &
             - int5(ri_79,temp79b,norm79(:,2),f(:,OP_DRZ),g(:,OP_DZ )) &
             + int5(ri_79,temp79b,norm79(:,2),f(:,OP_DZ ),g(:,OP_DRZ)) &
             - int5(ri_79,temp79b,norm79(:,2),f(:,OP_DR ),g(:,OP_DZZ))
        if(itor.eq.1) then
           temp79c = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
           temp = temp &
                + int4(ri2_79,temp79b,temp79c,norm79(:,1))
        endif
     else
        temp79a = h(:,OP_1)*e(:,OP_GS)  &
             + h(:,OP_DZ)*e(:,OP_DZ) + h(:,OP_DR)*e(:,OP_DR)
        temp79b = f(:,OP_DZP)*g(:,OP_DR ) - f(:,OP_DRP)*g(:,OP_DZ ) &
             +    f(:,OP_DZ )*g(:,OP_DRP) - f(:,OP_DR )*g(:,OP_DZP)
        temp79c = e(:,OP_DZ)*g(:,OP_DZP) + e(:,OP_DR)*g(:,OP_DRP)
        temp79d = h(:,OP_DP)*(e(:,OP_DZ)*f(:,OP_DR )-e(:,OP_DR)*f(:,OP_DZ )) &
             +    h(:,OP_1 )*(e(:,OP_DZ)*f(:,OP_DRP)-e(:,OP_DR)*f(:,OP_DZP))
        temp79e = h(:,OP_DP)*f(:,OP_GS) + h(:,OP_1)*f(:,OP_GSP) &
             + h(:,OP_DZP)*f(:,OP_DZ ) + h(:,OP_DRP)*f(:,OP_DR ) &
             + h(:,OP_DZ )*f(:,OP_DZP) + h(:,OP_DR )*f(:,OP_DRP)
     
        temp = int3(ri_79,temp79a,temp79b)               &
             + int4(ri_79,f(:,OP_DR),h(:,OP_DZ),temp79c) &
             - int4(ri_79,f(:,OP_DZ),h(:,OP_DR),temp79c) &
             - int3(ri_79,g(:,OP_GS),temp79d)            &
             - int4(ri_79,e(:,OP_DZ),g(:,OP_DR),temp79e) &
             + int4(ri_79,e(:,OP_DR),g(:,OP_DZ),temp79e)
        temp = -temp
     endif
  end select

  v1upsib = temp
#else
  v1upsib = 0.
#endif

end function v1upsib


! V1ubb 
! =====
vectype function v1ubb(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  temp = 0.
  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        if(itor.eq.1) then
           temp = 2.* &
                (   int5(ri3_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
                -   int5(ri3_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1)) &
                -2.*int5(ri4_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1 ),h(:,OP_1)))
        endif
        
#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp + &
             (int5(ri4_79,e(:,OP_DZ),f(:,OP_DZPP),g(:,OP_1),h(:,OP_1)) &
             +int5(ri4_79,e(:,OP_DR),f(:,OP_DRPP),g(:,OP_1),h(:,OP_1)))
        if(itor.eq.1) then
           temp = temp + &
                2.*int5(ri5_79,e(:,OP_1),f(:,OP_DRPP),g(:,OP_1),h(:,OP_1))
        endif
#endif
     end if

  case(1)


     if(surface_int) then
        temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
        temp = int4(e(:,OP_1),temp79a,norm79(:,2),h(:,OP_DR)) &
             - int4(e(:,OP_1),temp79a,norm79(:,1),h(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
        temp79a = e(:,OP_1)*h(:,OP_DP)
        temp = temp &
             + int5(ri2_79,temp79a,norm79(:,1),f(:,OP_DR),g(:,OP_DP)) &
             + int5(ri2_79,temp79a,norm79(:,2),f(:,OP_DZ),g(:,OP_DP)) &
             + int5(ri2_79,temp79a,norm79(:,1),f(:,OP_DRP),g(:,OP_1)) &
             + int5(ri2_79,temp79a,norm79(:,2),f(:,OP_DZP),g(:,OP_1))
#endif
     else
#if defined(USE3D) || defined(USECOMPLEX)
        temp79a = &
             (e(:,OP_DZ)*f(:,OP_DZPP) + e(:,OP_DR)*f(:,OP_DRPP))*g(:,OP_1) &
        + 2.*(e(:,OP_DZ)*f(:,OP_DZP) + e(:,OP_DR)*f(:,OP_DRP))*g(:,OP_DP) &
        +    (e(:,OP_DZ)*f(:,OP_DZ) + e(:,OP_DR)*f(:,OP_DR))*g(:,OP_DPP)
        temp = int3(ri2_79,temp79a,h(:,OP_1))
#endif
     end if


  end select
  v1ubb = temp
  return
end function v1ubb


! V1up
! ====
vectype function v1up(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = 0.
     end if

  case(1)
     if(itor.eq.0) then
        temp = 0.
     else
        if(surface_int) then
           temp = 2.* &
                (int5(r_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),norm79(:,2)) &
                -int5(r_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),norm79(:,2))) &
                +4.*gam*int4(e(:,OP_1),f(:,OP_DZ),g(:,OP_1),norm79(:,2))
        else
           temp = 2.* &
                (int4(r_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ)) &
                -int4(r_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DR))) &
                -4.*gam*int3(e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1))
        endif
     end if
  end select

  v1up = temp
  return
end function v1up



! V1vpsipsi
! =========
vectype function v1vpsipsi(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp79a = g(:,OP_DZ)* &
             ( h(:,OP_DZ)*e(:,OP_DRZ) + h(:,OP_DR)*e(:,OP_DRR)  &
             - e(:,OP_DZ)*h(:,OP_DRZ) - e(:,OP_DR)*h(:,OP_DRR)) &
             + g(:,OP_DR)* &
             ( e(:,OP_DR)*h(:,OP_DRZ) + e(:,OP_DZ)*h(:,OP_DZZ)  &
             - h(:,OP_DR)*e(:,OP_DRZ) - h(:,OP_DZ)*e(:,OP_DZZ))
        temp = -int3(ri3_79,f(:,OP_DP),temp79a)

        if(itor.eq.1) then
           temp79a = e(:,OP_DZ)* &
                (g(:,OP_DZ)*h(:,OP_DZ ) + g(:,OP_DR)*h(:,OP_DR)) &
             + 2.*e(:,OP_DR)*(g(:,OP_DZ)*h(:,OP_DR ))                         &
             +    e(:,OP_1) * &
                (g(:,OP_DR)*h(:,OP_DZZ) - g(:,OP_DZ)*h(:,OP_DRR))
           temp = temp &
                -2.*int3(ri4_79,f(:,OP_DP),temp79a) &
                -2.*int5(ri5_79,e(:,OP_1),f(:,OP_DP),g(:,OP_DZ),h(:,OP_DR))
        endif
     end if

  case(1)
     if(surface_int) then
        temp = 0.
     else
        temp79a = e(:,OP_DZ)*h(:,OP_DZP) + e(:,OP_DR)*h(:,OP_DRP)
        temp79b = f(:,OP_DP)*g(:,OP_GS) + f(:,OP_1)*g(:,OP_GSP) &
             + f(:,OP_DZP)*g(:,OP_DZ ) + f(:,OP_DRP)*g(:,OP_DR ) &
             + f(:,OP_DZ )*g(:,OP_DZP) + f(:,OP_DR )*g(:,OP_DRP)
        temp79c = &
             f(:,OP_DP)*(e(:,OP_DZ)*h(:,OP_DR ) - e(:,OP_DR)*h(:,OP_DZ )) &
          +  f(:,OP_1 )*(e(:,OP_DZ)*h(:,OP_DRP) - e(:,OP_DR)*h(:,OP_DZP))
        temp = int4(ri_79,g(:,OP_DZ),f(:,OP_DR),temp79a) &
             - int4(ri_79,g(:,OP_DR),f(:,OP_DZ),temp79a) &
             + int4(ri_79,e(:,OP_DZ),h(:,OP_DR),temp79b) &
             - int4(ri_79,e(:,OP_DR),h(:,OP_DZ),temp79b) &
             + int3(ri_79,g(:,OP_GS),temp79c)
        temp= -temp
     end if
  end select

  v1vpsipsi = temp
#else
  v1vpsipsi = 0.
#endif
  return
end function v1vpsipsi



! V1vpsib
! =======
vectype function v1vpsib(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  temp = 0.
  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        if(itor.eq.1) then
           temp = temp + &
             2.*(   int5(ri3_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
                -   int5(ri3_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1)) &
                +2.*int5(ri4_79,e(:,OP_DZ),f(:,OP_1 ),g(:,OP_DZ),h(:,OP_1)))
        endif
#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp - &
             (int5(ri4_79,e(:,OP_DZ),f(:,OP_DPP),g(:,OP_DZ),h(:,OP_1)) &
             +int5(ri4_79,e(:,OP_DR),f(:,OP_DPP),g(:,OP_DR),h(:,OP_1)))
        if(itor.eq.1) then
           temp = temp - &
                2.*int5(ri5_79,e(:,OP_1),f(:,OP_DPP),g(:,OP_DR),h(:,OP_1))
        endif
#endif
     end if

  case(1)

     if(surface_int) then
        temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
        temp = int4(e(:,OP_1),temp79a,norm79(:,2),h(:,OP_DR)) &
             - int4(e(:,OP_1),temp79a,norm79(:,1),h(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
        temp79a = e(:,OP_1)*h(:,OP_DP)
        temp = temp &
             - int5(ri2_79,temp79a,norm79(:,1),g(:,OP_DR),f(:,OP_DP)) &
             - int5(ri2_79,temp79a,norm79(:,2),g(:,OP_DZ),f(:,OP_DP)) &
             - int5(ri2_79,temp79a,norm79(:,1),g(:,OP_DRP),f(:,OP_1)) &
             - int5(ri2_79,temp79a,norm79(:,2),g(:,OP_DZP),f(:,OP_1))
#endif
     else
        temp = 0.

#if defined(USE3D) || defined(USECOMPLEX)
        temp79a = &
             f(:,OP_DPP)*(e(:,OP_DZ)*g(:,OP_DZ) + e(:,OP_DR)*g(:,OP_DR)) &
             +2.*f(:,OP_DP)*(e(:,OP_DZ)*g(:,OP_DZP) + e(:,OP_DR)*g(:,OP_DRP)) &
             +   f(:,OP_1)*(e(:,OP_DZ)*g(:,OP_DZPP) + e(:,OP_DR)*g(:,OP_DRPP))
        temp = -int3(ri2_79,h(:,OP_1),temp79a)
#endif
     end if
  end select

  v1vpsib = temp
  return
end function v1vpsib


! V1vp
! ====
vectype function v1vp(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = 0.
     end if
     
  case(1)

     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
     if(itor.eq.1) then
        if(surface_int) then
           temp = -2.* &
                (    int4(e(:,OP_1),f(:,OP_1),g(:,OP_DP),norm79(:,2)) &
                +gam*int4(e(:,OP_1),f(:,OP_DP),g(:,OP_1),norm79(:,2)))
        else
           temp79a =  f(:,OP_1 )*g(:,OP_DP) &
                + gam*f(:,OP_DP)*g(:,OP_1 )
           temp = 2.*int2(e(:,OP_DZ),temp79a)
        end if
     end if
#endif

  end select

  v1vp = temp
  return
end function v1vp



! V1chipsipsi
! ===========
vectype function v1chipsipsi(e,f,g,h)

  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e, f, g, h
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp79b = f(:,OP_DRZ)*g(:,OP_DZ ) + f(:,OP_DRR)*g(:,OP_DR ) &
             +    f(:,OP_DZ )*g(:,OP_DRZ) + f(:,OP_DR )*g(:,OP_DRR)
        temp79c = f(:,OP_DZZ)*g(:,OP_DZ ) + f(:,OP_DRZ)*g(:,OP_DR ) &
             +    f(:,OP_DZ )*g(:,OP_DZZ) + f(:,OP_DR )*g(:,OP_DRZ)
        
        temp79e = e(:,OP_DRZ)*h(:,OP_DR ) - e(:,OP_DRR)*h(:,OP_DZ ) &
             +    e(:,OP_DZ )*h(:,OP_DRR) - e(:,OP_DR )*h(:,OP_DRZ)
        temp79f = e(:,OP_DZZ)*h(:,OP_DR ) - e(:,OP_DRZ)*h(:,OP_DZ ) &
             +    e(:,OP_DZ )*h(:,OP_DRZ) - e(:,OP_DR )*h(:,OP_DZZ)

        temp = int3(ri_79,temp79b,temp79e) &
             + int3(ri_79,temp79c,temp79f) &
             + int4(ri_79,e(:,OP_DR),temp79c,h(:,OP_GS)) &
             - int4(ri_79,e(:,OP_DZ),temp79b,h(:,OP_GS))
     
        if(itor.eq.1) then
           temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
           
           temp = temp + 2.* &
                (int4(ri2_79,e(:,OP_1 ),temp79c,h(:,OP_GS )) &
                -int4(ri2_79,e(:,OP_1 ),temp79b,h(:,OP_DRZ)) &
                -int4(ri2_79,e(:,OP_1 ),temp79c,h(:,OP_DZZ)) &
                -int4(ri2_79,e(:,OP_DR),temp79b,h(:,OP_DZ )) &
                -int4(ri2_79,e(:,OP_DZ),temp79c,h(:,OP_DZ )) &
                +int4(ri2_79,e(:,OP_DZ),temp79b,h(:,OP_DR )) &
                -int4(ri2_79,e(:,OP_DR),temp79b,h(:,OP_DZ )))
        endif
     end if

  case(1)
     if(surface_int) then
        temp79a = e(:,OP_DR)*h(:,OP_DZ) - e(:,OP_DZ)*h(:,OP_DR)
        temp = int5(ri3_79,temp79a,norm79(:,1),f(:,OP_DRZ),g(:,OP_DZ )) &
             + int5(ri3_79,temp79a,norm79(:,1),f(:,OP_DRR),g(:,OP_DR )) &
             + int5(ri3_79,temp79a,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRZ)) &
             + int5(ri3_79,temp79a,norm79(:,1),f(:,OP_DR ),g(:,OP_DRR)) &
             + int5(ri3_79,temp79a,norm79(:,2),f(:,OP_DZZ),g(:,OP_DZ )) &
             + int5(ri3_79,temp79a,norm79(:,2),f(:,OP_DRZ),g(:,OP_DR )) &
             + int5(ri3_79,temp79a,norm79(:,2),f(:,OP_DZ ),g(:,OP_DZZ)) &
             + int5(ri3_79,temp79a,norm79(:,2),f(:,OP_DR ),g(:,OP_DRZ))
        if(itor.eq.1) then
           temp = temp - 2.* &
                (int5(ri4_79,temp79a,norm79(:,1),f(:,OP_DZ),g(:,OP_DZ)) &
                +int5(ri4_79,temp79a,norm79(:,1),f(:,OP_DR),g(:,OP_DR)))
        endif
     else
        temp79b = f(:,OP_DRZ)*g(:,OP_DZ ) + f(:,OP_DRR)*g(:,OP_DR ) &
             +    f(:,OP_DZ )*g(:,OP_DRZ) + f(:,OP_DR )*g(:,OP_DRR)
        temp79c = f(:,OP_DZZ)*g(:,OP_DZ ) + f(:,OP_DRZ)*g(:,OP_DR ) &
             +    f(:,OP_DZ )*g(:,OP_DZZ) + f(:,OP_DR )*g(:,OP_DRZ)
        
        temp79e = e(:,OP_DRZ)*h(:,OP_DR ) - e(:,OP_DRR)*h(:,OP_DZ ) &
             +    e(:,OP_DZ )*h(:,OP_DRR) - e(:,OP_DR )*h(:,OP_DRZ)
        temp79f = e(:,OP_DZZ)*h(:,OP_DR ) - e(:,OP_DRZ)*h(:,OP_DZ ) &
             +    e(:,OP_DZ )*h(:,OP_DRZ) - e(:,OP_DR )*h(:,OP_DZZ)


        temp = int3(ri3_79,temp79b,temp79e) &
             + int3(ri3_79,temp79c,temp79f) &
             + int4(ri3_79,e(:,OP_DR),temp79c,h(:,OP_GS)) &
             - int4(ri3_79,e(:,OP_DZ),temp79b,h(:,OP_GS))
     
        if(itor.eq.1) then
           temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
           temp79d = e(:,OP_DZ)*h(:,OP_DR) - e(:,OP_DR)*h(:,OP_DZ)
           
           temp = temp + &
                (2.*int4(ri4_79,e(:,OP_DZ),temp79a,h(:,OP_GS)) &
                -2.*int3(ri4_79,temp79e,temp79a) &
                +   int3(ri4_79,temp79d,temp79b) &
                -2.*int3(ri5_79,temp79d,temp79a))
        endif
     end if
  end select

  v1chipsipsi = temp
  return
end function v1chipsipsi


! V1chipsib
! =========
vectype function v1chipsib(e,f,g,h)
  
  use basic
  use arrays
  use m3dc1_nint

  implicit none
  
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e, f, g, h
  vectype :: temp

   temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
   select case(ivform)
   case(0)
      if(surface_int) then
         temp = 0.
      else
         temp79a = f(:,OP_DRP )*e(:,OP_DZZ)-f(:,OP_DZP )*e(:,OP_DRZ) &
              -    f(:,OP_DRRP)*e(:,OP_DR )-f(:,OP_DRZP)*e(:,OP_DZ )
         temp79b = f(:,OP_DZP )*e(:,OP_DRR)-f(:,OP_DRP )*e(:,OP_DRZ) &
              -    f(:,OP_DRZP)*e(:,OP_DR )-f(:,OP_DZZP)*e(:,OP_DZ )

         temp = temp +   &
              int4(ri2_79,g(:,OP_DR),temp79a,h(:,OP_1)) &
              +int4(ri2_79,g(:,OP_DZ),temp79b,h(:,OP_1))

         if(itor.eq.1) then
            temp = temp &
             -   int5(ri3_79,f(:,OP_DRP), e(:,OP_DR),g(:,OP_DR),h(:,OP_1))  &
             -2.*int5(ri3_79,f(:,OP_DRRP),e(:,OP_1), g(:,OP_DR),h(:,OP_1))  &
             -2.*int5(ri4_79,f(:,OP_DRP), e(:,OP_1), g(:,OP_DR),h(:,OP_1))  &
             +   int5(ri3_79,f(:,OP_DZP), e(:,OP_DR),g(:,OP_DZ),h(:,OP_1))  &
             -2.*int5(ri3_79,f(:,OP_DRP), e(:,OP_DZ),g(:,OP_DZ),h(:,OP_1))  &
             -2.*int5(ri3_79,f(:,OP_DRZP),e(:,OP_1), g(:,OP_DZ),h(:,OP_1))
         endif
      end if

  case(1)
     if(surface_int) then
        temp79a = f(:,OP_DRP)*g(:,OP_DR ) + f(:,OP_DZP)*g(:,OP_DZ ) &
             +    f(:,OP_DR )*g(:,OP_DRP) + f(:,OP_DZ )*g(:,OP_DZP)
        temp79b = e(:,OP_1)*h(:,OP_DP)
        temp = &
             - int5(ri4_79,norm79(:,1),e(:,OP_DR),temp79a,h(:,OP_1)) &
             - int5(ri4_79,norm79(:,2),e(:,OP_DZ),temp79a,h(:,OP_1)) &
             - int5(ri4_79,norm79(:,1),temp79b,f(:,OP_DRR),g(:,OP_DR )) &
             - int5(ri4_79,norm79(:,1),temp79b,f(:,OP_DRZ),g(:,OP_DZ )) &
             - int5(ri4_79,norm79(:,1),temp79b,f(:,OP_DR ),g(:,OP_DRR)) &
             - int5(ri4_79,norm79(:,1),temp79b,f(:,OP_DZ ),g(:,OP_DRZ)) &
             - int5(ri4_79,norm79(:,2),temp79b,f(:,OP_DRZ),g(:,OP_DR )) &
             - int5(ri4_79,norm79(:,2),temp79b,f(:,OP_DZZ),g(:,OP_DZ )) &
             - int5(ri4_79,norm79(:,2),temp79b,f(:,OP_DR ),g(:,OP_DRZ)) &
             - int5(ri4_79,norm79(:,2),temp79b,f(:,OP_DZ ),g(:,OP_DZZ))
        if(itor.eq.1) then
           temp = temp + 2.* &
                (int5(ri5_79,temp79b,norm79(:,1),f(:,OP_DZ),g(:,OP_DZ)) &
                +int5(ri5_79,temp79b,norm79(:,1),f(:,OP_DR),g(:,OP_DR)))
        endif
     else
        temp79a = h(:,OP_DZP)*f(:,OP_DR ) - h(:,OP_DRP)*f(:,OP_DZ ) &
             +    h(:,OP_DZ )*f(:,OP_DRP) - h(:,OP_DR )*f(:,OP_DZP)
        temp79b = (e(:,OP_DZ)*f(:,OP_DZ )+e(:,OP_DR)*f(:,OP_DR ))*h(:,OP_DP) &
             +    (e(:,OP_DZ)*f(:,OP_DZP)+e(:,OP_DR)*f(:,OP_DRP))*h(:,OP_1 )
        temp79c = f(:,OP_DZP)*g(:,OP_DZ ) + f(:,OP_DRP)*g(:,OP_DR ) &
             +    f(:,OP_DZ )*g(:,OP_DZP) + f(:,OP_DR )*g(:,OP_DRP)
        temp79d = h(:,OP_1)*f(:,OP_GS) &
             + h(:,OP_DZ)*f(:,OP_DZ) + h(:,OP_DR)*f(:,OP_DR)
        if(itor.eq.1) then
           temp79a = temp79a + 4.*ri_79* &
                (h(:,OP_DP)*f(:,OP_DZ) + h(:,OP_1)*f(:,OP_DZP))
           temp79d = temp79d - 2.*ri_79*h(:,OP_1)*f(:,OP_DR)
        endif
     
        temp = int4(ri4_79,e(:,OP_DZ),g(:,OP_DR),temp79a)  &
             - int4(ri4_79,e(:,OP_DR),g(:,OP_DZ),temp79a)  &
             - int3(ri4_79,g(:,OP_GS),temp79b)             &
             - int4(ri4_79,e(:,OP_GS),h(:,OP_1 ),temp79c)  &
             - int4(ri4_79,e(:,OP_DZ),h(:,OP_DZ),temp79c)  &
             - int4(ri4_79,e(:,OP_DR),h(:,OP_DR),temp79c)  &
             + int4(ri4_79,e(:,OP_DZ),g(:,OP_DZP),temp79d) &
             + int4(ri4_79,e(:,OP_DR),g(:,OP_DRP),temp79d)
        temp = -temp
     end if
  end select
#endif

  v1chipsib = temp
  return
end function v1chipsib


! V1chibb
! =======
vectype function v1chibb(e,f,g,h)

  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e, f, g, h

  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else

        temp = 0.
        if(itor.eq.1) then
           temp = temp -2.* &
                (int5(ri2_79,e(:,OP_DZ),g(:,OP_1 ),f(:,OP_GS),h(:,OP_1)) &
                +int5(ri2_79,e(:,OP_DZ),g(:,OP_DZ),f(:,OP_DZ),h(:,OP_1)) &
                +int5(ri2_79,e(:,OP_DZ),g(:,OP_DR),f(:,OP_DR),h(:,OP_1)))
        endif

#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp &
             + int5(ri3_79,f(:,OP_DZPP),e(:,OP_DR),g(:,OP_1),h(:,OP_1)) &
             - int5(ri3_79,f(:,OP_DRPP),e(:,OP_DZ),g(:,OP_1),h(:,OP_1))
        
        if(itor.eq.1) then
           temp = temp + &
                2.*int5(ri4_79,f(:,OP_DZPP),e(:,OP_1),g(:,OP_1),h(:,OP_1))
        endif
#endif
     end if

  case(1)
     if(surface_int) then
        temp79a = e(:,OP_1)*(norm79(:,1)*h(:,OP_DZ) - norm79(:,2)*h(:,OP_DR))
        temp = int4(ri3_79,temp79a,f(:,OP_GS),g(:,OP_1 )) &
             + int4(ri3_79,temp79a,f(:,OP_DZ),g(:,OP_DZ)) &
             + int4(ri3_79,temp79a,f(:,OP_DR),g(:,OP_DR)) 
        if(itor.eq.1) then
           temp = temp &
                - 2.*int4(ri4_79,temp79a,f(:,OP_DR),g(:,OP_1))
        endif
#if defined(USE3D) || defined(USECOMPLEX)
        temp79b = e(:,OP_1)*h(:,OP_DP)
        temp = temp &
             + int5(ri5_79,temp79b,g(:,OP_1 ),norm79(:,1),f(:,OP_DZP)) &
             - int5(ri5_79,temp79b,g(:,OP_1 ),norm79(:,2),f(:,OP_DRP)) &
             + int5(ri5_79,temp79b,g(:,OP_DP),norm79(:,1),f(:,OP_DZ )) &
             - int5(ri5_79,temp79b,g(:,OP_DP),norm79(:,2),f(:,OP_DR ))
#endif
     else
        temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
        temp79a = &
             (e(:,OP_DZ)*f(:,OP_DR) - e(:,OP_DR)*f(:,OP_DZ))*g(:,OP_DPP) &
        + 2.*(e(:,OP_DZ)*f(:,OP_DRP) - e(:,OP_DR)*f(:,OP_DZP))*g(:,OP_DP) &
        +    (e(:,OP_DZ)*f(:,OP_DRPP) - e(:,OP_DR)*f(:,OP_DZPP))*g(:,OP_1)
        temp = -int3(ri5_79,temp79a,h(:,OP_1))
#endif
     end if

  end select
  
  v1chibb = temp
  return
end function v1chibb


! V1chip
! ======
vectype function v1chip(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = 0.
     end if
     
  case(1)
     if(itor.eq.0) then
        temp = 0.
     else
        temp79a = gam*f(:,OP_GS)*g(:,OP_1) + &
             f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)

        if(surface_int) then
           temp = -2.*int4(ri2_79,e(:,OP_1),temp79a,norm79(:,2))
        else
           temp = 2.*int3(ri2_79,e(:,OP_DZ),temp79a)
        end if
     end if
  end select

  v1chip = temp
  return
end function v1chip


! V1ngrav
! =======
vectype function v1ngrav(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f
  vectype :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v1ngrav = 0.
     return
  endif

  if(surface_int) then
     temp = 0.
  else
     temp = gravz*int3( r_79,e(:,OP_1),f(:,OP_DR)) &
          - gravr*int3(ri_79,e(:,OP_1),f(:,OP_DZ))
  end if

  v1ngrav = temp
  return
end function v1ngrav


! V1ungrav
! ========
vectype function v1ungrav(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v1ungrav = 0.
     return
  endif

  if(surface_int) then
     temp = 0.
  else
     temp79a = f(:,OP_DR)*g(:,OP_DZ) - f(:,OP_DZ)*g(:,OP_DR)
  
     temp = gravz*int2(       e(:,OP_DR),temp79a) &
          - gravr*int3(ri2_79,e(:,OP_DZ),temp79a)
     
     if(itor.eq.1) &
          temp = temp + 2.*gravz*int3(ri_79,e(:,OP_1),temp79a)
  end if

  v1ungrav = temp
  return
end function v1ungrav


! V1chingrav
! ==========
vectype function v1chingrav(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v1chingrav = 0.
     return
  endif

  if(surface_int) then
     temp = 0.
  else
     temp79a = r_79*(f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR) &
          + g(:,OP_1)*f(:,OP_LP))

     temp = gravz*int2(       e(:,OP_DR),temp79a) &
          - gravr*int3(ri2_79,e(:,OP_DZ),temp79a)

     if(itor.eq.1) &
          temp = temp + 2.*gravz*int3(ri_79,e(:,OP_1),temp79a)
  endif

  v1chingrav = temp
  return
end function v1chingrav


! V1ndenmgrav
! ===========
vectype function v1ndenmgrav(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f
  real, intent(in) :: g
  vectype :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v1ndenmgrav = 0.
     return
  endif

  if(surface_int) then
     temp = 0.
  else
     temp79a = -g*r_79*f(:,OP_LP)

     temp = gravz*int2(       e(:,OP_DR),temp79a) &
          - gravr*int3(ri2_79,e(:,OP_DZ),temp79a)

     if(itor.eq.1) &
          temp = temp + 2.*gravz*int3(ri_79,e(:,OP_1),temp79a)
  end if

  v1ndenmgrav = temp
  return
end function v1ndenmgrav


! V1us
! ====
vectype function v1us(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(idens.eq.0 .or. nosig.eq.1) then
     v1us = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) + denm*nt79(:,OP_LP)

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = &
             - int3(e(:,OP_DZ),f(:,OP_DZ),temp79a) &
             - int3(e(:,OP_DR),f(:,OP_DR),temp79a)

        if(itor.eq.1) then 
           temp = temp - 2.*int4(ri_79,e(:,OP_1),f(:,OP_DR),temp79a)
        endif
     end if

  case(1)
     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp = 0.
        else
           temp =  &
                + int5(r2_79,e(:,OP_1),norm79(:,1),f(:,OP_DR),temp79a) &
                + int5(r2_79,e(:,OP_1),norm79(:,2),f(:,OP_DZ),temp79a)
        endif
     else
        temp = -int4(r2_79,e(:,OP_DZ),f(:,OP_DZ),temp79a) &
             -  int4(r2_79,e(:,OP_DR),f(:,OP_DR),temp79a)
     end if
  end select

  v1us = temp
  return
end function v1us


! V1chis
! ======
vectype function v1chis(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(idens.eq.0 .or. nosig.eq.1) then
     v1chis = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) + denm*nt79(:,OP_LP)


  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int4(r_79,e(:,OP_DZ),f(:,OP_DR),temp79a) &
             - int4(r_79,e(:,OP_DR),f(:,OP_DZ),temp79a)

        if(itor.eq.1) then 
           temp = temp - 2.*int3(e(:,OP_1),f(:,OP_DZ),temp79a)
        endif
     end if

  case(1)
     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp = 0.
        else
           temp = &
                + int5(ri_79,e(:,OP_1),norm79(:,1),f(:,OP_DZ),temp79a) &
                - int5(ri_79,e(:,OP_1),norm79(:,2),f(:,OP_DR),temp79a)
        endif
     else
        temp = int4(ri3_79,e(:,OP_DZ),f(:,OP_DR),temp79a) &
             - int4(ri3_79,e(:,OP_DR),f(:,OP_DZ),temp79a)
     endif
  end select

  v1chis = temp
  return
end function v1chis


! V1psif
! ======
vectype function v1psif(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
  select case (ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = &
             - int3(e(:,OP_DZ),f(:,OP_GS),g(:,OP_DZP)) &
             - int3(e(:,OP_DR),f(:,OP_GS),g(:,OP_DRP))
        if(itor.eq.1) then
           temp = temp - 2.* & 
                int4(ri_79,e(:,OP_1),f(:,OP_GS),g(:,OP_DRP))
        endif
     end if

  case(1)

     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp = 0.
        else
           temp =  &
                + int4(e(:,OP_1),f(:,OP_GS),norm79(:,1),g(:,OP_DRP)) &
                + int4(e(:,OP_1),f(:,OP_GS),norm79(:,2),g(:,OP_DZP))
        endif
     else
        temp = &
             - int3(e(:,OP_DZ),f(:,OP_GS),g(:,OP_DZP)) &
             - int3(e(:,OP_DR),f(:,OP_GS),g(:,OP_DRP))
     end if
  end select
  
#else
  temp = 0.
#endif
  v1psif = temp
  return
end function v1psif


! V1bf
! ====
vectype function v1bf(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
#if defined(USE3D) || defined(USECOMPLEX)
  select case (ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = &
             + int4(ri_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DRPP)) &
             - int4(ri_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DZPP))
        if(itor.eq.1) then
           temp = temp - 2.* & 
                int4(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZPP))
        endif
     end if

  case(1)
     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp = 0.
        else
           temp = &
                + int5(ri_79,e(:,OP_1),f(:,OP_1),norm79(:,1),g(:,OP_DZPP)) &
                - int5(ri_79,e(:,OP_1),f(:,OP_1),norm79(:,2),g(:,OP_DRPP))
        end if
     else
        temp = &
             + int4(ri_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DRPP)) &
             - int4(ri_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DZPP))
     end if
  end select
#else
  temp = 0.
#endif
  v1bf = temp
  return
end function v1bf


! V1p
! ===
vectype function v1p(e,f)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e, f

  vectype :: temp
  temp = 0.

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = 0.
     end if

  case(1)
     if(surface_int) then
        if(inoslip_pol.eq.1) then
           temp = 0.
        else
           temp = &
             + int4(r_79,e(:,OP_1),norm79(:,1),f(:,OP_DZ)) &
             - int4(r_79,e(:,OP_1),norm79(:,2),f(:,OP_DR))
        endif
     else
        temp = int3(r_79,e(:,OP_DZ),f(:,OP_DR)) &
             - int3(r_79,e(:,OP_DR),f(:,OP_DZ))
     end if
  end select

  v1p = temp
  return
end function v1p



!============================================================================
! V2 TERMS
!============================================================================


! V2vn
! ====
vectype function v2vn(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(surface_int) then
     v2vn = 0.
     return
  endif

  select case(ivform)
  case(0)
     temp = int3(e(:,OP_1),f(:,OP_1),g(:,OP_1))
  case(1)
     temp = int4(r2_79,e(:,OP_1),f(:,OP_1),g(:,OP_1))
  end select

  v2vn = temp
  return
end function v2vn


! V2umu
! =====
vectype function v2umu(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
  select case(ivform)
  case(0)
     ! not imeplemented
     temp = 0.

  case(1)
     if(surface_int) then
        temp = int5(r_79,e(:,OP_1),norm79(:,2),f(:,OP_DRP),g(:,OP_1)) &
             - int5(r_79,e(:,OP_1),norm79(:,1),f(:,OP_DZP),g(:,OP_1))
     else
        temp = int4(r_79,e(:,OP_DR),f(:,OP_DZP),g(:,OP_1)) &
             - int4(r_79,e(:,OP_DZ),f(:,OP_DRP),g(:,OP_1))
        if(itor.eq.1) then
           temp = temp &
                + 2.*int3(e(:,OP_1),f(:,OP_DZP),g(:,OP_1)) &
                - 4.*int3(e(:,OP_1),f(:,OP_DZP),h(:,OP_1))
        endif
     end if
  end select
#endif
  v2umu = temp
  return
end function v2umu


! V2vmu
! =====
vectype function v2vmu(e,f,g,h,i)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int3(e(:,OP_1),f(:,OP_GS),g(:,OP_1)) &
             + int3(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ)) &
             + int3(e(:,OP_1),f(:,OP_DR),g(:,OP_DR))
     
        if(itor.eq.1) then
           temp = temp - 2.*int4(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR))
        endif

#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp + 2.*int4(ri2_79,e(:,OP_1),f(:,OP_DPP),h(:,OP_1))
#endif

        ! hyperviscous
        if(hypv.ne.0.) then
           if(ihypamu.eq.1) then
              temp79a = e(:,OP_GS)*g(:,OP_1) + &
                   e(:,OP_DZ)*g(:,OP_DZ) + e(:,OP_DR)*g(:,OP_DR)
              if(itor.eq.1) temp79a = temp79a + 4.*ri_79*e(:,OP_DR)*g(:,OP_1)
           else
              temp79a = e(:,OP_GS)
              if(itor.eq.1) temp79a = temp79a + 4.*ri_79*e(:,OP_DR)
           endif
           temp = temp - int3(temp79a,f(:,OP_GS),i(:,OP_1))  
        endif
     end if

  case(1)
     if(surface_int) then
        temp = int5(r2_79,e(:,OP_1),norm79(:,1),f(:,OP_DR),g(:,OP_1)) &
             + int5(r2_79,e(:,OP_1),norm79(:,2),f(:,OP_DZ),g(:,OP_1))
     else
        temp = -int4(r2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1)) &
             -  int4(r2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_1))
     
#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp + 2.*int3(e(:,OP_1),f(:,OP_DPP),h(:,OP_1))
#endif

        ! hyperviscous
        if(hypv.ne.0.) then
           if(ihypamu.eq.1) then
              temp79a = e(:,OP_GS)*g(:,OP_1) + &
                   e(:,OP_DZ)*g(:,OP_DZ) + e(:,OP_DR)*g(:,OP_DR)
              if(itor.eq.1) temp79a = temp79a + 4.*ri_79*e(:,OP_DR)*g(:,OP_1)
           else
              temp79a = e(:,OP_GS)
              if(itor.eq.1) temp79a = temp79a + 4.*ri_79*e(:,OP_DR)
           endif
        
           temp = temp - int4(r2_79,temp79a,f(:,OP_GS),i(:,OP_1))
           if(itor.eq.1) then
              temp = temp - 4.*int4(r_79,temp79a,f(:,OP_DR),i(:,OP_1))
           endif
        end if
     end if
  end select

  v2vmu = temp
  return
end function v2vmu


! V2chimu
! =======
vectype function v2chimu(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
  select case(ivform)
  case(0)
     ! not implemented
     temp = 0.

  case(1)
     if(surface_int) then
        temp = int5(ri2_79,e(:,OP_1),norm79(:,1),f(:,OP_DRP),g(:,OP_1)) &
             + int5(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_DZP),g(:,OP_1))
     else
        temp79a = h(:,OP_1) - g(:,OP_1)
        temp = &
             - int4(ri2_79,e(:,OP_DZ),f(:,OP_DZP),g(:,OP_1)) &
             - int4(ri2_79,e(:,OP_DR),f(:,OP_DRP),g(:,OP_1)) &
             + 2.*int4(ri2_79,e(:,OP_1),f(:,OP_GSP),temp79a)
        if(itor.eq.1) then
           temp = temp &
                +2.*int4(ri3_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_1))
        endif
     end if
  end select
#endif
  v2chimu = temp
  return
end function v2chimu



! V2vun
! =====
vectype function v2vun(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(inertia.eq.0 .or. surface_int) then
     v2vun = 0.
     return
  end if


  select case(ivform)
  case(0)
     temp = int5(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1)) &
          - int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1))
  case(1)
     temp = int5(r3_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1)) &
          - int5(r3_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1))

     if(itor.eq.1) then
        temp = temp + &
             2.*int5(r2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),h(:,OP_1))
     end if

  end select

  v2vun = temp
  return
end function v2vun


! V2vvn
! =====
vectype function v2vvn(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(inertia.eq.0 .or. surface_int) then
     v2vvn = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  select case(ivform)
  case(0)
     temp = -int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DP),h(:,OP_1))
  case(1)
     temp = -int5(r2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DP),h(:,OP_1))
  end select
#else
  temp = 0.
#endif
  v2vvn = temp
end function v2vvn

! V2up
! ====
vectype function v2up(e,f,g)
  use basic
  use arrays
  use m3dc1_nint
  
  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e, f, g

  vectype :: temp
  temp = 0.

  if(surface_int) then
     v2up = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  select case(ivform)
  case(0)
     temp = & 
          -  int4(ri_79,e(:,OP_1),f(:,OP_DZP),g(:,OP_DR)) &
          +  int4(ri_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_DZ)) &
          -  int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DRP)) &
          +  int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZP))

  case(1)
     temp = int4(r_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_DZ)) &
          - int4(r_79,e(:,OP_1),f(:,OP_DZP),g(:,OP_DR)) &
          + int4(r_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZP)) &
          - int4(r_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DRP))
     if(itor.eq.1) then
        temp = temp - 2.*gam* &
             (int3(e(:,OP_1),f(:,OP_DZP),g(:,OP_1)) &
             +int3(e(:,OP_1),f(:,OP_DZ),g(:,OP_DP)))
     endif
  end select
#endif

  v2up = temp
  return
end function v2up

! V2vp
! ====
vectype function v2vp(e,f,g)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e, f, g

  vectype :: temp

  if(surface_int) then
     v2vp = 0.
     return
  end if

  temp = 0.

#if defined(USE3D) || defined(USECOMPLEX)
  select case(ivform)
  case(0)
     temp =          int4(ri2_79,e(:,OP_1),g(:,OP_DPP),f(:,OP_1)) &
          + (1.+gam)*int4(ri2_79,e(:,OP_1),g(:,OP_DP),f(:,OP_DP)) &
          + gam*     int4(ri2_79,e(:,OP_1),g(:,OP_1),f(:,OP_DPP))
  case(1)
     temp =          int3(e(:,OP_1),g(:,OP_DPP),f(:,OP_1)) &
          + (1.+gam)*int3(e(:,OP_1),g(:,OP_DP),f(:,OP_DP)) &
          + gam*     int3(e(:,OP_1),g(:,OP_1),f(:,OP_DPP))
  end select
#endif

  v2vp = temp
  return
end function v2vp


! V2chip
! ======
vectype function v2chip(e,f,g)

  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e, f, g
  vectype :: temp

  if(surface_int) then
     v2chip = 0.
     return
  end if

  temp = 0.

#if defined(USE3D) || defined(USECOMPLEX)
  select case (ivform)
  case(0)
     temp =     int3(e(:,OP_1),f(:,OP_DRP),g(:,OP_DR ))    &
          +     int3(e(:,OP_1),f(:,OP_DZP),g(:,OP_DZ ))    &
          +     int3(e(:,OP_1),f(:,OP_DR ),g(:,OP_DRP))    &
          +     int3(e(:,OP_1),f(:,OP_DZ ),g(:,OP_DZP))    &
          + gam*int3(e(:,OP_1),f(:,OP_LPP),g(:,OP_1 ))     &
          + gam*int3(e(:,OP_1),f(:,OP_LP ),g(:,OP_DP))

  case(1)
     temp =     int4(ri2_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_DR))    &
              + int4(ri2_79,e(:,OP_1),f(:,OP_DZP),g(:,OP_DZ))    &
          + gam*int4(ri2_79,e(:,OP_1),f(:,OP_GSP),g(:,OP_1 ))    &
              + int4(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DRP))    &
              + int4(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZP))    &
          + gam*int4(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_DP))
  end select
#endif

  v2chip = temp
  return
end function v2chip

! V2p
! ===
vectype function v2p(e,f)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e, f

  vectype :: temp

  if(surface_int) then
     v2p = 0.
     return
  end if

  temp = 0.

#if defined(USE3D) || defined(USECOMPLEX)
  ! same for both ivforms
  temp = -int2(e(:,OP_1),f(:,OP_DP))
#endif

  v2p = temp
  return
end function v2p


! V2psipsi
! ========
vectype function v2psipsi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(surface_int) then
     v2psipsi = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  ! same for both ivforms
  temp = - &
       (int4(ri2_79,e(:,OP_1),f(:,OP_DZP),g(:,OP_DZ)) &
       +int4(ri2_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_DR)))
#else
  temp = 0.
#endif
  v2psipsi = temp
  return
end function v2psipsi


! V2psib
! ======
vectype function v2psib(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(surface_int) then
     v2psib = 0.
     return
  end if

  temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ)) &
       - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR))

  v2psib = temp
  return
end function v2psib


! V2vpsipsi
! =========
vectype function v2vpsipsi(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  ! [nu,psi(2)]
  temp79a = e(:,OP_DZ)*h(:,OP_DR) - e(:,OP_DR)*h(:,OP_DZ)

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int4(ri2_79,f(:,OP_DR),g(:,OP_DZ),temp79a) &
             - int4(ri2_79,f(:,OP_DZ),g(:,OP_DR),temp79a)
        if(itor.eq.1) then
           temp = temp - 2.*int4(ri3_79,f(:,OP_1),g(:,OP_DZ),temp79a)
        endif

#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp + &
             (int5(ri4_79,e(:,OP_1),f(:,OP_DPP),g(:,OP_DZ),h(:,OP_DZ)) &
             +int5(ri4_79,e(:,OP_1),f(:,OP_DPP),g(:,OP_DR),h(:,OP_DR)))
#endif
     end if

  case(1)
     if(surface_int) then
        temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
        temp = int4(e(:,OP_1),temp79a,norm79(:,2),h(:,OP_DR)) &
             - int4(e(:,OP_1),temp79a,norm79(:,1),h(:,OP_DZ))
     else
        temp = int3(f(:,OP_DR),g(:,OP_DZ),temp79a) &
             - int3(f(:,OP_DZ),g(:,OP_DR),temp79a)

#if defined(USE3D) || defined(USECOMPLEX)
        temp79b = &
              f(:,OP_DPP)*(g(:,OP_DZ )*h(:,OP_DZ ) + g(:,OP_DR )*h(:,OP_DR )) &
          + 2.*f(:,OP_DP)*(g(:,OP_DZP)*h(:,OP_DZ ) + g(:,OP_DRP)*h(:,OP_DR )) &
          +   f(:,OP_DP )*(g(:,OP_DZ )*h(:,OP_DZP) + g(:,OP_DR )*h(:,OP_DRP)) &
          + f(:,OP_1 )*(g(:,OP_DZPP)*h(:,OP_DZ ) + g(:,OP_DRPP)*h(:,OP_DR )) &
          + f(:,OP_1 )*(g(:,OP_DZP )*h(:,OP_DZP) + g(:,OP_DRP )*h(:,OP_DRP))
        temp = temp + int3(ri2_79,e(:,OP_1),temp79b)
#endif
     end if

  end select

  v2vpsipsi = temp
  return
end function v2vpsipsi


! V2vpsib
! =======
vectype function v2vpsib(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp79a = f(:,OP_DP)*(g(:,OP_DZ )*h(:,OP_DR)-g(:,OP_DR )*h(:,OP_DZ)) &
             +    f(:,OP_1 )*(g(:,OP_DZP)*h(:,OP_DR)-g(:,OP_DRP)*h(:,OP_DZ))
        temp = int3(ri3_79,e(:,OP_1),temp79a) &
             + int5(ri3_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_DP)) &
             - int5(ri3_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_DP))
        if(itor.eq.1) then
           temp = temp + 2.* &
                int5(ri4_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),h(:,OP_DP))
        endif
     end if

  case(1)
     if(surface_int) then
        temp = 0.
     else
        temp79a = f(:,OP_DP)*(g(:,OP_DZ )*h(:,OP_DR)-g(:,OP_DR )*h(:,OP_DZ)) &
             +    f(:,OP_1 )*(g(:,OP_DZP)*h(:,OP_DR)-g(:,OP_DRP)*h(:,OP_DZ))
        temp = int3(ri_79,e(:,OP_1),temp79a)
     end if
  end select
  v2vpsib = temp
#else
  v2vpsib = 0.
#endif
  return
end function v2vpsib



! V2upsipsi
! =========
vectype function v2upsipsi(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)  
  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp79a = h(:,OP_DZ)* &
             (f(:,OP_DZZP)*g(:,OP_DR )-f(:,OP_DRZP)*g(:,OP_DZ )  &
             +f(:,OP_DZP )*g(:,OP_DRZ)-f(:,OP_DRP )*g(:,OP_DZZ)) &
             +    h(:,OP_DR)* &
             (f(:,OP_DRZP)*g(:,OP_DR )-f(:,OP_DRRP)*g(:,OP_DZ )  &
             +f(:,OP_DZP )*g(:,OP_DRR)-f(:,OP_DRP )*g(:,OP_DRZ))
        temp = -int3(ri3_79,e(:,OP_1),temp79a)

        if(itor.eq.1) then
           temp79a = (f(:,OP_DZP)*g(:,OP_DR)-f(:,OP_DRP)*g(:,OP_DZ))*h(:,OP_DR)
           temp = temp + int3(ri4_79,e(:,OP_1),temp79a)
        endif
     end if
 
  case(1)
     temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
     temp79b = f(:,OP_DZP)*g(:,OP_DR ) - f(:,OP_DRP)*g(:,OP_DZ ) &
          +    f(:,OP_DZ )*g(:,OP_DRP) - f(:,OP_DR )*g(:,OP_DZP)

     if(surface_int) then
        temp = &
             - int5(ri_79,e(:,OP_1),temp79a,norm79(:,1),h(:,OP_DRP)) &
             - int5(ri_79,e(:,OP_1),temp79a,norm79(:,2),h(:,OP_DZP)) &
             - int5(ri_79,e(:,OP_1),temp79b,norm79(:,1),h(:,OP_DR)) &
             - int5(ri_79,e(:,OP_1),temp79b,norm79(:,2),h(:,OP_DZ))
     else
        temp79c = e(:,OP_1)*h(:,OP_GS) &
             +    e(:,OP_DZ)*h(:,OP_DZ) + e(:,OP_DR)*h(:,OP_DR)
        temp79d = e(:,OP_1)*h(:,OP_GSP) &
             +    e(:,OP_DZ)*h(:,OP_DZP) + e(:,OP_DR)*h(:,OP_DRP)
        temp = int3(ri_79,temp79a,temp79d) &
             + int3(ri_79,temp79b,temp79c)
     end if
  end select
#else
  temp = 0.
#endif

  v2upsipsi = temp
  return
end function v2upsipsi



! V2upsib
! =======
vectype function v2upsib(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  temp79a = e(:,OP_DZ)*f(:,OP_DR) - e(:,OP_DR)*f(:,OP_DZ)
  select case(ivform)
  case(0)

     if(surface_int) then
        temp = 0.
     else
        temp = int4(ri2_79,temp79a,g(:,OP_DR),h(:,OP_DZ)) &
             - int4(ri2_79,temp79a,g(:,OP_DZ),h(:,OP_DR))
        if(itor.eq.1) then
           temp = temp + 2.* &
                (int5(ri3_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
                -int5(ri3_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)))
        endif

#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp - &
             (int5(ri4_79,e(:,OP_1),f(:,OP_DZPP),g(:,OP_DZ),h(:,OP_1)) &
             +int5(ri4_79,e(:,OP_1),f(:,OP_DRPP),g(:,OP_DR),h(:,OP_1)))
#endif
     end if

  case(1)
     if(surface_int) then
        temp79a = e(:,OP_1)*(h(:,OP_DZ)*g(:,OP_DR) - h(:,OP_DR)*g(:,OP_DZ))
        temp = int3(temp79a,norm79(:,1),f(:,OP_DZ)) &
             - int3(temp79a,norm79(:,2),f(:,OP_DR))
     else
        temp = (int3(temp79a,g(:,OP_DR),h(:,OP_DZ)) &
             -  int3(temp79a,g(:,OP_DZ),h(:,OP_DR)))

#if defined(USE3D) || defined(USECOMPLEX)
        temp79b = &
          2.*(f(:,OP_DZP)*g(:,OP_DZ ) + f(:,OP_DRP)*g(:,OP_DR ))*h(:,OP_DP ) &
          +  (f(:,OP_DZ )*g(:,OP_DZP) + f(:,OP_DR )*g(:,OP_DRP))*h(:,OP_DP ) &
          +  (f(:,OP_DZ )*g(:,OP_DZ ) + f(:,OP_DR )*g(:,OP_DR ))*h(:,OP_DPP) &
          +  (f(:,OP_DZPP)*g(:,OP_DZ ) + f(:,OP_DRPP)*g(:,OP_DR ))*h(:,OP_1) &
          +  (f(:,OP_DZP )*g(:,OP_DZP) + f(:,OP_DRP )*g(:,OP_DRP))*h(:,OP_1)
        temp = temp - int3(ri2_79,e(:,OP_1),temp79b)
#endif
     end if

  end select

  v2upsib = temp
  return
end function v2upsib


! V2ubb
! =====
vectype function v2ubb(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = &
             (int5(ri3_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_DZ),h(:,OP_1)) &
             -int5(ri3_79,e(:,OP_1),f(:,OP_DZP),g(:,OP_DR),h(:,OP_1)))
        if(itor.eq.1) then
           temp = temp - 2.* &
                int5(ri4_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_1),h(:,OP_DP))
        endif
     end if

  case(1)
     if(surface_int) then
        temp = 0.
     else
        temp79a = h(:,OP_DP)*(f(:,OP_DR)*g(:,OP_DZ) - f(:,OP_DZ)*g(:,OP_DR)) &
             +    h(:,OP_1)*(f(:,OP_DRP)*g(:,OP_DZ) - f(:,OP_DZP)*g(:,OP_DR))
        temp = int3(ri_79,e(:,OP_1),temp79a)
     end if
  end select

#else
  temp = 0.
#endif
  v2ubb = temp
  return
end function v2ubb



! v2upsisb2
! ========
vectype function v2upsisb2(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

  v2upsisb2 = 0.
  return
end function v2upsisb2


! v2ubsb1
! =======
vectype function v2ubsb1(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

  v2ubsb1 = 0.
  return
end function v2ubsb1


! v2chipsipsi
! ===========
vectype function v2chipsipsi(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = temp  &
             +int5(ri2_79,f(:,OP_DRRP),e(:,OP_1),g(:,OP_DR ),h(:,OP_DR)) &
             +int5(ri2_79,f(:,OP_DRP ),e(:,OP_1),g(:,OP_DRR),h(:,OP_DR)) &
             +int5(ri2_79,f(:,OP_DRZP),e(:,OP_1),g(:,OP_DZ ),h(:,OP_DR)) &
             +int5(ri2_79,f(:,OP_DZP ),e(:,OP_1),g(:,OP_DRZ),h(:,OP_DR)) &
             +int5(ri2_79,f(:,OP_DRZP),e(:,OP_1),g(:,OP_DR ),h(:,OP_DZ)) &
             +int5(ri2_79,f(:,OP_DRP ),e(:,OP_1),g(:,OP_DRZ),h(:,OP_DZ)) &
             +int5(ri2_79,f(:,OP_DZZP),e(:,OP_1),g(:,OP_DZ ),h(:,OP_DZ)) &
             +int5(ri2_79,f(:,OP_DZP ),e(:,OP_1),g(:,OP_DZZ),h(:,OP_DZ))
     end if

  case(1)
     temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
     temp79b = f(:,OP_DZP)*g(:,OP_DZ ) + f(:,OP_DRP)*g(:,OP_DR ) &
          +    f(:,OP_DZ )*g(:,OP_DZP) + f(:,OP_DR )*g(:,OP_DRP)

     if(surface_int) then
        temp = int5(ri4_79,e(:,OP_1),temp79a,norm79(:,1),h(:,OP_DRP)) &
             + int5(ri4_79,e(:,OP_1),temp79a,norm79(:,2),h(:,OP_DZP)) &
             + int5(ri4_79,e(:,OP_1),temp79b,norm79(:,1),h(:,OP_DR )) &
             + int5(ri4_79,e(:,OP_1),temp79b,norm79(:,2),h(:,OP_DZ ))
     else
        temp79c = e(:,OP_1)*h(:,OP_GS) &
             + e(:,OP_DZ)*h(:,OP_DZ) + e(:,OP_DR)*h(:,OP_DR)
        temp79d = e(:,OP_1)*h(:,OP_GSP) &
             + e(:,OP_DZ)*h(:,OP_DZP) + e(:,OP_DR)*h(:,OP_DRP)
        
        temp = -int3(ri4_79,temp79a,temp79d) &
               -int3(ri4_79,temp79b,temp79c)
     end if
  end select
#endif

  v2chipsipsi = temp
  return
end function v2chipsipsi


! v2chipsib
! =========
vectype function v2chipsib(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp


  temp79a = h(:,OP_1 )*f(:,OP_GS) &
       +    h(:,OP_DZ)*f(:,OP_DZ) + h(:,OP_DR)*f(:,OP_DR)

  temp79b = e(:,OP_DR)*h(:,OP_DZ) - e(:,OP_DZ)*h(:,OP_DR)

  temp79c = e(:,OP_DZ)*g(:,OP_DR) - e(:,OP_DR)*g(:,OP_DZ)
  
  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int3(ri_79,temp79a,temp79c) &
             + int4(ri_79,temp79b,f(:,OP_DZ),g(:,OP_DZ)) &
             + int4(ri_79,temp79b,f(:,OP_DR),g(:,OP_DR)) 
#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp  &
             -int5(ri3_79,f(:,OP_DZPP),e(:,OP_1),g(:,OP_DR),h(:,OP_1)) &
             +int5(ri3_79,f(:,OP_DRPP),e(:,OP_1),g(:,OP_DZ),h(:,OP_1))
#endif
     end if
     
  case(1)

     if(surface_int) then
        temp79a = e(:,OP_1)* &
             (norm79(:,1)*f(:,OP_DR) + norm79(:,2)*f(:,OP_DZ))
        temp79b = e(:,OP_1)*h(:,OP_1)* &
             (norm79(:,1)*g(:,OP_DZ) - norm79(:,2)*g(:,OP_DR))
        temp = int4(ri3_79,temp79a,g(:,OP_DZ),h(:,OP_DR)) &
             - int4(ri3_79,temp79a,g(:,OP_DR),h(:,OP_DZ)) &
             + int3(ri3_79,temp79b,f(:,OP_GS))
        if(itor.eq.1) then
           temp = temp - 2.*int3(ri4_79,temp79b,f(:,OP_DR))
        endif
     else
        temp = int3(ri3_79,temp79a,temp79c) &
             + int4(ri3_79,temp79b,f(:,OP_DZ),g(:,OP_DZ)) &
             + int4(ri3_79,temp79b,f(:,OP_DR),g(:,OP_DR))

        if(itor.eq.1) then
           temp = temp - &
                2.*int4(ri4_79,temp79c,f(:,OP_DR),h(:,OP_1))
        endif
#if defined(USE3D) || defined(USECOMPLEX)
        temp79d = &
         2.*(f(:,OP_DZP)*g(:,OP_DR ) - f(:,OP_DRP)*g(:,OP_DZ ))*h(:,OP_DP ) &
         +  (f(:,OP_DZ )*g(:,OP_DRP) - f(:,OP_DR )*g(:,OP_DZP))*h(:,OP_DP ) &
         +  (f(:,OP_DZ )*g(:,OP_DR ) - f(:,OP_DR )*g(:,OP_DZ ))*h(:,OP_DPP) &
         +  (f(:,OP_DZPP)*g(:,OP_DR ) - f(:,OP_DRPP)*g(:,OP_DZ ))*h(:,OP_1) &
         +  (f(:,OP_DZP )*g(:,OP_DRP) - f(:,OP_DRP )*g(:,OP_DZP))*h(:,OP_1)
        temp = temp - int3(ri5_79,e(:,OP_1),temp79d)
#endif
     end if

  end select

  v2chipsib = temp
  return
end function v2chipsib


! v2chibb
! =======
vectype function v2chibb(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
  select case (ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp =  &
             + int5(ri2_79,f(:,OP_DRP),e(:,OP_1),g(:,OP_DR),h(:,OP_1)) &
             + int5(ri2_79,f(:,OP_DZP),e(:,OP_1),g(:,OP_DZ),h(:,OP_1)) &
             - int5(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_1),h(:,OP_DP))
     end if

  case(1)
     if(surface_int) then
        temp = 0.
     else
        temp = int5(ri4_79,e(:,OP_1),g(:,OP_1),h(:,OP_DR),f(:,OP_DRP)) &
             + int5(ri4_79,e(:,OP_1),g(:,OP_1),h(:,OP_DZ),f(:,OP_DZP)) &
             + int5(ri4_79,e(:,OP_1),g(:,OP_DP),h(:,OP_DR),f(:,OP_DR)) &
             + int5(ri4_79,e(:,OP_1),g(:,OP_DP),h(:,OP_DZ),f(:,OP_DZ))
     end if
  end select
#endif

  v2chibb = temp
  return
end function v2chibb



! v2vchin
! =======
vectype function v2vchin(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(inertia.eq.0 .or. surface_int) then
     v2vchin = 0.
     return
  end if


  select case(ivform)
  case(0)
     temp =-int4(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
          - int4(e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1))
  case(1)
     temp =-int4(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
          - int4(e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1))
     if(itor.eq.1) then
        temp = temp &
             - 2.*int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),h(:,OP_1))
     endif
  end select

  v2vchin = temp
  return
end function v2vchin


! v2chibsb1
! =========
vectype function v2chibsb1(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

  v2chibsb1 = 0.
  return
end function v2chibsb1

! v2psisb2
! ========
vectype function v2psisb2(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ)) &
          - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR))
  endif

  v2psisb2 = temp
  return
end function v2psisb2

! v2bsb1
! ======
vectype function v2bsb1(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR)) &
          - int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ))
  end if

  v2bsb1 = temp
  return
end function v2bsb1


! V2vs
! ====
vectype function v2vs(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(idens.eq.0 .or. nosig.eq.1 .or. surface_int) then
     v2vs = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) + denm*nt79(:,OP_LP)

  select case(ivform)
  case(0)
     temp = -int3(e(:,OP_1),f(:,OP_1),temp79a)
  case(1)
     temp = -int4(r2_79,e(:,OP_1),f(:,OP_1),temp79a)
  end select

  v2vs = temp
  return
end function v2vs


! V2psif1
! =======
vectype function v2psif1(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(surface_int) then
     v2psif1 = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  ! same for both ivforms
  temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZPP)) &
       - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DRPP))
#else
  temp = 0.
#endif
  v2psif1 = temp
  return
end function v2psif1

! V2psif2
! =======
vectype function v2psif2(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(surface_int) then
     v2psif2 = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  ! same for both ivforms
  temp = int4(ri_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_DZP)) &
       - int4(ri_79,e(:,OP_1),f(:,OP_DZP),g(:,OP_DRP))
#else
  temp = 0.
#endif
  v2psif2 = temp
  return
end function v2psif2



! V2bf
! ====
vectype function v2bf(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
  
  if(surface_int) then
     v2bf = 0.
     return
  endif
#if defined(USE3D) || defined(USECOMPLEX)
  ! same for both ivforms
  temp = &
       - int3(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZP)) &
       - int3(e(:,OP_1),f(:,OP_DR),g(:,OP_DRP))
#else
  temp = 0.
#endif
  v2bf = temp
  return
end function v2bf

! V2ff
! ====
vectype function v2ff(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
  
  if(surface_int) then
     v2ff = 0.
     return
  endif
#if defined(USE3D) || defined(USECOMPLEX)
  ! same for both ivforms
  temp = &
       - int3(e(:,OP_1),f(:,OP_DZPP),g(:,OP_DZP)) &
       - int3(e(:,OP_1),f(:,OP_DRPP),g(:,OP_DRP))
#else
  temp = 0.
#endif
  v2ff = temp
  return
end function v2ff




!==============================================================================
! V3 TERMS
!==============================================================================

! V3chin
! ======
vectype function v3chin(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = &
             - int3(e(:,OP_DR),f(:,OP_DR),g(:,OP_1)) &
             - int3(e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1))
     end if

  case(1)
     if(surface_int) then
!!$        temp = int5(ri4_79,e(:,OP_1),g(:,OP_1),norm79(:,1),f(:,OP_DR)) &
!!$             + int5(ri4_79,e(:,OP_1),g(:,OP_1),norm79(:,2),f(:,OP_DZ))
        temp = 0.
     else
        temp = - int4(ri4_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_1)) &
               - int4(ri4_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1))
     end if
  end select

  v3chin = temp
  return
end function v3chin



! V3chimu
! =======
vectype function v3chimu(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = 2.*int3(e(:,OP_LP),f(:,OP_LP),h(:,OP_1 )) &
             +    int3(e(:,OP_LP),f(:,OP_DZ),g(:,OP_DZ)) &
             +    int3(e(:,OP_LP),f(:,OP_DR),g(:,OP_DR)) &
             +    int3(e(:,OP_DZ),f(:,OP_LP),g(:,OP_DZ)) &
             +    int3(e(:,OP_DR),f(:,OP_LP),g(:,OP_DR)) &
             +    int3(e(:,OP_DZ),f(:,OP_DZ),g(:,OP_LP)) &
             +    int3(e(:,OP_DR),f(:,OP_DR),g(:,OP_LP))
     end if

  case(1)
     if(surface_int) then
        temp = 2.* &
             (int5(ri4_79,g(:,OP_1),e(:,OP_DR),norm79(:,1),f(:,OP_DZZ)) &
             -int5(ri4_79,g(:,OP_1),e(:,OP_DR),norm79(:,2),f(:,OP_DRZ)) &
             -int5(ri4_79,g(:,OP_1),e(:,OP_DZ),norm79(:,1),f(:,OP_DRZ)) &
             +int5(ri4_79,g(:,OP_1),e(:,OP_DZ),norm79(:,2),f(:,OP_DRR)) &
             -int5(ri4_79,h(:,OP_1),norm79(:,1),e(:,OP_DR),f(:,OP_GS)) &
             -int5(ri4_79,h(:,OP_1),norm79(:,2),e(:,OP_DZ),f(:,OP_GS)))
        if(itor.eq.1) then
           temp = temp + 2.* &
                (int5(ri5_79,g(:,OP_1),e(:,OP_DZ),norm79(:,1),f(:,OP_DZ)) &
                -int5(ri5_79,g(:,OP_1),e(:,OP_DZ),norm79(:,2),f(:,OP_DR)) &
                +int5(ri5_79,g(:,OP_1),e(:,OP_DR),norm79(:,1),f(:,OP_DR)) &
                +int5(ri5_79,g(:,OP_1),e(:,OP_DR),norm79(:,2),f(:,OP_DZ)) &
                +int5(ri5_79,g(:,OP_1),e(:,OP_1),norm79(:,1),f(:,OP_DZZ)) &
                -int5(ri5_79,g(:,OP_1),e(:,OP_1),norm79(:,2),f(:,OP_DRZ)) &
                +2.*int5(ri6_79,g(:,OP_1),e(:,OP_1),norm79(:,2),f(:,OP_DZ)))
        endif

#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp &
             + int5(ri6_79,e(:,OP_1),g(:,OP_1),norm79(:,1),f(:,OP_DRPP)) &
             + int5(ri6_79,e(:,OP_1),g(:,OP_1),norm79(:,2),f(:,OP_DZPP))
#endif
     else
        temp79a = e(:,OP_DRR)
        temp79b = f(:,OP_DRR)
        temp79c = e(:,OP_DRZ)
        temp79d = f(:,OP_DRZ)
        if(itor.eq.1) then
           temp79a = temp79a - 2.*ri_79*e(:,OP_DR)
           temp79b = temp79b - 2.*ri_79*f(:,OP_DR)
           temp79c = temp79c -    ri_79*e(:,OP_DZ)
           temp79d = temp79d -    ri_79*f(:,OP_DZ)
        endif
        temp = 2.* &
             (int4(ri4_79,e(:,OP_DZZ),f(:,OP_DZZ),g(:,OP_1)) &
             +int4(ri4_79,temp79a,temp79b,g(:,OP_1)) &
             +2.*int4(ri4_79,temp79c,temp79d,g(:,OP_1)) & 
             +int4(ri4_79,e(:,OP_GS),f(:,OP_GS),h(:,OP_1)) &
             -int4(ri4_79,e(:,OP_GS),f(:,OP_GS),g(:,OP_1)))
        if(itor.eq.1) then
           temp = temp &
                + 2.*int4(ri6_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_1))
        endif
#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp - &
             (int4(ri6_79,e(:,OP_DZ),f(:,OP_DZPP),g(:,OP_1)) &
             +int4(ri6_79,e(:,OP_DR),f(:,OP_DRPP),g(:,OP_1)))
#endif
     end if
  end select

  v3chimu = temp
  return
end function v3chimu


! V3umu
! =====
vectype function v3umu(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int4(ri_79,e(:,OP_LP),f(:,OP_DR),g(:,OP_DZ)) &
             - int4(ri_79,e(:,OP_LP),f(:,OP_DZ),g(:,OP_DR)) &
             + int4(ri_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DR)) &
             - int4(ri_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DZ)) &
             + int4(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_LP)) &
             - int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_LP))
     end if

  case(1)
     temp79c = f(:,OP_DZZ) - f(:,OP_DRR)
     if(itor.eq.1) temp79c = temp79c - ri_79*f(:,OP_DR)

     if(surface_int) then
        temp = int5(ri_79,g(:,OP_1),norm79(:,1),e(:,OP_DZ),temp79c) &
             + int5(ri_79,g(:,OP_1),norm79(:,2),e(:,OP_DR),temp79c) &
             + 2.* &
             (int5(ri_79,g(:,OP_1),norm79(:,1),e(:,OP_DR),f(:,OP_DRZ)) &
             -int5(ri_79,g(:,OP_1),norm79(:,2),e(:,OP_DZ),f(:,OP_DRZ))) &
             + int5(ri_79,g(:,OP_1),norm79(:,1),e(:,OP_DZ),f(:,OP_LP)) &
             - int5(ri_79,g(:,OP_1),norm79(:,2),e(:,OP_DR),f(:,OP_LP))
        if(itor.eq.1) then
           temp79a = h(:,OP_1) - g(:,OP_1)
           temp = temp &
                + 2.*int5(ri2_79,g(:,OP_1),norm79(:,1),e(:,OP_DR),f(:,OP_DZ)) &
                + 2.*int5(ri2_79,g(:,OP_1),norm79(:,2),e(:,OP_1),f(:,OP_LP)) &
                + 4.* &
                (int5(ri2_79,temp79a,norm79(:,1),e(:,OP_DR),f(:,OP_DZ)) &
                +int5(ri2_79,temp79a,norm79(:,2),e(:,OP_DZ),f(:,OP_DZ))) &
                - 4.* &
                (int5(ri2_79,e(:,OP_1),norm79(:,1),f(:,OP_DRZ),h(:,OP_1)) &
                +int5(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_DZZ),h(:,OP_1)))
        endif

#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp &
             + int5(ri3_79,e(:,OP_1),norm79(:,2),f(:,OP_DRPP),g(:,OP_1)) &
             - int5(ri3_79,e(:,OP_1),norm79(:,1),f(:,OP_DZPP),g(:,OP_1))
#endif        

     else
        temp79a = e(:,OP_DZZ) - e(:,OP_DRR)
        temp79b = e(:,OP_DRZ)
        if(itor.eq.1) then
           temp79a = temp79a + 2.*ri_79*e(:,OP_DR)
           temp79b = temp79b -    ri_79*e(:,OP_DZ)
        endif
        
        temp = 2.* &
             (int4(ri_79,f(:,OP_DRZ),temp79a,g(:,OP_1)) &
             -int4(ri_79,temp79b,temp79c,g(:,OP_1)))
        if(itor.eq.1) then
           temp = temp - 2.* &
                (int4(ri2_79,e(:,OP_DRR),f(:,OP_DZ),g(:,OP_1)) &
                -int4(ri3_79,e(:,OP_DR ),f(:,OP_DZ),g(:,OP_1)))
           
           temp79d = g(:,OP_1)-h(:,OP_1)
           temp = temp + 4.*int4(ri2_79,e(:,OP_GS),f(:,OP_DZ),temp79d)
        endif
        
#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp &
             + int4(ri3_79,e(:,OP_DR),f(:,OP_DZPP),g(:,OP_1)) &
             - int4(ri3_79,e(:,OP_DZ),f(:,OP_DRPP),g(:,OP_1))
#endif
     end if
  end select

  v3umu = temp
  return
end function v3umu


! V3vmu
! =====
vectype function v3vmu(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
  select case(ivform)
  case(0)
     ! not implemented
     temp = 0.

  case(1)
     if(surface_int) then
        temp79a = h(:,OP_1) - g(:,OP_1)
        temp = -2.* &
             (int5(ri2_79,norm79(:,1),e(:,OP_DR),f(:,OP_DP),temp79a) &
             +int5(ri2_79,norm79(:,2),e(:,OP_DZ),f(:,OP_DP),temp79a)) &
             + 2.* &
             (int5(ri2_79,e(:,OP_1),norm79(:,1),f(:,OP_DRP),temp79a) &
             +int5(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_DZP),temp79a)) &
             + int5(ri2_79,e(:,OP_1),g(:,OP_1),norm79(:,1),f(:,OP_DRP)) &
             + int5(ri2_79,e(:,OP_1),g(:,OP_1),norm79(:,2),f(:,OP_DZP))

        if(itor.eq.1) then
           temp = temp &
                - 2.*int5(ri3_79,e(:,OP_1),g(:,OP_1),norm79(:,1),f(:,OP_DP))
        endif
     else
        temp = &
             - int4(ri2_79,e(:,OP_DZ),f(:,OP_DZP),g(:,OP_1)) &
             - int4(ri2_79,e(:,OP_DR),f(:,OP_DRP),g(:,OP_1)) &
             + 2.*int4(ri2_79,e(:,OP_GS),f(:,OP_DP),h(:,OP_1)) &
             - 2.*int4(ri2_79,e(:,OP_GS),f(:,OP_DP),g(:,OP_1))
        if(itor.eq.1) then
           temp = temp &
                + 2.*int4(ri3_79,e(:,OP_DR),f(:,OP_DP),g(:,OP_1))
        endif
     end if
  end select
#endif
  v3vmu = temp
  return
end function v3vmu

! V3un
! ====
vectype function v3un(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ)) &
             - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR)) 
     end if
  case(1)
     if(surface_int) then
!!$        temp = int5(ri_79,e(:,OP_1),g(:,OP_1),norm79(:,2),f(:,OP_DR)) &
!!$             - int5(ri_79,e(:,OP_1),g(:,OP_1),norm79(:,1),f(:,OP_DZ))
        temp = 0.
     else
        temp = int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_1)) &
             - int4(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_1)) 
     end if
  end select

  v3un = temp
  return
end function v3un


! V3p
! ===
vectype function v3p(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = -int2(e(:,OP_1),f(:,OP_LP))
     end if
     
  case(1)
     if(surface_int) then
!!$        temp = &
!!$             - int4(ri2_79,e(:,OP_1),norm79(:,1),f(:,OP_DR)) &
!!$             - int4(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_DZ))
        temp = 0.
     else
        temp = int3(ri2_79,e(:,OP_DZ),f(:,OP_DZ)) &
             + int3(ri2_79,e(:,OP_DR),f(:,OP_DR))
     end if
  end select

  v3p = temp
  return
end function v3p



! V3up
! ====
vectype function v3up(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int4(ri_79,e(:,OP_LP),f(:,OP_DR),g(:,OP_DZ)) &
             - int4(ri_79,e(:,OP_LP),f(:,OP_DZ),g(:,OP_DR))
     end if
  case(1)
     if(surface_int) then
        temp79a = norm79(:,1)*e(:,OP_DR) + norm79(:,2)*e(:,OP_DZ)
        temp = int4(ri_79,temp79a,f(:,OP_DZ),g(:,OP_DR)) &
             - int4(ri_79,temp79a,f(:,OP_DR),g(:,OP_DZ)) &
             + int5(ri_79,e(:,OP_1),norm79(:,1),f(:,OP_DRR),g(:,OP_DZ )) &
             - int5(ri_79,e(:,OP_1),norm79(:,1),f(:,OP_DRZ),g(:,OP_DR )) &
             + int5(ri_79,e(:,OP_1),norm79(:,1),f(:,OP_DR ),g(:,OP_DRZ)) &
             - int5(ri_79,e(:,OP_1),norm79(:,1),f(:,OP_DZ ),g(:,OP_DRR)) &
             + int5(ri_79,e(:,OP_1),norm79(:,2),f(:,OP_DRZ),g(:,OP_DZ )) &
             - int5(ri_79,e(:,OP_1),norm79(:,2),f(:,OP_DZZ),g(:,OP_DR )) &
             + int5(ri_79,e(:,OP_1),norm79(:,2),f(:,OP_DR ),g(:,OP_DZZ)) &
             - int5(ri_79,e(:,OP_1),norm79(:,2),f(:,OP_DZ ),g(:,OP_DRZ))
        if(itor.eq.1) then
           temp79b = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
           temp = temp + 2.*gam* &
                (int4(ri2_79,temp79a,f(:,OP_DZ),g(:,OP_1)) &
                -int5(ri2_79,e(:,OP_1),norm79(:,1),f(:,OP_DRZ),g(:,OP_1 )) &
                -int5(ri2_79,e(:,OP_1),norm79(:,1),f(:,OP_DZ ),g(:,OP_DR)) &
                -int5(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_DZZ),g(:,OP_1 )) &
                -int5(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_DZ ),g(:,OP_DZ))) &
                -int4(ri2_79,e(:,OP_1),norm79(:,1),temp79b)
        endif
     else
        temp = int4(ri_79,e(:,OP_GS),f(:,OP_DR),g(:,OP_DZ)) &
             - int4(ri_79,e(:,OP_GS),f(:,OP_DZ),g(:,OP_DR))
        if(itor.eq.1) then
           temp = temp - 2.*gam* &
                int4(ri2_79,e(:,OP_GS),f(:,OP_DZ),g(:,OP_1))
        endif
     end if
  end select

  v3up = temp
  
  return
end function v3up


! V3vp
! ====
vectype function v3vp(e,f,g)

  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e, f, g

  vectype :: temp
  temp = 0.

#if defined(USE3D) || defined(USECOMPLEX)
  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = gam*int4(ri2_79,e(:,OP_LP),f(:,OP_DP),g(:,OP_1))
     end if

  case(1)
     if(surface_int) then
        temp79a = norm79(:,1)*e(:,OP_DR) + norm79(:,2)*e(:,OP_DZ)
        temp = &
             - int4(ri2_79,temp79a,f(:,OP_1 ),g(:,OP_DP)) &
             - gam*int4(ri2_79,temp79a,f(:,OP_DP),g(:,OP_1 )) &
             + int5(ri2_79,e(:,OP_1),norm79(:,1),f(:,OP_DR),g(:,OP_DP )) &
             + int5(ri2_79,e(:,OP_1),norm79(:,1),f(:,OP_1 ),g(:,OP_DRP)) &
             + int5(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_DZ),g(:,OP_DP )) &
             + int5(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_1 ),g(:,OP_DZP)) &
             + gam * &
             (int5(ri2_79,e(:,OP_1),norm79(:,1),f(:,OP_DRP),g(:,OP_1 )) &
             +int5(ri2_79,e(:,OP_1),norm79(:,1),f(:,OP_DP ),g(:,OP_DR)) &
             +int5(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_DZP),g(:,OP_1 )) &
             +int5(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_DP ),g(:,OP_DZ)))
     else
        temp = gam*int4(ri2_79,e(:,OP_GS),f(:,OP_DP),g(:,OP_1)) &
             + int4(ri2_79,e(:,OP_GS),f(:,OP_1),g(:,OP_DP))
     end if
  end select
#endif

  v3vp = temp
  return
end function v3vp



! V3chip
! ======
vectype function v3chip(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int3(e(:,OP_LP),f(:,OP_DZ),g(:,OP_DZ)) &
             + int3(e(:,OP_LP),f(:,OP_DR),g(:,OP_DR)) &
             + gam*int3(e(:,OP_LP),f(:,OP_LP),g(:,OP_1))
     end if

  case(1)
     if(surface_int) then
        temp79a = norm79(:,1)*e(:,OP_DR) + norm79(:,2)*e(:,OP_DZ)
        temp = &
             - int4(ri4_79,temp79a,f(:,OP_DZ),g(:,OP_DZ)) &
             - int4(ri4_79,temp79a,f(:,OP_DR),g(:,OP_DR)) &
             - gam*int4(ri4_79,temp79a,f(:,OP_GS),g(:,OP_1))  &
             + int5(ri4_79,e(:,OP_1),norm79(:,1),f(:,OP_DRZ),g(:,OP_DZ )) &
             + int5(ri4_79,e(:,OP_1),norm79(:,1),f(:,OP_DRR),g(:,OP_DR )) &
             + int5(ri4_79,e(:,OP_1),norm79(:,1),f(:,OP_DZ ),g(:,OP_DRZ)) &
             + int5(ri4_79,e(:,OP_1),norm79(:,1),f(:,OP_DR ),g(:,OP_DRR)) &
             + int5(ri4_79,e(:,OP_1),norm79(:,2),f(:,OP_DZZ),g(:,OP_DZ )) &
             + int5(ri4_79,e(:,OP_1),norm79(:,2),f(:,OP_DRZ),g(:,OP_DR )) &
             + int5(ri4_79,e(:,OP_1),norm79(:,2),f(:,OP_DZ ),g(:,OP_DZZ)) &
             + int5(ri4_79,e(:,OP_1),norm79(:,2),f(:,OP_DR ),g(:,OP_DRZ))
        if(itor.eq.1) then
           temp = temp - 2.* &
                (int5(ri5_79,e(:,OP_1),norm79(:,1),f(:,OP_DZ),g(:,OP_DZ)) &
                +int5(ri5_79,e(:,OP_1),norm79(:,1),f(:,OP_DR),g(:,OP_DR)))
        endif
        !! missing some terms
     else
        temp = int4(ri4_79,e(:,OP_GS),f(:,OP_DZ),g(:,OP_DZ)) &
             + int4(ri4_79,e(:,OP_GS),f(:,OP_DR),g(:,OP_DR)) &
             + gam*int4(ri4_79,e(:,OP_GS),f(:,OP_GS),g(:,OP_1))
     end if
  end select

  v3chip = temp

  return
end function v3chip



! V3psipsi
! ========
vectype function v3psipsi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_GS)) &
             + int4(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_GS))
     end if

  case(1)
     if(surface_int) then
!!$        temp = &
!!$             - int5(ri4_79,e(:,OP_1),norm79(:,1),f(:,OP_DR),g(:,OP_GS)) &
!!$             - int5(ri4_79,e(:,OP_1),norm79(:,2),f(:,OP_DZ),g(:,OP_GS))
        temp = 0.
     else
        temp = int4(ri4_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_GS)) &
             + int4(ri4_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_GS))
     end if

  end select

  v3psipsi = temp

  return
end function v3psipsi


! V3psib
! ======
vectype function v3psib(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = - &
             (int4(ri3_79,e(:,OP_DZ),f(:,OP_DRP),g(:,OP_1)) &
             -int4(ri3_79,e(:,OP_DR),f(:,OP_DZP),g(:,OP_1)))
     end if

  case(1)
     if(surface_int) then
!!$        temp = int5(ri5_79,e(:,OP_1),norm79(:,2),f(:,OP_DRP),g(:,OP_1)) &
!!$             - int5(ri5_79,e(:,OP_1),norm79(:,1),f(:,OP_DZP),g(:,OP_1))
        temp = 0.
     else
        temp = - &
             (int4(ri5_79,e(:,OP_DZ),f(:,OP_DRP),g(:,OP_1)) &
             -int4(ri5_79,e(:,OP_DR),f(:,OP_DZP),g(:,OP_1)))
     end if
  end select
#else
  temp = 0.
#endif
  v3psib = temp
  return
end function v3psib


! V3psif
! ======
vectype function v3psif(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp
#if defined(USE3D) || defined(USECOMPLEX)
  select case(ivform)
  case(0)
     temp = int4(ri_79,e(:,OP_DZ),g(:,OP_DRP),f(:,OP_GS)) &
          - int4(ri_79,e(:,OP_DR),g(:,OP_DZP),f(:,OP_GS))

  case(1)
     if(surface_int) then
!!$        temp = int5(ri3_79,e(:,OP_1),f(:,OP_GS),norm79(:,1),g(:,OP_DZP)) &
!!$             - int5(ri3_79,e(:,OP_1),f(:,OP_GS),norm79(:,2),g(:,OP_DRP))
        temp = 0.
     else
        temp = int4(ri3_79,e(:,OP_DZ),g(:,OP_DRP),f(:,OP_GS)) &
             - int4(ri3_79,e(:,OP_DR),g(:,OP_DZP),f(:,OP_GS))
     end if
  end select
#else
  temp = 0.
#endif

  v3psif = temp

  return
end function v3psif


! V3bb
! ====
vectype function v3bb(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int4(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_1)) &
             + int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1))
     end if

  case(1)
     if(surface_int) then
!!$        temp = &
!!$             - int5(ri4_79,e(:,OP_1),norm79(:,1),f(:,OP_DR),g(:,OP_1)) &
!!$             - int5(ri4_79,e(:,OP_1),norm79(:,2),f(:,OP_DZ),g(:,OP_1))
        temp = 0.
     else
        temp = int4(ri4_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_1)) &
             + int4(ri4_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1))
     end if
  end select

  v3bb = temp
  return
end function v3bb


! V3bf
! ====
vectype function v3bf(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = - &
             (int4(ri2_79,e(:,OP_1),g(:,OP_DZPP),f(:,OP_DZ)) &
             +int4(ri2_79,e(:,OP_1),g(:,OP_DRPP),f(:,OP_DR)))
     end if

  case(1)
     if(surface_int) then
!!$        temp = &
!!$             - int5(ri4_79,e(:,OP_1),norm79(:,1),g(:,OP_DRPP),f(:,OP_1)) &
!!$             - int5(ri4_79,e(:,OP_1),norm79(:,2),g(:,OP_DZPP),f(:,OP_1))
        temp = 0.
     else
        temp = int4(ri4_79,e(:,OP_DR),g(:,OP_DRPP),f(:,OP_1)) &
             + int4(ri4_79,e(:,OP_DZ),g(:,OP_DZPP),f(:,OP_1))
     end if
  end select
#else
  temp = 0.
#endif

  v3bf = temp
end function v3bf


! V3psisb1
! ========
vectype function v3psisb1(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DZ)) &
          + int4(ri2_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DR)) &
          + int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_GS)) &
          + int4(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_GS))
  end if

  v3psisb1 = temp

end function v3psisb1


! V3bsb2
! ======
vectype function v3bsb2(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1 )) &
          + int4(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_1 )) &
          + int4(ri2_79,e(:,OP_DZ),f(:,OP_1 ),g(:,OP_DZ)) &
          + int4(ri2_79,e(:,OP_DR),f(:,OP_1 ),g(:,OP_DR))
  end if

  v3bsb2 = temp
  return
end function v3bsb2


! V3upsipsi
! =========
vectype function v3upsipsi(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e, f, g, h

  vectype :: temp

  if(surface_int) then
     select case(ivform)
     case(0)
        temp = 0.
     case(1)
        temp79a = e(:,OP_DZ)*h(:,OP_DZ) + e(:,OP_DR)*h(:,OP_DR)
        temp = int5(ri3_79,temp79a,norm79(:,1),f(:,OP_DRZ),g(:,OP_DR )) &
             - int5(ri3_79,temp79a,norm79(:,1),f(:,OP_DRR),g(:,OP_DZ )) &
             + int5(ri3_79,temp79a,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRR)) &
             - int5(ri3_79,temp79a,norm79(:,1),f(:,OP_DR ),g(:,OP_DRZ)) &
             + int5(ri3_79,temp79a,norm79(:,2),f(:,OP_DZZ),g(:,OP_DR )) &
             - int5(ri3_79,temp79a,norm79(:,2),f(:,OP_DRZ),g(:,OP_DZ )) &
             + int5(ri3_79,temp79a,norm79(:,2),f(:,OP_DZ ),g(:,OP_DRZ)) &
             - int5(ri3_79,temp79a,norm79(:,2),f(:,OP_DR ),g(:,OP_DZZ))
        if(itor.eq.1) then
           temp = temp &
                + int5(ri4_79,temp79a,norm79(:,1),f(:,OP_DZ),g(:,OP_DR)) &
                - int5(ri4_79,temp79a,norm79(:,1),f(:,OP_DR),g(:,OP_DZ))
        endif
     end select
  else
  
  ! [f,g],r
  temp79b = f(:,OP_DRZ)*g(:,OP_DR ) - f(:,OP_DRR)*g(:,OP_DZ) &
       +    f(:,OP_DZ )*g(:,OP_DRR) - f(:,OP_DR )*g(:,OP_DRZ)
  
  ! [f,g],z
  temp79c = f(:,OP_DZZ)*g(:,OP_DR ) - f(:,OP_DRZ)*g(:,OP_DZ) &
       +    f(:,OP_DZ )*g(:,OP_DRZ) - f(:,OP_DR )*g(:,OP_DZZ)
  
  temp = int4(ri3_79,e(:,OP_DZ ),temp79c,h(:,OP_GS )) &
       + int4(ri3_79,e(:,OP_DR ),temp79b,h(:,OP_GS )) &
       - int4(ri3_79,e(:,OP_DZZ),temp79c,h(:,OP_DZ )) &
       - int4(ri3_79,e(:,OP_DR ),temp79c,h(:,OP_DRZ)) &
       - int4(ri3_79,e(:,OP_DZ ),temp79c,h(:,OP_DZZ)) &
       - int4(ri3_79,e(:,OP_DRZ),temp79c,h(:,OP_DR )) &
       - int4(ri3_79,e(:,OP_DRZ),temp79b,h(:,OP_DZ )) &
       - int4(ri3_79,e(:,OP_DR ),temp79b,h(:,OP_DRR)) &
       - int4(ri3_79,e(:,OP_DZ ),temp79b,h(:,OP_DRZ)) &
       - int4(ri3_79,e(:,OP_DRR),temp79b,h(:,OP_DR ))
    
  if(itor.eq.1) then
     ! [f,g]
     temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
  
     select case(ivform)
     case(0)
        temp = temp &
             - int4(ri4_79,e(:,OP_DR ),temp79a,h(:,OP_GS )) &
             + int4(ri4_79,e(:,OP_DRZ),temp79a,h(:,OP_DZ )) &
             + int4(ri4_79,e(:,OP_DZ ),temp79a,h(:,OP_DRZ)) &
             + int4(ri4_79,e(:,OP_DRR),temp79a,h(:,OP_DR )) &
             + int4(ri4_79,e(:,OP_DR ),temp79a,h(:,OP_DRR))

     case(1)

        temp79d = e(:,OP_DZ)*h(:,OP_DZ) + e(:,OP_DR)*h(:,OP_DR)
        
        temp = temp &
             +2.*int3(ri4_79,temp79b,temp79d) &
             -   int4(ri4_79,e(:,OP_DRZ),temp79a,h(:,OP_DZ )) &
             -   int4(ri4_79,e(:,OP_DZ ),temp79a,h(:,OP_DRZ)) &
             -   int4(ri4_79,e(:,OP_DRR),temp79a,h(:,OP_DR )) &
             -   int4(ri4_79,e(:,OP_DR ),temp79a,h(:,OP_DRR)) &
             +   int4(ri4_79,e(:,OP_DR),temp79a,h(:,OP_GS)) &
             +2.*int3(ri5_79,temp79d,temp79a)
     end select
  endif

  end if

  v3upsipsi = temp
  return
end function v3upsipsi


! V3upsib
! =======
vectype function v3upsib(e,f,g,h)

  use basic
  use arrays
  use m3dc1_nint

  implicit none
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

  vectype :: temp
  temp = 0.

#if defined(USE3D) || defined(USECOMPLEX) 
  select case (ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = temp    &
             +int5(ri4_79,f(:,OP_DZZP),e(:,OP_DR ),g(:,OP_DR),h(:,OP_1)) &
             -int5(ri4_79,f(:,OP_DRZP),e(:,OP_DZ ),g(:,OP_DR),h(:,OP_1)) &
             -int5(ri4_79,f(:,OP_DRP ),e(:,OP_DRR),g(:,OP_DR),h(:,OP_1)) &
             -int5(ri4_79,f(:,OP_DZP ),e(:,OP_DRZ),g(:,OP_DR),h(:,OP_1)) &
             +int5(ri4_79,f(:,OP_DRRP),e(:,OP_DZ ),g(:,OP_DZ),h(:,OP_1)) &
             -int5(ri4_79,f(:,OP_DRZP),e(:,OP_DR ),g(:,OP_DZ),h(:,OP_1)) &
             -int5(ri4_79,f(:,OP_DZP ),e(:,OP_DZZ),g(:,OP_DZ),h(:,OP_1)) &
             -int5(ri4_79,f(:,OP_DRP ),e(:,OP_DRZ),g(:,OP_DZ),h(:,OP_1))
        
        if(itor.eq.1) then
           temp = temp    &
                -int5(ri5_79,f(:,OP_DRP),e(:,OP_DR),g(:,OP_DR),h(:,OP_1)) &
                -int5(ri5_79,f(:,OP_DRP),e(:,OP_DZ),g(:,OP_DZ),h(:,OP_1))
        end if
     end if

  case(1)
     if(surface_int) then
        temp79a = h(:,OP_1)*(norm79(:,1)*e(:,OP_DZ) - norm79(:,2)*e(:,OP_DR))
        temp = int4(ri4_79,temp79a,f(:,OP_DRP),g(:,OP_DZ )) &
             - int4(ri4_79,temp79a,f(:,OP_DZP),g(:,OP_DR )) &
             + int4(ri4_79,temp79a,f(:,OP_DR ),g(:,OP_DZP)) &
             - int4(ri4_79,temp79a,f(:,OP_DZ ),g(:,OP_DRP))
     else
        temp79a = e(:,OP_DZ)*g(:,OP_DRP) - e(:,OP_DR)*g(:,OP_DZP)
        temp79b = h(:,OP_DZ)*f(:,OP_DR) - h(:,OP_DR)*f(:,OP_DZ)
        temp79c = f(:,OP_DZP)*g(:,OP_DR) - f(:,OP_DRP)*g(:,OP_DZ) &
             +    f(:,OP_DZ)*g(:,OP_DRP) - f(:,OP_DR)*g(:,OP_DZP)
        temp79d = e(:,OP_DZ)*h(:,OP_DR) - e(:,OP_DR)*h(:,OP_DZ)
        temp79e = h(:,OP_DP)*f(:,OP_GS) + h(:,OP_1)*f(:,OP_GSP) &
             + f(:,OP_DZP)*h(:,OP_DZ ) + f(:,OP_DRP)*h(:,OP_DR ) &
             + f(:,OP_DZ )*h(:,OP_DZP) + f(:,OP_DR )*h(:,OP_DRP)
        temp79f = h(:,OP_DP)*(e(:,OP_DZ)*f(:,OP_DZ )+e(:,OP_DR)*f(:,OP_DR )) &
             +    h(:,OP_1 )*(e(:,OP_DZ)*f(:,OP_DZP)+e(:,OP_DR)*f(:,OP_DRP))
        if(itor.eq.1) then 
           temp79d = temp79d - 4.*ri_79*e(:,OP_DZ)*h(:,OP_1)
        endif
        temp = int3(ri4_79,temp79a,temp79b) &
             + int3(ri4_79,temp79c,temp79d) &
             + int4(ri4_79,temp79e,e(:,OP_DZ),g(:,OP_DZ)) &
             + int4(ri4_79,temp79e,e(:,OP_DR),g(:,OP_DR)) &
             + int3(ri4_79,temp79f,g(:,OP_GS))
     end if
  end select
#endif
  v3upsib = temp
  return
end function v3upsib



! V3ubb
! =====
vectype function v3ubb(e,f,g,h)

  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int5(ri3_79,e(:,OP_GS),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1)) &
             - int5(ri3_79,e(:,OP_GS),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1))

        if(itor.eq.1) then
           temp = temp + 2.* &
                int5(ri4_79,e(:,OP_GS),f(:,OP_DZ),g(:,OP_1),h(:,OP_1))
        endif
#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp -   &
             int5(ri5_79,f(:,OP_DRPP),e(:,OP_DZ),g(:,OP_1),h(:,OP_1)) &
             +int5(ri5_79,f(:,OP_DZPP),e(:,OP_DR),g(:,OP_1),h(:,OP_1))
#endif
     end if

  case(1)
     if(surface_int) then
        temp79a = h(:,OP_1)*(norm79(:,1)*e(:,OP_DR) + norm79(:,2)*e(:,OP_DZ))
        if(itor.eq.1) then
           temp79a = temp79a &
                - 2.*ri_79*norm79(:,1)*e(:,OP_1)*h(:,OP_1)
        end if
        temp = int4(ri3_79,temp79a,f(:,OP_DZ),g(:,OP_DR)) &
             - int4(ri3_79,temp79a,f(:,OP_DR),g(:,OP_DZ))
     else
        temp79a = h(:,OP_1)*(g(:,OP_DZ)*f(:,OP_DR) - g(:,OP_DR)*f(:,OP_DZ))

        temp = int3(ri3_79,e(:,OP_GS),temp79a)

!  scj removed 4/1/2011
!        if(itor.eq.1) then
!           temp = temp - &
!                2.*int3(ri4_79,e(:,OP_DR),temp79a)
!        endif

#if defined(USE3D) || defined(USECOMPLEX)
        temp79b = (e(:,OP_DZ)*f(:,OP_DR)-e(:,OP_DR)*f(:,OP_DZ))*g(:,OP_DPP) &
             + 2.*(e(:,OP_DZ)*f(:,OP_DRP)-e(:,OP_DR)*f(:,OP_DZP))*g(:,OP_DP) &
             +    (e(:,OP_DZ)*f(:,OP_DRPP)-e(:,OP_DR)*f(:,OP_DZPP))*g(:,OP_1)
        temp = temp - int3(ri5_79,temp79b,h(:,OP_1))
#endif
     end if

  end select


  v3ubb = temp
  return
end function v3ubb


! v3vpsipsi
! =========
vectype function v3vpsipsi(e,f,g,h)
!
!  e trial
!  f lin
!  g psi
!  h psi
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = temp  &
             +int5(ri2_79,f(:,OP_DP),e(:,OP_DR ),g(:,OP_DRR),h(:,OP_DR)) &
             +int5(ri2_79,f(:,OP_DP),e(:,OP_DRR),g(:,OP_DR ),h(:,OP_DR)) &
             +int5(ri2_79,f(:,OP_DP),e(:,OP_DZ ),g(:,OP_DRZ),h(:,OP_DR)) &
             +int5(ri2_79,f(:,OP_DP),e(:,OP_DRZ),g(:,OP_DZ ),h(:,OP_DR)) &
             +int5(ri2_79,f(:,OP_DP),e(:,OP_DR ),g(:,OP_DRZ),h(:,OP_DZ)) &
             +int5(ri2_79,f(:,OP_DP),e(:,OP_DRZ),g(:,OP_DR ),h(:,OP_DZ)) &
             +int5(ri2_79,f(:,OP_DP),e(:,OP_DZ ),g(:,OP_DZZ),h(:,OP_DZ)) &
             +int5(ri2_79,f(:,OP_DP),e(:,OP_DZZ),g(:,OP_DZ ),h(:,OP_DZ)) &
             -int5(ri2_79,f(:,OP_DP),e(:,OP_DR ),g(:,OP_GS ),h(:,OP_DR)) &
             -int5(ri2_79,f(:,OP_DP),e(:,OP_DZ ),g(:,OP_GS ),h(:,OP_DZ))
     end if

  case(1)
     if(surface_int) then
        temp = 0.
     else
        temp79a = e(:,OP_DZ)*h(:,OP_DRP) - e(:,OP_DR)*h(:,OP_DZP)
        temp79b = g(:,OP_DZ)*f(:,OP_DR) - g(:,OP_DR)*f(:,OP_DZ)
        temp79e = f(:,OP_DP)*g(:,OP_GS) + f(:,OP_1)*g(:,OP_GSP) &
             + g(:,OP_DZP)*f(:,OP_DZ ) + g(:,OP_DRP)*f(:,OP_DR ) &
             + g(:,OP_DZ )*f(:,OP_DZP) + g(:,OP_DR )*f(:,OP_DRP)
        temp79f = f(:,OP_DP)*(e(:,OP_DZ)*g(:,OP_DZ )+e(:,OP_DR)*g(:,OP_DR )) &
             +    f(:,OP_1 )*(e(:,OP_DZ)*g(:,OP_DZP)+e(:,OP_DR)*g(:,OP_DRP))

        temp = int3(ri4_79,temp79a,temp79b) &
             - int4(ri4_79,temp79e,e(:,OP_DZ),h(:,OP_DZ)) &
             - int4(ri4_79,temp79e,e(:,OP_DR),h(:,OP_DR)) &
             - int3(ri4_79,temp79f,h(:,OP_GS))
     end if
  end select
#endif

  v3vpsipsi = temp
  return
end function v3vpsipsi



! v3vpsib
! =======
vectype function v3vpsib(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int5(ri3_79,e(:,OP_GS),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1)) &
             - int5(ri3_79,e(:,OP_GS),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1))

        if(itor.eq.1) then
           temp = temp - 2.* &
                int5(ri4_79,e(:,OP_GS),f(:,OP_1),g(:,OP_DZ),h(:,OP_1))
        endif
        
#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp    &
             +int5(ri5_79,f(:,OP_DPP),e(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
             -int5(ri5_79,f(:,OP_DPP),e(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) 
#endif
     end if
     
  case(1)
     if(surface_int) then
        temp79a = norm79(:,1)*e(:,OP_DR) + norm79(:,2)*e(:,OP_DZ)
        if(itor.eq.1) then
           temp79a = temp79a &
                - 2.*ri_79*norm79(:,1)*e(:,OP_1)
        end if
        temp = int5(ri3_79,temp79a,f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
             - int5(ri3_79,temp79a,f(:,OP_DR),g(:,OP_DZ),h(:,OP_1))
     else
        temp79a = h(:,OP_1)*(g(:,OP_DZ)*f(:,OP_DR) - g(:,OP_DR)*f(:,OP_DZ))

        temp = int3(ri3_79,e(:,OP_GS),temp79a)

!   scj removed 4/1/2011        
!        if(itor.eq.1) then
!           temp = temp - &
!                2.*int3(ri4_79,e(:,OP_DR),temp79a)
!        endif

#if defined(USE3D) || defined(USECOMPLEX)
        temp79b = f(:,OP_DPP)*(e(:,OP_DZ)*g(:,OP_DR)-e(:,OP_DR)*g(:,OP_DZ)) &
             + 2.*f(:,OP_DP)*(e(:,OP_DZ)*g(:,OP_DRP)-e(:,OP_DR)*g(:,OP_DZP)) &
             +    f(:,OP_1)*(e(:,OP_DZ)*g(:,OP_DRPP)-e(:,OP_DR)*g(:,OP_DZPP))
        temp = temp + int3(ri5_79,temp79b,h(:,OP_1))
#endif
     end if

  end select

  v3vpsib = temp
  return
end function v3vpsib


! V3vbb
! =====
vectype function v3vbb(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp


#if defined(USE3D) || defined(USECOMPLEX)

  select case(ivform)
  case(0)
     ! not implemented
     temp = 0.
  case(1)
     if(surface_int) then
        temp = 0.
     else
        temp = 0.
     end if
  end select
#else
  temp = 0.
#endif

  v3vbb = temp
  return
end function v3vbb


! V3chipsipsi
! ===========
vectype function v3chipsipsi(e,f,g,h)

  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e, f, g, h
  vectype :: temp

  ! <f,g>,r
  temp79b = f(:,OP_DRR)*g(:,OP_DR ) + f(:,OP_DRZ)*g(:,OP_DZ ) &
       +    f(:,OP_DR )*g(:,OP_DRR) + f(:,OP_DZ )*g(:,OP_DRZ)
  
  ! <f,g>,z
  temp79c = f(:,OP_DRZ)*g(:,OP_DR ) + f(:,OP_DZZ)*g(:,OP_DZ ) &
       +    f(:,OP_DR )*g(:,OP_DRZ) + f(:,OP_DZ )*g(:,OP_DZZ)

  ! <e,h>,r
  temp79e = e(:,OP_DRR)*h(:,OP_DR ) + e(:,OP_DRZ)*h(:,OP_DZ ) &
       +    e(:,OP_DR )*h(:,OP_DRR) + e(:,OP_DZ )*h(:,OP_DRZ)
  
  ! <e,h>,z
  temp79f = e(:,OP_DRZ)*h(:,OP_DR ) + e(:,OP_DZZ)*h(:,OP_DZ ) &
       +    e(:,OP_DR )*h(:,OP_DRZ) + e(:,OP_DZ )*h(:,OP_DZZ)
  
  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int3(ri2_79,temp79b,temp79e) &
             + int3(ri2_79,temp79c,temp79f) &
             - int4(ri2_79,e(:,OP_DZ),temp79c,h(:,OP_GS)) &
             - int4(ri2_79,e(:,OP_DR),temp79b,h(:,OP_GS))
     end if

  case(1)
     if(surface_int) then
        temp79a = e(:,OP_DZ)*h(:,OP_DZ) + e(:,OP_DR)*h(:,OP_DR)
        temp = &
             - int5(ri6_79,temp79a,norm79(:,1),f(:,OP_DRZ),g(:,OP_DZ )) &
             - int5(ri6_79,temp79a,norm79(:,1),f(:,OP_DRR),g(:,OP_DR )) &
             - int5(ri6_79,temp79a,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRZ)) &
             - int5(ri6_79,temp79a,norm79(:,1),f(:,OP_DR ),g(:,OP_DRR)) &
             - int5(ri6_79,temp79a,norm79(:,2),f(:,OP_DZZ),g(:,OP_DZ )) &
             - int5(ri6_79,temp79a,norm79(:,2),f(:,OP_DRZ),g(:,OP_DR )) &
             - int5(ri6_79,temp79a,norm79(:,2),f(:,OP_DZ ),g(:,OP_DZZ)) &
             - int5(ri6_79,temp79a,norm79(:,2),f(:,OP_DR ),g(:,OP_DRZ))
        if(itor.eq.1) then
           temp = temp + 2.* &
                (int5(ri7_79,temp79a,f(:,OP_DZ),g(:,OP_DZ),norm79(:,1)) &
                +int5(ri7_79,temp79a,f(:,OP_DR),g(:,OP_DR),norm79(:,1)))
        endif
     else
        temp = int3(ri6_79,temp79b,temp79e) &
             + int3(ri6_79,temp79c,temp79f) &
             - int4(ri6_79,e(:,OP_DZ),temp79c,h(:,OP_GS)) &
             - int4(ri6_79,e(:,OP_DR),temp79b,h(:,OP_GS))

        if(itor.eq.1) then
           temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
           temp79d = e(:,OP_DZ)*h(:,OP_DZ) + e(:,OP_DR)*h(:,OP_DR)
           
           temp = temp &
                -2.*int3(ri7_79,temp79b,temp79d) &
                -2.*int3(ri7_79,temp79a,temp79e) &
                +2.*int4(ri7_79,e(:,OP_DR),temp79a,h(:,OP_GS)) &
                +4.*int3(ri8_79,temp79a,temp79d)
        endif
     end if
  end select

  v3chipsipsi = temp
  return
end function v3chipsipsi


! V3chipsib
! =========
vectype function v3chipsib(e,f,g,h)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

  vectype :: temp
  temp = 0.

#if defined(USE3D) || defined(USECOMPLEX)
  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = temp    &
             +int5(ri3_79,f(:,OP_DRRP),e(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
             +int5(ri3_79,f(:,OP_DRP),e(:,OP_DRZ),g(:,OP_DR),h(:,OP_1)) &
             -int5(ri3_79,f(:,OP_DRZP),e(:,OP_DR),g(:,OP_DR),h(:,OP_1)) &
             -int5(ri3_79,f(:,OP_DZP),e(:,OP_DRR),g(:,OP_DR),h(:,OP_1)) &
             +int5(ri3_79,f(:,OP_DRZP),e(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
             +int5(ri3_79,f(:,OP_DRP),e(:,OP_DZZ),g(:,OP_DZ),h(:,OP_1)) &
             -int5(ri3_79,f(:,OP_DZZP),e(:,OP_DR),g(:,OP_DZ),h(:,OP_1)) &
             -int5(ri3_79,f(:,OP_DZP),e(:,OP_DRZ),g(:,OP_DZ),h(:,OP_1))
     
        if(itor.eq.1) then
           temp = temp    &
                -int5(ri4_79,f(:,OP_DZP),e(:,OP_DR),g(:,OP_DR),h(:,OP_1)) &
                +int5(ri4_79,f(:,OP_DRP),e(:,OP_DZ),g(:,OP_DR),h(:,OP_1))
        endif
     end if

  case(1)
     if(surface_int) then
        temp79a = h(:,OP_1)*(norm79(:,1)*e(:,OP_DZ) - norm79(:,2)*e(:,OP_DR))
        temp = int4(ri7_79,temp79a,f(:,OP_DZP),g(:,OP_DZ )) &
             + int4(ri7_79,temp79a,f(:,OP_DRP),g(:,OP_DR )) &
             + int4(ri7_79,temp79a,f(:,OP_DZ ),g(:,OP_DZP)) &
             + int4(ri7_79,temp79a,f(:,OP_DR ),g(:,OP_DRP))
     else
        temp79a = h(:,OP_1)*f(:,OP_GS) &
             + h(:,OP_DZ)*f(:,OP_DZ) + h(:,OP_DR)*f(:,OP_DR)
        temp79b = e(:,OP_DZ)*h(:,OP_DR) - e(:,OP_DR)*h(:,OP_DZ)
        temp79c = f(:,OP_DZP)*h(:,OP_DR ) - f(:,OP_DRP)*h(:,OP_DZ ) &
             +    f(:,OP_DZ )*h(:,OP_DRP) - f(:,OP_DR )*h(:,OP_DZP)
        temp79d = (f(:,OP_DZ )*e(:,OP_DR)-f(:,OP_DR )*e(:,OP_DZ))*h(:,OP_DP) &
             +    (f(:,OP_DZP)*e(:,OP_DR)-f(:,OP_DRP)*e(:,OP_DZ))*h(:,OP_1 )
        if(itor.eq.1) then
           temp79a = temp79a - 2.*ri_79*f(:,OP_DR)*h(:,OP_1)
           temp79b = temp79b - 4.*ri_79*e(:,OP_DZ)*h(:,OP_1)
           temp79c = temp79c - 4.*ri_79* &
                (f(:,OP_DZP)*h(:,OP_1) + f(:,OP_DZ)*h(:,OP_DP))
        endif
        temp = int4(ri7_79,e(:,OP_DZ),g(:,OP_DRP),temp79a) &
             - int4(ri7_79,e(:,OP_DR),g(:,OP_DZP),temp79a) &
             - int4(ri7_79,temp79b,f(:,OP_DZP),g(:,OP_DZ)) &
             - int4(ri7_79,temp79b,f(:,OP_DRP),g(:,OP_DR)) &
             - int4(ri7_79,temp79b,f(:,OP_DZ),g(:,OP_DZP)) &
             - int4(ri7_79,temp79b,f(:,OP_DR),g(:,OP_DRP)) &
             + int4(ri7_79,temp79c,e(:,OP_DZ),g(:,OP_DZ)) &
             + int4(ri7_79,temp79c,e(:,OP_DR),g(:,OP_DR)) &
             + int3(ri7_79,g(:,OP_GS),temp79d)
     end if
  end select
#endif
  v3chipsib = temp
  return
end function v3chipsib


! V3chibb
! =======
vectype function v3chibb(e,f,g,h)

  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e, f, g, h

  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int5(ri2_79,e(:,OP_GS),f(:,OP_GS),g(:,OP_1 ),h(:,OP_1)) &
             + int5(ri2_79,e(:,OP_GS),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
             + int5(ri2_79,e(:,OP_GS),f(:,OP_DR),g(:,OP_DR),h(:,OP_1))

#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp    &
             - int5(ri4_79,e(:,OP_DR),f(:,OP_DRPP),g(:,OP_1),h(:,OP_1)) &
             - int5(ri4_79,e(:,OP_DZ),f(:,OP_DZPP),g(:,OP_1),h(:,OP_1))
#endif
     end if
     
  case(1)
     if(surface_int) then
        temp79a = h(:,OP_1)*(norm79(:,1)*e(:,OP_DR) + norm79(:,2)*e(:,OP_DZ))
!        if(itor.eq.1) then
!           temp79a = temp79a &
!                - 2.*ri_79*norm79(:,1)*e(:,OP_1)*h(:,OP_1)
!        end if
        temp = &
             - int4(ri6_79,temp79a,f(:,OP_GS),g(:,OP_1)) &
             - int4(ri6_79,temp79a,f(:,OP_DZ),g(:,OP_DZ)) &
             - int4(ri6_79,temp79a,f(:,OP_DR),g(:,OP_DR))
        if(itor.eq.1) then
           temp = temp + 2.*int4(ri7_79,temp79a,f(:,OP_DR),g(:,OP_1))
        endif
     else
        temp79a = g(:,OP_1)*f(:,OP_GS) &
             + g(:,OP_DZ)*f(:,OP_DZ) + g(:,OP_DR)*f(:,OP_DR)
        
        if(itor.eq.1) then
           temp79a = temp79a - &
                2.*ri_79*f(:,OP_DR)*g(:,OP_1)
        end if
        
        temp = int4(ri6_79,e(:,OP_GS),temp79a,h(:,OP_1))

!    scj removed 4/2/2011        
!        if(itor.eq.1) then
!           temp = temp - &
!                2.*int4(ri7_79,e(:,OP_DR),temp79a,h(:,OP_1))
!        endif
        
#if defined(USE3D) || defined(USECOMPLEX)
        temp79b = &
             (e(:,OP_DZ)*f(:,OP_DZPP) + e(:,OP_DR)*f(:,OP_DRPP))*g(:,OP_1  ) &
        + 2.*(e(:,OP_DZ)*f(:,OP_DZP ) + e(:,OP_DR)*f(:,OP_DRP ))*g(:,OP_DP ) &
        +    (e(:,OP_DZ)*f(:,OP_DZ  ) + e(:,OP_DR)*f(:,OP_DR  ))*g(:,OP_DPP)
        temp = temp - int3(ri8_79,h(:,OP_1),temp79b)
#endif
     end if

  end select

  v3chibb = temp

  return
end function v3chibb


! V3uun
! =====
vectype function v3uun(e,f,g,h)

  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e, f, g, h

  vectype :: temp

  if(inertia.eq.0 .or. surface_int) then
     v3uun = 0.
     return
  end if

  select case(ivform)
  case(0)

     ! Revised 5/23/08 NMF
     temp79a = e(:,OP_DZ)*g(:,OP_DZ) + e(:,OP_DR)*g(:,OP_DR)
     temp79b = e(:,OP_DZ)*h(:,OP_DZ) + e(:,OP_DR)*h(:,OP_DR) &
          + e(:,OP_LP)*h(:,OP_1)
     temp79c = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)

     temp = - int4(ri2_79,f(:,OP_GS),temp79a,h(:,OP_1)) &
          -.5*int3(ri2_79,temp79b,temp79c)

  case(1)
     temp79a = e(:,OP_DZ)*(f(:,OP_DR)*g(:,OP_DRZ) - f(:,OP_DZ)*g(:,OP_DRR)) &
          +    e(:,OP_DR)*(f(:,OP_DZ)*g(:,OP_DRZ) - f(:,OP_DR)*g(:,OP_DZZ))

     temp = int2(temp79a,h(:,OP_1))

     if(itor.eq.1) then
        temp79a = g(:,OP_DZ)*(e(:,OP_DR)*f(:,OP_DZ) - e(:,OP_DZ)*f(:,OP_DR))
        temp = temp + int3(ri_79,temp79a,h(:,OP_1))
     endif
  end select
     

  v3uun = temp
  return
end function v3uun


! V3uvn
! =====
vectype function v3uvn(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0 .or. surface_int) then
     v3uvn = 0.
     return
  end if

  select case(ivform)
  case(0)
     ! not implemented
     temp = 0.
  case(1)
     temp = int5(ri_79,e(:,OP_DZ),f(:,OP_DRP),g(:,OP_1),h(:,OP_1)) &
          - int5(ri_79,e(:,OP_DR),f(:,OP_DZP),g(:,OP_1),h(:,OP_1))
  end select
#else
  temp = 0.
#endif
  v3uvn = temp
  return
end function v3uvn



! V3vvn
! =====
vectype function v3vvn(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e, f, g, h
  vectype :: temp

  if(inertia.eq.0 .or. surface_int) then
     v3vvn = 0.
     return
  end if

  if(itor.eq.0) then
     temp = 0.
  else
     select case(ivform)
     case(0)
        temp = -int5(ri3_79,e(:,OP_DR),f(:,OP_1),g(:,OP_1),h(:,OP_1))
     case(1)
        temp = -int5(ri_79,e(:,OP_DR),f(:,OP_1),g(:,OP_1),h(:,OP_1))
     end select
  endif

  v3vvn = temp
  return
end function v3vvn


! V3uchin
! =======
vectype function v3uchin(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e, f, g, h
  vectype :: temp

  if(inertia.eq.0 .or. surface_int) then
     v3uchin = 0.
     return
  end if


  select case(ivform)
  case(0)
     temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)

     temp = int4(ri_79,e(:,OP_LP),temp79a,h(:,OP_1 )) &
          + int5(ri_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DR),h(:,OP_1)) & 
          - int5(ri_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DZ),h(:,OP_1)) & 
          + int4(ri_79,e(:,OP_DZ),temp79a,h(:,OP_DZ)) & 
          + int4(ri_79,e(:,OP_DR),temp79a,h(:,OP_DR))

  case(1)
     temp79a = e(:,OP_DZ)*(f(:,OP_DR)*g(:,OP_DZZ) + f(:,OP_DRZ)*g(:,OP_DZ)  &
                          -f(:,OP_DZ)*g(:,OP_DRZ) + f(:,OP_DRR)*g(:,OP_DR)) &
              +e(:,OP_DR)*(f(:,OP_DZ)*g(:,OP_DRR) + f(:,OP_DRZ)*g(:,OP_DR)  &
                          -f(:,OP_DR)*g(:,OP_DRZ) + f(:,OP_DZZ)*g(:,OP_DZ))
     temp = int3(ri3_79,temp79a,h(:,OP_1))

     if(itor.eq.1) then
        temp79a = e(:,OP_DZ)*(f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)) &
             +    f(:,OP_DZ)*(e(:,OP_DZ)*g(:,OP_DZ) + e(:,OP_DR)*g(:,OP_DR))

        temp = temp + int3(ri4_79,temp79a,h(:,OP_1))
     end if
  end select
  
  v3uchin = temp
  return
end function v3uchin


! V3vchin
! =======
vectype function v3vchin(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(inertia.eq.0 .or. surface_int) then
     v3vchin = 0.
     return
  end if

  select case(ivform)
  case(0)
     ! not implemented
     temp = 0.
  case(1)
     temp = int5(ri4_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DZP),h(:,OP_1)) &
          + int5(ri4_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DRP),h(:,OP_1))
  end select
#else
  temp = 0.
#endif
  v3vchin = temp
  return
end function v3vchin



! V3chichin
! =========
vectype function v3chichin(e,f,g,h)

  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e, f, g, h
  vectype :: temp

  if(inertia.eq.0 .or. surface_int) then
     v3chichin = 0.
     return
  end if

  select case(ivform)
  case(0)
     ! Revised 5/23/08 nmf
     temp79a = e(:,OP_DZ)*h(:,OP_DZ) + e(:,OP_DR)*h(:,OP_DR) &
          + e(:,OP_LP)*h(:,OP_1)
     temp79b = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
     temp = -0.5*int2(temp79a,temp79b)

  case(1)
     temp79a = e(:,OP_DZ)*(f(:,OP_DZ)*g(:,OP_DZZ) + f(:,OP_DR)*g(:,OP_DRZ)) &
          +    e(:,OP_DR)*(f(:,OP_DZ)*g(:,OP_DRZ) + f(:,OP_DR)*g(:,OP_DRR))

     temp = int3(ri6_79,temp79a,h(:,OP_1))

     if(itor.eq.1) then
        temp79a = g(:,OP_DR)*(e(:,OP_DZ)*f(:,OP_DZ) + e(:,OP_DR)*f(:,OP_DR))

        temp = temp - 2.*int3(ri7_79,temp79a,h(:,OP_1))
     endif

  end select

  v3chichin = temp
  return
end function v3chichin


! V3ngrav
! =======
vectype function v3ngrav(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f
  vectype :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v3ngrav = 0.
     return
  endif

  if(surface_int) then
     temp = 0.
  else
     temp = gravz*int2(       e(:,OP_DZ),f(:,OP_1)) & 
          + gravr*int3(ri2_79,e(:,OP_DR),f(:,OP_1)) 
  end if

  v3ngrav = temp
  return
end function v3ngrav


! V3ungrav
! ========
vectype function v3ungrav(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v3ungrav = 0.
     return
  endif

  if(surface_int) then
     temp = 0.
  else
     temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
     
     temp = gravz*int3( ri_79,e(:,OP_DZ),temp79a) &
          + gravr*int3(ri3_79,e(:,OP_DR),temp79a)
  end if

  v3ungrav = temp
  return
   end function v3ungrav


! V3chingrav
! ==========
vectype function v3chingrav(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v3chingrav = 0.
     return
  endif

  if(surface_int) then
     temp = 0.
  else
     temp79a = -(f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR) &
          + f(:,OP_LP)*g(:,OP_1))

     temp = gravz*int2(       e(:,OP_DZ),temp79a) &
          + gravr*int3(ri2_79,e(:,OP_DR),temp79a)
  end if

  v3chingrav = temp
  return
end function v3chingrav


! V3ndenmgrav
! ===========
vectype function v3ndenmgrav(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f
  real, intent(in) :: g
  vectype :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     v3ndenmgrav = 0.
     return
  endif

  if(surface_int) then
     temp = 0.
  else
     temp = gravz*int2(       e(:,OP_DZ),f(:,OP_LP)) &
          + gravr*int3(ri2_79,e(:,OP_DR),f(:,OP_LP))
  end if

  v3ndenmgrav = g*temp
  return
end function v3ndenmgrav


! V3us
! ====
vectype function v3us(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(idens.eq.0 .or. nosig.eq.1) then
     v3us = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) + denm*nt79(:,OP_LP)

  if(surface_int) then
     temp = 0.
  else
     temp = int4(ri_79,e(:,OP_DZ),f(:,OP_DR),temp79a) &
          -int4(ri_79,e(:,OP_DR),f(:,OP_DZ),temp79a)
  end if

  v3us = temp
  return
end function v3us


! V3chis
! ======
vectype function v3chis(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(idens.eq.0 .or. nosig.eq.1) then
     v3chis = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = g(:,OP_1) + denm*nt79(:,OP_LP)

  select case(ivform)
  case(0) 
     if(surface_int) then
        temp = 0.
     else
        temp = int3(e(:,OP_DZ),f(:,OP_DZ),temp79a) &
             + int3(e(:,OP_DR),f(:,OP_DR),temp79a)
     end if

  case(1)
     if(surface_int) then
        temp = 0.
     else
        temp = int4(ri4_79,e(:,OP_DZ),f(:,OP_DZ),temp79a) &
             + int4(ri4_79,e(:,OP_DR),f(:,OP_DR),temp79a)
     end if
  end select

  v3chis = temp
  return
end function v3chis


!==============================================================================
! B1 TERMS
!==============================================================================


! B1psi
! =====
vectype function b1psi(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f
  vectype :: temp

!!$  if(jadv.eq.0) then
!!$     temp = int2(e(:,OP_1),f(:,OP_1))
!!$  else
!!$!     temp79a = e(:,OP_DZ)*f(:,OP_DZ) + e(:,OP_DR)*f(:,OP_DR)
!!$!     temp = -int2(ri2_79,temp79a) 
!!$!...changed 7/22/08    scj
!!$     temp = -int3(ri2_79,e(:,OP_DZ),f(:,OP_DZ)) &
!!$            -int3(ri2_79,e(:,OP_DR),f(:,OP_DR))
!!$  endif

  if(surface_int) then
     if(jadv.eq.0) then
        temp = 0.
     else
        temp = int4(ri2_79,e(:,OP_1),norm79(:,1),f(:,OP_DR)) &
             + int4(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_DZ)) &
             - int4(ri2_79,f(:,OP_1),norm79(:,1),e(:,OP_DR)) &
             - int4(ri2_79,f(:,OP_1),norm79(:,2),e(:,OP_DZ))
     end if
  else
     if(jadv.eq.0) then
        temp79a = e(:,OP_1)
     else
        temp79a = ri2_79*e(:,OP_GS)
     endif

     temp = int2(temp79a,f(:,OP_1))
  end if

  b1psi = temp
  return
end function b1psi


! B1psiu
! ======
vectype function b1psiu(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        if(jadv.eq.0) then
           temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ)) &
                - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR))
        else
           temp = int4(ri3_79,e(:,OP_DR),f(:,OP_DRZ),g(:,OP_DR )) &
                + int4(ri3_79,e(:,OP_DR),f(:,OP_DZ ),g(:,OP_DRR)) &
                - int4(ri3_79,e(:,OP_DR),f(:,OP_DRR),g(:,OP_DZ )) &
                - int4(ri3_79,e(:,OP_DR),f(:,OP_DR ),g(:,OP_DRZ)) &
                + int4(ri3_79,e(:,OP_DZ),f(:,OP_DZZ),g(:,OP_DR )) &
                + int4(ri3_79,e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_DRZ)) &
                - int4(ri3_79,e(:,OP_DZ),f(:,OP_DRZ),g(:,OP_DZ )) &
                - int4(ri3_79,e(:,OP_DZ),f(:,OP_DR ),g(:,OP_DZZ)) 
           if(itor.eq.1) then
              temp = temp &
                   - int4(ri4_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DR)) &
                   + int4(ri4_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DZ)) 
           endif
        endif
     end if

  case(1)

     ! surface terms
     if(surface_int) then
        if(jadv.eq.0) then
           temp = 0.
        else
           temp79a = f(:,OP_DR)*g(:,OP_DZ) - f(:,OP_DZ)*g(:,OP_DR)
           temp79b = f(:,OP_DRR)*g(:,OP_DZ ) - f(:,OP_DRZ)*g(:,OP_DR ) &
                +    f(:,OP_DR )*g(:,OP_DRZ) - f(:,OP_DZ )*g(:,OP_DRR)
           temp79c = f(:,OP_DRZ)*g(:,OP_DZ ) - f(:,OP_DZZ)*g(:,OP_DR ) &
                +    f(:,OP_DR )*g(:,OP_DZZ) - f(:,OP_DZ )*g(:,OP_DRZ)

           temp = int4(ri_79,e(:,OP_1 ),norm79(:,1),temp79b) &
                + int4(ri_79,e(:,OP_1 ),norm79(:,2),temp79c) &
                - int4(ri_79,e(:,OP_DR),norm79(:,1),temp79a) &
                - int4(ri_79,e(:,OP_DZ),norm79(:,2),temp79a)
           if(itor.eq.1) then
              temp = temp + 2.*int4(ri2_79,e(:,OP_1),norm79(:,1),temp79a)
           end if
        endif

     ! volume terms
     else
        if(jadv.eq.0) then
           temp79a = e(:,OP_1)
        else
           temp79a = ri2_79*e(:,OP_GS)
        endif

        temp = int4(r_79,temp79a,f(:,OP_DR),g(:,OP_DZ)) &
             - int4(r_79,temp79a,f(:,OP_DZ),g(:,OP_DR))
     end if
  end select
  
  b1psiu = temp
  return
end function b1psiu


! B1psiv
! ======
vectype function b1psiv(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
    temp = 0.
  else
     select case(ivform)
     case(0)
        if(surface_int) then
           temp = 0.
        else
           temp = int4(ri4_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DP)) &
                + int4(ri4_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DP)) &
                + int4(ri4_79,e(:,OP_DR),f(:,OP_DRP),g(:,OP_1)) &
                + int4(ri4_79,e(:,OP_DZ),f(:,OP_DZP),g(:,OP_1))
        end if
     case(1)
        if(surface_int) then
           if(inoslip_tor.eq.1) then
              temp = 0.
           else
           temp = &
                - int5(ri2_79,e(:,OP_1),norm79(:,1),f(:,OP_DRP),g(:,OP_1 )) &
                - int5(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_DZP),g(:,OP_1 )) &
                - int5(ri2_79,e(:,OP_1),norm79(:,1),f(:,OP_DR ),g(:,OP_DP)) &
                - int5(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_DZ ),g(:,OP_DP))
           endif
        else
           temp = int4(ri2_79,e(:,OP_DR),f(:,OP_DRP),g(:,OP_1 )) &
                + int4(ri2_79,e(:,OP_DZ),f(:,OP_DZP),g(:,OP_1 )) &
                + int4(ri2_79,e(:,OP_DR),f(:,OP_DR ),g(:,OP_DP)) &
                + int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_DP))
        endif
     end select
  endif
#else
  temp = 0
#endif

  b1psiv = temp
  return
end function b1psiv


! B1psid
! ======
vectype function b1psid(e,f,g)
  use basic
  use m3dc1_nint

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(mass_ratio.eq.0. .or. dbf.eq.0.) then
     b1psid = 0.
     return
  endif

  if(surface_int) then
     temp = 0.
  else
     temp = int3(e(:,OP_1),f(:,OP_GS),g(:,OP_1))
  endif

  b1psid = temp*me_mi*mass_ratio*dbf**2
  return
end function b1psid


! B1bu
! ====
vectype function b1bu(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     temp = 0.
  else
     select case(ivform)
     case(0)
        if(surface_int) then
           temp = 0.
        else
           temp = -(int4(ri4_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DZP)) &
                +   int4(ri4_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DRP)) &
                +   int4(ri4_79,e(:,OP_DZ),f(:,OP_DP),g(:,OP_DZ)) &
                +   int4(ri4_79,e(:,OP_DR),f(:,OP_DP),g(:,OP_DR))) 
        end if

     case(1)
        if(surface_int) then
           if(inoslip_pol.eq.1) then
              temp = 0.
           else
           temp = int5(ri2_79,e(:,OP_1),norm79(:,1),g(:,OP_DRP),f(:,OP_1 )) &
                + int5(ri2_79,e(:,OP_1),norm79(:,2),g(:,OP_DZP),f(:,OP_1 )) &
                + int5(ri2_79,e(:,OP_1),norm79(:,1),g(:,OP_DR ),f(:,OP_DP)) &
                + int5(ri2_79,e(:,OP_1),norm79(:,2),g(:,OP_DZ ),f(:,OP_DP))
           endif
        else
           temp = -(int4(ri2_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DZP)) &
                  + int4(ri2_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DRP)) &
                  + int4(ri2_79,e(:,OP_DZ),f(:,OP_DP),g(:,OP_DZ)) &
                  + int4(ri2_79,e(:,OP_DR),f(:,OP_DP),g(:,OP_DR)))
        end if
     end select
  endif

#else
  temp  = 0.
#endif
  b1bu = temp
end function b1bu


! B1bv
! ====
vectype function b1bv(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  b1bv = 0.
end function b1bv


! B1psichi
! ========
vectype function b1psichi(e,f,g)

  use basic
  use m3dc1_nint
  
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        if(jadv.eq.0) then
           temp = -int3(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ)) &
                -  int3(e(:,OP_1),f(:,OP_DR),g(:,OP_DR))
        else
           temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_DZZ),g(:,OP_DZ )) &
                + int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_DZZ)) &
                + int4(ri2_79,e(:,OP_DZ),f(:,OP_DRZ),g(:,OP_DR )) &
                + int4(ri2_79,e(:,OP_DZ),f(:,OP_DR ),g(:,OP_DRZ)) &
                + int4(ri2_79,e(:,OP_DR),f(:,OP_DRZ),g(:,OP_DZ )) &
                + int4(ri2_79,e(:,OP_DR),f(:,OP_DZ ),g(:,OP_DRZ)) &
                + int4(ri2_79,e(:,OP_DR),f(:,OP_DRR),g(:,OP_DR )) &
                + int4(ri2_79,e(:,OP_DR),f(:,OP_DR ),g(:,OP_DRR))
        endif
     end if

  case(1)
     if(jadv.eq.0) then
        if(surface_int) then
           temp = 0.
        else
           temp = -int4(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ)) &
                  -int4(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR))
        end if
     else
        if(surface_int) then
!!$           temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
!!$           temp = int4(ri4_79,temp79a,norm79(:,1),e(:,OP_DR)) &
!!$                + int4(ri4_79,temp79a,norm79(:,2),e(:,OP_DZ)) &
!!$                - int5(ri4_79,norm79(:,1),f(:,OP_DRZ),g(:,OP_DZ ),e(:,OP_1)) &
!!$                - int5(ri4_79,norm79(:,1),f(:,OP_DRR),g(:,OP_DR ),e(:,OP_1)) &
!!$                - int5(ri4_79,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRZ),e(:,OP_1)) &
!!$                - int5(ri4_79,norm79(:,1),f(:,OP_DR ),g(:,OP_DRR),e(:,OP_1)) &
!!$                - int5(ri4_79,norm79(:,2),f(:,OP_DZZ),g(:,OP_DZ ),e(:,OP_1)) &
!!$                - int5(ri4_79,norm79(:,2),f(:,OP_DRZ),g(:,OP_DR ),e(:,OP_1)) &
!!$                - int5(ri4_79,norm79(:,2),f(:,OP_DZ ),g(:,OP_DZZ),e(:,OP_1)) &
!!$                - int5(ri4_79,norm79(:,2),f(:,OP_DR ),g(:,OP_DRZ),e(:,OP_1))
!!$           if(itor.eq.1) then
!!$              temp = temp + 2.*int4(ri5_79,temp79a,norm79(:,1),e(:,OP_1))
!!$           endif
           temp = 0.
        else
           temp = int4(ri4_79,e(:,OP_DZ),f(:,OP_DZZ),g(:,OP_DZ )) &
                + int4(ri4_79,e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_DZZ)) &
                + int4(ri4_79,e(:,OP_DZ),f(:,OP_DRZ),g(:,OP_DR )) &
                + int4(ri4_79,e(:,OP_DZ),f(:,OP_DR ),g(:,OP_DRZ)) &
                + int4(ri4_79,e(:,OP_DR),f(:,OP_DRZ),g(:,OP_DZ )) &
                + int4(ri4_79,e(:,OP_DR),f(:,OP_DZ ),g(:,OP_DRZ)) &
                + int4(ri4_79,e(:,OP_DR),f(:,OP_DRR),g(:,OP_DR )) &
                + int4(ri4_79,e(:,OP_DR),f(:,OP_DR ),g(:,OP_DRR))
           if(itor.eq.1) then
              temp = temp - 2.* &
                   (int4(ri5_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DR)) &
                   +int4(ri5_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DZ)))
           endif
        end if
     endif

  end select

  b1psichi = temp
  return
end function b1psichi


! B1bchi
! ======
vectype function b1bchi(e,f,g)

  use basic
  use m3dc1_nint
  
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     temp = 0.
  else
     select case (ivform)
     case(0)
        if(surface_int) then
           temp = 0.
        else
           temp = int4(ri3_79,e(:,OP_DZ),g(:,OP_DRP),f(:,OP_1))  &
                - int4(ri3_79,e(:,OP_DR),g(:,OP_DZP),f(:,OP_1))  &
                + int4(ri3_79,e(:,OP_DZ),g(:,OP_DR),f(:,OP_DP))  &
                - int4(ri3_79,e(:,OP_DR),g(:,OP_DZ),f(:,OP_DP))
        end if

     case(1)
        if(surface_int) then
           if(inoslip_pol.eq.1) then
              temp = 0.
           else
           temp = int5(ri5_79,e(:,OP_1),f(:,OP_1 ),norm79(:,1),g(:,OP_DZP)) &
                - int5(ri5_79,e(:,OP_1),f(:,OP_1 ),norm79(:,2),g(:,OP_DRP)) &
                + int5(ri5_79,e(:,OP_1),f(:,OP_DP),norm79(:,1),g(:,OP_DZ )) &
                - int5(ri5_79,e(:,OP_1),f(:,OP_DP),norm79(:,2),g(:,OP_DR ))
           endif
        else
           temp = int4(ri5_79,e(:,OP_DZ),g(:,OP_DRP),f(:,OP_1))  &
                - int4(ri5_79,e(:,OP_DR),g(:,OP_DZP),f(:,OP_1))  &
                + int4(ri5_79,e(:,OP_DZ),g(:,OP_DR),f(:,OP_DP))  &
                - int4(ri5_79,e(:,OP_DR),g(:,OP_DZ),f(:,OP_DP))
        end if
     end select
  endif

#else
  temp = 0.
#endif
  b1bchi = temp
  return
end function b1bchi


! B1psieta
! ========
vectype function b1psieta(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(jadv.eq.0) then

     if(surface_int) then
        temp = 0.
     else
        temp = int3(e(:,OP_1),f(:,OP_GS),g(:,OP_1))

        if(hypf.ne.0.) then
           if(ihypeta.eq.1) then
              temp79a = e(:,OP_1)*g(:,OP_LP) + e(:,OP_LP)*g(:,OP_1) &
                   + 2.*(e(:,OP_DZ)*g(:,OP_DZ) + e(:,OP_DR)*g(:,OP_DR))
              temp = temp - int3(temp79a,f(:,OP_GS),h(:,OP_1))
              if(itor.eq.1) temp = temp + &
                   2.*int4(ri_79,temp79a,f(:,OP_DR),h(:,OP_1))
           else
              temp = temp - int3(e(:,OP_LP),f(:,OP_GS),h(:,OP_1))

              if(itor.eq.1) then
                 temp = temp - 2.*int4(ri_79,e(:,OP_DR),f(:,OP_GS),h(:,OP_1))
              endif

#if defined(USE3D) || defined(USECOMPLEX)
              temp = temp &
                   - 2.*int4(ri2_79,e(:,OP_1),f(:,OP_GSPP),h(:,OP_1)) &
                   - int4(ri2_79,e(:,OP_DZ),f(:,OP_DZPP),h(:,OP_1)) &
                   - int4(ri2_79,e(:,OP_DR),f(:,OP_DRPP),h(:,OP_1))
              if(itor.eq.1) then
                 temp = temp + 2.*int4(ri3_79,e(:,OP_1),f(:,OP_DRPP),h(:,OP_1))
              endif
#endif
           end if
        end if
     end if
  else
     if(surface_int) then
        if(inocurrent_norm.eq.1) then
           temp = 0.
        else
           temp = int5(ri2_79,e(:,OP_1),f(:,OP_GS),norm79(:,1),g(:,OP_DR)) &
                + int5(ri2_79,e(:,OP_1),f(:,OP_GS),norm79(:,2),g(:,OP_DZ)) &
                - int5(ri2_79,g(:,OP_1),f(:,OP_GS),norm79(:,1),e(:,OP_DR)) &
                - int5(ri2_79,g(:,OP_1),f(:,OP_GS),norm79(:,2),e(:,OP_DZ))
        endif

#if defined(USE3D) || defined(USECOMPLEX)
        if(inocurrent_norm.eq.1) then
           temp = temp
        else
           temp = temp &
                + int5(ri4_79,e(:,OP_1),norm79(:,1),f(:,OP_DRPP),g(:,OP_1)) &
                + int5(ri4_79,e(:,OP_1),norm79(:,2),f(:,OP_DZPP),g(:,OP_1))
        endif
#endif
     else
        temp = int4(ri2_79,g(:,OP_1),e(:,OP_GS),f(:,OP_GS))

#if defined(USE3D) || defined(USECOMPLEX)
        temp = temp - &
             (int4(ri4_79,e(:,OP_DZ),f(:,OP_DZPP),g(:,OP_1)) &
             +int4(ri4_79,e(:,OP_DR),f(:,OP_DRPP),g(:,OP_1)))
#endif
     end if
  endif

  b1psieta = temp
  return
end function b1psieta


! B1beta
! ======
vectype function b1beta(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     temp = 0.
  else 
     if(surface_int) then
        if(inocurrent_norm.eq.1) then
           temp = 0.
        else
           temp = int5(ri3_79,e(:,OP_1),g(:,OP_1),norm79(:,2),f(:,OP_DRP)) &
                - int5(ri3_79,e(:,OP_1),g(:,OP_1),norm79(:,1),f(:,OP_DZP))
        endif
     else
        temp = int4(ri3_79,e(:,OP_DR),f(:,OP_DZP),g(:,OP_1)) &
             - int4(ri3_79,e(:,OP_DZ),f(:,OP_DRP),g(:,OP_1))

        if(ihypeta.eq.0) then
           temp = temp - 2.*int4(ri2_79,e(:,OP_1),f(:,OP_DZP),h(:,OP_1))
        endif
     end if
  endif
#else
  temp = 0.
#endif

  b1beta = temp
  return
end function b1beta



! B1psipsid
! =========
vectype function b1psipsid(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itwofluid.eq.0) then
     b1psipsid = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     if(surface_int) then
        temp = 0.
     else
        temp = int5(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZP),h(:,OP_1)) &
             + int5(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DRP),h(:,OP_1))
     endif
  else
     if(surface_int) then
        temp79a = f(:,OP_DZ)*g(:,OP_DZP) + f(:,OP_DR)*g(:,OP_DRP)
        temp79b = e(:,OP_1)*h(:,OP_1)
        temp79c = e(:,OP_1)*h(:,OP_DP)
        temp = int5(ri4_79,e(:,OP_1),temp79a,norm79(:,1),h(:,OP_DR)) &
             + int5(ri4_79,e(:,OP_1),temp79a,norm79(:,2),h(:,OP_DZ)) &
             + int5(ri4_79,temp79b,norm79(:,1),f(:,OP_DRZ),g(:,OP_DZP )) &
             + int5(ri4_79,temp79b,norm79(:,1),f(:,OP_DRR),g(:,OP_DRP )) &
             + int5(ri4_79,temp79b,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRZP)) &
             + int5(ri4_79,temp79b,norm79(:,1),f(:,OP_DR ),g(:,OP_DRRP)) &
             + int5(ri4_79,temp79b,norm79(:,2),f(:,OP_DZZ),g(:,OP_DZP )) &
             + int5(ri4_79,temp79b,norm79(:,2),f(:,OP_DRZ),g(:,OP_DRP )) &
             + int5(ri4_79,temp79b,norm79(:,2),f(:,OP_DZ ),g(:,OP_DZZP)) &
             + int5(ri4_79,temp79b,norm79(:,2),f(:,OP_DR ),g(:,OP_DRZP)) &
             - int5(ri4_79,temp79a,norm79(:,1),e(:,OP_DR),h(:,OP_1)) &
             - int5(ri4_79,temp79a,norm79(:,2),e(:,OP_DZ),h(:,OP_1)) &
             - int5(ri4_79,temp79b,f(:,OP_GSP),norm79(:,1),g(:,OP_DR )) &
             - int5(ri4_79,temp79b,f(:,OP_GSP),norm79(:,2),g(:,OP_DZ )) &
             - int5(ri4_79,temp79b,f(:,OP_GS ),norm79(:,1),g(:,OP_DRP)) &
             - int5(ri4_79,temp79b,f(:,OP_GS ),norm79(:,2),g(:,OP_DZP)) &
             - int5(ri4_79,temp79c,f(:,OP_GS ),norm79(:,1),g(:,OP_DR )) &
             - int5(ri4_79,temp79c,f(:,OP_GS ),norm79(:,2),g(:,OP_DZ ))
        if(itor.eq.1) then
           temp = temp &
                - 2.*int5(ri5_79,e(:,OP_1),temp79a,norm79(:,1),h(:,OP_1))
        endif
     else
        temp = int5(ri4_79,e(:,OP_GS),f(:,OP_DZ),g(:,OP_DZP),h(:,OP_1)) &
             + int5(ri4_79,e(:,OP_GS),f(:,OP_DR),g(:,OP_DRP),h(:,OP_1)) &
             + int5(ri4_79,e(:,OP_DZ),f(:,OP_DZP),g(:,OP_GS ),h(:,OP_1 )) &
             + int5(ri4_79,e(:,OP_DR),f(:,OP_DRP),g(:,OP_GS ),h(:,OP_1 )) &
             + int5(ri4_79,e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_GSP),h(:,OP_1 )) &
             + int5(ri4_79,e(:,OP_DR),f(:,OP_DR ),g(:,OP_GSP),h(:,OP_1 )) &
             + int5(ri4_79,e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_GS ),h(:,OP_DP)) &
             + int5(ri4_79,e(:,OP_DR),f(:,OP_DR ),g(:,OP_GS ),h(:,OP_DP))
     end if
  endif
#else
  temp = 0.
#endif
  b1psipsid = temp
  return
  
end function b1psipsid


! B1psibd1
! ========
vectype function b1psibd1(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itwofluid.eq.0) then
     b1psibd1 = 0.
     return
  end if

  if(jadv.eq.0) then
     if(surface_int) then
        temp = 0.
     else
        temp = int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
             - int5(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1))
     endif

  else
     if(surface_int) then
        temp79a = f(:,OP_DZ)*g(:,OP_DR) - f(:,OP_DR)*g(:,OP_DZ)
        temp79b = e(:,OP_1)*h(:,OP_1)
        temp = &
             + int5(ri3_79,e(:,OP_1),temp79a,norm79(:,1),h(:,OP_DR)) &
             + int5(ri3_79,e(:,OP_1),temp79a,norm79(:,2),h(:,OP_DZ)) &
             - int5(ri3_79,temp79b,norm79(:,1),f(:,OP_DRZ),g(:,OP_DR )) &
             + int5(ri3_79,temp79b,norm79(:,1),f(:,OP_DRR),g(:,OP_DZ )) &
             - int5(ri3_79,temp79b,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRR)) &
             + int5(ri3_79,temp79b,norm79(:,1),f(:,OP_DR ),g(:,OP_DRZ)) &
             - int5(ri3_79,temp79b,norm79(:,2),f(:,OP_DZZ),g(:,OP_DR )) &
             + int5(ri3_79,temp79b,norm79(:,2),f(:,OP_DRZ),g(:,OP_DZ )) &
             - int5(ri3_79,temp79b,norm79(:,2),f(:,OP_DZ ),g(:,OP_DRZ)) &
             + int5(ri3_79,temp79b,norm79(:,2),f(:,OP_DR ),g(:,OP_DZZ)) &
             - int5(ri3_79,temp79a,norm79(:,1),e(:,OP_DR),h(:,OP_1)) &
             - int5(ri3_79,temp79a,norm79(:,2),e(:,OP_DZ),h(:,OP_1))
        if(itor.eq.1) then
           temp = temp &
                - int5(ri4_79,e(:,OP_1),temp79a,norm79(:,1),h(:,OP_1))
        endif
     else
        temp = int5(ri3_79,e(:,OP_GS),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
             - int5(ri3_79,e(:,OP_GS),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1))
     end if
  endif

  b1psibd1 = temp
  return
end function b1psibd1

! B1psibd2
! ========
vectype function b1psibd2(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itwofluid.eq.0) then
     b1psibd2 = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     temp = 0.
  else
     if(surface_int) then
        temp79b = e(:,OP_1)*h(:,OP_1)
        temp79c = e(:,OP_1)*h(:,OP_DP)
        temp = &
             + int5(ri5_79,temp79b,norm79(:,2),f(:,OP_DRPP),g(:,OP_1 )) &
             - int5(ri5_79,temp79b,norm79(:,1),f(:,OP_DZPP),g(:,OP_1 )) &
             + int5(ri5_79,temp79b,norm79(:,2),f(:,OP_DRP ),g(:,OP_DP)) &
             - int5(ri5_79,temp79b,norm79(:,1),f(:,OP_DZP ),g(:,OP_DP)) &
             + int5(ri5_79,temp79c,norm79(:,2),f(:,OP_DRP ),g(:,OP_1 )) &
             - int5(ri5_79,temp79c,norm79(:,1),f(:,OP_DZP ),g(:,OP_1 ))
     else
        temp = &
             + int5(ri5_79,e(:,OP_DR),f(:,OP_DZPP),g(:,OP_1 ),h(:,OP_1 )) &
             - int5(ri5_79,e(:,OP_DZ),f(:,OP_DRPP),g(:,OP_1 ),h(:,OP_1 )) &
             + int5(ri5_79,e(:,OP_DR),f(:,OP_DZP ),g(:,OP_DP),h(:,OP_1 )) &
             - int5(ri5_79,e(:,OP_DZ),f(:,OP_DRP ),g(:,OP_DP),h(:,OP_1 )) &
             + int5(ri5_79,e(:,OP_DR),f(:,OP_DZP ),g(:,OP_1 ),h(:,OP_DP)) &
             - int5(ri5_79,e(:,OP_DZ),f(:,OP_DRP ),g(:,OP_1 ),h(:,OP_DP))
     end if
  endif
#else
  temp = 0.
#endif

  b1psibd2 = temp
  return
end function b1psibd2



! B1psifd1
! ========
vectype function b1psifd1(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itwofluid.eq.0) then
     b1psifd1 = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     if(surface_int) then
        temp = 0.
     else
        temp = int5(ri_79,e(:,OP_1),f(:,OP_DZP),g(:,OP_DRP ),h(:,OP_1)) &
             - int5(ri_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_DZP ),h(:,OP_1))
     end if
  else
     if(surface_int) then
        temp79b = e(:,OP_1)*h(:,OP_1)
        temp79c = e(:,OP_1)*h(:,OP_DP)
        temp = &
             - int5(ri3_79,temp79b,norm79(:,2),g(:,OP_DRPP),f(:,OP_GS )) &
             + int5(ri3_79,temp79b,norm79(:,1),g(:,OP_DZPP),f(:,OP_GS )) &
             - int5(ri3_79,temp79b,norm79(:,2),g(:,OP_DRP ),f(:,OP_GSP)) &
             + int5(ri3_79,temp79b,norm79(:,1),g(:,OP_DZP ),f(:,OP_GSP)) &
             - int5(ri3_79,temp79c,norm79(:,2),g(:,OP_DRP ),f(:,OP_GS )) &
             + int5(ri3_79,temp79c,norm79(:,1),g(:,OP_DZP ),f(:,OP_GS ))
     else
        temp = int5(ri3_79,e(:,OP_GS),f(:,OP_DZP),g(:,OP_DRP ),h(:,OP_1)) &
             - int5(ri3_79,e(:,OP_GS),f(:,OP_DRP),g(:,OP_DZP ),h(:,OP_1)) &
             + int5(ri3_79,e(:,OP_DZ),f(:,OP_GSP),g(:,OP_DRP ),h(:,OP_1 )) &
             - int5(ri3_79,e(:,OP_DR),f(:,OP_GSP),g(:,OP_DZP ),h(:,OP_1 )) &
             + int5(ri3_79,e(:,OP_DZ),f(:,OP_GS ),g(:,OP_DRPP),h(:,OP_1 )) &
             - int5(ri3_79,e(:,OP_DR),f(:,OP_GS ),g(:,OP_DZPP),h(:,OP_1 )) &
             + int5(ri3_79,e(:,OP_DZ),f(:,OP_GS ),g(:,OP_DRP ),h(:,OP_DP)) &
             - int5(ri3_79,e(:,OP_DR),f(:,OP_GS ),g(:,OP_DZP ),h(:,OP_DP))
     endif
  endif
  b1psifd1 = temp
#else
  b1psifd1 = 0.
#endif
  return
end function b1psifd1

! B1psifd2
! ========
vectype function b1psifd2(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itwofluid.eq.0) then
     b1psifd2 = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     if(surface_int) then
        temp = 0.
     else
        temp = int5(ri_79,e(:,OP_1),f(:,OP_DZ ),g(:,OP_DRPP),h(:,OP_1)) &
             - int5(ri_79,e(:,OP_1),f(:,OP_DR ),g(:,OP_DZPP),h(:,OP_1))
     end if
  else
     if(surface_int) then
        temp79a = f(:,OP_DZ)*g(:,OP_DRP) - f(:,OP_DR)*g(:,OP_DZP)
        temp79b = e(:,OP_1)*h(:,OP_1)
        temp = int5(ri3_79,e(:,OP_1),temp79a,norm79(:,1),h(:,OP_DR)) &
             + int5(ri3_79,e(:,OP_1),temp79a,norm79(:,2),h(:,OP_DZ)) &
             - int5(ri3_79,temp79b,norm79(:,1),f(:,OP_DRZ),g(:,OP_DRP )) &
             + int5(ri3_79,temp79b,norm79(:,1),f(:,OP_DRR),g(:,OP_DZP )) &
             - int5(ri3_79,temp79b,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRRP)) &
             + int5(ri3_79,temp79b,norm79(:,1),f(:,OP_DR ),g(:,OP_DRZP)) &
             - int5(ri3_79,temp79b,norm79(:,2),f(:,OP_DZZ),g(:,OP_DRP )) &
             + int5(ri3_79,temp79b,norm79(:,2),f(:,OP_DRZ),g(:,OP_DZP )) &
             - int5(ri3_79,temp79b,norm79(:,2),f(:,OP_DZ ),g(:,OP_DRZP)) &
             + int5(ri3_79,temp79b,norm79(:,2),f(:,OP_DR ),g(:,OP_DZZP)) &
             - int5(ri3_79,temp79a,norm79(:,1),e(:,OP_DR),h(:,OP_1))     &
             - int5(ri3_79,temp79a,norm79(:,2),e(:,OP_DZ),h(:,OP_1))
        if(itor.eq.1) then
           temp = temp &
                + int5(ri4_79,e(:,OP_1),temp79a,norm79(:,1),h(:,OP_1))
        endif
     else
        temp = int5(ri3_79,e(:,OP_GS),f(:,OP_DZ ),g(:,OP_DRPP),h(:,OP_1)) &
             - int5(ri3_79,e(:,OP_GS),f(:,OP_DR ),g(:,OP_DZPP),h(:,OP_1))
     endif
  endif
  b1psifd2 = temp
#else
  b1psifd2 = 0.
#endif
  return
end function b1psifd2




! B1bbd
! =====
vectype function b1bbd(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itwofluid.eq.0) then
     b1bbd = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     if(surface_int) then
        temp = 0.
     else
        temp = 0.
     end if
  else
     if(surface_int) then
        temp79a = ri4_79*e(:,OP_1)
        temp = &
             - int5(temp79a,f(:,OP_DP),norm79(:,1),g(:,OP_DR ),h(:,OP_1 )) &
             - int5(temp79a,f(:,OP_DP),norm79(:,2),g(:,OP_DZ ),h(:,OP_1 )) &
             - int5(temp79a,f(:,OP_1 ),norm79(:,1),g(:,OP_DRP),h(:,OP_1 )) &
             - int5(temp79a,f(:,OP_1 ),norm79(:,2),g(:,OP_DZP),h(:,OP_1 )) &
             - int5(temp79a,f(:,OP_1 ),norm79(:,1),g(:,OP_DR ),h(:,OP_DP)) &
             - int5(temp79a,f(:,OP_1 ),norm79(:,2),g(:,OP_DZ ),h(:,OP_DP))
     else
        temp = int5(ri4_79,e(:,OP_DZ),f(:,OP_DZP),g(:,OP_1 ),h(:,OP_1 )) &
             + int5(ri4_79,e(:,OP_DR),f(:,OP_DRP),g(:,OP_1 ),h(:,OP_1 )) &
             + int5(ri4_79,e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_DP),h(:,OP_1 )) &
             + int5(ri4_79,e(:,OP_DR),f(:,OP_DR ),g(:,OP_DP),h(:,OP_1 )) &
             + int5(ri4_79,e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_1 ),h(:,OP_DP)) &
             + int5(ri4_79,e(:,OP_DR),f(:,OP_DR ),g(:,OP_1 ),h(:,OP_DP))
     endif
  endif
  b1bbd = temp
#else
  b1bbd = 0.
#endif
  return
end function b1bbd


! B1bfd1
! ======
vectype function b1bfd1(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itwofluid.eq.0) then
     b1bfd1 = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     if(surface_int) then
        temp = 0.
     else
        temp = int4(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZP),h(:,OP_1)) &
             + int4(e(:,OP_1),f(:,OP_DR),g(:,OP_DRP),h(:,OP_1))
     end if
  else
     if(surface_int) then
        temp79a = f(:,OP_DZ)*g(:,OP_DZP) + f(:,OP_DR)*g(:,OP_DRP)
        temp79b = e(:,OP_1)*h(:,OP_1)
        temp = int5(ri4_79,e(:,OP_1),temp79a,norm79(:,1),h(:,OP_DR)) &
             + int5(ri4_79,e(:,OP_1),temp79a,norm79(:,2),h(:,OP_DZ)) &
             + int5(ri4_79,temp79b,norm79(:,1),f(:,OP_DRZ),g(:,OP_DZP )) &
             + int5(ri4_79,temp79b,norm79(:,1),f(:,OP_DRR),g(:,OP_DRP )) &
             + int5(ri4_79,temp79b,norm79(:,1),f(:,OP_DZ ),g(:,OP_DRZP)) &
             + int5(ri4_79,temp79b,norm79(:,1),f(:,OP_DR ),g(:,OP_DRRP)) &
             + int5(ri4_79,temp79b,norm79(:,2),f(:,OP_DZZ),g(:,OP_DZP )) &
             + int5(ri4_79,temp79b,norm79(:,2),f(:,OP_DRZ),g(:,OP_DRP )) &
             + int5(ri4_79,temp79b,norm79(:,2),f(:,OP_DZ ),g(:,OP_DZZP)) &
             + int5(ri4_79,temp79b,norm79(:,2),f(:,OP_DR ),g(:,OP_DRZP)) &
             - int5(ri4_79,temp79a,norm79(:,1),e(:,OP_DR),h(:,OP_1)) &
             - int5(ri4_79,temp79a,norm79(:,2),e(:,OP_DZ),h(:,OP_1))

        if(itor.eq.1) then
           temp = temp &
                - 2.*int5(ri5_79,e(:,OP_1),temp79a,norm79(:,1),h(:,OP_1))
        endif
     else
        temp = int5(ri2_79,e(:,OP_GS),f(:,OP_DZ),g(:,OP_DZP),h(:,OP_1)) &
             + int5(ri2_79,e(:,OP_GS),f(:,OP_DR),g(:,OP_DRP),h(:,OP_1))
     endif
  endif
  b1bfd1 = temp
#else
  b1bfd1 = 0.
#endif
  return
end function b1bfd1

! B1bfd2
! ======
vectype function b1bfd2(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itwofluid.eq.0) then
     b1bfd2 = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     temp = 0.
  else
     if(surface_int) then
        temp79b = e(:,OP_1)*h(:,OP_1)
        temp79c = e(:,OP_1)*h(:,OP_DP)
        temp = &
             - int5(ri4_79,temp79b,f(:,OP_DP),norm79(:,1),g(:,OP_DRPP)) &
             - int5(ri4_79,temp79b,f(:,OP_DP),norm79(:,2),g(:,OP_DZPP)) &
             - int5(ri4_79,temp79c,f(:,OP_1 ),norm79(:,1),g(:,OP_DRPP)) &
             - int5(ri4_79,temp79c,f(:,OP_1 ),norm79(:,2),g(:,OP_DZPP))

     else
        temp = int5(ri4_79,e(:,OP_DZ),g(:,OP_DZPP),f(:,OP_DP),h(:,OP_1 ))  &
             + int5(ri4_79,e(:,OP_DR),g(:,OP_DRPP),f(:,OP_DP),h(:,OP_1 ))  &
             + int5(ri4_79,e(:,OP_DZ),g(:,OP_DZPP),f(:,OP_1 ),h(:,OP_DP))  &
             + int5(ri4_79,e(:,OP_DR),g(:,OP_DRPP),f(:,OP_1 ),h(:,OP_DP))

#ifdef USECOMPLEX
        ! f''' term hack
        temp = temp + rfac* &
             (int5(ri4_79,e(:,OP_DZ),g(:,OP_DZPP),f(:,OP_1),h(:,OP_1))  &
             +int5(ri4_79,e(:,OP_DR),g(:,OP_DRPP),f(:,OP_1),h(:,OP_1)))
#endif
     endif
  endif
  b1bfd2 = temp
#else
  b1bfd2 = 0.
#endif
  return
end function b1bfd2



! B1ped
! =====
vectype function b1ped(e,f,g)

  use basic
  use m3dc1_nint

  implicit none
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(itwofluid.eq.0) then
     b1ped = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.0) then
     if(surface_int) then
        temp = 0.
     else
        temp = int3(e(:,OP_1),f(:,OP_DP),g(:,OP_1))
     end if
  else
     if(surface_int) then
        temp = int5(ri2_79,e(:,OP_1),f(:,OP_DP),norm79(:,1),g(:,OP_DR)) &
             + int5(ri2_79,e(:,OP_1),f(:,OP_DP),norm79(:,2),g(:,OP_DZ)) &
             - int5(ri2_79,e(:,OP_1),g(:,OP_DP),norm79(:,1),f(:,OP_DR)) &
             - int5(ri2_79,e(:,OP_1),g(:,OP_DP),norm79(:,2),f(:,OP_DZ)) &
             - int5(ri2_79,g(:,OP_1),f(:,OP_DP),norm79(:,1),e(:,OP_DR)) &
             - int5(ri2_79,g(:,OP_1),f(:,OP_DP),norm79(:,2),e(:,OP_DZ))
     else
        temp = int4(ri2_79,e(:,OP_GS),f(:,OP_DP),g(:,OP_1)) &
             + int4(ri2_79,e(:,OP_DZ),f(:,OP_DZP),g(:,OP_1 )) &
             + int4(ri2_79,e(:,OP_DR),f(:,OP_DRP),g(:,OP_1 )) &
             + int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_DP)) &
             + int4(ri2_79,e(:,OP_DR),f(:,OP_DR ),g(:,OP_DP))
     end if
  endif
#else
  temp = 0.
#endif

  b1ped = temp
  return
end function b1ped


! B1feta
! ======
vectype function b1feta(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)

  if(jadv.eq.0) then
     temp = 0.
  else
     if(surface_int) then
        temp = 0.
     else
#ifdef USECOMPLEX
        temp = rfac*int4(ri3_79,e(:,OP_DR),f(:,OP_DZPP),g(:,OP_1)) &
             - rfac*int4(ri3_79,e(:,OP_DZ),f(:,OP_DRPP),g(:,OP_1))
#else
        temp = -int4(ri3_79,e(:,OP_DRP),f(:,OP_DZPP),g(:,OP_1)) &
             +  int4(ri3_79,e(:,OP_DZP),f(:,OP_DRPP),g(:,OP_1))
#endif
     end if
  end if

  b1feta = temp
#else
  b1feta = 0.
#endif
end function b1feta


! b1fu
! ====
vectype function b1fu(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)

  select case(ivform)
  case(0)
     if(jadv.eq.0) then
        if(surface_int) then
           temp = 0.
        else
           temp = &
                - int3(e(:,OP_1),f(:,OP_DZP),g(:,OP_DZ)) &
                - int3(e(:,OP_1),f(:,OP_DRP),g(:,OP_DR))
        endif
     else
        if(surface_int) then
           temp = 0.
        else
           temp = &
                - int4(ri2_79,e(:,OP_GS),f(:,OP_DZP),g(:,OP_DZ)) &
                - int4(ri2_79,e(:,OP_GS),f(:,OP_DRP),g(:,OP_DR))
        endif
     endif
     
  case(1)

     if(jadv.eq.0) then
        if(surface_int) then
           temp = 0.
        else
           temp = &
                - int4(r2_79,e(:,OP_1),f(:,OP_DZP),g(:,OP_DZ)) &
                - int4(r2_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_DR))
        endif
     else
        if(surface_int) then
           temp79a = f(:,OP_DZP)*g(:,OP_DZ) + f(:,OP_DRP)*g(:,OP_DR)
           temp = int3(temp79a,norm79(:,1),e(:,OP_DR)) &
                + int3(temp79a,norm79(:,1),e(:,OP_DZ)) &
                - int4(e(:,OP_1),norm79(:,1),f(:,OP_DRZP),g(:,OP_DZ )) &
                - int4(e(:,OP_1),norm79(:,1),f(:,OP_DRRP),g(:,OP_DR )) &
                - int4(e(:,OP_1),norm79(:,1),f(:,OP_DZP ),g(:,OP_DRZ)) &
                - int4(e(:,OP_1),norm79(:,1),f(:,OP_DRP ),g(:,OP_DRR)) &
                - int4(e(:,OP_1),norm79(:,2),f(:,OP_DZZP),g(:,OP_DZ )) &
                - int4(e(:,OP_1),norm79(:,2),f(:,OP_DRZP),g(:,OP_DR )) &
                - int4(e(:,OP_1),norm79(:,2),f(:,OP_DZP ),g(:,OP_DZZ)) &
                - int4(e(:,OP_1),norm79(:,2),f(:,OP_DRP ),g(:,OP_DRZ))
           if(itor.eq.1) then
              temp = temp &
                   - 2.*int4(ri_79,e(:,OP_1),temp79a,norm79(:,1))
           endif
           
        else
           temp = &
                - int3(e(:,OP_GS),f(:,OP_DZP),g(:,OP_DZ)) &
                - int3(e(:,OP_GS),f(:,OP_DRP),g(:,OP_DR))
        endif
     endif
  end select

  b1fu = temp
#else
  b1fu = 0.
#endif
end function b1fu


! b1fv
! ====
vectype function b1fv(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)

  select case(ivform)
  case(0)
     temp = 0.
     ! not implemented

  case(1)
     if(jadv.eq.0) then
        temp = 0.
     else
        if(surface_int) then
           if(inoslip_tor.eq.1) then
              temp = 0.
           else
           temp = int5(ri_79,e(:,OP_1),norm79(:,1),f(:,OP_DZPP),g(:,OP_1 )) &
                - int5(ri_79,e(:,OP_1),norm79(:,2),f(:,OP_DRPP),g(:,OP_1 )) &
                + int5(ri_79,e(:,OP_1),norm79(:,1),f(:,OP_DZP ),g(:,OP_DP)) &
                - int5(ri_79,e(:,OP_1),norm79(:,2),f(:,OP_DRP ),g(:,OP_DP))
           endif
        else
           temp = int4(ri_79,e(:,OP_DZ),f(:,OP_DRPP),g(:,OP_1)) &
                - int4(ri_79,e(:,OP_DR),f(:,OP_DZPP),g(:,OP_1)) &
                + int4(ri_79,e(:,OP_DZ),f(:,OP_DRP ),g(:,OP_DP)) &
                - int4(ri_79,e(:,OP_DR),f(:,OP_DZP ),g(:,OP_DP))
        endif
     endif
  end select

  b1fv = temp
#else
  b1fv = 0.
#endif
end function b1fv


! b1fchi
! ======
vectype function b1fchi(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)

  select case(ivform)
  case(0)
     if(jadv.eq.0) then
        if(surface_int) then
           temp = 0.
        else
           temp = int4(r_79,e(:,OP_1),f(:,OP_DZP),g(:,OP_DR)) &
                - int4(r_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_DZ))
        endif
     else
        if(surface_int) then
           temp = 0.
        else
           temp = int4(ri_79,e(:,OP_GS),f(:,OP_DZP),g(:,OP_DR)) &
                - int4(ri_79,e(:,OP_GS),f(:,OP_DRP),g(:,OP_DZ))
        endif
     endif

  case(1)

     if(jadv.eq.0) then
        if(surface_int) then
           temp = 0.
        else
           temp = int4(ri_79,e(:,OP_1),f(:,OP_DZP),g(:,OP_DR)) &
                - int4(ri_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_DZ))
        endif
     else
        if(surface_int) then
           temp79a = f(:,OP_DRP)*g(:,OP_DZ) - f(:,OP_DZP)*g(:,OP_DR)
           temp = int4(ri3_79,temp79a,norm79(:,1),e(:,OP_DR)) &
                + int4(ri3_79,temp79a,norm79(:,2),e(:,OP_DZ)) &
                + int5(ri3_79,e(:,OP_1),norm79(:,1),f(:,OP_DRZP),g(:,OP_DR )) &
                - int5(ri3_79,e(:,OP_1),norm79(:,1),f(:,OP_DRRP),g(:,OP_DZ )) &
                + int5(ri3_79,e(:,OP_1),norm79(:,1),f(:,OP_DZP ),g(:,OP_DRR)) &
                - int5(ri3_79,e(:,OP_1),norm79(:,1),f(:,OP_DRP ),g(:,OP_DRZ)) &
                + int5(ri3_79,e(:,OP_1),norm79(:,2),f(:,OP_DZZP),g(:,OP_DR )) &
                - int5(ri3_79,e(:,OP_1),norm79(:,2),f(:,OP_DRZP),g(:,OP_DZ )) &
                + int5(ri3_79,e(:,OP_1),norm79(:,2),f(:,OP_DZP ),g(:,OP_DRZ)) &
                - int5(ri3_79,e(:,OP_1),norm79(:,2),f(:,OP_DRP ),g(:,OP_DZZ))
           if(itor.eq.1) then
              temp = temp &
                   + 2.*int4(ri4_79,e(:,OP_1),temp79a,norm79(:,1))
           endif
           
        else
           temp = int4(ri3_79,e(:,OP_GS),f(:,OP_DZP),g(:,OP_DR)) &
                - int4(ri3_79,e(:,OP_GS),f(:,OP_DRP),g(:,OP_DZ))
        endif
     endif
  end select

  b1fchi = temp
#else
  b1fchi = 0.
#endif
end function b1fchi


! B1e
! ===
vectype function b1e(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(jadv.eq.1) then
     temp = 0.
  else
     if(surface_int) then
        temp = 0.
     else
        temp = -int2(e(:,OP_1),f(:,OP_DP))
     end if
  endif
#else
  temp = 0.
#endif

  b1e = temp
  return
end function b1e


!==============================================================================
! B2 TERMS
!==============================================================================

! B2b
! ===
vectype function b2b(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f
  vectype :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = int3(ri2_79,e(:,OP_1),f(:,OP_1))
  end if

  b2b = temp
  return
end function b2b


! B2psieta
! ========
vectype function b2psieta(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     if(inocurrent_pol.eq.1) then
        temp = 0.
     else
        temp = int5(ri3_79,e(:,OP_1),norm79(:,1),f(:,OP_DZP),g(:,OP_1)) &
             - int5(ri3_79,e(:,OP_1),norm79(:,2),f(:,OP_DRP),g(:,OP_1))
     endif
  else
     temp = int4(ri3_79,e(:,OP_DZ),f(:,OP_DRP),g(:,OP_1)) &
          - int4(ri3_79,e(:,OP_DR),f(:,OP_DZP),g(:,OP_1))

     if(hypi.ne.0) then 
        if(ihypeta.eq.0) then
           temp79a = e(:,OP_DZZ) - e(:,OP_DRR)
           if(itor.eq.1) temp79a = temp79a +    ri_79*e(:,OP_DR)
           temp79b = e(:,OP_DRZ)
           if(itor.eq.1) temp79b = temp79b +    ri_79*e(:,OP_DZ)
           temp79c = e(:,OP_DRZ)
           if(itor.eq.1) temp79c = temp79c - 2.*ri_79*e(:,OP_DZ)
           
           temp = temp + 2.*&
                (int4(ri3_79,temp79a,f(:,OP_DRZP),h(:,OP_1)) &
                -int4(ri3_79,temp79b,f(:,OP_DZZP),h(:,OP_1)) &
                +int4(ri3_79,temp79c,f(:,OP_DRRP),h(:,OP_1)))
           
           if(itor.eq.1) then
              temp = temp - 2.* &
                   (   int4(ri4_79,temp79a,f(:,OP_DZP),h(:,OP_1)) &
                   +2.*int4(ri4_79,temp79c,f(:,OP_DRP),h(:,OP_1)))
           endif
        endif
     end if
  end if
#else
  temp = 0.
#endif
  b2psieta = temp
  return
end function b2psieta


! B2psimue
! ========
vectype function b2psimue(e,f,g)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp
  
  if(surface_int) then
     temp = 0.
  else
     temp = &
          - int4(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ)) &
          - int4(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR)) &
          - int4(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_1 ))
  end if

  b2psimue = temp
  return
end function b2psimue



! B2beta
! ======
vectype function b2beta(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(surface_int) then
     if(inocurrent_pol.eq.1) then
        temp = 0.
     else
        temp = int5(ri2_79,e(:,OP_1),norm79(:,1),f(:,OP_DR),g(:,OP_1)) &
             + int5(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_DZ),g(:,OP_1))
     end if
  else
     temp = &
          - int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1)) &
          - int4(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_1)) 

     if(hypi.ne.0.) then
        if(ihypeta.eq.1) then
           temp79a = (e(:,OP_GS)*g(:,OP_1) + &
                e(:,OP_DZ)*g(:,OP_DZ) + e(:,OP_DR)*g(:,OP_DR))
           temp = temp - int4(ri2_79,temp79a,f(:,OP_GS),h(:,OP_1))
        else
           temp79a = e(:,OP_DZZ) - e(:,OP_DRR)
           if(itor.eq.1) temp79a = temp79a + ri_79*e(:,OP_DR)
           temp79b = e(:,OP_DRZ)
           if(itor.eq.1) temp79b = temp79b - ri_79*e(:,OP_DZ)

           temp = temp &
                - int4(ri2_79,temp79a,f(:,OP_DZZ),h(:,OP_1)) &
                + int4(ri2_79,temp79a,f(:,OP_DRR),h(:,OP_1)) &
                - 2.*int4(ri2_79,e(:,OP_DRZ),f(:,OP_DRZ),h(:,OP_1)) &
                - 2.*int4(ri2_79,temp79b,f(:,OP_DRZ),h(:,OP_1))

           if(itor.eq.1) then
              temp = temp &
                   - int4(ri3_79,temp79a,f(:,OP_DR),h(:,OP_1)) &
                   - 2.*int4(ri3_79,e(:,OP_DRZ),f(:,OP_DZ),h(:,OP_1)) &
                   + 4.*int4(ri3_79,temp79b,f(:,OP_DZ),h(:,OP_1))
           endif

#if defined(USE3D) || defined(USECOMPLEX)
           temp = temp &
                + int4(ri4_79,e(:,OP_DZ),f(:,OP_DZPP),h(:,OP_1)) &
                + int4(ri4_79,e(:,OP_DR),f(:,OP_DRPP),h(:,OP_1))
#endif

        end if
     endif
  end if

  b2beta = temp
  return
end function b2beta


! B2feta
! ======
vectype function b2feta(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     if(inocurrent_pol.eq.1) then
        temp = 0.
     else
        temp = int5(ri2_79,e(:,OP_1),norm79(:,1),f(:,OP_DRPP),g(:,OP_1)) &
             + int5(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_DZPP),g(:,OP_1))
     end if
  else
     temp = &
          - int4(ri2_79,e(:,OP_DZ),f(:,OP_DZPP),g(:,OP_1)) &
          - int4(ri2_79,e(:,OP_DR),f(:,OP_DRPP),g(:,OP_1))

     if(hypi.ne.0.) then
        if(ihypeta.eq.0) then
           temp79a = e(:,OP_DZZ) - e(:,OP_DRR)
           if(itor.eq.1) temp79a = temp79a + ri_79*e(:,OP_DR)
           temp79b = e(:,OP_DRZ)
           if(itor.eq.1) temp79b = temp79b - ri_79*e(:,OP_DZ)

           temp = temp &
                - int4(ri2_79,temp79a,f(:,OP_DZZPP),h(:,OP_1)) &
                + int4(ri2_79,temp79a,f(:,OP_DRRPP),h(:,OP_1)) &
                - 2.*int4(ri2_79,e(:,OP_DRZ),f(:,OP_DRZPP),h(:,OP_1)) &
                - 2.*int4(ri2_79,temp79b,f(:,OP_DRZPP),h(:,OP_1))

           if(itor.eq.1) then
              temp = temp &
                   - int4(ri3_79,temp79a,f(:,OP_DRPP),h(:,OP_1)) &
                   - 2.*int4(ri3_79,e(:,OP_DRZ),f(:,OP_DZPP),h(:,OP_1)) &
                   + 4.*int4(ri3_79,temp79b,f(:,OP_DZPP),h(:,OP_1))
           endif
        endif
     endif
  end if

  b2feta = temp
#else
  b2feta = 0.
#endif
end function b2feta


! B2bu
! ====
vectype function b2bu(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
        
     if(surface_int) then
        temp = 0.
     else
        temp = int4(ri3_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DR)) &
             - int4(ri3_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DZ))
     end if
        
  case(1)
     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp = 0.
        else
           temp = int5(ri_79,e(:,OP_1),f(:,OP_1),norm79(:,1),g(:,OP_DZ)) &
                - int5(ri_79,e(:,OP_1),f(:,OP_1),norm79(:,2),g(:,OP_DR))
        endif
     else
        temp = int4(ri_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DR)) &
             - int4(ri_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DZ))
     endif
  end select

  b2bu = temp
  return
end function b2bu


! B2bchi
! ======
vectype function b2bchi(e,f,g)

  use basic
  use m3dc1_nint

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DZ)) &
             + int4(ri2_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DR))
     end if

  case(1)
     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp = 0.
        else
           temp = &
                - int5(ri4_79,e(:,OP_1),f(:,OP_1),norm79(:,1),g(:,OP_DR)) &
                - int5(ri4_79,e(:,OP_1),f(:,OP_1),norm79(:,2),g(:,OP_DZ))
        endif
     else
        temp = int4(ri4_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DZ)) &
             + int4(ri4_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DR))        
     end if
  end select

  b2bchi = temp
  return
end function b2bchi


! B2bd
! ====
vectype function b2bd(e,f,g)
  use basic
  use m3dc1_nint

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(mass_ratio.eq.0. .or. dbf.eq.0.) then
     b2bd = 0.
     return
  endif

  if(surface_int) then
     temp = 0.
  else
     temp = - &
          (int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1)) &
          +int4(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_1)))
  end if

  b2bd = temp*me_mi*mass_ratio*dbf**2
  return
end function b2bd



! B2psiv
! ======
vectype function b2psiv(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int4(ri3_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_1)) &
             - int4(ri3_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_1))
     end if
     
  case(1)
     if(surface_int) then
        if(inoslip_tor.eq.1) then
           temp = 0.
        else
           temp = int5(ri_79,e(:,OP_1),g(:,OP_1),norm79(:,2),f(:,OP_DR)) &
                - int5(ri_79,e(:,OP_1),g(:,OP_1),norm79(:,1),f(:,OP_DZ))
        endif
     else
        temp = int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_1)) &
             - int4(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_1))
     endif
  end select

  b2psiv = temp
  return
end function b2psiv


! B2fv
! ====
vectype function b2fv(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  select case(ivform)
  case(0)
     if(surface_int) then
        if(inoslip_tor.eq.1) then
           temp = 0.
        else
           temp = int5(ri2_79,e(:,OP_1),norm79(:,2),f(:,OP_DRP),g(:,OP_1)) &
                - int5(ri2_79,e(:,OP_1),norm79(:,1),f(:,OP_DZP),g(:,OP_1))
        endif
     else
        temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_DZP),g(:,OP_1)) &
             + int4(ri2_79,e(:,OP_DR),f(:,OP_DRP),g(:,OP_1))
     endif
     
  case(1)
     if(surface_int) then
        if(inoslip_tor.eq.1) then
           temp = 0.
        else
           temp = int4(e(:,OP_1),norm79(:,2),f(:,OP_DRP),g(:,OP_1)) &
                - int4(e(:,OP_1),norm79(:,1),f(:,OP_DZP),g(:,OP_1))
        endif
     else
        temp = int3(e(:,OP_DZ),f(:,OP_DZP),g(:,OP_1)) &
             + int3(e(:,OP_DR),f(:,OP_DRP),g(:,OP_1))
     endif
  end select
#else
  temp = 0.
#endif

  b2fv = temp
  return
end function b2fv


! B2psipsid
! =========
vectype function b2psipsid(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itwofluid.eq.0) then
     b2psipsid = 0.
     return
  end if

  if(surface_int) then
     if(inocurrent_tor.eq.1) then
        temp = 0.
     else
        temp79a = e(:,OP_1)*f(:,OP_GS)*h(:,OP_1)
        temp = int4(ri3_79,temp79a,norm79(:,2),g(:,OP_DR)) &
             - int4(ri3_79,temp79a,norm79(:,1),g(:,OP_DZ))
     end if
  else
     temp = int5(ri3_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DZ),h(:,OP_1)) &
          - int5(ri3_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DR),h(:,OP_1))
  endif

  b2psipsid = temp
  return 
end function b2psipsid


! B2psibd
! =======
vectype function b2psibd(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itwofluid.eq.0) then
     b2psibd = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     if(inocurrent_norm.eq.1) then
        temp = 0.
     else
        temp79a = e(:,OP_1)*g(:,OP_1)*h(:,OP_1)
        temp = int4(ri4_79,temp79a,norm79(:,1),f(:,OP_DRP)) &
             + int4(ri4_79,temp79a,norm79(:,2),f(:,OP_DZP))
     endif
  else
     temp = &
          -(int5(ri4_79,e(:,OP_DZ),f(:,OP_DZP),g(:,OP_1),h(:,OP_1)) &
           +int5(ri4_79,e(:,OP_DR),f(:,OP_DRP),g(:,OP_1),h(:,OP_1)))
  end if

  b2psibd = temp
#else
  b2psibd = 0.
#endif
  return 
end function b2psibd


! B2bbd
! =====
vectype function b2bbd(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itwofluid.eq.0) then
     b2bbd = 0.
     return
  end if

  if(surface_int) then
     if(inocurrent_norm.eq.1) then
        temp = 0.
     else
        temp79a = e(:,OP_1)*g(:,OP_1)*h(:,OP_1)
        temp = int4(ri3_79,temp79a,norm79(:,2),f(:,OP_DR)) &
             - int4(ri3_79,temp79a,norm79(:,1),f(:,OP_DZ))
     endif
  else
     temp = int5(ri3_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_1),h(:,OP_1)) &
          - int5(ri3_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_1),h(:,OP_1))
  endif
  
  b2bbd = temp
  return 
end function b2bbd



! B2ped
! =====
vectype function b2ped(e,f,g)

  use basic
  use m3dc1_nint

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(itwofluid.eq.0) then
     b2ped = 0.
     return
  end if

  if(surface_int) then
     temp = int5(ri_79,e(:,OP_1),norm79(:,2),f(:,OP_DR),g(:,OP_1)) &
          - int5(ri_79,e(:,OP_1),norm79(:,1),f(:,OP_DZ),g(:,OP_1))
  else
     temp = int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_1)) &
          - int4(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_1))
  end if

  b2ped = temp
  return
end function b2ped


! B2psifd
! =======
vectype function b2psifd(e,f,g,h)

  use basic
  use m3dc1_nint

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itwofluid.eq.0) then
     b2psifd = 0.
     return
  end if
  
#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     if(inocurrent_tor.eq.1) then
        temp = 0.
     else
        temp79a = e(:,OP_1)*f(:,OP_GS)*h(:,OP_1)
        temp = &
             - int4(ri2_79,temp79a,norm79(:,1),g(:,OP_DRP)) &
             - int4(ri2_79,temp79a,norm79(:,2),g(:,OP_DZP))
     endif
  else
     temp = int5(ri2_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DZP),h(:,OP_1)) &
          + int5(ri2_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DRP),h(:,OP_1))
  end if
  b2psifd = temp
#else
  b2psifd = 0.
#endif
  return
end function b2psifd


! B2bfd
! =====
vectype function b2bfd(e,f,g,h)

  use basic
  use m3dc1_nint

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itwofluid.eq.0) then
     b2bfd = 0.
     return
  end if
  
#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     if(inocurrent_norm.eq.1) then
        temp = 0.
     else 
        temp79a = e(:,OP_1)*f(:,OP_1)*h(:,OP_1)
        temp = int4(ri3_79,temp79a,norm79(:,2),g(:,OP_DRPP)) &
             - int4(ri3_79,temp79a,norm79(:,1),g(:,OP_DZPP))
     endif
  else
     temp = - &
          (int5(ri3_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DRPP),h(:,OP_1)) &
          -int5(ri3_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DZPP),h(:,OP_1)))
  end if
  
  b2bfd = temp
#else
  b2bfd = 0.
#endif
  return
end function b2bfd


!=============================================================================
! B3 TERMS
!=============================================================================

! B3pe
! ====
vectype function b3pe(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f
  vectype :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = int2(e(:,OP_1),f(:,OP_1))
  end if

  b3pe = temp
  return
end function b3pe

! B3q
! ===
vectype function b3q(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f
  vectype :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = int2(e(:,OP_1),f(:,OP_1))
  end if

  b3q = temp
  return
end function b3q

! B3psipsieta
! ===========
vectype function b3psipsieta(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  
  vectype :: temp

  if(gam.le.1. .or. surface_int) then
     temp = 0.
  else
     temp = (gam-1.)* &
           int5(ri2_79,e(:,OP_1),f(:,OP_GS), g(:,OP_GS), h(:,OP_1))   
#if defined(USE3D) || defined(USECOMPLEX)
     temp = temp + (gam-1)*   &
           (int5(ri4_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_DRP),h(:,OP_1))   &
         +  int5(ri4_79,e(:,OP_1),f(:,OP_DZP),g(:,OP_DZP),h(:,OP_1)))
#endif
  end if


  b3psipsieta = temp
  
  return
end function b3psipsieta

! B3psibeta
! ===========
vectype function b3psibeta(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  
  vectype :: temp

  if(gam.le.1. .or. surface_int) then
     temp = 0.
  else
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
     temp = 2.*(gam-1.)* &
          (int5(ri3_79,e(:,OP_1),f(:,OP_DZP),g(:,OP_DR),h(:,OP_1))  &
          -int5(ri3_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_DZ),h(:,OP_1)))
#endif
  end if


  b3psibeta = temp
  
  return
end function b3psibeta

! B3psifeta
! ===========
vectype function b3psifeta(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  
  vectype :: temp

  if(gam.le.1. .or. surface_int) then
     temp = 0.
  else
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
     temp = 2.*(gam-1.)* &
          (int5(ri3_79,e(:,OP_1),f(:,OP_DZP),g(:,OP_DRPP),h(:,OP_1))  &
          -int5(ri3_79,e(:,OP_1),f(:,OP_DRP),g(:,OP_DZPP),h(:,OP_1)))
#endif
  end if


  b3psifeta = temp
  
  return
end function b3psifeta


! B3bbeta
! =======
vectype function b3bbeta(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  
  vectype :: temp

  if(gam.le.1. .or. surface_int) then
     temp = 0.
  else 
     temp = (gam-1.)* &
          (int5(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
          +int5(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1)))
  end if

  b3bbeta = temp
  
  return
end function b3bbeta

! B3bfeta
! =======
vectype function b3bfeta(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  
  vectype :: temp

  if(gam.le.1. .or. surface_int) then
     temp = 0.
  else 
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
     temp = (gam-1.)* &
          (int5(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZPP),h(:,OP_1)) &
          +int5(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DRPP),h(:,OP_1)))
#endif
  end if

  b3bfeta = temp
  
  return
end function b3bfeta
! B3ffeta
! =======
vectype function b3ffeta(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  
  vectype :: temp

  if(gam.le.1. .or. surface_int) then
     temp = 0.
  else 
     temp = 0.
#if defined(USE3D) || defined(USECOMPLEX)
     temp = (gam-1.)* &
          (int5(ri2_79,e(:,OP_1),f(:,OP_DZPP),g(:,OP_DZPP),h(:,OP_1)) &
          +int5(ri2_79,e(:,OP_1),f(:,OP_DRPP),g(:,OP_DRPP),h(:,OP_1)))
#endif
  end if

  b3ffeta = temp
  
  return
end function b3ffeta



! B3pepsid
! ========
vectype function b3pepsid(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itwofluid.eq.0 .or. surface_int) then
     b3pepsid = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  temp = int5(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZP),h(:,OP_1)) &
       + int5(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DRP),h(:,OP_1)) &
       - int5(ri2_79,e(:,OP_1),f(:,OP_DP),g(:,OP_GS),h(:,OP_1)) &
       + gam* &
       (int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZP),h(:,OP_DZ)) &
       +int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DRP),h(:,OP_DR)) &
       -int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_GS),h(:,OP_DP)))
#else
  temp = 0.
#endif

  b3pepsid = temp
  
  return
end function b3pepsid


! B3pebd
! ======
vectype function b3pebd(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itwofluid.eq.0 .or. surface_int) then
     b3pebd = 0.
     return
  end if

  temp = int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
        -int5(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1)) &
       + gam* &
       (int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),h(:,OP_DZ)) &
       -int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),h(:,OP_DR)))

  b3pebd = temp
  
  return
end function b3pebd


! B3pefd
! ======
vectype function b3pefd(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(itwofluid.eq.0 .or. surface_int) then
     b3pefd = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  temp = int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DRPP),h(:,OP_1)) &
        -int5(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZPP),h(:,OP_1)) &
       + gam* &
       (int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DRPP),h(:,OP_DZ)) &
       -int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZPP),h(:,OP_DR)))
#else
  temp = 0.
#endif

  b3pefd = temp
  
  return
end function b3pefd



! B3pedkappa
! ==========
vectype function b3pedkappa(e,f,g,h,i)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i
  vectype :: temp

  if(gam.le.1.) then
     b3pedkappa = 0.
     return
  end if

  if(surface_int) then
     temp = int5(e(:,OP_1),norm79(:,1),f(:,OP_DR),g(:,OP_1 ),h(:,OP_1)) &
          + int5(e(:,OP_1),norm79(:,2),f(:,OP_DZ),g(:,OP_1 ),h(:,OP_1)) &
          + int5(e(:,OP_1),norm79(:,1),f(:,OP_1 ),g(:,OP_DR),h(:,OP_1)) &
          + int5(e(:,OP_1),norm79(:,2),f(:,OP_1 ),g(:,OP_DZ),h(:,OP_1))
  else
     temp = &
          - int4(e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1 ),h(:,OP_1)) &
          - int4(e(:,OP_DR),f(:,OP_DR),g(:,OP_1 ),h(:,OP_1)) &
          - int4(e(:,OP_DZ),f(:,OP_1 ),g(:,OP_DZ),h(:,OP_1)) &
          - int4(e(:,OP_DR),f(:,OP_1 ),g(:,OP_DR),h(:,OP_1))
  
#if defined(USE3D) || defined(USECOMPLEX)
     temp = temp +                       &
          int5(ri2_79,e(:,OP_1),f(:,OP_DPP),g(:,OP_1),h(:,OP_1))
#endif
     if(hypp.ne.0.) then
        ! Laplacian[f g]
        temp79a = f(:,OP_LP)*g(:,OP_1) + f(:,OP_1)*g(:,OP_LP) &
             + 2.*(f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR))

        if(ihypkappa.eq.1) then        
           temp = temp - &
                (int4(temp79a,e(:,OP_LP),h(:,OP_1 ),i(:,OP_1)) &
                +int4(temp79a,e(:,OP_DZ),h(:,OP_DZ),i(:,OP_1)) &
                +int4(temp79a,e(:,OP_DR),h(:,OP_DR),i(:,OP_1)))
        else
           temp = temp - int3(temp79a,e(:,OP_LP),i(:,OP_1))
        endif
     endif
  end if

  b3pedkappa = (gam-1.)*temp  
  return
end function b3pedkappa



!============================================================================
! N1 TERMS
!============================================================================

! N1n
! ===
vectype function n1n(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f
  vectype :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = int2(e(:,OP_1),f(:,OP_1))
  end if

  n1n = temp
  return
end function n1n


! N1ndenm
! =======
vectype function n1ndenm(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,h
  real, intent(in) :: g
  vectype :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = g*int2(e(:,OP_1),f(:,OP_LP))

#if defined(USE3D) || defined(USECOMPLEX)
     temp = temp + g*int3(ri2_79,e(:,OP_1),f(:,OP_DPP))
#endif

     if(hypp.ne.0.) then
        if(ihypkappa.eq.1) then
           temp = temp - g*int3(e(:,OP_LP),f(:,OP_LP),h(:,OP_1))
        else
           temp = temp - int3(e(:,OP_LP),f(:,OP_LP),h(:,OP_1))
        endif
     endif
  end if

  n1ndenm = temp
  return
end function n1ndenm



! N1nu
! ====
vectype function n1nu(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp        

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int4(ri_79, e(:,OP_1),f(:,OP_DR),g(:,OP_DZ)) &
             - int4(ri_79, e(:,OP_1),f(:,OP_DZ),g(:,OP_DR))
     end if
     
  case(1)
     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp = 0.
        else
           temp = int5(r_79,e(:,OP_1),f(:,OP_1),norm79(:,2),g(:,OP_DR)) &
                - int5(r_79,e(:,OP_1),f(:,OP_1),norm79(:,1),g(:,OP_DZ))
        endif
     else
        temp = int4(r_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DR)) &
             - int4(r_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DZ))
     endif
  end select

  n1nu = temp
  return
end function n1nu


! N1nv
! ====
vectype function n1nv(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = -int4(ri2_79,e(:,OP_1),f(:,OP_1 ),g(:,OP_DP)) &
             -  int4(ri2_79,e(:,OP_1),f(:,OP_DP),g(:,OP_1 ))
     end if

  case(1)
     if(surface_int) then
        temp = 0.
     else
        temp = -int3(e(:,OP_1),f(:,OP_1 ),g(:,OP_DP)) &
             -  int3(e(:,OP_1),f(:,OP_DP),g(:,OP_1 ))
     end if
  end select

  n1nv = temp
#else
  n1nv = 0.
#endif
  return
end function n1nv

! N1nchi
! ======
vectype function n1nchi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int3(e(:,OP_DZ),f(:,OP_1),g(:,OP_DZ)) &
             + int3(e(:,OP_DR),f(:,OP_1),g(:,OP_DR))
     end if

  case(1)
     if(surface_int) then
        if(inonormalflow.eq.1) then
           temp = 0.
        else
           temp = int5(ri2_79,e(:,OP_1),f(:,OP_1),norm79(:,1),g(:,OP_DR)) &
                + int5(ri2_79,e(:,OP_1),f(:,OP_1),norm79(:,2),g(:,OP_DZ))
        endif
     else
        temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DZ)) &
             + int4(ri2_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DR))
     end if
  end select

  n1nchi = temp
  return
end function n1nchi


! N1s
! ===
vectype function n1s(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f
  vectype :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = int2(e(:,OP_1),f(:,OP_1))
  end if

  n1s = temp
  return
end function n1s



!============================================================================
! P1 TERMS
!============================================================================

! P1pu
! ====
vectype function p1pu(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ)) &
             - int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR))
     end if

  case(1)
     if(surface_int) then
!!$        temp = int5(r_79,e(:,OP_1),f(:,OP_1),norm79(:,1),g(:,OP_DZ)) &
!!$             - int5(r_79,e(:,OP_1),f(:,OP_1),norm79(:,2),g(:,OP_DR))
        temp = 0.
     else
        temp = int4(r_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DR)) &
             - int4(r_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DZ))

        if(itor.eq.1) then
           temp = temp + &
                2.*(gam-1.)*int3(e(:,OP_1),f(:,OP_1),g(:,OP_DZ))
        endif
     end if
  end select

  p1pu = temp

  return
end function p1pu


! P1pv
! ====
vectype function p1pv(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = - int4(ri2_79,e(:,OP_1),f(:,OP_DP),g(:,OP_1)) &
             - gam*int4(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_DP))
     endif
  case(1)
     if(surface_int) then
        temp = 0.
     else
        temp = - int3(e(:,OP_1),f(:,OP_DP),g(:,OP_1)) &
             - gam*int3(e(:,OP_1),f(:,OP_1),g(:,OP_DP))
     endif
  end select
#else
  temp = 0.
#endif

  p1pv = temp

  return
end function p1pv


! P1pchi
! ======
vectype function p1pchi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = gam* &
             (int3(e(:,OP_DZ),f(:,OP_1),g(:,OP_DZ))  &
             +int3(e(:,OP_DR),f(:,OP_1),g(:,OP_DR))) &
             +(gam-1.)* &
             (int3(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ))  &
             +int3(e(:,OP_1),f(:,OP_DR),g(:,OP_DR)))
     end if
     
  case(1)
     if(surface_int) then
!!$        temp = &
!!$             - int5(ri2_79,e(:,OP_1),f(:,OP_1),norm79(:,1),g(:,OP_DR)) &
!!$             - int5(ri2_79,e(:,OP_1),f(:,OP_1),norm79(:,2),g(:,OP_DZ))
        temp = 0.
     else
        temp79a = e(:,OP_DR)*g(:,OP_DR) + e(:,OP_DZ)*g(:,OP_DZ) &
             - (gam-1.)*e(:,OP_1)*g(:,OP_GS)
     
        temp = int3(ri2_79,temp79a,f(:,OP_1))
     endif
  end select

  p1pchi = temp

  return
end function p1pchi


! P1psipsikappar
! ==============
vectype function p1psipsikappar(e,f,g,h,i,j,k)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j,k
  vectype :: temp

  if(gam.le.1.) then
     p1psipsikappar = 0.
     return
  end if

  if(surface_int) then
     temp79a = ri2_79*k(:,OP_1)*j(:,OP_1)*e(:,OP_1)* &
          (norm79(:,1)*f(:,OP_DZ) - norm79(:,2)*f(:,OP_DR))
     temp = int4(temp79a,g(:,OP_DZ),h(:,OP_DR),i(:,OP_1 )) &
          - int4(temp79a,g(:,OP_DR),h(:,OP_DZ),i(:,OP_1 )) &
          + int4(temp79a,g(:,OP_DZ),h(:,OP_1 ),i(:,OP_DR)) &
          - int4(temp79a,g(:,OP_DR),h(:,OP_1 ),i(:,OP_DZ))
  else
     temp79a = k(:,OP_1)*ri2_79* &
          (e(:,OP_DZ)*f(:,OP_DR) - e(:,OP_DR)*f(:,OP_DZ))*j(:,OP_1)

     temp = int4(temp79a,g(:,OP_DZ),h(:,OP_DR),i(:,OP_1 )) &
          - int4(temp79a,g(:,OP_DR),h(:,OP_DZ),i(:,OP_1 )) &
          + int4(temp79a,g(:,OP_DZ),h(:,OP_1 ),i(:,OP_DR)) &
          - int4(temp79a,g(:,OP_DR),h(:,OP_1 ),i(:,OP_DZ))
  end if

  p1psipsikappar = (gam - 1.) * temp
  return
end function p1psipsikappar

! P1psipsipnkappar
! ================
vectype function p1psipsipnkappar(e,f,g,h,i,fac1)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i
  integer, intent(in) :: fac1
  vectype :: temp

  if(gam.le.1. .or. fac1.eq.0) then
     p1psipsipnkappar = 0.
     return
  end if

  temp79a = -ri2_79*kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)

  ! [T,psi]*n = [n p/n^2,psi]*n
  temp79b = ni79(:,OP_1)* &
       (i(:,OP_1)*(h(:,OP_DZ)*g(:,OP_DR)-h(:,OP_DR)*g(:,OP_DZ)) &
       +h(:,OP_1)*(i(:,OP_DZ)*g(:,OP_DR)-i(:,OP_DR)*g(:,OP_DZ))) &
       + 2.*h(:,OP_1)*i(:,OP_1)* &
       (ni79(:,OP_DZ)*g(:,OP_DR) - ni79(:,OP_DR)*g(:,OP_DZ))

  if(surface_int) then
     temp = int5(temp79a,temp79b,e(:,OP_1),norm79(:,1),f(:,OP_DZ)) &
          - int5(temp79a,temp79b,e(:,OP_1),norm79(:,2),f(:,OP_DR))
  else
     temp = int4(temp79a,temp79b,e(:,OP_DZ),f(:,OP_DR)) &
          - int4(temp79a,temp79b,e(:,OP_DR),f(:,OP_DZ))
  end if

  p1psipsipnkappar = (gam - 1.) * temp
  return
end function p1psipsipnkappar


! P1psibkappar
! ============
vectype function p1psibkappar(e,f,g,h,i,j,k)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j,k
  vectype :: temp

  if(gam.le.1.) then
     p1psibkappar = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     temp79a = -ri3_79*k(:,OP_1)*j(:,OP_1)*e(:,OP_1)*g(:,OP_1)* &
          (norm79(:,1)*f(:,OP_DZ) - norm79(:,2)*f(:,OP_DR))

     temp = int3(temp79a,h(:,OP_DP),i(:,OP_1 )) &
          + int3(temp79a,h(:,OP_1 ),i(:,OP_DP))
  else
     temp79a = -k(:,OP_1)*ri3_79*g(:,OP_1)* &
          (e(:,OP_DZ)*f(:,OP_DR) - e(:,OP_DR)*f(:,OP_DZ))*j(:,OP_1)  

     temp79b = f(:,OP_DR)*(h(:,OP_DZ)*i(:,OP_1) + h(:,OP_1)*i(:,OP_DZ)) &
          -    f(:,OP_DZ)*(h(:,OP_DR)*i(:,OP_1) + h(:,OP_1)*i(:,OP_DR))

     temp79c = f(:,OP_DRP)*(h(:,OP_DZ )*i(:,OP_1 ) + h(:,OP_1 )*i(:,OP_DZ )) &
          -    f(:,OP_DZP)*(h(:,OP_DR )*i(:,OP_1 ) + h(:,OP_1 )*i(:,OP_DR )) &
          +    f(:,OP_DR )*(h(:,OP_DZP)*i(:,OP_1 ) + h(:,OP_DP)*i(:,OP_DZ )) &
          -    f(:,OP_DZ )*(h(:,OP_DRP)*i(:,OP_1 ) + h(:,OP_DP)*i(:,OP_DR )) &
          +    f(:,OP_DR )*(h(:,OP_DZ )*i(:,OP_DP) + h(:,OP_1 )*i(:,OP_DZP)) &
          -    f(:,OP_DZ )*(h(:,OP_DR )*i(:,OP_DP) + h(:,OP_1 )*i(:,OP_DRP))

     temp79d = temp79c*g(:,OP_1 )*j(:,OP_1 )*k(:,OP_1 ) &
          +    temp79b*g(:,OP_DP)*j(:,OP_1 )*k(:,OP_1 ) &
          +    temp79b*g(:,OP_1 )*j(:,OP_DP)*k(:,OP_1 ) &
          +    temp79b*g(:,OP_1 )*j(:,OP_1 )*k(:,OP_DP)

     temp = int3(temp79a,h(:,OP_DP),i(:,OP_1 )) &
          + int3(temp79a,h(:,OP_1 ),i(:,OP_DP)) &
          + int3(ri3_79,e(:,OP_1),temp79d)
  end if
#else
  temp = 0.
#endif

  p1psibkappar = (gam - 1.) * temp
  return
end function p1psibkappar

! P1psibpnkappar
! ==============
vectype function p1psibpnkappar(e,f,g,h,i,fac1,fac2)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i
  integer, intent(in) :: fac1, fac2
  vectype :: temp

  if(gam.le.1.) then
     p1psibpnkappar = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*g(:,OP_1)

  ! n*dT/dphi
  temp79e = ni79(:,OP_1)*(h(:,OP_1)*i(:,OP_DP) + h(:,OP_DP)*i(:,OP_1)) &
       + 2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DP)

  if(surface_int) then
     temp = fac1*int5(ri3_79,temp79a,temp79e,norm79(:,2),f(:,OP_DR)) &
          - fac1*int5(ri3_79,temp79a,temp79e,norm79(:,1),f(:,OP_DZ))
  else
     ! d(temp79a)/dphi
     temp79b = kar79(:,OP_DP)*b2i79(:,OP_1 )*ni79(:,OP_1 )*g(:,OP_1 ) &
          +    kar79(:,OP_1 )*b2i79(:,OP_DP)*ni79(:,OP_1 )*g(:,OP_1 ) &
          +    kar79(:,OP_1 )*b2i79(:,OP_1 )*ni79(:,OP_DP)*g(:,OP_1 ) &
          +    kar79(:,OP_1 )*b2i79(:,OP_1 )*ni79(:,OP_1 )*g(:,OP_DP)

     ! [T,psi]*n = [n p/n^2,psi]*n
     temp79c = ni79(:,OP_1)* &
          (i(:,OP_1)*(h(:,OP_DZ)*f(:,OP_DR)-h(:,OP_DR)*f(:,OP_DZ)) &
          +h(:,OP_1)*(i(:,OP_DZ)*f(:,OP_DR)-i(:,OP_DR)*f(:,OP_DZ))) &
          + 2.*h(:,OP_1)*i(:,OP_1)* &
          (ni79(:,OP_DZ)*f(:,OP_DR) - ni79(:,OP_DR)*f(:,OP_DZ))

     ! d(temp79c)/dphi
     temp79d = ni79(:,OP_DP)* &
          (i(:,OP_1)*(h(:,OP_DZ)*f(:,OP_DR)-h(:,OP_DR)*f(:,OP_DZ)) &
          +h(:,OP_1)*(i(:,OP_DZ)*f(:,OP_DR)-i(:,OP_DR)*f(:,OP_DZ))) &
          + 2.*h(:,OP_1)*i(:,OP_1)* &
          (ni79(:,OP_DZP)*f(:,OP_DR) - ni79(:,OP_DRP)*f(:,OP_DZ)) &
          +    ni79(:,OP_1)* &
          (i(:,OP_1)*(h(:,OP_DZ)*f(:,OP_DRP)-h(:,OP_DR)*f(:,OP_DZP)) &
          +h(:,OP_1)*(i(:,OP_DZ)*f(:,OP_DRP)-i(:,OP_DR)*f(:,OP_DZP))) &
          + 2.*h(:,OP_1)*i(:,OP_1)* &
          (ni79(:,OP_DZ)*f(:,OP_DRP) - ni79(:,OP_DR)*f(:,OP_DZP)) &
          +    ni79(:,OP_1)* &
          (i(:,OP_1 )*(h(:,OP_DZP)*f(:,OP_DR)-h(:,OP_DRP)*f(:,OP_DZ)) &
          +h(:,OP_DP)*(i(:,OP_DZ )*f(:,OP_DR)-i(:,OP_DR )*f(:,OP_DZ))) &
          + 2.*h(:,OP_DP)*i(:,OP_1)* &
          (ni79(:,OP_DZ)*f(:,OP_DR) - ni79(:,OP_DR)*f(:,OP_DZ)) &
          +    ni79(:,OP_1)* &
          (i(:,OP_DP)*(h(:,OP_DZ )*f(:,OP_DR)-h(:,OP_DR )*f(:,OP_DZ)) &
          +h(:,OP_1 )*(i(:,OP_DZP)*f(:,OP_DR)-i(:,OP_DRP)*f(:,OP_DZ))) &
          + 2.*h(:,OP_1)*i(:,OP_DP)* &
          (ni79(:,OP_DZ)*f(:,OP_DR) - ni79(:,OP_DR)*f(:,OP_DZ))

     temp = fac2*int4(ri3_79,e(:,OP_1),temp79a,temp79d) &
          + fac2*int4(ri3_79,e(:,OP_1),temp79b,temp79c) &
          + fac1*int5(ri3_79,e(:,OP_DR),temp79a,f(:,OP_DZ),temp79e) &
          - fac1*int5(ri3_79,e(:,OP_DZ),temp79a,f(:,OP_DR),temp79e)
  end if
#else
  temp = 0.
#endif

  p1psibpnkappar = (gam - 1.) * temp
  return
end function p1psibpnkappar


! P1bbkappar
! ==========
vectype function p1bbkappar(e,f,g,h,i,j,k)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j,k
  vectype :: temp

  if(gam.le.1.) then
     p1bbkappar = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     temp = 0.
  else
     temp79a = h(:,OP_DP)*i(:,OP_1) + h(:,OP_1)*i(:,OP_DP)
     temp79b = h(:,OP_DPP)*i(:,OP_1  ) &
          + 2.*h(:,OP_DP )*i(:,OP_DP ) &
          +    h(:,OP_1  )*i(:,OP_DPP)

     temp79c = f(:,OP_DP)*g(:,OP_1 )*temp79a*j(:,OP_1 )*k(:,OP_1 ) &
          +    f(:,OP_1 )*g(:,OP_DP)*temp79a*j(:,OP_1 )*k(:,OP_1 ) &
          +    f(:,OP_1 )*g(:,OP_1 )*temp79b*j(:,OP_1 )*k(:,OP_1 ) &
          +    f(:,OP_1 )*g(:,OP_1 )*temp79a*j(:,OP_DP)*k(:,OP_1 ) &
          +    f(:,OP_1 )*g(:,OP_1 )*temp79a*j(:,OP_1 )*k(:,OP_DP)

     temp = int3(ri4_79,e(:,OP_1),temp79c)
  end if
#else
  temp = 0.
#endif

  p1bbkappar = (gam - 1.) * temp
  return
end function p1bbkappar


! P1bbpnkappar
! ===========
vectype function p1bbpnkappar(e,f,g,h,i,fac1)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i
  integer, intent(in) :: fac1
  vectype :: temp
  

  if(gam.le.1. .or. fac1.eq.0) then
     p1bbpnkappar = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     temp = 0.
  else
     temp79a = kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1)*g(:,OP_1)

     ! d(temp79a)/dphi
     temp79b = kar79(:,OP_DP)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1)*g(:,OP_1) &
          +    kar79(:,OP_1)*b2i79(:,OP_DP)*ni79(:,OP_1)*f(:,OP_1)*g(:,OP_1) &
          +    kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_DP)*f(:,OP_1)*g(:,OP_1) &
          +    kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_DP)*g(:,OP_1) &
          +    kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1)*g(:,OP_DP)

     ! n*dT/dphi
     temp79c = ni79(:,OP_1)*(h(:,OP_1)*i(:,OP_DP) + h(:,OP_DP)*i(:,OP_1)) &
          + 2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DP)
     ! d(temp79c)/dphi
     temp79d = ni79(:,OP_DP)*(h(:,OP_1)*i(:,OP_DP) + h(:,OP_DP)*i(:,OP_1)) &
          + 2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DPP) &
          +    ni79(:,OP_1)*(h(:,OP_DP)*i(:,OP_DP) + h(:,OP_DPP)*i(:,OP_1)) &
          + 2.*h(:,OP_DP)*i(:,OP_1)*ni79(:,OP_DP) &
          +    ni79(:,OP_1)*(h(:,OP_1)*i(:,OP_DPP) + h(:,OP_DP)*i(:,OP_DP)) &
          + 2.*h(:,OP_1)*i(:,OP_DP)*ni79(:,OP_DP)

     temp = int4(ri4_79,e(:,OP_1),temp79a,temp79d) &
          + int4(ri4_79,e(:,OP_1),temp79b,temp79c)
  end if
#else
  temp = 0.
#endif

  p1bbpnkappar = (gam - 1.) * temp
  return
end function p1bbpnkappar



! P1psifkappar
! ============
vectype function p1psifkappar(e,f,g,h,i,j,k)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j,k
  vectype :: temp

  if(gam.le.1.) then
     p1psifkappar = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     temp79a = k(:,OP_1)*ri_79*e(:,OP_1)* &
          (norm79(:,2)*f(:,OP_DR) - norm79(:,1)*f(:,OP_DZ))*j(:,OP_1)
     temp79b = -k(:,OP_1)*ri_79*e(:,OP_1)* &
          (norm79(:,2)*g(:,OP_DZP) + norm79(:,1)*g(:,OP_DRP))*j(:,OP_1)

     temp = int4(temp79a,g(:,OP_DZP),h(:,OP_DZ),i(:,OP_1 )) &
          + int4(temp79a,g(:,OP_DRP),h(:,OP_DR),i(:,OP_1 )) &
          + int4(temp79a,g(:,OP_DZP),h(:,OP_1 ),i(:,OP_DZ)) &
          + int4(temp79a,g(:,OP_DRP),h(:,OP_1 ),i(:,OP_DR)) &
          + int4(temp79b,f(:,OP_DR ),h(:,OP_DZ),i(:,OP_1 )) &
          - int4(temp79b,f(:,OP_DZ ),h(:,OP_DR),i(:,OP_1 )) &
          + int4(temp79b,f(:,OP_DR ),h(:,OP_1 ),i(:,OP_DZ)) &
          - int4(temp79b,f(:,OP_DZ ),h(:,OP_1 ),i(:,OP_DR))
  else
     temp79a = k(:,OP_1)*ri_79* &
          (e(:,OP_DZ)*f(:,OP_DR) - e(:,OP_DR)*f(:,OP_DZ))*j(:,OP_1)
     temp79b = k(:,OP_1)*ri_79* &
          (e(:,OP_DZ)*g(:,OP_DZP) + e(:,OP_DR)*g(:,OP_DRP))*j(:,OP_1)

     temp = int4(temp79a,g(:,OP_DZP),h(:,OP_DZ),i(:,OP_1 )) &
          + int4(temp79a,g(:,OP_DRP),h(:,OP_DR),i(:,OP_1 )) &
          + int4(temp79a,g(:,OP_DZP),h(:,OP_1 ),i(:,OP_DZ)) &
          + int4(temp79a,g(:,OP_DRP),h(:,OP_1 ),i(:,OP_DR)) &
          + int4(temp79b,f(:,OP_DR ),h(:,OP_DZ),i(:,OP_1 )) &
          - int4(temp79b,f(:,OP_DZ ),h(:,OP_DR),i(:,OP_1 )) &
          + int4(temp79b,f(:,OP_DR ),h(:,OP_1 ),i(:,OP_DZ)) &
          - int4(temp79b,f(:,OP_DZ ),h(:,OP_1 ),i(:,OP_DR))
  end if
#else
  temp = 0.
#endif

  p1psifkappar = (gam - 1.) * temp
  return
end function p1psifkappar

! P1psifpnkappar
! ==============
vectype function p1psifpnkappar(e,f,g,h,i,fac1,fac2)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i
  integer, intent(in) :: fac1, fac2
  vectype :: temp

  if(gam.le.1.) then
     p1psifpnkappar = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = ri_79*kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)

  ! n*<T,f'>
  temp79b = ni79(:,OP_1)*h(:,OP_1)* &
       (i(:,OP_DR)*g(:,OP_DRP) + i(:,OP_DZ)*g(:,OP_DZP)) &
       + ni79(:,OP_1)*i(:,OP_1)* &
       (h(:,OP_DR)*g(:,OP_DRP) + h(:,OP_DZ)*g(:,OP_DZP)) &
       + 2.*h(:,OP_1)*i(:,OP_1)* &
       (ni79(:,OP_DR)*g(:,OP_DRP) + ni79(:,OP_DZ)*g(:,OP_DZP))
  ! n*[T,psi]
  temp79c = ni79(:,OP_1)*h(:,OP_1)* &
       (i(:,OP_DZ)*f(:,OP_DR) - i(:,OP_DR)*f(:,OP_DZ)) &
       + ni79(:,OP_1)*i(:,OP_1)* &
       (h(:,OP_DZ)*f(:,OP_DR) - h(:,OP_DR)*f(:,OP_DZ)) &
       + 2.*h(:,OP_1)*i(:,OP_1)* &
       (ni79(:,OP_DZ)*f(:,OP_DR) - ni79(:,OP_DR)*f(:,OP_DZ))
  
  if(surface_int) then
     temp = fac1*int5(temp79a,temp79b,e(:,OP_1),norm79(:,1),f(:,OP_DZ)) &
          - fac1*int5(temp79a,temp79b,e(:,OP_1),norm79(:,2),f(:,OP_DR)) &
          - fac2*int5(temp79a,temp79c,e(:,OP_1),norm79(:,1),g(:,OP_DRP)) &
          - fac2*int5(temp79a,temp79c,e(:,OP_1),norm79(:,2),g(:,OP_DZP))
  else
     temp = fac1*int4(temp79a,e(:,OP_DZ),f(:,OP_DR),temp79b) &
          - fac1*int4(temp79a,e(:,OP_DR),f(:,OP_DZ),temp79b) &
          + fac2*int4(temp79a,e(:,OP_DR),g(:,OP_DRP),temp79c) &
          + fac2*int4(temp79a,e(:,OP_DZ),g(:,OP_DZP),temp79c)
  end if
#else
  temp = 0.
#endif

  p1psifpnkappar = (gam - 1.) * temp
  return
end function p1psifpnkappar

! P1qpsikappar
! ==============
vectype function p1qpsikappar(e,f,g,i,k)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,i,k
  vectype :: temp

  if(gam.le.1.) then
     p1qpsikappar = 0.
     return
  end if

  temp79a = ri_79*k(:,OP_1)*i(:,OP_1)
  
  if(surface_int) then
     temp = 0.
  else
     temp = int4(temp79a,g(:,OP_DZ),e(:,OP_DR),f(:,OP_1)) &
          - int4(temp79a,g(:,OP_DR),e(:,OP_DZ),f(:,OP_1)) 
  end if

  p1qpsikappar = (gam - 1.) * temp
  return
end function p1qpsikappar

! P1bfkappar
! ==========
vectype function p1bfkappar(e,f,g,h,i,j,k)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j,k
  vectype :: temp

  if(gam.le.1.) then
     p1bfkappar = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     temp79a = -ri2_79*k(:,OP_1)*j(:,OP_1)*e(:,OP_1)*f(:,OP_1)* &
          (norm79(:,1)*g(:,OP_DRP) + norm79(:,2)*g(:,OP_DZP))

     temp = int3(temp79a,h(:,OP_DP),i(:,OP_1 )) &
          + int3(temp79a,h(:,OP_1 ),i(:,OP_DP))
  else
     temp79a = k(:,OP_1)*ri2_79*f(:,OP_1)* &
          (e(:,OP_DZ)*g(:,OP_DZP) + e(:,OP_DR)*g(:,OP_DRP))*j(:,OP_1)

     temp79b = g(:,OP_DZP)*(h(:,OP_DZ)*i(:,OP_1) + h(:,OP_1)*i(:,OP_DZ)) &
          +    g(:,OP_DRP)*(h(:,OP_DR)*i(:,OP_1) + h(:,OP_1)*i(:,OP_DR))

     temp79c = g(:,OP_DZPP)*(h(:,OP_DZ )*i(:,OP_1 ) + h(:,OP_1 )*i(:,OP_DZ )) &
          +    g(:,OP_DRPP)*(h(:,OP_DR )*i(:,OP_1 ) + h(:,OP_1 )*i(:,OP_DR )) &
          +    g(:,OP_DZP )*(h(:,OP_DZP)*i(:,OP_1 ) + h(:,OP_DP)*i(:,OP_DZ )) &
          +    g(:,OP_DRP )*(h(:,OP_DRP)*i(:,OP_1 ) + h(:,OP_DP)*i(:,OP_DR )) &
          +    g(:,OP_DZP )*(h(:,OP_DZ )*i(:,OP_DP) + h(:,OP_1 )*i(:,OP_DZP)) &
          +    g(:,OP_DRP )*(h(:,OP_DR )*i(:,OP_DP) + h(:,OP_1 )*i(:,OP_DRP))

     temp79d = temp79c*f(:,OP_1 )*j(:,OP_1 )*k(:,OP_1 ) &
          +    temp79b*f(:,OP_DP)*j(:,OP_1 )*k(:,OP_1 ) &
          +    temp79b*f(:,OP_1 )*j(:,OP_DP)*k(:,OP_1 ) &
          +    temp79b*f(:,OP_1 )*j(:,OP_1 )*k(:,OP_DP)

     temp = int3(temp79a,h(:,OP_DP),i(:,OP_1 )) &
          + int3(temp79a,h(:,OP_1 ),i(:,OP_DP)) &
          - int3(ri2_79,e(:,OP_1),temp79d)
  end if
#else
  temp = 0.
#endif

  p1bfkappar = (gam - 1.) * temp
  return
end function p1bfkappar


! P1qbkappar
! ==========
vectype function p1qbkappar(e,f,g,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,i,j
  vectype :: temp

  if(gam.le.1.) then
     p1qbkappar = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     temp79a = 0.
  else
     temp79a =  ri2_79*i(:,OP_1)*j(:,OP_1)*g(:,OP_1)


     temp = -int3(temp79a,e(:,OP_DP),f(:,OP_1 )) 
  end if
#else
  temp = 0.
#endif

  p1qbkappar = (gam - 1.) * temp
  return
end function p1qbkappar

! ==========
vectype function p1bfpnkappar(e,f,g,h,i,fac1,fac2)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i
  integer, intent(in) :: fac1, fac2
  vectype :: temp

  if(gam.le.1.) then
     p1bfpnkappar = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  temp79a = ri2_79*kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1)

  ! n*dT/dphi = n*d(n p/n^2)/dphi
  temp79e = ni79(:,OP_1)*(h(:,OP_1)*i(:,OP_DP) + h(:,OP_DP)*i(:,OP_1)) &
       + 2.*h(:,OP_1)*i(:,OP_1)*ni79(:,OP_DP)

  if(surface_int) then
     temp = -fac2*int5(temp79a,temp79e,e(:,OP_1),norm79(:,1),g(:,OP_DRP)) &
          -  fac2*int5(temp79a,temp79e,e(:,OP_1),norm79(:,2),g(:,OP_DZP))
  else
     ! d(temp79a)/dphi
     temp79b = ri2_79 * &
          (kar79(:,OP_DP)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_1) &
          +kar79(:,OP_1)*b2i79(:,OP_DP)*ni79(:,OP_1)*f(:,OP_1) &
          +kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_DP)*f(:,OP_1) &
          +kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)*f(:,OP_DP))

     ! n*<T, f'> = n*<n p/n^2, f'>
     temp79c = ni79(:,OP_1)*h(:,OP_1)* &
          (i(:,OP_DR)*g(:,OP_DRP) + i(:,OP_DZ)*g(:,OP_DZP)) &
          + ni79(:,OP_1)*i(:,OP_1)* &
          (h(:,OP_DR)*g(:,OP_DRP) + h(:,OP_DZ)*g(:,OP_DZP)) &
          + 2.*h(:,OP_1)*i(:,OP_1)* &
          (ni79(:,OP_DR)*g(:,OP_DRP) + ni79(:,OP_DZ)*g(:,OP_DZP))

     ! d(temp79c)/dphi
     temp79d = ni79(:,OP_DP)*h(:,OP_1)* &
          (i(:,OP_DR)*g(:,OP_DRP) + i(:,OP_DZ)*g(:,OP_DZP)) &
          + ni79(:,OP_DP)*i(:,OP_1)* &
          (h(:,OP_DR)*g(:,OP_DRP) + h(:,OP_DZ)*g(:,OP_DZP)) &
          + 2.*h(:,OP_1)*i(:,OP_1)* &
          (ni79(:,OP_DRP)*g(:,OP_DRP) + ni79(:,OP_DZP)*g(:,OP_DZP)) &
          +    ni79(:,OP_1)*h(:,OP_1)* &
          (i(:,OP_DR)*g(:,OP_DRPP) + i(:,OP_DZ)*g(:,OP_DZPP)) &
          + ni79(:,OP_1)*i(:,OP_1)* &
          (h(:,OP_DR)*g(:,OP_DRPP) + h(:,OP_DZ)*g(:,OP_DZPP)) &
          + 2.*h(:,OP_1)*i(:,OP_1)* &
          (ni79(:,OP_DR)*g(:,OP_DRPP) + ni79(:,OP_DZ)*g(:,OP_DZPP)) &
          +    ni79(:,OP_1)*h(:,OP_DP)* &
          (i(:,OP_DR)*g(:,OP_DRP) + i(:,OP_DZ)*g(:,OP_DZP)) &
          + ni79(:,OP_1)*i(:,OP_1)* &
          (h(:,OP_DRP)*g(:,OP_DRP) + h(:,OP_DZP)*g(:,OP_DZP)) &
          + 2.*h(:,OP_DP)*i(:,OP_1)* &
          (ni79(:,OP_DR)*g(:,OP_DRP) + ni79(:,OP_DZ)*g(:,OP_DZP)) &
          +    ni79(:,OP_1)*h(:,OP_1)* &
          (i(:,OP_DRP)*g(:,OP_DRP) + i(:,OP_DZP)*g(:,OP_DZP)) &
          + ni79(:,OP_1)*i(:,OP_DP)* &
          (h(:,OP_DR)*g(:,OP_DRP) + h(:,OP_DZ)*g(:,OP_DZP)) &
          + 2.*h(:,OP_1)*i(:,OP_DP)* &
          (ni79(:,OP_DR)*g(:,OP_DRP) + ni79(:,OP_DZ)*g(:,OP_DZP))

     temp = fac2*int4(temp79a,e(:,OP_DR),g(:,OP_DRP),temp79e) &
          + fac2*int4(temp79a,e(:,OP_DZ),g(:,OP_DZP),temp79e) &
          - fac1*int3(e(:,OP_1),temp79a,temp79d) &
          - fac1*int3(e(:,OP_1),temp79b,temp79c)
  end if
#else
  temp = 0.
#endif

  p1bfpnkappar = (gam - 1.) * temp
  return
end function p1bfpnkappar


! P1ffkappar
! ==========
vectype function p1ffkappar(e,f,g,h,i,j,k)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j,k
  vectype :: temp

  if(gam.le.1.) then
     p1ffkappar = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     temp79a =  k(:,OP_1)*e(:,OP_1)* &
          (norm79(:,2)*f(:,OP_DZP) + norm79(:,1)*f(:,OP_DRP))*j(:,OP_1)

     temp = int4(temp79a,g(:,OP_DZP),h(:,OP_DZ),i(:,OP_1 )) &
          + int4(temp79a,g(:,OP_DRP),h(:,OP_DR),i(:,OP_1 )) &
          + int4(temp79a,g(:,OP_DZP),h(:,OP_1 ),i(:,OP_DZ)) &
          + int4(temp79a,g(:,OP_DRP),h(:,OP_1 ),i(:,OP_DR)) 
  else
     temp79a = - k(:,OP_1)*                                            &
          (e(:,OP_DZ)*f(:,OP_DZP) + e(:,OP_DR)*f(:,OP_DRP))*j(:,OP_1)

     temp = int4(temp79a,g(:,OP_DZP),h(:,OP_DZ),i(:,OP_1 )) &
          + int4(temp79a,g(:,OP_DRP),h(:,OP_DR),i(:,OP_1 )) &
          + int4(temp79a,g(:,OP_DZP),h(:,OP_1 ),i(:,OP_DZ)) &
          + int4(temp79a,g(:,OP_DRP),h(:,OP_1 ),i(:,OP_DR)) 
  end if
#else
  temp = 0.
#endif

  p1ffkappar = (gam - 1.) * temp
  return
end function p1ffkappar


! P1ffpnkappar
! ============
vectype function p1ffpnkappar(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i
  vectype :: temp

  if(gam.le.1.) then
     p1ffpnkappar = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     temp = 0.
  else
     temp79a = kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)

     ! <T, f'>*n = <n p/n^2, f'>*n
     temp79c = ni79(:,OP_1)*h(:,OP_1)* &
          (i(:,OP_DR)*g(:,OP_DRP) + i(:,OP_DZ)*g(:,OP_DZP)) &
          + ni79(:,OP_1)*i(:,OP_1)* &
          (h(:,OP_DR)*g(:,OP_DRP) + h(:,OP_DZ)*g(:,OP_DZP)) &
          + 2.*h(:,OP_1)*i(:,OP_1)* &
          (ni79(:,OP_DR)*g(:,OP_DRP) + ni79(:,OP_DZ)*g(:,OP_DZP))

     temp = -int4(temp79a,e(:,OP_DR),f(:,OP_DRP),temp79c) &
          -  int4(temp79a,e(:,OP_DZ),f(:,OP_DZP),temp79c)
  end if
#else
  temp = 0.
#endif

  p1ffpnkappar = (gam - 1.) * temp
  return
end function p1ffpnkappar


! P1qfkappar
! ============
vectype function p1qfkappar(e,f,g,h,i)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i
  vectype :: temp

  if(gam.le.1.) then
     p1qfkappar = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     temp = 0.
  else
     temp79a = i(:,OP_1)*h(:,OP_1)


     temp = int4(temp79a,e(:,OP_DR),g(:,OP_DRP),f(:,OP_1)) &
          + int4(temp79a,e(:,OP_DZ),g(:,OP_DZP),f(:,OP_1))
  end if
#else
  temp = 0.
#endif

  p1qfkappar = (gam - 1.) * temp
  return
end function p1qfkappar

! P1kappax
! ========
vectype function p1kappax(e,f,g,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,i,j
  vectype :: temp

  if(gam.le.1.) then
     p1kappax = 0.
     return
  end if

  if(surface_int) then
     temp = 0.
  else
     temp79a = ri_79*i(:,OP_1)*(e(:,OP_DZ)*f(:,OP_DR) - e(:,OP_DR)*f(:,OP_DZ))

     temp = int3(g(:,OP_1),i(:,OP_1),temp79a) 
  endif

  p1kappax = (gam - 1.) * temp
  return
end function p1kappax



! P1uus
! =====
vectype function p1uus(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(idens.eq.0 .or. gam.le.1. .or. surface_int) then
     p1uus = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = h(:,OP_1) + denm*nt79(:,OP_LP)

  select case(ivform)
  case(0)
     temp = 0.5*(gam-1.)* &
          (int5(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),temp79a) &
          +int5(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),temp79a))
  case(1)
     temp = 0.5*(gam-1.)* &
          (int5(r2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),temp79a) &
          +int5(r2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),temp79a))
  end select

  p1uus = temp
  return
end function p1uus


! P1vvs
! =====
vectype function p1vvs(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(idens.eq.0 .or. gam.le.1. .or. surface_int) then
     p1vvs = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = h(:,OP_1) + denm*nt79(:,OP_LP)

  select case(ivform)
  case(0)
     temp = 0.5*(gam-1.)* &
          int5(ri2_79,e(:,OP_1),f(:,OP_1),g(:,OP_1),temp79a)
  case(1)
     temp = 0.5*(gam-1.)* &
          int5(r2_79,e(:,OP_1),f(:,OP_1),g(:,OP_1),temp79a)
  end select

  p1vvs = temp
  return
end function p1vvs


! P1chichis
! =========
vectype function p1chichis(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(idens.eq.0 .or. gam.le.1. .or. surface_int) then
     p1chichis = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = h(:,OP_1) + denm*nt79(:,OP_LP)

  select case(ivform)
  case(0)
     temp = 0.5*(gam-1.)* &
          (int4(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),temp79a) &
          +int4(e(:,OP_1),f(:,OP_DR),g(:,OP_DR),temp79a))

  case(1)
     temp = 0.5*(gam-1.)* &
          (int5(ri4_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),temp79a) &
          +int5(ri4_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),temp79a))
  end select
     
  p1chichis = temp
  return
end function p1chichis


! P1uchis
! =======
vectype function p1uchis(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(idens.eq.0 .or. gam.le.1. .or. surface_int) then
     p1uchis = 0.
     return
  endif

  ! add in density diffusion explicitly
  temp79a = h(:,OP_1) + denm*nt79(:,OP_LP)

  temp = -(gam-1.)* & 
       (int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),temp79a) &
       -int5(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),temp79a))

  p1uchis = temp
  return
end function p1uchis

!
! Extra diffusion to model upstream differencing
vectype function p1uspu(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp 

  select case(ivform)
  case(0)
!
!....needs to be defined
     if(surface_int) then
        temp = 0.
     else
        temp = 0.
     endif

  case(1)
     if(surface_int) then
!....needs to be defined
        temp = 0.
     else
        temp79a = abs(g(:,OP_DZ))
        temp79b = abs(g(:,OP_DR))
!
        temp = - int4(r_79,e(:,OP_DR),f(:,OP_DR),temp79a)  &
               - int4(r_79,e(:,OP_DZ),f(:,OP_DZ),temp79b)
     end if
  end select

  p1uspu = temp

  return
end function p1uspu

!
! Extra diffusion to model upstream differencing
vectype function p1uspchi(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp 

  select case(ivform)
  case(0)
!
!....needs to be defined
     if(surface_int) then
        temp = 0.
     else
        temp = 0.
     endif

  case(1)
     if(surface_int) then
!....needs to be defined
        temp = 0.
     else
        temp79a = abs(g(:,OP_DR))
        temp79b = abs(g(:,OP_DZ))
!
        temp = - int4(ri2_79,e(:,OP_DR),f(:,OP_DR),temp79a)  &
               - int4(ri2_79,e(:,OP_DZ),f(:,OP_DZ),temp79b)
     end if
  end select

  p1uspchi = temp

  return
end function p1uspchi

! Extra diffusion to model upstream differencing
vectype function p1uspv(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g

  vectype :: temp 

#if defined(USE3D) || defined(USECOMPLEX)
  select case(ivform)
  case(0)
!
!....needs to be defined
     if(surface_int) then
        temp = 0.
     else
        temp = 0.
     endif

  case(1)
     if(surface_int) then
!....needs to be defined
        temp = 0.
     else
        temp79a = abs(g(:,OP_1))
!
        temp =  int3(e(:,OP_1),f(:,OP_DPP),temp79a)  
     end if
  end select

  p1uspv = temp
#else

  p1uspv = 0.
#endif

  return
end function p1uspv


!======================================================================
! Parallel Viscous Terms
!======================================================================

! PVS1
! ====
subroutine PVS1(i,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i
  vectype, intent(out), dimension(MAX_PTS) :: o

  select case(ivform)
  case(0)
     o = 0.

  case(1)

     o = 0.
     if(itor.eq.1) then
        o = o - (1./3.)*i(:,OP_DZ)
     endif
  end select

end subroutine PVS1


! PVS1psipsi
! ==========
subroutine PVS1psipsi(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i,j,k
  vectype, intent(out), dimension(MAX_PTS) :: o


  select case(ivform)
  case(0)
     o = ri_79* &
          (i(:,OP_DRZ) * (j(:,OP_DR)*k(:,OP_DR) - j(:,OP_DZ)*k(:,OP_DZ)) &
          +j(:,OP_DR)*k(:,OP_DZ)*(i(:,OP_DZZ) - i(:,OP_DRR)))
        
     if(itor.eq.1) then
        o = o + ri2_79*(j(:,OP_DZ) * &
             (i(:,OP_DZ)*k(:,OP_DZ) + i(:,OP_DR)*k(:,OP_DR)))
     end if
        
     o = o * ri2_79*b2i79(:,OP_1)

  case(1)
     o = r_79* &
          (j(:,OP_DR)*k(:,OP_DZ)*(i(:,OP_DZZ) - i(:,OP_DRR)) &
          -(j(:,OP_DZ)*k(:,OP_DZ) - j(:,OP_DR)*k(:,OP_DR))*i(:,OP_DRZ))
     if(itor.eq.1) then
        o = o + j(:,OP_DR)* &
             (i(:,OP_DZ)*k(:,OP_DR) - i(:,OP_DR)*k(:,OP_DZ))
     endif

     o = o * ri2_79*b2i79(:,OP_1)
  end select

end subroutine PVS1psipsi

! PVS1psib
! ========
subroutine PVS1psib(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i,j,k
  vectype, intent(out), dimension(MAX_PTS) :: o

  select case(ivform)
  case(0)
     o = 0.

  case(1)

     o = 0.
#if defined(USE3D) || defined(USECOMPLEX)
     o = o + k(:,OP_1)* &
          (j(:,OP_DZ)*i(:,OP_DZP) + j(:,OP_DR)*i(:,OP_DRP))
     o = o * ri2_79*b2i79(:,OP_1)
#endif     
  end select

end subroutine PVS1psib


! PVS1bb
! ======
subroutine PVS1bb(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i, j, k 
  vectype, intent(out), dimension(MAX_PTS) :: o


  select case(ivform)
  case(0)
     o = 0.
     if(itor.eq.1) then
        o = o - ri2_79*(i(:,OP_DZ)*j(:,OP_1)*k(:,OP_1))
     end if
        
     o = o * ri2_79*b2i79(:,OP_1)

  case(1)
     o = 0.

  end select

end subroutine PVS1bb


! PVS2
! ====
subroutine PVS2(i,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i
  vectype, intent(out), dimension(MAX_PTS) :: o

  select case(ivform)
  case(0)
     o = 0.

  case(1)

#if defined(USE3D) || defined(USECOMPLEX)
     o = -i(:,OP_DP)/3.
#else
     o = 0.
#endif
  end select

end subroutine PVS2


! PVS2psipsi
! ==========
subroutine PVS2psipsi(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i,j,k
  vectype, intent(out), dimension(MAX_PTS) :: o

  o = 0.
end subroutine PVS2psipsi


! PVS2psib
! ========
subroutine PVS2psib(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i, j, k
  vectype, intent(out), dimension(MAX_PTS) :: o

  select case(ivform)
  case(0)
     o = ri_79*k(:,OP_1)* &
          (i(:,OP_DZ)*j(:,OP_DR) - i(:,OP_DR)*j(:,OP_DZ))
        
     if(itor.eq.1) then
        o = o + 2.*ri2_79*k(:,OP_1)*j(:,OP_DZ)*i(:,OP_1)
     endif

     o = o * ri2_79*b2i79(:,OP_1)

  case(1)
     o = r_79*k(:,OP_1)* &
          (i(:,OP_DZ)*j(:,OP_DR) - i(:,OP_DR)*j(:,OP_DZ))

     o = o * ri2_79*b2i79(:,OP_1)
  end select

end subroutine PVS2psib

! PVS2bb
! ======
subroutine PVS2bb(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i, j, k
  vectype, intent(out), dimension(MAX_PTS) :: o

  select case(ivform)
  case(0)
     o = 0.

  case(1)
#if defined(USE3D) || defined(USECOMPLEX)
     o = j(:,OP_1)*k(:,OP_1)*i(:,OP_DP)
     o = o * ri2_79*b2i79(:,OP_1)
#else
     o = 0.
#endif
  end select

end subroutine PVS2bb


! PVS3
! ====
subroutine PVS3(i,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i
  vectype, intent(out), dimension(MAX_PTS) :: o

  select case(ivform)
  case(0)
     o = -(1./3.)*i(:,OP_LP)

  case(1)
     o = -(1./3.)*ri2_79*i(:,OP_GS)

     if(itor.eq.1) then
        o = o + (1./3.)*ri3_79*i(:,OP_DR)
     end if
  end select
 
end subroutine PVS3


! PVS3psipsi
! ==========
subroutine PVS3psipsi(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i,j,k
  vectype, intent(out), dimension(MAX_PTS) :: o

  select case(ivform)
  case(0)
     o = i(:,OP_DZZ)*j(:,OP_DR)*k(:,OP_DR) &
       + i(:,OP_DRR)*j(:,OP_DZ)*k(:,OP_DZ) &
          - 2.*i(:,OP_DRZ)*j(:,OP_DR)*k(:,OP_DZ)
  
     o = o * ri2_79*b2i79(:,OP_1)

  case(1)
     o = -ri2_79* &
          (j(:,OP_DZ)*k(:,OP_DZ)*i(:,OP_DZZ) &
          +j(:,OP_DR)*k(:,OP_DR)*i(:,OP_DRR) &
          +2.*j(:,OP_DZ)*k(:,OP_DR)*i(:,OP_DRZ) &
          -i(:,OP_GS)*(j(:,OP_DZ)*k(:,OP_DZ) + j(:,OP_DR)*k(:,OP_DR))) 
  
     if(itor.eq.1) then
        o = o + 2.*ri3_79*j(:,OP_DZ) * &
             (i(:,OP_DZ)*k(:,OP_DR) - i(:,OP_DR)*k(:,OP_DZ))
        
     endif

     o = o * ri2_79*b2i79(:,OP_1)
  end select
 
end subroutine PVS3psipsi


! PVS3psib
! ========
subroutine PVS3psib(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i,j,k
  vectype, intent(out), dimension(MAX_PTS) :: o

  select case(ivform)
  case(0)
     o = 0.

  case(1)
     
#if defined(USE3D) || defined(USECOMPLEX)
     o = ri3_79*k(:,OP_1) * &
          (i(:,OP_DZP)*j(:,OP_DR) - i(:,OP_DRP)*j(:,OP_DZ))
     o = o * ri2_79*b2i79(:,OP_1)
#else
     o = 0.
#endif

  end select
 
end subroutine PVS3psib


! PVS3bb
! ======
subroutine PVS3bb(i,j,k,o)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: i,j,k
  vectype, intent(out), dimension(MAX_PTS) :: o

  select case(ivform)
  case(0)
     o = 0.
     if(itor.eq.1) then
        o = o + ri_79*i(:,OP_DR)*j(:,OP_1)*k(:,OP_1)
     endif

     o = o * ri2_79*b2i79(:,OP_1)

  case(1)
     o = 0.
  end select
 
end subroutine PVS3bb



! PVV1
! ====
subroutine PVV1(e,o)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e
  vectype, intent(out), dimension(MAX_PTS) :: o

  select case(ivform)
  case(0)
     if(surface_int) then
        o = 0.
     else
        o =   (e(:,OP_DZZ) - e(:,OP_DRR))*pst79(:,OP_DR)*pst79(:,OP_DZ) &
             + e(:,OP_DRZ)*(pst79(:,OP_DR)**2 - pst79(:,OP_DZ)**2)
 
        if(itor.eq.1) then
           o = o + ri_79* &
                (pst79(:,OP_DZ)* &
                (e(:,OP_DZ)*pst79(:,OP_DZ) + e(:,OP_DR)*pst79(:,OP_DR)) &
                -2.*(e(:,OP_DZ)*(pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2) &
                +2.*e(:,OP_DR)*pst79(:,OP_DR)*pst79(:,OP_DZ)) &
                -e(:,OP_DZ)*bzt79(:,OP_1)**2)
        endif
        
        o = 3.*ri_79*b2i79(:,OP_1)*o
     endif

  case(1)
     if(surface_int) then
        o = 3.*ri_79*b2i79(:,OP_1)* &
             (e(:,OP_DZ)*pst79(:,OP_DZ) + e(:,OP_DR)*pst79(:,OP_DR))* &
             (norm79(:,1)*pst79(:,OP_DZ) - norm79(:,2)*pst79(:,OP_DR)) &
             - r_79*(norm79(:,1)*e(:,OP_DZ) - norm79(:,2)*e(:,OP_DR))
!        o = 0.
     else
        o =   (e(:,OP_DZZ) - e(:,OP_DRR))*pst79(:,OP_DR)*pst79(:,OP_DZ) &
             + e(:,OP_DRZ)*(pst79(:,OP_DR)**2 - pst79(:,OP_DZ)**2)
        
        if(itor.eq.1) then
           o = o - ri_79* &
                (pst79(:,OP_DZ)* &
                (e(:,OP_DZ)*pst79(:,OP_DZ) + e(:,OP_DR)*pst79(:,OP_DR)) &
                +e(:,OP_DZ)*bzt79(:,OP_1)**2)
        endif
        
        o = 3.*ri_79*b2i79(:,OP_1)*o
        
        if(itor.eq.1) then
           o = o + 2.*e(:,OP_DZ)
        endif
        
#if defined(USECOMPLEX)
     o = o - rfac*3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1) * &
          (e(:,OP_DZ)*pst79(:,OP_DZ) + e(:,OP_DR)*pst79(:,OP_DR))
#endif
     end if
  end select

  o = -o
end subroutine  PVV1

! PVV2
! ====
subroutine PVV2(e,o)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e
  vectype, intent(out), dimension(MAX_PTS) :: o

  select case(ivform)
  case(1)
     if(surface_int) then
!!$     o = 3.*ri_79*b2i79(:,OP_1)*bzt79(:,OP_1)*e(:,OP_1)* &
!!$          (norm79(:,1)*pst79(:,OP_DZ) - norm79(:,2)*pst79(:,OP_DR))
        o = 0.
     else
        o = ri_79*(e(:,OP_DZ)*pst79(:,OP_DR) - e(:,OP_DR)*pst79(:,OP_DZ))

        o = -3.*b2i79(:,OP_1)*bzt79(:,OP_1)*o

#if defined(USECOMPLEX)
        ! This term is a total derivative in phi
        ! and therefore integrates out of the 3D case

        o = o - rfac*e(:,OP_1) * &
             (1. - 3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2)
#endif
     end if
  end select
end subroutine  PVV2


! PVV3
! ====
subroutine PVV3(e,o)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e
  vectype, intent(out), dimension(MAX_PTS) :: o

  select case(ivform)
  case(0)
     if(surface_int) then
        o = 0.
     else
        o = (1. - 3.*ri2_79*b2i79(:,OP_1)* &
             (pst79(:,OP_DR)**2 + pst79(:,OP_DZ)**2)) * e(:,OP_LP) &
             + 3.*ri2_79*b2i79(:,OP_1)* &
             (   e(:,OP_DZZ)*pst79(:,OP_DZ)**2 &
             +   e(:,OP_DRR)*pst79(:,OP_DR)**2 &
             +2.*e(:,OP_DRZ)*pst79(:,OP_DZ)*pst79(:,OP_DR))
        
        if(itor.eq.1) then
           o = o + 3.*ri3_79*b2i79(:,OP_1)*e(:,OP_DR)* &
                (pst79(:,OP_DZ)**2 + pst79(:,OP_DR)**2 - bzt79(:,OP_1)**2)
        endif
     endif

  case(1)
     if(surface_int) then
        o = 0.
     else
        o =      e(:,OP_DZZ)*pst79(:,OP_DR)**2 &
             +   e(:,OP_DRR)*pst79(:,OP_DZ)**2 &
             -2.*e(:,OP_DRZ)*pst79(:,OP_DZ)*pst79(:,OP_DR)
        
        if(itor.eq.1) then
           o = o + 2.*ri_79*pst79(:,OP_DZ)* &
                (e(:,OP_DZ)*pst79(:,OP_DR) - e(:,OP_DR)*pst79(:,OP_DZ)) &
                + ri_79*e(:,OP_DR)*bzt79(:,OP_1)**2
        endif
        
        o = 3.*ri4_79*b2i79(:,OP_1)*o - ri2_79*e(:,OP_GS)
        
#if defined(USECOMPLEX)
     o = o - rfac*3.*ri5_79*b2i79(:,OP_1)*bzt79(:,OP_1) * &
          (e(:,OP_DZ)*pst79(:,OP_DR) - e(:,OP_DR)*pst79(:,OP_DZ))
#endif
     end if
  end select
     
end subroutine  PVV3


! P1vip
! =====
vectype function P1vip(e)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e

  if(amupar.eq.0) then
     P1vip = 0.
     return
  endif

  call PVS1      (pht79,temp79b)
  call PVS1psipsi(pht79,pst79,pst79,temp79c)
  call PVS1psib  (pht79,pst79,bzt79,temp79d)
  call PVS1bb    (pht79,bzt79,bzt79,temp79e)
  temp79a = temp79b + temp79c + temp79d + temp79e

  if(numvar.ge.2) then
     call PVS2      (vzt79,temp79b)
     call PVS2psipsi(vzt79,pst79,pst79,temp79c)
     call PVS2psib  (vzt79,pst79,bzt79,temp79d)
     call PVS2bb    (vzt79,bzt79,bzt79,temp79e)
     temp79a = temp79a + temp79b + temp79c + temp79d + temp79e
  endif

  if(numvar.ge.3) then
     call PVS3      (cht79,temp79b)
     call PVS3psipsi(cht79,pst79,pst79,temp79c)
     call PVS3psib  (cht79,pst79,bzt79,temp79d)
     call PVS3bb    (cht79,bzt79,bzt79,temp79e)
     temp79a = temp79a + temp79b + temp79c + temp79d + temp79e
  endif

  P1vip = 3.*int4(e(:,OP_1),vip79(:,OP_1),temp79a,temp79a)
end function P1vip


!======================================================================
! Gyroviscous terms
!======================================================================

! g1u
! ===
vectype function g1u(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  if(surface_int) then
     g1u = 0.
     return
  end if

  select case(ivform)
  case(0)
     temp79b = e(:,OP_DZZ) - e(:,OP_DRR)
     if(itor.eq.1) temp79b = temp79b - 3.0*ri_79*e(:,OP_DR)
     temp79c = f(:,OP_DRZ)
     if(itor.eq.1) temp79c = temp79c - 0.5*ri_79*f(:,OP_DZ)
     temp79d = f(:,OP_DZZ) - f(:,OP_DRR)
     if(itor.eq.1) temp79d = temp79d +     ri_79*f(:,OP_DR)
     temp79e = e(:,OP_DRZ)
     if(itor.eq.1) temp79e = temp79e + 1.5*ri_79*e(:,OP_DZ)
     
     temp79a = &
          ri_79*b2i79(:,OP_1)*bzt79(:,OP_1) * &
          (1.+1.5*ri2_79*b2i79(:,OP_1)*(pst79(:,OP_DZ)**2 + pst79(:,OP_DR)**2)) * &
          (temp79b*temp79c - temp79d*temp79e)

     if(itor.eq.1) then 
        temp79a = temp79a &
             +2.25*ri4_79*b2i79(:,OP_1)**2*bzt79(:,OP_1) * &
             ((pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2) * &
              (e(:,OP_DZ)*temp79d - f(:,OP_DZ)*temp79b) &
             +4.*pst79(:,OP_DR)*pst79(:,OP_DZ) * &
              (e(:,OP_DZ)*temp79c - f(:,OP_DZ)*temp79e))
     endif

     g1u = int2(pit79(:,OP_1),temp79a)

  case(1)

     temp79b = e(:,OP_DZZ) - e(:,OP_DRR)
     if(itor.eq.1) temp79b = temp79b - ri_79*e(:,OP_DR)
     temp79c = f(:,OP_DZZ) - f(:,OP_DRR)
     if(itor.eq.1) temp79c = temp79c - ri_79*f(:,OP_DR)
     temp79d = e(:,OP_DRZ)
     if(itor.eq.1) temp79d = temp79d + ri_79*e(:,OP_DZ)

     temp79a = 2.*r_79*(5. - 3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2) &
          * bzt79(:,OP_1)*(temp79b*f(:,OP_DRZ) - e(:,OP_DRZ)*temp79c) &
          - 12.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DR)*pst79(:,OP_DZ) &
          * bzt79(:,OP_1)*(e(:,OP_DRZ)*f(:,OP_DZ) - e(:,OP_DZ)*f(:,OP_DRZ))

     if(itor.eq.1) then
        temp79a = temp79a &
             + 2.*(1. + 3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DR)**2) &
             *bzt79(:,OP_1)*(temp79b*f(:,OP_DZ) - e(:,OP_DZ)*temp79c)
     end if

#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = temp79a &
          + (1. - 3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2) &
          * (temp79b*(pst79(:,OP_DZ)*f(:,OP_DZP)-pst79(:,OP_DR)*f(:,OP_DRP)) &
            + 2.*(pst79(:,OP_DZ)*f(:,OP_DRP)*e(:,OP_DRZ) &
                 +pst79(:,OP_DR)*f(:,OP_DZP)*temp79d)) &
          + 3.*ri2_79*b2i79(:,OP_1)* &
            (pst79(:,OP_DZ)*f(:,OP_DZP) + pst79(:,OP_DR)*f(:,OP_DRP))* &
            (temp79b*(pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2) &
            +2.*pst79(:,OP_DZ)*pst79(:,OP_DR)*(e(:,OP_DRZ)+temp79d))

     if(itor.eq.1) then
        temp79a = temp79a + 2.*ri_79*e(:,OP_DZ)* &
             (pst79(:,OP_DZ)*f(:,OP_DRP) - pst79(:,OP_DR)*f(:,OP_DZP))
     end if

     temp79c = f(:,OP_DZZP) - f(:,OP_DRRP)
     if(itor.eq.1) temp79c = temp79c - ri_79*f(:,OP_DRP)
     temp79d = f(:,OP_DRZP)
     if(itor.eq.1) temp79d = temp79d + ri_79*f(:,OP_DZP)

     temp79a = temp79a + 2.* &
          (temp79c* &
          ((1.-3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DR)**2) &
          *e(:,OP_DR)*pst79(:,OP_DR) &
          +(1.-3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DZ)**2) &
          *e(:,OP_DZ)*pst79(:,OP_DZ)) &
          +(1.+3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DR)**2) &
          *e(:,OP_DR)*pst79(:,OP_DZ)*temp79d &
          +(1.+3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DZ)**2) &
          *e(:,OP_DZ)*pst79(:,OP_DR)*f(:,OP_DRZP) &
          -3.*ri2_79*b2i79(:,OP_1) &
          *((bzt79(:,OP_1)**2 - pst79(:,OP_DZ)**2) &
          * e(:,OP_DZ)*pst79(:,OP_DR)*temp79d &
          + (bzt79(:,OP_1)**2 - pst79(:,OP_DR)**2) &
          * e(:,OP_DR)*pst79(:,OP_DZ)*f(:,OP_DRZP)) &
          +(1.-3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2)*bzt79(:,OP_1) &
          *ri_79*(e(:,OP_DZ)*f(:,OP_DRPP) - e(:,OP_DR)*f(:,OP_DZPP)))
     
#endif

     g1u = 0.25*int3(temp79a,pit79(:,OP_1),b2i79(:,OP_1))

  end select
  
  return
end function g1u

! g1v
! ===
vectype function g1v(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  if(surface_int) then
     g1v = 0.
     return
  end if

  select case(ivform)
  case(0)
     temp79b = f(:,OP_DR)
     if(itor.eq.1) temp79b = temp79b - 2.0*ri_79*f(:,OP_1)
     temp79c = e(:,OP_DRR) - e(:,OP_DZZ)
     if(itor.eq.1) temp79c = temp79c + 3.0*ri_79*e(:,OP_DR)
     temp79d = e(:,OP_DRZ)
     if(itor.eq.1) temp79d = temp79d + 3.0*ri_79*e(:,OP_DZ)

     temp79a = &
          0.25*ri_79*b2i79(:,OP_1)*(1.-3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2) * &
          ((pst79(:,OP_DR)*f(:,OP_DZ) + pst79(:,OP_DZ)*temp79b)*temp79c &
          +2.*f(:,OP_DZ )*pst79(:,OP_DZ)*temp79d &
          -2.*e(:,OP_DRZ)*pst79(:,OP_DR)*temp79b)

     temp79b = pst79(:,OP_DZ)*f(:,OP_DR) - pst79(:,OP_DR)*f(:,OP_DZ)
     if(itor.eq.1) temp79b = temp79b - 2.*ri_79*f(:,OP_1)*pst79(:,OP_DZ)
     temp79c = e(:,OP_DZZ)
     if(itor.eq.1) temp79c = temp79c - 2.0*ri_79*e(:,OP_DR)
     temp79d = e(:,OP_DRZ)
     if(itor.eq.1) temp79d = temp79d + 1.5*ri_79*e(:,OP_DZ)
     temp79e = e(:,OP_DRR)
     if(itor.eq.1) temp79e = temp79e +     ri_79*e(:,OP_DR)
     

     temp79a = temp79a &
          -0.75*ri3_79*b2i79(:,OP_1)**2*temp79b * &
          (e(:,OP_GS)*(pst79(:,OP_DZ)**2 + pst79(:,OP_DR)**2) - &
           2.*(   temp79c*pst79(:,OP_DR)**2 &
              -2.*temp79d*pst79(:,OP_DR)*pst79(:,OP_DZ) &
              +   temp79e*pst79(:,OP_DZ)**2))

     if(itor.eq.1) then
        temp79a = temp79a &
             + 4.5*ri4_79*b2i79(:,OP_1)**2*bzt79(:,OP_1)**2*e(:,OP_DZ) * &
             (pst79(:,OP_DZ)*f(:,OP_DZ) + pst79(:,OP_DR)*f(:,OP_DR) &
             -2.*ri_79*f(:,OP_1)*pst79(:,OP_DR))
     endif

     g1v = int2(pit79(:,OP_1),temp79a)

  case(1)

     temp79b = e(:,OP_DZZ) - e(:,OP_DRR)
     if(itor.eq.1) temp79b = temp79b - ri_79*e(:,OP_DR)
     temp79c = e(:,OP_DRZ)
     if(itor.eq.1) temp79c = temp79c + ri_79*e(:,OP_DZ)
     temp79d = r_79*(f(:,OP_DZ)*pst79(:,OP_DR) - f(:,OP_DR)*pst79(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
     temp79d = temp79d + 2.*bzt79(:,OP_1)*f(:,OP_DP)
#endif

     temp79a = -r_79*(1.-3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2) &
          * (2.*e(:,OP_DRZ)* &
             (pst79(:,OP_DR)*f(:,OP_DR) - pst79(:,OP_DZ)*f(:,OP_DZ)) &
            +temp79b* &
             (pst79(:,OP_DZ)*f(:,OP_DR) + pst79(:,OP_DR)*f(:,OP_DZ))) &
          + 3.*ri2_79*b2i79(:,OP_1)*temp79d &
          * ((pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2)*temp79b &
            +2.*pst79(:,OP_DZ)*pst79(:,OP_DR)*(e(:,OP_DRZ)+temp79c))

     if(itor.eq.1) then
        temp79a = temp79a + 2.*e(:,OP_DZ)* &
             (pst79(:,OP_DZ)*f(:,OP_DZ) &
             +3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2 &
             *pst79(:,OP_DR)*f(:,OP_DR))
     endif

#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = temp79a + 2.* &
          ((1.-3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2) &
          *bzt79(:,OP_1)*(e(:,OP_DZ)*f(:,OP_DZP) + e(:,OP_DR)*f(:,OP_DRP)) &
          +(1.+3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2) &
          *ri_79*f(:,OP_DPP) &
          *(e(:,OP_DZ)*pst79(:,OP_DR) - e(:,OP_DR)*pst79(:,OP_DZ)))
#endif

     g1v = 0.25*int3(temp79a,pit79(:,OP_1),b2i79(:,OP_1))
  end select
       
  return
end function g1v

! g1chi
! =====
vectype function g1chi(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  if(surface_int) then
     g1chi = 0.
     return
  end if

  select case(ivform)
  case(0)
     temp79b = e(:,OP_DRZ)
     if(itor.eq.1) temp79b = temp79b + 3.0*ri_79*e(:,OP_DZ)
     temp79c = e(:,OP_DRR) - e(:,OP_DZZ)
     if(itor.eq.1) temp79c = temp79c + 3.0*ri_79*e(:,OP_DR)
     temp79d = f(:,OP_DZZ)
     if(itor.eq.1) temp79d = temp79d -     ri_79*f(:,OP_DR)
     temp79e = f(:,OP_DRR)
     if(itor.eq.1) temp79e = temp79e -     ri_79*f(:,OP_DR)

     temp79a = &
          0.5*b2i79(:,OP_1)*bzt79(:,OP_1) * &
           (temp79c*(f(:,OP_DRR) - f(:,OP_DZZ)) &
           +2.*(e(:,OP_DRZ) + temp79b)*f(:,OP_DRZ)) &
          +1.5*ri2_79*b2i79(:,OP_1)**2*bzt79(:,OP_1) * &
           (temp79c*(f(:,OP_GS)*(pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2) &
           -f(:,OP_DZZ)*pst79(:,OP_DZ)**2 &
           +f(:,OP_DRR)*pst79(:,OP_DR)**2) &
           +2.*(f(:,OP_DRZ)*(pst79(:,OP_DZ)**2*e(:,OP_DRZ) &
                            +pst79(:,OP_DR)**2*temp79b) &
               -pst79(:,OP_DR)*pst79(:,OP_DZ)* &
                (e(:,OP_DRZ)*temp79d + temp79b*temp79e)))
       
     g1chi = int2(pit79(:,OP_1),temp79a)

  case(1)

     temp79b = e(:,OP_DZZ) - e(:,OP_DRR)
     if(itor.eq.1) temp79b = temp79b - ri_79*e(:,OP_DR)
     temp79c = e(:,OP_DRZ)
     if(itor.eq.1) temp79c = temp79c + ri_79*e(:,OP_DZ)
     temp79d = f(:,OP_DZZ)
     if(itor.eq.1) temp79d = temp79d - ri_79*f(:,OP_DR)
     temp79e = f(:,OP_DRR)
     if(itor.eq.1) temp79e = temp79e - 3.*ri_79*f(:,OP_DR)
     temp79f = f(:,OP_DRZ)
     if(itor.eq.1) temp79f = temp79f - ri_79*f(:,OP_DZ)

     temp79a = 2.*ri2_79*bzt79(:,OP_1) * &
          (temp79b*(temp79d - temp79e) + 2.*temp79f*(e(:,OP_DRZ) + temp79c)) &
          + 12.*ri4_79*b2i79(:,OP_1)*bzt79(:,OP_1) * &
          (temp79f*(pst79(:,OP_DZ)**2*e(:,OP_DRZ)+pst79(:,OP_DR)**2*temp79c) &
          -pst79(:,OP_DZ)*pst79(:,OP_DR)* &
          (temp79e*temp79c + temp79d*e(:,OP_DRZ))) &
          + 6.*ri4_79*b2i79(:,OP_1)*bzt79(:,OP_1)*temp79b* &
          (pst79(:,OP_DR)**2*temp79d - pst79(:,OP_DZ)**2*temp79e)

#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = temp79a &
          + ri3_79*(1. - 3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2) * &
          (2.*(pst79(:,OP_DZ)*f(:,OP_DZP)*e(:,OP_DRZ) &
              -pst79(:,OP_DR)*f(:,OP_DRP)*temp79c) &
          -temp79b*(pst79(:,OP_DZ)*f(:,OP_DRP) + pst79(:,OP_DR)*f(:,OP_DZP))) &
          + 3.*ri5_79*b2i79(:,OP_1) * &
          (f(:,OP_DZP)*pst79(:,OP_DR) - f(:,OP_DRP)*pst79(:,OP_DZ)) * &
          (temp79b*(pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2) &
          +2.*pst79(:,OP_DZ)*pst79(:,OP_DR)*(e(:,OP_DRZ) + temp79c))

     if(itor.eq.1) then
        temp79a = temp79a + 2.*ri4_79*e(:,OP_DZ) * &
             (pst79(:,OP_DZ)*f(:,OP_DZP) + pst79(:,OP_DR)*f(:,OP_DRP))
     endif

     temp79d = f(:,OP_DZZP)
     if(itor.eq.1) temp79d = temp79d - ri_79*f(:,OP_DRP)
     temp79e = f(:,OP_DRRP)
     if(itor.eq.1) temp79e = temp79e - 3.*ri_79*f(:,OP_DRP)
     temp79f = f(:,OP_DRZP)
     if(itor.eq.1) temp79f = temp79f - ri_79*f(:,OP_DZP)

     temp79a = temp79a + 2.* &
          (-ri3_79*(1.+3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DZ)**2) &
          *e(:,OP_DZ)*pst79(:,OP_DR)*temp79e &
          + ri3_79*(1.+3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DR)**2) &
          *e(:,OP_DR)*pst79(:,OP_DZ)*temp79d &
          + 2.*ri3_79*temp79f &
          *((1.-3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DZ)**2) &
          * e(:,OP_DZ)*pst79(:,OP_DZ) &
          - (1.-3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DR)**2) &
          * e(:,OP_DR)*pst79(:,OP_DR)) &
          +3.*ri5_79*b2i79(:,OP_1) &
          *((bzt79(:,OP_1)**2 - pst79(:,OP_DR)**2) &
          * e(:,OP_DR)*pst79(:,OP_DZ)*temp79e &
          - (bzt79(:,OP_1)**2 - pst79(:,OP_DZ)**2) &
          * e(:,OP_DZ)*pst79(:,OP_DR)*temp79d) &
          +ri4_79*(1.-3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2) &
          *bzt79(:,OP_1)*(e(:,OP_DZ)*f(:,OP_DZPP) + e(:,OP_DR)*f(:,OP_DRPP)))
#endif

     g1chi = 0.25*int3(temp79a,pit79(:,OP_1),b2i79(:,OP_1))
  end select

  return
end function g1chi


! g2u
! ===
vectype function g2u(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  if(surface_int) then
     g2u = 0.
     return
  end if

  select case(ivform)
  case(0)
     temp79b = f(:,OP_DZZ) - f(:,OP_DRR)
     if(itor.eq.1) temp79b = temp79b +     ri_79*f(:,OP_DR)
     temp79c = f(:,OP_DRZ)
     if(itor.eq.1) temp79c = temp79c +     ri_79*f(:,OP_DZ)
     temp79d = f(:,OP_DRZ)
     if(itor.eq.1) temp79d = temp79d - 2.0*ri_79*f(:,OP_DZ)
     temp79e = f(:,OP_DRZ)
     if(itor.eq.1) temp79e = temp79e - 0.5*ri_79*f(:,OP_DZ)

     temp79a = 0.25*ri_79*b2i79(:,OP_1)* &
          (1.-3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2)* &
           ((e(:,OP_DR)*pst79(:,OP_DZ) + e(:,OP_DZ)*pst79(:,OP_DR))*temp79b &
           -2.*e(:,OP_DZ)*pst79(:,OP_DZ)*temp79c &
           +2.*e(:,OP_DR)*pst79(:,OP_DR)*temp79d) &
          + 0.75*ri3_79*b2i79(:,OP_1)**2 * &
           (pst79(:,OP_DZ)*e(:,OP_DR) - pst79(:,OP_DR)*e(:,OP_DZ)) * &
           (   temp79b*(pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2) &
           +4.*temp79e* pst79(:,OP_DR)    * pst79(:,OP_DZ)   )

     if(itor.eq.1) then
        temp79a = temp79a &
             -4.5*ri4_79*b2i79(:,OP_1)**2*bzt79(:,OP_1)**2 * &
              f(:,OP_DZ)*(e(:,OP_DZ)*pst79(:,OP_DZ) + &
                          e(:,OP_DR)*pst79(:,OP_DR))
     endif

     g2u = int2(pit79(:,OP_1),temp79a)

  case(1)


     temp79b = f(:,OP_DZZ) - f(:,OP_DRR)
     if(itor.eq.1) temp79b = temp79b - ri_79*f(:,OP_DR)
     temp79c = f(:,OP_DRZ)
     if(itor.eq.1) temp79c = temp79c + ri_79*f(:,OP_DZ)
     
     temp79a = 2.*r_79*temp79b* &
          ((1.-3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DZ)**2) &
          *pst79(:,OP_DZ)*e(:,OP_DR) &
          +(1.-3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DR)**2) &
          *pst79(:,OP_DR)*e(:,OP_DZ)) &
          + 2.*r_79* &
          ((1.+3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DR)**2) &
          *e(:,OP_DZ)*pst79(:,OP_DZ)*temp79c &
          -(1.+3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DZ)**2) &
          *e(:,OP_DR)*pst79(:,OP_DR)*f(:,OP_DRZ)) &
          + 6.*ri_79*b2i79(:,OP_1)* &
          ((bzt79(:,OP_1)**2 - pst79(:,OP_DZ)**2) &
          *e(:,OP_DR)*pst79(:,OP_DR)*temp79c &
          -(bzt79(:,OP_1)**2 - pst79(:,OP_DR)**2) &
          *e(:,OP_DZ)*pst79(:,OP_DZ)*f(:,OP_DRZ))

#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = temp79a &
          - 2.*(1. - 3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2) &
          *bzt79(:,OP_1)*(e(:,OP_DZ)*f(:,OP_DZP) + e(:,OP_DR)*f(:,OP_DRP))

     temp79b = f(:,OP_DZZP) - f(:,OP_DRRP)
     if(itor.eq.1) temp79b = temp79b - ri_79*f(:,OP_DRP)
     temp79c = f(:,OP_DRZP)
     if(itor.eq.1) temp79c = temp79c + ri_79*f(:,OP_DZP)

     temp79a = temp79a - 2.*e(:,OP_1)* &
          (3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1) &
          *(pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2)*temp79b &
          +2.*pst79(:,OP_DR)*pst79(:,OP_DZ)*(f(:,OP_DRZP) + temp79c) &
          -(1.+3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2) &
          *ri_79*(f(:,OP_DZPP)*pst79(:,OP_DR) - f(:,OP_DRPP)*pst79(:,OP_DZ)))

#endif

     g2u = -0.25*int3(temp79a,pit79(:,OP_1),b2i79(:,OP_1))
  end select

  return
end function g2u


! g2v
! ===
vectype function g2v(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  if(surface_int) then
     g2v = 0.
     return
  end if

  select case(ivform)
  case(0)
     temp79b = e(:,OP_DZ)*f(:,OP_DR) - e(:,OP_DR)*f(:,OP_DZ)
     if(itor.eq.1) temp79b = temp79b - 2.*ri_79*e(:,OP_DZ)*f(:,OP_1)
     
!!$     temp79a = 0.25*ri_79*b2i79(:,OP_1)*bzt79(:,OP_1)*temp79b* &
!!$          (1.-3.*ri2_79*b2i79(:,OP_1)* &
!!$          (pst79(:,OP_DZ)**2 + pst79(:,OP_DR)**2 - bzt79(:,OP_1)**2))

     temp79a = ri_79*b2i79(:,OP_1)*bzt79(:,OP_1)*temp79b* &
          (1. - 1.5*ri2_79*b2i79(:,OP_1)*(pst79(:,OP_DR)**2 + pst79(:,OP_DZ)**2))
     g2v = int2(pit79(:,OP_1),temp79a)
  case(1)

     temp79a = 2.*r_79*(1. - 3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2) &
          *bzt79(:,OP_1)*(e(:,OP_DZ)*f(:,OP_DR)-e(:,OP_DR)*f(:,OP_DZ))

#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = temp79a &
          - 2.*(1. + 3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2) &
          *(e(:,OP_DZ)*pst79(:,OP_DZ) + e(:,OP_DR)*pst79(:,OP_DR)) &
          *f(:,OP_DP)

     temp79a = temp79a - 2.*e(:,OP_1)* &
          (1.+3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2) &
          *(pst79(:,OP_DZ)*f(:,OP_DZP) + pst79(:,OP_DR)*f(:,OP_DRP))
#endif

     g2v = -0.25*int3(temp79a,pit79(:,OP_1),b2i79(:,OP_1))
  end select

  return
end function g2v


! g2chi
! =====
vectype function g2chi(e,f)

  use basic
  use m3dc1_nint

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  if(surface_int) then
     g2chi = 0.
     return
  end if

  select case(ivform)
  case(0)
     temp79b = f(:,OP_DZZ)
     if(itor.eq.1) temp79b = temp79b - ri_79*f(:,OP_DR)
     temp79c = f(:,OP_DRR)
     if(itor.eq.1) temp79c = temp79c - ri_79*f(:,OP_DR)

     temp79a = &
          -0.5*b2i79(:,OP_1)* &
           (e(:,OP_DZ)*temp79b*pst79(:,OP_DZ) &
           +e(:,OP_DR)*temp79c*pst79(:,OP_DR) &
           +f(:,OP_DRZ)*(e(:,OP_DR)*pst79(:,OP_DZ) + e(:,OP_DZ)*pst79(:,OP_DR))) &
          +1.5*ri2_79*b2i79(:,OP_1)**2 * &
           (pst79(:,OP_DZ)*e(:,OP_DR) - pst79(:,OP_DR)*e(:,OP_DZ)) * &
           (pst79(:,OP_DZ)*pst79(:,OP_DR)*(f(:,OP_DZZ) - f(:,OP_DRR)) &
           -(pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2)*f(:,OP_DRZ)) &
          -1.5*ri2_79*b2i79(:,OP_1)**2*bzt79(:,OP_1)**2 * &
           (e(:,OP_DZ)*temp79c*pst79(:,OP_DZ) &
           +e(:,OP_DR)*temp79b*pst79(:,OP_DR) &
           -f(:,OP_DRZ)*(e(:,OP_DR)*pst79(:,OP_DZ) + &
                         e(:,OP_DZ)*pst79(:,OP_DR)))

     g2chi = int2(pit79(:,OP_1),temp79a)
  case(1)

     temp79b = f(:,OP_DRR)
     if(itor.eq.1) temp79b = temp79b - 3.*ri_79*f(:,OP_DR)
     temp79c = f(:,OP_DZZ)
     if(itor.eq.1) temp79c = temp79c -    ri_79*f(:,OP_DR)
     temp79d = f(:,OP_DRZ)
     if(itor.eq.1) temp79d = temp79d -    ri_79*f(:,OP_DZ)

     temp79a = 2.*ri2_79* &
          ((1. + 3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DZ)**2) &
          *e(:,OP_DR)*pst79(:,OP_DR)*temp79b  &
          +(1. + 3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DR)**2) &
          *e(:,OP_DZ)*pst79(:,OP_DZ)*temp79c) &
          - 4.*ri2_79*temp79d* &
          ((1. - 3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DZ)**2) &
          *e(:,OP_DR)*pst79(:,OP_DZ) &
          +(1. - 3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DR)**2) &
          *e(:,OP_DZ)*pst79(:,OP_DR)) &
          + 6.*ri4_79*b2i79(:,OP_1)* &
          ((bzt79(:,OP_1)**2 - pst79(:,OP_DZ)**2) &
          *e(:,OP_DR)*pst79(:,OP_DR)*temp79c &
          +(bzt79(:,OP_1)**2 - pst79(:,OP_DR)**2) &
          *e(:,OP_DZ)*pst79(:,OP_DZ)*temp79b)

#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = temp79a &
          + 2.*ri3_79*(1. - 3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2) &
          *bzt79(:,OP_1)*(e(:,OP_DZ)*f(:,OP_DRP) - e(:,OP_DR)*f(:,OP_DZP))

     temp79b = f(:,OP_DZZP) - f(:,OP_DRRP)
     if(itor.eq.1) temp79b = temp79b + 2.*ri_79*f(:,OP_DRP)
     temp79c = f(:,OP_DRZP)
     if(itor.eq.1) temp79c = temp79c -    ri_79*f(:,OP_DZP)

     temp79a = temp79a - 2.*e(:,OP_1)* &
          (6.*ri5_79*b2i79(:,OP_1)*bzt79(:,OP_1) &
          *(pst79(:,OP_DR)*pst79(:,OP_DZ)*temp79b &
          + (pst79(:,OP_DR)**2 - pst79(:,OP_DZ)**2)*temp79c) &
          +ri4_79*(1.+3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2) &
          *(pst79(:,OP_DZ)*f(:,OP_DZPP) + pst79(:,OP_DR)*f(:,OP_DRPP)))
#endif

     g2chi = -0.25*int3(temp79a,pit79(:,OP_1),b2i79(:,OP_1))          
  end select
     
  return
end function g2chi


! g3u
! ===
vectype function g3u(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  if(surface_int) then
     g3u = 0.
     return
  end if

  select case(ivform)
  case(0)
     if(itor.eq.1) then
        temp79b = ri_79*e(:,OP_DR)
        temp79d = ri_79*f(:,OP_DR)
        temp79e = ri_79*f(:,OP_DZ)
     else
        temp79b = 0.
        temp79d = 0.
        temp79e = 0.
     endif

     temp79a = &
          0.5*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1) * &
           ((e(:,OP_DZZ) - e(:,OP_DRR))* &
            (f(:,OP_DZZ) - f(:,OP_DRR) + temp79d) &
           +4.*e(:,OP_DRZ)*(f(:,OP_DRZ) - 0.5*temp79e)) &
          +1.5*ri4_79*b2i79(:,OP_1)**2*bzt79(:,OP_1) * &
           ((f(:,OP_DZZ) - f(:,OP_DRR) + temp79d) * &
            (pst79(:,OP_DR)**2*(e(:,OP_DZZ) - temp79b) &
            -pst79(:,OP_DZ)**2*(e(:,OP_DRR) - temp79b)) &
           +2.*(e(:,OP_DRZ)*pst79(:,OP_DR)**2*(f(:,OP_DRZ) +    temp79e) &
               +e(:,OP_DRZ)*pst79(:,OP_DZ)**2*(f(:,OP_DRZ) - 2.*temp79e) &
               -pst79(:,OP_DR)*pst79(:,OP_DZ)* &
                ( f(:,OP_DRZ)           *(e(:,OP_DRR) - temp79b    ) &
                +(f(:,OP_DRZ) - temp79e)*(e(:,OP_DZZ) - temp79b    ) &
                + temp79e               *(e(:,OP_DRR) - e(:,OP_DZZ)))))
     
     g3u = int2(pit79(:,OP_1),temp79a)
  case(1)
          
     temp79b = e(:,OP_DRR)
     if(itor.eq.1) temp79b = temp79b - 2.*ri_79*e(:,OP_DR)
     temp79c = e(:,OP_DRZ)
     if(itor.eq.1) temp79c = temp79c -    ri_79*e(:,OP_DZ)
     temp79d = f(:,OP_DZZ) - f(:,OP_DRR)
     if(itor.eq.1) temp79d = temp79d -    ri_79*f(:,OP_DR)
     temp79e = f(:,OP_DRZ)
     if(itor.eq.1) temp79e = temp79e +    ri_79*f(:,OP_DZ)

     temp79a = bzt79(:,OP_1) * &
          ((1. + 3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DR)**2) &
          *(2.*temp79c*temp79e     + e(:,OP_DZZ)*temp79d) & 
          +(1. + 3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DZ)**2) &
          *(2.*temp79c*f(:,OP_DRZ) - temp79b    *temp79d) &
          -6.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DR)*pst79(:,OP_DZ) &
          *(e(:,OP_DZZ)*f(:,OP_DRZ) + temp79b*temp79e))

     if(itor.eq.1) then
        temp79a = temp79a + 3.*ri3_79*b2i79(:,OP_1)*bzt79(:,OP_1)*e(:,OP_DR)* &
             ((pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2)*temp79d &
             +2.*pst79(:,OP_DZ)*pst79(:,OP_DR)*(f(:,OP_DRZ) + temp79e))
     endif

#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = temp79a + 3.*ri3_79*b2i79(:,OP_1)* &
          (pst79(:,OP_DZ)*f(:,OP_DZP) + pst79(:,OP_DR)*f(:,OP_DRP)) * &
          ((pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2)*temp79c &
          +pst79(:,OP_DZ)*pst79(:,OP_DR)*(temp79b - e(:,OP_DZZ))) &
          + ri_79*(1. - 3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2)*temp79c* &
          (pst79(:,OP_DZ)*f(:,OP_DZP) - pst79(:,OP_DR)*f(:,OP_DRP)) &
          + 3.*ri3_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2*(e(:,OP_DZZ)-temp79b)* &
          (pst79(:,OP_DR)*f(:,OP_DZP) + pst79(:,OP_DZ)*f(:,OP_DRP)) &
          - ri_79*(1. + 3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2)* &
          (pst79(:,OP_DZ)*f(:,OP_DRP)*e(:,OP_DZZ) &
          -pst79(:,OP_DR)*f(:,OP_DZP)*temp79b)

     if(itor.eq.1) then
        temp79a = temp79a + ri2_79*e(:,OP_DR)* &
             (1. + 3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2)* &
             (pst79(:,OP_DZ)*f(:,OP_DRP) - pst79(:,OP_DR)*f(:,OP_DZP))
     endif

     temp79d = f(:,OP_DZZP) - f(:,OP_DRRP)
     if(itor.eq.1) temp79d = temp79d -    ri_79*f(:,OP_DRP)
     temp79e = f(:,OP_DRZP)
     if(itor.eq.1) temp79e = temp79e +    ri_79*f(:,OP_DZP)

     temp79a = temp79a - ri_79* &
          (temp79d &
          *((1.-3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DZ)**2) &
          * e(:,OP_DR)*pst79(:,OP_DZ) &
          + (1.-3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DR)**2) &
          * e(:,OP_DZ)*pst79(:,OP_DR)) &
          +(1.+3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DR)**2) &
          *e(:,OP_DZ)*pst79(:,OP_DZ)*temp79e &
          -(1.+3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DZ)**2) &
          *e(:,OP_DR)*pst79(:,OP_DR)*f(:,OP_DRZP) &
          +3.*ri2_79*b2i79(:,OP_1) &
          *((bzt79(:,OP_1)**2 - pst79(:,OP_DZ)**2) &
          * e(:,OP_DR)*pst79(:,OP_DR)*temp79e &
          - (bzt79(:,OP_1)**2 - pst79(:,OP_DR)**2) &
          * e(:,OP_DZ)*pst79(:,OP_DZ)*f(:,OP_DRZP)) &
          -ri_79*(1.-3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2) &
          *bzt79(:,OP_1)*(e(:,OP_DZ)*f(:,OP_DZPP) + e(:,OP_DR)*f(:,OP_DRPP)))
#endif
     
     g3u = 0.5*int4(ri2_79,temp79a,pit79(:,OP_1),b2i79(:,OP_1))

  end select

  return
end function g3u


! g3v
! ===
vectype function g3v(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  if(surface_int) then
     g3v = 0.
     return
  end if

  select case(ivform)
  case(0)
     if(itor.eq.1) then
        temp79b = ri_79*e(:,OP_DR)
        temp79c = ri_79*f(:,OP_1)
     else
        temp79b = 0.
        temp79c = 0.
     endif
     
     temp79a = &
          -0.5*ri2_79*b2i79(:,OP_1) * &
          (pst79(:,OP_DZ)*(e(:,OP_DZZ) - temp79b)* f(:,OP_DZ) &
          +pst79(:,OP_DR)*(e(:,OP_DRR) - temp79b)*(f(:,OP_DR) - 2.*temp79c) &
          +e(:,OP_DRZ)*(pst79(:,OP_DZ)*(f(:,OP_DR) - 2.*temp79c) &
                       +pst79(:,OP_DR)* f(:,OP_DZ))) &
          -1.5*ri4_79*b2i79(:,OP_1)**2 * &
          (f(:,OP_DZ)*pst79(:,OP_DR) - f(:,OP_DR)*pst79(:,OP_DZ) &
          +2.*temp79c*pst79(:,OP_DZ)) * &
          (pst79(:,OP_DZ)*pst79(:,OP_DR)*(e(:,OP_DZZ) - e(:,OP_DRR)) &
          -e(:,OP_DRZ)*(pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2)) &
          +1.5*ri4_79*b2i79(:,OP_1)**2*bzt79(:,OP_1)**2 * &
          (e(:,OP_DRZ)*(pst79(:,OP_DZ)*(f(:,OP_DR)-2.*temp79c) &
                       +pst79(:,OP_DR)* f(:,OP_DZ)) &
          -(e(:,OP_DZZ) - temp79b)*(f(:,OP_DR)-2.*temp79c)*pst79(:,OP_DR) &
          -(e(:,OP_DRR) - temp79b)* f(:,OP_DZ)            *pst79(:,OP_DZ))

     g3v = int2(pit79(:,OP_1),temp79a)     
  case(1)
        
     temp79b = e(:,OP_DRZ)
     if(itor.eq.1) temp79b = temp79b - ri_79*e(:,OP_DZ)
     temp79c = e(:,OP_DRR)
     if(itor.eq.1) temp79c = temp79c - 2.*ri_79*e(:,OP_DR)
     temp79d = ri2_79*(f(:,OP_DZ)*pst79(:,OP_DR) - f(:,OP_DR)*pst79(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
     temp79d = temp79d + 2.*ri3_79*bzt79(:,OP_1)*f(:,OP_DP)
#endif

     temp79a = &
          (1. - 3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2)* &
          (pst79(:,OP_DR)*f(:,OP_DR)*e(:,OP_DZZ) &
          +pst79(:,OP_DZ)*f(:,OP_DZ)*temp79c &
          -(pst79(:,OP_DZ)*f(:,OP_DR) + pst79(:,OP_DR)*f(:,OP_DZ))*temp79b) &
          + 3.*b2i79(:,OP_1)*temp79d* &
          ((pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2)*temp79b &
          +pst79(:,OP_DZ)*pst79(:,OP_DR)*(temp79c - e(:,OP_DZZ))) &
          - (e(:,OP_DZZ) + temp79c) * &
          (pst79(:,OP_DZ)*f(:,OP_DZ) + pst79(:,OP_DR)*f(:,OP_DR))

     if(itor.eq.1) then
        temp79a = temp79a + ri_79*e(:,OP_DR)* &
             (1. + 3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2)* &
             (pst79(:,OP_DZ)*f(:,OP_DZ) + pst79(:,OP_DR)*f(:,OP_DR))
     endif

#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = temp79a - ri_79* &
          ((1.-3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2) &
          *bzt79(:,OP_1)*(e(:,OP_DZ)*f(:,OP_DRP) - e(:,OP_DR)*f(:,OP_DZP)) &
          -(1.+3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2)*ri_79 &
          *f(:,OP_DPP)*(e(:,OP_DZ)*pst79(:,OP_DZ) + e(:,OP_DR)*pst79(:,OP_DR)))
#endif

     g3v = 0.5*int4(ri2_79,temp79a,pit79(:,OP_1),b2i79(:,OP_1))
  end select

  return
end function g3v


! g3chi
! =====
vectype function g3chi(e,f)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f

  if(surface_int) then
     g3chi = 0.
     return
  end if

  select case(ivform)
  case(0)
     if(itor.eq.1) then
        temp79b = ri_79*e(:,OP_DR)
        temp79c = ri_79*f(:,OP_DR)
     else
        temp79b = 0.
        temp79c = 0.
     endif

     temp79a = &
          -ri_79*b2i79(:,OP_1)*bzt79(:,OP_1)* &
           ((e(:,OP_DZZ) - e(:,OP_DRR))*f(:,OP_DRZ) &
           -(f(:,OP_DZZ) - f(:,OP_DRR))*e(:,OP_DRZ)) &
          -3.*ri3_79*b2i79(:,OP_1)**2*bzt79(:,OP_1) * &
           (pst79(:,OP_DR)**2 * &
            ((e(:,OP_DZZ) - temp79b)*f(:,OP_DRZ) &
            -(f(:,OP_DZZ) - temp79c)*e(:,OP_DRZ)) &
           +pst79(:,OP_DZ)**2 * &
            ((f(:,OP_DRR) - temp79c)*e(:,OP_DRZ) &
            -(e(:,OP_DRR) - temp79b)*f(:,OP_DRZ)) &
           +pst79(:,OP_DZ)*pst79(:,OP_DR) * &
            ((e(:,OP_DRR) - temp79b)*(f(:,OP_DZZ) - temp79c) &
            -(e(:,OP_DZZ) - temp79b)*(f(:,OP_DRR) - temp79c)))

     g3chi = int2(pit79(:,OP_1),temp79a)     
  case(1)

     temp79b = e(:,OP_DRR)
     if(itor.eq.1) temp79b = temp79b - 2.*ri_79*e(:,OP_DR)
     temp79c = e(:,OP_DRZ)
     if(itor.eq.1) temp79c = temp79c -    ri_79*e(:,OP_DZ)
     temp79d = f(:,OP_DZZ)
     if(itor.eq.1) temp79d = temp79d -    ri_79*f(:,OP_DR)
     temp79e = f(:,OP_DRZ)
     if(itor.eq.1) temp79e = temp79e -    ri_79*f(:,OP_DZ)
     temp79f = f(:,OP_DRR)
     if(itor.eq.1) temp79f = temp79f - 3.*ri_79*f(:,OP_DR)

     temp79a = 2.*ri3_79*bzt79(:,OP_1) * &
          ((1. + 3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DR)**2) &
          *(temp79c*temp79d - e(:,OP_DZZ)*temp79e) &
          +(1. + 3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DZ)**2) &
          *(temp79b*temp79e - temp79c*temp79f)) &
          + 6.*ri5_79*b2i79(:,OP_1)*bzt79(:,OP_1) * &
          pst79(:,OP_DR)*pst79(:,OP_DZ)*(e(:,OP_DZZ)*temp79f - temp79b*temp79d)

     if(itor.eq.1) then
        temp79a = temp79a + 6.*ri6_79*b2i79(:,OP_1)*bzt79(:,OP_1)*e(:,OP_DR)* &
             (pst79(:,OP_DZ)*pst79(:,OP_DR)*(temp79d-temp79f) &
             -(pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2)*temp79e)
     endif

#if defined(USE3D) || defined(USECOMPLEX)
     temp79a = temp79a &
          + ri4_79*(1. - 3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2)* &
          (pst79(:,OP_DR)*f(:,OP_DRP)*e(:,OP_DZZ) &
          +pst79(:,OP_DZ)*f(:,OP_DZP)*temp79b &
          -(pst79(:,OP_DZ)*f(:,OP_DRP) + pst79(:,OP_DR)*f(:,OP_DZP))*temp79c) &
          - ri4_79*(temp79b + e(:,OP_DZZ)) * &
          (pst79(:,OP_DZ)*f(:,OP_DZP) + pst79(:,OP_DR)*f(:,OP_DRP)) &
          + 3.*ri6_79*b2i79(:,OP_1) * &
          (pst79(:,OP_DZ)*pst79(:,OP_DR)*(temp79b - e(:,OP_DZZ)) &
          +(pst79(:,OP_DZ)**2 - pst79(:,OP_DR)**2)*temp79c)* &
          (f(:,OP_DZP)*pst79(:,OP_DR) - f(:,OP_DRP)*pst79(:,OP_DZ))

     if(itor.eq.1) then
        temp79a = temp79a + ri5_79*e(:,OP_DR) * &
             (1. + 3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2) * &
             (pst79(:,OP_DZ)*f(:,OP_DZP) + pst79(:,OP_DR)*f(:,OP_DRP))
     end if

     temp79d = f(:,OP_DZZP)
     if(itor.eq.1) temp79d = temp79d -    ri_79*f(:,OP_DRP)
     temp79e = f(:,OP_DRZP)
     if(itor.eq.1) temp79e = temp79e -    ri_79*f(:,OP_DZP)
     temp79f = f(:,OP_DRRP)
     if(itor.eq.1) temp79f = temp79f - 3.*ri_79*f(:,OP_DRP)

     temp79a = temp79a - ri_79* &
          (ri3_79* &
          ((1.+3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DZ)**2) &
          * e(:,OP_DR)*pst79(:,OP_DR)*temp79f &
          +(1.+3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DR)**2) &
          * e(:,OP_DZ)*pst79(:,OP_DZ)*temp79d) &
          -2.*ri3_79*temp79e* &
          ((1.-3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DZ)**2) &
          * e(:,OP_DR)*pst79(:,OP_DZ) &
          +(1.-3.*ri2_79*b2i79(:,OP_1)*pst79(:,OP_DR)**2) &
          * e(:,OP_DZ)*pst79(:,OP_DR)) &
          +3.*ri5_79*b2i79(:,OP_1)* &
          ((bzt79(:,OP_1)**2 - pst79(:,OP_DZ)**2) &
          * e(:,OP_DR)*pst79(:,OP_DR)*temp79d &
          +(bzt79(:,OP_1)**2 - pst79(:,OP_DR)**2) &
          * e(:,OP_DZ)*pst79(:,OP_DZ)*temp79f) &
          +(1.-3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2)*ri4_79*bzt79(:,OP_1) &
          *(e(:,OP_DZ)*f(:,OP_DRPP) - e(:,OP_DR)*f(:,OP_DZPP)))
          
#endif

     g3chi = 0.5*int4(ri2_79,temp79a,pit79(:,OP_1),b2i79(:,OP_1))
  end select

  return
end function g3chi





! ==============================================================
! Ohmic heating terms
! ==============================================================

! qpsipsieta
! ==========
vectype function qpsipsieta(e)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e
  vectype :: temp

  if(hypf.eq.0. .or. jadv.eq.1 .or. surface_int) then
     qpsipsieta = 0.
     return
  end if

  temp79a = ri2_79* &
       (jt79(:,OP_DZ)**2 + jt79(:,OP_DR)**2)
  if(itor.eq.1) then
     temp79a = temp79a - 2.*jt79(:,OP_1)* &
          (ri3_79*jt79(:,OP_DR) - ri4_79*jt79(:,OP_1))
  endif

  if(ihypeta.eq.1) then
     temp = hypf*int4(e(:,OP_1),eta79(:,OP_1),temp79a,sz79(:,OP_1))
  else
     temp = hypf*int3(e(:,OP_1),temp79a,sz79(:,OP_1))
  endif


  temp79a = jt79(:,OP_DZ)*ni79(:,OP_DZ) + jt79(:,OP_DR)*ni79(:,OP_DR)
  if(itor.eq.1) then
     temp79a = temp79a - ri_79*jt79(:,OP_1)*ni79(:,OP_DR)
  endif
  temp79a = temp79a * ri2_79*nt79(:,OP_1)*jt79(:,OP_1)
     
  if(ihypeta.eq.1) then
     temp = temp + hypf*int4(e(:,OP_1),eta79(:,OP_1),temp79a,sz79(:,OP_1))
  else
     temp = temp + hypf*int3(e(:,OP_1),temp79a,sz79(:,OP_1))
  endif

  qpsipsieta = temp
  return
end function qpsipsieta

! qbbeta
! ======
vectype function qbbeta(e)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e
  vectype :: temp

  if(hypi.eq.0. .or. surface_int) then
     qbbeta = 0.
     return
  end if

  temp79a = ri2_79*(bzt79(:,OP_GS)*bzt79(:,OP_GS) &
       + 2.*(bzt79(:,OP_DRZ)**2 - bzt79(:,OP_DRR)*bzt79(:,OP_DZZ)))
  
  if(itor.eq.1) then 
     temp79a = temp79a + 2.* &
          (ri3_79*bzt79(:,OP_DZZ)*bzt79(:,OP_DR) &
          -ri3_79*bzt79(:,OP_DRZ)*bzt79(:,OP_DZ) &
          +ri4_79*bzt79(:,OP_DZ )*bzt79(:,OP_DZ))
  endif

  if(ihypeta.eq.1) then
     temp = hypi*int4(e(:,OP_1),eta79(:,OP_1),temp79a,sz79(:,OP_1))
  else
     temp = hypi*int3(e(:,OP_1),temp79a,sz79(:,OP_1))
  endif

  temp79a = &
       bzt79(:,OP_DZ)*bzt79(:,OP_DZZ)*ni79(:,OP_DZ) &
       + bzt79(:,OP_DR)*bzt79(:,OP_DRZ)*ni79(:,OP_DZ) &
       + bzt79(:,OP_DZ)*bzt79(:,OP_DRZ)*ni79(:,OP_DR) &
       + bzt79(:,OP_DR)*bzt79(:,OP_DRR)*ni79(:,OP_DR)
  if(itor.eq.1) then
     temp79a = temp79a - ri_79*ni79(:,OP_DR)* &
          (bzt79(:,OP_DZ)**2 + bzt79(:,OP_DR)**2)
  endif
  temp79a = temp79a * ri2_79*nt79(:,OP_1)
     
  if(ihypeta.eq.1) then
     temp = temp + hypi*int4(e(:,OP_1),eta79(:,OP_1),temp79a,sz79(:,OP_1))
  else
     temp = temp + hypi*int3(e(:,OP_1),temp79a,sz79(:,OP_1))
  endif

  qbbeta = temp
  return
end function qbbeta


! ==============================================================
! Viscous heating terms
! ==============================================================

! quumu
! =====
vectype function quumu(e,f,g,h,i)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = -int5(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_GS),h(:,OP_1)) &
             + 4.* int5(ri2_79,e(:,OP_1),f(:,OP_DRR),g(:,OP_DZZ),h(:,OP_1)) &
             - 4.* int5(ri2_79,e(:,OP_1),f(:,OP_DRZ),g(:,OP_DRZ),h(:,OP_1))

        if(itor.eq.1) then
           temp = temp &
                + 4.*int5(ri3_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1)) &
                - 4.*int5(ri3_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
                + 4.*int5(ri4_79,e(:,OP_1 ),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1))
        endif
        
        if(hypc.ne.0.) then
           temp79a = h(:,OP_1)*i(:,OP_1)
           temp = temp - &
                (int5(ri2_79,e(:,OP_1),vot79(:,OP_DZ),vot79(:,OP_DZ),temp79a) &
                +int5(ri2_79,e(:,OP_1),vot79(:,OP_DR),vot79(:,OP_DR),temp79a))
        endif
     end if
     
  case(1)
        
     ! Not yet implemented
     temp = 0.

  end select

  quumu = temp
  return
end function quumu


! quchimu
! =======
vectype function quchimu(e,f,g,h,i,j)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i,j

  quchimu = 0.
  return
end function quchimu


! qvvmu
! =====
vectype function qvvmu(e,f,g,h,i)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = - &
             (int5(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
             +int5(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1)))

        if(itor.eq.1) then
           temp = temp &
                + 4.*int5(ri3_79,e(:,OP_1),f(:,OP_DR),g(:,OP_1),h(:,OP_1)) &
                - 4.*int5(ri4_79,e(:,OP_1),f(:,OP_1 ),g(:,OP_1),h(:,OP_1))
        endif
  
        if(hypv.ne.0.) then
           temp79a = h(:,OP_1)*i(:,OP_1)
           temp = temp - int5(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_GS),temp79a)
        endif
     end if

  case(1)
     if(surface_int) then
        temp = 0.
     else
        temp = - &
             (int5(r2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1)) &
             +int5(r2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1)))
        
        if(hypv.ne.0.) then
           temp79a = h(:,OP_1)*i(:,OP_1)
           temp = temp - int5(r2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_GS),temp79a)
        endif
     end if
        
  end select

  qvvmu = temp
  return
end function qvvmu


! qchichimu
! =========
vectype function qchichimu(e,f,g,h,i)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,i
  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = 4.*int5(ri_79,e(:,OP_1),h(:,OP_1),f(:,OP_DZZ),g(:,OP_DRZ)) &
             - 4.*int5(ri_79,e(:,OP_1),h(:,OP_1),f(:,OP_DRZ),g(:,OP_DZZ)) &
             + 4.*int5(ri_79,e(:,OP_1),h(:,OP_1),f(:,OP_DRZ),g(:,OP_DRR)) &
             - 4.*int5(ri_79,e(:,OP_1),h(:,OP_1),f(:,OP_DRR),g(:,OP_DRZ))

        if(itor.eq.1) then
           temp = temp &
                + 4.*int5(ri2_79,e(:,OP_DR),h(:,OP_1),f(:,OP_DZ),g(:,OP_DR)) &
                - 4.*int5(ri2_79,e(:,OP_DZ),h(:,OP_1),f(:,OP_DR),g(:,OP_DR))
        endif
  
        if(hypc.ne.0.) then
           temp79a = h(:,OP_1)*i(:,OP_1)
           temp = temp - 2.* &
                (int4(e(:,OP_1),cot79(:,OP_DZ),cot79(:,OP_DZ),temp79a) &
                +int4(e(:,OP_1),cot79(:,OP_DR),cot79(:,OP_DR),temp79a))
        endif
     end if

  case(1)
     ! Not yet implemented
     temp =0.
     
  end select

  qchichimu = temp
  return
end function qchichimu


!======================================================================
! ENERGY
!======================================================================

#ifdef USECOMPLEX
#define CONJUGATE(x) conjg(x)
#else
#define CONJUGATE(x) x
#endif

! Poloidal magnetic
! -----------------
real function energy_mp()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(linear.eq.1) then
     temp = .5* &
          (int3(ri2_79,ps179(:,OP_DZ),CONJUGATE(ps179(:,OP_DZ))) &
          +int3(ri2_79,ps179(:,OP_DR),CONJUGATE(ps179(:,OP_DR))))
#ifdef USECOMPLEX
     temp = temp + .5* &
          (int2(bf179(:,OP_DZP),CONJUGATE(bf179(:,OP_DZP))) &
          +int2(bf179(:,OP_DRP),CONJUGATE(bf179(:,OP_DRP))) &
          +int3(ri_79,ps179(:,OP_DZ),CONJUGATE(bf179(:,OP_DRP))) &
          -int3(ri_79,ps179(:,OP_DR),CONJUGATE(bf179(:,OP_DZP))) &
          +int3(ri_79,CONJUGATE(ps179(:,OP_DZ)),bf179(:,OP_DRP)) &
          -int3(ri_79,CONJUGATE(ps179(:,OP_DR)),bf179(:,OP_DZP)))
#endif
  else
!    nonlinear:   subtract off equilibrium piece
     temp = .5* &
          (int3(ri2_79,pst79(:,OP_DZ),pst79(:,OP_DZ)) &
          +int3(ri2_79,pst79(:,OP_DR),pst79(:,OP_DR))) &
          - .5* &
          (int3(ri2_79,ps079(:,OP_DZ),ps079(:,OP_DZ)) &
          +int3(ri2_79,ps079(:,OP_DR),ps079(:,OP_DR)))
#if defined(USE3D)
     temp = temp   &
          + .5* &
          (int2(bft79(:,OP_DZP),bft79(:,OP_DZP)) &
          +int2(bft79(:,OP_DRP),bft79(:,OP_DRP)) &
          +2.*int3(ri_79,pst79(:,OP_DZ),bft79(:,OP_DRP)) &
          -2.*int3(ri_79,pst79(:,OP_DR),bft79(:,OP_DRP)) )
#endif

#ifdef USECOMPLEX
     temp = temp + .5* &
          (int2(bft79(:,OP_DZP),CONJUGATE(bft79(:,OP_DZP))) &
          +int2(bft79(:,OP_DRP),CONJUGATE(bft79(:,OP_DRP))) &
          +int3(ri_79,pst79(:,OP_DZ),CONJUGATE(bft79(:,OP_DRP))) &
          -int3(ri_79,pst79(:,OP_DR),CONJUGATE(bft79(:,OP_DRP))) &
          +int3(ri_79,CONJUGATE(pst79(:,OP_DZ)),bft79(:,OP_DRP)) &
          -int3(ri_79,CONJUGATE(pst79(:,OP_DR)),bft79(:,OP_DZP)))
#endif
  endif

  energy_mp = temp
  return
end function energy_mp


! Toroidal magnetic
! -----------------
real function energy_mt()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(linear.eq.1) then
     temp = .5*int3(ri2_79,bz179(:,OP_1),CONJUGATE(bz179(:,OP_1)))
  else
!....nonlinear:  subtract off equilibrium piece
     temp = .5*int3(ri2_79,bzt79(:,OP_1),bzt79(:,OP_1))   &
          - .5*int3(ri2_79,bz079(:,OP_1),bz079(:,OP_1))
  endif

  energy_mt = temp
  return
end function energy_mt


! Pressure
! --------
real function energy_p()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(gam.le.1.) then 
     temp = 0.
  else
     if(linear.eq.1) then
        temp = int1(p179) / (gam - 1.)
     else
!.......nonlinear: subtract off equilibrium piece
        temp = (int1(pt79) - int1(p079))/ (gam - 1.)
     endif
  endif

  energy_p = temp
  return
end function energy_p



! Poloidal kinetic
! ----------------
real function energy_kp()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  select case(ivform)
  case(0)
     if(linear.eq.1) then
        temp = .5* &
             (int4(ri2_79,ph179(:,OP_DZ),CONJUGATE(ph179(:,OP_DZ)),n079(:,OP_1)) &
             +int4(ri2_79,ph179(:,OP_DR),CONJUGATE(ph179(:,OP_DR)),n079(:,OP_1)))
     else
        temp = .5* &
             (int4(ri2_79,pht79(:,OP_DZ),CONJUGATE(pht79(:,OP_DZ)),nt79(:,OP_1)) &
             +int4(ri2_79,pht79(:,OP_DR),CONJUGATE(pht79(:,OP_DR)),nt79(:,OP_1)))
     endif

  case(1)
     if(linear.eq.1) then
        temp = .5* &
             (int4(r2_79,ph179(:,OP_DZ),CONJUGATE(ph179(:,OP_DZ)),n079(:,OP_1)) &
             +int4(r2_79,ph179(:,OP_DR),CONJUGATE(ph179(:,OP_DR)),n079(:,OP_1)))
     else
        temp = .5* &
             (int4(r2_79,pht79(:,OP_DZ),pht79(:,OP_DZ),nt79(:,OP_1)) &
             +int4(r2_79,pht79(:,OP_DR),pht79(:,OP_DR),nt79(:,OP_1)))
     endif
  end select

  energy_kp = temp
  return
end function energy_kp


! Toroidal kinetic
! ----------------
real function energy_kt()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  select case(ivform)
  case(0)
     if(linear.eq.1) then
        temp = .5*int4(ri2_79,vz179(:,OP_1),CONJUGATE(vz179(:,OP_1)),n079(:,OP_1)) 
     else
        temp = .5*int4(ri2_79,vzt79(:,OP_1),CONJUGATE(vzt79(:,OP_1)),nt79(:,OP_1))
     endif

  case(1)
     if(linear.eq.1) then
        temp = .5*int4(r2_79,vz179(:,OP_1),CONJUGATE(vz179(:,OP_1)),n079(:,OP_1))
     else
        temp = .5*int4(r2_79,vzt79(:,OP_1),vzt79(:,OP_1),nt79(:,OP_1))
     endif
  end select

  energy_kt = temp
  return
end function energy_kt


! Compressional kinetic
! ---------------------
real function energy_k3()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  select case(ivform)
  case(0)
     if(linear.eq.1) then
        temp = .5* &
          (int3(ch179(:,OP_DZ),CONJUGATE(ch179(:,OP_DZ)),n079(:,OP_1)) &
          +int3(ch179(:,OP_DR),CONJUGATE(ch179(:,OP_DR)),n079(:,OP_1)) &
          +int4(ri_79,ch179(:,OP_DZ),CONJUGATE(ph179(:,OP_DR)),n079(:,OP_1)) &
          -int4(ri_79,ch179(:,OP_DR),CONJUGATE(ph179(:,OP_DZ)),n079(:,OP_1)) &
          +int4(ri_79,CONJUGATE(ch179(:,OP_DZ)),ph179(:,OP_DR),n079(:,OP_1)) &
          -int4(ri_79,CONJUGATE(ch179(:,OP_DR)),ph179(:,OP_DZ),n079(:,OP_1)))
     else
        temp = .5* &
          (int3(cht79(:,OP_DZ),CONJUGATE(cht79(:,OP_DZ)),nt79(:,OP_1)) &
          +int3(cht79(:,OP_DR),CONJUGATE(cht79(:,OP_DR)),nt79(:,OP_1)) &
          +int4(ri_79,cht79(:,OP_DZ),CONJUGATE(pht79(:,OP_DR)),nt79(:,OP_1)) &
          -int4(ri_79,cht79(:,OP_DR),CONJUGATE(pht79(:,OP_DZ)),nt79(:,OP_1)) &
          +int4(ri_79,CONJUGATE(cht79(:,OP_DZ)),pht79(:,OP_DR),nt79(:,OP_1)) &
          -int4(ri_79,CONJUGATE(cht79(:,OP_DR)),pht79(:,OP_DZ),nt79(:,OP_1)))
     endif
     
  case(1)
     if(linear.eq.1) then
        temp = .5* &
          (int4(ri4_79,ch179(:,OP_DZ),CONJUGATE(ch179(:,OP_DZ)),n079(:,OP_1)) &
          +int4(ri4_79,ch179(:,OP_DR),CONJUGATE(ch179(:,OP_DR)),n079(:,OP_1)) &
          +int4(ri_79,ch179(:,OP_DZ),CONJUGATE(ph179(:,OP_DR)),n079(:,OP_1)) &
          -int4(ri_79,ch179(:,OP_DR),CONJUGATE(ph179(:,OP_DZ)),n079(:,OP_1)) &
          +int4(ri_79,CONJUGATE(ch179(:,OP_DZ)),ph179(:,OP_DR),n079(:,OP_1)) &
          -int4(ri_79,CONJUGATE(ch179(:,OP_DR)),ph179(:,OP_DZ),n079(:,OP_1)))
     else
        temp = .5* &
          (int4(ri4_79,cht79(:,OP_DZ),cht79(:,OP_DZ),nt79(:,OP_1)) &
          +int4(ri4_79,cht79(:,OP_DR),cht79(:,OP_DR),nt79(:,OP_1)) &
          +int4(ri_79,cht79(:,OP_DZ),pht79(:,OP_DR),nt79(:,OP_1)) &
          -int4(ri_79,cht79(:,OP_DR),pht79(:,OP_DZ),nt79(:,OP_1)) &
          +int4(ri_79,cht79(:,OP_DZ),pht79(:,OP_DR),nt79(:,OP_1)) &
          -int4(ri_79,cht79(:,OP_DR),pht79(:,OP_DZ),nt79(:,OP_1)))
     endif
  end select

  energy_k3 = temp
  return
end function energy_k3


! Poloidal resistive
! ------------------
real function energy_mpd()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(linear.eq.1) then
     temp = -int4(ri2_79,ps179(:,OP_GS),CONJUGATE(ps179(:,OP_GS)),eta79(:,OP_1))
#ifdef USECOMPLEX
     temp = temp - &
          (int4(ri4_79,ps179(:,OP_DZP),CONJUGATE(ps179(:,OP_DZP)),eta79(:,OP_1)) &
          +int4(ri4_79,ps179(:,OP_DRP),CONJUGATE(ps179(:,OP_DRP)),eta79(:,OP_1)))
#endif
  else
     temp = -int4(ri2_79,pst79(:,OP_GS),CONJUGATE(pst79(:,OP_GS)),eta79(:,OP_1))
  endif

  energy_mpd = temp
  return
end function energy_mpd


! Toroidal resistive
! ------------------
real function energy_mtd()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(linear.eq.1) then
     temp = - &
          (int4(ri2_79,bz179(:,OP_DZ),CONJUGATE(bz179(:,OP_DZ)),eta79(:,OP_1))&
          +int4(ri2_79,bz179(:,OP_DR),CONJUGATE(bz179(:,OP_DR)),eta79(:,OP_1)))
  else
     temp = - &
          (int4(ri2_79,bzt79(:,OP_DZ),CONJUGATE(bzt79(:,OP_DZ)),eta79(:,OP_1))&
          +int4(ri2_79,bzt79(:,OP_DR),CONJUGATE(bzt79(:,OP_DR)),eta79(:,OP_1)))
  end if

  energy_mtd = temp
  return
end function energy_mtd


! Poloidal viscous
! ----------------
real function energy_kpd()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(linear.eq.1) then
     temp = -int4(ri2_79,ph179(:,OP_GS),CONJUGATE(ph179(:,OP_GS)),vis79(:,OP_1))
#ifdef USECOMPLEX
     temp = temp - &
          (int4(ri4_79,ph179(:,OP_DZP),CONJUGATE(ph179(:,OP_DZP)),vis79(:,OP_1)) &
          +int4(ri4_79,ph179(:,OP_DRP),CONJUGATE(ph179(:,OP_DRP)),vis79(:,OP_1)))
#endif
  else
     temp = -int4(ri2_79,pht79(:,OP_GS),CONJUGATE(pht79(:,OP_GS)),vis79(:,OP_1))
  endif

  energy_kpd = temp
  return
end function energy_kpd


! Toroidal viscous
! ----------------
real function energy_ktd()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  select case(ivform)
  case(0)
     if(linear.eq.1) then
        temp = - &
             (int4(ri2_79,vz179(:,OP_DZ),CONJUGATE(vz179(:,OP_DZ)),vis79(:,OP_1)) &
             +int4(ri2_79,vz179(:,OP_DR),CONJUGATE(vz179(:,OP_DR)),vis79(:,OP_1)))
     else
        temp = - &
             (int4(ri2_79,vzt79(:,OP_DZ),CONJUGATE(vzt79(:,OP_DZ)),vis79(:,OP_1)) &
             +int4(ri2_79,vzt79(:,OP_DR),CONJUGATE(vzt79(:,OP_DR)),vis79(:,OP_1)))
     endif
  case(1)
     if(linear.eq.1) then
        temp = - &
             (int4(r2_79,vz179(:,OP_DZ),CONJUGATE(vz179(:,OP_DZ)),vis79(:,OP_1)) &
             +int4(r2_79,vz179(:,OP_DR),CONJUGATE(vz179(:,OP_DR)),vis79(:,OP_1)))
     else
        temp = - &
             (int4(r2_79,vzt79(:,OP_DZ),CONJUGATE(vzt79(:,OP_DZ)),vis79(:,OP_1)) &
             +int4(r2_79,vzt79(:,OP_DR),CONJUGATE(vzt79(:,OP_DR)),vis79(:,OP_1)))
     endif
  end select

  energy_ktd = temp
  return
end function energy_ktd

! Compressional viscous
! ---------------------
real function energy_k3d()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(linear.eq.1) then
     temp = - 2.*int3(ch179(:,OP_LP),CONJUGATE(ch179(:,OP_LP)),vic79(:,OP_1))
  else
     temp = - 2.*int3(cht79(:,OP_LP),CONJUGATE(cht79(:,OP_LP)),vic79(:,OP_1))
  endif

  energy_k3d = temp
  return
end function energy_k3d


! Poloidal hyper-viscous
! ----------------------
real function energy_kph()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  temp = - hypc* &
       (int5(ri2_79,vot79(:,OP_DZ),CONJUGATE(vot79(:,OP_DZ)),vis79(:,OP_1),sz79(:,OP_1)) &
       +int5(ri2_79,vot79(:,OP_DR),CONJUGATE(vot79(:,OP_DR)),vis79(:,OP_1),sz79(:,OP_1)))

  energy_kph = temp
  return
end function energy_kph


! Toroidal hyper-viscous
! ----------------------
real function energy_kth()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  select case(ivform)
  case(0)
     temp = -hypv*int5(ri2_79,vzt79(:,OP_GS),CONJUGATE(vzt79(:,OP_GS)),vis79(:,OP_1),sz79(:,OP_1))
  case(1)
     temp = -hypv*int5(r2_79,vzt79(:,OP_GS),CONJUGATE(vzt79(:,OP_GS)),vis79(:,OP_1),sz79(:,OP_1))
  end select

  energy_kth = temp
  return
end function energy_kth

! Compressional hyper-viscous
! ---------------------------
real function energy_k3h()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  temp = -2.*hypc* &
       (int4(cot79(:,OP_DZ),CONJUGATE(cot79(:,OP_DZ)),vic79(:,OP_1),sz79(:,OP_1)) &
       +int4(cot79(:,OP_DR),CONJUGATE(cot79(:,OP_DR)),vic79(:,OP_1),sz79(:,OP_1)))

  energy_k3h = temp
  return
end function energy_k3h



!======================================================================
! FLUXES
!======================================================================


! Pressure convection
! -------------------
real function flux_pressure()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if((numvar.lt.3 .and. ipres.eq.0) .or. gam.eq.1. .or. gam.eq.0.) then
     flux_pressure = 0.
     return
  endif

  select case(ivform)
  case(0)
     temp = 0.5* &
          (int4(ri_79,pt79(:,OP_1),norm79(:,2),pht79(:,OP_DR)) &
          -int4(ri_79,pt79(:,OP_1),norm79(:,1),pht79(:,OP_DZ)) &
          +int3(pt79(:,OP_1),norm79(:,1),cht79(:,OP_DR)) &
          +int3(pt79(:,OP_1),norm79(:,2),cht79(:,OP_DZ)))

  case(1)
     temp = 0.5* &
          (int4(r_79,pt79(:,OP_1),norm79(:,2),pht79(:,OP_DR)) &
          -int4(r_79,pt79(:,OP_1),norm79(:,1),pht79(:,OP_DZ)) &
          +int4(ri2_79,pt79(:,OP_1),norm79(:,1),cht79(:,OP_DR)) &
          +int4(ri2_79,pt79(:,OP_1),norm79(:,2),cht79(:,OP_DZ)))

  end select

  if(dbf .ne. 0.) then
     temp = temp + dbf* &
          (int5(ri_79,pet79(:,OP_1),ni79(:,OP_1),norm79(:,1),bzt79(:,OP_DZ)) &
          -int5(ri_79,pet79(:,OP_1),ni79(:,OP_1),norm79(:,2),bzt79(:,OP_DR)))
  endif

  flux_pressure = gam*real(temp)/(gam-1.)

  return
end function flux_pressure


! Kinetic Energy Convection
! -------------------------
real function flux_ke()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  select case(ivform)
  case(0)
     temp79a = ri2_79*(pht79(:,OP_DZ)**2 + pht79(:,OP_DR)**2) &
          +    ri2_79*vzt79(:,OP_1)**2 &
          +           (cht79(:,OP_DZ)**2 + cht79(:,OP_DR)**2) &
          + 2.*ri_79* &
          (cht79(:,OP_DZ)*pht79(:,OP_DR)-cht79(:,OP_DR)*pht79(:,OP_DZ))
     
     temp = 0.5* &
          (int5(ri_79,nt79(:,OP_1),norm79(:,2),pht79(:,OP_DR),temp79a) &
          -int5(ri_79,nt79(:,OP_1),norm79(:,1),pht79(:,OP_DZ),temp79a) &
          +int4(nt79(:,OP_1),norm79(:,1),cht79(:,OP_DR),temp79a) &
          +int4(nt79(:,OP_1),norm79(:,2),cht79(:,OP_DZ),temp79a))
  case(1)
     temp79a = r2_79*(pht79(:,OP_DZ)**2 + pht79(:,OP_DR)**2) &
          +    r2_79*vzt79(:,OP_1)**2 &
          +   ri4_79*(cht79(:,OP_DZ)**2 + cht79(:,OP_DR)**2) &
          + 2.*ri_79* &
          (cht79(:,OP_DZ)*pht79(:,OP_DR)-cht79(:,OP_DR)*pht79(:,OP_DZ))

     temp = 0.5* &
          (int5(r_79,nt79(:,OP_1),norm79(:,2),pht79(:,OP_DR),temp79a) &
          -int5(r_79,nt79(:,OP_1),norm79(:,1),pht79(:,OP_DZ),temp79a) &
          +int5(ri2_79,nt79(:,OP_1),norm79(:,1),cht79(:,OP_DR),temp79a) &
          +int5(ri2_79,nt79(:,OP_1),norm79(:,2),cht79(:,OP_DZ),temp79a))
  end select

  flux_ke = real(temp)
  return
end function flux_ke


! Poynting flux
! -------------
real function flux_poynting()

  use math
  use basic
  use m3dc1_nint

  implicit none
  vectype :: temp

  temp = -vloop/twopi * &
       (int3(ri2_79,norm79(:,1),pst79(:,OP_DR)) &
       +int3(ri2_79,norm79(:,1),pst79(:,OP_DZ)))

  flux_poynting = real(temp)
  return
end function flux_poynting


! Heat flux
! ---------
real function flux_heat()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(numvar.lt.3 .and. ipres.eq.0) then
     flux_heat = 0.
     return
  endif

  temp = int4(kap79(:,1),norm79(:,1),pt79(:,OP_DR),ni79(:,OP_1)) &
       + int4(kap79(:,1),norm79(:,2),pt79(:,OP_DZ),ni79(:,OP_1)) &
       + int4(kap79(:,1),norm79(:,1),pt79(:,OP_1),ni79(:,OP_DR)) &
       + int4(kap79(:,1),norm79(:,2),pt79(:,OP_1),ni79(:,OP_DZ))

  if(kappar.ne.0.) then
     temp79a = ni79(:,OP_1)* &
          (pt79(:,OP_DZ)*pst79(:,OP_DR) - pt79(:,OP_DR)*pst79(:,OP_DZ)) &
          +    pt79(:,OP_1)* &
          (ni79(:,OP_DZ)*pst79(:,OP_DR) - ni79(:,OP_DR)*pst79(:,OP_DZ))
     temp79b = norm79(:,1)*pst79(:,OP_DZ) - norm79(:,2)*pst79(:,OP_DR)
     temp = temp &
          + int5(ri2_79,kar79(:,OP_1),b2i79(:,OP_1),temp79a,temp79b)
  endif

  flux_heat = real(temp)
  return
end function flux_heat


! Grav_pot
! --------
real function grav_pot()

  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(gravr.eq.0. .and. gravz.eq.0.) then
     grav_pot = 0.
     return
  endif

  temp = gravr*int3(ri3_79,pht79(:,OP_DZ),nt79(:,OP_1)) &
       - gravz*int3( ri_79,pht79(:,OP_DR),nt79(:,OP_1))
     
  if(numvar.ge.3) then
     temp = -gravr*int3(ri2_79,cht79(:,OP_DR),nt79(:,OP_1)) &
          -gravz*int2(       cht79(:,OP_DZ),nt79(:,OP_1))
  endif

  grav_pot = real(temp)
  return
end function grav_pot


!======================================================================
! Toroidal (angular) momentum
!======================================================================

! torque_em
! ~~~~~~~~~
vectype function torque_em()
  use m3dc1_nint

  implicit none

  vectype :: temp

  temp = int4(ri_79,bzt79(:,OP_1),norm79(:,2),pst79(:,OP_DR)) &
       - int4(ri_79,bzt79(:,OP_1),norm79(:,1),pst79(:,OP_DZ))
  
  torque_em = temp
end function torque_em

! torque_sol
! ~~~~~~~~~~
vectype function torque_sol()
  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  select case(ivform)
  case(0)
     temp = int5(ri_79,nt79(:,OP_1),vzt79(:,OP_1),norm79(:,1),pht79(:,OP_DZ)) &
          - int5(ri_79,nt79(:,OP_1),vzt79(:,OP_1),norm79(:,2),pht79(:,OP_DR))
  case(1)
     temp = int5(r3_79,nt79(:,OP_1),vzt79(:,OP_1),norm79(:,1),pht79(:,OP_DZ)) &
          - int5(r3_79,nt79(:,OP_1),vzt79(:,OP_1),norm79(:,2),pht79(:,OP_DR))
  end select
  
  torque_sol = temp
end function torque_sol

! torque_com
! ~~~~~~~~~~
vectype function torque_com()
  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(numvar.lt.3) then
     torque_com = 0.
     return
  end if

  select case(ivform)
  case(0)
     temp = &
          - int4(nt79(:,OP_1),vzt79(:,OP_1),norm79(:,1),cht79(:,OP_DR)) &
          - int4(nt79(:,OP_1),vzt79(:,OP_1),norm79(:,2),cht79(:,OP_DZ))
  case(1)
     temp = &
          - int4(nt79(:,OP_1),vzt79(:,OP_1),norm79(:,1),cht79(:,OP_DR)) &
          - int4(nt79(:,OP_1),vzt79(:,OP_1),norm79(:,2),cht79(:,OP_DZ))     
  end select

  torque_com = temp
end function torque_com


! torque_visc
! ~~~~~~~~~~~
vectype function torque_visc()
  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp

  select case(ivform)
  case(0)
     temp = int3(vis79(:,OP_1),norm79(:,1),vzt79(:,OP_DR)) &
          + int3(vis79(:,OP_1),norm79(:,2),vzt79(:,OP_DZ))
  case(1)
     temp = int4(r2_79,vis79(:,OP_1),norm79(:,1),vzt79(:,OP_DR)) &
          + int4(r2_79,vis79(:,OP_1),norm79(:,2),vzt79(:,OP_DZ))
  end select

  torque_visc = temp
end function torque_visc


! torque_parvisc
! ~~~~~~~~~~~~~~
vectype function torque_parvisc()
  use basic
  use m3dc1_nint

  implicit none

  vectype :: temp
  
  if(amupar.eq.0.) then
     torque_parvisc = 0.
     return
  endif
       
  call PVS1(pht79,temp79a)

  if(numvar.ge.2) then
     call PVS2(vzt79,temp79b)
     temp79a = temp79a + temp79b
  endif
  
  if(numvar.ge.3) then
     call PVS3(cht79,temp79c)
     temp79a = temp79a + temp79c
  endif

  temp79d = 3.*vip79(:,OP_1)*bzt79(:,OP_1)*b2i79(:,OP_1)*temp79a
  temp = int4(ri_79,temp79d,norm79(:,2),pst79(:,OP_DR)) &
       - int4(ri_79,temp79d,norm79(:,1),pst79(:,OP_DZ))

  torque_parvisc = temp
end function torque_parvisc



! torque_gyro
! ~~~~~~~~~~~
vectype function torque_gyro()
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype :: temp

  if(gyro.eq.0 .or. numvar.lt.2) then
     torque_gyro = 0.
     return
  endif

  temp79a = 0.25*dbf*b2i79(:,OP_1)
  temp79b = temp79a*(1. - 3.*ri2_79*b2i79(:,OP_1)*bzt79(:,OP_1)**2)

  select case(ivform)
  case(0)
     temp = 0.

  case(1)
     temp79c = norm79(:,1)*pst79(:,OP_DZ) + norm79(:,2)*pst79(:,OP_DR)
     temp79d = norm79(:,1)*pst79(:,OP_DR) - norm79(:,2)*pst79(:,OP_DZ)
     temp79e = 3.*temp79a*b2i79(:,OP_1)* &
          (norm79(:,1)*pst79(:,OP_DZ) - norm79(:,2)*pst79(:,OP_DR))
     temp79f = pht79(:,OP_DZZ) - pht79(:,OP_DRR)
     if(itor.eq.1) then
        temp79f = temp79f - ri_79*pht79(:,OP_DR)
     endif

     ! U contribution
     temp = int4(r_79,temp79b,temp79c,temp79f) &
        +2.*int4(r_79,temp79b,temp79d,pht79(:,OP_DRZ)) &
        +   int5(ri_79,temp79e,pst79(:,OP_DZ),pst79(:,OP_DZ),temp79f) &
        -   int5(ri_79,temp79e,pst79(:,OP_DR),pst79(:,OP_DR),temp79f) &
        +4.*int5(ri_79,temp79e,pst79(:,OP_DR),pst79(:,OP_DZ),pht79(:,OP_DRZ))

     if(itor.eq.1) then
        temp = temp &
           + 2.*int4(temp79b,norm79(:,1),pst79(:,OP_DR),pht79(:,OP_DZ)) &
           - 2.*int4(temp79a,pht79(:,OP_DZ),norm79(:,1),pst79(:,OP_DR)) &
           - 2.*int4(temp79a,pht79(:,OP_DZ),norm79(:,2),pst79(:,OP_DZ)) &
           + 2.*int5(ri2_79,temp79e,pst79(:,OP_DR),pst79(:,OP_DZ),pht79(:,OP_DZ))
     endif
     
     ! omega contribution
     temp = temp + 2.* &
          (int5(r_79,temp79b,bzt79(:,OP_1),norm79(:,1),vzt79(:,OP_DZ)) &
          -int5(r_79,temp79b,bzt79(:,OP_1),norm79(:,2),vzt79(:,OP_DR)))

     temp79f = cht79(:,OP_DZZ) - cht79(:,OP_DRR)
     if(itor.eq.1) then
        temp79f = temp79f + 2.*ri_79*cht79(:,OP_DR)
     endif

     ! chi constribution
     temp = temp &
          +    int4(ri2_79,temp79b,temp79d,temp79f) &
          - 2.*int4(ri2_79,temp79b,temp79c,cht79(:,OP_DRZ)) &
          +    int5(ri2_79,temp79b,cht79(:,OP_GS),norm79(:,1),pst79(:,OP_DR)) &
          +    int5(ri2_79,temp79b,cht79(:,OP_GS),norm79(:,2),pst79(:,OP_DZ)) &
          - 2.*int4(ri2_79,cht79(:,OP_GS),norm79(:,1),pst79(:,OP_DR)) &
          - 2.*int4(ri2_79,cht79(:,OP_GS),norm79(:,2),pst79(:,OP_DZ)) &
          + 2.*int5(ri4_79,temp79e,pst79(:,OP_DR),pst79(:,OP_DZ),temp79f) &
          - 2.*int5(ri4_79,temp79e,pst79(:,OP_DZ),pst79(:,OP_DZ),cht79(:,OP_DRZ)) &
          + 2.*int5(ri4_79,temp79e,pst79(:,OP_DR),pst79(:,OP_DR),cht79(:,OP_DRZ))

     if(itor.eq.1) then
        temp = temp &
             + 2.*int4(ri3_79,temp79b,temp79c,cht79(:,OP_DZ)) &
             + 3.*int5(ri3_79,temp79b,cht79(:,OP_DR),norm79(:,1),pst79(:,OP_DR)) &
             + 3.*int5(ri3_79,temp79b,cht79(:,OP_DR),norm79(:,2),pst79(:,OP_DZ)) &
             - 6.*int4(ri3_79,cht79(:,OP_DR),norm79(:,1),pst79(:,OP_DR)) &
             - 6.*int4(ri3_79,cht79(:,OP_DR),norm79(:,2),pst79(:,OP_DZ)) &
             + 2.*int5(ri5_79,temp79e,pst79(:,OP_DZ),pst79(:,OP_DZ),cht79(:,OP_DZ)) &
             - 2.*int5(ri5_79,temp79e,pst79(:,OP_DR),pst79(:,OP_DR),cht79(:,OP_DZ))
     endif
  end select

  torque_gyro = 0.
end function torque_gyro


! torque_denm
! ~~~~~~~~~~~
vectype function torque_denm()
  use basic
  use m3dc1_nint

  implicit none

  torque_denm = 0.

end function torque_denm





vectype function tepsipsikappar(e,f,g,h,j,k)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,j,k
  vectype :: temp

  if(gam.le.1.) then
     tepsipsikappar = 0.
     return
  end if

  if(surface_int) then
     temp79a = ri2_79*k(:,OP_1)*j(:,OP_1)*e(:,OP_1)* &
          (norm79(:,1)*f(:,OP_DZ) - norm79(:,2)*f(:,OP_DR))
     temp = int3(temp79a,g(:,OP_DZ),h(:,OP_DR)) &
          - int3(temp79a,g(:,OP_DR),h(:,OP_DZ))

  else
     temp79a = k(:,OP_1)*ri2_79* &
          (e(:,OP_DZ)*f(:,OP_DR) - e(:,OP_DR)*f(:,OP_DZ))*j(:,OP_1)

     temp = int3(temp79a,g(:,OP_DZ),h(:,OP_DR)) &
          - int3(temp79a,g(:,OP_DR),h(:,OP_DZ))

  end if

  tepsipsikappar = (gam - 1.) * temp
  return
end function tepsipsikappar
vectype function tepsibkappar(e,f,g,h,j,k)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,j,k
  vectype :: temp

  if(gam.le.1.) then
     tepsibkappar = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     temp79a = -ri3_79*k(:,OP_1)*j(:,OP_1)*e(:,OP_1)*g(:,OP_1)* &
          (norm79(:,1)*f(:,OP_DZ) - norm79(:,2)*f(:,OP_DR))

     temp = int2(temp79a,h(:,OP_DP))
  else
     temp79a = -k(:,OP_1)*ri3_79*g(:,OP_1)* &
          (e(:,OP_DZ)*f(:,OP_DR) - e(:,OP_DR)*f(:,OP_DZ))*j(:,OP_1)  

     temp79b = f(:,OP_DR)*(h(:,OP_DZ) ) &
          -    f(:,OP_DZ)*(h(:,OP_DR) )

     temp79d = temp79b*g(:,OP_1 )*j(:,OP_1 )*k(:,OP_1)

     temp = int2(temp79a,h(:,OP_DP)) &
          - int3(ri3_79,e(:,OP_DP),temp79d)
  end if
#else
  temp = 0.
#endif

  tepsibkappar = (gam - 1.) * temp
  return
end function tepsibkappar
vectype function tebbkappar(e,f,g,h,j,k)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,j,k
  vectype :: temp

  if(gam.le.1.) then
     tebbkappar = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     temp = 0.
  else
     temp79a = h(:,OP_DP)

     temp79c = f(:,OP_1)*g(:,OP_1 )*temp79a*j(:,OP_1 )*k(:,OP_1 )

     temp = -int3(ri4_79,e(:,OP_DP),temp79c)
  end if
#else
  temp = 0.
#endif

  tebbkappar = (gam - 1.) * temp
  return
end function tebbkappar
vectype function tepsifkappar(e,f,g,h,j,k)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,j,k
  vectype :: temp

  if(gam.le.1.) then
     tepsifkappar = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     temp79a = k(:,OP_1)*ri_79*e(:,OP_1)* &
          (norm79(:,2)*f(:,OP_DR) - norm79(:,1)*f(:,OP_DZ))*j(:,OP_1)
     temp79b = -k(:,OP_1)*ri_79*e(:,OP_1)* &
          (norm79(:,2)*g(:,OP_DZP) + norm79(:,1)*g(:,OP_DRP))*j(:,OP_1)

     temp = int3(temp79a,g(:,OP_DZP),h(:,OP_DZ)) &
          + int3(temp79a,g(:,OP_DRP),h(:,OP_DR)) &
          + int3(temp79b,f(:,OP_DR ),h(:,OP_DZ)) &
          - int3(temp79b,f(:,OP_DZ ),h(:,OP_DR))
  else
     temp79a = k(:,OP_1)*ri_79* &
          (e(:,OP_DZ)*f(:,OP_DR) - e(:,OP_DR)*f(:,OP_DZ))*j(:,OP_1)
     temp79b = k(:,OP_1)*ri_79* &
          (e(:,OP_DZ)*g(:,OP_DZP) + e(:,OP_DR)*g(:,OP_DRP))*j(:,OP_1)

     temp = int3(temp79a,g(:,OP_DZP),h(:,OP_DZ)) &
          + int3(temp79a,g(:,OP_DRP),h(:,OP_DR)) &
          + int3(temp79b,f(:,OP_DR ),h(:,OP_DZ)) &
          - int3(temp79b,f(:,OP_DZ ),h(:,OP_DR))
  end if
#else
  temp = 0.
#endif

  tepsifkappar = (gam - 1.) * temp
  return
end function tepsifkappar
vectype function tebfkappar(e,f,g,h,j,k)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,j,k
  vectype :: temp

  if(gam.le.1.) then
     tebfkappar = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     temp79a = -ri2_79*k(:,OP_1)*j(:,OP_1)*e(:,OP_1)*f(:,OP_1)* &
          (norm79(:,1)*g(:,OP_DRP) + norm79(:,2)*g(:,OP_DZP))

     temp = int2(temp79a,h(:,OP_DP))
  else
     temp79a = k(:,OP_1)*ri2_79*f(:,OP_1)* &
          (e(:,OP_DZ)*g(:,OP_DZP) + e(:,OP_DR)*g(:,OP_DRP))*j(:,OP_1)

     temp79b = g(:,OP_DZP)*(h(:,OP_DZ) )&
          +    g(:,OP_DRP)*(h(:,OP_DR) )

     temp79d = temp79b*f(:,OP_1 )*j(:,OP_1 )*k(:,OP_1 )

     temp = int2(temp79a,h(:,OP_DP)) &
          + int3(ri2_79,e(:,OP_DP),temp79d)
  end if
#else
  temp = 0.
#endif

  tebfkappar = (gam - 1.) * temp
  return
end function tebfkappar
vectype function teffkappar(e,f,g,h,j,k)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h,j,k
  vectype :: temp

  if(gam.le.1.) then
     teffkappar = 0.
     return
  end if

#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     temp79a =  k(:,OP_1)*e(:,OP_1)* &
          (norm79(:,2)*f(:,OP_DZP) + norm79(:,1)*f(:,OP_DRP))*j(:,OP_1)

     temp = int3(temp79a,g(:,OP_DZP),h(:,OP_DZ)) &
          + int3(temp79a,g(:,OP_DRP),h(:,OP_DR))
  else
     temp79a = - k(:,OP_1)*                                            &
          (e(:,OP_DZ)*f(:,OP_DZP) + e(:,OP_DR)*f(:,OP_DRP))*j(:,OP_1)

     temp = int3(temp79a,g(:,OP_DZP),h(:,OP_DZ)) &
          + int3(temp79a,g(:,OP_DRP),h(:,OP_DR))
  end if
#else
  temp = 0.
#endif

  teffkappar = (gam - 1.) * temp
  return
end function teffkappar
vectype function b3peeta(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = int3(e(:,OP_1),f(:,OP_1),g(:,OP_1))
  end if

  b3peeta = temp
  return
end function b3peeta

vectype function q1ppsi(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp


  if(surface_int) then
     temp = 0.
  else
     temp = int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1)) &
          - int5(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1))
  end if

  q1ppsi =   temp
  return
end function q1ppsi

vectype function q1pb(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp


#if defined(USE3D) || defined(USECOMPLEX)
  if(surface_int) then
     temp = 0.
  else
     temp = int5(ri2_79,e(:,OP_1),f(:,OP_DP),g(:,OP_1),h(:,OP_1))
  end if
#else
  temp = 0.
#endif

  q1pb =   temp
  return
end function q1pb

vectype function q1pf(e,f,g,h)
  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp


#if defined(USE3D) || defined(USECOMPLEX)

  if(surface_int) then
     temp = 0.
  else

     temp = - int4(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZP),h(:,OP_1)) &
            - int4(e(:,OP_1),f(:,OP_DR),g(:,OP_DRP),h(:,OP_1))
  end if
#else
  temp = 0.
#endif

  q1pf =   temp
  return
end function q1pf
vectype function t3tneta(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
  vectype :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = int3(e(:,OP_1),f(:,OP_1),g(:,OP_1))
  end if

  t3tneta = temp
  return
end function t3tneta
vectype function t3tn(e,f,g)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g
  vectype :: temp

  if(surface_int) then
     temp = 0.
  else
     temp = int2(e(:,OP_1),f(:,OP_1))
  end if

  t3tn = temp
  return
end function t3tn
vectype function t3tnu(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

  vectype :: temp

  select case(ivform)
  case(0)
!==> This needs to be redone
     if(surface_int) then
        temp = 0.
     else
        temp = int4(ri_79,e(:,OP_1),f(:,OP_DR),h(:,OP_DZ)) &
             - int4(ri_79,e(:,OP_1),f(:,OP_DZ),h(:,OP_DR))
     end if

  case(1)
     if(surface_int) then
!!$        temp = int5(r_79,e(:,OP_1),f(:,OP_1),norm79(:,1),h(:,OP_DZ)) &
!!$             - int5(r_79,e(:,OP_1),f(:,OP_1),norm79(:,2),h(:,OP_DR))
        temp = 0.
     else
        temp = int5(r_79,e(:,OP_1),f(:,OP_DR),g(:,OP_1),h(:,OP_DZ)) &
             - int5(r_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_1),h(:,OP_DR))

        if(itor.eq.1) then
           temp = temp + &
                2.*(gam-1.)*int4(e(:,OP_1),f(:,OP_1),g(:,OP_1),h(:,OP_DZ))
        endif
     end if
  end select

  t3tnu = temp

  return
end function t3tnu
vectype function t3tnv(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

  vectype :: temp

#if defined(USE3D) || defined(USECOMPLEX)
  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
!==> This needs to be checked
        temp = - int4(ri2_79,e(:,OP_1),f(:,OP_DP),h(:,OP_1)) &
             - gam*int4(ri2_79,e(:,OP_1),f(:,OP_1),h(:,OP_DP))
     endif
  case(1)
     if(surface_int) then
        temp = 0.
     else
        temp = - int4(e(:,OP_1),f(:,OP_DP),g(:,OP_1),h(:,OP_1)) &
             - (gam-1.)*int4(e(:,OP_1),f(:,OP_1),g(:,OP_1),h(:,OP_DP))
     endif
  end select
#else
  temp = 0.
#endif

  t3tnv = temp

  return
end function t3tnv
vectype function t3tnchi(e,f,g,h)

  use basic
  use m3dc1_nint

  implicit none

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h

  vectype :: temp

  select case(ivform)
  case(0)
     if(surface_int) then
        temp = 0.
     else
        temp = gam* &
             (int3(e(:,OP_DZ),f(:,OP_1),h(:,OP_DZ))  &
             +int3(e(:,OP_DR),f(:,OP_1),h(:,OP_DR))) &
             +(gam-1.)* &
             (int3(e(:,OP_1),f(:,OP_DZ),h(:,OP_DZ))  &
             +int3(e(:,OP_1),f(:,OP_DR),h(:,OP_DR)))
     end if
     
  case(1)
     if(surface_int) then
!!$        temp = &
!!$             - int5(ri2_79,e(:,OP_1),f(:,OP_1),norm79(:,1),h(:,OP_DR)) &
!!$             - int5(ri2_79,e(:,OP_1),f(:,OP_1),norm79(:,2),h(:,OP_DZ))
        temp = 0.
     else
        temp79a = -(f(:,OP_DR)*h(:,OP_DR) + f(:,OP_DZ)*h(:,OP_DZ))  &
                  - (gam-1.)*f(:,OP_1)*h(:,OP_GS)
     
        temp = int4(ri2_79,temp79a,e(:,OP_1),g(:,OP_1))
     endif
  end select

  t3tnchi = temp

  return
end function t3tnchi
end module metricterms_new
