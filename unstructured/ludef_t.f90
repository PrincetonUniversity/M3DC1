!======================================================================
! Vorticity Equation
!======================================================================
subroutine vorticity_lin(trial, lin, ssterm, ddterm, q_bf, advfield)
  
  use basic
  use arrays
  use m3dc1_nint
  use metricterms_new

  implicit none

  vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: trial, lin 
  vectype, dimension(num_fields), intent(out) :: ssterm, ddterm
  vectype, intent(out) :: q_bf
  integer, intent(in) :: advfield   ! if advfield = 1, eliminate rrterm by
                                    ! using analytic form of advanced field
  vectype :: temp
  real :: ththm

  select case(imp_mod)
  case(0)
     ththm = (1.-thimp*bdf)*thimp
  case(1)
     ththm = -bdf*thimp**2
  case(2)
     ththm = -bdf*thimp
  end select

  ssterm = 0.
  ddterm = 0.
  q_bf = 0.

  if(istatic.eq.1) then
     if(.not.surface_int) then
        temp = int2(trial,lin)
        ssterm(u_g) = temp
        ddterm(u_g) = temp
     endif
     return
  endif

  ! NUMVAR = 1
  ! ~~~~~~~~~~
  temp = v1un(trial,lin,nt79)
  ssterm(u_g) = ssterm(u_g) + temp
  ddterm(u_g) = ddterm(u_g) + temp*bdf
  
  temp = v1umu(trial,lin,vis79,vic79) &
       + v1us (trial,lin,sig79)
  ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
  ddterm(u_g) = ddterm(u_g) + (1.-thimp*bdf)*dt*temp

  if(linear.eq.0) then 
     temp = v1uun(trial,lin,ph179,nt79) &
          + v1uun(trial,ph179,lin,nt79)
     ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
     ddterm(u_g) = ddterm(u_g) + (.5-thimp*bdf)*dt*temp

     temp = v1uvn(trial,lin,vz179,nt79)
     ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
     ddterm(u_g) = ddterm(u_g) + (.5-thimp*bdf)*dt*temp

     temp = v1uchin(trial,lin,ch179,nt79)       
     ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
     ddterm(u_g) = ddterm(u_g) + (.5-thimp*bdf)*dt*temp
  endif
  
  if(gyro.eq.1) then
     temp = g1u(trial,lin)*dbf
     ssterm(u_g) = ssterm(u_g) +     thimp     *dt*temp
     ddterm(u_g) = ddterm(u_g) - (1.-thimp*bdf)*dt*temp
  endif

  if(advfield.eq.1) then
     temp = v1upsipsi(trial,lin,pst79,pst79) &
          + v1upsib  (trial,lin,pst79,bzt79) &
          + v1ubb    (trial,lin,bzt79,bzt79) &
          + v1up     (trial,lin,pt79) &
          + v1ungrav (trial,lin,nt79)
     ssterm(u_g) = ssterm(u_g) - thimp*thimp*dt*dt*temp
     ddterm(u_g) = ddterm(u_g) +       ththm*dt*dt*temp

     ddterm(psi_g) = ddterm(psi_g) + dt* &
          (v1psipsi(trial,lin,pss79) &
          +v1psipsi(trial,pss79,lin) &
          +v1psib  (trial,lin,bzs79))

     if(numvar.ge.3 .or. ipres.eq.1) then
        ddterm(p_g) = ddterm(p_g) + dt*    &
             v1p(trial,lin)
     endif
  else
     if(linear.eq.0) then
        temp = v1psipsi(trial,lin,ps179) &
             + v1psipsi(trial,ps179,lin)
        ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
        ddterm(psi_g) = ddterm(psi_g) + (.5-thimp*bdf)*dt*temp

        temp = v1psib(trial,lin,bz179)
        ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
        if(numvar.eq.1) temp = 2.*temp
        ddterm(psi_g) = ddterm(psi_g) + (.5-thimp*bdf)*dt*temp
     endif

     if(numvar.ge.3 .or. ipres.eq.1) then
        temp = v1p(trial,lin)
        ssterm(p_g) = ssterm(p_g) -     thimp     *dt*temp
        ddterm(p_g) = ddterm(p_g) + (1.-thimp*bdf)*dt*temp
     endif
  endif
    
  if(eqsubtract.eq.1) then
     temp = v1uun  (trial,lin,ph079,nt79) &    
          + v1uun  (trial,ph079,lin,nt79) &
          + v1uvn  (trial,lin,vz079,nt79) &
          + v1uchin(trial,lin,ch079,nt79)
     ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
     ddterm(u_g) = ddterm(u_g) + (1.-thimp*bdf)*dt*temp
     
     if(advfield.eq.0) then
        temp = v1psipsi(trial,lin,ps079) &
             + v1psipsi(trial,ps079,lin) &
             + v1psib  (trial,lin,bz079)
        ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
        ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp
     endif
  endif
  
  
  ! NUMVAR = 2
  ! ~~~~~~~~~~
  if(numvar.ge.2) then
     temp = v1vmu(trial,lin,vis79,vic79)
     ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
     ddterm(vz_g) = ddterm(vz_g) + (1.-thimp*bdf)*dt*temp

     if(linear.eq.0) then        
        temp = v1vvn(trial,lin,vz179,nt79) &
             + v1vvn(trial,vz179,lin,nt79) &
             + v1uvn(trial,ph179,lin,nt79)
        ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
        ddterm(vz_g) = ddterm(vz_g) + (.5-thimp*bdf)*dt*temp

        temp = v1vchin(trial,lin,ch179,nt79)       
        ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
        ddterm(vz_g) = ddterm(vz_g) + (.5-thimp*bdf)*dt*temp
     endif
     
     if(gyro.eq.1) then
        temp = g1v(trial,lin)*dbf
        ssterm(vz_g) = ssterm(vz_g) +     thimp     *dt*temp
        ddterm(vz_g) = ddterm(vz_g) - (1.-thimp*bdf)*dt*temp
     endif

     if(advfield.eq.1) then
        temp = v1vpsipsi(trial,lin,pst79,pst79) &
             + v1vpsib  (trial,lin,pst79,bzt79) &
             + v1vp     (trial,lin,pt79)
        ssterm(vz_g) = ssterm(vz_g) - thimp*thimp*dt*dt*temp
        ddterm(vz_g) = ddterm(vz_g) +       ththm*dt*dt*temp

        ddterm(bz_g) = ddterm(bz_g) + dt* &
             (v1psib(trial,pss79,lin) &
             +v1bb  (trial,lin,bzs79) &
             +v1bb  (trial,bzs79,lin))
     else
        if(linear.eq.0) then 
           temp = v1psib(trial,ps179,lin) &
                + v1bb  (trial,lin,bz179) &
                + v1bb  (trial,bz179,lin)
           ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
           ddterm(bz_g) = ddterm(bz_g) + (.5-thimp*bdf)*dt*temp
        endif
     endif
          
     if(eqsubtract.eq.1) then
        temp = v1uvn  (trial,ph079,lin,nt79) & 
             + v1vvn  (trial,lin,vz079,nt79) &
             + v1vvn  (trial,vz079,lin,nt79) &
             + v1vchin(trial,lin,ch079,nt79)
        ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
        ddterm(vz_g) = ddterm(vz_g) + (1.-thimp*bdf)*dt*temp

        if(advfield.eq.0) then
           temp = v1psib(trial,ps079,lin) &
                + v1bb  (trial,lin,bz079) &
                + v1bb  (trial,bz079,lin)
           ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
           ddterm(bz_g) = ddterm(bz_g) + (1.-thimp*bdf)*dt*temp
        endif
     endif   !  on numvar .ne. 1

     if(i3d.eq.1) then
        temp = v1psif(trial,lin,bft79)
        ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
        ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

        temp = v1bf(trial,lin,bft79)
        ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
        ddterm(bz_g) = ddterm(bz_g) + (1.-thimp*bdf)*dt*temp

        if(eqsubtract.eq.1) then
           q_bf = q_bf + dt* &
                (v1psif(trial,ps079,lin) &
                +v1bf  (trial,bz079,lin))
        endif
     endif
  endif
  
  ! NUMVAR = 3
  ! ~~~~~~~~~~
  if(numvar.ge.3) then
     
     if(linear.eq.0) then 
        temp = v1uchin  (trial,ph179,lin,nt79) &
             + v1vchin  (trial,vz179,lin,nt79) &
             + v1chichin(trial,lin,ch179,nt79) &
             + v1chichin(trial,ch179,lin,nt79)
        ssterm(chi_g) = ssterm(chi_g) -     thimp     *dt*temp
        ddterm(chi_g) = ddterm(chi_g) + (.5-thimp*bdf)*dt*temp
     endif
     
     temp = v1chin(trial,lin,nt79)*chiiner
     ssterm(chi_g) = ssterm(chi_g) + temp
     ddterm(chi_g) = ddterm(chi_g) + temp*bdf
     
     temp = v1chimu(trial,lin,vis79,vic79) &
          + v1chis (trial,lin,sig79)
     ssterm(chi_g) = ssterm(chi_g) -     thimp     *dt*temp
     ddterm(chi_g) = ddterm(chi_g) + (1.-thimp*bdf)*dt*temp

     if(gyro.eq.1) then
        temp = g1chi(trial,lin)*dbf
        ssterm(chi_g) = ssterm(chi_g) +     thimp     *dt*temp
        ddterm(chi_g) = ddterm(chi_g) - (1.-thimp*bdf)*dt*temp
     endif

     if(advfield.eq.1) then
        temp = v1chipsipsi(trial,lin,pst79,pst79) &
             + v1chipsib  (trial,lin,pst79,bzt79) &
             + v1chibb    (trial,lin,bzt79,bzt79) &
             + v1chip     (trial,lin,pt79)        &
             + v1chingrav (trial,lin,nt79)
        ssterm(chi_g) = ssterm(chi_g) - thimp*thimp*dt*dt*temp
        ddterm(chi_g) = ddterm(chi_g) +       ththm*dt*dt*temp
     endif
     
     if(eqsubtract.eq.1) then        
        temp = v1uchin  (trial,ph079,lin,nt79) &
             + v1vchin  (trial,vz079,lin,nt79) &
             + v1chichin(trial,lin,ch079,nt79) &
             + v1chichin(trial,ch079,lin,nt79)
        ssterm(chi_g) = ssterm(chi_g) -     thimp     *dt*temp
        ddterm(chi_g) = ddterm(chi_g) + (1.-thimp*bdf)*dt*temp
     endif
  endif

  if(idens.eq.1) then
     if(advfield.eq.1) then
        ddterm(den_g) = ddterm(den_g) + dt* &
             v1ngrav(trial,lin)
     else
        temp = v1ngrav(trial,lin)
        ssterm(den_g) = ssterm(den_g) -     thimp     *dt*temp
        ddterm(den_g) = ddterm(den_g) + (1.-thimp*bdf)*dt*temp
     endif

     if(eqsubtract.eq.1) then
        ddterm(den_g) = ddterm(den_g) + dt* &
             (v1uun    (trial,ph079,ph079,lin) &
             +v1uvn    (trial,ph079,vz079,lin) &
             +v1vvn    (trial,vz079,vz079,lin) &
             +v1uchin  (trial,vz079,ch079,lin) &
             +v1vchin  (trial,vz079,ch079,lin) &
             +v1chichin(trial,ch079,ch079,lin))
     endif
  endif


  ! Parallel Viscosity
  if(amupar.ne.0.) then
     call PVV1(trial,temp79b)

     call PVS1(lin,temp79c)
     temp = int3(vip79(:,OP_1),temp79b,temp79c)
     ssterm(u_g) = ssterm(u_g) +     thimp     *dt*temp
     ddterm(u_g) = ddterm(u_g) - (1.-thimp*bdf)*dt*temp

     if(numvar.ge.2) then
        call PVS2(lin,temp79c)
        temp = int3(vip79(:,OP_1),temp79b,temp79c)
        ssterm(vz_g) = ssterm(vz_g) +     thimp     *dt*temp
        ddterm(vz_g) = ddterm(vz_g) - (1.-thimp*bdf)*dt*temp
     endif

     if(numvar.ge.3) then
        call PVS3(lin,temp79c)
        temp = int3(vip79(:,OP_1),temp79b,temp79c)
        ssterm(chi_g) = ssterm(chi_g) +     thimp     *dt*temp
        ddterm(chi_g) = ddterm(chi_g) - (1.-thimp*bdf)*dt*temp
     endif
  endif

end subroutine vorticity_lin 


subroutine vorticity_nolin(trial, r4term)

  use basic
  use m3dc1_nint
  use metricterms_new

  implicit none

  vectype, intent(in), dimension(MAX_PTS, OP_NUM)  :: trial
  vectype, intent(out) :: r4term

  r4term = 0.

  if(linear.eq.1) return

  ! density terms
  ! ~~~~~~~~~~~~~
  if(idens.eq.1 .and. eqsubtract.eq.1) then
     r4term = r4term + dt* &
          v1us     (trial,ph079,sig79)
             
     if(numvar.ge.3) then
        r4term = r4term + dt* &
             v1chis   (trial,ch079,sig79)
     endif
  endif
    
end subroutine vorticity_nolin


!======================================================================
! Axial Velocity Equation
!======================================================================
subroutine axial_vel_lin(trial, lin, ssterm, ddterm, q_bf, advfield)
  
  use basic
  use arrays
  use m3dc1_nint
  use metricterms_new

  implicit none

  vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: trial, lin 
  vectype, dimension(num_fields), intent(out) :: ssterm, ddterm
  vectype, intent(out) :: q_bf

  integer, intent(in) :: advfield
  vectype :: temp
  real :: ththm

  vectype, dimension(MAX_PTS, OP_NUM) :: hv
  hv = hypv*sz79

  select case(imp_mod)
  case(0)
     ththm = (1.-thimp*bdf)*thimp
  case(1)
     ththm = -bdf*thimp**2
  case(2)
     ththm = -bdf*thimp
  end select

  ssterm = 0.
  ddterm = 0.
  q_bf = 0.

  if(numvar.lt.2) return

  if(istatic.eq.1) then
     if(.not.surface_int) then
        temp = int2(trial,lin)
        ssterm(vz_g) = temp
        ddterm(vz_g) = temp
     endif
     return
  endif

  temp = v2umu(trial,lin,vis79,vic79)
  ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
  ddterm(u_g) = ddterm(u_g) + (1.-thimp*bdf)*dt*temp
         
  if(linear.eq.0) then 
     temp = v2vun(trial,vz179,lin,nt79)
     ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
     ddterm(u_g) = ddterm(u_g) + (.5-thimp*bdf)*dt*temp

     temp = v2vun(trial,lin,ph179,nt79) &
          + v2vvn(trial,lin,vz179,nt79) &
          + v2vvn(trial,vz179,lin,nt79)
     ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
     ddterm(vz_g) = ddterm(vz_g) + (.5-thimp*bdf)*dt*temp

     temp = v2vchin(trial,lin,ch179,nt79)
     ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
     ddterm(vz_g) = ddterm(vz_g) + (.5-thimp*bdf)*dt*temp
  endif
          
  temp = v2vn(trial,lin,nt79)
  ssterm(vz_g) = ssterm(vz_g) + temp
  ddterm(vz_g) = ddterm(vz_g) + temp*bdf
     
  temp = v2vmu  (trial,lin,vis79,vic79,hv) &
       + v2vs   (trial,lin,sig79)
  ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
  ddterm(vz_g) = ddterm(vz_g) + (1.-thimp*bdf)*dt*temp

  if(gyro.eq.1) then
     temp = g2u(trial,lin)*dbf
     ssterm(u_g) = ssterm(u_g) +     thimp     *dt*temp
     ddterm(u_g) = ddterm(u_g) - (1.-thimp*bdf)*dt*temp
     
     temp = g2v(trial,lin)*dbf
     ssterm(vz_g) = ssterm(vz_g) +     thimp     *dt*temp
     ddterm(vz_g) = ddterm(vz_g) - (1.-thimp*bdf)*dt*temp
  endif
          
  if(advfield.eq.1) then 
     temp = v2upsipsi(trial,lin,pst79,pst79) &
          + v2upsib  (trial,lin,pst79,bzt79) &
          + v2ubb    (trial,lin,bzt79,bzt79) &
          + v2up     (trial,lin,pt79)
     ssterm(u_g) = ssterm(u_g) - thimp*thimp*dt*dt*temp
     ddterm(u_g) = ddterm(u_g) +       ththm*dt*dt*temp

     temp = v2vpsipsi(trial,lin,pst79,pst79) &
          + v2vpsib  (trial,lin,pst79,bzt79) &
          + v2vp     (trial,lin,pt79)
     ssterm(vz_g) = ssterm(vz_g) - thimp*thimp*dt*dt*temp
     ddterm(vz_g) = ddterm(vz_g) +       ththm*dt*dt*temp

     ddterm(psi_g) = ddterm(psi_g) + dt* &
          (v2psipsi(trial,lin,pss79) & 
          +v2psipsi(trial,pss79,lin) & 
          +v2psib  (trial,lin,bzs79))
     
     ddterm(bz_g) = ddterm(bz_g) + dt* &
          (v2psib(trial,pss79,lin))

     if(numvar.ge.3 .or. ipres.eq.1) then
        ddterm(p_g) = ddterm(p_g) + dt* &
             v2p(trial,lin)
     endif
  else
     if(linear.eq.0) then
        temp = v2psipsi(trial,lin,ps179) &
             + v2psipsi(trial,ps179,lin) &
             + v2psib  (trial,lin,bz179)
        ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
        ddterm(psi_g) = ddterm(psi_g) + (.5-thimp*bdf)*dt*temp
        
        temp = v2psib(trial,ps179,lin)
        ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
        ddterm(bz_g) = ddterm(bz_g) + (.5-thimp*bdf)*dt*temp
     endif

     if(numvar.ge.3 .or. ipres.eq.1) then
        temp = v2p(trial,lin)
        ssterm(p_g) = ssterm(p_g) -     thimp     *dt*temp
        ddterm(p_g) = ddterm(p_g) + (1.-thimp*bdf)*dt*temp
     endif
  endif

  if(eqsubtract.eq.1) then
     
     temp = v2vun(trial,vz079,lin,nt79)
     ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
     ddterm(u_g) = ddterm(u_g) + (1.-thimp*bdf)*dt*temp
              
     temp = v2vun  (trial,lin,ph079,nt79) &
          + v2vvn  (trial,lin,vz079,nt79) &
          + v2vvn  (trial,vz079,lin,nt79) &
          + v2vchin(trial,lin,ch079,nt79)
     ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
     ddterm(vz_g) = ddterm(vz_g) + (1.-thimp*bdf)*dt*temp

     if(advfield.eq.0) then
        temp = v2psipsi(trial,lin,ps079) &
             + v2psipsi(trial,ps079,lin) &
             + v2psib  (trial,lin,bz079)
        ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
        ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

        temp = v2psib(trial,ps079,lin)
        ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
        ddterm(bz_g) = ddterm(bz_g) + (1.-thimp*bdf)*dt*temp
     end if
  endif

  if(i3d.eq.1) then
     temp = v2psif(trial,lin,bft79)
     ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

     temp = v2bf(trial,lin,bft79)
     ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
     ddterm(bz_g) = ddterm(bz_g) + (1.-thimp*bdf)*dt*temp

     if(eqsubtract.eq.1) then
        q_bf = q_bf + dt* &
             (v2psif(trial,ps079,lin) &
             +v2bf  (trial,bz079,lin) &
             +v2ff  (trial,bf079,lin) &
             +v2ff  (trial,lin,bf079))
     endif
  endif

        
  ! NUMVAR = 3
  ! ~~~~~~~~~~
  if(numvar.ge.3) then

     if(linear.eq.0) then 
        temp = v2vchin(trial,vz179,lin,nt79)
        ssterm(chi_g) = ssterm(chi_g) -     thimp     *dt*temp
        ddterm(chi_g) = ddterm(chi_g) + (.5-thimp*bdf)*dt*temp
     endif

     temp = v2chimu  (trial,lin,vis79,vic79)
     ssterm(chi_g) = ssterm(chi_g) -     thimp     *dt*temp
     ddterm(chi_g) = ddterm(chi_g) + (1.-thimp*bdf)*dt*temp
           
     if(gyro.eq.1) then
        temp = g2chi(trial,lin)*dbf
        ssterm(chi_g) = ssterm(chi_g) +     thimp     *dt*temp
        ddterm(chi_g) = ddterm(chi_g) - (1.-thimp*bdf)*dt*temp
     endif

     if(advfield.eq.1) then
        temp = v2chipsipsi(trial,lin,pst79,pst79) &
             + v2chipsib  (trial,lin,pst79,bzt79) &
             + v2chibb    (trial,lin,bzt79,bzt79) &
             + v2chip     (trial,lin,pt79)
        ssterm(chi_g) = ssterm(chi_g) - thimp*thimp*dt*dt*temp
        ddterm(chi_g) = ddterm(chi_g) +       ththm*dt*dt*temp
     end if
           
     if(eqsubtract.eq.1) then             
        temp = v2vchin(trial,vz079,lin,nt79)
        ssterm(chi_g) = ssterm(chi_g) -     thimp     *dt*temp
        ddterm(chi_g) = ddterm(chi_g) + (1.-thimp*bdf)*dt*temp
     endif
  endif

  if(idens.eq.1) then
     if(eqsubtract.eq.1) then
        ddterm(den_g) = ddterm(den_g) + dt* &
             (v2vun  (trial,vz079,ph079,lin) &
             +v2vvn  (trial,vz079,vz079,lin) &
             +v2vchin(trial,vz079,ch079,lin))
     endif
  endif


  ! Parallel Viscosity
  if(amupar.ne.0.) then
     call PVV2(trial,temp79b)

     call PVS1(lin,temp79c)
     temp = int3(vip79(:,OP_1),temp79b,temp79c)
     ssterm(u_g) = ssterm(u_g) +     thimp     *dt*temp
     ddterm(u_g) = ddterm(u_g) - (1.-thimp*bdf)*dt*temp

     if(numvar.ge.2) then
        call PVS2(lin,temp79c)
        temp = int3(vip79(:,OP_1),temp79b,temp79c)
        ssterm(vz_g) = ssterm(vz_g) +     thimp     *dt*temp
        ddterm(vz_g) = ddterm(vz_g) - (1.-thimp*bdf)*dt*temp
     endif

     if(numvar.ge.3) then
        call PVS3(lin,temp79c)
        temp = int3(vip79(:,OP_1),temp79b,temp79c)
        ssterm(chi_g) = ssterm(chi_g) +     thimp     *dt*temp
        ddterm(chi_g) = ddterm(chi_g) - (1.-thimp*bdf)*dt*temp
     endif
  endif

end subroutine axial_vel_lin

subroutine axial_vel_nolin(trial, r4term)

  use basic
  use m3dc1_nint
  use metricterms_new

  implicit none

  vectype, intent(in), dimension(MAX_PTS, OP_NUM)  :: trial
  vectype, intent(out) :: r4term
  
  r4term = 0.

  if(numvar.lt.2 .or. linear.eq.1) return

  ! density terms
  ! ~~~~~~~~~~~~~
  if(idens.eq.1 .and. eqsubtract.eq.1) then
     r4term = r4term + dt* &
          (v2vun(trial,vz079,ph079,n179) &
          +v2vs (trial,vz079,sig79))
  endif

end subroutine axial_vel_nolin


!======================================================================
! Compression Equation
!======================================================================
subroutine compression_lin(trial, lin, ssterm, ddterm, q_bf, advfield)
  
  use basic
  use arrays
  use m3dc1_nint
  use metricterms_new

  implicit none

  vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: trial, lin 
  vectype, dimension(num_fields), intent(out) :: ssterm, ddterm
  vectype, intent(out) :: q_bf
  integer, intent(in) :: advfield

  vectype :: temp
  real :: ththm

  select case(imp_mod)
  case(0)
     ththm = (1.-thimp*bdf)*thimp
  case(1)
     ththm = -bdf*thimp**2
  case(2)
     ththm = -bdf*thimp
  end select

  ssterm = 0.
  ddterm = 0.
  q_bf = 0. 
                     
  if(numvar.lt.3) return

  if(istatic.eq.1) then
     if(.not.surface_int) then
        temp = int2(trial,lin)
        ssterm(chi_g) = temp
        ddterm(chi_g) = temp
     endif
     return
  endif

  ! regularize the chi equation
  if(inoslip_pol.eq.0 .and. (.not.surface_int)) then
     temp = -regular*int2(trial(:,OP_1),lin(:,OP_1))
     ssterm(chi_g) = ssterm(chi_g) + temp
     ddterm(chi_g) = ddterm(chi_g) + temp*bdf
  end if
         
  temp = v3un(trial,lin,nt79)
  ssterm(u_g) = ssterm(u_g) + temp
  ddterm(u_g) = ddterm(u_g) + temp*bdf
  
  temp = v3chin(trial,lin,nt79)*chiiner
  ssterm(chi_g) = ssterm(chi_g) + temp
  ddterm(chi_g) = ddterm(chi_g) + temp*bdf

  if(linear.eq.0) then 
     temp = v3uun  (trial,lin,ph179,nt79) &
          + v3uun  (trial,ph179,lin,nt79) &
          + v3uvn  (trial,lin,vz179,nt79) &
          + v3uchin(trial,lin,ch179,nt79)
     ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
     ddterm(u_g) = ddterm(u_g) + (.5-thimp*bdf)*dt*temp

     temp = v3uvn  (trial,ph179,lin,nt79) &
          + v3vvn  (trial,lin,vz179,nt79) &
          + v3vvn  (trial,vz179,lin,nt79) &
          + v3vchin(trial,lin,ch179,nt79)
     ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
     ddterm(vz_g) = ddterm(vz_g) + (.5-thimp*bdf)*dt*temp

     temp = v3uchin  (trial,ph179,lin,nt79) &
          + v3vchin  (trial,vz179,lin,nt79) &
          + v3chichin(trial,lin,ch179,nt79) &
          + v3chichin(trial,ch179,lin,nt79)  
     ssterm(chi_g) = ssterm(chi_g) -     thimp     *dt*temp
     ddterm(chi_g) = ddterm(chi_g) + (.5-thimp*bdf)*dt*temp
  endif
           
  temp = v3umu(trial,lin,vis79,vic79) &
       + v3us (trial,lin,sig79)
  ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
  ddterm(u_g) = ddterm(u_g) + (1.-thimp*bdf)*dt*temp

  temp = v3vmu(trial,lin,vis79,vic79)
  ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
  ddterm(vz_g) = ddterm(vz_g) + (1.-thimp*bdf)*dt*temp
                     
  temp = v3chimu(trial,lin,vis79,vic79) &
       + v3chis (trial,lin,sig79)
  ssterm(chi_g) = ssterm(chi_g) -     thimp     *dt*temp
  ddterm(chi_g) = ddterm(chi_g) + (1.-thimp*bdf)*dt*temp
  
  if(gyro.eq.1) then    
     temp = g3u(trial,lin)*dbf
     ssterm(u_g) = ssterm(u_g) +     thimp     *dt*temp
     ddterm(u_g) = ddterm(u_g) - (1.-thimp*bdf)*dt*temp
     
     temp = g3v(trial,lin)*dbf
     ssterm(vz_g) = ssterm(vz_g) +     thimp     *dt*temp
     ddterm(vz_g) = ddterm(vz_g) - (1.-thimp*bdf)*dt*temp
     
     temp = g3chi(trial,lin)*dbf
     ssterm(chi_g) = ssterm(chi_g) +     thimp     *dt*temp
     ddterm(chi_g) = ddterm(chi_g) - (1.-thimp*bdf)*dt*temp
  endif
  
  if(advfield.eq.1) then
     temp = v3up     (trial,lin,pt79)        &
          + v3upsipsi(trial,lin,pst79,pst79) &
          + v3upsib  (trial,lin,pst79,bzt79) &
          + v3ubb    (trial,lin,bzt79,bzt79) & 
          + v3ungrav (trial,lin,nt79)
     ssterm(u_g) = ssterm(u_g) - thimp*thimp*dt*dt*temp
     ddterm(u_g) = ddterm(u_g) +       ththm*dt*dt*temp

     temp = v3vp     (trial,lin,pt79)        &
          + v3vpsipsi(trial,lin,pst79,pst79) &
          + v3vpsib  (trial,lin,pst79,bzt79) &
          + v3vbb    (trial,lin,bzt79,bzt79)
     ssterm(vz_g) = ssterm(vz_g) - thimp*thimp*dt*dt*temp
     ddterm(vz_g) = ddterm(vz_g) +       ththm*dt*dt*temp

     temp = v3chip     (trial,lin,pt79)        &
          + v3chipsipsi(trial,lin,pst79,pst79) &
          + v3chipsib  (trial,lin,pst79,bzt79) &
          + v3chibb    (trial,lin,bzt79,bzt79) &
          + v3chingrav (trial,lin,nt79)
     ssterm(chi_g) = ssterm(chi_g) - thimp*thimp*dt*dt*temp
     ddterm(chi_g) = ddterm(chi_g) +       ththm*dt*dt*temp

     ddterm(psi_g) = ddterm(psi_g) + dt*  &
          (v3psipsi(trial,lin,pss79)      & 
          +v3psipsi(trial,pss79,lin)      &
          +v3psib  (trial,lin,bzs79))
     
     ddterm(bz_g) = ddterm(bz_g) + dt*  &
          (v3psib(trial,pss79,lin)      &
          +v3bb  (trial,lin,bzs79)      &  
          +v3bb  (trial,bzs79,lin))
  
     ddterm(p_g) = ddterm(p_g) + dt*  &
          v3p(trial,lin)
  else
     if(linear.eq.0) then 
        temp = v3psipsi(trial,lin,ps179) &
             + v3psipsi(trial,ps179,lin) &
             + v3psib  (trial,lin,bz179)
        ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
        ddterm(psi_g) = ddterm(psi_g) + (.5-thimp*bdf)*dt*temp
        
        temp = v3psib(trial,ps179,lin) &
             + v3bb  (trial,lin,bz179) &
             + v3bb  (trial,bz179,lin)
        ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
        ddterm(bz_g) = ddterm(bz_g) + (.5-thimp*bdf)*dt*temp
     endif

     temp = v3p(trial,lin) 
     ssterm(p_g) = ssterm(p_g) -     thimp     *dt*temp
     ddterm(p_g) = ddterm(p_g) + (1.-thimp*bdf)*dt*temp
  endif
  
  if(eqsubtract.eq.1) then             
     temp = v3uun  (trial,lin,ph079,nt79) &
          + v3uun  (trial,ph079,lin,nt79) &
          + v3uvn  (trial,lin,vz079,nt79) &
          + v3uchin(trial,lin,ch079,nt79)
     ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
     ddterm(u_g) = ddterm(u_g) + (1.-thimp*bdf)*dt*temp
     
     temp = v3uvn  (trial,ph079,lin,nt79) &
          + v3vvn  (trial,lin,vz079,nt79) &
          + v3vvn  (trial,vz079,lin,nt79) &
          + v3vchin(trial,lin,ch079,nt79)
     ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
     ddterm(vz_g) = ddterm(vz_g) + (1.-thimp*bdf)*dt*temp
     
     temp = v3uchin  (trial,ph079,lin,nt79) &
          + v3vchin  (trial,vz079,lin,nt79) &
          + v3chichin(trial,lin,ch079,nt79) &
          + v3chichin(trial,ch079,lin,nt79)
     ssterm(chi_g) = ssterm(chi_g) -     thimp     *dt*temp
     ddterm(chi_g) = ddterm(chi_g) + (1.-thimp*bdf)*dt*temp

     if(advfield.eq.0) then
        temp = v3psipsi(trial,lin,ps079) &
             + v3psipsi(trial,ps079,lin) &
             + v3psib  (trial,lin,bz079)
        ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
        ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

        temp = v3psib(trial,ps079,lin) &
             + v3bb  (trial,lin,bz079) &
             + v3bb  (trial,bz079,lin)
        ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
        ddterm(bz_g) = ddterm(bz_g) + (1.-thimp*bdf)*dt*temp
     end if
  endif

  if(idens.eq.1) then
     if(advfield.eq.1) then 
        ddterm(den_g) = ddterm(den_g) + dt* &
             v3ngrav(trial,lin)
     else
        temp = v3ngrav(trial,lin)
        ssterm(den_g) = ssterm(den_g) -     thimp     *dt*temp
        ddterm(den_g) = ddterm(den_g) + (1.-thimp*bdf)*dt*temp
     endif

     if(eqsubtract.eq.1) then 
        ddterm(den_g) = ddterm(den_g) + dt* &
             (v3uun    (trial,ph079,ph079,lin) &
             +v3uvn    (trial,ph079,vz079,lin) &
             +v3uchin  (trial,ph079,ch079,lin) &
             +v3vvn    (trial,vz079,vz079,lin) &
             +v3vchin  (trial,vz079,ch079,lin) &
             +v3chichin(trial,ch079,ch079,lin))
     endif
  endif

  ! Parallel Viscosity
  if(amupar.ne.0.) then
     call PVV3(trial,temp79b)

     call PVS1(lin,temp79c)
     temp = int3(vip79(:,OP_1),temp79b,temp79c)
     ssterm(u_g) = ssterm(u_g) +     thimp     *dt*temp
     ddterm(u_g) = ddterm(u_g) - (1.-thimp*bdf)*dt*temp

     if(numvar.ge.2) then
        call PVS2(lin,temp79c)
        temp = int3(vip79(:,OP_1),temp79b,temp79c)
        ssterm(vz_g) = ssterm(vz_g) +     thimp     *dt*temp
        ddterm(vz_g) = ddterm(vz_g) - (1.-thimp*bdf)*dt*temp
     endif

     if(numvar.ge.3) then
        call PVS3(lin,temp79c)
        temp = int3(vip79(:,OP_1),temp79b,temp79c)
        ssterm(chi_g) = ssterm(chi_g) +     thimp     *dt*temp
        ddterm(chi_g) = ddterm(chi_g) - (1.-thimp*bdf)*dt*temp
     endif
  endif


  ! 3D terms
  if(i3d.eq.1) then
     temp = v3psif(trial,lin,bft79)
     ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

     temp = v3bf(trial,lin,bft79)
     ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
     ddterm(bz_g) = ddterm(bz_g) + (1.-thimp*bdf)*dt*temp

     if(eqsubtract.eq.1) then
        ! The v3bf term causes problems in the isplitstep=0 case
        ! for some unknown reason.
        q_bf = q_bf + dt* &
             (v3psif(trial,ps079,lin) &
             +v3bf  (trial,bz079,lin))
     endif
  endif

end subroutine compression_lin

subroutine compression_nolin(trial, r4term)

  use basic
  use m3dc1_nint
  use metricterms_new

  implicit none

  vectype, intent(in), dimension(MAX_PTS, OP_NUM)  :: trial
  vectype, intent(out) :: r4term
  
  r4term = 0.

  if(numvar.lt.3 .or. linear.eq.1) return

  ! density terms
  ! ~~~~~~~~~~~~~
  if(idens.eq.1 .and. eqsubtract.eq.1) then
     r4term = r4term + dt* &
          (v3us     (trial,ph079,sig79) &
          +v3chis   (trial,ch079,sig79))

  endif

end subroutine compression_nolin


!======================================================================
! Flux Equation
!======================================================================
subroutine flux_lin(trial, lin, ssterm, ddterm, q_ni, q_bf, r_e, q_e)
  
  use basic
  use arrays
  use m3dc1_nint
  use metricterms_new
  use harned_mikic_mod

  implicit none

  vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: trial, lin 
  vectype, dimension(num_fields), intent(out) :: ssterm, ddterm
  vectype, intent(out) :: q_ni(2), q_bf, r_e, q_e
  vectype :: temp, temp2
  real :: thimpb, thimpf, thimpe

  vectype, dimension(MAX_PTS, OP_NUM) :: hf
  hf = hypf*sz79


  if(imp_mod.eq.0) then
     thimpb = thimp
  else
     thimpb = 1.
  endif

  thimpf = thimp
  thimpe = 1.

  ssterm = 0.
  ddterm = 0.
  q_ni = 0.
  q_bf = 0.
  r_e = 0.
  q_e = 0.

  if(iestatic.eq.1) then
     if(.not.surface_int) then
        temp = int2(trial,lin)
        ssterm(psi_g) = temp
        ddterm(psi_g) = temp
     endif
     return
  end if

  temp = b1psi (trial,lin) &
       - b1psid(trial,lin,ni79)   ! electron mass term
  ssterm(psi_g) = ssterm(psi_g) + temp
  ddterm(psi_g) = ddterm(psi_g) + temp*bdf

  temp = b1psieta(trial,lin,eta79,hf)
  ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
  ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

  temp = b1psiu  (trial,lin,pht79) &
       + b1psiv  (trial,lin,vzt79) &
       + b1psichi(trial,lin,cht79)
  ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
  ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

  if(linear.eq.0) then 
     temp = (b1psipsid(trial,lin,ps179,ni79) &
            +b1psipsid(trial,ps179,lin,ni79))*dbf
     ssterm(psi_g) = ssterm(psi_g) -     thimpf     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (.5-thimpf*bdf)*dt*temp

     temp = b1psibd(trial,lin,bz179,ni79)*dbf
     ssterm(psi_g) = ssterm(psi_g) -     thimpf     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (.5-thimpf*bdf)*dt*temp

     temp = b1psiu(trial,ps179,lin)
     ssterm(u_g) = ssterm(u_g) - thimpb*dt*temp
     ddterm(u_g) = ddterm(u_g) - thimpb*dt*temp*bdf

     temp = b1bu(trial,bz179,lin)
     ssterm(u_g) = ssterm(u_g) - thimpb*dt*temp
     if(numvar.eq.1) then
        ddterm(u_g) = ddterm(u_g) - (1.-2.*thimpb)*dt*temp*bdf
     else
        ddterm(u_g) = ddterm(u_g) - thimpb*dt*temp*bdf
     endif
  endif

  if(eqsubtract.eq.1) then
     temp = (b1psipsid(trial,lin,ps079,ni79) &
            +b1psipsid(trial,ps079,lin,ni79) &
            +b1psibd  (trial,lin,bz079,ni79))*dbf 
     ssterm(psi_g) = ssterm(psi_g) -     thimpf     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (1.-thimpf*bdf)*dt*temp

     temp = b1psiu(trial,ps079,lin) &
          + b1bu  (trial,bz079,lin)
     ssterm(u_g) = ssterm(u_g) -     thimpb     *dt*temp
     ddterm(u_g) = ddterm(u_g) + (1.-thimpb*bdf)*dt*temp
  endif

  call b1harnedmikic(trial,lin,temp,temp2)
  ssterm(psi_g) = ssterm(psi_g) - thimpf**2*dt*dt*temp
  ddterm(psi_g) = ddterm(psi_g) - thimpf**2*dt*dt*temp
  if(numvar.ge.2) then 
     ssterm(bz_g) = ssterm(bz_g) - thimpf**2*dt*dt*temp2
     ddterm(bz_g) = ddterm(bz_g) - thimpf**2*dt*dt*temp2
  endif


  ! NUMVAR = 2
  ! ~~~~~~~~~~
  if(numvar.ge.2) then
     if(linear.eq.0) then 
        temp = b1psibd(trial,ps179,lin,ni79)*dbf &
             + b1bbd  (trial,bz179,lin,ni79)*dbf &
             + b1bbd  (trial,lin,bz179,ni79)*dbf
        ssterm(bz_g) = ssterm(bz_g) -     thimpf     *dt*temp
        ddterm(bz_g) = ddterm(bz_g) + (.5-thimpf*bdf)*dt*temp
        
        temp = b1psiv(trial,ps179,lin) &
             + b1bv  (trial,bz179,lin)  
        ssterm(vz_g) = ssterm(vz_g) - thimpb*dt*temp
        ddterm(vz_g) = ddterm(vz_g) - thimpb*dt*temp*bdf
     endif

     temp = b1beta(trial,lin,eta79,hf)
     ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
     ddterm(bz_g) = ddterm(bz_g) + (1.-thimp*bdf)*dt*temp

     temp = b1bu(trial,lin,pht79) &
          + b1bv(trial,lin,vzt79)
     ssterm(bz_g) = ssterm(bz_g) -     thimpb     *dt*temp
     ddterm(bz_g) = ddterm(bz_g) + (1.-thimpb*bdf)*dt*temp
                      
     if(eqsubtract.eq.1) then             
        temp = b1psibd(trial,ps079,lin,ni79)*dbf &
             + b1bbd  (trial,bz079,lin,ni79)*dbf & 
             + b1bbd  (trial,lin,bz079,ni79)*dbf   
        ssterm(bz_g) = ssterm(bz_g) -     thimpf     *dt*temp
        ddterm(bz_g) = ddterm(bz_g) + (1.-thimpf*bdf)*dt*temp

        temp = b1psiv(trial,ps079,lin) &
             + b1bv  (trial,bz079,lin)
        ssterm(vz_g) = ssterm(vz_g) -     thimpb     *dt*temp
        ddterm(vz_g) = ddterm(vz_g) + (1.-thimpb*bdf)*dt*temp
     end if
  end if

  if(i3d.eq.1 .and. numvar.ge.2) then
     temp = b1fu(trial,bft79,lin)
     ssterm(u_g) = ssterm(u_g) -     thimpb     *dt*temp
     ddterm(u_g) = ddterm(u_g) + (1.-thimpb*bdf)*dt*temp

     temp = b1fv(trial,bft79,lin)
     ssterm(vz_g) = ssterm(vz_g) -     thimpb     *dt*temp
     ddterm(vz_g) = ddterm(vz_g) + (1.-thimpb*bdf)*dt*temp

     temp = b1psifd(trial,lin,bft79,ni79)*dbf
     ssterm(psi_g) = ssterm(psi_g) -     thimpf     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (1.-thimpf*bdf)*dt*temp

     temp = b1bfd(trial,lin,bft79,ni79)*dbf
     ssterm(bz_g) = ssterm(bz_g) -     thimpf     *dt*temp
     ddterm(bz_g) = ddterm(bz_g) + (1.-thimpf*bdf)*dt*temp

     if(eqsubtract.eq.1) then
        q_bf = q_bf + dt* &
             (b1fu   (trial,lin,ch079) &
             +b1fv   (trial,lin,vz079) &
             +b1psifd(trial,ps079,lin,ni79)*dbf &
             +b1bfd  (trial,bz079,lin,ni79)*dbf)
     endif
  endif


  ! NUMVAR = 3
  ! ~~~~~~~~~~
  if(numvar.ge.3) then
     temp = b1ped(trial,lin,ni79)*dbf*pefac
     ssterm(pe_g) = ssterm(pe_g) -     thimpf     *dt*temp
     ddterm(pe_g) = ddterm(pe_g) + (1.-thimpf*bdf)*dt*temp

     if(linear.eq.0) then 
        temp = b1psichi(trial,ps179,lin) &
             + b1bchi  (trial,bz179,lin)
        ssterm(chi_g) = ssterm(chi_g) - thimpb*dt*temp
        ddterm(chi_g) = ddterm(chi_g) - thimpb*dt*temp*bdf
     endif
     
     if(eqsubtract.eq.1) then
        temp = b1psichi(trial,ps079,lin) &
             + b1bchi  (trial,bz079,lin)
        ssterm(chi_g) = ssterm(chi_g) -     thimpb     *dt*temp
        ddterm(chi_g) = ddterm(chi_g) + (1.-thimpb*bdf)*dt*temp
     endif

     if(i3d.eq.1) then 
        temp = b1fchi(trial,bft79,lin)
        ssterm(chi_g) = ssterm(chi_g) -     thimpb     *dt*temp
        ddterm(chi_g) = ddterm(chi_g) + (1.-thimpb*bdf)*dt*temp

        if(eqsubtract.eq.1) then
           q_bf = q_bf + dt* &
                (b1fchi(trial,lin,ch079))
        endif
     endif
  end if

  ! density terms
  ! ~~~~~~~~~~~~~
  if(idens.eq.1 .and. eqsubtract.eq.1) then
     q_ni(1) = q_ni(1) + dt* &
          (b1psipsid(trial,ps079,ps079,lin)*dbf &
          +b1psibd  (trial,ps079,bz079,lin)*dbf &
          +b1bbd    (trial,bz079,bz079,lin)*dbf &
          +b1ped    (trial,pe079,lin)*dbf)
  endif

  ! electrostatic potential
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  temp = b1e(trial,lin)
  r_e = r_e       - thimpe     *dt*temp
  q_e = q_e + (1. - thimpe*bdf)*dt*temp
        
end subroutine flux_lin

subroutine flux_nolin(trial, r4term)

  use math
  use basic
  use m3dc1_nint
  use metricterms_new

  implicit none

  vectype, intent(in), dimension(MAX_PTS, OP_NUM)  :: trial
  vectype, intent(out) :: r4term
  
  r4term = 0.

  if(igauge.eq.1 .or. linear.eq.1) then
     r4term = r4term - dt* &
          vloop*int1(trial)/twopi
  endif
 
end subroutine flux_nolin


!======================================================================
! Axial Field Equation
!======================================================================
subroutine axial_field_lin(trial, lin, ssterm, ddterm, q_ni, q_bf)
  
  use basic
  use arrays
  use m3dc1_nint
  use metricterms_new
  use harned_mikic_mod

  implicit none

  vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: trial, lin 
  vectype, dimension(num_fields), intent(out) :: ssterm, ddterm
  vectype, intent(out) :: q_ni(2), q_bf
  vectype :: temp, temp2
  real :: thimpb, thimpf

  vectype, dimension(MAX_PTS, OP_NUM) :: hi
  hi = hypi*sz79

  if(imp_mod.eq.0) then
     thimpb = thimp
  else
     thimpb = 1.
  endif

  thimpf = thimp

  ssterm = 0.
  ddterm = 0.
  q_ni = 0.
  q_bf = 0.

  if(iestatic.eq.1) then
     if(.not.surface_int) then
        temp = int2(trial,lin)
        ssterm(bz_g) = temp
        ddterm(bz_g) = temp
     endif
     return
  end if


  if(numvar.lt.2) return          

  temp = b2psieta(trial,lin,eta79,hi)
  ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
  ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

!  if(ibootstrap.gt.0) then
!     temp = b2psimue(trial,lin,)*dbf
!     ssterm(psi_g) = ssterm(psi_g) -     thimpf     *dt*temp
!     ddterm(psi_g) = ddterm(psi_g) + (1.-thimpf*bdf)*dt*temp
!  endif
         
  temp = b2psiv(trial,lin,vzt79)
  ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
  ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

  if(linear.eq.0) then 
     temp = b2psipsid(trial,lin,ps179,ni79)*dbf &
          + b2psipsid(trial,ps179,lin,ni79)*dbf &
          + b2psibd  (trial,lin,bz179,ni79)*dbf
     ssterm(psi_g) = ssterm(psi_g) -     thimpf     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (.5-thimpf*bdf)*dt*temp

     temp = b2psibd(trial,ps179,lin,ni79)*dbf &
          + b2bbd  (trial,lin,bz179,ni79)*dbf &
          + b2bbd  (trial,bz179,lin,ni79)*dbf
     ssterm(bz_g) = ssterm(bz_g) -     thimpf     *dt*temp
     ddterm(bz_g) = ddterm(bz_g) + (.5-thimpf*bdf)*dt*temp

     temp = b2bu  (trial,bz179,lin)
     ssterm(u_g) = ssterm(u_g) - thimpb*dt*temp
     ddterm(u_g) = ddterm(u_g) - thimpb*dt*temp*bdf

     temp = b2psiv(trial,ps179,lin)
     ssterm(vz_g) = ssterm(vz_g) - thimpb*dt*temp
     ddterm(vz_g) = ddterm(vz_g) - thimpb*dt*temp*bdf
  endif

  temp = b2b (trial,lin) &
       - b2bd(trial,lin,ni79)   ! electron mass term
  ssterm(bz_g) = ssterm(bz_g) + temp
  ddterm(bz_g) = ddterm(bz_g) + temp*bdf

  temp = b2beta(trial,lin,eta79,hi)
  ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
  ddterm(bz_g) = ddterm(bz_g) + (1.-thimp*bdf)*dt*temp

  temp = b2bu(trial,lin,pht79)
  ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
  ddterm(bz_g) = ddterm(bz_g) + (1.-thimp*bdf)*dt*temp

  call b2harnedmikic(trial,lin,temp,temp2)
  ssterm(psi_g) = ssterm(psi_g) - thimpf**2*dt*dt*temp
  ddterm(psi_g) = ddterm(psi_g) - thimpf**2*dt*dt*temp
  ssterm(bz_g) = ssterm(bz_g) - thimpf**2*dt*dt*temp2
  ddterm(bz_g) = ddterm(bz_g) - thimpf**2*dt*dt*temp2
 

  if(eqsubtract.eq.1) then
     temp = b2psipsid(trial,lin,ps079,ni79)*dbf &
          + b2psipsid(trial,ps079,lin,ni79)*dbf &
          + b2psibd  (trial,lin,bz079,ni79)*dbf
     ssterm(psi_g) = ssterm(psi_g) -     thimpf     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (1.-thimpf*bdf)*dt*temp

     temp = b2psibd(trial,ps079,lin,ni79)*dbf &
          + b2bbd  (trial,lin,bz079,ni79)*dbf &
          + b2bbd  (trial,bz079,lin,ni79)*dbf
     ssterm(bz_g) = ssterm(bz_g) -     thimpf     *dt*temp
     ddterm(bz_g) = ddterm(bz_g) + (1.-thimpf*bdf)*dt*temp

     temp = b2bu  (trial,bz079,lin)
     ssterm(u_g) = ssterm(u_g) -     thimpb     *dt*temp
     ddterm(u_g) = ddterm(u_g) + (1.-thimpb*bdf)*dt*temp
     
     temp = b2psiv(trial,ps079,lin)
     ssterm(vz_g) = ssterm(vz_g) -     thimpb     *dt*temp
     ddterm(vz_g) = ddterm(vz_g) + (1.-thimpb*bdf)*dt*temp
  endif

  if(i3d.eq.1) then
     temp = b2fv(trial,bft79,lin)
     ssterm(vz_g) = ssterm(vz_g) -     thimpb *dt*temp
     ddterm(vz_g) = ddterm(vz_g) + (1.-thimpb)*dt*temp*bdf

     temp = b2psifd(trial,lin,bft79,ni79)*dbf
     ssterm(psi_g) = ssterm(psi_g) -     thimpf     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (1.-thimpf*bdf)*dt*temp

     temp = b2bfd(trial,lin,bft79,ni79)*dbf
     ssterm(bz_g) = ssterm(bz_g) -     thimpf     *dt*temp
     ddterm(bz_g) = ddterm(bz_g) + (1.-thimpf*bdf)*dt*temp

     q_bf = q_bf + dt* &
          b2feta(trial,lin,eta79,hi)

     if(eqsubtract.eq.1) then
        q_bf = q_bf + dt* &
             (b2psifd(trial,ps079,lin,ni79)*dbf & 
             +b2bfd  (trial,bz079,lin,ni79)*dbf &
             +b2fv   (trial,lin,vz079))
     endif
  endif


  ! NUMVAR = 3
  ! ~~~~~~~~~~
  if(numvar.ge.3) then
     temp = b2bchi(trial,lin,cht79)                  
     ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
     ddterm(bz_g) = ddterm(bz_g) + (1.-thimp*bdf)*dt*temp

     temp = b2ped(trial,lin,ni79)*dbf*pefac
     ssterm(pe_g) = ssterm(pe_g) -     thimpf     *dt*temp
     ddterm(pe_g) = ddterm(pe_g) + (1.-thimpf*bdf)*dt*temp

     if(linear.eq.0) then
        temp = b2bchi(trial,bz179,lin)                  
        ssterm(chi_g) = ssterm(chi_g) - thimpb*dt*temp
        ddterm(chi_g) = ddterm(chi_g) - thimpb*dt*temp*bdf
     endif
             
     if(eqsubtract.eq.1) then
        temp = b2bchi(trial,bz079,lin)
        ssterm(chi_g) = ssterm(chi_g) -     thimpb     *dt*temp
        ddterm(chi_g) = ddterm(chi_g) + (1.-thimpb*bdf)*dt*temp
     endif
  end if


  ! density terms
  ! ~~~~~~~~~~~~~
  if(idens.eq.1 .and. eqsubtract.eq.1) then
     q_ni(1) = q_ni(1) + dt* &
          (b2psipsid(trial,ps079,ps079,lin)*dbf &
          +b2psibd  (trial,ps079,bz079,lin)*dbf &
          +b2bbd    (trial,bz079,bz079,lin)*dbf &
          +b2ped    (trial,pe079,lin)*dbf*pefac)
  endif

end subroutine axial_field_lin


subroutine axial_field_nolin(trial, r4term)

  use basic
  use m3dc1_nint
  use metricterms_new

  implicit none

  vectype, intent(in), dimension(MAX_PTS, OP_NUM)  :: trial
  vectype, intent(out) :: r4term
  
  r4term = 0.

end subroutine axial_field_nolin


!======================================================================
! Electron Pressure Equation
!======================================================================
subroutine electron_pressure_lin(trial, lin, ssterm, ddterm, q_ni, q_bf)
  
  use basic
  use arrays
  use m3dc1_nint
  use metricterms_new

  implicit none

  vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: trial, lin 
  vectype, dimension(num_fields), intent(out) :: ssterm, ddterm
  vectype, intent(out) :: q_ni(2), q_bf
  vectype :: temp
  real :: thimpb, thimpf

  vectype, dimension(MAX_PTS, OP_NUM) :: hp
  hp = hypp*sz79


  if(imp_mod.eq.0) then
     thimpb = thimp
  else
     thimpb = 1.
  endif

  thimpf = thimp

  ssterm = 0.
  ddterm = 0.
  q_ni = 0.
  q_bf = 0.

  if(numvar.lt.3) return   

  temp = b3pe(trial,lin)
  ssterm(pe_g) = ssterm(pe_g) + temp
  ddterm(pe_g) = ddterm(pe_g) + temp*bdf

  if(linear.eq.0) then
    temp = b3psipsieta(trial,lin,ps179,eta79) &
         + b3psipsieta(trial,ps179,lin,eta79)
    ssterm(psi_g) = ssterm(psi_g) -     thimp_ohm     *dt*temp
    ddterm(psi_g) = ddterm(psi_g) + (.5-thimp_ohm*bdf)*dt*temp

    temp = b3bbeta(trial,lin,bz179,eta79) &
         + b3bbeta(trial,bz179,lin,eta79)
    ssterm(bz_g) = ssterm(bz_g) -     thimp_ohm     *dt*temp
    ddterm(bz_g) = ddterm(bz_g) + (.5-thimp_ohm*bdf)*dt*temp

    temp = b3pepsid(trial,pe179,lin,ni79)*dbf*pefac
    ssterm(psi_g) = ssterm(psi_g) -     thimpf     *dt*temp
    ddterm(psi_g) = ddterm(psi_g) + (.5-thimpf*bdf)*dt*temp
  
    temp = b3pebd(trial,pe179,lin,ni79)*dbf*pefac
    ssterm(bz_g) = ssterm(bz_g) -     thimpf     *dt*temp
    ddterm(bz_g) = ddterm(bz_g) + (.5-thimpf*bdf)*dt*temp

    temp = b3pepsid(trial,lin,ps179,ni79)*dbf*pefac & 
         + b3pebd  (trial,lin,bz179,ni79)*dbf*pefac
    ssterm(pe_g) = ssterm(pe_g) -     thimpf     *dt*temp
    ddterm(pe_g) = ddterm(pe_g) + (.5-thimpf*bdf)*dt*temp

    temp = p1pu(trial,pe179,lin)
    ssterm(u_g) = ssterm(u_g) - thimpb*dt*temp
    ddterm(u_g) = ddterm(u_g) - thimpb*dt*temp*bdf

    temp = p1pv(trial,pe179,lin)
    ssterm(vz_g) = ssterm(vz_g) - thimpb*dt*temp
    ddterm(vz_g) = ddterm(vz_g) - thimpb*dt*temp*bdf

    temp = p1pchi(trial,pe179,lin)
    ssterm(chi_g) = ssterm(chi_g) - thimpb*dt*temp
    ddterm(chi_g) = ddterm(chi_g) - thimpb*dt*temp*bdf

    if(ipres.eq.0) then
       temp = p1uus  (trial,lin,ph179,sig79) &
            + p1uus  (trial,ph179,lin,sig79) &
            + p1uchis(trial,lin,ch179,sig79) 
       ssterm(u_g) = ssterm(u_g) -     thimpb     *dt*temp
       ddterm(u_g) = ddterm(u_g) + (.5-thimpb*bdf)*dt*temp
       
       temp = p1vvs  (trial,lin,vz179,sig79) &
            + p1vvs  (trial,vz179,lin,sig79)
       ssterm(vz_g) = ssterm(vz_g) -     thimpb     *dt*temp
       ddterm(vz_g) = ddterm(vz_g) + (.5-thimpb*bdf)*dt*temp
       
       temp = p1chichis(trial,lin,ch179,sig79) &
            + p1chichis(trial,ch179,lin,sig79) &
            + p1uchis  (trial,ph179,lin,sig79) 
       ssterm(chi_g) = ssterm(chi_g) -     thimpb     *dt*temp
       ddterm(chi_g) = ddterm(chi_g) + (.5-thimpb*bdf)*dt*temp
    endif
  endif

  temp = p1pu  (trial,lin,pht79) & 
       + p1pv  (trial,lin,vzt79) &
       + p1pchi(trial,lin,cht79) &
       + b3pedkappa(trial,lin,ni79,kap79,hp)
  ssterm(pe_g) = ssterm(pe_g) -     thimp     *dt*temp
  ddterm(pe_g) = ddterm(pe_g) + (1.-thimp*bdf)*dt*temp

  ! Anisotropic Heat Flux
  if(kappar.ne.0.) then
     ! Assumes no contribution from equilibrium f

     if(linear.eq.0) then
        temp = p1psipsikappar(trial,lin,ps179,pe179,ni79,b2i79,kar79) &
             + p1psipsikappar(trial,ps179,lin,pe179,ni79,b2i79,kar79) &
             + p1psibkappar  (trial,lin,bz179,pe179,ni79,b2i79,kar79)
        ssterm(psi_g) = ssterm(psi_g) -          thimp     *dt*temp
        ddterm(psi_g) = ddterm(psi_g) + (1./3. - thimp*bdf)*dt*temp

        temp = p1psibkappar  (trial,ps179,lin,pe179,ni79,b2i79,kar79) &
             + p1bbkappar    (trial,lin,bz179,pe179,ni79,b2i79,kar79) &
             + p1bbkappar    (trial,bz179,lin,pe179,ni79,b2i79,kar79)
        ssterm(bz_g) = ssterm(bz_g) -          thimp     *dt*temp
        ddterm(bz_g) = ddterm(bz_g) + (1./3. - thimp*bdf)*dt*temp
        
        temp = p1psipsikappar(trial,ps179,ps179,lin,ni79,b2i79,kar79) &
             + p1psibkappar  (trial,ps179,bz179,lin,ni79,b2i79,kar79) &
             + p1bbkappar    (trial,bz179,bz179,lin,ni79,b2i79,kar79)
        ssterm(pe_g) = ssterm(pe_g) -          thimp     *dt*temp
        ddterm(pe_g) = ddterm(pe_g) + (1./3. - thimp*bdf)*dt*temp
     endif

     if(eqsubtract.eq.1) then 
        temp = p1psipsikappar(trial,lin,ps079,pe079,ni79,b2i79,kar79) &
             + p1psipsikappar(trial,ps079,lin,pe079,ni79,b2i79,kar79) &
             + p1psibkappar  (trial,lin,bz079,pe079,ni79,b2i79,kar79)
        ssterm(psi_g) = ssterm(psi_g) -       thimp     *dt*temp
        ddterm(psi_g) = ddterm(psi_g) + (1. - thimp*bdf)*dt*temp

        temp = p1psibkappar  (trial,ps079,lin,pe079,ni79,b2i79,kar79) &
             + p1bbkappar    (trial,lin,bz079,pe079,ni79,b2i79,kar79) &
             + p1bbkappar    (trial,bz079,lin,pe079,ni79,b2i79,kar79)
        ssterm(bz_g) = ssterm(bz_g) -       thimp     *dt*temp
        ddterm(bz_g) = ddterm(bz_g) + (1. - thimp*bdf)*dt*temp

        temp = p1psipsikappar(trial,ps079,ps079,lin,ni79,b2i79,kar79) &
             + p1psibkappar  (trial,ps079,bz079,lin,ni79,b2i79,kar79) &
             + p1bbkappar    (trial,bz079,bz079,lin,ni79,b2i79,kar79)
        ssterm(pe_g) = ssterm(pe_g) -       thimp     *dt*temp
        ddterm(pe_g) = ddterm(pe_g) + (1. - thimp*bdf)*dt*temp      

        if(linear.eq.0) then 
           temp = p1psipsikappar(trial,lin,ps179,pe079,ni79,b2i79,kar79) &
                + p1psipsikappar(trial,lin,ps079,pe179,ni79,b2i79,kar79) &
                + p1psipsikappar(trial,ps179,lin,pe079,ni79,b2i79,kar79) &
                + p1psipsikappar(trial,ps079,lin,pe179,ni79,b2i79,kar79) &
                + p1psibkappar  (trial,lin,bz179,pe079,ni79,b2i79,kar79) &
                + p1psibkappar  (trial,lin,bz179,pe179,ni79,b2i79,kar79)
           ssterm(psi_g) = ssterm(psi_g) -          thimp     *dt*temp
           ddterm(psi_g) = ddterm(psi_g) + (1./2. - thimp*bdf)*dt*temp

           temp = p1psibkappar  (trial,ps179,lin,pe079,ni79,b2i79,kar79) &
                + p1psibkappar  (trial,ps079,lin,pe179,ni79,b2i79,kar79) &
                + p1bbkappar    (trial,lin,bz179,pe079,ni79,b2i79,kar79) &
                + p1bbkappar    (trial,lin,bz079,pe179,ni79,b2i79,kar79) &
                + p1bbkappar    (trial,bz179,lin,pe079,ni79,b2i79,kar79) &
                + p1bbkappar    (trial,bz079,lin,pe179,ni79,b2i79,kar79)
           ssterm(bz_g) = ssterm(bz_g) -          thimp     *dt*temp
           ddterm(bz_g) = ddterm(bz_g) + (1./2. - thimp*bdf)*dt*temp
           
           temp = p1psipsikappar(trial,ps179,ps079,lin,ni79,b2i79,kar79) &
                + p1psipsikappar(trial,ps079,ps179,lin,ni79,b2i79,kar79) &
                + p1psibkappar  (trial,ps179,bz079,lin,ni79,b2i79,kar79) &
                + p1psibkappar  (trial,ps079,bz179,lin,ni79,b2i79,kar79) &
                + p1bbkappar    (trial,bz179,bz079,lin,ni79,b2i79,kar79) &
                + p1bbkappar    (trial,bz079,bz179,lin,ni79,b2i79,kar79)
           ssterm(pe_g) = ssterm(pe_g) -          thimp     *dt*temp
           ddterm(pe_g) = ddterm(pe_g) + (1./2. - thimp*bdf)*dt*temp      
        endif
     endif   

     if(i3d.eq.1) then
        if(eqsubtract.eq.1) then
           q_bf = q_bf + dt* &
                (p1psifkappar(trial,ps079,lin,pe079,ni79,b2i79,kar79) &
                +p1bfkappar  (trial,bz079,lin,pe079,ni79,b2i79,kar79))
        endif
     endif
  endif

  ! Cross-field Heat Flux
  if(kappax.ne.0.) then
     temp = p1kappax(trial,pe179,lin,ni79,kax79)
     ssterm(bz_g) = ssterm(bz_g) - thimp*dt*temp
     ddterm(bz_g) = ddterm(bz_g) - thimp*dt*temp*bdf

     temp = p1kappax(trial,lin,bzt79,ni79,kax79)
     ssterm(pe_g) = ssterm(pe_g) -     thimp     *dt*temp
     ddterm(pe_g) = ddterm(pe_g) + (1.-thimp*bdf)*dt*temp
  endif
              
  if(eqsubtract.eq.1) then
     temp = b3psipsieta(trial,lin,ps079,eta79) &
          + b3psipsieta(trial,ps079,lin,eta79)
     ssterm(psi_g) = ssterm(psi_g) -     thimp_ohm     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (1.-thimp_ohm*bdf)*dt*temp
     
     temp = b3bbeta(trial,lin,bz079,eta79) &
          + b3bbeta(trial,bz079,lin,eta79)
     ssterm(bz_g) = ssterm(bz_g) -     thimp_ohm     *dt*temp
     ddterm(bz_g) = ddterm(bz_g) + (1.-thimp_ohm*bdf)*dt*temp

     temp = b3pepsid(trial,pe079,lin,ni79)*dbf*pefac
     ssterm(psi_g) = ssterm(psi_g) -     thimpf     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (1.-thimpf*bdf)*dt*temp

     temp = b3pebd(trial,pe079,lin,ni79)*dbf*pefac
     ssterm(bz_g) = ssterm(bz_g) -     thimpf     *dt*temp
     ddterm(bz_g) = ddterm(bz_g) + (1.-thimpf*bdf)*dt*temp

     temp = b3pepsid(trial,lin,ps079,ni79)*dbf*pefac &
          + b3pebd  (trial,lin,bz079,ni79)*dbf*pefac
     ssterm(pe_g) = ssterm(pe_g) -     thimpf     *dt*temp
     ddterm(pe_g) = ddterm(pe_g) + (1.-thimpf*bdf)*dt*temp
     
     temp = p1pu(trial,pe079,lin)
     ssterm(u_g) = ssterm(u_g) -     thimpb     *dt*temp
     ddterm(u_g) = ddterm(u_g) + (1.-thimpb*bdf)*dt*temp
     
     temp = p1pv(trial,pe079,lin)
     ssterm(vz_g) = ssterm(vz_g) -     thimpb     *dt*temp
     ddterm(vz_g) = ddterm(vz_g) + (1.-thimpb*bdf)*dt*temp
     
     temp = p1pchi(trial,pe079,lin)                
     ssterm(chi_g) = ssterm(chi_g) -     thimpb     *dt*temp
     ddterm(chi_g) = ddterm(chi_g) + (1.-thimpb*bdf)*dt*temp
     
     if(ipres.eq.0) then
        temp = p1uus  (trial,lin,ph079,sig79) &
             + p1uus  (trial,ph079,lin,sig79) &
             + p1uchis(trial,lin,ch079,sig79) 
        ssterm(u_g) = ssterm(u_g) -     thimpb     *dt*temp
        ddterm(u_g) = ddterm(u_g) + (1.-thimpb*bdf)*dt*temp
        
        temp = p1vvs  (trial,lin,vz079,sig79) &
             + p1vvs  (trial,vz079,lin,sig79)
        ssterm(vz_g) = ssterm(vz_g) -     thimpb     *dt*temp
        ddterm(vz_g) = ddterm(vz_g) + (1.-thimpb*bdf)*dt*temp
        
        temp = p1chichis(trial,lin,ch079,sig79) &
             + p1chichis(trial,ch079,lin,sig79) &
             + p1uchis  (trial,ph079,lin,sig79) 
        ssterm(chi_g) = ssterm(chi_g) -     thimpb     *dt*temp
        ddterm(chi_g) = ddterm(chi_g) + (1.-thimpb*bdf)*dt*temp
     endif

  endif


  if(i3d.eq.1) then
     temp = b3pefd(trial,lin,bft79,ni79)*dbf*pefac
     ssterm(pe_g) = ssterm(pe_g) -     thimpf     *dt*temp
     ddterm(pe_g) = ddterm(pe_g) + (1.-thimpf*bdf)*dt*temp

     if(eqsubtract.eq.1) then
        q_bf = q_bf + dt* &
             (b3pefd(trial,pe079,lin,ni79)*dbf*pefac)
     endif
  endif


  ! density terms
  ! ~~~~~~~~~~~~~
  if(idens.eq.1 .and. eqsubtract.eq.1) then
     q_ni(1) = q_ni(1) + dt* &
          (b3pepsid(trial,pe079,ps079,lin)*dbf*pefac &
          +b3pebd  (trial,pe079,bz079,lin)*dbf*pefac &
          +b3pedkappa(trial,pe079,lin,kap79,hp) &
          +p1psipsikappar(trial,ps079,ps079,pe079,lin,b2i79,kar79) &
          +p1psibkappar  (trial,ps079,bz079,pe079,lin,b2i79,kar79) &
          +p1bbkappar    (trial,bz079,bz079,pe079,lin,b2i79,kar79))          
  endif

  if(eqsubtract.eq.1) then
     q_ni(2) = q_ni(2) + dt* &
          (b3pedkappa(trial,pe079,lin,kap79,hp) &
          +p1psipsikappar(trial,ps079,ps079,pe079,ni79,lin,kar79) &
          +p1psibkappar  (trial,ps079,bz079,pe079,ni79,lin,kar79) &
          +p1bbkappar    (trial,bz079,bz079,pe079,ni79,lin,kar79))
  endif

end subroutine electron_pressure_lin

subroutine electron_pressure_nolin(trial, r4term)

  use basic
  use m3dc1_nint
  use metricterms_new

  implicit none

  vectype, intent(in), dimension(MAX_PTS, OP_NUM)  :: trial
  vectype, intent(out) :: r4term

  vectype, dimension(MAX_PTS, OP_NUM) :: hv, hc
  hv = hypv*sz79
  hc = hypc*sz79
  
  r4term = 0.

  if(numvar.lt.3 .or. linear.eq.1) return

  ! source terms
  ! ~~~~~~~~~~~~
  if(gam.ne.1.) then
     ! hyper-ohmic heating
     r4term = r4term + dbf*dt*(gam-1.)* &
          (qpsipsieta(trial) &
          +qbbeta    (trial))

     ! viscous heating
     if(ipres.eq.0) then
        r4term = r4term - dt*(gam-1.)* &
             (quumu    (trial,pht79,pht79,vis79,      hc) &
             +qvvmu    (trial,vzt79,vzt79,vis79,      hv) &
             +quchimu  (trial,pht79,cht79,vis79,vic79,hc) &
             +qchichimu(trial,cht79,cht79,      vic79,hc) &
             +p1vip    (trial))
     endif
  end if

  ! density source terms
  ! ~~~~~~~~~~~~~~~~~~~~
  if(idens.eq.1 .and. eqsubtract.eq.1) then
     if(ipres.eq.0) then
        r4term = r4term + dt* &
             (p1uus    (trial,ph079,ph079,sig79) &
             +p1vvs    (trial,vz079,vz079,sig79) &
             +p1chichis(trial,ch079,ch079,sig79) &
             +p1uchis  (trial,ph079,ch079,sig79))
     endif
  endif

end subroutine electron_pressure_nolin


!======================================================================
! ludefall
! --------
!
! Clears, populates, and finalizes all matrices for the implicit 
! time advance.  Does not insert boundary conditions, or finalize the
! s* matrices.
!
!======================================================================
subroutine ludefall(ivel_def, idens_def, ipres_def, ifield_def)

  use p_data
  use mesh_mod
  use basic
  use arrays
  use sparse
  use m3dc1_nint
  use diagnostics
  use boundary_conditions
  use time_step
  use matrix_mod

  implicit none

  integer, intent(in) :: ivel_def   ! populate velocity advance matrices
  integer, intent(in) :: idens_def  ! populate density advance matrices
  integer, intent(in) :: ipres_def  ! populate pressure advance matrices
  integer, intent(in) :: ifield_def ! populate field advance matrices 

  integer :: itri, numelms
  integer :: def_fields

  real :: tstart, tend, tfield, telm, tsizefield, tfinalize
  logical :: is_edge(3)  ! is inode on boundary
  real :: n(2,3)
  integer :: iedge, idim(3)

  tfield = 0.
  telm = 0.
  tsizefield = 0.
  tfinalize = 0.

  numelms = local_elements()

  ! initialize matrices
  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, " initializing matrices..."

  ! Clear matrices
  select case(isplitstep)
  case(0)
     call clear_mat(s1_mat)
     call clear_mat(d1_mat)
     if(eqsubtract.eq.1) call clear_mat(q42_mat)
     if(i3d.eq.1) call clear_mat(o1_mat)
     q4_vec = 0.

  case(1)

     if(ivel_def.eq.1) then
        call clear_mat(s1_mat)
        call clear_mat(d1_mat)
        call clear_mat(q1_mat)
        call clear_mat(r14_mat)
        if(i3d.eq.1) call clear_mat(o1_mat)
        r4_vec = 0.
     end if

     if(ifield_def.eq.1) then
        call clear_mat(s2_mat)
        call clear_mat(d2_mat)
        call clear_mat(r2_mat)
        call clear_mat(q2_mat)
        if(eqsubtract.eq.1) call clear_mat(q42_mat)
        if(i3d.eq.1) call clear_mat(o2_mat)
        q4_vec = 0.
     end if

     if(idens_def.eq.1) then
        call clear_mat(s8_mat)
        call clear_mat(d8_mat)
        call clear_mat(r8_mat)
        call clear_mat(q8_mat)
        qn4_vec = 0.
     endif

     if(ipres_def.eq.1) then
        call clear_mat(s9_mat)
        call clear_mat(d9_mat)
        call clear_mat(r9_mat)
        call clear_mat(q9_mat)
        qp4_vec = 0.
     endif
  end select

  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, " populating matrices..."
    
  ! Specify which fields will be used in matrix population
  def_fields = FIELD_PSI + FIELD_I + FIELD_PHI + FIELD_ETA + FIELD_MU &
             + FIELD_N + FIELD_NI
  if(numvar.ge.2) def_fields = def_fields + FIELD_V

  if(numvar.ge.3) then
     def_fields = def_fields + &
          FIELD_CHI + FIELD_PE + FIELD_B2I + FIELD_J + FIELD_P + FIELD_KAP
  else
     if(ipres.gt.0) then
        def_fields = def_fields + FIELD_P
     endif
  endif

  if(idens_def.eq.1) then
     if(ipellet.eq.1 .or. ionization.eq.1 .or. isink.gt.0) &
          def_fields = def_fields + FIELD_SIG
  endif

  if(gyro.eq.1 .or. amupar.ne.0) then
     if(numvar.lt.3) def_fields = def_fields  + FIELD_PE + FIELD_B2I
  endif

  if(hyperc.ne.0. .and. numvar.ge.3) then
     def_fields = def_fields + FIELD_VOR + FIELD_COM
  endif


  if(integrator.eq.1 .and. ntime.gt.1) then
     bdf = 2.
  else
     bdf = 1.
  endif

  ! Loop over elements
  do itri=1,numelms

     ! calculate the field values and derivatives at the sampling points
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call define_element_quadrature(itri, int_pts_main, int_pts_tor)
     call define_fields(itri, def_fields, 1, linear)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        tfield = tfield + tend - tstart
     endif
     
     ! add element's contribution to matrices
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     if(ivel_def.eq.1) call ludefvel_n(itri)
     if(ifield_def.eq.1) call ludefphi_n(itri)
     if(idens_def.eq.1) call ludefden_n(itri)
     if(ipres_def.eq.1) call ludefpres_n(itri)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        telm = telm + tend - tstart
     endif

     if(isurface.eq.0) cycle
     if(.not.(ifixedb.eq.1 .and. ivform.eq.1)) cycle

     ! add surface terms
     call boundary_edge(itri, is_edge, n, idim)
     
     do iedge=1,3
        if(.not.is_edge(iedge)) cycle

        call define_boundary_quadrature(itri, iedge, 5, n, idim)
        call define_fields(itri, def_fields, 1, linear)

        if(ivel_def.eq.1) call ludefvel_n(itri)
        if(ifield_def.eq.1) call ludefphi_n(itri)
        if(idens_def.eq.1) call ludefden_n(itri)
!        if(ipres_def.eq.1) call ludefpres_n(itri)
     end do
  end do

  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, " finalizing matrices..."

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

  ! since a proc is contributing values to parts of the vector
  ! it does not own, we call sumsharedppplvecvals so that these values
  ! get summed up for all values shared by multiple procs
  ! and then update these values

  ! Finalize matrices for multiplication
  select case(isplitstep)
  case(0)
     call flush(s1_mat)
     call finalize(d1_mat)
     if(eqsubtract.eq.1) call finalize(q42_mat)
     if(i3d.eq.1) call finalize(o1_mat)
     call sum_shared(q4_vec)

  case(1)
     if(ivel_def.eq.1) then
        call flush(s1_mat)
        call finalize(d1_mat)
        call finalize(q1_mat)
        call finalize(r14_mat)
        if(i3d.eq.1) call finalize(o1_mat)
        call sum_shared(r4_vec)
     end if
     
     if(ifield_def.eq.1) then
        call flush(s2_mat)
        call finalize(d2_mat)
        call finalize(r2_mat)
        call finalize(q2_mat)
        if(i3d.eq.1) call finalize(o2_mat)
        if(eqsubtract.eq.1) call finalize(q42_mat)
        call sum_shared(q4_vec)
     end if
     
     if(idens_def.eq.1) then
        call flush(s8_mat)
        call finalize(d8_mat)
        call finalize(q8_mat)
        call finalize(r8_mat)
        call sum_shared(qn4_vec)
     endif ! on idens_def.eq.1
     
     
     if(ipres_def.eq.1) then
        call flush(s9_mat)
        call finalize(d9_mat)
        call finalize(q9_mat)
        call finalize(r9_mat)
        call sum_shared(qp4_vec)
     endif ! on ipres_def.eq.1
  end select

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     tfinalize = tfinalize + tend - tstart
  endif

  if(myrank.eq.0 .and. itimer.eq.1 .and. iprint.ge.1) then
     print *, '  Time spent defining fields: ', tfield
     print *, '  Time spent interpolating size field: ', tsizefield
     print *, '  Time spent calculating elements: ', telm
     print *, '  Time spent finalizing arrays: ', tfinalize
  endif

end subroutine ludefall


!======================================================================
! ludefvel_n
! ----------
!
! populates the matrices for implicit velocity advance
!
! itri: index of finite element
!======================================================================
subroutine ludefvel_n(itri)

  use basic
  use mesh_mod
  use m3dc1_nint
  use arrays
  use sparse
  use time_step

  implicit none

  integer, intent(in) :: itri

  integer :: i, j

  vectype, dimension(dofs_per_element,dofs_per_element,num_fields) :: ss, dd
  vectype, dimension(dofs_per_element,dofs_per_element) :: q_bf
  vectype, dimension(dofs_per_element) :: r4

  type(matrix_type), pointer :: vv1, vv0, vb1, vb0, vn1, vn0, vf0
  type(vector_type), pointer :: vsource
  integer :: advfield
  integer :: pp_i
  integer, dimension(3) :: ieq
  integer :: k
  integer, dimension(dofs_per_element) :: imask

  if(isplitstep.eq.1) then
     vv1 => s1_mat
     vv0 => d1_mat
     vb0 => q1_mat
     vn0 => r14_mat
     vf0 => o1_mat
     vsource => r4_vec
     pp_i = pe_i
  else
     vv1 => s1_mat
     vv0 => d1_mat
     vb1 => s1_mat
     vb0 => d1_mat
     vn1 => s1_mat
     vn0 => d1_mat
     vf0 => o1_mat
     vsource => q4_vec
     pp_i = p_i
  endif

  if(isplitstep.eq.1 .and. iestatic.eq.0) then
     advfield = 1 
  else 
     advfield = 0
  endif

  ieq(1) =   u_i
  ieq(2) =  vz_i
  ieq(3) = chi_i

  do k=1,numvar
     ss = 0.
     dd = 0.
     q_bf = 0.
     r4 = 0.

     select case(k)
     case(1)
        call get_vor_mask(itri, imask)
     case(2)
        call get_vz_mask(itri, imask)
     case(3)
        call get_chi_mask(itri, imask)
     end select

     do i=1,dofs_per_element
        if(imask(i).eq.0) then
           ss(i,:,:) = 0.
           dd(i,:,:) = 0.
           q_bf(i,:) = 0.
           r4(i) = 0.
           cycle
        endif

        do j=1,dofs_per_element
           select case(k)
           case(1)
              call vorticity_lin(mu79(:,:,i),nu79(:,:,j), &
                   ss(i,j,:),dd(i,j,:),q_bf(i,j),advfield)
           case(2)
              call axial_vel_lin(mu79(:,:,i),nu79(:,:,j), &
                   ss(i,j,:),dd(i,j,:),q_bf(i,j),advfield)
           case(3)
              call compression_lin(mu79(:,:,i),nu79(:,:,j), &
                   ss(i,j,:),dd(i,j,:),q_bf(i,j),advfield)             
           end select
        end do

        select case(k)
        case(1)
           call vorticity_nolin(mu79(:,:,i),r4(i))
        case(2)
           call axial_vel_nolin(mu79(:,:,i),r4(i))
        case(3)
           call compression_nolin(mu79(:,:,i),r4(i))
        end select
     end do

     call insert_block(vv1,itri,ieq(k),  u_i,ss(:,:,  u_g),MAT_ADD)
     call insert_block(vv0,itri,ieq(k),  u_i,dd(:,:,  u_g),MAT_ADD)
     call insert_block(vb0,itri,ieq(k),psi_i,dd(:,:,psi_g),MAT_ADD)
     if(numvar.ge.2) then
        call insert_block(vv1,itri,ieq(k), vz_i,ss(:,:,vz_g),MAT_ADD)
        call insert_block(vv0,itri,ieq(k), vz_i,dd(:,:,vz_g),MAT_ADD)
        call insert_block(vb0,itri,ieq(k), bz_i,dd(:,:,bz_g),MAT_ADD)
     endif
     if(numvar.ge.3) then
        call insert_block(vv1,itri,ieq(k),chi_i,ss(:,:,chi_g),MAT_ADD)
        call insert_block(vv0,itri,ieq(k),chi_i,dd(:,:,chi_g),MAT_ADD)
        call insert_block(vb0,itri,ieq(k), pp_i,dd(:,:,  p_g),MAT_ADD)
     endif
     if(idens.eq.1) then
        call insert_block(vn0,itri,ieq(k),den_i,dd(:,:,den_g),MAT_ADD)
     endif
     if(i3d.eq.1) then
        call insert_block(vf0,itri,ieq(k),bf_i,q_bf(:,:),MAT_ADD)
     endif
     if(isplitstep.eq.0) then
        call insert_block(vb1,itri,ieq(k),psi_i,ss(:,:,psi_g),MAT_ADD)
        if(numvar.ge.2) &
             call insert_block(vb1,itri,ieq(k), bz_i,ss(:,:, bz_g),MAT_ADD)
        if(numvar.ge.3) &
             call insert_block(vb1,itri,ieq(k), pp_i,ss(:,:,  p_g),MAT_ADD)
        if(idens.eq.1) &
             call insert_block(vn1,itri,ieq(k),den_i,ss(:,:,den_g),MAT_ADD)
     endif

     call vector_insert_block(vsource,itri,ieq(k),r4,VEC_ADD)
  end do

end subroutine ludefvel_n


!======================================================================
! ludefphi_n
! ----------
!
! populates the matrices for implicit field advance
!
! itri: index of finite element
!======================================================================
subroutine ludefphi_n(itri)
  use basic
  use mesh_mod
  use m3dc1_nint
  use arrays
  use sparse
  use electrostatic_potential
  use vacuum_interface
  use time_step

  implicit none

  integer, intent(in) :: itri

  integer :: i, j, k
  
  vectype, dimension(dofs_per_element,dofs_per_element,num_fields) :: ss, dd
  vectype, dimension(dofs_per_element,dofs_per_element) :: q_bf, r_e, q_e
  vectype, dimension(dofs_per_element) :: q4
  vectype, dimension(dofs_per_element,dofs_per_element,2) :: q_ni

  type(matrix_type), pointer :: bb1, bb0, bv1, bv0, bbf, bni
  type(vector_type), pointer :: bsource
  integer :: ieq(4)
  integer :: maxk
  integer :: imask(dofs_per_element)

  if(isplitstep.eq.1) then
     bb1 => s2_mat
     bb0 => d2_mat
     bv1 => r2_mat
     bv0 => q2_mat
     bbf => o2_mat
     bni => q42_mat
     bsource => q4_vec
  else
     bb1 => s1_mat
     bb0 => d1_mat
     bv1 => s1_mat
     bv0 => d1_mat
     bbf => o1_mat
     bni => q42_mat
     bsource => q4_vec
  endif

  ieq(1) = psi_i
  ieq(2) =  bz_i
  ieq(3) =  pe_i
  ieq(4) =   e_i

  if(jadv.eq.0 .and. i3d.eq.1) then
     maxk = numvar + 1
  else 
     maxk = numvar
  endif

  imask = 1

  do k=1,maxk
     ss = 0.
     dd = 0.
     q_ni = 0.
     q_bf = 0.
     r_e = 0.
     q_e = 0.
     q4 = 0.

     select case(k)
     case(1)
        call get_flux_mask(itri, imask)
     case(2)
        call get_bz_mask(itri, imask)
     case(3)
        call get_p_mask(itri, imask)
     case(4)
        imask = 1
     end select
     
     do i=1,dofs_per_element
        if(imask(i).eq.0) then
           ss(i,:,:) = 0.
           dd(i,:,:) = 0.
           q_ni(i,:,:) = 0.
           q_bf(i,:) = 0.
           r_e(i,:) = 0.
           q_e(i,:) = 0.
           q4(i) = 0.
           cycle
        endif

        do j=1,dofs_per_element
           select case(k)
           case(1)
              if(.not.surface_int) then
                 call flux_lin(mu79(:,:,i),nu79(:,:,j), &
                      ss(i,j,:),dd(i,j,:),q_ni(i,j,:),q_bf(i,j), &
                      r_e(i,j),q_e(i,j))
              end if
              
           case(2)
              if(.not.surface_int) then
                 call axial_field_lin(mu79(:,:,i),nu79(:,:,j), &
                      ss(i,j,:),dd(i,j,:),q_ni(i,j,:),q_bf(i,j))
              endif

           case(3)
              call electron_pressure_lin(mu79(:,:,i),nu79(:,:,j), &
                   ss(i,j,:),dd(i,j,:),q_ni(i,j,:),q_bf(i,j))
           case(4)
              call potential_lin(mu79(:,:,i),nu79(:,:,j), &
                   ss(i,j,:),dd(i,j,:),q_ni(i,j,1),q_bf(i,j),r_e(i,j))
           end select
        end do

        select case(k)
        case(1)
           call flux_nolin(mu79(:,:,i),q4(i))
        case(2)
           call axial_field_nolin(mu79(:,:,i),q4(i))
        case(3)
           call electron_pressure_nolin(mu79(:,:,i),q4(i))
        case(4)
           q4(i) = 0.
        end select
     end do
    
     call insert_block(bb1,itri,ieq(k),psi_i,ss(:,:,psi_g),MAT_ADD)
     call insert_block(bb0,itri,ieq(k),psi_i,dd(:,:,psi_g),MAT_ADD)
     call insert_block(bv1,itri,ieq(k),  u_i,ss(:,:,  u_g),MAT_ADD)
     call insert_block(bv0,itri,ieq(k),  u_i,dd(:,:,  u_g),MAT_ADD)
     if(numvar.ge.2) then
        call insert_block(bb1,itri,ieq(k), bz_i,ss(:,:,bz_g),MAT_ADD)
        call insert_block(bb0,itri,ieq(k), bz_i,dd(:,:,bz_g),MAT_ADD)
        call insert_block(bv1,itri,ieq(k), vz_i,ss(:,:,vz_g),MAT_ADD)
        call insert_block(bv0,itri,ieq(k), vz_i,dd(:,:,vz_g),MAT_ADD)
     endif
     if(numvar.ge.3) then
        call insert_block(bb1,itri,ieq(k), pe_i,ss(:,:, pe_g),MAT_ADD)
        call insert_block(bb0,itri,ieq(k), pe_i,dd(:,:, pe_g),MAT_ADD)
        call insert_block(bv1,itri,ieq(k),chi_i,ss(:,:,chi_g),MAT_ADD)
        call insert_block(bv0,itri,ieq(k),chi_i,dd(:,:,chi_g),MAT_ADD)
     endif
     if(eqsubtract.eq.1) then
        if(idens.eq.1) then
           call insert_block(bni,itri,ieq(k),1,q_ni(:,:,1),MAT_ADD)
        endif
        call insert_block(bni,itri,ieq(k),2,q_ni(:,:,2),MAT_ADD)        
     endif
     if(i3d.eq.1 .and. numvar.ge.2) then
        call insert_block(bbf,itri,ieq(k),bf_i,q_bf(:,:),MAT_ADD)
     endif

     call vector_insert_block(bsource,itri,ieq(k),q4,VEC_ADD)
  end do
end subroutine ludefphi_n


!======================================================================
! ludefden_n
! ----------
!
! populates the matrices for implicit density advance
!
! itri: index of finite element
!======================================================================
subroutine ludefden_n(itri)

  use basic
  use m3dc1_nint
  use arrays
  use sparse
  use metricterms_new
  use time_step

  implicit none

  integer, intent(in) :: itri

  integer :: i, j
  vectype, dimension(dofs_per_element, dofs_per_element) :: ssterm, ddterm
  vectype, dimension(dofs_per_element, dofs_per_element, 3) :: rrterm, qqterm
  vectype, dimension(dofs_per_element) :: oterm

  vectype :: temp

  type(matrix_type), pointer :: nn1, nn0, nv1, nv0
  type(vector_type), pointer :: nsource
  real :: thimpb
  integer :: imask(dofs_per_element)

  vectype, dimension(MAX_PTS,OP_NUM) :: hp
  hp = hypp*sz79

  if(imp_mod.eq.0) then
     thimpb = thimp
  else
     thimpb = 1.
  endif

  if(isplitstep.eq.1) then
     nn1 => s8_mat
     nn0 => d8_mat
     nv1 => r8_mat
     nv0 => q8_mat
     nsource => qn4_vec
  else
     nn1 => s1_mat
     nn0 => d1_mat
     nv1 => s1_mat
     nv0 => d1_mat
     nsource => q4_vec
  endif

  ssterm = 0.
  ddterm = 0.
  rrterm = 0.
  qqterm = 0.
  oterm = 0.

  call get_den_mask(itri, imask)

  do i=1,dofs_per_element
     if(imask(i).eq.0) then
        ssterm(i,:) = 0.
        ddterm(i,:) = 0.
        rrterm(i,:,:) = 0.
        qqterm(i,:,:) = 0.
        oterm(i) = 0.
        cycle
     endif

     do j=1,dofs_per_element

        ! NUMVAR = 1
        ! ~~~~~~~~~~
        temp = n1n(mu79(:,:,i),nu79(:,:,j))
        ssterm(i,j) = ssterm(i,j) + temp    
        ddterm(i,j) = ddterm(i,j) + temp*bdf
        
        temp = n1ndenm(mu79(:,:,i),nu79(:,:,j),denm,hp) &
             + n1nu   (mu79(:,:,i),nu79(:,:,j),pht79)
        ssterm(i,j) = ssterm(i,j) -     thimp     *dt*temp
        ddterm(i,j) = ddterm(i,j) + (1.-thimp*bdf)*dt*temp

        if(linear.eq.0) then
           temp = n1nu(mu79(:,:,i),n179,nu79(:,:,j))
           rrterm(i,j,1) = rrterm(i,j,1) + thimpb*dt*temp
           qqterm(i,j,1) = qqterm(i,j,1) - thimpb*dt*temp*bdf
        endif

        if(eqsubtract.eq.1) then
           temp = n1nu  (mu79(:,:,i),n079,nu79(:,:,j))
           rrterm(i,j,1) = rrterm(i,j,1) +     thimpb     *dt*temp
           qqterm(i,j,1) = qqterm(i,j,1) + (1.-thimpb*bdf)*dt*temp
        endif

#ifdef USECOMPLEX
        ! NUMVAR = 2
        ! ~~~~~~~~~~
        if(numvar.ge.2) then
           temp = n1nv(mu79(:,:,i),nu79(:,:,j),vzt79)
           ssterm(i,j) = ssterm(i,j) -     thimp     *dt*temp
           ddterm(i,j) = ddterm(i,j) + (1.-thimp*bdf)*dt*temp
           
           if(linear.eq.0) then 
              temp = n1nv(mu79(:,:,i),n179,nu79(:,:,j))
              rrterm(i,j,2) = rrterm(i,j,2) + thimpb*dt*temp
              qqterm(i,j,2) = qqterm(i,j,2) - thimpb*dt*temp*bdf
           endif

           if(eqsubtract.eq.1) then
              temp = n1nv(mu79(:,:,i),n079,nu79(:,:,j))
              rrterm(i,j,2) = rrterm(i,j,2) +     thimpb     *dt*temp
              qqterm(i,j,2) = qqterm(i,j,2) + (1.-thimpb*bdf)*dt*temp
           endif
        endif
#endif

        ! NUMVAR = 3
        ! ~~~~~~~~~~
        if(numvar.ge.3) then
           temp = n1nchi(mu79(:,:,i),nu79(:,:,j),cht79)
           ssterm(i,j) = ssterm(i,j) -     thimp     *dt*temp
           ddterm(i,j) = ddterm(i,j) + (1.-thimp*bdf)*dt*temp
           
           if(linear.eq.0) then
              temp = n1nchi(mu79(:,:,i),n179,nu79(:,:,j))
              rrterm(i,j,3) = rrterm(i,j,3) + thimpb*dt*temp
              qqterm(i,j,3) = qqterm(i,j,3) - thimpb*dt*temp*bdf
           endif

           if(eqsubtract.eq.1) then
              temp = n1nchi(mu79(:,:,i),n079,nu79(:,:,j))
              rrterm(i,j,3) = rrterm(i,j,3) +     thimpb     *dt*temp
              qqterm(i,j,3) = qqterm(i,j,3) + (1.-thimpb*bdf)*dt*temp
           endif
        endif

     enddo                     ! on j
     
     oterm(i) = dt*n1s(mu79(:,:,i),sig79)

     if(eqsubtract.eq.1) then
        oterm(i) = oterm(i) + dt*n1ndenm(mu79(:,:,i),n079,denm,hp)
     endif
  enddo                     ! on i
     
  if(isplitstep.eq.0) rrterm = -rrterm

  call insert_block(nn1,itri,den_i,den_i,ssterm,MAT_ADD)
  call insert_block(nn0,itri,den_i,den_i,ddterm,MAT_ADD)
  call insert_block(nv1,itri,den_i,u_i,rrterm(:,:,1),MAT_ADD)
  call insert_block(nv0,itri,den_i,u_i,qqterm(:,:,1),MAT_ADD)
  if(numvar.ge.2) then
     call insert_block(nv1,itri,den_i,vz_i,rrterm(:,:,2),MAT_ADD)
     call insert_block(nv0,itri,den_i,vz_i,qqterm(:,:,2),MAT_ADD)
  endif
  if(numvar.ge.3) then
     call insert_block(nv1,itri,den_i,chi_i,rrterm(:,:,3),MAT_ADD)
     call insert_block(nv0,itri,den_i,chi_i,qqterm(:,:,3),MAT_ADD)
  endif

  call vector_insert_block(nsource,itri,den_i,oterm,VEC_ADD)
end subroutine ludefden_n



!======================================================================
! ludefpres_n
! -----------
!
! populates the matrices for implicit pressure advance
!
! itri: index of finite element
!======================================================================
subroutine ludefpres_n(itri)

  use basic
  use m3dc1_nint
  use arrays
  use sparse
  use metricterms_new
  use time_step

  implicit none

  integer, intent(in) :: itri

  integer :: i, j
  vectype, dimension(dofs_per_element,dofs_per_element) :: ssterm, ddterm
  vectype, dimension(dofs_per_element,dofs_per_element,3) :: rrterm, qqterm
  vectype, dimension(dofs_per_element) :: oterm


  type(matrix_type), pointer :: pp1, pp0, pv1, pv0
  type(vector_type), pointer :: psource

  vectype :: temp
  real :: thimpb

  vectype, dimension(MAX_PTS,OP_NUM) :: hp
  hp = hypp*sz79

  if(imp_mod.eq.0) then
     thimpb = thimp
  else
     thimpb = 1.
  endif

  if(isplitstep.eq.1) then
     pp1 => s9_mat
     pp0 => d9_mat
     pv1 => r9_mat
     pv0 => q9_mat
     psource => qp4_vec
  else
     pp1 => s1_mat
     pp0 => d1_mat
     pv1 => s1_mat
     pv0 => d1_mat
     psource => q4_vec
  endif

  ssterm = 0.
  ddterm = 0.
  rrterm = 0.
  qqterm = 0.
  
  do i=1,dofs_per_element
     do j=1,dofs_per_element
     
        ! NUMVAR = 1
        ! ~~~~~~~~~~
        temp = int2(mu79(:,:,i),nu79(:,:,j))
        ssterm(i,j) = ssterm(i,j) + temp
        ddterm(i,j) = ddterm(i,j) + temp*bdf
        
        temp = p1pu   (mu79(:,:,i),nu79(:,:,j),pht79)
        ssterm(i,j) = ssterm(i,j) -     thimp     *dt*temp
        ddterm(i,j) = ddterm(i,j) + (1.-thimp*bdf)*dt*temp

        if(linear.eq.0) then 
           temp = p1pu(mu79(:,:,i),p179,nu79(:,:,j))
           rrterm(i,j,1) = rrterm(i,j,1) + thimpb*dt*temp
           qqterm(i,j,1) = qqterm(i,j,1) - thimpb*dt*temp*bdf
        endif

        if(eqsubtract.eq.1) then
           temp = p1pu  (mu79(:,:,i),p079,nu79(:,:,j))
           rrterm(i,j,1) = rrterm(i,j,1) +     thimpb     *dt*temp
           qqterm(i,j,1) = qqterm(i,j,1) + (1.-thimpb*bdf)*dt*temp
        endif

        ! NUMVAR = 2
        ! ~~~~~~~~~~
        if(numvar.ge.2) then
          temp = p1pv   (mu79(:,:,i),nu79(:,:,j),vzt79)
          ssterm(i,j) = ssterm(i,j) -     thimp     *dt*temp
          ddterm(i,j) = ddterm(i,j) + (1.-thimp*bdf)*dt*temp

          if(linear.eq.0) then 
             temp = p1pv(mu79(:,:,i),p179,nu79(:,:,j))
             rrterm(i,j,2) = rrterm(i,j,2) + thimpb*dt*temp
             qqterm(i,j,2) = qqterm(i,j,2) - thimpb*dt*temp*bdf
          endif

          if(eqsubtract.eq.1) then
             temp = p1pv  (mu79(:,:,i),p079,nu79(:,:,j))
             rrterm(i,j,2) = rrterm(i,j,2) +     thimpb     *dt*temp
             qqterm(i,j,2) = qqterm(i,j,2) + (1.-thimpb*bdf)*dt*temp
          endif
        endif

        ! NUMVAR = 3
        ! ~~~~~~~~~~
        if(numvar.ge.3) then
           temp = p1pchi    (mu79(:,:,i),nu79(:,:,j),cht79) &
                + b3pedkappa(mu79(:,:,i),nu79(:,:,j),ni79,kap79,hp) &
                + p1psipsikappar(mu79(:,:,i),pst79,pst79,nu79(:,:,j),ni79,b2i79,kar79) &
                + p1kappax  (mu79(:,:,i),nu79(:,:,j),bzt79,ni79,kax79)
           ssterm(i,j) = ssterm(i,j) -     thimp     *dt*temp
           ddterm(i,j) = ddterm(i,j) + (1.-thimp*bdf)*dt*temp
           
           if(linear.eq.0) then 
              temp = p1pchi(mu79(:,:,i),p179,nu79(:,:,j))
              rrterm(i,j,3) = rrterm(i,j,3) + thimpb*dt*temp
              qqterm(i,j,3) = qqterm(i,j,3) - thimpb*dt*temp*bdf
              
              temp = p1uus  (mu79(:,:,i),nu79(:,:,j),ph179,sig79) &
                   + p1uus  (mu79(:,:,i),ph179,nu79(:,:,j),sig79) &
                   + p1uchis(mu79(:,:,i),nu79(:,:,j),ch179,sig79) 
              rrterm(i,j,1) = rrterm(i,j,1) +     thimpb     *dt*temp
              qqterm(i,j,1) = qqterm(i,j,1) + (.5-thimpb*bdf)*dt*temp

              temp = p1vvs  (mu79(:,:,i),nu79(:,:,j),vz179,sig79) &
                   + p1vvs  (mu79(:,:,i),vz179,nu79(:,:,j),sig79)
              rrterm(i,j,2) = rrterm(i,j,2) +     thimpb     *dt*temp
              qqterm(i,j,2) = qqterm(i,j,2) + (.5-thimpb*bdf)*dt*temp
              
              temp = p1chichis(mu79(:,:,i),nu79(:,:,j),ph179,sig79) &
                   + p1chichis(mu79(:,:,i),ph179,nu79(:,:,j),sig79) &
                   + p1uchis  (mu79(:,:,i),ph079,nu79(:,:,j),sig79) 
              rrterm(i,j,3) = rrterm(i,j,3) +     thimpb     *dt*temp
              qqterm(i,j,3) = qqterm(i,j,3) + (.5-thimpb*bdf)*dt*temp
           endif

           if(eqsubtract.eq.1) then
              temp = p1pchi(mu79(:,:,i),p079,nu79(:,:,j))
              rrterm(i,j,3) = rrterm(i,j,3) +     thimpb     *dt*temp
              qqterm(i,j,3) = qqterm(i,j,3) + (1.-thimpb*bdf)*dt*temp

              temp = p1uus  (mu79(:,:,i),nu79(:,:,j),ph079,sig79) &
                   + p1uus  (mu79(:,:,i),ph079,nu79(:,:,j),sig79) &
                   + p1uchis(mu79(:,:,i),nu79(:,:,j),ch079,sig79) 
              rrterm(i,j,1) = rrterm(i,j,1) +     thimpb     *dt*temp
              qqterm(i,j,1) = qqterm(i,j,1) + (1.-thimpb*bdf)*dt*temp

              temp = p1vvs  (mu79(:,:,i),nu79(:,:,j),vz079,sig79) &
                   + p1vvs  (mu79(:,:,i),vz079,nu79(:,:,j),sig79)
              rrterm(i,j,2) = rrterm(i,j,2) +     thimpb     *dt*temp
              qqterm(i,j,2) = qqterm(i,j,2) + (1.-thimpb*bdf)*dt*temp

              temp = p1chichis(mu79(:,:,i),nu79(:,:,j),ph079,sig79) &
                   + p1chichis(mu79(:,:,i),ph079,nu79(:,:,j),sig79) &
                   + p1uchis  (mu79(:,:,i),ph079,nu79(:,:,j),sig79) 
              rrterm(i,j,3) = rrterm(i,j,3) +     thimpb     *dt*temp
              qqterm(i,j,3) = qqterm(i,j,3) + (1.-thimpb*bdf)*dt*temp
           endif
        endif
     enddo                     ! on j

     oterm(i) = 0.
     if(linear.eq.0) then
        oterm(i) = oterm(i) + dt* &
             (b3psipsieta(mu79(:,:,i),pst79,pst79,eta79))
      
        if(numvar.ge.2) then
           oterm(i) = oterm(i) + dt* &
                (b3bbeta(mu79(:,:,i),bzt79,bzt79,eta79))
        endif
        
        if(numvar.ge.3) then
           oterm(i) = oterm(i) + dt*dbf* &
                (b3pebd(mu79(:,:,i),pet79,bzt79,ni79))
        endif

        if(eqsubtract.eq.1) then
           oterm(i) = oterm(i) + dt* &
                (b3pedkappa(mu79(:,:,i),p079,ni79,kap79,hp) &
                +p1psipsikappar(mu79(:,:,i),ps079,ps079,p079,ni79,b2i79,kar79) &
                +p1kappax  (mu79(:,:,i),pe079,bz079,ni79,kax79) &
                +p1uus    (mu79(:,:,i),ph079,ph079,sig79) &
                +p1vvs    (mu79(:,:,i),vz079,vz079,sig79) &
                +p1chichis(mu79(:,:,i),ch079,ch079,sig79) &
                +p1uchis  (mu79(:,:,i),ph079,ch079,sig79))
           
           oterm(i) = oterm(i) + dt* &
                (p1pu   (mu79(:,:,i),p079,ph079))
           
           if(numvar.ge.3) then
              oterm(i) = oterm(i) + dt* &
                   (p1pchi(mu79(:,:,i),p079,ch079))
           endif
        endif
        
     endif  ! on linear.eq.0
  enddo                     ! on i

  if(isplitstep.eq.0) rrterm = -rrterm

  call insert_block(pp1,itri,den_i,den_i,ssterm,MAT_ADD)
  call insert_block(pp0,itri,den_i,den_i,ddterm,MAT_ADD)
  call insert_block(pv1,itri,den_i,u_i,rrterm(:,:,1),MAT_ADD)
  call insert_block(pv0,itri,den_i,u_i,qqterm(:,:,1),MAT_ADD)
  if(numvar.ge.2) then
     call insert_block(pv1,itri,den_i,vz_i,rrterm(:,:,2),MAT_ADD)
     call insert_block(pv0,itri,den_i,vz_i,qqterm(:,:,2),MAT_ADD)
  endif
  if(numvar.ge.3) then
     call insert_block(pv1,itri,den_i,chi_i,rrterm(:,:,3),MAT_ADD)
     call insert_block(pv0,itri,den_i,chi_i,qqterm(:,:,3),MAT_ADD)
  endif

  call vector_insert_block(psource,itri,den_i,oterm,VEC_ADD)

end subroutine ludefpres_n
