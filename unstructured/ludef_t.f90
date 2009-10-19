!======================================================================
! Vorticity Equation
!======================================================================
subroutine vorticity_lin(trial, lin, ssterm, ddterm, q_bf, advfield)
  
  use basic
  use arrays
  use nintegrate_mod
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

  ! NUMVAR = 1
  ! ~~~~~~~~~~
  temp = v1un(trial,lin,nt79)
  ssterm(u_g) = ssterm(u_g) + temp
  ddterm(u_g) = ddterm(u_g) + temp*bdf
  
  temp = v1umu(trial,lin,vis79,vic79) &
       + v1us (trial,lin,sig79)
  ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
  ddterm(u_g) = ddterm(u_g) + (1.-thimp*bdf)*dt*temp

  temp = v1uun(trial,lin,ph179,nt79) &
       + v1uun(trial,ph179,lin,nt79)
  ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
  ddterm(u_g) = ddterm(u_g) + (.5-thimp*bdf)*dt*temp

  if(idens.eq.1) then
     ddterm(den_g) = ddterm(den_g) + dt* &
          v1ngrav(trial,lin)
  endif
  
  if(gyro.eq.1) then
     temp = g1u(trial,lin)*dbf
     ssterm(u_g) = ssterm(u_g) +     thimp     *dt*temp
     ddterm(u_g) = ddterm(u_g) - (1.-thimp*bdf)*dt*temp
  endif

  if(advfield.eq.1) then
     temp = v1upsipsi(trial,lin,pst79,pst79) &
          + v1up     (trial,lin,pt79) &
          + v1ungrav (trial,lin,nt79)
     ssterm(u_g) = ssterm(u_g) - thimp*thimp*dt*dt*temp
     ddterm(u_g) = ddterm(u_g) +       ththm*dt*dt*temp

     ddterm(psi_g) = ddterm(psi_g)  + dt* &
          (v1psipsi(trial,lin,pss79)  &
          +v1psipsi(trial,pss79,lin))

     if(numvar.ge.3 .or. ipres.eq.1) then
        ddterm(p_g) = ddterm(p_g) + dt*    &
             v1p(trial,lin)
     endif
  
     if(isources.eq.1) then
        ddterm(psi_g) = ddterm(psi_g) + thimp*dt*dt* &
             (v1psisb1 (trial,lin,sb179))
     endif
  else
     temp = v1psipsi(trial,lin,ps179) &
          + v1psipsi(trial,ps179,lin)
     ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (.5-thimp*bdf)*dt*temp

     if(idens.eq.1) then
        temp = v1ngrav(trial,lin)
        ssterm(den_g) = ssterm(den_g) -     thimp     *dt*temp
        ddterm(den_g) = ddterm(den_g) + (1.-thimp*bdf)*dt*temp
     endif

     if(numvar.ge.3 .or. ipres.eq.1) then
        temp = v1p(trial,lin)
        ssterm(p_g) = ssterm(p_g) -     thimp     *dt*temp
        ddterm(p_g) = ddterm(p_g) + (1.-thimp*bdf)*dt*temp
     endif
  endif
    
  if(eqsubtract.eq.1) then
     temp = v1uun(trial,lin,ph079,nt79) &    
          + v1uun(trial,ph079,lin,nt79)
     ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
     ddterm(u_g) = ddterm(u_g) + (1.-thimp*bdf)*dt*temp

     if(idens.eq.1) then
        ddterm(den_g) = ddterm(den_g) + dt* &
             v1uun(trial,ph079,ph079,lin)
     endif
     
     if(advfield.eq.1) then
        ddterm(psi_g) = ddterm(psi_g) + thimp*dt*dt* &
             (v1upsipsi(trial,ph079,lin,pss79) &
             +v1upsipsi(trial,ph079,pss79,lin))

        if(idens.eq.1) then
           ddterm(den_g) = ddterm(den_g) + thimp*dt*dt*  &
                (v1ungrav   (trial,ph079,lin)) !    &
!               +v1ndenmgrav(trial,lin, denm))
        endif

        if(numvar.ge.3 .or. ipres.eq.1) then
           ddterm(p_g) = ddterm(p_g) + thimp*dt*dt* &
                v1up(trial,ph079,lin)
        endif
     else
        temp = v1psipsi(trial,lin,ps079) &
             + v1psipsi(trial,ps079,lin)
        ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
        ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp
     endif
  endif
  
  
  ! NUMVAR = 2
  ! ~~~~~~~~~~
  if(numvar.eq.1) then
     ! These terms exist since F0=bzero when numvar=1
     if(advfield.eq.1) then
        temp = v1upsib(trial,lin,pst79,bzt79) &
             + v1ubb  (trial,lin,bzt79,bzt79)
        ssterm(u_g) = ssterm(u_g) - thimp*thimp*dt*dt*temp
        ddterm(u_g) = ddterm(u_g) +       ththm*dt*dt*temp
        
        ddterm(psi_g) = ddterm(psi_g) + dt* &
             (v1psib(trial,lin,bzt79))
     else
        temp = v1psib(trial,lin,bzt79)
        ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
        ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp
     endif

  else       ! on numvar.eq.1
     temp = v1vmu(trial,lin,vis79,vic79)
     ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
     ddterm(vz_g) = ddterm(vz_g) + (1.-thimp*bdf)*dt*temp

     temp = v1uvn(trial,lin,vz179,nt79)
     ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
     ddterm(u_g) = ddterm(u_g) + (.5-thimp*bdf)*dt*temp

     temp = v1vvn(trial,lin,vz179,nt79) &
          + v1vvn(trial,vz179,lin,nt79) &
          + v1uvn(trial,ph179,lin,nt79)
     ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
     ddterm(vz_g) = ddterm(vz_g) + (.5-thimp*bdf)*dt*temp
     
     if(gyro.eq.1) then
        temp = g1v(trial,lin)*dbf
        ssterm(vz_g) = ssterm(vz_g) +     thimp     *dt*temp
        ddterm(vz_g) = ddterm(vz_g) - (1.-thimp*bdf)*dt*temp
     endif

     if(advfield.eq.1) then
        temp = v1upsib(trial,lin,pst79,bzt79) &
             + v1ubb  (trial,lin,bzt79,bzt79)
        ssterm(u_g) = ssterm(u_g) - thimp*thimp*dt*dt*temp
        ddterm(u_g) = ddterm(u_g) +       ththm*dt*dt*temp

        temp = v1vpsipsi(trial,lin,pst79,pst79) &
             + v1vpsib  (trial,lin,pst79,bzt79) &
             + v1vp     (trial,lin,pt79)
        ssterm(vz_g) = ssterm(vz_g) - thimp*thimp*dt*dt*temp
        ddterm(vz_g) = ddterm(vz_g) +       ththm*dt*dt*temp

        ddterm(psi_g) = ddterm(psi_g) + dt* &
             (v1psib(trial,lin,bzs79))

        ddterm(bz_g) = ddterm(bz_g) + dt* &
             (v1psib(trial,pss79,lin) &
             +v1bb  (trial,lin,bzs79) &
             +v1bb  (trial,bzs79,lin))

        if(isources.eq.1) then
           ddterm(bz_g) = ddterm(bz_g) + thimp*dt*dt* &
                (v1bsb2(trial,lin,sb279))
        endif
     else
        temp = v1psib(trial,lin,bz179)
        ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
        ddterm(psi_g) = ddterm(psi_g) + (.5-thimp*bdf)*dt*temp

        temp = v1psib(trial,ps179,lin) &
             + v1bb  (trial,lin,bz179) &
             + v1bb  (trial,bz179,lin)
        ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
        ddterm(bz_g) = ddterm(bz_g) + (.5-thimp*bdf)*dt*temp
     endif
          
     if(eqsubtract.eq.1) then
        temp = v1uvn(trial,lin,vz079,nt79)
        ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
        ddterm(u_g) = ddterm(u_g) + (1.-thimp*bdf)*dt*temp

        temp = v1vvn(trial,lin,vz079,nt79) &
             + v1vvn(trial,vz079,lin,nt79) &
             + v1uvn(trial,ph079,lin,nt79)
        ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
        ddterm(vz_g) = ddterm(vz_g) + (1.-thimp*bdf)*dt*temp

        if(idens.eq.1) then
           ddterm(den_g) = ddterm(den_g) + dt*       &
                (v1uvn(trial,ph079,vz079,lin) &
                +v1vvn(trial,vz079,vz079,lin))
        endif

        if(advfield.eq.1) then
           ddterm(psi_g) = ddterm(psi_g) + thimp*dt*dt* &
                (v1upsib  (trial,ph079,lin,bzs79) &
                +v1vpsipsi(trial,vz079,lin,pss79) &
                +v1vpsipsi(trial,vz079,pss79,lin) &
                +v1vpsib  (trial,vz079,lin,bzs79))
        
           ddterm(bz_g) = ddterm(bz_g) + thimp*dt*dt* &
                (v1upsib(trial,ph079,pss79,lin) &
                +v1ubb  (trial,ph079,lin,bzs79) &
                +v1ubb  (trial,ph079,bzs79,lin) &
                +v1vpsib(trial,vz079,pss79,lin))

           if(numvar.ge.3. .or. ipres.eq.1) then
              ddterm(p_g) = ddterm(p_g) +  thimp*dt*dt* &
                   v1vp  (trial,vz079,lin)
           endif
        else
           temp = v1psib(trial,lin,bz079)
           ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
           ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

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
     
     temp = v1uchin(trial,lin,ch179,nt79)       
     ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
     ddterm(u_g) = ddterm(u_g) + (.5-thimp*bdf)*dt*temp

     temp = v1vchin(trial,lin,ch179,nt79)       
     ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
     ddterm(vz_g) = ddterm(vz_g) + (.5-thimp*bdf)*dt*temp
     
     temp = v1chin(trial,lin,nt79)*chiiner
     ssterm(chi_g) = ssterm(chi_g) + temp
     ddterm(chi_g) = ddterm(chi_g) + temp*bdf
     
     temp = v1chimu(trial,lin,vis79,vic79) &
          + v1chis (trial,lin,sig79)
     ssterm(chi_g) = ssterm(chi_g) -     thimp     *dt*temp
     ddterm(chi_g) = ddterm(chi_g) + (1.-thimp*bdf)*dt*temp
     
     temp = v1uchin  (trial,ph179,lin,nt79) &
          + v1vchin  (trial,vz179,lin,nt79) &
          + v1chichin(trial,lin,ch179,nt79) &
          + v1chichin(trial,ch179,lin,nt79)
     ssterm(chi_g) = ssterm(chi_g) -     thimp     *dt*temp
     ddterm(chi_g) = ddterm(chi_g) + (.5-thimp*bdf)*dt*temp
          
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
        temp = v1uchin(trial,lin,ch079,nt79)
        ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
        ddterm(u_g) = ddterm(u_g) + (1.-thimp*bdf)*dt*temp

        temp = v1vchin(trial,lin,ch079,nt79)
        ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
        ddterm(vz_g) = ddterm(vz_g) + (1.-thimp*bdf)*dt*temp
        
        temp = v1uchin  (trial,ph079,lin,nt79) &
             + v1vchin  (trial,vz079,lin,nt79) &
             + v1chichin(trial,lin,ch079,nt79) &
             + v1chichin(trial,ch079,lin,nt79)
        ssterm(chi_g) = ssterm(chi_g) -     thimp     *dt*temp
        ddterm(chi_g) = ddterm(chi_g) + (1.-thimp*bdf)*dt*temp

        if(idens.eq.1) then
           ddterm(den_g) = ddterm(den_g) + dt* &
                (v1uchin  (trial,vz079,ch079,lin) &
                +v1vchin  (trial,vz079,ch079,lin) &
                +v1chichin(trial,ch079,ch079,lin))
        endif
        
        if(advfield.eq.1) then
           ddterm(psi_g) = ddterm(psi_g) + thimp*dt*dt* &
                (v1chipsipsi(trial,ch079,lin,pss79) &
                +v1chipsipsi(trial,ch079,pss79,lin) &
                +v1chipsib  (trial,ch079,lin,bzs79))
        
           ddterm(bz_g) = ddterm(bz_g) + thimp*dt*dt* &
                (v1chipsib(trial,ch079,pss79,lin) &
                +v1chibb  (trial,ch079,lin,bzs79) &
                +v1chibb  (trial,ch079,bzs79,lin))

           ddterm(p_g) = ddterm(p_g) +  thimp*dt*dt* &
                v1chip(trial,ch079,lin)

           if(idens.eq.1) then
              ddterm(den_g) = ddterm(den_g) + thimp*dt*dt* &
                   v1chingrav(trial,ch079,lin)
           endif
        else
        endif
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
  use nintegrate_mod
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
  use nintegrate_mod
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

  temp = v2umu(trial,lin,vis79,vic79)
  ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
  ddterm(u_g) = ddterm(u_g) + (1.-thimp*bdf)*dt*temp
         
  temp = v2vun(trial,vz179,lin,nt79)
  ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
  ddterm(u_g) = ddterm(u_g) + (.5-thimp*bdf)*dt*temp
          
  temp = v2vn(trial,lin,nt79)
  ssterm(vz_g) = ssterm(vz_g) + temp
  ddterm(vz_g) = ddterm(vz_g) + temp*bdf
     
  temp = v2vmu  (trial,lin,vis79,vic79,hv) &
       + v2vs   (trial,lin,sig79)
  ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
  ddterm(vz_g) = ddterm(vz_g) + (1.-thimp*bdf)*dt*temp

  temp = v2vun(trial,lin,ph179,nt79) &
       + v2vvn(trial,lin,vz179,nt79) &
       + v2vvn(trial,vz179,lin,nt79)
  ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
  ddterm(vz_g) = ddterm(vz_g) + (.5-thimp*bdf)*dt*temp

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
     
     if(isources.eq.1) then
        ddterm(psi_g) = ddterm(psi_g) + thimp*dt*dt* &
             (v2psisb2(trial,lin,sb279))
        ddterm(bz_g) = ddterm(bz_g) + thimp*dt*dt* &
             (v2bsb1(trial,lin,sb179))
     endif
  else
     temp = v2psipsi(trial,lin,ps179) &
          + v2psipsi(trial,ps179,lin) &
          + v2psib  (trial,lin,bz179)
     ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (.5-thimp*bdf)*dt*temp

     temp = v2psib(trial,ps179,lin)
     ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
     ddterm(bz_g) = ddterm(bz_g) + (.5-thimp*bdf)*dt*temp

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
              
     temp = v2vun(trial,lin,ph079,nt79) &
          + v2vvn(trial,lin,vz079,nt79) &
          + v2vvn(trial,vz079,lin,nt79)
     ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
     ddterm(vz_g) = ddterm(vz_g) + (1.-thimp*bdf)*dt*temp

     if(advfield.eq.1) then
        ddterm(psi_g) = ddterm(psi_g) + thimp*dt*dt* &
             (v2upsipsi(trial,ph079,lin,pss79) &
             +v2upsipsi(trial,ph079,pss79,lin) &
             +v2upsib  (trial,ph079,lin,bzs79) &
             +v2vpsipsi(trial,vz079,lin,pss79) &
             +v2vpsipsi(trial,vz079,pss79,lin) &
             +v2vpsib  (trial,vz079,lin,bzs79))
        
        ddterm(bz_g) = ddterm(bz_g) + thimp*dt*dt* &
             (v2upsib(trial,ph079,pss79,lin) &
             +v2ubb  (trial,ph079,bzs79,lin) &
             +v2ubb  (trial,ph079,lin,bzs79) &
             +v2vpsib(trial,vz079,pss79,lin))

        if(numvar.ge.3 .or. ipres.eq.1) then
           ddterm(p_g) = ddterm(p_g) + thimp*dt*dt* &
                (v2up  (trial,ph079,lin) &
                +v2vp  (trial,vz079,lin))
        endif
     else
        temp = v2psipsi(trial,lin,ps079) &
             + v2psipsi(trial,ps079,lin) &
             + v2psib  (trial,lin,bz079)
        ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
        ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

        temp = v2psib(trial,ps079,lin)
        ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
        ddterm(bz_g) = ddterm(bz_g) + (1.-thimp*bdf)*dt*temp

        if(idens.eq.1) then
           ddterm(den_g) = ddterm(den_g) + dt* &
                (v2vun(trial,vz079,ph079,lin) &
                +v2vvn(trial,vz079,vz079,lin))
        endif
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
             +v2bf  (trial,bz079,lin))
     endif
  endif

        
  ! NUMVAR = 3
  ! ~~~~~~~~~~
  if(numvar.ge.3) then

     temp = v2vchin(trial,lin,ch179,nt79)
     ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
     ddterm(vz_g) = ddterm(vz_g) + (.5-thimp*bdf)*dt*temp

     temp = v2chimu  (trial,lin,vis79,vic79)
     ssterm(chi_g) = ssterm(chi_g) -     thimp     *dt*temp
     ddterm(chi_g) = ddterm(chi_g) + (1.-thimp*bdf)*dt*temp
           
     temp = v2vchin(trial,vz179,lin,nt79)
     ssterm(chi_g) = ssterm(chi_g) -     thimp     *dt*temp
     ddterm(chi_g) = ddterm(chi_g) + (.5-thimp*bdf)*dt*temp

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
        temp = v2vchin(trial,lin,ch079,nt79)
        ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
        ddterm(vz_g) = ddterm(vz_g) + (1.-thimp*bdf)*dt*temp
              
        temp = v2vchin(trial,vz079,lin,nt79)
        ssterm(chi_g) = ssterm(chi_g) -     thimp     *dt*temp
        ddterm(chi_g) = ddterm(chi_g) + (1.-thimp*bdf)*dt*temp

        if(advfield.eq.1) then
           ddterm(psi_g) = ddterm(psi_g) + thimp*dt*dt* &
                (v2chipsipsi(trial,ch079,lin,pss79) &
                +v2chipsipsi(trial,ch079,pss79,lin) &
                +v2chipsib  (trial,ch079,lin,bzs79))
              
           ddterm(bz_g) = ddterm(bz_g) + thimp*dt*dt* &
                (v2chipsib(trial,ch079,pss79,lin) &
                +v2chibb  (trial,ch079,lin,bzs79) &
                +v2chibb  (trial,ch079,bzs79,lin))

           ddterm(p_g) = ddterm(p_g) + thimp*dt*dt* &
                v2chip(trial,ch079,lin)

        else
           if(idens.eq.1) then
              ddterm(den_g) = ddterm(den_g) + dt* &
                   v2vchin(trial,vz079,ch079,lin)
           endif
        endif
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
  use nintegrate_mod
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
  use nintegrate_mod
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

  temp = v3uun  (trial,lin,ph179,nt79) &
       + v3uun  (trial,ph179,lin,nt79) &
       + v3uvn  (trial,lin,vz179,nt79) &
       + v3uchin(trial,lin,ch179,nt79)
  ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
  ddterm(u_g) = ddterm(u_g) + (.5-thimp*bdf)*dt*temp
           
  temp = v3umu(trial,lin,vis79,vic79) &
       + v3us (trial,lin,sig79)
  ssterm(u_g) = ssterm(u_g) -     thimp     *dt*temp
  ddterm(u_g) = ddterm(u_g) + (1.-thimp*bdf)*dt*temp

  temp = v3vmu(trial,lin,vis79,vic79)
  ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
  ddterm(vz_g) = ddterm(vz_g) + (1.-thimp*bdf)*dt*temp
                     
  temp = v3uvn  (trial,ph179,lin,nt79) &
       + v3vvn  (trial,lin,vz179,nt79) &
       + v3vvn  (trial,vz179,lin,nt79) &
       + v3vchin(trial,lin,ch179,nt79)
  ssterm(vz_g) = ssterm(vz_g) -     thimp     *dt*temp
  ddterm(vz_g) = ddterm(vz_g) + (.5-thimp*bdf)*dt*temp

  temp = v3chimu(trial,lin,vis79,vic79) &
       + v3chis (trial,lin,sig79)
  ssterm(chi_g) = ssterm(chi_g) -     thimp     *dt*temp
  ddterm(chi_g) = ddterm(chi_g) + (1.-thimp*bdf)*dt*temp
  
  temp = v3uchin  (trial,ph179,lin,nt79) &
       + v3vchin  (trial,vz179,lin,nt79) &
       + v3chichin(trial,lin,ch179,nt79) &
       + v3chichin(trial,ch179,lin,nt79)  
  ssterm(chi_g) = ssterm(chi_g) -     thimp     *dt*temp
  ddterm(chi_g) = ddterm(chi_g) + (.5-thimp*bdf)*dt*temp

  if(idens.eq.1) then
     ddterm(den_g) = ddterm(den_g) + dt* &
          v3ngrav(trial,lin)
  endif
    
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

     if(isources.eq.1) then
        ddterm(psi_g) = ddterm(psi_g) +    &
             thimp*dt*dt*                  &
             (v3psisb1(trial,lin,sb179))
        
        ddterm(bz_g) = ddterm(bz_g) +      &
             thimp*dt*dt*                  &
             (v3bsb2(trial,lin,sb279))
     end if
  else
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

     temp = v3p(trial,lin) 
     ssterm(p_g) = ssterm(p_g) -     thimp     *dt*temp
     ddterm(p_g) = ddterm(p_g) + (1.-thimp*bdf)*dt*temp

     if(idens.eq.1) then
        temp = v3ngrav(trial,lin)
        ssterm(den_g) = ssterm(den_g) -     thimp     *dt*temp
        ddterm(den_g) = ddterm(den_g) + (1.-thimp*bdf)*dt*temp
     endif
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

     if(idens.eq.1) then
        ddterm(den_g) = ddterm(den_g) + dt* &
             (v3uun    (trial,ph079,ph079,lin) &
             +v3uchin  (trial,ph079,ch079,lin) &
             +v3vvn    (trial,vz079,vz079,lin) &
             +v3chichin(trial,ch079,ch079,lin))
     endif

     if(advfield.eq.1) then
        ddterm(psi_g) = ddterm(psi_g) + thimp*dt*dt* &
             (v3upsipsi  (trial,ph079,lin,pss79) &
             +v3upsipsi  (trial,ph079,pss79,lin) &
             +v3upsib    (trial,ph079,lin,bzs79) &
             +v3vpsipsi  (trial,vz079,lin,pss79) &
             +v3vpsipsi  (trial,vz079,pss79,lin) &
             +v3vpsib    (trial,vz079,lin,bzs79) &
             +v3chipsipsi(trial,ch079,lin,pss79) &
             +v3chipsipsi(trial,ch079,pss79,lin) &
             +v3chipsib  (trial,ch079,lin,bzs79))
     
        ddterm(bz_g) = ddterm(bz_g) + thimp*dt*dt* &
             (v3upsib  (trial,ph079,pss79,lin) &
             +v3ubb    (trial,ph079,lin,bzs79) &
             +v3ubb    (trial,ph079,bzs79,lin) &
             +v3vpsib  (trial,vz079,pss79,lin) &
             +v3vbb    (trial,vz079,lin,bzs79) &
             +v3vbb    (trial,vz079,bzs79,lin) &
             +v3chipsib(trial,ch079,pss79,lin) &
             +v3chibb  (trial,ch079,lin,bzs79) &
             +v3chibb  (trial,ch079,bzs79,lin))
     
        ddterm(p_g) = ddterm(p_g) + thimp*dt*dt* &
             (v3up  (trial,ph079,lin) &
             +v3vp  (trial,vz079,lin) &
             +v3chip(trial,ch079,lin))
        
        if(idens.eq.1) then
           ddterm(den_g) = ddterm(den_g) + thimp*dt*dt* &
                (v3ungrav   (trial,ph079,lin) &
                +v3chingrav (trial,ch079,lin)) ! &
!                +v3ndenmgrav(trial,n179, denm))
        endif

     else
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

        if(idens.eq.1) then
           ddterm(den_g) = ddterm(den_g) + dt*             &
                (v3uun    (trial,ph079,ph079,lin)  &
                +v3vvn    (trial,vz079,vz079,lin)  &
                +v3uchin  (trial,ph079,ch079,lin)  &
                +v3chichin(trial,ch079,ch079,lin))
        end if
     end if
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
        q_bf = q_bf + dt* &
             (v3psif(trial,ps079,lin) &
             +v3bf  (trial,bz079,lin))
     endif
  endif

end subroutine compression_lin

subroutine compression_nolin(trial, r4term)

  use basic
  use nintegrate_mod
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
subroutine flux_lin(trial, lin, ssterm, ddterm, q_ni, q_bf)
  
  use basic
  use arrays
  use nintegrate_mod
  use metricterms_new

  implicit none

  vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: trial, lin 
  vectype, dimension(num_fields), intent(out) :: ssterm, ddterm
  vectype, intent(out) :: q_ni, q_bf
  vectype :: temp
  real :: thimpb

  vectype, dimension(MAX_PTS, OP_NUM) :: hf
  hf = hypf*sz79


  if(imp_mod.eq.0) then
     thimpb = thimp
  else
     thimpb = 1.
  endif

  ssterm = 0.
  ddterm = 0.
  q_ni = 0.
  q_bf = 0.

  temp = b1psi (trial,lin) &
       - b1psid(trial,lin,ni79)   ! electron mass term
  ssterm(psi_g) = ssterm(psi_g) + temp
  ddterm(psi_g) = ddterm(psi_g) + temp*bdf

  temp = b1psieta(trial,lin,eta79,hf)
  ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
  ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

  temp = b1psiu(trial,lin,pht79)
  ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
  ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

  temp = (b1psipsid(trial,lin,ps179,ni79) &
         +b1psipsid(trial,ps179,lin,ni79))*dbf
  ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
  ddterm(psi_g) = ddterm(psi_g) + (.5-thimp*bdf)*dt*temp

  temp = b1psiu(trial,ps179,lin)
  ssterm(u_g) = ssterm(u_g) - thimpb*dt*temp
  ddterm(u_g) = ddterm(u_g) - thimpb*dt*temp*bdf

  if(eqsubtract.eq.1) then
     temp = (b1psipsid(trial,lin,ps079,ni79) &
            +b1psipsid(trial,ps079,lin,ni79))*dbf
     ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

     temp = b1psiu(trial,ps079,lin)
     ssterm(u_g) = ssterm(u_g) -     thimpb     *dt*temp
     ddterm(u_g) = ddterm(u_g) + (1.-thimpb*bdf)*dt*temp
  endif

  ! NUMVAR = 2
  ! ~~~~~~~~~~
  if(numvar.eq.1) then
     ! Term due to F0 = bzero when numvar=1
     temp = b1psibd(trial,lin,bzt79,ni79)*dbf
     ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

     temp = b1bu(trial,bzt79,lin)
     ssterm(u_g) = ssterm(u_g) -     thimpb     *dt*temp
     ddterm(u_g) = ddterm(u_g) + (1.-thimpb*bdf)*dt*temp
  else
     temp = b1psibd(trial,lin,bz179,ni79)*dbf
     ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (.5-thimp*bdf)*dt*temp

     temp = b1beta(trial,lin,eta79)
     ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
     ddterm(bz_g) = ddterm(bz_g) + (1.-thimp*bdf)*dt*temp

     temp = b1psiv  (trial,lin,vzt79)
     ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

     temp = b1bu(trial,lin,pht79) &
          + b1bv(trial,lin,vzt79)
     ssterm(bz_g) = ssterm(bz_g) -     thimpb     *dt*temp
     ddterm(bz_g) = ddterm(bz_g) + (1.-thimpb*bdf)*dt*temp
           
     temp = b1psibd(trial,ps179,lin,ni79)*dbf &
          + b1bbd  (trial,bz179,lin,ni79)*dbf &
          + b1bbd  (trial,lin,bz179,ni79)*dbf
     ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
     ddterm(bz_g) = ddterm(bz_g) + (.5-thimp*bdf)*dt*temp

     temp = b1psiv(trial,ps179,lin)
     ssterm(vz_g) = ssterm(vz_g) - thimpb*dt*temp
     ddterm(vz_g) = ddterm(vz_g) - thimpb*dt*temp*bdf

     temp = b1bu(trial,bz179,lin)
     ssterm(u_g) = ssterm(u_g) - thimpb*dt*temp
     ddterm(u_g) = ddterm(u_g) - thimpb*dt*temp*bdf

     temp = b1bv(trial,bz179,lin)
     ssterm(vz_g) = ssterm(vz_g) - thimpb*dt*temp
     ddterm(vz_g) = ddterm(vz_g) - thimpb*dt*temp*bdf
           
     if(eqsubtract.eq.1) then
        temp = b1psibd(trial,lin,bz079,ni79)*dbf
        ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
        ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp
              
        temp = b1psibd(trial,ps079,lin,ni79)*dbf &
             + b1bbd  (trial,bz079,lin,ni79)*dbf &
             + b1bbd  (trial,lin,bz079,ni79)*dbf
        ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
        ddterm(bz_g) = ddterm(bz_g) + (1.-thimp*bdf)*dt*temp

        temp = b1psiv(trial,ps079,lin)
        ssterm(vz_g) = ssterm(vz_g) -     thimpb     *dt*temp
        ddterm(vz_g) = ddterm(vz_g) + (1.-thimpb*bdf)*dt*temp

        temp = b1bu(trial,bz079,lin)
        ssterm(u_g) = ssterm(u_g) -     thimpb     *dt*temp
        ddterm(u_g) = ddterm(u_g) + (1.-thimpb*bdf)*dt*temp

        temp = b1bv(trial,bz079,lin)
        ssterm(vz_g) = ssterm(vz_g) -     thimpb     *dt*temp
        ddterm(vz_g) = ddterm(vz_g) + (1.-thimpb*bdf)*dt*temp
     end if
  end if

  if(i3d.eq.1 .and. numvar.ge.2) then
     temp = b1psifd(trial,lin,bft79,ni79)*dbf
     ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

     temp = b1bfd(trial,lin,bft79,ni79)*dbf
     ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
     ddterm(bz_g) = ddterm(bz_g) + (1.-thimp*bdf)*dt*temp

     if(eqsubtract.eq.1) then
        q_bf = q_bf + dt* &
             (b1psifd(trial,ps079,lin,ni79)*dbf &
             +b1bfd  (trial,bz079,lin,ni79)*dbf)
     endif
  endif


  ! NUMVAR = 3
  ! ~~~~~~~~~~
  if(numvar.ge.3) then
     temp = b1ped(trial,lin,ni79)*dbf*pefac
     ssterm(pe_g) = ssterm(pe_g) -     thimp     *dt*temp
     ddterm(pe_g) = ddterm(pe_g) + (1.-thimp*bdf)*dt*temp

     temp = b1psichi(trial,lin,cht79)                
     ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

     temp = b1psichi(trial,ps179,lin)                
     ssterm(chi_g) = ssterm(chi_g) - thimpb*dt*temp
     ddterm(chi_g) = ddterm(chi_g) - thimpb*dt*temp*bdf

     temp = b1bchi(trial,bz179,lin)
     ssterm(chi_g) = ssterm(chi_g) - thimpb*dt*temp
     ddterm(chi_g) = ddterm(chi_g) - thimpb*dt*temp*bdf
     
     if(eqsubtract.eq.1) then
        temp = b1psichi(trial,ps079,lin)
        ssterm(chi_g) = ssterm(chi_g) -     thimpb     *dt*temp
        ddterm(chi_g) = ddterm(chi_g) + (1.-thimpb*bdf)*dt*temp

        temp = b1bchi(trial,bz079,lin)
        ssterm(chi_g) = ssterm(chi_g) -     thimpb     *dt*temp
        ddterm(chi_g) = ddterm(chi_g) + (1.-thimpb*bdf)*dt*temp
     endif
  end if

  ! density terms
  ! ~~~~~~~~~~~~~
  if(idens.eq.1 .and. eqsubtract.eq.1) then
     q_ni = q_ni + dt* &
          (b1psipsid(trial,ps079,ps079,lin)*dbf)

     if(numvar.eq.1) then
        ! Term due to F0 = bzero when numvar=1
        q_ni = q_ni + dt* &
             (b1psibd(trial,ps079,bzt79,lin)*dbf)
     endif

     if(numvar.ge.2) then
        q_ni = q_ni + dt* &
             (b1psibd(trial,ps079,bz079,lin)*dbf &
             +b1bbd  (trial,bz079,bz079,lin)*dbf)
     endif
     if(numvar.ge.3) then
        q_ni = q_ni + dt* &
             (b1ped(trial,pe079,lin)*dbf)
     endif
  endif
        
end subroutine flux_lin

subroutine flux_nolin(trial, r4term)

  use basic
  use nintegrate_mod
  use metricterms_new

  implicit none

  vectype, intent(in), dimension(MAX_PTS, OP_NUM)  :: trial
  vectype, intent(out) :: r4term
  
  r4term = 0.

  if(igauge.eq.1 .or. linear.eq.1) then
     r4term = r4term - dt* &
          vloop*int1(trial)/(2.*pi)
  endif
 
end subroutine flux_nolin


!======================================================================
! Axial Field Equation
!======================================================================
subroutine axial_field_lin(trial, lin, ssterm, ddterm, q_ni, q_bf)
  
  use basic
  use arrays
  use nintegrate_mod
  use metricterms_new

  implicit none

  vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: trial, lin 
  vectype, dimension(num_fields), intent(out) :: ssterm, ddterm
  vectype, intent(out) :: q_ni, q_bf
  vectype :: temp
  real :: thimpb

  vectype, dimension(MAX_PTS, OP_NUM) :: hi
  hi = hypi*sz79

  if(imp_mod.eq.0) then
     thimpb = thimp
  else
     thimpb = 1.
  endif

  ssterm = 0.
  ddterm = 0.
  q_ni = 0.
  q_bf = 0.

  if(numvar.lt.2) return          

  temp = b2psieta(trial,lin,eta79,hi)
  ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
  ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

  if(ibootstrap.gt.0) then
     temp = b2psimue(trial,lin,visc_e)*dbf
     ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp
  endif
         
  temp = b2psiv(trial,lin,vzt79)
  ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
  ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

  temp = b2psipsid(trial,lin,ps179,ni79)*dbf &
       + b2psipsid(trial,ps179,lin,ni79)*dbf &
       + b2psibd  (trial,lin,bz179,ni79)*dbf
  ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
  ddterm(psi_g) = ddterm(psi_g) + (.5-thimp*bdf)*dt*temp

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

  temp = b2psibd(trial,ps179,lin,ni79)*dbf &
       + b2bbd  (trial,lin,bz179,ni79)*dbf &
       + b2bbd  (trial,bz179,lin,ni79)*dbf
  ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
  ddterm(bz_g) = ddterm(bz_g) + (.5-thimp*bdf)*dt*temp

  temp = b2bu  (trial,bz179,lin)
  ssterm(u_g) = ssterm(u_g) - thimpb*dt*temp
  ddterm(u_g) = ddterm(u_g) - thimpb*dt*temp*bdf

  temp = b2psiv(trial,ps179,lin)
  ssterm(vz_g) = ssterm(vz_g) - thimpb*dt*temp
  ddterm(vz_g) = ddterm(vz_g) - thimpb*dt*temp*bdf

  temp = b2harnedmikic_psi(trial,lin)
  ssterm(psi_g) = ssterm(psi_g) + dt*dt*temp
  ddterm(psi_g) = ddterm(psi_g) + dt*dt*temp

  temp = b2harnedmikic_b(trial,lin)
  ssterm(bz_g) = ssterm(bz_g) + dt*dt*temp
  ddterm(bz_g) = ddterm(bz_g) + dt*dt*temp

  if(eqsubtract.eq.1) then
     temp = b2psipsid(trial,lin,ps079,ni79)*dbf &
          + b2psipsid(trial,ps079,lin,ni79)*dbf &
          + b2psibd  (trial,lin,bz079,ni79)*dbf
     ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

     temp = b2psibd(trial,ps079,lin,ni79)*dbf &
          + b2bbd  (trial,lin,bz079,ni79)*dbf &
          + b2bbd  (trial,bz079,lin,ni79)*dbf
     ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
     ddterm(bz_g) = ddterm(bz_g) + (1.-thimp*bdf)*dt*temp

     temp = b2bu  (trial,bz079,lin)
     ssterm(u_g) = ssterm(u_g) -     thimpb     *dt*temp
     ddterm(u_g) = ddterm(u_g) + (1.-thimpb*bdf)*dt*temp
     
     temp = b2psiv(trial,ps079,lin)
     ssterm(vz_g) = ssterm(vz_g) -     thimpb     *dt*temp
     ddterm(vz_g) = ddterm(vz_g) + (1.-thimpb*bdf)*dt*temp
  endif

  if(i3d.eq.1) then
     temp = b2fv(trial,bft79,lin)
     ssterm(vz_g) = ssterm(vz_g) - thimpb*dt*temp
     ddterm(vz_g) = ddterm(vz_g) - thimpb*dt*temp*bdf

     temp = b2psifd(trial,lin,bft79,ni79)*dbf
     ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

     temp = b2bfd(trial,lin,bft79,ni79)*dbf
     ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
     ddterm(bz_g) = ddterm(bz_g) + (1.-thimp*bdf)*dt*temp

     q_bf = q_bf + dt* &
          b2feta(trial,lin,eta79)

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
     ssterm(pe_g) = ssterm(pe_g) -     thimp     *dt*temp
     ddterm(pe_g) = ddterm(pe_g) + (1.-thimp*bdf)*dt*temp

     temp = b2bchi(trial,bz179,lin)                  
     ssterm(chi_g) = ssterm(chi_g) - thimpb*dt*temp
     ddterm(chi_g) = ddterm(chi_g) - thimpb*dt*temp*bdf
             
     if(eqsubtract.eq.1) then
        temp = b2bchi(trial,bz079,lin)
        ssterm(chi_g) = ssterm(chi_g) -     thimpb     *dt*temp
        ddterm(chi_g) = ddterm(chi_g) + (1.-thimpb*bdf)*dt*temp
     endif
  end if


  ! density terms
  ! ~~~~~~~~~~~~~
  if(idens.eq.1 .and. eqsubtract.eq.1) then
     q_ni = q_ni + dt* &
          (b2psipsid(trial,ps079,ps079,lin)*dbf &
          +b2psibd  (trial,ps079,bz079,lin)*dbf &
          +b2bbd    (trial,bz079,bz079,lin)*dbf)

     if(numvar.ge.3) then
        q_ni = q_ni + dt* &
             (b2ped (trial,pe079,lin)*dbf*pefac)
     endif
  endif

end subroutine axial_field_lin


subroutine axial_field_nolin(trial, r4term)

  use basic
  use nintegrate_mod
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
  use nintegrate_mod
  use metricterms_new

  implicit none

  vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: trial, lin 
  vectype, dimension(num_fields), intent(out) :: ssterm, ddterm
  vectype, intent(out) :: q_ni, q_bf
  vectype :: temp
  real :: thimpb

  vectype, dimension(MAX_PTS, OP_NUM) :: hp
  hp = hypp*sz79


  if(imp_mod.eq.0) then
     thimpb = thimp
  else
     thimpb = 1.
  endif

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
  endif

  temp = b3pepsid(trial,pe179,lin,ni79)*dbf*pefac
  ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
  ddterm(psi_g) = ddterm(psi_g) + (.5-thimp*bdf)*dt*temp
  
  temp = b3pebd(trial,pe179,lin,ni79)*dbf*pefac
  ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
  ddterm(bz_g) = ddterm(bz_g) + (.5-thimp*bdf)*dt*temp

  temp = b3pepsid(trial,lin,ps179,ni79)*dbf*pefac & 
       + b3pebd  (trial,lin,bz179,ni79)*dbf*pefac
  ssterm(pe_g) = ssterm(pe_g) -     thimp     *dt*temp
  ddterm(pe_g) = ddterm(pe_g) + (.5-thimp*bdf)*dt*temp

  temp = p1pu  (trial,lin,pht79) & 
       + p1pv  (trial,lin,vzt79) &
       + p1pchi(trial,lin,cht79) &
       + b3pedkappa(trial,lin,ni79,kap79,hp)
  ssterm(pe_g) = ssterm(pe_g) -     thimp     *dt*temp
  ddterm(pe_g) = ddterm(pe_g) + (1.-thimp*bdf)*dt*temp

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

  ! Anisotropic Heat Flux
  if(kappar.ne.0.) then
     ! Assumes no contribution from equilibrium f

     if(linear.eq.0) then
        temp = p1psipsikappar(trial,lin,ps179,pe179,ni79,b2i79,kar79) &
             + p1psipsikappar(trial,ps179,lin,pe179,ni79,b2i79,kar79) &
             + p1psibkappar  (trial,lin,bz179,pe179,ni79,b2i79,kar79)
        ssterm(psi_g) = ssterm(psi_g) -          thimp     *dt*temp
        ddterm(psi_g) = ddterm(psi_g) + (1./3. - thimp*bdf)*dt*temp

        temp = p1psibkappar  (trial,ps179,lin,pe179,ni79,b2i79,kar79)
        ssterm(bz_g) = ssterm(bz_g) -          thimp     *dt*temp
        ddterm(bz_g) = ddterm(bz_g) + (1./3. - thimp*bdf)*dt*temp
        
        temp = p1psipsikappar(trial,ps179,ps179,lin,ni79,b2i79,kar79) &
             + p1psibkappar  (trial,ps179,bz179,lin,ni79,b2i79,kar79)
        ssterm(pe_g) = ssterm(pe_g) -          thimp     *dt*temp
        ddterm(pe_g) = ddterm(pe_g) + (1./3. - thimp*bdf)*dt*temp
     endif

     if(eqsubtract.eq.1) then 
        temp = p1psipsikappar(trial,lin,ps079,pe079,ni79,b2i79,kar79) &
             + p1psipsikappar(trial,ps079,lin,pe079,ni79,b2i79,kar79) &
             + p1psibkappar  (trial,lin,bz079,pe079,ni79,b2i79,kar79)
        ssterm(psi_g) = ssterm(psi_g) -       thimp     *dt*temp
        ddterm(psi_g) = ddterm(psi_g) + (1. - thimp*bdf)*dt*temp

        temp = p1psibkappar  (trial,ps079,lin,pe079,ni79,b2i79,kar79)
        ssterm(bz_g) = ssterm(bz_g) -       thimp     *dt*temp
        ddterm(bz_g) = ddterm(bz_g) + (1. - thimp*bdf)*dt*temp

        temp = p1psipsikappar(trial,ps079,ps079,lin,ni79,b2i79,kar79) &
             + p1psibkappar  (trial,ps079,bz079,lin,ni79,b2i79,kar79)
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
                + p1psibkappar  (trial,ps079,lin,pe179,ni79,b2i79,kar79)
           ssterm(bz_g) = ssterm(bz_g) -          thimp     *dt*temp
           ddterm(bz_g) = ddterm(bz_g) + (1./2. - thimp*bdf)*dt*temp
           
           temp = p1psipsikappar(trial,ps179,ps079,lin,ni79,b2i79,kar79) &
                + p1psipsikappar(trial,ps079,ps179,lin,ni79,b2i79,kar79) &
                + p1psibkappar  (trial,ps179,bz079,lin,ni79,b2i79,kar79) &
                + p1psibkappar  (trial,ps079,bz179,lin,ni79,b2i79,kar79)
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
     ssterm(psi_g) = ssterm(psi_g) -     thimp     *dt*temp
     ddterm(psi_g) = ddterm(psi_g) + (1.-thimp*bdf)*dt*temp

     temp = b3pebd(trial,pe079,lin,ni79)*dbf*pefac
     ssterm(bz_g) = ssterm(bz_g) -     thimp     *dt*temp
     ddterm(bz_g) = ddterm(bz_g) + (1.-thimp*bdf)*dt*temp

     temp = b3pepsid(trial,lin,ps079,ni79)*dbf*pefac &
          + b3pebd  (trial,lin,bz079,ni79)*dbf*pefac
     ssterm(pe_g) = ssterm(pe_g) -     thimp     *dt*temp
     ddterm(pe_g) = ddterm(pe_g) + (1.-thimp*bdf)*dt*temp
     
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
     ssterm(pe_g) = ssterm(pe_g) -     thimp     *dt*temp
     ddterm(pe_g) = ddterm(pe_g) + (1.-thimp*bdf)*dt*temp

     if(eqsubtract.eq.1) then
        q_bf = q_bf + dt* &
             (b3pefd(trial,pe079,lin,ni79)*dbf*pefac)
     endif
  endif


  ! density terms
  ! ~~~~~~~~~~~~~
  if(idens.eq.1 .and. eqsubtract.eq.1) then
     q_ni = q_ni + dt* &
          (b3pepsid(trial,pe079,ps079,lin)*dbf*pefac &
          +b3pebd  (trial,pe079,bz079,lin)*dbf*pefac &
          +p1psipsikappar(trial,ps079,ps079,pe079,lin,b2i79,kar79) &
          +p1psibkappar  (trial,ps079,bz079,pe079,lin,b2i79,kar79))
  endif

end subroutine electron_pressure_lin

subroutine electron_pressure_nolin(trial, r4term)

  use basic
  use nintegrate_mod
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
  use t_data
  use basic
  use arrays
  use sparse
  use nintegrate_mod
  use diagnostics

  implicit none

! include 'mpif.h'
#include "finclude/petsc.h"

  integer, intent(in) :: ivel_def   ! populate velocity advance matrices
  integer, intent(in) :: idens_def  ! populate density advance matrices
  integer, intent(in) :: ipres_def  ! populate pressure advance matrices
  integer, intent(in) :: ifield_def ! populate field advance matrices 

  integer :: itri, numelms
  integer :: def_fields
  real :: x, z, xmin, zmin, factor

  real :: tstart, tend, tfield, telm, tsizefield, tfinalize
  logical :: is_edge(3)  ! is inode on boundary
  real :: n(2,3)
  integer :: iedge, idim(3)


  double precision cogcoords(3)

  integer :: ier
  PetscTruth :: flg_petsc, flg_solve2, flg_solve1 

  tfield = 0.
  telm = 0.
  tsizefield = 0.
  tfinalize = 0.

  call numfac(numelms)
  call getmincoord(xmin,zmin)

  ! initialize matrices
  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, " initializing matrices..."

  ! default linear solver superlu cj-april-09-2008 
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-ipetsc', flg_petsc,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve2', flg_solve2,ier)
  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-solve1', flg_solve1,ier)

  ! Clear matrices
  select case(isplitstep)
  case(0)
     if(flg_petsc.eq.PETSC_TRUE) then
        call zeropetscmatrix(s1matrix_sm,icomplex, vecsize_vel)
        if(myrank.eq.0 .and. iprint.ge.2) &
        print *, "	ludef_t_ludefall zeropetscmatrix", s1matrix_sm 
     else
        call zerosuperlumatrix(s1matrix_sm,icomplex, vecsize_vel)
        if(myrank.eq.0 .and. iprint.ge.2) &
        print *, "	ludef_t_ludefall zerosuperlumatrix", s1matrix_sm
     endif
     call zeromultiplymatrix(d1matrix_sm,icomplex,vecsize_vel)
     if(idens_def.eq.1 .and. eqsubtract.eq.1) &
          call zeromultiplymatrix(q42matrix_sm,icomplex,vecsize_vel)
#ifdef USECOMPLEX
     call zeromultiplymatrix(o1matrix_sm,icomplex,vecsize_vel)
#endif
     q4 = 0.

  case(1)

     if(ivel_def.eq.1) then
        if(flg_petsc.eq.PETSC_TRUE) then
           call zeropetscmatrix(s1matrix_sm,icomplex,vecsize_vel)
           if(myrank.eq.0 .and. iprint.ge.2) &
           print *, "	ludef_t_ludefall zeropetscmatrix", s1matrix_sm 
        else
           call zerosuperlumatrix(s1matrix_sm,icomplex,vecsize_vel)
           if(myrank.eq.0 .and. iprint.ge.2) &
           print *, "	ludef_t_ludefall zerosuperlumatrix", s1matrix_sm
        endif
        call zeromultiplymatrix(d1matrix_sm,icomplex,vecsize_vel)
        call zeromultiplymatrix(q1matrix_sm,icomplex,vecsize_phi)
        call zeromultiplymatrix(r14matrix_sm,icomplex,vecsize_phi)
#ifdef USECOMPLEX
        call zeromultiplymatrix(o1matrix_sm,icomplex,vecsize_phi)
#endif   
        r4 = 0.
     end if

     if(ifield_def.eq.1) then
        if(flg_petsc.eq.PETSC_TRUE) then
           call zeropetscmatrix(s2matrix_sm,icomplex,vecsize_phi)
           if(myrank.eq.0 .and. iprint.ge.2) &
           print *, "	ludef_t_ludefall zeropetscmatrix", s2matrix_sm 
        else
           call zerosuperlumatrix(s2matrix_sm,icomplex,vecsize_phi)
           if(myrank.eq.0 .and. iprint.ge.2) &
           print *, "	ludef_t_ludefall zerosuperlumatrix", s2matrix_sm
        endif 
        call zeromultiplymatrix(d2matrix_sm,icomplex,vecsize_phi)
        call zeromultiplymatrix(r2matrix_sm,icomplex,vecsize_vel)
        call zeromultiplymatrix(q2matrix_sm,icomplex,vecsize_vel)
        if(idens_def.eq.1 .and. eqsubtract.eq.1) &
             call zeromultiplymatrix(q42matrix_sm,icomplex,vecsize_vel)
#ifdef USECOMPLEX
        call zeromultiplymatrix(o2matrix_sm,icomplex,vecsize_phi)
#endif
        q4 = 0.
     end if

     if(idens_def.eq.1) then
        if(flg_petsc.eq.PETSC_TRUE) then
           call zeropetscmatrix(s8matrix_sm,icomplex,vecsize_n)
           if(myrank.eq.0 .and. iprint.ge.2) &
           print *, "	ludef_t_ludefall zeropetscmatrix", s8matrix_sm 
        else 
           call zerosuperlumatrix(s8matrix_sm,icomplex,vecsize_n)
           if(myrank.eq.0 .and. iprint.ge.2) &
           print *, "	ludef_t_ludefall zerosuperlumatrix", s8matrix_sm
        endif 
        call zeromultiplymatrix(d8matrix_sm,icomplex,vecsize_n)
        call zeromultiplymatrix(q8matrix_sm,icomplex,vecsize_vel)
        call zeromultiplymatrix(r8matrix_sm,icomplex,vecsize_vel)
        qn4 = 0.
     endif
     if(ipres_def.eq.1) then
        if(flg_petsc.eq.PETSC_TRUE) then
           call zeropetscmatrix(s9matrix_sm,icomplex,vecsize_p)
           if(myrank.eq.0 .and. iprint.ge.2) &
           print *, "	ludef_t_ludefall zeropetscmatrix", s9matrix_sm 
        else
           call zerosuperlumatrix(s9matrix_sm,icomplex,vecsize_p)
           if(myrank.eq.0 .and. iprint.ge.2) &
           print *, "	ludef_t_ludefall zerosuperlumatrix", s9matrix_sm
        endif
        call zeromultiplymatrix(d9matrix_sm,icomplex,vecsize_p)
        call zeromultiplymatrix(q9matrix_sm,icomplex,vecsize_vel)
        call zeromultiplymatrix(r9matrix_sm,icomplex,vecsize_vel)
        qp4 = 0.
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

  if(isources.eq.1) def_fields = def_fields + FIELD_SRC
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

     call cogfac(itri,cogcoords)
     x = cogcoords(1)-xmin
     z = cogcoords(2)-zmin
   
     if(imask.eq.1) then
        call mask(x,z,factor)
        dbf = db*factor
     else
        dbf = db
     endif

     ! calculate the field values and derivatives at the sampling points
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call define_triangle_quadrature(itri, int_pts_main)
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

        call define_edge_quadrature(itri, iedge, 5, n, idim)
        call define_fields(itri, def_fields, 1, linear)

        if(ivel_def.eq.1) call ludefvel_n(itri)
        if(ifield_def.eq.1) call ludefphi_n(itri)
!!$        if(idens_def.eq.1) call ludefden_n(itri)
!!$        if(ipres_def.eq.1) call ludefpres_n(itri)
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
     call finalizematrix(d1matrix_sm)
     if(idens_def.eq.1 .and. eqsubtract.eq.1) &
          call finalizematrix(q42matrix_sm)
#ifdef USECOMPLEX
     call finalizematrix(o1matrix_sm)
#endif   
     call sumsharedppplvecvals(q4)

  case(1)
     if(ivel_def.eq.1) then
        call finalizematrix(d1matrix_sm)
        call finalizematrix(q1matrix_sm)
        call finalizematrix(r14matrix_sm)
#ifdef USECOMPLEX
        call finalizematrix(o1matrix_sm)
#endif   
        call sumsharedppplvecvals(r4)
     end if
     
     if(ifield_def.eq.1) then
        call finalizematrix(d2matrix_sm)
        call finalizematrix(r2matrix_sm)
        call finalizematrix(q2matrix_sm)
#ifdef USECOMPLEX
        call finalizematrix(o2matrix_sm)
#endif   
        if(idens_def.eq.1 .and. eqsubtract.eq.1) &
             call finalizematrix(q42matrix_sm)
        call sumsharedppplvecvals(q4)
     end if
     
     if(idens_def.eq.1) then
        call finalizematrix(d8matrix_sm)
        call finalizematrix(q8matrix_sm)
        call finalizematrix(r8matrix_sm)
        call sumsharedppplvecvals(qn4)
     endif ! on idens_def.eq.1
     
     
     if(ipres_def.eq.1) then
        call finalizematrix(d9matrix_sm)
        call finalizematrix(q9matrix_sm)
        call finalizematrix(r9matrix_sm)
        call sumsharedppplvecvals(qp4)
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
  use t_data
  use nintegrate_mod
  use arrays
  use sparse

  implicit none

  integer, intent(in) :: itri

  integer :: i, ii, iii, j, jj, jjj, ip, iv, jp, jv
  integer :: ib_vel, ib_phi, jb_vel, jb_phi, iendplusone

  vectype :: temp

  vectype, dimension(num_fields,num_fields) :: ss, dd
  vectype, dimension(num_fields) :: q_bf

  integer :: vv1, vv0, vb1, vb0, vn1, vn0, vf0
  integer :: advfield
  integer :: pp_off
  vectype, pointer :: vsource(:)

  if(isplitstep.eq.1) then
     vv1 = s1matrix_sm
     vv0 = d1matrix_sm
     vb0 = q1matrix_sm
     vn0 = r14matrix_sm
     vf0 = o1matrix_sm
     vsource => r4
     pp_off = pe_off
  else
     vv1 = s1matrix_sm
     vv0 = d1matrix_sm
     vb1 = s1matrix_sm
     vb0 = d1matrix_sm
     vn1 = s1matrix_sm
     vn0 = d1matrix_sm
     vf0 = o1matrix_sm
     vsource => q4
     pp_off = p_off
  endif

  if(isplitstep.eq.1 .and. iestatic.eq.0) then
     advfield = 1 
  else 
     advfield = 0
  endif

  do iii=1,3
  call entdofs(vecsize_phi,  ist(itri,iii)+1, 0, ib_phi, iendplusone)
  call entdofs(vecsize_vel,  ist(itri,iii)+1, 0, ib_vel, iendplusone)

  do ii=1,6
     i = (iii-1)*6 + ii
     iv = ib_vel + ii - 1
     ip = ib_phi + ii - 1

     do jjj=1,3
     call entdofs(vecsize_phi,  ist(itri,jjj)+1, 0, jb_phi, iendplusone)
     call entdofs(vecsize_vel,  ist(itri,jjj)+1, 0, jb_vel, iendplusone)
     do jj=1,6
        j = (jjj-1)*6 + jj
        jv = jb_vel + jj - 1
        jp = jb_phi + jj - 1

        if(istatic.eq.1) then
           ss = 0.
           dd = 0.
           temp = int2(g79(:,OP_1,i),g79(:,OP_1,j))
           ss(u_g,u_g) = temp; ss(vz_g, vz_g) = temp; ss(chi_g,chi_g) = temp
           dd(u_g,u_g) = temp; dd(vz_g, vz_g) = temp; dd(chi_g,chi_g) = temp
        else 
           call vorticity_lin(g79(:,:,i),g79(:,:,j), &
                ss(u_g,:),dd(u_g,:),q_bf(u_g),advfield)

           if(surface_int) then
              ss(vz_g,:) = 0.
              dd(vz_g,:) = 0.
              q_bf(vz_g) = 0.
           else
              if(numvar.ge.2) then
                 call axial_vel_lin(g79(:,:,i),g79(:,:,j), &
                      ss(vz_g,:),dd(vz_g,:),q_bf(vz_g),advfield)
              endif
           end if
           if(numvar.ge.3) then
              call compression_lin(g79(:,:,i),g79(:,:,j), &
                   ss(chi_g,:),dd(chi_g,:),q_bf(chi_g),advfield)
           endif
        endif

        if(iestatic.eq.1 .and. isplitstep.eq.1) then 
           ss(:,psi_g) = 0; ss(:,bz_g) = 0; ss(:,p_g) = 0
           dd(:,psi_g) = 0; dd(:,bz_g) = 0; dd(:,p_g) = 0
        end if

        call insval(vv1,ss(  u_g,  u_g),icomplex,iv+  u_off,jv+  u_off,1)
        call insval(vv0,dd(  u_g,  u_g),icomplex,iv+  u_off,jv+  u_off,1)
        call insval(vb0,dd(  u_g,psi_g),icomplex,ip+  u_off,jp+psi_off,1)
        if(idens.eq.1) &
             call insval(vn0,dd(  u_g,den_g),icomplex,iv+  u_off,jv+den_off,1)
        if(isplitstep.eq.0) then
           call insval(vb1,ss(  u_g, psi_g),icomplex,ip+  u_off,jp+psi_off,1)
           if(idens.eq.1) &
                call insval(vn1,ss(  u_g,den_g),icomplex,iv+  u_off,jv+den_off,1)
        endif
        if(i3d.eq.1) then
           call insval(vf0,q_bf(u_g),icomplex,iv+u_off,jv+bf_off,1)
        endif
        if(numvar.ge.2) then
           call insval(vv1,ss(  u_g, vz_g),icomplex,iv+  u_off,jv+ vz_off,1)
           call insval(vv1,ss( vz_g,  u_g),icomplex,iv+ vz_off,jv+  u_off,1)
           call insval(vv1,ss( vz_g, vz_g),icomplex,iv+ vz_off,jv+ vz_off,1)
           call insval(vv0,dd(  u_g, vz_g),icomplex,iv+  u_off,jv+ vz_off,1)
           call insval(vv0,dd( vz_g,  u_g),icomplex,iv+ vz_off,jv+  u_off,1)
           call insval(vv0,dd( vz_g, vz_g),icomplex,iv+ vz_off,jv+ vz_off,1)
           call insval(vb0,dd(  u_g, bz_g),icomplex,ip+  u_off,jp+ bz_off,1)
           call insval(vb0,dd( vz_g,psi_g),icomplex,ip+ vz_off,jp+psi_off,1)
           call insval(vb0,dd( vz_g, bz_g),icomplex,ip+ vz_off,jp+ bz_off,1)
           if(idens.eq.1) &
                call insval(vn0,dd( vz_g,den_g),icomplex,iv+vz_off,jv+den_off,1)
           if(isplitstep.eq.0) then
              call insval(vb1,ss(  u_g, bz_g),icomplex,ip+  u_off,jp+ bz_off,1)
              call insval(vb1,ss( vz_g,psi_g),icomplex,ip+ vz_off,jp+psi_off,1)
              call insval(vb1,ss( vz_g, bz_g),icomplex,ip+ vz_off,jp+ bz_off,1)
              if(idens.eq.1) &
                   call insval(vn1,ss( vz_g,den_g),icomplex,iv+vz_off,jv+den_off,1)
           end if
           if(i3d.eq.1) then
              call insval(vf0,q_bf(vz_g),icomplex,iv+vz_off,jv+bf_off,1)
           endif
        endif
        if(numvar.ge.3) then
           call insval(vv1,ss(  u_g,chi_g),icomplex,iv+  u_off,jv+chi_off,1)
           call insval(vv1,ss( vz_g,chi_g),icomplex,iv+ vz_off,jv+chi_off,1)
           call insval(vv1,ss(chi_g,  u_g),icomplex,iv+chi_off,jv+  u_off,1)
           call insval(vv1,ss(chi_g, vz_g),icomplex,iv+chi_off,jv+ vz_off,1)
           call insval(vv1,ss(chi_g,chi_g),icomplex,iv+chi_off,jv+chi_off,1)
           call insval(vv0,dd(  u_g,chi_g),icomplex,iv+  u_off,jv+chi_off,1)
           call insval(vv0,dd( vz_g,chi_g),icomplex,iv+ vz_off,jv+chi_off,1)
           call insval(vv0,dd(chi_g,  u_g),icomplex,iv+chi_off,jv+  u_off,1)
           call insval(vv0,dd(chi_g, vz_g),icomplex,iv+chi_off,jv+ vz_off,1)
           call insval(vv0,dd(chi_g,chi_g),icomplex,iv+chi_off,jv+chi_off,1)
           call insval(vb0,dd(  u_g,  p_g),icomplex,ip+  u_off,jp+ pp_off,1)
           call insval(vb0,dd( vz_g,  p_g),icomplex,ip+ vz_off,jp+ pp_off,1)
           call insval(vb0,dd(chi_g,psi_g),icomplex,ip+chi_off,jp+psi_off,1)
           call insval(vb0,dd(chi_g, bz_g),icomplex,ip+chi_off,jp+ bz_off,1)
           call insval(vb0,dd(chi_g,  p_g),icomplex,ip+chi_off,jp+ pp_off,1)
           if(idens.eq.1) &
                call insval(vn0,dd(chi_g,den_g),icomplex,iv+chi_off,jv+den_off,1)
           if(isplitstep.eq.0) then
              call insval(vb1,ss(  u_g,  p_g),icomplex,ip+  u_off,jp+ pp_off,1)
              call insval(vb1,ss( vz_g,  p_g),icomplex,ip+ vz_off,jp+ pp_off,1)
              call insval(vb1,ss(chi_g,psi_g),icomplex,ip+chi_off,jp+psi_off,1)
              call insval(vb1,ss(chi_g, bz_g),icomplex,ip+chi_off,jp+ bz_off,1)
              call insval(vb1,ss(chi_g,  p_g),icomplex,ip+chi_off,jp+ pp_off,1)
              if(idens.eq.1) &
                   call insval(vn1,ss(chi_g,den_g),icomplex,iv+chi_off,jv+den_off,1)
           endif
           if(i3d.eq.1) then
              call insval(vf0,q_bf(chi_g),icomplex,iv+chi_off,jv+bf_off,1)
           endif
        endif
     enddo               ! on j
     enddo

     ! Definition of R4
     ! ================
     if(istatic.eq.1) then
        vsource(iv+u_off) = 0.
        if(numvar.ge.2) vsource(iv+vz_off) = 0.
        if(numvar.ge.3) vsource(iv+chi_off) = 0.
     else if(linear.eq.0) then  
        call vorticity_nolin(g79(:,:,i),temp)
        vsource(iv+  u_off) = vsource(iv+  u_off) + temp
        if(numvar.ge.2) then
           call axial_vel_nolin(g79(:,:,i),temp)
           vsource(iv+ vz_off) = vsource(iv+ vz_off) + temp
        endif
        if(numvar.ge.3) then
           call compression_nolin(g79(:,:,i),temp)
           vsource(iv+chi_off) = vsource(iv+chi_off) + temp
        endif
     end if
  enddo                  ! on i
  enddo

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
  use t_data
  use nintegrate_mod
  use arrays
  use sparse

  implicit none

  integer, intent(in) :: itri

  integer :: i, ii, iii, j, jj, jjj, ip, iv, jp, jv
  integer :: ib_vel, ib_phi, jb_vel, jb_phi, iendplusone
  
  vectype, dimension(num_fields,num_fields) :: ss, dd
  vectype, dimension(num_fields) :: q_ni, q_bf
  vectype :: temp

  integer :: bb1, bb0, bv1, bv0, bbf, bni
  vectype, pointer :: bsource(:)

  if(isplitstep.eq.1) then
     bb1 = s2matrix_sm
     bb0 = d2matrix_sm
     bv1 = r2matrix_sm
     bv0 = q2matrix_sm
     bbf = o2matrix_sm
     bni = q42matrix_sm
     bsource => q4
  else
     bb1 = s1matrix_sm
     bb0 = d1matrix_sm
     bv1 = s1matrix_sm
     bv0 = d1matrix_sm
     bbf = o1matrix_sm
     bni = q42matrix_sm
     bsource => q4
  endif

  do iii=1,3
  call entdofs(vecsize_phi,  ist(itri,iii)+1, 0, ib_phi, iendplusone)
  call entdofs(vecsize_vel,  ist(itri,iii)+1, 0, ib_vel, iendplusone)

  do ii=1,6
     i = (iii-1)*6 + ii
     iv = ib_vel + ii - 1
     ip = ib_phi + ii - 1

     do jjj=1,3
     call entdofs(vecsize_phi,  ist(itri,jjj)+1, 0, jb_phi, iendplusone)
     call entdofs(vecsize_vel,  ist(itri,jjj)+1, 0, jb_vel, iendplusone)
     do jj=1,6
        j = (jjj-1)*6 + jj
        jv = jb_vel + jj - 1
        jp = jb_phi + jj - 1

        if(surface_int) then
           ss(psi_g,:) = 0.
           dd(psi_g,:) = 0.
           q_ni(psi_g) = 0.
           q_bf(psi_g) = 0.
        else
           call flux_lin(g79(:,:,i),g79(:,:,j), &
                ss(psi_g,:),dd(psi_g,:),q_ni(psi_g),q_bf(psi_g))
        endif
        if(numvar.ge.2) then
           call axial_field_lin(g79(:,:,i),g79(:,:,j), &
                ss(bz_g,:),dd(bz_g,:),q_ni(bz_g),q_bf(bz_g))
        endif
        if(numvar.ge.3) then
           call electron_pressure_lin(g79(:,:,i),g79(:,:,j), &
                ss(pe_g,:),dd(pe_g,:),q_ni(pe_g),q_bf(pe_g))
        endif

        call insval(bb1,ss(psi_g,psi_g),icomplex,ip+psi_off,jp+psi_off,1)
        call insval(bb0,dd(psi_g,psi_g),icomplex,ip+psi_off,jp+psi_off,1)
        call insval(bv1,ss(psi_g,  u_g),icomplex,iv+psi_off,jv+  u_off,1)
        call insval(bv0,dd(psi_g,  u_g),icomplex,iv+psi_off,jv+  u_off,1)
        if(i3d.eq.1) then
           call insval(bbf,q_bf(psi_g),icomplex,ip+psi_off,jp+bf_off,1)
        endif
        if(numvar.ge.2) then
           call insval(bb1,ss(psi_g, bz_g),icomplex,ip+psi_off,jp+ bz_off,1)
           call insval(bb1,ss( bz_g,psi_g),icomplex,ip+ bz_off,jp+psi_off,1)
           call insval(bb1,ss( bz_g, bz_g),icomplex,ip+ bz_off,jp+ bz_off,1)
           call insval(bb0,dd(psi_g, bz_g),icomplex,ip+psi_off,jp+ bz_off,1)
           call insval(bb0,dd( bz_g,psi_g),icomplex,ip+ bz_off,jp+psi_off,1)
           call insval(bb0,dd( bz_g, bz_g),icomplex,ip+ bz_off,jp+ bz_off,1)
           call insval(bv1,ss(psi_g, vz_g),icomplex,iv+psi_off,jv+ vz_off,1)
           call insval(bv1,ss( bz_g,  u_g),icomplex,iv+ bz_off,jv+  u_off,1)
           call insval(bv1,ss( bz_g, vz_g),icomplex,iv+ bz_off,jv+ vz_off,1)
           call insval(bv0,dd(psi_g, vz_g),icomplex,iv+psi_off,jv+ vz_off,1)
           call insval(bv0,dd( bz_g,  u_g),icomplex,iv+ bz_off,jv+  u_off,1)
           call insval(bv0,dd( bz_g, vz_g),icomplex,iv+ bz_off,jv+ vz_off,1)
           if(i3d.eq.1) then
              call insval(bbf,q_bf(bz_g),icomplex,ip+bz_off,jp+bf_off,1)
           endif
        endif
        if(numvar .eq. 3) then        
           call insval(bb1,ss(psi_g, pe_g),icomplex,ip+psi_off,jp+ pe_off,1)
           call insval(bb1,ss( bz_g, pe_g),icomplex,ip+ bz_off,jp+ pe_off,1)
           call insval(bb1,ss( pe_g, pe_g),icomplex,ip+ pe_off,jp+ pe_off,1)
           call insval(bb1,ss( pe_g,psi_g),icomplex,ip+ pe_off,jp+psi_off,1)
           call insval(bb1,ss( pe_g, bz_g),icomplex,ip+ pe_off,jp+ bz_off,1)
           call insval(bb0,dd(psi_g, pe_g),icomplex,ip+psi_off,jp+ pe_off,1)
           call insval(bb0,dd( bz_g, pe_g),icomplex,ip+ bz_off,jp+ pe_off,1)
           call insval(bb0,dd( pe_g, pe_g),icomplex,ip+ pe_off,jp+ pe_off,1)
           call insval(bb0,dd( pe_g,psi_g),icomplex,ip+ pe_off,jp+psi_off,1)
           call insval(bb0,dd( pe_g, bz_g),icomplex,ip+ pe_off,jp+ bz_off,1)
           call insval(bv1,ss(psi_g,chi_g),icomplex,iv+psi_off,jv+chi_off,1)
           call insval(bv1,ss( bz_g,chi_g),icomplex,iv+ bz_off,jv+chi_off,1)
           call insval(bv1,ss( pe_g,chi_g),icomplex,iv+ pe_off,jv+chi_off,1)
           call insval(bv1,ss( pe_g,  u_g),icomplex,iv+ pe_off,jv+  u_off,1)
           call insval(bv1,ss( pe_g, vz_g),icomplex,iv+ pe_off,jv+ vz_off,1)
           call insval(bv0,dd(psi_g,chi_g),icomplex,iv+psi_off,jv+chi_off,1)
           call insval(bv0,dd( bz_g,chi_g),icomplex,iv+ bz_off,jv+chi_off,1)
           call insval(bv0,dd( pe_g,chi_g),icomplex,iv+ pe_off,jv+chi_off,1)
           call insval(bv0,dd( pe_g,  u_g),icomplex,iv+ pe_off,jv+  u_off,1)
           call insval(bv0,dd( pe_g, vz_g),icomplex,iv+ pe_off,jv+ vz_off,1)
           if(i3d.eq.1) then
              call insval(bbf,q_bf(pe_g),icomplex,ip+pe_off,jp+bf_off,1)
           endif
           if(idens.eq.1 .and. eqsubtract.eq.1) then
              call insval(bni,q_ni(psi_g),icomplex,ip+psi_off,jp,1)
              if(numvar.ge.2) &
                   call insval(bni,q_ni(bz_g),icomplex,ip+bz_off,jp,1)
              if(numvar.ge.3) &
                   call insval(bni,q_ni(pe_g),icomplex,ip+pe_off,jp,1)
           endif
        endif
       
     enddo ! on j
     enddo

     call flux_nolin(g79(:,:,i),temp)
     bsource(ip+psi_off) = bsource(ip+psi_off) + temp
     if(numvar.ge.2) then
        call axial_field_nolin(g79(:,:,i),temp)
        bsource(ip+ bz_off) = bsource(ip+ bz_off) + temp
     endif
     if(numvar.ge.3) then
        call electron_pressure_nolin(g79(:,:,i),temp)
        bsource(ip+ pe_off) = bsource(ip+ pe_off) + temp
     endif

  enddo ! on i
  enddo
     

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
  use nintegrate_mod
  use arrays
  use sparse
  use metricterms_new

  implicit none

  integer, intent(in) :: itri

  integer :: i, i1, ione, j, j1, jone
  vectype :: ssterm, ddterm
  vectype, dimension(3) :: rrterm, qqterm

  vectype :: temp

  integer :: nn1, nn0, nv1, nv0
  vectype, pointer :: nsource(:)
  real :: thimpb

  vectype, dimension(MAX_PTS,OP_NUM) :: hp
  hp = hypp*sz79

  if(imp_mod.eq.0) then
     thimpb = thimp
  else
     thimpb = 1.
  endif

  if(isplitstep.eq.1) then
     nn1 = s8matrix_sm
     nn0 = d8matrix_sm
     nv1 = r8matrix_sm
     nv0 = q8matrix_sm
     nsource => qn4
  else
     nn1 = s1matrix_sm
     nn0 = d1matrix_sm
     nv1 = s1matrix_sm
     nv0 = d1matrix_sm
     nsource => q4
  endif

  do i=1,18
     if(isplitstep.eq.1) then
        ione = isval1(itri,i)
     else
        ione = isvaln(itri,i)
     endif
     i1 = isvaln(itri,i)

     do j=1,18         
        ssterm = 0.
        ddterm = 0.
        rrterm = 0.
        qqterm = 0.

        if(isplitstep.eq.1) then
           jone = isval1(itri,j)
        else
           jone = isvaln(itri,j)
        endif
        j1 = isvaln(itri,j)

        ! NUMVAR = 1
        ! ~~~~~~~~~~
        temp = n1n(g79(:,:,i),g79(:,:,j))
        ssterm = ssterm + temp    
        ddterm = ddterm + temp*bdf
        
        temp = n1ndenm(g79(:,:,i),g79(:,:,j),denm,hp) &
             + n1nu   (g79(:,:,i),g79(:,:,j),pht79)
        ssterm = ssterm -     thimp     *dt*temp
        ddterm = ddterm + (1.-thimp*bdf)*dt*temp

        temp = n1nu(g79(:,:,i),n179,g79(:,:,j))
        rrterm(1) = rrterm(1) + thimpb*dt*temp
        qqterm(1) = qqterm(1) - thimpb*dt*temp*bdf

        if(eqsubtract.eq.1) then
           temp = n1nu  (g79(:,:,i),n079,g79(:,:,j))
           rrterm(1) = rrterm(1) +     thimpb     *dt*temp
           qqterm(1) = qqterm(1) + (1.-thimpb*bdf)*dt*temp
        endif

#ifdef USECOMPLEX
        ! NUMVAR = 2
        ! ~~~~~~~~~~
        if(numvar.ge.2) then
           temp = n1nv(g79(:,:,i),g79(:,:,j),vzt79)
           ssterm = ssterm -     thimp     *dt*temp
           ddterm = ddterm + (1.-thimp*bdf)*dt*temp
           
           temp = n1nv(g79(:,:,i),n179,g79(:,:,j))
           rrterm(2) = rrterm(2) + thimpb*dt*temp
           qqterm(2) = qqterm(2) - thimpb*dt*temp*bdf

           if(eqsubtract.eq.1) then
              temp = n1nv(g79(:,:,i),n079,g79(:,:,j))
              rrterm(2) = rrterm(2) +     thimpb     *dt*temp
              qqterm(2) = qqterm(2) + (1.-thimpb*bdf)*dt*temp
           endif
        endif
#endif

        ! NUMVAR = 3
        ! ~~~~~~~~~~
        if(numvar.ge.3) then
           temp = n1nchi(g79(:,:,i),g79(:,:,j),cht79)
           ssterm = ssterm -     thimp     *dt*temp
           ddterm = ddterm + (1.-thimp*bdf)*dt*temp
           
           temp = n1nchi(g79(:,:,i),n179,g79(:,:,j))
           rrterm(3) = rrterm(3) + thimpb*dt*temp
           qqterm(3) = qqterm(3) - thimpb*dt*temp*bdf

           if(eqsubtract.eq.1) then
              temp = n1nchi(g79(:,:,i),n079,g79(:,:,j))
              rrterm(3) = rrterm(3) +     thimpb     *dt*temp
              qqterm(3) = qqterm(3) + (1.-thimpb*bdf)*dt*temp
           endif
        endif

        if(isplitstep.eq.0) rrterm = -rrterm

        call insval(nn1, ssterm, icomplex, ione+den_off, jone+den_off, 1)
        call insval(nn0, ddterm, icomplex, ione+den_off, jone+den_off, 1)
        call insval(nv1, rrterm(1), icomplex, i1+den_off,j1+  u_off, 1)
        call insval(nv0, qqterm(1), icomplex, i1+den_off,j1+  u_off, 1)
        if(numvar.ge.2) then
           call insval(nv1,rrterm(2), icomplex, i1+den_off,j1+vz_off,1)
           call insval(nv0,qqterm(2), icomplex, i1+den_off,j1+vz_off,1)
        endif
        if(numvar.ge.3) then
           call insval(nv1,rrterm(3), icomplex, i1+den_off,j1+chi_off,1)
           call insval(nv0,qqterm(3), icomplex, i1+den_off,j1+chi_off,1)
        endif
     enddo                     ! on j

     nsource(ione+den_off) = nsource(ione+den_off) + dt* &
             n1s(g79(:,:,i),sig79)

     if(eqsubtract.eq.1) then
        nsource(ione+den_off) = nsource(ione+den_off) + dt* &
             (n1ndenm(g79(:,:,i),n079,denm,hp))
     endif
     
  enddo                     ! on i

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
  use nintegrate_mod
  use arrays
  use sparse
  use metricterms_new

  implicit none

  integer, intent(in) :: itri

  integer :: i, i1, ione, j, j1, jone
  vectype :: ssterm, ddterm
  vectype, dimension(3) :: rrterm, qqterm

  integer :: pp1, pp0, pv1, pv0
  vectype, pointer :: psource(:)

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
     pp1 = s9matrix_sm
     pp0 = d9matrix_sm
     pv1 = r9matrix_sm
     pv0 = q9matrix_sm
     psource => qp4
  else
     pp1 = s1matrix_sm
     pp0 = d1matrix_sm
     pv1 = s1matrix_sm
     pv0 = d1matrix_sm
     psource => q4
  endif

  do i=1,18
     if(isplitstep.eq.1) then
        ione = isval1(itri,i)
     else
        ione = isvaln(itri,i)
     endif
     i1 = isvaln(itri,i)
     
     do j=1,18         
        ssterm = 0.
        ddterm = 0.
        rrterm = 0.
        qqterm = 0.

        if(isplitstep.eq.1) then
           jone = isval1(itri,j)
        else
           jone = isvaln(itri,j)
        endif
        j1 = isvaln(itri,j)


        ! NUMVAR = 1
        ! ~~~~~~~~~~
        temp = int2(g79(:,:,i),g79(:,:,j))
        ssterm = ssterm + temp
        ddterm = ddterm + temp*bdf
        
        temp = p1pu   (g79(:,:,i),g79(:,:,j),pht79)
        ssterm = ssterm -     thimp     *dt*temp
        ddterm = ddterm + (1.-thimp*bdf)*dt*temp

        temp = p1pu(g79(:,:,i),p179,g79(:,:,j))
        rrterm(1) = rrterm(1) + thimpb*dt*temp
        qqterm(1) = qqterm(1) - thimpb*dt*temp*bdf

        if(eqsubtract.eq.1) then
           temp = p1pu  (g79(:,:,i),p079,g79(:,:,j))
           rrterm(1) = rrterm(1) +     thimpb     *dt*temp
           qqterm(1) = qqterm(1) + (1.-thimpb*bdf)*dt*temp
        endif

        ! NUMVAR = 2
        ! ~~~~~~~~~~
        if(numvar.ge.2) then
          temp = p1pv   (g79(:,:,i),g79(:,:,j),vzt79)
          ssterm = ssterm -     thimp     *dt*temp
          ddterm = ddterm + (1.-thimp*bdf)*dt*temp

          temp = p1pv(g79(:,:,i),p179,g79(:,:,j))
          rrterm(2) = rrterm(2) + thimpb*dt*temp
          qqterm(2) = qqterm(2) - thimpb*dt*temp*bdf

          if(eqsubtract.eq.1) then
             temp = p1pv  (g79(:,:,i),p079,g79(:,:,j))
             rrterm(2) = rrterm(2) +     thimpb     *dt*temp
             qqterm(2) = qqterm(2) + (1.-thimpb*bdf)*dt*temp
          endif
        endif

        ! NUMVAR = 3
        ! ~~~~~~~~~~
        if(numvar.ge.3) then
           temp = p1pchi    (g79(:,:,i),g79(:,:,j),cht79) &
                + b3pedkappa(g79(:,:,i),g79(:,:,j),ni79,kap79,hp) &
                + p1psipsikappar(g79(:,:,i),pst79,pst79,g79(:,:,j),ni79,b2i79,kar79) &
                + p1kappax  (g79(:,:,i),g79(:,:,j),bzt79,ni79,kax79)
           ssterm = ssterm -     thimp     *dt*temp
           ddterm = ddterm + (1.-thimp*bdf)*dt*temp
           
           temp = p1pchi(g79(:,:,i),p179,g79(:,:,j))
           rrterm(3) = rrterm(3) + thimpb*dt*temp
           qqterm(3) = qqterm(3) - thimpb*dt*temp*bdf

           temp = p1uus  (g79(:,:,i),g79(:,:,j),ph179,sig79) &
                + p1uus  (g79(:,:,i),ph179,g79(:,:,j),sig79) &
                + p1uchis(g79(:,:,i),g79(:,:,j),ch179,sig79) 
           rrterm(1) = rrterm(1) +     thimpb     *dt*temp
           qqterm(1) = qqterm(1) + (.5-thimpb*bdf)*dt*temp

           temp = p1vvs  (g79(:,:,i),g79(:,:,j),vz179,sig79) &
                + p1vvs  (g79(:,:,i),vz179,g79(:,:,j),sig79)
           rrterm(2) = rrterm(2) +     thimpb     *dt*temp
           qqterm(2) = qqterm(2) + (.5-thimpb*bdf)*dt*temp

           temp = p1chichis(g79(:,:,i),g79(:,:,j),ph179,sig79) &
                + p1chichis(g79(:,:,i),ph179,g79(:,:,j),sig79) &
                + p1uchis  (g79(:,:,i),ph079,g79(:,:,j),sig79) 
           rrterm(3) = rrterm(3) +     thimpb     *dt*temp
           qqterm(3) = qqterm(3) + (.5-thimpb*bdf)*dt*temp

           if(eqsubtract.eq.1) then
              temp = p1pchi(g79(:,:,i),p079,g79(:,:,j))
              rrterm(3) = rrterm(3) +     thimpb     *dt*temp
              qqterm(3) = qqterm(3) + (1.-thimpb*bdf)*dt*temp

              temp = p1uus  (g79(:,:,i),g79(:,:,j),ph079,sig79) &
                   + p1uus  (g79(:,:,i),ph079,g79(:,:,j),sig79) &
                   + p1uchis(g79(:,:,i),g79(:,:,j),ch079,sig79) 
              rrterm(1) = rrterm(1) +     thimpb     *dt*temp
              qqterm(1) = qqterm(1) + (1.-thimpb*bdf)*dt*temp

              temp = p1vvs  (g79(:,:,i),g79(:,:,j),vz079,sig79) &
                   + p1vvs  (g79(:,:,i),vz079,g79(:,:,j),sig79)
              rrterm(2) = rrterm(2) +     thimpb     *dt*temp
              qqterm(2) = qqterm(2) + (1.-thimpb*bdf)*dt*temp

              temp = p1chichis(g79(:,:,i),g79(:,:,j),ph079,sig79) &
                   + p1chichis(g79(:,:,i),ph079,g79(:,:,j),sig79) &
                   + p1uchis  (g79(:,:,i),ph079,g79(:,:,j),sig79) 
              rrterm(3) = rrterm(3) +     thimpb     *dt*temp
              qqterm(3) = qqterm(3) + (1.-thimpb*bdf)*dt*temp
           endif
        endif

        if(isplitstep.eq.0) rrterm = -rrterm

        call insval(pp1, ssterm, icomplex, ione+p_off, jone+p_off, 1)
        call insval(pp0, ddterm, icomplex, ione+p_off, jone+p_off, 1)
        call insval(pv1, rrterm(1), icomplex, i1+p_off, j1+u_off, 1)
        call insval(pv0, qqterm(1), icomplex, i1+p_off, j1+u_off, 1)
        if(numvar.ge.2) then
           call insval(pv1,rrterm(2), icomplex, i1+p_off,j1+vz_off,1)
           call insval(pv0,qqterm(2), icomplex, i1+p_off,j1+vz_off,1)
        endif
        if(numvar.ge.3) then
           call insval(pv1,rrterm(3), icomplex, i1+p_off,j1+chi_off,1)
           call insval(pv0,qqterm(3), icomplex, i1+p_off,j1+chi_off,1)
        endif

     enddo                     ! on j

     if(linear.eq.0) then
        psource(ione) = psource(ione) + dt* &
             (b3psipsieta(g79(:,:,i),pst79,pst79,eta79))
      
        if(numvar.ge.2) then
           psource(ione) = psource(ione) + dt* &
                (b3bbeta(g79(:,:,i),bzt79,bzt79,eta79))
        endif
        
        if(numvar.ge.3) then
           psource(ione) = psource(ione) + dt*dbf* &
                (b3pebd(g79(:,:,i),pet79,bzt79,ni79))
        endif

        if(eqsubtract.eq.1) then
           psource(ione) = psource(ione) + dt* &
                (b3pedkappa(g79(:,:,i),p079,ni79,kap79,hp) &
                +p1psipsikappar(g79(:,:,i),ps079,ps079,p079,ni79,b2i79,kar79) &
                +p1kappax  (g79(:,:,i),pe079,bz079,ni79,kax79) &
                +p1uus    (g79(:,:,i),ph079,ph079,sig79) &
                +p1vvs    (g79(:,:,i),vz079,vz079,sig79) &
                +p1chichis(g79(:,:,i),ch079,ch079,sig79) &
                +p1uchis  (g79(:,:,i),ph079,ch079,sig79))
           
           psource(ione) = psource(ione) + dt* &
                (p1pu   (g79(:,:,i),p079,ph079))
           
           if(numvar.ge.3) then
              psource(ione) = psource(ione) + dt* &
                   (p1pchi(g79(:,:,i),p079,ch079))
           endif
           
        endif
        
     endif  ! on linear.eq.0

  enddo                     ! on i
end subroutine ludefpres_n
