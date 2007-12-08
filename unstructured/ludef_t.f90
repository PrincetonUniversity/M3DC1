!======================================================================
! Vorticity Equation
!======================================================================
subroutine vorticity_lin(trial, lin, ssterm, ddterm, rrterm, qqterm, advfield)
  
  use basic
  use nintegrate_mod
#ifdef NEW_VELOCITY
  use metricterms_new
#else
  use metricterms_n
#endif

  implicit none

  vectype, dimension(79, OP_NUM), intent(in) :: trial, lin 
  vectype, dimension(3), intent(out) :: ssterm, ddterm
  vectype, dimension(4), intent(out) :: rrterm, qqterm
  integer, intent(in) :: advfield   ! if advfield = 1, eliminate rrterm by
                                    ! using analytic form of advanced field
  real :: parvisc
  real :: temp
  real :: ththm

  if(imp_mod.eq.1) then
     ththm = -bdf*thimp**2
  else
     ththm = (1.-thimp*bdf)*thimp
  endif

  ssterm = 0.
  ddterm = 0.
  rrterm = 0.
  qqterm = 0.


  ! NUMVAR = 1
  ! ~~~~~~~~~~
  temp = v1un(trial,lin,nt79)
  ssterm(1) = ssterm(1) + temp
  ddterm(1) = ddterm(1) + temp*bdf
  
  temp = v1umu(trial,lin,vis79) &
       + v1us (trial,lin,sig79)
  ssterm(1) = ssterm(1) -     thimp     *dt*temp
  ddterm(1) = ddterm(1) + (1.-thimp*bdf)*dt*temp
  
  temp = v1uun(trial,lin,ph179,nt79) &
       + v1uun(trial,ph179,lin,nt79)
  ssterm(1) = ssterm(1) -     thimp     *dt*temp
  ddterm(1) = ddterm(1) + (.5-thimp*bdf)*dt*temp

  if(gravr.ne.0 .or. gravz.ne.0) then          
     qqterm(4) = qqterm(4) + dt* &
          v1ngrav(trial,lin)
  endif
  
  if(gyro.eq.1) then
     temp = g1u(trial,lin)*dbf
     ssterm(1) = ssterm(1) +     thimp     *dt*temp
     ddterm(1) = ddterm(1) - (1.-thimp*bdf)*dt*temp
  endif

  if(advfield.eq.1) then
     temp = v1upsipsi(trial,lin,pst79,pst79) 
     if(gravr.ne.0 .or. gravz.ne.0) then          
        temp = temp + v1ungrav(trial,lin,nt79)
     endif
     ssterm(1) = ssterm(1) - thimp*thimp*dt*dt*temp
     ddterm(1) = ddterm(1) +       ththm*dt*dt*temp

     qqterm(1) = qqterm(1)  + dt* &
          (v1psipsi(trial,lin,pss79)  &
          +v1psipsi(trial,pss79,lin))
  
     if(isources.eq.1) then
        qqterm(1) = qqterm(1) + thimp*dt*dt* &
             (v1psisb1 (trial,lin,sb179))
     endif
  else
     temp = v1psipsi(trial,lin,ps179) &
          + v1psipsi(trial,ps179,lin)
     rrterm(1) = rrterm(1) +     thimp     *dt*temp
     qqterm(1) = qqterm(1) + (.5-thimp*bdf)*dt*temp

     if(idens.eq.1 .and. (gravr.ne.0 .or. gravz.ne.0)) then
        temp = v1ngrav(trial,lin)
        rrterm(4) = rrterm(4) +     thimp     *dt*temp
        qqterm(4) = qqterm(4) + (1.-thimp*bdf)*dt*temp
     endif
  endif
    
  if(linear.eq.1 .or. eqsubtract.eq.1) then
     temp = v1uun(trial,lin,ph079,nt79) &    
          + v1uun(trial,ph079,lin,nt79)
     ssterm(1) = ssterm(1) -     thimp     *dt*temp
     ddterm(1) = ddterm(1) + (1.-thimp*bdf)*dt*temp
     
     if(advfield.eq.1) then
        qqterm(1) = qqterm(1) + thimp*dt*dt* &
             (v1upsipsi(trial,ph079,lin,pss79) &
             +v1upsipsi(trial,ph079,pss79,lin))

        if(idens.eq.1) then
           qqterm(4) = qqterm(4) + dt*              &
                (v1un     (trial,ph079,lin)         &
                +v1uun    (trial,ph079,ph079,lin))
           if(gravr.ne.0 .or. gravz.ne.0) then 
              qqterm(4) = qqterm(4) + thimp*dt*dt*  &
                   (v1ungrav   (trial,ph079,lin)    &
                   +v1ndenmgrav(trial,lin, denm))
           endif
        endif
     else
        temp = v1psipsi(trial,lin,ps079) &
             + v1psipsi(trial,ps079,lin)
        rrterm(1) = rrterm(1) +     thimp     *dt*temp
        qqterm(1) = qqterm(1) + (1.-thimp*bdf)*dt*temp

        if(idens.eq.1) then
           qqterm(4) = qqterm(4) + dt* &
                v1uun(trial,ph079,ph079,lin)
        endif
     endif
  endif
  
  
  ! NUMVAR = 2
  ! ~~~~~~~~~~
  if(numvar.eq.1) then
     ! These terms exist since F0=bzero when numvar=1
     if(advfield.eq.1) then
        temp = v1upsib(trial,lin,pst79,bzt79) &
             + v1ubb  (trial,lin,bzt79,bzt79)
        ssterm(1) = ssterm(1) - thimp*thimp*dt*dt*temp
        ddterm(1) = ddterm(1) +       ththm*dt*dt*temp
        
        qqterm(1) = qqterm(1) + dt* &
             (v1psib  (trial,lin,bzt79))
     else
        temp = v1psib(trial,lin,bzt79)
        rrterm(1) = rrterm(1) +     thimp     *dt*temp
        qqterm(1) = qqterm(1) + (1.-thimp*bdf)*dt*temp
     endif

  else
          
     temp = v1vvn(trial,lin,vz179,nt79) &
          + v1vvn(trial,vz179,lin,nt79)
     ssterm(2) = ssterm(2) -     thimp     *dt*temp
     ddterm(2) = ddterm(2) + (.5-thimp*bdf)*dt*temp
     
     if(gyro.eq.1) then
        temp = g1v(trial,lin)*dbf
        ssterm(2) = ssterm(2) +     thimp     *dt*temp
        ddterm(2) = ddterm(2) - (1.-thimp*bdf)*dt*temp
     endif

     if(advfield.eq.1) then
        temp = v1ubb  (trial,lin,bzt79,bzt79)
        ssterm(1) = ssterm(1) - thimp*thimp*dt*dt*temp
        ddterm(1) = ddterm(1) +       ththm*dt*dt*temp

        temp = v1vpsib(trial,lin,pst79,bzt79)
        ssterm(2) = ssterm(2) - thimp*thimp*dt*dt*temp
        ddterm(2) = ddterm(2) +       ththm*dt*dt*temp

        qqterm(2) = qqterm(2) + dt* &
             (v1bb     (trial,lin,bzs79)  &
             +v1bb     (trial,bzs79,lin))
     
        if(isources.eq.1) then
           qqterm(2) = qqterm(2) + thimp*dt*dt* &
                (v1bsb2   (trial,lin,sb279))
        endif
     else
        temp = v1bb(trial,lin,bz179) &
             + v1bb(trial,bz179,lin)
        rrterm(2) = rrterm(2) +     thimp     *dt*temp
        qqterm(2) = qqterm(2) + (.5-thimp*bdf)*dt*temp
     endif
          
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        
        temp = v1vvn(trial,lin,vz079,nt79) &
             + v1vvn(trial,vz079,lin,nt79)
        ssterm(2) = ssterm(2) -     thimp     *dt*temp
        ddterm(2) = ddterm(2) + (1.-thimp*bdf)*dt*temp

        if(advfield.eq.1) then
           qqterm(1) = qqterm(1) + thimp*dt*dt* &
                (v1vpsib  (trial,vz079,lin,bzs79))
        
           qqterm(2) = qqterm(2) + thimp*dt*dt* &
                (v1ubb    (trial,ph079,lin,bzs79) &
                +v1ubb    (trial,ph079,bzs79,lin) &
                +v1vpsib  (trial,vz079,pss79,lin))

           if(idens.eq.1) then
              qqterm(4) = qqterm(4) + dt*       &
                   v1vvn(trial,vz079,vz079,lin)
           endif
        else
           temp = v1bb(trial,lin,bz079) &
                + v1bb(trial,bz079,lin)
           rrterm(2) = rrterm(2) +     thimp     *dt*temp
           qqterm(2) = qqterm(2) + (1.-thimp*bdf)*dt*temp

           if(idens.eq.1) then
              qqterm(4) = qqterm(4) + dt* &
                   v1vvn(trial,vz079,vz079,lin)
           endif
        endif
     endif
  endif
  
  ! NUMVAR = 3
  ! ~~~~~~~~~~
  if(numvar.ge.3) then
     
     temp = v1uchin(trial,lin,ch179,nt79)       
     ssterm(1) = ssterm(1) -     thimp     *dt*temp
     ddterm(1) = ddterm(1) + (.5-thimp*bdf)*dt*temp
     
     temp = v1chin(trial,lin,nt79)
     ssterm(3) = ssterm(3) + temp
     ddterm(3) = ddterm(3) + temp*bdf
     
     temp = v1chimu(trial,lin,vis79) &
          + v1chis (trial,lin,sig79)
     ssterm(3) = ssterm(3) -     thimp     *dt*temp
     ddterm(3) = ddterm(3) + (1.-thimp*bdf)*dt*temp
     
     temp = v1uchin  (trial,ph179,lin,nt79) &
          + v1chichin(trial,lin,ch179,nt79) &
          + v1chichin(trial,ch179,lin,nt79)
     ssterm(3) = ssterm(3) -     thimp     *dt*temp
     ddterm(3) = ddterm(3) + (.5-thimp*bdf)*dt*temp
          
     if(gyro.eq.1) then
        temp = g1chi(trial,lin)*dbf
        ssterm(3) = ssterm(3) +     thimp     *dt*temp
        ddterm(3) = ddterm(3) - (1.-thimp*bdf)*dt*temp
     endif

     if(advfield.eq.1) then
        temp = v1chipsipsi(trial,lin,pst79,pst79)  &
             + v1chibb    (trial,lin,bzt79,bzt79)
        if(gravr.ne.0 .or. gravz.ne.0) then          
           temp = temp &
                + v1chingrav(trial,lin,nt79)
        endif
        ssterm(3) = ssterm(3) - thimp*thimp*dt*dt*temp
        ddterm(3) = ddterm(3) +       ththm*dt*dt*temp
     endif
     
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        temp = v1uchin(trial,lin,ch079,nt79)
        ssterm(1) = ssterm(1) -     thimp     *dt*temp
        ddterm(1) = ddterm(1) + (1.-thimp*bdf)*dt*temp
        
        temp = v1uchin(trial,ph079,lin,nt79)
        ssterm(3) = ssterm(3) -     thimp     *dt*temp
        ddterm(3) = ddterm(3) + (1.-thimp*bdf)*dt*temp
        
        if(advfield.eq.1) then
           qqterm(1) = qqterm(1) + thimp*dt*dt* &
                (v1chipsipsi(trial,ch079,lin,pss79) &
                +v1chipsipsi(trial,ch079,pss79,lin))
        
           qqterm(2) = qqterm(2) + thimp*dt*dt* &
                (v1chibb(trial,ch079,lin,bzs79) &
                +v1chibb(trial,ch079,bzs79,lin))

           if(idens.eq.1) then
              qqterm(4) = qqterm(4) + dt* &
                   (v1uchin  (trial,ph079,ch079,lin) &
                   +v1chichin(trial,ch079,ch079,lin))
              if(gravr.ne.0 .or. gravz.ne.0) then          
                 qqterm(4) = qqterm(4) + thimp*dt*dt* &
                      v1chingrav (trial,ch079,lin)
              endif
           endif
        else
           if(idens.eq.1) then
              qqterm(4) = qqterm(4) + dt* &
                   (v1uchin  (trial,vz079,ch079,lin) &
                   +v1chichin(trial,ch079,ch079,lin))
           endif
        endif
     endif
  endif

  ! Parallel Viscosity
  if(amupar.ne.0) then
     call PVS1(lin,temp79b)
     parvisc = int2(vic79(:,OP_1),temp79b,weight_79,79)

     temp = PVV1(trial)*parvisc
     ssterm(1) = ssterm(1) +     thimp     *dt*temp
     ddterm(1) = ddterm(1) - (1.-thimp*bdf)*dt*temp

     if(numvar.ge.2) then
        temp = PVV2(trial)*parvisc
        ssterm(2) = ssterm(2) +     thimp     *dt*temp
        ddterm(2) = ddterm(2) - (1.-thimp*bdf)*dt*temp
     endif

     if(numvar.ge.3) then
        temp = PVV3(trial)*parvisc
        ssterm(3) = ssterm(3) +     thimp     *dt*temp
        ddterm(3) = ddterm(3) - (1.-thimp*bdf)*dt*temp
     endif
  endif

end subroutine vorticity_lin 

subroutine vorticity_nolin(trial, r4term)

  use basic
  use nintegrate_mod
#ifdef NEW_VELOCITY
  use metricterms_new
#else
  use metricterms_n
#endif

  implicit none

  vectype, intent(in), dimension(79, OP_NUM)  :: trial
  vectype, intent(out) :: r4term
  vectype :: temp

  r4term = 0.

  ! density terms
  ! ~~~~~~~~~~~~~
  if(idens.eq.1 .and. (linear.eq.1 .or. eqsubtract.eq.1)) then
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
subroutine axial_vel_lin(trial, lin, ssterm, ddterm, rrterm, qqterm, advfield)
  
  use basic
  use nintegrate_mod
  use metricterms_new

  implicit none

  vectype, dimension(79, OP_NUM), intent(in) :: trial, lin 
  vectype, dimension(3), intent(out) :: ssterm, ddterm
  vectype, dimension(4), intent(out) :: rrterm, qqterm
  integer, intent(in) :: advfield

  real :: temp, parvisc
  real :: ththm

  if(imp_mod.eq.1) then
     ththm = -bdf*thimp**2
  else
     ththm = (1.-thimp*bdf)*thimp
  endif

  ssterm = 0.
  ddterm = 0.
  rrterm = 0.
  qqterm = 0.

  if(numvar.lt.2) return
         
  temp = v2vun(trial,vz179,lin,nt79)
  ssterm(1) = ssterm(1) -     thimp     *dt*temp
  ddterm(1) = ddterm(1) + (.5-thimp*bdf)*dt*temp
          
  temp = v2vn(trial,lin,nt79)
  ssterm(2) = ssterm(2) + temp
  ddterm(2) = ddterm(2) + temp*bdf
     
  temp = v2vmu  (trial,lin,vis79,hypv*sz79) &
       + v2vs   (trial,lin,sig79)
  ssterm(2) = ssterm(2) -     thimp     *dt*temp
  ddterm(2) = ddterm(2) + (1.-thimp*bdf)*dt*temp

  temp = v2vun(trial,lin,ph179,nt79)
  ssterm(2) = ssterm(2) -      thimp     *dt*temp
  ddterm(2) = ddterm(2) + (0.5-thimp*bdf)*dt*temp

  if(gyro.eq.1) then
     temp = g2u(trial,lin)*dbf
     ssterm(1) = ssterm(1) +     thimp     *dt*temp
     ddterm(1) = ddterm(1) - (1.-thimp*bdf)*dt*temp
     
     temp = g2v(trial,lin)*dbf
     ssterm(2) = ssterm(2) +     thimp     *dt*temp
     ddterm(2) = ddterm(2) - (1.-thimp*bdf)*dt*temp
  endif
          
  if(advfield.eq.1) then 
     temp = v2upsib(trial,lin,pst79,bzt79)
     ssterm(1) = ssterm(1) - thimp*thimp*dt*dt*temp
     ddterm(1) = ddterm(1) +       ththm*dt*dt*temp

     temp = v2vpsipsi(trial,lin,pst79,pst79)
     ssterm(2) = ssterm(2) - thimp*thimp*dt*dt*temp
     ddterm(2) = ddterm(2) +       ththm*dt*dt*temp

     qqterm(1) = qqterm(1) + dt* &
          (v2psib(trial,lin,bzs79))
     
     qqterm(2) = qqterm(2) + dt* &
          (v2psib(trial,pss79,lin))
     
     if(isources.eq.1) then
        qqterm(1) = qqterm(1) + thimp*dt*dt* &
             (v2psisb2(trial,lin,sb279))
        qqterm(2) = qqterm(2) + thimp*dt*dt* &
             (v2bsb1(trial,lin,sb179))
     endif
  else
     temp = v2psib(trial,lin,bz179)
     rrterm(1) = rrterm(1) +     thimp     *dt*temp
     qqterm(1) = qqterm(1) + (.5-thimp*bdf)*dt*temp

     temp = v2psib(trial,ps179,lin)
     rrterm(2) = rrterm(2) +     thimp     *dt*temp
     qqterm(2) = qqterm(2) + (.5-thimp*bdf)*dt*temp
  endif

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     
     temp = v2vun(trial,vz079,lin,nt79)
     ssterm(1) = ssterm(1) -     thimp     *dt*temp
     ddterm(1) = ddterm(1) + (1.-thimp*bdf)*dt*temp
              
     temp = v2vun(trial,lin,ph079,nt79)
     ssterm(2) = ssterm(2) -     thimp     *dt*temp
     ddterm(2) = ddterm(2) + (1.-thimp*bdf)*dt*temp

     if(advfield.eq.1) then
        qqterm(1) = qqterm(1) + thimp*dt*dt* &
             (v2upsib  (trial,ph079,lin,bzs79) &
             +v2vpsipsi(trial,vz079,lin,pss79) &
             +v2vpsipsi(trial,vz079,pss79,lin))
        
        qqterm(2) = qqterm(2) + thimp*dt*dt* &
             (v2upsib  (trial,ph079,pss79,lin))        
     else
        temp = v2psib(trial,lin,bz079)
        rrterm(1) = rrterm(1) +     thimp     *dt*temp
        qqterm(1) = qqterm(1) + (1.-thimp*bdf)*dt*temp

        temp = v2psib(trial,ps079,lin)
        rrterm(2) = rrterm(2) +     thimp     *dt*temp
        qqterm(2) = qqterm(2) + (1.-thimp*bdf)*dt*temp

        if(idens.eq.1) then
           qqterm(4) = qqterm(4) + dt* &
                v2vun(trial,vz079,ph079,lin)
        endif
     end if
  endif
        
  ! NUMVAR = 3
  ! ~~~~~~~~~~
  if(numvar.ge.3) then

     temp = v2vchin(trial,lin,ch179,nt79)
     ssterm(2) = ssterm(2) -     thimp     *dt*temp
     ddterm(2) = ddterm(2) + (.5-thimp*bdf)*dt*temp
           
     temp = v2vchin(trial,vz179,lin,nt79)
     ssterm(3) = ssterm(3) -     thimp     *dt*temp
     ddterm(3) = ddterm(3) + (.5-thimp*bdf)*dt*temp

     if(gyro.eq.1) then
        temp = g2chi(trial,lin)*dbf
        ssterm(3) = ssterm(3) +     thimp     *dt*temp
        ddterm(3) = ddterm(3) - (1.-thimp*bdf)*dt*temp
     endif

     if(advfield.eq.1) then
        temp = v2chipsib(trial,lin,pst79,bzt79) 
        ssterm(3) = ssterm(3) - thimp*thimp*dt*dt*temp
        ddterm(3) = ddterm(3) +       ththm*dt*dt*temp
     end if
           
     if(linear.eq.1 .or. eqsubtract.eq.1) then              
        temp = v2vchin(trial,lin,ch079,nt79)
        ssterm(2) = ssterm(2) -     thimp     *dt*temp
        ddterm(2) = ddterm(2) + (1.-thimp*bdf)*dt*temp
              
        temp = v2vchin(trial,vz079,lin,nt79)
        ssterm(3) = ssterm(3) -     thimp     *dt*temp
        ddterm(3) = ddterm(3) + (1.-thimp*bdf)*dt*temp

        if(advfield.eq.1) then
           qqterm(1) = qqterm(1) + thimp*dt*dt* &
                (v2chipsib(trial,ch079,lin,bzs79))
              
           qqterm(2) = qqterm(2) + thimp*dt*dt* &
                (v2chipsib(trial,ch079,pss79,lin))
        else
           if(idens.eq.1) then
              qqterm(4) = qqterm(4) + dt* &
                   v2vchin(trial,vz079,ch079,lin)
           endif
        endif
     endif
  endif

  ! Parallel Viscosity
  if(amupar.ne.0) then
     call PVS2(lin,temp79b)
     parvisc = int2(vic79(:,OP_1),temp79b,weight_79,79)

     temp = PVV1(trial)*parvisc
     ssterm(1) = ssterm(1) +     thimp     *dt*temp
     ddterm(1) = ddterm(1) - (1.-thimp*bdf)*dt*temp

     if(numvar.ge.2) then
        temp = PVV2(trial)*parvisc
        ssterm(2) = ssterm(2) +     thimp     *dt*temp
        ddterm(2) = ddterm(2) - (1.-thimp*bdf)*dt*temp
     endif

     if(numvar.ge.3) then
        temp = PVV3(trial)*parvisc
        ssterm(3) = ssterm(3) +     thimp     *dt*temp
        ddterm(3) = ddterm(3) - (1.-thimp*bdf)*dt*temp
     endif
  endif

end subroutine axial_vel_lin

subroutine axial_vel_nolin(trial, r4term)

  use basic
  use nintegrate_mod
  use metricterms_new

  implicit none

  vectype, intent(in), dimension(79, OP_NUM)  :: trial
  vectype, intent(out) :: r4term
  vectype :: temp
  
  r4term = 0.

  if(numvar.lt.2) return

  ! density terms
  ! ~~~~~~~~~~~~~
  if(idens.eq.1 .and. (linear.eq.1 .or. eqsubtract.eq.1)) then
     r4term = r4term + dt* &
          (v2vun(trial,vz079,ph079,n179) &
          +v2vs (trial,vz079,sig79))
  endif

end subroutine axial_vel_nolin


!======================================================================
! Compression Equation
!======================================================================
subroutine compression_lin(trial, lin, ssterm, ddterm, rrterm, qqterm, advfield)
  
  use basic
  use nintegrate_mod
  use metricterms_new

  implicit none

  vectype, dimension(79, OP_NUM), intent(in) :: trial, lin 
  vectype, dimension(3), intent(out) :: ssterm, ddterm
  vectype, dimension(4), intent(out) :: rrterm, qqterm
  integer, intent(in) :: advfield

  real :: temp, parvisc
  real :: ththm

  if(imp_mod.eq.1) then
     ththm = -bdf*thimp**2
  else
     ththm = (1.-thimp*bdf)*thimp
  endif

  ssterm = 0.
  ddterm = 0.
  rrterm = 0.
  qqterm = 0.

  if(numvar.lt.3) return

  ! regularize the chi equation
  temp = -regular*int2(trial(:,OP_1),lin(:,OP_1),weight_79,79)
  ssterm(3) = ssterm(3) + temp
  ddterm(3) = ddterm(3) + temp*bdf
         
  temp = v3un(trial,lin,nt79)
  ssterm(1) = ssterm(1) + temp
  ddterm(1) = ddterm(1) + temp*bdf

  temp = v3uun  (trial,lin,ph179,nt79) &
       + v3uun  (trial,ph179,lin,nt79) &
       + v3uchin(trial,lin,ch179,nt79)
  ssterm(1) = ssterm(1) -     thimp     *dt*temp
  ddterm(1) = ddterm(1) + (.5-thimp*bdf)*dt*temp
           
  temp = v3umu   (trial,lin,vis79,vic79) &
       + v3us    (trial,lin,sig79)
  ssterm(1) = ssterm(1) -     thimp     *dt*temp
  ddterm(1) = ddterm(1) + (1.-thimp*bdf)*dt*temp
                     
  temp = v3vvn(trial,lin,vz179,nt79) &
       + v3vvn(trial,vz179,lin,nt79)
  ssterm(2) = ssterm(2) -     thimp     *dt*temp
  ddterm(2) = ddterm(2) + (.5-thimp*bdf)*dt*temp
  
  temp = v3chin(trial,lin,nt79)
  ssterm(3) = ssterm(3) + temp
  ddterm(3) = ddterm(3) + temp*bdf

  temp = v3chimu(trial,lin,vis79,vic79) &
       + v3chis (trial,lin,sig79)
  ssterm(3) = ssterm(3) -     thimp     *dt*temp
  ddterm(3) = ddterm(3) + (1.-thimp*bdf)*dt*temp
  
  temp = v3uchin  (trial,ph179,lin,nt79) &
       + v3chichin(trial,lin,ch179,nt79) &
       + v3chichin(trial,ch179,lin,nt79)  
  ssterm(3) = ssterm(3) -     thimp     *dt*temp
  ddterm(3) = ddterm(3) + (.5-thimp*bdf)*dt*temp

  if(gravr.ne.0 .or. gravz.ne.0) then          
     qqterm(4) = qqterm(4) + dt* &
          v3ngrav(trial,lin)
  endif
    
  if(gyro.eq.1) then    
     temp = g3u(trial,lin)*dbf
     ssterm(1) = ssterm(1) +     thimp     *dt*temp
     ddterm(1) = ddterm(1) - (1.-thimp*bdf)*dt*temp
     
     temp = g3v(trial,lin)*dbf
     ssterm(2) = ssterm(2) +     thimp     *dt*temp
     ddterm(2) = ddterm(2) - (1.-thimp*bdf)*dt*temp
     
     temp = g3chi(trial,lin)*dbf
     ssterm(3) = ssterm(3) +     thimp     *dt*temp
     ddterm(3) = ddterm(3) - (1.-thimp*bdf)*dt*temp
  endif
  
  if(advfield.eq.1) then
     temp = v3up     (trial,lin,pt79)        &
          + v3upsipsi(trial,lin,pst79,pst79) &
          + v3ubb    (trial,lin,bzt79,bzt79) 
     if(gravr.ne.0 .or. gravz.ne.0) then
        temp = temp + thimp*dt* &
             v3ungrav(trial,lin,nt79)
     endif
     ssterm(1) = ssterm(1) - thimp*thimp*dt*dt*temp
     ddterm(1) = ddterm(1) +       ththm*dt*dt*temp

     temp = v3vpsib(trial,lin,pst79,bzt79)   
     ssterm(2) = ssterm(2) - thimp*thimp*dt*dt*temp
     ddterm(2) = ddterm(2) +       ththm*dt*dt*temp

     temp = v3chip     (trial,lin,pt79)        &
          + v3chipsipsi(trial,lin,pst79,pst79) &
          + v3chibb    (trial,lin,bzt79,bzt79)
     if(gravr.ne.0 .or. gravz.ne.0) then
        temp = temp + &
             v3chingrav(trial,lin,nt79)
     endif
     ssterm(3) = ssterm(3) - thimp*thimp*dt*dt*temp
     ddterm(3) = ddterm(3) +       ththm*dt*dt*temp

     qqterm(1) = qqterm(1) +            &
          dt*                           &
          (v3psipsi(trial,lin,pss79)    & 
          +v3psipsi(trial,pss79,lin))
     
     qqterm(2) = qqterm(2) + dt*        &
          (v3bb(trial,lin,bzs79)        &     
          +v3bb(trial,bzs79,lin)) 
  
     qqterm(3) = qqterm(3) + dt*        &
          v3p(trial,lin)

     if(isources.eq.1) then
        qqterm(1) = qqterm(1) +            &
             thimp*dt*dt*                  &
             (v3psisb1(trial,lin,sb179))
        
        qqterm(2) = qqterm(2) +            &
             thimp*dt*dt*                  &
             (v3bsb2(trial,lin,sb279))
     end if
  else
     temp = v3psipsi(trial,lin,ps179) &
          + v3psipsi(trial,ps179,lin)
     rrterm(1) = rrterm(1) +     thimp     *dt*temp
     qqterm(1) = qqterm(1) + (.5-thimp*bdf)*dt*temp

     temp = v3bb(trial,lin,bz179) &
          + v3bb(trial,bz179,lin)
     rrterm(2) = rrterm(2) +     thimp     *dt*temp
     qqterm(2) = qqterm(2) + (.5-thimp*bdf)*dt*temp

     temp = v3p(trial,lin) 
     rrterm(3) = rrterm(3) +     thimp     *dt*temp
     qqterm(3) = qqterm(3) + (1.-thimp*bdf)*dt*temp

     if(idens.eq.1) then
        temp = v3ngrav(trial,lin)
        rrterm(4) = rrterm(4) +     thimp     *dt*temp
        qqterm(4) = qqterm(4) + (1.-thimp*bdf)*dt*temp
     endif
  endif
  
  if(linear.eq.1 .or. eqsubtract.eq.1) then             
     temp = v3uun  (trial,lin,ph079,nt79) &
          + v3uun  (trial,ph079,lin,nt79) &
          + v3uchin(trial,lin,ch079,nt79)
     ssterm(1) = ssterm(1) -     thimp     *dt*temp
     ddterm(1) = ddterm(1) + (1.-thimp*bdf)*dt*temp
     
     temp = v3vvn(trial,lin,vz079,nt79) &
          + v3vvn(trial,vz079,lin,nt79)
     ssterm(2) = ssterm(2) -     thimp     *dt*temp
     ddterm(2) = ddterm(2) + (1.-thimp*bdf)*dt*temp
     
     temp = v3uchin  (trial,ph079,lin,nt79) &
          + v3chichin(trial,lin,ch079,nt79) &
          + v3chichin(trial,ch079,lin,nt79)
     ssterm(3) = ssterm(3) -     thimp     *dt*temp
     ddterm(3) = ddterm(3) + (1.-thimp*bdf)*dt*temp

     if(advfield.eq.1) then
        qqterm(1) = qqterm(1) + thimp*dt*dt* &
             (v3upsipsi  (trial,ph079,lin,pss79) &
             +v3upsipsi  (trial,ph079,pss79,lin) &
             +v3vpsib    (trial,vz079,lin,bzs79) &
             +v3chipsipsi(trial,ch079,lin,pss79) &
             +v3chipsipsi(trial,ch079,pss79,lin))
     
        qqterm(2) = qqterm(2) + thimp*dt*dt* &
             (v3ubb  (trial,ph079,lin,bzs79) &
             +v3ubb  (trial,ph079,bzs79,lin) &
             +v3vpsib(trial,vz079,pss79,lin) &
             +v3chibb(trial,ch079,lin,bzs79) &
             +v3chibb(trial,ch079,bzs79,lin))
     
        qqterm(3) = qqterm(3) + thimp*dt*dt* &
             (v3up  (trial,ph079,lin) &
             +v3chip(trial,ch079,lin))
     else
        temp = v3psipsi(trial,lin,ps079) &
             + v3psipsi(trial,ps079,lin)
        rrterm(1) = rrterm(1) +     thimp     *dt*temp
        qqterm(1) = qqterm(1) + (1.-thimp*bdf)*dt*temp

        temp = v3bb(trial,lin,bz079) &
             + v3bb(trial,bz079,lin)
        rrterm(2) = rrterm(2) +     thimp     *dt*temp
        qqterm(2) = qqterm(2) + (1.-thimp*bdf)*dt*temp

        if(idens.eq.1) then
           qqterm(4) = qqterm(4) + dt*             &
                (v3uun    (trial,ph079,ph079,lin)  &
                +v3vvn    (trial,vz079,vz079,lin)  &
                +v3uchin  (trial,ph079,ch079,lin)  &
                +v3chichin(trial,ch079,ch079,lin))
        end if
     end if
  endif

  ! Parallel Viscosity
  if(amupar.ne.0) then
     call PVS3(lin,temp79b)
     parvisc = int2(vic79(:,OP_1),temp79b,weight_79,79)

     temp = PVV1(trial)*parvisc
     ssterm(1) = ssterm(1) +     thimp     *dt*temp
     ddterm(1) = ddterm(1) - (1.-thimp*bdf)*dt*temp

     if(numvar.ge.2) then
        temp = PVV2(trial)*parvisc
        ssterm(2) = ssterm(2) +     thimp     *dt*temp
        ddterm(2) = ddterm(2) - (1.-thimp*bdf)*dt*temp
     endif

     if(numvar.ge.3) then
        temp = PVV3(trial)*parvisc
        ssterm(3) = ssterm(3) +     thimp     *dt*temp
        ddterm(3) = ddterm(3) - (1.-thimp*bdf)*dt*temp
     endif
  endif

end subroutine compression_lin

subroutine compression_nolin(trial, r4term)

  use basic
  use nintegrate_mod
  use metricterms_new

  implicit none

  vectype, intent(in), dimension(79, OP_NUM)  :: trial
  vectype, intent(out) :: r4term
  vectype :: temp
  
  r4term = 0.

  if(numvar.lt.3) return

  ! density terms
  ! ~~~~~~~~~~~~~
  if(idens.eq.1 .and. (linear.eq.1 .or. eqsubtract.eq.1)) then
     if(gravr.ne.0 .or. gravz.ne.0) then          
        r4term = r4term + thimp*dt*dt* &
             (v3ungrav   (trial,ph079,n179) &
             +v3chingrav (trial,ch079,n179) &
             +v3ndenmgrav(trial,n179, denm))
     endif
     r4term = r4term + dt* &
          (v3uun    (trial,ph079,ph079,n179) &
          +v3uchin  (trial,ph079,ch079,n179) &
          +v3vvn    (trial,vz079,vz079,n179) &
          +v3chichin(trial,ch079,ch079,n179) &
          +v3us     (trial,ph079,sig79) &
          +v3chis   (trial,ch079,sig79))
  endif

end subroutine compression_nolin


!======================================================================
! Flux Equation
!======================================================================
subroutine flux_lin(trial, lin, ssterm, ddterm, rrterm, qqterm)
  
  use basic
  use nintegrate_mod
  use metricterms_new

  implicit none

  vectype, dimension(79, OP_NUM), intent(in) :: trial, lin 
  vectype, dimension(3), intent(out) :: ssterm, ddterm, rrterm, qqterm
  real :: temp, thimpb

  if(imp_mod.eq.1) then
     thimpb = 1.
  else
     thimpb = thimp
  endif

  ssterm = 0.
  ddterm = 0.
  rrterm = 0.
  qqterm = 0.

  temp = b1psi(trial,lin)
  ssterm(1) = ssterm(1) + temp
  ddterm(1) = ddterm(1) + temp*bdf

  temp = b1psieta(trial,lin,eta79,hypf*sz79) &
       + b1psiu  (trial,lin,pht79)
  ssterm(1) = ssterm(1) -     thimpb     *dt*temp
  ddterm(1) = ddterm(1) + (1.-thimpb*bdf)*dt*temp

  temp = (b1psipsid(trial,lin,ps179,ni79) &
         +b1psipsid(trial,ps179,lin,ni79))*dbf
  ssterm(1) = ssterm(1) -     thimpb     *dt*temp
  ddterm(1) = ddterm(1) + (.5-thimpb*bdf)*dt*temp

  temp = b1psiu(trial,ps179,lin)
  rrterm(1) = rrterm(1) + thimpb*dt*temp
  qqterm(1) = qqterm(1) - thimpb*dt*temp*bdf

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     temp = (b1psipsid(trial,lin,ps079,ni79) &
            +b1psipsid(trial,ps079,lin,ni79))*dbf
     ssterm(1) = ssterm(1) -     thimpb     *dt*temp
     ddterm(1) = ddterm(1) + (1.-thimpb*bdf)*dt*temp

     temp = b1psiu(trial,ps079,lin)
     rrterm(1) = rrterm(1) +     thimpb     *dt*temp
     qqterm(1) = qqterm(1) + (1.-thimpb*bdf)*dt*temp
  endif

  ! NUMVAR = 2
  ! ~~~~~~~~~~
  if(numvar.eq.1) then
     ! Term due to F0 = bzero when numvar=1
     temp = b1psibd(trial,lin,bzt79,ni79)*dbf
     ssterm(1) = ssterm(1) -     thimpb     *dt*temp
     ddterm(1) = ddterm(1) + (1.-thimpb*bdf)*dt*temp

     temp = b1bu(trial,bzt79,lin) 
     rrterm(1) = rrterm(1) + thimpb*dt*temp
     qqterm(1) = qqterm(1) - thimpb*dt*temp*bdf           
  else           
     temp = b1psibd(trial,lin,bz179,ni79)*dbf
     ssterm(1) = ssterm(1) -     thimpb     *dt*temp
     ddterm(1) = ddterm(1) + (.5-thimpb*bdf)*dt*temp
           
     temp = b1psibd(trial,ps179,lin,ni79)*dbf              
     ssterm(2) = ssterm(2) -     thimpb     *dt*temp
     ddterm(2) = ddterm(2) + (.5-thimpb*bdf)*dt*temp
           
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        temp = b1psibd(trial,lin,bz079,ni79)*dbf
        ssterm(1) = ssterm(1) -     thimpb     *dt*temp
        ddterm(1) = ddterm(1) + (1.-thimpb*bdf)*dt*temp
              
        temp = b1psibd(trial,ps079,lin,ni79)*dbf
        ssterm(2) = ssterm(2) -     thimpb     *dt*temp
        ddterm(2) = ddterm(2) + (1.-thimpb*bdf)*dt*temp
     end if
  end if

  ! NUMVAR = 3
  ! ~~~~~~~~~~
  if(numvar.ge.3) then

     temp = b1psichi(trial,lin,cht79)                
     ssterm(1) = ssterm(1) -     thimpb     *dt*temp
     ddterm(1) = ddterm(1) + (1.-thimpb*bdf)*dt*temp

     temp = b1psichi(trial,ps179,lin)                
     rrterm(3) = rrterm(3) + thimpb*dt*temp
     qqterm(3) = qqterm(3) - thimpb*dt*temp*bdf
     
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        temp = b1psichi(trial,ps079,lin)
        rrterm(3) = rrterm(3) +     thimpb     *dt*temp
        qqterm(3) = qqterm(3) + (1.-thimpb*bdf)*dt*temp
     endif
  end if
        
end subroutine flux_lin

subroutine flux_nolin(trial, r4term)

  use basic
  use nintegrate_mod
  use metricterms_new

  implicit none

  vectype, intent(in), dimension(79, OP_NUM)  :: trial
  vectype, intent(out) :: r4term
  vectype :: temp
  
  r4term = 0.

  if(igauge.eq.1) then
     r4term = r4term - dt* &
          vloop*int1(trial,weight_79,79)/(2.*pi)
  endif

  ! density terms
  ! ~~~~~~~~~~~~~
  if(idens.eq.1 .and. (linear.eq.1 .or. eqsubtract.eq.1)) then
     r4term = r4term + dt* &
          (b1psipsid(trial,ps079,ps079,ni79)*dbf)

     if(numvar.eq.1) then
        ! Term due to F0 = bzero when numvar=1
        r4term = r4term + dt* &
             (b1psibd  (trial,ps079,bzt79,ni79)*dbf)
     else
        r4term = r4term + dt* &
             (b1psibd  (trial,ps079,bz079,ni79)*dbf)
     endif
  endif
  
end subroutine flux_nolin


!======================================================================
! Axial Field Equation
!======================================================================
subroutine axial_field_lin(trial, lin, ssterm, ddterm, rrterm, qqterm)
  
  use basic
  use nintegrate_mod
#ifdef NEW_VELOCITY
  use metricterms_new
#else
  use metricterms_n
#endif

  implicit none

  vectype, dimension(79, OP_NUM), intent(in) :: trial, lin 
  vectype, dimension(3), intent(out) :: ssterm, ddterm, rrterm, qqterm
  vectype :: temp
  real :: thimpb

  if(imp_mod.eq.1) then
     thimpb = 1.
  else
     thimpb = thimp
  endif

  ssterm = 0.
  ddterm = 0.
  rrterm = 0.
  qqterm = 0.

  if(numvar.lt.2) return          
          
  temp = b2psiv   (trial,lin,vzt79)
  ssterm(1) = ssterm(1) -     thimpb     *dt*temp
  ddterm(1) = ddterm(1) + (1.-thimpb*bdf)*dt*temp

  temp = b2psipsid(trial,lin,ps179,ni79)*dbf &
       + b2psipsid(trial,ps179,lin,ni79)*dbf
  ssterm(1) = ssterm(1) -     thimpb     *dt*temp
  ddterm(1) = ddterm(1) + (.5-thimpb*bdf)*dt*temp

  temp = b2b(trial,lin)
  ssterm(2) = ssterm(2) + temp
  ddterm(2) = ddterm(2) + temp*bdf

  temp = b2beta(trial,lin,eta79,hypi*sz79) &
       + b2bu  (trial,lin,pht79)
  ssterm(2) = ssterm(2) -     thimpb     *dt*temp
  ddterm(2) = ddterm(2) + (1.-thimpb*bdf)*dt*temp

  temp = b2bbd(trial,lin,bz179,ni79)*dbf &
       + b2bbd(trial,bz179,lin,ni79)*dbf
  ssterm(2) = ssterm(2) -     thimpb     *dt*temp
  ddterm(2) = ddterm(2) + (.5-thimpb*bdf)*dt*temp

  temp = b2bu  (trial,bz179,lin)
  rrterm(1) = rrterm(1) + thimpb*dt*temp
  qqterm(1) = qqterm(1) - thimpb*dt*temp*bdf

  temp = b2psiv(trial,ps179,lin)
  rrterm(2) = rrterm(2) + thimpb*dt*temp
  qqterm(2) = qqterm(2) - thimpb*dt*temp*bdf

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     temp = b2psipsid(trial,lin,ps079,ni79)*dbf &
          + b2psipsid(trial,ps079,lin,ni79)*dbf
     ssterm(1) = ssterm(1) -     thimpb     *dt*temp
     ddterm(1) = ddterm(1) + (1.-thimpb*bdf)*dt*temp

     temp = b2bbd (trial,lin,bz079,ni79)*dbf &
          + b2bbd (trial,bz079,lin,ni79)*dbf
     ssterm(2) = ssterm(2) -     thimpb     *dt*temp
     ddterm(2) = ddterm(2) + (1.-thimpb*bdf)*dt*temp

     temp = b2bu  (trial,bz079,lin)
     rrterm(1) = rrterm(1) +     thimpb     *dt*temp
     qqterm(1) = qqterm(1) + (1.-thimpb*bdf)*dt*temp
     
     temp = b2psiv(trial,ps079,lin)
     rrterm(2) = rrterm(2) +     thimpb     *dt*temp
     qqterm(2) = qqterm(2) + (1.-thimpb*bdf)*dt*temp
  endif

  ! NUMVAR = 3
  ! ~~~~~~~~~~
  if(numvar.ge.3) then

     temp = b2bchi(trial,lin,cht79)                  
     ssterm(2) = ssterm(2) -     thimpb     *dt*temp
     ddterm(2) = ddterm(2) + (1.-thimpb*bdf)*dt*temp

     temp = b2ped(trial,lin,ni79)*dbf*pefac
     ssterm(3) = ssterm(3) -     thimpb     *dt*temp
     ddterm(3) = ddterm(3) + (1.-thimpb*bdf)*dt*temp

     temp = b2bchi(trial,bz179,lin)                  
     rrterm(3) = rrterm(3) + thimpb*dt*temp
     qqterm(3) = qqterm(3) - thimpb*dt*temp*bdf
             
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        temp = b2bchi(trial,bz079,lin)
        rrterm(3) = rrterm(3) +     thimpb     *dt*temp
        qqterm(3) = qqterm(3) + (1.-thimpb*bdf)*dt*temp
     endif
  end if

end subroutine axial_field_lin


subroutine axial_field_nolin(trial, r4term)

  use basic
  use nintegrate_mod
#ifdef NEW_VELOCITY
  use metricterms_new
#else
  use metricterms_n
#endif

  implicit none

  vectype, intent(in), dimension(79, OP_NUM)  :: trial
  vectype, intent(out) :: r4term
  vectype :: temp
  
  r4term = 0.

  if(numvar .lt. 2) return

  ! density terms
  ! ~~~~~~~~~~~~~
  if(idens.eq.1 .and. (linear.eq.1 .or. eqsubtract.eq.1)) then
     r4term = r4term + dt* &
          (b2psipsid(trial,ps079,ps079,ni79)*dbf &
          +b2bbd    (trial,bz079,bz079,ni79)*dbf)

     if(numvar.ge.3) then
        r4term = r4term + dt* &
             (b2ped (trial,pe079,ni79)*dbf*pefac)
     endif
  endif

end subroutine axial_field_nolin


!======================================================================
! Electron Pressure Equation
!======================================================================
subroutine electron_pressure_lin(trial, lin, ssterm, ddterm, rrterm, qqterm)
  
  use basic
  use nintegrate_mod
#ifdef NEW_VELOCITY
  use metricterms_new
#else
  use metricterms_n
#endif

  implicit none

  vectype, dimension(79, OP_NUM), intent(in) :: trial, lin 
  vectype, dimension(3), intent(out) :: ssterm, ddterm, rrterm, qqterm
  vectype :: temp
  real :: thimpb

  if(imp_mod.eq.1) then
     thimpb = 1.
  else
     thimpb = thimp
  endif

  ssterm = 0.
  ddterm = 0.
  rrterm = 0.
  qqterm = 0.

  if(numvar.lt.3) return   

  temp = b3psipsieta(trial,lin,ps179,eta79) &
       + b3psipsieta(trial,ps179,lin,eta79)
  ssterm(1) = ssterm(1) -     thimp_ohm     *dt*temp
  ddterm(1) = ddterm(1) + (.5-thimp_ohm*bdf)*dt*temp

  temp = b3pebd(trial,pe179,lin,ni79)*dbf*pefac
  ssterm(2) = ssterm(2) - thimpb*dt*temp
  ddterm(2) = ddterm(2) - thimpb*dt*temp*bdf

  temp = b3bbeta(trial,lin,bz179,eta79) &
       + b3bbeta(trial,bz179,lin,eta79)
  ssterm(2) = ssterm(2) -     thimp_ohm     *dt*temp
  ddterm(2) = ddterm(2) + (.5-thimp_ohm*bdf)*dt*temp

  temp = b3pe(trial,lin)
  ssterm(3) = ssterm(3) + temp
  ddterm(3) = ddterm(3) + temp*bdf

  temp = b3pebd(trial,lin,bzt79,ni79)*dbf*pefac &
       + p1pu  (trial,lin,pht79)                & 
       + p1pchi(trial,lin,cht79)                & 
       + b3pedkappa(trial,lin,ni79,kap79,hypp*sz79)
  ssterm(3) = ssterm(3) -     thimpb     *dt*temp
  ddterm(3) = ddterm(3) + (1.-thimpb*bdf)*dt*temp

  temp = p1pu(trial,pe179,lin)
  rrterm(1) = rrterm(1) + thimpb*dt*temp
  qqterm(1) = qqterm(1) - thimpb*dt*temp*bdf

  temp = p1pchi(trial,pe179,lin)
  rrterm(3) = rrterm(3) + thimpb*dt*temp
  qqterm(3) = qqterm(3) - thimpb*dt*temp*bdf

  if(ipres.eq.0) then
     temp = p1uus  (trial,lin,ph179,sig79) &
          + p1uus  (trial,ph179,lin,sig79) &
          + p1uchis(trial,lin,ch179,sig79) 
     rrterm(1) = rrterm(1) +     thimpb     *dt*temp
     qqterm(1) = qqterm(1) + (.5-thimpb*bdf)*dt*temp

     temp = p1vvs  (trial,lin,vz179,sig79) &
          + p1vvs  (trial,vz179,lin,sig79)
     rrterm(2) = rrterm(2) +     thimpb     *dt*temp
     qqterm(2) = qqterm(2) + (.5-thimpb*bdf)*dt*temp

     temp = p1chichis(trial,lin,ch179,sig79) &
          + p1chichis(trial,ch179,lin,sig79) &
          + p1uchis  (trial,ph179,lin,sig79) 
     rrterm(3) = rrterm(3) +     thimpb     *dt*temp
     qqterm(3) = qqterm(3) + (.5-thimpb*bdf)*dt*temp
  endif

  ! Anisotropic Heat Flux
  if(kappar.ne.0.) then
     temp = (p1kappar(trial,lin,pst79,pet79,ni79,b2i79,kar79) &
          +  p1kappar(trial,pst79,lin,pet79,ni79,b2i79,kar79))
     ssterm(1) = ssterm(1) - thimpb*dt*temp
     ddterm(1) = ddterm(1) - thimpb*dt*temp*bdf
     
     temp = p1kappar(trial,pst79,pst79,lin,ni79,b2i79,kar79)
     ssterm(3) = ssterm(3) -     thimpb     *dt*temp
     ddterm(3) = ddterm(3) + (1.-thimpb*bdf)*dt*temp

     if(linear.eq.1 .or. eqsubtract.eq.1) then
        temp = (p1kappar(trial,lin,pss79,pe079,ni79,b2i79,kar79) &
             +  p1kappar(trial,pss79,lin,pe079,ni79,b2i79,kar79))
        ddterm(1) = ddterm(1) + dt*temp
     endif
  endif

  ! Cross-field Heat Flux
  if(kappax.ne.0.) then
     temp = p1kappax(trial,pe179,lin,ni79,kax79)
     ssterm(2) = ssterm(2) - thimpb*dt*temp
     ddterm(2) = ddterm(2) - thimpb*dt*temp*bdf

     temp = p1kappax(trial,lin,bzt79,ni79,kax79)
     ssterm(3) = ssterm(3) -     thimpb     *dt*temp
     ddterm(3) = ddterm(3) + (1.-thimpb*bdf)*dt*temp

     if(linear.eq.1 .or. eqsubtract.eq.1) then
        temp = p1kappax(trial,pe079,lin,ni79,kax79)
        ssterm(2) = ssterm(2) -     thimpb     *dt*temp
        ddterm(2) = ddterm(2) + (1.-thimpb*bdf)*dt*temp
     endif
  endif
              
  if(linear.eq.1 .or. eqsubtract.eq.1) then
     temp = b3psipsieta(trial,lin,ps079,eta79) &
          + b3psipsieta(trial,ps079,lin,eta79)
     ssterm(1) = ssterm(1) -     thimp_ohm     *dt*temp
     ddterm(1) = ddterm(1) + (1.-thimp_ohm*bdf)*dt*temp

     temp = b3bbeta(trial,lin,bz079,eta79) &
          + b3bbeta(trial,bz079,lin,eta79)
     ssterm(2) = ssterm(2) -     thimp_ohm     *dt*temp
     ddterm(2) = ddterm(2) + (1.-thimp_ohm*bdf)*dt*temp
     
     temp = b3pebd(trial,pe079,lin,ni79)*dbf*pefac
     ssterm(2) = ssterm(2) -     thimpb     *dt*temp
     ddterm(2) = ddterm(2) + (1.-thimpb*bdf)*dt*temp
     
     temp = p1pu(trial,pe079,lin)
     rrterm(1) = rrterm(1) +     thimpb     *dt*temp
     qqterm(1) = qqterm(1) + (1.-thimpb*bdf)*dt*temp
     
     temp = p1pchi(trial,pe079,lin)                
     rrterm(3) = rrterm(3) +     thimpb     *dt*temp
     qqterm(3) = qqterm(3) + (1.-thimpb*bdf)*dt*temp
     
     if(ipres.eq.0) then
        temp = p1uus  (trial,lin,ph079,sig79) &
             + p1uus  (trial,ph079,lin,sig79) &
             + p1uchis(trial,lin,ch079,sig79) 
        rrterm(1) = rrterm(1) +     thimpb     *dt*temp
        qqterm(1) = qqterm(1) + (1.-thimpb*bdf)*dt*temp
        
        temp = p1vvs  (trial,lin,vz079,sig79) &
             + p1vvs  (trial,vz079,lin,sig79)
        rrterm(2) = rrterm(2) +     thimpb     *dt*temp
        qqterm(2) = qqterm(2) + (1.-thimpb*bdf)*dt*temp
        
        temp = p1chichis(trial,lin,ch079,sig79) &
             + p1chichis(trial,ch079,lin,sig79) &
             + p1uchis  (trial,ph079,lin,sig79) 
        rrterm(3) = rrterm(3) +     thimpb     *dt*temp
        qqterm(3) = qqterm(3) + (1.-thimpb*bdf)*dt*temp
     endif

  endif

end subroutine electron_pressure_lin

subroutine electron_pressure_nolin(trial, r4term)

  use basic
  use nintegrate_mod
#ifdef NEW_VELOCITY
  use metricterms_new
#else
  use metricterms_n
#endif

  implicit none

  vectype, intent(in), dimension(79, OP_NUM)  :: trial
  vectype, intent(out) :: r4term
  real :: temp
  
  r4term = 0.

  if(numvar .lt. 3) return

  ! source terms
  ! ~~~~~~~~~~~~
  ! ohmic heating
  r4term = r4term + dt*(gam-1.)* &
       (qpsipsieta(trial) &
       +qbbeta    (trial))

  ! viscous heating
  if(ipres.eq.0) then
     r4term = r4term - dt*(gam-1.)* &
          (quumu    (trial,pht79,pht79,vis79,      hypc*sz79) &
          +qvvmu    (trial,vzt79,vzt79,vis79,      hypv*sz79) &
          +quchimu  (trial,pht79,cht79,vis79,vic79,hypc*sz79) &
          +qchichimu(trial,cht79,cht79,      vic79,hypc*sz79) &
          +p1vip    (trial))
  endif

  ! density terms
  ! ~~~~~~~~~~~~~
  if(idens.eq.1 .and. (linear.eq.1 .or. eqsubtract.eq.1)) then
     r4term = r4term + dt* &
          (b3pebd(trial,pe079,bz079,ni79)*dbf*pefac &
          +p1kappar(trial,ps079,ps079,pe079,ni79,b2i79,kar79))
     if(ipres.eq.1) then
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
subroutine ludefall

  use p_data
  use t_data
  use basic
  use arrays
  use sparse
  use nintegrate_mod
  use diagnostics

  implicit none

  include 'mpif.h'


  integer :: itri, numelms, i
  integer :: def_fields
  real :: x, z, xmin, zmin, factor

  real :: tstart, tend, tfield, telm, tsizefield, tfinalize

  double precision cogcoords(3)

  tfield = 0.
  telm = 0.
  tsizefield = 0.

  call numfac(numelms)
  call getmincoord(xmin,zmin)

  if(isplitstep.eq.1) then
     ! form field matrices
     call zerosuperluarray(s1matrix_sm,vecsize)
     call zeromultiplyarray(d1matrix_sm,vecsize)
     call zeromultiplyarray(q1matrix_sm,vecsize)
     call zerosuperluarray(s2matrix_sm,vecsize)
     call zeromultiplyarray(d2matrix_sm,vecsize)
     call zeromultiplyarray(r2matrix_sm,vecsize)
     call zeromultiplyarray(q2matrix_sm,vecsize)
     call zeromultiplyarray(r14matrix_sm,vecsize)
     if(iresolve.eq.1) then
        call zerosuperluarray(s10matrix_sm,vecsize)
        call zeromultiplyarray(d10matrix_sm,vecsize)
        call zeromultiplyarray(r10matrix_sm,vecsize)
        call zeromultiplyarray(q10matrix_sm,vecsize)
     endif
     if(idens.eq.1) then
        call zerosuperluarray(s8matrix_sm,vecsize1)
        call zeromultiplyarray(d8matrix_sm,vecsize1)
        call zeromultiplyarray(q8matrix_sm,vecsize)
        call zeromultiplyarray(r8matrix_sm,vecsize)
     endif
     if(ipres.eq.1) then
        call zerosuperluarray(s9matrix_sm,vecsize1)
        call zeromultiplyarray(d9matrix_sm,vecsize1)
        call zeromultiplyarray(q9matrix_sm,vecsize)
        call zeromultiplyarray(r9matrix_sm,vecsize)
     endif
     r4 = 0.
     if(idens.eq.1) qn4 = 0.
     if(ipres.eq.1) qp4 = 0.
  else
     call zerosuperluarray(s1matrix_sm,vecsize)
     call zeromultiplyarray(d1matrix_sm,vecsize)
  endif

  q4 = 0.
    
  ! Determine which fields need to be calculated
  def_fields = FIELD_PSI + FIELD_I + FIELD_PHI + FIELD_ETA + FIELD_MU
  if(numvar.ge.2) def_fields = def_fields + FIELD_V
  if(numvar.ge.3) def_fields = def_fields + &
       FIELD_CHI + FIELD_PE + FIELD_B2I + FIELD_J + FIELD_P + FIELD_KAP
  if(idens.eq.1) then
     if(ipellet.eq.1 .or. ionization.eq.1) def_fields = def_fields + FIELD_SIG
     def_fields = def_fields + FIELD_N + FIELD_NI
  endif

  if(gyro.eq.1) then
     if(numvar.lt.3) def_fields = def_fields + FIELD_P + FIELD_PE + FIELD_B2I
  endif

  if(hyperc.ne.0 .and. numvar.ge.3) then
     def_fields = def_fields + FIELD_VOR + FIELD_COM
  endif

  if(isources.eq.1) def_fields = def_fields + FIELD_SRC


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
     call define_fields_79(itri, def_fields)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        tfield = tfield + tend - tstart
     endif
     
     ! add element's contribution to matrices
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call ludefvel_n(itri)
     call ludefphi_n(itri)
     if(idens.eq.1) call ludefden_n(itri)
     if(ipres.eq.1) call ludefpres_n(itri)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        telm = telm + tend - tstart
     endif

  end do

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

  ! since a proc is contributing values to parts of the vector
  ! it does not own, we call sumshareddofs so that these values
  ! get summed up for all values shared by multiple procs
  ! and then update these values
  call sumshareddofs(q4)
  call finalizearray(d1matrix_sm)

  if(isplitstep.eq.1) then
     call sumshareddofs(r4)
     if(idens.eq.1) call sumshareddofs(qn4)

     ! Finalize matrices for multiplication
     call finalizearray(q1matrix_sm)
     
     call finalizearray(d2matrix_sm)
     call finalizearray(r2matrix_sm)
     call finalizearray(q2matrix_sm)
     call finalizearray(r14matrix_sm)
     
     if(iresolve.eq.1) then
        call finalizearray(d10matrix_sm)
        call finalizearray(q10matrix_sm)
        call finalizearray(r10matrix_sm)
     endif

     if(idens.eq.1) then
        call finalizearray(d8matrix_sm)
        call finalizearray(q8matrix_sm)
        call finalizearray(r8matrix_sm)
     endif ! on idens.eq.1
     
     
     if(ipres.eq.1) then
        call finalizearray(d9matrix_sm)
        call finalizearray(q9matrix_sm)
        call finalizearray(r9matrix_sm)
     endif ! on ipres.eq.1
  endif

  if(myrank.eq.0 .and. itimer.eq.1) then
     call second(tend)
     tfinalize = tfinalize + tend - tstart
  endif

  if(myrank.eq.0 .and. iprint.ge.1) print *, "Done defining matrices."

  if(myrank.eq.0 .and. itimer.eq.1) then
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
  use nintegrate_mod
  use arrays
  use sparse

  implicit none

  integer, intent(in) :: itri

  integer :: i, i1, j, j1
  vectype, dimension(3,3) :: ssterm, ddterm
  vectype, dimension(3,4) :: rrterm, qqterm
  vectype :: temp

  integer :: vv1, vv0, vb1, vb0, vn1, vn0
  vectype, pointer :: vsource(:)

  if(isplitstep.eq.1) then
     vv1 = s1matrix_sm
     vv0 = d1matrix_sm
     vb0 = q1matrix_sm
     vn0 = r14matrix_sm
     vsource => r4
  else
     vv1 = s1matrix_sm
     vv0 = d1matrix_sm
     vb1 = s1matrix_sm
     vb0 = d1matrix_sm
     vn1 = s1matrix_sm
     vn0 = d1matrix_sm
     vsource => q4
  endif

  do i=1,18

     i1 = isvaln(itri,i)

     do j=1,18
        j1 = isvaln(itri,j)

        call vorticity_lin(g79(:,:,i),g79(:,:,j), &
             ssterm(1,:),ddterm(1,:),rrterm(1,:),qqterm(1,:),isplitstep)
        if(numvar.ge.2) then
           call axial_vel_lin(g79(:,:,i),g79(:,:,j), &
                ssterm(2,:),ddterm(2,:),rrterm(2,:),qqterm(2,:),isplitstep)
        endif
        if(numvar.ge.3) then
           call compression_lin(g79(:,:,i),g79(:,:,j), &
                ssterm(3,:),ddterm(3,:),rrterm(3,:),qqterm(3,:),isplitstep)
        endif

        if(isplitstep.eq.0) rrterm = -rrterm
           
        call insertval(vv1,ssterm(1,1),i1+phi_off,j1+phi_off,1)
        call insertval(vv0,ddterm(1,1),i1+phi_off,j1+phi_off,1)
        call insertval(vb0,qqterm(1,1),i1+phi_off,j1+psi_off,1)
        if(idens.eq.1) call insertval(vn0,qqterm(1,4),i1+phi_off,j1+den_off,1)
        if(isplitstep.eq.0) then
           call insertval(vb1,rrterm(1,1),i1+phi_off,j1+psi_off,1)
           if(idens.eq.1) call insertval(vn1,rrterm(1,4),i1+phi_off,j1+den_off,1)
        endif
        if(numvar.ge.2) then
           call insertval(vv1,ssterm(1,2),i1+phi_off,j1+ vz_off,1)
           call insertval(vv1,ssterm(2,1),i1+ vz_off,j1+phi_off,1)
           call insertval(vv1,ssterm(2,2),i1+ vz_off,j1+ vz_off,1)
           call insertval(vv0,ddterm(1,2),i1+phi_off,j1+ vz_off,1)
           call insertval(vv0,ddterm(2,1),i1+ vz_off,j1+phi_off,1)
           call insertval(vv0,ddterm(2,2),i1+ vz_off,j1+ vz_off,1)
           call insertval(vb0,qqterm(1,2),i1+phi_off,j1+ bz_off,1)
           call insertval(vb0,qqterm(2,1),i1+ vz_off,j1+psi_off,1)
           call insertval(vb0,qqterm(2,2),i1+ vz_off,j1+ bz_off,1)
           if(idens.eq.1) call insertval(vn0,qqterm(2,4),i1+vz_off,j1+den_off,1)
           if(isplitstep.eq.0) then
              call insertval(vb1,rrterm(1,2),i1+phi_off,j1+ bz_off,1)
              call insertval(vb1,rrterm(2,1),i1+ vz_off,j1+psi_off,1)
              call insertval(vb1,rrterm(2,2),i1+ vz_off,j1+ bz_off,1)
              if(idens.eq.1) call insertval(vn1,rrterm(2,4),i1+vz_off,j1+den_off,1)
           end if
        endif
        if(numvar.ge.3) then
           call insertval(vv1,ssterm(1,3),i1+phi_off,j1+chi_off,1)
           call insertval(vv1,ssterm(2,3),i1+ vz_off,j1+chi_off,1)
           call insertval(vv1,ssterm(3,1),i1+chi_off,j1+phi_off,1)
           call insertval(vv1,ssterm(3,2),i1+chi_off,j1+ vz_off,1)
           call insertval(vv1,ssterm(3,3),i1+chi_off,j1+chi_off,1)
           call insertval(vv0,ddterm(1,3),i1+phi_off,j1+chi_off,1)
           call insertval(vv0,ddterm(2,3),i1+ vz_off,j1+chi_off,1)
           call insertval(vv0,ddterm(3,1),i1+chi_off,j1+phi_off,1)
           call insertval(vv0,ddterm(3,2),i1+chi_off,j1+ vz_off,1)
           call insertval(vv0,ddterm(3,3),i1+chi_off,j1+chi_off,1)
           call insertval(vb0,qqterm(1,3),i1+phi_off,j1+ pe_off,1)
           call insertval(vb0,qqterm(2,3),i1+ vz_off,j1+ pe_off,1)
           call insertval(vb0,qqterm(3,1),i1+chi_off,j1+psi_off,1)
           call insertval(vb0,qqterm(3,2),i1+chi_off,j1+ bz_off,1)
           call insertval(vb0,qqterm(3,3),i1+chi_off,j1+ pe_off,1)
           if(idens.eq.1) call insertval(vn0,qqterm(3,4),i1+chi_off,j1+den_off,1)
           if(isplitstep.eq.0) then
              call insertval(vb1,rrterm(1,3),i1+phi_off,j1+ pe_off,1)
              call insertval(vb1,rrterm(2,3),i1+ vz_off,j1+ pe_off,1)
              call insertval(vb1,rrterm(3,1),i1+chi_off,j1+psi_off,1)
              call insertval(vb1,rrterm(3,2),i1+chi_off,j1+ bz_off,1)
              call insertval(vb1,rrterm(3,3),i1+chi_off,j1+ pe_off,1)
              if(idens.eq.1) call insertval(vn1,rrterm(3,4),i1+chi_off,j1+den_off,1)
           endif
        endif
        
        if(iresolve.eq.1) then
           call vorticity_lin(g79(:,:,i),g79(:,:,j), &
                ssterm(1,:),ddterm(1,:),rrterm(1,:),qqterm(1,:),0)
           call axial_vel_lin(g79(:,:,i),g79(:,:,j), &
                ssterm(2,:),ddterm(2,:),rrterm(2,:),qqterm(2,:),0)
           call compression_lin(g79(:,:,i),g79(:,:,j), &
                ssterm(3,:),ddterm(3,:),rrterm(3,:),qqterm(3,:),0)
           
           call insertval(s10matrix_sm,ssterm(1,1),i1,j1,1)
           call insertval(d10matrix_sm,ddterm(1,1),i1,j1,1)
           call insertval(r10matrix_sm,rrterm(1,1),i1,j1,1)
           call insertval(q10matrix_sm,qqterm(1,1),i1,j1,1)
           if(numvar.ge.2) then
              call insertval(s10matrix_sm,ssterm(1,2),i1  ,j1+6,1)
              call insertval(s10matrix_sm,ssterm(2,1),i1+6,j1  ,1)
              call insertval(s10matrix_sm,ssterm(2,2),i1+6,j1+6,1)
              call insertval(d10matrix_sm,ddterm(1,2),i1  ,j1+6,1)
              call insertval(d10matrix_sm,ddterm(2,1),i1+6,j1  ,1)
              call insertval(d10matrix_sm,ddterm(2,2),i1+6,j1+6,1)
              call insertval(r10matrix_sm,rrterm(1,2),i1  ,j1+6,1)
              call insertval(r10matrix_sm,rrterm(2,1),i1+6,j1  ,1)
              call insertval(r10matrix_sm,rrterm(2,2),i1+6,j1+6,1)
              call insertval(q10matrix_sm,qqterm(1,2),i1  ,j1+6,1)
              call insertval(q10matrix_sm,qqterm(2,1),i1+6,j1  ,1)
              call insertval(q10matrix_sm,qqterm(2,2),i1+6,j1+6,1)
           endif
           if(numvar.ge.3) then
              call insertval(s10matrix_sm,ssterm(1,3),i1,   j1+12,1)
              call insertval(s10matrix_sm,ssterm(2,3),i1+6, j1+12,1)
              call insertval(s10matrix_sm,ssterm(3,1),i1+12,j1,   1)
              call insertval(s10matrix_sm,ssterm(3,2),i1+12,j1+6, 1)
              call insertval(s10matrix_sm,ssterm(3,3),i1+12,j1+12,1)
              call insertval(d10matrix_sm,ddterm(1,3),i1,   j1+12,1)
              call insertval(d10matrix_sm,ddterm(2,3),i1+6, j1+12,1)
              call insertval(d10matrix_sm,ddterm(3,1),i1+12,j1,   1)
              call insertval(d10matrix_sm,ddterm(3,2),i1+12,j1+6, 1)
              call insertval(d10matrix_sm,ddterm(3,3),i1+12,j1+12,1)
              call insertval(r10matrix_sm,rrterm(1,3),i1,   j1+12,1)
              call insertval(r10matrix_sm,rrterm(2,3),i1+6, j1+12,1)
              call insertval(r10matrix_sm,rrterm(3,1),i1+12,j1,   1)
              call insertval(r10matrix_sm,rrterm(3,2),i1+12,j1+6, 1)
              call insertval(r10matrix_sm,rrterm(3,3),i1+12,j1+12,1)
              call insertval(q10matrix_sm,qqterm(1,3),i1,   j1+12,1)
              call insertval(q10matrix_sm,qqterm(2,3),i1+6, j1+12,1)
              call insertval(q10matrix_sm,qqterm(3,1),i1+12,j1,   1)
              call insertval(q10matrix_sm,qqterm(3,2),i1+12,j1+6, 1)
              call insertval(q10matrix_sm,qqterm(3,3),i1+12,j1+12,1)
           endif
           
        endif
     enddo               ! on j

     ! Definition of R4
     ! ================
     call vorticity_nolin(g79(:,:,i),temp)
     vsource(i1+phi_off) = vsource(i1+phi_off) + temp
     if(numvar.ge.2) then
        call axial_vel_nolin(g79(:,:,i),temp)
        vsource(i1+ vz_off) = vsource(i1+ vz_off) + temp
     endif
     if(numvar.ge.3) then
        call compression_nolin(g79(:,:,i),temp)
        vsource(i1+chi_off) = vsource(i1+chi_off) + temp
     endif
  enddo                  ! on i

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
  use nintegrate_mod
  use arrays
  use sparse

  implicit none

  integer, intent(in) :: itri

  integer :: i, i1, j, j1
  vectype, dimension(3,3) :: ssterm, ddterm, rrterm, qqterm
  vectype :: temp

  integer :: bb1, bb0, bv1, bv0
  vectype, pointer :: bsource(:)

  if(isplitstep.eq.1) then
     bb1 = s2matrix_sm
     bb0 = d2matrix_sm
     bv1 = r2matrix_sm
     bv0 = q2matrix_sm
     bsource => q4
  else
     bb1 = s1matrix_sm
     bb0 = d1matrix_sm
     bv1 = s1matrix_sm
     bv0 = d1matrix_sm
     bsource => q4
  endif

  do i=1,18

     i1 = isvaln(itri,i)

     do j=1,18

        j1 = isvaln(itri,j)

        call flux_lin(g79(:,:,i),g79(:,:,j), &
             ssterm(1,:),ddterm(1,:),rrterm(1,:),qqterm(1,:))
        if(numvar.ge.2) then
           call axial_field_lin(g79(:,:,i),g79(:,:,j), &
                ssterm(2,:),ddterm(2,:),rrterm(2,:),qqterm(2,:))
        endif
        if(numvar.ge.3) then
           call electron_pressure_lin(g79(:,:,i),g79(:,:,j), &
                ssterm(3,:),ddterm(3,:),rrterm(3,:),qqterm(3,:))
        endif

        if(isplitstep.eq.0) rrterm = -rrterm
        
        call insertval(bb1,ssterm(1,1),i1+psi_off,j1+psi_off,1)
        call insertval(bb0,ddterm(1,1),i1+psi_off,j1+psi_off,1)
        call insertval(bv1,rrterm(1,1),i1+psi_off,j1+phi_off,1)
        call insertval(bv0,qqterm(1,1),i1+psi_off,j1+phi_off,1)
        if(numvar.ge.2) then
           call insertval(bb1,ssterm(1,2),i1+psi_off,j1+ bz_off,1)
           call insertval(bb1,ssterm(2,1),i1+ bz_off,j1+psi_off,1)
           call insertval(bb1,ssterm(2,2),i1+ bz_off,j1+ bz_off,1)
           call insertval(bb0,ddterm(1,2),i1+psi_off,j1+ bz_off,1)
           call insertval(bb0,ddterm(2,1),i1+ bz_off,j1+psi_off,1)
           call insertval(bb0,ddterm(2,2),i1+ bz_off,j1+ bz_off,1)
           call insertval(bv1,rrterm(1,2),i1+psi_off,j1+ vz_off,1)
           call insertval(bv1,rrterm(2,1),i1+ bz_off,j1+phi_off,1)
           call insertval(bv1,rrterm(2,2),i1+ bz_off,j1+ vz_off,1)
           call insertval(bv0,qqterm(1,2),i1+psi_off,j1+ vz_off,1)
           call insertval(bv0,qqterm(2,1),i1+ bz_off,j1+phi_off,1)
           call insertval(bv0,qqterm(2,2),i1+ bz_off,j1+ vz_off,1)
        endif
        if(numvar .eq. 3) then        
           call insertval(bb1,ssterm(1,3),i1+psi_off,j1+ pe_off,1)
           call insertval(bb1,ssterm(2,3),i1+ bz_off,j1+ pe_off,1)
           call insertval(bb1,ssterm(3,3),i1+ pe_off,j1+ pe_off,1)
           call insertval(bb1,ssterm(3,1),i1+ pe_off,j1+psi_off,1)
           call insertval(bb1,ssterm(3,2),i1+ pe_off,j1+ bz_off,1)
           call insertval(bb0,ddterm(1,3),i1+psi_off,j1+ pe_off,1)
           call insertval(bb0,ddterm(2,3),i1+ bz_off,j1+ pe_off,1)
           call insertval(bb0,ddterm(3,3),i1+ pe_off,j1+ pe_off,1)
           call insertval(bb0,ddterm(3,1),i1+ pe_off,j1+psi_off,1)
           call insertval(bb0,ddterm(3,2),i1+ pe_off,j1+ bz_off,1)
           call insertval(bv1,rrterm(1,3),i1+psi_off,j1+chi_off,1)
           call insertval(bv1,rrterm(2,3),i1+ bz_off,j1+chi_off,1)
           call insertval(bv1,rrterm(3,3),i1+ pe_off,j1+chi_off,1)
           call insertval(bv1,rrterm(3,1),i1+ pe_off,j1+phi_off,1)
           call insertval(bv1,rrterm(3,2),i1+ pe_off,j1+ vz_off,1)
           call insertval(bv0,qqterm(1,3),i1+psi_off,j1+chi_off,1)
           call insertval(bv0,qqterm(2,3),i1+ bz_off,j1+chi_off,1)
           call insertval(bv0,qqterm(3,3),i1+ pe_off,j1+chi_off,1)
           call insertval(bv0,qqterm(3,1),i1+ pe_off,j1+phi_off,1)
           call insertval(bv0,qqterm(3,2),i1+ pe_off,j1+ vz_off,1)
        endif
       
     enddo ! on j


     call flux_nolin(g79(:,:,i),temp)
     bsource(i1+psi_off) = bsource(i1+psi_off) + temp
     if(numvar.ge.2) then
        call axial_field_nolin(g79(:,:,i),temp)
        bsource(i1+ bz_off) = bsource(i1+ bz_off) + temp
     endif
     if(numvar.ge.3) then
        call electron_pressure_nolin(g79(:,:,i),temp)
        bsource(i1+ pe_off) = bsource(i1+ pe_off) + temp
     endif

  enddo ! on i

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

#ifdef NEW_VELOCITY
  use metricterms_new
#else
  use metricterms_n
#endif


  implicit none

  integer, intent(in) :: itri

  integer :: i, i1, ione, j, j1, jone
  vectype :: ssterm, ddterm
  vectype, dimension(3) :: rrterm, qqterm
  vectype :: temp

  integer :: nn1, nn0, nv1, nv0
  vectype, pointer :: nsource(:)

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
        
        temp = n1ndenm(g79(:,:,i),g79(:,:,j),denm,hypp*sz79) &
             + n1nu   (g79(:,:,i),g79(:,:,j),pht79)
        ssterm = ssterm -     thimp     *dt*temp
        ddterm = ddterm + (1.-thimp*bdf)*dt*temp

        temp = n1nu(g79(:,:,i),n179,g79(:,:,j))
        rrterm(1) = rrterm(1) + thimp*dt*temp
        qqterm(1) = qqterm(1) - thimp*dt*temp*bdf

        if(linear.eq.1 .or. eqsubtract.eq.1) then
           temp = n1nu  (g79(:,:,i),n079,g79(:,:,j))
           rrterm(1) = rrterm(1) +     thimp     *dt*temp
           qqterm(1) = qqterm(1) + (1.-thimp*bdf)*dt*temp
        endif

        ! NUMVAR = 3
        ! ~~~~~~~~~~
        if(numvar.ge.3) then
           temp = n1nchi(g79(:,:,i),g79(:,:,j),cht79)
           ssterm = ssterm -     thimp     *dt*temp
           ddterm = ddterm + (1.-thimp*bdf)*dt*temp
           
           temp = n1nchi(g79(:,:,i),n179,g79(:,:,j))
           rrterm(3) = rrterm(3) + thimp*dt*temp
           qqterm(3) = qqterm(3) - thimp*dt*temp*bdf

           if(linear.eq.1 .or. eqsubtract.eq.1) then
              temp = n1nchi(g79(:,:,i),n079,g79(:,:,j))
              rrterm(3) = rrterm(3) +     thimp     *dt*temp
              qqterm(3) = qqterm(3) + (1.-thimp*bdf)*dt*temp
           endif
        endif

        if(isplitstep.eq.0) rrterm = -rrterm

        call insertval(nn1, ssterm, ione+den_off, jone+den_off, 1)
        call insertval(nn0, ddterm, ione+den_off, jone+den_off, 1)
        call insertval(nv1, rrterm(1), i1+den_off,j1+phi_off, 1)
        call insertval(nv0, qqterm(1), i1+den_off,j1+phi_off, 1)
        if(numvar.ge.2) then
           call insertval(nv1,rrterm(2), i1+den_off,j1+vz_off,1)
           call insertval(nv0,qqterm(2), i1+den_off,j1+vz_off,1)
        endif
        if(numvar.ge.3) then
           call insertval(nv1,rrterm(3), i1+den_off,j1+chi_off,1)
           call insertval(nv0,qqterm(3), i1+den_off,j1+chi_off,1)
        endif
     enddo                     ! on j

     nsource(ione+den_off) = nsource(ione+den_off) + dt* &
             n1s(g79(:,:,i),sig79)

     if(linear.eq.1 .or. eqsubtract.eq.1) then
        nsource(ione+den_off) = nsource(ione+den_off) + dt* &
             (n1ndenm(g79(:,:,i),n079,denm,hypp*sz79))
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

#ifdef NEW_VELOCITY
  use metricterms_new
#else
  use metricterms_n
#endif


  implicit none

  integer, intent(in) :: itri

  integer :: i, i1, ione, j, j1, jone
  vectype :: ssterm, ddterm
  vectype, dimension(3) :: rrterm, qqterm
  vectype :: temp

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
        temp = int2(g79(:,:,i),g79(:,:,j),weight_79,79)
        ssterm = ssterm + temp
        ddterm = ddterm + temp*bdf
        
        temp = p1pu   (g79(:,:,i),g79(:,:,j),pht79)
        ssterm = ssterm -     thimp     *dt*temp
        ddterm = ddterm + (1.-thimp*bdf)*dt*temp

        temp = p1pu(g79(:,:,i),p179,g79(:,:,j))
        rrterm(1) = rrterm(1) + thimp*dt*temp
        qqterm(1) = qqterm(1) - thimp*dt*temp*bdf

        if(linear.eq.1 .or. eqsubtract.eq.1) then
           temp = p1pu  (g79(:,:,i),p079,g79(:,:,j))
           rrterm(1) = rrterm(1) +     thimp     *dt*temp
           qqterm(1) = qqterm(1) + (1.-thimp*bdf)*dt*temp
        endif

        ! NUMVAR = 3
        ! ~~~~~~~~~~
        if(numvar.ge.3) then
           temp = p1pchi    (g79(:,:,i),g79(:,:,j),cht79) &
                + b3pedkappa(g79(:,:,i),g79(:,:,j),ni79,kap79,hypp*sz79) &
                + p1kappar  (g79(:,:,i),pst79,pst79,g79(:,:,j),ni79,b2i79,kar79) &
                + p1kappax  (g79(:,:,i),g79(:,:,j),bzt79,ni79,kax79)
           ssterm = ssterm -     thimp     *dt*temp
           ddterm = ddterm + (1.-thimp*bdf)*dt*temp
           
           temp = p1pchi(g79(:,:,i),p179,g79(:,:,j))
           rrterm(3) = rrterm(3) + thimp*dt*temp
           qqterm(3) = qqterm(3) - thimp*dt*temp*bdf

           temp = p1uus  (g79(:,:,i),g79(:,:,j),ph179,sig79) &
                + p1uus  (g79(:,:,i),ph179,g79(:,:,j),sig79) &
                + p1uchis(g79(:,:,i),g79(:,:,j),ch179,sig79) 
           rrterm(1) = rrterm(1) +     thimp     *dt*temp
           qqterm(1) = qqterm(1) + (.5-thimp*bdf)*dt*temp

           temp = p1vvs  (g79(:,:,i),g79(:,:,j),vz179,sig79) &
                + p1vvs  (g79(:,:,i),vz179,g79(:,:,j),sig79)
           rrterm(2) = rrterm(2) +     thimp     *dt*temp
           qqterm(2) = qqterm(2) + (.5-thimp*bdf)*dt*temp

           temp = p1chichis(g79(:,:,i),g79(:,:,j),ph179,sig79) &
                + p1chichis(g79(:,:,i),ph179,g79(:,:,j),sig79) &
                + p1uchis  (g79(:,:,i),ph079,g79(:,:,j),sig79) 
           rrterm(3) = rrterm(3) +     thimp     *dt*temp
           qqterm(3) = qqterm(3) + (.5-thimp*bdf)*dt*temp

           if(linear.eq.1 .or. eqsubtract.eq.1) then
              temp = p1pchi(g79(:,:,i),p079,g79(:,:,j))
              rrterm(3) = rrterm(3) +     thimp     *dt*temp
              qqterm(3) = qqterm(3) + (1.-thimp*bdf)*dt*temp

              temp = p1uus  (g79(:,:,i),g79(:,:,j),ph079,sig79) &
                   + p1uus  (g79(:,:,i),ph079,g79(:,:,j),sig79) &
                   + p1uchis(g79(:,:,i),g79(:,:,j),ch079,sig79) 
              rrterm(1) = rrterm(1) +     thimp     *dt*temp
              qqterm(1) = qqterm(1) + (1.-thimp*bdf)*dt*temp

              temp = p1vvs  (g79(:,:,i),g79(:,:,j),vz079,sig79) &
                   + p1vvs  (g79(:,:,i),vz079,g79(:,:,j),sig79)
              rrterm(2) = rrterm(2) +     thimp     *dt*temp
              qqterm(2) = qqterm(2) + (1.-thimp*bdf)*dt*temp

              temp = p1chichis(g79(:,:,i),g79(:,:,j),ph079,sig79) &
                   + p1chichis(g79(:,:,i),ph079,g79(:,:,j),sig79) &
                   + p1uchis  (g79(:,:,i),ph079,g79(:,:,j),sig79) 
              rrterm(3) = rrterm(3) +     thimp     *dt*temp
              qqterm(3) = qqterm(3) + (1.-thimp*bdf)*dt*temp
           endif
        endif

        if(isplitstep.eq.0) rrterm = -rrterm
        
        call insertval(s9matrix_sm, ssterm, ione, jone, 1)
        call insertval(d9matrix_sm, ddterm, ione, jone, 1)
        call insertval(r9matrix_sm, rrterm(1), i1, j1, 1)
        call insertval(q9matrix_sm, qqterm(1), i1, j1, 1)
        if(numvar.ge.2) then
           call insertval(r9matrix_sm,rrterm(2), i1,j1+6,1)
           call insertval(q9matrix_sm,qqterm(2), i1,j1+6,1)
        endif
        if(numvar.ge.3) then
           call insertval(r9matrix_sm,rrterm(3), i1,j1+12,1)
           call insertval(q9matrix_sm,qqterm(3), i1,j1+12,1)
        endif

     enddo                     ! on j

     qp4(ione) = qp4(ione) + dt* &
          (b3psipsieta(g79(:,:,i),pst79,pst79,eta79))

     if(numvar.ge.2) then
        qp4(ione) = qp4(ione) + dt* &
             (b3bbeta(g79(:,:,i),bzt79,bzt79,eta79))
     endif

     if(numvar.ge.3) then
        qp4(ione) = qp4(ione) + dt* &
             (b3pebd(g79(:,:,i),pet79,bzt79,ni79))
     endif

     if(linear.eq.1 .or. eqsubtract.eq.1) then
                
        qp4(ione) = qp4(ione) + dt* &
             (b3pedkappa(g79(:,:,i),p079,ni79,kap79,hypp*sz79) &
             +p1kappar  (g79(:,:,i),ps079,ps079,p079,ni79,b2i79,kar79) &
             +p1kappax  (g79(:,:,i),pe079,bz079,ni79,kax79) &
             +p1uus    (g79(:,:,i),ph079,ph079,sig79) &
             +p1vvs    (g79(:,:,i),vz079,vz079,sig79) &
             +p1chichis(g79(:,:,i),ch079,ch079,sig79) &
             +p1uchis  (g79(:,:,i),ph079,ch079,sig79))
     endif

     ! EQUILIBRIUM TERMS
     if(linear.eq.0 .and. eqsubtract.eq.1) then 

        qp4(ione) = qp4(ione) + dt* &
             (p1pu   (g79(:,:,i),p079,ph079))
       
        if(numvar.ge.3) then
           qp4(ione) = qp4(ione) + dt* &
                (p1pchi(g79(:,:,i),p079,ch079))
        endif
     endif

  enddo                     ! on i
end subroutine ludefpres_n
