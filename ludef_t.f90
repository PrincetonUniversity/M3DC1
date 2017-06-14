!======================================================================
! Vorticity Equation
!======================================================================
subroutine vorticity_lin(trialx, lin, ssterm, ddterm, r_bf, q_bf, advfield, &
     izone)
  
  use basic
  use arrays
  use m3dc1_nint
  use metricterms_new
  use two_fluid

  implicit none

  vectype, dimension(MAX_PTS, OP_NUM, dofs_per_element), intent(in) :: trialx
  vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: lin 
  vectype, dimension(dofs_per_element,num_fields), intent(out) :: ssterm, ddterm
  vectype, dimension(dofs_per_element), intent(out) :: r_bf, q_bf
  integer, intent(in) :: advfield   ! if advfield = 1, eliminate rrterm by
                                    ! using analytic form of advanced field
  integer, intent(in) :: izone
  vectype :: temp
  vectype, dimension(dofs_per_element) :: tempx
  real :: ththm, nv, thimp_bf
  integer :: i

  vectype :: freq_fac

  vectype, dimension(MAX_PTS, OP_NUM) :: trial

  if(numvar.eq.1) then
     nv = 1.
  else
     nv = .5
  end if

  if(imp_bf.eq.1 .and. isplitstep.eq.0) then
     thimp_bf = thimp
  else
     thimp_bf = 0.
  end if

  select case(imp_mod)
  case(0)
     ththm = (1.-thimp*bdf)*thimp
  case(1)
     ththm = -bdf*thimp**2
  case(2)
     ththm = -bdf*thimp
  case default
     ththm = (1.-thimp*bdf)*thimp
  end select

  if(itime_independent.eq.1) then
#ifdef USECOMPLEX
     freq_fac = (0,1)*frequency*dt
#else
     freq_fac = 0.
#endif
  else
     freq_fac = 1.
  end if

  ssterm = 0.
  ddterm = 0.
  q_bf = 0.
  r_bf = 0.

  if(istatic.eq.1) then
     if(.not.surface_int) then
        tempx = intx2(trialx(:,OP_1,:),lin(:,OP_1))
        ssterm(:,u_g) = tempx
        ddterm(:,u_g) = tempx
     endif
     return
  endif

  if(izone.ne.1) then
     if(inonormalflow.eq.1 .and. .not.surface_int) then
        tempx = intx2(trialx(:,OP_1,:),lin(:,OP_1))
     else
        tempx = v1un(trialx,lin,nt79)
     end if
     ssterm(:,u_g) = ssterm(:,u_g) + tempx
     if(numvar.ge.3) then
        if(inoslip_pol.eq.0) then
           tempx = v1chin(trialx,lin,nt79)*chiiner
           ssterm(:,chi_g) = ssterm(:,chi_g) + tempx
        end if
     end if
     return
  end if

  ! Regularization term
  ! ~~~~~~~~~~~~~~~~~~~
  if(inonormalflow.eq.0 .and. (.not.surface_int)) then
     tempx = -regular*intx2(trialx(:,OP_1,:),lin(:,OP_1))
     ssterm(:,u_g) = ssterm(:,u_g) + tempx
     ddterm(:,u_g) = ddterm(:,u_g) + tempx*bdf
  end if


  ! Time Derivative
  ! ~~~~~~~~~~~~~~~
  tempx = v1un(trialx,lin,nt79)*freq_fac
  ssterm(:,u_g) = ssterm(:,u_g) + tempx
  if(itime_independent.eq.0) ddterm(:,u_g) = ddterm(:,u_g) + tempx*bdf
  
  if(numvar.ge.3) then
     tempx = v1chin(trialx,lin,nt79)*chiiner*freq_fac
     ssterm(:,chi_g) = ssterm(:,chi_g) + tempx
     if(itime_independent.eq.0) ddterm(:,chi_g) = ddterm(:,chi_g) + tempx*bdf
  end if


  ! Viscosity
  ! ~~~~~~~~~
  tempx = v1umu(trialx,lin,vis79,vic79) &
       +  v1us (trialx,lin,sig79)
  ssterm(:,u_g) = ssterm(:,u_g) -     thimp     *dt*tempx
  ddterm(:,u_g) = ddterm(:,u_g) + (1.-thimp*bdf)*dt*tempx
  if(numvar.ge.2) then
     tempx = v1vmu(trialx,lin,vis79,vic79)
     ssterm(:,vz_g) = ssterm(:,vz_g) -     thimp     *dt*tempx
     ddterm(:,vz_g) = ddterm(:,vz_g) + (1.-thimp*bdf)*dt*tempx
  end if
  if(numvar.ge.3) then
     tempx = v1chimu(trialx,lin,vis79,vic79) &
          +  v1chis (trialx,lin,sig79)
     ssterm(:,chi_g) = ssterm(:,chi_g) -     thimp     *dt*tempx
     ddterm(:,chi_g) = ddterm(:,chi_g) + (1.-thimp*bdf)*dt*tempx
  end if


  ! Advection
  ! ~~~~~~~~~
  if(linear.eq.0) then 
     tempx = v1uun(trialx,lin,ph179,nt79) &
          +  v1uun(trialx,ph179,lin,nt79)
     ssterm(:,u_g) = ssterm(:,u_g) -     thimp     *dt*tempx
     ddterm(:,u_g) = ddterm(:,u_g) + (.5-thimp*bdf)*dt*tempx

     tempx = v1uvn(trialx,lin,vz179,nt79)
     ssterm(:,u_g) = ssterm(:,u_g) -     thimp     *dt*tempx
     ddterm(:,u_g) = ddterm(:,u_g) + (.5-thimp*bdf)*dt*tempx

     tempx = v1uchin(trialx,lin,ch179,nt79)       
     ssterm(:,u_g) = ssterm(:,u_g) -     thimp     *dt*tempx
     ddterm(:,u_g) = ddterm(:,u_g) + (.5-thimp*bdf)*dt*tempx

     if(numvar.ge.2) then
        tempx = v1vvn(trialx,lin,vz179,nt79) &
             +  v1vvn(trialx,vz179,lin,nt79) &
             +  v1uvn(trialx,ph179,lin,nt79)
        ssterm(:,vz_g) = ssterm(:,vz_g) -     thimp     *dt*tempx
        ddterm(:,vz_g) = ddterm(:,vz_g) + (.5-thimp*bdf)*dt*tempx

        tempx = v1vchin(trialx,lin,ch179,nt79)       
        ssterm(:,vz_g) = ssterm(:,vz_g) -     thimp     *dt*tempx
        ddterm(:,vz_g) = ddterm(:,vz_g) + (.5-thimp*bdf)*dt*tempx
     end if

     if(numvar.ge.3) then
        tempx = v1uchin  (trialx,ph179,lin,nt79) &
             +  v1vchin  (trialx,vz179,lin,nt79) &
             +  v1chichin(trialx,lin,ch179,nt79) &
             +  v1chichin(trialx,ch179,lin,nt79)
        ssterm(:,chi_g) = ssterm(:,chi_g) -     thimp     *dt*tempx
        ddterm(:,chi_g) = ddterm(:,chi_g) + (.5-thimp*bdf)*dt*tempx
      end if
  end if

  if(eqsubtract.eq.1) then
     tempx = v1uun  (trialx,lin,ph079,nt79) &    
          +  v1uun  (trialx,ph079,lin,nt79) &
          +  v1uvn  (trialx,lin,vz079,nt79) &
          +  v1uchin(trialx,lin,ch079,nt79)
     ssterm(:,u_g) = ssterm(:,u_g) -     thimp     *dt*tempx
     ddterm(:,u_g) = ddterm(:,u_g) + (1.-thimp*bdf)*dt*tempx

     if(numvar.ge.2) then
        tempx = v1uvn  (trialx,ph079,lin,nt79) & 
             +  v1vvn  (trialx,lin,vz079,nt79) &
             +  v1vvn  (trialx,vz079,lin,nt79) &
             +  v1vchin(trialx,lin,ch079,nt79)
        ssterm(:,vz_g) = ssterm(:,vz_g) -     thimp     *dt*tempx
        ddterm(:,vz_g) = ddterm(:,vz_g) + (1.-thimp*bdf)*dt*tempx
     end if

     if(numvar.ge.3) then
        tempx = v1uchin  (trialx,ph079,lin,nt79) &
             +  v1vchin  (trialx,vz079,lin,nt79) &
             +  v1chichin(trialx,lin,ch079,nt79) &
             +  v1chichin(trialx,ch079,lin,nt79)
        ssterm(:,chi_g) = ssterm(:,chi_g) -     thimp     *dt*tempx
        ddterm(:,chi_g) = ddterm(:,chi_g) + (1.-thimp*bdf)*dt*tempx
     end if

     if(idens.eq.1) then
        ddterm(:,den_g) = ddterm(:,den_g) + dt* &
             (v1uun    (trialx,ph079,ph079,lin) &
             +v1uvn    (trialx,ph079,vz079,lin) &
             +v1vvn    (trialx,vz079,vz079,lin) &
             +v1uchin  (trialx,vz079,ch079,lin) &
             +v1vchin  (trialx,vz079,ch079,lin) &
             +v1chichin(trialx,ch079,ch079,lin))
     end if
  end if


  ! Grad(p)
  ! ~~~~~~~
  if(numvar.ge.3 .or. ipres.eq.1) then
     ! Split time-step
     if(advfield.eq.1) then
        ddterm(:,p_g) = ddterm(:,p_g) + dt*v1p(trialx,lin)

        ! "Parabolization" terms
        tempx = v1up(trialx,lin,pt79)
        ssterm(:,u_g) = ssterm(:,u_g) - thimp*thimp*dt*dt*tempx
        ddterm(:,u_g) = ddterm(:,u_g) +       ththm*dt*dt*tempx

        if(numvar.ge.2) then
           tempx = v1vp(trialx,lin,pt79)
           ssterm(:,vz_g) = ssterm(:,vz_g) - thimp*thimp*dt*dt*tempx
           ddterm(:,vz_g) = ddterm(:,vz_g) +       ththm*dt*dt*tempx
        endif

        if(numvar.ge.3) then
           tempx = v1chip(trialx,lin,pt79)
           ssterm(:,chi_g) = ssterm(:,chi_g) - thimp*thimp*dt*dt*tempx
           ddterm(:,chi_g) = ddterm(:,chi_g) +       ththm*dt*dt*tempx
        end if

     ! Unsplit time-step
     else
        if    (kinetic.le.1) then
        tempx = v1p(trialx,lin)
        ssterm(:,p_g) = ssterm(:,p_g) -     thimp     *dt*tempx
        ddterm(:,p_g) = ddterm(:,p_g) + (1.-thimp*bdf)*dt*tempx
 
        ! CGL (anisotropic pressure)
        elseif(kinetic.gt.1) then
           tempx = v1par(trialx,lin)    &
                + v1parb2ipsipsi(trialx,lin,b2i79,pstx79,pstx79)  &
                + v1parb2ipsib(trialx,lin,b2i79,pstx79,bztx79)
           ssterm(:,pe_g) = ssterm(:,pe_g) -     thimp     *dt*tempx
           ddterm(:,pe_g) = ddterm(:,pe_g) + (1.-thimp*bdf)*dt*tempx
           ssterm(:,p_g)  = ssterm(:,p_g)  +     thimp     *dt*tempx
           ddterm(:,p_g)  = ddterm(:,p_g)  - (1.-thimp*bdf)*dt*tempx
        endif
     end if
  end if

  ! Gravity
  ! ~~~~~~~
  if(idens.eq.1) then
     ! Split time-step
     if(advfield.eq.1) then
        ddterm(:,den_g) = ddterm(:,den_g) + dt* &
             v1ngrav(trialx,lin)

        ! parabolization terms
        tempx = v1ungrav (trialx,lin,nt79)
        ssterm(:,u_g) = ssterm(:,u_g) - thimp*thimp*dt*dt*tempx
        ddterm(:,u_g) = ddterm(:,u_g) +       ththm*dt*dt*tempx
     
        if(numvar.ge.3) then
           tempx = v1chingrav (trialx,lin,nt79)
           ssterm(:,chi_g) = ssterm(:,chi_g) - thimp*thimp*dt*dt*tempx
           ddterm(:,chi_g) = ddterm(:,chi_g) +       ththm*dt*dt*tempx
        end if
        
     ! Unsplit time-step
     else
        tempx = v1ngrav(trialx,lin)
        ssterm(:,den_g) = ssterm(:,den_g) -     thimp     *dt*tempx
        ddterm(:,den_g) = ddterm(:,den_g) + (1.-thimp*bdf)*dt*tempx
     end if
  end if


  ! JxB
  ! ~~~
  ! Split time-step
  if(advfield.eq.1) then
     ddterm(:,psi_g) = ddterm(:,psi_g) + dt* &
          (v1psipsi(trialx,lin,pss79) &
          +v1psipsi(trialx,pss79,lin) &
          +v1psib  (trialx,lin,bzs79))

     if(numvar.ge.2) then
        ddterm(:,bz_g) = ddterm(:,bz_g) + dt* &
             (v1psib(trialx,pss79,lin) &
             +v1bb  (trialx,lin,bzs79) &
             +v1bb  (trialx,bzs79,lin))
     end if

     if(use_external_fields .and. linear.eq.0) then
        ddterm(:,psi_g) = ddterm(:,psi_g) + dt* &
             (v1psipsi(trialx,psx79,lin) &
             +v1psib  (trialx,lin,bzx79))

        if(numvar.ge.2) then
           ddterm(:,bz_g) = ddterm(:,bz_g) + dt* &
                v1bb  (trialx,lin,bzx79)
        end if
     end if
 
  ! Unsplit time-step
  else
     if(linear.eq.0) then
        tempx = v1psipsi(trialx,lin,ps179) &
             +  v1psipsi(trialx,ps179,lin)
        ssterm(:,psi_g) = ssterm(:,psi_g) -     thimp     *dt*tempx
        ddterm(:,psi_g) = ddterm(:,psi_g) + (.5-thimp*bdf)*dt*tempx

        tempx = v1psib(trialx,lin,bz179)
        ssterm(:,psi_g) = ssterm(:,psi_g) -     thimp     *dt*tempx
        ddterm(:,psi_g) = ddterm(:,psi_g) + (nv-thimp*bdf)*dt*tempx

        if(numvar.ge.2) then
           tempx = v1psib(trialx,ps179,lin) &
                +  v1bb  (trialx,lin,bz179) &
                +  v1bb  (trialx,bz179,lin)
           ssterm(:,bz_g) = ssterm(:,bz_g) -     thimp     *dt*tempx
           ddterm(:,bz_g) = ddterm(:,bz_g) + (.5-thimp*bdf)*dt*tempx
        end if
     end if

     if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
        tempx = v1psipsi(trialx,lin,ps079) &
             +  v1psipsi(trialx,ps079,lin) &
             +  v1psib  (trialx,lin,bz079)
        ssterm(:,psi_g) = ssterm(:,psi_g) -     thimp     *dt*tempx
        ddterm(:,psi_g) = ddterm(:,psi_g) + (1.-thimp*bdf)*dt*tempx

        if(numvar.ge.2) then
           tempx = v1psib(trialx,ps079,lin) &
                +  v1bb  (trialx,lin,bz079) &
                +  v1bb  (trialx,bz079,lin)
           ssterm(:,bz_g) = ssterm(:,bz_g) -     thimp     *dt*tempx
           ddterm(:,bz_g) = ddterm(:,bz_g) + (1.-thimp*bdf)*dt*tempx
        end if
     end if

     if(use_external_fields .and. linear.eq.0) then
        tempx = v1psipsi(trialx,psx79,lin) &
             +  v1psib  (trialx,lin,bzx79)
        ssterm(:,psi_g) = ssterm(:,psi_g) -     thimp     *dt*tempx
        ddterm(:,psi_g) = ddterm(:,psi_g) + (1.-thimp*bdf)*dt*tempx

        if(numvar.ge.2) then
           tempx = v1bb(trialx,lin,bzx79)
           ssterm(:,bz_g) = ssterm(:,bz_g) -     thimp     *dt*tempx
           ddterm(:,bz_g) = ddterm(:,bz_g) + (1.-thimp*bdf)*dt*tempx
        end if
     end if
  end if

  if(i3d.eq.1 .and. numvar.ge.2) then
     if(linear.eq.0) then
        tempx = v1psif(trialx,lin,bf179)
        ssterm(:,psi_g) = ssterm(:,psi_g) -     thimp_bf     *dt*tempx
        ddterm(:,psi_g) = ddterm(:,psi_g) + (.5-thimp_bf*bdf)*dt*tempx

        tempx = v1bf(trialx,lin,bf179)
        ssterm(:,bz_g) = ssterm(:,bz_g) -     thimp_bf     *dt*tempx
        ddterm(:,bz_g) = ddterm(:,bz_g) + (.5-thimp_bf*bdf)*dt*tempx

        tempx = v1psif(trialx,ps179,lin) &
             + v1bf  (trialx,bz179,lin)
        r_bf = r_bf -     thimp_bf     *dt*tempx
        q_bf = q_bf + (.5-thimp_bf*bdf)*dt*tempx
     end if
     if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
        tempx = v1psif(trialx,lin,bf079)
        ssterm(:,psi_g) = ssterm(:,psi_g) -     thimp_bf     *dt*tempx
        ddterm(:,psi_g) = ddterm(:,psi_g) + (1.-thimp_bf*bdf)*dt*tempx

        tempx = v1bf(trialx,lin,bf079)
        ssterm(:,bz_g) = ssterm(:,bz_g) -     thimp_bf     *dt*tempx
        ddterm(:,bz_g) = ddterm(:,bz_g) + (1.-thimp_bf*bdf)*dt*tempx

        tempx = v1psif(trialx,ps079,lin) &
             + v1bf  (trialx,bz079,lin)
        r_bf = r_bf -     thimp_bf     *dt*tempx
        q_bf = q_bf + (1.-thimp_bf*bdf)*dt*tempx
     end if
     if(use_external_fields .and. linear.eq.0) then
        tempx = v1psif(trialx,lin,bfx79)
        ssterm(:,psi_g) = ssterm(:,psi_g) -     thimp_bf     *dt*tempx
        ddterm(:,psi_g) = ddterm(:,psi_g) + (1.-thimp_bf*bdf)*dt*tempx

        tempx = v1bf  (trialx,bzx79,lin)
        r_bf = r_bf -     thimp_bf     *dt*tempx
        q_bf = q_bf + (1.-thimp_bf*bdf)*dt*tempx
     end if
  end if



  do i=1, dofs_per_element
     trial = trialx(:,:,i)


  ! JxB
  ! ~~~
  ! Split time-step
  if(advfield.eq.1) then
     ! parabolization terms
     temp = v1upsipsi(trial,lin,pst79,pst79) &
          + v1upsib  (trial,lin,pst79,bzt79) &
          + v1ubb    (trial,lin,bzt79,bzt79)
     ssterm(i,u_g) = ssterm(i,u_g) - thimp*thimp*dt*dt*temp
     ddterm(i,u_g) = ddterm(i,u_g) +       ththm*dt*dt*temp

     ! two-fluid contribution
     if(db .gt. 0 .and. itwofluid.eq.2) then
        temp = v1hupsi(trial,lin,pst79) & 
             + v1hub  (trial,lin,bzt79)
        if(i3d.eq.1) temp = temp + v1huf(trial,lin,bf179)
        ssterm(i,u_g) = ssterm(i,u_g) + db*thimp*dt*temp
        ddterm(i,u_g) = ddterm(i,u_g) + db*thimp*dt*temp
     endif

     if(numvar.ge.2) then
        temp = v1vpsipsi(trial,lin,pst79,pst79) &
             + v1vpsib  (trial,lin,pst79,bzt79)
        ssterm(i,vz_g) = ssterm(i,vz_g) - thimp*thimp*dt*dt*temp
        ddterm(i,vz_g) = ddterm(i,vz_g) +       ththm*dt*dt*temp

       ! two-fluid contribution
       if(db .gt. 0 .and. itwofluid.eq.2) then
          temp = v1hvpsi(trial,lin,pst79) & 
               + v1hvb  (trial,lin,bzt79)
          if(i3d.eq.1) temp = temp + v1hvf(trial,lin,bf179)
          ssterm(i,vz_g) = ssterm(i,vz_g) + db*thimp*dt*temp
          ddterm(i,vz_g) = ddterm(i,vz_g) + db*thimp*dt*temp
       endif
     endif

     if(numvar.ge.3) then
       temp = v1chipsipsi(trial,lin,pst79,pst79) &
            + v1chipsib  (trial,lin,pst79,bzt79) &
            + v1chibb    (trial,lin,bzt79,bzt79)
       ssterm(i,chi_g) = ssterm(i,chi_g) - thimp*thimp*dt*dt*temp
       ddterm(i,chi_g) = ddterm(i,chi_g) +       ththm*dt*dt*temp

       ! two-fluid contribution
       if(db .gt. 0 .and. itwofluid.eq.2) then
          temp = v1hchipsi(trial,lin,pst79) & 
               + v1hchib  (trial,lin,bzt79)
          if(i3d.eq.1) temp = temp + v1hchif(trial,lin,bf179)
          ssterm(i,chi_g) = ssterm(i,chi_g) + db*thimp*dt*temp
          ddterm(i,chi_g) = ddterm(i,chi_g) + db*thimp*dt*temp
       endif

     endif
  end if


  ! Gyroviscosity  
  ! ~~~~~~~~~~~~~
  if(gyro.eq.1 .and. db.ne.0.) then
     temp = g1u(trial,lin)*db
     ssterm(i,u_g) = ssterm(i,u_g) -     thimp     *dt*temp
     ddterm(i,u_g) = ddterm(i,u_g) + (1.-thimp*bdf)*dt*temp

     if(numvar.ge.2) then
        temp = g1v(trial,lin)*db
        ssterm(i,vz_g) = ssterm(i,vz_g) -     thimp     *dt*temp
        ddterm(i,vz_g) = ddterm(i,vz_g) + (1.-thimp*bdf)*dt*temp
     end if

     if(numvar.ge.3) then
        temp = g1chi(trial,lin)*db
        ssterm(i,chi_g) = ssterm(i,chi_g) -     thimp     *dt*temp
        ddterm(i,chi_g) = ddterm(i,chi_g) + (1.-thimp*bdf)*dt*temp
     end if
  end if

  ! Parallel Viscosity
  ! ~~~~~~~~~~~~~~~~~~
  if(amupar.ne.0.) then
     call PVV1(trial,temp79f)

     call PVS1      (lin,temp79b)
     call PVS1psipsi(lin,pst79,pst79,temp79c)
     call PVS1psib  (lin,pst79,bzt79,temp79d)
     call PVS1bb    (lin,bzt79,bzt79,temp79e)
     temp79a = temp79b + temp79c + temp79d + temp79e

     temp = int3(vip79(:,OP_1),temp79a,temp79f)
     ssterm(i,u_g) = ssterm(i,u_g) -     thimp     *dt*temp
     ddterm(i,u_g) = ddterm(i,u_g) + (1.-thimp*bdf)*dt*temp

     if(numvar.ge.2) then      
        call PVS2      (lin,temp79b)
        call PVS2psipsi(lin,pst79,pst79,temp79c)
        call PVS2psib  (lin,pst79,bzt79,temp79d)
        call PVS2bb    (lin,bzt79,bzt79,temp79e)
        temp79a = temp79b + temp79c + temp79d + temp79e

        temp = int3(vip79(:,OP_1),temp79a,temp79f)
        ssterm(i,vz_g) = ssterm(i,vz_g) -     thimp     *dt*temp
        ddterm(i,vz_g) = ddterm(i,vz_g) + (1.-thimp*bdf)*dt*temp
     endif

     if(numvar.ge.3) then
        call PVS3      (lin,temp79b)
        call PVS3psipsi(lin,pst79,pst79,temp79c)
        call PVS3psib  (lin,pst79,bzt79,temp79d)
        call PVS3bb    (lin,bzt79,bzt79,temp79e)
        temp79a = temp79b + temp79c + temp79d + temp79e

        temp = int3(vip79(:,OP_1),temp79a,temp79f)
        ssterm(i,chi_g) = ssterm(i,chi_g) -     thimp     *dt*temp
        ddterm(i,chi_g) = ddterm(i,chi_g) + (1.-thimp*bdf)*dt*temp
     endif
  endif

  ! Poloidal force term (Lucca)
  if(ipforce.ge.1) then
     temp = v1psiforce(trial,lin,for79)
     ddterm(i,psi_g) = ddterm(i,psi_g) + dt*temp
  endif


  end do
end subroutine vorticity_lin 


subroutine vorticity_nolin(trialx, r4term)

  use basic
  use m3dc1_nint
  use metricterms_new
  use particles

  implicit none

  vectype, intent(in), dimension(MAX_PTS, OP_NUM, dofs_per_element) :: trialx
  vectype, intent(out), dimension(dofs_per_element) :: r4term

  r4term = 0.

  ! JxB_ext
  ! ~~~~~~~
  if(use_external_fields .and. (eqsubtract.eq.1 .or. icsubtract.eq.1)) then 
     r4term = r4term + dt* &
          (v1psipsi(trialx,psx79,ps079) &
          +v1psib  (trialx,ps079,bzx79))
        
     if(numvar.ge.2) then
        r4term = r4term + dt* &
             v1bb  (trialx,bz079,bzx79)
     end if

     if(i3d.eq.1 .and. numvar.ge.2) then
        r4term = r4term + dt* &
             (v1psif(trialx,ps079,bfx79) &
             +v1bf  (trialx,bzx79,bf079))
     end if
  end if

#ifdef USEPARTICLES
  ! kinetic terms
  ! ~~~~~~~~~~~~~
  if(kinetic .eq. 1) then
     r4term = r4term + dt* &
                 (v1par(trialx,ppar79)    &
                + v1parb2ipsipsi(trialx,ppar79,b2i79,pstx79,pstx79)  &
                + v1parb2ipsib(trialx,ppar79,b2i79,pstx79,bztx79)    &
                - v1par(trialx,pper79)    &
                - v1parb2ipsipsi(trialx,pper79,b2i79,pstx79,pstx79)  &
                - v1parb2ipsib(trialx,pper79,b2i79,pstx79,bztx79))
  endif
#endif

  if(linear.eq.1) return


  ! density terms
  ! ~~~~~~~~~~~~~
  if(idens.eq.1 .and. eqsubtract.eq.1) then
     r4term = r4term + dt* &
          v1us(trialx,ph079,sig79)
             
     if(numvar.ge.3) then
        r4term = r4term + dt* &
             v1chis(trialx,ch079,sig79)
     endif
  endif

end subroutine vorticity_nolin


!======================================================================
! Axial Velocity Equation
!======================================================================
subroutine axial_vel_lin(trialx, lin, ssterm, ddterm, r_bf, q_bf, advfield, &
     izone)
  
  use basic
  use arrays
  use m3dc1_nint
  use metricterms_new
  use two_fluid

  implicit none

  vectype, dimension(MAX_PTS, OP_NUM, dofs_per_element), intent(in) :: trialx
  vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: lin 
  vectype, dimension(dofs_per_element,num_fields), intent(out) :: ssterm, ddterm
  vectype, dimension(dofs_per_element), intent(out) :: r_bf, q_bf

  integer, intent(in) :: advfield
  integer, intent(in) :: izone
  vectype :: temp
  vectype, dimension(dofs_per_element) :: tempx
  real :: ththm, thimp_bf
  vectype :: freq_fac

  integer :: i
  vectype, dimension(MAX_PTS, OP_NUM) :: trial

  select case(imp_mod)
  case(0)
     ththm = (1.-thimp*bdf)*thimp
  case(1)
     ththm = -bdf*thimp**2
  case(2)
     ththm = -bdf*thimp
  case default
     ththm = (1.-thimp*bdf)*thimp
  end select

  if(imp_bf.eq.1 .and. isplitstep.eq.0) then
     thimp_bf = thimp
  else
     thimp_bf = 0.
  end if

  if(itime_independent.eq.1) then
#ifdef USECOMPLEX
     freq_fac = (0,1)*frequency*dt
#else
     freq_fac = 0.
#endif
  else
     freq_fac = 1.
  end if

  ssterm = 0.
  ddterm = 0.
  q_bf = 0.
  r_bf = 0.

  if(numvar.lt.2) return

  if(istatic.eq.1 .or. izone.ne.1) then
     if(.not.surface_int) then
        tempx = intx2(trialx(:,OP_1,:),lin(:,OP_1))
        ssterm(:,vz_g) = tempx
        ddterm(:,vz_g) = tempx
     endif
     return
  endif

  ! Time Derivative
  ! ~~~~~~~~~~~~~~~
  tempx = v2vn(trialx,lin,nt79)*freq_fac
  ssterm(:,vz_g) = ssterm(:,vz_g) + tempx
  if(itime_independent.eq.0) ddterm(:,vz_g) = ddterm(:,vz_g) + tempx*bdf


  ! Viscosity
  ! ~~~~~~~~~
  tempx = v2umu(trialx,lin,vis79,vic79)
  ssterm(:,u_g) = ssterm(:,u_g) -     thimp     *dt*tempx
  ddterm(:,u_g) = ddterm(:,u_g) + (1.-thimp*bdf)*dt*tempx

  tempx = v2vmu(trialx,lin,vis79,vic79) &
       +  v2vs (trialx,lin,sig79)
  ssterm(:,vz_g) = ssterm(:,vz_g) -     thimp     *dt*tempx
  ddterm(:,vz_g) = ddterm(:,vz_g) + (1.-thimp*bdf)*dt*tempx

  if(numvar.ge.3) then
     tempx = v2chimu(trialx,lin,vis79,vic79)
     ssterm(:,chi_g) = ssterm(:,chi_g) -     thimp     *dt*tempx
     ddterm(:,chi_g) = ddterm(:,chi_g) + (1.-thimp*bdf)*dt*tempx
  end if


  ! Advection
  ! ~~~~~~~~~
  if(linear.eq.0) then 
     tempx = v2vun(trialx,vz179,lin,nt79)
     ssterm(:,u_g) = ssterm(:,u_g) -     thimp     *dt*tempx
     ddterm(:,u_g) = ddterm(:,u_g) + (.5-thimp*bdf)*dt*tempx

     tempx = v2vun(trialx,lin,ph179,nt79) &
          +  v2vvn(trialx,lin,vz179,nt79) &
          +  v2vvn(trialx,vz179,lin,nt79)
     ssterm(:,vz_g) = ssterm(:,vz_g) -     thimp     *dt*tempx
     ddterm(:,vz_g) = ddterm(:,vz_g) + (.5-thimp*bdf)*dt*tempx

     tempx = v2vchin(trialx,lin,ch179,nt79)
     ssterm(:,vz_g) = ssterm(:,vz_g) -     thimp     *dt*tempx
     ddterm(:,vz_g) = ddterm(:,vz_g) + (.5-thimp*bdf)*dt*tempx

     if(numvar.ge.3) then
        tempx = v2vchin(trialx,vz179,lin,nt79)
        ssterm(:,chi_g) = ssterm(:,chi_g) -     thimp     *dt*tempx
        ddterm(:,chi_g) = ddterm(:,chi_g) + (.5-thimp*bdf)*dt*tempx
     end if
  end if
  
  if(eqsubtract.eq.1) then
     tempx = v2vun(trialx,vz079,lin,nt79)
     ssterm(:,u_g) = ssterm(:,u_g) -     thimp     *dt*tempx
     ddterm(:,u_g) = ddterm(:,u_g) + (1.-thimp*bdf)*dt*tempx
              
     tempx = v2vun  (trialx,lin,ph079,nt79) &
          +  v2vvn  (trialx,lin,vz079,nt79) &
          +  v2vvn  (trialx,vz079,lin,nt79) &
          +  v2vchin(trialx,lin,ch079,nt79)
     ssterm(:,vz_g) = ssterm(:,vz_g) -     thimp     *dt*tempx
     ddterm(:,vz_g) = ddterm(:,vz_g) + (1.-thimp*bdf)*dt*tempx

     if(numvar.ge.3) then
        tempx = v2vchin(trialx,vz079,lin,nt79)
        ssterm(:,chi_g) = ssterm(:,chi_g) -     thimp     *dt*tempx
        ddterm(:,chi_g) = ddterm(:,chi_g) + (1.-thimp*bdf)*dt*tempx
     end if

     if(idens.eq.1) then
        ddterm(:,den_g) = ddterm(:,den_g) + dt* &
             (v2vun  (trialx,vz079,ph079,lin) &
             +v2vvn  (trialx,vz079,vz079,lin) &
             +v2vchin(trialx,vz079,ch079,lin))
     end if
  end if


  ! JxB
  ! ~~~
  ! Split time-step
  if(advfield.eq.1) then 
     ddterm(:,psi_g) = ddterm(:,psi_g) + dt* &
          (v2psipsi(trialx,lin,pss79) & 
          +v2psipsi(trialx,pss79,lin) & 
          +v2psib  (trialx,lin,bzs79))
     
     ddterm(:,bz_g) = ddterm(:,bz_g) + dt* &
          (v2psib(trialx,pss79,lin))

     if(use_external_fields .and. linear.eq.0) then
        ddterm(:,psi_g) = ddterm(:,psi_g) + dt* &
             (v2psipsi(trialx,lin,psx79))
     
        ddterm(:,bz_g) = ddterm(:,bz_g) + dt* &
             (v2psib(trialx,psx79,lin))
     end if

  ! Unsplit time-step
  else
     if(linear.eq.0) then
        tempx = v2psipsi(trialx,lin,ps179) &
             +  v2psipsi(trialx,ps179,lin) &
             +  v2psib  (trialx,lin,bz179)
        ssterm(:,psi_g) = ssterm(:,psi_g) -     thimp     *dt*tempx
        ddterm(:,psi_g) = ddterm(:,psi_g) + (.5-thimp*bdf)*dt*tempx
        
        tempx = v2psib(trialx,ps179,lin)
        ssterm(:,bz_g) = ssterm(:,bz_g) -     thimp     *dt*tempx
        ddterm(:,bz_g) = ddterm(:,bz_g) + (.5-thimp*bdf)*dt*tempx
     end if

     if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
        tempx = v2psipsi(trialx,lin,ps079) &
             +  v2psipsi(trialx,ps079,lin) &
             +  v2psib  (trialx,lin,bz079)
        ssterm(:,psi_g) = ssterm(:,psi_g) -     thimp     *dt*tempx
        ddterm(:,psi_g) = ddterm(:,psi_g) + (1.-thimp*bdf)*dt*tempx

        tempx = v2psib(trialx,ps079,lin)
        ssterm(:,bz_g) = ssterm(:,bz_g) -     thimp     *dt*tempx
        ddterm(:,bz_g) = ddterm(:,bz_g) + (1.-thimp*bdf)*dt*tempx
     end if

     if(use_external_fields .and. linear.eq.0) then
        tempx = v2psipsi(trialx,lin,psx79)
        ssterm(:,psi_g) = ssterm(:,psi_g) -     thimp     *dt*tempx
        ddterm(:,psi_g) = ddterm(:,psi_g) + (1.-thimp*bdf)*dt*tempx

        tempx = v2psib(trialx,psx79,lin)
        ssterm(:,bz_g) = ssterm(:,bz_g) -     thimp     *dt*tempx
        ddterm(:,bz_g) = ddterm(:,bz_g) + (1.-thimp*bdf)*dt*tempx
     end if
  end if

  if(i3d.eq.1) then
     if(linear.eq.0) then
        tempx = v2psif1(trialx,lin,bf179) &
             +  v2psif2(trialx,lin,bf179)
        ssterm(:,psi_g) = ssterm(:,psi_g) -     thimp_bf     *dt*tempx
        ddterm(:,psi_g) = ddterm(:,psi_g) + (.5-thimp_bf*bdf)*dt*tempx

        tempx = v2bf(trialx,lin,bf179)
        ssterm(:,bz_g) = ssterm(:,bz_g) -     thimp_bf     *dt*tempx
        ddterm(:,bz_g) = ddterm(:,bz_g) + (.5-thimp_bf*bdf)*dt*tempx

        tempx = v2psif1(trialx,ps179,lin) &
             +  v2psif2(trialx,ps179,lin) &
             +  v2bf   (trialx,bz179,lin) &
             +  v2ff   (trialx,bf179,lin) &
             +  v2ff   (trialx,lin,bf179)
        r_bf = r_bf -     thimp_bf     *dt*tempx
        q_bf = q_bf + (.5-thimp_bf*bdf)*dt*tempx
     end if

     if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
        tempx = v2psif1(trialx,lin,bf079) &
             +  v2psif2(trialx,lin,bf079)
        ssterm(:,psi_g) = ssterm(:,psi_g) -     thimp_bf     *dt*tempx
        ddterm(:,psi_g) = ddterm(:,psi_g) + (1.-thimp_bf*bdf)*dt*tempx

        tempx = v2bf(trialx,lin,bf079)
        ssterm(:,bz_g) = ssterm(:,bz_g) -     thimp_bf     *dt*tempx
        ddterm(:,bz_g) = ddterm(:,bz_g) + (1.-thimp_bf*bdf)*dt*tempx

        tempx = v2psif1(trialx,ps079,lin) &
             +  v2psif2(trialx,ps079,lin) &
             +  v2bf   (trialx,bz079,lin) &
             +  v2ff   (trialx,bf079,lin) &
             +  v2ff   (trialx,lin,bf079)
        r_bf = r_bf -     thimp_bf     *dt*tempx
        q_bf = q_bf + (1.-thimp_bf*bdf)*dt*tempx
     end if

     if(use_external_fields .and. linear.eq.0) then
        tempx = v2psif2(trialx,lin,bfx79)
        ssterm(:,psi_g) = ssterm(:,psi_g) -     thimp_bf     *dt*tempx
        ddterm(:,psi_g) = ddterm(:,psi_g) + (1.-thimp_bf*bdf)*dt*tempx

        tempx = v2bf(trialx,lin,bfx79)
        ssterm(:,bz_g) = ssterm(:,bz_g) -     thimp_bf     *dt*tempx
        ddterm(:,bz_g) = ddterm(:,bz_g) + (1.-thimp_bf*bdf)*dt*tempx

        tempx = v2psif1(trialx,psx79,lin) &
             +  v2ff   (trialx,lin,bfx79)
        r_bf = r_bf -     thimp_bf     *dt*tempx
        q_bf = q_bf + (1.-thimp_bf*bdf)*dt*tempx
     end if
  end if


  ! Grad(p)
  ! ~~~~~~~
  if(numvar.ge.3 .or. ipres.eq.1) then
     ! Split time-step
     if(advfield.eq.1) then
        ddterm(:,p_g) = ddterm(:,p_g) + dt* &
             v2p(trialx,lin)

        ! parabolization terms
        tempx = v2up(trialx,lin,pt79)
        ssterm(:,u_g) = ssterm(:,u_g) - thimp*thimp*dt*dt*tempx
        ddterm(:,u_g) = ddterm(:,u_g) +       ththm*dt*dt*tempx

        tempx = v2vp(trialx,lin,pt79)
        ssterm(:,vz_g) = ssterm(:,vz_g) - thimp*thimp*dt*dt*tempx
        ddterm(:,vz_g) = ddterm(:,vz_g) +       ththm*dt*dt*tempx

        if(numvar.ge.3) then
           tempx = v2chip(trialx,lin,pt79)
           ssterm(:,chi_g) = ssterm(:,chi_g) - thimp*thimp*dt*dt*tempx
           ddterm(:,chi_g) = ddterm(:,chi_g) +       ththm*dt*dt*tempx
        end if

     ! Unsplit time-step
     else
        if(kinetic.le.1) then
           tempx = v2p(trialx,lin)
           ssterm(:,p_g) = ssterm(:,p_g) -     thimp     *dt*tempx
           ddterm(:,p_g) = ddterm(:,p_g) + (1.-thimp*bdf)*dt*tempx
        elseif(kinetic.eq.3) then  ! full CGL model
           tempx = v2parpb2ipsipsi(trialx,lin,b2i79,pstx79,pstx79)   &
                -  v2parpb2ipsib  (trialx,lin,b2i79,pstx79,bztx79)
           ssterm(:,p_g) = ssterm(:,p_g) -      thimp*dt*tempx
           ddterm(:,p_g) = ddterm(:,p_g) + (1.-thimp*bdf)*dt*tempx
           tempx = v2parpb2ibb  (trialx,lin,b2i79,bztx79,bztx79)       &
                +  v2parpb2ipsib(trialx,lin,b2i79,pstx79,bztx79)
           ssterm(:,pe_g) = ssterm(:,pe_g) -      thimp*dt*tempx
           ddterm(:,pe_g) = ddterm(:,pe_g) + (1.-thimp*bdf)*dt*tempx
        endif
     end if
  end if



  do i=1, dofs_per_element
     trial = trialx(:,:,i)

  ! incompressible constraint for CGL (kinetic.eq.2)
  if(kinetic.eq.2) then
     temp = incupsi(trial,lin,pstx79)
     ssterm(i,u_g) = temp
     if(numvar.ge.2) then
        temp = incvb(trial,lin,bztx79)
        ssterm(i,vz_g) = temp
        if(numvar.ge.3) then
           temp = incchipsi(trial,lin,pstx79)
           ssterm(i,chi_g) = temp
        endif
     endif
     cycle
  endif


  ! JxB
  ! ~~~
  ! Split time-step
  if(advfield.eq.1) then 

     ! parabolization terms
     temp = v2upsipsi(trial,lin,pst79,pst79) &
          + v2upsib  (trial,lin,pst79,bzt79) &
          + v2ubb    (trial,lin,bzt79,bzt79)
     ssterm(i,u_g) = ssterm(i,u_g) - thimp*thimp*dt*dt*temp
     ddterm(i,u_g) = ddterm(i,u_g) +       ththm*dt*dt*temp

     ! two-fluid contribution
     if(db .gt. 0 .and. itwofluid.eq.2) then
        temp = v2hupsi(trial,lin,pst79) & 
             + v2hub  (trial,lin,bzt79)
        if(i3d.eq.1) temp = temp + v2huf(trial,lin,bf179)
        ssterm(i,u_g) = ssterm(i,u_g) + db*thimp*dt*temp
        ddterm(i,u_g) = ddterm(i,u_g) + db*thimp*dt*temp
     endif

     temp = v2vpsipsi(trial,lin,pst79,pst79) &
          + v2vpsib  (trial,lin,pst79,bzt79)
     ssterm(i,vz_g) = ssterm(i,vz_g) - thimp*thimp*dt*dt*temp
     ddterm(i,vz_g) = ddterm(i,vz_g) +       ththm*dt*dt*temp

     ! two-fluid contribution
     if(db .gt. 0 .and. itwofluid.eq.2) then
        temp = v2hvpsi(trial,lin,pst79) & 
             + v2hvb  (trial,lin,bzt79)
        if(i3d.eq.1) temp = temp + v2hvf(trial,lin,bf179)
        ssterm(i,vz_g) = ssterm(i,vz_g) + db*thimp*dt*temp
        ddterm(i,vz_g) = ddterm(i,vz_g) + db*thimp*dt*temp
     endif

     if(numvar.ge.3) then
        temp = v2chipsipsi(trial,lin,pst79,pst79) &
             + v2chipsib  (trial,lin,pst79,bzt79) &
             + v2chibb    (trial,lin,bzt79,bzt79)
        ssterm(i,chi_g) = ssterm(i,chi_g) - thimp*thimp*dt*dt*temp
        ddterm(i,chi_g) = ddterm(i,chi_g) +       ththm*dt*dt*temp

        ! two-fluid contribution
        if(db .gt. 0 .and. itwofluid.eq.2) then
           temp = v2hchipsi(trial,lin,pst79) & 
                + v2hchib  (trial,lin,bzt79)
           if(i3d.eq.1) temp = temp + v2hchif(trial,lin,bf179)
           ssterm(i,chi_g) = ssterm(i,chi_g) + db*thimp*dt*temp
           ddterm(i,chi_g) = ddterm(i,chi_g) + db*thimp*dt*temp
        endif
     end if
  end if



  ! Gyroviscosity
  ! ~~~~~~~~~~~~~
  if(gyro.eq.1 .and. db.ne.0.) then
     temp = g2u(trial,lin)*db
     ssterm(i,u_g) = ssterm(i,u_g) -     thimp     *dt*temp
     ddterm(i,u_g) = ddterm(i,u_g) + (1.-thimp*bdf)*dt*temp
     
     temp = g2v(trial,lin)*db
     ssterm(i,vz_g) = ssterm(i,vz_g) -     thimp     *dt*temp
     ddterm(i,vz_g) = ddterm(i,vz_g) + (1.-thimp*bdf)*dt*temp

     if(numvar.ge.3) then
        temp = g2chi(trial,lin)*db
        ssterm(i,chi_g) = ssterm(i,chi_g) -     thimp     *dt*temp
        ddterm(i,chi_g) = ddterm(i,chi_g) + (1.-thimp*bdf)*dt*temp
     end if
  end if


  ! Parallel Viscosity
  ! ~~~~~~~~~~~~~~~~~~
  if(amupar.ne.0.) then
     call PVV2(trial,temp79f)

     call PVS1      (lin,temp79b)
     call PVS1psipsi(lin,pst79,pst79,temp79c)
     call PVS1psib  (lin,pst79,bzt79,temp79d)
     call PVS1bb    (lin,bzt79,bzt79,temp79e)
     temp79a = temp79b + temp79c + temp79d + temp79e

     temp = int3(vip79(:,OP_1),temp79a,temp79f)
     ssterm(i,u_g) = ssterm(i,u_g) -     thimp     *dt*temp
     ddterm(i,u_g) = ddterm(i,u_g) + (1.-thimp*bdf)*dt*temp

     if(numvar.ge.2) then      
        call PVS2      (lin,temp79b)
        call PVS2psipsi(lin,pst79,pst79,temp79c)
        call PVS2psib  (lin,pst79,bzt79,temp79d)
        call PVS2bb    (lin,bzt79,bzt79,temp79e)
        temp79a = temp79b + temp79c + temp79d + temp79e

        temp = int3(vip79(:,OP_1),temp79a,temp79f)
        ssterm(i,vz_g) = ssterm(i,vz_g) -     thimp     *dt*temp
        ddterm(i,vz_g) = ddterm(i,vz_g) + (1.-thimp*bdf)*dt*temp
     endif

     if(numvar.ge.3) then
        call PVS3      (lin,temp79b)
        call PVS3psipsi(lin,pst79,pst79,temp79c)
        call PVS3psib  (lin,pst79,bzt79,temp79d)
        call PVS3bb    (lin,bzt79,bzt79,temp79e)
        temp79a = temp79b + temp79c + temp79d + temp79e

        temp = int3(vip79(:,OP_1),temp79a,temp79f)
        ssterm(i,chi_g) = ssterm(i,chi_g) -     thimp     *dt*temp
        ddterm(i,chi_g) = ddterm(i,chi_g) + (1.-thimp*bdf)*dt*temp
     endif
  endif

  end do

end subroutine axial_vel_lin

subroutine axial_vel_nolin(trialx, r4term)

  use basic
  use m3dc1_nint
  use metricterms_new
  use transport_coefficients
  use particles

  implicit none

  vectype, intent(in), dimension(MAX_PTS, OP_NUM, dofs_per_element) :: trialx
  vectype, intent(out), dimension(dofs_per_element) :: r4term
  
  r4term = 0.

  if(numvar.lt.2) return

  ! JxB_ext
  ! ~~~~~~~
  if(use_external_fields .and. (eqsubtract.eq.1 .or. icsubtract.eq.1)) then 
     r4term = r4term + dt* &
          (v2psipsi(trialx,ps079,psx79) &
          +v2psib  (trialx,psx79,bz079))
     
     if(i3d.eq.1) then
        r4term = r4term + dt* &
             (v2psif1(trialx,psx79,bf079) &
             +v2psif2(trialx,ps079,bfx79) &
             +v2bf   (trialx,bz079,bfx79) &
             +v2ff   (trialx,bf079,bfx79))
     end if
  end if


#ifdef USEPARTICLES
  ! kinetic terms
  ! ~~~~~~~~~~~~~
  if(kinetic .eq. 1) then
     r4term = r4term + dt* &
                 (v2parpb2ipsipsi(trialx,pper79,b2i79,pstx79,pstx79)   &
                - v2parpb2ipsib  (trialx,pper79,b2i79,pstx79,bztx79)   &
                + v2parpb2ibb    (trialx,ppar79,b2i79,bztx79,bztx79)   &
                + v2parpb2ipsib  (trialx,ppar79,b2i79,pstx79,bztx79))
  endif
#endif

  ! density terms
  ! ~~~~~~~~~~~~~
  if(idens.eq.1 .and. eqsubtract.eq.1) then
     r4term = r4term + dt* &
          (v2vun(trialx,vz079,ph079,n179) &
          +v2vs (trialx,vz079,sig79))
  endif

  if(momentum_source) then
     r4term = r4term + dt*intx3(trialx(:,OP_1,:),r_79,fy79(:,OP_1))
  end if

end subroutine axial_vel_nolin


!======================================================================
! Compression Equation
!======================================================================
subroutine compression_lin(trialx, lin, ssterm, ddterm, r_bf, q_bf, advfield, &
     izone)
  
  use basic
  use arrays
  use m3dc1_nint
  use metricterms_new
  use two_fluid

  implicit none

  vectype, dimension(MAX_PTS, OP_NUM, dofs_per_element), intent(in) :: trialx 
  vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: lin
  vectype, dimension(dofs_per_element, num_fields), intent(out) :: ssterm, ddterm
  vectype, dimension(dofs_per_element), intent(out) :: r_bf, q_bf
  integer, intent(in) :: advfield
  integer, intent(in) :: izone

  vectype :: temp
  real :: ththm, thimp_bf

  vectype :: freq_fac

  integer :: i
  vectype, dimension(MAX_PTS, OP_NUM) :: trial
  vectype, dimension(dofs_per_element) :: tempx

  select case(imp_mod)
  case(0)
     ththm = (1.-thimp*bdf)*thimp
  case(1)
     ththm = -bdf*thimp**2
  case(2)
     ththm = -bdf*thimp
  case default
     ththm = (1.-thimp*bdf)*thimp
  end select

  if(imp_bf.eq.1 .and. isplitstep.eq.0) then
     thimp_bf = thimp
  else
     thimp_bf = 0.
  end if

  if(itime_independent.eq.1) then
#ifdef USECOMPLEX
     freq_fac = (0,1)*frequency*dt
#else
     freq_fac = 0.
#endif
  else
     freq_fac = 1.
  end if

  ssterm = 0.
  ddterm = 0.
  q_bf = 0. 
  r_bf = 0.
                     
  if(numvar.lt.3) return

  if(istatic.eq.1) then
     if(.not.surface_int) then
        tempx = intx2(trialx(:,OP_1,:),lin(:,OP_1))
        ssterm(:,chi_g) = tempx
        ddterm(:,chi_g) = tempx
     endif
     return
  endif
        
  if(istatic.eq.3) then    !   zero out chi only
     if(.not.surface_int) then
        tempx = intx2(trialx(:,OP_1,:),lin(:,OP_1))
        ssterm(:,chi_g) = tempx
        ddterm(:,chi_g) = 0.
     endif
     return
  endif

  ! Regularization term
  ! ~~~~~~~~~~~~~~~~~~~
  if(inoslip_pol.eq.0 .and. (.not.surface_int)) then
     tempx = -regular*intx2(trialx(:,OP_1,:),lin(:,OP_1))
     ssterm(:,chi_g) = ssterm(:,chi_g) + tempx
     ddterm(:,chi_g) = ddterm(:,chi_g) + tempx*bdf
  end if

  if(izone.ne.1) then 
     if(inonormalflow.eq.0) then
        tempx = v3un(trialx,lin,nt79)
        ssterm(:,u_g) = ssterm(:,u_g) + tempx
     end if
     
     if(inoslip_pol.eq.1 .and. .not.surface_int) then
        tempx = intx2(trialx(:,OP_1,:),lin(:,OP_1))
     else
        tempx = v3chin(trialx,lin,nt79)
     end if
     ssterm(:,chi_g) = ssterm(:,chi_g) + tempx
     return
  end if

  ! Time Derivatives
  ! ~~~~~~~~~~~~~~~~
  tempx = v3un(trialx,lin,nt79)*freq_fac
  ssterm(:,u_g) = ssterm(:,u_g) + tempx
  if(itime_independent.eq.0) ddterm(:,u_g) = ddterm(:,u_g) + tempx*bdf
     
  tempx = v3chin(trialx,lin,nt79)*chiiner*freq_fac
  ssterm(:,chi_g) = ssterm(:,chi_g) + tempx
  if(itime_independent.eq.0) ddterm(:,chi_g) = ddterm(:,chi_g) + tempx*bdf


  ! Viscosity
  ! ~~~~~~~~~
  tempx = v3umu(trialx,lin,vis79,vic79) &
       +  v3us (trialx,lin,sig79)
  ssterm(:,u_g) = ssterm(:,u_g) -     thimp     *dt*tempx
  ddterm(:,u_g) = ddterm(:,u_g) + (1.-thimp*bdf)*dt*tempx

  tempx = v3vmu(trialx,lin,vis79,vic79)
  ssterm(:,vz_g) = ssterm(:,vz_g) -     thimp     *dt*tempx
  ddterm(:,vz_g) = ddterm(:,vz_g) + (1.-thimp*bdf)*dt*tempx
                     
  tempx = v3chimu(trialx,lin,vis79,vic79) &
       +  v3chis (trialx,lin,sig79)
  ssterm(:,chi_g) = ssterm(:,chi_g) -     thimp     *dt*tempx
  ddterm(:,chi_g) = ddterm(:,chi_g) + (1.-thimp*bdf)*dt*tempx


  ! Advection
  ! ~~~~~~~~~
  if(linear.eq.0) then 
     tempx = v3uun  (trialx,lin,ph179,nt79) &
          + v3uun  (trialx,ph179,lin,nt79) &
          + v3uvn  (trialx,lin,vz179,nt79) &
          + v3uchin(trialx,lin,ch179,nt79)
     ssterm(:,u_g) = ssterm(:,u_g) -     thimp     *dt*tempx
     ddterm(:,u_g) = ddterm(:,u_g) + (.5-thimp*bdf)*dt*tempx

     tempx = v3uvn  (trialx,ph179,lin,nt79) &
          + v3vvn  (trialx,lin,vz179,nt79) &
          + v3vvn  (trialx,vz179,lin,nt79) &
          + v3vchin(trialx,lin,ch179,nt79)
     ssterm(:,vz_g) = ssterm(:,vz_g) -     thimp     *dt*tempx
     ddterm(:,vz_g) = ddterm(:,vz_g) + (.5-thimp*bdf)*dt*tempx

     tempx = v3uchin  (trialx,ph179,lin,nt79) &
          + v3vchin  (trialx,vz179,lin,nt79) &
          + v3chichin(trialx,lin,ch179,nt79) &
          + v3chichin(trialx,ch179,lin,nt79)  
     ssterm(:,chi_g) = ssterm(:,chi_g) -     thimp     *dt*tempx
     ddterm(:,chi_g) = ddterm(:,chi_g) + (.5-thimp*bdf)*dt*tempx
  endif

  if(eqsubtract.eq.1) then             
     tempx = v3uun  (trialx,lin,ph079,nt79) &
          + v3uun  (trialx,ph079,lin,nt79) &
          + v3uvn  (trialx,lin,vz079,nt79) &
          + v3uchin(trialx,lin,ch079,nt79)
     ssterm(:,u_g) = ssterm(:,u_g) -     thimp     *dt*tempx
     ddterm(:,u_g) = ddterm(:,u_g) + (1.-thimp*bdf)*dt*tempx
     
     tempx = v3uvn  (trialx,ph079,lin,nt79) &
          + v3vvn  (trialx,lin,vz079,nt79) &
          + v3vvn  (trialx,vz079,lin,nt79) &
          + v3vchin(trialx,lin,ch079,nt79)
     ssterm(:,vz_g) = ssterm(:,vz_g) -     thimp     *dt*tempx
     ddterm(:,vz_g) = ddterm(:,vz_g) + (1.-thimp*bdf)*dt*tempx
     
     tempx = v3uchin  (trialx,ph079,lin,nt79) &
          + v3vchin  (trialx,vz079,lin,nt79) &
          + v3chichin(trialx,lin,ch079,nt79) &
          + v3chichin(trialx,ch079,lin,nt79)
     ssterm(:,chi_g) = ssterm(:,chi_g) -     thimp     *dt*tempx
     ddterm(:,chi_g) = ddterm(:,chi_g) + (1.-thimp*bdf)*dt*tempx

     if(idens.eq.1) then 
        ddterm(:,den_g) = ddterm(:,den_g) + dt* &
             (v3uun    (trialx,ph079,ph079,lin) &
             +v3uvn    (trialx,ph079,vz079,lin) &
             +v3uchin  (trialx,ph079,ch079,lin) &
             +v3vvn    (trialx,vz079,vz079,lin) &
             +v3vchin  (trialx,vz079,ch079,lin) &
             +v3chichin(trialx,ch079,ch079,lin))
     end if
  end if


  ! Grad(p)
  ! ~~~~~~~
  ! Split time-step
  if(advfield.eq.1) then
     ddterm(:,p_g) = ddterm(:,p_g) + dt*  &
          v3p(trialx,lin)

     ! parabolization terms
     tempx = v3up     (trialx,lin,pt79)
     ssterm(:,u_g) = ssterm(:,u_g) - thimp*thimp*dt*dt*tempx
     ddterm(:,u_g) = ddterm(:,u_g) +       ththm*dt*dt*tempx

     tempx = v3vp     (trialx,lin,pt79)
     ssterm(:,vz_g) = ssterm(:,vz_g) - thimp*thimp*dt*dt*tempx
     ddterm(:,vz_g) = ddterm(:,vz_g) +       ththm*dt*dt*tempx

     tempx = v3chip     (trialx,lin,pt79)
     ssterm(:,chi_g) = ssterm(:,chi_g) - thimp*thimp*dt*dt*tempx
     ddterm(:,chi_g) = ddterm(:,chi_g) +       ththm*dt*dt*tempx

  ! Unsplit time-step
  else
     if    (kinetic.le.1) then
        tempx = v3p(trialx,lin) 
        ssterm(:,p_g) = ssterm(:,p_g) -     thimp     *dt*tempx
        ddterm(:,p_g) = ddterm(:,p_g) + (1.-thimp*bdf)*dt*tempx
 
        ! CGL (anisotropic pressure)
     elseif(kinetic.gt.1) then
           tempx = v3par(trialx,lin)    &
                + v3parb2ipsipsi(trialx,lin,b2i79,pstx79,pstx79)  &
                + v3parb2ipsib(trialx,lin,b2i79,pstx79,bztx79)
           ssterm(:,pe_g) = ssterm(:,pe_g) -     thimp     *dt*tempx
           ddterm(:,pe_g) = ddterm(:,pe_g) + (1.-thimp*bdf)*dt*tempx
           ssterm(:,p_g)  = ssterm(:,p_g)  +     thimp     *dt*tempx
           ddterm(:,p_g)  = ddterm(:,p_g)  - (1.-thimp*bdf)*dt*tempx
     endif
  endif


  ! JxB
  ! ~~~
  ! Split time-step
  if(advfield.eq.1) then
     ddterm(:,psi_g) = ddterm(:,psi_g) + dt*  &
          (v3psipsi(trialx,lin,pss79)      & 
          +v3psipsi(trialx,pss79,lin)      &
          +v3psib  (trialx,lin,bzs79))
     
     ddterm(:,bz_g) = ddterm(:,bz_g) + dt*  &
          (v3psib(trialx,pss79,lin)      &
          +v3bb  (trialx,lin,bzs79)      &  
          +v3bb  (trialx,bzs79,lin))

     if(use_external_fields .and. linear.eq.0) then
        ddterm(:,psi_g) = ddterm(:,psi_g) + dt*  &
             (v3psipsi(trialx,psx79,lin)      &
             +v3psib  (trialx,lin,bzx79))
     
        ddterm(:,bz_g) = ddterm(:,bz_g) + dt*  &
             (v3bb  (trialx,lin,bzx79))
     end if

  ! Unsplit time-step
  else
     if(linear.eq.0) then 
        tempx = v3psipsi(trialx,lin,ps179) &
             + v3psipsi(trialx,ps179,lin) &
             + v3psib  (trialx,lin,bz179)
        ssterm(:,psi_g) = ssterm(:,psi_g) -     thimp     *dt*tempx
        ddterm(:,psi_g) = ddterm(:,psi_g) + (.5-thimp*bdf)*dt*tempx
        
        tempx = v3psib(trialx,ps179,lin) &
             + v3bb  (trialx,lin,bz179) &
             + v3bb  (trialx,bz179,lin)
        ssterm(:,bz_g) = ssterm(:,bz_g) -     thimp     *dt*tempx
        ddterm(:,bz_g) = ddterm(:,bz_g) + (.5-thimp*bdf)*dt*tempx
     end if
     if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
        tempx = v3psipsi(trialx,lin,ps079) &
             + v3psipsi(trialx,ps079,lin) &
             + v3psib  (trialx,lin,bz079)
        ssterm(:,psi_g) = ssterm(:,psi_g) -     thimp     *dt*tempx
        ddterm(:,psi_g) = ddterm(:,psi_g) + (1.-thimp*bdf)*dt*tempx

        tempx = v3psib(trialx,ps079,lin) &
             + v3bb  (trialx,lin,bz079) &
             + v3bb  (trialx,bz079,lin)
        ssterm(:,bz_g) = ssterm(:,bz_g) -     thimp     *dt*tempx
        ddterm(:,bz_g) = ddterm(:,bz_g) + (1.-thimp*bdf)*dt*tempx
     end if
     if(use_external_fields .and. linear.eq.0) then
        tempx = v3psipsi(trialx,psx79,lin) &
             + v3psib  (trialx,lin,bzx79)
        ssterm(:,psi_g) = ssterm(:,psi_g) -     thimp     *dt*tempx
        ddterm(:,psi_g) = ddterm(:,psi_g) + (1.-thimp*bdf)*dt*tempx

        tempx = v3bb  (trialx,lin,bzx79)
        ssterm(:,bz_g) = ssterm(:,bz_g) -     thimp     *dt*tempx
        ddterm(:,bz_g) = ddterm(:,bz_g) + (1.-thimp*bdf)*dt*tempx
     end if
  end if

  if(i3d.eq.1) then
     if(linear.eq.0) then
        tempx = v3psif(trialx,lin,bf179)
        ssterm(:,psi_g) = ssterm(:,psi_g) -     thimp_bf     *dt*tempx
        ddterm(:,psi_g) = ddterm(:,psi_g) + (.5-thimp_bf*bdf)*dt*tempx
        
        tempx = v3bf  (trialx,lin,bf179)
        ssterm(:,bz_g) = ssterm(:,bz_g) -     thimp_bf     *dt*tempx
        ddterm(:,bz_g) = ddterm(:,bz_g) + (.5-thimp_bf*bdf)*dt*tempx
        
        tempx = v3psif(trialx,ps179,lin) &
             + v3bf  (trialx,bz179,lin)
        r_bf = r_bf -     thimp_bf     *dt*tempx
        q_bf = q_bf + (.5-thimp_bf*bdf)*dt*tempx
     end if
     if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
        tempx = v3psif(trialx,lin,bf079)
        ssterm(:,psi_g) = ssterm(:,psi_g) -     thimp_bf     *dt*tempx
        ddterm(:,psi_g) = ddterm(:,psi_g) + (1.-thimp_bf*bdf)*dt*tempx
        
        tempx = v3bf  (trialx,lin,bf079)
        ssterm(:,bz_g) = ssterm(:,bz_g) -     thimp_bf     *dt*tempx
        ddterm(:,bz_g) = ddterm(:,bz_g) + (1.-thimp_bf*bdf)*dt*tempx        

        ! The following terms cause problems in isplitstep=0 for some reason
        ! unless imp_bf = 1
        tempx = v3psif(trialx,ps079,lin) &
             + v3bf  (trialx,bz079,lin)
        r_bf = r_bf -     thimp_bf     *dt*tempx
        q_bf = q_bf + (1.-thimp_bf*bdf)*dt*tempx
     end if
     if(use_external_fields .and. linear.eq.0) then
        tempx = v3psif(trialx,lin,bfx79)
        ssterm(:,psi_g) = ssterm(:,psi_g) -     thimp_bf     *dt*tempx
        ddterm(:,psi_g) = ddterm(:,psi_g) + (1.-thimp_bf*bdf)*dt*tempx
        
        tempx = v3bf  (trialx,bzx79,lin)
        r_bf = r_bf -     thimp_bf     *dt*tempx
        q_bf = q_bf + (1.-thimp_bf*bdf)*dt*tempx
     end if
  endif



  do i=1, dofs_per_element
     trial = trialx(:,:,i)


  ! JxB
  ! ~~~
  ! Split time-step
  if(advfield.eq.1) then

     ! parabolization terms
     temp = v3upsipsi(trial,lin,pst79,pst79) &
          + v3upsib  (trial,lin,pst79,bzt79) &
          + v3ubb    (trial,lin,bzt79,bzt79)
     ssterm(i,u_g) = ssterm(i,u_g) - thimp*thimp*dt*dt*temp
     ddterm(i,u_g) = ddterm(i,u_g) +       ththm*dt*dt*temp

     ! two-fluid contribution
     if(db .gt. 0 .and. itwofluid.gt.1) then
        temp = v3hupsi(trial,lin,pst79) & 
             + v3hub  (trial,lin,bzt79)
        if(i3d.eq.1) temp = temp + v3huf(trial,lin,bf179)
        ssterm(i,u_g) = ssterm(i,u_g) + db*thimp*dt*temp
        ddterm(i,u_g) = ddterm(i,u_g) + db*thimp*dt*temp
     endif

     temp = v3vpsipsi(trial,lin,pst79,pst79) &
          + v3vpsib  (trial,lin,pst79,bzt79) &
          + v3vbb    (trial,lin,bzt79,bzt79)
     ssterm(i,vz_g) = ssterm(i,vz_g) - thimp*thimp*dt*dt*temp
     ddterm(i,vz_g) = ddterm(i,vz_g) +       ththm*dt*dt*temp

     ! two-fluid contribution
     if(db .gt. 0 .and. itwofluid.gt.1) then
        temp = v3hvpsi(trial,lin,pst79) & 
             + v3hvb  (trial,lin,bzt79)
        if(i3d.eq.1) temp = temp + v3hvf(trial,lin,bf179)
        ssterm(i,vz_g) = ssterm(i,vz_g) + db*thimp*dt*temp
        ddterm(i,vz_g) = ddterm(i,vz_g) + db*thimp*dt*temp
     endif


     temp = v3chipsipsi(trial,lin,pst79,pst79) &
          + v3chipsib  (trial,lin,pst79,bzt79) &
          + v3chibb    (trial,lin,bzt79,bzt79)
     ssterm(i,chi_g) = ssterm(i,chi_g) - thimp*thimp*dt*dt*temp
     ddterm(i,chi_g) = ddterm(i,chi_g) +       ththm*dt*dt*temp

     ! two-fluid contribution
     if(db .gt. 0 .and. itwofluid.gt.1) then
        temp = v3hchipsi(trial,lin,pst79) & 
             + v3hchib  (trial,lin,bzt79)
        if(i3d.eq.1) temp = temp + v3hchif(trial,lin,bf179)
        ssterm(i,chi_g) = ssterm(i,chi_g) + db*thimp*dt*temp
        ddterm(i,chi_g) = ddterm(i,chi_g) + db*thimp*dt*temp
     endif
  endif


  ! Gyroviscisity
  ! ~~~~~~~~~~~~~
  if(gyro.eq.1 .and. db.ne.0.) then    
     temp = g3u(trial,lin)*db
     ssterm(i,u_g) = ssterm(i,u_g) -     thimp     *dt*temp
     ddterm(i,u_g) = ddterm(i,u_g) + (1.-thimp*bdf)*dt*temp
     
     temp = g3v(trial,lin)*db
     ssterm(i,vz_g) = ssterm(i,vz_g) -     thimp     *dt*temp
     ddterm(i,vz_g) = ddterm(i,vz_g) + (1.-thimp*bdf)*dt*temp
     
     temp = g3chi(trial,lin)*db
     ssterm(i,chi_g) = ssterm(i,chi_g) -     thimp     *dt*temp
     ddterm(i,chi_g) = ddterm(i,chi_g) + (1.-thimp*bdf)*dt*temp
  endif
    

  ! Parallel Viscosity
  ! ~~~~~~~~~~~~~~~~~~
  if(amupar.ne.0.) then
     call PVV3(trial,temp79f)

     call PVS1      (lin,temp79b)
     call PVS1psipsi(lin,pst79,pst79,temp79c)
     call PVS1psib  (lin,pst79,bzt79,temp79d)
     call PVS1bb    (lin,bzt79,bzt79,temp79e)
     temp79a = temp79b + temp79c + temp79d + temp79e

     temp = int3(vip79(:,OP_1),temp79a,temp79f)
     ssterm(i,u_g) = ssterm(i,u_g) -     thimp     *dt*temp
     ddterm(i,u_g) = ddterm(i,u_g) + (1.-thimp*bdf)*dt*temp

     if(numvar.ge.2) then      
        call PVS2      (lin,temp79b)
        call PVS2psipsi(lin,pst79,pst79,temp79c)
        call PVS2psib  (lin,pst79,bzt79,temp79d)
        call PVS2bb    (lin,bzt79,bzt79,temp79e)
        temp79a = temp79b + temp79c + temp79d + temp79e

        temp = int3(vip79(:,OP_1),temp79a,temp79f)
        ssterm(i,vz_g) = ssterm(i,vz_g) -     thimp     *dt*temp
        ddterm(i,vz_g) = ddterm(i,vz_g) + (1.-thimp*bdf)*dt*temp
     endif

     if(numvar.ge.3) then
        call PVS3      (lin,temp79b)
        call PVS3psipsi(lin,pst79,pst79,temp79c)
        call PVS3psib  (lin,pst79,bzt79,temp79d)
        call PVS3bb    (lin,bzt79,bzt79,temp79e)
        temp79a = temp79b + temp79c + temp79d + temp79e

        temp = int3(vip79(:,OP_1),temp79a,temp79f)
        ssterm(i,chi_g) = ssterm(i,chi_g) -     thimp     *dt*temp
        ddterm(i,chi_g) = ddterm(i,chi_g) + (1.-thimp*bdf)*dt*temp
     endif
  endif


  ! Gravity
  ! ~~~~~~~
  if(idens.eq.1) then
     ! Split time-step
     if(advfield.eq.1) then 
        ddterm(i,den_g) = ddterm(i,den_g) + dt* &
             v3ngrav(trial,lin)

        ! parabolization terms
        temp = v3ungrav (trial,lin,nt79)
        ssterm(i,u_g) = ssterm(i,u_g) - thimp*thimp*dt*dt*temp
        ddterm(i,u_g) = ddterm(i,u_g) +       ththm*dt*dt*temp

        temp = v3chingrav (trial,lin,nt79)
        ssterm(i,chi_g) = ssterm(i,chi_g) - thimp*thimp*dt*dt*temp
        ddterm(i,chi_g) = ddterm(i,chi_g) +       ththm*dt*dt*temp

     ! Unsplit time-step
     else
        temp = v3ngrav(trial,lin)
        ssterm(i,den_g) = ssterm(i,den_g) -     thimp     *dt*temp
        ddterm(i,den_g) = ddterm(i,den_g) + (1.-thimp*bdf)*dt*temp
     endif
  endif

  ! Poloidal force term (Lucca)
  ! ~~~~~~~~~~~~~~~~~
  if(ipforce.ge.1) then
     temp = v3psiforce(trial,lin,for79)
     ddterm(i,psi_g) = ddterm(i,psi_g) + dt*temp
  endif

  end do
end subroutine compression_lin

subroutine compression_nolin(trialx, r4term)

  use basic
  use m3dc1_nint
  use metricterms_new
  use particles

  implicit none

  vectype, intent(in), dimension(MAX_PTS, OP_NUM, dofs_per_element) :: trialx
  vectype, intent(out), dimension(dofs_per_element) :: r4term

  r4term = 0.

  if(numvar.lt.3) return

  ! density terms
  ! ~~~~~~~~~~~~~
  if(idens.eq.1 .and. eqsubtract.eq.1) then
     r4term = r4term + dt* &
          (v3us     (trialx,ph079,sig79) &
          +v3chis   (trialx,ch079,sig79))
  endif


  ! JxB_ext
  ! ~~~~~~~
  if(use_external_fields .and. (eqsubtract.eq.1 .or. icsubtract.eq.1)) then 
     r4term = r4term + dt* &
          (v3psipsi(trialx,psx79,ps079) &
          +v3psib  (trialx,ps079,bzx79) &
          +v3bb    (trialx,bz079,bzx79))
     
     if(i3d.eq.1) then
        r4term = r4term + dt* &
             (v3psif(trialx,ps079,bfx79) &
             +v3bf  (trialx,bzx79,bf079))
     end if
  endif


#ifdef USEPARTICLES
  ! kinetic terms
  ! ~~~~~~~~~~~~~
  if(kinetic .eq. 1) then
     r4term = r4term + dt* &
                 (v3par(trialx,ppar79)    &
                + v3parb2ipsipsi(trialx,ppar79,b2i79,pstx79,pstx79)  &
                + v3parb2ipsib(trialx,ppar79,b2i79,pstx79,bztx79)    &
                - v3par(trialx,pper79)    &
                - v3parb2ipsipsi(trialx,pper79,b2i79,pstx79,pstx79)  &
                - v3parb2ipsib(trialx,pper79,b2i79,pstx79,bztx79))
  endif
#endif

end subroutine compression_nolin


!======================================================================
! Flux Equation
!======================================================================
subroutine flux_lin(trialx, lin, ssterm, ddterm, q_ni, r_bf, q_bf, izone)
  
  use basic
  use arrays
  use m3dc1_nint
  use metricterms_new
  use two_fluid
  use harned_mikic_mod
  use bootstrap

  implicit none

  vectype, dimension(MAX_PTS, OP_NUM, dofs_per_element), intent(in) :: trialx
  vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: lin 
  vectype, dimension(dofs_per_element,num_fields), intent(out) :: ssterm, ddterm
  vectype, intent(out), dimension(dofs_per_element) :: r_bf, q_bf 
  vectype, intent(out), dimension(dofs_per_element,2) :: q_ni
  integer, intent(in) :: izone

  vectype, dimension(dofs_per_element) :: tempx
  vectype, dimension(MAX_PTS, OP_NUM) :: trial
  integer :: i

  vectype :: temp, temp2
  real :: thimpb, thimpf, thimpe, thimp_bf, nv

  vectype :: freq_fac

  thimpf = thimp
  thimpe = 1.

  if(imp_mod.eq.0) then
     thimpb = thimp
  else
     thimpb = 1.
  endif

  if(imp_bf.eq.1) then
     thimp_bf = thimpf
  else
     thimp_bf = 0.
  end if

  if(numvar.eq.1) then
     nv = 1.
  else
     nv = 0.5
  endif

  if(itime_independent.eq.1) then
#ifdef USECOMPLEX
     freq_fac = (0,1)*frequency*dt
#else
     freq_fac = 0.
#endif
  else
     freq_fac = 1.
  end if

  ssterm = 0.
  ddterm = 0.
  q_ni = 0.
  r_bf = 0.
  q_bf = 0.


  if(iestatic.eq.1) then
     if(.not.surface_int) then
        tempx = intx2(trialx(:,OP_1,:),lin(:,OP_1))
        ssterm(:,psi_g) = tempx
        ddterm(:,psi_g) = tempx
     endif
     return
  end if

  ! Regularization term
  ! ~~~~~~~~~~~~~~~~~~~
  if(iconst_bn.eq.0 .and. (.not.surface_int)) then
     tempx = -regular*intx2(trialx(:,OP_1,:),lin(:,OP_1))
     ssterm(:,psi_g) = ssterm(:,psi_g) + tempx
     if(itime_independent.eq.0) ddterm(:,psi_g) = ddterm(:,psi_g) + tempx*bdf
  end if

  ! Resistive and Hyper Terms
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~
  tempx = b1psieta(trialx,lin,eta79,vz079,eta_mod.eq.1)
  ssterm(:,psi_g) = ssterm(:,psi_g) -     thimp     *dt*tempx
  ddterm(:,psi_g) = ddterm(:,psi_g) + (1.-thimp*bdf)*dt*tempx
  ! implicit hyperresistivity
  if(jadv.eq.1 .and. imp_hyper.eq.1) then
     tempx = b1jeta(trialx,lin,eta79)
     ssterm(:,e_g) = ssterm(:,e_g) - dt*tempx
  endif

  if(numvar.ge.2) then
     tempx = b1beta(trialx,lin,eta79)
     ssterm(:,bz_g) = ssterm(:,bz_g) -     thimp     *dt*tempx
     ddterm(:,bz_g) = ddterm(:,bz_g) + (1.-thimp*bdf)*dt*tempx
     ! implicit hyperresistivity
     if(jadv.eq.1 .and. imp_hyper.eq.2) then
        tempx = b1bj(trialx,bzt79,lin) + b1psij(trialx,pst79,lin)
        ssterm(:,e_g) = ssterm(:,e_g) - dt*tempx
     endif

     if(i3d.eq.1) then
        tempx = b1feta(trialx,lin,eta79)
        r_bf = r_bf -     thimp_bf     *dt*tempx
        q_bf = q_bf + (1.-thimp_bf*bdf)*dt*tempx
        ! implicit hyperrestivity
        if(jadv.eq.1 .and. imp_hyper.eq.2) then
           tempx = b1fj(trialx,bft79,lin)
           ssterm(:,e_g) = ssterm(:,e_g) - dt*tempx
        endif
     end if
  endif


  ! Zone 3: eta J = 0.
  if(izone.eq.3) return


  ! Time Derivatives
  ! ~~~~~~~~~~~~~~~~
  tempx = (b1psi (trialx,lin) &
       -  b1psid(trialx,lin,ni79)) &  ! electron mass term
       * freq_fac
  ssterm(:,psi_g) = ssterm(:,psi_g) + tempx
  if(itime_independent.eq.0) ddterm(:,psi_g) = ddterm(:,psi_g) + tempx*bdf

  ! Zone 2: E = eta J
  if(izone.ne.1) return


  ! VxB
  ! ~~~
  if(linear.eq.0) then 
     tempx = b1psiu  (trialx,lin,ph179) &
          + b1psiv  (trialx,lin,vz179) &
          + b1psichi(trialx,lin,ch179)
     ssterm(:,psi_g) = ssterm(:,psi_g) -     thimpb     *dt*tempx
     ddterm(:,psi_g) = ddterm(:,psi_g) + (.5-thimpb*bdf)*dt*tempx

     tempx = b1psiu  (trialx,ps179,lin)
     ssterm(:,u_g) = ssterm(:,u_g) -     thimpb     *dt*tempx
     ddterm(:,u_g) = ddterm(:,u_g) + (.5-thimpb*bdf)*dt*tempx

     tempx = b1bu    (trialx,bz179,lin)
     ssterm(:,u_g) = ssterm(:,u_g) -     thimp     *dt*tempx
     ddterm(:,u_g) = ddterm(:,u_g) + (nv-thimp*bdf)*dt*tempx

     if(numvar.ge.2) then
        tempx = b1bu  (trialx,lin,ph179) &
             + b1bv  (trialx,lin,vz179) &
             + b1bchi(trialx,lin,ch179)
        ssterm(:,bz_g) = ssterm(:,bz_g) -     thimpb     *dt*tempx
        ddterm(:,bz_g) = ddterm(:,bz_g) + (.5-thimpb*bdf)*dt*tempx

        tempx = b1psiv(trialx,ps179,lin) &
             + b1bv  (trialx,bz179,lin)  
        ssterm(:,vz_g) = ssterm(:,vz_g) -     thimpb     *dt*tempx
        ddterm(:,vz_g) = ddterm(:,vz_g) + (.5-thimpb*bdf)*dt*tempx
     end if

     if(numvar.ge.3) then
        tempx = b1psichi(trialx,ps179,lin) &
             + b1bchi  (trialx,bz179,lin)
        ssterm(:,chi_g) = ssterm(:,chi_g) -     thimpb     *dt*tempx
        ddterm(:,chi_g) = ddterm(:,chi_g) + (.5-thimpb*bdf)*dt*tempx
     endif
  end if

  if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
     tempx = b1psiu  (trialx,lin,ph079) &
          + b1psiv  (trialx,lin,vz079) &
          + b1psichi(trialx,lin,ch079)
     ssterm(:,psi_g) = ssterm(:,psi_g) -     thimpb     *dt*tempx
     ddterm(:,psi_g) = ddterm(:,psi_g) + (1.-thimpb*bdf)*dt*tempx

     tempx = b1psiu  (trialx,ps079,lin) &
          + b1bu    (trialx,bz079,lin)
     ssterm(:,u_g) = ssterm(:,u_g) -     thimpb     *dt*tempx
     ddterm(:,u_g) = ddterm(:,u_g) + (1.-thimpb*bdf)*dt*tempx

     if(numvar.ge.2) then
        tempx = b1bu  (trialx,lin,ph079) &
             + b1bv  (trialx,lin,vz079) &
             + b1bchi(trialx,lin,ch079)
        ssterm(:,bz_g) = ssterm(:,bz_g) -     thimpb     *dt*tempx
        ddterm(:,bz_g) = ddterm(:,bz_g) + (1.-thimpb*bdf)*dt*tempx

        tempx = b1psiv(trialx,ps079,lin) &
             + b1bv  (trialx,bz079,lin)
        ssterm(:,vz_g) = ssterm(:,vz_g) -     thimpb     *dt*tempx
        ddterm(:,vz_g) = ddterm(:,vz_g) + (1.-thimpb*bdf)*dt*tempx
     end if

     if(numvar.ge.3) then
        tempx = b1psichi(trialx,ps079,lin) &
             + b1bchi  (trialx,bz079,lin)
        ssterm(:,chi_g) = ssterm(:,chi_g) -     thimpb     *dt*tempx
        ddterm(:,chi_g) = ddterm(:,chi_g) + (1.-thimpb*bdf)*dt*tempx
     endif
  end if

  if(use_external_fields .and. linear.eq.0) then
     tempx = b1psiu  (trialx,psx79,lin) &
          + b1bu    (trialx,bzx79,lin)
     ssterm(:,u_g) = ssterm(:,u_g) -     thimpb     *dt*tempx
     ddterm(:,u_g) = ddterm(:,u_g) + (1.-thimpb*bdf)*dt*tempx

     if(numvar.ge.2) then
        tempx = b1psiv(trialx,psx79,lin) &
             + b1bv  (trialx,bzx79,lin)
        ssterm(:,vz_g) = ssterm(:,vz_g) -     thimpb     *dt*tempx
        ddterm(:,vz_g) = ddterm(:,vz_g) + (1.-thimpb*bdf)*dt*tempx
     end if

     if(numvar.ge.3) then
        tempx = b1psichi(trialx,psx79,lin) &
             + b1bchi  (trialx,bzx79,lin)
        ssterm(:,chi_g) = ssterm(:,chi_g) -     thimpb     *dt*tempx
        ddterm(:,chi_g) = ddterm(:,chi_g) + (1.-thimpb*bdf)*dt*tempx
     endif
  end if


  if(i3d.eq.1 .and. numvar.ge.2) then
     if(linear.eq.0) then
        tempx = b1fu(trialx,bf179,lin)
        ssterm(:,u_g) = ssterm(:,u_g) -     thimp_bf     *dt*tempx
        ddterm(:,u_g) = ddterm(:,u_g) + (.5-thimp_bf*bdf)*dt*tempx
     
        tempx = b1fv(trialx,bf179,lin)
        ssterm(:,vz_g) = ssterm(:,vz_g) -     thimp_bf     *dt*tempx
        ddterm(:,vz_g) = ddterm(:,vz_g) + (.5-thimp_bf*bdf)*dt*tempx

        if(numvar.ge.3) then
           tempx = b1fchi(trialx,bf179,lin)
           ssterm(:,chi_g) = ssterm(:,chi_g) -     thimp_bf     *dt*tempx
           ddterm(:,chi_g) = ddterm(:,chi_g) + (.5-thimp_bf*bdf)*dt*tempx
        end if
        
        tempx = b1fu  (trialx,lin,ph179) &
             + b1fv  (trialx,lin,vz179) &
             + b1fchi(trialx,lin,ch179)
        r_bf = r_bf -     thimp_bf     *dt*tempx
        q_bf = q_bf + (.5-thimp_bf*bdf)*dt*tempx
     end if
     if(eqsubtract.eq.1) then
        tempx = b1fu(trialx,bf079,lin)
        ssterm(:,u_g) = ssterm(:,u_g) -     thimp_bf     *dt*tempx
        ddterm(:,u_g) = ddterm(:,u_g) + (1.-thimp_bf*bdf)*dt*tempx
     
        tempx = b1fv(trialx,bf079,lin)
        ssterm(:,vz_g) = ssterm(:,vz_g) -     thimp_bf     *dt*tempx
        ddterm(:,vz_g) = ddterm(:,vz_g) + (1.-thimp_bf*bdf)*dt*tempx

        if(numvar.ge.3) then
           tempx = b1fchi(trialx,bf079,lin)
           ssterm(:,chi_g) = ssterm(:,chi_g) -     thimp_bf     *dt*tempx
           ddterm(:,chi_g) = ddterm(:,chi_g) + (1.-thimp_bf*bdf)*dt*tempx
        end if
        
        tempx = b1fu  (trialx,lin,ph079) &
             + b1fv  (trialx,lin,vz079) &
             + b1fchi(trialx,lin,ch079)
        r_bf = r_bf -     thimp_bf     *dt*tempx
        q_bf = q_bf + (1.-thimp_bf*bdf)*dt*tempx
     end if
     if(use_external_fields .and. linear.eq.0) then
        tempx = b1fu(trialx,bfx79,lin)
        ssterm(:,u_g) = ssterm(:,u_g) -     thimp_bf     *dt*tempx
        ddterm(:,u_g) = ddterm(:,u_g) + (1.-thimp_bf*bdf)*dt*tempx
     
        tempx = b1fv(trialx,bfx79,lin)
        ssterm(:,vz_g) = ssterm(:,vz_g) -     thimp_bf     *dt*tempx
        ddterm(:,vz_g) = ddterm(:,vz_g) + (1.-thimp_bf*bdf)*dt*tempx

        if(numvar.ge.3) then
           tempx = b1fchi(trialx,bfx79,lin)
           ssterm(:,chi_g) = ssterm(:,chi_g) -     thimp_bf     *dt*tempx
           ddterm(:,chi_g) = ddterm(:,chi_g) + (1.-thimp_bf*bdf)*dt*tempx
        end if
     end if
  end if


  do i=1, dofs_per_element
     trial = trialx(:,:,i)


  if(jadv.eq.0) then
  ! electrostatic potential
  ! ~~~~~~~~~~~~~~~~~~~~~~~
    temp = b1e(trial,lin)
    ssterm(i,e_g) = ssterm(i,e_g)       - thimpe     *dt*temp
    ddterm(i,e_g) = ddterm(i,e_g) + (1. - thimpe*bdf)*dt*temp
  endif

  select case (itwofluid)

  case(-1) !  electron form with extra density linearization
  ! JxB
  ! ~~~
  if(db.ne.0.) then
     if(linear.eq.0) then
        temp = b1psipsin(trial,lin,ps179,nt79)*db &
             + b1psipsin(trial,ps179,lin,nt79)*db
        ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
        ddterm(i,psi_g) = ddterm(i,psi_g) + (.5-thimpf*bdf)*dt*temp
        
        temp = b1psibn1(trial,lin,bz179,nt79)*db &
             + b1psibn2(trial,lin,bz179,nt79)*db
        ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
        ddterm(i,psi_g) = ddterm(i,psi_g) + (.5-thimpf*bdf)*dt*temp

        if(numvar.ge.2) then
           temp = b1psibn1(trial,ps179,lin,nt79)*db &
                + b1psibn2(trial,ps179,lin,nt79)*db &
                + b1bbn   (trial,bz179,lin,nt79)*db &
                + b1bbn   (trial,lin,bz179,nt79)*db
           ssterm(i,bz_g) = ssterm(i,bz_g) -     thimpf     *dt*temp
           ddterm(i,bz_g) = ddterm(i,bz_g) + (.5-thimpf*bdf)*dt*temp
        end if
     end if

     if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
        temp = b1psipsin(trial,lin,ps079,nt79)*db &
             + b1psipsin(trial,ps079,lin,nt79)*db &
             + b1psibn1 (trial,lin,bz079,nt79)*db &
             + b1psibn2 (trial,lin,bz079,nt79)*db 
        ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
        ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimpf*bdf)*dt*temp

        if(numvar.ge.2) then
           temp = b1psibn1(trial,ps079,lin,nt79)*db &
                + b1psibn2(trial,ps079,lin,nt79)*db &
                + b1bbn   (trial,bz079,lin,nt79)*db & 
                + b1bbn   (trial,lin,bz079,nt79)*db   
           ssterm(i,bz_g) = ssterm(i,bz_g) -     thimpf     *dt*temp
           ddterm(i,bz_g) = ddterm(i,bz_g) + (1.-thimpf*bdf)*dt*temp
        end if
     end if

     if(use_external_fields .and. linear.eq.0) then
        temp = b1psipsin(trial,psx79,lin,nt79)*db &
             + b1psibn2 (trial,lin,bzx79,nt79)*db 
        ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
        ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimpf*bdf)*dt*temp

        if(numvar.ge.2) then
           temp = b1psibn1(trial,psx79,lin,nt79)*db &
                + b1bbn   (trial,lin,bzx79,nt79)*db
           ssterm(i,bz_g) = ssterm(i,bz_g) -     thimpf     *dt*temp
           ddterm(i,bz_g) = ddterm(i,bz_g) + (1.-thimpf*bdf)*dt*temp
        end if
     end if

     if(i3d.eq.1 .and. numvar.ge.2) then
        if(linear.eq.0) then
           temp = b1psifn1(trial,lin,bf179,nt79)*db &
                + b1psifn2(trial,lin,bf179,nt79)*db
           ssterm(i,psi_g) = ssterm(i,psi_g) -     thimp_bf     *dt*temp
           ddterm(i,psi_g) = ddterm(i,psi_g) + (.5-thimp_bf*bdf)*dt*temp

           temp = b1bfn1  (trial,lin,bf179,nt79)*db &
                + b1bfn2  (trial,lin,bf179,nt79)*db
           ssterm(i,bz_g) = ssterm(i,bz_g) -     thimp_bf     *dt*temp
           ddterm(i,bz_g) = ddterm(i,bz_g) + (.5-thimp_bf*bdf)*dt*temp

           temp = b1psifn1(trial,ps179,lin,nt79)*db &
                + b1psifn2(trial,ps179,lin,nt79)*db &
                + b1bfn1  (trial,bz179,lin,nt79)*db &
                + b1bfn2  (trial,bz179,lin,nt79)*db
           r_bf(i) = r_bf(i) -     thimp_bf     *dt*temp
           q_bf(i) = q_bf(i) + (.5-thimp_bf*bdf)*dt*temp
        end if
        if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
           temp = b1psifn1(trial,lin,bf079,nt79)*db &
                + b1psifn2(trial,lin,bf079,nt79)*db
           ssterm(i,psi_g) = ssterm(i,psi_g) -     thimp_bf     *dt*temp
           ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimp_bf*bdf)*dt*temp

           temp = b1bfn1  (trial,lin,bf079,nt79)*db &
                + b1bfn2  (trial,lin,bf079,nt79)*db
           ssterm(i,bz_g) = ssterm(i,bz_g) -     thimp_bf     *dt*temp
           ddterm(i,bz_g) = ddterm(i,bz_g) + (1.-thimp_bf*bdf)*dt*temp

           temp = b1psifn1(trial,ps079,lin,nt79)*db &
                + b1psifn2(trial,ps079,lin,nt79)*db &
                + b1bfn1  (trial,bz079,lin,nt79)*db &
                + b1bfn2  (trial,bz079,lin,nt79)*db
           r_bf(i) = r_bf(i) -     thimp_bf     *dt*temp
           q_bf(i) = q_bf(i) + (1.-thimp_bf*bdf)*dt*temp
        end if
        if(use_external_fields .and. linear.eq.0) then
           temp = b1psifn1(trial,lin,bfx79,nt79)*db
           ssterm(i,psi_g) = ssterm(i,psi_g) -     thimp_bf     *dt*temp
           ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimp_bf*bdf)*dt*temp

           temp = b1bfn1  (trial,lin,bfx79,nt79)*db
           ssterm(i,bz_g) = ssterm(i,bz_g) -     thimp_bf     *dt*temp
           ddterm(i,bz_g) = ddterm(i,bz_g) + (1.-thimp_bf*bdf)*dt*temp

           temp = b1psifn2(trial,psx79,lin,nt79)*db &
                + b1bfn2  (trial,bzx79,lin,nt79)*db
           r_bf(i) = r_bf(i) -     thimp_bf     *dt*temp
           q_bf(i) = q_bf(i) + (1.-thimp_bf*bdf)*dt*temp
        end if
     end if

     if(idens.eq.1) then
        temp = b1psipsin(trial,pstx79,pst79,lin)*db &
             + b1psibn1 (trial,pstx79,bzt79,lin)*db &
             + b1psibn2 (trial,pst79,bztx79,lin)*db &
             + b1bbn    (trial,bzt79,bztx79,lin)*db
        if(i3d.eq.1 .and. numvar.ge.2) then 
           temp = temp &
                +b1psifn1(trial,pst79,bftx79,lin)*db &
                +b1psifn2(trial,pstx79,bft79,lin)*db &
                +b1bfn1  (trial,bzt79,bftx79,lin)*db &
                +b1bfn2  (trial,bztx79,bft79,lin)*db
        end if
        ssterm(i,den_g) = ssterm(i,den_g) + thimpf*dt*temp
        ddterm(i,den_g) = ddterm(i,den_g) + thimpf*dt*temp

        if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
           temp = b1psipsin(trial,ps079,ps079,lin)*db &
                + b1psibn1 (trial,ps079,bz079,lin)*db &
                + b1psibn2 (trial,ps079,bz079,lin)*db &
                + b1bbn    (trial,bz079,bz079,lin)*db
           if(i3d.eq.1 .and. numvar.ge.2) then 
              temp = temp &
                   +b1psifn1(trial,ps079,bf079,lin)*db &
                   +b1psifn2(trial,ps079,bf079,lin)*db &
                   +b1bfn1  (trial,bz079,bf079,lin)*db &
                   +b1bfn2  (trial,bz079,bf079,lin)*db
           end if
           if(linear.eq.1) temp = -temp
           ddterm(i,den_g) = ddterm(i,den_g) + dt*temp
        end if
     end if
  end if

  ! Grad(p_e)
  ! ~~~~~~~~~
  if(db.ne.0.) then
     if(numvar.ge.3) then
        temp = b1pen(trial,lin,nt79)*db
        ssterm(i,pe_g) = ssterm(i,pe_g) -     thimpf     *dt*temp
        ddterm(i,pe_g) = ddterm(i,pe_g) + (1.-thimpf*bdf)*dt*temp
     end if

     if(idens.eq.1) then
        temp = b1pen(trial,pet79,lin)*db
        ssterm(i,den_g) = ssterm(i,den_g) + thimpf*dt*temp
        ddterm(i,den_g) = ddterm(i,den_g) + thimpf*dt*temp

        if(eqsubtract.eq.1) then 
           temp = b1pen(trial,pe079,lin)*db
           if(linear.eq.1) temp = -temp
           ddterm(i,den_g) = ddterm(i,den_g) + dt*temp
        end if
     end if
  end if  


  case(1)     ! electron form
  ! JxB
  ! ~~~
  if(db.ne.0.) then
     if(linear.eq.0) then
        temp = b1psipsid(trial,lin,ps179,ni79)*db &
             + b1psipsid(trial,ps179,lin,ni79)*db
        ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
        ddterm(i,psi_g) = ddterm(i,psi_g) + (.5-thimpf*bdf)*dt*temp
        
        temp = b1psibd1(trial,lin,bz179,ni79)*db &
             + b1psibd2(trial,lin,bz179,ni79)*db
        ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
        ddterm(i,psi_g) = ddterm(i,psi_g) + (.5-thimpf*bdf)*dt*temp

        if(numvar.ge.2) then
           temp = b1psibd1(trial,ps179,lin,ni79)*db &
                + b1psibd2(trial,ps179,lin,ni79)*db &
                + b1bbd   (trial,bz179,lin,ni79)*db &
                + b1bbd   (trial,lin,bz179,ni79)*db
           ssterm(i,bz_g) = ssterm(i,bz_g) -     thimpf     *dt*temp
           ddterm(i,bz_g) = ddterm(i,bz_g) + (.5-thimpf*bdf)*dt*temp
        end if
     end if

     if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
        temp = b1psipsid(trial,lin,ps079,ni79)*db &
             + b1psipsid(trial,ps079,lin,ni79)*db &
             + b1psibd1 (trial,lin,bz079,ni79)*db &
             + b1psibd2 (trial,lin,bz079,ni79)*db 
        ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
        ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimpf*bdf)*dt*temp

        if(numvar.ge.2) then
           temp = b1psibd1(trial,ps079,lin,ni79)*db &
                + b1psibd2(trial,ps079,lin,ni79)*db &
                + b1bbd   (trial,bz079,lin,ni79)*db & 
                + b1bbd   (trial,lin,bz079,ni79)*db   
           ssterm(i,bz_g) = ssterm(i,bz_g) -     thimpf     *dt*temp
           ddterm(i,bz_g) = ddterm(i,bz_g) + (1.-thimpf*bdf)*dt*temp
        end if
     end if

     if(use_external_fields .and. linear.eq.0) then
        temp = b1psipsid(trial,psx79,lin,ni79)*db &
             + b1psibd2 (trial,lin,bzx79,ni79)*db 
        ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
        ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimpf*bdf)*dt*temp

        if(numvar.ge.2) then
           temp = b1psibd1(trial,psx79,lin,ni79)*db &
                + b1bbd   (trial,lin,bzx79,ni79)*db
           ssterm(i,bz_g) = ssterm(i,bz_g) -     thimpf     *dt*temp
           ddterm(i,bz_g) = ddterm(i,bz_g) + (1.-thimpf*bdf)*dt*temp
        end if
     end if

     if(i3d.eq.1 .and. numvar.ge.2) then
        if(linear.eq.0) then
           temp = b1psifd1(trial,lin,bf179,ni79)*db &
                + b1psifd2(trial,lin,bf179,ni79)*db
           ssterm(i,psi_g) = ssterm(i,psi_g) -     thimp_bf     *dt*temp
           ddterm(i,psi_g) = ddterm(i,psi_g) + (.5-thimp_bf*bdf)*dt*temp

           temp = b1bfd1  (trial,lin,bf179,ni79)*db &
                + b1bfd2  (trial,lin,bf179,ni79)*db
           ssterm(i,bz_g) = ssterm(i,bz_g) -     thimp_bf     *dt*temp
           ddterm(i,bz_g) = ddterm(i,bz_g) + (.5-thimp_bf*bdf)*dt*temp

           temp = b1psifd1(trial,ps179,lin,ni79)*db &
                + b1psifd2(trial,ps179,lin,ni79)*db &
                + b1bfd1  (trial,bz179,lin,ni79)*db &
                + b1bfd2  (trial,bz179,lin,ni79)*db
           r_bf(i) = r_bf(i) -     thimp_bf     *dt*temp
           q_bf(i) = q_bf(i) + (.5-thimp_bf*bdf)*dt*temp
        end if
        if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
           temp = b1psifd1(trial,lin,bf079,ni79)*db &
                + b1psifd2(trial,lin,bf079,ni79)*db
           ssterm(i,psi_g) = ssterm(i,psi_g) -     thimp_bf     *dt*temp
           ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimp_bf*bdf)*dt*temp

           temp = b1bfd1  (trial,lin,bf079,ni79)*db &
                + b1bfd2  (trial,lin,bf079,ni79)*db
           ssterm(i,bz_g) = ssterm(i,bz_g) -     thimp_bf     *dt*temp
           ddterm(i,bz_g) = ddterm(i,bz_g) + (1.-thimp_bf*bdf)*dt*temp

           temp = b1psifd1(trial,ps079,lin,ni79)*db &
                + b1psifd2(trial,ps079,lin,ni79)*db &
                + b1bfd1  (trial,bz079,lin,ni79)*db &
                + b1bfd2  (trial,bz079,lin,ni79)*db
           r_bf(i) = r_bf(i) -     thimp_bf     *dt*temp
           q_bf(i) = q_bf(i) + (1.-thimp_bf*bdf)*dt*temp
        end if
        if(use_external_fields .and. linear.eq.0) then
           temp = b1psifd1(trial,lin,bfx79,ni79)*db
           ssterm(i,psi_g) = ssterm(i,psi_g) -     thimp_bf     *dt*temp
           ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimp_bf*bdf)*dt*temp

           temp = b1bfd1  (trial,lin,bfx79,ni79)*db
           ssterm(i,bz_g) = ssterm(i,bz_g) -     thimp_bf     *dt*temp
           ddterm(i,bz_g) = ddterm(i,bz_g) + (1.-thimp_bf*bdf)*dt*temp

           temp = b1psifd2(trial,psx79,lin,ni79)*db &
                + b1bfd2  (trial,bzx79,lin,ni79)*db
           r_bf(i) = r_bf(i) -     thimp_bf     *dt*temp
           q_bf(i) = q_bf(i) + (1.-thimp_bf*bdf)*dt*temp
        end if
     end if

     if(idens.eq.1 .and. (eqsubtract.eq.1 .or. icsubtract.eq.1)) then
        q_ni(i,1) = q_ni(i,1) + db*dt* &
             (b1psipsid (trial,ps079,ps079,lin) &
             +b1psibd1  (trial,ps079,bz079,lin) &
             +b1psibd2  (trial,ps079,bz079,lin) &
             +b1bbd     (trial,bz079,bz079,lin))
     endif
  end if

  ! Grad(p_e)
  ! ~~~~~~~~~
  if(db.ne.0.) then
     if(numvar.ge.3) then
        temp = b1ped(trial,lin,ni79)*db
        ssterm(i,pe_g) = ssterm(i,pe_g) -     thimpf     *dt*temp
        ddterm(i,pe_g) = ddterm(i,pe_g) + (1.-thimpf*bdf)*dt*temp
     end if

     if(idens.eq.1 .and. eqsubtract.eq.1) then
        q_ni(i,1) = q_ni(i,1) + db*dt* &
             (b1ped(trial,pe079,lin))
     endif
  end if  

  ! Harned-Mikic Term
  ! ~~~~~~~~~~~~~~~~~
  call b1harnedmikic(trial,lin,temp,temp2)
  ssterm(i,psi_g) = ssterm(i,psi_g) - thimpf**2*dt*dt*temp
  ddterm(i,psi_g) = ddterm(i,psi_g) - thimpf**2*dt*dt*temp
  if(numvar.ge.2) then 
     ssterm(i,bz_g) = ssterm(i,bz_g) - thimpf**2*dt*dt*temp2
     ddterm(i,bz_g) = ddterm(i,bz_g) - thimpf**2*dt*dt*temp2
  endif

  case(2)  ! ion form
!

  !  velocity time derivatives
  !  ~~~~~~~~~~~~~~~~~~~~~~~~
  if(db.ne.0.) then
     temp = b1vzdot(trial,lin)
     ssterm(i,vz_g) = ssterm(i,vz_g) + db*temp
     ddterm(i,vz_g) = ddterm(i,vz_g) + db*temp*bdf

     if(numvar.ge.3) then
        temp = b1chidot(trial,lin)
        ssterm(i,chi_g) = ssterm(i,chi_g) + db*temp
        ddterm(i,chi_g) = ddterm(i,chi_g) + db*temp*bdf
     endif
 
  endif

  ! -Grad(p_i)
  ! ~~~~~~~~~
  if(db.ne.0.) then
     if(numvar.ge.3) then
        temp = b1ped(trial,lin,ni79)*db
        ssterm(i,p_g)  = ssterm(i,p_g)  +     thimpf     *dt*temp
        ddterm(i,p_g)  = ddterm(i,p_g)  - (1.-thimpf*bdf)*dt*temp
        ssterm(i,pe_g) = ssterm(i,pe_g) -     thimpf     *dt*temp
        ddterm(i,pe_g) = ddterm(i,pe_g) + (1.-thimpf*bdf)*dt*temp
     end if

     if(idens.eq.1 .and. eqsubtract.eq.1) then
        q_ni(i,1) = q_ni(i,1) - db*dt* &
             (b1ped(trial,p079,lin) - b1ped(trial,pe079,lin))
     endif
  end if  

  ! Convective derivative terms
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(db.ne.0) then
     temp = .5*b1vphi(trial,lin,pht79)*db
     ssterm(i,vz_g) = ssterm(i,vz_g) -     thimpf     *dt*temp
     ddterm(i,vz_g) = ddterm(i,vz_g) + (1.-thimpf*bdf)*dt*temp

     temp = .5*b1vphi(trial,vzt79,lin)*db
     ssterm(i,u_g) = ssterm(i,u_g) -     thimpf     *dt*temp
     ddterm(i,u_g) = ddterm(i,u_g) + (1.-thimpf*bdf)*dt*temp

     if(numvar.ge.3) then
        temp = .5*b1vchi(trial,lin,cht79)*db
        ssterm(i,vz_g) = ssterm(i,vz_g) -     thimpf     *dt*temp
        ddterm(i,vz_g) = ddterm(i,vz_g) + (1.-thimpf*bdf)*dt*temp

        temp = .5*b1vchi(trial,vzt79,lin)*db
        ssterm(i,chi_g) = ssterm(i,chi_g) -     thimpf     *dt*temp
        ddterm(i,chi_g) = ddterm(i,chi_g) + (1.-thimpf*bdf)*dt*temp
     endif
  endif

! Viscosity terms
! ~~~~~~~~~~~~~~~
  if(db.ne.0) then
     temp = b1vmun(trial,lin,vis79,ni79)*db
     ssterm(i,vz_g) = ssterm(i,vz_g) -     thimpf     *dt*temp
     ddterm(i,vz_g) = ddterm(i,vz_g) + (1.-thimpf*bdf)*dt*temp
  endif
!
  case(3)  !  only include parallel electron pressure gradient
  ! Grad(p_e)
  ! ~~~~~~~~~
  if(db.ne.0.) then
     if(numvar.ge.3) then
        temp = b1psi2bfpe(trial,pst79,pst79,bzt79,bft79,lin)*db
        ssterm(i,pe_g) = ssterm(i,pe_g) -     thimpf     *dt*temp
        ddterm(i,pe_g) = ddterm(i,pe_g) + (1.-thimpf*bdf)*dt*temp

        if(eqsubtract.eq.1) then
          temp = b1psi2bfpe(trial,lin,pst79,bzt79,bft79,pe079)*db
          ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
          ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimpf*bdf)*dt*temp
        endif
     end if
  end if  

  end select        

  if(ibootstrap_model.eq.1) then 
     call bootstrap_flux(trial, lin, ssterm(i,:), ddterm(i,:), r_bf(i), q_bf(i), thimpf, thimp_bf)
  end if

  end do
end subroutine flux_lin


subroutine flux_nolin(trialx, r4term)

  use math
  use basic
  use m3dc1_nint
  use metricterms_new

  implicit none

  vectype, intent(in), dimension(MAX_PTS, OP_NUM, dofs_per_element)  :: trialx
  vectype, intent(out), dimension(dofs_per_element) :: r4term

  vectype, dimension(MAX_PTS, OP_NUM)  :: trial
  integer :: i
  
  r4term = 0.

  if(igauge.eq.1 .or. linear.eq.1) then
     r4term = r4term - dt* &
          vloop*intx1(trialx(:,OP_1,:))/twopi
  endif

  if(icd_source.gt.0) then
     if(jadv.eq.0) then
        r4term = r4term - &
             dt*intx3(trialx(:,OP_1,:),eta79(:,OP_1),cd79(:,OP_1))
     else
        r4term = r4term - &
             dt*intx4(trialx(:,OP_GS,:),ri2_79,eta79(:,OP_1),cd79(:,OP_1))
     endif
  endif

  if(use_external_fields .and. (eqsubtract.eq.1 .or. icsubtract.eq.1)) then
     ! VxB
     ! ~~~
     r4term = r4term + dt* &
          (b1psiu  (trialx,psx79,ph079) &
          +b1psiv  (trialx,psx79,vz079) &
          +b1psichi(trialx,psx79,ch079) &
          +b1bu    (trialx,bzx79,ph079) &
          +b1bv    (trialx,bzx79,vz079) &
          +b1bchi  (trialx,bzx79,ch079))

     if(i3d.eq.1 .and. numvar.ge.2) then
        r4term = r4term + dt* &
             (b1fu  (trialx,bfx79,ph079) &
             +b1fv  (trialx,bfx79,vz079) &
             +b1fchi(trialx,bfx79,ch079))
     end if
  end if


  do i=1, dofs_per_element
     trial = trialx(:,:,i)

  if(use_external_fields .and. (eqsubtract.eq.1 .or. icsubtract.eq.1)) then 
     ! JxB
     ! ~~~
     if(db.ne.0.) then
        r4term(i) = r4term(i) + db*dt* &
             (b1psipsid(trial,psx79,ps079,ni79) &
             +b1psibd1 (trial,psx79,bz079,ni79) &
             +b1psibd2 (trial,ps079,bzx79,ni79) &
             +b1bbd    (trial,bz079,bzx79,ni79))

        if(i3d.eq.1 .and. numvar.ge.2) then
           r4term(i) = r4term(i) + db*dt* &
                (b1psifd1(trial,ps079,bfx79,ni79) &
                +b1psifd2(trial,psx79,bf079,ni79) &
                +b1bfd1  (trial,bz079,bfx79,ni79) &
                +b1bfd2  (trial,bzx79,bf079,ni79))
        end if
     end if
  end if
  
  end do
 
end subroutine flux_nolin


!======================================================================
! Axial Field Equation
!======================================================================
subroutine axial_field_lin(trialx, lin, ssterm, ddterm, q_ni, r_bf, q_bf, &
     izone)
  
  use basic
  use arrays
  use m3dc1_nint
  use metricterms_new
  use two_fluid
  use harned_mikic_mod
  use bootstrap

  implicit none

  vectype, dimension(MAX_PTS, OP_NUM, dofs_per_element), intent(in) :: trialx
  vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: lin
  vectype, dimension(dofs_per_element,num_fields), intent(out) :: ssterm, ddterm

  vectype, intent(out), dimension(dofs_per_element,2) :: q_ni
  vectype, intent(out), dimension(dofs_per_element) :: r_bf, q_bf
  integer, intent(in) :: izone

  vectype, dimension(MAX_PTS, OP_NUM) :: trial
  vectype, dimension(dofs_per_element) :: tempx
  integer :: i

  vectype :: temp, temp2
  real :: thimpb, thimpf, thimp_bf

  vectype :: freq_fac

  thimpf = thimp

  if(imp_mod.eq.0) then
     thimpb = thimp
  else
     thimpb = 1.
  endif

  if(imp_bf.eq.1) then
     thimp_bf = thimpf
  else
     thimp_bf = 0.
  end if

  if(itime_independent.eq.1) then
#ifdef USECOMPLEX
     freq_fac = (0,1)*frequency*dt
#else
     freq_fac = 0.
#endif
  else
     freq_fac = 1.
  end if

  ssterm = 0.
  ddterm = 0.
  q_ni = 0.
  r_bf = 0.
  q_bf = 0.

  if(iestatic.eq.1) then
     if(.not.surface_int) then
        tempx = intx2(trialx(:,OP_1,:),lin(:,OP_1))
        ssterm(:,bz_g) = tempx
        ddterm(:,bz_g) = tempx
     endif
     return
  end if

  if(numvar.lt.2) return


  ! Resistive and Hyper Terms
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~
  tempx = b2psieta(trialx,lin,eta79)
  ssterm(:,psi_g) = ssterm(:,psi_g) -     thimp     *dt*tempx
  ddterm(:,psi_g) = ddterm(:,psi_g) + (1.-thimp*bdf)*dt*tempx
  
  tempx = b2beta(trialx,lin,eta79,vz079)
  ssterm(:,bz_g) = ssterm(:,bz_g) -     thimp     *dt*tempx
  ddterm(:,bz_g) = ddterm(:,bz_g) + (1.-thimp*bdf)*dt*tempx

  if(i3d.eq.1) then
     tempx = b2feta(trialx,lin,eta79)
     r_bf = r_bf -     thimp_bf     *dt*tempx
     q_bf = q_bf + (1.-thimp_bf*bdf)*dt*tempx

     if(imp_hyper.eq.2) then
        tempx = b2fj(trialx,bft79,lin)
        ssterm(:,e_g) = ssterm(:,e_g) - dt*tempx
     endif
  end if

  if(imp_hyper.eq.2) then
     tempx = b2psij(trialx,pst79,lin)
     ssterm(:,e_g) = ssterm(:,e_g) - dt*tempx
  endif

  if(izone.eq.3) return


  ! Time Derivative
  ! ~~~~~~~~~~~~~~~
  tempx = (b2b (trialx,lin) &
       -  b2bd(trialx,lin,ni79)) &   ! electron mass term
       * freq_fac
  ssterm(:,bz_g) = ssterm(:,bz_g) + tempx
  if(itime_independent.eq.0) ddterm(:,bz_g) = ddterm(:,bz_g) + tempx*bdf

  if(izone.ne.1) return


  ! VxB
  ! ~~~
  if(linear.eq.0) then
     tempx = b2psiv(trialx,lin,vz179)
     ssterm(:,psi_g) = ssterm(:,psi_g) -     thimpb     *dt*tempx
     ddterm(:,psi_g) = ddterm(:,psi_g) + (.5-thimpb*bdf)*dt*tempx

     tempx = b2bu  (trialx,lin,ph179)  &
          + b2bchi(trialx,lin,ch179)
     ssterm(:,bz_g) = ssterm(:,bz_g) -     thimpb     *dt*tempx
     ddterm(:,bz_g) = ddterm(:,bz_g) + (.5-thimpb*bdf)*dt*tempx

     tempx = b2bu  (trialx,bz179,lin)
     ssterm(:,u_g) = ssterm(:,u_g) -     thimpb     *dt*tempx
     ddterm(:,u_g) = ddterm(:,u_g) + (.5-thimpb*bdf)*dt*tempx

     tempx = b2psiv(trialx,ps179,lin)
     ssterm(:,vz_g) = ssterm(:,vz_g) -     thimpb     *dt*tempx
     ddterm(:,vz_g) = ddterm(:,vz_g) + (.5-thimpb*bdf)*dt*tempx

     if(numvar.ge.3) then
        tempx = b2bchi(trialx,bz179,lin)                  
        ssterm(:,chi_g) = ssterm(:,chi_g)      -thimpb *dt*tempx
        ddterm(:,chi_g) = ddterm(:,chi_g) + (.5-thimpb)*dt*tempx*bdf
     end if
  end if

  if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
     tempx = b2psiv(trialx,lin,vz079)
     ssterm(:,psi_g) = ssterm(:,psi_g) -     thimpb     *dt*tempx
     ddterm(:,psi_g) = ddterm(:,psi_g) + (1.-thimpb*bdf)*dt*tempx

     tempx = b2bu  (trialx,lin,ph079) &
          + b2bchi(trialx,lin,ch079)
     ssterm(:,bz_g) = ssterm(:,bz_g) -     thimpb     *dt*tempx
     ddterm(:,bz_g) = ddterm(:,bz_g) + (1.-thimpb*bdf)*dt*tempx

     tempx = b2bu  (trialx,bz079,lin)
     ssterm(:,u_g) = ssterm(:,u_g) -     thimpb     *dt*tempx
     ddterm(:,u_g) = ddterm(:,u_g) + (1.-thimpb*bdf)*dt*tempx
    
     tempx = b2psiv(trialx,ps079,lin)
     ssterm(:,vz_g) = ssterm(:,vz_g) -     thimpb     *dt*tempx
     ddterm(:,vz_g) = ddterm(:,vz_g) + (1.-thimpb*bdf)*dt*tempx

     if(numvar.ge.3) then
        tempx = b2bchi(trialx,bz079,lin)
        ssterm(:,chi_g) = ssterm(:,chi_g) -     thimpb     *dt*tempx
        ddterm(:,chi_g) = ddterm(:,chi_g) + (1.-thimpb*bdf)*dt*tempx
     end if
  endif

  if(use_external_fields .and. linear.eq.0) then
     tempx = b2bu  (trialx,bzx79,lin)
     ssterm(:,u_g) = ssterm(:,u_g) -     thimpb     *dt*tempx
     ddterm(:,u_g) = ddterm(:,u_g) + (1.-thimpb*bdf)*dt*tempx
    
     tempx = b2psiv(trialx,psx79,lin)
     ssterm(:,vz_g) = ssterm(:,vz_g) -     thimpb     *dt*tempx
     ddterm(:,vz_g) = ddterm(:,vz_g) + (1.-thimpb*bdf)*dt*tempx

     if(numvar.ge.3) then
        tempx = b2bchi(trialx,bzx79,lin)
        ssterm(:,chi_g) = ssterm(:,chi_g) -     thimpb     *dt*tempx
        ddterm(:,chi_g) = ddterm(:,chi_g) + (1.-thimpb*bdf)*dt*tempx
     end if
  endif

  if(i3d.eq.1) then
     if(linear.eq.0) then
        tempx = b2fv(trialx,bf179,lin)
        ssterm(:,vz_g) = ssterm(:,vz_g) -     thimp_bf *dt*tempx
        ddterm(:,vz_g) = ddterm(:,vz_g) + (.5-thimp_bf)*dt*tempx*bdf

        tempx = b2fv(trialx,lin,vz179)
        r_bf = r_bf -     thimp_bf *dt*tempx
        q_bf = q_bf + (.5-thimp_bf)*dt*tempx*bdf
     end if
     if(eqsubtract.eq.1) then
        tempx = b2fv(trialx,bf079,lin)
        ssterm(:,vz_g) = ssterm(:,vz_g) -     thimp_bf *dt*tempx
        ddterm(:,vz_g) = ddterm(:,vz_g) + (1.-thimp_bf)*dt*tempx*bdf

        tempx = b2fv(trialx,lin,vz079)
        r_bf = r_bf -     thimp_bf *dt*tempx
        q_bf = q_bf + (1.-thimp_bf)*dt*tempx*bdf
     end if
     if(use_external_fields .and. linear.eq.0) then
        tempx = b2fv(trialx,bfx79,lin)
        ssterm(:,vz_g) = ssterm(:,vz_g) -     thimp_bf *dt*tempx
        ddterm(:,vz_g) = ddterm(:,vz_g) + (1.-thimp_bf)*dt*tempx*bdf
     end if
  end if



  do i=1, dofs_per_element
     trial = trialx(:,:,i)



  select case(itwofluid)

  case(-1)
  ! JxB
  ! ~~~
  if(db.ne.0.) then
     if(linear.eq.0) then
        temp = b2psipsin(trial,lin,ps179,nt79)*db &
             + b2psipsin(trial,ps179,lin,nt79)*db &
             + b2psibn  (trial,lin,bz179,nt79)*db
        ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
        ddterm(i,psi_g) = ddterm(i,psi_g) + (.5-thimpf*bdf)*dt*temp

        temp = b2psibn(trial,ps179,lin,nt79)*db &
             + b2bbn  (trial,lin,bz179,nt79)*db &
             + b2bbn  (trial,bz179,lin,nt79)*db
        ssterm(i,bz_g) = ssterm(i,bz_g) -     thimpf     *dt*temp
        ddterm(i,bz_g) = ddterm(i,bz_g) + (.5-thimpf*bdf)*dt*temp
     end if
     if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
        temp = b2psipsin(trial,lin,ps079,nt79)*db &
             + b2psipsin(trial,ps079,lin,nt79)*db &
             + b2psibn  (trial,lin,bz079,nt79)*db
        ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
        ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimpf*bdf)*dt*temp

        temp = b2psibn(trial,ps079,lin,nt79)*db &
             + b2bbn  (trial,lin,bz079,nt79)*db &
             + b2bbn  (trial,bz079,lin,nt79)*db
        ssterm(i,bz_g) = ssterm(i,bz_g) -     thimpf     *dt*temp
        ddterm(i,bz_g) = ddterm(i,bz_g) + (1.-thimpf*bdf)*dt*temp
     end if
     if(use_external_fields .and. linear.eq.0) then
        temp = b2psipsin(trial,lin,psx79,nt79)*db &
             + b2psibn  (trial,lin,bzx79,nt79)*db
        ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
        ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimpf*bdf)*dt*temp

        temp = b2bbn  (trial,lin,bzx79,nt79)*db
        ssterm(i,bz_g) = ssterm(i,bz_g) -     thimpf     *dt*temp
        ddterm(i,bz_g) = ddterm(i,bz_g) + (1.-thimpf*bdf)*dt*temp
     end if

     if(i3d.eq.1) then
        if(linear.eq.0) then
           temp = b2psifn(trial,lin,bf179,nt79)*db
           ssterm(i,psi_g) = ssterm(i,psi_g) -     thimp_bf     *dt*temp
           ddterm(i,psi_g) = ddterm(i,psi_g) + (.5-thimp_bf*bdf)*dt*temp
        
           temp = b2bfn  (trial,lin,bf179,nt79)*db
           ssterm(i,bz_g) = ssterm(i,bz_g) -     thimp_bf     *dt*temp
           ddterm(i,bz_g) = ddterm(i,bz_g) + (.5-thimp_bf*bdf)*dt*temp
        
           temp = b2psifn(trial,ps179,lin,nt79)*db &
                + b2bfn  (trial,bz179,lin,nt79)*db
           r_bf(i) = r_bf(i) -     thimp_bf     *dt*temp
           q_bf(i) = q_bf(i) + (.5-thimp_bf*bdf)*dt*temp
        end if
        if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
           temp = b2psifn(trial,lin,bf079,nt79)*db
           ssterm(i,psi_g) = ssterm(i,psi_g) -     thimp_bf     *dt*temp
           ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimp_bf*bdf)*dt*temp
        
           temp = b2bfn  (trial,lin,bf079,nt79)*db
           ssterm(i,bz_g) = ssterm(i,bz_g) -     thimp_bf     *dt*temp
           ddterm(i,bz_g) = ddterm(i,bz_g) + (1.-thimp_bf*bdf)*dt*temp
        
           temp = b2psifn(trial,ps079,lin,nt79)*db &
                + b2bfn  (trial,bz079,lin,nt79)*db
           r_bf(i) = r_bf(i) -     thimp_bf     *dt*temp
           q_bf(i) = q_bf(i) + (1.-thimp_bf*bdf)*dt*temp
        end if
        if(use_external_fields .and. linear.eq.0) then
           temp = b2psifn(trial,lin,bfx79,nt79)*db
           ssterm(i,psi_g) = ssterm(i,psi_g) -     thimp_bf     *dt*temp
           ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimp_bf*bdf)*dt*temp
        
           temp = b2bfn  (trial,bzx79,lin,nt79)*db
           r_bf(i) = r_bf(i) -     thimp_bf     *dt*temp
           q_bf(i) = q_bf(i) + (1.-thimp_bf*bdf)*dt*temp
        end if
     end if

     if(idens.eq.1) then
        temp = b2psipsin(trial,pst79,pstx79,lin)*db &
             + b2psibn  (trial,pst79,bztx79,lin)*db &
             + b2bbn    (trial,bzt79,bztx79,lin)*db
        if(i3d.eq.1 .and. numvar.ge.2) then 
           temp = temp &
                + b2psifn(trial,pst79,bftx79,lin)*db &
                + b2bfn  (trial,bztx79,bft79,lin)*db
        end if
        ssterm(i,den_g) = ssterm(i,den_g) + thimpf*dt*temp
        ddterm(i,den_g) = ddterm(i,den_g) + thimpf*dt*temp

        if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
           temp = b2psipsin(trial,ps079,ps079,lin)*db &
                + b2psibn  (trial,ps079,bz079,lin)*db &
                + b2bbn    (trial,bz079,bz079,lin)*db
           if(i3d.eq.1 .and. numvar.ge.2) then 
              temp = temp &
                   +b2psifn(trial,ps079,bf079,lin)*db &
                   +b2bfn  (trial,bz079,bf079,lin)*db
           end if
           if(linear.eq.1) temp = -temp
           ddterm(i,den_g) = ddterm(i,den_g) + dt*temp
        end if
     end if
  end if

  ! Grad(pe)
  ! ~~~~~~~~
  if(db.ne.0 .and. numvar.ge.3) then
     temp = b2pen(trial,lin,nt79)*db
     ssterm(i,pe_g) = ssterm(i,pe_g) -     thimpf     *dt*temp
     ddterm(i,pe_g) = ddterm(i,pe_g) + (1.-thimpf*bdf)*dt*temp

     if(idens.eq.1) then
        temp = b2pen(trial,pet79,lin)*db
        ssterm(i,den_g) = ssterm(i,den_g) + thimpf*dt*temp
        ddterm(i,den_g) = ddterm(i,den_g) + thimpf*dt*temp

        if(eqsubtract.eq.1) then 
           temp = b2pen(trial,pe079,lin)*db
           if(linear.eq.1) temp = -temp
           ddterm(i,den_g) = ddterm(i,den_g) + dt*temp
        end if
     end if
  end if

!  if(ibootstrap.gt.0) then
!     temp = b2psimue(trial,lin,)*db
!     ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
!     ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimpf*bdf)*dt*temp
!  endif


  case(1)  ! electron form
  ! JxB
  ! ~~~
  if(db.ne.0.) then
     if(linear.eq.0) then
        temp = b2psipsid(trial,lin,ps179,ni79)*db &
             + b2psipsid(trial,ps179,lin,ni79)*db &
             + b2psibd  (trial,lin,bz179,ni79)*db
        ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
        ddterm(i,psi_g) = ddterm(i,psi_g) + (.5-thimpf*bdf)*dt*temp

        temp = b2psibd(trial,ps179,lin,ni79)*db &
             + b2bbd  (trial,lin,bz179,ni79)*db &
             + b2bbd  (trial,bz179,lin,ni79)*db
        ssterm(i,bz_g) = ssterm(i,bz_g) -     thimpf     *dt*temp
        ddterm(i,bz_g) = ddterm(i,bz_g) + (.5-thimpf*bdf)*dt*temp
     end if
     if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
        temp = b2psipsid(trial,lin,ps079,ni79)*db &
             + b2psipsid(trial,ps079,lin,ni79)*db &
             + b2psibd  (trial,lin,bz079,ni79)*db
        ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
        ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimpf*bdf)*dt*temp

        temp = b2psibd(trial,ps079,lin,ni79)*db &
             + b2bbd  (trial,lin,bz079,ni79)*db &
             + b2bbd  (trial,bz079,lin,ni79)*db
        ssterm(i,bz_g) = ssterm(i,bz_g) -     thimpf     *dt*temp
        ddterm(i,bz_g) = ddterm(i,bz_g) + (1.-thimpf*bdf)*dt*temp
     end if
     if(use_external_fields .and. linear.eq.0) then
        temp = b2psipsid(trial,lin,psx79,ni79)*db &
             + b2psibd  (trial,lin,bzx79,ni79)*db
        ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
        ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimpf*bdf)*dt*temp

        temp = b2bbd  (trial,lin,bzx79,ni79)*db
        ssterm(i,bz_g) = ssterm(i,bz_g) -     thimpf     *dt*temp
        ddterm(i,bz_g) = ddterm(i,bz_g) + (1.-thimpf*bdf)*dt*temp
     end if

     if(i3d.eq.1) then
        if(linear.eq.0) then
           temp = b2psifd(trial,lin,bf179,ni79)*db
           ssterm(i,psi_g) = ssterm(i,psi_g) -     thimp_bf     *dt*temp
           ddterm(i,psi_g) = ddterm(i,psi_g) + (.5-thimp_bf*bdf)*dt*temp
        
           temp = b2bfd  (trial,lin,bf179,ni79)*db
           ssterm(i,bz_g) = ssterm(i,bz_g) -     thimp_bf     *dt*temp
           ddterm(i,bz_g) = ddterm(i,bz_g) + (.5-thimp_bf*bdf)*dt*temp
        
           temp = b2psifd(trial,ps179,lin,ni79)*db &
                + b2bfd  (trial,bz179,lin,ni79)*db
           r_bf(i) = r_bf(i) -     thimp_bf     *dt*temp
           q_bf(i) = q_bf(i) + (.5-thimp_bf*bdf)*dt*temp
        end if
        if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
           temp = b2psifd(trial,lin,bf079,ni79)*db
           ssterm(i,psi_g) = ssterm(i,psi_g) -     thimp_bf     *dt*temp
           ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimp_bf*bdf)*dt*temp
        
           temp = b2bfd  (trial,lin,bf079,ni79)*db
           ssterm(i,bz_g) = ssterm(i,bz_g) -     thimp_bf     *dt*temp
           ddterm(i,bz_g) = ddterm(i,bz_g) + (1.-thimp_bf*bdf)*dt*temp
        
           temp = b2psifd(trial,ps079,lin,ni79)*db &
                + b2bfd  (trial,bz079,lin,ni79)*db
           r_bf(i) = r_bf(i) -     thimp_bf     *dt*temp
           q_bf(i) = q_bf(i) + (1.-thimp_bf*bdf)*dt*temp
        end if
        if(use_external_fields .and. linear.eq.0) then
           temp = b2psifd(trial,lin,bfx79,ni79)*db
           ssterm(i,psi_g) = ssterm(i,psi_g) -     thimp_bf     *dt*temp
           ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimp_bf*bdf)*dt*temp
        
           temp = b2bfd  (trial,bzx79,lin,ni79)*db
           r_bf(i) = r_bf(i) -     thimp_bf     *dt*temp
           q_bf(i) = q_bf(i) + (1.-thimp_bf*bdf)*dt*temp
        end if
     end if
     if(idens.eq.1 .and. (eqsubtract.eq.1 .or. icsubtract.eq.1)) then
        q_ni(i,1) = q_ni(i,1) + db*dt* &
             (b2psipsid(trial,ps079,ps079,lin) &
             +b2psibd  (trial,ps079,bz079,lin) &
             +b2bbd    (trial,bz079,bz079,lin))
     endif

  end if

  ! Grad(pe)
  ! ~~~~~~~~
  if(db.ne.0 .and. numvar.ge.3) then
     temp = b2ped(trial,lin,ni79)*db
     ssterm(i,pe_g) = ssterm(i,pe_g) -     thimpf     *dt*temp
     ddterm(i,pe_g) = ddterm(i,pe_g) + (1.-thimpf*bdf)*dt*temp

     if(idens.eq.1 .and. eqsubtract.eq.1) then
        q_ni(i,1) = q_ni(i,1) + db*dt* &
             (b2ped(trial,pe079,lin))
     end if
  end if

!  if(ibootstrap.gt.0) then
!     temp = b2psimue(trial,lin,)*db
!     ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
!     ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimpf*bdf)*dt*temp
!  endif


  call b2harnedmikic(trial,lin,temp,temp2)
  ssterm(i,psi_g) = ssterm(i,psi_g) - thimpf**2*dt*dt*temp
  ddterm(i,psi_g) = ddterm(i,psi_g) - thimpf**2*dt*dt*temp
  ssterm(i,bz_g) = ssterm(i,bz_g) - thimpf**2*dt*dt*temp2
  ddterm(i,bz_g) = ddterm(i,bz_g) - thimpf**2*dt*dt*temp2

  case(2)  ! ion form

!

  !  velocity time derivatives
  !  ~~~~~~~~~~~~~~~~~~~~~~~~
  if(db.ne.0.) then
     temp = b2phidot(trial,lin)
     ssterm(i,u_g) = ssterm(i,u_g) + db*temp
     ddterm(i,u_g) = ddterm(i,u_g) + db*temp*bdf

     if(numvar.ge.3) then
        temp = b2chidot(trial,lin)
        ssterm(i,chi_g) = ssterm(i,chi_g) + db*temp
        ddterm(i,chi_g) = ddterm(i,chi_g) + db*temp*bdf
     endif
 
  endif

  ! -Grad(pi)
  ! ~~~~~~~~
  if(db.ne.0 .and. numvar.ge.3) then
     temp = b2ped(trial,lin,ni79)*db
     ssterm(i,p_g)  = ssterm(i,p_g)  +     thimpf     *dt*temp
     ddterm(i,p_g)  = ddterm(i,p_g)  - (1.-thimpf*bdf)*dt*temp
     ssterm(i,pe_g) = ssterm(i,pe_g) -     thimpf     *dt*temp
     ddterm(i,pe_g) = ddterm(i,pe_g) + (1.-thimpf*bdf)*dt*temp

     if(idens.eq.1 .and. eqsubtract.eq.1) then
        q_ni(i,1) = q_ni(i,1) - db*dt* &
             (b2ped(trial,p079,lin) - b2ped(trial,pe079,lin))
     end if
  endif

  ! convective derivative terms
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(db.ne.0) then
     temp = 0.5*(b2uu(trial,lin,pht79) + b2uu(trial,pht79,lin))*db
     ssterm(i,u_g) = ssterm(i,u_g) - thimpf*dt*temp
     ddterm(i,u_g) = ddterm(i,u_g) + (1.-thimpf*bdf)*dt*temp
     if(numvar.ge.3) then
        temp = b2uchi(trial,lin,cht79)*db
        ssterm(i,u_g) = ssterm(i,u_g) - thimpf*dt*temp
        ddterm(i,u_g) = ddterm(i,u_g) + (1.-thimpf*bdf)*dt*temp

        temp = b2uchi(trial,pht79,lin)*db
        ssterm(i,chi_g) = ssterm(i,chi_g) - thimpf*dt*temp
        ddterm(i,chi_g) = ddterm(i,chi_g) + (1.-thimpf*bdf)*dt*temp

     endif
  endif

  ! Viscosity terms
  ! ~~~~~~~~~~~~~~~
  if(db.ne.0) then

     temp = b2phimun(trial,lin,vis79,ni79)*db
     ssterm(i,u_g) = ssterm(i,u_g) - thimpf*dt*temp
     ddterm(i,u_g) = ddterm(i,u_g) + (1.-thimpf*bdf)*dt*temp
     if(numvar.ge.3) then
        temp = b2chimun(trial,lin,vis79,ni79)*db
        ssterm(i,chi_g) = ssterm(i,chi_g) - thimpf*dt*temp
        ddterm(i,chi_g) = ddterm(i,chi_g) + (1.-thimpf*bdf)*dt*temp
     endif

  endif
!
  case(3)  !  only include parallel electron gradient
  ! Grad(p_e)
  ! ~~~~~~~~~
  if(db.ne.0.) then
     if(numvar.ge.3) then
        temp = b2psi2bfpe(trial,pst79,pst79,bzt79,bft79,lin)*db
        ssterm(i,pe_g) = ssterm(i,pe_g) -     thimpf     *dt*temp
        ddterm(i,pe_g) = ddterm(i,pe_g) + (1.-thimpf*bdf)*dt*temp
        if(eqsubtract.eq.1) then
           temp = b2psi2bfpe(trial,lin,pst79,bzt79,bft79,pe079)*db
           ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
           ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimpf*bdf)*dt*temp
        endif
     end if
  end if  
  end select

  if(ibootstrap_model.eq.1) then 
     call bootstrap_axial_field(trial, lin, ssterm(i,:), ddterm(i,:), r_bf(i), q_bf(i), thimpf, thimp_bf)
  end if

  end do

end subroutine axial_field_lin


subroutine axial_field_nolin(trialx, r4term)

  use basic
  use m3dc1_nint
  use metricterms_new

  implicit none

  vectype, intent(in), dimension(MAX_PTS, OP_NUM, dofs_per_element) :: trialx
  vectype, intent(out), dimension(dofs_per_element) :: r4term
  
  vectype, dimension(MAX_PTS, OP_NUM)  :: trial
  integer :: i

  r4term = 0.

  if(numvar.lt.2) return

  if(use_external_fields .and. (eqsubtract.eq.1 .or. icsubtract.eq.1)) then
     ! VxB
     ! ~~~
     r4term = r4term + dt* &
          (b2psiv(trialx,psx79,vz079) &
          +b2bu  (trialx,bzx79,ph079) &
          +b2bchi(trialx,bzx79,ch079))

     if(i3d.eq.1) then
        r4term = r4term + dt*b2fv(trialx,bfx79,vz079)
     end if
  end if

  do i=1, dofs_per_element
     trial = trialx(:,:,i)

  if(use_external_fields .and. (eqsubtract.eq.1 .or. icsubtract.eq.1)) then
     ! JxB
     ! ~~~
     if(db.ne.0.) then
        r4term(i) = r4term(i) + db*dt* &
             (b2psipsid(trial,ps079,psx79,ni79) &
             +b2psibd  (trial,ps079,bzx79,ni79) &
             +b2bbd    (trial,bz079,bzx79,ni79))
        
        if(i3d.eq.1) then
           r4term(i) = r4term(i) + db*dt* &
                (b2psifd(trial,ps079,bfx79,ni79) &
                +b2bfd  (trial,bzx79,bf079,ni79))
        end if
     end if
  end if

  end do

end subroutine axial_field_nolin


!======================================================================
! Pressure Equation
!======================================================================
subroutine pressure_lin(trialx, lin, ssterm, ddterm, q_ni, r_bf, q_bf,&
     total_pressure, thimpf, izone)
  
  use basic
  use math
  use arrays
  use m3dc1_nint
  use metricterms_new
  use bootstrap

  implicit none

  vectype, dimension(MAX_PTS, OP_NUM, dofs_per_element), intent(in) :: trialx
  vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: lin
  vectype, dimension(dofs_per_element,num_fields), intent(out) :: ssterm, ddterm
  vectype, intent(out), dimension(dofs_per_element, 2) :: q_ni
  vectype, intent(out), dimension(dofs_per_element) :: r_bf, q_bf
  logical, intent(in) :: total_pressure
  integer, intent(in) :: izone
  real, intent(in) :: thimpf

  vectype :: freq_fac

  vectype, dimension(MAX_PTS, OP_NUM) :: trial
  vectype, dimension(dofs_per_element) :: tempx
  integer :: i

  vectype, dimension(MAX_PTS, OP_NUM) :: pp079, pp179, ppt79
  vectype :: temp
  real :: thimpb, thimp_bf, nv, coefeq
  integer :: pp_g
  coefeq =  3853.*(n0_norm/1.e14)*(l0_norm/100.)**2     ! 3853 = (gam-1)*3 *mu_0*e^2*[1.e20]/M_i   MKS units

  if(total_pressure) then
     pp079 = p079
     pp179 = p179
     ppt79 = pt79
     pp_g = p_g
  else
     pp079 = pe079
     pp179 = pe179
     ppt79 = pet79
     pp_g = pe_g
  end if

  if(imp_mod.eq.0) then
     thimpb = thimp
  else
     thimpb = 1.
  endif

  if(imp_bf.eq.1) then
     thimp_bf = thimpf
  else
     thimp_bf = 0.
  end if

  if(numvar.eq.1) then
     nv = .5
  else
     nv = 1./3.
  end if

  if(itime_independent.eq.1) then
#ifdef USECOMPLEX
     freq_fac = (0,1)*frequency*dt
#else
     freq_fac = 0.
#endif
  else
     freq_fac = 1.
  end if

  ssterm = 0.
  ddterm = 0.
  q_ni = 0.
  r_bf = 0.
  q_bf = 0.

  if(izone.ne.1) then
     tempx = b3pe(trialx,lin)
     ssterm(:,pp_g) = ssterm(:,pp_g) + tempx
     ddterm(:,pp_g) = ddterm(:,pp_g) + tempx*bdf
     return
  end if

  ! Time Derivative
  ! ~~~~~~~~~~~~~~~
  tempx = b3pe(trialx,lin)*freq_fac
  ssterm(:,pp_g) = ssterm(:,pp_g) + tempx
  if(itime_independent.eq.0) ddterm(:,pp_g) = ddterm(:,pp_g) + tempx*bdf

  ! special to peg pressure for itaylor=27
  if(iheat_sink.eq.1 .and. itaylor.eq.27) then
      tempx = b3pe27(trialx,lin)
      ssterm(:,pp_g) = ssterm(:,pp_g) +     thimpb     *dt*(gam-1.)*coolrate*tempx
      ddterm(:,pp_g) = ddterm(:,pp_g) - (1.-thimpb*bdf)*dt*(gam-1.)*coolrate*tempx
  end if


  ! Perpendicular Heat Flux
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  tempx = b3pedkappa(trialx,lin,ni79,kap79,vz079)
  ssterm(:,pp_g) = ssterm(:,pp_g) -     thimp     *dt*tempx
  ddterm(:,pp_g) = ddterm(:,pp_g) + (1.-thimp*bdf)*dt*tempx
  if(eqsubtract.eq.1) then
     if(idens.eq.1) then
        q_ni(:,1) = q_ni(:,1) + dt* &
             (b3pedkappa(trialx,pp079,lin,kap79,vz079))
     end if
  end if


  ! Equipartition
  ! ~~~~~~~~~~~~~
!...this presently only works for ipressplit=1

 if((.not. total_pressure) .and. ipres.eq.1 .and. ipressplit.eq.1) then
    if(linear.eq.0 .and. eqsubtract.eq.0) then
       tempx = dt*coefeq*b3peeta(trialx,lin,eta79)
       ssterm(:,pe_g) = ssterm(:,pe_g) +  2.*tempx        !*thimpf
!       ddterm(:,pp_g) = ddterm(:,pp_g) +  2.*(1.-thimpf)*tempx
       ssterm(:,p_g)  = ssterm(:,p_g)  -     tempx        !*thimpf
!       ddterm(:,p_g)  = ddterm(:,p_g)  -  (1.-thimpf)*tempx
    endif
 endif


  ! Ohmic Heating
  ! ~~~~~~~~~~~~~
  if(linear.eq.0) then
     tempx = b3psipsieta(trialx,lin,ps179,eta79) &
          + b3psipsieta(trialx,ps179,lin,eta79)
     ssterm(:,psi_g) = ssterm(:,psi_g) -     thimpf     *dt*tempx
     ddterm(:,psi_g) = ddterm(:,psi_g) + (.5-thimpf*bdf)*dt*tempx

     if(numvar.ge.2) then
        tempx = b3bbeta(trialx,lin,bz179,eta79) &
             + b3bbeta(trialx,bz179,lin,eta79)
        ssterm(:,bz_g) = ssterm(:,bz_g) -     thimpf     *dt*tempx
        ddterm(:,bz_g) = ddterm(:,bz_g) + (.5-thimpf*bdf)*dt*tempx

        tempx = b3psibeta(trialx,lin,bz179,eta79) 
        ssterm(:,psi_g) = ssterm(:,psi_g) -     thimpf     *dt*tempx
        ddterm(:,psi_g) = ddterm(:,psi_g) + (.5-thimpf*bdf)*dt*tempx

        tempx = b3psibeta(trialx,ps179,lin,eta79)
        ssterm(:,bz_g) = ssterm(:,bz_g) -     thimpf     *dt*tempx
        ddterm(:,bz_g) = ddterm(:,bz_g) + (.5-thimpf*bdf)*dt*tempx
        if(i3d .eq. 1) then
           tempx = b3psifeta(trialx,lin,bf179,eta79)
           ssterm(:,psi_g) = ssterm(:,psi_g) -     thimpf     *dt*tempx
           ddterm(:,psi_g) = ddterm(:,psi_g) + (.5-thimpf*bdf)*dt*tempx

           tempx = b3psifeta(trialx,ps179,lin,eta79)
           r_bf = r_bf -     thimp_bf     *dt*tempx
           q_bf = q_bf + (.5-thimp_bf*bdf)*dt*tempx

           tempx = b3bfeta(trialx,lin,bf179,eta79)
           ssterm(:,bz_g) = ssterm(:,bz_g) -     thimpf     *dt*tempx
           ddterm(:,bz_g) = ddterm(:,bz_g) + (.5-thimpf*bdf)*dt*tempx

           tempx = b3bfeta(trialx,bz179,lin,eta79)
           r_bf = r_bf -     thimp_bf     *dt*tempx
           q_bf = q_bf + (.5-thimp_bf*bdf)*dt*tempx

           tempx = b3ffeta(trialx,lin,bf179,eta79)      &
                + b3ffeta(trialx,bf179,lin,eta79)
           r_bf = r_bf -     thimp_bf     *dt*tempx
           q_bf = q_bf + (.5-thimp_bf*bdf)*dt*tempx
        endif

     end if
  end if
  if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
     tempx = b3psipsieta(trialx,lin,ps079,eta79) &
          + b3psipsieta(trialx,ps079,lin,eta79)
     ssterm(:,psi_g) = ssterm(:,psi_g) -     thimpf     *dt*tempx
     ddterm(:,psi_g) = ddterm(:,psi_g) + (1.-thimpf*bdf)*dt*tempx
     
     if(numvar.ge.2) then
        tempx = b3bbeta(trialx,lin,bz079,eta79) &
             + b3bbeta(trialx,bz079,lin,eta79)
        ssterm(:,bz_g) = ssterm(:,bz_g) -     thimpf     *dt*tempx
        ddterm(:,bz_g) = ddterm(:,bz_g) + (1.-thimpf*bdf)*dt*tempx

        tempx = b3psibeta(trialx,lin,bz079,eta79)
        ssterm(:,psi_g) = ssterm(:,psi_g) -     thimpf     *dt*tempx
        ddterm(:,psi_g) = ddterm(:,psi_g) + (1.-thimpf*bdf)*dt*tempx

        if(i3d .eq. 1) then
           tempx = b3bfeta(trialx,bz079,lin,eta79)
           r_bf = r_bf -     thimp_bf     *dt*tempx
           q_bf = q_bf + (1.-thimp_bf*bdf)*dt*tempx
        endif
     end if
  endif


  ! Pressure Advection
  ! ~~~~~~~~~~~~~~~~~~
  if(linear.eq.0) then
     tempx = p1pu(trialx,pp179,lin)
     ssterm(:,u_g) = ssterm(:,u_g) -     thimpb     *dt*tempx
     ddterm(:,u_g) = ddterm(:,u_g) + (.5-thimpb*bdf)*dt*tempx

     if(numvar.ge.2) then
        tempx = p1pv(trialx,pp179,lin)
        ssterm(:,vz_g) = ssterm(:,vz_g) -     thimpb     *dt*tempx
        ddterm(:,vz_g) = ddterm(:,vz_g) + (.5-thimpb*bdf)*dt*tempx
     end if

     if(numvar.ge.3) then
        tempx = p1pchi(trialx,pp179,lin)
        ssterm(:,chi_g) = ssterm(:,chi_g) -     thimpb     *dt*tempx
        ddterm(:,chi_g) = ddterm(:,chi_g) + (.5-thimpb*bdf)*dt*tempx
     end if

     tempx = p1pu  (trialx,lin,ph179) & 
          + p1pv  (trialx,lin,vz179) &
          + p1pchi(trialx,lin,ch179)
     ssterm(:,pp_g) = ssterm(:,pp_g) -     thimp     *dt*tempx
     ddterm(:,pp_g) = ddterm(:,pp_g) + (.5-thimp*bdf)*dt*tempx
  endif

  if(eqsubtract.eq.1) then
     if(kinetic.eq.0) then
        tempx = p1pu(trialx,pp079,lin)
        ssterm(:,u_g) = ssterm(:,u_g) -     thimpb     *dt*tempx
        ddterm(:,u_g) = ddterm(:,u_g) + (1.-thimpb*bdf)*dt*tempx
     
        if(numvar.ge.2) then
           tempx = p1pv(trialx,pp079,lin)
           ssterm(:,vz_g) = ssterm(:,vz_g) -     thimpb     *dt*tempx
           ddterm(:,vz_g) = ddterm(:,vz_g) + (1.-thimpb*bdf)*dt*tempx
        end if
     
        if(numvar.ge.3) then
           tempx = p1pchi(trialx,pp079,lin)                
           ssterm(:,chi_g) = ssterm(:,chi_g) -     thimpb     *dt*tempx
           ddterm(:,chi_g) = ddterm(:,chi_g) + (1.-thimpb*bdf)*dt*tempx
        end if

        tempx = p1pu  (trialx,lin,ph079) & 
             + p1pv  (trialx,lin,vz079) &
             + p1pchi(trialx,lin,ch079)
        ssterm(:,pp_g) = ssterm(:,pp_g) -     thimp     *dt*tempx
        ddterm(:,pp_g) = ddterm(:,pp_g) + (1.-thimp*bdf)*dt*tempx

     else   ! on kinetic.eq.0
        if(total_pressure) then
           tempx = pperpu(trialx,ppt79,lin)  &
                + pperpupsipsib2(trialx,ppt79,lin,pstx79,pstx79,b2i79)  &
                + pperpubbb2(trialx,ppt79,lin,bztx79,bztx79,b2i79)
           ssterm(:,u_g) = ssterm(:,u_g) -     thimpb     *dt*tempx
           ddterm(:,u_g) = ddterm(:,u_g) + (1.-thimpb*bdf)*dt*tempx
     
           if(numvar.ge.2) then
              tempx = pperpvpsibb2(trialx,ppt79,lin,pstx79,bztx79,b2i79)  &
                   + pperpvbbb2(trialx,ppt79,lin,bztx79,bztx79,b2i79)
              ssterm(:,vz_g) = ssterm(:,vz_g) -     thimpb     *dt*tempx
              ddterm(:,vz_g) = ddterm(:,vz_g) + (1.-thimpb*bdf)*dt*tempx
           end if
     
           if(numvar.ge.3) then
              tempx = pperpchi(trialx,ppt79,lin)  &
                   + pperpchipsipsib2(trialx,ppt79,lin,pstx79,pstx79,b2i79)  &
                   + pperpchipsibb2(trialx,ppt79,lin,pstx79,bztx79,b2i79)  &
                   + pperpchibbb2(trialx,ppt79,lin,bztx79,bztx79,b2i79)
              ssterm(:,chi_g) = ssterm(:,chi_g) -     thimpb     *dt*tempx
              ddterm(:,chi_g) = ddterm(:,chi_g) + (1.-thimpb*bdf)*dt*tempx
           end if
        else   !on total_pressure
           tempx = pparpu(trialx,ppt79,lin)  &
                + pparpupsipsib2(trialx,ppt79,lin,pstx79,pstx79,b2i79)  &
                + pparpubbb2(trialx,ppt79,lin,bztx79,bztx79,b2i79)
           ssterm(:,u_g) = ssterm(:,u_g) -     thimpb     *dt*tempx
           ddterm(:,u_g) = ddterm(:,u_g) + (1.-thimpb*bdf)*dt*tempx
     
           if(numvar.ge.2) then
              tempx = pparpvpsibb2(trialx,ppt79,lin,pstx79,bztx79,b2i79)  &
                   + pparpvbbb2(trialx,ppt79,lin,bztx79,bztx79,b2i79)
              ssterm(:,vz_g) = ssterm(:,vz_g) -     thimpb     *dt*tempx
              ddterm(:,vz_g) = ddterm(:,vz_g) + (1.-thimpb*bdf)*dt*tempx
           end if
     
           if(numvar.ge.3) then
              tempx = pparpchi(trialx,ppt79,lin)  &
                   + pparpchipsipsib2(trialx,ppt79,lin,pstx79,pstx79,b2i79)  &
                   + pparpchipsibb2(trialx,ppt79,lin,pstx79,bztx79,b2i79)  &
                   + pparpchibbb2(trialx,ppt79,lin,bztx79,bztx79,b2i79)
              ssterm(:,chi_g) = ssterm(:,chi_g) -     thimpb     *dt*tempx
              ddterm(:,chi_g) = ddterm(:,chi_g) + (1.-thimpb*bdf)*dt*tempx
           end if
        endif
     endif
  endif



  do i=1, dofs_per_element
     trial = trialx(:,:,i)


  
  ! Electron Pressure Advection
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(db.ne.0.) then
     if(linear.eq.0) then
        temp = b3pepsid(trial,pe179,lin,ni79)*db
        ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
        ddterm(i,psi_g) = ddterm(i,psi_g) + (.5-thimpf*bdf)*dt*temp
        
        if(numvar.ge.2) then
           temp = b3pebd(trial,pe179,lin,ni79)*db
           ssterm(i,bz_g) = ssterm(i,bz_g) -     thimpf     *dt*temp
           ddterm(i,bz_g) = ddterm(i,bz_g) + (.5-thimpf*bdf)*dt*temp
        end if

        if(numvar.ge.3) then
           temp = b3pepsid(trial,lin,ps179,ni79)*db & 
                + b3pebd  (trial,lin,bz179,ni79)*db
           ssterm(i,pe_g) = ssterm(i,pe_g) -     thimpf     *dt*temp
           ddterm(i,pe_g) = ddterm(i,pe_g) + (.5-thimpf*bdf)*dt*temp
        end if
     end if
     if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
        temp = b3pepsid(trial,pe079,lin,ni79)*db
        ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
        ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimpf*bdf)*dt*temp

        if(numvar.ge.2) then
           temp = b3pebd(trial,pe079,lin,ni79)*db
           ssterm(i,bz_g) = ssterm(i,bz_g) -     thimpf     *dt*temp
           ddterm(i,bz_g) = ddterm(i,bz_g) + (1.-thimpf*bdf)*dt*temp
        end if

        if(numvar.ge.3) then
           temp = b3pepsid(trial,lin,ps079,ni79)*db &
                + b3pebd  (trial,lin,bz079,ni79)*db
           ssterm(i,pe_g) = ssterm(i,pe_g) -     thimpf     *dt*temp
           ddterm(i,pe_g) = ddterm(i,pe_g) + (1.-thimpf*bdf)*dt*temp
        end if
     end if
     if(i3d.eq.1 .and. numvar.ge.2) then
        if(linear.eq.0) then
           if(numvar.ge.3) then
              temp = b3pefd(trial,lin,bf179,ni79)*db
              ssterm(i,pe_g) = ssterm(i,pe_g) -     thimp_bf     *dt*temp
              ddterm(i,pe_g) = ddterm(i,pe_g) + (.5-thimp_bf*bdf)*dt*temp
           end if

           temp = b3pefd(trial,pe179,lin,ni79)*db
           r_bf(i) = r_bf(i) -     thimp_bf     *dt*temp
           q_bf(i) = q_bf(i) + (.5-thimp_bf*bdf)*dt*temp
        end if
        if(eqsubtract.eq.1) then
           if(numvar.ge.3) then
              temp = b3pefd(trial,lin,bf079,ni79)*db
              ssterm(i,pe_g) = ssterm(i,pe_g) -     thimp_bf     *dt*temp
              ddterm(i,pe_g) = ddterm(i,pe_g) + (1.-thimp_bf*bdf)*dt*temp
           end if

           temp = b3pefd(trial,pe079,lin,ni79)*db
           r_bf(i) = r_bf(i) -     thimp_bf     *dt*temp
           q_bf(i) = q_bf(i) + (1.-thimp_bf*bdf)*dt*temp
        end if
     end if
     if(idens.eq.1 .and. (eqsubtract.eq.1 .or. icsubtract.eq.1)) then
        q_ni(i,1) = q_ni(i,1) + db*dt* &
             (b3pepsid(trial,pe079,ps079,lin) &
             +b3pebd  (trial,pe079,bz079,lin))
     end if
  end if


  ! Source Cooling
  ! ~~~~~~~~~~~~~~
  if(total_pressure) then
     if(linear.eq.0) then
        temp = p1uus  (trial,lin,ph179,sig79) &
             + p1uus  (trial,ph179,lin,sig79) &
             + p1uchis(trial,lin,ch179,sig79) 
        ssterm(i,u_g) = ssterm(i,u_g) -     thimpb     *dt*temp
        ddterm(i,u_g) = ddterm(i,u_g) + (.5-thimpb*bdf)*dt*temp
     
        if(numvar.ge.2) then
           temp = p1vvs  (trial,lin,vz179,sig79) &
                + p1vvs  (trial,vz179,lin,sig79)
           ssterm(i,vz_g) = ssterm(i,vz_g) -     thimpb     *dt*temp
           ddterm(i,vz_g) = ddterm(i,vz_g) + (.5-thimpb*bdf)*dt*temp
        end if
        
        if(numvar.ge.3) then
           temp = p1chichis(trial,lin,ch179,sig79) &
                + p1chichis(trial,ch179,lin,sig79) &
                + p1uchis  (trial,ph179,lin,sig79) 
           ssterm(i,chi_g) = ssterm(i,chi_g) -     thimpb     *dt*temp
           ddterm(i,chi_g) = ddterm(i,chi_g) + (.5-thimpb*bdf)*dt*temp
        end if
     end if
     if(eqsubtract.eq.1) then
        temp = p1uus  (trial,lin,ph079,sig79) &
             + p1uus  (trial,ph079,lin,sig79) &
             + p1uchis(trial,lin,ch079,sig79) 
        ssterm(i,u_g) = ssterm(i,u_g) -     thimpb     *dt*temp
        ddterm(i,u_g) = ddterm(i,u_g) + (1.-thimpb*bdf)*dt*temp

        if(numvar.ge.2) then
           temp = p1vvs  (trial,lin,vz079,sig79) &
                + p1vvs  (trial,vz079,lin,sig79)
           ssterm(i,vz_g) = ssterm(i,vz_g) -     thimpb     *dt*temp
           ddterm(i,vz_g) = ddterm(i,vz_g) + (1.-thimpb*bdf)*dt*temp
        end if
        
        if(numvar.ge.3) then
           temp = p1chichis(trial,lin,ch079,sig79) &
                + p1chichis(trial,ch079,lin,sig79) &
                + p1uchis  (trial,ph079,lin,sig79) 
           ssterm(i,chi_g) = ssterm(i,chi_g) -     thimpb     *dt*temp
           ddterm(i,chi_g) = ddterm(i,chi_g) + (1.-thimpb*bdf)*dt*temp
        end if
     end if
  end if


  ! Gradient-dependent heat flux
  if(kappag.ne.0) then
     if(linear.eq.0) then
        temp = b3pedkappag(trial,lin,pp179,pp179) &
             + b3pedkappag(trial,pp179,lin,pp179) &
             + b3pedkappag(trial,pp179,pp179,lin)
        ssterm(i,pp_g) = ssterm(i,pp_g) -        thimp     *dt*temp
        ddterm(i,pp_g) = ddterm(i,pp_g) + (1./3.-thimp*bdf)*dt*temp
        
        if(eqsubtract.eq.1) then
           temp = b3pedkappag(trial,lin,pp179,pp079) &
                + b3pedkappag(trial,pp179,lin,pp079) &
                + b3pedkappag(trial,pp179,pp079,lin) &
                + b3pedkappag(trial,lin,pp079,pp179) &
                + b3pedkappag(trial,pp079,lin,pp179) &
                + b3pedkappag(trial,pp079,pp179,lin)
           ssterm(i,pp_g) = ssterm(i,pp_g) -        thimp     *dt*temp
           ddterm(i,pp_g) = ddterm(i,pp_g) + (1./2.-thimp*bdf)*dt*temp     
        end if
     end if

     if(eqsubtract.eq.1) then
        temp = b3pedkappag(trial,lin,pp079,pp079) &
             + b3pedkappag(trial,pp079,lin,pp079) &
             + b3pedkappag(trial,pp079,pp079,lin)
        ssterm(i,pp_g) = ssterm(i,pp_g) -     thimp     *dt*temp
        ddterm(i,pp_g) = ddterm(i,pp_g) + (1.-thimp*bdf)*dt*temp     
     end if
  end if

  ! Parallel Heat Flux
  ! ~~~~~~~~~~~~~~~~~~
  if(kappar.ne.0.) then
        
    if(linear.eq.0) then
!
!....note:  for ikappar_ni .eq.0, the metric calls have been changed from p1.... to te.... 
!    and ,ni79 has been removed from the calls due to numerical instability 
!    for linear.eq.0  (1/25/2012 scj)

       if(ikappar_ni.eq.0) then
          temp = tepsipsikappar(trial,pstx79,pstx79,lin,b2i79,kar79) &
               + tepsibkappar  (trial,pstx79,bztx79,lin,b2i79,kar79) &
               + tebbkappar    (trial,bztx79,bztx79,lin,b2i79,kar79)
          ssterm(i,pp_g) = ssterm(i,pp_g) -          thimp     *dt*temp
          ddterm(i,pp_g) = ddterm(i,pp_g) + (1.    - thimp*bdf)*dt*temp
          if(i3d.eq.1 .and. numvar.ge.2) then
             temp = tepsifkappar(trial,pstx79,bftx79,lin,b2i79,kar79) &
                  + tebfkappar  (trial,bztx79,bftx79,lin,b2i79,kar79) &
                  + teffkappar  (trial,bftx79,bftx79,lin,b2i79,kar79)
             ssterm(i,pp_g) = ssterm(i,pp_g) -          thimp     *dt*temp
             ddterm(i,pp_g) = ddterm(i,pp_g) + (1.    - thimp*bdf)*dt*temp
          endif

       if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
          temp = tepsipsikappar(trial,lin,ps179,pp079,b2i79,kar79) &
               + tepsipsikappar(trial,ps179,lin,pp079,b2i79,kar79) &
               + tepsibkappar  (trial,lin,bz179,pp079,b2i79,kar79)
          ssterm(i,psi_g) = ssterm(i,psi_g) -          thimp     *dt*temp
          ddterm(i,psi_g) = ddterm(i,psi_g) + (1. - .5*thimp*bdf)*dt*temp

          temp = tepsibkappar  (trial,ps179,lin,pp079,b2i79,kar79) &
               + tebbkappar    (trial,lin,bz179,pp079,b2i79,kar79) &
               + tebbkappar    (trial,bz179,lin,pp079,b2i79,kar79)
          ssterm(i,bz_g) = ssterm(i,bz_g) -          thimp     *dt*temp
          ddterm(i,bz_g) = ddterm(i,bz_g) + (1. - .5*thimp*bdf)*dt*temp

          temp = tepsipsikappar(trial,lin,ps079,pp079,b2i79,kar79) &
               + tepsipsikappar(trial,ps079,lin,pp079,b2i79,kar79) &
               + tepsibkappar  (trial,lin,bz079,pp079,b2i79,kar79)
          ssterm(i,psi_g) = ssterm(i,psi_g) -       thimp     *dt*temp
          ddterm(i,psi_g) = ddterm(i,psi_g) + (1. - thimp*bdf)*dt*temp

          temp = tepsibkappar  (trial,ps079,lin,pp079,b2i79,kar79) &
               + tebbkappar    (trial,lin,bz079,pp079,b2i79,kar79) &
               + tebbkappar    (trial,bz079,lin,pp079,b2i79,kar79)
          ssterm(i,bz_g) = ssterm(i,bz_g) -       thimp     *dt*temp
          ddterm(i,bz_g) = ddterm(i,bz_g) + (1. - thimp*bdf)*dt*temp
        
          if(i3d.eq.1 .and. numvar.ge.2) then
             temp = tepsifkappar(trial,lin,bf179,pp079,b2i79,kar79)
             ssterm(i,psi_g) = ssterm(i,psi_g) -          thimp     *dt*temp
             ddterm(i,psi_g) = ddterm(i,psi_g) + (1. - .5*thimp*bdf)*dt*temp

             temp = tebfkappar(trial,lin,bf179,pp079,b2i79,kar79)
             ssterm(i,bz_g) = ssterm(i,bz_g) -          thimp     *dt*temp
             ddterm(i,bz_g) = ddterm(i,bz_g) + (1. - .5*thimp*bdf)*dt*temp

             temp = teffkappar(trial,lin,bf179,pp079,b2i79,kar79) &
                  + teffkappar(trial,bf179,lin,pp079,b2i79,kar79)
             r_bf(i) = r_bf(i) -          thimp_bf     *dt*temp
             q_bf(i) = q_bf(i) + (1. - .5*thimp_bf*bdf)*dt*temp

             temp = tepsifkappar(trial,lin,bf079,pp079,b2i79,kar79)
             ssterm(i,psi_g) = ssterm(i,psi_g) -       thimp     *dt*temp
             ddterm(i,psi_g) = ddterm(i,psi_g) + (1. - thimp*bdf)*dt*temp
             
             temp = tebfkappar(trial,lin,bf079,pp079,b2i79,kar79)
             ssterm(i,bz_g) = ssterm(i,bz_g) -       thimp     *dt*temp
             ddterm(i,bz_g) = ddterm(i,bz_g) + (1. - thimp*bdf)*dt*temp

             temp = teffkappar(trial,lin,bf079,pp079,b2i79,kar79) &
                  + teffkappar(trial,bf079,lin,pp079,b2i79,kar79)
             r_bf(i) = r_bf(i) -       thimp_bf     *dt*temp
             q_bf(i) = q_bf(i) + (1. - thimp_bf*bdf)*dt*temp
          end if
       endif

       if(use_external_fields) then
          temp = tepsipsikappar(trial,lin,psx79,pp079,b2i79,kar79) &
               + tepsipsikappar(trial,psx79,lin,pp079,b2i79,kar79) &
               + tepsibkappar  (trial,lin,bzx79,pp079,b2i79,kar79)
          ssterm(i,psi_g) = ssterm(i,psi_g) -       thimp     *dt*temp
          ddterm(i,psi_g) = ddterm(i,psi_g) + (1. - thimp*bdf)*dt*temp

          temp = tepsibkappar  (trial,psx79,lin,pp079,b2i79,kar79) &
               + tebbkappar    (trial,lin,bzx79,pp079,b2i79,kar79) &
               + tebbkappar    (trial,bzx79,lin,pp079,b2i79,kar79)
          ssterm(i,bz_g) = ssterm(i,bz_g) -       thimp     *dt*temp
          ddterm(i,bz_g) = ddterm(i,bz_g) + (1. - thimp*bdf)*dt*temp

          temp = tepsifkappar(trial,lin,bfx79,pp079,b2i79,kar79)
          ssterm(i,psi_g) = ssterm(i,psi_g) -       thimp     *dt*temp
          ddterm(i,psi_g) = ddterm(i,psi_g) + (1. - thimp*bdf)*dt*temp
             
          temp = tebfkappar(trial,lin,bfx79,pp079,b2i79,kar79)
          ssterm(i,bz_g) = ssterm(i,bz_g) -       thimp     *dt*temp
          ddterm(i,bz_g) = ddterm(i,bz_g) + (1. - thimp*bdf)*dt*temp

          temp = teffkappar(trial,lin,bfx79,pp079,b2i79,kar79) &
               + teffkappar(trial,bfx79,lin,pp079,b2i79,kar79)
          r_bf(i) = r_bf(i) -       thimp_bf     *dt*temp
          q_bf(i) = q_bf(i) + (1. - thimp_bf*bdf)*dt*temp
       end if

       else ! on ikappar_ni
          temp = p1psipsikappar(trial,pstx79,pstx79,lin,ni79,b2i79,kar79) &
               + p1psibkappar  (trial,pstx79,bztx79,lin,ni79,b2i79,kar79) &
               + p1bbkappar    (trial,bztx79,bztx79,lin,ni79,b2i79,kar79)
          ssterm(i,pp_g) = ssterm(i,pp_g) -       thimp     *dt*temp
          ddterm(i,pp_g) = ddterm(i,pp_g) + (1. - thimp*bdf)*dt*temp
          if(i3d.eq.1 .and. numvar.ge.2) then
             temp = p1psifkappar(trial,pstx79,bftx79,lin,ni79,b2i79,kar79) &
                  + p1bfkappar  (trial,bztx79,bftx79,lin,ni79,b2i79,kar79) &
                  + p1ffkappar  (trial,bftx79,bftx79,lin,ni79,b2i79,kar79)
             ssterm(i,pp_g) = ssterm(i,pp_g) -       thimp     *dt*temp
             ddterm(i,pp_g) = ddterm(i,pp_g) + (1. - thimp*bdf)*dt*temp
          endif

          temp = p1psipsikappar(trial,lin,pstx79,ppt79,ni79,b2i79,kar79) &
               + p1psipsikappar(trial,pstx79,lin,ppt79,ni79,b2i79,kar79) &
               + p1psibkappar  (trial,lin,bztx79,ppt79,ni79,b2i79,kar79)
          ssterm(i,psi_g) = ssterm(i,psi_g) - thimp    *dt*temp
          ddterm(i,psi_g) = ddterm(i,psi_g) - thimp*bdf*dt*temp
          if(i3d.eq.1 .and. numvar.ge.2) then
             temp = p1psifkappar(trial,lin,bftx79,ppt79,ni79,b2i79,kar79)
             ssterm(i,psi_g) = ssterm(i,psi_g) - thimp    *dt*temp
             ddterm(i,psi_g) = ddterm(i,psi_g) - thimp*bdf*dt*temp
          endif 

          temp = p1psibkappar(trial,pstx79,lin,ppt79,ni79,b2i79,kar79) &
               + p1bbkappar  (trial,lin,bztx79,ppt79,ni79,b2i79,kar79) &
               + p1bbkappar  (trial,bztx79,lin,ppt79,ni79,b2i79,kar79)
          ssterm(i,bz_g) = ssterm(i,bz_g) - thimp    *dt*temp
          ddterm(i,bz_g) = ddterm(i,bz_g) - thimp*bdf*dt*temp
          if(i3d.eq.1 .and. numvar.ge.2) then
             temp = p1bfkappar(trial,lin,bftx79,ppt79,ni79,b2i79,kar79)
             ssterm(i,bz_g) = ssterm(i,bz_g) - thimp    *dt*temp
             ddterm(i,bz_g) = ddterm(i,bz_g) - thimp*bdf*dt*temp
          endif

          if(i3d.eq.1 .and. numvar.ge.2) then
             temp = p1ffkappar(trial,lin,bftx79,ppt79,ni79,b2i79,kar79) &
                  + p1ffkappar(trial,bftx79,lin,ppt79,ni79,b2i79,kar79)
             r_bf(i) = r_bf(i) - thimp_bf    *dt*temp
             q_bf(i) = q_bf(i) - thimp_bf*bdf*dt*temp
          end if

          if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
             temp = p1psipsikappar(trial,lin,ps179,pp079,ni79,b2i79,kar79) &
                  + p1psipsikappar(trial,ps179,lin,pp079,ni79,b2i79,kar79) &
                  + p1psibkappar  (trial,lin,bz179,pp079,ni79,b2i79,kar79)
             ddterm(i,psi_g) = ddterm(i,psi_g) + .5*dt*temp

             temp = p1psibkappar  (trial,ps179,lin,pp079,ni79,b2i79,kar79) &
                  + p1bbkappar    (trial,lin,bz179,pp079,ni79,b2i79,kar79) &
                  + p1bbkappar    (trial,bz179,lin,pp079,ni79,b2i79,kar79)
             ddterm(i,bz_g) = ddterm(i,bz_g) + .5*dt*temp

             if(i3d.eq.1 .and. numvar.ge.2) then
                temp = p1psifkappar(trial,lin,bf179,pp079,ni79,b2i79,kar79)
                ddterm(i,psi_g) = ddterm(i,psi_g) + .5*dt*temp
                
                temp = p1bfkappar(trial,lin,bf179,pp079,ni79,b2i79,kar79)
                ddterm(i,bz_g) = ddterm(i,bz_g) + .5*dt*temp
                
                temp = p1psifkappar(trial,ps179,lin,pp079,ni79,b2i79,kar79) &
                     + p1bfkappar  (trial,bz179,lin,pp079,ni79,b2i79,kar79) &
                     + p1ffkappar  (trial,lin,bf179,pp079,ni79,b2i79,kar79) &
                     + p1ffkappar  (trial,bf179,lin,pp079,ni79,b2i79,kar79)
                q_bf(i) = q_bf(i) + .5*dt*temp
             end if

             temp = p1psipsikappar(trial,lin,ps079,pp079,ni79,b2i79,kar79) &
                  + p1psipsikappar(trial,ps079,lin,pp079,ni79,b2i79,kar79) &
                  + p1psibkappar  (trial,lin,bz079,pp079,ni79,b2i79,kar79)
             ddterm(i,psi_g) = ddterm(i,psi_g) + dt*temp
             
             temp = p1psibkappar  (trial,ps079,lin,pp079,ni79,b2i79,kar79) &
                  + p1bbkappar    (trial,lin,bz079,pp079,ni79,b2i79,kar79) &
                  + p1bbkappar    (trial,bz079,lin,pp079,ni79,b2i79,kar79)
             ddterm(i,bz_g) = ddterm(i,bz_g) + dt*temp
             
             if(i3d.eq.1 .and. numvar.ge.2) then               
                temp = p1psifkappar(trial,lin,bf079,pp079,ni79,b2i79,kar79)
                ddterm(i,psi_g) = ddterm(i,psi_g) + dt*temp
                
                temp = p1bfkappar(trial,lin,bf079,pp079,ni79,b2i79,kar79)
                ddterm(i,bz_g) = ddterm(i,bz_g) + dt*temp
                
                temp = p1psifkappar(trial,ps079,lin,pp079,ni79,b2i79,kar79) &
                     + p1bfkappar  (trial,bz079,lin,pp079,ni79,b2i79,kar79) &
                     + p1ffkappar  (trial,lin,bf079,pp079,ni79,b2i79,kar79) &
                     + p1ffkappar  (trial,bf079,lin,pp079,ni79,b2i79,kar79)
                q_bf(i) = q_bf(i) + dt*temp
             end if
          endif
          
          if(use_external_fields) then
             temp = p1psipsikappar(trial,lin,psx79,pp079,ni79,b2i79,kar79) &
                  + p1psipsikappar(trial,psx79,lin,pp079,ni79,b2i79,kar79) &
                  + p1psibkappar  (trial,lin,bzx79,pp079,ni79,b2i79,kar79)
             ddterm(i,psi_g) = ddterm(i,psi_g) + dt*temp
             
             temp = p1psibkappar  (trial,psx79,lin,pp079,ni79,b2i79,kar79) &
                  + p1bbkappar    (trial,lin,bzx79,pp079,ni79,b2i79,kar79) &
                  + p1bbkappar    (trial,bzx79,lin,pp079,ni79,b2i79,kar79)
             ddterm(i,bz_g) = ddterm(i,bz_g) + dt*temp
             
             if(i3d.eq.1 .and. numvar.ge.2) then               
                temp = p1psifkappar(trial,lin,bfx79,pp079,ni79,b2i79,kar79)
                ddterm(i,psi_g) = ddterm(i,psi_g) + dt*temp
                
                temp = p1bfkappar(trial,lin,bfx79,pp079,ni79,b2i79,kar79)
                ddterm(i,bz_g) = ddterm(i,bz_g) + dt*temp
                
                temp = p1psifkappar(trial,psx79,lin,pp079,ni79,b2i79,kar79) &
                     + p1bfkappar  (trial,bzx79,lin,pp079,ni79,b2i79,kar79) &
                     + p1ffkappar  (trial,lin,bfx79,pp079,ni79,b2i79,kar79) &
                     + p1ffkappar  (trial,bfx79,lin,pp079,ni79,b2i79,kar79)
                q_bf(i) = q_bf(i) + dt*temp
             end if
          end if
       endif
       
    else  ! on linear.eq.0
       ! Assumes no contribution from equilibrium f

       temp = p1psipsipnkappar(trial,lin,ps079,pp079,n079,ieq_bdotgradt) &
            + p1psipsipnkappar(trial,ps079,lin,pp079,n079,1) &
            + p1psibpnkappar  (trial,lin,bz079,pp079,n079,ieq_bdotgradt,1)
       ssterm(i,psi_g) = ssterm(i,psi_g) -       thimpf     *dt*temp
       ddterm(i,psi_g) = ddterm(i,psi_g) + (1. - thimpf*bdf)*dt*temp

       if(numvar.ge.2) then
          temp = p1psibpnkappar(trial,ps079,lin,pp079,n079,1,ieq_bdotgradt) &
               + p1bbpnkappar  (trial,lin,bz079,pp079,n079,ieq_bdotgradt) &
               + p1bbpnkappar  (trial,bz079,lin,pp079,n079,1)
          ssterm(i,bz_g) = ssterm(i,bz_g) -       thimpf     *dt*temp
          ddterm(i,bz_g) = ddterm(i,bz_g) + (1. - thimpf*bdf)*dt*temp
       end if

       temp = p1psipsipnkappar(trial,ps079,ps079,lin,n079,1) &
            + p1psibpnkappar  (trial,ps079,bz079,lin,n079,1,1) &
            + p1bbpnkappar    (trial,bz079,bz079,lin,n079,1)
       ssterm(i,pp_g) = ssterm(i,pp_g) -       thimp     *dt*temp
       ddterm(i,pp_g) = ddterm(i,pp_g) + (1. - thimp*bdf)*dt*temp

       ! this term has the opposite sign because n comes in with -1 power
       if(idens.eq.1) then
          temp = p1psipsipnkappar(trial,ps079,ps079,p079,lin,1) &
               + p1psibpnkappar  (trial,ps079,bz079,p079,lin,1,1) &
               + p1bbpnkappar    (trial,bz079,bz079,p079,lin,1)
          ssterm(i,den_g) = ssterm(i,den_g) +       thimpf     *dt*temp
          ddterm(i,den_g) = ddterm(i,den_g) - (1. - thimpf*bdf)*dt*temp
       end if

       if(i3d.eq.1 .and. numvar.ge.2) then
          temp = p1psifpnkappar(trial,ps079,lin,pp079,n079,1,ieq_bdotgradt) &
               + p1bfpnkappar  (trial,bz079,lin,pp079,n079,1,ieq_bdotgradt)
          r_bf(i) = r_bf(i) -       thimp_bf     *dt*temp
          q_bf(i) = q_bf(i) + (1. - thimp_bf*bdf)*dt*temp
       endif
    endif  ! on linear.eq.0
  endif  ! on kappar.ne.0

  ! Cross-field Heat Flux
  ! ~~~~~~~~~~~~~~~~~~~~~
  if(kappax.ne.0.) then
     temp = p1kappax(trial,pp179,lin,ni79,kax79)
     ssterm(i,bz_g) = ssterm(i,bz_g) - thimp*dt*temp
     ddterm(i,bz_g) = ddterm(i,bz_g) - thimp*dt*temp*bdf

     temp = p1kappax(trial,lin,bzt79,ni79,kax79)
     ssterm(i,pp_g) = ssterm(i,pp_g) -     thimp     *dt*temp
     ddterm(i,pp_g) = ddterm(i,pp_g) + (1.-thimp*bdf)*dt*temp
  endif

  if(ibootstrap_model.eq.1) then
     call bootstrap_pressure(trial, lin, ssterm, ddterm, pp_g, thimp)
  end if

  end do
end subroutine pressure_lin

!======================================================================
! Temperature Equation
!======================================================================
subroutine temperature_lin(trialx, lin, ssterm, ddterm, q_ni, r_bf, q_bf,&
     electron_temperature, thimpf, izone)
  
  use basic
  use math
  use arrays
  use m3dc1_nint
  use metricterms_new

  implicit none

  vectype, dimension(MAX_PTS, OP_NUM, dofs_per_element), intent(in) :: trialx
  vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: lin
  vectype, dimension(dofs_per_element, num_fields), intent(out) :: ssterm, ddterm
  vectype, intent(out), dimension(dofs_per_element,2) :: q_ni
  vectype, intent(out), dimension(dofs_per_element) :: r_bf, q_bf
  logical, intent(in) :: electron_temperature
  integer, intent(in) :: izone
  real, intent(in) :: thimpf

  vectype :: freq_fac

  vectype, dimension(MAX_PTS, OP_NUM) :: trial
  vectype, dimension(dofs_per_element) :: tempx
  integer :: i

  vectype, dimension(MAX_PTS, OP_NUM) :: pp079, pp179
  vectype :: temp
  real :: thimpb, thimp_bf, nv, coefeq, ohfac
  integer :: pp_g
  coefeq =  3853.*(n0_norm/1.e14)*(l0_norm/100.)**2     ! 3853 = (gam-1)*3 *mu_0*e^2*[1.e20]/M_i   MKS units

  if(electron_temperature) then
     pp079 = te079
     pp179 = te179
     pp_g = te_g
  else
     pp079 = ti079
     pp179 = ti179
     pp_g = ti_g
  end if

  if(imp_mod.eq.0) then
     thimpb = thimp
  else
     thimpb = 1.
  endif

  if(imp_bf.eq.1) then
     thimp_bf = thimpf
  else
     thimp_bf = 0.
  end if

  if(numvar.eq.1) then
     nv = .5
  else
     nv = 1./3.
  end if

  if(itime_independent.eq.1) then
#ifdef USECOMPLEX
     freq_fac = (0,1)*frequency*dt
#else
     freq_fac = 0.
#endif
  else
     freq_fac = 1.
  end if

  ohfac = 1.
  if(ipres.eq.0) ohfac = 0.5
!...note:  all OH terms should be multiplied by 0.5 for ipres=0   (since we assume T_e = T_i


  ssterm = 0.
  ddterm = 0.
  q_ni = 0.
  r_bf = 0.
  q_bf = 0.

  if(izone.ne.1) then
     tempx = t3tn(trialx,lin,nt79)
     ssterm(:,pp_g) = ssterm(:,pp_g) + tempx
     ddterm(:,pp_g) = ddterm(:,pp_g) + tempx*bdf
     return
  end if


  if(electron_temperature) then
  ! Ohmic Heating
  ! ~~~~~~~~~~~~~
    if(linear.eq.0) then
       tempx = b3psipsieta(trialx,lin,ps179,eta79) &
            + b3psipsieta(trialx,ps179,lin,eta79)
       ssterm(:,psi_g) = ssterm(:,psi_g) -     thimpf     *dt*tempx*ohfac
       ddterm(:,psi_g) = ddterm(:,psi_g) + (.5-thimpf*bdf)*dt*tempx*ohfac

       if(numvar.ge.2) then
          tempx = b3bbeta(trialx,lin,bz179,eta79) &
               + b3bbeta(trialx,bz179,lin,eta79)
          ssterm(:,bz_g) = ssterm(:,bz_g) -     thimpf     *dt*tempx*ohfac
          ddterm(:,bz_g) = ddterm(:,bz_g) + (.5-thimpf*bdf)*dt*tempx*ohfac

          tempx = b3psibeta(trialx,lin,bz179,eta79) 
          ssterm(:,psi_g) = ssterm(:,psi_g) -     thimpf     *dt*tempx*ohfac
          ddterm(:,psi_g) = ddterm(:,psi_g) + (.5-thimpf*bdf)*dt*tempx*ohfac

          tempx = b3psibeta(trialx,ps179,lin,eta79)
          ssterm(:,bz_g) = ssterm(:,bz_g) -     thimpf     *dt*tempx*ohfac
          ddterm(:,bz_g) = ddterm(:,bz_g) + (.5-thimpf*bdf)*dt*tempx*ohfac
          if(i3d .eq. 1) then
             tempx = b3psifeta(trialx,lin,bf179,eta79)
             ssterm(:,psi_g) = ssterm(:,psi_g) -     thimpf     *dt*tempx*ohfac
             ddterm(:,psi_g) = ddterm(:,psi_g) + (.5-thimpf*bdf)*dt*tempx*ohfac

             tempx = b3psifeta(trialx,ps179,lin,eta79)
             r_bf = r_bf -     thimp_bf     *dt*tempx*ohfac
             q_bf = q_bf + (.5-thimp_bf*bdf)*dt*tempx*ohfac

             tempx = b3bfeta(trialx,lin,bf179,eta79)
             ssterm(:,bz_g) = ssterm(:,bz_g) -     thimpf     *dt*tempx*ohfac
             ddterm(:,bz_g) = ddterm(:,bz_g) + (.5-thimpf*bdf)*dt*tempx*ohfac

             tempx = b3bfeta(trialx,bz179,lin,eta79)
             r_bf = r_bf -     thimp_bf     *dt*tempx*ohfac
             q_bf = q_bf + (.5-thimp_bf*bdf)*dt*tempx*ohfac

             tempx = b3ffeta(trialx,lin,bf179,eta79)      &
                  + b3ffeta(trialx,bf179,lin,eta79)
             r_bf = r_bf -     thimp_bf     *dt*tempx*ohfac
             q_bf = q_bf + (.5-thimp_bf*bdf)*dt*tempx*ohfac
          endif

       end if
    end if
    if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
       tempx = b3psipsieta(trialx,lin,ps079,eta79) &
            + b3psipsieta(trialx,ps079,lin,eta79)
       ssterm(:,psi_g) = ssterm(:,psi_g) -     thimpf     *dt*tempx*ohfac
       ddterm(:,psi_g) = ddterm(:,psi_g) + (1.-thimpf*bdf)*dt*tempx*ohfac
     
       if(numvar.ge.2) then
          tempx = b3bbeta(trialx,lin,bz079,eta79) &
               + b3bbeta(trialx,bz079,lin,eta79)
          ssterm(:,bz_g) = ssterm(:,bz_g) -     thimpf     *dt*tempx*ohfac
          ddterm(:,bz_g) = ddterm(:,bz_g) + (1.-thimpf*bdf)*dt*tempx*ohfac

          tempx = b3psibeta(trialx,lin,bz079,eta79)
          ssterm(:,psi_g) = ssterm(:,psi_g) -     thimpf     *dt*tempx*ohfac
          ddterm(:,psi_g) = ddterm(:,psi_g) + (1.-thimpf*bdf)*dt*tempx*ohfac

          if(i3d .eq. 1) then
             tempx = b3bfeta(trialx,bz079,lin,eta79)
             r_bf = r_bf -     thimp_bf     *dt*tempx*ohfac
             q_bf = q_bf + (1.-thimp_bf*bdf)*dt*tempx*ohfac
          endif
       end if
    endif
  endif


  do i=1, dofs_per_element
     trial = trialx(:,:,i)


  ! Time Derivative
  ! ~~~~~~~~~~~~~~~
!
! NOTE:  iadiabat=1 is correct form;   adiabat=0 is for backwards compatibility (6/2/16)
  if(iadiabat.eq.1) then
     temp = t3tn(trial,lin,nt79)*freq_fac
     ssterm(i,pp_g) = ssterm(i,pp_g) + temp
     if(itime_independent.eq.0) ddterm(i,pp_g) = ddterm(i,pp_g) + temp*bdf
  else
     temp = t3t(trial,lin)*freq_fac
     ssterm(i,pp_g) = ssterm(i,pp_g) + temp
     if(itime_independent.eq.0) ddterm(i,pp_g) = ddterm(i,pp_g) + temp*bdf
  endif

!

  ! Equipartition
  ! ~~~~~~~~~~~
   if(ipres.eq.1 .and. ipressplit.eq.1) then
      if(linear.eq.0 .and. eqsubtract.eq.0) then
         temp = dt*coefeq*t3tneta(trial,lin,nt79,eta79)
         if(electron_temperature) then
            ssterm(i,pp_g) = ssterm(i,pp_g) + temp       
            ssterm(i,ti_g) = ssterm(i,ti_g) - temp
         else
            ssterm(i,pp_g) = ssterm(i,pp_g) + temp       
            ssterm(i,te_g) = ssterm(i,te_g) - temp
         endif
     endif       
  endif

  ! Temperature Advection
  ! ~~~~~~~~~~~~~~~~~~
  if(no_vdg_T .eq. 0) then   ! debug
     if(linear.eq.0) then
        temp = t3tnu(trial,pp179,nt79,lin)
        ssterm(i,u_g) = ssterm(i,u_g) -     thimpb     *dt*temp
        ddterm(i,u_g) = ddterm(i,u_g) + (.5-thimpb*bdf)*dt*temp

        if(numvar.ge.2) then
           temp = t3tnv(trial,pp179,nt79,lin)
           ssterm(i,vz_g) = ssterm(i,vz_g) -     thimpb     *dt*temp
           ddterm(i,vz_g) = ddterm(i,vz_g) + (.5-thimpb*bdf)*dt*temp
        end if

        if(numvar.ge.3) then
           temp = t3tnchi(trial,pp179,nt79,lin)
           ssterm(i,chi_g) = ssterm(i,chi_g) -     thimpb     *dt*temp
           ddterm(i,chi_g) = ddterm(i,chi_g) + (.5-thimpb*bdf)*dt*temp
        end if

        temp = t3tnu  (trial,lin,nt79,ph179) &
             + t3tnv  (trial,lin,nt79,vz179) &
             + t3tnchi(trial,lin,nt79,ch179)
        ssterm(i,pp_g) = ssterm(i,pp_g) -     thimp     *dt*temp
        ddterm(i,pp_g) = ddterm(i,pp_g) + (.5-thimp*bdf)*dt*temp
     endif
     if(eqsubtract.eq.1) then
        temp = t3tnu(trial,pp079,nt79,lin)
        ssterm(i,u_g) = ssterm(i,u_g) -     thimpb     *dt*temp
        ddterm(i,u_g) = ddterm(i,u_g) + (1.-thimpb*bdf)*dt*temp
     
        if(numvar.ge.2) then
           temp = t3tnv(trial,pp079,nt79,lin)
           ssterm(i,vz_g) = ssterm(i,vz_g) -     thimpb     *dt*temp
           ddterm(i,vz_g) = ddterm(i,vz_g) + (1.-thimpb*bdf)*dt*temp
        end if
     
        if(numvar.ge.3) then
           temp = t3tnchi(trial,pp079,nt79,lin)
           ssterm(i,chi_g) = ssterm(i,chi_g) -     thimpb     *dt*temp
           ddterm(i,chi_g) = ddterm(i,chi_g) + (1.-thimpb*bdf)*dt*temp
        end if

        temp = t3tnu  (trial,lin,nt79,ph079) &
             + t3tnv  (trial,lin,nt79,vz079) &
             + t3tnchi(trial,lin,nt79,ch079)
        ssterm(i,pp_g) = ssterm(i,pp_g) -     thimp     *dt*temp
        ddterm(i,pp_g) = ddterm(i,pp_g) + (1.-thimp*bdf)*dt*temp
     endif
  endif  ! on no_vdg_T .eq. 0

  
  ! Electron Pressure Advection
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(db.ne.0.) then
     if(linear.eq.0) then
        temp = b3pepsid(trial,pe179,lin,ni79)*db
        ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
        ddterm(i,psi_g) = ddterm(i,psi_g) + (.5-thimpf*bdf)*dt*temp
        
        if(numvar.ge.2) then
           temp = b3pebd(trial,pe179,lin,ni79)*db
           ssterm(i,bz_g) = ssterm(i,bz_g) -     thimpf     *dt*temp
           ddterm(i,bz_g) = ddterm(i,bz_g) + (.5-thimpf*bdf)*dt*temp
        end if

        if(numvar.ge.3) then
           temp = b3pepsid(trial,lin,ps179,ni79)*db & 
                + b3pebd  (trial,lin,bz179,ni79)*db
           ssterm(i,pe_g) = ssterm(i,pe_g) -     thimpf     *dt*temp
           ddterm(i,pe_g) = ddterm(i,pe_g) + (.5-thimpf*bdf)*dt*temp
        end if
     end if
     if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
        temp = b3pepsid(trial,pe079,lin,ni79)*db
        ssterm(i,psi_g) = ssterm(i,psi_g) -     thimpf     *dt*temp
        ddterm(i,psi_g) = ddterm(i,psi_g) + (1.-thimpf*bdf)*dt*temp

        if(numvar.ge.2) then
           temp = b3pebd(trial,pe079,lin,ni79)*db
           ssterm(i,bz_g) = ssterm(i,bz_g) -     thimpf     *dt*temp
           ddterm(i,bz_g) = ddterm(i,bz_g) + (1.-thimpf*bdf)*dt*temp
        end if

        if(numvar.ge.3) then
           temp = b3pepsid(trial,lin,ps079,ni79)*db &
                + b3pebd  (trial,lin,bz079,ni79)*db
           ssterm(i,pe_g) = ssterm(i,pe_g) -     thimpf     *dt*temp
           ddterm(i,pe_g) = ddterm(i,pe_g) + (1.-thimpf*bdf)*dt*temp
        end if
     end if
     if(i3d.eq.1 .and. numvar.ge.2) then
        if(linear.eq.0) then
           if(numvar.ge.3) then
              temp = b3pefd(trial,lin,bf179,ni79)*db
              ssterm(i,pe_g) = ssterm(i,pe_g) -     thimp_bf     *dt*temp
              ddterm(i,pe_g) = ddterm(i,pe_g) + (.5-thimp_bf*bdf)*dt*temp
           end if

           temp = b3pefd(trial,pe179,lin,ni79)*db
           r_bf(i) = r_bf(i) -     thimp_bf     *dt*temp
           q_bf(i) = q_bf(i) + (.5-thimp_bf*bdf)*dt*temp
        end if
        if(eqsubtract.eq.1) then
           if(numvar.ge.3) then
              temp = b3pefd(trial,lin,bf079,ni79)*db
              ssterm(i,pe_g) = ssterm(i,pe_g) -     thimp_bf     *dt*temp
              ddterm(i,pe_g) = ddterm(i,pe_g) + (1.-thimp_bf*bdf)*dt*temp
           end if

           temp = b3pefd(trial,pe079,lin,ni79)*db
           r_bf(i) = r_bf(i) -     thimp_bf     *dt*temp
           q_bf(i) = q_bf(i) + (1.-thimp_bf*bdf)*dt*temp
        end if
     end if
     if(idens.eq.1 .and. (eqsubtract.eq.1 .or. icsubtract.eq.1)) then
        q_ni(i,1) = q_ni(i,1) + db*dt* &
             (b3pepsid(trial,pe079,ps079,lin) &
             +b3pebd  (trial,pe079,bz079,lin))
     end if
  end if




     ! Perpendicular Heat Flux
     ! ~~~~~~~~~~~~~~~~~~~~~~~
     temp = b3tekappa(trial,lin,kap79,vz079)
     ssterm(i,pp_g) = ssterm(i,pp_g) -     thimp     *dt*temp
     ddterm(i,pp_g) = ddterm(i,pp_g) + (1.-thimp*bdf)*dt*temp

  ! Parallel Heat Flux
  ! ~~~~~~~~~~~~~~~~~~
  if(kappar.ne.0.) then
        
    if(linear.eq.0) then
!
          temp = tepsipsikappar(trial,pstx79,pstx79,lin,b2i79,kar79) &
               + tepsibkappar  (trial,pstx79,bztx79,lin,b2i79,kar79) &
               + tebbkappar    (trial,bztx79,bztx79,lin,b2i79,kar79)
          ssterm(i,pp_g) = ssterm(i,pp_g) -          thimp     *dt*temp
          ddterm(i,pp_g) = ddterm(i,pp_g) + (1.    - thimp*bdf)*dt*temp
          if(i3d.eq.1 .and. numvar.ge.2) then
             temp = tepsifkappar(trial,pstx79,bftx79,lin,b2i79,kar79) &
                  + tebfkappar  (trial,bztx79,bftx79,lin,b2i79,kar79) &
                  + teffkappar  (trial,bftx79,bftx79,lin,b2i79,kar79)
             ssterm(i,pp_g) = ssterm(i,pp_g) -          thimp     *dt*temp
             ddterm(i,pp_g) = ddterm(i,pp_g) + (1.    - thimp*bdf)*dt*temp
          endif

       if(eqsubtract.eq.1 .or. icsubtract.eq.1) then
          temp = tepsipsikappar(trial,lin,ps179,pp079,b2i79,kar79) &
               + tepsipsikappar(trial,ps179,lin,pp079,b2i79,kar79) &
               + tepsibkappar  (trial,lin,bz179,pp079,b2i79,kar79)
          ssterm(i,psi_g) = ssterm(i,psi_g) -          thimpf     *dt*temp*.5
          ddterm(i,psi_g) = ddterm(i,psi_g) + (1. -    thimpf*bdf)*dt*temp*.5

          temp = tepsibkappar  (trial,ps179,lin,pp079,b2i79,kar79) &
               + tebbkappar    (trial,lin,bz179,pp079,b2i79,kar79) &
               + tebbkappar    (trial,bz179,lin,pp079,b2i79,kar79)
          ssterm(i,bz_g) = ssterm(i,bz_g) -          thimpf     *dt*temp*.5
          ddterm(i,bz_g) = ddterm(i,bz_g) + (1. -    thimpf*bdf)*dt*temp*.5

          temp = tepsipsikappar(trial,lin,ps079,pp079,b2i79,kar79) &
               + tepsipsikappar(trial,ps079,lin,pp079,b2i79,kar79) &
               + tepsibkappar  (trial,lin,bz079,pp079,b2i79,kar79)
          ssterm(i,psi_g) = ssterm(i,psi_g) -       thimpf     *dt*temp
          ddterm(i,psi_g) = ddterm(i,psi_g) + (1. - thimpf*bdf)*dt*temp

          temp = tepsibkappar  (trial,ps079,lin,pp079,b2i79,kar79) &
               + tebbkappar    (trial,lin,bz079,pp079,b2i79,kar79) &
               + tebbkappar    (trial,bz079,lin,pp079,b2i79,kar79)
          ssterm(i,bz_g) = ssterm(i,bz_g) -       thimpf     *dt*temp
          ddterm(i,bz_g) = ddterm(i,bz_g) + (1. - thimpf*bdf)*dt*temp
        
          if(i3d.eq.1 .and. numvar.ge.2) then
             temp = tepsifkappar(trial,lin,bf179,pp079,b2i79,kar79)
             ssterm(i,psi_g) = ssterm(i,psi_g) -          thimpf     *dt*temp*.5
             ddterm(i,psi_g) = ddterm(i,psi_g) + (1. -    thimpf*bdf)*dt*temp*.5

             temp = tebfkappar(trial,lin,bf179,pp079,b2i79,kar79)
             ssterm(i,bz_g) = ssterm(i,bz_g) -          thimpf     *dt*temp*.5
             ddterm(i,bz_g) = ddterm(i,bz_g) + (1. -    thimpf*bdf)*dt*temp*.5

             temp = teffkappar(trial,lin,bf179,pp079,b2i79,kar79) &
                  + teffkappar(trial,bf179,lin,pp079,b2i79,kar79) &
                + tepsifkappar(trial,ps179,lin,pp079,b2i79,kar79) &
                  + tebfkappar(trial,bz179,lin,pp079,b2i79,kar79)
             r_bf(i) = r_bf(i) -          thimp_bf     *dt*temp*.5
             q_bf(i) = q_bf(i) + (1. -    thimp_bf*bdf)*dt*temp*.5

             temp = tepsifkappar(trial,lin,bf079,pp079,b2i79,kar79)
             ssterm(i,psi_g) = ssterm(i,psi_g) -       thimpf     *dt*temp
             ddterm(i,psi_g) = ddterm(i,psi_g) + (1. - thimpf*bdf)*dt*temp
             
             temp = tebfkappar(trial,lin,bf079,pp079,b2i79,kar79)
             ssterm(i,bz_g) = ssterm(i,bz_g) -       thimpf     *dt*temp
             ddterm(i,bz_g) = ddterm(i,bz_g) + (1. - thimpf*bdf)*dt*temp

             temp = teffkappar(trial,lin,bf079,pp079,b2i79,kar79) &
                  + teffkappar(trial,bf079,lin,pp079,b2i79,kar79) &
                + tepsifkappar(trial,ps079,lin,pp079,b2i79,kar79) &
                  + tebfkappar(trial,bz079,lin,pp079,b2i79,kar79)
             r_bf(i) = r_bf(i) -       thimp_bf     *dt*temp
             q_bf(i) = q_bf(i) + (1. - thimp_bf*bdf)*dt*temp
          end if
       endif

       if(use_external_fields) then
          temp = tepsipsikappar(trial,lin,psx79,pp079,b2i79,kar79) &
               + tepsipsikappar(trial,psx79,lin,pp079,b2i79,kar79) &
               + tepsibkappar  (trial,lin,bzx79,pp079,b2i79,kar79)
          ssterm(i,psi_g) = ssterm(i,psi_g) -       thimpf     *dt*temp
          ddterm(i,psi_g) = ddterm(i,psi_g) + (1. - thimpf*bdf)*dt*temp

          temp = tepsibkappar  (trial,psx79,lin,pp079,b2i79,kar79) &
               + tebbkappar    (trial,lin,bzx79,pp079,b2i79,kar79) &
               + tebbkappar    (trial,bzx79,lin,pp079,b2i79,kar79)
          ssterm(i,bz_g) = ssterm(i,bz_g) -       thimpf     *dt*temp
          ddterm(i,bz_g) = ddterm(i,bz_g) + (1. - thimpf*bdf)*dt*temp

          temp = tepsifkappar(trial,lin,bfx79,pp079,b2i79,kar79)
          ssterm(i,psi_g) = ssterm(i,psi_g) -       thimpf     *dt*temp
          ddterm(i,psi_g) = ddterm(i,psi_g) + (1. - thimpf*bdf)*dt*temp
             
          temp = tebfkappar(trial,lin,bfx79,pp079,b2i79,kar79)
          ssterm(i,bz_g) = ssterm(i,bz_g) -       thimpf     *dt*temp
          ddterm(i,bz_g) = ddterm(i,bz_g) + (1. - thimpf*bdf)*dt*temp

          temp = teffkappar(trial,lin,bfx79,pp079,b2i79,kar79) &
               + teffkappar(trial,bfx79,lin,pp079,b2i79,kar79)
          r_bf(i) = r_bf(i) -       thimp_bf     *dt*temp
          q_bf(i) = q_bf(i) + (1. - thimp_bf*bdf)*dt*temp
       end if


    else  ! on linear.eq.0
       ! Assumes no contribution from equilibrium f

        temp = tepsipsikappar(trial,lin,ps079,pp079,b2i79,kar79) &
            + tepsipsikappar(trial,ps079,lin,pp079,b2i79,kar79) &
             + tepsibkapparl  (trial,lin,bz079,pp079,b2i79,kar79)
        ssterm(i,psi_g) = ssterm(i,psi_g) -       thimpf     *dt*temp
        ddterm(i,psi_g) = ddterm(i,psi_g) + (1. - thimpf*bdf)*dt*temp

       if(numvar.ge.2) then
          temp = tepsibkapparl(trial,ps079,lin,pp079,b2i79,kar79) &
               + tebbkapparl  (trial,lin,bz079,pp079,b2i79,kar79) &
               + tebbkapparl  (trial,bz079,lin,pp079,b2i79,kar79)
          ssterm(i,bz_g) = ssterm(i,bz_g) -       thimpf     *dt*temp
          ddterm(i,bz_g) = ddterm(i,bz_g) + (1. - thimpf*bdf)*dt*temp
       end if

        temp = tepsipsikappar(trial,ps079,ps079,lin,b2i79,kar79) &
             + tepsibkapparl  (trial,ps079,bz079,lin,b2i79,kar79) &
             + tebbkapparl    (trial,bz079,bz079,lin,b2i79,kar79)
        ssterm(i,pp_g) = ssterm(i,pp_g) -       thimp     *dt*temp
        ddterm(i,pp_g) = ddterm(i,pp_g) + (1. - thimp*bdf)*dt*temp
 !

       if(i3d.eq.1 .and. numvar.ge.2) then
          temp = tepsifkappar(trial,ps079,lin,pp079,b2i79,kar79) &
               + tebfkapparl (trial,bz079,lin,pp079,b2i79,kar79)
          r_bf(i) = r_bf(i) -       thimp_bf     *dt*temp
          q_bf(i) = q_bf(i) + (1. - thimp_bf*bdf)*dt*temp
       endif

    endif  ! on linear.eq.0
  endif  ! on kappar.ne.0
 

  ! Cross-field Heat Flux
  ! ~~~~~~~~~~~~~~~~~~~~~
  if(kappax.ne.0.) then
     temp = p1kappax(trial,pp179,lin,ni79,kax79)
     ssterm(i,bz_g) = ssterm(i,bz_g) - thimp*dt*temp
     ddterm(i,bz_g) = ddterm(i,bz_g) - thimp*dt*temp*bdf

     temp = p1kappax(trial,lin,bzt79,ni79,kax79)
     ssterm(i,pp_g) = ssterm(i,pp_g) -     thimp     *dt*temp
     ddterm(i,pp_g) = ddterm(i,pp_g) + (1.-thimp*bdf)*dt*temp
  endif

  ! terms due to time-dependent density
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(idens.eq.1 .and. iadiabat.eq.1) then
     temp = t3tndenm(trial,lin,nt79,denm) + t3ts(trial,lin,sig79)
     ssterm(i,pp_g) = ssterm(i,pp_g) -     thimp     *dt*temp
     ddterm(i,pp_g) = ddterm(i,pp_g) + (1.-thimp*bdf)*dt*temp
  endif

  end do
end subroutine temperature_lin


subroutine pressure_nolin(trialx, r4term, total_pressure)

  use basic
  use m3dc1_nint
  use metricterms_new

  implicit none

  vectype, intent(in), dimension(MAX_PTS, OP_NUM, dofs_per_element) :: trialx
  vectype, intent(out), dimension(dofs_per_element) :: r4term

  vectype, dimension(MAX_PTS, OP_NUM)  :: trial
  integer :: i

  vectype, dimension(MAX_PTS, OP_NUM) :: pp079
  logical, intent(in) :: total_pressure

  if(itemp.eq.0) then
     if(total_pressure) then
        pp079 = p079
     else
        pp079 = pe079
     end if
  else
     if(total_pressure) then
        pp079 = ti079
     else
        pp079 = te079
     endif
  endif
  
  r4term = 0.

  do i=1, dofs_per_element
     trial = trialx(:,:,i)

  ! Contribution from external fields
  if(use_external_fields .and. (eqsubtract.eq.1 .or. icsubtract.eq.1)) then
     if(kappar.ne.0.) then
        r4term(i) = r4term(i) + dt* &
             (p1psipsipnkappar(trial,psx79,ps079,pp079,n079,ieq_bdotgradt) &
             +p1psipsipnkappar(trial,ps079,psx79,pp079,n079,1) &
             +p1psibpnkappar  (trial,psx79,bz079,pp079,n079,ieq_bdotgradt,1) &
             +p1psibpnkappar  (trial,ps079,psx79,pp079,n079,1,ieq_bdotgradt) &
             +p1bbpnkappar    (trial,bzx79,bz079,pp079,n079,ieq_bdotgradt) &
             +p1bbpnkappar    (trial,bz079,bzx79,pp079,n079,1))

        if(i3d.eq.1 .and. numvar.ge.2) then
           r4term(i) = r4term(i) + dt* &
                (p1psifpnkappar(trial,ps079,bfx79,pp079,n079,1,ieq_bdotgradt) &
                +p1bfpnkappar  (trial,bz079,bfx79,pp079,n079,1,ieq_bdotgradt))
        endif
     endif ! on kappar.ne.0
  endif  

  if(linear.eq.1) cycle

  ! source terms
  ! ~~~~~~~~~~~~
  if(gam.ne.1.) then
     ! heat source
     r4term(i) = r4term(i) + dt*(gam-1.)*(b3q(trial,q79) + b3q(trial,rad79))

     ! hyper-ohmic heating
     if(db.ne.0.) then 
        r4term(i) = r4term(i) + db*dt*(gam-1.)* &
             (qpsipsieta(trial) &
             +qbbeta    (trial))
     endif

     ! viscous heating
     if(total_pressure) then
        if(eqsubtract.eq.1) then
           r4term(i) = r4term(i) - dt*(gam-1.)* &
                (quumu    (trial,ph079,ph179,vis79      ) &
                +qvvmu    (trial,vz079,vz179,vis79      ) &
                +quchimu  (trial,ph079,ch179,vis79,vic79) &
                +qchichimu(trial,ch079,ch179,      vic79))
           r4term(i) = r4term(i) - dt*(gam-1.)* &
                (quumu    (trial,ph179,ph079,vis79      ) &
                +qvvmu    (trial,vz179,vz079,vis79      ) &
                +quchimu  (trial,ph179,ch079,vis79,vic79) &
                +qchichimu(trial,ch179,ch079,      vic79))
           r4term(i) = r4term(i) - dt*(gam-1.)* &
                (p1vip    (trial))
        else
           r4term(i) = r4term(i) - dt*(gam-1.)* &
                (quumu    (trial,pht79,pht79,vis79      ) &
                +qvvmu    (trial,vzt79,vzt79,vis79      ) &
                +quchimu  (trial,pht79,cht79,vis79,vic79) &
                +qchichimu(trial,cht79,cht79,      vic79) &
                +p1vip    (trial))
        end if
     endif
  end if

  ! density source terms
  ! ~~~~~~~~~~~~~~~~~~~~
  if(idens.eq.1 .and. eqsubtract.eq.1) then
     if(total_pressure) then
        r4term(i) = r4term(i) + dt* &
             (p1uus    (trial,ph079,ph079,sig79) &
             +p1vvs    (trial,vz079,vz079,sig79) &
             +p1chichis(trial,ch079,ch079,sig79) &
             +p1uchis  (trial,ph079,ch079,sig79))
     endif
     if(itemp.eq.1 .and. iadiabat.eq.1) then
        r4term(i) = r4term(i) + dt* &
             (t3tndenm(trial,pp079,n179,denm) &
             +t3ts(trial,pp079,sig79))
     endif
  endif

  end do
end subroutine pressure_nolin


!======================================================================
! bf_equation
!======================================================================
subroutine bf_equation_lin(trial, lin, ssterm, ddterm, r_bf, q_bf)
  
  use basic
  use arrays
  use m3dc1_nint
  use metricterms_new

  implicit none

  vectype, dimension(MAX_PTS, OP_NUM, dofs_per_element), intent(in) :: trial
  vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: lin 
  vectype, dimension(dofs_per_element,num_fields), intent(out) :: ssterm, ddterm
  vectype, dimension(dofs_per_element), intent(out) :: r_bf, q_bf

  ssterm = 0.
  ddterm = 0.
  r_bf = 0.
  q_bf = 0.

  r_bf = - intx3(trial(:,OP_1,:),r2_79,lin(:,OP_LP))
  if(ifbound.eq.2) then 
     r_bf = r_bf - regular*intx3(trial(:,OP_1,:),r2_79,lin(:,OP_1))
  end if
  ssterm(:,bz_g) = intx2(trial(:,OP_1,:),lin(:,OP_1 ))
end subroutine bf_equation_lin

!======================================================================
! bf_equation
!======================================================================
subroutine bf_equation_nolin(trial, r4term)
  
  use basic
  use arrays
  use m3dc1_nint
  use metricterms_new

  implicit none

  vectype, intent(in), dimension(MAX_PTS, OP_NUM, dofs_per_element)  :: trial
  vectype, intent(out), dimension(dofs_per_element) :: r4term

  if(eqsubtract.eq.1) then
     r4term = 0.
     return
  end if

  if(itor.eq.0) then
     temp79a = bzero
  else
     temp79a = bzero*rzero
  end if
  r4term = intx2(trial(:,OP_1,:),temp79a)
end subroutine bf_equation_nolin


!======================================================================
! j_equation
!======================================================================
subroutine j_equation_lin(trial, lin, ssterm, ddterm, r_bf, q_bf)
  
  use basic
  use arrays
  use m3dc1_nint
  use metricterms_new

  implicit none

  vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: trial, lin 
  vectype, dimension(num_fields), intent(out) :: ssterm, ddterm
  vectype, intent(out) :: r_bf, q_bf

  vectype :: temp

  ssterm = 0.
  ddterm = 0.
  r_bf = 0
  q_bf = 0

  ssterm(e_g) =     int2(trial(:,OP_1),lin(:,OP_1 ))
  select case(imp_hyper)
  case(1)  ! note e_g = - delstar(psi)

    ssterm(psi_g) =   int2(trial(:,OP_1),lin(:,OP_GS))
  
  case(2)  ! note e_g = J dot B / B**2

    temp = j1b2ipsib(trial,b2i79,lin,bzt79)
    ssterm(psi_g) = ssterm(psi_g) - temp

    temp = j1b2ibpsi(trial,b2i79,pst79,lin)
    ssterm(bz_g) = ssterm(bz_g) - temp

#if defined(USE3D) || defined(USECOMPLEX)

    temp = j1b2ipsif(trial,b2i79,pst79,lin)
    r_bf = r_bf - temp

    temp = j1b2ifb(trial,b2i79,lin,bzt79)
    r_bf = r_bf - temp

    temp = j1b2iff(trial,b2i79,bft79,lin)
    r_bf = r_bf - temp

    temp = j1b2ipsipsi(trial,b2i79,lin,pst79)
    ssterm(psi_g) = ssterm(psi_g) - temp

#endif

  end select
end subroutine j_equation_lin



!======================================================================
! ludefall
! --------
!
! Clears, populates, and finalizes all matrices for the implicit 
! time advance.  Does not insert boundary conditions, or finalize the
! s* matrices.
!
!======================================================================
subroutine ludefall(ivel_def, idens_def, ipres_def, ipressplit_def,  ifield_def)

  use mesh_mod
  use basic
  use arrays
  use sparse
  use m3dc1_nint
  use diagnostics
  use boundary_conditions
  use time_step
  use matrix_mod
  use transport_coefficients
  use gyroviscosity
  use runaway_mod

  implicit none

  integer, intent(in) :: ivel_def   ! populate velocity advance matrices
  integer, intent(in) :: idens_def  ! populate density advance matrices
  integer, intent(in) :: ipres_def  ! populate pressure advance matrices
  integer, intent(in) :: ipressplit_def  ! also populate pressure advance matrices
  integer, intent(in) :: ifield_def ! populate field advance matrices 

  integer :: itri, numelms, izone
  integer :: def_fields

  real :: tstart, tend, tfield, telm, tsizefield, tfinalize
  integer :: is_edge(3)  ! is inode on boundary
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

  call clear_matrices

  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, " populating matrices..."
    
  ! Specify which fields will be used in matrix population
  def_fields = FIELD_PSI + FIELD_I + FIELD_P &
       + FIELD_PHI + FIELD_V + FIELD_CHI &
       + FIELD_ETA + FIELD_MU &
       + FIELD_N + FIELD_NI

  if(numvar.ge.3 .or. ipres.eq.1) then
     def_fields = def_fields + FIELD_KAP
  end if

  if(itemp.gt.0) then
     def_fields = def_fields + FIELD_TE 
     if(ipres.gt.0) then
        def_fields = def_fields + FIELD_TI
     endif
  endif

  if(idens_def.eq.1) then
     if(density_source) def_fields = def_fields + FIELD_SIG
  endif
  if(momentum_source) def_fields = def_fields + FIELD_F
!...poloidal momentum source
  if(ipforce.gt.0) def_fields = def_fields + FIELD_PF
  if(heat_source) def_fields = def_fields + FIELD_Q
  if(icd_source.gt.0) def_fields = def_fields + FIELD_CD
  if(rad_source) def_fields = def_fields + FIELD_RAD

  if(gyro.eq.1 .or. amupar.ne.0 .or. kappar.ne.0 .or. kinetic.ne.0) then
     def_fields = def_fields + FIELD_B2I
  endif

  if(numvar.ge.3 .or. ipres.eq.1) then
     if(hyper.eq.0.) def_fields = def_fields + FIELD_J
     if(hyperc.ne.0.) def_fields = def_fields + FIELD_VOR + FIELD_COM
  end if

  if(kinetic.gt.0) then
     def_fields = def_fields + FIELD_KIN
  endif

  if(irunaway.gt.0) then
     def_fields = def_fields + FIELD_RE
  end if

  if(integrator.eq.1 .and. ntime.gt.1) then
     bdf = 2.
  else
     bdf = 1.
  endif

  ! Loop over elements
  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, " begin loop over elements"
  do itri=1,numelms

     call get_zone(itri, izone)

     ! calculate the field values and derivatives at the sampling points
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     call define_element_quadrature(itri, int_pts_main, int_pts_tor)
     call define_fields(itri, def_fields, 1, linear)
     if(gyro.eq.1) call gyro_common
     if(irunaway.gt.0) call eval_runaway(itri,izone)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        tfield = tfield + tend - tstart
     endif
     
     ! add element's contribution to matrices
     if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)
     if(ivel_def.eq.1) call ludefvel_n(itri)
     if(ifield_def.eq.1) call ludefphi_n(itri)
     if(idens_def.eq.1) call ludefden_n(itri)
     if(ipres_def.eq.1 .or. ipressplit_def.eq.1) call ludefpres_n(itri)
     if(myrank.eq.0 .and. itimer.eq.1) then
        call second(tend)
        telm = telm + tend - tstart
     endif

     if(isurface.eq.0) cycle
     if(.not.(nonrect.eq.1 .and. ivform.eq.1)) cycle

     ! add surface terms
     call boundary_edge(itri, is_edge, n, idim)
     
     do iedge=1,3
        if(is_edge(iedge).eq.0) cycle

        call define_boundary_quadrature(itri, iedge, 5, 5, n, idim)
        call define_fields(itri, def_fields, 1, linear)

!        write(*,'(A,8F10.4)') 'EDGE: ', x_79(1), z_79(1), x_79(5), z_79(5), &
!           norm79(3,1), norm79(3,2)

        if(gyro.eq.1) call gyro_common

        if(ivel_def.eq.1) call ludefvel_n(itri)
        if(ifield_def.eq.1) call ludefphi_n(itri)
        if(idens_def.eq.1) call ludefden_n(itri)
        if(ipres_def.eq.1 .or. ipressplit_def.eq.1) call ludefpres_n(itri)
     end do
  end do

  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, " finalizing matrices..."

  if(myrank.eq.0 .and. itimer.eq.1) call second(tstart)

  ! since a proc is contributing values to parts of the vector
  ! it does not own, we call sumsharedppplvecvals so that these values
  ! get summed up for all values shared by multiple procs
  ! and then update these values

  call finalize_matrices

  if(myrank.eq.0 .and. iprint.ge.1) &
       print *, " done finalizing matrices."

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
  use model
  use boundary_conditions

  implicit none

  integer, intent(in) :: itri

  integer :: i, j

  vectype, dimension(dofs_per_element,dofs_per_element,num_fields) :: ss, dd
  vectype, dimension(dofs_per_element,dofs_per_element) :: r_bf, q_bf
  vectype, dimension(dofs_per_element) :: r4

  type(matrix_type), pointer :: vv1, vb1, vn1, vf1, vp1
  type(matrix_type), pointer :: vv0, vb0, vn0, vf0, vp0
  type(vector_type), pointer :: vsource
  integer :: advfield
  integer :: pp_i
  integer, dimension(3) :: ieq
  integer :: k, izone
  integer, dimension(dofs_per_element) :: imask

  call get_zone(itri, izone)

  if(isplitstep.ge.1) then
     vv1 => s1_mat
     vv0 => d1_mat
     vb0 => q1_mat
     vn0 => r14_mat
     if(imp_bf.eq.1) then
        vf1 => q1_mat
        vf0 => q1_mat
     else
        vf0 => o1_mat
     end if
     vsource => r4_vec

     ! In splitstep case, pe_i spot in phivec is occupied by total pressure
     ! in velocity advance for numvar==3
     if((ipres.eq.1 .and. numvar.lt.3) .or. ipressplit.eq.1) then
        vp0 => p1_mat
        pp_i = p_i
     else
        vp0 => q1_mat
        pp_i = pe_i
     end if    
  else
     vv1 => s1_mat
     vv0 => d1_mat
     vb1 => s1_mat
     vb0 => d1_mat
     vp1 => s1_mat
     vp0 => d1_mat
     vn1 => s1_mat
     vn0 => d1_mat
     vsource => q4_vec
     pp_i = p_i
     if(imp_bf.eq.1) then
        vf1 => s1_mat
        vf0 => d1_mat
     else
        vf0 => o1_mat
     end if
  endif

  if(isplitstep.ge.1 .and. iestatic.eq.0) then
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
     r_bf = 0.
     q_bf = 0.
     r4 = 0.

     select case(k)
     case(1)
        do j=1,dofs_per_element
           call vorticity_lin(mu79,nu79(:,:,j), &
                ss(:,j,:),dd(:,j,:),r_bf(:,j),q_bf(:,j),advfield,izone)
        end do
     case(2)
        do j=1,dofs_per_element
           call axial_vel_lin(mu79,nu79(:,:,j), &
                ss(:,j,:),dd(:,j,:),r_bf(:,j),q_bf(:,j),advfield,izone)
        end do
     case(3)
        do j=1,dofs_per_element
           call compression_lin(mu79,nu79(:,:,j), &
                ss(:,j,:),dd(:,j,:),r_bf(:,j),q_bf(:,j),advfield,izone)
        end do
     end select

     if(izone.eq.1) then 
        select case(k)
        case(1)
           call vorticity_nolin(mu79,r4)
        case(2)
           call axial_vel_nolin(mu79,r4)
        case(3)
           call compression_nolin(mu79,r4)
        end select
     end if


     if(idifv .gt. 0) dd(:,:,  u_g) = dd(:,:,  u_g) - ss(:,:,  u_g)


     ! Zero-out rows that will be used for boundary conditions
     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     select case(k)
     case(1)
        call get_vor_mask(itri, imask)
     case(2)
        call get_vz_mask(itri, imask)
     case(3)
        call get_chi_mask(itri, imask)
     end select
     do i=1, num_fields
        call apply_boundary_mask(itri, 0, ss(:,:,i), imask)
        call apply_boundary_mask(itri, 0, dd(:,:,i), imask)
     end do
     call apply_boundary_mask(itri, 0, r_bf, imask)
     call apply_boundary_mask(itri, 0, q_bf, imask)
     call apply_boundary_mask_vec(itri, 0, r4, imask)


     ! Insert values into matrix
     ! ~~~~~~~~~~~~~~~~~~~~~~~~~
     call insert_block(vv1,itri,ieq(k),  u_i,ss(:,:,  u_g),MAT_ADD)
     call insert_block(vv0,itri,ieq(k),  u_i,dd(:,:,  u_g),MAT_ADD)
     call insert_block(vb0,itri,ieq(k),psi_i,dd(:,:,psi_g),MAT_ADD)
     if(numvar.ge.2) then
        if(idifv .gt. 0) dd(:,:,  vz_g) = dd(:,:,  vz_g) - ss(:,:,  vz_g)
        call insert_block(vv1,itri,ieq(k), vz_i,ss(:,:,vz_g),MAT_ADD)
        call insert_block(vv0,itri,ieq(k), vz_i,dd(:,:,vz_g),MAT_ADD)
        call insert_block(vb0,itri,ieq(k), bz_i,dd(:,:,bz_g),MAT_ADD)
     endif
     if(numvar.ge.3) then
        if(idifv .gt. 0) dd(:,:,  chi_g) = dd(:,:,  chi_g) - ss(:,:,  chi_g)
        call insert_block(vv1,itri,ieq(k),chi_i,ss(:,:,chi_g),MAT_ADD)
        call insert_block(vv0,itri,ieq(k),chi_i,dd(:,:,chi_g),MAT_ADD)
     endif
     if(numvar.ge.3 .or. ipres.eq.1) then
        call insert_block(vp0,itri,ieq(k), pp_i,dd(:,:,  p_g),MAT_ADD)
     end if
     if(idens.eq.1) then
        call insert_block(vn0,itri,ieq(k),den_i,dd(:,:,den_g),MAT_ADD)
     endif
     if(i3d.eq.1 .and. numvar.ge.2) then
        call insert_block(vf0,itri,ieq(k),bf_i,q_bf(:,:),MAT_ADD)
        if(imp_bf.eq.1) then
           if(isplitstep.ge.1) r_bf = -r_bf
           call insert_block(vf1,itri,ieq(k),bf_i,r_bf(:,:),MAT_ADD)
        end if
     endif
     if(isplitstep.eq.0) then
        call insert_block(vb1,itri,ieq(k),psi_i,ss(:,:,psi_g),MAT_ADD)
        if(numvar.ge.2) &
             call insert_block(vb1,itri,ieq(k), bz_i,ss(:,:, bz_g),MAT_ADD)
        if(numvar.ge.3 .or. ipres.eq.1) &
             call insert_block(vp1,itri,ieq(k), pp_i,ss(:,:,  p_g),MAT_ADD)
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
  use time_step
  use model
  use boundary_conditions

  implicit none

  integer, intent(in) :: itri

  integer :: i, j, k, izone
  
  vectype, dimension(dofs_per_element,dofs_per_element,num_fields) :: ss, dd
  vectype, dimension(dofs_per_element,dofs_per_element) :: r_bf, q_bf
  vectype, dimension(dofs_per_element) :: q4
  vectype, dimension(dofs_per_element,dofs_per_element,2) :: q_ni

  type(matrix_type), pointer :: bb1, bb0, bv1, bv0, bf0, bf1, bn1, bn0
  type(vector_type), pointer :: bsource
  integer :: ieq(5)
  integer :: maxk, pp_i, ppe_i
  integer :: imask(dofs_per_element)

  call get_zone(itri, izone)

  ! pp_i is the column multiplying p
  ! ppe_i is the column multiplying pe

  if(isplitstep.ge.1) then
     bb1 => s2_mat
     bb0 => d2_mat
     bv1 => r2_mat
     bv0 => q2_mat
     bn1 => r42_mat
     bn0 => q42_mat
     if(imp_bf.eq.1) then
        bf1 => s2_mat
        bf0 => d2_mat
     else
        bf0 => o2_mat
     end if
     bsource => q4_vec
     pp_i = pe_i
     ppe_i = pe_i
     maxk = vecsize_phi
  else
     bb1 => s1_mat
     bb0 => d1_mat
     bv1 => s1_mat
     bv0 => d1_mat
     bn1 => s1_mat
     bn0 => d1_mat
     if(imp_bf.eq.1) then
        bf1 => s1_mat
        bf0 => d1_mat
     else
        bf0 => o1_mat
     end if
     bsource => q4_vec
     pp_i = p_i
     if(ipres.eq.1 .and. numvar.ge.3) then
        ppe_i = pe_i
     else
        ! If electron pressure equation is not included,
        ! terms that should multiply pe will instead multiply p
        ! (we add a factor of pefac to these terms below)
        ppe_i = p_i
     end if
     maxk = numvar
     if(imp_bf.eq.1 .and. numvar.ge.2) maxk = maxk + 1
!!$     if(ipressplit.ne.0 .and. numvar.ge.3) maxk = maxk - 1   !  is this generally valid?
     if((jadv.eq.0).or.(jadv.eq.1 .and. imp_hyper.ge.1)) &
          maxk = maxk + 1
  endif

  ieq(1) = psi_i
  ieq(2) =  bz_i
  if(ipressplit.eq.0 .and. numvar.ge.3) then
     ieq(3) = ppe_i
  end if

  ! add bf equation and e equations:
  ! NOTE:  e=electrostatic potential for jadv=0, del_star_psi for jadv=1
  if(imp_bf.eq.1) then
     if((jadv.eq.0).or.(jadv.eq.1 .and. imp_hyper.ge.1)) then
        ieq(maxk-1) = bf_i
        ieq(maxk)   = e_i
     else
        ieq(maxk) = bf_i
     endif
  else
     if((jadv.eq.0).or.(jadv.eq.1 .and. imp_hyper.ge.1)) then
        ieq(maxk) = e_i
     endif
  endif

  do k=1,maxk
     ss = 0.
     dd = 0.
     q_ni = 0.
     r_bf = 0.
     q_bf = 0.
     q4 = 0.
    
     do j=1,dofs_per_element
        if     (ieq(k).eq.psi_i) then
           if(.not.surface_int) then
              call flux_lin(mu79,nu79(:,:,j), &
                   ss(:,j,:),dd(:,j,:),q_ni(:,j,:),r_bf(:,j),q_bf(:,j),izone)
           end if
        else if(ieq(k).eq.bz_i .and. numvar.ge.2) then
           call axial_field_lin(mu79,nu79(:,:,j), &
                ss(:,j,:),dd(:,j,:),q_ni(:,j,:),r_bf(:,j),q_bf(:,j), &
                izone)
        else if(ieq(k).eq.ppe_i .and. ipressplit.eq.0 .and. numvar.ge.3) then
              ! if ipres==0, this is the total pressure equation
              ! if ipres==1, this is the electron pressure equation
           call pressure_lin(mu79,nu79(:,:,j), &
                   ss(:,j,:),dd(:,j,:),q_ni(:,j,:),r_bf(:,j),q_bf(:,j), &
                   ipres.eq.0, thimp, izone)
        else if(ieq(k).eq.bf_i .and. imp_bf.eq.1) then
           call bf_equation_lin(mu79,nu79(:,:,j), &
                ss(:,j,:),dd(:,j,:),r_bf(:,j),q_bf(:,j))
        else if(ieq(k).eq.e_i) then
           if(jadv.eq.0) then
              do i=1,dofs_per_element
                 call potential_lin(mu79(:,:,i),nu79(:,:,j), &
                      ss(i,j,:),dd(i,j,:),q_ni(i,j,1),r_bf(i,j),q_bf(i,j))
              end do
           else   !jadv.eq.1
              do i=1,dofs_per_element
                 call j_equation_lin(mu79(:,:,i),nu79(:,:,j), &
                      ss(i,j,:),dd(i,j,:),r_bf(i,j),q_bf(i,j))   
              end do
           endif
           
        else
           print *, 'error!'
           call safestop(32)
        end if
     end do

     if(izone.eq.1) then
        if(ieq(k).eq.psi_i) then
           call flux_nolin(mu79,q4)
        else if(ieq(k).eq.bz_i .and. numvar.ge.2) then
           call axial_field_nolin(mu79,q4)
        else if(ieq(k).eq.ppe_i .and. ipressplit.eq.0 .and. numvar.ge.3) then
           call pressure_nolin(mu79,q4,ipres.eq.0)
        else if(ieq(k).eq.bf_i .and. imp_bf.eq.1) then
           call bf_equation_nolin(mu79,q4)
        end if
     end if
    
     ! If electron pressure equation is not included,
     ! terms that should multiply pe will instead multiply p
     ! so we multiply these terms my pefac
     if(.not.(ipres.eq.1 .and. numvar.ge.3) .and. ipressplit.eq.0) then
        ss(:,:,pe_g) = ss(:,:,pe_g)*pefac
        dd(:,:,pe_g) = dd(:,:,pe_g)*pefac
     end if


     ! Zero-out rows that will be used for boundary conditions
     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     if     (ieq(k).eq.psi_i) then
        call get_flux_mask(itri, imask)
     else if(ieq(k).eq.bz_i .and. numvar.ge.2) then
        call get_bz_mask(itri, imask)
     else if(ieq(k).eq.ppe_i .and. ipressplit.eq.0 .and. numvar.ge.3) then
        call get_pres_mask(itri, imask)
     else if(ieq(k).eq.bf_i .and. imp_bf.eq.1) then
        call get_bf_mask(itri, imask)
     else if(ieq(k).eq.e_i) then
        if(jadv.eq.0) then
           call get_esp_mask(itri, imask)
        else
           call get_j_mask(itri, imask)
        endif
     else
        if(myrank.eq.0) print *, 'error in ludefphi_n', k, ieq(k), maxk
        call safestop(31)
     end if
     do i=1, num_fields
        call apply_boundary_mask(itri, 0, ss(:,:,i), imask)
        call apply_boundary_mask(itri, 0, dd(:,:,i), imask)
     end do
     do i=1, 2
        call apply_boundary_mask(itri, 0, q_ni(:,:,i), imask)
     end do
     call apply_boundary_mask(itri, 0, r_bf, imask)
     call apply_boundary_mask(itri, 0, q_bf, imask)
     call apply_boundary_mask_vec(itri, 0, q4, imask)


     ! Insert values into matrix
     ! ~~~~~~~~~~~~~~~~~~~~~~~~~
     if(idiff .gt. 0) dd(:,:,psi_g) = dd(:,:,psi_g) - ss(:,:,psi_g)
     call insert_block(bb1,itri,ieq(k),psi_i,ss(:,:,psi_g),MAT_ADD)
     call insert_block(bb0,itri,ieq(k),psi_i,dd(:,:,psi_g),MAT_ADD)
     call insert_block(bv1,itri,ieq(k),  u_i,ss(:,:,  u_g),MAT_ADD)
     call insert_block(bv0,itri,ieq(k),  u_i,dd(:,:,  u_g),MAT_ADD)
     if((jadv.eq.0) .or. (jadv.eq.1 .and. imp_hyper.ge.1)) then
        if(idiff .gt. 0) dd(:,:,e_g)  = dd(:,:,e_g)  - ss(:,:,e_g)
        call insert_block(bb1,itri,ieq(k),e_i,ss(:,:,e_g),MAT_ADD)
        call insert_block(bb0,itri,ieq(k),e_i,dd(:,:,e_g),MAT_ADD)     
     endif
     if(numvar.ge.2) then
        if(idiff .gt. 0) dd(:,:,bz_g) = dd(:,:,bz_g) - ss(:,:,bz_g)
        call insert_block(bb1,itri,ieq(k), bz_i,ss(:,:,bz_g),MAT_ADD)
        call insert_block(bb0,itri,ieq(k), bz_i,dd(:,:,bz_g),MAT_ADD)
        call insert_block(bv1,itri,ieq(k), vz_i,ss(:,:,vz_g),MAT_ADD)
        call insert_block(bv0,itri,ieq(k), vz_i,dd(:,:,vz_g),MAT_ADD)
     endif
     if(numvar.ge.3) then
        ! if ipres==0, total pressure equation is in pe_i slot
        if(ipressplit.eq.0) then
!           if(ipres.eq.0) then
              if(idiff .gt. 0) dd(:,:,p_g) = dd(:,:,p_g) - ss(:,:,p_g)
              call insert_block(bb1,itri,ieq(k), pp_i,ss(:,:,  p_g),MAT_ADD)
              call insert_block(bb0,itri,ieq(k), pp_i,dd(:,:,  p_g),MAT_ADD)
!           end if
           if(idiff .gt. 0) dd(:,:,pe_g) = dd(:,:,pe_g) - ss(:,:,pe_g)
           call insert_block(bb1,itri,ieq(k), ppe_i,ss(:,:, pe_g),MAT_ADD)
           call insert_block(bb0,itri,ieq(k), ppe_i,dd(:,:, pe_g),MAT_ADD)
        end if
        call insert_block(bv1,itri,ieq(k),chi_i,ss(:,:,chi_g),MAT_ADD)
        call insert_block(bv0,itri,ieq(k),chi_i,dd(:,:,chi_g),MAT_ADD)
     endif
     if(idens.eq.1) then
        call insert_block(bn1,itri,ieq(k),den_i,ss(:,:,den_g),MAT_ADD)
        call insert_block(bn0,itri,ieq(k),den_i,dd(:,:,den_g),MAT_ADD)
     endif
     if(i3d.eq.1 .and. numvar.ge.2) then
        if(imp_bf.eq.1 .and. idiff .gt. 0) q_bf(:,:) = q_bf(:,:) - r_bf(:,:)
        call insert_block(bf0,itri,ieq(k),bf_i,q_bf(:,:),MAT_ADD)
        if(imp_bf.eq.1) then  
           call insert_block(bf1,itri,ieq(k),bf_i,r_bf(:,:),MAT_ADD)
        end if
     endif

     call vector_insert_block(bsource,itri,ieq(k),q4,VEC_ADD)
  end do
end subroutine ludefphi_n

!======================================================================
! ludefpres_n
! -----------
!
! populates the matrices for pressure advance
!
! itri: index of finite element
!======================================================================
subroutine ludefpres_n(itri)
  use basic
  use mesh_mod
  use m3dc1_nint
  use arrays
  use sparse
  use electrostatic_potential
  use time_step
  use model
  use boundary_conditions

  implicit none

  integer, intent(in) :: itri

  integer :: i, j, k, izone
  
  vectype, dimension(dofs_per_element,dofs_per_element,num_fields) :: ss, dd
  vectype, dimension(dofs_per_element,dofs_per_element) :: r_bf, q_bf
  vectype, dimension(dofs_per_element) :: q4
  vectype, dimension(dofs_per_element,dofs_per_element,2) :: q_ni

  type(matrix_type), pointer :: pp1, pp0, pv1, pv0, pb0, pb1, pf0, pn1, pn0
  type(vector_type), pointer :: psource
  integer :: imask(dofs_per_element)
  integer :: ieq(2)
  integer :: maxk
  real :: thimpf

  call get_zone(itri, izone)

  if(isplitstep.ge.1) then
     pp1 => s9_mat
     pp0 => d9_mat
     pv1 => r9_mat
     pv0 => q9_mat
     pb0 => o9_mat
     pn1 => rp42_mat
     pn0 => qp42_mat
     if(imp_bf.eq.1) then
        pf0 => o9_mat
     else
        pf0 => o3_mat
     endif
     psource => qp4_vec
     thimpf = 0.
  else
     pp1 => s1_mat
     pp0 => d1_mat
     pv1 => s1_mat
     pv0 => d1_mat
     pb1 => s1_mat
     pb0 => d1_mat
     pn1 => s1_mat
     pn0 => d1_mat
     if(imp_bf.eq.1) then
        pf0 => d1_mat
     else
        pf0 => o1_mat
     endif
     psource => q4_vec
     thimpf = thimp
  endif
  select case(ipressplit)
  case(0)
     ieq(1) = p_i
     maxk = 1
  case(1)
     select case(imode)
     case(1)
        maxk = 1
        ieq(1) = p_i
     case(2)
        maxk = 1
        ieq(1) = te_i
     case(3)
        maxk = 2
        ieq(1) = p_i
        ieq(2) = pe_i
     case(4)
        maxk = 2
        ieq(1) = te_i
        ieq(2) = ti_i
     end select
  end select

  do k=1,maxk
     ss = 0.
     dd = 0.
     q_ni = 0.
     r_bf = 0.
     q_bf = 0.
     q4 = 0.

     do j=1,dofs_per_element
        if(ipressplit.eq.0) then
           call pressure_lin(mu79,nu79(:,:,j), &
                ss(:,j,:),dd(:,j,:),q_ni(:,j,:),r_bf(:,j),q_bf(:,j), &
                .true., thimpf, izone)
        else ! ipressplit=1
           select case(imode)
           case(1)
              call pressure_lin(mu79,nu79(:,:,j), &
                   ss(:,j,:),dd(:,j,:),q_ni(:,j,:),r_bf(:,j),q_bf(:,j), &
                   .true., thimpf, izone)
           case(2)
              call temperature_lin(mu79,nu79(:,:,j), &
                   ss(:,j,:),dd(:,j,:),q_ni(:,j,:),r_bf(:,j),q_bf(:,j), &
                   .true., thimpf, izone)
           case(3)
              if(k.eq.1) call pressure_lin(mu79,nu79(:,:,j), &
                   ss(:,j,:),dd(:,j,:),q_ni(:,j,:),r_bf(:,j),q_bf(:,j), &
                   .true., thimpf, izone)
              if(k.eq.2) call pressure_lin(mu79,nu79(:,:,j), &
                   ss(:,j,:),dd(:,j,:),q_ni(:,j,:),r_bf(:,j),q_bf(:,j), &
                   .false., thimpf, izone)
           case(4)
              if(k.eq.1) call temperature_lin(mu79,nu79(:,:,j), &
                   ss(:,j,:),dd(:,j,:),q_ni(:,j,:),r_bf(:,j),q_bf(:,j), &
                   .true., thimpf, izone)
              if(k.eq.2) call temperature_lin(mu79,nu79(:,:,j), &
                   ss(:,j,:),dd(:,j,:),q_ni(:,j,:),r_bf(:,j),q_bf(:,j), &
                   .false., thimpf, izone)
           end select
        endif ! on ipressplit
     end do  ! on j

     if(izone.eq.1) then
        if(ipressplit.eq.0) then
           call pressure_nolin(mu79,q4,.true.)
        endif
        if(ipressplit.eq.1) then
           if(imode.eq.1 .or. (imode.eq.3 .and. k.eq.1) &
                .or. (imode.eq.4 .and. k.eq.2))  then
              call pressure_nolin(mu79,q4,.true.)
           else
              call pressure_nolin(mu79,q4,.false.)
           endif
        end if
     endif  ! ipressplit

     ! Zero-out rows that will be used for boundary conditions
     call get_pres_mask(itri, imask)
     do i=1, num_fields
        call apply_boundary_mask(itri, 0, ss(:,:,i), imask)
        call apply_boundary_mask(itri, 0, dd(:,:,i), imask)
     end do
     do i=1, 2
        call apply_boundary_mask(itri, 0, q_ni(:,:,i), imask)
     end do
     call apply_boundary_mask(itri, 0, r_bf, imask)
     call apply_boundary_mask(itri, 0, q_bf, imask)
     call apply_boundary_mask_vec(itri, 0, q4, imask)


     if(ipressplit.eq.0) then
        if(idiff .gt. 0) dd(:,:, p_g) = dd(:,:, p_g) - ss(:,:, p_g)
        call insert_block(pp1,itri,ieq(k),  p_i,ss(:,:,  p_g),MAT_ADD)
        call insert_block(pp0,itri,ieq(k),  p_i,dd(:,:,  p_g),MAT_ADD)
        call insert_block(pv1,itri,ieq(k),  u_i,ss(:,:,  u_g),MAT_ADD)
        call insert_block(pv0,itri,ieq(k),  u_i,dd(:,:,  u_g),MAT_ADD)
        call insert_block(pb0,itri,ieq(k),psi_i,dd(:,:,psi_g),MAT_ADD)
        if(numvar.ge.2) then
           call insert_block(pv1,itri,ieq(k),vz_i,ss(:,:,vz_g),MAT_ADD)
           call insert_block(pv0,itri,ieq(k),vz_i,dd(:,:,vz_g),MAT_ADD)
           call insert_block(pb0,itri,ieq(k),bz_i,dd(:,:,bz_g),MAT_ADD)
        endif
        if(numvar.ge.3) then
           call insert_block(pv1,itri,ieq(k),chi_i,ss(:,:,chi_g),MAT_ADD)
           call insert_block(pv0,itri,ieq(k),chi_i,dd(:,:,chi_g),MAT_ADD)
           call insert_block(pb0,itri,ieq(k), pe_i,dd(:,:, pe_g),MAT_ADD)
        endif

        if(isplitstep.eq.0) then
           call insert_block(pb1,itri,ieq(k),psi_i,ss(:,:,psi_g),MAT_ADD)
           if(numvar.ge.2) call insert_block(pb1,itri,ieq(k),bz_i,ss(:,:,bz_g),MAT_ADD)
           if(numvar.ge.3) call insert_block(pb1,itri,ieq(k),pe_i,ss(:,:,pe_g),MAT_ADD)
        end if
     else ! ipressplit=1

        select case(imode)
        case(1)
           if(idiff .gt. 0) dd(:,:, p_g) = dd(:,:, p_g) - ss(:,:, p_g)
           call insert_block(pp1,itri,ieq(k),  p_i,ss(:,:,  p_g),MAT_ADD)
           call insert_block(pp0,itri,ieq(k),  p_i,dd(:,:,  p_g),MAT_ADD)
        case(2)
           if(idiff .gt. 0) dd(:,:, te_g) = dd(:,:, te_g) - ss(:,:, te_g)
           call insert_block(pp1,itri,ieq(k),  te_i,ss(:,:,  te_g),MAT_ADD)
           call insert_block(pp0,itri,ieq(k),  te_i,dd(:,:,  te_g),MAT_ADD)
        case(3)
           if(idiff .gt. 0) then
              dd(:,:, p_g) = dd(:,:, p_g) - ss(:,:, p_g)
              dd(:,:, pe_g) = dd(:,:, pe_g) - ss(:,:, pe_g)
           endif
           call insert_block(pp1,itri,ieq(k),  p_i,ss(:,:,  p_g),MAT_ADD)
           call insert_block(pp0,itri,ieq(k),  p_i,dd(:,:,  p_g),MAT_ADD)
           call insert_block(pp1,itri,ieq(k), pe_i,ss(:,:, pe_g),MAT_ADD)
           call insert_block(pp0,itri,ieq(k), pe_i,dd(:,:, pe_g),MAT_ADD)
        case(4)
           if(idiff .gt. 0) then
              dd(:,:, te_g) = dd(:,:, te_g) - ss(:,:, te_g)
              dd(:,:, ti_g) = dd(:,:, ti_g) - ss(:,:, ti_g)
           endif
           call insert_block(pp1,itri,ieq(k), te_i,ss(:,:, te_g),MAT_ADD)
           call insert_block(pp0,itri,ieq(k), te_i,dd(:,:, te_g),MAT_ADD)
           call insert_block(pp1,itri,ieq(k), ti_i,ss(:,:, ti_g),MAT_ADD)
           call insert_block(pp0,itri,ieq(k), ti_i,dd(:,:, ti_g),MAT_ADD)
        end select
        call insert_block(pv1,itri,ieq(k),  u_i,ss(:,:,  u_g),MAT_ADD)
        call insert_block(pv0,itri,ieq(k),  u_i,dd(:,:,  u_g),MAT_ADD)
        call insert_block(pb0,itri,ieq(k),psi_i, dd(:,:,psi_g),MAT_ADD)
!...the following line should not be needed as ss(:,:, psi_g) should be zero
        call insert_block(pb0,itri,ieq(k),psi_i,-ss(:,:,psi_g),MAT_ADD)
        call insert_block(pv1,itri,ieq(k), vz_i,ss(:,:, vz_g),MAT_ADD)
        call insert_block(pv0,itri,ieq(k), vz_i,dd(:,:, vz_g),MAT_ADD)
        call insert_block(pb0,itri,ieq(k), bz_i, dd(:,:, bz_g),MAT_ADD)
!...the following line should not be needed as ss(:,:, bz_g) should be zero
        call insert_block(pb0,itri,ieq(k), bz_i,-ss(:,:, bz_g),MAT_ADD)
        call insert_block(pv1,itri,ieq(k),chi_i,ss(:,:,chi_g),MAT_ADD)
        call insert_block(pv0,itri,ieq(k),chi_i,dd(:,:,chi_g),MAT_ADD)
     endif  ! on ipressplit

     call vector_insert_block(psource,itri,ieq(k),q4,VEC_ADD)

     if(idens.eq.1) then
        call insert_block(pn1,itri,ieq(k),den_i,ss(:,:,den_g),MAT_ADD)
        call insert_block(pn0,itri,ieq(k),den_i,dd(:,:,den_g),MAT_ADD)
     endif
     if(i3d.eq.1 .and. numvar.ge.2) then
        call insert_block(pf0,itri,ieq(k),bf_i, q_bf(:,:),MAT_ADD)
        call insert_block(pf0,itri,ieq(k),bf_i,-r_bf(:,:),MAT_ADD)
     endif

  enddo ! on k
end subroutine ludefpres_n



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
  use model
  use boundary_conditions

  implicit none

  integer, intent(in) :: itri

  integer :: i, j, izone
  vectype, dimension(dofs_per_element, dofs_per_element) :: ssterm, ddterm
  vectype, dimension(dofs_per_element, dofs_per_element, 3) :: rrterm, qqterm
  vectype, dimension(dofs_per_element) :: oterm

  vectype, dimension(dofs_per_element) :: tempx
  vectype, dimension(dofs_per_element,dofs_per_element) :: tempxx

  vectype :: freq_fac

  type(matrix_type), pointer :: nn1, nn0, nv1, nv0
  type(vector_type), pointer :: nsource
  real :: thimpb
  integer :: imask(dofs_per_element)

  call get_zone(itri, izone)

  if(imp_mod.eq.0) then
     thimpb = thimp
  else
     thimpb = 1.
  endif

  if(itime_independent.eq.1) then
#ifdef USECOMPLEX
     freq_fac = (0,1)*frequency*dt
#else
     freq_fac = 0.
#endif
  else
     freq_fac = 1.
  end if

  if(isplitstep.ge.1) then
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

  if(izone.ne.1) then
     tempxx = n1n(mu79,nu79)
     ssterm = ssterm + tempxx
     ddterm = ddterm + tempxx*bdf
     goto 400
  end if

  ! NUMVAR = 1
  ! ~~~~~~~~~~
  tempxx = n1n(mu79,nu79)*freq_fac
  ssterm = ssterm + tempxx
  if(itime_independent.eq.0) ddterm = ddterm + tempxx*bdf

  do j=1,dofs_per_element
     
     tempx = n1ndenm(mu79,nu79(:,:,j),denm,vz079) &
          +  n1nu   (mu79,nu79(:,:,j),pht79)
     ssterm(:,j) = ssterm(:,j) -     thimp     *dt*tempx
     ddterm(:,j) = ddterm(:,j) + (1.-thimp*bdf)*dt*tempx

     if(linear.eq.0) then
        tempx = n1nu(mu79,n179,nu79(:,:,j))
        rrterm(:,j,1) = rrterm(:,j,1) + thimpb*dt*tempx
        qqterm(:,j,1) = qqterm(:,j,1) - thimpb*dt*tempx*bdf
     endif

     if(eqsubtract.eq.1) then
        tempx = n1nu(mu79,n079,nu79(:,:,j))
        rrterm(:,j,1) = rrterm(:,j,1) +     thimpb     *dt*tempx
        qqterm(:,j,1) = qqterm(:,j,1) + (1.-thimpb*bdf)*dt*tempx
     endif

#if defined(USECOMPLEX) || defined(USE3D)
     ! NUMVAR = 2
     ! ~~~~~~~~~~
     if(numvar.ge.2) then
        tempx = n1nv(mu79,nu79(:,:,j),vzt79)
        ssterm(:,j) = ssterm(:,j) -     thimp     *dt*tempx
        ddterm(:,j) = ddterm(:,j) + (1.-thimp*bdf)*dt*tempx
        
        if(linear.eq.0) then 
           tempx = n1nv(mu79,n179,nu79(:,:,j))
           rrterm(:,j,2) = rrterm(:,j,2) + thimpb*dt*tempx
           qqterm(:,j,2) = qqterm(:,j,2) - thimpb*dt*tempx*bdf
        endif
        
        if(eqsubtract.eq.1) then
           tempx = n1nv(mu79,n079,nu79(:,:,j))
           rrterm(:,j,2) = rrterm(:,j,2) +     thimpb     *dt*tempx
           qqterm(:,j,2) = qqterm(:,j,2) + (1.-thimpb*bdf)*dt*tempx
        endif
     endif
#endif

     ! NUMVAR = 3
     ! ~~~~~~~~~~
     if(numvar.ge.3) then
        tempx = n1nchi(mu79,nu79(:,:,j),cht79)
        ssterm(:,j) = ssterm(:,j) -     thimp     *dt*tempx
        ddterm(:,j) = ddterm(:,j) + (1.-thimp*bdf)*dt*tempx
        
        if(linear.eq.0) then
           tempx = n1nchi(mu79,n179,nu79(:,:,j))
           rrterm(:,j,3) = rrterm(:,j,3) + thimpb*dt*tempx
           qqterm(:,j,3) = qqterm(:,j,3) - thimpb*dt*tempx*bdf
        endif
        
        if(eqsubtract.eq.1) then
           tempx = n1nchi(mu79,n079,nu79(:,:,j))
           rrterm(:,j,3) = rrterm(:,j,3) +     thimpb     *dt*tempx
           qqterm(:,j,3) = qqterm(:,j,3) + (1.-thimpb*bdf)*dt*tempx
        endif
     endif
  end do

  ! Source term
  ! ~~~~~~~~~~~~ 
  oterm = dt*n1s(mu79,sig79)



400 continue
  
  if(isplitstep.eq.0) rrterm = -rrterm
  if(idiff .gt. 0) ddterm = ddterm - ssterm


  ! Zero-out rows that will be used for boundary conditions
  call get_den_mask(itri, imask)
  call apply_boundary_mask(itri, 0, ssterm, imask)
  call apply_boundary_mask(itri, 0, ddterm, imask)
  do i=1, 3
     call apply_boundary_mask(itri, 0, rrterm(:,:,i), imask)
     call apply_boundary_mask(itri, 0, qqterm(:,:,i), imask)
  end do
  call apply_boundary_mask_vec(itri, 0, oterm, imask)


  ! Insert data into matrices
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


