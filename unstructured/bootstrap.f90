! Module containing terms 
!
module bootstrap
  
  implicit none

  integer :: ibootstrap_model
  ! 1 : add -eta*J_BS term to Ohm's law
  !     where J_BS = alpha I <p, psi> B
  !     In Sauter model, alpha = L31 / (p B^2 |Grad(psi)|^2)

  real :: bootstrap_alpha

contains

  subroutine bootstrap_flux(trial, lin, ssterm, ddterm, r_bf, q_bf, thimpf, thimp_bf)
    use basic
    use arrays
    use m3dc1_nint

    implicit none

    vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: trial, lin 
    vectype, dimension(num_fields), intent(inout) :: ssterm, ddterm
    vectype, intent(out) :: r_bf, q_bf
    real, intent(in) :: thimpf, thimp_bf

    vectype :: temp
    
    if(numvar.eq.1) then
       temp = bs_b1psipsib(trial,lin,pst79,bzt79) &
            + bs_b1psipsib(trial,pst79,lin,bzt79)
       ssterm(psi_g) = ssterm(psi_g) -          thimpf     *dt*temp
       ddterm(psi_g) = ddterm(psi_g) + (1./2. - thimpf*bdf)*dt*temp
       
       temp = bs_b1psibb  (trial,lin,bzt79,bzt79) &
            + bs_b1psibf  (trial,lin,bzt79,bft79)
       ssterm(psi_g) = ssterm(psi_g) -          thimpf     *dt*temp
       ddterm(psi_g) = ddterm(psi_g) + (1. - thimpf*bdf)*dt*temp
    else
       temp = bs_b1psipsib(trial,lin,pst79,bzt79) &
            + bs_b1psipsib(trial,pst79,lin,bzt79) &
            + bs_b1psibb  (trial,lin,bzt79,bzt79) &
            + bs_b1psibf  (trial,lin,bzt79,bft79)
       ssterm(psi_g) = ssterm(psi_g) -          thimpf     *dt*temp
       ddterm(psi_g) = ddterm(psi_g) + (1./3. - thimpf*bdf)*dt*temp
       
       temp = bs_b1psipsib(trial,pst79,pst79,lin) &
            + bs_b1psibb  (trial,pst79,lin,bzt79) &
            + bs_b1psibb  (trial,pst79,bzt79,lin) &
            + bs_b1psibf  (trial,pst79,lin,bft79)
       ssterm(bz_g) = ssterm(bz_g) -          thimpf     *dt*temp
       ddterm(bz_g) = ddterm(bz_g) + (1./3. - thimpf*bdf)*dt*temp
       
       temp = bs_b1psibf  (trial,pst79,bzt79,lin)
       r_bf = r_bf -          thimp_bf     *dt*temp
       q_bf = q_bf + (1./3. - thimp_bf*bdf)*dt*temp
       
       if(eqsubtract.eq.1) then
          temp = bs_b1psipsib(trial,lin,pst79,bz079) &
               + bs_b1psipsib(trial,lin,ps079,bzt79) &
               + bs_b1psipsib(trial,lin,ps079,bz079)*2. &
               + bs_b1psipsib(trial,pst79,lin,bz079) &
               + bs_b1psipsib(trial,ps079,lin,bzt79) &
               + bs_b1psipsib(trial,ps079,lin,bz079)*2. &
               + bs_b1psibb  (trial,lin,bzt79,bz079) &
               + bs_b1psibb  (trial,lin,bz079,bzt79) &
               + bs_b1psibb  (trial,lin,bz079,bz079)*2. &
               + bs_b1psibf  (trial,lin,bzt79,bf079) &
               + bs_b1psibf  (trial,lin,bz079,bft79) &
               + bs_b1psibf  (trial,lin,bz079,bf079)*2.
          ddterm(psi_g) = ddterm(psi_g) + (1./6.)*dt*temp
          
          temp = bs_b1psipsib(trial,pst79,ps079,lin) &
               + bs_b1psipsib(trial,ps079,pst79,lin) &
               + bs_b1psipsib(trial,ps079,ps079,lin)*2. &
               + bs_b1psibb  (trial,pst79,lin,bz079) &
               + bs_b1psibb  (trial,ps079,lin,bzt79) &
               + bs_b1psibb  (trial,ps079,lin,bz079)*2. &
               + bs_b1psibb  (trial,pst79,bz079,lin) &
               + bs_b1psibb  (trial,ps079,bzt79,lin) &
               + bs_b1psibb  (trial,ps079,bz079,lin)*2. &
               + bs_b1psibf  (trial,pst79,lin,bf079) &
               + bs_b1psibf  (trial,ps079,lin,bft79) &
               + bs_b1psibf  (trial,ps079,lin,bf079)*2.
          ddterm(bz_g) = ddterm(bz_g) + (1./6.)*dt*temp
          
          temp = bs_b1psibf  (trial,pst79,bz079,lin) &
               + bs_b1psibf  (trial,ps079,bzt79,lin) &
               + bs_b1psibf  (trial,ps079,bz079,lin)*2.
          q_bf = q_bf + (1./6.)*dt*temp
       end if
    end if
  end subroutine bootstrap_flux

  subroutine bootstrap_axial_field(trial, lin, ssterm, ddterm, r_bf, q_bf, thimpf, thimp_bf)
    use basic
    use arrays
    use m3dc1_nint

    implicit none

    vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: trial, lin 
    vectype, dimension(num_fields), intent(inout) :: ssterm, ddterm
    vectype, intent(out) :: r_bf, q_bf
    real, intent(in) :: thimpf, thimp_bf

    vectype :: temp

    temp = bs_b2psipsib(trial,lin,pst79,bzt79) &
         + bs_b2psipsib(trial,pst79,lin,bzt79) &
         + bs_b2psibf  (trial,lin,bzt79,bft79)
    ssterm(psi_g) = ssterm(psi_g) -          thimpf     *dt*temp
    ddterm(psi_g) = ddterm(psi_g) + (1./3. - thimpf*bdf)*dt*temp
    
    temp = bs_b2psipsib(trial,pst79,pst79,lin) &
         + bs_b2psibf  (trial,pst79,lin,bft79)
    ssterm(bz_g) = ssterm(bz_g) -          thimpf     *dt*temp
    ddterm(bz_g) = ddterm(bz_g) + (1./3. - thimpf*bdf)*dt*temp
    
    temp = bs_b2psibf  (trial,pst79,bzt79,lin)
    r_bf = r_bf -          thimp_bf     *dt*temp
    q_bf = q_bf + (1./3. - thimp_bf*bdf)*dt*temp
    
    if(eqsubtract.eq.1) then
       temp = bs_b2psipsib(trial,lin,pst79,bz079) &
            + bs_b2psipsib(trial,lin,ps079,bzt79) &
            + bs_b2psipsib(trial,lin,ps079,bz079)*2. &
            + bs_b2psipsib(trial,pst79,lin,bz079) &
            + bs_b2psipsib(trial,ps079,lin,bzt79) &
            + bs_b2psipsib(trial,ps079,lin,bz079)*2. &
            + bs_b2psibf  (trial,lin,bzt79,bf079) &
            + bs_b2psibf  (trial,lin,bz079,bft79) &
            + bs_b2psibf  (trial,lin,bz079,bf079)*2.
       ddterm(psi_g) = ddterm(psi_g) + (1./6.)*dt*temp
       
       temp = bs_b2psipsib(trial,pst79,ps079,lin) &
            + bs_b2psipsib(trial,ps079,pst79,lin) &
            + bs_b2psipsib(trial,ps079,ps079,lin)*2. &
            + bs_b2psibf  (trial,pst79,lin,bf079) &
            + bs_b2psibf  (trial,ps079,lin,bft79) &
            + bs_b2psibf  (trial,ps079,lin,bf079)*2.
       ddterm(bz_g) = ddterm(bz_g) + (1./6.)*dt*temp
       
       temp = bs_b2psibf  (trial,pst79,bz079,lin) &
            + bs_b2psibf  (trial,ps079,bzt79,lin) &
            + bs_b2psibf  (trial,ps079,bz079,lin)*2.
       q_bf = q_bf + (1./6.)*dt*temp
    end if
  end subroutine bootstrap_axial_field

  subroutine bootstrap_pressure(trial, lin, ssterm, ddterm, pp_g, thimpi)
    use basic
    use arrays
    use m3dc1_nint

    implicit none

    vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: trial, lin 
    vectype, dimension(num_fields), intent(inout) :: ssterm, ddterm
    integer, intent(in) :: pp_g
    real, intent(in) :: thimpi

    vectype :: temp

    temp = bs_b3pe(trial,lin)
    ssterm(pp_g) = ssterm(pp_g) -       thimpi     *dt*temp
    ddterm(pp_g) = ddterm(pp_g) + (1. - thimpi*bdf)*dt*temp
  end subroutine bootstrap_pressure

  
  ! B1psipsib
  ! =========
  vectype function bs_b1psipsib(e,f,g,h)
    use basic
    use m3dc1_nint

    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
    vectype :: temp

    temp = 0.

#if defined(USECOMPLEX) || defined(USE3D)
    if(jadv.eq.0) then
       temp = 0.
    else 
       if(surface_int) then
          temp = 0.
       else
          temp79a = eta79(:,OP_1)*h(:,OP_1)* &
               (pt79(:,OP_DZ)*f(:,OP_DZ) + pt79(:,OP_DR)*f(:,OP_DR))

          temp = int4(ri3_79,e(:,OP_DRP),temp79a,g(:,OP_DZ)) &
               - int4(ri3_79,e(:,OP_DZP),temp79a,g(:,OP_DR))

#ifdef USECOMPLEX
          temp = temp - rfac* &
               (int4(ri3_79,temp79a,e(:,OP_DR),g(:,OP_DZ)) &
               -int4(ri3_79,temp79a,e(:,OP_DZ),g(:,OP_DR))) 
#endif

       endif
    endif
#endif


    bs_b1psipsib = bootstrap_alpha*temp
  end function bs_b1psipsib

 
  ! B1psibb
  ! =======
  vectype function bs_b1psibb(e,f,g,h)
    use basic
    use m3dc1_nint

    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
    vectype :: temp

    if(jadv.eq.0) then
       if(surface_int) then
          temp = 0.
       else
          temp79a = eta79(:,OP_1)*g(:,OP_1)* &
               (pt79(:,OP_DZ)*f(:,OP_DZ) + pt79(:,OP_DR)*f(:,OP_DR))
          temp = int3(e(:,OP_1),temp79a,h(:,OP_1))
       end if
    else 
       if(surface_int) then
          temp = 0.
       else
          temp79a = eta79(:,OP_1)*g(:,OP_1)* &
               (pt79(:,OP_DZ)*f(:,OP_DZ) + pt79(:,OP_DR)*f(:,OP_DR))

          temp = int4(ri2_79,e(:,OP_GS),temp79a,h(:,OP_1))

          if(itor.eq.1) then
             temp = temp - 2.*int4(ri3_79,e(:,OP_DR),temp79a,h(:,OP_1))
          end if
       endif
    endif

    bs_b1psibb = bootstrap_alpha*temp
  end function bs_b1psibb


  ! B1psibf
  ! =======
  vectype function bs_b1psibf(e,f,g,h)
    use basic
    use m3dc1_nint

    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
    vectype :: temp

    temp = 0.

#ifdef defined(USECOMPLEX) || defined(USE3D)
    if(jadv.eq.0) then
       temp = 0.
    else 
       if(surface_int) then
          temp = 0.
       else
          temp79a = eta79(:,OP_1)*g(:,OP_1)* &
               (pt79(:,OP_DZ)*f(:,OP_DZ) + pt79(:,OP_DR)*f(:,OP_DR))

          temp = int4(ri2_79,e(:,OP_DRP),temp79a,h(:,OP_DRP) &
               + int4(ri2_79,e(:,OP_DZP),temp79a,h(:,OP_DZP))

#ifdef USECOMPLEX
          temp = temp - rfac* &
               (int4(ri2_79,e(:,OP_DR),temp79a,h(:,OP_DRP)) &
               +int4(rie_79,e(:,OP_DZ),temp79a,h(:,OP_DZP))
#endif

       endif
    endif
#endif

    bs_b1psibf = bootstrap_alpha*temp
  end function bs_b1psibf


  ! B2psipsib
  ! =========
  vectype function bs_b2psipsib(e,f,g,h)
    use basic
    use m3dc1_nint

    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
    vectype :: temp

    temp79a = eta79(:,OP_1)*h(:,OP_1)* &
         (pt79(:,OP_DZ)*f(:,OP_DZ) + pt79(:,OP_DR)*f(:,OP_DR))

    temp = int4(ri2_79,e(:,OP_DR),g(:,OP_DR),temp79a) &
         + int4(ri2_79,e(:,OP_DZ),g(:,OP_DZ),temp79a)

    bs_b2psipsib = bootstrap_alpha*temp
  end function bs_b2psipsib

  ! B2psibf
  ! =======
  vectype function bs_b2psibf(e,f,g,h)
    use basic
    use m3dc1_nint

    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f,g,h
    vectype :: temp

#ifdef defined(USECOMPLEX) || defined(USE3D)
    temp79a = eta79(:,OP_1)*g(:,OP_1)* &
         (pt79(:,OP_DZ)*f(:,OP_DZ) + pt79(:,OP_DR)*f(:,OP_DR))


    temp = int4(ri3_79,e(:,OP_DZ),h(:,OP_DRP),temp79a) &
         - int4(ri3_79,e(:,OP_DR),h(:,OP_DZP),temp79a)
#else 
    temp = 0.
#endif

    bs_b2psibf = bootstrap_alpha*temp
  end function bs_b2psibf

  ! B3pe
  ! ====
  vectype function bs_b3pe(e,f)
    use basic
    use m3dc1_nint

    implicit none

    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: e,f
    vectype :: temp

    temp79a = e(:,OP_1)*eta79(:,OP_1)*bzt79(:,OP_1)* &
         (f(:,OP_DZ)*pst79(:,OP_DZ) + f(:,OP_DR)*pst79(:,OP_DR))

    ! J.B
    temp79b = -ri2_79*bzt79(:,OP_1)*pst79(:,OP_GS) &
         + ri2_79*(bzt79(:,OP_DZ)*pst79(:,OP_DZ) + bzt79(:,OP_DR)*pst79(:,OP_DR))
#if defined(USE3D) || defined(USECOMPLEX)
    temp79b = temp79b &
         + ri2_79*(bft79(:,OP_DZPP)*pst79(:,OP_DZ)  + bft79(:,OP_DRPP)*pst79(:,OP_DR)) &
         + ri_79 *(bzt79(:,OP_DZ)*  bft79(:,OP_DRP) - bzt79(:,OP_DR)*  bft79(:,OP_DZP)) &
         + ri_79 *(bft79(:,OP_DZPP)*bft79(:,OP_DRP) - bft79(:,OP_DRPP)*bft79(:,OP_DZP)) &
         - ri3_79*(pst79(:,OP_DZ)*  pst79(:,OP_DRP) - pst79(:,OP_DR)*  pst79(:,OP_DZP)) &
         - ri2_79*(bft79(:,OP_DZP)* pst79(:,OP_DZP) + bft79(:,OP_DRP)* pst79(:,OP_DRP))
#endif

    temp = int2(temp79a, temp79b)

    bs_b3pe = -(gam-1.)*bootstrap_alpha*temp
  end function bs_b3pe


end module bootstrap
