module runaway_mod

  use basic
  use arrays
  use m3dc1_nint
  use field

  implicit none

  type(field_type), private :: dnre_field1,dndt_field,f_field

  type(field_type), private :: dnre_field2,jre_field,ere_field

  type(field_type), private :: depar_field,ecrit_field
    
  real, private, parameter :: eps0 = 8.854187817D-12 ![F/m]
  real, private, parameter :: c = 299792458.D0 ![m/s]
  real, private, parameter :: ec = 1.60217662D-19 ![C]
  real, private, parameter :: hbar = 1.05457180D-34
  real, private, parameter :: me = 9.10938356D-31 ![kg]
  real, private, parameter :: va = 1542350 ![m/s]

  real, dimension(MAX_PTS) :: re_j79

contains

  subroutine runaway_deallocate()
    implicit none
    if(irunaway.eq.0) return

    call destroy_field(dnre_field1)
    call destroy_field(dnre_field2)
    call destroy_field(depar_field)
    call destroy_field(ere_field)
    call destroy_field(ecrit_field)
    call destroy_field(dndt_field)
    call destroy_field(jre_field)
    call destroy_field(f_field)
  end subroutine runaway_deallocate
    
  subroutine runaway_init()
    use basic
    use math
    
    implicit none
    if(irunaway.eq.0) return

    print *, 'Estimated Ecrit for runaways = ', &
         ec**3*(n0_norm*1e6)*17./(4.*pi*eps0**2*me*c**2), ' V/m'

    call create_field(dnre_field1)
    call create_field(dnre_field2)
    call create_field(depar_field)
    call create_field(ere_field)
    call create_field(ecrit_field)
    call create_field(dndt_field)
    call create_field(jre_field)
    call create_field(f_field)

    dnre_field1 = 0.
    dnre_field2 = 0.
    depar_field = 0.
    ere_field = 0.
    ecrit_field = 0.
    dndt_field = 0.
    jre_field = 0.
    f_field = 0.
  end subroutine runaway_init
  
  elemental subroutine runaway_current(nre,epar,Temp,Dens,&
                                       Zeff,eta,Ecrit,re_79,&
                                       re_j79,re_epar,dndt,&
                                       mr,bz,bi,ri,z)
    use math
    use basic

    implicit none

    real, intent(in) :: nre  ! [1/m^3]
    real, intent(in) :: epar,eta
    real, intent(in) :: Temp ! [eV]
    real, intent(in) :: Dens ! [1/m^3]
    real, intent(in) :: Zeff ! [1]
    real, intent(in) :: bz
    real, intent(in) :: ri,z
    real, intent(in) :: bi
    real, intent(out) :: re_79,dndt,re_j79,re_epar,Ecrit
    real :: Clog,x,nu,vth,esign,teval,jpar,nrel,a,r,Ed, &
            f,dt_si,tmp,nretmp,Dens1,sd,sa,nra,gamma,tau
    integer, intent(in) :: mr
    integer :: l, nl
    
    dndt = 0.
    sd = 0.
    sa = 0.
    dt_si = dt * t0_norm
    nrel = nre
    nl = 10
    tmp = 0.
    nretmp = 0.
    Dens1 = Dens 
    
    if(mr.ne.0) then 
      re_79 = 0.D0
      re_j79 = 0.D0 
      return
    endif
             
!   based on formula in Stahl, et al, PRL 114 115002 (2015)
    teval = max(1.,Temp)   ! note: this sets a minimum temperature of 1 eV for runaway production
    esign = sign(1., -epar)
       jpar = nre!*ec*cre*va
       nra = nre/ec/cre/va 
       re_epar = epar! - 1.0*abs(eta*jpar/b1/bz/ri)
       if (Temp.ge.teval .and. mr.eq.0) then
          Clog = 14.78D0-0.5*log(Dens1/1.d20)+log(Temp/1.d3)
          Ecrit = ec**3*Dens1*Clog/(4*pi*eps0**2*me*c**2) 
          vth = sqrt(2*ec*Temp/me)
          nu = Dens1*ec**4*Clog/(4*pi*eps0**2*me**2*vth**3)
          tau = me*c/ec/Ecrit
          a = sqrt((1/ri-xmag)**2+(z-zmag)**2)
          r = sqrt(1/ri**2+z**2)
          gamma = 1/(1+1.46*sqrt(a/r)+1.72*a/r) 
          x = (abs(re_epar)*ec*Temp)/(Ecrit*me*c**2)
          Ed = abs(re_epar)/Ecrit
          if(Ed < 1) Ed = 1
          if(abs(re_epar).gt.Ecrit) then
              sd = Dens1*nu*x**(-3.D0*(1.D0+Zeff)/1.6D1) &
                 *exp(-1.D0/(4*x)-sqrt((1.D0+Zeff)/x))
              !write(0,*) re_epar, Ecrit
              if (re_epar*nre<0) then
                 sa = nra/tau/Clog*sqrt(pi*gamma/3/(Zeff+5))*(Ed-1)* &
                    1/sqrt(1-1/Ed+4*pi*(Zeff+1)**2/3/gamma/(Zeff+5)/(Ed**2+4/gamma**2-1)) 
              else
                 sa = 0.
              endif
              dndt = (sd*esign + sa)*cre*ec*va
                 nrel = nre + dndt*dt_si
          else
              nrel = nre 
              dndt = 0.
          end if

       else
           nrel = nre
       endif
    
       re_79 = nrel
       re_j79 = nre
  
       if (abs(re_j79) .ge. abs(ri*bz)) then
            dndt = 0.
       endif

  end subroutine runaway_current

  subroutine eval_runaway(itri,izone)
    use basic
    use m3dc1_nint
    use electric_field
    use diagnostics
    use math
    implicit none

    integer, intent(in) :: itri,izone
    vectype, dimension(MAX_PTS) :: epar, te, ne, nre, eta
    vectype, dimension(MAX_PTS) :: dndt, re_79, re_j79, &
                                   re_epar, ecrit, bz, bi  
    vectype, dimension(dofs_per_element) :: dofs
    integer, dimension(MAX_PTS) :: mr, tmp


    if(irunaway.eq.0 .or. izone.ne.1) return
#ifdef USECOMPLEX
    re_79 = 0.
    re_j79 = 0.
    return
#else
    call electric_field_par(0, epar, izone)

    dofs = intx2(mu79(:,:,OP_1),epar)
    call vector_insert_block(ere_field%vec,itri,1,dofs,VEC_ADD)
    
    ! convert to SI units
    epar = epar*e0_norm*c*1e-4
    te = tet79(:,OP_1)*(p0_norm/n0_norm)/1.6022e-12
    ne = net79(:,OP_1)*n0_norm*1e6
    nre = nre179(:,OP_1)*j0_norm/(c*1e-3)!n0_norm*1e6
    eta = eta79(:,OP_1)*e0_norm/j0_norm*c**2*1e-7
    call magnetic_region(pst79(:,OP_1),&
         pst79(:,OP_DR),pst79(:,OP_DZ),x_79,Z_79,mr)
    bz = pst79(:,OP_GS)/((c*1e-3)/j0_norm)
    bi = bi79(:,OP_1)
    
    call runaway_current(nre,epar,te,ne,z_ion,eta,ecrit,&
                         re_79,re_j79,re_epar,dndt,mr,bz,bi,ri_79,z_79)

    ! convert back to normalized units
    dndt = dndt*t0_norm/(j0_norm/c/1e-3)

    re_79 = re_79*(c*1e-3)/j0_norm!/(n0_norm*1e6)

    re_j79 = re_j79*(c*1e-3)/j0_norm
   
    re_j79 = re_j79+dndt*dt

    re_epar = re_epar/(e0_norm*c*1e-4)

    ecrit = ecrit/(e0_norm*c*1e-4)
    
    dofs = intx2(mu79(:,:,OP_1),dndt) 
    call vector_insert_block(dndt_field%vec,itri,1,dofs,VEC_ADD)

    dofs = intx2(mu79(:,:,OP_1),re_79)
    call vector_insert_block(dnre_field1%vec,itri,1,dofs,VEC_ADD)

    dofs = intx2(mu79(:,:,OP_1),re_epar)
    call vector_insert_block(depar_field%vec,itri,1,dofs,VEC_ADD)

    dofs = intx2(mu79(:,:,OP_1),ecrit)
    call vector_insert_block(ecrit_field%vec,itri,1,dofs,VEC_ADD)

    dofs = intx2(mu79(:,:,OP_1),re_j79)
    call vector_insert_block(jre_field%vec,itri,1,dofs,VEC_ADD)



#endif

  end subroutine eval_runaway

  subroutine runaway_advance
    use basic
    use newvar_mod
    use m3dc1_nint
    implicit none

    integer :: ier

    if(irunaway.eq.0) return

    call newvar_solve(jre_field%vec, mass_mat_lhs)

    nre_field(1) = jre_field
    dnre_field2 = 0.
    dnre_field1 = 0.
    jre_field = 0.

  end subroutine runaway_advance

  subroutine smooth_runaway
    use basic
    use m3dc1_nint
    use newvar_mod
    use diagnostics
    use math
    implicit none

    vectype, dimension(MAX_PTS) :: jre
    vectype, dimension(MAX_PTS, OP_NUM) :: jre79
    vectype, dimension(dofs_per_element) :: dofs
    integer :: numelms, itri
    numelms = local_elements()    
    jre_field = 0.     
    do itri=1, numelms
       call define_element_quadrature(itri,int_pts_main,int_pts_tor)
       call define_fields(itri,0,1,0)
       call eval_ops(itri, nre_field(1), jre79,rfac)

       where(real(jre79(:,OP_1)) .ge. 0.0)
               jre79(:,OP_1) = 0.0
       endwhere

       jre = jre79(:,OP_1)
       dofs = intx2(mu79(:,:,OP_1),jre)
       call vector_insert_block(jre_field%vec,itri,1,dofs,VEC_ADD)
    enddo

    call newvar_solve(jre_field%vec, mass_mat_lhs)
    nre_field(1) = jre_field
    jre_field= 0.




  end subroutine



     
end module runaway_mod
