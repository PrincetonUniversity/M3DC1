module runaway_mod

  use basic
  use arrays
  use m3dc1_nint
  use field

  implicit none

  type(field_type), private :: dnre_field
    
  real, private, parameter :: eps0 = 8.854187817D-12 ![F/m]
  real, private, parameter :: c = 299792458.D0 ![m/s]
  real, private, parameter :: ec = 1.60217662D-19 ![C]
  real, private, parameter :: hbar = 1.05457180D-34
  real, private, parameter :: me = 9.10938356D-31 ![kg]

  real, dimension(MAX_PTS) :: re_j79

contains

  subroutine runaway_deallocate()
    implicit none
    if(irunaway.eq.0) return

    call destroy_field(nre_field)
    call destroy_field(dnre_field)
  end subroutine runaway_deallocate
    
  subroutine runaway_init()
    use basic
    use math
    
    implicit none
    if(irunaway.eq.0) return

    print *, 'Estimated Ecrit for runaways = ', &
         ec**3*(n0_norm*1e6)*17./(4.*pi*eps0**2*me*c**2), ' V/m'

    call create_field(nre_field)
    call create_field(dnre_field)

    nre_field = 0.
    dnre_field = 0.
  end subroutine runaway_init
  
  elemental subroutine runaway_current(nre,epar,Temp,Dens,Zeff,jpar,dndt)
    use math

    implicit none

    real, intent(in) :: nre  ! [1/m^3]
    real, intent(in) :: epar
    real, intent(in) :: Temp ! [eV]
    real, intent(in) :: Dens ! [1/m^3]
    real, intent(in) :: Zeff ! [1]
    real, intent(out) :: jpar,dndt
    real :: Clog,x,Ecrit,nu,vth,jsign,teval


    teval = max(10.,Temp)
    jsign = sign(1.,Epar)
    Clog = 14.9D0-0.5*log(dens/1.d20)+log(teval/1.d3)
    Ecrit = ec**3*Dens*Clog/(4*pi*eps0**2*me*c**2) 
    vth = sqrt(2*ec*teval/me)
    nu = dens*ec**4*clog/(4*pi*eps0**2*me**2*vth**3) 
    x = (abs(Epar)*ec*teval)/(Ecrit*me*c**2)
    if(abs(Epar).gt.Ecrit) then
       dndt = dens*nu*x**(-3.D0*(1.D0+Zeff)/1.6D1) &
            *exp(-1.D0/(4*x)-sqrt((1.D0+Zeff)/x))
    else
       dndt = 0.D0
    end if
    jpar = -jsign*ec*nre
  end subroutine runaway_current

  subroutine eval_runaway(itri,izone)
    use basic
    use m3dc1_nint
    use electric_field
    implicit none

    integer, intent(in) :: itri,izone
    vectype, dimension(MAX_PTS) :: epar, te, ne, nre
    vectype, dimension(MAX_PTS) :: dndt
    vectype, dimension(dofs_per_element) :: dofs

    if(irunaway.eq.0 .or. izone.ne.1) return
#ifdef USECOMPLEX
    re_j79 = 0.
    return
#else
    call electric_field_par(0, epar, izone)
    
    ! convert to SI units
    epar = epar*e0_norm*c*1e-4
    te = tet79(:,OP_1)*(p0_norm/n0_norm)/1.6022e-12
    ne = net79(:,OP_1)*n0_norm*1e6
    nre = nre79(:,OP_1)*n0_norm*1e6

    call runaway_current(nre,epar,te,ne,z_ion,re_j79,dndt)

    ! convert back to normalized units
    dndt = dndt*t0_norm/(n0_norm*1e6)
    re_j79 = re_j79/(c*1e-3)/j0_norm

    dofs = intx2(mu79(:,:,OP_1),nre79(:,OP_1)) + &
         intx2(mu79(:,:,OP_1),dndt)*dt

    call vector_insert_block(dnre_field%vec,itri,1,dofs,VEC_SET)
#endif

  end subroutine eval_runaway

  subroutine runaway_advance
    use basic
    use newvar_mod

    implicit none

    integer :: ier

    if(irunaway.eq.0) return

    if(iprint.ge.1) print *, 'Solving RE advance'
    call sum_shared(dnre_field%vec)
    call newsolve(mass_mat_lhs%mat,dnre_field%vec,ier)
    nre_field = dnre_field
    dnre_field = 0.
  end subroutine runaway_advance

end module runaway_mod
     

