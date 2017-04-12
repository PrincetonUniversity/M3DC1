!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
!> Module for estimating the runaway current
!<
module runaway_mod

  implicit none
  
  private
  real*8, parameter :: pi = 4.D0*atan(1.D0) 
  real*8, parameter :: eps0 = 8.854187817D-12 ![F/m]
  real*8, parameter :: c = 299792458.D0 ![m/s]
  real*8, parameter :: ec = 1.60217662D-19 ![C]
  real*8, parameter :: hbar = 1.05457180D-34
  real*8, parameter :: me = 9.10938356D-31 ![kg]
  
  type runaway_data
     real*8 :: rdens ! runaway density [1/m^3]
     real*8 :: jsign ! sign of the parallel current
  end type runaway_data

  type(runaway_data), allocatable, dimension(:) :: runaway_status

  public :: runaway_test
  public :: runaway_deallocate
  public :: runaway_init
  public :: runaway_current
  
contains

  subroutine runaway_deallocate()
    if(allocated(runaway_status)) deallocate(runaway_status)
  end subroutine runaway_deallocate
    
  subroutine runaway_init(dim)
    integer, intent(in) :: dim
    integer :: i
    call runaway_deallocate()
    allocate(runaway_status(dim))
    do i=1,dim
       runaway_status(i)%rdens=0.D0
    end do
  end subroutine runaway_init
  
  subroutine runaway_current(dt,dim,ER,EP,EZ,BR,BP,BZ,Temp,Dens,Zeff,JR,JP,JZ)
    integer, intent(in) :: dim
    real*8, intent(in), dimension(dim) :: ER,EP,EZ,BR,BP,BZ
    real*8, intent(in), dimension(dim) :: Temp ! [eV]
    real*8, intent(in), dimension(dim) :: Dens ! [1/m^3]
    real*8, intent(in), dimension(dim) :: Zeff ! [1]
    real*8, intent(in) :: dt   ! [s]
    real*8, intent(out), dimension(dim) :: JR,JP,JZ
    integer :: i
    real*8 :: Bmagn,dndt, Epar,Clog,x,Ecrit,nu,vth
    do i=1,dim
       Bmagn=sqrt(BR(i)**2+BP(i)**2+BZ(i)**2)
       Epar=(ER(i)*BR(i)+EP(i)*BP(i)+EZ(i)*BZ(i))/Bmagn
       runaway_status(i)%jsign=sign(1.D0,Epar)
       Epar=abs(Epar)
       Clog=14.9D0-0.5*log(dens(i)/1.d20)+log(temp(i)/1.d3)
       Ecrit = ec**3*Dens(i)*Clog/(4*pi*eps0**2*me*c**2) 
       vth=sqrt(2*ec*Temp(i)/me)
       nu=dens(i)*ec**4*clog/(4*pi*eps0**2*me**2*vth**3) 
       x=(Epar*ec*temp(i))/(Ecrit*me*c**2)
       if(Epar>Ecrit) then
          dndt=dens(i)*nu*x**(-3.D0*(1.D0+Zeff(i))/1.6D1)*exp(-1.D0/(4*x)-sqrt((1.D0+Zeff(i))/x))
       else
          dndt=0.D0
       end if
       runaway_status(i)%rdens=runaway_status(i)%rdens+dt*dndt
       JR(i)=-runaway_status(i)%jsign*BR(i)*ec*runaway_status(i)%rdens/Bmagn
       JP(i)=-runaway_status(i)%jsign*BP(i)*ec*runaway_status(i)%rdens/Bmagn
       JZ(i)=-runaway_status(i)%jsign*BZ(i)*ec*runaway_status(i)%rdens/Bmagn
    end do
  end subroutine runaway_current

  subroutine runaway_test()
    integer, parameter :: dim=1
    real*8, dimension(dim) :: ER,EP,EZ,BR,BP,BZ,JR,JP,JZ,Temp,Dens,Zeff
    real*8 :: dt
    integer :: nt,it
    ER=0.D0;    EP=5.D-2;    EZ=0.D0
    BR=0.D0;    BP=1.D0;    BZ=0.D0
    Temp=1000.D0
    Dens=1.D19
    Zeff=1.D0
    JR=0.D0
    JP=0.D0
    JZ=0.D0
    nt=10
    dt=1.D-2
    call runaway_init(dim)
    do it=1,nt
       call runaway_current(dt,dim,ER,EP,EZ,BR,BP,BZ,Temp,Dens,Zeff,JR,JP,JZ)
       write(*,*) 'current: ',JR,JP,JZ, ' [A/m^2]' 
    end do
    call runaway_deallocate()
  end subroutine runaway_test

end module runaway_mod
     

