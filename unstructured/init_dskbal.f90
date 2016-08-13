!==============================================================================
! Dskbal_eq
! ~~~~~~~~~
!
! Loads dskbal equilibrium
!==============================================================================
module dskbal_eq

  implicit none

contains

subroutine dskbal_init()
  use math
  use basic
  use arrays
  use dskbal
  use gradshafranov
  use newvar_mod
  use sparse
  use diagnostics

  implicit none

  integer :: i, inode, numnodes, icounter_tt
  real :: dp, a(4), minden, minte, minti

  print *, "dskbal_init called"

  call load_dskbal
  p_bal = p_bal*amu0
  pprime_bal = pprime_bal*amu0

  ! find "vacuum" values
  minden = 0.
  minte = 0.
  minti = 0.
  do i=1,npsi_bal
     if(ne_bal(i) .le. 0.) exit
     minden = ne_bal(i)
     minte = te_bal(i)
     minti = ti_bal(i)
  end do

  ! replace zeroed-out values with "vacuum" values
  do i=1,npsi_bal
     if(ne_bal(i) .le. 0.) ne_bal(i) = minden
     if(ti_bal(i) .le. 0.) ti_bal(i) = minti
     if(te_bal(i) .le. 0.) te_bal(i) = minte
     if(p_bal(i) .le. 0.) p_bal(i) = minden*(minte + minti)*1.6022e-12*amu0/10.
  end do

  xmag = x_bal(1,1)
  zmag = z_bal(1,1)
  rzero = xmag
  bzero = f_bal(npsi_bal)/rzero

  ifixedb = 1

  if(iread_dskbal.eq.2) then
     call default_profiles
  else
     call create_profile(npsi_bal,p_bal,pprime_bal,f_bal,ffprime_bal,psi_bal)
  end if
  
  ! initial plasma current filament
  call deltafun(xmag,zmag,tcuro,jphi_field)

  call gradshafranov_solve
  call gradshafranov_per  

  ! set density profile
  numnodes = owned_nodes()
  do icounter_tt=1,numnodes
     inode = nodes_owned(icounter_tt)
     call get_local_vals(inode)
     do i=1, npsi_bal-1
        if((psi_bal(i+1)-psi_bal(1))/(psi_bal(npsi_bal)-psi_bal(1)) &
             .gt. (real(psi0_l(1))-psimin)/(psilim-psimin)) exit
     end do

     call cubic_interpolation_coeffs(ne_bal,npsi_bal,i,a)

     dp = ((real(psi0_l(1))-psimin) &
          *(psi_bal(npsi_bal) - psi_bal(1))/(psilim - psimin) &
          -(psi_bal(i)-psi_bal(1))) / (psi_bal(i+1) - psi_bal(i))

     den0_l(1) = a(1) + a(2)*dp + a(3)*dp**2 + a(4)*dp**3
     den0_l(2) = (a(2) + 2.*a(3)*dp + 3.*a(4)*dp**2)*psi0_l(2)
     den0_l(3) = (a(2) + 2.*a(3)*dp + 3.*a(4)*dp**2)*psi0_l(3)
     den0_l(4) = (2.*a(3) + 6.*a(4)*dp)*psi0_l(4)
     den0_l(5) = (2.*a(3) + 6.*a(4)*dp)*psi0_l(5)
     den0_l(6) = (2.*a(3) + 6.*a(4)*dp)*psi0_l(6)
     den0_l = den0_l / n0_norm
     call set_local_vals(inode)
  enddo

  call unload_dskbal

end subroutine dskbal_init
  
end module dskbal_eq
