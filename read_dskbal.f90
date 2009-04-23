module dskbal

  implicit none

  character*80 :: heading
  character*80, parameter :: heading1 = '&end'

  integer :: npsi_bal, nthet_bal

  real, allocatable :: psi_bal(:), pprime_bal(:), f_bal(:), ffprime_bal(:)
  real, allocatable :: chipsi_bal(:), q_bal(:), x_bal(:,:), z_bal(:,:)
  real, allocatable :: ne_bal(:), neprime_bal(:)
  real, allocatable :: te_bal(:), teprime_bal(:), ti_bal(:), tiprime_bal(:)
  real, allocatable :: p_bal(:)

contains

subroutine load_dskbal
  
  implicit none

  integer, parameter :: ifile = 32

  integer :: i, j

  open(unit=ifile,file='dskbal',status='old')

  do i=1,100
     read(ifile,*,err=100,end=100) heading
     if(heading.eq.heading1) exit
  end do

  read(ifile,*) heading
  read(ifile,'(2I5)') npsi_bal, nthet_bal

  allocate(psi_bal(npsi_bal), pprime_bal(npsi_bal), f_bal(npsi_bal), &
       ffprime_bal(npsi_bal), chipsi_bal(npsi_bal), q_bal(npsi_bal), &
       x_bal(nthet_bal,npsi_bal), z_bal(nthet_bal,npsi_bal), &
       ne_bal(npsi_bal), neprime_bal(npsi_bal), te_bal(npsi_bal), &
       teprime_bal(npsi_bal), ti_bal(npsi_bal), tiprime_bal(npsi_bal), &
       p_bal(npsi_bal))

  read(ifile,*) heading
  read(ifile,'(1p5e20.12)') (psi_bal(j),j=1,npsi_bal)

  read(ifile,*) heading
  read(ifile,'(1p5e20.12)') (pprime_bal(j),j=1,npsi_bal)

  read(ifile,*) heading
  read(ifile,'(1p5e20.12)') (f_bal(j),j=1,npsi_bal)

  read(ifile,*) heading
  read(ifile,'(1p5e20.12)') (ffprime_bal(j),j=1,npsi_bal)

  read(ifile,*) heading
  read(ifile,'(1p5e20.12)') (chipsi_bal(j),j=1,npsi_bal)

  read(ifile,*) heading
  read(ifile,'(1p5e20.12)') (q_bal(j),j=1,npsi_bal)

  read(ifile,*) heading
  do i=1,nthet_bal
     read(ifile,'(1p5e20.12)') (x_bal(i,j),j=1,npsi_bal)
  end do

  read(ifile,*) heading
  do i=1,nthet_bal
     read(ifile,'(1p5e20.12)') (z_bal(i,j),j=1,npsi_bal)
  end do

  read(ifile,*) heading
  read(ifile,'(1p5e20.12)') (ne_bal(j),j=1,npsi_bal)

  read(ifile,*) heading
  read(ifile,'(1p5e20.12)') (neprime_bal(j),j=1,npsi_bal)

  read(ifile,*) heading
  read(ifile,'(1p5e20.12)') (te_bal(j),j=1,npsi_bal)

  read(ifile,*) heading
  read(ifile,'(1p5e20.12)') (teprime_bal(j),j=1,npsi_bal)

  read(ifile,*) heading
  read(ifile,'(1p5e20.12)') (ti_bal(j),j=1,npsi_bal)

  read(ifile,*) heading
  read(ifile,'(1p5e20.12)') (tiprime_bal(j),j=1,npsi_bal)

  ! p is not provided; calculate from density and temperature profiles
  p_bal = ne_bal*(te_bal + ti_bal)*1.6022e-12

  ! convert to SI
  x_bal = x_bal / 100.
  z_bal = z_bal / 100.

  psi_bal = psi_bal / 1e8
  f_bal = f_bal / 1e6
  p_bal = p_bal / 10.
  ffprime_bal = ffprime_bal / 1e6 / 1e6 * 1e8
  pprime_bal = pprime_bal / 10. * 1e8
  

  goto 101
100 print *, 'Error reading file.'

101 close(ifile)

end subroutine load_dskbal

subroutine unload_dskbal
  implicit none

  deallocate(psi_bal, pprime_bal, f_bal, ffprime_bal, &
       chipsi_bal, q_bal, x_bal, z_bal, &
       ne_bal, neprime_bal, te_bal, teprime_bal, &
       ti_bal, tiprime_bal, p_bal)
end subroutine unload_dskbal

end module dskbal
