program read_test

  implicit none

  integer :: i, itime, ierr
  real(8) :: R, Phi, Z
  real(8) :: Br, Bphi, Bz

  print *, 'loading...'
  ! time = 1 for total field, 0 for vacuum field, -1 for equilibrium only 
  itime = 1
  call m3dc1_load_file(itime, ierr)

  if(ierr.ne.0) then
     print *, 'Error loding C1.h5 file'
     return
  end if

  print *, 'evaluating...'
  Phi = 0.
  Z = 0.
  do i=1, 10
     R = 1.1 + i/10.
     call m3dc1_get_field(R, Phi, Z, Br, Bphi, Bz)
     write(*, '("Field at (", 3F6.2, ") is ", 3F12.4)') R, Phi, Z, Br, Bphi, Bz
  end do

  print *, 'unloading..'
  call m3dc1_unload_file()

end program read_test
