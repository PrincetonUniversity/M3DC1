program m3dc1_fortran_test

  implicit none

  integer :: i, ierr
  real :: r, phi, z, br, bphi, bz, n, p
  integer :: n_handle, p_handle

  ! Open the file
  ! (note that we're calling a c function 
  !  so the string must be null-terminated)
  call m3dc1_open_file("C1.h5"//char(0), ierr)
  if(ierr .ne. 0) then
     print *, 'Error opening file.'
     stop
  end if

  call m3dc1_load_magnetic_field(1, ierr)
  if(ierr .ne. 0) then
     print *, 'Error loading magnetic field.'
     goto 100
  end if

  call m3dc1_load_field("n"//char(0), 1, n_handle, ierr)
  if(ierr.ne.0) then
     print *, 'Error loading n field.'
     goto 10
  end if

  do i=1, 10
     r = 1.2 + 0.1*i
     phi = 0.
     z = 0.
     call m3dc1_eval_magnetic_field(r, phi, z, br, bphi, bz)
     if(ierr .ne. 0) then
        print *, 'Error evaluating magnetic field.'
     end if
     
     call m3dc1_eval_field(n_handle, r, phi, z, n, ierr)
     if(ierr .ne. 0) then 
        print *, 'Error evaluating n field.'
     end if

     write(*, '(3f12.6)') br, bphi, bz, n, p
  end do


10 continue
  call m3dc1_unload_magnetic_field(1, ierr)
  if(ierr .ne. 0) then
     print *, 'Error unloading magnetic field.'
     stop
  end if
  
100 continue
  call m3dc1_close_file()

end program m3dc1_fortran_test
