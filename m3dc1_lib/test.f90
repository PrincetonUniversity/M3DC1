program m3dc1_fortran_test

  implicit none

  integer :: i, ierr
  real :: r, phi, z, br, bphi, bz, n, p
  integer :: n_handle, p_handle
  integer :: timeslice = 1

  ! Open the file
  ! (note that we're calling a c function 
  !  so the string must be null-terminated)
  call m3dc1_open_file("C1.h5"//char(0), ierr)
  if(ierr .ne. 0) then
     print *, 'Error opening file.'
     stop
  end if

  call m3dc1_load_magnetic_field(timeslice, ierr)
  if(ierr .ne. 0) then
     print *, 'Error loading magnetic field.'
     goto 100
  end if

  call m3dc1_load_field("den"//char(0), timeslice, n_handle, ierr)
  if(ierr.ne.0) then
     print *, 'Error loading density field.'
     goto 100
  end if

  call m3dc1_load_field("P"//char(0), timeslice, p_handle, ierr)
  if(ierr.ne.0) then
     print *, 'Error loading pressure field.'
     goto 100
  end if

  do i=1, 12
     r = 1.0 + 0.1*i
     phi = 0.
     z = 0.
     call m3dc1_eval_magnetic_field(r, phi, z, br, bphi, bz, ierr)
     if(ierr .ne. 0) then
        print *, 'Error evaluating magnetic field.', ierr
     end if
     
     call m3dc1_eval_field(n_handle, r, phi, z, n, ierr)
     if(ierr .ne. 0) then 
        print *, 'Error evaluating density field.'
     end if

     call m3dc1_eval_field(p_handle, r, phi, z, p, ierr)
     if(ierr .ne. 0) then 
        print *, 'Error evaluating pressure field.'
     end if

     write(*, '(5f12.6)') br, bphi, bz, n, p
  end do
  
100 continue
  call m3dc1_close_file()

end program m3dc1_fortran_test
