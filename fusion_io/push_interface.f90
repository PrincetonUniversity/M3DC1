module fio_push
  use fusion_io

  implicit none

  type :: field_type
     real :: phi
     real, dimension(3) :: gradphi, a
     real, dimension(3,3) :: grada, ginv
     real, dimension(3,3,3) :: gmat1
  end type field_type

  integer, private :: isrc
  integer, private, dimension(2) :: iphi, ie, ia, igrada

contains

  subroutine fio_push_initialize(itype, filename, linfac)
    implicit none

    integer, intent(in) :: itype
    character(len=*), intent(in) :: filename
    real, intent(in) :: linfac

    integer :: ierr

    call fio_open_source_f(itype, trim(filename), isrc, ierr)
    if(ierr.ne.0) return

     ! Set options appropriate to this source
     call fio_get_options_f(isrc, ierr)
     call fio_set_real_option_f(FIO_LINEAR_SCALE, linfac, ierr)

     ! read vacuum fields
     call fio_set_int_option_f(FIO_TIMESLICE, 0, ierr)
     call fio_get_field_f(isrc, FIO_SCALAR_POTENTIAL, iphi(1), ierr);
     call fio_get_field_f(isrc, FIO_ELECTRIC_FIELD, ie(1), ierr);
     call fio_get_field_f(isrc, FIO_VECTOR_POTENTIAL, ia(1), ierr);
     call fio_get_field_f(isrc, FIO_GRAD_VECTOR_POTENTIAL, igrada(1), ierr);

     ! read plasma response fields
     call fio_set_int_option_f(FIO_TIMESLICE, 1, ierr)
     call fio_get_field_f(isrc, FIO_SCALAR_POTENTIAL, iphi(2), ierr);
     call fio_get_field_f(isrc, FIO_ELECTRIC_FIELD, ie(2), ierr);
     call fio_get_field_f(isrc, FIO_VECTOR_POTENTIAL, ia(2), ierr);
     call fio_get_field_f(isrc, FIO_GRAD_VECTOR_POTENTIAL, igrada(2), ierr);
   end subroutine fio_push_initialize


   subroutine fio_push_finalize()
     implicit none
     integer :: i, ierr
     
     do i=1, 2
        call fio_close_field_f(iphi(i), ierr)
        call fio_close_field_f(ie(i), ierr)
        call fio_close_field_f(ia(i), ierr)
        call fio_close_field_f(igrada(i), ierr)
     end do
   end subroutine fio_push_finalize
   

   subroutine fio_push_field_eval(t,q,field)
     implicit none

     real, intent(in) :: t
     real, dimension(3), intent(in) :: q
     type(field_type), intent(out) :: field

     integer :: it, ierr
     real :: phi(1)

     if(t.eq.0.) then
        it = 1
     else
        it = 2
     end if
     
     call fio_eval_field_f(iphi(it), q, phi, ierr)
     field%phi = phi(1)
     call fio_eval_field_f(ie(it), q, field%gradphi, ierr)
     field%gradphi = -field%gradphi
     call fio_eval_field_f(ia(it), q, field%a, ierr)
     call fio_eval_field_f(igrada(it), q, field%grada, ierr)

     ! Result is in (R, Phi, Z) coordinates
     ! metric tensor = diag(1, r^2, 1)

     ! inverse of metric tensor
     field%ginv = 0.
     field%ginv(1,1) = 1.
     field%ginv(2,2) = 1./q(1)**2
     field%ginv(3,3) = 1.
     
     ! gradient of metric tensor
     field%gmat1 = 0.
     field%gmat1(1,2,2) = 2.*q(1)    
   end subroutine fio_push_field_eval
end module fio_push
