module jsolver

  implicit none

  real, allocatable :: p_jsv(:), ppxx_jsv(:), q_jsv(:), qp_jsv(:)
  real, allocatable :: gxx_jsv(:), gpx_jsv(:), fb_jsv(:), fbp_jsv(:)
  real, allocatable :: f_jsv(:), fp_jsv(:), psival_jsv(:)
  real, allocatable :: x_jsv(:,:), z_jsv(:,:), aj3_jsv(:,:), aj_jsv(:,:)
  real, allocatable :: te_jsv(:), ane_jsv(:), ti_jsv(:)
  real, allocatable :: avez_jsv(:), zeffa_jsv(:)

  integer :: nz1, nz2, nthe, npsi_jsv, isym
  real :: dr, dt, xzero_jsv, p0_jsv, plimst, pminst

contains

subroutine load_jsolver
  
  implicit none

  integer, parameter :: ifile = 33
  integer :: i,j

  open(unit=ifile,file='fixed',status='old')

  read(ifile,1098) nz1,nz2,nthe,npsi_jsv,isym,dr,dt
  read(ifile,1099) xzero_jsv,p0_jsv,plimst,pminst

1098 format(5i5,1p2e12.4)
1099 format(1p5e16.8)

  allocate(p_jsv(nz2), ppxx_jsv(nz2), q_jsv(nz2), qp_jsv(nz2), &
       gxx_jsv(nz2), gpx_jsv(nz2), fb_jsv(nz2), fbp_jsv(nz2), &
       f_jsv(nz2), fp_jsv(nz2), psival_jsv(nz2), &
       x_jsv(nz1,nz2), z_jsv(nz1,nz2), aj3_jsv(nz1,nz2), aj_jsv(nz1,nz2), &
       te_jsv(nz2), ane_jsv(nz2), ti_jsv(nz2), &
       avez_jsv(nz2), zeffa_jsv(nz2))

  read(ifile,1099) (p_jsv(j),j=1,nz2)
  read(ifile,1099) (ppxx_jsv(j),j=1,nz2)
  read(ifile,1099) (q_jsv(j),j=1,nz2)
  read(ifile,1099) (qp_jsv(j),j=1,nz2)
  read(ifile,1099) (gxx_jsv(j),j=1,nz2)
  read(ifile,1099) (gpx_jsv(j),j=1,nz2)
  read(ifile,1099) (fb_jsv(j),j=1,nz2)
  read(ifile,1099) (fbp_jsv(j),j=1,nz2)
  read(ifile,1099) (f_jsv(j),j=1,nz2)
  read(ifile,1099) (fp_jsv(j),j=1,nz2)
  read(ifile,1099) (psival_jsv(j),j=1,nz2)
  read(ifile,1099) ((x_jsv(i,j),i=1,nz1),j=1,nz2)
  read(ifile,1099) ((z_jsv(i,j),i=1,nz1),j=1,nz2)
  read(ifile,1099) ((aj3_jsv(i,j),i=1,nz1),j=1,nz2)
  read(ifile,1099) ((aj_jsv(i,j),i=1,nz1),j=1,nz2)
  read(ifile,1099) (te_jsv(j),j=1,nz2)
  read(ifile,1099) (ane_jsv(j),j=1,nz2)
  read(ifile,1099) (ti_jsv(j),j=1,nz2)
  read(ifile,1099) (avez_jsv(j),j=1,nz2)
  read(ifile,1099) (zeffa_jsv(j),j=1,nz2)

  close(ifile)

end subroutine load_jsolver


subroutine unload_jsolver
  implicit none

  deallocate(p_jsv, ppxx_jsv, q_jsv, qp_jsv, &
       gxx_jsv, gpx_jsv, fb_jsv, fbp_jsv, &
       f_jsv, fp_jsv, psival_jsv, &
       x_jsv, z_jsv, aj3_jsv, aj_jsv, &
       te_jsv, ane_jsv, ti_jsv, &
       avez_jsv, zeffa_jsv)

end subroutine unload_jsolver


end module jsolver


#ifdef READ_JSOLVER

Program readjsolver

  use jsolver
  use polar

  implicit none

  call load_jsolver

  call write_polar(x_jsv(2:nthe,1:npsi_jsv), z_jsv(2:nthe,1:npsi_jsv), &
       nthe-1, npsi_jsv, 0)

  call unload_jsolver

end Program readjsolver

#endif

