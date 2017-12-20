!======================================================================
! polar
!
! This module contains the subroutine write_polar
! which takes an array of x and z locations
! and writes a file in POLAR format.
!======================================================================
module polar

  implicit none

contains

  subroutine write_polar(x, z, m, n, itranspose)

    implicit none

    integer, intent(in) :: itranspose
    integer, intent(in) :: m,n
    real, intent(in), dimension(m,n) :: x, z
    
    integer, parameter :: ifile = 74
    integer :: i, j, nt, np
    real, allocatable :: norm(:,:), curv(:)

    if(itranspose.eq.0) then
       nt = m
       np = n
    else
       nt = n
       np = m
    endif

    print *, 'Writing POLAR file with '
    print *, ' ntheta, npsi = ', nt, np
    print *, ' itranspose = ', itranspose

    allocate(norm(2,nt),curv(nt))

    open(unit=ifile,file="POLAR",status="unknown")
    write(ifile,'(2i5)') np,nt

    if(itranspose.eq.0) then
       do j=1,np-1
          do i=1,nt
             write(ifile,'(1p2e18.10)') x(i,j),z(i,j)
          enddo
       enddo
       call calc_norms(x(:,np),z(:,np),nt,norm,curv)
       do i=1,nt
          write(ifile,'(1p5e18.10)') &
               x(i,np),z(i,np),norm(1,i),norm(2,i),curv(i)
       end do
    else
       do j=1,np-1
          do i=1,nt
             write(ifile,'(1p2e18.10)') x(j,i),z(j,i)
          enddo
       enddo
       call calc_norms(x(np,:),z(np,:),nt,norm,curv)
       do i=1,nt
          write(ifile,'(1p5e18.10)') &
               x(np,i),z(np,i),norm(1,i),norm(2,i),curv(i)
       end do
    end if
    
    close(ifile)

    deallocate(norm, curv)
  
  end subroutine write_polar

subroutine calc_norms(x, z, n, norm, curv)

  implicit none

  integer, intent(in) :: n
  real, intent(in), dimension(n) :: x, z
  real, intent(out), dimension(2,n) :: norm
  real, intent(out), dimension(n) :: curv

  integer :: i, i1, i2
  real :: dx1, dx2, dz1, dz2, l, t1, t2, t3
  real :: l1(n), l2(n), norm1(2), norm2(2)

  ! calculate normals
  do i=1,n
     i1 = mod(i+n-2,n)+1
     i2 = mod(i,n)+1
     dx1 = x(i) - x(i1)
     dx2 = x(i2) - x(i)
     dz1 = z(i) - z(i1)
     dz2 = z(i2) - z(i)

     l1(i) = 1./sqrt(dx1**2 + dz1**2)
     l2(i) = 1./sqrt(dx2**2 + dz2**2)
     norm1(1) =  dz1*l1(i)
     norm1(2) = -dx1*l1(i)
     norm2(1) =  dz2*l2(i)
     norm2(2) = -dx2*l2(i)

     ! perform weigted average of adjacent edge normals
     norm(:,i) = (l1(i)*norm1 + l2(i)*norm2)/(l1(i)+l2(i))

     ! normalize normal
     l = sqrt(norm(1,i)**2 + norm(2,i)**2)
     norm(:,i) = norm(:,i)/l
  enddo

  ! calculate curvature
  do i=1,n
     i1 = mod(i+n-2,n)+1
     i2 = mod(i,n)+1

     t1 = atan2(norm(2,i1),norm(1,i1))
     t2 = atan2(norm(2,i ),norm(1,i ))
     t3 = atan2(norm(2,i2),norm(1,i2))
     
     ! perform weigted average of dtheta/dl
     curv(i) = ((t2-t1)*l1(i)**2 + (t3-t2)*l2(i)**2)/(l1(i) + l2(i))
  end do
  

end subroutine calc_norms


end module polar
