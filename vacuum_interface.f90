module vacuum_interface

  implicit none

  integer :: nodes
  integer, allocatable :: global_id(:)
  integer, allocatable :: local_id(:)
  complex, allocatable :: zgrbth(:,:), zgrbph(:,:), zgrbthp(:,:), zgrbphp(:,:)
  real, allocatable :: xnode(:), znode(:)

contains 

  subroutine load_boundary_nodes(ierr)
    implicit none
    
    integer, intent(out) :: ierr

    character(len=*), parameter :: filename = 'ordered.points'
    integer, parameter :: ifile = 1
    integer :: i
    real :: dummy

    ierr = 0

    open(unit=ifile,file=filename,action='read',status='unknown')

    ! read total number of nodes
    read(ifile, '(I8)') nodes

    allocate(global_id(nodes))
    allocate(local_id(nodes))
    allocate(xnode(nodes+1))
    allocate(znode(nodes+1))

    ! write node entry for first node
    do i=1, nodes
       read(ifile, '(I8,2f12.6)') global_id(i), xnode(i),znode(i)
       call globalidnod(global_id(i),local_id(i))
    end do

    close(ifile)
    
  end subroutine load_boundary_nodes


  subroutine load_response_matrix(ierr)
    use basic
    implicit none

    integer, intent(out) :: ierr
    
    character(len=*), parameter :: filename = 'RESPONSE-M3DC1'
    integer, parameter :: ifile = 2
    integer :: idum, i, j, m, mmax
    real :: dtheta, ka, thetai, thetaj, fac, bessk, besskp, grate, temp

    ierr = 0
!
!...check if analytic test problem
    if(itaylor.ne.10) then

    open(unit=ifile,file=filename,action='read',status='unknown')

    read(ifile, '(5x, "Response Matrix in theta space for M3D_C1" )')

    read(ifile, '(/,1x, "No of internal points used = ", i5 )') idum 
    write(*, '(/,1x, "No of internal points used = ", i5 )') idum 

    read(ifile, '(/,1x, "Min Fourier mode used = ", i5 )') idum
    write(*, '(/,1x, "Min Fourier mode used = ", i5 )') idum

    read(ifile, '(  1x, "Max Fourier mode used = ", i5 )' ) idum
    write(*, '(  1x, "Max Fourier mode used = ", i5 )' ) idum

    read(ifile, '(/,1x, "Toroidal mode no. = ", i5 )' ) idum
    write(*, '(/,1x, "Toroidal mode no. = ", i5 )' ) idum

    read(ifile, &
         '(/,1x, "Number of surface points on the closed domain = ", i5 )' ) &
         idum

    if(idum .ne. nodes+1) then 
       print *, 'Error: nodes in response file different from nodes in ordered.points', idum-1, nodes
       ierr = 1
       goto 100
    endif

    allocate(zgrbph (nodes+1,nodes+1),zgrbth (nodes+1,nodes+1))
    allocate(zgrbphp(nodes+1,nodes+1),zgrbthp(nodes+1,nodes+1))

    read(ifile, & 
         '(/,1x, "B_theta response Matrix, Rth(obs,srce):" )' )

    do i=1, nodes+1
       read(ifile, '(/, 1x, "i_obs = ", i5 )' ) idum
       read(ifile, '( (1x, 8es14.6) )' ) (zgrbth(idum,j), j=1, nodes+1)
    END DO
    
    read(ifile, '(/,1x, "IMAG. B_phi response Matrix, Rph(obs,src):" )')
    
    do i=1, nodes+1
       read(ifile, '(/, 1x, "i_obs = ", i5 )' ) idum
       read(ifile, '( (1x, 8es14.6) )' ) (zgrbph(idum,j), j=1, nodes+1)
    end do

    dtheta = 2.*3.14159265358979323846/nodes
!
!...NOTE:  equivalent to multiplying by dtheta and dividing by 2 pi
    zgrbth = zgrbth/nodes
    zgrbph = zgrbph/nodes


    ! calculate derivatives (wrt i)
    zgrbthp(1,:) = 0.5*(zgrbth(2,:) - zgrbth(nodes+1,:))/dtheta
    zgrbphp(1,:) = 0.5*(zgrbph(2,:) - zgrbph(nodes+1,:))/dtheta
    do i=2,nodes
       zgrbthp(i,:) = 0.5*(zgrbth(i+1,:) - zgrbth(i-1,:))/dtheta
       zgrbphp(i,:) = 0.5*(zgrbph(i+1,:) - zgrbph(i-1,:))/dtheta
    end do
    
100 continue

    close(ifile)
!
    else      ! on itaylor.eq.10
      allocate(zgrbph (nodes+1,nodes+1),zgrbth (nodes+1,nodes+1))
      allocate(zgrbphp(nodes+1,nodes+1),zgrbthp(nodes+1,nodes+1))
!...  define matrices with analytic formula
      xnode(nodes+1) = xnode(1)
      znode(nodes+1) = znode(1)
      zgrbph = (0.,0.)
      ka = aminor*ntor/xzero
      mmax = 30
      do m=1,mmax
        fac = 2.*m/(nodes*ka)*bessk(m,ka)/besskp(m,ka)
        do i=1,nodes+1
        do j=1,nodes+1
          thetai = atan2(znode(i),xnode(i)-xzero)
          thetaj = atan2(znode(j),xnode(j)-xzero)
          zgrbth(i,j) = zgrbth(i,j) + fac*sin(m*(thetai-thetaj))
        enddo
        enddo
      enddo

      ka = aminor*ntor/xzero
      grate = (eta_wall/(delta_wall*aminor))*(mpol + (mpol**2*bessk(mpol,ka))/(ka*besskp(mpol,ka)))
      write(*,1001) grate
      write(98,1001) grate, bessk(mpol,ka), besskp(mpol,ka), xzero
 1001 format(" Analytic Decay rate ",1p4e12.4)
    endif      ! on itaylor.eq.10

!      do i=1,nodes+1
!
!.....debug:   change write to read
!        read(97,1002) (zgrbth(i,j),j=1,nodes+1)
!        do j=1,nodes+1
!          zgrbthp(i,j) = 0.
!          zgrbph(i,j) = 0.
!          zgrbphp(i,j) = 0.
!        enddo
!      enddo
 1002 format(1p6e12.4)
    
  end subroutine load_response_matrix

  subroutine load_vacuum_data(ierr)
    implicit none

    integer, intent(out) :: ierr

    call load_boundary_nodes(ierr)
    if(ierr .ne. 0) return
    print *, 'boundary nodes loaded'
    
    call load_response_matrix(ierr)
    if(ierr .ne. 0) return
    print *, 'response matrix loaded'
      
  end subroutine load_vacuum_data

  subroutine unload_vacuum_data
    implicit none

    print *, 'unloading vacuum data'
    if(allocated(global_id)) deallocate(global_id)
    if(allocated(local_id)) deallocate(local_id)
    if(allocated(zgrbth)) deallocate(zgrbth)
    if(allocated(zgrbph)) deallocate(zgrbph)
    if(allocated(zgrbthp)) deallocate(zgrbthp)
    if(allocated(zgrbphp)) deallocate(zgrbphp)

  end subroutine unload_vacuum_data


end module vacuum_interface
