module vacuum_interface

  implicit none

  integer :: nodes
  integer, allocatable :: global_id(:)
  integer, allocatable :: local_id(:)
  real, allocatable :: normal(:,:)
  real, allocatable :: zgrbth(:,:), zgrbph(:,:)
  real, allocatable :: zgrbthp(:,:), zgrbphp(:,:)

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
    allocate(normal(2,nodes))

    ! write node entry for first node
    do i=1, nodes
       read(ifile, '(I8,4f12.8)') global_id(i), dummy, dummy, &
            normal(1,i), normal(2,i)
       call globalidnod(global_id(i),local_id(i))
    end do

    close(ifile)
    
  end subroutine load_boundary_nodes


  subroutine load_response_matrix(ierr)
    implicit none

    integer, intent(out) :: ierr
    
    character(len=*), parameter :: filename = 'RESPONSE-M3DC1'
    integer, parameter :: ifile = 2
    integer :: idum, i, j
    real :: dtheta

    ierr = 0

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
    zgrbth = zgrbth*dtheta
    zgrbph = zgrbph*dtheta

    ! calculate derivatives (wrt i)
    zgrbthp(1,:) = 0.5*(zgrbth(2,:) - zgrbth(nodes+1,:))/dtheta
    zgrbphp(1,:) = 0.5*(zgrbph(2,:) - zgrbph(nodes+1,:))/dtheta
    do i=2,nodes
       zgrbthp(i,:) = 0.5*(zgrbth(i+1,:) - zgrbth(i-1,:))/dtheta
       zgrbphp(i,:) = 0.5*(zgrbph(i+1,:) - zgrbph(i-1,:))/dtheta
    end do
    
100 continue

    close(ifile)

    
  end subroutine load_response_matrix

  subroutine load_vacuum_data(ierr)
    implicit none

    integer, intent(out) :: ierr

    call load_boundary_nodes(ierr)
    if(ierr .ne. 0) return
    
    call load_response_matrix(ierr)
    if(ierr .ne. 0) return
      
  end subroutine load_vacuum_data

  subroutine unload_vacuum_data
    implicit none

    print *, 'unloading vacuum data'
    if(allocated(global_id)) deallocate(global_id)
    if(allocated(local_id)) deallocate(local_id)
    if(allocated(normal)) deallocate(normal)
    if(allocated(zgrbth)) deallocate(zgrbth)
    if(allocated(zgrbph)) deallocate(zgrbph)
    if(allocated(zgrbthp)) deallocate(zgrbthp)
    if(allocated(zgrbphp)) deallocate(zgrbphp)

  end subroutine unload_vacuum_data


end module vacuum_interface
