Program m3dc1_skeleton
  use mesh
  use m3dc1_data
  use matdef
  use etimes
  use omp_lib

  implicit none

  include 'mpif.h'

  integer :: ierr
  real :: tstart, tend
  real :: tempin(2), tempout(2)
  integer :: ithread, isize
  intrinsic cpu_time

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, maxrank, ierr)

  !2019.april.11
  !the initial size is 64(mesh_size) on 8 cpu, 
  !now mesh_size=64*isize for weaking scaling profiling
  isize=1

  if(myrank.eq.0) print *, 'loading mesh...'
  call load_mesh(isize)

  if(myrank.eq.0) print *, 'allocating data...'
  call allocate_data

  if(myrank.eq.0) print *, 'calculating mesh data...'
  call tridef

  if(myrank.eq.0) print *, 'setting initial conditions...'
  call initialize_fields

!!$OMP PARALLEL
  !ithread = omp_get_thread_num()
  !if(myrank.eq.0 .and. ithread.eq.0) &
  !   print *, 'NUM_THREADS = ', omp_get_num_threads()
!!$OMP END PARALLEL

  if(myrank.eq.0) print *, 'calculating matrix elements...'
  call get_etime(tstart)
  call matdefall
  call get_etime(tend)

  print *, 'rank, elapsed time: ', myrank, tend - tstart

  if(maxrank.gt.1) then 
     tempin(1) = current
     tempin(2) = volume
     call MPI_Allreduce(tempin, tempout, 2, MPI_DOUBLE_PRECISION, &
          MPI_SUM, MPI_COMM_WORLD, ierr)
     current = tempout(1)
     volume  = tempout(2)
  end if
  if(myrank.eq.0) write(*, '(A,2G15.8)') 'current, volume = ', current, volume

  call safestop(0)

end Program m3dc1_skeleton

subroutine safestop(out)
  use mesh
  use m3dc1_data

  implicit none

  integer, intent(in) :: out

  integer :: ierr

  if(myrank.eq.0) print *, 'deallocating data...'
  call deallocate_data

  if(myrank.eq.0) print *, 'unloading mesh...'
  call unload_mesh
  call MPI_Finalize(ierr)

  if(myrank.eq.0) print *, 'Stopped with output ', out
  stop
end subroutine safestop

!============================================================
! rotation
! ~~~~~~~~
! calculates the rotation matrix rot given angle theta
!============================================================
subroutine rotation(rot,ndim,theta)
  implicit none

  integer, intent(in) :: ndim
  real, intent(in) :: theta
  real, intent(out) :: rot(ndim,*)

  integer :: i, j
  real :: r1(6,6), co, sn

  co = cos(theta)
  sn = sin(theta)
  do i=1,6
     do j=1,6
        r1(i,j) = 0.
     enddo
  enddo

  r1(1,1) = 1.

  r1(2,2) = co
  r1(2,3) = sn

  r1(3,2) = -sn
  r1(3,3) = co

  r1(4,4) = co**2
  r1(4,5) = 2.*sn*co
  r1(4,6) = sn**2

  r1(5,4) = -sn*co
  r1(5,5) = co**2-sn**2
  r1(5,6) = sn*co

  r1(6,4) = sn**2
  r1(6,5) = -2.*sn*co
  r1(6,6) = co**2

  do i=1,18
     do j=1,18
        rot(i,j) = 0.
     enddo
  enddo

  do i=1,6
     do j=1,6
        rot(i,j)       = r1(i,j)
        rot(i+6,j+6)   = r1(i,j)
        rot(i+12,j+12) = r1(i,j)
     enddo
  enddo

  return
end subroutine rotation


  !============================================================
  ! tridef
  ! ~~~~~~
  ! populates the *tri arrays
  !============================================================
  subroutine tridef
    use math
    use mesh

    implicit none

    include 'mpif.h'
  
    type(element_data) :: d
    integer :: itri, i, j, k, ii, jj, numelms, numnodes, ndofs
    real, dimension(coeffs_per_tri,coeffs_per_tri) :: ti 
    real, dimension(dofs_per_tri, dofs_per_tri) :: rot, newrot
    real :: sum, theta
    real :: norm(2), curv, x, z
    integer :: inode(nodes_per_element)
    logical :: is_boundary
    integer :: izone, izonedim

    numelms = local_elements()
    numnodes = local_nodes()
    ndofs = numnodes*dofs_per_node

    ! start the loop over triangles within a rectangular region
    do itri=1,numelms

       ! define a,b,c and theta
       call get_element_data(itri,d)
       
       ! define the Inverse Transformation Matrix that enforces the 
       ! condition that the normal slope between triangles has only 
       ! cubic variation
       call tmatrix(ti,coeffs_per_tri,d%a,d%b,d%c)
       
       ! calculate the rotation matrix rot
       theta = atan2(d%sn,d%co)
       call rotation(rot,dofs_per_tri,theta)
       
       newrot = 0.
       call get_element_nodes(itri, inode)

       do i=1, 3
          call boundary_node(inode(i), &
               is_boundary, izone, izonedim, norm, curv, x, z)

          k = (i-1)*6 + 1
          if(is_boundary) then
             newrot(k  ,k  ) = 1.
             newrot(k+1,k+1) =  norm(1)
             newrot(k+1,k+2) =  norm(2)
             newrot(k+1,k+3) =  curv*norm(2)**2
             newrot(k+1,k+4) = -curv*norm(1)*norm(2)
             newrot(k+1,k+5) =  curv*norm(1)**2
             newrot(k+2,k+1) = -norm(2)
             newrot(k+2,k+2) =  norm(1)
             newrot(k+2,k+3) =  2.*curv*norm(1)*norm(2)
             newrot(k+2,k+4) = -curv*(norm(1)**2 - norm(2)**2) 
             newrot(k+2,k+5) = -2.*curv*norm(1)*norm(2)
             newrot(k+3,k+3) =  norm(1)**2 
             newrot(k+3,k+4) =  norm(1)*norm(2)
             newrot(k+3,k+5) =  norm(2)**2
             newrot(k+4,k+3) = -2.*norm(1)*norm(2)
             newrot(k+4,k+4) =  norm(1)**2 - norm(2)**2
             newrot(k+4,k+5) =  2.*norm(1)*norm(2)
             newrot(k+5,k+3) =  norm(2)**2
             newrot(k+5,k+4) = -norm(1)*norm(2)
             newrot(k+5,k+5) =  norm(1)**2
          else
             do j=1, 6
                newrot(k+j-1,k+j-1) = 1.
             end do
          end if
       end do     

       ! form the matrix g by multiplying ti and rot
       do k=1, coeffs_per_tri
          do j=1, dofs_per_tri
             sum = 0.
             do ii = 1, dofs_per_tri
                sum = sum + ti(k,ii)*rot(ii,j)
             enddo
             gtri_old(k,j,itri) = sum

             sum = 0.
             do ii = 1, dofs_per_tri
                do jj=1, dofs_per_tri
                   sum = sum + newrot(j,jj)*ti(k,ii)*rot(ii,jj)
                end do
             enddo
             gtri(k,j,itri) = sum
          enddo
       enddo

       htri(1,1,itri) = 1.
#ifdef USE3D
       htri(2,1,itri) = 0.
       htri(3,1,itri) =-3./d%d**2
       htri(4,1,itri) = 2./d%d**3

       htri(1,2,itri) = 0.
       htri(2,2,itri) = 1.
       htri(3,2,itri) =-2./d%d
       htri(4,2,itri) = 1./d%d**2

       htri(1,3,itri) = 0.
       htri(2,3,itri) = 0.
       htri(3,3,itri) = 3./d%d**2
       htri(4,3,itri) =-2./d%d**3

       htri(1,4,itri) = 0.
       htri(2,4,itri) = 0.
       htri(3,4,itri) =-1./d%d
       htri(4,4,itri) = 1./d%d**2
#endif
    end do

  end subroutine tridef

