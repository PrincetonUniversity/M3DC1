module eqdsk

  character*10 :: name(6)      

  integer :: nw, nh, nbbbs, limitr
  real :: rdim, zdim, rcentr, rleft, zmid
  real :: rmaxis, zmaxis, simag, sibry, bcentr
  real :: current

  real, allocatable :: psirz(:,:),fpol(:),press(:),ffprim(:), &
       pprime(:),qpsi(:),rbbbs(:),zbbbs(:), &
       rlim(:),zlim(:)

  integer, private :: size, rank

contains

subroutine load_eqdsk

  implicit none

  include 'mpif.h'

  integer, parameter :: neqdsk = 20
  integer :: idum, i, j, ierr
  real :: xdum, ffp2, pp2, dpsi
  integer :: ibuff(2)
  real :: rbuff(20)

  call MPI_comm_size(MPI_COMM_WORLD,size,ierr)
  call MPI_comm_rank(MPI_COMM_WORLD,rank,ierr)

  if(rank.eq.0) then
     print *, 'Reading EQDSK g-file'
     open(unit=neqdsk,file='geqdsk',status='old')

2000 format (6a8,3i4)
2020 format (5e16.9)
2022 format (2i5)

     read (neqdsk,2000) (name(i),i=1,6),idum,nw,nh

     print *, ' Resolution: ', nw, nh
     ibuff(1) = nw
     ibuff(2) = nh
  endif

  call mpi_bcast(ibuff, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  nw = ibuff(1)
  nh = ibuff(2)

  allocate(psirz(nw,nh),fpol(nw),press(nw),ffprim(nw),pprime(nw),qpsi(nw))

  if(rank.eq.0) then
     read (neqdsk,2020) rdim,zdim,rcentr,rleft,zmid
     read (neqdsk,2020) rmaxis,zmaxis,simag,sibry,bcentr
     read (neqdsk,2020) current,xdum,xdum,xdum,xdum
     read (neqdsk,2020) xdum,xdum,xdum,xdum,xdum
     print *, 'reading fpol'
     read (neqdsk,2020) (fpol(i),i=1,nw)
     print *, 'reading press'
     read (neqdsk,2020) (press(i),i=1,nw)
     print *, 'reading ffprim'
     read (neqdsk,2020) (ffprim(i),i=1,nw)
     print *, 'reading pprim'
     read (neqdsk,2020) (pprime(i),i=1,nw)
     print *, 'reading psirz'
     read (neqdsk,2020) ((psirz(i,j),i=1,nw),j=1,nh)
     print *, 'reading qpsi'
     read (neqdsk,2020) (qpsi(i),i=1,nw)
     print *, 'reading nbbbs, limitr'
!!$  read (neqdsk,2022) nbbbs,limitr
!!$
!!$  print *, 'nbbbs, limitr', nbbbs, limitr
!!$  allocate(rbbbs(nbbbs),zbbbs(nbbbs),rlim(limitr),zlim(limitr))
!!$  read (neqdsk,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)

     print *, name
     write(*,'(A,2e12.4)') " Magnetic axis: ", rmaxis, zmaxis
     write(*,'(A,1e12.4)') " Toroidal field on axis: ", bcentr
     write(*,'(A,1e12.4)') " Plasma current: ", current
     write(*,'(A,2e12.4)') " Poloidal flux at axis, boundary", simag, sibry
     write(*,'(A,4e12.4)') " Bounding box:", rleft, zmid-zdim/2., rdim, zdim

!  call read_out1(neqdsk)

     close(neqdsk)

     if(current.lt.0) then
        print *, 'Current is negative; flipping flux'
        simag = -simag
        sibry = -sibry
        psirz = -psirz
        pprime = -pprime
        ffprim = -ffprim
     end if
  endif

  rbuff(1)  = rdim
  rbuff(2)  = zdim
  rbuff(3)  = rcentr
  rbuff(4)  = rleft
  rbuff(5)  = zmid
  rbuff(6)  = rmaxis
  rbuff(7)  = zmaxis
  rbuff(8)  = simag
  rbuff(9)  = sibry
  rbuff(10) = bcentr
  rbuff(11) = current
  call mpi_bcast(rbuff, 11, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  rdim    = rbuff(1)
  zdim    = rbuff(2)
  rcentr  = rbuff(3)
  rleft   = rbuff(4)
  zmid    = rbuff(5)
  rmaxis  = rbuff(6)
  zmaxis  = rbuff(7)
  simag   = rbuff(8)
  sibry   = rbuff(9)
  bcentr  = rbuff(10)
  current = rbuff(11)
  call mpi_bcast(fpol,   nw, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(press,  nw, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(ffprim, nw, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(pprime, nw, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(psirz,  nw*nh, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(qpsi,   nw, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

         if(rank.eq.0) then
           open(unit=78,file="eqdsk-out",status="unknown")
           write(78,2010) sibry,simag,current
2010       format("sibry,simag,current =",1p3e12.4,/,  &
                "  l   press       pprime      fpol        ffprim      ffp2        pprime2     qpsi")
                                                                                                        
           dpsi = (sibry-simag)/(nw-1.)
           do i=1,nw
              if(i.gt.1 .and. i.lt.nw)  then
                ffp2 = fpol(i)*(fpol(i+1)-fpol(i-1))/(2*dpsi)
                pp2 = (press(i+1) - press(i-1))/(2*dpsi)
              else
                ffp2 = 0
                pp2 = 0
              endif
              write(78,2011) i,press(i),pprime(i),fpol(i),ffprim(i),ffp2,pp2,qpsi(i)
           enddo
2011       format(i3,1p7e12.4)
           close(78)
        endif

  if(rank.eq.0) then
     print *, 'Done reading EQDSK g-file.'
  endif

end subroutine load_eqdsk

subroutine unload_eqdsk
  implicit none
  
  deallocate(psirz,fpol,press,ffprim,pprime,qpsi)
!!$  deallocate(rbbbs,zbbbs,rlim,zlim)
  
end subroutine unload_eqdsk


subroutine read_out1(neqdsk)

  implicit none

  integer, intent(in) :: neqdsk

  integer :: ishot, itime, icurrt, nbdry, mxiter, nxiter, limitr
  integer :: iconvr, ibunmn, nqpsi, npress
  real :: betap0, rzero, qenp, enp, emp, plasma, btor, fwtcur
  real :: error, rcentr
  real :: expmp2(76), coils(44),brsp(487),rbdry(110),zbdry(110), fwtsi(44)
  real :: xlim(160), ylim(160), pressr(132), rpress(132), sigpre(132)
  real :: pres(65), qpsi(65), pressw(65)
  namelist /out1/ &
       ishot, itime, betap0, rzero, qenp, enp, emp, plasma, expmp2, &
       coils, btor, rcentr, brsp, icurrt, rbdry, zbdry, nbdry, &
       fwtsi, fwtcur, mxiter, nxiter, limitr, xlim, ylim, error, &
       iconvr, ibunmn, pressr, rpress, qpsi, pressw, pres, nqpsi, npress, &
       sigpre

  print *, "Reading namelist..."

  read(neqdsk,nml=out1)

  print *, "Shot = ", ishot

end subroutine read_out1

end module
