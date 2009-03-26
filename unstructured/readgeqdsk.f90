module eqdsk

  character*10 :: name(6)      

  integer :: nw, nh, nbbbs, limitr
  real :: rdim, zdim, rcentr, rleft, zmid
  real :: rmaxis, zmaxis, simag, sibry, bcentr
  real :: current

  real, allocatable :: psirz(:,:),fpol(:),press(:),ffprim(:), &
       pprime(:),qpsi(:),rbbbs(:),zbbbs(:), &
       rlim(:),zlim(:)

contains

subroutine load_eqdsk

  implicit none

  integer, parameter :: neqdsk = 20
  integer :: idum, i, j
  real :: xdum

  print *, 'reading geqdsk'

  open(unit=neqdsk,file='geqdsk',status='old')

2000 format (6a8,3i4)
2020 format (5e16.9)
2022 format (2i5)

  read (neqdsk,2000) (name(i),i=1,6),idum,nw,nh
  allocate(psirz(nw,nh),fpol(nw),press(nw),ffprim(nw),pprime(nw),qpsi(nw))

  read (neqdsk,2020) rdim,zdim,rcentr,rleft,zmid
  read (neqdsk,2020) rmaxis,zmaxis,simag,sibry,bcentr
  read (neqdsk,2020) current,simag,xdum,rmaxis,xdum
  read (neqdsk,2020) zmaxis,xdum,sibry,xdum,xdum
  read (neqdsk,2020) (fpol(i),i=1,nw)
  read (neqdsk,2020) (press(i),i=1,nw)
  read (neqdsk,2020) (ffprim(i),i=1,nw)
  read (neqdsk,2020) (pprime(i),i=1,nw)
  read (neqdsk,2020) ((psirz(i,j),i=1,nw),j=1,nh)
  read (neqdsk,2020) (qpsi(i),i=1,nw)
  read (neqdsk,2022) nbbbs,limitr

  allocate(rbbbs(nbbbs),zbbbs(nbbbs),rlim(limitr),zlim(limitr))
  read (neqdsk,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)


  print *, name
  write(*,'(A,2e12.4)') "Magnetic axis: ", rmaxis, zmaxis
  write(*,'(A,1e12.4)') "Toroidal field on axis: ", bcentr
  write(*,'(A,1e12.4)') "Plasma current: ", current
  write(*,'(A,2e12.4)') "Poloidal flux at axis, boundary", simag, sibry

!  call read_out1(neqdsk)
  
  close(neqdsk)  

end subroutine load_eqdsk

subroutine unload_eqdsk
  implicit none
  
  deallocate(psirz,fpol,press,ffprim,pprime,qpsi, &
       rbbbs,zbbbs,rlim,zlim)
  
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
