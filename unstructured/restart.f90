subroutine wrrestart
  use mesh_mod
  use basic
  use arrays
  use diagnostics
  use time_step
  use pellet

  implicit none

#ifdef USESCOREC

  integer :: mmnn18, j1, numnodes, numelms
  integer :: ndofs
  integer, save :: ifirstrs = 1
  character (len=30) :: fname, oldfname

  call createfilename(fname, oldfname)
  numnodes = local_nodes()
  numelms = local_elements()
  mmnn18 = 0
  call numdofs(num_fields, ndofs)
  
  if(ifirstrs .ne. 1) call rename(fname, oldfname)
  ifirstrs = 0

  open(56,file=fname,form='unformatted',status='replace', action='write')
  ! first put in information to check on run information
  write(56) numnodes
  write(56) numelms
  write(56) mmnn18
  write(56) numvar
  write(56) iper
  write(56) jper 
  write(56) myrank
  write(56) maxrank
  write(56) eqsubtract
  write(56) linear
  write(56) icomplex
  
  do j1=1,ndofs 
     write(56) field_vec%data(j1)
  enddo
  do j1=1,ndofs 
     write(56) field0_vec%data(j1)
  enddo

  write(56) ntime,time,dt
  write(56) totcur0,tflux0,gbound,ptot,vloop,   &
          i_control%err_i, i_control%err_p_old, n_control%err_i, n_control%err_p_old
  write(56) psimin,psilim,psibound
  write(56) xnull, znull
  write(56) xmag, zmag
  call numdofs(1, ndofs)
  do j1=1,ndofs 
     write(56) bf_field(1)%vec%data(j1)
  enddo
  do j1=1,ndofs 
     write(56) bf_field(0)%vec%data(j1)
  enddo
  write(56) pellet_x, pellet_phi, pellet_z, &
       pellet_velx, pellet_velphi, pellet_velz, pellet_var
  write(56) version

  close(56)

#endif

end subroutine wrrestart

!============================================================
subroutine rdrestart
  use mesh_mod
  use basic
  use arrays
  use diagnostics
  use time_step
  use pellet

  implicit none
  
#ifdef USESCOREC

  integer :: j1, numnodes, inumnodes
  integer :: inumelms, immnn18, inumvar, iiper, ijper, imyrank
  integer :: imaxrank, numelms, ieqsubtract, ilinear, icomp
  character (len=30) :: fname, oldfname
  integer :: ndofs
  integer :: iversion

  call createfilename(fname, oldfname)
  call numdofs(num_fields, ndofs)
  numnodes = local_nodes()
  numelms = local_elements()

  open(56,file=fname,form='unformatted',status='unknown')
  read(56) inumnodes
  read(56) inumelms
  read(56) immnn18
  read(56) inumvar
  read(56) iiper
  read(56) ijper 
  read(56) imyrank
  read(56) imaxrank
  read(56) ieqsubtract
  read(56) ilinear
  read(56) icomp

  if(inumnodes .ne. numnodes .or. inumelms .ne. numelms .or. &
       iiper .ne. iper .or. ijper .ne. jper .or. &
       imyrank .ne. myrank .or. imaxrank .ne. maxrank) then
     write(*,*) 'Restart file information does not match!'
     close(56)
     if(inumnodes .ne. numnodes) then
        write(*,*) 'numnodes ',inumnodes, numnodes, myrank 
     endif
     if(inumelms .ne. numelms) then
        write(*,*) 'numelms ',inumnodes, numnodes, myrank 
     endif
     if(iiper .ne. iper) then
        write(*,*) 'iper',iiper, iper, myrank
     endif
     if(ijper .ne. jper) then
        write(*,*) 'jper',ijper, jper, myrank
     endif
     if(imyrank .ne. myrank) then
        write(*,*) 'myrank',imyrank,myrank
     endif
     if(imaxrank .ne. maxrank) then
        write(*,*) 'maxrank',imaxrank, maxrank, myrank
     endif
     call safestop(2)
  endif

  do j1=1,ndofs 
     read(56) field_vec%data(j1)
  enddo
  do j1=1,ndofs 
     read(56) field0_vec%data(j1)
  enddo
                         
  ! If we are running a linear simulation, but the restart file was
  ! a nonlinear simulation, make the restart data be the equilibrium
  if(linear.eq.1 .and. ilinear.eq.0) then
     field0_vec%data = field0_vec%data + field_vec%data
     field_vec%data = 0.
  endif

  read(56) ntime,time,dt
  read(56) totcur0,tflux0,gbound,ptot,vloop,   &
          i_control%err_i, i_control%err_p_old, n_control%err_i, n_control%err_p_old
  read(56,END=1199) psimin,psilim,psibound
  read(56,END=1199) xnull,znull
  read(56,END=1199) xmag,zmag

  call numdofs(1, ndofs)
  do j1=1,ndofs 
     read(56,END=1199) bf_field(1)%vec%data(j1)
  enddo
  do j1=1,ndofs 
     read(56,END=1199) bf_field(0)%vec%data(j1)
  enddo

  read(56, END=1199) pellet_x, pellet_phi, pellet_z, &
       pellet_velx, pellet_velphi, pellet_velz, pellet_var

  read(56, END=1199) iversion

  goto 1200
1199 if(myrank.eq.0) &
          print *, 'Warning: reading from a previous restart version'
1200 close(56)

#endif
end subroutine rdrestart

!============================================================
subroutine createfilename(filename, oldfilename)
  implicit none

  include 'mpif.h'
  character (len=30) :: filename, oldfilename
  character (len=5) :: charprocnum
  integer :: myrank, j, ier, i

  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ier)
  if (ier /= 0) then
     print *,'Error in MPI_Comm_rank:',ier
     call safestop(1)
  endif
                                ! initialize the SUPERLU process grid
  do j=1,5
     i = (myrank / 10**(5-j)) - (myrank / 10**(6-j)) * 10
     select case(i)
     case(0)
        charprocnum(j:j) = '0'
     case(1)
        charprocnum(j:j) = '1'
     case(2)
        charprocnum(j:j) = '2'
     case(3)
        charprocnum(j:j) = '3'
     case(4)
        charprocnum(j:j) = '4'
     case(5)
        charprocnum(j:j) = '5'
     case(6)
        charprocnum(j:j) = '6'
     case(7)
        charprocnum(j:j) = '7'
     case(8)
        charprocnum(j:j) = '8'
     case(9)
        charprocnum(j:j) = '9'
     end select
  enddo

  filename = 'C1restart'//charprocnum
  oldfilename = 'C1restarto'//charprocnum
  
  return
end subroutine createfilename
!============================================================
subroutine wrrestartglobal
  use mesh_mod
  use basic
  use arrays
  use diagnostics
  implicit none

#ifdef USESCOREC
  include 'mpif.h'
  integer :: ifirstrs2, numdofsglobal, mpitag
  integer :: ierr, numentsglobal(4), i, numnodes, j, globalid
  integer :: begindofnumber, enddofnumberplusone
  integer :: recvinfo(1)
  integer :: status(MPI_STATUS_SIZE)
  data ifirstrs2/1/

  if(ifirstrs2 .ne. 1 .and. myrank .eq. 0) then
     call rename('C1restart', 'C1restarto')
  endif

  ifirstrs2 = 0
  mpitag = 9070
  recvinfo(1) = 1
  ierr = 0

  if(myrank .ne. 0) then
     ! unaligned access warnings on PPPL sgi machines 
     ! is from communication here
     call MPI_RECV(recvinfo, 1, MPI_INTEGER, myrank-1,  &
          mpitag, MPI_COMM_WORLD, status, ierr)
     if(ierr .ne. MPI_SUCCESS) then
        write(*,*) 'MPI recv problem in restart', ierr
        call safestop(828)
     endif

#ifdef _AIX
     print *, "must use wrrestart on ibm instead of wrrestartglobal"
     call safestop(9229)
#else
     open(56,file='C1restart',form='unformatted',status='old',&
          position='append')
#endif
     endfile 56
  else
     if(ifirstrs2 .ne. 1) call rename('C1restart', 'C1restarto')
     ifirstrs2 = 0
     open(56,file='C1restart',form='unformatted',status='replace', &
          action='write')
     ! first put in information to check on run information
     call numglobaldofs(num_fields, numdofsglobal)
     call numglobalents(numentsglobal(1),numentsglobal(2), &
          numentsglobal(3), numentsglobal(4))
     write(56) numentsglobal(1)
     write(56) numentsglobal(2)
     write(56) numentsglobal(3)
     write(56) numdofsglobal
     write(56) numvar
     write(56) iper
     write(56) jper   
     write(56) ntime,time,dt
     write(56) totcur0,tflux0,gbound,ptot,vloop,   &
          i_control%err_i, i_control%err_p_old, n_control%err_i, n_control%err_p_old
     write(56) icomplex
  endif

  numnodes = local_nodes()
  do j=1,numnodes
     call entglobalid(j, 0, globalid)
     write(56) globalid
     call entdofs(num_fields, j, 0, begindofnumber, &
          enddofnumberplusone)
     do i=begindofnumber,enddofnumberplusone-1
        write(56) field_vec%data(i)
     enddo
     do i=begindofnumber,enddofnumberplusone-1
        write(56) field0_vec%data(i)
     enddo
  enddo
  i = -1
  !  below specifies that this is the last line
  if(myrank .eq. maxrank-1) write(56) i

  if(myrank .ne. maxrank-1) then
     call MPI_SEND(recvinfo, 1, MPI_INTEGER, myrank+1, &
          mpitag, MPI_COMM_WORLD, ierr)
     if(ierr .ne. MPI_SUCCESS) then
        write(*,*) 'MPI send problem in restart'
        call safestop(827)
     endif
  endif

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  if(myrank.eq.0) then
     write(56) psimin,psilim,psibound,xnull,znull
     write(56) xmag,zmag
  endif
  close(56)
      
1201 format("* * * restart file written at cycle,time = ",i5,1pe12.4,    &
          "  * * *")
#endif
end subroutine wrrestartglobal
!============================================================
subroutine rdrestartglobal
  use mesh_mod
  use basic
  use arrays
  use diagnostics
  
  implicit none

#ifdef USESCOREC

  include 'mpif.h'

  integer :: numdofsglobal, numentsglobal(4), i, numnodes, globalid
  integer :: begindofnumber, enddofnumberplusone, in1, in2, in3, nodeid
  integer :: icomp
  vectype :: tempvariable
  real :: realtemp

  open(56,file='C1restart',form='unformatted',status='old')
      
  ! first put in information to check on run information
  call numglobaldofs(num_fields, numdofsglobal)
  call numglobalents(numentsglobal(1),numentsglobal(2), &
       numentsglobal(3), numentsglobal(4))
  read(56) in1!, in2, in3
  read(56) in2
  read(56) in3
  write(*,*) 'in',in1, in2, in3
  if(in1 .ne. numentsglobal(1) .or. in2 .ne. numentsglobal(2) .or. &
       in3 .ne. numentsglobal(3)) then
     write(*,*) 'Mesh does not match for rdrestartglobal', &
          in1, in2, in3, &
          numentsglobal(1), numentsglobal(2), numentsglobal(3)
     call safestop(567)
  endif
  read(56) numdofsglobal
  read(56) numvar
  read(56) iper
  read(56) jper   
  read(56) ntime,time,dt
  read(56) totcur0,tflux0,gbound,ptot,vloop,   &
          i_control%err_i, i_control%err_p_old, n_control%err_i, n_control%err_p_old
  read(56) icomp
  write(*,*) 'in restartglobal time, dts are ', ntime, time, dt
  numnodes = local_nodes()
  read(56) globalid
  do while(globalid .ne. -1)
     call globalidnod(globalid, nodeid)
     if(nodeid .ne. -1) then ! this node exists on this proc
        
        if(icomp.eq.0 .and. icomplex.eq.1) then
           call entdofs(num_fields, nodeid, 0, begindofnumber, &
                enddofnumberplusone)
           do i=begindofnumber,enddofnumberplusone-1
              read(56) realtemp
              field_vec%data(i) = realtemp
           enddo
           do i=begindofnumber,enddofnumberplusone-1
              read(56) realtemp
              field0_vec%data(i) = realtemp
           enddo
        else
           call entdofs(num_fields, nodeid, 0, begindofnumber, &
                enddofnumberplusone)
           do i=begindofnumber,enddofnumberplusone-1
              read(56) field_vec%data(i)
           enddo
           do i=begindofnumber,enddofnumberplusone-1
              read(56) field0_vec%data(i)
           enddo
        endif
     else ! node is not on this proc
        if(icomp.eq.0 .and. icomplex.eq.1) then
           ! the number of dof values at this node
           do i=1,num_fields*6*3 
              read(56) realtemp
           enddo
        else
           ! the number of dof values at this node
           do i=1,num_fields*6*3 
              read(56) tempvariable
           enddo
        endif
     endif
     read(56) globalid
  enddo
  
  read(56,END=1199) psimin,psilim,psibound
  write(*,*) 'psimin,psilim,psibound',psimin,psilim,psibound
  read(56,END=1199) xnull,znull
  read(56,END=1199) xmag,zmag
  goto 1200
1199 if(myrank.eq.0) &
          print *, 'Warning: reading from a previous restart version'
1200 close(56)

#endif
end subroutine rdrestartglobal



subroutine wrrestart_adios
  use mesh_mod
  use basic
  use arrays
  use diagnostics
  use time_step
  use pellet

  implicit none

   include 'mpif.h'

!cj aug05-2011 #ifdef USESCOREC 
#ifdef USEADIOS
  integer :: mmnn18, j1, numnodes, numelms
  integer :: ndofs_1, ndofs_2, i, j
  integer, save :: ifirstrs = 1
  character (len=30) :: fname, oldfname

  ! ADIOS variables declarations for matching gwrite_restart_c1.fh 
  integer             :: comm, ierr
  integer                 :: adios_err
  integer*8               :: adios_groupsize, adios_totalsize
  integer*8               :: adios_handle 
  real, allocatable :: tmp_field_vec(:), tmp_field0_vec(:)
  real, allocatable :: tmp_bf_field_1(:), tmp_bf_field_0(:)
  real, allocatable :: tmp_psi_ext(:), tmp_bz_ext(:), tmp_bf_ext(:)

  integer :: useext

  if(ntime.lt.2) return
  fname="restart.bp"
  oldfname="restarto.bp"
  numnodes = local_nodes()
  numelms = local_elements()

  if(ifirstrs .ne. 1) call rename(fname, oldfname)

  call numdofs(num_fields, ndofs_1)
  allocate(tmp_field_vec(ndofs_1))
  allocate(tmp_field0_vec(ndofs_1))
  tmp_field_vec(1:ndofs_1) = field_vec%data(1:ndofs_1)
  tmp_field0_vec(1:ndofs_1) = field0_vec%data(1:ndofs_1)

  call numdofs(1, ndofs_2)
  allocate(tmp_bf_field_1(ndofs_2))
  allocate(tmp_bf_field_0(ndofs_2))
  tmp_bf_field_1(1:ndofs_2)= bf_field(1)%vec%data(1:ndofs_2)
  tmp_bf_field_0(1:ndofs_2) = bf_field(0)%vec%data(1:ndofs_2)

  allocate(tmp_psi_ext(ndofs_2))
  allocate(tmp_bz_ext(ndofs_2))
  allocate(tmp_bf_ext(ndofs_2))
  if(use_external_fields) then
     tmp_psi_ext(1:ndofs_2)= psi_ext%vec%data(1:ndofs_2)
     tmp_bz_ext(1:ndofs_2) = bz_ext%vec%data(1:ndofs_2)
     tmp_bf_ext(1:ndofs_2) = bf_ext%vec%data(1:ndofs_2)
     useext = 1
  else
     useext = 0
  end if


    call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr) 
    if( ifirstrs .eq. 1 ) then
#ifdef USECOMPLEX
    call adios_init ("m3dc1_cplx.xml", adios_err)
#else
    call adios_init ("m3dc1.xml", adios_err)
#endif
    endif
    call adios_open (adios_handle, "restart", fname, "w", comm, adios_err)
#ifdef USECOMPLEX
#include "gwrite_restart_c1_cplx.fh" 
#else
#include "gwrite_restart_c1.fh" 
#endif
    call adios_close (adios_handle, adios_err)
    call MPI_Barrier (comm, ierr)
    if( (ntimemax-(ntime-ntime0)) .lt. ntimepr ) then
    call adios_finalize (myrank, adios_err)
    endif

  if(myrank.eq.0) &
      write(*,*) "OUTPUT: wrrestart_adios groupsize totalsize", &
                  adios_groupsize, adios_totalsize

  deallocate(tmp_field_vec, tmp_field0_vec, tmp_bf_field_1, tmp_bf_field_0, &
       tmp_psi_ext, tmp_bz_ext, tmp_bf_ext)

  ifirstrs = 0
#endif
end subroutine wrrestart_adios

!============================================================
subroutine rdrestart_adios
  use mesh_mod
  use basic
  use arrays
  use diagnostics
  use time_step
  use pellet
  implicit none
  
   include 'mpif.h'

!cj aug05-2011 #ifdef USESCOREC
#ifdef USEADIOS

  integer :: j1, numnodes, inumnodes
  integer :: inumelms, immnn18, inumvar, iiper, ijper, imyrank
  integer :: imaxrank, numelms, ieqsubtract, ilinear, icomp
  character (len=30) :: fname, oldfname
  integer :: ndofs1, ndofs2, ndofs_1, ndofs_2

  ! ADIOS variables declarations for matching gread_restart_c11.fh &  gread_restart_c12.fh
  integer             :: comm, ierr, ier, itmp
  integer                 :: adios_err
  integer*8               :: adios_groupsize, adios_totalsize
  integer*8               :: adios_handle, adios_buf_size
  real, allocatable :: tmp_field_vec(:), tmp_field0_vec(:)
  real, allocatable :: tmp_bf_field_1(:), tmp_bf_field_0(:)
  real, allocatable :: tmp_psi_ext(:), tmp_bz_ext(:), tmp_bf_ext(:)

  !!! ADIOS variables for reading 
  integer*8                    :: fh, gh          ! file handler and group handler
  character(10), dimension(1)  :: gnamelist       ! expect 1 group only in restart file
  character(30), dimension(50) :: vnamelist, adios_varlist   ! list of all vars 
  character(30), dimension(10) :: anamelist, adios_attrlist  ! list of all attributes 
  integer                      :: tstart, ntsteps ! timesteps available
  integer*8                    :: start, vstart(2), readsize, vreadsize(2)
  integer                      :: i, gcnt, vcnt, acnt, ts, lastts
  integer*8                    :: read_bytes ! bytes read by adios_read_var()
  integer                      :: vartype, ndim, timedim ! adios_inq_var()
  integer*8, dimension(2)      :: dims                   ! adios_inq_var()
  integer                      :: elemsize               ! double complex or double
  integer :: iuseext
  integer :: iversion


  fname="restart.bp"
  oldfname="restarto.bp"
  numnodes = local_nodes()
  numelms = local_elements() 
  call numdofs(num_fields, ndofs1)
  call numdofs(1, ndofs2)

    call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
!#ifdef USECOMPLEX
!    call adios_init ("m3dc1_cplx.xml", adios_err)
!#else
!    call adios_init ("m3dc1.xml", adios_err)
!#endif
!    call adios_open (adios_handle, "restart", fname, "r", comm, adios_err)
!#include "gread_restart_c11.fh" 
!    call adios_close (adios_handle, adios_err)
!    ! Note, we have to close to perform the reading of the variables above

    !! set read method, default=0 the BP reader, 1=BP subfiles with better performance
    call adios_set_read_method (0, adios_err);
    call adios_fopen (fh,fname,MPI_COMM_WORLD,gcnt,adios_err)
    if (adios_err .ne. 0) call safestop(2)
    call adios_inq_file (fh,vcnt,acnt,tstart,ntsteps,gnamelist,adios_err)

!cj    write (*,'("Number of timesteps : ",i0," starting from ",i0)') ntsteps, tstart
!cj    write (*,'("Number of groups : ",i0)') gcnt

!cj    do i=1,gcnt
!cj        write (*,"(i5, a, a)") i,")  ", trim(gnamelist(i))
!cj    enddo

    !! Have to open a group from the file to read anything
    call adios_gopen (fh, gh, gnamelist(1), vcnt, acnt, ierr)
    ! not necessary to inquire if you know what to read
    call adios_inq_group(gh, vnamelist, anamelist, ts, lastts, ierr)

!cj    write (*,'("Number of variables in group ",a,": ",i0)') trim(gnamelist(1)), vcnt
!cj    do i=1,vcnt
!cj        write (*,"(i5, a, a)") i,")  ", trim(vnamelist(i))
!cj    enddo

!cj    write (*,'("Number of attributes in group ",a,": ",i0)') trim(gnamelist(1)), acnt
!cj    do i=1,acnt
!cj        write (*,"(i5, a, a)") i,")  ", trim(anamelist(i))
!cj    enddo

    ! Read in scalars (offset = 0, read size = 1 "in all dimensions")
    start = 0
    readsize = 1 ! number of elements of a type, not bytes!

    ! These scalars are varying over processors so we use adios_read_local_var() to 
    !  read the myrank'th value of the scalar from the file
    call adios_read_local_var (gh, "numnodes",   myrank, start, readsize, inumnodes, read_bytes)
    call adios_read_local_var (gh, "numelms",    myrank, start, readsize, inumelms, read_bytes)
    call adios_read_local_var (gh, "mmnn18",     myrank, start, readsize, immnn18, read_bytes)
    call adios_read_local_var (gh, "numvar",     myrank, start, readsize, inumvar, read_bytes)
    call adios_read_local_var (gh, "iper",       myrank, start, readsize, iiper, read_bytes)
    call adios_read_local_var (gh, "jper",       myrank, start, readsize, ijper, read_bytes)
    call adios_read_local_var (gh, "myrank",     myrank, start, readsize, imyrank, read_bytes)
    call adios_read_local_var (gh, "maxrank",    myrank, start, readsize, imaxrank, read_bytes)
    call adios_read_local_var (gh, "eqsubtract", myrank, start, readsize, ieqsubtract, read_bytes)
    call adios_read_local_var (gh, "useext",     myrank, start, readsize, iuseext, read_bytes)
    call adios_read_local_var (gh, "linear",     myrank, start, readsize, ilinear, read_bytes)
    call adios_read_local_var (gh, "icomplex",   myrank, start, readsize, icomp, read_bytes)
    call adios_read_local_var (gh, "ntime",      myrank, start, readsize, ntime, read_bytes)
    call adios_read_local_var (gh, "time",       myrank, start, readsize, time, read_bytes)
    call adios_read_local_var (gh, "dt",       myrank, start, readsize, dt, read_bytes)
    call adios_read_local_var (gh, "totcur0",    myrank, start, readsize, totcur0, read_bytes)
    call adios_read_local_var (gh, "tflux0",     myrank, start, readsize, tflux0, read_bytes)
    call adios_read_local_var (gh, "gbound",     myrank, start, readsize, gbound, read_bytes)
    call adios_read_local_var (gh, "ptot",       myrank, start, readsize, ptot, read_bytes)
    call adios_read_local_var (gh, "vloop",      myrank, start, readsize, vloop, read_bytes)
    call adios_read_local_var (gh, "i_control%err_i",      myrank, start, readsize, i_control%err_i, read_bytes)
    call adios_read_local_var (gh, "i_control%err_p_old",      myrank, start, readsize, i_control%err_p_old, read_bytes)
    call adios_read_local_var (gh, "n_control%err_i",      myrank, start, readsize, n_control%err_i, read_bytes)
    call adios_read_local_var (gh, "n_control%err_p_old",      myrank, start, readsize, n_control%err_p_old, read_bytes)
    call adios_read_local_var (gh, "psimin",     myrank, start, readsize, psimin, read_bytes)
    call adios_read_local_var (gh, "psilim",     myrank, start, readsize, psilim, read_bytes)
    call adios_read_local_var (gh, "psibound",     myrank, start, readsize, psibound, read_bytes)
    call adios_read_local_var (gh, "xnull",      myrank, start, readsize, xnull, read_bytes)
    call adios_read_local_var (gh, "znull",      myrank, start, readsize, znull, read_bytes)
    call adios_read_local_var (gh, "xmag",       myrank, start, readsize, xmag, read_bytes)
    call adios_read_local_var (gh, "zmag",       myrank, start, readsize, zmag, read_bytes)
    call adios_read_local_var (gh, "ndofs_1",    myrank, start, readsize, ndofs_1, read_bytes)
    call adios_read_local_var (gh, "ndofs_2",    myrank, start, readsize, ndofs_2, read_bytes)
    call adios_read_local_var (gh, "pellet_x",    myrank, start, readsize, pellet_x, read_bytes)
    call adios_read_local_var (gh, "pellet_phi",    myrank, start, readsize, pellet_phi, read_bytes)
    call adios_read_local_var (gh, "pellet_z",    myrank, start, readsize, pellet_z, read_bytes)
    call adios_read_local_var (gh, "pellet_velx",    myrank, start, readsize, pellet_velx, read_bytes)
    call adios_read_local_var (gh, "pellet_velphi",    myrank, start, readsize, pellet_velphi, read_bytes)
    call adios_read_local_var (gh, "pellet_velz",    myrank, start, readsize, pellet_velz, read_bytes)
    call adios_read_local_var (gh, "pellet_var",    myrank, start, readsize, pellet_var, read_bytes)
    call adios_read_local_var (gh, "version",    myrank, start, readsize, iversion, read_bytes)

  if(inumnodes .ne. numnodes .or. inumelms .ne. numelms .or. &
     numvar .ne. inumvar .or. &
     iiper .ne. iper .or. ijper .ne. jper .or. &
     imyrank .ne. myrank .or. imaxrank .ne. maxrank .or. &
     ndofs1 .ne. ndofs_1 .or. ndofs2.ne. ndofs_2) then
     write(*,*) 'Restart file information does not match!'

     if(inumnodes.ne.numnodes) print *, 'numnodes ',inumnodes,numnodes,myrank 
     if(inumelms .ne. numelms) print *, 'numelms ',inumelms,numelms,myrank 
     if(iiper .ne. iper) write(*,*) 'iper',iiper, iper, myrank
     if(ijper .ne. jper) write(*,*) 'jper',ijper, jper, myrank
     if(imyrank .ne. myrank) write(*,*) 'myrank',imyrank,myrank
     if(ndofs1 .ne. ndofs_1) write(*,*) 'ndofs 1',myrank, ndofs1, ndofs_1
     if(ndofs2 .ne. ndofs_2) write(*,*) 'ndofs 2',myrank, ndofs2, ndofs_2
     ierr = 1
  else
     ierr = 0
  end if

  ! check for errors
  call mpi_allreduce(ierr,itmp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ier)
  if(itmp.ne.0) call safestop(6)

    ! Allocate space for the arrays tmp_
  allocate(tmp_field_vec(ndofs_1))
  allocate(tmp_field0_vec(ndofs_1)) 
  allocate(tmp_bf_field_1(ndofs_2))
  allocate(tmp_bf_field_0(ndofs_2)) 
  allocate(tmp_psi_ext(ndofs_2))
  allocate(tmp_bz_ext(ndofs_2))
  allocate(tmp_bf_ext(ndofs_2))

    ! Check types of the array variables
    call adios_inq_var (gh, "tmp_field_vec", vartype, ndim, dims, timedim, adios_err)
#ifdef USECOMPLEX 
    elemsize = 16 ! sizeof double complex
    if (vartype .ne. 11) then  ! 11 = double complex in adios
#else
    elemsize = 8 ! sizeof double
    if (vartype .ne. 6) then  ! 6 = double in adios
#endif
         write(*,*) 'Type mismatch (complex vs real) of variable tmp_field_vec in restart file!'
         call safestop(2)
    endif

    

    ! Read in arrays organized as a global array: offset = (0, rank), read size = (*,1)
    start    = 0
!   vstart(2)    = myrank
    readsize = ndofs_1 ! number of elements of a type, not bytes!
!   vreadsize(2) = 1
    
    call adios_read_local_var (gh, "tmp_field_vec", myrank, start, readsize, tmp_field_vec, read_bytes)
    if (read_bytes .ne. elemsize*ndofs_1) then
         write(*,*) 'Size mismatch at reading tmp_field_vec!', read_bytes, elemsize*ndofs_1
         call safestop(2)
    endif
    
    call adios_read_local_var (gh, "tmp_field0_vec", myrank, start, readsize, tmp_field0_vec, read_bytes)
    if (read_bytes .ne. elemsize*ndofs_1) then 
         write(*,*) 'Size mismatch at reading tmp_field0_vec!', read_bytes, elemsize*ndofs_1
         call safestop(2)
    endif
    
    readsize = ndofs_2 ! number of elements of a type, not bytes!
    
    call adios_read_local_var (gh, "tmp_bf_field_1", myrank, start, readsize, tmp_bf_field_1, read_bytes)
    if (read_bytes .ne. elemsize*ndofs_2) then  
         write(*,*) 'Size mismatch at reading tmp_bf_field_1!', read_bytes, elemsize*ndofs_2
         call safestop(2)
    endif
    
    call adios_read_local_var (gh, "tmp_bf_field_0", myrank, start, readsize, tmp_bf_field_0, read_bytes)
    if (read_bytes .ne. elemsize*ndofs_2) then 
         write(*,*) 'Size mismatch at reading tmp_bf_field_1!', read_bytes, elemsize*ndofs_2
         call safestop(2)
    endif

    call adios_read_local_var (gh, "tmp_psi_ext", myrank, start, readsize, tmp_psi_ext, read_bytes)
    if (read_bytes .ne. elemsize*ndofs_2) then 
         write(*,*) 'Size mismatch at reading tmp_psi_ext!', read_bytes, elemsize*ndofs_2
         call safestop(2)
    endif

    call adios_read_local_var (gh, "tmp_bz_ext", myrank, start, readsize, tmp_bz_ext, read_bytes)
    if (read_bytes .ne. elemsize*ndofs_2) then 
         write(*,*) 'Size mismatch at reading tmp_bz_ext!', read_bytes, elemsize*ndofs_2
         call safestop(2)
    endif

    call adios_read_local_var (gh, "tmp_bf_ext", myrank, start, readsize, tmp_bf_ext, read_bytes)
    if (read_bytes .ne. elemsize*ndofs_2) then 
         write(*,*) 'Size mismatch at reading tmp_bf_ext!', read_bytes, elemsize*ndofs_2
         call safestop(2)
    endif
    
    !! Close group and file
    call MPI_Barrier (comm, ierr)
    call adios_gclose(gh, ierr)
    call adios_fclose(fh, ierr)

!    call adios_open (adios_handle, "restart", fname, "r", comm, adios_err)

!#ifdef USECOMPLEX
!#include "gread_restart_c12_cplx.fh" 
!#else
!#include "gread_restart_c12.fh" 
!#endif
!    call adios_close (adios_handle, adios_err) 
!    call MPI_Barrier (comm, ierr) 
!    call adios_finalize (myrank, adios_err)

      field_vec%data(1:ndofs_1) = tmp_field_vec(1:ndofs_1)
      field0_vec%data(1:ndofs_1) = tmp_field0_vec(1:ndofs_1)

  ! If we are running a linear simulation, but the restart file was
  ! a nonlinear simulation, make the restart data be the equilibrium
  if(linear.eq.1 .and. ilinear.eq.0) then
     field0_vec%data = field0_vec%data + field_vec%data
     field_vec%data = 0.
  endif 

  bf_field(1)%vec%data(1:ndofs_2) = tmp_bf_field_1(1:ndofs_2)
  bf_field(0)%vec%data(1:ndofs_2) = tmp_bf_field_0(1:ndofs_2) 

  if(iuseext.eq.1) then 
     use_external_fields = .true.
     call create_field(psi_ext)
     call create_field(bz_ext)
     call create_field(bf_ext)
     psi_ext%vec%data(1:ndofs_2) = tmp_psi_ext(1:ndofs_2)
     bz_ext%vec%data(1:ndofs_2) = tmp_bz_ext(1:ndofs_2)
     bf_ext%vec%data(1:ndofs_2) = tmp_bf_ext(1:ndofs_2)
  end if


!  if(myrank.eq.0) &
!      write(*,*) "INPUT: rdrestart_adios buf_size groupsize totalsize", &
!                 adios_buf_size, adios_groupsize, adios_totalsize

  deallocate(tmp_field_vec, tmp_field0_vec, tmp_bf_field_1, tmp_bf_field_0, &
       tmp_psi_ext, tmp_bz_ext, tmp_bf_ext) 
#endif
end subroutine rdrestart_adios
