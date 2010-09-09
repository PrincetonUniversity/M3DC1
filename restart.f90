subroutine wrrestart
  use p_data
  use mesh_mod
  use basic
  use arrays
  use diagnostics
  use time_step

  implicit none

#ifdef USESCOREC

  integer :: mmnn18, j1, numnodes, numelms
  integer :: ndofs
  integer, save :: ifirstrs = 1
  character (len=30) :: fname, oldfname

  call createfilename(fname, oldfname)
  numnodes = local_nodes()
  numelms = local_elements()
  call numdofs(vecsize_phi, mmnn18)
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

  write(56) ntime,time
  write(56) totcur0,tflux0,gbound,ptot,vloop
  write(56) psimin, psilim
  write(56) xnull, znull
  write(56) xmag, zmag
  call numdofs(1, ndofs)
  do j1=1,ndofs 
     write(56) bf_field(1)%vec%data(j1)
  enddo
  do j1=1,ndofs 
     write(56) bf_field(0)%vec%data(j1)
  enddo


  close(56)

#endif

end subroutine wrrestart

!============================================================
subroutine rdrestart
  use p_data
  use mesh_mod
  use basic
  use arrays
  use diagnostics
  use time_step
  implicit none
  
#ifdef USESCOREC

  integer :: mmnn18, j1, numnodes, inumnodes
  integer :: inumelms, immnn18, inumvar, iiper, ijper, imyrank
  integer :: imaxrank, numelms, ieqsubtract, ilinear, icomp
  character (len=30) :: fname, oldfname
  integer :: ndofs

  call createfilename(fname, oldfname)
  call numdofs(vecsize_phi, mmnn18)
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

  read(56) ntime,time
  read(56) totcur0,tflux0,gbound,ptot,vloop
  read(56,END=1199) psimin,psilim
  read(56,END=1199) xnull,znull
  read(56,END=1199) xmag,zmag

  call numdofs(1, ndofs)
  do j1=1,ndofs 
     read(56,END=1199) bf_field(1)%vec%data(j1)
  enddo
  do j1=1,ndofs 
     read(56,END=1199) bf_field(0)%vec%data(j1)
  enddo

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
  use p_data
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
     open(56,file='C1restart',form='unformatted',status='old', ACCESS='append')
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
     write(56) ntime,time
     write(56) totcur0,tflux0,gbound,ptot,vloop
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
     write(56) psimin,psilim,xnull,znull
     write(56) xmag,zmag
  endif
  close(56)
      
1201 format("* * * restart file written at cycle,time = ",i5,1pe12.4,    &
          "  * * *")
#endif
end subroutine wrrestartglobal
!============================================================
subroutine rdrestartglobal
  use p_data
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
  read(56) ntime,time
  read(56) totcur0,tflux0,gbound,ptot,vloop
  read(56) icomp
  write(*,*) 'in restartglobal times are ', ntime, time
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
  
  read(56,END=1199) psimin,psilim
  write(*,*) 'psimin,psilim',psimin,psilim
  read(56,END=1199) xnull,znull
  read(56,END=1199) xmag,zmag
  goto 1200
1199 if(myrank.eq.0) &
          print *, 'Warning: reading from a previous restart version'
1200 close(56)

#endif
end subroutine rdrestartglobal

