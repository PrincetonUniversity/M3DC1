      subroutine wrrestart(timein, ntimein)
      use p_data
      use t_data
      use basic
      use arrays
      use diagnostics

      implicit none
      integer mmnn18, j1, mmnn6, ifirstrs, numnodes, numelms, ntimein
      integer :: ndofs
      real timein
      data ifirstrs/1/
      character (len=30) :: fname, oldfname
      integer :: comp

#ifdef USECOMPLEX
      comp = 1
#else
      comp = 0
#endif

      call createfilename(fname, oldfname)
      call numnod(numnodes)
      call numfac(numelms)
      call numdofs(vecsize, mmnn18)
      call numdofs(num_fields, ndofs)

      if(ifirstrs .ne. 1) call rename(fname, oldfname)
      ifirstrs = 0
      open(56,file=fname,form='unformatted',status='unknown')
c     first put in information to check on run information
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
      write(56) comp

      do j1=1,ndofs 
         write(56) field(j1)
      enddo
      do j1=1,ndofs 
         write(56) field0(j1)
      enddo
      do j1=1,ndofs 
         write(56) fieldi(j1)
      enddo

      do j1=1,mmnn18
         write(56) phiold(j1)
      enddo

      if(isplitstep.eq.1) then
         do j1=1,mmnn18
            write(56) velold(j1)
         enddo
         if(idens.eq.1) then
            call numdofs(1, mmnn6)
            do j1=1,mmnn6
               write(56) denold(j1)
            enddo
         endif
         if(ipres.eq.1) then 
            do j1=1,mmnn6       
               write(56) presold(j1)
            enddo
         endif
      endif
!
      write(56) ntimein,timein
      write(56) totcur0,tflux0,gbound,facd,fbound,ptot,vloop
      close(56)

      if(maxrank.eq.1) then
         write(9,1201) ntime,time
 1201   format("* * * restart file written at cycle,time = ",i5,1pe12.4,    
     &        "  * * *")
      endif

      return
      end
!============================================================
      subroutine rdrestart
      use p_data
      use t_data
      use basic
      use arrays
      use diagnostics
      implicit none
      integer mmnn18, j1, mmnn6, itime, igr, numnodes, inumnodes, 
     1     inumelms, immnn18, inumvar, iiper, ijper, imyrank,
     2     imaxrank, numelms, ieqsubtract, ilinear
      character (len=30) :: fname, oldfname
      real, allocatable :: temp(:)
      integer :: comp, icomp
      integer :: ndofs
#ifdef USECOMPLEX
      comp = 1
#else
      comp = 0
#endif

      call createfilename(fname, oldfname)
      call numdofs(numvar, mmnn18)
      call numdofs(num_fields, ndofs)
      call numnod(numnodes)
      call numfac(numelms)

      open(56,file=fname,form='unformatted',status='old')
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

      if(inumnodes .ne. numnodes .or. inumelms .ne. numelms .or.
     1     immnn18 .ne. mmnn18 .or. inumvar .ne. numvar .or.
     2     iiper .ne. iper .or. ijper .ne. jper .or.
     3     imyrank .ne. myrank .or. imaxrank .ne. maxrank) then
         write(*,*) 'Restart file information does not match!'
         close(56)
         if(inumnodes .ne. numnodes) then
            write(*,*) 'numnodes ',inumnodes, numnodes, myrank 
         endif
         if(inumelms .ne. numelms) then
            write(*,*) 'numelms ',inumnodes, numnodes, myrank 
         endif
         if(immnn18 .ne. mmnn18) then
            write(*,*) 'mmnn18',immnn18,mmnn18, myrank
         endif
         if(inumvar .ne. numvar) then
            write(*,*) 'numvar',inumvar,numvar, myrank
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

      ! If we are using complex vectors but the restart file was 
      ! written using real vectors, we must create a real buffer to 
      ! read input.
      if(icomp.eq.0 .and. comp.eq.1) then
         allocate(temp(ndofs))

         do j1=1,ndofs 
            read(56) temp
         enddo
         field = temp
         do j1=1,ndofs 
            read(56) temp(j1)
         enddo
         field0 = temp
         do j1=1,ndofs 
            read(56) temp
         enddo
         fieldi = temp(j1)

         do j1=1,mmnn18
            read(56) temp(j1)
         enddo
         phiold = temp(1:mmnn18)
         
         if(isplitstep.eq.1) then
            do j1=1,mmnn18
               read(56) temp(j1)
            enddo
            velold = temp(1:mmnn18)
            if(idens.eq.1) then
               call numdofs(1, mmnn6)
               do j1=1,mmnn6
                  read(56) temp(j1)
               enddo
               denold = temp(1:mmnn6)
            endif
            denold = temp(1:mmnn6)
            if(ipres.eq.1) then 
               do j1=1,mmnn6       
                  read(56) temp(j1)
               enddo
               presold = temp(1:mmnn6)
            endif
         endif

         deallocate(temp)

      ! Otherwise, just read in what was written.
      else 

         do j1=1,ndofs 
            read(56) field(j1)
         enddo
         do j1=1,ndofs 
            read(56) field0(j1)
         enddo
         do j1=1,ndofs 
            read(56) fieldi(j1)
         enddo
         
         do j1=1,mmnn18
            read(56) phiold(j1)
         enddo
         
         if(isplitstep.eq.1) then
            do j1=1,mmnn18
               read(56) velold(j1)
            enddo
            if(idens.eq.1) then
               call numdofs(1, mmnn6)
               do j1=1,mmnn6
                  read(56) denold(j1)
               enddo
            endif
            if(ipres.eq.1) then 
               do j1=1,mmnn6       
                  read(56) presold(j1)
               enddo
            endif
         endif


      endif
                         
      ! If we are running a linear simulation, but the restart file was
      ! a nonlinear simulation, make the restart data be the equilibrium
      if(linear.eq.1 .and. ilinear.eq.0) then
         field0 = field0 + field
         field = 0
      endif

      read(56) ntimer,timer
      read(56) totcur0,tflux0,gbound,facd,fbound,ptot,vloop
      close(56)
      write(*,1201) ntimer,timer
      if(maxrank.eq.1) then
         write(9,1201) ntimer,timer
 1201    format("* * * restart file read at cycle,time = ",i5,1pe12.4,     &
     &        "  * * *")
      endif
      return
      end
!============================================================
      subroutine createfilename(filename, oldfilename)
      implicit none
      include 'mpif.h'
      character (len=30) :: filename, oldfilename
      character (len=5) :: charprocnum
      integer myrank, j, ier, i

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
      end
c$$$!============================================================
c$$$      subroutine wrrestartglobal(timein, ntimein)
c$$$      use p_data
c$$$      use t_data
c$$$      use basic
c$$$      use arrays
c$$$      use diagnostics
c$$$      implicit none
c$$$      include 'mpif.h'
c$$$      integer ifirstrs2, numdofsglobal, mpitag, status, 
c$$$     &     ierr, request, numentsglobal(4), i, numnodes, j, globalid,
c$$$     &     begindofnumber, enddofnumberplusone, ntimein
c$$$      real timein
c$$$      integer recvinfo(8)
c$$$      data ifirstrs2/1/
c$$$
c$$$      if(ifirstrs2 .ne. 1 .and. myrank .eq. 0) then
c$$$         call rename('C1restart', 'C1restarto')
c$$$      endif
c$$$      ifirstrs2 = 0
c$$$      mpitag = 9070
c$$$      do i=1,8
c$$$         recvinfo(i) = i
c$$$      enddo
c$$$      if(myrank .ne. 0) then
c$$$         ! unaligned access warnings on PPPL sgi machines is from communication here
c$$$         call MPI_Recv(recvinfo(2), 4, MPI_INTEGER, myrank-1, 
c$$$     &        mpitag, MPI_COMM_WORLD, status, ierr)
c$$$         if(ierr .ne. MPI_SUCCESS) then
c$$$            write(*,*) 'MPI recv problem in restart'
c$$$            call safestop(828)
c$$$         endif         
c$$$         call MPI_Wait(request, status, ierr)
c$$$         if(ierr .ne. MPI_SUCCESS) then
c$$$            write(*,*) 'MPI wait problem in restart'
c$$$            call safestop(829)
c$$$         endif         
c$$$#ifdef _AIX
c$$$	 print *, "must use wrrestart on ibm instead of wrrestartglobal"
c$$$ 	 call safestop(9229)
c$$$#else
c$$$         open(56,file='C1restart',form='unformatted',status='old',
c$$$     &        ACCESS='append')
c$$$#endif
c$$$         endfile 56
c$$$      else
c$$$         if(ifirstrs2 .ne. 1) call rename('C1restart', 'C1restarto')
c$$$         ifirstrs2 = 0
c$$$         open(56,file='C1restart',form='unformatted',status='unknown')
c$$$c     first put in information to check on run information
c$$$         call numglobaldofs(numvar, numdofsglobal)
c$$$         call numglobalents(numentsglobal(1),numentsglobal(2),
c$$$     &        numentsglobal(3), numentsglobal(4))
c$$$         write(56) numentsglobal(1)!, numentsglobal(2), numentsglobal(3)
c$$$         write(56) numentsglobal(2)
c$$$         write(56) numentsglobal(3)
c$$$         write(56) numdofsglobal
c$$$         write(56) numvar
c$$$         write(56) iper
c$$$         write(56) jper   
c$$$         write(56) ntimein,timein
c$$$         write(56) totcur0,tflux0,gbound,facd,fbound,ptot,vloop
c$$$      endif
c$$$
c$$$      call numnod(numnodes)
c$$$      do j=1,numnodes
c$$$         call entglobalid(j, 0, globalid)
c$$$         write(56) globalid
c$$$         call entdofs(numvar, j, 0, begindofnumber, 
c$$$     &        enddofnumberplusone)
c$$$         do i=begindofnumber,enddofnumberplusone-1
c$$$            write(56) vel(i)
c$$$         enddo
c$$$         do i=begindofnumber,enddofnumberplusone-1
c$$$            write(56) vel0(i)
c$$$         enddo
c$$$         do i=begindofnumber,enddofnumberplusone-1
c$$$            write(56) vels(i)
c$$$         enddo
c$$$         do i=begindofnumber,enddofnumberplusone-1
c$$$            write(56) velold(i)
c$$$         enddo
c$$$         do i=begindofnumber,enddofnumberplusone-1
c$$$            write(56) phi(i)
c$$$         enddo
c$$$         do i=begindofnumber,enddofnumberplusone-1
c$$$            write(56) phi0(i)
c$$$         enddo
c$$$         do i=begindofnumber,enddofnumberplusone-1
c$$$            write(56) phis(i)
c$$$         enddo
c$$$         do i=begindofnumber,enddofnumberplusone-1
c$$$            write(56) phiold(i)
c$$$         enddo
c$$$         call entdofs(1, j, 0, begindofnumber, 
c$$$     &        enddofnumberplusone)
c$$$         do i=begindofnumber,enddofnumberplusone-1
c$$$            write(56) jphi(i)
c$$$         enddo         
c$$$         do i=begindofnumber,enddofnumberplusone-1
c$$$            write(56) vor(i)
c$$$         enddo         
c$$$         do i=begindofnumber,enddofnumberplusone-1
c$$$            write(56) com(i)
c$$$         enddo         
c$$$         do i=begindofnumber,enddofnumberplusone-1
c$$$            write(56) den(i)
c$$$         enddo         
c$$$         do i=begindofnumber,enddofnumberplusone-1
c$$$            write(56) den0(i)
c$$$         enddo         
c$$$         do i=begindofnumber,enddofnumberplusone-1
c$$$            write(56) denold(i)
c$$$         enddo         
c$$$      enddo 
c$$$      i = -1
c$$$c     below specifies that this is the last line
c$$$      if(myrank .eq. maxrank-1) write(56) i
c$$$
c$$$      close(56)
c$$$      if(myrank .ne. maxrank-1) then
c$$$         call MPI_Send(recvinfo(2), 4, MPI_INTEGER, myrank+1,
c$$$     &        mpitag, MPI_COMM_WORLD, ierr)
c$$$         if(ierr .ne. MPI_SUCCESS) then
c$$$            write(*,*) 'MPI send problem in restart'
c$$$            call safestop(827)
c$$$         endif
c$$$      endif
c$$$!      write(9,1201) ntime,time !9 is FILE__C1NEW
c$$$      call MPI_Barrier(MPI_COMM_WORLD, i)
c$$$      
c$$$ 1201 format("* * * restart file written at cycle,time = ",i5,1pe12.4,    &
c$$$     &     "  * * *")
c$$$
c$$$      return
c$$$      end
c$$$!============================================================
c$$$      subroutine rdrestartglobal
c$$$      use p_data
c$$$      use t_data
c$$$      use basic
c$$$      use arrays
c$$$      use diagnostics
c$$$      implicit none
c$$$      include 'mpif.h'
c$$$      integer numdofsglobal, numentsglobal(4), i, numnodes, j, globalid,
c$$$     &     begindofnumber, enddofnumberplusone, in1, in2, in3, nodeid
c$$$      real tempvariable
c$$$
c$$$      open(56,file='C1restart',form='unformatted',status='old')
c$$$      
c$$$c     first put in information to check on run information
c$$$      call numglobaldofs(numvar, numdofsglobal)
c$$$      call numglobalents(numentsglobal(1),numentsglobal(2),
c$$$     &     numentsglobal(3), numentsglobal(4))
c$$$      read(56) in1!, in2, in3
c$$$      read(56) in2
c$$$      read(56) in3
c$$$      write(*,*) 'in',in1, in2, in3
c$$$      if(in1 .ne. numentsglobal(1) .or. in2 .ne. numentsglobal(2) .or.
c$$$     &     in3 .ne. numentsglobal(3)) then
c$$$         write(*,*) 'Mesh does not match for rdrestartglobal', in1, in2, in3,
c$$$     &        numentsglobal(1), numentsglobal(2), numentsglobal(3)
c$$$         call safestop(567)
c$$$      endif
c$$$      read(56) numdofsglobal
c$$$      read(56) numvar
c$$$      read(56) iper
c$$$      read(56) jper   
c$$$      read(56) ntimer,timer
c$$$      read(56) totcur0,tflux0,gbound,facd,fbound,ptot,vloop
c$$$      write(*,*) 'in restartglobal times are ', ntimer, timer
c$$$      call numnod(numnodes)
c$$$      read(56) globalid
c$$$      do while(globalid .ne. -1)
c$$$         call globalidnod(globalid, nodeid)
c$$$         if(nodeid .ne. -1) then ! this node exists on this proc
c$$$            call entdofs(numvar, nodeid, 0, begindofnumber, 
c$$$     &           enddofnumberplusone)
c$$$            do i=begindofnumber,enddofnumberplusone-1
c$$$               read(56) vel(i)
c$$$            enddo
c$$$            do i=begindofnumber,enddofnumberplusone-1
c$$$               read(56) vel0(i)
c$$$            enddo
c$$$            do i=begindofnumber,enddofnumberplusone-1
c$$$               read(56) vels(i)
c$$$            enddo
c$$$            do i=begindofnumber,enddofnumberplusone-1
c$$$               read(56) velold(i)
c$$$            enddo
c$$$            do i=begindofnumber,enddofnumberplusone-1
c$$$               read(56) phi(i)
c$$$            enddo
c$$$            do i=begindofnumber,enddofnumberplusone-1
c$$$               read(56) phi0(i)
c$$$            enddo
c$$$            do i=begindofnumber,enddofnumberplusone-1
c$$$               read(56) phis(i)
c$$$            enddo
c$$$            do i=begindofnumber,enddofnumberplusone-1
c$$$               read(56) phiold(i)
c$$$            enddo
c$$$            call entdofs(1, nodeid, 0, begindofnumber, 
c$$$     &           enddofnumberplusone)
c$$$            do i=begindofnumber,enddofnumberplusone-1
c$$$               read(56) jphi(i)
c$$$            enddo         
c$$$            do i=begindofnumber,enddofnumberplusone-1
c$$$               read(56) vor(i)
c$$$            enddo         
c$$$            do i=begindofnumber,enddofnumberplusone-1
c$$$               read(56) com(i)
c$$$            enddo         
c$$$            do i=begindofnumber,enddofnumberplusone-1
c$$$               read(56) den(i)
c$$$            enddo         
c$$$            do i=begindofnumber,enddofnumberplusone-1
c$$$               read(56) den0(i)
c$$$            enddo         
c$$$            do i=begindofnumber,enddofnumberplusone-1
c$$$               read(56) denold(i)
c$$$            enddo         
c$$$         else ! node is not on this proc
c$$$            do i=1,numvar*6*8 + 6*6 ! the number of dof values at this node
c$$$               read(56) tempvariable
c$$$            enddo
c$$$         endif
c$$$         read(56) globalid
c$$$      enddo 
c$$$      
c$$$      close(56)
c$$$      return
c$$$      end
c$$$!============================================================
