      subroutine wrrestart
      use p_data
      use t_data
      use basic
      use arrays

      implicit none
      integer mmnn18, j1, mmnn6, ifirstrs, numnodes, numelms 
      real ekinemagt
      data ifirstrs/1/
      character (len=30) :: fname, oldfname

      call createfilename(fname, oldfname)
      call numnod(numnodes)
      call numfac(numelms)
      call numdofs(numvar, mmnn18)

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
      do j1=1,mmnn18 
         write(56) vel(j1)
      enddo
      do j1=1,mmnn18       
         write(56) vel0(j1)
      enddo
      do j1=1,mmnn18 
         write(56) vels(j1)
      enddo
      do j1=1,mmnn18 
         write(56) velold(j1)
      enddo
      do j1=1,mmnn18 
         write(56) phi(j1)
      enddo
      do j1=1,mmnn18 
         write(56) phi0(j1)
      enddo
      do j1=1,mmnn18 
         write(56) phis(j1)
      enddo
      do j1=1,mmnn18 
         write(56) phiold(j1)
      enddo
!
      call numdofs(1, mmnn6)
      do j1=1,mmnn6       
         write(56) den(j1)
      enddo
      do j1=1,mmnn6       
         write(56) den0(j1)
      enddo
      do j1=1,mmnn6       
         write(56) denold(j1)
      enddo
!
      write(56) ntime,time
      write(56) totcur0,tflux0,gbound,facd
      close(56)

      if(maxrank.eq.1) then
         write(9,1201) ntime,time
 1201    format("* * * restart file written at cycle,time = ",i5,1pe12.4,    &
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
      implicit none
      integer mmnn18, j1, mmnn6, itime, igr, numnodes, inumnodes, 
     1     inumelms, immnn18, inumvar, iiper, ijper, imyrank,
     2     imaxrank, numelms
      character (len=30) :: fname, oldfname
      real mytemp(10000)

      call createfilename(fname, oldfname)
      call numdofs(numvar, mmnn18)
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
      do j1=1,mmnn18
         read(56) vel(j1)
      enddo
      do j1=1,mmnn18
         read(56) vel0(j1)
      enddo
      do j1=1,mmnn18
         read(56) vels(j1)
      enddo
      do j1=1,mmnn18
         read(56) velold(j1)
      enddo
      do j1=1,mmnn18
         read(56) phi(j1)
      enddo
      do j1=1,mmnn18
         read(56) phi0(j1)
      enddo
      do j1=1,mmnn18
         read(56) phis(j1)
      enddo
      do j1=1,mmnn18
         read(56) phiold(j1)
      enddo

      call numdofs(1, mmnn6)
      do j1=1,mmnn6
         read(56) den(j1)
      enddo
      do j1=1,mmnn6
         read(56) den0(j1)
      enddo
      do j1=1,mmnn6
         read(56) denold(j1)
      enddo

!
      read(56) ntimer,timer
      read(56) totcur0,tflux0,gbound,facd
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
