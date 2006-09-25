      subroutine wrrestart
      use p_data
      use t_data
      use basic
      use arrays
      implicit none
      integer mmnn18, j1, mmnn6, ifirstrs
      real ekinemagt
      data ifirstrs/1/
      
      call numdofs(numvar, mmnn18)

      if(ifirstrs .ne. 1) call rename('C1restartout','C1restartouto')
      ifirstrs = 0
      open(56,file='C1restartout',form='formatted',status='unknown')
      do j1=1,mmnn18 
         write(56,1204) j1,vel(j1)
      enddo
      do j1=1,mmnn18       
         write(56,1204) j1, vel0(j1)
      enddo
      do j1=1,mmnn18 
         write(56,1204) j1, vels(j1)
      enddo
      do j1=1,mmnn18 
         write(56,1204) j1, phi(j1)
      enddo
      do j1=1,mmnn18 
         write(56,1204) j1, phi0(j1)
      enddo
      do j1=1,mmnn18 
         write(56,1204) j1, phis(j1)
      enddo
!
      call numdofs(1, mmnn6)
      do j1=1,mmnn6     
         write(56,1204) j1, jphi(j1)
      enddo
      do j1=1,mmnn6       
         write(56,1204) j1, vor(j1)
      enddo
      do j1=1,mmnn6       
         write(56,1204) j1, com(j1)
      enddo
      do j1=1,mmnn6       
         write(56,1204) j1, den(j1)
      enddo
      do j1=1,mmnn6       
         write(56,1204) j1, den0(j1)
      enddo
!
      write(56, 1204) ntime,time
!      write(56) ((graphit(itime,igr),igr=1,maxplots),itime=1,ntime)
      write(56,120) ekin,emag,ekind,emagd,ekinp,emagp,ekinpd,emagpd,
     &     ekint, emagt,ekintd,emagtd,ekinph,ekinth,emagph,emagth,          
     &     ekin3,ekin3d,ekin3h,emag3,emag3d,emag3h,tflux,gbound
      close(56)
      write(*,*) 'emagph is now ',emagph
      write(*,1201) ntime,time
      write(9,1201) ntime,time
 1201 format("* * * restart file written at cycle,time = ",i5,1pe12.4,    &
     &     "  * * *")
 1204 format(I5, 1X, E12.4)
 120  format(24(e12.4, 1x))

      return
      end
!============================================================
      subroutine rdrestart
      use p_data
      use t_data
      use basic
      use arrays
      implicit none
      integer mmnn18, j1, mmnn6, itime, igr, numnodes
      
      call numdofs(numvar, mmnn18)
      call safestop(111)
      open(56,file='restart',form='unformatted',status='old')
      read(56) (vel(j1),j1=1,mmnn18)
      read(56) (vel0(j1),j1=1,mmnn18)
      read(56) (vels(j1),j1=1,mmnn18)
      read(56) (phi(j1),j1=1,mmnn18)
      read(56) (phi0(j1),j1=1,mmnn18)
      read(56) (phis(j1),j1=1,mmnn18)
!
!....DEBUG
!     if(myrank.le.0) then
!     do jj=1,m
!     write(98,1002) jj
!1002 format("jj =",i3)
!     write(98,1001) (vel(j1),j1=jj*18-5,jj*18)
!1001 format(1p6e12.4)
!     enddo
!     endif
!     call safestop(2)
!
!
      call numdofs(1, mmnn6)
      read(56) (jphi(j1),j1=1,mmnn6)
      read(56) (vor(j1),j1=1,mmnn6)
      read(56) (com(j1),j1=1,mmnn6)
      read(56) (den(j1),j1=1,mmnn6)
      read(56) (den0(j1),j1=1,mmnn6)
!
      read(56) ntimer,timer
      read(56) ((graphit(itime,igr),igr=1,maxplots),itime=1,ntimer)
      read(56) ekin,emag,ekind,emagd,ekinp,emagp,ekinpd,emagpd,ekint,   &
     &          emagt,ekintd,emagtd,ekinph,ekinth,emagph,emagth,        &
     &          ekin3,ekin3d,ekin3h,emag3,emag3d,emag3h,tflux,gbound
      close(56)
      write(*,1201) ntimer,timer
      write(9,1201) ntimer,timer
 1201 format("* * * restart file read at cycle,time = ",i5,1pe12.4,     &
     &     "  * * *")
      return
      end
!============================================================

