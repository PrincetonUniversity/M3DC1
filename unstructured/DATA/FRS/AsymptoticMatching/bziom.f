c...........................................................................
c     Source code for a ZIO emulator. Uses Fortran Input/Output extensions 
c     DEC-ALPHA version
c     J. Manickam 6.29.95
c     Includes entries for 
c     ZOP
c     ZCL
c     ZWR
c     ZRD
c.........................
c     Usage Notes:
c     Uses unformatted read/writes with a recordtype 'stream'. This has
c     the free field format similar to direct access. However it treats the
c     file as a sequential file. There is a function 'seek' which allows a
c     preposition of the read so that this can be made to look like ZIO
c     The file name can be specified as in a typical OPEN call
c...........................................................................
c     
      subroutine zop(ioc,name,nsize,idisk,icode,ilab)
      character*(8) name
c
      open(unit=ioc,file=name,form='unformatted',
     $     status='unknown')
      return
      end
c...................................................................
      subroutine zcl(ioc,ierr)
c...................................................................
c     
      close(ioc)
      return
      end
c...................................................................
      subroutine zwr(ioc,a,nwords,nadres,lgivup,irr)
c...................................................................
c     
      dimension a(1)
      integer*4 iocs,noffsets,fseek,ftell
      nbytes = 8
      iocs = ioc
c     check on location
c     ncurr = ftell(iocs)
      ncurr = 0
      noffsets = nadres * nbytes - ncurr
      ierr = 0
      write(ioc)(a(i),i=1,nwords)
c     
c     Error checks
c     
      if(ierr .eq. 0) return
c     write error code
      write(6,100)ierr      
 100  format(' Error in ZWR, error code = ', i4, 
     $        ' Check man 3f perror ')
c      stop
      end
c...................................................................
      subroutine zrd(ioc,a,nwords,nadres,lgivup,irr)
c...................................................................
c     
      dimension a(1)
      integer*4 iocs,noffsets,fseek,ftell
      nbytes = 8
      iocs = ioc
c     check on location
c     ncurr = ftell(iocs)
      ncurr = 0
      noffsets = nadres * nbytes - ncurr
      ierr = 0
      read(ioc)(a(i),i=1,nwords)
c     
c     Error checks
c     
      if(ierr .eq. 0) return

c        write error code
         write(6,100)ierr
 100     format(' Error in ZWR, error code = ', i4, 
     $        ' Check man 3f perror ')
         stop
         end

	subroutine skipeof(iva,iva1)
	return
	end


	subroutine userinfo(nuser,nacct,ndrop,nsfx)
	return
	end

	subroutine  close(iun)
	close(iun)
	return
	end
