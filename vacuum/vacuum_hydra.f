c.............................................................
      subroutine date_time ( date_array, datex, timex,nout1,nout2 )
c..............................................................
c   
c$$$  character*len fdate, d
c$$$  external fdate
c$$$  d = fdate( )
c$$$
c$$$Description
c$$$  The fdate function returns the current date and time as a 24-character
c$$$  string in the format described in ctime(3).  Neither `newline' nor NULL is
c$$$  included.
c$$$
c$$$  The fdate function can be called as a function or as a subroutine.  If
c$$$  called as a function, the calling routine must define the type and length.
c$$$  For example:
c
       character*24   fdate
       character*(*) date_array
       character*(*) datex,timex
       external       fdate
       write(nout1,*) fdate()
       write(nout2,*) fdate()

       date_array = fdate()

c$$$  If fdate is called as a subroutine, the call is as follows:
c$$$
c$$$       character*24   d
c$$$       external       fdate
c$$$       call fdate(d)
c
      return
      end
c
c...........................................................
      subroutine clock ( ntim )
c...........................................................
c
      character*(10) logo(2),ntim,ndate,mach,nsfx
c
      return
      end
c
c............................................................
      subroutine shella
c...........................................................
c
c     
      call system ( 'rm -f modovmc' )
      call system ( 'rm -f pestotv' )
      call system ( 'rm -f vdata' )
      call system ( 'rm -f vac.cgm' )
      call system ( 'rm -f vacadj' )
      call system ( 'rm -f vacdcon' )
      call system ( 'rm -f vacgato' )
c
      return
      end
c
c............................................................
      subroutine shellb
c...........................................................
c
      return
      end
c
c..........................................................
      subroutine gelima  ( copmat,nfm,uvpr,nfm1,jmax1,jmax6,uvp0,nfm2,
     $     wrki,waa,nfm3,wbb,nfm5,ifail )
c..........................................................
c
      call f04aef ( copmat,nfm,uvpr,nfm1,jmax1,jmax6,uvp0,nfm2,wrki,
     $     waa,nfm3,wbb,nfm5,ifail )
c
      return
      end
c
c.........................................................
      subroutine gelimb ( copmat,nfm,uvpwr,nfm1,jmax1,jmax6,uvpw0,
     $     nfm2,wrki,ifail )
c...........................................................
c
      call f04aaf ( copmat,nfm,uvpwr,nfm1,jmax1,jmax6,uvpw0,nfm2,wrki,
     $     ifail )
c
      return
      end

c..................................................................
      SUBROUTINE evalef ( msiz,zvec,nd,zbvc,ndb,zeps1,zalfar,zalfai,
     $     zbeta,lzmatv,zwk,ndzz,izter,izfail )
c..................................................................

      CALL f02bjf( msiz,zvec,nd,zbvc,ndb,zeps1,zalfar,zalfai,
     $     zbeta,lzmatv,zwk,ndzz,izter,izfail )

      END SUBROUTINE evalef

c...............................................
      subroutine timedate(ntim,ndat,mach,nsfx)
c................................................
c
        return
        end
c
c...................................................
        subroutine userinfo(nuser,nacct,ndrop,nsfx)
c..................................................
c
        return
        end
c
c...............................
        subroutine close(iun)
c...............................
c
        close(iun)
c
        return
        end
                                                                             
