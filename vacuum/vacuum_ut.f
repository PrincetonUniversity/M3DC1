c
c......................................................................
      subroutine errmes ( nout,mesage )
c......................................................................
c
      character*(*) mesage
c
      write(5,100)mesage
      write(nout,100)mesage
  100 format(" error in subroutine",1x, a10/)
c
c      call exit
      stop
      end
c
c...................................................
c
c...............................................................
ccc
      subroutine wrtout ( nsrf, vec, wvec, m1,m2 )
c
      dimension vec(*)
      character*(*) wvec
c
      write ( 16, 100 ) nsrf, wvec
  100 format ( //,1x, 'nsrf = ',i3, 2x, a, / )
c
      write ( 16, 200 ) ( vec(i), i = m1,m2 )
  200 format ( 1x, 1p10e13.5 )
c
      return
      end
ccc
c..............................................................
      subroutine vecwrt ( mm, vec, wvec, m1,m2, nout1,nout2 )
c...............................................................
c
      dimension vec(*)
      character*(*) wvec
c
      if ( nout1 .ne. 0 ) then
         write ( nout1, 100 ) wvec, mm
         write ( nout1, 200 ) ( vec(i), i = m1,m2 )
      end if
c
      if ( nout2 .ne. 0 ) then
         write ( nout2, 100 ) wvec, mm
         write ( nout2, 200 ) ( vec(i), i = m1,m2 )
      end if
c
  100 format ( /,1x, "VECTOR ELEMENTS of ", a, 2x, "mm = ",i3 )
  200 format ( 1p10e12.4 )
c
      return
      end
ccc
c...............................................................
      subroutine writg1 ( string, seps, iseps, nout ) 
c................................................................
ccc
      dimension nout(*)
      character*(*) string
      character*(*) seps
ccc
      do io = 1, 3
         if ( nout(io) .gt. 0 )then
            if ( iseps .eq. 0 )
     $           write ( nout(io), '(/,a)' ) string
            if ( iseps .eq. 1 )
     $           write ( nout(io), '(/, a,/, 5x, a)' ) seps, string
            if ( iseps .eq. 2 )
     $           write ( nout(io), '(/, 5x,a, /,a )' ) string, seps
            if ( iseps .eq. 3 )
     $           write ( nout(io), '( /, a ,/, 5x,a,/,a )' )
     $           seps, string, seps
         end if
      enddo
ccc
      return
      end
ccc
      subroutine errnc ( mesage, nrcode, nout )
ccc
      character*(*) mesage
      integer nout(*)
ccc
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) mesage, nrcode
         end if
      enddo
ccc
 100  format ( "nrcode in ", a16, " = ", i5 )
ccc
      return
      end
ccc
      subroutine nwopn ( cdfid, file, nrw, rcode, nout )
ccc
      integer cdfid, rcode, nout(*)
      character*(*) file, nrw
ccc
      if ( rcode .gt. 0 ) call errnc ( "ncopn", rcode, nout )
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) cdfid, file, nrw
         end if
      enddo
ccc
 100  format ( /,"cdfid = ", i3, /,
     $     "ncopn called for file, ",a, " ---nrw = ", a )
ccc
      return
      end
ccc
      subroutine nwinq ( cdfid, ndims, nvars, ngatts, recid, rcode,
     $     nout )
ccc
      integer cdfid, rcode, recid, nout(*)
ccc
      if ( rcode .gt. 0 ) call errnc ( "ncinq", rcode, nout )
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) cdfid, ndims, nvars,
     $           ngatts, recid
         end if
      enddo
ccc
 100  format ( /, "ncinq: cdfid = ", i3, /,
     $     "number of dimensions, ndims is ...", i5, /,
     $     "number of variables, nvars is ..." i5, /,
     $     "number of global attributes, ngatts is ...", i5, /,
     $     "ID of the unlimited dimension, recid is ...", i5 )
ccc
c     return
      end
ccc
      subroutine nwdid ( cdfid, dimnam, ndid, rcode, nout ) 
ccc
      integer cdfid, rcode, nout(*)
      character*(*) dimnam
ccc      
      if ( rcode .gt. 0 ) call errnc ( "ncdid", rcode, nout )
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) cdfid, dimnam, ndid
         end if
      enddo
ccc
 100  format ( /,"ncdid: cdfid = ", i3,
     $     " dimension name is ...", a, /, "dimension id is ...",i5 )
ccc
      return
      end
ccc
      subroutine nwdinq ( cdfid, dimid, dimnam, dimsiz, rcode,
     $     nout )
ccc
      integer cdfid, rcode, dimid, dimsiz, nout(*)
      character*(*) dimnam
ccc      
      if ( rcode .gt. 0 ) call errnc ( "ncdinq", rcode, nout )
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) cdfid, dimnam, dimid, dimsiz
         end if
      enddo
ccc
 100  format ( /,"ncdinq: cdfid = ", i3,
     $     ".  Dimension name, id, and size = ",/ ,a, 2i5 )
ccc
      return
      end
ccc
      subroutine nwvid ( cdfid, varnam, vid, rcode, nout )
ccc
      integer cdfid, rcode, vid, nout(*)
      character*(*) varnam
ccc      
      if ( rcode .gt. 0 ) call errnc ( "ncvid", rcode, nout )
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) cdfid, varnam, vid
         end if
      enddo
ccc
 100  format ( /,"ncvid: cdfid = ", i3, " varnam= ",a ,/,
     $     "variable id is ...", i5 )
ccc
      return
      end
ccc
      subroutine nwvinq ( cdfid, vid, varnam, vartyp, nvdims,
     $     vdims, nvatts, rcode, nout )
ccc
      integer cdfid, rcode, vid, vartyp, vdims(*), nout(*)
      character*(*) varnam
      character*8 fmt
      character*72 fmt1
      character*72 fmt2
      character*72 fmtn(10)
ccc
      call vfrmt2 ( nvdims, "i5", fmt )
      fmt1 = '/,"ncvinq: cdfid = ", i3, ".  Vid= ",i5,/,'
      fmt2 = '( "Vector of Dimension ID''s is ...",'//fmt//')'
      fmtn(1) = '('//fmt1
      fmtn(2) = '"variable type and name is ...", i5,'
      fmtn(3) = '20x, ".....  ", a16, /,'
      fmtn(4) = '"no. of dimensions and attributes are ... ", 2i5 )'
ccc      
      if ( rcode .gt. 0 ) call errnc ( "ncvinq", rcode, nout )
ccc
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), fmtn ) cdfid, vid, vartyp, varnam,
     $        nvdims, nvatts
            if ( nvdims .gt. 0 ) then
               write ( nout(io), fmt2 ) (vdims(j),j=1,nvdims)
            end if
         end if
      enddo
ccc
      return
      end
ccc
      subroutine nwvg1i ( cdfid, vid, vindx,nd, val1, rcode,
     $     nout )
ccc
      integer cdfid, rcode, vid, vindx(*), val1, nout(*)
      character*8 fmt
      character*72 fmt1
      character*72 fmtn(10)
ccc
      if ( rcode .gt. 0 ) call errnc ( "ncvg1i", rcode, nout )
ccc
      if (  nd .eq. 0 ) nd = 1
ccc
      call vfrmt2 ( nd, "i5", fmt )
      fmt1 = '/,"ncvgli: cdfid = ", i3, " vid= ",i5,/,'
      fmtn(1) = '('//fmt1
      fmtn(2) = '"index of the integer value is ...",'
     $     //fmt//',/,'
      fmtn(3) = '"value of the integer is ...",i5 )'
ccc      
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), fmtn ) cdfid, vid,
     $           (vindx(i),i=1,nd), val1
         end if
      enddo
ccc
      return
      end
ccc
      subroutine nwvgt1 ( cdfid, vid, vindx,nd, val1, rcode,
     $     nout )
ccc
      integer cdfid, rcode, vid, vindx(*), nout(*)
      character*8 fmt
      character*72 fmt1
      character*72 fmtn(10)
ccc
      if ( rcode .gt. 0 ) call errnc ( "ncvgt1", rcode, nout )
ccc
      if ( nd .eq. 0 ) nd = 1
ccc
      call vfrmt2 ( nd, "i5", fmt )
      fmt1 = '/,"ncvgt1: cdfid = ", i3, " vid= ",i5,/,'
      fmtn(1) = '('//fmt1
      fmtn(2) = '"index of the data value is ...",'
     $     //fmt//',/,'
      fmtn(3) = '"value of the data is ...",1pe12.4 )'
ccc   
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), fmtn ) cdfid, vid,
     $           (vindx(i),i=1,nd), val1
         end if
      enddo
ccc
      return
      end
ccc
      subroutine nwvg1c ( cdfid, vid, vindx,nd, chval, rcode,
     $     nout )
ccc
      integer cdfid, rcode, vid, vindx(*), nout(*)
      character chval
      character*8 fmt
      character*72 fmt1
      character*72 fmtn(10)
ccc
      call vfrmt2 ( nd, "i5", fmt )
      fmt1 = '/,"ncvg1c: cdfid = ", i3, " vid= ",i5,/,'
      fmtn(1) = '('//fmt1
      fmtn(2) = '"index of the character value is ...",'
     $     //fmt//',/,'
      fmtn(3) = '"value of the character is ... ",a )'
ccc
      if ( rcode .gt. 0 ) call errnc ( "ncvg1c", rcode, nout )
ccc
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), fmtn ) cdfid, vid,
     $           (vindx(j),j=1,nd), chval
         end if
      enddo
ccc      
      return
      end
ccc
      subroutine nwvgti ( cdfid, vid, start, count,nd, val0, vals,
     $     rcode, nd1,nd2, nout )
ccc
ccc..... Repack only 2D matrix system for now.
ccc
      integer cdfid, rcode, vid, start(*), count(*),
     $     val0(*), vals(nd1,nd2), nout(*)
ccc
      do i2 = 1, count(2)
         do i1 = 1, count(1)
            vals(i1,i2) = val0( i1+(i2-1)*count(1) ) 
         enddo
      enddo      
ccc      
      if ( rcode .gt. 0 ) call errnc ( "ncvgt", rcode, nout )
ccc
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) cdfid, vid,
     $           (start(i),i=1,5), (count(j),j=1,5)
         end if
      enddo
ccc
 100  format ( /,"ncvgt1: cdfid = ", i3, " vid= ",i5,/,
     $     "start, first corner of hyperslab is ...",5i5,/,
     $     "count, edge lengths of the hyperslab is ...", 5i5 )
ccc
      call writis ( "Hyperslab of integers is..", vals, nd1,nd2,
     $     count, nd, nout )
ccc
      return
      end
ccc
      subroutine nwvgt ( cdfid, vid, start,count,nd, val0, vals,
     $     nd1,nd2, rcode, nout )
ccc
ccc.....  Repack only 1D and 2D matrix system for now.
ccc
      integer cdfid, rcode, vid, start(*), count(*), nout(*)
      character*8 fmt
      character*72 fmt1
      character*72 fmtn(10)
      dimension vals(nd1,nd2), val0(*)
ccc
      if ( nd .eq. 1 ) count(2) = 1
      do i2 = 1, count(2)
         do i1 = 1, count(1)
            vals(i1,i2) = val0( i1+(i2-1)*count(1) )
         enddo
      enddo
ccc
      call vfrmt2 ( nd, "i5", fmt )
      fmt1 = '/,"ncvgt: cdfid = ", i3, " vid= ",i5,/,'
      fmtn(1) = '('//fmt1
      fmtn(2) = '"start, first corner of hyperslab is ...",'
     $     //fmt//',/,'
      fmtn(3) = '"count, edge lengths of the hyperslab is ...",'
     $     //fmt//')'
ccc      
      if ( rcode .gt. 0 ) call errnc ( "ncvgt", rcode, nout )
ccc
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), fmtn ) cdfid, vid,
     $           (start(i),i=1,nd), (count(j),j=1,nd)
         end if
      enddo
ccc
      call writvs ( "Hyperslab of Values is..", vals, nd1,nd2,
     $     count, nd, nout )
ccc
      return
      end
ccc
      subroutine nwvgtc ( cdfid, vid, start,count,nd, string, lenstr,
     $     rcode, nout )
ccc
      integer cdfid, rcode, vid, start(*), count(*), nout(*)
      character*(*) string
      character*8 fmt
      character*72 fmt1
      character*72 fmtn(10)
ccc
      call vfrmt2 ( nd, "i5", fmt )
      fmt1 = '/,"ncvgtc: cdfid = ", i3, " vid= ",i5,/,'
      fmtn(1) = '('//fmt1
      fmtn(2) = '"start, first corner of hyperslab is ...",'
     $     //fmt//',/,'
      fmtn(3) = '"count, edge lengths of the hyperslab is ...",'
     $     //fmt//',/,'
      fmtn(4) = '"string, character string is ...",/, a )'
ccc      
      if ( rcode .gt. 0 ) call errnc ( "ncvgtc", rcode, nout )
ccc
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), fmtn ) cdfid, vid,
     $           (start(i),i=1,nd),(count(j),j=1,nd), string
         end if
      enddo
ccc
      return
      end
ccc
      subroutine nwainq ( cdfid, vid, attnam, attype, attlen, rcode,
     $     nout )
ccc
      integer cdfid, rcode, vid, attype, attlen, nout(*)
      character*(*) attnam
ccc
      if ( rcode .gt. 0 ) call errnc ( "ncainq", rcode, nout ) 
ccc
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) cdfid, vid, attnam,
     $           attype, attlen
         end if
      enddo
ccc
 100  format ( /,"ncainq: cdfid = ", i3, ". Vid = ",i5,/,
     $     "attribute's name, type, and length: ",/,"... ", a,
     $     "... ", 2i5 )
ccc
      return
      end
ccc
      subroutine nwagtc ( cdfid, vid, attnam, astrng, attlen, maxa1,
     $     rcode, nout )
ccc
      integer cdfid, rcode,vid, attlen, nout(*)
      character*(*) attnam, astrng
      character*8 fmt
      character*72 fmt1
      character*72 fmtn(10)
ccc
      call vfrmt3 ( attlen, "a", fmt )
ccc      fmt = 'a'
      fmt1 = '/,"ncagtc: cdfid = ", i3,/,'
      fmtn(1) = '('//fmt1
      fmtn(2) = '"attribute''s name, Vid, length, and value: ",/,'
      fmtn(3) = '"... ", a," ...", 2i5,/,'
      fmtn(4) = '">> ",'//fmt//' " <<" )'
ccc
      if ( rcode .gt. 0 ) call errnc ( "ncatgc", rcode, nout )
ccc
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
ccc            write ( nout(io), '("format, fmt = ", a)' ) fmt
            if ( attlen .gt. maxa1 ) then
               write ( nout(io),50 ) attnam, attlen, maxa1
               return
            else
               write ( nout(io), fmtn ) cdfid, attnam, vid, attlen,
     $              astrng
            end if
         end if
      enddo
ccc
 50   format ( /, "*** ", a, "'s length too long: ", i4, " >",i4 )
 100  format ( /,"ncagtc: cdfid =  ", i3,
     $     " attribute's Vid, name and value  = ", i5,/,
     $     "... ", a," ...",/,
     $     ">> ", a, " <<",/,
     $     "attribute length is ... ", i5 )
ccc
      return
      end
ccc
      subroutine nwagt ( cdfid, vid, attnam, attval, attlen,
     $     rcode, nout )
ccc
      integer cdfid, rcode,vid, attlen, nout(*)
      character*(*) attnam
      character*8 fmt
      character*72 fmt1
      character*72 fmtn(10)
ccc
      call vfrmt2 ( attlen, "f13.5", fmt )
      fmt1 = '/,"ncagt: cdfid = ", i3,'
      fmtn(1) = '('//fmt1
      fmtn(2) = '" attribute''s name, Vid, length, and value: ",/,'
      fmtn(3) = '"... ", a," ...", 2i5,/,'
      fmtn(4) = '">> ",'//fmt//' " <<" )'
ccc
      if ( rcode .gt. 0 ) call errnc ( "ncatgc", rcode, nout )
ccc
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
               write ( nout(io), fmtn ) cdfid, attnam, vid, attlen,
     $           attval
         end if
      enddo
ccc
 100  format ( /,"ncagt: cdfid =  ", i3,
     $     " attribute's id, name and value  = ", i5,/,
     $     "... ", a," ...",/,
     $     ">> ", 1pe12.4, " <<",/,
     $     "attribute length is ... ", i5 )
ccc
      return
      end
ccc
      subroutine nwanam ( cdfid, vid, n, attnam, rcode, nout )
ccc
      integer cdfid, rcode,vid, nout(*)
      character*(*) attnam
ccc
      if ( rcode .gt. 0 ) call errnc ( "nwanam", rcode, nout )
ccc
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) cdfid, vid, n, attnam
         end if
      enddo
ccc
 100  format ( /,"ncanam: cdfid = ", i3,
     $     ".  Attribute's variable ID is ... ", i5,/,
     $     "attribute's no. and name is ... ",
     $     i4," -- ", a )
ccc
      return
      end
ccc
ccc
ccc
      subroutine writn1 ( mesage, n, nout )
ccc
      character*(*) mesage
      integer nout(*)
ccc
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) mesage, n
         end if
      enddo
ccc
 100  format ( a, i5 ) 
c
      return
      end
ccc
      subroutine writnn ( mesage, n, nd, nout )
ccc
      character*(*) mesage
      dimension n(*), nout(*)
ccc
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) mesage, (n(i), i=1,nd)
         end if
      enddo
ccc
 100  format ( a, 5i5 ) 
c
      return
      end
ccc
      subroutine writa1 ( mesage, string1, nout )
ccc
      character*(*) mesage
      character*(*) string1
      integer nout(*)
ccc
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) mesage, string1
         end if
      enddo
ccc
 100  format ( a, a ) 
c
      return
      end
ccc
      subroutine writv1 ( mesage, a, nout )
ccc
      character*(*) mesage
      integer nout(*)
ccc
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 100 ) mesage, a
         end if
      enddo
ccc
 100  format ( a, 1pe12.4 ) 
c
      return
      end
ccc
      subroutine writis ( mesage, a, nd1, nd2, nw, nd, nout )
ccc
ccc... writes out only the matrix for the first two indices for now.
      character*(*) mesage
      integer a(nd1,nd2), nout(*), nw(*)
ccc
  10  format (/, a )
  50  format ( 10(1x,10i5,/) )
ccc
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), 10 ) mesage
            do 100 j1 = 1, nw(1)
               write ( nout(io), 50 ) ( a(j1,j2),j2 = 1, nw(2) )
 100        continue
         end if
      enddo
ccc
      return
      end
ccc
      subroutine writvs ( mesage, a, nd1, nd2, nw, nd, nout )
ccc
ccc..... writes out only the matrix for the first two indices for now.
ccc
      character*(*) mesage
      dimension a(nd1,nd2), nw(*)
      integer nout(*)
ccc
  10  format (/, a )
  50  format ( 10(1x,5e13.5,/) )
ccc
      nw1 = nw(1)
      nw2 = nw(2)
      if ( nd .eq. 1 ) then
         nw1 = nw(2)
         nw2 = nw(1)
      end if
      do io = 1, 3
         if ( nout(io) .gt. 0 ) then
            write ( nout(io), '(/,a)' ) mesage
            do 100 j1 = 1, nw1
               if ( nd .ne. 1 ) then
               write ( nout(io), 50 ) ( a(j1,j2),j2 = 1, nw2 )
               else
                  write ( nout(io), 50 ) ( a(j11,1), j11 = 1, nw2 )
               end if
 100        continue
         end if
      enddo
ccc
      return
      end
ccc
      subroutine vfrmt1 ( n, field1, fmt )
ccc
      character*(*) fmt, field1
      character*8 nn
ccc
      write ( nn, 10 ) n
 10   format ( i4 )
ccc
      do while ( index(nn,' ') .eq. 1 )
         nn = nn(2:)
      enddo
ccc
      fmt = '('//nn(1:index(nn, ' ')-1)//field1//')'
ccc
      return
      end
ccc
      subroutine vfrmt2 ( n, field1, fmt )
ccc
      character*(*) fmt, field1
      character*8 nn
ccc
      write ( nn, 10 ) n
 10   format ( i4 )
ccc
      do while ( index(nn,' ') .eq. 1 )
         nn = nn(2:)
      enddo
ccc
      fmt = nn(1:index(nn, ' ')-1)//field1
ccc
      return
      end
ccc
      subroutine vfrmt3 ( n, field1, fmt )
ccc
      character*(*) fmt, field1
      character*8 nn
ccc
      write ( nn, 10 ) n
 10   format ( i4 )
ccc
      do while ( index(nn,' ') .eq. 1 )
         nn = nn(2:)
      enddo
ccc
      fmt = field1//nn(1:index(nn, ' ')-1)
ccc
      return
      end
c
c........................................................
      SUBROUTINE bounds ( x1,z1, n1,n2, xmin,xmax, zmin,zmax )
c.........................................................

      DIMENSION x1(*), z1(*)

      IF ( n1 /= n2 ) THEN
         xmin = MINVAL ( x1(n1:n2) )
         xmax = MAXVAL ( x1(n1:n2) )
         zmin = MINVAL ( z1(n1:n2) )
         zmax = MAXVAL ( z1(n1:n2) )
      ELSE
         xmin = x1(n1)
         xmax = x1(n1)
         zmin = z1(n1)
         zmax = z1(n1)
      END IF
c
c$$$      xmax = -1.0e30
c$$$      zmax = -1.0e30
c$$$      xmin =  1.0e30
c$$$      zmin =  1.0e30
c$$$c
c$$$      do 100 i = n1, n2
c$$$c
c$$$      if ( xmin .gt. x1(i) ) xmin = x1(i)
c$$$      if ( zmin .gt. z1(i) ) zmin = z1(i)
c$$$      if ( xmax .lt. x1(i) ) xmax = x1(i)
c$$$      if ( zmax .lt. z1(i) ) zmax = z1(i)
c$$$c
c$$$  100 continue
c
      RETURN
      END SUBROUTINE bounds
c     
c........................................................
      SUBROUTINE boundsi ( x1,z1, n1,n2, xmin,xmax, zmin,zmax,
     $     ixn,ixx,izn,izx )
c.......................................................
     
      REAL, DIMENSION(*) :: x1, z1
      INTEGER, DIMENSION(1) :: izxn, izxx, izzn, izzx

      xmin = MINVAL ( x1(n1:n2) )
      xmax = MAXVAL ( x1(n1:n2) )
      zmin = MINVAL ( z1(n1:n2) )
      zmax = MAXVAL ( z1(n1:n2) )

      izxn = MINLOC ( x1(n1:n2) )
      izxx = MAXLOC ( x1(n1:n2) )
      izzn = MINLOC ( z1(n1:n2) )
      izzx = MAXLOC ( z1(n1:n2) )

      ixn = izxn(1)
      ixx = izxx(1)
      izn = izzn(1)
      izx = izzx(1)

c$$$      xmax = -1.0e30
c$$$      zmax = -1.0e30
c$$$      xmin =  1.0e30
c$$$      zmin =  1.0e30
c$$$c     
c$$$      do 100 i = n1, n2
c$$$c     
c$$$         if ( xmin .gt. x1(i) ) then 
c$$$            xmin = x1(i)
c$$$            ixn = i
c$$$         end if
c$$$         if ( zmin .gt. z1(i) ) then
c$$$            zmin = z1(i)
c$$$            izn = i
c$$$         end if
c$$$         if ( xmax .lt. x1(i) ) then
c$$$            xmax = x1(i)
c$$$            ixx = i
c$$$         end if
c$$$         if ( zmax .lt. z1(i) ) then
c$$$            zmax = z1(i)
c$$$            izx = i
c$$$         end if
c$$$c     
c$$$ 100  continue
c     
      RETURN
      END SUBROUTINE boundsi
c     
c...........................................................               
      subroutine macopy ( ain,ndin, aout,ndout )
c..........................................................
c
      dimension ain(ndin,ndin), aout(ndout,ndout)
c
      do 20 i = 1, ndin
         do 10 j = 1, ndin
            aout(i,j) = ain(i,j)
 10      continue
 20   continue
c
      return
      end
c
c.......................................................
      subroutine matwrt ( a, maxj1, maxj2, jmax1, jmax6, label )
c...................................................................
c     
c     writes out the elements of matrix a.  label must be
c     20 letters long.
c     
c     integer label ( 2 )
c     
      character*(*) label
      integer out
      dimension a ( maxj1, maxj2 )
c     
      out = 23
c     
      write ( out, 10 ) label
 10   format ( //, 5x, "matrix elements of  ", a )
c     
      do 30 j1 = 1, jmax1
         write ( out, 20 ) j1, ( a(j1,j2),j2 = 1, jmax6 )
 20      format ( i4,10(1x,1p10e11.4,/) )
 30   continue
c     
      return
      end
c     
c..............................................................
      SUBROUTINE vacasym ( vacin, nd,nnj,label, nout1,nout2 )
c...............................................................
c     
c     include 'vacuum1.inc'
c     
      character*(*) label
      dimension vacin(nd,nd)
c     
c.......calculate asymmetries...
c     
      vacrd = 0.0
      vacrso = 0.0
      vacrsd = 0.0
c     
      write ( nout1, 111 ) label
      write ( nout2, 111 ) label
 111  format (/, 1x,"       Asymmetries in ", a )
c     
      do 110 l1 = 1, nnj
         vacrsd = vacrsd + vacin(l1,l1)
         do 108 l2 = 1, l1
            vacrd = vacrd + abs ( vacin(l1,l2) - vacin(l2,l1) )
            if ( l1 .ne. l2 )
     .           vacrso = vacrso + abs ( vacin(l1,l2)+vacin(l2,l1) )
 108     continue
 110  continue
c
      asymv2 = vacrd / ( vacrso + vacrsd )
c
      if ( abs(vacrso) .lt. 1.0e-20 ) then
         vacrso = 1.0
         write (nout1, '("Caution: sum of off-diagonals .lt. 1.e-20")' )
         write (nout2, '("Caution: sum of off-diagonals .lt. 1.e-20")' )
      end if
c
      asymv1 = vacrd / vacrso
c     
      write ( nout1,112 ) vacrsd, vacrso, vacrd, asymv1,asymv2
      write ( nout2,112 ) vacrsd, vacrso, vacrd, asymv1,asymv2
 112  format(  /,1x," sum of diagonals, sd =            ",1pe12.4,
     $     /,1x," sum of off-diagonals, sod =       ",1pe12.4,
     .     /,1x," difference in off diagonals, dod =",1pe12.4,
     .     //,1x," dod/sod =                         ",1pe12.4,
     $     /,1x," dod/(sd+sod) =                    ",1pe12.4,/ )
c     
      RETURN
      END

c..............................................................
      SUBROUTINE vacasyma ( vacin, nd,nnj, zf,label, nout1,nout2 )
c...............................................................
c     
c     Calculate asymmetries of VACIN, normalized to the sum of the 
c     abs(diagonals) and also to the sum of abs(diag) and abs (off diag).
c     zf is  1.0 if VACIN is supposed to be symmetric. Input
c     zf is -1.0 if VACIN is supposed to be antisymmetric.Input

      CHARACTER*(*) label
      DIMENSION vacin(nd,nd)

c.......calculate asymmetries...
c     
      vacrd = 0.0
      vacrso = 0.0
      vacrsd = 0.0
c     
      write ( nout1, 111 ) label, zf
      write ( nout2, 111 ) label, zf
 111  FORMAT (/, 1x,"       Asymmetries in ", a, " -- ZF = ", f5.2 )
c     
      DO 110 l1 = 1, nnj
         vacrsd = vacrsd + abs( vacin(l1,l1) )
         DO 108 l2 = 1, l1
            vacrd = vacrd + abs ( vacin(l1,l2) - zf*vacin(l2,l1) )
            IF ( l1 .ne. l2 )
     .           vacrso = vacrso + abs ( vacin(l1,l2) )
     $           + abs ( vacin(l2,l1) )
 108     CONTINUE
 110  CONTINUE
c
      asymv2 = vacrd / ( vacrso + vacrsd )
c
      IF ( abs(vacrso) .lt. 1.0e-20 ) THEN
         vacrso = 1.0
         WRITE (nout1, '("Caution: sum of off-diagonals .lt. 1.e-20")' )
         WRITE (nout2, '("Caution: sum of off-diagonals .lt. 1.e-20")' )
      END IF
c
      asymv1 = vacrd / vacrso
c     
      WRITE ( nout1,112 ) vacrsd, vacrso, vacrd, asymv1,asymv2
      WRITE ( nout2,112 ) vacrsd, vacrso, vacrd, asymv1,asymv2
 112  FORMAT(  /,1x," sum of diagonals, sd =            ",1pe12.4,
     $     /,1x," sum of off-diagonals, sod =       ",1pe12.4,
     .     /,1x," difference in off diagonals, dod =",1pe12.4,
     .     //,1x," dod/sod =                         ",1pe12.4,
     $     /,1x," dod/(sd+sod) =                    ",1pe12.4,/ )
c     
      RETURN
      END

c..............................................................
      SUBROUTINE masym0 ( vacin, nd,nnj,n1,n2,n10,n20, zf,
     $     label, nout1,nout2 )
c...............................................................

      CHARACTER*(*) label
      DIMENSION vacin(nd,nd)

c     Calculate asymmetries of VACIN, normalized to the sum of the 
c     abs(diagonals) and also to the sum of abs(diag) and abs (off diag).
c     It ignores  lla = n10 and llb = n20 elements.
c     zf is  1.0 if VACIN is supposed to be symmetric. Input
c     zf is -1.0 if VACIN is supposed to be antisymmetric.Input

c     
      vacrd = 0.0
      vacrso = 0.0
      vacrsd = 0.0
c     
      write ( nout1, 111 ) label, zf, n10, n20
      write ( nout2, 111 ) label, zf, n10, n20
 111  format (/, 1x,"  Asymmetries in ", a, " -- ZF = ", f5.2,
     $     " -- Ignores la,lb = ", 2i3 )
c     
      do 110 l1 = 1, nnj
         lla = l1 - 1 + n1
         if ( lla .eq. n10 ) go to 110
         vacrsd = vacrsd + abs( vacin(l1,l1) )
         do 108 l2 = 1, l1
            llb =  l2 - 1 + n2
            if ( llb .eq. n20 ) go to 108
            vacrd = vacrd + abs ( vacin(l1,l2) - zf*vacin(l2,l1) )
            if ( l1 .ne. l2 )
     .           vacrso = vacrso + abs ( vacin(l1,l2) )
     $           + abs( vacin(l2,l1) )
 108     continue
 110  continue
c     
      asymv2 = vacrd / ( vacrso + vacrsd )
c
      if ( abs(vacrso) .lt. 1.0e-20 ) then
         vacrso = 1.0
         write (nout1, '("Caution: sum of off-diagonals .lt. 1.e-20")' )
         write (nout2, '("Caution: sum of off-diagonals .lt. 1.e-20")' )
      end if
c
      asymv1 = vacrd / vacrso
c     
      write ( nout1,112 ) vacrsd, vacrso, vacrd, asymv1,asymv2
      write ( nout2,112 ) vacrsd, vacrso, vacrd, asymv1,asymv2
 112  format(  /,1x," sum of diagonals, sd =            ",1pe12.4,
     $     /,1x," sum of off-diagonals, sod =       ",1pe12.4,
     .     /,1x," difference in off diagonals, dod =",1pe12.4,
     .     //,1x," dod/sod =                         ",1pe12.4,
     $     /,1x," dod/(sd+sod) =                    ",1pe12.4,/ )
c     
      RETURN
      END
c
c..............................................................
      subroutine vacasymi ( vacin, nd,nnj,label, nout1,nout2 )
c...............................................................
c     
      include 'vacuum1.inc'
c     
      character*(*) label
      dimension vacin(nd,nd), wrkrd(nfm),wrkrso(nfm),wrkrsd(nfm),
     $     zsymv1(nfm),zsymv2(nfm)
c     
c.......calculate asymmetries of nested sub-matrices...
c     
      do 10 j = 1, nnj
         wrkrd(j) = 0.0
         wrkrso(j) = 0.0
         wrkrsd(j) = 0.0
 10   continue
c     
      write ( nout1, 111 ) label
      write ( nout2, 111 ) label
 111  format (/, 1x,"  Asymmetries in nested sub-matrices of ", a )
c     
      ncnt = 0
      do 110 l1 = 1, nnj/2 + 1
         ja = l1
         jb = nnj - l1 + 1
         jab = jb - ja + 1
         if ( jab .lt. 2 ) go to 110
         do 108 jin1 = 1, jab
            ja1 = ja + jin1 - 1
            wrkrsd(l1) = wrkrsd(l1) + vacin(ja1,ja1)
            do 106 jin2 = 1, jin1
               ja2 = ja + jin2 - 1
               wrkrd(l1) = wrkrd(l1) +
     $              abs ( vacin(ja1,ja2) - vacin(ja2,ja1) )
               if ( ja1 .ne. ja2 )
     $              wrkrso(l1) = wrkrso(l1) +
     $              abs ( vacin(ja1,ja2)+vacin(ja2,ja1) )
 106        continue
 108     continue
         ncnt = ncnt + 1
 110  continue
c     
      do 150 j = 1, ncnt
         zsymv1(j) = wrkrd(j) / wrkrso(j)
         zsymv2(j) = wrkrd(j) / ( wrkrso(j) + wrkrsd(j) )
 150  continue
c     
      write ( nout1,'("ncnt = ",i3)') ncnt
      write ( nout2,'("ncnt = ",i3)') ncnt
      write ( nout1,221 ) ( wrkrsd(i), i = 1, ncnt )
      write ( nout1,222 ) ( wrkrso(i), i = 1, ncnt )
      write ( nout1,223 ) ( wrkrd(i), i = 1, ncnt )
      write ( nout1,224 ) ( zsymv1(i), i = 1, ncnt )
      write ( nout1,225 ) ( zsymv2(i), i = 1, ncnt )
c     
      write ( nout2,221 ) ( wrkrsd(i), i = 1, ncnt )
      write ( nout2,222 ) ( wrkrso(i), i = 1, ncnt )
      write ( nout2,223 ) ( wrkrd(i), i = 1, ncnt )
      write ( nout2,224 ) ( zsymv1(i), i = 1, ncnt )
      write ( nout2,225 ) ( zsymv2(i), i = 1, ncnt )
c     
 221  format(/,1x," sum of diagonals, sd =            ",/,(1p8e12.4))
 222  format(1x," sum of off-diagonals, sod =       ",/,(1p8e12.4))
 223  format(1x," difference in off diagonals, dod =",/,(1p8e12.4))
 224  format(1x," dod/sod =                         ",/,(1p8e12.4))
 225  format(1x," dod/(sd+sod) =                    ",/,(1p8e12.4))
c     
      return
      end
c
c$$$c......................................................................
c$$$      subroutine timer ( nout,mesage )
c$$$c......................................................................
c$$$c
c$$$c     dimension mesage(2)
c$$$      character*(*) mesage
c$$$ccc
c$$$      return
c$$$ccc
c$$$ccc      call timeused(i1,i2,i3)
c$$$      x1 = i1 / 1000000.
c$$$      x2 = i2 / 1000000.
c$$$      x3 = i3 / 1000000.
c$$$      xdif1 = x1 - x1old
c$$$      xdif2 = x2 - x2old
c$$$      xdif3 = x3 - x3old
c$$$      x1old = x1
c$$$      x2old = x2
c$$$      x3old = x3
c$$$      write(nout,200)
c$$$      write ( iotty,200 )
c$$$      write(nout,100)mesage,xdif1,x1,xdif2,x2,xdif3,x3
c$$$      write(iotty,100)mesage,xdif1,x1,xdif2,x2,xdif3,x3
c$$$      write(nout,200)
c$$$      write ( iotty,200 )
c$$$  100 format
c$$$     .   (1x," time check at",2x,a20,/," this segment",3x,"cpu......",
c$$$     .    1pe10.3,5x," cumulative",3x,"cpu......",1pe10.3,/,
c$$$     .    16x,"i/o......",1pe10.3,19x,"i/o......",1pe10.3,/,
c$$$     .    16x,"sys......",1pe10.3,19x,"sys......",1pe10.3)
c$$$  200 format(" ...................................................")
c$$$      return
c$$$      end
c$$$c
c......................................................................
      SUBROUTINE timer ( nout,nout2,mesage )
c......................................................................

c     DIMENSION mesage(2)
      CHARACTER*(*) mesage

      CALL cpu_time ( x1 )

      x2 = 0.0
      x3 = 0.0
      x1old = 0.0
      x2old = 0.0
      x3old = 0.0

      xdif1 = x1 - x1old
      xdif2 = x2 - x2old
      xdif3 = x3 - x3old
      x1old = x1
      x2old = x2
      x3old = x3
      
      WRITE (nout,200)
      IF ( nout2 /= 0 ) WRITE ( nout2,200 )
      WRITE (nout,100) mesage,xdif1,x1,xdif2,x2,xdif3,x3
      IF ( nout2 /= 0 )
     $     WRITE (nout2,100) mesage,xdif1,x1,xdif2,x2,xdif3,x3
      WRITE (nout,200)
      IF ( nout2 /= 0 ) WRITE ( nout2,200 )

  100 FORMAT
     .   (1x," Time Check at",2x,a20,/," This Segment",3x,"cpu sec...",
     .    1pe10.3,2x,"Cumulative",3x,"cpu-sec.......",1pe10.3,/,
     .    16x,"i/o.......",1pe10.3,19x,"i/o.......",1pe10.3,/,
     .    16x,"sys.......",1pe10.3,19x,"sys.......",1pe10.3)
  200 FORMAT(" ...................................................")

      RETURN
      END

c...........................................................
      subroutine shftpi ( xin,fin, xout,fout, n )
c...........................................................
c
      dimension xin(*), fin(*), xout(*), fout(*)
c
      nm = n-1
      n2 = nm/2
      do 10 i = 1, n
         im = i-1 + n2
         ig = mod(im,nm) + 1
         xout(ig) = xin(i)
         if ( i .gt. n2 ) xout(ig) = xin(i) - xin(n)
         fout(ig) = fin(i)
 10   continue
c
      fout(n) = fout(1)
      xout(n) = xin(n2+1)
c
      return
      end
c
c...........................................................
      subroutine shflen ( fin, fout, n, m )
c...........................................................
c
c... Shift the periodic length vector "fin" to start at 
c    fin(m). Stores the result in fout". fout(1) = 0.
c
      dimension fin(*), fout(*)
c
      do  i = 1, n - m + 1
         fout(i) = fin(m+i-1) - fin(m)
      end do
c
      do i = 1, m
         fout(n-m+i) = fout(n-m+1) + fin(i) - fin(1)
      end do
c
      return
      end
c
c.........................................................
      subroutine shftgr ( a, mth )
c........................................................
c
      dimension a(*)
c
      mthh = mth / 2
c
      do i = 1, mthh 
         za = a(i)
         a(i) = a(i+mthh)
         a(i+mthh) = za
      end do
c
      a(mth+1) = a(1)
c
      return
      end
c
c.........................................................
      subroutine shfgrb ( a, b, mth )
c........................................................
c
      dimension a(*), b(*)
c
      mthh = mth / 2
c
      do i = 1, mthh 
         za = a(i)
         b(i) = a(i+mthh)
         b(i+mthh) = za
      end do
c
      b(mth+1) = b(1)
c
      return
      end
c
c............................................................
      subroutine diff5 ( fin, h, df5 )
c...........................................................
c     
      dimension fin(*)
c     
      fm2 = fin(1)
      fm1 = fin(2)
      fp1 = fin(4)
      fp2 = fin(5)
c     
      df5 = 2.0 * ( fp1-fm1 -(fp2-fm2)/8.0 ) / (3.0*h)
c     
      return
      end
c
c     ..........................................................
      subroutine matwrtn ( a, maxj1,maxj2,l1,l2,jmax1,jmax6,jn1,jn2,
     $     label, nout1,nout2 )
c...............................................................
c     
c     writes out the elements of matrix a. jn1,jn2 roughly equally spaced
c     points.
c     
      character*(*) label
      dimension a ( maxj1, maxj2 ), jw2(101)
c
c.....      nwrt2 is the number of fields across the page
c           nwrt1 the number down the page. Not activated now!
c
      nwrt2 = 9
      nwrt1 = jmax1
c
      if ( nout1 .ne. 0 )
     $     write ( nout1, 10 ) label, jn1, jn2
      if ( nout2 .ne. 0 )
     $     write ( nout2, 10 ) label, jn1, jn2
 10   format ( /, 5x, "matrix elements of  ", a, 2i4  )
c     
c...(jmax1d,jmax2d) is the submatrix to be written. Full matrix now.
      jmax1d = jmax1
      jmax6d = jmax6
c$$$      jmax1d = jmax1 - l1 + 1
c$$$      jmax6d = jmax6 - l2 + 1

      jd1 = jmax1d / jn1
      jd2 = jmax6d / jn2
      if ( jd1 .lt. 1 ) jd1 = 1
      if ( jd2 .lt. 1 ) jd2 = 1
c$$$      njd1 = min ( jmax1,(nwrt1-1)*jd1 + 1 )
      njd2 = min ( jmax6d,(nwrt2-1)*jd2 + 1 )
      njd1 = min ( jmax1d,(nwrt1-1)*jd1 + 1 )
c
      if ( jn2 .gt. 101 ) go to 11
      j2m = 0
      do j2 = 1, jmax6d, jd2
         j2m = j2m + 1
         jw2(j2m) = j2 -1 + l2
      end do
      j2m = min(nwrt2,j2m)
      if ( nout1. ne. 0 )
     $     write ( nout1, '(10i11)' ) ( jw2(i),i = 1, j2m )
      if ( nout2. ne. 0 )
     $     write ( nout2, '(10i11)' ) ( jw2(i),i = 1, j2m )
 11   continue
c     
      do 30 j1 = 1, jmax1, jd1
         jw = j1 - 1 + l1
         if (  nout1 .ne. 0 )
     $        write ( nout1, 20 ) jw, ( a(j1,j2),j2 = 1,njd2,jd2 )
         if (  nout2 .ne. 0 )
     $        write ( nout2, 20 ) jw, ( a(j1,j2),j2 = 1,njd2,jd2 )
 20      format ( i4,10(1x,1p10e11.4,/) )
 30   continue
c     
      return
      end
c
c     ..........................................................
      subroutine matwrts ( a, maxj1,maxj2,l1,l2,jmax1,jmax6,jn1,jn2,
     $     jmin1,jmin2,label, nout1,nout2 )
c...............................................................
c     
c     writes out the elements of submatrix of  matrix a. 
c     jn1,jn2 roughly equally spaced points. 
c     Submatrix starts from (jmin1,jmin2) to (jmax1,jmax6)
c     So matrix a is offset from upper left corner. 
c     (l1,l2) is the starting labels of the full matrix. The 
c     value of the written label uses (jmin1,jmin2). (l1.l2) can be
c     optionally set it to (1,1).
c     
      character*(*) label
      dimension a ( maxj1, maxj2 ), jw2(101)
c
c.....      nwrt2 is the number of fields across the page
c           nwrt1 the number down the page.
c
      nwrt2 = 9
      nwrt1 = jmax1 - jmin1 + 1
c
      if ( nout1 .ne. 0 )
     $     write ( nout1, 10 ) label, jn1, jn2
      if ( nout2 .ne. 0 )
     $     write ( nout2, 10 ) label, jn1, jn2
 10   format ( /, 5x, "matrix elements of  ", a, 2i4  )
c     
c...(jmax1d,jmax2d) is the submatrix to be written.
      jmax1d = jmax1 - jmin1 + 1
      jmax6d = jmax6 - jmin2 + 1
c$$$      jmax1d = jmax1 - l1 + 1
c$$$      jmax6d = jmax6 - l2 + 1

      jd1 = jmax1d / jn1
      jd2 = jmax6d / jn2
      if ( jd1 .lt. 1 ) jd1 = 1
      if ( jd2 .lt. 1 ) jd2 = 1
c$$$      njd1 = min ( jmax1,(nwrt1-1)*jd1 + 1 )
      njd2 = min ( jmax6d,(nwrt2-1)*jd2 + 1 )
      njd1 = min ( jmax1d,(nwrt1-1)*jd1 + 1 )
c
c  Get vector jw2 of labels. Starts at l2.
      if ( jn2 .gt. 101 ) go to 11
      j2m = 0
      do j2 = 1, njd2, jd2
         j2m = j2m + 1
c         jw2(j2m) = j2 - 1 + l2
         jw2(j2m) = j2 - 1 + l2 + jmin2 - 1 
      end do
      j2m = min(nwrt2,j2m)
      if ( nout1. ne. 0 )
     $     write ( nout1, '(10i11)' ) ( jw2(i),i = 1, j2m )
      if ( nout2. ne. 0 )
     $     write ( nout2, '(10i11)' ) ( jw2(i),i = 1, j2m )
 11   continue
c     
c   Label jw starts at l. Submatrix starts at (jmin1.jmin2).
      do 30 j1 = 1, njd1, jd1
c         jw = j1 - 1 + l1
         jw = j1 - 1 + l1 + jmin1 - 1
         if (  nout1 .ne. 0 )
     $        write ( nout1, 20 )
     $        jw, ( a(j1-1+jmin1,j2-1+jmin2),j2 = 1,njd2,jd2 )
         if (  nout2 .ne. 0 )
     $        write ( nout2, 20 )
     $        jw, ( a(j1-1+jmin1,j2-1+jmin2),j2 = 1,njd2,jd2 )
 20      format ( i4,10(1x,1p10e11.4,/) )
 30   continue
c     
      return
      end
c
c     ..........................................................
      subroutine matwrtnbak ( a, maxj1,maxj2,l1,l2,jmax1,jmax6,jn1,jn2,
     $     label, nout1,nout2 )
c...............................................................
c     
c     writes out the elements of matrix a. jn1,jn2 roughly equally spaced
c     points.
c     
      character*(*) label
      dimension a ( maxj1, maxj2 ), jw2(101)
c
c.....      nwrt2 is the number of fields across the page
c           nwrt1 the number down the page. Not activated now!
c
      nwrt2 = 9
      nwrt1 = jmax1
c
      if ( nout1 .ne. 0 )
     $     write ( nout1, 10 ) label, jn1, jn2
      if ( nout2 .ne. 0 )
     $     write ( nout2, 10 ) label, jn1, jn2
 10   format ( /, 5x, "matrix elements of  ", a, 2i4  )
c     
      jmax1d = jmax1 - l1 + 1
      jmax6d = jmax6 - l2 + 1
      jd1 = jmax1d / jn1
      jd2 = jmax6d / jn2
      if ( jd1 .lt. 1 ) jd1 = 1
      if ( jd2 .lt. 1 ) jd2 = 1
c$$$      njd1 = min ( jmax1,(nwrt1-1)*jd1 + 1 )
      njd2 = min ( jmax6d,(nwrt2-1)*jd2 + 1 )
      njd1 = min ( jmax1d,(nwrt1-1)*jd1 + 1 )
c
      if ( jn2 .gt. 101 ) go to 11
      j2m = 0
      do j2 = l2, jmax6, jd2
         j2m = j2m + 1
         jw2(j2m) = j2 
      end do
      j2m = min(nwrt2,j2m)
      if ( nout1. ne. 0 )
     $     write ( nout1, '(10i11)' ) ( jw2(i),i = 1, j2m )
      if ( nout2. ne. 0 )
     $     write ( nout2, '(10i11)' ) ( jw2(i),i = 1, j2m )
 11   continue
c     
      do 30 j1 = l1, jmax1, jd1
         jw = j1
         if (  nout1 .ne. 0 )
     $        write ( nout1, 20 ) jw, ( a(j1,j2+l2-1),j2 = 1,njd2,jd2 )
         if (  nout2 .ne. 0 )
     $        write ( nout2, 20 ) jw, ( a(j1,j2+l2-1),j2 = 1,njd2,jd2 )
 20      format ( i4,10(1x,1p10e11.4,/) )
 30   continue
c     
      return
      end
c
c.........................................................     
      SUBROUTINE matwrtx ( a,ndi,ndj, i1,i2, j1,j2, label,
     $     nout1,nout2 )
c ...............................................................

      REAL, DIMENSION(ndi,ndj) :: a
      INTEGER, DIMENSION(ndi) :: iaxis
      INTEGER, DIMENSION(ndj) :: jaxis
      CHARACTER(*) label

      ni = i2 - i1 + 1
      nj = j2 - j1 + 1

      iaxis(1:ni) = (/ (i,i=i1, i2) /)
      jaxis(1:nj) = (/ (j,j=j1, j2) /)

      IF ( nout1 /= 0 ) then
         WRITE ( nout1, '(/,5x,"Matrix of ", a, 2i5)' ) label, i1,i2
         WRITE ( nout1, '(i6,9i11)' ) ( j, j = j1, j2 )
         DO i = i1, i2
            WRITE (nout1, '(10es11.4)' ) ( a(i,j), j = j1, j2 )
         END DO
      END IF

      IF ( nout2 /= 0 ) then
         WRITE ( nout2, '(/,5x,"Matrix of ", a, 2i5)' ) label, i1,i2
         WRITE ( nout2, '(i6,9i11)' ) ( j, j = j1, j2 )
         DO i = i1, i2
            WRITE (nout2, '(10es11.4)' ) ( a(i,j), j = j1, j2 )
         END DO
      END IF

      END SUBROUTINE matwrtx

c...............................................................
      subroutine matwrt9 ( a, maxj1,maxj2,l1,l2, label, nout1,nout2 )
c...............................................................
c     
c     writes out at most the innermost 9 elements of matrix a.
c     
      character*(*) label
      dimension a ( maxj1, maxj2 ), jw2(101)
c     
      if ( nout1 .ne. 0 )
     $     write ( nout1, 10 ) label, l1, l2
      if ( nout2 .ne. 0 )
     $     write ( nout2, 10 ) label, l1, l2
 10   format ( /, 5x, "matrix elements of  ", a, 2i4  )
c     
      ln4 = max ( l1, -4 )
      lx4 = min ( l2, +4 )
      jnx = lx4 - ln4 + 1
c     
      do j2 = 1, jnx
         jw2(j2) = j2 - 1 + ln4
      end do
      if ( nout1. ne. 0 )
     $     write ( nout1, '(10i11)' ) ( jw2(i),i = 1, jnx )
      if ( nout2. ne. 0 )
     $     write ( nout2, '(10i11)' ) ( jw2(i),i = 1, jnx )
 11   continue
c     
      js = max ( -4-l1, 0 ) + 1
c     
      lw = 0
      do 30 j1 = js, js+jnx-1
c$$$  jw = ln + j1 + 1 
         lw = lw + 1
         if (  nout1 .ne. 0 )
     $        write ( nout1, 20 ) jw2(lw),(a(j1,j2),j2 = js, js+jnx-1)
         if (  nout2 .ne. 0 )
     $        write ( nout2, 20 ) jw2(Lw),(a(j1,j2),j2 = js, js+jnx-1)
 20      format ( i4,10(1x,1p10e11.4,/) )
 30   continue
c     
      return
      end
c     ..........................................................
      subroutine mawrtnn ( a, maxj1,maxj2,maxn, l1,l2,n1,
     $     jmax1,jmax6,nmax, jn1,jn2, label, nout1,nout2 )
c...............................................................
c     
c     writes out the elements of matrix a(l1,l2,n) for each n.
c      jn1,jn2 roughly equally spaced in l1,l2.
c     
      character*(*) label
      dimension a ( maxj1, maxj2, maxn )
c     
      do 50 jn = 1, nmax
         nw = n1 + jn - 1
c     
         if ( nout1 .ne. 0 )
     $        write ( nout1, 10 ) label, jn1, jn2, nw
         if ( nout2 .ne. 0 )
     $        write ( nout2, 10 ) label, jn1, jn2, nw
 10      format ( /, 5x, "matrix elements of  ", a, 2i4, " n = ", i4  )
c     
         jd1 = jmax1 / jn1
         jd2 = jmax6 / jn2
         if ( jd1 .lt. 1 ) jd1 = 1
         if ( jd2 .lt. 1 ) jd2 = 1
c     
         do 30 j1 = 1, jmax1,jd1
            jw = j1 -1 + l1
            if (  nout1 .ne. 0 ) write ( nout1, 20 ) jw,
     $           ( a(j1,j2,jn),j2 = 1, jmax6,jd2 )
            if (  nout2 .ne. 0 ) write ( nout2, 20 ) jw,
     $           ( a(j1,j2,jn),j2 = 1, jmax6,jd2 )
 20         format ( i4,10(1x,1p10e11.4,/) )
 30      continue
 50   continue
c     
      return
      end
c     
c.................................................................
      subroutine mtrans ( a, nd, n )
c.................................................................
c
c... Transpose a square  matrix, a, into itself. a is overwritten.
c
      dimension a(nd,nd)
c
      do 200 i = 1, n
         do 100 j = 1, i
            z = a(i,j)
            a(i,j) = a(j,i)
            a(j,i) = z
 100     continue
 200  continue
c
      return
      end
c
c.................................................................
      subroutine mtransr ( a, nd1,nd2, n1,n2, b )
c.................................................................
c
c... Transpose a rectangular matrix, a, into b. a is not overwrittten.
c
      dimension a(nd1,nd2), b(nd2,nd1)
c
      do 200 i = 1, n1
         do 100 j = 1, n2
            b(j,i) = a(i,j)
 100     continue
 200  continue
c
      return
      end
c
c..............................................................
      subroutine scalmul ( a, n1,n2, f )
c...............................................................
c
      dimension a(n1,n2)
c
      do i = 1, n1
         do j = 1, n2
            a(i,j) = f * a(i,j)
         end do
      end do
c
      return
      end
c
c..............................................................
      subroutine scalmul0 ( a, n1,n2,ln,lx, f )
c...............................................................
c
      dimension a(n1,n2)
c
c....Multiply all but l=0 rows of a by f.
c
      do j = 1, n2
         if ( j .eq. (-ln + 1)) go to 100
         do i = 1, n1
            a(i,j) = f * a(i,j)
         end do
 100     continue
      end do
c
      return
      end
c
c..............................................................
      subroutine vecmula ( a, n1,n2,ln,lx, fv )
c...............................................................
c
      dimension a(n1,n2), fv(*)
c
c.... Multiply rows of a by the vector fv
c
      do j = 1, n2
         do i = 1, n1
            a(i,j) = fv(j) * a(i,j)
         end do
 100     continue
      end do
c
      return
      end
c
c............................................................
      subroutine shiftmat ( a, nd1,nd2, l1,l2, ldl,ldu )
c............................................................
c
c.... Shifts the elements of a upwards and leftwards by ldl,ldu.
c     Need only l1,l2 elements.
c
      dimension a(nd1,nd2)
c
      do i = 1, l2
         do j = 1, l1
            a(j,i) = a(j+ldu,i+ldl)
         end do
      end do
c
      return
      end
c
c............................................................
      SUBROUTINE reversevecc ( a,b, mth1 )
c............................................................

c... Changes the direction of the indexing of an array a, outputs b

      COMPLEX, DIMENSION(*) :: a, b

      b(1:mth1) = (/ (a(mth1+1-j),j=1,mth1) /)

      RETURN
      END SUBROUTINE reversevecc

c............................................................
      SUBROUTINE reversematc ( a, nd1,nd2, mth )
c............................................................

c... Changes the direction of the indexing of both rows and columns
c    in a complex matrix. Matrix overwritten

      COMPLEX, DIMENSION(nd1,nd2) :: a
      
      mth1 = mth + 1
      mth2 = mth1 + 1
c.. First the rows
      DO i = 1, mth1
         a(i,1:mth1) = (/ (a(i,mth1+1-j),j=1,mth1) /)
      END DO

c.. Now the columns
      DO j = 1, mth1
         a(1:mth1,j) = (/ (a(mth1+1-i,j),i=1,mth1) /)
      END DO

c.. Put in the data at the mth2 points
      DO i = 1, mth1
         a(i,mth2) = a(i,2)
      END DO

      DO j = 1, mth2
         a(mth2,j) = a(2,j)
      END DO

      RETURN
      END SUBROUTINE reversematc

c......................................................
      SUBROUTINE transc ( a, b, nd1,nd2, mth,mthout )
c......................................................

c... Interpolate a complex matrix from a mth to a mthout grid.
c  Uses TRANS an interpolater on a real array.

      COMPLEX, DIMENSION(nd1,nd2) :: a, b
      REAL, DIMENSION(:), ALLOCATABLE :: zgri, zgro, zgrop

      mth2 = mth + 2
      mthout2 = mthout + 2
      ALLOCATE ( zgri(mth2+3), zgro(mth2+3), zgrop(mth2+3) )

c... Interpolate columns to mthout points

      DO i = 1, mth2
         zgri(1:mth2) = REAL( a(1:mth2, i) )
         CALL trans ( zgri,mth, zgro,mthout )
         zgri(1:mth2) = AIMAG( a(1:mth2, i) )
         CALL trans ( zgri,mth, zgrop,mthout )
         b(1:mthout2,i) = CMPLX( zgro(1:mthout2),zgrop(1:mthout2) )
      END DO

c... Now the rows

      DO i = 1, mthout2
         zgri(1:mth2) = REAL( b(i, 1:mth2) )
         CALL trans ( zgri,mth, zgro,mthout )
         zgri(1:mth2) = AIMAG( b(i, 1:mth2) )
         CALL trans ( zgri,mth, zgrop,mthout )
         b(i,1:mthout2) = CMPLX( zgro(1:mthout2),zgrop(1:mthout2) )
      END DO

      DEALLOCATE ( zgri, zgro, zgrop )

      RETURN
      END SUBROUTINE transc
