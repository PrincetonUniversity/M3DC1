c.....................................................................
      SUBROUTINE drawc6(xx,arr,xx2,arr2,xx3,arr3,xx4,arr4,xx5,arr5,
     $     xx6,arr6,nl,nr,nl3,nr3,nl4,nr4,nl5,nr5,nl6,nr6,
     $     yla,xla,xmx,xma,zma,
     $     plrad,if,im,idot, jobid, ishape,a,b,aw,bw,cw,dw,tw,
     $     abulg, bbulg, tbulg, wcentr, pcircm, wcircm,
     $     pelong, pdness, wmaj,wlrad, welong, wdness, chgt,chgta,chgtb)
c.......................................................................
c
c     To plot a function arr of xx, from xx(nl) to xx(nr).
c     Called by wwall in module vacuum_vac.f
c
      dimension arr(*), xx(*), arr2(*), xx2(*), arr3(*), xx3(*),
     $     arr4(*), xx4(*), arr5(*), xx5(*),arr6(*), xx6(*),
     $     chgt(*), chgta(*), chgtb(*) 
c
      character*(*) yla, xla, jobid
      character*(80) stringv(5), string1
c
c     find minimum and maximum of arr and xx...

      CALL bounds ( arr,xx, nl,nr, amin,amax, bmin,bmax )
      CALL bounds ( arr2,xx2, nl,nr, amin2,amax2, bmin2,bmax2 )
      CALL bounds ( arr3,xx3, nl3,nr3, amin3,amax3, bmin3,bmax3 )
      CALL bounds ( arr4,xx4, nl4,nr4, amin4,amax4, bmin4,bmax4 )
      CALL bounds ( arr5,xx5, nl5,nr5, amin5,amax5, bmin5,bmax5 )
      CALL bounds ( arr6,xx6, nl6,nr6, amin6,amax6, bmin6,bmax6 )

c$$$      amin = 1.0e20
c$$$      amax = -1.0e20
c$$$      bmin = 1.0e20
c$$$      bmax = -1.0e20
c$$$      amin2 =  1.0e20
c$$$      amax2 = -1.0e20
c$$$      bmin2 =  1.0e20
c$$$      bmax2 = -1.0e20
c$$$      amin3 =  1.0e20
c$$$      amax3 = -1.0e20
c$$$      bmin3 =  1.0e20
c$$$      bmax3 = -1.0e20
c$$$      amin4 =  1.0e20
c$$$      amax4 = -1.0e20
c$$$      bmin4 =  1.0e20
c$$$      bmax4 = -1.0e20
c$$$      amin5 =  1.0e20
c$$$      amax5 = -1.0e20
c$$$      bmin5 =  1.0e20
c$$$      bmax5 = -1.0e20
c$$$
c$$$c
c$$$      do  nn= nl, nr
c$$$         if ( arr(nn) .lt. amin )  amin = arr(nn)
c$$$         if ( arr(nn) .gt. amax )  amax = arr(nn)
c$$$         if ( xx(nn) .lt. bmin ) bmin = xx(nn)
c$$$         if ( xx(nn) .gt. bmax ) bmax = xx(nn)
c$$$c
c$$$         if ( arr2(nn) .lt. amin2 ) amin2 = arr2(nn)
c$$$         if ( arr2(nn) .gt. amax2 ) amax2 = arr2(nn)
c$$$         if ( xx2(nn) .lt. bmin2 ) bmin2 = xx2(nn)
c$$$         if ( xx2(nn) .gt. bmax2 ) bmax2 = xx2(nn)
c$$$      end do
c$$$c
c$$$      do nn = nl3, nr3
c$$$         if ( arr3(nn) .lt. amin3 ) amin3 = arr3(nn)
c$$$         if ( arr3(nn) .gt. amax3 ) amax3 = arr3(nn)
c$$$         if ( xx3(nn) .lt. bmin3 ) bmin3 = xx3(nn)
c$$$         if ( xx3(nn) .gt. bmax3 ) bmax3 = xx3(nn)
c$$$      end do
c$$$c
c$$$      do nn = nl4, nr4
c$$$         if ( arr4(nn) .lt. amin4 ) amin4 = arr4(nn)
c$$$         if ( arr4(nn) .gt. amax4 ) amax4 = arr4(nn)
c$$$         if ( xx4(nn) .lt. bmin4 ) bmin4 = xx4(nn)
c$$$         if ( xx4(nn) .gt. bmax4 ) bmax4 = xx4(nn)
c$$$      end do
c$$$
c$$$      do nn = nl5, nr5
c$$$         if ( arr5(nn) .lt. amin5 ) amin5 = arr5(nn)
c$$$         if ( arr5(nn) .gt. amax5 ) amax5 = arr5(nn)
c$$$         if ( xx5(nn) .lt. bmin5 ) bmin5 = xx5(nn)
c$$$         if ( xx5(nn) .gt. bmax5 ) bmax5 = xx5(nn)
c$$$      end do

c     write ( 23,101 ) amin,amax, bmin,bmax, amin2,amax2, bmin2,bmax2
 101  format ( 1x,"amin,amax, bmin,bmax, amin2,amax2, bmin2,bmax2 = ",
     .     /, 1p8e12.5 )
c
      amin = amin1 ( amin,amin2 )
      bmin = amin1 ( bmin,bmin2 )
      amax = amax1 ( amax,amax2 )
      bmax = amax1 ( bmax,bmax2 )
c
      amin = amin1 ( amin,amin3 )
      bmin = amin1 ( bmin,bmin3 )
      amax = amax1 ( amax,amax3 )
      bmax = amax1 ( bmax,bmax3 )
c
      amin = amin1 ( amin,amin4 )
      bmin = amin1 ( bmin,bmin4 )
      amax = amax1 ( amax,amax4 )
      bmax = amax1 ( bmax,bmax4 )

      amin = amin1 ( amin,amin5 )
      bmin = amin1 ( bmin,bmin5 )
      amax = amax1 ( amax,amax5 )
      bmax = amax1 ( bmax,bmax5 )

      amin = amin1 ( amin,amin6 )
      bmin = amin1 ( bmin,bmin6 )
      amax = amax1 ( amax,amax6 )
      bmax = amax1 ( bmax,bmax6 )

c     write ( 23,101 ) amin,amax, bmin,bmax, amin2,amax2, bmin2,bmax2
c
c.....fix it so that grid is undistorted...
c
      dela = amax - amin
      delb = bmax - bmin
c
      if ( dela .gt. delb ) go to 4
c
      b2 = bmax + 0.1*delb
      b1 = bmin - 0.1*delb
      a2 = ( amax + amin + b2 - b1 ) / 2.0
      a1 = ( amax + amin - b2 + b1 ) / 2.0
c
      go to 6
c
    4 continue
c
      a2 = amax + 0.1*dela
      a1 = amin - 0.1*dela
      b2 = ( bmax + bmin + a2 - a1 ) / 2.0
      b1 = ( bmax + bmin - a2 + a1 ) / 2.0
c
    6 continue
c
c     write (23,105 ) a1,a2, b1,b2
 105  format ( 1x, "a1, a2, b1, b2 = ",/, 1p4e12.5 )
c
      call maps ( b1,b2, a1,a2, .114,.828,.285,1.0 )
c
c     write (23,105 ) a1,a2, b1,b2
c
      xl = xx(nl)
      xr = xx(nr)
      dx = ( xr-xl ) / 50.0
      call pointc ( 'x', xmx, 0.0, 1, -1, -1, 0., 0. )
      if ( a .lt. -100.0 )
     $     call pointc ( 'x', -xmx, 0.0, 1, -1, -1, 0., 0. )
      call pointc ( '0', xma, zma, 1, -1, -1, 0., 0. )
      if ( a .lt. -100.0 )
     $     call pointc ( '0', -xma, zma, 1, -1, -1, 0., 0. )
      call pointc ( '*', wcentr, 0.0, 1, -1, -1, 0., 0. )
      call pointc ( 'X', chgt(1),  chgt(2), 1, -1, -1, 0., 0. )
      call pointc ( 'X', chgt(3),  chgt(4), 1, -1, -1, 0., 0. )
      call pointc ( 'X', chgta(1),  chgta(2), 1, -1, -1, 0., 0. )
      call pointc ( 'X', chgta(3),  chgta(4), 1, -1, -1, 0., 0. )
      call pointc ( 'X', chgtb(1),  chgtb(2), 1, -1, -1, 0., 0. )
      call pointc ( 'X', chgtb(3),  chgtb(4), 1, -1, -1, 0., 0. )
      call pointc ( '+', xx6(nl6),  arr6(nl6), 1, -1, -1, 0., 0. )
      call pointc ( '+', xx6(nr6),  arr6(nr6), 1, -1, -1, 0., 0. )

c
      if ( idot .eq. 0 ) then
c
         call point ( xx(nl),arr(nl) )
c
         do i = nl, nr
            x = xx(i)
            y = arr(i)
            call vector ( x, y )
         end do
c
         call point ( xx2(nl),arr2(nl) )
c
         do i = nl, nr
            x = xx2(i)
            y = arr2(i)
            call vector ( x,y )
         end do
c
         if ( a .lt. -100.0 ) then
            call point ( -xx2(nl),arr2(nl) )
c
            do i = nl, nr
               x = -xx2(i)
               y = arr2(i)
               call vector ( x, y )
            end do
c
         end if
c
      nlr = nr3 - nl3 + 1
      call points ( xx3, arr3, nlr, -1,-1,0.,0. )

      nlr = nr4 - nl4 + 1
      call points ( xx4, arr4, nlr, -1,-1,0.,0. )
c
      nlr = nr5 - nl5 + 1
      call points ( xx5, arr5, nlr, -1,-1,0.,0. )
c
      nlr = nr6 - nl6 + 1
      call points ( xx6, arr6, nlr, -1,-1,0.,0. )

c
c$$$
c$$$         call point ( xx3(nl3),arr3(nl3) )
c$$$c
c$$$         do i = nl3, nr3
c$$$            x = xx3(i)
c$$$            y = arr3(i)
c$$$            call vector ( x,y )
c$$$         end do
c
         go to 60
c
      end if
c
 15   continue
c
      go to 50
c
c$$$      call setpch ( 0, 0, -1,-1 )
c$$$      do i = nl, nr
c$$$         call pointc ( 1hx, xx(i), arr(i), 1, -1, -1, 0., 0. )
c$$$      end do
c$$$c
c$$$      do i = nl, nr
c$$$         call pointc ( 1h0, xx2(i), arr2(i), 1, -1, -1, 0., 0. )
c$$$      end do
c$$$c
c$$$      do i = nl3, nr3
c$$$         call pointc ( 1h*, xx3(i), arr3(i), 1, -1, -1, 0., 0. )
c$$$      end do
c
 50   continue
c
      nlr = nr - nl + 1
      call points ( xx, arr, nlr, -1,-1,0.,0. )
      call points ( xx2, arr2, nlr, -1,-1,0.,0. )
c$$$      nlr = nr3 - nl3 + 1
c$$$      call points ( xx3, arr3, nlr, -1,-1,0.,0. )
c
 60   continue
c
      call dders(1)
c     ab = a1 - 0.1 * ( a2-a1 )
c     xm = 0.5 * ( b1+b2 )
c     call setlch ( xm,ab,0,2,0,-1 )
c     call gtext ( xla, -1, 0 )
      xm = b1 - 0.1 * ( b2-b1 )
c     ab = 0.5 * ( a1 + a2 )
      ab =  a1 + 0.4 * ( a2-a1 )
      string1 = xla//'-'//yla//' plane'
      call setlch ( xm,ab,0,2,1,-1 )
      call gtext ( string1, -1, 0 )
      xm = b1 + 0.05 * ( b2-b1 )
      ab = a2 - 0.05 * ( a2 - a1 )
      call setlch ( xm,ab, 0,1,0,-1 )
      call gtext ( jobid, -1, 0 )
      xm = b1 - (b2-b1) / 10.0
      ab = a1 - (a2-a1) / 8.0
c
      call setlch ( xm,ab, 0,1,0,-1 )
c
      write( 6,'("Circumference of Plasma, pl = ",1pe11.3)')
     $     pcircm
      write(23,'("Circumference of Plasma, pl = ",1pe11.3)')
     $     pcircm
c
      write ( stringv, 3001 ) plrad, xmx, xma, pcircm
 3001 format ( 1x, "plrad= ", f7.3, 1x, "xmaj-[x]= ",f7.3, 1x,
     $     "xma-[0]= ", f7.3, " Pcircum= ", f7.3 )
      call wrtstr ( stringv, 1 )
c
      write ( stringv, 3002) ishape, wcentr, a, b, aw, bw, cw, dw, tw,
     $     abulg, bbulg, tbulg
 3002 format ( 1x, "ishape= ", i3, " wcentr-[*]=",f6.3,
     $     " a= ",f8.3, " b= ",f7.3,/,
     $     1x, "aw= ",f8.3," bw= ",f8.3,
     $     " cw= ",f8.3," dw= ",f8.3," tw= ",f8.3,/,
     $     1x, "abulg= "f8.3, " bbulg= ",f8.3, " tbulg= ", f8.3 )
      call wrtstr ( stringv, 3 )
c
      write ( stringv, 3003) pelong, pdness, wcircm,
     $     wmaj,wlrad, welong, wdness,
     $     chgt(1), chgt(2), chgt(3), chgt(4)
 3003 format ( 1x, "pelong= ",f8.3, " pdness= ",f8.3,
     $     " wcircm= ",f8.3,/,
     $     " wmaj= ",f8.3, " wlrad = ",f8.3, " welong= ",f8.3,
     $     " wdness= ",f8.3, /,
     $     " Coil Tips = ", 2f8.3, " and ", 2f8.3 )
      call wrtstr ( stringv, 3 )
c
      if ( if .eq. 1 ) call framep(jobid, ff)
c
      RETURN
      END
c.....................................................................
      subroutine drawc5(xx,arr,xx2,arr2,xx3,arr3,xx4,arr4,xx5,arr5,
     $     nl,nr,nl3,nr3,nl4,nr4,nl5,nr5,
     $     yla,xla,xmx,xma,zma,
     $     plrad,if,im,idot, jobid, ishape,a,b,aw,bw,cw,dw,tw,
     $     abulg, bbulg, tbulg, wcentr, pcircm, wcircm,
     $     pelong, pdness, wmaj,wlrad, welong, wdness, chgt,chgta,chgtb)
c.......................................................................
c
c     To plot a function arr of xx, from xx(nl) to xx(nr).
c     Called by wwall in module vacuum_vac.f
c
      dimension arr(*), xx(*), arr2(*), xx2(*), arr3(*), xx3(*),
     $     arr4(*), xx4(*), arr5(*), xx5(*),
     $     chgt(*), chgta(*), chgtb(*) 
c
      character*(*) yla, xla, jobid
      character*(80) stringv(5), string1
c
c     find minimum and maximum of arr and xx...

      CALL bounds ( arr,xx, nl,nr, amin,amax, bmin,bmax )
      CALL bounds ( arr2,xx2, nl,nr, amin2,amax2, bmin2,bmax2 )
      CALL bounds ( arr3,xx3, nl3,nr3, amin3,amax3, bmin3,bmax3 )
      CALL bounds ( arr4,xx4, nl4,nr4, amin4,amax4, bmin4,bmax4 )
      CALL bounds ( arr5,xx5, nl5,nr5, amin5,amax5, bmin5,bmax5 )

c$$$      amin = 1.0e20
c$$$      amax = -1.0e20
c$$$      bmin = 1.0e20
c$$$      bmax = -1.0e20
c$$$      amin2 =  1.0e20
c$$$      amax2 = -1.0e20
c$$$      bmin2 =  1.0e20
c$$$      bmax2 = -1.0e20
c$$$      amin3 =  1.0e20
c$$$      amax3 = -1.0e20
c$$$      bmin3 =  1.0e20
c$$$      bmax3 = -1.0e20
c$$$      amin4 =  1.0e20
c$$$      amax4 = -1.0e20
c$$$      bmin4 =  1.0e20
c$$$      bmax4 = -1.0e20
c$$$      amin5 =  1.0e20
c$$$      amax5 = -1.0e20
c$$$      bmin5 =  1.0e20
c$$$      bmax5 = -1.0e20
c$$$
c$$$c
c$$$      do  nn= nl, nr
c$$$         if ( arr(nn) .lt. amin )  amin = arr(nn)
c$$$         if ( arr(nn) .gt. amax )  amax = arr(nn)
c$$$         if ( xx(nn) .lt. bmin ) bmin = xx(nn)
c$$$         if ( xx(nn) .gt. bmax ) bmax = xx(nn)
c$$$c
c$$$         if ( arr2(nn) .lt. amin2 ) amin2 = arr2(nn)
c$$$         if ( arr2(nn) .gt. amax2 ) amax2 = arr2(nn)
c$$$         if ( xx2(nn) .lt. bmin2 ) bmin2 = xx2(nn)
c$$$         if ( xx2(nn) .gt. bmax2 ) bmax2 = xx2(nn)
c$$$      end do
c$$$c
c$$$      do nn = nl3, nr3
c$$$         if ( arr3(nn) .lt. amin3 ) amin3 = arr3(nn)
c$$$         if ( arr3(nn) .gt. amax3 ) amax3 = arr3(nn)
c$$$         if ( xx3(nn) .lt. bmin3 ) bmin3 = xx3(nn)
c$$$         if ( xx3(nn) .gt. bmax3 ) bmax3 = xx3(nn)
c$$$      end do
c$$$c
c$$$      do nn = nl4, nr4
c$$$         if ( arr4(nn) .lt. amin4 ) amin4 = arr4(nn)
c$$$         if ( arr4(nn) .gt. amax4 ) amax4 = arr4(nn)
c$$$         if ( xx4(nn) .lt. bmin4 ) bmin4 = xx4(nn)
c$$$         if ( xx4(nn) .gt. bmax4 ) bmax4 = xx4(nn)
c$$$      end do
c$$$
c$$$      do nn = nl5, nr5
c$$$         if ( arr5(nn) .lt. amin5 ) amin5 = arr5(nn)
c$$$         if ( arr5(nn) .gt. amax5 ) amax5 = arr5(nn)
c$$$         if ( xx5(nn) .lt. bmin5 ) bmin5 = xx5(nn)
c$$$         if ( xx5(nn) .gt. bmax5 ) bmax5 = xx5(nn)
c$$$      end do

c     write ( 23,101 ) amin,amax, bmin,bmax, amin2,amax2, bmin2,bmax2
 101  format ( 1x,"amin,amax, bmin,bmax, amin2,amax2, bmin2,bmax2 = ",
     .     /, 1p8e12.5 )
c
      amin = amin1 ( amin,amin2 )
      bmin = amin1 ( bmin,bmin2 )
      amax = amax1 ( amax,amax2 )
      bmax = amax1 ( bmax,bmax2 )
c
      amin = amin1 ( amin,amin3 )
      bmin = amin1 ( bmin,bmin3 )
      amax = amax1 ( amax,amax3 )
      bmax = amax1 ( bmax,bmax3 )
c
      amin = amin1 ( amin,amin4 )
      bmin = amin1 ( bmin,bmin4 )
      amax = amax1 ( amax,amax4 )
      bmax = amax1 ( bmax,bmax4 )

      amin = amin1 ( amin,amin5 )
      bmin = amin1 ( bmin,bmin5 )
      amax = amax1 ( amax,amax5 )
      bmax = amax1 ( bmax,bmax5 )

c     write ( 23,101 ) amin,amax, bmin,bmax, amin2,amax2, bmin2,bmax2
c
c.....fix it so that grid is undistorted...
c
      dela = amax - amin
      delb = bmax - bmin
c
      if ( dela .gt. delb ) go to 4
c
      b2 = bmax + 0.1*delb
      b1 = bmin - 0.1*delb
      a2 = ( amax + amin + b2 - b1 ) / 2.0
      a1 = ( amax + amin - b2 + b1 ) / 2.0
c
      go to 6
c
    4 continue
c
      a2 = amax + 0.1*dela
      a1 = amin - 0.1*dela
      b2 = ( bmax + bmin + a2 - a1 ) / 2.0
      b1 = ( bmax + bmin - a2 + a1 ) / 2.0
c
    6 continue
c
c     write (23,105 ) a1,a2, b1,b2
 105  format ( 1x, "a1, a2, b1, b2 = ",/, 1p4e12.5 )
c
      call maps ( b1,b2, a1,a2, .114,.828,.285,1.0 )
c
c     write (23,105 ) a1,a2, b1,b2
c
      xl = xx(nl)
      xr = xx(nr)
      dx = ( xr-xl ) / 50.0
      call pointc ( 'x', xmx, 0.0, 1, -1, -1, 0., 0. )
      if ( a .lt. -100.0 )
     $     call pointc ( 'x', -xmx, 0.0, 1, -1, -1, 0., 0. )
      call pointc ( '0', xma, zma, 1, -1, -1, 0., 0. )
      if ( a .lt. -100.0 )
     $     call pointc ( '0', -xma, zma, 1, -1, -1, 0., 0. )
      call pointc ( '*', wcentr, 0.0, 1, -1, -1, 0., 0. )
      call pointc ( 'X', chgt(1),  chgt(2), 1, -1, -1, 0., 0. )
      call pointc ( 'X', chgt(3),  chgt(4), 1, -1, -1, 0., 0. )
      call pointc ( 'X', chgta(1),  chgta(2), 1, -1, -1, 0., 0. )
      call pointc ( 'X', chgta(3),  chgta(4), 1, -1, -1, 0., 0. )
      call pointc ( 'X', chgtb(1),  chgtb(2), 1, -1, -1, 0., 0. )
      call pointc ( 'X', chgtb(3),  chgtb(4), 1, -1, -1, 0., 0. )
      call pointc ( '+', xx5(nl5),  arr5(nl5), 1, -1, -1, 0., 0. )
      call pointc ( '+', xx5(nr5),  arr5(nr5), 1, -1, -1, 0., 0. )

c
      if ( idot .eq. 0 ) then
c
         call point ( xx(nl),arr(nl) )
c
         do i = nl, nr
            x = xx(i)
            y = arr(i)
            call vector ( x, y )
         end do
c
         call point ( xx2(nl),arr2(nl) )
c
         do i = nl, nr
            x = xx2(i)
            y = arr2(i)
            call vector ( x,y )
         end do
c
         if ( a .lt. -100.0 ) then
            call point ( -xx2(nl),arr2(nl) )
c
            do i = nl, nr
               x = -xx2(i)
               y = arr2(i)
               call vector ( x, y )
            end do
c
         end if
c
      nlr = nr3 - nl3 + 1
      call points ( xx3, arr3, nlr, -1,-1,0.,0. )

      nlr = nr4 - nl4 + 1
      call points ( xx4, arr4, nlr, -1,-1,0.,0. )
c
      nlr = nr5 - nl5 + 1
      call points ( xx5, arr5, nlr, -1,-1,0.,0. )
c
c$$$
c$$$         call point ( xx3(nl3),arr3(nl3) )
c$$$c
c$$$         do i = nl3, nr3
c$$$            x = xx3(i)
c$$$            y = arr3(i)
c$$$            call vector ( x,y )
c$$$         end do
c
         go to 60
c
      end if
c
 15   continue
c
      go to 50
c
c$$$      call setpch ( 0, 0, -1,-1 )
c$$$      do i = nl, nr
c$$$         call pointc ( 1hx, xx(i), arr(i), 1, -1, -1, 0., 0. )
c$$$      end do
c$$$c
c$$$      do i = nl, nr
c$$$         call pointc ( 1h0, xx2(i), arr2(i), 1, -1, -1, 0., 0. )
c$$$      end do
c$$$c
c$$$      do i = nl3, nr3
c$$$         call pointc ( 1h*, xx3(i), arr3(i), 1, -1, -1, 0., 0. )
c$$$      end do
c
 50   continue
c
      nlr = nr - nl + 1
      call points ( xx, arr, nlr, -1,-1,0.,0. )
      call points ( xx2, arr2, nlr, -1,-1,0.,0. )
      nlr = nr3 - nl3 + 1
      call points ( xx3, arr3, nlr, -1,-1,0.,0. )
c
 60   continue
c
      call dders(1)
c     ab = a1 - 0.1 * ( a2-a1 )
c     xm = 0.5 * ( b1+b2 )
c     call setlch ( xm,ab,0,2,0,-1 )
c     call gtext ( xla, -1, 0 )
      xm = b1 - 0.1 * ( b2-b1 )
c     ab = 0.5 * ( a1 + a2 )
      ab =  a1 + 0.4 * ( a2-a1 )
      string1 = xla//'-'//yla//' plane'
      call setlch ( xm,ab,0,2,1,-1 )
      call gtext ( string1, -1, 0 )
      xm = b1 + 0.05 * ( b2-b1 )
      ab = a2 - 0.05 * ( a2 - a1 )
      call setlch ( xm,ab, 0,1,0,-1 )
      call gtext ( jobid, -1, 0 )
      xm = b1 - (b2-b1) / 10.0
      ab = a1 - (a2-a1) / 8.0
c
      call setlch ( xm,ab, 0,1,0,-1 )
c
      write( 6,'("Circumference of Plasma, pl = ",1pe11.3)')
     $     pcircm
      write(23,'("Circumference of Plasma, pl = ",1pe11.3)')
     $     pcircm
c
      write ( stringv, 3001 ) plrad, xmx, xma, pcircm
 3001 format ( 1x, "plrad= ", f7.3, 1x, "xmaj-[x]= ",f7.3, 1x,
     $     "xma-[0]= ", f7.3, " Pcircum= ", f7.3 )
      call wrtstr ( stringv, 1 )
c
      write ( stringv, 3002) ishape, wcentr, a, b, aw, bw, cw, dw, tw,
     $     abulg, bbulg, tbulg
 3002 format ( 1x, "ishape= ", i3, " wcentr-[*]=",f6.3,
     $     " a= ",f8.3, " b= ",f7.3,/,
     $     1x, "aw= ",f8.3," bw= ",f8.3,
     $     " cw= ",f8.3," dw= ",f8.3," tw= ",f8.3,/,
     $     1x, "abulg= "f8.3, " bbulg= ",f8.3, " tbulg= ", f8.3 )
      call wrtstr ( stringv, 3 )
c
      write ( stringv, 3003) pelong, pdness, wcircm,
     $     wmaj,wlrad, welong, wdness,
     $     chgt(1), chgt(2), chgt(3), chgt(4)
 3003 format ( 1x, "pelong= ",f8.3, " pdness= ",f8.3,
     $     " wcircm= ",f8.3,/,
     $     " wmaj= ",f8.3, " wlrad = ",f8.3, " welong= ",f8.3,
     $     " wdness= ",f8.3, /,
     $     " Coil Tips = ", 2f8.3, " and ", 2f8.3 )
      call wrtstr ( stringv, 3 )
c
      if ( if .eq. 1 ) call framep(jobid, ff)
c
      return
      end
c................................................................
      subroutine drawc4(xx,arr,xx2,arr2,xx3,arr3,xx4,arr4,
     $     nl,nr,nl3,nr3,nl4,nr4,
     $     yla,xla,xmx,xma,zma,
     $     plrad,if,im,idot, jobid, ishape,a,b,aw,bw,cw,dw,tw,
     $     abulg, bbulg, tbulg, wcentr, pcircm, wcircm,
     $     pelong, pdness, wmaj,wlrad, welong, wdness, chgt )
c.......................................................................
c
c     To plot a function arr of xx, from xx(nl) to xx(nr).
c     Called by wwall in module vacuum_vac.f
c
      dimension arr(*), xx(*), arr2(*), xx2(*), arr3(*), xx3(*), chgt(*)
      dimension arr4(*),xx4(*)
c
      character*(*) yla, xla, jobid
      character*(80) stringv(5), string1
c
c     find minimum and maximum of arr and xx...
      amin = 1.0e20
      amax = -1.0e20
      bmin = 1.0e20
      bmax = -1.0e20
      amin2 =  1.0e20
      amax2 = -1.0e20
      bmin2 =  1.0e20
      bmax2 = -1.0e20
      amin3 =  1.0e20
      amax3 = -1.0e20
      bmin3 =  1.0e20
      bmax3 = -1.0e20
      amin4 =  1.0e20
      amax4 = -1.0e20
      bmin4 =  1.0e20
      bmax4 = -1.0e20
c
      do  nn= nl, nr
         if ( arr(nn) .lt. amin )  amin = arr(nn)
         if ( arr(nn) .gt. amax )  amax = arr(nn)
         if ( xx(nn) .lt. bmin ) bmin = xx(nn)
         if ( xx(nn) .gt. bmax ) bmax = xx(nn)
c
         if ( arr2(nn) .lt. amin2 ) amin2 = arr2(nn)
         if ( arr2(nn) .gt. amax2 ) amax2 = arr2(nn)
         if ( xx2(nn) .lt. bmin2 ) bmin2 = xx2(nn)
         if ( xx2(nn) .gt. bmax2 ) bmax2 = xx2(nn)
      end do
c
      do nn = nl3, nr3
         if ( arr3(nn) .lt. amin3 ) amin3 = arr3(nn)
         if ( arr3(nn) .gt. amax3 ) amax3 = arr3(nn)
         if ( xx3(nn) .lt. bmin3 ) bmin3 = xx3(nn)
         if ( xx3(nn) .gt. bmax3 ) bmax3 = xx3(nn)
      end do
c
      do nn = nl4, nr4
         if ( arr4(nn) .lt. amin4 ) amin4 = arr4(nn)
         if ( arr4(nn) .gt. amax4 ) amax4 = arr4(nn)
         if ( xx4(nn) .lt. bmin4 ) bmin4 = xx4(nn)
         if ( xx4(nn) .gt. bmax4 ) bmax4 = xx4(nn)
      end do
c     write ( 23,101 ) amin,amax, bmin,bmax, amin2,amax2, bmin2,bmax2
 101  format ( 1x,"amin,amax, bmin,bmax, amin2,amax2, bmin2,bmax2 = ",
     .     /, 1p8e12.5 )
c
      amin = amin1 ( amin,amin2 )
      bmin = amin1 ( bmin,bmin2 )
      amax = amax1 ( amax,amax2 )
      bmax = amax1 ( bmax,bmax2 )
c
      amin = amin1 ( amin,amin3 )
      bmin = amin1 ( bmin,bmin3 )
      amax = amax1 ( amax,amax3 )
      bmax = amax1 ( bmax,bmax3 )
c
      amin = amin1 ( amin,amin4 )
      bmin = amin1 ( bmin,bmin4 )
      amax = amax1 ( amax,amax4 )
      bmax = amax1 ( bmax,bmax4 )
c     write ( 23,101 ) amin,amax, bmin,bmax, amin2,amax2, bmin2,bmax2
c
c.....fix it so that grid is undistorted...
c
      dela = amax - amin
      delb = bmax - bmin
c
      if ( dela .gt. delb ) go to 4
c
      b2 = bmax + 0.1*delb
      b1 = bmin - 0.1*delb
      a2 = ( amax + amin + b2 - b1 ) / 2.0
      a1 = ( amax + amin - b2 + b1 ) / 2.0
c
      go to 6
c
    4 continue
c
      a2 = amax + 0.1*dela
      a1 = amin - 0.1*dela
      b2 = ( bmax + bmin + a2 - a1 ) / 2.0
      b1 = ( bmax + bmin - a2 + a1 ) / 2.0
c
    6 continue
c
c     write (23,105 ) a1,a2, b1,b2
 105  format ( 1x, "a1, a2, b1, b2 = ",/, 1p4e12.5 )
c
      call maps ( b1,b2, a1,a2, .114,.828,.285,1.0 )
c
c     write (23,105 ) a1,a2, b1,b2
c
      xl = xx(nl)
      xr = xx(nr)
      dx = ( xr-xl ) / 50.0
      call pointc ( 'x', xmx, 0.0, 1, -1, -1, 0., 0. )
      if ( a .lt. -100.0 )
     $     call pointc ( 'x', -xmx, 0.0, 1, -1, -1, 0., 0. )
      call pointc ( '0', xma, zma, 1, -1, -1, 0., 0. )
      if ( a .lt. -100.0 )
     $     call pointc ( '0', -xma, zma, 1, -1, -1, 0., 0. )
      call pointc ( '*', wcentr, 0.0, 1, -1, -1, 0., 0. )
      call pointc ( 'X', chgt(1),  chgt(2), 1, -1, -1, 0., 0. )
      call pointc ( 'X', chgt(3),  chgt(4), 1, -1, -1, 0., 0. )
c
      if ( idot .eq. 0 ) then
c
         call point ( xx(nl),arr(nl) )
c
         do i = nl, nr
            x = xx(i)
            y = arr(i)
            call vector ( x, y )
         end do
c
         call point ( xx2(nl),arr2(nl) )
c
         do i = nl, nr
            x = xx2(i)
            y = arr2(i)
            call vector ( x,y )
         end do
c
         if ( a .lt. -100.0 ) then
            call point ( -xx2(nl),arr2(nl) )
c
            do i = nl, nr
               x = -xx2(i)
               y = arr2(i)
               call vector ( x, y )
            end do
c
         end if
c
      nlr = nr3 - nl3 + 1
      call points ( xx3, arr3, nlr, -1,-1,0.,0. )
c
      nlr = nr4 - nl4 + 1
      call points ( xx4, arr4, nlr, -1,-1,0.,0. )
c
c$$$
c$$$         call point ( xx3(nl3),arr3(nl3) )
c$$$c
c$$$         do i = nl3, nr3
c$$$            x = xx3(i)
c$$$            y = arr3(i)
c$$$            call vector ( x,y )
c$$$         end do
c
         go to 60
c
      end if
c
 15   continue
c
      go to 50
c
c$$$      call setpch ( 0, 0, -1,-1 )
c$$$      do i = nl, nr
c$$$         call pointc ( 1hx, xx(i), arr(i), 1, -1, -1, 0., 0. )
c$$$      end do
c$$$c
c$$$      do i = nl, nr
c$$$         call pointc ( 1h0, xx2(i), arr2(i), 1, -1, -1, 0., 0. )
c$$$      end do
c$$$c
c$$$      do i = nl3, nr3
c$$$         call pointc ( 1h*, xx3(i), arr3(i), 1, -1, -1, 0., 0. )
c$$$      end do
c
 50   continue
c
      nlr = nr - nl + 1
      call points ( xx, arr, nlr, -1,-1,0.,0. )
      call points ( xx2, arr2, nlr, -1,-1,0.,0. )
      nlr = nr3 - nl3 + 1
      call points ( xx3, arr3, nlr, -1,-1,0.,0. )
c
 60   continue
c
      call dders(1)
c     ab = a1 - 0.1 * ( a2-a1 )
c     xm = 0.5 * ( b1+b2 )
c     call setlch ( xm,ab,0,2,0,-1 )
c     call gtext ( xla, -1, 0 )
      xm = b1 - 0.1 * ( b2-b1 )
c     ab = 0.5 * ( a1 + a2 )
      ab =  a1 + 0.4 * ( a2-a1 )
      string1 = xla//'-'//yla//' plane'
      call setlch ( xm,ab,0,2,1,-1 )
      call gtext ( string1, -1, 0 )
      xm = b1 + 0.05 * ( b2-b1 )
      ab = a2 - 0.05 * ( a2 - a1 )
      call setlch ( xm,ab, 0,1,0,-1 )
      call gtext ( jobid, -1, 0 )
      xm = b1 - (b2-b1) / 10.0
      ab = a1 - (a2-a1) / 8.0
c
      call setlch ( xm,ab, 0,1,0,-1 )
c
      write( 6,'("Circumference of Plasma, pl = ",1pe11.3)')
     $     pcircm
      write(23,'("Circumference of Plasma, pl = ",1pe11.3)')
     $     pcircm
c
      write ( stringv, 3001 ) plrad, xmx, xma, pcircm
 3001 format ( 1x, "plrad= ", f7.3, 1x, "xmaj-[x]= ",f7.3, 1x,
     $     "xma-[0]= ", f7.3, " Pcircum= ", f7.3 )
      call wrtstr ( stringv, 1 )
c
      write ( stringv, 3002) ishape, wcentr, a, b, aw, bw, cw, dw, tw,
     $     abulg, bbulg, tbulg
 3002 format ( 1x, "ishape= ", i3, " wcentr-[*]=",f6.3,
     $     " a= ",f8.3, " b= ",f7.3,/,
     $     1x, "aw= ",f8.3," bw= ",f8.3,
     $     " cw= ",f8.3," dw= ",f8.3," tw= ",f8.3,/,
     $     1x, "abulg= "f8.3, " bbulg= ",f8.3, " tbulg= ", f8.3 )
      call wrtstr ( stringv, 3 )
c
      write ( stringv, 3003) pelong, pdness, wcircm,
     $     wmaj,wlrad, welong, wdness,
     $     chgt(1), chgt(2), chgt(3), chgt(4)
 3003 format ( 1x, "pelong= ",f8.3, " pdness= ",f8.3,
     $     " wcircm= ",f8.3,/,
     $     " wmaj= ",f8.3, " wlrad = ",f8.3, " welong= ",f8.3,
     $     " wdness= ",f8.3, /,
     $     " Coil Tips = ", 2f8.3, " and ", 2f8.3 )
      call wrtstr ( stringv, 3 )
c
      if ( if .eq. 1 ) call framep(jobid, ff)
c
      return
      end
c.......................................................................
      subroutine drawc3(xx,arr,xx2,arr2,xx3,arr3,nl,nr,nl3,nr3,
     $     yla,xla,xmx,xma,zma,
     $     plrad,if,im,idot, jobid, ishape,a,b,aw,bw,cw,dw,tw,
     $     abulg, bbulg, tbulg, wcentr, pcircm, wcircm,
     $     pelong, pdness, wmaj,wlrad, welong, wdness, chgt )
c.......................................................................
c
c     To plot a function arr of xx, from xx(nl) to xx(nr).
c     Called by wwall in module vacuum_vac.f
c
      dimension arr(*), xx(*), arr2(*), xx2(*), arr3(*), xx3(*), chgt(*)
c
      character*(*) yla, xla, jobid
      character*(80) stringv(5), string1
c
c     find minimum and maximum of arr and xx...
      amin = 1.0e20
      amax = -1.0e20
      bmin = 1.0e20
      bmax = -1.0e20
      amin2 =  1.0e20
      amax2 = -1.0e20
      bmin2 =  1.0e20
      bmax2 = -1.0e20
      amin3 =  1.0e20
      amax3 = -1.0e20
      bmin3 =  1.0e20
      bmax3 = -1.0e20
c
      do  nn= nl, nr
         if ( arr(nn) .lt. amin )  amin = arr(nn)
         if ( arr(nn) .gt. amax )  amax = arr(nn)
         if ( xx(nn) .lt. bmin ) bmin = xx(nn)
         if ( xx(nn) .gt. bmax ) bmax = xx(nn)
c
         if ( arr2(nn) .lt. amin2 ) amin2 = arr2(nn)
         if ( arr2(nn) .gt. amax2 ) amax2 = arr2(nn)
         if ( xx2(nn) .lt. bmin2 ) bmin2 = xx2(nn)
         if ( xx2(nn) .gt. bmax2 ) bmax2 = xx2(nn)
      end do
c
      do nn = nl3, nr3
         if ( arr3(nn) .lt. amin3 ) amin3 = arr3(nn)
         if ( arr3(nn) .gt. amax3 ) amax3 = arr3(nn)
         if ( xx3(nn) .lt. bmin3 ) bmin3 = xx3(nn)
         if ( xx3(nn) .gt. bmax3 ) bmax3 = xx3(nn)
      end do
c
c     write ( 23,101 ) amin,amax, bmin,bmax, amin2,amax2, bmin2,bmax2
 101  format ( 1x,"amin,amax, bmin,bmax, amin2,amax2, bmin2,bmax2 = ",
     .     /, 1p8e12.5 )
c
      amin = amin1 ( amin,amin2 )
      bmin = amin1 ( bmin,bmin2 )
      amax = amax1 ( amax,amax2 )
      bmax = amax1 ( bmax,bmax2 )
c
      amin = amin1 ( amin,amin3 )
      bmin = amin1 ( bmin,bmin3 )
      amax = amax1 ( amax,amax3 )
      bmax = amax1 ( bmax,bmax3 )
c
c     write ( 23,101 ) amin,amax, bmin,bmax, amin2,amax2, bmin2,bmax2
c
c.....fix it so that grid is undistorted...
c
      dela = amax - amin
      delb = bmax - bmin
c
      if ( dela .gt. delb ) go to 4
c
      b2 = bmax + 0.1*delb
      b1 = bmin - 0.1*delb
      a2 = ( amax + amin + b2 - b1 ) / 2.0
      a1 = ( amax + amin - b2 + b1 ) / 2.0
c
      go to 6
c
    4 continue
c
      a2 = amax + 0.1*dela
      a1 = amin - 0.1*dela
      b2 = ( bmax + bmin + a2 - a1 ) / 2.0
      b1 = ( bmax + bmin - a2 + a1 ) / 2.0
c
    6 continue
c
c     write (23,105 ) a1,a2, b1,b2
 105  format ( 1x, "a1, a2, b1, b2 = ",/, 1p4e12.5 )
c
      call maps ( b1,b2, a1,a2, .114,.828,.285,1.0 )
c
c     write (23,105 ) a1,a2, b1,b2
c
      xl = xx(nl)
      xr = xx(nr)
      dx = ( xr-xl ) / 50.0
      call pointc ( 'x', xmx, 0.0, 1, -1, -1, 0., 0. )
      if ( a .lt. -100.0 )
     $     call pointc ( 'x', -xmx, 0.0, 1, -1, -1, 0., 0. )
      call pointc ( '0', xma, zma, 1, -1, -1, 0., 0. )
      if ( a .lt. -100.0 )
     $     call pointc ( '0', -xma, zma, 1, -1, -1, 0., 0. )
      call pointc ( '*', wcentr, 0.0, 1, -1, -1, 0., 0. )
      call pointc ( 'X', chgt(1),  chgt(2), 1, -1, -1, 0., 0. )
      call pointc ( 'X', chgt(3),  chgt(4), 1, -1, -1, 0., 0. )
c
      if ( idot .eq. 0 ) then
c
         call point ( xx(nl),arr(nl) )
c
         do i = nl, nr
            x = xx(i)
            y = arr(i)
            call vector ( x, y )
         end do
c
         call point ( xx2(nl),arr2(nl) )
c
         do i = nl, nr
            x = xx2(i)
            y = arr2(i)
            call vector ( x,y )
         end do
c
         if ( a .lt. -100.0 ) then
            call point ( -xx2(nl),arr2(nl) )
c
            do i = nl, nr
               x = -xx2(i)
               y = arr2(i)
               call vector ( x, y )
            end do
c
         end if
c
      nlr = nr3 - nl3 + 1
      call points ( xx3, arr3, nlr, -1,-1,0.,0. )
c$$$
c$$$         call point ( xx3(nl3),arr3(nl3) )
c$$$c
c$$$         do i = nl3, nr3
c$$$            x = xx3(i)
c$$$            y = arr3(i)
c$$$            call vector ( x,y )
c$$$         end do
c
         go to 60
c
      end if
c
 15   continue
c
      go to 50
c
c$$$      call setpch ( 0, 0, -1,-1 )
c$$$      do i = nl, nr
c$$$         call pointc ( 1hx, xx(i), arr(i), 1, -1, -1, 0., 0. )
c$$$      end do
c$$$c
c$$$      do i = nl, nr
c$$$         call pointc ( 1h0, xx2(i), arr2(i), 1, -1, -1, 0., 0. )
c$$$      end do
c$$$c
c$$$      do i = nl3, nr3
c$$$         call pointc ( 1h*, xx3(i), arr3(i), 1, -1, -1, 0., 0. )
c$$$      end do
c
 50   continue
c
      nlr = nr - nl + 1
      call points ( xx, arr, nlr, -1,-1,0.,0. )
      call points ( xx2, arr2, nlr, -1,-1,0.,0. )
      nlr = nr3 - nl3 + 1
      call points ( xx3, arr3, nlr, -1,-1,0.,0. )
c
 60   continue
c
      call dders(1)
c     ab = a1 - 0.1 * ( a2-a1 )
c     xm = 0.5 * ( b1+b2 )
c     call setlch ( xm,ab,0,2,0,-1 )
c     call gtext ( xla, -1, 0 )
      xm = b1 - 0.1 * ( b2-b1 )
c     ab = 0.5 * ( a1 + a2 )
      ab =  a1 + 0.4 * ( a2-a1 )
      string1 = xla//'-'//yla//' plane'
      call setlch ( xm,ab,0,2,1,-1 )
      call gtext ( string1, -1, 0 )
      xm = b1 + 0.05 * ( b2-b1 )
      ab = a2 - 0.05 * ( a2 - a1 )
      call setlch ( xm,ab, 0,1,0,-1 )
      call gtext ( jobid, -1, 0 )
      xm = b1 - (b2-b1) / 10.0
      ab = a1 - (a2-a1) / 8.0
c
      call setlch ( xm,ab, 0,1,0,-1 )
c
      write( 6,'("Circumference of Plasma, pl = ",1pe11.3)')
     $     pcircm
      write(23,'("Circumference of Plasma, pl = ",1pe11.3)')
     $     pcircm
c
      write ( stringv, 3001 ) plrad, xmx, xma, pcircm
 3001 format ( 1x, "plrad= ", f7.3, 1x, "xmaj-[x]= ",f7.3, 1x,
     $     "xma-[0]= ", f7.3, " Pcircum= ", f7.3 )
      call wrtstr ( stringv, 1 )
c
      write ( stringv, 3002) ishape, wcentr, a, b, aw, bw, cw, dw, tw,
     $     abulg, bbulg, tbulg
 3002 format ( 1x, "ishape= ", i3, " wcentr-[*]=",f6.3,
     $     " a= ",f8.3, " b= ",f7.3,/,
     $     1x, "aw= ",f8.3," bw= ",f8.3,
     $     " cw= ",f8.3," dw= ",f8.3," tw= ",f8.3,/,
     $     1x, "abulg= "f8.3, " bbulg= ",f8.3, " tbulg= ", f8.3 )
      call wrtstr ( stringv, 3 )
c
      write ( stringv, 3003) pelong, pdness, wcircm,
     $     wmaj,wlrad, welong, wdness,
     $     chgt(1), chgt(2), chgt(3), chgt(4)
 3003 format ( 1x, "pelong= ",f8.3, " pdness= ",f8.3,
     $     " wcircm= ",f8.3,/,
     $     " wmaj= ",f8.3, " wlrad = ",f8.3, " welong= ",f8.3,
     $     " wdness= ",f8.3, /,
     $     " Coil Tips = ", 2f8.3, " and ", 2f8.3 )
      call wrtstr ( stringv, 3 )
c
      if ( if .eq. 1 ) call framep(jobid, ff)
c
      return
      end
c
c....................................................................
c.......................................................................
      subroutine drawc2(xx,arr,xx2,arr2,nl,nr,yla,xla,xmx,xma,zma,
     $     plrad,if,im,idot, jobid, ishape,a,b,aw,bw,cw,dw,tw,
     $     abulg, bbulg, tbulg, wcentr, pcircm, wcircm,
     $     pelong, pdness, wmaj,wlrad, welong, wdness )
c.......................................................................
c
c         To plot a function arr of xx, from xx(nl) to xx(nr).
c         Called by wwall in module vacuum_vac.f
c
      dimension arr(*), xx(*), arr2(*), xx2(*)
c
      character*(*) yla, xla, jobid
      character*(80) stringv(5), string1
c
c         find minimum and maximum of arr and xx...
      amin = 1.0e20
      amax = -1.0e20
      bmin = 1.0e20
      bmax = -1.0e20
      amin2 =  1.0e20
      amax2 = -1.0e20
      bmin2 =  1.0e20
      bmax2 = -1.0e20
      do 5 nn= nl, nr
      if ( arr(nn) .lt. amin )  amin = arr(nn)
      if ( arr(nn) .gt. amax )  amax = arr(nn)
      if ( xx(nn) .lt. bmin ) bmin = xx(nn)
      if ( xx(nn) .gt. bmax ) bmax = xx(nn)
c
      if ( arr2(nn) .lt. amin2 ) amin2 = arr2(nn)
      if ( arr2(nn) .gt. amax2 ) amax2 = arr2(nn)
      if ( xx2(nn) .lt. bmin2 ) bmin2 = xx2(nn)
      if ( xx2(nn) .gt. bmax2 ) bmax2 = xx2(nn)
c
    5 continue
c
c
c     write ( 23,101 ) amin,amax, bmin,bmax, amin2,amax2, bmin2,bmax2
 101  format ( 1x,"amin,amax, bmin,bmax, amin2,amax2, bmin2,bmax2 = ",
     .      /, 1p8e12.5 )
c
      amin = amin1 ( amin,amin2 )
      bmin = amin1 ( bmin,bmin2 )
      amax = amax1 ( amax,amax2 )
      bmax = amax1 ( bmax,bmax2 )
c
c     write ( 23,101 ) amin,amax, bmin,bmax, amin2,amax2, bmin2,bmax2
c
c..... fix it so that grid is undistorted...
c
      dela = amax - amin
      delb = bmax - bmin
c
      if ( dela .gt. delb ) go to 4
c
      b2 = bmax + 0.1*delb
      b1 = bmin - 0.1*delb
      a2 = ( amax + amin + b2 - b1 ) / 2.0
      a1 = ( amax + amin - b2 + b1 ) / 2.0
c
      go to 6
c
    4 continue
c
      a2 = amax + 0.1*dela
      a1 = amin - 0.1*dela
      b2 = ( bmax + bmin + a2 - a1 ) / 2.0
      b1 = ( bmax + bmin - a2 + a1 ) / 2.0
c
    6 continue
c
c     write (23,105 ) a1,a2, b1,b2
 105  format ( 1x, "a1, a2, b1, b2 = ",/, 1p4e12.5 )
c
      call maps ( b1,b2, a1,a2, .114,.828,.285,1.0 )
c
c     write (23,105 ) a1,a2, b1,b2
c
      xl = xx(nl)
      xr = xx(nr)
      dx = ( xr-xl ) / 50.0
      call pointc ( 'x', xmx, 0.0, 1, -1, -1, 0., 0. )
      if ( a .lt. -100.0 )
     $     call pointc ( 'x', -xmx, 0.0, 1, -1, -1, 0., 0. )
      call pointc ( '0', xma, zma, 1, -1, -1, 0., 0. )
      if ( a .lt. -100.0 )
     $     call pointc ( '0', -xma, zma, 1, -1, -1, 0., 0. )
      call pointc ( '*', wcentr, 0.0, 1, -1, -1, 0., 0. )
c
      if ( idot .eq. 0 ) then
c
         call point ( xx(nl),arr(nl) )
c
         do 10 i = nl, nr
            x = xx(i)
            y = arr(i)
            call vector ( x, y )
c
 10      continue
c
         call point ( xx2(nl),arr2(nl) )
c
         do 12 i = nl, nr
c
            x = xx2(i)
            y = arr2(i)
            call vector ( x,y )
c
 12      continue
c
         if ( a .lt. -100.0 ) then
            call point ( -xx2(nl),arr2(nl) )
c
            do 13 i = nl, nr
               x = -xx2(i)
               y = arr2(i)
               call vector ( x, y )
 13         continue
         end if
c
         go to 60
c
      end if
c
 15   continue
c
      go to 50
c
c$$$      call setpch ( 0, 0, -1,-1 )
c$$$      do 20 i = nl, nr
c$$$         call pointc ( 1hx, xx(i), arr(i), 1, -1, -1, 0., 0. )
c$$$ 20   continue
c$$$c
c$$$      do 22 i = nl, nr
c$$$c
c$$$         call pointc ( 1h0, xx2(i), arr2(i), 1, -1, -1, 0., 0. )
c$$$c
c$$$ 22   continue
c$$$c
 50   continue
c
      nlr = nr - nl + 1
      call points ( xx, arr, nlr, -1,-1,0.,0. )
      call points ( xx2, arr2, nlr, -1,-1,0.,0. )
c
 60   continue
c
      call dders(1)
c      ab = a1 - 0.1 * ( a2-a1 )
c      xm = 0.5 * ( b1+b2 )
c      call setlch ( xm,ab,0,2,0,-1 )
c      call gtext ( xla, -1, 0 )
      xm = b1 - 0.1 * ( b2-b1 )
c      ab = 0.5 * ( a1 + a2 )
      ab =  a1 + 0.4 * ( a2-a1 )
      string1 = xla//'-'//yla//' plane'
      call setlch ( xm,ab,0,2,1,-1 )
      call gtext ( string1, -1, 0 )
      xm = b1 + 0.05 * ( b2-b1 )
      ab = a2 - 0.05 * ( a2 - a1 )
      call setlch ( xm,ab, 0,1,0,-1 )
      call gtext ( jobid, -1, 0 )
      xm = b1 - (b2-b1) / 10.0
      ab = a1 - (a2-a1) / 8.0
c
      call setlch ( xm,ab, 0,1,0,-1 )
c
         write( 6,'("Circumference of Plasma, pl = ",1pe11.3)')
     $        pcircm
         write(23,'("Circumference of Plasma, pl = ",1pe11.3)')
     $        pcircm
c
      write ( stringv, 3001 ) plrad, xmx, xma, pcircm
 3001 format ( 1x, "plrad= ", f7.3, 1x, "xmaj-[x]= ",f7.3, 1x,
     $     "xma-[0]= ", f7.3, " Pcircum= ", f7.3 )
      call wrtstr ( stringv, 1 )
c
      write ( stringv, 3002) ishape, wcentr, a, b, aw, bw, cw, dw, tw,
     $     abulg, bbulg, tbulg
 3002 format ( 1x, "ishape= ", i3, " wcentr-[*]=",f6.3,
     $     " a= ",f8.3, " b= ",f7.3,/,
     $     1x, "aw= ",f8.3," bw= ",f8.3,
     $     " cw= ",f8.3," dw= ",f8.3," tw= ",f8.3,/,
     $     1x, "abulg= "f8.3, " bbulg= ",f8.3, " tbulg= ", f8.3 )
      call wrtstr ( stringv, 3 )
c
      write ( stringv, 3003) pelong, pdness, wcircm,
     $     wmaj,wlrad, welong, wdness
 3003 format ( 1x, "pelong= ",f8.3, " pdness= ",f8.3,
     $     " wcircm= ",f8.3,/,
     $     " wmaj= ",f8.3, " wlrad = ",f8.3, " welong= ",f8.3,
     $     " wdness= ",f8.3 )
      call wrtstr ( stringv, 3 )
c
      if ( if .eq. 1 ) call framep(jobid, ff)
c
      return
      end
c
c....................................................................
      subroutine drawc1(xx,arr,nl,nr,yla,xla,xmx,xma,zma,
     $     plrad,if,im,idot, jobid, ishape,a,b,aw,bw,cw,dw,tw,
     $     abulg, bbulg, tbulg, wcentr, pcircm,
     $     pelong, pdness, wmaj,wlrad, welong, wdness )
c....................................................................
c
c         To plot a function arr of xx, from xx(nl) to xx(nr).
c         Called by wwall in module vacuum_vac.f
c
      dimension arr(*), xx(*)
c
      character*(*) yla, xla, jobid
      character*(80) stringv(5), string1
c
c         find minimum and maximum of arr and xx...
      amin = 1.0e20
      amax = -1.0e20
      bmin = 1.0e20
      bmax = -1.0e20
      amin2 =  1.0e20
      amax2 = -1.0e20
      bmin2 =  1.0e20
      bmax2 = -1.0e20
      do 5 nn= nl, nr
      if ( arr(nn) .lt. amin )  amin = arr(nn)
      if ( arr(nn) .gt. amax )  amax = arr(nn)
      if ( xx(nn) .lt. bmin ) bmin = xx(nn)
      if ( xx(nn) .gt. bmax ) bmax = xx(nn)
c
    5 continue
c
c
c     write ( 23,101 ) amin,amax, bmin,bmax
 101  format ( 1x,"amin,amax, bmin,bmax = ",
     .      /, 1p4e12.5 )
c
c..... fix it so that grid is undistorted...
c
      dela = amax - amin
      delb = bmax - bmin
c
      if ( dela .gt. delb ) go to 4
c
      b2 = bmax + 0.1*delb
      b1 = bmin - 0.1*delb
      a2 = ( amax + amin + b2 - b1 ) / 2.0
      a1 = ( amax + amin - b2 + b1 ) / 2.0
c
      go to 6
c
    4 continue
c
      a2 = amax + 0.1*dela
      a1 = amin - 0.1*dela
      b2 = ( bmax + bmin + a2 - a1 ) / 2.0
      b1 = ( bmax + bmin - a2 + a1 ) / 2.0
c
    6 continue
c
c     write (23,105 ) a1,a2, b1,b2
 105  format ( 1x, "a1, a2, b1, b2 = ",/, 1p4e12.5 )
c
      call maps ( b1,b2, a1,a2, .114,.828,.285,1.0 )
c
c     write (23,105 ) a1,a2, b1,b2
c
      xl = xx(nl)
      xr = xx(nr)
      dx = ( xr-xl ) / 50.0
      call pointc ( 'x', xmx, 0.0, 1, -1, -1, 0., 0. )
      if ( a .lt. -100.0 )
     $     call pointc ( 'x', -xmx, 0.0, 1, -1, -1, 0., 0. )
      call pointc ( '0', xma, zma, 1, -1, -1, 0., 0. )
      if ( a .lt. -100.0 )
     $     call pointc ( '0', -xma, zma, 1, -1, -1, 0., 0. )
      call pointc ( '*', wcentr, 0.0, 1, -1, -1, 0., 0. )
c
      if ( idot .eq. 0 ) then
c
         call point ( xx(nl),arr(nl) )
c
         do 10 i = nl, nr
            x = xx(i)
            y = arr(i)
            call vector ( x, y )
c
 10      continue
c
         go to 60
c
      end if
c
 15   continue
c
      go to 50
c
c$$$      call setpch ( 0, 0, -1,-1 )
c$$$      do 20 i = nl, nr
c$$$         call pointc ( 1hx, xx(i), arr(i), 1, -1, -1, 0., 0. )
c$$$ 20   continue
c
 50   continue
c
      nlr = nr - nl + 1
      call points ( xx, arr, nlr, -1,-1,0.,0. )
c
 60   continue
c
      call dders(1)
c      ab = a1 - 0.1 * ( a2-a1 )
c      xm = 0.5 * ( b1+b2 )
c      call setlch ( xm,ab,0,2,0,-1 )
c      call gtext ( xla, -1, 0 )
      xm = b1 - 0.1 * ( b2-b1 )
c      ab = 0.5 * ( a1 + a2 )
      ab =  a1 + 0.4 * ( a2-a1 )
      string1 = xla//'-'//yla//' plane'
      call setlch ( xm,ab,0,2,1,-1 )
      call gtext ( string1, -1, 0 )
      xm = b1 + 0.05 * ( b2-b1 )
      ab = a2 - 0.05 * ( a2 - a1 )
      call setlch ( xm,ab, 0,1,0,-1 )
      call gtext ( jobid, -1, 0 )
      xm = b1 - (b2-b1) / 10.0
      ab = a1 - (a2-a1) / 8.0
c
      call setlch ( xm,ab, 0,1,0,-1 )
c
         write( 6,'("Circumference of Plasma, pl = ",1pe11.3)')
     $        pcircm
         write(23,'("Circumference of Plasma, pl = ",1pe11.3)')
     $        pcircm
c
      write ( stringv, 3001 ) plrad, xmx, xma, pcircm
 3001 format ( 1x, "plrad= ", f7.3, 1x, "xmaj-[x]= ",f7.3, 1x,
     $     "xma-[0]= ", f7.3, " Pcircum= ", f7.3 )
      call wrtstr ( stringv, 1 )
c
      write ( stringv, 3002) ishape, wcentr, a, b, aw, bw, cw, dw, tw,
     $     abulg, bbulg, tbulg
 3002 format ( 1x, "ishape= ", i3, " wcentr-[*]=",f6.3,
     $     " a= ",f8.3, " b= ",f7.3,/,
     $     1x, "aw= ",f8.3," bw= ",f8.3,
     $     " cw= ",f8.3," dw= ",f8.3," tw= ",f8.3,/,
     $     1x, "abulg= "f8.3, " bbulg= ",f8.3, " tbulg= ", f8.3 )
      call wrtstr ( stringv, 3 )
c
      write ( stringv, 3003) pelong, pdness, wmaj,wlrad, welong, wdness
 3003 format ( 1x, "pelong= ",f8.3, " pdness= ",f8.3,/,
     $     " wmaj= ",f8.3, " wlrad = ",f8.3, " welong= ",f8.3,
     $     " wdness= ",f8.3 )
      call wrtstr ( stringv, 2 )
c
      if ( if .eq. 1 ) call framep(jobid, ff)

c$$$      call atpoint ( "LEAVING DRAWC1", "nlr", nlr, "pcircm",pcircm,
c$$$     $     iotty, outmod )
c
      return
      end
c
c......................................................
c
c
      subroutine wrtgr1 ( jobid, xw, yw, cmax, n, q,
     $     omsq, lfm,xil,xlab,jmax1 )
c
      real n
      dimension lfm(*), xil(*)
c
      character*(80) stringv(11), string1(5)
      character*(*) jobid, xlab
c
      call setlch ( xw, yw, 1,0,0,0 )
c
      write ( string1, '(a)' ) jobid
      call wrtstr ( string1, 1 )
c
      write ( string1, 1002 ) cmax
      call wrtstr ( string1, 1 )
c
      write( stringv, 1001) q, n, omsq, xlab
      call wrtstr ( stringv, 4)
c
      do j = 1, jmax1
         write ( string1, 1003 ) lfm(j), xil(j)
         call wrtstr ( string1, 1 )
      end do
c
 1002 format(1x, "cmx =", 1pe9.3 )
 1001 format(1x,"qs =",f5.2,/,2x,"n =",f4.1,/
     .     1x,"o=",e12.5,/,3x,"l",3x, a )
 1003 format ( i4, 1pe11.3 )
c
      return
      end
c
c.....................................................
      subroutine dart(ax0,ay0,vx,vy,zdum)
      data pi,r,s, epszer/ 3.14159265,.2,.15, 1.0e-30/
c
c     all passed parameters must be real*4
c     velocities are assumed scaled to coordinates
c----------------------------------------------------------
c     Pass in ax0, ay0 and define ax, ay such that arrow is
c     centered.  Put in epszer. M.Chance
c
      ax = ax0 - vx/2.0
      ay = ay0 - vy/2.0
c----------------------------------------------------------
      bx = ax + vx
      by = ay + vy
c      if ( abs(vx) .le. epszer ) vx = sign ( epszer, vx )
      theta = atan2m(vy,vx)
      psi = theta + pi/2.
      ab = sqrt(vx**2 + vy**2)
      rab = r*ab
      sab = s*ab
      delx = sab*cos(psi)
      dely = sab*sin(psi)
      xbar = ax + (1.-r)*vx
      ybar = ay + (1.-r)*vy
      cx = xbar - delx
      cy = ybar - dely
      dx = xbar + delx
      dy = ybar + dely
      call setcrt(ax,ay)
      call vector(bx,by)
      call vector(cx,cy)
      call setcrt(dx,dy)
      call vector(bx,by)
      return
      end
c$$$
c$$$c
c$$$c
c$$$      subroutine dart ( x0, y0, vx, vy, h0 )
c$$$c
c$$$      x = x0 - vx
c$$$      y = y0 - vy
c$$$      u = x0 + vx
c$$$      v = y0 + vy
c$$$      h = h0 * sqrt ( sqrt(vx**2+vy**2) )
c$$$c
c$$$c     call plotv(x,y,u,v,h)
c$$$c
c$$$      return
c$$$      end
c
c.......................................................
c
c
      subroutine pospl1 ( x,y,jmin,jmax, ipos, label, l )
c
      dimension x(*), y(*)
c
      character*(*) label
      character*(80) string1(5)
c
c.. Zeps is related to the precision of the plot package. For single
c    precision use 
      zeps = 1.0e-6

      CALL bounds ( x,y, jmin,jmax, xmin,xmax, ymin,ymax )

      jd = jmax - jmin + 1

      if ( abs(ymax-ymin) .le. zeps )
     $     ymax = ymin + 1.0
      if ( abs(xmax-xmin) .le. zeps )
     $     xmax = xmin + 1.0
      dx = xmax - xmin
      dy = ymax - ymin
      xw = xmin + dx/20.0
      yw = ymax + dy/20.0
c
      go to (110,120,130,140), ipos
c
  110 call mapg(xmin,xmax, ymin,ymax, 0.12,0.42, 0.70,0.97)
      go to 200
  120 call mapg(xmin,xmax, ymin,ymax, 0.5,0.8, 0.70,0.97)
      go to 200
  130 call mapg(xmin,xmax, ymin,ymax, 0.12,0.42, 0.27,0.54)
      go to 200
  140 call mapg(xmin,xmax, ymin,ymax, 0.5,0.8, 0.27,0.54)
c
  200 continue
c
      call trace ( x, y, jd, 1,1,0.,0. )
      call setlch ( xw,yw, 0,1,0,-1 )

      lenz = LEN_TRIM(label)
c
      write ( string1, 300 ) label(1:lenz), l
      call wrtstr ( string1, 1 )
  300 format ( a, ' z=',i5 )
c
      return
      end
c
c........................................................
      SUBROUTINE pospl2 ( x,y,y2, jmin,jmax, ipos, label, l )
c............................................................
      dimension x(*), y(*), y2(*)
c
      character*(*) label
      character*(80) string1(5)
c
         CALL bounds ( x,y, jmin,jmax, xmin,xmax, ymin,ymax )
         CALL bounds ( x,y2, jmin,jmax, xmin,xmax, y2mn,y2mx )

c$$$      ymax = -1.0e40
c$$$      ymin =  1.0e40
c$$$      y2mx = -1.0e40
c$$$      y2mn =  1.0e40
c
      jd = jmax - jmin + 1
c
c$$$      do 50 j = jmin, jmax
c$$$c
c$$$      if ( y(j) .gt. ymax ) ymax = y(j)
c$$$      if ( y(j) .lt. ymin ) ymin = y(j)
c$$$      if ( y2(j) .gt. y2mx ) y2mx = y2(j)
c$$$      if ( y2(j) .lt. y2mn ) y2mn = y2(j)
c$$$c
c$$$   50 continue
c
c$$$      xmin = x(jmin)
c$$$      xmax = x(jmax)
      ymax = amax1 ( ymax,y2mx )
      ymin = amin1 ( ymin,y2mn )
      if ( ymax .eq. ymin ) ymax = ymin + 1.0
      dx = xmax - xmin
      dy = ymax - ymin
      xw = xmin + dx/20.0
      yw = ymax + dy/20.0
c
      go to (110,120,130,140), ipos
c
  110 call mapg(xmin,xmax, ymin,ymax, 0.12,0.42, 0.70,0.97)
      go to 200
  120 call mapg(xmin,xmax, ymin,ymax, 0.5,0.8, 0.70,0.97)
      go to 200
  130 call mapg(xmin,xmax, ymin,ymax, 0.12,0.42, 0.27,0.54)
      go to 200
  140 call mapg(xmin,xmax, ymin,ymax, 0.5,0.8, 0.27,0.54)
c
  200 continue
c
      call trace ( x, y, jd, 1,1,0.,0. )
      call tracep ( x, y2, jd, 4, 1,1 )
      call setlch ( xw,yw, 0,1,0,-1 )
c
      lenz = LEN_TRIM( label )
      write ( string1, 300 ) label(1:lenz), l
      call wrtstr ( string1, 1 )
  300 format ( a, 1x, i5 )
c
      return
      end
c
c........................................................
c
      SUBROUTINE pospl2s ( x,y,y2, jmin,jmax, ipos, label, l )
c...................................................
c
      DIMENSION x(*), y(*), y2(*)
c
      CHARACTER*(*) label
      CHARACTER*(80) string1(5)
c

      CALL bounds ( x,y, jmin,jmax, xmin,xmax, ymin,ymax )
      CALL bounds ( x,y2, jmin,jmax, xmin,xmax, y2mn,y2mx )

c$$$      ymax = -1.0e40
c$$$      ymin =  1.0e40
c$$$      y2mx = -1.0e40
c$$$      y2mn =  1.0e40
c
      jd = jmax - jmin + 1
c
c$$$      DO  j = jmin, jmax
c$$$c
c$$$         if ( y(j) .gt. ymax ) ymax = y(j)
c$$$         if ( y(j) .lt. ymin ) ymin = y(j)
c$$$         if ( y2(j) .gt. y2mx ) y2mx = y2(j)
c$$$         if ( y2(j) .lt. y2mn ) y2mn = y2(j)
c$$$c
c$$$      END DO
c
c$$$      xmin = x(jmin)
c$$$      xmax = x(jmax)
      ymax = amax1 ( ymax,y2mx )
      ymin = amin1 ( ymin,y2mn )
      IF ( ymax == ymin ) ymax = ymin + 1.0
      dx = xmax - xmin
      dy = ymax - ymin
      xw = xmin - dx/18.0
      yw = ymax + dy/20.0
c
      GO TO (110,120,130,140), ipos
c
  110 CALL mapg(xmin,xmax, ymin,ymax, 0.114,0.374, 0.80,0.97)
      go to 200
  120 CALL mapg(xmin,xmax, ymin,ymax, 0.454,0.714, 0.80,0.97)
      go to 200
  130 CALL mapg(xmin,xmax, ymin,ymax, 0.12,0.42, 0.53,0.70)
      go to 200
  140 CALL mapg(xmin,xmax, ymin,ymax, 0.5,0.8, 0.53,0.70)
c
  200 CONTINUE
c
      CALL trace ( x, y, jd, 1,1,0.,0. )
      CALL tracep ( x, y2, jd, 4, 1,1 )
      CALL setlch ( xw,yw, 0,1,0,-1 )
c
      lenz = LEN_TRIM(label)
      write ( string1, 300 ) label(1:lenz), l
      CALL wrtstr ( string1, 1 )
  300 FORMAT ( a, i4 )
c
      RETURN
      END
c
c........................................................
c
      SUBROUTINE pospl2sc ( x,y,y2, jmin,jmax, ipos, label,
     $     zxpt,zzpt,npt, l )
c...................................................
c
      DIMENSION x(*), y(*), y2(*), zxpt(*), zzpt(*)
c
      CHARACTER*(*) label
      CHARACTER*(80) string1(5)
c

      CALL bounds ( x,y, jmin,jmax, xmin,xmax, ymin,ymax )
      CALL bounds ( x,y2, jmin,jmax, xmin,xmax, y2mn,y2mx )
c$$$
c$$$      ymax = -1.0e40
c$$$      ymin =  1.0e40
c$$$      y2mx = -1.0e40
c$$$      y2mn =  1.0e40
c
      jd = jmax - jmin + 1
c
c$$$      DO  j = jmin, jmax
c$$$c
c$$$         if ( y(j) .gt. ymax ) ymax = y(j)
c$$$         if ( y(j) .lt. ymin ) ymin = y(j)
c$$$         if ( y2(j) .gt. y2mx ) y2mx = y2(j)
c$$$         if ( y2(j) .lt. y2mn ) y2mn = y2(j)
c$$$c
c$$$      END DO
c$$$c
c$$$      xmin = x(jmin)
c$$$      xmax = x(jmax)

      ymax = amax1 ( ymax,y2mx )
      ymin = amin1 ( ymin,y2mn )
      IF ( ymax == ymin ) ymax = ymin + 1.0
      dx = xmax - xmin
      dy = ymax - ymin
      xw = xmin - dx/18.0
      yw = ymax + dy/20.0
c
      GO TO (110,120,130,140), ipos
c
  110 CALL mapg(xmin,xmax, ymin,ymax, 0.10,0.40, 0.80,0.97)
      go to 200
  120 CALL mapg(xmin,xmax, ymin,ymax, 0.50,0.80, 0.80,0.97)
      go to 200
  130 CALL mapg(xmin,xmax, ymin,ymax, 0.10,0.40, 0.53,0.70)
      go to 200
  140 CALL mapg(xmin,xmax, ymin,ymax, 0.50,0.80, 0.53,0.70)
c
  200 CONTINUE
c
      CALL trace ( x, y, jd, 1,1,0.,0. )
      CALL tracep ( x, y2, jd, 4, 1,1 )
      DO ipt = 1, npt
         CALL pointc ( 'X', zxpt(ipt),zzpt(ipt), 1, -1, -1, 0., 0. )
      END DO

      CALL setlch ( xw,yw, 0,1,0,-1 )

      lenz = LEN_TRIM(label)

      WRITE ( string1, 300 ) label(1:lenz), l
      CALL wrtstr ( string1, 1 )
  300 FORMAT ( a, i4 )
c

      xw = xmin + dx/10
      yw = ymin - dy/5.0

      CALL setlch ( xw, yw, 1,1,0,-1 )

      WRITE ( string1, '()' )
      CALL wrtstr ( string1, 1 )

      WRITE ( string1, '( 5x,"K phase")' )
      CALL wrtstr ( string1, 1 )

      DO it = 1, npt
         WRITE ( string1, '( 1x, es13.5 )' )
     $        zzpt(it)
      CALL wrtstr ( string1, 1 )
      END DO

      RETURN
      END
c
c...................................................
c
c............................................................
      subroutine posplv2 ( x,y,y2,zwkx,zwky, jmin,jmax,
     $     ipos, label, l )
c............................................................
c
      dimension x(*), y(*), y2(*), zwkx(*), zwky(*)
c
      character*(*) label
      character*(80) string1(5)
c

      CALL bounds ( x,y, jmin,jmax, xmin,xmax, ymin,ymax )
      CALL bounds ( x,y2, jmin,jmax, xmin,xmax, y2mn,y2mx )

c$$$      ymax = -1.0e40
c$$$      ymin =  1.0e40
c$$$      y2mx = -1.0e40
c$$$      y2mn =  1.0e40
c
      jd = jmax - jmin + 1
c
c$$$      do j = jmin, jmax
c$$$         zwkx(j-jmin+1) = x(j)
c$$$         zwky(j-jmin+1) = y(j)
c$$$c
c$$$         if ( y(j) .gt. ymax ) ymax = y(j)
c$$$         if ( y(j) .lt. ymin ) ymin = y(j)
c$$$         if ( y2(j) .gt. y2mx ) y2mx = y2(j)
c$$$         if ( y2(j) .lt. y2mn ) y2mn = y2(j)
c$$$      end do
c$$$c
c$$$      xmin = x(jmin)
c$$$      xmax = x(jmax)

      ymax = amax1 ( ymax,y2mx )
      ymin = amin1 ( ymin,y2mn )
      if ( ymax .eq. ymin ) ymax = ymin + 1.0
      dx = xmax - xmin
      dy = ymax - ymin
      xw = xmin + dx/20.0
      yw = ymax + dy/20.0
c
      go to (110,120,130,140), ipos
c
 110  call mapg(xmin,xmax, ymin,ymax, 0.12,0.42, 0.70,0.97)
      go to 200
 120  call mapg(xmin,xmax, ymin,ymax, 0.5,0.8, 0.70,0.97)
      go to 200
 130  call mapg(xmin,xmax, ymin,ymax, 0.12,0.42, 0.27,0.54)
      go to 200
 140  call mapg(xmin,xmax, ymin,ymax, 0.5,0.8, 0.27,0.54)
c
 200  continue
c
c      call pointc ( 'x', zwkx, zwky, jd, -1,-1,0.,0. )
      call trace ( zwkx, zwky, jd, 1,1,0.,0. )
c
      do j = jmin, jmax
         zwky(j-jmin+1) = y2(j)
      end do
c
c      call pointc ( 'o', zwkx, zwky, jd, -1,-1,0.,0. )
      call tracep (  zwkx, zwky, jd, 4, 1,1 )
      call setlch ( xw,yw, 0,1,0,-1 )
c
      lenz = LEN_TRIM(label)

      write ( string1, 300 ) label(1:lenz), jd
      call wrtstr ( string1, 1 )
 300  format ( a, ' N=',i5 )

      RETURN
      END
c
c...................................................
c
