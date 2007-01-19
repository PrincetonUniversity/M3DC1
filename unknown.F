      subroutine local(xi,zi,b,theta,itri)
      use basic
      implicit none
      integer itri, iz, ix, nodeids(4)
      real xi, zi, theta, b, small, co, sn, xmin, zmin
      dimension xi(*),zi(*)
      double precision coords(3)
      small = 1.e-8
!
!     calculates the coefficients xi and zi for region ll in the local
!     coordinate expansion:  x = Sum{ xi(i)*si**mi(i)*eta**ni(i)}; {i=1,3}
!                            z = Sum{ zi(i)*si**mi(i)*eta**ni(i)}; {i=1,3}

      call getmincoord(xmin,zmin)
      call nodfac(itri,nodeids)
      call xyznod(nodeids(2),coords)
      
      co = cos(theta)
      sn = sin(theta)
      xi(1) = coords(1)
      xi(2) =  co
      xi(3) = -sn
      call xyznod(nodeids(1),coords)
      zi(1) = coords(2)
      zi(2) =  sn
      zi(3) =  co
      return
      end
!============================================================
      subroutine whattri(x,z,itri,xref,zref)
      use basic
      implicit none
      integer itri, nodeids(4)
      real x, z, xref, zref
      double precision coords(3), x2, z2

      x2 = x
      z2 = z

      call usesearchstructure(x2,z2,itri)

      if(itri.lt.0) return

      call nodfac(itri,nodeids)
      call xyznod(nodeids(1), coords)
      xref = coords(1)
      zref = coords(2)

!      if(itri.eq.1) then
!         write(*,*) "x1,y1=", xref, zref
!      endif
!
      return
      end
!============================================================
      subroutine axis(dum,xsep,zsep,itime)
      use basic
      use p_data
      use t_data
      implicit none
      integer nsep, itrit, nmax,jrect, irect, itri, inews, iii, ii,
     &     i, index, k, inewt, itrinew, isep, itime, 
     &     numelms, nodeids(4), ibegin, iendplusone
      real sumx, sumz, x1, z1, wlocal, dum, x, z, theta, b, co, 
     &     sn, si, eta, sum, avector, sum1, sum2, sum3, sum4, sum5, 
     &     term1, term2, term3, term4, term5, pt, pt1, pt2, p11, 
     &     p22, p12, denom, sinew, etanew, xnew, znew, xsep, zsep, deex
      dimension xsep(*),zsep(*)
      dimension dum(*),wlocal(20),avector(20)
      real xcord, zcord
      double precision coords(3)
!
!     itime=0 for first call
!          =1 for subsequent calls
!     
!     locates the magnetic axis and separatrix locations
!
!.....evaluate the solution to get the value [ans] at one point (x,z)
!      write(*,*) 'did not check on how deex was changed by acbauer!'
!      call safestop(435) 

      call numfac(numelms)

!     loop over all triangles
      nsep = 0

      do itrit=1,numelms

         itri = itrit

         call nodfac(itri,nodeids)
         inews = 0
 300     inews = inews+1
         call safestop(111)
c     below is not tested by acbauer
         x1 = coords(1)
         z1 = coords(2)
!
!.....construct the 18 vector corresponding to this triangle
!.....calculate the index and local coordinates for this triangle
         sumx = 0.
         sumz = 0.
         do iii=1,3 
           call entdofs(numvar, ist(itri,iii)+1, 0, ibegin, iendplusone)
            do ii=1,6
               i = (iii-1)*6 + ii
               index = ibegin+ii-1
               wlocal(i) = dum(index)
            enddo
            index = ist(itri,iii) + 1
            call xyznod(index,coords)
            sumx = sumx + coords(1)
            sumz = sumz + coords(2)
         enddo
         if(inews.le.1) then
            x = sumx/3.
            z = sumz/3.
         endif

         ! calculate local coordinates
         theta = ttri(itri)
         b = btri(itri)
         co = cos(theta)
         sn = sin(theta)
         si  = (x-x1)*co + (z-z1)*sn - b
         eta =-(x-x1)*sn + (z-z1)*co

         !     calculate the function value corresponding to this point
         do i=1,20
            sum = 0.
            do k=1,18
               sum = sum + gtri(i,k,itri)*wlocal(k)
            enddo
            avector(i) = sum
         enddo
         inewt = 0
 301     continue
         inewt = inewt + 1
!
!.....evaluate the polynomial and second derivative
         sum = 0.
         sum1 = 0.
         sum2 = 0.
         sum3 = 0.
         sum4 = 0.
         sum5 = 0.
         do i=1,20
            sum = sum + avector(i)*si**mi(i)*eta**ni(i)
            term1 = 0.
            if(mi(i).ge.1) term1 = mi(i)*si**(mi(i)-1)*eta**ni(i)
            term2 = 0.
        if(ni(i).ge.1) term2 = ni(i)*si**mi(i)*eta**(ni(i)-1)
        term3 = 0.
        if(mi(i).ge.2) term3 = mi(i)*(mi(i)-1)*si**(mi(i)-2)*eta**ni(i)
        term4 = 0.
        if(ni(i).ge.2) term4 = ni(i)*(ni(i)-1)*si**mi(i)*eta**(ni(i)-2)
        term5 = 0.
        if(ni(i)*mi(i) .ge. 1)                                          &
     &                 term5 = mi(i)*ni(i)*si**(mi(i)-1)*eta**(ni(i)-1)
        sum1 = sum1 + avector(i)*term1
        sum2 = sum2 + avector(i)*term2
        sum3 = sum3 + avector(i)*term3
        sum4 = sum4 + avector(i)*term4
        sum5 = sum5 + avector(i)*term5
      enddo
      pt  = sum
      pt1  = sum1
      pt2  = sum2
      p11 = sum3
      p22 = sum4
      p12 = sum5
!
      denom =  p22*p11 - p12**2
      sinew = si -  ( p22*pt1 - p12*pt2)/denom
      etanew= eta - (-p12*pt1 + p11*pt2)/denom
      xnew = x1 + co*(b+sinew) - sn*etanew
      znew = z1 + sn*(b+sinew) + co*etanew
!
!.....determine if this new point is nearby
      call whattri(xnew,znew,itrinew,x1,z1)
      call getdeex(itri,deex)
      if( (xnew-x)**2 + (znew-z)**2 .le. 2*deex**2                      &
     &   .and. xnew .gt. 0 .and. xnew.lt.4.                             &
     &   .and. znew .gt. 0 .and. znew.lt.4) then
            x = xnew
            z = znew
            itri = itrinew
            irect = (x/deex) + 1
            jrect = (z/deex) + 1
            if(inews.le.4 )go to 300
            if(inews.eq.5) then
!
              if(nsep.eq.0) then
                nsep = 1
                xsep(nsep) = x
                zsep(nsep) = z
              else
                do isep = 1,nsep
                if(sqrt((x-xsep(isep))**2+(z-zsep(isep))**2)            &
     &           .le. .001*deex) go to 400
                enddo
                nsep = nsep + 1
                xsep(nsep) = x
                zsep(nsep) = z
 400            continue
              endif
!
            endif
 3001       format(1p2e12.4)
      endif
!
      enddo
!
!.....end of loop over all triangles
!     return
      end
