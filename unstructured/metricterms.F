c=============================================
c FOURTH ORDER TENSORS - in alphabetical order
c=============================================
      subroutine calcc1term2(itri, c1term2, fintl, dofs1, dofs2,
     &     ind1, ind2, inddof1, inddof2)
      use basic
      use t_data
      implicit none
      real c1term2(18,18), fintl(-6:maxi,-6:maxi),c1term2p(20,20), 
     &     dofs1(20), dofs2(20), sum
      integer itri, msum, nsum, ind1, ind2, inddof1, inddof2
      integer, pointer :: i1, i2, idof1, idof2, p, q, r, s
      include 'set.h'

      allocate(p, q, r, s)
      call set4(i1, i2, idof1, idof2, p, q, r, s, ind1, ind2, inddof1)

      c1term2p = 0.
      do p=1,20
         do q=1,20
            do r=1,20
               do s=1,20
                  msum = mi(p)+mi(q)+mi(r)+mi(s)
                  nsum = ni(p)+ni(q)+ni(r)+ni(s)
                  sum =                                            
     &                 (mi(r)*(mi(r)-1)*((mi(s)*ni(q)-mi(q)*ni(s))              
     &                 *((mi(q)+mi(s)-1)*ni(p)-(ni(q)+ni(s)-1)*mi(p))         
     &                 -(mi(s)*ni(p)-mi(p)*ni(s))              
     &                 *((mi(p)+mi(s)-1)*ni(q)-(ni(p)+ni(s)-1)*mi(q)))        
     &                 -mi(q)*(mi(q)-1)*(mi(s)*ni(p)-mi(p)*ni(s))             
     &                 *(mi(r)*(ni(p)+ni(s)-1)-ni(r)*(mi(p)+mi(s)-1))     
     &                 +2*(mi(s)*ni(p)-mi(p)*ni(s))*mi(r)*mi(q)           
     &                 *((mi(r)-1)*ni(q)-ni(r)*(mi(q)-1)))
     &                 *fintl(msum-4,nsum-2)
                                !
     &                 +(ni(r)*(ni(r)-1)*((mi(s)*ni(q)-mi(q)*ni(s))       
     &                 *((mi(q)+mi(s)-1)*ni(p)-(ni(q)+ni(s)-1)*mi(p))     
     &                 -(mi(s)*ni(p)-mi(p)*ni(s))              
     &                 *((mi(p)+mi(s)-1)*ni(q)-(ni(p)+ni(s)-1)*mi(q)))    
     &                 -ni(q)*(ni(q)-1)*(mi(s)*ni(p)-mi(p)*ni(s))              
     &                 *(mi(r)*(ni(p)+ni(s)-1)-ni(r)*(mi(p)+mi(s)-1))     
     &                 +2*(mi(s)*ni(p)-mi(p)*ni(s))*ni(r)*ni(q)           
     &                 *(mi(r)*(ni(q)-1)-(ni(r)-1)*mi(q)))
     &                 *fintl(msum-2,nsum-4)
                  c1term2p(i1,i2) = c1term2p(i1,i2) 
     &                 + sum*dofs1(idof1)*dofs2(idof2)
               enddo
            enddo
         enddo
      enddo
      call transform2(c1term2p, c1term2, itri)
      deallocate(p, q, r, s)

      return
      end
c============================================
      subroutine calcv1terml(itri, v1terml, fintl, dofs1, dofs2,
     &     ind1, ind2, inddof1, inddof2)
      use basic
      use t_data
      implicit none
      real v1terml(18,18), fintl(-6:maxi,-6:maxi),v1termlp(20,20), 
     &     dofs1(20), dofs2(20), tmp
      integer itri, msum, nsum, ind1, ind2, inddof1, inddof2
      integer, pointer :: i1, i2, idof1, idof2, p, q, r, s
      include 'set.h'
      allocate(p, q, r, s)

      call set4(i1, i2, idof1, idof2, p, q, r, s, ind1, ind2, inddof1)
      
      v1termlp = 0.
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tmp = mi(q)*(mi(r)*((mi(q)-1)*(ni(s)+ni(p))-ni(q)*mi(s))       
     &             - ni(r)*mi(p)*(mi(q)-1))*fintl(msum-3,nsum-1)         &
     &       +ni(q)*(ni(r)*(mi(q)*ni(s)-(ni(q)-1)*(mi(s)+mi(p)))        &
     &             + mi(r)*ni(p)*(ni(q)-1))*fintl(msum-1,nsum-3)
              v1termlp(i1,i2) = v1termlp(i1,i2) +
     &             tmp*dofs1(idof1)*dofs2(idof2)
            enddo
          enddo
        enddo
      enddo

      call transform2(v1termlp, v1terml, itri)

      deallocate(p, q, r, s)
      return
      end
c============================================
      subroutine calcv2terml(itri, v2terml, fintl, dofs1, dofs2,
     &     ind1, ind2, inddof1, inddof2)
      use basic
      use t_data
      implicit none
      real v2terml(18,18), fintl(-6:maxi,-6:maxi),v2termlp(20,20), 
     &     dofs1(20), dofs2(20), tmp
      integer itri, msum, nsum, ind1, ind2, inddof1, inddof2
      integer, pointer :: i1, i2, idof1, idof2, p, q, r, s
      include 'set.h'
      allocate(p, q, r, s)

      call set4(i1, i2, idof1, idof2, p, q, r, s, ind1, ind2, inddof1)
      
      v2termlp = 0.
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tmp =  mi(r)*(mi(q)*((mi(r)-1)*(ni(s)+ni(p))-ni(r)*mi(s))        &
     &             - ni(q)*mi(p)*(mi(r)-1))*fintl(msum-3,nsum-1)         &
     &       +ni(r)*(ni(q)*(mi(r)*ni(s)-(ni(r)-1)*(mi(s)+mi(p)))        &
     &             + mi(q)*ni(p)*(ni(r)-1))*fintl(msum-1,nsum-3)
              v2termlp(i1,i2) = v2termlp(i1,i2) + 
     &             tmp*dofs1(idof1)*dofs2(idof2)
            enddo
          enddo
        enddo
      enddo

      call transform2(v2termlp, v2terml, itri)

      deallocate(p, q, r, s)
      return
      end

c============================================
c THIRD ORDER TENSORS - in alphabetical order
c============================================

!============================================================
      subroutine calcg0term(itri, g0term, fintl, dofs, ind1, 
     &     ind2, inddof)
      use basic
      use t_data
      implicit none
      real g0term(18,18), fintl(-6:maxi,-6:maxi),g0termp(20,20), 
     &     dofs(20), tmp
      integer itri, msum, nsum, ind1, ind2, inddof
      integer, pointer :: i1, i2, idof, p, q, r
      include 'set.h'
      allocate(p, q, r)

      call set3(i1, i2, idof, p, q, r, ind1, ind2)

      g0termp = 0.

      do p=1,20
         do q=1,20
            do r=1,20
               msum = mi(q)+mi(r)+mi(p)
               nsum = ni(q)+ni(r)+ni(p)
               tmp = - (mi(p)*ni(r)-mi(r)*ni(p))                 
     &              *(mi(q)*(mi(q)-1)*fintl(msum-3,nsum-1) 
     &              + ni(q)*(ni(q)-1)*fintl(msum-1,nsum-3))
               g0termp(i1,i2) = g0termp(i1,i2) + tmp*dofs(idof)
            enddo
         enddo
      enddo

      call transform2(g0termp, g0term, itri)
      
      deallocate(p,q,r)
      return
      end
c================================================
      subroutine calcg12erm(itri, g12erm, fintl, dofs, ind1, 
     &     ind2, inddof)
      use basic
      use t_data
      implicit none
      
      real g12erm(18,18), fintl(-6:maxi,-6:maxi),g12ermp(20,20),
     &     dofs(20), tmp
      integer itri,msum, nsum, ind1, ind2, inddof
      integer, pointer :: i1, i2, idof, p, q, r
      include 'set.h'
      allocate(p, q, r)

      call set3(i1, i2, idof, p, q, r, ind1, ind2)
      
      g12ermp = 0.
      do p=1,20
         do q=1,20
            do r=1,20
               msum = mi(q)+mi(r)+mi(p)
               nsum = ni(q)+ni(r)+ni(p)
               tmp = 
     &        mi(q)*(mi(q)-1)*(mi(p)+mi(r))*(mi(p)+mi(r)-1)             &
     &                 *fintl(msum-4,nsum)                               &
     &      +(mi(q)*(mi(q)-1)*(ni(p)+ni(r))*(ni(p)+ni(r)-1)             &
     &       +ni(q)*(ni(q)-1)*(mi(p)+mi(r))*(mi(p)+mi(r)-1))            &
     &                 *fintl(msum-2,nsum-2)                             &
     &       +ni(q)*(ni(q)-1)*(ni(p)+ni(r))*(ni(p)+ni(r)-1)             &
     &                 *fintl(msum,nsum-4)
               g12ermp(i1,i2) = g12ermp(i1,i2) + tmp*dofs(idof)
            enddo
         enddo
      enddo

      call transform2(g12ermp, g12erm, itri)
               
      deallocate(p, q, r)

      return
      end
c================================================
      subroutine calcg2term(itri, g2term, fintl, dofs, ind1, 
     &     ind2, inddof)
      use basic
      use t_data
      implicit none
      
      real g2term(18,18),fintl(-6:maxi,-6:maxi),g2termp(20,20),dofs(20)
      integer itri,msum, nsum, ind1, ind2, inddof
      integer, pointer :: i1, i2, idof, p, q, r
      include 'set.h'
      allocate(p, q, r)

      call set3(i1, i2, idof, p, q, r, ind1, ind2)
      g2termp = 0.
      do p=1,20
         do q=1,20
            do r=1,20
               msum = mi(q)+mi(r)+mi(p)
               nsum = ni(q)+ni(r)+ni(p)
               g2termp(i1,i2) =  g2termp(i1,i2) -
     &              (mi(p)*mi(q)*fintl(msum-2,nsum)
     &              +ni(p)*ni(q)*fintl(msum,nsum-2))*dofs(idof)
            enddo
         enddo
      enddo

      call transform2(g2termp, g2term, itri)

      deallocate(p, q, r)
      
      return
      end
c================================================
      subroutine calck0term(itri, k0term, fintl, dofs, ind1, 
     &     ind2, inddof)
      use basic
      use t_data
      implicit none
      real k0term(18,18), fintl(-6:maxi,-6:maxi),k0termp(20,20),
     &     dofs(*)
      integer msum, nsum, ind1, ind2, inddof, itri
      integer, pointer :: i1, i2, idof, p, q, r
      include 'set.h'
      allocate(p, q, r)

      call set3(i1, i2, idof, p, q, r, ind1, ind2)

      k0termp = 0.
      do p=1,20
         do q=1,20
            do r=1,20
               msum = mi(q)+mi(r)+mi(p)
               nsum = ni(q)+ni(r)+ni(p)
c...........poisson bracket of 2 scalars
               k0termp(i1,i2)= k0termp(i1,i2) + dofs(idof)*
     &              (mi(q)*ni(r)-mi(r)*ni(q))*fintl(msum-1,nsum-1)
            enddo
         enddo
      enddo

      call transform2(k0termp, k0term, itri)

      deallocate(p, q, r)
      
      return 
      end
c============================================
c SECOND ORDER TENSORS - in alphabetical order
c============================================
      subroutine calcaterm(itri, aterm, fintl)
      use basic
      use t_data
      implicit none
      integer itri, p, q, i, j, k, l, msum, nsum
      real fintl(-6:maxi,-6:maxi), aterm(18,18), 
     &     amatrix(20,20)
      
!.....begin local integrations of metric elements
      do q=1,20
         do p=1,20
            msum = mi(p)+mi(q)
            nsum = ni(p)+ni(q)
!.......del-sq operator
            amatrix(p,q) =  -(mi(p)*mi(q)*fintl(msum-2,nsum) +               & 
     &           ni(p)*ni(q)*fintl(msum,nsum-2))
         enddo
      enddo
!......transform from local coordinates to global coordinates
      call transform2(amatrix, aterm, itri)

      return
      end
c================================================
      subroutine calcbterm(itri, bterm, fintl)
      use basic
      use t_data
      implicit none
      integer itri, p, q, i, j, k, l, msum, nsum
      real fintl(-6:maxi,-6:maxi), bterm(18,18),  
     &     bmatrix(20,20)
      
!.....begin local integrations of metric elements
      do q=1,20
         do p=1,20
            msum = mi(p)+mi(q)
            nsum = ni(p)+ni(q)
!.......del-fourth operator
        bmatrix(p,q)=(mi(p)*(mi(p)-1)*mi(q)*(mi(q)-1)*fintl(msum-4,nsum)&
     &              + ni(p)*(ni(p)-1)*ni(q)*(ni(q)-1)*fintl(msum,nsum-4)&
     &    + (mi(p)*(mi(p)-1)*ni(q)*(ni(q)-1)                            &
     &    +  mi(q)*(mi(q)-1)*ni(p)*(ni(p)-1))*fintl(msum-2,nsum-2))
         enddo
      enddo
!......transform from local coordinates to global coordinates
      call transform2(bmatrix, bterm, itri)

      return
      end
c================================================
      subroutine calcdterm(itri, dterm, fintl)
      use basic
      use t_data
      implicit none
      integer itri, p, q, i, j, k, l, msum, nsum
      real fintl(-6:maxi,-6:maxi), dterm(18,18),  
     &     dmatrix(20,20)
      
!.....begin local integrations of metric elements
      do q=1,20
         do p=1,20
            msum = mi(p)+mi(q)
            nsum = ni(p)+ni(q)
!.......mass matrix
            dmatrix(p,q) = fintl(msum,nsum)
         enddo
      enddo
!......transform from local coordinates to global coordinates
      call transform2(dmatrix, dterm, itri)

      return
      end
c================================================
      subroutine calcsterm(itri, sterm, fintl)
      use basic
      use t_data
      implicit none
      integer itri, p, q, i, j, k, l, msum, nsum
      real fintl(-6:maxi,-6:maxi), sterm(18,18),  
     &     qmatrix(20,20)
      
!.....begin local integrations of metric elements
      do q=1,20
         do p=1,20
            msum = mi(p)+mi(q)
            nsum = ni(p)+ni(q)
!.......del-star operator
        qmatrix(p,q) =  (mi(p)+mi(q)-1)*mi(p)*fintl(msum-2,nsum)         &
     &                + (ni(p)+ni(q)-1)*ni(p)*fintl(msum,nsum-2)

         enddo
      enddo
!......transform from local coordinates to global coordinates
      call transform2(qmatrix, sterm, itri)

      return
      end
c============================================
c FIRST ORDER TENSORS - in alphabetical order
c================================================
      subroutine calcd2term(itri, d2term, fintl)
      use basic
      use t_data
      implicit none
      integer itri, p, q, i, j, k, l, msum, nsum
      real fintl(-6:maxi,-6:maxi), d2term(18),  
     &     d2matrix(20)
      
!.....begin local integrations of metric elements
      do p=1,20
         d2matrix(p) = fintl(mi(p),ni(p))
      enddo
!......transform from local coordinates to global coordinates
      call transform1(d2matrix, d2term, itri)

      return
      end
c================================================
      subroutine set4(p1, p2, p3, p4, i1, i2, i3, i4,  
     &     ind1, ind2, ind3)
      implicit none
c     subroutine sets p1, p2, p3, and p4 to point
c     to i1, i2, i3, i4 depending on ind1, ind2 and ind3
c     i1, i2, i3 and i4 must already be allocated
      integer, pointer :: p1, p2, p3, p4, i1, i2, i3, i4
      integer ind1, ind2, ind3

      if(ind1 .eq. 1) then
         p1 => i1
         if(ind2 .eq. 2) then
            p2 => i2
            if(ind3 .eq. 3) then
               p3 => i3
               p4 => i4
            else
               p3 => i4
               p4 => i3
            endif
         else 
            p2 => i3
            p3 => i2
         endif
      else if(ind1 .eq. 2) then
         p1 => i2
         if(ind2 .eq. 1) then
            p2 => i1
            p3 => i3
         else
            p2 => i3
            p3 => i1
         endif
      else if(ind1 .eq. 3) then
         p1 => i3
         if(ind2 .eq. 1) then
            p2 => i1
            p3 => i2
         else
            p2 => i2
            p3 => i1
         endif
      else
         call safestop(5532)
      endif

      end subroutine set4
c=========================c===========================
      subroutine set3(p1, p2, p3, i1, i2, i3, ind1, ind2)
      implicit none
c     subroutine sets p1, p2, and p3 to point
c     to i1, i2, i3 depending on ind1 and ind2
c     i1, i2 and i3 must already be allocated
      integer, pointer :: p1, p2, p3, i1, i2, i3
      integer ind1, ind2

      if(ind1 .eq. 1) then
         p1 => i1
         if(ind2 .eq. 2) then
            p2 => i2
            p3 => i3
         else
            p2 => i3
            p3 => i2
         endif
      else if(ind1 .eq. 2) then
         p1 => i2
         if(ind2 .eq. 1) then
            p2 => i1
            p3 => i3
         else
            p2 => i3
            p3 => i1
         endif
      else if(ind1 .eq. 3) then
         p1 => i3
         if(ind2 .eq. 1) then
            p2 => i1
            p3 => i2
         else
            p2 => i2
            p3 => i1
         endif
      else
         call safestop(5532)
      endif

      end subroutine set3
c=========================
c transforms second order tensor from local to global values
      subroutine transform2(termp, term, itri)
      use basic
      use t_data
      implicit none
      integer itri, i, j, k, l
      real termp(20,20), term(18,18), sum, gk

c......transform from local coordinates to global coordinates 
c......for a second order tensor
      term = 0.
      do i=1,18
         do j=1,18
            sum = 0.
            do k = 1,20
               gk = gtri(k,i,itri)
               do l = 1,20
                  sum = sum + termp(k,l)*gk*gtri(l,j,itri)
               enddo
            enddo
            term(i,j) = term(i,j) + sum
         enddo   
      enddo      

      return
      end
c=========================
c transforms first order tensor from local to global values
      subroutine transform1(termp, term, itri)
      use basic
      use t_data
      implicit none
      integer itri, i, k
      real termp(20), term(18)

c......transform from local coordinates to global coordinates 
c......for a first order tensor
      do i=1,18
         term(i) = 0.
         do k = 1,20
            term(i) = term(i) + termp(k)*gtri(k,i,itri)
         enddo   
      enddo      

      return
      end
c=========================
      subroutine calcfint(fint, maxi, atri, btri, ctri)
      implicit none
      integer maxi, i, j
      real fint(-6:maxi,-6:maxi), f, atri,btri,ctri
      
      do i=-6,maxi
         do j=-6,maxi
            fint(i,j) = 0.
            if(i.ge.0 .and.j.ge.0)  then
               fint(i,j) = f(i,j,atri,btri,ctri)
            endif
         enddo
      enddo

      return
      end

         
      
      
      
