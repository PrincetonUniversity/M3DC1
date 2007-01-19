      subroutine errorcalc(numvar, phi, itype)
      implicit none
      integer numelms, ipolyorder, nodeids(4), i, ilength, 
     &     numpoints, m, n, numvar, itype, iodd, j, i2, itri
      parameter(ilength=20)
      real error, jac, weights(ilength), points(2,ilength), 
     &     elmerror, deex, deez, x, z, truvalue, elmtruvalue,
     &     simvalue, elmsimvalue, svalue2, phi(*), svalue, tvalue,
     &     weight
      double precision coords(9)

      call compload()
      ipolyorder = 5            ! polynomial order to be integrated

      error = 0.
      simvalue = 0.
      truvalue = 0.
      call integrationpoints(ipolyorder, numpoints, points, weights, 
     &     ilength)
      call getparams(n, m, deex, deez)
      do iodd=1,2
         do i=1,m-1
            do j=1,m-1
               coords(1) = deex*(i-1.)
               coords(2) = deez*(j-1.)
               if(iodd .eq. 1) then
                  coords(4) = deex*i
                  coords(5) = deez*(j-1.)
                  coords(7) = deex*i
                  coords(8) = deez*j
               else
                  coords(4) = deex*i
                  coords(5) = deez*j
                  coords(7) = deex*(i-1.)
                  coords(8) = deez*j
               endif
               call jacobian(coords, jac)
c     integration points loop
               elmerror = 0.
               elmtruvalue = 0.
               elmsimvalue = 0.
               do i2=1, numpoints
                  call parametric2global(coords, points(1,i2), 
     &                 points(2,i2), x, z)
                  call compevaluate(x, z, tvalue, itype)
                  call evaluate(x, z, svalue, svalue2, phi, itype,
     &                 numvar)
                  weight = weights(i2)*jac*0.5
                  elmerror = elmerror + (tvalue-svalue)*
     &                 (tvalue-svalue)*weight
                  elmtruvalue = elmtruvalue + tvalue*tvalue*weight
                  elmsimvalue = elmsimvalue + svalue*svalue*weight
               enddo            ! do i2=1,numpoints
               error = error + elmerror
               simvalue = simvalue + elmsimvalue
               truvalue = truvalue + elmtruvalue
            enddo
         enddo
      enddo

      simvalue = sqrt(simvalue)
      truvalue = sqrt(truvalue)
      error = sqrt(error)
      write(*,*) 'the L2 norm of the simulation value is',simvalue
      write(*,*) 'the L2 norm of the true value is',truvalue
      write(*,*) 'the L2 norm of the error and relative error are',
     &     error,(error/truvalue)

c      call comptecplot()
      call compfree()
      
      return
      end

      subroutine integrationpoints(ipolyorder, numpoints, points, 
     &     weights, ilength)
      implicit none
      integer maxnumpoints
      parameter(maxnumpoints=7)
      integer ipolyorder, numpoints, ilength,
     &     lnumpoints(0:maxnumpoints), i
      real points(2,ilength), weights(ilength), onethird, onesixth, 
     &     twothirds
      data lnumpoints /1,1,3,4,6,7,12,13/

      if(ipolyorder .lt. 0 .or. ipolyorder .gt. maxnumpoints) then
         write(*,*) 'wrong integration order ', ipolyorder
         call safestop(444)   
      endif
      numpoints = lnumpoints(ipolyorder)

      if(numpoints .gt. ilength) then
         write(*,*) 'points and weights arrays need to be larger'
         call safestop(445)
      endif

      onethird = 1./3.
      twothirds = 2./3.
      onesixth = 1./6.

      if(numpoints .eq. 1) then
         points(1,1) = onethird
         points(2,1) = onethird
         weights(1) = 1.
      else if(numpoints .eq. 3) then
         points(1,1) = onesixth
         points(2,1) = onesixth
         points(1,2) = twothirds
         points(2,2) = onesixth
         points(1,3) = onesixth
         points(2,3) = twothirds
         weights(1) = onethird
         weights(2) = onethird
         weights(3) = onethird
      else if(numpoints .eq. 4) then
         points(1,1) = onethird
         points(2,1) = onethird
         points(1,2) = .6
         points(2,2) = .2
         points(1,3) = .2
         points(2,3) = .6
         points(1,4) = .2
         points(2,4) = .2
         weights(1) = -0.5625
         weights(2) = 1.5625*onethird
         weights(3) = 1.5625*onethird
         weights(4) = 1.5625*onethird
      else if(numpoints .eq. 6) then
         points(1,1) = 0.816847572980459
         points(2,1) = 0.091576213509771
         weights(1) = 0.109951743655322
         points(1,2) = 0.091576213509771
         points(2,2) = 0.816847572980459
         weights(2) = 0.109951743655322
         points(1,3) = 0.091576213509771
         points(2,3) = 0.091576213509771
         weights(3) = 0.109951743655322
         points(1,4) = 0.108103018168070
         points(2,4) = 0.445948490915965
         weights(4) = 0.223381589678011
         points(1,5) = 0.445948490915965
         points(2,5) = 0.108103018168070
         weights(5) = 0.223381589678011
         points(1,6) = 0.445948490915965
         points(2,6) = 0.445948490915965
         weights(6) = 0.223381589678011
      else if(numpoints .eq. 7) then
         points(1,1) = onethird
         points(2,1) = onethird
         weights(1) = 0.225000000000000
         points(1,2) = 0.797426985353087
         points(2,2) = 0.101286507323456
         weights(2) = 0.125939180544827
         points(1,3) = 0.101286507323456
         points(2,3) = 0.797426985353087
         weights(3) = 0.125939180544827
         points(1,4) = 0.101286507323456
         points(2,4) = 0.101286507323456
         weights(4) = 0.125939180544827
         points(1,5) = 0.470142064105115
         points(2,5) = 0.059715871789770
         weights(5) = 0.132394152788506
         points(1,6) = 0.059715871789770
         points(2,6) = 0.470142064105115
         weights(6) = 0.132394152788506
         points(1,7) = 0.470142064105115
         points(2,7) = 0.470142064105115
         weights(7) = 0.132394152788506
      else if(numpoints .eq. 12) then
         points(1,1) = 0.873821971016996
         points(2,1) = 0.063089014491502
         weights(1) = 0.050844906370207
         points(1,2) = 0.063089014491502
         points(2,2) = 0.873821971016996
         weights(2) = 0.050844906370207
         points(1,3) = 0.063089014491502
         points(2,3) = 0.063089014491502
         weights(3) = 0.050844906370207
         points(1,4) = 0.501426509658179
         points(2,4) = 0.249286745170910
         weights(4) = 0.116786275726379
         points(1,5) = 0.249286745170910
         points(2,5) = 0.501426509658179
         weights(5) = 0.116786275726379
         points(1,6) = 0.249286745170910
         points(2,6) = 0.249286745170910
         weights(6) = 0.116786275726379
         points(1,7) = 0.636502499121399
         points(2,7) = 0.310352451033785
         weights(7) = 0.082851075618374
         points(1,8) = 0.310352451033785
         points(2,8) = 0.636502499121399
         weights(8) = 0.082851075618374
         points(1,9) = 0.636502499121399
         points(2,9) = 0.053145049844816
         weights(9) = 0.082851075618374
         points(1,10) = 0.310352451033785
         points(2,10) = 0.053145049844816
         weights(10) = 0.082851075618374
         points(1,11) = 0.053145049844816
         points(2,11) = 0.310352451033785
         weights(11) = 0.082851075618374
         points(1,12) = 0.053145049844816
         points(2,12) = 0.636502499121399
         weights(12) = 0.082851075618374
      else if(numpoints .eq. 13) then
         points(1,1) = 0.333333333333333
         points(2,1) = 0.333333333333333
         weights(1) = -0.149570044467682
         points(1,2) = 0.479308067841920
         points(2,2) = 0.260345966079040
         weights(2) = 0.175615257433208
         points(1,3) = 0.260345966079040
         points(2,3) = 0.479308067841920
         weights(3) = 0.175615257433208
         points(1,4) = 0.260345966079040
         points(2,4) = 0.260345966079040
         weights(4) = 0.175615257433208
         points(1,5) = 0.869739794195568
         points(2,5) = 0.065130102902216
         weights(5) = 0.053347235608838
         points(1,6) = 0.065130102902216
         points(2,6) = 0.869739794195568
         weights(6) = 0.053347235608838
         points(1,7) = 0.065130102902216
         points(2,7) = 0.065130102902216
         weights(7) = 0.053347235608838
         points(1,8) = 0.048690315425316
         points(2,8) = 0.312865496004874
         weights(8) = 0.077113760890257
         points(1,9) = 0.312865496004874
         points(2,9) = 0.048690315425316
         weights(9) = 0.077113760890257
         points(1,10) = 0.638444188569810
         points(2,10) = 0.048690315425316
         weights(10) = 0.077113760890257
         points(1,11) = 0.048690315425316
         points(2,11) = 0.638444188569810
         weights(11) = 0.077113760890257
         points(1,12) = 0.312865496004874
         points(2,12) = 0.638444188569810
         weights(12) = 0.077113760890257
         points(1,13) = 0.638444188569810
         points(2,13) = 0.312865496004874
         weights(13) = 0.077113760890257
      endif

      return 
      end

      
      subroutine jacobian(coords, jac)
      implicit none
      
      double precision coords(9)
      real v1(2), v2(2), jac

      v1(1) = coords(4) - coords(1)
      v1(2) = coords(5) - coords(2)
      v2(1) = coords(7) - coords(1)
      v2(2) = coords(8) - coords(2)
      
      jac = v1(1)*v2(2) - v1(2)*v2(1)

      return
      end

      subroutine parametric2global(coords, xl, zl, xg, zg)
c     takes local (xl,zl) coords and maps to global (xg,zg)
      implicit none
      integer i
      real xl, zl, xg, zg, w
      double precision coords(9)

      w = 1.-xl-zl
      if(w .lt. 0.) w = 0.
c$$$      xg = w*coords(0) + xl*coords(4) + zl*coords(7)
c$$$      zg = w*coords(1) + xl*coords(5) + zl*coords(8)

      xg = w*coords(1) + xl*coords(5) + zl*coords(8)
      zg = w*coords(2) + xl*coords(6) + zl*coords(9)


      return
      end
