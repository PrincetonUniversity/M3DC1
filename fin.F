      subroutine rinvdef(rinv)

      use basic

      implicit none
      integer i, numnodes, ibegin, iendplusone
      real, intent(out) :: rinv(*)
      double precision coords(3)

      real :: x

      call numnod(numnodes)
      do i=1,numnodes
         call xyznod(i, coords)

         x = coords(1) + xzero

         call entdofs(1, i, 0, ibegin, iendplusone)
         rinv(ibegin) =  1./x
         rinv(ibegin+1) = -1./x**2
         rinv(ibegin+3) =  2./x**3
      enddo

      return
      end
!===================================
      subroutine safestop(iarg)

      use basic
      use hdf5_output

      implicit none
      
      integer, intent(in) :: iarg

      include 'mpif.h'

      integer :: ier
      
      ! finalize hdf5
      print *, "finalizing hdf5..."
      call hdf5_finalize(ier)
      print *, "done."

      print *, "finalizing SCOREC software..."
      call clearscorecdata()
      print *, "done."

      call sludexit
      call MPI_Finalize(ier)
      if (ier /= 0) print *,'Error terminating MPI:',ier

      close(UNIT=9, STATUS='KEEP')
      write(*,*) "stopped at", iarg
      stop
      end
!==========================================
      subroutine MPE_DECOMP1D(nproc,mrows,myrank,jfirst,jlast)
      implicit none
      integer nproc, mrows, myrank, jfirst, jlast, ilast, ifirst,
     1     mleft, nleft, i, ifrac
      ilast = 0
      ifirst = 1
      mleft = mrows
      nleft = nproc + 1
      do i=0,myrank
        mleft = mleft - (ilast-ifirst+1)
        nleft = nleft - 1
        ifirst = ilast + 1
        ifrac = mleft/nleft
        if(nleft*ifrac.eq.mleft) then
          ilast = ifirst + ifrac - 1
        else
          ilast = ifirst + ifrac 
        endif
      enddo
      jfirst = ifirst
      jlast  = ilast
      return
      end
!============================================================
      subroutine input
      use basic
      use arrays
      use sparse
      implicit none
!
!.....set itest .ne. 0 for a particular test problem
!     itest = 1 for tilting columns (nonlinear, n=m=21)
!     itest = 2 for Taylor problem  (n=m=31, dt=2  )
!     itest = 3 for tilting columns (linear)
!     itest = 4 for Taylor problem  (Fitzpatrick base case)
!     itest = 5 for tilting columns - two fluid (linear)
!
!.....switch for second order time advance:
      isecondorder = 0
!
!.....linear parameter  0-nonlinear,  1-linear
      linear = 0
!
!.....density advance parameter: 0 no advance,  1-advance density
      idens = 1
!
!.....restart parameter 0-no restart, 1 restart
      irestart = 0
!
!.....associated restart parameter: 0-normal, 1-start time from zero
      istart = 0
!
!.....setup parameter 0-read from disk, 1 compute
!      isetup = 1
!
!.....boundary parameter  0-nonperiodic, 1-periodic
!      iper = 0
!      jper = 0
!
!.....switch for Taylor Problem
      itaylor = 0
!
!.....time integrator: 0 = Crank-Nicholson, 1 = BDF2
      integrator = 0
!
!.....CN implicitness parameter (thimp > 0.5 for stability)
      thimp = 0.5
!
!.....hyper-resistivity coefficients
      hyper = 0.5000
      hyperi= 0.5000
      hyperv= 0.50
      hyperp= 0.50
!
!.....viscosity, resistivity, and heat conduction
      amu = .00010
      etar = .00010
      eta0 = 0.
      kappa = .001
!
!.....loop voltage
      vloop = 0.
!
!......toroidal magnetic field
      bzero = 1.
!
      facw = 1.
!.....multiplies the diffusion [test] terms 
      facd = 1.
!     NOTE:  if facw=1., facd is zeroed after the first cycle
!
!......2-fluid coefficients
      cb = 0.000
      db = 0.00
!
!     masking switch (0 no mask,  1 mask hyper terms at boundary
      imask = 0
!
!     source terms in velcity advance (1 = include, 0 = exclude)
      isources = 1
!
!......maximum in taylor expansion loop (5 is maximum)
      maxs = 5
!.....how many cycles to skip before inverting matrix
      nskip = 1
!
!.....timestep
      dt = 0.1
!.....ntimemax is max time cycles, nimepr is cycles between print cycle
!
      ntimemax = 20
      ntimepr   = 5
!
!.....itor=0 for cylinder, 1 for torus
      itor=0
!
!.....parameters needed for toroidal plasma
      tcuro = 1.
      xmag = 12.1
      zmag = 0.
      xlim = 10.5
      zlim = 0.
!
!.....pressure = p0*(1 + p1*psi + p2*psi**2 + ...)
!                where psi is the normalized poloidal flux
      p0 = 0.01      
      p1 = -1.
      p2 = -2.
      pedge = 0.
!
!     derivative of current wrt normalized poloidal flux
      djdpsi = 0.0
!
!....number of velocity and field variables
      numvar = 3
!
!.....parameters for taylor problem (itaylor=1)
      eps = .01
!
!.....gravitational acceleration
      grav = 0.
!
!.....1 = conserve toroidal flux
      iconstflux = 0
!
      if(itest.eq.0) then
        open(5,file='C1input',form='formatted',status='old')
        read(5,nml=inputnl)
        if(amuc.eq.0) amuc = amu
        go to 100
      endif
      go to (1,2,3,4,5,6,7,8,9,10), itest
!
 1    continue
!.....itest=1 tilting columns (nonlinear)
      itaylor = 0
      dt = 0.05
      facd = 1.
      etar = .0010
      amu = .001
      kappa = .001
      hyper = 0.01
      hyperi= 0.01
      hyperv = 0.01
      hyperc = 0.01
      hyperp= 0.50
      cb = 1.
      db = .0001
      ntimemax = 100
      ntimepr = 1
      thimp = 0.55
!      isetup = 0
      imask = 1
      numvar = 2
      bzero = 1.
      go to 100
!
 2    continue
!.....itest=2 Taylor problem, high resistivity
      itaylor = 1
      dt = 0.5
      facd = 0.
      etar = .0001
      amu  = .00010
      kappa = .0001
      hyper = 1.00
      hyperi = 1.00
      hyperv = 1.0
      hyperc = 1.0
      hyperp= 0.50
      cb = 1.
      db = 0.0
      ntimemax = 80
      ntimepr   = 10
      thimp = 0.6
!      isetup = 1
      numvar = 2
      imask = 0
      bzero = 0.
      go to 100
 3    continue
!.....itest=3  tilting columns--linear
      linear = 1
      itaylor = 0
      dt = 0.10
      facd = 1.
      etar = .0001
      amu = .0001
      kappa = .0001
      hyper = 0.
      hyperi = 0.
      hyperv = 0.0
      hyperc = 0.0
      cb = 1.0
      db = 0.
      ntimemax = 80
      ntimepr = 5
      thimp = 0.80
!      isetup = 0
      imask = 0
!
      go to 100
 4    continue
!.....itest=4 Taylor problem..Fitzpatrick base case with eta=mu=1.e-5
      itaylor = 1
      dt = 2.0
      facd = 0.
      etar = .000010
      amu  = .0000100
      kappa = .00001
      hyper = 0.
      hyperi = 0.
      hyperv = 0.0
      hyperc = 0.0
      hyperp= 0.50
      cb = 0.
      db = 0.
      ntimemax = 40
      ntimepr   = 5
      thimp = 0.60
!      isetup = 1
      go to 100
 5    continue
!.....itest=5 tilting columns--linear-two fluid
      linear = 1
      numvar = 2
      itaylor = 0
      dt = 0.05
      facd = 1.
      etar = .001
      amu = .001
      kappa = .001
      hyper = 1.00
      hyperi = 1.00
      hyperv = 4.0
      hyperc = 4.0
      hyperp= 1.0
      cb = 1.0
      db = 0.2
      ntimemax = 81
      ntimepr = 1
      thimp = 0.51
!      isetup = 0
      irestart = 0  !0=start from scratch;  1=restart
      imask = 0
      bzero = 1.
      beta = 0.0
!
      go to 100
 6    continue
!.....itest=6 Taylor problem..eta=mu=1.e-5, with two-fluid terms
      itaylor = 1
      dt = 2.0
      facd = 0.
      etar = 0.00010
      amu  = 0.0000100
      kappa= 0.00001
      hyper = .0010
      hyperi = .0010
      hyperv = .0010
      hyperc = .0010
      hyperp= 0.50
      cb = 1.0
      db = 0.
      ntimemax = 20
      ntimepr   = 10
      thimp = 0.6
!      isetup = 0
      imask = 1
      go to 100
!
 7    continue
!.....itest=7 tilting columns--non-linear-two fluid
      itaylor = 0
      dt = 0.025
      facd = 1.
      etar = .0010
      amu = .0010
      kappa = .001
      hyper =  1.
      hyperi = 1.
      hyperv = 1.
      hyperc = 1.
      hyperp= 0.50
      cb = 1.0
      db = .4
      ntimemax = 5
      ntimepr = 5
      thimp = 0.55
!      isetup = 0
      imask = 0
      numvar = 2
      bzero = 0.001
      go to 100
!
    8 continue
!.....itest=8 Taylor state
      itaylor = 2
      dt = 0.18
      facd = 1.
      etar = .0001
      amu = .0001
      kappa = .0001
      hyper = .00
      hyperi = 0.0
      hyperv = 0.0
      hyperc = 0.0
      hyperp= 0.50
      cb = 0.3
      db = 0.
      ntimemax = 10
      ntimepr = 2
      thimp = 0.6
!      isetup = 0
      imask = 0
      go to 100
!
    9 continue
!.....itest=9 GEM reconnection
      linear = 0
      itaylor = 3
      dt = .25
      facd = 0.
      etar = .0010
      amu = .0010
      kappa = .001
      hyper = 1.0  !hyper-resistivity coefficient
      hyperi = 1.0 !hyper-resistivity coefficient for I
      hyperv = 2.0 !hyper-viscosity coefficient for V
      hyperc = 0.0 !hyper-viscosity coefficient for chi
      hyperp= 0.50
      cb = 1.0
      db = 0.0001
      ntimemax = 150
      ntimepr = 50
      thimp = .60
!      isetup = 1
      imask = 1
      bzero = 0.
      numvar = 2
      eps = 0.1
      go to 100
!
   10 continue
!.....itest=10 GEM reconnection...(two fluid)
      linear = 0
      itaylor = 3
      dt = 0.25    !time step size
      facd = 0.
      etar = .0010
      amu = .0010
      kappa = .001
      hyper = 2.00  !hyper-resistivity coefficient
      hyperi = 2.00 !hyper-resistivity coefficient for I
      hyperv = 4.0 !hyper-viscosity coefficient for V
      hyperc = 4.0 !hyper-viscosity coefficient for chi
      hyperp= 1.00
      cb = 1.0
      db = 1.0
      ntimemax = 50  !Number of time steps
      ntimepr = 10
      thimp = 0.6
!      isetup = 0  !0=already exists; 1=recalc
      irestart = 1  !0=start from scratch;  1=restart
      imask = 1
      bzero = 0.0
      numvar = 3
      eps = 0.1
      go to 100
 100  continue
!
      call createdofnumbering(numvar1_numbering, iper, jper, 
     &     6, 0, 0, 0, maxdofs1)
      call createdofnumbering(numvar2_numbering, iper, jper, 
     &     12, 0, 0, 0, maxdofs2)
      if(numvar .ge. 3) then
         call createdofnumbering(numvar3_numbering, iper, jper, 
     &        18, 0, 0, 0, maxdofs3)
      endif

      if(myrank.eq.0) write(*,nml=inputnl)
      return
      end
!============================================================
      function f(m,n,a,b,c)
!
!.....f calculates the double integral over a triangle
!     with base a + b and height c  [see ref 2, figure 1]
!
!     f = Int{si**m*eta**n}d(si)d(eta)
!
      implicit none
      real f, a, b, c, anum, denom, fact
      integer m,n

      anum = c**(n+1)*(a**(m+1)-(-b)**(m+1))*fact(m)*fact(n)
      denom = fact(m+n+2)
      f = anum/denom
      return
      end function f
!============================================================
      subroutine rotation(rot,ndim,theta)
      implicit none
      integer ndim, i, j
      real rot(ndim,*),r1(6,6), co, theta, sn

!
!     calculate the rotation matrix R from Table 2 of Ref 2
!
      co = cos(theta)
      sn = sin(theta)
      do i=1,6
        do j=1,6
          r1(i,j) = 0
        enddo
      enddo
!
      r1(1,1) = 1.
!
      r1(2,2) = co
      r1(2,3) = sn
!
      r1(3,2) = -sn
      r1(3,3) = co
!
      r1(4,4) = co**2
      r1(4,5) = 2.*sn*co
      r1(4,6) = sn**2
!
      r1(5,4) = -sn*co
      r1(5,5) = co**2-sn**2
      r1(5,6) = sn*co
!
      r1(6,4) = sn**2
      r1(6,5) = -2*sn*co
      r1(6,6) = co**2
! 
      do i=1,18
        do j=1,18
          rot(i,j) = 0
        enddo
      enddo
!
      do i=1,6
        do j=1,6
          rot(i,j)       = r1(i,j)
          rot(i+6,j+6)   = r1(i,j)
          rot(i+12,j+12) = r1(i,j)
        enddo
      enddo
!
      return
      end
!============================================================

      subroutine tmatrix(ti,ndim,a,b,c)
      implicit none
      integer ndim, i, j, ierrorchk, ifail, info1, info2
      real a, b, c, danaly, det, percent, diff
      real t(20,20),ti(ndim,*),wkspce(9400),ipiv(20)
      ierrorchk = 0

!   define the 20 x 20 Transformation Matrix that enforces the condition that
!   the nomal slope between triangles has only cubic variation..
!...first initialize to zero
!
!
      do i=1,20
        do j=1,20
          t(i,j) = 0.
        enddo
      enddo
!
!...Table 1 of Ref. [2]
!
      t(1,1)   = 1.
      t(1,2)   = -b
      t(1,4)   = b**2
      t(1,7)   = -b**3
      t(1,11)  = b**4
      t(1,16)  = -b**5
!
      t(2,2)   = 1
      t(2,4)   = -2*b
      t(2,7)   = 3*b**2
      t(2,11)  = -4*b**3
      t(2,16)  = 5*b**4
!
      t(3,3)   = 1
      t(3,5)   = -b
      t(3,8)   = b**2
      t(3,12)  = -b**3
!
      t(4,4)   = 2.
      t(4,7)   = -6.*b
      t(4,11)  = 12*b**2
      t(4,16)  = -20*b**3
!
      t(5,5)   = 1.
      t(5,8)   = -2.*b
      t(5,12)  = 3*b**2
!
      t(6,6)   = 2.
      t(6,9)   = -2*b
      t(6,13)  = 2*b**2
      t(6,17)  = -2*b**3
!
      t(7,1)   = 1.
      t(7,2)   = a
      t(7,4)   = a**2
      t(7,7)   = a**3
      t(7,11)  = a**4
      t(7,16)  = a**5
!
      t(8,2)   = 1.
      t(8,4)   = 2*a
      t(8,7)   = 3*a**2
      t(8,11)  = 4*a**3
      t(8,16)  = 5*a**4
!
      t(9,3)   = 1.
      t(9,5)   = a
      t(9,8)   = a**2
      t(9,12)  = a**3
!
      t(10,4)  = 2
      t(10,7)  = 6*a
      t(10,11) = 12*a**2
      t(10,16) = 20*a**3
!
      t(11,5)  = 1.
      t(11,8)  = 2.*a
      t(11,12) = 3*a**2
!
      t(12,6)  = 2.
      t(12,9)  = 2*a
      t(12,13) = 2*a**2
      t(12,17) = 2*a**3
!
      t(13,1)  = 1
      t(13,3)  = c
      t(13,6)  = c**2
      t(13,10) = c**3
      t(13,15) = c**4
      t(13,20) = c**5
!
      t(14,2)  = 1.
      t(14,5)  = c
      t(14,9)  = c**2
      t(14,14) = c**3
      t(14,19) = c**4
!
      t(15,3)  = 1.
      t(15,6)  = 2*c
      t(15,10) = 3*c**2
      t(15,15) = 4*c**3
      t(15,20) = 5*c**4
!
      t(16,4)  = 2.
      t(16,8)  = 2*c
      t(16,13) = 2*c**2
      t(16,18) = 2*c**3
!
      t(17,5)  = 1.
      t(17,9)  = 2*c
      t(17,14) = 3*c**2
      t(17,19) = 4*c**3
!
      t(18,6)  = 2.
      t(18,10) = 6*c
      t(18,15) = 12*c**2
      t(18,20) = 20*c**3
!
      t(19,16) = 5*a**4*c
      t(19,17) = 3*a**2*c**3 - 2*a**4*c
      t(19,18) = -2*a*c**4+3*a**3*c**2
      t(19,19) = c**5-4*a**2*c**3
      t(19,20) = 5*a*c**4
!
      t(20,16) = 5*b**4*c
      t(20,17) = 3*b**2*c**3 - 2*b**4*c
      t(20,18) = 2*b*c**4 - 3*b**3*c**2
      t(20,19) = c**5 - 4*b**2*c**3
      t(20,20) = -5*b*c**4
!
      if(ierrorchk.eq.0) go to 100
!
!.....analytic formula for determinant
      danaly = -64*(a+b)**17*c**20*(a**2+c**2)*(b**2+c**2)
!
!.....calculate determinant using nag
      ifail = 0
      do i=1,20
        do j=1,20
          ti(i,j) = t(i,j)
        enddo
      enddo
      det = 0.
!     call f03aaf(ti,20,20,det,wkspce,ifail)
!
      diff = det - danaly
      percent = 100* diff / danaly
!
      if(abs(percent) .gt. 1.e-12)                                      &
     &write(*,1001) percent
 1001 format("percent error in determinant",1pe12.4)
!
 100  continue
!
!
!.....calculate the inverse of t using NAG library routines
      info1 = 0
      info2 = 0
      do i=1,20
        do j=1,20
          ti(i,j) = t(i,j)
        enddo
      enddo
      call f07adf(20,20,ti,20,ipiv,info1)
      call f07ajf(20,ti,20,ipiv,wkspce,400,info2)
      if(info1.ne.0.or.info2.ne.0)write(*,1002) info1,info2
 1002 format(3i5)
!
! error checking

      return
      end
!============================================================
      subroutine offset(itri)
      use basic
      use t_data
      implicit none
      integer itri, nodeids(4)
      call nodfac(itri, nodeids)
      ist(itri,1) = nodeids(1)-1
      ist(itri,2) = nodeids(2)-1
      ist(itri,3) = nodeids(3)-1    
      
  !    write(*,*) 'itri ',itri, ist(itri,1), ist(itri,2), ist(itri,3) 
      return
!
      end
!============================================================
      function fact(n)
! acbauer - my implict none changes here caused problems, probably from real double precision integer things
!
!.....fact calculates the factorial of n...an integer
      ans = 1.
      if (n.le.1) go to 100
      do i=1,n
        ans = ans*i
      enddo
 100  fact = ans
      return
      end
!============================================================
      subroutine tridef
      use p_data
      use t_data
      use basic
      use nintegrate_mod

      implicit none

      integer itri, j, k, q, numelms
      real ti(20,20),rot(18,18), sum

      call numfac(numelms)
!      
!     start the loop over triangles within a rectangular region
      do 510 itri=1,numelms
!
!
!     define a,b,c and theta
! acbauer -- might want to take this out later 
      call getelmparams(itri,atri(itri),btri(itri),
     &        ctri(itri),ttri(itri))
!      
!
!   define the Inverse Transformation Matrix that enforces the 
!   condition that the normal slope between triangles has only cubic variation
!
      call tmatrix(ti,20,atri(itri),btri(itri),ctri(itri))
!
!
!     calculate the rotation matrix rot
      call rotation(rot,18,ttri(itri))
!
!     form the matrix g by multiplying ti and rot
      do k=1,20
         do j=1,18
            sum = 0
            do q = 1,18
               sum = sum + ti(k,q)*rot(q,j)
            enddo
            gtri(k,j,itri) = sum
         enddo
      enddo
!
      call offset(itri)

510      continue
!.....end of loop over itri
!
!
      return
      end







