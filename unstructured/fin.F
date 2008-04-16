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

      call finalizesolvers
      call PetscFinalize(ier)
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

      implicit none

      if(myrank.eq.0 .and. iprint.eq.1) print *, " entering input.."
!
!.....linear parameter  0-nonlinear,  1-linear
      linear = 0
!
!.....density advance parameter: 0 no advance,  1-advance density
      idens = 1
!
!.....electron pressure advance parameter: 0 no advance, 1-advance pe
      ipres = 0
!
!.....use current density equation instead of poloidal flux equation
      jadv = 0
!
!.....restart parameter 0-no restart, 1 restart
      irestart = 0
!
!.....associated restart parameter: 0-normal, 1-start time from zero
      istart = 0
!
!.....maximum KE before solution is re-scaled in linear case
      max_ke = 1.
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
!.....do timestep splitting
      isplitstep = 1
!
!.....iterate field solve
      iteratephi = 0
!
!.....iterate field solve
      irecalc_eta = 0
!
!.....modifications to implicitization
      imp_mod = 1
      igauge = 0
!
!.....form of toroidal velocity
      ivform = 0
!
!.....number of grad-shafranov iterations
      igs = 80
!
!.....CN implicitness parameter (thimp > 0.5 for stability)
      thimp = 0.5
      thimp_ohm = -1.
!
!.....hyper-resistivity coefficients
      deex = 1.
      hyper = 0.0
      hyperi= 0.0
      hyperv= 0.0
      hyperp= 0.0
!
!.....scale hyper-resistivity with dx**ihypdx
      ihypdx = 2
!
!.....scale hyper-resistivity with eta
      ihypeta = 1
!
!.....don't evolve resistivity
      iconst_eta = 0
!
!.....scale kappar with kappa
      ikapscale = 0
!
!.....add upwinding-like diffusion to pressure and density advance
      iupwind = 0
!
!.....viscosity, resistivity, and heat conduction
      amu = 0.0
      amupar = 0.0
      etar = 0.0
      eta0 = 0.0
      kappat = 0.0
      kappar = 0.0
      kappa0 = 0.0
      kappah = 0.0
      kappax = 0.0
!
!.....ratio of specific heats
      gam = 5./3.
!
!.....loop voltage
      vloop = 0.
!
!......toroidal magnetic field
      bzero = 1.
!
!.....initial toroidal velocity
      vzero = 0.
!
!.....boundary conditions
!     com_bc: 1 = sets laplacian(chi) = 0 on boundary
      com_bc = 1
      iconst_t = 1
      inograd_t = 0
      inonormalflow = 1
      inoslip_pol = 0
      inoslip_tor = 1
      inostress_tor = 0
      inocurrent_pol = 1
      inocurrent_tor = 1
!
!.....density source
      ipellet = 0
      pellet_x = xmag
      pellet_z = zmag
      pellet_rate = denm 
      pellet_var = 1.
!
      ionization = 0
      ionization_rate = denm
      ionization_temp = 0.01
      ionization_depth = 0.01
!
      nosig = 0
!
!.....ion skin depth
      db = 0.0
!
!     masking switch (0 no mask,  1 mask hyper terms at boundary
      imask = 0
!
      amu_edge = 0.
!
      isources = 0
!
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
      idevice = 0
!
!.....parameters needed for toroidal plasma
      tcuro = 1.
      xmag = xzero + 1.
      zmag = 0.
      xlim = xzero
      zlim = 0.
!
!..... toroidal current controller parameters
      tcur = tcuro
      control_p = 0.
      control_i = 0.
      control_d = 0.
!
!.....pressure = p0*(1 + p1*psi + p2*psi**2 + ...)
!                where psi is the normalized poloidal flux
      p0 = 0.01      
      p1 = -1.
      p2 = -2.
      pedge = 0.
      expn = 0.
      q0 = 1.
      th_gs = 0.5
!
!.....divertor information
      divertors = 0
      xdiv = xmag
      zdiv = zzero-1.
      divcur = 0.1
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
      gravr = 0.
      gravz = 0.
!
!.....1 = don't advance velocity
      istatic = 0
      iestatic = 0
!
!.....1 = conserve toroidal flux
      iconstflux = 0
!
!.....Toroidal mode number
      ntor = 0
!
!.....1 = output/input global restart files
      iglobalout = 0
      iglobalin = 0
!
!.....1 = adapt mesh upon initialization
      iadapt = 0
!
!
      open(5,file='C1input',form='formatted',status='old')
      read(5,nml=inputnl)
      if(amuc.eq.0.) amuc = amu
      if(thimp_ohm.lt.0) thimp_ohm = thimp
!      
      if(myrank.eq.0) write(*,nml=inputnl)

      if(myrank.eq.0 .and. iprint.eq.1) print *, " exiting input.."

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
      use arrays

      implicit none

      integer :: itri, i, j, k, ii, iii, numelms
      integer :: ibegin, iendplusone
      integer :: ibegin2, iendplusone2
      integer :: ibeginn, iendplusonen
      integer, dimension(4) :: nodeids
      real :: ti(20,20),rot(18,18), sum
      real :: f

      if(myrank.eq.1 .and. iprint.ge.1) print *, " Entering tridef..."

      call numfac(numelms)
!      
!     start the loop over triangles within a rectangular region
      do itri=1,numelms

         call nodfac(itri, nodeids)
         ist(itri,1) = nodeids(1)-1
         ist(itri,2) = nodeids(2)-1
         ist(itri,3) = nodeids(3)-1    
!
!        define a,b,c and theta
         call getelmparams(itri,atri(itri),btri(itri),
     &        ctri(itri),ttri(itri))  

!        define the Inverse Transformation Matrix that enforces the 
!        condition that the normal slope between triangles has only cubic variation
         call tmatrix(ti,20,atri(itri),btri(itri),ctri(itri))

!        calculate the rotation matrix rot
         call rotation(rot,18,ttri(itri))
!
!        form the matrix g by multiplying ti and rot
         do k=1,20
            do j=1,18
               sum = 0
               do ii = 1,18
                  sum = sum + ti(k,ii)*rot(ii,j)
               enddo
               gtri(k,j,itri) = sum
            enddo
         enddo

!        calculate fint (used for analytic integrations)
         do i=-6,maxi
            do j=-6,maxi
               fint(i,j) = 0
               if(i.ge.0 .and.j.ge.0)                                          &
     &              fint(i,j) = f(i,j,atri(itri),btri(itri),ctri(itri))
            enddo
         enddo

!        populate isval* arrays
         do iii=1,3
            call entdofs(1, ist(itri,iii)+1, 0, ibegin,  iendplusone )
            call entdofs(2, ist(itri,iii)+1, 0, ibegin2, iendplusone2)
            call entdofs(vecsize_vel,  ist(itri,iii)+1, 0,                  &
     &           ibeginn, iendplusonen)
            do ii=1,6
               i = (iii-1)*6 + ii
               isval1(itri,i) = ibegin +ii-1 ! 6*ist(itri,iii)+ii
               isval2(itri,i) = ibegin2+ii-1 ! 12*ist(itri,iii)+ii
               isvaln(itri,i) = ibeginn+ii-1 ! numvar*6*ist(itri,iii)+ii
            enddo
         enddo

      end do
!
      if(myrank.eq.1 .and. iprint.ge.1) print *, " Exiting tridef."

      return
      end







