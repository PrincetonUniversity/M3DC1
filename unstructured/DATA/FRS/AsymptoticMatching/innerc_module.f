      module innerc_module
      contains
      subroutine innerc(q,e,f,g,h,k,de,do,ifail,epsd,rmatch,iplot,yo,ye,xary)
c**************************************************************
c
c      author: s.c. jardin, ppl,march,82.
c      modified:   complex q  ray grimm, 30/3/82.
c      modified:  corrected complex q, leqt3f  (sj)   17/5/82
c
c
c      q is growth rate
c      e,f,g,h,k are surface quantities
c      de,do are ratios(returned) for even psi and odd psi solutions
c      ifail is returned 0 if a convergent answer is obtained
c      epsd is error tolorance (default is 1.e-4)
c
      real k
      integer iplot
      optional iplot,yo,ye,xary
c++++++++++++++++++
      complex de,do,dey,doy,hv,gv,gp,gpp,hsq,q,qsq,sq,
     1        sum,sum1,sum2,sum3,suma,sumb,a,b,c,dd2,dd3,
     2        t1,t2,wk,ev,vec,psia,psiap,aa,cc,dd,
     3        bquad,cquad,disc,s,sqdisc,denom,
     4        term1,term2,term3,xv,ratio,psim,sim,upsm,
     5        psi,si,ups,biglc,biglic,ex1,ex2,ex3,yo,ye
          
c+++++++++++++++++++++++++++rg 30/3/82
      dimension psi(1610),si(1610),ups(1610),
     1  a(3,3,1610),b(3,3,1610),c(3,3,1610),ev(3,3,1610),
     2  t1(3,3),t2(3,3),xary(1610),yo(1610,3),ye(1610,3)
      complex sc1,sc2,ai,sc1i,sc2i,cfac1,cfac2
      complex tl1c1,tl2c1,tl3c1,tp1c1,tp2c1,tp3c1
      complex tl1c2,tl2c2,tl3c2,tp1c2,tp2c2,tp3c2
c
      dimension hv(3,3),gv(3,3),dd2(3,3),dd3(3,3),dd3i(3,3),wk(9),
     1          d4(3,3),aa(3,3),psia(4,3),psiap(4,3),xv(3),
     2          amp(202),nas(2),vec(3,202),cc(3,3),dd(3,3),
     3          gp(3,3),gpp(3,3),hsq(3,3)
c
      dimension sol(2,6),sola(2,6),bci(3),doy(10),dey(10)
c
      de = 0.
      do = 0.
c
c.....choose matching point
c      rmatch = 3.
      if(epsd.le.0) epsd = .0001
      z1 = 4.*abs(g+k*f)
      z2 = 4.*abs(g-k*e)
      z3 = 4.*abs(k*h)
      zmax = rmatch**2
      if(z1.gt.zmax) zmax = z1
      if(z2.gt.zmax) zmax = z2
      if(z3.gt.zmax) zmax = z3
      qa = cabs(q)
      if(qa.lt.1.) bigl = (qa**0.25)*rmatch
      if(qa.ge.1.) bigl = qa*(zmax**0.5)
      biglsv = bigl
      imatch = 0
      go to 5
    4 bigl = 2*bigl
      imatch = imatch + 1
      if(imatch.gt.5) go to 3
    5 continue
      bigli = 1./bigl
      biglic = cmplx(bigli,0.)
      biglc = cmplx(bigl,0.)
c
c.....define constants
      nmaxa = 20
      gmke = g-k*e
      gpkf = g+k*f
      qsq = q*q
      ds = e+f+h
c
c.....start loop on step size
      inval = 0
    1 inval = inval + 1
      nval = 25*2**(inval-1)
      ifail = 0
      n = nval
c
      nsave = (n/2)*2 + 1
      n = nsave
c
c
c
c
      nm1 = n - 1
      del = bigl/(n-1)
      del2 = del*del
      del3 = del2*del
      del4 = del2*del2
c
      hv(1,1) = 0.
      hv(1,2) = 0.
      hv(1,3) = h
      hv(2,1) = -h/qsq
      hv(2,2) = 0.
      hv(2,3) = 0.
      hv(3,1) = h*k*q
      hv(3,2) = 0.
      hv(3,3) = 0.
c
      gv(1,1) = q
      gv(1,3) = 0.
      gv(2,3) = -(e+f)/qsq
      gv(3,2) = -(g-k*e)*q
c
      gp(1,1) = 0.
      gp(1,2) = -q
      gp(1,3) = 0.
      gp(2,1) = -1./q
      gp(2,3) = 0.
      gp(3,1) = -1./q
      gp(3,2) = 0.
c
      gpp(1,1) = 0.
      gpp(1,2) = 0.
      gpp(1,3) = 0.
      gpp(2,1) = 0.
      gpp(2,2) = 2./q
      gpp(2,3) = 0.
      gpp(3,1) = 0.
      gpp(3,2) = 0.
      gpp(3,3) = 2./q
c
      hsq(1,1) = q*k*h**2
      hsq(1,2) = 0.
      hsq(1,3) = 0.
      hsq(2,1) = 0.
      hsq(2,2) = 0.
      hsq(2,3) = -h**2/qsq
      hsq(3,1) = 0.
      hsq(3,2) = 0.
      hsq(3,3) = q*k*h**2
c
      do 302 ioddx=1,2
      iodd = ioddx-1
      do 99 l=1,n
      x = (l-1)*del
      xary(l) = x
c
      gp(3,3) = 2.*x/q
      gp(2,2) = 2.*x/q
      gv(3,3) = (g+k*f)*q + x**2/q
      gv(3,1) = -x/q
      gv(2,1) = -x/q
      gv(2,2) = x**2/q
      gv(1,2) = -q*x
c
      if(l.eq.n) go to 98
c
      do 21 i=1,3
      do 22 j=1,3
      sum1 = 0.
      do 23 kk=1,3
   23 sum1 = sum1 + gv(i,kk)*hv(kk,j) - hv(i,kk)*gv(kk,j)
     1            -hsq(i,kk)*hv(kk,j)
c
      a(i,j,l) = -.5*del*hv(i,j) - (del3/24.)*(2.*gp(i,j)+sum1)
c
      sum2 = 0.
      do 24 kk=1,3
   24 sum2 = sum2 - hv(i,kk)*gp(kk,j) - hsq(i,kk)*gv(kk,j)
     1            + gv(i,kk)*gv(kk,j)
c
      b(i,j,l) = -del2*gv(i,j) - (del4/12.)*(gpp(i,j)+sum2)
c
      c(i,j,l) = .5*del*hv(i,j) + (del3/24.)*(2.*gp(i,j)+sum1)
   22 continue
      a(i,i,l) = a(i,i,l) + 1.
      b(i,i,l) = b(i,i,l) - 2.
      c(i,i,l) = c(i,i,l) + 1.
   21 continue
      go to 99
   98 continue
c
      do 11 i=1,3
      do 10 j=1,3
      sum = 0.
      do 12 kk=1,3
      sum = sum + gv(i,kk)*hv(kk,j)
     1          + hv(i,kk)*gv(kk,j)
     2          + hv(i,kk)*hsq(kk,j)
   12 continue
      dd3(i,j) = -.5*del*hv(i,j)
     1         + (1./6.)*del2*(hsq(i,j)+gv(i,j))
     2         - (1./24.)*del3*(2.*gp(i,j)+sum)
      if(i.eq.j) dd3(i,j) = dd3(i,j) + 1.
      sum2 = 0.
      sum3 = 0.
      do 13 kk=1,3
      sum2 = sum2 + hv(i,kk)*gv(kk,j)
      sum3 = sum3 + hv(i,kk)*gp(kk,j)
     1            + hsq(i,kk)*gv(kk,j)
     2            + gv(i,kk)*gv(kk,j)
   13 continue
      dd2(i,j) = .5*del*gv(i,j) - (1./6.)*del2*(gp(i,j)+sum2)
     1         + (1./24.)*del3*(gpp(i,j)+sum3)
      if(i.eq.j) dd2(i,j) = dd2(i,j) + 1./del
   10 continue
   11 continue
c
   99 continue
c
c
c     solve coupled second order difference equations
c
c
c     left hand boundary condition
c
      do 51 i=1,3
      do 51 j=1,3
   51 t1(i,j) = -b(i,j,1)
      if(iodd.eq.0) go to 20
c
c     iodd=1 . . . even solution
c
      do 52 i=1,3
      t2(i,1) = a(i,1,1)+c(i,1,1)
      t2(i,2) = a(i,2,1)-c(i,2,1)
   52 t2(i,3) = a(i,3,1)-c(i,3,1)
      go to 30
   20 continue
c
c     iodd=0 . . . odd solution
c
      do 53 i=1,3
      t2(i,1) = a(i,1,1)-c(i,1,1)
      t2(i,2) = a(i,2,1)+c(i,2,1)
   53 t2(i,3) = a(i,3,1)+c(i,3,1)
   30 continue
c
      call leqt1f(t1,3,3,3,t2,0,wk,ier)
      do 54 i=1,3
      do 54 j=1,3
   54 ev(i,j,1) = t2(i,j)
      do 32 l=2,nm1
      do 34 i=1,3
      do 34 j=1,3
      sum1 = 0
      do 33 kk=1,3
   33 sum1 = sum1 + c(i,kk,l)*ev(kk,j,l-1)
      t1(i,j) = b(i,j,l) + sum1
   34 t2(i,j) = a(i,j,l)
      call leqt3fc(t1,t2,3)
      do 35 i=1,3
      do 35 j=1,3
   35 ev(i,j,l) = -t2(i,j)
   32 continue
c
      do 43 i=1,3
      do 42 j=1,3
      dd2(i,j) = dd2(i,j) - ev(i,j,n-1)/del
   42 continue
   43 continue
      call leqt1f(dd3,3,3,3,dd2,0,wk,ier)
c
c.....power-like asymtpotic solutions
      if(ioddx.eq.2) go to 301
      if(inval.gt.1) go to 301
      do 300 l=1,2
      s = -.5
      term = .5*sqrt(1.-4.*ds)
      if(l.eq.1) s = s+term
      if(l.eq.2) s = s-term
      do 320 i=1,3
      vec(i,1) = 0
      vec(i,2) = 1.
  320 continue
      psia(l,1) = biglc**(s+1)
      psia(l,2) = biglc**s
      nn = 0
      psia(l,3) = biglc**s
      psiap(l,1) = (s+1)*biglc**s
      psiap(l,2) = s*biglc**(s-1)
      psiap(l,3) = s*biglc**(s-1)
c
      do 401 nn=1,nmaxa
      sm2 = s-2*(nn-1)-2
      sm1 = s-2*(nn-1)-1
      sp1 = s-2*(nn-1)+1
      ss  = s-2*(nn-1)
      sp2 = s-2*(nn-1)+2
      ab = (h+sm2)*sm1
      bb = e+f-h*sm2
      bbmab = -bb-ab
      cc(1,1) = -qsq*bb*k*h*sp1/bbmab
      cc(1,2) = qsq*(ss*sm1 + bb*gmke)/bbmab
      cc(1,3) = -qsq*bb*gpkf/bbmab
      cc(2,1) = h*sp1*(1./q - qsq*bb*k/bbmab)
      cc(2,2) = qsq*(ss*sm1 + bb*gmke)/bbmab
      cc(2,3) = (e+f)/q - bb*qsq*gpkf/bbmab
      cc(3,1) = qsq*ab*k*h*sp1/bbmab
      cc(3,2) = qsq*(ss*sm1 - ab*gmke)/bbmab
      cc(3,3) = qsq*ab*gpkf/bbmab
c
      dd(1,1) = 0
      dd(1,2) = 0
      dd(1,3) = q*bb*sp2*sp1/bbmab
      dd(2,1) = 0
      dd(2,2) = q*sp2*sp1
      dd(2,3) = q*bb*sp2*sp1/bbmab
      dd(3,1) = 0
      dd(3,2) = 0
      dd(3,3) =-q*ab*sp2*sp1/bbmab
c
      do 310 i=1,3
      suma = 0
      sumb = 0
      do 305 j=1,3
      suma = suma + cc(i,j)*vec(j,nn+1)
      sumb = sumb + dd(i,j)*vec(j,nn)
  305 continue
      vec(i,nn+2) = suma + sumb
  310 continue
  401 continue
      nasy = nmaxa
      ampmin = 1.e20
      do 402 nn=1,nmaxa
      term1 = cabs(vec(1,nn+2)*biglc**(s-2*nn+1))
      term2 = cabs(vec(2,nn+2)*biglc**(s-2*nn))
      term3 = cabs(vec(3,nn+2)*biglc**(s-2*nn))
      amp(nn) = term1+term2+term3
      if(amp(nn) .gt. ampmin) go to 402
      ampmin = amp(nn)
      nasy = nn
      nas(l) = nn
  402 continue
      if(ampmin .gt. 0.1*epsd) go to 4
      do 400 nn=1,nasy
      psia(l,1) = psia(l,1) + vec(1,nn+2)*biglc**(s-2*nn+1)
      psia(l,2) = psia(l,2) + vec(2,nn+2)*biglc**(s-2*nn)
      psia(l,3) = psia(l,3) + vec(3,nn+2)*biglc**(s-2*nn)
      psiap(l,1) = psiap(l,1)+vec(1,nn+2)*(s-2*nn+1)*biglc**(s-2*nn)
      psiap(l,2) = psiap(l,2)+vec(2,nn+2)*(s-2*nn)*biglc**(s-2*nn-1)
      psiap(l,3) = psiap(l,3)+vec(3,nn+2)*(s-2*nn)*biglc**(s-2*nn-1)
  400 continue
  300 continue
c
c.....small exponential asymptotic solutions
      dr = e+f+h**2
      sq = csqrt(q)
      bquad = q*sq*(gmke+1.+k*dr)
      cquad = q**3*gpkf + (e+f) - dr*gmke
      disc = bquad**2 - 4.*cquad
c     irank = 6
c     if(disc.lt.0) go to 700
c     if(disc.gt.0) go to 451
c     irank = 5
c     sqdisc = 0
c     go to 452
  451 continue
      sqdisc = .25*csqrt(disc)
  452 continue
c
      do 600 l=3,4
      s = -.5-.25*bquad
      if(l.eq.3) s = s+sqdisc
      if(l.eq.4) s = s-sqdisc
      denom = dr/q - h*sq
      term1 = (-ds*q-2.*s*h*q)/denom
      term2 = 1.
      term3 = (qsq+sq*(2*s+1-h))/denom
c
c     factor out biglc**s*exp(-biglc**2/(2*sq))
      psia(l,1) = term1*bigli
      psia(l,2) = term2
      psia(l,3) = term3
      psiap(l,1) = term1*((s-1.)*biglic**2-1./sq)
      psiap(l,2) = term2*(s*bigli-bigl/sq)
      psiap(l,3) = term3*(s*bigli-bigl/sq)
  600 continue
      go to 701
  700 continue
      ai = cmplx(0.,1.)
      sqdisc = .25*csqrt(-disc)
      sreal = -.5-.25*bquad
      sc1i = cmplx(0,aimag(sqdisc))
      sc1 = cmplx(sreal,aimag(sqdisc))
      sc2i = cmplx(0,aimag(-sqdisc))
      sc2 = cmplx(sreal,aimag(-sqdisc))
      cfac1 = biglc**sc1i
      cfac2 = biglc**sc2i
      denom = dr/q-h*sq
      tl1c1 = ((-ds*q-2.*sc1*h*q)/denom)*cfac1*bigli
      tl1c2 = ((-ds*q-2.*sc2*h*q)/denom)*cfac2*bigli
      tl2c1 = cfac1
      tl2c2 = cfac2
      tl3c1 = ((qsq+sq*(2.*sc1+1-h))/denom)*cfac1
      tl3c2 = ((qsq*sq*(2.*sc2+1-h))/denom)*cfac2
      tp1c1 = tl1c1*((sc1-1.)*bigli - bigl/sq)
      tp1c2 = tl1c2*((sc2-1.)*bigli - bigl/sq)
      tp2c1 = tl2c1*(sc1*bigli - bigl/sq)
      tp2c2 = tl2c2*(sc2*bigli - bigl/sq)
      tp3c1 = tl3c1*(sc1*bigli - bigl/sq)
      tp3c2 = tl3c2*(sc2*bigli - bigl/sq)
      psia(3,1) = real(tl1c1)
      psia(3,2) = real(tl2c1)
      psia(3,3) = real(tl3c1)
      psia(4,1) = aimag(tl1c1)
      psia(4,2) = aimag(tl2c1)
      psia(4,3) = aimag(tl3c1)
      psiap(3,1) = real(tp1c1)
      psiap(3,2) = real(tp2c1)
      psiap(3,3) = real(tp3c1)
      psiap(4,1) = aimag(tp1c1)
      psiap(4,2) = aimag(tp2c1)
      psiap(4,3) = aimag(tp3c1)
  701 continue
  301 continue
c
c.....apply matching condition
      do 45 i=1,3
      do 40 j=1,3
      sum = 0.
      do 47 kk=1,3
   47 sum = sum + dd2(i,kk)*psia(j+1,kk)
      aa(i,j) = -psiap(j+1,i) + sum
   40 continue
      sum = 0
      do 41 kk=1,3
   41 sum = sum + dd2(i,kk)*psia(1,kk)
      xv(i) = psiap(1,i) - sum
   45 continue
      ier = 0
      call leqt1f(aa,1,3,3,xv,0,wk,ier)
      ratio = 1./xv(1)
      psim = psia(1,1)+xv(1)*psia(2,1)+xv(2)*psia(3,1)+xv(3)*psia(4,1)
      sim  = psia(1,2)+xv(1)*psia(2,2)+xv(2)*psia(3,2)+xv(3)*psia(4,2)
      upsm = psia(1,3)+xv(1)*psia(2,3)+xv(2)*psia(3,3)+xv(3)*psia(4,3)
c
c
      if(ioddx.eq.1) do = ratio
      if(ioddx.eq.2) de = ratio

      if(iplot.eq.0) go to 603
      go to(601,602),ioddx
 601   continue
      yo(n,1) = psim*ratio
      yo(n,2) = sim*ratio
      yo(n,3) = upsm*ratio
      do 48 ll=2,n
      l = n+1-ll
      do 49 i=1,3
      sum = 0
      do 46 j=1,3
 46       sum = sum + ev(i,j,l)*yo(l+1,j)
 49        yo(l,i) = sum
 48         continue
      go to 603
 602   continue
      ye(n,1) = psim*ratio
      ye(n,2) = sim*ratio
      ye(n,3) = upsm*ratio
      do 68 ll=2,n
      l = n+1-ll
      do 67 i=1,3
      sum = 0
      do 66 j=1,3
 66       sum = sum + ev(i,j,l)*ye(l+1,j)
 67        ye(l,i) = sum
 68         continue
 603         continue


  302 continue
      doy(inval) = do
      dey(inval) = de
c     write(59,1113) n,bigl,ampmin,de,do
 1113 format(i5,4e16.8)
c
c.....check for convergence
      if(inval.eq.1) go to 1
      if(cabs(doy(inval)-doy(inval-1)).gt.
     1   epsd*cabs(doy(inval))) ifail = 1
      if(cabs(dey(inval)-dey(inval-1)).gt.
     1   epsd*cabs(dey(inval))) ifail = ifail+2
      if(ifail.gt.0 .and. inval.le.6) go to 1
c
c      apply richardson's extrapolation
      do = (16.*doy(inval)-doy(inval-1))/15.
      de = (16.*dey(inval)-dey(inval-1))/15.
c     define boundary conditions at right endpoint
c
c      psi(n) = psim
c      si(n) = sim
c      ups(n) = upsm
c
c     step solution from right to left
c
    3 continue
c$$$      open(unit=10,file='xary.txt')
c$$$      do 62 i=1,1610
c$$$         write(10,*) xary(i)
c$$$ 62   continue
c$$$      open(unit=11,file='si.txt')
c$$$      do 63 i=1,1610
c$$$         write(11,*) si(i)
c$$$ 63      continue
c$$$         stop
      return
      end subroutine
      subroutine leqt1f (a,m,n,ia,b,idgt,wkarea,ier)
c
c-leqt1f--------s-------library 3---------------------------------------
c
c   function            - linear equation solution - full storage
c                           mode - space economizer solution.
c   usage               - call leqt1f (a,m,n,ia,b,idgt,wkarea,ier)
c   parameters   a      - input matrix of dimension n by n containing
c                           the coefficient matrix of the equation
c                           ax = b.
c                         on output, a is replaced by the lu
c                           decomposition of a rowwise permutation of
c                           a.
c                m      - number of right-hand sides.(input)
c                n      - order of a and number of rows in b.(input)
c                ia     - number of rows in the dimension statement
c                           for a and b in the calling program. (input)
c                b      - input matrix of dimension n by m containing
c                           right-hand sides of the equation ax = b.
c                         on output, the n by m solution x replaces b.
c                idgt   - input option.
c                         if idgt is greater than 0, the elements of
c                           a and b are assumed to be correct to idgt
c                           decimal digits and the routine performs
c                           an accuracy test.
c                         if idgt equals zero, the accuracy test is
c                           bypassed.
c                wkarea - work area of dimension greater than or equal
c                           to n.
c                ier    - error parameter
c                         terminal error = 128+n.
c                           n = 1 indicates that a is algorithmically
c                             singular. (see the chapter l prelude).
c                         warning error = 32+n.
c                           n = 2 indicates that the accuracy test
c                                 failed.
c                                 the computed solution may be in error
c                                 by more than can be accounted for by
c                                 the uncertainty of the data.
c                                 this warning can be produced only if
c                                 idgt is greater than 0 on input.
c                                 see chapter l prelude for further
c                                 discussion.
c   precision           - single
c   reqd. imsl routines - ludatf,luelmf,uertst
c   language            - fortran
c-----------------------------------------------------------------------
c   latest revision     - august 15, 1973
c
      implicit complex (a-h,o-z)
      dimension          a(ia,1),b(ia,1),wkarea(1)
c                                  initialize ier
      ier=0
c                                  decompose a
      call ludatf (a,a,n,ia,idgt,d1,d2,wkarea,wkarea,wa,ier)
      if (ier .gt. 128) go to 9000
c                                  call routine luelmf (forward and
c                                  backward substitutions)
      do 10 j=1,m
         call luelmf (a,b(1,j),wkarea,n,ia,b(1,j))
   10 continue
      if (ier .eq. 0) go to 9005
 9000 continue
      call uertst (ier,6hleqt1f)
 9005 return
      end subroutine
      subroutine ludatf (a,lu,n,ia,idgt,d1,d2,ipvt,equil,wa,ier)
c
c-ludatf--------s-------library 3---------------------------------------
c
c   function            - l-u decomposition by the crout algorithm
c                           with optional accuracy test.
c   usage               - call ludatf(a,lu,n,ia,idgt,d1,d2,ipvt,
c                                     equil,wa,ier)
c   parameters   a      - input matrix of dimension n by n containing
c                           the matrix to be decomposed
c                lu     - real output matrix of dimension n by n
c                           containing the l-u decomposition of a
c                           rowwise permutation of the input matrix.
c                           for a description of the format of lu, see
c                           example.
c                n      - input scalar containing the order of the
c                           matrix a.
c                ia     - input scalar containing the row dimension of
c                           matrices a and lu in the calling program.
c                idgt   - input option.
c                           if idgt is greater than zero, the non-zero
c                           elements of a are assumed to be correct to
c                           idgt decimal places.  ludatf performs an
c                           accuracy test to determine if the computed
c                           decomposition is the exact decomposition
c                           of a matrix which differs from the given one
c                           by less than its uncertainty.
c                         if idgt is equal to zero, the accuracy test is
c                           bypassed.
c                d1     - output scalar containing one of the two
c                           components of the determinant. see
c                           description of parameter d2, below.
c                d2     - output scalar containing one of the
c                           two components of the determinant. the
c                           determinant may be evaluated as (d1)(2**d2)
c                ipvt   - output vector of length n containing the
c                           permutation indices. see document
c                           (algorithm).
c                equil  - output vector of length n containing
c                           reciprocals of the absolute values of
c                           the largest (in absolute value) element
c                           in each row.
c                wa     - accuracy test parameter, output only if
c                           idgt is greater than zero.
c                           see element documentation for details.
c                ier    - error parameter
c                         terminal error=128+n
c                           n = 1 indicates that matrix a is
c                                 algorithmically singular. (see the
c                                 chapter l prelude).
c                         warning error=32+n
c                           n = 2 indicates that the accuracy test
c                                 failed.
c                                 the computed solution may be in error
c                                 by more than can be accounted for by
c                                 the uncertainty of the data.
c                                 this warning can be produced only if
c                                 idgt is greater than 0 on input.
c                                 see chapter l prelude for further
c                                 discussion.
c   precision           - single
c   reqd. imsl routines - uertst
c   language            - fortran
c-----------------------------------------------------------------------
c   latest revision     - august 15, 1973
c
      implicit complex (a-h,o-z)
      dimension          a(ia,1),lu(ia,1),ipvt(1),equil(1)
      complex               lu,ipvt,equil
      data               zero,one,four,sixtn,sixth/0.0,1.,4.,16.,.0625/
c                                  initialization
      ier = 0
      rn = n
      wrel = 0.0
      d1 = 1.0
      d2 = 0.0
      biga = 0.0
      do 10 i=1,n
         big = 0.0
         do 5 j=1,n
            p = a(i,j)
            lu(i,j) = p
            p = cabs(p)
            if (cabs(p) .gt. cabs(big)) big = p
    5    continue
         if (cabs(big) .gt. cabs(biga)) biga = big
         if (cabs(big) .eq. 0.0) go to 110
         equil(i) = 1.0/big
   10 continue
      do 105 j=1,n
         jm1 = j-1
         if (jm1 .lt. 1) go to 40
c                                  compute u(i,j), i=1,...,j-1
         do 35 i=1,jm1
            sum = lu(i,j)
            im1 = i-1
            if (idgt .eq. 0) go to 25
c                                  with accuracy test
            ai = cabs(sum)
            wi = 0.0
            if (im1 .lt. 1) go to 20
            do 15 k=1,im1
               t = lu(i,k)*lu(k,j)
               sum = sum-t
               wi = wi+cabs(t)
   15       continue
            lu(i,j) = sum
   20       wi = wi+cabs(sum)
            if (cabs(ai) .eq. 0.0) ai = biga
            test = wi/ai
            if (cabs(test) .gt. cabs(wrel)) wrel = test
            go to 35
c                                  without accuracy
   25       if (im1 .lt. 1) go to 35
            do 30 k=1,im1
               sum = sum-lu(i,k)*lu(k,j)
   30       continue
            lu(i,j) = sum
   35    continue
   40    p = 0.0
c                                  compute u(j,j) and l(i,j), i=j+1,...,
         do 70 i=j,n
            sum = lu(i,j)
            if (idgt .eq. 0) go to 55
c                                  with accuracy test
            ai = cabs(sum)
            wi = 0.0
            if (jm1 .lt. 1) go to 50
            do 45 k=1,jm1
               t = lu(i,k)*lu(k,j)
               sum = sum-t
               wi = wi+cabs(t)
   45       continue
            lu(i,j) = sum
   50       wi = wi+cabs(sum)
            if (cabs(ai) .eq. 0.0) ai = biga
            test = wi/ai
            if (cabs(test) .gt. cabs(wrel)) wrel = test
            go to 65
c                                  without accuracy test
   55       if (jm1 .lt. 1) go to 65
            do 60 k=1,jm1
               sum = sum-lu(i,k)*lu(k,j)
   60       continue
            lu(i,j) = sum
   65       q = equil(i)*cabs(sum)
            if (cabs(p) .ge. cabs(q)) go to 70
            p = q
            imax = i
   70    continue
c                                  test for algorithmic singularity
         if (rn+p .eq. rn) go to 110
         if (j .eq. imax) go to 80
c                                  interchange rows j and imax
         d1 = -d1
         do 75 k=1,n
            p = lu(imax,k)
            lu(imax,k) = lu(j,k)
            lu(j,k) = p
   75    continue
         equil(imax) = equil(j)
   80    ipvt(j) = imax
         d1 = d1*lu(j,j)
   85    if (cabs(d1) .le. 1.0) go to 90
         d1 = d1*sixth
         d2 = d2+4.0
         go to 85
   90    if (cabs(d1) .ge. cabs(sixth)) go to 95
         d1 = d1*sixtn
         d2 = d2-4.0
         go to 90
   95    continue
         jp1 = j+1
         if (jp1 .gt. n) go to 105
c                                  divide by pivot element u(j,j)
         p = lu(j,j)
         do 100 i=jp1,n
            lu(i,j) = lu(i,j)/p
  100    continue
  105 continue
c                                  perform accuracy test
      if (idgt .eq. 0) go to 9005
      p = 3*n+3
      wa = p*wrel
      if (wa+10.0**(-idgt) .ne. wa) go to 9005
      ier = 34
      go to 9000
c                                  algorithmic singularity
  110 ier = 129
      d1 = 0.0
      d2 = 0.0
 9000 continue
c                                  print error
      call uertst(ier,6hludatf)
 9005 return
      end subroutine
      subroutine luelmf (a,b,ipvt,n,ia,x)
c
c-luelmf--------s-------library 3---------------------------------------
c
c   function            - elimination part of solution of ax=b -
c                           full storage mode
c   usage               - call luelmf (a,b,ipvt,n,ia,x)
c   parameters   a      - the result, lu, computed in the subroutine
c                           *ludatf*, where l is a lower triangular
c                           matrix with ones on the mai diagonal. u is
c                           upper triangular. l and u are stored as a
c                           single matrix a, and the unit diagonal of
c                           l is not stored
c                b      - b is a vector of length n on the right hand
c                           side of the equation ax=b
c                ipvt   - the permutation matrix returned from the
c                           subroutine *ludatf*, stored as an n length
c                           vector
c                n      - order of a and number of rows in b
c                ia     - number of rows in the dimension statement
c                           for a in the calling program.
c                x      - the result x
c   precision           - single
c   language            - fortran
c-----------------------------------------------------------------------
c   latest revision     - april 11,1975
c
      implicit complex (a-h,o-z)
      dimension          a(ia,1),b(1),ipvt(1),x(1)
      complex ipvt
c                                  solve ly = b for y
      do 5 i=1,n
    5 x(i) = b(i)
      iw = 0
      do 20 i=1,n
         ip = ipvt(i)
         sum = x(ip)
         x(ip) = x(i)
         if (iw .eq. 0) go to 15
         im1 = i-1
         do 10 j=iw,im1
            sum = sum-a(i,j)*x(j)
   10    continue
         go to 20
   15    if (sum .ne. 0.) iw = i
   20 x(i) = sum
c                                  solve ux = y for x
      do 30 ib=1,n
         i = n+1-ib
         ip1 = i+1
         sum = x(i)
         if (ip1 .gt. n) go to 30
         do 25 j=ip1,n
            sum = sum-a(i,j)*x(j)
   25   continue
   30 x(i) = sum/a(i,i)
      return
      end subroutine
      subroutine uertst (ier,name)
c
c-uertst----------------library 3---------------------------------------
c
c   function            - error message generation
c   usage               - call uertst(ier,name)
c   parameters   ier    - error parameter. type + n  where
c                           type= 128 implies terminal error
c                                  64 implies warning with fix
c                                  32 implies warning
c                           n   = error code relevant to calling routine
c                name   - input scalar containing the name of the
c                           calling routine as a 6-character literal
c                           string.
c   language            - fortran
c-----------------------------------------------------------------------
c   latest revision     - august 1, 1973
c
      character*6 name
      dimension          ityp(2,4),ibit(4)
      character*8 ityp
      integer            warn,warf,term,printr
      equivalence        (ibit(1),warn),(ibit(2),warf),(ibit(3),term)
      data     ityp       /8hwarning   ,8h          ,
     *                    8h warn(wi  ,8hth fix)   ,
     *                    8hterminal  ,8h          ,
     *                    8hnon-def   ,8hd         /,
     *         ibit      / 32,64,128,0/
      data               printr/59/
      ier2=ier
      if (ier2 .ge. warn) go to 5
c                                  non-defined
      ier1=4
      go to 20
   5  if (ier2 .lt. term) go to 10
c                                  terminal
      ier1=3
      go to 20
  10  if (ier2 .lt. warf) go to 15
c                                  warning(with fix)
      ier1=2
      go to 20
c                                  warning
  15  ier1=1
c                                  extract *n*
  20  ier2=ier2-ibit(ier1)
c                                  print error message
c     write (printr,25) (ityp(i,ier1),i=1,2),name,ier2,ier
   25 format(26h *** i m s l(uertst) ***  ,2a10,4x,a6,4x,i2,
     1   8h (ier = ,i3,1h))
      return
      end subroutine
      subroutine leqt3fc(c,b,n)
      implicit complex(a-h,o-z)
      dimension c(9),b(3,n)
c changed b(3,3) into b(3,n)
      t2 = c(1)*c(5)-c(4)*c(2)
      t4 = c(2)*c(6)-c(5)*c(3)
      t5 = c(1)*c(8)-c(7)*c(2)
      t6 = c(2)*c(6)-c(5)*c(3)
      t7 = c(2)*c(9)-c(8)*c(3)
      t8 = c(1)*c(5)-c(4)*c(2)
      do 10 i=1,n
      t1 = c(2)*b(3,i)-b(2,i)*c(3)
      t3 = c(1)*b(2,i)-b(1,i)*c(2)
      x3=(-t1*t2+t3*t4) / (t5*t6-t7*t8)
      x2=(-b(2,i)*c(3)+b(3,i)*c(2)-t7*x3)/t4
      b(1,i)=(-c(4)*x2-c(7)*x3+b(1,i))/c(1)
      b(2,i)=x2
      b(3,i)=x3
   10 continue
      return
      end subroutine

      end module innerc_module
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccc end of innerc subroutines cccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
