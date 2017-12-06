      module modtest
      contains
      subroutine innerc(q,e,f,g,h,k,de,do,ifail,epsd,rmatch)
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
c++++++++++++++++++
      complex de,do,dey,doy,hv,gv,gp,gpp,hsq,q,qsq,sq,
     1        sum,sum1,sum2,sum3,suma,sumb,a,b,c,dd2,dd3,
     2        t1,t2,wk,ev,vec,psia,psiap,aa,cc,dd,
     3        bquad,cquad,disc,s,sqdisc,denom,
     4        term1,term2,term3,xv,ratio,psim,sim,upsm,
     5        psi,si,ups,biglc,biglic,ex1,ex2,ex3
c+++++++++++++++++++++++++++rg 30/3/82
      dimension psi(1610),si(1610),ups(1610),
     1  a(3,3,1610),b(3,3,1610),c(3,3,1610),ev(3,3,1610),
     2  t1(3,3),t2(3,3),xary(1610)
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
      psi(n) = psim
      si(n) = sim
      ups(n) = upsm
c
c     step solution from right to left
c
    3 continue
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
c                           matrix with ones on the main diagonal. u is
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccc end of innerc subroutines cccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccc Array utilities cccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine linspace(xmin,xmax,x,n)
      real xmin,xmax
      real, dimension(n)::x
      integer  i,n
      
      do i=1,n
       x(i) = (xmax-xmin) * real(i-1) / real(n-1) + xmin
      end do
      end subroutine linspace

      subroutine linspaceint(xmin,xmax,x,n)
      integer xmin,xmax
      integer, dimension(n)::x
      integer  i,n
      
      do i=1,n
       x(i) = int((xmax-xmin) * real(i-1) / real(n-1) + xmin)
      end do
      end subroutine linspaceint

      subroutine logspace(xmin,xmax,x,n)
      integer n
      real xmin,xmax
      real, dimension(n)::x
      call linspace(log10(xmin),log10(xmax),x,n)
      x = 10.**x
      end subroutine logspace


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccc   subroutine equi  cccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine equi(size,r,p,dp,q,dq,Bz,Bp,k,i_system,alpha,q0,Bp0)
      implicit none
      integer size
      integer i,i_system
      INTEGER IER, PGBEG,ifail
      real k,q0,q_edge,alpha_q,p0,a,Bz0,step,er,Bz_edge,beta
      real p,q,dp,dq,integral1,integral2,Bz,Bp
      real Bp0
      real alpha
      real r
      dimension r(size),p(size),q(size),dp(size),dq(size),
     1integral1(size),integral2(size),Bz(size),Bp(size)
      
      a=1.
      step=a/real(size-1)
      call linspace(0.,1.,r,size)
c      r = (/((i*step),i=0,size-1)/)

c$$$ccccccc Tokamak ccccccccc   
c$$$      if(i_system==0)then
c$$$      k=1./3.
c$$$      q0=0.8
c$$$      q_edge=3.2
c$$$      alpha_q=q_edge/q0-1.
c$$$      
c$$$c     a should be kept to one for the time being
c$$$      beta=0.05
c$$$      p0=0.5*beta
c$$$c      Bz0=1.
c$$$      Bz_edge=1
c$$$      
c$$$      integral1=0.
c$$$      integral2=0.
c$$$     
c$$$      q=q0*(1+alpha_q*(r/a)**2)
c$$$      p=p0*(1-(r/a)**2)
c$$$      dq=q0*2.*alpha_q*r/a**2
c$$$      dp=p0*(-2.)*r/a**2
c$$$      ifail=0
c$$$      do 10 i=4,size
c$$$c      call D01GAF(r(1:i),(2.*r(1:i)-r(1:i)**2*dq(1:i)/q(1:i))/(q(1:i)**2
c$$$c     1/k**2+r(1:i)**2 ),i,integral1(i),er,ifail)
c$$$      call D01GAF(r(size-i+1:size),-(2.*r(size-i+1:size)-r(size-i+1:size
c$$$     1)**2*dq(size-i+1:size)/q(size-i+1:size))/(q(size-i+1:size)**2
c$$$     2/k**2+r(size-i+1:size)**2 ),i,integral1(size-i+1),er,ifail)
c$$$c      write(*,*) 'integral1 =', integral1(i)
c$$$c       write(*,*)'er =',er
c$$$10    continue
c$$$
c$$$      do 11 i=4,size
c$$$c      call D01GAF(r(1:i),exp(2.*integral1(1:i))*(-(1/k)**2*q(1:i)**2*
c$$$c     1dp(1:i))/0.5/(q(1:i)**2/k**2+r(1:i)**2),i,integral2(i),er,ifail)
c$$$      call D01GAF(r(size-i+1:size),-exp(2.*integral1(size-i+1:size))*(-
c$$$     1(1/k)**2*q(size-i+1:size)**2*
c$$$     2dp(size-i+1:size))/0.5/(q(size-i+1:size)**2/k**2+r(size-i+1:size)
c$$$     3**2),i,integral2(size-i+1),er,ifail)
c$$$c      write(*,*) 'integral2 =', integral2(i)
c$$$c        write(*,*)'ifail =',ifail
c$$$ 11   continue
c$$$      
c$$$c      Bz=exp(-integral1)*sqrt(Bz0**2+integral2)
c$$$      Bz=exp(-integral1)*sqrt(Bz_edge**2+integral2)
c$$$      Bp=Bz*r*k/q
c$$$
c$$$c      do 12 i=1,size
c$$$c      write(*,*) 'q=', q(i)
c$$$c      write(*,*) 'Bz2 =', Bz2(i)
c$$$c 12   continue
c$$$     
c$$$
c$$$c     
c$$$c.....initialize plotting
c$$$      call ncarcgm(1,'Bz.cgm')
c$$$c....make plot
c$$$      call mapg(r(1),r(size),0.,2.,.1,.5,.7,1.)
c$$$      call trace(r,Bz,size,-1,-1,0.,0.)
c$$$c      call setlch(xary(1)-1.0,yary(1)-0.1,0,2,1,-1)
c$$$c      write(s100,1002)
c$$$c      call gtext(s100,80,0)
c$$$c 1002 format(" test plot")
c$$$c
c$$$c
c$$$      call frame(0)
c$$$c.....finalize plot file
c$$$      call plote
c$$$
c$$$c.....initialize plotting
c$$$      call ncarcgm(1,'Bp.cgm')
c$$$c....make plot
c$$$      call mapg(r(1),r(size),0.,0.1,.1,.5,.7,1.)
c$$$      call trace(r,Bp,size,-1,-1,0.,0.)
c$$$c      call setlch(xary(1)-1.0,yary(1)-0.1,0,2,1,-1)
c$$$c      write(s100,1002)
c$$$c      call gtext(s100,80,0)
c$$$c 1002 format(" test plot")
c$$$c
c$$$c
c$$$      call frame(0)
c$$$c.....finalize plot file
c$$$      call plote
c$$$    
c$$$      end if

ccccc Spheromak (DeLucia)ccccccccc

      if(i_system==1)then
c      Bp0=1.
c      k=1./2.
c      q0=0.4
c      q_edge=3.2
c      alpha=0.7
      
c     a should be kept to one for the time being
c      beta=0.05
c      p0=0.5*beta

       
      integral1=0.
      integral2=0.
     
      q=q0*(1-(r/a)**2)
c      p=p0*(1-(r/a)**2)
      dq=-q0*2.*r/a**2
c      dp=p0*(-2.)*r/a**2
      ifail=0
c       write(*,*) 'bug'
      do 20 i=4,size
c          write(*,*) 'i=',i
c          write(*,*) 'r(i)=',r(i)

      call D01GAF(r(1:i),(2.*r(1:i)+q/k**2*dq
     1 -alpha/8./k**2*r(1:i)*dq**2 )/(q(1:i)**2/k**2
     2 +r(1:i)**2 ),i,integral1(i),er,ifail)
c      call D01GAF(r(size-i+1:size),-(2.*r(size-i+1:size)-r(size-i+1:size
c     1)**2*dq(size-i+1:size)/q(size-i+1:size))/(q(size-i+1:size)**2
c     2/k**2+r(size-i+1:size)**2 ),i,integral1(size-i+1),er,ifail)
c      write(*,*) 'integral1 =', integral1(i)
c       write(*,*)'er =',er
 20   continue
     
      Bp=Bp0*r*exp(-integral1)
      Bz=Bp0*q/k*exp(-integral1)
      
      do 21 i=4,size
      call D01GAF(r(1:i),r(1:i)*Bz(1:i)**2*(dq(1:i)/q(1:i))**2,
     1 i,integral2(i),er,ifail)
c      call D01GAF(r(size-i+1:size),-exp(2.*integral1(size-i+1:size))*(-
c     1(1/k)**2*q(size-i+1:size)**2*
c     2dp(size-i+1:size))/0.5/(q(size-i+1:size)**2/k**2+r(size-i+1:size)
c     3**2),i,integral2(size-i+1),er,ifail)
c      write(*,*) 'integral2 =', integral2(i)
c        write(*,*)'ifail =',ifail
 21   continue
      p=-alpha/8.*integral2
      dp=-alpha*r*Bz**2/8.*(dq/q)**2
      p(size)=p(size-1)
      write(*,*) 'Bzedge=',Bz(size)
      write(*,*) 'Bz0=',Bz(1)
      beta=abs(2.*p(size)/Bz(1)**2)
      write(*,*) 'beta=',beta
      
c      Bz=exp(-integral1)*sqrt(Bz0**2+integral2)
c      Bz=exp(-integral1)*sqrt(Bz_edge**2+integral2)
c      Bp=Bz*r*k/q
   

c.....initialize plotting
      call ncarcgm(1,'q.cgm')
c....make plot
      call mapg(r(1),r(size),minval(q),maxval(q),.1,.5,.7,1.)
      call trace(r,q,size,-1,-1,0.,0.)
c      call setlch(xary(1)-1.0,yary(1)-0.1,0,2,1,-1)
c      write(s100,1002)
c      call gtext(s100,80,0)
c 1002 format(" test plot")
c
c
      call frame(0)
c.....finalize plot file
      call plote

   
c.....initialize plotting
      call ncarcgm(1,'p.cgm')
c....make plot
      call mapg(r(1),r(size),minval(p),maxval(p),.1,.5,.7,1.)
      call trace(r,p,size,-1,-1,0.,0.)
c      call setlch(xary(1)-1.0,yary(1)-0.1,0,2,1,-1)
c      write(s100,1002)
c      call gtext(s100,80,0)
c 1002 format(" test plot")
c
c
      call frame(0)
c.....finalize plot file
      call plote



c.....initialize plotting
      call ncarcgm(1,'Bz.cgm')
c....make plot
      call mapg(r(1),r(size),0.,3.,.1,.5,.7,1.)
      call trace(r,Bz,size,-1,-1,0.,0.)
c      call setlch(xary(1)-1.0,yary(1)-0.1,0,2,1,-1)
c      write(s100,1002)
c      call gtext(s100,80,0)
c 1002 format(" test plot")
c
c
      call frame(0)
c.....finalize plot file
      call plote

c.....initialize plotting
      call ncarcgm(1,'Bp.cgm')
c....make plot
      call mapg(r(1),r(size),0.,3.,.1,.5,.7,1.)
      call trace(r,Bp,size,-1,-1,0.,0.)
c      call setlch(xary(1)-1.0,yary(1)-0.1,0,2,1,-1)
c      write(s100,1002)
c      call gtext(s100,80,0)
c 1002 format(" test plot")
c
c
      call frame(0)
c.....finalize plot file
      call plote

      end if

      
      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccc   subroutine rational_surf cccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rational_surf(m,n,size,r,p,dp,q,dq,Bz,Bp,k,E,F,G,H,K0
     1 ,ifail,dqc,Bpc,rc,ta,ic)
      implicit none
      integer m,n,size,ic,ifail
      real, dimension(size):: r,p,q,dp,dq,Bz,Bp
      real k,qc,rc,dpc,dqc,Bzc,Bpc,pc,E,F,G,H,K0,ta
      integer, dimension(1)::ic_array,val_array
      real tab
      qc=real(m)/real(n)
     

c      val_array=minval(abs(q-qc))
c      write(*,*) 'min =',val_array(1)
c      write(*,*) 'min =',minval(abs(q-qc))
c      write(*,*) 'qc =',qc
      if (minval(abs(q-qc)).ge.0.01)  then
         ifail=1
c         write(*,*) 'outside'
      return
      end if
      ic_array=minloc(abs(q-qc))
      ic=ic_array(1)
      write(*,*) 'ic =',ic
      rc=r(ic)
       write(*,*) 'rc =',rc
      Bpc=Bp(ic)
      Bzc=Bz(ic)
      ta=rc/sqrt(Bzc**2+Bpc**2)
c*sqrt(4*3.14159265)
      write(*,*) 'Bc =',sqrt(Bzc**2+Bpc**2)
      write(*,*) 'Bzc =',Bzc
      write(*,*) 'Bpc =',Bpc
      write(*,*) 'ta_b =',ta
      tab=rc/Bpc
      write(*,*) 'ta_bp=',tab
      write(*,*) 'ta_bz',rc/Bzc
      ta=rc/bzc*sqrt(4.*3.14159265)
      pc=p(ic)
      dpc=dp(ic)
      dqc=dq(ic)
      F=dpc**2*rc**2/(Bpc**4/k**2*dqc**2)
      E = -2.*dpc*rc/(Bpc**2/k**2*dqc**2) - F
      K0=1./F
      G = 3./5. * (Bzc**2+Bpc**2)/pc
c       G = 3./5. * Bpc**2/pc
      H=0.
      

      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccc outerfunc cccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine outerfunc(rin,vecin,vecout)
      implicit none
      real rin,rc
      real, dimension(2)::vecin,vecout
      real, dimension(10000)::rb,fr,gr,rbb,frb,grb,rbb2
      real, dimension(ninterp*(ninterp+1)/2)::Cf,Cg
      real, dimension(1)::Cfedge,Cgedge
      integer, dimension(1)::icb
      integer ic,last,n,np,ninterp
      common /fg_block/ fr,gr,rc,ic,rb
      common /size_interp/ ninterp
      n=ninterp
      frb=fr
      grb=gr
      rbb=rb
      rbb2=rb
c                           write(*,*) 'r(icbfirst)',rbb(2)
      icb=minloc(abs(rbb-rin))
      
c$$$      write(*,*) '----------------------------'

c      write(*,*) 'rin=',rin
c                     write(*,*) 'icbfirst=',icb(1)
c                     write(*,*) 'r(icbfirst)',rbb(icb(1))
c$$$      write(*,*) 'er=',vecin(1)
c$$$      write(*,*) 'ksir=',vecin(2)
c      icb=minloc(abs(rbb-rin))
      if(rbb(icb(1)).lt.rin .or. rin==1.)then 
         icb(1)=icb(1)-1
         end if
c               write(*,*) 'icb=',icb(1)
c .or. icb(1)==10000-2 .or. icb(1)==10000-3  .or. icb(1)==3 .or. icb(1)==4
      if (icb(1)==10000-1  .or. icb(1)==2 )  then
         last=1
         call E01AAF(rbb(icb(1):icb(1)+1), frb(icb(1):icb(1)+1), Cfedge, 2, 1, 1, rin)
         vecout(1)=vecin(2)/Cfedge(last)
c$$$         write(*,*) 'f=',Cfedge(last)
         call E01AAF(rbb2(icb(1):icb(1)+1), grb(icb(1):icb(1)+1), Cgedge, 2,1, 1, rin)
         vecout(2)=vecin(1)*Cgedge(last)
c       write(*,*) 'g=',Cgedge(last)

c$$$      write(*,*) 'der=',vecout(1)
c$$$      write(*,*) 'dksir=',vecout(2)
c$$$      write(*,*) 'this was the edge'
      else

      
c      write(*,*) rb(101)
c nint must be even
      
      last=n*(n+1)/2
      np=(n+1)/2
c      last=(10000-1)*10000/2
      
c      write(*,*) icb(1)
      
      call E01AAF(rbb((icb(1)-(np-1)):(icb(1)+np)), frb((icb(1)-(np-1)):(icb(1)+np)), Cf,n+1, last, n, rin)
c      write(*,*) 'fr=',C(last)
      vecout(1)=vecin(2)/Cf(last)
c$$$      write(*,*) 'f=',Cf(last)
c$$$      
      call E01AAF(rbb2((icb(1)-(np-1)):(icb(1)+np)), grb((icb(1)-(np-1)):(icb(1)+np)), Cg, n+1, last, n, rin)
c      Cg(last)=gr(icb(1))
      vecout(2)=vecin(1)*Cg(last)
c$$$       write(*,*) 'g=',Cg(last)
c$$$            write(*,*) 'der=',vecout(1)
c$$$      write(*,*) 'dksir=',vecout(2)
      end if
      end subroutine

c$$$      subroutine outputfunc(xsol,y)
c$$$      real xsol
c$$$      real, dimension(2)::y
c$$$      xsol=0.
c$$$c$$$      write(*,*) 'r=',xsol
c$$$c$$$      write(*,*) 'er=',y(1)
c$$$      end subroutine
c$$$
c$$$      function gfunc(x,y)
c$$$      real x
c$$$      real, dimension(2)::y
c$$$c      real rc
c$$$c      common /fg_block/ rc
c$$$c      gfunc=x-rc
c$$$      gfunc=y(1)-1.186127220166149E-04
c$$$      end function




cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccc   subroutine outer cccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine outer(m,n,size,r,p,dp,q,dq,Bz,Bp,k,order,
     1 d_right,d_left,DI,ifail,icb,iendleft,iendright,istartleft,alpha,prec)
      implicit none
      real D02EJW,D02EJX,D02EJY
      external D02EJW,D02EJX,D02EJY
      integer size,m,n,order,ic,ifail,fact,icb
      real, dimension(10000):: r,p,q,dp,dq,Bz,Bp,fr,gr,rb,grbis
      real, dimension(ic-1):: er_l_left_array, er_s_left_array
      real, dimension(10000-ic)::er_l_right_array,er_s_right_array
      real, dimension(10000)::ksir,er,ksil,el
      real k,qc,rc,step,sigma_l,sigma_s,d_left,d_right
      real d_prime
      real er_l_left,er_s_left,er_l_right,er_s_right
      real DI
      real mr,nr
c      integer, dimension(1)::ic_array
      real, dimension(14)::dfrc,erest,dgrc,errors_dfrc,errors_dgrc
      real,dimension(order)::er_l,er_s
      real, dimension(order+1)::ksi_l,ksi_s
      real eps1_l,eps2_l,eps1_s,eps2_s
      real eps,er0,el0,der0,del0,ksir0,ksil0
      real z
      real ksi_l_right_shoot, ksi_s_right_shoot, ksi_l_left_shoot,ksi_s_left_shoot
      real er_l_right_shoot, er_s_right_shoot, er_l_left_shoot,er_s_left_shoot
      real der_l_right_shoot, der_s_right_shoot, der_l_left_shoot,der_s_left_shoot
      real ksi1_l, ksi2_l, ksi3_l, ksi1_s, ksi2_s, ksi3_s
      real alpha,dqc,ta,d_right_shoot,d_left_shoot
      integer istartleft,iendright,iendleft,i,j
      integer iw
      real zright,zleft,d_prime_shoot
      real, dimension(2)::vecin,vecout,vec_r0,vec_l0
      real, dimension(40)::W
      real, dimension(14*2+50)::W2
      real tol,rightstart,rightend,leftstart,leftend
      real der_s_left_bis,der_l_left_bis,der_s_right_bis,der_l_right_bis
      real er_s_left_bis,er_l_left_bis,er_s_right_bis,er_l_right_bis
      real Cr
      integer ninterp
      real d_left_shoot_test
      real matchright, matchleft
      real prec
c      real c1,c2,c3,c4
      common /fg_block/ fr,gr,rc,ic,rb
c      ,rc,ic,rb
      common /size_interp/ ninterp
      ninterp=3
c      common /plot_block/ c1,c2,c3,c4
      step=1./real(size-1)
c supposes a=1
      qc=real(m)/real(n)
      ic=icb
      rb=r
      mr=real(m)
      nr=real(n)
c      ic_array=minloc(abs(q-qc))
c      write(*,*) 'minloc =',ic_array
c      write(*,*) 'minloc =',ic_array(1)
c      ic=ic_array(1)
      rc=r(icb)
      dqc=dq(icb)

c$$$c.....initialize plotting
c$$$      call ncarcgm(1,'dq.cgm')
c$$$c....make plot
c$$$      call mapg(0.,1.,minval(dq),maxval(dq),.1,1.,.1,1.)
c$$$c      call agsetf('X/LOGARITHMIC.',1.) 
c$$$      ta=2.
c$$$      call trace(r,dq,size
c$$$     1 ,-1,-1,0.,0.)
c$$$      call frame(0)
c$$$c.....finalize plot file
c$$$      call plote
c$$$      write(*,*) 'dqc',dqc
c$$$      write(*,*) 'mr',mr
c$$$      write(*,*) 'nr',nr
c$$$      write(*,*) 'k',k
c$$$      write(*,*) 'alpha',alpha
c$$$      do 10 i=1,size
c$$$         write(*,*) '------------'
c$$$         write(*,*) 'r=',r(i)
c$$$         write(*,*) 'q=',q(i)
c$$$         write(*,*) 'dq=',dq(i)
c$$$c         write(*,*) 'fr',fr(i) 
c$$$c         write(*,*) 'gr=',gr(i) 
c$$$c         if (i==50) stop
c$$$c         if(abs(x-(rc-(2*i-1)*step)).le.step*0.01)funfr=fr(ic-(2*i-1))
c$$$ 10   continue
c      write(*,*) q
      
      fr = r * Bp**2 * (mr-nr*q)**2 / (mr**2+nr**2*k**2*r**2)
c$$$      gr = 2.*nr**2*k**2*r**2*dp/(mr**2+nr**2*k**2*r**2) + 1./r * Bp**2 * 
c$$$     1(mr-nr*q)**2*(nr**2*k**2*r**2+mr**2-1)/(mr**2+nr**2*k**2*r**2)+
c$$$     22.* nr**2*k**2*r * Bp**2 * (nr**2*q**2-mr**2)/ 
c$$$     3 (mr**2+nr**2*k**2*r**2)**2

c$$$      gr=nr**2*k**2*r*Bp**2/(mr**2+nr**2*k**2*r**2)*(-alpha/4.*(dq/k**2)**2+(1.+(mr**2-1.)/(nr**2*k**2*r**2))*(mr-nr*q)**2+2.*(nr**2*q**2-mr**2)/(mr**2+nr**2*k**2*r**2) )
c$$$      gr(1)=gr(2)
c$$$      gr(size)=gr(size-1)
      gr=2.*nr**2*k**2*r**2/(k**2*nr**2*r**2+mr**2)*dp + Bp**2 / r *(mr-nr*q)**2*(k**2*nr**2*r**2+m**2-1.)/(k**2*r**2*nr**2+mr**2)
     A + 2.*k**2*nr**2*r*Bp**2*(n**2*q**2-m**2)/(k**2*nr**2*r**2+m**2)**2
c$$$
c$$$      gr=grbis  
      gr(1)=gr(2)
      gr(size)=gr(size-1)
c$$$      write(*,*) 'fr(2000)=',fr(2000)
c$$$      call outerfunc((r(2000)+r(4000))/2.,vecin,vecout)
c$$$      write(*,*) 'fr(4000)=',fr(4000)

c.....initialize plotting
      call ncarcgm(1,'fr.cgm')
c....make plot
      call mapg(0.,1.,minval(fr),maxval(fr),.1,1.,.1,1.)
c      call agsetf('X/LOGARITHMIC.',1.) 
      ta=2.
      call trace(r,fr,size
     1 ,-1,-1,0.,0.)
      call frame(0)
c.....finalize plot file
      call plote



c.....initialize plotting
      call ncarcgm(1,'gr.cgm')
c....make plot
      call mapg(0.,1.,minval(gr),maxval(gr),.1,1.,.1,1.)
c      call agsetf('X/LOGARITHMIC.',1.)
      ta=2.
      call trace(r,gr,size,-1,-1,0.,0.)
      call colora('green')
      call trace(r,grbis,size,-1,-1,0.,0.)
      call frame(0)
c.....finalize plot file
      call plote


c$$$      write(*,*) 'test'
c$$$      do 10 i=1,size
c$$$         write(*,*) '------------'
c$$$         write(*,*) 'r=',r(i)
c$$$         write(*,*) 'fr',fr(i) 
c$$$         write(*,*) 'gr=',gr(i) 
c$$$c         if (i==50) stop
c$$$c         if(abs(x-(rc-(2*i-1)*step)).le.step*0.01)funfr=fr(ic-(2*i-1))
c$$$ 10   continue

c      write(*,*) 'frc =', funfr(rc)

c      stop
            write(*,*) 'g0=',gr(icb)
      call D04AAF(rc,14,step,dfrc,erest,funfr,ifail)
      errors_dfrc=erest/abs(dfrc)
      call D04AAF(rc,14,step,dgrc,erest,fungr,ifail)
      errors_dgrc=erest/abs(dgrc)
      
      do 12 i=1,14
      write(*,*) '------- i=',i   
      write(*,*) 'error_dfrc =', errors_dfrc(i)
      write(*,*) 'error_dgrc =', errors_dgrc(i)
      write(*,*) 'dfrc =', dfrc(i)
      write(*,*) 'dgrc =', dgrc(i)

 12   continue
      fact=1
      do 33 i=2,order+2
         fact=fact*i
         dfrc(i)=dfrc(i)/real(fact)
         dgrc(i)=dgrc(i)/real(fact)
 33   continue
      write(*,*) 'g0=',gr(icb)
      write(*,*) 'f2=',dfrc(2)
c      DI=-0.25-gr(ic)/(dfrc(2))
c      DI=E+F+H-0.25
      write(*,*) 'DI=alpha-1/4',DI
      write(*,*) 'DI=-0.25-g0/f2',-0.25-gr(icb)/(dfrc(2))
      if (-DI.ge.0.) then
         sigma_l=-0.5-sqrt(-DI)
         sigma_s=-0.5+sqrt(-DI)

      else
         write(*,*) 'complex sigma; Suydam criterion not satisfied'
         ifail=1
         return
      end if
      

      write(*,*) 'sigma_l=',sigma_l
      write(*,*) 'sigma_s=',sigma_s
      write(*,*) 'sigma_s-sigma_l=',sigma_s-sigma_l
      er_l(1)=-(dfrc(3)*sigma_l*(sigma_l+2)-dgrc(1))/(dfrc(2)
     1*(sigma_l+1)*(sigma_l+2)-gr(ic))
      er_s(1)=-(dfrc(3)*sigma_s*(sigma_s+2)-dgrc(1))/(dfrc(2)
     1*(sigma_s+1)*(sigma_s+2)-gr(ic))
      do 13 i=2,order
         er_l(i)=(dfrc(i+2)*(sigma_l)
     1 *(sigma_l+i+1)-dgrc(i))
         er_s(i)=(dfrc(i+2)*(sigma_s)
     1 *(sigma_s+i+1)-dgrc(i))
         do 14 j=1,i-1
            er_l(i)=er_l(i)+er_l(j)*(dfrc(i+2-j)*(sigma_l+j)
     1 *(sigma_l+i+1)-dgrc(i-j))
            er_s(i)=er_s(i)+er_s(j)*(dfrc(i+2-j)*(sigma_s+j)
     1 *(sigma_s+i+1)-dgrc(i-j))
 14   continue
      er_l(i)=-er_l(i)/(dfrc(2)*(sigma_l+i)*(sigma_l+i+1)-gr(ic))
      er_s(i)=-er_s(i)/(dfrc(2)*(sigma_s+i)*(sigma_s+i+1)-gr(ic))
 13   continue

c better precision for order=2
      eps1_l=-sigma_l/2.*(sigma_l*(sigma_l+2)*dfrc(3)-dgrc(1))/gr(ic)
      eps2_l=-(sigma_l*(sigma_l+3)*dfrc(4)-dgrc(2)+eps1_l*((sigma_l+1)*(sigma_l+3)*dfrc(3)-dgrc(1)))/((sigma_l+2)*(sigma_l+3)*dfrc(2)-gr(ic))

      eps1_s=-sigma_s/2.*(sigma_s*(sigma_s+2)*dfrc(3)-dgrc(1))/gr(ic)
      eps2_s=-(sigma_s*(sigma_s+3)*dfrc(4)-dgrc(2)+eps1_s*((sigma_s+1)*(sigma_s+3)*dfrc(3)-dgrc(1)))/((sigma_s+2)*(sigma_s+3)*dfrc(2)-gr(ic))

      ksi1_l=gr(ic)/(sigma_l+1.)
      ksi2_l=0.5*(dgrc(1)-sigma_l**2*dfrc(3))
      ksi3_l=sigma_l*dfrc(4)+(sigma_l+1)*dfrc(3)*eps1_l+(sigma_l+2)*dfrc(2)*eps2_l

      ksi1_s=gr(ic)/(sigma_s+1)
      ksi2_s=0.5*(dgrc(1)-sigma_s**2*dfrc(3))
      ksi3_s=sigma_s*dfrc(4)+(sigma_s+1)*dfrc(3)*eps1_s+(sigma_s+2)*dfrc(2)*eps2_s

      write(*,*) 'ksi1_l=',ksi1_l
      write(*,*) 'bis', sigma_l*dfrc(2)
      write(*,*) 'ksi2_l=',ksi2_l
      write(*,*) 'ksi3_l=',ksi3_l
      write(*,*) 'ksi1_s=',ksi1_s
      write(*,*) 'ksi2_s=',ksi2_s
      write(*,*) 'ksi3_s=',ksi3_s


      er_l(1)=eps1_l
      er_s(1)=eps1_s
      er_l(2)=eps2_l
      er_s(2)=eps2_s

      matchleft=real(step*iendleft)
      matchright=real(step*iendright)
      

      ksi_l(1)=ksi1_l
      ksi_l(2)=ksi2_l
      ksi_l(3)=ksi3_l

      ksi_s(1)=ksi1_s
      ksi_s(2)=ksi2_s
      ksi_s(3)=ksi3_s


c$$$
c$$$      er_l_right_shoot=1.
c$$$      er_s_right_shoot=1.
c$$$      er_l_left_shoot=1.
c$$$      er_s_left_shoot=1.
      

      er_s_right_bis=(matchright)**(sigma_s)
      er_l_right_bis=(matchright)**(sigma_l)      

      er_s_left_bis=(matchleft)**(sigma_s)
      er_l_left_bis=(matchleft)**(sigma_l)

         write(*,*) 'er_l_right_zero_order =',(matchright)**(sigma_s)
         write(*,*) 'er_s_right_zero_order =',(matchright)**(sigma_l)     
         write(*,*) 'er_l_left_zero_order =',(matchleft)**(sigma_s)
         write(*,*) 'er_s_left_zero_order =',(matchleft)**(sigma_l)
c$$$      der_l_right_shoot=0.
c$$$c      der_l_right_bis=sigma_l*(step
c$$$      der_s_right_shoot=0.
c$$$      der_l_left_shoot=0.
c$$$      der_s_left_shoot=0.

      der_s_right_bis=sigma_s*(matchright)**(sigma_s-1.)
      der_l_right_bis=sigma_l*(matchright)**(sigma_l-1.)      

      der_s_left_bis=-sigma_s*(matchleft)**(sigma_s-1.)
      der_l_left_bis=-sigma_l*(matchleft)**(sigma_l-1.)

               write(*,*) 'der_l_right_zero_order =',sigma_s*(matchright)**(sigma_s-1.)
         write(*,*) 'der_s_right_zero_order =',sigma_l*(matchright)**(sigma_l-1.)        
         write(*,*) 'der_l_left_zero_order =',-sigma_s*(matchleft)**(sigma_s-1.)
         write(*,*) 'der_s_left_zero_order =',-sigma_l*(matchleft)**(sigma_l-1.)


      ksi_l_right_shoot=0.
      ksi_s_right_shoot=0.
      ksi_l_left_shoot=0.
      ksi_s_left_shoot=0.

      do 17 i=1,order
         write(*,*) 'order =',i
c         fact=fact*i

c$$$         er_l_right_shoot=er_l_right_shoot+er_l(i)*(matchright)**i
c$$$         er_s_right_shoot=er_s_right_shoot+er_s(i)*(matchright)**i
c$$$         er_l_left_shoot=er_l_left_shoot+er_l(i)*(-matchleft)**i
c$$$         er_s_left_shoot=er_s_left_shoot+er_s(i)*(-matchleft)**i

c$$$         write(*,*) 'der_l_right_shoot =',er_l(i)*real(i)*(matchright)**(i-1)
c$$$         write(*,*) 'der_s_right_shoot =',er_s(i)*real(i)*(matchright)**(i-1)
c$$$         write(*,*) 'der_l_left_shoot =',er_l(i)*real(i)*(-matchleft)**(i-1)
c$$$         write(*,*) 'der_s_left_shoot =',er_s(i)*real(i)*(-matchleft)**(i-1)

c$$$         der_l_right_shoot=der_l_right_shoot+er_l(i)*real(i)*(matchright)**(i-1)
c$$$         der_s_right_shoot=der_s_right_shoot+er_s(i)*real(i)*(matchright)**(i-1)
c$$$         der_l_left_shoot=der_l_left_shoot+er_l(i)*real(i)*(-matchleft)**(i-1)
c$$$         der_s_left_shoot=der_s_left_shoot+er_s(i)*real(i)*(-matchleft)**(i-1)


         er_l_right_bis=er_l_right_bis+er_l(i)*(matchright)**(sigma_l+i)
         er_s_right_bis=er_s_right_bis+er_s(i)*(matchright)**(sigma_s+i)
         er_l_left_bis=er_l_left_bis+er_l(i)*(-1.)**i*(matchleft)**(sigma_l+i)
         er_s_left_bis=er_s_left_bis+er_s(i)*(-1.)**i*(matchleft)**(sigma_s+i)

         write(*,*) 'er_l_right_shoot_bis =',er_l(i)*(matchright)**(sigma_l+i)
         write(*,*) 'er_s_right_shoot_bis =',er_s(i)*(matchright)**(sigma_s+i)
         write(*,*) 'er_l_left_shoot_bis =',er_l(i)*(-1.)**i*(matchleft)**(sigma_l+i)
         write(*,*) 'er_s_left_shoot_bis =',er_s(i)*(-1.)**i*(matchleft)**(sigma_s+i)


         der_s_right_bis=der_s_right_bis+er_s(i)*(sigma_s+i)*(matchright)**(sigma_s+i-1)
         der_l_right_bis=der_l_right_bis+er_l(i)*(sigma_l+i)*(matchright)**(sigma_l+i-1)
         der_s_left_bis=der_s_left_bis-(-1.)**(i)*er_s(i)*(sigma_s+i)*(matchleft)**(sigma_s+i-1)
         der_l_left_bis=der_l_left_bis-(-1.)**(i)*er_l(i)*(sigma_l+i)*(matchleft)**(sigma_l+i-1)

         write(*,*) 'der_l_right_shoot_bis =',er_s(i)*(sigma_s+i)*(matchright)**(sigma_s+i-1)
         write(*,*) 'der_s_right_shoot_bis =',er_l(i)*(sigma_l+i)*(matchright)**(sigma_l+i-1)
         write(*,*) 'der_l_left_shoot_bis =',(-1.)**(i)*er_s(i)*(sigma_s+i)*(matchleft)**(sigma_s+i-1)
         write(*,*) 'der_s_left_shoot_bis =',(-1.)**(i)*er_l(i)*(sigma_l+i)*(matchleft)**(sigma_l+i-1)
c         er_l_right_shoot_array=er_l_right_shoot+er_l(i)*(r(ic+1:size)-rc)**i
c         er_s_right_shoot_array=er_s_right_shoot+er_s(i)*(r(ic+1:size)-rc)**i
c         er_l_left_shoot_array=er_l_left_shoot+er_l(i)*(r(1:ic-1)-rc)**i
c         er_s_left_shoot_array=er_s_left_shoot+er_s(i)*(r(1:ic-1)-rc)**i
 17   continue

      do 16 i=1,3
         write(*,*) 'i=',i
         write(*,*) 'ksi_l_right_shoot',ksi_l(i)*(matchright)**i
         write(*,*) 'ksi_s_right_shoot',ksi_s(i)*(matchright)**i
         write(*,*) 'ksi_l_left_shoot',ksi_l(i)*(-matchleft)**i
         write(*,*) 'ksi_s_left_shoot',ksi_s(i)*(-matchleft)**i

         ksi_l_right_shoot=ksi_l_right_shoot+ksi_l(i)*(matchright)**i
         ksi_s_right_shoot=ksi_s_right_shoot+ksi_s(i)*(matchright)**i
         ksi_l_left_shoot=ksi_l_left_shoot+ksi_l(i)*(-matchleft)**i
         ksi_s_left_shoot=ksi_s_left_shoot+ksi_s(i)*(-matchleft)**i
 16       continue


c      write(*,*) 'rc =',rc

c$$$      der_l_right_shoot=sigma_l*(step*real(iendright))**(sigma_l-1.)*er_l_right_shoot+(real(matchright))**sigma_l*der_l_right_shoot
c$$$      der_s_right_shoot=sigma_s*(step*real(iendright))**(sigma_s-1.)*er_s_right_shoot+(real(matchright))**sigma_s*der_s_right_shoot
c$$$      der_l_left_shoot=-sigma_l*(step*real(iendleft))**(sigma_l-1.)*er_l_left_shoot+(real(matchright))**sigma_l*der_l_left_shoot
c$$$      der_s_left_shoot=-sigma_s*(step*real(iendleft))**(sigma_s-1.)*er_s_left_shoot+(real(matchright))**sigma_s*der_s_left_shoot 
c$$$   
c$$$
c$$$
c$$$      er_l_right_shoot=(matchright)**sigma_l*er_l_right_shoot
c$$$      er_s_right_shoot=(matchright)**sigma_s*er_s_right_shoot
c$$$      er_l_left_shoot=(matchleft)**sigma_l*er_l_left_shoot
c$$$      er_s_left_shoot=(matchleft)**sigma_s*er_s_left_shoot

      ksi_l_right_shoot=(matchright)**sigma_l*ksi_l_right_shoot
      ksi_s_right_shoot=(matchright)**sigma_s*ksi_s_right_shoot
      ksi_l_left_shoot=(matchleft)**sigma_l*ksi_l_left_shoot
      ksi_s_left_shoot=(matchleft)**sigma_s*ksi_s_left_shoot
 
c      d_right_shoot=1.
c      d_left_shoot=1.
c      do 99 iendleft=1,20
c         do 98 iendright=1,5
c         iendright=iendleft
c      istartleft=-1
c min is 1
c      iendleft=30
c default is 1
c      iendright=30
c default is 1

      er_s_right_shoot=er_s_right_bis
      er_l_right_shoot=er_l_right_bis
      er_s_left_shoot=er_s_left_bis
      er_l_left_shoot=er_l_left_bis
      
c$$$      der_s_right_shoot=der_s_right_bis
c$$$      der_l_right_shoot=der_l_right_bis
c$$$      der_s_left_shoot=der_s_left_bis
c$$$      der_l_left_shoot=der_l_left_bis

      der_s_right_shoot=ksi_s_right_shoot/fr(ic+iendright)
      der_l_right_shoot=ksi_l_right_shoot/fr(ic+iendright)
      der_s_left_shoot=ksi_s_left_shoot/fr(ic-iendleft)
      der_l_left_shoot=ksi_l_left_shoot/fr(ic-iendleft)

      write(*,*) '-------------'
      write(*,*) 'iendright=',iendright
      write(*,*) 'iendleft=',iendleft


c$$$      write(*,*) 'er_l_right_shoot=',er_l_right_shoot
c$$$      write(*,*) 'er_s_right_shoot=',er_s_right_shoot
c$$$      write(*,*) 'er_l_left_shoot=',er_l_left_shoot
c$$$      write(*,*) 'er_s_left_shoot=',er_s_left_shoot

      write(*,*) 'er_l_right_bis=',er_l_right_bis
      write(*,*) 'er_s_right_bis=',er_s_right_bis

      write(*,*) 'er_l_left_bis=',er_l_left_bis
      write(*,*) 'er_s_left_bis=',er_s_left_bis

c$$$      write(*,*) 'der_l_right_shoot=',der_l_right_shoot
c$$$      write(*,*) 'der_s_right_shoot=',der_s_right_shoot
c$$$      write(*,*) 'der_l_left_shoot=',der_l_left_shoot
c$$$      write(*,*) 'der_s_left_shoot=',der_s_left_shoot

      write(*,*) 'der_l_right_bis=',der_l_right_bis
      write(*,*) 'der_s_right_bis=',der_s_right_bis

      write(*,*) 'der_l_left_bis=',der_l_left_bis
      write(*,*) 'der_s_left_bis=',der_s_left_bis
      
      write(*,*) 'der_l_right_bisbis=',ksi_l_right_shoot/fr(ic+iendright)
      write(*,*) 'der_s_right_bisbis=',ksi_s_right_shoot/fr(ic+iendright)

      write(*,*) 'der_l_left_bisbis=',ksi_l_left_shoot/fr(ic-iendleft)
      write(*,*) 'der_s_left_bisbis=',ksi_s_left_shoot/fr(ic-iendleft)
      write(*,*) 'fr_right',fr(ic+iendright)
      write(*,*) 'fr_left',fr(ic-iendleft)
c$$$      write(*,*) 'der_l_right_shootbis=',ksi_l_right_shoot/fr(ic+iendright)
c$$$      write(*,*) 'der_s_right_shootbis=',ksi_s_right_shoot/fr(ic+iendright)
c$$$      write(*,*) 'der_l_left_shootbis=',ksi_l_left_shoot/fr(ic-iendleft)
c$$$      write(*,*) 'der_s_left_shootbis=',ksi_s_left_shoot/fr(ic-iendleft)


ccccccccccccccccccccccccccccccccccccccccc
cccccccccccccc RUNGE KUTTA cccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccc      

ccccc RIGHT REGION cccccc 
      ifail=0
      vec_r0(1)=0.
      vec_r0(2)=fr(size)*(-1.e-1)
      iw=14*2+50
      rightstart=1.
      rightend=rc+iendright*step
      tol=1.e-8
      call D02EJF(rightstart,rightend,2, vec_r0,outerfunc,D02EJY,tol,    'D', D02EJX,D02EJW,W2,iw,ifail)
c      call D02BJF(1.,rc+iendright*step,2, vec_r0,outerfunc,0.01,'D', D02BJX,D02BJW,W,ifail)
c      write(*,*) 'ouf'
      write(*,*) 'zright_RK=',vec_r0(1)/(vec_r0(2)/fr(ic+iendright))
      write(*,*) 'er=',vec_r0(1)
      write(*,*) 'der=',vec_r0(2)/fr(ic+iendright)
      write(*,*) 'fr=',fr(ic+iendright)
      zright=vec_r0(1)/(vec_r0(2))*fr(ic+iendright)
c      zrightinv=
ccccc LEFT REGION cccccc
      ifail=0
      tol=1.e-8
      iw=14*2+50
      if(m==1)then
         vec_l0(1)=-1.
         vec_l0(2)=0.
      else
         vec_l0(1)=(istartleft*step)**(m-1)
         vec_l0(2)=real(m-1)*(istartleft*step)**(m-2)*fr(istartleft+1)
      end if
      leftstart=istartleft*step
      leftend=rc-iendleft*step
      call D02EJF(leftstart,leftend,2, vec_l0,outerfunc,D02EJY,tol,        'D', D02EJX,D02EJW,W2,iw,ifail)


      write(*,*) 'zleft_RK=',vec_l0(1)/(vec_l0(2)/fr(ic-iendleft))
      write(*,*) 'er=',vec_l0(1)
      write(*,*) 'der=',vec_l0(2)/fr(ic-iendleft)
      write(*,*) 'fr=',fr(ic-iendleft)
      zleft=vec_l0(1)/(vec_l0(2)/fr(ic-iendleft))
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccc






      d_right_shoot=(zright*der_l_right_shoot-er_l_right_shoot)/(er_s_right_shoot-zright*der_s_right_shoot)
      write(*,*) 'zright*der_l_right=',zright*der_l_right_shoot
      write(*,*) 'zright*der_s_right=',zright*der_s_right_shoot
      write(*,*) 'numerator=',(zright*der_l_right_shoot-er_l_right_shoot)
      write(*,*) 'denominator=',(er_s_right_shoot-zright*der_s_right_shoot)
      d_left_shoot=(zleft*der_l_left_shoot-er_l_left_shoot)/(er_s_left_shoot-zleft*der_s_left_shoot)
      d_prime_shoot=d_right_shoot+d_left_shoot
      write(*,*) 'd_prime_shoot =',d_prime_shoot
      d_right=d_right_shoot
      d_left=d_left_shoot
      write(*,*) 'd_right =',d_right
  
      write(*,*) 'd_left =',d_left
      write(*,*) 'W0right',er_l_right_shoot*fr(ic+iendright)*der_s_right_shoot-fr(ic+iendright)*der_l_right_shoot*er_s_right_shoot
      write(*,*) 'W0bisright',er_l_right_shoot*ksi_s_right_shoot-ksi_l_right_shoot* er_s_right_shoot

            write(*,*) 'W0left',er_l_left_shoot*fr(ic-iendleft)*der_s_left_shoot-fr(ic-iendleft)*der_l_left_shoot*er_s_left_shoot
      write(*,*) 'W0bisleft',er_l_left_shoot*ksi_s_left_shoot-ksi_l_left_shoot* er_s_left_shoot

      write(*,*) 'Wexact',(2.*sigma_s+1.)*dfrc(2)
      prec=abs(er_l_right_shoot*fr(ic+iendright)*der_s_right_shoot-fr(ic+iendright)*der_l_right_shoot*er_s_right_shoot-(2.*sigma_s+1.)*dfrc(2))
c      write(*,*) 'd_right_shoot2',
c$$$      write(*,*) 'der/er_right=',1/zright
c$$$      write(*,*) 'der/er_left=',1/zleft
c$$$      write(*,*) 'dprime=',1/zright-1/zleft
      write(*,*) '***********************************************'
c$$$      der_l_right_shoot=ksi_l_right_shoot/fr(ic+iendright)
c$$$      der_s_right_shoot=ksi_s_right_shoot/fr(ic+iendright)
c$$$      der_l_left_shoot=ksi_l_left_shoot/fr(ic-iendleft)
c$$$      der_s_left_shoot=ksi_s_left_shoot/fr(ic-iendleft)


      




      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccc   function funfr   cccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function funfr(x)
      implicit none
      real x,funfr
      real step
      real rc
      integer size,ic,i
      common /size_block/ size
      real,dimension(10000) ::r,fr,gr,rb
      common /fg_block/ fr,gr,rc,ic,rb
      step=1./real(size-1)
      r=rb
c supposes a=1
      if (abs(x-rc).le.step*0.01) then
         funfr=fr(ic) 
c         write(*,*) 'i=0'
         return
      end if

      do 10 i=1,10
         if(abs(x-(rc+(2*i-1)*step)).le.step*0.01) then 
            funfr=fr(ic+(2*i-1))
c            write(*,*) 'i_plus =',i
            return
          end if

          if(abs(x-(rc-(2*i-1)*step)).le.step*0.01) then 
            funfr=fr(ic-(2*i-1))
c            write(*,*) 'funfr =',funfr
c            write(*,*) 'i_minus =',i
            return
          end if

c         if(abs(x-(rc-(2*i-1)*step)).le.step*0.01)funfr=fr(ic-(2*i-1))
 10   continue
      write(*,*) 'oups:',x
      
      end function
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccc   function fungr   cccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function fungr(x)
      implicit none
      real x,fungr
      real step
      real rc
      integer size,ic,i
      common /size_block/ size
      real,dimension(10000) ::r,fr,gr,rb
      common /fg_block/ fr,gr,rc,ic,rb
      step=1./real(size-1)
      r=rb
c supposes a=1
      if (abs(x-rc).le.step*0.01) then
         fungr=gr(ic) 
c         write(*,*) 'i_g=0'
         return
      end if

      do 10 i=1,10
         if(abs(x-(rc+(2*i-1)*step)).le.step*0.01) then 
            fungr=gr(ic+(2*i-1))
c            write(*,*) 'i_plus_g =',i
            return
          end if

          if(abs(x-(rc-(2*i-1)*step)).le.step*0.01) then 
            fungr=gr(ic-(2*i-1))
c            write(*,*) 'fungr =',fungr
c            write(*,*) 'i_minus_g =',i
            return
          end if

c         if(abs(x-(rc-(2*i-1)*step)).le.step*0.01)fungr=gr(ic-(2*i-1))
 10   continue
      write(*,*) 'oups:',x
       
      end function




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccc dpfunc ccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      function dpfunc(nq0)
       implicit none
       real dpfunc
      integer size,m,n,order,size_plot,istartleft,iendright,iendleft
      real step
      real d_right,d_left
      integer ifail,ic
   
      real E,F,G,H,K0
      real rmatch
      real, dimension(10000):: p,q,dp,dq,Bz,Bp
      real, dimension(10000):: r
      real, dimension(100,3)::gamma_dless_array,
     1 gamma_array
c      real, dimension(100)::eta_array,det_array,nq0_array,d_prime_shoot_array
c      real, dimension(400)::d_right_array,d_left_array,i_array,gamma_array2
      real, dimension(100)::dp_array,prec_array
      real, dimension(100)::Cr_array
      real, dimension(100)::iend_array
      real k
      real gamma_dless
      real Bpc,dqc,rc
      real eta
      real rho
      real Lr,gammar,DI,Lr_DI,DIbis
      real ta
      complex de,do
      real q0,nq0,alpha
      integer i
      real prec
      real Bp0
c      common /delta_block/ d_right,d_left,E,F,G,H,K0,Lr_DI,Lr
      
c      gamma=0.01
      
c      epsd = 0.0
ccccccccccccccccccc
c Fig 6a) na/R=3
ccccccccccccccccccc
      k=3./4.
      size=10000
      step=1./(size-1.)
      alpha=0.5
      m=1
      n=4
      rho=1
      order=3
      Bp0=2.
c      nq0=2.8
            write(*,*) '----------------------'
         write(*,*) 'm= ',m,'; n= ',n
               write(*,*) 'nq0 =',nq0
c      call linspace(1.3,1.6,nq0_array,size_plot)
c      do 67 i=1,size_plot
c         write(*,*) 'nqo =',nq0_array(i)
         q0=nq0/real(n)  
          write(*,*) 'qo =',q0
      call equi(size,r,p,dp,q,dq,Bz,Bp,k,1,alpha,q0,Bp0)
      
c      do 10  n=1,10
c         do 20 m=n+1,11
      
c      eta=1.0e-5
c      size_plot=100
c      step=1.0e-4/real(size_plot)

c      eta_array = (/((i*step),i=1,size_plot)/)
  
      ifail=0

      
     

      call rational_surf(m,n,size,r,p,dp,q,dq,Bz,Bp,k,E,F,G,H,K0,ifail,
     1 dqc,Bpc,rc,ta,ic)
            write(*,*) 'DI=E+F+H-0.25',E+F+H-0.25
c       DI=E+F+H-0.25
      DI=(alpha-1.)/4.
      if (ifail==1) then
         
         write(*,*) 'Rational surface is outside'
c         cycle
      end if
c      call linspace(33.,50.,iend_array,10)
      ifail=0
      istartleft=100
c      do 10 i=1,10
c      iendright=int(iend_array(i))
c      iendleft=int(iend_array(i))
      iendright=140
      iendleft=140
c      call logspace(1.e-8,1.e-1,Cr_array,100)
c      do 33 i=1,100
c         iendleft=(i-1)*2+2
c         iendright=(i-1)*2+2
c         iend_array(i)=iendleft
c       write(*,*) i
c      iendleft=20
c      iendright=20
     
c      DI=(alpha-1.)/4.
      call outer(m,n,size,r,p,dp,q,dq,Bz,Bp,k,order,d_right,d_left,DI
     1 ,ifail,ic,iendleft,iendright,istartleft,alpha,prec)
      dpfunc=d_right+d_left
c      dp_array(i)=dpfunc
      
c 33   continue
      
c$$$
c$$$      call ncarcgm(1,'prec.cgm')
c$$$c....make plot
c$$$      call mapg(minval(iend_array*step),maxval(iend_array*step),minval(prec_array),maxval(prec_array),.1,1.,.1,1.)
c$$$c      call agsetf('X/LOGARITHMIC.',1.)
c$$$c      ta=2.
c$$$      call colora('green')
c$$$      call trace(iend_array*step,prec_array,100,-1,-1,0.,0.)
c$$$      call frame(0)
c$$$c.....finalize plot file
c$$$      call plote
c$$$c$$$
c$$$c$$$
c$$$      call ncarcgm(1,'iend_d_prime.cgm')
c$$$c....make plot
c$$$      call mapg(minval(iend_array*step),maxval(iend_array*step),minval(dp_array),maxval(dp_array),.1,1.,.1,1.)
c$$$c      call agsetf('X/LOGARITHMIC.',1.)
c$$$c      ta=2.
c$$$      call colora('green')
c$$$      call trace(iend_array*step,dp_array,100,-1,-1,0.,0.)
c$$$      call frame(0)
c$$$c.....finalize plot file
c$$$      call plote

      

c$$$      call ncarcgm(1,'Cr.cgm')
c$$$c....make plot
c$$$      call mapg(minval(log10(Cr_array)),maxval(log10(Cr_array)),minval(dp_array),maxval(dp_array),.1,1.,.1,1.)
c$$$c      call agsetf('X/LOGARITHMIC.',1.)
c$$$c      ta=2.
c$$$      call colora('green')
c$$$      call trace(log10(Cr_array),dp_array,100,-1,-1,0.,0.)
c$$$      call frame(0)
c$$$c.....finalize plot file
c$$$      call plote


      end function
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccc  subroutine main_fig6 ccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine main_fig6(size,size_plot) 
      implicit none
      integer size,m,n,order,size_plot,i
      real step
      real d_right,d_left
      integer ifail,ic,j
      real E,F,G,H,K0
      real rmatch
      real, dimension(size):: p,q,dp,dq,Bz,Bp
      real, dimension(size):: r
      real, dimension(size_plot,3)::gamma_dless_array,
     1 gamma_array
      real, dimension(size_plot)::eta_array,det_array,nq0_array,d_prime_shoot_array
      real, dimension(400)::d_right_array,d_left_array,i_array,gamma_array2
      real k
      real gamma_dless
      real Bpc,dqc,rc
      real eta
      real rho
      real Lr,gammar,DI,Lr_DI,DIbis
      real ta
      complex de,do
      real nq0s
      real alpha,gamlim,gamma
c      common /dp_block/ m,n,size,r,p,dp,q,dq,Bz,Bp,k,order,iendleft,iendright,istartleft,alpha
c      common /delta_block/ d_right,d_left,E,F,G,H,K0,Lr_DI,Lr
c      write(*,*) dpfunc(1.3)
c      call C05ADF(1.1,2.,1.e-6,0.,dpfunc,nq0s,ifail)
c      write(*,*) 'nq0s =',nq0s
      

c      gamma=0.01
      
c      epsd = 0.0
c$$$      alpha=0.5
c$$$      m=1
c$$$      n=4
c$$$      rho=1
c$$$      order=2
c$$$            write(*,*) '----------------------'
c$$$         write(*,*) 'm= ',m,'; n= ',n
c$$$      
ccccccccccccccccccc
c Fig 6a) na/R=3
ccccccccccccccccccc
      call linspace(1.4,2.2,nq0_array,size_plot)
      do 67 i=1,size_plot
         d_prime_shoot_array(i)=dpfunc(nq0_array(i))
c         stop
 67      continue
         
c      d_prime_shoot_array(1)=dpfunc(7.)
c      stop
      call ncarcgm(1,'d_prime.cgm')
c....make plot
      call mapg(minval(nq0_array),maxval(nq0_array),minval(d_prime_shoot_array),maxval(d_prime_shoot_array),.1,1.,.1,1.)
c      call agsetf('X/LOGARITHMIC.',1.)
c      ta=2.
      call colora('green')
      call trace(nq0_array,d_prime_shoot_array,size_plot,-1,-1,0.,0.)
      call frame(0)
c.....finalize plot file
      call plote

      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccc  subroutine main_fig7 ccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine main_fig7(size,size_plot) 
      integer size,m,n,order,size_plot
      real step
      real d_right,d_left
      integer ifail,ic
   
      real E,F,G,H,K0
      real rmatch
      real, dimension(size):: p,q,dp,dq,Bz,Bp
      real, dimension(size):: r
      real, dimension(size_plot,1)::gamma_dless_array,
     1 gamma_array
      real, dimension(1)::ta_array
      real, dimension(size_plot)::eta_array,det_array
      real, dimension(400)::d_right_array,d_left_array,i_array,gamma_array2
      real k
      real gamma_dless
      real Bpc,dqc,rc
      real eta
      real rho
      real Lr,gammar,DI,Lr_DI,DIbis
      real ta
      real prec
      real Bp0
      complex de,do
      
      common /delta_block/ d_right,d_left,E,F,G,H,K0,Lr_DI,Lr

cccccccccccccccccccccccccccccc
cccccc fig7, solid, alpha=0.5
cccccccccccccccccccccccccccccc
      
c      gamma=0.01
      
c      epsd = 0.0
      alpha=0.5
      Bp0=2.
c      Bp0=1.
      do 77 j=1,1
c         if (j==1) Bp0=1.
c Bp0=1.
c alpha=0.5
c         if (j==2) Bp0=1.75
c Bp0=0.85
c alpha=0.1
c         if (j==3) Bp0=2.
c Bp0=0.65
c alpha=0.04
c         if (j==4) alpha=0.02
c         if (j==5) alpha=0.01
         k=1./2.
         
      call equi(size,r,p,dp,q,dq,Bz,Bp,k,1,alpha,0.4,Bp0)
c      do 10  n=1,10
c         do 20 m=n+1,11
      m=1
      n=4
      rho=1
       order=2
c      eta=1.0e-5
c      size_plot=100
c      step=1.0e-4/real(size_plot)

c      eta_array = (/((i*step),i=1,size_plot)/)
  
      ifail=0

      
     
      write(*,*) '----------------------'
         write(*,*) 'm= ',m,'; n= ',n
      
      call rational_surf(m,n,size,r,p,dp,q,dq,Bz,Bp,k,E,F,G,H,K0,ifail,
     1 dqc,Bpc,rc,ta_array(j),ic)
      
      if (ifail==1) then
         
         write(*,*) 'Rational surface is outside'
c         cycle
      end if
      ifail=0
      istartleft=100
      iendright=140
      iendleft=140
c      do 33 i=1,400
c         iendleft=i+70
c         iendright=i+70
c       write(*,*) i
c      iendleft=20
c      iendright=20
      DI=E+F+H-0.25
      call outer(m,n,size,r,p,dp,q,dq,Bz,Bp,k,order,d_right,d_left,DI
     1 ,ifail,ic,iendleft,iendright,istartleft,alpha,prec)

     
      
c      if (ifail==1) cycle
      if (ifail==1) stop
c$$$      eta=1.e-8/ta
c$$$      gammar=(eta*n**2*Bpc**2*dqc**2/(rho*rc**2))**(1./3.)
c$$$      Lr=(rho*eta**2*rc**2/(n**2*Bpc**2*dqc**2))**(1./6.)
c$$$      Lr_DI=Lr**(-2.*sqrt(-DI))
      
      call logspace(1.e-9,1.e-3,eta_array,size_plot)
c      call linspace(-1.e-5,0.25,gamma_dless_array,size_plot)
      do 30 i=1,size_plot
c      gamma_dless=0.001
      ifail=0
c      DIbis=-1./4.+E+F+H
      eta=eta_array(i)
     
       gammar=(eta*real(n)**2*Bpc**2*dqc**2/(rho*rc**2))**(1./3.)

      Lr=(rho*eta**2*rc**2/(real(n)**2*Bpc**2*dqc**2))**(1./6.)
      Lr_DI=Lr**(-2.*sqrt(-DI))
      
c      write(*,*) 'DI =',DI
c       write(*,*) 'DIbis =',DIbis
     
    
c      det_array(i) = det(gamma_dless_array(i,j))
      
   






      write(*,*) '-----------------------------------'
      write(*,*) 'eta =',eta
      write(*,*) 'alpha=',alpha
c      call C05AJF(gamma_dless,0.0001,0.0,det,500,ifail)
      write(*,*) '--- biss ---'
c$$$      if (j==1) gamlim=0.2
c$$$      if (j==2) gamlim=0.07
c$$$      if (j==3) gamlim=0.03
c$$$      if (j==4) gamlim=0.02
c$$$      if (j==5) gamlim=0.015
      gamlim=0.2
      call C05ADF(1.e-9,gamlim,1.e-12,0.,det,gamma_dless,ifail)
      gamma=gammar*gamma_dless
c      gamma_bis=gammar*gamma_dless_bis
      write(*,*) 'gamma dimensionless =',gamma_dless
      write(*,*) 'gamma =',gamma
c      write(*,*) 'gammabis dimensionless =',gamma_dless_bis
c      write(*,*) 'gammabis =',gamma_bis
      gamma_dless_array(i,j)=gamma_dless
      gamma_array(i,j)=gamma


 30   continue


c            write(*,*) '************************************'
c      call C05AJF(gamma_dless,1.e-6,0.0,det,500,ifail)
c      write(*,*) '--- biss ---'
c      call C05ADF(1.e-5,0.06,1.e-6,0.,det,gamma_dless_bis,ifail)
c      gamma=gammar*gamma_dless
c      gamma_bis=gammar*gamma_dless_bis
c      write(*,*) 'gamma dimensionless =',gamma_dless
c      write(*,*) 'gamma =',gamma
c      write(*,*) 'gammabis dimensionless =',gamma_dless_bis
c      write(*,*) 'gammabis =',gamma_bis
 


      write(*,*) 'eta =',eta
      write(*,*) 'log(resistivity) =',log(eta*ta_array(j))
c      write(*,*) 'ln(gamma) =',log(gamma)
      write(*,*) 'd_right =',d_right
      write(*,*) 'd_left =',d_left
            write(*,*) 'DI =',DI
c       write(*,*) 'DIbis =',DIbis
       write(*,*) 'ta =',ta_array(j)
       write(*,*) 'order =',order
       write(*,*) 'lambda =',-2.*sqrt(-DI)+2*order-2
c$$$cccccccc PLOTTING cccccccccccc
c$$$c.....initialize plotting
c$$$      call ncarcgm(1,'det.cgm')
c$$$c....make plot
c$$$      call mapg(minval(gamma_dless_array(:,j)),maxval(gamma_dless_array(:,j)),-0.1,0.1,.1,1.,.1,1.)
c$$$c      call agsetf('X/LOGARITHMIC.',1.)
c$$$c      do 98 i=1,10
c$$$c      ta=1.5+i*0.1
c$$$      call colora('green')
c$$$      call trace(gamma_dless_array(:,j),det_array/maxval(abs(det_array)) ,size_plot,-1,-1,0.,0.)
c$$$
c$$$c 98   continue 
c$$$      call frame(0)
c$$$c.....finalize plot file
c$$$      call plote

       
 77    continue
c 33     continue

c$$$      call ncarcgm(1,'d_right.cgm')
c$$$c....make plot
c$$$      call mapg(minval(i_array),maxval(i_array),-10., 0.,.1,1.,.1,1.)
c$$$c      call agsetf('X/LOGARITHMIC.',1.)
c$$$c      ta=2.
c$$$      call trace(i_array,d_right_array,400,-1,-1,0.,0.)
c$$$      call frame(0)
c$$$c.....finalize plot file
c$$$      call plote
c$$$
c$$$c.....initialize plotting
c$$$      call ncarcgm(1,'d_left.cgm')
c$$$c....make plot
c$$$      call mapg(minval(i_array),maxval(i_array),-20.,0.,.1,1.,.1,1.)
c$$$c      call agsetf('X/LOGARITHMIC.',1.)
c$$$c      ta=2.
c$$$      call trace(i_array,d_left_array,400,-1,-1,0.,0.)
c$$$      call frame(0)
c$$$c.....finalize plot file
c$$$      call plote
c$$$
c$$$c.....initialize plotting
c$$$      call ncarcgm(1,'gamma_startleft.cgm')
c$$$c....make plot
c$$$      call mapg(minval(i_array),maxval(i_array),0.0015,0.0020,.1,1.,.1,1.)
c$$$c      call agsetf('X/LOGARITHMIC.',1.)
c$$$c      ta=2.
c$$$      call trace(i_array,gamma_array2,400,-1,-1,0.,0.)
c$$$      call frame(0)
c$$$c.....finalize plot file
c$$$      call plote


c        stop


c.....initialize plotting
      call ncarcgm(1,'gamma_fig7.cgm')
c....make plot
      call mapg(-9.,-3.,-20.,0.,.1,1.,.1,1.)
c      call agsetf('X/LOGARITHMIC.',1.)
c      do 98 i=1,10
c      ta=1.5+i*0.1
c      ta_array(1)=2.3
      do 97 j=1,1
c         if(j==2) call colora('green')
c         if(j==3) call colora('red')  
      call trace(log10(eta_array),log(gamma_array(:,j)) ,size_plot,-1,-1,0.,0.)
 97   continue
c 98   continue 
      call frame(0)
c.....finalize plot file
      call plote

c$$$c.....initialize plotting
c$$$      call ncarcgm(1,'gamma_dless_fig7.cgm')
c$$$c....make plot
c$$$      call mapsll(eta_array(1),eta_array(size_plot),0.1,
c$$$     A 0.2,.1,1.,.1,1.)
c$$$c      call mapg(eta_array(1),eta_array(size_plot),5.e-2,6.e-2,.1,.5,.7,1
c$$$c     1 .)
c$$$      call trace(eta_array,gamma_dless_array,size_plot,-1,-1,0.,0.)
c$$$c      call setlch(xary(1)-1.0,yary(1)-0.1,0,2,1,-1)
c$$$c      write(s100,1002)
c$$$c      call gtext(s100,80,0)
c$$$c 1002 format(" test plot")^
c$$$c
c$$$c
c$$$      call frame(0)
c$$$c.....finalize plot file
c$$$      call plote

c 34   continue
c 33   continue
      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccc   function det     cccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function det(gamma)
      implicit none
      real det,gamma
      real d_right,d_left
      complex de,do
      real der,dor
      real E,F,G,H,K0
      real epsd
      real rmatch
      integer ifail,n
      real eta,rc,dqc,Bpc,Lr_DI,Lr
      common /delta_block/ d_right,d_left,E,F,G,H,K0,Lr_DI,Lr
      if (gamma.le.0) then
         write(*,*) 'negative gamma =',gamma
c         stop
      end if

      rmatch=8
c      epsd=1.e-3
      ifail=0
      write(*,*) '------------------------------------'
       write(*,*) 'gamma_iter =',gamma
      call innerc(cmplx(gamma),E,F,G,H,K0,de,do,ifail,epsd
     1 ,rmatch)
      if (ifail.ne.0)then
         write(*,*) 'innerc did not converge','ifail = ',ifail
         det=0.
         return
c         stop
         
      end if
      der=real(de*Lr_DI)
      dor=real(do*Lr_DI)
c      der=real(de)
c      dor=real(do)
c      write(*,*) 'Lr_DI=',Lr_DI
c      write(*,*) 'der =',der
c      write(*,*) 'dor =',dor

      det=-(d_right-der)*(d_left-dor)-(d_right-dor)*(d_left-der)
      write(*,*) 'det=',det
c      write(*,*) 'gamma_iter =',gamma
      end function


      subroutine testode()
      real D02EJW,D02EJX,D02EJY
      external D02EJW,D02EJX,D02EJY
      integer ifail,iw
      real, dimension(100)::W2 
      real,dimension(1):: y
      real x,xend,tol
      iw=100
      y=0.
      ifail=0
      x=0.
      xend=2.
      tol=1.e-4
      call D02EJF(x,xend,1, y,odefun,D02EJY,tol,    'D', D02EJX,D02EJW,W2,iw,ifail)
      write(*,*) ifail
      write(*,*) x
      write(*,*) y
      end subroutine

      subroutine odefun(x,yin,yout)
      real x
      real,dimension(1)::yin,yout
      yout=x
      write(*,*) x
      end subroutine

      end module modtest 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccc   main program     cccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program moddemo
      
      use modtest
      integer size
      common /size_block/ size
      size=10000
      call main_fig7(size,100)
c     call testode()
      end program moddemo
