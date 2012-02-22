c
      subroutine spl1d1(n,x,f,w,iop,ij,a,b,c)
c     where n= number of points in the interpolation
c           x= origin of table of independent variable
c           f= origin of table of dependent variable
c           w= an array of dimension n which contains the calculated
c              second derivatives upon return
c           iop= an array of dimension 2 which contains combinations of
c                the integers 1 thru 5 used to specify the boundary
c                conditions
c           ij= spacing in the f and w tables
c           a,b,c= arrays of dimension n used for temporary storage
c
c$$$      dimension iop(*),x(*),f(*),w(*),a(*),b(*),c(*), comm(6)
c    comm made real*8 in order to hold 8 characters in ibm-land
c$$$      data  comm          /8hspl1d1 n,8h less th,8han 4. re,8hsults in,
c$$$     18hcorrect.,8h        /

c...Modernize  comm(6) in above:
      DIMENSION iop(*),x(*),f(*),w(*),a(*),b(*),c(*)
      CHARACTER(LEN=50) :: comm
      DATA  comm / "spl1d1 n less than 4. results incorrect." /

      data zz,oz,tz,sz/0.0e0,1.0e0,3.0e0,6.0e0/
c.. Initialize w = 0.0, since it's not.  Chance 09/28/2008      
      w(1:n) = 0.0
      k=n-1
      a(2)=-(x(2)-x(1))/sz
      b(2)=(x(3)-x(1))/tz
      w(ij+1)=(f(2*ij+1)-f(ij+1))/(x(3)-x(2))-(f(ij+1)-f(1))
     1/(x(2)-x(1))
      if (n-3)3,4,3
    3 do 10 i=3,k
      m=(i-1)*ij+1
      j1=m+ij
      j2=m-ij
      con=(x(i+1)-x(i-1))/tz
      don=(x(i)-x(i-1))/sz
      b(i)=con-(don**2)/b(i-1)
      e=(f(j1)-f(m))/(x(i+1)-x(i))-(f(m)-f(j2))/
     1(x(i)-x(i-1))
      w(m)=e-(don*w(j2))/b(i-1)
   10 a(i)=-(don*a(i-1))/b(i-1)
    4 k1=(n-2)*ij+1
      c(n-1)=-((x(n)-x(n-1))/sz)/b(n-1)
      w(k1)=w(k1)/b(n-1)
      a(n-1)=a(n-1)/b(n-1)
      k2=k-1
      if (n-3)7,8,7
    7 do 20 i=2,k2
      j=n-i
      con=(x(j+1)-x(j))/sz
      a(j)=(a(j)-con*a(j+1))/b(j)
      c(j)=-(con*c(j+1))/b(j)
      k3=(j-1)*ij+1
      m=k3+ij
   20 w(k3)=(w(k3)-con*w(m))/b(j)
    8 k4=(n-1)*ij+1
      if (iop(1)-5) 201,200,201
  201 c1=w(1)
      if (iop(2)-5) 203,202,203
  203 c2=w(k4)
      go to 205
  200 if (n-4)300,302,302
  302 a1=x(1)-x(2)
      a2=x(1)-x(3)
      a3=x(1)-x(4)
      a4=x(2)-x(3)
      a5=x(2)-x(4)
      a6=x(3)-x(4)
      w(1)=f(1)*(oz/a1+oz/a2+oz/a3)-a2*a3*f(ij+1)/(a1*a4*a5)+
     1 a1*a3*f(2*ij+1)/(a2*a4*a6 )-a1*a2*f(3*ij+1)/(a3*a5*a6)
      go to 201
  202 if (n-4)300,303,303
  303 b1=x(n)-x(n-3)
      b2=x(n)-x(n-2)
      b3=x(n)-x(n-1)
      b4=x(n-1)-x(n-3)
      b5=x(n-1)-x(n-2)
      b6=x(n-2)-x(n-3)
      l1=k4-ij
      l2=l1-ij
      l3=l2-ij
      w(k4)=-b2*b3*f(l3)/(b6*b4*b1)+b1*b3*f(l2)/(b6*b5*b2)
     1 -b1*b2*f(l1)/(b4*b5*b3)+f(k4)*(oz/b1+oz/b2+oz/b3)
      go to 203
c  cdc compiler apparently permits transfer into the range of a do-loop
c   corrected for ibm  compilers (br)
c 205 do 50 i=1,k
 205          i    =    1
 2051 continue
      m=(i-1)*ij+1
      go to 60
   70 if (i-1)80,50,80
   80 w(1)=w(1)-bob*w(m)
      w(k4)=w(k4)-bill*w(m)
      a(1)=a(1)-bob*a(i)
      a(n)=a(n)-bill*a(i)
      c(1)=c(1)-bob*c(i)
      c(n)=c(n)-bill*c(i)
c     see   note at statement label 250
   50 continue
      i=i+1
      if ( i .le. k )   go to 2051
c  50 continue
      go to 100
   60 mk=iop(1)
      go to (62,64,66,68,66),mk
   62 if (i-1)71,63,71
   63 a(1)=-oz
      c(1)=zz
      go to 500
   71 bob=zz
      go to 500
   64 if (i-1)73,76,73
   76 a(1)=-oz
      c(1)=zz
      w(1)=zz
      go to 500
   73 if (i-2)81,81,82
   81 bob=-c1
      go to 500
   82 bob=zz
      go to 500
   66 if (i-1)83,84,83
   84 a(1)=-(x(2)-x(1))/tz
      c(1)=zz
      w(1)=-c1+(f(ij+1)-f(1))/(x(2)-x(1))
      go to 500
   83 if (i-2)85,85,86
   85 bob=(x(2)-x(1))/sz
      go to 500
   86 bob=zz
      go to 500
   68 if (i-1)87,88,87
   88 a(1)=-oz
      c(1)=oz
      w(1)=zz
      go to 500
   87 bob=zz
  500 ml=iop(2)
      go to (120,130,140,150,140),ml
  120 if (i-1)121,122,121
  122 a(n)=zz
      c(n)=-oz
      go to 70
  121 bill=zz
      go to 70
  130 if (i-1)131,132,131
  132 a(n)=zz
      c(n)=-oz
      w(k4)=zz
      go to 70
  131 if (i-k)134,133,134
  133 bill=-c2
      go to 70
  134 bill=zz
      go to 70
  140 if (i-1)141,142,141
  142 a(n)=zz
      c(n)=(x(n-1)-x(n))/tz
      w(k4)=c2-(f(k4)-f(k1))/(x(n)-x(n-1))
      go to 70
  141 if (i-k)143,144,143
  144 bill=(x(n)-x(n-1))/sz
      go to 70
  143 bill=zz
      go to 70
  150 if (i-1)151,152,151
  152 a(n)=zz
      c(n)=(x(n-1)+x(1)-x(n)-x(2))/tz
      w(k4)=(f(ij+1)-f(1))/(x(2)-x(1))-(f(k4)-f(k1))/(x(n)-x(n-1))
      go to 70
  151 if (i-2)153,154,153
  154 bill=(x(2)-x(1))/sz
      go to 70
  153 if (i-k)155,156,155
  156 bill=(x(n)-x(n-1))/sz
      go to 70
  155 bill=zz
      go to 70
  100 con=a(1)*c(n)-c(1)*a(n)
      d1=-w(1)
      d2=-w(k4)
      w(1)=(d1*c(n)-c(1)*d2)/con
      w(k4)=(a(1)*d2-d1*a(n))/con
      do 110 i=2,k
      m=(i-1)*ij+1
  110 w(m)=w(m)+a(i)*w(1)+c(i)*w(k4)
      go to 305
  300 call labrt(1,comm,1)
  305 return
      end
c
c................................................................
c
      subroutine spl1d2(n,x,f,w,ij,y,tab)
      data wz,sz/2.0e0,6.0e0/
c     where n= number of points in the interpolation
c           x= origin of table of the independent variable
c           f= origin of table of the dependent variable
c           w= origin of table of second derivatives as calculated by
c              spl1d1
c           ij= spacing in the tables f and w
c           y= the point at which interpolation is desired
c           tab= an array of dimension 3 which contains the function
c                value, first derivative, and second derivative at y
c
      dimension x(*),f(*),w(*),tab(3)
c
      mflag = 0
c
c     locate y in the x table
c
      if(y-x(1))10,10,20
   10 i=1
      go to 30
   20 if(y-x(n))15,40,40
   40 i=n-1
      go to 30
   15 call search(y,x,n,i,mflag)
   30 mi=(i-1)*ij+1
      k1=mi+ij
      flk=x(i+1)-x(i)
c
c     calculate f(y)
c
      a=(w(mi)*(x(i+1)-y)**3+w(k1)*(y-x(i))**3)/(sz*flk)
      b=(f(k1)/flk-w(k1)*flk/sz)*(y-x(i))
      c=(f(mi)/flk-flk*w(mi)/sz)*(x(i+1)-y)
      tab(1)=a+b+c
c
c     calculate the first derivative at y
c
      a=(w(k1)*(y-x(i))**2-w(mi)*(x(i+1)-y)**2)/(wz*flk)
      b=(f(k1)-f(mi))/flk
      c=flk*(w(mi)-w(k1))/sz
      tab(2)=a+b+c
c
c     calculate the second derivative at y
c
      tab(3)=(w(mi)*(x(i+1)-y)+w(k1)*(y-x(i)))/flk
      return
      end
c
c...........................................................
c
      subroutine labrt(isw,lhol,inx)
c      dimension lhol(8)
      character(*) :: lhol
      logical ps, ts
      data np/10/,ps/.true./,ts/.false./
      if((isw.eq.0).or.(isw.gt.5))return
      go to ( 1,2,3,4,5 ), isw
    1 if ( ps .and. (np .gt. 0) )   write ( 3, 27 )   lhol, inx
   27 format(1h0,9x,a,3x,z4)
      np=np-1
c      if ( ts )   call exit
      if ( ts )   stop
      return
    2 ps=.false.
      return
    3 ps=.true.
      np=inx
      return
    4 ts=.true.
      return
    5 ts=.false.
      return
      end
c
c..............................................................
c
c.................................................................
      subroutine search(xbar,x,n,i,mflag)
c.................................................................
c  double precision version
c   the array x is really a double precision real array
c   indexing changed
c      dimension com1(8)
      CHARACTER(132) ::  com1
c     integer xbar,x(1)
      dimension x(*)          !!  Changed from 1 to * MAy 9 2011. MSC
c$$$      data com1/ 'search  ','xbar is ','outside ','range of',' table',
c$$$     $     " ", " ", " "/
      DATA com1 / "Search: xbar is otside range of table. " /
      mflag=0
      i=n
      if(xbar.eq.x(n))return
      i=1
      if(n.le.1) return
      ixbar=sign(1.0,xbar)
      ix1=sign(1.0,x(1))
      ixn=sign(1.0,x(n))
      do 5 k=1,n
      j=i+i
      if(j.ge.n) go to 6
  5   i=j
 6    k=i
      mflag = 1
      do 115  l=2,n
      a=x(l-1)
      b=x(l)
      if(sign(1.0,a)-sign(1.0,b)) 7,113,8
 113   if(a-b)7,115,8
 115    continue
  7    j=1
      if(ixbar.lt.ix1.or.(ixbar.eq.ix1.and.xbar.lt.x(1)).or.ixbar
     1.gt.ixn.or.(ixbar.eq.ixn.and.xbar.gt.x(n))) go to 16
      go to 10
  8    j=2
      if(ixbar.lt.ixn.or.(ixbar.eq.ixn.and.xbar.lt.x(n)).or.ixbar
     1.gt.ix1.or.(ixbar.eq.ix1.and.xbar.gt.x(1))) go to 16
   10 k=k/2
       a=x(i)
      go to (11,20),j
  11   if(ixbar-sign(1.0,a)) 111,1111,2111
 1111 if(xbar-a)111,14,2111
 2111 b=x(i+1)
      if(ixbar-sign(1.0,b)) 2112,2113,12
 2113   if(xbar.ge.b) go to 12
 2112   return
 111  i = i-k
      go to 13
   12 i = i+k
   13 if (i.lt.n) go to 10
      k=k/2
      go to 111
   14 mflag=0
      return
   16 call labrt(1,com1,1)
      mflag=2
      return
  20   if(ixbar-sign(1.0,a) ) 2120,2121,111
 2121  if(xbar-a) 2120,14,111
 2120 b=x(i+1)
       if(ixbar-sign(1.0,b)) 12,2122,2112
 2122  if(xbar-b) 12,12,2112
 8888    return
      end
c............................................................
c
      subroutine searchx(xbar,x,n,i,mflag)
c  double precision version
c   the array x is really a double precision real array
c   indexing changed
      CHARACTER(132) ::  com1
c      dimension com1(8)
      integer xbar,x(*)           !!  Changed from 1 to * MAy 9 2011. MSC
      if(.false.) go to 8888
      DATA com1 / "Searchx: xbar is otside range of table. " /
c$$$
c$$$      data com1/ "search  ","xbar is ","outside ","range of"," table",
c$$$     $     " ", " ", " "/
      mflag=0
      i=n
      if(xbar.eq.x(n))return
      i=1
      if(n.le.1) return
      ixbar=isign(1,xbar)
      ix1=isign(1,x(1))
      ixn=isign(1,x(n))
      do 5 k=1,n
      j=i+i
      if(j.ge.n) go to 6
  5   i=j
 6    k=i
      mflag = 1
      do 115  l=2,n
      ia=x(l-1)
      ib=x(l)
      if(isign(1,ia)-isign(1,ib)) 7,113,8
 113   if(ia-ib)7,115,8
 115    continue
  7    j=1
      if(ixbar.lt.ix1.or.(ixbar.eq.ix1.and.xbar.lt.x(1)).or.ixbar
     1.gt.ixn.or.(ixbar.eq.ixn.and.xbar.gt.x(n))) go to 16
      go to 10
  8    j=2
      if(ixbar.lt.ixn.or.(ixbar.eq.ixn.and.xbar.lt.x(n)).or.ixbar
     1.gt.ix1.or.(ixbar.eq.ix1.and.xbar.gt.x(1))) go to 16
   10 k=k/2
       ia=x(i)
      go to (11,20),j
  11   if(ixbar-isign(1,ia)) 111,1111,2111
 1111 if(xbar-ia)111,14,2111
 2111 ib=x(i+1)
      if(ixbar-isign(1,ib)) 2112,2113,12
 2113   if(xbar.ge.ib) go to 12
 2112   return
 111  i = i-k
      go to 13
   12 i = i+k
   13 if (i.lt.n) go to 10
      k=k/2
      go to 111
   14 mflag=0
      return
   16 call labrt(1,com1,1)
      mflag=2
      return
  20   if(ixbar-isign(1,ia) ) 2120,2121,111
 2121  if(xbar-ia) 2120,14,111
 2120 ib=x(i+1)
       if(ixbar-isign(1,ib)) 12,2122,2112
 2122  if(xbar-ib) 12,12,2112
 8888    return
      end

c............................................
      SUBROUTINE green
c.................................................................

      include 'vacuum1.inc'
      include 'vacuum3.inc'

      COMMON / grvac / aval0, phex0, zphfc0, iscs, isco,
     $     j1com,j2com

c.....Bval not divided by twopi like aval is.  Take care ov this
c     in subroutine kernel. 6-27-94

      pm = 0.
      pn = 0.
      pp = 0.
      aleg0 = 0.
      aleg1 = 0.
c     
      sqpi  = sqrt(pye)
      pii  = two / pye
      gam  = sqpi
c     n is floating point in common
      nloc  = n + 0.1e0
      xs2  = xs * xs
      xt2  = xt * xt
      xp2  = xs2 + xt2
      xm2  = xt2 - xs2
      zm2  = ( zt - zs )**2
      r14  = xm2*xm2 + zm2*zm2 + two*xp2*zm2
      r1sq = sqrt( r14 )
      r1 = sqrt( r1sq )
      s  = (xp2 + zm2 )/r1sq   ! This is the standard argument of P(n)
c     
c     use upwards recurrence relations...
     
      call aleg ( s,nloc, pm,pn,pp, aleg0,aleg1 )
c     
c     calculate gamma(half-n)
      kloc=0
      ak=zero
      if ( nloc .eq. 0 )  go to 10
    5 kloc = kloc+1
      ak = float(kloc)
      ak02 = half-ak
      gam = gam / ak02
      if ( kloc .ne. nloc )  go to 5
c     
c     now pp contains p(n+1) and pn contains p(n)
c     
 10   gg  = -two * sqpi * gam / r1
      bval  = -gg*pn
      aval1 = ( n*(xs2+xt2+zm2)*(xt2-xs2-zm2)+xt2*(xm2+zm2))*pn
      aval2 = two*xt*xs*(xm2-zm2)*pp
      aval3 = ztp*(aval1+aval2) / xt
      aval4 = ( two*n+one)*(xp2+zm2)*pn+four*xt*xs*pp
      aval5 = xtp*(zt-zs)*aval4
      aval6 =(aval3-aval5) / ( xt*r14 )
      aval = - xt2*aval6 * gg / twopi
      aval0 = ztp*(two*xs*(zm2-xm2)*aleg1 - xt*(xm2+zm2)*aleg0)
      aval0 = aval0 + xtp*(zt-zs)*(four*xt*xs*aleg1+(xp2+zm2)*aleg0)
      aval0 = -aval0*xt / (r14*r1)
      return
      end
c     
c.....................................................
      SUBROUTINE aleg(x,nloc,pm,pn,pp, aleg0,aleg1 )
c.....................................................
c     
c     subroutine to calculate half integral legendre functions.
c     uses upwards recurrence relations starting from elliptic
c     integrals evaluated using Bulirsch's algorithm
c     these expressions are very bad for large values of nloc.
c     zkisq is ths the 1 - k**2 in Elliptic integeral parlance.     

c     This modified from the old aleg subroutine to use the 
c     Bulirsch algorithms for the Elliptic functions. 
c     The new integral representation of the Legendre function is used
c     here for n*rhohat >= 0.1

c     Reference: JCP 221 (2007) 330-348

      PARAMETER ( pye=3.141592653589793, pii=2.0/pye, sqpi=SQRT(pye),
     $            sqtwo=SQRT(2.0), half=0.5 )

c...  Sum of ak_i = pi/2. Sum of ae_i = pi/2 - 1.0

!::::::::::::::::::::::::::::::::::::::::::::::::::::::

!     This stuff for Gaussian Itegration:

c$$$      REAL, DIMENSION(8):: tgaus, wgaus 
      REAL, DIMENSION(32):: tg32, wg32, xg32
      REAL, DIMENSION(5):: xu, xl
c$$$      REAL, DIMENSION(10):: cfac, wksp

!.... Weights and abscissae for 32 points gaussian quadrature.

      wg32(1)  =  0.007018610009470096600 
      wg32(2)  =  0.016274394730905670605
      wg32(3)  =  0.025392065309262059456
      wg32(4)  =  0.034273862913021433103
      wg32(5)  =  0.042835898022226680657
      wg32(6)  =  0.050998059262376176196
      wg32(7)  =  0.058684093478535547145
      wg32(8)  =  0.065822222776361846838
      wg32(9)  =  0.072345794108848506225
      wg32(10) =  0.078193895787070306472
      wg32(11) =  0.083311924226946755222
      wg32(12) =  0.087652093004403811143
      wg32(13) =  0.091173878695763884713
      wg32(14) =  0.093844399080804565639
      wg32(15) =  0.095638720079274859419
      wg32(16) =  0.096540088514727800567
      
      DO i = 1, 16
         wg32(16+i) = wg32(17-i)
      END DO

      xg32(1:16) = (/ -0.997263861849481563545,
     $ 0.985611511545268335400, 
     $ 0.964762255587506430774,
     $ 0.934906075937739689171, 
     $ 0.896321155766052123965, 
     $ 0.849367613732569970134, 
     $ 0.794483795967942406963, 
     $ 0.732182118740289680387, 
     $ 0.663044266930215200975, 
     $ 0.587715757240762329041, 
     $ 0.506899908932229390024, 
     $ 0.421351276130635345364, 
     $ 0.331868602282127649780, 
     $ 0.239287362252137074545, 
     $ 0.144471961582796493485, 
     $ 0.048307665687738316235 /)

!     xg32(17:32) = - (/ xg32(16:1) /)

      DO i = 1, 16
         xg32(16+i) = - xg32(17-i)
      END DO

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c$$$      pii = 2.0 / pye
c$$$      sqpi = SQRT ( pye )
c$$$      sqtwo = SQRT(2.0)
c$$$c      init = init + 1
c$$$      half = 0.5

      gam = sqpi
      xxq = x*x
      ysq = xxq - 1.0
      y = SQRT( ysq )
      w = x+y

      rhohatsq = 1.0 / ( 2.0 * y*w )
      rhohat = SQRT (rhohatsq)

      zk1i = w              
      zk1 = 1.0/w           ! This is k1 = SQRT(1-k**2) = SQRT(m_1)
      zk1sq = zk1**2        ! This is m_1
      zk1sqrt = SQRT(zk1)   ! This is m_1^(1/4)
      zk1sqrti = SQRT(zk1i) ! This is m_1^(-1/4)

      errbu = 1.0e-8
      CALL ek3 ( zk1sq, ierbu, errbu, 10, elipk, elipe, convbu, kcbu )

      pn = pii * zk1sqrt * elipk
      pnp = pii * zk1sqrti * elipe

      aleg0 = pn

      pp = (  pnp - x*pn ) / (2.0*y)

      aleg1 = pp

c... Use Gaussian Integration if ...

      IF ( nloc*rhohat >= 0.1 ) GO TO 100

c      GO TO 100
c 85    CONTINUE

      kloc=0
      ak = 0.0

      IF ( nloc == 0 )  GO TO 10

    5 kloc=kloc+1
      ak = FLOAT(kloc)
      ak02 = 0.5 - ak
      pm = pn
      pn = pp
      pp = -2.0*ak*x*pn/y - ak02*ak02*pm
      gam = gam / ak02
      IF ( kloc /= nloc )  GO TO 5

 10   CONTINUE

      GO TO 500

 100  CONTINUE

c...  use Gauss integration of the new integral representation 
c     if n*rhohat >= 0.1
c...  The integration is done in nng segments [xl(ing),xu(ing)]. 
c     Each stored in gint.

      ngauss = 32
      nng = 1
      xl(1) = 0.0
      xu(1) = 5.0

      gint = 0.0
      gintp = 0.0

      DO 165 ing = 1, nng

!.....xl, xu are the lower and upper limits of the gaussian integration
!     The integration is done in nng sections
!     This will calculate P(n) and P(n+1) together. 
!        variables for P(n+1) will usually have p appended.

         agaus = half*( xu(ing)+xl(ing) )
         bgaus = half*( xu(ing)-xl(ing) )
         
         tg32(1:32) = agaus + xg32(1:32) * bgaus

         ginti = 0.0
         gintip = 0.0

         DO ig = 1, ngauss
            tg0 = tg32(ig)
            tg02 = tg0**2
            tg1  = tg02 / (2.0*nloc)
            tg1p = tg02 / (2.0*nloc+2.0)
            sinhtg1  = SINH(tg1)
            sinhtg1p = SINH(tg1p)
            sinhtg12  = sinhtg1  * sinhtg1
            sinhtg12p = sinhtg1p * sinhtg1p
            dnom  = x * sinhtg12  + sinhtg1 *SQRT(1.0 + sinhtg12)
            dnomp = x * sinhtg12p + sinhtg1p*SQRT(1.0 + sinhtg12p)
            dnom  = SQRT(dnom)
            dnomp = SQRT(dnomp)
            anumr = tg0 * EXP(-tg02)
            ginti  = ginti  + wg32(ig)*anumr / dnom
            gintip = gintip + wg32(ig)*anumr / dnomp
         END DO                 ! 32 point Gaussian
         
         ginti  = bgaus * ginti
         gintip = bgaus * gintip
         gint  = gint  + ginti
         gintp = gintp + gintip

 165  CONTINUE                  !  Gaussian integration segments

c... Now calculate the coeficients for the Legendre functions.

      pcoef = SQRT ( (x-1.0)/(x+1.0) )
      twopi = 2.0 * pye

c.. gamn is  Gamma[1/2-n]
c   gamp is  Gamma[1/2-(n+1)]

      gamn = sqpi
      gamp = - 2.0 * sqpi

      IF ( nloc /= 0 ) THEN

         gamn = sqpi /
     $        PRODUCT( (/ ( -(i-1)-0.5, i = 1, nloc ) /) )
         gamp = - gamn / (nloc+0.5)

c$$$         nwrt = 1
c$$$         IF ( nwrt == 1 ) THEN
c$$$         WRITE (6, '("nloc, gamn, gamp = ", i3, 2es12.4)' )
c$$$     $        nloc, gamn, gamp
c$$$         nwrt = nwrt + 1
c$$$         END IF
         
      END IF

      gint  = sqtwo * pcoef**nloc * gint / (nloc*sqpi*gamn)
      gintp = sqtwo * pcoef**(nloc+1) * gintp / ((nloc+1.0)*sqpi*gamp)
      pn = gint  ! P(n)
      pp = gintp  ! P(n+1)

 500  CONTINUE
c$$$
c$$$         nwrt = 1
c$$$         IF ( nwrt == 1 ) THEN
c$$$            WRITE (23, '("nloc, x, rhohat, pn, pp = ", i3, 4es12.4)' )
c$$$     $        nloc, x, rhohat, pn, pp
c$$$         nwrt = nwrt + 1
c$$$         END IF

      RETURN
      END
c     
!...................................................
       SUBROUTINE ek3(eta,ier,error,maxit,cel1,cel2,convg, kounter)
!..................................................

!  Compute the complete elliptic integral of first and second kind
!      cel(kc,p,a,b).  
!  Bulirsch's method. Numerical Recipes, modified by Turnbull to 
!    calculate both K and E simultaneously.

!  Returns cel1 = K, cel2 = E.
!  Precision is error**2, 

!  eta, the complementary parameter, (1 - k^2), is the square of 
!          the argument kc
!  p   is 1
!  a   is 1
!  b   is 1 for the first kind and b is eta( = kc**2) for the second kind

       PARAMETER (pi=3.1415926535897932385 , pi2 = pi/2.0)


       pp     = 1.0
       aa     = 1.0
       bb1    = 1.0
       bb2    = ABS(eta)

       ier    = 0
       IF(eta .LE. 0.0  .OR.  eta > 1.0) THEN
          IF(eta < 0.0) ier   = 1
          IF(eta == 0.0) ier   = 2
          IF(eta > 1.0) ier   = 3
          cel1  = 0.0
          cel2  = 0.0
          RETURN
       END IF


       qcval  = SQRT(ABS(eta))
       aval0  = aa
       bval1  = bb1
       bval2  = bb2
       pval0  = pp

       eval   = qcval
       emval  = 1.0


       IF(pval0 > 0.0) THEN
          pval  = SQRT(pval0)
          aval1 = aval0
          aval2 = aval0
          bval1 = bval1/pval
          bval2 = bval2/pval

       else
          fval  = qcval*qcval
          tval  = 1.0  - fval
          gval  = 1.0  - pval0
          fval  = fval - pval0
          qval1 = tval*(bval1 - aval0*pval0)
          qval2 = tval*(bval2 - aval0*pval0)

          pval  = SQRT(fval/gval)
          aval1 = (aval0 - bval1) / gval
          aval2 = (aval0 - bval2) / gval
          bval1 =  aval1*pval - qval1/(gval*gval*pval)
          bval2 =  aval2*pval - qval2/(gval*gval*pval)
       END IF


       kounter = 0
 100   CONTINUE
       kounter = kounter + 1

       hval1  = aval1
       hval2  = aval2
       aval1  = aval1 + bval1/pval
       aval2  = aval2 + bval2/pval
       rval   = eval/pval
       bval1  = bval1 + hval1*rval
       bval1  = bval1 + bval1
       bval2  = bval2 + hval2*rval
       bval2  = bval2 + bval2
       pval   = rval + pval

       sval   = emval
       emval  = qcval + emval

       IF (ABS(sval-qcval) > sval*error) THEN
          qcval  = SQRT(eval)
          qcval  = qcval + qcval
          eval   = qcval*emval
          GO TO 100
       END IF

       IF(sval /= 0.0) snorm = sval*sval
       IF(sval == 0.0) snorm = 1.0
       convg  = (sval-qcval)*(sval-qcval) / snorm
       IF ( convg <= 1.0e-100 ) convg = 1.0e-100
       cnvlog = ALOG10(ABS(convg))
       logcnv = IFIX(cnvlog)

       IF (kounter > maxit) THEN
          IF(logcnv < 0) ier   = logcnv
          IF(logcnv >= 0) ier   = -1
          cel1  = pi2*(bval1 + aval1*emval) / (emval*(emval+pval))
          cel2  = pi2*(bval2 + aval2*emval) / (emval*(emval+pval))
          RETURN
       END IF


       cel1  = pi2*(bval1 + aval1*emval) / (emval*(emval+pval))
       cel2  = pi2*(bval2 + aval2*emval) / (emval*(emval+pval))

       RETURN
       END


c....................................................
c     
      subroutine alegold(x,nloc,pm,pn,pp, aleg0,aleg1 )
c     
c     subroutine to calculate half integral legendre functions.
c     uses upwards recurrence relations starting from elliptic
c     integrals,evaluated from hasting's formulae.
c     these expressions are very bad for large values of nloc.
c     
      data pye/3.141592653589793/, init/0/,
     .     ak0 /1.38629436112/,
     .     ak1 /0.09666344259/,
     .     ak2 /0.03590092383/,
     .     ak3 /0.03742563713/,
     .     ak4 /0.01451196212/,
     .     bk0 /0.5/,
     .     bk1 /0.12498593597/,
     .     bk2 /0.06880248576/,
     .     bk3 /0.03328355346/,
     .     bk4 /0.00441787012/,
     .     ae1 /0.44325141463/,
     .     ae2 /0.0626060122/,
     .     ae3 /0.04757383546/,
     .     ae4 /0.01736506451/,
     .     be1 /0.2499836831/,
     .     be2 /0.09200180037/,
     .     be3 /0.04069697526/,
     .     be4 /0.00526449639/
c     
c... Sum of ak_i = pi/2. Sum of ae_i = pi/2 - 1.0
c
c$$$      if ( init .gt. 0 )  go to 110
      pii = 2.0 / pye
      sqpi = sqrt ( pye )
      init = init + 1
c$$$ 110  continue
c     
      gam = sqpi
      xxq = x*x
      s = (x-1.0)/(x+1.0)
      ysq = xxq-1.0
      y = sqrt(ysq)
      w = x+y
      v = 2.0*y/w
      x1 = 2.0 / (x+1.0)
      x2 = x1*x1
      x3 = x2*x1
      x4 = x3*x1
c     
      elipk = ak0+ak1*x1+ak2*x2+ak3*x3+ak4*x4
     .     - (bk0+bk1*x1+bk2*x2+bk3*x3+bk4*x4)*alog(x1)
c     
      pn = pii*sqrt(2.0/(x+1.0))*elipk
c     
      aleg0 = pn
c     
      x1 = 1.0 / w**2
      x2=x1*x1
      x3 = x2*x1
      x4 = x3*x1
c     
      elipe=1.0
      if(abs(x1) .gt. 1.0e-6)
     .     elipe=1.0+ae1*x1+ae2*x2+ae3*x3+ae4*x4
     .     - (be1*x1+be2*x2+be3*x3+be4*x4)*alog(x1)
c     
      pp = (pii*sqrt(w)*elipe-x*pn)/(2.0*y)
c     
      aleg1 = pp
c     
      kloc=0
      ak = 0.0
      if ( nloc .eq. 0 )  go to 10
    5 kloc=kloc+1
      ak = float(kloc)
      ak02 = 0.5 - ak
      pm = pn
      pn = pp
      pp = -2.0*ak*x*pn/y - ak02*ak02*pm
      gam = gam / ak02
      if ( kloc .ne. nloc )  go to 5
 10   continue
c     
      RETURN
      END

c...................................................................
c     
      SUBROUTINE trans ( vecin,mthin, vecout,mth )
c..................................................................     
      DIMENSION vecin(*), vecout(*)
c     
      vecin(mthin+1) = vecin(1)
      vecin(mthin+2) = vecin(2)
c     
      if ( mth .ne. mthin ) go to 20
c     
      do 10 i = 1, mth + 2
         vecout(i) = vecin(i)
 10   continue
c     
      return
c     
 20   continue
c     
      do 50 i = 1, mth
         ai = i-1
         x = ai / mth
         iop = 1
         call lagpe4 ( vecin,mthin, x, vecout(i), df, 1 )
 50   continue
c     
      mth1 = mth + 1
      mth2 = mth1 + 1
      vecout(mth1) = vecout(1)
      vecout(mth2) = vecout(2)
c     
      RETURN
      END
c     
c...................................................................
      SUBROUTINE transdx ( vecin,mthin, vecout,mth, dx0,dx1 )
c...................................................................
c
c..Starts new origin of data at dx1/mth + dx0/mthin and Interpolates from 
c   mthin to mth points.      VECIN preserved. 
c... Note. Using lagpe4 with abscissa's range = [0,1]
c     
      REAL, DIMENSION(*) ::  vecin, vecout
c     
      vecin(mthin+1) = vecin(1)
      vecin(mthin+2) = vecin(2)
c     
      if ( (mth .ne. mthin) .or. (abs(dx0) .gt. 1.e-6) ) go to 20
c     
      do 10 i = 1, mth + 2
         vecout(i) = vecin(i)
 10   continue
c     
      return
c     
 20   continue
c     
      do 50 i = 1, mth
         ai = i-1 + dx1
         x = ai / mth + dx0 / mthin
         iop = 1
         call lagpe4 ( vecin,mthin, x, vecout(i), df, 1 )
 50   continue
c     
      mth1 = mth + 1
      mth2 = mth1 + 1
      vecout(mth1) = vecout(1)
      vecout(mth2) = vecout(2)
c     
      RETURN
      END
c     
c...................................................................
      SUBROUTINE transdxc ( vecin,mthin, vecout,mth, dx0,dx1 )
c...................................................................
c
c..Starts new origin of data at dx1/mth + dx0/mthin and Interpolates from 
c   mthin to mth points.      VECIN preserved. 
c... Note. Using lagpe4 with abscissa's range = [0,1]
c     
      COMPLEX, DIMENSION(*) ::  vecin, vecout
      COMPLEX :: df
c     
      vecin(mthin+1) = vecin(1)
      vecin(mthin+2) = vecin(2)
c     
      IF ( (mth /= mthin) .OR. (ABS(dx0) > 1.e-6) ) GO TO 20
c     
      DO 10 i = 1, mth + 2
         vecout(i) = vecin(i)
 10   CONTINUE
c     
      RETURN
c     
 20   CONTINUE
c     
      DO 50 i = 1, mth
         ai = i-1 + dx1
         x = ai / mth + dx0 / mthin
         iop = 1
         CALL lagpe4c ( vecin,mthin, x, vecout(i), df, 1 )
 50   CONTINUE
c     
      mth1 = mth + 1
      mth2 = mth1 + 1
      vecout(mth1) = vecout(1)
      vecout(mth2) = vecout(2)
C     
      RETURN
      END

c...................................................................
      SUBROUTINE transdxx ( vecin,mthin, vecout,mth, dx0 )
c...................................................................
c
c..Shifts data by dx0/mthin and Interpolates from mthin to mth points.
c       VECIN overwritten and VECOUT used for temporary strorage.
c... Note. Using lagpe4 with abscissa's range = [0,1]

      dimension vecin(*), vecout(*)
c     
      vecin(mthin+1) = vecin(1)
      vecin(mthin+2) = vecin(2)
c     
      if ( (mth .ne. mthin) .or. (abs(dx0) .gt. 1.e-6) ) go to 20
c     
      do i = 1, mth + 2
         vecout(i) = vecin(i)
      end do
c     
      return
c     
 20   continue
c     
      do i = 1, mth
         ai = i-1
         x = ai / mth + dx0 / mthin
         iop = 1
         call lagpe4 ( vecin,mthin, x, vecout(i), df, 1 )
      end do
c     
      mth1 = mth + 1
      mth2 = mth1 + 1
      vecout(mth1) = vecout(1)
      vecout(mth2) = vecout(2)
c
      do i = 1, mth2
         vecin(i) = vecout(i)
      end do
c     
      RETURN
      END
      
c.....................................................
      subroutine smooth0 ( g, n )
c....................................................
c     
      dimension g(*)
c     
      nm = n - 1
c     
      do 1 j = 1, nm
    1 g(j) = g(j) + g(j+1)
c     
      do 2 j = 1, nm
    2 g(n+1-j) = g(n+1-j) + g(n-j)
c     
      g(1) = 2.0 * g(1)
      g(n) = 2.0 * g(n)
c     
      do 3 j = 1, n
    3 g(j) = g(j) / 4.0
c     
      return
      end
c     

c........................................................
      SUBROUTINE smooth1 ( ain, aout, n )
c.........................................................
c     
      REAL, DIMENSION (n) :: zain, za
      DIMENSION ain(*), aout(*)
      
      nm = n - 1

c     Save ain from being overwritten. Use zain for working.

      zain(1:n) = ain(1:n)

      DO j = 1, nm
         za(j) = zain(j)
         zain(j) = zain(j) + zain(j+1)
      END DO

      za(n) = zain(n)
    
      DO j = 1, nm
         za(n+1-j) = za(n+1-j) + za(n-j)
      END DO
     
c......periodicity....

      za(1) = za(n)
      zain(n) = zain(1)
     
      DO j = 1, n
         aout(j) = ( zain(j)+za(j) ) / 4.0
      END DO
    
      END SUBROUTINE smooth1
     
c.....................................................
      subroutine smooth ( g, n )
c.........................................................
c     
      REAL, DIMENSION (n) :: h
      DIMENSION g(*)
c     
      nm = n - 1
c     
      do 1 j = 1, nm
         h(j) = g(j)
    1 g(j) = g(j) + g(j+1)
      h(n) = g(n)
c     
      do 2 j = 1, nm
    2 h(n+1-j) = h(n+1-j) + h(n-j)
c     
c......periodicity....
      h(1) = h(n)
      g(n) = g(1)
c     
      do 3 j = 1, n
    3 g(j) = ( g(j)+h(j) ) / 4.0
c     
      return
      end
c     
c.......................................................
      SUBROUTINE msc_smooths1y ( fin, fout, npts, ns )
c.......................................................

c... Smooth n array, fin, using smooth1 (which smooths in both direction).
c    Does this ns times. npts is no. of points on open interval.

      DIMENSION fin(*), fout(*)
      REAL, DIMENSION(:), ALLOCATABLE :: zfsm

      ALLOCATE ( zfsm(npts+5) )

      zfsm(1:npts+1)  = fin (1:npts+1)

      DO n = 1, ns
         CALL smooth1(zfsm, fout ,npts+1) 
         zfsm(1:npts+1) = fout(1:npts+1)
      END DO

      DEALLOCATE ( zfsm )

      END SUBROUTINE msc_smooths1y

c.................................................

      subroutine lagp ( ax, af, m, nl, x, f, df, iop, iper )
c...................................................
c     
c......written by m.s. chance           . .................
c     
      dimension ax(*), af(*)
c     
c.....ax is array of independent variable, increasing with index.
c     af is array of dependent variable.
c     m is significant dimension of ax and af.
c     nl is desired order for the interpolation.
c     x is the value of x at which f is desired.
c     iop = 0    value of f alone.
c     iop = 1    value 0f df alone.
c     iop = 2    value of f and df.
c     
      zero = 0.0
      one = 1.0
      dax = ax(m) - ax(1)
c     
c     calculate index range of the interpolation.
c     
      in = 1
      do 20 i = 1, m
         if ( ax(i) .gt. x ) go to 25
 20   continue
 25   continue
      in = i - 1
c     
      if ( (in .eq. 0) .or. (in .eq. m) ) go to 30
      if ( ax(in+1)-x .lt. x-ax(in) ) in = in + 1
c     
 30   continue
c     
c     interpolation will be centered at in except near boundaries of
c     array.  the left and right index are respectively,
c     
      inmm = (nl-0.1) / 2.0
      inpp = (nl+0.1) / 2.0
      nll = in - inmm
      nlr = in + inpp
c     
c.....even nl can be off, so adjust...
c     
      if ( ( (nl/2)*2 .eq. nl) .and. ( ax(in) .gt. x ) ) then
         nll = nll - 1
         nlr = nlr - 1
      end if
c     
      if ( (nll.ge.1) .and. (nlr.le.m) ) go to 35
c     
      if ( iper .eq. 1 ) go to 35
c     
      if ( nlr .le. m ) go to 33
      nlr = m
      nll = nlr - nl + 1
      go to 35
c     
 33   continue
      nll = 1
      nlr = nl
c     
 35   continue
c     
c     now interpolate.
c     
      if ( iop .eq. 1 ) go to 120
c     
      f = zero
c     
      do 100 i0 = nll, nlr
c     
         call shft ( i0,i, axi, ax, m, dax )
c     
         alag = one
c     
         do 50 j0 = nll, nlr
c     
            call shft ( j0,j, axj, ax, m, dax )
c     
            if ( i0 .eq. j0 ) go to 50
            alag = alag * ( x-axj ) / ( axi - axj )
 50      continue
c     
         f = f + alag * af(i)
c     
 100  continue
c     
c     end of f.
c     
      if ( iop .eq. 0 ) return
c     
 120  continue
c     
      df = zero
c     
      do 400 i0 = nll, nlr
c     
         call shft ( i0, i, axi, ax, m, dax )
c     
         slag = zero
         do 300 id0 = nll, nlr
c     
            call shft ( id0, id, axid, ax, m, dax )
c     
            if ( id0 .eq. i0 ) go to 300
c     
            alag = one
            do 200 j0 = nll, nlr
c     
               call shft ( j0, j, axj, ax, m, dax )
c     
               if ( j0 .eq. i0 ) go to 200
               if ( j0 .eq. id0 ) go to 160
c     
               alag = alag*( x-axj ) / ( axi-axj )
               go to 200
c     
 160           continue
               alag = alag / ( axi-axid )
c     
 200        continue
            slag = slag + alag
 300     continue
c     
         df = df + slag * af(i)
c     
 400  continue
c     
      return
      end
c
c........................................................
c
         SUBROUTINE shft ( i0,i, axi, ax, m, dax )
c.......................................................

c.. Shifts the index i0 and ax by dax, if i0 falls outside of
c    [1,m] Used in LAGP for periodic functions.
c
         dimension ax(*)
c
         if ( (i0 .ge. 1) .and. (i0 .le. m) ) then
            i = i0
            axi = ax(i)
         else
c               
            if ( i0 .gt. m ) then
               i = mod(i0,m) + 1
               axi = ax(i) + dax
            end if
c
            if ( i0 .lt. 1 ) then
               i = m + i0 - 1
               axi = ax(i) - dax
            end if
         end if
c
         return
         end    
c     
c....................................................................
c     
      SUBROUTINE lagpe4 ( f0,m, x,f, df, iop )
c...................................................................     
c... 4 point lagrange interpolation. Equally spaced abscissa. Abscissa's
c         range is [0,1]

      REAL, DIMENSION(*) :: f0
c     
      h = 1.0 / m
      f0(m+1) = f0(1)
      f0(m+2) = f0(2)
      f0(m+3) = f0(3)
c     
      m0 = x / h
      x0 = m0 * h
      p = (x-x0) / h
      pm1 = p - 1.0
      pm2 = p - 2.0
      pp1 = p + 1.0
      pp2 = p + 2.0
c     
      am1 = - p*pm1*pm2 / 6.0
      a00 =  pm1*pp1*pm2 / 2.0
      ap1 = - p*pp1*pm2 / 2.0
      ap2 =  p*pp1*pm1 / 6.0
c     
      IF ( m0 == 0 ) fm1 = f0(m)
      IF ( m0 >= 1 ) fm1 = f0(m0)
      f00 = f0(m0+1)
      fp1 = f0(m0+2)
      fp2 = f0(m0+3)
c     
      f = am1*fm1 + a00*f00 + ap1*fp1 + ap2*fp2
C     
      RETURN
      END
c....................................................................
c     
      SUBROUTINE lagpe4c ( f0,m, x,f, df, iop )
c...................................................................     
c... 4 point lagrange interpolation. Equally spaced abscissa. Abscissa's
c         range is [0,1]

      COMPLEX, DIMENSION(*) :: f0
      COMPLEX :: f, df
c     
      h = 1.0 / m
      f0(m+1) = f0(1)
      f0(m+2) = f0(2)
      f0(m+3) = f0(3)
c     
      m0 = x / h
      x0 = m0 * h
      p = (x-x0) / h
      pm1 = p - 1.0
      pm2 = p - 2.0
      pp1 = p + 1.0
      pp2 = p + 2.0
c     
      am1 = - p*pm1*pm2 / 6.0
      a00 =  pm1*pp1*pm2 / 2.0
      ap1 = - p*pp1*pm2 / 2.0
      ap2 =  p*pp1*pm1 / 6.0
c     
      IF ( m0 == 0 ) fm1 = f0(m)
      IF ( m0 >= 1 ) fm1 = f0(m0)
      f00 = f0(m0+1)
      fp1 = f0(m0+2)
      fp2 = f0(m0+3)
c     
      f = am1*fm1 + a00*f00 + ap1*fp1 + ap2*fp2
C     
      RETURN
      END
c           
c.........................................................
      subroutine lagpe5 ( f0, fw, m, r, x, f, df, iop, iper )
c.........................................................
c     
c...  Interpolates and find first derivative of a periodic function using 
c     5 point lagrange formula. the range of the function is extended
c     iper points on the left and right into a new work function fw.
c     The independent variable is also defined from x-h*iper, h=grid size.
c     Call intset first to extend range using entry point.                
c                                        .   M.S. Chance
c     
c     f0 = input function on equally spaced parameter.
c     fw = work vector of dimension at least m+2*iper.
c     m = number of points on the open domain.
c     r = range on the closed interval. m+1 
c     iop = 0 function only
c     iop = 1 derivative only
c     iop = 2 both function and derivative.
c     iper = no. of points to extend range. ( >= 3, say )
c     
      dimension f0(*), fw(*), a(5), d(5)
c     
      go to 20
c     
      entry intset ( f0, fw, m, r, iper )
c     
      h = r / m
c     
      do i = 1, m + 1 
         fw(i+iper) = f0(i)
      end do
c     
      do i = 1, iper
         fw(i) = f0(m-iper+i)
         fw(m+iper+i) = f0(i)
      end do
c     
      return
c     
 20   continue
c     
      xw = x + h*iper
c     
      m0 = ( xw + 0.5*h ) / h
      x0 = m0 * h
c     
      p = (xw-x0) / h
c     
      if ( iop .eq. 1 ) go to 100
c     
      pm1 = p - 1.0
      pm2 = p - 2.0
      pp1 = p + 1.0
      pp2 = p + 2.0
      p2m1 = pm1*pp1
      p2m4 = pm2*pp2
c     
      a(1) =   p2m1*p*pm2 / 24.0
      a(2) = - pm1*p*p2m4 / 6.0
      a(3) =   p2m1*p2m4  / 4.0
      a(4) = - pp1*p*p2m4 / 6.0
      a(5) =   p2m1*p*pp2 / 24.0
c
      f = 0.0
c     
      do i = 1, 5
         f = f + a(i) * fw(m0-2+i)
      end do
c     
      if ( iop .eq. 0 ) return
c     
 100  continue
c     
c...  Derivative.
c     
      p2 = p*p
      p3 = p2*p
c     
      d(1) =   ( 2.0*p3 - 3.0*p2 - p + 1.0 ) / 12.0
      d(2) = - ( 4.0*p3 - 3.0*p2 - 8.0*p + 4.0 ) / 6.0
      d(3) =   ( 2.0*p3 - 5.0*p ) / 2.0
      d(4) = - ( 4.0*p3 + 3.0*p2 - 8.0*p - 4.0 ) / 6.0
      d(5) =   ( 2.0*p3 + 3.0*p2 - p - 1.0 ) / 12.0
c     
      df = 0.0
c     
      do i = 1, 5
         df = df + d(i) * fw(m0-2+i)
      end do
c     
      df = df / h
c     
c$$$  write ( 16,'(i4, 1p6e12.4)')
c$$$  $     m0, x, xw, f0(m0-iper+1),fw(m0+1), f, df
c     
      return
      end
c.........................................................
      subroutine lag ( ax, af, m, nl, x, f, df, iop )
c.........................................................
c     
      dimension ax(*), af(*)
c     
c.....ax is array of independent variable, increasing with index.
c     af is array of dependent variable.
c     m is significant dimension of ax and af on a closed range.
c     nl is desired order for the interpolation.
c     x is the value of x at which f is desired.
c     iop = 0    value of f alone.
c     iop = 1    value 0f df alone.
c     iop = 2    value of f and df.
c     
      zero = 0.0
      one = 1.0
c     
c     calculate index range of the interpolation.
c     
      in = 1
      do 20 i = 1, m
         if ( ax(i) .gt. x ) go to 25        ! ge to gt changed 01-24-03
 20   continue
 25   continue
      in = i - 1
c     
      if ( in .eq. m ) go to 30
      if ( ax(in+1)-x .lt. x-ax(in) ) in = in + 1
c     
 30   continue
c     
c     interpolation will be centered at in except near boundaries of
c     array.  the left and right index are respectively,
c     
      inmm = (nl-0.1) / 2.0
      inpp = (nl+0.1) / 2.0
      nll = in - inmm
      nlr = in + inpp
c     
c.......even nl can be off, so adjust..
c     
      if ( ((nl/2)*2 .eq. nl) .and. (ax(in) .gt. x) ) then
         nll = nll - 1
         nlr = nlr - 1                 ! nlr + 1 to nlr -1 changed 01-24-03
      end if
c     
      if ( (nll.ge.1) .and. (nlr.le.m) ) go to 35
c     
      if ( nlr .le. m ) go to 33
      nlr = m
      nll = nlr - nl + 1
      go to 35
c     
 33   continue
      nll = 1
      nlr = nl
c     
 35   continue
c     
c     now interpolate.
c     
      if ( iop .eq. 1 ) go to 120
c     
      f = zero
c     
      do 100 i = nll, nlr
c     
         alag = one
c     
         do 50 j = nll, nlr
            if ( i .eq. j ) go to 50
            alag = alag * ( x-ax(j) ) / ( ax(i) - ax(j) )
 50      continue
c     
         f = f + alag * af(i)
c     
 100  continue
c     
c     end of f.
c     
      if ( iop .eq. 0 ) return
c     
 120  continue
c     
      df = zero
c     
      do 400 i = nll, nlr
c     
         slag = zero
         do 300 id = nll, nlr
c     
            if ( id .eq. i ) go to 300
c     
            alag = one
            do 200 j = nll, nlr
c     
               if ( j .eq. i ) go to 200
               if ( j .eq. id ) go to 160
c     
               alag = alag*( x-ax(j) ) / ( ax(i)-ax(j) )
               go to 200
c     
 160           continue
               alag = alag / ( ax(i)-ax(id) )
c     
 200        continue
            slag = slag + alag
 300     continue
c     
         df = df + slag * af(i)
c     
 400  continue
c     
      return
      end
c
c...............................................................
c
c
c     ..................................................................
c
c        description of parameters
c           a - original matrix (symmetric), destroyed in computation.
c               resultant eigenvalues are developed in diagonal of
c               matrix a in descending order.
c           r - resultant matrix of eigenvectors (stored columnwise,
c               in same sequence as eigenvalues)
c           n - order of matrices a and r
c           mv- input code
c                   0   compute eigenvalues and eigenvectors
c                   1   compute eigenvalues only (r need not be
c                       dimensioned but must still appear in calling
c                       sequence)
c
c        remarks
c           original matrix a must be real symmetric (storage mode=1)
c           matrix a cannot be in the same location as matrix r
c
c        subroutines and function subprograms required
c           none
c
c        method
c           diagonalization method originated by jacobi and adapted
c           by von neumann for large computers as found in "mathematical
c           methods for digital computers", edited by a. ralston and
c           h.s. wilf, john wiley and sons, new york, 1962, chapter 7
c
c     ..................................................................
c
      subroutine eigen(a,r,n,mv)
c     lcm  (eigen)
      dimension a(*),r(*)
c
c        ...............................................................
c
c        if a double precision version of this routine is desired, the
c        c in column 1 should be removed from the double precision
c        statement which follows.
c
c
c        the c must also be removed from double precision statements
c        appearing in other routines used in conjunction with this
c        routine.
c
c        the double precision version of this subroutine must also
c        contain double precision fortran functions.  sqrt in statements
c        40, 68, 75, and 78 must be changed to dsqrt.  abs in statement
c        62 must be changed to abs. the constant in statement 5 should
c        be changed to 1.0d-12.
c
c        ...............................................................
c
c        generate identity matrix
c
    5 range=1.0e-12
      if(mv-1) 10,25,10
   10 iq=-n
      do 20 j=1,n
      iq=iq+n
      do 20 i=1,n
      ij=iq+i
      r(ij)=0.0
      if(i-j) 20,15,20
   15 r(ij)=1.0
   20 continue
c
c        compute initial and final norms (anorm and anormx)
c
   25 anorm=0.0
      do 35 i=1,n
      do 35 j=i,n
      if(i-j) 30,35,30
   30 ia=i+(j*j-j)/2
      anorm=anorm+a(ia)*a(ia)
   35 continue
      if(anorm) 165,165,40
   40 anorm=1.414*sqrt(anorm)
      anrmx=anorm*range/float(n)
c
c        initialize indicators and compute threshold, thr
c
      ind=0
      thr=anorm
   45 thr=thr/float(n)
   50 l=1
   55 m=l+1
c
c        compute sin and cos
c
   60 mq=(m*m-m)/2
      lq=(l*l-l)/2
      lm=l+mq
   62 if( abs(a(lm))-thr) 130,65,65
   65 ind=1
      ll=l+lq
      mm=m+mq
      x=0.5*(a(ll)-a(mm))
   68 y=-a(lm)/ sqrt(a(lm)*a(lm)+x*x)
      if(x) 70,75,75
   70 y=-y
   75 continue
      yp = 1.0 - y*y
      if ( abs(yp) .lt. 1.0e-10 ) yp = 0.0
      sinx=y/ sqrt(2.0*(1.0+sqrt(yp+0.0)))
      sinx2=sinx*sinx
   78 cosx= sqrt(1.0-sinx2)
      cosx2=cosx*cosx
      sincs =sinx*cosx
c
c        rotate l and m columns
c
      ilq=n*(l-1)
      imq=n*(m-1)
      do 125 i=1,n
      iq=(i*i-i)/2
      if(i-l) 80,115,80
   80 if(i-m) 85,115,90
   85 im=i+mq
      go to 95
   90 im=m+iq
   95 if(i-l) 100,105,105
  100 il=i+lq
      go to 110
  105 il=l+iq
  110 x=a(il)*cosx-a(im)*sinx
      a(im)=a(il)*sinx+a(im)*cosx
      a(il)=x
  115 if(mv-1) 120,125,120
  120 ilr=ilq+i
      imr=imq+i
      x=r(ilr)*cosx-r(imr)*sinx
      r(imr)=r(ilr)*sinx+r(imr)*cosx
      r(ilr)=x
  125 continue
      x=2.0*a(lm)*sincs
      y=a(ll)*cosx2+a(mm)*sinx2-x
      x=a(ll)*sinx2+a(mm)*cosx2+x
      a(lm)=(a(ll)-a(mm))*sincs+a(lm)*(cosx2-sinx2)
      a(ll)=y
      a(mm)=x
c
c        tests for completion
c
c        test for m = last column
c
  130 if(m-n) 135,140,135
  135 m=m+1
      go to 60
c
c        test for l = second from last column
c
  140 if(l-(n-1)) 145,150,145
  145 l=l+1
      go to 55
  150 if(ind-1) 160,155,160
  155 ind=0
      go to 50
c
c        compare threshold with final norm
c
  160 if(thr-anrmx) 165,165,45
c
c        sort eigenvalues and eigenvectors
c
  165 iq=-n
      do 185 i=1,n
      iq=iq+n
      ll=i+(i*i-i)/2
      jq=n*(i-2)
      do 185 j=i,n
      jq=jq+n
      mm=j+(j*j-j)/2
      if(a(ll)-a(mm)) 170,185,185
  170 x=a(ll)
      a(ll)=a(mm)
      a(mm)=x
      if(mv-1) 175,185,175
  175 do 180 k=1,n
      ilr=iq+k
      imr=jq+k
      x=r(ilr)
      r(ilr)=r(imr)
  180 r(imr)=x
  185 continue
      return
      end
c     
c..........................................................
c     
c     
c     ..................................................................
c     
      subroutine mult ( a, b, c, ndim, ndim2, m, l )
c     
      dimension a(ndim,ndim), b(ndim,ndim2), c(ndim,ndim2)
c     
      do 40 j1 = 1, m
         do 40 j2 = 1, l
c     
            sum = 0.0
            do 20 i = 1, m
               sum = sum + a(j1,i) * b(i,j2)
 20         continue
c     
            c(j1,j2) = sum
 40      continue
c     
         return
         end
c     
c.........................................................
c     
c     
      subroutine matmul0 ( a, b, nda,ndb, n, c,ndc )
c     
      dimension a(nda,nda), b(ndb,ndb), c(ndc,ndc)
c     
      do 100 i = 1, n
         do 100 j = 1, n
c     
            sum = 0.0
            do 60 k = 1, n
               sum = sum + a(i,k)*b(k,j)
 60         continue
c     
            c(i,j) = sum
c     
 100     continue
c     
         return
         end
c.........................................................
c     
c     
      subroutine matmul3 ( a, b, nda,ndb, la,lb,lc, c,ndc )
c     
      dimension a(nda,nda), b(ndb,ndb), c(ndc,ndc)
c     
      do 100 i = 1, la
         do 100 j = 1, lc
c     
            sum = 0.0
            do 60 k = 1, lb
               sum = sum + a(i,k)*b(k,j)
 60         continue
c     
            c(i,j) = sum
c     
 100     continue
c     
         return
         end
c     
c...............................................................     
      subroutine dvmat ( a,f, b, nd1,nd2, ni,nj, lrl )
c...............................................................
c     
c.. Direct Multiply matrix B by vector A
c   If lrl = 1, left multiply: A*B
c   If lrl = 2, right multiply: B*A.
c.
      dimension a(*), b(nd1,nd2)
c     
      do i = 1, ni
         do  j = 1, nj
            if ( lrl .eq. 1 ) b(i,j) =  a(i)*b(i,j)
            if ( lrl .eq. 2 ) b(i,j) =  a(j)*b(i,j)
         end do
      end do
c     
         return
         end
c............................................................
c     
c     
      subroutine indef4 ( f, fin, dx, n1, n2, defint, iend )
c     
c     calculates the indefinite and definite integral of f from an
c     equally spaced mesh, dx, from n1 to n2. puts indefinite integral
c     in fin and the definite integral in defint.  if at least one extra
c     data point is available at both ends of the range then set iend=0.
c     m. chance
c     
      dimension f(*), fin(*)
c     
c.....take care of end points.
c     
      dendf = 0.0
      fin(n1) = 0.0
      ind = 1
c     
      if ( iend .eq. 1 ) then
         fin(n1+1) = dx * ( 9.0*f(n1) + 19.0*f(n1+1) - 5.0*f(n1+2)
     $        + f(n1+3) ) / 24.0
         dendf = dx * ( 9.0*f(n2) + 19.0*f(n2-1) - 5.0*f(n2-2)
     $        + f(n2-3) ) / 24.0
         ind = 2
      end if
c     
      do 100 i = n1+ind, n2-ind+1
         ip1 = i+1
         im1 = i-1
         im2 = i-2
         fin(i) = fin(im1) + dx * ( -f(im2) + 13.0*f(im1) + 13.0*f(i)
     $        -f(ip1) ) / 24.0
 100  continue
c     
      fin(n2) = fin(n2-ind+1) + dendf
      defint = fin(n2)
c     
      return
      end
c     
c.................................................................
c     
c     
c..........................................................
      SUBROUTINE difspl ( nsp,the,xin,xout )
c............................................................
c     
      include 'vacuum1.inc'
      dimension  the(*), xin(*), xout(*)
      dimension xpp(nths)
      dimension iop(2), ww1(nths),ww2(nths),ww3(nths),tab(3)
c     
      ns1 = nsp + 1
c     
      iop(1) = 4
      iop(2) = 4
c      dth = the(ns1) / nsp      ! New line since June-3-2003

      call spl1d1(ns1,the,xin,xpp,iop,1,ww1,ww2,ww3)
c     
      do 20 i = 1, ns1
         theta = (i-1) * dth
         call spl1d2 ( ns1,the, xin,xpp, 1, theta, tab )
         xout(i) = tab(2)
 20   continue
c     
      return
      end
c     
c..........................................................
      subroutine difspl2 ( nsp,the,xin,xout,xpp )
c............................................................
c     
      include 'vacuum1.inc'
      dimension  the(*), xin(*), xout(*)
      dimension xpp(nths)
      dimension iop(2), ww1(nths),ww2(nths),ww3(nths),tab(3)
c     
      ns1 = nsp + 1
c     
      iop(1) = 4
      iop(2) = 4
      call spl1d1(ns1,the,xin,xpp,iop,1,ww1,ww2,ww3)
c     
      do 20 i = 1, ns1
         theta = (i-1) * dth
         call spl1d2 ( ns1,the, xin,xpp, 1, theta, tab )
         xout(i) = tab(2)
 20   continue

      xout(nsp+2) = xout(2)
      xpp(nsp+2) = xpp(2)
c     
      return
      end

c..........................................................
      SUBROUTINE difsplt ( nsp,the,xin,xout )
c............................................................
c     
c.. This version calculates zdth internally

      include 'vacuum1.inc'
      dimension  the(*), xin(*), xout(*)
      dimension xpp(nths)
      dimension iop(2), ww1(nths),ww2(nths),ww3(nths),tab(3)
c     
      ns1 = nsp + 1
c     
      iop(1) = 4
      iop(2) = 4

       zdth = the(ns1) / nsp      ! New line since June-3-2003

      call spl1d1(ns1,the,xin,xpp,iop,1,ww1,ww2,ww3)
c     
      do 20 i = 1, ns1
         theta = (i-1) * zdth
         call spl1d2 ( ns1,the, xin,xpp, 1, theta, tab )
         xout(i) = tab(2)
 20   continue
c     
      return
      end
c     
c..............................................................
      subroutine fouran ( gij, gil, cs, m00,l00 )
c..............................................................
c
c... Integrates the second index of GIJ over the JMAX1 functions CS(J,L)
c        and stores results in GIL(M00+I,L00+L)
c    Does not contain the DTH or the 2PI factor.
c
      include 'vacuum1.inc'
c
      real nq
      dimension gij(nths,nths), gil(nths2,nfm2), cs(nths,nfm)
c
c.....the entire cos(lth) matrix is 1-d stored before sin(lth) matrix.
c     
      jmax1 = lmax(1) - lmin(1) + 1
c
      do l1 = 1, jmax1
         do i = 1, mth
            gil(m00+i,l00+l1) = 0.0
         end do
      end do
c
      do 135 l1 = 1, jmax1
c     
         ll = l1 - 1 + lmin(1)
         do 130 j = 1, mth
c     
            do 140 i = 1, mth
               gil(m00+i,l00+l1) = gil(m00+i,l00+l1) +
     $              cs(j,l1)*gij(i,j)
 140        continue
c     
 130     continue
 135  continue
c     
      return
      end
c
c..............................................................
      subroutine foranv ( gil, gll, cs, m00,l00, fff )
c..............................................................
c
c... Integrates the first index of GIL(M00+I,L00+L) over the JMAX1 
c        functions CS(I,L) and stores results in GLL
c....Note that GIL is dimensioned NTHS2,NFM2
c    Contains the DTH and the 2PI factor.
c
c.....The twopi represents the FACTPI needed as explained in VACCAL.
c     fff is a factor for NORM modifications, etc.
c
      include 'vacuum1.inc'
c
      real nq
      dimension gil(nths2,nfm2), gll(nfm,nfm), cs(nths,nfm)
c
c.....the entire cos(lth) matrix is 1-d stored before sin(lth) matrix.
c     
      jmax1 = lmax(1) - lmin(1) + 1
c     
      do l1 = 1, jmax1
         do l2 = 1, jmax1
            gll(l1,l2) = 0.0
         end do
      end do
c
      do l1 = 1, jmax1
         do l2 = 1, jmax1
            do i = 1, mth
               gll(l2,l1) = gll(l2,l1) +
     $              dth*cs(i,l2)*gil(m00+i,l00+l1)*twopi * fff
            end do
         end do
      end do
c     
      return
      end
c
c..............................................................
      subroutine foranv2 ( gil,nid1,nid2,li, gll,nld1,nld2,
     $     cs,ncd1,ncd2,lc, nth, m00,l00, rnge,fff )
c..............................................................
c
c... Integrates the first index of GIL(M00+I,L00+L) over the 
c        LI functions CS(I,L)  and stores results in GLL
c    Contains the DTH factor.
c    fff is a factor for NORM modifications, etc.
c
      include 'vacuum1.inc'
c
c.....The twopi represents the FACTPI needed as explained in VACCAL.
c     fff is a factor for NORM modifications, etc.
c
      real nq
      dimension gil(nid1,nid2), gll(nld1,nld2), cs(ncd1,ncd2)
c
c.....the entire cos(lth) matrix is 1-d stored before sin(lth) matrix.
c     
      zdth = rnge / nth
c     
      write ( 6, '(1x,"li,lc,nth,m00,l00, fff = ",/, 5i4, 1pe13.5)' )
     $     li,lc,nth,m00,l00, fff
c
      do l1 = 1, li
         do l2 = 1, lc
            gll(l1,l2) = 0.0
         end do
      end do
c
      do l1 = 1, li
         do l2 = 1, lc
            do i = 1, nth
               gll(l2,l1) = gll(l2,l1) +
     $              cs(i,l2)*gil(m00+i,l00+l1)*fff*zdth
            end do
         end do
      end do
c     
      return
      end
c
c
c..............................................................
      subroutine foran1 ( gil, gll, cs, m00,l00 )
c..............................................................
c
c... Integrates the first index of GIL(M00+I,L00+L) over the JMAX1 
c        functions CS(I,L) and stores results in GLL
c    Contains the DTH and the 2PI factor.
c....Note that GIL is dimensioned only to NTHS,NFM so M00 and L00 must be 
c      both carefully set to 0.
c    Contains the DTH and the 2PI factor.
c
c.....The twopi represents the FACTPI needed as explained in VACCAL.
c
      include 'vacuum1.inc'
c
      real nq
      dimension gil(nths,nfm), gll(nfm,nfm), cs(nths,nfm)
c
c.....the entire cos(lth) matrix is 1-d stored before sin(lth) matrix.
c     
      jmax1 = lmax(1) - lmin(1) + 1
c     
      do l1 = 1, jmax1
         do l2 = 1, jmax1
            gll(l1,l2) = 0.0
         end do
      end do
c
      do l1 = 1, jmax1
         do l2 = 1, jmax1
            do i = 1, mth
               gll(l2,l1) = gll(l2,l1) +
     $              dth*cs(i,l2)*gil(m00+i,l00+l1)*twopi
            end do
         end do
      end do
c     
      return
      end
c
c..............................................................
      subroutine foura2 ( gij,m01,m02, gil, m00 )
c..............................................................
c
c... Integrates the second index of GIJ(M01+I,M02+J) with COS(ELTHNQ)
c       and SIN(ELTHNQ). Stores the former in GIL(M00+I,L) and the 
c       latter in GIL(M00+I,JMAX1+L)
c    Does not contain the DTH or the 2PI factor.
c
      include 'vacuum1.inc'
c
      real nq
      dimension gij(nths2,nths2), gil(nths,nths)
c
      q = qa1
      nq = n*q
      jmax1 = lmax(1) - lmin(1) + 1
c
c.....the entire cos(lth) matrix is 1-d stored before sin(lth) matrix.
c     
      do l1 = 1, jmax1
         do i = 1, mth
            gil(m00+i,l1) = 0.0
            gil(m00+i,jmax1+l1) = 0.0
         end do
      end do
c
      do 135 l1 = 1, jmax1
c     
         ll = l1 - 1 + lmin(1)
         do 130 j = 1, mth
c     
            theta = (j-1) * dth
            elth = ll * theta
            elthnq = elth + nq*delta(j)
            sinlth = sin(elthnq)
            coslth = cos(elthnq)
c     
            do 140 i = 1, mth
               gil(m00+i,l1) = gil(m00+i,l1) + coslth*gij(i+m01,j+m02)
               gil(m00+i,jmax1+l1) = gil(m00+i,jmax1+l1)
     $              + sinlth*gij(i+m01,j+m02)
 140        continue
c     
 130     continue
 135  continue
c     
      return
      end
c     
c.............................................................
      subroutine fanal ( fth, nt, flc,fls, l1,l2, pi,ddt )
c.............................................................
c
c  Fourier analyse fth(nt) into flc(l) and fls(l). DDT is the shift 
c  the theta coordinate. Used when quantities are shifted in indices,
c  like from GATO's input.
c
      dimension fth(*), flc(*), fls(*)
c
      dt0 = 1.0 / nt
      dt = 2.0 * pi * dt0
      fth(nt+1) = fth(1)
      nl = l2 - l1 + 1
c
      do l = 1, nl
         ll = l1 - 1 + l
         flc(l) = 0.0
         fls(l) = 0.0
         do i = 1, nt
            th = dt*(i-1 + ddt)
            flc(l) = flc(l) + fth(i)*cos(ll*th)
            fls(l) = fls(l) + fth(i)*sin(ll*th)
         end do
         flc(l) = dt0*flc(l)
         fls(l) = dt0*fls(l)
      end do
c
      return
      end
c
c..............................................................
      subroutine fanal0 ( gi,ndi1,ndi2,mi1, goc,gos,ndo1,ndo2 )
c..............................................................
c
      include 'vacuum1.inc'
c
      real nq
      dimension gi(ndi1,ndi2), goc(ndo1,ndo2), gos(ndo1,ndo2)
c     
c......Fourier analyse gi matrix. Cosine and Sine parts.
c  Offset in gi by mi1(input).
c
      q = qa1
      nq = n*q
      jmax1 = lmax(1) - lmin(1) + 1
c
         do  l1 = 1, jmax1
            do l2 = 1, jmax1
               goc(l2,l1) = 0.0
               gos(l2,l1) = 0.0
            end do
         end do
c     
         do l1 = 1, jmax1
            do l2 = 1, jmax1
               ll2 = l2 - 1 + lmin(1)
c     
               do i = 1, mth
                  elth = ll2*(i-1)*dth 
                  goc(l1,l2) = goc(l1,l2)
     $                 + cos(elth) * gi(mi1+i,l1)
                  gos(l1,l2) = gos(l1,l2)
     $                 + sin(elth) * gi(mi1+i,l1)
               end do
c     
               goc(l1,l2) = goc(l1,l2) * dth / pye
               gos(l1,l2) = gos(l1,l2) * dth / pye
c     
            end do
         end do
c
      return
      end
c..............................................................
      subroutine fanal1 ( gi,ndi1,ndi2,mi1, gor,goi,ndo1,ndo2 )
c..............................................................
c
      include 'vacuum1.inc'
c
      real nq
      dimension gi(ndi1,ndi2), gor(ndo1,ndo2), goi(ndo1,ndo2)
c     
c......fourier analyse gi matrix. Real and imaginary parts.
c  Offset in gi by gi(mi1(input),jmax1)
c  Note that indices of input is (obs,srce), output is (srce,obs)
c
      q = qa1
      nq = n*q
      jmax1 = lmax(1) - lmin(1) + 1
c
c.....the entire cos(lth) matrix is 1-d stored before sin(lth) matrix.
c     
         do  l1 = 1, jmax1
            do l2 = 1, jmax1
               gor(l2,l1) = 0.0
               goi(l2,l1) = 0.0
            end do
         end do
c     
         do l1 = 1, jmax1
            do l2 = 1, jmax1
               ll2 = l2 - 1 + lmin(1)
c     
               do i = 1, mth
                  elth = ll2*(i-1)*dth 
                  gor(l1,l2) = gor(l1,l2)
     $                 + cos(elth) * gi(mi1+i,l1)
     $                 + sin(elth) * gi(mi1+i,jmax1+l1)
                  goi(l1,l2) = goi(l1,l2)
     $                 + cos(elth) * gi(mi1+i,jmax1+l1)
     $                 - sin(elth) * gi(mi1+i,l1)
               end do
c     
               gor(l1,l2) = gor(l1,l2) * dth / twopi
               goi(l1,l2) = goi(l1,l2) * dth / twopi
c     
            end do
         end do
c
      return
      end
c     
c..............................................................
      SUBROUTINE fanal1d ( gi,ndi1,ndi2,mi1, gor,goi,ndo1,ndo2 )
c..............................................................
c
      INCLUDE 'vacuum1.inc'
      INCLUDE 'vacuum5.inc'
c
      REAL nq
      DIMENSION gi(ndi1,ndi2), gor(ndo1,ndo2), goi(ndo1,ndo2)
c     
c......fourier analyse gi matrix. Real and imaginary parts.
c      This has the coordinate function, DELTA, included.
c  Offset in gi by gi(mi1(input),jmax1)
c  Note that indices of input is (obs,srce), output is (srce,obs)
c
      q = qa1
      nq = n*q
      jmax1 = lmax(1) - lmin(1) + 1
c
c.....the entire cos(lth) matrix is 1-d stored before sin(lth) matrix.
c     
         DO  l1 = 1, jmax1
            DO l2 = 1, jmax1
               gor(l2,l1) = 0.0
               goi(l2,l1) = 0.0
            END DO
         END DO
C     
         DO l1 = 1, jmax1
            DO l2 = 1, jmax1
               ll2 = l2 - 1 + lmin(1)
c     
               DO i = 1, mth
                  elth = ll2*(i-1)*dth 
                  gor(l1,l2) = gor(l1,l2)
     $                 + cslth(i,l2) * gi(mi1+i,l1)
     $                 + snlth(i,l2) * gi(mi1+i,jmax1+l1)
                  goi(l1,l2) = goi(l1,l2)
     $                 + cslth(i,l2) * gi(mi1+i,jmax1+l1)
     $                 - snlth(i,l2) * gi(mi1+i,l1)
               END DO
c     
               gor(l1,l2) = gor(l1,l2) * dth / twopi
               goi(l1,l2) = goi(l1,l2) * dth / twopi
c     
            END DO
         END DO
C
      RETURN
      END
c     

c..............................................................
      subroutine fanal2 ( gi1,gi2,ndi1,ndi2,mi1, gor,goi,ndo1,ndo2 )
c..............................................................
c
      include 'vacuum1.inc'
c
      real nq
      dimension gi1(ndi1,ndi2), gi2(ndi1,ndi2),
     $     gor(ndo1,ndo2), goi(ndo1,ndo2)
c     
c......Fourier analyse gi matrices: gi1, gi2 are cosine and sine parts..
c  Real and imaginary parts in go1, go2.
c
      q = qa1
      nq = n*q
      jmax1 = lmax(1) - lmin(1) + 1
c
         do  l1 = 1, jmax1
            do l2 = 1, jmax1
               gor(l2,l1) = 0.0
               goi(l2,l1) = 0.0
            end do
         end do
c     
         do l1 = 1, jmax1
            do l2 = 1, jmax1
               ll2 = l2 - 1 + lmin(1)
c     
               do i = 1, mth
                  elth = ll2*(i-1)*dth 
                  gor(l1,l2) = gor(l1,l2)
     $                 + cos(elth) * gi1(mi1+i,l1)
     $                 + sin(elth) * gi2(mi1+i,l1)
                  goi(l1,l2) = goi(l1,l2)
     $                 + cos(elth) * gi2(mi1+i,l1)
     $                 - sin(elth) * gi1(mi1+i,l1)
               end do
c     
               gor(l1,l2) = gor(l1,l2) * dth / twopi
               goi(l1,l2) = goi(l1,l2) * dth / twopi
c     
            end do
         end do
c
      return
      end
c     
c*****************************************************************
c................................................................
      subroutine felang ( gij, gil, cs, m00,l00 )
c................................................................
c 
      data izcal / 0 /
      include 'vacuum1.inc'
c     
      real nq
      dimension gij(nths,nths), gil(nths2,nfm2), cs(*)
c     
      izcal = izcal + 1
      nzwrt = 8
      nzwdel = mth / nzwrt
c
c Normalize finite elements to Sqrt(MFEL). Put in integration weights
c for convenience.
      znorm = sqrt( float(mfel) )
      zwt1 = 0.5 * znorm
      zwt2 = 1.0 * znorm
      zws1 = 1.0 * znorm /3.0
      zws2 = 2.0 * zws1
      zws4 = 2.0 * zws2
c     
      nzdel = ndfel
      nzdel1 = nzdel + 1
c     
      q = qa1
      nq = n*q
c
      do l1 = 1, mfel
         do i = 1, mth
            gil(m00+i,l00+l1) = 0.0
         end do
      end do
c     
c..   Loop over observer points.
            do i = 1, mth
c     
c..   Loop over finite elements.
      do l1 = 1, mfel
               izwrt = 0
               if ( (i .eq. (i/nzwdel * nzwdel + 1))
     $        .and. ( l1 .eq. (l1/8 * 8 +1 ))
     $        .and. izcal .eq. 1 ) then
                  izwrt = 1
c$$$                  write ( outmod, '(".......................")' )
               end if
c     
         mzl = (l1-1) * nzdel + 1
         mzr = mzl + nzdel
c     
c..   Loop over source integral.
         do j = 1, nzdel1
c     
            jth0 = mzl + j - 2
            jth = mod(jth0,mth)
            if ( jth .lt. 0 ) jth = mth + jth
            jth1 = jth + 1
c
               if ( nzdel .eq. 1 ) then
                  gil(m00+i,l00+l1) = gil(m00+i,l00+l1) +
     $                 zwt1 * gij(i,jth1) * cs(jth1)
               else
                  zwt = zws2
                  if ( (j/2)*2 .eq. j ) zwt = zws4
                  if ( (j .eq. 1) .or. (j .eq. nzdel1) ) zwt = zws1
                  gil(m00+i,l00+l1) = gil(m00+i,l00+l1) +
     $                 zwt * gij(i,jth1) * cs(jth1)
               end if
c$$$               if ( izwrt .eq. 1 )
c$$$     $              write ( outmod, '("Iobs, Mel, Jsce, JTH1, Zwt = ",
c$$$     $              4i4, f9.3)' ) i, l1, j, jth1, zwt
            end do
         end do
      end do
c     
      return
      end
c................................................................
      subroutine felanv ( gil, gll, cs, m00,l00 )
c................................................................
c     
c.....The twopi represents the FACTPI needed as explained in VACCAL.c     
c
      include 'vacuum1.inc'
c     
      data izcal / 0 /
      real nq
      dimension gil(nths2,nfm2), gll(nfm,nfm), cs(*)
c     
c Normalize finite elements to Sqrt(MFEL). Put in integration weights
c for convenience.
c     
      izcal = izcal + 1
c
      znorm = sqrt( float(mfel) )
      zwt1 = 0.5 * dth * znorm
      zwt2 = dth * znorm
      zws1 = dth * znorm /3.0
      zws2 = 2.0 * zws1 
      zws4 = 2.0 * zws2 
c     
      nzdel = ndfel
      nzdel1 = nzdel + 1
c     
      q = qa1
      nq = n*q
c     
      do l1 = 1, mfel
         do l2 = 1, mfel
            gll(l1,l2) = 0.0
         end do
      end do
c     
c..   Loop over original source finite elements.
            do l2 = 1, mfel
c
c..   Loop to integrate over observer finite elements.
      do l1 = 1, mfel
               izwrt = 0
               if ((l2 .eq. (l2/8 * 8 + 1))
     $              .and. (l1 .eq. (l1/8 * 8 +1))
     $              .and. (izcal .eq. 1)) then
                  izwrt = 1
c                  write ( outmod, '("..........................")' )
               end if
c     
         mzl = (l1-1) * nzdel + 1
         mzr = mzl + nzdel
c     
c..   Loop over obvserver integral.
         do j = 1, nzdel1
c     
            jth0 = mzl + j - 2
            jth = mod(jth0,mth)
            if ( jth .lt. 0 ) jth = mth + jth
            jth1 = jth + 1
               if ( nzdel .eq. 1 ) then
                  gll(l1,l2) = gll(l1,l2) +
     $                 zwt1 * gil(m00+jth1,l00+l2) * cs(jth1) * twopi
               else
                  zwt = zws2
                  if ( (j/2)*2 .eq. j ) zwt = zws4
                  if ( (j .eq. 1) .or. (j .eq. nzdel1) ) zwt = zws1
                  if ( jth1 .ne. mth1 ) then
                     gll(l1,l2) = gll(l1,l2) +
     $                    zwt * gil(m00+jth1,l00+l2) * cs(jth1) * twopi
                  else 
                     gll(l1,l2) = gll(l1,l2)
     $                    + zwt * gil(m00+1,l00+l2) * cs(jth1) * twopi
                  end if
               end if
c$$$               if ( izwrt .eq. 1 ) 
c$$$     $              write ( outmod, '("L2, L1, Jobs, Jth1, zwt = ",
c$$$     $              4i4, f9.3 )' ) l2, l1, j, jth1, zwt
            end do
         end do
      end do
c     
      return
      end
c     
c................................................................
      subroutine felang3 ( gij,nd1,nd2, gil, m00,l00 )
c................................................................
c     
c.... Finite element alaready built in so just take gij into gil
c     with znorm
c... !!!!! ZNORM SHOULD PROBABLY BE 1/SQRT(mfel) IF SOURCE IS 
c      MEANT TO BE B(theta_k) RATHER THAN  B_k !!!!!!
C     OR ZNORM = 1 IF SOURCE IS MEANT TO BE B_k. THIS MAY ACCOUNT FOR 
C     THE FACTOR OF mfel IN THE VACUUM MATRIX.
c     
      include 'vacuum1.inc'
c     
      real nq
      dimension gij(nd1,nd2), gil(nths2,nfm2)
c     
c     Normalize finite elements to Sqrt(MFEL). Put in integration weights
c     for convenience.
      znorm = sqrt( float(mfel) )
c     
c..   Loop over observer points.
      do i = 1, mth
c     
c..   Loop over finite elements.
         do l1 = 1, mfel
c     
            gil(m00+i,l00+l1) =  znorm * gij(i,l1) 
c     
         end do
      end do
c     
      return
      end
c................................................................
      subroutine felanv3 ( gil, gll, cs, m00,l00 )
c................................................................
c     
c.....The twopi represents the FACTPI needed as explained in VACCAL.c     
c     
      include 'vacuum1.inc'
c     
      data izcal / 0 /
      real nq
      dimension gil(nths2,nfm2), gll(nfm,nfm), cs(*)
c     
c     Normalize finite elements to Sqrt(MFEL). Put in integration weights
c     for convenience.
c     
      izcal = izcal + 1
      tpdth = twopi * dth
c     
      znorm = sqrt( float(mfel) )
c     
      nzdel = ndfel
      nzdel1 = nzdel + 1
c     
      q = qa1
      nq = n*q
c     
      do l1 = 1, mfel
         do l2 = 1, mfel
            gll(l1,l2) = 0.0
         end do
      end do
c     
c..   Loop over original source finite elements.
      do l2 = 1, mfel
c     
c..   Loop to integrate over observer finite elements.
         do l1 = 1, mfel
c
            gll(l1,l2) = gll(l1,l2) +
     $           znorm * gil(m00+l1,l00+l2) * cs(l1) * tpdth
c
         end do
      end do
c     
      return
      end


