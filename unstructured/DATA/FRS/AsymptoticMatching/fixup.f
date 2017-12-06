      subroutine fixup
      entry reginit
c     entry walld
      entry rawmeas
      entry regler
      entry vforce
      entry restcon
      entry vforcepl
      entry fpplot
c     entry hyper
c     entry lingin1
c     entry lingin2
c     entry lingin3
      entry appvolto
      entry spdpbx
      entry spdtftr
c     entry spdd3d
      entry spdpbxm
      entry spdpbxn
      entry fedtsc
      entry spdasdex
      entry missionc
      entry missionw
      entry getlscrestart
      entry putlscrestart
      entry growth
c     entry lingplt
c     entry lingpl2
c     entry lingpl3
      entry colorc
c     entry balloon
      entry lsc
c     entry tridiag
      entry d02bae
c     for double precision nag use d02baf
      entry clock
      return
      end
      subroutine second(tcpu)
      real*4 etime
      external etime
      dimension tarray(2)
      tcpu = etime(tarray)
      return
      end
      subroutine tremain(secleft)
c
c.....returns the number of seconds left before
c.....program times out
      secleft = 1000.
      return
      end
      subroutine s15aee(a1,n1)
      return
      end
      subroutine f04aae(a1,n1,a2,n2,n3,n4,a3,n5,a4,n6)
      call       f04aaf(a1,n1,a2,n2,n3,n4,a3,n5,a4,n6)
      return
      end
ckuma st at end
      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.

        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY

        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue

      return
      END
