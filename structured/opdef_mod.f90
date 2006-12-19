module opdef_mod

implicit none
contains

!==============================================================
subroutine ireaddef
  use basic

  integer :: itype

  if(myrank.eq.0 .and. iprint.gt.0) then
     print *, "Entering ireaddef."
  endif


  iread = 0
  ntensor = 79

  do itype=1,ntensor
!
! possibly skip all gyroviscous terms:
     if(gyro.ne.1) then
        if((itype.ge.51).and.(itype.le.72)) go to 400
     endif
! possibly skip all anisotropic head conduction
     if(chipar.eq.0.) then
        if((itype.ge.73).and.(itype.le.79)) go to 400
     endif
!
! the following are not used
     if(itype.eq.34) go to 400
     if(itype.eq.35) go to 400
     if(itype.eq.36) go to 400
     if(itype.eq.49) go to 400
     if(itype.eq.50) go to 400
     if(numvar.le.1) then
        if(itype.eq. 1) go to 400
        if(itype.eq. 2) go to 400 
        if(itype.eq. 3) go to 400
        if(itype.eq. 4) go to 400
!!$        if(itype.eq. 5) go to 400
        if(itype.eq. 6) go to 400
        if(itype.eq. 9) go to 400
        if(itype.eq.10) go to 400
!
        if(itype.eq.17) go to 400
        if(itype.eq.18) go to 400
        if(itype.eq.21) go to 400
        if(itype.eq.22) go to 400
!
        if(itype.eq.47) go to 400
        if(itype.eq.48) go to 400
        if(itype.eq.51) go to 400
        if(itype.eq.52) go to 400
        if(itype.eq.53) go to 400
        if(itype.eq.55) go to 400
        if(itype.eq.56) go to 400
        if(itype.eq.57) go to 400
        if(itype.eq.61) go to 400
        if(itype.eq.66) go to 400
        if(itype.eq.67) go to 400
        if(itype.eq.68) go to 400
     endif
     if(numvar.le.2) then
        if(itype.eq. 8) go to 400
!
        if(itype.eq.15) go to 400
        if(itype.eq.16) go to 400
        if(itype.eq.19) go to 400
        if(itype.eq.20) go to 400
        if(itype.eq.23) go to 400
        if(itype.eq.24) go to 400
        if(itype.eq.25) go to 400
        if(itype.eq.26) go to 400
        if(itype.eq.27) go to 400
        if(itype.eq.28) go to 400
        if(itype.eq.29) go to 400
        if(itype.eq.30) go to 400
        if(itype.eq.31) go to 400
        if(itype.eq.32) go to 400
!
        if(itype.eq.35) go to 400
        if(itype.eq.37) go to 400
        if(itype.eq.38) go to 400
        if(itype.eq.39) go to 400
        if(itype.eq.40) go to 400
        if(itype.eq.41) go to 400
        if(itype.eq.42) go to 400
        if(itype.eq.43) go to 400
        if(itype.eq.44) go to 400
        if(itype.eq.45) go to 400
        if(itype.eq.46) go to 400
        if(itype.eq.58) go to 400
        if(itype.eq.59) go to 400
        if(itype.eq.60) go to 400
        if(itype.eq.62) go to 400
        if(itype.eq.63) go to 400
        if(itype.eq.64) go to 400
        if(itype.eq.65) go to 400
        if(itype.eq.69) go to 400
        if(itype.eq.70) go to 400
        if(itype.eq.71) go to 400
        if(itype.eq.72) go to 400
        if(itype.eq.73) go to 400
        if(itype.eq.74) go to 400
        if(itype.eq.75) go to 400
        if(itype.eq.76) go to 400
        if(itype.eq.77) go to 400
        if(itype.eq.78) go to 400
        if(itype.eq.79) go to 400
     endif
     !
     iread(itype) = 1

400  continue
  end do

  if(myrank.eq.0 .and. iprint.gt.0) then
     print *, "Exiting ireaddef."
  endif


end subroutine ireaddef

!==========================================================
subroutine opdef1
  use p_data
  use t_data
  use basic
  use arrays
  
  implicit none
    
  integer :: iodd, i, j, k, l, ll,  msum, nsum, itri
  real :: rinvjps
  integer :: ii, iii, jj, jjj, jp, jjp, jjjp
  real :: terma, termb, termd, termd2, termq, termx, termy
  real :: gfac, gk, gl, fackj, sum, co, sn
  real :: sumg0, sumg2, sumg3, sumg4, sumg5, sumg6, sumg7, sumg8, sumg9
  real :: sug10, sug11, sug12, sug13, sug14, sug15, sug16            
  real :: sumh3, sumh5, sumk0, sumk1, sumk2, sumh2
  real :: sumx0, sumy0, sumx1, sumy1, sumx2, sumy2
  real :: f

  if(myrank.eq.0 .and. iprint.gt.0) then
     print *, "Entering opdef1."
  endif
  
  do iodd=1,ioddmx
  
  co = cos(ttri(iodd))
  sn = sin(ttri(iodd))

  !  calculate matrix elements of the integration matrix fint
  do i=-6,maxi
     do j=-6,maxi
        fint(i,j) = 0
        if(i.ge.0 .and.j.ge.0) &           
             fint(i,j) = f(i,j,atri(iodd),btri(iodd),ctri(iodd))
     enddo
  enddo
     
  ! begin local integrations of metric elements
  do p=1,20
     d2matrix(p) = fint(mi(p),ni(p))
     do q=1,20
        msum = mi(p)+mi(q)
        nsum = ni(p)+ni(q)
        
        ! mass matrix
        dmatrix(p,q) = fint(msum,nsum)
        
        ! del-sq operator
        amatrix(p,q) =  -(mi(p)*mi(q)*fint(msum-2,nsum) + & 
             ni(p)*ni(q)*fint(msum,nsum-2))
        
        ! del-star operator
        qmatrix(p,q) =  (mi(p)+mi(q)-1)*mi(p)*fint(msum-2,nsum) &
             + (ni(p)+ni(q)-1)*ni(p)*fint(msum,nsum-2)
        
        ! del-fourth operator
        bmatrix(p,q) =(mi(p)*(mi(p)-1)*mi(q)*(mi(q)-1)*fint(msum-4,nsum)    &
                     + ni(p)*(ni(p)-1)*ni(q)*(ni(q)-1)*fint(msum,nsum-4)    &
                    + (mi(p)*(mi(p)-1)*ni(q)*(ni(q)-1)                      &
                    +  mi(q)*(mi(q)-1)*ni(p)*(ni(p)-1))*fint(msum-2,nsum-2))
           
        ! x derivative
        xmatrix(p,q) = co*mi(q)*fint(msum-1,nsum) &
                      -sn*ni(q)*fint(msum,nsum-1)
        ! y derivative
        ymatrix(p,q) = sn*mi(q)*fint(msum-1,nsum) &
                      +co*ni(q)*fint(msum,nsum-1)
           
        do r=1,20
           msum = mi(q)+mi(r)+mi(p)
           nsum = ni(q)+ni(r)+ni(p)
              
           ! poisson bracket involving del-sq (note sign)
           termg0(p,q,r) =  -(mi(p)*ni(r)-mi(r)*ni(p))             &
                *(mi(q)*(mi(q)-1)*fint(msum-3,nsum-1)              &
                + ni(q)*(ni(q)-1)*fint(msum-1,nsum-3))
              
           termg2(p,q,r) =  -(mi(p)*mi(q)*fint(msum-2,nsum)        &
                             +ni(p)*ni(q)*fint(msum,nsum-2))
              
           termg3(p,q,r) =                                                   &
                mi(p)*mi(q)*(mi(q)-1)*(mi(p)+mi(r)-1) *fint(msum-4,nsum)     &
                +( ni(p)*mi(q)*(mi(q)-1)*(ni(p)+ni(r)-1)                     &
                  +mi(p)*ni(q)*(ni(q)-1)*(mi(p)+mi(r)-1))*fint(msum-2,nsum-2)&
                  +ni(p)*ni(q)*(ni(q)-1)*(ni(p)+ni(r)-1)*fint(msum,nsum-4)

           termg4(p,q,r)=mi(q)*(mi(q)-1)*mi(p)*mi(r) *fint(msum-4,nsum)  &
                      + (mi(q)*(mi(q)-1)*ni(p)*ni(r)                     & 
                        +ni(q)*(ni(q)-1)*mi(p)*mi(r))*fint(msum-2,nsum-2)&
                        +ni(q)*(ni(q)-1)*ni(p)*ni(r) *fint(msum,nsum-4)

           termg5(p,q,r) = mi(p)*mi(r)*fint(msum-2,nsum)                &
                          +ni(p)*ni(r)*fint(msum,nsum-2)

           termg6(p,q,r) =                                               &
                - mi(p)*mi(r)*mi(q)*(mi(r)+mi(q)-2) *fint(msum-4,nsum)   &
                -(mi(p)*ni(r)*ni(q)*(mi(r)+mi(q)  )                      &
                 +ni(p)*mi(r)*mi(q)*(ni(r)+ni(q)  ))*fint(msum-2,nsum-2) &
                - ni(p)*ni(r)*ni(q)*(ni(r)+ni(q)-2) *fint(msum  ,nsum-4) &
             +gam*(mi(q)*(mi(q)-1)*mi(p)*(mi(p)-1) *fint(msum-4,nsum)    &
                +(mi(q)*(mi(q)-1)*ni(p)*(ni(p)-1)                        &
                 +mi(p)*(mi(p)-1)*ni(q)*(ni(q)-1))*fint(msum-2,nsum-2)   &
                 +ni(q)*(ni(q)-1)*ni(p)*(ni(p)-1) *fint(msum  ,nsum-4))

           termg7(p,q,r) =  -(mi(r)*ni(q)-mi(q)*ni(r))               &
                *(mi(p)*(mi(r)+mi(q)-1)*fint(msum-3,nsum-1)          &
                 +ni(p)*(ni(r)+ni(q)-1)*fint(msum-1,nsum-3))

           termg8(p,q,r) =                                          &
                -(mi(q)*(mi(q)-1)*mi(p)*mi(r) *fint(msum-4,nsum)    &
               + (mi(q)*(mi(q)-1)*ni(p)*ni(r)                       &
                 +ni(q)*(ni(q)-1)*mi(p)*mi(r))*fint(msum-2,nsum-2)  &
               +  ni(q)*(ni(q)-1)*ni(p)*ni(r) *fint(msum  ,nsum-4))

           termg9(p,q,r) =                                           &
                (-(mi(p)*ni(r)-mi(r)*ni(p))*mi(q)*(mi(q)-1)          &
                -2*mi(q)*ni(q)*mi(r)*(mi(r)-1)                       &
                +2*mi(r)*ni(r)*mi(q)*(mi(q)-1))*fint(msum-3,nsum-1)  &
               +(-(mi(p)*ni(r)-mi(r)*ni(p))*ni(q)*(ni(q)-1)          &
                +2*mi(q)*ni(q)*ni(r)*(ni(r)-1)                       &
                -2*mi(r)*ni(r)*ni(q)*(ni(q)-1))*fint(msum-1,nsum-3)    

           terg10(p,q,r) =    0.5*(                                &
                mi(p)*(mi(p)-1)*mi(r)*mi(q)*fint(msum-4,nsum)      &
                +(mi(p)*(mi(p)-1)*ni(r)*ni(q)                      &
                + ni(p)*(ni(p)-1)*mi(r)*mi(q))*fint(msum-2,nsum-2) &
                + ni(p)*(ni(p)-1)*ni(r)*ni(q)*fint(msum,  nsum-4))

           terg11(p,q,r) =  2.0*(                                     &
                mi(q)*ni(q)*mi(r)*ni(r)                               &
                -mi(q)*(mi(q)-1)*ni(r)*(ni(r)-1))*fint(msum-2,nsum-2)
!
!...........changed on 12/15/06...scj
            terg12(p,q,r) =                                          &
                  mi(q)*(mi(q)-1)*mi(r)*(mi(r)-1)*fint(msum-4,nsum)  &
                 +ni(q)*(ni(q)-1)*ni(r)*(ni(r)-1)*fint(msum,nsum-4)  &
                 +2*mi(q)*ni(q)*mi(r)*ni(r)*fint(msum-2,nsum-2)
               
!           terg12(p,q,r) =                                          &
!                  mi(q)*(mi(q)-1)*(mi(p)+mi(r))*(mi(p)+mi(r)-1)     &
!                *fint(msum-4,nsum)                                  &
!                +(mi(q)*(mi(q)-1)*(ni(p)+ni(r))*(ni(p)+ni(r)-1)     &
!                + ni(q)*(ni(q)-1)*(mi(p)+mi(r))*(mi(p)+mi(r)-1))    &
!                *fint(msum-2,nsum-2)                                &
!                + ni(q)*(ni(q)-1)*(ni(p)+ni(r))*(ni(p)+ni(r)-1)     &
!                *fint(msum,nsum-4)

           terg13(p,q,r) =                                        &
                 ((mi(q)+mi(r))*(mi(q)+mi(r)-1)*mi(p)*(mi(p)-1))  &
                *fint(msum-4,nsum)                                &
               +(((mi(q)+mi(r))*(mi(q)+mi(r)-1)*ni(p)*(ni(p)-1))  &
                +((ni(q)+ni(r))*(ni(q)+ni(r)-1)*mi(p)*(mi(p)-1))) &
                *fint(msum-2,nsum-2)                              &
                +((ni(q)+ni(r))*(ni(q)+ni(r)-1)*ni(p)*(ni(p)-1))  &
                *fint(msum,nsum-4)

           terg14(p,q,r) =                                        &
                   mi(q)*(mi(q)-1)*mi(r)*(mi(r)-1)                &
                *fint(msum-4,nsum)                                &
                +(-ni(q)*(ni(q)-1)*mi(r)*(mi(r)-1)                &
                -  mi(q)*(mi(q)-1)*ni(r)*(ni(r)-1)                &
                + 4*mi(q)*ni(q)*mi(r)*ni(r))*fint(msum-2,nsum-2)  &
                +  ni(q)*(ni(q)-1)*ni(r)*(ni(r)-1)                &
                *fint(msum,nsum-4)

           terg15(p,q,r) =                                            &
                 2.*mi(q)*(mi(q)-1)*mi(r)*(mi(r)-1)*fint(msum-4,nsum) &
                +4.*mi(q)*ni(q)*mi(r)*ni(r)*fint(msum-2,nsum-2)       &
                +2.*ni(q)*(ni(q)-1)*ni(r)*(ni(r)-1)*fint(msum,nsum-4)

           terg16(p,q,r) =                                            &
                 4.*mi(q)*mi(r)*(ni(q)*(mi(r)-1)-ni(r)*(mi(q)-1))     &
                *fint(msum-3,nsum-1)                                  &
                -4.*ni(q)*ni(r)*(mi(q)*(ni(r)-1)-mi(r)*(ni(q)-1))     &
                *fint(msum-1,nsum-3)

           termh3(p,q,r) =                                            &
                -mi(p)*mi(q)*mi(r)*(mi(r)+mi(q)-2)*fint(msum-4,nsum)  &
               -(mi(q)*ni(r)*(mi(p)*(ni(r)-1)+ni(p)*(mi(q)-1))        &
                +ni(q)*mi(r)*(ni(p)*(mi(r)-1)+mi(p)*(ni(q)-1)))       &
                *fint(msum-2,nsum-2)                                  &
                -ni(p)*ni(q)*ni(r)*(ni(r)+ni(q)-2)*fint(msum,nsum-4)

           termh5(p,q,r) = -mi(p)*(mi(r)+mi(q))*fint(msum-2,nsum)   &
                           -ni(p)*(ni(r)+ni(q))*fint(msum,nsum-2)

           ! Grad-Shafranov operator
           termh(p,q,r) =-(mi(p)*mi(q)*fint(msum-2,nsum)         &
                         + ni(p)*ni(q)*fint(msum  ,nsum-2))

           ! poisson bracket of 2 scalars
           termk0(p,q,r) =(mi(q)*ni(r)-mi(r)*ni(q))*fint(msum-1,nsum-1)

           termk1(p,q,r) =  fint(msum,nsum)

           termk2(p,q,r) = -(mi(q)*(mi(q)-1)*fint(msum-2,nsum)   &
                             +ni(q)*(ni(q)-1)*fint(msum,nsum-2))

           ! x, y derivative
           termx0(p,q,r) = co*mi(q)*fint(msum-1,nsum)             &
                          -sn*ni(q)*fint(msum,nsum-1)
           termy0(p,q,r) = co*ni(q)*fint(msum,nsum-1)             &
                          +sn*mi(q)*fint(msum-1,nsum)

           ! x, y derivative of poisson bracket
           termx1(p,q,r) = (mi(q)*ni(r) - mi(r)*ni(q))        &
                *(co*mi(p)*fint(msum-2,nsum-1)   &
                 -sn*ni(p)*fint(msum-1,nsum-2))
           termy1(p,q,r) = (mi(q)*ni(r) - mi(r)*ni(q))        &
                *(sn*mi(p)*fint(msum-2,nsum-1)   &
                 +co*ni(p)*fint(msum-1,nsum-2))
           termx2(p,q,r) = -mi(q)*(mi(r)+mi(q)-1)             &
                *(co*mi(p)*fint(msum-3,nsum)     &
                 -sn*ni(p)*fint(msum-2,nsum-1))  &
                           -ni(q)*(ni(r)+ni(q)-1)             &
                *(co*mi(p)*fint(msum-1,nsum-2)   &
                 -sn*ni(p)*fint(msum,nsum-3))
           termy2(p,q,r) = -mi(q)*(mi(r)+mi(q)-1)             &
                *(sn*mi(p)*fint(msum-3,nsum)     &
                 +co*ni(p)*fint(msum-2,nsum-1))  &
                           -ni(q)*(ni(r)+ni(q)-1)             &
                *(sn*mi(p)*fint(msum-1,nsum-2)   &
                 +co*ni(p)*fint(msum,nsum-3))
        enddo
     enddo
  enddo
  
  ! Transform from local coordinates to global coordinates
  do i=1,18
     termd2 = 0.
     do k=1,20
        termd2 = termd2 + d2matrix(k)*gtri(k,i,iodd)
     enddo
     d2term(iodd,i) = termd2
     
     do j=1,18
        
        termd = 0
        termb = 0
        terma = 0
        termq = 0
        termx = 0
        termy = 0
        do k=1,20
           
           do l=1,20
              fackj =  gtri(k,i,iodd)*gtri(l,j,iodd)
              termd =  termd + dmatrix(k,l)*fackj
              termb =  termb + bmatrix(k,l)*fackj
              terma =  terma + amatrix(k,l)*fackj
              termq =  termq + qmatrix(k,l)*fackj
              termx =  termx + xmatrix(k,l)*fackj
              termy =  termy + ymatrix(k,l)*fackj
           enddo
        enddo
        aterm(iodd,i,j) = terma
        dterm(iodd,i,j) = termd
        bterm(iodd,i,j) = termb
        sterm(iodd,i,j) = termq
        xterm(iodd,i,j) = termx
        yterm(iodd,i,j) = termy
        
        do jp=1,18
           sumg0 = 0.
           sumg2 = 0.
           sumg3 = 0.
           sumg4 = 0.
           sumg5 = 0.
           sumg6 = 0.
           sumg7 = 0.
           sumg8 = 0.
           sumg9 = 0.
           sug10 = 0.
           sug11 = 0.
           sug12 = 0.
           sug13 = 0.
           sug14 = 0.
           sug15 = 0.
           sug16 = 0.
           sumh3 = 0.
           sumh5 = 0.
           sumk0 = 0.
           sumk1 = 0.
           sumk2 = 0.
           sumx0 = 0.
           sumy0 = 0.
           sumx1 = 0.
           sumx2 = 0.
           sumy1 = 0.
           sumy2 = 0.
           do p = 1,20
              gk = gtri(p,jp,iodd)
              do k = 1,20
                 gl = gtri(k,i,iodd)
                 do l = 1,20 
                    gfac = gk*gl*gtri(l,j,iodd)
                    sumg0 = sumg0 + termg0(k,l,p)*gfac
                    sumg2 = sumg2 + termg2(k,l,p)*gfac
                    sumg3 = sumg3 + termg3(k,l,p)*gfac
                    sumg4 = sumg4 + termg4(k,l,p)*gfac
                    sumg5 = sumg5 + termg5(k,l,p)*gfac
                    sumg6 = sumg6 + termg6(k,l,p)*gfac
                    sumg7 = sumg7 + termg7(k,l,p)*gfac
                    sumg8 = sumg8 + termg8(k,l,p)*gfac
                    sumg9 = sumg9 + termg9(k,l,p)*gfac
                    sug10 = sug10 + terg10(k,l,p)*gfac
                    sug11 = sug11 + terg11(k,l,p)*gfac
                    sug12 = sug12 + terg12(k,l,p)*gfac
                    sug13 = sug13 + terg13(k,l,p)*gfac
                    sug14 = sug14 + terg14(k,l,p)*gfac
                    sug15 = sug15 + terg15(k,l,p)*gfac
                    sug16 = sug16 + terg16(k,l,p)*gfac
                    sumh3 = sumh3 + termh3(k,l,p)*gfac
                    sumh5 = sumh5 + termh5(k,l,p)*gfac
                    sumk0 = sumk0 + termk0(k,l,p)*gfac
                    sumk1 = sumk1 + termk1(k,l,p)*gfac
                    sumk2 = sumk2 + termk2(k,l,p)*gfac
                    sumx0 = sumx0 + termx0(k,l,p)*gfac
                    sumy0 = sumy0 + termy0(k,l,p)*gfac
                    sumx1 = sumx1 + termx1(k,l,p)*gfac
                    sumy1 = sumy1 + termy1(k,l,p)*gfac
                    sumx2 = sumx2 + termx2(k,l,p)*gfac
                    sumy2 = sumy2 + termy2(k,l,p)*gfac
                 enddo
              enddo
           enddo
           
           g0term(iodd,i,j,jp) = g0term(iodd,i,j,jp) + sumg0
           g2term(iodd,i,j,jp) = g2term(iodd,i,j,jp) + sumg2
           g3term(iodd,i,j,jp) = g3term(iodd,i,j,jp) + sumg3
           g4term(iodd,i,j,jp) = g4term(iodd,i,j,jp) + sumg4
           g5term(iodd,i,j,jp) = g5term(iodd,i,j,jp) + sumg5
           g6term(iodd,i,j,jp) = g6term(iodd,i,j,jp) + sumg6
           g7term(iodd,i,j,jp) = g7term(iodd,i,j,jp) + sumg7
           g8term(iodd,i,j,jp) = g8term(iodd,i,j,jp) + sumg8
           g9term(iodd,i,j,jp) = g9term(iodd,i,j,jp) + sumg9
           g10erm(iodd,i,j,jp) = g10erm(iodd,i,j,jp) + sug10
           g11erm(iodd,i,j,jp) = g11erm(iodd,i,j,jp) + sug11
           g12erm(iodd,i,j,jp) = g12erm(iodd,i,j,jp) + sug12
           g13erm(iodd,i,j,jp) = g13erm(iodd,i,j,jp) + sug13
           g14erm(iodd,i,j,jp) = g14erm(iodd,i,j,jp) + sug14
           g15erm(iodd,i,j,jp) = g15erm(iodd,i,j,jp) + sug15              
           g16erm(iodd,i,j,jp) = g16erm(iodd,i,j,jp) + sug16
           h3term(iodd,i,j,jp) = h3term(iodd,i,j,jp) + sumh3
           h5term(iodd,i,j,jp) = h5term(iodd,i,j,jp) + sumh5
           k0term(iodd,i,j,jp) = k0term(iodd,i,j,jp) + sumk0
           k1term(iodd,i,j,jp) = k1term(iodd,i,j,jp) + sumk1
           k2term(iodd,i,j,jp) = k2term(iodd,i,j,jp) + sumk2
           x0term(iodd,i,j,jp) = x0term(iodd,i,j,jp) + sumx0
           y0term(iodd,i,j,jp) = y0term(iodd,i,j,jp) + sumy0
           x1term(iodd,i,j,jp) = x1term(iodd,i,j,jp) + sumx1
           y1term(iodd,i,j,jp) = y1term(iodd,i,j,jp) + sumy1
           x2term(iodd,i,j,jp) = x2term(iodd,i,j,jp) + sumx2
           y2term(iodd,i,j,jp) = y2term(iodd,i,j,jp) + sumy2
           
        enddo ! on jp
        
        if(itor.ne.0) then
        do p=1,20
           sumh2 = 0
           do k=1,20
              do l=1,20
                 sumh2 = sumh2 + &
                      termh(k,l,p)*gtri(k,i,iodd)*gtri(l,j,iodd)
              enddo ! on l
           enddo ! on k
           termh2(iodd,i,j,p) = sumh2
        enddo ! on p
        endif! on itor
        
     enddo ! on j
  enddo ! on i
  
  enddo ! on iodd

  if(itor.ne.0) then
     itri = 0
     do iodd=1,ioddmx
        do ll=1,nreg
           itri = itri + 1
           
           do iii=1,3
              do ii=1,6
                 i = (iii-1)*6 + ii
                 
                 do jjj=1,3
                    do jj=1,6            
                       j = (jjj-1)*6 + jj
                       
                       sum = 0.
                       do jjjp=1,3
                          do jjp=1,6            
                             jp = (jjjp-1)*6 + jjp
                             rinvjps = rinv(6*ist(itri,jjjp)+jjp)
                             do p=1,20
                                gfac = gtri(p,jp,iodd)
                                sum = sum + gfac*termh2(iodd,i,j,p)*rinvjps
                             enddo ! on p
                          enddo ! on jjp
                       enddo ! on jjjp
                       
                       hterm(itri,i,j) = hterm(itri,i,j) + sum
                    enddo !on jj
                 enddo ! on jjj
                 
              enddo ! on ii
           enddo ! on iii
           
        enddo ! on ll
     enddo !on iodd
  endif

  if(myrank.eq.0 .and. iprint.gt.0) then
     print *, "Exiting opdef1."
  endif
  
  return
end subroutine opdef1
!============================================================
subroutine opdef2
  use p_data
  use t_data
  use basic
  use arrays
  implicit none
!
!.....number of 4th rank tensors computed
!
!    index   name      description
!     1      U3   
!     2      U4
!     3      U1
!     4      U2
!     5      V7
!     6      V8
!     7      V9
!     8      V10
!     9      U5
!     10     U6
!     11     V1
!     12     V2
!     13     C1
!     14     C2
!     15     C3
!     16     C4
!     17     C5
!     18     C6
!     19     C7
!     20     C8
!     21     C9
!     22     C10
!     23     C11
!     24     C12
!     25     C13
!     26     C14
!     27     C15
!     28     C16
!     29     C17
!     30     C18
!     31     C19
!     32     C20
!
!
!     37     V3
!     38     V4
!     39     V11
!     40     V12
!     41     V13
!     42     V14
!     43     V15
!     44     V16
!     45     V5
!     46     V6
!     47     U7
!     48     U8
!     49     not used
!     50     not used
!     51     GMU
!     52     GW2
!     53     GW3
!     54     GV11a
!     55     GV12a
!     56     GV12b
!     57     GV12c
!     58     GV13a
!     59     GV13b
!     60     GV13c
!     61     GV22a
!     62     GV23a
!     63     GV23b
!     64     GV33b 
!     65     GV33c 
!     66     GV12ax
!     67     GV12bx
!     68     GV12cx
!     69     GV13bx
!     70     GV13cx
!     71     GV23ax
!     72     GV23bx
!     73     T1
!     74     T2
!     75     T3
!     76     T4
!     77     T5
!     78     W2
!     79     W3


  integer :: iodd, i, j, k, l, itype, msum, nsum
  real :: co, sn, co2, sn2, sum2, gq, gr
  real :: f
  allocate (tensor(20,20,20,20),tterm(18,18,18,18))

  do iodd=1,ioddmx
     
  co = cos(ttri(iodd))
  sn = sin(ttri(iodd))
  co2 = cos(2.*ttri(iodd))
  sn2 = sin(2.*ttri(iodd))
!
!     calculate matrix elements of the integration matrix fint
  do i=-6,maxi
     do j=-6,maxi
        fint(i,j) = 0
        if(i.ge.0 .and.j.ge.0)                                      &
             fint(i,j) = f(i,j,atri(iodd),btri(iodd),ctri(iodd))
     enddo
  enddo

  do itype=1,ntensor
     
     if(iprint.eq.1) write(*,*) "myrank,iodd,itype,iread",            &
          myrank,iodd,itype,iread(itype)
     if(iread(itype).ne.1) go to 400
         
      ! begin local integrations of metric elements for 
      ! 4th rank tensor terms
      
     if(itype.gt.ntensor) go to 200
     go to(101,102,103,104,105, 106,107,108,109,110,                   &
           111,112,113,114,115, 116,117,118,119,120,                   &
           121,122,123,124,125, 126,127,128,129,130,                   &
           131,132,133,134,135, 136,137,138,139,140,                   &
           141,142,143,144,145, 146,147,148,149,150,                   &
           151,152,153,154,155, 156,157,158,159,160,                   &
           161,162,163,164,165, 166,167,168,169,170,                   &
           171,172,173,174,175, 176,177,178,179),itype
!
!     U3
  101 continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = (mi(q)*ni(r)-mi(r)*ni(q))               &
                                *fint(msum-1,nsum-1)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     U4 = U3 (q <-> r)
  102 continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = (mi(r)*ni(q)-mi(q)*ni(r))               &
                               *fint(msum-1,nsum-1)
            enddo
          enddo
        enddo
      enddo
!
      go to 200
!
!     U1
  103 continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                           &
              -(mi(p)*ni(r)-mi(r)*ni(p))                                  &
                 *(mi(q)*(mi(q)-1)*fint(msum-3,nsum-1)                    &
                  +ni(q)*(ni(q)-1)*fint(msum-1,nsum-3))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     U2 = U1 (q <-> r)
 104  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                           &
              -(mi(p)*ni(q)-mi(q)*ni(p))                                  &
                 *(mi(r)*(mi(r)-1)*fint(msum-3,nsum-1)                    &
                  +ni(r)*(ni(r)-1)*fint(msum-1,nsum-3))
            enddo
          enddo
        enddo
      enddo
!
      go to 200
!
!     V7
 105  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = (mi(q)*ni(r)-mi(r)*ni(q)) &
                                * fint(msum-1,nsum-1)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     V8 = V7 (q <-> r)
 106  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = (mi(r)*ni(q)-mi(q)*ni(r)) &
                                * fint(msum-1,nsum-1)
            enddo
          enddo
        enddo
      enddo
!
      go to 200
!
!     V9
 107  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = mi(q)*mi(r)*fint(msum-2,nsum)           &
                              + ni(q)*ni(r)*fint(msum,nsum-2)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     V10 = V9(q <-> r)
 108  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = mi(r)*mi(q)*fint(msum-2,nsum)           &
                              + ni(r)*ni(q)*fint(msum,nsum-2)
            enddo
          enddo
        enddo
      enddo
!
      go to 200
!
!     U5
 109  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s)=((mi(q)*ni(r)-mi(r)*ni(q))                 &
                 + gam*(mi(s)*ni(r)-mi(r)*ni(s)))*fint(msum-1,nsum-1)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     U6 = U5(q <-> r)
 110  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s)=((mi(r)*ni(q)-mi(q)*ni(r))                 &
                 + gam*(mi(s)*ni(q)-mi(q)*ni(s)))*fint(msum-1,nsum-1)
            enddo
          enddo
        enddo
      enddo
!
      go to 200
!
!     V1
 111  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s)=                                          &
              mi(q)*(mi(r)*((mi(q)-1)*(ni(s)+ni(p))-ni(q)*mi(s))        &
                   - ni(r)*mi(p)*(mi(q)-1))*fint(msum-3,nsum-1)         &
             +ni(q)*(ni(r)*(mi(q)*ni(s)-(ni(q)-1)*(mi(s)+mi(p)))        &
                   + mi(r)*ni(p)*(ni(q)-1))*fint(msum-1,nsum-3)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     V2   (note:   V2 = V1(q <-> r)
 112  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s)=                                          &
              mi(r)*(mi(q)*((mi(r)-1)*(ni(s)+ni(p))-ni(r)*mi(s))        &
                   - ni(q)*mi(p)*(mi(r)-1))*fint(msum-3,nsum-1)         &
             +ni(r)*(ni(q)*(mi(r)*ni(s)-(ni(r)-1)*(mi(s)+mi(p)))        &
                   + mi(q)*ni(p)*(ni(r)-1))*fint(msum-1,nsum-3)
            enddo
          enddo
        enddo
      enddo
!
      go to 200
!
!      C1
  113  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = (mi(s)*ni(q)-mi(q)*ni(s))                      &
                   *((mi(r)*(mi(r)-1)                                          &
                   *((mi(q)+mi(s)-1)*ni(p)-(ni(q)+ni(s)-1)*mi(p))              &
                   +(mi(p)*ni(r)-mi(r)*ni(p))*(mi(s)+mi(q)-1)*(mi(p)+mi(r)-1)) &
                   *fint(msum-4,nsum-2)                                        &
                   +(ni(r)*(ni(r)-1)                                           &
                   *((mi(q)+mi(s)-1)*ni(p)-(ni(q)+ni(s)-1)*mi(p))              &
                   +(mi(p)*ni(r)-mi(r)*ni(p))*(ni(s)+ni(q)-1)*(ni(p)+ni(r)-1)) &
                   *fint(msum-2,nsum-4))
            enddo
          enddo
        enddo
      enddo
       go to 200
!
!     C2
  114  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = (mi(s)*ni(r)-mi(r)*ni(s))                      &
                   *((mi(q)*(mi(q)-1)                                          &
                   *((mi(r)+mi(s)-1)*ni(p)-(ni(r)+ni(s)-1)*mi(p))              &
                   +(mi(p)*ni(q)-mi(q)*ni(p))*(mi(s)+mi(r)-1)*(mi(p)+mi(q)-1)) &
                   *fint(msum-4,nsum-2)                                        &
                   +(ni(q)*(ni(q)-1)                                           &
                   *((mi(r)+mi(s)-1)*ni(p)-(ni(r)+ni(s)-1)*mi(p))              &
                   +(mi(p)*ni(q)-mi(q)*ni(p))*(ni(s)+ni(r)-1)*(ni(p)+ni(q)-1)) &
                   *fint(msum-2,nsum-4))
            enddo
          enddo
        enddo
      enddo
!
       go to 200
!
!     C3
  115  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =   mi(r)*mi(q)                                 &
              *((mi(p)*ni(s)-mi(s)*ni(p))*((mi(r)+mi(q)-2)*(mi(p)+mi(s)-1))   &
                +mi(s)*(mi(s)-1)*((mi(r)+mi(q)-2)*ni(p)-mi(p)*(ni(r)+ni(q)))) &
                                 *fint(msum-5,nsum-1)                         &
              + (mi(r)*mi(q)                                                  &
              *((mi(p)*ni(s)-mi(s)*ni(p))*((ni(r)+ni(q))*(ni(p)+ni(s)-1))     &
                +ni(s)*(ni(s)-1)*((mi(r)+mi(q)-2)*ni(p)-mi(p)*(ni(r)+ni(q)))) &   
                +ni(r)*ni(q)                                                  &
              *((mi(p)*ni(s)-mi(s)*ni(p))*((mi(r)+mi(q))*(mi(p)+mi(s)-1))     &
                +mi(s)*(mi(s)-1)*((mi(r)+mi(q))*ni(p)-mi(p)*(ni(r)+ni(q)-2))))&
                                 *fint(msum-3,nsum-3)                         &
              + ni(r)*ni(q)                                                   &
              *((mi(p)*ni(s)-mi(s)*ni(p))*((ni(r)+ni(q)-2)*(ni(p)+ni(s)-1))   &
                +ni(s)*(ni(s)-1)*((mi(r)+mi(q))*ni(p)-mi(p)*(ni(r)+ni(q)-2))) &
                                 *fint(msum-1,nsum-5)                         
            enddo
          enddo
        enddo
      enddo
       go to 200
!
!     C4
 116  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =   mi(q)*mi(r)                                 &
              *((mi(p)*ni(s)-mi(s)*ni(p))*((mi(q)+mi(r)-2)*(mi(p)+mi(s)-1))   &
                +mi(s)*(mi(s)-1)*((mi(q)+mi(r)-2)*ni(p)-mi(p)*(ni(q)+ni(r)))) &
                                 *fint(msum-5,nsum-1)                         &
              + (mi(q)*mi(r)                                                  &
              *((mi(p)*ni(s)-mi(s)*ni(p))*((ni(q)+ni(r))*(ni(p)+ni(s)-1))     &
                +ni(s)*(ni(s)-1)*((mi(q)+mi(r)-2)*ni(p)-mi(p)*(ni(q)+ni(r)))) &
                +ni(q)*ni(r)                                                  &
              *((mi(p)*ni(s)-mi(s)*ni(p))*((mi(q)+mi(r))*(mi(p)+mi(s)-1))     &
                +mi(s)*(mi(s)-1)*((mi(q)+mi(r))*ni(p)-mi(p)*(ni(q)+ni(r)-2))))&
                                 *fint(msum-3,nsum-3)                         &
              + ni(q)*ni(r)                                                   &
              *((mi(p)*ni(s)-mi(s)*ni(p))*((ni(q)+ni(r)-2)*(ni(p)+ni(s)-1))   &
                +ni(s)*(ni(s)-1)*((mi(q)+mi(r))*ni(p)-mi(p)*(ni(q)+ni(r)-2))) &
                                 *fint(msum-1,nsum-5)                         
            enddo
          enddo
        enddo
      enddo
!
      go to 200
!
!     C5
 117  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         &
            ((mi(q)*ni(r)-mi(r)*ni(q))*(mi(p)*ni(s)-mi(s)*ni(p))        &
            -(mi(q)*ni(s)-mi(s)*ni(q))*(mi(p)*ni(r)-mi(r)*ni(p)))       &
                              *fint(msum-2,nsum-2)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     C6
 118  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         &
            ((mi(r)*ni(q)-mi(q)*ni(r))*(mi(p)*ni(s)-mi(s)*ni(p))        &
            -(mi(r)*ni(s)-mi(s)*ni(r))*(mi(p)*ni(q)-mi(q)*ni(p)))       &
                              *fint(msum-2,nsum-2)
            enddo
          enddo
        enddo
      enddo
!
      go to 200
!
!     C7
 119  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = -mi(q)*(mi(r)*(mi(p)*ni(s)-mi(s)*ni(p)) &
               + (mi(r)*ni(p) - mi(p)*ni(r))*(mi(s)+mi(q)-1))           &
                              *fint(msum-3,nsum-1)                      &
                                -ni(q)*(ni(r)*(mi(p)*ni(s)-mi(s)*ni(p)) &
               + (mi(r)*ni(p) - mi(p)*ni(r))*(ni(s)+ni(q)-1))           &
                              *fint(msum-1,nsum-3)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     C8
 120  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = -mi(r)*(mi(q)*(mi(p)*ni(s)-mi(s)*ni(p)) &
               + (mi(q)*ni(p) - mi(p)*ni(q))*(mi(s)+mi(r)-1))           &
                              *fint(msum-3,nsum-1)                      &
                                -ni(r)*(ni(q)*(mi(p)*ni(s)-mi(s)*ni(p)) &
               + (mi(q)*ni(p) - mi(p)*ni(q))*(ni(s)+ni(r)-1))           &
                              *fint(msum-1,nsum-3)
            enddo
          enddo
        enddo
      enddo
!
      go to 200
!
!     C9
 121  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = -0.5*((mi(r)*ni(p)-mi(p)*ni(r))         &
                                     *(mi(s)*ni(q)-mi(q)*ni(s))         &
                                     +(mi(s)*ni(p)-mi(p)*ni(s))         &
                                     *(mi(r)*ni(q)-mi(q)*ni(r)))        &
                              *fint(msum-2,nsum-2)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     C10
 122  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = -0.5*((mi(q)*ni(p)-mi(p)*ni(q))         &
                                     *(mi(s)*ni(r)-mi(r)*ni(s))         &
                                     +(mi(s)*ni(p)-mi(p)*ni(s))         &
                                     *(mi(q)*ni(r)-mi(r)*ni(q)))        &
                              *fint(msum-2,nsum-2)
            enddo
          enddo
        enddo
      enddo
!
      go to 200
!
!     C11
 123  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = -(mi(s)*ni(q)-mi(q)*ni(s)) *                           &
                   (mi(p)*mi(r)*(mi(p)-1)*(mi(s)+mi(q)-1)*fint(msum-5,nsum-1) +        &
                    ni(p)*ni(r)*(ni(p)-1)*(ni(s)+ni(q)-1)*fint(msum-1,nsum-5) +        &
                   ((ni(p)*ni(r)*(mi(p)+mi(r))-mi(p)*ni(r)*(ni(r)-1))*(mi(s)+mi(q)-1)  &
                   +(mi(p)*mi(r)*(ni(p)+ni(r))-ni(p)*mi(r)*(mi(r)-1))*(ni(s)+ni(q)-1)) &
                   * fint(msum-3,nsum-3))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     C12
 124  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = -(mi(s)*ni(r)-mi(r)*ni(s)) *                           &
                   (mi(p)*mi(q)*(mi(p)-1)*(mi(s)+mi(r)-1)*fint(msum-5,nsum-1) +        &
                    ni(p)*ni(q)*(ni(p)-1)*(ni(s)+ni(r)-1)*fint(msum-1,nsum-5) +        &
                   ((ni(p)*ni(q)*(mi(p)+mi(q))-mi(p)*ni(q)*(ni(q)-1))*(mi(s)+mi(r)-1)  &
                   +(mi(p)*mi(q)*(ni(p)+ni(q))-ni(p)*mi(q)*(mi(q)-1))*(ni(s)+ni(r)-1)) &
                   * fint(msum-3,nsum-3))
            enddo
          enddo
        enddo
      enddo
!
      go to 200
!
!     C13
 125  continue
       do p=1,20
         do q=1,20
           do r=1,20
             do s=1,20
               msum = mi(p)+mi(q)+mi(r)+mi(s)
               nsum = ni(p)+ni(q)+ni(r)+ni(s)
                tensor(p,q,r,s) =                                                                    &
                       -(mi(p)-1)*mi(s)*mi(p)*mi(r)*mi(q)*(mi(r)+mi(q)-2)*fint(msum-6,nsum)          &
                +(((mi(s)-1)*ni(r)*ni(q)-mi(r)*mi(q)*(ni(p)+ni(s)))*mi(s)*mi(p)*(ni(r)+ni(q))        &
                  +((ni(s)-ni(p)-1)*mi(p)-ni(p)*mi(s))*ni(s)*mi(r)*mi(q)*(mi(r)+mi(q)-2)             &
               +((mi(s)-1)*ni(p)*mi(r)*mi(q)-ni(r)*ni(q)*mi(p)*(mi(p)+mi(s)-2))*mi(s)*(mi(r)+mi(q))) &
                                       *fint(msum-4,nsum-2)                                          &
                +(((ni(s)-1)*mi(r)*mi(q)-ni(r)*ni(q)*(mi(p)+mi(s)))*ni(s)*ni(p)*(mi(r)+mi(q))        &
                  +((mi(s)-mi(p)-1)*ni(p)-mi(p)*ni(s))*mi(s)*ni(r)*ni(q)*(ni(r)+ni(q)-2)             &
               +((ni(s)-1)*mi(p)*ni(r)*ni(q)-mi(r)*mi(q)*ni(p)*(ni(p)+ni(s)-2))*ni(s)*(ni(r)+ni(q))) &
                                       *fint(msum-2,nsum-4)                                          &
                  -(ni(p)-1)*ni(s)*ni(p)*ni(r)*ni(q)*(ni(r)+ni(q)-2)*fint(msum,nsum-6)                
             enddo
           enddo
         enddo
       enddo
      go to 200
!
!     C14
 126  continue
       do p=1,20
         do q=1,20
           do r=1,20
             do s=1,20
               msum = mi(p)+mi(q)+mi(r)+mi(s)
               nsum = ni(p)+ni(q)+ni(r)+ni(s)
                tensor(p,q,r,s) =                                                                    &
                       -(mi(p)-1)*mi(s)*mi(p)*mi(q)*mi(r)*(mi(q)+mi(r)-2)*fint(msum-6,nsum)          &
                +(((mi(s)-1)*ni(q)*ni(r)-mi(q)*mi(r)*(ni(p)+ni(s)))*mi(s)*mi(p)*(ni(q)+ni(r))        &
                  +((ni(s)-ni(p)-1)*mi(p)-ni(p)*mi(s))*ni(s)*mi(q)*mi(r)*(mi(q)+mi(r)-2)             &
               +((mi(s)-1)*ni(p)*mi(q)*mi(r)-ni(q)*ni(r)*mi(p)*(mi(p)+mi(s)-2))*mi(s)*(mi(q)+mi(r))) &
                                       *fint(msum-4,nsum-2)                                          &
                +(((ni(s)-1)*mi(q)*mi(r)-ni(q)*ni(r)*(mi(p)+mi(s)))*ni(s)*ni(p)*(mi(q)+mi(r))        &
                  +((mi(s)-mi(p)-1)*ni(p)-mi(p)*ni(s))*mi(s)*ni(q)*ni(r)*(ni(q)+ni(r)-2)             &
               +((ni(s)-1)*mi(p)*ni(q)*ni(r)-mi(q)*mi(r)*ni(p)*(ni(p)+ni(s)-2))*ni(s)*(ni(q)+ni(r))) &
                                       *fint(msum-2,nsum-4)                                          &
                  -(ni(p)-1)*ni(s)*ni(p)*ni(q)*ni(r)*(ni(q)+ni(r)-2)*fint(msum,nsum-6)
             enddo
           enddo
         enddo
       enddo
!
      go to 200
!
!     C15
 127  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =-(mi(s)*ni(q)-ni(s)*mi(q))               &
               *( mi(p)*(mi(p)-1)*fint(msum-3,nsum-1)                   &
                 +ni(p)*(ni(p)-1)*fint(msum-1,nsum-3))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     C16
 128  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =-(mi(s)*ni(r)-ni(s)*mi(r))               &
               *( mi(p)*(mi(p)-1)*fint(msum-3,nsum-1)                   &
                 +ni(p)*(ni(p)-1)*fint(msum-1,nsum-3))
            enddo
          enddo
        enddo
      enddo
!
      go to 200
!
!      C17
 129  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         &
               -mi(p)*mi(q)*(mi(p)-1)*(mi(q)+mi(s)-1)*fint(msum-4,nsum) &
           -(mi(p)*ni(q)*(mi(p)-1)*(ni(q)+ni(s)-1)                      &
            +ni(p)*mi(q)*(ni(p)-1)*(mi(q)+mi(s)-1))*fint(msum-2,nsum-2) &
               -ni(p)*ni(q)*(ni(p)-1)*(ni(q)+ni(s)-1)*fint(msum,nsum-4)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     C18
 130  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         &
               -mi(p)*mi(r)*(mi(p)-1)*(mi(r)+mi(s)-1)*fint(msum-4,nsum) &
           -(mi(p)*ni(r)*(mi(p)-1)*(ni(r)+ni(s)-1)                      &
            +ni(p)*mi(r)*(ni(p)-1)*(mi(r)+mi(s)-1))*fint(msum-2,nsum-2) &
               -ni(p)*ni(r)*(ni(p)-1)*(ni(r)+ni(s)-1)*fint(msum,nsum-4)
            enddo
          enddo
        enddo
      enddo
!
      go to 200
!
!     C19
 131  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =-(mi(r)*ni(q)-mi(q)*ni(r))               &
                     *(mi(p)*(mi(p)-1)*fint(msum-3,nsum-1)              &
                      +ni(p)*(ni(p)-1)*fint(msum-1,nsum-3))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     C20
 132  continue
!
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =-(mi(q)*ni(r)-mi(r)*ni(q))               &
                     *(mi(p)*(mi(p)-1)*fint(msum-3,nsum-1)              &
                      +ni(p)*(ni(p)-1)*fint(msum-1,nsum-3))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
 133  continue
      go to 200
!
!     NOT USED
 134  continue
!
      go to 200
!
 135  continue
      go to 200
!
!     NOT USED
 136  continue
!
      go to 200
!
!     V3
 137  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         &
               (-mi(r)*mi(p)*mi(q)*(mi(q)-1))*fint(msum-4,nsum)         &
              +(mi(r)*ni(q)*((mi(r)+mi(q)-1)*ni(s)                      &
                   -(ni(q)-1)*(mi(s)+mi(p)) - mi(s)*ni(r))              &
               +ni(r)*mi(q)*((ni(r)+ni(q)-1)*mi(s)                      &
                   -(mi(q)-1)*(ni(s)+ni(p)) - ni(s)*mi(r)))             &
                                             *fint(msum-2,nsum-2)       &
              +(-ni(r)*ni(p)*ni(q)*(ni(q)-1))*fint(msum,nsum-4)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     V4 = V3(q <-> r)
 138  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         &
               (-mi(q)*mi(p)*mi(r)*(mi(r)-1))*fint(msum-4,nsum)         &
              +(mi(q)*ni(r)*((mi(q)+mi(r)-1)*ni(s)                      &
                   -(ni(r)-1)*(mi(s)+mi(p)) - mi(s)*ni(q))              &
               +ni(q)*mi(r)*((ni(q)+ni(r)-1)*mi(s)                      &
                   -(mi(r)-1)*(ni(s)+ni(p)) - ni(s)*mi(q)))             &
                                             *fint(msum-2,nsum-2)       &
              +(-ni(q)*ni(p)*ni(r)*(ni(r)-1))*fint(msum,nsum-4)
            enddo
          enddo
        enddo
      enddo
!
      go to 200
!
!     V11
 139  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         &
              mi(q)*((mi(q)-1)*(mi(r)*(ni(p)+ni(s)+2*ni(r))             &
              -mi(p)*ni(r))+mi(r)*(mi(s)*(ni(r)-ni(q))                  &
                   -(mi(r)-1)*(ni(s)+2*ni(q))))*fint(msum-3,nsum-1)     &
             -ni(q)*((ni(q)-1)*(ni(r)*(mi(p)+mi(s)+2*mi(r))             &
              -ni(p)*mi(r))+ni(r)*(ni(s)*(mi(r)-mi(q))                  &
                   -(ni(r)-1)*(mi(s)+2*mi(q))))*fint(msum-1,nsum-3)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     V12 = V11 (q <-> r)
 140  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         &
              mi(r)*((mi(r)-1)*(mi(q)*(ni(p)+ni(s)+2*ni(q))             &
              -mi(p)*ni(q))+mi(q)*(mi(s)*(ni(q)-ni(r))                  &
                   -(mi(q)-1)*(ni(s)+2*ni(r))))*fint(msum-3,nsum-1)     &
             -ni(r)*((ni(r)-1)*(ni(q)*(mi(p)+mi(s)+2*mi(q))             &
              -ni(p)*mi(q))+ni(q)*(ni(s)*(mi(q)-mi(r))                  &
                   -(ni(q)-1)*(mi(s)+2*mi(r))))*fint(msum-1,nsum-3)
            enddo
          enddo
        enddo
      enddo
!
      go to 200
!
!     V13
 141  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         &
              .5*mi(p)*mi(q)*mi(r)*(mi(p)+mi(s)-1)*fint(msum-4,nsum)    &
            +(.5*mi(p)*ni(q)*ni(r)*(mi(p)+mi(s)-1)                      &
            + .5*ni(p)*mi(q)*mi(r)*(ni(p)+ni(s)-1))*fint(msum-2,nsum-2) &
            + .5*ni(p)*ni(q)*ni(r)*(ni(p)+ni(s)-1)*fint(msum,nsum-4) 
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     V14 = V13 (q <-> r)
 142  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         &
     &        .5*mi(p)*mi(r)*mi(q)*(mi(p)+mi(s)-1)*fint(msum-4,nsum)    &
     &      +(.5*mi(p)*ni(r)*ni(q)*(mi(p)+mi(s)-1)                      &
     &      + .5*ni(p)*mi(r)*mi(q)*(ni(p)+ni(s)-1))*fint(msum-2,nsum-2) &
     &      + .5*ni(p)*ni(r)*ni(q)*(ni(p)+ni(s)-1)*fint(msum,nsum-4)
            enddo
          enddo
        enddo
      enddo
!
      go to 200
!
!     V15
 143  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         &
                (2.*mi(q)*ni(r)*(ni(q)*mi(r) - (mi(q)-1)*(ni(r)-1))     &
                +ni(q)*mi(r)*(mi(s)*ni(r)-(mi(r)-1)*ni(s))              &
                +mi(q)*ni(r)*(mi(r)*ni(s)-(ni(r)-1)*mi(s)))             &
                              *fint(msum-2,nsum-2)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     V16 = V15(q <-> r)
 144  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         &
                (2.*mi(r)*ni(q)*(ni(r)*mi(q) - (mi(r)-1)*(ni(q)-1))     &
                +ni(r)*mi(q)*(mi(s)*ni(q)-(mi(q)-1)*ni(s))              &
                +mi(r)*ni(q)*(mi(q)*ni(s)-(ni(q)-1)*mi(s)))             &
                              *fint(msum-2,nsum-2)
            enddo
          enddo
        enddo
      enddo
!
      go to 200
!
!     V5
 145  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = 0.5*(mi(s)*ni(p)-mi(p)*ni(s))           &
                *(mi(r)*mi(q)*fint(msum-3,nsum-1)                       &
                 +ni(r)*ni(q)*fint(msum-1,nsum-3))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     V6
 146  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = 0.5*(mi(s)*ni(p)-mi(p)*ni(s))           &
                *(mi(q)*mi(r)*fint(msum-3,nsum-1)                       &
                 +ni(q)*ni(r)*fint(msum-1,nsum-3))
            enddo
          enddo
        enddo
      enddo
!
      go to 200
!
!     U7
 147  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =  (mi(s)*ni(r)-mi(r)*ni(s))              &
                                                *fint(msum-1,nsum-1)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     U8 = U7(q <->r)
 148  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =  (mi(s)*ni(q)-mi(q)*ni(s))              &
                                                *fint(msum-1,nsum-1)         
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     not used
 149  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         &
                 ((mi(q)*ni(r)-mi(r)*ni(q))*(mi(p)*ni(s)-mi(s)*ni(p))   &
                 +(mi(s)*ni(r)-mi(r)*ni(s))*(mi(p)*ni(q)-mi(q)*ni(p)))  &
                          *fint(msum-2,nsum-2)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     not used
 150  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         &
               (mi(q)*(mi(q)-1)*((mi(s)*ni(r)-mi(r)*ni(s))              &
                 *((mi(r)+mi(s)-1)*ni(p)-(ni(r)+ni(s)-1)*mi(p))         &
                                -(mi(s)*ni(p)-mi(p)*ni(s))              &
                 *((mi(p)+mi(s)-1)*ni(r)-(ni(p)+ni(s)-1)*mi(r)))        &
               -mi(r)*(mi(r)-1)*(mi(s)*ni(p)-mi(p)*ni(s))               &
                     *(mi(q)*(ni(p)+ni(s)-1)-ni(q)*(mi(p)+mi(s)-1))     &
               +2*(mi(s)*ni(p)-mi(p)*ni(s))*mi(q)*mi(r)                 &
                 *((mi(q)-1)*ni(r)-ni(q)*(mi(r)-1)))*fint(msum-4,nsum-2)&
!
              +(ni(q)*(ni(q)-1)*((mi(s)*ni(r)-mi(r)*ni(s))              &
                 *((mi(r)+mi(s)-1)*ni(p)-(ni(r)+ni(s)-1)*mi(p))         &
                                -(mi(s)*ni(p)-mi(p)*ni(s))              &
                 *((mi(p)+mi(s)-1)*ni(r)-(ni(p)+ni(s)-1)*mi(r)))        &
               -ni(r)*(ni(r)-1)*(mi(s)*ni(p)-mi(p)*ni(s))               &
                     *(mi(q)*(ni(p)+ni(s)-1)-ni(q)*(mi(p)+mi(s)-1))     &
               +2*(mi(s)*ni(p)-mi(p)*ni(s))*ni(q)*ni(r)                 &
                 *(mi(q)*(ni(r)-1)-(ni(q)-1)*mi(r)))*fint(msum-2,nsum-4)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     GMU
 151  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = fint(msum,nsum)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     GW2
 152  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         &
                   co2*(mi(q)*mi(r)*fint(msum-2,nsum)                   &
                       -ni(q)*ni(r)*fint(msum,nsum-2))                  &
                  -sn2*(mi(q)*ni(r)+mi(r)*ni(q))*fint(msum-1,nsum-1)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     GW3
 153  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         &
                   sn2*(mi(q)*mi(r)*fint(msum-2,nsum)                   &
                       -ni(q)*ni(r)*fint(msum,nsum-2))                  &
                  +co2*(mi(q)*ni(r)+mi(r)*ni(q))*fint(msum-1,nsum-1)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     GV11a
 154  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         & 
                   mi(p)*mi(q)*((mi(q)-1)*ni(p)-(mi(p)-1)*ni(q))        &
                   *fint(msum-3,nsum-1)                                 &
                  -ni(p)*ni(q)*((ni(q)-1)*mi(p)-(ni(p)-1)*mi(q))        &
                   *fint(msum-1,nsum-3)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     GV12a
 155  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = 0.5*                                    & 
                   ((0.5*mi(p)*(mi(p)-1)*(mi(r)*ni(q)+mi(q)*ni(r))      &
                        -mi(p)*mi(q)*mi(r)*ni(p))*fint(msum-3,nsum-1)   &
                   -(0.5*ni(p)*(ni(p)-1)*(ni(r)*mi(q)+ni(q)*mi(r))      &
                        -ni(p)*ni(q)*ni(r)*mi(p))*fint(msum-1,nsum-3))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     GV12b
 156  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         &
                   (mi(r)*ni(q)-mi(q)*ni(r))                            &
                   *(0.5*co2*(mi(p)*(mi(p)-1)*fint(msum-3,nsum-1)       &
                             -ni(p)*(ni(p)-1)*fint(msum-1,nsum-3))      &
                    -sn2*mi(p)*ni(p)*fint(msum-2,nsum-2))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     GV12c
 157  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         &
                   (mi(r)*ni(q)-mi(q)*ni(r))                            &
                   *(0.5*sn2*(mi(p)*(mi(p)-1)*fint(msum-3,nsum-1)       &
                             -ni(p)*(ni(p)-1)*fint(msum-1,nsum-3))      &
                    +co2*mi(p)*ni(p)*fint(msum-2,nsum-2))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     GV13a
 158  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = -0.5*                                   & 
                   (mi(p)*(mi(p)-1)*mi(q)*(mi(q)-1)*fint(msum-4,nsum)   &
                   +ni(p)*(ni(p)-1)*ni(q)*(ni(q)-1)*fint(msum,nsum-4)   &
                   +(4.*mi(p)*mi(q)*ni(p)*ni(q)                         &
                    -(mi(p)*(mi(p)-1)*ni(q)*(ni(q)-1)                   &
                     +ni(p)*(ni(p)-1)*mi(q)*(mi(q)-1)))                 &
                    *fint(msum-2,nsum-2))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     GV13b
 159  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         & 
                   0.5*co2*(mi(p)*(mi(p)-1)*                            &
                            (mi(q)*(mi(q)-1)*fint(msum-4,nsum)          &
                            +ni(q)*(ni(q)-1)*fint(msum-2,nsum-2))       &
                           -ni(p)*(ni(p)-1)*                            &
                            (mi(q)*(mi(q)-1)*fint(msum-2,nsum-2)        &
                            +ni(q)*(ni(q)-1)*fint(msum,nsum-4)))        &
                  -sn2*mi(p)*ni(p)*(mi(q)*(mi(q)-1)*fint(msum-3,nsum-1) &
                                   +ni(q)*(ni(q)-1)*fint(msum-1,nsum-3))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     GV13c
 160  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         &
                   0.5*sn2*(mi(p)*(mi(p)-1)*                            &
                            (mi(q)*(mi(q)-1)*fint(msum-4,nsum)          &
                            +ni(q)*(ni(q)-1)*fint(msum-2,nsum-2))       &
                           -ni(p)*(ni(p)-1)*                            &
                            (mi(q)*(mi(q)-1)*fint(msum-2,nsum-2)        &
                            +ni(q)*(ni(q)-1)*fint(msum,nsum-4)))        &
                  +co2*mi(p)*ni(p)*(mi(q)*(mi(q)-1)*fint(msum-3,nsum-1) &
                                   +ni(q)*(ni(q)-1)*fint(msum-1,nsum-3))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     GV22a
 161  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = 0.25*(mi(p)*ni(q)-mi(q)*ni(p))          &
                                *fint(msum-1,nsum-1)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     GV23a
 162  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = -0.5*                                   & 
                  (mi(p)*mi(q)*(mi(q)-1)*mi(r)*fint(msum-4,nsum)        &
                  +ni(p)*ni(q)*(ni(q)-1)*ni(r)*fint(msum,nsum-4)        &
                  +mi(q)*ni(q)*(mi(r)*ni(p)+mi(p)*ni(r))                &
                   *fint(msum-2,nsum-2))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     GV23b
 163  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = 0.5*                                    & 
                   (mi(q)*ni(q)*(mi(r)*ni(p)+mi(p)*ni(r))               &
                   -(mi(p)*mi(r)*ni(q)*(ni(q)-1)                        &
                    +ni(p)*ni(r)*mi(q)*(mi(q)-1)))                      &
                   *fint(msum-2,nsum-2)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     GV33b
 164  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         & 
                  co2*(mi(p)*mi(q)*((mi(p)-1)*ni(q)-(mi(q)-1)*ni(p))    &
                       *fint(msum-3,nsum-1)                             &
                      +ni(p)*ni(q)*((ni(p)-1)*mi(q)-(ni(q)-1)*mi(p))    &
                       *fint(msum-1,nsum-3))                            &
                 -sn2*(mi(p)*(mi(p)-1)*ni(q)*(ni(q)-1)                  &
                      -ni(p)*(ni(p)-1)*mi(q)*(mi(q)-1))                 &
                  *fint(msum-2,nsum-2)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     GV33c
 165  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         & 
                  sn2*(mi(p)*mi(q)*((mi(p)-1)*ni(q)-(mi(q)-1)*ni(p))    &
                       *fint(msum-3,nsum-1)                             &
                      +ni(p)*ni(q)*((ni(p)-1)*mi(q)-(ni(q)-1)*mi(p))    &
                       *fint(msum-1,nsum-3))                            &
                 +co2*(mi(p)*(mi(p)-1)*ni(q)*(ni(q)-1)                  &
                      -ni(p)*(ni(p)-1)*mi(q)*(mi(q)-1))                 &
                  *fint(msum-2,nsum-2)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     GV12ax
 166  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = 0.5*                                    & 
                   ((0.5*mi(q)*(mi(q)-1)*(mi(r)*ni(p)+mi(p)*ni(r))      &
                        -mi(q)*mi(p)*mi(r)*ni(q))*fint(msum-3,nsum-1)   &
                   -(0.5*ni(q)*(ni(q)-1)*(ni(r)*mi(p)+ni(p)*mi(r))      &
                        -ni(q)*ni(p)*ni(r)*mi(q))*fint(msum-1,nsum-3))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     GV12bx
 167  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         &
                   (mi(r)*ni(p)-mi(p)*ni(r))                            &
                   *(0.5*co2*(mi(q)*(mi(q)-1)*fint(msum-3,nsum-1)       &
                             -ni(q)*(ni(q)-1)*fint(msum-1,nsum-3))      &
                    -sn2*mi(q)*ni(q)*fint(msum-2,nsum-2))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     GV12cx
 168  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         &
                   (mi(r)*ni(p)-mi(p)*ni(r))                            &
                   *(0.5*sn2*(mi(q)*(mi(q)-1)*fint(msum-3,nsum-1)       &
                             -ni(q)*(ni(q)-1)*fint(msum-1,nsum-3))      &
                    +co2*mi(q)*ni(q)*fint(msum-2,nsum-2))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     GV13bx
 169  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         & 
                   0.5*co2*(mi(q)*(mi(q)-1)*                            &
                            (mi(p)*(mi(p)-1)*fint(msum-4,nsum)          &
                            +ni(p)*(ni(p)-1)*fint(msum-2,nsum-2))       &
                           -ni(q)*(ni(q)-1)*                            &
                            (mi(p)*(mi(p)-1)*fint(msum-2,nsum-2)        &
                            +ni(p)*(ni(p)-1)*fint(msum,nsum-4)))        &
                  -sn2*mi(q)*ni(q)*(mi(p)*(mi(p)-1)*fint(msum-3,nsum-1) &
                                   +ni(p)*(ni(p)-1)*fint(msum-1,nsum-3))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     GV13cx
 170  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) =                                         & 
                   0.5*sn2*(mi(q)*(mi(q)-1)*                            &
                            (mi(p)*(mi(p)-1)*fint(msum-4,nsum)          &
                            +ni(p)*(ni(p)-1)*fint(msum-2,nsum-2))       &
                           -ni(q)*(ni(q)-1)*                            &
                            (mi(p)*(mi(p)-1)*fint(msum-2,nsum-2)        &
                            +ni(p)*(ni(p)-1)*fint(msum,nsum-4)))        &
                  +co2*mi(q)*ni(q)*(mi(p)*(mi(p)-1)*fint(msum-3,nsum-1) &
                                   +ni(p)*(ni(p)-1)*fint(msum-1,nsum-3))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     GV23ax
 171  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = -0.5*                                   & 
                  (mi(q)*mi(p)*(mi(p)-1)*mi(r)*fint(msum-4,nsum)        &
                  +ni(q)*ni(p)*(ni(p)-1)*ni(r)*fint(msum,nsum-4)        &
                  +mi(p)*ni(p)*(mi(r)*ni(q)+mi(q)*ni(r))                &
                   *fint(msum-2,nsum-2))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     GV23bx
 172  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = 0.5*                                    & 
                   (mi(p)*ni(p)*(mi(r)*ni(q)+mi(q)*ni(r))               &
                   -(mi(q)*mi(r)*ni(p)*(ni(p)-1)                        &
                    +ni(q)*ni(r)*mi(p)*(mi(p)-1)))                      &
                   *fint(msum-2,nsum-2)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     T1
 173  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = (mi(p)*ni(r)-mi(r)*ni(p))              &
                   *fint(msum-1,nsum-1)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     T2
 174  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = (mi(p)*ni(r)-mi(r)*ni(p))              &
                   *(co*mi(q)*fint(msum-2,nsum-1)                      &
                    -sn*ni(q)*fint(msum-1,nsum-2))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     T3
 175  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = (mi(p)*ni(r)-mi(r)*ni(p))              &
                   *(co*ni(q)*fint(msum-1,nsum-2)                      &
                    +sn*mi(q)*fint(msum-2,nsum-1))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!     T4
 176  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = (mi(p)*ni(q)-mi(q)*ni(p))              &
                   *(co*mi(r)*fint(msum-2,nsum-1)                      &
                    -sn*ni(r)*fint(msum-1,nsum-2))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     T5
 177  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = (mi(p)*ni(q)-mi(q)*ni(p))              &
                   *(co*ni(r)*fint(msum-1,nsum-2)                      &
                    +sn*mi(r)*fint(msum-2,nsum-1))
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     W2
 178  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = co*ni(q)*fint(msum,nsum-1)             &
                               +sn*mi(q)*fint(msum-1,nsum)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
!     W3
 179  continue
      do p=1,20
        do q=1,20
          do r=1,20
            do s=1,20
              msum = mi(p)+mi(q)+mi(r)+mi(s)
              nsum = ni(p)+ni(q)+ni(r)+ni(s)
              tensor(p,q,r,s) = co*mi(q)*fint(msum-1,nsum)             &
                               -sn*ni(q)*fint(msum,nsum-1)
            enddo
          enddo
        enddo
      enddo
      go to 200
!
!
 200  continue

      ! Transform from local coordinates to global coordinates
!
        tterm = 0.
!      general case...NO SYMMETRY
      do i=1,18
        do j=1,18
            do k=1,18
              do l=1,18
                sum2 = 0.
                do q = 1,20
                  gq = gtri(q,j,iodd)
                  do r = 1,20
                    gr = gtri(r,k,iodd)
                    do p = 1,20
                      do s = 1,20
                        sum2 = sum2 + tensor(p,q,r,s)*gq*gr*             &
                                gtri(p,i,iodd)*gtri(s,l,iodd)               
                      enddo
                    enddo
                  enddo
                enddo
                tterm(i,j,k,l) = tterm(i,j,k,l) + sum2
                enddo
              enddo
            enddo
          enddo

          ! open scratch disks for use in storing metric elements
          if(iodd.eq.1) then

             ! open file in submit directory for reuse of metric terms
             if(myrank.eq.0) then
                open(itype+200,file=filename(itype), form='unformatted',     &
                     status='unknown')
                rewind(itype+200)
             endif

             open(itype+100,form='unformatted',status='scratch')
             rewind(itype+100)
          endif

          ! save tensor to disk
          do i=1,18
             do j=1,18

                ! write a record for this iodd,i,j value
                write(itype+100) iodd,i,j,((tterm(i,j,k,l),l=1,18),k=1,18)
                if(myrank.eq.0) then
                   write(itype+200) iodd,i,j,((tterm(i,j,k,l),l=1,18),k=1,18)
                endif

             enddo
          enddo
          go to 400
399       iread(itype) = 0

400    enddo ! on itype

    enddo ! on iodd

    deallocate(tensor, tterm)
    if(myrank.eq.0) then
       do itype=1,ntensor
          if(iread(itype).ne.0) then
             close(UNIT=itype+200,STATUS='KEEP')
          endif
       enddo
    endif

    return
  end subroutine opdef2

end module opdef_mod
