cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccc Array utilities cccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module modtest 
      contains
      subroutine linspace(xmin,xmax,x,n)
      real xmin,xmax
      real, dimension(n)::x
      integer  i,n
      
      do i=1,n
       x(i) = (xmax-xmin) * real(i-1) / real(n-1) + xmin
      end do
      end subroutine linspace
cccccccccccccccccccccccccccccccccccccccccccc
c$$$      subroutine linspaceint(xmin,xmax,x,n)
c$$$      integer xmin,xmax
c$$$      integer, dimension(n)::x
c$$$      integer  i,n
c$$$      
c$$$      do i=1,n
c$$$       x(i) = int((xmax-xmin) * real(i-1) / real(n-1) + xmin)
c$$$      end do
c$$$      end subroutine linspaceint
ccccccccccccccccccccccccccccccccccccccccccccc
      subroutine logspace(xmin,xmax,x,n)
      integer n
      real xmin,xmax
      real, dimension(n)::x
      call linspace(log10(xmin),log10(xmax),x,n)
      x = 10.**x
      end subroutine logspace
cccccccccccccccccccccccccccccccccccccccccccc
      subroutine myplot(x,y,n,name)
      real, dimension(n)::x,y
      integer n
      character(*) name

      call ncarcgm(1,name)

      call mapg(minval(x),maxval(x),minval(y),maxval(y),.1,1.,.1,1.)
      call trace(x,y,n,-1,-1,0.,0.)

      call frame(0)

      call plote

      end subroutine
cccccccccccccccccccccccccccccccccccccccccccc
      subroutine mysemilogyplot(x,y,n,name)
      real, dimension(n)::x,y
      integer n
      character(*) name

      call ncarcgm(1,name)

      call mapgsl(minval(x),maxval(x),minval(y),maxval(y),.1,1.,.1,1.)
      call trace(x,y,n,-1,-1,0.,0.)

      call frame(0)

      call plote

      end subroutine
cccccccccccccccccccccccccccccccccccccccccccc
      subroutine mywrite(x,name)
      real x
      character(*) name
      write(*,*) name,x
      end subroutine




cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccsubroutine equi  cccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine equi(size,r,p,dp,q,dq,d2q,Bz,Bp,k,i_system,q0,q_edge,beta,i_solver,r0,rs,m,n,is)
      implicit none
      integer size
      integer i,i_system,i_solver,is,m,n
      integer, dimension(1)::is_array
      INTEGER IER, PGBEG,ifail
      real k,q0,q_edge,alpha_q,p0,a,Bz0,step,er,Bz_edge,beta,rs,integ,Bz_edgebis
      real, dimension(size):: r,p,q,dp,dq,integral1,integral2,Bz,Bp,dBz,jz,dBp,fx,integral1bis,Bzbis,integral2bis,Ir,dq_an,d2q
      real Bp0,r0,inc1,inc2
      real qa, betaq, Lq
      real j0
      real qmin,alpha_rev,rmin
      
c$$$  dimension r(size),p(size),q(size),dp(size),dq(size),
c$$$  1integral1(size),integral2(size),Bz(size),Bp(size)
      
      a=1.
      step=a/real(size-1)
      call linspace(0.,1.,r,size)
c     r = (/((i*step),i=0,size-1)/)

cccccccTokamak ccccccccc 


      if(i_system==0)then



cccccccccccccccccccccccccccccc
cccc  Bp-solver Furth ccccccccccccc
ccccccccccccccccccccccccccccc
         if (i_solver==1) then
            is_array=minloc(abs(r-rs))
            is=is_array(1)
            write(*,*) 'is =',is
            
c     Bz_edge=1.
c     r0=1./2.
            Bz_edge=1.
            p0=0.5*beta
            p=p0*(1-(r/a)**2)
            dp=p0*(-2.)*r/a**2
cccc  Peaked model
            fx=(r/r0)/(1+(r/r0)**2)

            

            Ir=(-1./2.)*((1.+(r/r0)**2)**(-2)-(1.+(1./r0)**2)**(-2))
            Bp0=sqrt((Bz_edge**2 -2.*p(is))/((real(m)/real(n)/k*fx(is)/rs )**2+2.*Ir(is)))
            write(*,*) 'Bp0 =',Bp0
            Bp=Bp0*(r/r0)/(1+(r/r0)**2)
            dBp=((1./r0)/(1.+(r/r0)**2)+(r/r0)*(-1)/(1.+(r/r0)**2)**2*2.*r/r0**2)*Bp0
             jz=2.*r0*(r0**2/(r0**2+r**2)**2)*Bp0
c$$$            integral1=0.
c$$$            integral1bis=0.
c$$$  do 33 i=4,size
c$$$  call D01GAF(r(size-i+1:size),-dp-Bp*jz,
c$$$  1           i,integral1(size-i+1),er,ifail)
c$$$  c     write(*,*) 'integral1 =', integral1(i)
c$$$  c     write(*,*)'er =',er
c$$$  33      continue
c$$$            integ=0.
c$$$            do 34 i=1,size-1
c$$$               integ=integ+((-dp(size-i)-Bp(size-i)*jz(size-i))+(-dp(size-i+1)-Bp(size-i+1)*jz(size-i+1)))*step/2.
c$$$               integral1bis(size-i)=integ
c$$$ 34         continue
            
c$$$  call ncarcgm(1,'integs.cgm')
c$$$  
c$$$  call mapg(minval(r),maxval(r),minval(integral1),maxval(integral1),.1,1.,.1,1.)
c$$$  c      call agsetf('X/LOGARITHMIC.',1.)
c$$$  c      ta=2.
c$$$  call colora('green')
c$$$  call trace(r,integral1bis,size,-1,-1,0.,0.)
c$$$  call colora('red')
c$$$  call trace(r,integral1,size,-1,-1,0.,0.)
c$$$  call frame(0)
c$$$  c.....finalize plot file
c$$$  call plote
c$$$  
c$$$            integral1=integral1bis
            
c            Bz_edge=sqrt( (real(m)/k*Bp0*fx(is)/r(is))**2   + 2*integral1bis(is)        )
c     Bz_edgebis=sqrt( (real(m)/k*Bp0*fx(is)/r(is))**2   + 2*integral1bis(is)        )
            
            write(*,*) 'Bzedge=',Bz_edge
c     stop
            Bz=sqrt(Bz_edge**2-2.*p-2.*Bp0**2*Ir )
c     Bzbis=sqrt(Bz_edgebis**2-2*integral1bis)
            q=r*Bz*k/Bp
            q(1)=q(2)
            dq=Bz*k/Bp+k*r/Bp*(-Bp*jz-dp)/Bz-r*k*Bz/Bp**2*dBp
            dq(1)=dq(2)

            

            call mywrite(q(1),'q0=')
            call mywrite(q(size),'q_edge=')
c            call myplot(r,q,size,'q__Furth.cgm')
         end if

ccccccccccccccccccccccccccccccc
ccccccc q-solver FTZ cccccccccccccc
ccccccccccccccccccccccccccccccc  
         if (i_solver==0) then

           

            
            p0=0.5*beta
            p=p0*(1-(r/a)**2)
            dp=p0*(-2.)*r/a**2
            
            qa=3.5
            betaq=-0.3
            Lq=0.1
            
            q=qa/(2.-r**2)+betaq*(r-rs)*exp(-(r-rs)**2/Lq**2)
            dq=qa*2.*r/(2.-r**2)**2 + betaq * exp(-(r-rs)**2/Lq**2) * (1. - 2.* (r-rs)**2/Lq**2 )

c            call myplot(r(2:size),-dp(2:size)*(q(2:size)/dq(2:size))**2/r(2:size),size-1,'Ds.cgm')

      call ncarcgm(1,'Ds.cgm')
c.... make plot
      call mapg(minval(r),maxval(r),-0.01,0.5,.1,1.,.1,1.)
        call trace(r,-2*dp*(q/dq)**2/r,size,-1,-1,0.,0.)
 

      call frame(0)
c.....finalize plot file
      call plote

            
            
c$$$            integral1=0.
c$$$            integral1bis=0.
c$$$            integral2=0.
c$$$            

 
            ifail=0
c$$$            do 10 i=4,size
c$$$c     call D01GAF(r(1:i),(2.*r(1:i)-r(1:i)**2*dq(1:i)/q(1:i))/(q(1:i)**2
c$$$c     1/k**2+r(1:i)**2 ),i,integral1(i),er,ifail)
c$$$               call D01GAF(r(size-i+1:size),-(2.*r(size-i+1:size)-r(size-i+1:size
c$$$     1              )**2*dq(size-i+1:size)/q(size-i+1:size))/(q(size-i+1:size)**2
c$$$     2              /k**2+r(size-i+1:size)**2 ),i,integral1(size-i+1),er,ifail)
c$$$c     write(*,*) 'integral1 =', integral1(i)
c$$$c     write(*,*)'er =',er
c$$$ 10         continue
           
            integ=0
            do 35 i=1,size-1
               inc1= - (2.*r(size-i+1)-r(size-i+1
     1              )**2*dq(size-i+1)/q(size-i+1))/(q(size-i+1)**2
     2              /k**2+r(size-i+1)**2) 
               inc2= - (2.*r(size-i)-r(size-i
     1              )**2*dq(size-i)/q(size-i))/(q(size-i)**2
     2              /k**2+r(size-i)**2) 

               integ=integ + (inc1+inc2)/2.*step
               integral1(size-i)=integ
 35         continue

 
            
c$$$            do 11 i=4,size
c$$$c     call D01GAF(r(1:i),exp(2.*integral1(1:i))*(-(1/k)**2*q(1:i)**2*
c$$$c     1dp(1:i))/0.5/(q(1:i)**2/k**2+r(1:i)**2),i,integral2(i),er,ifail)
c$$$               call D01GAF(r(size-i+1:size),-exp(2.*integral1(size-i+1:size))*(-
c$$$     1              (1/k)**2*q(size-i+1:size)**2*
c$$$     2              dp(size-i+1:size))/0.5/(q(size-i+1:size)**2/k**2+r(size-i+1:size)
c$$$     3              **2),i,integral2(size-i+1),er,ifail)
c$$$c     write(*,*) 'integral2 =', integral2(i)
c$$$c     write(*,*)'ifail =',ifail
c$$$ 11         continue
            integ=0.
            do 36 i=1,size-1
               inc1=-exp(2.*integral1(size-i+1))*(-
     1              (1/k)**2*q(size-i+1)**2*
     2              dp(size-i+1))/0.5/(q(size-i+1)**2/k**2+r(size-i+1)
     3              **2) 
                inc2=-exp(2.*integral1(size-i))*(-
     1              (1/k)**2*q(size-i)**2*
     2              dp(size-i))/0.5/(q(size-i)**2/k**2+r(size-i)
     3              **2) 


               integ=integ + (inc1+inc2)/2.*step
               integral2bis(size-i)=integ
 36         continue


            Bz_edge=1.
            
c     Bz=exp(-integral1)*sqrt(Bz0**2+integral2)
            Bz=exp(-integral1)*sqrt(Bz_edge**2+integral2bis)
            Bp=Bz*r*k/q

            dBz=(-dp-Bz**2*k**2/q*(2.*r/q-r**2/q**2*dq ))/(Bz*(1+r**2*k**2/q**2 ))
            jz=k/r*(2.*r*Bz/q+r**2/q*dBz-r**2*Bz/q**2*dq )

c     do 12 i=1,size
c     write(*,*) 'q=', q(i)
c     write(*,*) 'Bz2 =', Bz2(i)
c     12   continue

         end if

ccccccccccccccccccccccccccccccccccccccccccccccccc
cccc  j=1-r**2 Solver (FTZ Tearing) ccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccc
         if (i_solver==2) then
            is_array=minloc(abs(r-rs))
            is=is_array(1)
            write(*,*) 'is =',is
            
c     Bz_edge=1.
c     r0=1./2.
            Bz_edge=1.
            p=0.
            dp=0.

            j0 = Bz_edge* 4. * k / 2.8
            write(*,*) 'j0 =',j0

            Bp = j0 * (r/2. - r**3/4.)
            dBp = j0 * (1./2. - 3. * r**2 / 4.)

            Bz = sqrt(Bz_edge**2 - 2. * j0**2 * (r**2/4. - 3. * r**4/16. + r**6/24. - 5./48.))
            dBz = - j0**2 / Bz * (r/2. - r**3/4.)*(1.-r**2)

            q=Bz*k/j0/(1./2.-r**2/4.)
            dq=k/j0*(dBz/(1./2.-r**2/4.)+Bz*r/2./(1./2.-r**2/4.)**2)


            
         end if

ccccccccccccccccccccccccccccccccccccccccccccccccc
cccc  j=1-r**2 Solver (FTZ Tearing) ccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccc
         if (i_solver==3) then
            is_array=minloc(abs(r-rs))
            is=is_array(1)
            write(*,*) 'is =',is

            Bz_edge=1.
            p=0.
            dp=0.

            qmin=2.
            alpha_rev=3.
            rmin=0.6
            q=qmin*((alpha_rev-1.)*(r/rmin)**4-2.*(alpha_rev-1.)*(r/rmin)**2+alpha_rev)
            dq=qmin*((alpha_rev-1.)*4.*r**3/rmin**4-2.*(alpha_rev-1.)*2.*r/rmin**2)
            d2q=qmin*((alpha_rev-1.)*12.*r**2/rmin**4-4.*(alpha_rev-1.)/rmin**2)

            ifail=0

           
            integ=0
            do 45 i=1,size-1
               inc1= - (2.*r(size-i+1)-r(size-i+1
     1              )**2*dq(size-i+1)/q(size-i+1))/(q(size-i+1)**2
     2              /k**2+r(size-i+1)**2) 
               inc2= - (2.*r(size-i)-r(size-i
     1              )**2*dq(size-i)/q(size-i))/(q(size-i)**2
     2              /k**2+r(size-i)**2) 

               integ=integ + (inc1+inc2)/2.*step
               integral1(size-i)=integ
 45         continue

 

            integ=0.
            do 46 i=1,size-1
               inc1=-exp(2.*integral1(size-i+1))*(-
     1              (1/k)**2*q(size-i+1)**2*
     2              dp(size-i+1))/0.5/(q(size-i+1)**2/k**2+r(size-i+1)
     3              **2) 
                inc2=-exp(2.*integral1(size-i))*(-
     1              (1/k)**2*q(size-i)**2*
     2              dp(size-i))/0.5/(q(size-i)**2/k**2+r(size-i)
     3              **2) 


               integ=integ + (inc1+inc2)/2.*step
               integral2bis(size-i)=integ
 46         continue


     
            
c     Bz=exp(-integral1)*sqrt(Bz0**2+integral2)
            Bz=exp(-integral1)*sqrt(Bz_edge**2+integral2bis)
            Bp=Bz*r*k/q

            dBz=(-dp-Bz**2*k**2/q*(2.*r/q-r**2/q**2*dq ))/(Bz*(1+r**2*k**2/q**2 ))
            jz=k/r*(2.*r*Bz/q+r**2/q*dBz-r**2*Bz/q**2*dq )

            

            
         end if


         write(*,*) 'Plot Bz'
         call myplot(r,Bz,size,'Bz.cgm')
         write(*,*) 'Plot dBz'
         call myplot(r,dBz,size,'dBz.cgm')
         write(*,*) 'Plot Bp'
         call myplot(r,Bp,size,'Bp.cgm')
         write(*,*) 'Plot p'
         call myplot(r,p,size,'p.cgm')
         write(*,*) 'Plot q'
         call myplot(r,q,size,'q.cgm')
         write(*,*) 'Plot dq'
         call myplot(r,dq,size,'dq.cgm')
         write(*,*) 'Plot jz'
         call myplot(r,jz,size,'jz.cgm')
         write(*,*) 'Plot jzBp'
         call myplot(r,-jz*Bp,size,'jzBp.cgm')
         

c$$$      open(unit=10,file='jz.txt')
c$$$      do 62 i=1,size
c$$$         write(10,*) jz(i)
c$$$ 62   continue
c$$$
c$$$      open(unit=11,file='r.txt')
c$$$      do 63 i=1,size
c$$$         write(11,*) r(i)
c$$$ 63   continue

      end if



      
      end subroutine



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccsubroutine main_tok_rev ccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine main_tok_rev(size,size_plot,size_plot2) 
      implicit none
      integer size,m,n,order,size_plot,size_plot2,i,j,ii
      integer istartleft,iendleft,iendright
      real step
      real d_right,d_left
      integer ifail,ic,i_solver,is
      
      real epsd,rmatch
      real E,F,G,H,K0
      real, dimension(size):: p,q,dp,dq,Bz,Bp,Bz2,Bp2,alpha_DL_array,d2q
c      real, dimension(size):: p2,q2,dp2,dq2,Bz2,Bp2
      real, dimension(size):: r
      real, dimension(size_plot)::gamma_dless_array,
     1     gamma_array, gamma_dless_GGJ_array,gamma_dless_array2, gamma_array2, gamma_dless_GGJ_array2
      real, dimension(1)::ta_array
      real, dimension(size_plot)::eta_array,det_array,eta_array2
c      real, dimension(400)::d_right_array,d_left_array,i_array,gamma_array2
      real k
      real gamma_dless
      real Bpc,dqc,rc
      real eta
      real rho
      real Lr,gammar,DI,Lr_DI,DIbis
      real ta
      real prec,prec_left
      real Bp0
      complex de,do
      real beta, q0,q_edge
      real rs,xs 
      real r0
      real, dimension(10,size_plot2)::dp_array,prec_array,prec_left_array,dr_array,dl_array
      real, dimension(10)::Cr_array
c     real, dimension(10)::iend_array
      real gamma, gamlim
      real pi,tol,T1,T2

      real, dimension(10)::tol_rf_array,epsd_array,rmatch_array
      real, dimension(size_plot)::Bp0_array,xs_array,Lr_array,gammar_array,Lr_array2,gammar_array2
      real, dimension(size_plot,size_plot2)::dp2_array
      real, dimension(size_plot2)::beta_array
      integer, dimension(size_plot2)::iend_array
      
      real,dimension(2)::prec_min_array,prec_max_array

      real, dimension(7)::coeff
      real, dimension(2,6)::sol
      real, dimension(14)::work
      real a,b,gamfun14, gamfun34, Dr, d_prime
      real erl1,ers1
      real psi1l,psi1s,eps1l,eps1s,delta0l,delta0s
      real sigma_l,sigma_s
      
      real, dimension(size_plot):: accl, accs
      integer solmark
      real Qmin, Qmax
      real alpha_DL

      integer, dimension(1)::maxalphaloc_array
      integer maxalphaloc
      common /delta_block/ d_right,d_left,E,F,G,H,K0,Lr_DI,Lr,epsd,rmatch
      rho=1.

      k=1./10.
c$$$      q0=1.470898932131725
c$$$      q_edge=7.350522933203288
      call linspace(0.9,1.5,xs_array,size_plot)
c      call logspace(0.0,6.e-3,beta_array,size_plot2)
c$$$      iend_array(1)=1
c$$$      iend_array(2)=2
c$$$      iend_array(3)=3
c      beta=0.001

c$$$      open(unit=4,file='dprime_beta2_1e3.txt')
c$$$      open(unit=5,file='beta2_1e3.txt')

     
      do 11 j=1,1
c     beta=0.5
         order=3
         step=1./(size-1.)
         r0=0.5
c     i_solver=0 : q-solver
c     i_solver=1 : Bp-solver
         i_solver=3
c     call linspace(0.03,0.06,Bp0_array,size_plot)

c$$$      beta=beta_array(j)
         beta=0.05

         do 10 i=1,1
            
c     Bp0=Bp0_array(i)
c            Bp0=0.035
            m=2
            n=1
            xs=xs_array(i)
            rs=r0*xs
            if (i_solver==0) rs=0.5
            
            call mywrite(beta,'beta= ')
c            call mywrite(Bp0,'Bp0=')
            call equi(size,r,p,dp,q,dq,d2q,Bz,Bp,k,0,q0,q_edge,beta,i_solver,r0,rs,m,n,is)
c$$$            i_solver=1

ccccccccc Comparison equilibrium with Furth Bp solver and q solver  ccccccccccccc          
c$$$            i_solver=0
c$$$            call equi(size,r,p,dp,q,dq,Bz,Bp,k,0,q0,q_edge,beta,i_solver,r0,rs,m,n,is)
c$$$            
c$$$      call ncarcgm(1,'Bzcomp.cgm')
c$$$c.... make plot
c$$$      call mapg(minval(r),maxval(r),0.99,1.01,.1,1.,.1,1.)
c$$$        call trace(r,Bz,size,-1,-1,0.,0.)
c$$$         call colora('red')
c$$$         call trace(r,Bz2,size,-1,-1,0.,0.)
c$$$
c$$$      call frame(0)
c$$$c.....finalize plot file
c$$$      call plote
c$$$
c$$$      call ncarcgm(1,'Bpcomp.cgm')
c$$$c.... make plot
c$$$      call mapg(minval(r),maxval(r),minval(Bp),0.025,.1,1.,.1,1.)
c$$$        call trace(r,Bp,size,-1,-1,0.,0.)
c$$$         call colora('red')
c$$$         call trace(r,Bp2,size,-1,-1,0.,0.)
c$$$
c$$$      call frame(0)
c$$$c.....finalize plot file
c$$$      call plote
c$$$            stop
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            write(*,*) '----------------------'
            write(*,*) 'm= ',m,'; n= ',n
            
            ifail=0
            call rational_surf(m,n,size,r,p,dp,q,dq,Bz,Bp,k,E,F,G,H,K0,ifail,
     1           dqc,Bpc,rc,ic)
            if (ifail==1) then
               
               write(*,*) 'Rational surface is outside'
               stop

            end if   



            

            
c     xs_array(i)=rc/r0
c$$$            DI=E+F+H-0.25
c$$$c            write(*,*) 'lambda =',-2.*sqrt(-DI)+2*order-2
c$$$            if (-DI.ge.0.) then
c$$$               sigma_l=-0.5-sqrt(-DI)
c$$$               sigma_s=-0.5+sqrt(-DI)
c$$$
c$$$            else
c$$$               write(*,*)                                                    'complex sigma; Suydam criterion not satisfied'
c$$$               ifail=1
c$$$               stop
c$$$            end if

            ifail=0
c$$$            istartleft=(1000001-1)/1000
            istartleft=10
            tol=1.e-12

cccccccccccc without iend study cccccccccccccccc
            iendleft=1
            iendright=1

            call outer(m,n,size,r,p,dp,q,dq,d2q,Bz,Bp,k,order,d_right,d_left,DI
     1           ,ifail,ic,iendleft,iendright,istartleft,prec,prec_left,tol,erl1,ers1)
            stop
c 05/17
c$$$            d_right = -397.9688373840144
c$$$            d_left = 396.6712643462110

c old
c$$$            d_right = -571.7752398772411
c$$$            d_left = 577.3055066092531

c beta=0.05        
            d_right = -84.59992171676191
            d_left = 89.46808711672807

            d_prime=d_right+d_left
            call mywrite(q(1),'q0=')
            call mywrite(q(size),'q_edge=')
            write(*,*) 'd_prime=',d_prime
            write(*,*) '----------------------------------------'
            
c$$$            write(4,*) d_prime
c$$$            write(5,*) beta
            
            


cccccccccccccc iend study cccccccccccccccccccc
c$$$            open(unit=4,file='iend.txt')
c$$$            open(unit=5,file='d_prime.txt')
c$$$            open(unit=7,file='prec_right.txt')
c$$$            open(unit=8,file='prec_left.txt')
c$$$            
c$$$            do 13 ii=1,5
c$$$            
c$$$            iendleft=ii*200
c$$$            iendright=ii*200
c$$$            write(4,*) ii*200
c$$$            
c$$$            call outer(m,n,size,r,p,dp,q,dq,Bz,Bp,k,order,d_right,d_left,DI
c$$$     1           ,ifail,ic,iendleft,iendright,istartleft,prec,prec_left,tol)
c$$$            dp2_array(i,j)=d_right+d_left
c$$$            write(*,*) 'd_prime=',d_right+d_left
c$$$            d_prime=d_right+d_left
c$$$            write(5,*) d_prime
c$$$            write(7,*) prec
c$$$            write(8,*) prec_left
c$$$
c$$$ 13         continue

cccccccccccccc tol study cccccccccccccccccccc
c$$$            open(unit=4,file='tol_5.txt')
c$$$            open(unit=5,file='d_prime_tol_5.txt')
c$$$            open(unit=7,file='time_5.txt')
c$$$c            open(unit=7,file='prec_right.txt')
c$$$c            open(unit=8,file='prec_left.txt')
c$$$            
c$$$            do 13 ii=1,4
c$$$               if (ii==1) tol=1.e-7
c$$$               if (ii==2) tol=1.e-8
c$$$               if (ii==3) tol=1.e-9
c$$$               if (ii==4) tol=1.e-10
c$$$c               if (ii==5) tol=1.e-11
c$$$c               if (ii==6) tol=1.e-12              
c$$$            write(4,*) tol
c$$$            CALL CPU_TIME(T1)
c$$$            call outer(m,n,size,r,p,dp,q,dq,Bz,Bp,k,order,d_right,d_left,DI
c$$$     1           ,ifail,ic,iendleft,iendright,istartleft,prec,prec_left,tol)
c$$$            CALL CPU_TIME(T2)
c$$$            write(*,*) 'tol =',tol
c$$$            write(*,*) 'time =',T2-T1
c$$$            write(7,*) 'time =',T2-T1
c$$$            dp2_array(i,j)=d_right+d_left
c$$$            write(*,*) 'd_prime=',d_right+d_left
c$$$            d_prime=d_right+d_left
c$$$            write(5,*) d_prime
c$$$c            write(7,*) prec
c$$$c            write(8,*) prec_left
c$$$
c$$$ 13         continue

c           stop
 10      continue
 11   continue
      
c$$$      call ncarcgm(1,'fig1.cgm')
c$$$c.... make plot
c$$$      call mapg(minval(xs_array),maxval(xs_array),minval(dp2_array*r0),maxval(dp2_array*r0),.1,1.,.1,1.)
c$$$c     call agsetf('X/LOGARITHMIC.',1.)
c$$$c     ta=2.
c$$$      
c$$$      do 13 j=1,1
c$$$         if (j==1) call colora('green')
c$$$         if (j==2) call colora('yellow')
c$$$         if (j==3) call colora('red')
c$$$         call trace(xs_array,dp2_array(:,j)*r0,size_plot,-1,-1,0.,0.)
c$$$ 13   continue
c$$$      call frame(0)
c$$$c.....finalize plot file
c$$$      call plote

cccccccccccccccccccccccccccccccccccccccccccccccc
cccccccc Find gamma code cccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccc

c$$$      eta=1.e-5
c$$$
c$$$
c$$$      gammar=(eta*real(n)**2*Bpc**2*dqc**2/(rho*rc**2))**(1./3.)
c$$$
c$$$      Lr=(rho*eta**2*rc**2/(real(n)**2*Bpc**2*dqc**2))**(1./6.)
c$$$      Lr_DI=Lr**(-2.*sqrt(-DI))
c$$$
c$$$
c$$$      call linspace(1.e-3,0.9,gamma_dless_array,size_plot)
c$$$      do 99 i=1,size_plot
c$$$         det_array(i)=det(gamma_dless_array(i))
c$$$ 99   continue
c$$$c      call myplot(gamma_dless_array,det_array/maxval(det_array),size_plot,'det.cgm')
c$$$         
c$$$         call ncarcgm(1,'det.cgm')
c$$$
c$$$         call mapg(minval(gamma_dless_array),maxval(gamma_dless_array),-1.,1. ,.1,1.,.1,1.)
c$$$
c$$$
c$$$         call trace(gamma_dless_array, det_array/maxval(abs(det_array)),size_plot,-1,-1,0.,0.)
c$$$
c$$$
c$$$         call frame(0)
c$$$
c$$$         call plote
c$$$
c$$$         open(unit=1,file='det.txt')
c$$$         open(unit=2, file='Q.txt')
c$$$         
c$$$         do 33 i=1,size_plot
c$$$            write(1,*) det_array(i)/maxval(abs(det_array))
c$$$            write(2,*) gamma_dless_array(i)
c$$$ 33         continue
c$$$
c$$$
c$$$         stop

cccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccc tol of rootfinder study ccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$
c$$$      eta=1.e-6
c$$$     
c$$$       gammar=(eta*real(n)**2*Bpc**2*dqc**2/(rho*rc**2))**(1./3.)
c$$$       gammar_array2(i)=gammar
c$$$      Lr=(rho*eta**2*rc**2/(real(n)**2*Bpc**2*dqc**2))**(1./6.)
c$$$      Lr_DI=Lr**(-2.*sqrt(-DI))
c$$$
c$$$      write(*,*) '-----------------------------------'
c$$$      write(*,*) 'eta =',eta
c$$$
c$$$      gamlim=0.36
c$$$      call logspace(1.e-10,3.806487718378616E-16,tol_rf_array,10)
c$$$
c$$$      open(unit=4,file='tol_rf.txt')
c$$$      open(unit=5,file='gamma_rf.txt')
c$$$
c$$$      do 59 i=1,10
c$$$               write(*,*) '------------tol_rf =',tol_rf_array(i)
c$$$      call C05ADF(5.e-3,gamlim,tol_rf_array(i),0.,det,gamma_dless,ifail)  
c$$$      gamma=gammar*gamma_dless
c$$$      write(4,*) tol_rf_array(i)
c$$$
c$$$      write(5,*) gamma
c$$$ 59   continue
c$$$
c$$$      stop
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccc epsd study ccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$
c$$$      eta=1.e-6
c$$$     
c$$$       gammar=(eta*real(n)**2*Bpc**2*dqc**2/(rho*rc**2))**(1./3.)
c$$$       gammar_array2(i)=gammar
c$$$      Lr=(rho*eta**2*rc**2/(real(n)**2*Bpc**2*dqc**2))**(1./6.)
c$$$      Lr_DI=Lr**(-2.*sqrt(-DI))
c$$$
c$$$      write(*,*) '-----------------------------------'
c$$$      write(*,*) 'eta =',eta
c$$$
c$$$      gamlim=0.36
c$$$      call logspace(1.e-3, 2.782559402207126E-05,epsd_array,10)
c$$$
c$$$      open(unit=4,file='epsd.txt')
c$$$      open(unit=5,file='gamma_epsd.txt')
c$$$
c$$$      do 59 i=1,10
c$$$               write(*,*) '------------epsd =',epsd_array(i)
c$$$               epsd=epsd_array(i)
c$$$      call C05ADF(5.e-3,gamlim,1.e-15,0.,det,gamma_dless,ifail)  
c$$$      gamma=gammar*gamma_dless
c$$$      write(4,*) epsd_array(i)
c$$$
c$$$      write(5,*) gamma
c$$$ 59   continue
c$$$
c$$$      stop
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccc rmatch study ccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$
c$$$      eta=1.e-6
c$$$     
c$$$       gammar=(eta*real(n)**2*Bpc**2*dqc**2/(rho*rc**2))**(1./3.)
c$$$       gammar_array2(i)=gammar
c$$$      Lr=(rho*eta**2*rc**2/(real(n)**2*Bpc**2*dqc**2))**(1./6.)
c$$$      Lr_DI=Lr**(-2.*sqrt(-DI))
c$$$
c$$$      write(*,*) '-----------------------------------'
c$$$      write(*,*) 'eta =',eta
c$$$
c$$$      gamlim=0.36
c$$$      call linspace(3., 13.,rmatch_array,10)
c$$$
c$$$      open(unit=4,file='rmatch.txt')
c$$$      open(unit=5,file='gamma_rmatch.txt')
c$$$
c$$$      do 59 i=1,10
c$$$               write(*,*) '------------rmatch =',rmatch_array(i)
c$$$               rmatch=rmatch_array(i)
c$$$      call C05ADF(5.e-3,gamlim,1.e-15,0.,det,gamma_dless,ifail)  
c$$$      gamma=gammar*gamma_dless
c$$$      write(4,*) rmatch_array(i)
c$$$
c$$$      write(5,*) gamma
c$$$ 59   continue
c$$$
c$$$      stop
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc


c       call logspace(0.5e-8,1.e-3,eta_array,size_plot)






c$$$       call logspace(1.e-12,3.e-6,eta_array2,size_plot)
c$$$
c$$$
c$$$c      call linspace(0.101,0.9,gamma_dless_array,size_plot)
c$$$c      eta=1.e-6
c$$$
c$$$
c$$$      do 30 i=1,size_plot
c$$$
c$$$      ifail=0
c$$$ 
c$$$
c$$$
c$$$c$$$      eta=eta_array(i)
c$$$c$$$     
c$$$c$$$       gammar=(eta*real(n)**2*Bpc**2*dqc**2/(rho*rc**2))**(1./3.)
c$$$c$$$       gammar_array(i)=gammar
c$$$c$$$      Lr=(rho*eta**2*rc**2/(real(n)**2*Bpc**2*dqc**2))**(1./6.)
c$$$c$$$      Lr_DI=Lr**(-2.*sqrt(-DI))
c$$$c$$$      Lr_array(i)=Lr
c$$$c$$$c      write(*,*) 'DI =',DI
c$$$c$$$c       write(*,*) 'DIbis =',DIbis
c$$$c$$$     
c$$$c$$$    
c$$$c$$$c      det_array(i) = det(gamma_dless_array(i))
c$$$c$$$      
c$$$c$$$   
c$$$c$$$
c$$$c$$$
c$$$c$$$
c$$$c$$$
c$$$c$$$
c$$$c$$$
c$$$c$$$      write(*,*) '-----------------------------------'
c$$$c$$$      write(*,*) 'eta =',eta
c$$$c$$$c      write(*,*) 'alpha=',alpha
c$$$c$$$c      call C05AJF(gamma_dless,0.0001,0.0,det,500,ifail)
c$$$c$$$c      write(*,*) '--- biss ---'
c$$$c$$$c$$$      if (j==1) gamlim=0.2
c$$$c$$$c$$$      if (j==2) gamlim=0.07
c$$$c$$$c$$$      if (j==3) gamlim=0.03
c$$$c$$$c$$$      if (j==4) gamlim=0.02
c$$$c$$$c$$$      if (j==5) gamlim=0.015
c$$$c$$$
c$$$c$$$
c$$$c$$$ccccc for beta = 0.01
c$$$c$$$c      gamlim=0.35
c$$$c$$$ccccc for beta= 0.1
c$$$c$$$c      gamlim=0.75
c$$$c$$$cccc for beta=0.001 and eta<=1.e-8
c$$$c$$$c$$$      gamlim=0.09
c$$$c$$$c$$$      call C05ADF(1.e-3,gamlim,1.e-15,0.,det,gamma_dless,ifail)
c$$$c$$$ccccc for beta=0.001
c$$$c$$$
c$$$c$$$      call C05ADF(0.101,0.6,1.e-15,0.,det,gamma_dless,ifail)
c$$$c$$$
c$$$c$$$      gamma=gammar*gamma_dless
c$$$c$$$      write(*,*) 'gamma dimensionless =',gamma_dless
c$$$c$$$      write(*,*) 'gamma =',gamma
c$$$c$$$
c$$$c$$$      gamma_dless_array(i)=gamma_dless
c$$$c$$$      gamma_array(i)=gamma
c$$$
c$$$      eta=eta_array2(i)
c$$$     
c$$$      gammar=(eta*real(n)**2*Bpc**2*dqc**2/(rho*rc**2))**(1./3.)
c$$$      gammar_array2(i)=gammar
c$$$      Lr=(rho*eta**2*rc**2/(real(n)**2*Bpc**2*dqc**2))**(1./6.)
c$$$      Lr_DI=Lr**(-2.*sqrt(-DI))
c$$$      Lr_array2(i)=Lr
c$$$
c$$$
c$$$
c$$$
c$$$c      Lr_array(i)=Lr
c$$$c      write(*,*) 'DI =',DI
c$$$c       write(*,*) 'DIbis =',DIbis
c$$$     
c$$$    
c$$$c      det_array(i) = det(gamma_dless_array(i))
c$$$      
c$$$   
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$      write(*,*) '-----------------------------------'
c$$$      write(*,*) 'eta =',eta
c$$$
c$$$c      gamlim=0.14
c$$$c$$$      Qmin=3.e-3
c$$$c$$$      Qmax=4.8e-3
c$$$      Qmin=0.01
c$$$      Qmax=0.3
c$$$      call C05ADF(Qmin,Qmax,1.e-15,0.,det,gamma_dless,ifail)  
c$$$      gamma=gammar*gamma_dless
c$$$      gamma_dless_array2(i)=gamma_dless
c$$$      gamma_array2(i)=gamma
c$$$      write(*,*) 'gamma dimensionless =',gamma_dless
c$$$      write(*,*) 'gamma =',gamma
c$$$

c$$$
c$$$cccccccccc Accuracy AM ccccccccccccccccc
c$$$      psi1l=(0.5*(1.-2.*sigma_l))*(-gamma_dless**2*sigma_l*(sigma_l-1)-delta0l*(E+F-H*(sigma_l-2)) )
c$$$      psi1s=(0.5*(1.-2.*sigma_s))*(-gamma_dless**2*sigma_s*(sigma_s-1)-delta0s*(E+F-H*(sigma_s-2)) )
c$$$
c$$$      delta0l=-gamma_dless**2*(G+K*F)+gamma_dless**2*(G-K*E)-gamma_dless**2*K*H*(sigma_l+1)
c$$$      delta0s=-gamma_dless**2*(G+K*F)+gamma_dless**2*(G-K*E)-gamma_dless**2*K*H*(sigma_s+1)
c$$$
c$$$      eps1l=psi1l+(H+gamma_dless)*(sigma_l+1)+(E+F)/gamma_dless
c$$$      eps1s=psi1s+(H+gamma_dless)*(sigma_s+1)+(E+F)/gamma_dless
c$$$      
c$$$      accl(i)=Lr*abs(eps1l)**(0.5)*erl1
c$$$      accs(i)=Lr*abs(eps1s)**(0.5)*ers1
c$$$
c$$$ccccccccccccccccccccccccccccccccccccc
c$$$ 30   continue
c$$$
c$$$ 
c$$$c$$$      call myplot(gamma_dless_array,det_array/maxval(abs(det_array)),size_plot                 ,'det_beta-3.cgm')
c$$$
c$$$
c$$$c            write(*,*) '************************************'
c$$$c      call C05AJF(gamma_dless,1.e-6,0.0,det,500,ifail)
c$$$c      write(*,*) '--- biss ---'
c$$$c      call C05ADF(1.e-5,0.06,1.e-6,0.,det,gamma_dless_bis,ifail)
c$$$c      gamma=gammar*gamma_dless
c$$$c      gamma_bis=gammar*gamma_dless_bis
c$$$c      write(*,*) 'gamma dimensionless =',gamma_dless
c$$$c      write(*,*) 'gamma =',gamma
c$$$c      write(*,*) 'gammabis dimensionless =',gamma_dless_bis
c$$$c      write(*,*) 'gammabis =',gamma_bis
c$$$ 
c$$$
c$$$
c$$$c      write(*,*) 'eta =',eta
c$$$c      write(*,*) 'log(resistivity) =',log(eta*ta_array(j))
c$$$c      write(*,*) 'ln(gamma) =',log(gamma)
c$$$      write(*,*) 'd_right =',d_right
c$$$      write(*,*) 'd_left =',d_left
c$$$            write(*,*) 'DI =',DI
c$$$c       write(*,*) 'DIbis =',DIbis
c$$$c       write(*,*) 'ta =',ta_array(j)
c$$$       write(*,*) 'order =',order
c$$$       write(*,*) 'lambda =',-2.*sqrt(-DI)+2*order-2
c$$$       write(*,*) 'Dr =',E+F
c$$$       pi=3.14159265358979323846264338327950288419716939937510
c$$$       gamfun34=1.2254167024651776451290983033628905268512392481080706112301
c$$$       gamfun14=3.6256099082219083119306851558676720029951676828800654674333
c$$$       Dr=E+F
c$$$c$$$      do 50 i=1,size_plot
c$$$c$$$          a=2*pi**2*rc**2*(1./k)/Lr_array(i)*gamfun34/gamfun14
c$$$c$$$          b=a*pi*Dr/4.
c$$$c$$$          coeff(1)=a
c$$$c$$$          coeff(6)=-d_prime
c$$$c$$$          coeff(7)=-b
c$$$c$$$          solmark=0
c$$$c$$$          call C02AGF(coeff,6,.true.,sol,work,ifail)
c$$$c$$$          do 51 j=1,6
c$$$c$$$             if(sol(2,j)==0.) then
c$$$c$$$                if(sol(1,j).ge.0.) then
c$$$c$$$                   gamma_dless_GGJ_array(i)=sol(1,j)**4
c$$$c$$$                   if (solmark==0) then
c$$$c$$$                      solmark=1
c$$$c$$$                   else
c$$$c$$$                      write(*,*) 'multiple real positive growth rate'
c$$$c$$$                      stop
c$$$c$$$                   end if
c$$$c$$$                end if
c$$$c$$$             end if
c$$$c$$$ 51          continue
c$$$c$$$          
c$$$c$$$         call mywrite(gamma_array(i)/(gamma_dless_GGJ_array(i)*gammar_array(i)),'mygamma/GGJgamma=')
c$$$c$$$ 50      continue
c$$$
c$$$      do 52 i=1,size_plot
c$$$          a=2*pi**2*rc**2*(1./k)/Lr_array2(i)*gamfun34/gamfun14
c$$$          b=a*pi*Dr/4.
c$$$          coeff(1)=a
c$$$          coeff(6)=-d_prime
c$$$          coeff(7)=-b
c$$$          solmark=0
c$$$          call C02AGF(coeff,6,.true.,sol,work,ifail)
c$$$          do 53 j=1,6
c$$$             if(sol(2,j)==0.) then
c$$$                if(sol(1,j).ge.0.) then
c$$$                   gamma_dless_GGJ_array2(i)=sol(1,j)**4
c$$$                   if (solmark==0) then
c$$$                      solmark=1
c$$$                   else
c$$$                      write(*,*) 'multiple real positive growth rate'
c$$$                      stop
c$$$                   end if
c$$$                end if
c$$$             end if
c$$$ 53       continue
c$$$          
c$$$c         call mywrite(gamma_array(i)/(gamma_dless_GGJ_array(i)*gammar_array(i)),'mygamma/GGJgamma=')
c$$$ 52   continue
c$$$      call myplot(eta_array2,accl,size_plot,'accl.cgm')
c$$$      call myplot(eta_array2,accs,size_plot,'accs.cgm')
c$$$c$$$      open(unit=1,file='eta_final.txt')
c$$$c$$$      do 55 i=1,size_plot
c$$$c$$$         write(1,*) eta_array2(i)
c$$$c$$$ 55   continue
c$$$c$$$
c$$$c$$$      open(unit=2,file='gamma_final.txt')
c$$$c$$$      do 56 i=1,size_plot
c$$$c$$$         write(2,*) gamma_array2(i)
c$$$c$$$ 56   continue
c$$$c$$$
c$$$c$$$      open(unit=3,file='gamma_GGJ_final.txt')
c$$$c$$$      do 57 i=1,size_plot
c$$$c$$$         write(3,*) gamma_dless_GGJ_array2(i)*gammar_array2(i)
c$$$c$$$ 57   continue
c$$$
c$$$c.....initialize plotting
c$$$      call ncarcgm(1,'acc.cgm')
c$$$c....make plot
c$$$      call mapgll(minval(eta_array2),maxval(eta_array2),minval(accl),maxval(accl),.1,1.,.1,1.)
c$$$
c$$$      call trace(eta_array2,accl,size_plot,-1,-1,0.,0.)
c$$$
c$$$      call colora('green')
c$$$
c$$$      call trace(eta_array2,accs,size_plot,-1,-1,0.,0.)
c$$$
c$$$      call frame(0)
c$$$c.....finalize plot file
c$$$      call plote
c$$$
c$$$c.....initialize plotting
c$$$      call ncarcgm(1,'gamma_eta.cgm')
c$$$c....make plot
c$$$      call mapgll(minval(eta_array2),maxval(eta_array2),minval(gamma_array2),maxval(gamma_dless_GGJ_array2*gammar_array2),.1,1.,.1,1.)
c$$$c      call agsetf('X/LOGARITHMIC.',1.)
c$$$c      do 98 i=1,10
c$$$c      ta=1.5+i*0.1
c$$$c      ta_array(1)=2.3
c$$$      do 97 j=1,1
c$$$c         if(j==2) call colora('green')
c$$$c         if(j==3) call colora('red')  
c$$$
c$$$c      call trace(eta_array,gamma_array,size_plot,-1,-1,0.,0.)
c$$$      call trace(eta_array2,gamma_array2,size_plot,-1,-1,0.,0.)
c$$$cccccccccc equation 90 of GGJ (approx of eq 88 of GGJ at small eta)
c$$$      call colora('green')
c$$$c      call trace(eta_array,((pi*(E+F)/4.)**(2./3.)*(eta_array*real(n)**2*Bpc**2*dqc**2/(rho*rc**2))**(1./3.)) ,size_plot,-1,-1,0.,0.)
c$$$      call trace(eta_array2,((pi*(E+F)/4.)**(2./3.)*(eta_array2*real(n)**2*Bpc**2*dqc**2/(rho*rc**2))**(1./3.)) ,size_plot,-1,-1,0.,0.)
c$$$cccccccccc equation 93 of GGJ
c$$$      if(d_prime.gt.0.) then
c$$$      call colora('blue')
c$$$c      call trace(eta_array,(d_prime/(2*pi**2*rc**2*(1./k)/Lr_array*gamfun34/gamfun14))**(4./5.)
c$$$c     A*(eta_array*real(n)**2*Bpc**2*dqc**2/(rho*rc**2))**(1./3.) ,size_plot,-1,-1,0.,0.)
c$$$      call trace(eta_array2,(d_prime/(2*pi**2*rc**2*(1./k)/Lr_array2*gamfun34/gamfun14))**(4./5.)
c$$$     A *(eta_array2*real(n)**2*Bpc**2*dqc**2/(rho*rc**2))**(1./3.) ,size_plot,-1,-1,0.,0.)
c$$$      endif
c$$$
c$$$cccccccccc equation 88 of GGJ
c$$$      call colora('red')
c$$$c      call trace(eta_array,(gamma_dless_GGJ_array*gammar_array) ,size_plot,-1,-1,0.,0.)      
c$$$      call trace(eta_array2,(gamma_dless_GGJ_array2*gammar_array2) ,size_plot,-1,-1,0.,0.)   
c$$$ 97   continue
c$$$c 98   continue 
c$$$
c$$$      call tginit(1)
c$$$      call setpch(0,2,-1,1)
c$$$      call colora('white')
c$$$      call pointc('*',(/1.e-6,2.5e-6,5.e-6,7.5e-6,1.e-5,2.5e-5,5.e-5,7.5e-5,         1.e-4,2.5e-4/),
c$$$     A (/4.65e-4,8.38e-4,12.18e-4,15.0075e-4,16.56e-4,23.335e-4,27.94e-4,31.83e-4,29.91e-4,25.71e-4 /),10,1,1,0.,0.)
c$$$
c$$$      call frame(0)
c$$$c.....finalize plot file
c$$$      call plote
c$$$
c$$$
c$$$c$$$      open(unit=1,file='eta_C2.txt')
c$$$c$$$c$$$      do 61 i=1,size_plot
c$$$c$$$c$$$         write(1,*) eta_array(i)
c$$$c$$$c$$$ 61   continue
c$$$c$$$      do 62 i=1,size_plot
c$$$c$$$         write(1,*) eta_array2(i)
c$$$c$$$ 62   continue  
c$$$c$$$
c$$$c$$$      open(unit=3,file='gamma_C2.txt')
c$$$c$$$c$$$      do 65 i=1,size_plot
c$$$c$$$c$$$         write(3,*) gamma_array(i)
c$$$c$$$c$$$ 65   continue
c$$$c$$$      do 66 i=1,size_plot
c$$$c$$$         write(3,*) gamma_array2(i)
c$$$c$$$ 66   continue
c$$$c$$$    
c$$$c$$$      open(unit=2,file='gamma_inter_C2.txt')
c$$$c$$$c$$$      do 63 i=1,size_plot
c$$$c$$$c$$$         write(2,*) ((pi*(E+F)/4.)**(2./3.)*(eta_array(i)*real(n)**2*Bpc**2*dqc**2/(rho*rc**2))**(1./3.))
c$$$c$$$c$$$ 63   continue
c$$$c$$$      do 64 i=1,size_plot
c$$$c$$$         write(2,*) ((pi*(E+F)/4.)**(2./3.)*(eta_array2(i)*real(n)**2*Bpc**2*dqc**2/(rho*rc**2))**(1./3.))
c$$$c$$$ 64   continue
c$$$c$$$
c$$$c$$$      open(unit=4,file='gamma_tear_C2.txt')
c$$$c$$$c$$$      do 67 i=1,size_plot
c$$$c$$$c$$$         write(4,*) (d_prime/(2*pi**2*rc**2*(1./k)/Lr_array(i)*gamfun34/gamfun14))**(4./5.)
c$$$c$$$c$$$     A*(eta_array(i)*real(n)**2*Bpc**2*dqc**2/(rho*rc**2))**(1./3.)
c$$$c$$$c$$$ 67   continue
c$$$c$$$      do 68 i=1,size_plot
c$$$c$$$         write(4,*) (d_prime/(2*pi**2*rc**2*(1./k)/Lr_array2(i)*gamfun34/gamfun14))**(4./5.)
c$$$c$$$     A*(eta_array2(i)*real(n)**2*Bpc**2*dqc**2/(rho*rc**2))**(1./3.)
c$$$c$$$ 68   continue
c$$$c$$$
c$$$c$$$      open(unit=5,file='gamma_an_C2.txt')
c$$$c$$$c$$$      do 69 i=1,size_plot
c$$$c$$$c$$$         write(5,*) gamma_dless_GGJ_array(i)*gammar_array(i)
c$$$c$$$c$$$ 69   continue
c$$$c$$$      do 70 i=1,size_plot
c$$$c$$$         write(5,*) gamma_dless_GGJ_array2(i)*gammar_array2(i)
c$$$c$$$ 70      continue

      end subroutine



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccc subroutine rational_surf cccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rational_surf(m,n,size,r,p,dp,q,dq,Bz,Bp,k,E,F,G,H,K0
     1     ,ifail,dqc,Bpc,rc,ic)
      implicit none
      integer m,n,size,ic,ifail
      real, dimension(size):: r,p,q,dp,dq,Bz,Bp
      real k,qc,rc,dpc,dqc,Bzc,Bpc,pc,E,F,G,H,K0,ta
      integer, dimension(1)::ic_array,val_array
      real tab
      qc=real(m)/real(n)
      

c     val_array=minval(abs(q-qc))
c     write(*,*) 'min =',val_array(1)
c     write(*,*) 'min =',minval(abs(q-qc))
c     write(*,*) 'qc =',qc
      if (minval(abs(q-qc)).ge.0.01)  then
         ifail=1
c     write(*,*) 'outside'
         return
      end if
      ic_array=minloc(abs(q-qc))
      ic=ic_array(1)
      write(*,*) 'ic =',ic
      rc=r(ic)
      write(*,*) 'rc =',rc
      Bpc=Bp(ic)
      Bzc=Bz(ic)


      write(*,*) 'Bc =',sqrt(Bzc**2+Bpc**2)
      write(*,*) 'Bzc =',Bzc
      write(*,*) 'Bpc =',Bpc

      pc=p(ic)
      dpc=dp(ic)
      dqc=dq(ic)
      F=dpc**2*rc**2/(Bpc**4/k**2*dqc**2)
      E = -2.*dpc*rc/(Bpc**2/k**2*dqc**2) - F
      K0=1./F
      if(pc.ne.0) G = 3./5. * (Bzc**2+Bpc**2)/pc
      
c     G = 3./5. * Bpc**2/pc
      H=0.
      

      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccouterfunc cccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine outerfunc(rin,vecin,vecout)
      implicit none
      real rin,rc
      real, dimension(2)::vecin,vecout
      real, dimension(1000001)::rb,fr,gr,rbb,frb,grb,rbb2
      real, dimension(ninterp*(ninterp+1)/2)::Cf,Cg
      real, dimension(1)::Cfedge,Cgedge
      integer, dimension(1)::icb
      real, dimension(2)::tempx,tempy
      integer ic,last,n,np,ninterp
      common /fg_block/ fr,gr,rc,ic,rb
      common /size_interp/ ninterp
      n=ninterp
      frb=fr
      grb=gr
      rbb=rb
      rbb2=rb
c     write(*,*) 'r(icbfirst)',rbb(2)
      icb=minloc(abs(rbb-rin))
      
c      write(*,*) '----------------------------'
c$$$
c      write(*,*) 'rin=',rin
c$$$      write(*,*) 'icbfirst=',icb(1)
c$$$      write(*,*) 'r(icbfirst)',rbb(icb(1))
c      write(*,*) 'er=',vecin(1)
c      write(*,*) 'ksir=',vecin(2)
c     icb=minloc(abs(rbb-rin))
      if(rbb(icb(1)).lt.rin .or. rin==1.)then 
         icb(1)=icb(1)-1
      end if
c$$$      write(*,*) 'icb=',icb(1)
c     .or. icb(1)==1000001-2 .or. icb(1)==1000001-3  .or. icb(1)==3 .or. icb(1)==4
      if (icb(1)==1000001-1  .or. icb(1)==2 )  then
         last=1
         call E01AAF(rbb(icb(1):icb(1)+1), frb(icb(1):icb(1)+1), Cfedge, 2, 1, 1, rin)
         vecout(1)=vecin(2)/Cfedge(last)
c$$$  write(*,*) 'f=',Cfedge(last)
         call E01AAF(rbb2(icb(1):icb(1)+1), grb(icb(1):icb(1)+1), Cgedge, 2,1, 1, rin)
         vecout(2)=vecin(1)*Cgedge(last)
c     write(*,*) 'g=',Cgedge(last)

c$$$  write(*,*) 'der=',vecout(1)
c$$$  write(*,*) 'dksir=',vecout(2)
c$$$  write(*,*) 'this was the edge'

                 
c$$$         call ncarcgm(1,'interp.cgm')
c$$$
c$$$         call mapg(minval(rb((icb(1)):(icb(1)+1))),maxval(rb((icb(1)):(icb(1)+1)))+0.000001,minval(gr((icb(1)):(icb(1)+1))),
c$$$     A        maxval(gr((icb(1)):(icb(1)+1))),.1,1.,.1,1.)
c$$$
c$$$         call colora('green')
c$$$         call trace(rb((icb(1)):(icb(1)+1)), gr((icb(1)):(icb(1)+1)),2,-1,-1,0.,0.)
c$$$         call colora('red')
c$$$         tempx(1)=rin
c$$$         tempx(2)=rin
c$$$         tempy(1)=0.
c$$$         tempy(2)=Cgedge(last)
c$$$         
c$$$         call trace(tempx,tempy,2,-1,-1,0.,0.)
c$$$         call frame(0)
c$$$
c$$$         call plote
c$$$         
c$$$         stop
         
      else

         
c     write(*,*) rb(101)
c     nint must be even
         
         last=n*(n+1)/2
         np=(n+1)/2
c     last=(1000001-1)*1000001/2
         
c     write(*,*) icb(1)
         
         call E01AAF(rbb((icb(1)-(np-1)):(icb(1)+np)), frb((icb(1)-(np-1)):(icb(1)+np)), Cf,n+1, last, n, rin)
c     write(*,*) 'fr=',C(last)
         vecout(1)=vecin(2)/Cf(last)
c$$$  write(*,*) 'f=',Cf(last)
c$$$  
         call E01AAF(rbb2((icb(1)-(np-1)):(icb(1)+np)), grb((icb(1)-(np-1)):(icb(1)+np)), Cg, n+1, last, n, rin)
c     Cg(last)=gr(icb(1))
         vecout(2)=vecin(1)*Cg(last)





c$$$  write(*,*) 'g=',Cg(last)
c$$$  write(*,*) 'der=',vecout(1)
c$$$  write(*,*) 'dksir=',vecout(2)
      end if
      end subroutine

      subroutine outputRKfuncleft(xsol,y)
      real xsol, interval, zleft, d_left_shoot
      real, dimension(2)::y
      integer n
      real,dimension(3)::er_l,er_s
      real, dimension(3)::ksi_l,ksi_s
      real ksi_l_right, ksi_s_right, ksi_l_left,ksi_s_left
      real er_l_right, er_s_right, er_l_left,er_s_left
      real der_l_right, der_s_right, der_l_left,der_s_left
      real matchleft
      real rc
      integer ic
      integer i0p
      real c1, c2, f2, f3
      real, dimension(1000001):: fr,gr,rb
      common /ksi_er_block/ ksi_l, ksi_s, er_l, er_s, sigma_s, sigma_l, i0p, c1, c2, f2, f3
      common /fg_block/ fr,gr,rc,ic,rb
      interval=1.e-6
            if (xsol.le.rc-0.002) then
         xsol=rc-0.001
         return
         endif

      write(11,*) xsol
      write(23,*) y(1)
      write(*,*) 'r = ',xsol

     
      zleft=y(1)/(y(2)/funfr(xsol))

      matchleft=rc-xsol

      if (i0p .eq. 0) then

      er_s_left=(matchleft)**(sigma_s)
      er_l_left=(matchleft)**(sigma_l)



      do 17 i=1,3


         er_l_left=er_l_left+er_l(i)*(-1.)**i*(matchleft)**(sigma_l+i)
         er_s_left=er_s_left+er_s(i)*(-1.)**i*(matchleft)**(sigma_s+i)

 17   continue

      ksi_l_left=0.
      ksi_s_left=0.
      do 16 i=1,3

         ksi_l_left=ksi_l_left+ksi_l(i)*(-matchleft)**i
         ksi_s_left=ksi_s_left+ksi_s(i)*(-matchleft)**i
 16   continue

      ksi_l_left=(matchleft)**sigma_l*ksi_l_left
      ksi_s_left=(matchleft)**sigma_s*ksi_s_left
      

      der_s_left=ksi_s_left/funfr(xsol)
      der_l_left=ksi_l_left/funfr(xsol)

      endif


ccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccc beta=0 cccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (i0p .eq. 1) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccc  er_l = er_s ; ksi_l = ksi_s => We always use the _s ones
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccc Large solution cccccccccccccccccccccccccc
         er_l_left = 1./matchleft * (1. + er_s(1) * (-matchleft) + ( c2 + er_s(2))*matchleft**2 + c1 * (-matchleft) * log(matchleft ) * (1. +  er_s(1) * (-matchleft) +  er_s(2)*matchleft**2))

         ksi_l_left = - f2 / matchleft  * (-matchleft - 2.*er_s(1)*matchleft**2 - ( 3.*er_s(2)-er_s(1)**2 + f3/f2 * er_s(1)) * (-matchleft)**3 - c1 * (-matchleft)**3 * log (matchleft) * (er_s(1) + (f3/f2*er_s(1)+2.*er_s(2))*(-matchleft)))

ccccccccc Small solution ccccccccccccccccccccccccccc
      er_s_left=(matchleft)**(sigma_s)

      write(*,*) 'er_s_left_zero_order =',(matchleft)**(sigma_s)

      do 18 i=1,2
         write(*,*) 'order =',i

         write(*,*) 'er_s_left =',er_s(i)*(-1.)**i*(matchleft)**(sigma_s+i)

         write(*,*) 'er_s_left_conv =',er_s(i)*(-1.)**i*(matchleft)**(sigma_s+i)/er_s_left

         er_s_left=er_s_left+er_s(i)*(-1.)**i*(matchleft)**(sigma_s+i)

 18      continue

      ksi_s_left=0.
      do 19 i=1,3
         write(*,*) 'i=',i

 
         write(*,*) 'ksi_s_left',ksi_s(i)*(-matchleft)**i
         if(i.ge.3) then


            write(*,*) 'ksi_s_left_conv',ksi_s(i)*(-matchleft)**i/ksi_s_left
         end if

  
         ksi_s_left=ksi_s_left+ksi_s(i)*(-matchleft)**i
 19      continue

      ksi_s_left=(matchleft)**sigma_s*ksi_s_left
 
      der_s_left=ksi_s_left/funfr(xsol)
      der_l_left=ksi_l_left/funfr(xsol)
 
         endif

      d_left=(zleft*der_l_left-er_l_left)/(er_s_left-zleft*der_s_left)

      write(10,*) d_left

      xsol=xsol+interval


      end subroutine

      subroutine outputRKfuncright(xsol,y)
      real xsol, interval, zleft, d_left_shoot
      real, dimension(2)::y
      integer n
      real,dimension(3)::er_l,er_s
      real, dimension(3)::ksi_l,ksi_s
      real ksi_l_right, ksi_s_right, ksi_l_left,ksi_s_left
      real er_l_right, er_s_right, er_l_left,er_s_left
      real der_l_right, der_s_right, der_l_left,der_s_left
      real matchright
      real rc
      integer ic
      real, dimension(1000001):: fr,gr,rb
      integer i0p
      real c1,c2,f2,f3
      common /ksi_er_block/ ksi_l, ksi_s, er_l, er_s, sigma_s, sigma_l, i0p, c1, c2, f2, f3
      common /fg_block/ fr,gr,rc,ic,rb
      interval=1.e-6
      if (xsol.eq.1.) then
         xsol=rc+0.001
         return
         endif

      write(13,*) xsol
      write(22,*) y(1)
      write(*,*) 'r = ',xsol
      
     
      zright=y(1)/(y(2)/funfr(xsol))

      matchright=xsol-rc

      if (i0p .eq. 0) then

      er_s_right=(matchright)**(sigma_s)
      er_l_right=(matchright)**(sigma_l)      


c$$$
c$$$      write(*,*) 'er_l_right_zero_order =',(matchright)**(sigma_s)
c$$$      write(*,*) 'er_s_right_zero_order =',(matchright)**(sigma_l)     
c$$$


      do 17 i=1,3
c$$$         write(*,*) 'order =',i
c$$$c     fact=fact*i
c$$$
c$$$
c$$$         write(*,*) 'er_l_right =',er_l(i)*(matchright)**(sigma_l+i)
c$$$         write(*,*) 'er_s_right =',er_s(i)*(matchright)**(sigma_s+i)
c$$$
c$$$
c$$$         write(*,*) 'er_l_right_conv =',er_l(i)*(matchright)**(sigma_l+i)/er_l_right
c$$$         write(*,*) 'er_s_right_conv =',er_s(i)*(matchright)**(sigma_s+i)/ er_s_right


         er_l_right=er_l_right+er_l(i)*(matchright)**(sigma_l+i)
         er_s_right=er_s_right+er_s(i)*(matchright)**(sigma_s+i)



 17   continue

      ksi_l_right=0.
      ksi_s_right=0.
      ksi_l_left=0.
      ksi_s_left=0.
      do 16 i=1,3
c$$$         write(*,*) 'i=',i
c$$$         write(*,*) 'ksi_l_right',ksi_l(i)*(matchright)**i
c$$$         write(*,*) 'ksi_s_right',ksi_s(i)*(matchright)**i
c$$$
c$$$         if(i.ne.1) then
c$$$            write(*,*) 'ksi_l_right_conv',ksi_l(i)*(matchright)**i/ksi_l_right
c$$$            write(*,*) 'ksi_s_right_conv',ksi_s(i)*(matchright)**i/ ksi_s_right
c$$$
c$$$         end if
         ksi_l_right=ksi_l_right+ksi_l(i)*(matchright)**i
         ksi_s_right=ksi_s_right+ksi_s(i)*(matchright)**i

 16   continue



      ksi_l_right=(matchright)**sigma_l*ksi_l_right
      ksi_s_right=(matchright)**sigma_s*ksi_s_right

      


      der_s_right=ksi_s_right/funfr(xsol)
      der_l_right=ksi_l_right/funfr(xsol)







      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccc beta=0 cccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (i0p .eq. 1) then
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccc  er_l = er_s ; ksi_l = ksi_s => We always use the _s ones
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccc Large solution cccccccccccccccccccccccccc
         er_l_right = 1./matchright * (1. + er_s(1) * (matchright) + ( c2 + er_s(2))*matchright**2 + c1 * (matchright) * log(matchright ) * (1. +  er_s(1) * (matchright) +  er_s(2)*matchright**2))

         ksi_l_right = - f2 / matchright  * (matchright - 2.*er_s(1)*matchright**2 - ( 3.*er_s(2)-er_s(1)**2 + f3/f2 * er_s(1)) * (matchright)**3 - c1 * (matchright)**3 * log (matchright) * (er_s(1) + (f3/f2*er_s(1)+2.*er_s(2))*(matchright)))

ccccccccc Small solution ccccccccccccccccccccccccccc
      er_s_right=(matchright)**(sigma_s)    

      write(*,*) 'er_s_right_zero_order =',(matchright)**(sigma_s)

      do 18 i=1,2
         write(*,*) 'order =',i
c     fact=fact*i

         write(*,*) 'er_s_right =',er_s(i)*(-1.)**i*(matchright)**(sigma_s+i)
         write(*,*) 'er_s_right_conv =',er_s(i)*(-1.)**i*(matchright)**(sigma_s+i)/er_s_right
    
         er_s_right=er_s_right+er_s(i)*(-1.)**i*(matchright)**(sigma_s+i)


 18   continue

      ksi_s_right=0.
      do 19 i=1,3
         write(*,*) 'i=',i

      
         write(*,*) 'ksi_s_right',ksi_s(i)*(matchright)**i
         if(i.ge.3) then

            write(*,*) 'ksi_s_right_conv',ksi_s(i)*(matchright)**i/ksi_s_right
         end if
  
         ksi_s_right=ksi_s_right+ksi_s(i)*(matchright)**i
 19   continue
  
      ksi_s_right=(matchright)**sigma_s*ksi_s_right
      
      der_s_right=ksi_s_right/funfr(xsol)
      der_l_right=ksi_l_right/funfr(xsol)
         endif

      d_right=(zright*der_l_right-er_l_right)/(er_s_right-zright*der_s_right)
      write(12,*) d_right

      xsol=xsol-interval




c$$$      write(*,*) '-----------------------------'
c$$$      write(*,*) 'r=',xsol
c$$$      write(*,*) 'er=',y(1)
      end subroutine



c$$$  
c$$$  function gfunc(x,y)
c$$$  real x
c$$$  real, dimension(2)::y
c$$$  c      real rc
c$$$  c      common /fg_block/ rc
c$$$  c      gfunc=x-rc
c$$$  gfunc=y(1)-1.186127220166149E-04
c$$$  end function




cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccsubroutine outer cccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine outer(m,n,size,r,p,dp,q,dq,d2q,Bz,Bp,k,order,
     1     d_right,d_left,DI,ifail,icb,iendleft,iendright,istartleft,prec,precleft,tol,erl1,ers1)
      implicit none
      real D02EJW,D02EJX,D02EJY
      external D02EJW,D02EJX,D02EJY
      integer size,m,n,order,ic,ifail,fact,icb
      real, dimension(1000001):: r,p,q,dp,dq,Bz,Bp,fr,gr,rb,grbis,d2q,Bzb
      real, dimension(ic-1):: er_l_left_array, er_s_left_array
      real, dimension(1000001-ic)::er_l_right_array,er_s_right_array
      real, dimension(1000001)::ksir,er,ksil,el
      real k,qc,rc,step,sigma_l,sigma_s,d_left,d_right
      real d_prime

      real DI
      real mr,nr
c     integer, dimension(1)::ic_array
      real, dimension(14)::dfrc,erest,dgrc,errors_dfrc,errors_dgrc,Bzrc,errors_Bzrc
      real,dimension(3)::er_l,er_s
      real, dimension(3)::ksi_l,ksi_s
      real eps1_l,eps2_l,eps1_s,eps2_s
      real eps,er0,el0,der0,del0,ksir0,ksil0
      real z
      real ksi_l_right, ksi_s_right, ksi_l_left,ksi_s_left
      real er_l_right, er_s_right, er_l_left,er_s_left
      real der_l_right, der_s_right, der_l_left,der_s_left
      real ksi1_l, ksi2_l, ksi3_l, ksi1_s, ksi2_s, ksi3_s
      real alpha,dqc,ta
      integer istartleft,iendright,iendleft,i,j
      integer iw
      real zright,zleft
      real, dimension(2)::vecin,vecout,vec_r0,vec_l0
      real, dimension(40)::W
      real, dimension(14*2+50)::W2
      real tol,rightstart,rightend,leftstart,leftend
c      real der_s_left_bis,der_l_left_bis,der_s_right_bis,der_l_right_bis
c      real er_s_left_bis,er_l_left_bis,er_s_right_bis,er_l_right_bis
      real Cr
      integer ninterp
      real d_left_test
      real matchright, matchleft
      real prec,precleft
      real erl1,ers1
      integer i0p
      real c1,c2,f2,f3
      real f4,g2
c     real c1,c2,c3,c4
      common /fg_block/ fr,gr,rc,ic,rb,Bzb
 
c     ,rc,ic,rb
      common /size_interp/ ninterp
      common /ksi_er_block/ ksi_l, ksi_s, er_l, er_s, sigma_s, sigma_l, i0p, c1, c2, f2, f3
      ninterp=3
c     common /plot_block/ c1,c2,c3,c4
      step=1./real(size-1)
c     supposes a=1
      qc=real(m)/real(n)
      ic=icb
      rb=r
      mr=real(m)
      nr=real(n)

      Bzb=Bz
      rc=r(icb)
      dqc=dq(icb)
      i0p=0
      if (dp(icb).eq.0.) i0p=1

c$$$  do 10 i=1,size
c$$$  write(*,*) '------------'
c$$$  write(*,*) 'r=',r(i)
c$$$  write(*,*) 'q=',q(i)
c$$$  write(*,*) 'dq=',dq(i)
c$$$  c         write(*,*) 'fr',fr(i) 
c$$$  c         write(*,*) 'gr=',gr(i) 
c$$$  c         if (i==50) stop
c$$$  c         if(abs(x-(rc-(2*i-1)*step)).le.step*0.01)funfr=fr(ic-(2*i-1))
c$$$  10   continue
c     write(*,*) q
      
      fr = r * Bp**2 * (mr-nr*q)**2 / (mr**2+nr**2*k**2*r**2)
c$$$  gr = 2.*nr**2*k**2*r**2*dp/(mr**2+nr**2*k**2*r**2) + 1./r * Bp**2 * 
c$$$  1(mr-nr*q)**2*(nr**2*k**2*r**2+mr**2-1)/(mr**2+nr**2*k**2*r**2)+
c$$$  22.* nr**2*k**2*r * Bp**2 * (nr**2*q**2-mr**2)/ 
c$$$  3 (mr**2+nr**2*k**2*r**2)**2

c$$$  gr=nr**2*k**2*r*Bp**2/(mr**2+nr**2*k**2*r**2)*(-alpha/4.*(dq/k**2)**2+(1.+(mr**2-1.)/(nr**2*k**2*r**2))*(mr-nr*q)**2+2.*(nr**2*q**2-mr**2)/(mr**2+nr**2*k**2*r**2) )
c$$$  gr(1)=gr(2)
c$$$  gr(size)=gr(size-1)
      gr=2.*nr**2*k**2*r**2/(k**2*nr**2*r**2+mr**2)*dp + Bp**2 / r *(mr-nr*q)**2*(k**2*nr**2*r**2+m**2-1.)/(k**2*r**2*nr**2+mr**2)
     A     + 2.*k**2*nr**2*r*Bp**2*(n**2*q**2-m**2)/(k**2*nr**2*r**2+m**2)**2
c$$$  
c$$$  gr=grbis  
      gr(1)=gr(2)
      gr(size)=gr(size-1)+gr(size-1)-gr(size-2)
c$$$  write(*,*) 'fr(2000)=',fr(2000)
c$$$  call outerfunc((r(2000)+r(4000))/2.,vecin,vecout)
c$$$  write(*,*) 'fr(4000)=',fr(4000)

c.....initialize plotting
      call ncarcgm(1,'fr.cgm')
c.... make plot
      call mapg(0.,1.,minval(fr),maxval(fr),.1,1.,.1,1.)
c     call agsetf('X/LOGARITHMIC.',1.) 
      ta=2.
      call trace(r,fr,size
     1     ,-1,-1,0.,0.)
      call frame(0)
c.....finalize plot file
      call plote



c.....initialize plotting
      call ncarcgm(1,'gr.cgm')
c.... make plot
      call mapg(0.,1.,minval(gr),maxval(gr),.1,1.,.1,1.)
c     call agsetf('X/LOGARITHMIC.',1.)
      ta=2.
      call trace(r,gr,size,-1,-1,0.,0.)

      call frame(0)
c.....finalize plot file
      call plote


c$$$  write(*,*) 'test'
c$$$  do 10 i=1,size
c$$$  write(*,*) '------------'
c$$$  write(*,*) 'r=',r(i)
c$$$  write(*,*) 'fr',fr(i) 
c$$$  write(*,*) 'gr=',gr(i) 
c$$$  c         if (i==50) stop
c$$$  c         if(abs(x-(rc-(2*i-1)*step)).le.step*0.01)funfr=fr(ic-(2*i-1))
c$$$  10   continue

c     write(*,*) 'frc =', funfr(rc)

c     stop
      write(*,*) 'g0=',gr(icb)
      call D04AAF(rc,14,1.e-4,dfrc,erest,funfr,ifail)
      errors_dfrc=erest/abs(dfrc)
      call D04AAF(rc,14,1.e-4,dgrc,erest,fungr,ifail)
      errors_dgrc=erest/abs(dgrc)
      call D04AAF(rc,14,1.e-4,Bzrc,erest,funBz,ifail)
      errors_Bzrc=erest/abs(Bzrc)
      
      do 12 i=1,14
         write(*,*) '------- i=',i   
         write(*,*) 'error_dfrc =', errors_dfrc(i)
         write(*,*) 'error_dgrc =', errors_dgrc(i)
         write(*,*) 'error_Bzrc =', errors_Bzrc(i)
         write(*,*) 'dfrc =', dfrc(i)
         write(*,*) 'dgrc =', dgrc(i)
         write(*,*) 'Bzrc =', Bzrc(i)

 12   continue
      fact=1
      do 33 i=2,order+2
         fact=fact*i
         dfrc(i)=dfrc(i)/real(fact)
         dgrc(i)=dgrc(i)/real(fact)
c$$$         write(*,*) 'dfrc(',i,')=',dfrc(i)
c$$$         write(*,*) 'dgrc(',i,')=',dgrc(i)
 33   continue

      matchleft=real(step*iendleft)
      matchright=real(step*iendright)


cccccccccccccccccccccccccccccccccccccccccccccccc
cc Calculation of the exponents cccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccc
      
      f4=rc*Bp(icb)**2.*nr**2/(mr**2+nr**2*k**2*rc**2)*(d2q(icb))**2/6.
      g2=2.*nr**2*k**2*rc*Bp(icb)**2/(n**2*k**2*rc**2+mr**2)**2*nr**2*qc*d2q(icb)
      write(*,*) 'f4=',f4
      write(*,*) 'f4num=',dfrc(4)
      write(*,*) 'g2=',g2
      write(*,*) 'g2num=',dgrc(2)
      sigma_l=-3./2.-sqrt(9./4.-g2/f4)
      sigma_s=-3./2.+sqrt(9./4.-g2/f4)
      write(*,*) 'sigma_l=',sigma_l
      write(*,*) 'sigma_s=',sigma_s

      write(*,*) '2+G2/F4=',2.+g2/f4
      write(*,*) '-A12*A21b/A22 =',1./k*qc*(-1./rc**2*Bzrc(1)+1./rc*Bzrc(2))/Bp(icb)/(d2q(icb)/2.)/(nr/mr)
      sigma_l=-3./2.-sqrt(1./4.+1./k*qc*(-1./rc**2*Bzrc(1)+1./rc*Bzrc(2))/Bp(icb)/(d2q(icb)/2.)/(nr/mr))
      sigma_s=-3./2.+sqrt(1./4.+1./k*qc*(-1./rc**2*Bzrc(1)+1./rc*Bzrc(2))/Bp(icb)/(d2q(icb)/2.)/(nr/mr))
      write(*,*) 'sigma_l_FOREST=',sigma_l
      write(*,*) 'sigma_s_FOREST=',sigma_s
      
      stop
c$$$cccccccccccccccccccccccccccccccccccccccccccccccc
c$$$cccccccc     FINITE PRESSURE   ccccccccccccccc
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$
c$$$      if (i0p==0) then
c$$$
c$$$      if (-DI.ge.0.) then
c$$$         sigma_l=-0.5-sqrt(-DI)
c$$$         sigma_s=-0.5+sqrt(-DI)
c$$$
c$$$      else
c$$$         write(*,*) 'complex sigma; Suydam criterion not satisfied'
c$$$         ifail=1
c$$$         return
c$$$      end if
c$$$
c$$$c Suydam criterion must be fulfilled only at left endpoints ==> not necessary to check that      
c$$$c$$$      if(-dp(size)-r(size)/8.*Bz(size)**2*(dq(size)/q(size))**2 .ge. 0) then
c$$$c$$$         write(*,*) 'Suydam criterion not satisfied at the edge'
c$$$c$$$         stop
c$$$c$$$      end if
c$$$
c$$$         
c$$$
c$$$      write(*,*) 'sigma_l=',sigma_l
c$$$      write(*,*) 'sigma_s=',sigma_s
c$$$      write(*,*) 'sigma_s-sigma_l=',sigma_s-sigma_l
c$$$      er_l(1)=-(dfrc(3)*sigma_l*(sigma_l+2)-dgrc(1))/(dfrc(2)
c$$$     1     *(sigma_l+1)*(sigma_l+2)-gr(ic))
c$$$      er_s(1)=-(dfrc(3)*sigma_s*(sigma_s+2)-dgrc(1))/(dfrc(2)
c$$$     1     *(sigma_s+1)*(sigma_s+2)-gr(ic))
c$$$      do 13 i=2,order
c$$$         er_l(i)=(dfrc(i+2)*(sigma_l)
c$$$     1        *(sigma_l+i+1)-dgrc(i))
c$$$         er_s(i)=(dfrc(i+2)*(sigma_s)
c$$$     1        *(sigma_s+i+1)-dgrc(i))
c$$$         do 14 j=1,i-1
c$$$            er_l(i)=er_l(i)+er_l(j)*(dfrc(i+2-j)*(sigma_l+j)
c$$$     1           *(sigma_l+i+1)-dgrc(i-j))
c$$$            er_s(i)=er_s(i)+er_s(j)*(dfrc(i+2-j)*(sigma_s+j)
c$$$     1           *(sigma_s+i+1)-dgrc(i-j))
c$$$ 14      continue
c$$$         er_l(i)=-er_l(i)/(dfrc(2)*(sigma_l+i)*(sigma_l+i+1)-gr(ic))
c$$$         er_s(i)=-er_s(i)/(dfrc(2)*(sigma_s+i)*(sigma_s+i+1)-gr(ic))
c$$$ 13   continue
c$$$
c$$$c     better precision for order=2
c$$$      eps1_l=-sigma_l/2.*(sigma_l*(sigma_l+2)*dfrc(3)-dgrc(1))/gr(ic)
c$$$      eps2_l=-(sigma_l*(sigma_l+3)*dfrc(4)-dgrc(2)+eps1_l*((sigma_l+1)*(sigma_l+3)*dfrc(3)-dgrc(1)))/((sigma_l+2)*(sigma_l+3)*dfrc(2)-gr(ic))
c$$$
c$$$      eps1_s=-sigma_s/2.*(sigma_s*(sigma_s+2)*dfrc(3)-dgrc(1))/gr(ic)
c$$$      eps2_s=-(sigma_s*(sigma_s+3)*dfrc(4)-dgrc(2)+eps1_s*((sigma_s+1)*(sigma_s+3)*dfrc(3)-dgrc(1)))/((sigma_s+2)*(sigma_s+3)*dfrc(2)-gr(ic))
c$$$
c$$$      ksi1_l=gr(ic)/(sigma_l+1.)
c$$$      ksi2_l=0.5*(dgrc(1)-sigma_l**2*dfrc(3))
c$$$      ksi3_l=sigma_l*dfrc(4)+(sigma_l+1)*dfrc(3)*eps1_l+(sigma_l+2)*dfrc(2)*eps2_l
c$$$
c$$$      ksi1_s=gr(ic)/(sigma_s+1)
c$$$      ksi2_s=0.5*(dgrc(1)-sigma_s**2*dfrc(3))
c$$$      ksi3_s=sigma_s*dfrc(4)+(sigma_s+1)*dfrc(3)*eps1_s+(sigma_s+2)*dfrc(2)*eps2_s
c$$$
c$$$      write(*,*) 'ksi1_l=',ksi1_l
c$$$      write(*,*) 'bis', sigma_l*dfrc(2)
c$$$      write(*,*) 'ksi2_l=',ksi2_l
c$$$      write(*,*) 'ksi3_l=',ksi3_l
c$$$      write(*,*) 'ksi1_s=',ksi1_s
c$$$      write(*,*) 'ksi2_s=',ksi2_s
c$$$      write(*,*) 'ksi3_s=',ksi3_s
c$$$
c$$$
c$$$      er_l(1)=eps1_l
c$$$      er_s(1)=eps1_s
c$$$      er_l(2)=eps2_l
c$$$      er_s(2)=eps2_s
c$$$
c$$$      erl1=er_l(1)
c$$$      ers1=er_s(1)
c$$$      
c$$$      
c$$$      return
c$$$      
c$$$
c$$$      ksi_l(1)=ksi1_l
c$$$      ksi_l(2)=ksi2_l
c$$$      ksi_l(3)=ksi3_l
c$$$
c$$$      ksi_s(1)=ksi1_s
c$$$      ksi_s(2)=ksi2_s
c$$$      ksi_s(3)=ksi3_s
c$$$
c$$$
c$$$      
c$$$
c$$$      er_s_right=(matchright)**(sigma_s)
c$$$      er_l_right=(matchright)**(sigma_l)      
c$$$
c$$$      er_s_left=(matchleft)**(sigma_s)
c$$$      er_l_left=(matchleft)**(sigma_l)
c$$$
c$$$      write(*,*) 'er_l_right_zero_order =',(matchright)**(sigma_s)
c$$$      write(*,*) 'er_s_right_zero_order =',(matchright)**(sigma_l)     
c$$$      write(*,*) 'er_l_left_zero_order =',(matchleft)**(sigma_s)
c$$$      write(*,*) 'er_s_left_zero_order =',(matchleft)**(sigma_l)
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$      do 17 i=1,order
c$$$         write(*,*) 'order =',i
c$$$c     fact=fact*i
c$$$
c$$$
c$$$         write(*,*) 'er_l_right =',er_l(i)*(matchright)**(sigma_l+i)
c$$$         write(*,*) 'er_s_right =',er_s(i)*(matchright)**(sigma_s+i)
c$$$         write(*,*) 'er_l_left =',er_l(i)*(-1.)**i*(matchleft)**(sigma_l+i)
c$$$         write(*,*) 'er_s_left =',er_s(i)*(-1.)**i*(matchleft)**(sigma_s+i)
c$$$
c$$$         write(*,*) 'er_l_right_conv =',er_l(i)*(matchright)**(sigma_l+i)/er_l_right
c$$$         write(*,*) 'er_s_right_conv =',er_s(i)*(matchright)**(sigma_s+i)/ er_s_right
c$$$         write(*,*) 'er_l_left_conv =',er_l(i)*(-1.)**i*(matchleft)**(sigma_l+i)/er_l_left
c$$$         write(*,*) 'er_s_left_conv =',er_s(i)*(-1.)**i*(matchleft)**(sigma_s+i)/er_s_left
c$$$
c$$$
c$$$
c$$$
c$$$         er_l_right=er_l_right+er_l(i)*(matchright)**(sigma_l+i)
c$$$         er_s_right=er_s_right+er_s(i)*(matchright)**(sigma_s+i)
c$$$         er_l_left=er_l_left+er_l(i)*(-1.)**i*(matchleft)**(sigma_l+i)
c$$$         er_s_left=er_s_left+er_s(i)*(-1.)**i*(matchleft)**(sigma_s+i)
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$ 17   continue
c$$$
c$$$      ksi_l_right=0.
c$$$      ksi_s_right=0.
c$$$      ksi_l_left=0.
c$$$      ksi_s_left=0.
c$$$      do 16 i=1,3
c$$$         write(*,*) 'i=',i
c$$$         write(*,*) 'ksi_l_right',ksi_l(i)*(matchright)**i
c$$$         write(*,*) 'ksi_s_right',ksi_s(i)*(matchright)**i
c$$$         write(*,*) 'ksi_l_left',ksi_l(i)*(-matchleft)**i
c$$$         write(*,*) 'ksi_s_left',ksi_s(i)*(-matchleft)**i
c$$$         if(i.ne.1) then
c$$$            write(*,*) 'ksi_l_right_conv',ksi_l(i)*(matchright)**i/ksi_l_right
c$$$            write(*,*) 'ksi_s_right_conv',ksi_s(i)*(matchright)**i/ ksi_s_right
c$$$            write(*,*) 'ksi_l_left_conv',ksi_l(i)*(-matchleft)**i/ksi_l_left
c$$$            write(*,*) 'ksi_s_left_conv',ksi_s(i)*(-matchleft)**i/ksi_s_left
c$$$         end if
c$$$         ksi_l_right=ksi_l_right+ksi_l(i)*(matchright)**i
c$$$         ksi_s_right=ksi_s_right+ksi_s(i)*(matchright)**i
c$$$         ksi_l_left=ksi_l_left+ksi_l(i)*(-matchleft)**i
c$$$         ksi_s_left=ksi_s_left+ksi_s(i)*(-matchleft)**i
c$$$ 16   continue
c$$$
c$$$
c$$$
c$$$      ksi_l_right=(matchright)**sigma_l*ksi_l_right
c$$$      ksi_s_right=(matchright)**sigma_s*ksi_s_right
c$$$      ksi_l_left=(matchleft)**sigma_l*ksi_l_left
c$$$      ksi_s_left=(matchleft)**sigma_s*ksi_s_left
c$$$      
c$$$
c$$$
c$$$      der_s_right=ksi_s_right/fr(ic+iendright)
c$$$      der_l_right=ksi_l_right/fr(ic+iendright)
c$$$      der_s_left=ksi_s_left/fr(ic-iendleft)
c$$$      der_l_left=ksi_l_left/fr(ic-iendleft)
c$$$      endif
c$$$
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c--------------- ZERO PRESSURE ------------------
c$$$c***********************************************
c$$$c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$
c$$$      if (i0p==1) then
c$$$         sigma_s=0.
c$$$         sigma_l=-1.
c$$$      f2=dfrc(2)
c$$$      f3=dfrc(3)
c$$$      er_s(1) = dgrc(1) / 2. / dfrc(2)
c$$$      er_s(2) = 1. / 6. / dfrc(2) * (dgrc(2) - dgrc(1) / 2. / dfrc(2) * (3. *dfrc(3) - dgrc(1)))
c$$$      ksi_s(1)=0.
c$$$      ksi_s(2)=0.5*dgrc(1)
c$$$      ksi_s(3)=1./3. * (dgrc(2)+0.5*dgrc(1)**2/dfrc(2))
c$$$      c1=f3/f2 + 2. * er_s(1)
c$$$      c2=dfrc(4)/dfrc(2)+er_s(1)**2+2.*er_s(2)+2.*f3/f2*er_s(1)-c1**2
c$$$
c$$$ccccccccccccccccccccccccccccccccccccccccccc
c$$$c---------------LEFT--------------------cccc
c$$$ccccccccccccccccccccccccccccccccccccccccccc
c$$$
c$$$ccccccccc Large solution cccccccccccccccccccccccccc
c$$$         er_l_left = 1./matchleft * (1. + er_s(1) * (-matchleft) + ( c2 + er_s(2))*matchleft**2 + c1 * (-matchleft) * log(matchleft ) * (1. +  er_s(1) * (-matchleft) +  er_s(2)*matchleft**2))
c$$$
c$$$         ksi_l_left = - f2 / matchleft  * (-matchleft - 2.*er_s(1)*matchleft**2 - ( 3.*er_s(2)-er_s(1)**2 + f3/f2 * er_s(1)) * (-matchleft)**3 - c1 * (-matchleft)**3 * log (matchleft) * (er_s(1) + (f3/f2*er_s(1)+2.*er_s(2))*(-matchleft)))
c$$$
c$$$ccccccccc Small solution ccccccccccccccccccccccccccc
c$$$      er_s_left=(matchleft)**(sigma_s)
c$$$
c$$$      write(*,*) 'er_s_left_zero_order =',(matchleft)**(sigma_s)
c$$$
c$$$      do 18 i=1,2
c$$$         write(*,*) 'order =',i
c$$$c     fact=fact*i
c$$$
c$$$         write(*,*) 'er_s_left =',er_s(i)*(-1.)**i*(matchleft)**(sigma_s+i)
c$$$
c$$$         write(*,*) 'er_s_left_conv =',er_s(i)*(-1.)**i*(matchleft)**(sigma_s+i)/er_s_left
c$$$
c$$$         er_s_left=er_s_left+er_s(i)*(-1.)**i*(matchleft)**(sigma_s+i)
c$$$ 18   continue
c$$$
c$$$   
c$$$      ksi_s_left=0.
c$$$      do 19 i=1,3
c$$$         write(*,*) 'i=',i
c$$$
c$$$c$$$         write(*,*) 'ksi_l_left',ksi_l(i)*(-matchleft)**i
c$$$         write(*,*) 'ksi_s_left',ksi_s(i)*(-matchleft)**i
c$$$         if(i.ne.1) then
c$$$
c$$$c$$$            write(*,*) 'ksi_l_left_conv',ksi_l(i)*(-matchleft)**i/ksi_l_left
c$$$            write(*,*) 'ksi_s_left_conv',ksi_s(i)*(-matchleft)**i/ksi_s_left
c$$$         end if
c$$$
c$$$  
c$$$         ksi_s_left=ksi_s_left+ksi_s(i)*(-matchleft)**i
c$$$ 19   continue
c$$$
c$$$      ksi_s_left=(matchleft)**sigma_s*ksi_s_left
c$$$
c$$$      der_s_left=ksi_s_left/fr(ic-iendleft)
c$$$      der_l_left=ksi_l_left/fr(ic-iendleft)
c$$$
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c---------------RIGHT-------------------------
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$
c$$$
c$$$ccccccccc Large solution cccccccccccccccccccccccccc
c$$$         er_l_right = 1./matchright * (1. + er_s(1) * (matchright) + ( c2 + er_s(2))*matchright**2 + c1 * (matchright) * log(matchright ) * (1. +  er_s(1) * (matchright) +  er_s(2)*matchright**2))
c$$$
c$$$         ksi_l_right = - f2 / matchright  * (matchright - 2.*er_s(1)*matchright**2 - ( 3.*er_s(2)-er_s(1)**2 + f3/f2 * er_s(1)) * (matchright)**3 - c1 * (matchright)**3 * log (matchright) * (er_s(1) + (f3/f2*er_s(1)+2.*er_s(2))*(matchright)))
c$$$
c$$$ccccccccc Small solution ccccccccccccccccccccccccccc
c$$$      er_s_right=(matchright)**(sigma_s)    
c$$$
c$$$      write(*,*) 'er_s_right_zero_order =',(matchright)**(sigma_s)
c$$$
c$$$      do 20 i=1,2
c$$$         write(*,*) 'order =',i
c$$$c     fact=fact*i
c$$$
c$$$         write(*,*) 'er_s_right =',er_s(i)*(-1.)**i*(matchright)**(sigma_s+i)
c$$$         write(*,*) 'er_s_right_conv =',er_s(i)*(-1.)**i*(matchright)**(sigma_s+i)/er_s_right
c$$$    
c$$$         er_s_right=er_s_right+er_s(i)*(-1.)**i*(matchright)**(sigma_s+i)
c$$$
c$$$
c$$$ 20   continue
c$$$
c$$$      ksi_s_right=0.
c$$$      do 21 i=1,3
c$$$         write(*,*) 'i=',i
c$$$
c$$$        
c$$$         write(*,*) 'ksi_s_right',ksi_s(i)*(matchright)**i
c$$$         if(i.ne.1) then
c$$$
c$$$         
c$$$            write(*,*) 'ksi_s_right_conv',ksi_s(i)*(matchright)**i/ksi_s_right
c$$$         end if
c$$$  
c$$$         ksi_s_right=ksi_s_right+ksi_s(i)*(matchright)**i
c$$$ 21   continue
c$$$  
c$$$      ksi_s_right=(matchright)**sigma_s*ksi_s_right
c$$$      
c$$$      der_s_right=ksi_s_right/fr(ic+iendright)
c$$$      der_l_right=ksi_l_right/fr(ic+iendright)
c$$$
c$$$
c$$$      endif
c$$$
c$$$
c$$$      write(*,*) '-------------'
c$$$      write(*,*) 'iendright=',iendright
c$$$      write(*,*) 'iendleft=',iendleft
c$$$
c$$$
c$$$
c$$$
c$$$      write(*,*) 'er_l_right=',er_l_right
c$$$      write(*,*) 'er_s_right=',er_s_right
c$$$
c$$$      write(*,*) 'er_l_left=',er_l_left
c$$$      write(*,*) 'er_s_left=',er_s_left
c$$$
c$$$
c$$$      
c$$$      write(*,*) 'der_l_right=',der_l_right
c$$$      write(*,*) 'der_s_right=',der_s_right
c$$$
c$$$      write(*,*) 'der_l_left=',der_l_left
c$$$      write(*,*) 'der_s_left=',der_s_left
c$$$
c$$$      write(*,*) 'fr_right',fr(ic+iendright)
c$$$      write(*,*) 'fr_left',fr(ic-iendleft)

      
 
ccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccRUNGE KUTTA cccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccc RIGHT REGION cccccc 
      open(unit=12,file='dright.txt')
      open(unit=13, file='rright.txt')
      open(unit=22, file='eright.txt')
      ifail=0
      vec_r0(1)=0.
      vec_r0(2)=fr(size)*(-1.e-2)
      iw=14*2+50
      rightstart=1.
      rightend=rc+iendright*step
c      tol=1.e-9
      call D02EJF(rightstart,rightend,2, vec_r0,outerfunc,D02EJY,tol,    'D', outputRKfuncright,D02EJW,W2,iw,ifail)
      
      if (ifail.ne.0) then
         write(*,*) 'RK failed on right'
         stop
         end if

c     call D02BJF(1.,rc+iendright*step,2, vec_r0,outerfunc,0.01,'D', D02BJX,D02BJW,W,ifail)
c     write(*,*) 'ouf'
      write(*,*) 'zright_RK=',vec_r0(1)/(vec_r0(2)/fr(ic+iendright))
      write(*,*) 'er=',vec_r0(1)
      write(*,*) 'der=',vec_r0(2)/fr(ic+iendright)
      write(*,*) 'fr=',fr(ic+iendright)
      zright=vec_r0(1)/(vec_r0(2))*fr(ic+iendright)

c$$$ccccccc Suydam check ccccccccccccc
c$$$      ifail=0
c$$$      vec_r0(1)=er_s_right_bis
c$$$      vec_r0(2)=ksi_s_right
c$$$      rightstart=1.
c$$$      rightend=rc+iendright*step
c$$$      tol=1.e-8
c$$$      write(*,*) '----------- Suydam check right'
c$$$      call D02EJF(rightend,rightstart,2, vec_r0,outerfunc,D02EJY,tol,    'D', D02EJX,D02EJW,W2,iw,ifail)


      

c     zrightinv=
ccccc LEFT REGION cccccc
      open(unit=10,file='dleft.txt')
      open(unit=11, file='rleft.txt')
      open(unit=23, file='eleft.txt')
c$$$      do 62 i=1,size
c$$$         write(10,*) jz(i)
c$$$ 62   continue
      ifail=0
c      tol=1.e-9
      iw=14*2+50
      if(m==1)then
         vec_l0(1)=-1.
         vec_l0(2)=0.
      else
         vec_l0(1)=(istartleft*step)**(m-1)*1.e-1
         vec_l0(2)=real(m-1)*(istartleft*step)**(m-2)*fr(istartleft+1)*1.e-1
      end if
      leftstart=istartleft*step
      leftend=rc-iendleft*step
c$$$      write(*,*) '----------- Suydam check left'
      call D02EJF(leftstart,leftend,2, vec_l0,outerfunc,D02EJY,tol,        'D', outputRKfuncleft,D02EJW,W2,iw,ifail)
            if (ifail.ne.0) then
         write(*,*) 'RK failed on right'
         stop
         end if

      write(*,*) 'zleft_RK=',vec_l0(1)/(vec_l0(2)/fr(ic-iendleft))
      write(*,*) 'er=',vec_l0(1)
      write(*,*) 'der=',vec_l0(2)/fr(ic-iendleft)
      write(*,*) 'fr=',fr(ic-iendleft)
      zleft=vec_l0(1)/(vec_l0(2)/fr(ic-iendleft))
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccc





c$$$
c$$$      d_right=(zright*der_l_right-er_l_right)/(er_s_right-zright*der_s_right)
c$$$      write(*,*) 'zright*der_l_right=',zright*der_l_right
c$$$      write(*,*) 'zright*der_s_right=',zright*der_s_right
c$$$      write(*,*) 'numerator=',(zright*der_l_right-er_l_right)
c$$$      write(*,*) 'denominator=',(er_s_right-zright*der_s_right)
c$$$      d_left=(zleft*der_l_left-er_l_left)/(er_s_left-zleft*der_s_left)
c$$$      d_prime=d_right+d_left
c$$$      write(*,*) 'd_prime =',d_prime
c$$$      d_right=d_right
c$$$      d_left=d_left
c$$$      write(*,*) 'd_right =',d_right
c$$$      
c$$$      write(*,*) 'd_left =',d_left
c$$$      write(*,*) 'W0right',er_l_right*fr(ic+iendright)*der_s_right-fr(ic+iendright)*der_l_right*er_s_right
c$$$
c$$$
c$$$      write(*,*) 'W0left',er_l_left*fr(ic-iendleft)*der_s_left-fr(ic-iendleft)*der_l_left*er_s_left
c$$$
c$$$
c$$$      write(*,*) 'Wexact',(2.*sigma_s+1.)*dfrc(2)
c$$$      prec=abs(er_l_right*fr(ic+iendright)*der_s_right-fr(ic+iendright)*der_l_right*er_s_right-(2.*sigma_s+1.)*dfrc(2))/abs((2.*sigma_s+1.)*dfrc(2))
c$$$      precleft=abs(er_l_left*fr(ic-iendleft)*der_s_left-fr(ic-iendleft)*der_l_left*er_s_left+(2.*sigma_s+1.)*dfrc(2))/abs((2.*sigma_s+1.)*dfrc(2))


      write(*,*) '***********************************************'


      
      



      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccfunction funfr   cccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function funfr(x)
      implicit none
      real x,funfr
      real step
      real rc
      real, dimension(6)::Cf
      integer size,ic,i,n,last,np
      integer, dimension(1) ::icb
      common /size_block/ size
      real,dimension(1000001) ::r,fr,gr,rb,frb,rbb,rbb2,Bzb
      common /fg_block/ fr,gr,rc,ic,rb,Bzb

     

     
      n=3
      frb=fr
      rbb=rb
      rbb2=rb
c     write(*,*) 'r(icbfirst)',rbb(2)
      icb=minloc(abs(rbb-x))

         
         last=n*(n+1)/2
         np=(n+1)/2

         call E01AAF(rbb((icb(1)-(np-1)):(icb(1)+np)), frb((icb(1)-(np-1)):(icb(1)+np)), Cf,n+1, last, n, x)

         funfr=Cf(last)

      end function
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccfunction fungr   cccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function fungr(x)
      implicit none
      real x,fungr
      real step
      real rc
      real, dimension(6)::Cf
      integer size,ic,i,n,last,np
      common /size_block/ size
      integer, dimension(1) ::icb
      real,dimension(1000001) ::r,fr,gr,rb,grb,rbb,rbb2,Bzb
      common /fg_block/ fr,gr,rc,ic,rb,Bzb


     
      n=3
      grb=gr
      rbb=rb
      rbb2=rb
c     write(*,*) 'r(icbfirst)',rbb(2)
      icb=minloc(abs(rbb-x))

         
         last=n*(n+1)/2
         np=(n+1)/2

         call E01AAF(rbb((icb(1)-(np-1)):(icb(1)+np)), grb((icb(1)-(np-1)):(icb(1)+np)), Cf,n+1, last, n, x)

         fungr=Cf(last)
      end function

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccfunction funBz   cccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function funBz(x)
      implicit none
      real x,funBz
      real step
      real rc
      real, dimension(6)::Cf
      integer size,ic,i,n,last,np
      common /size_block/ size
      integer, dimension(1) ::icb
      real,dimension(1000001) ::r,Bzb,rb,grb,rbb,rbb2,fr,gr
      common /fg_block/ fr,gr,rc,ic,rb,Bzb


     
      n=3
      
      rbb=rb
      rbb2=rb
c     write(*,*) 'r(icbfirst)',rbb(2)
      icb=minloc(abs(rbb-x))

         
         last=n*(n+1)/2
         np=(n+1)/2

         call E01AAF(rbb((icb(1)-(np-1)):(icb(1)+np)), Bzb((icb(1)-(np-1)):(icb(1)+np)), Cf,n+1, last, n, x)

         funBz=Cf(last)
      end function




      


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccfunction det     cccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function det(gamma)
      use innerc_module
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
      common /delta_block/ d_right,d_left,E,F,G,H,K0,Lr_DI,Lr,epsd,rmatch
      if (gamma.le.0) then
         write(*,*) 'negative gamma =',gamma
c     stop
      end if

      rmatch=8
      ifail=0
      write(*,*) '------------------------------------'
      write(*,*) 'gamma_iter =',gamma
      call innerc(cmplx(gamma),E,F,G,H,K0,de,do,ifail,epsd
     1     ,rmatch)
      if (ifail.ne.0)then
         write(*,*) 'innerc did not converge','ifail = ',ifail
         stop
         det=0.
         return
c     stop
         
      end if
      der=real(de*Lr_DI)
      dor=real(do*Lr_DI)
c     der=real(de)
c     dor=real(do)
c     write(*,*) 'Lr_DI=',Lr_DI
c     write(*,*) 'der =',der
c     write(*,*) 'dor =',dor

      det=-(d_right-der)*(d_left-dor)-(d_right-dor)*(d_left-der)
      write(*,*) 'det=',det
c     write(*,*) 'gamma_iter =',gamma
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
ccccccccccccccccccccccmain program     cccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program moddemo
      
      use modtest
      integer size
      common /size_block/ size
      size=1000001
      call main_tok_rev(size,100,5)
c     call testode()
      end program moddemo
