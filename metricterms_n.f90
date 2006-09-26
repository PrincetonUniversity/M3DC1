!===============================================================================
! V1 TERMS
!===============================================================================


! V1umu  (velocity diffusion)
! ===========================
subroutine v1umu(itri,i,j,ssterm, ddterm, rrterm)

  use basic
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm

  real :: temp

  temp = int2(g25(:,OP_GS,i),g25(:,OP_LP,j),weight_25,25)
  temp = -amu*temp

  ssterm(1,1) = ssterm(1,1) -     thimp *dt*temp
  ddterm(1,1) = ddterm(1,1) + (1.-thimp)*dt*temp

end subroutine v1umu


! V1un (vorticity advection)
! ===========================
subroutine v1un(itri,i,j,ssterm, ddterm, rrterm)

  use basic
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm

  real :: temp

  if(idens.eq.0) then
     temp = int2(g25(:,OP_DR,i),g25(:,OP_DR,j),weight_25,25) &
          + int2(g25(:,OP_DZ,i),g25(:,OP_DZ,j),weight_25,25)
  else
     temp = int3(g79(:,OP_DR,i),g79(:,OP_DR,j),nt79(:,OP_1),weight_79,79) &
          + int3(g79(:,OP_DZ,i),g79(:,OP_DZ,j),nt79(:,OP_1),weight_79,79)
  endif

  ssterm(1,1) = ssterm(1,1) + temp
  ddterm(1,1) = ddterm(1,1) + temp  
end subroutine v1un


! V1chin
! ======
real function v1chin(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g

  real :: temp

  if(idens.eq.0) then
     temp = 0.
  else
     temp = int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79) &
          - int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79)
  endif

  v1chin = temp
  return
end function v1chin


! V1uun
! =====
subroutine v1uun(itri,i,j,ssterm, ddterm, rrterm)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm

  integer :: i1
  real :: temp
           
  if(idens.eq.0) then

     ! Linearize in u(1)
     ! ~~~~~~~~~~~~~~~~~
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        temp = int4(r_79,g79(:,OP_DR,i),g79(:,OP_DZ,j),ph079(:,OP_LP),weight_79,79) &
             - int4(r_79,g79(:,OP_DZ,i),g79(:,OP_DR,j),ph079(:,OP_LP),weight_79,79)
        if(itor.eq.1) then
           temp = temp &
             + int3(g79(:,OP_1,i ),g79(:,OP_DZ, j),ph079(:,OP_LP),weight_79,79) &
             + int3(g79(:,OP_1,i ),g79(:,OP_DZZ,j),ph079(:,OP_DZ),weight_79,79) &
             + int3(g79(:,OP_1,i ),g79(:,OP_DRZ,j),ph079(:,OP_DR),weight_79,79) &
             + int3(g79(:,OP_DZ,i),g79(:,OP_DZ, j),ph079(:,OP_DZ),weight_79,79) &
             + int3(g79(:,OP_DR,i),g79(:,OP_DZ, j),ph079(:,OP_DR),weight_79,79)
        endif
        ssterm(1,1) = ssterm(1,1) -     thimp *dt*temp
        ddterm(1,1) = ddterm(1,1) + (1.-thimp)*dt*temp
     endif

     temp = int4(r_79,g79(:,OP_DR,i),g79(:,OP_DZ,j),ph179(:,OP_LP),weight_79,79) &
          - int4(r_79,g79(:,OP_DZ,i),g79(:,OP_DR,j),ph179(:,OP_LP),weight_79,79)
     if(itor.eq.1) then
        temp = temp &
          + int3(g79(:,OP_1, i),g79(:,OP_DZ, j),ph179(:,OP_LP),weight_79,79) &
          + int3(g79(:,OP_1, i),g79(:,OP_DZZ,j),ph179(:,OP_DZ),weight_79,79) &
          + int3(g79(:,OP_1, i),g79(:,OP_DRZ,j),ph179(:,OP_DR),weight_79,79) &
          + int3(g79(:,OP_DZ,i),g79(:,OP_DZ, j),ph179(:,OP_DZ),weight_79,79) &
          + int3(g79(:,OP_DR,i),g79(:,OP_DZ, j),ph179(:,OP_DR),weight_79,79)
     endif
     ssterm(1,1) = ssterm(1,1) -      thimp *dt*temp
     ddterm(1,1) = ddterm(1,1) + (0.5-thimp)*dt*temp
             
     ! Linearize in u(2)
     ! ~~~~~~~~~~~~~~~~~
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        temp = int4(r_79,g79(:,OP_DR,i),ph079(:,OP_DZ),g79(:,OP_LP,j),weight_79,79) &
             - int4(r_79,g79(:,OP_DZ,i),ph079(:,OP_DR),g79(:,OP_LP,j),weight_79,79)
        if(itor.eq.1) then
           temp = temp &
             + int3(g79(:,OP_1, i),ph079(:,OP_DZ ),g79(:,OP_LP,j),weight_79,79) &
             + int3(g79(:,OP_1, i),ph079(:,OP_DZZ),g79(:,OP_DZ,j),weight_79,79) &
             + int3(g79(:,OP_1, i),ph079(:,OP_DRZ),g79(:,OP_DR,j),weight_79,79) &
             + int3(g79(:,OP_DZ,i),ph079(:,OP_DZ ),g79(:,OP_DZ,j),weight_79,79) &
             + int3(g79(:,OP_DR,i),ph079(:,OP_DZ ),g79(:,OP_DR,j),weight_79,79)
        endif
        ssterm(1,1) = ssterm(1,1) -     thimp *dt*temp
        ddterm(1,1) = ddterm(1,1) + (1.-thimp)*dt*temp
     endif
              
     temp = int4(r_79,g79(:,OP_DR,i),ph179(:,OP_DZ),g79(:,OP_LP,j),weight_79,79) &
          - int4(r_79,g79(:,OP_DZ,i),ph179(:,OP_DR),g79(:,OP_LP,j),weight_79,79)
     if(itor.eq.1) then
        temp = temp &
          + int3(g79(:,OP_1, i),ph179(:,OP_DZ),g79(:,OP_LP,j),weight_79,79) &
          + int3(g79(:,OP_1, i),ph179(:,OP_DZZ),g79(:,OP_DZ,j),weight_79,79) &
          + int3(g79(:,OP_1, i),ph179(:,OP_DRZ),g79(:,OP_DR,j),weight_79,79) &
          + int3(g79(:,OP_DZ,i),ph179(:,OP_DZ),g79(:,OP_DZ,j),weight_79,79) &
          + int3(g79(:,OP_DR,i),ph179(:,OP_DZ),g79(:,OP_DR,j),weight_79,79)
     endif
     ssterm(1,1) = ssterm(1,1) -      thimp *dt*temp
     ddterm(1,1) = ddterm(1,1) + (0.5-thimp)*dt*temp
              
  else

     ! Linearize in u(1)
     ! ~~~~~~~~~~~~~~~~~
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        temp = -int5(r_79,g79(:,OP_1, i),g79(:,OP_DZ,j),ph079(:,OP_DRR),nt79(:,OP_DR),weight_79,79) &
             +  int5(r_79,g79(:,OP_1, i),g79(:,OP_DR,j),ph079(:,OP_DRZ),nt79(:,OP_DR),weight_79,79) &
             -  int5(r_79,g79(:,OP_1, i),g79(:,OP_DZ,j),ph079(:,OP_DRZ),nt79(:,OP_DZ),weight_79,79) &
             +  int5(r_79,g79(:,OP_1, i),g79(:,OP_DR,j),ph079(:,OP_DZZ),nt79(:,OP_DZ),weight_79,79) &
             +  int5(r_79,g79(:,OP_1, i),g79(:,OP_DZ,j),ph079(:,OP_LP ),nt79(:,OP_DR),weight_79,79) &
             -  int5(r_79,g79(:,OP_1, i),g79(:,OP_DR,j),ph079(:,OP_LP ),nt79(:,OP_DZ),weight_79,79) &
             +  int5(r_79,g79(:,OP_DR,i),g79(:,OP_DZ,j),ph079(:,OP_LP ),nt79(:,OP_1),weight_79,79) &
             -  int5(r_79,g79(:,OP_DZ,i),g79(:,OP_DR,j),ph079(:,OP_LP ),nt79(:,OP_1),weight_79,79)
        if(itor.eq.1) then
           temp = temp &
             +  int4(g79(:,OP_1 ,i),g79(:,OP_DZ, j),ph079(:,OP_LP),nt79(:,OP_1),weight_79,79) &
             +  int4(g79(:,OP_1 ,i),g79(:,OP_DZZ,j),ph079(:,OP_DR),nt79(:,OP_1),weight_79,79) &
             +  int4(g79(:,OP_1 ,i),g79(:,OP_DRZ,j),ph079(:,OP_DZ),nt79(:,OP_1),weight_79,79) &
             +  int4(g79(:,OP_DZ,i),g79(:,OP_DZ, j),ph079(:,OP_DR),nt79(:,OP_1),weight_79,79) &
             +  int4(g79(:,OP_DZ,i),g79(:,OP_DR, j),ph079(:,OP_DZ),nt79(:,OP_1),weight_79,79)
        endif
        ssterm(1,1) = ssterm(1,1) -     thimp *dt*temp
        ddterm(1,1) = ddterm(1,1) + (1.-thimp)*dt*temp
     endif

     temp = -int5(r_79,g79(:,OP_1, i),g79(:,OP_DZ,j),ph179(:,OP_DRR),nt79(:,OP_DR),weight_79,79) &
          +  int5(r_79,g79(:,OP_1, i),g79(:,OP_DR,j),ph179(:,OP_DRZ),nt79(:,OP_DR),weight_79,79) &
          -  int5(r_79,g79(:,OP_1, i),g79(:,OP_DZ,j),ph179(:,OP_DRZ),nt79(:,OP_DZ),weight_79,79) &
          +  int5(r_79,g79(:,OP_1, i),g79(:,OP_DR,j),ph179(:,OP_DZZ),nt79(:,OP_DZ),weight_79,79) &
          +  int5(r_79,g79(:,OP_1, i),g79(:,OP_DZ,j),ph179(:,OP_LP ),nt79(:,OP_DR),weight_79,79) &
          -  int5(r_79,g79(:,OP_1, i),g79(:,OP_DR,j),ph179(:,OP_LP ),nt79(:,OP_DZ),weight_79,79) &
          +  int5(r_79,g79(:,OP_DR,i),g79(:,OP_DZ,j),ph179(:,OP_LP ),nt79(:,OP_1),weight_79,79) &
          -  int5(r_79,g79(:,OP_DZ,i),g79(:,OP_DR,j),ph179(:,OP_LP ),nt79(:,OP_1),weight_79,79)
     if(itor.eq.1) then
        temp = temp &
          +  int4(g79(:,OP_1 ,i),g79(:,OP_DZ, j),ph179(:,OP_LP),nt79(:,OP_1),weight_79,79) &
          +  int4(g79(:,OP_1 ,i),g79(:,OP_DZZ,j),ph179(:,OP_DR),nt79(:,OP_1),weight_79,79) &
          +  int4(g79(:,OP_1 ,i),g79(:,OP_DRZ,j),ph179(:,OP_DZ),nt79(:,OP_1),weight_79,79) &
          +  int4(g79(:,OP_DZ,i),g79(:,OP_DZ, j),ph179(:,OP_DR),nt79(:,OP_1),weight_79,79) &
          +  int4(g79(:,OP_DZ,i),g79(:,OP_DR, j),ph179(:,OP_DZ),nt79(:,OP_1),weight_79,79)
     endif
     ssterm(1,1) = ssterm(1,1) -      thimp *dt*temp
     ddterm(1,1) = ddterm(1,1) + (0.5-thimp)*dt*temp

     ! Linearize in u(2)
     ! ~~~~~~~~~~~~~~~~~
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        temp = -int5(r_79,g79(:,OP_1, i),ph079(:,OP_DZ),g79(:,OP_DRR,j),nt79(:,OP_DR),weight_79,79) &
             +  int5(r_79,g79(:,OP_1, i),ph079(:,OP_DR),g79(:,OP_DRZ,j),nt79(:,OP_DR),weight_79,79) &
             -  int5(r_79,g79(:,OP_1, i),ph079(:,OP_DZ),g79(:,OP_DRZ,j),nt79(:,OP_DZ),weight_79,79) &
             +  int5(r_79,g79(:,OP_1, i),ph079(:,OP_DR),g79(:,OP_DZZ,j),nt79(:,OP_DZ),weight_79,79) &
             +  int5(r_79,g79(:,OP_1, i),ph079(:,OP_DZ),g79(:,OP_LP, j),nt79(:,OP_DR),weight_79,79) &
             -  int5(r_79,g79(:,OP_1, i),ph079(:,OP_DR),g79(:,OP_LP, j),nt79(:,OP_DZ),weight_79,79) &
             +  int5(r_79,g79(:,OP_DR,i),ph079(:,OP_DZ),g79(:,OP_LP, j),nt79(:,OP_1),weight_79,79) &
             -  int5(r_79,g79(:,OP_DZ,i),ph079(:,OP_DR),g79(:,OP_LP, j),nt79(:,OP_1),weight_79,79)
        if(itor.eq.1) then
           temp = temp &
             +  int4(g79(:,OP_1 ,i),ph079(:,OP_DZ ),g79(:,OP_LP,j),nt79(:,OP_1),weight_79,79) &
             +  int4(g79(:,OP_1 ,i),ph079(:,OP_DZZ),g79(:,OP_DR,j),nt79(:,OP_1),weight_79,79) &
             +  int4(g79(:,OP_1 ,i),ph079(:,OP_DRZ),g79(:,OP_DZ,j),nt79(:,OP_1),weight_79,79) &
             +  int4(g79(:,OP_DZ,i),ph079(:,OP_DZ ),g79(:,OP_DR,j),nt79(:,OP_1),weight_79,79) &
             +  int4(g79(:,OP_DZ,i),ph079(:,OP_DR ),g79(:,OP_DZ,j),nt79(:,OP_1),weight_79,79)
        endif
        ssterm(1,1) = ssterm(1,1) -     thimp *dt*temp
        ddterm(1,1) = ddterm(1,1) + (1.-thimp)*dt*temp
     endif

     temp = -int5(r_79,g79(:,OP_1, i),ph179(:,OP_DZ),g79(:,OP_DRR,j),nt79(:,OP_DR),weight_79,79) &
          +  int5(r_79,g79(:,OP_1, i),ph179(:,OP_DR),g79(:,OP_DRZ,j),nt79(:,OP_DR),weight_79,79) &
          -  int5(r_79,g79(:,OP_1, i),ph179(:,OP_DZ),g79(:,OP_DRZ,j),nt79(:,OP_DZ),weight_79,79) &
          +  int5(r_79,g79(:,OP_1, i),ph179(:,OP_DR),g79(:,OP_DZZ,j),nt79(:,OP_DZ),weight_79,79) &
          +  int5(r_79,g79(:,OP_1, i),ph179(:,OP_DZ),g79(:,OP_LP, j),nt79(:,OP_DR),weight_79,79) &
          -  int5(r_79,g79(:,OP_1, i),ph179(:,OP_DR),g79(:,OP_LP, j),nt79(:,OP_DZ),weight_79,79) &
          +  int5(r_79,g79(:,OP_DR,i),ph179(:,OP_DZ),g79(:,OP_LP, j),nt79(:,OP_1),weight_79,79) &
          -  int5(r_79,g79(:,OP_DZ,i),ph179(:,OP_DR),g79(:,OP_LP, j),nt79(:,OP_1),weight_79,79)
     if(itor.eq.1) then
        temp = temp &
          +  int4(g79(:,OP_1 ,i),ph179(:,OP_DZ ),g79(:,OP_LP,j),nt79(:,OP_1),weight_79,79) &
          +  int4(g79(:,OP_1 ,i),ph179(:,OP_DZZ),g79(:,OP_DR,j),nt79(:,OP_1),weight_79,79) &
          +  int4(g79(:,OP_1 ,i),ph179(:,OP_DRZ),g79(:,OP_DZ,j),nt79(:,OP_1),weight_79,79) &
          +  int4(g79(:,OP_DZ,i),ph179(:,OP_DZ ),g79(:,OP_DR,j),nt79(:,OP_1),weight_79,79) &
          +  int4(g79(:,OP_DZ,i),ph179(:,OP_DR ),g79(:,OP_DZ,j),nt79(:,OP_1),weight_79,79)
     endif
     ssterm(1,1) = ssterm(1,1) -      thimp *dt*temp
     ddterm(1,1) = ddterm(1,1) + (0.5-thimp)*dt*temp

     ! Linearize in n
     ! ~~~~~~~~~~~~~~
     if(linear.eq.1 .or. eqsubtract.eq.1) then
        i1 = isvaln(itri,i)

        temp = -int5(r_79,g79(:,OP_1, i),ph079(:,OP_DZ),ph079(:,OP_DRR),nt79(:,OP_DR),weight_79,79) &
             +  int5(r_79,g79(:,OP_1, i),ph079(:,OP_DR),ph079(:,OP_DRZ),nt79(:,OP_DR),weight_79,79) &
             -  int5(r_79,g79(:,OP_1, i),ph079(:,OP_DZ),ph079(:,OP_DRZ),nt79(:,OP_DZ),weight_79,79) &
             +  int5(r_79,g79(:,OP_1, i),ph079(:,OP_DR),ph079(:,OP_DZZ),nt79(:,OP_DZ),weight_79,79) &
             +  int5(r_79,g79(:,OP_1, i),ph079(:,OP_DZ),ph079(:,OP_LP ),nt79(:,OP_DR),weight_79,79) &
             -  int5(r_79,g79(:,OP_1, i),ph079(:,OP_DR),ph079(:,OP_LP ),nt79(:,OP_DZ),weight_79,79) &
             +  int5(r_79,g79(:,OP_DR,i),ph079(:,OP_DZ),ph079(:,OP_LP ),nt79(:,OP_1),weight_79,79) &
             -  int5(r_79,g79(:,OP_DZ,i),ph079(:,OP_DR),ph079(:,OP_LP ),nt79(:,OP_1),weight_79,79)
        if(itor.eq.1) then
           temp = temp &
             +  int4(g79(:,OP_1 ,i),ph079(:,OP_DZ ),ph079(:,OP_LP),nt79(:,OP_1),weight_79,79) &
             +  int4(g79(:,OP_1 ,i),ph079(:,OP_DZZ),ph079(:,OP_DR),nt79(:,OP_1),weight_79,79) &
             +  int4(g79(:,OP_1 ,i),ph079(:,OP_DRZ),ph079(:,OP_DZ),nt79(:,OP_1),weight_79,79) &
             +  int4(g79(:,OP_DZ,i),ph079(:,OP_DZ ),ph079(:,OP_DR),nt79(:,OP_1),weight_79,79) &
             +  int4(g79(:,OP_DZ,i),ph079(:,OP_DR ),ph079(:,OP_DZ),nt79(:,OP_1),weight_79,79)
        endif
        r4(i1) = r4(i1) + dt*temp
     endif

  endif

end subroutine v1uun


! v1uchin
! ~~~~~~~
real function v1uchin(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h

  real :: temp

  if(idens.eq.0) then
     temp =-int3(e(:,OP_DZ),f(:,OP_LP),g(:,OP_DZ),weight_79,79) &
          - int3(e(:,OP_DR),f(:,OP_LP),g(:,OP_DR),weight_79,79)
  else 
     temp =-int4(e(:,OP_DZ),f(:,OP_LP ),g(:,OP_DZ ),h(:,OP_1 ),weight_79,79) &
          - int4(e(:,OP_DR),f(:,OP_LP ),g(:,OP_DR ),h(:,OP_1 ),weight_79,79) &
          - int4(e(:,OP_1 ),f(:,OP_LP ),g(:,OP_DZ ),h(:,OP_DZ),weight_79,79) &
          - int4(e(:,OP_1 ),f(:,OP_LP ),g(:,OP_DR ),h(:,OP_DR),weight_79,79) &
          + int4(e(:,OP_1 ),f(:,OP_DRR),g(:,OP_DR ),h(:,OP_DR),weight_79,79) &
          + int4(e(:,OP_1 ),f(:,OP_DRZ),g(:,OP_DR ),h(:,OP_DZ),weight_79,79) &
          + int4(e(:,OP_1 ),f(:,OP_DRZ),g(:,OP_DZ ),h(:,OP_DR),weight_79,79) &
          + int4(e(:,OP_1 ),f(:,OP_DZZ),g(:,OP_DZ ),h(:,OP_DZ),weight_79,79) &
          + int4(e(:,OP_1 ),f(:,OP_DZ ),g(:,OP_DRR),h(:,OP_DZ),weight_79,79) &
          - int4(e(:,OP_1 ),f(:,OP_DZ ),g(:,OP_DRZ),h(:,OP_DR),weight_79,79) &
          - int4(e(:,OP_1 ),f(:,OP_DR ),g(:,OP_DRZ),h(:,OP_DZ),weight_79,79) &
          + int4(e(:,OP_1 ),f(:,OP_DR ),g(:,OP_DZZ),h(:,OP_DR),weight_79,79)
     if(itor.eq.1) then
        temp = temp &
             + int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_DZ),weight_79,79) &
             + int5(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_DR),weight_79,79)
     endif
  endif

  v1uchin = temp
end function v1uchin


! v1chichin
! =========
real function v1chichin(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h

  real :: temp

  if(idens.eq.0) then
     temp = 0
  else
     temp = -0.5* &
          (int5(ri_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_DR),weight_79,79) &
          +int5(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DR),h(:,OP_DR),weight_79,79) &
          -int5(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_DZ),weight_79,79) &
          -int5(ri_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DR),h(:,OP_DZ),weight_79,79))
  endif
  
  v1chichin = temp

  return
end function v1chichin



! V1psisb1
! ========
subroutine v1psisb1(itri,i,j,ssterm, ddterm, rrterm)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm

  integer :: i1
  real :: temp

  temp = int4(ri3_79,g79(:,OP_DZ,i),g79(:,OP_GS,j),sb179(:,OP_DR),weight_79,79) &
       - int4(ri3_79,g79(:,OP_DR,i),g79(:,OP_GS,j),sb179(:,OP_DZ),weight_79,79) &
       + int4(ri3_79,g79(:,OP_DZ,i),g79(:,OP_DR,j),sb179(:,OP_GS),weight_79,79) &
       - int4(ri3_79,g79(:,OP_DR,i),g79(:,OP_DZ,j),sb179(:,OP_GS),weight_79,79)

  rrterm(1,1) = rrterm(1,1) + thimp *dt*dt*temp

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     i1 = isvaln(itri,i)

     temp = int4(ri3_79,g79(:,OP_DZ,i),ph079(:,OP_GS),sb179(:,OP_DR),weight_79,79) &
          - int4(ri3_79,g79(:,OP_DR,i),ph079(:,OP_GS),sb179(:,OP_DZ),weight_79,79) &
          + int4(ri3_79,g79(:,OP_DZ,i),ph079(:,OP_DR),sb179(:,OP_GS),weight_79,79) &
          - int4(ri3_79,g79(:,OP_DR,i),ph079(:,OP_DZ),sb179(:,OP_GS),weight_79,79)

     r4(i1) = r4(i1) + thimp*dt*dt*temp
  end if

end subroutine v1psisb1


! V1psipsi
! ========
subroutine v1psipsi(itri,i,j,ssterm, ddterm, rrterm)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm

  integer :: i1
  real :: temp
          
  temp = int4(ri3_79,g79(:,OP_DZ,i),g79(:,OP_GS,j),ps179(:,OP_DR),weight_79,79) &
       - int4(ri3_79,g79(:,OP_DR,i),g79(:,OP_GS,j),ps179(:,OP_DZ),weight_79,79) &
       + int4(ri3_79,g79(:,OP_DZ,i),ps179(:,OP_GS),g79(:,OP_DR,j),weight_79,79) &
       - int4(ri3_79,g79(:,OP_DR,i),ps179(:,OP_GS),g79(:,OP_DZ,j),weight_79,79)
  rrterm(1,1) = rrterm(1,1) + dt*temp/2.

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     i1 = isval1(itri,i)

     temp = int4(ri3_79,g79(:,OP_DZ,i),g79(:,OP_GS,j),ps079(:,OP_DR),weight_79,79) &
          - int4(ri3_79,g79(:,OP_DR,i),g79(:,OP_GS,j),ps079(:,OP_DZ),weight_79,79) &
          + int4(ri3_79,g79(:,OP_DZ,i),ps079(:,OP_GS),g79(:,OP_DR,j),weight_79,79) &
          - int4(ri3_79,g79(:,OP_DR,i),ps079(:,OP_GS),g79(:,OP_DZ,j),weight_79,79)
     rrterm(1,1) = rrterm(1,1) + dt*temp

     ! EQUILIBRIUM TERM
     temp = int4(ri3_79,g79(:,OP_DZ,i),ps079(:,OP_GS),ps079(:,OP_DR),weight_79,79) &
          - int4(ri3_79,g79(:,OP_DR,i),ps079(:,OP_GS),ps079(:,OP_DZ),weight_79,79)
     r4(i1) = r4(i1) + dt*temp
  endif
  

end subroutine v1psipsi


! V1upsipsi
! =========
subroutine v1upsipsi(itri,i,j,ssterm, ddterm, rrterm)

  use basic
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm

  real :: temp

  ! Linearize in u
  ! ~~~~~~~~~~~~~~
           
  ! |psi(1), u|
  temp79a = g79(:,OP_DR, j)*pst79(:,OP_DZ ) - g79(:,OP_DZ, j)*pst79(:,OP_DR)
  ! |psi(1), u|,r
  temp79b = g79(:,OP_DRR,j)*pst79(:,OP_DZ ) - g79(:,OP_DRZ,j)*pst79(:,OP_DR ) &
       +    g79(:,OP_DR, j)*pst79(:,OP_DRZ) - g79(:,OP_DZ, j)*pst79(:,OP_DRR)
  ! |psi(1), u|,z
  temp79c = g79(:,OP_DRZ,j)*pst79(:,OP_DZ ) - g79(:,OP_DZZ,j)*pst79(:,OP_DR ) &
       +    g79(:,OP_DR, j)*pst79(:,OP_DZZ) - g79(:,OP_DZ, j)*pst79(:,OP_DRZ)

  ! |psi(2), nu|
  temp79d = g79(:,OP_DR, i)*pst79(:,OP_DZ ) - g79(:,OP_DZ, i)*pst79(:,OP_DR)
  ! |psi(2), nu|,r
  temp79e = g79(:,OP_DRR,i)*pst79(:,OP_DZ ) - g79(:,OP_DRZ,i)*pst79(:,OP_DR ) &
       +    g79(:,OP_DR, i)*pst79(:,OP_DRZ) - g79(:,OP_DZ, i)*pst79(:,OP_DRR)
  ! |psi(2), nu|,z
  temp79f = g79(:,OP_DRZ,i)*pst79(:,OP_DZ ) - g79(:,OP_DZZ,i)*pst79(:,OP_DR ) &
       +    g79(:,OP_DR, i)*pst79(:,OP_DZZ) - g79(:,OP_DZ, i)*pst79(:,OP_DRZ)
  
  temp = int4(ri2_79,g79(:,OP_DR,i),temp79c,pst79(:,OP_GS),weight_79,79) &
       - int4(ri2_79,g79(:,OP_DZ,i),temp79b,pst79(:,OP_GS),weight_79,79) &
       - int3(ri2_79,temp79b,temp79e,weight_79,79) &
       - int3(ri2_79,temp79c,temp79f,weight_79,79)
  if(itor.eq.1) then
     temp = temp &
       + int3(ri3_79,temp79b,temp79d,weight_79,79) &
       - int3(ri3_79,temp79a,temp79e,weight_79,79) &
       - int4(ri3_79,g79(:,OP_DZ,i),temp79a,pst79(:,OP_GS),weight_79,79) &
       + int3(ri4_79,temp79a,temp79d,weight_79,79)
  endif
  ssterm(1,1) = ssterm(1,1) - thimp*    thimp *dt*dt*temp
  ddterm(1,1) = ddterm(1,1) + thimp*(1.-thimp)*dt*dt*temp

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     ! Linearize in psi(1)
     ! ~~~~~~~~~~~~~~~~~~~
     tm79 = ps079 + ps179/2.
     
     ! |psi(1), u|   
     temp79a = ph079(:,OP_DR)*g79(:,OP_DZ,j) - ph079(:,OP_DZ)*g79(:,OP_DR,j)
     ! |psi(1), u|,r
     temp79b = ph079(:,OP_DRR)*g79(:,OP_DZ, j) - ph079(:,OP_DRZ)*g79(:,OP_DR, j) &
          +    ph079(:,OP_DR )*g79(:,OP_DRZ,j) - ph079(:,OP_DZ )*g79(:,OP_DRR,j)
     ! |psi(1), u|,z
     temp79c = ph079(:,OP_DRZ)*g79(:,OP_DZ, j) - ph079(:,OP_DZZ)*g79(:,OP_DR, j) &
          +    ph079(:,OP_DR )*g79(:,OP_DZZ,j) - ph079(:,OP_DZ )*g79(:,OP_DRZ,j)
     
     ! |psi(2), nu|
     temp79d = g79(:,OP_DR,i)*tm79(:,OP_DZ) - g79(:,OP_DZ,i)*tm79(:,OP_DR)
     ! |psi(2), nu|,r
     temp79e = g79(:,OP_DRR,i)*tm79(:,OP_DZ ) - g79(:,OP_DRZ,i)*tm79(:,OP_DR ) &
          +    g79(:,OP_DR, i)*tm79(:,OP_DRZ) - g79(:,OP_DZ, i)*tm79(:,OP_DRR)
     ! |psi(2), nu|,z
     temp79f = g79(:,OP_DRZ,i)*tm79(:,OP_DZ ) - g79(:,OP_DZZ,i)*tm79(:,OP_DR ) &
          +    g79(:,OP_DR, i)*tm79(:,OP_DZZ) - g79(:,OP_DZ, i)*tm79(:,OP_DRZ)
     
     temp = int4(ri2_79,g79(:,OP_DR,i),temp79c,tm79(:,OP_GS),weight_79,79) &
          - int4(ri2_79,g79(:,OP_DZ,i),temp79b,tm79(:,OP_GS),weight_79,79) &
          - int3(ri2_79,temp79b,temp79e,weight_79,79) &
          - int3(ri2_79,temp79c,temp79f,weight_79,79)
     if(itor.eq.1) then
        temp = temp &
          + int3(ri3_79,temp79b,temp79d,weight_79,79) &
          - int3(ri3_79,temp79a,temp79e,weight_79,79) &
          - int4(ri3_79,g79(:,OP_DZ,i),temp79a,tm79(:,OP_GS),weight_79,79) &
          + int3(ri4_79,temp79a,temp79d,weight_79,79)
     endif
     rrterm(1,1) = rrterm(1,1) + thimp*dt*dt*temp
     

     ! Linearize in psi(2)
     ! ~~~~~~~~~~~~~~~~~~~

     ! |psi(1), u|
     temp79a = ph079(:,OP_DR)*tm79(:,OP_DZ) - ph079(:,OP_DZ)*tm79(:,OP_DR)
     ! |psi(1), u|,r
     temp79b = ph079(:,OP_DRR)*tm79(:,OP_DZ ) - ph079(:,OP_DRZ)*tm79(:,OP_DR ) &
          +    ph079(:,OP_DR )*tm79(:,OP_DRZ) - ph079(:,OP_DZ )*tm79(:,OP_DRR)
     ! |psi(1), u|,z
     temp79c = ph079(:,OP_DRZ)*tm79(:,OP_DZ ) - ph079(:,OP_DZZ)*tm79(:,OP_DR ) &
          +    ph079(:,OP_DR )*tm79(:,OP_DZZ) - ph079(:,OP_DZ )*tm79(:,OP_DRZ)
     
     ! |psi(2), nu|
     temp79d = g79(:,OP_DR,i)*g79(:,OP_DZ, j) - g79(:,OP_DZ,i)*g79(:,OP_DR,j)
     ! |psi(2), nu|,r
     temp79e = g79(:,OP_DRR,i)*g79(:,OP_DZ, j) - g79(:,OP_DRZ,i)*g79(:,OP_DR, j) &
          +    g79(:,OP_DR, i)*g79(:,OP_DRZ,j) - g79(:,OP_DZ, i)*g79(:,OP_DRR,j)
     ! |psi(2), nu|,z
     temp79f = g79(:,OP_DRZ,i)*g79(:,OP_DZ, j) - g79(:,OP_DZZ,i)*g79(:,OP_DR, j) &
          +    g79(:,OP_DR, i)*g79(:,OP_DZZ,j) - g79(:,OP_DZ, i)*g79(:,OP_DRZ,j)
     
     temp = int4(ri2_79,g79(:,OP_DR,i),temp79c,g79(:,OP_GS,j),weight_79,79) &
          - int4(ri2_79,g79(:,OP_DZ,i),temp79b,g79(:,OP_GS,j),weight_79,79) &
          - int3(ri2_79,temp79b,temp79e,weight_79,79) &
          - int3(ri2_79,temp79c,temp79f,weight_79,79)
     if(itor.eq.1) then
        temp = temp &
          + int3(ri3_79,temp79b,temp79d,weight_79,79) &
          - int3(ri3_79,temp79a,temp79e,weight_79,79) &
          - int4(ri3_79,g79(:,OP_DZ,i),temp79a,g79(:,OP_GS,j),weight_79,79) &
          + int3(ri4_79,temp79a,temp79d,weight_79,79)
     endif
     rrterm(1,1) = rrterm(1,1) + thimp*dt*dt*temp
  endif
end subroutine v1upsipsi


! V1vpsib
! =========================================================================
subroutine v1vpsib(itri,i,j,ssterm, ddterm, rrterm)

  use basic
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm

  real :: temp

  if(itor.eq.0) return

  ! Linearize in v
  ! ~~~~~~~~~~~~~~
  temp = int5(ri6_79,g79(:,OP_DZ,i),g79(:,OP_DZ,j),pst79(:,OP_DR),bzt79(:,OP_1),weight_79,79) &
       - int5(ri6_79,g79(:,OP_DZ,i),g79(:,OP_DR,j),pst79(:,OP_DZ),bzt79(:,OP_1),weight_79,79) &
       + int5(ri7_79,g79(:,OP_DZ,i),g79(:,OP_1,j),pst79(:,OP_DZ),bzt79(:,OP_1),weight_79,79)

  ssterm(1,2) = ssterm(1,2) - thimp*    thimp *dt*dt*temp
  ddterm(1,2) = ddterm(1,2) + thimp*(1.-thimp)*dt*dt*temp

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     ! Linearize in psi
     ! ~~~~~~~~~~~~~~~~
     temp = int5(ri6_79,g79(:,OP_DZ,i),vz079(:,OP_DZ),g79(:,OP_DR,j),bz179(:,OP_1),weight_79,79)/2. &
          - int5(ri6_79,g79(:,OP_DZ,i),vz079(:,OP_DR),g79(:,OP_DZ,j),bz179(:,OP_1),weight_79,79)/2. &
          + int5(ri6_79,g79(:,OP_DZ,i),vz079(:,OP_DZ),g79(:,OP_DR,j),bz079(:,OP_1),weight_79,79) &
          - int5(ri6_79,g79(:,OP_DZ,i),vz079(:,OP_DR),g79(:,OP_DZ,j),bz079(:,OP_1),weight_79,79) &
          + int5(ri7_79,g79(:,OP_DZ,i),vz079(:,OP_1),g79(:,OP_DZ,j),bz179(:,OP_1),weight_79,79)/2. &
          + int5(ri7_79,g79(:,OP_DZ,i),vz079(:,OP_1),g79(:,OP_DZ,j),bz079(:,OP_1),weight_79,79)

     rrterm(1,1) = rrterm(1,1) + thimp*dt*dt*temp

     ! Linearize in bz
     ! ~~~~~~~~~~~~~~~
     temp = int5(ri6_79,g79(:,OP_DZ,i),vz079(:,OP_DZ),ps179(:,OP_DR),g79(:,OP_1,j),weight_79,79)/2. &
          - int5(ri6_79,g79(:,OP_DZ,i),vz079(:,OP_DR),ps179(:,OP_DZ),g79(:,OP_1,j),weight_79,79)/2. &
          + int5(ri6_79,g79(:,OP_DZ,i),vz079(:,OP_DZ),ps079(:,OP_DR),g79(:,OP_1,j),weight_79,79) &
          - int5(ri6_79,g79(:,OP_DZ,i),vz079(:,OP_DR),ps079(:,OP_DZ),g79(:,OP_1,j),weight_79,79) &
          + int5(ri7_79,g79(:,OP_DZ,i),vz079(:,OP_1),ps179(:,OP_DZ),g79(:,OP_1,j),weight_79,79)/2. &
          + int5(ri7_79,g79(:,OP_DZ,i),vz079(:,OP_1),ps079(:,OP_DZ),g79(:,OP_1,j),weight_79,79)

     rrterm(1,2) = rrterm(1,2) + thimp*dt*dt*temp

  endif


end subroutine v1vpsib


! V1bb 
! ====
subroutine v1bb(itri,i,j,ssterm, ddterm, rrterm)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm

  integer :: i1
  real :: temp

  if(itor.eq.0) return

  temp = -2.* &
       (int4(ri4_79,g79(:,OP_1,i),g79(:,OP_1,j),bz179(:,OP_DZ),weight_79,79) &
       +int4(ri4_79,g79(:,OP_1,i),bz179(:,OP_1),g79(:,OP_DZ,j),weight_79,79))

  rrterm(1,2) = rrterm(1,2) + dt*temp/2.


  if(linear.eq.1 .or. eqsubtract.eq.1) then
     temp = -2.* &
          (int4(ri4_79,g79(:,OP_1,i),g79(:,OP_1,j),bz079(:,OP_DZ),weight_79,79) &
          +int4(ri4_79,g79(:,OP_1,i),bz079(:,OP_1),g79(:,OP_DZ,j),weight_79,79))
     
     rrterm(1,2) = rrterm(1,2) + dt*temp


     ! EQUILIBRIUM TERM
     i1 = isvaln(itri,i)

     temp = -2.*int4(ri4_79,g79(:,OP_1,i),bz079(:,OP_1),bz079(:,OP_DZ),weight_79,79)

     r4(i1) = r4(i1) + dt*temp
  endif
end subroutine v1bb


! V1ubb 
! =====
real function v1ubb(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp

  if(itor.eq.0) then
     temp = 0.
  else
     temp = 2.* &
          (int5(ri5_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DR),g(:,OP_1),weight_79,79) &
          -int5(ri5_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ),g(:,OP_1),weight_79,79))
  endif
  
  v1ubb = temp
  return
end function v1ubb


! V1chipsipsi
! ===========
real function v1chipsipsi(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  real :: temp

  temp79a = f(:,OP_DZ)*g(:,OP_DZ) + f(:,OP_DR)*g(:,OP_DR)
  temp79b = f(:,OP_DRZ)*g(:,OP_DZ ) + f(:,OP_DRR)*g(:,OP_DR ) &
       +    f(:,OP_DZ )*g(:,OP_DRZ) + f(:,OP_DR )*g(:,OP_DRR)
  temp79c = f(:,OP_DZZ)*g(:,OP_DZ ) + f(:,OP_DRZ)*g(:,OP_DR ) &
       +    f(:,OP_DZ )*g(:,OP_DZZ) + f(:,OP_DR )*g(:,OP_DRZ)

  temp = int4(ri2_79,h(:,OP_GS ),e(:,OP_DR ),temp79c,weight_79,79) &
       - int4(ri2_79,h(:,OP_GS ),e(:,OP_DZ ),temp79b,weight_79,79) &
       + int4(ri3_79,e(:,OP_DZZ),h(:,OP_DR ),temp79c,weight_79,79) &
       - int4(ri3_79,e(:,OP_DRZ),h(:,OP_DZ ),temp79c,weight_79,79) &
       + int4(ri3_79,e(:,OP_DZ ),h(:,OP_DRZ),temp79c,weight_79,79) &
       - int4(ri3_79,e(:,OP_DR ),h(:,OP_DZZ),temp79c,weight_79,79) &
       + int4(ri3_79,e(:,OP_DRZ),h(:,OP_DR ),temp79b,weight_79,79) &
       - int4(ri3_79,e(:,OP_DRR),h(:,OP_DZ ),temp79b,weight_79,79) &
       + int4(ri3_79,e(:,OP_DZ ),h(:,OP_DRR),temp79b,weight_79,79) &
       - int4(ri3_79,e(:,OP_DR ),h(:,OP_DRZ),temp79b,weight_79,79)

  if(itor.eq.1) then
     temp = temp &
          + int4(ri4_79,e(:,OP_DR),h(:,OP_DZ),temp79b,weight_79,79) &
          - int4(ri4_79,e(:,OP_DZ),h(:,OP_DR),temp79b,weight_79,79) &
          + int4(ri4_79,e(:,OP_DR),h(:,OP_DZ),temp79b,weight_79,79) &
          - int4(ri4_79,e(:,OP_DZ),h(:,OP_DR),temp79b,weight_79,79)
  endif

  v1chipsipsi = temp
  return
end function v1chipsipsi



! V1chibb
! =======
real function v1chibb(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  real :: temp

  if(itor.eq.0) then
     v1chibb = 0
     return
  endif

  temp = -2.* &
       (int5(ri6_79,e(:,OP_DZ),g(:,OP_1 ),f(:,OP_GS),h(:,OP_1),weight_79,79) &
       +int5(ri6_79,e(:,OP_DZ),g(:,OP_DZ),f(:,OP_DZ),h(:,OP_1),weight_79,79) &
       +int5(ri6_79,e(:,OP_DZ),g(:,OP_DR),f(:,OP_DR),h(:,OP_1),weight_79,79))
  
  v1chibb = temp
  return
end function v1chibb



! V1bsb2
! ======
subroutine v1bsb2(itri,i,j,ssterm, ddterm, rrterm)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm

  integer :: i2
  real :: temp

  if(itor.eq.0) return

  temp = 2.*int4(ri4_79,g79(:,OP_DZ,i),g79(:,OP_1,j),sb279(:,OP_1),weight_79,79)

  rrterm(1,1) = rrterm(1,1) + thimp*dt*dt*temp

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     i2 = isvaln(itri,i) + 6
     
     temp = 2.*int4(ri4_79,g79(:,OP_DZ,i),ph079(:,OP_1),sb279(:,OP_1),weight_79,79)

     r4(i2) = r4(i2) + thimp*dt*dt*temp
  end if

end subroutine v1bsb2


!===============================================================================
! V2 TERMS
!===============================================================================


! V2v
! ===========================
subroutine v2v(itri,i,j,ssterm, ddterm, rrterm)

  use basic
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm

  real :: temp

  if(idens.eq.0) then
     temp = int2(g25(:,OP_1,i),g25(:,OP_1,j),weight_25,25)
  else
     temp = int3(g79(:,OP_1,i),g79(:,OP_1,j),nt79(:,itri),weight_79,79)
  endif

  ssterm(2,2) = ssterm(2,2) + temp
  ddterm(2,2) = ddterm(2,2) + temp

end subroutine v2v


! V2vmu
! ===========================
subroutine v2vmu(itri,i,j,ssterm, ddterm, rrterm)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm

  integer :: i2
  real :: temp

  temp = -int2(g25(:,OP_DZ,i),g25(:,OP_DZ,j),weight_25,25) &
       -  int2(g25(:,OP_DR,i),g25(:,OP_DR,j),weight_25,25)
  if(itor.eq.1) then
     temp = temp - int3(ri3_79,g79(:,OP_1,i),g79(:,OP_1,j),weight_79,79)
  endif


  ssterm(2,2) = ssterm(2,2) -     thimp *dt*temp*amu
  ddterm(2,2) = ddterm(2,2) + (1.-thimp)*dt*temp*amu


  ! EQUILIBRIUM TERM
  if(linear.eq.1 .or. eqsubtract.eq.1) then
     i2 = isvaln(itri,i) + 6

     temp = -int2(g25(:,OP_DZ,i),vz025(:,OP_DZ),weight_25,25) &
          -  int2(g25(:,OP_DR,i),vz025(:,OP_DR),weight_25,25)
     if(itor.eq.1) then
        temp = temp - int3(ri3_79,g79(:,OP_1,i),vz025(:,OP_1),weight_79,79)
     endif

     r4(i2) = r4(i2) + dt*temp*amu
  endif

end subroutine v2vmu

! V2vun
! =====
subroutine v2vun(itri,i,j,ssterm, ddterm, rrterm)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm

  integer :: i2
  real :: temp

  if(idens.eq.0) then
     temp = int3(g79(:,OP_1,i),g79(:,OP_DR,j),ph179(:,OP_DZ),weight_79,79) &
          - int3(g79(:,OP_1,i),g79(:,OP_DZ,j),ph179(:,OP_DR),weight_79,79)
  else
     temp = int4(g79(:,OP_1,i),g79(:,OP_DR,j),ph179(:,OP_DZ),nt79(:,OP_1),weight_79,79) &
          - int4(g79(:,OP_1,i),g79(:,OP_DZ,j),ph179(:,OP_DR),nt79(:,OP_1),weight_79,79)
  endif

  ssterm(2,2) = ssterm(2,2) -      thimp *dt*temp
  ddterm(2,2) = ddterm(2,2) + (0.5-thimp)*dt*temp

  if(idens.eq.0) then
     temp = int3(g79(:,OP_1,i),vz179(:,OP_DR),g79(:,OP_DZ,j),weight_79,79) &
          - int3(g79(:,OP_1,i),vz179(:,OP_DZ),g79(:,OP_DR,j),weight_79,79)
  else
     temp = int4(g79(:,OP_1,i),vz179(:,OP_DR),g79(:,OP_DZ,j),nt79(:,OP_1),weight_79,79) &
          - int4(g79(:,OP_1,i),vz179(:,OP_DZ),g79(:,OP_DR,j),nt79(:,OP_1),weight_79,79)
  endif

  ssterm(2,1) = ssterm(2,1) -      thimp *dt*temp
  ddterm(2,1) = ddterm(2,1) + (0.5-thimp)*dt*temp
  
  if(linear.eq.1 .or. eqsubtract.eq.1) then
     if(idens.eq.0) then
        temp = int3(g79(:,OP_1,i),g79(:,OP_DR,j),ph079(:,OP_DZ),weight_79,79) &
             - int3(g79(:,OP_1,i),g79(:,OP_DZ,j),ph079(:,OP_DR),weight_79,79)
     else
        temp = int4(g79(:,OP_1,i),g79(:,OP_DR,j),ph079(:,OP_DZ),nt79(:,OP_1),weight_79,79) &
             - int4(g79(:,OP_1,i),g79(:,OP_DZ,j),ph079(:,OP_DR),nt79(:,OP_1),weight_79,79)
     endif

     ssterm(2,2) = ssterm(2,2) -     thimp *dt*temp
     ddterm(2,2) = ddterm(2,2) + (1.-thimp)*dt*temp

     if(idens.eq.0) then
        temp = int3(g79(:,OP_1,i),vz079(:,OP_DR),g79(:,OP_DZ,j),weight_79,79) &
             - int3(g79(:,OP_1,i),vz079(:,OP_DZ),g79(:,OP_DR,j),weight_79,79)
     else
        temp = int4(g79(:,OP_1,i),vz079(:,OP_DR),g79(:,OP_DZ,j),nt79(:,OP_1),weight_79,79) &
             - int4(g79(:,OP_1,i),vz079(:,OP_DZ),g79(:,OP_DR,j),nt79(:,OP_1),weight_79,79)
     endif

     ssterm(2,1) = ssterm(2,1) -     thimp *dt*temp
     ddterm(2,1) = ddterm(2,1) + (1.-thimp)*dt*temp

     i2 = isvaln(itri,i) + 6

     if(idens.eq.1) then
        temp = int4(g79(:,OP_1,i),vz079(:,OP_DR),ph079(:,OP_DZ),n179(:,OP_1),weight_79,79) &
             - int4(g79(:,OP_1,i),vz079(:,OP_DZ),ph079(:,OP_DR),n179(:,OP_1),weight_79,79)

        r4(i2) = r4(i2) + dt*temp
     endif

     ! EQUILIBRIUM TERM
     if(idens.eq.0) then
        temp = int3(g79(:,OP_1,i),vz079(:,OP_DR),ph079(:,OP_DZ),weight_79,79) &
             - int3(g79(:,OP_1,i),vz079(:,OP_DZ),ph079(:,OP_DR),weight_79,79)
     else
        temp = int4(g79(:,OP_1,i),vz079(:,OP_DR),ph079(:,OP_DZ),n079(:,OP_1),weight_79,79) &
             - int4(g79(:,OP_1,i),vz079(:,OP_DZ),ph079(:,OP_DR),n079(:,OP_1),weight_79,79)
     endif

     r4(i2) = r4(i2) + dt*temp
  endif

end subroutine v2vun


! V2psib
! ===========================
subroutine v2psib(itri,i,j,ssterm, ddterm, rrterm)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm

  integer :: i2
  real :: temp

  ! Linearize in psi
  ! ~~~~~~~~~~~~~~~~
  temp = 0.5* &
       (int4(ri2_79,g79(:,OP_1,i),g79(:,OP_DR,j),bz179(:,OP_DZ),weight_79,79) &
       -int4(ri2_79,g79(:,OP_1,i),g79(:,OP_DZ,j),bz179(:,OP_DR),weight_79,79))

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     temp = temp + &
          (int4(ri2_79,g79(:,OP_1,i),g79(:,OP_DR,j),bz079(:,OP_DZ),weight_79,79) &
          -int4(ri2_79,g79(:,OP_1,i),g79(:,OP_DZ,j),bz079(:,OP_DR),weight_79,79))
  endif

  rrterm(2,1) = rrterm(2,1) + dt*temp

  ! Linearize in bz
  ! ~~~~~~~~~~~~~~~
  temp = 0.5* &
       (int4(ri2_79,g79(:,OP_1,i),ps179(:,OP_DR),g79(:,OP_DZ,j),weight_79,79) &
       -int4(ri2_79,g79(:,OP_1,i),ps179(:,OP_DZ),g79(:,OP_DR,j),weight_79,79))

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     temp = temp + &
          (int4(ri2_79,g79(:,OP_1,i),ps079(:,OP_DR),g79(:,OP_DZ,j),weight_79,79) &
          -int4(ri2_79,g79(:,OP_1,i),ps079(:,OP_DZ),g79(:,OP_DR,j),weight_79,79))
  endif

  rrterm(2,2) = rrterm(2,2) + dt*temp

  ! EQUILIBRIUM TERM
  if(linear.eq.1 .or. eqsubtract.eq.1) then
     i2 = isvaln(itri,i) + 6

     temp = int4(ri2_79,g79(:,OP_1,i),ps079(:,OP_DR),bz079(:,OP_DZ),weight_79,79) &
          - int4(ri2_79,g79(:,OP_1,i),ps079(:,OP_DZ),bz079(:,OP_DR),weight_79,79)

     r4(i2) = r4(i2) + dt*temp
  endif  

end subroutine v2psib


! V2psisb1
! ========
subroutine v2psisb1(itri,i,j,ssterm, ddterm, rrterm)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm

  integer :: i2
  real :: temp

  temp = int4(ri2_79,g79(:,OP_1,i),sb179(:,OP_DR),g79(:,OP_DZ,j),weight_79,79) &
       - int4(ri2_79,g79(:,OP_1,i),sb179(:,OP_DZ),g79(:,OP_DR,j),weight_79,79)

  rrterm(2,2) = rrterm(2,2) + thimp*dt*dt*temp

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     i2 = isvaln(itri,i) + 6

     temp = int4(ri2_79,g79(:,OP_1,i),sb179(:,OP_DR),bz079(:,OP_DZ),weight_79,79) &
          - int4(ri2_79,g79(:,OP_1,i),sb179(:,OP_DZ),bz079(:,OP_DR),weight_79,79)

     r4(i2) = r4(i2) + thimp*dt*dt*temp
  endif  

end subroutine v2psisb1


! V2psisb2
! ========
subroutine v2psisb2(itri,i,j,ssterm, ddterm, rrterm)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm

  integer :: i2
  real :: temp

  temp = int4(ri2_79,g79(:,OP_1,i),g79(:,OP_DR,j),sb279(:,OP_DZ),weight_79,79) &
       - int4(ri2_79,g79(:,OP_1,i),g79(:,OP_DZ,j),sb279(:,OP_DR),weight_79,79)

  rrterm(2,1) = rrterm(2,1) + thimp*dt*dt*temp

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     i2 = isvaln(itri,i) + 6

     temp = int4(ri2_79,g79(:,OP_1,i),ps079(:,OP_DR),sb279(:,OP_DZ),weight_79,79) &
          - int4(ri2_79,g79(:,OP_1,i),ps079(:,OP_DZ),sb279(:,OP_DR),weight_79,79)

     r4(i2) = r4(i2) + thimp*dt*dt*temp
  endif  

end subroutine v2psisb2


! V2vpsipsi
! =========
subroutine v2vpsipsi(itri,i,j,ssterm, ddterm, rrterm)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm

  integer :: i2
  real :: temp


  ! Linearize in v
  ! ~~~~~~~~~~~~~~

  ! [nu,psi(2)] + nu*psi(2)_z/r
  temp79a = g79(:,OP_DZ,i)*pst79(:,OP_DR) - g79(:,OP_DR,i)*pst79(:,OP_DZ)
  if(itor.eq.1) temp79a = temp79a + ri_79*g79(:,OP_1,i)*pst79(:,OP_DZ)

  ! [v,psi(1)] + v*psi(1)_z/r
  temp79b = g79(:,OP_DZ,j)*pst79(:,OP_DR) - g79(:,OP_DR,j)*pst79(:,OP_DZ)
  if(itor.eq.1) temp79b = temp79b + ri_79*g79(:,OP_1,j)*pst79(:,OP_DZ)

  temp = -int3(ri4_79, temp79a, temp79b, weight_79, 79)

  ssterm(2,2) = ssterm(2,2) - thimp*    thimp *dt*dt*temp
  ddterm(2,2) = ddterm(2,2) + thimp*(1.-thimp)*dt*dt*temp


  if(linear.eq.1 .or. eqsubtract.eq.1) then

     ! Linearize in psi(1)
     ! ~~~~~~~~~~~~~~~~~~~
     
     ! [nu,psi(2)] + nu*psi(2)_z/r
     temp79a = (g79(:,OP_DZ,i)*ps179(:,OP_DR) - g79(:,OP_DR,i)*ps179(:,OP_DZ))/2. &
          + g79(:,OP_DZ,i)*ps079(:,OP_DR) - g79(:,OP_DR,i)*ps079(:,OP_DZ)
     if(itor.eq.1) temp79a = temp79a + ri_79*g79(:,OP_1,i)*(ps079(:,OP_DZ)+ps179(:,OP_DZ)/2.)

     ! [v,psi(1)] + v*psi(1)_z/r
     temp79b = vz079(:,OP_DZ)*g79(:,OP_DR,j) - vz079(:,OP_DR)*g79(:,OP_DZ,j)
     if(itor.eq.1) temp79b = temp79b + ri_79*vz079(:,OP_1)*g79(:,OP_DZ,j)

     temp = -int3(ri4_79, temp79a, temp79b, weight_79, 79)
  
     rrterm(2,1) = rrterm(2,1) + thimp*dt*dt*temp
  

     ! Linearize in psi(2)
     ! ~~~~~~~~~~~~~~~~~~~
     
     ! [nu,psi(2)] + nu*psi(2)_z/r
     temp79a = g79(:,OP_DZ,i)*g79(:,OP_DR,j) - g79(:,OP_DR,i)*g79(:,OP_DZ,j)
     if(itor.eq.1) temp79a = temp79a + ri_79*g79(:,OP_1,i)*g79(:,OP_DZ,j)

     ! [v,psi(1)] + v*psi(1)_z/r
     temp79b = (vz079(:,OP_DZ)*ps179(:,OP_DR) - vz079(:,OP_DR)*ps179(:,OP_DZ))/2. &
          + g79(:,OP_DZ,i)*ps079(:,OP_DR) - g79(:,OP_DR,i)*ps079(:,OP_DZ)
     if(itor.eq.1) temp79b = temp79b + ri_79*vz079(:,OP_1)*(ps079(:,OP_DZ)+ps179(:,OP_DZ)/2.)

     temp = -int3(ri4_79, temp79a, temp79b, weight_79, 79)
  
     rrterm(2,1) = rrterm(2,1) + thimp*dt*dt*temp

     ! EQUILIBRIUM TERM
     i2 = isvaln(itri, i) + 6
        
     ! [nu,psi(2)] + nu*psi(2)_z/r
     temp79a = g79(:,OP_DZ,i)*ps079(:,OP_DR) - g79(:,OP_DR,i)*ps079(:,OP_DZ)
     if(itor.eq.1) temp79a = temp79a + ri_79*g79(:,OP_1,i)*ps079(:,OP_DZ)
     
     ! [v,psi(1)] + v*psi(1)_z/r
     temp79b = vz079(:,OP_DZ)*ps079(:,OP_DR) - vz079(:,OP_DR)*ps079(:,OP_DZ)
     if(itor.eq.1) temp79b = temp79b + ri_79*vz079(:,OP_1)*ps079(:,OP_DZ)
     
     temp = -int3(ri4_79, temp79a, temp79b, weight_79, 79)
     
     r4(i2) = r4(i2) + thimp*dt*dt*temp
  endif

end subroutine v2vpsipsi


! V2upsib
! =========
subroutine v2upsib(itri,i,j,ssterm, ddterm, rrterm)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm

  integer :: i2
  real :: temp


  ! Linearize in v
  ! ~~~~~~~~~~~~~~

  ! [nu,bz] + nu*bz_z/r
  temp79a = g79(:,OP_DZ,i)*bzt79(:,OP_DR) - g79(:,OP_DR,i)*bzt79(:,OP_DZ)
  if(itor.eq.1) temp79a = temp79a + ri_79*g79(:,OP_1,i)*bzt79(:,OP_DZ)

  ! [nu,psi] + nu*psi_z/r
  temp79b = g79(:,OP_DZ,i)*pst79(:,OP_DR) - g79(:,OP_DR,i)*pst79(:,OP_DZ)
  if(itor.eq.1) temp79b = temp79b + ri_79*g79(:,OP_1,i)*pst79(:,OP_DZ)


  temp = int4(ri_79, g79(:,OP_DZ,j), pst79(:,OP_DR), temp79a, weight_79, 79)  &
       - int4(ri_79, g79(:,OP_DR,j), pst79(:,OP_DZ), temp79a, weight_79, 79)  &
       + int4(ri3_79, g79(:,OP_DR,j), bzt79(:,OP_DZ), temp79b, weight_79, 79) &
       - int4(ri3_79, g79(:,OP_DZ,j), bzt79(:,OP_DR), temp79b, weight_79, 79)

  ssterm(2,1) = ssterm(2,1) - thimp*    thimp *dt*dt*temp
  ddterm(2,1) = ddterm(2,1) + thimp*(1.-thimp)*dt*dt*temp


  if(linear.eq.1 .or. eqsubtract.eq.1) then

     ! Linearize in psi
     ! ~~~~~~~~~~~~~~~~
     tm79 = bz079 + bz179/2.

     ! [nu,bz] + nu*bz_z/r
     temp79a = g79(:,OP_DZ,i)*tm79(:,OP_DR) - g79(:,OP_DR,i)*tm79(:,OP_DZ)
     if(itor.eq.1) temp79a = temp79a + ri_79*g79(:,OP_1,i)*tm79(:,OP_DZ)

     ! [nu,psi] + nu*psi_z/r
     temp79b = g79(:,OP_DZ,i)*g79(:,OP_DR,j) - g79(:,OP_DR,i)*g79(:,OP_DZ,j)
     if(itor.eq.1) temp79b = temp79b + ri_79*g79(:,OP_1,i)*g79(:,OP_DZ,j)

     temp = int4(ri_79, vz079(:,OP_DZ), g79(:,OP_DR,j), temp79a, weight_79, 79) &
          - int4(ri_79, vz079(:,OP_DR), g79(:,OP_DZ,j), temp79a, weight_79, 79) &
          + int4(ri3_79, vz079(:,OP_DR), tm79(:,OP_DZ), temp79b, weight_79, 79) &
          - int4(ri3_79, vz079(:,OP_DZ), tm79(:,OP_DR), temp79b, weight_79, 79)

     rrterm(2,1) = rrterm(2,1) + thimp**dt*dt*temp


     ! Linearize in bz
     ! ~~~~~~~~~~~~~~~
     tm79 = ps079 + ps179/2.

     ! [nu,bz] + nu*bz_z/r
     temp79a = g79(:,OP_DZ,i)*g79(:,OP_DR,j) - g79(:,OP_DR,i)*g79(:,OP_DZ,j)
     if(itor.eq.1) temp79a = temp79a + ri_79*g79(:,OP_1,i)*g79(:,OP_DZ,j)

     ! [nu,psi] + nu*psi_z/r
     temp79b = g79(:,OP_DZ,i)*tm79(:,OP_DR) - g79(:,OP_DR,i)*tm79(:,OP_DZ)
     if(itor.eq.1) temp79b = temp79b + ri_79*g79(:,OP_1,i)*tm79(:,OP_DZ)

     temp = int4(ri_79, vz079(:,OP_DZ), tm79(:,OP_DR), temp79a, weight_79, 79)   &
          - int4(ri_79, vz079(:,OP_DR), tm79(:,OP_DZ), temp79a, weight_79, 79)   &
          + int4(ri3_79, vz079(:,OP_DR), g79(:,OP_DZ,j), temp79b, weight_79, 79) &
          - int4(ri3_79, vz079(:,OP_DZ), g79(:,OP_DR,j), temp79b, weight_79, 79)

     rrterm(2,2) = rrterm(2,2) + thimp**dt*dt*temp


     ! EQUILIBRIUM TERM

     i2 = isvaln(itri,i)

     ! [nu,bz] + nu*bz_z/r
     temp79a = g79(:,OP_DZ,i)*bz079(:,OP_DR) - g79(:,OP_DR,i)*bz079(:,OP_DZ)
     if(itor.eq.1) temp79a = temp79a + ri_79*g79(:,OP_1,i)*bz079(:,OP_DZ)

     ! [nu,psi] + nu*psi_z/r
     temp79b = g79(:,OP_DZ,i)*ps079(:,OP_DR) - g79(:,OP_DR,i)*ps079(:,OP_DZ)
     if(itor.eq.1) temp79b = temp79b + ri_79*g79(:,OP_1,i)*ps079(:,OP_DZ)

     temp = int4(ri_79, vz079(:,OP_DZ), ps079(:,OP_DR), temp79a, weight_79, 79)  &
          - int4(ri_79, vz079(:,OP_DR), ps079(:,OP_DZ), temp79a, weight_79, 79)  &
          + int4(ri3_79, vz079(:,OP_DR), bz079(:,OP_DZ), temp79b, weight_79, 79) &
          - int4(ri3_79, vz079(:,OP_DZ), bz079(:,OP_DR), temp79b, weight_79, 79)

     r4(i2) = r4(i2) + dt*temp
  endif

end subroutine v2upsib


! v2upsisb2
! =========
real function v2upsisb2(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp

  temp = int5(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DR),h(:,OP_DZ),weight_79,79) &
       - int5(ri_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DZ),h(:,OP_DZ),weight_79,79) &
       - int5(ri_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DR),h(:,OP_DR),weight_79,79) &
       + int5(ri_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ),h(:,OP_DR),weight_79,79)
  if(itor.eq.1) then
     temp = temp &
          + int5(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_DZ),weight_79,79) &
          - int5(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_DZ),weight_79,79)
  endif

  v2upsisb2 = temp
  return
end function v2upsisb2


! v2chipsib
! =========
real function v2chipsib(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp


  temp79a = h(:,OP_1 )*f(:,OP_GS) &
       +    h(:,OP_DZ)*f(:,OP_DZ) + h(:,OP_DR)*f(:,OP_DR)
  temp79b = e(:,OP_DZ)*g(:,OP_DR) - e(:,OP_DR)*g(:,OP_DZ)

  temp = int5(ri2_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_DZ),weight_79,79) &
       + int5(ri2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DR),h(:,OP_DZ),weight_79,79) &
       - int5(ri2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_DR),weight_79,79) &
       - int5(ri2_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DR),h(:,OP_DR),weight_79,79) &
       + int3(ri4_79,temp79a,temp79b,weight_79,79)
  if(itor.eq.1) then
     temp = temp &
          - int5(ri3_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_DZ),weight_79,79) &
          - int5(ri3_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_DZ),weight_79,79) &
          + int4(ri5_79,temp79a,e(:,OP_1),g(:,OP_DZ),weight_79,79)
  endif

  v2chipsib = temp
  return
end function v2chipsib


! v2vchin
! =======
real function v2vchin(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp

  if(idens.eq.0) then
     temp =-int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
          - int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),weight_79,79)
  else
     temp =-int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          - int5(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1),weight_79,79)
  endif

  v2vchin = temp
  return
end function v2vchin



!===============================================================================
! V3 TERMS
!===============================================================================

! V3chin
! ======
real function v3chin(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  
  real :: temp

  if(idens.eq.0) then
     temp = - int2(e(:,OP_DR),f(:,OP_DR),weight_79,79) &
          -   int2(e(:,OP_DZ),f(:,OP_DZ),weight_79,79)
  else
     temp = - int3(e(:,OP_DR),f(:,OP_DR),g(:,OP_1),weight_79,79) &
          -   int3(e(:,OP_DZ),f(:,OP_DZ),g(:,OP_1),weight_79,79)
  endif

  v3chin = temp
  return
end function v3chin



! V3chimu
! =======
real function v3chimu(e,f)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real :: temp

  temp = int2(e(:,OP_LP),f(:,OP_LP),weight_79,79)
  if(itor.eq.1) then
     temp = temp + 2.*int3(ri_79,e(:,OP_DR),f(:,OP_LP),weight_79,79)
  endif

  v3chimu = temp*amu
  return
end function v3chimu


! V3umu
! =====
subroutine v3umu(itri,i,j,ssterm, ddterm, rrterm)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm

  integer :: i3
  real :: temp

  if(itor.eq.0) return

  temp = -4.*int3(ri_79,g79(:,OP_DZ,i),g79(:,OP_DR,j),weight_79,79)

  ssterm(3,1) = ssterm(3,1) -     thimp *dt*temp*amu
  ddterm(3,1) = ddterm(3,1) + (1.-thimp)*dt*temp*amu

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     ! EQUILIBRIUM TERM
     i3 = isvaln(itri,i) + 12
     
     temp = -4.*int3(ri_79,g79(:,OP_DZ,i),ph079(:,OP_DR),weight_79,79)
  
     r4(i3) = r4(i3) + dt*temp*amu
  endif

end subroutine v3umu


! V3un
! ====
subroutine v3un(itri,i,j,ssterm, ddterm, rrterm)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm

  integer :: i3
  real :: temp

  if(idens.eq.0) return

  temp = int4(r_79,g79(:,OP_1,i),g79(:,OP_DR,j),nt79(:,OP_DZ),weight_79,79) &
       - int4(r_79,g79(:,OP_1,i),g79(:,OP_DZ,j),nt79(:,OP_DR),weight_79,79) 

  ssterm(3,1) = ssterm(3,1) + temp
  ddterm(3,1) = ddterm(3,1) + temp

end subroutine v3un


! V3p
! ===
real function v3p(e,f)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f
  real :: temp

  temp = int2(e(:,OP_DZ),f(:,OP_DZ),weight_79,79) &
       + int2(e(:,OP_DR),f(:,OP_DR),weight_79,79)

  if(itor.eq.1) then
     temp = temp + 2.*int3(ri_79,e(:,OP_1),f(:,OP_DR),weight_79,79)
  endif

  v3p = temp
  return
end function v3p



! V3up
! ====
real function v3up(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g

  real :: temp

  temp = int4(r_79,e(:,OP_LP),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
       - int4(r_79,e(:,OP_LP),f(:,OP_DZ),g(:,OP_DR),weight_79,79)

  if(itor.eq.1) then
     temp = temp + 2.* &
          (int3(e(:,OP_DR),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
          -int3(e(:,OP_DR),f(:,OP_DZ),g(:,OP_DR),weight_79,79) &
          -2.*gam* &
          (int3(e(:,OP_LP),f(:,OP_DZ),g(:,OP_1),weight_79,79) &
          +int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_1),weight_79,79)))
  endif

  v3up = temp
  
  return
end function v3up


! V3chip
! ======
real function v3chip(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp

  temp = int3(e(:,OP_LP),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
       + int3(e(:,OP_LP),f(:,OP_DR),g(:,OP_DR),weight_79,79) &
       + gam*int3(e(:,OP_LP),f(:,OP_LP),g(:,OP_1),weight_79,79)
  if(itor.eq.1) then
     temp = temp + 2.* &
          (int4(ri_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
          +int4(ri_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DR),weight_79,79) &
          +gam*int4(ri_79,e(:,OP_DR),f(:,OP_LP),g(:,OP_1),weight_79,79))
  endif

  v3chip = temp

  return
end function v3chip



! V3psipsi
! ========
real function v3psipsi(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g

  real :: temp

  temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DZ),weight_79,79) &
       + int4(ri2_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DR),weight_79,79)

  if(itor.eq.1) then
     temp = temp + 2.*int4(ri3_79,e(:,OP_1),f(:,OP_GS),g(:,OP_DR),weight_79,79)
  endif

  v3psipsi = temp

  return
end function v3psipsi


! V3bb
! ====
real function v3bb(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g

  real :: temp

  temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_1),g(:,OP_DZ),weight_79,79) &
       + int4(ri2_79,e(:,OP_DR),f(:,OP_1),g(:,OP_DR),weight_79,79)

  if(itor.eq.1) then
     temp = temp + 2.* &
          int4(ri3_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),weight_79,79)
  endif

  v3bb = temp
  return
end function v3bb


! V3psisb1
! ========
real function v3psisb1(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  real :: temp

  temp = int4(ri2_79,e(:,OP_DZ),f(:,OP_GS),g(:,OP_DZ),weight_79,79) &
       + int4(ri2_79,e(:,OP_DR),f(:,OP_GS),g(:,OP_DR),weight_79,79)

  if(itor.eq.1) then
     temp = temp + 2.* &
          int4(ri3_79,e(:,OP_1),f(:,OP_GS),g(:,OP_DR),weight_79,79)
  endif

  v3psisb1 = temp

end function v3psisb1


! V3bsb2
! ======
real function v3bsb2(e,f,g)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g

  real :: temp

  temp = -int4(ri2_79,e(:,OP_LP),f(:,OP_1),g(:,OP_1),weight_79,79)

  if(itor.eq.1) then
     temp = temp + 4.*int4(ri4_79,e(:,OP_1),f(:,OP_1),g(:,OP_1),weight_79,79)
  endif

  v3bsb2 = temp
  return
end function v3bsb2


! V3upsipsi
! =========
real function v3upsipsi(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  real :: temp

  ! Linearize in u
  ! ~~~~~~~~~~~~~~

  ! [b,c]
  temp79a = f(:,OP_DZ )*g(:,OP_DR ) - f(:,OP_DR )*g(:,OP_DZ)

  ! [b,c],r
  temp79b = f(:,OP_DRZ)*g(:,OP_DR ) - f(:,OP_DRR)*g(:,OP_DZ) &
       +    f(:,OP_DZ )*g(:,OP_DRR) - f(:,OP_DR )*g(:,OP_DRZ)

  ! [b,c],z
  temp79c = f(:,OP_DZZ)*g(:,OP_DR ) - f(:,OP_DRZ)*g(:,OP_DZ) &
       +    f(:,OP_DZ )*g(:,OP_DRZ) - f(:,OP_DR )*g(:,OP_DZZ)

  temp =-int4(ri_79,temp79b,e(:,OP_DRZ),h(:,OP_DRZ),weight_79,79) & 
       - int4(ri_79,temp79b,e(:,OP_DRR),h(:,OP_DRR),weight_79,79) &
       - int4(ri_79,temp79c,e(:,OP_DZZ),h(:,OP_DZZ),weight_79,79) & 
       - int4(ri_79,temp79c,e(:,OP_DRZ),h(:,OP_DRZ),weight_79,79) &
       + int4(ri_79,e(:,OP_DR),temp79b,h(:,OP_GS),weight_79,79) &
       + int4(ri_79,e(:,OP_DZ),temp79c,h(:,OP_GS),weight_79,79)

  if(itor.eq.1) then
     temp = temp &
          + int4(ri2_79,e(:,OP_DR ),temp79a,h(:,OP_GS ),weight_79,79) &
          - int4(ri2_79,e(:,OP_DRZ),temp79a,h(:,OP_DZ ),weight_79,79) &
          - int4(ri2_79,e(:,OP_DRR),temp79a,h(:,OP_DR ),weight_79,79) &
          - int4(ri2_79,e(:,OP_DZ ),temp79a,h(:,OP_DRZ),weight_79,79) &
          - int4(ri2_79,e(:,OP_DR ),temp79a,h(:,OP_DRR),weight_79,79) &
          - 2.*int4(ri2_79,e(:,OP_1  ),temp79b,h(:,OP_DRR),weight_79,79) &
          - 2.*int4(ri2_79,e(:,OP_1  ),temp79c,h(:,OP_DRZ),weight_79,79) &
          - 2.*int4(ri2_79,e(:,OP_DR ),temp79b,h(:,OP_DR ),weight_79,79) &
          - 2.*int4(ri2_79,e(:,OP_DZ ),temp79c,h(:,OP_DR ),weight_79,79) &
          + 2.*int4(ri2_79,e(:,OP_1  ),temp79b,h(:,OP_GS ),weight_79,79) &
          + 2.*int4(ri3_79,e(:,OP_1  ),temp79b,h(:,OP_DR ),weight_79,79) &
          - 2.*int4(ri3_79,e(:,OP_DR ),temp79a,h(:,OP_DR ),weight_79,79) &
          - 2.*int4(ri3_79,e(:,OP_1  ),temp79a,h(:,OP_DRR),weight_79,79) &
          + 2.*int4(ri3_79,e(:,OP_1  ),temp79a,h(:,OP_GS ),weight_79,79) &
          + 2.*int4(ri4_79,e(:,OP_1  ),temp79a,h(:,OP_DR ),weight_79,79)
  endif

  v3upsipsi = temp
  return
end function v3upsipsi


! V3ubb
! =====
real function v3ubb(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h

  real :: temp

  temp = int5(ri3_79,e(:,OP_LP),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
       - int5(ri3_79,e(:,OP_LP),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79)

  if(itor.eq.1) then
     temp = temp - 4.* &
          (int5(ri5_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          -int5(ri5_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79))
  endif

  v3ubb = temp
  return
end function v3ubb


! v3vpsib
! =======
real function v3vpsib(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  real :: temp

  temp = int5(ri4_79,e(:,OP_LP),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
       - int5(ri4_79,e(:,OP_LP),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79)
  if(itor.eq.1) then
     temp = temp &
          - int5(ri5_79,e(:,OP_LP),f(:,OP_1),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          -4.* &
          (int5(ri6_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          -int5(ri6_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
          -int5(ri7_79,e(:,OP_1),f(:,OP_1 ),g(:,OP_DZ),h(:,OP_1),weight_79,79))
  endif

  v3vpsib = temp
  return
end function v3vpsib


! V3chipsipsi
! ===========
real function v3chipsipsi(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  real :: temp

  ! <f,g>
  temp79a = f(:,OP_DR )*g(:,OP_DR ) + f(:,OP_DZ )*g(:,OP_DZ )

  ! <f,g>,r
  temp79b = f(:,OP_DRR)*g(:,OP_DR ) + f(:,OP_DZ )*g(:,OP_DRZ) &
       +    f(:,OP_DRR)*g(:,OP_DR ) + f(:,OP_DZ )*g(:,OP_DRZ)

  ! <f,g>,z
  temp79c = f(:,OP_DRZ)*g(:,OP_DR ) + f(:,OP_DZ )*g(:,OP_DZZ) &
       +    f(:,OP_DRZ)*g(:,OP_DR ) + f(:,OP_DZ )*g(:,OP_DZZ)

  temp = int4(ri2_79,temp79b,e(:,OP_DRZ),h(:,OP_DZ ),weight_79,79) &
       + int4(ri2_79,temp79b,e(:,OP_DRR),h(:,OP_DR ),weight_79,79) &
       + int4(ri2_79,temp79b,e(:,OP_DZ ),h(:,OP_DRZ),weight_79,79) &
       + int4(ri2_79,temp79b,e(:,OP_DR ),h(:,OP_DRR),weight_79,79) &
       + int4(ri2_79,temp79c,e(:,OP_DZZ),h(:,OP_DZ ),weight_79,79) &
       + int4(ri2_79,temp79c,e(:,OP_DRZ),h(:,OP_DR ),weight_79,79) &
       + int4(ri2_79,temp79c,e(:,OP_DZ ),h(:,OP_DZZ),weight_79,79) &
       + int4(ri2_79,temp79c,e(:,OP_DR ),h(:,OP_DRZ),weight_79,79) &
       - int4(ri2_79,temp79c,e(:,OP_DZ ),h(:,OP_GS ),weight_79,79) &
       - int4(ri2_79,temp79b,e(:,OP_DR ),h(:,OP_GS ),weight_79,79)
  if(itor.eq.1) then
     temp = temp + 2.* &
          (int4(ri3_79,temp79c,e(:,OP_1  ),h(:,OP_DRZ),weight_79,79) &
          +int4(ri3_79,temp79b,e(:,OP_1  ),h(:,OP_DRR),weight_79,79) &
          +int4(ri3_79,temp79c,e(:,OP_DRZ),h(:,OP_1  ),weight_79,79) &
          +int4(ri3_79,temp79b,e(:,OP_DRR),h(:,OP_1  ),weight_79,79) &
          -int4(ri3_79,temp79b,e(:,OP_1  ),h(:,OP_GS ),weight_79,79) &
          -int4(ri4_79,temp79b,e(:,OP_1  ),h(:,OP_DR ),weight_79,79))
  endif

  v3chipsipsi = temp

  return
end function v3chipsipsi


! V3chibb
! =======
real function v3chibb(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  real :: temp

  temp = int5(ri4_79,e(:,OP_LP),f(:,OP_GS),g(:,OP_1 ),h(:,OP_1),weight_79,79) &
       + int5(ri4_79,e(:,OP_LP),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
       + int5(ri4_79,e(:,OP_LP),f(:,OP_DR),g(:,OP_DR),h(:,OP_1),weight_79,79)
  if(itor.eq.1) then
     temp = temp - 4.*&
          (int5(ri6_79,e(:,OP_1),f(:,OP_GS),g(:,OP_1 ),h(:,OP_1),weight_79,79) &
          +int5(ri6_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          +int5(ri6_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),h(:,OP_1),weight_79,79))
  endif

  v3chibb = temp

  return
end function v3chibb


! V3uun
! =====
real function v3uun(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  real :: temp

  if(idens.eq.0) then
     temp = int4(r2_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DRZ),weight_79,79) &
          - int4(r2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DRR),weight_79,79) &
          + int4(r2_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DRZ),weight_79,79) &
          - int4(r2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DZZ),weight_79,79)
     if(itor.eq.1) then
        temp = temp &
             + int4(r_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ ),weight_79,79) &
             - int4(r_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DZ ),weight_79,79) &
             - int4(r_79,e(:,OP_DZ),g(:,OP_DR),f(:,OP_DZ ),weight_79,79) &
             + int4(r_79,e(:,OP_DR),g(:,OP_DZ),f(:,OP_DZ ),weight_79,79) &
             + int4(r_79,e(:,OP_1 ),f(:,OP_DZ),g(:,OP_DRZ),weight_79,79) &
             - int4(r_79,e(:,OP_1 ),f(:,OP_DR),g(:,OP_DZZ),weight_79,79)
     endif
  else
     temp = int5(r2_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DRZ),h(:,OP_1),weight_79,79) &
          - int5(r2_79,e(:,OP_DZ),f(:,OP_DZ),g(:,OP_DRR),h(:,OP_1),weight_79,79) &
          + int5(r2_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DRZ),h(:,OP_1),weight_79,79) &
          - int5(r2_79,e(:,OP_DR),f(:,OP_DR),g(:,OP_DZZ),h(:,OP_1),weight_79,79)
     if(itor.eq.1) then
        temp = temp &
             + int5(r_79,e(:,OP_DZ),f(:,OP_DR),g(:,OP_DZ ),h(:,OP_1 ),weight_79,79) &
             - int5(r_79,e(:,OP_DR),f(:,OP_DZ),g(:,OP_DZ ),h(:,OP_1 ),weight_79,79) &
             - int5(r_79,e(:,OP_DZ),g(:,OP_DR),f(:,OP_DZ ),h(:,OP_1 ),weight_79,79) &
             + int5(r_79,e(:,OP_DR),g(:,OP_DZ),f(:,OP_DZ ),h(:,OP_1 ),weight_79,79) &
             + int5(r_79,e(:,OP_1 ),f(:,OP_DZ),g(:,OP_DRZ),h(:,OP_1 ),weight_79,79) &
             - int5(r_79,e(:,OP_1 ),f(:,OP_DR),g(:,OP_DZZ),h(:,OP_1 ),weight_79,79) &
             - int5(r_79,e(:,OP_1 ),f(:,OP_DZ),g(:,OP_DZ ),h(:,OP_DR),weight_79,79) &
             + int5(r_79,e(:,OP_1 ),f(:,OP_DR),g(:,OP_DZ ),h(:,OP_DZ),weight_79,79)
     endif
  endif

  v3uun = temp
  return
end function v3uun


! V3uchin
! =======
real function v3uchin(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  real :: temp

  if(idens.eq.0) then
     temp = 2.* &
          (int3(e(:,OP_1),f(:,OP_DR ),g(:,OP_DRZ),weight_79,79) &
          -int3(e(:,OP_1),f(:,OP_DZZ),g(:,OP_DZ ),weight_79,79) &
          -int3(e(:,OP_1),f(:,OP_DRZ),g(:,OP_DR ),weight_79,79) &
          -int3(e(:,OP_1),f(:,OP_DZ ),g(:,OP_DRR),weight_79,79))
     if(itor.eq.1) then
        temp = temp &
             -2.*int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
     endif
  else
     temp = 2.* &
          (int4(e(:,OP_1),f(:,OP_DR ),g(:,OP_DRZ),h(:,OP_1),weight_79,79) &
          -int4(e(:,OP_1),f(:,OP_DZZ),g(:,OP_DZ ),h(:,OP_1),weight_79,79) &
          -int4(e(:,OP_1),f(:,OP_DRZ),g(:,OP_DR ),h(:,OP_1),weight_79,79) &
          -int4(e(:,OP_1),f(:,OP_DZ ),g(:,OP_DRR),h(:,OP_1),weight_79,79))
     if(itor.eq.1) then
        temp = temp &
             -2.*int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79)
     endif
  endif

  v3uchin = temp
  return
end function v3uchin



! V3chichin
! =========
real function v3chichin(e,f,g,h)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e, f, g, h

  real :: temp

  if(idens.eq.0) then
     temp = int3(e(:,OP_DZ),f(:,OP_DZZ),g(:,OP_DZ ),weight_79,79) &
          + int3(e(:,OP_DZ),f(:,OP_DRZ),g(:,OP_DR ),weight_79,79) &
          + int3(e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_DZZ),weight_79,79) &
          + int3(e(:,OP_DZ),f(:,OP_DR ),g(:,OP_DRZ),weight_79,79) &
          + int3(e(:,OP_DR),f(:,OP_DRZ),g(:,OP_DZ ),weight_79,79) &
          + int3(e(:,OP_DR),f(:,OP_DRR),g(:,OP_DR ),weight_79,79) &
          + int3(e(:,OP_DR),f(:,OP_DZ ),g(:,OP_DRZ),weight_79,79) &
          + int3(e(:,OP_DR),f(:,OP_DR ),g(:,OP_DRR),weight_79,79)
     if(itor.eq.1) then
        temp = temp &
          - int4(ri_79,e(:,OP_DR),f(:,OP_DZ ),g(:,OP_DZ ),weight_79,79) &
          - int4(ri_79,e(:,OP_DR),f(:,OP_DR ),g(:,OP_DR ),weight_79,79) &
          + int4(ri_79,e(:,OP_1 ),f(:,OP_DRZ),g(:,OP_DZ ),weight_79,79) &
          + int4(ri_79,e(:,OP_1 ),f(:,OP_DRR),g(:,OP_DR ),weight_79,79) &
          + int4(ri_79,e(:,OP_1 ),f(:,OP_DZ ),g(:,OP_DRZ),weight_79,79) &
          + int4(ri_79,e(:,OP_1 ),f(:,OP_DR ),g(:,OP_DRR),weight_79,79)
     endif
  else
     temp = int4(e(:,OP_DZ),f(:,OP_DZZ),g(:,OP_DZ ),h(:,OP_1),weight_79,79) &
          + int4(e(:,OP_DZ),f(:,OP_DRZ),g(:,OP_DR ),h(:,OP_1),weight_79,79) &
          + int4(e(:,OP_DZ),f(:,OP_DZ ),g(:,OP_DZZ),h(:,OP_1),weight_79,79) &
          + int4(e(:,OP_DZ),f(:,OP_DR ),g(:,OP_DRZ),h(:,OP_1),weight_79,79) &
          + int4(e(:,OP_DR),f(:,OP_DRZ),g(:,OP_DZ ),h(:,OP_1),weight_79,79) &
          + int4(e(:,OP_DR),f(:,OP_DRR),g(:,OP_DR ),h(:,OP_1),weight_79,79) &
          + int4(e(:,OP_DR),f(:,OP_DZ ),g(:,OP_DRZ),h(:,OP_1),weight_79,79) &
          + int4(e(:,OP_DR),f(:,OP_DR ),g(:,OP_DRR),h(:,OP_1),weight_79,79)
     if(itor.eq.1) then
        temp = temp &
          - int5(ri_79,e(:,OP_DR),f(:,OP_DZ ),g(:,OP_DZ ),h(:,OP_1 ),weight_79,79) &
          - int5(ri_79,e(:,OP_DR),f(:,OP_DR ),g(:,OP_DR ),h(:,OP_1 ),weight_79,79) &
          + int5(ri_79,e(:,OP_1 ),f(:,OP_DRZ),g(:,OP_DZ ),h(:,OP_1 ),weight_79,79) &
          + int5(ri_79,e(:,OP_1 ),f(:,OP_DRR),g(:,OP_DR ),h(:,OP_1 ),weight_79,79) &
          + int5(ri_79,e(:,OP_1 ),f(:,OP_DZ ),g(:,OP_DRZ),h(:,OP_1 ),weight_79,79) &
          + int5(ri_79,e(:,OP_1 ),f(:,OP_DR ),g(:,OP_DRR),h(:,OP_1 ),weight_79,79) &
          - int5(ri_79,e(:,OP_DR),f(:,OP_DZ ),g(:,OP_DZ ),h(:,OP_DR),weight_79,79) &
          - int5(ri_79,e(:,OP_DR),f(:,OP_DR ),g(:,OP_DR ),h(:,OP_DR),weight_79,79)
     endif
  endif

  temp = temp / 2.
  v3chichin = temp

  return
end function v3chichin

!===============================================================================
! B1 TERMS
!===============================================================================


! B1psi
! =====
subroutine b1psi(itri,i,j,ssterm, ddterm, rrterm, qqterm)

  use basic
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm, qqterm

  real :: temp

  temp = int2(g25(:,OP_1,i),g25(:,OP_1,j),weight_25,25)
  ssterm(1,1) = ssterm(1,1) + temp
  ddterm(1,1) = ddterm(1,1) + temp

end subroutine b1psi


! B1psiu (flux advection)
! =======================
subroutine b1psiu(itri,i,j,ssterm, ddterm, rrterm, qqterm)

  use basic
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm, qqterm

  real :: temp

  temp = int4(r_79,g79(:,OP_1,i),g79(:,OP_DR,j),pht79(:,OP_DZ),weight_79,79) &
       - int4(r_79,g79(:,OP_1,i),g79(:,OP_DZ,j),pht79(:,OP_DR),weight_79,79)

  ssterm(1,1) = ssterm(1,1) -     thimp *dt*temp
  ddterm(1,1) = ddterm(1,1) + (1.-thimp)*dt*temp
  
  temp = int4(r_79,g79(:,OP_1,i),ps179(:,OP_DR),g79(:,OP_DZ,j),weight_79,79) &
       - int4(r_79,g79(:,OP_1,i),ps179(:,OP_DZ),g79(:,OP_DR,j),weight_79,79)
  rrterm(1,1) = rrterm(1,1) + thimp*dt*temp
  qqterm(1,1) = qqterm(1,1) - thimp*dt*temp 

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     temp = int4(r_79,g79(:,OP_1,i),ps079(:,OP_DR),g79(:,OP_DZ,j),weight_79,79) &
          - int4(r_79,g79(:,OP_1,i),ps079(:,OP_DZ),g79(:,OP_DR,j),weight_79,79)
     rrterm(1,1) = rrterm(1,1) +     thimp *dt*temp
     qqterm(1,1) = qqterm(1,1) + (1.-thimp)*dt*temp
  endif

end subroutine b1psiu


! B1psieta
! ========
subroutine b1psieta(itri,i,j,ssterm, ddterm, rrterm, qqterm)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm, qqterm

  integer :: i1
  real :: temp

  temp =-int2(g25(:,OP_DR,i),g25(:,OP_DR,j),weight_25,25) &
       - int2(g25(:,OP_DZ,i),g25(:,OP_DZ,j),weight_25,25)
  if(itor.eq.1) then
     temp = temp - 2.*int3(ri_79,g79(:,OP_1,i),g79(:,OP_DR,j),weight_79,79)
  endif

  ssterm(1,1) = ssterm(1,1) - dt*    thimp *temp*etar
  ddterm(1,1) = ddterm(1,1) + dt*(1.-thimp)*temp*etar

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     i1 = isvaln(itri,i)

     temp =-int2(g25(:,OP_DR,i),ps025(:,OP_DR),weight_25,25) &
          - int2(g25(:,OP_DZ,i),ps025(:,OP_DZ),weight_25,25)
     if(itor.eq.1) then
        temp = temp - 2.*int3(ri_79,g79(:,OP_1,i),ps079(:,OP_DR),weight_79,79)
     endif

     q4(i1) = q4(i1) + dt*temp*etar
  endif

end subroutine b1psieta


! B1psibd
! =======
subroutine b1psibd(itri,i,j,ssterm, ddterm, rrterm, qqterm)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm, qqterm

  integer :: i2
  real :: temp

  ! Linearize in psi
  ! ~~~~~~~~~~~~~~~~
  if(idens.eq.0) then
     temp = int4(ri_79,g79(:,OP_1,i),g79(:,OP_DZ,j),bz179(:,OP_DR),weight_79,79) &
          - int4(ri_79,g79(:,OP_1,i),g79(:,OP_DR,j),bz179(:,OP_DZ),weight_79,79)
  else
     temp = int5(ri_79,g79(:,OP_1,i),g79(:,OP_DZ,j),bz179(:,OP_DR),ni79(:,OP_1),weight_79,79) &
          - int5(ri_79,g79(:,OP_1,i),g79(:,OP_DR,j),bz179(:,OP_DZ),ni79(:,OP_1),weight_79,79)
  endif

  ssterm(1,1) = ssterm(1,1) -      thimp *dt*temp*dbf
  ddterm(1,1) = ddterm(1,1) + (0.5-thimp)*dt*temp*dbf

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     if(idens.eq.0) then
        temp = int4(ri_79,g79(:,OP_1,i),g79(:,OP_DZ,j),bz079(:,OP_DR),weight_79,79) &
             - int4(ri_79,g79(:,OP_1,i),g79(:,OP_DR,j),bz079(:,OP_DZ),weight_79,79)
     else
        temp = int5(ri_79,g79(:,OP_1,i),g79(:,OP_DZ,j),bz079(:,OP_DR),ni79(:,OP_1),weight_79,79) &
             - int5(ri_79,g79(:,OP_1,i),g79(:,OP_DR,j),bz079(:,OP_DZ),ni79(:,OP_1),weight_79,79)
     endif   

     ssterm(1,1) = ssterm(1,1) -     thimp *dt*temp*dbf
     ddterm(1,1) = ddterm(1,1) + (1.-thimp)*dt*temp*dbf
  endif


  ! Linearize in bz
  ! ~~~~~~~~~~~~~~~
  if(idens.eq.0) then
     temp = int4(ri_79,g79(:,OP_1,i),ps179(:,OP_DZ),g79(:,OP_DR,j),weight_79,79) &
          - int4(ri_79,g79(:,OP_1,i),ps179(:,OP_DR),g79(:,OP_DZ,j),weight_79,79)
  else
     temp = int5(ri_79,g79(:,OP_1,i),ps179(:,OP_DZ),g79(:,OP_DR,j),ni79(:,OP_1),weight_79,79) &
          - int5(ri_79,g79(:,OP_1,i),ps179(:,OP_DR),g79(:,OP_DZ,j),ni79(:,OP_1),weight_79,79)
  endif

  ssterm(1,2) = ssterm(1,2) -      thimp *dt*temp*dbf
  ddterm(1,2) = ddterm(1,2) + (0.5-thimp)*dt*temp*dbf

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     if(idens.eq.0) then
        temp = int4(ri_79,g79(:,OP_1,i),ps179(:,OP_DZ),g79(:,OP_DR,j),weight_79,79) &
             - int4(ri_79,g79(:,OP_1,i),ps179(:,OP_DR),g79(:,OP_DZ,j),weight_79,79)
     else
        temp = int5(ri_79,g79(:,OP_1,i),ps179(:,OP_DZ),g79(:,OP_DR,j),ni79(:,OP_1),weight_79,79) &
             - int5(ri_79,g79(:,OP_1,i),ps179(:,OP_DR),g79(:,OP_DZ,j),ni79(:,OP_1),weight_79,79)
     endif

     ssterm(1,2) = ssterm(1,2) -     thimp *dt*temp*dbf
     ddterm(1,2) = ddterm(1,2) + (1.-thimp)*dt*temp*dbf
  endif

  ! Linearize in 1/n + Equilibrium term
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(linear.eq.1 .or. eqsubtract.eq.1) then
     i2 = isvaln(itri,i) + 6

     temp = int5(ri_79,g79(:,OP_1,i),ps179(:,OP_DZ),g79(:,OP_DR,j),ni79(:,OP_1),weight_79,79) &
          - int5(ri_79,g79(:,OP_1,i),ps179(:,OP_DR),g79(:,OP_DZ,j),ni79(:,OP_1),weight_79,79)
     q4(i2) = q4(i2) + dt*temp*dbf

  endif
end subroutine b1psibd


!===============================================================================
! B2 TERMS
!===============================================================================

! B2b
! ===
subroutine b2b(itri,i,j,ssterm, ddterm, rrterm, qqterm)

  use basic
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm, qqterm

  real :: temp

  temp = int2(g25(:,OP_1,i),g25(:,OP_1,j),weight_25,25)

  ssterm(2,2) = ssterm(2,2) + temp
  ddterm(2,2) = ddterm(2,2) + temp

end subroutine b2b


! B2beta
! ======
subroutine b2beta(itri,i,j,ssterm, ddterm, rrterm, qqterm)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm, qqterm

  integer :: i2
  real :: temp

  temp = -int3(ri2_79,g79(:,OP_DZ,i),g79(:,OP_DZ,j),weight_79,79) &
       -  int3(ri2_79,g79(:,OP_DR,i),g79(:,OP_DR,j),weight_79,79)

  ssterm(2,2) = ssterm(2,2) -     thimp *dt*temp*etar
  ddterm(2,2) = ddterm(2,2) + (1.-thimp)*dt*temp*etar

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     i2 = isvaln(itri,i) + 6

     temp = -int3(ri2_79,g79(:,OP_DZ,i),bz079(:,OP_DZ),weight_79,79) &
          -  int3(ri2_79,g79(:,OP_DR,i),bz079(:,OP_DR),weight_79,79)

     q4(i2) = q4(i2) + dt*temp*etar
  endif

end subroutine b2beta


! B2bu
! ====
subroutine b2bu(itri,i,j,ssterm, ddterm, rrterm, qqterm)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm, qqterm

  integer :: i2
  real :: temp

  temp = int4(ri_79,g79(:,OP_1,i),g79(:,OP_DR,j),pht79(:,OP_DZ),weight_79,79) &
       - int4(ri_79,g79(:,OP_1,i),g79(:,OP_DZ,j),pht79(:,OP_DR),weight_79,79)

  ssterm(2,2) = ssterm(2,2) -     thimp *dt*temp
  ddterm(2,2) = ddterm(2,2) + (1.-thimp)*dt*temp

  temp = int4(ri_79,g79(:,OP_1,i),bzt79(:,OP_DR),g79(:,OP_DZ,j),weight_79,79) &
       - int4(ri_79,g79(:,OP_1,i),bzt79(:,OP_DZ),g79(:,OP_DR,j),weight_79,79)

  rrterm(2,1) = rrterm(2,1) + thimp*dt*temp
  qqterm(2,1) = qqterm(2,1) - thimp*dt*temp

  if(linear.eq.1 .or. eqsubtract.eq.1) then
       temp = int4(ri_79,g79(:,OP_1,i),bz079(:,OP_DR),g79(:,OP_DZ,j),weight_79,79) &
            - int4(ri_79,g79(:,OP_1,i),bz079(:,OP_DZ),g79(:,OP_DR,j),weight_79,79)

       qqterm(2,1) = qqterm(2,1) + dt*temp

       ! EQUILIBRIUM TERM
       i2 = isvaln(itri,i) + 6

       temp = int4(ri_79,g79(:,OP_1,i),bz079(:,OP_DR),ph079(:,OP_DZ),weight_79,79) &
            - int4(ri_79,g79(:,OP_1,i),bz079(:,OP_DZ),ph079(:,OP_DR),weight_79,79)

       r4(i2) = r4(i2) + dt*temp
  endif

end subroutine b2bu


! B2psiv
! ======
subroutine b2psiv(itri,i,j,ssterm, ddterm, rrterm, qqterm)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm, qqterm

  integer :: i2
  real :: temp

  ! Linearize in psi
  ! ~~~~~~~~~~~~~~~~
  temp = int4(ri2_79,g79(:,OP_1,i),g79(:,OP_DR,j),vzt79(:,OP_DZ),weight_79,79) &
       - int4(ri2_79,g79(:,OP_1,i),g79(:,OP_DZ,j),vzt79(:,OP_DR),weight_79,79)
  if(itor.eq.1) then
     temp = temp + int4(ri3_79,g79(:,OP_1,i),g79(:,OP_DZ,j),vzt79(:,OP_1),weight_79,79)
  endif

  ssterm(2,1) = ssterm(2,1) -     thimp *dt*temp
  ddterm(2,1) = ddterm(2,1) + (1.-thimp)*dt*temp

  ! Linearize in v
  ! ~~~~~~~~~~~~~~
  temp = int4(ri2_79,g79(:,OP_1,i),pst79(:,OP_DR),g79(:,OP_DZ,j),weight_79,79) &
       - int4(ri2_79,g79(:,OP_1,i),pst79(:,OP_DZ),g79(:,OP_DR,j),weight_79,79)
  if(itor.eq.1) then
     temp = temp + int4(ri3_79,g79(:,OP_1,i),pst79(:,OP_DZ),g79(:,OP_1,j),weight_79,79)
  endif
  rrterm(2,2) = rrterm(2,2) + thimp*dt*temp
  qqterm(2,2) = qqterm(2,2) - thimp*dt*temp

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     temp = int4(ri2_79,g79(:,OP_1,i),ps079(:,OP_DR),g79(:,OP_DZ,j),weight_79,79) &
          - int4(ri2_79,g79(:,OP_1,i),ps079(:,OP_DZ),g79(:,OP_DR,j),weight_79,79)
     if(itor.eq.1) then
        temp = temp + int4(ri3_79,g79(:,OP_1,i),ps079(:,OP_DZ),g79(:,OP_1,j),weight_79,79)
     endif
     qqterm(2,2) = qqterm(2,2) + dt*temp

     ! EQUILIBRIUM TERM
     i2 = isvaln(itri,i)
     temp = int4(ri2_79,g79(:,OP_1,i),ps079(:,OP_DR),vz079(:,OP_DZ),weight_79,79) &
          - int4(ri2_79,g79(:,OP_1,i),ps079(:,OP_DZ),vz079(:,OP_DR),weight_79,79)
     if(itor.eq.1) then
        temp = temp + int4(ri3_79,g79(:,OP_1,i),ps079(:,OP_DZ),vz079(:,OP_1),weight_79,79)
     endif
     r4(i2) = r4(i2) + dt*temp
  endif

end subroutine b2psiv


! B2psipsid
! =========
subroutine b2psipsid(itri,i,j,ssterm, ddterm, rrterm, qqterm)

  use basic
  use arrays
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm, qqterm

  integer :: i2
  real :: temp

  if(idens.eq.0) then
     temp = int4(ri3_79,g79(:,OP_DR,i),g79(:,OP_DZ,j),ps179(:,OP_GS),weight_79,79) &
          - int4(ri3_79,g79(:,OP_DZ,i),g79(:,OP_DR,j),ps179(:,OP_GS),weight_79,79) &
          + int4(ri3_79,g79(:,OP_DR,i),ps179(:,OP_DZ),g79(:,OP_GS,j),weight_79,79) &
          - int4(ri3_79,g79(:,OP_DZ,i),ps179(:,OP_DR),g79(:,OP_GS,j),weight_79,79)
  else
     temp = int5(ri3_79,g79(:,OP_DR,i),g79(:,OP_DZ,j),ps179(:,OP_GS),ni79(:,OP_1),weight_79,79) &
          - int5(ri3_79,g79(:,OP_DZ,i),g79(:,OP_DR,j),ps179(:,OP_GS),ni79(:,OP_1),weight_79,79) &
          + int5(ri3_79,g79(:,OP_DR,i),ps179(:,OP_DZ),g79(:,OP_GS,j),ni79(:,OP_1),weight_79,79) &
          - int5(ri3_79,g79(:,OP_DZ,i),ps179(:,OP_DR),g79(:,OP_GS,j),ni79(:,OP_1),weight_79,79)
  endif

  ssterm(2,1) = ssterm(2,1) -      thimp *dt*temp*dbf
  ddterm(2,1) = ddterm(2,1) + (0.5-thimp)*dt*temp*dbf

  
  if(linear.eq.1 .or. eqsubtract.eq.1) then
     if(idens.eq.0) then
        temp = int4(ri3_79,g79(:,OP_DR,i),g79(:,OP_DZ,j),ps079(:,OP_GS),weight_79,79) &
             - int4(ri3_79,g79(:,OP_DZ,i),g79(:,OP_DR,j),ps079(:,OP_GS),weight_79,79) &
             + int4(ri3_79,g79(:,OP_DR,i),ps079(:,OP_DZ),g79(:,OP_GS,j),weight_79,79) &
             - int4(ri3_79,g79(:,OP_DZ,i),ps079(:,OP_DR),g79(:,OP_GS,j),weight_79,79)
     else
        temp = int5(ri3_79,g79(:,OP_DR,i),g79(:,OP_DZ,j),ps079(:,OP_GS),ni79(:,OP_1),weight_79,79) &
             - int5(ri3_79,g79(:,OP_DZ,i),g79(:,OP_DR,j),ps079(:,OP_GS),ni79(:,OP_1),weight_79,79) &
             + int5(ri3_79,g79(:,OP_DR,i),ps079(:,OP_DZ),g79(:,OP_GS,j),ni79(:,OP_1),weight_79,79) &
             - int5(ri3_79,g79(:,OP_DZ,i),ps079(:,OP_DR),g79(:,OP_GS,j),ni79(:,OP_1),weight_79,79)
     endif
  
     ssterm(2,1) = ssterm(2,1) -     thimp *dt*temp*dbf
     ddterm(2,1) = ddterm(2,1) + (1.-thimp)*dt*temp*dbf

     ! Linearize in ni
     ! ~~~~~~~~~~~~~~~
     if(idens.eq.1) then
        i2 = isvaln(itri,i)

        temp = int5(ri3_79,g79(:,OP_DR,i),ps079(:,OP_DZ),ps079(:,OP_GS),ni79(:,OP_1),weight_79,79) &
             - int5(ri3_79,g79(:,OP_DZ,i),ps079(:,OP_DR),ps079(:,OP_GS),ni79(:,OP_1),weight_79,79)

        r4(i2) = r4(i2) + dt*temp*dbf
     endif
  endif
end subroutine b2psipsid


!===============================================================================
! B3 TERMS
!===============================================================================

! B3pe
! ====
subroutine b3pe(itri,i,j,ssterm, ddterm, rrterm, qqterm)

  use basic
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout), dimension(3,3) :: ssterm, ddterm, rrterm, qqterm

  real :: temp

  temp = int2(g25(:,OP_1,i),g25(:,OP_1,j),weight_25,25)

  ssterm(3,3) = ssterm(3,3) + temp
  ddterm(3,3) = ddterm(3,3) + temp

end subroutine b3pe


! B3psipsieta
! ~~~~~~~~~~~
real function b3psipsieta(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  
  real :: temp

  if(itor.eq.0) then
     b3psipsieta = 0.
     return
  endif

  temp = (gam-1.)*etar*int4(ri2_79,e(:,OP_1),f(:,OP_GS),g(:,OP_GS),weight_79,79)

  b3psipsieta = temp
  
  return
end function


! B3bbeta
! ~~~~~~~
real function b3bbeta(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g
  
  real :: temp

  if(itor.eq.0) then
     b3bbeta = 0.
     return
  endif

  temp = (gam-1.)*etar* &
       (int4(ri2_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),weight_79,79) &
       +int4(ri2_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DR),weight_79,79))

  b3bbeta = temp
  
  return
end function


! B3pebd
! ~~~~~~
real function b3pebd(e,f,g,h)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g,h
  
  real :: temp

  if(idens.eq.0) then
     temp = int4(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79) &
          - int4(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79)
  else 
     temp = int5(ri_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),h(:,OP_1),weight_79,79) &
          - int5(ri_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),h(:,OP_1),weight_79,79) &
          + gam* &
          (int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DR),h(:,OP_DZ),weight_79,79) &
          -int5(ri_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),h(:,OP_DR),weight_79,79))
  endif

  b3pebd = temp
  
  return
end function



!===============================================================================
! N1 TERMS
!===============================================================================

! N1n
! ===
subroutine n1n(itri,i,j,ssterm, ddterm, rrterm, qqterm)

  use basic
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout) :: ssterm, ddterm
  real, intent(inout), dimension(3) :: rrterm, qqterm

  real :: temp

  temp = int2(g25(:,OP_1,i),g25(:,OP_1,j),weight_25,25)

  ssterm = ssterm + temp
  ddterm = ddterm + temp

end subroutine n1n



! N1nu  (velocity diffusion)
! =========================
subroutine n1nu(itri,i,j,ssterm, ddterm, rrterm, qqterm)

  use basic
  use nintegrate_mod

  implicit none

  integer, intent(in) :: itri, i, j
  real, intent(inout) :: ssterm, ddterm
  real, intent(inout), dimension(3) :: rrterm, qqterm

  real :: temp

  ! Linearize in n
  ! ~~~~~~~~~~~~~~

  temp = int4(r_79, g79(:,OP_1,i), g79(:,OP_DR,j), pht79(:,OP_DZ), weight_79,79) &
       - int4(r_79, g79(:,OP_1,i), g79(:,OP_DZ,j), pht79(:,OP_DR), weight_79,79)
  if(itor.eq.1) then
     temp = temp + 2.*int3(g79(:,OP_1,i), g79(:,OP_1,j), pht79(:,OP_DZ), weight_79,79)
  endif

  ssterm = ssterm -     thimp *dt*temp
  ddterm = ddterm + (1.-thimp)*dt*temp


  ! Linearize in u
  ! ~~~~~~~~~~~~~~

  temp = int4(r_79, g79(:,OP_1,i), n179(:,OP_DR), g79(:,OP_DZ,j), weight_79,79) &
       - int4(r_79, g79(:,OP_1,i), n179(:,OP_DZ), g79(:,OP_DR,j), weight_79,79)
  if(itor.eq.1) then
     temp = temp + 2.*int3(g79(:,OP_1,i), n179(:,OP_1), g79(:,OP_DZ,j), weight_79,79)
  endif

  rrterm(1) = rrterm(1) + thimp*dt*temp
  qqterm(1) = qqterm(1) - thimp*dt*temp

  if(linear.eq.1 .or. eqsubtract.eq.1) then
     temp = int4(r_79, g79(:,OP_1,i), n079(:,OP_DR), g79(:,OP_DZ,j), weight_79,79) &
          - int4(r_79, g79(:,OP_1,i), n079(:,OP_DZ), g79(:,OP_DR,j), weight_79,79)
     if(itor.eq.1) then
        temp = temp + 2.*int3(g79(:,OP_1,i), n079(:,OP_1), g79(:,OP_DZ,j), weight_79,79)
     endif
     
     rrterm(1) = rrterm(1) +     thimp *dt*temp
     qqterm(1) = qqterm(1) + (1.-thimp)*dt*temp
  endif


end subroutine n1nu


!===============================================================================
! P1 TERMS
!===============================================================================

! P1pu
! ~~~~
real function p1pu(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g

  real :: temp

  temp = int4(r_79,e(:,OP_1),f(:,OP_DR),g(:,OP_DZ),weight_79,79) &
       - int4(r_79,e(:,OP_1),f(:,OP_DZ),g(:,OP_DR),weight_79,79)
  if(itor.eq.1) then
     temp = temp + 2.*gam*int4(r_79,e(:,OP_1),f(:,OP_1),g(:,OP_DZ),weight_79,79)
  endif

  p1pu = temp

  return
end function p1pu


! P1pchi
! ~~~~~~
real function p1pchi(e,f,g)

  use basic
  use nintegrate_mod

  implicit none

  real, intent(in), dimension(79,OP_NUM) :: e,f,g

  real :: temp

  temp = gam* &
       (int3(e(:,OP_DZ),f(:,OP_1),g(:,OP_DZ),weight_79,79)  &
       +int3(e(:,OP_DR),f(:,OP_1),g(:,OP_DR),weight_79,79)) &
       +(gam-1.)* &
       (int3(e(:,OP_1),f(:,OP_DZ),g(:,OP_DZ),weight_79,79)  &
       +int3(e(:,OP_1),f(:,OP_DR),g(:,OP_DR),weight_79,79))

  p1pchi = temp

  return
end function p1pchi
