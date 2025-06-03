! Module containing terms 
!
module bootstrap
  
  implicit none
  integer, private, parameter :: dp = selected_real_kind(15,307)
  integer :: ibootstrap_model
  real :: ibootstrap_regular
  integer :: ibootstrap_map_te  ! Switch for bootstrap current coefficients mapping either as a function of (0) Psi_normal or (1) Te 
  ! 1 : add -eta*J_BS term to Ohm's law
  !     where J_BS = jbscommon * B
    
    !ibootstrap_model=1: Sauter & Angioni (1999) 

    !A = d lnp  /d psi    = (1/R psi_z + f'_r) p_z    + (1/R psi_r - f'_z) p_r
    !B = d lnTe /d psi    = (1/R psi_z + f'_r) Te_z/Te + (1/R psi_r - f'_z) Te_r/Te
    !C = d lnTi /d psi    = (1/R psi_z + f'_r) Ti_z/Te + (1/R psi_r - f'_z) Ti_r/Te
    !tempD =  L31 (A) + L32 Pe (B) + L34  alpha (P-Pe) (C)

    !ibootstrap_model=2: Redl et al (2021) 
    !tempD =  p L31 (A) + (L31+L32) Pe (B) + (L31+L34alpha)  (P-Pe) (C)
    !A    = (1/R psi_z + f'_r) nt_z/nt    + (1/R psi_r - f'_z) nt_r/nt


    

    !jbscommon   = -  F / <B^2>   [ L31 (dlnp/dpsi) + L32 Pe (dlnTe/dpsi) + L34 Pe alpha (dlnTi/dpsi) ]
    ! but d/dpsi = - del phi . (B x del)/||Bp||^2
    !     d/dpsi = - 1/||Bp||^2 1/R ( (1/R psi_z + f'_r) d/dz + (1/R psi_r - f'_z) d/dr)
    ! which gives
    !jbscommon   =  tempD 1/|Bp|^2  1 / <B^2> 1/R F


    ! L31,L32,L34,alpha,1/ <B^2>,dtedpsit,G are input as a function of psi or Te
    !via an external file named ProfileJBSCoeff_Psi_L31_32_34_alpha_B2_dtedpsit_G OR 
    ! ProfileJBSCoeff_Te_L31_32_34_alpha_B2_dtedpsit_G
    !bootsrap_alpha = amplification factor  
  real :: bootstrap_alpha

contains


!calculation of the bootstrap term in the pressure eq remains the same since it is linearized only in pe
subroutine bootstrap_pressure(trial, lin, ssterm, ddterm, pp_g, thimpi)
  use basic
  use arrays
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element,MAX_PTS, OP_NUM), intent(in) :: trial
  vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: lin
  vectype, dimension(dofs_per_element,num_fields), intent(inout) :: ssterm, ddterm
  integer, intent(in) :: pp_g
  real, intent(in) :: thimpi

  vectype, dimension(dofs_per_element) :: temp

  temp = bs_b3pe(trial,lin)
  ssterm(:,pp_g) = ssterm(:,pp_g) -       thimpi     *dt*temp
  ddterm(:,pp_g) = ddterm(:,pp_g) + (1. - thimpi*bdf)*dt*temp
end subroutine bootstrap_pressure




  
  ! B3pe
  ! ====
function bs_b3pe(e,f)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element) :: bs_b3pe
  vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(MAX_PTS) :: tempDD, tempAA, tempBB, tempCC, iBpsq



   
     !ibootstrap_model=1,3: Sauter & Angioni (1999) ! equivalent to 1 but a simplified version
    !A = d lnp  /d psi    = (1/R psi_z + f'_r) p_z    + (1/R psi_r - f'_z) p_r
    !B = d lnTe /d psi    = (1/R psi_z + f'_r) Te_z/Te + (1/R psi_r - f'_z) Te_r/Te
    !C = d lnTi /d psi    = (1/R psi_z + f'_r) Ti_z/Te + (1/R psi_r - f'_z) Ti_r/Te
    !tempD =  L31 (A) + L32 Pe (B) + L34  alpha (P-Pe) (C)
   
     !ibootstrap_model=2,4: Redl et al (2021) ! equivalent to 2 but a simplified version
     !tempD =  p L31 (A) + (L31+L32) Pe (B) + (L31+L34alpha)  (P-Pe) (C)
     !A    = (1/R psi_z + f'_r) nt_z/nt    + (1/R psi_r - f'_z) nt_r/nt
   
   
    !jbscommon   = -  F / <B^2>   [ L31 (dlnp/dpsi) + L32 Pe (dlnTe/dpsi) + L34 Pe alpha (dlnTi/dpsi) ]
    ! but d/dpsi = - del phi . (B x del)/||Bp||^2
    !     d/dpsi = - 1/||Bp||^2 1/R ( (1/R psi_z + f'_r) d/dz + (1/R psi_r - f'_z) d/dr)
    ! which gives
    !jbscommon   =  tempD 1/|Bp|^2  1 / <B^2> 1/R F
   
    !temp79a =  jbscommon * bootsrap_alpha = tempD 1/|Bp|^2  1 / <B^2>* 1/R * F * bootsrap_alpha
    


!linearizing in pe

  tempBB = (pst79(:,OP_DZ)*ri_79 + bfpt79(:,OP_DR))*tet79(:,OP_DZ)/tet79(:,OP_1) &
          + (pst79(:,OP_DR)*ri_79 - bfpt79(:,OP_DZ))*tet79(:,OP_DR)/tet79(:,OP_1)

  tempCC =  (pst79(:,OP_DZ)*ri_79 + bfpt79(:,OP_DR))*tit79(:,OP_DZ)/tit79(:,OP_1) &
          + (pst79(:,OP_DR)*ri_79 - bfpt79(:,OP_DZ))*tit79(:,OP_DR)/tit79(:,OP_1)

  if(ibootstrap_model.eq.1 .or. ibootstrap_model.eq.3)then !Sauter & Angioni (1999) 
      tempAA = (pst79(:,OP_DZ)*ri_79 + bfpt79(:,OP_DR))*pt79(:,OP_DZ) &
              + (pst79(:,OP_DR)*ri_79 - bfpt79(:,OP_DZ))*pt79(:,OP_DR)

      tempDD = jbsl3179(:,OP_1)*(tempAA) + &
               jbsl3279(:,OP_1)*f(:,OP_1)*(tempBB) + &
               jbsl3479(:,OP_1)*jbsalpha79(:,OP_1)*(pt79(:,OP_1)-f(:,OP_1))*(tempCC)
  else if (ibootstrap_model.eq.2 .or. ibootstrap_model.eq.4)then !Redl et al (2021) 
      tempAA = (pst79(:,OP_DZ)*ri_79 + bfpt79(:,OP_DR))*nt79(:,OP_DZ)/nt79(:,OP_1)*pt79(:,OP_1) &
               +  (pst79(:,OP_DR)*ri_79 - bfpt79(:,OP_DZ))*nt79(:,OP_DR)/nt79(:,OP_1)*pt79(:,OP_1)

      tempDD = jbsl3179(:,OP_1)*(tempAA) + &
               (jbsl3179(:,OP_1)+jbsl3279(:,OP_1))*f(:,OP_1)*(tempBB) + &
               (jbsl3179(:,OP_1)+jbsl3479(:,OP_1)*jbsalpha79(:,OP_1))*(pt79(:,OP_1)-f(:,OP_1))*(tempCC)
  end if
 
  iBpsq(:) = 1./((-ri_79*pst79(:,OP_DZ)-bfpt79(:,OP_DR))*(-ri_79*pst79(:,OP_DZ)-bfpt79(:,OP_DR)) &
              + (ri_79*pst79(:,OP_DR)-bfpt79(:,OP_DZ))*(ri_79*pst79(:,OP_DR)-bfpt79(:,OP_DZ)))

  temp79a=tempDD*jbsfluxavg_iBsq_B79(:,OP_1)*iBpsq(:)*ri_79*bzt79(:,OP_1)*bootstrap_alpha


  ! J.B
  temp79b = -ri2_79*bzt79(:,OP_1)*pst79(:,OP_GS) &
       + ri2_79*(bzt79(:,OP_DZ)*pst79(:,OP_DZ) + bzt79(:,OP_DR)*pst79(:,OP_DR))
#if defined(USE3D) || defined(USECOMPLEX)
  temp79b = temp79b &
       + ri2_79*(bfpt79(:,OP_DZP)*pst79(:,OP_DZ)  + bfpt79(:,OP_DRP)*pst79(:,OP_DR)) &
       + ri_79 *(bzt79(:,OP_DZ)*  bfpt79(:,OP_DR) - bzt79(:,OP_DR)*  bfpt79(:,OP_DZ)) &
       + ri_79 *(bfpt79(:,OP_DZP)*bfpt79(:,OP_DR) - bfpt79(:,OP_DRP)*bfpt79(:,OP_DZ)) &
       - ri3_79*(pst79(:,OP_DZ)*  pst79(:,OP_DRP) - pst79(:,OP_DR)*  pst79(:,OP_DZP)) &
       - ri2_79*(bfpt79(:,OP_DZ)* pst79(:,OP_DZP) + bfpt79(:,OP_DR)* pst79(:,OP_DRP))
#endif

  temp = intx3(e(:,:,OP_1),temp79a, temp79b)

  bs_b3pe = -(gam-1.)*temp

end function bs_b3pe




        ! B1psifpsib
   ! =========
   function bs_b1psifpsib(e,f,g,h,i)
    !mu,psi,f',psi,F
     !Sauter & Angioni (1999) 
    !A = d lnp  /d psi    = (1/R psi_z + f'_r) p_z    + (1/R psi_r - f'_z) p_r
    !B = d lnTe /d psi    = (1/R psi_z + f'_r) Te_z/Te + (1/R psi_r - f'_z) Te_r/Te
    !C = d lnTi /d psi    = (1/R psi_z + f'_r) Ti_z/Te + (1/R psi_r - f'_z) Ti_r/Te
    !tempD =  L31 (A) + L32 Pe (B) + L34  alpha (P-Pe) (C)

    !Redl et al (2021)  
    !tempD =  p L31 (A) + (L31+L32) Pe (B) + (L31+L34alpha)  (P-Pe) (C)
    !A    = (1/R psi_z + f'_r) nt_z/nt    + (1/R psi_r - f'_z) nt_r/nt


    !jbscommon=  tempD 1/|Bp|^2  1 / <B^2> 1/R F

    !temp = eta jbscommon ( (-1/r^3) (mu'_z psi_r - mu'_r psi_z) )
    !bootsrap_alpha = amplification factor  

    use basic
    use m3dc1_nint

    implicit none
    
    vectype, dimension(dofs_per_element) :: bs_b1psifpsib
    vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
    vectype, dimension(MAX_PTS) :: tempDD, tempAA, tempBB, tempCC, iBpsq
    vectype, dimension(dofs_per_element) :: temp
    
    temp = 0.

#if defined(USECOMPLEX) || defined(USE3D)
    if(jadv.eq.0) then
       temp = 0.
    else 
       if(surface_int) then
        temp = 0.
       else          

          tempBB = (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*tet79(:,OP_DZ)/tet79(:,OP_1) &
                 + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*tet79(:,OP_DR)/tet79(:,OP_1)
        
          tempCC =  (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*tit79(:,OP_DZ)/tit79(:,OP_1) &
                 + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*tit79(:,OP_DR)/tit79(:,OP_1)

          if(ibootstrap_model.eq.1)then !Sauter & Angioni (1999)
            tempAA = (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*pt79(:,OP_DZ) &
                 + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*pt79(:,OP_DR)

            tempDD = jbsl3179(:,OP_1)*(tempAA) + &
                     jbsl3279(:,OP_1)*pet79(:,OP_1)*(tempBB) + &
                     jbsl3479(:,OP_1)*jbsalpha79(:,OP_1)*(pt79(:,OP_1)-pet79(:,OP_1))*(tempCC)
          else if (ibootstrap_model.eq.2)then !Redl et al (2021)
            tempAA = (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*nt79(:,OP_DZ)/nt79(:,OP_1)*pt79(:,OP_1) &
                   + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*nt79(:,OP_DR)/nt79(:,OP_1)*pt79(:,OP_1)

            tempDD = jbsl3179(:,OP_1)*(tempAA) + &
                     (jbsl3179(:,OP_1)+jbsl3279(:,OP_1))*pet79(:,OP_1)*(tempBB) + &
                     (jbsl3179(:,OP_1)+jbsl3479(:,OP_1)*jbsalpha79(:,OP_1))*(pt79(:,OP_1)-pet79(:,OP_1))*(tempCC)
          end if

          iBpsq(:) = 1./((-ri_79*pst79(:,OP_DZ)-bfpt79(:,OP_DR))*(-ri_79*pst79(:,OP_DZ)-bfpt79(:,OP_DR)) &
              + (ri_79*pst79(:,OP_DR)-bfpt79(:,OP_DZ))*(ri_79*pst79(:,OP_DR)-bfpt79(:,OP_DZ)))

          temp79a = jbsfluxavg_iBsq_B79(:,OP_1)*iBpsq(:)*eta79(:,OP_1)*i(:,OP_1)* tempDD 

          temp = intx4(e(:,:,OP_DRP),ri4_79,h(:,OP_DZ),temp79a) &
               - intx4(e(:,:,OP_DZP),ri4_79,h(:,OP_DR),temp79a)

#ifdef USECOMPLEX
          temp = temp - rfac* &
                (intx4(e(:,:,OP_DR),ri4_79,h(:,OP_DZ),temp79a) &
               - intx4(e(:,:,OP_DZ),ri4_79,h(:,OP_DR),temp79a)) 
#endif

       endif
    endif
#endif

    bs_b1psifpsib = bootstrap_alpha*temp



end function bs_b1psifpsib


  ! B1psifbb
! =======
function bs_b1psifbb(e,f,g,h,i)
   !mu,psi,f',F,F

    !Sauter & Angioni (1999) 
    !A = d lnp  /d psi    = (1/R psi_z + f'_r) p_z    + (1/R psi_r - f'_z) p_r
    !B = d lnTe /d psi    = (1/R psi_z + f'_r) Te_z/Te + (1/R psi_r - f'_z) Te_r/Te
    !C = d lnTi /d psi    = (1/R psi_z + f'_r) Ti_z/Te + (1/R psi_r - f'_z) Ti_r/Te
    !tempD =  L31 (A) + L32 Pe (B) + L34  alpha (P-Pe) (C)

    !Redl et al (2021)  
    !tempD =  p L31 (A) + (L31+L32) Pe (B) + (L31+L34alpha)  (P-Pe) (C)
    !A    = (1/R psi_z + f'_r) nt_z/nt    + (1/R psi_r - f'_z) nt_r/nt


    !jbscommon=  tempD 1/|Bp|^2  1 / <B^2> 1/R F

    !temp = eta jbscommon (1/r^2 F (del* mu ))
    !bootsrap_alpha = amplification factor  
    use basic
    use m3dc1_nint

    implicit none
    
    vectype, dimension(dofs_per_element) :: bs_b1psifbb
    vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
    vectype, dimension(MAX_PTS) :: tempDD, tempAA, tempBB, tempCC, iBpsq
    vectype, dimension(dofs_per_element) :: temp

    temp = 0.

    if(jadv.eq.0) then
       if(surface_int) then
          temp = 0.
       else

          tempBB = (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*tet79(:,OP_DZ)/tet79(:,OP_1) &
                 + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*tet79(:,OP_DR)/tet79(:,OP_1)
        
          tempCC =  (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*tit79(:,OP_DZ)/tit79(:,OP_1) &
                 + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*tit79(:,OP_DR)/tit79(:,OP_1)
                 
          if(ibootstrap_model.eq.1)then !Sauter & Angioni (1999)
            tempAA = (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*pt79(:,OP_DZ) &
                 + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*pt79(:,OP_DR)

            tempDD = jbsl3179(:,OP_1)*(tempAA) + &
                     jbsl3279(:,OP_1)*pet79(:,OP_1)*(tempBB) + &
                     jbsl3479(:,OP_1)*jbsalpha79(:,OP_1)*(pt79(:,OP_1)-pet79(:,OP_1))*(tempCC)
          else if (ibootstrap_model.eq.2)then !Redl et al (2021)
            tempAA = (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*nt79(:,OP_DZ)/nt79(:,OP_1)*pt79(:,OP_1) &
                   + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*nt79(:,OP_DR)/nt79(:,OP_1)*pt79(:,OP_1)

            tempDD = jbsl3179(:,OP_1)*(tempAA) + &
                     (jbsl3179(:,OP_1)+jbsl3279(:,OP_1))*pet79(:,OP_1)*(tempBB) + &
                     (jbsl3179(:,OP_1)+jbsl3479(:,OP_1)*jbsalpha79(:,OP_1))*(pt79(:,OP_1)-pet79(:,OP_1))*(tempCC)
          end if

          iBpsq(:) = 1./((-ri_79*pst79(:,OP_DZ)-bfpt79(:,OP_DR))*(-ri_79*pst79(:,OP_DZ)-bfpt79(:,OP_DR)) &
              + (ri_79*pst79(:,OP_DR)-bfpt79(:,OP_DZ))*(ri_79*pst79(:,OP_DR)-bfpt79(:,OP_DZ)))

          temp79a = jbsfluxavg_iBsq_B79(:,OP_1)*iBpsq(:)*eta79(:,OP_1)*i(:,OP_1)* tempDD 
          temp = intx4(e(:,:,OP_1),ri_79,temp79a,h(:,OP_1))

       end if
    else 
       if(surface_int) then
          temp = 0.
       else

          tempBB = (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*tet79(:,OP_DZ)/tet79(:,OP_1) &
                 + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*tet79(:,OP_DR)/tet79(:,OP_1)
        
          tempCC =  (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*tit79(:,OP_DZ)/tit79(:,OP_1) &
                 + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*tit79(:,OP_DR)/tit79(:,OP_1)
                 
          if(ibootstrap_model.eq.1)then !Sauter & Angioni (1999)
            tempAA = (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*pt79(:,OP_DZ) &
                 + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*pt79(:,OP_DR)

            tempDD = jbsl3179(:,OP_1)*(tempAA) + &
                     jbsl3279(:,OP_1)*pet79(:,OP_1)*(tempBB) + &
                     jbsl3479(:,OP_1)*jbsalpha79(:,OP_1)*(pt79(:,OP_1)-pet79(:,OP_1))*(tempCC)
          else if (ibootstrap_model.eq.2)then !Redl et al (2021)
            tempAA = (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*nt79(:,OP_DZ)/nt79(:,OP_1)*pt79(:,OP_1) &
                   + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*nt79(:,OP_DR)/nt79(:,OP_1)*pt79(:,OP_1)

            tempDD = jbsl3179(:,OP_1)*(tempAA) + &
                     (jbsl3179(:,OP_1)+jbsl3279(:,OP_1))*pet79(:,OP_1)*(tempBB) + &
                     (jbsl3179(:,OP_1)+jbsl3479(:,OP_1)*jbsalpha79(:,OP_1))*(pt79(:,OP_1)-pet79(:,OP_1))*(tempCC)
          end if

          iBpsq(:) = 1./((-ri_79*pst79(:,OP_DZ)-bfpt79(:,OP_DR))*(-ri_79*pst79(:,OP_DZ)-bfpt79(:,OP_DR)) &
              + (ri_79*pst79(:,OP_DR)-bfpt79(:,OP_DZ))*(ri_79*pst79(:,OP_DR)-bfpt79(:,OP_DZ)))

          temp79a = jbsfluxavg_iBsq_B79(:,OP_1)*iBpsq(:)*eta79(:,OP_1)*i(:,OP_1)* tempDD 
          temp = intx4(e(:,:,OP_GS),ri3_79,temp79a,h(:,OP_1))
 
          if(itor.eq.1) then
             temp = temp - 2.*intx4(e(:,:,OP_DR),ri4_79,temp79a,h(:,OP_1))
          end if
       endif
    endif

    bs_b1psifbb = bootstrap_alpha*temp
  end function bs_b1psifbb


 ! B1psiffb
! =======
  function bs_b1psiffb(e,f,g,h,i)
    !mu,psi,f',f',F

    !Sauter & Angioni (1999) 
    !A = d lnp  /d psi    = (1/R psi_z + f'_r) p_z    + (1/R psi_r - f'_z) p_r
    !B = d lnTe /d psi    = (1/R psi_z + f'_r) Te_z/Te + (1/R psi_r - f'_z) Te_r/Te
    !C = d lnTi /d psi    = (1/R psi_z + f'_r) Ti_z/Te + (1/R psi_r - f'_z) Ti_r/Te
    !tempD =  L31 (A) + L32 Pe (B) + L34  alpha (P-Pe) (C)

    !Redl et al (2021)  
    !tempD =  p L31 (A) + (L31+L32) Pe (B) + (L31+L34alpha)  (P-Pe) (C)
    !A    = (1/R psi_z + f'_r) nt_z/nt    + (1/R psi_r - f'_z) nt_r/nt


    !jbscommon=  tempD 1/|Bp|^2  1 / <B^2> 1/R F

    !temp = eta jbscommon (1/r^2 (psi_z mu'_z + psi_r mu'_r))
    !bootsrap_alpha = amplification factor  
    use basic
    use m3dc1_nint

    implicit none
    
    vectype, dimension(dofs_per_element) :: bs_b1psiffb
    vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
    vectype, dimension(MAX_PTS) :: tempDD, tempAA, tempBB, tempCC, iBpsq
    vectype, dimension(dofs_per_element) :: temp


    temp = 0.

#if defined(USECOMPLEX) || defined(USE3D)
    if(jadv.eq.0) then
       temp = 0.
    else 
       if(surface_int) then
          temp = 0.
       else

        tempBB = (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*tet79(:,OP_DZ)/tet79(:,OP_1) &
                 + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*tet79(:,OP_DR)/tet79(:,OP_1)
        
        tempCC =  (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*tit79(:,OP_DZ)/tit79(:,OP_1) &
                 + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*tit79(:,OP_DR)/tit79(:,OP_1)

        if(ibootstrap_model.eq.1)then !Sauter & Angioni (1999)
            tempAA = (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*pt79(:,OP_DZ) &
                 + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*pt79(:,OP_DR)

            tempDD = jbsl3179(:,OP_1)*(tempAA) + &
                     jbsl3279(:,OP_1)*pet79(:,OP_1)*(tempBB) + &
                     jbsl3479(:,OP_1)*jbsalpha79(:,OP_1)*(pt79(:,OP_1)-pet79(:,OP_1))*(tempCC)
        else if (ibootstrap_model.eq.2)then !Redl et al (2021)
            tempAA = (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*nt79(:,OP_DZ)/nt79(:,OP_1)*pt79(:,OP_1) &
                   + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*nt79(:,OP_DR)/nt79(:,OP_1)*pt79(:,OP_1)

            tempDD = jbsl3179(:,OP_1)*(tempAA) + &
                     (jbsl3179(:,OP_1)+jbsl3279(:,OP_1))*pet79(:,OP_1)*(tempBB) + &
                     (jbsl3179(:,OP_1)+jbsl3479(:,OP_1)*jbsalpha79(:,OP_1))*(pt79(:,OP_1)-pet79(:,OP_1))*(tempCC)
        end if

        iBpsq(:) = 1./((-ri_79*pst79(:,OP_DZ)-bfpt79(:,OP_DR))*(-ri_79*pst79(:,OP_DZ)-bfpt79(:,OP_DR)) &
              + (ri_79*pst79(:,OP_DR)-bfpt79(:,OP_DZ))*(ri_79*pst79(:,OP_DR)-bfpt79(:,OP_DZ)))
        
        temp79a = jbsfluxavg_iBsq_B79(:,OP_1)*iBpsq(:)*eta79(:,OP_1)*i(:,OP_1)* tempDD 

        temp = intx4(e(:,:,OP_DRP),ri3_79,temp79a,h(:,OP_DR)) &
             + intx4(e(:,:,OP_DZP),ri3_79,temp79a,h(:,OP_DZ))

#ifdef USECOMPLEX
        temp = temp - rfac* &
               (intx4(e(:,:,OP_DR),ri3_79,temp79a,h(:,OP_DR)) &
               +intx4(e(:,:,OP_DZ),ri3_79,temp79a,h(:,OP_DZ)))
#endif

       endif
    endif
#endif

    bs_b1psiffb = bootstrap_alpha*temp
  end function bs_b1psiffb



      ! B2psifpsib
   ! =========
  function bs_b2psifpsib(e,f,g,h,i)
    !mu,psi,f',psi,F

    !Sauter & Angioni (1999) 
    !A = d lnp  /d psi    = (1/R psi_z + f'_r) p_z    + (1/R psi_r - f'_z) p_r
    !B = d lnTe /d psi    = (1/R psi_z + f'_r) Te_z/Te + (1/R psi_r - f'_z) Te_r/Te
    !C = d lnTi /d psi    = (1/R psi_z + f'_r) Ti_z/Te + (1/R psi_r - f'_z) Ti_r/Te
    !tempD =  L31 (A) + L32 Pe (B) + L34  alpha (P-Pe) (C)

    !Redl et al (2021)  
    !tempD =  p L31 (A) + (L31+L32) Pe (B) + (L31+L34alpha)  (P-Pe) (C)
    !A    = (1/R psi_z + f'_r) nt_z/nt    + (1/R psi_r - f'_z) nt_r/nt


    !jbscommon=  tempD 1/|Bp|^2  1 / <B^2> 1/R F

    !temp = eta jbscommon (1/r^2 (psi_z mu_z + psi_r mu_r))
    !bootsrap_alpha = amplification factor  

    use basic
    use m3dc1_nint

    implicit none

    vectype, dimension(dofs_per_element) :: bs_b2psifpsib
    vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
    vectype, dimension(MAX_PTS) :: tempDD, tempAA, tempBB, tempCC, iBpsq
    vectype, dimension(dofs_per_element) :: temp


    temp = 0.

    tempBB = (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*tet79(:,OP_DZ)/tet79(:,OP_1) &
                 + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*tet79(:,OP_DR)/tet79(:,OP_1)
        
    tempCC =  (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*tit79(:,OP_DZ)/tit79(:,OP_1) &
                 + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*tit79(:,OP_DR)/tit79(:,OP_1)

    if(ibootstrap_model.eq.1)then !Sauter & Angioni (1999)
            tempAA = (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*pt79(:,OP_DZ) &
                 + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*pt79(:,OP_DR)

            tempDD = jbsl3179(:,OP_1)*(tempAA) + &
                     jbsl3279(:,OP_1)*pet79(:,OP_1)*(tempBB) + &
                     jbsl3479(:,OP_1)*jbsalpha79(:,OP_1)*(pt79(:,OP_1)-pet79(:,OP_1))*(tempCC)
    else if (ibootstrap_model.eq.2)then !Redl et al (2021)
            tempAA = (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*nt79(:,OP_DZ)/nt79(:,OP_1)*pt79(:,OP_1) &
                   + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*nt79(:,OP_DR)/nt79(:,OP_1)*pt79(:,OP_1)

            tempDD = jbsl3179(:,OP_1)*(tempAA) + &
                     (jbsl3179(:,OP_1)+jbsl3279(:,OP_1))*pet79(:,OP_1)*(tempBB) + &
                     (jbsl3179(:,OP_1)+jbsl3479(:,OP_1)*jbsalpha79(:,OP_1))*(pt79(:,OP_1)-pet79(:,OP_1))*(tempCC)
    end if

    iBpsq(:) = 1./((-ri_79*pst79(:,OP_DZ)-bfpt79(:,OP_DR))*(-ri_79*pst79(:,OP_DZ)-bfpt79(:,OP_DR)) &
              + (ri_79*pst79(:,OP_DR)-bfpt79(:,OP_DZ))*(ri_79*pst79(:,OP_DR)-bfpt79(:,OP_DZ)))

    temp79a = jbsfluxavg_iBsq_B79(:,OP_1)*iBpsq(:)*eta79(:,OP_1)*i(:,OP_1)* tempDD 

    temp = intx4(e(:,:,OP_DZ),ri3_79,h(:,OP_DZ),temp79a) &
         + intx4(e(:,:,OP_DR),ri3_79,h(:,OP_DR),temp79a)

    bs_b2psifpsib = bootstrap_alpha*temp
  end function bs_b2psifpsib




! B2psiffb
 ! =======
 function bs_b2psiffb(e,f,g,h,i)

    !mu,psi,f',f',F
    !Sauter & Angioni (1999) 
    !A = d lnp  /d psi    = (1/R psi_z + f'_r) p_z    + (1/R psi_r - f'_z) p_r
    !B = d lnTe /d psi    = (1/R psi_z + f'_r) Te_z/Te + (1/R psi_r - f'_z) Te_r/Te
    !C = d lnTi /d psi    = (1/R psi_z + f'_r) Ti_z/Te + (1/R psi_r - f'_z) Ti_r/Te
    !tempD =  L31 (A) + L32 Pe (B) + L34  alpha (P-Pe) (C)

    !Redl et al (2021)  
    !tempD =  p L31 (A) + (L31+L32) Pe (B) + (L31+L34alpha)  (P-Pe) (C)
    !A    = (1/R psi_z + f'_r) nt_z/nt    + (1/R psi_r - f'_z) nt_r/nt
    
    !jbscommon=  tempD 1/|Bp|^2  1 / <B^2> 1/R F

    !temp = eta jbscommon (1/r   ( mu_z f'_r - mu_r f'_z))
    !bootsrap_alpha = amplification factor  

    use basic
    use m3dc1_nint

    implicit none

    vectype, dimension(dofs_per_element) :: bs_b2psiffb
    vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
    vectype, dimension(MAX_PTS) :: tempDD, tempAA, tempBB, tempCC, iBpsq
    vectype, dimension(dofs_per_element) :: temp

    temp = 0.

#if defined(USECOMPLEX) || defined(USE3D)
    tempBB = (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*tet79(:,OP_DZ)/tet79(:,OP_1) &
            + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*tet79(:,OP_DR)/tet79(:,OP_1)

    tempCC =  (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*tit79(:,OP_DZ)/tit79(:,OP_1) &
            + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*tit79(:,OP_DR)/tit79(:,OP_1)

    if(ibootstrap_model.eq.1)then !Sauter & Angioni (1999) 
        tempAA = (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*pt79(:,OP_DZ) &
                + (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*pt79(:,OP_DR)

        tempDD = jbsl3179(:,OP_1)*(tempAA) + &
                 jbsl3279(:,OP_1)*pet79(:,OP_1)*(tempBB) + &
                 jbsl3479(:,OP_1)*jbsalpha79(:,OP_1)*(pt79(:,OP_1)-pet79(:,OP_1))*(tempCC)
    else if (ibootstrap_model.eq.2)then !Redl et al (2021) 
        tempAA = (f(:,OP_DZ)*ri_79 + g(:,OP_DR))*nt79(:,OP_DZ)/nt79(:,OP_1)*pt79(:,OP_1) &
                 +  (f(:,OP_DR)*ri_79 - g(:,OP_DZ))*nt79(:,OP_DR)/nt79(:,OP_1)*pt79(:,OP_1)

        tempDD = jbsl3179(:,OP_1)*(tempAA) + &
                 (jbsl3179(:,OP_1)+jbsl3279(:,OP_1))*pet79(:,OP_1)*(tempBB) + &
                 (jbsl3179(:,OP_1)+jbsl3479(:,OP_1)*jbsalpha79(:,OP_1))*(pt79(:,OP_1)-pet79(:,OP_1))*(tempCC)
    end if

    iBpsq(:) = 1./((-ri_79*pst79(:,OP_DZ)-bfpt79(:,OP_DR))*(-ri_79*pst79(:,OP_DZ)-bfpt79(:,OP_DR)) &
              + (ri_79*pst79(:,OP_DR)-bfpt79(:,OP_DZ))*(ri_79*pst79(:,OP_DR)-bfpt79(:,OP_DZ)))

    temp79a = jbsfluxavg_iBsq_B79(:,OP_1)*iBpsq(:)*eta79(:,OP_1)*i(:,OP_1)*tempDD 

    temp = intx4(e(:,:,OP_DZ),ri2_79,h(:,OP_DR),temp79a) &
         - intx4(e(:,:,OP_DR),ri2_79,h(:,OP_DZ),temp79a)
#else 
    temp = 0.
#endif 
    bs_b2psiffb = bootstrap_alpha*temp


  end function bs_b2psiffb





  subroutine bootstrap_flux(trial, lin, ssterm, ddterm, r_bf, q_bf, thimpf, thimp_bf)
    use basic
    use arrays
    use m3dc1_nint

    implicit none

    vectype, dimension(dofs_per_element, MAX_PTS, OP_NUM), intent(in) :: trial
    vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: lin 
    vectype, dimension(dofs_per_element, num_fields), intent(inout) :: ssterm, ddterm
    vectype, dimension(dofs_per_element), intent(out) :: r_bf, q_bf
    real, intent(in) :: thimpf, thimp_bf

    vectype, dimension(dofs_per_element) :: temp
    
    if(numvar.eq.1) then
       ! linearizing in psi  
       temp = bs_b1psifpsib(trial,lin,bfpt79,pst79,bzt79) &
            + bs_b1psifpsib(trial,pst79,bfpt79,lin,bzt79) 
       

        ssterm(:,psi_g) = ssterm(:,psi_g) -          thimpf     *dt*temp
        ddterm(:,psi_g) = ddterm(:,psi_g) + (1./2. - thimpf*bdf)*dt*temp
       
        ! linearizing in psi
       temp = bs_b1psifbb(trial,lin,bfpt79,bzt79,bzt79) &
            + bs_b1psiffb(trial,lin,bfpt79,bfpt79,bzt79) 

       ssterm(:,psi_g) = ssterm(:,psi_g) -          thimpf     *dt*temp
       ddterm(:,psi_g) = ddterm(:,psi_g) + (1. - thimpf*bdf)*dt*temp
    else     
        ! linearizing in psi
       temp = bs_b1psifpsib(trial,lin,bfpt79,pst79,bzt79) &
            + bs_b1psifpsib(trial,pst79,bfpt79,lin,bzt79) &
            + bs_b1psifbb(trial,lin,bfpt79,bzt79,bzt79) &
            + bs_b1psiffb(trial,lin,bfpt79,bfpt79,bzt79) 

       ssterm(:,psi_g) = ssterm(:,psi_g) -          thimpf     *dt*temp
       ddterm(:,psi_g) = ddterm(:,psi_g) + (3./12. - thimpf*bdf)*dt*temp
  
         ! linearizing in b  
       temp = bs_b1psifpsib(trial,pst79,bfpt79,pst79,lin) &
            + bs_b1psifbb(trial,pst79,bfpt79,lin,bzt79) &
            + bs_b1psifbb(trial,pst79,bfpt79,bzt79,lin) &
            + bs_b1psiffb(trial,pst79,bfpt79,bfpt79,lin) 

       
       ssterm(:,bz_g) = ssterm(:,bz_g) -          thimpf     *dt*temp
       ddterm(:,bz_g) = ddterm(:,bz_g) + (3./12. - thimpf*bdf)*dt*temp
       

     ! linearizing in f'  
       temp = bs_b1psifpsib(trial,pst79,lin,pst79,bzt79) &
            + bs_b1psifbb(trial,pst79,lin,bzt79,bzt79) &
            + bs_b1psiffb(trial,pst79,lin,bfpt79,bzt79) &
            + bs_b1psiffb(trial,pst79,bfpt79,lin,bzt79) 
       r_bf = r_bf -          thimp_bf     *dt*temp
       q_bf = q_bf + (3./12. - thimp_bf*bdf)*dt*temp            
       
       if(eqsubtract.eq.1) then
        !(at+a0)(bt+b0)(ct+c0)=(at bt + at b0 + a0 bt + a0 b0)(ct + c0)=
       != (at bt ct + at b0 ct + a0 bt ct + a0 b0 ct) + (at bt c0 + at b0 c0 + a0 bt c0 + a0 b0 c0)

        temp = bs_b1psifpsib(trial,lin,bfpt79,pst79,bz079) & !at bt c0  
             + bs_b1psifpsib(trial,lin,bfp079,pst79,bz079) & !a0 bt c0
             + bs_b1psifpsib(trial,lin,bfpt79,ps079,bz079) & !at b0 c0
             + bs_b1psifpsib(trial,lin,bfp079,ps079,bz079)*6 & !a0 b0 c0
             + bs_b1psifpsib(trial,lin,bfp079,pst79,bzt79) & !a0 bt ct
             + bs_b1psifpsib(trial,lin,bfpt79,ps079,bzt79) & !at b0 ct 
             + bs_b1psifpsib(trial,lin,bfp079,ps079,bzt79) & !a0 b0 ct
        !--------------------------------------------------------!
             + bs_b1psifpsib(trial,pst79,bfpt79,lin,bz079) & !at bt c0  
             + bs_b1psifpsib(trial,ps079,bfpt79,lin,bz079) & !a0 bt c0
             + bs_b1psifpsib(trial,pst79,bfp079,lin,bz079) & !at b0 c0
             + bs_b1psifpsib(trial,ps079,bfp079,lin,bz079)*6 & !a0 b0 c0
             + bs_b1psifpsib(trial,ps079,bfpt79,lin,bzt79) & !a0 bt ct
             + bs_b1psifpsib(trial,pst79,bfp079,lin,bzt79) & !at b0 ct 
             + bs_b1psifpsib(trial,ps079,bfp079,lin,bzt79) & !a0 b0 ct
         !--------------------------------------------------------!
             + bs_b1psifbb  (trial,lin,bfpt79,bzt79,bz079) & !at bt c0  
             + bs_b1psifbb  (trial,lin,bfp079,bzt79,bz079) & !a0 bt c0
             + bs_b1psifbb  (trial,lin,bfpt79,bz079,bz079) & !at b0 c0
             + bs_b1psifbb  (trial,lin,bfp079,bz079,bz079)*6 & !a0 b0 c0
             + bs_b1psifbb  (trial,lin,bfp079,bzt79,bzt79) & !a0 bt ct
             + bs_b1psifbb  (trial,lin,bfpt79,bz079,bzt79) & !at b0 ct   
             + bs_b1psifbb  (trial,lin,bfp079,bz079,bzt79) &  !a0 b0 ct
         !--------------------------------------------------------! 
             + bs_b1psiffb  (trial,lin,bfpt79,bfpt79,bz079) & !at bt c0  
             + bs_b1psiffb  (trial,lin,bfp079,bfpt79,bz079) & !a0 bt c0
             + bs_b1psiffb  (trial,lin,bfpt79,bfp079,bz079) & !at b0 c0
             + bs_b1psiffb  (trial,lin,bfp079,bfp079,bz079)*6 & !a0 b0 c0
             + bs_b1psiffb  (trial,lin,bfp079,bfpt79,bzt79) & !a0 bt ct
             + bs_b1psiffb  (trial,lin,bfpt79,bfp079,bzt79) & !at b0 ct   
             + bs_b1psiffb  (trial,lin,bfp079,bfp079,bzt79)   !a0 b0 ct
         !--------------------------------------------------------!  
        ddterm(:,psi_g) = ddterm(:,psi_g) + (1./12.)*dt*temp      
          

        temp = bs_b1psifpsib(trial,pst79,bfpt79,ps079,lin) & !at bt c0  
             + bs_b1psifpsib(trial,ps079,bfpt79,ps079,lin) & !a0 bt c0
             + bs_b1psifpsib(trial,pst79,bfp079,ps079,lin) & !at b0 c0
             + bs_b1psifpsib(trial,ps079,bfp079,ps079,lin)*6 & !a0 b0 c0
             + bs_b1psifpsib(trial,ps079,bfpt79,pst79,lin) & !a0 bt ct
             + bs_b1psifpsib(trial,pst79,bfp079,pst79,lin) & !at b0 ct 
             + bs_b1psifpsib(trial,ps079,bfp079,pst79,lin) & !a0 b0 ct
        !--------------------------------------------------------!
             + bs_b1psifbb  (trial,pst79,bfpt79,lin,bz079) & !at bt c0  
             + bs_b1psifbb  (trial,ps079,bfpt79,lin,bz079) & !a0 bt c0
             + bs_b1psifbb  (trial,pst79,bfp079,lin,bz079) & !at b0 c0
             + bs_b1psifbb  (trial,ps079,bfp079,lin,bz079)*6 & !a0 b0 c0
             + bs_b1psifbb  (trial,ps079,bfpt79,lin,bzt79) & !a0 bt ct
             + bs_b1psifbb  (trial,pst79,bfp079,lin,bzt79) & !at b0 ct   
             + bs_b1psifbb  (trial,ps079,bfp079,lin,bzt79) &  !a0 b0 ct
         !--------------------------------------------------------!
             + bs_b1psifbb  (trial,pst79,bfpt79,bz079,lin) & !at bt c0  
             + bs_b1psifbb  (trial,ps079,bfpt79,bz079,lin) & !a0 bt c0
             + bs_b1psifbb  (trial,pst79,bfp079,bz079,lin) & !at b0 c0
             + bs_b1psifbb  (trial,ps079,bfp079,bz079,lin)*6 & !a0 b0 c0
             + bs_b1psifbb  (trial,ps079,bfpt79,bzt79,lin) & !a0 bt ct
             + bs_b1psifbb  (trial,pst79,bfp079,bzt79,lin) & !at b0 ct   
             + bs_b1psifbb  (trial,ps079,bfp079,bzt79,lin) &  !a0 b0 ct
         !--------------------------------------------------------! 
             + bs_b1psiffb  (trial,pst79,bfpt79,bfp079,lin) & !at bt c0  
             + bs_b1psiffb  (trial,ps079,bfpt79,bfp079,lin) & !a0 bt c0
             + bs_b1psiffb  (trial,pst79,bfp079,bfp079,lin) & !at b0 c0
             + bs_b1psiffb  (trial,ps079,bfp079,bfp079,lin)*6 & !a0 b0 c0
             + bs_b1psiffb  (trial,ps079,bfpt79,bfpt79,lin) & !a0 bt ct
             + bs_b1psiffb  (trial,pst79,bfp079,bfpt79,lin) & !at b0 ct   
             + bs_b1psiffb  (trial,ps079,bfp079,bfpt79,lin)   !a0 b0 ct
         !--------------------------------------------------------!  
        ddterm(:,bz_g) = ddterm(:,bz_g) + (1./12.)*dt*temp    

        temp = bs_b1psifpsib(trial,pst79,lin,pst79,bz079) & !at bt c0  
             + bs_b1psifpsib(trial,ps079,lin,pst79,bz079) & !a0 bt c0
             + bs_b1psifpsib(trial,pst79,lin,ps079,bz079) & !at b0 c0
             + bs_b1psifpsib(trial,ps079,lin,ps079,bz079)*6 & !a0 b0 c0
             + bs_b1psifpsib(trial,ps079,lin,pst79,bzt79) & !a0 bt ct
             + bs_b1psifpsib(trial,pst79,lin,ps079,bzt79) & !at b0 ct 
             + bs_b1psifpsib(trial,ps079,lin,ps079,bzt79) & !a0 b0 ct
        !--------------------------------------------------------!
             + bs_b1psifbb  (trial,pst79,lin,pst79,bz079) & !at bt c0  
             + bs_b1psifbb  (trial,ps079,lin,pst79,bz079) & !a0 bt c0
             + bs_b1psifbb  (trial,pst79,lin,ps079,bz079) & !at b0 c0
             + bs_b1psifbb  (trial,ps079,lin,ps079,bz079)*6 & !a0 b0 c0
             + bs_b1psifbb  (trial,ps079,lin,pst79,bzt79) & !a0 bt ct
             + bs_b1psifbb  (trial,pst79,lin,ps079,bzt79) & !at b0 ct   
             + bs_b1psifbb  (trial,ps079,lin,ps079,bzt79) &  !a0 b0 ct
         !--------------------------------------------------------!
             + bs_b1psiffb  (trial,pst79,lin,bfpt79,bz079) & !at bt c0  
             + bs_b1psiffb  (trial,ps079,lin,bfpt79,bz079) & !a0 bt c0
             + bs_b1psiffb  (trial,pst79,lin,bfp079,bz079) & !at b0 c0
             + bs_b1psiffb  (trial,ps079,lin,bfp079,bz079)*6 & !a0 b0 c0
             + bs_b1psiffb  (trial,ps079,lin,bfpt79,bzt79) & !a0 bt ct
             + bs_b1psiffb  (trial,pst79,lin,bfp079,bzt79) & !at b0 ct   
             + bs_b1psiffb  (trial,ps079,lin,bfp079,bzt79) &  !a0 b0 ct
         !--------------------------------------------------------! 
             + bs_b1psiffb  (trial,pst79,bfpt79,lin,bz079) & !at bt c0  
             + bs_b1psiffb  (trial,ps079,bfpt79,lin,bz079) & !a0 bt c0
             + bs_b1psiffb  (trial,pst79,bfp079,lin,bz079) & !at b0 c0
             + bs_b1psiffb  (trial,ps079,bfp079,lin,bz079)*6 & !a0 b0 c0
             + bs_b1psiffb  (trial,ps079,bfpt79,lin,bzt79) & !a0 bt ct
             + bs_b1psiffb  (trial,pst79,bfp079,lin,bzt79) & !at b0 ct   
             + bs_b1psiffb  (trial,ps079,bfp079,lin,bzt79)   !a0 b0 ct
         !--------------------------------------------------------!      
         q_bf = q_bf + (1./12.)*dt*temp    
       end if
    end if
  end subroutine bootstrap_flux







   subroutine bootstrap_axial_field(trial, lin, ssterm, ddterm, r_bf, q_bf, thimpf, thimp_bf)
    use basic
    use arrays
    use m3dc1_nint

    implicit none

    vectype, dimension(dofs_per_element, MAX_PTS, OP_NUM), intent(in) :: trial
    vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: lin
    vectype, dimension(dofs_per_element, num_fields), intent(inout) :: ssterm, ddterm
    vectype, dimension(dofs_per_element), intent(out) :: r_bf, q_bf
    real, intent(in) :: thimpf, thimp_bf

    vectype, dimension(dofs_per_element) :: temp


  ! linearizing in psi
    temp = bs_b2psifpsib(trial,lin,bfpt79,pst79,bzt79) &
         + bs_b2psifpsib(trial,pst79,bfpt79,lin,bzt79) &
         + bs_b2psiffb  (trial,lin,bfpt79,bfpt79,bzt79)
    ssterm(:,psi_g) = ssterm(:,psi_g) -          thimpf     *dt*temp
    ddterm(:,psi_g) = ddterm(:,psi_g) + (3./12. - thimpf*bdf)*dt*temp
    

  ! linearizing in F
    temp = bs_b2psifpsib(trial,pst79,bfpt79,pst79,lin) &
         + bs_b2psiffb  (trial,pst79,bfpt79,bfpt79,lin)
    ssterm(:,bz_g) = ssterm(:,bz_g) -          thimpf     *dt*temp
    ddterm(:,bz_g) = ddterm(:,bz_g) + (3./12. - thimpf*bdf)*dt*temp
    
  ! linearizing in f'
    temp = bs_b2psifpsib(trial,pst79,lin,pst79,bzt79) &
         + bs_b2psiffb  (trial,pst79,lin,bfpt79,bzt79) &
         + bs_b2psiffb  (trial,pst79,bfpt79,lin,bzt79)

    r_bf = r_bf -          thimp_bf     *dt*temp
    q_bf = q_bf + (3./12. - thimp_bf*bdf)*dt*temp

    if(eqsubtract.eq.1) then

       !(at+a0)(bt+b0)(ct+c0)=(at bt + at b0 + a0 bt + a0 b0)(ct + c0)=
       != (at bt ct + at b0 ct + a0 bt ct + a0 b0 ct) + (at bt c0 + at b0 c0 + a0 bt c0 + a0 b0 c0)

       temp = bs_b2psifpsib(trial,lin,bfpt79,pst79,bz079) & !at bt c0
            + bs_b2psifpsib(trial,lin,bfp079,pst79,bz079) & !a0 bt c0
            + bs_b2psifpsib(trial,lin,bfpt79,ps079,bz079) & !at b0 c0
            + bs_b2psifpsib(trial,lin,bfp079,ps079,bz079)*6 & !a0 b0 c0
            + bs_b2psifpsib(trial,lin,bfp079,pst79,bzt79) & !a0 bt ct
            + bs_b2psifpsib(trial,lin,bfpt79,ps079,bzt79) & !at b0 ct
            + bs_b2psifpsib(trial,lin,bfp079,ps079,bzt79) & !a0 b0 ct
            !--------------------------------------------------------!
            + bs_b2psifpsib(trial,pst79,bfpt79,lin,bz079) & !at bt c0
            + bs_b2psifpsib(trial,ps079,bfpt79,lin,bz079) & !a0 bt c0
            + bs_b2psifpsib(trial,pst79,bfp079,lin,bz079) & !at b0 c0
            + bs_b2psifpsib(trial,ps079,bfp079,lin,bz079)*6 & !a0 b0 c0
            + bs_b2psifpsib(trial,ps079,bfpt79,lin,bzt79) & !a0 bt ct
            + bs_b2psifpsib(trial,pst79,bfp079,lin,bzt79) & !at b0 ct
            + bs_b2psifpsib(trial,ps079,bfp079,lin,bzt79) & !a0 b0 ct
             !--------------------------------------------------------!
            + bs_b2psiffb  (trial,lin,bfpt79,bfpt79,bz079) & !at bt c0
            + bs_b2psiffb  (trial,lin,bfp079,bfpt79,bz079) & !a0 bt c0
            + bs_b2psiffb  (trial,lin,bfpt79,bfp079,bz079)*6 & !at b0 c0
            + bs_b2psiffb  (trial,lin,bfp079,bfp079,bz079) & !a0 b0 c0
            + bs_b2psiffb  (trial,lin,bfp079,bfpt79,bzt79) & !a0 bt ct
            + bs_b2psiffb  (trial,lin,bfpt79,bfp079,bzt79) & !at b0 ct
            + bs_b2psiffb  (trial,lin,bfp079,bfp079,bzt79)   !a0 b0 ct
             !--------------------------------------------------------!
       ddterm(:,psi_g) = ddterm(:,psi_g) + (1./12.)*dt*temp


        temp = bs_b2psifpsib(trial,pst79,bfpt79,ps079,lin) & !at bt c0
             + bs_b2psifpsib(trial,ps079,bfpt79,ps079,lin) & !a0 bt c0
             + bs_b2psifpsib(trial,pst79,bfp079,ps079,lin) & !at b0 c0
             + bs_b2psifpsib(trial,ps079,bfp079,ps079,lin)*6 & !a0 b0 c0
             + bs_b2psifpsib(trial,ps079,bfpt79,pst79,lin) & !a0 bt ct
             + bs_b2psifpsib(trial,pst79,bfp079,pst79,lin) & !at b0 ct
             + bs_b2psifpsib(trial,ps079,bfp079,pst79,lin) & !a0 b0 ct
             !--------------------------------------------------------!
             + bs_b2psiffb  (trial,pst79,bfpt79,bfp079,lin) & !at bt c0
             + bs_b2psiffb  (trial,ps079,bfpt79,bfp079,lin) & !a0 bt c0
             + bs_b2psiffb  (trial,pst79,bfp079,bfp079,lin) & !at b0 c0
             + bs_b2psiffb  (trial,ps079,bfp079,bfp079,lin)*6 & !a0 b0 c0
             + bs_b2psiffb  (trial,ps079,bfpt79,bfpt79,lin) & !a0 bt ct
             + bs_b2psiffb  (trial,pst79,bfp079,bfpt79,lin) & !at b0 ct
             + bs_b2psiffb  (trial,ps079,bfp079,bfpt79,lin)   !a0 b0 ct
       ddterm(:,bz_g) = ddterm(:,bz_g) + (1./12.)*dt*temp



       temp = bs_b2psifpsib(trial,pst79,lin,pst79,bz079) & !at bt c0
            + bs_b2psifpsib(trial,ps079,lin,pst79,bz079) & !a0 bt c0
            + bs_b2psifpsib(trial,pst79,lin,ps079,bz079) & !at b0 c0
            + bs_b2psifpsib(trial,ps079,lin,ps079,bz079)*6 & !a0 b0 c0
            + bs_b2psifpsib(trial,ps079,lin,pst79,bzt79) & !a0 bt ct
            + bs_b2psifpsib(trial,pst79,lin,ps079,bzt79) & !at b0 ct
            + bs_b2psifpsib(trial,ps079,lin,ps079,bzt79) & !a0 b0 ct
            !--------------------------------------------------------!
            + bs_b2psiffb  (trial,pst79,lin,bfpt79,bz079) & !at bt c0
            + bs_b2psiffb  (trial,ps079,lin,bfpt79,bz079) & !a0 bt c0
            + bs_b2psiffb  (trial,pst79,lin,bfp079,bz079) & !at b0 c0
            + bs_b2psiffb  (trial,ps079,lin,bfp079,bz079)*6 & !a0 b0 c0
            + bs_b2psiffb  (trial,ps079,lin,bfpt79,bzt79) & !a0 bt ct
            + bs_b2psiffb  (trial,pst79,lin,bfp079,bzt79) & !at b0 ct
            + bs_b2psiffb  (trial,ps079,lin,bfp079,bzt79) & !a0 b0 ct
            !--------------------------------------------------------!
            + bs_b2psiffb  (trial,pst79,bfpt79,lin,bz079) & !at bt c0
            + bs_b2psiffb  (trial,ps079,bfpt79,lin,bz079) & !a0 bt c0
            + bs_b2psiffb  (trial,pst79,bfp079,lin,bz079) & !at b0 c0
            + bs_b2psiffb  (trial,ps079,bfp079,lin,bz079)*6 & !a0 b0 c0
            + bs_b2psiffb  (trial,ps079,bfpt79,lin,bzt79) & !a0 bt ct
            + bs_b2psiffb  (trial,pst79,bfp079,lin,bzt79) & !at b0 ct
            + bs_b2psiffb  (trial,ps079,bfp079,lin,bzt79)   !a0 b0 ct
       q_bf = q_bf + (1./12.)*dt*temp
    end if
  end subroutine bootstrap_axial_field



!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
! See simplified version of the Bootstrap models below where the terms within jbscommon term 
! that contains psi,f',F are not linearized and computed in a separate subroutine  
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------  

   !calculating bootstrap current
  subroutine calculate_CommonTerm_Lambda(temp1,temp2,tempAA, tempBB, tempCC)
    !using da  /d psi    = 1/|Bp|^2 1/R [(1/R psi_z + f'_r) a_z     + (1/R psi_r - f'_z) a_r]

   !Sauter & Angioni (1999): ibootstrap_model=3: ! equivalent to 1 but a simplified version
    !temp3 = A = -  F L31 p d lnp  /d psi               
    !temp4 = B = -  F L32 Pe d lnTe /d psi              
    !temp5 = C = -  F L34  alpha (p-pe)d lnTi /d psit   
    !p d lnp  /d psi  = 1/|Bp|^2 1/R [(1/R psi_z + f'_r) p_z     + (1/R psi_r - f'_z) p_r]
    !d lnTe /d psi    = 1/|Bp|^2 1/R [(1/R psi_z + f'_r) Te_z/Te + (1/R psi_r - f'_z) Te_r/Te]
    !d lnTi /d psi    = 1/|Bp|^2 1/R [(1/R psi_z + f'_r) Ti_z/Te + (1/R psi_r - f'_z) Ti_r/Te]


   !Redl et al (2021): ibootstrap_model=4:! equivalent to 2 but a simplified version
    !temp3 = A = -F L31        p d lnn  /d psi         
    !temp4 = B = -F (L31+L32) Pe d lnTe /d psi           
    !temp5 = C = -F (L31+L34alpha) (p-pe)d lnTi /d psi  
    !d lnp  /d psi    = 1/|Bp|^2 1/R [(1/R psi_z + f'_r) nt_z/nt    + (1/R psi_r - f'_z) nt_r/nt]
    !d lnTe /d psi    = 1/|Bp|^2 1/R [(1/R psi_z + f'_r) Te_z/Te    + (1/R psi_r - f'_z) Te_r/Te]
    !d lnTi /d psi    = 1/|Bp|^2 1/R [(1/R psi_z + f'_r) Ti_z/Te    + (1/R psi_r - f'_z) Ti_r/Te]

    !tempD = <J.B> = [(A) +  (B) +  (C)]
   
   
   !jbscommon   = -  1 / <B^2>   [ (A) +  (B) +  (C) ]
    !temp1 =  <J.B>/<B^2>  * bootsrap_alpha   =  tempD 1 / <B^2> bootsrap_alpha
    !temp2 =  <J.B>                           =  tempD           bootsrap_alpha 

     use basic
     use m3dc1_nint
   
     implicit none
   
     vectype, dimension(MAX_PTS) :: tempDD, tempAA, tempBB, tempCC, temp1, temp2, iBpsq
   
     temp1 = 0.
     temp2 = 0.

     iBpsq(:) = 1./((-ri_79*pst79(:,OP_DZ)-bfpt79(:,OP_DR))*(-ri_79*pst79(:,OP_DZ)-bfpt79(:,OP_DR)) &
              + (ri_79*pst79(:,OP_DR)-bfpt79(:,OP_DZ))*(ri_79*pst79(:,OP_DR)-bfpt79(:,OP_DZ)))
   
     tempBB = (pst79(:,OP_DZ)*ri_79 + bfpt79(:,OP_DR))*tet79(:,OP_DZ)/tet79(:,OP_1) &
             + (pst79(:,OP_DR)*ri_79 - bfpt79(:,OP_DZ))*tet79(:,OP_DR)/tet79(:,OP_1)
   
     tempCC =  (pst79(:,OP_DZ)*ri_79 + bfpt79(:,OP_DR))*tit79(:,OP_DZ)/tit79(:,OP_1) &
             + (pst79(:,OP_DR)*ri_79 - bfpt79(:,OP_DZ))*tit79(:,OP_DR)/tit79(:,OP_1)
   
     if(ibootstrap_model.eq.1 .or. ibootstrap_model.eq.3)then !Sauter & Angioni (1999) 
         tempAA = (pst79(:,OP_DZ)*ri_79 + bfpt79(:,OP_DR))*pt79(:,OP_DZ) &
                 + (pst79(:,OP_DR)*ri_79 - bfpt79(:,OP_DZ))*pt79(:,OP_DR)
        
         !temp3 = A = -  F L31 p d lnp  /d psit  
         tempAA = iBpsq(:)*ri_79*bzt79(:,OP_1)*jbsl3179(:,OP_1)*tempAA

         !temp4 = B = -F (L31+L32) Pe d lnTe /d psit
         tempBB = iBpsq(:)*ri_79*bzt79(:,OP_1)*jbsl3279(:,OP_1)*pet79(:,OP_1)*(tempBB)

         !temp5 = C = -  F L34  alpha (p-pe)d lnTi /d psit 
         tempCC = iBpsq(:)*ri_79*bzt79(:,OP_1)*jbsl3479(:,OP_1)*jbsalpha79(:,OP_1)*(pt79(:,OP_1)-pet79(:,OP_1))*(tempCC)
       
         tempDD = (tempAA) + (tempBB) +(tempCC)

     else if (ibootstrap_model.eq.2 .or. ibootstrap_model.eq.4)then !Redl et al (2021) 
         tempAA = (pst79(:,OP_DZ)*ri_79 + bfpt79(:,OP_DR))*nt79(:,OP_DZ)/nt79(:,OP_1)*pt79(:,OP_1) &
                  +  (pst79(:,OP_DR)*ri_79 - bfpt79(:,OP_DZ))*nt79(:,OP_DR)/nt79(:,OP_1)*pt79(:,OP_1)

        !temp3 = A = -F L31        p d lnn  /d psi 
        tempAA = iBpsq(:)*ri_79*bzt79(:,OP_1)*jbsl3179(:,OP_1)*tempAA

        !temp4 = B = -F (L31+L32) Pe d lnTe /d psi
        tempBB = iBpsq(:)*ri_79*bzt79(:,OP_1)*(jbsl3179(:,OP_1)+jbsl3279(:,OP_1))*pet79(:,OP_1)*(tempBB)

        !temp5 = C = -F (L31+L34alpha) (p-pe)d lnTi /d psi
        tempCC = iBpsq(:)*ri_79*bzt79(:,OP_1)*(jbsl3179(:,OP_1)+jbsl3479(:,OP_1)*jbsalpha79(:,OP_1))* &
                (pt79(:,OP_1)-pet79(:,OP_1))*(tempCC)
        
        tempDD = (tempAA) + (tempBB) +(tempCC)
     end if
    
         

     if(ibootstrap_model.eq.5)then !Constant Lambda =1/ jbscommon   = -  F / <B^2>   [ L31 (dlnp/dpsi) + L32 Pe (dlnTe/dpsi) + L34 Pe alpha (dlnTi/dpsi) ] = 1
      temp1=1.0
      temp2=1.0
     else !if (ibootstrap_model = 1,2,3,4)
      temp1=tempDD*jbsfluxavg_iBsq_B79(:,OP_1)*bootstrap_alpha
      temp2=tempDD*bootstrap_alpha
      
     endif

   
    end subroutine calculate_CommonTerm_Lambda
   
   
  !calculating bootstrap current
  subroutine calculate_CommonTerm_Lambda_fordtedpsit(temp1,temp2,tempAA, tempBB, tempCC)
  !ibootstrap=2 using dte/dpsit: da/dpsit=da/dTe dTe/dpsit = (del a.del Te)/(|del Te|^2 + chi^2) dTe/dpsit


   !Sauter & Angioni (1999): ibootstrap_model=3: ! equivalent to 1 but a simplified version
    !temp3 = A = -  2piq F L31 p d lnp  /d psit               =  (del p.del Te)/(|del Te|^2 + chi^2) dTe/dpsit
    !temp4 = B = -  2piq F L32 Pe d lnTe /d psit              =  pe/Te  dTe/dpsit 
    !temp5 = C = -  2piq F L34  alpha (p-pe)d lnTi /d psit    =  (p-pe)/Ti (del Ti .del Te)/(|del Te|^2 + chi^2) dTe/dpsit
    

   !Redl et al (2021): ibootstrap_model=4:! equivalent to 2 but a simplified version
    !temp3 = A = -2pi Gbar / (iota - helicity_N) L31        p d lnn  /d psit            =  (ne_s Te_s + ni_s Ti_s)/ne (d lnne / d psit)) = (ne_s Te_s + ni_s Ti_s)/ne (del ne.del Te)/(|del Te|^2 + chi^2) dTe/dpsit
    !temp4 = B = -2pi Gbar / (iota - helicity_N) (L31+L32) Pe d lnTe /d psit            =  pe/Te  dTe/dpsit 
    !temp5 = C = -2pi Gbar / (iota - helicity_N) (L31+L34alpha) (p-pe)d lnTi /d psit    =  (p-pe)/Ti (del Ti .del Te)/(|del Te|^2 + chi^2) dTe/dpsit
    
    !tempD = <J.B> = [(A) +  (B) +  (C)]
   
   
   !jbscommon   = -  1 / <B^2>   [ (A) +  (B) +  (C) ]
    !temp1 =  <J.B>/<B^2>  * bootsrap_alpha   =  tempD 1 / <B^2> bootsrap_alpha
    !temp2 =  <J.B>                           =  tempD           bootsrap_alpha 

   
   
     
    
    use basic
    use m3dc1_nint
  
    implicit none
  
    vectype, dimension(MAX_PTS) :: tempDD, tempAA, tempBB, tempCC, temp1, temp2, iBpsq, temp_delmagTe, temp_delTe
    vectype, dimension(MAX_PTS) :: chisq,const1,adaptive_regularization
    integer :: i
    real(dp):: tempbeta,tempvar

    temp1 = 0.
    temp2 = 0.


    !(del a.del Te)=dadr dTedr + 1/r^2 dadphi dTedphi + dadz dTedz 
    !|del Te|^2= dTedr ^2+ 1/r^2 dTedphi^2 +dTedz^2 
    !(|del Te|^2)
    temp_delmagTe = tet79(:,OP_DR)*tet79(:,OP_DR)+ tet79(:,OP_DZ)*tet79(:,OP_DZ)
#if defined(USE3D) || defined(USECOMPLEX)
        if(itor.eq.1) temp_delmagTe  = temp_delmagTe + tet79(:,OP_DP)*tet79(:,OP_DP)*ri2_79
#endif


! Adaptive regularization term: grows larger when gradients are small
! adaptive_regularization = epsilon / (1.0 + alpha * grad_Te_magnitude)
! Compute the regularized expression
! regularized_expression = (dot(del_p, del_Te)) / (grad_Te_magnitude**2 + adaptive_regularization)

   do i = 1, MAX_PTS
      tempbeta=temp_delmagTe(i)
      chisq(i) = ibootstrap_regular / (1.0 + 1e-2 * tempbeta) 
    enddo
   
    tempBB = 1./(tet79(:,OP_1))
         
    
    tempCC=tit79(:,OP_DR)*tet79(:,OP_DR)+ tit79(:,OP_DZ)*tet79(:,OP_DZ)
#if defined(USE3D) || defined(USECOMPLEX)
    if(itor.eq.1) tempCC = tempCC + tit79(:,OP_DP)*tet79(:,OP_DP)*ri2_79
#endif
    tempCC = tempCC/tit79(:,OP_1)/(temp_delmagTe+chisq)
     
  
    if(ibootstrap_model.eq.1 .or. ibootstrap_model.eq.3)then !Sauter & Angioni (1999) 

      !dlnp  / dpsit    =  (del p.del Te)/(|del Te|^2 + chi^2) dTe/dpsit
      tempAA=pt79(:,OP_DR)*tet79(:,OP_DR)+ pt79(:,OP_DZ)*tet79(:,OP_DZ)
#if defined(USE3D) || defined(USECOMPLEX)
      if(itor.eq.1) tempAA = tempAA + pt79(:,OP_DP)*tet79(:,OP_DP)*ri2_79
#endif
      tempAA = tempAA/(temp_delmagTe+chisq) 

       !dnds_term = -2piq G  L31 ((p dlnp / dpsit)      
       tempAA = jbsfluxavg_G79(:,OP_1)*jbsl3179(:,OP_1)*jbs_dtedpsit79(:,OP_1)*(tempAA)

       !dTeds_term = -2piq G L32 (pe dlnTe / dpsit) 
       tempBB = jbsfluxavg_G79(:,OP_1)*jbsl3279(:,OP_1)*pet79(:,OP_1)*jbs_dtedpsit79(:,OP_1)*(tempBB)

       !dTids_term = -2piq G L34 (pe dlnTi / dpsit) 
       tempCC = jbsfluxavg_G79(:,OP_1)*jbsl3479(:,OP_1)*jbsalpha79(:,OP_1)*(pt79(:,OP_1)-pet79(:,OP_1))* &
               jbs_dtedpsit79(:,OP_1)*(tempCC)
       
       !jdotB = dnds_term + dTeds_term + dTids_term
       tempDD = (tempAA) + (tempBB) + (tempCC)



    else if (ibootstrap_model.eq.2 .or. ibootstrap_model.eq.4)then !Redl et al (2021) 

       !A = d lnn  /d psi    = p/n (del n.del Te)/(|del Te|^2 + chi^2) dTe/dpsit
       !A = pe d lnne  /d psi +  pi d lnni  /d psi    = (tene + tini)/ne (del ne.del Te)/(|del Te|^2 + chi^2) dTe/dpsit
      tempAA=net79(:,OP_DR)*tet79(:,OP_DR)+ net79(:,OP_DZ)*tet79(:,OP_DZ)
#if defined(USE3D) || defined(USECOMPLEX)
      if(itor.eq.1) tempAA = tempAA + net79(:,OP_DP)*tet79(:,OP_DP)*ri2_79
#endif

      tempAA = (net79(:,OP_1)*tet79(:,OP_1)+nt79(:,OP_1)*tit79(:,OP_1))/net79(:,OP_1) * tempAA/(temp_delmagTe+chisq)
        
 !       !dnds_term = -2pi Gbar / (iota - helicity_N)  L31 (ne_s Te_s + ni_s Ti_s)/ne (d lnne / d psit))
        tempAA = jbsfluxavg_G79(:,OP_1)*jbsl3179(:,OP_1)*jbs_dtedpsit79(:,OP_1)*(tempAA)
        
 !       !dTeds_term = -2pi Gbar / (iota - helicity_N) (L31 + L32) pe_s (d lnTe / d psit)
        tempBB = jbsfluxavg_G79(:,OP_1)*(jbsl3179(:,OP_1)+jbsl3279(:,OP_1))*pet79(:,OP_1)*jbs_dtedpsit79(:,OP_1)*(tempBB)
        
 !       !dTids_term = -2pi Gbar / (iota - helicity_N) (L31 + L34 * alpha) pi_s (d lnTi / d psit)
        tempCC = jbsfluxavg_G79(:,OP_1)*(jbsl3179(:,OP_1)+jbsl3479(:,OP_1)*jbsalpha79(:,OP_1))*(pt79(:,OP_1)-pet79(:,OP_1))* &
                jbs_dtedpsit79(:,OP_1)*(tempCC)

 !       !jdotB = dnds_term + dTeds_term + dTids_term
        tempDD = (tempAA) + (tempBB) + (tempCC)
    end if
   
    

    if(ibootstrap_model.eq.5)then 
     temp1=1.0
     temp2=1.0
    else !if (ibootstrap_model = 1,2,3,4)
     temp1=-tempDD*jbsfluxavg_iBsq_B79(:,OP_1)*bootstrap_alpha
     temp2=-tempDD*bootstrap_alpha
    endif

    
    
    
    
   end subroutine calculate_CommonTerm_Lambda_fordtedpsit
   
   
         ! B1psi
    ! =========
    function bs_b1psi(e,h)
      !temp = eta jbscommon ( (-1/r^3) (mu'_z psi_r - mu'_r psi_z) )
     use basic
     use m3dc1_nint
   
     implicit none
     
     vectype, dimension(dofs_per_element) :: bs_b1psi
     vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
     vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: h
     vectype, dimension(dofs_per_element) :: temp
     temp = 0.
#if defined(USECOMPLEX) || defined (USE3D)
         if(jadv.eq.0) then
            temp = 0.
         else 
            if(surface_int) then
             temp = 0.
            else          
              if (ibootstrap.eq.1)then
                call calculate_CommonTerm_Lambda(temp79a,temp79b,temp79c,temp79d,temp79e)
              else if (ibootstrap.eq.2)then
                call calculate_CommonTerm_Lambda_fordtedpsit(temp79a,temp79b,temp79c,temp79d,temp79e)
              endif
             temp = intx5(e(:,:,OP_DRP),ri3_79,h(:,OP_DZ),temp79a,eta79(:,OP_1)) &
                    - intx5(e(:,:,OP_DZP),ri3_79,h(:,OP_DR),temp79a,eta79(:,OP_1))
#ifdef USECOMPLEX
             temp = temp - rfac* &
                     (intx5(e(:,:,OP_DR),ri3_79,h(:,OP_DZ),temp79a,eta79(:,OP_1)) &
                    - intx5(e(:,:,OP_DZ),ri3_79,h(:,OP_DR),temp79a,eta79(:,OP_1))) 
       
#endif
            endif
         endif
     
#endif  
     bs_b1psi = temp
         
         
   end function bs_b1psi
   
   
   ! B1b
   ! =======
   function bs_b1b(e,h)
      !temp = eta jbscommon (1/r^2 F (del* mu ))
     use basic
     use m3dc1_nint
   
     implicit none
     
     vectype, dimension(dofs_per_element) :: bs_b1b
     vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
     vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: h
     vectype, dimension(dofs_per_element) :: temp
   
     temp = 0.
   
     if(jadv.eq.0) then
        if(surface_int) then
           temp = 0.
        else
   
          if (ibootstrap.eq.1)then
            call calculate_CommonTerm_Lambda(temp79a,temp79b,temp79c,temp79d,temp79e)
          else if (ibootstrap.eq.2)then
            call calculate_CommonTerm_Lambda_fordtedpsit(temp79a,temp79b,temp79c,temp79d,temp79e)
          endif
           temp = intx4(e(:,:,OP_1),temp79a,h(:,OP_1),eta79(:,OP_1))
   
        end if
     else 
        if(surface_int) then
           temp = 0.
        else
   
          if (ibootstrap.eq.1)then
            call calculate_CommonTerm_Lambda(temp79a,temp79b,temp79c,temp79d,temp79e)
          else if (ibootstrap.eq.2)then
            call calculate_CommonTerm_Lambda_fordtedpsit(temp79a,temp79b,temp79c,temp79d,temp79e)
          endif
           temp = intx5(e(:,:,OP_GS),ri2_79,temp79a,h(:,OP_1),eta79(:,OP_1))
   
           if(itor.eq.1) then
              temp = temp - 2.*intx5(e(:,:,OP_DR),ri3_79,temp79a,h(:,OP_1),eta79(:,OP_1))
           end if
        endif
     endif
   
     bs_b1b = temp
   end function bs_b1b
   
   
   ! B1f
   ! =======
   function bs_b1f(e,h)
      !temp = eta jbscommon (1/r^2 (f'_z mu'_z + f'_r mu'_r))
      
     use basic
     use m3dc1_nint
   
     implicit none
     
     vectype, dimension(dofs_per_element) :: bs_b1f
     vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
     vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: h
     vectype, dimension(dofs_per_element) :: temp
   
   
     temp = 0.
   
#if defined(USECOMPLEX) || defined(USE3D)
      if(jadv.eq.0) then
        temp = 0.
      else 
        if(surface_int) then
           temp = 0.
        else
   
          if (ibootstrap.eq.1)then
            call calculate_CommonTerm_Lambda(temp79a,temp79b,temp79c,temp79d,temp79e)
          else if (ibootstrap.eq.2)then
            call calculate_CommonTerm_Lambda_fordtedpsit(temp79a,temp79b,temp79c,temp79d,temp79e)
          endif
   
         temp = intx5(e(:,:,OP_DRP),ri2_79,temp79a,h(:,OP_DR),eta79(:,OP_1)) &
              + intx5(e(:,:,OP_DZP),ri2_79,temp79a,h(:,OP_DZ),eta79(:,OP_1))
#ifdef USECOMPLEX
         temp = temp - rfac* &
                (intx5(e(:,:,OP_DR),ri2_79,temp79a,h(:,OP_DR),eta79(:,OP_1)) &
                +intx5(e(:,:,OP_DZ),ri2_79,temp79a,h(:,OP_DZ),eta79(:,OP_1)))
#endif
   
        endif
      endif
#endif
   
     bs_b1f = temp
   end function bs_b1f
   
   
   
       ! B2psi
    ! =========
   function bs_b2psi(e,h)
     !temp = eta jbscommon (1/r^2 (psi_z mu_z + psi_r mu_r))
   
     use basic
     use m3dc1_nint
   
     implicit none
   
     vectype, dimension(dofs_per_element) :: bs_b2psi
     vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
     vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: h
     vectype, dimension(dofs_per_element) :: temp
   
   
     temp = 0.
   
     if (ibootstrap.eq.1)then
      call calculate_CommonTerm_Lambda(temp79a,temp79b,temp79c,temp79d,temp79e)
    else if (ibootstrap.eq.2)then
      call calculate_CommonTerm_Lambda_fordtedpsit(temp79a,temp79b,temp79c,temp79d,temp79e)
    endif
   
     temp = intx5(e(:,:,OP_DZ),ri2_79,h(:,OP_DZ),temp79a,eta79(:,OP_1)) &
          + intx5(e(:,:,OP_DR),ri2_79,h(:,OP_DR),temp79a,eta79(:,OP_1))
   
     bs_b2psi = temp
   end function bs_b2psi
   
   
   
   
   ! B2f
   ! =======
   function bs_b2f(e,h)
   
     !temp = eta jbscommon ( 1/r (mu_z f'_r - mu_r f'_z) )
   
     use basic
     use m3dc1_nint
   
     implicit none
   
     vectype, dimension(dofs_per_element) :: bs_b2f
     vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
     vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: h
     vectype, dimension(dofs_per_element) :: temp
   
     temp = 0.
   
#if defined(USECOMPLEX) || defined(USE3D)
  if (ibootstrap.eq.1)then
    call calculate_CommonTerm_Lambda(temp79a,temp79b,temp79c,temp79d,temp79e)
  else if (ibootstrap.eq.2)then
    call calculate_CommonTerm_Lambda_fordtedpsit(temp79a,temp79b,temp79c,temp79d,temp79e)
  endif
   
     temp = intx5(e(:,:,OP_DZ),ri_79,h(:,OP_DR),temp79a,eta79(:,OP_1)) &
          - intx5(e(:,:,OP_DR),ri_79,h(:,OP_DZ),temp79a,eta79(:,OP_1))
#else 
     temp = 0.
#endif 
     bs_b2f = temp
   
   
   end function bs_b2f
   
   
   
   
   
   subroutine bootstrap_flux_simplified(trial, lin, ssterm, ddterm, r_bf, q_bf, thimpf, thimp_bf)
     use basic
     use arrays
     use m3dc1_nint
   
     implicit none
   
     vectype, dimension(dofs_per_element, MAX_PTS, OP_NUM), intent(in) :: trial
     vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: lin 
     vectype, dimension(dofs_per_element, num_fields), intent(inout) :: ssterm, ddterm
     vectype, dimension(dofs_per_element), intent(out) :: r_bf, q_bf
     real, intent(in) :: thimpf, thimp_bf
   
     vectype, dimension(dofs_per_element) :: temp
     
     if(numvar.eq.1) then
        ! linearizing in psi  
        temp = bs_b1psi(trial,lin) 
        
         ssterm(:,psi_g) = ssterm(:,psi_g) -          thimpf     *dt*temp
         ddterm(:,psi_g) = ddterm(:,psi_g) + (1. - thimpf*bdf)*dt*temp
        
     else     
         ! linearizing in psi
        temp = bs_b1psi(trial,lin) 
   
        ssterm(:,psi_g) = ssterm(:,psi_g) -          thimpf     *dt*temp
        ddterm(:,psi_g) = ddterm(:,psi_g) + (1.    - thimpf*bdf)*dt*temp
   
          ! linearizing in b  
        temp = bs_b1b(trial,lin) 
        
        ssterm(:,bz_g) = ssterm(:,bz_g) -          thimpf     *dt*temp
        ddterm(:,bz_g) = ddterm(:,bz_g) + (1.    - thimpf*bdf)*dt*temp
        
   
      ! linearizing in f'  
        temp = bs_b1f(trial,lin) 
        r_bf = r_bf -          thimp_bf     *dt*temp
        q_bf = q_bf + (1. - thimp_bf*bdf)*dt*temp            
        
        if(eqsubtract.eq.1) then
         !at=(a1+a0)
         !temp = bs_b1psi(trial,lin) 
         !--------------------------------------------------------!
         !ddterm(:,psi_g) = ddterm(:,psi_g) + (1.)*dt*temp      
           
   
         !temp = bs_b1b(trial,lin) 
         !--------------------------------------------------------!  
         !ddterm(:,bz_g) = ddterm(:,bz_g) + (1.)*dt*temp    
   
         !temp = bs_b1f(trial,lin) 
         !--------------------------------------------------------!      
         !q_bf = q_bf + (1.)*dt*temp    
        end if
     end if
   end subroutine bootstrap_flux_simplified
   
   
   
   
   
   
   
    subroutine bootstrap_axial_field_simplified(trial, lin, ssterm, ddterm, r_bf, q_bf, thimpf, thimp_bf)
     use basic
     use arrays
     use m3dc1_nint
   
     implicit none
   
     vectype, dimension(dofs_per_element, MAX_PTS, OP_NUM), intent(in) :: trial
     vectype, dimension(MAX_PTS, OP_NUM), intent(in) :: lin
     vectype, dimension(dofs_per_element, num_fields), intent(inout) :: ssterm, ddterm
     vectype, dimension(dofs_per_element), intent(out) :: r_bf, q_bf
     real, intent(in) :: thimpf, thimp_bf
   
     vectype, dimension(dofs_per_element) :: temp
   
   
   ! linearizing in psi
     temp = bs_b2psi(trial,lin)
          
     ssterm(:,psi_g) = ssterm(:,psi_g) -          thimpf     *dt*temp
     ddterm(:,psi_g) = ddterm(:,psi_g) + (1.    - thimpf*bdf)*dt*temp
     
   
   ! linearizing in F - No F term
   !  temp = bs_b2b(trial,lin)
   !  ssterm(:,bz_g) = ssterm(:,bz_g) -          thimpf     *dt*temp
   !  ddterm(:,bz_g) = ddterm(:,bz_g) + (1.    - thimpf*bdf)*dt*temp
     
   ! linearizing in f'
     temp = bs_b2f(trial,lin)
   
     r_bf = r_bf -          thimp_bf     *dt*temp
     q_bf = q_bf + (1. - thimp_bf*bdf)*dt*temp
   
     if(eqsubtract.eq.1) then
   
        !at=a1+a0
   
      !  temp = bs_b2psi(trial,lin) 
      !  ddterm(:,psi_g) = ddterm(:,psi_g) + (1.)*dt*temp
   
      ! linearizing in F - No F term
      !  temp = bs_b2b(trial,lin) 
      !  ddterm(:,bz_g) = ddterm(:,bz_g) + (1.)*dt*temp
   
   
   
      !  temp = bs_b2f(trial,lin) 
      !  q_bf = q_bf + (1..)*dt*temp
     end if
   end subroutine bootstrap_axial_field_simplified


end module bootstrap


