
module temperature_plots

contains

subroutine advection(o)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(MAX_PTS), intent(out) :: o
 ! Advection
  ! ~~~~~~~~~~~~~~~~~~
     o =  - r_79*tet79(:,OP_DZ)*nt79(:,OP_1)*pht79(:,OP_DR) &
                + r_79*tet79(:,OP_DR)*nt79(:,OP_1)*pht79(:,OP_DZ)    
         
     if(itor.eq.1) then
        o = o + 2.*(gam-1.)*tet79(:,OP_1)*nt79(:,OP_1)*pht79(:,OP_DZ)
     endif

#if defined(USE3D) || defined(USECOMPLEX)
     if(numvar.ge.2) then
        o = o - tet79(:,OP_DP)*nt79(:,OP_1)*vzt79(:,OP_1) &
         - (gam-1.)*tet79(:,OP_1)*nt79(:,OP_1)*vzt79(:,OP_DP)
     end if
#endif

     if(numvar.ge.3) then
        temp79b = -(tet79(:,OP_DR)*cht79(:,OP_DR) &
                  + tet79(:,OP_DZ)*cht79(:,OP_DZ))  &
                  - (gam-1.)*tet79(:,OP_1)*cht79(:,OP_GS)
     
        o = o + ri2_79*temp79b*nt79(:,OP_1)
     end if

  ! Electron Temp Advection
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(db.ne.0.) then
#if defined(USE3D) || defined(USECOMPLEX)
  temp79b = ri2_79*tet79(:,OP_DZ)*pst79(:,OP_DZP)*ni79(:,OP_1) &
       + ri2_79*tet79(:,OP_DR)*pst79(:,OP_DRP)*ni79(:,OP_1) &
       - ri2_79*tet79(:,OP_DP)*pst79(:,OP_GS)*ni79(:,OP_1) &
       + gam* &
       (ri2_79*tet79(:,OP_1)*pst79(:,OP_DZP)*ni79(:,OP_DZ) &
       +ri2_79*tet79(:,OP_1)*pst79(:,OP_DRP)*ni79(:,OP_DR) &
       -ri2_79*tet79(:,OP_1)*pst79(:,OP_GS)*ni79(:,OP_DP))
  o = o + temp79b*db
#endif
        if(numvar.ge.2) then
           temp79b =  ri_79*pet79(:,OP_DZ)*bzt79(:,OP_DR)*ni79(:,OP_1) &
                 - ri_79*pet79(:,OP_DR)*bzt79(:,OP_DZ)*ni79(:,OP_1) &
                 + gam* &
                  (ri_79*pet79(:,OP_1)*bzt79(:,OP_DR)*ni79(:,OP_DZ) &
                 - ri_79*pet79(:,OP_1)*bzt79(:,OP_DZ)*ni79(:,OP_DR))
           o = o + temp79b*db
        end if
 
        if(i3d.eq.1) then
           if(numvar.ge.3) then
#if defined(USE3D) || defined(USECOMPLEX)
              o = o +  ri_79*pet79(:,OP_DZ)*bfp179(:,OP_DRP)*ni79(:,OP_1) &
                    - ri_79*pet79(:,OP_DR)*bfp179(:,OP_DZP)*ni79(:,OP_1) &
                  + gam* &
                    ( ri_79*pet79(:,OP_1)*bfp179(:,OP_DRP)*ni79(:,OP_DZ) &
                    - ri_79*pet79(:,OP_1)*bfp179(:,OP_DZP)*ni79(:,OP_DR))
#endif
             
           end if
       end if
  endif
 o = o * (1. + (1.-pefac)/pefac)
end subroutine advection
subroutine advection1(o)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(MAX_PTS), intent(out) :: o
 ! Advection
  ! ~~~~~~~~~~~~~~~~~~
     o = 0.
     o =  o - r_79*tet79(:,OP_DZ)*nt79(:,OP_1)*pht79(:,OP_DR) &
                + r_79*tet79(:,OP_DR)*nt79(:,OP_1)*pht79(:,OP_DZ)    
         
     if(itor.eq.1) then
        o = o + 2.*(gam-1.)*tet79(:,OP_1)*nt79(:,OP_1)*pht79(:,OP_DZ)
     endif


 o = o * (1. + (1.-pefac)/pefac)
end subroutine advection1

subroutine advection2(o)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(MAX_PTS), intent(out) :: o
 ! Advection
  ! ~~~~~~~~~~~~~~~~~~
     o = 0.

#if defined(USE3D) || defined(USECOMPLEX)
     if(numvar.ge.2) then
        o = o - tet79(:,OP_DP)*nt79(:,OP_1)*vzt79(:,OP_1) &
         - (gam-1.)*tet79(:,OP_1)*nt79(:,OP_1)*vzt79(:,OP_DP)
     end if
#endif



 o = o * (1. + (1.-pefac)/pefac)
end subroutine advection2

subroutine advection3(o)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(MAX_PTS), intent(out) :: o
 ! Advection
  ! ~~~~~~~~~~~~~~~~~~
     o = 0.

     if(numvar.ge.3) then
        temp79b = -(tet79(:,OP_DR)*cht79(:,OP_DR) &
                  + tet79(:,OP_DZ)*cht79(:,OP_DZ))  &
                  - (gam-1.)*tet79(:,OP_1)*cht79(:,OP_GS)
     
        o = o + ri2_79*temp79b*nt79(:,OP_1)
     end if

 o = o * (1. + (1.-pefac)/pefac)
end subroutine advection3

subroutine hf_perp(dofs)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element), intent(out) :: dofs

  ! Perpendicular Heat Flux
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  dofs =   -(gam-1)*    &
          ( intx3(mu79(:,:,OP_DZ),tet79(:,OP_DZ),kap79(:,OP_1)) &
          + intx3(mu79(:,:,OP_DR),tet79(:,OP_DR),kap79(:,OP_1)))
#if defined(USE3D) || defined(USECOMPLEX)
  dofs = dofs + &
       (gam-1)*intx4(mu79(:,:,OP_1),ri2_79,tet79(:,OP_DPP),kap79(:,OP_1))
#endif
  dofs = (1.+kappai_fac*(1.-pefac)/pefac)*dofs

end subroutine hf_perp

subroutine hf_par(dofs)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element), intent(out) :: dofs
!
  temp79b = kar79(:,OP_1)*ri2_79*b2i79(:,OP_1)

  dofs = (gam-1.)* &
       (intx5(mu79(:,:,OP_DZ),temp79b,pstx79(:,OP_DR),pstx79(:,OP_DZ),tet79(:,OP_DR)) &
       -intx5(mu79(:,:,OP_DR),temp79b,pstx79(:,OP_DZ),pstx79(:,OP_DZ),tet79(:,OP_DR)) &
       -intx5(mu79(:,:,OP_DZ),temp79b,pstx79(:,OP_DR),pstx79(:,OP_DR),tet79(:,OP_DZ)) &
       +intx5(mu79(:,:,OP_DR),temp79b,pstx79(:,OP_DZ),pstx79(:,OP_DR),tet79(:,OP_DZ)))

#if defined(USE3D) || defined(USECOMPLEX)
       temp79b = -kar79(:,OP_1)*ri3_79*bztx79(:,OP_1)*b2i79(:,OP_1)
       temp79c = pstx79(:,OP_DR)*tet79(:,OP_DZ) &
               - pstx79(:,OP_DZ)*tet79(:,OP_DR)
       temp79d = temp79c*bztx79(:,OP_1 )*b2i79(:,OP_1 )*kar79(:,OP_1)

       dofs = dofs + (gam - 1.)* &
            (intx4(mu79(:,:,OP_DZ),temp79b,pstx79(:,OP_DR),tet79(:,OP_DP)) &
            -intx4(mu79(:,:,OP_DR),temp79b,pstx79(:,OP_DZ),tet79(:,OP_DP)))
       
       dofs = dofs - (gam - 1.)*intx3(mu79(:,:,OP_DP),ri3_79,temp79d)

       temp79b = bztx79(:,OP_1)*bztx79(:,OP_1 )   &
            *tet79(:,OP_DP)*b2i79(:,OP_1 )*kar79(:,OP_1 )
       dofs = dofs - (gam - 1.)*intx3(mu79(:,:,OP_DP),ri4_79,temp79b)


       if(i3d.eq.1 .and. numvar.ge.2) then
          temp79b = kar79(:,OP_1)*ri_79*b2i79(:,OP_1)
          
          dofs = dofs + (gam -1.)*&
               (intx5(mu79(:,:,OP_DZ),pstx79(:,OP_DR),temp79b,bfptx79(:,OP_DZ),tet79(:,OP_DZ)) &
               -intx5(mu79(:,:,OP_DR),pstx79(:,OP_DZ),temp79b,bfptx79(:,OP_DZ),tet79(:,OP_DZ)) &
               +intx5(mu79(:,:,OP_DZ),pstx79(:,OP_DR),temp79b,bfptx79(:,OP_DR),tet79(:,OP_DR)) &
               -intx5(mu79(:,:,OP_DR),pstx79(:,OP_DZ),temp79b,bfptx79(:,OP_DR),tet79(:,OP_DR)) &
               +intx5(mu79(:,:,OP_DZ),bfptx79(:,OP_DZ),temp79b,pstx79(:,OP_DR ),tet79(:,OP_DZ)) &
               +intx5(mu79(:,:,OP_DR),bfptx79(:,OP_DR),temp79b,pstx79(:,OP_DR ),tet79(:,OP_DZ)) &
               -intx5(mu79(:,:,OP_DZ),bfptx79(:,OP_DZ),temp79b,pstx79(:,OP_DZ ),tet79(:,OP_DR)) &
               -intx5(mu79(:,:,OP_DR),bfptx79(:,OP_DR),temp79b,pstx79(:,OP_DZ ),tet79(:,OP_DR)))
          
          temp79b = kar79(:,OP_1)*ri2_79*bztx79(:,OP_1)*b2i79(:,OP_1)* &
               tet79(:,OP_DP)
          temp79c = bfptx79(:,OP_DZ)*tet79(:,OP_DZ) &
               + bfptx79(:,OP_DR)*tet79(:,OP_DR)          
          temp79d = temp79c*bztx79(:,OP_1 )*b2i79(:,OP_1 )*kar79(:,OP_1 )
          
          dofs = dofs + (gam - 1.)* &
               (intx3(mu79(:,:,OP_DZ),bfptx79(:,OP_DZ),temp79b) &
               +intx3(mu79(:,:,OP_DR),bfptx79(:,OP_DR),temp79b) &
               +intx3(mu79(:,:,OP_DP),ri2_79,temp79d))
          
          temp79b = - kar79(:,OP_1)*b2i79(:,OP_1)

          dofs = dofs + (gam - 1.)* &
               (intx5(mu79(:,:,OP_DZ),bfptx79(:,OP_DZ),temp79b,bfptx79(:,OP_DZ),tet79(:,OP_DZ)) &
               +intx5(mu79(:,:,OP_DR),bfptx79(:,OP_DR),temp79b,bfptx79(:,OP_DZ),tet79(:,OP_DZ)) &
               +intx5(mu79(:,:,OP_DZ),bfptx79(:,OP_DZ),temp79b,bfptx79(:,OP_DR),tet79(:,OP_DR)) &
               +intx5(mu79(:,:,OP_DR),bfptx79(:,OP_DR),temp79b,bfptx79(:,OP_DR),tet79(:,OP_DR)))
       endif
#endif

end subroutine hf_par

subroutine ohmic(o)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(MAX_PTS), intent(out) :: o

  ! Ohmic Heating
  ! ~~~~~~~~~~~~~
      o = (gam-1.)* &
              ri2_79*pst79(:,OP_GS)* pst79(:,OP_GS)* eta79(:,OP_1)   
#if defined(USE3D) || defined(USECOMPLEX)
       o = o + (gam-1)*   &
           (ri4_79*pst79(:,OP_DRP)*pst79(:,OP_DRP)*eta79(:,OP_1)   &
         +  ri4_79*pst79(:,OP_DZP)*pst79(:,OP_DZP)*eta79(:,OP_1))
#endif

       if(numvar.ge.2) then
         
         o = o + (gam-1.)* &
              (ri2_79*bzt79(:,OP_DZ)*bzt79(:,OP_DZ)*eta79(:,OP_1) &
              +ri2_79*bzt79(:,OP_DR)*bzt79(:,OP_DR)*eta79(:,OP_1))

#if defined(USE3D) || defined(USECOMPLEX)
         o = o + 2.*(gam-1.)* &
              (ri3_79*pst79(:,OP_DZP)*bzt79(:,OP_DR)*eta79(:,OP_1)  &
              -ri3_79*pst79(:,OP_DRP)*bzt79(:,OP_DZ)*eta79(:,OP_1))

          if(i3d .eq. 1) then

             o = o +  2.*(gam-1.)* &
             (ri3_79*pst79(:,OP_DZP)*bfpt79(:,OP_DRP)*eta79(:,OP_1)  &
             -ri3_79*pst79(:,OP_DRP)*bfpt79(:,OP_DZP)*eta79(:,OP_1))

             o = o + (gam-1.)* &
             (ri2_79*bzt79(:,OP_DZ)*bfpt79(:,OP_DZP)*eta79(:,OP_1) &
             +ri2_79*bzt79(:,OP_DR)*bfpt79(:,OP_DRP)*eta79(:,OP_1))

             o = o + (gam-1.)* &
               (ri2_79*bfpt79(:,OP_DZP)*bfpt79(:,OP_DZP)*eta79(:,OP_1) &
             +  ri2_79*bfpt79(:,OP_DRP)*bfpt79(:,OP_DRP)*eta79(:,OP_1))      
          endif
#endif
       end if

! ADD heat source if present
  if(heat_source) o = o + (gam-1.)*q79(:,OP_1)

! ADD rad source if present
  if(rad_source) o = o + (gam-1.)*totrad79(:,OP_1)


end subroutine ohmic
subroutine vpar_get(o)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(MAX_PTS), intent(out) :: o


!      bsquared
  temp79c = (pst79(:,OP_DR)**2 + pst79(:,OP_DZ)**2)*ri2_79    &
          + bzt79(:,OP_1)**2*ri2_79      
#if defined(USE3D) || defined(USECOMPLEX)
  temp79c = temp79c     &
         + 2.*ri_79*(bfpt79(:,OP_DR)*pst79(:,OP_DZ) - bfpt79(:,OP_DZ)*pst79(:,OP_DR))  &
         + bfpt79(:,OP_DR)**2 + bfpt79(:,OP_DZ)**2
#endif

!      v dot b
  temp79b = pht79(:,OP_DR)*pst79(:,OP_DR) +  pht79(:,OP_DZ)*pst79(:,OP_DZ)   &
          + bzt79(:,OP_1)*vzt79(:,OP_1)   &
          + ri3_79*(pst79(:,OP_DR)*cht79(:,OP_DZ) - pst79(:,OP_DZ)*cht79(:,OP_DR))  
 
#if defined(USE3D) || defined(USECOMPLEX)
  temp79b = temp79b     &
          + r_79*(bfpt79(:,OP_DR)*pht79(:,OP_DZ) - bfpt79(:,OP_DZ)*pht79(:,OP_DR))  &
          - ri2_79*( cht79(:,OP_DR)*bfpt79(:,OP_DR) +  cht79(:,OP_DZ)*bfpt79(:,OP_DZ))
#endif
  
  o = temp79b/sqrt(temp79c)

end subroutine vpar_get

subroutine f1vplot_sub(term)
  use basic
  use m3dc1_nint
  use metricterms_new
  
  implicit none
  vectype, intent(out), dimension(dofs_per_element) :: term
  vectype, dimension(dofs_per_element) :: temp

  temp = b2psiv(mu79,pst79,vzt79)
  term = temp

  temp = b2bu  (mu79,bzt79,pht79)  &
       + b2bchi(mu79,bzt79,cht79)
  term = term + temp

end subroutine f1vplot_sub

subroutine f1eplot_sub(term)
  use basic
  use m3dc1_nint
  use metricterms_new
  
  implicit none
  vectype, intent(out), dimension(dofs_per_element) :: term
  vectype, dimension(dofs_per_element) :: temp

  ! Resistive and Hyper Terms
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~
  temp = b2psieta(mu79,pst79,eta79)
  term = temp
  
  temp = b2beta(mu79,bzt79,eta79,vz079)
  term = term + temp

  if(i3d.eq.1) then
     temp = b2feta(mu79,bfpt79,eta79)
     term = term + temp
  end if

end subroutine f1eplot_sub

subroutine f2vplot_sub(term)
  use basic
  use m3dc1_nint
  use metricterms_new
  
  implicit none
  vectype, intent(out), dimension(dofs_per_element) :: term
  integer :: jadvs
  vectype, dimension(dofs_per_element) :: temp

  jadvs = jadv
  jadv = 1  ! only for evaluation of these functions

  temp = b1psiu  (mu79,pst79,pht79) &
       + b1psiv  (mu79,pst79,vzt79) &
       + b1psichi(mu79,pst79,cht79)
  term = temp
  
  temp =  b1bu  (mu79,bzt79,pht79) &
       + b1bv  (mu79,bzt79,vzt79) &
       + b1bchi(mu79,bzt79,cht79)
  term = term + temp
  
  if(i3d.eq.1 .and. numvar.ge.2) then
     temp = b1fu  (mu79,bfpt79,pht79)  &
          + b1fv  (mu79,bfpt79,vzt79)  &
          + b1fchi(mu79,bfpt79,cht79)
     term = term + temp
  end if
  jadv = jadvs

end subroutine f2vplot_sub

subroutine f2eplot_sub(term)
  use basic
  use m3dc1_nint
  use metricterms_new
  
  implicit none
  vectype, intent(out), dimension(dofs_per_element) :: term
  vectype, dimension(dofs_per_element) :: temp
  vectype, dimension(MAX_PTS, OP_NUM) ::  siw79

  if(iadiabat.eq.1) then
    siw79 = sie79 + sii79*(1.-pefac)/pefac
    temp =  t3tndenm(mu79,tet79,net79,denm79)  &
         +  t3ts(mu79,tet79,siw79)
  endif
  term = temp
end subroutine f2eplot_sub

subroutine f3vplot_sub(term)
  use basic
  use m3dc1_nint
  use metricterms_new
  
  implicit none
  vectype, intent(out), dimension(dofs_per_element) :: term
  vectype, dimension(dofs_per_element) :: temp

  temp = (gam-1)*(q_delta1(mu79,tit79)-q_delta1(mu79,tet79))
  term = temp
end subroutine f3vplot_sub

subroutine f3eplot_sub(term)
  use basic
  use m3dc1_nint
  use metricterms_new
  
  implicit none
  vectype, intent(out), dimension(dofs_per_element) :: term

  vectype, dimension(dofs_per_element) :: tempx

  ! Ohmic Heating
  ! ~~~~~~~~~~~~~
  tempx = b3psipsieta(mu79,pst79,pst79,eta79,nre179) &
       +  b3bbeta    (mu79,bzt79,bzt79,eta79) &
       +  b3psibeta  (mu79,pst79,bzt79,eta79,nre179) 
  term = tempx

  if(i3d .eq. 1) then
     tempx = b3psifeta(mu79,pst79,bfpt79,eta79,nre179) &
          +  b3bfeta  (mu79,bzt79,bfpt79,eta79,nre179) &
          +  b3ffeta  (mu79,bfpt79,bfpt79,eta79,nre179)   
     term = term + tempx
  endif

  ! Perpendicular Heat Flux
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  tempx = b3tekappa(mu79,tet79,kap79,vz079)
  term = term + tempx


  ! Source terms
  ! ~~~~~~~~~~~~
  tempx = (gam-1.)*b3q(mu79,q79)
  term = term + tempx


  ! Parallel Heat Flux
  ! ~~~~~~~~~~~~~~~~~~
  if(kappar.ne.0. .or. ikapparfunc.eq.2) then

     tempx = tepsipsikappar(mu79,pstx79,pstx79,tet79,b2i79,kar79) &
          + tepsibkappar  (mu79,pstx79,bztx79,tet79,b2i79,kar79) &
          + tebbkappar    (mu79,bztx79,bztx79,tet79,b2i79,kar79)
     term = term + tempx
     if(i3d.eq.1 .and. numvar.ge.2) then
        tempx = tepsifkappar(mu79,pstx79,bfptx79,tet79,b2i79,kar79) &
             + tebfkappar  (mu79,bztx79,bfptx79,tet79,b2i79,kar79) &
             + teffkappar  (mu79,bfptx79,bfptx79,tet79,b2i79,kar79)
        term = term + tempx
     endif
  endif
end subroutine f3eplot_sub

subroutine potential2(dofs)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(dofs_per_element), intent(out) :: dofs

  ! Electrical Potential due to resistive terms
  ! ~~~~~~~~~~~~~~~~~~~~~~~
  dofs = - ( intx4(mu79(:,:,OP_DZ),eta79(:,OP_1), ri_79,bzt79(:,OP_DR)) &
           - intx4(mu79(:,:,OP_DR),eta79(:,OP_1), ri_79,bzt79(:,OP_DZ)) )
#if defined(USE3D) || defined(USECOMPLEX)
  dofs = dofs - ( intx4(mu79(:,:,OP_DZ),eta79(:,OP_1), ri_79,bfpt79(:,OP_DRP)) &
                - intx4(mu79(:,:,OP_DR),eta79(:,OP_1), ri_79,bfpt79(:,OP_DZP)) &
                + intx4(mu79(:,:,OP_DR),eta79(:,OP_1),ri3_79,pst79(:,OP_DRP)) &
                + intx4(mu79(:,:,OP_DZ),eta79(:,OP_1),ri3_79,pst79(:,OP_DZP)) )
#endif
!  if(izone.ne.1) return
  dofs = dofs - intx3(mu79(:,:,OP_DR),pht79(:,OP_DR),bzt79(:,OP_1)) &
              - intx3(mu79(:,:,OP_DZ),pht79(:,OP_DZ),bzt79(:,OP_1)) &
              + intx3(mu79(:,:,OP_DR),pst79(:,OP_DR),vzt79(:,OP_1)) &       
              + intx3(mu79(:,:,OP_DZ),pst79(:,OP_DZ),vzt79(:,OP_1)) &
              + intx4(mu79(:,:,OP_DZ),bzt79(:,OP_1), ri3_79,cht79(:,OP_DR)) &
              - intx4(mu79(:,:,OP_DR),bzt79(:,OP_1), ri3_79,cht79(:,OP_DZ))

#if defined(USE3D) || defined(USECOMPLEX)
  dofs = dofs + intx4(mu79(:,:,OP_DZ),vzt79(:,OP_1), r_79,bfpt79(:,OP_DR)) &
       	      -	intx4(mu79(:,:,OP_DR),vzt79(:,OP_1), r_79,bfpt79(:,OP_DZ))

#endif






 

end subroutine potential2

end module temperature_plots


 

 
