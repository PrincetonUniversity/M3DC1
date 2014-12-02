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
              o = o +  ri_79*pet79(:,OP_DZ)*bf179(:,OP_DRPP)*ni79(:,OP_1) &
                    - ri_79*pet79(:,OP_DR)*bf179(:,OP_DZPP)*ni79(:,OP_1) &
                  + gam* &
                    ( ri_79*pet79(:,OP_1)*bf179(:,OP_DRPP)*ni79(:,OP_DZ) &
                    - ri_79*pet79(:,OP_1)*bf179(:,OP_DZPP)*ni79(:,OP_DR))
#endif
             
           end if
       end if
  endif
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

end subroutine advection3

subroutine hf_perp(i,o)
  use basic
  use m3dc1_nint

  implicit none

  integer, intent(in) :: i
  vectype, dimension(MAX_PTS), intent(out) :: o


     ! Perpendicular Heat Flux
     ! ~~~~~~~~~~~~~~~~~~~~~~~
     o =   -(gam-1)*    &
          ( mu79(:,OP_DZ,i)*tet79(:,OP_DZ)*kap79(:,OP_1) &
          + mu79(:,OP_DR,i)*tet79(:,OP_DR)*kap79(:,OP_1) )
#if defined(USE3D) || defined(USECOMPLEX)
     o = o + (gam-1)*ri2_79*mu79(:,OP_1,i)*tet79(:,OP_DPP)*kap79(:,OP_1)
#endif

end subroutine hf_perp

subroutine hf_par(i,o)
  use basic
  use m3dc1_nint

  implicit none

  integer, intent(in) :: i
  vectype, dimension(MAX_PTS), intent(out) :: o
!
          temp79b = kar79(:,OP_1)*ri2_79*        &
                  (mu79(:,OP_DZ,i)*pstx79(:,OP_DR) &
                 - mu79(:,OP_DR,i)*pstx79(:,OP_DZ))*b2i79(:,OP_1)

          o = (gam-1.)*(temp79b*pstx79(:,OP_DZ)*tet79(:,OP_DR) &
                         - temp79b*pstx79(:,OP_DR)*tet79(:,OP_DZ))

#if defined(USE3D) || defined(USECOMPLEX)
          temp79b = -kar79(:,OP_1)*ri3_79*bztx79(:,OP_1)* &
                    (mu79(:,OP_DZ,i)*pstx79(:,OP_DR) &
                   - mu79(:,OP_DR,i)*pstx79(:,OP_DZ))*b2i79(:,OP_1)
          temp79c = pstx79(:,OP_DR)*tet79(:,OP_DZ) &
                  - pstx79(:,OP_DZ)*tet79(:,OP_DR)
          temp79d = temp79c*bztx79(:,OP_1 )*b2i79(:,OP_1 )*kar79(:,OP_1)
          o = o + (gam - 1.)* (temp79b*tet79(:,OP_DP) &
                                - ri3_79*mu79(:,OP_DP,i)*temp79d)

          temp79b = bztx79(:,OP_1)*bztx79(:,OP_1 )   &
                                *tet79(:,OP_DP)*b2i79(:,OP_1 )*kar79(:,OP_1 )
          o = o -(gam-1)*(ri4_79*mu79(:,OP_DP,i)*temp79b)


          if(i3d.eq.1 .and. numvar.ge.2) then
             temp79b = kar79(:,OP_1)*ri_79* &
                     (mu79(:,OP_DZ,i)*pstx79(:,OP_DR) &
                    - mu79(:,OP_DR,i)*pstx79(:,OP_DZ))*b2i79(:,OP_1)
             temp79c = kar79(:,OP_1)*ri_79* &
                     (mu79(:,OP_DZ,i)*bftx79(:,OP_DZP) &
                    + mu79(:,OP_DR,i)*bftx79(:,OP_DRP))*b2i79(:,OP_1)

             o = o + (gam -1.)*(temp79b*bftx79(:,OP_DZP)*tet79(:,OP_DZ) &
                             + temp79b*bftx79(:,OP_DRP)*tet79(:,OP_DR) &
                             + temp79c*pstx79(:,OP_DR )*tet79(:,OP_DZ) &
                             - temp79c*pstx79(:,OP_DZ )*tet79(:,OP_DR))
       
              temp79b = kar79(:,OP_1)*ri2_79*bztx79(:,OP_1)* &
                    (mu79(:,OP_DZ,i)*bftx79(:,OP_DZP) &
                   + mu79(:,OP_DR,i)*bftx79(:,OP_DRP))*b2i79(:,OP_1)

              temp79c = bftx79(:,OP_DZP)*tet79(:,OP_DZ) &
                      + bftx79(:,OP_DRP)*tet79(:,OP_DR)

              temp79d = temp79c*bztx79(:,OP_1 )*b2i79(:,OP_1 )*kar79(:,OP_1 )

              o = o + (gam - 1.)*(temp79b*tet79(:,OP_DP) &
                                   + ri2_79*mu79(:,OP_DP,i)*temp79d)
           
              temp79b = - kar79(:,OP_1)*              &
                    (mu79(:,OP_DZ,i)*bftx79(:,OP_DZP) &
                   + mu79(:,OP_DR,i)*bftx79(:,OP_DRP))*b2i79(:,OP_1)

              o = o + (gam - 1.)*(temp79b*bftx79(:,OP_DZP)*tet79(:,OP_DZ) &
                            + temp79b*bftx79(:,OP_DRP)*tet79(:,OP_DR))
          endif

#endif

end subroutine hf_par

subroutine ohmic(o)
  use basic
  use m3dc1_nint

  implicit none

  vectype, dimension(MAX_PTS), intent(out) :: o
  real :: ohfac
  ohfac = 1.
  if(ipres.eq.0) ohfac = 0.5

  ! Ohmic Heating
  ! ~~~~~~~~~~~~~
      o = (gam-1.)*ohfac* &
              ri2_79*pst79(:,OP_GS)* pst79(:,OP_GS)* eta79(:,OP_1)   
#if defined(USE3D) || defined(USECOMPLEX)
       o = o + (gam-1)*ohfac*   &
           (ri4_79*pst79(:,OP_DRP)*pst79(:,OP_DRP)*eta79(:,OP_1)   &
         +  ri4_79*pst79(:,OP_DZP)*pst79(:,OP_DZP)*eta79(:,OP_1))
#endif

       if(numvar.ge.2) then
         
         o = o + (gam-1.)*ohfac* &
              (ri2_79*bzt79(:,OP_DZ)*bzt79(:,OP_DZ)*eta79(:,OP_1) &
              +ri2_79*bzt79(:,OP_DR)*bzt79(:,OP_DR)*eta79(:,OP_1))

#if defined(USE3D) || defined(USECOMPLEX)
         o = o + 2.*(gam-1.)*ohfac* &
              (ri3_79*pst79(:,OP_DZP)*bzt79(:,OP_DR)*eta79(:,OP_1)  &
              -ri3_79*pst79(:,OP_DRP)*bzt79(:,OP_DZ)*eta79(:,OP_1))

          if(i3d .eq. 1) then

             o = o +  2.*(gam-1.)*ohfac* &
             (ri3_79*pst79(:,OP_DZP)*bft79(:,OP_DRPP)*eta79(:,OP_1)  &
             -ri3_79*pst79(:,OP_DRP)*bft79(:,OP_DZPP)*eta79(:,OP_1))

             o = o + (gam-1.)*ohfac* &
             (ri2_79*bzt79(:,OP_DZ)*bft79(:,OP_DZPP)*eta79(:,OP_1) &
             +ri2_79*bzt79(:,OP_DR)*bft79(:,OP_DRPP)*eta79(:,OP_1))

             o = o + (gam-1.)*ohfac* &
               (ri2_79*bft79(:,OP_DZPP)*bft79(:,OP_DZPP)*eta79(:,OP_1) &
             +  ri2_79*bft79(:,OP_DRPP)*bft79(:,OP_DRPP)*eta79(:,OP_1))      
          endif
#endif
       end if

! ADD heat source if present
  if(heat_source) o = o + (gam-1.)*q79(:,OP_1)

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
         + 2.*ri_79*(bft79(:,OP_DRP)*pst79(:,OP_DZ) - bft79(:,OP_DZP)*pst79(:,OP_DR))  &
         + bft79(:,OP_DRP)**2 + bft79(:,OP_DZP)**2
#endif

!      v dot b
  temp79b = pht79(:,OP_DR)*pst79(:,OP_DR) +  pht79(:,OP_DZ)*pst79(:,OP_DZ)   &
          + bzt79(:,OP_1)*vzt79(:,OP_1)   &
          + ri3_79*(pst79(:,OP_DR)*cht79(:,OP_DZ) - pst79(:,OP_DZ)*cht79(:,OP_DR))  
 
#if defined(USE3D) || defined(USECOMPLEX)
  temp79b = temp79b     &
          + r_79*(bft79(:,OP_DRP)*pht79(:,OP_DZ) - bft79(:,OP_DZP)*pht79(:,OP_DR))  &
          - ri2_79*( cht79(:,OP_DR)*bft79(:,OP_DRP) +  cht79(:,OP_DZ)*bft79(:,OP_DZP))
#endif
  
  o = temp79b/sqrt(temp79c)

end subroutine vpar_get

subroutine f1vplot_sub(i,term)
  use basic
  use m3dc1_nint
  use metricterms_new
  
  implicit none
  vectype, intent(out) :: term
  integer, intent(in) :: i
  vectype :: temp

     temp = b2psiv(mu79(:,:,i),pst79,vzt79)
     term = temp

     temp = b2bu  (mu79(:,:,i),bzt79,pht79)  &
          + b2bchi(mu79(:,:,i),bzt79,cht79)
     term = term + temp

end subroutine f1vplot_sub

subroutine f1eplot_sub(i,term)
  use basic
  use m3dc1_nint
  use metricterms_new
  
  implicit none
  vectype, intent(out) :: term
  vectype :: temp
  vectype, dimension(MAX_PTS, OP_NUM) :: hf
  integer, intent(in) :: i

  ! Resistive and Hyper Terms
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~
  hf = hypi*sz79
  temp = b2psieta(mu79(:,:,i),pst79,eta79,hf)
  term = temp
  
  temp = b2beta(mu79(:,:,i),bzt79,eta79,hf)
  term = term + temp

  if(i3d.eq.1) then
     temp = b2feta(mu79(:,:,i),bft79,eta79,hf)
     term = term + temp
  end if

end subroutine f1eplot_sub

subroutine f2vplot_sub(i,term)
  use basic
  use m3dc1_nint
  use metricterms_new
  
  implicit none
  vectype, intent(out) :: term
  integer, intent(in) :: i
  vectype :: temp

     temp = b1psiu  (mu79(:,:,i),pst79,pht79) &
          + b1psiv  (mu79(:,:,i),pst79,vzt79) &
          + b1psichi(mu79(:,:,i),pst79,cht79)
     term = temp

     temp =  b1bu  (mu79(:,:,i),bzt79,pht79) &
           + b1bv  (mu79(:,:,i),bzt79,vzt79) &
           + b1bchi(mu79(:,:,i),bzt79,cht79)
     term = term + temp

 if(i3d.eq.1 .and. numvar.ge.2) then
        temp = b1fu  (mu79(:,:,i),bft79,pht79)  &
             + b1fv  (mu79(:,:,i),bft79,vzt79)  &
             + b1fchi(mu79(:,:,i),bft79,cht79)
        term = term + temp
 end if

end subroutine f2vplot_sub

subroutine f2eplot_sub(i,term)
  use basic
  use m3dc1_nint
  use metricterms_new
  
  implicit none
  vectype, intent(out) :: term
  integer, intent(in) :: i
  vectype, dimension(MAX_PTS, OP_NUM) :: hf
  vectype :: temp

      ! Resistive and Hyper Terms
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~
  hf = hypf*sz79
  temp = b1psieta(mu79(:,:,i),pst79,eta79,hf)

  term = temp

end subroutine f2eplot_sub

subroutine f3vplot_sub(i,term)
  use basic
  use m3dc1_nint
  use metricterms_new
  
  implicit none
  vectype, intent(out) :: term
  integer, intent(in) :: i
  vectype :: temp


        temp = t3tnu  (mu79(:,:,i),tet79,nt79,pht79) &
             + t3tnv  (mu79(:,:,i),tet79,nt79,vzt79) &
             + t3tnchi(mu79(:,:,i),tet79,nt79,cht79)
        term = temp
!

end subroutine f3vplot_sub

subroutine f3eplot_sub(i,term)
  use basic
  use m3dc1_nint
  use metricterms_new
  
  implicit none
  vectype, intent(out) :: term
  integer, intent(in) :: i
  vectype :: temp, ohfac
  vectype, dimension(MAX_PTS, OP_NUM) :: hf

  ohfac = 1.
  if(ipres.eq.0) ohfac = 0.5

 ! Ohmic Heating
       temp = b3psipsieta(mu79(:,:,i),pst79,pst79,eta79) &
            + b3bbeta    (mu79(:,:,i),bzt79,bzt79,eta79) &
            + b3psibeta  (mu79(:,:,i),pst79,bzt79,eta79) 
       term = temp*ohfac
          if(i3d .eq. 1) then
             temp = b3psifeta(mu79(:,:,i),pst79,bft79,eta79) &
                  + b3bfeta  (mu79(:,:,i),bzt79,bft79,eta79) &
                  + b3ffeta(mu79(:,:,i),bft79,bft79,eta79)   
             term = temp*ohfac
          endif

! Perpendicular Heat Flux
! ~~~~~~~~~~~~~~~~~~~~~~~
       hf = hypp*sz79
       temp = b3tekappa(mu79(:,:,i),tet79,kap79,hf)
       term = term + temp

! Parallel Heat Flux
! ~~~~~~~~~~~~~~~~~~
  if(kappar.ne.0.) then
        
          temp = tepsipsikappar(mu79(:,:,i),pstx79,pstx79,tet79,b2i79,kar79) &
               + tepsibkappar  (mu79(:,:,i),pstx79,bztx79,tet79,b2i79,kar79) &
               + tebbkappar    (mu79(:,:,i),bztx79,bztx79,tet79,b2i79,kar79)
          term = term + temp
          if(i3d.eq.1 .and. numvar.ge.2) then
             temp = tepsifkappar(mu79(:,:,i),pstx79,bftx79,tet79,b2i79,kar79) &
                  + tebfkappar  (mu79(:,:,i),bztx79,bftx79,tet79,b2i79,kar79) &
                  + teffkappar  (mu79(:,:,i),bftx79,bftx79,tet79,b2i79,kar79)
             term = term + temp
!source terms
             temp = (gam-1.)*b3q(mu79(:,:,i),q79)
             term = term + temp

          endif
  endif
end subroutine f3eplot_sub


end module temperature_plots


 

 
