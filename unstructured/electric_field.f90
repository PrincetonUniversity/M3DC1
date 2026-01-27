module electric_field

contains

subroutine electric_field_R(ilin,o,izone)
  use basic
  use m3dc1_nint

  implicit none

  integer, intent(in) :: ilin, izone
  vectype, dimension(MAX_PTS), intent(out) :: o

  o = 0.
  if(izone.eq.ZONE_VACUUM) return

  if(izone.eq.ZONE_CONDUCTOR) then
     temp79b = etaRZ79(:,OP_1)
  else
     temp79b = eta79(:,OP_1)
  end if
  
  ! eta J
  ! ~~~~~
  if(ilin.eq.1) then
     o = o - temp79b*bz179(:,OP_DZ)*ri_79
#if defined(USE3D) || defined(USECOMPLEX)
     o = o + temp79b* &
          (ri2_79*ps179(:,OP_DRP) - ri_79*bfp179(:,OP_DZP))
#endif
  else
     o = o - temp79b*bzt79(:,OP_DZ)*ri_79
#if defined(USE3D) || defined(USECOMPLEX)
     o = o + temp79b* &
          (ri2_79*pst79(:,OP_DRP) - ri_79*bfpt79(:,OP_DZP))
#endif
  end if
  if(izone.eq.ZONE_CONDUCTOR) return

  ! -VxB
  ! ~~~~
  if(ilin.eq.1) then
     o = o &
          + bz079(:,OP_1)*ph179(:,OP_DR) &
          + bz179(:,OP_1)*ph079(:,OP_DR) &
          - vz079(:,OP_1)*ps179(:,OP_DR) &
          - vz179(:,OP_1)*ps079(:,OP_DR) &
          + ri3_79*bz079(:,OP_1)*ch179(:,OP_DZ) &
          + ri3_79*bz179(:,OP_1)*ch079(:,OP_DZ)
     if(use_external_fields) then
        o = o + bzx79(:,OP_1)*ph079(:,OP_DR) &
             - vz079(:,OP_1)*psx79(:,OP_DR) &
             + ri3_79*bzx79(:,OP_1)*ch079(:,OP_DZ)
     end if
#if defined(USE3D) || defined(USECOMPLEX)
     o = o + vz079(:,OP_1)*r_79*bfp179(:,OP_DZ) &
          + vz179(:,OP_1)*r_79*bfp079(:,OP_DZ)
     if(use_external_fields) then
        o = o + vz079(:,OP_1)*r_79*bfpx79(:,OP_DZ)
     end if
#endif
  else
     o = o &
          + bztx79(:,OP_1)*pht79(:,OP_DR) &
          - vzt79(:,OP_1)*pstx79(:,OP_DR) &
          + ri3_79*bztx79(:,OP_1)*cht79(:,OP_DZ)
#if defined(USE3D) || defined(USECOMPLEX)
     o = o + vzt79(:,OP_1)*r_79*bfptx79(:,OP_DZ)
#endif
#ifdef USEPARTICLES
     o = o &
          - bz079(:,OP_1)*ph079(:,OP_DR) &
          + vz079(:,OP_1)*ps079(:,OP_DR) &
          - ri3_79*bz079(:,OP_1)*ch079(:,OP_DZ)
#if defined(USE3D) || defined(USECOMPLEX)
     o = o - vz079(:,OP_1)*r_79*bfp079(:,OP_DZ)
#endif
#endif
  end if

  !! JxB
  !! ~~~
  !if(db .ne. 0.) then
  !   if(ilin.eq.1) then 
  !      o = o - db*ni79(:,OP_1)*ri2_79* &
  !           (bz079(:,OP_1)*bz179(:,OP_DR) + ps079(:,OP_GS)*ps179(:,OP_DR) &
  !           +bz179(:,OP_1)*bz079(:,OP_DR) + ps179(:,OP_GS)*ps079(:,OP_DR))
  !      if(use_external_fields) then 
  !         o = o - db*ni79(:,OP_1)*ri2_79* &
  !              (ps079(:,OP_GS)*psx79(:,OP_DR) &
  !              +bzx79(:,OP_1) *bz079(:,OP_DR))
  !      end if

!#if defined(USE3D) || defined(USECOMPLEX)
  !      o = o + db*ni79(:,OP_1)* &
  !           (ri_79*ps079(:,OP_GS)*bfp179(:,OP_DZ) &
  !           +ri_79*ps179(:,OP_GS)*bfp079(:,OP_DZ) &
  !           -ri2_79*bz079(:,OP_1)*bfp179(:,OP_DRP) &
  !           -ri2_79*bz179(:,OP_1)*bfp079(:,OP_DRP) &
  !           -ri3_79*bz079(:,OP_1)*ps179(:,OP_DZP) &
  !           -ri3_79*bz179(:,OP_1)*ps079(:,OP_DZP))
  !      if(use_external_fields) then 
  !         o = o + db*ni79(:,OP_1)* &
  !              (ri_79*ps079(:,OP_GS)*bfpx79(:,OP_DZ) &
  !              -ri2_79*bzx79(:,OP_1)*bfp079(:,OP_DRP) &
  !              -ri3_79*bzx79(:,OP_1)*ps079(:,OP_DZP))
  !      end if
!#endif
  !   else
  !      o = o - db*ni79(:,OP_1)*ri2_79* &
  !           (bztx79(:,OP_1)*bzt79(:,OP_DR) + pst79(:,OP_GS)*pstx79(:,OP_DR))
!#if defined(USE3D) || defined(USECOMPLEX)
  !      o = o + db*ni79(:,OP_1)* &
  !           (ri_79*pst79(:,OP_GS)*bfptx79(:,OP_DZ) &
  !           -ri2_79*bztx79(:,OP_1)*bfpt79(:,OP_DRP) &
  !           -ri3_79*bztx79(:,OP_1)*pst79(:,OP_DZP))
!#endif
  !   end if
  !end if

  !! Grad(Pe)
  !! ~~~~~~~~
  !if(db .ne. 0.) then 
  !   if(ilin.eq.1) then
  !      o = o - db*ni79(:,OP_1)*pe179(:,OP_DR)
  !   else
  !      o = o - db*ni79(:,OP_1)*pet79(:,OP_DR)
  !   endif
  !end if

end subroutine electric_field_R

subroutine electric_field_Z(ilin,o,izone)
  use basic
  use m3dc1_nint

  implicit none

  integer, intent(in) :: ilin, izone
  vectype, dimension(MAX_PTS), intent(out) :: o

  o = 0.
  if(izone.eq.ZONE_VACUUM) return

  if(izone.eq.ZONE_CONDUCTOR) then
     temp79b = etaRZ79(:,OP_1)
  else
     temp79b = eta79(:,OP_1)
  end if
  
  ! eta J
  ! ~~~~~
  if(ilin.eq.1) then 
     o = o + temp79b*bz179(:,OP_DR)*ri_79
#if defined(USE3D) || defined(USECOMPLEX)
     o = o + temp79b* &
          (ri2_79*ps179(:,OP_DZP) + ri_79*bfp179(:,OP_DRP))
#endif
  else
     o = o + temp79b*bzt79(:,OP_DR)*ri_79
#if defined(USE3D) || defined(USECOMPLEX)
     o = o + temp79b* &
          (ri2_79*pst79(:,OP_DZP) + ri_79*bfpt79(:,OP_DRP))
#endif
  end if
  if(izone.eq.ZONE_CONDUCTOR) return

  ! -VxB
  ! ~~~~
  if(ilin.eq.1) then 
     o = o &
          + bz079(:,OP_1)*ph179(:,OP_DZ) &
          + bz179(:,OP_1)*ph079(:,OP_DZ) &
          - vz079(:,OP_1)*ps179(:,OP_DZ) &
          - vz179(:,OP_1)*ps079(:,OP_DZ) &
          - ri3_79*bz079(:,OP_1)*ch179(:,OP_DR) &
          - ri3_79*bz179(:,OP_1)*ch079(:,OP_DR)
     if(use_external_fields) then
        o = o + bzx79(:,OP_1)*ph079(:,OP_DZ) &
             - vz079(:,OP_1)*psx79(:,OP_DZ) &
             - ri3_79*bzx79(:,OP_1)*ch079(:,OP_DR)
     end if
#if defined(USE3D) || defined(USECOMPLEX)
     o = o &
          - vz079(:,OP_1)*r_79*bfp179(:,OP_DR) &
          - vz179(:,OP_1)*r_79*bfp079(:,OP_DR)
     if(use_external_fields) then
        o = o - vz079(:,OP_1)*r_79*bfpx79(:,OP_DR)
     end if
#endif
  else
     o = o &
          + bztx79(:,OP_1)*pht79(:,OP_DZ) &
          - vzt79(:,OP_1)*pstx79(:,OP_DZ) &
          - ri3_79*bztx79(:,OP_1)*cht79(:,OP_DR)
#if defined(USE3D) || defined(USECOMPLEX)
     o = o - vzt79(:,OP_1)*r_79*bfptx79(:,OP_DR)
#endif
#ifdef USEPARTICLES
     o = o &
          - bz079(:,OP_1)*ph079(:,OP_DZ) &
          + vz079(:,OP_1)*ps079(:,OP_DZ) &
          + ri3_79*bz079(:,OP_1)*ch079(:,OP_DR)
#if defined(USE3D) || defined(USECOMPLEX)
     o = o + vz079(:,OP_1)*r_79*bfp079(:,OP_DR)
#endif
#endif
  end if


  !! JxB
  !! ~~~
  !if(db .ne. 0.) then
  !   if(ilin.eq.1) then 
  !      o = o - db*ni79(:,OP_1)*ri2_79* &
  !           (bz079(:,OP_1)*bz179(:,OP_DZ) + ps079(:,OP_GS)*ps179(:,OP_DZ) &
  !           +bz179(:,OP_1)*bz079(:,OP_DZ) + ps179(:,OP_GS)*ps079(:,OP_DZ))
  !      if(use_external_fields) then 
  !         o = o - db*ni79(:,OP_1)*ri2_79* &
  !              (ps079(:,OP_GS)*psx79(:,OP_DZ) &
  !              +bzx79(:,OP_1)*bz079(:,OP_DZ))
  !      end if
!#if defined(USE3D) || defined(USECOMPLEX)
  !      o = o - db*ni79(:,OP_1)* &
  !           (ri_79*ps079(:,OP_GS)*bfp179(:,OP_DR) &
  !           +ri_79*ps179(:,OP_GS)*bfp079(:,OP_DR) &
  !           +ri2_79*bz079(:,OP_1)*bfp179(:,OP_DZP) &
  !           +ri2_79*bz179(:,OP_1)*bfp079(:,OP_DZP) &
  !           -ri3_79*bz079(:,OP_1)*ps179(:,OP_DRP) &
  !           -ri3_79*bz179(:,OP_1)*ps079(:,OP_DRP))
  !      if(use_external_fields) then
  !         o = o - db*ni79(:,OP_1)* &
  !              (ri_79*ps079(:,OP_GS)*bfpx79(:,OP_DR) &
  !              +ri2_79*bzx79(:,OP_1)*bfp079(:,OP_DZP) &
  !              -ri3_79*bzx79(:,OP_1)*ps079(:,OP_DRP))
  !      end if
!#endif
  !   else
  !      o = o - db*ni79(:,OP_1)*ri2_79* &
  !           (bztx79(:,OP_1)*bzt79(:,OP_DZ) + pst79(:,OP_GS)*pstx79(:,OP_DZ))
!#if defined(USE3D) || defined(USECOMPLEX)
  !      o = o - db*ni79(:,OP_1)* &
  !           (ri_79*pst79(:,OP_GS)*bfptx79(:,OP_DR) &
  !           +ri2_79*bztx79(:,OP_1)*bfpt79(:,OP_DZP) &
  !           -ri3_79*bztx79(:,OP_1)*pst79(:,OP_DRP))
!#endif
  !   end if
  !end if

  !! Grad(Pe)
  !! ~~~~~~~~
  !if(db .ne. 0.) then 
  !   if(ilin.eq.1) then
  !      o = o - db*ni79(:,OP_1)*pe179(:,OP_DZ)
  !   else
  !      o = o - db*ni79(:,OP_1)*pet79(:,OP_DZ)
  !   end if
  !end if

end subroutine electric_field_Z

subroutine electric_field_phi(ilin,o, izone)
  use basic
  use m3dc1_nint

  implicit none

  integer, intent(in) :: ilin, izone
  vectype, dimension(MAX_PTS), intent(out) :: o

  o = 0.
  if(izone.eq.ZONE_VACUUM) return

  ! eta J
  ! ~~~~~
  if(ilin.eq.1) then
     o = o - ri_79*eta79(:,OP_1)*ps179(:,OP_GS)
  else
     o = o - ri_79*eta79(:,OP_1)*pst79(:,OP_GS)
  end if
  if(izone.eq.ZONE_CONDUCTOR) return

  ! VxB
  ! ~~~
     if(ilin.eq.1) then
        o = o &
             +ps079(:,OP_DZ)*ph179(:,OP_DR)-ps079(:,OP_DR)*ph179(:,OP_DZ) &
             +ps179(:,OP_DZ)*ph079(:,OP_DR)-ps179(:,OP_DR)*ph079(:,OP_DZ) &
             + ri3_79* &
             (ps079(:,OP_DZ)*ch179(:,OP_DZ)+ps079(:,OP_DR)*ch179(:,OP_DR) &
             +ps179(:,OP_DZ)*ch079(:,OP_DZ)+ps179(:,OP_DR)*ch079(:,OP_DR))
        if(use_external_fields) then
           o =  o + &
                psx79(:,OP_DZ)*ph079(:,OP_DR)-psx79(:,OP_DR)*ph079(:,OP_DZ) &
                + ri3_79* &
                (psx79(:,OP_DZ)*ch079(:,OP_DZ)+psx79(:,OP_DR)*ch079(:,OP_DR))
        end if
#if defined(USE3D) || defined(USECOMPLEX) 
        o = o + r_79* &
             (ph079(:,OP_DZ)*bfp179(:,OP_DZ)+ph079(:,OP_DR)*bfp179(:,OP_DR)+ &
              ph179(:,OP_DZ)*bfp079(:,OP_DZ)+ph179(:,OP_DR)*bfp079(:,OP_DR)) &
             + ri2_79* &
             (ch079(:,OP_DZ)*bfp179(:,OP_DR)-ch079(:,OP_DR)*bfp179(:,OP_DZ)+ &
              ch179(:,OP_DZ)*bfp079(:,OP_DR)-ch179(:,OP_DR)*bfp079(:,OP_DZ))
        if(use_external_fields) then
           o = o + r_79* &
                (ph079(:,OP_DZ)*bfpx79(:,OP_DZ) &
                +ph079(:,OP_DR)*bfpx79(:,OP_DR)) &
                + ri2_79* &
                (ch079(:,OP_DZ)*bfpx79(:,OP_DZ) &
                +ch079(:,OP_DR)*bfpx79(:,OP_DR))
        end if
#endif
     else
        o = o &
             +pstx79(:,OP_DZ)*pht79(:,OP_DR)-pstx79(:,OP_DR)*pht79(:,OP_DZ) &
             + ri3_79* &
             (pstx79(:,OP_DZ)*cht79(:,OP_DZ)+pstx79(:,OP_DR)*cht79(:,OP_DR))
#if defined(USE3D) || defined(USECOMPLEX) 
        o = o + r_79* &
             (pht79(:,OP_DZ)*bfptx79(:,OP_DZ) &
             +pht79(:,OP_DR)*bfptx79(:,OP_DR)) &
             + ri2_79* &
             (cht79(:,OP_DZ)*bfptx79(:,OP_DR) &
             -cht79(:,OP_DR)*bfptx79(:,OP_DZ))
#endif
#ifdef USEPARTICLES
        o = o &
             -ps079(:,OP_DZ)*ph079(:,OP_DR)+ps079(:,OP_DR)*ph079(:,OP_DZ) &
             - ri3_79* &
             (ps079(:,OP_DZ)*ch079(:,OP_DZ)+ps079(:,OP_DR)*ch079(:,OP_DR))
#if defined(USE3D) || defined(USECOMPLEX) 
        o = o - r_79* &
             (ph079(:,OP_DZ)*bfp079(:,OP_DZ) &
             +ph079(:,OP_DR)*bfp079(:,OP_DR)) &
             - ri2_79* &
             (ch079(:,OP_DZ)*bfp079(:,OP_DR) &
             -ch079(:,OP_DR)*bfp079(:,OP_DZ))
#endif
#endif
     end if


  !! JxB
  !! ~~~
  !if(db .ne. 0.) then
  !   if(ilin.eq.1) then
  !      o = o + db*ni79(:,OP_1)*ri2_79* &
  !           (bz079(:,OP_DZ)*ps179(:,OP_DR) - bz079(:,OP_DR)*ps179(:,OP_DZ) &
  !           +bz179(:,OP_DZ)*ps079(:,OP_DR) - bz179(:,OP_DR)*ps079(:,OP_DZ))
  !      if(use_external_fields) then
  !         o = o + db*ni79(:,OP_1)*ri2_79* &
  !              (bz079(:,OP_DZ)*psx79(:,OP_DR) - bz079(:,OP_DR)*psx79(:,OP_DZ))
  !      end if
!#if defined(USE3D) || defined(USECOMPLEX)
  !      o = o + db*ni79(:,OP_1)* &
  !           (ri2_79*(bfp079(:,OP_DZP)*ps179(:,OP_DR) &
  !                   +bfp179(:,OP_DZP)*ps079(:,OP_DR) &
  !                   -bfp079(:,OP_DRP)*ps179(:,OP_DZ) &
  !                   -bfp179(:,OP_DRP)*ps079(:,OP_DZ)) &
  !           -ri_79*((bz079(:,OP_DZ)+bfp079(:,OP_DZP))*bfp179(:,OP_DZ) &
  !                  +(bz179(:,OP_DZ)+bfp179(:,OP_DZP))*bfp079(:,OP_DZ) &
  !                  +(bz079(:,OP_DR)+bfp079(:,OP_DRP))*bfp179(:,OP_DR) &
  !                  +(bz179(:,OP_DR)+bfp179(:,OP_DRP))*bfp079(:,OP_DR)) &
  !           -ri3_79*(ps079(:,OP_DZP)*ps179(:,OP_DZ) &
  !                   +ps179(:,OP_DZP)*ps079(:,OP_DZ) &
  !                   +ps079(:,OP_DRP)*ps179(:,OP_DR) &
  !                   +ps179(:,OP_DRP)*ps079(:,OP_DR)) &
  !           -ri2_79*(ps079(:,OP_DZP)*bfp179(:,OP_DR) &
  !                   +ps179(:,OP_DZP)*bfp079(:,OP_DR) &
  !                   -ps079(:,OP_DRP)*bfp179(:,OP_DZ) &
  !                   -ps179(:,OP_DRP)*bfp079(:,OP_DZ)))
  !      if(use_external_fields) then
  !         o = o + db*ni79(:,OP_1)* &
  !              (ri2_79*(bfp079(:,OP_DZP)*psx79(:,OP_DR) &
  !                      -bfp079(:,OP_DRP)*psx79(:,OP_DZ)) &
  !              -ri_79*((bz079(:,OP_DZ)+bfp079(:,OP_DZP))*bfpx79(:,OP_DZ) &
  !                     +(bz079(:,OP_DR)+bfp079(:,OP_DRP))*bfpx79(:,OP_DR)) &
  !              -ri3_79*(ps079(:,OP_DZP)*psx79(:,OP_DZ) &
  !                      +ps079(:,OP_DRP)*psx79(:,OP_DR)) &
  !              -ri2_79*(ps079(:,OP_DZP)*bfpx79(:,OP_DR) &
  !                      -ps079(:,OP_DRP)*bfpx79(:,OP_DZ)))
  !      end if
!#endif
  !   else
  !      o = o + db*ni79(:,OP_1)*ri2_79* &
  !           (bzt79(:,OP_DZ)*pstx79(:,OP_DR) - bzt79(:,OP_DR)*pstx79(:,OP_DZ))
!#if defined(USE3D) || defined(USECOMPLEX)
  !      o = o + db*ni79(:,OP_1)* &
  !           (ri2_79*(bfpt79(:,OP_DZP)*pstx79(:,OP_DR) &
  !                   -bfpt79(:,OP_DRP)*pstx79(:,OP_DZ)) &
  !           -ri_79*((bzt79(:,OP_DZ)+bfpt79(:,OP_DZP))*bfptx79(:,OP_DZ) &
  !                  +(bzt79(:,OP_DR)+bfpt79(:,OP_DRP))*bfptx79(:,OP_DR)) &
  !           -ri3_79*(pst79(:,OP_DZP)*pstx79(:,OP_DZ) &
  !                   +pst79(:,OP_DRP)*pstx79(:,OP_DR)) &
  !           -ri2_79*(pst79(:,OP_DZP)*bfptx79(:,OP_DR) &
  !                   -pst79(:,OP_DRP)*bfptx79(:,OP_DZ)))
!#endif
  !   end if
  !end if

  !! grad(Pe)
  !! ~~~~~~~~
!#if defined(USE3D) || defined(USECOMPLEX)
  !if(db .ne. 0.) then
  !   if(ilin.eq.1) then
  !      o = o - db*ni79(:,OP_1)*ri_79*pe179(:,OP_DP)
  !   else
  !      o = o - db*ni79(:,OP_1)*ri_79*pet79(:,OP_DP)
  !   end if
  !end if
!#endif
end subroutine electric_field_phi

subroutine electric_field_par(ilin,o, izone)
  use basic
  use m3dc1_nint

  implicit none

  integer, intent(in) :: ilin, izone
  vectype, dimension(MAX_PTS), intent(out) :: o
  vectype, dimension(MAX_PTS) :: b2, osign,ere

  o = 0.
  if(izone.eq.ZONE_VACUUM) return

  if(izone.eq.ZONE_CONDUCTOR) then
     temp79b = etaRZ79(:,OP_1)
  else
     temp79b = eta79(:,OP_1)
  end if

  ! eta J
  ! ~~~~~
  o = o - ri2_79*eta79(:,OP_1)*bztx79(:,OP_1)*pst79(:,OP_GS)
  o = o + ri2_79*temp79b*(bzt79(:,OP_DR)*pstx79(:,OP_DR) + bzt79(:,OP_DZ)*pstx79(:,OP_DZ))

#if defined(USE3D) || defined(USECOMPLEX)
  o = o + temp79b* &
       (ri2_79* &
       (bfpt79(:,OP_DRP)*pstx79(:,OP_DR) + bfpt79(:,OP_DZP)*pstx79(:,OP_DZ) &
       -(bfptx79(:,OP_DR)*pst79(:,OP_DRP) + bfptx79(:,OP_DZ)*pst79(:,OP_DZP))) &
       -ri3_79* &
       (pstx79(:,OP_DZ)*pst79(:,OP_DRP) - pstx79(:,OP_DR)*pst79(:,OP_DZP)) &
       +ri_79* &
       ((bzt79(:,OP_DZ)+bfpt79(:,OP_DZP))*bfptx79(:,OP_DR) &
       -(bzt79(:,OP_DR)+bfpt79(:,OP_DRP))*bfptx79(:,OP_DZ)))
#endif

  if(izone.eq.ZONE_PLASMA) then

     !! grad(Pe)
     !! ~~~~~~~~
     !if(db .ne. 0.) then
     !   o = o - db*ni79(:,OP_1)*ri_79* &
     !        (pet79(:,OP_DZ)*pstx79(:,OP_DR) - pet79(:,OP_DR)*pstx79(:,OP_DZ))

!#if defined(USE3D) || defined(USECOMPLEX)
     !   o = o - db*ni79(:,OP_1)* &
     !        (ri2_79*pet79(:,OP_DP)*bztx79(:,OP_1) &
     !        -(pet79(:,OP_DZ)*bfptx79(:,OP_DZ) + pet79(:,OP_DR)*bfptx79(:,OP_DR)))
!#endif
     !end if

  end if

  b2 = (bztx79(:,OP_1)**2 + pstx79(:,OP_DR)**2 + pstx79(:,OP_DZ)**2)*ri2_79
#if defined(USE3D) || defined(USECOMPLEX)
  b2 = b2 &
       + 2.*(pstx79(:,OP_DZ)*bfptx79(:,OP_DR) &
       -     pstx79(:,OP_DR)*bfptx79(:,OP_DZ)) * ri_79 &
       + bfptx79(:,OP_DZ)**2 + bfptx79(:,OP_DR)**2
#endif


  o = o / sqrt(b2)
#if defined(USECOMPLEX)
  return
#else  
  ! eta J_RA
  ! ~~~~~~~~
  if(irunaway .gt. 0) then
      !ere = abs(eta79(:,OP_1)*nre179(:,OP_1) )
      osign = sign(1.0,o)
      !o = osign*(abs(o) - ere)
      !where (sign(1.,o) .ne. osign) 
      !o = 0.
      !end where
      ere = eta79(:,OP_1)*nre179(:,OP_1)
      ! only include RE term when RE current has the same sign as total current
      where (sign(1.,ere) .ne. osign) 
         o = o + ere
      end where

  endif
#endif
 
end subroutine electric_field_par


subroutine electric_field_eta_j(ilin,o)
  use basic
  use m3dc1_nint

  implicit none

  integer, intent(in) :: ilin
  vectype, dimension(MAX_PTS), intent(out) :: o


  ! eta J
  ! ~~~~~
  if(ilin.eq.1) then
     o = - ri_79*eta79(:,OP_1)*ps179(:,OP_GS)
  else
     o = - ri_79*eta79(:,OP_1)*pst79(:,OP_GS)
  end if

end subroutine electric_field_eta_j

subroutine electric_field_psidot(ilin,o)
  use basic
  use m3dc1_nint

  implicit none

  integer, intent(in) :: ilin
  vectype, dimension(MAX_PTS), intent(out) :: o

  ! VxB
  ! ~~~
     if(ilin.eq.1) then
        o =   ps079(:,OP_DZ)*ph179(:,OP_DR)-ps079(:,OP_DR)*ph179(:,OP_DZ) &
             +ps179(:,OP_DZ)*ph079(:,OP_DR)-ps179(:,OP_DR)*ph079(:,OP_DZ) &
             + ri3_79* &
             (ps079(:,OP_DZ)*ch179(:,OP_DZ)+ps079(:,OP_DR)*ch179(:,OP_DR) &
             +ps179(:,OP_DZ)*ch079(:,OP_DZ)+ps179(:,OP_DR)*ch079(:,OP_DR))
        if(use_external_fields) then
           o =  o + &
                psx79(:,OP_DZ)*ph079(:,OP_DR)-psx79(:,OP_DR)*ph079(:,OP_DZ) &
                + ri3_79* &
                (psx79(:,OP_DZ)*ch079(:,OP_DZ)+psx79(:,OP_DR)*ch079(:,OP_DR))
        end if
#if defined(USE3D) || defined(USECOMPLEX) 
        o = o + r_79* &
             (ph079(:,OP_DZ)*bfp179(:,OP_DZ)+ph079(:,OP_DR)*bfp179(:,OP_DR))+&
             (ch179(:,OP_DZ)*bfp079(:,OP_DR)-ch179(:,OP_DR)*bfp079(:,OP_DZ))
        if(use_external_fields) then
           o = o + r_79* &
                (ph079(:,OP_DZ)*bfpx79(:,OP_DZ) &
                +ph079(:,OP_DR)*bfpx79(:,OP_DR))
        end if
#endif
     else
        o =   pstx79(:,OP_DZ)*pht79(:,OP_DR)-pstx79(:,OP_DR)*pht79(:,OP_DZ) &
             + ri3_79* &
             (pstx79(:,OP_DZ)*cht79(:,OP_DZ)+pstx79(:,OP_DR)*cht79(:,OP_DR))
#if defined(USE3D) || defined(USECOMPLEX) 
        o = o + r_79* &
             (pht79(:,OP_DZ)*bfptx79(:,OP_DZ) &
             +pht79(:,OP_DR)*bfptx79(:,OP_DR))+ &
             (cht79(:,OP_DZ)*bfptx79(:,OP_DR) &
             -cht79(:,OP_DR)*bfptx79(:,OP_DZ))
#endif
     end if


  ! eta J
  ! ~~~~~
  if(ilin.eq.1) then
     o = o - ri_79*eta79(:,OP_1)*ps179(:,OP_GS)
  else
     o = o - ri_79*eta79(:,OP_1)*pst79(:,OP_GS)
  end if

  !! JxB
  !! ~~~
  !if(db .ne. 0.) then
  !   if(ilin.eq.1) then
  !      o = o + db*ni79(:,OP_1)*ri2_79* &
  !           (bz079(:,OP_DZ)*ps179(:,OP_DR) - bz079(:,OP_DR)*ps179(:,OP_DZ) &
  !           +bz179(:,OP_DZ)*ps079(:,OP_DR) - bz179(:,OP_DR)*ps079(:,OP_DZ))
!#if defined(USE3D) || defined(USECOMPLEX)
  !      o = o + db*ni79(:,OP_1)* &
  !           (ri2_79*(bfp079(:,OP_DZP)*ps179(:,OP_DR) &
  !                   +bfp179(:,OP_DZP)*ps079(:,OP_DR) &
  !                   -bfp079(:,OP_DRP)*ps179(:,OP_DZ) &
  !                   -bfp179(:,OP_DRP)*ps079(:,OP_DZ)) &
  !           -ri_79*((bz079(:,OP_DZ)+bfp079(:,OP_DZP))*bfp179(:,OP_DZ) &
  !                  +(bz179(:,OP_DZ)+bfp179(:,OP_DZP))*bfp079(:,OP_DZ) &
  !                  +(bz079(:,OP_DR)+bfp079(:,OP_DRP))*bfp179(:,OP_DR) &
  !                  +(bz179(:,OP_DR)+bfp179(:,OP_DRP))*bfp079(:,OP_DR)) &
  !           -ri3_79*(ps079(:,OP_DZP)*ps179(:,OP_DZ) &
  !                   +ps179(:,OP_DZP)*ps079(:,OP_DZ) &
  !                   +ps079(:,OP_DRP)*ps179(:,OP_DR) &
  !                   +ps179(:,OP_DRP)*ps079(:,OP_DR)) &
  !           -ri_79*(ps079(:,OP_DZP)*bfp179(:,OP_DR) &
  !                  +ps179(:,OP_DZP)*bfp079(:,OP_DR) &
  !                  -ps079(:,OP_DRP)*bfp179(:,OP_DZ) &
  !                  -ps179(:,OP_DRP)*bfp079(:,OP_DZ)))
!#endif
  !   else
  !      o = o + db*ni79(:,OP_1)*ri2_79* &
  !           (bzt79(:,OP_DZ)*pst79(:,OP_DR) - bzt79(:,OP_DR)*pst79(:,OP_DZ))
!#if defined(USE3D) || defined(USECOMPLEX)
  !      o = o + db*ni79(:,OP_1)* &
  !           (ri2_79*(bfpt79(:,OP_DZP)*pst79(:,OP_DR) &
  !                   -bfpt79(:,OP_DRP)*pst79(:,OP_DZ)) &
  !           -ri_79*((bzt79(:,OP_DZ)+bfpt79(:,OP_DZP))*bfpt79(:,OP_DZ) &
  !                  +(bzt79(:,OP_DR)+bfpt79(:,OP_DRP))*bfpt79(:,OP_DR)) &
  !           -ri3_79*(pst79(:,OP_DZP)*pst79(:,OP_DZ) &
  !                   +pst79(:,OP_DRP)*pst79(:,OP_DR)) &
  !           -ri_79*(pst79(:,OP_DZP)*bfpt79(:,OP_DR) &
  !                  -pst79(:,OP_DRP)*bfpt79(:,OP_DZ)))
!#endif
  !   end if
  !end if

  !! grad(Pe)
  !! ~~~~~~~~
!#if defined(USE3D) || defined(USECOMPLEX)
  !if(db .ne. 0.) then
  !   if(ilin.eq.1) then
  !      o = o - db*ni79(:,OP_1)*ri_79*pe179(:,OP_DP)
  !   else
  !      o = o - db*ni79(:,OP_1)*ri_79*pet79(:,OP_DP)
  !   end if
  !end if
  ! electric potential
  ! ~~~~~~~~~~~~~~~~~~
  !if(jadv.eq.0 .and. i3d.eq.1) then
     !o = o + ri_79*es179(:,OP_DP)
  !endif
!#endif
end subroutine electric_field_psidot
subroutine electric_field_veldif(ilin,o)
  use basic
  use m3dc1_nint

  implicit none

  integer, intent(in) :: ilin
  vectype, dimension(MAX_PTS), intent(out) :: o

  select case(iveldif)
  case(0)   !default
  ! VxB
  ! ~~~
     if(ilin.eq.1) then
        o =   ps079(:,OP_DZ)*ph179(:,OP_DR)-ps079(:,OP_DR)*ph179(:,OP_DZ) &
             +ps179(:,OP_DZ)*ph079(:,OP_DR)-ps179(:,OP_DR)*ph079(:,OP_DZ) &
             + ri3_79* &
             (ps079(:,OP_DZ)*ch179(:,OP_DZ)+ps079(:,OP_DR)*ch179(:,OP_DR) &
             +ps179(:,OP_DZ)*ch079(:,OP_DZ)+ps179(:,OP_DR)*ch079(:,OP_DR))
        if(use_external_fields) then
           o =  o + &
                psx79(:,OP_DZ)*ph079(:,OP_DR)-psx79(:,OP_DR)*ph079(:,OP_DZ) &
                + ri3_79* &
                (psx79(:,OP_DZ)*ch079(:,OP_DZ)+psx79(:,OP_DR)*ch079(:,OP_DR))
        end if
#if defined(USE3D) || defined(USECOMPLEX) 
        o = o + r_79* &
             (ph079(:,OP_DZ)*bfp179(:,OP_DZ)+ph079(:,OP_DR)*bfp179(:,OP_DR))+&
             (ch179(:,OP_DZ)*bfp079(:,OP_DR)-ch179(:,OP_DR)*bfp079(:,OP_DZ))
        if(use_external_fields) then
           o = o + r_79* &
                (ph079(:,OP_DZ)*bfpx79(:,OP_DZ) &
                +ph079(:,OP_DR)*bfpx79(:,OP_DR))
        end if
#endif
     else
        o =   pstx79(:,OP_DZ)*pht79(:,OP_DR)-pstx79(:,OP_DR)*pht79(:,OP_DZ) &
             + ri3_79* &
             (pstx79(:,OP_DZ)*cht79(:,OP_DZ)+pstx79(:,OP_DR)*cht79(:,OP_DR))
#if defined(USE3D) || defined(USECOMPLEX) 
        o = o + r_79* &
             (pht79(:,OP_DZ)*bfptx79(:,OP_DZ) &
             +pht79(:,OP_DR)*bfptx79(:,OP_DR))+ &
      ri2_79*(cht79(:,OP_DZ)*bfptx79(:,OP_DR) &
             -cht79(:,OP_DR)*bfptx79(:,OP_DZ))
#endif
     end if


  ! electric potential
  ! ~~~~~~~~~~~~~~~~~~
#if defined(USE3D) || defined(USECOMPLEX) 
  if(jadv.eq.0 .and. i3d.eq.1) then
     o = o + ri_79*es179(:,OP_DP)
  endif
#endif
  case(1)   ! pht79 only
   o =   pstx79(:,OP_DZ)*pht79(:,OP_DR)-pstx79(:,OP_DR)*pht79(:,OP_DZ) 
            
#if defined(USE3D) || defined(USECOMPLEX) 
        o = o + r_79* &
             (pht79(:,OP_DZ)*bfptx79(:,OP_DZ) &
             +pht79(:,OP_DR)*bfptx79(:,OP_DR))
#endif
  
  case(2)  ! cht79 onl
   o =   ri3_79*(pstx79(:,OP_DZ)*cht79(:,OP_DZ)+pstx79(:,OP_DR)*cht79(:,OP_DR))
#if defined(USE3D) || defined(USECOMPLEX) 
        o = o + ri2_79*(cht79(:,OP_DZ)*bfptx79(:,OP_DR) &
                       -cht79(:,OP_DR)*bfptx79(:,OP_DZ))
#endif
 
  case(3)  ! elecric potential only
    o = 0.
  ! electric potential
  ! ~~~~~~~~~~~~~~~~~~
#if defined(USE3D) || defined(USECOMPLEX) 
  if(jadv.eq.0 .and. i3d.eq.1) then
     o = ri_79*es179(:,OP_DP)
  endif
#endif
  case(4)
#if defined(USE3D) || defined(USECOMPLEX) 
  o = ri_79*bztx79(:,OP_1)*pht79(:,OP_DP)
  if(jadv.eq.0 .and. i3d.eq.1) then
     o = o + ri_79*es179(:,OP_DP)
  endif
#endif
  case(5)
   o =   pstx79(:,OP_DZ)*pht79(:,OP_DR)-pstx79(:,OP_DR)*pht79(:,OP_DZ) 
            
#if defined(USE3D) || defined(USECOMPLEX) 
   o = o + r_79* &
             (pht79(:,OP_DZ)*bfptx79(:,OP_DZ) &
             +pht79(:,OP_DR)*bfptx79(:,OP_DR))
   o = o - ri_79*bztx79(:,OP_1)*pht79(:,OP_DP)
#endif

  case(6)
  o = 0
#if defined(USE3D) || defined(USECOMPLEX) 
        o = o + ri2_79*(cht79(:,OP_DZ)*bfptx79(:,OP_DR) &
                       -cht79(:,OP_DR)*bfptx79(:,OP_DZ))
#endif
  end select
end subroutine electric_field_veldif
subroutine ef_eta_jdb(ilin,o)
  use basic
  use m3dc1_nint

  implicit none

  integer, intent(in) :: ilin
  vectype, dimension(MAX_PTS), intent(out) :: o

  o = ri2_79*bzt79(:,OP_DR)*pst79(:,OP_DR)   &
    + ri2_79*bzt79(:,OP_DZ)*pst79(:,OP_DZ)   &
    - ri2_79*bzt79(:,OP_1)*pst79(:,OP_GS)

  ! ~~~~~~~~~~~~~~~~~~
#if defined(USE3D) || defined(USECOMPLEX) 
  if(jadv.eq.0 .and. i3d.eq.1) then

  o = o + ri2_79*bfpt79(:,OP_DRP)*pst79(:,OP_DR)  &
        + ri2_79*bfpt79(:,OP_DZP)*pst79(:,OP_DZ)  &
        - ri_79 *bfpt79(:,OP_DZ)*bzt79(:,OP_DR)   &
        + ri_79 *bfpt79(:,OP_DR)*bzt79(:,OP_DZ)   &
        - ri_79 *bfpt79(:,OP_DZ)*bfpt79(:,OP_DRP) &
        + ri_79 *bfpt79(:,OP_DR)*bfpt79(:,OP_DZP) &
        + ri3_79*pst79(:,OP_DZP)*pst79(:,OP_DR)   &
        - ri3_79*pst79(:,OP_DRP)*pst79(:,OP_DZ)   &
        - ri2_79*pst79(:,OP_DRP)*bfpt79(:,OP_DR)  &
        - ri2_79*pst79(:,OP_DZP)*bfpt79(:,OP_DZ) 

  endif
#endif
     o = o*eta79(:,OP_1)
end subroutine ef_eta_jdb
subroutine jdbobs_sub(o)
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


  o = ri2_79*bzt79(:,OP_DR)*pst79(:,OP_DR)   &
    + ri2_79*bzt79(:,OP_DZ)*pst79(:,OP_DZ)   &
    - ri2_79*bzt79(:,OP_1)*pst79(:,OP_GS)

  ! ~~~~~~~~~~~~~~~~~~
#if defined(USE3D) || defined(USECOMPLEX) 
  if(jadv.eq.0 .and. i3d.eq.1) then

  o = o + ri2_79*bfpt79(:,OP_DRP)*pst79(:,OP_DR)  &
        + ri2_79*bfpt79(:,OP_DZP)*pst79(:,OP_DZ)  &
        - ri_79 *bfpt79(:,OP_DZ)*bzt79(:,OP_DR)   &
        + ri_79 *bfpt79(:,OP_DR)*bzt79(:,OP_DZ)   &
        - ri_79 *bfpt79(:,OP_DZ)*bfpt79(:,OP_DRP) &
        + ri_79 *bfpt79(:,OP_DR)*bfpt79(:,OP_DZP) &
        + ri3_79*pst79(:,OP_DZP)*pst79(:,OP_DR)   &
        - ri3_79*pst79(:,OP_DRP)*pst79(:,OP_DZ)   &
        - ri2_79*pst79(:,OP_DRP)*bfpt79(:,OP_DR)  &
        - ri2_79*pst79(:,OP_DZP)*bfpt79(:,OP_DZ) 

  endif
#endif
     o = o/temp79c
end subroutine jdbobs_sub
subroutine ef_bdgp(ilin,o)
  use basic
  use m3dc1_nint

  implicit none

  integer, intent(in) :: ilin
  vectype, dimension(MAX_PTS), intent(out) :: o

  o = 0.
  select case(ibdgp)
  case(0)   !  full evaluation

  o = - ri_79*es179(:,OP_DZ)*pst79(:,OP_DR)  &
   +  ri_79*es179(:,OP_DR)*pst79(:,OP_DZ)

#if defined(USE3D) || defined(USECOMPLEX) 
  if(jadv.eq.0 .and. i3d.eq.1) then
     o = o + es179(:,OP_DR)*bfpt79(:,OP_DR)  &
           + es179(:,OP_DZ)*bfpt79(:,OP_DZ)
     o = o - ri2_79*bzt79(:,OP_1)*es179(:,OP_DP)
  endif
#endif

  case(1)  ! only psi term

  o = - ri_79*es179(:,OP_DZ)*pst79(:,OP_DR)  &
   +  ri_79*es179(:,OP_DR)*pst79(:,OP_DZ)

  case(2)  ! only f term

#if defined(USE3D) || defined(USECOMPLEX) 
  if(jadv.eq.0 .and. i3d.eq.1) then
     o = o + es179(:,OP_DR)*bfpt79(:,OP_DR)  &
           + es179(:,OP_DZ)*bfpt79(:,OP_DZ)
  endif
#endif
  case(3)  ! only F (bz) term

#if defined(USE3D) || defined(USECOMPLEX) 
  if(jadv.eq.0 .and. i3d.eq.1) then
     o = o - ri2_79*bzt79(:,OP_1)*es179(:,OP_DP)
  endif
#endif









  end select
end subroutine ef_bdgp
subroutine ef_vlbdgp(ilin,o)
  use basic
  use m3dc1_nint
  use math

  implicit none

  integer, intent(in) :: ilin
  vectype, dimension(MAX_PTS), intent(out) :: o


 o =  -vloop*ri2_79*bzt79(:,OP_1)/toroidal_period


end subroutine ef_vlbdgp
end module electric_field
