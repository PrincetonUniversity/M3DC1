module harned_mikic_mod

contains


subroutine b1harnedmikic(trial,lin,psiterm,bterm)
  use basic
  use m3dc1_nint

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: trial, lin
  vectype, intent(out) :: psiterm, bterm

  bterm = 0.
  psiterm = 0.

  if(itwofluid.eq.0 .or. surface_int .or. jadv.eq.1 .or. &
       dbf.eq.0. .or. harned_mikic.eq.0.) then
     return
  end if

  temp79a = trial(:,OP_DZ )*pst79(:,OP_DR ) - trial(:,OP_DR )*pst79(:,OP_DZ )
  temp79b = trial(:,OP_DZZ)*pst79(:,OP_DR ) - trial(:,OP_DRZ)*pst79(:,OP_DZ ) &
       +    trial(:,OP_DZ )*pst79(:,OP_DRZ) - trial(:,OP_DR )*pst79(:,OP_DZZ)
  temp79c = trial(:,OP_DRZ)*pst79(:,OP_DR ) - trial(:,OP_DRR)*pst79(:,OP_DZ ) &
       +    trial(:,OP_DZ )*pst79(:,OP_DRR) - trial(:,OP_DR )*pst79(:,OP_DRZ)

  temp79e = ni79(:,OP_1)**2

  temp79d = temp79b*pst79(:,OP_DR) - temp79c*pst79(:,OP_DZ)
  if(itor.eq.1) temp79d = temp79d - ri_79*temp79a*pst79(:,OP_DZ)

  psiterm = -int4(ri2_79,lin(:,OP_GS),temp79d,temp79e)


  temp79d = temp79b*lin(:,OP_DR) - temp79c*lin(:,OP_DZ)
  if(itor.eq.1) temp79d = temp79d - ri_79*temp79a*lin(:,OP_DZ)

  bterm = -int4(ri2_79,bzt79(:,OP_1),temp79d,temp79e)

#ifdef USECOMPLEX
  temp79d = temp79b*lin(:,OP_DZP) + temp79c*lin(:,OP_DRP)
  if(itor.eq.1) temp79d = temp79d + ri_79*temp79a*lin(:,OP_DRP)
  psiterm = psiterm - int4(ri3_79,bzt79(:,OP_1),temp79d,temp79e)

  temp79d = pst79(:,OP_DZ)* &
       (pst79(:,OP_DZZ)*lin(:,OP_DZPP ) + pst79(:,OP_DRZ)*lin(:,OP_DRPP )  &
       +pst79(:,OP_DZ )*lin(:,OP_DZZPP) + pst79(:,OP_DR )*lin(:,OP_DRZPP)) &
       +    pst79(:,OP_DR)* &
       (pst79(:,OP_DRZ)*lin(:,OP_DZPP ) + pst79(:,OP_DRR)*lin(:,OP_DRPP )  &
       +pst79(:,OP_DZ )*lin(:,OP_DRZPP) + pst79(:,OP_DR )*lin(:,OP_DRRPP))
  if(itor.eq.1) then
     temp79d = temp79d - 2.*ri_79*pst79(:,OP_DR) * &
          (pst79(:,OP_DZ)*lin(:,OP_DZPP) + pst79(:,OP_DR)*lin(:,OP_DRPP))
  endif
  psiterm = psiterm + int4(ri4_79,trial(:,OP_1),temp79d,temp79e)

  temp79d = pst79(:,OP_DZ)**2 + pst79(:,OP_DR)**2
  psiterm = psiterm &
       - int5(ri4_79,trial(:,OP_1),temp79d,lin(:,OP_GSPP),temp79e)

  temp79d = rfac* &
       (pst79(:,OP_DZ)*lin(:,OP_DRPP) - pst79(:,OP_DR)*lin(:,OP_DZPP))
  psiterm = psiterm &
       + int5(ri5_79,trial(:,OP_1),bzt79(:,OP_1),temp79d,temp79e)

  temp79d = pst79(:,OP_DZ)* &
       (lin(:,OP_DZZP)*pst79(:,OP_DR ) - lin(:,OP_DRZP)*pst79(:,OP_DZ )  &
       +lin(:,OP_DZP )*pst79(:,OP_DRZ) - lin(:,OP_DRP )*pst79(:,OP_DZZ)) &
       + pst79(:,OP_DR)* &
       (lin(:,OP_DRZP)*pst79(:,OP_DR ) - lin(:,OP_DRRP)*pst79(:,OP_DZ )  &
       +lin(:,OP_DZP )*pst79(:,OP_DRR) - lin(:,OP_DRP )*pst79(:,OP_DRZ))
  if(itor.eq.1) then 
     temp79d = temp79d - ri_79*pst79(:,OP_DR) * &
          (lin(:,OP_DZP)*pst79(:,OP_DR) - lin(:,OP_DRP)*pst79(:,OP_DZ))
  endif
  bterm = bterm - int4(ri3_79,trial(:,OP_1),temp79d,temp79e)

  temp79d = pst79(:,OP_DZ)*lin(:,OP_DZPP) + pst79(:,OP_DR)*lin(:,OP_DRPP)
  bterm = bterm &
       - int5(ri4_79,trial(:,OP_1),bzt79(:,OP_1),temp79d,temp79e)
#endif


  bterm = bterm*dbf**2*harned_mikic
  psiterm = psiterm*dbf**2*harned_mikic
 
end subroutine b1harnedmikic

subroutine b2harnedmikic(trial,lin,psiterm,bterm)

  use basic
  use m3dc1_nint

  vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: trial, lin
  vectype, intent(out) :: psiterm, bterm

  bterm = 0.
  psiterm = 0.

  if(itwofluid.eq.0 .or. surface_int .or. &
       dbf.eq.0. .or. harned_mikic.eq.0.) then
     return
  end if

  temp79e = ni79(:,OP_1)**2

  temp79a = trial(:,OP_DZ )*bzt79(:,OP_DR ) - trial(:,OP_DR )*bzt79(:,OP_DZ )
  temp79b = trial(:,OP_DZZ)*bzt79(:,OP_DR ) - trial(:,OP_DRZ)*bzt79(:,OP_DZ ) &
          + trial(:,OP_DZ )*bzt79(:,OP_DRZ) - trial(:,OP_DR )*bzt79(:,OP_DZZ)
  temp79c = trial(:,OP_DRZ)*bzt79(:,OP_DR ) - trial(:,OP_DRR)*bzt79(:,OP_DZ ) &
          + trial(:,OP_DZ )*bzt79(:,OP_DRR) - trial(:,OP_DR )*bzt79(:,OP_DRZ)

  temp79d = temp79b*pst79(:,OP_DR) - temp79c*pst79(:,OP_DZ)
  if(itor.eq.1) then
     temp79d = temp79d &
          +    ri_79*temp79a*pst79(:,OP_DZ) &
          - 2.*ri_79*trial(:,OP_DZ)* &
          (bzt79(:,OP_DZ)*pst79(:,OP_DR) - bzt79(:,OP_DR)*pst79(:,OP_DZ)) &
          - 2.*ri_79*bzt79(:,OP_1)* &
          (trial(:,OP_DZZ)*pst79(:,OP_DR) - trial(:,OP_DRZ)*pst79(:,OP_DZ)) &
          - 4.*ri2_79*bzt79(:,OP_1)*pst79(:,OP_DZ)*trial(:,OP_DZ)
  endif
  psiterm = psiterm - int4(ri4_79,lin(:,OP_GS),temp79d,temp79e)

  temp79d = temp79b*lin(:,OP_DR) - temp79c*lin(:,OP_DZ)
  if(itor.eq.1) then
     temp79d = temp79d &
          +    ri_79*temp79a*lin(:,OP_DZ) &
          - 2.*ri_79*trial(:,OP_DZ)* &
          (bzt79(:,OP_DZ)*lin(:,OP_DR) - bzt79(:,OP_DR)*lin(:,OP_DZ)) &
          - 2.*ri_79*bzt79(:,OP_1)* &
          (trial(:,OP_DZZ)*lin(:,OP_DR) - trial(:,OP_DRZ)*lin(:,OP_DZ)) &
          - 4.*ri2_79*bzt79(:,OP_1)*lin(:,OP_DZ)*trial(:,OP_DZ)
  endif
  bterm = bterm - int4(ri4_79,bzt79(:,OP_1),temp79d,temp79e)
  

  temp79a = trial(:,OP_DZ )*pst79(:,OP_DR ) - trial(:,OP_DR )*pst79(:,OP_DZ )
  temp79b = trial(:,OP_DZZ)*pst79(:,OP_DR ) - trial(:,OP_DRZ)*pst79(:,OP_DZ ) &
       +    trial(:,OP_DZ )*pst79(:,OP_DRZ) - trial(:,OP_DR )*pst79(:,OP_DZZ)
  temp79c = trial(:,OP_DRZ)*pst79(:,OP_DR ) - trial(:,OP_DRR)*pst79(:,OP_DZ ) &
       +    trial(:,OP_DZ )*pst79(:,OP_DRR) - trial(:,OP_DR )*pst79(:,OP_DRZ)

  if(itor.eq.1) then 
     temp79d = temp79b* &
          (pst79(:,OP_DZZ)*lin(:,OP_DR ) - pst79(:,OP_DRZ)*lin(:,OP_DZ )  &
          +pst79(:,OP_DZ )*lin(:,OP_DRZ) - pst79(:,OP_DR )*lin(:,OP_DZZ)) &
          + (temp79c - ri_79*temp79a)* &
          (pst79(:,OP_DRZ)*lin(:,OP_DR ) - pst79(:,OP_DRR)*lin(:,OP_DZ )  &
          +pst79(:,OP_DZ )*lin(:,OP_DRR) - pst79(:,OP_DR )*lin(:,OP_DRZ)  &
          - ri_79*(pst79(:,OP_DZ)*lin(:,OP_DR) - pst79(:,OP_DR)*lin(:,OP_DZ)))
  else
     temp79d = temp79b* &
          (pst79(:,OP_DZZ)*lin(:,OP_DR ) - pst79(:,OP_DRZ)*lin(:,OP_DZ )  &
          +pst79(:,OP_DZ )*lin(:,OP_DRZ) - pst79(:,OP_DR )*lin(:,OP_DZZ)) &
          +    temp79c* &
          (pst79(:,OP_DRZ)*lin(:,OP_DR ) - pst79(:,OP_DRR)*lin(:,OP_DZ )  &
          +pst79(:,OP_DZ )*lin(:,OP_DRR) - pst79(:,OP_DR )*lin(:,OP_DRZ))
  endif

  bterm = bterm + int3(ri4_79,temp79d,temp79e)


#ifdef USECOMPLEX
  if(itor.eq.1) then 
     temp79d = temp79b* &
          (pst79(:,OP_DZZ)*lin(:,OP_DZP ) + pst79(:,OP_DRZ)*lin(:,OP_DRP )  &
          +pst79(:,OP_DZ )*lin(:,OP_DZZP) + pst79(:,OP_DR )*lin(:,OP_DRZP)) &
          + (temp79c - ri_79*temp79a)* &
          (pst79(:,OP_DRZ)*lin(:,OP_DZP ) + pst79(:,OP_DRR)*lin(:,OP_DRP )  &
          +pst79(:,OP_DZ )*lin(:,OP_DRZP) + pst79(:,OP_DR )*lin(:,OP_DRRP)  &
          - 2.*ri_79* &
          (pst79(:,OP_DZ)*lin(:,OP_DZ) + pst79(:,OP_DR)*lin(:,OP_DRP)))
  else
     temp79d = temp79b* &
          (pst79(:,OP_DZZ)*lin(:,OP_DZP ) + pst79(:,OP_DRZ)*lin(:,OP_DRP )  &
          +pst79(:,OP_DZ )*lin(:,OP_DZZP) + pst79(:,OP_DR )*lin(:,OP_DRZP)) &
          +    temp79c* &
          (pst79(:,OP_DRZ)*lin(:,OP_DZP ) + pst79(:,OP_DRR)*lin(:,OP_DRP )  &
          +pst79(:,OP_DZ )*lin(:,OP_DRZP) + pst79(:,OP_DR )*lin(:,OP_DRRP))
  endif
  psiterm = psiterm + int3(ri5_79,temp79d,temp79e)

  temp79d = pst79(:,OP_DZ)*temp79b + pst79(:,OP_DR)*temp79c
  if(itor.eq.1) temp79d = temp79d - ri_79*temp79a*pst79(:,OP_DR)
  psiterm = psiterm - int4(ri5_79,lin(:,OP_GSP),temp79d,temp79e)

  temp79d = temp79b*lin(:,OP_DRPP) - temp79c*lin(:,OP_DZPP)
  if(itor.eq.1) temp79d = temp79d + ri_79*temp79a*lin(:,OP_DZPP)
  psiterm = psiterm + int4(ri6_79,bzt79(:,OP_1),temp79d,temp79e)

  temp79d = trial(:,OP_DZ)* &
       (pst79(:,OP_DZZ)*lin(:,OP_DZPP ) + pst79(:,OP_DRZ)*lin(:,OP_DRPP )  &
       +pst79(:,OP_DZ )*lin(:,OP_DZZPP) + pst79(:,OP_DR )*lin(:,OP_DRZPP)) &
       + trial(:,OP_DR)* &
       (pst79(:,OP_DRZ)*lin(:,OP_DZPP ) + pst79(:,OP_DRR)*lin(:,OP_DRPP )  &
       +pst79(:,OP_DZ )*lin(:,OP_DRZPP) + pst79(:,OP_DR )*lin(:,OP_DRRPP))
  if(itor.eq.1) then 
     temp79d = temp79d - 2.*ri_79*trial(:,OP_DR)* &
          (pst79(:,OP_DZ)*lin(:,OP_DZPP) + pst79(:,OP_DR)*lin(:,OP_DRPP))
  endif
  psiterm = psiterm - int4(ri6_79,bzt79(:,OP_1),temp79d,temp79e)

  temp79d = trial(:,OP_DZ)*pst79(:,OP_DZ) + trial(:,OP_DR)*pst79(:,OP_DR)
  psiterm = psiterm + int5(ri6_79,bzt79(:,OP_1),temp79d,lin(:,OP_GSPP),temp79e)

  temp79d = rfac* &
       (trial(:,OP_DZ)*lin(:,OP_DRPP) - trial(:,OP_DR)*lin(:,OP_DZPP))
  psiterm = psiterm - int5(ri7_79,bzt79(:,OP_1),bzt79(:,OP_1),temp79d,temp79e)



  temp79d = lin(:,OP_DZP)*temp79b + lin(:,OP_DRP)*temp79c
  if(itor.eq.1) temp79d = temp79d - ri_79*temp79a*lin(:,OP_DRP)
  bterm = bterm - int4(ri5_79,bz179(:,OP_1),temp79d,temp79e)

  temp79d = trial(:,OP_DZ)* &
       (pst79(:,OP_DZZ)*lin(:,OP_DRP ) - pst79(:,OP_DRZ)*lin(:,OP_DZP )  &
       +pst79(:,OP_DZ )*lin(:,OP_DRZP) - pst79(:,OP_DR )*lin(:,OP_DZZP)) &
       +    trial(:,OP_DR)* &
       (pst79(:,OP_DRZ)*lin(:,OP_DRP ) - pst79(:,OP_DRR)*lin(:,OP_DZP )  &
       +pst79(:,OP_DZ )*lin(:,OP_DRRP) - pst79(:,OP_DR )*lin(:,OP_DRZP))
  if(itor.eq.1) then 
     temp79d = temp79d - ri_79*trial(:,OP_DR)* &
          (pst79(:,OP_DZ)*lin(:,OP_DRP) - pst79(:,OP_DR)*lin(:,OP_DZP))
  endif
  bterm = bterm - int4(ri5_79,bzt79(:,OP_1),temp79d,temp79e)

  temp79d = trial(:,OP_DZ)*lin(:,OP_DZPP) + trial(:,OP_DR)*lin(:,OP_DRPP)
  bterm = bterm + int5(ri6_79,bzt79(:,OP_1),bzt79(:,OP_1),temp79d,temp79e)
#endif

  bterm = bterm*dbf**2*harned_mikic
  psiterm = psiterm*dbf**2*harned_mikic

end subroutine b2harnedmikic



end module harned_mikic_mod
