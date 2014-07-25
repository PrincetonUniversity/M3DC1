;======================================================
; is_in_tri
; ~~~~~~~~~
;
; returns 1 if the local coordinates "localp" descibes
; is within the triangle described by a, b, c;
; otherwise, returns 0
;======================================================
function is_in_tri, localp, a, b, c

   small = (a+b+c)*1e-4

   if(localp[1] lt 0. - small) then return, 0
   if(localp[1] gt c + small) then return, 0
   
   x = 1.-localp[1]/c
   if(localp[0] lt -b*x - small) then return, 0
   if(localp[0] gt a*x + small) then return, 0

   return, 1
end
