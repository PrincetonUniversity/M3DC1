;==============================================================
; eval
; ~~~~
;
; given avec field data "field" for element "elm", the value
; of the field at the local position "localpos"
;==============================================================
function eval, field, localpos, theta, elm, operation=op

   sz = size(field, /dim)
   if(sz[0] eq 80) then begin
       threed = 1
   endif else begin
       threed = 0
   endelse
  
   ;     0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9
   mi = [0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,3,2,1,0]
   ni = [0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,2,3,4,5]

   mi1 = mi - 1 > 0
   mi2 = mi - 2 > 0
   mi3 = mi - 3 > 0
   ni1 = ni - 1 > 0
   ni2 = ni - 2 > 0
   ni3 = ni - 3 > 0

   lp0 = [1., localpos[0], localpos[0]^2, localpos[0]^3, localpos[0]^4, localpos[0]^5]
   lp1 = [1., localpos[1], localpos[1]^2, localpos[1]^3, localpos[1]^4, localpos[1]^5]

   if(n_elements(op) eq 0) then op = 1

   co = cos(theta)
   sn = sin(theta)

   op2 = (op-1) / 10
   op1 = op - op2*10

   sum = 0.
   temp = 0.
   
   case op1 of
      1: temp = lp0[mi]*lp1[ni]
      2: begin
         temp =   co*mi*lp0[mi1]*lp1[ni] $
                - sn*ni*lp0[mi]*lp1[ni1]
      end
      3: begin
         temp =   sn*mi*lp0[mi1]*lp1[ni] $
                + co*ni*lp0[mi]*lp1[ni1]
      end
      4: begin
         temp =   co*co*mi*mi1*lp0[mi2]*lp1[ni] $
                + sn*sn*ni*ni1*lp0[mi]*lp1[ni2] $
                -2.*co*sn*ni*mi*lp0[mi1]*lp1[ni1]
      end
      5: begin
         temp =  co*sn*mi*mi1*lp0[mi2]*lp1[ni] $
               - co*sn*ni*ni1*lp0[mi]*lp1[ni2] $
                 +(co*co-sn*sn) $
                 *ni*mi*lp0[mi1]*lp1[ni1]
      end
      6: begin
         temp =  sn*sn*mi*(mi1)*lp0[mi2]*lp1[ni] $
               + co*co*ni*(ni1)*lp0[mi]*lp1[ni2] $
               + 2.*co*sn*ni*mi*lp0[mi1]*lp1[ni1]
      end
      7: begin
         temp =   mi*mi1*lp0[mi2]*lp1[ni] $
                + ni*ni1*lp1[ni2]*lp0[mi]
      end
      8: begin
         temp = $
            ( co*mi*mi1*mi2*lp0[mi3]*lp1[ni] $
            - sn*ni*ni1*ni2*lp1[ni3]*lp0[mi] $
            - sn*mi*mi1*lp0[mi2]*ni*lp1[ni1] $
            + co*mi*lp0[mi1]*ni*ni1*lp1[ni2])
      end
      9: begin
         temp = $
            ( sn*mi*mi1*mi2*lp0[mi3]*lp1[ni] $
            + co*ni*ni1*ni2*lp1[ni3]*lp0[mi] $
            + co*mi*mi1*lp0[mi2]*ni*lp1[ni1] $
            + sn*mi*lp0[mi1]*ni*ni1*lp1[ni2])
      end
   end

   case op2 of
      0: begin
         sum = field[0:19,elm]*temp
         if(threed eq 1) then begin
            sum = sum + temp* $
                  (field[20:39,elm]*localpos[2]   $
                  +field[40:59,elm]*localpos[2]^2 $
                  +field[60:79,elm]*localpos[2]^3)
         endif
      end
      1: begin
         if(threed eq 1) then begin
            sum = temp* $
                  (field[20:39,elm]   $
                  +field[40:59,elm]*localpos[2]*2. $
                  +field[60:79,elm]*localpos[2]^2*3.)
         endif else sum = 0.
      end
      2: begin
         if(threed eq 1) then begin
            sum = temp* $
                  (field[40:59,elm]*2. $
                  +field[60:79,elm]*localpos[2]*6.)
         endif else sum = 0.
      end
   end

   return, total(sum)
end
