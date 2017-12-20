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
   
   mi = [0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,3,2,1,0]
   ni = [0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,2,3,4,5]
   sum = 0.

   if(n_elements(op) eq 0) then op = 1

   co = cos(theta)
   sn = sin(theta)

   op2 = (op-1) / 10
   op1 = op - op2*10

   for p=0, 19 do begin
       temp = 0.
        case op1 of
        1: temp = localpos[0]^mi[p]*localpos[1]^ni[p]
        2: begin
            if(mi[p] ge 1) then $
              temp = temp $
              + co*mi[p]*localpos[0]^(mi[p]-1)*localpos[1]^ni[p]
            if(ni[p] ge 1) then $
              temp = temp $
              - sn*ni[p]*localpos[0]^mi[p]*localpos[1]^(ni[p]-1)
           end
        3: begin
            if(mi[p] ge 1) then $
              temp = temp $
              + sn*mi[p]*localpos[0]^(mi[p]-1)*localpos[1]^ni[p]
            if(ni[p] ge 1) then $
              temp = temp $
              + co*ni[p]*localpos[0]^mi[p]*localpos[1]^(ni[p]-1)
           end
        4: begin
            if(mi[p] ge 2) then $
              temp = temp $
              + co*co*mi[p]*(mi[p]-1)*localpos[0]^(mi[p]-2)*localpos[1]^ni[p]
            if(ni[p] ge 2) then $
              temp = temp $
              + sn*sn*ni[p]*(ni[p]-1)*localpos[0]^mi[p]*localpos[1]^(ni[p]-2)
            if(mi[p] ge 1 and ni[p] ge 1) then $
              temp = temp $
              -2.*co*sn*ni[p]*mi[p]*localpos[0]^(mi[p]-1)*localpos[1]^(ni[p]-1)
        end
        5: begin
            if(mi[p] ge 2) then $
              temp = temp $
              + co*sn*mi[p]*(mi[p]-1)*localpos[0]^(mi[p]-2)*localpos[1]^ni[p]
            if(ni[p] ge 2) then $
              temp = temp $
              - co*sn*ni[p]*(ni[p]-1)*localpos[0]^mi[p]*localpos[1]^(ni[p]-2)
            if(mi[p] ge 1 and ni[p] ge 1) then $
              temp = temp $
              +(co*co-sn*sn) $
              *ni[p]*mi[p]*localpos[0]^(mi[p]-1)*localpos[1]^(ni[p]-1)
           end
        6: begin
            if(mi[p] ge 2) then $
              temp = temp $
              + sn*sn*mi[p]*(mi[p]-1)*localpos[0]^(mi[p]-2)*localpos[1]^ni[p]
            if(ni[p] ge 2) then $
              temp = temp $
              + co*co*ni[p]*(ni[p]-1)*localpos[0]^mi[p]*localpos[1]^(ni[p]-2)
            if(mi[p] ge 1 and ni[p] ge 1) then $
              temp = temp $
              +2.*co*sn*ni[p]*mi[p]*localpos[0]^(mi[p]-1)*localpos[1]^(ni[p]-1)
           end
        7: begin
            if(mi[p] ge 2) then $
              temp = temp + $
              mi[p]*(mi[p]-1)*localpos[0]^(mi[p]-2)*localpos[1]^ni[p]
            if(ni[p] ge 2) then $
              temp = temp + $
              ni[p]*(ni[p]-1)*localpos[1]^(ni[p]-2)*localpos[0]^mi[p]
           end
        8: begin
            temp = temp + $
             ( co $
              *mi[p]*(mi[p]-1)*(mi[p]-2)*localpos[0]^(mi[p]-3>0) $
              *                          localpos[1]^ ni[p]      $
             - sn $
              *                          localpos[0]^ mi[p]      $
              *ni[p]*(ni[p]-1)*(ni[p]-2)*localpos[1]^(ni[p]-3>0) $
             - sn $
              *mi[p]*(mi[p]-1)*          localpos[0]^(mi[p]-2>0) $
              *ni[p]*                    localpos[1]^(ni[p]-1>0) $
             + co $
              *mi[p]*                    localpos[0]^(mi[p]-1>0) $
              *ni[p]*(ni[p]-1)*          localpos[1]^(ni[p]-2>0))
           end
        9: begin
            temp = temp + $
             ( sn $
              *mi[p]*(mi[p]-1)*(mi[p]-2)*localpos[0]^(mi[p]-3>0) $
              *                          localpos[1]^ ni[p]      $
             + co $
              *                          localpos[0]^ mi[p]      $
              *ni[p]*(ni[p]-1)*(ni[p]-2)*localpos[1]^(ni[p]-3>0) $
             + co $
              *mi[p]*(mi[p]-1)*          localpos[0]^(mi[p]-2>0) $
              *ni[p]*                    localpos[1]^(ni[p]-1>0) $
             + sn $
              *mi[p]*                    localpos[0]^(mi[p]-1>0) $
              *ni[p]*(ni[p]-1)*          localpos[1]^(ni[p]-2>0))
        end
       end


       case op2 of
           0: begin
               sum = sum + field[p,elm]*temp
               if(threed eq 1) then begin
                   sum = sum + temp* $
                     (field[p+20,elm]*localpos[2]   $
                     +field[p+40,elm]*localpos[2]^2 $
                     +field[p+60,elm]*localpos[2]^3)
               endif
           end
           1: begin
               if(threed eq 1) then begin
                   sum = sum + temp* $
                     (field[p+20,elm]   $
                     +field[p+40,elm]*localpos[2]*2. $
                     +field[p+60,elm]*localpos[2]^2*3.)
               endif
           end
           2: begin
               if(threed eq 1) then begin
                   sum = sum + temp* $
                     (field[p+40,elm]*2. $
                     +field[p+60,elm]*localpos[2]*6.)
               endif
           end
       end
   end

   return, sum
end
