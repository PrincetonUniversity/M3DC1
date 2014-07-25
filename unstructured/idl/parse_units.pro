;=====================================================================
; parse_units
; ~~~~~~~~~~~
;
; x is a vector containing the dimensions of
; [4pi, c, n0, vA0, B0, T0, i0, e0, tA0, L0, ePot]
; [  0, 1,  2,   3,  4,  5,  6,  7,   8,  9,   10]
; 
; output is a string containing units
;=====================================================================
function parse_units, x, cgs=cgs, mks=mks
   if(keyword_set(cgs)) then begin
       x[9] = x[9] - 3*x[2] + x[3] + x[1]
       x[8] = x[8]          - x[3] - x[1]
       x[0] = 0
       x[1] = 0
       x[2] = 0
       x[3] = 0
       u = ['!64!7p', '!8c', '!6cm', '!8v!DA!60!N', $
            '!6G', '!6eV', '!6statamps', '!6erg', '!6s', '!6cm', $
            '!6statvolts'] + '!X'
   endif else if(keyword_set(mks)) then begin
       x[9] = x[9] - 3*x[2] + x[3] + x[1]
       x[8] = x[8]          - x[3] - x[1]
       x[0] = 0
       x[1] = 0
       x[2] = 0
       x[3] = 0
       u = ['!64!7p', '!8c', '!6cm', '!8v!DA!60!N', $
            '!6T', '!6eV', '!6A', '!6J', '!6s', '!6m', $
            '!6V'] + '!X'
   endif else begin
       x[0] = x[0]   - x[5] - x[6] -   x[7]
       x[1] = x[1]          + x[6]
       x[4] = x[4] + 2*x[5] + x[6] + 2*x[7]
       x[9] = x[9]          + x[6] - 3*x[7]
       x[5] = 0
       x[6] = 0
       x[7] = 0

       u = ['!64!7p', '!8c', '!8n!D!60!N', '!8v!DA!60!N', $
            '!8B!D!60!N', '!6temp', '!6curr', $
            '!6energy', '!7s!DA!60!N', '!8L!D!60!N', $
            '!6potential'] + '!X'
   endelse
   units = ''

   nu = n_elements(x)

   if(max(x) gt 0) then pos=1 else pos=0
   if(min(x) lt 0) then neg=1 else neg=0

   is = 0
   sscript = '("!U!6",G0,"!N!X")'
   for i=0, nu-1 do begin
       if(x[i] gt 0) then begin 
           if(is eq 1) then units = units + ' '
           units = units + u[i]
           if(x[i] ne 1) then $
             units = units + string(format=sscript,x[i])
           is = 1
       endif
   end

   if(pos eq 1) then begin
       if(neg eq 1) then begin
           units = units + '!6/!X'
           is = 0
       end
       for i=0, nu-1 do begin
           if(x[i] lt  0) then begin
               if(is eq 1) then units = units + ' '
               units = units + u[i]
               if(x[i] ne -1) then $
                   units = units + string(format=sscript,-x[i])
               is = 1
           endif
       end
   endif else begin
       for i=0, nu-1 do begin
           if(x[i] lt 0) then begin
               if(is eq 1) then units = units + ' '
               units = units + u[i]
               units = units + string(format=sscript, x[i])
               is = 1
           endif
       end
   endelse

   return, units
end
