function read_gamma, filename=filename, _EXTRA=extra

   if(n_elements(filename) eq 0) then filename='C1.h5'

   n = n_elements(filename)
   gamma = fltarr(n)
   pt = 10

   for i=0, n-1 do begin
       ke = read_scalar('ke', file=filename[i], time=t, _EXTRA=extra)
       m = n_elements(ke)
       if(m lt pt) then begin
           gamma[i] = 0
       endif else begin
           g = deriv(t, alog(ke)) / 2.
           gamma[i] = median(g[m-pt:m-1])
       end
   end

   return, gamma
end
