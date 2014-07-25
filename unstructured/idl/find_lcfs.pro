; ========================================================
; find_lcfs
; ~~~~~~~~~
;
; returns the flux value of the last closed flux surface
; ========================================================
function find_lcfs, psi, x, z, axis=axis, xpoint=xpoint, $
               flux0=flux0, psilim=psilim, xlim=xlim, _EXTRA=extra

   if(n_elements(psi) eq 0 or n_elements(x) eq 0 or n_elements(z) eq 0) then $
     psi = read_field('psi',x,z,t,_EXTRA=extra,linear=0)

   nulls, psi, x, z, xpoint=xpoint, axis=axis, _EXTRA=extra

   ; flux at magnetic axis
   if(n_elements(axis) lt 2) then begin
       print, "Error: no magnetic axis"
       flux0 = max(psi[0,*,*])
   endif else begin
       flux0 = field_at_point(psi,x,z,axis[0],axis[1])
   endelse
   if(n_elements(axis) gt 2) then begin
       print, "Warning: there is more than one magnetic axis"
   endif
   print, "Flux on axis:", flux0

   ; If xlim is set, then use limiter values
   ; given in C1input file
   if(keyword_set(xlim)) then begin
       limiter = fltarr(2)
       limiter[0] = read_parameter('xlim', _EXTRA=extra)
       limiter[1] = read_parameter('zlim', _EXTRA=extra)
       print, "limiter at: ", limiter
   endif

   ; limiting value
   ; Find limiting flux by calculating outward normal derivative of
   ; the normalized flux.  If this derivative is negative, there is a
   ; limiter.
   if(n_elements(psilim) eq 0) then begin
       if(n_elements(limiter) eq 2) then begin
           xerr = min(x-limiter[0],xi,/absolute)
           zerr = min(z-limiter[1],zi,/absolute)
           psilim = psi[0,xi,zi]
       endif else begin
           print, ' No limiter provided, using wall.'
           sz = size(psi)
           
           psiz = dz(psi,z)
           psix = dx(psi,x)
           
           normal_mask = psi*0.
           normal_mask[0,      *,      0] = 1.
           normal_mask[0,      *,sz[3]-1] = 1.
           normal_mask[0,      0,      *] = 1.
           normal_mask[0,sz[2]-1,      *] = 1.
           
           xx = fltarr(1,sz[2],sz[3])
           zz = fltarr(1,sz[2],sz[3])
           for i=0,sz[3]-1 do xx[0,*,i] = x - axis[0]
           for i=0,sz[2]-1 do zz[0,i,*] = z - axis[1]
           normal_deriv = (psix*xx + psiz*zz)*normal_mask
           
           normal_deriv = normal_deriv lt 0
           
           psi_bound = psi*normal_deriv + (1-normal_deriv)*1e10
           
           psilim = min(psi_bound-flux0, i, /absolute)
           psilim = psi_bound[i]
       endelse
       
       print, "Flux at limiter", psilim
       
       
       ; flux at separatrix
       sz = size(xpoint)
       if(sz[0] gt 0 and (not keyword_set(xlim))) then begin
           psix = field_at_point(psi,x,z,xpoint[0],xpoint[1])
           print, "Flux at separatrix:", psix
           
           if(abs(psix-flux0) gt abs(psilim-flux0)) then begin
               print, "Plasma is limited."
           endif else begin
               print, "Plasma is diverted."
               psilim = psix
           endelse
       endif
   end
   
   return, psilim
end
