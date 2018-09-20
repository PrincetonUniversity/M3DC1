function path_at_flux, psi,x,z,t,flux,refine=refine,$
                       interval=interval, axis=axis, psilim=psilim, $
                       contiguous=contiguous, path_points=pts

   contour, psi[0,*,*], x, z, levels=flux, closed=0, $
     path_xy=xy, path_info=info, /path_data_coords, /overplot

   if(n_elements(xy) eq 0) then begin
       print, 'Error: no points at this flux value', flux
       return, 0
   end

   ; refine the surface found by contour with a single newton iteration
   if(keyword_set(refine)) then begin
       i = n_elements(x)*(xy[0,*]-min(x)) / (max(x)-min(x))
       j = n_elements(z)*(xy[1,*]-min(z)) / (max(z)-min(z))
       f = interpolate(reform(psi[0,*,*]),i,j)
       fx = interpolate(reform(dx(psi[0,*,*],x)),i,j)
       fz = interpolate(reform(dz(psi[0,*,*],z)),i,j)
       gf2 = fx^2 + fz^2
       l = (flux-f)/gf2
       xy[0,*] = xy[0,*] + l*fx
       xy[1,*] = xy[1,*] + l*fz
   endif

   ; find the biggest closed path
   if(n_elements(info) gt 1 and keyword_set(contiguous) $
      and n_elements(axis) ne 0) then begin
      for k=0, n_elements(info)-1 do begin
         path = xy[*,info[k].offset:info[k].offset+info[k].n-1]
         if(path_contains_point(path, axis)) then begin
            xy = path
            break
         end
      end
;       ibig = 0
;       nbig = 0
;       for k=0, n_elements(info)-1 do begin
;           if(info[k].n gt nbig) then begin
;               ibig = k
;               nbig = info[k].n
;           end
;       end
;       xy = xy[*,info[ibig].offset:info[ibig].offset+info[ibig].n-1]
   end

   if(n_elements(interval) ne 0) then begin
       if(interval eq 0) then begin
           spline_p, xy[0,*], xy[1,*], xp_new, zp_new
       endif else begin
           spline_p, xy[0,*], xy[1,*], xp_new, zp_new, interval=interval
       endelse
       xy = transpose([[xp_new],[zp_new]])
   endif

   if(n_elements(pts) ne 0) then begin
       ind = findgen(pts)*n_elements(xy[0,*])/pts
       oldxy = xy
       xy = fltarr(2,pts)
       xy[0,*] = interpolate(oldxy[0,*], ind)
       xy[1,*] = interpolate(oldxy[1,*], ind)
   end

   return, xy
end
