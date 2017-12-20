; ==============================================
; field_at_flux
;
; evaluates field at points where psi=flux
; ==============================================
function field_at_flux, field, psi, x, z, t, flux, theta=theta, angle=angle, $
                        xp=xp, zp=zp, integrate=integrate, axis=axis, $
                        refine=refine, interval=interval, $
                        psilim=psilim, contiguous=contiguous, _EXTRA=extra

   if(n_elements(field) le 1) then return, 0

   xy = path_at_flux(psi,x,z,t,flux,refine=refine, interval=interval, $
                    axis=axis,psilim=psilim,contiguous=contiguous)
   if(n_elements(xy) le 1) then return, 0

   i = n_elements(x)*(xy[0,*]-min(x)) / (max(x)-min(x))
   j = n_elements(z)*(xy[1,*]-min(z)) / (max(z)-min(z))

   test = interpolate(reform(field[0,*,*]),i,j)
   if(n_elements(theta) ne 0) then begin
       angle = interpolate(reform(theta[0,*,*]),i,j)
   endif else begin
       if(n_elements(axis) eq 0) then begin
           angle = findgen(n_elements(test))
       endif else begin
           angle = atan(xy[1,*]-axis[1],xy[0,*]-axis[0])
       endelse
   endelse

   ; re-order array
   mint = min(angle,p)
   p = p+1
   test = shift(test, -p)
   angle = shift(angle, -p)
   xp = shift(xy[0,*],-p)
   zp = shift(xy[1,*],-p)

   if(keyword_set(integrate)) then begin
       r = radius_matrix(x,z,t)
       jac = interpolate(r/a_bracket(theta,psi,x,z),i,j)
       jac = shift(jac, -p)
       dtheta = deriv(angle)
       test = total(test*dtheta*jac, /cumulative)
       d = min(angle, i,/abs)
       test = test - test[i]
   endif

   return, test
end
