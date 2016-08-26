;==============================================================================
; EXTEND_PROFILE
;
; Reads an ascii file with two columns (x and y)
; fits values to y = 0.5*c*[1+d*(1-x)]*{1-tanh[b*(x-a)]} + d
; and uses fit to extend x and y 
; Extended profiles are written to file with '.extended' suffixed
;
; psimax:   largest value of x to extrapolate to (default = 1.05)
; fitrange: domain of x used for fitting (default = [0.95, max(x)])
; minval:   value of d parameter in fit (if not set, d is free parameter)
;==============================================================================
pro tanhfit, x, a, f, pder
   t = tanh(a[1]*(x-a[0]))
   s = 0.5*a[2]*(1.+a[3]*(1.-x))
   f = s*(1.-t)

   if(n_params() ge 4) then begin
       pder = [[-s*(1.-t^2)*(-a[1])],$
               [-s*(1.-t^2)*(x-a[0])],$
               [0.5*(1.+a[3]*(1.-x))*(1.-t)],$
               [0.5*a[2]*(1.-x)*(1.-t)]]
   end
end

pro tanhfit2, x, a, f, pder
   t = tanh(a[1]*(x-a[0]))
   s = 0.5*a[2]*(1.+a[3]*(1.-x))
   f = s*(1.-t) + a[4]

   if(n_params() ge 4) then begin
       pder = [[-s*(1.-t^2)*(-a[1])],$
               [-s*(1.-t^2)*(x-a[0])],$
               [0.5*(1.+a[3]*(1.-x))*(1.-t)],$
               [0.5*a[2]*(1.-x)*(1.-t)],$
               [replicate(1.,n_elements(x))]]
   end
end

pro extend_profile, filein, psimax=psimax, fitrange=fitrange, minval=minval, $
                    smooth=sm, _EXTRA=extra
   if(n_elements(psimax) eq 0) then psimax = 1.05

   r = read_ascii(filein)

   x = reform(r.field1[0,*])
   y = reform(r.field1[1,*])
   
   if(n_elements(fitrange) eq 0) then begin
       fitrange = [0.95, max(x)]
   endif else if(n_elements(fitrange) eq 1) then begin
       fitrange = [fitrange, max(x)]
   end

   print, 'Fitting points in range ', fitrange
   if(psimax le max(x)) then begin
       print, 'Error: profile already extends to ', psimax
   end

   i = where(x ge fitrange[0] and x le fitrange[1], c)
   print, 'Fitting to ', c, ' points'
   if(c eq 0) then begin
       print, 'Error: no data points in range'
       return
   endif

   xf = x[i]
   yf = y[i]

   ; fit profile
   w = 1. / yf
   if(n_elements(minval) eq 0) then begin
       a = [0.98, 1./0.01, max(yf), 0., min(yf)]
       yfit = curvefit(xf, yf, w, a, sig, $
                       function_name='tanhfit2')
   endif else begin
       a = [0.98, 1./0.01, max(yf), 0.]
       yfit = curvefit(xf, yf-minval, w, a, sig, $
                       function_name='tanhfit')
   endelse

   print, 'Fit center: ', a[0]
   print, 'Fit width: ', 1./a[1]
   print, 'Fit height: ', a[2]
   if(n_elements(minval) eq 0) then begin
       print, 'Fit floor: ', a[4]
       if(a[4] le 0.) then begin
           print, 'Warning: fit asymptotes to negative values'
           print, 'Try setting minval = desired asymptotic value'
       endif
   endif
  
   ; extend profile
   print, 'Extending profile to ', psimax
   n = n_elements(x)
   dx = x[n-1] - x[n-2]
   m = (psimax - x[n-1])/dx
   newx = fltarr(n+m)
   newy = fltarr(n+m)
   newx[0:n-1] = x
   newx[n:n+m-1] = x[n-1] + (findgen(m)+1.)*dx
   newy[0:n-1] = y
   if(n_elements(minval) eq 0) then begin
       tanhfit2, newx[n:n+m-1], a, z
       newy[n:n+m-1] = z
   endif else begin
       tanhfit, newx[n:n+m-1], a, z
       newy[n:n+m-1] = z + minval
   endelse

   if(n_elements(sm) ne 0) then begin
      newy = smooth(newy, sm)
   end

   plot, newx, newy, xrange=[fitrange[0], psimax], _EXTRA=extra
   oplot, xf, yf, psym=4

   ; write data
   fileout = filein + '.extended'
   openw, ifile, fileout, /get_lun
   printf, ifile, format='(2F12.6)', transpose([[newx], [newy]])
   close, ifile
end
   
