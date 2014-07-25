;==================================================================
; flux_coord_field
; ~~~~~~~~~~~~~~~~
;==================================================================
function flux_coord_field, field, psi, x, z, t, slice=slice, area=area, i0=i0,$
                           fbins=fbins,  tbins=tbins, flux=flux, angle=angle, $
                           psirange=frange, nflux=nflux, qval=q, pest=pest, $
                           dV=dV, volume=volume, _EXTRA=extra, qflux=qflux

   if(n_elements(psi) eq 0) then begin
       linear = read_parameter('linear',_EXTRA=extra)
       if(linear eq 1) then begin
           psi = read_field('psi',x,z,t,slice=-1,_EXTRA=extra)
       endif else begin
           psi = read_field('psi',x,z,t,slice=slice,_EXTRA=extra)
       endelse
   endif

   sz = size(field)

   if(n_elements(fbins) eq 0) then fbins = sqrt(sz[2]*sz[3])
   if(n_elements(tbins) eq 0) then tbins = sqrt(sz[2]*sz[3])

   type = size(field, /type)
   if((type eq 6) or (type eq 9)) then begin
       result = complexarr(sz[1], fbins, tbins)
   endif else begin
       result = fltarr(sz[1], fbins, tbins)
   endelse
   flux = fltarr(sz[1], fbins)
   angle = fltarr(sz[1], tbins)
   area = fltarr(sz[1], fbins)
   volume = fltarr(sz[1], fbins)
   dV = fltarr(sz[1], fbins)

   psival = lcfs(psi,x,z,axis=axis,xpoint=xpoint,slice=slice,flux0=flux0, $
                 _EXTRA=extra)
   
   if(n_elements(range) eq 0) then begin
       ; if range not provided, use all flux within lcfs
       range = fltarr(sz[1],2)
       for k=0, sz[1]-1 do range[k,*] = [psival, flux0]
   endif else if(n_elements(range) eq 2) then begin
       oldrange = range
       range = fltarr(sz[1],2)
       for k=0, sz[1]-1 do range[k,*] = oldrange
   endif

   r = radius_matrix(x,z,t)

   bp = sqrt(s_bracket(psi,psi,x,z)/r^2)
   dpsidr = dx(psi,x)
   dpsidz = dz(psi,z)

   if(keyword_set(pest)) then begin
       linear = read_parameter('linear',_EXTRA=extra)
       if(n_elements(i0) le 1) then begin
           print, 'DBG: flux_coord_field reading field'
           if(linear eq 1) then begin
               i0 = read_field('i',x,z,t,slice=-1,_EXTRA=extra)
           endif else begin
               i0 = read_field('i',x,z,t,slice=slice,_EXTRA=extra)
           endelse
       endif
       bt = i0/r
       db = bt/(r*bp)
   endif

   print, 'binning with fbins, tbins=', fbins, tbins
   printed = 0
   left_handed = (flux0-psival)*mean(i0) lt 0 
   print, 'left_handed = ', left_handed

   q = fltarr(sz[1], fbins)

   for k=0, sz[1]-1 do begin

       angle[k,*] = 2.*!pi*findgen(tbins)/float(tbins) - !pi
       dpsi = float(range[k,1] - range[k,0])/float(fbins)

       for p=0, fbins-1 do begin
           flux[k,p] = range[k,1] - dpsi*(p+0.5)
       
           f = field_at_flux(field, psi, x, z, t, flux[k,p], $
                             angle=a, xp=xp, zp=zp, axis=axis, $
                             psilim=psival, /contiguous)

           if(n_elements(xp) le 2) then print, 'Too few points!'

           dx = deriv(xp)
           dz = deriv(zp)
           ds = sqrt(dx^2 + dz^2)
           area[k,p] = 2.*!pi*int_tabulated(findgen(n_elements(ds)),ds*xp)
           if(min(x) eq max(x)) then print, 'X ERROR!'
           if(min(z) eq max(z)) then print, 'Z ERROR!'
           ix = n_elements(x)*(xp - min(x))/(max(x) - min(x))
           iz = n_elements(z)*(zp - min(z))/(max(z) - min(z))
           h = interpolate(reform(bp[k,*,*]),ix,iz)
           dV[k,p] = 2.*!pi*int_tabulated(findgen(n_elements(ds)),ds/h)

           thimp = 0.
           if(keyword_set(pest)) then begin
               g = interpolate(reform(db[k,*,*]),ix,iz)

               ; determine position where angle changes sign
               if(p eq 0) then begin
                   ; if this is the first point, 
                   ; minimize the geometric angle
                   func = a
                   dum = min(func, i, /abs)
               endif else begin
                   ; otherwise, minimize the distance to the
                   ; point normal to the last surface
                   dp = [interpolate(reform(dpsidr[k,*,*]),ix,iz), $
                         interpolate(reform(dpsidz[k,*,*]),ix,iz)] $
                     / (xp*h) * (-dpsi*thimp)
                   func = (xp-dp[0]-p1[0])^2 + (zp-dp[1]-p1[1])^2
                   dum = min(func, i, /abs)
                   func = deriv(func)
               endelse
               da = deriv(func)
               if(da[i] eq 0) then print, 'DA ERROR!'
               index = i
               di = func[i]/da[i]
               if(abs(di) lt 0.5) then begin
                   index = index - di
               endif else begin
                   print, 'index correction too large: ', di, func[i], da[i]
               endelse
                   
               if(index lt 0 or index ge n_elements(xp)) then begin
                   print, 'Interpolation error ', index, p
                   index = i
               end

               p0 = [interpolate(xp,index), interpolate(zp,index)]
               ix0 = n_elements(x)*(p0[0] - min(x))/(max(x) - min(x))
               iz0 = n_elements(z)*(p0[1] - min(z))/(max(z) - min(z))

               dp = [interpolate(reform(dpsidr[k,*,*]),ix0,iz0), $
                     interpolate(reform(dpsidz[k,*,*]),ix0,iz0)] $
                 / (p0[0]*interpolate(reform(bp[k,*,*]),ix0,iz0))
               p1 = p0 + (-dpsi*(1.-thimp))*dp
;               print, 'p0, p1', p0, p1
                             
               ; calculate dt
               dt = ds*g
               if(not left_handed) then dt = -dt

               pest_angle = fltarr(n_elements(dt))
               pest_angle[0] = 0.
               for j=1, n_elements(dt)-1 do begin
                   pest_angle[j] = pest_angle[j-1] + (dt[j]+dt[j-1])/2.
               end
             
               ; center pest_angle to change sign where a changes sign
               pest_angle = pest_angle - interpolate(pest_angle,index)

               ; rescale pest_angle
               q[k,p] = (max(pest_angle)-min(pest_angle)) $
                 /(2.*!pi)
               if(left_handed) then q[k,p] = -q[k,p]
               pest_angle = pest_angle/abs(q[k,p])
;               qval = interpol(q, qflux, flux[k,p])
;               pest_angle = pest_angle/qval
               ; constrain pest angle to +/- pi
               pest_angle = clamp_and_shift(pest_angle, shift=count)
               f = shift(f,-count)

               tol = 0.0001
               while(abs(pest_angle[0]-pest_angle[1]) lt tol) do begin
                   pest_angle = pest_angle[1:n_elements(pest_angle)-1]
                   f = f[1:n_elements(f)-1]
               end
               while(abs(pest_angle[n_elements(pest_angle)-2] - $
                      pest_angle[n_elements(pest_angle)-1]) lt tol) do begin
                   pest_angle = pest_angle[0:n_elements(pest_angle)-2]
                   f = f[0:n_elements(f)-2]
               end

               result[k,p,*] = interpol(f,pest_angle,angle[k,*])

;                problems = where(result[k,p,*] eq 1./0., count)
;                if(count gt 0) then begin
;                    print, "result = ", reform(result[k,p,*])
;                    print, "f = ", reform(f)
;                    print, "pest_angle = ", reform(pest_angle)
;                    stop
;                end

           endif else begin
               result[k,p,*] = interpol(f,a,angle[k,*])
           endelse 
       end

       if(dpsi gt 0) then begin
           dV[k,*] = -dV[k,*]
       end

       for i=1, n_elements(flux)-1 do $
         volume[k,i] = volume[k,i-1] $
         + (dV[k,i]+dV[k,i-1])/2. * (flux[k,i] - flux[k,i-1])
   endfor

   nflux = (flux-flux0)/(psival-flux0)

   return, result
end
