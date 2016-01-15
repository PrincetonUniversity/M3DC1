;==================================================================
; flux_average_field
; ~~~~~~~~~~~~~~~~~~
;
; Computes the flux averages of a field at a given time.
; This function is only for internal use.  Users should instead
; call flux_average
;==================================================================
function flux_average_field, field, psi, x, z, t, bins=bins, flux=flux, $
                             area=area, volume=volume, dV=dV, psirange=range, $
                             integrate=integrate, r0=r0, surface_weight=sw, $
                             nflux=nflux, elongation=elongation, _EXTRA=extra, $
                             fc=fc

   sz = size(field)

   points = sqrt(sz[2]*sz[3])

   if(n_elements(bins) eq 0) then bins = fix(points)

   print, 'flux averaging with ', bins, ' bins'

   type = size(field, /type)
   if(type eq 6) then begin
       result = complexarr(sz[1], bins)
   endif else begin
       result = fltarr(sz[1], bins)
    endelse

  ; if flux coordinates structure is provided, use it
  if(isa(fc) and not isa(fc, 'Int')) then begin
     print, 'FLUX_AVERAGE_FIELD using provided fc'
     
     f = field_at_point(field, x, z, fc.r, fc.z)
     
     for j=0, fc.n-1 do begin
        result[0, j] = total(f[*,j]*fc.j[*,j])/total(fc.j[*,j])
     end

     if(keyword_set(integrate)) then begin
        val = reform(result[0,*])
        result[0,0] = val[0]*fc.dV[0]*(fc.psi[0] - fc.psi0)/2.
        for j=1, fc.n-1 do begin
           result[0,j] = result[0,j-1] + $
                         (val[j]*fc.dV[j]+val[j-1]*fc.dV[j-1])*(fc.psi[j]-fc.psi[j-1])/2.
        end
     endif

     flux = fc.psi
     nflux = fc.psi_norm
     dV = fc.dV
     area = fc.area
     volume = fc.v
     return, result
  endif


   flux = fltarr(sz[1], bins)
   dV = fltarr(sz[1], bins)
   area = fltarr(sz[1], bins)
   elongation = fltarr(sz[1], bins)

   psival = lcfs(psi, x, z, axis=axis, xpoint=xpoint, flux0=flux0,_EXTRA=extra)
   r0 = axis[0]

   if(n_elements(range) eq 0) then begin
       ; if range not provided, use all flux within lcfs
       range = fltarr(sz[1],2)
       for k=0, sz[1]-1 do range[k,*] = [psival, flux0]
   endif else if(n_elements(range) eq 2) then begin
       oldrange = range
       range = fltarr(sz[1],2)
       for k=0, sz[1]-1 do range[k,*] = oldrange
   endif

   ; remove divertor region from consideration
   div_mask = fltarr(n_elements(x), n_elements(z))
   if(n_elements(xpoint) ge 2) then begin
       if(xpoint[0] ne 0. or xpoint[1] ne 0.) then begin
           if(xpoint[1] lt axis[1]) then begin
               for i=0, n_elements(z)-1 do begin
                   if(z[i] ge xpoint[1]) then div_mask[*,i] = 1.
               end
               print, ' removing points below z=', xpoint[1]
           endif else begin
               for i=0, n_elements(z)-1 do begin
                   if(z[i] le xpoint[1]) then div_mask[*,i] = 1.
               end
               print, ' removing points above z=', xpoint[1]
           endelse
       endif else div_mask = 1.
   endif else begin
       div_mask = 1.
   endelse

   bp = sqrt(s_bracket(psi,psi,x,z))

   for k=0, sz[1]-1 do begin
       dpsi = float(range[k,1] - range[k,0])/float(bins)

       for p=0, bins-1 do begin
           fval = range[k,1] - dpsi*(p+1)
           flux[k,p] = fval + dpsi/2.

           faf = field_at_flux(field[k,*,*], psi[k,*,*], $
                               x, z, t, xp=xp, zp=zp, flux[k,p], /contiguous)
           bpf = field_at_flux(bp[k,*,*], psi[k,*,*], $
                               x, z, t, xp=xp, zp=zp, flux[k,p], /contiguous)

           if(dpsi gt 0.) then bpf = -bpf

           if(n_elements(bpf) lt 3) then continue

           dl = sqrt(deriv(xp)^2 + deriv(zp)^2)

           dV[k,p] = 2.*!pi*total(dl*xp/bpf)
           area[k,p] = 2.*!pi*total(dl*xp)
           if(keyword_set(integrate)) then begin
               result[k,p] = 2.*!pi*total(faf*xp*dl*(-dpsi)/bpf)
           end else begin
               if(keyword_set(sw)) then begin
                   ; weight by differential surface area
                   result[k,p] = total(faf*xp*dl)/total(xp*dl)
               endif else begin
                   ; weight by differential volume (standard definition)
                   result[k,p] = total(faf*xp*dl/bpf)/total(xp*dl/bpf)
               end
           endelse
           width = max(xp) - min(xp)
           height = max(zp) - min(zp)
           
           elongation[k,p] = height / width
       endfor
       if(keyword_set(integrate)) then begin
           val = reform(result[k,*])
           result[k,0] = val[0]/2.
           for i=1, n_elements(val)-1 do begin
               result[k,i] = result[k,i-1] + (val[i]+val[i-1])/2.
           end
       endif
   endfor

   nflux = (flux0-flux)/(flux0-psival)

   volume = fltarr(n_elements(flux))
   for i=1, n_elements(flux)-1 do $
     volume[i] = volume[i-1] + (dV[i]+dV[i-1])/2. * (flux[i] - flux[i-1])

   return, result
end
