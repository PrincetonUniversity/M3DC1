;==================================================================
; flux_average_field
; ~~~~~~~~~~~~~~~~~~
;
; Computes the flux averages of a field at a given time.
; This function is only for internal use.  Users should instead
; call flux_average
;==================================================================
function flux_average_field, field, psi, x, z, t, bins=bins, flux=flux, $
                             area=area, volume=volume, psirange=range, $
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
   endif else begin
      fc = flux_coordinates(_EXTRA=extra)
   endelse
   
   f = reform(field_at_point(field, x, z, fc.r, fc.z))
     
   for j=0, fc.n-1 do begin
      result[0, j] = total(f[*,j]*fc.j[*,j])/total(fc.j[*,j])
   end

   if(keyword_set(integrate)) then begin
      val = reform(result[0,*])
      result[0,0] = val[0]*fc.dV[0]*fc.psi_norm[0]/2.
      for j=1, fc.n-1 do begin
         result[0,j] = result[0,j-1] + $
                       (val[j]*fc.dV_dchi[j]+val[j-1]*fc.dV_dchi[j-1]) $
                       *(fc.psi_norm[j]-fc.psi_norm[j-1])/2.
      end
   endif

   flux = fc.psi
   nflux = fc.psi_norm
   dV = fc.dV_dchi
   area = fc.area
   volume = fc.v
   return, result
end
