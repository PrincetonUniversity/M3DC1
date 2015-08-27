function flux_coord_field, field, psi, x, z, t, slice=slice, area=area, i0=i0,$
                           fbins=fbins,  tbins=tbins, flux=flux, angle=angle, $
                           psirange=frange, nflux=nflux, qval=q, pest=pest, $
                           dV=dV, volume=volume, _EXTRA=extra, qflux=qflux, $
                           fc=fc

  if(n_elements(field) le 1) then begin
     field = read_field(field,x,z,t,slice=slice,_EXTRA=extra)
  end

  if(isa(fc)) then begin
     print, 'FLUX_COORD_FIELD reusing flux coordinate info'
  endif else begin
     fc = flux_coordinates(slice=slice,_EXTRA=extra, $
                           tbins=tbins, fbins=fbins, $
                           psi0=psi,i0=i0,/fast,x=x,z=z,pest=pest)
  end

  dV = fc.dV
  volume = fc.v
  flux = fc.psi
  nflux = fc.psi_norm
  angle = fc.theta
  q = fc.q
  area = fc.area

  sz = size(field)
  type = size(field, /type)
  if((type eq 6) or (type eq 9)) then begin
     result = complexarr(sz[1], fc.n, fc.m)
  endif else begin
     result = fltarr(sz[1], fc.n, fc.m)
  endelse

  result[0,*,*] = transpose(field_at_point(field, x, z, fc.r, fc.z))
  return, result
end
