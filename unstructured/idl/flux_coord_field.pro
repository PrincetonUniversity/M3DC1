function flux_coord_field, field, psi, x, z, t, slice=slice, area=area, i0=i0,$
                           fbins=fbins,  tbins=tbins, flux=flux, angle=angle, $
                           psirange=frange, nflux=nflux, qval=q, pest=pest, $
                           volume=volume, qflux=qflux, $
                           fc=fc, dpsi0_dx=psi0_r, dpsi0_dz=psi0_z, $
                           _EXTRA=extra

  if(n_elements(field) le 1) then begin
     field = read_field(field,x,z,t,slice=slice,_EXTRA=extra)
  end

  if(isa(fc)) then begin
     print, 'FLUX_COORD_FIELD reusing flux coordinate info'
  endif else begin
     print, 'FLUX_COORD_FIELD NOT reusing flux coordinate info'
     fc = flux_coordinates(slice=slice, $
                           tbins=tbins, fbins=fbins, $
                           dpsi0_dx=psi0_r, dpsi0_dz=psi0_z, $
                           psi0=psi,i0=i0,/fast,x=x,z=z,pest=pest, $
                           _EXTRA=extra)
  end

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

  for i=0, sz[1]-1 do begin
     result[i,*,*] = transpose(field_at_point(field[i,*,*], x, z, fc.r, fc.z))
  end
  return, result
end
