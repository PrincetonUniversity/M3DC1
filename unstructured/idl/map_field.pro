function map_field, field, ix, iy, mask=mask, outval=val
  sz = size(field)
  n = sz[2]
  result = field
  f = reform(field[0,*,*])
  if(n_elements(val) eq 0) then val=0.

  for i=0, n-1 do begin
     for j=0, n-1 do begin
        if(ix[0,i,j] eq -1) then begin
           result[0,i,j] = val
        endif else begin
           result[0,i,j] = interpolate(f, ix[0,i,j], iy[0,i,j])
        endelse
     end
  end

  mask = ix eq -1

  return, result
end
