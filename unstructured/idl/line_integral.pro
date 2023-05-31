; Compute the line integral of field inside the boundary
function line_integral, name, x, y, t, time=time, slice=slice, filename=filename, points=pts, cutx=cutx, cutz=cutz, _EXTRA=ex

  if(size(name, /type) eq 7) then begin
    field = read_field(name,x,y,t,slices=slice,filename=filename,points=pts,time=time, _EXTRA=ex)
  endif else begin
    field = name
  endelse
  R = radius_matrix(x,y,t)
  Z = z_matrix(x,y,t)
  
  B = get_boundary_path(filename=filename,points=pts)
  P = obj_new('IDLanROI',B[0,*],B[1,*])
  mask = P->ContainsPoints(R, Z)
  mask = reform(mask,[1,pts,pts]) gt 0
  
  field = mask * real_part(field)
  
  if(n_elements(cutx) gt 0) then begin
    data = cumtrapz(y, field_at_point(field, x, y, replicate(cutx, n_elements(x)), y))
  endif else if(n_elements(cutz) gt 0) then begin
    data = cumtrapz(x, field_at_point(field, x, y, x, replicate(cutz, n_elements(y))))
  end
  
  return, data[-1]
  
end
