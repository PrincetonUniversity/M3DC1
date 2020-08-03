; Compute the toroidal volume integral of field inside the boundary
; /core: volume inside LCFS instead of inside boundary
function volume_integral, name, x, y, t, slice=slice, filename=filename, points=pts, core=core,_EXTRA=ex

  i3d = read_parameter('3d', filename=filename)
  nplanes = read_parameter('nplanes', filename=filename)
  if(size(name, /type) eq 7) then begin
    if(i3d eq 0) then begin
      field = read_field(name,x,y,t,slices=slice,filename=filename,points=pts,_EXTRA=ex)
    endif else begin
      field = read_field(name,x,y,t,slices=slice,filename=filename,points=pts,taverage=5*nplanes,_EXTRA=ex)
    endelse
  endif else begin
    field = name
  endelse
  R = radius_matrix(x,y,t)
  Z = z_matrix(x,y,t)
  
  if keyword_set(core) then begin
    if(i3d eq 0) then begin
      psin = read_field('psi_norm',x,y,t,slices=slice,filename=filename,points=pts,_EXTRA=ex)
    endif else begin
      psin = read_field('psi_norm',x,y,t,slices=slice,filename=filename,points=pts,taverage=5*nplanes,_EXTRA=ex)
    endelse
    lcfs = find_lcfs(psi,x,y,xpoint=xp,filename=filename,points=pts,_EXTRA=ex)
    mask = psin le 1.0
    ; exclude private flux region
    if (xp[1] lt 0.) then mask = mask and (Z ge xp[1]) else mask = mask and (Z le xp[1])
  endif else begin
    B = get_boundary_path(filename=filename,points=pts)
    P = obj_new('IDLanROI',B[0,*],B[1,*])
    mask = P->ContainsPoints(R,Z)
    mask = reform(mask,[1,pts,pts]) gt 0
    R = reform(R,[1,pts,pts])
    Z = reform(Z,[1,pts,pts])
  endelse
  ; assumes uniform x and uniform y grids
  return, 2.0*!pi*(x[1]-x[0])*(y[1]-y[0])*total(mask*R*field,/nan)
  
end
