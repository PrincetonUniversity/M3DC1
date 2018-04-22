; Compute the toroidal volume integral of field inside the boundary
; /core: volume inside LCFS instead of inside boundary
function volume_integral, field, slice, x, y, t, filename=filename, points=pts, core=core,_EXTRA=ex

  if(size(field, /type) eq 7) then out = read_field(field,x,y,t,slices=slice,filename=filename,points=pts,_EXTRA=ex)
  R = radius_matrix(x,y,t)
  Z = z_matrix(x,y,t)
  
  if keyword_set(core) then begin
    psin = read_field('psi_norm',x,y,t,slices=slice,filename=filename,points=pts,_EXTRA=ex)
    lcfs = find_lcfs(psi,x,y,xpoint=xp,filename=filename,points=pts,_EXTRA=ex)
    mask = psin le 1.0
    ; exclude private flux region
    if (xp[1] lt 0.) then mask = mask and (Z ge xp[1]) else mask = mask and (Z le xp[1])
  endif else begin
    B = get_boundary_path(filename=filename,points=pts)
    mask = mask_bound(R,Z,B)
  endelse
    
  ; assumes uniform x and uniform y grids
  return, 2.0*!pi*(x[1]-x[0])*(y[1]-y[0])*total(mask*R*out,/nan)
  
end
