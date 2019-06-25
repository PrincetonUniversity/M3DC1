function get_shape, _EXTRA=extra

  shape = fltarr(2)

  psi = read_field('psi',r,z,t,_EXTRA=extra)
  psival = lcfs(psi,r,z,axis=axis,_EXTRA=extra)

  ; plot contour
  path = path_at_flux(psi, r, z, t, psival, /contiguous, axis=axis)

  plot, path[0,*], path[1,*], thick=!p.thick*1.5

  x = axis[0]
  y = axis[1]

  left = 0.
  right = 0.
  top = 0.
  bottom = 0.

  n= n_elements(path[0,*])

  ; find minor radius
  for i=0, n-1 do begin
     x0 = path[0,i]
     y0 = path[1,i]
     if(i lt n-1) then begin
        x1 = path[0,i+1]
        y1 = path[1,i+1]
     endif else begin
        x1 = path[0,0]
        y1 = path[1,0]
     endelse
     dx0 = x-x0
     dx1 = x-x1
     dy0 = y-y0
     dy1 = y-y1
     if(dy0*dy1 lt 0.)  then  begin
        xx = (x0*dy1 - x1*dy0)/(dy1-dy0)
        if(dx0 lt 0) then begin
           right = xx
        endif else begin
           left = xx
        end
     end
  end

  top  = max(path[1,*], i)
  top_x = path[0,i]
  bottom = min(path[1,*], i)
  bottom_x = path[0,i]
  
  a = (right-left)/2.
  r0 = (right+left)/2.
  du = (r0 - top_x)/a
  dl = (r0 - bottom_x)/a
  d = (du + dl)/2.
  b = (top-bottom)/2.

  shape = { a:a, b:b, r0:r0, z0:(bottom+top)/2., kappa:b/a, $
            delta_upper:du, delta_lower:dl, delta:d}

  print, 'r0 = ', shape.r0
  print, 'a = ', shape.a
  print, 'z0 = ', shape.z0
  print, 'b = ', shape.b
  print, 'kappa = ', shape.kappa
  print, 'trangularity = ', shape.delta
  print, 'delta upper = ', shape.delta_upper
  print, 'delta lower = ', shape.delta_lower

  oplot, shape.r0+shape.a*[1., -1., -shape.delta_upper, -shape.delta_lower], $
         shape.z0+[0., 0., shape.b, -shape.b], psym=6

  return, shape
end
