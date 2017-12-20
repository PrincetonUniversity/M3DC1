function find_local_extreme, f, x, z, dfdr=dfdr, dfdz=dfdz, guess=guess
  tol = 1e-5*(max(f) - min(f))
  maxits = 20

  if(n_elements(dfdr) eq 0) then dfdr = dx(f,x)
  if(n_elements(dfdz) eq 0) then dfdz = dz(f,z)
  dfdrr = dx(dfdr,x)
  dfdrz = (dx(dfdz,x) + dz(dfdr,z))/2.
  dfdzz = dz(dfdz,z)

  x0 = guess
  dx = [0., 0.]

  for i=0, maxits-1 do begin
     fr = field_at_point(dfdr,x,z,x0[0],x0[1])
     fz = field_at_point(dfdz,x,z,x0[0],x0[1])

     if(i gt 0) then begin
        frr = field_at_point(dfdrr,x,z,x0[0],x0[1])
        frz = field_at_point(dfdrz,x,z,x0[0],x0[1])
        fzz = field_at_point(dfdzz,x,z,x0[0],x0[1])
        
        det = frr*fzz - frz*frz
        
        if(det ne 0.) then begin
           dx[0] = (-fzz*fr + frz*fz) / det
           dx[1] = (-frr*fz + frz*fr) / det
           x0 = x0 + dx 
        endif else begin
           print, 'Error: det == 0'
           return, x0
        endelse
     end

     if((fr*fr + fz*fz) lt tol*tol) then return, x0
  end

  print, 'Error: find_local_extreme did not converge.'
  return, guess
end
