; find indices i0, j0 such that
; rfield[i0,j0] = r[i,j]
; zfield[i0,jz] = z[i,j]
; on a regular i,j grid
pro create_map, rfield, zfield, r=r, z=z, ix=i0, iy=j0
  i0 = rfield
  j0 = zfield

  sz = size(rfield)
  n = sz[2]

  i = findgen(n)
  rf = reform(rfield[0,*,*])
  zf = reform(zfield[0,*,*])
  rf_x = reform(dx(rfield, i))
  rf_z = reform(dz(rfield, i))
  zf_x = reform(dx(zfield, i))
  zf_z = reform(dz(zfield, i))

  rmin = min(rf)
  rmax = max(rf)
  zmin = min(zf)
  zmax = max(zf)

  
  r = findgen(n)*(rmax-rmin)/(n-1.) + rmin
  z = findgen(n)*(zmax-zmin)/(n-1.) + zmin

  tol = 1e-4
  its = 100
  converged = 0
  iguess = n/2.
  for i=0, n-1 do begin
     jguess = n/2.
     for j=0, n-1 do begin

        iguess_old = iguess
        jguess_old = jguess

        for k=0, its-1 do begin
           rval = interpolate(rf, iguess, jguess)
           zval = interpolate(zf, iguess, jguess)
           d2 = (rval - r[i])^2 + (zval - z[j])^2
           if(d2 le tol^2) then begin
              i0[0,i,j] = iguess
              j0[0,i,j] = jguess
              converged = converged + 1
              goto, next
           end

           drdi = interpolate(rf_x, iguess, jguess)
           drdj = interpolate(rf_z, iguess, jguess)
           dzdi = interpolate(zf_x, iguess, jguess)
           dzdj = interpolate(zf_z, iguess, jguess)
           
           dd2di = (rval - r[i])*drdi + (zval - z[j])*dzdi
           dd2dj = (rval - r[i])*drdj + (zval - z[j])*dzdj
           
           ; d2 + vec(dl) . grad(d2) = 0
           ; vec(dl) = dl * grad(d2)/|grad(d2)|
           ; dl = -d2 grad(d2)/ | grad(d2) |^2
           dl = -d2 / (dd2di^2 + dd2dj^2)
           di = dl*dd2di
           dj = dl*dd2dj
           
                                ; check to make sure we're
                                ; staying in bounds.  If not, move
                                ; 90% of the way to the boundary
           if(iguess+di lt 0.) then begin
              di = -iguess*0.99
           endif else if(iguess+di gt n-1) then begin
              di = (n-1-iguess)*0.99
           endif
           if(jguess+dj lt 0.) then begin
              dj = -jguess*0.99
           endif else if(jguess+dj gt n-1) then begin
              dj = (n-1-jguess)*0.99
           endif

           iguess = iguess + di
           jguess = jguess + dj
        end
;        print, 'diverged at', r[i], z[j]
        iguess = iguess_old
        jguess = jguess_old
        i0[0,i,j] = -1
        j0[0,i,j] = -1
next:
     end
  end

  for i=1, n-2 do begin
     for j=1, n-2 do begin
        if(i0[0,i,j] eq -1) then begin
           interp = 0
           if(i0[0,i-1,j] ne -1 and i0[0,i+1,j] ne -1) then begin
              i0[0,i,j] = (i0[0,i-1,j]+i0[0,i+1,j])/2.
              j0[0,i,j] = (j0[0,i-1,j]+j0[0,i+1,j])/2.
              interp = 1
           end
           if(i0[0,i,j-1] ne -1 and i0[0,i,j+1] ne -1) then begin
              if(interp eq 1) then begin
                 i0[0,i,j] = (i0[0,i,j-1]+i0[0,i,j+1])/2.
                 j0[0,i,j] = (j0[0,i,j-1]+j0[0,i,j+1])/2.
              endif else begin
                 i0[0,i,j] = (i0[0,i,j-1]+i0[0,i,j+1] $
                              +  i0[0,i-1,j]+i0[0,i+1,j])/4.
                 j0[0,i,j] = (j0[0,i,j-1]+j0[0,i,j+1] $
                              +  j0[0,i-1,j]+j0[0,i+1,j])/4.
              endelse
           end
        end
     end
  end
  
  print, 'converged = ', converged
end
