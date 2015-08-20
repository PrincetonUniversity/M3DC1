pro flux_coordinates, _EXTRA=extra, pest=pest, points=pts, rpath=rpath, zpath=zpath, $
                      flux=psi, theta=theta, nflux=psi_norm, q=q, jacobian=jac, $
                      plot=makeplot

  if(keyword_set(pest)) then begin
     print, 'Creating flux coordinates using PEST angle'
  endif else begin
     print, 'Creating flux coordinates using GEOMETRIC angle'
  end

;  if(n_elements(bmncdf) eq 0) then bmncdf='eq.nc'

  psi0 = read_field('psi',x,z,t,slice=-1,points=pts,_EXTRA=extra)
  psi0_r = read_field('psi',x,z,t,slice=-1,points=pts,_EXTRA=extra,op=2)
  psi0_z = read_field('psi',x,z,t,slice=-1,points=pts,_EXTRA=extra,op=3)
  i0 = read_field('I',x,z,t,slice=-1,points=pts,_EXTRA=extra)

  n = pts  ; number of flux surfaces
  m = pts  ; number of poloidal points

  rpath = fltarr(m,n)
  zpath = fltarr(m,n)
  jac = fltarr(m,n)

  psi_s = lcfs(psi0, x, z, $
          axis=axis, xpoint=xpoint, flux0=flux0, _EXTRA=extra)

  psi = (psi_s - flux0)*(findgen(n)+.5)/n + flux0
  psi_norm = (psi - flux0)/(psi_s - flux0)
  theta = 2.*!pi*findgen(m)/m

  tol = 1e-4
  maxits = 20
  for i=0, m-1 do begin
     rho = 0.001
     for j=0, n-1 do begin
        ; do newton iterations to find (R,Z) at (psi, theta)
        converged = 0
        for k=0, maxits-1 do begin
           rpath[i,j] = rho*cos(theta[i]) + axis[0]
           zpath[i,j] = rho*sin(theta[i]) + axis[1]
        
           psi_x = field_at_point(psi0, x, z, rpath[i,j], zpath[i,j])
           psi_r = field_at_point(psi0_r, x, z, rpath[i,j], zpath[i,j])
           psi_z = field_at_point(psi0_z, x, z, rpath[i,j], zpath[i,j])
        
           if(abs(psi_x - psi[j]) lt tol) then begin
              converged = 1
              break
           endif
           
           dpsi_drho = psi_r*cos(theta[i]) + psi_z*sin(theta[i])
           drho = (psi[j] - psi_x)/dpsi_drho
           if(rho + drho lt 0.) then begin
              rho = rho / 2.
           endif else begin
              rho = rho + drho
           endelse
        end
        if(converged eq 0) then begin
           print, 'Error at (psi_n, theta) = ', psi_norm[j], theta[i]
           print, 'Did not converge after, ', maxits, ' iterations.', $
                  rho, rpath[i,j], zpath[i,j]
           rpath[i,j] = 0.
           zpath[i,j] = 0.
           rho = 0.001
        end
     end
  end

  ; find pest angle
  if(psi_s gt flux0) then begin
     fac = -1.
  endif else begin
     fac = 1.
  endelse
  q = fltarr(n)
  theta_pest = fltarr(m,n)
  integrand = fltarr(m)
  gradpsi = sqrt(psi0_r^2 + psi0_z^2)
  for j=0, n-1 do begin
     rx = [rpath[m-1,j],rpath[*,j],rpath[0,j]]
     zx = [zpath[m-1,j],zpath[*,j],zpath[0,j]]
     dlx = sqrt(deriv(rx)^2 + deriv(zx)^2)
     dl = dlx[1:m]

     for i=0, m-1 do begin
        bp = fac*field_at_point(gradpsi,x,z,rpath[i,j],zpath[i,j]) $
             / rpath[i,j]
        ix = field_at_point(i0,x,z,rpath[i,j],zpath[i,j])
        integrand[i] = ix/(rpath[i,j]^2*bp)

        if(i eq 0) then begin
           theta_pest[i,j] = 0.
        endif else begin
           theta_pest[i,j] = theta_pest[i-1,j] $
                             +dl[i]*integrand[i]/2. $
                             +dl[i-1]*integrand[i-1]/2.
        end
     end
     q[j] = total(integrand*dl)/(2.*!pi)
     theta_pest[*,j] = theta_pest[*,j]/q[j]
  end

  if(keyword_set(pest)) then begin
     for j=0, n-1 do begin
        newm = interpol(findgen(m),theta_pest[*,j],theta)
        rpath[*,j] = interpolate(rpath[*,j],newm)
        zpath[*,j] = interpolate(zpath[*,j],newm)

        ; we have analytic expression for PEST Jacobian
        for i=0,m-1 do begin
           jac[i,j] = rpath[i,j]^2*q[j] / $
                 field_at_point(i0,x,z,rpath[i,j],zpath[i,j])
        end
     end
     
  endif else begin
     ; calculate jacobian
     dr_dpsi = rpath
     dz_dpsi = zpath
     dr_dtheta = rpath
     dz_dtheta = zpath
     for i=0, m-1 do begin
        dr_dpsi[i,*] = deriv(psi, rpath[i,*])
        dz_dpsi[i,*] = deriv(psi, zpath[i,*])
     end
     for j=0, n-1 do begin
        dr_dtheta[*,j] = deriv(theta, rpath[*,j])
        dz_dtheta[*,j] = deriv(theta, zpath[*,j])
     end
     jac = rpath*(dr_dtheta*dz_dpsi - dr_dpsi*dz_dtheta)
  endelse
  
  if(keyword_set(makeplot)) then begin

     window, 1
     ; plot jacobian
;     contour, alog10(abs(jac)), theta, psi_norm, $
;              xtitle='!7h!X', ytitle='!7W!X', $
;              xrange=[0,2*!pi], xstyle=1, /follow
     contour_and_legend, jac, theta, psi_norm, table=39, $
              xtitle='!7h!X', ytitle='!7W!X', /lines

     window, 0
     !p.multi = [0,2,1]

     ; plot coordinates
     plot, [min(rpath), max(rpath)], [min(zpath), max(zpath)], $
           /nodata, /iso, xtitle='!8R!X', ytitle='!8Z!X'
     for j=0, n-1, 10 do begin
        oplot, [rpath[*,j],rpath[0,j]], [zpath[*,j],zpath[0,j]]
     end
     oplot, [rpath[*,n-1],rpath[0,n-1]], [zpath[*,n-1],zpath[0,n-1]]
     for i=0, m-1, 10 do begin
        oplot, rpath[i,*], zpath[i,*]
     end

     ; plot q
     plot, psi_norm, abs(q), xtitle='!7W!X', ytitle='!8q!X'
     
     !p.multi=0
  end
end
