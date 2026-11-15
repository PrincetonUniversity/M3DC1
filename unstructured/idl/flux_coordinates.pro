;==================================================================
; flux_coordinates
; ~~~~~~~~~~~~~~~~
; This function calculates flux coordinates and returns a structure
; containing a description of the coordinate system 
; (psi_norm, theta, zeta)
; psi_norm increases outward from the axis, 
; zeta increases counter-clockwise looking down from above
; theta increases clockwise around magnetic axis
;
; INPUT
;  /pest     : use PEST angle (default is geometric angle)
;  /makeplot : plot some of the coordinate system data
;  /fast     : don't calculate things involving toroidal field
;  points    : number of points (radial and poloidal) to use 
; OUTPUT
;  .m        : number of poloidal points
;  .n        : number of radial points
;  
;==================================================================

function flux_coordinates, _EXTRA=extra, pest=pest, points=pts, $
                           plot=makeplot, fast=fast0, psi0=psi0, $
                           i0=i0, x=x, z=z, fbins=fbins, $
                           tbins=tbins, boozer=boozer, hamada=hamada, $
                           njac=njac, itor=itor, psin_range=psin_range, $
                           r0=r0, filename=filename, $
                           dpsi0_dx=psi0_r, dpsi0_dz=psi0_z
  if(n_elements(filename) eq 0) then begin
     fn = 'C1.h5'
  endif else begin
     fn = filename[0]
  end

  if(n_elements(fast0) eq 0) then fast0 = 0
  if(n_elements(itor) ne 1) then $
     itor = read_parameter('itor', file=fn, _EXTRA=extra)
  if(n_elements(r0) ne 1) then $
     r0 = read_parameter('rzero', file=fn, _EXTRA=extra)
  if(r0 eq 0) then begin
     print, 'Error: r0 = 0. Using r0 = 1'
     r0 = 1.
  end
  if(itor eq 0) then begin
     period = 2.*!pi*r0
  endif else begin
     period = 2.*!pi
  end
  help, r0, itor, period
  fast = fast0

  if(n_elements(pts) eq 0) then begin
     if(n_elements(x) eq 0 or n_elements(z) eq 0) then begin
        pts=200
     endif else begin
        pts=sqrt(n_elements(x)*n_elements(z))
     endelse
  end
  if(n_elements(fbins) eq 0) then fbins=pts
  if(n_elements(tbins) eq 0) then tbins=pts

  geo = 0
  if(keyword_set(pest)) then begin
     print, 'Creating flux coordinates using PEST angle'
     fast = 0
  endif else if(keyword_set(boozer)) then begin
     print, 'Creating flux coordinates using BOOZER angle'
     fast = 0
  endif else if(keyword_set(hamada)) then begin
     print, 'Creating flux coordinates using HAMADA angle'
     fast = 0
  endif else begin
     print, 'Creating flux coordinates using GEOMETRIC angle'
     geo = 1
  end
  print, 'FAST MODE: ', keyword_set(fast)
  print, 'Using FC resolution ', tbins, fbins

  if(n_elements(psi0) eq 0 or n_elements(x) eq 0 or n_elements(z) eq 0) then begin
     print, 'READING PSI IN FLUX_COORDINATES'
     psi0 = read_field('psi',x,z,t,points=pts,/equilibrium,$
                       filename=fn,_EXTRA=extra)
  end

  if(keyword_set(fast)) then begin
     if(n_elements(psi0_r) eq 0 or n_elements(psi0_z) eq 0) then begin
        psi0_r = dx(psi0, x)
        psi0_z = dz(psi0, z)
     end
  endif else begin
     if(n_elements(psi0_r) eq 0 or n_elements(psi0_z) eq 0) then begin
        psi0_r = read_field('psi',x,z,t,points=pts,/equilibrium,$
                            filename=fn,_EXTRA=extra,op=2)
        psi0_z = read_field('psi',x,z,t,points=pts,/equilibrium,$
                            filename=fn,_EXTRA=extra,op=3)
     end
        if(n_elements(i0) eq 0) then begin
           i0 = read_field('I',x,z,t,points=pts,/equilibrium,$
                           filename=fn,_EXTRA=extra)
        end
  endelse

  psi_s = lcfs(psi0, x, z, axis=axis, xpoint=xpoint, flux0=flux0, $
               /refine, filename=fn, _EXTRA=extra)

  print, 'Magnetic axis = ', axis
  print, 'flux0, psi_s = ', flux0, psi_s

  m = tbins  ; number of poloidal points
  n = fbins  ; number of radial points
  
  if(n_elementS(psin_range) eq 0) then begin
     psi = (psi_s - flux0)*(findgen(n)+.5)/n + flux0
     psin = (psi - flux0)/(psi_s - flux0)
  endif else begin
     dpsi = psi_s - flux0
     p0 = flux0 + (psi_s - flux0)*psin_range[0]
     p1 = flux0 + (psi_s - flux0)*psin_range[1]
     psi = (p1 - p0)*(findgen(n)+.5)/n + p0
     psin = (psi - flux0)/(psi_s - flux0)
  endelse
  dpsi_dpsin = psi_s - flux0
  print, 'dpsi_dpsin = ', dpsi_dpsin

  psin0 = (psi0 - flux0) / dpsi_dpsin
  psin0_r = psi0_r / dpsi_dpsin
  psin0_z = psi0_z / dpsi_dpsin

  ; Remove too large target psin that cannot be found along some poloidal angles
  ; Ensure target psin values form closed flux surfaces within a given rectangular boundary
  min_psin0_at_boundary= min([min([psin0[0,-1,*],psin0[0,0,*]]),min([psin0[0,*,0],psin0[0,*,-1]])])
  if(psin[-1] gt min_psin0_at_boundary) then begin
    psin = psin[0] + (0.999*min_psin0_at_boundary - psin[0]) *FINDGEN(n)/(n-1)
    psi = psin*(psi_s -flux0) +flux0
  endif

  theta = 2.*!pi*findgen(m)/m
  rpath = fltarr(m,n)
  zpath = fltarr(m,n)
  jac = fltarr(m,n)
  omega = fltarr(m,n)
  q = fltarr(n)
  dV = fltarr(n)
  V = fltarr(n)
  phi = fltarr(n)
  area = fltarr(n)
  current = fltarr(n)

  tol_psin = 1e-4
  tol_r = 1e-4*(max(x) - min(x))
  maxits = 20
  rho_old = 0.
  psin_old = 0.

  tot = long(0)

  ; find points on magnetic surfaces
  for i=0, m-1 do begin
     ; minus sign is because theta increases clockwise about axis
     co = cos(-theta[i])
     sn = sin(-theta[i])
     dpsin_drho = 0.
     max_drho = sqrt(((x[n_elements(x)-1] - x[0])*co)^2 + $
                     ((z[n_elements(z)-1] - z[0])*sn)^2) * 0.1
     for j=0, n-1 do begin
        ; do newton iterations to find (R,Z) at (psin, theta)
        converged = 0
        if(j eq 0 or dpsin_drho eq 0.) then begin
           old_psin_x=0.0
           rho = 0.01
           drho = rho
        endif else begin
           rho = rho + (psin[j]-psin[j-1])/dpsin_drho
           drho = (psin[j]-psin[j-1])/dpsin_drho
           old_psin_x = psin[j-1]
        endelse
        for k=0, maxits-1 do begin
           rpath[i,j] = rho*co + axis[0]
           zpath[i,j] = rho*sn + axis[1]
        
           psin_x = field_at_point(psin0, x, z, rpath[i,j], zpath[i,j])

           if(abs(psin_x - psin[j]) lt tol_psin) then begin
              converged = 1
              break
           endif

           if(keyword_set(fast)) then begin
              if((k eq 0 and j eq 0) $
                 or abs(psin_x - psin_old) lt tol_psin $
                 or abs(rho - rho_old) lt tol_r) then begin
                 psin_r = field_at_point(psin0_r, x, z, rpath[i,j], zpath[i,j])
                 psin_z = field_at_point(psin0_z, x, z, rpath[i,j], zpath[i,j])

                 dpsin_drho = psin_r*co + psin_z*sn
              endif else begin
                 dpsin_drho = (psin_x - psin_old) / (rho - rho_old)
              endelse
              rho_old = rho
              psin_old = psin_x

           endif else begin
              psin_r = field_at_point(psin0_r, x, z, rpath[i,j], zpath[i,j])
              psin_z = field_at_point(psin0_z, x, z, rpath[i,j], zpath[i,j])
                   
              dpsin_drho = psin_r*co + psin_z*sn
           end
           if(psin[j] lt 0.95) then begin
             ; Do usual Newton's method
             drho = (psin[j] - psin_x)/dpsin_drho
           endif else begin
             ; Do the Bisection method for a more reliable root finding
             if( (old_psin_x-psin[j])*(psin_x-psin[j])<0 ) then begin
                 drho=-0.5*drho
             endif
           endelse
           old_psin_x = psin_x

           if(rho + drho lt 0.) then begin
              rho = rho / 2.
           endif else if(drho gt max_drho) then begin
              rho = rho + max_drho
           endif else begin
              rho = rho + drho
           endelse
        end
        tot = tot + k
        if(converged eq 0) then begin
           print, 'Error at (psi_n, theta) = ', psin[j], theta[i]
           print, 'Did not converge after, ', maxits, ' iterations.', $
                  rho, rpath[i,j], zpath[i,j]
           rho = 0.01
        end
     end
  end

  print, 'Total Newton iterations ', tot
  print, 'Average Newton iterations ', float(tot)/float(long(m)*long(n))

;  find pest angle
  if(psi_s lt flux0) then begin
     print, 'Ip > 0'
     print, 'grad(psi) is inward'
     fac = -1.
  endif else begin
     print, 'Ip < 0'
     print, 'grad(psi) is outward'
     fac = 1.
  endelse
  theta_sfl = fltarr(m,n)
  f = fltarr(n)
  fjr2 = fltarr(m)
  dthetadl = fltarr(m)
  bp = fltarr(m)
  gradpsi = sqrt(psi0_r^2 + psi0_z^2)

  rp = rpath
  if(itor eq 0) then rp[*,*] = 1.

  for j=0, n-1 do begin
     rx = [rpath[m-1,j],rpath[*,j],rpath[0,j]]
     zx = [zpath[m-1,j],zpath[*,j],zpath[0,j]]
     drx = deriv(rx)
     dzx = deriv(zx)
     dr = drx[1:m]
     dz = dzx[1:m]
     dl = sqrt(dr^2 + dz^2)
     br = -field_at_point(psi0_z,x,z,rpath[*,j],zpath[*,j])/rp[*,j]
     bz =  field_at_point(psi0_r,x,z,rpath[*,j],zpath[*,j])/rp[*,j]
     
     for i=0, m-1 do begin

        bp[i] = field_at_point(gradpsi,x,z,rpath[i,j],zpath[i,j])/rp[i,j]

        if(not keyword_set(fast)) then begin
           ix = field_at_point(i0,x,z,rpath[i,j],zpath[i,j])
           ; dtheta/dl ~ 1./(Bp*Jac)
           if(keyword_set(pest)) then begin
              dthetadl[i] = -fac*ix/(rp[i,j]^2*bp[i])
           endif else if(keyword_set(boozer)) then begin
              dthetadl[i] = -fac*(bp[i]^2 + (ix/rp[i,j])^2) / bp[i]
           endif else if(keyword_set(hamada)) then begin
              dthetadl[i] = -fac*1./bp[i]
           endif
           fjr2[i] = -fac*ix/(rp[i,j]^2*bp[i])
           
           if(i eq 0) then begin
              theta_sfl[i,j] = 0.
           endif else begin
              theta_sfl[i,j] = theta_sfl[i-1,j] $
                                +dl[i]*dthetadl[i]/2. $
                                +dl[i-1]*dthetadl[i-1]/2.
           end
        end
     end
     
     current[j] = total(br*dr + bz*dz)
     area[j] = period*total(dl*rp[*,j])
     dV[j] = period*total(fac*dpsi_dpsin*dl/bp)
     if(j eq 0) then begin
        V[j] = dV[j]*psin[j]/2.
     endif else begin
        V[j] = V[j-1] + (dV[j]+dV[j-1])*(psin[j]-psin[j-1])/2.
     endelse

     ; calculate q
     if(not keyword_set(fast)) then begin
        ; normalize theta to 2 pi
        if(geo eq 0) then begin
           f[j] = total(dthetadl*dl)/(2.*!pi)
           theta_sfl[*,j] = theta_sfl[*,j]/f[j]
        end
        q[j] = total(fjr2*dl)/period
           
        ; calculate toroidal flux
        if(j eq 0) then begin
           phi[j] = -period*q[j]*(psi[j]-flux0)
        endif else begin
           phi[j] = phi[j-1] - period*(q[j]+q[j-1])*(psi[j]-psi[j-1])/2.
        endelse
     end
                
     ; sanity checks
     if(V[j] lt 0) then begin
        print, 'ERROR, volume is negative'
        print, j, V[j], dV[j], psi[j]-flux0, bp[j]
        return, 0
     end
     if(not keyword_set(fast)) then begin
        if(theta_sfl[0,j] gt theta_sfl[m-1,j]) then begin
           print, 'ERROR, theta_sfl is clockwise'
           return, 0
        end
        if(fac*q[j]*ix gt 0.) then begin
           print, 'ERROR, q has wrong sign'
           stop
           return, 0
        end
        if(phi[j]*ix lt 0.) then begin
           print, 'ERROR, phi has wrong sign'
           return, 0
        end
     end
  end
     
  if(geo eq 0) then begin
     for j=0, n-1 do begin
        ; interpolate fields to be evenly spaced in PEST angle
        newm = interpol(findgen(m),theta_sfl[*,j],theta)
        rpath[*,j] = interpolate(rpath[*,j],newm)
        zpath[*,j] = interpolate(zpath[*,j],newm)
        if(itor eq 1) then rp[*,j] = rpath[*,j]
        
        ; use analytic expression for Jacobian
        if(keyword_set(pest)) then begin
           for i=0,m-1 do begin
              jac[i,j] = -dpsi_dpsin*rp[i,j]^2*q[j] / $
                         field_at_point(i0,x,z,rpath[i,j],zpath[i,j])
           end
        endif else if(keyword_set(boozer)) then begin
           b2r2 = i0^2 + psi0_r^2 + psi0_z^2
           for i=0,m-1 do begin
              jac[i,j] = -dpsi_dpsin*f[j]*rp[i,j]^2 / $
                         field_at_point(b2r2,x,z,rpath[i,j],zpath[i,j])
           end
        endif else if(keyword_set(hamada)) then begin
           jac[*,j] = -dpsi_dpsin*f[j]
        end
     end
     
  endif 

  ;  Calculate Jacobian if requested (or if analytic expression isn't available)

  print, 'Calculating Jacobian numerically'
                                ; calculate jacobian
  dr_dpsi = rpath
  dz_dpsi = zpath
  dr_dtheta = rpath
  dz_dtheta = zpath
  for i=0, m-1 do begin
     dr_dpsi[i,*] = deriv(psin, rpath[i,*])
     dz_dpsi[i,*] = deriv(psin, zpath[i,*])
  end
  for j=0, n-1 do begin
     dr_dtheta[*,j] = deriv(theta, rpath[*,j])
     dz_dtheta[*,j] = deriv(theta, zpath[*,j])
  end
  jac_test = -(dr_dpsi*dz_dtheta - dr_dtheta*dz_dpsi)
  if(itor eq 1) then jac_test = jac_test*rpath

  if(mean(jac_test) lt 0) then begin
     print, 'ERROR: numerical jacobian is negative!'
     return, 0
  end
  if(geo eq 1 or keyword_set(njac)) then jac = jac_test

  if(mean(jac) lt 0) then begin
     print, 'ERROR: jacobian is negative!'
     return, 0
  end

  ; calculate deviation of toroidal angle from geometric toroidal angle
  ; omega = zeta - phi
  if(not keyword_set(geo) $
;  if(not keyword_set(pest) and not keyword_set(geo) $
     and not keyword_set(fast)) then begin
     ; q B.Grad(theta) = q*dpsi_dpsin/Jac = B.Grad(zeta)
     ; q B.Grad(theta_pest) = q*dpsi_dpsin/Jac_pest = B.Grad(phi) = I/R^2
     ; B.Grad(zeta) - B.Grad(phi) = 1/Jac - 1/Jac_pest = (1-Jac/Jac_pest)/Jac

     ; d(zeta - phi) = (1 - Jac/Jac_pest)/Jac dl
     ; dtheta = psi'/(Jac*Bp) dl
     ; d(zeta - phi) = (Bp/psi') * (1 - Jac/Jac_pest)

     for j=0, n-1 do begin
        for i=0, m-1 do begin
           jac_pest = -dpsi_dpsin*rp[i,j]^2*q[j] / $
                      field_at_point(i0,x,z,rpath[i,j],zpath[i,j])
           if(i eq 0) then begin
              omega[i,j] = 0.
           endif else begin
              dtheta = theta[i] - theta[i-1]
              omega[i,j] = omega[i-1,j] + dtheta* $
                                ((1. - jac[i  ,j]/jac_pest    )/2. $
                                +(1. - jac[i-1,j]/jac_pest_old)/2.)
           end
           jac_pest_old = jac_pest
        end
        omega[*,j] = omega[*,j]*q[j]
     end
  end

  ; Define flux coordinates structure
  ; psi = poloidal flux
  ; phi = toroidal flux
  fc = { m:m, n:n, r:rpath, z:zpath, r0:axis[0], z0:axis[1], omega:omega, $
         psi:psi, psi_norm:psin, flux_pol:-period*(psi-psi[0]), theta:theta, $
         j:jac, q:q, area:area, dV_dchi:dV, pest:keyword_set(pest), $
         boozer:keyword_set(boozer), hamada:keyword_set(hamada), $
         V:V, flux_tor:phi, phi_norm:phi/phi[n-1], rho:sqrt(phi/phi[n-1]), $
         period:period, itor:itor, current:current, dpsi_dchi:dpsi_dpsin}

  if(keyword_set(makeplot)) then begin

     ; plot coordinates
     window, 2
     plot, [min(rpath), max(rpath)], [min(zpath), max(zpath)], $
           /nodata, /iso, xtitle='!8R!X', ytitle='!8Z!X'
     for j=0, n-1, 10 do begin
        oplot, [rpath[*,j],rpath[0,j]], [zpath[*,j],zpath[0,j]], thick=1+(j eq 10)
     end
     oplot, [rpath[*,n-1],rpath[0,n-1]], [zpath[*,n-1],zpath[0,n-1]]
     for i=0, m-1, 10 do begin
        oplot, rpath[i,*], zpath[i,*], thick=1+(i eq 10)
     end

     window, 1
                                ; plot jacobian
;     contour, alog10(abs(jac)), theta, psin, $
;              xtitle='!7h!X', ytitle='!7W!X', $
;              xrange=[0,2*!pi], xstyle=1, /follow
     contour_and_legend, jac, theta, psin, table=39, $
                         xtitle='!7h!X', ytitle='!7W!X', /lines, $
                         _EXTRA=extra

     ; plot q
     window, 0
     !p.multi = [0,2,3]
     plot, psin, q, xtitle='!7W!X', ytitle='!8q!X'
     plot, psin, V, xtitle='!7W!X', ytitle='!8V!X'
     plot, psin, area, xtitle='!7W!X', ytitle='!8A!X'
     plot, psin, current, xtitle='!7W!X', ytitle='!8I!X'
     plot, theta, omega[*,20], xtitle='!7W!X', ytitle='!7x!X'
     plot, psin, phi, xtitle='!7W!X', ytitle='!7u!D!8T!N!X'
     !p.multi=0
  end

  return, fc
end
