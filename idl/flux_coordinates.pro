;==================================================================
; flux_coordinates
; ~~~~~~~~~~~~~~~~
; This function calculates flux coordinates and returns a structure
; containing a description of the coordinate system
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
                           njac=njac

  if(n_elements(fast0) eq 0) then fast0 = 0
  fast = fast0

  if(n_elements(pts) eq 0) then begin
     if(n_elements(x) eq 0 or n_elements(z) eq 0) then begin
        pts=200
     endif else begin
        pts=sqrt(n_elements(x)*n_elements(z))
     endelse
  end
  if(n_elements(fbins) eq 0) then fbins=200
  if(n_elements(tbins) eq 0) then tbins=200

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
     psi0 = read_field('psi',x,z,t,points=pts,/equilibrium,_EXTRA=extra)
  end

  if(keyword_set(fast)) then begin
     psi0_r = dx(psi0, x)
     psi0_z = dz(psi0, z)
  endif else begin
     psi0_r = read_field('psi',x,z,t,points=pts,/equilibrium,_EXTRA=extra,op=2)
     psi0_z = read_field('psi',x,z,t,points=pts,/equilibrium,_EXTRA=extra,op=3)
     if(n_elements(i0) eq 0) then begin
        i0 = read_field('I',x,z,t,points=pts,/equilibrium,_EXTRA=extra)
     end
  endelse

  psi_s = lcfs(psi0, x, z, $
          axis=axis, xpoint=xpoint, flux0=flux0, _EXTRA=extra)

  m = tbins  ; number of poloidal points
  n = fbins  ; number of radial points
  
  psi = (psi_s - flux0)*(findgen(n)+.5)/n + flux0
  psi_norm = (psi - flux0)/(psi_s - flux0)

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

  tol_psi = 1e-4*abs(psi_s - flux0)
  tol_r = 1e-4*(max(x) - min(x))
  maxits = 20
  rho_old = 0.
  psi_old = 0.

  tot = long(0)

  for i=0, m-1 do begin
     co = cos(theta[i])
     sn = sin(theta[i])
     for j=0, n-1 do begin
        ; do newton iterations to find (R,Z) at (psi, theta)
        converged = 0
        if(j eq 0) then begin
           rho = 0.001
        endif else begin
           rho = rho + (psi[j]-psi[j-1])/dpsi_drho
        endelse
        for k=0, maxits-1 do begin
           rpath[i,j] = rho*co + axis[0]
           zpath[i,j] = rho*sn + axis[1]
        
           psi_x = field_at_point(psi0, x, z, rpath[i,j], zpath[i,j])

           if(abs(psi_x - psi[j]) lt tol_psi) then begin
              converged = 1
              break
           endif

           if(keyword_set(fast)) then begin
              if((k eq 0 and j eq 0) $
                 or abs(psi_x - psi_old) lt tol_psi $
                 or abs(rho - rho_old) lt tol_r) then begin
                 psi_r = field_at_point(psi0_r, x, z, rpath[i,j], zpath[i,j])
                 psi_z = field_at_point(psi0_z, x, z, rpath[i,j], zpath[i,j])

                 dpsi_drho = psi_r*co + psi_z*sn
              endif else begin
                 dpsi_drho = (psi_x - psi_old) / (rho - rho_old)
              endelse
              rho_old = rho
              psi_old = psi_x

           endif else begin
              psi_r = field_at_point(psi0_r, x, z, rpath[i,j], zpath[i,j])
              psi_z = field_at_point(psi0_z, x, z, rpath[i,j], zpath[i,j])
                   
              dpsi_drho = psi_r*co + psi_z*sn
           end
           drho = (psi[j] - psi_x)/dpsi_drho
           if(rho + drho lt 0.) then begin
              rho = rho / 2.
           endif else begin
              rho = rho + drho
           endelse
        end
        tot = tot + k
        if(converged eq 0) then begin
           print, 'Error at (psi_n, theta) = ', psi_norm[j], theta[i]
           print, 'Did not converge after, ', maxits, ' iterations.', $
                  rho, rpath[i,j], zpath[i,j]
           rho = 0.001
        end
     end
  end

  print, 'Total Newton iterations ', tot
  print, 'Average Newton iterations ', float(tot)/float(long(m)*long(n))

  ; find pest angle
  if(psi_s lt flux0) then begin
     print, 'grad(psi) is inward'
     fac = -1.
  endif else begin
     print, 'grad(psi) is outward'
     fac = 1.
  endelse
  theta_sfl = fltarr(m,n)
  f = fltarr(n)
  fjr2 = fltarr(m)
  dthetadl = fltarr(m)
  bp = fltarr(m)
  gradpsi = sqrt(psi0_r^2 + psi0_z^2)
  for j=0, n-1 do begin
     rx = [rpath[m-1,j],rpath[*,j],rpath[0,j]]
     zx = [zpath[m-1,j],zpath[*,j],zpath[0,j]]
     dlx = sqrt(deriv(rx)^2 + deriv(zx)^2)
     dl = dlx[1:m]
     
     for i=0, m-1 do begin
        bp[i] = fac*field_at_point(gradpsi,x,z,rpath[i,j],zpath[i,j]) $
                / rpath[i,j]

        if(not keyword_set(fast)) then begin
           ix = field_at_point(i0,x,z,rpath[i,j],zpath[i,j])
           ; dtheta/dl = -1./(Bp*Jac)
           if(keyword_set(pest)) then begin
              dthetadl[i] = -ix/(rpath[i,j]^2*bp[i])
           endif else if(keyword_set(boozer)) then begin
              dthetadl[i] = -(bp[i]^2 + (ix/rpath[i,j])^2) / bp[i]
           endif else if(keyword_set(hamada)) then begin
              dthetadl[i] = -1./bp[i]
           endif
           fjr2[i] = -ix/(rpath[i,j]^2*bp[i])
           
           if(i eq 0) then begin
              theta_sfl[i,j] = 0.
           endif else begin
              theta_sfl[i,j] = theta_sfl[i-1,j] $
                                +dl[i]*dthetadl[i]/2. $
                                +dl[i-1]*dthetadl[i-1]/2.
           end
        end
     end
     
     area[j] = 2.*!pi*total(dl*rpath[*,j])
     dV[j] = 2.*!pi*total(dl/bp)
     if(j eq 0) then begin
        V[j] = dV[j]*(psi[j]-flux0)
     endif else begin
        V[j] = V[j-1] + (dV[j]+dV[j-1])*(psi[j]-psi[j-1])/2.
     endelse

     ; calculate q
     if(not keyword_set(fast)) then begin
        ; normalize theta to 2 pi
        if(geo eq 0) then begin
           f[j] = total(dthetadl*dl)/(2.*!pi)
           theta_sfl[*,j] = theta_sfl[*,j]/f[j]
        end
        q[j] = total(fjr2*dl)/(2.*!pi)
           
        ; calculate toroidal flux
        if(j eq 0) then begin
           phi[j] = -2.*!pi*q[j]*(psi[j]-flux0)
        endif else begin
           phi[j] = phi[j-1] - 2.*!pi*(q[j]+q[j-1])*(psi[j]-psi[j-1])/2.
        endelse
     end
                
     ; sanity checks
     if(V[j] lt 0) then begin
        print, 'ERROR, volume is negative'
        print, j, V[j], dV[j], psi[j]-flux0, bp[j]
;        return, 0
     end
     if(not keyword_set(fast)) then begin
        if(theta_sfl[0,j] gt theta_sfl[m-1,j]) then begin
           print, 'ERROR, theta_sfl is clockwise'
           return, 0
        end
        if(fac*q[j]*ix gt 0.) then begin
           print, 'ERROR, q has wrong sign'
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
        
        ; use analytic expression for Jacobian
        if(keyword_set(pest)) then begin
           for i=0,m-1 do begin
              jac[i,j] = rpath[i,j]^2*q[j] / $
                         field_at_point(i0,x,z,rpath[i,j],zpath[i,j])
           end
        endif else if(keyword_set(boozer)) then begin
           b2r2 = i0^2 + psi0_r^2 + psi0_z^2
           for i=0,m-1 do begin
              jac[i,j] = f[j]*rpath[i,j]^2 / $
                         field_at_point(b2r2,x,z,rpath[i,j],zpath[i,j])
           end
        endif else if(keyword_set(hamada)) then begin
           jac[*,j] = f[j]
        end
     end
     
  endif 

  ; Calculate Jacobian if requested (or if analytic expression isn't available)
  if(geo eq 1 or keyword_set(njac)) then begin
     print, 'Calculating Jacobian numerically'
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
  end

  ; calculate deviation of toroidal angle from geometric toroidal angle
  if(not keyword_set(pest) and not keyword_set(geo) $
     and not keyword_set(fast)) then begin
     for j=0, n-1 do begin
        for i=0, m-1 do begin
           jac_pest = rpath[i,j]^2*q[j] / $
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

  ; define flux coordinates structure
  fc = { m:m, n:n, r:rpath, z:zpath, r0:axis[0], z0:axis[1], omega:omega, $
         psi1:psi_s, psi0:flux0, psi:psi, psi_norm:psi_norm, theta:theta, $
         j:jac, q:q, area:area, dV:dV, pest:keyword_set(pest), $
         V:V, phi:phi, phi_norm:phi/phi[n-1], rho:sqrt(phi/phi[n-1]) }

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
;     contour, alog10(abs(jac)), theta, psi_norm, $
;              xtitle='!7h!X', ytitle='!7W!X', $
;              xrange=[0,2*!pi], xstyle=1, /follow
     contour_and_legend, jac, theta, psi_norm, table=39, $
                         xtitle='!7h!X', ytitle='!7W!X', /lines

     ; plot q
     window, 0
     !p.multi = [0,2,2]
     plot, psi_norm, q, xtitle='!7W!X', ytitle='!8q!X'
     plot, psi_norm, V, xtitle='!7W!X', ytitle='!8V!X'
;     plot, psi_norm, area, xtitle='!7W!X', ytitle='!8A!X'
     plot, psi_norm, omega, xtitle='!7W!X', ytitle='!7x!X'
     plot, psi_norm, phi, xtitle='!7W!X', ytitle='!7u!D!8T!N!X'
     !p.multi=0
  end

  return, fc
end
