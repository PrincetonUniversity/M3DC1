pro delta_W, filename=filename, _EXTRA=extra

  gamma = read_gamma(filename=filename,_EXTRA=extra)

  itor = read_parameter('itor',filename=filename,_EXTRA=extra)
  ntor = read_parameter('ntor',filename=filename,_EXTRA=extra)
  complex = read_parameter('icomplex',filename=filename,_EXTRA=extra)
  linear = read_parameter('linear',filename=filename,_EXTRA=extra)

  jr0 = read_field('jx',x,y,t,filename=filename,/equilibrium,$
                  _EXTRA=extra)
  jphi0 = read_field('jy',x,y,t,mesh=mesh,filename=filename,/equilibrium,$
                    _EXTRA=extra)
  jz0 = read_field('jz',x,y,t,mesh=mesh,filename=filename,/equilibrium,$
                  _EXTRA=extra)
  jr1 = read_field('jx',x,y,t,mesh=mesh,filename=filename, $
                   complex=complex,linear=linear,_EXTRA=extra)
  jphi1 = read_field('jy',x,y,t,mesh=mesh,filename=filename, $
                     complex=complex,linear=linear,_EXTRA=extra)
  jz1 = read_field('jz',x,y,t,mesh=mesh,filename=filename,$
                     complex=complex,linear=linear,_EXTRA=extra)
  br0 = read_field('bx',x,y,t,mesh=mesh,filename=filename,/equilibrium,$
                  _EXTRA=extra)
  bphi0 = read_field('by',x,y,t,mesh=mesh,filename=filename,/equilibrium,$
                    _EXTRA=extra)
  bz0 = read_field('bz',x,y,t,mesh=mesh,filename=filename,/equilibrium,$
                  _EXTRA=extra)
  br1 = read_field('bx',x,y,t,mesh=mesh,filename=filename,$
                  complex=complex,linear=linear,_EXTRA=extra)
  bphi1 = read_field('by',x,y,t,mesh=mesh,filename=filename,$
                  complex=complex,linear=linear,_EXTRA=extra)
  bz1 = read_field('bz',x,y,t,mesh=mesh,filename=filename,$
                  complex=complex,linear=linear,_EXTRA=extra)
  p1 = read_field('p',x,y,t,mesh=mesh,filename=filename,$
                  complex=complex,linear=linear,_EXTRA=extra)
  pr1 = read_field('p',x,y,t,mesh=mesh,filename=filename,$
                  complex=complex,linear=linear,_EXTRA=extra,op=2)
  pz1 = read_field('p',x,y,t,mesh=mesh,filename=filename,$
                  complex=complex,linear=linear,_EXTRA=extra,op=3)
  
  xir = read_field('xi_x',x,y,t,mesh=mesh,filename=filename, $
                  complex=complex,linear=linear,_EXTRA=extra)
  xiphi = read_field('xi_y',x,y,t,mesh=mesh,filename=filename, $
                  complex=complex,linear=linear,_EXTRA=extra)
  xiz = read_field('xi_z',x,y,t,mesh=mesh,filename=filename, $
                  complex=complex,linear=linear,_EXTRA=extra)

  den0 = read_field('den',x,y,t,mesh=mesh,filename=filename,/equilibrium,$
                   _EXTRA=extra)
  
  pphi1 = complex(0,ntor)*p1

  if(itor eq 1) then  begin
     r = radius_matrix(x,y,t)
  endif else r = 1.

  fr   = (jphi0*bz1 - jz0*bphi1 + jphi1*bz0 - jz1*bphi0 - pr1)
  fphi = (jz0*br1 - jr0*bz1 + jz1*br0 - jr1*bz0         - pphi1/r)
  fz   = (jr0*bphi1 - jphi0*br1 + jr1*bphi0 - jphi1*br0 - pz1)
  
  deltaW = -0.5*(conj(xir)*fr + conj(xiphi)*fphi + conj(xiz)*fz)

  ; factor of 1/2 from toroidal average
  Kx = 0.5*den0*(abs(xir)^2 + abs(xiphi)^2 + abs(xiz)^2)/2.  
  b2 = (abs(br1)^2 + abs(bphi1)^2 + abs(bz1)^2) / 2.

;  deltaW = read_field('delta_W',x,y,t,filename=filename, $
;                      mask=mask,_EXTRA=extra)
;  Kx = read_field('K',x,y,t,filename=filename, $
;                 mask=mask,_EXTRA=extra)
;  b2 = read_field('b2',x,y,t,filename=filename, $
;                  mask=mask,_EXTRA=extra)
  psi = read_field('psi',x,y,t,filename=filename,$
                   mask=mask,/equilibrium,_EXTRA=extra)

  psi1 = lcfs(flux0=psi0,filename=filename,_EXTRA=extra)
  psin = (psi - psi0) / (psi1 - psi0)

  Kx = Kx*(mask eq 0)
  deltaW = deltaW*(mask eq 0)
  b2 = b2*(mask eq 0)

  help, deltaW, x, y
  help, fr, fphi, fz, xir, xiphi, xiz
  contour_and_legend, real_part(deltaW), x, y

  dr = mean(deriv(x))
  dz = mean(deriv(y))

  K = 2.*!pi*total(Kx*r*dr*dz)
  dW = 2.*!pi*total(deltaW*r*dr*dz)
  dW_V = 2.*!pi*total((psin gt 1)*deltaW*r*dr*dz)
  dW_F = dW - dW_V

  dw_V_b2 = 2.*!pi*total((psin gt 1)*b2*r*dr*dz)

  contour_and_legend, real_part((psin gt 1)*deltaW), x, y
  
  print, 'K = ', K
  print, 'dW = ', dW
  print, 'dW_F = ', dW_F
  print, 'dW_V = ', dW_V
  print, 'dW_V from b^2 = ', dW_V_b2
  print, 'dW_F + dW_V_b2 = ', dW_F + dW_V_b2
  print, 'gamma = ', gamma
  print, 'gamma from dW / K = ', sqrt(-dW / K)
  print, 'gamma from (dW_F + dW_V_b2) / K = ', sqrt(-(dW_F + dW_V_b2) / K)
end
