pro write_geqdsk, eqfile=eqfile, $
                  psilim=psilim, points=pts, _EXTRA=extra
  
  if(n_elements(slice) eq 0) then begin
      slice = read_parameter('ntime', _EXTRA=extra) - 1
  end
  if(n_elements(eqfile) eq 0) then eqfile = 'geqdsk.out'
  get_normalizations, b0=b0, n0=n0, l0=l0, _EXTRA=extra

  bzero = read_parameter('bzero', _EXTRA=extra)
  rzero = read_parameter('rzero', _EXTRA=extra)

  ; calculate flux averages
  psi = read_field('psi',x,z,t,mesh=mesh,/equilibrium,_EXTRA=extra)
  psi_r = read_field('psi',x,z,t,mesh=mesh,/equilibrium,op=2,_EXTRA=extra)
  psi_z = read_field('psi',x,z,t,mesh=mesh,/equilibrium,op=3,_EXTRA=extra)
  psi_lp= read_field('psi',x,z,t,/equilibrium,op=7,_EXTRA=extra)
  p0 = read_field('p',x,z,t,/equilibrium,_EXTRA=extra)
  p0_r = read_field('p',x,z,t,/equilibrium,op=2,_EXTRA=extra)
  p0_z = read_field('p',x,z,t,/equilibrium,op=3,_EXTRA=extra)
  I0 = read_field('I',x,z,t,/equilibrium,_EXTRA=extra)
  I0_r = read_field('I',x,z,t,/equilibrium,op=2,_EXTRA=extra)
  I0_z = read_field('I',x,z,t,/equilibrium,op=3,_EXTRA=extra)
  r = radius_matrix(x,z,t)
  beta = r^2*2.*p0/(s_bracket(psi,psi,x,z) + I0^2)
  beta0 = mean(2.*p0*r^2/(bzero*rzero)^2)
  b2 = (s_bracket(psi,psi,x,z) + I0^2)/r^2
  dx = (max(x)-min(x))/(n_elements(x) - 1.)
  dz = (max(z)-min(z))/(n_elements(z) - 1.)
  jphi = psi_lp - psi_r/r
  tcur = read_scalar('ip',_EXTRA=extra,/mks)
  print, 'current = ', tcur[0]

  ; calculate magnetic axis and xpoint
  lcfs_psi = lcfs(psi,x,z, axis=axis, xpoint=xpoint, $
                  flux0=flux0, _EXTRA=extra)

  ; plot psi
  contour_and_legend, psi, x, z

  ; calculate wall points
  bound_xy = get_boundary_path(mesh=mesh, _EXTRA=extra)
  nwall = n_elements(reform(bound_xy[0,*]))
  rwall = fltarr(nwall)
  zwall = fltarr(nwall)
  rwall[*] = bound_xy[0,*]
  zwall[*] = bound_xy[1,*]

  ifixedb = read_parameter('ifixedb', _EXTRA=extra)
  print, 'ifixedb = ', ifixedb

  if(ifixedb eq 1) then begin
      psilim = 0.
      lcfs_psi = 0.
      
      ; for ifixedb eq 1, boundary points are same as wall points
      nlim = nwall
      rlim = rwall
      zlim = zwall

  endif else begin
      ; find points at psi = psilim to use as boundary points
      if(n_elements(psilim) eq 0) then begin
          ; use psilim = lcfs if psilim is not given
          psilim = lcfs_psi
      endif
      print, 'lcfs_psi = ', lcfs_psi
      print, 'psilim = ', psilim
      
      window, 0
      lcfs_xy = path_at_flux(psi,x,z,t,psilim)

      ; count only points on separatrix above the xpoint
      if(n_elements(xpoint) gt 1) then begin
          if(xpoint[0] ne 0 or xpoint[1] ne 0) then begin
              if(xpoint[1] lt axis[1]) then begin
                  lcfs_mask = lcfs_xy[1,*] gt xpoint[1]
              endif else begin
                  lcfs_mask = lcfs_xy[1,*] lt xpoint[1]
              endelse
          endif else begin
              lcfs_mask = fltarr(1, n_elements(lcfs_xy[1,*]))
              lcfs_mask[*] = 1
          endelse
      endif else begin
          lcfs_mask = fltarr(1, n_elements(lcfs_xy[1,*]))
          lcfs_mask[*] = 1
      endelse
      
      oplot, lcfs_xy[0,*], lcfs_xy[1,*]
      
      nlim = fix(total(lcfs_mask))
      rlim = fltarr(nlim)
      zlim = fltarr(nlim)
      j = 0 
      for i=0, n_elements(lcfs_mask)-1 do begin
          if(not lcfs_mask[0, i]) then continue
          rlim[j] = lcfs_xy[0,i]
          zlim[j] = lcfs_xy[1,i]
          j = j+1
      end
      print, 'lim points = ', nlim, j
  endelse

  ; reduce the boundary points if necessary
  while(nlim ge 500) do begin
      print, 'reducing lim points...'
      nlim = nlim / 2
      new_rlim = fltarr(nlim)
      new_zlim = fltarr(nlim)
      for k=0, nlim-1 do begin
          new_rlim[k] = rlim[2*k]
          new_zlim[k] = zlim[2*k]
      endfor
      rlim = new_rlim
      zlim = new_zlim
      print, 'new lim points = ', nlim
  end
      
  oplot, rlim, zlim, psym=4
  oplot, rwall, zwall, psym=6
 
  r2bp = psi_r^2 + psi_z^2

  ; calculate flux averages
  p = flux_average_field(p0,psi,x,z,flux=flux,nflux=nflux,_EXTRA=extra)
  pp = (p0_r*psi_r + p0_z*psi_z)/r2bp
  pprime = flux_average_field(pp,psi,x,z,nflux=nflux,flux=flux,_EXTRA=extra)
  I = flux_average_field(I0,psi,x,z,flux=flux,nflux=nflux,_EXTRA=extra)
  ffp = I0*(I0_r*psi_r + I0_z*psi_z)/r2bp
  ffprim = flux_average_field(ffp,psi,x,z,flux=flux,nflux=nflux,_EXTRA=extra)
  q = flux_average('q',slice=time,psi=psi,x=x,z=z,t=t,flux=flux,nflux=nflux,_EXTRA=extra)
;  q = smooth(q,5,/edge)
  jb = (I0_z*psi_r - I0_r*psi_z - jphi*I0)/r^2
  jdotb = flux_average_field(jb,psi,x,z,t,flux=flux,nflux=nflux,_EXTRA=extra)
  r2i = flux_average_field(1./r^2,psi,x,z,flux=flux,nflux=nflux,_EXTRA=extra)

  betacent = field_at_point(beta[0,*,*], x, z, axis[0], axis[1])

  ; to cgs............. to si
  c = 3e10
  p = p*b0^2/(4.*!pi)           / 10.
  p0 = p0*b0^2/(4.*!pi)         / 10.
  psi = psi*b0*l0^2             / 1e8
  flux = flux*b0*l0^2           / 1e8
  flux0 = flux0*b0*l0^2         / 1e8
  psilim=psilim*b0*l0^2         / 1e8
  pprime = pprime*b0^2/(4.*!pi) / 10.
  pprime = pprime/(b0*l0^2)     * 1e8
  I = I*b0*l0                   / (1e4*100.)
  ffprim = ffprim*(b0*l0)^2     / (1e4*100.)^2
  ffprim = ffprim/(b0*l0^2)     * 1e8
  bzero = bzero*b0              / 1e4
  b2 = b2*b0^2                  / (1e4)^2
;  tcur = tcur*b0*c*l0/(4.*!pi)  / 3e9
  r = r*l0                      / 100.
  x = x*l0                      / 100.
  z = z*l0                      / 100.
  axis[0] = axis[0]*l0          / 100.
  axis[1] = axis[1]*l0          / 100.
  rzero = rzero*l0              / 100.
  jdotb = jdotb*b0^2*c/(l0*4.*!pi) / (1e4*3e5)
  r2i = r2i/(l0^2)              * 100.^2

  nr = n_elements(flux)
  print, 'nr = ', nr

  psimin = flux0

  if(psimin gt psilim) then begin
      psi = -psi
      psilim = -psilim
      psimin = -psimin
      flux = -flux
      pprime = -pprime
      ffprim = -ffprim
  end

  name = ['M3DC1', '03/17/', '2013    ', $
          '#000000', '0000', '']
  idum = 3
  nr = n_elements(flux)
  nz = n_elements(z)
  rdim = max(x) - min(x)
  zdim = max(z) - min(z)
  xplas = 1.
  ccon = min(x)
  zmid = (max(z) + min(z))/2.
  rmag = axis[0]
  zmag = axis[1]
  zip = tcur[0]
  bcentr = bzero*rzero/rmag
  beta0 = beta0
  beta_n = 100.*(bzero*rzero/rmag)*beta0/(zip/1e6)
  xdum = 0.

  print, 'rdim, zdim', rdim, zdim
  print, 'rmag, zmag', rmag, zmag
  print, 'bcentr = ', bcentr
  print, 'beta0 = ', beta0
  print, 'betacent = ', betacent
  print, 'beta_n = ', beta_n
  print, 'zip = ', zip
  print, 'psimin, psilim = ', psimin, psilim
  print, 'min, max (flux) = ', min(flux), max(flux)


  ; set up formatting codes
  f2000 = '(6A8,3I4)'
  f2020 = '(5E16.9)'
  f2022 = '(2i5)' 


  ; output to eqdsk file
  print, 'outputting to eqdsk format...'
  file = 1
  
  openw, file, eqfile

  printf, file, format=f2000, name, idum, nr, nz
  printf, file, format=f2020, rdim, zdim, xplas, ccon, zmid
  printf, file, format=f2020, rmag, zmag, psimin, psilim, bcentr
  printf, file, format=f2020, zip, psimin, beta0, rmag, betacent
  printf, file, format=f2020, zmag, beta_n, psilim, xdum, xdum

  printf, file, format=f2020, I
  printf, file, format=f2020, p
  printf, file, format=f2020, ffprim
  printf, file, format=f2020, pprime
  printf, file, format=f2020, reform(psi[0,*,*])
  printf, file, format=f2020, q
  printf, file, format=f2022, nlim, nwall
  printf, file, format=f2020, transpose([[rlim],[zlim]])
  printf, file, format=f2020, transpose([[rwall],[zwall]])

  close, file

  ; output to jsolver
  print, 'outputting to jsolver format...'

  jsfile = 'jsfile'
  openw, file, jsfile  

  ncycle=  1
  isyms=  0
  ipest=  1
  kmax=nlim-1
  npsit = n_elements(flux)-1 ; remove final point to avoid NaN's
  
  times=  0.1140E-01
  xaxes=  axis[0]
  zmags=  axis[1]
  apls=   0.3656E+05
  betas=  0.5184E-02
  betaps= 0.1472E+00
  ali2s=  0.5458E+00
  qsaws=  0.5000E+00
  psimins=min(flux)
  psilims=max(flux)

  ; jsolver wants boundary points in clockwise direction
  if(ifixedb eq 1) then begin
      print, 'Reversing boundary points'
      rlim = reverse(rlim)
      zlim = reverse(zlim)
  end

  f2201 = '(20x,10a8)'
  f6100 = '(5i10)'
  f6101 = '(5e20.12)'

  ;  gzeros= R * B_T in m-T   (at vacuum)
  gzeros = bzero*rzero

  printf, file, format=f6100, ncycle,isyms,ipest,npsit,kmax
  printf, file, format=f6101, times,xaxes,zmags,gzeros,apls,betas,betaps, $
    ali2s,qsaws,psimins,psilims

  mu0 = (4.*!pi*1.e-7)
  ajpest2 = jdotb/(I*r2i)

  printf, file, format=f6101, mu0*p[0:npsit-1]
  printf, file, format=f6101, mu0*pprime[0:npsit-1]
  printf, file, format=f6101, -mu0*ajpest2[0:npsit-1]
  printf, file, format=f6101, flux[0:npsit-1] - flux[0]
  printf, file, format=f6101, rlim
  printf, file, format=f6101, zlim

  close, file

  window, 1
;  contour_and_legend, psi,x,z, /iso
;  loadct,12
;  oplot, rlim, zlim, color=color(1,3), thick=3.0
;  oplot, rwall, zwall, color=color(2,3), thick=3.0
 
;  plot, flux, pprime

  loadct,12
  !p.multi = [0,3,2]
  plot, nflux, mu0*p, title='p'
  plot, nflux, mu0*pprime, title="p'"
  oplot, nflux, mu0*deriv(flux,p), color=color(1,2), linestyle=2
  plot, nflux, I, title='f'
  plot, nflux, ffprim, title="f f'"
  oplot, nflux, I*deriv(flux,I), color=color(1,2), linestyle=2
  plot, nflux, q, title='q', yrange=[0,10]
  plot, nflux, mu0*ajpest2, title='ajpest2'

 
  !p.multi=0

  window, 0

end
