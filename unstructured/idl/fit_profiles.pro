pro fit_profiles, outfile=out, _EXTRA=extra

  n = 200
  grid_psin = findgen(n)/(n-1.)*(1.1)
  grid_p = fltarr(n)
  grid_te = fltarr(n)
  grid_ne = fltarr(n)
  grid_pe = fltarr(n)
  grid_omega = fltarr(n)

  ; fit profile
  !p.multi=[0,3,2]

  print, 'PRESSURE'
  efit = readg('geqdsk')
  efit_p = efit.pres
  efit_psin = findgen(n_elements(efit_p))/n_elements(efit_p)
  fit_modtanh, efit_psin, efit_p, grid_psin, grid_p, param=a_p, $
               chisq=chi_p, /plot
  print, 'Fit center: ', a_p[0]
  print, 'Fit width: ', 1./a_p[1]
  print, 'Fit height: ', a_p[2]

  print, 'ELECTRON TEMPERATURE'
  te_file = read_ascii('profile_te')
  te = reform(te_file.field1[1,*])
  te_psin = reform(te_file.field1[0,*])
  fit_modtanh, te_psin, te, grid_psin, grid_te, param=a_te, minval=1e-3, $
               chisq=chi_te, /plot
  print, 'Fit center: ', a_te[0]
  print, 'Fit width: ', 1./a_te[1]
  print, 'Fit height: ', a_te[2]

  print, 'ELECTRON PRESSURE'
  dene_file = read_ascii('profile_ne')
  dene = reform(dene_file.field1[1,*])
  dene_psin = reform(dene_file.field1[0,*])
  pe = interpol(dene, dene_psin, te_psin)*te
  pe_psin = te_psin
  fit_modtanh, pe_psin, pe, grid_psin, grid_pe, param=a_pe, $
               chisq=chi_pe, /plot
  print, 'Fit center: ', a_pe[0]
  print, 'Fit width: ', 1./a_pe[1]
  print, 'Fit height: ', a_pe[2]

  grid_pi = grid_p - grid_pe
  if(min(grid_pi) lt 0.) then begin
     print, 'ERROR: ion pressure < 0'
  endif else begin
     print, 'Ion pressure is positive everywhere.'
  end

  print, 'ELECTRON DENSITY' 
  grid_ne = grid_pe / grid_te
  ; enforce monotonicity on ne outside psi_norm = 0.98
  for i=0, n_elements(grid_ne)-1 do begin
     if(grid_psin[i] lt 0.98) then continue
     if(grid_ne[i] gt grid_ne[i-1]) then grid_ne[i] = grid_ne[i-1]
  end
  plot, dene_psin, dene
  oplot, grid_psin, grid_ne, color=color(1)
  

  print, 'OMEGA_ExB'
  w_file = read_ascii('profile_omega')
  w = reform(w_file.field1[1,*])
  w_psin = reform(w_file.field1[0,*])
  dpsi = efit.ssibry - efit.ssimag 
;  a_w = a_pe[0:5]
;  a_w[2:5] = a_w[2:5]/dpsi
;  print, 'a_w = ', a_w

  fit_modtanh, w_psin, w, grid_psin, grid_w, param=a_w, /deriv, $
               chisq=chi_w, /plot
;  mod_dtanhfit, grid_psin, a_w, grid_w
;  plot, w_psin, w
;  oplot, grid_psin, grid_w, color=color(1)


  !p.multi=0

  if(keyword_set(out)) then begin
     help, a_p, a_te, a_pe, a_w
     help, chi_p, chi_te, chi_pe, chi_w
     openw, ifile, out, /get_lun
     printf, ifile, format='(8F10.4)', chi_p, a_p
     printf, ifile, format='(7F10.4)', chi_te, a_te
     printf, ifile, format='(8F10.4)', chi_pe, a_pe
     printf, ifile, format='(7F10.4)', chi_w, a_w
     free_lun, ifile
  end

end
