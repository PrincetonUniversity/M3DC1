pro fit_profiles, _EXTRA=extra, error=error

  n = 200
  grid_psin = findgen(n)/(n-1.)*(1.1)
  grid_p = fltarr(n)
  grid_te = fltarr(n)
  grid_ne = fltarr(n)
  grid_pe = fltarr(n)
  grid_w = fltarr(n)
  chi_p = 0.
  chi_te = 0.
  chi_pe = 0.
  chi_ne = 0.
  chi_w = 0.

  itmax = 100
  pmin = 1000.     ; minimum pressure, in Pa
  pemin = 100.

  ; fit profile
  !p.multi=[0,3,2]

  if(~file_test('geqdsk',/read)) then begin
     print, 'Error: geqdsk does not exist'
     return
  endif else begin
     efit = readg('geqdsk')
     efit_p = efit.pres
     efit_psin = findgen(n_elements(efit_p))/n_elements(efit_p)
     fit_modtanh, efit_psin, efit_p, grid_psin, grid_p, param=a_p, $
                  chisq=chi_p, /plot, itmax=itmax, minval=pmin
  end

  if(~file_test('profile_te',/read)) then begin
     print, 'Error: profile_te does not exist'
     return
  endif else begin
     te_file = read_ascii('profile_te')
     te = reform(te_file.field1[1,*])
     te_psin = reform(te_file.field1[0,*])
     fitrange = where(te_psin lt 1.)
     fit_modtanh, te_psin[fitrange], te[fitrange], grid_psin, grid_te, param=a_te, minval=1e-3, $
                  chisq=chi_te, /plot, itmax=itmax
  end

  if(~file_test('profile_ne',/read)) then begin
     print, 'Error: profile_ne does not exist'
     return
  endif else begin
     dene_file = read_ascii('profile_ne')
     dene = reform(dene_file.field1[1,*])
     dene_psin = reform(dene_file.field1[0,*])
     fitrange = where(te_psin lt 1.)
     fit_modtanh, dene_psin[fitrange], dene[fitrange], grid_psin, grid_ne, param=a_ne, $
                  chisq=chi_ne, /plot, itmax=itmax

;     pe = interpol(dene, dene_psin, te_psin)*te
;     ; convert to Pa
;     pe = pe*1000.*1.60218e-19*1e20
;     pe_psin = te_psin
;     fit_modtanh, pe_psin, pe, grid_psin, grid_pe, param=a_pe, minval=pemin, $
;                  chisq=chi_pe, /plot, itmax=itmax
;     grid_ne = grid_pe/(1000.*1.60218e-19*1e20) / grid_te
                                ; enforce monotonicity on ne outside psi_norm = 0.98
;     for i=0, n_elements(grid_ne)-1 do begin
;        if(grid_psin[i] lt 0.98) then continue
;        if(grid_ne[i] gt grid_ne[i-1]) then grid_ne[i] = grid_ne[i-1]
;     end

;     plot, dene_psin, dene
;     oplot, grid_psin, grid_ne, color=color(1)
  end
  

;  print, 'OMEGA_ExB'
  if(~file_test('profile_omega',/read)) then begin
     print, 'Error: profile_omega does not exist'
     return
  endif else begin
     w_file = read_ascii('profile_omega')
     w = reform(w_file.field1[1,*])
     w_psin = reform(w_file.field1[0,*])
;     dpsi = efit.ssibry - efit.ssimag 
;  a_w = a_pe[0:5]
;  a_w[2:5] = a_w[2:5]/dpsi
;  print, 'a_w = ', a_w
     
     fit_modtanh, w_psin, w, grid_psin, grid_w, param=a_w, /deriv, $
                  chisq=chi_w, /plot, itmax=itmax

  endelse

  !p.multi=0


; Do Checks
  error = 0
  if(min(grid_p) lt 0.) then begin
     print, 'ERROR: pressure < 0'
     error = 1
  end
;  if(min(grid_pe) lt 0.) then begin
;     error = 2
;     print, 'ERROR: electron pressure < 0'
;  end
  if(min(grid_te) lt 0.) then begin
     error = 3
     print, 'ERROR: electron temperature < 0'
  end
  if(min(grid_p - grid_pe) lt 0.) then begin
     error = 4
     print, 'ERROR: ion pressure < 0'
  end
  if(min(grid_ne) lt 0.) then begin
     error = 5
     print, 'ERROR: density < 0'
  end

  
  openw, ifile, "fit_profiles.log", /get_lun
  printf, ifile, error
  printf, ifile, format='(7E12.4)', chi_p, a_p
  printf, ifile, format='(7E12.4)', chi_te, a_te
  printf, ifile, format='(8E12.4)', chi_ne, a_ne
  printf, ifile, format='(7E12.4)', chi_w, a_w
  free_lun, ifile

  openw, ifile, "profile_p.fit", /get_lun
  for i=0, n-1 do begin
     printf, ifile, grid_psin[i], grid_p[i]
  end
  free_lun, ifile

  openw, ifile, "profile_te.fit", /get_lun
  for i=0, n-1 do begin
     printf, ifile, grid_psin[i], grid_te[i]
  end
  free_lun, ifile

  openw, ifile, "profile_ne.fit", /get_lun
  for i=0, n-1 do begin
     printf, ifile, grid_psin[i], grid_ne[i]
  end
  free_lun, ifile

  openw, ifile, "profile_omega.fit", /get_lun
  for i=0, n-1 do begin
     printf, ifile, grid_psin[i], grid_w[i]
  end
  free_lun, ifile

end
