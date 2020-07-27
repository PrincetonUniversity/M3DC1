pro write_neo_input, q, _EXTRA=extra, out=outfile, nphi=nphi, $
                     scalefac=scalefac, phase=phase, sum=sum, points=pts, $
                     plotscale=plotscale, filename=filename, $
                     nonlinear=nonlinear, slice=slice, smooth=smooth

  if(n_elements(filename) eq 0) then filename='C1.h5'
  if(n_elements(nonlinear) eq 0) then nonlinear=1
  if(n_elements(pts) eq 0) then pts=600
  if(n_elements(bins) eq 0) then bins = 50
  ntheta = 400
  if(n_elements(smooth) eq 0) then smooth = 3
  if(n_elements(nphi) eq 0) then nphi = 16
  if(n_elements(plotscale) eq 0) then plotscale=10.
  if(n_elements(q) eq 0) then begin
     q = flux_average('q', flux=qflux, psi=psi0, x=x,z=z,t=t,fc=fc,points=pts,$
                      bins=bins, i0=i0, slice=-1, filename=filename[0])
  end
  nr = n_elements(q)
  
  print, 'A'

  theta = 2.*!pi*findgen(ntheta)/(ntheta) - !pi
  phi = 360.*findgen(nphi)/(nphi)  ; toroidal angle in degrees
  r0 = fltarr(nr,nphi,ntheta)
  z0 = fltarr(nr,nphi,ntheta)

  plot_perturbed_surface, q, filename=filename, $
                          xy_out=xy_out, theta=theta, phi=phi, $
                          flux=flux, nonlinear=nonlinear,  $
                          scalefac=scalefac, phase=phase, smooth=smooth, $
                          /sum, points=pts, plotscale=plotscale, slice=slice

  r0[*,*,*] = xy_out[0,*,0,*,*]
  z0[*,*,*] = xy_out[0,*,1,*,*]

  print, 'B'

  ne0 = flux_average('ne',flux=flux,psi=psi0,x=x,z=z,t=t,fc=fc,points=pts,$
                     bins=bins, i0=i0, slice=-1,/mks,filename=filename[0])
  ni0 = flux_average('den',flux=flux,psi=psi0,x=x,z=z,t=t,fc=fc,points=pts,$
                     bins=bins, i0=i0, slice=-1,/mks,filename=filename[0])
  Te0 = flux_average('Te',flux=flux,psi=psi0,x=x,z=z,t=t,fc=fc,points=pts,$
                     bins=bins, i0=i0, slice=-1,/mks,filename=filename[0])
  Ti0 = flux_average('Ti',flux=flux,psi=psi0,x=x,z=z,t=t,fc=fc,points=pts,$
                     bins=bins, i0=i0, slice=-1,/mks,filename=filename[0])

  print, 'C'

  if(n_elements(qflux) eq 0) then begin
     qflux = flux_at_q(q, psi=psi0,x=x,z=z,t=t,fc=fc, /unique)
  end
  ;; ne0_x = interpol(ne0,qflux,flux)
  ;; ni0_x = interpol(ni0,qflux,flux)
  ;; Te0_x = interpol(Te0,qflux,flux)
  ;; Ti0_x = interpol(Ti0,qflux,flux)

  ;; plot, qflux, Te0
  ;; if(n_elements(flux) gt 1) then begin
  ;;    oplot, flux, Te0_x, psym=4
  ;; endif
;  stop

  ion_mass = read_parameter('ion_mass', filename=filename[0])
  psilim = read_lcfs(axis=axis, xpoint=xpoint, flux0=flux0,/mks, $
                     filename=filename[0])

  if(n_elements(outfile) eq 0) then outfile='m3dc1_neo.nc'
  id = ncdf_create(outfile, /clobber)
  ncdf_attput, id, 'version', 1, /short, /global
  ncdf_attput, id, 'psi_0', flux0, /global
  ncdf_attput, id, 'psi_1', psilim, /global
  ncdf_attput, id, 'ion_mass', ion_mass , /global

  npsi = n_elements(flux)
  npsi_id = ncdf_dimdef(id, 'npsi', npsi); number of psi points for profiles
  nr_id = ncdf_dimdef(id, 'nr', nr)      ; number of radial points for surfaces
  np_id = ncdf_dimdef(id, 'np', ntheta)  ; number of poloidal points
  nt_id = ncdf_dimdef(id, 'nt', nphi)    ; number of toroidal points

  q_var = ncdf_vardef(id, 'q', [nr_id], /float)
  psi_var = ncdf_vardef(id, 'psi', [nr_id], /float)
  phi_var = ncdf_vardef(id, 'Phi', [nt_id], /float)
  psi_var = ncdf_vardef(id, 'psi0', [npsi_id], /float)
  Te_var = ncdf_vardef(id, 'Te0', [npsi_id], /float)
  Ti_var = ncdf_vardef(id, 'Ti0', [npsi_id], /float)
  ni_var = ncdf_vardef(id, 'ni0', [npsi_id], /float)
  ne_var = ncdf_vardef(id, 'ne0', [npsi_id], /float)
  r_var = ncdf_vardef(id, 'R', [nr_id,nt_id,np_id], /float)
  z_var = ncdf_vardef(id, 'Z', [nr_id,nt_id,np_id], /float)
  ncdf_control, id, /endef
  ncdf_varput, id, 'q', q
  ncdf_varput, id, 'psi', reform(qflux)
  ncdf_varput, id, 'Phi', phi
  ncdf_varput, id, 'psi0', reform(flux)
  ncdf_varput, id, 'Te0', reform(Te0)
  ncdf_varput, id, 'Ti0', reform(Ti0)
  ncdf_varput, id, 'ne0', reform(ne0)
  ncdf_varput, id, 'ni0', reform(ni0)
  ncdf_varput, id, 'R', r0
  ncdf_varput, id, 'Z', z0
  ncdf_close, id

end
