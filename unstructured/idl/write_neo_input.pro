pro write_neo_input, q, _EXTRA=extra, out=outfile

  if(n_elements(q) eq 0) then begin
     q = flux_average('q', flux=qflux, psi=psi0, x=x,z=z,t=t,fc=fc, $
                      bins=bin, i0=i0, slice=-1, _EXTRA=extra)
  end
  nr = n_elements(q)
  ntheta = 100
  nphi = 16
  
  theta = 2.*!pi*findgen(ntheta)/(ntheta) - !pi
  phi = 2.*!pi*findgen(nphi)/(nphi)
  r = fltarr(nr,nphi,ntheta)
  z = fltarr(nr,nphi,ntheta)

  for i=0, nphi-1 do begin
     plot_perturbed_surface, q, $
                             xy_out=xy_out, theta=theta, phi=phi[i]*180./!pi, $
                             flux=flux, scale=10, _EXTRA=extra, overplot=(i gt 0), noplot=(nr ge 10)
     r[*,i,*] = xy_out[0,*,0,*]
     z[*,i,*] = xy_out[0,*,1,*]
  end


  ne0 = flux_average('ne',flux=qflux,psi=psi0,x=x,z=z,t=t,fc=fc,$
                     bins=bins, i0=i0, slice=-1,/mks,_EXTRA=extra)
  ni0 = flux_average('den',flux=qflux,psi=psi0,x=x,z=z,t=t,fc=fc,$
                     bins=bins, i0=i0, slice=-1,/mks,_EXTRA=extra)
  Te0 = flux_average('Te',flux=qflux,psi=psi0,x=x,z=z,t=t,fc=fc,$
                     bins=bins, i0=i0, slice=-1,/mks,_EXTRA=extra)
  Ti0 = flux_average('Ti',flux=qflux,psi=psi0,x=x,z=z,t=t,fc=fc,$
                     bins=bins, i0=i0, slice=-1,/mks,_EXTRA=extra)
  
  ne0_x = interpol(ne0,qflux,flux)
  ni0_x = interpol(ni0,qflux,flux)
  Te0_x = interpol(Te0,qflux,flux)
  Ti0_x = interpol(Ti0,qflux,flux)

  plot, qflux, Te0
  if(n_elements(flux) gt 1) then begin
     oplot, flux, Te0_x, psym=4
  endif
;  stop

  ion_mass = read_parameter('ion_mass', _EXTRA=extra)

  if(n_elements(outfile) eq 0) then outfile='m3dc1_neo.nc'
  id = ncdf_create(outfile, /clobber)
  ncdf_attput, id, 'version', 1, /short, /global
  ncdf_attput, id, 'ion_mass', ion_mass , /global

  nr_id = ncdf_dimdef(id, 'nr', nr)      ; number of radial points
  np_id = ncdf_dimdef(id, 'np', ntheta)  ; number of poloidal points
  nt_id = ncdf_dimdef(id, 'nt', nphi)    ; number of toroidal points

  q_var = ncdf_vardef(id, 'q', [nr_id], /float)
  phi_var = ncdf_vardef(id, 'Phi', [nt_id], /float)
  psi_var = ncdf_vardef(id, 'psi0', [nr_id], /float)
  Te_var = ncdf_vardef(id, 'Te0', [nr_id], /float)
  Ti_var = ncdf_vardef(id, 'Ti0', [nr_id], /float)
  ni_var = ncdf_vardef(id, 'ni0', [nr_id], /float)
  ne_var = ncdf_vardef(id, 'ne0', [nr_id], /float)
  r_var = ncdf_vardef(id, 'R', [nr_id,nt_id,np_id], /float)
  z_var = ncdf_vardef(id, 'Z', [nr_id,nt_id,np_id], /float)
  ncdf_control, id, /endef
  ncdf_varput, id, 'q', q
  ncdf_varput, id, 'Phi', phi*180./!pi
  ncdf_varput, id, 'psi0', reform(flux)
  ncdf_varput, id, 'Te0', reform(Te0_x)
  ncdf_varput, id, 'Ti0', reform(Ti0_x)
  ncdf_varput, id, 'ne0', reform(ne0_x)
  ncdf_varput, id, 'ni0', reform(ni0_x)
  ncdf_varput, id, 'R', r
  ncdf_varput, id, 'Z', z
  ncdf_close, id

end
