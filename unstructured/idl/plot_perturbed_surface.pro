pro plot_perturbed_surface, q, scalefac=scalefac, points=pts, $
                            filename=filename, phi=phi0, _EXTRA=extra, $
                            out=out, color=color, names=names, $
                            overplot=overplot, xy_out=xy_out, theta=theta, $
                            ntheta=ntheta, axis=axis_out, flux=flux

   n = n_elements(filename)
   nr = n_elements(q)
   if(n_elements(theta) eq 0) then begin
      if(n_elements(ntheta) eq 0) then ntheta = 100
      theta = 2.*!pi*findgen(ntheta)/(ntheta) - !pi
   endif else begin
      if(n_elements(ntheta) eq 0) then ntheta = n_elements(theta)
   end

   xy_out = fltarr(n,nr,2,ntheta)
   axis_out = fltarr(n,2)
   flux = fltarr(n,nr)

   if(n gt 1) then begin
       c = shift(get_colors(), -1)
       for i=0, n-1 do begin
           plot_perturbed_surface, q, scalefac=scalefac, poinst=pts, $
             filename=filename[i], phi=phi0, color=c[i], overplot=(i ne 0), $
             _EXTRA=extra, xy_out=xy_out_tmp, axis=axis_tmp, flux=flux_tmp
           xy_out[i,*,*,*] = xy_out_tmp[0,*,*,*]
           axis_out[i,*] = axis_tmp[0,*]
           flux[i,*] = flux_tmp[0,*]
       end

       if(n_elements(names) eq 0) then names = filename
       plot_legend, names, color=c
       return
   end

   icomp = read_parameter('icomplex', filename=filename)
   if(n_elements(scalefac) eq 0) then scalefac=1.
   if(n_elements(scalefac) eq 1) then $
     scalefac=replicate(scalefac, n_elements(q))

   psi0 = read_field('psi',x,z,t,filename=filename, slice=-1, $
                     points=pts, _EXTRA=extra)
   psi0_r = read_field('psi',x,z,t,filename=filename, slice=-1, $
                       points=pts, _EXTRA=extra, op=2)
   psi0_z = read_field('psi',x,z,t, filename=filename, slice=-1, $
                       points=pts, _EXTRA=extra, op=3)
   zi = read_field('displacement',x,z,t,filename=filename, /linear, $
                   points=pts, slice=slice, complex=icomp, $
                   phi=phi0, _EXTRA=extra)
   
   if(n_elements(bins) eq 0) then bins = n_elements(x)
   fvals = flux_at_q(q, points=pts, filename=filename, $
                     normalized_flux=norm, bins=bins)
   print, fvals

   xhat = psi0_r/sqrt(psi0_r^2 + psi0_z^2)
   zhat = psi0_z/sqrt(psi0_r^2 + psi0_z^2)

   if(not keyword_set(overplot)) then begin
       plot, x, z, /nodata, /iso, _EXTRA=extra, $
         xtitle='!8R!6 (m)!X', ytitle='!8Z!6 (m)!X'
   end
   c = get_colors()
   c0 = c
   if(n_elements(fvals) eq 1) then begin
       l0 = 0
       c[1] = c0[3]
   endif else begin
       l0 = 1
   endelse

   psimin = read_lcfs(axis=axis, xpoint=xpoint, flux0=flux0, $
                    filename=filename, slice=slice)

   for k=0, n_elements(fvals)-1 do begin
       if(n_elements(color) ne 0) then begin
           cc = color
           cc0 = c[0]
       endif else begin
           cc = c[k+1]
           cc0 = c0[k+1]
       endelse
       xy = path_at_flux(psi0, x, z, t, fvals[k], axis=axis, /contiguous)
       xy_new = xy

       for j=0, n_elements(xy[0,*])-1 do begin
           dx = field_at_point(xhat, x, z, xy[0,j], xy[1,j])
           dz = field_at_point(zhat, x, z, xy[0,j], xy[1,j])
           dr = field_at_point(zi, x, z, xy[0,j], xy[1,j])
           xy_new[0,j] = xy[0,j] + dr*dx*scalefac[k]
           xy_new[1,j] = xy[1,j] + dr*dz*scalefac[k]
       end

       ; interpolate onto regular theta grid
       theta0 = atan(xy_new[1,*]-axis[1],xy_new[0,*]-axis[0])
       i = sort(theta0)
       xy_out[0,k,0,*] = interpol(xy_new[0,i],theta0[i],theta)
       xy_out[0,k,1,*] = interpol(xy_new[1,i],theta0[i],theta)

       oplot, xy[0,*], xy[1,*], linestyle=l0, color=cc0
       oplot, [xy[0,n_elements(xy[0,*])-1], xy[0,0]],  $
         [xy[1,n_elements(xy[0,*])-1], xy[1,0]], linestyle=l0, color=cc0

       oplot, xy_new[0,*], xy_new[1,*], color=cc
       oplot, [xy_new[0,n_elements(xy[0,*])-1], xy_new[0,0]], $
         [xy_new[1,n_elements(xy[0,*])-1], xy_new[1,0]], color=cc
   end

   if(n_elements(out) ne 0) then begin
       for k=0, n_elements(fvals)-1 do begin
           ; write equilibrium
           out0 = string(format='(A,"_q=",F0.3,"_scale=0.txt")',out,q[k])

           openw, ifile, out0, /get_lun
           printf, ifile, format='(2F15.5)', xy
           free_lun, ifile

           ; write perturbed surface
           out1 = string(format='(A,"_q=",F0.3,"_scale=",F0.3,".txt")',out, $
                         q[k], scalefac)
           openw, ifile, out1, /get_lun
           printf, ifile, format='(2F15.5)', xy_new
           free_lun, ifile
       end
   end
end
