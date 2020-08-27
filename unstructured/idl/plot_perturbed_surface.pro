; plot_perturbed_surface
; ~~~~~~~~~~~~~~~~~~~~~~
; returns the perturbed surfaces at specified values of q
;
; INPUTS
;  ntheta: number of poloidal points per purface
;  tpoints: number of toroidal points per surface
;  points: number of sampling points in the R,Z plane
;  /sum: add together the perturbed fields from the different files.
;        Equilibrium fields are taken from the first file.
;  phase: the toroidal phase factor to apply to each perturbed field
;  scalefac: the amplitude factor to apply to each perturbed field
; 
; IN / OUT
;  theta: if given as input, this specifies the poloidal angles
;         on output contains the values of the poloidal angles
;  phi:   if given as input, this specifies the toroidal angles
;         on output contains the values of the toroidal angles
;
; OUTPUTS
;  

pro plot_perturbed_surface, psi_norm, scalefac=scalefac, points=pts,  $
                            filename=filename, phi=phi0, _EXTRA=extra, $
                            out=out, color=color, names=names, $
                            overplot=overplot, xy_out=xy_out, theta=theta, $
                            ntheta=ntheta, axis=axis_out, flux=flux, $
                            noplot=noplot, nonlinear=nonlinear, sum_files=sum,$
                            smooth=smooth, plotscale=plotscale, phase=phase, $
                            tpoints=nphi, slice=slice, fc=fc

   if(n_elements(filename) eq 0) then filename='C1.h5'
   n = n_elements(filename)
   nr = n_elements(psi_norm)
   if(n_elements(theta) eq 0) then begin
      if(n_elements(ntheta) eq 0) then ntheta = 400
      theta = 2.*!pi*findgen(ntheta)/(ntheta) - !pi
   endif else begin
      if(n_elements(ntheta) eq 0) then ntheta = n_elements(theta)
   end
   if(n_elements(filename) gt 1) then begin
      itor = read_parameter('itor', filename=filename[0])
      rzero = read_parameter('rzero', filename=filename[0])
   endif else begin
      itor = read_parameter('itor', filename=filename)
      rzero = read_parameter('rzero', filename=filename)
   end
   if(itor eq 1) then begin
      period = 2.*!pi
   endif else begin
      period = 2.*!pi*rzero
   end
   if(n_elements(phi0) eq 0) then phi0=0.
   if(n_elements(nphi) eq 0) then begin
      nphi = n_elements(phi0)
   endif else begin
      phi0 = findgen(nphi)/nphi * period
      if(itor eq 1) then phi0 = phi0*180./!pi
   end
   if(n_elements(scalefac) eq 0) then scalefac=1.
   if(n_elements(scalefac) lt n) then scalefac = replicate(scalefac, n)
   if(n_elements(plotscale) eq 0) then plotscale=1.
   if(n_elements(phase) eq 0) then phase = 0.
   if(n_elements(phase) lt n) then phase = replicate(phase, n)

   nn = n
   if(keyword_set(sum)) then nn=1
   xy_out = fltarr(nn,nr,2,nphi,ntheta)
   axis_out = fltarr(nn,2)
   flux = fltarr(nn,nr)

   if(n gt 1 and not keyword_set(sum)) then begin
       c = shift(get_colors(), -1)

       for i=0, n-1 do begin
           plot_perturbed_surface, psi_norm, scalefac=scalefac[i], points=pts,$
                                   filename=filename[i], phi=phi0, color=c[i],$
                                   overplot=(i ne 0), xy_out=xy_out_tmp, $
                                   axis=axis_tmp, flux=flux_tmp, $
                                   noplot=noplot, phase=phase[i], $
                                   plotscale=plotscale, $
                                   ntheta=ntheta, smooth=smooth, $
                                   slice=slice, fc=fc
           xy_out[i,*,*,*,*] = xy_out_tmp[0,*,*,*,*]
           axis_out[i,*] = axis_tmp[0,*]
           flux[i,*] = flux_tmp[0,*]
        end

       if(n_elements(names) eq 0) then names = filename
       if(not keyword_set(noplot)) then plot_legend, names, color=c
       return
   end

   psi0 = read_field('psi',x,z,t,filename=filename[0], slice=-1, $
                     points=pts, _EXTRA=extra)
   psi0_r = read_field('psi',x,z,t,filename=filename[0], slice=-1, $
                       points=pts, _EXTRA=extra, op=2)
   psi0_z = read_field('psi',x,z,t, filename=filename[0], slice=-1, $
                       points=pts, _EXTRA=extra, op=3)
   
   if(keyword_set(nonlinear)) then begin

      Te = fltarr(nphi,pts,pts)

      Te0 = read_field('Te',x,z,t,filename=filename[0], slice=-1, $
                       points=pts, _EXTRA=extra)
      Te0_fa = flux_average_field(Te0, psi0, x, z, t, filename=filename[0], $
                                  flux=te_flux, points=pts, fc=fc)
      for i=0, nphi-1 do Te[i,*,*] = Te0

      for i=0, n-1 do begin
         icomp = read_parameter('icomplex', filename=filename[i])
         phi1 = phi0 - phase[i]
         Te1 = read_field('Te',x,z,t,filename=filename[i], /linear, $
                          points=pts, slice=slice, complex=icomp, $
                          phi=phi1)

         Te = Te + scalefac[i]*real_part(Te1)
      end
     
   endif else begin
      icomp = read_parameter('icomplex', filename=filename[0])
      phi1 = phi0 - phase[0]
      zi = read_field('displacement',x,z,t,filename=filename[0], /linear, $
                      points=pts, slice=slice, complex=icomp, $
                      phi=phi1) * scalefac[0]
      for i=0, n-1 do begin
         icomp = read_parameter('icomplex', filename=filename[i])
         phi1 = phi0 - phase[i]
         zi = zi + $
              read_field('displacement',x,z,t,filename=filename[i], /linear, $
                         points=pts, slice=slice, complex=icomp, $
                         phi=phi1) * scalefac[i]
      end
      zi = real_part(zi)
;      if(n_elements(smooth) ne 0) then begin
;         for k=0, nphi-1 do zi[k,*,*] = gauss_smooth(reform(zi[k,*,*]),smooth)
;      end
   end

   if(n_elements(bins) eq 0) then bins = n_elements(x)

   psilim = lcfs(file=filename[0], flux0=psimin)
   fvals = psimin + psi_norm*(psilim - psimin)
   flux[0,*] = fvals

   xhat = psi0_r/sqrt(psi0_r^2 + psi0_z^2)
   zhat = psi0_z/sqrt(psi0_r^2 + psi0_z^2)

   if(not keyword_set(noplot)) then begin
      if(not keyword_set(overplot)) then begin
         plot, x, z, /nodata, /iso, _EXTRA=extra, $
               xtitle='!8R!6 (m)!X', ytitle='!8Z!6 (m)!X', $
               subtitle='plot scale = x'+string(format='(G0.0)',plotscale)
      end
      c = get_colors()
      c0 = c
      if(n_elements(fvals) eq 1) then begin
         l0 = 0
         c[1] = c0[3]
      endif else begin
         l0 = 1
      endelse
   endif

   psimin = read_lcfs(axis=axis, xpoint=xpoint, flux0=flux0, $
                    filename=filename[0], slice=slice)

   nk = 1
   if(n_elements(fvals) gt 8) then nk = ceil(n_elements(fvals) / 8)

   xy0 = fltarr(nr,2,nphi,ntheta)

   for k=0, n_elements(fvals)-1 do begin
      print, 'Surface ', k+1, ' of ', n_elements(fvals)
      for m=0, nphi-1 do begin
         xy = path_at_flux(psi0, x, z, t, fvals[k], axis=axis, /contiguous, $
                          /refine)
         if(keyword_set(nonlinear)) then begin
            teval = interpol(Te0_fa,te_flux,fvals[k])
            xy_new = path_at_flux(Te[m,*,*], x, z, t, teval, axis=axis, $
                                  /contiguous, /refine)
         endif else begin
            xy_new = xy
            
            for j=0, n_elements(xy[0,*])-1 do begin
               dx = field_at_point(xhat, x, z, xy[0,j], xy[1,j])
               dz = field_at_point(zhat, x, z, xy[0,j], xy[1,j])
               dr = field_at_point(zi[m,*,*], x, z, xy[0,j], xy[1,j])
               xy_new[0,j] = xy[0,j] + dr*dx
               xy_new[1,j] = xy[1,j] + dr*dz
            end
         end
         
                                ; interpolate onto regular theta grid
         theta0 = atan(xy[1,*]-axis[1],xy[0,*]-axis[0])
         i = sort(theta0)
         xy0[k,0,m,*] = interpol(xy[0,i],theta0[i],theta)
         xy0[k,1,m,*] = interpol(xy[1,i],theta0[i],theta)
         
         theta0 = atan(xy_new[1,*]-axis[1],xy_new[0,*]-axis[0])
         i = sort(theta0)
         t0 = theta0[i]
         xy_out[0,k,0,m,*] = interpol(xy_new[0,i],t0,theta)
         xy_out[0,k,1,m,*] = interpol(xy_new[1,i],t0,theta)
      
         d = reform(xy_out[0,k,*,m,*]) - reform(xy0[k,*,m,*])
         if(n_elements(smooth) gt 0) then begin
            d[0,*] = gauss_smooth(d[0,*],smooth,/edge_wrap)
            d[1,*] = gauss_smooth(d[1,*],smooth,/edge_wrap)

            xy_out[0,k,*,m,*] = xy0[k,*,m,*] + d[*,*]
         end

         if(not keyword_set(noplot) and (k mod nk) eq 0) then begin
            kk = k / nk
            if(n_elements(color) ne 0) then begin
               cc = color
               cc0 = c[0]
            endif else begin
               cc = c[kk+1]
               cc0 = c0[kk+1]
            endelse

            oplot, xy0[k,0,m,*], xy0[k,1,m,*], linestyle=l0, color=cc0
            oplot, xy0[k,0,m,*] + d[0,*]*plotscale, $
                   xy0[k,1,m,*] + d[1,*]*plotscale, color=cc
         end
      end
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
           printf, ifile, format='(2F15.5)', xy_out[0,k]
           free_lun, ifile
       end
   end
end
