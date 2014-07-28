pro schaffer_plot, field, x,z,t, q=q, _EXTRA=extra, bins=bins, q_val=q_val, $
                   psi_val=psi_val, ntor=ntor, label=label, psi0=psi0, i0=i0, $
                   m_val=m_val, phase=phase, overplot=overplot, $
                   linestyle=linestyle, outfile=outfile, bmnfile=bmnfile, $
                   bmncdf=bmncdf, rhs=rhs, reverse_q=reverse_q

   print, 'Drawing schaffer plot'

   if(n_elements(psi0) eq 0) then begin
;       psi0 = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra)
       psi0 = read_field('psi',x,z,t,slice=-1,_EXTRA=extra)
   endif
   if(n_elements(i0) eq 0) then begin
;       i0   = read_field('i'  ,x,z,t,/equilibrium,_EXTRA=extra)
       i0   = read_field('i'  ,x,z,t,slice=-1,_EXTRA=extra)
   endif
   if(n_elements(ntor) eq 0) then begin
       ntor = read_parameter('ntor',_EXTRA=extra)
   endif

   r = radius_matrix(x,z,t)

   ; From Schaffer 2008
   ; Br_mn = [(2*pi)^2/S] * [J Br]_mn
   ; J = B_theta*q*R^3/(R_0*B_0)
   bp = sqrt(s_bracket(psi0,psi0,x,z))/r
   jac = r^3*bp/abs(i0)

   if(size(field, /type) eq 7) then begin
       field = read_field(field,x,z,t,/complex,_EXTRA=extra)
   endif

   a_r = flux_coord_field(real_part(field)*jac,psi0,x,z,t, $
                          flux=flux,angle=angle,qval=q, $
                          area=area,nflux=nflux,tbins=bins,fbins=bins, $
                          /pest,i0=i0, _EXTRA=extra)
   a_i = flux_coord_field(imaginary(field)*jac,psi0,x,z,t, $
                          flux=flux,angle=angle, qval=q, $
                          area=area,nflux=nflux,tbins=bins,fbins=bins, $
                          /pest,i0=i0, _EXTRA=extra)

   if(keyword_set(reverse_q)) then q = -q

   for i=0, n_elements(angle)-1 do begin
       a_r[0,*,i] = (2.*!pi)^2*a_r[0,*,i]*q/area
       a_i[0,*,i] = (2.*!pi)^2*a_i[0,*,i]*q/area
   end 

   a = complex(a_r, a_i)
   b = transpose(a,[0,2,1])
   c = fft(b, -1, dimension=2)

   ; shift frequency values so that most negative frequency comes first
   n = n_elements(angle)
   f = indgen(n)
   f[n/2+1] = n/2 + 1 - n + findgen((n-1)/2)
   m = shift(f,-(n/2+1))
   d = shift(c,0,-(n/2+1),0)

   if(n_elements(bmnfile) ne 0) then begin
       openw, ifile, /get_lun, bmnfile

       printf, ifile, format='(3I5)', ntor, n_elements(nflux), n_elements(m)
       printf, ifile, format='(1F13.6)', nflux
       printf, ifile, format='(1I5)', transpose(m)
       for i=0, n_elements(nflux)-1 do begin
           printf, ifile, format='(2F13.6)', $
             real_part(reform(d[0,*,i])), $
             imaginary(reform(d[0,*,i]))
       end

       free_lun, ifile
   end

   if(n_elements(bmncdf) ne 0) then begin
        omega_i = flux_average('v_omega',flux=qflux,psi=psi0,x=x,z=z,t=t,$
                               bins=bins, i0=i0, slice=-1,/mks,_EXTRA=extra)
        omega_e = flux_average('ve_omega',flux=qflux,psi=psi0,x=x,z=z,t=t,$
                               bins=bins, i0=i0, slice=-1,/mks,_EXTRA=extra)
        F = flux_average('I',flux=qflux,psi=psi0,x=x,z=z,t=t,$
                               bins=bins, i0=i0, slice=-1,/mks,_EXTRA=extra)
        p = flux_average('p',flux=qflux,psi=psi0,x=x,z=z,t=t,$
                               bins=bins, i0=i0, slice=-1,/mks,_EXTRA=extra)
        pe = flux_average('pe',flux=qflux,psi=psi0,x=x,z=z,t=t,$
                               bins=bins, i0=i0, slice=-1,/mks,_EXTRA=extra)
        den_e = flux_average('ne',flux=qflux,psi=psi0,x=x,z=z,t=t,$
                             bins=bins, i0=i0, slice=-1,/mks,_EXTRA=extra)

       rpath = fltarr(n_elements(m), n_elements(nflux))
       zpath = fltarr(n_elements(m), n_elements(nflux))
       read_nulls, axis=axis, _EXTRA=extra

       for i=0, n_elements(nflux)-1 do begin
           xy = path_at_flux(psi0,x,z,t,flux[i],/contiguous,$
                             path_points=n_elements(m))
           rpath[i,*] = xy[0,*]
           zpath[i,*] = xy[1,*]

                                ; do 3 Newton iterations
           for j=0, 3 do begin
               theta_geom = reform(atan(zpath[i,*]-axis[1],rpath[i,*]-axis[0]))
               dum = min(theta_geom, ind, /abs)
               if(dum eq 0 and ind eq 0) then break
               im = ind-1 mod n_elements(m)
               ip = ind+1 mod n_elements(m)
               if((dum lt 0 and theta_geom(ip) gt 0) or $
                  (dum gt 0 and theta_geom(ip) lt 0)) then begin
                   dtheta = theta_geom(ip) - dum
               endif else if((dum lt 0 and theta_geom(im) gt 0) or $
                             (dum gt 0 and theta_geom(im) lt 0)) then begin
                   dtheta = dum - theta_geom(im)
               endif else begin
                   print, 'Error!', theta_geom(im), dum, theta_geom(ip)
               endelse
               dind = -dum/dtheta
               index = findgen(n_elements(m)) + ind + dind
               wg = where(index ge n_elements(m), c)
               if(c gt 0) then index[wg] = index[wg] - n_elements(m)
               wg = where(index lt 0, c)
               if(c gt 0) then index[wg] = index[wg] + n_elements(m)
               rpath[i,*] = interpolate(rpath[i,*], index)
               zpath[i,*] = interpolate(zpath[i,*], index)
           end
       end
       
       plot, [1.0, 2.3], [-1.5, 1.5], /nodata
       for i=0, n_elements(nflux) - 1, 10 do begin
           oplot, rpath[*, i], zpath[*,i]
       end
;       stop

       bpval = fltarr(n_elements(m), n_elements(nflux))
       bpval = field_at_point(bp, x, z, rpath, zpath)
       
       id = ncdf_create(bmncdf, /clobber)
       ncdf_attput, id, 'ntor', ntor, /short, /global
       n_id = ncdf_dimdef(id, 'npsi', n_elements(nflux))
       m_id = ncdf_dimdef(id, 'mpol', n_elements(m))
       psi_var = ncdf_vardef(id, 'psi', [n_id], /float)
       flux_pol_var = ncdf_vardef(id, 'flux_pol', [n_id], /float)
       q_var = ncdf_vardef(id, 'q', [n_id], /float)
       p_var = ncdf_vardef(id, 'p', [n_id], /float)
       F_var = ncdf_vardef(id, 'F', [n_id], /float)
       pe_var = ncdf_vardef(id, 'pe', [n_id], /float)
       ne_var = ncdf_vardef(id, 'ne', [n_id], /float)
       omega_i_var = ncdf_vardef(id, 'omega_i', [n_id], /float)
       omega_e_var = ncdf_vardef(id, 'omega_e', [n_id], /float)
       m_var = ncdf_vardef(id, 'm', [m_id], /short)
       bmn_real_var = ncdf_vardef(id, 'bmn_real', [m_id,n_id], /float)
       bmn_imag_var = ncdf_vardef(id, 'bmn_imag', [m_id,n_id], /float)
       rpath_var = ncdf_vardef(id, 'rpath', [m_id,n_id], /float)
       zpath_var = ncdf_vardef(id, 'zpath', [m_id,n_id], /float)
       bp_var = ncdf_vardef(id, 'Bp', [m_id,n_id], /float)
       ncdf_control, id, /endef
       ncdf_varput, id, 'psi', reform(nflux[0,*])
       ncdf_varput, id, 'flux_pol', reform(flux[0,*])
       ncdf_varput, id, 'm', m
       ncdf_varput, id, 'q', abs(reform(q))
       ncdf_varput, id, 'p', reform(p)
       ncdf_varput, id, 'F', reform(F)
       ncdf_varput, id, 'pe', reform(pe)
       ncdf_varput, id, 'ne', reform(den_e)
       ncdf_varput, id, 'omega_e', reform(omega_e)
       ncdf_varput, id, 'omega_i', reform(omega_i)
       ncdf_varput, id, 'bmn_real', real_part(reform(d[0,*,*]))
       ncdf_varput, id, 'bmn_imag', imaginary(reform(d[0,*,*]))
       ncdf_varput, id, 'rpath', rpath
       ncdf_varput, id, 'zpath', zpath
       ncdf_varput, id, 'Bp', bpval
       ncdf_close, id
   end

   
   if(n_elements(m_val) ne 0) then begin

       q_val = abs(m_val/float(ntor))
       indices = interpol(findgen(n_elements(q)), abs(q), q_val)

       if(n_elements(linestyle) eq 0) then linestyle=0

       for i=0, n_elements(m_val)-1 do begin
           dum = min(m-m_val[i], j, /abs)

           c = get_colors()

           if(keyword_set(phase)) then begin
               data = atan(imaginary(d[0,j,*]), real_part(d[0,j,*]))
               ytitle='!6Phase!X'
           endif else begin
               data = abs(d[0,j,*])
               ytitle=label
           endelse
           if(i eq 0 and not keyword_set(overplot)) then begin
               plot, nflux, data, color=c[0], $
                 /nodata, _EXTRA=extra, xtitle='!7W!X', ytitle=ytitle
           endif

           oplot, nflux, data, color=c[i+1], linestyle=linestyle

           fv = interpolate(nflux,indices[i])
           oplot, [fv,fv], !y.crange, color=c[i+1], linestyle=1
       end

       if(n_elements(m_val) gt 1) then begin
           names = string(format='("!8m!6 = ",I0,"!X")', m_val)
           plot_legend, names, color=c[1:n_elements(m_val)], _EXTRA=extra
       endif

       return
   endif

   if (n_elements(psi_val) ne 0) then begin
       indices = interpol(findgen(n_elements(nflux)), nflux, psi_val)
   endif else if(n_elements(q_val) ne 0) then begin
       indices = interpol(findgen(n_elements(q)), abs(q), q_val)
   endif

   if(n_elements(indices) ne 0) then begin
       print, abs(q[fix(indices)])
       print, abs(q[fix(indices+1)])

       b = complexarr(n_elements(angle), n_elements(indices))
       for i=0, n_elements(angle)-1 do begin
           b[i,*] = interpolate(reform(a[0,*,i]), indices)
       end
;       b = interpolate(reform(a[0,*,*]), indices)
       
       c = fft(b, -1, dimension=1)
       dold = d
       if(n_elements(indices) eq 1) then begin
           d = shift(c,-(n/2+1))
       endif else begin
           d = shift(c,-(n/2+1),0)
       endelse

       col = fltarr(n_elements(indices))
       c = get_colors()
       for i=0, n_elements(col)-1 do begin
           col[i] = c[i mod n_elements(c)]
       end
       if(n_elements(outfile) ne 0) then begin
           openw, ifile, outfile, /get_lun
       end

       for i=0, n_elements(indices)-1 do begin
           dum = min(m-q_val[i]*ntor, j, /abs)
           dum = min(m+q_val[i]*ntor, k, /abs)

           print, 'q, Psi = ', interpolate(abs(q),indices[i]), $
             interpolate(nflux, indices[i])
           print, 'Resonant field: m (mag, phase) = ', m[j], abs(d[j,i]), $
             atan(imaginary(d[j,i]),real_part(d[j,i]))
           print, 'Resonant field: m (mag, phase) = ', m[k], abs(d[k,i]), $
             atan(imaginary(d[k,i]),real_part(d[k,i]))

           if(i eq 0) then begin
               plot, m, abs(d[*,i]), xrange=[-20,20], yrange=[0, max(abs(d))]
           endif else begin
               oplot, m, abs(d[*,i]), color=col[i]
           end
           
           if(n_elements(outfile) ne 0) then begin
               print, 'q[0] = ', q[0]
               if(q[0] gt 0) then mi = j else mi = k
               print, 'm[mi] = ', m[mi]
               printf, ifile, format='(I5,7F12.6)', m[mi], abs(d[mi,i]), $
                 atan(imaginary(d[mi,i]),real_part(d[mi,i])), $
                 interpolate(nflux, indices[i]), $
                 interpolate(abs(q), indices[i]), $
                 interpolate(deriv(nflux, abs(q)), indices[i]), $
                 interpolate(area, indices[i]), $
                 interpolate(deriv(nflux, flux), indices[i])
           end
       end
       if(n_elements(outfile) ne 0) then begin
           free_lun, ifile
       end
       return
   endif

   if(1 eq strcmp(!d.name, 'PS', /fold_case)) then begin
       xsize = 1.33
   endif else begin
       xsize = 1.
   endelse

   xtitle='!8m!X'
   ytitle='!9r!7W!X'

   y = sqrt(nflux)
;   y = nflux

   contour_and_legend, abs(d), m, y,  $
     table=39, xtitle=xtitle, ytitle=ytitle, $
     xrange=[-20,20], yrange=[0,1], /lines, c_thick=1, $
     ccolor=!d.table_size-1, label=label, $
     _EXTRA=extra, xsize=xsize

   oplot, ntor*q, y, linestyle=2, color=!d.table_size-1
end
