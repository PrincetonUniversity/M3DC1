pro schaffer_plot, field, x,z,t, q=q, _EXTRA=extra, bins=bins, q_val=q_val, $
                   psi_val=psi_val, ntor=ntor, label=label, psi0=psi0, i0=i0, $
                   m_val=m_val

   print, 'Drawing schaffer plot'

   if(n_elements(psi0) eq 0) then begin
       psi0 = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra)
   endif
   if(n_elements(i0) eq 0) then begin
       i0   = read_field('i'  ,x,z,t,/equilibrium,_EXTRA=extra)
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
                          flux=flux,angle=angle,q=q, $
                          area=area,nflux=nflux,tbins=bins,fbins=bins, $
                          /pest, _EXTRA=extra)
   a_i = flux_coord_field(imaginary(field)*jac,psi0,x,z,t, $
                          flux=flux,angle=angle, q=q,$
                          area=area,nflux=nflux,tbins=bins,fbins=bins, $
                          /pest, _EXTRA=extra)

;   plot, nflux, q, _EXTRA=extra
;   return

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
   
   if(n_elements(m_val) ne 0) then begin

       q_val = abs(m_val/ntor)
       indices = interpol(findgen(n_elements(q)), q, q_val)       

       for i=0, n_elements(m_val)-1 do begin
           dum = min(m-m_val[i], j, /abs)

           c = colors()
           if(i eq 0) then begin
               plot, nflux, abs(d[0,j,*]), color=c[0], /nodata, _EXTRA=extra
           endif

           oplot, nflux, abs(d[0,j,*]), color=c[i+1]

           fv = interpolate(nflux,indices[i])
           oplot, [fv,fv], !y.crange, color=c[i+1], linestyle=1
       end

       if(n_elements(m_val) gt 1) then begin
           plot_legend, m_val, color=c[1:n_elements(m_val)], _EXTRA=extra
       endif

       return
   endif

   if (n_elements(psi_val) ne 0) then begin
       indices = interpol(findgen(n_elements(nflux)), nflux, psi_val)
   endif else if(n_elements(q_val) ne 0) then begin
       indices = interpol(findgen(n_elements(q)), q, q_val)       
   endif

   if(n_elements(indices) ne 0) then begin
       print, q[fix(indices)]
       print, q[fix(indices+1)]

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

       col = colors()
       for i=0, n_elements(indices)-1 do begin
           dum = min(m-q_val[i]*ntor, j, /abs)
           dum = min(m+q_val[i]*ntor, k, /abs)

           print, 'q, Psi = ', interpolate(q,indices[i]), $
             interpolate(nflux, indices[i])
           print, 'Resonant field: m = ', m[j], abs(d[j,i]), $
             atan(imaginary(d[j,i]),real_part(d[j,i]))
           print, 'Resonant field: m = ', m[k], abs(d[k,i]), $
             atan(imaginary(d[k,i]),real_part(d[k,i]))

           if(i eq 0) then begin
               plot, m, abs(d[*,i]), xrange=[-20,20], yrange=[0, max(abs(d))]
           endif else begin
               oplot, m, abs(d[*,i]), color=col[i]
           end
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

   contour_and_legend, abs(d), m, sqrt(nflux),  $
     table=39, xtitle=xtitle, ytitle=ytitle, $
     xrange=[-20,20], yrange=[0,1], /lines, c_thick=1, $
     ccolor=!d.table_size-1, label=label, $
     _EXTRA=extra, xsize=xsize

   oplot, -ntor*q, sqrt(nflux), linestyle=2, color=!d.table_size-1
end

pro plot_br, _EXTRA=extra, bins=bins, q_val=q_val, $
             subtract_vacuum=subtract_vacuum, ntor=ntor, $
             plotbn=plotbn, slice=slice

   if(n_elements(ntor) eq 0) then begin
       ntor = read_parameter('ntor', _EXTRA=extra)
   endif
   print, 'ntor = ', ntor
   if(n_elements(slice) eq 0) then last=1

   psi0 = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra)
   i0   = read_field('i'  ,x,z,t,/equilibrium,_EXTRA=extra)

   bx = read_field('bx',x,z,t,last=last,slice=slice, $
                   /linear,_EXTRA=extra,/complex)
   by = read_field('by',x,z,t,last=last,slice=slice, $
                   /linear,_EXTRA=extra,/complex)

   if(keyword_set(subtract_vacuum)) then begin
       bx0 = read_field('bx',x,z,t,slice=0,/linear,_EXTRA=extra)
       by0 = read_field('by',x,z,t,slice=0,/linear,_EXTRA=extra)
       bx = bx - bx0
       by = by - by0
   endif   

   r = radius_matrix(x,z,t)
   y = z_matrix(x,z,t)

   bt = sqrt(s_bracket(psi0,psi0,x,z))/r

   br = -(bx*s_bracket(psi0,r,x,z) + by*s_bracket(psi0,y,x,z)) / (r*bt)

   ; convert to cgs
   get_normalizations, b0=b0_norm, n0=n0_norm, l0=l0_norm, _EXTRA=extra
   br = br*b0_norm

   schaffer_plot, br, x, z, t, _EXTRA=extra, psi0=psi0,i0=i0, q_val=q_val, $
     label='!8B!Dn!N!6 (!8G!6)!X', points=points, bins=bins, ntor=ntor
end

pro plot_jpar, _EXTRA=extra, bins=bins, q_val=q_val, $
             ntor=ntor, slice=slice

   if(n_elements(ntor) eq 0) then begin
       ntor = read_parameter('ntor', _EXTRA=extra)
   endif
   print, 'ntor = ', ntor
   if(n_elements(slice) eq 0) then last=1

   psi0 = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra)
   i0 = read_field('i',x,z,t,/equilibrium,_EXTRA=extra)

   psi0_r = read_field('psi',x,z,t,last=last,/equilibrium, $
                   /linear,_EXTRA=extra,/complex,op=2)
   psi0_z = read_field('psi',x,z,t,last=last,/equilibrium, $
                   /linear,_EXTRA=extra,/complex,op=3)
   psi_r = read_field('psi',x,z,t,last=last,slice=slice, $
                       /linear,_EXTRA=extra,/complex,op=2)
   psi_z = read_field('psi',x,z,t,last=last,slice=slice, $
                       /linear,_EXTRA=extra,/complex,op=3)
   psi_lp = read_field('psi',x,z,t,last=last,slice=slice, $
                       /linear,_EXTRA=extra,/complex,op=7)
   i_r = read_field('i',x,z,t,last=last,slice=slice, $
                   /linear,_EXTRA=extra,/complex,op=2)
   i_z = read_field('i',x,z,t,last=last,slice=slice, $
                   /linear,_EXTRA=extra,/complex,op=3)
   f_r = read_field('f',x,z,t,last=last,slice=slice, $
                   /linear,_EXTRA=extra,/complex,op=2)
   f_z = read_field('f',x,z,t,last=last,slice=slice, $
                   /linear,_EXTRA=extra,/complex,op=3)

   r = radius_matrix(x,z,t)

   ddpsi = complex(0.,1.)*ntor

   jB = -i0*(psi_lp - psi_r/r)/r^2 $
     + ((i_r+ddpsi^2*f_r)*psi0_r + (i_z+ddpsi^2*f_z)*psi0_z)/r^2 $
     + (ddpsi*psi_z*psi0_r - ddpsi*psi_r*psi0_z)/r^3

   B2 = (abs(psi0_r)^2 + abs(psi0_z)^2)/r^2 + abs(i0)^2/r^2

   jpar = jB/sqrt(B2)

   ; convert to cgs
   get_normalizations, b0=b0_norm, n0=n0_norm, l0=l0_norm, _EXTRA=extra
   jpar = jpar*b0_norm*(3.e10)/(4.*!pi*l0_norm)
   ; convert to mks
   jpar = jpar/(3.e5)

   schaffer_plot, jpar, x, z, t, _EXTRA=extra, psi0=psi0,i0=i0, $
     label='!8J!d!9#!N!6 (A/m!U2!N!N)!X', points=points, bins=bins, ntor=ntor
end


pro plot_vpar, _EXTRA=extra, bins=bins, q_val=q_val, $
             ntor=ntor, slice=slice

   if(n_elements(ntor) eq 0) then begin
       ntor = read_parameter('ntor', _EXTRA=extra)
   endif
   print, 'ntor = ', ntor
   if(n_elements(slice) eq 0) then last=1

   psi0 = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra)
   i0 = read_field('i',x,z,t,/equilibrium,_EXTRA=extra)

   psi0_r = read_field('psi',x,z,t,last=last,/equilibrium, $
                   /linear,_EXTRA=extra,/complex,op=2)
   psi0_z = read_field('psi',x,z,t,last=last,/equilibrium, $
                   /linear,_EXTRA=extra,/complex,op=3)
   u_r = read_field('phi',x,z,t,last=last,slice=slice, $
                       /linear,_EXTRA=extra,/complex,op=2)
   u_z = read_field('phi',x,z,t,last=last,slice=slice, $
                       /linear,_EXTRA=extra,/complex,op=3)
   chi_r = read_field('chi',x,z,t,last=last,slice=slice, $
                   /linear,_EXTRA=extra,/complex,op=2)
   chi_z = read_field('chi',x,z,t,last=last,slice=slice, $
                   /linear,_EXTRA=extra,/complex,op=3)
   omega = read_field('omega',x,z,t,last=last,slice=slice, $
                   /linear,_EXTRA=extra,/complex,op=2)

   r = radius_matrix(x,z,t)

   ddpsi = complex(0.,1.)*ntor

   vB = omega*i0 $
     + (u_r*psi0_r + u_z*psi0_z) $
     + (chi_z*psi0_r - chi_r*psi0_z)/r^3

   B2 = (abs(psi0_r)^2 + abs(psi0_z)^2)/r^2 + abs(i0)^2/r^2

   vpar = vB/sqrt(B2)

   ; convert to cgs
   get_normalizations, b0=b0_norm, n0=n0_norm, l0=l0_norm, _EXTRA=extra
   vpar = vpar*b0_norm/sqrt(4.*!pi*n0_norm*1.6726e-24)
   ; convert to mks
   vpar = vpar/100.

   schaffer_plot, vpar, x, z, t, _EXTRA=extra, psi0=psi0,i0=i0, $
     label='!8v!d!9#!N!6 (m/s!N)!X', points=points, bins=bins, ntor=ntor
end


pro plot_epar, _EXTRA=extra, bins=bins, q_val=q_val, $
             ntor=ntor, slice=slice

   if(n_elements(ntor) eq 0) then begin
       ntor = read_parameter('ntor', _EXTRA=extra)
   endif
   print, 'ntor = ', ntor
   if(n_elements(slice) eq 0) then last=1

   psi0_lp = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra,op=7)
   psi0_r = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra,op=2)
   psi0_z = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra,op=3)
   i0 = read_field('i',x,z,t,/equilibrium,_EXTRA=extra)
   i0_r = read_field('i',x,z,t,/equilibrium,_EXTRA=extra,op=2)
   i0_z = read_field('i',x,z,t,/equilibrium,_EXTRA=extra,op=3)
   eta = read_field('eta',x,z,t,/equilibrium,_EXTRA=extra)

   psi1_lp = read_field('psi',x,z,t,last=last, $
                     /linear,_EXTRA=extra,/complex,op=7)
   psi1_r = read_field('psi',x,z,t,last=last, $
                     /linear,_EXTRA=extra,/complex,op=2)
   psi1_z = read_field('psi',x,z,t,last=last, $
                     /linear,_EXTRA=extra,/complex,op=3)
   i1 = read_field('i',x,z,t,last=last,slice=slice, $
                   /linear,_EXTRA=extra,/complex)
   f1 = read_field('f',x,z,t,last=last,slice=slice, $
                   /linear,_EXTRA=extra,/complex)

   r = radius_matrix(x,z,t)

   jphi0 = psi0_lp - psi0_r/r
   jphi1 = psi1_lp - psi1_r/r


   ddpsi = complex(0.,1.)*ntor

   epar = (eta/r^2)* $
     (s_bracket(i1 + ddpsi*ddpsi*f1,psi0,x,z) + (i0_r*psi1_r + i0_z*psi1_z) $
      - i0*jphi1 - i1*jphi0 $
      - ddpsi*(psi0_z*psi1_r - psi0_r*psi1_z)/r $
      + r*a_bracket(i0,ddpsi*f1,x,z))

   B2 = (abs(psi0_r)^2 + abs(psi0_z)^2)/r^2 + abs(i0)^2/r^2

   epar = epar/sqrt(B2)

;   contour_and_legend, abs(epar), x, z, _EXTRA=extra
;   return

   ; convert to cgs
   get_normalizations, b0=b0_norm, n0=n0_norm, l0=l0_norm, _EXTRA=extra
   epar = epar*b0_norm^2/sqrt(4.*!pi*n0_norm*1.6726e-24)/3.e10
   ; convert to mks
   epar = epar*3.e4

   schaffer_plot, epar, x, z, t, _EXTRA=extra, psi0=psi0,i0=i0, $
     label='!8E!D!9#!N!6 (V/m!N)!X', points=points, bins=bins, ntor=ntor
end



pro integrate_j, df=df, _EXTRA=extra
   
   ntor = read_parameter('ntor', _EXTRA=extra)
   q0 = findgen(8)/float(ntor)

   fq = flux_at_q(q0, _EXTRA=extra)

   psi0  = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra)
   psix_r = read_field('psi',x,z,t,/last,/linear,op=2,_EXTRA=extra)
   psix_i = read_field('psi_i',x,z,/last,/linear,op=2,_EXTRA=extra)
   psilp_r = read_field('psi',x,z,t,/last,/linear,op=7,_EXTRA=extra)
   psilp_i = read_field('psi_i',x,z,t,/last,/linear,op=7,_EXTRA=extra)
;   jphi = read_field('jphi',x,z,t,_EXTRA=extra)

   r = radius_matrix(x,z,t)

   jz = -(complex(psilp_r, psilp_i) - complex(psix_r,psix_i)/r)/r

;   contour_and_legend, abs(jz), x, z
;   return

   psibound = lcfs(flux0=flux0, _EXTRA=extra)

   ndf = 500

   int_j = fltarr(n_elements(fq), ndf)
   df = (findgen(ndf)+1.)/3000.

   for i=0, n_elements(fq)-1 do begin   
       for k=0, ndf-1 do begin
           deltaf = abs((flux0 - psibound)) * df[k]
           
           fmin = fq[i] - deltaf
           fmax = fq[i] + deltaf
           j = ((psi0 ge fmin) and (psi0 le fmax)) * abs(jz)
           int_j[i, k] = total(j*r)*mean(deriv(x))*mean(deriv(z))
       end
       
       if(i eq 0) then begin
           plot, df, deriv(int_j[i,*]), color=color(i,n_elements(fq)), yrange=[0,.1]
       endif else begin
           oplot, df, deriv(int_j[i,*]), color=color(i,n_elements(fq))
       endelse
   end

   plot_legend, string(format='("m = ",I)', fix(q0*ntor)), $
     color=colors(n_elements(fq))
end
