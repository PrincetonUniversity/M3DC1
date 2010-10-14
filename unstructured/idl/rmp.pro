pro schaffer_plot, field, x,z,t, q=q, _EXTRA=extra, bins=bins, q_val=q_val, $
                   ntor=ntor, label=label, psi0=psi0, i0=i0

   print, 'Drawing schaffer plot'

   if(n_elements(psi0) eq 0) then begin
       psi0 = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra)
   endif
   if(n_elements(i0) eq 0) then begin
       i0   = read_field('i'  ,x,z,t,/equilibrium,_EXTRA=extra)
   endif

   r = radius_matrix(x,z,t)

   bt = sqrt(s_bracket(psi0,psi0,x,z))/r
   jac = r^3*bt/abs(i0)

   a_r = flux_coord_field(real_part(field)*jac,psi0,x,z,t, $
                          flux=flux,angle=angle,q=q, $
                          area=area,nflux=nflux,tbins=bins,fbins=bins, $
                          /pest, _EXTRA=extra)
   a_i = flux_coord_field(imaginary(field)*jac,psi0,x,z,t, $
                          flux=flux,angle=angle, q=q,$
                          area=area,nflux=nflux,tbins=bins,fbins=bins, $
                          /pest, _EXTRA=extra)

    for i=0, n_elements(angle)-1 do begin
        a_r[0,*,i] = a_r[0,*,i]*q/area
        a_i[0,*,i] = a_i[0,*,i]*q/area
    end

   a = complex(a_r, a_i)
   b = transpose(a,[0,2,1])
   c = sqrt(2.*!pi)*fft(b, -1, dimension=2)

   ; shift frequency values so that most negative frequency comes first
   n = n_elements(angle)
   f = indgen(n)
   f[n/2+1] = n/2 + 1 - n + findgen((n-1)/2)
   m = shift(f,-(n/2+1))
   d = shift(c,0,-(n/2+1),0)
   
   xtitle='!8m!X'
   ytitle='!9r!7W!X'

   if(n_elements(q_val) ne 0) then begin
       dum = min(q-q_val, i, /abs)
       print, 'Psi, q = ', nflux[i], q[0,i]
       plot, m[*], abs(d[0,*,i])
       return
   endif

   if(1 eq strcmp(!d.name, 'PS', /fold_case)) then begin
       xsize = 1.33
   endif else begin
       xsize = 1.
   endelse

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

   schaffer_plot, br, x, z, t, _EXTRA=extra, psi0=psi0,i0=i0, $
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
