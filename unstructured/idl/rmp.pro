pro plot_br, _EXTRA=extra, plotq=plotq, bins=bins, usepsi=usepsi, $
             subtract_vacuum=subtract_vacuum, vacuum=vacuum, ntor=ntor, $
             plotbn=plotbn

   if(n_elements(ntor) eq 0) then begin
       ntor = read_parameter('ntor', _EXTRA=extra)
   endif
   print, 'ntor = ', ntor

   ; convert result to Tesla / radian
   get_normalizations, b0=b0, n0=n0, l0=l0, _EXTRA=extra
   if(b0 eq 0.) then b0 = 1.e4
   if(l0 eq 0.) then l0 = 100.
   print, 'normalizations: b0, l0 = ', b0, l0

   psi0 = read_field('psi',x,z,t,slice=-1,_EXTRA=extra)

   If(keyword_set(vacuum)) then begin
       vacuum_field,psi1r,x,z,ntor
;       vacuum_field,psi1r2,x,z,-ntor
;       psi1r = (psi1r + psi1r2)/2.
       psi1i = psi1r*0.
   endif else begin
       bx = read_field('bx',x,z,t,/last,/linear,_EXTRA=extra, /complex)
       by = read_field('by',x,z,t,/last,/linear,_EXTRA=extra, /complex)
       i0 = read_field('i',x,z,t,slice=-1,_EXTRA=extra)
   endelse

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

   br_r = real_part(br)
   br_i = imaginary(br)

   if(keyword_set(usepsi)) then begin
       a_r = flux_coord_field(psi1r,psi0,x,z,t,flux=flux,angle=angle,q=q, $
                              area=area,nflux=nflux,tbins=bins,fbins=bins, $
                              _EXTRA=extra)
       a_i = flux_coord_field(psi1i,psi0,x,z,t,flux=flux,angle=angle,q=q, $
                              area=area, nflux=nflux,tbins=bins,fbins=bins, $
                              _EXTRA=extra)
   endif else begin
       if(keyword_set(plotbn)) then begin
           jac = 1.
           pest = 0
       endif else begin
           jac = r^3*bt/i0
           pest = 1
       endelse
       a_r = flux_coord_field(br_r*jac,psi0,x,z,t,flux=flux,angle=angle,q=q, $
                              area=area,nflux=nflux,tbins=bins,fbins=bins, $
                              pest=pest, _EXTRA=extra)
       a_i = flux_coord_field(br_i*jac,psi0,x,z,t,flux=flux,angle=angle, q=q,$
                              area=area,nflux=nflux,tbins=bins,fbins=bins, $
                              pest=pest, _EXTRA=extra)


       if(keyword_set(plotbn)) then begin
           dum = min(sqrt(nflux)-.97817, i, /abs)
           print, 'Psi, q = ', nflux[i], q[0,i]
           plot, angle, a_r[0,i,*], _EXTRA=extra
           return
       endif

       ; ignore flux surfaces where q > 10
       q = q*(q lt 10.)

       for i=0, n_elements(angle)-1 do begin
           a_r[0,*,i] = a_r[0,*,i]*q/area
           a_i[0,*,i] = a_i[0,*,i]*q/area
       end
   endelse
   a = complex(a_r, a_i)
   b = transpose(a,[0,2,1])
   c = sqrt(2.*!pi)*fft(b, -1, dimension=2)

   ; shift frequency values so that most negative frequency comes first
   n = n_elements(angle)
   f = indgen(n)
   f[n/2+1] = n/2 + 1 - n + findgen((n-1)/2)
   m = shift(f,-(n/2+1))
   d = shift(c,0,-(n/2+1),0)

   if(keyword_set(usepsi)) then begin
       ; B_r_{m,n} = psi_{m,n} * (2*pi)^2*m/S
       for i=0, n-1 do begin
           for j=0, n_elements(flux)-1 do begin
               d[0,i,j] = d[0,i,j]*(2.*!pi)^2*m[i]/area[j]
           end
       end
   endif
   
   xtitle='!6m!X'
   ytitle='!9r!7W!X'

   contour_and_legend, abs(d), m, sqrt(nflux),  $
     table=39, xtitle=xtitle, ytitle=ytitle, $
     xrange=[-15,15], yrange=[0,1], /lines, $
     ccolor=!d.table_size-1, label='!8B!Dn!N!6 (!8T!6)!X', $
     _EXTRA=extra

   if(keyword_set(plotq)) then begin
       ntor = 3
       oplot, ntor*q, sqrt(nflux), linestyle=2, color=!d.table_size-1
       oplot, -ntor*q, sqrt(nflux), linestyle=2, color=!d.table_size-1
   endif
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
