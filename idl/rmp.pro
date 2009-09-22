pro plot_br, _EXTRA=extra, plotq=plotq, b=b, usepsi=usepsi

   psi0 = read_field('psi',x,z,t,slice=-1,_EXTRA=extra)
   psi1r = read_field('psi',x,z,t,/last,/linear,_EXTRA=extra)
   psi1i = read_field('psi_i',x,z,t,/last,/linear,_EXTRA=extra)

   bp = sqrt(s_bracket(psi1r,psi1r,x,z) + s_bracket(psi1i,psi1i,x,z)) $
     / radius_matrix(x,z,t)
   br_r = a_bracket(psi0,psi1r,x,z) $
     / (sqrt(s_bracket(psi0,psi0,x,z))*radius_matrix(x,z,t))
   br_i = a_bracket(psi0,psi1i,x,z) $
     / (sqrt(s_bracket(psi0,psi0,x,z))*radius_matrix(x,z,t))

   if(keyword_set(plotq)) then begin
       q = flux_average('q',psi=psi,x=x,z=z,t=t,flux=flux2,slice=-1, $
                        nflux=nflux2, _EXTRA=extra)
   end

   ; shift frequency values so that most negative frequency comes first
   if(keyword_set(usepsi)) then begin
       a_r = flux_coord_field(psi1r,psi0,x,z,t,flux=flux,angle=angle, $
                              area=area,nflux=nflux,_EXTRA=extra)
       a_i = flux_coord_field(psi1i,psi0,x,z,t,flux=flux,angle=angle, $
                              area=area, nflux=nflux,_EXTRA=extra)
   endif else begin
       a_r = flux_coord_field(br_r,psi0,x,z,t,flux=flux,angle=angle, $
                              area=area,nflux=nflux,_EXTRA=extra)
       a_i = flux_coord_field(br_i,psi0,x,z,t,flux=flux,angle=angle, $
                              area=area,nflux=nflux,_EXTRA=extra)
   endelse

   a = complex(a_r, a_i)
   b = transpose(a,[0,2,1])
   c = fft(b, -1, dimension=2)

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

   ; convert result to Tesla / radian
   get_normalizations, b0=b0, n0=n0, l0=l0, _EXTRA=extra
   if(b0 eq 0.) then b0 = 1.e4

   b0 = b0/1.e4

   print, 'b0 = ', b0

   contour_and_legend, abs(d)*b0, m, sqrt(nflux),  $
     table=39, xtitle=xtitle, ytitle=ytitle, $
     xrange=[-15,15], yrange=[0,1], _EXTRA=extra

   if(keyword_set(plotq)) then begin
       ntor = 3
       oplot, ntor*q, sqrt(nflux2), linestyle=2, color=!d.table_size-1
       oplot, -ntor*q, sqrt(nflux2), linestyle=2, color=!d.table_size-1
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
