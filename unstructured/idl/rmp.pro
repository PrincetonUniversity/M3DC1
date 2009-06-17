pro plot_br, _EXTRA=extra, plotq=plotq

   get_normalizations, b0=b0, n0=n0, l0=l0, _EXTRA=extra

   if(b0 eq 0.) then b0 = 1e4

   psi0 = read_field('psi',x,z,t,slice=-1,_EXTRA=extra)
   psi1r = read_field('psi',x,z,t,/last,/linear,_EXTRA=extra)
   psi1i = read_field('psi_i',x,z,t,/last,/linear,_EXTRA=extra)
  
   if(keyword_set(plotq)) then begin
       q = flux_average('q',slice,psi=psi,x=x,z=z,t=t,flux=flux2,slice=-1, $
                        nflux=nflux2, _EXTRA=extra)
   end

   a_r = flux_coord_field(psi1r,psi0,x,z,t,flux=flux,angle=angle,area=area, $
                        nflux=nflux,_EXTRA=extra)
   a_i = flux_coord_field(psi1i,psi0,x,z,t,flux=flux,angle=angle,area=area, $
                        nflux=nflux,_EXTRA=extra)

   a = complex(a_r, a_i)

   b = transpose(a,[0,2,1])

   c = fft(b, -1, dimension=2)

   n = n_elements(angle)
   d = c
   if((n mod 2) eq 0) then begin
       d[0,0:n/2-2,*] = c[0,n/2+1:n-1,*]
       d[0,n/2-1,*] = c[0,0,*]
       d[0,n/2:n-1,*] = c[0,1:n/2,*]
       f = indgen(n) - n/2 + 1
   endif else begin
   endelse

    for i=0, n-1 do begin
        for j=0, n_elements(flux)-1 do begin
            d[0,i,j] = d[0,i,j]*f[i]/area[j]
        end
    end

   xtitle='!6m!X'
   ytitle='!9r!7W!X'

   contour_and_legend, abs(d)*b0, f, sqrt(nflux),  $
     color_table=39, xtitle=xtitle, ytitle=ytitle, $
     xrange=[-15,15], yrange=[0,1]

   if(keyword_set(plotq)) then begin
       ntor = 3
       oplot, ntor*q, sqrt(nflux2), linestyle=2, color=!d.table_size-1
       oplot, -ntor*q, sqrt(nflux2), linestyle=2, color=!d.table_size-1
   endif
end
