pro force_balance_x, _EXTRA=extra
   itor = read_parameter('itor', _EXTRA=extra)
   ntor = read_parameter('ntor', _EXTRA=extra)
   icomplex = read_parameter('icomplex', _EXTRA=extra)
   linear = read_parameter('linear', _EXTRA=extra)
   print, itor, ntor, icomplex, linear

   icomplex = 1

   psi0_r = read_field('psi',x,z,t,/equilibrium,complex=icomplex,_EXTRA=extra,$
                      op=2)
   psi0_z = read_field('psi',x,z,t,/equilibrium,complex=icomplex,_EXTRA=extra,$
                      op=3)
   g0 = read_field('i',x,z,t,/equilibrium,complex=icomplex,_EXTRA=extra)
   g0_r = read_field('i',x,z,t,/equilibrium,complex=icomplex,_EXTRA=extra,op=2)
   jphi0 = read_field('psi',x,z,t,/equilibrium,complex=icomplex,_EXTRA=extra,$
                      op=7)

   psi1_r = read_field('psi',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=2)
   psi1_z = read_field('psi',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=3)
   g1 = read_field('i',x,z,t,/linear,complex=icomplex,_EXTRA=extra)
   g1_r = read_field('i',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=2)
   if(icomplex eq 1) then begin
       f1_r = read_field('f',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=2)
       f1_z = read_field('f',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=3)
   endif
   p1_r = read_field('p',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=2)
   jphi1 = read_field('psi',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=7)


   ddphi = complex(0,1)*ntor

   if(itor eq 0) then begin
       rzero = read_parameter('rzero',_EXTRA=extra)
       ddphi = ddphi/rzero
       r = 1.
   endif else begin
       r = radius_matrix(x,z,t)
       jphi0 = jphi0 - psi0_r/r
       jphi1 = jphi1 - psi1_r/r
   endelse
   
   term1 = jphi0*ddphi*f1_z/r $
     -g0*ddphi*psi1_z/r^3 $
     -g0*(g1_r + ddphi*ddphi*f1_r)/r^2 $
     -g0_r*g1/r^2 $
     -(jphi0*psi1_r + jphi1*psi0_r)/r^2
   term2 = -p1_r
   
   term0 = term1+term2

   names = ['total', 'JxB', 'p']
   contour_and_legend, real_part([term0, term1, term2]), $
     x,z, title=names, range=[-0.01, 0.01]

   print, total(abs(term1)), total(abs(term2))
   print, total(abs(term0))/(total(abs(term1)) + total(abs(term2)))
end

pro force_balance_phi, _EXTRA=extra
   itor = read_parameter('itor', _EXTRA=extra)
   ntor = read_parameter('ntor', _EXTRA=extra)
   icomplex = read_parameter('icomplex', _EXTRA=extra)
   linear = read_parameter('linear', _EXTRA=extra)
   print, itor, ntor, icomplex, linear

   icomplex = 1

   psi0_r = read_field('psi',x,z,t,/equilibrium,complex=icomplex,_EXTRA=extra,$
                      op=2)
   psi0_z = read_field('psi',x,z,t,/equilibrium,complex=icomplex,_EXTRA=extra,$
                      op=3)
   g0_r = read_field('i',x,z,t,/equilibrium,complex=icomplex,_EXTRA=extra,op=2)
   g0_z = read_field('i',x,z,t,/equilibrium,complex=icomplex,_EXTRA=extra,op=3)


   mu = read_field('visc',x,z,t,/equilibrium,complex=icomplex, $
                     _EXTRA=extra)
   mu_c = read_field('visc_c',x,z,t,/equilibrium,complex=icomplex, $
                     _EXTRA=extra)
   w = read_field('omega',x,z,t,/linear,complex=icomplex,_EXTRA=extra)
   w_lp = read_field('omega',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=7)
   w_r = read_field('omega',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=2)
   
   psi1_r = read_field('psi',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=2)
   psi1_z = read_field('psi',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=3)
   g1_r = read_field('i',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=2)
   g1_z = read_field('i',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=3)
   if(icomplex eq 1) then begin
       f1_r = read_field('f',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=2)
       f1_z = read_field('f',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=3)
   endif
   p1 = read_field('p',x,z,t,/linear,complex=icomplex,_EXTRA=extra)

   ddphi = complex(0,1)*ntor

   if(itor eq 0) then begin
       rzero = read_parameter('rzero',_EXTRA=extra)
       ddphi = ddphi/rzero
       r = 1.
   endif else begin
       r = radius_matrix(x,z,t)
   endelse
   
   term1 = (g0_z*psi1_r - g0_r*psi1_z)/r $
     + ((g1_z+ddphi^2*f1_z)*psi0_r - (g1_r+ddphi^2*f1_r)*psi0_z)/r $
     - ddphi*(g0_z*f1_z + g0_r*f1_r) $
     - ddphi*(psi1_r*psi0_r + psi1_z*psi0_z)/r^2
   term2 = -ddphi*p1

   term3 = 2.*mu_c*ddphi^2*w + r^2*mu*(w_lp+3*w_r/r)
   
   term0 = term1+term2+term3

   names = ['total', 'JxB', 'p', 'visc']
   contour_and_legend, real_part([term0, term1, term2, term3]), $
     x,z, title=names, range=[-0.001, 0.001]

   print, total(abs(term1)), total(abs(term2))
   print, total(abs(term0))/(total(abs(term1)) + total(abs(term2)))
end


pro grad_shafranov, names=names, terms=terms, _EXTRA=extra
   psi_r  = read_field('psi', x, z, t, _EXTRA=extra, op=2)
   psi_z  = read_field('psi', x, z, t, _EXTRA=extra, op=3)
   psi_lp = read_field('psi', x, z, t, _EXTRA=extra, op=7)
   p_r = read_field('p', x, z, t, _EXTRA=extra, op=2)
   p_z = read_field('p', x, z, t, _EXTRA=extra, op=3)
   i   = read_field('i', x, z, t, _EXTRA=extra)
   i_r = read_field('i', x, z, t, _EXTRA=extra, op=2)
   i_z = read_field('i', x, z, t, _EXTRA=extra, op=3)
   den = read_field('den', x, z, t, _EXTRA=extra)
   omega = read_field('omega', x, z, t, _EXTRA=extra)
   r = radius_matrix(x,z,t)

   names = ['!7D!6*!7w!X', "!8R!U!62!N p'!X", "!8FF'!X", '!6v!9.G!6v!X']

   terms = fltarr(n_elements(names), n_elements(x), n_elements(z))

   terms[0,*,*] = psi_lp - psi_r/r
   terms[1,*,*] = r^2*(p_r*psi_r + p_z*psi_z) / (psi_r^2 + psi_z^2)
   terms[2,*,*] = i*(i_r*psi_r + i_z*psi_z) / (psi_r^2 + psi_z^2)
   terms[3,*,*] = -den*r^3*omega^2*psi_r / (psi_r^2 + psi_z^2)
end

pro force_balance_parallel, names=names, terms=terms, _EXTRA=extra
   psi_r  = read_field('psi', x, z, t, _EXTRA=extra, op=2)
   psi_z  = read_field('psi', x, z, t, _EXTRA=extra, op=3)
   p_r = read_field('p', x, z, t, _EXTRA=extra, op=2)
   p_z = read_field('p', x, z, t, _EXTRA=extra, op=3)
   den = read_field('den', x, z, t, _EXTRA=extra)
   omega = read_field('omega', x, z, t, _EXTRA=extra)
   r = radius_matrix(x,z,t)

   names = ["!6B!9.G!8p'!X", "!6B!9.!6v!9.G!6v!X"]

   terms = fltarr(n_elements(names), n_elements(x), n_elements(z))

   terms[0,*,*] = reform((p_z*psi_r - p_r*psi_z)/r)
   terms[1,*,*] = reform(psi_z*den*omega^2)
end


pro plot_equation, eqname, xrange=xrange, yrange=yrange,$
                   cosine=cosine, sine=sine, $
                   surface_average=surface_average, _EXTRA=extra
  call_procedure, eqname, names=names, terms=terms, _EXTRA=extra

  if(keyword_set(surface_average)) then begin
      psi0 = read_field('psi',x,z,t,_EXTRA=extra,/equilibrium)
      if(keyword_set(cosine)) then begin
          co = cos(read_field('theta',x,z,t,_EXTRA=extra,/equilibrium))
      endif else begin
          co = 1.
      end
      if(keyword_set(sine)) then begin
          sn = sin(read_field('theta',x,z,t,_EXTRA=extra,/equilibrium))
      endif else begin
          sn = 1.
      end

      c = shift(get_colors(),-1)
      if(n_elements(xrange) eq 0) then xrange = [0,1]
      if(n_elements(yrange) eq 0) then $
        yrange = [min(terms*co*sn), max(terms*co*sn)]
      plot, xrange, yrange, /nodata
      for i=0, n_elements(names)-1 do begin
          fa = flux_average_field(terms[i,*,*]*co*sn, psi0, x, z, t, $
                                  bins=bins, nflux=nflux, _EXTRA=extra)
          if(i eq 0) then begin
              tot = fa
          endif else begin
              tot = tot + fa
          end
          oplot, nflux, fa, color=c[i]
      end
      oplot, nflux, tot, linestyle=1

      plot_legend, names, color=c, _EXTRA=extra
  endif else begin
      contour_and_legend, terms, title=names
  endelse
end


pro pressure, kar=kar, _EXTRA=extra
   itor = read_parameter('itor', _EXTRA=extra)
   ntor = read_parameter('ntor', _EXTRA=extra)
   icomplex = read_parameter('icomplex', _EXTRA=extra)
   linear = read_parameter('linear', _EXTRA=extra)
   print, itor, ntor, icomplex, linear

   if(n_elements(kar) eq 0) then kar = 0.
   kap = read_parameter('kappat', _EXTRA=extra)
   print, 'kappa = ', kap
   gamma = 5./3.

   eta = read_field('eta',x,z,t,/equilibrium,_EXTRA=extra)
   psi0_lp = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra,op=7)
   psi0_r = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra,op=2)
   g0_r = read_field('i',x,z,t,/equilibrium,_EXTRA=extra,op=2)
   g0_z = read_field('i',x,z,t,/equilibrium,_EXTRA=extra,op=3)
   psi1_r = read_field('psi',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=2)
   psi1_z = read_field('psi',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=3)
   psi1_lp = read_field('psi',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=7)
   g1_r = read_field('i',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=2)
   g1_z = read_field('i',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=3)
   if(icomplex eq 1) then begin
       f1_r = read_field('f',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=2)
       f1_z = read_field('f',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=3)
   endif

   if(kar ne 0) then begin
       psi0 = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra)
       psi0_z = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra,$
                           op=3)
       g0 = read_field('i',x,z,t,/equilibrium,_EXTRA=extra)
   end
   n0 = read_field('den',x,z,t,/equilibrium,_EXTRA=extra)
   n0_r = read_field('den',x,z,t,/equilibrium,_EXTRA=extra,op=2)
   n0_z = read_field('den',x,z,t,/equilibrium,_EXTRA=extra,op=3)
   n0_lp = read_field('den',x,z,t,/equilibrium,_EXTRA=extra,op=7)
   
   p0 = read_field('p',x,z,t,/equilibrium,_EXTRA=extra)
   p0_r = read_field('p',x,z,t,/equilibrium,_EXTRA=extra,op=2)
   p0_z = read_field('p',x,z,t,/equilibrium,_EXTRA=extra,op=3)
   p0_lp = read_field('p',x,z,t,/equilibrium,_EXTRA=extra,op=7)

   w0 = read_field('omega',x,z,t,/equilibrium,_EXTRA=extra)

   u1_r = read_field('phi',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=2)
   u1_z = read_field('phi',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=3)
   w1 = read_field('omega',x,z,t,/linear,complex=icomplex,_EXTRA=extra)
   chi1_r = read_field('chi',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=2)
   chi1_z = read_field('chi',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=3)
   chi1_lp = read_field('chi',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=7)
   
   p1 = read_field('p',x,z,t,/linear,complex=icomplex,_EXTRA=extra)
   p1_r = read_field('p',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=2)
   p1_z = read_field('p',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=3)
   p1_lp = read_field('p',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=7)
   n1 = read_field('den',x,z,t,/linear,complex=icomplex,_EXTRA=extra)
   n1_r = read_field('den',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=2)
   n1_z = read_field('den',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=3)
   n1_lp = read_field('den',x,z,t,/linear,complex=icomplex,_EXTRA=extra,op=7)

   ddphi = complex(0,1)*ntor

   if(itor eq 0) then begin
       rzero = read_parameter('rzero',_EXTRA=extra)
       ddphi = ddphi/rzero
       r = 1.
   endif else begin
       r = radius_matrix(x,z,t)
   endelse
   
   term1 = r*(u1_z*p0_r - u1_r*p0_z) $
     - ddphi*p1*w0 - (p0_r*chi1_r + p0_z*chi1_z)/r^2

   term2 = -gamma*p0*(-2.*u1_z + ddphi*w1 + (chi1_lp - chi1_r/r)/r^2)

   divT = ddphi^2*(p1/n0 - n1*p0/n0^2)/r^2 $
     + (p1_lp+p1_r/r)/n0 $
     - 2.*(p1_r*n0_r + p1_z*n0_z + p0_r*n1_r + p0_z*n1_z)/n0^2 $
     - p1*(n0_lp+n0_r/r)/n0^2 - p0*(n1_lp+n1_r/r)/n0^2 $
     + 2.*p1*(n0_r^2 + n0_z^2)/n0^3 + 4.*p0*(n0_r*n1_r + n0_z*n1_z)/n0^3
   ; remove n1 terms
   divT = ddphi^2*(p1/n0)/r^2 $
     + (p1_lp+p1_r/r)/n0 $
     - 2.*(p1_r*n0_r + p1_z*n0_z)/n0^2 $
     - p1*(n0_lp+n0_r/r)/n0^2 $
     + 2.*p1*(n0_r^2 + n0_z^2)/n0^3

   term3 = (gamma-1.)*kap*divT

   if(kar ne 0) then begin
       B1dotT0 = (p0_z*psi1_r - p0_r*psi1_z)/(r*n0) $
         - (n0_z*psi1_r - n0_z*psi1_r)*p0/(r*n0^2) $
         - (p0_r*f1_r + p0_z*f1_z)/n0 $
         + (n0_r*f1_r + n0_z*f1_z)*p0/n0^2
       B0dotT1 = (p1_z*psi0_r - p1_r*psi0_z)/(r*n0) $
         - (n1_z*psi0_r - n1_r*psi0_z)*p0/(r*n0^2) $
         + (g0/r^2)*(ddphi*p1/n0 - ddphi*n1*p0/n0^2)
       B02 = (psi0_r^2 + g0^2)/r^2
       k = kar*(B1dotT0 + B0dotT1)/B02

       term4 = a_bracket(k, psi0, x, z)/r + g0*ddphi*k/r^2
   endif else begin
       term4 = term2*0.
   end

   term5 = 2.*(gamma-1.)*(eta/r^2)* $
     ((psi0_lp-psi0_r/r)*(psi1_lp-psi1_r/r) $
      + (g0_r*(g1_r+ddphi*ddphi*f1_r) + g0_z*(g1_z+ddphi*ddphi*f1_z)) $
      + (ddphi*psi1_z*g0_r - ddphi*psi1_r*g0_z)/r)
   
   term0 = term1+term2+term3+term4+term5

   names = ['total', 'conv', 'comp', 'k_perp', 'k_par', 'ohmic']
   contour_and_legend, real_part([term0, term1, term2, term3, term4, term5]), $
     x,z, title=names, range=[-1e-5, 1e-5]

   print, total(abs(term1)), total(abs(term2)), total(abs(term3))
   print, total(abs(term0))/(total(abs(term1)) + total(abs(term2)))
end
