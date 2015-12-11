;==================================================================
; flux_average
; ~~~~~~~~~~~~
;
; Computes the flux averages of a field at a given time
;
; input:
;   field: the field to flux-average.  This may be a string (the
;          name of the field) or the field values themselves.
;   psi:   the flux field
;   x:     r-coordinates
;   z:     z-coordinates
;   t:     t-coordinates;
;   bins:  number of bins into which to subdivide the flux
;   range: range of flux values
;
; output: 
;   flux:  flux value of each bin
;     dV:  differential volume of each flux surface
;   area:  surface area of each flux surface
;   name:  the formatted name of the field
; symbol:  the formatted symbol of the field
;  units:  the formatted units of the field
;==================================================================
function flux_average, field, psi=psi, i0=i0, x=x, z=z, t=t, r0=r0, $
                       flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
                       points=points, name=name, symbol=symbol, units=units, $
                       integrate=integrate, complex=complex, abs=abs, $
                       phase=phase, stotal=total, fac=fac, fc=fc, $
                       elongation=elongation, filename=filename, $
                       _EXTRA=extra

   type = size(field, /type)
   sz = size(field)

   if(n_elements(points) eq 0) then begin
       if(type ne 7) then points = sz[2] else points=200
   endif

   if(not isa(fc)) then begin
      itor = read_parameter('itor', filename=filename)
      r0 =read_parameter('rzero', filename=filename)
      fc = flux_coordinates(/fast,points=points,filename=filename,$
                            tbins=bins,fbins=bins,itor=itor,r0=r0, $
                            _EXTRA=extra)
      if(isa(fc, /int)) then begin
         print, 'Error calculating flux coordinates'
         return, 0
      end
   end

   if(type eq 7) then begin ; named field
       if (strcmp(field, 'Safety Factor', /fold_case) eq 1) or $
         (strcmp(field, 'q', /fold_case) eq 1) then begin

           flux_t = flux_average('flux_t', psi=psi, i0=i0, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, last=last, filename=filename, fc=fc, _EXTRA=extra)
           
           units = ''
           name = '!6Safety Factor!X'
           symbol = '!8q!X'

           return, abs(deriv(flux, flux_t))/fc.period

       endif else $
         if(strcmp(field, 'alpha', /fold_case) eq 1) then begin

           p = flux_average('p', psi=psi, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, last=last, filename=filename, fc=fc, _EXTRA=extra)

           pp = deriv(fc.psi, p)

           units = ''
           name = '!6Ballooning Parameter!X'
           symbol = '!7a!X'

           return, -fc.dV/(2.*!pi^2) * sqrt(abs(fc.V)/(2.*!pi^2*fc.r0)) * pp

       endif else $
         if(strcmp(field, 'shear', /fold_case) eq 1) then begin

           q = flux_average('q', psi=psi, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, last=last, filename=filename, fc=fc, _EXTRA=extra)
           
           dqdV = deriv(fc.V, q)

           units = ''
           name = '!6Magnetic Shear!X'
           symbol = '!8s!X'

           return, 2.*fc.V*dqdV/q
       endif else $
         if(strcmp(field, 'elongation', /fold_case) eq 1) then begin

           psi = flux_average('psi', psi=psi, x=x, z=z, t=t, $
                              r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
                              points=points, last=last, elongation=elongation, filename=filename, $
                              fc=fc, _EXTRA=extra)

           units = ''
           name = '!6Elongation!X'
           symbol = '!7j!X'

           return, elongation

       endif else $
         if(strcmp(field, 'dqdrho', /fold_case) eq 1) then begin

           q = flux_average('q', psi=psi, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, last=last, filename=filename, fc=fc, _EXTRA=extra)
           rho = flux_average('rho', psi=psi, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, last=last, filename=filename, fc=fc, _EXTRA=extra)
           
           units = make_units(l0=-1, filename=filename, _EXTRA=extra)
           name = '!8dq!3/!8d!7q!X'
           symbol = name

           return, deriv(rho, q)

       endif else $
         if(strcmp(field, 'lambda', /fold_case) eq 1) then begin
           ; this is m*lambda from Hegna, Callen Phys. Plasmas 1 (1994) p.2308

           I = read_field('I', x, z, t, points=points, /equilibrium, $
                          units=units, filename=filename, _EXTRA=extra)
           sigma = read_field('jpar_B', x, z, t, points=points, /equilibrium, $
                              units=units, filename=filename, _EXTRA=extra)
           sp = s_bracket(psi,sigma,x,z)/s_bracket(psi,psi,x,z)
           r = radius_matrix(x,z,t)

           q = flux_average('q', psi=psi, i0=i, x=x, z=z, t=t, fc=fc, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, /equilibrium, filename=filename, _EXTRA=extra)

           qp = deriv(flux,q)
           qrz = interpol(q, flux, psi)
           qprz = interpol(qp, flux, psi)
           jac = r^2*qrz/I

           dpsi = sqrt(s_bracket(psi,psi,x,z))
           dtheta = I/(r*qrz*dpsi)

           gpsi = flux_average(dpsi/jac, psi=psi, i0=i, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, filename=filename, fc=fc, _EXTRA=extra)
           gchi = flux_average(dtheta/jac, psi=psi, i0=i, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, filename=filename, fc=fc, _EXTRA=extra)
           gsp = flux_average(sp/jac, psi=psi, i0=i, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, filename=filename, fc=fc, _EXTRA=extra)
           gi = flux_average(I, psi=psi, i0=i, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, filename=filename, fc=fc, _EXTRA=extra)
           
           units = ''
           name = '!8m!7k!X'
           symbol = name          

           return, -gi*q*gsp/(2.*qp)*(1./(gpsi*gchi))

       endif else $
         if(strcmp(field, 'rho', /fold_case) eq 1) then begin

           flux_t = flux_average('flux_t', psi=psi, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, last=last, filename=filename, fc=fc, _EXTRA=extra)
           bzero = read_parameter('bzero', filename=filename, _EXTRA=extra)
           print, 'bzero = ', bzero          

           units = make_units(/l0, filename=filename, _EXTRA=extra)
           name = '!7q!X'
           symbol = name

           return, sqrt(flux_t/(!pi*bzero))

       endif else $
         if(strcmp(field, 'flux_t', /fold_case) eq 1) then begin
           print, 'DBG: flux_t reading field'

           if(n_elements(i0) le 1) then begin
               i0 = read_field('I', x, z, t, points=points, filename=filename, _EXTRA=extra)
           endif
           
           itor = read_parameter('itor', filename=filename)

           r = radius_matrix(x,z,t)

           if(itor eq 1) then begin
              field = i0/(fc.period*r^2)
           endif else begin
              field = i0/(fc.period)
           endelse

           units = ''
           name = '!6Toroidal Flux!X'
           symbol = '!7w!D!8t!N!X'

           integrate = 1

       endif else $
         if(strcmp(field, 'beta_pol', /fold_case) eq 1) then begin
           I = read_field('I', x, z, t, points=points, filename=filename, _EXTRA=extra)
           
           bzero = read_parameter('bzero',filename=filename, _EXTRA=extra)
           rzero = read_parameter('rzero',filename=filename, _EXTRA=extra)
           izero = bzero*rzero
           print, 'izero = ', izero
       
           r = radius_matrix(x,z,t)

           bpol2 = s_bracket(psi,psi,x,z)/r^2

           ii = flux_average_field(izero^2-i^2,psi,x,z,t, r0=r0, $
             flux=flux, area=area, dV=dV, bins=bins, filename=filename, $
                                   fc=fc, _EXTRA=extra)
           rr = flux_average_field(r,psi,x,z,t, r0=r0, $
             flux=flux, area=area, dV=dV, bins=bins, filename=filename, $
                                   fc=fc, _EXTRA=extra)
           bb = flux_average_field(bpol2,psi,x,z,t, r0=r0, $
             flux=flux, area=area, dV=dV, bins=bins, filename=filename, $
                                   fc=fc, _EXTRA=extra)
         
           symbol = '!7b!6!Dpol!N!X'
           units = ''
           name = '!6Poloidal Beta!X'

           return, 1.+0.5*ii/(bb*rr^2)
           
       endif else $
         if(strcmp(field, 'alpha2', /fold_case) eq 1) then begin
           q = flux_average('q',psi=psi,x=x,z=z,nflux=nflux, $
                            flux=flux, bins=bins, r0=r0, fc=fc, $
                            points=points, last=last, filename=filename, _EXTRA=extra)
 
           beta = read_field('beta',x,z,t,points=points,last=last,filename=filename, _EXTRA=extra)
           r = read_field('r',x,z,t,points=points,last=last,filename=filename, _EXTRA=extra)

           betap = s_bracket(beta,r,x,z)
           alpha = $
             flux_average_field(betap, psi, x, z, t, r0=r0, flux=flux, fc=fc, $
                              nflux=nflux, area=area, dV=dV, bins=bins, $
                              integrate=integrate, filename=filename, _EXTRA=extra)

           symbol = '!7a!X'
           units = ''
           name = '!7a!X'

           return, -q^2*fc.r0*alpha

       endif else $
         if(strcmp(field, 'kappa_implied', /fold_case) eq 1) then begin
           Q =  flux_average('heat_source', psi=psi, i0=i, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, filename=filename, _EXTRA=extra, fc=fc, /integrate)
           p = read_field('p',x,z,t,points=points,last=last,filename=filename, _EXTRA=extra)
           n = read_field('den',x,z,t,points=points,last=last,filename=filename, _EXTRA=extra)
           temp = p/n

           pprime = s_bracket(temp,psi,x,z)
           GradP = flux_average_field(pprime, psi, x, z, t, r0=r0, flux=flux, $
                                      nflux=nflux, area=area, dV=dV, fc=fc, $
                                      bins=bins, filename=filename, _EXTRA=extra)

           symbol = '!7j!X'
           d = dimensions(l0=2, t0=-1, n0=1)
           units = parse_units(d, _EXTRA=extra)
           name = '!7j!X'

           return, -Q/(dV*GradP)

        endif else $
          if(strcmp(field, 'amu_implied', /fold_case) eq 1) then begin
           Q =  flux_average('torque', psi=psi, i0=i, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, dV=dV, bins=bins, $
             points=points, filename=filename, fc=fc, _EXTRA=extra, /integrate)
           w= read_field('omega',x,z,t,points=points,last=last,filename=filename, _EXTRA=extra)
           r = radius_matrix(x,z,t)
           wprime = s_bracket(w*r,psi,x,z)
           GradW = flux_average_field(wprime, psi, x, z, t, r0=r0, flux=flux, $
                                      nflux=nflux, area=area, dV=dV, fc=fc, $
                                      bins=bins, filename=filename, _EXTRA=extra)

           symbol = '!7l!X'
           d = dimensions(/p0, /t0)
           units = parse_units(d, _EXTRA=extra)
           name = '!7l!X'

           return, -Q/(dV*GradW)

       ;; endif else $
          ;; if(strcmp(field, 'b2', /fold_case) eq 1) then begin

          ;;  bx = read_field('bx', x, z, t, points=points, complex=complex, $
          ;;                  symbol=symbol, units=units,dimensions=d,$
          ;;                  abs=abs, phase=phase, filename=filename, $
          ;;                  _EXTRA=extra)
          ;;  by = read_field('by', x, z, t, points=points, complex=complex, $
          ;;                  symbol=symbol, units=units,dimensions=d,$
          ;;                  abs=abs, phase=phase, filename=filename, $
          ;;                  _EXTRA=extra)
          ;;  bz = read_field('bz', x, z, t, points=points, complex=complex, $
          ;;                  symbol=symbol, units=units,dimensions=d,$
          ;;                  abs=abs, phase=phase, filename=filename, $
          ;;                  _EXTRA=extra)
          ;;  b2 = bx*conj(bx) + by*conj(by) + bz*conj(bz)

          ;;  b2_fa = flux_average_field(b2, psi, x, z, t, r0=r0, flux=flux, $
          ;;                             nflux=nflux, area=area, dV=dV, $
          ;;                             bins=bins, filename=filename, $
          ;;                             _EXTRA=extra, integrate=integrate)

          ;;  symbol = '!3|!6B!3|!6!U2!N!X'
          ;;  d = dimensions(b0 = 2)
          ;;  units = parse_units(d, _EXTRA=extra)
          ;;  name = '!3|!6B!3|!6!U2!N!X'

          ;;  return, b2_fa

       endif else begin
           field = read_field(field, x, z, t, points=points, complex=complex, $
                              symbol=symbol, units=units,dimensions=d,fac=fac,$
                              abs=abs, phase=phase, filename=filename, $
                              _EXTRA=extra)
           name = symbol
           if(keyword_set(total)) then begin
               d = d + dimensions(l0=2)
               units = parse_units(d, _EXTRA=extra)
           endif else if(keyword_set(integrate)) then begin
               d = d + dimensions(l0=3)
               units = parse_units(d, _EXTRA=extra)
           endif else begin
               symbol = '!12<!X' + symbol + '!12>!X'
           endelse
       endelse
   endif else begin
       name = ''
       symbol = ''
       units = make_units()
   endelse

   fa = flux_average_field(field, psi, x, z, t, r0=r0, flux=flux, $
                           nflux=nflux, area=area, dV=dV, bins=bins, $
                           integrate=integrate, surface_weight=total, $
                           elongation=elongation, fc=fc, $
                           filename=filename, _EXTRA=extra)

   if(keyword_set(total)) then begin
       fa = fa*area
   end

   return, fa
end
