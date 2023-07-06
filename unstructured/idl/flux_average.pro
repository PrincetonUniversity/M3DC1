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
;   area:  surface area of each flux surface
;   name:  the formatted name of the field
; symbol:  the formatted symbol of the field
;  units:  the formatted units of the field
;==================================================================
function flux_average, field, psi=psi, i0=i0, x=x, z=z, t=t, r0=r0, $
                       flux=flux, nflux=nflux, area=area, bins=bins, $
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
      fc = flux_coordinates(points=points,filename=filename,$
                            tbins=points,fbins=bins,itor=itor,r0=r0, $
                            psi0=psi,i0=i0,x=x,z=z,_EXTRA=extra)
      if(isa(fc, "Int")) then begin
         print, 'Error calculating flux coordinates'
         return, 0
      end
   end
   flux = fc.psi
   nflux = fc.psi_norm

   if(type eq 7) then begin ; named field
       if (strcmp(field, 'Safety Factor', /fold_case) eq 1) or $
         (strcmp(field, 'q', /fold_case) eq 1) then begin

           units = ''
           name = '!6Safety Factor!X'
           symbol = '!8q!X'

           return, abs(fc.q)

       endif else $
         if(strcmp(field, 'alpha', /fold_case) eq 1) then begin

           p = flux_average('p', psi=psi, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, bins=bins, $
             points=points, last=last, filename=filename, fc=fc, _EXTRA=extra)

           pp = deriv(fc.psi, p)
           dV = fc.dV_dchi / fc.dpsi_dchi

           units = ''
           name = '!6Ballooning Parameter!X'
           symbol = '!7a!X'

           ; From Miller, et al. 1998
           return, -dV/(2.*!pi^2) * sqrt(fc.V/(2.*!pi^2*fc.r0)) * pp

       endif else $
         if(strcmp(field, 'shear', /fold_case) eq 1) then begin

           q = flux_average('q', psi=psi, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, bins=bins, $
             points=points, last=last, filename=filename, fc=fc, _EXTRA=extra)
           
           dqdV = deriv(fc.V, q)

           units = ''
           name = '!6Magnetic Shear!X'
           symbol = '!8s!X'

           ; From Miller, et al. 1998
           return, 2.*fc.V*dqdV/q
       endif else $
         if(strcmp(field, 'elongation', /fold_case) eq 1) then begin

           psi = flux_average('psi', psi=psi, x=x, z=z, t=t, $
                              r0=r0, flux=flux, nflux=nflux, area=area, bins=bins, $
                              points=points, last=last, elongation=elongation, filename=filename, $
                              fc=fc, _EXTRA=extra)

           units = ''
           name = '!6Elongation!X'
           symbol = '!7j!X'

           return, elongation

       endif else $
         if(strcmp(field, 'dqdrho', /fold_case) eq 1) then begin

           q = flux_average('q', psi=psi, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, bins=bins, $
             points=points, last=last, filename=filename, fc=fc, _EXTRA=extra)
           
           units = make_units(l0=-1, filename=filename, _EXTRA=extra)
           name = '!8dq!3/!8d!7q!X'
           symbol = name

           return, deriv(fc.rho, q)

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
             r0=r0, flux=flux, nflux=nflux, area=area, bins=bins, $
             points=points, /equilibrium, filename=filename, _EXTRA=extra)

           qp = deriv(flux,q)
           qrz = interpol(q, flux, psi)
           qprz = interpol(qp, flux, psi)
           jac = r^2*qrz/I

           dpsi = sqrt(s_bracket(psi,psi,x,z))
           dtheta = I/(r*qrz*dpsi)

           gpsi = flux_average(dpsi/jac, psi=psi, i0=i, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, bins=bins, $
             points=points, filename=filename, fc=fc, _EXTRA=extra)
           gchi = flux_average(dtheta/jac, psi=psi, i0=i, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, bins=bins, $
             points=points, filename=filename, fc=fc, _EXTRA=extra)
           gsp = flux_average(sp/jac, psi=psi, i0=i, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, bins=bins, $
             points=points, filename=filename, fc=fc, _EXTRA=extra)
           gi = flux_average(I, psi=psi, i0=i, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, bins=bins, $
             points=points, filename=filename, fc=fc, _EXTRA=extra)
           
           units = ''
           name = '!8m!7k!X'
           symbol = name          

           return, -gi*q*gsp/(2.*qp)*(1./(gpsi*gchi))

       endif else $
         if(strcmp(field, 'rho', /fold_case) eq 1) then begin

           units = make_units(/l0, filename=filename, _EXTRA=extra)
           name = '!7q!X'
           symbol = name

           return, fc.rho

       endif else $
         if(strcmp(field, 'flux_t', /fold_case) eq 1) then begin

           units = ''
           name = '!6Toroidal Flux!X'
           symbol = '!7W!D!8t!N!X'

           return, fc.flux_tor

        endif else $
         if(strcmp(field, 'flux_p', /fold_case) eq 1) then begin

           units = ''
           name = '!6Poloidal Flux!X'
           symbol = '!7W!D!8p!N!X'

           return, fc.flux_pol

        endif else $
         if(strcmp(field, 'volume', /fold_case) eq 1) then begin

           units = make_units(l0=3,filename=filename,_EXTRA=extra)
           name = '!6Volume!X'
           symbol = '!8V!X'

           return, fc.V

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
             flux=flux, area=area, bins=bins, filename=filename, $
                                   fc=fc, _EXTRA=extra)
           rr = flux_average_field(r,psi,x,z,t, r0=r0, $
             flux=flux, area=area, bins=bins, filename=filename, $
                                   fc=fc, _EXTRA=extra)
           bb = flux_average_field(bpol2,psi,x,z,t, r0=r0, $
             flux=flux, area=area, bins=bins, filename=filename, $
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
                              nflux=nflux, area=area, bins=bins, $
                              integrate=integrate, filename=filename, _EXTRA=extra)

           symbol = '!7a!X'
           units = ''
           name = '!7a!X'

           return, -q^2*fc.r0*alpha

       endif else $
         if(strcmp(field, 'kappa_implied', /fold_case) eq 1) then begin

           Q_ext = read_field('heat_source', x, z, t, points=points, complex=complex, $
                              symbol=symbol, units=units,dimensions=d,fac=fac,$
                              abs=abs, phase=phase, filename=filename, $
                              _EXTRA=extra)
           eta = read_field('eta',x,z,t,points=points, $
                            last=last,filename=filename,_EXTRA=extra)
           Jx = read_field('jx',x,z,t,points=points, $
                          last=last,filename=filename,_EXTRA=extra)
           Jy = read_field('jy',x,z,t,points=points, $
                          last=last,filename=filename,_EXTRA=extra)
           Jz = read_field('jz',x,z,t,points=points, $
                          last=last,filename=filename,_EXTRA=extra)
           J2 = Jx^2 + Jy^2 + Jz^2

           te = read_field('Te',x,z,t,points=points,last=last,filename=filename, _EXTRA=extra)
           ti = read_field('Ti',x,z,t,points=points,last=last,filename=filename, _EXTRA=extra)
           temp = te + ti

           pprime = s_bracket(temp,psi,x,z)
           GradP = flux_average_field(pprime, psi, x, z, t, psi=psi, i0=i, $
                                      r0=r0, flux=flux, $
                                      nflux=nflux, area=area, fc=fc, $
                                      bins=bins, filename=filename, _EXTRA=extra)
           dV = fc.dV_dchi / fc.dpsi_dchi

           qq = Q_ext + eta*J2
           Q = flux_average_field(qq, psi, x, z, t, fc=fc, /integrate, psi=psi, i0=i, $
                                   r0=r0, flux=flux, nflux=nflux, area=area, bins=bins, $
                                   points=points, filename=filename, _EXTRA=extra)

           symbol = '!7j!X'
           d = dimensions(l0=2, t0=-1, n0=1)
           units = parse_units(d, _EXTRA=extra)
           name = '!7j!X'

           return, -Q/(dV*GradP)

       endif else $
          if(strcmp(field, 'amu_implied', /fold_case) eq 1) then begin
           Q =  flux_average('torque', psi=psi, i0=i, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, bins=bins, $
             points=points, filename=filename, fc=fc, _EXTRA=extra, /integrate)
           w= read_field('omega',x,z,t,points=points,last=last,filename=filename, _EXTRA=extra)
           r = radius_matrix(x,z,t)
           wprime = s_bracket(w*r,psi,x,z)
           GradW = flux_average_field(wprime, psi, x, z, t, r0=r0, flux=flux, $
                                      nflux=nflux, area=area, fc=fc, $
                                      bins=bins, filename=filename, _EXTRA=extra)

           dV = fc.dV_dchi / fc.dpsi_dchi

           symbol = '!7l!X'
           d = dimensions(/p0, /t0)
           units = parse_units(d, _EXTRA=extra)
           name = '!7l!X'

           return, -Q/(dV*GradW)

        endif else $
          if(strcmp(field, 'denm_implied', /fold_case) eq 1) then begin
           Q =  flux_average('sigma', psi=psi, i0=i, x=x, z=z, t=t, $
             r0=r0, flux=flux, nflux=nflux, area=area, bins=bins, $
             points=points, filename=filename, fc=fc, _EXTRA=extra, /integrate)
           w= read_field('den',x,z,t,points=points,last=last,filename=filename, _EXTRA=extra)
           r = radius_matrix(x,z,t)
           wprime = s_bracket(w,psi,x,z)
           GradW = flux_average_field(wprime, psi, x, z, t, r0=r0, flux=flux, $
                                      nflux=nflux, area=area, fc=fc, $
                                      bins=bins, filename=filename, _EXTRA=extra)

           dV = fc.dV_dchi / fc.dpsi_dchi

           symbol = '!8D!Dn!N!X'
           d = dimensions(l0=2, t0=-1)
           units = parse_units(d, _EXTRA=extra)
           name = '!8D!Dn!N!X'

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
          ;;                             nflux=nflux, area=area, $
          ;;                             bins=bins, filename=filename, $
          ;;                             _EXTRA=extra, integrate=integrate)

          ;;  symbol = '!3|!6B!3|!6!U2!N!X'
          ;;  d = dimensions(b0 = 2)
          ;;  units = parse_units(d, _EXTRA=extra)
          ;;  name = '!3|!6B!3|!6!U2!N!X'

          ;;  return, b2_fa
        endif else $
           if((strcmp(field, 'd_i', /fold_case) eq 1) or $
              (strcmp(field, 'd_r', /fold_case) eq 1)) then begin

           i = read_field('i',x,z,t,points=points, $
                          last=last,filename=filename,_EXTRA=extra)
           psi = read_field('psi',x,z,t,points=points, $
                            last=last,filename=filename,_EXTRA=extra)
           psi_r = read_field('psi',x,z,t,points=points, $
                            last=last,filename=filename,op=2,_EXTRA=extra)
           psi_z = read_field('psi',x,z,t,points=points, $
                            last=last,filename=filename,op=3,_EXTRA=extra)
           p_r = read_field('p',x,z,t,points=points, $
                            last=last,filename=filename,op=2,_EXTRA=extra)
           p_z = read_field('p',x,z,t,points=points, $
                            last=last,filename=filename,op=3,_EXTRA=extra)

           if(itor eq 1) then begin
              r = radius_matrix(x,z,t)
           endif else begin
              r = 1.
           end

           psi2 = psi_r^2 + psi_z^2
           b2 = psi2/r^2 + i^2/r^2
           pprime = (psi_r*p_r + psi_z*p_z)/psi2

           b2_av = $
              flux_average_field(b2, psi, x, z, t, psi=psi, i0=i, $
                                 r0=r0, flux=flux, $
                                 nflux=nflux, area=area, fc=fc, $
                                 bins=bins, filename=filename, $
                                 _EXTRA=extra)
           b2_psi2_av = $
              flux_average_field(b2/psi2, psi, x, z, t, psi=psi, i0=i, $
                                 r0=r0, flux=flux, $
                                 nflux=nflux, area=area, fc=fc, $
                                 bins=bins, filename=filename, $
                                 _EXTRA=extra)
           inv_b2psi2_av = $
              flux_average_field(1./(b2*psi2), psi, x, z, t, psi=psi, i0=i, $
                                 r0=r0, flux=flux, $
                                 nflux=nflux, area=area, fc=fc, $
                                 bins=bins, filename=filename, $
                                 _EXTRA=extra)
           inv_psi2_av = $
              flux_average_field(1./psi2, psi, x, z, t, psi=psi, i0=i, $
                                 r0=r0, flux=flux, $
                                 nflux=nflux, area=area, fc=fc, $
                                 bins=bins, filename=filename, $
                                 _EXTRA=extra)
           inv_b2_av = $
              flux_average_field(1./b2, psi, x, z, t, psi=psi, i0=i, $
                                 r0=r0, flux=flux, $
                                 nflux=nflux, area=area, fc=fc, $
                                 bins=bins, filename=filename, $
                                 _EXTRA=extra)
           pp = $
              flux_average_field(pprime, psi, x, z, t, psi=psi, i0=i, $
                                 r0=r0, flux=flux, $
                                 nflux=nflux, area=area, fc=fc, $
                                 bins=bins, filename=filename, $
                                 _EXTRA=extra)
           g = $
              flux_average_field(i, psi, x, z, t, psi=psi, i0=i, $
                                 r0=r0, flux=flux, $
                                 nflux=nflux, area=area, fc=fc, $
                                 bins=bins, filename=filename, $
                                 _EXTRA=extra)

           Vp_2pi2 = fc.dV_dchi / fc.dpsi_dchi / (2.*!pi)^2
           Vpp_2pi2 = deriv(fc.psi, Vp_2pi2)
           q = fc.q
           qp = deriv(fc.psi,q)

           ; here the first term of E and all of H have the opposite
           ; sign of the standard derivation because M3D-C1 defines
           ; B = grad(psi) x grad(phi) rather than grad(phi)xgrad(psi)
           ; but q is still defined as positive when IP and BT are in
           ; the same direction
           e = -pp*Vp_2pi2/(qp^2) * b2_psi2_av *$
               (qp*g / b2_av + Vpp_2pi2)
           f = (pp*Vp_2pi2/qp)^2 * $
               (g^2*(b2_psi2_av*inv_b2psi2_av - inv_psi2_av^2) $
                + b2_psi2_av*inv_b2_av)
           h = -g*pp*Vp_2pi2/qp * (inv_psi2_av - b2_psi2_av / b2_av)

           d = dimensions()
           units = parse_units(d, _EXTRA=extra)

           if(strcmp(field, 'd_i', /fold_case) eq 1) then begin
              symbol = '!8D!DI!N!X'
              name = '!6Ideal Interchange Criterion!X'
              return, e + f + h - 1./4.
           endif else begin
              symbol = '!8D!DR!N!X'
              name = '!6Resistive Interchange Criterion!X'
              return, e + f + h^2
           endelse

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
                           nflux=nflux, area=area, bins=bins, $
                           integrate=integrate, surface_weight=total, $
                           elongation=elongation, fc=fc, $
                           filename=filename, _EXTRA=extra)

   if(keyword_set(total)) then begin
       fa = fa*area
   end

   return, fa
end
