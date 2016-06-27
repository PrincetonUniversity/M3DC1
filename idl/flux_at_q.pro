function flux_at_q, qval, normalized_flux=norm, points=pts, $
                    q=q, flux=flux, psi=psi, x=x, z=z, t=t, _EXTRA=extra

  if(n_elements(q) eq 0 or n_elements(flux) eq 0) then begin
;     q = flux_average('q', flux=flux, nflux=nflux, /equilibrium, points=pts, $
;                      _EXTRA=extra, psi=psi, x=x, z=z, t=t)
     fc = flux_coordinates(_EXTRA=extra, points=pts, psi0=psi, x=x, z=z, $
                          /equilibrium)
     q = abs(fc.q)
     if(keyword_set(norm)) then flux=fc.psi_norm
  end
  
   dq_dpsi = deriv(flux,q)

   n = n_elements(qval)
   if(n eq 0) then return, 0

   fval = fltarr(n)
   for k=0, n-1 do begin
       ; make initial guess
       dum = min(q-qval[k],/abs,i)
       fval[k] = flux[i]
       
       ; perform newton iterations to refine result
       for j=0, 5 do begin
           dq = interpol(dq_dpsi,flux,fval[k]) 
           q0 = interpol(q,flux,fval[k])
           dpsi = (qval[k] - q0)/dq_dpsi[i]
           fval[k] = fval[k] + dpsi
           if((fval[k] gt max(flux)) or (fval[k] lt min(flux))) then begin
               print, 'flux_at_q: could find surface with q = ', $
                 qval[k]
               break
           endif
;           print, qval[k], q0, fval[k]
       end
   end

   return, fval
end
