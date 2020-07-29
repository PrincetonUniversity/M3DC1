function flux_at_q, qval, q, normalized_flux=norm, points=pts, $
                    q=qout, flux=flux, psi=psi, x=x, z=z, t=t, $
                    guess=guess, fc=fc, unique=unique, _EXTRA=extra

  if(n_elements(q) eq 0 or n_elements(flux) eq 0) then begin
     if(not isa(fc)) then begin
        fc = flux_coordinates(_EXTRA=extra, points=pts, psi0=psi, x=x, z=z, $
                              /equilibrium)
     end
     q = abs(fc.q)
     if(keyword_set(norm)) then flux=fc.psi_norm else flux=fc.psi
  end

  
   dq_dpsi = deriv(flux,q)

   n = n_elements(qval)
   if(n eq 0) then return, 0

   if(n gt 1) then begin
      for k=0, n-1 do begin
         temp_f = flux_at_q(qval[k], q, normalized_flux=norm, points=pts, $
                  q=temp_q, fc=fc, flux=flux, _EXTRA=extra)
         if(temp_f[0] eq 0) then continue
         
         if(n_elements(fval) eq 0) then begin
            fval = temp_f
            qout = temp_q
         endif else begin
            fval = [fval, temp_f]
            qout = [qout, temp_q]
         end
      end
      if(n_elements(fval) eq 0) then return, 0
      return, fval
   end
   
   ; do point-by-point search
   nv = 0
   for i=1, n_elements(q)-1 do begin
      if((q[i]-qval)*(q[i-1]-qval) gt 0) then continue

      f0 = flux[i]

      ; perform newton iterations to refine result
      for j=0, 5 do begin
         dq = interpol(dq_dpsi,flux,f0) 
         q0 = interpol(q,flux,f0)
         dpsi = (qval - q0)/dq
         f0 = f0 + dpsi
         if((f0 gt max(flux)) or (f0 lt min(flux))) then begin
            print, 'flux_at_q: error converging on surface with q = ', qval
            break
         endif
;           print, qval[k], q0, fval[k]
      end

      if(nv eq 0) then begin
         qout = interpol(q, flux, f0)
         fval = f0
      endif else begin
         qout = [qout, interpol(q, flux, f0)]
         fval = [fval, f0]
      end
      nv = nv + 1

      if(keyword_set(unique)) then break
   end


   if(nv eq 0) then begin
      print, 'flux_at_q: could not find surface with q = ', qval
      return, 0
   end
   return, fval
end
