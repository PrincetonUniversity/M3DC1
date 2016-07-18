pro read_bmncdf, file=filename, bmn=bmn, psi=psi, m=m, q=q, ntor=ntor, $
                 rho=rho, cur=cur

  if(n_elements(cur) eq 0) then cur=1.

  if(n_elements(filename) eq 1) then begin

     id = ncdf_open(filename)
     bmnr_id = ncdf_varid(id, "bmn_real")
     bmni_id = ncdf_varid(id, "bmn_imag")
     psi_id = ncdf_varid(id, "psi")
     m_id = ncdf_varid(id, "m")
     q_id = ncdf_varid(id, "q")
     flux_pol_id = ncdf_varid(id, "flux_pol")
       
     ncdf_attget, id, "ntor", ntor, /global
     ncdf_varget, id, bmnr_id, bmnr
     ncdf_varget, id, bmni_id, bmni
     ncdf_varget, id, psi_id, psi
     ncdf_varget, id, m_id, m
     ncdf_varget, id, q_id, q
     ncdf_varget, id, flux_pol_id, flux_pol
     
     ncdf_close, id
     
     bmn = complex(bmnr, bmni)*cur

     ; calculate rho
     dflux_tor = deriv(flux_pol)
     flux_tor = flux_pol
     flux_tor[0] = 0.
     for i=1, n_elements(flux_pol)-1 do begin
        flux_tor[i] = flux_tor[i-1] + $
                      (q[i-1]+q[i])*(dflux_tor[i-1] + dflux_tor[i])/4.
     end
     rho = sqrt(flux_tor / flux_tor(n_elements(flux_pol)-1))
     return
  endif else begin

     if(n_elements(cur) lt n_elements(filename)) then $
        cur = replicate(cur, n_elementS(filename))

     for i=0, n_elements(filename)-1 do begin
        if(i eq 0) then begin
           read_bmncdf, file=filename[i], bmn=bmn, psi=psi, m=m, q=q, $
                        ntor=ntor, rho=rho
        endif else begin
           read_bmncdf, file=filename[i], $
                        bmn=bmn1, psi=psi1, m=m1, q=q1, ntor=ntor1, $
                        rho=rho

           if(ntor1 ne ntor) then begin
              print, "Error: ntor doesn't match"
              continue
           end
           bmn = bmn + bmn1
        endelse
     end
  end
end
