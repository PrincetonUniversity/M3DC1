pro read_bmncdf, file=filename, bmn=bmn, psi=psi, m=m, q=q, ntor=ntor

  if(n_elements(filename) eq 1) then begin

     id = ncdf_open(filename)
     bmnr_id = ncdf_varid(id, "bmn_real")
     bmni_id = ncdf_varid(id, "bmn_imag")
     psi_id = ncdf_varid(id, "psi")
     m_id = ncdf_varid(id, "m")
     q_id = ncdf_varid(id, "q")
       
     ncdf_attget, id, "ntor", ntor, /global
     ncdf_varget, id, bmnr_id, bmnr
     ncdf_varget, id, bmni_id, bmni
     ncdf_varget, id, psi_id, psi
     ncdf_varget, id, m_id, m
     ncdf_varget, id, q_id, q
     
     ncdf_close, id
     
     bmn = complex(bmnr, bmni)
     return
  endif else begin
     
     for i=0, n_elements(filename)-1 do begin
        if(i eq 0) then begin
           read_bmncdf, file=filename[i], bmn=bmn, psi=psi, m=m, q=q, ntor=ntor
        endif else begin
           read_bmncdf, file=filename[i], $
                        bmn=bmn1, psi=psi1, m=m1, q=q1, ntor=ntor1

           if(ntor1 ne ntor) then begin
              print, "Error: ntor doesn't match"
              continue
           end
           bmn = bmn + bmn1
        endelse
     end
  end
end
