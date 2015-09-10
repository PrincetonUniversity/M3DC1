pro read_bmncdf, file=filename, _EXTRA=extra, bmn=bmn, psi=psi, m=m, q=q, $
                 ntor=ntor

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
end
