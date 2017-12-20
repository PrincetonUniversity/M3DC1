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
       
;       for i=0, n_elements(psi)-1 do begin
;           bmni[*,i] = bmni[*,i]*(-1)^abs(m)
;           bmnr[*,i] = bmnr[*,i]*(-1)^abs(m)
;       end
           
;       angle = atan(bmni, bmnr) + !pi
;       k = where(angle lt 0.)
;       angle[k] = angle[k] + 2.*!pi

;       bmn = fltarr(1, n_elements(m), n_elements(psi))
;       bmn[0,*,*] = sqrt(bmnr^2 + bmni^2)

;       dum = min(psi-0.89,i,/abs)
;       print, i



;       plot, m, angle[*,i], _EXTRA=extra
;       plot, m, bmnr[*,i], _EXTRA=extra
;       oplot, m, bmni[*,i], linestyle=1, _EXTRA=extra


;       contour_and_legend, bmn, m, psi, _EXTRA=extra, table=39
;       ct3
;       oplot, -ntor*q, psi, color=color(5)
end
