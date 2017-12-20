pro read_bmncdf, filename=filename, bmn=bmn, psi=psi, m=m, q=q, ntor=ntor, $
                 rho=rho, cur=cur, flux_pol=flux_pol, area=area, bpol=bpol, $
                 symbol=symbol, units=units, version=version, boozer=boozer, $
                 amn=amn, ip=ip, fpol=fpol, jacobian=jacobian, jmn=jmn

  if(n_elements(cur) eq 0) then cur=1.

  if(n_elements(filename) eq 1) then begin

     id = ncdf_open(filename)

     result = ncdf_attinq(id, /global, "version")
     if(result.length eq 0) then begin
        version = 0
     endif else begin
        ncdf_attget, id, /global, "version", version
     end

     bmnr_id = ncdf_varid(id, "bmn_real")
     bmni_id = ncdf_varid(id, "bmn_imag")
     if(version eq 0) then begin
        psi_id = ncdf_varid(id, "psi")
     endif else begin
        psi_id = ncdf_varid(id, "psi_norm")
     end

     symbol = "!8B!Dmn!N"
     units = "G"
     if(version ge 2) then begin
        ncdf_attget, id, "symbol", bytes, /global
        symbol = string(bytes)
        ncdf_attget, id, "units", bytes, /global
        units = string(bytes)
     end
     m_id = ncdf_varid(id, "m")
     q_id = ncdf_varid(id, "q")
     flux_pol_id = ncdf_varid(id, "flux_pol")
     area_id = ncdf_varid(id, "area")
     bp_id = ncdf_varid(id, "Bp")
     F_id = ncdf_varid(id, "F")
     ip_id = ncdf_varid(id, "current")
     jac_id = ncdf_varid(id, "jacobian")
     r_id = ncdf_varid(id, "rpath")
     z_id = ncdf_varid(id, "zpath")
       
     ncdf_attget, id, "ntor", ntor, /global
     ncdf_varget, id, bmnr_id, bmnr
     ncdf_varget, id, bmni_id, bmni
     ncdf_varget, id, psi_id, psi
     ncdf_varget, id, m_id, m
     ncdf_varget, id, q_id, q
     if(area_id ne -1) then  ncdf_varget, id, area_id, area
     ncdf_varget, id, bp_id, bpol
     ncdf_varget, id, flux_pol_id, flux_pol
     if(ip_id ne -1) then ncdf_varget, id, ip_id, ip
     ncdf_varget, id, jac_id, jacobian
     ncdf_varget, id, r_id, rpath
     ncdf_varget, id, z_id, zpath

     if(F_id ne -1) then ncdf_varget, id, F_id, fpol

     if(keyword_set(boozer)) then begin
        alphar_id = ncdf_varid(id, "alpha_real")
        alphai_id = ncdf_varid(id, "alpha_imag")
        ncdf_varget, id, alphar_id, amnr
        ncdf_varget, id, alphai_id, amni
        amn = complex(amnr, amni)*cur
     end

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

     ; calculate resonant current as defined by Callen
     ; Jmn = -(c/4pi) <|Grad(psi0)|> (d/dpsi0) Bmn / (i m)
     gpsipsi = fltarr(n_elements(psi))
     jmn = bmn

     for i=0, n_elements(m)-1 do begin
        jmn[i,*] = deriv(flux_pol,bmn[i,*])
     end

     for j=0, n_elements(psi)-1 do begin
        gpsipsi[j] = total(bpol[*,j]*rpath[*,j]*jacobian[*,j]) / $
                     total(jacobian[*,j])
        jmn[*,j] = -2.998e10/(4.*!pi) * gpsipsi[j] * jmn[*,j] $
                   / (complex(0,1)*m)
     end

     ; convert to MKS
     jmn = jmn / 2.998e5

     return
  endif else begin

     if(n_elements(cur) lt n_elements(filename)) then $
        cur = replicate(cur, n_elementS(filename))

     for i=0, n_elements(filename)-1 do begin
        if(i eq 0) then begin
           read_bmncdf, file=filename[i], bmn=bmn, psi=psi, m=m, q=q, $
                        ntor=ntor, rho=rho, area=area, bpol=bpol, $
                        symbol=symbol, units=units
        endif else begin
           read_bmncdf, file=filename[i], $
                        bmn=bmn1, psi=psi1, m=m1, q=q1, ntor=ntor1, $
                        rho=rho, area=area, bpol=bpol, $
                        symbol=symbol, units=units

           if(ntor1 ne ntor) then begin
              print, "Error: ntor doesn't match"
              continue
           end
           bmn = bmn + bmn1
        endelse
     end
  end
end
