function read_neo_input, filename
  id = ncdf_open(filename)

  npsi_id = ncdf_dimid(id, 'npsi')  ; number of psi points for profiles
  nr_id = ncdf_dimid(id, 'nr')      ; number of radial points for surfaces
  np_id = ncdf_dimid(id, 'np')      ; number of poloidal points
  nt_id = ncdf_dimid(id, 'nt')      ; number of toroidal points

  ncdf_diminq, id, npsi_id, dimstr, npsi
  ncdf_diminq, id, nr_id, dimstr, nr
  ncdf_diminq, id, np_id, dimstr, np
  ncdf_diminq, id, nt_id, dimstr, nt

  print, 'Dimensions'
  print, ' npsi = ', npsi
  print, ' nr = ', nr
  print, ' np = ', np
  print, ' nt = ', nt

  ncdf_attget, id, 'version', version, /global
  ncdf_attget, id, 'psi_0', psi_0, /global
  ncdf_attget, id, 'psi_1', psi_1, /global
  ncdf_attget, id, 'ion_mass', ion_mass, /global

  q_id = ncdf_varid(id, 'q')
  psi_id = ncdf_varid(id, 'psi')
  phi_id = ncdf_varid(id, 'Phi')
  psi0_id = ncdf_varid(id, 'psi0')
  te0_id = ncdf_varid(id, 'Te0')
  ti0_id = ncdf_varid(id, 'Ti0')
  ne0_id = ncdf_varid(id, 'ne0')
  ni0_id = ncdf_varid(id, 'ni0')
  r_id = ncdf_varid(id, 'R')
  z_id = ncdf_varid(id, 'Z')

  ncdf_varget, id, q_id, q
  ncdf_varget, id, psi_id, psi
  ncdf_varget, id, phi_id, Phi
  ncdf_varget, id, psi0_id, psi0
  ncdf_varget, id, te0_id, Te0
  ncdf_varget, id, ti0_id, Ti0
  ncdf_varget, id, ne0_id, ne0
  ncdf_varget, id, ni0_id, ni0
  ncdf_varget, id, r_id, R
  ncdf_varget, id, z_id, Z

  ncdf_close, id

  a = { npsi:npsi, nr:nr, np:np, nt:nt, $
        version:version, psi_0:psi_0, psi_1:psi_1, ion_mass:ion_mass, $
        q:q, psi:psi, Phi:Phi, psi0:psi0, Te0:Te0, Ti0:Ti0, ni0:ni0, ne0:ne0, $
        R:R, Z:Z }


  return, a
end
