;=========================================================
; read_quad_points
; ~~~~~~~~~~~~~~~~
;
; Returns the mesh data structure at a given time slice
;=========================================================
pro read_quad_points, filename=filename, slice=t, outfile=outfile

  if(n_elements(outfile) eq 0) then outfile='quad_pts.dat'
  
  mesh = read_mesh(filename=filename, slice=t)

  nelms = mesh.nelms._data
  npts = n_elements(mesh.quad_phi._data)

  sz = size(mesh.quad_r._data)
  nelms = sz[2]
  npts = sz[1]
  
  print, 'nelms = ', nelms
  print, 'npts = ', npts
  
  openw, ifile, outfile, /get_lun
  for i=0, nelms-1 do begin
     for j=0, npts-1 do begin
        printf, ifile, mesh.quad_r._data[j, i], mesh.quad_z._data[j, i], mesh.quad_phi._data[j, i]
     end
  end

  free_lun, ifile
end
