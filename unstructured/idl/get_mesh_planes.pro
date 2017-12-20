function get_mesh_planes, mesh, filename=filename, slice=slice
  if(n_params() eq 0) then begin
     mesh = read_mesh(filename=filename, slice=slice)
  end

  if(n_elements(mesh.elements._data[*,0]) lt 9) then begin
     print, n_elements(mesh.elements._data[*,0])
     return, 0.
  endif

  nelms = mesh.nelms._data 

  phi = float(mesh.elements._data[8,0])
  nplanes = 1

  for i=long(1), nelms-1 do begin
     phi0 = mesh.elements._data[8,i]
     newplane = 1
     for j=0, nplanes-1 do begin
        if(phi0 eq phi[j]) then begin
           newplane = 0
           break
        end
     end
     if(newplane eq 1) then begin
        phi = [phi, phi0]
        nplanes = nplanes+1
     end
  end
  
  return, phi[sort(phi)]
end
