pro test_mesh, filename, nplanes=nplanes, _EXTRA=extra
  mesh = read_mesh(filename=filename, _EXTRA=extra)
  nelms = mesh.nelms._data

  if(n_elements(nplanes) eq 0) then begin
     nplanes = read_parameter('nplanes', filename=filename)
     print, 'read nplanes.', nplanes
  end
  n = nelms/nplanes

  print, 'elms = ', nelms
  print, 'nplanes = ', nplanes
  print, 'elms_per_plane = ', n

  tol = 1e-6

  correct = long(0)
  wrong = lonarr(6)
  for i=long(0), n-1 do begin
     for k=1, nplanes-1 do begin
        for j=0, 5 do begin
           if(abs(mesh.elements._data[j,i]- $
                  mesh.elements._data[j,i+n*k]) gt tol) then begin
              wrong[j] = wrong[j] + 1
              print, 'misalignment of element ', i, ' at plane ', k
           endif else correct = correct +  1
        end
     end
  end
  for j=0, 5 do begin
     print, 'wrong ', j, ' = ', wrong[j]
  end
  print, 'correct = ', correct/6
end
