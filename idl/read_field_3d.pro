function read_field_3d, name, phi, x, z, t, points=points, tpoints=tpoints, $
                        symbol=symbol, units=u, ntor=ntor, _EXTRA=extra

  if(n_elements(points) eq 0) then points=200
  if(n_elements(tpoints) eq 0) then tpoints=16

  field = fltarr(tpoints, points, points)

  phi = fltarr(tpoints)

  for i=0, tpoints-2 do begin
     phi[i] = 360.*i/(tpoints-1)
     field[i,*,*] = reform(read_field(name,x,z,t,phi=phi[i],points=points,$
                                      symbol=symbol, units=u, _EXTRA=extra))
  end
  phi[tpoints-1] = 360.
  field[tpoints-1,*,*] = field[0,*,*]

  if(n_elements(ntor) eq 0) then begin
     return, field
  endif else begin
     fftfield = complexarr(tpoints,points,points)
     for i=0, points-1 do begin
        for j=0, points-1 do begin
           fftfield[*,i,j] = fft(field[*,i,j])
        end
     end
     lastfield = complexarr(1,points,points)
     lastfield = fftfield[ntor, *, *]
     return, lastfield
  end
end
