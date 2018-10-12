pro interpolate_boundary_points, infile, d=d, smooth=sm

  in = read_ascii(infile, data=1)
  rl = reform(in.field1[0,*])
  zl = reform(in.field1[1,*])
  if(n_elements(d) eq 0) then d = 0.1

  plot, rl, zl, psym=4

  n = n_elements(rl)

  for i=0, n-1 do begin
     if(i eq n-1) then begin
        rn = rl[0]
        zn = zl[0]
     endif else begin
        rn = rl[i+1]
        zn = zl[i+1]
     end
     dx = rn - rl[i]
     dz = zn - zl[i]
     dl = sqrt(dx^2 + dz^2)
     parts = floor(dl / d) + 1
     dd = dl / parts
     for m=0, parts-1 do begin
        r0 = rl[i] + m*dd*dx/dl
        z0 = zl[i] + m*dd*dz/dl
        if(m eq 0 and i eq 0) then begin
           r = r0
           z = z0
        endif else begin
           r = [r, r0]
           z = [z, z0]
        end
     end
  end

  p = n_elements(r)

  if(n_elements(sm) ne 0) then begin
     if(sm gt 0) then begin
        s =2*sm
        r = [r[p-s:p-1],r,r[0:s-1]]
        z = [z[p-s:p-1],z,z[0:s-1]]
        r = smooth(r, s)
        z = smooth(z, s)
        r = r[s:p-1+s]
        z = z[s:p-1+s]
     end
  end

  openw, ifile, 'boundary.dat', /get_lun
  printf, ifile, p
  for i=0, p-1 do printf, ifile, r[i], z[i]
  free_lun, ifile

  ct3
;  plot, r0, z0, /iso
;  oplot, r0, z0, psym=2
  oplot, r, z, color=color(1), linestyle=2
  oplot, r, z, color=color(1), linestyle=2, psym=2
end
