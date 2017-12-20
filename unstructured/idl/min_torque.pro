; find the relative phase and amplitude of file2
; that minimizes quasilinear torque
; with error field file1

function get_torque, br1, br2, bz1, bz2, $
                     jr1, jr2, jz1, jz2, $
                     fac, r, z, t
   n = n_elements(fac)

   torque = fltarr(n, n)
   area = (max(r) - min(r))*(max(z) - min(z))
   rr = radius_matrix(r,z,t)

   for i=0, n-1 do begin
       for j=0, n-1 do begin
           f = complex(fac[i],fac[j])
           br = br1 + f*br2
           bz = bz1 + f*bz2
           jr = jr1 + f*jr2
           jz = jz1 + f*jz2
           tphi = (conj(jz)*br + jz*conj(br)) $
             -    (conj(jr)*bz + jr*conj(bz))
           torque[i, j] = mean(rr*real_part(tphi))*area
       end
   end

   return, torque
end

function get_modb, br1, br2, bz1, bz2, bphi1, bphi2, $
                   fac, r, z, t
   n = n_elements(fac)

   modb = fltarr(n,n)
   area = (max(r) - min(r))*(max(z) - min(z))
   rr = radius_matrix(r,z,t)

   for i=0, n-1 do begin
       for j=0, n-1 do begin
           f = complex(fac[i],fac[j])
           br = br1 + f*br2
           bz = bz1 + f*bz2
           bphi = bphi1 + f*bphi2
           b2 = conj(br)*br + conj(bz)*bz + conj(bphi)*bphi
           modb[i, j] = mean(rr*real_part(b2))*area
       end
   end

   return, modb
end

function get_overlap, file1, file2, fac

   n = n_elements(fac)
   overlap = fltarr(n,n)

   for i=0, n-1 do begin
       print, i
       for j=0, n-1 do begin
           cur = [complex(1,0), complex(fac[i],fac[j])]
           overlap[i, j] = island_overlap([file1,file2],cur=cur,/sum_files)
       end
   end

   return, overlap
end

function get_sigma, file1, file2, fac, psi0

   n = n_elements(fac)
   overlap = fltarr(n,n)

   for i=0, n-1 do begin
       print, i
       for j=0, n-1 do begin
           cur = [complex(1,0), complex(fac[i],fac[j])]
           p = psi0
           overlap[i, j] = chi95([file1,file2],cur=cur,/sum_files,psi0=p)
       end
   end

   return, overlap
end

function get_bmn, bmn1, bmn2, psi, m, psi0, m0, fac
   n = n_elements(fac)
   bmn = fltarr(n,n)

   ipsi = interpol(findgen(n_elements(psi)), psi, psi0)
   if(ipsi lt 0 or ipsi ge n_elements(psi)) then begin
       print, 'Psi ', psi0, 'not found'
       return, bmn
   end
   print, psi0, ipsi, psi[fix(ipsi)]

   im = where(m eq m0, count)
   if(n_elements(count) ne 1) then begin
       print, 'Mode ', m0, ' not found'
       return, bmn
   end
   b1 = interpolate(bmn1, im, ipsi)
   b2 = interpolate(bmn2, im, ipsi)

   for i=0, n-1 do begin
       for j=0, n-1 do begin
           f = complex(fac[i],fac[j])
           bmn[i, j] = abs(b1 + f*b2)
       end
   end

   return, bmn
end

function get_bmntot, bmn1, bmn2, psi, m, psi0, m0, fac
   n = n_elements(fac)
   bmn = fltarr(n,n)

   for i=0, n-1 do begin
       for j=0, n-1 do begin
           f = complex(fac[i],fac[j])
           bmn[i, j] = total(abs(bmn1 + f*bmn2))
       end
   end

   return, bmn
end

pro min_torque, file1, file2, modb=modb, maxamp=maxamp, overlap=overlap, $
                sigma=sigma, mpol=mpol, psi0=psi0, q0=q0, bmntot=bmntot, $
                scale=scale, _EXTRA=extra

   if(n_elements(maxamp) eq 0) then maxamp=2.
   if(n_elements(psi0) eq 0) then psi0=0.95
   if(n_elements(scale) eq 0) then scale = 1.
   if(n_elements(scale) eq 1) then scale = [scale, scale]

   if(not keyword_set(overlap) and $
      not keyword_set(sigma) and $
      not keyword_set(bmntot) and $
      not keyword_set(mpol)) then begin
       ntor1 = read_parameter('ntor', file=file1)
       ntor2 = read_parameter('ntor', file=file2)
       if(ntor1 eq ntor2) then begin
           ntor = ntor1
       endif else begin
           print, 'Error: files have different ntor'
           return
       end
       
       extsub1 = read_parameter('extsubtract', file=file1)
       extsub2 = read_parameter('extsubtract', file=file1)
       if(extsub1 eq 1) then suf1 = '_plasma' else suf1 = ''
       if(extsub2 eq 1) then suf2 = '_plasma' else suf2 = ''
   end

   if(keyword_set(overlap)) then begin
   endif else if(keyword_set(sigma)) then begin
   endif else if(keyword_set(mpol) or keyword_set(bmntot)) then begin
       read_bmncdf, file=file1, bmn=bmn1, psi=psi1, m=m1, q=q1
       read_bmncdf, file=file2, bmn=bmn2, psi=psi2, m=m2, q=q2
       bmn1 = bmn1*scale[0]
       bmn2 = bmn2*scale[1]

       if(keyword_set(q0)) then begin
           psi0 = interpol(psi1, q1, q0)
           print, 'psi0 = ', psi0
       end

   endif else if(keyword_set(modb)) then begin
       slice=1
       br1 = read_field('bx',x,z,t,file=file1,slice=slice,/linear, $
                        /complex,/mks,edge_val=0.,_EXTRA=extra)*scale[0]
       br2 = read_field('bx',x,z,t,file=file2,slice=slice,/linear, $
                        /complex,/mks,edge_val=0.,_EXTRA=extra)*scale[1]
       bz1 = read_field('bz',x,z,t,file=file1,slice=slice,/linear, $
                        /complex,/mks,edge_val=0.,_EXTRA=extra)*scale[0]
       bz2 = read_field('bz',x,z,t,file=file2,slice=slice,/linear, $
                        /complex,/mks,edge_val=0.,_EXTRA=extra)*scale[1]
       bphi1 = read_field('by',x,z,t,file=file1,slice=slice,/linear, $
                          /complex,/mks,edge_val=0.,_EXTRA=extra)*scale[0]
       bphi2 = read_field('by',x,z,t,file=file2,slice=slice,/linear, $
                          /complex,/mks,edge_val=0.,_EXTRA=extra)*scale[1]

   endif else begin
       ; read fields from coils
       coil_br1 = read_field('bx',x,z,t,file=file1,slice=0,/linear, $
                             /complex,/mks,edge_val=0.,_EXTRA=extra)*scale[0]
       coil_br2 = read_field('bx',x,z,t,file=file2,slice=0,/linear, $
                             /complex,/mks,edge_val=0.,_EXTRA=extra)*scale[1]
       coil_bz1 = read_field('bz',x,z,t,file=file1,slice=0,/linear, $
                             /complex,/mks,edge_val=0.,_EXTRA=extra)*scale[0]
       coil_bz2 = read_field('bz',x,z,t,file=file2,slice=0,/linear, $
                             /complex,/mks,edge_val=0.,_EXTRA=extra)*scale[1]
       
       ; read currents from plasma
       plasma_jr1 = read_field('jr'+suf1,x,z,t,file=file1,slice=1,/linear, $
                               /complex,/mks,edge_val=0.,_EXTRA=extra)*scale[0]
       plasma_jr2 = read_field('jr'+suf2,x,z,t,file=file2,slice=1,/linear, $
                               /complex,/mks,edge_val=0.,_EXTRA=extra)*scale[1]
       plasma_jz1 = read_field('jz'+suf1,x,z,t,file=file1,slice=1,/linear, $
                               /complex,/mks,edge_val=0.,_EXTRA=extra)*scale[0]
       plasma_jz2 = read_field('jz'+suf2,x,z,t,file=file2,slice=1,/linear, $
                               /complex,/mks,edge_val=0.,_EXTRA=extra)*scale[1]
   end

   namp = 100
   fac = 2.*maxamp*findgen(namp)/(namp - 1) - maxamp

   if(keyword_set(overlap)) then begin
       val = get_overlap(file1, file2, fac*scale[0])
       label = 'Island Overlap Width (!7W!X))'
   endif else if(keyword_set(sigma)) then begin
       val = get_sigma(file1, file2, fac*scale[0], psi0)
       label = string(format='("Local Chirikov at !7W!X = ",G0)', psi0)
   endif else if(keyword_set(mpol)) then begin
       val = get_bmn(bmn1, bmn2, psi1, m1, psi0, mpol, fac)
       label = 'Bmn'
   endif else if(keyword_set(bmntot)) then begin
       val = get_bmntot(bmn1, bmn2, psi1, m1, psi0, mpol, fac)
       label = 'Bmn Total'
   endif else if(keyword_set(modb)) then begin
       val = get_modb(br1, br2, bz1, bz2, bphi1, bphi2, $
                      fac,x,z,t)
       val = val/(2.*!pi*4.e-7)
       label = 'Perturbed Magnetic Energy (J)'
   endif else begin
       val = get_torque(coil_br1, coil_br2, coil_bz1, coil_bz2, $
                        plasma_jr1, plasma_jr2, plasma_jz1, plasma_jz2, $
                        fac,x,z,t)
       label = 'J!9x!XB Torque (N!9.!Xm)!X'
   end

   minval = min(val, /abs, i)
   j = i mod namp
   k = i / namp

   plt = findgen(1,namp,namp)
   plt[0,*,*] = val

   opt_amp = sqrt(fac[j]^2 + fac[k]^2)
   opt_ang = 360.*atan(fac[k],fac[j])/!pi
   subtitle = string(format='("Optimum: ",F0.2," kA at ",F0.1,"!9%!X")', $
                     opt_amp, opt_ang)

   if(1 eq strcmp('X', !d.name)) then begin
       window, 0
       contour_and_legend, plt, fac, fac, label=label, subtitle=subtitle, $
         /lines, _EXTRA=extra
       oplot, [0,0], !y.crange
       oplot, !x.crange, [0,0]
       ct3
       oplot, [fac[j], fac[j]], [fac[k], fac[k]], psym=7
   end

   print, minval, val[j,k]

   ; find all points where |fac|=1

   theta = 2.*!pi*findgen(namp)/namp
   xi = (cos(theta) - min(fac))/(max(fac) - min(fac))*n_elements(fac)
   zi = (sin(theta) - min(fac))/(max(fac) - min(fac))*n_elements(fac)
   if(1 eq strcmp('X', !d.name)) then $
     oplot, interpolate(fac,xi), interpolate(fac,zi)

   vint = interpolate(val, xi, zi)
   if(1 eq strcmp('X', !d.name)) then window, 1
   if(1 eq strcmp('X', !d.name)) then print, 'Using X'

   plot, 180.*theta/!pi, vint, xrange=[0, 360], xstyle=1, $
     xtitle='Phase (deg)', ytitle=label, _EXTRA=extra
end
