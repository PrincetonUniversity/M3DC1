;===========================================================================
; FIELD_SPECTRUM
; ~~~~~~~~~~~~~~
; Calculates the poloidal Fourier components of field
;
; Input:
;  field:         field to be decomposed
;  x, z:          x and z coordinates of field indices
;  fc (optional): flux coordinate structure
;
; Output:
;  spectral components of field 
;  fc:  flux coordinate structure
;  m:   poloidal mode numbers
;  n:   toroidal mode numbers
;===========================================================================
function field_spectrum, field, x, z, psi0=psi0, i0=i0, fc=fc, m=m, n=n, $
                         ignore_jacobian=nojac, $
                         dpsi0_dx=psi0_r, dpsi0_dz=psi0_z,_EXTRA=extra

;  b = flux_coord_field_new(field, x, z, fc=fc, /pest, _EXTRA=extra)

  a = flux_coord_field(field, psi0, x, z, i0=i0, fc=fc, $
                       dpsi0_dx=psi0_r, dpsi0_dz=psi0_z, _EXTRA=extra)
  b = transpose(a,[0,2,1])

  sz = size(field)
  nn = sz[1]

  for i=0, nn-1 do begin
     b[i,*,*] = b[i,*,*]*(2.*!pi)^2 * exp(-complex(0,1)*fc.omega)
     ; this was previously complex(0, ntor), but not sure that's right...
  end
  
  if(keyword_set(ignore_jacobian)) then begin
     print, 'field_spectrum: Ignoring Jacobian'
  endif else begin
     print, 'field_spectrum: Including Jacobian'
     for i=0, nn-1 do begin
        b[i,*,*] = b[i,*,*]*fc.j
     end
  end

  c = fft(b, -1, dimension=2)

  ; shift frequency values so that most negative frequency comes first
  nm = fc.m
  f = indgen(nm)
  f[nm/2+1] = nm/2 + 1 - nm + findgen((nm-1)/2)
  m = shift(f,-(nm/2+1))
  d = shift(c,0,-(nm/2+1),0)

  ; Flip m because we want a ~ exp(i n phi - i m theta)
  m = -m
  d = conj(d)

  ; do toroidal fourier transform
  if(nn gt 1) then begin
     c = fft(d, -1, dimension=1)
;     f = indgen(nn)
;     f[nn/2+1] = nn/2 + 1 - nn + findgen((nn-1)/2)
;     n = -shift(f,-(nn/2+1))
;     d = shift(c,-(nn/2+1),0,0)
     nt = nn/2+1
     n = indgen(nt)
;     print, 'nyquist = ', n[nt-1]
     ; combine positive and negative toroidal mode numbers
     for i=1, (nn-1)/2 do begin
;        print, 'adding ', -f[i], -f[nn-i]
        d[i,*,*] = c[i,*,*] + c[nn-i,*,*]
     end
     d = d[0:nt-1,*,*]

  endif else begin
     n = read_parameter('ntor', _EXTRA=extra)
  end

  return, d
end
