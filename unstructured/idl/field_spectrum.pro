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
;===========================================================================
function field_spectrum, field, x, z, psi0=psi0, i0=i0, fc=fc, m=m, $
                         ignore_jacobian=nojac, _EXTRA=extra

;  b = flux_coord_field_new(field, x, z, fc=fc, /pest, _EXTRA=extra)

  a = flux_coord_field(field, psi0, x, z, i0=i0, fc=fc, _EXTRA=extra)
  b = transpose(a,[0,2,1])

  ntor = read_parameter('ntor', _EXTRA=extra)
  print, 'NTOR = ', ntor

  if(keyword_set(ignore_jacobian)) then begin
     print, 'field_spectrum: Ignoring Jacobian'
     b[0,*,*] = b[0,*,*]*(2.*!pi)^2 $
                * exp(complex(0,ntor)*fc.omega)
  endif else begin
     print, 'field_spectrum: Including Jacobian'
     b[0,*,*] = b[0,*,*]*(2.*!pi)^2*fc.j $
                * exp(complex(0,ntor)*fc.omega)
  end

  c = fft(b, -1, dimension=2)

  ; shift frequency values so that most negative frequency comes first
  n = fc.m
  f = indgen(n)
  f[n/2+1] = n/2 + 1 - n + findgen((n-1)/2)
  m = shift(f,-(n/2+1))
  d = shift(c,0,-(n/2+1),0)

  return, d
end
