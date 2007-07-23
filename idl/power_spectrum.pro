;==================================================================
; power_spectrum
; ~~~~~~~~~~~~~~
;
; Returns the power spectrum of signal "f".
; The frequency of each returned element is stored in "frequency".
; The total time interval may be specified with "t"
;==================================================================
function power_spectrum, f, frequency=frequency, t=t
  
  fft = fft(f)
  n = n_elements(fft)

  phi_old = abs(fft)^2
  phi = fltarr(n)
  frequency = fltarr(n)

  if(n_elements(t) eq 0) then t = n

  mid = n/2
  right = fix(n/2.-0.5)
  left = n-right-1

  phi[mid] = phi_old[0]
  frequency[mid] = 0.

  for i=1, left do begin
      phi[left-i] = phi_old[n-i]
      frequency[left-i] = -2.*3.14159*i/float(t)
  endfor
  for i=1, right do begin
      phi[mid+i] = phi_old[i]
      frequency[mid+i] = 2.*3.14159*i/float(t)
  endfor

  return, phi
end
