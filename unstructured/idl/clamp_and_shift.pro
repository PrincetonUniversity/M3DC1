function clamp_and_shift, vec, shift=n
  ; clamp
  new = vec
  n = 0
  for i=0, n_elements(vec)-1 do begin
      if(new[i] lt -!pi) then new[i] = new[i] + 2.*!pi
      if(new[i] ge !pi) then new[i] = new[i] - 2.*!pi
      if(i gt 0) then begin
          if(abs(last - new[i]) gt !pi) then n = i
      endif
      last = new[i]
  end
  new = shift(new, -n)
  return, new
end
