;=====================================================================
; dimensions
; ~~~~~~~~~~
;
; returs a vector with specified dimensions
;=====================================================================
function dimensions, energy=ener, eta=eta, j0=j, $
                     p0=pres, temperature=temp, potential=pot, $
                     l0=len, light=c, fourpi=fourpi, n0=den, $
                     v0=vel, t0=time, b0=mag, mu0=visc
  d = intarr(11)

  fp =    [1,0,0,0,0,0,0,0,0,0,0]
  c0 =    [0,1,0,0,0,0,0,0,0,0,0]
  n0 =    [0,0,1,0,0,0,0,0,0,0,0]
  v0 =    [0,0,0,1,0,0,0,0,0,0,0]
  b0 =    [0,0,0,0,1,0,0,0,0,0,0]
  temp0 = [0,0,0,0,0,1,0,0,0,0,0]
  i0 =    [0,0,0,0,0,0,1,0,0,0,0]
  e0 =    [0,0,0,0,0,0,0,1,0,0,0]
  t0 =    [0,0,0,0,0,0,0,0,1,0,0]
  l0 =    [0,0,0,0,0,0,0,0,0,1,0]
  pot0 =  [0,0,0,0,0,0,0,0,0,0,1]
  p0 = e0 - 3*l0
  
  if(keyword_set(fourpi)) then d = d + fp*fourpi
  if(keyword_set(light))  then d = d + c0*light
  if(keyword_set(den))    then d = d + n0*den
  if(keyword_set(vel))    then d = d + v0*vel
  if(keyword_set(mag))    then d = d + b0*mag
  if(keyword_set(time))   then d = d + t0*time
  if(keyword_set(len))    then d = d + l0*len
  if(keyword_set(temp))   then d = d + temp0*temp
  if(keyword_set(j))      then d = d + i0*j - 2*l0
  if(keyword_set(ener))   then d = d + ener*e0

  if(keyword_set(pot))  then d = d + pot*pot0
  if(keyword_set(eta))  then d = d +  eta*(2*l0-t0+fp-2*c0)
  if(keyword_set(pres)) then d = d + pres*(p0)
  if(keyword_set(visc)) then d = d + visc*(p0+t0)

  return, d
end
