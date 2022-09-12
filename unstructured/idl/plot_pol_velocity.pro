; ==================================================
; plot_pol_velocity
; ~~~~~~~~~~~~~~~~~
;
; makes a vector plot of the poloidal velocity
; ==================================================
pro plot_pol_velocity, time,  maxval=maxval, points=points, $
                       lcfs=lcfs, _EXTRA=extra

  if(n_elements(pts) eq 0) then pts=25

  nv = read_parameter('numvar', _EXTRA=extra)
  itor = read_parameter('itor', _EXTRA=extra)
  if(n_elements(itor) gt 1) then itor=itor[0]

  phi = read_field('phi', x, z, t, points=points, _EXTRA=extra, slice=time)
  if(n_elements(phi) le 1) then return

  if(itor eq 1) then r = radius_matrix(x,z,t) else r = 1.

;....modified for new form of velocity
  vx = -dz(phi,z)*r
  vz =  dx(phi,x)*r

  chi = read_field('chi', x, z, t, points=points, _EXTRA=extra, slice=time)
  vx = vx + dx(chi,x)/r^2
  vz = vz + dz(chi,z)/r^2

  bigvel = max(sqrt(vx^2 + vz^2))
  print, "maximum velocity: ", bigvel
  if(n_elements(maxval) ne 0) then begin
      length = bigvel/maxval
  endif else length=1

  if(n_elements(title) eq 0) then begin
      title = '!6Poloidal Flow!X'
      if(t gt 0) then begin
          title = title +  $
            string(FORMAT='("!6(!8t!6 = ",G0," !7s!D!8A!60!N)!X")', t)
      endif else begin
          title = title + $
            string(FORMAT='("!6(!8t!6 = ",G0,")!X")', t)
      endelse
  endif
  
  maxstr=string(format='("!6max(!8u!Dpol!N!6) = ",G0.3,"!X")',bigvel) + $
    '!6 ' + make_units(/v0) + '!X'

  velovect, reform(vx), reform(vz), x, z, length=length, _EXTRA=extra, $
    xtitle=make_label('!8R!X', /l0, _EXTRA=extra), $
    ytitle=make_label('!8Z!X', /l0, _EXTRA=extra), $
    title=title, subtitle=maxstr

  if(keyword_set(lcfs)) then plot_lcfs, /over, points=200, slice=time, _EXTRA=extra
end
