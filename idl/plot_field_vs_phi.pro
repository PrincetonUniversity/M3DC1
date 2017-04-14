pro plot_field_vs_phi, field, _EXTRA=extra, cutz=cutz, rrange=rrange, $
                       tpts=tpts, mesh=mesh, phirange=phirange

  if(n_elements(phirange) eq 0) then phirange = [0, 360.]

  ytitle = '!9P!6 (deg)!X'
  xtitle = '!8R!6 (m)!X'

  complex = read_parameter('complex',_EXTRA=extra)
  help, complex
  if(complex eq 1) then begin
     ntor = read_parameter('ntor',_EXTRA=extra)
     field = read_field(field,x,z,t,symbol=symbol,units=units,$
                        /complex,_EXTRA=extra)
     if(n_elements(tpts) eq 0) then tpts = 200
  endif else begin
     ff = read_field(field,x,z,t,symbol=symbol,units=units,$
                     _EXTRA=extra)
     if(n_elements(tpts) eq 0) then tpts = 20
  endelse
  n = n_elements(x)

  if(n_elements(cutz) eq 0) then begin
     cutz = 0.
  end
  if(n_elements(rrange) eq 0) then begin
     rrange = [min(x), max(x)]
  end

  r0 = rrange[0]
  r1 = rrange[1]

  z0 = cutz
  z1 = cutz

  rr = (r1 - r0)*findgen(n) / (n-1) + r0
  zz = (z1 - z0)*findgen(n) / (n-1) + z0
  phi = (phirange[1]-phirange[0])*findgen(tpts)/(tpts-1) + phirange[0]

  knx = fltarr(1,n,tpts)
  for j=0, tpts-1 do begin
     if(complex eq 1) then begin
        ff = real_part(field*exp(complex(0,ntor)*phi[j]*!pi/180.))
     endif else begin
        print, 'PHI = ', phi[j]
        ff = read_field(field,x,z,t,symbol=symbol,units=units,$
                        phi=phi[j],_EXTRA=extra)
     endelse
     kni = field_at_point(ff,x,z,rr,zz)

     knx[0,*,j] = kni[*]
  end
 
  label = symbol + ' ('+units +') at !8Z!6 = '+string(format='(g0)',cutz)+ '!X'
  contour_and_legend, knx, rr, phi, xtitle=xtitle, ytitle=ytitle, $
                      title=title, label=label
;  oplot, [rsep, rsep], !y.crange, linestyle=2

  if(keyword_set(mesh)) then begin
     ct3
     y = get_mesh_planes(_EXTRA=extra)
     for i=0, n_elements(y)-1 do begin
        oplot, !x.crange, [y[i], y[i]]*180./!pi, color=color(9,16)
     end
  end

end
