pro plot_field_vs_phi, field, _EXTRA=extra, cutx=cutx, cutz=cutz, $
                       tpts=tpts, mesh=mesh, phirange=phirange, $
                       slice=slice


  itor = read_parameter('itor',_EXTRA=extra)
  if(n_elements(phirange) eq 0) then begin
     if(itor eq 1) then begin
        phirange = [0, 360.]
     endif else begin
        rzero = read_parameter('rzero',_EXTRA=extra)
        phirange = [0, 2.*!pi*rzero]
     end
  end

  if(itor eq 1) then begin
     ytitle = '!9P!6 (deg)!X'
  endif else begin
     ytitle = '!8y!6 (m)!X'
  end

  complex = read_parameter('icomplex',_EXTRA=extra)

  if(slice eq -1) then begin
     ntor = 0.
     field = read_field(field,x,z,t,symbol=symbol,units=units,$
                        slice=-1,_EXTRA=extra)
     if(n_elements(tpts) eq 0) then tpts = 200
  endif else if(complex eq 1) then begin
     ntor = read_parameter('ntor',_EXTRA=extra)
     field = read_field(field,x,z,t,symbol=symbol,units=units,$
                        /complex,slice=slice,_EXTRA=extra)
     if(n_elements(tpts) eq 0) then tpts = 200
  endif else begin
     ff = read_field(field,x,z,t,symbol=symbol,units=units,$
                     _EXTRA=extra, slice=slice)
     if(n_elements(tpts) eq 0) then tpts = 20
  endelse
  n = n_elements(x)

  rrange = [min(x), max(x)]
  zrange = [min(z), max(z)]
  if(n_elements(cutz) ne 0) then begin
     r0 = rrange[0]
     r1 = rrange[1]

     z0 = cutz
     z1 = cutz
     cutv = cutz
     cut = 'Z'
     xtitle = '!8R!6 (m)!X'
  endif else if(n_elements(cutx) ne 0) then begin
     r0 = cutx
     r1 = cutx

     z0 = zrange[0]
     z1 = zrange[1]
     cutv = cutx
     cut = 'R'
     xtitle = '!8Z!6 (m)!X'
  end


  rr = (r1 - r0)*findgen(n) / (n-1) + r0
  zz = (z1 - z0)*findgen(n) / (n-1) + z0
  phi = (phirange[1]-phirange[0])*findgen(tpts)/(tpts-1) + phirange[0]

  knx = fltarr(1,n,tpts)
  for j=0, tpts-1 do begin
     if(slice eq -1) then begin
        ff = reform(field)
     endif else if(complex eq 1) then begin
        ff = real_part(field*exp(complex(0,ntor)*phi[j]*!pi/180.))
     endif else begin
        print, 'PHI = ', phi[j]
        ff = read_field(field,x,z,t,symbol=symbol,units=units,$
                        phi=phi[j],slice=slice,_EXTRA=extra)
     endelse
     kni = field_at_point(ff,x,z,rr,zz)

     knx[0,*,j] = kni[*]
  end

  if(strlen(strtrim(units)) ne 0) then begin
     u = ' (' + units +')'
  endif else begin
     u = ''
  end

  if(n_elements(cutz) ne 0) then begin
     v = rr
  endif else if(n_elements(cutx) ne 0) then begin
     v = zz
  end

  
  label = symbol + u + ' at !8'+cut+'!6 = '+string(format='(g0)',cutv)+ '!X'
  contour_and_legend, knx, v, phi, xtitle=xtitle, ytitle=ytitle, $
                      title=title, label=label, _EXTRA=extra
;  oplot, [rsep, rsep], !y.crange, linestyle=2

  if(keyword_set(mesh)) then begin
     ct3
     y = get_mesh_planes(_EXTRA=extra)
     for i=0, n_elements(y)-1 do begin
        oplot, !x.crange, [y[i], y[i]]*180./!pi, color=color(9,16)
     end
  end
end
