pro movie_field, name, maxtime,dt=dt, outfile=outfile,rate=rate, width=width, height=height, $
                 lcfs=lcfs, ext=ext, Nphi=Nphi, plot_pellet=plot_pellet, fac=fac, $
                 maxvel=maxvel, thick=thick, charsize=charsize, points=points, phi=phi, _EXTRA=ex

  if n_elements(width) eq 0 then width = 580
  if n_elements(height) eq 0 then height = 750
  if n_elements(phi) eq 0 then phi = 0.
  
  ; Available formats: avi flv gif matroska mjpeg mov mp4 swf wav webm
  if n_elements(ext) eq 0 then ext='mp4'
  
  if keyword_set(plot_pellet) then begin
    pr = read_scalar('pellet_r',time=pt,_EXTRA=ex)
    pz = read_scalar('pellet_z',time=pt,_EXTRA=ex)
  endif
  
  window, 0, xsize=width, ysize=height
  
  if n_elements(outfile) eq 0 then outfile=name+'.'+ext

  oVid = IDLffVideoWrite(outfile)
  if n_elements(rate) eq 0 then rate=5
  vidStream = oVid.AddVideoStream(width,height,rate)
  
  if n_elements(Nphi) ne 0 then begin
    Nf = Nphi + 1
    phis = 360.*findgen(Nf)/(Nf-1)
    slices = 0.*indgen(Nf) + maxtime
  endif else begin
    Nf = maxtime+1
    slices = indgen(Nf)
    phis = 0.*findgen(Nf) + phi
  endelse
  
  
  for i = 0,Nf-1,dt do begin
    phi = phis[i]
    slice = slices[i]
    plot_field, name, slice, time=time, fac=fac, thick=thick, charsize=charsize, $
      points=points, phi=phi, lcfs=lcfs, _EXTRA=ex
    
    if keyword_set(plot_pellet) then begin
      r = interpol(pr,pt,time)
      z = interpol(pz,pt,time)
      oplot,[r,r],[z,z],psym=4,symsize=1,thick=5,color=color(7)
    endif
    if n_elements(maxvel) ne 0 then begin
      print,'maxvel'
      plot_pol_velocity,slice,points=40,/over, maxval=maxvel,phi=phi,_EXTRA=ex
    endif
    if n_elements(Nphi) ne 0 then begin
      xyouts,0.55*width,0.91*height,'phi = '+strtrim(string(phi),2),charsize=1.5,color=color(0),/device
    endif else begin
      xyouts,0.55*width,0.91*height,'t = '+strtrim(string(time,format='(E10.3)'),2),charsize=1.5,color=color(0),/device
    endelse
    frame = tvrd(true=1)
    !NULL = oVid.Put(vidStream,frame)
  endfor
  oVid = 0
end
