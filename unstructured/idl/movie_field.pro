pro movie_field, name, maxtime, outfile=outfile,rate=rate, width=width, height=height, lcfs=lcfs, ext=ext, _EXTRA=ex

  if n_elements(width) eq 0 then width = 580
  if n_elements(height) eq 0 then height = 750
  
  ; Available formats: avi flv gif matroska mjpeg mov mp4 swf wav webm
  if n_elements(ext) eq 0 then ext='mp4'
  
  window, 0, xsize=width, ysize=height
  
  if n_elements(outfile) eq 0 then outfile=name+'.'+ext

  oVid = IDLffVideoWrite(outfile)
  if n_elements(rate) eq 0 then rate=5
  vidStream = oVid.AddVideoStream(width,height,rate)
  
  for i = 0,maxtime do begin
    plot_field,name,i,time=time,_EXTRA=ex
    if keyword_set(lcfs) then begin
      plot_lcfs, slice=0,_EXTRA=ex
    endif
    xyouts,0.57*width,0.91*height,'t = '+strtrim(string(time,format='(E10.3)'),2),charsize=1.5,color=0,/device
    frame = tvrd(true=1)
    !NULL = oVid.Put(vidStream,frame)
  endfor
  oVid = 0
end