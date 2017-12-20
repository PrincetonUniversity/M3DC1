pro plot_field, name, time, mpeg=mpeg, _EXTRA=extra

  if(n_elements(mpeg) eq 0) then mpeg = 'movie.mpeg'

  ; open mpeg object
  mpegid = mpeg_open([640,480],quality=100)

  for i=time[0],time[1] do begin
     plot_field, name, i, _EXTRA=extra
     image = tvrd(true=1)
               
     image[0,*,*] = rotate(reform(image[0,*,*]), 7)
     image[1,*,*] = rotate(reform(image[1,*,*]), 7)
     image[2,*,*] = rotate(reform(image[2,*,*]), 7)
     
     mpeg_put, mpegid, image=image, frame=(i-time[0])
  end

  print, 'Writing mpeg...'
  mpeg_save, mpegid, filename=mpeg
  mpeg_close, mpegid
  return
end
