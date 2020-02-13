function read_lcfs, axis=axis, xpoint=xpoint, flux0=flux0, $
                    filename=filename, slice=time, last=last, $
                    mks=mks, cgs=cgs, _EXTRA=extra

   s = read_scalars(filename=filename)
   if(n_elements(time) eq 0) then begin
       slice = 0
   endif else slice = time
   if(keyword_set(last)) then begin
       slice = read_parameter("ntime", filename=filename)-1
   endif

   t0 = get_slice_time(filename=filename, slice=slice)

   tmp = s.time._data[*] - t0[0]
   dum = min(tmp, i, /abs)
   print, 'slice time = ', t0
   print, 'time step time: ', s.time._data[i]
   print, 'time slice: ', i
   
   xpoint = fltarr(2)
   axis = fltarr(2)
   xpoint[0] = s.xnull._data[i]
   xpoint[1] = s.znull._data[i]
   axis[0] = s.xmag._data[i]
   axis[1] = s.zmag._data[i]
   
   flux0 = s.psimin._data[i]
   psilim = s.psi_lcfs._data[i]

   if(keyword_set(mks) or keyword_set(cgs)) then begin
      itor = read_parameter("itor", filename=filename)
      get_normalizations, filename=filename,b0=b0,n0=n0,l0=l0,ion=mi
      convert_units, xpoint, dimensions(/l0), b0, n0, l0, mi, cgs=cgs, mks=mks
      convert_units, axis, dimensions(/l0), b0, n0, l0, mi, cgs=cgs, mks=mks
      convert_units, flux0, dimensions(l0=1+itor, /b0), b0, n0, l0, mi, $
                     cgs=cgs, mks=mks
      convert_units, psilim, dimensions(l0=1+itor, /b0), b0, n0, l0, mi, $
                     cgs=cgs, mks=mks
   end

   return, psilim
end
