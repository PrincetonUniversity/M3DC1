pro plot_kspits, filename=filename, yrange=yrange
   if(n_elements(filename) eq 0) then filename = 'C1.h5'
   if(hdf5_file_test(filename) eq 0) then return

   file_id = h5f_open(filename)
   root_id = h5g_open(file_id, "/")
   data = h5_parse(root_id, "kspits", /read_data)
   h5g_close, root_id
   h5f_close, file_id
   kspits = data.KSPITS._DATA

   if(n_elements(yrange) eq 0) then begin
   kspits_minmax = minmax(kspits)
   yrange=[kspits_minmax[0], kspits_minmax[1]]
   end
   print, 'plot range = ', yrange[0], yrange[1]

   dimn = size(kspits, /dim)
   print, 'total number of linear solvers and timesteps = ', dimn

   if(n_elements(maxn) eq 0) then begin
   maxn = dimn[0]
   end
   ntimes = dimn[1]
   print, 'max number of linear solvers to be plotted = ', maxn
   
   x = fltarr(ntimes)
   tmp = fltarr(ntimes)
   
   for n=0, maxn-1 do begin
      for t=0, ntimes-1 do begin
         ind = n + t*dimn[0]
         tmp[t] = kspits[ind]
         x[t] = t
      endfor
      if(n lt 1) then begin
         plot, x, tmp, yrange=yrange, TITLE='KSPSolve iteration numbers for 5, 1, 17, 6', linestyle=0
      endif else begin
         oplot, x, tmp, linestyle=0
      endelse

      if(n eq 0) then begin
      isolver=5
      end
      if(n eq 1) then begin
      isolver=1
      end
      if(n eq 2) then begin
      isolver=17
      end
      if(n eq 3) then begin
      isolver=6
      end
      numberAsString = STRTRIM(isolver, 2)
      xyouts, x[ntimes/2], tmp[ntimes/2]+0.25, numberAsString
   endfor

end


