pro plot_hmn, filename=filename,  maxn=maxn, growth=growth, outfile=outfile,$
                yrange=yrange, smooth=sm, overplot=over, _EXTRA=extra, $
              ke=ke, me=me, xscale=xscale, labelx=labelx
   if(n_elements(labelx) eq 0) then labelx = 0.5
   if(n_elements(filename) eq 0) then filename = 'C1.h5'
   if(hdf5_file_test(filename) eq 0) then return

   ; read harmonics [N, ntimes]
   file_id = h5f_open(filename)
   root_id = h5g_open(file_id, "/")
   if(keyword_set(me)) then begin
      data = h5_parse(root_id, "bharmonics", /read_data)
      name = 'Magnetic Energy'
      kehmn = data.BHARMONICS._DATA
   endif else begin
      data = h5_parse(root_id, "keharmonics", /read_data)
      name = 'Kinetic Energy'
      kehmn = data.KEHARMONICS._DATA
   end
   h5g_close, root_id
   h5f_close, file_id

   dimn = size(kehmn, /dim)
   print, 'total number of Fourier harmonics and timesteps = ', dimn

   ; read times, timestep [ntimes]
   time = read_scalar('time', filename=filename, units=u, _EXTRA=extra)
   xtitle = '!8t !6(' + u + ')!X'

   ; write harmonics [N, ntimes] into "outfile"
      if(keyword_set(outfile)) then begin
         ;format=string(39B)+'(' + STRTRIM(1+dimn[0], 2) + 'E16.6)'+string(39B)
         format='(' + STRTRIM(1+dimn[0], 2) + 'E16.6)'
         print, format
         openw,ifile,outfile,/get_lun
         printf,ifile,format=format,[transpose(time),kehmn]
         free_lun, ifile
      endif

   ; get the maximum number of fourier harmonics to be plotted, default to dim[0]
   if(n_elements(maxn) eq 0) then maxn = dimn[0]

   ntimes = dimn[1]
   if(ntimes gt n_elements(time)) then ntimes=n_elements(time)
   print, 'max number of Fourier harmonics to be plotted = ', maxn, ntimes
   ke = fltarr(maxn, ntimes)
   grate=fltarr(maxn ,ntimes)
   for n=0, maxn-1 do begin
      ke[n,*] = kehmn[n:n+(ntimes-1)*dimn[0]:dimn[0]]
      grate[n,*] = deriv(time, alog(ke[n,*]))      ;  /2. removed 12/12/16 scj
   endfor

   ; plot range 1:ntimes
   if(keyword_set(growth)) then begin
      tmp = grate
      ytitle='!6Growth Rate!X'
   endif else begin
      tmp = ke
      ytitle='!6' + name + '!X'
   endelse

   ; smooth data
   if(n_elements(sm) ne 0) then begin
      for n=0,maxn-1 do tmp[n,*] = smooth(tmp[n,*], sm)
   end

   ; get plot's yrange, default to minmax(...)
   if(n_elements(yrange) eq 0) then yrange=[min(tmp), max(tmp)]

   if(n_elements(xscale) eq 0) then xscale = 1.

   c = get_colors(n)
   for n=0, maxn-1 do begin

      if(n lt 1 and not keyword_set(over)) then begin
         plot, time[1:ntimes-1]*xscale, tmp[n,1:ntimes-1], $
               xtitle=xtitle, ytitle=ytitle, yrange=yrange, $
               _EXTRA=extra
      endif else begin
         oplot, time*xscale, tmp[n,*], linestyle=0, color=c[n]
      endelse

      numberAsString = STRTRIM(n, 2)
      m = min([ntimes*labelx, ntimes-1])
      xyouts, time[m]*xscale, tmp[n,m], numberAsString, color=c[n]
   endfor
   plot_legend, string(format='("!8n!6=",I0,"!X")',indgen(maxn)), $
                color=c, _EXTRA=extra

end
