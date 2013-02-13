function minmax,array,NAN=nan, DIMEN=dimen, $
	SUBSCRIPT_MAX = subscript_max, SUBSCRIPT_MIN = subscript_min
;+
; NAME:
;      MINMAX
; PURPOSE:
;      Return a 2 element array giving the minimum and maximum of an array
; EXPLANATION:
;      Using MINMAX() is faster than doing a separate MAX and MIN.
;
; CALLING SEQUENCE:
;      value = minmax( array )
; INPUTS:
;      array - an IDL numeric scalar, vector or array.
;
; OUTPUTS:
;      value = a two element vector (if DIMEN is not supplied)
;            value[0] = minimum value of array
;            value[1] = maximum value of array
;
;            If the DIMEN keyword is supplied then value will be a 2 x N element
;            array where N is the number of elements in the specified
;            dimension
;              
; OPTIONAL INPUT KEYWORDS:
;      /NAN   - Set this keyword to cause the routine to check for occurrences
;            of the IEEE floating-point value NaN in the input data.  Elements 
;            with the value NaN are treated as missing data.
;
;      DIMEN - (V5.5 or later) integer (either 1 or 2) specifying which 
;            dimension of a 2-d array to  take the minimum and maximum.   Note
;            that DIMEN is only valid for a 2-d array, larger dimensions are 
;            not supported.
;
; OPTIONAL OUTPUT KEYWORDS:
;      SUBSCRIPT_MAX and SUBSCRIPT_MIN  Set either of these keywords to 
;            named variables to return the subscripts of the MIN and MAX
;	     values (V5.5 or later).
; EXAMPLE:
;     (1)  Print the minimum and maximum of an image array, im
; 
;            IDL> print, minmax( im )
;
;     (2) Given a 2-dimension array of (echelle) wavelengths w, print the
;         minimum and maximum of each order (requires V5.5 or later)
;
;         print,minmax(w,dimen=1)
;
; PROCEDURE:
;      The MIN function is used with the MAX keyword
;
; REVISION HISTORY:
;      Written W. Landsman                January, 1990
;      Converted to IDL V5.0   W. Landsman   September 1997
;      Added NaN keyword.      M. Buie       June 1998
;      Added DIMENSION keyword    W. Landsman  January 2002
;      Added SUBSCRIPT_MIN and SUBSCRIPT_MAX  BT Jan 2005
;      Check for IDL 5.5 or later, W. Thompson, 24-Feb-2005
;-
 On_error,2
 if !version.release ge '5.5' then begin
   if N_elements(DIMEN) GT 0 then begin
      amin = min(array, subscript_min, $
      	MAX = amax, NAN = nan, DIMEN = dimen, SUBSCRIPT_MAX = subscript_max) 
      return, transpose( [[amin], [amax] ])
   endif else  begin 
     amin = min( array, subscript_min, $
     	MAX = amax, NAN=nan, SUBSCRIPT_MAX = subscript_max)
     return, [ amin, amax ]
   endelse
 end else begin
   amin = min( array, MAX = amax, NAN=nan)
   return, [ amin, amax ]
 endelse
 end

pro plot_kehmn, filename=filename, xrange=xrange, yrange=yrange, maxn=maxn, ylog=ylog, growth=growth, outfile=outfile
   if(n_elements(filename) eq 0) then filename = 'C1.h5'
   if(hdf5_file_test(filename) eq 0) then return

   ; read harmonics [N, ntimes]
   file_id = h5f_open(filename)
   root_id = h5g_open(file_id, "/")
   data = h5_parse(root_id, "keharmonics", /read_data)
   h5g_close, root_id
   h5f_close, file_id
   kehmn = data.KEHARMONICS._DATA
   dimn = size(kehmn, /dim)
   print, 'total number of Fourier harmonics and timesteps = ', dimn

   ; read times, timestep [ntimes]
   s = read_scalars(filename=filename)
   time = s.time._data
   dt = s.dt._data

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
   if(n_elements(maxn) eq 0) then begin
   maxn = dimn[0]
   end
   ntimes = dimn[1]
   print, 'max number of Fourier harmonics to be plotted = ', maxn, ntimes
   
   ; if growth rate to be plotted
   grate=fltarr( maxn , (ntimes-1) )
   for n=0, maxn-1 do begin
      for t=0, ntimes-2 do begin
         ind = n + t*dimn[0]
         ind1 = n + (t+1)*dimn[0]
         ;print, n, t, ind, ind1, kehmn[ind] , kehmn[ind1]
         grate[n,t] = 2. / (kehmn[ind1] + kehmn[ind]) * (kehmn[ind1] - kehmn[ind]) / dt[t+1]
      endfor
   endfor

   ; get plot's yrange, default to minmax(...)
   if(n_elements(yrange) eq 0) then begin
      if(keyword_set(growth)) then begin
         grate_minmax = minmax(grate)
         yrange=[grate_minmax[0], grate_minmax[1]]
      endif else begin
         kehmn_minmax = minmax(kehmn)
         yrange=[kehmn_minmax[0], kehmn_minmax[1]]
      endelse
   end
   print, 'plot yrange = ', yrange[0], yrange[1]
   ; get plot's xrange, default to [time[1],time[ntimes]]
   if(n_elements(xrange) eq 0) then begin
         xrange=[time[1], time[dimn[1]-1]]
   end
   print, 'plot xrange = ', xrange[0], xrange[1]

   ; plot range 1:ntimes
   x = fltarr(ntimes-1)
   tmp = fltarr(ntimes-1)

   for n=0, maxn-1 do begin
      for t=0, ntimes-2 do begin
         if(keyword_set(growth)) then begin
            ind = n + t*dimn[0]
            tmp[t] = grate[ind]
            title='growth rate for each harmonics'
         endif else begin
            ind = n + (t+1)*dimn[0]
            tmp[t] = kehmn[ind]
            title='kinetic energy for each harmonics'
         endelse
            x[t] = time[t+1]
      endfor

      if(n lt 1) then begin
         if(keyword_set(ylog)) then begin
         plot, x, tmp, xrange=xrange, yrange=yrange, /ylog, TITLE=title, linestyle=0
         endif else begin
         plot, x, tmp, xrange=xrange, yrange=yrange, TITLE=title, linestyle=0
         endelse
      endif else begin
         oplot, x, tmp, linestyle=0
      endelse

      numberAsString = STRTRIM(n, 2)
      xyouts, x[ntimes/2], tmp[ntimes/2], numberAsString
   endfor

end

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


