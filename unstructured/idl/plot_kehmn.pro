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

pro plot_kehmn, filename=filename,  maxn=maxn, growth=growth, outfile=outfile,$
                yrange=yrange, smooth=sm, _EXTRA=extra
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
   print, 'max number of Fourier harmonics to be plotted = ', maxn, ntimes
   ke = fltarr(maxn, ntimes)
   grate=fltarr(maxn ,ntimes)
   for n=0, maxn-1 do begin
      ke[n,*] = kehmn[n:n+(ntimes-1)*dimn[0]:dimn[0]]
      grate[n,*] = deriv(time, alog(ke[n,*]))        ; /2. removed 12/12/16 scj
   endfor

   ; plot range 1:ntimes
   if(keyword_set(growth)) then begin
      tmp = grate
      ytitle='!6Growth Rate!X'
   endif else begin
      tmp = ke
      ytitle='!6Kinetic Energy!X'
   endelse

   ; smooth data
   if(n_elements(sm) ne 0) then begin
      for n=0,maxn-1 do tmp[n,*] = smooth(tmp[n,*], sm)
   end

   ; get plot's yrange, default to minmax(...)
   if(n_elements(yrange) eq 0) then begin
      mm = minmax(tmp[*,1:ntimes-1])
      yrange=[mm[0], mm[1]]
   end

   c = get_colors(n)
   for n=0, maxn-1 do begin

      if(n lt 1) then begin
         plot, time[1:ntimes-1], tmp[n,1:ntimes-1], xtitle=xtitle, ytitle=ytitle, yrange=yrange, $
               _EXTRA=extra
      endif else begin
         oplot, time, tmp[n,*], linestyle=0, color=c[n]
      endelse

      numberAsString = STRTRIM(n, 2)
      xyouts, time[ntimes/2], tmp[n,ntimes/2], numberAsString, color=c[n]
   endfor
   plot_legend, string(format='("!8n!6=",I0,"!X")',indgen(maxn)), $
                ylog=ylog, color=c

end

