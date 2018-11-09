pro read_nulls, axis=axis, xpoints=xpoint, _EXTRA=extra
   s = read_scalars(_EXTRA=extra)

   t0 = get_slice_time(_EXTRA=extra)

   dum = min(s.time._data - t0[0], i, /abs)

   xpoint = fltarr(2)
   axis = fltarr(2)
   xpoint[0] = s.xnull._data[i]
   xpoint[1] = s.znull._data[i]
   axis[0] = s.xmag._data[i]
   axis[1] = s.zmag._data[i]
   
   return
end

; ==============================================================
; find_nulls
; ----------
;
;  Finds field nulls.
;  xpoint = fltarr(2) : position of active x-point
;  axis   = fltarr(2) : position of magnetic axis
; ==============================================================
pro find_nulls, psi, x, z, axis=axis, xpoints=xpoint, _EXTRA=extra

   if(n_elements(time) eq 0) then time = 0
   
   field = s_bracket(psi,psi,x,z)
   d2 = dz(dz(psi,z),z)*dx(dx(psi,x),x)
   
   nulls = field lt mean(field)/10.

   sz = size(field)

   xpoints = 0
   xpoint = 0.
   axes = 0
   axis = 0.

   oldxflux = min(psi)
   oldaflux = min(psi)
   
   for i=0, sz[2]-1 do begin
       for j=0, sz[3]-1 do begin
           if(nulls[0,i,j] eq 0) then continue

           currentmin = field[0,i,j]
           currentpos = [i,j]
           foundlocalmin = 0

           ; find local minimum
           for m=i, sz[2]-1 do begin
               if(nulls[0,m,j] eq 0) then break
               for n=j+1, sz[3]-1 do begin
                   if(nulls[0,m,n] eq 0) then break

                   if(field[0,m,n] le currentmin) then begin
                       currentmin = field[0,m,n]
                       currentpos = [m,n]
                       foundlocalmin = 1
                   endif
                   nulls[0,m,n] = 0
               endfor
               for n=j-1, 0, -1 do begin
                   if(nulls[0,m,n] eq 0) then break

                   if(field[0,m,n] le currentmin) then begin
                       currentmin = field[0,m,n]
                       currentpos = [m,n]
                       foundlocalmin = 1
                   endif
                   nulls[0,m,n] = 0
               endfor
           endfor

           if (foundlocalmin eq 0) then continue

           ; throw out local minima on boundaries
           if (currentpos[0] eq 0) or (currentpos[0] eq sz[2]-1) then continue
           if (currentpos[1] eq 0) or (currentpos[1] eq sz[3]-1) then continue

           ; determine if point is an x-point or an axis and
           ; append the location index to the appropriate array
           if(d2[0,currentpos[0],currentpos[1]] lt 0) then begin
               if(psi[0,currentpos[0],currentpos[1]] gt oldxflux) then begin
                   xpoint = [x[currentpos[0]], z[currentpos[1]]]
                   oldxflux = psi[0,currentpos[0],currentpos[1]]
               end
           endif else begin
               if(psi[0,currentpos[0],currentpos[1]] gt oldxflux) then begin
                   axis = [x[currentpos[0]], z[currentpos[1]]]
                   oldaflux = psi[0,currentpos[0],currentpos[1]]
               endif
           endelse
       endfor
   endfor 

   if(n_elements(axis) lt 2) then begin
       print, 'Warning: cannot find magnetic axis'
   endif

   fieldr = dx(field,x)
   fieldz = dz(field,z)

   ; do iterative refinement on magnetic axis
   for k=1, 2 do begin
       pt = field_at_point(field, x, z, axis[0], axis[1])
       pt1 = field_at_point(fieldr, x, z, axis[0], axis[1])
       pt2 = field_at_point(fieldz, x, z, axis[0], axis[1])

       denom = pt1^2 + pt2^2           
       dx = pt*pt1/denom
       dz = pt*pt2/denom           

       axis[0] = axis[0] - dx
       axis[1] = axis[1] - dz
   end

   print, 'Found axis at ', axis[0], axis[1]


   if(n_elements(xpoint) ge 2) then begin
       print, 'x-point guess at ', xpoint[0], xpoint[1]
       ; do iterative refinement on x-point
       for k=1, 10 do begin
           pt = field_at_point(field, x, z, xpoint[0], xpoint[1])
           pt1 = field_at_point(fieldr, x, z, xpoint[0], xpoint[1])
           pt2 = field_at_point(fieldz, x, z, xpoint[0], xpoint[1])
           
           denom = pt1^2 + pt2^2           
           dx = pt*pt1/denom
           dz = pt*pt2/denom           
           
           xpoint[0] = xpoint[0] - dx
           xpoint[1] = xpoint[1] - dz

           if(xpoint[0] le min(x) or xpoint[0] ge max(x)) then begin
               xpoint = 0
               break
           endif
           if(xpoint[1] le min(z) or xpoint[1] ge max(z)) then begin
               xpoint = 0
               break
           endif           
       end
   endif
   if(n_elements(xpoint) eq 2) then begin
       print, 'Found X-point at ', xpoint[0], xpoint[1]
   endif else begin
       print, 'No X-point found'
   endelse
end

pro nulls, psi, x, z, axis=axis, xpoints=xpoint, _EXTRA=extra
    version = read_parameter('version', _EXTRA=extra)
    if(version ge 3) then begin
        read_nulls, axis=axis, xpoint=xpoint, _EXTRA=extra
    endif else begin
        find_nulls, psi, x, z, axis=axis, xpoint=xpoint, _EXTRA=extra
    endelse 
end

