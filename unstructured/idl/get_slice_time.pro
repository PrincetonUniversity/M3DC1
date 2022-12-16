;========================================================
; get_slice_time
; ~~~~~~~~~~~~~~
;
; Returns the physical time associated with time slice
;========================================================
function get_slice_time, filename=filename, slice=slice, _EXTRA=extra
   if(n_elements(filename) eq 0) then filename='C1.h5'
   if(n_elements(slice) eq 0) then begin
      nt = read_parameter('ntime', filename=filename)
      slice = indgen(nt)
   end

   n = n_elements(slice)

   t = fltarr(n)

   get_normalizations, filename=filename, b0=b0, n0=n0, l0=l0, $
                       ion_mass=mi, _EXTRA=extra

   for i=0, n-1 do begin
       file_id = h5f_open(filename)
       time_group_id = h5g_open(file_id, time_name(slice[i]))
       time_id = h5a_open_name(time_group_id, "time")

       tt = h5a_read(time_id)
       convert_units, tt, dimensions(/t0), $
                      b0, n0, l0, mi, _EXTRA=extra
       t[i] = tt

       h5a_close, time_id
       h5g_close, time_group_id
       h5f_close, file_id
   end

   return, t
end
