;==================================================
; read_parameter
; ~~~~~~~~~~~~~~
;
; returns the value associated with the attribute
; "name" in the root group of "filename"
;==================================================
function read_parameter, name, filename=filename, print=pr, mks=mks, cgs=cgs
   if(n_elements(filename) eq 0) then filename='C1.h5'

   if(n_elements(filename) gt 1) then begin
       attr = fltarr(n_elements(filename))
       for i=0, n_elements(filename)-1 do begin
           attr[i] = read_parameter(name,filename=filename[i],print=pr)
       end
       return, attr
   endif

   if(hdf5_file_test(filename) eq 0) then return, 0

   file_id = h5f_open(filename)
   root_id = h5g_open(file_id, "/")
   attr = read_attribute(root_id, name)
   h5g_close, root_id
   h5f_close, file_id

   if(keyword_set(cgs) or keyword_set(mks)) then begin
      itor = read_parameter('itor', filename=filename)
      symbol = field_data(name, units=d, itor=itor)
      d0 = d
      get_normalizations, filename=filename,b0=b0,n0=n0,l0=l0,zeff=zeff,ion=mi
      convert_units, attr, d0, b0, n0, l0, zeff, mi, cgs=cgs, mks=mks
   end

   if(keyword_set(pr)) then print, name, " = ", attr

   return, attr
end

