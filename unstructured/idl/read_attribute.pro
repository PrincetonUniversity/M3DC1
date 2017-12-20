;================================================
; read_attribute
; ~~~~~~~~~~~~~~
;
; returns the value associated with the attribute
; "name" in "group_id"
;================================================
function read_attribute, group_id, name
   nmembers = h5a_get_num_attrs(group_id)

   for i=0, nmembers-1 do begin
       attr_id = h5a_open_idx(group_id, i)
       if(strcmp(name, h5a_get_name(attr_id), /fold_case)) then begin
           result = h5a_read(attr_id)
           h5a_close, attr_id
           return, result
       endif
       h5a_close, attr_id
   end
   
   print, "Error: cannot find attribute ", name

   return, 0
end
