;==============================================
; hdf5_file_test
; ~~~~~~~~~~~~~~
; 
; returns 1 if "filename" is a valid hdf5 file;
; otherwise returns 0
; =============================================
function hdf5_file_test, filename
   if(n_elements(filename) ne 1) then begin

   endif else if(file_test(filename) eq 0) then begin

   endif else if(h5f_is_hdf5(filename) eq 0) then begin

   endif else return, 1

   print, "Error: ", filename, " is not a valid file."
   return, 0
end
