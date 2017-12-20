;=========================================================
; read_mesh
; ~~~~~~~~~
;
; Returns the mesh data structure at a given time slice
;=========================================================
function read_mesh, filename=filename, slice=t
   if(n_elements(filename) eq 0) then filename='C1.h5'
   if(n_elements(t) eq 0) then t = 0

   if(hdf5_file_test(filename) eq 0) then return, 0

   file_id = h5f_open(filename)

   time_group_id = h5g_open(file_id, time_name(t))
   mesh = h5_parse(time_group_id, 'mesh', /read_data)   

   h5g_close, time_group_id
   h5f_close, file_id

   return, mesh
end
