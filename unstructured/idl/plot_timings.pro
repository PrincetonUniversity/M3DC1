pro plot_timings, filename=filename, overplot=overplot, _EXTRA=extra

   if(n_elements(filename) eq 0) then filename = 'C1.h5'

   if(hdf5_file_test(filename) eq 0) then return

   file_id = h5f_open(filename)
   root_id = h5g_open(file_id, "/")
   timings = h5_parse(root_id, "timings", /read_data)
   h5g_close, root_id
   h5f_close, file_id

   v = read_parameter('version', filename=filename)

   t_solve = timings.t_solve_b._data + timings.t_solve_v._data + $
     timings.t_solve_n._data + timings.t_solve_p._data
   t_output = timings.t_output_cgm._data + timings.t_output_hdf5._data + $
     timings.t_output_reset._data


   loadct, 12
   c = get_colors()

   if(keyword_set(overplot)) then begin
       oplot, timings.t_onestep._data
   endif else begin
       plot, timings.t_onestep._data>0, title='!6Timings!3', $
         xtitle='!6Time Step!3', ytitle='!8t!6 (s)!3', _EXTRA=extra
   endelse
   oplot, timings.t_ludefall._data, linestyle=2, color=c[1]
   oplot, timings.t_sources._data, linestyle=1, color=c[2]
   oplot, timings.t_aux._data, linestyle=1, color=c[3]
   oplot, timings.t_smoother._data, linestyle=1, color=c[4]
   oplot, timings.t_mvm._data, linestyle=1, color=c[5]
   oplot, t_solve, linestyle=2, color=c[6]
   oplot, t_output, linestyle=2, color=c[7]
   if(v ge 20) then oplot, timings.t_kprad._data, linestyle=2, color=c[8]

   plot_legend, ['Onestep', 'ludefall', 'sources', 'aux', $
                 'smoother', 'mat vec mult', 'solve', 'output','kprad'], $
     linestyle=[0,2,1,1,1,1,2,2,2], color=c

end
