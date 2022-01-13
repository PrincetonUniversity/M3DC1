pro plot_wall_regions, filename=filename, slice=slice, overplot=overplot, $
  _EXTRA=extra

  if(n_elements(filename) eq 0) then filename = 'C1.h5'
  if(n_elements(slice) eq 0) then slice = 0

  file_id = h5f_open(filename)
  
  time_group_id = h5g_open(file_id, time_name(slice))
  att_id = h5a_open_name(time_group_id, "version")
  version = h5a_read(att_id)
  h5a_close, att_id

  print, 'version = ', version
  if(version lt 40) then begin
     print, 'Wall region data not found'
  endif else begin
     att_id = h5a_open_name(time_group_id, 'iwall_regions')
     nregions = h5a_read(att_id)
     h5a_close, att_id

     print, 'nregions = ', nregions

     for i=1, nregions do begin
        region_name = string(FORMAT='("wall_region_",I3.3)', i)
        region_id = h5g_open(time_group_id, region_name)
        att_id = h5a_open_name(region_id, "nplanes")
        nplanes = h5a_read(att_id)
        h5a_close, att_id

        print, 'nplanes = ', nplanes
        for j=1, nplanes do begin
           plane_name = string(FORMAT='("plane_",I3.3)', j)
           plane_id = h5g_open(region_id, plane_name)
           x_id = h5d_open(plane_id, "x")
           y_id = h5d_open(plane_id, "y")
           x = h5d_read(x_id)
           y = h5d_read(y_id)
           h5d_close, x_id
           h5d_close, y_id

           if(i eq 1 and j eq 1 and not keyword_set(overplot)) then begin
              plot, x, y, _EXTRA=extra
           endif else begin
              oplot, x, y, _EXTRA=extra
           end
           h5g_close, plane_id
        end
        
        h5g_close, region_id
     end
  end


  h5g_close, time_group_id
  h5f_close, file_id

  
end
