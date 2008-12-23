function get_fieldnames, descriptions=descriptions
   fieldnames = ['beta', 'chi', 'eta', 'den', 'I', 'jphi', $
                 'kappa', 'Omega', 'p', 'pe', 'phi', 'pi', 'psi', $
                 'T', 'Te', 'Ti', 'V', 'vz']
   return, fieldnames
end

function get_scalarnames, descriptions=descriptions
   fieldnames = ['Beta', 'Toroidal Current', 'Kinetic Energy', $
                 'Loop Voltage', 'Particles']
   return, fieldnames
end

pro c1view_event, ev
   widget_control, ev.top, get_uvalue=state
   widget_control, ev.id, get_uvalue=uval

   if(n_elements(uval) eq 0) then return

   print, uval

   case uval of
       'AUTO_CT' : begin
           state.auto_ct = 1-state.auto_ct
           ct_button = widget_info(ev.top, find_by_uname='SET_CT')
           widget_control, ct_button, sensitive=1-state.auto_ct
       end
       'DONE' : widget_control, ev.top, /destroy 
       'FIELD_OPTIONS' : begin
           print, ev.value
           case ev.value of 
               'Contours': state.lines = 1-state.lines 
               'Isotropic' : state.isotropic = 1 - state.isotropic
               'LCFS' : state.lcfs = 1 - state.lcfs
               'Linear' : state.lcfs = 1 - state.linear
               'Mesh' : state.mesh = 1 - state.mesh
           end
       end
       'FIELDNAME' : begin
           fieldnames = get_fieldnames()
           state.fieldname = fieldnames(ev.index)
       end
       'FILENAME' : state.filename = ev.value
       'GAMMA' : state.growth_rate = 1-state.growth_rate
       'PLOT_FIELD' : begin 
           if(state.to_file eq 1) then begin
               print, 'plotting to file', state.output_filename
               setplot, 'ps'
               device, filename=state.output_filename, /color, /encapsulated
           endif
           plot_field, state.fieldname, state.slice, $
             filename=state.filename, $
             isotropic=state.isotropic, $
             lcfs=state.lcfs, $
             linear=state.linear, $ 
             mesh=state.mesh, $
             lines=state.lines, $
             range=state.zrange, zlog=state.zlog, $
             noautoct=1-state.auto_ct 
           if(state.to_file eq 1) then begin
               device, /close
               setplot, 'x'
           endif
       end
       'PLOT_SCALAR' : begin
           if(state.to_file eq 1) then begin
               print, 'plotting to file', state.output_filename
               setplot, 'ps'
               device, filename=state.output_filename, /color, /encapsulated
           endif
           plot_scalar, state.scalarname, $
             filename=state.filename[0:state.max_files-1], $
             xrange=state.xrange, yrange=state.yrange, $
             xlog=state.xlog, ylog=state.ylog, $
             growth_rate=state.growth_rate
           if(state.to_file eq 1) then begin
               device, /close
               setplot, 'x'
           endif
       end
       'SCALARNAME' : begin
           scalarnames = get_scalarnames()
           state.scalarname = scalarnames(ev.index)
       end
       'SLICE' : begin
           widget_control, ev.id, get_value=value
           state.slice = value
       end
       'SET_CT' : xloadct
       'TO_FILE' : state.to_file = 1-state.to_file
       'TO_FILENAME' : begin
           widget_control, ev.id, get_value=filename
           state.output_filename = filename
           print, 'filename = ', filename
       end
       'XLOG' : state.xlog = 1-state.xlog
       'XMIN' : begin
           widget_control, ev.id, get_value=xmin
           state.xrange[0] = float(xmin)
       end
       'XMAX' : begin
           widget_control, ev.id, get_value=xmax
           state.xrange[1] = float(xmax)
       end
       'YLOG' : state.ylog = 1-state.ylog
       'YMIN' : begin
           widget_control, ev.id, get_value=ymin
           state.yrange[0] = float(ymin)
       end
       'YMAX' : begin
           widget_control, ev.id, get_value=ymax
           state.yrange[1] = float(ymax)
       end
       'ZLOG' : state.zlog = 1-state.zlog
       'ZMIN' : begin
           widget_control, ev.id, get_value=zmin
           state.zrange[0] = float(zmin)
       end
       'ZMAX' : begin
           widget_control, ev.id, get_value=zmax
           state.zrange[1] = float(zmax)
       end

       else : print, 'No event defined'
   end

   widget_control, ev.top, set_uval=state, bad_id=id
end


pro filename_event, ev
   widget_control, ev.id, get_uval=uval
   widget_control, ev.top, get_uval=state
   widget_control, ev.id, get_value=value

   
   state.filename[uval-1] = value

   for i=0, n_elements(state.filename)-1 do begin
       if(strlen(state.filename[i]) eq 0) then break
   end
   state.max_files = i
   widget_control, ev.top, set_uval=state, bad_id=id

   print, state.max_files, state.filename
end

pro c1view
   state = { filename:['C1.h5', '', ''], max_files:1, $
             label:['1', '2', '3'], $
             fieldname:'beta', slice:0, $
             scalarname: 'beta', growth_rate:0, color_table:-1, auto_ct:1, $
             points:50, isotropic:1, lcfs:0, mesh:0, lines:0, linear:0, $
             xrange:[0.,0.], yrange:[0.,0.], zrange:[0.,0.], $
             xlog:0, ylog:0, zlog:0, $
             to_file:0, output_filename:'plot.eps' }

   base = widget_base(/row, /base_align_top)


   ; Field Plot Tab
   ; ~~~~~~~~~~~~~~
   tab = widget_tab(base, value='Field')
   base_field = widget_base(tab, title='Field Plot', /column)


   field_timeslice = cw_field(base_field, /all_events, uvalue='SLICE', $
                              /integer, value=0, title='Time Slice')

   list_field = widget_list(base_field, value=get_fieldnames(), ysize=5, $
                            uval='FIELDNAME')
   widget_control, list_field, set_list_select=0

   options = ['Contours', 'Isotropic', 'LCFS', 'Linear', 'Mesh']
   value = [state.lines,state.isotropic, state.lcfs, state.linear, state.mesh]
   group_plot_options = cw_bgroup(base_field, options, frame=1, $
                                  uvalue='FIELD_OPTIONS', column=2, $
                                  /nonexclusive, label_top='Plot Options', $
                                 set_value=value, /return_name)

   tab_zrange = widget_tab(base_field, value='zrange')
   base_zrange = widget_base(tab_zrange, /row, title='z Range')
   text_zmin = cw_field(base_zrange, value=state.zrange[0], /all_events, $
                           uvalue='ZMIN', xsize=6, title='', /float)
   text_zmin = cw_field(base_zrange, value=state.zrange[1], /all_events, $
                           uvalue='ZMAX', xsize=6, title='', /float)
   button_zlog = cw_bgroup(base_zrange, 'Log', uvalue='ZLOG', /nonexclusive)

   button_plot_field = widget_button(base_field, value='Plot', $
                                     uvalue='PLOT_FIELD')


   ; Scalar Plot Tab
   ; ~~~~~~~~~~~~~~~
   base_scalar = widget_base(tab, title='Scalar Plot', /column)



   scalar_list = widget_list(base_scalar, value=get_scalarnames(), ysize=5, $
                             uval='SCALARNAME')
   widget_control, scalar_list, set_list_select=0
   button_gamma = cw_bgroup(base_scalar, 'Growth Rate', uvalue='GAMMA', $
                            /nonexclusive)
   button_plot_scalar = widget_button(base_scalar, value='Plot', $
                                     uvalue='PLOT_SCALAR')


   ; Plot Options Tab
   ; ~~~~~~~~~~~~~~~~
   tab_plot = widget_tab(base, value='Plot')
   base_plot = widget_base(tab_plot, title='Plot Options', /column)

   tab_xrange = widget_tab(base_plot, value='xrange')
   base_xrange = widget_base(tab_xrange, /row, title='x Range')
   text_xmin = cw_field(base_xrange, value=state.xrange[0], /all_events, $
                           uvalue='XMIN', xsize=6, title='', /float)
   text_xmin = cw_field(base_xrange, value=state.xrange[1], /all_events, $
                           uvalue='XMAX', xsize=6, title='', /float)
   button_xlog = cw_bgroup(base_xrange, 'Log', uvalue='XLOG', /nonexclusive)
   tab_yrange = widget_tab(base_plot, value='yrange')
   base_yrange = widget_base(tab_yrange, /row, title='y Range')
   text_ymin = cw_field(base_yrange, value=state.yrange[0], /all_events, $
                           uvalue='YMIN', xsize=6, title='', /float)
   text_ymin = cw_field(base_yrange, value=state.yrange[1], /all_events, $
                           uvalue='YMAX', xsize=6, title='', /float)
   button_ylog = cw_bgroup(base_yrange, 'Log', uvalue='YLOG', /nonexclusive)


   tab_ct = widget_tab(base_plot, value='ct')
   base_ct = widget_base(tab_ct, /row, title='Color Table')
   button_ct = widget_button(base_ct, value='Set color table', $
                             uvalue='SET_CT', uname='SET_CT', $
                            sensitive=(1-state.auto_ct))
   button_ct_auto = cw_bgroup(base_ct, 'Auto', uvalue='AUTO_CT', $
                              /nonexclusive, set_value=state.auto_ct)

   button_to_file = cw_bgroup(base_plot, 'Output to file', $ 
                              uvalue='TO_FILE', /nonexclusive)
   field_to_file = cw_field(base_plot, value='plot.eps', $
                         uvalue='TO_FILENAME', /all_events)

   ; Filenames Tab
   ; ~~~~~~~~~~~~~
   base_filenames = widget_base(tab_plot, title='Filenames',column=2)

   label_filenames = widget_label(base_filenames, value='Filenames')
   text_file1 = widget_text(base_filenames, /all_events, uvalue=1, $
                            value=state.filename[0], /editable, $
                           event_pro='filename_event')
   text_file2 = widget_text(base_filenames, /all_events, uvalue=2, $
                            value=state.filename[1], /editable, $
                           event_pro='filename_event')
   text_file3 = widget_text(base_filenames, /all_events, uvalue=3, $
                            value=state.filename[2], /editable, $
                           event_pro='filename_event')
   label_labels = widget_label(base_filenames, value='Labels')
   text_label1 = widget_text(base_filenames, /all_events, uvalue=1, $
                             value='', /editable)
   text_label2 = widget_text(base_filenames, /all_events, uvalue=2, $
                             value='', /editable)
   text_label3 = widget_text(base_filenames, /all_events, uvalue=3, $
                             value='', /editable)


   button_done = widget_button(base, value='Done', uvalue='DONE')

   draw = widget_draw(base, xsize=500, ysize=400)

   widget_control, base, set_uvalue = state
   widget_control, base, /realize
   xmanager, 'c1view', base, /no_block
end

