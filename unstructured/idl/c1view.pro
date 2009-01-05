function get_fieldnames
   fieldnames = ['beta', 'chi', 'eta', 'den', 'I', 'jphi', $
                 'kappa', 'Omega', 'p', 'pe', 'phi', 'pi', 'psi', $
                 'T', 'Te', 'Ti', 'V', 'vz']
   return, fieldnames
end

function get_scalarnames
   fieldnames = ['Beta', 'Toroidal Current', 'Kinetic Energy', $
                 'Loop Voltage', 'Particles']
   return, fieldnames
end

function get_conservation_names
   fieldnames = ['Energy', 'Flux', 'Angular Momentum', 'Particle Flux']
   return, fieldnames
end


pro c1view_event, ev
   widget_control, ev.top, get_uvalue=state
   widget_control, ev.id, get_uvalue=uval

   if(n_elements(uval) eq 0) then return

   case uval of
       'CHANGE' : state.change = 1-state.change
       'CONS_NAME' : begin
           name = get_conservation_names()
           state.conservation_name = name(ev.index)
       end
       'CT_OPT' : begin
           state.ct_opt = ev.value
           if(ev.value eq 2) then loadct,0
           ct_button = widget_info(ev.top, find_by_uname='SET_CT')
           widget_control, ct_button, sensitive=(state.ct_opt eq 0)
       end
       'DONE' : widget_control, ev.top, /destroy 
       'FIELD_OPTIONS' : begin
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
       'XTICKS' : begin
           widget_control, ev.id, get_value=xticks
           state.xticks = float(xticks)
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
       'YTICKS' : begin
           widget_control, ev.id, get_value=yticks
           state.yticks = float(yticks)
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


pro plot_event, ev
   widget_control, ev.top, get_uvalue=state
   widget_control, ev.id, get_uvalue=uval

   if(n_elements(uval) eq 0) then return

   for i=0, n_elements(state.label)-1 do begin
       if(strlen(state.label[i]) eq 0) then $
         state.label[i] = state.filename[i]
   end

   if(state.to_file eq 1) then begin
       print, 'plotting to file', state.output_filename
       setplot, 'ps'
       device, filename=state.output_filename, /color, /encapsulated
   endif

   case uval of
       'PLOT_CONS' : begin
           plot_energy, state.conservation_name, $
             filename=state.filename[0], $
             diff=state.change, $
             xrange=state.xrange, yrange=state.yrange, $
             xlog=state.xlog, ylog=state.ylog
       end
       'PLOT_FIELD' : begin 
           plot_field, state.fieldname, state.slice, $
             filename=state.filename, $
             isotropic=state.isotropic, $
             lcfs=state.lcfs, $
             linear=state.linear, $ 
             mesh=state.mesh, $
             lines=state.lines, $
             range=state.zrange, zlog=state.zlog, $
             noautoct=(state.ct_opt eq 0), $
             xticks=state.xticks, yticks=state.yticks
       end
       'PLOT_SCALAR' : begin
           plot_scalar, state.scalarname, $
             filename=state.filename[0:state.max_files-1], $
             xrange=state.xrange, yrange=state.yrange, $
             xlog=state.xlog, ylog=state.ylog, $
             growth_rate=state.growth_rate, $
             bw=(state.ct_opt eq 2), $
             names=state.label[0:state.max_files-1]             
       end
   end

   if(state.to_file eq 1) then begin
       device, /close
       setplot, 'x'
   endif

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
end

pro label_event, ev
   widget_control, ev.id, get_uval=uval
   widget_control, ev.top, get_uval=state
   widget_control, ev.id, get_value=value
   
   state.label[uval-1] = value

   widget_control, ev.top, set_uval=state, bad_id=id
end


pro c1view
   state = { filename:['C1.h5', '', '', '', ''], max_files:1, $
             label:['', '', '', '', ''], $
             fieldname:'beta', slice:0, $
             scalarname:'beta', growth_rate:0, $
             conservation_name:'energy', change:0, $
             color_table:-1, ct_opt:1, $
             points:50, isotropic:1, lcfs:0, mesh:0, lines:0, linear:0, $
             monochrome:0, $
             xrange:[0.,0.], yrange:[0.,0.], zrange:[0.,0.], $
             xlog:0, ylog:0, zlog:0, xticks:0, yticks:0, $
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
                                     uvalue='PLOT_FIELD', $
                                     event_pro='plot_event')


   ; Scalar Plot Tab
   ; ~~~~~~~~~~~~~~~
   base_scalar = widget_base(tab, title='Scalar Plot', /column)



   scalar_list = widget_list(base_scalar, value=get_scalarnames(), ysize=5, $
                             uval='SCALARNAME')
   widget_control, scalar_list, set_list_select=0
   button_gamma = cw_bgroup(base_scalar, 'Growth Rate', uvalue='GAMMA', $
                            set_value=state.growth_rate, /nonexclusive)
   button_plot_scalar = widget_button(base_scalar, value='Plot', $
                                      uvalue='PLOT_SCALAR', $
                                      event_pro='plot_event')

   ; Convervation Laws Tab
   ; ~~~~~~~~~~~~~~~~~~~~~
   base_cons = widget_base(tab, title='Conservation Laws', /column)
   scalar_list = widget_list(base_cons, value=get_conservation_names(), $
                             ysize=5, uval='CONS_NAME')
   widget_control, scalar_list, set_list_select=0
   button_options = cw_bgroup(base_cons, 'Change', uvalue='CHANGE', $
                              set_value=state.change, /nonexclusive)
   button_plot_scalar = widget_button(base_cons, value='Plot', $
                                      uvalue='PLOT_CONS', $
                                      event_pro='plot_event')



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
   text_xticks = cw_field(base_xrange, title='Ticks', value=state.xticks, $
                          /all_events, uvalue='XTICKS', xsize=6, /integer)

   tab_yrange = widget_tab(base_plot, value='yrange')
   base_yrange = widget_base(tab_yrange, /row, title='y Range')
   text_ymin = cw_field(base_yrange, value=state.yrange[0], /all_events, $
                           uvalue='YMIN', xsize=6, title='', /float)
   text_ymin = cw_field(base_yrange, value=state.yrange[1], /all_events, $
                           uvalue='YMAX', xsize=6, title='', /float)
   button_ylog = cw_bgroup(base_yrange, 'Log', uvalue='YLOG', /nonexclusive)
   text_yticks = cw_field(base_yrange, title='Ticks', value=state.yticks, $
                          /all_events, uvalue='YTICKS', xsize=6, /integer)



   tab_ct = widget_tab(base_plot, value='ct')
   base_ct = widget_base(tab_ct, /row, title='Color Table')
   button_ct = widget_button(base_ct, value='Set color table', $
                             uvalue='SET_CT', uname='SET_CT', $
                            sensitive=(state.ct_opt eq 0))
   button_ct_auto = cw_bgroup(base_ct, ['User', 'Auto', 'Mono'], $
                              set_value=state.ct_opt, /row, $
                              uvalue='CT_OPT', /exclusive)

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
   text_file4 = widget_text(base_filenames, /all_events, uvalue=4, $
                            value=state.filename[3], /editable, $
                            event_pro='filename_event')
   text_file5 = widget_text(base_filenames, /all_events, uvalue=5, $
                            value=state.filename[4], /editable, $
                            event_pro='filename_event')

   label_labels = widget_label(base_filenames, value='Labels')
   text_label1 = widget_text(base_filenames, /all_events, uvalue=1, $
                             value=state.label[0], /editable, $
                             event_pro='label_event')
   text_label2 = widget_text(base_filenames, /all_events, uvalue=2, $
                             value=state.label[1], /editable, $
                             event_pro='label_event')
   text_label3 = widget_text(base_filenames, /all_events, uvalue=3, $
                             value=state.label[2], /editable, $
                             event_pro='label_event')
   text_label4 = widget_text(base_filenames, /all_events, uvalue=4, $
                             value=state.label[3], /editable, $
                             event_pro='label_event')
   text_label5 = widget_text(base_filenames, /all_events, uvalue=5, $
                             value=state.label[4], /editable, $
                             event_pro='label_event')


   button_done = widget_button(base, value='Done', uvalue='DONE')

   draw = widget_draw(base, xsize=600, ysize=400)

   widget_control, base, set_uvalue = state
   widget_control, base, /realize
   xmanager, 'c1view', base, /no_block
end

