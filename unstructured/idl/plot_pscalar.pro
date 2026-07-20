;=====================================================
; read_scalars
; ~~~~~~~~~~~~
;
; returns a structure populated with the data from
; the "scalars" group of "filename"
;=====================================================
function read_pscalars, filename=filename

   ; 1. Declare the universal common block (named 'global_h5_data')
   ;    We map the shared variable slot to our internal name 'pscalars'
   COMMON global_h5_data, pscalars
   
   ; If N_ELEMENTS is 0, the variable is undefined (file hasn't been read)
   IF N_ELEMENTS(pscalars) NE 0 THEN BEGIN
       PRINT, 'pscalars already exists in memory! Skipping read.'
       RETURN, pscalars
   ENDIF ELSE BEGIN
       PRINT, 'pscalars does not exist yet. Reading file...'
   ENDELSE


   if(n_elements(filename) eq 0) then filename='C1.h5'

   if(hdf5_file_test(filename) eq 0) then return, 0

;    ; check if this is a primitive field
;    file_id = h5f_open(filename)
;    time_group_id = h5g_open(file_id, time_name(time))
;    mesh = h5_parse(time_group_id, 'mesh', /read_data)

   file_id = h5f_open(filename)

;    root_id = h5g_open(file_id, "/")
;    pscalars = h5_parse(root_id, "particle_tracing", /read_data)
;    h5g_close, root_id

   ptrace_id = h5g_open(file_id, "particle_tracing")
   pscalars = h5_parse(ptrace_id, "particle_tracing_data", /read_data)
   h5g_close, ptrace_id
   
   h5f_close, file_id

   return, pscalars
end


function read_pscalar, scalarname, filename=filename, title=title, $
                      symbol=symbol, units=units, time=time, final=final, $
                      integrate=integrate, ipellet=ipellet, pid=pid, substep=substep, _EXTRA=extra

   ; 1. Declare the same common block. 
   ;    The first variable position links directly to 'pscalars' from above!
   COMMON global_h5_data, s
   COMMON global_pdata, pdata, pdata_substeps, nsubsteps_trace


   if(n_elements(scalarname) eq 0) then begin
       print, "Error: no scalar name provided"
       return, 0
   end

   if(n_elements(filename) eq 0) then filename='C1.h5'
   
   if(n_elements(ipellet) eq 0) then ipellet=0

   if(n_elements(scalarname) eq 0) then scalarname = "r"

   if(n_elements(pid) eq 0) then pid=0
   
   if(n_elements(filename) gt 1) then begin
       data = fltarr(n_elements(filename))
       for i=0, n_elements(filename)-1 do begin
           data[i] = read_pscalar(scalarname, filename=filename[i], $
                                 title=title, symbol=symbol, units=units, $
                                 time=time, /final)
       end
       return, data
   endif


   if( (n_elements(pdata) eq 0) ) then begin
       PRINT, 'Generating pdata and pdata_substeps...'

       ; ;;; C1.h5 read test
   ;    s = read_scalars(filename=filename)
   ;    n = tag_names(s)
   ;    smatch = where(strcmp(n, scalarname, /fold_case) eq 1,scount)
   ;    print, 'tag_names = ', n
   ;    print,"e_kp_size = ",size(s.E_KP._data)
   ;    print,"times_size = ",size(s.time._data)
   ;    print,"time = ",s.time._data
   
      s = read_pscalars(filename=filename)
      n = tag_names(s)
      smatch = where(strcmp(n, scalarname, /fold_case) eq 1,scount)
      ;print, 'tag_names = ', n
   
      time = s.time._data
      ; print,"times_size = ",size(time,/DIMENSIONS)
      ; print,"time = ",time
      dims = size(time,/DIMENSIONS)
      ntimemax0 = dims[0]
      ;print,"ntimemax = ",ntimemax
   
   
   ; ;;; test output in time slice format  
   ;    pdata = s.pdata._data
   ;    print, "pdata_size = ", SIZE(pdata,/DIMENSIONS)
   ;    print,"pdata1 = ",pdata[*,0]
   ;    ;print,"pdata13 = ",pdata[*,12]
   
   
   ;;; output all traced particle information
      all_traced_particle_data = s.all_traced_particle._data
      npart_trace = s.npart_trace._data[0]
      pdims_trace = s.pdims_trace._data[0]
      pdata_trace_id = s.pdata_trace_id._data

      partloss_count = s.partloss_count._data
      ;print, "partloss_size = ", SIZE(partloss_count,/DIMENSIONS)
   
   ;   print, "all_traced_particle_size = ", SIZE(all_traced_particle_data,/DIMENSIONS)
      print, "number of traced particles = ",npart_trace
   
   ;;; particle tracing for MHD steps only.
      dims = SIZE(all_traced_particle_data,/DIMENSIONS)
      print,'dims_pdata = ',dims
      ntimemax = dims[1]
      if(ntimemax NE ntimemax0) then print, 'Warning: ntimemax is not consistent with C1input, simulation may be interrupted.'
      print,"ntimemax = ",ntimemax
      pdims = pdims_trace ;dims[0]/npart_trace
      print,"traced particle data dims = ",pdims
      pdata = fltarr(pdims,npart_trace,ntimemax)
      FOR k = 0, ntimemax-1, 1 DO BEGIN
          pdata[*,*,k] = REFORM(all_traced_particle_data[*,k],[pdims,npart_trace])
      ENDFOR
   
   ;;; particle tracing with substeps
      nsubsteps_trace = s.nsubsteps_trace._data[0]
      ;HELP, nsubsteps_trace   ;;; <--- Add this temporarily
      if keyword_set(nsubsteps_trace>0) then begin
          print,'substep trajectory is exported for nsubsteps =',nsubsteps_trace
          all_traced_particle_substeps_data = s.all_traced_particle_substeps._data
        
          dims_substeps = SIZE(all_traced_particle_substeps_data,/DIMENSIONS)
          print,'dims_pdata_substeps = ',dims_substeps
          ;ntimemax = dims[1]
          ;print,"ntimemax = ",ntimemax
          ;pdims = pdims_trace ;dims_substeps[0]/npart_trace/nsubsteps_trace
          ;print,"traced particle data dims = ",pdims
          pdata_substeps = fltarr(pdims,npart_trace,nsubsteps_trace*ntimemax)
          FOR k = 0, ntimemax-1, 1 DO BEGIN
              ktmp = k*nsubsteps_trace
              pdata_substeps[*,*,ktmp:ktmp+nsubsteps_trace-1] = $
                  REFORM(all_traced_particle_substeps_data[*,k],[pdims,npart_trace,nsubsteps_trace])
          ENDFOR
          ;;;; test output:
          ;print,'all_traced_particle_sub_data[0:15,istep=1] =',all_traced_particle_substeps_data[0:15,1]
          ;print,'pdata_substeps[*,0,(istep=1)] = ',pdata_substeps[*,0,nsubsteps_trace]
      endif

   ENDIF ELSE BEGIN
       PRINT, 'Found pdata in memory!'
       if keyword_set(nsubsteps_trace>0) then print, "Found pdata_substeps in memory!"
   
   endelse   ; if(pdata and pdata_substeps not exist)

;;; extract data for specific particle pid.
   select_id = pid ; id-1

   ; ;;; test output:
   ; rdata_substeps = reform(pdata_substeps[0,select_id,*])
   ; print,"rdata_substeps_size",size(rdata_substeps,/DIMENSIONS)
   ; print,"rdata_substeps[490:510] = ",rdata_substeps[490:510]
   ; print,"rdata_substeps[990:1010] = ",rdata_substeps[990:1010]
   ; 
   ; rdata = reform(pdata[0,select_id,*])
   ; print,"rdata[0:4] = ",rdata[0:4]

   
   if(scalarname eq "r") then begin
       ; use reform to remove all dims of size 1
       rdata = reform(pdata[0,select_id,*])
       ;print,"rdata_size",size(rdata,/DIMENSIONS)
       ;print,"rdata = ",rdata
       if keyword_set(substep) and keyword_set(nsubsteps_trace>0) then rdata = reform(pdata_substeps[0,select_id,*])
       print,"reading R data"
       return, rdata
   endif

   if(scalarname eq "theta") then begin
       ; use reform to remove all dims of size 1
       tdata = reform(pdata[1,select_id,*])
       if keyword_set(substep) and keyword_set(nsubsteps_trace>0) then tdata = reform(pdata_substeps[1,select_id,*])
       print,"reading Theta data"
       return, tdata
   endif

   if(scalarname eq "z") then begin
       ; use reform to remove all dims of size 1
       zdata = reform(pdata[2,select_id,*])
       if keyword_set(substep) and keyword_set(nsubsteps_trace>0) then zdata = reform(pdata_substeps[2,select_id,*])
       print,"reading Z data"
       return, zdata
   endif

   if(scalarname eq "partloss_count") then begin
       ; use reform to remove all dims of size 1
       time = s.time._data
       ; print,"times_size = ",size(time,/DIMENSIONS)
       partloss_count = s.partloss_count._data
       ; print, "partloss_size = ", SIZE(partloss_count,/DIMENSIONS)
       
       time = time - time[0]

       ;if keyword_set(substep) and keyword_set(nsubsteps_trace>0) then zdata = reform(pdata_substeps[2,select_id,*])
       
       print,"reading partloss_count data"
       return, {time: time, data: partloss_count}
   endif

;     print,"plot RZ trajectory of traced particle, #",select_id
;    title_str = 'RZ trajectory of traced particle, #'+STRTRIM(STRING(floor(select_id+1)),2)
;    ; PSYM sets marker style, COLOR sets marker color
;    PLOT, rdata, zdata, TITLE=title_str, XTITLE='R axis', YTITLE='Z axis';, PSYM=4, COLOR="red"; , LINESTYLE=1
 

; ;;; test output in time history format
;    test_particle_data = s.test_particle._data
;    print,"test_particle_size = ",SIZE(test_particle_data,/DIMENSIONS)
;    ntimemax = N_ELEMENTS(time)
;    ;print,"time count: ",ntimemax
;    npt = min([10,ntimemax])-1
;    ;print,"npt=",npt
;    print,"test_particle[0:3,0:10] = ",test_particle_data[0:3,0:npt]
; 
;    gid_test = test_particle_data[0,0]
;    rdata = test_particle_data[1,*]
;    tdata = test_particle_data[2,*]
;    zdata = test_particle_data[3,*]
;    print,"rdata = ",rdata
; 
;    title_str = 'RZ trajectory of test particle, id = '+STRTRIM(STRING(floor(gid_test)),2)
;    ; PSYM sets marker style, COLOR sets marker color
;    PLOT, rdata, zdata, TITLE=title_str, XTITLE='R axis', YTITLE='Z axis';, PSYM=4, COLOR="red"; , LINESTYLE=1
 
end


pro plot_pscalar, scalarname, x, filename=filename, names=names, $
                 overplot=overplot, difference=diff, $
                 ylog=ylog, xlog=xlog, absolute_value=absolute, $
                 power_spectrum=pspec, per_length=per_length, $
                 growth_rate=growth_rate, bw=bw, nolegend=nolegend, $
                 cgs=cgs,mks=mks,linestyle=ls, color=co, outfile=outfile, $
                 smooth=sm, compensate_renorm=comp, integrate=integrate, $
                 xscale=xscale, ipellet=ipellet, factor=fac, versus=versus, $
                 xabs=xabs, pid=pid, substep=substep, _EXTRA=extra

  if(n_elements(filename) eq 0) then filename='C1.h5'
  if(n_elements(xscale) eq 0) then xscale=1.
  if(n_elements(fac) eq 0) then fac=1.

  if(n_elements(names) eq 0) then names=filename

  if(n_elements(scalarname) eq 0) then scalarname = "rz"

  if(n_elements(pid) eq 0) then pid=1
  
  ; idl index starts from 0
  pid = pid-1
  
  
  nfiles = n_elements(filename)
;   data = read_pscalar(scalarname, filename=filename, time=time, ipellet=ipellet, $
;                      title=title, symbol=symbol, units=units, cgs=cgs, mks=mks, integrate=integrate)

      rdata = read_pscalar("r", filename=filename, pid=pid, substep=substep, _EXTRA=extra)
      rtitle = "R axis"
      tdata = read_pscalar("theta", filename=filename, pid=pid, substep=substep, _EXTRA=extra)
      ttitle = "Theta axis"
      zdata = read_pscalar("z", filename=filename, pid=pid, substep=substep, _EXTRA=extra)
      ztitle = "Z axis"
 
  ; Open the file for writing (write mode)
  unit = 13
  OPENW, unit, 'test_particle_trajectory.out'
  
  ; Loop through the data and write each row
  FOR i = 0, N_ELEMENTS(rdata) - 1 DO BEGIN
      PRINTF, unit, FORMAT='(E16.8, E16.8, E16.8)', rdata[i], tdata[i], zdata[i]
  END
  
  ; Close the file
  CLOSE, unit

  if(scalarname eq "rz") then begin
      xdata = rdata ;read_pscalar("r", filename=filename, pid=pid, _EXTRA=extra)
      xtitle = "R axis"
      ydata = zdata ;read_pscalar("z", filename=filename, pid=pid, _EXTRA=extra)
      ytitle = "Z axis"

      print,"Trajectory of traced particle, #",pid+1
      title_str = 'Trajectory of traced particle, #'+STRTRIM(STRING(floor(pid+1)),2)
      plot_dims = 2
  endif

  if(scalarname eq "rtheta") then begin
      ;rdata = read_pscalar("r", filename=filename, pid=pid, _EXTRA=extra)
      rtitle = "R axis"
      ;tdata = read_pscalar("theta", filename=filename, pid=pid, _EXTRA=extra)
      ttitle = "Theta axis"

      ; Convert cylindrical coordinates (r, theta, z) to Cartesian (x, y, z)
      xdata = rdata * COS(tdata)
      xtitle = "X axis"
      ;print,"xdata_size",size(xdata,/DIMENSIONS)
      ;print,"xdata = ",xdata
      ydata = rdata * SIN(tdata)
      ytitle = "Y axis"

      print,"Trajectory of traced particle, #",pid+1
      title_str = 'Trajectory of traced particle, #'+STRTRIM(STRING(floor(pid+1)),2)
      plot_dims = 2
;       title_str = 'R-Phi Trajectory of traced particle, #'+STRTRIM(STRING(floor(pid+1)),2)
;       ; PSYM sets marker style, COLOR sets marker color
;       pplot = POLARPLOT(xdata, ydata, TITLE=title_str, XTITLE=xtitle, YTITLE=ytitle) ;, PSYM=4, COLOR="red"; , LINESTYLE=1
  endif

  if(scalarname eq "thetaz") then begin
      xdata = tdata ;read_pscalar("theta", filename=filename, pid=pid, _EXTRA=extra)
      xtitle = "Theta axis"
      ydata = zdata ;read_pscalar("z", filename=filename, pid=pid, _EXTRA=extra)
      ytitle = "Z axis"

      print,"Trajectory of traced particle, #",pid+1
      title_str = 'Trajectory of traced particle, #'+STRTRIM(STRING(floor(pid+1)),2)
      plot_dims = 2
  endif

  if(scalarname eq "rtz") then begin   
      ;rdata = read_pscalar("r", filename=filename, pid=pid, _EXTRA=extra)
      rtitle = "R axis"
      ;tdata = read_pscalar("theta", filename=filename, pid=pid, _EXTRA=extra)
      ttitle = "Theta axis"
      ;zdata = read_pscalar("z", filename=filename, pid=pid, _EXTRA=extra)
      ztitle = "Z axis"
   
      ; Convert cylindrical coordinates (r, theta, z) to Cartesian (x, y, z)
      xdata = rdata * COS(tdata)
      xtitle = "X axis"
      ;print,"xdata_size",size(xdata,/DIMENSIONS)
      ;print,"xdata = ",xdata
      ydata = rdata * SIN(tdata)
      ytitle = "Y axis"

      plot_dims = 3
  endif

  if(scalarname eq "partloss") then begin

      plot_data = read_pscalar("partloss_count", filename=filename, pid=pid, substep=substep, _EXTRA=extra)
      ; print,"plot_data_size",size(plot_data,/DIMENSIONS)
    
      xdata = plot_data.time ;read_pscalar("r", filename=filename, pid=pid, _EXTRA=extra)
      xtitle = "time"
      ydata = plot_data.data ;read_pscalar("z", filename=filename, pid=pid, _EXTRA=extra)
      ytitle = "Particle loss count"

      print,"Particle loss history plotted."
      title_str = "Particle loss history"
      plot_dims = 2
  endif

  if(plot_dims eq 3) then begin  

;       x = FINDGEN(100) / 10.0
;       y = SIN(x)
;       z = COS(x)
; 
;       my3DPlot = PLOT3D(x, y, z)

      print,"Real Trajectory of traced particle, #",pid+1
      title_str = '3D Trajectory of traced particle, #'+STRTRIM(STRING(floor(pid+1)),2)
      ; Create a 3D plot of the curve
      my3DPlot = PLOT3D(xdata, ydata, zdata) ;, TITLE=title_str, XTITLE=xtitle, YTITLE=ytitle, ZTITLE=ztitle )
  endif else if(plot_dims eq 2) then begin
   ;   s = read_pscalars(filename=filename)
   ;   n = tag_names(s)
   ;   smatch = where(strcmp(n, scalarname, /fold_case) eq 1,scount)
   ;   print, 'tag_names = ', n
    
   ;  xtitle = make_label('!6Frequency!X', t0=-1, cgs=cgs, mks=mks, _EXTRA=extra)
   
   ;  data = power_spectrum(data, frequency=tdata, t=max(time))
   
      ; print,"Trajectory of traced particle, #",pid+1
      ; title_str = 'Trajectory of traced particle, #'+STRTRIM(STRING(floor(pid+1)),2)
      ; PSYM sets marker style, COLOR sets marker color
      PLOT, xdata, ydata, TITLE=title_str, XTITLE=xtitle, YTITLE=ytitle, PSYM=3, LINESTYLE=0, _EXTRA=extra, ylog=ylog, xlog=xlog ;, PSYM=4, COLOR="red"; , LINESTYLE=1
      ;PLOT, xdata, ydata, TITLE=title_str, XTITLE=xtitle, YTITLE=ytitle, LINESTYLE=0 ;, PSYM=4, COLOR="red"; , LINESTYLE=1
      
      ; 2. Mark the START point (index 0) with a large circle/diamond
      ; PSYM=4 is a diamond, SYMSIZE increases the size, COLOR=255 assumes standard RGB red
      OPLOT, [xdata[0]], [ydata[0]], PSYM=4, SYMSIZE=2.5, COLOR='FF0000'x
      
      ; 3. Mark the LAST point (index -1) with a large X
      ; Using XYOUTS lets you print an actual "X" character exactly on the coordinates
      OPLOT, [xdata[-1]], [ydata[-1]], PSYM=7, SYMSIZE=2.5, COLOR='FF0000'x
      ;XYOUTS, xdata[-1], ydata[-1], 'X', ALIGNMENT=0.5, CHARSIZE=2.0, COLOR='0000FF'x
  endif

;
end




