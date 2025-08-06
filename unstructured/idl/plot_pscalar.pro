;=====================================================
; read_scalars
; ~~~~~~~~~~~~
;
; returns a structure populated with the data from
; the "scalars" group of "filename"
;=====================================================
function read_pscalars, filename=filename
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
                      integrate=integrate, ipellet=ipellet, pid=pid, _EXTRA=extra

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
   ntimemax = dims[0]
   ;print,"ntimemax = ",ntimemax


; ;;; test output in time slice format  
;    pdata = s.pdata._data
;    print, "pdata_size = ", SIZE(pdata,/DIMENSIONS)
;    print,"pdata1 = ",pdata[*,0]
;    ;print,"pdata13 = ",pdata[*,12]


;;; output all traced particle information
   all_traced_particle_data = s.all_traced_particle._data
   npart_trace = s.npart_trace._data
   pdims_trace = s.pdims_trace._data
   pdata_trace_id = s.pdata_trace_id._data
;   print, "all_traced_particle_size = ", SIZE(all_traced_particle_data,/DIMENSIONS)
   print, "number of traced particles = ",npart_trace

   dims = SIZE(all_traced_particle_data,/DIMENSIONS)
   ;ntimemax = dims[1]
   ;print,"ntimemax = ",ntimemax
   pdims = dims[0]/npart_trace
   print,"traced particle data dims = ",pdims
   pdata = fltarr(pdims,npart_trace,ntimemax)
   FOR k = 0, ntimemax-1, 1 DO BEGIN
       pdata[*,*,k] = REFORM(all_traced_particle_data[*,k],[pdims,npart_trace])
   ENDFOR


   select_id = pid ; id-1
   
   if(scalarname eq "r") then begin
       ; use reform to remove all dims of size 1
       rdata = reform(pdata[0,select_id,*])
       ;print,"rdata_size",size(rdata,/DIMENSIONS)
       ;print,"rdata = ",rdata
       print,"reading R data"
       return, rdata
   endif

   if(scalarname eq "theta") then begin
       ; use reform to remove all dims of size 1
       tdata = reform(pdata[1,select_id,*])
       print,"reading Theta data"
       return, tdata
   endif

   if(scalarname eq "z") then begin
       ; use reform to remove all dims of size 1
       zdata = reform(pdata[2,select_id,*])
       print,"reading Z data"
       return, zdata
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
                 xabs=xabs, pid=pid, _EXTRA=extra

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

      rdata = read_pscalar("r", filename=filename, pid=pid, _EXTRA=extra)
      rtitle = "R axis"
      tdata = read_pscalar("theta", filename=filename, pid=pid, _EXTRA=extra)
      ttitle = "Theta axis"
      zdata = read_pscalar("z", filename=filename, pid=pid, _EXTRA=extra)
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
      xdata = read_pscalar("r", filename=filename, pid=pid, _EXTRA=extra)
      xtitle = "R axis"
      ydata = read_pscalar("z", filename=filename, pid=pid, _EXTRA=extra)
      ytitle = "Z axis"

      plot_dims = 2
  endif

  if(scalarname eq "rtheta") then begin
      rdata = read_pscalar("r", filename=filename, pid=pid, _EXTRA=extra)
      rtitle = "R axis"
      tdata = read_pscalar("theta", filename=filename, pid=pid, _EXTRA=extra)
      ttitle = "Theta axis"

      ; Convert cylindrical coordinates (r, theta, z) to Cartesian (x, y, z)
      xdata = rdata * COS(tdata)
      xtitle = "X axis"
      ;print,"xdata_size",size(xdata,/DIMENSIONS)
      ;print,"xdata = ",xdata
      ydata = rdata * SIN(tdata)
      ytitle = "Y axis"

      plot_dims = 2
;       title_str = 'R-Phi Trajectory of traced particle, #'+STRTRIM(STRING(floor(pid+1)),2)
;       ; PSYM sets marker style, COLOR sets marker color
;       pplot = POLARPLOT(xdata, ydata, TITLE=title_str, XTITLE=xtitle, YTITLE=ytitle) ;, PSYM=4, COLOR="red"; , LINESTYLE=1
  endif

  if(scalarname eq "thetaz") then begin
      xdata = read_pscalar("theta", filename=filename, pid=pid, _EXTRA=extra)
      xtitle = "Theta axis"
      ydata = read_pscalar("z", filename=filename, pid=pid, _EXTRA=extra)
      ytitle = "Z axis"

      plot_dims = 2
  endif

  if(scalarname eq "rtz") then begin   
      rdata = read_pscalar("r", filename=filename, pid=pid, _EXTRA=extra)
      rtitle = "R axis"
      tdata = read_pscalar("theta", filename=filename, pid=pid, _EXTRA=extra)
      ttitle = "Theta axis"
      zdata = read_pscalar("z", filename=filename, pid=pid, _EXTRA=extra)
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
   
      print,"Trajectory of traced particle, #",pid+1
      title_str = 'Trajectory of traced particle, #'+STRTRIM(STRING(floor(pid+1)),2)
      ; PSYM sets marker style, COLOR sets marker color
      PLOT, xdata, ydata, TITLE=title_str, XTITLE=xtitle, YTITLE=ytitle, PSYM=3, LINESTYLE=0 ;, PSYM=4, COLOR="red"; , LINESTYLE=1
  endif

;
end




