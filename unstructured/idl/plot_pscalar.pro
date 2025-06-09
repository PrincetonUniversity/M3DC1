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
                      integrate=integrate, ipellet=ipellet, _EXTRA=extra

   if(n_elements(scalarname) eq 0) then begin
       print, "Error: no scalar name provided"
       return, 0
   end

   if(n_elements(filename) eq 0) then filename='C1.h5'
   
   if(n_elements(ipellet) eq 0) then ipellet=0

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
   print, 'tag_names = ', n

   time = s.time._data
   ; print,"times_size = ",size(time,/DIMENSIONS)
   ; print,"time = ",time
   dims = size(time,/DIMENSIONS)
   ntimemax = dims[0]
   print,"ntimemax = ",ntimemax


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


   plot_id = 3 ; id-1
   ; use reform to remove all dims of size 1
   rdata = reform(pdata[0,plot_id,*])
   tdata = reform(pdata[1,plot_id,*])
   zdata = reform(pdata[2,plot_id,*])
   ;print,"rdata_size",size(rdata,/DIMENSIONS)
   ;print,"rdata = ",rdata

   print,"plot RZ trajectory of traced particle, #",plot_id
   title_str = 'RZ trajectory of traced particle, #'+STRTRIM(STRING(floor(plot_id+1)),2)
   ; PSYM sets marker style, COLOR sets marker color
   PLOT, rdata, zdata, TITLE=title_str, XTITLE='R axis', YTITLE='Z axis';, PSYM=4, COLOR="red"; , LINESTYLE=1
 

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
                 xabs=xabs, _EXTRA=extra

  if(n_elements(filename) eq 0) then filename='C1.h5'
  if(n_elements(xscale) eq 0) then xscale=1.
  if(n_elements(fac) eq 0) then fac=1.

  if(n_elements(names) eq 0) then names=filename

  nfiles = n_elements(filename)

;   data = read_pscalar(scalarname, filename=filename, time=time, ipellet=ipellet, $
;                      title=title, symbol=symbol, units=units, cgs=cgs, mks=mks, integrate=integrate)
  scalarname = "rz"
 
  data = read_pscalar(scalarname, filename=filename, _EXTRA=extra)
;   s = read_pscalars(filename=filename)
;   n = tag_names(s)
;   smatch = where(strcmp(n, scalarname, /fold_case) eq 1,scount)
;   print, 'tag_names = ', n
 
;  xtitle = make_label('!6Frequency!X', t0=-1, cgs=cgs, mks=mks, _EXTRA=extra)

;  data = power_spectrum(data, frequency=tdata, t=max(time))
 
end




