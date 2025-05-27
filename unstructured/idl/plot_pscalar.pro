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

   s = read_scalars(filename=filename)
   n = tag_names(s)
   smatch = where(strcmp(n, scalarname, /fold_case) eq 1,scount)
   print, 'tag_names = ', n
   print,"e_kp_size = ",size(s.E_KP._data)
   print,"times_size = ",size(s.time._data)
   print,"time = ",s.time._data

   s = read_pscalars(filename=filename)
   n = tag_names(s)
   smatch = where(strcmp(n, scalarname, /fold_case) eq 1,scount)
   print, 'tag_names = ', n

   time = s.time._data
   print,"times_size = ",size(time,/DIMENSIONS)
   print,"time = ",time
  
   pdata = s.pdata._data
   print, "pdata_size = ", SIZE(pdata,/DIMENSIONS)
   print,"pdata1 = ",pdata[*,0]
   ;print,"pdata13 = ",pdata[*,12]

   test_particle_data = s.test_particle._data
   print,"test_particle_size = ",SIZE(test_particle_data,/DIMENSIONS)
   ntimemax = N_ELEMENTS(time)
   ;print,"time count: ",ntimemax
   npt = min([10,ntimemax])-1
   ;print,"npt=",npt
   print,"test_particle[0:3,0:10] = ",test_particle_data[0:3,0:npt]

   gid_test = test_particle_data[0,0]
   rdata = test_particle_data[1,*]
   tdata = test_particle_data[2,*]
   zdata = test_particle_data[3,*]
   print,"rdata = ",rdata

   title_str = 'RZ trajectory of test particle, id = '+STRTRIM(STRING(floor(gid_test)),2)
   ; PSYM sets marker style, COLOR sets marker color
   PLOT, rdata, zdata, TITLE=title_str, XTITLE='R axis', YTITLE='Z axis';, PSYM=4, COLOR="red"; , LINESTYLE=1
 
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




