;+ 
; NAME: 
;	EFIT_READ
;
; PURPOSE:
;
;	Handles EFIT retrieval logic common to both READA and READG
;
; CATEGORY: 
;
;	DIII-D development 
;
; CALLING SEQUENCE: 
;
;	data = EFIT_READ(arg1 [,arg2] ,TYPE=type [,MODE=mode] [,RUNID=runid] 
;	                      [,INFO=info] [,SOURCE=source] 
;			      [,EXACT_TIME=exact_time]
;			      [,VERBOSE=verbose] [,SERVER=server]
;			      [,DEBUG=debug] [,STATUS=status] )
;
; INPUT PARAMETERS: 
;
;	arg1:  	Either a string specifying the filename to read, or a long 
;               integer specifying the shot number (if the latter, arg2 must 
;               be present) to read.
;
;	arg2:	Float or double specifying the EFIT time to read and return.  
;               If arg1 specifies the shot number, arg2 must be used.
;
; OPTIONAL INPUT PARAMETERS: 
;
;	none 
;
; KEYWORDS: 
;
; 	TYPE: *mandatory* Either "a" or "g" specifying whether AEQDSK or
;	GEQDSK data is required.
;
;	(all the rest are optional)
;
;	
;	MODE:  If "FILE", will restrict EFIT_READ to retrieving EFIT data from
;	       files only, not from MDSplus.  If "MDSPLUS", will restrict 
;	       EFIT_READ to retrieving EFIT data from MDSplus only, not files.
;	       If not specified, EFIT_READ will first attempt to retrieve the 
;	       data from a file, and then from MDSplus.
;	
;	RUNID:  EFIT "run ID" to use in MDSplus.  This defaults to "EFIT01" - 
;		the non-MSE automatic control room EFIT.
;	
;	INFO:  A structure with the following form:
;	
;		{mode:'', file:'', shot:0l, time:0.0d0, runid:''}
;	
;	       If specified as an input to EFIT_READ, INFO will superceed the 
;	       arguments specified in arg1, arg2, and the keyword values of 
;	       MODE and RUNID.
;	
;	       INFO is also returned from EFIT_READ to indicate the values it 
;	       used to find the EFIT.
;	
;	SOURCE:  Either "FILE" or "MDSPLUS" - specifies the data source from
;	         where the EFIT data were retrieved.
;	
;	EXACT_TIME:  If set, forces EFIT_READ to match the time specified, rather
;		     than using the default behavior of returning the nearest
;		     time.  
;
;	VERBOSE:  If set, EFIT_READ will print out informational messages on its
;	          progress.
;	
;       SERVER:  If set to a string containing a valid IP address for
;                an MDSplus data server, EFIT_READ will read the EFITs
;                from the specified server instead of the default for DIII-D.
;	
;	DEBUG:  If set, EFIT_READ will print out additional debugging 
;		information, as well as turn off all error handling.  This 
;	        will allow EFIT_READ to crash if there is an unexpected error.
;	
;	STATUS:  TRUE if EFIT_READ was able to retrieve the data successfully,
;                FALSE if not. This information is also provided in the 
;		 output of the function (see below).
;
;
; OUTPUTS: 
;
;	structure containing the AEQDSK, GEQDSK or MEASUREMENTS data 
;	retrieved for the shot and time specified.
;
;	Note that the structure contains the tag ERROR, which is 0 if the
;	data was read successfully, and 1 if not.
;
; COMMON BLOCKS: 
;
;       None.   If GEQDSK or MEASUREMENTS data is requested, the routines
;	READG_MDSREAD or READM_MDSREAD are called.  These have common blocks
;
;		COMMON EFIT_READG_CACHE,info_cache,data_cache
;		COMMON EFIT_READM_CACHE,info_cache,data_cache
;
; SIDE EFFECTS: 
;
;	Calls back into functions stored in READA.PRO, READG.PRO, or READM.PRO 
;	to retrieve the AEQDSK, GEQDSK, or MEASUREMENTS data.  These functions 
;	are only automatically compiled if READA, READG, or READM are compiled
;	before the call to EFIT_READ.
;
; RESTRICTIONS:
;
;	Requires that READA.PRO, READG.PRO, and READM.PRO be compiled first.
;
; PROCEDURE: 
;
;	EFIT_READ retrieves data from an EFIT run for a particular shot and
;	timeslice.  It is intended to be called only from READA, READG, or
;	READM, not directly by the user.  This is because it calls back into
;	functions stored in READA.PRO, READG.PRO, and READM.PRO.  These
;	functions will only be automatically compiled if READA (or READG or
;	READM) is called or compiled before the call to EFIT_READ.
;	
;	EFIT_READ is used to save duplication of the logic used to retrieve
;	the AEQDSK, GEQDSK, or MEASUREMENTS data.  It handles all of the 
;	functions common to these retrievals.  Functions specific to the type 
;	of retrieval are handled by the functions and procedures in READA.PRO,
;	READG.PRO, or READM.PRO.
;	
;	- If arg1 specifies a file, EFIT_READ attempts to determine the 6 digit
;	  shot number and time from the filename, assuming it has the format
;	  .../XSSSSSS.TTTTT_TTT.  _TTT is optional - used for specifying
;	  sub-millisecond timeslices.  If it cannot, it will still attempt to 
;	  read the file specified, but if the file attempt fails, the MDSplus
;	  attempt will also fail (see below). "X" is the file type - either "a"
;	  "g" or "m".  NOTE THAT if arg1 specifies a file, EFIT_READ will act 
;	  as if the EXACT_TIME keyword is set - that is it will not attempt to 
;	  find the nearest time if it cannot find the exact time.
;	
;	- If arg1 specifies the shot and arg2 the time, the filename
;	  XSSSSSS.TTTTT_TTT is composed, where X is the file type specified by
;	  TYPE, SSSSSS is the 6 digit shot number and TTTTT_TTT is the time, 
;	  optionally using the _TTT for sub-ms timeslices.
;	
;	- If the filename contains a directory specification, EFIT_READ will 
;         look for the file specified in the place specified.  If it does not,
;         EFIT_READ will search the following locations (in order) for the file:
;	
;	  1) The current directory ./  (VMS: [])
;	  2) The subdirectory ./shotSSSSSS  (VMS: [.shotSSSSSS])
;	  3) The subdirectory ./shotSSSSS  (VMS: [.shotSSSSS])    
;					(SSSSS = 5 digit shot number)
;	
;	- If the file is found in one of these places, an attempt is made to 
;	  read it.
;
;	- If the file is not found in one of these three places, and the 
;	  keyword EXACT_TIME is *NOT* set, the same locations are searched 
;	  for an AEQDSK, GEQDSK, or MEASUREMENTS file with a time *nearest* 
;	  the time specified.  
;	
;	- If that attempt fails, or if te file is not found, EFIT_READ will 
;         attempt to read the data from MDSplus, using the shot number and 
;         time specified (or determined from the filename).  Data from the 
;	  time *nearest* the time specified will be returned if the MDSplus 
;	  read attempt is successful (unless the keyword EXACT_TIME is set).
;	
;	- If the value of the keyword MODE is "MDSPLUS", EFIT_READ will not
;	  attempt to read the data from a file, instead proceeding directly to
;	  the MDSplus read attempt. 
;	
;	- If the value ofthe keyword MODE is "FILE", EFIT_READ will not attempt
;	  to read the data from MDSplus if the file read attempt is 
;         unsuccessful.
;	
;	- Any other value of MODE will have no effect on the read logic.
;	
; EASE OF USE: Can be used with existing documentation
;
; OPERATING SYSTEMS:  HP-UX, OSF/Unix, OpenVMS, MacOS
;
; EXTERNAL CALLS:  MDSplus
;
; RESPONSIBLE PERSON: Jeff Schachter
;
; CODE TYPE: modeling, analysis, control  
;
; CODE SUBJECT:  handling, equilibrium
;
; DATE OF LAST MODIFICATION:  6/04/01
;
; MODIFICATION HISTORY:
;
;	Version 1.0: Released by Jeff Schachter 98.03.19
;	   98.03.25: EFIT_READ_FILESEARCH: only do searchpath for file if file
;                    does not already contain directory
;	   98.04.21: Debugged for VMS by Bill Davis
;	Version 2.0: Jeff Schachter 98.04.27
;                    Standardized behavior in finding nearest time, and added 
;		     keyword EXACT_TIME 
;	Version 3.0: Jeff Schachter 98.05.08
;		     Added more functionality with respect to finding nearest 
;		     time
;       Version 3.1: For workstations without MDSplus installed,
;                    calling with MODE="FILE" will now work
;       Version 3.2: Jeff Schachter 1998.10.06
;                   - modify calls to MDSplus functions so that this
;                     procedure works with both client/server and
;                     native access
;       1999.11.02   Jeff Schachter - moved readg_mdsreadall to
;                                     efit_read_mdsreadall (general
;                                     code)
;       1999.11.03   Jeff Schachter - trim strings returned from
;                                     MDSplus so that structure tags
;                                     are set ok
;       1999.11.17   Jeff Schachter - fix bug in handling RUNID
;                                     specified as null string
;       2000.08.01   Jeff Schachter - make test of mode case
;                                     insensitive
;       2000.09.01   Jeff S.: fix check of RUNID.  Move from setinfo
;                    to efit_read_mds so that connection to server is already
;                    established before call to mdsvalue('EFITTREE()') 
;                    is made.
;       2001.06.04   Jeff S.: add SERVER keyword.
;       2002.02.15   Ted Terpstra - corrected error of Porters where
;                    g=readg(106005, 3000) returned g.lim=[89,2] when it
;                    should have been                     [2,89].
;                    Coincidentally, 2*89 = 178 the number of times in shot.
;                    This forced a transpose. 
;                    Now only look at sz[0:nsz-3], where sz = size(data). 
;                    Change in efit_read_mdsreadall.
;	2003.12.09   Q.Peng, make sure efit_read_timematch returns a scalar
;		     rather than an array of size 1 (fix for IDL 6)
;	2007.06.27   Q.P. Use nodes themselves instead of their children's
;                    parents to determine NODE_NAME, FULLPATH for M trees
;                    as not all nodes under M have child "label".
;-

pro efit_read_message,verbose,msg
  if (verbose) then print,msg
end

function efit_read_error
  return,{error:1}
end

function efit_read_timerange,time,time_range
  t0 = time - time_range[0]	
  if (n_elements(time_range) eq 2) then begin
    t1 = time + time_range[1]
  endif else t1 = time + time_range[0]
  return,{t0:temporary(t0), t1:temporary(t1)}
end

function efit_read_timematch,timereq,time_range,times
  tr = efit_read_timerange(timereq,time_range)
  itemp = where(times ge tr.t0 and times le tr.t1, ntemp)
  if (ntemp gt 0) then begin
    delta = min(abs(times[itemp] - timereq),ix)
    itime = itemp[ix]
  endif else itime = -1
  return,itime[0]
end


;-------------------------------------------------------------------------------------

function efit_read_filenearest,path,type,shot,time,time_range,count=count,debug=debug
  forward_function efit_shotfiles
  info = efit_shotfiles('.',subdir=path,shot=shot,types=type)
  if (keyword_set(debug)) then help,path,shot,type,info,/str

  if (info.nshots eq 1) then begin

    if (keyword_set(time_range)) then begin
      tr = efit_read_timerange(time,time_range)
      itime = efit_read_timematch(time,time_range,(*info.ptrs[0]).dep.times)

    endif else begin

      dummy = min(abs((*info.ptrs[0]).dep.times - time),itime)

    endelse

    if (itime ge 0) then begin
      fileFound = FILE_SEARCH((*info.ptrs[0]).dep.files[itime],count=count)
    endif else begin
      count = 0
      fileFound = ''
    endelse

    ptr_free,info.ptrs	 ; only if info.shots > 0 will info.ptrs be defined.

  endif else begin
    count = 0
    fileFound = ''
  endelse
  return,fileFound[0]

end
  
function efit_read_filesearch,info,exact_time=exact_time,time_range=time_range,$
				time_indep=time_indep,debug=debug,status=status

  ;====== find file in search path

  shot6str = 'shot'+string(info.shot,format='(i6.6)')
  shotstr = 'shot'+strtrim(info.shot,2)

  case (!VERSION.OS_FAMILY) of 
    'unix' : begin
      ; if file starts with . or / then directory is specified, otherwise, add search path
      if ((strpos(info.file,'/') ne 0) and (strpos(info.file,'.') ne 0) ) then begin 
        searchpath = ['./', './'+shot6str+'/',  './'+shotstr+'/']	
      endif else searchpath = [''] 					
    end
    'vms' : begin
      ;if file contains [ or : then directory is specified, otherwise, add searchpath
      if ( (strpos(info.file,'[') EQ -1) and (strpos(info.file,':') EQ -1) ) then begin 
	searchpath = ['[]', '[.'+shot6str+']', '[.'+shotstr+']']
      endif else searchpath = ['']
    end
    'MacOS' : begin
;NOTE: FILL IN FILE SEARCH PATH INFO FOR MAC LATER
        searchpath = [''] 
    end

    else : message, !VERSION.OS_FAMILY+" operating system not supported"
  endcase


  ;==== first look for exact file in any of the specified subdirs

  count = 0
  i = 0

  while (count eq 0 and i lt n_elements(searchpath)) do begin
    fileFound = (FILE_SEARCH(searchpath[i]+info.file, count=count))[0]
    if (keyword_set(debug)) then print,"Searched for: "+searchpath[i]+info.file+ $
        "   Found: "+fileFound
    i = i + 1
  endwhile

  status = (count eq 1)

  ;==== if not found, and time independent files are allowed, look for time in independent file
  if (not(status) and keyword_set(time_indep)) then begin
    filetmp = 'm'+string(info.shot,format='(i6.6)')+'.nc'
    i = 0
    while (count eq 0 and i lt n_elements(searchpath)) do begin
      fileFound = (FILE_SEARCH(searchpath[i]+filetmp, count=count))[0]
      if (keyword_set(debug)) then print,"Searched for: "+searchpath[i]+filetmp+ $
          "   Found: "+fileFound
      i = i + 1
    endwhile
    status = (count eq 1)
  endif

  ;==== if still not found, find closest matching time from files

  if (not(status) and not(keyword_set(exact_time))) then begin
    i = 0
    while (count eq 0 and i lt n_elements(searchpath)) do begin
      fileFound = efit_read_filenearest(searchpath[i],info.type,info.shot,info.time,time_range,count=count,debug=debug)
      if (keyword_set(debug)) then begin
        print,"Searching: "+searchpath[i]+" for nearest match to: "+strtrim(info.shot,2)+" "+strtrim(info.time,2)
        print,"FOUND: ",fileFound
      endif
      i = i + 1
    endwhile
    status = (count eq 1)
  endif

  return,fileFound

end

;-------------------------------------------------------------------------------------

function efit_read_setinfo, arg1, arg2, type=type, mode=mode, runid=runid, exact_time=exact_time, $
					verbose=verbose, debug=debug, status=status

  forward_function efit_filename_parse

  ;--- Initialize INFO structure

  info={mode:'', type:type, file:'', shot:0l, time:0.0d0, runid:'EFIT01'}  ; default to EFIT01

  size1 = size(arg1)
  size2 = size(arg2)

  ;--- Set MODE and RUNID

  if (n_elements(mode) gt 0) then info.mode = strupcase(mode)

  ;--- Set SHOT, TIME, FILE 

  if (size1(n_elements(size1)-2) eq 7) then begin      ; arg1 is a string

    info.file = arg1
    ; parse filename to get shot and time if can
    fileparse=efit_filename_parse(info.file) 
    info.shot=fileparse.shots(0)  ; will be -1 if can't parse the shot number from the filename
    info.time=fileparse.times(0)  ; will be -1 if can't parse the time from the filename
    if (info.time eq -1.) then if (keyword_set(arg2)) then info.time = double(arg2)
  
    exact_time = 1  ; force finding exact time!

    status = 1

  endif else begin

    if (keyword_set(arg2)) then begin 			; arg2 must be present (SHOT,TIME)

      info.shot = long(arg1)
      info.time = double(arg2)
      info.file = type + string(info.shot,format='(i6.6)') + '.' + string(info.time,format='(i5.5)')

      ; allow sub-millisecond times g099999.12345_123

      delta = info.time - long(info.time)
      if (delta gt 0) then begin
	s=strtrim(delta,2)
	info.file = info.file + '_' + strmid(s,strpos(s,'.')+1,strlen(s))
      endif

      status = 1

    endif else begin

      status = 0

    endelse

  endelse

  if (n_elements(mode) gt 0) then stat=(strupcase(strtrim(mode,2)) ne 'FILE') else stat=1
  if (n_elements(runid) eq 1) then begin
    info.runid = runid 
  endif else begin 
    ;; 2001.06.07 JMS - scratch trees have 8 digit shot numbers (at DIII-D)
    if (info.shot gt 999999) then info.runid='EFIT' else info.runid='EFIT01'
  endelse

  ;--- Done setup

  if (keyword_set(debug)) then begin
    efit_read_message,1,'EFIT_READ_SETINFO: Information set: '
    help,status
    help,info,/str
  endif


  return,info

end
    

;-------------------------------------------------------------------------------------

function efit_read_mds_time,efittime,timereq,exact_time=exact_time,time_range=time_range, $
                            maxdt=maxdt,verbose=verbose,debug=debug,status=status

  ;==== FIND NEAREST TIME, or EXACT TIME IF REQUESTED

  ;== first look for exact match
  itime = where(efittime eq timereq,nmatch)
  status = (nmatch eq 1)

  ;== if no exact match, what to do depends on keywords
  if (nmatch eq 0) then begin
    case (1) of 

      ;=  return failure
      keyword_set(exact_time) : status = 0  

      ;= look for nearest only within specified range
      keyword_set(time_range) : begin	
	itime = efit_read_timematch(timereq,time_range,efittime)
	if (keyword_set(debug)) then print,'Found time: ',itime
	status = (itime ge 0)
      end

      ;= look for nearest within maxdt range
      keyword_set(maxdt) : begin
	delta = min(abs(efittime - timereq),itime)
	if (delta le maxdt) then status = 1 else status = 0
      end

      ;= look for nearest anywhere
      else : begin
        delta = min(abs(efittime - timereq),itime) 
        status = 1
      end

    endcase
  endif
  return,itime[0]

end

function efit_read_mdsreadall,info,verbose=verbose,debug=debug,status=status

  forward_function efit_read_error, efit_read_mds_tags, mdsvalue, efit_read_getmdstime
  quiet = 1 - (keyword_set(debug))

  time = efit_read_getmdstime(info.type,debug=debug,status=status)
  code_version = mdsvalue('\TOP.RESULTS:CODE_VERSION',quiet=quiet)
  if (status) then begin
    s={shot:info.shot, time:time, code_version:code_version, error:0}
    if (status) then begin
      siginfo = efit_read_mds_tags(info.type,debug=debug)
      status = siginfo.status
      if (status) then begin
        for i=0,n_elements(siginfo.parents)-1 do begin
          if (keyword_set(debug)) then print,siginfo.parents[i],' --- ',siginfo.tags[i]
          data=mdsvalue('_s = '+siginfo.parents[i],quiet=quiet,status=statread)
          if (statread) then begin
            ndims = size(data,/n_dimensions) 
            if (ndims gt 1) then begin
              ;;; if last dimension is not time, transpose array.
              ;;; Hopefully only 2D arrays will need transposing.
              ;;; If a 3D array (say for e.g. psirz[t,r,z]) is
              ;;; transposed, then you end up with psirz[z,r,t]
              ;;; which is still not good.
              sz = size(data)

;;;;          tbt 20020215  fixes g=readg(106005, 3000) which had g.lim=[89,2]
;;;;          Last two elements of size are type and total elements in array.
	      nsz = N_elements(sz)
              sz = sz[0:nsz-3]
;;;;---------------------------------------

              itimedim = where(sz eq n_elements(time),ntimedim)

;;;;          Make sure lim is never transposed. tbt 20020215
              If siginfo.tags[i] Eq 'lim' Then ntimedim = 0

              case (ntimedim) of
                0 : 
                1 : if (itimedim[0] ne sz[0]) then data=temporary(transpose(data))
                else : begin
                  unitslast = mdsvalue('UNITS(DIM_OF(_s,$))',ndims-1,quiet=quiet,status=statunits)
                  if (statunits) then begin
                    if (unitslast ne 's' and unitslast ne 'ms') then data = temporary(transpose(data)) 
                  endif
                end
              endcase
            endif
            s=create_struct(temporary(s),siginfo.tags[i],data)
          endif
        endfor
      endif else s.error=1
    endif else s.error=1
  endif else s.error=1
  return,s

end
function efit_read_mds_tags,type,debug=debug

  quiet = (1-keyword_set(debug))

  if (strupcase(type) eq 'M') then begin
    checkHeader = '\EFIT_MEASUREMENTS'
    fallbackHeader = '\TOP.MEASUREMENTS'
    child = 'LABEL'
  endif else begin
    checkHeader = '\EFIT_'+type+'EQDSK'
    fallbackHeader = '\TOP.RESULTS.'+type+'EQDSK'
    child = 'READA_NAME'
  endelse

  headerPath = mdsvalue('GETNCI($,"MINPATH")',checkHeader,/quiet,status=statheader)
  if (not(statheader)) then headerPath = fallbackHeader

  expression = headerPath + ":*"
  if (strupcase(type) eq 'M') then $  ; not all nodes in M have "child"
  nids = mdsvalue('GETNCI($,"NID_NUMBER")',expression,quiet=quiet,status=status) else $
  nids = mdsvalue('GETNCI($,"NID_NUMBER")',expression+":"+child,quiet=quiet,status=status)

  if (status) then begin
    if (strupcase(type) eq 'M') then begin
      parents = strtrim(mdsvalue('GETNCI($,"FULLPATH")',nids,quiet=quiet,status=status),2)
      nodes = strtrim(mdsvalue('GETNCI($,"NODE_NAME")',nids,quiet=quiet,status=status),2)
    endif else begin
      parents = strtrim(mdsvalue('GETNCI(GETNCI($,"PARENT"),"FULLPATH")',nids,quiet=quiet,status=status),2)
      nodes = strtrim(mdsvalue('GETNCI(GETNCI($,"PARENT"),"NODE_NAME")',nids,quiet=quiet,status=status),2)
    end
    if (status) then begin

      ;=== MEASUREMENTS signals do not have "READA_NAME" nodes underneath

      if (strupcase(type) eq 'M') then begin
	tags = strtrim(nodes,2)
      endif else begin
	tags = strtrim(mdsvalue('GETNCI($,"RECORD")',nids,quiet=quiet,status=status),2)
      endelse

      ;=== all three arrays - PARENTS, NODES, and TAGS, are aligned
      ;=== This is guaranteed by MDSplus

      ;=== substitute node name for tag name if READA_NAME is empty

      i = where(strcompress(tags,/remove_all) eq '',n)
      if (n gt 0) then tags[i] = nodes[i] 

      d = {nodes:temporary(nodes), parents:temporary(parents), tags:temporary(tags), status:status} 
    endif else d = {status:status}

  endif else begin ; Do not look for children (C-Mod)

    nids = mdsvalue('GETNCI($,"NID_NUMBER")',expression,quiet=quiet,status=status)
    if (status) then begin
      parents = strtrim(mdsvalue('GETNCI($,"FULLPATH")',nids,quiet=quiet,status=status1),2)
      nodes = strtrim(mdsvalue('GETNCI($,"NODE_NAME")',nids,quiet=quiet,status=status2),2)
      status = status1 and status2
      if (status) then begin
        ind=where(nodes NE 'CASE') ; illegal tag name from CMOD
        d = {nodes:nodes[ind], parents:parents[ind], tags:nodes[ind], status:status}
      endif else d = {status:status}
    endif else d={status:status}
  endelse

  return,d

end


;-------------------------------------------------------------------------------------
function efit_read_file,info,exact_time=exact_time,time_range=time_range, $
			     verbose=verbose, debug=debug, status=status

  forward_function reada_file, readg_file, readm_file

  info.file = strtrim(info.file,2)
  file = efit_read_filesearch(info, exact_time=exact_time, time_range=time_range, $
				    time_indep=(strupcase(info.type) eq 'M'), $
	                            debug=debug, status=status)

  if (status) then begin

    info.file = file
    efit_read_message,verbose,'EFIT_READ: reading file '+file

    ;===== readX_file functions in file readX.pro
    case (strupcase(info.type)) of 
      'A' : data = reada_file(info, verbose=verbose, debug=debug, status=status)
      'G' : data = readg_file(info, verbose=verbose, debug=debug, status=status)
      'M' : data = readm_file(info, verbose=verbose, debug=debug, status=status)
      else : begin
	status = 0 
	efit_read_message,verbose,'EFIT_READ: Unrecognized file type: '+info.type
	if (keyword_set(debug)) then stop
      end
    endcase

  endif else begin
    data = efit_read_error()
    efit_read_message,verbose,"EFIT_READ: Could not find file: "+info.file
  endelse

  return,data
end

;-------------------------------------------------------------------------------------

function efit_read_getmdstime,type,debug=debug,status=status
  quiet = (1-keyword_set(debug))
  expression = '\'+type+'TIME'
  time = mdsvalue(expression,quiet=quiet,status=status)
  if (not(status)) then begin
    time = mdsvalue('\EFIT_AEQDSK:TIME',quiet=quiet,status=status) ; C-Mod trees
  endif
  return,time
end

;-------------------------------------------------------------------------------------

function efit_read_mds,info,exact_time=exact_time,time_range=time_range, maxdt=maxdt, $
			     verbose=verbose, debug=debug, status=status, quiet=quiet  

  forward_function reada_mdsread, readg_mdsread, readm_mdsread
 
  ;;--- was a tree specified as a tag (ie. not "EFITnn") 
  if (strpos(strupcase(info.runid),'EFIT') ne 0) then begin
    tree = mdsvalue('EFITTREE($,$)',info.runid,info.shot,/quiet,status=stat)
    if (stat) then info.runid = tree ; if error, leave set to bogus runid so efit_read_mds will generate error
    efit_read_message,verbose,'EFIT_READ_SETINFO: Using EFIT runid: --->'+info.runid+'<---'
  endif

  mdsopen,info.runid,info.shot,quiet=quiet,status=status

  if (status) then begin
    efittime = efit_read_getmdstime(info.type,debug=debug,status=status)
    if (status) then begin
      itime = efit_read_mds_time(efittime,info.time,exact_time=exact_time,time_range=time_range,$
                                                  verbose=verbose,debug=debug,status=status,maxdt=maxdt)

      if (status) then begin

	;readX_mdsread functions in file readX.pro

	case (strupcase(info.type)) of
	  'A' : data = reada_mdsread(info, itime, verbose=verbose,debug=debug,status=status)
	  'G' : data = readg_mdsread(info, itime, verbose=verbose,debug=debug,status=status)  
	  'M' : data = readm_mdsread(info, itime, verbose=verbose,debug=debug,status=status)
	  else: begin
	    status = 0
	    data = efit_read_error()
	    efit_read_message,verbose,'EFIT_READ: Unrecognized file type: '+info.type
            if (keyword_set(debug)) then stop
    	  end
	endcase

      endif else begin
	data = efit_read_error()
   	if not(keyword_set(quiet)) then print,"EFIT_READ: Could not find suitable time in MDSplus for "$
                                              +strtrim(info.shot,2)+" "+strtrim(info.time,2)
      endelse

    endif else begin
      data = efit_read_error()
      efit_read_message,verbose,"EFIT_READ: Could not get MDSplus timebase \"+info.type+"TIME"
    endelse
  endif else begin
    data = efit_read_error()
    efit_read_message,verbose,"EFIT_READ:  Could not open tree "+info.runid+" "+strtrim(info.shot,2)
  endelse
	
  return,data
end

;-------------------------------------------------------------------------------------

function efit_read, arg1, arg2, type=type, mode=mode, runid=runid, info=info, source=source, $
		  		exact_time=exact_time, time_range=time_range, quiet=quiet, $
                                server=server, maxdt=maxdt, verbose=verbose, debug=debug, status=status

forward_function mdsplus_setup

  ;====== initialize keywords

  if (not(keyword_set(type))) then type='a' ; default to reading A file if no type specified
  if (not(keyword_set(verbose))) then verbose = 0
  if (keyword_set(debug)) then verbose = 1   ; if called in debug mode, then will always print messages


  ;====== Set up INFO structure if not passed in.
  ;====== Will contain shot, time, file, and EFIT run.

  calltype='READ'+strupcase(type)  ; used in efit_read_message messages

  if (not(keyword_set(info))) then begin
    if (not(keyword_set(arg1))) then begin
      status = 0
      efit_read_message,1,calltype+": Either call "+calltype+"(shot,time) or "+calltype+"(filename)"
      if (keyword_set(debug)) then stop
    endif else begin
      info=efit_read_setinfo(arg1,arg2,type=type,mode=mode,runid=runid,exact_time=exact_time, $
				       verbose=verbose,debug=debug,status=status)
      ;EXACT_TIME will be set to 1 if efit_read_setinfo discovers that ARG1 is a string (ie. specifies a file)
    endelse
  endif else status=1  ; no error checking on info if passed in!


  ;====== if able to assemble INFO structure, attempt to read from file, then MDSplus

  if (status) then begin  

    ;====== if mode is not restricted to MDSPLUS

    if (info.mode ne 'MDSPLUS') then begin
      data = efit_read_file(info,exact_time=exact_time,time_range=time_range, $
			     verbose=verbose, debug=debug, status=status)
    endif else begin

      status = 0  ; pretend that FILE read attempt was unsuccessful
      efit_read_message,verbose,"EFIT_READ: MDSplus mode selected"
      
    endelse
    
    
    ;====== If found file, then done.  Otherwise, try MDSplus if allowed.

    if (status) then begin

      source = 'FILE'
      efit_read_message,verbose,'EFIT_READ: Read from file '+info.file+' successful'

    endif else begin

      if (info.mode ne 'FILE') then begin
        
        status = mdsplus_setup(server=server)
        if (not(status)) then message,'Could not initialize MDSplus'
 
        data = efit_read_mds(info,exact_time=exact_time,time_range=time_range, quiet=quiet, $
                             verbose=verbose, debug=debug, status=status, maxdt=maxdt)
        
        if (status) then begin
          source = 'MDSPLUS'
        endif else begin
          efit_read_message,verbose,'EFIT_READ: MDSplus read unsuccessful'
        endelse
        
      endif else begin

        efit_read_message,verbose,'EFIT_READ: Error reading data from file '+info.file
	if (keyword_set(debug)) then stop
 
      endelse

    endelse


  endif else begin

    efit_read_message,verbose,"EFIT_READ: Either call "+calltype+"(shot,time) or "+calltype+"(filename)"
    if (keyword_set(debug)) then stop

  endelse

 
  if (status) then begin
    tmp_datatime = float(string(data.time,format='(F8.3)'))
    tmp_infotime = float(string(info.time,format='(F8.3)'))
    if (tmp_datatime ne tmp_infotime) then begin
      msg = 'WARNING: time requested = '+strtrim(info.time,2)+' does not match.  Time returned = '+strtrim(data.time,2)
      if not(keyword_set(quiet)) then print,msg
    endif
  endif else data = efit_read_error()

  return,data

end



