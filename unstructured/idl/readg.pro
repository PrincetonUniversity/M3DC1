;+ 
; NAME: 
;	READG
;
; PURPOSE: 
;
;	Retrieves GEQDSK data from file or MDSplus
;
; CATEGORY: 
;
;	DIII-D development 
;
; CALLING SEQUENCE: 
;
;	a = READG(arg1 [,arg2] [,MODE=mode] [,RUNID=runid] [,INFO=info] 
;                              [,SOURCE=source] [,EXACT_TIME=exact_time]
;			       [,VERBOSE=verbose] [,SERVER=server]
;                              [,DEBUG=debug] [,STATUS=status] )
;
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
;	(all are optional)
;	
;	MODE:  If "FILE", will restrict READG to retrieving EFIT data from
;	       files only, not from MDSplus.  If "MDSPLUS", will restrict 
;	       READG to retrieving EFIT data from MDSplus only, not files.  If
;	       not specified, READG will first attempt to retrieve the data 
;	       from a file, and then from MDSplus.
;	
;	RUNID:  EFIT "run ID" to use in MDSplus.  This defaults to "EFIT01" - 
;		the non-MSE automatic control room EFIT.
;	
;	INFO:  A structure with the following form:
;	
;		{mode:'', file:'', shot:0l, time:0.0d0, runid:''}
;	
;	       If specified as an input to READG, INFO will superceed the 
;	       arguments specified in arg1, arg2, and the keyword values of 
;              MODE and RUNID.
;	
;	       INFO is also returned from READG to indicate the values it used
;              to find the EFIT.
;	
;	SOURCE:  Either "FILE" or "MDSPLUS" - specifies the data source from
;	         where the EFIT data were retrieved.  Note that this information
;		 is also available in the returned structure as G.SOURCE.
;	
;	EXACT_TIME:  If set, forces READG to match the time specified, rather
;		     than using the default behavior of returning the nearest
;		     time.  
;
;	VERBOSE:  If set, READG will print out informational messages on its 
;                 progress.
;
;       SERVER:  If set to a string containing a valid IP address for
;                an MDSplus data server, READG will read the EFITs
;                from the specified server instead of the default for DIII-D.
;	
;	DEBUG:  If set, READG will print out additional debugging information,
;	        as well as turn off all error handling.  This will allow READG 
;	        to crash if there is an unexpected error.
;	
;	STATUS:  TRUE if READG was able to retrieve the data successfully, 
;                FALSE if not.  This information is also provided in the 
;		 output of the function (see below).
;
; OUTPUTS: 
;
;	Structure containing the GEQDSK data retrieved for the shot and
;	time specified.
;
;	Note that the structure contains the tag ERROR, which is 0 if the
;	data was read successfully, and 1 if not.
;
;	The returned structure also contains the tag SOURCE, a string that
;	describes the source from which the data was obtained (MDSplus or
;	File, which shot, EFIT run, and time).
;
; COMMON BLOCKS: 
;
;       COMMON EFIT_READG_CACHE,info_cache,data_cache
;
;	This common block caches GEQDSK data read from MDSplus.  READG reads 
;	the data for the entire time history, and then subscripts it at the 
;	time of interest. Subsequent references to data from the same shot and
;	EFIT run but different time will retrieve the data from the cache 
;	rather than reading it again from MDSplus.
;
; SIDE EFFECTS: 
;
;	Calls function EFIT_READ to handle read logic - as this logic is the 
;	same as for READA.
;
; RESTRICTIONS:
;
;	None.
;
; PROCEDURE: 
;
;	READG retrieves GEQDSK data from an EFIT run for a particular shot and
;	timeslice.  READG uses the following logic to locate the EFIT:
;	
;	- If arg1 specifies a file, READG attempts to determine the 6 digit
;	  shot number and time from the filename, assuming it has the format
;	  .../gSSSSSS.TTTTT_TTT.  _TTT is optional - used for specifying
;	  sub-millisecond timeslices.  If it cannot, it will still attempt to 
;	  read the file specified, but if the file attempt fails, the MDSplus
;	  attempt will also fail (see below).  NOTE THAT if arg1 specifies a
;	  file, READG will act as if the EXACT_TIME keyword is set - that is
;	  it will not attempt to find the nearest time if it cannot find the
;	  exact time.
;	
;	- If arg1 specifies the shot and arg2 the time, the filename
;	  gSSSSSS.TTTTT_TTT is composed, where SSSSSS is the 6 digit shot 
;	  number and TTTTT_TTT is the time, optionally using the _TTT for 
;	  sub-ms timeslices.
;	
;	- If the filename contains a directory specification, READG will look 
;	  for the file specified in the place specified.  If it does not, 
;	  READG will search  the following locations (in order) for the file:
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
;	  for an GEQDSK file with a time *nearest* the time specified.  
;	
;	- If the read attempt fails, or if the file is not found, READG will 
;	  attempt to read the data from MDSplus, using the shot number and 
;	  time specified (or determined from the filename).  Data from the 
;	  time *nearest* the time specified will be returned if the MDSplus 
;	  read attempt is successful (unless the keyword EXACT_TIME is set).
;	
;	- If the value of the keyword MODE is "MDSPLUS", READG will not
;	  attempt to read the data from a file, instead proceeding directly to
;	  the MDSplus read attempt. 
;	
;	- If the value ofthe keyword MODE is "FILE", READG will not attempt to
;	  read the data from MDSplus if the file read attempt is unsuccessful.
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
; DATE OF LAST MODIFICATION: 06/04/01
;
; MODIFICATION HISTORY:
;
;	Version 1.0: Released by Jeff Schachter 98.03.19
;	   98.03.26: bunch of bug fixes
;	Version 2.0: Standardized behavior in finding nearest time, and added keyword
;		     EXACT_TIME (Jeff Schachter, 98.04.27)
;	Version 2.1: 98.05.05 - Jeff Schachter - fix bug with unformatted reads
;	Version 2.2: 98.05.16 - Jeff Schachter -  added M.SOURCE to indicate from where
;		     data was obtained
;       Version 3.0: 98.08.21 - Jeff Schachter - return unpadded
;                               arrays when read from a file
;          98.08.25: fix bug with calculation of g.z
;       Version 3.1: 98.08.27 - Jeff Schachter - bdry array must be
;                               padded because number of bdry points
;                               fluctuates as a function of time
;       Version 3.2: Jeff Schachter 98.10.05
;                    - removed log of user path if /link/idl not present
;       Version 3.3: Jeff Schachter 1998.10.06
;                    - modify calls to MDSplus functions so that this
;                      procedure works with both client/server and
;                      native access
;	10-22-98 J.Ferron - added checking for program that creates the file.
;	11-06-98 Q.Peng - fixed a bug that failed files with 0 nbdry
;	12-17-98 Q.Peng - get time string by removing preceding and following
;		 chars instead of cutting 4 digits so that it can read 5-digit
;		 time correctly.
;	01-28-99 Q.Peng - if error when converting the time read from g file,
;		 fill the g.time with info.time - special fix when time field 
;		 in g file is ****.
; 	03-16-99 Q.Peng - moved the version checking for rotation and mass 
;		 from 970702 to 970916.
;	03-17-99 Q.Peng - fixed a bug so that it does not depend on delimitor
;		 when parsing time. Recognize 'EFIT*' instead of 'EFITD'.
;	05-06-99 Q.Peng - changed the logic for reading bdry and limiter 
;		 section so that it defaults to read as oppose to not to.
;		 It can now handle ONETWO eqdsk in addtion to
;		 efit,mfit,fixb.x.
;       1999.05.20: Jeff Schachter - ignore new G signal \PSIN (do not 
;                put in structure or try to subscript)
;	07-12-99 Q.Peng - pad bdry only when author_program is EFIT/MFIT so
;	         that creating G structure array still needs no extra effort. 
;	         But declare bdrypad dynamically if from other codes so if
;	         nbdry is large., e.g,ONETWO could be up to 3000, it can
;	         handle without increasing the size for general cases.
;		 Also reduce bdrymax from 600 to 400 (same as that in EFIT).
;	08-16-99 Q.Peng - added checking for author_program TEQ whose eqdsk
;	         stops after psirz. It does not have qpsi, bdry, limiter etc.
;	10-26-99 Q.Peng - Catch I/O error when reading the last part epoten
;		 to avoid crashing when a file does not have this info but 
;		 has namelists attached to the end.
;       10-27-99 J.Schachter - switch time from LONG to FLOAT (for
;                              sub-millisecond EFITs)
;       1999.11.02   Jeff Schachter - moved readg_mdsreadall to
;                                     efit_read_mdsreadall (general code)
; 	12-07-99 Q.P. Y2K compliant - Accommodate 4-digit year without 
;		breaking 2-digit year for version number
;       2000.07.28 J. Schachter - handle single timeslice EFITs in
;                                 readg_mdssub
;       2000.09.01 J. Schachter - make more robust against missing
;                                 data in MDSplus
;       2001.06.04 J. Schachter - add SERVER keyword
;	2002.07.03 Q. Peng - add NSTX
;	2003.11.04 Q. Peng - IDL6.0 returns (struct.tag)[i] as a size 1 array
;		rather than scalar as before. Convert to scalar explicitly. 
;	2005.05.18 Q. Peng - if failed to get the shot number from file 
;               (e.g. eqdsk generted by none-efit), use info.shot instead
;       2005.07.06 Q. Peng - catch and skip errors when creating the structure
;               from data_cache to skip special tags such as rgrid,zgrid from CMOD.
;               These tags are already dealt with separately.
;	2006.01.20 Ed Lazarus - process TRANSP(TRXPL) equilibrium header 
;               separately, incorporated by QP
;       2006.12.19 Q.P.- for machines other than D3D,CMOD,NSTX, add r,z,lim
;                  tags assuming the tree structure is similar to that of D3D.
;       2007.06.14 Q.P - handle ONETWO EQDSK header separately per Prater's request.
;                      - in readg_mdssub, when there is only one time slice
;                        convert size-1 arrays to scalars to avoid side-effects.

;-	


;-------------------------------------------------------------------------------------
; FILE SPECIFIC CODE
;-------------------------------------------------------------------------------------

function readg_file,info,verbose=verbose,debug=debug,status=status

  forward_function efit_read_error, efit_read_filesearch

  ;====== ONE ERROR HANDLER FOR ALL I/O ERRORS
  if (not(keyword_set(debug))) then catch,err else err=0
  if (err ne 0) then begin
    catch,/cancel
    efit_read_message,1,'READG_FILE: Error reading gfile '+info.file+': '+!ERR_STRING
    if (keyword_set(lun)) then free_lun,lun
    status = 0
    return,efit_read_error()
  endif

  ;====== open and read file

  limmax  = 300
  bdrymax = 3000  ; increased by request (L.Lao) to 3000 to handle efit x513 
    
;  g={GEQDSK,shot:long(0),time:long(0),error:long(0),source:info.file,$
;              ecase:strarr(6),mw:long(0),mh:long(0),xdim:float(0),$
;              zdim:float(0),rzero:float(0),rgrid1:float(0),zmid:float(0),$
;              rmaxis:float(0),zmaxis:float(0),ssimag:float(0),ssibry:float(0),$
;              bcentr:float(0),cpasma:float(0),$
;              fpol:fltarr(129),pres:fltarr(129),ffprim:fltarr(129),$
;              pprime:fltarr(129),psirz:fltarr(129,129),qpsi:fltarr(129),$
;              nbdry:long(0),limitr:long(0),bdry:fltarr(2,bdrymax),$
;              lim:fltarr(2,limmax),R:fltarr(129),Z:fltarr(129),$
;    	  rhovn:fltarr(129),epoten:fltarr(129)}
    
;  g.error = 1 ; will be reset later


;JMS 98.04.27.... get shot and time from geqdsk file rather than from info structure
;;;  g.shot = info.shot
;;;  g.time = info.time
    
    
  ;====== Create the variables for the first read.
    
  ecase = strarr(6)
  idum  = long(0)
  mw    = long(0)
  mh    = long(0)
  notes = ''    ; remaing of the 1st line, useful in EQDSK by codes like ONETWO

  ;====== Open the file for formatted reads.  
    
  openr,lun,info.file,/get_lun

    
  ; Read the first group of variables assuming that the file is formatted.
  ; If there is an error jump ahead to try to read the file assuming that
  ; it is unformatted.
    
  fileformatted = 0 ; initialize to "unformatted" - if first read is successful,
		    ; this will be switched to 1 - formatted

  on_ioerror,READG_UNFORMATTED

  readf, lun, ecase, idum, mw, mh, notes, format='(6a8,3i4,a)'

  ;==== if formatted read attempt successful, will reach here

  fileformatted = 1 
    
READG_UNFORMATTED:

  on_ioerror,NULL  ;==== turn off special error handling.  Now will be handled by error handler at top of procedure.

  if (not(fileformatted)) then begin    

    efit_read_message,verbose,'READG_FILE: File is not formatted.  Attempting unformatted read.'

    ; Close the file then reopen it for unformatted reads.

    close,lun
    openr,lun,info.file,/f77_unformatted,/segmented

    ; Recreate the variables for the first read.  MUST do for 
    ; unformatted read to work ok.

    ecase = strarr(6)
    ecase(*) = '          '
    bdum = bytarr(10,6)
    idum = long(0)
    mw = long(0)
    mh = long(0)
          
    ;readu,lun,ecase,idum,mw,mh	; this is commented out in old READG
    
    readu,lun,bdum,idum,mw,mh	; this is used in old READG
    ecase = string(bdum)

  endif else if (keyword_set(debug)) then print,'READG_FILE: File is formatted.'
    
  ;====  Now we have determined whether the file is formatted or unformatted.
  ;====  Procede to get the remainder of the data.

  ;====  Error check on data read so far

  ; At this point we should have the values of mw and mh.  A typical error
  ; seems to be that these are 0 for some as yet undetermined reason.
  
  if ( (mw eq 0) or (mh eq 0) or (mw gt 1000) or (mh gt 1000) ) then begin
    msg = "Problem reading G file.  MW = "+strtrim(mw,2)+"  MH = "+strtrim(mh,2)
    if (keyword_set(debug)) then begin
	print,msg
	stop
    endif else message,msg  ; generate error for handler
  endif
    
; Find out which code wrote this eqdsk file.
; Possible author programs are 
; EFITx - general efit codes, x being D,J,?
; MFITD - mfit
; feqdsk- fixb.x (a part of the toq suite of programs which labels itself 
;         feqdsk in ecase).
; g     - ONETWO wrteqdsk routine
; TEQ   - TEQ code from LLNL
; TRXPL - TRANSP equilibrium
; ONETWO- ONETWO EQDSK (header is different)
;
  IF strpos(notes,'ONETWO') GE 0 THEN author_program = 'ONETWO' ELSE $
  author_program = StrCompress(ecase(0),/Remove_All)

  ;==== proceed with read

  ;
  ; create the second group of variables.
  ;
    
  xdim=float(0)
  zdim=float(0)
  rzero=float(0)
  rgrid1=float(0)
  zmid = float(0)
  rmaxis = float(0)
  zmaxis = float(0)
  ssimag = float(0)
  ssibry = float(0)
  bcentr = float(0)   
  cpasma = float(0)
  xdum = float(0)
  xdum1 = float(0)
  fpol = fltarr(mw)
  pres = fltarr(mw)
  ffprim = fltarr(mw)
  pprime = fltarr(mw)
  psirz = fltarr(mw,mh)
  qpsi = fltarr(mw)
  
  ;
  ; Read the second group of variables.
  ;
    
  if (fileformatted) then begin
    readf,lun,xdim,zdim,rzero,rgrid1,zmid,format='(5e16.9)'
    readf,lun,rmaxis,zmaxis,ssimag,ssibry,bcentr,format='(5e16.9)'
    readf,lun,cpasma,ssimag,xdum,rmaxis,xdum,format='(5e16.9)'
    readf,lun,zmaxis,xdum,ssibry,xdum,xdum,format='(5e16.9)'
    readf,lun,fpol,format='(5e16.9)'
    readf,lun,pres,format='(5e16.9)'
    readf,lun,ffprim,format='(5e16.9)'
    readf,lun,pprime,format='(5e16.9)'
    readf,lun,psirz,format='(5e16.9)'
    if author_program NE "TEQ" then $
    readf,lun,qpsi,format='(5e16.9)'
  endif else begin
    readu,lun,xdim,zdim,rzero,rgrid1,zmid
    readu,lun,rmaxis,zmaxis,ssimag,ssibry,bcentr
    readu,lun,cpasma,ssimag,xdum,rmaxis,xdum
    readu,lun,zmaxis,xdum,ssibry,xdum,xdum
    readu,lun,fpol
    readu,lun,pres
    readu,lun,ffprim
    readu,lun,pprime
    readu,lun,psirz
    if author_program NE "TEQ" then $
    readu,lun,qpsi
  endelse
    
  ;
  ; Find the efit version number that was used to create the file.
  ; If it is large enough, then boundary and limiter data is present
  ; in the file so it can be read.
  ;
  ; Create the third group of variables.  These are needed even if the
  ; data cannot be read from the file.
  ;
    
  nbdry  = long(0)
  limitr = long(0)
  bdry = 0
  lim = 0

  ; Y2K compliant - Accommodate 4-digit year without breaking 2-digit 
  ; year for version number. QP 12-7-99
  
  ;nvernum = long(strmid(ecase(3-1),2-1,2)+strmid(ecase(2-1),4-1,2)+$
  ;	strmid(ecase(2-1),7-1,2))
  if author_program eq 'ONETWO' then nvernum = 0 else $
  nvernum = long(strcompress(strmid(ecase(3-1),2-1,4),/Remove)+$
	strmid(ecase(2-1),4-1,2)+strmid(ecase(2-1),7-1,2))
    
  if (keyword_set(debug)) then print,'READG_FILE: Version: '+strtrim(nvernum,2)+' '+info.file
    
  if (NOT ((StrMid(author_program,0,4) EQ "EFIT" $
	    OR author_program EQ 'MFITD' ) AND (nvernum lt 870520)) AND $
      NOT ((author_program EQ "feqdsk") AND (nvernum lt 980709)) AND $
      NOT (author_program EQ "TEQ") ) then begin
    ;
    ; Read the third group of variables.
    ;
    if (fileformatted) then begin
      readf,lun,nbdry,limitr,format='(2i5)'
    endif else begin
      readu,lun,nbdry,limitr
    endelse
    ;
    ; Make the fourth group of variables.
    ;
    if nbdry ne 0 then   bdry = fltarr(2,nbdry)
    lim = fltarr(2,limitr)
    ;
    ; Read the fourth group of variables.
    ;
    chardm = ''
    if (fileformatted) then begin
      if nbdry ne 0 then readf,lun,bdry,format='(5e16.9)' else  readf,lun, chardm
      readf,lun,lim,format='(5e16.9)'
    endif else begin
      readu,lun,bdry
      readu,lun,lim
    endelse  
  endif
    
  ; 
  ; If efit version number is >= 970702, read electrostatic potential funciton
  ; after skipping rotation information and ion mass density profile.
  ; The version for this option seems no earlier than 970916. -QP 03-16-99
  ;
    
  rhovn=fltarr(mw)
  epoten=fltarr(mw)
  
  if((StrMid(author_program,0,4) EQ "EFIT") AND (nvernum ge 970916)) then BEGIN
;     the feqdsk "author program" does not contain this information.
    ; 
    ; Read and skip rotation information if available
    ;
    kvtor=long(0)
    rvtor=float(0)
    nmass=long(0)
    if (fileformatted) then begin
      readf,lun,kvtor,rvtor,nmass,format='(i5,e16.9,i5)' 
    endif else readu,lun,kvtor,rvtor,nmass  
    if (kvtor gt 0) then begin
      pressw=fltarr(mw)
      pwprim=fltarr(mw)
      if (fileformatted) then begin
        readf,lun,pressw,format='(5e16.9)'
        readf,lun,pwprim,format='(5e16.9)'
      endif else begin
        readu,lun,pressw
        readu,lun,pwprim
      endelse
    endif
  ;
  ; Read ion mass density profile if available
  ;
    if (nmass gt 0) then begin
      dmion=fltarr(mw)
      if (fileformatted) then begin
        readf,lun,dmion,format='(5e16.9)' 
      endif else readu,lun,dmion
    endif
  ;
  ; Read electrostatic potential function if available
  ; Error may occur if namelists attached to the end are mistaken as 
  ; the data. Skip instead crash in this case. QP 10-26-99
  ;

on_ioerror,nodata	
    keecur=long(0)
    if (fileformatted) then begin
      readf,lun,rhovn,format='(5e16.9)'
      readf,lun,keecur,format='(i5)'
      if(keecur gt 0) then readf,lun,epoten,format='(5e16.9)'
    endif else begin
      readu,lun,rhovn
      readu,lun,keecur
      if(keecur gt 0) then readu,lun,epoten
    endelse  
  endif
nodata:
on_ioerror,null
    
      ;
      ; NOTE: at this point the header information is not implemented on 
      ; non-VAX computers.
      ;
      ;header = string(replicate(32b,42))
      ;print,'READG:  header ="'+header+'"'
      ;if(fileformatted) then begin
      ;   readf,lun,header,format='(a43)'
      ;endif else begin
      ;   readu,lun,header
      ;endelse
      ;print,'READG:  header ="'+header+'"'
      ;
      ; Copy the variables into the structure.
      ;
;  g.ecase = ecase
;  g.mw = mw
;  g.mh = mh
;  g.xdim=xdim
;  g.zdim=zdim
;  g.rzero=rzero
;  g.rgrid1=rgrid1
;  g.zmid=zmid
;  g.rmaxis=rmaxis
;  g.zmaxis=zmaxis
;  g.ssimag=ssimag 
;  g.ssibry=ssibry 
;  g.bcentr=bcentr 
;  g.cpasma=cpasma 
;  g.fpol=fpol
;  g.pres=pres 
;  g.ffprim=ffprim 
;  g.pprime=pprime 
;  for i=0,mh-1 do g.psirz(0:mw-1,i) = psirz(0:mw-1,i)
;  g.qpsi=qpsi 
;  g.nbdry=nbdry 
;  g.limitr=limitr 
;  if (nbdry gt 0) then begin
;    if (nbdry le bdrymax) then for i=0,nbdry-1 do g.bdry(0:1,i)=bdry(0:1,i)
;  endif
;  if (limitr gt 0) then begin
;    if (limitr le limmax) then for i=0,limitr-1 do g.lim(0:1,i)=lim(0:1,i)
;  endif
;  g.rhovn=rhovn
;  g.epoten=epoten
  ;
  ;************** added g.r,g.z calc and changed signs on flux--B. Rice
  ;
  dR = xdim/(mw-1)
  dz = zdim/(mh-1)
;  for i=0,mw-1 do begin
;    g.r(i) = rgrid1+i*dR
;  endfor
;  for i=0,mh-1 do begin
;    g.z(i) = zmid-0.5*zdim+i*dz
;  endfor
  ;
  ; ip_sign= g.cpasma/abs(g.cpasma)
  ; g.ssimag= -g.ssimag*ip_sign ;change signs here so Bz comes out with right sign
  ; g.ssibry= -g.ssibry*ip_sign
  ; g.psirz= -g.psirz*ip_sign 
  ;
  ; Get the actual shot number and time from the G file data.
  ;
  if author_program eq 'ONETWO' then begin
   ; ONETWO equilibrium
     shotuse = long(ecase(0));
     nntime = ecase(1);
     if info.time eq 0 then info.time = float(nntime)
  endif else if strpos(author_program,'TRXPL') lt 0 then begin
   ; regular EFIT equilibrium
     nine = strmid(ecase(4-1),3-1,1)
     if(nine ne ' ') then begin
       nshot6=strmid(ecase(4-1),3-1,6)
;      g.shot = long(nshot6)
       shotuse = long(nshot6)
     endif else begin
       nshot5=strmid(ecase(4-1),4-1,5)
;      g.shot = long(nshot5)
       shotuse = long(nshot5)
     endelse
;    050518 QP - if failed to figure out the shot number, use info.shot instead
     if shotuse eq 0 then shotuse = info.shot 
;    nntime = strmid(ecase(5-1),3-1,4)
     nntime = str_sep(ecase[4],',')	; remove preceding ","
     nntime = nntime[N_Elements(nntime)-1]
     nntime = (str_sep(nntime,'ms'))[0]	; remove following "ms"
  endif else begin
   ; process TRANSP(TRXPL) equilibrium separately
     alla = ecase(0) & for ipl=1,5 do alla = alla+ecase(ipl)
     alla = str_sep(alla," ")
     transpid = alla[2]
     transptime=float(alla(6))
     transptime = long(transptime*1000)
     alphabet = ['a','b','c','d','e','f',$
                 'g','h','i','j','k','l','m','n','o','p',$
                 'q','r','s','t','u','v','w','x','y','z']
     alphabet = [alphabet,'A','B','C','D','E','F',$
                 'G','H','I','J','K','L','M','N','O','P',$
                 'Q','R','S','T','U','V','W','X','Y','Z']
     l1 = strlen(transpid)
     l2 = strarr(l1)
     for i = 0,l1-1 do l2(i) = strmid(transpid,i,1)
     l3 = intarr(l1)
     for i = 0,l1-1 do l3(i) = where (l2(i) eq alphabet)
     l4 = where(l3 ne -1)
     nshot6 = strmid(transpid,0,l4)
     shotuse = long(nshot6)
     nntime = string(transptime)
     alphabet = 0
     alla = 0
     l1 = 0 & l2 = 0 & l3 = 0 & l4 = 0
  endelse

;;;=========================
;;; DETERMINE TIME EITHER FROM G FILE OR FROM FILENAME
;;; G FILE DOES *NOT* STORE SUB-MILLISECOND TIME
;;; - revised by Jeff Schachter 1999.10.27
  gtimeFile = 0.
  gtime = info.time 			; use time from info - overwrite next.
  on_ioerror,bad_time			; catch convertion error for time
  gtimeFile = float(nntime)	        ; overwrite - use time read from file  ; change from long to float 99.10.27
bad_time:				; error - usually '****'
  on_ioerror,null			; turn off the special handling.
  if (gtime ne gtimeFile and (long(gtime) eq gtime)) then begin
    message,/info,'Trouble determining exact time.  Filename gives: '+strtrim(gtime,2)+'  but file contains: '+strtrim(gtimeFile,2)
  endif

  ;;; G structure returned will use gtime, so it will be using the
  ;;; time determined from the file name.
;;;=========================
  
      ;
      ; The routine executed successfully if the boundary or limiter arrays
      ; did not overflow.
      ;
;  if( (nbdry le bdrymax) and (limitr le limmax) ) then g.error = 0
      ;
      ; close the file
      ;

  bdrypad = 0
  if (StrMid(author_program,1,3) EQ "FIT") $
  then bdrypad = fltarr(2,bdrymax) $	; pad only for EFIT/MFIT g files.
  else if nbdry ne 0 then bdrypad = fltarr(2,nbdry) ; dynamic for other codes.
  if n_elements(bdry) ne 0 then if nbdry ne 0 then $
     bdrypad[*,0:nbdry-1] = bdry

  code_version = strtrim(ecase[1]+ecase[2],2)
  g={shot:shotuse, time:gtime, code_version:code_version, error:0l, source:info.file,$
              ecase:ecase, mw:mw, mh:mh, xdim:xdim, $ 
              zdim:zdim, rzero:rzero, rgrid1:rgrid1, zmid:zmid, $
              rmaxis:rmaxis, zmaxis:zmaxis, ssimag:ssimag, ssibry:ssibry, $
              bcentr:bcentr, cpasma:cpasma, $
              fpol:fpol, pres:pres, ffprim:ffprim, $
              pprime:pprime, psirz:psirz, qpsi:qpsi, $
              nbdry:nbdry, limitr:limitr, bdry:bdrypad, $
              lim:lim, r:rgrid1+findgen(mw)*dR, z:zmid-0.5*zdim+findgen(mh)*dz, $
              rhovn:rhovn, epoten:epoten}
   

  free_lun,lun
  status = 1

  return,g


end

;-------------------------------------------------------------------------------------
;  MDSplus SPECIFIC CODE
;-------------------------------------------------------------------------------------



function readg_mdssub,shot,itime,runid,status=status

  common efit_readg_cache,info_cache,data_cache


  source = 'MDSplus, shot = '+strtrim(shot,2)+', run = '+runid+', time = '+strtrim(data_cache.time[itime[0]],2)

  ;;; 2007JUN15 QP
  ;;; Instead of returning the whole structure when there is only one time slice
  ;;; create a new struct and convert arrays of size 1 to scalars since
  ;;; size-1-array in place of scalar causes side-effect in postprocessing.
  if (n_elements(data_cache.time) eq 1) then begin
    data = create_struct('tree',runid,'source',source)
    tags=tag_names(data_cache)
    for i=0,n_tags(data_cache)-1 do begin
	if (size(data_cache.(i),/n_dim) eq 1) and $     ; test for size-1-array
	   (size(data_cache.(i),/n_ele) eq 1) then $
	data = create_struct(temporary(data),tags[i],(data_cache.(i))[0]) else $
	data = create_struct(temporary(data),tags[i],data_cache.(i))
    end
    return,data
  endif

  data={shot:shot, time:data_cache.time[itime], tree:runid, error:0, source:source}

  tags=tag_names(data_cache)


  i = where(tags eq 'LIM',n)
  machine = mdsvalue('MACHINE()')
  case machine of
  'D3D': begin
    ;;; D3D
    dummy = where(tags eq 'R' or tags eq 'Z',n)
    if (n ne 2) then begin
      data.error = 1
      status = 0
      return,data
    endif
    data = create_struct(data,'lim',data_cache.lim,'r',data_cache.r,'z',data_cache.z)
   end
   'CMOD': begin
    ;;; CMOD
    s = size(data_cache.rbbbs)
    bdry = fltarr(2,s[1],s[2])
    bdry[0,*,*] = data_cache.rbbbs
    bdry[1,*,*] = data_cache.zbbbs
    data = create_struct(data,'lim',transpose([[data_cache.xlim], [data_cache.ylim]]), $
                              'r',data_cache.rgrid, $
                              'z',data_cache.zgrid);, $
                              ;'bdry',temporary(bdry), $
                              ;'rgrid1',data_cache.rgrid[0])
   end
   'NSTX': begin
    s = size(data_cache.rbdry)
    bdry = fltarr(2,s[1])
    bdry[0,*] = (data_cache.rbdry)[*,itime]
    bdry[1,*] = (data_cache.zbdry)[*,itime]
    data = create_struct(data,$
	'lim',transpose([[data_cache.xlim[*,0]],[data_cache.ylim[*,0]]]),$
	'r',(data_cache.r)[*,itime],$
	'z',(data_cache.z)[*,itime],$
	'bdry',temporary(bdry))
   end
   else: begin
    ; other machine, assume same structure as D3D
    dummy = where(tags eq 'R' or tags eq 'Z',n)
    if (n ne 2) then begin
      message,machine+' - no R, Z in MDSplus',/info
      data.error = 1
      status = 0
      return,data
    endif
    data = create_struct(data,'lim',data_cache.lim,'r',data_cache.r,'z',data_cache.z)
   end
   endcase


  if (where(tags eq 'NBDRY'))[0] eq -1 then begin ; NSTX has both NBDRY, NBBBS
    i = where(tags eq 'NBBBS',n)
    if (n eq 1) then tags[i] = 'NBDRY' ; rename CMOD for compatibility with D3D
  endif

  ix=where(tags ne 'SHOT' and tags ne 'TIME' and tags ne 'R' and tags ne 'Z' $
	and tags ne 'LIM' and tags ne 'ERROR' and tags ne 'PSIN', nx)

  if (nx gt 0) then begin
    for i=0,nx-1 do begin
      j=ix[i]
      ; 050706 QP - catch and skip problematic tags such as rgrid, zgrid from CMOD
      ; if unexpected tags are skipped, should print and check !ERR_STRING
      if (not(keyword_set(debug))) then catch,err else err=0
      if (err ne 0) then begin
	 err = 0
         efit_read_message,1,'skipping '+tags[j]  ;+': '+!ERR_STRING
      	 data = temporary(datacopy)
	 goto,skip_tag
      endif
      datacopy = data
      case ((size(data_cache.(j)))[0]) of 
	0 : data=create_struct(temporary(data),tags[j],data_cache.(j))
        1 : data=create_struct(temporary(data),tags[j],(data_cache.(j)[itime])[0])
        2 : data=create_struct(temporary(data),tags[j],data_cache.(j)[*,itime])
        3 : data=create_struct(temporary(data),tags[j],data_cache.(j)[*,*,itime])
      endcase
skip_tag:
      catch,/cancel
    endfor
    status  = 1
  endif else begin
    status = 0
    data.error = 1
  endelse

  return,data
end

function readg_mdsread,info,itime,verbose=verbose,debug=debug,status=status

  common efit_readg_cache,info_cache,data_cache
  forward_function efit_read_error, efit_read_mdsreadall

  ;====== initialize cache

  if (not(keyword_set(info_cache))) then begin
    efit_read_message,verbose,'READG_MDSREAD: Initializing info_cache'
    info_cache={shot:0l, runid:''}
  endif


  ;====== if data already cached

  if (keyword_set(data_cache) and (info.shot eq info_cache.shot) $
	                      and (info.runid eq info_cache.runid)) then begin

    data = readg_mdssub(info.shot,itime,info.runid,status=status)
    efit_read_message,verbose,'READG_MDSREAD: MDSplus data was cached already'
    
  endif else begin

    ;====== read data and cache it

    data_cache=efit_read_mdsreadall(info,verbose=verbose,debug=debug,status=status) 

    if (status) then begin

      info_cache=info
      data=readg_mdssub(info.shot,itime,info.runid,status=status)
      if (status) then $
          efit_read_message,verbose,'READG_MDSREAD: MDSplus read from '+info.runid+' '+strtrim(info.shot,2)+' successful.'

    endif else begin

      data_cache = '' ; unset data_cache variable
      data = efit_read_error()
      efit_read_message,verbose,'READG_MDSREAD: Error reading MDSplus GEQDSK data '+strtrim(info.shot,2)+' '+info.runid+' '+strtrim(info.time,2)
      if (keyword_set(debug)) then stop

    endelse

  endelse

  return,data
end

      
;-------------------------------------------------------------------------------------

function readg, arg1, arg2, mode=mode, runid=runid, info=info, source=source, $
			    exact_time=exact_time, time_range=time_range, $
                            server=server, quiet=quiet, maxdt=maxdt, $
			    verbose=verbose, debug=debug, status=status 

  forward_function efit_read, efitde_read

  if keyword_set(quiet) then verbose=0

  if n_elements(runid) gt 0 then begin
    if strupcase(strtrim(strjoin(runid),2)) eq 'EFITDE' then begin
      g = efitde_read(arg1,type='g',debug=debug,status=status,server=server)
      return,g
    endif
  endif
  
  g = efit_read(arg1, arg2, type='g', mode=mode, runid=runid, $
                info=info, source=source, $
                exact_time=exact_time, time_range=time_range, $
                server=server, maxdt=maxdt, quiet=quiet, $
                verbose=verbose, debug=debug, status=status)
  
  return,g

end

