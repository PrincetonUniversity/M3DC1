;2000.01.07 - allow any file extension FNNNNNN.TTTTT* and convert only 
;             those of form FNNNNNN.TTTTT_TTT to sub-millisecond form
;             - Jeff Schachter

function efit_filename_parse_pad,array
  i = where(array eq 0,n)
  if (n gt 0) then array[i]=32b
  return,array
end

function efit_filename_parse,files

  ; purpose:  to get shot numbers and times from list of EFIT files.
  ; acceptable filenames are of form x012345.01234 or y987654.3210_123
 
  ; rather than looping over the list of files and using strmid etc, 
  ; I convert the string to a byte array and take advantage of the fact
  ; that when converted back to a string, the byte array is truncated at
  ; the first occurance of 0b for each row;

  ; for example,

  ; if files=['a012345.01234', 'b999999.99999_444'], then 
  ; byte(files) is a 17 x 2 array, where 17 is the length of the longest
  ; element in files.  Shorter elements are padded with 32b (space), which
  ; can be removed later (globally - no loop needed) with strtrim(x,2) 

  br=efit_filename_parse_pad(temporary(reverse(byte(files))))

  case (!VERSION.OS_FAMILY) of
   
    'unix' : begin

      ; First, strip off everything except what follows last / 
      ; This should be the filename itself.
      ; If no / in filename, everything's ok

      i=where(br eq (byte('/'))[0],n)  
      if (n gt 0) then begin
        br[i]=0b
      endif 

    end

    'vms' : begin
;      message,'VMS not yet supported'

      i = where(br eq (byte(']'))[0],n)
      if (n gt 0) then begin
        br[i]=0b
      endif 
      i = where(br eq (byte(':'))[0],n)
      if (n gt 0) then begin
        br[i]=0b
      endif 
    end

    ; must handle files having directory in name, or just having logical
    ; disk.  For example, both DISK:[DIRECTORY]X012345.01234 and 
    ; LOGICAL:Y999999.99999 must be handled.

    'MacOS' : begin

      ; First, strip off everything except what follows last / 
      ; This should be the filename itself.
      ; If no / in filename, everything's ok

      i=where(br eq (byte(':'))[0],n)  
      if (n gt 0) then begin
        br[i]=0b
      endif 

    end


  endcase


  ; remainder of operations are independent of operating system

  br = efit_filename_parse_pad(temporary(reverse(byte(string(br)))))  ; string(br) cuts strings where br=0b

;@@@ change for any file extension @@@;  ;; now replace _ (found in file extension) with . so can get time when take float of it
;@@@ change for any file extension @@@;  
;@@@ change for any file extension @@@;  i = where(br eq (byte('_'))[0], n)
;@@@ change for any file extension @@@;  if (n gt 0) then br(i) = (byte('.'))[0]


  filenames = strtrim(string(temporary(br)),2)

  ; now replace "nc" with -1 - signals time independent files
  i = where(strmid(filenames,8,2) eq 'nc',n)
  if (n gt 0) then begin
    br = byte(filenames)
    br[8,i] = (byte('-'))[0]
    br[9,i] = (byte('1'))[0]
    filenames = strtrim(string(temporary(br)),2)
  endif
    

  ;@@@ change for any file extension @@@;  REPLACE _ with . only if have _[0-9][0-9][0-9]
  ;if _ is present, it will be the 14th character in the string X000000.00000_000
  test = strmid(filenames,14,1) ; first char after _
  i=where(strpos(filenames,'_') eq 13 and test ge '0' and test le '9',n) ; have at least _[0-9]
  if (n gt 0) then begin
    test = strmid(filenames[i],14,max(strlen(filenames[i])))
    j = where(string(long(test),format='(i3.3)') eq test,n) ; find only _[0-9][0-9][0-9]
    if (n gt 0) then begin
      br = byte(filenames)
      br[13,i[j]] = 46b ; set to '.'
      filenames = strtrim(string(temporary(br)),2)
    endif
  endif
  ;@@@ change for any file extension @@@;  END CHANGE

  ; now fish out files that aren't of the form X000000.00000
  ; *** this is not a strict test... bogus filenames can slip by. ***
  ; Result will be that the long() operation on the shots and the double() on the times
  ; both return 0.

  test = byte(strmid(filenames,1,1))	; first look for files whose 2nd char isn't a number
  i = where(test lt 47 or test gt 57,n)
  if (n gt 0) then filenames(i) = "x-00001.-0001" ; so shot and time are returned as -1
  return,{shots:long(strmid(filenames,1,6)), $ 
;@@@ change for any file extension @@@;  times:double(strmid(filenames,8,max(strlen(filenames)))), $
	  times:strmid(filenames,8,max(strlen(filenames))), $
	  types:strmid(filenames,0,1), $
	  shottimes:strtrim(strmid(filenames,1,max(strlen(filenames))),2)}

end
