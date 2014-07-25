;======================================================================
; make_units
; ~~~~~~~~~~
;
; creates a string given the specified unit dimensions
;======================================================================
function make_units, _EXTRA=extra
   return, parse_units(dimensions(_EXTRA=extra), _EXTRA=extra)
end
