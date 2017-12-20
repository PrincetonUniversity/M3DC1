function make_label, s, d, _EXTRA=extra
   if(n_elements(d) eq 0) then d = dimensions(_EXTRA=extra)
   if(max(d,/abs) eq 0.) then return, s
   return, s + ' !6(!X' + parse_units(d, _EXTRA=extra) + '!6)!X'
end
