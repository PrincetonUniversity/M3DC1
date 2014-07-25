; Returns the name of timeslice t
function time_name, t
   if(t lt 0) then begin
       label = "equilibrium"
   endif else begin
       label = string(FORMAT='("time_",I3.3)', t)
   endelse
   return, label
end
