function get_colors, maxcolors

   if(n_elements(maxcolors) eq 0) then maxcolors=10
   mc = max([maxcolors,10])

   c = indgen(mc) * (!d.table_size-25) / mc + 25

   if (1 EQ strcmp(!d.name, 'PS')) then begin
       c[0] = 0
   endif else begin
       c[0] = !d.table_size-1
   endelse 
   
   return, c
end
