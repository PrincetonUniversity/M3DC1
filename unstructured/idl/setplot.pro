pro setplot, p
   if(1 eq strcmp(p, 'ps', /fold_case)) then begin
       !p.charthick=5
       !p.thick=5
       !p.charsize=1.5
   endif else if (1 eq strcmp(p, 'x', /fold_case)) then begin
       !p.charthick=1
       !p.thick=1
       !p.charsize=1.5
   endif

   set_plot, p
end
