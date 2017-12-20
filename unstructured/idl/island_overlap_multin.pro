; filename[i,j] : first index is toroidal mode number, second is 
; current spans rows
function island_overlap_multin, filename, current=cur, plot=plot0, ntor=ntor, $
                                phi0=phi0
   if(n_elements(phi0) eq 0) then phi0 = 0.
   nmodes = n_elements(filename[*,0])
   nrows = n_elements(filename[0,*])

   width = fltarr(nmodes, 20)
   psi = fltarr(nmodes, 20)
   q = fltarr(nmodes, 20)
   nres = intarr(nmodes)

   if(keyword_set(plot0)) then plot, [0.0, 1], [0, 5], /nodata
   for i=0, nmodes-1 do begin
       scur = cur*exp(-complex(0.,1.)*ntor[i]*phi0)
       w = island_widths(reform(filename[i,*]), $
                         psin=psin, current=scur, q=qval, /sum)
       if((n_elements(w) eq 1)) then if(w eq 0) then continue

       nres[i] = n_elements(psin[*])
       width[i,0:nres[i]-1] = w[0,*]
       psi[i,0:nres[i]-1] = psin[0,*]
       q[i,0:nres[i]-1] = qval[0,*]

       if(keyword_set(plot0)) then begin
           for j=0, nres[i]-1 do begin
               oplot, [psi[i,j]-width[i,j]/2., psi[i,j]+width[i,j]/2.], $
                 [q[i,j], q[i,j]]+(i-1)/20., color=color(i)
               oplot, [psi[i,j], psi[i,j]], !y.crange, color=color(i), $
                 linestyle=1
           end
       end
   end

   ; determine island overlap width
   pos = 1.
   extended = 1
   while(extended eq 1) do begin
       extended = 0
       for i=0, nmodes-1 do begin
           for j=0, nres[i]-1 do begin
               pmax = psi[i,j]+width[i,j]/2.
               pmin = psi[i,j]-width[i,j]/2.
               if((pmax ge pos) and (pmin lt pos)) then begin
                   extended=1
                   pos=pmin
               end
           end
       end
   end
   
   if(keyword_set(plot0)) then oplot, [pos,pos], !y.crange

   if((1.-pos) gt 1.) then stop

   return, 1.-pos
end
