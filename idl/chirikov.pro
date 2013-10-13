function island_widths, filename, psin=psin, current=cur, q=q, $
                        sum_files=sum_files

    result = read_bmn(filename, m, bmn, phase, psin=psin, qval=q, $
                      qprime=qprime, area=area, psiprime=psiprime, $
                      sum_files=sum_files, factor=cur)

    if(result eq 1) then return, 0

    bmn = bmn/1e4               ; convert to tesla
    width = bmn*0.

    n = n_elements(bmn[*,0])

    for i=0, n-1 do begin
        width[i,*] = (2./!pi)*sqrt((area[i,*]*bmn[i,*]/abs(m*psiprime[i,*])) $
                                   *abs(q[i,*]/qprime[i,*]))
    end

    j = where(width ne width, count)
    if(count ne 0) then stop

    return, width
end

function chirikov, filename, psin=psin, psimid=psimid, current=cur, $
                   sum_files=sum_files
   width = island_widths(filename, psin=psin, current=cur, sum_files=sum_files)

   chi = fltarr(n_elements(width[*,0]), n_elements(psin[0,*])-1)
   psimid = fltarr(n_elements(psin[*,0]), n_elements(psin[0,*])-1)

   for j=0, n_elements(width[0,*])-2 do begin
       chi[*,j] = (width[*,j+1]+width[*,j])/2. $
         / (psin[*,j+1]-psin[*,j])
       psimid[*,j] = (psin[*,j+1]+psin[*,j])/2.
   end

   return, chi
end

function chi95, filename, psi0=psi0, current=cur, sum_files=sum_files

   if(n_elements(psi0) eq 0) then psi0 = 0.95
   if(n_elements(psi0) lt n_elements(filename)) then $
     psi0 = replicate(psi0, n_elements(filename))

   chi = chirikov(filename, psin=psin, psimid=psimid, current=cur, $
                  sum_files=sum_files)

   m = n_elements(chi[0,*])
;   if(keyword_set(sum_files)) then n = 1 else n =
;   n_elements(filename)
   n = n_elements(chi[*,0])
   c95 = fltarr(n)
   for i=0, n-1 do begin
       c95[i] = interpol(chi[i,0:m-1], psimid[i,0:m-1], psi0[i])
   end

   return, c95
end

function island_overlap, filename, current=cur, sum_files=sum_files
   width = island_widths(filename, psin=psin,current=cur,sum_files=sum_files,$
                        q=q)

   n = n_elements(psin[*,0])
   overlap = fltarr(n)

   for j=0, n-1 do begin
       lpl = 1.
       for i=n_elements(psin[j,*])-1, 0, -1 do begin
           pr = psin[j,i] + width[j,i]/2.
           pl = psin[j,i] - width[j,i]/2.

           if(pr lt lpl) then break
           lpl = pl
       end
       overlap[j] = 1.-lpl
   end

   return, overlap
end

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
