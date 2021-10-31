function island_widths, filename, psin=psin, current=cur, q=q, $
                        sum_files=sum_files, netcdf=netcdf, plot=pl, $
                        ntor=ntor, _EXTRA=extra

    result = read_bmn(filename, m, bmn, phase, psin=psin, qval=q, ntor=ntor, $
                      qprime=qprime, area=area, psiprime=psiprime, $
                      sum_files=sum_files, factor=cur, netcdf=netcdf)

    if(result eq 1) then return, 0

    bmn = bmn/1e4               ; convert to tesla
    width = bmn*0.

    n = n_elements(bmn[*,0])

    if(keyword_set(pl)) then begin
       plot, psin[0,*], abs(q), _EXTRA=extra
    end

    for i=0, n-1 do begin
        width[i,*] = (2./!pi)*sqrt((area[i,*]*bmn[i,*]/abs(m*psiprime[i,*])) $
                                   *abs(q[i,*]/qprime[i,*]))
        if(keyword_set(pl)) then begin
           for j=0, n_elements(q[i,*])-1 do begin
              oplot, [-.5,.5]*width[i,j] + psin[i,j], [1.,1.]*abs(q[i,j])
           end
        end
    end

    j = where(width ne width, count)
    if(count ne 0) then stop

    return, width
end
