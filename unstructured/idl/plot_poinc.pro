pro plot_poinc, folder, overplot=overplot, color=color, skip=skip, start=start

  Ns = n_elements(file_search(folder+'/out*'))
  if n_elements(skip) eq 0 then skip = 1
  if n_elements(start) eq 0 then start = 0

  col = get_colors(Ns+2)
  if n_elements(color) ne 0 then begin
    csize = size(color)
    if (csize[1] eq 7) and (color eq 'k' or color eq 'w') then begin
      tvlct,r,g,b,/get
      if color eq 'k' then set_python_ct,100 else set_python_ct,100,/flip_bw
      color = color(0)
    endif
    col = 0.*col + color
  endif
  
  if not keyword_set(overplot) then begin
    Xmax = -1e20
    Zmax = -1e20
    Xmin = 1e20
    Zmin = 1e20
    for i=start,Ns-1,skip do begin
      P = read_ascii(folder+'/out'+strtrim(string(i),2))
      X = P.field1[1,*]
      Z = P.field1[2,*]
      if max(X) gt Xmax then Xmax = max(X)
      if max(Z) gt Zmax then Zmax = max(Z)
      if min(X) lt Xmin then Xmin = min(X)
      if min(Z) lt Zmin then Zmin = min(Z)
    endfor
    plot,[Xmin,Xmax],[Zmin,Zmax],/nodata
  endif
  
  for i=start,Ns-1,skip do begin
    P = read_ascii(folder+'/out'+strtrim(string(i),2))
    X = P.field1[1,*]
    Z = P.field1[2,*]
    oplot,X,Z,psym=3,color=col[i+1],symsize=0.1
  endfor
  
  if n_elements(g) ne 0 then tvlct,r,g,b
  
end