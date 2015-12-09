pro grad_shafranov, _EXTRA=extra, nterms=nterms, terms=term, names=names, $
               title=title, x=x, z=z
  title = 'Grad Shafranov Equation'
  nterms = 3
  names = ['del*(psi)',"R^2 p'","FF'"]
  
  psi_r = read_field('psi',x,z,t,_EXTRA=extra,/equilibrium,op=2)
  psi_z = read_field('psi',x,z,t,_EXTRA=extra,/equilibrium,op=3)
  p_r = read_field('p',x,z,t,_EXTRA=extra,/equilibrium,op=2)
  p_z = read_field('p',x,z,t,_EXTRA=extra,/equilibrium,op=3)
  i_r = read_field('i',x,z,t,_EXTRA=extra,/equilibrium,op=2)
  i_z = read_field('i',x,z,t,_EXTRA=extra,/equilibrium,op=3)
  i = read_field('i',x,z,t,_EXTRA=extra,/equilibrium)
  jy = read_field('jy',x,z,t,_EXTRA=extra,/equilibrium)

  r = radius_matrix(x,z,t)

  d = size(i,/dim)
  term = complexarr(nterms,d[1],d[2])
     
  ; dV / dt
  print, 'defining del*(psi)'
  term[0,*,*] = -r*jy

  ; -JxB
  print, "defining p'"
  term[1,*,*] = r^2 * (p_r*psi_r + p_z*psi_z) / (psi_r^2 + psi_z^2)

  ; Grad(p)
  print, "defining FF'"
  term[2,*,*] = i*(i_r*psi_r + i_z*psi_z) / (psi_r^2 + psi_z^2)
end
