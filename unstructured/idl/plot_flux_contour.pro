pro plot_flux_contour, fval, _EXTRA=extra
   n = n_elements(fval)

   psi = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra)

   loadct, 12
   contour, psi, x, z, closed=0, levels=fval(sort(fval)), color=color(6,10), $
     _EXTRA=extra
end
