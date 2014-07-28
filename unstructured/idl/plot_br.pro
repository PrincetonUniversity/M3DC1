pro plot_br, _EXTRA=extra, bins=bins, q_val=q_val, $
             subtract_vacuum=subtract_vacuum, ntor=ntor, $
             plotbn=plotbn, slice=slice, extsubtract=extsubtract, $
             overplot=overplot

   if(n_elements(ntor) eq 0) then begin
       ntor = read_parameter('ntor', _EXTRA=extra)
   endif
   if(n_elements(extsubtract) eq 0) then begin
       extsubtract = read_parameter('extsubtract' ,_EXTRA=extra)
   end
   print, 'ntor = ', ntor
   if(n_elements(slice) eq 0) then last=1

;   psi0 = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra)
;   psi0_r = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra,op=2)
;   psi0_z = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra,op=3)
;   i0   = read_field('i'  ,x,z,t,/equilibrium,_EXTRA=extra)

   threed = read_parameter('3d', _EXTRA=extra)
   print, '3D = ', threed

   psi0 = read_field('psi',x,z,t,slice=-1,_EXTRA=extra)
   psi0_r = read_field('psi',x,z,t,slice=-1,_EXTRA=extra,op=2)
   psi0_z = read_field('psi',x,z,t,slice=-1,_EXTRA=extra,op=3)
   i0   = read_field('i'  ,x,z,t,slice=-1,_EXTRA=extra)

   if(threed eq 1) then begin
      bx = read_field_3d('bx',phi,x,z,t,last=last,slice=slice, $
                      /linear,_EXTRA=extra, ntor=ntor)
      bz = read_field_3d('bz',phi,x,z,t,last=last,slice=slice, $
                      /linear,_EXTRA=extra, ntor=ntor)
   endif else begin
      bx = read_field('bx',x,z,t,last=last,slice=slice, $
                      /linear,_EXTRA=extra,/complex)
      bz = read_field('bz',x,z,t,last=last,slice=slice, $
                      /linear,_EXTRA=extra,/complex)
   endelse

   if(keyword_set(subtract_vacuum)) then begin
       if(extsubtract eq 0) then begin
           bx0 = read_field('bx',x,z,t,slice=0,/linear,_EXTRA=extra)
           bz0 = read_field('bz',x,z,t,slice=0,/linear,_EXTRA=extra)
           bx = bx - bx0
           bz = bz - bz0
       end
   endif


   r = radius_matrix(x,z,t)
   y = z_matrix(x,z,t)

   br = -(bx*psi0_r + bz*psi0_z) / sqrt(psi0_r^2 + psi0_z^2)

   ; convert to cgs
   get_normalizations, b0=b0_norm, n0=n0_norm, l0=l0_norm, _EXTRA=extra
   br = br*b0_norm

   schaffer_plot, br, x, z, t, _EXTRA=extra, psi0=psi0,i0=i0, q_val=q_val, $
     label='!8B!Dn!N!6 (G)!X', points=points, bins=bins, ntor=ntor, $
     overplot=overplot
end
