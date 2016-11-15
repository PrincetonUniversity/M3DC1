pro plot_br, _EXTRA=extra, bins=bins, q_val=q_val, $
             subtract_vacuum=subtract_vacuum, ntor=ntor, $
             slice=slice, extsubtract=extsubtract, $
             overplot=overplot, filename=filename, scale=scale, $
             linfac=linfac, sum=sum, vacfield=vacfield, velocity=velocity

   if(n_elements(filename) eq 0) then filename='C1.h5'
   if(n_elements(ntor) eq 0) then begin
       ntor = read_parameter('ntor', _EXTRA=extra, filename=filename[0])
   endif
   if(n_elements(extsubtract) eq 0) then begin
       extsubtract = read_parameter('extsubtract' ,_EXTRA=extra,$
                                    filename=filename[0])
   end
   print, 'ntor = ', ntor
   if(n_elements(slice) eq 0) then last=1
   if(n_elements(scale) eq 0) then scale = 1.
   if(n_elements(scale) eq 1 and n_elements(filename) gt 1) then $
     scale = replicate(scale, n_elements(filename))

;   psi0 = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra)
;   psi0_r = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra,op=2)
;   psi0_z = read_field('psi',x,z,t,/equilibrium,_EXTRA=extra,op=3)
;   i0   = read_field('i'  ,x,z,t,/equilibrium,_EXTRA=extra)

   threed = read_parameter('3d', _EXTRA=extra, filename=filename[0])
   print, '3D = ', threed

   if(n_elements(filename) gt 1 and n_elements(linfac) eq 0 and $
      n_elements(scale) gt 0) then linfac=scale

   if(threed eq 1) then begin
      psi0 = read_field('psi',x,z,t,slice=-1,_EXTRA=extra, $
                           filename=filename[0])
      psi0_r = read_field('psi',x,z,t,slice=-1,_EXTRA=extra,op=2, $
                             filename=filename[0])
      psi0_z = read_field('psi',x,z,t,slice=-1,_EXTRA=extra,op=3, $
                             filename=filename[0])
      i0   = read_field('i'  ,x,z,t,slice=-1,_EXTRA=extra, $
                           filename=filename[0])

      bx = read_field_3d('bx',phi,x,z,t,last=last,slice=slice, $
                         /linear,_EXTRA=extra, ntor=ntor, $
                         filename=filename,sum=sum)
      bz = read_field_3d('bz',phi,x,z,t,last=last,slice=slice, $
                         /linear,_EXTRA=extra, ntor=ntor, $
                         filename=filename,sum=sum)
   endif else begin
      psi0 = read_field('psi',x,z,t,slice=-1,_EXTRA=extra, $
                        filename=filename[0])
      psi0_r = read_field('psi',x,z,t,slice=-1,_EXTRA=extra,op=2, $
                          filename=filename[0])
      psi0_z = read_field('psi',x,z,t,slice=-1,_EXTRA=extra,op=3, $
                          filename=filename[0])
      i0   = read_field('i'  ,x,z,t,slice=-1,_EXTRA=extra, $
                        filename=filename[0])

      if(keyword_set(velocity)) then begin
         bx = read_field('vx',x,z,t,last=last,slice=slice, $
                         /linear,_EXTRA=extra,/complex, $
                         filename=filename,linfac=linfac,sum=sum)
         bz = read_field('vz',x,z,t,last=last,slice=slice, $
                         /linear,_EXTRA=extra,/complex, $
                         filename=filename,linfac=linfac,sum=sum)         
      endif else begin
         bx = read_field('bx',x,z,t,last=last,slice=slice, $
                         /linear,_EXTRA=extra,/complex, $
                         filename=filename,linfac=linfac,sum=sum)
         bz = read_field('bz',x,z,t,last=last,slice=slice, $
                         /linear,_EXTRA=extra,/complex, $
                         filename=filename,linfac=linfac,sum=sum)
      end
   endelse

   if(keyword_set(subtract_vacuum)) then begin
       if(extsubtract eq 0) then begin
           bx0 = read_field('bx',x,z,t,slice=0,/linear,$
                            linfac=linfac,sum=sum,_EXTRA=extra, $
                            filename=filename)
           bz0 = read_field('bz',x,z,t,slice=0,/linear,$
                            linfac=linfac,sum=sum,_EXTRA=extra, $
                            filename=filename)
           bx = bx - bx0
           bz = bz - bz0
       end
   endif


   y = z_matrix(x,z,t)

   ; calculate B.grad(psi)
   br = (bx*psi0_r + bz*psi0_z)

   if(keyword_set(vacfield)) then begin
      print, 'Using VACFIELD integrand: B.Grad(psi)/|Grad(psi)|'
      br = br / sqrt(psi0_r^2 + psi0_z^2)
   end
   
   ; convert to cgs
   get_normalizations, b0=b0_norm, n0=n0_norm, l0=l0_norm, _EXTRA=extra, $
     filename=filename[0]
   br = br*b0_norm

   schaffer_plot, br, x, z, t, psi0=psi0,i0=i0, q_val=q_val, $
                  label='!8B!Dn!N!6 (G)!X', points=points, bins=bins, ntor=ntor, $
                  overplot=overplot, filename=filename[0], _EXTRA=extra, $
                  ignore_jacobian=vacfield
end
