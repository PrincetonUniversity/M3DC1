pro plot_br, bins=bins, q_val=q_val, $
             subtract_vacuum=subtract_vacuum, ntor=ntor, $
             slice=slice, extsubtract=extsubtract, $
             overplot=overplot, filename=filename, scale=scale, $
             linfac=linfac, sum=sum, vacfield=vacfield, velocity=velocity, $
             jpar=jpar, tpoints=tpts,_EXTRA=extra

   if(n_elements(filename) eq 0) then filename='C1.h5'
   if(n_elements(ntor) eq 0) then begin
       ntor = read_parameter('ntor', _EXTRA=extra, filename=filename[0])
   endif
   if(n_elements(extsubtract) eq 0) then begin
       extsubtract = read_parameter('extsubtract' ,_EXTRA=extra,$
                                    filename=filename[0])
   end
   itor = read_parameter('itor', _EXTRA=extra, filename=filename[0])

   print, 'ntor = ', ntor
   if(n_elements(slice) eq 0) then last=1
   if(n_elements(scale) eq 0) then scale = 1.
   if(n_elements(scale) eq 1 and n_elements(filename) gt 1) then $
     scale = replicate(scale, n_elements(filename))

   threed = read_parameter('3d', _EXTRA=extra, filename=filename[0])
   version = read_parameter('version', _EXTRA=extra, filename=filename[0])
   print, '3D = ', threed

   if(n_elements(filename) gt 1 and n_elements(linfac) eq 0 and $
      n_elements(scale) gt 0) then linfac=scale

   if(threed eq 1) then begin
      psi = read_field('psi',x,z,t,slice=slice,tpoints=tpts, $
                        filename=filename[0],_EXTRA=extra)
      psi_r = read_field('psi',x,z,t,slice=slice,_EXTRA=extra,op=2, $
                             filename=filename[0],tpoints=tpts)
      psi_z = read_field('psi',x,z,t,slice=slice,_EXTRA=extra,op=3, $
                         filename=filename[0],tpoints=tpts)
      print, 'version = ', version
      if(version lt 35) then begin
         print, 'reading f'
         f_rp = read_field('f',x,z,t,slice=slice,_EXTRA=extra,op=12, $
                           filename=filename[0],tpoints=tpts)
         f_zp = read_field('f',x,z,t,slice=slice,_EXTRA=extra,op=13, $
                           filename=filename[0],tpoints=tpts)
      endif else begin
         print, 'reading fp'
         f_rp = read_field('fp',x,z,t,slice=slice,_EXTRA=extra,op=2, $
                           filename=filename[0],tpoints=tpts)
         f_zp = read_field('fp',x,z,t,slice=slice,_EXTRA=extra,op=3, $
                           filename=filename[0],tpoints=tpts)
      end
      
      i   = read_field('i'  ,x,z,t,slice=slice,_EXTRA=extra, $
                           filename=filename[0],tpoints=tpts)

      psi0 = fltarr(1,n_elements(x),n_elements(z))
      psi0_r = fltarr(1,n_elements(x),n_elements(z))
      psi0_z = fltarr(1,n_elements(x),n_elements(z))
      i0 = fltarr(1,n_elements(x),n_elements(z))

      psi0[0,*,*] = mean(psi, dimension=1)
      psi0_r[0,*,*] = mean(psi_r, dimension=1)
      psi0_z[0,*,*] = mean(psi_z, dimension=1)
      i0[0,*,*] = mean(i, dimension=1)

      r = psi

      for j=0, tpts-1 do begin
         psi[j,*,*] = psi[j,*,*] - psi0[0,*,*]
         psi_r[j,*,*] = psi_r[j,*,*] - psi0_r[0,*,*]
         psi_z[j,*,*] = psi_z[j,*,*] - psi0_z[0,*,*]
         r[j,*,*] = radius_matrix(x,z,t)
      end
      if(itor eq 0) then r=1.

      bx = -psi_z/r - f_rp
      bz =  psi_r/r - f_zp
      
   endif else begin
      psi0 = read_field('psi',x,z,t,slice=-1,_EXTRA=extra, $
                        filename=filename[0])
      i0   = read_field('i'  ,x,z,t,slice=-1,_EXTRA=extra, $
                        filename=filename[0])
      if(not keyword_set(jpar)) then begin
         psi0_r = read_field('psi',x,z,t,slice=-1,_EXTRA=extra,op=2, $
                             filename=filename[0])
         psi0_z = read_field('psi',x,z,t,slice=-1,_EXTRA=extra,op=3, $
                             filename=filename[0])
      end

      if(keyword_set(velocity)) then begin
         bx = read_field('vx',x,z,t,last=last,slice=slice, $
                         /linear,_EXTRA=extra,/complex,/mks, $
                         filename=filename,linfac=linfac,sum=sum)
         bz = read_field('vz',x,z,t,last=last,slice=slice, $
                         /linear,_EXTRA=extra,/complex,/mks, $
                         filename=filename,linfac=linfac,sum=sum)         
      endif else if(keyword_set(jpar)) then begin
         br = read_field('jpar',x,z,t,last=last,slice=slice, $
                         /linear,_EXTRA=extra,/complex,/mks, $
                         filename=filename,linfac=linfac,sum=sum)
      endif else begin
         bx = read_field('bx',x,z,t,last=last,slice=slice, $
                         /linear,_EXTRA=extra,/complex,/mks, $
                         filename=filename,linfac=linfac,sum=sum)
         bz = read_field('bz',x,z,t,last=last,slice=slice, $
                         /linear,_EXTRA=extra,/complex,/mks, $
                         filename=filename,linfac=linfac,sum=sum)
      end
   endelse

   if(keyword_set(subtract_vacuum)) then begin
       if(extsubtract eq 0) then begin
           bx0 = read_field('bx',x,z,t,slice=0,/linear,$
                            linfac=linfac,sum=sum,_EXTRA=extra, $
                            filename=filename,tpoints=tpts)
           bz0 = read_field('bz',x,z,t,slice=0,/linear,$
                            linfac=linfac,sum=sum,_EXTRA=extra, $
                            filename=filename,tpoints=tpts)
           bx = bx - bx0
           bz = bz - bz0
       end
   endif


   y = z_matrix(x,z,t)

   ; calculate B.grad(psi)
   if(not keyword_set(jpar)) then begin
      br = bx
      sz = size(br)
      for k=0, sz[1]-1 do begin
         br[k,*,*] = (bx[k,*,*]*psi0_r + bz[k,*,*]*psi0_z)
      end
   end

   if(keyword_set(vacfield)) then begin
      print, 'Using VACFIELD integrand: B.Grad(psi)/|Grad(psi)|'
      br = br / sqrt(psi0_r^2 + psi0_z^2)
   end
   
   if(keyword_set(velocity)) then begin
      symbol='!8V!Dn!N'
      units='!6m/s!X'
   endif else if(keyword_set(jpar)) then begin
      symbol='!8J!D!3||!N!6 (A/m!U2!N)'
      units='A/m!U2!N'
   endif else begin
      symbol='!8B!Dn!N!6'
      units='G'
      br = br*1e4
   end

   schaffer_plot, br, x, z, t, psi0=psi0,i0=i0, q_val=q_val, $
                  symbol=symbol, points=points, bins=bins, ntor=ntor, $
                  overplot=overplot, filename=filename[0], $
                  ignore_jacobian=vacfield, units=units, $
                  dpsi0_dx=psi0_r, dpsi0_dz=psi0_z, _EXTRA=extra
end
