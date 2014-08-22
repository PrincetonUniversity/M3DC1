function read_b_at_points, r0, z0, a0, slice=slice, filename=filename, $
                           plasma=plasma, _EXTRA=extra

   n = n_elements(r0)
   if(n_elements(z0) ne n or n_elements(a0) ne n) then begin
       print, 'r0, z0, and a0 must have same dimensions'
       return, 0
   end
   
   complex = read_parameter('icomplex',file=filename)
   extsub = read_parameter('extsubtract',file=filename)
   linear = read_parameter('eqsubtract',file=filename)
   ntor = read_parameter('ntor',file=filename)

   psi_r = read_field('psi',x,z,t,_EXTRA=extra,slice=slice,complex=complex,$
                      op=2,linear=linear,filename=filename)
   psi_z = read_field('psi',x,z,t,_EXTRA=extra,slice=slice,complex=complex,$
                      op=3,linear=linear,filename=filename)
   if(complex eq 1) then begin
       f_r = read_field('f',x,z,t,_EXTRA=extra,slice=slice,complex=complex,$
                        op=2,linear=linear,filename=filename)
       f_z = read_field('f',x,z,t,_EXTRA=extra,slice=slice,complex=complex,$
                        op=3,linear=linear,filename=filename)
       f_rp = complex(0,ntor)*f_r
       f_zp = complex(0,ntor)*f_z
   end

   if(keyword_set(plasma)) then begin
       psi_r = psi_r - $
         read_field('psi',x,z,t,_EXTRA=extra,slice=0,complex=complex,$
                    op=2,linear=linear,filename=filename)
       psi_z = psi_z - $
         read_field('psi',x,z,t,_EXTRA=extra,slice=0,complex=complex,$
                          op=3,linear=linear,filename=filename)
       if(complex eq 1) then begin
           f_r = read_field('f',x,z,t,_EXTRA=extra,slice=0,complex=complex,$
                            op=2,linear=linear,filename=filename)
           f_z = read_field('f',x,z,t,_EXTRA=extra,slice=0,complex=complex,$
                            op=3,linear=linear,filename=filename)
           f_rp = f_rp - complex(0,ntor)*f_r
           f_zp = f_zp - complex(0,ntor)*f_z
       end
   end

   r = radius_matrix(x,z,t)
       
   br = -psi_z / r - f_rp
   bz =  psi_r / r - f_zp

   ir = n_elements(x)*(r0 - x[0]) / (x[n_elements(x)-1] - x[0])
   iz = n_elements(z)*(z0 - z[0]) / (z[n_elements(z)-1] - z[0])
   br_int = interpolate(reform(br),ir,iz)
   bz_int = interpolate(reform(bz),ir,iz)
   return, br_int*cos(a0) + bz_int*sin(a0) 
end
