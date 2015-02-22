pro continuity, _EXTRA=extra, nterms=nterms, terms=term, names=names, $
               title=title, x=x, z=z
  title = 'Continuity Equation'
  nterms = 4
  names = ['dn/dt','v.Grad(n)','n Div(v)','-D del^2(n)']
  
  complex = read_parameter('icomplex',_EXTRA=extra)
  itor = read_parameter('itor',_EXTRA=extra)
  ntor = read_parameter('ntor',_EXTRA=extra)
  rzero = read_parameter('rzero',_EXTRA=extra)
  numvar = read_parameter('numvar',_EXTRA=extra)
  ivform = read_parameter('ivform',_EXTRA=extra)
  denm = read_parameter('denm',_EXTRA=extra)

  print, 'numvar = ', numvar
  print, 'itor = ', itor

  if(complex eq 1) then begin
     gamma = read_gamma(_EXTRA=extra)

     n0 = read_field('den',x,z,t,_EXTRA=extra,/equilibrium)
     n0_r = read_field('den',x,z,t,_EXTRA=extra,/equilibrium,op=2)
     n0_z = read_field('den',x,z,t,_EXTRA=extra,/equilibrium,op=3)
     w0 = read_field('omega',x,z,t,_EXTRA=extra,/equilibrium)

     n1  = read_field('den',x,z,t,_EXTRA=extra,/linear,/complex)
     n1_lp  = read_field('den',x,z,t,_EXTRA=extra,/linear,/complex,op=7)
     if(itor eq 1) then begin
        n1_r  = read_field('den',x,z,t,_EXTRA=extra,/linear,/complex,op=2)
     end

     u1_r  = read_field('phi',x,z,t,_EXTRA=extra,/linear,/complex,op=2)
     u1_z  = read_field('phi',x,z,t,_EXTRA=extra,/linear,/complex,op=3)
     if(numvar ge 2) then begin
        w1 = read_field('omega',x,z,t,_EXTRA=extra,/linear,/complex)
     endif else begin
        w1 = n1*0.
     end
       
     if(numvar eq 3) then begin
        chi1_r = read_field('chi',x,z,t,_EXTRA=extra,/linear,/complex,op=2)
        chi1_z = read_field('chi',x,z,t,_EXTRA=extra,/linear,/complex,op=3)
        chi1_lp = read_field('chi',x,z,t,_EXTRA=extra,/linear,/complex,op=7)
     endif else begin
        chi1_r = n1*0.
        chi1_z = n1*0.
        chi1_lp = n1*0.
     end
     
     if(itor eq 1) then begin
        r = radius_matrix(x,z,t)
        rfac = complex(0,ntor)
        n1_lp = n1_lp - n1_r/r
     endif else begin
        r = 1.
        rfac = complex(0,ntor)/rzero
     end

     d = size(n0,/dim)
     term = complexarr(nterms,d[1],d[2])
     
     ; dn / dt
     print, 'defining dn/dt'
     term[0,*,*] = n1 * gamma[0]

     ; v.grad(n)
     print, 'defining v.Grad(n)'
     if(ivform eq 1) then begin
        term[1,*,*] = r*(n0_z*u1_r - n0_r*u1_z) $
                      + w0*rfac*n1 $
                      + (n0_r*chi1_r + n0_z*chi1_z) / r^2
     endif else begin
        term[1,*,*] = (n0_z*u1_r - n0_r*u1_z)/r $
                      + w0*rfac*n1 $
                      + n0_r*chi1_r + n0_z*chi1_z
     endelse

     ; n Div(v)
     print, 'defining n Div(v)'
     if(ivform eq 1) then begin
        term[2,*,*] = n0*(rfac*w1 + chi1_lp/r^2 $
                      + (itor eq 1)*(-2.*u1_z - 2.*chi1_r/r^3))
     endif else begin
        term[2,*,*] = n0*(rfac*w1 + chi1_lp)
     end

     ; -del^2 (n)
     print, 'defining del^2(n)'
     term[3,*,*] = -denm*n1_lp


     print, 'gamma = ', reform(gamma[0])
  end  

end
