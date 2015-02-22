pro vorticity, _EXTRA=extra, nterms=nterms, terms=term, names=names, $
               title=title, x=x, z=z
  
  title = 'Vorticity Equation'
  nterms = 3
  names = ['dV/dt','-JxB','Grad(p)']
  
  complex = read_parameter('icomplex',_EXTRA=extra)
  itor = read_parameter('itor',_EXTRA=extra)
  ntor = read_parameter('ntor',_EXTRA=extra)
  rzero = read_parameter('rzero',_EXTRA=extra)
  numvar = read_parameter('numvar',_EXTRA=extra)

  if(complex eq 1) then begin
     gamma = read_gamma(_EXTRA=extra)

     psi0 = read_field('psi',x,z,t,_EXTRA=extra,/equilibrium)
     dspsi0 = read_field('jphi',x,z,t,_EXTRA=extra,/equilibrium)
     i0 = read_field('i',x,z,t,_EXTRA=extra,/equilibrium)
     den0 = read_field('den',x,z,t,_EXTRA=extra,/equilibrium)

     psi1 = read_field('psi',x,z,t,_EXTRA=extra,/linear,/complex)
;     dspsi1 = read_field('jphi',x,z,t,_EXTRA=extra,/linear,/complex)
     u1 = read_field('phi',x,z,t,_EXTRA=extra,/linear,/complex)
     psi1_lp = read_field('psi',x,z,t,_EXTRA=extra,/linear,/complex,op=7)
     if(itor eq 1) then begin
        psi1_r = read_field('psi',x,z,t,_EXTRA=extra,/linear,/complex,op=2)
     end
       
     if(numvar ge 2) then begin
        i1 = read_field('i',x,z,t,_EXTRA=extra,/linear,/complex)
        f1 = read_field('f',x,z,t,_EXTRA=extra,/linear,/complex)
     endif else begin
        i1 = psi1*0.
        f1 = psi1*0.
     end
     if(numvar eq 3) then begin
        x1 = read_field('chi',x,z,t,_EXTRA=extra,/linear,/complex)
        p1_z = read_field('p',x,z,t,_EXTRA=extra,/linear,/complex,op=3)
     endif else begin
        x1 = psi1*0.
        p1_z = psi1*0.
     end
     
     if(itor eq 1) then begin
        r = radius_matrix(x,z,t)
        rfac = complex(0,ntor)
        dspsi1 = psi1_lp - psi1_r/r
     endif else begin
        r = 1.
        rfac = complex(0,ntor)/rzero
        dspsi1 = psi1_lp
     end

     d = size(psi0,/dim)
     term = complexarr(nterms,d[1],d[2])
     
     ; dV / dt
     print, 'defining dv/dt'
     term[0,*,*] = (r^2*den0*laplacian(u1,x,z,tor=itor) $
                +s_bracket(r^2*den0,u1,x,z) $
                -a_bracket(den0,x1,x,z)/r) * gamma[0]

     ; -JxB
     print, 'defining JxB'
     term[1,*,*] = a_bracket(psi0,dspsi1,x,z)/r $
                 + a_bracket(psi1,dspsi0,x,z)/r $
                 + a_bracket(rfac^2*f1,i0,x,z)/r $
                 - i0*rfac*dspsi1/r^2 $
                 - s_bracket(i0,rfac*psi1,x,z)/r^2 $
                 + dspsi0*delstar(rfac*f1,x,z,tor=itor) $
                 + s_bracket(rfac*f1,r^2*dspsi0,x,z)/r^2

     ; Grad(p)
     print, 'defining grad(p)'
     term[2,*,*] = (itor eq 1)*2.*p1_z

     print, 'gamma = ', reform(gamma[0])
  end
end
