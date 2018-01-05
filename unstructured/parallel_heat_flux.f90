module parallel_heat_flux
  implicit none

contains

  function kappar_p(e,f)
    use basic
    use m3dc1_nint

    implicit none

    vectype, dimension(dofs_per_element) :: kappar_p
    vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f
    vectype, dimension(dofs_per_element) :: temp

    if(gam.le.1) then 
       kappar_p = 0.
       return
    end if

    temp79a = -kar79(:,OP_1)*b2i79(:,OP_1)

    ! for ikappar_ni==0, don't include 1/n term
    if(ikappar_ni.eq.1) temp79a = temp79a*ni79(:,OP_1)

    ! B . Grad(p)
    ! [ p, psi] / R + 1/R^2 F p' - <f', p>
    temp79b = ri_79*(f(:,OP_DZ)*pstx79(:,OP_DR) - f(:,OP_DR)*pstx79(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
    temp79b = temp79b &
         + ri2_79*bztx79(:,OP_1)*f(:,OP_DP) &
         - (bftx79(:,OP_DZP)*f(:,OP_DZ) + bftx79(:,OP_DRP)*f(:,OP_DR))
#endif

    if(surface_int.eq.1) then
       temp = 0.
    else
       temp = intx5(e(:,:,OP_DZ),ri_79,temp79a,temp79b,pstx79(:,OP_DR)) &
            - intx5(e(:,:,OP_DR),ri_79,temp79a,temp79b,pstx79(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
       temp = temp &
            - intx5(e(:,:,OP_DP),ri2_79,temp79a,temp79b,bztx79(:,OP_1)) &
            - intx4(e(:,:,OP_DZ),temp79a,temp79b,bftx79(:,OP_DZP)) &
            - intx4(e(:,:,OP_DR),temp79a,temp79b,bftx79(:,OP_DRP))
#endif
    end if
    
    kappar_p = (gam - 1.) * temp
  end function kappar_p


  function kappar_pn(e,f,g)
    use basic
    use m3dc1_nint

    implicit none

    vectype, dimension(dofs_per_element) :: kappar_pn
    vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f, g
    vectype, dimension(dofs_per_element) :: temp

    if(gam.le.1) then 
       kappar_pn = 0.
       return
    end if
    
    temp79a = kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)**2*f(:,OP_1)

    ! B . Grad(n)
    ! [ n, psi] / R + 1/R^2 F n' - <f', n>
    temp79b = ri_79*(g(:,OP_DZ)*pstx79(:,OP_DR) - g(:,OP_DR)*pstx79(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
    temp79b = temp79b &
         + ri2_79*bztx79(:,OP_1)*g(:,OP_DP) &
         - (bftx79(:,OP_DZP)*g(:,OP_DZ) + bftx79(:,OP_DRP)*g(:,OP_DR))
#endif

    if(surface_int.eq.1) then
       temp = 0.
    else
       temp = intx5(e(:,:,OP_DZ),ri_79,temp79a,temp79b,pstx79(:,OP_DR)) &
            - intx5(e(:,:,OP_DR),ri_79,temp79a,temp79b,pstx79(:,OP_DZ))
#if defined(USE3D) || defined(USECOMPLEX)
       temp = temp &
            + intx5(e(:,:,OP_DP),ri2_79,temp79a,temp79b,bztx79(:,OP_1)) &
            - intx4(e(:,:,OP_DZ),temp79a,temp79b,bftx79(:,OP_DZP)) &
            - intx4(e(:,:,OP_DR),temp79a,temp79b,bftx79(:,OP_DRP))
#endif
    end if
    
    kappar_pn = (gam - 1.) * temp
  end function kappar_pn



  function kappar_psipsip(e,f,g,h)
    use basic
    use m3dc1_nint

    implicit none

    vectype, dimension(dofs_per_element) :: kappar_psipsip
    vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h
    vectype, dimension(dofs_per_element) :: temp

    if(gam.le.1) then 
       kappar_psipsip = 0.
       return
    end if
    
    temp79a = -ri2_79*kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)

    ! [p, psi]
    temp79b = h(:,OP_DZ)*g(:,OP_DR) - h(:,OP_DR)*g(:,OP_DZ)

    if(surface_int.eq.1) then
       temp = 0.
    else
       temp = intx4(e(:,:,OP_DZ),temp79a,temp79b,f(:,OP_DR)) &
            - intx4(e(:,:,OP_DR),temp79a,temp79b,f(:,OP_DZ))
    end if
    
    kappar_psipsip = (gam - 1.) * temp
  end function kappar_psipsip


  function kappar_psipsipn(e,f,g,h,i)
    use basic
    use m3dc1_nint

    implicit none

    vectype, dimension(dofs_per_element) :: kappar_psipsipn
    vectype, intent(in), dimension(dofs_per_element,MAX_PTS,OP_NUM) :: e
    vectype, intent(in), dimension(MAX_PTS,OP_NUM) :: f,g,h,i
    vectype, dimension(dofs_per_element) :: temp

    if(gam.le.1) then 
       kappar_psipsipn = 0.
       return
    end if
    
    temp79a = ri2_79*kar79(:,OP_1)*b2i79(:,OP_1)*ni79(:,OP_1)**2*h(:,OP_1)

    ! [n, psi]
    temp79b = i(:,OP_DZ)*g(:,OP_DR) - i(:,OP_DR)*g(:,OP_DZ)

    if(surface_int.eq.1) then
       temp = 0.
    else
       temp = intx4(e(:,:,OP_DZ),temp79a,temp79b,f(:,OP_DR)) &
            - intx4(e(:,:,OP_DR),temp79a,temp79b,f(:,OP_DZ))
    end if
    
    kappar_psipsipn = (gam - 1.) * temp
  end function kappar_psipsipn

end module parallel_heat_flux
