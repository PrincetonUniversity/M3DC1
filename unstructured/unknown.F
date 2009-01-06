      subroutine local(xi,zi,b,theta,itri)
      use basic
      implicit none
      integer itri, nodeids(4)
      real xi, zi, theta, b, small, co, sn
      dimension xi(*),zi(*)
      double precision coords(3)
      small = 1.e-8
!
!     calculates the coefficients xi and zi for region ll in the local
!     coordinate expansion:  x = Sum{ xi(i)*si**mi(i)*eta**ni(i)}; {i=1,3}
!                            z = Sum{ zi(i)*si**mi(i)*eta**ni(i)}; {i=1,3}

      call nodfac(itri,nodeids)
      call xyznod(nodeids(2),coords)
      
      co = cos(theta)
      sn = sin(theta)
      xi(1) = coords(1)
      xi(2) =  co
      xi(3) = -sn
      call xyznod(nodeids(1),coords)
      zi(1) = coords(2)
      zi(2) =  sn
      zi(3) =  co
      return
      end
!============================================================
      subroutine nodcoord(inode,x,z)
      use basic
      implicit none
      integer, intent(in) :: inode
      real, intent(out) :: x, z
      double precision :: coords(3)
      real :: xmin, zmin

      call xyznod(inode,coords)
      
      select case(nonrect)
      case(0)
         call getmincoord(xmin,zmin)
         x = coords(1) - xmin + xzero
         z = coords(2) - zmin + zzero
      case(1)
         x = coords(1)
         z = coords(2)
      end select
      end
!============================================================
      subroutine whattri(x,z,itri,xref,zref)
      use basic
      implicit none
      integer itri, nodeids(4)
      real x, z, xref, zref
      double precision coords(3), x2, z2


      select case(nonrect)
      case(0)
         x2 = x - xzero
         z2 = z - zzero
      case(1)
         x2 = x
         z2 = z
      end select

      call usesearchstructure(x2,z2,itri)

      if(itri.lt.0) return

      call nodfac(itri,nodeids)
      call xyznod(nodeids(1), coords)
      
      select case(nonrect)
      case(0)
         xref = coords(1) + xzero
         zref = coords(2) + zzero
      case(1)
         xref = coords(1)
         zref = coords(2)
      end select

      return
      end
