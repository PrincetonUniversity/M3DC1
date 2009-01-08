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

      real, intent(in) :: x, z
      integer, intent(out) :: itri
      real, intent(out) :: xref, zref
      double precision :: x2, z2
      integer :: nodeids(4)
      real :: xmin, zmin


      call getmincoord(xmin,zmin)
      select case(nonrect)
      case(0)
         x2 = x - xzero + xmin
         z2 = z - zzero - zmin
      case(1)
         x2 = x
         z2 = z
      end select

      call usesearchstructure(x2,z2,itri)

      if(itri.lt.0) return

      call nodfac(itri,nodeids)
      call nodcoord(nodeids(1), xref, zref)

      return
      end
