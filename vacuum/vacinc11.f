      MODULE vacuum11_mod

      REAL, DIMENSION(:), ALLOCATABLE :: xobs, zobs

      CONTAINS
c................................................................
      SUBROUTINE vacuum11_alloc ( ndim )
c...............................................................

      ALLOCATE ( xobs(ndim), zobs(ndim) )

      xobs = 0.0
      zobs = 0.0

      RETURN
      END SUBROUTINE vacuum11_alloc
c................................................................
      SUBROUTINE vacuum11_dealloc
c...............................................................

      DEALLOCATE ( xobs, zobs )

      RETURN
      END SUBROUTINE vacuum11_dealloc
c...............................................................
      END MODULE vacuum11_mod
