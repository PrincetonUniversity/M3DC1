      MODULE vacuum13_mod

      REAL, DIMENSION(:), ALLOCATABLE :: xloop, zloop

      CONTAINS
c................................................................
      SUBROUTINE vacuum13_alloc ( ndim )
c...............................................................

      ALLOCATE ( xloop(ndim), zloop(ndim) )

      xloop = 0.0
      zloop = 0.0

      RETURN
      END SUBROUTINE vacuum13_alloc
c................................................................
      SUBROUTINE vacuum13_dealloc
c...............................................................

      DEALLOCATE ( xloop, zloop )

      RETURN
      END SUBROUTINE vacuum13_dealloc
c...............................................................
      END MODULE vacuum13_mod
