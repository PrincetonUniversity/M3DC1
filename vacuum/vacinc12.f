      MODULE vacuum12_mod

      REAL, DIMENSION(:), ALLOCATABLE :: xobs, zobs

      CONTAINS
c................................................................
      SUBROUTINE vacuum12_alloc (mx, mz )
c...............................................................

      ALLOCATE ( xobs(mx), zobs(mz) )

      xobs = 0.0
      zobs = 0.0

      RETURN
      END SUBROUTINE vacuum12_alloc
c................................................................
      SUBROUTINE vacuum12_dealloc
c...............................................................

      DEALLOCATE ( xobs, zobs )

      RETURN
      END SUBROUTINE vacuum12_dealloc
c...............................................................
      END MODULE vacuum12_mod
