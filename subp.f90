      subroutine f07adf(n1,n2,ti,n3,ipiv,info1)
      call DGETRF(n1,n2,ti,n3,ipiv,info1)
      return
      end
      subroutine f07ajf(n1,ti,n2,ipiv,wkspce,n3,info2)
      call DGETRI(n1,ti,n2,ipiv,wkspce,n3,info2)
      return
      end
      function s17aef(arg,idum)
      s17aef = dbesj0(arg)
      return
      end
      function s17aff(arg,idum)
      s17aff = dbesj1(arg)
      return
      end
