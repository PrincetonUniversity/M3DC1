       interface
       subroutine set3(p1, p2, p3, i1, i2, i3, ind1, ind2)
	    integer, pointer :: p1, p2, p3, i1, i2, i3
	 integer :: ind1, ind2
	 end subroutine set3
	 end interface

	 interface
	 subroutine set4(p1, p2, p3, p4, i1, i2, i3, i4, ind1, ind2, ind3)
	    integer, pointer :: p1, p2, p3, p4, i1, i2, i3, i4
	 integer :: ind1, ind2, ind3
	 end subroutine set4
	 end interface
