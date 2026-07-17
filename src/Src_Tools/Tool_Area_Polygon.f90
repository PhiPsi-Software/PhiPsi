!-----------------------------------------------------------
! Brief: Computes the signed area of a 2D polygon from its vertex coordinates.
!
! Parameters:
!   Input:  x, y - arrays of vertex x/y coordinates
!   Input:  nb   - number of vertices (closing vertex may repeat)
!   Output: area - signed area of the polygon
!
! Notes:   Uses the shoelace formula; tolerates polygons whose first and last
!          vertices coincide and self-intersecting curves.
!-----------------------------------------------------------

subroutine Tool_Area_Polygon(x,y,nb,area)

! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-04  Time: 12:24:06
use Global_Float_Type
implicit none
integer,intent(in)        :: nb
real(kind=FT),intent(in)  :: x(nb),y(nb)
real(kind=FT),intent(out) :: area
integer :: i, n, nm1
real(kind=FT) :: a


n = nb


if((abs(x(1)-x(n)) <=Tol_15) .and. (abs(y(1)-y(n))<=Tol_15))then
    n = n - 1
endif

select case (n)
case (:2)
    area = ZR
case (3)
    area=HLF*((x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1)))
case default
    nm1 = n - 1
    a = x(1)*(y(2) - y(n)) + x(n)*(y(1) - y(nm1))
    do i = 2, nm1
        a = a + x(i)*(y(i+1) - y(i-1))
    end do
    area = HLF*a
end select

return
end subroutine Tool_Area_Polygon
