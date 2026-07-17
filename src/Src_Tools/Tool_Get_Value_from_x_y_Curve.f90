!-----------------------------------------------------------
! Brief: Linearly interpolate a y value from a discrete (x,y) curve at a given x.
!
! Parameters:
!   Input:  Curve_x,Curve_y - tabulated curve coordinates
!   Input:  num_Point - number of tabulated points
!   Input:  x - query abscissa
!   Output: y - interpolated ordinate at x
!
! Notes:   Out-of-range queries are clamped to the first or last
!   ordinate. A zero-width segment guards against division by zero.
!-----------------------------------------------------------

subroutine Tool_Get_Value_from_x_y_Curve(Curve_x,Curve_y, num_Point,x,y)
!     Given an x value, obtain the corresponding y value from the xy curve, linearly.

use Global_Float_Type
implicit none
integer,intent(in)::  num_Point
real(kind=FT),intent(in):: Curve_x(num_Point), Curve_y(num_Point),x
real(kind=FT),intent(out):: y

integer :: i
integer location_x
real(kind=FT) L_x,R_x,L_y,R_y
real(kind=FT) :: temp

do i =1,num_Point-1
    if(x>=Curve_x(i) .and. x<=Curve_x(i+1))then
        location_x = i
        exit
    end if
end do

L_x = Curve_x(location_x)
R_x = Curve_x(location_x+1)
L_y = Curve_y(location_x)
R_y = Curve_y(location_x+1)

temp = R_x - L_x
if(R_x==L_x)then
    temp = 1.0D30
endif
y = L_y + (x-L_x)*(R_y-L_y)/temp

if (x < minval(Curve_x))then
    y = Curve_y(1)
endif
if (x > maxval(Curve_x))then
    y = Curve_y(num_Point)
endif

return 
end subroutine Tool_Get_Value_from_x_y_Curve                          
