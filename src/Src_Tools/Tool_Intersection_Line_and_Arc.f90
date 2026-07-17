!-----------------------------------------------------------
! Brief: Compute the intersection points of a line segment and an arc
!
! Parameters:
!   Input:  A, B            - Endpoints of the line segment
!   Input:  c_Arc_Crack_Coor - 11-element arc descriptor (center, radius,
!                              start/end radians, endpoints)
!   Output: num_Inter       - Number of arc-line intersections
!   Output: Inter           - Intersection coordinates (2,2)
!   Output: Yes_Cross       - True if any arc-line intersection exists
!
! Notes:   First finds line-circle intersections, then keeps only those
!          that lie within the arc's angular range.
!-----------------------------------------------------------

subroutine Tool_Intersection_Line_and_Arc( A,B, c_Arc_Crack_Coor, &
num_Inter,Inter,Yes_Cross)
!     Calculate the intersection points of line segment AB and the arc segment (Added on 2017-07-16)
!     A: Coordinates of point A of line segment AB;
!     B: Coordinates of point B of line segment AB;
!     c_Arc_Crack_Coor(1:10),   !10 elements of the arc segment: o_x, o_y, r, Radian_Start, Radian_End, Radian, Point_Start_x, Point_Start_y, Point_End_x, Point_End_y
!     ---------------------------
!     The calculation approach is to first find the intersection points between the line segment and the circle on which the arc lies. If there are intersection points, then determine whether each intersection point lies within the range of the arc.

use Global_Float_Type
use Global_Elem_Area_Vol

implicit none
real(kind=FT),intent(in)::A(2),B(2),c_Arc_Crack_Coor(11)
integer,intent(out)::num_Inter
real(kind=FT),intent(out)::Inter(2,2)
logical,intent(out)::Yes_Cross

real(kind=FT) o_x,o_y
real(kind=FT) :: r
integer c_num_Inter,State
real(kind=FT) c_Inter(2,2)
integer :: i_Inter
logical :: c_Yes_ON

Yes_Cross = .False.

o_x               = c_Arc_Crack_Coor(1)
o_y               = c_Arc_Crack_Coor(2)
r                 = c_Arc_Crack_Coor(4)

c_num_Inter = 0
call Tool_Intersection_Line_and_Circle(o_x,o_y, r,A,B, c_num_Inter,State, c_Inter)

num_Inter = 0
do i_Inter = 1,c_num_Inter
    call Tool_Yes_On_Arc(c_Inter(i_Inter,1),c_Inter(i_Inter,2), c_Arc_Crack_Coor,c_Yes_ON)
    if(c_Yes_ON)then
        Yes_Cross = .True.
        num_Inter = num_Inter +1
        Inter(num_Inter,1:2) = c_Inter(i_Inter,1:2)
    endif
enddo

return 
end subroutine Tool_Intersection_Line_and_Arc                         
