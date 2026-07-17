!-----------------------------------------------------------
! Brief: Offset a single point by delta along the normal to AB.
!
! Parameters:
!   Input:  C              - point to be offset (2D)
!           Line_AB        - reference line segment (2x2)
!           offset_delta   - signed offset distance
!   Output: Offsetted_C_Up   - point shifted to the left of AB
!           Offsetted_C_Down - point shifted to the right of AB
!
! Notes:   If C coincides with an endpoint of Line_AB, the
!          offset is performed along the line direction first
!          before the normal shift to avoid degenerate cases.
!-----------------------------------------------------------

subroutine Cal_Offseted_Single_Point(C, Line_AB, offset_delta, Offsetted_C_Up, Offsetted_C_Down)
!=========================================================================
! Offset the point C along the up and down normal of line_AB by the 
! increment of offset_delta.
!=========================================================================
!2025-12-26.
use Global_Float_Type
implicit none

real(kind=FT), dimension(2), intent(in)    :: C
real(kind=FT), dimension(2,2), intent(in)  :: Line_AB
real(kind=FT), intent(in)                  :: offset_delta

real(kind=FT), dimension(2), intent(out)   :: Offsetted_C_Up
real(kind=FT), dimension(2), intent(out)   :: Offsetted_C_Down

real(kind=FT) :: a_x, a_y, b_x, b_y
real(kind=FT) :: theta
real(kind=FT), dimension(2) :: C_local
real(kind=FT), parameter :: eps = 1.0d-10

a_x = Line_AB(1, 1)
a_y = Line_AB(1, 2)
b_x = Line_AB(2, 1)
b_y = Line_AB(2, 2)

theta = atan2(b_y - a_y, b_x - a_x)

C_local = C

if (abs(C(1) - b_x) < eps .and. abs(C(2) - b_y) < eps) then
    C_local(1) = C_local(1) - offset_delta * cos(theta)
    C_local(2) = C_local(2) - offset_delta * sin(theta)

else if (abs(C(1) - a_x) < eps .and. abs(C(2) - a_y) < eps) then
    C_local(1) = C_local(1) + offset_delta * cos(theta)
    C_local(2) = C_local(2) + offset_delta * sin(theta)
end if

Offsetted_C_Up(1) = C_local(1) - offset_delta * sin(theta)
Offsetted_C_Up(2) = C_local(2) + offset_delta * cos(theta)

Offsetted_C_Down(1) = C_local(1) + offset_delta * sin(theta)
Offsetted_C_Down(2) = C_local(2) - offset_delta * cos(theta)
    
end subroutine Cal_Offseted_Single_Point