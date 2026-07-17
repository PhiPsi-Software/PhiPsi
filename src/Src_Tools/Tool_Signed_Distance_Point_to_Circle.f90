!-----------------------------------------------------------
! Brief: Compute the signed distance from a point to a 2D circle.
!
! Parameters:
!   Input:  x0,y0     - Circle center coordinates
!   Input:  r         - Circle radius
!   Input:  Point_C   - Query point (2-vector)
!   Output: S_Distance- Signed distance; negative inside, positive
!                        outside, zero on the circle
!
! Notes:   A relative tolerance is applied to declare the point
!   to lie on the circle when the distance to the boundary is
!   numerically negligible.
!-----------------------------------------------------------

subroutine Tool_Signed_Distance_Point_to_Circle(x0,y0,r,Point_C,S_Distance)
!     Calculate Symbol Distance
!     This function calculates the signed distance from the Point_C to the circle.
!     Negative inside the circle, positive outside the circle, zero on the circle

use Global_Float_Type      
implicit none
real(kind=FT),intent(in)::x0,y0,r,Point_C(2)
real(kind=FT),intent(out)::S_Distance

real(kind=FT) c_Dis

c_Dis = sqrt((x0-Point_C(1))**2 + (y0-Point_C(2))**2)

if(c_Dis < r)then
    S_Distance = -abs(c_Dis-r)
elseif(c_Dis > R)then
    S_Distance = abs(c_Dis-r)
endif

if(abs(c_Dis-R)<1.0D-8*r)then
    S_Distance = ZR
endif

return 
end SUBROUTINE Tool_Signed_Distance_Point_to_Circle                        
