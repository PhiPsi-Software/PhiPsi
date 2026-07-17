!-----------------------------------------------------------
! Brief: Tests whether a 2D point is inside, on, or outside a tilted (rotated) ellipse.
!
! Parameters:
!   Input:  point(2) - coordinates of the test point
!   Input:  x0,y0 - center of the ellipse
!   Input:  a,b - semi-major and semi-minor axes
!   Input:  theta - rotation angle in degrees
!   Input:  Tol - boundary tolerance
!   Output: Return_Statu - 0 on boundary, 1 inside, 2 outside
!
! Notes:   Uses the implicit rotated-ellipse equation; theta is converted
!          to radians internally.
!-----------------------------------------------------------

subroutine Tool_Yes_Point_in_Oblique_Ellipse(point, x0,y0,a,b,theta, Return_Statu,Tol)
!     Determine whether a point is inside, outside, or on a tilted ellipse.
!     Return_Status=0 indicates on the ellipse; =1 indicates inside; =2 indicates outside.
!     See the equation of an inclined ellipse in \Theoretical_Documents\011 How to Draw the Equation of Any Inclined Ellipse (Baidu Wenku)
!     2020-08-9.

use Global_Float_Type
use Global_Inclusion
use Global_Common

implicit none
real(kind=FT),intent(in)::point(2),x0,y0,a,b,theta
real(kind=FT),intent(in):: Tol
integer,intent(out):: Return_Statu
real(kind=FT) x,y,temp1,temp2,c,c_theta

c_theta = theta*pi/180.0D0
c = sqrt(a**2-b**2)
x = point(1)
y = point(2)
temp1 =(a**2-c**2*cos(c_theta)**2)*(x-x0)**2 + (a**2-c**2*sin(c_theta)**2)*(y-y0)**2 - &
c**2*sin(TWO*c_theta)*(x-x0)*(y-y0)
temp2 = (a**2)*(b**2)

if((temp1-temp2)>=Tol)then
    Return_Statu =2
elseif((temp1-temp2)<=-Tol)then
    Return_Statu =1
else 
    Return_Statu =0
endif

return 
end subroutine Tool_Yes_Point_in_Oblique_Ellipse                          
