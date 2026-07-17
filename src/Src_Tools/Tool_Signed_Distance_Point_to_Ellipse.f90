!-----------------------------------------------------------
! Brief: Compute the signed distance from a point to a tilted ellipse.
!
! Parameters:
!   Input:  x0,y0     - Ellipse center coordinates
!   Input:  a,b       - Semi-axes of the ellipse
!   Input:  theta     - Orientation angle of the ellipse
!   Input:  Point_C   - Query point (2-vector)
!   Output: S_Distance- Signed distance: -1 inside, +1 outside, 0 on
!
! Notes:   Uses the in/out classifier of the oblique ellipse
!   tool and reports a unit-signed scalar.
!-----------------------------------------------------------

subroutine Tool_Signed_Distance_Point_to_Ellipse(x0,y0,a,b,theta,Point_C,S_Distance)
!     Calculate the signed distance from a point to a tilted ellipse
!     Negative inside the ellipse, positive outside the ellipse, zero on the ellipse.
!     2020-08-09.

use Global_Float_Type      
implicit none
real(kind=FT),intent(in)::x0,y0,a,b,theta,Point_C(2)
real(kind=FT),intent(out)::S_Distance

real(kind=FT) c_Dis,Tol
integer Return_Statu

Tol = 1.0D-10

call Tool_Yes_Point_in_Oblique_Ellipse(Point_C,x0,y0,a,b,theta,Return_Statu,Tol)
if(Return_Statu==1)then
    S_Distance =-ONE
elseif(Return_Statu==2)then
    S_Distance =ONE
else
    S_Distance =ZR
endif
return 
end SUBROUTINE Tool_Signed_Distance_Point_to_Ellipse                        
