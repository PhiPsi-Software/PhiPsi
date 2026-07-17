!-----------------------------------------------------------
! Brief: Compute intersection of a line segment with a tilted 2D ellipse
!
! Parameters:
!   Input:  x0, y0   - Ellipse center
!   Input:  a, b     - Semi-axes of the ellipse
!   Input:  theta    - Rotation angle of the ellipse
!   Input:  Seg_A, Seg_B - Endpoints of the line segment
!   Output: num_Inter - Number of intersections (0 or 1)
!   Output: Inter     - Intersection point coordinates
!
! Notes:   Assumes one endpoint is inside and the other outside the
!          ellipse; uses an iterative bisection method on the parameter.
!-----------------------------------------------------------

subroutine Tool_Intersection_Line_and_Oblique_Ellipse( x0,y0,a,b,theta,Seg_A,Seg_B, num_Inter,Inter)
!     Calculate the intersection points of line segment AB and the ellipse, assuming that one endpoint of the line segment is inside the ellipse and the other is outside.
!     2020-08-09.

use Global_Float_Type      
use Global_Common
implicit none
real(kind=FT),intent(in):: x0,y0,a,b,theta,Seg_A(2),Seg_B(2)
integer,intent(out)::num_Inter
real(kind=FT),intent(out)::Inter(2)
real(kind=FT) :: Tol
integer Return_Statu_A,Return_Statu_B,Return_Statu_P

real(kind=FT) in_P(2),out_P(2),Possi_P(2)

integer :: i_Try


Tol = 1.0D-10
call Tool_Yes_Point_in_Oblique_Ellipse(Seg_A, x0,y0,a,b,theta, &
Return_Statu_A,Tol)
call Tool_Yes_Point_in_Oblique_Ellipse(Seg_B, x0,y0,a,b,theta, &
Return_Statu_B,Tol)
if((Return_Statu_A+Return_Statu_B)/=3)then
    num_Inter = 0
    return
endif



if(Return_Statu_A ==1)then
    in_P(1:2)  =Seg_A(1:2)
    out_P(1:2) =Seg_B(1:2)
else
    in_P(1:2)  =Seg_B(1:2)
    out_P(1:2) =Seg_A(1:2)
endif
do i_Try = 1,5000
    Possi_P(1:2) = (in_P(1:2) + out_P(1:2))/TWO
    call Tool_Yes_Point_in_Oblique_Ellipse(Possi_P, &
    x0,y0,a,b,theta, Return_Statu_P,Tol)
    if(Return_Statu_P==1)then   
        in_P = Possi_P
    elseif(Return_Statu_P==2)then
        out_P = Possi_P
    elseif(Return_Statu_P==0)then
        Inter(1:2) =Possi_P
        num_Inter = 1
        exit
    endif
enddo



return 
end subroutine Tool_Intersection_Line_and_Oblique_Ellipse  
