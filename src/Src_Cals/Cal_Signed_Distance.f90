!-----------------------------------------------------------
! Brief: Signed perpendicular distance from a point to a 2D line.
!
! Parameters:
!   Input:  Line_AB   - endpoints of the line ((2,2): x,y per row)
!           Point_C   - query point coordinates (x,y)
!   Output: S_Distance - signed distance (det / |AB|)
!
! Notes:   Sign follows the 2x2 determinant orientation of
!          (B-A, C-A), so positive/negative indicates which side
!          of the directed line the point lies on.
!-----------------------------------------------------------

subroutine Cal_Signed_Distance(Line_AB,Point_C,S_Distance)
!     Calculate Symbol Distance
!     This function calculates the signed distance from the Point_C to the Line_AB.
use Global_Float_Type      
implicit none
real(kind=FT),intent(in)::Line_AB(2,2),Point_C(2)
real(kind=FT),intent(out)::S_Distance

real(kind=FT) tem_1(2,2), tem_2(2)
real(kind=FT) tem_Det,tem_Norm


tem_1(1,:) = Line_AB(2,:)-  Line_AB(1,:)
tem_1(2,:) = Point_C     -  Line_AB(1,:)


tem_2(:)   = Line_AB(2,:)-Line_AB(1,:)


tem_Det = tem_1(1,1)*tem_1(2,2) - tem_1(2,1)*tem_1(1,2)

call Vector_Norm2(2,tem_2,tem_Norm) 

S_Distance = tem_Det / tem_Norm

return 
end subroutine Cal_Signed_Distance                          
