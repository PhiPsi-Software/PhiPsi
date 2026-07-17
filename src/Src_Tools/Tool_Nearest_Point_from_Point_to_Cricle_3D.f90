!-----------------------------------------------------------
! Brief: Find the closest point on a 3D circle to a given query point
!
! Parameters:
!   Input:  In_Point - Query point in 3D
!   Input:  Center   - Circle center
!   Input:  r        - Circle radius
!   Input:  n        - Unit normal of the circle's plane
!   Output: Out_Point - Closest point lying on the circle
!
! Notes:   Projects the query vector onto the plane and rescales to the
!          radius; works for any plane orientation.
!-----------------------------------------------------------

subroutine Tool_Nearest_Point_from_Point_to_Cricle_3D(In_Point,Center,r,n,Out_Point)
!     Calculate the point on the circle that is closest to the vertex In_Point.
!     Ref: \theory_documents\027 Distance to Circles in 3D-2021-10-30.pdf Eq.(2)
!     2021-10-30.

use Global_Float_Type

implicit none
real(kind=FT),intent(in)::In_Point(3),Center(3),r,n(3)
real(kind=FT),intent(out)::Out_Point(3)
real(kind=FT) delta(3),tem(3),nor_tem

delta = In_Point-Center
tem = delta - dot_product(n,delta)*n
call Vector_Norm2(3,tem,nor_tem)   
Out_Point = Center + r*tem/nor_tem


return 
end SUBROUTINE Tool_Nearest_Point_from_Point_to_Cricle_3D                         
