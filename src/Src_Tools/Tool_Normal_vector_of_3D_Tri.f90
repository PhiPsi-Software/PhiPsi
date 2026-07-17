!-----------------------------------------------------------
! Brief: Calculate the unit outward normal of a 3D triangle.
!
! Parameters:
!   Input:  Tri_P1 - First vertex of the triangle (3-vector)
!   Input:  Tri_P2 - Second vertex of the triangle (3-vector)
!   Input:  Tri_P3 - Third vertex of the triangle (3-vector)
!   Output: Np     - Unit normal vector of the triangle (3-vector)
!
! Notes:   Uses cross product of edge vectors and normalizes the
!   result to a unit length vector.
!-----------------------------------------------------------

subroutine Tool_Normal_vector_of_3D_Tri( Tri_P1,Tri_P2,Tri_P3,Np)
!     Calculate the outward normal vector of a spatial triangle.
!     Modified on 2024-02-22.

!.......................
! Variable declaration.
!.......................
use Global_Float_Type     
implicit none
real(kind=FT),intent(in)::Tri_P1(3),Tri_P2(3),Tri_P3(3)
real(kind=FT),intent(out):: Np(3)
real(kind=FT) P1P2(3),P1P3(3)

Np(1:3) = ZR

P1P2 = Tri_P2 - Tri_P1
P1P3 = Tri_P3 - Tri_P1
call Vector_Cross_Product_3(P1P2,P1P3,Np) 

call Vector_Normalize(3,Np(1:3))   


return 
end SUBROUTINE Tool_Normal_vector_of_3D_Tri                   
