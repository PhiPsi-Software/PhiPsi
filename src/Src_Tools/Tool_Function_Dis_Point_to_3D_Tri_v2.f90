!-----------------------------------------------------------
! Brief: Unsigned distance from a point to a 3D triangle
!
! Parameters:
!   Input:  Point(3)    - query point coordinates
!   Input:  Tri_P1(3)   - first triangle vertex
!   Input:  Tri_P2(3)   - second triangle vertex
!   Input:  Tri_P3(3)   - third triangle vertex
!   Returns: real - unsigned perpendicular distance
!
! Notes:   Uses the simplified Jones 1995 formula; not signed.
!-----------------------------------------------------------

function Tool_Function_Dis_Point_to_3D_Tri_v2(Point, Tri_P1,Tri_P2,Tri_P3)
!     Calculating the distance from a point to a triangle in space, Reference: Jones_1995_3D Distance from a Point to a Triangle
!     Tool_Function_Dis_Point_to_3D_Tri_v2 is a simplified version of Tool_Dis_Point_to_3D_Tri, only calculating the unsigned distance.
!     2022-05-05. NEWFTU2022050502.

!......................
! Variable Declaration
!......................
use Global_Float_Type     
implicit none
real(kind=FT),intent(in)::Point(3), Tri_P1(3),Tri_P2(3),Tri_P3(3)
real(kind=FT) P1P2(3),P1P3(3),P1P0(3),Np(3),cosa
real(kind=FT) abs_P1P0,abs_Np
real(kind=FT) :: Tool_Function_Dis_Point_to_3D_Tri_v2

Tool_Function_Dis_Point_to_3D_Tri_v2   =  ZR

P1P2 = Tri_P2 - Tri_P1
P1P3 = Tri_P3 - Tri_P1


call Vector_Cross_Product_3(P1P2,P1P3,Np)  

P1P0 = Point - Tri_P1
abs_P1P0 = sqrt(P1P0(1)**2 + P1P0(2)**2 + P1P0(3)**2)
abs_Np   = sqrt(Np(1)**2   + Np(2)**2   + Np(3)**2)
cosa = dot_product(P1P0,Np)/abs_P1P0/abs_Np 


Tool_Function_Dis_Point_to_3D_Tri_v2 = abs(abs_P1P0 * cosa)


end function  Tool_Function_Dis_Point_to_3D_Tri_v2                 
