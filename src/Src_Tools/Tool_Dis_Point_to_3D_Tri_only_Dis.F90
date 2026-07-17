!-----------------------------------------------------------
! Brief: Lightweight variant of Tool_Dis_Point_to_3D_Tri that returns only the signed distance.
!
! Parameters:
!   Input:  Point, Tri_P1, Tri_P2, Tri_P3 - query point and triangle vertices
!   Output: Distance - signed distance to the triangle plane
!
! Notes:   Skips the perpendicular foot and inside/edge detection; cheaper when
!          only the distance value is needed (e.g., inside large inner loops).
!-----------------------------------------------------------

subroutine Tool_Dis_Point_to_3D_Tri_only_Dis(Point,Tri_P1,Tri_P2,Tri_P3, Distance)
! Calculating the distance from a point to a triangle in space, Reference: Jones_1995_3D Distance
! from a Point to a Triangle
! PER indicates the coordinates of the foot of the perpendicular.
! Regarding the determination of the distance sign: it is positive if it is consistent with the
! orientation from the origin O to the plane (c_Vector_o_Orient(1:3)) (this method of judgment was
! later abandoned)
!      
! Tool_Dis_Point_to_3D_Tri.f modified. 2023-01-14.
! Only compute the signed distance compared to Tool_Dis_Point_to_3D_Tri.f.
!      2023-01-14.

!......................
! Variable Declaration
!......................
use Global_Float_Type     
use Global_Common
implicit none
real(kind=FT),intent(in)::Point(3),Tri_P1(3),Tri_P2(3),Tri_P3(3)
real(kind=FT),intent(out):: Distance
real(kind=FT) P1P2(3),P1P3(3),P1P0(3),Np(3)
real(kind=FT) abs_P1P0,abs_Np
real(kind=FT) B_i(3),C_i(3),X_i(3)
real(kind=FT) :: s

Distance   =  ZR

P1P2 = Tri_P2 - Tri_P1
P1P3 = Tri_P3 - Tri_P1
Np(1) = P1P2(2) * P1P3(3) - P1P2(3) * P1P3(2)
Np(2) = P1P2(3) * P1P3(1) - P1P2(1) * P1P3(3)
Np(3) = P1P2(1) * P1P3(2) - P1P2(2) * P1P3(1)

P1P0     = Point - Tri_P1
abs_P1P0 = sqrt(P1P0(1)**2 + P1P0(2)**2 + P1P0(3)**2)
abs_Np   = sqrt(Np(1)**2   + Np(2)**2   + Np(3)**2)

if (abs_Np <= Tol_15) then
    Distance = ZR
    return
endif



B_i = Tri_P2-Tri_P1
C_i = Tri_P3-Tri_P1
X_i = Point-Tri_P1


s = dot_product(P1P0, Np)
Distance = s / abs_Np




if (isnan(Distance)) then
    print *, '    Error-2023061501 :: Distance is NAN in Tool_Dis_Point_to_3D_Tri_only_Dis.f90!'
    call Warning_Message('S',Keywords_Blank)
endif
return 
end SUBROUTINE Tool_Dis_Point_to_3D_Tri_only_Dis                  
