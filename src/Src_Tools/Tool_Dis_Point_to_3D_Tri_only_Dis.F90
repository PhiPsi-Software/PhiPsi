!     ================================================= !
!             ____  _       _   ____  _____   _         !
!            |  _ \| |     |_| |  _ \|  ___| |_|        !
!            | |_) | |___   _  | |_) | |___   _         !
!            |  _ /|  _  | | | |  _ /|___  | | |        !
!            | |   | | | | | | | |    ___| | | |        !
!            |_|   |_| |_| |_| |_|   |_____| |_|        !
!     ================================================= !
!     PhiPsi:     a general-purpose computational       !
!                 mechanics program written in Fortran. !
!     Website:    http://phipsi.top                     !
!     Author:     Shi Fang, Huaiyin Institute of        !
!                 Technology, Huaian, JiangSu, China    !
!     Email:      shifang@hyit.edu.cn                   !
!     ------------------------------------------------- !
!     Please cite the following papers:                 !
!     (1)Shi F., Lin C. Modeling fluid-driven           !
!        propagation of 3D complex crossing fractures   !
!        with the extended finite element method.       !
!        Computers and Geotechnics, 2024, 172, 106482.  !
!     (2)Shi F., Wang D., Li H. An XFEM-based approach  !
!        for 3D hydraulic fracturing simulation         !
!        considering crack front segmentation. Journal  !
!        of Petroleum Science and Engineering, 2022,    !
!        214, 110518.                                   !
!     (3)Shi F., Wang D., Yang Q. An XFEM-based         !
!        numerical strategy to model three-dimensional  !
!        fracture propagation regarding crack front     !
!        segmentation. Theoretical and Applied Fracture !
!        Mechanics, 2022, 118, 103250.                  !
!     (4)Shi F., Liu J. A fully coupled hydromechanical !
!        XFEM model for the simulation of 3D non-planar !
!        fluid-driven fracture propagation. Computers   !
!        and Geotechnics, 2021, 132: 103971.            !
!     (5)Shi F., Wang X.L., Liu C., Liu H., Wu H.A. An  !
!        XFEM-based method with reduction technique     !
!        for modeling hydraulic fracture propagation    !
!        in formations containing frictional natural    !
!        fractures. Engineering Fracture Mechanics,     !
!        2017, 173: 64-90.                              !
!     ------------------------------------------------- !
 
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




#ifndef Silverfrost
if (isnan(Distance)) then
    print *, '    Error-2023061501 :: Distance is NAN in Tool_Dis_Point_to_3D_Tri_only_Dis.f90!'
    call Warning_Message('S',Keywords_Blank)
endif
#endif  
return 
end SUBROUTINE Tool_Dis_Point_to_3D_Tri_only_Dis                  
