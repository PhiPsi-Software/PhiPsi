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
 
subroutine Tool_Yes_Point_on_3D_Triangle_with_Tol(P,A,B,C,Yes_on,Tol)
! Determine whether a point in space lies on the plane of a spatial triangle.
! Algorithm: If point p is inside the triangle, then the original triangle ABC is divided into three
! triangles: ABP, ACP, and BCP.
! The sum of the areas of the three triangles equals the area of triangle ABC.
!
! Tol is the area-to-volume tolerance.
!
!2023-08-13. NEWFTU2023081301.
!

!......................
! Variable Declaration
!......................
use Global_Float_Type
implicit none
real(kind=FT),intent(in)::P(3),A(3),B(3),C(3),Tol
logical,intent(out):: Yes_on
real(kind=FT) coor_x_max,coor_x_min,coor_y_max,coor_y_min,coor_z_max,coor_z_min
real(kind=FT) Area_Tri,Area_1,Area_2,Area_3
Yes_on = .False.

coor_x_max = max(A(1),B(1),C(1))
coor_x_min = min(A(1),B(1),C(1))
coor_y_max = max(A(2),B(2),C(2))
coor_y_min = min(A(2),B(2),C(2))
coor_z_max = max(A(3),B(3),C(3))
coor_z_min = min(A(3),B(3),C(3))


call Tool_Area_Tri_3D(A,B,C,Area_Tri)
call Tool_Area_Tri_3D(A,B,P,Area_1)
call Tool_Area_Tri_3D(A,C,P,Area_2)
call Tool_Area_Tri_3D(B,C,P,Area_3)


if(abs(Area_1+Area_2+Area_3-Area_Tri)/Area_Tri<=Tol)then
  Yes_on = .True.
endif

return 
end SUBROUTINE Tool_Yes_Point_on_3D_Triangle_with_Tol              
