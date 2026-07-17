!-----------------------------------------------------------
! Brief: Tests whether a 3D point lies on the plane of a triangle, with a user-supplied tolerance.
!
! Parameters:
!   Input:  P(3),A(3),B(3),C(3) - test point and three triangle vertices
!   Input:  Tol - relative area tolerance
!   Output: Yes_on - .true. if P is coplanar with the triangle
!
! Notes:   Same area-decomposition strategy as the non-tolerance variant but
!          applies a relative area-ratio tolerance to be scale-independent.
!-----------------------------------------------------------

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
