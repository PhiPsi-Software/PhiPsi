!-----------------------------------------------------------
! Brief: Tests whether a 3D point lies inside or on a quadrilateral pyramid.
!
! Parameters:
!   Input:  Point(3) - XYZ coordinates of the test point
!   Input:  A,B,C,D,E(3) - apex and base vertex coordinates of the pyramid
!   Output: Yes_in - .true. if the point is strictly inside the pyramid
!   Output: Yes_on - .true. if the point lies on the pyramid boundary
!
! Notes:   Splits the pyramid into two tetrahedra (A-BCE and A-CDE) and
!          delegates each test to Tool_Yes_Point_in_3D_Tetrahedron. Exits
!          early on first positive match for performance.
!-----------------------------------------------------------

subroutine Tool_Yes_Point_in_3D_Pyramid(Point, A,B,C,D,E, Yes_in,Yes_on)
!     Determine whether a point in space is inside a spatial quadrilateral pyramid (split it into two spatial tetrahedrons to determine, A-BCE, A-CDE)

!......................
! Variable Declaration
!......................
use Global_Float_Type
implicit none
real(kind=FT),intent(in)::Point(3), A(3),B(3),C(3),D(3),E(3)
logical,intent(out):: Yes_in,Yes_on
logical c_Yes_in_1,c_Yes_on_1,c_Yes_in_2,c_Yes_on_2
real(kind=FT) coor_x_max,coor_x_min,coor_y_max,coor_y_min
real(kind=FT) coor_z_max,coor_z_min

Yes_in = .False.
Yes_on = .False.

coor_x_max = max(A(1),B(1),C(1),D(1),E(1))
if(Point(1) > (coor_x_max+Tol_8))then
    return
endif
coor_x_min = min(A(1),B(1),C(1),D(1),E(1))
if(Point(1) < (coor_x_min-Tol_8))then
    return
endif

coor_y_max = max(A(2),B(2),C(2),D(2),E(2))
if(Point(2) > (coor_y_max+Tol_8))then
    return
endif
coor_y_min = min(A(2),B(2),C(2),D(2),E(2))
if(Point(2) < (coor_y_min-Tol_8))then
    return
endif

coor_z_max = max(A(3),B(3),C(3),D(3),E(3))
if(Point(3) > (coor_z_max+Tol_8))then
    return
endif
coor_z_min = min(A(3),B(3),C(3),D(3),E(3))
if(Point(3) < (coor_z_min-Tol_8))then
    return
endif


call Tool_Yes_Point_in_3D_Tetrahedron(Point, A,B,C,E, c_Yes_in_1,c_Yes_on_1)
if (c_Yes_in_1) then
    Yes_in = .True.
    return
endif
if (c_Yes_on_1) then
    Yes_on = .True.
    return
endif

call Tool_Yes_Point_in_3D_Tetrahedron(Point, A,C,D,E, c_Yes_in_2,c_Yes_on_2)
if (c_Yes_in_2) then
    Yes_in = .True.
    return
endif
if (c_Yes_on_2) then
    Yes_on = .True.
    return
endif


return 
end subroutine Tool_Yes_Point_in_3D_Pyramid                
