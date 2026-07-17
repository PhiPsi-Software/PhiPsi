!-----------------------------------------------------------
! Brief: Test whether a 3D point is inside a 3D planar polygon (area-sum method)
!
! Parameters:
!   Input:  num_Point       - number of polygon vertices
!   Input:  Polygon        - polygon vertex coordinates
!   Input:  Point(3)       - test point coordinates
!   Output: Yes_Point_in   - true if point is inside the polygon
!
! Notes:   Compares the sum of sub-triangle areas to the full polygon area.
!-----------------------------------------------------------

subroutine Tool_Yes_Point_in_3D_Plane_Polygon(num_Point,Polygon,Point, Yes_Point_in)
! Check whether the point is inside the spatial plane polygon.
! IMPROV2022122003.
! Ref: https://blog.csdn.net/tjcwt2011/article/details/102737081 or Ref: \theory_documents\043
! Determining Whether a Point is Inside a Polygon_2022-12-20.pdf

! (1) Area and discriminant method: Determine whether the sum of the areas of the triangles formed
! by the target point and each edge of the polygon is equal to the area of the polygon. If equal,
! the point is inside the polygon.
! (2) Angle sum and discrimination method: Determine whether the sum of the angles between the
! target point and all the edges is 360 degrees; if it is 360 degrees, the point is inside the
! polygon.
! (3) Ray Casting Method: Draw a ray starting from the target point and count the number of
! intersections this ray has with all the edges of the polygon. If there is an odd number of
! intersections, the point is inside; if there is an even number of intersections, the point is
! outside.

! Here, method (1) is used.

use Global_Float_Type 
implicit none
integer,intent(in) :: num_Point
real(kind=FT),intent(in) :: Point(3),Polygon(num_Point,3)
logical,intent(out) :: Yes_Point_in
real(kind=FT) Area_Polygon,Sum_Area,c_Area
real(kind=FT) Point_1(3),Point_2(3),Point_3(3)
integer num_lines,i_Lines

Yes_Point_in = .False.

call Tool_Area_3D_Plane_Polygon(num_Point,Polygon(1:num_Point,1:3),Area_Polygon)

Sum_Area  = ZR
num_lines = num_Point
Point_1  = Point
do i_Lines = 1,num_lines
    Point_2 = Polygon(i_Lines,1:3)
    if(i_Lines<num_lines) then
        Point_3 = Polygon(i_Lines+1,1:3)
    else
        Point_3 = Polygon(1,1:3)
    endif
    call Tool_Area_Tri_3D(Point_1,Point_2,Point_3,c_Area)
    Sum_Area = Sum_Area + c_Area
enddo


if(abs(Sum_Area-Area_Polygon)<=Tol_9) then
    Yes_Point_in = .True.
    return
endif

return 
end SUBROUTINE Tool_Yes_Point_in_3D_Plane_Polygon                     
