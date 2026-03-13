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
 
subroutine Tool_Yes_Point_in_3D_Plane_Polygon(num_Point,Polygon,Point,&
                                              Yes_Point_in)
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
