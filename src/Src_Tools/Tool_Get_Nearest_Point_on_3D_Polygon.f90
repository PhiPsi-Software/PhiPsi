!-----------------------------------------------------------
! Brief: Find the point on a 3D polygon closest to a given query point.
!
! Parameters:
!   Input:  num_Point_Polygon - number of polygon vertices
!   Input:  Polygon - vertex coordinates of the polygon
!   Input:  Input_Point(3) - query point
!   Output: Nearest_Point(3) - closest point on the polygon
!   Output: Point_Pre_with_Gap(3) - previous vertex at gap distance
!   Output: Point_Next_with_Gap(3) - next vertex at gap distance
!
! Notes:   Performs a coarse vertex search then refines on the
!   two adjacent polygon edges; the gap outputs are vertex
!   offsets following the polygon's rotation order.
!-----------------------------------------------------------

subroutine Tool_Get_Nearest_Point_on_3D_Polygon(num_Point_Polygon,Polygon,Input_Point,&
Nearest_Point,Point_Pre_with_Gap,Point_Next_with_Gap) 
! Get the point on a 3D polygon closest to a given point. The point is not necessarily limited to
! the polygon's vertices.
! 2023-08-12.

!
! Point_Pre_with_Gap(3) refers to the previous point at a distance of Gap from the nearest point
! (following the rotation order of Polygons).
! Point_Next_with_Gap(3) is the next point at a distance of Gap from the nearest point (following
! the rotation order of the polygons).
!

use Global_Float_Type
use Global_Common

implicit none
integer,intent(in)::num_Point_Polygon
real(kind=FT),intent(in) :: Polygon(num_Point_Polygon,3)
real(kind=FT),intent(in) :: Input_Point(3)       
real(kind=FT),intent(out):: Nearest_Point(3),Point_Pre_with_Gap(3),Point_Next_with_Gap(3)  
integer Nearest_Point_Num,Nearest_Point_Pre_Num,Nearest_Point_Next_Num
real(kind=FT) :: Line_AB(3,2), Line_CD(3,2)
real(kind=FT) nearestPoint_to_AB(3),nearestPoint_to_CD(3)
real(kind=FT) distance_to_AB,distance_to_CD
real(kind=FT) length_AB,length_CD,c_length
real(kind=FT) Tool_Function_2Point_Dis_3D
real(kind=FT) :: Line_AB_23(2,3), Line_CD_23(2,3)


call Tool_Get_Nearest_Point_3D(num_Point_Polygon,Polygon(1:num_Point_Polygon,1:3),Input_Point,Nearest_Point_Num) 

if(Nearest_Point_Num==1) then
    Nearest_Point_Pre_Num = num_Point_Polygon
else
    Nearest_Point_Pre_Num = Nearest_Point_Num - 1
endif

if(Nearest_Point_Num == num_Point_Polygon) then
    Nearest_Point_Next_Num = 1
else
    Nearest_Point_Next_Num = Nearest_Point_Num + 1
endif



Line_AB(1:3,1) = Polygon(Nearest_Point_Pre_Num,1:3)
Line_AB(1:3,2) = Polygon(Nearest_Point_Num,1:3)
Line_CD(1:3,1) = Polygon(Nearest_Point_Num,1:3)
Line_CD(1:3,2) = Polygon(Nearest_Point_Next_Num,1:3)


length_AB = sqrt( (Line_AB(1,1) - Line_AB(1,2)) * (Line_AB(1,1) - Line_AB(1,2)) &
+ (Line_AB(2,1) - Line_AB(2,2)) * (Line_AB(2,1) - Line_AB(2,2)) &
+ (Line_AB(3,1) - Line_AB(3,2)) * (Line_AB(3,1) - Line_AB(3,2)) )
length_CD = sqrt( (Line_CD(1,1) - Line_CD(1,2)) * (Line_CD(1,1) - Line_CD(1,2)) &
+ (Line_CD(2,1) - Line_CD(2,2)) * (Line_CD(2,1) - Line_CD(2,2)) &
+ (Line_CD(3,1) - Line_CD(3,2)) * (Line_CD(3,1) - Line_CD(3,2)) )



Line_AB_23(1,1:3) = Line_AB(1:3,1)
Line_AB_23(2,1:3) = Line_AB(1:3,2)
Line_CD_23(1,1:3) = Line_CD(1:3,1)
Line_CD_23(2,1:3) = Line_CD(1:3,2)
call Tool_Nearest_Point_from_Point_to_Segment_3D(Input_Point, Line_AB_23, nearestPoint_to_AB, distance_to_AB)
call Tool_Nearest_Point_from_Point_to_Segment_3D(Input_Point, Line_CD_23, nearestPoint_to_CD, distance_to_CD)



if (abs(distance_to_AB - distance_to_CD) <= Tol_11) then
    Nearest_Point(1:3) = Polygon(Nearest_Point_Num,1:3)
    call Tool_Offset_Point_A_to_Point_B_3D(Nearest_Point(1:3), Line_AB(1:3,1), length_AB/TWO, Point_Pre_with_Gap)
    call Tool_Offset_Point_A_to_Point_B_3D(Nearest_Point(1:3), Line_CD(1:3,2), length_CD/TWO, Point_Next_with_Gap)

elseif (distance_to_AB <= distance_to_CD) then
    Nearest_Point(1:3) = nearestPoint_to_AB(1:3)
    c_length = Tool_Function_2Point_Dis_3D(Nearest_Point(1:3), Line_AB(1:3,1))
    call Tool_Offset_Point_A_to_Point_B_3D(Nearest_Point(1:3), Line_AB(1:3,1), c_length/TWO, Point_Pre_with_Gap)
    c_length = Tool_Function_2Point_Dis_3D(Nearest_Point(1:3), Line_AB(1:3,2))
    call Tool_Offset_Point_A_to_Point_B_3D(Nearest_Point(1:3), Line_AB(1:3,2), c_length/TWO, Point_Next_with_Gap)

else
    Nearest_Point(1:3) = nearestPoint_to_CD(1:3)
    c_length = Tool_Function_2Point_Dis_3D(Nearest_Point(1:3), Line_CD(1:3,1))
    call Tool_Offset_Point_A_to_Point_B_3D(Nearest_Point(1:3), Line_CD(1:3,1), c_length/TWO, Point_Pre_with_Gap)
    c_length = Tool_Function_2Point_Dis_3D(Nearest_Point(1:3), Line_CD(1:3,2))
    call Tool_Offset_Point_A_to_Point_B_3D(Nearest_Point(1:3), Line_CD(1:3,2), c_length/TWO, Point_Next_with_Gap)
endif

return 
end SUBROUTINE Tool_Get_Nearest_Point_on_3D_Polygon                  
