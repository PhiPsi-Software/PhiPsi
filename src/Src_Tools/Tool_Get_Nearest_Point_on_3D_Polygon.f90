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
real(kind=FT) Line_AB(2,3),Line_CD(2,3)
real(kind=FT) nearestPoint_to_AB(3),nearestPoint_to_CD(3)
real(kind=FT) distance_to_AB,distance_to_CD
real(kind=FT) length_AB,length_CD,c_length
real(kind=FT) Tool_Function_2Point_Dis_3D


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


Line_AB(1,1:3) = Polygon(Nearest_Point_Pre_Num,1:3)
Line_AB(2,1:3) = Polygon(Nearest_Point_Num,1:3)

Line_CD(1,1:3) = Polygon(Nearest_Point_Num,1:3)
Line_CD(2,1:3) = Polygon(Nearest_Point_Next_Num,1:3)

length_AB = Tool_Function_2Point_Dis_3D(Line_AB(1,1:3),Line_AB(2,1:3))
length_CD = Tool_Function_2Point_Dis_3D(Line_CD(1,1:3),Line_CD(2,1:3))


call Tool_Nearest_Point_from_Point_to_Segment_3D(Input_Point,Line_AB,nearestPoint_to_AB,distance_to_AB)
call Tool_Nearest_Point_from_Point_to_Segment_3D(Input_Point,Line_AB,nearestPoint_to_CD,distance_to_CD)



if(abs(distance_to_AB - distance_to_CD)<=Tol_11)then
    Nearest_Point(1:3) = Polygon(Nearest_Point_Num,1:3)
    call Tool_Offset_Point_A_to_Point_B_3D(Nearest_Point(1:3),Line_AB(1,1:3),length_AB/TWO,Point_Pre_with_Gap)
    call Tool_Offset_Point_A_to_Point_B_3D(Nearest_Point(1:3),Line_CD(2,1:3),length_AB/TWO,Point_Next_with_Gap)
elseif(distance_to_AB<=distance_to_CD)then
    Nearest_Point(1:3) = nearestPoint_to_AB(1:3)
    c_length = Tool_Function_2Point_Dis_3D(Nearest_Point(1:3),Line_AB(1,1:3))
    call Tool_Offset_Point_A_to_Point_B_3D(Nearest_Point(1:3),Line_AB(1,1:3),c_length/TWO,Point_Pre_with_Gap)
    c_length = Tool_Function_2Point_Dis_3D(Nearest_Point(1:3),Line_AB(2,1:3))
    call Tool_Offset_Point_A_to_Point_B_3D(Nearest_Point(1:3),Line_AB(2,1:3),c_length/TWO,Point_Next_with_Gap)
else
    Nearest_Point(1:3) = nearestPoint_to_CD(1:3)
    c_length = Tool_Function_2Point_Dis_3D(Nearest_Point(1:3),Line_CD(1,1:3))
    call Tool_Offset_Point_A_to_Point_B_3D(Nearest_Point(1:3),Line_CD(1,1:3),c_length/TWO,Point_Pre_with_Gap)
    c_length = Tool_Function_2Point_Dis_3D(Nearest_Point(1:3),Line_CD(2,1:3))
    call Tool_Offset_Point_A_to_Point_B_3D(Nearest_Point(1:3),Line_CD(2,1:3),c_length/TWO,Point_Next_with_Gap)
endif

return 
end SUBROUTINE Tool_Get_Nearest_Point_on_3D_Polygon                  
