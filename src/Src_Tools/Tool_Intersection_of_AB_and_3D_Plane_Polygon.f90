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
 
subroutine Tool_Intersection_of_AB_and_3D_Plane_Polygon(A,B, &
                  num_Edges, Plane_Polygon, &
                  Yes_Inter,InterSection_P) 
! Used to calculate whether a spatial line and a spatial polygon have an intersection point.
! Main idea: Split it into triangles for evaluation.
! 2023-01-09.

use Global_Float_Type
use Global_Elem_Area_Vol

implicit none
real(kind=FT),intent(in)::A(3),B(3)
integer,intent(in)::num_Edges
real(kind=FT),intent(in):: Plane_Polygon(num_Edges,3)
real(kind=FT),intent(out)::InterSection_P(3)
logical,intent(out)::Yes_Inter


real(kind=FT) Tool_Function_2Point_Dis_3D
logical c_Yes_Inter
real(kind=FT) c_InterSection_P(3),DIS_AB
integer i_Tri
real(kind=FT) c_Tri(3,3)
real(kind=FT) AB_min_x,AB_max_x,AB_min_y,AB_max_y,AB_min_z,AB_max_z
real(kind=FT) Polygon_min_x,Polygon_max_x,Polygon_min_y,Polygon_max_y,Polygon_min_z,Polygon_max_z


Yes_Inter = .False.
InterSection_P(1:3)=ZR 


DIS_AB = Tool_Function_2Point_Dis_3D(A,B)
if (DIS_AB <=Tol_11)then
  return
endif

Polygon_min_x = minval(Plane_Polygon(:,1))
AB_max_x =  max(A(1),B(1))
if(AB_max_x < Polygon_min_x) return

Polygon_min_y = minval(Plane_Polygon(:,2))
AB_max_y =  max(A(2),B(2))
if(AB_max_y < Polygon_min_y) return

Polygon_min_z = minval(Plane_Polygon(:,3))
AB_max_z =  max(A(3),B(3))
if(AB_max_z < Polygon_min_z) return

Polygon_max_x = maxval(Plane_Polygon(:,1))
AB_min_x =  min(A(1),B(1))
if(AB_min_x > Polygon_max_x) return

Polygon_max_y = maxval(Plane_Polygon(:,2))
AB_min_y =  min(A(2),B(2))
if(AB_min_y > Polygon_max_y) return

Polygon_max_z = maxval(Plane_Polygon(:,3))
AB_min_z =  min(A(3),B(3))
if(AB_min_z > Polygon_max_z) return

do i_Tri=1, num_Edges-2
  c_Tri(1,1:3) = Plane_Polygon(1,1:3)
  c_Tri(2,1:3) = Plane_Polygon(i_Tri+1,1:3)
  c_Tri(3,1:3) = Plane_Polygon(i_Tri+2,1:3)
  call Tool_Intersection_of_AB_and_Triangle_3D(A,B,c_Tri(1,1:3),c_Tri(2,1:3),c_Tri(3,1:3),c_Yes_Inter,c_InterSection_P) 
  if(c_Yes_Inter)then
      Yes_Inter = .True.
      InterSection_P(1:3)=c_InterSection_P
      return
  endif
enddo

return 
end SUBROUTINE Tool_Intersection_of_AB_and_3D_Plane_Polygon                       
