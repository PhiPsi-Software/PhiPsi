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
 
subroutine Tool_Nearest_Point_from_Point_to_Segment_3D(point,line,nearestPoint,distance)
! Calculate the closest point from a point to a 3D line segment.
!ref:
!https://math.stackexchange.com/questions/1905533/find-perpendicular-distance-from-point-to-line-in-3d
!2023-08-12.

use Global_Float_Type

implicit none
real(kind=FT),intent(in) :: point(3)
real(kind=FT),intent(in) :: line(2,3)
real(kind=FT),intent(out) :: nearestPoint(3),distance
real(kind=FT) :: P1(3),P2(3),t

real(kind=FT) :: P2_P1(3)
real(kind=FT) :: dis_to_P1,dis_to_P2
real(kind=FT) :: direction_vector(3),v(3)

logical Yes

P1(1:3) = line(1,1:3)
P2(1:3) = line(2,1:3)

P2_P1 = P2-P1

direction_vector = (P2-P1)/(sqrt(P2_P1(1)**2 + P2_P1(2)**2 + P2_P1(3)**2))
v = Point-P1

t = dot_product(v,direction_vector)

nearestPoint = P1+t*direction_vector
distance = sqrt((point(1)-nearestPoint(1))**2 +&
                (point(2)-nearestPoint(2))**2 +&
                (point(3)-nearestPoint(3))**2)

dis_to_P1 = sqrt((point(1)-P1(1))**2 +(point(2)-P1(2))**2 +(point(3)-P1(3))**2)
dis_to_P2 = sqrt((point(1)-P2(1))**2 +(point(2)-P2(2))**2 +(point(3)-P2(3))**2)

call Tool_Yes_Point_on_Line_Segment_3D(P1,P2,nearestPoint,Yes)

if(Yes .eqv. .false.) then
    if(dis_to_P1<=dis_to_P2)then
        nearestPoint = P1
        distance = dis_to_P1
    else
        nearestPoint = P2
        distance = dis_to_P2
    endif
endif 

return 
end SUBROUTINE Tool_Nearest_Point_from_Point_to_Segment_3D                         
