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
 
subroutine Tool_Delete_Collinear_Points_3D(input_points, input_point_num, output_points, output_point_num)
! Used to remove collinear points in a space. The first and last points are not deleted.
!Fang Shi.
!2024-05-03.
! Algorithm: If a point and the points on its left and right do not lie on the same line, then keep
! the middle point.
!
use Global_Float_Type
implicit none
integer, intent(in) :: input_point_num
integer :: output_point_num
real(kind=FT), dimension(input_point_num, 1:3), intent(in) :: input_points
real(kind=FT), dimension(input_point_num, 1:3), intent(out) :: output_points
integer :: i, c_count
logical :: Yes_Three_Points_Collinear
    
real(kind=FT) P1(3),P2(3),P3(3)

output_point_num = 1
c_count = 1


output_points(1, :) = input_points(1, :)

do i = 2, input_point_num-1


    P1(1:3) = input_points(i-1,1:3)
    P2(1:3) = input_points(i,1:3)
    P3(1:3) = input_points(i+1,1:3)
    
    call Tool_Yes_Three_Points_Collinear(p1, p2, p3, Yes_Three_Points_Collinear)
    
    
    if(Yes_Three_Points_Collinear .eqv. .False.) then
        c_count = c_count + 1
        output_points(c_count, :) = input_points(i, :)
    endif
    
end do

output_points(c_count+1, :) = input_points(input_point_num, :)

output_point_num = c_count + 1
end subroutine Tool_Delete_Collinear_Points_3D
