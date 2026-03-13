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
 
SUBROUTINE Tool_Fit_Least_Square_Straight_Line(x, y, n, k, b)
! The least squares fit of the data is a straight line.
!
! 2023-08-22.
!
use Global_Float_Type
implicit none
integer, intent(in) :: n
real(kind=FT), intent(in) :: x(n), y(n)
real(kind=FT), intent(out) :: k, b
real(kind=FT) :: sum_x, sum_y, sum_xy, sum_x2
integer :: i

sum_x = ZR
sum_y = ZR
sum_xy = ZR
sum_x2 = ZR

do i = 1, n
    sum_x = sum_x + x(i)
    sum_y = sum_y + y(i)
    sum_xy = sum_xy + x(i) * y(i)
    sum_x2 = sum_x2 + x(i) * x(i)
end do

k = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x)
b = (sum_y - k * sum_x) / n

RETURN
END SUBROUTINE Tool_Fit_Least_Square_Straight_Line
