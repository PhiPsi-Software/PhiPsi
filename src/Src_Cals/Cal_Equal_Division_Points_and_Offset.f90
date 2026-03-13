!     ...........................................
!             ____  _       _   ____  _____   _        
!            |  _ \| |     |_| |  _ \|  ___| |_|       
!            | |_) | |___   _  | |_) | |___   _        
!            |  _ /|  _  | | | |  _ /|___  | | |       
!            | |   | | | | | | | |    ___| | | |       
!            |_|   |_| |_| |_| |_|   |_____| |_|       
!     ...........................................
!     PhiPsi:     a general-purpose computational      
!                 mechanics program written in Fortran.
!     Website:    http://phipsi.top                    
!     Author:     Fang Shi  
!     Contact me: shifang@hyit.edu.cn     

subroutine Cal_Equal_Division_Points_and_Offset(Num_Diversion, Line_AB, &
offset_delta,Div_Points, Offsetted_D_P_Up, &
    Offsetted_D_P_Down)
!=========================================================================
! Get the equal diversion points of line AB, then offset them by a small delta.
! Diversion points are arranged from A to B.
! Note: Key_Include_Endpoint is fixed to 0
! 
! Input:
!   Num_Diversion        - Number of divisions
!   Line_AB(2,2)         - Line segment [a_x a_y; b_x b_y]
!   offset_delta         - Offset distance
!
! Output:
!   Div_Points(Num_Diversion-1, 2)          - Division points
!   Offsetted_D_P_Up(Num_Diversion-1, 2)    - Offset points to left
!   Offsetted_D_P_Down(Num_Diversion-1, 2)  - Offset points to right
!=========================================================================
use Global_Float_Type
implicit none

integer, intent(in) :: Num_Diversion
real(kind=FT), intent(in) :: Line_AB(2,2)
real(kind=FT), intent(in) :: offset_delta

real(kind=FT), intent(out) :: Div_Points(Num_Diversion-1, 2)
real(kind=FT), intent(out) :: Offsetted_D_P_Up(Num_Diversion-1, 2)
real(kind=FT), intent(out) :: Offsetted_D_P_Down(Num_Diversion-1, 2)

integer :: i, num_points
real(kind=FT) :: a_x, a_y, b_x, b_y
real(kind=FT) :: theta

a_x = Line_AB(1, 1)
a_y = Line_AB(1, 2)
b_x = Line_AB(2, 1)
b_y = Line_AB(2, 2)

num_points = Num_Diversion - 1

do i = 1, num_points
    Div_Points(i, 1) = (real(i, FT) * b_x + real(Num_Diversion - i, FT) * a_x) / real(Num_Diversion, FT)
    Div_Points(i, 2) = (real(i, FT) * b_y + real(Num_Diversion - i, FT) * a_y) / real(Num_Diversion, FT)
end do

theta = atan2(b_y - a_y, b_x - a_x)

Div_Points(1, 1) = Div_Points(1, 1) + offset_delta * cos(theta)
Div_Points(1, 2) = Div_Points(1, 2) + offset_delta * sin(theta)

Div_Points(num_points, 1) = Div_Points(num_points, 1) - offset_delta * cos(theta)
Div_Points(num_points, 2) = Div_Points(num_points, 2) - offset_delta * sin(theta)

do i = 1, num_points
    Offsetted_D_P_Up(i, 1) = Div_Points(i, 1) - offset_delta * sin(theta)
    Offsetted_D_P_Up(i, 2) = Div_Points(i, 2) + offset_delta * cos(theta)
end do

do i = 1, num_points
    Offsetted_D_P_Down(i, 1) = Div_Points(i, 1) + offset_delta * sin(theta)
    Offsetted_D_P_Down(i, 2) = Div_Points(i, 2) - offset_delta * cos(theta)
end do

return
end subroutine Cal_Equal_Division_Points_and_Offset