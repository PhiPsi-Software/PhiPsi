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


subroutine Cal_3D_SIFs_IIM_q_Line_Integral(i_Crack, i_Vertex, L_front, q_line_integral)
! Calculate the Q-letter line integral along the front of the crack.
! 2026-01-28.

use Global_Float_Type

implicit none
integer, intent(in) :: i_Crack    
integer, intent(in) :: i_Vertex    
real(kind=FT), intent(in) :: L_front 
real(kind=FT), intent(out) :: q_line_integral
integer :: i_seg, num_segments
real(kind=FT) :: s_start, s_end, s_mid, delta_s
real(kind=FT) :: q_start, q_end
real(kind=FT) :: integral_segment

q_line_integral = ZR

num_segments = 100  
delta_s = (FOU * L_front) / real(num_segments, FT)

do i_seg = 1, num_segments
    
    s_start = -TWO * L_front + real(i_seg-1, FT) * delta_s
    s_end   = -TWO * L_front + real(i_seg,   FT) * delta_s
    s_mid   =  ZP5 * (s_start + s_end)
    
    q_start = Compute_Q_Front_Component(s_start, L_front)
    q_end   = Compute_Q_Front_Component(s_end,   L_front)
    
    integral_segment = ZP5 * (q_start + q_end) * delta_s
    
    q_line_integral = q_line_integral + integral_segment
    
end do

contains

function Compute_Q_Front_Component(s, L_f) result(q_front)
    implicit none
    real(kind=FT), intent(in) :: s
    real(kind=FT), intent(in) :: L_f
    real(kind=FT) :: q_front
    if (abs(s) <= L_f) then
        q_front = 1.0D0
    else if (abs(s) >= 2.0D0 * L_f) then
        q_front = 0.0D0
    else
        q_front = (2.0D0 * L_f - abs(s)) / L_f
    end if
end function Compute_Q_Front_Component
    
    
end subroutine Cal_3D_SIFs_IIM_q_Line_Integral


