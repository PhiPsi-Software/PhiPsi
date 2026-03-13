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


subroutine Cal_3D_SIFs_IIM_q_Function(Elem_ID, q_nodal, N, dN_dx, dN_dy, dN_dz, &
                                  GP_coords, Tip_Coords, e1, e2, e3, &
                                  R_in, R_out, L_f, q_val, grad_q)
! PURPOSE: Evaluate q-function and its gradient at Gauss point.
! 2026-01-28.
!
! METHOD:
!   1. Interpolate q from nodal values (shape functions)
!   2. Compute gradient from nodal values and shape function derivatives
!   3. Can optionally recompute directly from GP position for accuracy
use Global_Float_Type
implicit none
integer, intent(in) :: Elem_ID
real(kind=FT), intent(in) :: q_nodal(:,:), N(8), dN_dx(8), dN_dy(8), dN_dz(8)
real(kind=FT), intent(in) :: GP_coords(3), Tip_Coords(3)
real(kind=FT), intent(in) :: e1(3), e2(3), e3(3)
real(kind=FT), intent(in) :: R_in, R_out, L_f
real(kind=FT), intent(out) :: q_val, grad_q(3)

real(kind=FT) :: Vec_to_GP(3), r_radial, s_along_front
real(kind=FT) :: Vec_radial(3), e_radial(3)
real(kind=FT) :: q_r, q_s, dq_r, dq_s,dist


 Vec_to_GP = GP_coords - Tip_Coords
 
 if (L_f > Tol_12) then
     s_along_front = dot_product(Vec_to_GP, e3)
     Vec_radial = Vec_to_GP - s_along_front * e3
     r_radial = norm2(Vec_radial)
     
     if (r_radial > Tol_10) then
         e_radial = Vec_radial / r_radial
     else
         e_radial = e1
     endif
     
     if (r_radial <= R_in) then
         q_r = ONE
         dq_r = ZR
     else if (r_radial >= R_out) then
         q_r = ZR
         dq_r = ZR
     else
         q_r = (R_out - r_radial) / (R_out - R_in)
         dq_r = -ONE / (R_out - R_in)
     endif
     
     if (abs(s_along_front) <= L_f) then
         q_s = ONE
         dq_s = ZR
     else if (abs(s_along_front) >= TWO*L_f) then
         q_s = ZR
         dq_s = ZR
     else
         q_s = (TWO*L_f - abs(s_along_front)) / L_f
         dq_s = -sign(ONE, s_along_front) / L_f
     endif
     q_val = q_r * q_s
     grad_q = dq_r * q_s * e_radial + q_r * dq_s * e3
    
 else
    dist = norm2(Vec_to_GP)
    if (dist <= R_in) then
        q_val = ONE
        grad_q = ZR  
    else if (dist >= R_out) then
        q_val = ZR
        grad_q = ZR
    else
        q_val = (R_out - dist) / (R_out - R_in)
        grad_q = -Vec_to_GP / (dist * (R_out - R_in))
    endif
 endif
end subroutine Cal_3D_SIFs_IIM_q_Function



