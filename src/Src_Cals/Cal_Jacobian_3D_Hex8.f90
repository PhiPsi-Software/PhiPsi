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
 

subroutine Cal_Jacobian_3D_Hex8(X, Y, Z, dN_dxi, dN_deta, dN_dzeta, J_mat, det_J)
!2026-01-28.

use Global_Float_Type      
implicit none

real(kind=FT), intent(in) :: X(8), Y(8), Z(8)
real(kind=FT), intent(in) :: dN_dxi(8), dN_deta(8), dN_dzeta(8)
real(kind=FT), intent(out) :: J_mat(3,3), det_J

J_mat(1,1) = sum(dN_dxi * X)
J_mat(1,2) = sum(dN_dxi * Y)
J_mat(1,3) = sum(dN_dxi * Z)
J_mat(2,1) = sum(dN_deta * X)
J_mat(2,2) = sum(dN_deta * Y)
J_mat(2,3) = sum(dN_deta * Z)
J_mat(3,1) = sum(dN_dzeta * X)
J_mat(3,2) = sum(dN_dzeta * Y)
J_mat(3,3) = sum(dN_dzeta * Z)

det_J = J_mat(1,1)*(J_mat(2,2)*J_mat(3,3) - J_mat(2,3)*J_mat(3,2)) - &
        J_mat(1,2)*(J_mat(2,1)*J_mat(3,3) - J_mat(2,3)*J_mat(3,1)) + &
        J_mat(1,3)*(J_mat(2,1)*J_mat(3,2) - J_mat(2,2)*J_mat(3,1))
    
end subroutine Cal_Jacobian_3D_Hex8
