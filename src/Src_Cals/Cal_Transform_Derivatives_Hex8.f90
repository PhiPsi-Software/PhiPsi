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
 

subroutine Cal_Transform_Derivatives_Hex8(dN_dxi, dN_deta, dN_dzeta, inv_J, dN_dx, dN_dy, dN_dz)
!This subroutine computes the spatial derivatives of the shape functions for an 8-node hexahedral
!(Hex8)
!finite element with respect to the global Cartesian coordinates (x, y, z). 
!2026-01-28.

use Global_Float_Type      
implicit none
real(kind=FT), intent(in) :: dN_dxi(8), dN_deta(8), dN_dzeta(8), inv_J(3,3)
real(kind=FT), intent(out) :: dN_dx(8), dN_dy(8), dN_dz(8)
integer :: i

do i = 1, 8
    dN_dx(i) = inv_J(1,1)*dN_dxi(i) + inv_J(1,2)*dN_deta(i) + inv_J(1,3)*dN_dzeta(i)
    dN_dy(i) = inv_J(2,1)*dN_dxi(i) + inv_J(2,2)*dN_deta(i) + inv_J(2,3)*dN_dzeta(i)
    dN_dz(i) = inv_J(3,1)*dN_dxi(i) + inv_J(3,2)*dN_deta(i) + inv_J(3,3)*dN_dzeta(i)
enddo
    
end subroutine Cal_Transform_Derivatives_Hex8
