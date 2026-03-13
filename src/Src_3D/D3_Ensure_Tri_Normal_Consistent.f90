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

subroutine D3_Ensure_Tri_Normal_Consistent(i_C, node1, node2, node3, refNormal)
! Ensure triangle node ordering gives consistent normal direction
! If the normal is opposite to refNormal, swap node2 and node3.
!
! 2026-02-14.
!
use Global_Float_Type
use Global_Common
use Global_Crack_3D
implicit none

integer, intent(in) :: i_C
integer, intent(inout) :: node1, node2, node3
real(kind=FT), intent(in) :: refNormal(3)

real(kind=FT) p1(3), p2(3), p3(3)
real(kind=FT) v1(3), v2(3), n(3)
real(kind=FT) dotv
integer :: tmp

p1 = Crack3D_Meshed_Node(i_C)%row(node1,1:3)
p2 = Crack3D_Meshed_Node(i_C)%row(node2,1:3)
p3 = Crack3D_Meshed_Node(i_C)%row(node3,1:3)

v1 = p2 - p1
v2 = p3 - p1
call Vector_Cross_Product_3(v1, v2, n)

dotv = n(1)*refNormal(1) + n(2)*refNormal(2) + n(3)*refNormal(3)

if(dotv < ZR) then
    tmp = node2
    node2 = node3
    node3 = tmp
endif

end subroutine D3_Ensure_Tri_Normal_Consistent