!-----------------------------------------------------------
! Brief: Swap two triangle vertices if the computed normal is inverted.
!
! Parameters:
!   Input:  i_C      - crack index
!           refNormal - reference normal vector (3 components)
!   In/Out: node1, node2, node3 - triangle vertex indices; node2/node3
!                                  swapped when dot product < 0
!
! Notes:   Computes (p2-p1) x (p3-p1) and compares against refNormal.
!-----------------------------------------------------------

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