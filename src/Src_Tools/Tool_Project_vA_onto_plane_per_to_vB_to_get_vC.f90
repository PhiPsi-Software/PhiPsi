!-----------------------------------------------------------
! Brief: Project a 3D vector onto the plane perpendicular to another.
!
! Parameters:
!   Input:  vA - Vector to be projected (3-vector)
!   Input:  vB - Normal defining the projection plane (3-vector)
!   Output: vC - Projection of vA onto the plane perpendicular to vB
!
! Notes:   vC = vA - (vA.vB / vB.vB) * vB; returns vA unchanged if
!   vB is the zero vector (length below tolerance).
!-----------------------------------------------------------

subroutine Tool_Project_vA_onto_plane_per_to_vB_to_get_vC(vA, vB, vC)
! Project vector A onto the plane perpendicular to vector B to obtain vector C.
! 
! Theory:
! The projection of vA onto vB is: proj_vB(vA) = (vA·vB / |vB|^2) * vB
!   The component of vA perpendicular to vB is: vC = vA - proj_vB(vA)
!   This vC lies in the plane perpendicular to vB.
!
! Input:
!   vA(3)  - Vector A to be projected
!   vB(3)  - Vector B (normal to the plane)
! Output:
!   vC(3)  - Projected vector C (lying in the plane perpendicular to vB)
!
! Created: 2026-02-01. NEWFTU-2026020102.
!

!......................
! Variable Declaration
!......................
use Global_Float_Type

implicit none
real(kind=FT), intent(in)  :: vA(3)
real(kind=FT), intent(in)  :: vB(3)
real(kind=FT), intent(out) :: vC(3)
real(kind=FT) :: dot_AB, dot_BB
real(kind=FT) :: proj_vB(3)
real(kind=FT) :: vB_magnitude

dot_AB = vA(1)*vB(1) + vA(2)*vB(2) + vA(3)*vB(3)
dot_BB = vB(1)*vB(1) + vB(2)*vB(2) + vB(3)*vB(3)

vB_magnitude = sqrt(dot_BB)
if (vB_magnitude < Tol_15) then
    vC(1:3) = vA(1:3)
    return
endif

proj_vB(1:3) = (dot_AB / dot_BB) * vB(1:3)

vC(1:3) = vA(1:3) - proj_vB(1:3)

return
end subroutine Tool_Project_vA_onto_plane_per_to_vB_to_get_vC