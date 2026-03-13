 
subroutine Tool_Project_vA_onto_plane_per_to_vB_to_get_vC(vA, vB, vC)

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