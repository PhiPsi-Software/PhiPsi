!-----------------------------------------------------------
! Brief: Computes the von Mises equivalent stress for a full 3D stress tensor.
!
! Parameters:
!   Input:  stress_x,stress_y,stress_z - normal stress components
!   Input:  stress_xy,stress_yz,stress_xz - shear stress components
!   Output: stress_VM - von Mises equivalent stress
!
! Notes:   Applies the standard Mises formula directly to the six independent
!          components of the 3D stress tensor.
!-----------------------------------------------------------

subroutine Tool_von_Mises_3D(stress_x,stress_y,stress_z, stress_xy,stress_yz,stress_xz, stress_VM)
! Calculate the Mises stress for a three-dimensional problem.

use Global_Float_Type
use Global_Common

implicit none
real(kind=FT),intent(in)::stress_x,stress_y,stress_z
real(kind=FT),intent(in)::stress_xy,stress_yz,stress_xz
real(kind=FT),intent(out)::stress_VM


stress_VM = sqrt(((stress_x-stress_y)**2 + (stress_y-stress_z)**2 + (stress_x-stress_z)**2 + &
SIX*(stress_xy**2+stress_yz**2+stress_xz**2))/TWO)

return
end SUBROUTINE Tool_von_Mises_3D
