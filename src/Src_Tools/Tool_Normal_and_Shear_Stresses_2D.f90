!-----------------------------------------------------------
! Brief: Compute normal and shear stress on an inclined 2D section
!
! Parameters:
!   Input:  S_xx, S_yy, S_xy - In-plane stress components
!   Input:  theta_stress     - Section inclination angle
!   Output: S_normal         - Normal stress on the inclined plane
!   Output: S_shear          - Shear stress on the inclined plane
!
! Notes:   Applies the standard stress-rotation formulas for plane stress;
!          the principal-stress expression in the source comments is documentation only.
!-----------------------------------------------------------

subroutine Tool_Normal_and_Shear_Stresses_2D(S_xx,S_yy,S_xy,theta_stress,S_normal,S_shear)
!     This function calculates the normal stress and shear stress on a cross-section.
use Global_Float_Type
use Global_Common

implicit none
real(kind=FT),intent(in)::S_xx,S_yy,S_xy,theta_stress
real(kind=FT),intent(out)::S_normal,S_shear


S_normal=(S_xx+S_yy)/TWO + (S_xx-S_yy)/TWO*cos(TWO*theta_stress)-S_xy*sin(TWO*theta_stress)
S_shear =(S_xx-S_yy)/TWO*sin(TWO*theta_stress)+S_xy*cos(TWO*theta_stress)      
RETURN
END SUBROUTINE Tool_Normal_and_Shear_Stresses_2D