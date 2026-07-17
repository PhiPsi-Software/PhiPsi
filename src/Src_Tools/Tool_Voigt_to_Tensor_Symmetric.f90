!-----------------------------------------------------------
! Brief: Convert a 6-component Voigt engineering-strain vector to a symmetric 3x3 tensor
!
! Parameters:
!   Input:  Voigt(6)  - Voigt engineering strain vector
!   Input:  Tensor(3,3) - placeholder for symmetric 3x3 tensor
!   Output: Tensor(3,3) - symmetric 3x3 strain tensor
!
! Notes:   Multiplies shear entries by 0.5 to undo engineering-strain factor.
!-----------------------------------------------------------

subroutine Tool_Voigt_to_Tensor_Symmetric(Voigt, Tensor)
!This subroutine converts a symmetric second-order tensor stored in Voigt notation 
!(a 6-component vector) into its full 3��3 matrix representation. 
!2026-01-28.

use Global_Float_Type

implicit none
real(kind=FT), intent(in) :: Voigt(6)
real(kind=FT), intent(out) :: Tensor(3,3)

Tensor(1,1) = Voigt(1)
Tensor(2,2) = Voigt(2)
Tensor(3,3) = Voigt(3)

Tensor(1,2) = HLF * Voigt(4)
Tensor(2,1) = Tensor(1,2)

Tensor(2,3) = HLF * Voigt(5)
Tensor(3,2) = Tensor(2,3)

Tensor(1,3) = HLF * Voigt(6)
Tensor(3,1) = Tensor(1,3)
end subroutine Tool_Voigt_to_Tensor_Symmetric

