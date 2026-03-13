 

subroutine Tool_Voigt_to_Tensor_Symmetric(Voigt, Tensor)

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

