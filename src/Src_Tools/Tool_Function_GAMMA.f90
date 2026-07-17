!-----------------------------------------------------------
! Brief: Wrapper for the intrinsic Gamma function
!
! Parameters:
!   Input:  X  - argument of Gamma
!   Output: GA - computed Gamma(X)
!
! Notes:   Thin wrapper around the compiler GAMMA intrinsic.
!-----------------------------------------------------------

subroutine Tool_Function_GAMMA(X,GA)
!     ================================
!     Purpose: Compute gamma function
!     ================================
use Global_Float_Type
real(kind=FT),intent(in)::X
real(kind=FT),intent(out)::GA

GA = GAMMA(X)


return
end SUBROUTINE Tool_Function_GAMMA
