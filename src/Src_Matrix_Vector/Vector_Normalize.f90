!-----------------------------------------------------------
! Brief: Normalize a real vector in place to unit Euclidean length.
!
! Parameters:
!   Input:  n      - vector length
!   In/Out: Vector - vector to be normalized; replaced by Vector/||Vector||
!
! Notes:   Guards against division by zero by using a small tolerance
!          when the input norm is below the zero threshold.
!-----------------------------------------------------------

subroutine Vector_Normalize(n,Vector)   
!     Vector Normalization.

use Global_Float_Type
use Global_Common  

implicit none
integer,intent(in):: n
real(kind=FT),intent(inout):: Vector(n)
real(kind=FT) :: Norm_2

call Vector_Norm2(n,Vector,Norm_2)   
if (Norm_2==ZR)then
    Norm_2 = Tol_30
endif

Vector = Vector/Norm_2

return
END subroutine Vector_Normalize



