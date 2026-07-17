!-----------------------------------------------------------
! Brief: Compute the entrywise 2-norm of a square matrix.
!
! Parameters:
!   Input:  n      - order of the matrix
!           Matrix - input matrix (n,n)
!   Output: Norm_2 - sqrt of sum of squared entries (Frobenius norm)
!-----------------------------------------------------------

subroutine Matrix_Norm2(n,Matrix,Norm_2)   
!     The matrix's 2-norm (square root of the sum of squares).

use Global_Float_Type
implicit none
integer,intent(in):: n
real(kind=FT),intent(in):: Matrix(n,n)
real(kind=FT),intent(out)::Norm_2

Norm_2  = sqrt(sum(Matrix(1:n,1:n)**2))   

return
END subroutine Matrix_Norm2



