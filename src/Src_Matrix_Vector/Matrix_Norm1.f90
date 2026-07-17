!-----------------------------------------------------------
! Brief: Compute the entrywise 1-norm of a square matrix.
!
! Parameters:
!   Input:  n      - order of the matrix
!           Matrix - input matrix (n,n)
!   Output: Norm_1 - sum of absolute values of all entries
!-----------------------------------------------------------

subroutine Matrix_Norm1(n,Matrix,Norm_1)   
!     The 1-norm of a matrix (the sum of absolute values).
!     2022-07-09.

use Global_Float_Type
implicit none
integer,intent(in):: n
real(kind=FT),intent(in):: Matrix(n,n)
real(kind=FT),intent(out)::Norm_1

Norm_1  = sum(abs(Matrix(1:n,1:n)))   


return
END subroutine Matrix_Norm1



