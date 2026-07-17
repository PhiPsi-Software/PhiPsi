!-----------------------------------------------------------
! Brief: Compute the inverse of a 2x2 matrix by closed-form formula.
!
! Parameters:
!   Input:  Matrix_A  - input 2x2 matrix
!   Output: Matrix_invA - inverse matrix (2,2)
!
! Notes:   Division by the determinant is performed inline; assumes
!          the input matrix is non-singular.
!-----------------------------------------------------------

subroutine Matrix_Inverse_2x2(Matrix_A,Matrix_invA)   
!     Find the inverse matrix of a 2x2 matrix A.

use Global_Float_Type  
implicit none
real(kind=FT),intent(in)::Matrix_A(2,2)
real(kind=FT),intent(out)::Matrix_invA(2,2)

real(kind=FT) :: tem

tem = Matrix_A(1,2)*Matrix_A(2,1) - Matrix_A(1,1)*Matrix_A(2,2) 

Matrix_invA(1,1) = -Matrix_A(2,2) /tem
Matrix_invA(1,2) =  Matrix_A(1,2)/tem
Matrix_invA(2,1) =  Matrix_A(2,1)/tem
Matrix_invA(2,2) = -Matrix_A(1,1)/tem      

return
END subroutine Matrix_Inverse_2x2



