!-----------------------------------------------------------
! Brief: Compute the inverse of a 3x3 matrix by the cofactor formula.
!
! Parameters:
!   Input:  Matrix_A  - input 3x3 matrix
!   Output: Matrix_invA - inverse matrix (3,3)
!
! Notes:   Computes the cofactor matrix and divides by the
!          determinant; assumes the matrix is non-singular.
!-----------------------------------------------------------

SUBROUTINE Matrix_Inverse_3x3(Matrix_A,Matrix_invA)   
! Find the inverse matrix of the 3x3 matrix A.
! https://ardoris.wordpress.com/2008/07/18/general-formula-for-the-inverse-of-a-3x3-matrix/    

use Global_Float_Type
implicit none
real(kind=FT),intent(in)::Matrix_A(3,3)
real(kind=FT),intent(out)::Matrix_invA(3,3)

real(kind=FT) tem

tem =  Matrix_A(1,1)*(Matrix_A(2,2)*Matrix_A(3,3) -Matrix_A(2,3)*Matrix_A(3,2))- &
Matrix_A(1,2)*(Matrix_A(2,1)*Matrix_A(3,3) -Matrix_A(2,3)*Matrix_A(3,1))+ &
Matrix_A(1,3)*(Matrix_A(2,1)*Matrix_A(3,2) -Matrix_A(2,2)*Matrix_A(3,1))

Matrix_invA(1,1) = Matrix_A(2,2)*Matrix_A(3,3) -Matrix_A(2,3)*Matrix_A(3,2)
Matrix_invA(2,1) = Matrix_A(2,3)*Matrix_A(3,1) -Matrix_A(2,1)*Matrix_A(3,3)
Matrix_invA(3,1) = Matrix_A(2,1)*Matrix_A(3,2) -Matrix_A(2,2)*Matrix_A(3,1)
Matrix_invA(1,2) = Matrix_A(1,3)*Matrix_A(3,2) -Matrix_A(1,2)*Matrix_A(3,3)
Matrix_invA(2,2) = Matrix_A(1,1)*Matrix_A(3,3) -Matrix_A(1,3)*Matrix_A(3,1)
Matrix_invA(3,2) = Matrix_A(1,2)*Matrix_A(3,1) -Matrix_A(1,1)*Matrix_A(3,2)
Matrix_invA(1,3) = Matrix_A(1,2)*Matrix_A(2,3) -Matrix_A(1,3)*Matrix_A(2,2)
Matrix_invA(2,3) = Matrix_A(1,3)*Matrix_A(2,1) -Matrix_A(1,1)*Matrix_A(2,3)
Matrix_invA(3,3) = Matrix_A(1,1)*Matrix_A(2,2) -Matrix_A(1,2)*Matrix_A(2,1)
Matrix_invA = Matrix_invA/tem

return
END SUBROUTINE Matrix_Inverse_3x3



