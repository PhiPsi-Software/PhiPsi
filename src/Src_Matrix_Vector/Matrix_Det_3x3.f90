!-----------------------------------------------------------
! Brief: Compute the determinant of a 3x3 matrix by the Sarrus formula.
!
! Parameters:
!   Input:  Matrix - input 3x3 matrix
!   Output: det    - determinant value
!-----------------------------------------------------------

subroutine Matrix_Det_3x3(Matrix,det)   
!     Find the determinant of the matrix.
use Global_Float_Type       
implicit none
real(kind=FT),intent(in)::Matrix(3,3)
real(kind=FT),intent(out)::det

det = Matrix(1,1)*Matrix(2,2)*Matrix(3,3) + Matrix(1,2)*Matrix(2,3)*Matrix(3,1) + Matrix(2,1)*Matrix(3,2)*Matrix(1,3) - &
Matrix(1,1)*Matrix(2,3)*Matrix(3,2) - Matrix(1,2)*Matrix(2,1)*Matrix(3,3) - Matrix(1,3)*Matrix(2,2)*Matrix(3,1)

return
END subroutine Matrix_Det_3x3



