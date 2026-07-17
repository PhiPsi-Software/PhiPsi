!-----------------------------------------------------------
! Brief: Compute the determinant of a 2x2 matrix by direct formula.
!
! Parameters:
!   Input:  Matrix - input 2x2 matrix
!   Output: det    - determinant (a11*a22 - a12*a21)
!-----------------------------------------------------------

subroutine Matrix_Det_2x2(Matrix,det)   
!     Find the determinant of a square matrix.
use Global_Float_Type       
implicit none
real(kind=FT),intent(in)::Matrix(2,2)
real(kind=FT),intent(out)::det

det = Matrix(1,1)*Matrix(2,2) - Matrix(1,2)*Matrix(2,1)

return
END subroutine Matrix_Det_2x2



