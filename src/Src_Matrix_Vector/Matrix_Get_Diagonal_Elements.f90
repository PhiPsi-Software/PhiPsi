!-----------------------------------------------------------
! Brief: Extract the diagonal entries of a square matrix into a vector.
!
! Parameters:
!   Input:  n        - order of the matrix
!           Matrix   - input square matrix (n,n)
!   Output: vector_D - extracted diagonal values (size n)
!-----------------------------------------------------------

subroutine Matrix_Get_Diagonal_Elements(n,Matrix,vector_D)   
!     Extract the diagonal elements of the matrix and store them in vector_D.
!     2020-02-11.

use Global_Float_Type       
implicit none
integer,intent(in)::n
real(kind=FT),intent(in)::Matrix(n,n)
real(kind=FT),intent(out)::vector_D(n)
integer :: i

do i = 1, n
    vector_D(i) = Matrix(i,i)
end do

return
END subroutine Matrix_Get_Diagonal_Elements



