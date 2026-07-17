!-----------------------------------------------------------
! Brief: Compute the outer product of two vectors into a matrix.
!
! Parameters:
!   Input:  Vector1   - first vector of length n1
!   Input:  n1        - length of Vector1
!   Input:  Vector2   - second vector of length n2
!   Input:  n2        - length of Vector2
!   Output: Matrix_OUT - resulting n1-by-n2 matrix
!
! Notes:   Matrix_OUT(i,j) = Vector1(i) * Vector2(j);
!          useful for rank-1 updates in matrix assembly.
!-----------------------------------------------------------

subroutine Vectors_Multi(Vector1,n1,Vector2,n2,Matrix_OUT)   
!     When two vectors are multiplied, a matrix is obtained. For example, (3*1) * (1*4) results in a 3*4 matrix.

use Global_Float_Type
implicit none
integer,intent(in):: n1,n2
real(kind=FT),intent(in)::Vector1(n1),Vector2(n2)
real(kind=FT),intent(out)::Matrix_OUT(n1,n2)

integer i,j

do i=1,n1
    do j=1,n2
        Matrix_OUT(i,j) = Vector1(i)*Vector2(j)
    end do
end do

return
END subroutine Vectors_Multi



