!-----------------------------------------------------------
! Brief: Remove duplicate rows from an integer matrix.
!
! Parameters:
!   Input:  m              - number of matrix rows
!           n              - number of matrix columns
!           m_Finish       - last row index to consider
!           Matrix         - integer input matrix (m x n)
!   Output: Uniqued_Matrix - de-duplicated matrix (m x n)
!           Uniqued_m      - number of unique rows kept
!           Uni_Mat_Count  - occurrence count of each kept row
!
! Notes:   Integer exact comparison; only the first m_Finish rows of
!          Matrix are scanned.
!-----------------------------------------------------------

subroutine Matrix_Unique_Row_Int(m,n,m_Finish,Matrix, Uniqued_Matrix,Uniqued_m, Uni_Mat_Count)
!     Operation target: from the 1st row to the m_Finish-th row of the matrix
!     Remove duplicate rows from the original matrix and store them in Uniqued_Matrix as integers.
!     m, n: Original matrix dimensions
!     Uniqued_Matrix: New Matrix
!     Uniqued_m: Number of usable data rows in the new matrix
!     Uni_Mat_Count: Number of times each row in the new matrix is repeated

use Global_Float_Type
implicit none
integer,intent(in):: m,n,m_Finish
integer,intent(in):: Matrix(m,n)
integer,intent(out):: Uniqued_m
integer,intent(out):: Uniqued_Matrix(m,n)
integer,intent(out):: Uni_Mat_Count(m)

integer :: tem(n)
integer i,j,k

Uniqued_Matrix(1:m,1:n) = 0

Uni_Mat_Count(1:m)=1

k = 1
Uniqued_Matrix(1,:) = Matrix(1,:)

outer: do i=2,m_Finish
do j=1,k
    tem = Uniqued_Matrix(j,:) - Matrix(i,:)
    if ((maxval(tem).eq. 0).and. (minval(tem).eq. 0)) then
        Uni_Mat_Count(j) = Uni_Mat_Count(j)+1
        cycle outer
    end if
end do
k = k + 1
Uniqued_Matrix(k,:) = Matrix(i,:)
end do outer

Uniqued_m = k            

return
END subroutine Matrix_Unique_Row_Int



