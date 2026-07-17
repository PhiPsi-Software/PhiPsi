!-----------------------------------------------------------
! Brief: Detect whether a real matrix has any duplicate rows.
!
! Parameters:
!   Input:  m       - number of matrix rows
!           n       - number of matrix columns
!           Matrix  - real input matrix (m x n)
!   Output: Yes_Dup - .true. if any two rows match within Tol_11
!           num_Dup - number of duplicate rows (m - Uniqued_m)
!
! Notes:   Built atop the row-uniqueness scan used by
!          Matrix_Unique_Row_Dou.
!-----------------------------------------------------------

subroutine Matrix_Yes_Duplicated_Row_Dou(m,n,Matrix,Yes_Dup, num_Dup)
!     Determine whether the matrix has duplicate rows.
!     Written based on Matrix_Unique_Row_Dou.f.
!     Added on 2022-04-23.

use Global_Float_Type
implicit none
integer,intent(in):: m,n
real(kind=FT),intent(in):: Matrix(m,n)
logical,intent(out):: Yes_Dup
integer,intent(out):: num_Dup
integer Uniqued_m,Uni_Mat_Count(m)
real(kind=FT) Uniqued_Matrix(m,n),tem(n)
integer i,j,k

Uniqued_Matrix(1:m,1:n) = ZR
Yes_Dup = .False.
Uni_Mat_Count(1:m)=1

k = 1
Uniqued_Matrix(1,:) = Matrix(1,:)

outer: do i=2,m
do j=1,k
    tem = Uniqued_Matrix(j,:) - Matrix(i,:)
    if ((abs(maxval(tem))<=Tol_11) .and. (abs(minval(tem))<=Tol_11)) then
        Uni_Mat_Count(j) = Uni_Mat_Count(j)+1
        Yes_Dup = .True.
        cycle outer
    end if
end do
k = k + 1
Uniqued_Matrix(k,:) = Matrix(i,:)
end do outer

Uniqued_m = k    

num_Dup = m-Uniqued_m

return
END subroutine Matrix_Yes_Duplicated_Row_Dou



