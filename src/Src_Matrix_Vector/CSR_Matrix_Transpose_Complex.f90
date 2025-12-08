 
subroutine CSR_Matrix_Transpose_Complex(m, n, nnz, values, col_ind, row_ptr, values_t, col_ind_t, row_ptr_t)
use Global_Float_Type  
implicit none

integer, intent(in) :: m, n
integer, intent(in) :: nnz
complex(kind=FT), intent(in) :: values(nnz)
integer, intent(in) :: col_ind(nnz)
integer, intent(in) :: row_ptr(m+1)

complex(kind=FT), intent(out) :: values_t(nnz)
integer, intent(out) :: col_ind_t(nnz)
integer, intent(out) :: row_ptr_t(n+1)

integer :: k
integer i
integer i_NZZ,j, count
integer, allocatable :: row_counts(:)

row_ptr_t = 0

allocate(row_counts(n))
row_counts = 0
do i_NZZ = 1, nnz
    row_counts(col_ind(i_NZZ)) = row_counts(col_ind(i_NZZ)) + 1
end do

row_ptr_t(1) = 1
do i = 1, n
    row_ptr_t(i+1) = row_ptr_t(i) + row_counts(i)
end do

row_counts = 0

do i = 1, m
    do j = row_ptr(i), row_ptr(i+1) - 1
        k = col_ind(j)
        count = row_ptr_t(k) + row_counts(k)
        values_t(count) = values(j)
        col_ind_t(count) = i
        row_counts(k) = row_counts(k) + 1
    end do
end do

deallocate(row_counts)
    
end subroutine CSR_Matrix_Transpose_Complex