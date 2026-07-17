!-----------------------------------------------------------
! Brief: Sort an n x 3 real matrix and sum rows with equal (i,j).
!
! Parameters:
!   Input:  n     - number of input rows
!           M     - real input matrix (n x 3); cols 1,2 = indices,
!                   col 3 = value
!   Output: newM  - real output matrix (n x 3) without duplicate
!                  (i,j) pairs; col 3 holds the summed value
!           new_n - number of valid rows in newM
!
! Notes:   Sorts by (col1, col2) using Matrix_n_x_2_Quick_Sort_Int
!          and validates with Matrix_n_x_3_Check_Sort. Rows with
!          identical (i,j) pairs are merged by accumulating col 3.
!-----------------------------------------------------------

subroutine Matrix_n_x_3_Combine_Same_i_j_Lines(n, M, newM, new_n)
!
! Sort an n-row, 3-column matrix.
! The first and second columns represent positions, and the third column represents values.
! The result is sorted in ascending order by the first column, and on this basis, the second column
! is sorted in ascending order.
!
! 2024-09.
!
use Global_Float_Type 
use module_INTERFACE_Matrix_n_x_3_Check_Sort

implicit none
integer, intent(in) :: n
real(kind=FT), intent(in) :: M(n, 3)
real(kind=FT), intent(out) :: newM(n, 3)
integer, intent(out) :: new_n
integer :: i
real(kind=FT) Sorted_M(n,3)
real(kind=FT) Sorted_M_2(n,3)
integer Matrix_M1_M2(n,2)
integer idx(n)
LOGICAL check_matrix
integer count_error
integer i_Try

Sorted_M = M
do i_Try = 1,100
    Matrix_M1_M2(1:n,1) =  INT(Sorted_M(1:n,1))
    Matrix_M1_M2(1:n,2) =  INT(Sorted_M(1:n,2))
    idx = [(i, i=1,n)]

    call Matrix_n_x_2_Quick_Sort_Int(Matrix_M1_M2(1:n,1:2),n,idx(1:n))  


    do i = 1, n
        Sorted_M_2(i,1:3) = Sorted_M(idx(i),1:3) 
    enddo

    Sorted_M  = Sorted_M_2

    call Matrix_n_x_3_Check_Sort(Sorted_M, n, check_matrix,count_error)

    if(check_matrix)then
        exit
    endif    
enddo  

if (check_matrix .eqv. .false.) then
    write ( *, '(a)' ) 'ERROR :: Failed to sort the matrix! In Matrix_n_x_3_Combine_Same_i_j_Lines.f90!'
    call Warning_Message('S',' ')
endif


new_n = 0
newM = 0.0D0
new_n = 1
newM(new_n, :) = Sorted_M(1, :)
do i = 2, n
    if (Sorted_M(i, 1) == newM(new_n, 1) .and. Sorted_M(i, 2) == newM(new_n, 2)) then
        newM(new_n, 3) = newM(new_n, 3) + Sorted_M(i, 3)
    else
        new_n = new_n + 1
        newM(new_n, :) = Sorted_M(i, :)
    end if
end do

end subroutine Matrix_n_x_3_Combine_Same_i_j_Lines
