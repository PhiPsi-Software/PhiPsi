 
subroutine Matrix_n_x_3_Check_Sort_Int(Matrix, n, check_matrix, count_error)
implicit none
integer, intent(in) :: Matrix(n,3)
integer, intent(in) :: n
logical, intent(out) :: check_matrix
integer, intent(out):: count_error

integer :: i

check_matrix = .true.

do i = 2, n
    if (Matrix(i, 1) < Matrix(i-1, 1)) then
        check_matrix = .false.
        count_error = n
        return
    end if
end do


count_error =0

do i = 2, n
    if (Matrix(i, 1) == Matrix(i-1, 1) .and. Matrix(i, 2) < Matrix(i-1, 2)) then
        check_matrix = .false.
        count_error = count_error +1
        
    end if
end do
end subroutine Matrix_n_x_3_Check_Sort_Int
    
    
    
 