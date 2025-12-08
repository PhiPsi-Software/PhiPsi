 
subroutine Matrix_n_x_3_Check_Sort(Matrix, n, check_matrix, count_error)
    use Global_Float_Type
    implicit none
    real(kind=FT), intent(in) :: Matrix(n,3)
    integer, intent(in) :: n
    logical, intent(out) :: check_matrix
    integer, intent(out) :: count_error

    integer :: i

    check_matrix = .true.

    do i = 2, n
        if (int(Matrix(i, 1)) < int(Matrix(i-1, 1))) then
            check_matrix = .false.
            count_error = n
            return
        end if
    end do
    
    count_error =0
    
    do i = 2, n
        if (int(Matrix(i, 1)) == int(Matrix(i-1, 1)) .and. int(Matrix(i, 2)) < int(Matrix(i-1, 2))) then
            check_matrix = .false.
            count_error = count_error +1
            
        end if
    end do
end subroutine Matrix_n_x_3_Check_Sort
    
    
    
 