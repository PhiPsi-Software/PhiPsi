!-----------------------------------------------------------
! Brief: Verify a real n x 3 matrix is sorted by (col1, col2).
!
! Parameters:
!   Input:  n           - number of matrix rows
!           Matrix      - real matrix to inspect (n x 3)
!   Output: check_matrix - .true. if ascending; .false. otherwise
!           count_error  - number of secondary-key violations
!                            (set to n if primary key fails)
!
! Notes:   Column 1 must be non-decreasing; within ties, column 2
!          must be non-decreasing. Uses INT-truncated comparison.
!-----------------------------------------------------------

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
    
    
    
 