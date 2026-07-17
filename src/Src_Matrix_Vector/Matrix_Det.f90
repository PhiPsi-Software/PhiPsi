!-----------------------------------------------------------
! Brief: Compute the determinant of a square matrix via Gaussian elimination.
!
! Parameters:
!   Input:  n      - order of the matrix
!           Matrix - input square matrix (n,n)
!   Output: det    - determinant value (0 if the matrix is singular)
!
! Notes:   Uses row pivoting; sign of det flips with each pivot row
!          swap. Low efficiency, intended for small matrices only.
!-----------------------------------------------------------

subroutine Matrix_Det(n,Matrix,det)   
! Find the determinant of the matrix.
! Note: The calculation efficiency is very low. 2022-08-18.

use Global_Float_Type       
implicit none
integer :: n
real(kind=FT),intent(in)::Matrix(n,n)
real(kind=FT),intent(out)::det

real(kind=FT) m,temp
real(kind=FT) temp_Matrix(n,n)

integer :: i, j, k, l

logical :: DetExists = .true.
l = 1

temp_Matrix = Matrix

do k = 1, n-1
    if (temp_Matrix(k,k) .eq. ZR) then 
        DetExists = .false.            
        do i = k+1, n 
            if (temp_Matrix(i,k) /= ZR) then
                do j = 1, n 
                    temp = temp_Matrix(i,j)
                    temp_Matrix(i,j)= temp_Matrix(k,j)
                    temp_Matrix(k,j) = temp
                end do
                DetExists = .true.
                l=-l
                exit
            end if
        end do
        if (DetExists .eqv. .false.) then
            det = 0
            return
        end if
    end if
    do j = k+1, n
        m = temp_Matrix(j,k)/temp_Matrix(k,k)
        do i = k+1, n
            temp_Matrix(j,i) = temp_Matrix(j,i) - m*temp_Matrix(k,i)
        end do
    end do
end do  

det = l
do i = 1, n
    det = det * temp_Matrix(i,i)
end do

return
END subroutine Matrix_Det



