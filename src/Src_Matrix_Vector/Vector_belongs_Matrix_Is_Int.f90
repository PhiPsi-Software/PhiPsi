!-----------------------------------------------------------
! Brief: Test whether an integer vector matches any row of a matrix.
!
! Parameters:
!   Input:  m        - number of rows in the matrix
!   Input:  n        - number of columns
!   Input:  Matrix   - m-by-n integer matrix
!   Input:  Vector   - candidate vector of length n
!   Output: Location - row index where match was found (0 if none)
!   Output: Yes     - .true. if Vector matches some row
!
! Notes:   Exact equality test using ALL(Vector_tem .eq. Vector);
!          scans rows and returns on the first match.
!-----------------------------------------------------------

subroutine Vector_belongs_Matrix_Is_Int(m,n,Matrix,Vector, Location,Yes)
!     Check whether the vector belongs to the matrix.  
use Global_Float_Type     
implicit none

integer,intent(in)::m,n
integer,intent(in)::Matrix(m,n)
integer,intent(in)::Vector(n)
integer,intent(out)::Location
logical,intent(out)::Yes
integer :: i
integer :: Vector_tem(n)

Location = 0
Yes = .False.
do i=1,m
    Vector_tem(1:n) = Matrix(i,1:n)
    if (ALL(Vector_tem .eq. Vector)) then
        Yes = .True.
        Location = i
        exit
    end if

end do

return 
end subroutine Vector_belongs_Matrix_Is_Int



