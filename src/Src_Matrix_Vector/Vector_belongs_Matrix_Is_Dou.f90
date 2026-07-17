!-----------------------------------------------------------
! Brief: Test whether a double precision vector matches any row of a matrix.
!
! Parameters:
!   Input:  m        - number of rows in the matrix
!   Input:  n        - number of columns
!   Input:  Matrix   - m-by-n double precision matrix
!   Input:  Vector   - candidate vector of length n
!   Output: Location - row index where match was found (0 if none)
!   Output: Yes     - .true. if Vector matches some row
!
! Notes:   Comparison uses Tol_11 tolerance; element-wise loop
!          with early exit on first match.
!-----------------------------------------------------------

subroutine Vector_belongs_Matrix_Is_Dou(m,n,Matrix,Vector,Location,Yes)   
!Check whether the vector belongs to the matrix, double precision.        
!implicit none
!integer,intent(in)::m,n  ! m is the number of rows in the matrix, and n is the number of columns in the matrix
!real(kind=FT),intent(in)::Vector(n)!,Vector_tem(n)
!
!Yes = .False.
!
!do i=1,m
!  if ((abs(MaxVal(Matrix(i,1:n)-Vector)) < Tol_11).and.  &
!      (abs(MinVal(Matrix(i,1:n)-Vector)) < Tol_11)) then
!      Yes = .True.
!      Location = i
!      exit
!  end if
!end do


!IMPROV-2026040802.
use Global_Float_Type
implicit none

integer,       intent(in)  :: m, n
real(kind=FT), intent(in)  :: Matrix(m,n)
real(kind=FT), intent(in)  :: Vector(n)
integer,       intent(out) :: Location
logical,       intent(out) :: Yes

integer :: i, j
logical :: same_row

Yes = .False.
Location = 0

do i = 1, m
    same_row = .True.
    do j = 1, n
        if (abs(Matrix(i,j) - Vector(j)) >= Tol_11) then
            same_row = .False.
            exit
        end if
    end do

    if (same_row) then
        Yes = .True.
        Location = i
        exit
    end if
end do

return 
end subroutine Vector_belongs_Matrix_Is_Dou



