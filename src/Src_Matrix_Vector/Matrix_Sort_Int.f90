!-----------------------------------------------------------
! Brief: Sort each row of an integer matrix in ascending order.
!
! Parameters:
!   Input:  m      - number of matrix rows
!           n      - number of matrix columns
!   In/Out: Matrix - integer matrix (m x n); each row sorted
!                    independently in-place
!
! Notes:   Per-row insertion sort.
!-----------------------------------------------------------

subroutine Matrix_Sort_Int(m,n,Matrix)   
!     Matrix sorting, integers, sort each row.
use Global_Float_Type    
implicit none
integer i,j,m,n
integer :: p
integer Matrix(m,n),a

do j=2,n
    do p=1,m
        a=Matrix(p,j)
        do i=j-1,1,-1
            if (Matrix(p,i).le.a) goto 10
            Matrix(p,i+1)=Matrix(p,i)
        end do
        i=0
        10         Matrix(p,i+1)=a
    end do
end do

return
END subroutine Matrix_Sort_Int



