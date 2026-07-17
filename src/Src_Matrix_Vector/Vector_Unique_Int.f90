!-----------------------------------------------------------
! Brief: Remove duplicate entries from an integer vector.
!
! Parameters:
!   Input:  n           - declared length of the input array
!   Input:  n_Op_F      - effective length to process (Vector(1:n_Op_F))
!   Input:  Vector      - integer source array
!   Output: Uniqued_Vec - integer output array with duplicates removed
!   Output: Uniqued_n   - number of unique entries returned
!
! Notes:   Same algorithm as the double precision version;
!          output buffer is zero-initialised before scanning.
!-----------------------------------------------------------

subroutine Vector_Unique_Int(n,n_Op_F,Vector, Uniqued_Vec,Uniqued_n)
!     Remove duplicate elements from the array, generate a new array, and store it in Uniqued_Vec, integer type
!     n: length of the original array
!     n_Op_F: The end operation position of the original array, that is, the Vector(1:n_Op_F) of the processed array
!     Uniqued_Vec: New Array
!     Uniqued_n: Length of usable data in the new array
use Global_Float_Type
implicit none
integer,intent(in):: n,n_Op_F
integer,intent(in):: Vector(n)
integer,intent(out)::Uniqued_Vec(n)
integer,intent(out)::Uniqued_n

integer i,j,k
Uniqued_Vec(1:n) = 0
k = 1
Uniqued_Vec(1) = Vector(1)

outer: do i=2,n_Op_F
do j=1,k
    if (Uniqued_Vec(j) .eq. Vector(i)) then
        cycle outer
    end if
end do
k = k + 1
Uniqued_Vec(k) = Vector(i)
end do outer

Uniqued_n   = k

return
END subroutine Vector_Unique_Int



