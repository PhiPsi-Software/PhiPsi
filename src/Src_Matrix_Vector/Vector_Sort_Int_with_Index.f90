!-----------------------------------------------------------
! Brief: Sort an integer vector and produce a permutation index.
!
! Parameters:
!   Input:  n       - length of the array
!   In/Out: Vector - integer array to be sorted in place
!   Output: V_Index - integer array giving source positions
!
! Notes:   Inserts then uses Vector_Location_Int to recover the
!          original positions; a simple self-contained index sort.
!-----------------------------------------------------------

subroutine Vector_Sort_Int_with_Index(n,Vector,V_Index)   
!     One-dimensional array sorting, integer.
!     The V_Index array stores the positions of the sorted elements in the original array.
!     2022-07-27.

use Global_Float_Type     
implicit none
integer,intent(in)::n
integer,intent(inout):: Vector(n)
integer,intent(out):: V_Index(n)
integer :: Vector_Old(n)
integer i,j,a
logical :: c_Yes_In

Vector_Old = Vector

do j=2, n
    a=Vector(j)
    do i=j-1,1,-1
        if (Vector(i).le.a) goto 10
        Vector(i+1)=Vector(i)


    end do
    i=0

    10     Vector(i+1)=a

end do


do i=1,n
    call Vector_Location_Int(n,Vector_Old,Vector(i), V_Index(i),c_Yes_In)
enddo


return
END subroutine Vector_Sort_Int_with_Index



