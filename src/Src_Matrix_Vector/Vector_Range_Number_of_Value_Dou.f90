!-----------------------------------------------------------
! Brief: Find the index of the interval containing a given real value.
!
! Parameters:
!   Input:  value_to_be_find - scalar to locate
!           n                - number of values (interval endpoints)
!           values           - ascending-ordered endpoint vector
!   Output: Yes_In           - true if value falls within [values(1), values(n)]
!           Range_Num        - index i of the interval [values(i), values(i+1)]
!-----------------------------------------------------------

subroutine Vector_Range_Number_of_Value_Dou(value_to_be_find,n, values,Yes_In, Range_Num)
!     Get the position of the value in the vector (which segment it is in). The vector should be ensured to be in ascending order.
!     For example, value_to_be_find: 7.0 is at position 2 in values: 1.0, 5.0, 9.0, 11.1. Range_Num = 2.
!     2022-07-06.

use Global_Float_Type
implicit none

integer,intent(in)::n
real(kind=FT),intent(in)::value_to_be_find,values(n)
integer,intent(out)::Range_Num
logical,intent(out)::Yes_In
integer :: i

Yes_In = .False.

do i=1,n-1
    if(value_to_be_find >= values(i)   .and. value_to_be_find <= values(i+1)  )then
        Range_Num = i
        Yes_In = .True.
        return
    endif
end do

return 
end subroutine Vector_Range_Number_of_Value_Dou



