!-----------------------------------------------------------
! Brief: Find a value's position in an integer vector via linear search.
!
! Parameters:
!   Input:  n        - vector length
!           Vector   - integer vector to search
!           Variable - value to find
!   Output: location - 1-based index of the first match (0 if absent)
!-----------------------------------------------------------

subroutine Vector_Location_Int_v2(n,Vector,Variable,location)   
!     Check if a variable is in a vector, and if it is, return its position.
!     2022-08-16.

use Global_Float_Type          
implicit none
integer,intent(in)::n,Variable
integer, intent(in) :: Vector(n)
integer,intent(out)::location
integer :: i

location = 0
do i=1,n
    if (Vector(i) == Variable) then
        location = i
        return
    end if
end do  

return
END subroutine Vector_Location_Int_v2




