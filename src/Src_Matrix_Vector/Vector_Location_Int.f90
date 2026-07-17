!-----------------------------------------------------------
! Brief: Locate a value in an integer vector by linear search.
!
! Parameters:
!   Input:  n        - vector length
!           Vector   - integer vector to search
!           Variable - value to find
!   Output: location - 1-based index of the first match (0 if absent)
!           Yes_In   - true if a match was found
!-----------------------------------------------------------

subroutine Vector_Location_Int(n,Vector,Variable,location,Yes_In)   
!     Check if Variable is in Vector; if it is, then Yes_In = .True. and return the position.

use Global_Float_Type          
implicit none
integer,intent(in)::n,Vector(n),Variable
integer,intent(out)::location
logical,intent(out):: Yes_In
integer :: i

Yes_In =.False.
location = 0
do i=1,n
    if (Vector(i) == Variable) then
        location = i
        Yes_In = .True.
        exit
    end if
end do  

return
END subroutine Vector_Location_Int



