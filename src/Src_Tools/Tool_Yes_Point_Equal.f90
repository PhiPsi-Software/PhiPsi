!-----------------------------------------------------------
! Brief: Test whether two 2D points coincide within tolerance
!
! Parameters:
!   Input:  A(2),B(2) - point coordinates
!   Output: Yes_Equal - true if both components are within Tol_11
!
! Notes:   Component-wise distance test with Tol_11.
!-----------------------------------------------------------

subroutine Tool_Yes_Point_Equal(A,B,Yes_Equal)
!     Determine whether two points are the same (within a certain margin of error)
use Global_Float_Type

implicit none
real(kind=FT),intent(in):: A(2),B(2)
logical,intent(out):: Yes_Equal

Yes_Equal = .False.

if ((abs(A(1)-B(1)) <=Tol_11) .and.  (abs(A(2)-B(2)) <=Tol_11)) then
    Yes_Equal = .True. 
end if

return 
end SUBROUTINE Tool_Yes_Point_Equal                       
