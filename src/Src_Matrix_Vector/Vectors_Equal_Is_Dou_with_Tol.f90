!-----------------------------------------------------------
! Brief: Test double precision vector equality with a tolerance.
!
! Parameters:
!   Input:  Vector_A - first double precision vector
!   Input:  Vector_B - second double precision vector
!   Input:  n        - length of the vectors
!   Output: Yes      - .true. if all element differences within tolerance
!
! Notes:   Per-element comparison uses Tol_11 absolute tolerance;
!          returns at the first element exceeding the tolerance.
!-----------------------------------------------------------

subroutine Vectors_Equal_Is_Dou_with_Tol(Vector_A,Vector_B,n,Yes)   
!     Compare whether two vectors are the same, double precision.
!     Modified based on Vectors_Equal_Is_Dou.f, with added tolerance.
!     2022-04-29.

use Global_Float_Type
implicit none
integer,intent(in):: n
real(kind=FT),intent(in):: Vector_A(n),Vector_B(n)
logical,intent(out)::Yes

integer :: i

Yes = .False.

do i=1,n
    if(abs(Vector_A(i)-Vector_B(i))>Tol_11) then
        return
    end if
end do

Yes = .True.

return
END subroutine Vectors_Equal_Is_Dou_with_Tol



