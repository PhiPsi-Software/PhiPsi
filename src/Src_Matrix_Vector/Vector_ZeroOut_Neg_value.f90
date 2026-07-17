!-----------------------------------------------------------
! Brief: Clamp negative entries in a vector to zero.
!
! Parameters:
!   Input:  n      - length of the array
!   In/Out: Vector - real array whose negative entries are replaced
!
! Notes:   Uses Fortran WHERE construct; relies on ZR zero
!          reference value from Global_Float_Type.
!-----------------------------------------------------------

subroutine Vector_ZeroOut_Neg_value(Vector,n)   
!     Clamp negative vector values to zero.
!      
use Global_Float_Type     
implicit none
integer,intent(in):: n
real(kind=FT),intent(inout)::Vector(n)

where (Vector<ZR)
    Vector=ZR
end where

return 
end subroutine Vector_ZeroOut_Neg_value



