!-----------------------------------------------------------
! Brief: Compute the Euclidean (L2) norm of a real vector.
!
! Parameters:
!   Input:  n      - vector length
!           Vector - real vector of length n
!   Output: Norm_2 - Euclidean norm sqrt(sum(Vector^2))
!-----------------------------------------------------------

subroutine Vector_Norm2(n,Vector,Norm_2)   
!     The Euclidean norm of a vector (square root of the sum of squares).
use Global_Float_Type
implicit none
integer,intent(in):: n
real(kind=FT),intent(in):: Vector(n)
real(kind=FT),intent(out)::Norm_2

real(kind=FT) :: tem

tem = sum(Vector(1:n)**2)
Norm_2  = sqrt(tem)

return
END subroutine Vector_Norm2



