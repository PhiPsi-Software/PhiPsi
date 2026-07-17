!-----------------------------------------------------------
! Brief: Compute the sign of a real number (-1, 0, or +1).
!
! Parameters:
!   Input:  Variable - real value to evaluate
!   Output: Sign_O   - +1 if Variable > 0, -1 if < 0, 0 otherwise
!-----------------------------------------------------------

subroutine Cal_Sign(Variable,Sign_O)

use Global_Float_Type     
implicit none

real(kind=FT),intent(in)::Variable
real(kind=FT),intent(out)::Sign_O

if(Variable.eq.ZR) then
    Sign_O = ZR
elseif (Variable > ZR) then
    Sign_O = ONE
elseif (Variable < ZR) then
    Sign_O = -ONE      
end if

return 
end subroutine Cal_Sign                
