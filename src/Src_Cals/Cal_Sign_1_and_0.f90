!-----------------------------------------------------------
! Brief: Return 1 when variable is positive, else 0.
!
! Parameters:
!   Input:  Variable - real value to evaluate
!   Output: Sign_O   - 1.0 if Variable > 0, 0.0 otherwise
!-----------------------------------------------------------

subroutine Cal_Sign_1_and_0(Variable,Sign_O)

use Global_Float_Type     
implicit none

real(kind=FT),intent(in)::Variable
real(kind=FT),intent(out)::Sign_O

if(Variable.eq.ZR) then
    Sign_O = ZR
elseif (Variable > ZR) then
    Sign_O = ONE
elseif (Variable < ZR) then
    Sign_O = ZR      
end if

return 
end subroutine Cal_Sign_1_and_0               
