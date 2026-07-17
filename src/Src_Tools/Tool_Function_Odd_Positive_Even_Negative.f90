!-----------------------------------------------------------
! Brief: Sign function based on parity of an integer
!
! Parameters:
!   Input:  c_number - integer to inspect
!   Returns: real - +1 for odd, -1 for even
!-----------------------------------------------------------

function Tool_Function_Odd_Positive_Even_Negative(c_number)
! Odd numbers are positive (1), and even numbers are negative (-1).
! 2022-08-01.

use Global_Float_Type
implicit none
integer,intent(in)::c_number
real(kind=FT) :: Tool_Function_Odd_Positive_Even_Negative

if (mod(c_number,2)==0) then
    Tool_Function_Odd_Positive_Even_Negative =  -ONE
else
    Tool_Function_Odd_Positive_Even_Negative =   ONE
endif


return 
end function Tool_Function_Odd_Positive_Even_Negative                      
