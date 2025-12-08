 
function Tool_Function_Odd_Positive_Even_Negative(c_number)

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
