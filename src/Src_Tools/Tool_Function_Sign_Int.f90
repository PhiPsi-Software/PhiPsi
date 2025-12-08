 
pure integer function Tool_Function_Sign_Int(c_value)

use Global_Float_Type
implicit none
real(kind=FT),intent(in)::c_value


if(c_value > Tol_15) then
    Tool_Function_Sign_Int = 1
    return
elseif (c_value < -Tol_15) then
    Tool_Function_Sign_Int = -1
    return
else
    Tool_Function_Sign_Int = 0
    return
endif

end function Tool_Function_Sign_Int                        
