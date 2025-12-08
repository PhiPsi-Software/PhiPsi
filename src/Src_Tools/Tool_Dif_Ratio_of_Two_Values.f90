 
subroutine Tool_Dif_Ratio_of_Two_Values(Value1,Value2,Method_Key,&
                                        Threshold,Logical_Dif)

use Global_Float_Type 
implicit none
integer,intent(in) :: Method_Key
real(kind=FT),intent(in) :: Value1,Value2,Threshold
logical,intent(out) :: Logical_Dif

Logical_Dif = .False.

if(Method_Key==1) then
    if(abs(abs(Value1)-abs(Value2))<=Threshold)then
        Logical_Dif = .True.
    endif
endif

return 
end SUBROUTINE Tool_Dif_Ratio_of_Two_Values
              
