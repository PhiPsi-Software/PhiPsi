!-----------------------------------------------------------
! Brief: Sets a logical flag indicating whether two scalar values differ within a threshold.
!
! Parameters:
!   Input:  Value1, Value2 - values to compare
!   Input:  Method_Key     - currently only method 1 (absolute-value comparison)
!   Input:  Threshold      - allowed difference
!   Output: Logical_Dif    - true if |abs(V1)-abs(V2)| <= Threshold
!
! Notes:   Designed for use in convergence tests where the magnitude matters
!          more than the sign of the residuals.
!-----------------------------------------------------------

subroutine Tool_Dif_Ratio_of_Two_Values(Value1,Value2,Method_Key, Threshold,Logical_Dif)
! Calculate the ratio of the difference between two values. NEWFTU2023022701.
!2023-02-27.
! Method_Key=1, compare absolute values.

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
              
