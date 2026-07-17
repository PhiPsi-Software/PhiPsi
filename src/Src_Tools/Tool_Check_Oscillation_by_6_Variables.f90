!-----------------------------------------------------------
! Brief: Detects whether a 6-sample sequence oscillates between two or three repeated values.
!
! Parameters:
!   Input:  Input_Var - 6-element array of samples
!   Output: Yes_Oscill - true if pattern A,b,A,b,A,b or A,B,C,A,B,C is detected
!
! Notes:   Applies two tolerance-driven rules: equal groups within Yita1 and
!          inter-group difference exceeding Yita2 to flag oscillation.
!-----------------------------------------------------------

subroutine Tool_Check_Oscillation_by_6_Variables(Input_Var, Yes_Oscill)
!     This subroutine check if the input six variables vary in an oscillating form.
!     In other wrods, check if the six variables arrange like this:
!     A,b,A,b,A,b  or  A,B,C,A,B,C

!     The rule for determining whether there is volatility is:
!     First case: for: A, b, A, b, A, b
!     (1) 1, 3, and 5 are basically the same (with an error not exceeding Yita1, e.g., 5%), and 2, 4, and 6 are basically the same (with an error not exceeding Yita1, e.g., 5%)
!     (2) On the premise of satisfying rule (1), 1, 3, and 5 should have significant differences from 2, 4, and 6, with a difference greater than Yita2 (e.g., 10%)
!      Specifically, compare the sum of 1, 3, and 5 with the sum of 2, 4, and 6.
!     Second case: for: A, B, C, A, B, C
!     (1) 1 and 4 are basically the same (with an error not exceeding Yita1, e.g., 5%), and 2 and 5 are basically the same (with an error not exceeding Yita1, e.g., 5%).
!     , 3, and 6 are basically the same (with an error not exceeding Yita1 (e.g., 5%))
!     (2) On the premise of meeting rule (1), there must be a significant difference between at least two of A, B, and C.
use Global_Float_Type      
implicit none
real(kind=FT),intent(in):: Input_Var(6)
logical,intent(out):: Yes_Oscill

real(kind=FT) Error_1_3_5,Error_2_4_6
real(kind=FT) Yita1,Yita2
real(kind=FT) sum_1_3_5,sum_2_4_6,Error_sum
real(kind=FT) Error_1_4,Error_2_5,Error_3_6
real(kind=FT) Error_1_2,Error_2_3,Error_1_3

Yes_Oscill = .False.
Yita1 = 5.0D-2;
Yita2 = 10.0D-2;

call Tool_Error_Between_3_Variables(Input_Var(1), Input_Var(3), Input_Var(5), Error_1_3_5)
call Tool_Error_Between_3_Variables(Input_Var(2), Input_Var(4), Input_Var(6), Error_2_4_6)
if(Error_1_3_5 <= Yita1  .and. Error_2_4_6 <= Yita1 )then
    sum_1_3_5 = Input_Var(1) + Input_Var(3) + Input_Var(5)
    sum_2_4_6 = Input_Var(2) + Input_Var(4) + Input_Var(6)
    Error_sum = abs(sum_1_3_5- sum_2_4_6 ) / max(sum_1_3_5, sum_2_4_6 )
    if(Error_sum >=Yita2 )then
        Yes_Oscill = .True.
        return
    endif
endif

Error_1_4 = abs(Input_Var(1)- Input_Var(4) ) / max(Input_Var(1), Input_Var(4) )
Error_2_5 = abs(Input_Var(2)- Input_Var(5) ) / max(Input_Var(2), Input_Var(5) )
Error_3_6 = abs(Input_Var(3)- Input_Var(6) ) / max(Input_Var(3), Input_Var(6) )
if(Error_1_4 <= Yita1  .and. Error_2_5 <= Yita1  .and. Error_3_6 <= Yita1 )    then
    Error_1_2 = abs(Input_Var(1)- Input_Var(2) ) / max(Input_Var(1), Input_Var(2) )
    Error_2_3 = abs(Input_Var(2)- Input_Var(3) ) / max(Input_Var(2), Input_Var(3) )
    Error_1_3 = abs(Input_Var(1)- Input_Var(3) ) / max(Input_Var(1), Input_Var(3) )
    if(Error_1_2 >=Yita2  .or. Error_2_3 >=Yita2  .or. Error_1_3 >=Yita2 )then
        Yes_Oscill = .True.
        return
    endif
endif

return 
end subroutine Tool_Check_Oscillation_by_6_Variables                         
