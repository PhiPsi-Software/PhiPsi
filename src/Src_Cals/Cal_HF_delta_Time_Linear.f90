!-----------------------------------------------------------
! Brief: Estimate the HF time increment from mass balance.
!
! Parameters:
!   Input:  Counter           - current iteration counter
!           Last_Cr_CalP_Aper - previous-step apertures per calc pt
!           total_Time        - current accumulated time
!   Output: delta_Time        - estimated time step
!
! Notes:   Computes the volume change of each HF unit (using the
!          rectangular or trapezoidal rule, controlled by
!          Key_Cal_deltaTime) and divides by the injection rate
!          to obtain a mass-consistent delta_Time.
!-----------------------------------------------------------

subroutine Cal_HF_delta_Time_Linear(Counter,Last_Cr_CalP_Aper, total_Time,delta_Time)

! Calculate the delta_T of each water injection crack system.
! For linear HF elements.
use Global_Float_Type
use Global_Crack
use Global_Crack_Common
use Global_HF
use Global_Material

implicit none

integer,intent(in)::Counter
real(kind=FT),intent(in):: Last_Cr_CalP_Aper(Max_Num_Cr, &
Max_Num_Cr_CalP)
real(kind=FT),intent(in):: total_Time
real(kind=FT),intent(out)::delta_Time

integer i_C,i_HF_Elem,Num_Div_Points
real(kind=FT) :: Total_delta_Vol
real(kind=FT) delta_w_1,delta_w_2,delta_w,Inject_Q


Total_delta_Vol = ZR
do i_C = 1,num_Crack
    if (Cracks_HF_State(i_C) == 1) then   
        Num_Div_Points = Cracks_CalP_Num(i_C)
        do i_HF_Elem = 1,Num_Div_Points-1
            delta_w_1 =   Cracks_CalP_Aper(i_C,i_HF_Elem)- Last_Cr_CalP_Aper(i_C,i_HF_Elem)
            delta_w_2 =   Cracks_CalP_Aper(i_C,i_HF_Elem+1)- Last_Cr_CalP_Aper(i_C,i_HF_Elem+1)
            if (Key_Cal_deltaTime == 1) then
                delta_w = delta_w_1  
            elseif (Key_Cal_deltaTime == 2) then
                delta_w = HLF*(delta_w_1 + delta_w_2)
            end if
            Total_delta_Vol = Total_delta_Vol + delta_w * Cracks_HF_Ele_L(i_C,i_HF_Elem)
        end do
    end if
end do
call Tool_Get_Value_from_x_y_Curve(Inject_Q_Time,Inject_Q_Val, 200,total_Time,Inject_Q)

if(Key_Symm_HF==0)then
    delta_Time = Total_delta_Vol/Inject_Q
elseif(Key_Symm_HF==1)then
    delta_Time = Total_delta_Vol/(Inject_Q*HLF)
endif
return 
end subroutine Cal_HF_delta_Time_Linear               
