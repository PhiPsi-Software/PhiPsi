!-----------------------------------------------------------
! Brief: Compute the arc length from an injection point to a calculation point
!
! Parameters:
!   Input:  i_Crack  - Crack index
!   Input:  CalP_x, CalP_y - Coordinates of the calculation point (informational)
!   Input:  i_CalP   - Calculation point index along the crack
!   Output: Length_Out - Cumulative length from the injection point to i_CalP
!
! Notes:   Walks the calculation points in both directions from the local
!          injection point; supports both full and symmetric-half HF models.
!-----------------------------------------------------------

subroutine Tool_Length_CalP_to_InjP(i_Crack,Length_Out, CalP_x,CalP_y,i_CalP)
!     Calculate the distance from the point to the water injection point
!     Calculate the point number based on the water injection point (this method is very convenient)
use Global_Float_Type      
use Global_Common
use Global_Crack
use Global_Crack_Common
use Global_HF

implicit none
integer,intent(in)::i_Crack,i_CalP
real(kind=FT),intent(in)::CalP_x,CalP_y
real(kind=FT),intent(out)::Length_Out
real(kind=FT) :: delta_L
integer :: i
real(kind=FT) CalP_L_x,CalP_L_y,CalP_R_x,CalP_R_y

if(i_Crack>num_Crack)then
    print *,'     Error:: wrong crack number for ' // 'Tool_Length_of_Crack.f!'
    call Warning_Message('S',Keywords_Blank) 
endif

Length_Out = ZR

if(Key_Symm_HF==0)then
    if(i_CalP>CalP_num_InjP_Local)then
        do i=CalP_num_InjP_Local,i_CalP-1
            CalP_L_x=Cracks_CalP_Coors(i_Crack,i,1)
            CalP_L_y=Cracks_CalP_Coors(i_Crack,i,2)
            CalP_R_x=Cracks_CalP_Coors(i_Crack,i+1,1)
            CalP_R_y=Cracks_CalP_Coors(i_Crack,i+1,2)
            delta_L=sqrt((CalP_L_x-CalP_R_x)**2 + (CalP_L_y-CalP_R_y)**2)
            Length_Out = Length_Out + delta_L
        end do
    elseif(i_CalP<CalP_num_InjP_Local)then
        do i=i_CalP,CalP_num_InjP_Local-1
            CalP_L_x=Cracks_CalP_Coors(i_Crack,i,1)
            CalP_L_y=Cracks_CalP_Coors(i_Crack,i,2)
            CalP_R_x=Cracks_CalP_Coors(i_Crack,i+1,1)
            CalP_R_y=Cracks_CalP_Coors(i_Crack,i+1,2)
            delta_L=sqrt((CalP_L_x-CalP_R_x)**2 + (CalP_L_y-CalP_R_y)**2)
            Length_Out = Length_Out + delta_L
        end do
    else
        Length_Out=ZR
    end if
end if
if(Key_Symm_HF==1)then
    if(i_CalP>0)then
        do i=1,i_CalP-1
            CalP_L_x=Cracks_CalP_Coors(i_Crack,i,1)
            CalP_L_y=Cracks_CalP_Coors(i_Crack,i,2)
            CalP_R_x=Cracks_CalP_Coors(i_Crack,i+1,1)
            CalP_R_y=Cracks_CalP_Coors(i_Crack,i+1,2)
            delta_L=sqrt((CalP_L_x-CalP_R_x)**2 + (CalP_L_y-CalP_R_y)**2)
            Length_Out = Length_Out + delta_L
        end do
    else
        Length_Out=ZR
    end if
end if

return 
end subroutine Tool_Length_CalP_to_InjP                         
