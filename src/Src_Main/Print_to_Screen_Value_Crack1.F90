!-----------------------------------------------------------
! Brief: Print aperture, pressure, or flux of crack 1 to screen.
!
! Parameters:
!   Input:  num_Space - leading spaces for indentation (unused logic)
!           c_String  - selector ('WIDTH', 'PRESS', or 'FLUXX')
!
! Notes:   Selects the requested crack-1 quantity and prints the maximum
!          value in the current unit system (SI or mm-ton-MPa).
!-----------------------------------------------------------

subroutine Print_to_Screen_Value_Crack1(num_Space,c_String)
use Global_Float_Type
use Global_Common 
use Global_Crack
implicit none
integer,intent(in):: num_Space
character(len=5) ::  c_String
real(kind=FT),allocatable::tem_vector_real(:)

if(c_String=='WIDTH')then
    if(Key_Unit_System==1)then
        print *,'    Max aperture of crack 1(mm):', maxval(Cracks_CalP_Aper(1,1:Cracks_CalP_Num(1))) *1000.0D0
    elseif(Key_Unit_System==2)then
        print *,'    Max aperture of crack 1(mm):', maxval(Cracks_CalP_Aper(1,1:Cracks_CalP_Num(1)))
    endif
elseif(c_String=='PRESS')then    
    if(Key_Unit_System==1)then
        print *,'    Max pressure of crack 1(MPa):', maxval(Cracks_CalP_Pres(1,1:Cracks_CalP_Num(1))) /1.0D6
    elseif(Key_Unit_System==2)then
        print *,'    Max pressure of crack 1(MPa):', maxval(Cracks_CalP_Pres(1,1:Cracks_CalP_Num(1)))
    endif
elseif(c_String=='FLUXX')then  
    if(Key_Unit_System==1)then
        print *,'    Max flux of crack 1(m*2/s):', maxval(Cracks_CalP_Quan(1,1:Cracks_CalP_Num(1)))
    elseif(Key_Unit_System==2)then
        print *,'    Max pressure of crack 1(m*2/s):', maxval(Cracks_CalP_Quan(1,1:Cracks_CalP_Num(1)))/1.0D6
    endif
else
    print *, '    Error :: wrong input par for ' //'Print_to_Screen_Value_Crack1'
    call Warning_Message('S',Keywords_Blank) 
endif

return 
end subroutine Print_to_Screen_Value_Crack1                
