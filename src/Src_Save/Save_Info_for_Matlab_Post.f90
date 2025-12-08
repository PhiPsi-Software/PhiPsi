 
SUBROUTINE Save_Info_for_Matlab_Post
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Filename
use Global_DISP
use Global_POST
use Global_Elem_Area_Vol

implicit none

character(200) c_File_name_1
logical alive

if (Key_Save_Nothing==1) return 

print *,'    Saving information for matlab post processing...'
c_File_name_1   =  trim(Full_Pathname)//'.post'

inquire(file=c_File_name_1, exist=alive)  
if(alive.EQV..True.)then
    OPEN  (UNIT=105, FILE=c_File_name_1, STATUS='OLD')
    CLOSE (UNIT=105, STATUS='DELETE')
endif

open(301,file=c_File_name_1,status='unknown',position='append',action='write')
write(301,'(A88)') 'Analysis type | Crack-tip type | Data format | Key_H_Value | Key_Hole_Value | Ave_Elem_L'

write(301, '(6X,I3,10X,I3,13X,I3,13X,I3,13X,I3,5X,F12.5)')       &
                  Key_Analysis_Type,Key_TipEnrich,        &
                  Key_Data_Format,Key_Heaviside_Value,Key_Hole_Value,Ave_Elem_L
close(301)

RETURN
END SUBROUTINE Save_Info_for_Matlab_Post
