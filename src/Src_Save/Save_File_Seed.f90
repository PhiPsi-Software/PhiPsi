 
SUBROUTINE Save_File_Seed(c_Seed)

use Global_Common
use Global_Filename
use Global_POST

implicit none
integer,intent(in):: c_Seed
character(200) c_File_name_1

if (Key_Save_Nothing==1) return 

c_File_name_1   =  trim(Full_Pathname)//'.seed'   

open(301,file=c_File_name_1,status='unknown',action='write')
write(301,*) '----  Randonmly generated seed ----'  
write(301, '(I11)') Seed 
close(301)   

RETURN
END SUBROUTINE Save_File_Seed
