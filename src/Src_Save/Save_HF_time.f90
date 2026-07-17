!-----------------------------------------------------------
! Brief: Appends HF analysis time-history row to the .hftm log.
!
! Parameters:
!   Input: imf          - macro step index
!   Input: ifra         - fracture step index
!   Input: Counter_Iter - cumulative iteration count
!   Input: c_Time       - current physical time
!
! Notes:   Writes header and clears existing file at (1,1,1).
!-----------------------------------------------------------

subroutine Save_HF_time(imf,ifra,Counter_Iter,c_Time)
!     Save time vs. iteration step data for HF analysis.

use Global_Float_Type      
use Global_Common
use Global_Model
use Global_Filename
use Global_DISP
use Global_POST

implicit none

integer,intent(in)::imf,ifra,Counter_Iter
real(kind=FT),intent(in)::c_Time
character(200) c_File_name_1
logical :: alive

if (Key_Save_Nothing==1) return 

c_File_name_1   =  trim(Full_Pathname)//'.hftm'   

if(imf==1.and. ifra==1 .and. Counter_Iter==1)then
    inquire(file=c_File_name_1, exist=alive) 
    if(alive.eqv..True.)then
        open  (UNIT=105, FILE=c_File_name_1, STATUS='OLD') 
        close (UNIT=105, STATUS='DELETE')
    endif
endif

open(301,file=c_File_name_1,status='unknown', position='append',action='write')
if(imf==1.and.ifra==1.and.Counter_Iter==1)then
    write(301,*) '    imf   |   ifra   | total_ter|   time   '   
endif   
write(301, '(3I10,1F18.5)') imf,ifra,Counter_Iter,c_Time
close(301)    

return
END subroutine Save_HF_time
