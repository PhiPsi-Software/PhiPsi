!-----------------------------------------------------------
! Brief: Append the current load factor and crack-opening displacement.
!
! Parameters:
!   Input:  isub     - current sub-step index
!   Input:  c_Lambda - current load factor
!   Input:  c_crack  - crack identifier whose COD is recorded
!
! Notes:   Computes COD as the maximum aperture along crack 1 and
!          appends one row to the .fccu Force-COD curve file per step.
!-----------------------------------------------------------

subroutine Save_File_F_COD_Curve(isub,c_Lambda,c_crack)
!     Save the load-COD curve of a specific fracture at current step
!     into the fccu file.
!     2022-10-15.

use Global_Float_Type
use Global_Common
use Global_Crack
use Global_Model
use Global_Filename
use Global_DISP
use Global_POST

implicit none

integer,intent(in):: isub,c_crack
real(kind=FT),intent(in):: c_Lambda
integer i,j
character(200) c_File_name_1
logical :: alive
real(kind=FT) :: COD

if (Key_Save_Nothing==1) return 

COD =   maxval(Cracks_CalP_Aper(1,1:Cracks_CalP_Num(1)))


print *,'    Saving Force-COD curve file...'
c_File_name_1   =  trim(Full_Pathname)//'.fccu'   

if(isub==1)then
    inquire(file=c_File_name_1, exist=alive)  
    if(alive.eqv..True.)then
        open  (UNIT=105, FILE=c_File_name_1, STATUS='OLD') 
        close (UNIT=105, STATUS='DELETE')
    endif
endif

open(301,file=c_File_name_1,status='unknown', position='append',action='write')
if(isub==1)then
    write(301,*) ' crack  |  isub  |   lambda   |   COD'  
    write(301, '(2I7,2F15.7)')  c_crack,0,ZR,ZR
else
    write(301, '(2I7,2F15.7)')c_crack,isub,c_Lambda,COD
endif   

close(301)   

return
END subroutine Save_File_F_COD_Curve
