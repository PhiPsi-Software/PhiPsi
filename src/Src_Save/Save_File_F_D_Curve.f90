!-----------------------------------------------------------
! Brief: Append a force-displacement data point for a specific node.
!
! Parameters:
!   Input:  isub        - current sub-step index
!   Input:  c_Lambda    - current load factor
!   Input:  c_node_num  - target node whose displacement magnitude is saved
!
! Notes:   Computes the nodal displacement magnitude from the global
!          DISP vector and appends it to the .fdcu curve file each step.
!-----------------------------------------------------------

subroutine Save_File_F_D_Curve(isub,c_Lambda,c_node_num)
!     Save the load-displacement curve of specific nodes at current step
!     into the fdcu file.

use Global_Float_Type
use Global_Common
use Global_Crack
use Global_Model
use Global_Filename
use Global_DISP
use Global_POST

implicit none

integer,intent(in):: isub,c_node_num
real(kind=FT),intent(in):: c_Lambda
integer i,j
character(200) c_File_name_1
logical :: alive
real(kind=FT) :: dis_of_node

if (Key_Save_Nothing==1) return

dis_of_node = sqrt(DISP(2*c_node_num-1)**2 + DISP(2*c_node_num)**2)


print *,'    Saving Force-Displacement curve file...'
c_File_name_1   =  trim(Full_Pathname)//'.fdcu'   

if(isub==1)then
    inquire(file=c_File_name_1, exist=alive)
    if(alive.eqv..True.)then
        open  (UNIT=105, FILE=c_File_name_1, STATUS='OLD') 
        close (UNIT=105, STATUS='DELETE')
    endif
endif

open(301,file=c_File_name_1,status='unknown', position='append',action='write')
if(isub==1)then
    write(301,*) ' node |  isub | lambda | dis_of_node'  
    write(301, '(2I7,2E15.7)')c_node_num,0,ZR,ZR
else
    write(301, '(2I7,2E15.7)')c_node_num,isub, c_Lambda,dis_of_node
endif   

close(301)   



return
END subroutine Save_File_F_D_Curve
