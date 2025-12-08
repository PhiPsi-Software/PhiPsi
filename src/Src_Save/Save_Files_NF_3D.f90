 
SUBROUTINE Save_Files_NF_3D(isub,Save_Type)

use Global_Float_Type
use Global_Common
use Global_Crack_Common
use Global_Crack_3D
use Global_Model
use Global_Filename
use Global_POST

implicit none
integer isub,Save_Type

integer i,j
character(200) Filename_1,Filename_2,Filename_3
character(5) temp    


if (Key_Save_Nothing==1) return 

selectcase (Save_Type)

case(1)
    print *,'    Saving coors files of 3D NFs...'
    write(temp,'(I5)') isub
    Filename_1 = trim(Full_Pathname)//'.nfcx'
    Filename_2 = trim(Full_Pathname)//'.nfcy'
    Filename_3 = trim(Full_Pathname)//'.nfcz'  
    open(101,file=Filename_1,status='unknown')     
    open(102,file=Filename_2,status='unknown')     
    open(103,file=Filename_3,status='unknown')   
    do i=1,num_Rand_Na_Crack
        write(101,'(100E20.12)')  (Na_Crack3D_Coor(i,j,1),j=1,Each_NaCr3D_Poi_Num(i))
        write(102,'(100E20.12)')  (Na_Crack3D_Coor(i,j,2),j=1,Each_NaCr3D_Poi_Num(i))
        write(103,'(100E20.12)')  (Na_Crack3D_Coor(i,j,3),j=1,Each_NaCr3D_Poi_Num(i))
    end do
    close(101)      
    close(102)   
    close(103)   
case(2)
end select
      
RETURN
END SUBROUTINE Save_Files_NF_3D
