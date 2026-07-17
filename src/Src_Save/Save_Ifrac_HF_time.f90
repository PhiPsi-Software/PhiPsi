!-----------------------------------------------------------
! Brief: Saves per-step HF (hydraulic-fracture) time in seconds.
!
! Parameters:
!   Input: ifra_HF_time - array of HF times per fracture step
!   Input: Num_Frac     - number of fracture steps
!
! Notes:   Output file is <Full_Pathname>.ihft; only positive entries written.
!-----------------------------------------------------------

subroutine Save_Ifrac_HF_time(ifra_HF_time,Num_Frac)
!     Save fracturing time (in seconds) for each fracture step
!     upon completion of all fracture calculations.

use Global_Float_Type      
use Global_Common
use Global_Filename
use Global_POST

implicit none

real(kind=FT),intent(in):: ifra_HF_time(Num_Frac)
integer,intent(in):: Num_Frac
integer :: i
character(200) c_File_name_1

if (Key_Save_Nothing==1) return

print *,'    Saving HF time (s) of each time step...'
c_File_name_1   =  trim(Full_Pathname)//'.ihft'   
open(201,file=c_File_name_1,status='unknown')   
write(201,*) '    HF time of each time step' 
do i=1,Num_Frac
    if(ifra_HF_time(i)>ZR)then
        write(201, '(F12.5)') ifra_HF_time(i)
    endif
end do
close(201)  

return
END subroutine Save_Ifrac_HF_time
