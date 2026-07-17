!-----------------------------------------------------------
! Brief: Saves per-step linear-solver CPU time for fracture iterations.
!
! Parameters:
!   Input: ifra_CPU_time - linear-solver CPU time per fracture step
!   Input: Num_Frac      - number of fracture steps
!
! Notes:   Output file is <Full_Pathname>.icpt_LS for paper-data extraction.
!-----------------------------------------------------------

subroutine Save_Ifrac_LS_CPU_time(ifra_CPU_time,Num_Frac)
!     After all fracture steps are calculated, save the time consumed
!     by the linear solver for each fracture step (for extracting paper data).

use Global_Float_Type      
use Global_Common
use Global_Filename
use Global_POST

implicit none

real(kind=FT),intent(in):: ifra_CPU_time(Num_Frac)
integer,intent(in):: Num_Frac
integer :: i
character(200) c_File_name_1

if (Key_Save_Nothing==1) return 

print *,'    Saving LS CPU time (s) of each time step...'
c_File_name_1   =  trim(Full_Pathname)//'.icpt_LS'   
open(201,file=c_File_name_1,status='unknown')   
write(201,*) '    LS CPU time (s) of each fracturing step:' 
do i=1,Num_Frac
    write(201, '(I10,E15.5)') i,ifra_CPU_time(i)
end do
close(201)    

return
END subroutine Save_Ifrac_LS_CPU_time
