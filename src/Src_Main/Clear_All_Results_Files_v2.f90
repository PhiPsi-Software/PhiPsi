!-----------------------------------------------------------
! Brief: Delete result files in the work directory.
!
! Parameters:
!
! Notes:   Invokes a Python helper script based on the host
!          operating system (Linux, macOS, or Windows).
!-----------------------------------------------------------

SUBROUTINE Clear_All_Results_Files_v2
!Delete all result files in the work directory before start phipsi.
!2022-07-25.

use global_common

print *,'    Clear result files in work directory using Python 3...'
!call system ('cd .\&&phipsi_results_file_deleter.py')   

if(Operation_System_Type==1)then
    call system ('./python_tools/phipsi_results_file_deleter.py') 
elseif(Operation_System_Type==2)then
    call system ('./python_tools/phipsi_results_file_deleter.py') 
elseif(Operation_System_Type==3)then
    call system ('cd .\python_tools\&&phipsi_results_file_deleter.py')    
elseif(Operation_System_Type==0)then
    !do nothing.
endif

RETURN
END SUBROUTINE Clear_All_Results_Files_v2
