 
SUBROUTINE Save_current_folder_dat_file

use Global_Common
use Global_Filename

implicit none

character(200) c_File_name_1,c_File_name_2

if (Key_Clear_All==0) return 

print *,'    Saving current_folder.dat file...'

if(Operation_System_Type==1)then
    c_File_name_1   =  trim(ADJUSTL(PhiPsi_Current_Directory))//'/python_tools/current_folder.dat'   
    open(301,file=c_File_name_1,status='unknown') 
    write(301, *) trim(ADJUSTL(work_directory))//'/'
    close(301)
    
    c_File_name_2   =  trim(ADJUSTL(PhiPsi_Current_Directory))//'/current_folder.dat'   
    open(302,file=c_File_name_2,status='unknown') 
    write(302, *) trim(ADJUSTL(work_directory))//'/'
    close(302)   
elseif(Operation_System_Type==2)then
    c_File_name_1   =  trim(ADJUSTL(PhiPsi_Current_Directory))//'/python_tools/current_folder.dat'   
    open(301,file=c_File_name_1,status='unknown') 
    write(301, *) trim(ADJUSTL(work_directory))//'/'
    close(301)   
    
    c_File_name_2   =  trim(ADJUSTL(PhiPsi_Current_Directory))//'/current_folder.dat'   
    open(302,file=c_File_name_2,status='unknown') 
    write(302, *) trim(ADJUSTL(work_directory))//'/'
    close(302)   
elseif(Operation_System_Type==3)then
    c_File_name_1   =  trim(ADJUSTL(PhiPsi_Current_Directory))//'\python_tools\current_folder.dat'   
    open(301,file=c_File_name_1,status='unknown') 
    write(301, *) trim(ADJUSTL(work_directory))//'\'
    close(301)   
    
    c_File_name_2   =  trim(ADJUSTL(PhiPsi_Current_Directory))//'\current_folder.dat'   
    open(302,file=c_File_name_2,status='unknown') 
    write(302, *) trim(ADJUSTL(work_directory))//'\'
    close(302) 
elseif(Operation_System_Type==0)then
endif

 

RETURN
END SUBROUTINE Save_current_folder_dat_file
