!     ================================================= !
!             ____  _       _   ____  _____   _         !
!            |  _ \| |     |_| |  _ \|  ___| |_|        !
!            | |_) | |___   _  | |_) | |___   _         !
!            |  _ /|  _  | | | |  _ /|___  | | |        !
!            | |   | | | | | | | |    ___| | | |        !
!            |_|   |_| |_| |_| |_|   |_____| |_|        !
!     ================================================= !
!     PhiPsi:     a general-purpose computational       !
!                 mechanics program written in Fortran. !
!     Website:    http://phipsi.top                     !
!     Author:     Shi Fang, Huaiyin Institute of        !
!                 Technology, Huaian, JiangSu, China    !
!     Email:      shifang@hyit.edu.cn                   !
!     ------------------------------------------------- !
!     Please cite the following papers:                 !
!     (1)Shi F., Lin C. Modeling fluid-driven           !
!        propagation of 3D complex crossing fractures   !
!        with the extended finite element method.       !
!        Computers and Geotechnics, 2024, 172, 106482.  !
!     (2)Shi F., Wang D., Li H. An XFEM-based approach  !
!        for 3D hydraulic fracturing simulation         !
!        considering crack front segmentation. Journal  !
!        of Petroleum Science and Engineering, 2022,    !
!        214, 110518.                                   !
!     (3)Shi F., Wang D., Yang Q. An XFEM-based         !
!        numerical strategy to model three-dimensional  !
!        fracture propagation regarding crack front     !
!        segmentation. Theoretical and Applied Fracture !
!        Mechanics, 2022, 118, 103250.                  !
!     (4)Shi F., Liu J. A fully coupled hydromechanical !
!        XFEM model for the simulation of 3D non-planar !
!        fluid-driven fracture propagation. Computers   !
!        and Geotechnics, 2021, 132: 103971.            !
!     (5)Shi F., Wang X.L., Liu C., Liu H., Wu H.A. An  !
!        XFEM-based method with reduction technique     !
!        for modeling hydraulic fracture propagation    !
!        in formations containing frictional natural    !
!        fractures. Engineering Fracture Mechanics,     !
!        2017, 173: 64-90.                              !
!     ------------------------------------------------- !
 
SUBROUTINE Save_current_folder_dat_file
! Save the current working path to the current_folder.dat file.
!2022-07-25.

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
