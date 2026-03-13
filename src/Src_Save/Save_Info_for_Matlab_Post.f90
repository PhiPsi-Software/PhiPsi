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
 
SUBROUTINE Save_Info_for_Matlab_Post
!
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Filename
use Global_DISP
use Global_POST
use Global_Elem_Area_Vol

implicit none

character(200) c_File_name_1
logical alive

if (Key_Save_Nothing==1) return 

print *,'    Saving information for matlab post processing...'
c_File_name_1   =  trim(Full_Pathname)//'.post'

inquire(file=c_File_name_1, exist=alive)  
if(alive.EQV..True.)then
    OPEN  (UNIT=105, FILE=c_File_name_1, STATUS='OLD')
    CLOSE (UNIT=105, STATUS='DELETE')
endif

open(301,file=c_File_name_1,status='unknown',position='append',action='write')
write(301,'(A88)') 'Analysis type | Crack-tip type | Data format | Key_H_Value | Key_Hole_Value | Ave_Elem_L'

write(301, '(6X,I3,10X,I3,13X,I3,13X,I3,13X,I3,5X,F12.5)')       &
                  Key_Analysis_Type,Key_TipEnrich,        &
                  Key_Data_Format,Key_Heaviside_Value,Key_Hole_Value,Ave_Elem_L
close(301)

RETURN
END SUBROUTINE Save_Info_for_Matlab_Post
