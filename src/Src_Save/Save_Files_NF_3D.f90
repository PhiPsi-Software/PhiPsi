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
 
SUBROUTINE Save_Files_NF_3D(isub,Save_Type)
! Save 3D Natural Fractures. NEWFTU2023010902.
! Save_Type=1: Saves vertex coords of fractures (quadrilateral/polygonal only).

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
