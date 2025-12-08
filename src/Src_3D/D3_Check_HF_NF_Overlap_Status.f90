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
 
SUBROUTINE D3_Check_HF_NF_Overlap_Status(HF_Coor_Ranges,NF_Coor_Ranges,Yes_Overlap)
! Used to detect the coordinate coincidence of hydraulic fractures and natural fractures.
! Only supports rectangular or polygonal natural fractures.
! 2023-01-09.

!*****************************
! Read public variable module
!*****************************
use Global_Float_Type

!**********************
! Variable Declaration
!**********************
implicit none
logical,intent(out)::Yes_Overlap
real(kind=FT),intent(in)::HF_Coor_Ranges(3,2),NF_Coor_Ranges(3,2)
logical x_Logical_Yes,y_Logical_Yes,z_Logical_Yes

!-----------
! Detection
!-----------
Yes_Overlap = .False.
call Tool_Yes_Two_Ranges_Overlapped_Double(HF_Coor_Ranges(1,1:2),NF_Coor_Ranges(1,1:2),x_Logical_Yes)  
if (x_Logical_Yes .eqv. .False.) then
    return
endif
call Tool_Yes_Two_Ranges_Overlapped_Double(HF_Coor_Ranges(2,1:2),NF_Coor_Ranges(2,1:2),y_Logical_Yes)  
if (y_Logical_Yes .eqv. .False.) then
    return
endif
call Tool_Yes_Two_Ranges_Overlapped_Double(HF_Coor_Ranges(3,1:2),NF_Coor_Ranges(3,1:2),z_Logical_Yes)  
if (z_Logical_Yes .eqv. .False.) then
    return
endif  
Yes_Overlap = .True.


RETURN
END SUBROUTINE D3_Check_HF_NF_Overlap_Status
