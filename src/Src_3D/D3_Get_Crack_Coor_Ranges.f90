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
 
SUBROUTINE D3_Get_Crack_Coor_Ranges(i_C,Coor_Ranges)
! Obtain the coordinate range of discrete fracture surfaces for the specified 3D fractures.
! 2023-01-09.

!-----------------------------
! Read public variable module
!-----------------------------
use Global_Float_Type
use Global_Common
use Global_Crack_Common
use Global_Crack_3D   
      
!---------------------------
! Variable Type Declaration
!---------------------------
implicit none
integer i_C
real(kind=FT),intent(out):: Coor_Ranges(3,2)
integer c_mesh_ndoes_num

!--------------------------
! Extract coordinate range
!--------------------------
c_mesh_ndoes_num = Crack3D_Meshed_Node_num(i_C)  
! x range
Coor_Ranges(1,1) =  minval(Crack3D_Meshed_Node(i_C)%row(1:c_mesh_ndoes_num,1))
Coor_Ranges(1,2) =  maxval(Crack3D_Meshed_Node(i_C)%row(1:c_mesh_ndoes_num,1))     
! y-range
Coor_Ranges(2,1) =  minval(Crack3D_Meshed_Node(i_C)%row(1:c_mesh_ndoes_num,2))
Coor_Ranges(2,2) =  maxval(Crack3D_Meshed_Node(i_C)%row(1:c_mesh_ndoes_num,2))   
! z range
Coor_Ranges(3,1) =  minval(Crack3D_Meshed_Node(i_C)%row(1:c_mesh_ndoes_num,3))
Coor_Ranges(3,2) =  maxval(Crack3D_Meshed_Node(i_C)%row(1:c_mesh_ndoes_num,3))           

RETURN
END SUBROUTINE D3_Get_Crack_Coor_Ranges
