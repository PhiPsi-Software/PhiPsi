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
 
SUBROUTINE D3_Get_Cracks_Centers(Cracks_Coor_Centers)
! Obtain the central positions of each discrete fracture surface. NEWFTU2022050501.
! 2022-05-05.

!-----------------------------
! Read public variable module
!-----------------------------
use Global_Float_Type
use Global_Common
use Global_Crack_Common
use Global_Crack_3D
use omp_lib      
      
!---------------------------
! Variable Type Declaration
!---------------------------
implicit none
integer i_C
real(kind=FT),intent(out):: Cracks_Coor_Centers(num_Crack,3)
integer c_mesh_ndoes_num

!-------------
! Crack cycle
!-------------
Cracks_Coor_Centers(1:num_Crack,1:3) = ZR
! IMPROV2022072204. OpenMP Parallelization.
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,c_mesh_ndoes_num) &
!$OMP SCHEDULE(static) 
do i_C =1,num_Crack
  c_mesh_ndoes_num = Crack3D_Meshed_Node_num(i_C)  
  Cracks_Coor_Centers(i_C,1) = sum(Crack3D_Meshed_Node(i_C)%row(1:c_mesh_ndoes_num,1))/c_mesh_ndoes_num    
  Cracks_Coor_Centers(i_C,2) = sum(Crack3D_Meshed_Node(i_C)%row(1:c_mesh_ndoes_num,2))/c_mesh_ndoes_num    
  Cracks_Coor_Centers(i_C,3) = sum(Crack3D_Meshed_Node(i_C)%row(1:c_mesh_ndoes_num,3))/c_mesh_ndoes_num    
enddo
!$omp end parallel do       

RETURN
END SUBROUTINE D3_Get_Cracks_Centers
