!-----------------------------------------------------------
! Brief: Allocate and fill bounding ranges for every 3D crack.
!
! Parameters:
!
! Notes:   (Re)allocates the global Crack_Coor_Range(num_Crack,3,2)
!          array and fills it with min/max of meshed node coordinates.
!          OpenMP-parallelised over cracks.
!-----------------------------------------------------------

SUBROUTINE D3_Get_Cracks_Coor_Ranges
! Obtain the coordinate ranges of each discrete fracture surface. NEWFTU2022050401.
! 2022-05-04.

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
integer c_mesh_ndoes_num

!--------------------------------
! Allocate Memory and Initialize
!--------------------------------
if (allocated(Crack_Coor_Range)) then
    deallocate(Crack_Coor_Range)
endif
allocate(Crack_Coor_Range(num_Crack,3,2))
Crack_Coor_Range(1:num_Crack,1:3,1:2) = ZR

!-------------
! Crack cycle
!-------------
! IMPROV2022072202. OpenMP Parallelization.
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,c_mesh_ndoes_num) SCHEDULE(static) 
do i_C =1,num_Crack
    c_mesh_ndoes_num = Crack3D_Meshed_Node_num(i_C)  
    ! x range
    Crack_Coor_Range(i_C,1,1) =  minval(Crack3D_Meshed_Node(i_C)%row(1:c_mesh_ndoes_num,1))
    Crack_Coor_Range(i_C,1,2) =  maxval(Crack3D_Meshed_Node(i_C)%row(1:c_mesh_ndoes_num,1))     
    ! y-range
    Crack_Coor_Range(i_C,2,1) =  minval(Crack3D_Meshed_Node(i_C)%row(1:c_mesh_ndoes_num,2))
    Crack_Coor_Range(i_C,2,2) =  maxval(Crack3D_Meshed_Node(i_C)%row(1:c_mesh_ndoes_num,2))   
    ! z range
    Crack_Coor_Range(i_C,3,1) =  minval(Crack3D_Meshed_Node(i_C)%row(1:c_mesh_ndoes_num,3))
    Crack_Coor_Range(i_C,3,2) =  maxval(Crack3D_Meshed_Node(i_C)%row(1:c_mesh_ndoes_num,3))        
enddo
!$omp end parallel do       

RETURN
END SUBROUTINE D3_Get_Cracks_Coor_Ranges
