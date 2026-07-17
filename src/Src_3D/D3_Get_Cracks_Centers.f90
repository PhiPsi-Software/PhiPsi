!-----------------------------------------------------------
! Brief: Compute centroids of all 3D crack meshed surfaces.
!
! Parameters:
!   Output: Cracks_Coor_Centers - per-crack centroids (num_Crack x 3)
!
! Notes:   Parallelised loop calling D3_Get_Crack_Center for every crack.
!-----------------------------------------------------------

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
