!-----------------------------------------------------------
! Brief: Compute centroid of a 3D crack's discrete meshed surface.
!
! Parameters:
!   Input:  i_C                - crack index
!   Output: Cracks_Coor_Center - centroid coordinates (x, y, z)
!
! Notes:   Average of Crack3D_Meshed_Node(i_C)%row(:,1:3) over the
!          meshed node count Crack3D_Meshed_Node_num(i_C).
!-----------------------------------------------------------

SUBROUTINE D3_Get_Crack_Center(i_C,Cracks_Coor_Center)
! Obtain the center position of the discrete fracture surface for a given crack. NEWFTU2022092902.
! 2022-09-29.

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
integer,intent(in):: i_C
real(kind=FT),intent(out):: Cracks_Coor_Center(3)
integer c_mesh_ndoes_num

!-------------
! Crack cycle
!-------------
Cracks_Coor_Center(1:3) = ZR
c_mesh_ndoes_num = Crack3D_Meshed_Node_num(i_C)  
Cracks_Coor_Center(1) =  sum(Crack3D_Meshed_Node(i_C)%row(1:c_mesh_ndoes_num,1))/c_mesh_ndoes_num    
Cracks_Coor_Center(2) =  sum(Crack3D_Meshed_Node(i_C)%row(1:c_mesh_ndoes_num,2))/c_mesh_ndoes_num    
Cracks_Coor_Center(3) =  sum(Crack3D_Meshed_Node(i_C)%row(1:c_mesh_ndoes_num,3))/c_mesh_ndoes_num    
    

RETURN
END SUBROUTINE D3_Get_Crack_Center
