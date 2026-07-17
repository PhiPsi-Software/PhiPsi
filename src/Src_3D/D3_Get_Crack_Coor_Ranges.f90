!-----------------------------------------------------------
! Brief: Compute axis-aligned bounding range of one 3D crack mesh.
!
! Parameters:
!   Input:  i_C         - crack index
!   Output: Coor_Ranges - min/max for x, y, z (3 x 2)
!
! Notes:   Uses minval/maxval over Crack3D_Meshed_Node(i_C)%row.
!-----------------------------------------------------------

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
