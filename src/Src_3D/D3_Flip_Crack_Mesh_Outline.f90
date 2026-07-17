!-----------------------------------------------------------
! Brief: Reverse the node order of a crack surface outline.
!
! Parameters:
!   Input:  i_C - crack index
!
! Notes:   Reverses Crack3D_Meshed_Outline(i_C) so that endpoint nodes
!          of every boundary line are swapped; reverses the implicit
!          traversal direction of the outline.
!-----------------------------------------------------------

subroutine D3_Flip_Crack_Mesh_Outline(i_C)
!     Outline of the ordering and other information of discrete nodes at the front edge vertices of flipped 3D crack surfaces.
!     2022-07-13. Ref: My PhiPsi Development Notebook V1-P41.

!     ----------------------------
!     Read public variable module
!     ----------------------------
use Global_Float_Type
use Global_Common
use Global_Crack_3D

!     --------------------------
!     Variable Type Declaration
!     --------------------------
implicit none
integer,intent(in)::i_C
integer :: i_out,Outline_num,c_Max_N_Node_3D
integer,allocatable::Temp(:,:)

c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)

!     ---------------------
!     Main program section
!     ---------------------
allocate(Temp(c_Max_N_Node_3D,4))
Temp(1:c_Max_N_Node_3D,1:4)= Crack3D_Meshed_Outline(i_C)%row(1:c_Max_N_Node_3D,1:4)

Outline_num = Crack3D_Meshed_Outline_num(i_C)

do i_out =1,Outline_num
Crack3D_Meshed_Outline(i_C)%row(i_out,1)= Temp(Outline_num-i_out+1,2)
Crack3D_Meshed_Outline(i_C)%row(i_out,2)= Temp(Outline_num-i_out+1,1)
Crack3D_Meshed_Outline(i_C)%row(i_out,3)= Temp(Outline_num-i_out+1,3)
enddo

deallocate(Temp)

return
end subroutine D3_Flip_Crack_Mesh_Outline
