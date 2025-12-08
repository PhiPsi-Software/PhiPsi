 
subroutine Cal_Crack_Avarage_Radius(i_Crack,Crack_Radius)

use Global_Float_Type
use Global_Crack_Common
use Global_Crack_3D
use Global_Model
use Global_Elem_Area_Vol
use Global_Common
use Global_HF
use Global_Post
use Global_Filename
use omp_lib
use Global_Crack
use Global_Cal_Ele_Num_by_Coors_3D  

implicit none
integer,intent(in)::i_Crack
real(kind=FT),intent(out)::Crack_Radius
real(kind=FT) center_x,center_y,center_z,c_x_front,c_y_front,c_z_front
integer Num_CrMesh_Outlines,i_Out_Node,c_Mesh_Node
real(kind=FT) tem_dis,dis_Sum
real(kind=FT) Tool_Function_2Point_Dis_3D

Crack_Radius = ZR
dis_Sum = ZR

center_x = Crack3D_Cir_Coor(i_Crack,1) 
center_y = Crack3D_Cir_Coor(i_Crack,2) 
center_z = Crack3D_Cir_Coor(i_Crack,3) 
Num_CrMesh_Outlines = Crack3D_Meshed_Outline_num(i_Crack)
do i_Out_Node = 1,Num_CrMesh_Outlines
  c_Mesh_Node = Crack3D_Meshed_Outline(i_Crack)%row(i_Out_Node,1)
  c_x_front  = Crack3D_Meshed_Node(i_Crack)%row(c_Mesh_Node,1) 
  c_y_front   = Crack3D_Meshed_Node(i_Crack)%row(c_Mesh_Node,2) 
  c_z_front   = Crack3D_Meshed_Node(i_Crack)%row(c_Mesh_Node,3)
  tem_dis = Tool_Function_2Point_Dis_3D([center_x,center_y,center_z],[c_x_front,c_y_front,c_z_front])
  dis_Sum = dis_Sum + tem_dis
enddo

Crack_Radius = dis_Sum/Num_CrMesh_Outlines
   
return 
end SUBROUTINE Cal_Crack_Avarage_Radius             
