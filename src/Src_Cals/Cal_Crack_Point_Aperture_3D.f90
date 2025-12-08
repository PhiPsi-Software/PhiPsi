 
subroutine Cal_Crack_Point_Aperture_3D(c_DISP,Crack_Number,Crack_Point,Crack_Point_Aperture,&
                                       Fluid_element_number,Fluid_node_number,&
                                       Crack_mesh_ele_number,Crack_mesh_node_number)

use Global_Float_Type
use Global_Common
use Global_Crack_Common
use Global_Crack_3D
use Global_Model
use Global_Elem_Area_Vol
use Global_DISP
use Global_Inter_Tool_Cal_Dis_Point_to_3D_Quad
use Global_INTERFACE_Tool_ThetaX_ThetaY_ThetaZ_3D_rotation

implicit none

real(kind=FT),intent(in)::c_DISP(Total_FD)
integer,intent(in)::Crack_Number
real(kind=FT),intent(in)::Crack_Point(3)
real(kind=FT),intent(out)::Crack_Point_Aperture(3)
integer Fluid_element_number,Fluid_node_number,Crack_mesh_ele_number,Crack_mesh_node_number

integer Solid_Element
real(kind=FT) c_mesh_element_points(3,3)
real(kind=FT) Nor_Vector(3)
integer c_mesh_node(3)
real(kind=FT) c_Kesi,c_Yita,c_Zeta
integer i_Crack_Ele
logical Yes_on
real(kind=FT) F_funtion_r
real(kind=FT) Point_Coor_Local(3)
integer Ele_Num_Cache
real(kind=FT) tem_vector(3)
real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8)
real(kind=FT) N(3,24)
real(kind=FT) c_N(8)
integer i_N,NODES_iE(8),Num_Enri_Node
integer i_F
real(kind=FT) c_T_Matrx(3,3)
integer Num_CrMesh_Outlines,i_Out_Node
integer c_Outline_Node
real(kind=FT) Nearest_Point(3),Point_Pre_with_Gap(3),Point_Next_with_Gap(3)  
real(kind=FT),ALLOCATABLE::c_Polygon(:,:)
real(kind=FT) Vector_Crack_x(3),Vector_Crack_y(3),Vector_Crack_z(3)
real(kind=FT) ThetaX,ThetaY,ThetaZ
real(kind=FT) Tip_vector_x(3),Tip_vector_y(3),Tip_vector_z(3)
real(kind=FT) Tool_Function_2Point_Dis_3D
logical yes_found_Nor_Vector
integer P1_CalP_num,P2_CalP_num,P3_CalP_num
integer i_FluidEle
real(kind=FT) c_P1(3),c_P2(3),c_P3(3)
real(kind=FT) Yes_Point_on_3D_Triangle_Tol


Crack_Point_Aperture(1:3) = ZR
Yes_Point_on_3D_Triangle_Tol = 1.0D-3

Num_CrMesh_Outlines = Crack3D_Meshed_Outline_num(Crack_Number)
allocate(c_Polygon(Num_CrMesh_Outlines,3))
do i_Out_Node = 1,Num_CrMesh_Outlines
  c_Outline_Node = Crack3D_Meshed_Outline(Crack_Number)%row(i_Out_Node,1)
  c_Polygon(i_Out_Node,1:3) = Crack3D_Meshed_Node(Crack_Number)%row(c_Outline_Node,1:3) 
enddo
call Tool_Get_Nearest_Point_on_3D_Polygon&
          (Num_CrMesh_Outlines,c_Polygon(1:Num_CrMesh_Outlines,1:3),&
           Crack_Point,Nearest_Point,Point_Pre_with_Gap,Point_Next_with_Gap) 
deallocate(c_Polygon)           
if(Tool_Function_2Point_Dis_3D(Nearest_Point,Crack_Point) <Tol_11) then
    return
endif        
Vector_Crack_x(1:3) = [ONE,ZR,ZR]
Vector_Crack_y(1:3) = [ZR,ONE,ZR]
Vector_Crack_z(1:3) = [ZR,ZR,ONE]
Tip_vector_x = Nearest_Point - Crack_Point
call Vector_Normalize(3,Tip_vector_x)   
Tip_vector_z = -(Point_Next_with_Gap - Point_Pre_with_Gap)
call Vector_Normalize(3,Tip_vector_z) 
call Vector_Cross_Product_3(Tip_vector_x,Tip_vector_z,Tip_vector_y)  
Tip_vector_y = -Tip_vector_y
call Vector_Normalize(3,Tip_vector_y) 
call Tool_ThetaX_ThetaY_ThetaZ_3D_rotation(1,                       &
          Tip_vector_x,&
          Tip_vector_y,&
          Tip_vector_z,&
          Vector_Crack_x(1:3),Vector_Crack_y(1:3),Vector_Crack_z(1:3),&
          ThetaX,ThetaY,ThetaZ,c_T_Matrx(1:3,1:3))       
    
Ele_Num_Cache = 1
call Cal_Ele_Num_by_Coors_3D(Crack_Point(1),Crack_Point(2),Crack_Point(3),Ele_Num_Cache,Solid_Element) 
if(Solid_Element==0)then
  if (Key_Stop_Outside_Crack ==1) then
      print *,'       WARN-2023081201 :: Solid_Element=0 in Cal_Crack_Point_Aperture_3D.f90!'
      print *,'                          The coordinates of point:',Crack_Point(1:3)   
      print *,'                          Crack ',Crack_Number, ' will stop propagation!'
      Crack_Type_Status_3D(Crack_Number,3) = 0
      return
  endif
  if (Key_Allow_3D_Outside_Crack==1) then
      return
  endif  
  print *,'       ERROR-2023081202 :: Solid_Element=0 in Cal_Crack_Point_Aperture_3D.f90!'
  call Warning_Message('S',Keywords_Blank)    
endif 

NODES_iE = G_NN(1:8,Solid_Element)

c_X_NODES = G_X_NODES(1:8,Solid_Element)
c_Y_NODES = G_Y_NODES(1:8,Solid_Element)
c_Z_NODES = G_Z_NODES(1:8,Solid_Element)

call Cal_KesiYita_by_Coor_3D(Crack_Point,Solid_Element,c_Kesi,c_Yita,c_Zeta)


call Cal_N_3D(c_Kesi,c_Yita,c_Zeta,N)
c_N(1:8) = N(1,1:24:3)

yes_found_Nor_Vector = .False.
if(Fluid_element_number>0)then
    Nor_Vector(1:3) = Cracks_FluidEle_Vector_3D(Crack_Number)%row(Fluid_element_number,1:3)    
    yes_found_Nor_Vector = .True.
    goto 100
endif
if(Crack_mesh_ele_number>0)then
    Nor_Vector(1:3) = Crack3D_Meshed_Ele_Nor_Vector(Crack_Number)%row(Crack_mesh_ele_number,1:3)
    yes_found_Nor_Vector = .True.
    goto 100
endif
if(Fluid_node_number>0)then
    Nor_Vector(1:3) = Cracks_CalP_Orient_3D(Crack_Number)%row(Fluid_node_number,1:3)
    yes_found_Nor_Vector = .True.
    goto 100
endif
if(Crack_mesh_node_number>0)then
    Nor_Vector(1:3) = Crack3D_Meshed_Node_Nor_Vector(Crack_Number)%row(Crack_mesh_node_number,1:3)
    yes_found_Nor_Vector = .True.
    goto 100
endif
if(Crack_mesh_ele_number==0) then
    do i_Crack_Ele = 1,Crack3D_Meshed_Ele_num(Crack_Number) 
        c_mesh_node(1:3) = Crack3D_Meshed_Ele(Crack_Number)%row(i_Crack_Ele,1:3) 
        c_mesh_element_points(1,1:3) = Crack3D_Meshed_Node(Crack_Number)%row(c_mesh_node(1),1:3)
        c_mesh_element_points(2,1:3) = Crack3D_Meshed_Node(Crack_Number)%row(c_mesh_node(2),1:3)
        c_mesh_element_points(3,1:3) = Crack3D_Meshed_Node(Crack_Number)%row(c_mesh_node(3),1:3)
        call Tool_Yes_Point_on_3D_Triangle_with_Tol(Crack_Point,&
                c_mesh_element_points(1,1:3),c_mesh_element_points(2,1:3),c_mesh_element_points(3,1:3),&
                Yes_on,Yes_Point_on_3D_Triangle_Tol)
        if(Yes_on) then
            Nor_Vector(1:3) = Crack3D_Meshed_Ele_Nor_Vector(Crack_Number)%row(i_Crack_Ele,1:3)
            yes_found_Nor_Vector = .True.
            goto 100
        endif
    enddo
endif
if(Fluid_element_number==0) then
    do i_FluidEle = 1,Cracks_FluidEle_num_3D(Crack_Number)
        P1_CalP_num = Cracks_FluidEle_CalP_3D(Crack_Number)%row(i_FluidEle,1)     
        P2_CalP_num = Cracks_FluidEle_CalP_3D(Crack_Number)%row(i_FluidEle,2)  
        P3_CalP_num = Cracks_FluidEle_CalP_3D(Crack_Number)%row(i_FluidEle,3)
        c_P1(1:3)   = Cracks_CalP_Coors_3D(Crack_Number)%row(P1_CalP_num,1:3)
        c_P2(1:3)   = Cracks_CalP_Coors_3D(Crack_Number)%row(P2_CalP_num,1:3)
        c_P3(1:3)   = Cracks_CalP_Coors_3D(Crack_Number)%row(P3_CalP_num,1:3)
        
        call Tool_Yes_Point_on_3D_Triangle_with_Tol(Crack_Point,c_P1,c_P2,c_P3,Yes_on,Yes_Point_on_3D_Triangle_Tol)
        if(Yes_on) then
            Nor_Vector(1:3) = Cracks_FluidEle_Vector_3D(Crack_Number)%row(i_FluidEle,1:3)   
            yes_found_Nor_Vector = .True.
            goto 100
        endif
    enddo
endif

100 continue 

if(yes_found_Nor_Vector .eqv. .False.) then
  print *,"    ERROR-2023081203 :: point is not on the crack surface, in Cal_Crack_Point_Aperture_3D.f90!"
  print *,'                        Crack_Number:',Crack_Number  
  print *,'                        Crack_Point:',Crack_Point(1:3)   
  call Warning_Message('S',Keywords_Blank)  
endif

if(sum(c_POS_3D(NODES_iE,Crack_Number))==0) then
  print *,"     ERROR-2023081201 :: in Cal_Crack_Point_Aperture_3D.f90."
  print *,"                         sum(c_POS_3D(NODES_iE,Crack_Number))==0!"
  call Warning_Message('S',Keywords_Blank)  
endif

do i_N = 1,8
  if (Enriched_Node_Type_3D(NODES_iE(i_N),Crack_Number) ==2) then
      Num_Enri_Node = c_POS_3D(NODES_iE(i_N),Crack_Number)
      if ((Elem_Type_3D(Solid_Element,Crack_Number) == 2)) then
          Crack_Point_Aperture(1) = Crack_Point_Aperture(1) + c_N(i_N)*c_DISP(Num_Enri_Node*3-2)
          Crack_Point_Aperture(2) = Crack_Point_Aperture(2) + c_N(i_N)*c_DISP(Num_Enri_Node*3-1)
          Crack_Point_Aperture(3) = Crack_Point_Aperture(3) + c_N(i_N)*c_DISP(Num_Enri_Node*3)
      end if
  elseif (Enriched_Node_Type_3D(NODES_iE(i_N),Crack_Number).eq.3)then

  elseif (Enriched_Node_Type_3D(NODES_iE(i_N),Crack_Number).eq.1)then
      
      Point_Coor_Local = MATMUL(c_T_Matrx,Crack_Point(1:3) - Nearest_Point(1:3))
      F_funtion_r = sqrt(Point_Coor_Local(1)**2 + Point_Coor_Local(2)**2)
      
      
      do i_F =1,1
          Num_Enri_Node=c_POS_3D(NODES_iE(i_N),Crack_Number)+i_F-1 
          Crack_Point_Aperture(1) = Crack_Point_Aperture(1) + sqrt(F_funtion_r)*c_N(i_N)*c_DISP(Num_Enri_Node*3-2)
          Crack_Point_Aperture(2) = Crack_Point_Aperture(2) + sqrt(F_funtion_r)*c_N(i_N)*c_DISP(Num_Enri_Node*3-1)
          Crack_Point_Aperture(3) = Crack_Point_Aperture(3) + sqrt(F_funtion_r)*c_N(i_N)*c_DISP(Num_Enri_Node*3)
      end do    
      
  end if
end do

tem_vector(1:3) = Crack_Point_Aperture(1:3)
Crack_Point_Aperture(1) = TWO*(Nor_Vector(1)*tem_vector(1))
Crack_Point_Aperture(2) = TWO*(Nor_Vector(2)*tem_vector(2))
Crack_Point_Aperture(3) = TWO*(Nor_Vector(3)*tem_vector(3))


return 
end SUBROUTINE Cal_Crack_Point_Aperture_3D      
