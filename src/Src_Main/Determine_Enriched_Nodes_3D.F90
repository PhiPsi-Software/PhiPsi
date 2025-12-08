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
 
SUBROUTINE Determine_Enriched_Nodes_3D(ifra,isub)
! This subroutine is used to determine the enriched nodes (3D problem).
      
!..................................
! Explanation of Variable Meanings
!..................................
! Dis_Node_to_FS(num_Node,num_Crack) ! Symbolic distance from nodes to each crack surface, FS
! stands for Fracture Surface
! Elem_Type_3D(Num_Elem,num_crack);  ! Element type: 1--Tip element;
!                                         !               2--Heaviside enriched element;
!                                         !               3--Fully cracked element with kink point
!                                         !               4--Junction elements
! 
! 
      
!.............................
! Read public variable module
!.............................
use Global_Float_Type
use Global_Common   
use Global_Filename
use Global_Model
use Global_XFEM_Elements
use Global_Elem_Area_Vol
use Global_Crack_Common
use Global_Crack_3D
use Global_HF
use Global_Inter_Tool_Cal_Dis_Point_to_3D_Quad
use omp_lib
use Global_Ragged_Array_Real_Classs
use Global_Cal_Ele_Num_by_Coors_3D  
use Global_INTERFACE_Vector_Normalize
use Global_INTERFACE_Tool_ThetaX_ThetaY_ThetaZ_3D_rotation
use Global_INTERFACE_Vector_Location_Int_v2
use Global_INTERFACE_D3_Get_Signed_Dis_to_Crack_Mesh
use Global_INTERFACE_D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane
!
!...........................
! Variable Type Declaration
!...........................
implicit none
integer,intent(in)::ifra,isub
integer i_C,i_Node,i_E 
real(kind=FT) c_Node_Coor(3)
integer c_NN(8)
real(kind=FT) c_max,c_min
real(kind=FT) Vol_Enrich
integer iEl_contour
real(kind=FT) Tri123_Distance,Tri123_PER(3)
logical Tri123_Yes_PER_in,Tri123_Yes_PER_on
integer i_Crack_Ele
integer Crack_Node1,Crack_Node2,Crack_Node3
real(kind=FT) Point1(3),Point2(3),Point3(3)
real(kind=FT) c_Min_Signed_Dis
integer num_Vertex,i_V,c_Mesh_Node
integer in_Ele_Num,mat_num
real(kind=FT) c_x,c_y,c_z
integer c_Mesh_Nex_Node,c_Mesh_Pre_Node
real(kind=FT) c_InterSection_P(3)
logical c_Yes_Inter
real(kind=FT) c_Point_1(3),c_Point_2(3)
real(kind=FT) Vector_Crack_x(3),Vector_Crack_y(3),Vector_Crack_z(3)
real(kind=FT) ThetaX,ThetaY,ThetaZ,T_Matrx(3,3)
integer c_V,c_V_Pre,c_V_Nex
real(kind=FT) Vertex_1(3),Vertex_2(3),Vertex_1_Pre(3),Vertex_2_Nex(3)
real(kind=FT) Vertex_3(3),Vertex_3_Nex(3)
real(kind=FT) c_Inter_Point_A(3),c_Inter_Point_B(3)
integer c_Vertex_1,c_Vertex_2
integer num_tip_enrich,ref_elem
logical added_tip_el(Num_Elem)
logical added_tip_node(Num_Elem,8)
integer added_tip_node_refe(Num_Elem,8)
integer num_Inter_Outline
integer num_H_enrich
integer XFEM_counter,FEM_counter
real(kind=FT) c_Four_Points(4,3)
integer j_C
logical Crack_Coor_Overlap_Status(num_Crack,num_Crack)
integer i_Cr_Node1,i_Cr_Node2,i_Cr_Node3
integer j_Cr_Node1,j_Cr_Node2,j_Cr_Node3
integer j_Crack_Ele
real(kind=FT) c_Tri_1(3,3),c_Tri_2(3,3)
real(kind=FT) c_Inter_Points(10,3)
integer c_Elem
integer Tem_Elem_Type(Num_Elem,num_Crack)
                                            ! (Used for Heaviside and Junction element type determination when crack tip enrichment is not
                                            ! considered)
real(kind=FT) Cracks_Coor_Centers(num_Crack,3)
real(kind=FT) new_A(3),delta_L
logical Yes_Found_Min_Signed_Dis
real(kind=FT) Check_Ball_R,c_n_Vector(3)
type(Ragged_Array_2D),allocatable::Ragged_Array_2D_n_Vectors(:)
real(kind=FT) c_Dis_Node_to_FS
real(kind=FT) c_Dis_Node_to_FS_v2
logical c_Yes_Side_on
integer c_Side_Tri_Number
integer c_Num_Eles,c_Ele
real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8)
integer i_G  
real(kind=FT),ALLOCATABLE::kesi_Enr(:),yita_Enr(:),zeta_Enr(:),weight_Enr(:)     
real(kind=FT) N(3,24),Global_coor_Gauss(3)
real(kind=FT) c_PER_Node_to_FS(3)
real(kind=FT) ,ALLOCATABLE::c_Distance(:,:)
logical c_Yes_Node_PER_in_FS
integer Nodes_to_be_Checked(Num_Node),c_Count_Nodes,c_Node
logical c_Logical_Positive,c_Logical_Negative
real(kind=FT) c_Sign
integer Num_Gauss_Cubes,Tem_Num_Gauss_3D
type(Ragged_Array_1D),allocatable::RaggedArray1D_Dis_Node_to_FS_v2(:)
real(kind=FT) c_Signed_Dis_v2
logical c_Yes_Parallel
integer c_Cr_Location
integer n_Vector_Sign_Int_1(8),n_Vector_Sign_Int_2(8)
integer Tool_Function_Sign_Int
integer Heaviside_Enrich_Type_Int_1,Heaviside_Enrich_Type_Int_2
integer c_Elem_num_Related_Cracks
integer c_Tri_Edge_Num
integer Check_Mesh_Node_1,Check_Mesh_Node_2
logical Mesh_Node_1_on_Outline,Mesh_Node_2_on_Outline
integer c_num_Outlines,c_location
integer,ALLOCATABLE::tem_Solid_El_Vertex_num(:,:,:)
integer Vertex_1_num,Vertex_2_num,Vertex_3_num
real(kind=FT) c_Eight_Values(8)
integer Ele_Num_Cache
integer c_Logical_Positive_Num,c_Logical_Negative_Num
integer c_count
integer c_num_Points
logical c_Yes_Overlaped_Ele
logical c_Logical_Dif
real(kind=FT) Vol_Enrich_Eles(Num_Elem)
integer Vol_Enrich_Eles_Flag(Num_Elem)
integer num_XFEM_Elem_Flag(Num_Elem)
901  FORMAT(5X,'Number of elements with enriched nodes is ',I8)         
902  FORMAT(5X,'Number of elements without enriched nodes is ',I8)           
1001 FORMAT(5X,'Enrichment ratio is ',F8.4,'%')     
1901 FORMAT(5X,'Caution :: Duplicate Heaviside element ',I8,' for cracks ',I4, ' and ',I4 ,' found!')

#ifndef Silverfrost       
!......................................................
! Allocate memory space. 2022-07-20. IMPROV2022072001.
!......................................................
IF(ALLOCATED(Dis_Node_to_FS)) DEALLOCATE(Dis_Node_to_FS)  
ALLOCATE(Dis_Node_to_FS(Num_Node))
IF(ALLOCATED(Elem_Type_3D)) DEALLOCATE(Elem_Type_3D)
ALLOCATE(Elem_Type_3D(Num_Elem,num_Crack))
IF(ALLOCATED(Enriched_Node_Type_3D)) DEALLOCATE(Enriched_Node_Type_3D)   
ALLOCATE(Enriched_Node_Type_3D(Num_Node,num_Crack))
IF(ALLOCATED(Ele_Num_Tip_Enriched_Node_3D)) DEALLOCATE(Ele_Num_Tip_Enriched_Node_3D)
ALLOCATE(Ele_Num_Tip_Enriched_Node_3D(Num_Node)) 
IF(ALLOCATED(Enriched_Node_Crack_n_Vector_3D)) DEALLOCATE(Enriched_Node_Crack_n_Vector_3D)      
ALLOCATE(Enriched_Node_Crack_n_Vector_3D(Num_Node))
! -------------High memory usage (fixed using RaggedArray1D_Dis_Node_to_FS_v2)-------------
ALLOCATE(RaggedArray1D_Dis_Node_to_FS_v2(Num_Node))
! -------------Low memory usage-------------
IF(ALLOCATED(Solid_El_num_Crs)) DEALLOCATE(Solid_El_num_Crs)
ALLOCATE(Solid_El_num_Crs(Num_Elem))
IF(ALLOCATED(Solid_El_Crs)) DEALLOCATE(Solid_El_Crs) 
ALLOCATE(Solid_El_Crs(Num_Elem,Solid_El_Max_num_Crs))
! ----High memory usage (resolved with staggered arrays, IMPROV2022091701.)-----
if (Key_Junction_Enrich==1) then
    IF(ALLOCATED(Jun_Ele_Negative_Cr_Num_3D)) DEALLOCATE(Jun_Ele_Negative_Cr_Num_3D)    
    ALLOCATE(Jun_Ele_Negative_Cr_Num_3D(Num_Elem))
    IF(ALLOCATED(Node_Jun_elem_3D)) DEALLOCATE(Node_Jun_elem_3D)    
    ALLOCATE(Node_Jun_elem_3D(Num_Node))
    IF(ALLOCATED(Coors_Junction_3D)) DEALLOCATE(Coors_Junction_3D)    
    ALLOCATE(Coors_Junction_3D(Num_Elem))
    
endif
! ---Extremely low memory usage---
IF(ALLOCATED(Vector_o_Orient)) DEALLOCATE(Vector_o_Orient)
ALLOCATE(Vector_o_Orient(num_Crack,3))
IF(ALLOCATED(Sign_o_Orient)) DEALLOCATE(Sign_o_Orient)
ALLOCATE(Sign_o_Orient(num_Crack))     
! ---Very low memory usage (very low memory usage)---
IF(ALLOCATED(Elem_num_Related_Cracks)) DEALLOCATE(Elem_num_Related_Cracks)    
ALLOCATE(Elem_num_Related_Cracks(Num_Elem))
! ---Very low memory usage (very low memory usage)---
IF(ALLOCATED(Elem_Related_Cracks)) DEALLOCATE(Elem_Related_Cracks)    
ALLOCATE(Elem_Related_Cracks(Num_Elem,Key_Ele_Max_Related_Cracks))
! ---Low Memory Usage---
IF(ALLOCATED(tem_Solid_El_Vertex_num)) DEALLOCATE(tem_Solid_El_Vertex_num)   
ALLOCATE(tem_Solid_El_Vertex_num(Num_Elem,Solid_El_Max_num_Crs,20))
!.................................
! Related variable initialization
!.................................
Elem_Type_3D(1:Num_Elem,1:num_Crack)          = 0
Tem_Elem_Type(1:Num_Elem,1:num_Crack)         = 0
Enriched_Node_Type_3D(1:Num_Node,1:num_Crack) = 0  
Vector_o_Orient(1:num_Crack,1:3)           = ZR
Sign_o_Orient(1:num_Crack)                 = 0
allocate(Ragged_Array_2D_n_Vectors(Num_Node))
Elem_num_Related_Cracks(1:Num_Elem)            = 0
Elem_Related_Cracks(1:Num_Elem,1:Key_Ele_Max_Related_Cracks)  = 0
!----------  IMPROV2022081502.
Solid_El_num_Crs(1:Num_Elem)  = 0
Solid_El_Crs(1:Num_Elem,1:Solid_El_Max_num_Crs) = 0


tem_Solid_El_Vertex_num(1:Num_Elem,1:Solid_El_Max_num_Crs,1:20)   = 0

!...........................................
!                                         .
! Step 1: Data Preparation.               .
!                                         .
!...........................................
!######################################################################################
! Step 0.1: Obtain the coordinate range of each discrete fracture surface and store it
!          Crack_Coor_Range(num_Crack,3,2)
!          NEWFTU2022050401.
!######################################################################################
! call D3_Get_Cracks_Coor_Ranges(Crack_Coor_Range)   !Already parallelized with OpenMP.
call D3_Get_Cracks_Coor_Ranges


!###############################################################################
! Step 0.2: Determine whether the coordinate ranges between each crack overlap.
! Store into Crack_Coor_Overlap_Status(1:num_Crack,1:num_Crack)
!          NEWFTU2022050402.
!###############################################################################
! call D3_Get_Cracks_Coor_Overlap_Status(Crack_Coor_Range, Crack_Coor_Overlap_Status) !Already
! parallelized with OpenMP.
call D3_Get_Cracks_Coor_Overlap_Status(Crack_Coor_Overlap_Status)

!##########################################################################################
! Step 0.3: Obtain the center positions of each discrete fracture plane. NEWFTU2022050501.
!##########################################################################################
call D3_Get_Cracks_Centers(Cracks_Coor_Centers)

!##################################################################################################
! Step 0.4: Generate the variable object of the crack enhancement element related disparity array.
! 2022-09-04. IMPROV2022090401.
!##################################################################################################
IF(ALLOCATED(Solid_El_Vertex_Num)) DEALLOCATE(Solid_El_Vertex_Num)    
IF(ALLOCATED(Solid_El_Vertex_Coor)) DEALLOCATE(Solid_El_Vertex_Coor)    
IF(ALLOCATED(Solid_El_Vertex_Nor_Vec)) DEALLOCATE(Solid_El_Vertex_Nor_Vec)    
IF(ALLOCATED(Solid_El_Vertex_x_Vec)) DEALLOCATE(Solid_El_Vertex_x_Vec)    
IF(ALLOCATED(Solid_El_Vertex_y_Vec)) DEALLOCATE(Solid_El_Vertex_y_Vec)    
IF(ALLOCATED(Solid_El_Vertex_z_Vec)) DEALLOCATE(Solid_El_Vertex_z_Vec)    
IF(ALLOCATED(Solid_El_Pre_Vertex_Coor)) DEALLOCATE(Solid_El_Pre_Vertex_Coor)    
IF(ALLOCATED(Solid_El_Pre_Vertex_Nor_Vec)) DEALLOCATE(Solid_El_Pre_Vertex_Nor_Vec)    
IF(ALLOCATED(Solid_El_Pre_Vertex_x_Vec)) DEALLOCATE(Solid_El_Pre_Vertex_x_Vec)    
IF(ALLOCATED(Solid_El_Pre_Vertex_y_Vec)) DEALLOCATE(Solid_El_Pre_Vertex_y_Vec)    
IF(ALLOCATED(Solid_El_Pre_Vertex_z_Vec)) DEALLOCATE(Solid_El_Pre_Vertex_z_Vec)    
IF(ALLOCATED(Solid_El_Nex_Vertex_Coor)) DEALLOCATE(Solid_El_Nex_Vertex_Coor)    
IF(ALLOCATED(Solid_El_Nex_Vertex_Nor_Vec)) DEALLOCATE(Solid_El_Nex_Vertex_Nor_Vec)  
IF(ALLOCATED(Solid_El_Nex_Vertex_x_Vec)) DEALLOCATE(Solid_El_Nex_Vertex_x_Vec)    
IF(ALLOCATED(Solid_El_Nex_Vertex_y_Vec)) DEALLOCATE(Solid_El_Nex_Vertex_y_Vec)    
IF(ALLOCATED(Solid_El_Nex_Vertex_z_Vec)) DEALLOCATE(Solid_El_Nex_Vertex_z_Vec)    
IF(ALLOCATED(Solid_El_Tip_BaseLine)) DEALLOCATE(Solid_El_Tip_BaseLine)    
IF(ALLOCATED(Solid_El_Tip_BaseLine_Nor_Vec)) DEALLOCATE(Solid_El_Tip_BaseLine_Nor_Vec)    
IF(ALLOCATED(Solid_El_Tip_BaseLine_x_Vec)) DEALLOCATE(Solid_El_Tip_BaseLine_x_Vec)    
IF(ALLOCATED(Solid_El_Tip_BaseLine_y_Vec)) DEALLOCATE(Solid_El_Tip_BaseLine_y_Vec)   
IF(ALLOCATED(Solid_El_Tip_BaseLine_z_Vec)) DEALLOCATE(Solid_El_Tip_BaseLine_z_Vec) 
IF(ALLOCATED(Solid_El_Tip_BaseLine_T_theta)) DEALLOCATE(Solid_El_Tip_BaseLine_T_theta)   
IF(ALLOCATED(Solid_El_Tip_BaseLine_T_Matrix)) DEALLOCATE(Solid_El_Tip_BaseLine_T_Matrix) 
! Does not use memory at this time:
allocate(Solid_El_Vertex_Num(num_Elem))
allocate(Solid_El_Vertex_Coor(num_Elem))
allocate(Solid_El_Vertex_Nor_Vec(num_Elem))
allocate(Solid_El_Vertex_x_Vec(num_Elem))
allocate(Solid_El_Vertex_y_Vec(num_Elem))
allocate(Solid_El_Vertex_z_Vec(num_Elem))
allocate(Solid_El_Pre_Vertex_Coor(num_Elem))
allocate(Solid_El_Pre_Vertex_Nor_Vec(num_Elem))
allocate(Solid_El_Pre_Vertex_x_Vec(num_Elem))
allocate(Solid_El_Pre_Vertex_y_Vec(num_Elem))
allocate(Solid_El_Pre_Vertex_z_Vec(num_Elem))
allocate(Solid_El_Nex_Vertex_Coor(num_Elem))
allocate(Solid_El_Nex_Vertex_Nor_Vec(num_Elem))
allocate(Solid_El_Nex_Vertex_x_Vec(num_Elem))
allocate(Solid_El_Nex_Vertex_y_Vec(num_Elem))
allocate(Solid_El_Nex_Vertex_z_Vec(num_Elem))
allocate(Solid_El_Tip_BaseLine(num_Elem))
allocate(Solid_El_Tip_BaseLine_Nor_Vec(num_Elem))
allocate(Solid_El_Tip_BaseLine_x_Vec(num_Elem))
allocate(Solid_El_Tip_BaseLine_y_Vec(num_Elem))
allocate(Solid_El_Tip_BaseLine_z_Vec(num_Elem))
allocate(Solid_El_Tip_BaseLine_T_theta(num_Elem))
allocate(Solid_El_Tip_BaseLine_T_Matrix(num_Elem))
  
!##################################################################################
! Calculate the relevant information of the vertices 
! on the boundaries of discrete fracture surfaces contained in each solid element.
! Including the adjacent previous and next vertices (up to 3 vertices)
!##################################################################################
! There may be data conflicts, so OpenMP is not advisable. 2022-08-19. In addition, the computation
! in this part is not large. 2022-10-06.
do i_C=1,num_Crack
  num_Vertex = Crack3D_Meshed_Outline_num(i_C)
  do i_V =1,num_Vertex
      ! Discrete node number at the crack surface boundary
      c_Mesh_Node = Crack3D_Meshed_Outline(i_C)%row(i_V,1)   
      ! A discrete node number on the edge of the crack surface
      if(i_V==1) then
        c_Mesh_Pre_Node=Crack3D_Meshed_Outline(i_C)%row(num_Vertex,1)
        c_V_Pre = num_Vertex
      else
        c_Mesh_Pre_Node=Crack3D_Meshed_Outline(i_C)%row(i_V-1,1)
        c_V_Pre = i_V-1
      endif
      ! A discrete node number beneath the crack surface boundary
      if(i_V==num_Vertex) then
        c_Mesh_Nex_Node=Crack3D_Meshed_Outline(i_C)%row(1,1)
        c_V_Nex = 1
      else
        c_Mesh_Nex_Node=Crack3D_Meshed_Outline(i_C)%row(i_V,2)
        c_V_Nex = i_V+1
      endif              
      in_Ele_Num  = Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(c_Mesh_Node)
      
      ! If the discrete nodes on the crack surface are outside the model, proceed to the next iteration,
      ! 2020-03-09
      if(in_Ele_Num==0)then
          cycle
      else
          ! Allocate memory for variables related to crack tip enhancement corresponding to this element and
          ! initialize them. 2022-09-04. IMPROV2022090401.
          call D3_Allocate_Ele_Memory(in_Ele_Num,1)   
      endif
      
      !****************************************
      ! Reduce memory usage. IMPROV2022081502.
      !****************************************
      ! Check whether this element has already been marked for the crack. 2022-08-15. IMPROV2022081502.
      call Vector_Location_Int_v2(Solid_El_Max_num_Crs,Solid_El_Crs(in_Ele_Num,1:Solid_El_Max_num_Crs),i_C,c_Cr_Location)
          
      ! If this crack was not previously marked in the element, mark it. 2022-08-15.
      if(c_Cr_Location==0 .or. Solid_El_num_Crs(in_Ele_Num)==0)then 
          ! Update the number of cracks associated with each solid element. 2022-08-15. IMPROV2022081502.
          Solid_El_num_Crs(in_Ele_Num) = Solid_El_num_Crs(in_Ele_Num) + 1
          if(Solid_El_num_Crs(in_Ele_Num) > Solid_El_Max_num_Crs) then
              print *,'    ERROR-2022081501 :: Solid_El_num_Crs(in_Ele_Num)>Solid_El_Max_num_Crs!'
              print *,'                        In Determine_Enriched_Nodes_3D.f90!'
              print *,'                        Try to increase Solid_El_Max_num_Crs!'
              call Warning_Message('S',Keywords_Blank)  
          endif
          
          ! Update the crack number related to each solid element. 2022-08-15. IMPROV2022081502.
          Solid_El_Crs(in_Ele_Num,Solid_El_num_Crs(in_Ele_Num)) = i_C
          
          ! Number of discrete fracture surface boundary vertices contained in each solid element
          Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num)) = &
                          Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num))  +1  
          !2022-12-22.                
          if(Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num)) > 20) then
              print *,'    ERROR-2022112301 :: Solid_El_Vertex_Num(in_Ele_Num)>20!'
              print *,'                        In Determine_Enriched_Nodes_3D.f90!'
              print *,'                        Increase size 3 of tem_Solid_El_Vertex_num!'
              call Warning_Message('S',Keywords_Blank)  
          endif             
          
          
          ! The vertex numbers of cracks in each solid element i_C. 2022-08-26.
          tem_Solid_El_Vertex_num(in_Ele_Num,Solid_El_num_Crs(in_Ele_Num),&
                    Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num))) =i_V
      ! If this element has previously marked this crack, increase the number of vertices. 2022-08-15.
      elseif(c_Cr_Location>=1)then
          ! Number of discrete fracture surface boundary vertices contained in each solid element
          Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num)) = &
                    Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num))  +1
                    
          !2022-12-22.           
          if(Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num)) > 20) then
              print *,'    ERROR-2022112302 :: Solid_El_Vertex_Num(in_Ele_Num)>20!'
              print *,'                        In Determine_Enriched_Nodes_3D.f90!'
              print *,'                        Increase size 3 of tem_Solid_El_Vertex_num!'
              call Warning_Message('S',Keywords_Blank)  
          endif    
          
          ! The vertex numbers of cracks i_C contained in each solid element. 2022-08-26.
          tem_Solid_El_Vertex_num(in_Ele_Num,Solid_El_num_Crs(in_Ele_Num),&
                    Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num))) =i_V                          
      endif       
      
      !**********************************************************
      ! Current information on the edge nodes of discrete cracks
      !**********************************************************
      ! Coordinates of discrete fracture edge nodes
      c_x  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,1) 
      c_y  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,2) 
      c_z  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,3) 
      ! Coordinates of the vertices on the boundaries of discrete fracture surfaces contained in each
      ! solid element (up to 3)
      Solid_El_Vertex_Coor(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num),&
                           Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num)),1:3)= [c_x,c_y,c_z]  
      ! The outward normal vectors at the vertices of the discrete fracture plane boundaries contained in
      ! each solid element (up to 3 vertices)
      Solid_El_Vertex_Nor_Vec(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num),&
            Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num)),1:3)=Crack3D_Meshed_Node_Nor_Vector(i_C)%row(i_V,1:3)     
      ! Local x-axis vectors of the discrete fracture face boundary vertices contained in each solid
      ! element (up to 3 vertices)
      Solid_El_Vertex_x_Vec(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num),&
            Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num)),1:3)=Crack3D_Meshed_Vertex_x_Vector(i_C)%row(i_V,1:3)
      ! Local y-axis vectors of the discrete fracture face boundary vertices contained in each solid
      ! element (up to 3 vertices)
      Solid_El_Vertex_y_Vec(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num),&
            Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num)),1:3)=Crack3D_Meshed_Vertex_y_Vector(i_C)%row(i_V,1:3)
      ! Local z-axis vectors of the discrete fracture surface boundary vertices contained in each solid
      ! element (up to 3 vertices)
      Solid_El_Vertex_z_Vec(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num),&
            Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num)),1:3)=Crack3D_Meshed_Vertex_z_Vector(i_C)%row(i_V,1:3)     
      !******************************************************************
      ! Information related to the previous discrete fracture edge nodes
      !******************************************************************
      ! Coordinates of discrete fracture edge nodes
      c_x  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Pre_Node,1) 
      c_y  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Pre_Node,2) 
      c_z  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Pre_Node,3)  
      ! Each solid element contains the coordinates of the vertices of discrete fracture surface
      ! boundaries (up to 3)
      Solid_El_Pre_Vertex_Coor(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num),&
                     Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num)),1:3)= [c_x,c_y,c_z]  
      ! The outward normal vectors at the vertices of the discrete fracture plane boundaries contained in
      ! each solid element (up to 3 vertices)
      Solid_El_Pre_Vertex_Nor_Vec(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num),&
            Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num)),1:3)=&
                 Crack3D_Meshed_Node_Nor_Vector(i_C)%row(c_V_Pre,1:3)     
      ! Local x-axis vectors of the discrete fracture surface boundary vertices contained in each solid
      ! element (up to 3 vertices)
      Solid_El_Pre_Vertex_x_Vec(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num),&
            Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num)),1:3)=&
                 Crack3D_Meshed_Vertex_x_Vector(i_C)%row(c_V_Pre,1:3)
      ! Local y-axis vectors of the discrete fracture face boundary vertices contained in each solid
      ! element (up to 3 vertices)
      Solid_El_Pre_Vertex_y_Vec(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num),&
            Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num)),1:3)=&
                 Crack3D_Meshed_Vertex_y_Vector(i_C)%row(c_V_Pre,1:3)
      ! Local z-axis vectors of the discrete fracture surface boundary vertices contained in each solid
      ! element (up to 3 vertices)
      Solid_El_Pre_Vertex_z_Vec(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num),&
             Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num)),1:3) =&
                Crack3D_Meshed_Vertex_z_Vector(i_C)%row(c_V_Pre,1:3)  


      !**********************************************************
      ! Information related to the next discrete crack edge node
      !**********************************************************
      ! Coordinates of discrete fracture edge nodes
      c_x  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Nex_Node,1) 
      c_y  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Nex_Node,2) 
      c_z  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Nex_Node,3)  
      ! Each solid element contains the coordinates of the vertices of discrete fracture surface
      ! boundaries (up to 3)
      Solid_El_Nex_Vertex_Coor(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num),&
            Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num)),1:3)= [c_x,c_y,c_z]  
      ! The outward normal vectors at the vertices of the discrete fracture plane boundaries contained in
      ! each solid element (up to 3 vertices)
      Solid_El_Nex_Vertex_Nor_Vec(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num),&
            Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num)),1:3) =&
                Crack3D_Meshed_Node_Nor_Vector(i_C)%row(c_V_Nex,1:3)     
      ! Local x-axis vectors of the discrete fracture surface boundary vertices contained in each solid
      ! element (up to 3 vertices)
      Solid_El_Nex_Vertex_x_Vec(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num),&
            Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num)),1:3) =&
                Crack3D_Meshed_Vertex_x_Vector(i_C)%row(c_V_Nex,1:3)
      ! Local y-axis vectors of the discrete fracture face boundary vertices contained in each solid
      ! element (up to 3 vertices)
      Solid_El_Nex_Vertex_y_Vec(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num),&
            Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num)),1:3) =&
                Crack3D_Meshed_Vertex_y_Vector(i_C)%row(c_V_Nex,1:3)
      ! Local z-axis vectors of the discrete fracture face boundary vertices contained in each solid
      ! element (up to 3 vertices)
      Solid_El_Nex_Vertex_z_Vec(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num),&
            Solid_El_Vertex_Num(in_Ele_Num)%row(Solid_El_num_Crs(in_Ele_Num)),1:3) =&
                Crack3D_Meshed_Vertex_z_Vector(i_C)%row(c_V_Nex,1:3)        
  enddo
enddo


!#########################################################
! Take different measures based on the number of vertices
!#########################################################
!BUGFIX2023021901.
do i_E = 1,Num_Elem   
    if(Solid_El_num_Crs(i_E)==0)then
        cycle
    endif
      
    c_NN  = G_NN(1:8,i_E)
    do i_C=1,num_Crack
      ! Check whether this element is related to this crack. IMPROV2022081502.
      call Vector_Location_Int_v2(Solid_El_Max_num_Crs,Solid_El_Crs(i_E,1:Solid_El_Max_num_Crs),i_C,c_Cr_Location)
      !c_Cr_Location =
      !minloc(Solid_El_Crs(i_E,1:Solid_El_num_Crs(i_E)),1,mask=(Solid_El_Crs(i_E,1:Solid_El_num_Crs(i_E))==i_C))
      
      ! If this element is not related to this crack, skip it. IMPROV2022081502.
      if(c_Cr_Location==0)then
          cycle
      endif
      
      num_Vertex = Crack3D_Meshed_Outline_num(i_C)
      
      !*****************************************************************************
      ! CASE 1: Only contains a single discrete fracture surface vertex, 2020-01-01
      !*****************************************************************************
      if(Solid_El_Vertex_Num(i_E)%row(c_Cr_Location)==1) then
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Crack Tip Enhancement Baseline
          ! Store into: Solid_El_Tip_BaseLine(i_E,1:2,1:3)
          ! Algorithm description: If there is one discrete point, the line formed by the discrete
          ! point and the previous discrete point intersects with the solid.
          ! The intersection points of the element, and the intersections of the line formed by a
          ! discrete point and the next discrete point with the solid element,
          ! Together, they form the baseline. The calculation of the crack tip enhancement function
          ! will be carried out on the baseline.
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Calculate the first coordinate point of the baseline
          c_Point_1(1:3) =  Solid_El_Vertex_Coor(i_E)%row(c_Cr_Location,1,1:3)
          c_Point_2(1:3) =  Solid_El_Pre_Vertex_Coor(i_E)%row(c_Cr_Location,1,1:3)
          call Tool_Intersection_of_AB_and_Brick_Ele(c_Point_1,c_Point_2,i_E,c_Yes_Inter,c_InterSection_P) 
          !BUGFIX2022092701. 
          if(c_Yes_Inter) then
              Solid_El_Tip_BaseLine(i_E)%row(c_Cr_Location,1,1:3) =c_InterSection_P
          else
              Solid_El_Tip_BaseLine(i_E)%row(c_Cr_Location,1,1:3) =c_Point_2(1:3) 
          endif
          ! Calculate the second coordinate point of the baseline
          c_Point_1(1:3) =  Solid_El_Vertex_Coor(i_E)%row(c_Cr_Location,1,1:3)
          c_Point_2(1:3) =  Solid_El_Nex_Vertex_Coor(i_E)%row(c_Cr_Location,1,1:3)
          call Tool_Intersection_of_AB_and_Brick_Ele(c_Point_1,c_Point_2,i_E,c_Yes_Inter,c_InterSection_P) 
          !BUGFIX2022092701.
          if(c_Yes_Inter) then
              Solid_El_Tip_BaseLine(i_E)%row(c_Cr_Location,2,1:3) =c_InterSection_P
          else
              Solid_El_Tip_BaseLine(i_E)%row(c_Cr_Location,2,1:3) =c_Point_2(1:3)
          endif
          
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Calculate the outward normal vector of the crack tip enhancement baseline
          ! Store into: Solid_El_Tip_BaseLine_Nor_Vec(1:Num_Elem, 1:3)
          ! Algorithm description: If there is one discrete point, the discrete 
          ! point is connected to the previous discrete point and the next discrete point.
          ! The average of the vectors in the local coordinate systems of the three, 
          ! used as the value to be determined
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          Solid_El_Tip_BaseLine_Nor_Vec(i_E)%row(c_Cr_Location,1:3) = &
                    (Solid_El_Vertex_Nor_Vec(i_E)%row(c_Cr_Location,1,1:3) +     &
                     Solid_El_Pre_Vertex_Nor_Vec(i_E)%row(c_Cr_Location,1,1:3) + &
                     Solid_El_Nex_Vertex_Nor_Vec(i_E)%row(c_Cr_Location,1,1:3))/THR
          call Vector_Normalize(3,Solid_El_Tip_BaseLine_Nor_Vec(i_E)%row(c_Cr_Location,1:3))
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Calculate the local coordinate system vectors of the crack tip enhancement baseline
          ! Store into: Solid_El_Tip_BaseLine_x_Vec(1:Num_Elem,1:3)
          !      Solid_El_Tip_BaseLine_y_Vec(1:Num_Elem,1:3)
          !      Solid_El_Tip_BaseLine_z_Vec(1:Num_Elem,1:3)
          ! Algorithm description: If there is one discrete point, the discrete 
          ! point is connected to the previous discrete point and the next discrete point.
          ! The average of the vectors in the local coordinate systems of the 
          ! three, used as the value to be determined
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          Solid_El_Tip_BaseLine_x_Vec(i_E)%row(c_Cr_Location,1:3) = &
                 (Solid_El_Vertex_x_Vec(i_E)%row(c_Cr_Location,1,1:3)     + &
                  Solid_El_Pre_Vertex_x_Vec(i_E)%row(c_Cr_Location,1,1:3) + &
                  Solid_El_Nex_Vertex_x_Vec(i_E)%row(c_Cr_Location,1,1:3))/THR       
          Solid_El_Tip_BaseLine_y_Vec(i_E)%row(c_Cr_Location,1:3) = &
                 (Solid_El_Vertex_y_Vec(i_E)%row(c_Cr_Location,1,1:3)      + &
                  Solid_El_Pre_Vertex_y_Vec(i_E)%row(c_Cr_Location,1,1:3)  + &
                  Solid_El_Nex_Vertex_y_Vec(i_E)%row(c_Cr_Location,1,1:3))/THR    
          Solid_El_Tip_BaseLine_z_Vec(i_E)%row(c_Cr_Location,1:3) = &
                 (Solid_El_Vertex_z_Vec(i_E)%row(c_Cr_Location,1,1:3)      + &
                  Solid_El_Pre_Vertex_z_Vec(i_E)%row(c_Cr_Location,1,1:3)  + &
                  Solid_El_Nex_Vertex_z_Vec(i_E)%row(c_Cr_Location,1,1:3))/THR    
          call Vector_Normalize(3,Solid_El_Tip_BaseLine_x_Vec(i_E)%row(c_Cr_Location,1:3))
          call Vector_Normalize(3,Solid_El_Tip_BaseLine_y_Vec(i_E)%row(c_Cr_Location,1:3))
          call Vector_Normalize(3,Solid_El_Tip_BaseLine_z_Vec(i_E)%row(c_Cr_Location,1:3))
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Calculate rotation angle and rotation matrix
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          Vector_Crack_x(1:3) = [ONE,ZR,ZR]
          Vector_Crack_y(1:3) = [ZR,ONE,ZR]
          Vector_Crack_z(1:3) = [ZR,ZR,ONE]
          ! Determine ThetaX, ThetaY, and ThetaZ so that all three sets of transformations above are
          ! satisfied.
          call Tool_ThetaX_ThetaY_ThetaZ_3D_rotation(1,                       &
                      Solid_El_Tip_BaseLine_x_Vec(i_E)%row(c_Cr_Location,1:3),&
                      Solid_El_Tip_BaseLine_y_Vec(i_E)%row(c_Cr_Location,1:3),&
                      Solid_El_Tip_BaseLine_z_Vec(i_E)%row(c_Cr_Location,1:3),&
                      Vector_Crack_x(1:3),Vector_Crack_y(1:3),Vector_Crack_z(1:3),&
                      ThetaX,ThetaY,ThetaZ,T_Matrx(1:3,1:3))
          Solid_El_Tip_BaseLine_T_Matrix(i_E)%row(c_Cr_Location,1:3,1:3)=T_Matrx(1:3,1:3)
          Solid_El_Tip_BaseLine_T_theta(i_E)%row(c_Cr_Location,1:3)     =[ThetaX,ThetaY,ThetaZ]
           
          ! If crack tip enhancement is allowed
          if (Key_TipEnrich>=1) then
              ! If the number of cracks associated with this enhancement element is less than
              ! Key_Ele_Max_Related_Cracks, and the current crack has not yet been associated, then execute
              ! IMPROV2022122701.
              if ((Elem_num_Related_Cracks(i_E) <Key_Ele_Max_Related_Cracks) .and. &
                 (.not. any(Elem_Related_Cracks(i_E,1:Key_Ele_Max_Related_Cracks)==i_C)))then   
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  ! Update Elem_num_Related_Cracks(i_E) and Elem_Related_Cracks(i_E,1:)
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  Elem_num_Related_Cracks(i_E)  = Elem_num_Related_Cracks(i_E)  + 1
                  Elem_Related_Cracks(i_E,Elem_num_Related_Cracks(i_E)) = i_C
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  ! Determine crack tip enhancement element
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  Elem_Type_3D(i_E,i_C)   = 1
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  ! Identify crack tip enhancement node
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  Enriched_Node_Type_3D(c_NN,i_C) = 1
                  
                  ! Ele_Num_Tip_Enriched_Node_3D(c_NN,i_C) = i_E ! Enhanced element number corresponding to the
                  ! enriched node (reference element number)
                  ! Using a staggered array, IMPROV2022091804. There may be potential data conflicts, so a CRITICAL
                  ! directive is added.
                  do i_Node =1,8
                      if(.not. allocated(Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node))%row)) then
                          allocate(Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node))%row(num_Crack))
                          Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node))%row(1:num_Crack) = 0
                      endif
                      Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node))%row(i_C) = i_E 
                  enddo
              endif
          endif
          Tem_Elem_Type(i_E,i_C) = 1
      !****************************************************************
      ! CASE 2: Contains 2 discrete crack surface vertices, 2020-01-06
      !****************************************************************
      elseif(Solid_El_Vertex_Num(i_E)%row(c_Cr_Location)==2) then    
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Calculation of crack tip enhancement baseline
          ! Store into: Solid_El_Tip_BaseLine(i_E,1:2,1:3)
          ! Algorithm description: If there are 2 discrete points, the straight
          ! line formed by discrete point 1 and the previous discrete point of 
          ! 1 intersects the solid.
          ! The intersection of the element, and the intersection of the line
          ! formed by discrete point 2 and the next discrete point of 2 
          ! with the solid element,
          ! Together, they form the baseline. The calculation of the crack 
          ! tip enhancement function will be carried out on the baseline.
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~              
          Vertex_1 = Solid_El_Vertex_Coor(i_E)%row(c_Cr_Location,1,1:3)
          Vertex_2 = Solid_El_Vertex_Coor(i_E)%row(c_Cr_Location,2,1:3)
          Vertex_1_Pre = Solid_El_Pre_Vertex_Coor(i_E)%row(c_Cr_Location,1,1:3)
          Vertex_2_Nex = Solid_El_Nex_Vertex_Coor(i_E)%row(c_Cr_Location,2,1:3)
          
          !......................................................................................
          ! Special case handling. Ref: My PhiPsi Development Notebook V1-P76. BUGFIX2022082601.
          !......................................................................................
          ! The vertex numbers of cracks in each solid element i_C. 2022-08-26.
          Vertex_1_num = tem_Solid_El_Vertex_num(i_E,c_Cr_Location,1) 
          Vertex_2_num = tem_Solid_El_Vertex_num(i_E,c_Cr_Location,2) 
          ! Case (1)-(b) in Ref: My PhiPsi Development Notebook V1-P76.
          if (Vertex_1_num == 1 .and. Vertex_2_num ==num_Vertex) then
              ! Change Vertex_1 to the point corresponding to the num_Vertex-th discrete point on the boundary
              ! line of the discrete fracture surface.
              Vertex_1 = Crack3D_Meshed_Node(i_C)%row(Crack3D_Meshed_Outline(i_C)%row(num_Vertex,1),1:3)
              ! Change Vertex_2 to the point corresponding to the first discrete point on the boundary line of the
              ! discrete fracture surface.
              Vertex_2 = Crack3D_Meshed_Node(i_C)%row(Crack3D_Meshed_Outline(i_C)%row(1,1),1:3)
              ! Set the num_Vertex-1 discrete point of the discrete fracture surface boundary line as the point
              ! corresponding to Vertex_1_Pre.
              Vertex_1_Pre = Crack3D_Meshed_Node(i_C)%row(Crack3D_Meshed_Outline(i_C)%row(num_Vertex-1,1),1:3)
              ! Determine the second discrete point on the boundary line of the discrete fracture surface as the
              ! point corresponding to Vertex_2_Nex.
              Vertex_2_Nex = Crack3D_Meshed_Node(i_C)%row(Crack3D_Meshed_Outline(i_C)%row(2,1),1:3)
          elseif(Vertex_1_num == num_Vertex .and. Vertex_2_num ==1)then
              print *,'    ERROR-2022082602 :: Unforseen case occurred!'
              print *,'                        in Determine_Enriched_Nodes_3D.f90!'
              call Warning_Message('S',Keywords_Blank)  
          endif
          ! Calculate the first coordinate point of the baseline
          call Tool_Intersection_of_AB_and_Brick_Ele(Vertex_1,Vertex_1_Pre,i_E,c_Yes_Inter,c_InterSection_P) 
          Solid_El_Tip_BaseLine(i_E)%row(c_Cr_Location,1,1:3) =c_InterSection_P
          ! Calculate the second coordinate point of the baseline
          call Tool_Intersection_of_AB_and_Brick_Ele(Vertex_2,Vertex_2_Nex,i_E,c_Yes_Inter,c_InterSection_P) 
          Solid_El_Tip_BaseLine(i_E)%row(c_Cr_Location,2,1:3) =c_InterSection_P             

          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Calculate the outward normal vector of the crack tip enhancement baseline.
          ! BUGFIX2022081501.
          ! Store into: Solid_El_Tip_BaseLine_Nor_Vec(1:Num_Elem, 1:3)
          ! Algorithm description: If there are 2 discrete points, the average of the 2 discrete
          ! points is used as the value to be calculated.
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          Solid_El_Tip_BaseLine_Nor_Vec(i_E)%row(c_Cr_Location,1:3) = &
                       (Solid_El_Vertex_Nor_Vec(i_E)%row(c_Cr_Location,1,1:3) +  &
                        Solid_El_Vertex_Nor_Vec(i_E)%row(c_Cr_Location,2,1:3))/TWO
            
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Calculate the local coordinate system vectors of the crack tip enhancement baseline
          ! Algorithm description: If there are 2 discrete points, then discrete point 1 and
          ! discrete point 2
          ! The average of the vectors in the local coordinate systems of the two, used as the value
          ! to be determined
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          Solid_El_Tip_BaseLine_x_Vec(i_E)%row(c_Cr_Location,1:3) = &
                (Solid_El_Vertex_x_Vec(i_E)%row(c_Cr_Location,1,1:3)  +   &
                 Solid_El_Vertex_x_Vec(i_E)%row(c_Cr_Location,2,1:3))/TWO
          Solid_El_Tip_BaseLine_y_Vec(i_E)%row(c_Cr_Location,1:3) = &
                (Solid_El_Vertex_y_Vec(i_E)%row(c_Cr_Location,1,1:3)  +   &
                 Solid_El_Vertex_y_Vec(i_E)%row(c_Cr_Location,2,1:3))/TWO
          Solid_El_Tip_BaseLine_z_Vec(i_E)%row(c_Cr_Location,1:3) =&
                (Solid_El_Vertex_z_Vec(i_E)%row(c_Cr_Location,1,1:3)  +   &
                 Solid_El_Vertex_z_Vec(i_E)%row(c_Cr_Location,2,1:3))/TWO   
          call Vector_Normalize(3,Solid_El_Tip_BaseLine_x_Vec(i_E)%row(c_Cr_Location,1:3))
          call Vector_Normalize(3,Solid_El_Tip_BaseLine_y_Vec(i_E)%row(c_Cr_Location,1:3))
          call Vector_Normalize(3,Solid_El_Tip_BaseLine_z_Vec(i_E)%row(c_Cr_Location,1:3))
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Calculate rotation angle and rotation matrix
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          Vector_Crack_x(1:3) = [ONE,ZR,ZR]
          Vector_Crack_y(1:3) = [ZR,ONE,ZR]
          Vector_Crack_z(1:3) = [ZR,ZR,ONE]
          call Tool_ThetaX_ThetaY_ThetaZ_3D_rotation(1,            &
                      Solid_El_Tip_BaseLine_x_Vec(i_E)%row(c_Cr_Location,1:3),    &
                      Solid_El_Tip_BaseLine_y_Vec(i_E)%row(c_Cr_Location,1:3),    & 
                      Solid_El_Tip_BaseLine_z_Vec(i_E)%row(c_Cr_Location,1:3),    & 
                      Vector_Crack_x(1:3),Vector_Crack_y(1:3),Vector_Crack_z(1:3),& 
                      ThetaX,ThetaY,ThetaZ,T_Matrx(1:3,1:3))
          Solid_El_Tip_BaseLine_T_Matrix(i_E)%row(c_Cr_Location,1:3,1:3)= T_Matrx(1:3,1:3)
          Solid_El_Tip_BaseLine_T_theta(i_E)%row(c_Cr_Location,1:3)  =[ThetaX,ThetaY,ThetaZ]
          
          ! If crack tip enhancement is allowed
          if (Key_TipEnrich>=1) then     
              ! If the number of cracks associated with this enhancement element is less than
              ! Key_Ele_Max_Related_Cracks, and the current crack has not yet been associated, then execute
              ! IMPROV2022122701.
              if ((Elem_num_Related_Cracks(i_E) <Key_Ele_Max_Related_Cracks) .and.&
                  (.not. any(Elem_Related_Cracks(i_E,1:Key_Ele_Max_Related_Cracks)==i_C)))then   
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  ! Update Elem_num_Related_Cracks(i_E) and Elem_Related_Cracks(i_E,1:4)
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  Elem_num_Related_Cracks(i_E)  = Elem_num_Related_Cracks(i_E)  + 1
                  Elem_Related_Cracks(i_E,Elem_num_Related_Cracks(i_E)) = i_C
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  ! Determine crack tip enhancement element
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  Elem_Type_3D(i_E,i_C)   = 1 
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  ! Identify crack tip enhancement node
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  Enriched_Node_Type_3D(c_NN,i_C) = 1    
                   
                  ! Ele_Num_Tip_Enriched_Node_3D(c_NN,i_C) = i_E  ! Enhanced element 
                  ! number corresponding to the enriched node (reference element number)
                  ! Using a staggered array, IMPROV2022091804. There may be potential 
                  ! data conflicts, so a CRITICAL directive is added.
                  do i_Node =1,8
                      if(.not. allocated(Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node))%row)) then
                          allocate(Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node))%row(num_Crack))
                          Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node))%row(1:num_Crack) = 0
                      endif
                      Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node))%row(i_C) = i_E 
                  enddo      
              endif
          endif
          Tem_Elem_Type(i_E,i_C) = 1
          
      !*********************************************************************************************
      ! CASE 3: Contains 0 discrete fracture plane vertices, 2022-08-03. *
      ! Explanation: (1) In most cases at this time, the crack surface outline 
      ! line segment passes completely through the element.
      ! Therefore, there are no crack surfaces at the vertices within the element.
      ! (2) The outline line segment of the crack surface will intersect the element at two points.
      ! (3) Mark as crack enhancement element. *
      !*********************************************************************************************
      elseif(Solid_El_Vertex_Num(i_E)%row(c_Cr_Location)==0) then    
          ! Note: Cases with 0 discrete fracture surface vertices will be handled separately in Step 3.1.
          

          
      !****************************************************************
      ! CASE 4: Contains 3 discrete crack surface vertices, 2020-01-06
      !****************************************************************
      elseif(Solid_El_Vertex_Num(i_E)%row(c_Cr_Location)==3) then
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Calculation of crack tip enhancement baseline
          ! Store into: Solid_El_Tip_BaseLine(i_E,1:2,1:3)
          ! Algorithm description: If there are three discrete points, the 
          ! line formed by discrete point 1 and the previous discrete point 
          ! of 1 intersects with the solid
          ! The intersection of the element, and the intersection of the 
          ! line formed by discrete point 3 and the next discrete point 
          ! after 3 with the solid element,
          ! Together, they form the baseline. The calculation of the crack
          ! tip enhancement function will be carried out on the baseline.
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~              
          Vertex_1 = Solid_El_Vertex_Coor(i_E)%row(c_Cr_Location,1,1:3)
          Vertex_3 = Solid_El_Vertex_Coor(i_E)%row(c_Cr_Location,3,1:3)
          Vertex_1_Pre = Solid_El_Pre_Vertex_Coor(i_E)%row(c_Cr_Location,1,1:3)
          Vertex_3_Nex = Solid_El_Nex_Vertex_Coor(i_E)%row(c_Cr_Location,3,1:3)  

          !....................................................................
          ! Special case handling. Ref: My PhiPsi Development Notebook V1-P76. 
          ! BUGFIX2022082601.
          !....................................................................
          ! The vertex numbers of cracks i_C contained in each solid element. 2022-08-26.
          ! Note that Vertex_1_num, Vertex_2_num, and Vertex_3_num must be arranged in ascending order.
          Vertex_1_num = tem_Solid_El_Vertex_num(i_E,c_Cr_Location,1) 
          Vertex_2_num = tem_Solid_El_Vertex_num(i_E,c_Cr_Location,2) 
          Vertex_3_num = tem_Solid_El_Vertex_num(i_E,c_Cr_Location,3) 
          ! Case (2)-(a) in Ref: My PhiPsi Development Notebook V1-P76.
          if (Vertex_1_num == 1 .and. Vertex_2_num ==2 .and. Vertex_3_num ==num_Vertex) then
              ! Change Vertex_1 to the point corresponding to the num_Vertex-th discrete point on the boundary
              ! line of the discrete fracture surface.
              Vertex_1 = Crack3D_Meshed_Node(i_C)%row(Crack3D_Meshed_Outline(i_C)%row(num_Vertex,1),1:3)   
              ! Change Vertex_3 to the point corresponding to the second discrete point on the boundary line of
              ! the discrete fracture surface.
              Vertex_3 = Crack3D_Meshed_Node(i_C)%row(Crack3D_Meshed_Outline(i_C)%row(2,1),1:3)              
              ! Determine the num_Vertex-1 discrete point of the discrete fracture surface boundary line as the
              ! point corresponding to Vertex_1_Pre.
              Vertex_1_Pre = Crack3D_Meshed_Node(i_C)%row(Crack3D_Meshed_Outline(i_C)%row(num_Vertex-1,1),1:3)
              ! Determine the three discrete points on the boundary line of the discrete fracture surface as the
              ! points corresponding to Vertex_3_Nex.
              Vertex_3_Nex = Crack3D_Meshed_Node(i_C)%row(Crack3D_Meshed_Outline(i_C)%row(3,1),1:3)
          ! Case (2)-(b) in Ref: My PhiPsi Development Notebook V1-P76.
          elseif(Vertex_1_num == 1 .and. (Vertex_2_num ==num_Vertex-1) .and. (Vertex_3_num ==num_Vertex)) then
              ! Change Vertex_1 to the point corresponding to the num_Vertex-th discrete point on the boundary
              ! line of the discrete fracture surface.
              Vertex_1 = Crack3D_Meshed_Node(i_C)%row(Crack3D_Meshed_Outline(i_C)%row(num_Vertex-1,1),1:3)           
              ! Vertex_3 remains stationary.
              !Vertex_3 = Vertex_1              
              ! Set the (num_Vertex-2)th discrete point on the boundary line of the discrete fracture surface as
              ! the point corresponding to Vertex_1_Pre.
              Vertex_1_Pre = Crack3D_Meshed_Node(i_C)%row(Crack3D_Meshed_Outline(i_C)%row(num_Vertex-2,1),1:3)
              ! Determine the second discrete point on the boundary line of the discrete fracture surface as the
              ! point corresponding to Vertex_3_Nex.
              Vertex_3_Nex = Crack3D_Meshed_Node(i_C)%row(Crack3D_Meshed_Outline(i_C)%row(2,1),1:3)
          endif
          
          
          ! Calculate the first coordinate point of the baseline
          call Tool_Intersection_of_AB_and_Brick_Ele(Vertex_1,Vertex_1_Pre,i_E,c_Yes_Inter,c_InterSection_P) 
          Solid_El_Tip_BaseLine(i_E)%row(c_Cr_Location,1,1:3) =c_InterSection_P
          ! Calculate the second coordinate point of the baseline
          call Tool_Intersection_of_AB_and_Brick_Ele(Vertex_3,Vertex_3_Nex,i_E,c_Yes_Inter,c_InterSection_P) 
          Solid_El_Tip_BaseLine(i_E)%row(c_Cr_Location,2,1:3) =c_InterSection_P             
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Calculate the outward normal vector of the crack tip enhancement
          ! baseline. BUGFIX2022081501.
          ! Store into: Solid_El_Tip_BaseLine_Nor_Vec(1:Num_Elem, 1:3)
          ! Algorithm description: If there are 3 discrete points, the average 
          ! of the 3 discrete points is used as the value to be calculated.
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          Solid_El_Tip_BaseLine_Nor_Vec(i_E)%row(c_Cr_Location,1:3) =&
                (Solid_El_Vertex_Nor_Vec(i_E)%row(c_Cr_Location,1,1:3) +  &
                 Solid_El_Vertex_Nor_Vec(i_E)%row(c_Cr_Location,2,1:3) +  &
                 Solid_El_Vertex_Nor_Vec(i_E)%row(c_Cr_Location,3,1:3))/THR
            
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Calculate the local coordinate system vectors of the crack tip 
          ! enhancement baseline
          ! Algorithm description: If there are 3 discrete points, then 
          ! discrete point 1 and discrete points 2 and 3
          ! The average of the vectors in the local coordinate systems of
          ! the three, used as the value to be determined
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          Solid_El_Tip_BaseLine_x_Vec(i_E)%row(c_Cr_Location,1:3) =&
                    (Solid_El_Vertex_x_Vec(i_E)%row(c_Cr_Location,1,1:3)  + &
                     Solid_El_Vertex_x_Vec(i_E)%row(c_Cr_Location,2,1:3)  + &
                     Solid_El_Vertex_x_Vec(i_E)%row(c_Cr_Location,3,1:3))/THR
          Solid_El_Tip_BaseLine_y_Vec(i_E)%row(c_Cr_Location,1:3) = &
                    (Solid_El_Vertex_y_Vec(i_E)%row(c_Cr_Location,1,1:3)  + &
                    Solid_El_Vertex_y_Vec(i_E)%row(c_Cr_Location,2,1:3)   + &
                    Solid_El_Vertex_y_Vec(i_E)%row(c_Cr_Location,3,1:3))/THR
          Solid_El_Tip_BaseLine_z_Vec(i_E)%row(c_Cr_Location,1:3) = &
                    (Solid_El_Vertex_z_Vec(i_E)%row(c_Cr_Location,1,1:3)  + &
                     Solid_El_Vertex_z_Vec(i_E)%row(c_Cr_Location,2,1:3)  + &
                     Solid_El_Vertex_z_Vec(i_E)%row(c_Cr_Location,3,1:3))/THR   
          call Vector_Normalize(3,Solid_El_Tip_BaseLine_x_Vec(i_E)%row(c_Cr_Location,1:3))
          call Vector_Normalize(3,Solid_El_Tip_BaseLine_y_Vec(i_E)%row(c_Cr_Location,1:3))
          call Vector_Normalize(3,Solid_El_Tip_BaseLine_z_Vec(i_E)%row(c_Cr_Location,1:3))
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Calculate the rotation angle and rotation matrix
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          Vector_Crack_x(1:3) = [ONE,ZR,ZR]
          Vector_Crack_y(1:3) = [ZR,ONE,ZR]
          Vector_Crack_z(1:3) = [ZR,ZR,ONE]
          ! Determine ThetaX, ThetaY, and ThetaZ so that all three sets of transformations above are
          ! satisfied.
          call Tool_ThetaX_ThetaY_ThetaZ_3D_rotation(1,            &
                      Solid_El_Tip_BaseLine_x_Vec(i_E)%row(c_Cr_Location,1:3),     &
                      Solid_El_Tip_BaseLine_y_Vec(i_E)%row(c_Cr_Location,1:3),     & 
                      Solid_El_Tip_BaseLine_z_Vec(i_E)%row(c_Cr_Location,1:3),     & 
                      Vector_Crack_x(1:3),Vector_Crack_y(1:3),Vector_Crack_z(1:3), & 
                      ThetaX,ThetaY,ThetaZ,T_Matrx(1:3,1:3))
          Solid_El_Tip_BaseLine_T_Matrix(i_E)%row(c_Cr_Location,1:3,1:3)= T_Matrx(1:3,1:3)
          Solid_El_Tip_BaseLine_T_theta(i_E)%row(c_Cr_Location,1:3)  =[ThetaX,ThetaY,ThetaZ]
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! Given crack tip enhancement
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          ! If crack tip enhancement is allowed
          if (Key_TipEnrich>=1) then    
              ! If the number of cracks associated with this enhancement element is less 
              ! than Key_Ele_Max_Related_Cracks, and the current crack has not yet been 
              ! associated, then execute IMPROV2022122701.
              if ((Elem_num_Related_Cracks(i_E) <Key_Ele_Max_Related_Cracks) .and.&
                  (.not. any(Elem_Related_Cracks(i_E,1:Key_Ele_Max_Related_Cracks)==i_C)))then   
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  ! Update Elem_num_Related_Cracks(i_E) and Elem_Related_Cracks(i_E,1:4)
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  Elem_num_Related_Cracks(i_E)  = Elem_num_Related_Cracks(i_E)  + 1
                  Elem_Related_Cracks(i_E,Elem_num_Related_Cracks(i_E)) = i_C
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  ! Determine crack tip enhancement element
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  Elem_Type_3D(i_E,i_C)   = 1 
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  ! Identify crack tip enhancement node
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  Enriched_Node_Type_3D(c_NN,i_C) = 1    
                  
                  ! Using a staggered array, IMPROV2022091804. There may be potential data conflicts, so a CRITICAL
                  ! directive is added.
                  do i_Node =1,8
                      if(.not. allocated(Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node))%row)) then
                          allocate(Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node))%row(num_Crack))
                          Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node))%row(1:num_Crack) = 0
                      endif
                      Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node))%row(i_C) = i_E 
                  enddo 
              endif
          endif
          Tem_Elem_Type(i_E,i_C) = 1
      !*******************************************************************
      ! CASE 5: Contains 4 discrete fracture surface vertices, 2022-08-26
      !*******************************************************************
      elseif(Solid_El_Vertex_Num(i_E)%row(c_Cr_Location)==4) then
          print *,'    Warn-2022102301 :: Crack tip element ', i_E ,' with 4 vertexes found!'
          print *,'                       In Determine_Enriched_Nodes_3D.f90!'     
      !*********************************************************************
      ! CASE 6: Contains 5 or more discrete crack face vertices, 2022-10-23
      !*********************************************************************
      elseif(Solid_El_Vertex_Num(i_E)%row(c_Cr_Location)>=5) then    
          print *,'    Warn-2022102302 :: Crack tip element ', i_E ,' with 5 vertexes found!'
          print *,'                       In Determine_Enriched_Nodes_3D.f90'             
      endif
    enddo
enddo
!!$omp end parallel do          

      
!.........................................................................
!                                                                       .
!                                                                       .
! Step 2: Calculate the position of each node relative to each fracture, 
! i.e., the signed distance (positive or negative), as well as the foot 
! coordinates.
! Used to determine Heaviside enrichment elements and nodes.
!                                                                       .
!                                                                       .
!.........................................................................
! IMPROV2022072205. OpenMP Parallelization.
!$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_C,i_Node,c_Node_Coor,c_Min_Signed_Dis,Check_Ball_R, &
!$OMP                  Yes_Found_Min_Signed_Dis,c_n_Vector,c_Yes_Node_PER_in_FS,c_PER_Node_to_FS,&
!$OMP                  c_Dis_Node_to_FS,c_Dis_Node_to_FS_v2)&
!$OMP SCHEDULE(static)       
!do i_C=1,num_Crack
do i_Node = 1,num_Node
  !*********************************
  ! New code for discrete fractures
  !*********************************
  !do i_Node = 1,num_Node
  do i_C=1,num_Crack
      c_Node_Coor = Coor(i_Node,1:3)
      c_Min_Signed_Dis  = Con_Big_20
      ! Calculate the signed distance from the current node to the fracture plane and determine the foot
      ! position
      ! Check_Ball_R = 3.0D0*Ave_Elem_L !Only calculate the discrete triangles of the crack surface within
      ! the detection
      Check_Ball_R  = 3.0D0*Node_Max_L(i_Node)
      
      ! If out-of-plane extension.
      if(Key_InPlane_Growth== 0 ) then
          call D3_Get_Signed_Dis_to_Crack_Mesh(c_Node_Coor,i_C,Check_Ball_R,  &
                    c_Dis_Node_to_FS,c_Dis_Node_to_FS_v2,                     &
                    c_Yes_Node_PER_in_FS,c_PER_Node_to_FS(1:3),               &
                    Yes_Found_Min_Signed_Dis,c_n_Vector)
                    
          !Dis_Node_to_FS(i_Node,i_C)= c_Dis_Node_to_FS          
          if(.not. allocated(Dis_Node_to_FS(i_Node)%row))then
            allocate(Dis_Node_to_FS(i_Node)%row(num_Crack))
            Dis_Node_to_FS(i_Node)%row(1:num_Crack) = ZR
          endif
          Dis_Node_to_FS(i_Node)%row(i_C) =  c_Dis_Node_to_FS  
          
          !Dis_Node_to_FS_v2(i_Node,i_C) = c_Dis_Node_to_FS_v2
          
          ! Using jagged arrays. 2022-09-15. IMPROV2022091502.
          if(abs(c_Dis_Node_to_FS_v2) > Tol_10) then
              if(.not. allocated(RaggedArray1D_Dis_Node_to_FS_v2(i_Node)%row))then  
                  allocate(RaggedArray1D_Dis_Node_to_FS_v2(i_Node)%row(num_Crack))
                  RaggedArray1D_Dis_Node_to_FS_v2(i_Node)%row(1:num_Crack) = ZR
                  !RaggedArray1D_Dis_Node_to_FS_v2(i_Node)%row(i_C) = c_Dis_Node_to_FS_v2
              endif
              RaggedArray1D_Dis_Node_to_FS_v2(i_Node)%row(i_C) = c_Dis_Node_to_FS_v2
          endif
      ! For in-plane extension, there is no need to find the minimum symbol distance, resulting in low
      ! resource consumption. 2022-09-08. NEWFTU2022090801.
      elseif(Key_InPlane_Growth== 1) then
          call D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane_Growth(c_Node_Coor,i_C,&
                                   Check_Ball_R,c_Dis_Node_to_FS ,                &
                                   c_Yes_Node_PER_in_FS,c_PER_Node_to_FS(1:3),    &
                                   Yes_Found_Min_Signed_Dis,c_n_Vector)
                                   
          !Dis_Node_to_FS(i_Node,i_C) =  c_Dis_Node_to_FS    
          if(.not. allocated(Dis_Node_to_FS(i_Node)%row))then
              allocate(Dis_Node_to_FS(i_Node)%row(num_Crack))
              Dis_Node_to_FS(i_Node)%row(1:num_Crack) = ZR
          endif
          Dis_Node_to_FS(i_Node)%row(i_C) =  c_Dis_Node_to_FS  
          
          !if(i_Node==2802) then
          !endif
          
          !Dis_Node_to_FS_v2(i_Node,i_C) = Dis_Node_to_FS(i_Node,i_C)
          
          ! Using jagged arrays. 2022-09-15. IMPROV2022091502.
          if(abs(c_Dis_Node_to_FS ) > Tol_10) then
              if(.not. allocated(RaggedArray1D_Dis_Node_to_FS_v2(i_Node)%row))then  
                  allocate(RaggedArray1D_Dis_Node_to_FS_v2(i_Node)%row(num_Crack))
                  RaggedArray1D_Dis_Node_to_FS_v2(i_Node)%row(1:num_Crack)= ZR
                  !RaggedArray1D_Dis_Node_to_FS_v2(i_Node)%row(i_C) = Dis_Node_to_FS(i_Node,i_C) 
              endif
              RaggedArray1D_Dis_Node_to_FS_v2(i_Node)%row(i_C) = c_Dis_Node_to_FS
          endif
          
      endif
      
      ! n_Vectors(i_Node,i_C,1:3) = c_n_Vector                        !Save the n_Vector for each node, 2022-05-12.

      if(sum(abs(c_n_Vector))>Tol_10) then
          if(.not. allocated(Ragged_Array_2D_n_Vectors(i_Node)%row))then  
            allocate(Ragged_Array_2D_n_Vectors(i_Node)%row(num_Crack,3))
            ! Ragged_Array_2D_n_Vectors(i_Node)%row(1:num_Crack,1:3) = ZR !Major bug fix, before the bug fix.
            ! BUGFIX2022102101.
            Ragged_Array_2D_n_Vectors(i_Node)%row(i_C,1:3) = c_n_Vector
          ! BUGFIX2022122301. On 2022-12-23, the following two key lines of code were added.
          else
            Ragged_Array_2D_n_Vectors(i_Node)%row(i_C,1:3) = c_n_Vector
          endif
      endif
  enddo
enddo
!$omp end parallel do        
      
!..........................................................................
! Step 3: Determine the type of enriched Heaviside element and the crack tip
!         enriched element without vertices.
! If the nodes of this element include both positive and negative nodes, there are two cases:
! (1) If the element is not intersected by the boundary of a 3D crack, it 
!     is a Heaviside-enriched element;
! (2) If the element is intersected by the boundary of a 3D crack, it is 
!     considered a crack-tip enhanced unit (the intersection point is taken
!     as the reference line).
! This part is not processed, see Step 3.2 for details.
!        ---------------                                                  .
! Note: The crack-tip enhancement elements determined in this section are 
!       referred to as Type 2 crack-tip enhancement elements (excluding discrete cracks).
!..........................................................................
! IMPROV2022072205. OpenMP parallelization. Extensive testing required.
!$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_C,i_E,c_NN,c_max,c_min,c_Yes_Inter,c_Inter_Point_A,c_Inter_Point_B,  &
!$OMP                  c_Vertex_1,c_Vertex_2,num_Inter_Outline,c_Four_Points,c_Eight_Values)      &
!$OMP SCHEDULE(STATIC)   
!do i_C=1,num_Crack
!  do i_E = 1,Num_Elem  
  
!BUGFIX2023021901.
do i_E = 1,Num_Elem  
  do i_C=1,num_Crack  
     c_NN  = G_NN(:,i_E)
     
     ! Using jagged arrays. 2022-09-18. IMPROV2022091802.
     c_Eight_Values(1:8) = ZR
     do i_Node = 1,8
         if(allocated(Dis_Node_to_FS(c_NN(i_Node))%row))then 
             c_Eight_Values(i_Node) = Dis_Node_to_FS(c_NN(i_Node))%row(i_C)
         endif
     enddo
     c_max = maxval(c_Eight_Values(1:8))     
     c_min = minval(c_Eight_Values(1:8)) 
        
     ! If the element's nodes include both positive and negative nodes
     !if( c_max > ZR .and. c_min < ZR)then
     if( c_max > Tol_11 .and. c_min < -Tol_11)then
       ! First, ensure that the element has not been enhanced (relative to i_C)
       !if(Elem_Type_3D(i_E,i_C) ==0)then
       if(Elem_Type_3D(i_E,i_C) ==0 .and. Tem_Elem_Type(i_E,i_C) ==0)then
         ! Determine whether the 3D element is intersected by the outer boundary of a 3D discrete fracture, 
         ! and calculate the coordinates of the intersection points
         ! Note: It is a boundary line that runs straight through.
         c_Yes_Inter =.False.
         ! HOTCOD2023021802. Hotfix code (fixed). 2023-02-18.
         call Tool_Yes_Ele_Intersected_by_3D_Crack_Outline(i_E,i_C,c_Yes_Inter,c_Inter_Point_A,c_Inter_Point_B, &
                        c_Vertex_1,c_Vertex_2,num_Inter_Outline,c_Four_Points)

         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! If the element is not intersected by an individual 3D crack
         ! boundary, it is a Heaviside-enriched element.
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         if(c_Yes_Inter .eqv. .False.)then
              ! If the number of cracks associated with this enhancement element is less than
              ! Key_Ele_Max_Related_Cracks,
              ! and the current crack has not yet been associated, then execute IMPROV2022122701.
              if ((Elem_num_Related_Cracks(i_E) <Key_Ele_Max_Related_Cracks) .and. &
                  (.not. any(Elem_Related_Cracks(i_E,1:Key_Ele_Max_Related_Cracks)==i_C)))then   
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  ! Update Elem_num_Related_Cracks(i_E) and Elem_Related_Cracks(i_E,1:4)
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  Elem_num_Related_Cracks(i_E)  = Elem_num_Related_Cracks(i_E)  + 1
                  Elem_Related_Cracks(i_E,Elem_num_Related_Cracks(i_E)) = i_C
                  ! Update unit type.
                  Elem_Type_3D(i_E,i_C)  = 2 
                  Tem_Elem_Type(i_E,i_C) = 2
              endif
         endif
       endif
     endif
  enddo
enddo
!$omp end parallel do  


!...........................................................................................
!                                                                                          .
!                                                                                          .
! Step 3.1: Identify crack tip enrichment elements without Vertex. 2022-08-04. BUGFIX2022080401.
! If the nodes of this element include both positive and negative nodes, there are two cases:
! (1) If the element is intersected by the boundary of a 3D crack, it is considered a
!     crack-tip enriched element (the intersection point is taken as the reference line);
! (2) Note: When calculating the nodal symbol distance, be sure not to consider whether 
!     the foot of the perpendicular is on the crack surface.
! Use Dis_Node_to_FS_v2(c_NN, i_C) instead of Dis_Node_to_FS(c_NN, i_C)!
! Because if the element node is on the fracture surface, the signed distance is 0,
! and the perpendicular foot of other nodes may not be on the element surface at this time.
! The symbol distance of all eight unit nodes may be 0, which does not match the actual 
! situation, causing errors in the detection of crack tip enhancement elements.
!        ---------------                                                                   .
! Note: The crack-tip enhancement elements determined in this section are referred to as 
!       Type 2 crack-tip enhancement elements (excluding discrete cracks).
!                                                                                          .
!                                                                                          .
!...........................................................................................
!BUGFIX2023021901.
do i_E = 1,Num_Elem      
  do i_C=1,num_Crack
     c_NN  = G_NN(1:8,i_E)
     ! Using jagged arrays. 2022-09-15. IMPROV2022091502.
     c_Eight_Values(1:8) = ZR
     do i_Node = 1,8
         if(allocated(RaggedArray1D_Dis_Node_to_FS_v2(c_NN(i_Node))%row))then 
             c_Eight_Values(i_Node) = RaggedArray1D_Dis_Node_to_FS_v2(c_NN(i_Node))%row(i_C)
         endif
     enddo
     c_max = maxval(c_Eight_Values(1:8))     
     c_min = minval(c_Eight_Values(1:8)) 
     
     ! If the element's nodes include both positive and negative nodes
     if(c_max > Tol_11 .and. c_min < -Tol_11)then
       ! First, ensure that the element has not been enhanced (relative to i_C)
       !if(Elem_Type_3D(i_E,i_C) ==0)then
       if(Elem_Type_3D(i_E,i_C) ==0 .and. Tem_Elem_Type(i_E,i_C) ==0  )then
         ! Determine whether the 3D element is intersected by the outer boundary of 
         ! a 3D discrete fracture, and calculate the coordinates of the intersection points
         ! Note: It is a boundary line that runs straight through.
         c_Yes_Inter =.False.
         call Tool_Yes_Ele_Intersected_by_3D_Crack_Outline(i_E,i_C,c_Yes_Inter,c_Inter_Point_A,c_Inter_Point_B, &
                        c_Vertex_1,c_Vertex_2,num_Inter_Outline,c_Four_Points)
         
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! If the element is crossed by the boundary of a 
         ! 3D crack, and only one boundary line passes through,
         ! then it becomes a crack tip enhancement element
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         !if (c_Yes_Inter .eqv. .True.)then
         if ((c_Yes_Inter .eqv. .True.) .and.  num_Inter_Outline==1)then
           ! Check whether this element has already been marked for this crack before. 2022-08-15.
           call Vector_Location_Int_v2(Solid_El_Max_num_Crs,Solid_El_Crs(i_E,1:Solid_El_Max_num_Crs),i_C,c_Cr_Location)
           !c_Cr_Location =
           !minloc(Solid_El_Crs(i_E,1:Solid_El_num_Crs(i_E)),1,mask=(Solid_El_Crs(i_E,1:Solid_El_num_Crs(i_E))==i_C))
           ! If the element had not previously marked the crack. 2022-08-15.
           !if(c_Cr_Location==0)then
           if(c_Cr_Location==0 .or. Solid_El_num_Crs(i_E)==0)then
               ! Allocate memory if the i_E ragged array object is not assigned. 2022-09-04.
               if(.not. ALLOCATED(Solid_El_Vertex_Num(i_E)%row)) then
                   ! Allocate memory for the variables related to crack tip enhancement for this element and initialize
                   ! them. 2022-09-04. IMPROV2022090401.
                   call D3_Allocate_Ele_Memory(i_E,1)   
               endif                    
               ! Update the number of cracks associated with each solid element. 2022-08-15. IMPROV2022081502.
               Solid_El_num_Crs(i_E) = Solid_El_num_Crs(i_E) + 1
               
               !2022-11-23.
               if(Solid_El_num_Crs(i_E) > Solid_El_Max_num_Crs) then
                  print *,'    ERROR-2022112301 :: Solid_El_num_Crs(i_E)>Solid_El_Max_num_Crs!'
                  print *,'                        in Determine_Enriched_Nodes_3D.f90!'
                  print *,'                        Try to increase Solid_El_Max_num_Crs!'
                  call Warning_Message('S',Keywords_Blank)  
               endif
              
               ! Update the crack number related to each solid element. 2022-08-15. IMPROV2022081502.
               Solid_El_Crs(i_E,Solid_El_num_Crs(i_E)) = i_C
              
               ! The first coordinate point of the baseline
               Solid_El_Tip_BaseLine(i_E)%row(Solid_El_num_Crs(i_E),1,1:3)=c_Inter_Point_A
               ! Calculate the second coordinate point of the baseline
               Solid_El_Tip_BaseLine(i_E)%row(Solid_El_num_Crs(i_E),2,1:3)=c_Inter_Point_B                  
               ! Calculate the local coordinate system vectors of the crack tip enhancement baseline
               Solid_El_Tip_BaseLine_x_Vec(i_E)%row(Solid_El_num_Crs(i_E),1:3) = &
                                                     (Crack3D_Meshed_Vertex_x_Vector(i_C)%row(c_Vertex_1,1:3)+ &
                                                      Crack3D_Meshed_Vertex_x_Vector(i_C)%row(c_Vertex_2,1:3))/TWO
               Solid_El_Tip_BaseLine_y_Vec(i_E)%row(Solid_El_num_Crs(i_E),1:3) = &
                                                     (Crack3D_Meshed_Vertex_y_Vector(i_C)%row(c_Vertex_1,1:3)+ &
                                                      Crack3D_Meshed_Vertex_y_Vector(i_C)%row(c_Vertex_2,1:3))/TWO
               Solid_El_Tip_BaseLine_z_Vec(i_E)%row(Solid_El_num_Crs(i_E),1:3) = &
                                                     (Crack3D_Meshed_Vertex_z_Vector(i_C)%row(c_Vertex_1,1:3)+ &
                                                      Crack3D_Meshed_Vertex_z_Vector(i_C)%row(c_Vertex_2,1:3))/TWO 
               call Vector_Normalize(3,Solid_El_Tip_BaseLine_x_Vec(i_E)%row(Solid_El_num_Crs(i_E),1:3))
               call Vector_Normalize(3,Solid_El_Tip_BaseLine_y_Vec(i_E)%row(Solid_El_num_Crs(i_E),1:3))
               call Vector_Normalize(3,Solid_El_Tip_BaseLine_z_Vec(i_E)%row(Solid_El_num_Crs(i_E),1:3))
               ! Calculate the normal vector of the crack tip enhancement baseline. 2022-08-15. BUGFIX2022081501.
               Solid_El_Tip_BaseLine_Nor_Vec(i_E)%row(Solid_El_num_Crs(i_E),1:3) = &
                                  (Crack3D_Meshed_Node_Nor_Vector(i_C)%row(c_Vertex_1,1:3)+ &
                                   Crack3D_Meshed_Node_Nor_Vector(i_C)%row(c_Vertex_2,1:3))/TWO
               call Vector_Normalize(3,Solid_El_Tip_BaseLine_Nor_Vec(i_E)%row(Solid_El_num_Crs(i_E),1:3))
               ! Calculate rotation angle and rotation matrix
               Vector_Crack_x(1:3) = [ONE,ZR,ZR]
               Vector_Crack_y(1:3) = [ZR,ONE,ZR]
               Vector_Crack_z(1:3) = [ZR,ZR,ONE]
               ! Determine ThetaX, ThetaY, ThetaZ so that all three sets of transformations above are satisfied.
               call Tool_ThetaX_ThetaY_ThetaZ_3D_rotation(1,                            &
                      Solid_El_Tip_BaseLine_x_Vec(i_E)%row(Solid_El_num_Crs(i_E),1:3),  &
                      Solid_El_Tip_BaseLine_y_Vec(i_E)%row(Solid_El_num_Crs(i_E),1:3),  &
                      Solid_El_Tip_BaseLine_z_Vec(i_E)%row(Solid_El_num_Crs(i_E),1:3),  &             
                      Vector_Crack_x(1:3),Vector_Crack_y(1:3),Vector_Crack_z(1:3),      &
                      ThetaX,ThetaY,ThetaZ,T_Matrx(1:3,1:3))
               Solid_El_Tip_BaseLine_T_Matrix(i_E)%row(Solid_El_num_Crs(i_E),1:3,1:3)=T_Matrx(1:3,1:3)
               Solid_El_Tip_BaseLine_T_theta(i_E)%row(Solid_El_num_Crs(i_E),1:3) =[ThetaX,ThetaY,ThetaZ]
           endif
           
           ! If crack tip enhancement is allowed
           if (Key_TipEnrich>=1) then  
              ! If the number of cracks associated with the enhancement element is less than 4,
              ! and the current crack has not yet been associated, then execute IMPROV2022122701.
              if ((Elem_num_Related_Cracks(i_E) <Key_Ele_Max_Related_Cracks) .and.&
                  (.not. any(Elem_Related_Cracks(i_E,1:Key_Ele_Max_Related_Cracks)==i_C)))then   
                   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   ! Update Elem_num_Related_Cracks(i_E) and Elem_Related_Cracks(i_E,1:4)
                   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                   Elem_num_Related_Cracks(i_E)  = Elem_num_Related_Cracks(i_E)  + 1
                   Elem_Related_Cracks(i_E,Elem_num_Related_Cracks(i_E)) = i_C
                  
                   ! Determine crack tip enhancement element
                   Elem_Type_3D(i_E,i_C)   = 1 
                   Enriched_Node_Type_3D(c_NN,i_C) = 1    
                   
                   ! Ele_Num_Tip_Enriched_Node_3D(c_NN,i_C) = i_E  ! Enhanced element number 
                   ! corresponding to the enriched node (reference element number)
                   ! Using a staggered array, IMPROV2022091804. There may be potential 
                   ! data conflicts, so a CRITICAL directive is added.
                   do i_Node =1,8
                      if(.not. allocated(Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node))%row)) then
                          allocate(Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node))%row(num_Crack))
                          Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node))%row(1:num_Crack) = 0
                      endif
                      Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node))%row(i_C) = i_E 
                   enddo        
               endif
           ! If crack tip enhancement is not allowed, no enhancement will be performed. 2022-05-07.
           ! BUGFIX2022050701.
           elseif(Key_TipEnrich==0) then 
               !Do Nothing.
           endif
           Tem_Elem_Type(i_E,i_C) = 1
         endif
       endif
     endif
  enddo
enddo

!..........................................................................
! Step 4: To prevent the stiffness matrix from becoming singular, remove 
!         (or prompt a warning for) the Heaviside enriched elements of 
!          overlapping cracks.
!         2022-08-24. BUGFIX2022082401. 
! Overlap refers to Heaviside enhancement elements that have the same 
! enhancement configuration. It is not limited to this.
! Completely overlapping; if located in the same unit, also delete 
! the former (the one corresponding to the small crack number).
! Ref: My PhiPsi Development Notebook V1-P66.
! Ref: My PhiPsi Development Notebook V1-P69.
!.......................................................................... 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,c_Elem_num_Related_Cracks,  &
!$OMP        c_NN,i_C,n_Vector_Sign_Int_1,n_Vector_Sign_Int_2,            &
!$OMP        Heaviside_Enrich_Type_Int_1,Heaviside_Enrich_Type_Int_2,     &
!$OMP        j_C,i_Node,c_Yes_Overlaped_Ele,c_Logical_Dif)                &
!$OMP SCHEDULE(static)       
do i_E = 1,Num_Elem  
    ! Number of Heaviside-enhanced cracks in the current element.
    c_Elem_num_Related_Cracks = count(Elem_Type_3D(i_E,1:num_Crack)==2)
    ! If the current element is related to two or more cracks, further assessment is required.
    if (c_Elem_num_Related_Cracks >= 2)then
        c_NN  = G_NN(:,i_E)
        do i_C=1,num_Crack  
            ! If it is a Heaviside enhancement element
            if(Elem_Type_3D(i_E,i_C) == 2)then
               !===============================================================
               ! Enhanced configuration integer code for Heaviside enhancement
               ! elements, used for near overlap or overlap deletion
               !===============================================================
               ! The integer sign (1, -1, 0) of the symbolic distance of the current 
               ! Heaviside enriched element's 8 nodes relative to crack i_C
               n_Vector_Sign_Int_1 = 0
               do i_Node = 1,8
                   !n_Vector_Sign_Int_1(i_Node) = Tool_Function_Sign_Int(Dis_Node_to_FS(c_NN(i_Node),i_C))
                   n_Vector_Sign_Int_1(i_Node) = Tool_Function_Sign_Int(Dis_Node_to_FS(c_NN(i_Node))%row(i_C))
               enddo
               ! Current Heaviside enhanced element enhancement configuration integer code
               Heaviside_Enrich_Type_Int_1 = 10000000*n_Vector_Sign_Int_1(1) +&
                                              1000000*n_Vector_Sign_Int_1(2) +&
                                               100000*n_Vector_Sign_Int_1(3) +&
                                                10000*n_Vector_Sign_Int_1(4) +&
                                                 1000*n_Vector_Sign_Int_1(5) +&
                                                  100*n_Vector_Sign_Int_1(6) +&
                                                   10*n_Vector_Sign_Int_1(7) +&
                                                    1*n_Vector_Sign_Int_1(8)    
               ! Level 2 Fracture Cycle
               do j_C = i_C,num_Crack 
                   if(j_C /= i_C) then
                       ! If the current Heaviside enhanced element is also a Heaviside enhanced element relative to j_C
                       if (Elem_Type_3D(i_E,j_C) == 2) then
                           !===============================================================
                           ! Enhanced configuration integer code for Heaviside enhancement 
                           ! elements, used for near overlap or overlap removal
                           !===============================================================
                           ! The integer sign (1, -1, 0) of the symbolic distance of the 
                           ! current Heaviside enriched element's 8 nodes relative to crack i_C
                           n_Vector_Sign_Int_2 = 0
                           do i_Node = 1,8                        
                               !n_Vector_Sign_Int_2(i_Node) = Tool_Function_Sign_Int(Dis_Node_to_FS(c_NN(i_Node),j_C))
                               n_Vector_Sign_Int_2(i_Node) = Tool_Function_Sign_Int(Dis_Node_to_FS(c_NN(i_Node))%row(j_C))
                           enddo
                           ! Current Heaviside enhanced element's enhanced configuration integer code
                           Heaviside_Enrich_Type_Int_2 = &
                                             10000000*n_Vector_Sign_Int_2(1) +&
                                              1000000*n_Vector_Sign_Int_2(2) +&
                                               100000*n_Vector_Sign_Int_2(3) +&
                                                10000*n_Vector_Sign_Int_2(4) +&
                                                 1000*n_Vector_Sign_Int_2(5) +&
                                                  100*n_Vector_Sign_Int_2(6) +&
                                                   10*n_Vector_Sign_Int_2(7) +&
                                                    1*n_Vector_Sign_Int_2(8)   
                                                   

                           ! IMPROV2023022701. Improved overlapping crack detection.
                           !IIIIIIIIIIIIIIIIIIIIIIIIIII
                           ! If the configuration is the same.
                           !IIIIIIIIIIIIIIIIIIIIIIIIIII
                           if(Heaviside_Enrich_Type_Int_2==  Heaviside_Enrich_Type_Int_1) then
                               ! Further calculation.
                               c_Yes_Overlaped_Ele = .True.
                               do i_Node = 1,8
                                   call Tool_Dif_Ratio_of_Two_Values(Dis_Node_to_FS(c_NN(i_Node))%row(i_C),&
                                                                     Dis_Node_to_FS(c_NN(i_Node))%row(j_C),1,&
                                                                     Elem_Ave_L(i_E)*0.05D0,&
                                                                     c_Logical_Dif)
                                   if(c_Logical_Dif.eqv..False.) then
                                       c_Yes_Overlaped_Ele = .False.
                                       exit
                                   endif
                               enddo
                               if(c_Yes_Overlaped_Ele .eqv. .True.) then
                                   !/////////////////////////
                                   ! Terminal output prompt.
                                   !/////////////////////////
                                   write(*,1901) i_E,j_C,i_C   
                                   !$OMP CRITICAL
                                   Elem_Type_3D(i_E,j_C) = 0  
                                   Enriched_Node_Type_3D(c_NN(1:8),j_C) = 0
                                   !$OMP END CRITICAL
                               endif
                           endif
                           !IIIIIIIIIIIIIIIIIIIIIIIIIII
                           ! If the configuration is opposite.
                           !IIIIIIIIIIIIIIIIIIIIIIIIIII
                           if(Heaviside_Enrich_Type_Int_2== -Heaviside_Enrich_Type_Int_1) then
                               ! Further calculation.
                               c_Yes_Overlaped_Ele = .True.
                               do i_Node = 1,8
                                   !c_Threshold = Tool_Function_2Point_Dis_3D(Coor(Num_Node,1:3),B)
                                   call Tool_Dif_Ratio_of_Two_Values(Dis_Node_to_FS(c_NN(i_Node))%row(i_C),&
                                                                     Dis_Node_to_FS(c_NN(i_Node))%row(j_C),1,&
                                                                     Elem_Ave_L(i_E)*0.05D0,&
                                                                     c_Logical_Dif)
                                   if(c_Logical_Dif.eqv..False.) then
                                       c_Yes_Overlaped_Ele = .False.
                                       exit
                                   endif
                               enddo
                               if(c_Yes_Overlaped_Ele .eqv. .True.) then
                                   !/////////////////////////
                                   ! Terminal output prompt.
                                   !/////////////////////////
                                   write(*,1901) i_E,j_C,i_C   
                                   !$OMP CRITICAL
                                   Elem_Type_3D(i_E,j_C) = 0
                                   Enriched_Node_Type_3D(c_NN(1:8),j_C) = 0
                                   !$OMP END CRITICAL
                               endif
                           endif
                       endif
                   endif
               enddo
            endif
        enddo
    endif
enddo
!$OMP END PARALLEL DO  



!.............................................................
!                                                           .
! Step 5: All nodes of the Heaviside enhanced elements are 
!         Heaviside enhanced nodes.
!                                                           .
!.............................................................
! IMPROV2022072205. OpenMP parallelization. Extensive testing required.
!$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_C,i_E,i_Node) &
!$OMP SCHEDULE(static)  
!do i_C=1,num_Crack   
!  do i_E = 1,Num_Elem     
  
!BUGFIX2023021901.   
do i_E = 1,Num_Elem 
  do i_C=1,num_Crack
     if(Elem_Type_3D(i_E,i_C) == 2)then
       do i_Node=1,8
          ! If the node has not been enhanced
          if (Enriched_Node_Type_3D(G_NN(i_Node,i_E),i_C) ==0)then
              Enriched_Node_Type_3D(G_NN(i_Node,i_E),i_C) = 2 
              
              ! BUGFIX2022091801. Due to potential data conflicts, a CRITICAL directive has been added here.
              !$OMP CRITICAL   
              
              ! Enriched_Node_Crack_n_Vector_3D: The outward normal vector of the 
              ! crack surface corresponding to the enriched node
              ! (used for penalty function treatment of compression-shear cracks).
              ! Using the staggered array. 2022-09-18. IMPROV2022091801.
              if(.not. allocated(Enriched_Node_Crack_n_Vector_3D(G_NN(i_Node,i_E))%row))then  
                  allocate(Enriched_Node_Crack_n_Vector_3D(G_NN(i_Node,i_E))%row(num_Crack,3))
                  Enriched_Node_Crack_n_Vector_3D(G_NN(i_Node,i_E))%row(1:num_Crack,1:3) = ZR
              endif
                  
              
              ! Major bug fix, before the bug fix. BUGFIX2022102102.
              !if(allocated(Ragged_Array_2D_n_Vectors(i_Node)%row)) then
              !    Enriched_Node_Crack_n_Vector_3D(G_NN(i_Node,i_E))%row(i_C,1:3) =
              !    Ragged_Array_2D_n_Vectors(i_Node)%row(i_C,1:3)
              !endif
              
              ! Major bug fix. BUGFIX2022102102.
              if(allocated(Ragged_Array_2D_n_Vectors(G_NN(i_Node,i_E))%row)) then 
                  Enriched_Node_Crack_n_Vector_3D(G_NN(i_Node,i_E))%row(i_C,1:3) = &
                                         Ragged_Array_2D_n_Vectors(G_NN(i_Node,i_E))%row(i_C,1:3) 
              endif
              
              !$OMP END CRITICAL
          endif
       enddo
     endif
  enddo
enddo
!$omp end parallel do  

!.................................................................
!                                                             .
! Step 5.1: Calculate the Enriched_Node_Crack_n_Vector_3D for 
!           the crack tip enrichment element.
! Used for penalty function handling of compression-shear cracks.
!          2023-08-24. IMPROV2023082401.                      .
!.................................................................
! Enriched_Node_Crack_n_Vector_3D: The outward normal vector of 
! the crack surface corresponding to the enriched node 
! (used for the penalty function treatment of compression-shear cracks).
if(Key_TipEnrich/=0) then
    !$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_C,i_E,i_Node) &
    !$OMP SCHEDULE(static)   
    do i_E = 1,Num_Elem 
        do i_C=1,num_Crack
             if(Elem_Type_3D(i_E,i_C) == 1)then
               do i_Node=1,8
                  ! BUGFIX2022091801. Due to potential data conflicts, a CRITICAL directive has been added here.
                  !$OMP CRITICAL   
                  ! Using the staggered array. 2022-09-18. IMPROV2022091801.
                  if(.not. allocated(Enriched_Node_Crack_n_Vector_3D(G_NN(i_Node,i_E))%row))then  
                      allocate(Enriched_Node_Crack_n_Vector_3D(G_NN(i_Node,i_E))%row(num_Crack,3))
                      Enriched_Node_Crack_n_Vector_3D(G_NN(i_Node,i_E))%row(1:num_Crack,1:3) = ZR
                  endif
                  if(allocated(Ragged_Array_2D_n_Vectors(G_NN(i_Node,i_E))%row)) then 
                      Enriched_Node_Crack_n_Vector_3D(G_NN(i_Node,i_E))%row(i_C,1:3) = &
                                             Ragged_Array_2D_n_Vectors(G_NN(i_Node,i_E))%row(i_C,1:3) 
                  endif
                  !$OMP END CRITICAL
               enddo
             endif
        enddo 
    enddo   
    !$omp end parallel do  
endif

!....................................................................
!                                                                  .
!                                                                  .
! Step 6: To prevent the stiffness matrix from being singular, 
!         remove the Heaviside-enhanced nodes that do not meet
!         the criteria.
!         2022-07-06. IMPROV2022070601.                            .
!                                                                  .
!                                                                  .
!....................................................................
! Depending on the points system:
if(allocated(kesi_Enr)) deallocate(kesi_Enr)
if(allocated(yita_Enr)) deallocate(yita_Enr)
if(allocated(zeta_Enr)) deallocate(zeta_Enr)
if(allocated(weight_Enr)) deallocate(weight_Enr)      
select case(Key_Integral_Sol)
! Fixed points
case(2)
  allocate(kesi_Enr(Num_Gau_Points_3D))
  allocate(yita_Enr(Num_Gau_Points_3D))
  allocate(zeta_Enr(Num_Gau_Points_3D))
  allocate(weight_Enr(Num_Gau_Points_3D))
  call Cal_Gauss_Points_3D_8nodes(Num_Gau_Points_3D,kesi_Enr,yita_Enr,zeta_Enr,weight_Enr)  
  Tem_Num_Gauss_3D = Num_Gau_Points_3D
! Fixed block integration. NEWFTU2022072901.
case(3)
  Num_Gauss_Cubes = Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube
  allocate(kesi_Enr(Num_Gauss_Cubes))
  allocate(yita_Enr(Num_Gauss_Cubes))
  allocate(zeta_Enr(Num_Gauss_Cubes))
  allocate(weight_Enr(Num_Gauss_Cubes))
  call Cal_Gauss_Points_3D_for_SUBCUBES(Num_Sub_3D_Cubes,Num_Gau_Points_3D_Cube,kesi_Enr,yita_Enr,zeta_Enr,weight_Enr)   
  Tem_Num_Gauss_3D = Num_Gauss_Cubes
! Integration by Parts (Temporarily Abandoned). 2022-07-27.
case(4)
  goto 1120
endselect

! First, obtain the list of nodes to be checked.
! NEWFTU2022093002. OpenMP parallelization. Test passed.
Nodes_to_be_Checked = 0
c_Count_Nodes = 0   
!$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_C,i_Node)   &     
!$OMP           SCHEDULE(static)        
do i_Node = 1,Num_Node
  do i_C=1,num_Crack
      if(Enriched_Node_Type_3D(i_Node,i_C) == 2) then
          !$OMP CRITICAL 
          c_Count_Nodes = c_Count_Nodes +1
          Nodes_to_be_Checked(c_Count_Nodes) = i_Node
          !$OMP END CRITICAL 
          exit
      endif
  enddo
enddo
!$omp end parallel do  
      

! For in-plane extension, calculate the boundary line area of each crack. 
! Used for OPTION 4 (i.e., NEWFTU202302202) in subsequent loops. NEWFTU202302202.
if(Key_InPlane_Growth== 1) then
    if (allocated(Cracks_Outline_Area)) then
        deallocate(Cracks_Outline_Area)
    endif
    allocate(Cracks_Outline_Area(num_Crack))
    Cracks_Outline_Area(1:num_Crack) = ZR
    !$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_C,c_num_Points) &
    !$OMP SCHEDULE(static)   
    do i_C=1,num_Crack
        c_num_Points = Crack3D_Meshed_Outline_num(i_C)   
        call Tool_Area_3D_Plane_Polygon(c_num_Points,&
             Crack3D_Meshed_Node(i_C)%row(Crack3D_Meshed_Outline(i_C)%row(1:c_num_Points,1),1:3),&
             Cracks_Outline_Area(i_C))
    enddo
    !$omp end parallel do 
endif

! Node loop.
! Core code.
! Note, this part is very time-consuming for out-of-plane extensions. 2023-07-17.
!$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_C,i_Node,c_Num_Eles,i_E,c_Ele,c_NN,c_X_NODES,c_Y_NODES,c_Z_NODES,i_G,&
!$OMP              N,Global_coor_Gauss,c_Distance,c_Node,Check_Ball_R,c_Yes_Node_PER_in_FS,                      &
!$OMP              c_PER_Node_to_FS,Yes_Found_Min_Signed_Dis,c_n_Vector,c_Sign,c_Signed_Dis_v2,                  &
!$OMP              c_Logical_Positive,c_Logical_Negative,c_Logical_Positive_Num,c_Logical_Negative_Num      )    &
!$OMP        SCHEDULE(static)   
do i_Node = 1,c_Count_Nodes
    c_Node = Nodes_to_be_Checked(i_Node)
    
    !c_Logical_Positive = .False.
    !c_Logical_Negative = .False.
    
    do i_C=1,num_Crack
        c_Logical_Positive_Num = 0
        c_Logical_Negative_Num = 0
        c_Logical_Positive = .False.
        c_Logical_Negative = .False.
        !Logical_Positive_XFEM = .False. 
        !Logical_Negative_XFEM = .False. 
        !Logical_XFEM = .False. 
        
        ! If relative to i_C, it is a Heaviside-enhanced node.
        if(Enriched_Node_Type_3D(c_Node,i_C) == 2) then
              ! The number of elements surrounding the current node.
              c_Num_Eles = num_Node_Elements(c_Node) 
              ! Cell cycle around the current node.
              allocate(c_Distance(c_Num_Eles,Tem_Num_Gauss_3D))
              c_Distance(1:c_Num_Eles,1:Tem_Num_Gauss_3D) = ZR
              do i_E = 1,c_Num_Eles
              
                !if(Elem_Type_3D(i_E,i_C)/=2) exit !2022-09-26.
                
                !c_Ele = Node_Elements(c_Node,i_E)
                c_Ele = Node_Elements_3D(c_Node)%row(i_E)
                
                ! If it is a potential split-tip enhancement element (and not a Junction enhancement element), skip
                ! it. BUGFIX2022092701.
                !if(Tem_Elem_Type(c_Ele,i_C)==1 .and. Elem_Type_3D(c_Ele,i_C)/=4) exit
                
                c_NN(1:8)    = G_NN(1:8,c_Ele)
                c_X_NODES = G_X_NODES(1:8,c_Ele)
                c_Y_NODES = G_Y_NODES(1:8,c_Ele)  
                c_Z_NODES = G_Z_NODES(1:8,c_Ele)   
                ! Gauss point loop
                !do i_G = 1,Num_Gau_Points_3D
                do i_G = 1,Tem_Num_Gauss_3D
                  call Cal_N_3D(kesi_Enr(i_G),yita_Enr(i_G),zeta_Enr(i_G),N)
                  !Global coordinates of the gauss point.
                  Global_coor_Gauss(1) = DOT_PRODUCT(N(1,1:24:3),c_X_NODES(1:8))
                  Global_coor_Gauss(2) = DOT_PRODUCT(N(1,1:24:3),c_Y_NODES(1:8))      
                  Global_coor_Gauss(3) = DOT_PRODUCT(N(1,1:24:3),c_Z_NODES(1:8))      
                  !****************************************************************************
                  ! Calculate the signed distance from the current Gauss point to the fracture
                  ! plane.
                  !****************************************************************************
                  ! Check_Ball_R = 3.0D0*Ave_Elem_L !Only calculates the discrete triangular surfaces of cracks within
                  ! the detection IMPROV2022050501.
                  Check_Ball_R  = THR*Elem_Max_L(c_Ele)
                   
                  !call D3_Get_Signed_Dis_to_Crack_Mesh(Global_coor_Gauss,i_C,Check_Ball_R,c_Distance(i_E,i_G),   &
                  !                  c_Signed_Dis_v2,                                                             &
                  !                  c_Yes_Node_PER_in_FS,c_PER_Node_to_FS(1:3),Yes_Found_Min_Signed_Dis,c_n_Vector)
                  ! If out-of-plane extension.
                  if(Key_InPlane_Growth== 0 ) then
                      ! Resource consumption is enormous. 2023-07-17. Hot code.
                      call D3_Get_Signed_Dis_to_Crack_Mesh(Global_coor_Gauss,i_C,Check_Ball_R,c_Distance(i_E,i_G),   &
                                        c_Signed_Dis_v2,c_Yes_Node_PER_in_FS,c_PER_Node_to_FS(1:3),&
                                        Yes_Found_Min_Signed_Dis,c_n_Vector)    
                  ! For in-plane extension, there is no need to find the minimum symbol distance, 
                  ! resulting in low resource consumption. 2022-09-08. NEWFTU2022090801.
                  elseif(Key_InPlane_Growth== 1) then
                      !-------------------------------------------------------------------
                      ! OPTION 1. The most stable, but the calculation speed is not high.
                      !-------------------------------------------------------------------
                      !HOTCOD2023021801.
                      call D3_Get_Signed_Dis_to_Crack_Mesh_for_InPlane_Growth(Global_coor_Gauss,i_C,&
                                               Check_Ball_R,c_Distance(i_E,i_G),       &
                                               c_Yes_Node_PER_in_FS,c_PER_Node_to_FS(1:3),    &
                                               Yes_Found_Min_Signed_Dis,c_n_Vector)
                                        
                  endif
                  
                  call Cal_Sign(c_Distance(i_E,i_G),c_Sign)
                  
                  !///////////////////////////////////////////////////////////
                  !/////////////////         OPTION -2      /////////////////
                  !///////////////////////////////////////////////////////////
                  ! 2022-09-26. BUGFIX2022092601. It should be ensured that elements 
                  ! containing both positive and negative Gauss points are enhanced elements.
                  if(Elem_Type_3D(c_Ele,i_C)>=1) then
                      if(c_Sign>=  (ONE -Tol_6)) then
                          c_Logical_Positive = .True.
                          c_Logical_Positive_Num = c_Logical_Positive_Num + 1
                      endif
                      if(c_Sign<= (-ONE +Tol_6)) then 
                          c_Logical_Negative = .True.    
                          c_Logical_Negative_Num = c_Logical_Negative_Num + 1
                      endif
                  endif
                  
                  !BUGFIX2022092602.
                  if(c_Logical_Positive .and. c_Logical_Negative .and. &
                        c_Logical_Positive_Num>=Must_Gauss_Number_3D .and. &
                        c_Logical_Negative_Num>=Must_Gauss_Number_3D) then
                      exit
                  endif
                enddo
                ! If both 1 and -1 signs exist at the same time, exit the loop.
                !if(c_Logical_Positive .and. c_Logical_Negative) then
                
                !BUGFIX2022092602.
                if(c_Logical_Positive .and. c_Logical_Negative .and. &
                   c_Logical_Positive_Num>=Must_Gauss_Number_3D .and. &
                   c_Logical_Negative_Num>=Must_Gauss_Number_3D) then 
                   
                      exit
                endif                
              enddo
              
              deallocate(c_Distance)
              ! If both 1 and -1 signs exist at the same time, exit the loop.
              !if(c_Logical_Positive .and. c_Logical_Negative) then
              
              !BUGFIX2022092602.
              if(c_Logical_Positive .and. c_Logical_Negative .and. &
                 c_Logical_Positive_Num>=Must_Gauss_Number_3D .and. &
                 c_Logical_Negative_Num>=Must_Gauss_Number_3D) then
              else
                  !//////////////////////////////
                  ! Terminal output improvement.
                  !//////////////////////////////
                  !write(*,903) c_Node,i_C    
                  
                  Enriched_Node_Type_3D(c_Node,i_C) = 0
              endif              
              !if(Cal_Sign(Variable,Sign_O))
        endif
    enddo
enddo
!$omp end parallel do       

      
1120 CONTINUE
 
!....................................................................
!                                                                  .
!                                                                  .
! Step 7: Multiple tip fractures enhancement, 2020-02-11.
!                                                                  .
!                                                                  .
!....................................................................
!call Cal_Ele_Num_by_Coors_3D(0.75D0,0.75D0,0.48D0,c_OUT_Elem)
if(Key_Multi_TipEnrNode==1)then
    do i_C=1,num_Crack
    added_tip_el(1:Num_Elem) = .False.
    added_tip_node(1:Num_Elem,1:Num_Node) = .False.
    added_tip_node_refe(1:Num_Elem,1:Num_Node) =0        
    !******************
    !   Step1: finding
    !******************
    do i_E = 1,Num_Elem  
         c_NN  = G_NN(:,i_E)
         ! If the element has not been enhanced
         if(Elem_Type_3D(i_E,i_C) == 0)then
           ! If the element contains four or more crack tip enrichment nodes and does not contain Heaviside
           ! enrichment nodes
           num_tip_enrich=count(Enriched_Node_Type_3D(c_NN(1:8),i_C)==1) 
           num_H_enrich  =count(Enriched_Node_Type_3D(c_NN(1:8),i_C)==2)      
           if(num_tip_enrich>=4 .and. num_H_enrich <=0)then
             !~~~~~~~~~~~~~~~~~~
             !Get the ref_elem.
             !~~~~~~~~~~~~~~~~~~
             do i_Node=1,8
              if(Enriched_Node_Type_3D(c_NN(i_Node),i_C)==1)then
                !ref_elem=Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node),i_C)
                ref_elem=Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node))%row(i_C)
                exit
              endif
             enddo
             !~~~~~~~~~~~~~~~~~~~~~~~~~~
             !Save temporary variables.
             !~~~~~~~~~~~~~~~~~~~~~~~~~~
             
             added_tip_el(i_E)    = .True.
             do i_Node=1,8
               ! If the node has not been enhanced
               if(Enriched_Node_Type_3D(c_NN(i_Node),i_C) ==0)then
                added_tip_node(i_E,i_Node) = .True.
                added_tip_node_refe(i_E,i_Node) = ref_elem
               endif
             enddo
           endif
         endif
    enddo
    !*******************
    !   Step2: updating
    !*******************
    do i_E = 1,Num_Elem  
         c_NN  = G_NN(:,i_E)
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ! Bug fix: Removed the tear tip enhancement element marker
         ! You only need to obtain the additional crack tip enhancement nodes, 
         ! without needing to mark them as crack tip enhanced elements.
         ! The reasons are as follows:
         ! If it is marked as a crack tip enhanced element, then the ref_elem
         ! number obtained through the following statement is incorrect.
         !     if ((Elem_Type_3D(in_Elem,i_C).eq.1 )) then
         !         ref_elem=in_Elem
         !     else
         ! !Number of the enriched element corresponding to the enriched node (reference element number)
         ! ref_elem = Ele_Num_Tip_Enriched_Node_3D(NODES_iE(i_N), i_C)  
         !     endif
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         do i_Node=1,8
             if(added_tip_node(i_E,i_Node) .eqv. .True.) then
               Enriched_Node_Type_3D(c_NN(i_Node),i_C) =1
               ref_elem = added_tip_node_refe(i_E,i_Node) 
               
               ! Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node),i_C) = ref_elem ! Reference element number for the crack
               ! tip enriched node
                if(.not. allocated(Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node))%row)) then
                    allocate(Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node))%row(num_Crack))
                    Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node))%row(1:num_Crack) = 0
                endif
                ! The outward normal vector of the crack surface corresponding to the enhanced node (used for
                ! the penalty function treatment of compression-shear cracks).
                Ele_Num_Tip_Enriched_Node_3D(c_NN(i_Node))%row(i_C) =ref_elem
                Enriched_Node_Crack_n_Vector_3D(G_NN(i_Node,ref_elem))%row(i_C,1:3)= &
                          Ragged_Array_2D_n_Vectors(G_NN(i_Node,ref_elem))%row(i_C,1:3) 
!                Enriched_Node_Crack_n_Vector_3D(G_NN(i_Node,i_E))%row(i_C,1:3) = &
!                                         Ragged_Array_2D_n_Vectors(G_NN(i_Node,i_E))%row(i_C,1:3) 
             endif
         enddo
    enddo
    enddo
endif



!.....................................................................
!                                                                   . 
!                                                                   .
! Step 8: Identify the Junction enhancement node, NEWFTU2022050301. .
!                                                                   .
!                                                                   .
!.....................................................................
!!###################################################
! OPTION 1 - Determined by the intersection points of 
! discrete triangles on the fracture surface
!!###################################################
! IMPROV2022072205. OpenMP Parallelization.
! If you enable Junction enhancement.
if (Key_Junction_Enrich==1) then
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C,j_C,i_Crack_Ele,i_Cr_Node1,i_Cr_Node2,i_Cr_Node3,c_Tri_1,  &
  !$OMP              j_Crack_Ele,j_Cr_Node1,j_Cr_Node2,j_Cr_Node3,c_Tri_2,c_Yes_Side_on,c_Inter_Points,    &
  !$OMP              c_Side_Tri_Number,c_Elem,c_NN,delta_L,new_A,c_Tri_Edge_Num,                           &
  !$OMP              Check_Mesh_Node_1,Check_Mesh_Node_2,Mesh_Node_1_on_Outline,Mesh_Node_2_on_Outline,    &
  !$OMP              c_num_Outlines,c_location,i_Node,Ele_Num_Cache)                                       &
  !$OMP SCHEDULE(static)       
  do i_C=1,num_Crack
    Ele_Num_Cache = 1
    !do j_C=1,num_Crack
    do j_C=i_C,num_Crack
    ! If there is an overlap between the coordinate ranges of the two cracks.
    if(Crack_Coor_Overlap_Status(i_C,j_C) .eqv. .True.) then
        ! num_i_j_Inter = 0 !The total number of intersection points of these two cracks
        ! i_C Discrete Triangular Loop of Cracks
        do i_Crack_Ele =1,Crack3D_Meshed_Ele_num(i_C)
          i_Cr_Node1 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,1)
          i_Cr_Node2 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,2)
          i_Cr_Node3 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,3)
          c_Tri_1(1,1:3) = Crack3D_Meshed_Node(i_C)%row(i_Cr_Node1,1:3)
          c_Tri_1(2,1:3) = Crack3D_Meshed_Node(i_C)%row(i_Cr_Node2,1:3)
          c_Tri_1(3,1:3) = Crack3D_Meshed_Node(i_C)%row(i_Cr_Node3,1:3)
          ! Discrete triangular loop of j_C cracks
          do j_Crack_Ele =1,Crack3D_Meshed_Ele_num(j_C)
            j_Cr_Node1 = Crack3D_Meshed_Ele(j_C)%row(j_Crack_Ele,1)
            j_Cr_Node2 = Crack3D_Meshed_Ele(j_C)%row(j_Crack_Ele,2)
            j_Cr_Node3 = Crack3D_Meshed_Ele(j_C)%row(j_Crack_Ele,3)
            c_Tri_2(1,1:3) =Crack3D_Meshed_Node(j_C)%row(j_Cr_Node1,1:3)
            c_Tri_2(2,1:3) =Crack3D_Meshed_Node(j_C)%row(j_Cr_Node2,1:3)
            c_Tri_2(3,1:3) =Crack3D_Meshed_Node(j_C)%row(j_Cr_Node3,1:3)
            
            ! The following code segment has a bug, which BUGFIX2022053001 has fixed.
            ! Check whether the discrete fracture surface triangles of the two cracks intersect.
            ! call
            ! Tool_Intersections_of_Two_Triangles_3D(c_Tri_1,c_Tri_2,c_Logical_Inter,c_num_Inters,c_Inter_Points,c_Logical_Parallel)
            ! If there is an intersection.
            !if(c_Logical_Inter .eqv. .True.)then           
            
            !BUGFIX2022053001.
            ! Check the discrete crack surface triangles of the two cracks: if an edge of one triangle lies on
            ! another triangle.
            call Tool_Yes_Triangle_Egde_on_3D_Triangle(c_Tri_1,c_Tri_2,c_Yes_Side_on,c_Inter_Points(1,1:3), &
                      c_Side_Tri_Number,c_Tri_Edge_Num)    
            ! Check whether the two triangles are completely parallel (with the outward normal directions
            ! aligned). 2022-08-05. BUGFIX2022080501.
            call Tool_Yes_Two_Triangles_Parallel(c_Tri_1,c_Tri_2,c_Yes_Parallel)    
            
            ! If one of the sides of a triangle lies on another triangle.
            !if(c_Yes_Side_on .eqv. .True.)then  
            
            ! If one side of a triangle lies on another triangle, and the two triangles are not parallel.
            ! 2022-08-05. BUGFIX2022080501.
            if((c_Yes_Side_on .eqv. .True.) .and. (c_Yes_Parallel .eqv. .False.))then
                !////////////////////////////////////////////////////////////////////////////////////
                ! Check whether the discrete points of the two fracture surfaces corresponding to 
                ! c_Tri_Edge_Num are on the boundary line Outlines of the discrete fracture
                ! surfaces.
                ! A detection condition has been added: the intersection line between an active
                ! crack
                ! and a passive crack must lie on the boundary line of the discrete crack surface of
                ! the active crack.
                ! 2022-08-25. BUGFIX2022082501.
                !////////////////////////////////////////////////////////////////////////////////////
                if(c_Side_Tri_Number==1)then
                    if (c_Tri_Edge_Num ==1) then
                        Check_Mesh_Node_1 = i_Cr_Node1
                        Check_Mesh_Node_2 = i_Cr_Node2
                    elseif(c_Tri_Edge_Num ==2) then
                        Check_Mesh_Node_1 = i_Cr_Node2
                        Check_Mesh_Node_2 = i_Cr_Node3
                    elseif(c_Tri_Edge_Num ==3) then
                        Check_Mesh_Node_1 = i_Cr_Node3
                        Check_Mesh_Node_2 = i_Cr_Node1
                    else
                        print *,'    Unexpected case in Determine_Enriched_Nodes_3D.f90!'
                    endif
                    c_num_Outlines = Crack3D_Meshed_Outline_num(i_C)
                    ! Check whether Check_Mesh_Node_1 is in crack i_C (active crack)
                    call Vector_Location_Int(c_num_Outlines,Crack3D_Meshed_Outline(i_C)%row(1:c_num_Outlines,1),&
                                             Check_Mesh_Node_1,c_location,Mesh_Node_1_on_Outline)   
                    call Vector_Location_Int(c_num_Outlines,Crack3D_Meshed_Outline(i_C)%row(1:c_num_Outlines,1),&
                                             Check_Mesh_Node_2,c_location,Mesh_Node_2_on_Outline)      
                elseif(c_Side_Tri_Number==2)then
                    if (c_Tri_Edge_Num ==1) then
                        Check_Mesh_Node_1 = j_Cr_Node1
                        Check_Mesh_Node_2 = j_Cr_Node2
                    elseif(c_Tri_Edge_Num ==2) then
                        Check_Mesh_Node_1 = j_Cr_Node2
                        Check_Mesh_Node_2 = j_Cr_Node3
                    elseif(c_Tri_Edge_Num ==3) then
                        Check_Mesh_Node_1 = j_Cr_Node3
                        Check_Mesh_Node_2 = j_Cr_Node1
                    else
                        print *,'    Unexpected case in Determine_Enriched_Nodes_3D.f90!'
                    endif
                    c_num_Outlines = Crack3D_Meshed_Outline_num(j_C)
                    ! Check whether Check_Mesh_Node_1 is in the crack j_C (active crack)
                    call Vector_Location_Int(c_num_Outlines,Crack3D_Meshed_Outline(j_C)%row(1:c_num_Outlines,1),&
                                             Check_Mesh_Node_1,c_location,Mesh_Node_1_on_Outline)   
                    call Vector_Location_Int(c_num_Outlines,Crack3D_Meshed_Outline(j_C)%row(1:c_num_Outlines,1),&
                                             Check_Mesh_Node_2,c_location,Mesh_Node_2_on_Outline)                                            
                endif
                
                !//////////////////////////////////////////////////////////////////////////////////
                ! If both points are on the boundary line outlines of the discrete fracture plane, 
                ! they are marked as Junction enhancement elements.
                ! 2022-08-25. BUGFIX2022082501.
                !//////////////////////////////////////////////////////////////////////////////////
                if(Mesh_Node_1_on_Outline .and. Mesh_Node_2_on_Outline) then
                  ! Get the cell number of the intersection
                  call Cal_Ele_Num_by_Coors_3D(c_Inter_Points(1,1),c_Inter_Points(1,2),c_Inter_Points(1,3),Ele_Num_Cache,c_Elem) 
                  ! Obtain active and passive crack numbers:
                  ! If the element is a crack-tip enriched element relative to i_C and a Heaviside enriched
                  ! element relative to j_C, then the element is a Junction enriched element.
                  if(Tem_Elem_Type(c_Elem,i_C)==1 .and. Tem_Elem_Type(c_Elem,j_C)==2 .and. &
                     c_Side_Tri_Number==1 .and. Elem_Type_3D(c_Elem,i_C)/=4) then
                    Elem_Type_3D(c_Elem,i_C) = 4
                    c_NN  = G_NN(:,c_Elem)
                    Enriched_Node_Type_3D(c_NN,i_C) = 3
                    do i_Node =1,8
                      if(allocated(Ragged_Array_2D_n_Vectors(c_NN(i_Node))%row)) then
                          !Enriched_Node_Crack_n_Vector_3D(c_NN(i_Node),i_C,1:3) = &
                          !                   Ragged_Array_2D_n_Vectors(c_NN(i_Node))%row(i_C,1:3) 
                          
                          ! Using the staggered array. 2022-09-18. IMPROV2022091801.
                          if(.not. allocated(Enriched_Node_Crack_n_Vector_3D(c_NN(i_Node))%row))then  
                              allocate(Enriched_Node_Crack_n_Vector_3D(c_NN(i_Node))%row(num_Crack,3))
                              Enriched_Node_Crack_n_Vector_3D(c_NN(i_Node))%row(1:num_Crack,1:3) = ZR
                          endif
                          Enriched_Node_Crack_n_Vector_3D(c_NN(i_Node))%row(i_C,1:3) = &
                                             Ragged_Array_2D_n_Vectors(c_NN(i_Node))%row(i_C,1:3)     
                      endif
                    enddo
                                                                        ! Note: For Junction enhancement elements, the n_Vectors of some nodes may be 0, 2022-05-12.
                    ! Jun_Ele_Negative_Cr_Num_3D(c_Elem, i_C) = j_C !Corresponding passive crack number (main crack
                    ! number)
                    if(.not. allocated(Jun_Ele_Negative_Cr_Num_3D(c_Elem)%row))then  
                        allocate(Jun_Ele_Negative_Cr_Num_3D(c_Elem)%row(num_Crack))
                        Jun_Ele_Negative_Cr_Num_3D(c_Elem)%row(1:num_Crack) = 0
                    endif
                    Jun_Ele_Negative_Cr_Num_3D(c_Elem)%row(i_C)         = j_C 
                    
                    ! Node_Jun_elem_3D(c_NN,i_C) = c_Elem !Junction enhanced node corresponding to the Junction element
                    ! number
                    do i_Node =1,8
                      if(.not. allocated(Node_Jun_elem_3D(c_NN(i_Node))%row)) then
                          allocate(Node_Jun_elem_3D(c_NN(i_Node))%row(num_Crack))
                          Node_Jun_elem_3D(c_NN(i_Node))%row(1:num_Crack) = 0
                      endif
                      Node_Jun_elem_3D(c_NN(i_Node))%row(i_C)  = c_Elem                      
                    enddo
                    
                    ! Intersection of active cracks (secondary cracks) and element boundaries (does not necessarily have
                    ! to be an intersection;
                    ! the purpose is to determine on which side of the main crack (passive crack) the secondary crack is
                    ! located)
                    ! The approach used here is to offset a certain distance from intersection 1 toward the centroid of
                    ! the secondary fracture.
                    !delta_L = Ave_Elem_L/TWO
                    delta_L = Elem_Max_L(c_Elem)/TWO
                    call Tool_Offset_Point_A_to_Point_B_3D(c_Inter_Points(1,1:3),Cracks_Coor_Centers(i_C,1:3),delta_L,new_A)   
                    
                    !Coors_Junction_3D(c_Elem,i_C,1:3) = new_A
                    if(.not. allocated(Coors_Junction_3D(c_Elem)%row))then  
                        allocate(Coors_Junction_3D(c_Elem)%row(num_Crack,3))
                        Coors_Junction_3D(c_Elem)%row(1:num_Crack,1:3) = 0
                    endif
                    Coors_Junction_3D(c_Elem)%row(i_C,1:3)  = new_A
                  endif
                  ! If the element is a crack-tip enriched element relative to j_C and a Heaviside enriched element
                  ! relative to i_C,
                  ! then the element is a Junction enriched element.
                  if(Tem_Elem_Type(c_Elem,j_C)==1 .and. Tem_Elem_Type(c_Elem,i_C)==2 .and.  &
                     c_Side_Tri_Number==2         .and.  &
                     Elem_Type_3D(c_Elem,j_C)/=4) then           
                    Elem_Type_3D(c_Elem,j_C) = 4
                    c_NN  = G_NN(:,c_Elem)
                    Enriched_Node_Type_3D(c_NN,j_C) = 3
                    do i_Node =1,8
                      if(allocated(Ragged_Array_2D_n_Vectors(c_NN(i_Node))%row)) then
                          !Enriched_Node_Crack_n_Vector_3D(c_NN(i_Node),j_C,1:3) = &
                          !                   Ragged_Array_2D_n_Vectors(c_NN(i_Node))%row(j_C,1:3) 
                          
                          ! Using the staggered array. 2022-09-18. IMPROV2022091801.
                          if(.not. allocated(Enriched_Node_Crack_n_Vector_3D(c_NN(i_Node))%row))then  
                              allocate(Enriched_Node_Crack_n_Vector_3D(c_NN(i_Node))%row(num_Crack,3))
                              Enriched_Node_Crack_n_Vector_3D(c_NN(i_Node))%row(1:num_Crack,1:3) = ZR
                          endif
                          Enriched_Node_Crack_n_Vector_3D(c_NN(i_Node))%row(j_C,1:3) = &
                                             Ragged_Array_2D_n_Vectors(c_NN(i_Node))%row(j_C,1:3)                             
                      endif
                    enddo                                                                        
                                                                        
                                                                      
                    ! Jun_Ele_Negative_Cr_Num_3D(c_Elem,j_C) = i_C !Corresponding passive crack number (main crack
                    ! number)
                    if(.not. allocated(Jun_Ele_Negative_Cr_Num_3D(c_Elem)%row))then  
                        allocate(Jun_Ele_Negative_Cr_Num_3D(c_Elem)%row(num_Crack))
                        Jun_Ele_Negative_Cr_Num_3D(c_Elem)%row(1:num_Crack) = 0
                    endif
                    Jun_Ele_Negative_Cr_Num_3D(c_Elem)%row(j_C)         = i_C 
                    
                    ! Node_Jun_elem_3D(c_NN,j_C) = c_Elem !Junction enhanced node corresponding to the Junction element
                    ! number
                    do i_Node =1,8
                      if(.not. allocated(Node_Jun_elem_3D(c_NN(i_Node))%row)) then
                          allocate(Node_Jun_elem_3D(c_NN(i_Node))%row(num_Crack))
                          Node_Jun_elem_3D(c_NN(i_Node))%row(1:num_Crack) = 0
                      endif
                      Node_Jun_elem_3D(c_NN(i_Node))%row(j_C)  = c_Elem                      
                    enddo
                    
                    ! Intersection of active cracks (secondary cracks) and unit boundaries (does not necessarily have 
                    ! to be an intersection; the purpose is to determine on which side of 
                    ! the main crack (passive crack) the secondary crack is located)
                    ! The approach used here is to offset a certain distance from intersection 1 toward the 
                    ! centroid of the secondary fracture.
                    !delta_L = Ave_Elem_L/TWO
                    delta_L = Elem_Max_L(c_Elem)/TWO
                    call Tool_Offset_Point_A_to_Point_B_3D(c_Inter_Points(1,1:3),Cracks_Coor_Centers(j_C,1:3),delta_L,new_A)    
                    
                    !Coors_Junction_3D(c_Elem,j_C,1:3) = new_A       
                    if(.not. allocated(Coors_Junction_3D(c_Elem)%row))then  
                        allocate(Coors_Junction_3D(c_Elem)%row(num_Crack,3))
                        Coors_Junction_3D(c_Elem)%row(1:num_Crack,1:3) = 0
                    endif
                    Coors_Junction_3D(c_Elem)%row(j_C,1:3)  = new_A
                  endif
                  
                  ! For passive cracks, if the element is a crack tip enriched element, Heaviside enrichment is not
                  ! applied.
                  !Yes_Tip_El= any(Tem_Elem_Type(c_Elem,:)==1)
                endif
            endif
          enddo
        enddo
    endif
    enddo
  enddo
  !$omp end parallel do         
endif

!................................................................
!                                                      .
!                                                      .
! Step 10: Determine the average unit area of the enhanced unit.
!                                                      .
!                                                      .
!................................................................
! Cycle between elements
! NEWFTU2022093002. OpenMP Parallelization.
Vol_Enrich = ZR
Vol_Enrich_Eles(1:Num_Elem)      = ZR
Vol_Enrich_Eles_Flag(1:Num_Elem) = 0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E) SCHEDULE(static) 
do i_E = 1,Num_Elem
    if (maxval(Elem_Type_3D(i_E,1:num_crack)) /=0 .or.minval(Elem_Type_3D(i_E,1:num_crack)) /=0) then
        Vol_Enrich_Eles_Flag(i_E)  = 1
        Vol_Enrich_Eles(i_E) = Elem_Vol(i_E)
    end if
end do
!$omp end parallel do   
iEl_contour =  sum(Vol_Enrich_Eles_Flag(1:Num_Elem))   
Vol_Enrich  =  sum(Vol_Enrich_Eles(1:Num_Elem))
Ave_Elem_Vol_Enrich  = Vol_Enrich/iEl_contour
Ave_Elem_L_Enrich    = Ave_Elem_Vol_Enrich**(1.0D0/3.0D0)
Ave_Elem_Area_Enrich = Ave_Elem_L_Enrich**2
! For the local encryption process, obtain the feature length of the enhanced unit before local
! encryption (2021-08-22)
if (key_local_mesh_refine>=0  .and.  (Flag_Local_Refined.eqv..false.))then
      Ave_Elem_L_Enrich_Unlocalrefined = Ave_Elem_L_Enrich 
endif

!.........................................................
!                                                       .
! Step 11: Determine the list of FEM elements and the 
! list of XFEM enhanced elements, 2022-04-16.
! Note: The so-called enhanced unit refers to a unit 
! that contains enhanced nodes.
!                                                       .
!.........................................................
num_XFEM_Elem = 0
num_XFEM_Elem_Flag(1:Num_Elem) = 0
! IMPROV2022072205. OpenMP Parallelization.
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,c_NN) SCHEDULE(static)           
do i_E = 1,Num_Elem
    c_NN    = G_NN(:,i_E)
    if(sum(abs(Enriched_Node_Type_3D(c_NN,1:num_Crack))).gt.0)then               
        num_XFEM_Elem_Flag(i_E) = 1
    endif
enddo
!$omp end parallel do  
num_XFEM_Elem = sum(num_XFEM_Elem_Flag(1:Num_Elem))

if (allocated(XFEM_Elem_List))DEALLOCATE(XFEM_Elem_List)
if (allocated(FEM_Elem_List))DEALLOCATE(FEM_Elem_List)

! Saved Elem_Location_old. 2022-06-24.
! if (isub>=2) then !BUGFIX2023071801 before the fix.
if (isub>=2 .and. allocated(Elem_Location))then
    if (allocated(Elem_Location_old))DEALLOCATE(Elem_Location_old)
    allocate(Elem_Location_Old(num_Elem,2)) 
    Elem_Location_old = 0
    Elem_Location_old = Elem_Location
endif

allocate(XFEM_Elem_List(num_XFEM_Elem))
num_FEM_Elem = Num_Elem-num_XFEM_Elem

write(*,901) num_XFEM_Elem
write(*,902) num_FEM_Elem
write(*,1001) dble(num_XFEM_Elem)/dble(num_Elem)*dble(100)   

allocate(FEM_Elem_List(num_FEM_Elem))

! Saved Elem_XFEM_Flag_Old. 2022-06-24.
if (isub>=2)then
    !BUGFIX2023012302. 2023-01-23.
    if (allocated(Elem_XFEM_Flag_Old))then
        DEALLOCATE(Elem_XFEM_Flag_Old)
    endif
    allocate(Elem_XFEM_Flag_Old(num_Elem))
    if(allocated(Elem_XFEM_Flag)) then
        Elem_XFEM_Flag_Old = 0
        Elem_XFEM_Flag_Old = Elem_XFEM_Flag    
    else
        Elem_XFEM_Flag_Old = 0
    endif
endif      

if(allocated(Elem_XFEM_Flag))DEALLOCATE(Elem_XFEM_Flag)
allocate(Elem_XFEM_Flag(num_Elem))

if (allocated(Elem_Location)) DEALLOCATE(Elem_Location)
allocate(Elem_Location(num_Elem,2)) 
Elem_Location(1:num_Elem,1:2) = 0

XFEM_counter = 0
FEM_counter  = 0
! This does not require (and is not suitable for) OpenMP parallel computing. 
! Using OpenMP will result in the elements of the generated list being in a 
! different order each time. 2022-09-30.
do i_E = 1,Num_Elem
  c_NN    = G_NN(:,i_E)
  if(sum(abs(Enriched_Node_Type_3D(c_NN,1:num_Crack))).gt.0)then           
      XFEM_counter  = XFEM_counter  + 1
      XFEM_Elem_List(XFEM_counter) = i_E
      Elem_XFEM_Flag(i_E) = 1
      
      !2022-10-25. IMPROV2022102502.
      Elem_Location(i_E,1)  = XFEM_counter
  else              
      FEM_counter  = FEM_counter  + 1
      FEM_Elem_List(FEM_counter) = i_E  
      Elem_XFEM_Flag(i_E) = 0
      
      !2022-10-25. IMPROV2022102502.
      Elem_Location(i_E,2)  = FEM_counter
  end if
enddo

! If the previous step used XFEM elements and the current step uses
! FEM elements, it is marked as a degraded element. IMPROV2023021501. 2023-02-15.
if (isub>=2)then
    Num_Rollbacked_FEM_Elements = count((Elem_XFEM_Flag_Old(:)-Elem_XFEM_Flag(:))==1)
    if (Num_Rollbacked_FEM_Elements >0) then
        if(allocated(Rollbacked_FEM_Elements))DEALLOCATE(Rollbacked_FEM_Elements)
        allocate(Rollbacked_FEM_Elements(Num_Rollbacked_FEM_Elements)) 
        c_count =0
        do i_E = 1,Num_Elem
            if (Elem_XFEM_Flag_Old(i_E)==1 .and. Elem_XFEM_Flag(i_E)==0) then
                c_count = c_count + 1
                Rollbacked_FEM_Elements(c_count) = i_E
            endif
        enddo 
        print *,'    WARN-2023021501 :: XFEM --> FEM rollbacked elements found!'
        print *,'                       In Determine_Enriched_Nodes_3D.f90!'
        print *,'                       Num_Rollbacked_FEM_Elements :',Num_Rollbacked_FEM_Elements
        print *,'                       First Rollbacked_FEM_Element:',Rollbacked_FEM_Elements(1)
    endif
endif

!..................................................................
!                                                                .
! Step 13: If a crack has no reinforcement nodes, terminate the
! program. 2022-09-27.
!                                                                .
!..................................................................
!IMPROV2022092701.
! The cause of this phenomenon: it is caused either by BUGFIX2022092601 
! or by deleting too many duplicate nodes.
906  FORMAT(5X,'ERROR-2022092706 :: crack ',I8,' has no enriched node!')     
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_C)  
do i_C=1,num_Crack    
    if(sum(Enriched_Node_Type_3D(1:Num_Node,i_C))==0) then
        write(*,906) i_C
        if (Key_Check_and_Adjust_Cracks_3D/=5)then
            print *, '    Set *Key_Check_and_Adjust_Cracks_3D to 5 and try again!'
        endif
        call Warning_Message('S',Keywords_Blank)
    endif
enddo  
!$omp end parallel do


!..................................................................
!                                                                .
! Step 14: Output the relevant information of 
! Elem_num_Related_Cracks(:). 2023-02-25.
!                                                                .
!..................................................................
Max_Related_Cracks = maxval(Elem_num_Related_Cracks(:))
print *,'    Max number of related cracks of elements:',Max_Related_Cracks
if (Max_Related_Cracks>=2) then
    print *,'    Coresppongding element number:',maxloc(Elem_num_Related_Cracks(:),1)
endif

!....................................
!                                  .
!Step 15:  Clear memory.           .
!                                  .
!....................................
deallocate(Ragged_Array_2D_n_Vectors)
deallocate(RaggedArray1D_Dis_Node_to_FS_v2)
IF(ALLOCATED(tem_Solid_El_Vertex_num)) DEALLOCATE(tem_Solid_El_Vertex_num)   

#endif


#ifdef Silverfrost
print *,'    ERROR :: Silverfrost compiler failed to compile Determine_Enriched_Nodes_3D.F90!'
print *,'             In Determine_Enriched_Nodes_3D.F90.'
call Warning_Message('S',Keywords_Blank)
#endif  
       
RETURN
END SUBROUTINE Determine_Enriched_Nodes_3D
