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
 
SUBROUTINE D3_Mesh_Initial_Crack
! Generate a mesh for the initial three-dimensional crack (a rectangular 
! crack surface defined by 4 points, a circular crack surface, or an 
! elliptical crack surface).
! (1) Rectangle:
!                 line 3
!              4 ------- 3
!              |         |
!        line 4|         |line 2
!              |         |
!              1 ------- 2      
!                 line 1
! (2) Circular or oval
! (3) Polygon
!
! real(kind=FT) Crack3D_Meshed_Node(Max_Num_Cr_3D, Max_N_Node_3D, 3) !3D crack node coordinates
! after discretization, each crack consists of up to 1000 points
! integer Cr3D_Meshed_Node_in_Ele_Num(Max_Num_Cr_3D, Max_N_Node_3D) ! Element number of the 3D
! fracture nodes after discretization
! real(kind=FT) Cr3D_Meshed_Node_in_Ele_Local(Max_Num_Cr_3D, Max_N_Node_3D, 3)! Local coordinates of
! the 3D crack nodes within the element after discretization
! integer Crack3D_Meshed_Node_num(Max_Num_Cr_3D)                           !Number of 3D crack nodes after discretization
! integer Crack3D_Meshed_Ele(Max_Num_Cr_3D, Max_N_Node_3D, 3) !3D crack element numbers after
! discretization, each crack composed of up to 1000 points
! integer Crack3D_Meshed_Ele_Attri(Max_Num_Cr_3D, Max_N_Node_3D, 5) !3D crack element characteristic
! parameters after discretization (perimeter, area, etc.)
! real(kind=FT) Crack3D_Meshed_Ele_Nor_Vector(Max_Num_Cr_3D, Max_N_Node_3D, 3) !Normal vectors of 3D
! crack elements after discretization
! real(kind=FT) Crack3D_Meshed_Node_Nor_Vector(Max_Num_Cr_3D,Max_N_Node_3D,3) !3D crack node outward
! normal vector after discretization
! integer Crack3D_Meshed_Ele_num(Max_Num_Cr_3D)                            !Number of 3D crack elements after discretization
! integer Crack3D_Meshed_Outline(Max_Num_Cr_3D, Max_N_Node_3D, 4) !3D crack outer boundary after
! discretization
! !Data 1 is the first point on the boundary line of the crack front
! !Data 2 is the second point on the boundary line of the crack front
! !Data 3 corresponds to the discrete fracture element number
! !Data 4 is used to mark whether the two points of this boundary line are allowed to extend,
! extending in very small steps (2021-08-20)
! integer Crack3D_Meshed_Outline_num(Max_Num_Cr_3D) !Number of 3D crack boundary lines after
! discretization
! real(kind=FT) Crack3D_Meshed_Vertex_x_Vector(Max_Num_Cr_3D, Max_N_Node_3D,3) !Local x-axis vector
! of 3D crack boundary points after discretization
! real(kind=FT) Crack3D_Meshed_Vertex_y_Vector(Max_Num_Cr_3D, Max_N_Node_3D, 3) !Local y-axis vector
! of 3D crack boundary points after discretization
! real(kind=FT) Crack3D_Meshed_Vertex_z_Vector(Max_Num_Cr_3D, Max_N_Node_3D, 3) !Local z-axis vector
! of 3D crack boundary points after discretization


!-----------------------------
! Read public variable module
!-----------------------------
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack_Common
use Global_Crack_3D
use Global_Cal_Ele_Num_by_Coors_3D  

!---------------------------
! Variable Type Declaration
!---------------------------
implicit none
integer num_Division_x,num_Division_y
integer i_C,i_Node,i_D,j_D
real(kind=FT) CR_P1(3),CR_P2(3),CR_P3(3),CR_P4(3)
real(kind=FT) Points_P1P2(200,3)
real(kind=FT) Points_P2P3(200,3)
real(kind=FT) Points_P3P4(200,3)
real(kind=FT) Points_P4P1(200,3)
real(kind=FT) Points_Edge24_Lines(200,200,3)
real(kind=FT) tem_Points(200,3)
real(kind=FT) c_P1(3),c_P2(3)
integer Num_Div_P_P1P2,Num_Div_P_P2P3,Num_Div_P_P3P4,Num_Div_P_P4P1
integer tem_Num_Div
integer c_Count
integer node_1,node_2,node_3,node_4
real(kind=FT) length_line_P1P2,length_line_P2P3
real(kind=FT) Tool_Function_2Point_Dis_3D
integer num_of_nodes,num_of_ele
integer in_Elem_num,i_Cr_Node
real(kind=FT) c_Kesi,c_Yita,c_Zeta
real(kind=FT) cir_center(3),cir_vec(3),radius
real(kind=FT) a(3),b(3)
integer i_theta,num_divison,i_Ele
real(kind=FT) theta,norm_a,norm_b,circumference
real(kind=FT),ALLOCATABLE:: c_xyz(:,:)
integer Crack_Node1,Crack_Node2,Crack_Node3
real(kind=FT) Point1(3),Point2(3),Point3(3)
real(kind=FT) radius_a,radius_b
integer i_Crack_Ele
real(kind=FT) Vector_V1(3),Vector_V2(3),Vector_Normal(3)
integer i_Crack_Node
integer c_count_Ele
real(kind=FT) c_All_Nor_Vector(3)
integer num_Cr_Nodes
real(kind=FT) centroid_x,centroid_y,centroid_z
real(kind=FT) leg
real(kind=FT) Cricle_Elem_L,CR_Center(3),Mesh_Size
integer Circle_Elem,CR_Center_El,i_P
integer num_Cr_Poi,i_Edge,num_Division
integer Ele_Num_Cache
integer c_Max_N_Node_3D  
integer max_outline_num
integer c_outline_num
real(kind=FT),ALLOCATABLE::input_points(:,:),output_points(:,:)
integer input_point_num, output_point_num
  
if (num_Crack==0) return    

#ifndef Silverfrost    

Points_Edge24_Lines(1:200,1:200,1:3) = ZR
  
1001 FORMAT(5X,'Mesh size of the initial crack ', I3,' is ',F8.4, ' m')  

!--------------------
! Each fracture loop
!--------------------
do i_C=1,num_Crack
    Ele_Num_Cache = 1
    ! If the crack has already been discretized previously, proceed to the next loop. IMPROV2022053101.
    if(Cracks_Initial_Meshed(i_C)==1)then
        cycle
    endif
    
    ! If no memory has been allocated for the discrete fracture surfaces before, reallocate it.
    ! 2022-09-03.
    call D3_Allocate_Crack_Memory(i_C,1,0)

    !##################################
    ! Case1: Rectangular Crack Surface
    !##################################
    !if(sum(abs(Crack3D_Coor(i_C,1:4,1:3)))>Tol_20)then        
    if(sum(abs(Crack3D_Coor(i_C,1:4,1:3)))>Tol_20 .and.Each_Cr_Poi_Num(i_C) == 4)then
      CR_P1 = Crack3D_Coor(i_C,1,1:3)
      CR_P2 = Crack3D_Coor(i_C,2,1:3)
      CR_P3 = Crack3D_Coor(i_C,3,1:3)
      CR_P4 = Crack3D_Coor(i_C,4,1:3)
      
      CR_Center = (CR_P1 + CR_P2 + CR_P3 + CR_P4)/FOU
      call Cal_Ele_Num_by_Coors_3D(CR_Center(1),CR_Center(2),CR_Center(3),Ele_Num_Cache,CR_Center_El)
      if(CR_Center_El==0) then
           print *,'    ERROR-2022110105 :: CR_Center_El=0 in D3_Mesh_Initial_Crack.f90!'
           call Warning_Message('S',Keywords_Blank)
      endif 
      
      
      if(CR_Center_El==0)then
          Mesh_Size = Ave_Elem_L
          print *, '    WARNING :: initial crack center outside the model!'
      else
          Mesh_Size = Elem_Vol(CR_Center_El)**(ONE/THR)
      endif
      
      
      !//////////////////
      ! Terminal output.
      !//////////////////
      length_line_P1P2= Tool_Function_2Point_Dis_3D(CR_P1,CR_P2)
      length_line_P2P3= Tool_Function_2Point_Dis_3D(CR_P2,CR_P3)
      num_Division_x = int(ceiling(length_line_P1P2/Mesh_Size))
      num_Division_y = int(ceiling(length_line_P2P3/Mesh_Size))
      ! Error message. 2022-09-05.
      if(num_Division_x>200 )then
          print *,'    ERROR-2022090501 :: num_Division_x>200, in error in D3_Mesh_Initial_Crack.f90!'
          call Warning_Message('S',Keywords_Blank) 
      endif
      if(num_Division_y>200 )then
          print *,'    ERROR-2022090501 :: num_Division_y>200, in error in D3_Mesh_Initial_Crack.f90!'
          call Warning_Message('S',Keywords_Blank) 
      endif     
      
      !...............
      ! Division line
      !...............
      !Divide line P1-P2
      call Cal_Equal_Division_Points_3D(CR_P1,CR_P2,num_Division_x,.True.,Points_P1P2,Num_Div_P_P1P2)
      !Divide line P2-P3
      call Cal_Equal_Division_Points_3D(CR_P2,CR_P3,num_Division_y,.True.,Points_P2P3,Num_Div_P_P2P3)   
      !Divide line P3-P4
      call Cal_Equal_Division_Points_3D(CR_P3,CR_P4,num_Division_x,.True.,Points_P3P4,Num_Div_P_P3P4)   
      !Divide line P4-P1
      call Cal_Equal_Division_Points_3D(CR_P4,CR_P1,num_Division_y,.True.,Points_P4P1,Num_Div_P_P4P1)   
      
      !...............................................
      ! Individual sectioning of the connection lines
      ! between Node 2 and Node 4
      !...............................................
      do i_D = 1,num_Division_y-1
          c_P1 = Points_P2P3(i_D+1,1:3)
          c_P2 = Points_P4P1(num_Division_y+1-i_D,1:3)
          call Cal_Equal_Division_Points_3D(c_P1,c_P2,num_Division_x,.True.,tem_Points,tem_Num_Div)  
          Points_Edge24_Lines(i_D,1:num_Division_x+1,1:3)=tem_Points(1:num_Division_x+1,1:3)
          
      enddo       
      !..............................
      ! Node Coordinates and Numbers
      !..............................
      c_Count = 0
      do i_D = 1,num_Division_y+1
          do j_D = 1,num_Division_x+1
              c_Count = c_Count+1
              
              c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
              if(c_Count > c_Max_N_Node_3D) then
                  ! Expand memory. NEWFTU2022110501.
                  call D3_Allocate_Crack_Memory(i_C,1,1)
                  !c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
              endif  
        
              if (i_D ==1) then
                  Crack3D_Meshed_Node(i_C)%row(c_Count,1:3) = Points_P1P2(j_D,1:3)
              elseif (i_D ==num_Division_y+1) then
                  Crack3D_Meshed_Node(i_C)%row(c_Count,1:3) = Points_P3P4(num_Division_x+2-j_D,1:3)
              else
                  Crack3D_Meshed_Node(i_C)%row(c_Count,1:3) = Points_Edge24_Lines(i_D-1,num_Division_x+2-j_D,1:3)
              endif 
          enddo
      enddo
      num_of_nodes = c_Count
      Crack3D_Meshed_Node_num(i_C)  = num_of_nodes
      ! During detection, the number of division nodes cannot exceed 5000, because
      ! Save_Files_Crack_3D(isub) can save a maximum of 5000 crack division nodes.
      if(c_Count>5000) then
          print *, "    Error :: too much nodes(>5000), error in D3_Mesh_Initial_Crack.f!" 
          call Warning_Message('S',Keywords_Blank) 
      endif
      
      !.....................
      ! element Node Number
      !.....................
      c_Count = 0
      do i_D = 1,num_Division_y
          do j_D = 1,num_Division_x
              node_1 = (i_D-1)*(num_Division_x+1)  +  j_D  
              node_2 = (i_D-1)*(num_Division_x+1)  +  j_D +1
              node_3 = i_D*(num_Division_x+1) +1 +  j_D
              node_4 = i_D*(num_Division_x+1) + j_D 
              c_Count = c_Count+1
              
              c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
              if(c_Count > c_Max_N_Node_3D) then
                  ! Expand memory. NEWFTU2022110501.
                  call D3_Allocate_Crack_Memory(i_C,1,1)
                  !c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
              endif 
              
              Crack3D_Meshed_Ele(i_C)%row(c_Count,1) =  node_1
              Crack3D_Meshed_Ele(i_C)%row(c_Count,2) =  node_2
              Crack3D_Meshed_Ele(i_C)%row(c_Count,3) =  node_4
              c_Count = c_Count+1
             
              c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
              if(c_Count > c_Max_N_Node_3D) then
                  ! Expand memory. NEWFTU2022110501.
                  call D3_Allocate_Crack_Memory(i_C,1,1)
                  !c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
              endif 
              
              Crack3D_Meshed_Ele(i_C)%row(c_Count,1) =  node_2
              Crack3D_Meshed_Ele(i_C)%row(c_Count,2) =  node_3
              Crack3D_Meshed_Ele(i_C)%row(c_Count,3) =  node_4   
          enddo
      enddo
      num_of_ele = c_Count
      Crack3D_Meshed_Ele_num(i_C)  = c_Count
      !.............................................................................................
      ! Determine the outer boundary of the fracture surface discretization elements, and store the 
      ! results in the following two global variables:
      ! integer Crack3D_Meshed_Outline(Max_Num_Cr, Max_N_Node_3D, 3) !3D crack outer boundary after
      ! discretization (data 3 corresponds to the element number)
      ! integer Crack3D_Meshed_Outline_num(Max_Num_Cr) !Number of 3D crack boundary lines after
      ! discretization
      !.............................................................................................
      call D3_Find_Crack_Mesh_Outline(i_C,num_of_ele,num_of_nodes,Crack3D_Meshed_Ele(i_C)%row(1:num_of_ele,1:3))
      !............................................................................................
      ! Determine the element number where the discrete node of the crack surface is located and
      ! store it in Cr3D_Meshed_Node_in_Ele_Num(i_C, i_Cr_Node)
      ! If the element number is 0, it indicates that the crack surface node is outside the model.
      ! Determine the local coordinates of the elements where the crack surface discrete nodes are
      ! located, and store them in Cr3D_Meshed_Node_in_Ele_Local(Max_Num_Cr,1000,3)
      !............................................................................................
      do i_Cr_Node =1,num_of_nodes
          call Cal_Ele_Num_by_Coors_3D(   &
                       Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,1),  &
                       Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,2), &
                       Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,3),Ele_Num_Cache,in_Elem_num)    
          !if(in_Elem_num==0) then 
          if(in_Elem_num==0 .and. Key_Allow_3D_Outside_Crack==0) then
               print *,'    ERROR-2022110104 :: in_Elem_num=0 in D3_Mesh_Initial_Crack.f90!'
               call Warning_Message('S',Keywords_Blank)
          endif 
          Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(i_Cr_Node)= in_Elem_num
          if (in_Elem_num>0)then
              call Cal_KesiYita_by_Coor_3D(Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,1:3),in_Elem_num,c_Kesi,c_Yita,c_Zeta)
              Cr3D_Meshed_Node_in_Ele_Local(i_C)%row(i_Cr_Node,1:3)=[c_Kesi,c_Yita,c_Zeta] 
          endif
      enddo
      
    !##################################################################
    ! Case 1.1: Polygonal crack surface. NEWFTU2022061601. 2022-06-16.
    !##################################################################
    elseif(sum(abs(Crack3D_Coor(i_C,1:4,1:3)))>Tol_20 .and.  Each_Cr_Poi_Num(i_C) >= 5)then
      num_Cr_Poi = Each_Cr_Poi_Num(i_C)
      CR_Center(1) =sum(Crack3D_Coor(i_C,1:num_Cr_Poi,1))/num_Cr_Poi
      CR_Center(2) =sum(Crack3D_Coor(i_C,1:num_Cr_Poi,2))/num_Cr_Poi
      CR_Center(3) =sum(Crack3D_Coor(i_C,1:num_Cr_Poi,3))/num_Cr_Poi
      call Cal_Ele_Num_by_Coors_3D(CR_Center(1),CR_Center(2),CR_Center(3),Ele_Num_Cache,CR_Center_El)
      if(CR_Center_El==0) then
           print *,'    ERROR-2022110103 :: CR_Center_El=0 in D3_Mesh_Initial_Crack.f90!'
           call Warning_Message('S',Keywords_Blank)
      endif
      !Mesh_Size = Ave_Elem_L/1.5
      
      if(CR_Center_El==0)then
          Mesh_Size = Ave_Elem_L
          print *, '    WARNING :: initial crack center outside the model!'
      else
          Mesh_Size = Elem_Vol(CR_Center_El)**(ONE/THR)
      endif
      
      
      write(*,1001) i_C,Mesh_Size
      
      !%%%%%%%%%%%%%%%%%%%%%%%%
      ! Loop through each side
      !%%%%%%%%%%%%%%%%%%%%%%%%
      c_Count = 0
      do i_Edge= 1,num_Cr_Poi
          if(i_Edge/=num_Cr_Poi) then
            CR_P1 = Crack3D_Coor(i_C,i_Edge,1:3)
            CR_P2 = Crack3D_Coor(i_C,i_Edge+1,1:3)
          else
            CR_P1 = Crack3D_Coor(i_C,i_Edge,1:3)
            CR_P2 = Crack3D_Coor(i_C,1,1:3)
          endif
          
          length_line_P1P2= Tool_Function_2Point_Dis_3D(CR_P1,CR_P2)
          num_Division = int(ceiling(length_line_P1P2/Mesh_Size))
          
          ! The first point of the current edge
          c_Count = c_Count +1
          Crack3D_Meshed_Node(i_C)%row(c_Count,1:3) = CR_P1
          
          !Divide line P1-P2
          if(num_Division>=2) then
              call Cal_Equal_Division_Points_3D(CR_P1,CR_P2,num_Division,.False., Points_P1P2,Num_Div_P_P1P2)
              
              do i_P = 1,Num_Div_P_P1P2
                  c_Count = c_Count + 1
                                  
                  c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
                  if(c_Count > c_Max_N_Node_3D) then
                      ! Expand memory. NEWFTU2022110501.
                      call D3_Allocate_Crack_Memory(i_C,1,1)
                      !c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
                  endif 
              
                  Crack3D_Meshed_Node(i_C)%row(c_Count,1:3) =Points_P1P2(i_P,1:3)
              enddo
          endif 
          
      enddo
      ! Add the centroid
      c_Count = c_Count +1
      Crack3D_Meshed_Node(i_C)%row(c_Count,1:3) = CR_Center
      
      ! Total number of discrete fracture surface points
      num_of_nodes = c_Count
      Crack3D_Meshed_Node_num(i_C)  = num_of_nodes
      
      
      ! During detection, the number of division nodes cannot exceed 5000, 
      ! because Save_Files_Crack_3D(isub) can save a maximum of 5000 crack division nodes.
      if(c_Count>5000) then
          print *, "    Error :: too much nodes(>5000), error in D3_Mesh_Initial_Crack.f!" 
          call Warning_Message('S',Keywords_Blank) 
      endif  
      
      
      !.....................
      ! element Node Number
      !.....................
      c_Count = 0
      do i_Ele = 1,num_of_nodes -1
          c_Count = c_Count+1
           
          c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
          if(c_Count > c_Max_N_Node_3D) then
              ! Expand memory. NEWFTU2022110501.
              call D3_Allocate_Crack_Memory(i_C,1,1)
              !c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
          endif 
          
          node_1 = i_Ele
          if(i_Ele/=(num_of_nodes -1))then
              node_2 = i_Ele + 1
          else
              node_2 = 1
          endif
          node_3 = num_of_nodes
          
          Crack3D_Meshed_Ele(i_C)%row(c_Count,1) =  node_1
          Crack3D_Meshed_Ele(i_C)%row(c_Count,2) =  node_2
          Crack3D_Meshed_Ele(i_C)%row(c_Count,3) =  node_3
          
      enddo
      num_of_ele = c_Count
      Crack3D_Meshed_Ele_num(i_C)  = c_Count      
      
      !.............................................................................................
      ! Determine the outer boundary of the fracture surface discretization elements, and store 
      ! sthe results in the following two global variables:
      ! integer Crack3D_Meshed_Outline(Max_Num_Cr, Max_N_Node_3D, 3) !3D crack outer boundary after
      ! discretization (data 3 corresponds to the element number)
      ! integer Crack3D_Meshed_Outline_num(Max_Num_Cr) !Number of 3D crack boundary lines after
      ! discretization
      !.............................................................................................
      call D3_Find_Crack_Mesh_Outline(i_C,num_of_ele,num_of_nodes,Crack3D_Meshed_Ele(i_C)%row(1:num_of_ele,1:3))
      !............................................................................................
      ! Determine the element number where the discrete node of the crack surface is located and 
      ! store it in Cr3D_Meshed_Node_in_Ele_Num(i_C, i_Cr_Node)
      ! If the element number is 0, it indicates that the crack surface node is outside the model.
      ! Determine the local coordinates of the elements where the crack face discrete nodes are 
      ! located, and store them in Cr3D_Meshed_Node_in_Ele_Local(Max_Num_Cr,1000,3)
      !............................................................................................
      do i_Cr_Node =1,num_of_nodes
          call Cal_Ele_Num_by_Coors_3D(  &
                        Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,1), &
                        Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,2), &
                        Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,3),Ele_Num_Cache,in_Elem_num)      
           if(in_Elem_num==0) then
               print *,'    ERROR-2022110102 :: in_Elem_num=0 in D3_Mesh_Initial_Crack.f90!'
               call Warning_Message('S',Keywords_Blank)
           endif                        
          Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(i_Cr_Node) = in_Elem_num
          if (in_Elem_num>0)then
              call Cal_KesiYita_by_Coor_3D(Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,1:3),in_Elem_num,c_Kesi,c_Yita,c_Zeta)
              Cr3D_Meshed_Node_in_Ele_Local(i_C)%row(i_Cr_Node,1:3)=[c_Kesi,c_Yita,c_Zeta] 
          endif
      enddo
      
      
    !################################
    ! Case 2: Circular crack surface
    !################################
    elseif(sum(abs(Crack3D_Cir_Coor(i_C,1:7)))>Tol_20)then 
       cir_center(1:3)= Crack3D_Cir_Coor(i_C,1:3)
       ! Length of the element containing the center
       call Cal_Ele_Num_by_Coors_3D(cir_center(1),cir_center(2),cir_center(3),Ele_Num_Cache,Circle_Elem)
       if(Circle_Elem==0) then
           print *,'    ERROR-2022110101 :: Circle_Elem=0 in D3_Mesh_Initial_Crack.f90!'
           print *,'                        Crack number:',i_C
           print *,'                        cir_center:',cir_center(1:3)
           print *,"                        Maybe outside the model!"
           call Warning_Message('S',Keywords_Blank)
       endif
       Cricle_Elem_L = (Elem_Vol(Circle_Elem))**(ONE/THR)
       ! Normal vector of a circle
       cir_vec(1:3)= Crack3D_Cir_Coor(i_C,4:6)
       ! Normalize the vector
       leg=sqrt(sum(cir_vec**2))
       cir_vec = cir_vec/leg
       
       radius = Crack3D_Cir_Coor(i_C,7)
       circumference =  TWO*pi*radius
       num_divison =  int(circumference/Cricle_Elem_L)
       
       !2024-03-13.
       if(num_divison<=2)then
           print *,'    ERROR-2024031301 :: illegal num_divison in D3_Mesh_Initial_Crack.f90!'
           print *,'                        Num_divison:',num_divison
           print *,'                        Radius:',radius
           print *,'                        Circumference:',circumference
           print *,'                        Cricle_Elem_L:',Cricle_Elem_L           
           print *,"                        Maybe caused by too small cracks compared to element size!"
           call Warning_Message('S',Keywords_Blank)
       endif
       
       allocate(c_xyz(num_divison,3))
       c_xyz(1:num_divison,1:3) = ZR
       ! Discretize the circle
       ! Reference: http://blog.sina.com.cn/s/blog_622fbc040102wt9o.html
       do i_theta = 1,num_divison
           theta = (i_theta-1)*TWO*pi/num_divison
           ! a = cross(n, [1 0 0]); % cross n with i to get the vector a
           call Vector_Cross_Product_3(cir_vec,[ONE,ZR,ZR],a)   
           ! if ~any(a) % If a is a zero vector, cross n with j
           !   a=cross(n,[0 1 0]);
           !end                
           if(sum(abs(a))<=Tol_20) then
              call Vector_Cross_Product_3(cir_vec,[ZR,ONE,ZR],a)   
           endif
           ! b = cross(n, a)! Calculate the b vector
           call Vector_Cross_Product_3(cir_vec,a,b)   
           ! a = a / norm(a)! Normalize the vector a
           ! b = b / norm(b)! Normalizing the b vector
           call Vector_Norm2(3,a,norm_a)   
           call Vector_Norm2(3,b,norm_b)  
           a=a/norm_a
           b=b/norm_b
 
           c_xyz(i_theta,1)=cir_center(1)+radius*a(1)*cos(theta)+radius*b(1)*sin(theta)
           c_xyz(i_theta,2)=cir_center(2)+radius*a(2)*cos(theta)+radius*b(2)*sin(theta)
           c_xyz(i_theta,3)=cir_center(3)+radius*a(3)*cos(theta)+radius*b(3)*sin(theta)
      enddo
      ! Number of 3D fracture nodes after discretization
      num_of_nodes = num_divison + 1   
           
      c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
      if(num_of_nodes > c_Max_N_Node_3D) then
          ! Expand memory. NEWFTU2022110501.
          call D3_Allocate_Crack_Memory(i_C,1,1)
          !c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
      endif 
              
      Crack3D_Meshed_Node_num(i_C) = num_of_nodes
      ! Number of 3D crack elements after discretization
      num_of_ele = num_divison  
      Crack3D_Meshed_Ele_num(i_C)  = num_of_ele        
      ! Coordinates of 3D fracture nodes after discretization, with each fracture consisting of up to
      ! 1,000 points
      Crack3D_Meshed_Node(i_C)%row(1,1:3) = cir_center(1:3)  
      do i_node = 1,num_divison
          Crack3D_Meshed_Node(i_C)%row(i_node+1,1:3)=c_xyz(i_node,1:3)  
      enddo
      deallocate(c_xyz)
                
      ! Determine the element number where the discrete node of the crack surface 
      ! is located and store it in Cr3D_Meshed_Node_in_Ele_Num(i_C, i_Cr_Node)
      ! If the element number is 0, it indicates that the crack surface node 
      ! is outside the model.
      do i_Cr_Node =1,num_divison + 1  
          call Cal_Ele_Num_by_Coors_3D(   &
                        Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,1), &
                        Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,2), &
                        Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,3), &
                        Ele_Num_Cache,in_Elem_num)      
          if(in_Elem_num==0) then
               print *,'    WARN-2022110106 :: in_Elem_num=0 in D3_Mesh_Initial_Crack.f90!'
               print *,'                        Crack number:',i_C
               print *,"                        Coors:",Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,1:3)
               print *,"                        Maybe outside the model!"
               !call Warning_Message('S',Keywords_Blank)
          endif                         
          Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(i_Cr_Node) = in_Elem_num
          ! Determine the local coordinates of the elements containing the discrete nodes
          ! on the crack surface, and store them in Cr3D_Meshed_Node_in_Ele_Local(Max_Num_Cr,1000,3)
          if (in_Elem_num>0)then
              call Cal_KesiYita_by_Coor_3D(Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,1:3),in_Elem_num,c_Kesi,c_Yita,c_Zeta)
              Cr3D_Meshed_Node_in_Ele_Local(i_C)%row(i_Cr_Node,1:3)=[c_Kesi,c_Yita,c_Zeta] 
          endif
      enddo
      ! Node numbering of discrete fracture elements
      do i_Ele =1,num_of_ele-1 
          node_1 = 1
          node_2 = i_Ele + 1
          node_3 = i_Ele + 2
          Crack3D_Meshed_Ele(i_C)%row(i_Ele,1:3)=[node_1,node_2,node_3]  
      enddo          
      node_1 = 1
      node_2 = num_divison + 1
      node_3 = 2
      Crack3D_Meshed_Ele(i_C)%row(num_of_ele,1:3)=[node_1,node_2,node_3]  
      ! Determine the outer boundary of the fracture surface discretization elements, 
      ! and store the results in the following two global variables:
      ! integer Crack3D_Meshed_Outline(Max_Num_Cr, Max_N_Node_3D, 4) !3D crack outer boundary after
      ! discretization (data 3 corresponds to the discrete crack element number)
      ! integer Crack3D_Meshed_Outline_num(Max_Num_Cr) !Number of 3D crack boundary lines after
      ! discretization
      call D3_Find_Crack_Mesh_Outline(i_C,num_of_ele,num_of_nodes,Crack3D_Meshed_Ele(i_C)%row(1:num_of_ele,1:3))          
    !#####################################
    ! Case 3: Elliptical Fracture Surface
    !#####################################
    elseif(sum(abs(Crack3D_Ellip_Coor(i_C,1:8)))>Tol_20)then 
       cir_center(1:3)= Crack3D_Ellip_Coor(i_C,1:3)
       cir_vec(1:3)= Crack3D_Ellip_Coor(i_C,4:6)
       ! Length of the element containing the ellipse
       call Cal_Ele_Num_by_Coors_3D(cir_center(1),cir_center(2),cir_center(3),Ele_Num_Cache,Circle_Elem)
       if(Circle_Elem==0) then
           print *,'    WARN-2022110107 :: Circle_Elem=0 in D3_Mesh_Initial_Crack.f90!'
           print *,'                        Crack number:',i_C
           print *,'                        cir_center:',cir_center(1:3)
           print *,"                        Maybe outside the model!"
           !call Warning_Message('S',Keywords_Blank)
       endif          
       Cricle_Elem_L = (Elem_Vol(Circle_Elem))**(ONE/THR)
       ! Normal of the circle?
       ! Normalize the vector
       leg=sqrt(sum(cir_vec**2))
       cir_vec = cir_vec/leg
       
       radius_a = Crack3D_Ellip_Coor(i_C,7)
       radius_b = Crack3D_Ellip_Coor(i_C,8)    
       circumference =  TWO*pi*radius_b + FOU*(radius_a-radius_b)
       num_divison =  int(circumference/Cricle_Elem_L)
       allocate(c_xyz(num_divison,3))
       c_xyz(1:num_divison,1:3) = ZR           
       ! Discretize the circle
       ! Reference: http://blog.sina.com.cn/s/blog_622fbc040102wt9o.html
       !         https://blog.csdn.net/qq_43539817/article/details/83475392
       do i_theta = 1,num_divison
           theta = (i_theta-1)*TWO*pi/num_divison
          ! a = cross(n, [1 0 0]); % cross n with i to get the vector a
          call Vector_Cross_Product_3(cir_vec,[ONE,ZR,ZR],a)   
          ! if ~any(a) % If a is a zero vector, cross n with j
          !   a=cross(n,[0 1 0]);
          !end                
          if(sum(abs(a))<=Tol_20) then
             call Vector_Cross_Product_3(cir_vec,[ZR,ONE,ZR],a)   
          endif
          ! b = cross(n, a)! Calculate the b vector
          call Vector_Cross_Product_3(cir_vec,a,b)   
          ! a = a / norm(a)! Normalize the vector a
          ! b = b / norm(b)! Normalizing the b vector
          call Vector_Norm2(3,a,norm_a)   
          call Vector_Norm2(3,b,norm_b)  
          a=a/norm_a
          b=b/norm_b
          c_xyz(i_theta,1)=cir_center(1)+radius_a*a(1)*cos(theta)+radius_b*b(1)*sin(theta)
          c_xyz(i_theta,2)=cir_center(2)+radius_a*a(2)*cos(theta)+radius_b*b(2)*sin(theta)
          c_xyz(i_theta,3)=cir_center(3)+radius_a*a(3)*cos(theta)+radius_b*b(3)*sin(theta)
      enddo
      ! Number of 3D fracture nodes after discretization
      num_of_nodes = num_divison + 1  
      
      c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
      if(num_of_nodes > c_Max_N_Node_3D) then
          ! Expand memory. NEWFTU2022110501.
          call D3_Allocate_Crack_Memory(i_C,1,1)
          !c_Max_N_Node_3D = size(Crack3D_Meshed_Node(i_C)%row,1)
      endif 
      
      Crack3D_Meshed_Node_num(i_C) = num_of_nodes
      ! Number of 3D crack elements after discretization
      num_of_ele = num_divison  
      Crack3D_Meshed_Ele_num(i_C)  = num_of_ele        
      ! Coordinates of 3D fracture nodes after discretization, with each fracture consisting of up to
      ! 1,000 points
      Crack3D_Meshed_Node(i_C)%row(1,1:3) = cir_center(1:3)  
      do i_node = 1,num_divison
          Crack3D_Meshed_Node(i_C)%row(i_node+1,1:3)=c_xyz(i_node,1:3)
      enddo
      deallocate(c_xyz)
      ! Determine the element number where the discrete node of the crack surface is 
      ! located and store it in Cr3D_Meshed_Node_in_Ele_Num(i_C, i_Cr_Node)
      ! If the element number is 0, it indicates that the crack surface node is 
      ! outside the model.
      do i_Cr_Node =1,num_divison + 1  
          call Cal_Ele_Num_by_Coors_3D(                            & 
                        Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,1), &
                        Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,2), &
                        Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,3), &
                        Ele_Num_Cache,in_Elem_num)          
          if(in_Elem_num==0) then
               print *,'    ERROR-2022110108 :: in_Elem_num=0 in D3_Mesh_Initial_Crack.f90!'
               call Warning_Message('S',Keywords_Blank)
          endif                          
          Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(i_Cr_Node) = in_Elem_num
          ! Determine the local coordinates of the elements containing the discrete nodes on the crack
          ! surface, and store them in Cr3D_Meshed_Node_in_Ele_Local(Max_Num_Cr,1000,3)
          if (in_Elem_num>0)then
              call Cal_KesiYita_by_Coor_3D(Crack3D_Meshed_Node(i_C)%row(i_Cr_Node,1:3),in_Elem_num,c_Kesi,c_Yita,c_Zeta)
              Cr3D_Meshed_Node_in_Ele_Local(i_C)%row(i_Cr_Node,1:3)=[c_Kesi,c_Yita,c_Zeta] 
          endif
      enddo
      ! Node numbering of discrete fracture elements
      do i_Ele =1,num_of_ele-1 
          node_1 = 1
          node_2 = i_Ele + 1
          node_3 = i_Ele + 2
          Crack3D_Meshed_Ele(i_C)%row(i_Ele,1:3) =[node_1,node_2,node_3]  
      enddo          
      node_1 = 1
      node_2 = num_divison + 1
      node_3 = 2
      Crack3D_Meshed_Ele(i_C)%row(num_of_ele,1:3) =[node_1,node_2,node_3]  
      ! Determine the outer boundary of the crack surface discrete element, 
      ! and store the result in the following global variable:
      ! integer Crack3D_Meshed_Outline(Max_Num_Cr, Max_N_Node_3D, 3) !3D crack outer boundary after
      ! discretization (data 3 corresponds to the element number)
      ! integer Crack3D_Meshed_Outline_num(Max_Num_Cr) !Number of 3D crack boundary lines after
      ! discretization
      call D3_Find_Crack_Mesh_Outline(i_C,num_of_ele,num_of_nodes,Crack3D_Meshed_Ele(i_C)%row(1:num_of_ele,1:3))            
    endif
enddo
  
  
!-----------------------------------------------------------------------
! Calculate the properties of discrete crack elements and store them in 
! Crack3D_Meshed_Ele_Attri(i_Crack_Ele, 1:5)
!-----------------------------------------------------------------------
do i_C =1,num_Crack  
  ! If the crack has already been discretized previously, proceed to the next loop. IMPROV2022053101.
  if(Cracks_Initial_Meshed(i_C)==1)then
      cycle
  endif    
  do i_Crack_Ele =1,Crack3D_Meshed_Ele_num(i_C)
      Crack_Node1 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,1)
      Crack_Node2 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,2)
      Crack_Node3 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,3)
      Point1(1:3) =Crack3D_Meshed_Node(i_C)%row(Crack_Node1,1:3)
      Point2(1:3) =Crack3D_Meshed_Node(i_C)%row(Crack_Node2,1:3)
      Point3(1:3) =Crack3D_Meshed_Node(i_C)%row(Crack_Node3,1:3)
      ! Perimeter of a discrete element (spatial triangle)
      Crack3D_Meshed_Ele_Attri(i_C)%row(i_Crack_Ele,1)  = Tool_Function_2Point_Dis_3D(Point1,Point2) + &
                                                     Tool_Function_2Point_Dis_3D(Point2,Point3) + &
                                                     Tool_Function_2Point_Dis_3D(Point1,Point3)
      ! Area of discrete elements (spatial triangles)
      ! ABC(1,1:3) = Point1
      ! ABC(2,1:3) = Point2
      ! ABC(3,1:3) = Point3
      call Tool_Area_Tri_3D(Point1,Point2,Point3,Crack3D_Meshed_Ele_Attri(i_C)%row(i_Crack_Ele,2))
      ! Outer normal vector Crack3D_Meshed_Ele_Nor_Vector(Max_Num_Cr_3D, Max_N_Node_3D, 3)
      Vector_V1 = Point2- Point1
      Vector_V2 = Point3- Point1
      call Vector_Normalize(3,Vector_V1)  
      call Vector_Normalize(3,Vector_V2)  
      call Vector_Cross_Product_3(Vector_V1,Vector_V2,Vector_Normal) 
      call Vector_Normalize(3,Vector_Normal)  
      Crack3D_Meshed_Ele_Nor_Vector(i_C)%row(i_Crack_Ele,1:3)=Vector_Normal
  enddo
enddo
  
  
!----------------------------------------------------------------------
! Outer normal vectors of 3D crack nodes after discretization
!        Crack3D_Meshed_Node_Nor_Vector(Max_Num_Cr_3D,Max_N_Node_3D,3) 
!----------------------------------------------------------------------

do i_C =1,num_Crack    
  ! If the crack has already been discretized previously, proceed to the next loop. IMPROV2022053101.
  if(Cracks_Initial_Meshed(i_C)==1)then
      cycle
  endif     
  Crack3D_Meshed_Node_Nor_Vector(i_C)%row(1:Max_N_Node_3D(i_C),1:3) = ZR          
  do i_Crack_Node =1,Crack3D_Meshed_Node_num(i_C) 
      c_count_Ele = 0
      c_All_Nor_Vector = [ZR,ZR,ZR]
      do i_Crack_Ele =1,Crack3D_Meshed_Ele_num(i_C)
          Crack_Node1=Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,1)
          Crack_Node2=Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,2)
          Crack_Node3=Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,3)
          if(i_Crack_Node == Crack_Node1 .or.i_Crack_Node == Crack_Node2 .or. i_Crack_Node == Crack_Node3 )then
              c_All_Nor_Vector  = c_All_Nor_Vector  +Crack3D_Meshed_Ele_Nor_Vector(i_C)%row(i_Crack_Ele,1:3)
              c_count_Ele = c_count_Ele +1
          endif
      enddo
      Crack3D_Meshed_Node_Nor_Vector(i_C)%row(i_Crack_Node,1:3)  = c_All_Nor_Vector/c_count_Ele 
      call Vector_Normalize(3,Crack3D_Meshed_Node_Nor_Vector(i_C)%row(i_Crack_Node,1:3))
  enddo
enddo

!----------------------------------------------------------
! Centroid of the 3D fracture surface after discretization
!----------------------------------------------------------
do i_C =1,num_Crack  
    ! If the crack has already been discretized previously, proceed to the next loop. IMPROV2022053101.
    if(Cracks_Initial_Meshed(i_C)==1)then
        cycle
    endif      
    num_Cr_Nodes = Crack3D_Meshed_Node_num(i_C)
    centroid_x=sum(Crack3D_Meshed_Node(i_C)%row(1:num_Cr_Nodes,1))/dble(num_Cr_Nodes)
    centroid_y=sum(Crack3D_Meshed_Node(i_C)%row(1:num_Cr_Nodes,2))/dble(num_Cr_Nodes)
    centroid_z=sum(Crack3D_Meshed_Node(i_C)%row(1:num_Cr_Nodes,3))/dble(num_Cr_Nodes)
    Crack3D_Centroid(i_C,1) = centroid_x
    Crack3D_Centroid(i_C,2) = centroid_y
    Crack3D_Centroid(i_C,3) = centroid_z
enddo      
  
!----------------------------------------------------------------------
! Check whether the crack has been previously initialized as discrete. 
! IMPROV2022053101.
!----------------------------------------------------------------------
do i_C =1,num_Crack
    ! If the crack has not been discretized before, mark it as discretized.
    if(Cracks_Initial_Meshed(i_C)==0)then
        Cracks_Initial_Meshed(i_C)=1
    endif  
enddo   
  

!----------------------------------------------------------------
! Boundary Simplification. IMPROV202450203.
! For points located on the same straight line, delete the part.
!----------------------------------------------------------------
if(Key_InPlane_Growth==1)then

    if(Key_Scheme_Signed_Dis_InPlane==3)then
        max_outline_num = maxval(Crack3D_Meshed_Outline_num(1:num_Crack))
        
        if(allocated(simplified_crack_outline_num)) deallocate(simplified_crack_outline_num)
        allocate(simplified_crack_outline_num(num_Crack))
        simplified_crack_outline_num(1:num_Crack) =0
        
        if(allocated(simplified_crack_outline_points)) deallocate(simplified_crack_outline_points)
        allocate(simplified_crack_outline_points(num_Crack,max_outline_num,3))
        simplified_crack_outline_points(1:num_Crack,1:max_outline_num,1:3) =ZR
        
        do i_C =1,num_Crack
        
            c_outline_num = Crack3D_Meshed_Outline_num(i_C)
            
            if (c_outline_num<=10) then
                cycle
            endif
            
            ALLOCATE(input_points(c_outline_num,3))
            ALLOCATE(output_points(c_outline_num,3))
            
            do i_P =1,c_outline_num
                input_points(i_P,1:3) = Crack3D_Meshed_Node(i_C)%row(Crack3D_Meshed_Outline(i_C)%row(i_P,1),1:3) 
            enddo
            
            ! Delete collinear points.
            ! Be careful not to delete the first point and the last point.
            call Tool_Delete_Collinear_Points_3D(input_points(1:c_outline_num,1:3),c_outline_num,&
                             output_points(1:c_outline_num,1:3),output_point_num)
                
            
            simplified_crack_outline_num(i_C) = output_point_num
            simplified_crack_outline_points(i_C,1:output_point_num,1:3) =output_points(1:output_point_num,1:3)
            
            DEALLOCATE(input_points)
            DEALLOCATE(output_points)
        enddo
    endif  
endif  
#endif


#ifdef Silverfrost
print *,'    ERROR :: Silverfrost compiler failed to compile D3_Mesh_Initial_Crack.F90!'
print *,'             In D3_Mesh_Initial_Crack.F90.'
call Warning_Message('S',Keywords_Blank)
#endif  
RETURN
END SUBROUTINE D3_Mesh_Initial_Crack
