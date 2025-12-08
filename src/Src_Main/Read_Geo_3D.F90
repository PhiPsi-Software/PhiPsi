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
 
SUBROUTINE Read_Geo_3D
  
!-----------------------------
! Read public variable module
!-----------------------------
use Global_Float_Type
use Global_Common
use Global_Filename
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack_Common
use Global_Crack_3D
use Global_HF
use Global_Contact
use Global_POST
use Global_Stress
use omp_lib

!---------------------------
! Variable Type Declaration
!---------------------------
implicit none
LOGICAL ALIVE
character*200 temp_name
integer Tool_Count_Lines

!--------------------
! Temporary variable
!--------------------
real(kind=FT),ALLOCATABLE::Temp_DATA(:,:)
logical Flag_Blank
integer,ALLOCATABLE::tem1(:,:)
integer,ALLOCATABLE::tem2(:)
integer,ALLOCATABLE::tem3(:)
integer,ALLOCATABLE::tem4(:)
integer i,j
integer all_num_Outline
integer,ALLOCATABLE::All_Outline(:,:)
integer N1,N2,N3,N4,N5,N6,N7,N8,NN(8)
real(kind=FT) Vol
integer num_Outline
character(200) c_File_name_1
character(200) temp_string
integer Num_CkElem
real(kind=FT) EdgeL(12),Tool_Function_2Point_Dis_3D
integer i_E,i_count,c_E,i_Edge,c_node_1,c_node_2,tem_num_edges
real(kind=FT) c_1_x,c_1_y,c_1_z,c_2_x,c_2_y,c_2_z
real(kind=FT),ALLOCATABLE::Elem_c_L_min(:),Elem_c_L_max(:)
integer c_NN(8),c_Node
integer n
integer c_node_1_num_Ele,c_node_2_num_Ele
integer,ALLOCATABLE::tem_vector1(:),tem_vector2(:)
integer temp_variable
integer Uniqued_Vec_8(8)
real(kind=FT) P1(3),P2(3)
integer tem_Elem_Node(8),tem_Uniqued_Vec_8(8)
#ifdef Silverfrost  
real(kind=FT),ALLOCATABLE::tem_vector_real(:) 
#endif
integer i_Mat,c_Mat_Type


! If Key_Cpp_Call_Fortran_Lib=1, the file reading is skipped. 2023-03-23.
if (Key_Cpp_Call_Fortran_Lib==1) then 
    goto 600
endif

!------------------------------------------------------------------------
!     
! STEP 1: Check if the node coordinate file exists; if it does, read it.
!
!------------------------------------------------------------------------
print *, "    Trying to read nodal files..." 
temp_name = trim(trim(Full_Pathname)//'.node')
inquire(file=temp_name, exist=alive)
if(alive.EQV..FALSE.)then
  print *, "    Error :: Can not find nodal files!" 
  print *, "             Check whether the following file exists or not:"
  print *, "             ",temp_name       
  call Warning_Message('S',Keywords_Blank) 
else
  ! Check whether the node file is a 3D problem file
  open(101,file=temp_name,status='old')
      read(101,'(A)') temp_string
      !read(101,*) temp_value(3)
  close(101)
  
  ! Count file lines
  Num_Node = Tool_Count_Lines(temp_name)
  IF(ALLOCATED(Coor)) DEALLOCATE(Coor)
  ALLOCATE(Coor(Num_Node,3))
  ALLOCATE(Temp_DATA(Num_Node,3))
  Call Tool_Read_File(temp_name,"node",Num_Node,3,Temp_DATA,Flag_Blank)
  Coor = Temp_DATA
  DEALLOCATE(Temp_DATA)
endif

!-----------------------------------------------------------------------   
!   
! STEP 2: Check if the element number file exists; if it does, read it.
!
!-----------------------------------------------------------------------
print *, "    Trying to read element files..." 
temp_name = trim(trim(Full_Pathname)//'.elem')
inquire(file=temp_name, exist=alive)  
if(alive.EQV..FALSE.)then
  print *, "    Error :: Can not find element files!!!" 
  call Warning_Message('S',Keywords_Blank) 
else
  Num_Elem = Tool_Count_Lines(temp_name)
  IF(ALLOCATED(Elem_Node)) DEALLOCATE(Elem_Node)
  ALLOCATE(Elem_Node(Num_Elem,8))
  IF(ALLOCATED(Elem_Mat)) DEALLOCATE(Elem_Mat)
  ALLOCATE(Elem_Mat(Num_Elem)) 
  ALLOCATE(Temp_DATA(Num_Elem,9))
  Call Tool_Read_File(temp_name,"elem",Num_Elem,9,Temp_DATA,Flag_Blank)
  Elem_Node = int(Temp_DATA(:,1:8))
  Elem_Mat  = int(Temp_DATA(:,9))
  num_of_Material = MaxVal(Elem_Mat)
  DEALLOCATE(Temp_DATA)
  !endif
endif  


!-------------------------------------------------------------------
!
! STEP 3: Check if the x boundary constraint file exists, and if it 
!         does, read it.
!
!-------------------------------------------------------------------
print *, "    Trying to read boux files..." 
temp_name = trim(trim(Full_Pathname)//'.boux')
inquire(file=temp_name, exist=alive)  
if(alive.EQV..FALSE.)then
  print *, "    Warning :: Can not find boux files!!!" 
else
  Num_Bou_x = Tool_Count_Lines(temp_name) 
  IF(ALLOCATED(Bou_x)) DEALLOCATE(Bou_x)
  ALLOCATE( Bou_x(Num_Bou_x))
  ALLOCATE( Temp_DATA(Num_Bou_x,1))
  Call Tool_Read_File(temp_name,"boux",Num_Bou_x,1,Temp_DATA,Flag_Blank)
  Bou_x  = int(Temp_DATA(:,1))
  DEALLOCATE(Temp_DATA)
endif  

!--------------------------------------------------------------------------------- 
!     
! STEP 4: Check if the y boundary constraint file exists, and if it does, read it
!
!---------------------------------------------------------------------------------
print *, "    Trying to read bouy files..." 
temp_name = trim(trim(Full_Pathname)//'.bouy')
inquire(file=temp_name, exist=alive)  
if(alive.EQV..FALSE.)then
  print *, "    Warning :: Can not find bouy files!!!" 
else
  Num_Bou_y = Tool_Count_Lines(temp_name) 
  IF(ALLOCATED(Bou_y)) DEALLOCATE(Bou_y)
  ALLOCATE( Bou_y(Num_Bou_y))
  ALLOCATE( Temp_DATA(Num_Bou_y,1))
  Call Tool_Read_File(temp_name,"bouy",Num_Bou_y,1,Temp_DATA,Flag_Blank)
  Bou_y  = int(Temp_DATA(:,1))
  DEALLOCATE(Temp_DATA)
endif 

!---------------------------------------------------------------------------------  
!    
! STEP 5: Check if the z boundary constraint file exists, and if it does, read it
!
!---------------------------------------------------------------------------------
print *, "    Trying to read bouz files..." 
temp_name = trim(trim(Full_Pathname)//'.bouz')
inquire(file=temp_name, exist=alive)  
if(alive.EQV..FALSE.)then
  print *, "    Warning :: Can not find bouz files!!!" 
else
  Num_Bou_z = Tool_Count_Lines(temp_name) 
  IF(ALLOCATED(Bou_z)) DEALLOCATE(Bou_z)
  ALLOCATE( Bou_z(Num_Bou_z))
  ALLOCATE( Temp_DATA(Num_Bou_z,1))
  Call Tool_Read_File(temp_name,"bouz",Num_Bou_z,1,Temp_DATA,Flag_Blank)
  Bou_z  = int(Temp_DATA(:,1))
  DEALLOCATE(Temp_DATA)
endif 

!-----------------------------------------------------------------------------------------
!      
! STEP 6: Check whether the load file in the x direction exists, and if it does, read it.
!
!-----------------------------------------------------------------------------------------
print *, "    Trying to read focx files..." 
temp_name = trim(trim(Full_Pathname)//'.focx')
inquire(file=temp_name, exist=alive)  
if(alive.EQV..FALSE.)then
  print *, "    Warning :: Can not find focx files!!!" 
else
  Num_Foc_x = Tool_Count_Lines(temp_name) 
  IF(ALLOCATED(Foc_x)) DEALLOCATE(Foc_x)
  ALLOCATE( Foc_x(Num_Foc_x,2))
  ALLOCATE( Temp_DATA(Num_Foc_x,2))
  Call Tool_Read_File(temp_name,"focx",Num_Foc_x,2,Temp_DATA,Flag_Blank)
  Foc_x  = Temp_DATA(:,:)
  DEALLOCATE(Temp_DATA)
endif  

!-----------------------------------------------------------------------------------------
!    
! STEP 7: Check whether the load file in the y direction exists, and if it does, read it.
!
!-----------------------------------------------------------------------------------------
print *, "    Trying to read focy files..." 
temp_name = trim(trim(Full_Pathname)//'.focy')
inquire(file=temp_name, exist=alive)  
if(alive.EQV..FALSE.)then
  print *, "    Warning :: Can not find focy files!!!" 
else
  Num_Foc_y = Tool_Count_Lines(temp_name) 
  IF(ALLOCATED(Foc_y)) DEALLOCATE(Foc_y)
  ALLOCATE( Foc_y(Num_Foc_y,2))
  ALLOCATE( Temp_DATA(Num_Foc_y,2))
  Call Tool_Read_File(temp_name,"focy",Num_Foc_y,2,Temp_DATA,Flag_Blank)
  Foc_y  = Temp_DATA(:,:)
  DEALLOCATE(Temp_DATA)
endif 

!--------------------------------------------------------------------------------
!   
! STEP 8: Check if the load file in the z direction exists. If it does, read it.
!
!--------------------------------------------------------------------------------
print *, "    Trying to read focz files..." 
temp_name = trim(trim(Full_Pathname)//'.focz')
inquire(file=temp_name, exist=alive)  
if(alive.EQV..FALSE.)then
  print *, "    Warning :: Can not find focz files!!!" 
else
  Num_Foc_z = Tool_Count_Lines(temp_name) 
  IF(ALLOCATED(Foc_z)) DEALLOCATE(Foc_z)
  ALLOCATE( Foc_z(Num_Foc_z,2))
  ALLOCATE( Temp_DATA(Num_Foc_z,2))
  Call Tool_Read_File(temp_name,"focz",Num_Foc_z,2,Temp_DATA,Flag_Blank)
  Foc_z  = Temp_DATA(:,:)
  DEALLOCATE(Temp_DATA)
endif 

600 continue

!---------------------------------------------------------------------------
!
! STEP 9: Store the coordinates of the 8 nodes of the element and calculate
! the average area of the element, etc.
! And other data extraction.
!
!---------------------------------------------------------------------------
print *, "    Get node and element data..." 
! Obtain the maximum node number difference in the element to calculate the maximum half-bandwidth.
! Ref: A First Course in the Finite Element Method, 5th edition, B.4.
Max_Diff_Elem_Num = maxval(maxval(Elem_Node(:,1:8),2)- minval(Elem_Node(:,1:8),2))
Max_Half_Band_Width = 3*(Max_Diff_Elem_Num + 1)
!
! Coordinate extrema
Max_X_Coor = maxval(Coor(:,1))
Min_X_Coor = minval(Coor(:,1))
Max_Y_Coor = maxval(Coor(:,2))
Min_Y_Coor = minval(Coor(:,2))
Max_Z_Coor = maxval(Coor(:,3))
Min_Z_Coor = minval(Coor(:,3))
! Model center coordinates. 2022-11-21. IMPROV2022112102.
Model_Center_x = (Max_X_Coor + Min_X_Coor)/TWO 
Model_Center_y = (Max_Y_Coor + Min_Y_Coor)/TWO
Model_Center_z = (Max_Z_Coor + Min_Z_Coor)/TWO
!Model coordinates range.
Model_X_Range = Max_X_Coor-Min_X_Coor
Model_Y_Range = Max_Y_Coor-Min_Y_Coor
Model_Z_Range = Max_Z_Coor-Min_Z_Coor
Max_Model_Range=max(Model_X_Range,Model_Y_Range,Model_Z_Range)
Min_Model_Range=min(Model_X_Range,Model_Y_Range,Model_Z_Range)
IF(ALLOCATED(G_NN)) DEALLOCATE(G_NN)
ALLOCATE(G_NN(8,Num_Elem))
IF(ALLOCATED(G_X_NODES)) DEALLOCATE(G_X_NODES)
ALLOCATE(G_X_NODES(8,Num_Elem))
IF(ALLOCATED(G_Y_NODES)) DEALLOCATE(G_Y_NODES)
ALLOCATE(G_Y_NODES(8,Num_Elem))
IF(ALLOCATED(G_Z_NODES)) DEALLOCATE(G_Z_NODES)
ALLOCATE(G_Z_NODES(8,Num_Elem))
IF(ALLOCATED(Elem_Vol)) DEALLOCATE(Elem_Vol)
ALLOCATE(Elem_Vol(Num_Elem))
IF(ALLOCATED(Elem_Max_L)) DEALLOCATE(Elem_Max_L)
ALLOCATE(Elem_Max_L(Num_Elem))
IF(ALLOCATED(Elem_Min_L)) DEALLOCATE(Elem_Min_L)
ALLOCATE(Elem_Min_L(Num_Elem))
IF(ALLOCATED(Elem_Ave_L)) DEALLOCATE(Elem_Ave_L)
ALLOCATE(Elem_Ave_L(Num_Elem))
IF(ALLOCATED(Node_Max_L)) DEALLOCATE(Node_Max_L)
ALLOCATE(Node_Max_L(Num_Node))
IF(ALLOCATED(Elem_Centroid)) DEALLOCATE(Elem_Centroid)
ALLOCATE(Elem_Centroid(Num_Elem,3))

IF(ALLOCATED(EleGaus_yes_FEM_asemd)) then
  DEALLOCATE(EleGaus_yes_FEM_asemd)
endif

! ALLOCATE(EleGaus_yes_FEM_asemd(Num_Elem, Num_Gau_Points_3D)) ! Used to mark whether the FEM
! element stiffness at this Gauss point for this element has already been assembled into the global
! stiffness
Num_Max_3D_gauss = max(Num_Gau_Points_3D,Num_Gau_Points_3D_MC)
if(Key_Integral_Sol  ==3 ) then
 Num_Max_3D_gauss = max(Num_Max_3D_gauss,Num_Sub_3D_Cubes*Num_Gau_Points_3D_Cube) 
endif
ALLOCATE(EleGaus_yes_FEM_asemd(Num_Elem,Num_Max_3D_gauss))


IF(ALLOCATED(x_max_Elements)) DEALLOCATE(x_max_Elements)
ALLOCATE(x_max_Elements(Num_Elem)) 
IF(ALLOCATED(x_min_Elements)) DEALLOCATE(x_min_Elements)
ALLOCATE(x_min_Elements(Num_Elem)) 
IF(ALLOCATED(y_max_Elements)) DEALLOCATE(y_max_Elements)
ALLOCATE(y_max_Elements(Num_Elem)) 
IF(ALLOCATED(y_min_Elements)) DEALLOCATE(y_min_Elements)
ALLOCATE(y_min_Elements(Num_Elem))
IF(ALLOCATED(z_max_Elements)) DEALLOCATE(z_max_Elements)
ALLOCATE(z_max_Elements(Num_Elem)) 
IF(ALLOCATED(z_min_Elements)) DEALLOCATE(z_min_Elements)
ALLOCATE(z_min_Elements(Num_Elem)) 
if(Key_Post_Elements_Gauss_Num== 1) then
IF(ALLOCATED(Elements_Gauss_Num)) DEALLOCATE(Elements_Gauss_Num)
    ALLOCATE(Elements_Gauss_Num(Num_Elem))     
endif

EleGaus_yes_FEM_asemd(1:Num_Elem,1:Num_Max_3D_gauss) =.false.
ALLOCATE(Elem_c_L_min(Num_Elem))
ALLOCATE(Elem_c_L_max(Num_Elem))      
!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,N1,N2,N3,N4,N5,N6,N7,N8,NN,EdgeL)    
!BUGFIX2022092401.
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,N1,N2,N3,N4,N5,N6,N7,N8,NN,EdgeL,Vol,P1,P2)   
do i=1,Num_Elem
  N1  = Elem_Node(i,1)                                                
  N2  = Elem_Node(i,2)                                              
  N3  = Elem_Node(i,3)                                             
  N4  = Elem_Node(i,4)            
  N5  = Elem_Node(i,5)   
  N6  = Elem_Node(i,6)   
  N7  = Elem_Node(i,7)   
  N8  = Elem_Node(i,8)   
  NN  = [N1,N2,N3,N4,N5,N6,N7,N8]                                                
  G_NN(1:8,i)      = NN
  G_X_NODES(1:8,i) = Coor(NN,1)
  G_Y_NODES(1:8,i) = Coor(NN,2) 
  G_Z_NODES(1:8,i) = Coor(NN,3) 
  ! The range of x, y, z coordinates for each element
  x_max_Elements(i) = maxval(G_X_NODES(1:8,i))
  x_min_Elements(i) = minval(G_X_NODES(1:8,i))
  y_max_Elements(i) = maxval(G_Y_NODES(1:8,i))
  y_min_Elements(i) = minval(G_Y_NODES(1:8,i))
  z_max_Elements(i) = maxval(G_Z_NODES(1:8,i))
  z_min_Elements(i) = minval(G_Z_NODES(1:8,i))
  ! element side length, obtaining the minimum element side length for explicit dynamic analysis
  ! (2021-08-23)
  P1 = Coor(N1,1:3); P2 = Coor(N2,1:3)
  EdgeL(1)=Tool_Function_2Point_Dis_3D(P1,P2)
  P1 = Coor(N2,1:3); P2 = Coor(N3,1:3)
  EdgeL(2)=Tool_Function_2Point_Dis_3D(P1,P2)
  P1 = Coor(N3,1:3); P2 = Coor(N4,1:3)
  EdgeL(3)=Tool_Function_2Point_Dis_3D(P1,P2)
  P1 = Coor(N4,1:3); P2 = Coor(N1,1:3)
  EdgeL(4)=Tool_Function_2Point_Dis_3D(P1,P2)
  P1 = Coor(N1,1:3); P2 = Coor(N5,1:3)
  EdgeL(5)=Tool_Function_2Point_Dis_3D(P1,P2)
  P1 = Coor(N2,1:3); P2 = Coor(N6,1:3)
  EdgeL(6)=Tool_Function_2Point_Dis_3D(P1,P2)
  P1 = Coor(N3,1:3); P2 = Coor(N7,1:3)
  EdgeL(7)=Tool_Function_2Point_Dis_3D(P1,P2)
  P1 = Coor(N4,1:3); P2 = Coor(N8,1:3)
  EdgeL(8)=Tool_Function_2Point_Dis_3D(P1,P2)
  P1 = Coor(N5,1:3); P2 = Coor(N6,1:3)
  EdgeL(9)=Tool_Function_2Point_Dis_3D(P1,P2)
  P1 = Coor(N6,1:3); P2 = Coor(N7,1:3)
  EdgeL(10)=Tool_Function_2Point_Dis_3D(P1,P2)
  P1 = Coor(N7,1:3); P2 = Coor(N8,1:3)
  EdgeL(11)=Tool_Function_2Point_Dis_3D(P1,P2)
  P1 = Coor(N5,1:3); P2 = Coor(N8,1:3)
  EdgeL(12)=Tool_Function_2Point_Dis_3D(P1,P2)
  ! Maximum element edge length. IMPROV2023022202.
  Elem_Max_L(i) = maxval(EdgeL(1:12))
  ! Minimum element edge length. IMPROV2023022702.
  Elem_Min_L(i) = minval(EdgeL(1:12))
  ! Average element edge length. IMPROV2023022702.
  Elem_Ave_L(i) = sum(EdgeL(1:12))/12.0D0
  Elem_c_L_min(i) = minval(EdgeL(1:12))
  Elem_c_L_max(i) = maxval(EdgeL(1:12))
  !Volume of element.
  call Cal_Vol_8nodes(i,Coor(NN,1),Coor(NN,2),Coor(NN,3),Vol)
  Elem_Vol(i) = Vol
  !Centroid of element.
  call Cal_Coor_by_KesiYita_3D(ZR,ZR,ZR,G_X_NODES(1:8,i),G_Y_NODES(1:8,i),G_Z_NODES(1:8,i), &
                     Elem_Centroid(i,1),Elem_Centroid(i,2),Elem_Centroid(i,3))
end do
!$OMP END PARALLEL DO

Min_Ele_Edge_Length = minval(Elem_c_L_min)
Max_Ele_Edge_Length = maxval(Elem_c_L_max)
DEALLOCATE(Elem_c_L_min)
DEALLOCATE(Elem_c_L_max)

Max_Elem_Vol =MaxVal(Elem_Vol)
Min_Elem_Vol =MinVal(Elem_Vol)
Ave_Elem_Vol =sum(Elem_Vol(1:Num_Elem))/size(Elem_Vol(1:Num_Elem))
!Ave_Elem_L   = sqrt(sqrt(Ave_Elem_Vol)) !Bug fixed, 2021-04-16
Ave_Elem_L   =Ave_Elem_Vol**(ONE/THR)
Ave_Elem_Area = Ave_Elem_L**2

!---------------------------------------------------------------------------------------
!
! STEP 9.1: Number of independent nodes for each element. Used for degenerate elements.
!           IMPROV2023061402.
! 
!---------------------------------------------------------------------------------------
ALLOCATE(Elem_Uniqued_Nodes(Num_Elem))
Elem_Uniqued_Nodes = 0
ALLOCATE(Yes_Degenarated_Elem(Num_Elem))
Yes_Degenarated_Elem = .False.
Uniqued_Vec_8(1:8) = 0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,Uniqued_Vec_8,tem_Elem_Node,tem_Uniqued_Vec_8)   
do i=1,Num_Elem 
    tem_Elem_Node(1:8)     = Elem_Node(i,1:8)
    tem_Uniqued_Vec_8(1:8) = Uniqued_Vec_8(1:8)
    call Vector_Unique_Int(8,8,tem_Elem_Node,tem_Uniqued_Vec_8,Elem_Uniqued_Nodes(i))   
    if (Elem_Uniqued_Nodes(i)<8)then
        Yes_Degenarated_Elem(i) = .True.
    endif
enddo
!$OMP END PARALLEL DO
Num_Degenarated_Elems = count(Yes_Degenarated_Elem .eqv. .True.)
print *, '    Number of degenarated elements: ', Num_Degenarated_Elems
print '(A,F8.3)', '     Percent of degenarated elements (%):', dble(Num_Degenarated_Elems)/dble(Num_Elem)*TEN_2
!--------------------------------------------------------------------
!
! STEP 10: The elements around each node and the number of elements.
!
!--------------------------------------------------------------------
print *, '    Get surounding elements of nodes...'   
IF(ALLOCATED(num_Node_Elements)) DEALLOCATE(num_Node_Elements)
ALLOCATE(num_Node_Elements(Num_Node))
num_Node_Elements(1:Num_Node) = 0  

IF(ALLOCATED(Node_Elements_3D)) DEALLOCATE(Node_Elements_3D)
ALLOCATE(Node_Elements_3D(Num_Node))
 
!STEP 1.

!STEP 1.1.
!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,c_NN,c_Node)   
! Not suitable for parallel computing. There are data conflicts.
do j=1,Num_Elem                  
  c_NN    = G_NN(1:8,j)
  do i=1,8                
      c_Node = c_NN(i)
      num_Node_Elements(c_Node) = num_Node_Elements(c_Node) +1
  end do
end do   
!!$omp end parallel do

!STEP 1.2. 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)   
do i = 1,Num_Node    
    allocate(Node_Elements_3D(i)%row(num_Node_Elements(i)))
    Node_Elements_3D(i)%row(1:num_Node_Elements(i)) = 0
enddo      
!$omp end parallel do

!STEP 1.3.
num_Node_Elements(1:Num_Node) = 0  
do j=1,Num_Elem                  
  c_NN    = G_NN(1:8,j)
  do i=1,8                
      c_Node = c_NN(i)
      num_Node_Elements(c_Node) = num_Node_Elements(c_Node) +1
      Node_Elements_3D(c_Node)%row(num_Node_Elements(c_Node)) = j
  end do
end do         
  
! STEP 2. Remove duplicates.
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,n,tem_vector1,tem_vector2,temp_variable)    
do i=1,Num_Node  
    n = num_Node_Elements(i)
    allocate(tem_vector1(n),tem_vector2(n))
    tem_vector1(1:n) = Node_Elements_3D(i)%row(1:n)
    tem_vector2(1:n) = 0
    call Vector_Unique_Int(n,n,tem_vector1(1:n),tem_vector2(1:n),temp_variable)
    Node_Elements_3D(i)%row(1:n) = tem_vector2(1:n)
    num_Node_Elements(i) = temp_variable
    deallocate(tem_vector1,tem_vector2)
enddo
!$omp end parallel do


print *, '    Checking num_Node_Elements (<50)...'
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)     
do i=1,Num_Node
  ! If the number of elements around a node exceeds 50, a warning will pop up (under normal
  ! circumstances, this is unlikely to happen). 2023-06-14.
  if(num_Node_Elements(i)> 50)then
    print *, '    Warning-2023061401 :: num_Node_Elements(i)> 50 in Read_Geo_3D.f!'
    print *, '    Node number:', i
  endif
end do         
!$omp end parallel do  


!---------------------------------------------------------------------------
!
! STEP 10.1: Maximum length of elements around each node. IMPROV2023022203.
!
!---------------------------------------------------------------------------
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)  
do i=1,Num_Node 
    !Node_Max_L(i) = maxval(Elem_Max_L(Node_Elements_3D(i,1:num_Node_Elements(i))))
    
#ifndef Silverfrost
    Node_Max_L(i) = maxval(Elem_Max_L(Node_Elements_3D(i)%row(1:num_Node_Elements(i))))
#endif
#ifdef Silverfrost  
    ! For the Silverfrost compiler, maxval only applies to the entire array and is not suitable for
    ! obtaining the extreme values of parts of the array. 2024-09-11.
    allocate(tem_vector_real(num_Node_Elements(i)))
    tem_vector_real = Elem_Max_L(Node_Elements_3D(i)%row(1:num_Node_Elements(i)))
    Node_Max_L(i) = maxval(tem_vector_real)
    deallocate(tem_vector_real)
#endif
enddo
!$omp end parallel do  

!----------------------------------------
!
! STEP 11: The 12 edges of each element.
!
!----------------------------------------
print *, "    Finding edges of elements..." 
IF(ALLOCATED(Element_Edges)) DEALLOCATE(Element_Edges) 
ALLOCATE(Element_Edges(12,2,Num_Elem))
Element_Edges =0 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)         
do i=1,Num_Elem
  Element_Edges(1,1,i) = Elem_Node(i,1)
  Element_Edges(1,2,i) = Elem_Node(i,2)
  Element_Edges(2,1,i) = Elem_Node(i,2)
  Element_Edges(2,2,i) = Elem_Node(i,3)
  Element_Edges(3,1,i) = Elem_Node(i,3)
  Element_Edges(3,2,i) = Elem_Node(i,4)
  Element_Edges(4,1,i) = Elem_Node(i,4)
  Element_Edges(4,2,i) = Elem_Node(i,1)
  !----
  Element_Edges(5,1,i) = Elem_Node(i,1)
  Element_Edges(5,2,i) = Elem_Node(i,5)
  Element_Edges(6,1,i) = Elem_Node(i,2)
  Element_Edges(6,2,i) = Elem_Node(i,6)
  Element_Edges(7,1,i) = Elem_Node(i,3)
  Element_Edges(7,2,i) = Elem_Node(i,7)
  Element_Edges(8,1,i) = Elem_Node(i,4)
  Element_Edges(8,2,i) = Elem_Node(i,8)
  !----
  Element_Edges(9,1,i)  = Elem_Node(i,5)
  Element_Edges(9,2,i)  = Elem_Node(i,6)
  Element_Edges(10,1,i) = Elem_Node(i,6)
  Element_Edges(10,2,i) = Elem_Node(i,7)
  Element_Edges(11,1,i) = Elem_Node(i,7)
  Element_Edges(11,2,i) = Elem_Node(i,8)
  Element_Edges(12,1,i) = Elem_Node(i,8)
  Element_Edges(12,2,i) = Elem_Node(i,5)
end do 
!$omp end parallel do    

!-------------------------------------------------------------------------------------------------
!   
! STEP 12: If it is a block model, determine the number of boundary cells and the boundary cells,
! 2022-04-13
!
!-------------------------------------------------------------------------------------------------
if(Key_Block_Model == 1)then 
    Num_Elem_Block_Bou = 0
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E) schedule(static)   &
    !$OMP reduction(+:Num_Elem_Block_Bou)   
    do i_E = 1,Num_Elem
    if ((abs(x_max_Elements(i_E)  - Max_X_Coor) <= Tol_20 ) .or.    &
            (abs(x_min_Elements(i_E)  - Min_X_Coor) <= Tol_20 ) .or.  & 
             (abs(y_max_Elements(i_E)  - Max_Y_Coor) <= Tol_20 ) .or. &
             (abs(y_min_Elements(i_E)  - Min_Y_Coor) <= Tol_20 ) .or. &  
             (abs(z_max_Elements(i_E)  - Max_Z_Coor) <= Tol_20 ) .or. &
             (abs(z_min_Elements(i_E)  - Min_Z_Coor) <= Tol_20 ))then 
          Num_Elem_Block_Bou = Num_Elem_Block_Bou +1
    endif
    enddo
    !$OMP END PARALLEL DO

    ALLOCATE(Elems_Block_Bou(Num_Elem_Block_Bou))  
    i_count = 0
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E) schedule(static)
    do i_E = 1,Num_Elem
    if ((abs(x_max_Elements(i_E)  - Max_X_Coor) <= Tol_20) .or.      &
             (abs(x_min_Elements(i_E)  - Min_X_Coor) <= Tol_20) .or.   & 
             (abs(y_max_Elements(i_E)  - Max_Y_Coor) <= Tol_20) .or.   &
             (abs(y_min_Elements(i_E)  - Min_Y_Coor) <= Tol_20) .or.   &
             (abs(z_max_Elements(i_E)  - Max_Z_Coor) <= Tol_20) .or.   &
             (abs(z_min_Elements(i_E)  - Min_Z_Coor) <= Tol_20))then 
          !$OMP Critical
          i_count = i_count +1
          Elems_Block_Bou(i_count) = i_E
          !$OMP end Critical
    endif
    enddo    
    !$OMP END PARALLEL DO        
endif


!------------------------------------------------------------------------------------------------
!    
! STEP 13: Determine the outer boundaries and external surfaces of all 3D elements of the model.
!
!------------------------------------------------------------------------------------------------
print *, "    Get *.outl and *.outa files (preparing)..." 
!********************************
!                              *
!   Key_Block_Model == 0       *
!                              *
!********************************
if(Key_Block_Model == 0)then
  !IMPROV2022100202. 2022-10-02. 
  ! First, find the boundary cells.
  Num_Elem_Block_Bou = 0
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E ,c_NN,j,c_Node) schedule(static)   &
  !$OMP reduction(+:Num_Elem_Block_Bou)    
  do i_E =1,Num_Elem  
    c_NN = G_NN(1:8,i_E)
    do j=1,8 
        c_Node = c_NN(j)
        ! If a node in the current cell has fewer than 7 neighboring cells, the current cell may be a
        ! boundary cell.
        if (num_Node_Elements(c_Node)<8)then
            Num_Elem_Block_Bou = Num_Elem_Block_Bou + 1
            exit
        endif
    enddo
  enddo
  !$OMP END PARALLEL DO
      
  ALLOCATE(Elems_Block_Bou(Num_Elem_Block_Bou))  
  
  i_count = 0
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,c_NN,j,c_Node) schedule(static) 
  do i_E =1,Num_Elem  
    c_NN    = G_NN(1:8,i_E)
    do j=1,8 
        c_Node = c_NN(j)
        ! If a node in the current cell has fewer than 7 neighboring cells, the current cell may be a
        ! boundary cell.
        if (num_Node_Elements(c_Node)<8)then
            !$OMP Critical
            i_count = i_count +1
            Elems_Block_Bou(i_count) = i_E   
            !$OMP end Critical
            exit
        endif
    enddo
  enddo
  !$OMP END PARALLEL DO   

  ! Get the number of edges that meet the condition: tem_num_edges
  Num_CkElem = Num_Elem_Block_Bou
  tem_num_edges = 0
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,c_E,i_Edge,c_node_1,c_node_2,&
  !$OMP                      c_node_1_num_Ele,c_node_2_num_Ele) &
  !$OMP reduction(+:tem_num_edges)  
  do i_E = 1,Num_CkElem
      c_E = Elems_Block_Bou(i_E)
      do i_Edge = 1,12
          c_node_1 = Elem_Node(c_E,Ele_3D_Edges_Node(i_Edge,1)) 
          c_node_2 = Elem_Node(c_E,Ele_3D_Edges_Node(i_Edge,2)) 
          ! Number of elements around c_node_1.
          c_node_1_num_Ele = num_Node_Elements(c_node_1)
          ! The number of elements around c_node_2.
          c_node_2_num_Ele = num_Node_Elements(c_node_2)
          ! Each of the two nodes on the boundary edge has no more than 4 adjacent elements. 2022-10-02.
          if(c_node_1_num_Ele<=4 .and. c_node_2_num_Ele <=4) then
              !!$OMP CRITICAL
              tem_num_edges = tem_num_edges +1
              !!$OMP END CRITICAL
          endif
      enddo
  enddo  
  !$OMP END PARALLEL DO   
  
  ! Create a temporary array based on the number of edges counted
  ALLOCATE(tem1(tem_num_edges,2))
  ALLOCATE(tem2(tem_num_edges))
  
  ! Store the edges of the element into the array tem1
  Num_CkElem = Num_Elem_Block_Bou
  tem_num_edges = 0
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,c_E,i_Edge,c_node_1,c_node_2,&
  !$OMP                      c_node_1_num_Ele,c_node_2_num_Ele) &
  !$OMP  schedule(static)   
  do i_E = 1,Num_CkElem
      c_E = Elems_Block_Bou(i_E)
      do i_Edge = 1,12
          c_node_1 = Elem_Node(c_E,Ele_3D_Edges_Node(i_Edge,1)) 
          c_node_2 = Elem_Node(c_E,Ele_3D_Edges_Node(i_Edge,2)) 
          ! Number of elements around c_node_1.
          c_node_1_num_Ele = num_Node_Elements(c_node_1)
          ! The number of elements around c_node_2.
          c_node_2_num_Ele = num_Node_Elements(c_node_2)
          ! Each of the two nodes on the boundary edge has no more than 4 adjacent elements. 2022-10-02.
          if(c_node_1_num_Ele<=4 .and. c_node_2_num_Ele <=4) then
              !$OMP CRITICAL
              tem_num_edges = tem_num_edges +1
              ! Store the edges of the element into an array
              tem1(tem_num_edges,1:2) = [c_node_1,c_node_2]   
              !$OMP END CRITICAL
          endif
      enddo
  enddo    
  !$OMP END PARALLEL DO   
  
!********************************
!                              *
!   Key_Block_Model == 1       *
!                              *
!********************************
elseif(Key_Block_Model == 1)then
  Num_CkElem = Num_Elem_Block_Bou
  !................................................................................
  ! Count the number of cell edges located on the boundary faces of the cube model
  !................................................................................
  i_count = 0
  
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,c_E,i_Edge,c_node_1, &
  !$OMP                             c_node_2,c_1_x,c_1_y,c_1_z, &
  !$OMP                              c_2_x,c_2_y,c_2_z) &
  !$OMP reduction(+:i_count)  
  do i_E = 1,Num_CkElem
      c_E = Elems_Block_Bou(i_E)
      do i_Edge = 1,12
          c_node_1 = Elem_Node(c_E,Ele_3D_Edges_Node(i_Edge,1)) 
          c_node_2 = Elem_Node(c_E,Ele_3D_Edges_Node(i_Edge,2)) 
          c_1_x = Coor(c_node_1,1)
          c_1_y = Coor(c_node_1,2)
          c_1_z = Coor(c_node_1,3)
          c_2_x = Coor(c_node_2,1)
          c_2_y = Coor(c_node_2,2)
          c_2_z = Coor(c_node_2,3)
          if ((abs(c_1_x-Max_X_Coor)<=Tol_20 .and.            &
                 abs(c_2_x-Max_X_Coor)<=Tol_20)       .or.    &
                 (abs(c_1_x-Min_X_Coor)<=Tol_20 .and.         &
                  abs(c_2_x-Min_X_Coor)<=Tol_20)       .or.   &    
                 (abs(c_1_y-Max_Y_Coor)<=Tol_20 .and.         &  
                  abs(c_2_y-Max_Y_Coor)<=Tol_20)       .or.   &
                 (abs(c_1_y-Min_Y_Coor)<=Tol_20 .and.         &
                  abs(c_2_y-Min_Y_Coor)<=Tol_20)       .or.   &  
                 (abs(c_1_z-Max_Z_Coor)<=Tol_20 .and.         &
                  abs(c_2_z-Max_Z_Coor)<=Tol_20)       .or.   &
                 (abs(c_1_z-Min_Z_Coor)<=Tol_20 .and.         &
                  abs(c_2_z-Min_Z_Coor)<=Tol_20)    )then 
                  i_count = i_count +1
          endif                  
      enddo
  enddo
  !$OMP END PARALLEL DO            
  
  !...............................................................
  ! Create a temporary array based on the counted number of edges
  !...............................................................
  ALLOCATE(tem1(i_count,2))
  ALLOCATE(tem2(i_count))
  !....................................................
  ! Store the edges of the element into the array tem1
  !....................................................
  tem_num_edges = 0
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E,c_E,i_Edge,c_node_1, &
  !$OMP                             c_node_2,c_1_x,c_1_y,c_1_z, &
  !$OMP                              c_2_x,c_2_y,c_2_z)   
  do i_E = 1,Num_CkElem
      c_E = Elems_Block_Bou(i_E)
      do i_Edge = 1,12
          c_node_1 = Elem_Node(c_E,Ele_3D_Edges_Node(i_Edge,1)) 
          c_node_2 = Elem_Node(c_E,Ele_3D_Edges_Node(i_Edge,2)) 
          c_1_x = Coor(c_node_1,1)
          c_1_y = Coor(c_node_1,2)
          c_1_z = Coor(c_node_1,3)
          c_2_x = Coor(c_node_2,1)
          c_2_y = Coor(c_node_2,2)
          c_2_z = Coor(c_node_2,3)
          if ((abs(c_1_x-Max_X_Coor)<=Tol_20 .and.            &
                  abs(c_2_x-Max_X_Coor)<=Tol_20)       .or.   &
                 (abs(c_1_x-Min_X_Coor)<=Tol_20 .and.         &
                  abs(c_2_x-Min_X_Coor)<=Tol_20)       .or.   &    
                 (abs(c_1_y-Max_Y_Coor)<=Tol_20 .and.         &
                  abs(c_2_y-Max_Y_Coor)<=Tol_20)       .or.   &
                 (abs(c_1_y-Min_Y_Coor)<=Tol_20 .and.         &
                  abs(c_2_y-Min_Y_Coor)<=Tol_20)       .or.   &   
                 (abs(c_1_z-Max_Z_Coor)<=Tol_20 .and.         &
                  abs(c_2_z-Max_Z_Coor)<=Tol_20)       .or.   &
                 (abs(c_1_z-Min_Z_Coor)<=Tol_20 .and.         &
                  abs(c_2_z-Min_Z_Coor)<=Tol_20)    )then 
                  !$OMP CRITICAL
                  tem_num_edges = tem_num_edges +1
                  ! Store the edges of the element into an array
                  tem1(tem_num_edges,1:2) = [c_node_1,c_node_2]
                  !$OMP END CRITICAL
          endif                  
      enddo
  enddo
  !$OMP END PARALLEL DO    
endif





print *, "    Get *.outl and *.outa files (sorting)..." 
! Sort
call Matrix_Sort_Int(tem_num_edges,2,tem1)

!call Tool_Get_Current_Time(current_data,date_time,S_time)
print *, "    Get *.outl and *.outa files (counting)..." 
call Matrix_Count_Row_Int_m_x_2(tem_num_edges,tem1,tem2,all_num_Outline,1)


! Outer Boundary Extraction
! Extract edges that appear only once, including both outer and inner edges
print *, "    Get *.outl and *.outa files (extarcting outl)..." 
ALLOCATE(All_Outline(all_num_Outline,2))
j=0
do i=1,tem_num_edges
  if (tem2(i).eq.1)then
      j=j+1
      All_Outline(j,1:2) = tem1(i,1:2)
  end if
end do
num_Outline = j
IF(ALLOCATED(Outline)) DEALLOCATE(Outline)
ALLOCATE(Outline(num_Outline,2))
do i=1,num_Outline
  Outline(i,1:2) = All_Outline(i,1:2) 
end do
DEALLOCATE(All_Outline)

! Save 3D model outer boundary file outl
if(Key_Save_Nothing/= 1) then
  print *, "    Saving outl file..." 
  c_File_name_1   =  trim(Full_Pathname)//'.outl'
  if (Key_Data_Format==1) then
      open(101,file=c_File_name_1,status='unknown')  
          do i=1,num_Outline
              write(101, '(2I10)') (Outline(i,1:2))
          end do
      close(101)
  elseif(Key_Data_Format==2)then
      open(101,file=c_File_name_1,status='unknown',form='unformatted',access='stream')  
          do i=1,num_Outline
              write(101) (Outline(i,1:2))
          end do
      close(101)          
  endif        
endif

! Exterior Surface Extraction
! Extract edges that appear only twice, including both outer and inner edges
print *, "    Get *.outl and *.outa files (extarcting outa)..." 
ALLOCATE(All_Outline(Num_CkElem*12,2))
All_Outline(1:Num_CkElem*12,1:2) = 0
j=0
do i=1,tem_num_edges
  if (tem2(i).eq.2)then
      j=j+1
      All_Outline(j,1:2) = tem1(i,1:2)
  end if
end do
num_Outline = j
IF(ALLOCATED(OutArea)) DEALLOCATE(OutArea)
ALLOCATE(OutArea(num_Outline,2))
do i=1,num_Outline
  OutArea(i,1:2) = All_Outline(i,1:2) 
end do
DEALLOCATE(All_Outline)
! Save 3D model outer surface file outa
if(Key_Save_Nothing/= 1) then
  print *, "    Saving outa file..." 
  c_File_name_1   =  trim(Full_Pathname)//'.outa'
  if (Key_Data_Format==1) then
      open(101,file=c_File_name_1,status='unknown')  
      do i=1,num_Outline
          write(101, '(2I10)') (OutArea(i,1:2))
      end do
      close(101)   
  elseif(Key_Data_Format==2)then
      open(101,file=c_File_name_1,status='unknown',form='unformatted',access='stream') 
      do i=1,num_Outline
          write(101) (OutArea(i,1:2))
      end do    
      close(101)
  endif 
endif


if(allocated(tem1)) DEALLOCATE(tem1)
if(allocated(tem2)) DEALLOCATE(tem2)     

! External Surface Node Number Extraction (2021-09-08)
! First, convert OutArea(:,2) into a vector and store it in the temporary variable tem3.
ALLOCATE(tem3(2*num_Outline))
ALLOCATE(tem4(2*num_Outline))
tem3(1:num_Outline)               = OutArea(1:num_Outline,1)
tem3(num_Outline+1:num_Outline*2) = OutArea(1:num_Outline,2)
! Delete duplicate node numbers
call Vector_Unique_Int(2*num_Outline,2*num_Outline,tem3(1:2*num_Outline),tem4,Num_Surface_Nodes)   
! Allocate memory for global variables
if (.not. ALLOCATED(Surface_Nodes)) then
  ALLOCATE(Surface_Nodes(Num_Surface_Nodes))
endif
Surface_Nodes(1:Num_Surface_Nodes) = 0
Surface_Nodes(1:Num_Surface_Nodes) =tem4(1:Num_Surface_Nodes)
! Sorting (Quick Sort)
call Vector_QuickSort_Int(Surface_Nodes,1,Num_Surface_Nodes)
DEALLOCATE(tem3)  
DEALLOCATE(tem4)  
! Save the 3D model exterior surface node number file outn
if(Key_Save_Nothing/= 1) then
  print *, "    Saving outn file..." 
  c_File_name_1   =  trim(Full_Pathname)//'.outn'
  if (Key_Data_Format==1) then
      open(101,file=c_File_name_1,status='unknown')  
          do i=1,Num_Surface_Nodes
              write(101, '(I10)') (Surface_Nodes(i))
          end do
      close(101) 
  elseif(Key_Data_Format==2)then
      open(101,file=c_File_name_1,status='unknown',form='unformatted',access='stream')  
          do i=1,Num_Surface_Nodes
              write(101) (Surface_Nodes(i))
          end do
      close(101)
  endif   
endif

!----------------------------------------------------------------------------------------------------
! STEP 15: Divide the 3D model into 8 regions according to the coordinates and save the cell numbers
! of the 8 regions.
! Convenient for finding the element number based on coordinates.
! For details on the region division, see: PhiPsi Record Notebook, v1, P84.
! It is necessary to ensure that there is overlap between the different areas.
!IMPROV2022112101. 
!----------------------------------------------------------------------------------------------------
print *, "    Setting elements domains..."    
call D3_Set_Elements_Domains
 
!----------------------------------------
!
! STEP 16: Calculate the half bandwidth.
!
!----------------------------------------

!----------------------------------------------------------------------
!
! STEP 17: Determine the dimensions for the relevant public variables.
!
!----------------------------------------------------------------------
!2022-06-12.
!Max_Num_Crack = Max_Num_Cr_3D
IF(ALLOCATED(num_GP_Elem)) DEALLOCATE(num_GP_Elem)    
ALLOCATE(num_GP_Elem(Num_Elem))
IF(ALLOCATED(Ele_GP_Start_Num)) DEALLOCATE(Ele_GP_Start_Num)    
ALLOCATE(Ele_GP_Start_Num(Num_Elem))

!------------------------------------------------------------------------------------------------
!   
! STEP 18: Generate the variable object array related to discrete fracture surfaces. 2022-09-04.
! IMPROV2022090401.
!
!------------------------------------------------------------------------------------------------
allocate(Crack3D_Meshed_Node(Max_Num_Cr_3D))
allocate(Crack3D_Meshed_Ele(Max_Num_Cr_3D))
allocate(Crack3D_Meshed_Node_Value(Max_Num_Cr_3D))
allocate(Cr3D_Meshed_Node_in_Ele_Num(Max_Num_Cr_3D))
allocate(Cr3D_Meshed_Node_in_Ele_Local(Max_Num_Cr_3D))
allocate(Crack3D_Meshed_Ele_Attri(Max_Num_Cr_3D))
allocate(Crack3D_Meshed_Ele_Nor_Vector(Max_Num_Cr_3D))
allocate(Crack3D_Meshed_Node_Nor_Vector(Max_Num_Cr_3D))
allocate(Crack3D_Meshed_Outline(Max_Num_Cr_3D))
allocate(Crack3D_Meshed_Vertex_x_Vector(Max_Num_Cr_3D))
allocate(Crack3D_Meshed_Vertex_y_Vector(Max_Num_Cr_3D))
allocate(Crack3D_Meshed_Vertex_z_Vector(Max_Num_Cr_3D))
allocate(Crack3D_Meshed_Vertex_T_Matrx(Max_Num_Cr_3D))
allocate(Crack3D_Meshed_Vertex_T_Theta(Max_Num_Cr_3D))
allocate(Crack3D_Vector_S1(Max_Num_Cr_3D))
allocate(Cracks_FluidEle_CalP_3D(Max_Num_Cr_3D))
allocate(Cracks_FluidEle_Glo_CalP_3D(Max_Num_Cr_3D))
allocate(Cracks_FluidEle_num_CalP_3D(Max_Num_Cr_3D))
allocate(Cracks_FluidEle_EleNum_3D(Max_Num_Cr_3D))
allocate(Cracks_FluidEle_Area_3D(Max_Num_Cr_3D))
allocate(Cracks_FluidEle_Centroid_3D(Max_Num_Cr_3D))
allocate(Cracks_FluidEle_LCS_x_3D(Max_Num_Cr_3D))
allocate(Cracks_FluidEle_LCS_y_3D(Max_Num_Cr_3D))
allocate(Cracks_FluidEle_LCS_z_3D(Max_Num_Cr_3D))
allocate(Cracks_FluidEle_LCS_T_3D(Max_Num_Cr_3D))
allocate(Cracks_FluidEle_Vector_3D(Max_Num_Cr_3D))
allocate(Cracks_FluidEle_Aper_3D(Max_Num_Cr_3D))
allocate(Cracks_CalP_Coors_3D(Max_Num_Cr_3D))
allocate(Cracks_CalP_Orient_3D(Max_Num_Cr_3D))
allocate(Cracks_CalP_MeshedEl_3D(Max_Num_Cr_3D))
allocate(Cracks_CalP_Elem_3D(Max_Num_Cr_3D))
allocate(Cracks_CalP_Aper_3D(Max_Num_Cr_3D))
allocate(Cracks_CalP_UpDis_3D(Max_Num_Cr_3D))
allocate(Cracks_CalP_LowDis_3D(Max_Num_Cr_3D))
allocate(Cracks_CalP_Pres_3D(Max_Num_Cr_3D))
allocate(Cracks_CalP_Tractions_3D(Max_Num_Cr_3D))
allocate(Cracks_CalP_Pgra_3D(Max_Num_Cr_3D))
allocate(Cracks_CalP_Velo_3D(Max_Num_Cr_3D))
allocate(Cracks_CalP_Quan_3D(Max_Num_Cr_3D))
allocate(Cracks_CalP_Conc_3D(Max_Num_Cr_3D))
allocate(Cracks_CalP_Remo_Strs_3D(Max_Num_Cr_3D))

!----------------------------------------------------------------------------------
!
! STEP 19: Obtain the element list based on the material number. NEWFTU2023012401.
!
!----------------------------------------------------------------------------------
if (Key_XA  /= 2) then 
    print *, "    Getting element mat list..."    
    allocate(List_Elements_Mat(num_of_Material))
    do i_Mat =1,num_of_Material
        Elements_Num_Mat(i_Mat)  = count(Elem_Mat(1:Num_Elem)==i_Mat)
        allocate(List_Elements_Mat(i_Mat)%row(Elements_Num_Mat(i_Mat)))
    enddo
    Elements_Num_Mat(1:num_of_Material) = 0
    do i_E =1,Num_Elem  
        c_Mat_Type = Elem_Mat(i_E)
        Elements_Num_Mat(c_Mat_Type) = Elements_Num_Mat(c_Mat_Type) +1
        List_Elements_Mat(c_Mat_Type)%row(Elements_Num_Mat(c_Mat_Type)) = i_E
    enddo
endif

!-------------------------------------------------------------------------------------------------
!
! STEP 20: Assign the initial temperature and the current temperature of the element. 2023-03-13.
! IMPROV2023031301.
!
!-------------------------------------------------------------------------------------------------
if (Key_Thermal_Stress==1) then 
    print *, "    Setting element temperature..."  
    IF(ALLOCATED(Elem_Initial_T)) DEALLOCATE(Elem_Initial_T)
    ALLOCATE(Elem_Initial_T(Num_Elem))
    IF(ALLOCATED(Elem_Current_T)) DEALLOCATE(Elem_Current_T)
    ALLOCATE(Elem_Current_T(Num_Elem))
    IF(ALLOCATED(Elem_T_for_Stress)) DEALLOCATE(Elem_T_for_Stress)
    ALLOCATE(Elem_T_for_Stress(Num_Elem))
    Elem_Initial_T(1:Num_Elem)    = ZR
    Elem_Current_T(1:Num_Elem)    = ZR
    Elem_T_for_Stress(1:Num_Elem) = ZR
    
    ! If XA calls the Fortran lib, then exit.
    if(Key_Cpp_Call_Fortran_Lib==1) goto 2001
    
    ! Set according to the material number.
    if(Key_Initial_Temperature==1)then
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i_E) 
        do i_E =1,Num_Elem  
            Elem_Initial_T(i_E) = Thermal_Str_Temper( Elem_Mat(i_E))
        enddo
        !$omp end parallel do  
        Elem_Current_T = Elem_Initial_T
    ! Read from an external file.
    elseif(Key_Initial_Temperature==2)then
        !TOBEDN2023031301.
    endif
    ! Elem_T_for_Stress, the actual temperature of the element involved in the calculation.
    if(Key_Scheme_Thermal_Stress == 1) then
        Elem_T_for_Stress = Elem_Current_T
    elseif(Key_Scheme_Thermal_Stress == 2) then
        Elem_T_for_Stress = Elem_Current_T - Elem_Initial_T
    endif
    
    2001 continue 
endif

!-----------------------------------------------------------------------------------------
!
! STEP 21: Initial pore pressure, current pore pressure, Biot coefficient of the element.
! 2023-03-19. IMPROV2023031901.
!
!-----------------------------------------------------------------------------------------
if (Key_PoreP==1) then 
    print *, "    Setting element pore pressure..."  
    IF(ALLOCATED(Elem_Initial_PoreP)) DEALLOCATE(Elem_Initial_PoreP)
    ALLOCATE(Elem_Initial_PoreP(Num_Elem))
    IF(ALLOCATED(Elem_Current_PoreP)) DEALLOCATE(Elem_Current_PoreP)
    ALLOCATE(Elem_Current_PoreP(Num_Elem))
    IF(ALLOCATED(Elem_Biots)) DEALLOCATE(Elem_Biots)
    ALLOCATE(Elem_Biots(Num_Elem))    
    Elem_Initial_PoreP(1:Num_Elem)    = ZR
    Elem_Current_PoreP(1:Num_Elem)    = ZR
    Elem_Biots(1:Num_Elem)            = ZR
    
    Elem_Initial_PoreP = Initial_PoreP
    Elem_Current_PoreP = Elem_Initial_PoreP
    Elem_Biots         = Initial_Biot_Coeff
endif

!------------------------------------------------------------------------------------------
!
! STEP 22: The element's elastic modulus, Poisson's ratio, specific heat coefficient, etc.
! (activated only when Key_XA=2). 2023-03-19. IMPROV2023031903.
!         Elem_E_XA(:),Elem_Mu_XA(:),Elem_TEC_XA(:),Elem_KIc_XA(:),Elem_St_XA(:)
! D matrix of the element: Elem_D_XA(Num_Elem,6,6)
!
!------------------------------------------------------------------------------------------
if (Key_XA  == 2) then 
    print *, "    Setting element values for Key_XA = 2..."  
    IF(ALLOCATED(Elem_E_XA)) DEALLOCATE(Elem_E_XA)
    ALLOCATE(Elem_E_XA(Num_Elem))
    IF(ALLOCATED(Elem_Mu_XA)) DEALLOCATE(Elem_Mu_XA)
    ALLOCATE(Elem_Mu_XA(Num_Elem))
    IF(ALLOCATED(Elem_TEC_XA)) DEALLOCATE(Elem_TEC_XA)
    ALLOCATE(Elem_TEC_XA(Num_Elem))
    IF(ALLOCATED(Elem_KIc_XA)) DEALLOCATE(Elem_KIc_XA)
    ALLOCATE(Elem_KIc_XA(Num_Elem))
    IF(ALLOCATED(Elem_St_XA)) DEALLOCATE(Elem_St_XA)
    ALLOCATE(Elem_St_XA(Num_Elem))
    
    ! Given the initial values. These are basically set arbitrarily. Here, they are temporarily set to
    ! the parameters of Material 1 (which may not actually be provided).
    ! In fact, the relevant parameters are passed in by C.
    Elem_E_XA = Material_Para(1,1)
    Elem_Mu_XA = Material_Para(1,2)
    Elem_St_XA = Material_Para(1,5)
    Elem_KIc_XA = Material_Para(1,6)
    Elem_TEC_XA = Material_Para(1,8)
    
    !IMPROV2023031904.
    IF(ALLOCATED(Elem_D_XA)) DEALLOCATE(Elem_D_XA)
    ALLOCATE(Elem_D_XA(Num_Elem,6,6))
    Elem_D_XA = ZR
endif

!---------------------------------------------------------------------------------------
!
! STEP 23: Allocate memory for the initial stress variables of the element. 2023-03-24.
! NEWFTU2023032401.
!
!---------------------------------------------------------------------------------------
if(Key_Cpp_Call_Fortran_Lib==1) then
    print *, "    Allocate memory for Key_Cpp_Call_Fortran_Lib =1..." 
    IF(ALLOCATED(XA_Ele_InSitu_S1_Vector)) DEALLOCATE(XA_Ele_InSitu_S1_Vector)
    ALLOCATE(XA_Ele_InSitu_S1_Vector(Num_Elem,3))
    IF(ALLOCATED(XA_Ele_InSitu_S2_Vector)) DEALLOCATE(XA_Ele_InSitu_S2_Vector)
    ALLOCATE(XA_Ele_InSitu_S2_Vector(Num_Elem,3))
    IF(ALLOCATED(XA_Ele_InSitu_S3_Vector)) DEALLOCATE(XA_Ele_InSitu_S3_Vector)
    ALLOCATE(XA_Ele_InSitu_S3_Vector(Num_Elem,3))
    IF(ALLOCATED(XA_Ele_InSitu_S1_S2_S3)) DEALLOCATE(XA_Ele_InSitu_S1_S2_S3)
    ALLOCATE(XA_Ele_InSitu_S1_S2_S3(Num_Elem,3))
    XA_Ele_InSitu_S1_Vector = ZR
    XA_Ele_InSitu_S2_Vector = ZR
    XA_Ele_InSitu_S3_Vector = ZR
    XA_Ele_InSitu_S1_S2_S3  = ZR
    IF(ALLOCATED(XA_Ele_Stress)) DEALLOCATE(XA_Ele_Stress)
    ALLOCATE(XA_Ele_Stress(Num_Elem,6))
endif

RETURN
END SUBROUTINE Read_Geo_3D