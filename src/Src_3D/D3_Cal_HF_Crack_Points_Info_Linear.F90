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
 
subroutine D3_Cal_HF_Crack_Points_Info_Linear(isub)
! Calculate the fluid element calculation point numbers for cracks, the coordinates 
! and directions of the calculation points, the crack segment numbers they belong to,
! and the element numbers they are in, for hydraulic fracturing linear HF elements.
! In addition, this subroutine also needs to specify the type of calculation point.
! The keyword Key_Kink_Point is used to control whether the kink points of crack
! segments are used as calculation points.
!
! Fix bug on 2020-01-01.
!

! The core variables involved are explained as follows:
! integer Cracks_FluidEle_CalP_3D(num_Crack, Max_N_CalP_3D, 5) !Fluid element calculation point
! number for each crack
! integer Cracks_FluidEle_Glo_CalP_3D(num_Crack, Max_N_CalP_3D, 5) !Global numbering of fluid
! element calculation points for each crack (used for assembling the Q matrix)
! integer Cracks_FluidEle_num_CalP_3D(num_Crack, Max_N_CalP_3D) !Number of calculation points for
! each fluid element in each crack
! integer Cracks_FluidEle_num_3D(num_Crack)                             !Number of fluid elements for each crack
! real(kind=FT) Cracks_CalP_Coors_3D(num_Crack, Max_N_CalP_3D, 3) !Calculate the coordinates of
! points for each crack (some of which do not participate in the calculation)
! integer Cracks_CalP_Num_3D(num_Crack) !Number of calculation points for each crack (some do not
! participate in the calculation)
! integer Cracks_Real_CalP_Num_3D(1:num_Crack) !Number of calculation points for each crack
! (participating in the calculation)
! integer Cracks_Real_CalPs_3D(num_Crack, Max_N_CalP_3D) !Number of calculation points for each
! crack (participating in the calculation)
! real(kind=FT) Cracks_CalP_Orient_3D(num_Crack, Max_N_CalP_3D, 3) !Calculate the orientation of
! each crack at calculation points (fluid nodes) (element outward normal vector of the crack
! surface)
! integer Cracks_FluidEle_EleNum_3D(i_C, Max_N_CalP_3D) !The solid element number corresponding to
! the fluid element
! integer Cracks_CalP_MeshedEl_3D(num_Crack, Max_N_CalP_3D) !Crack discrete element number
! corresponding to each calculation point of the crack
! integer Cracks_CalP_Elem_3D(num_Crack, Max_N_CalP_3D, 2) !Element number and edge number of each
! calculation point for each crack
! real(kind=FT) Cracks_FluidEle_Area_3D(1:num_Crack,1:Max_N_CalP_3D) !Area of fluid elements for
! each crack
! real(kind=FT) Cracks_FluidEle_Centroid_3D(1:num_Crack,1:Max_N_CalP_3D) !Centroid of fluid elements
! for each crack
! real(kind=FT) Cracks_FluidEle_Vector_3D(1:num_Crack,1:Max_N_CalP_3D,1:3) !Average outward normal
! vector of fluid elements for each crack
! real(kind=FT) Cracks_FluidEle_LCS_x_3D(1:num_Crack,1:Max_N_CalP_3D,1:3) !Local coordinate system
! x-axis of the centroid positions of fluid elements for each crack
! real(kind=FT) Cracks_FluidEle_LCS_y_3D(1:num_Crack,1:Max_N_CalP_3D,1:3) !y-axis of the local
! coordinate system for the centroid position of each crack's fluid element
! real(kind=FT) Cracks_FluidEle_LCS_z_3D(1:num_Crack,1:Max_N_CalP_3D,1:3) !Local coordinate system
! z-axis (normal vector) of the centroid position of each crack's fluid element
! real(kind=FT) Cracks_FluidEle_LCS_T_3D(1:num_Crack,1:Max_N_CalP_3D,1:3,1:3)!Transformation matrix
! of the local coordinate system at the centroid position of the fluid element for each crack

! Save the local information corresponding to the global fluid node numbers: 
! including fracture numbers, fluid element numbers, and the fluid node 
! numbers corresponding to the fluid elements. 2022-06-04.
!Cracks_FluidEle_CalP_Glo_Info(Glo_location,1) = i_C
!Cracks_FluidEle_CalP_Glo_Info(Glo_location,2) = i_FluidEle
!Cracks_FluidEle_CalP_Glo_Info(Glo_location,3) = i_CalP
      
! The in-situ stress at the globally numbered fluid nodes (perpendicular to the fracture surface), 
! that is, the initial in-situ stress in the direction normal to the fracture surface. 2022-06-04.

! Local variable
! integer temp_CalP_Near_Model_Bou(num_Crack,Max_N_CalP_3D) ! Used to mark whether fluid element
! nodes are on the model boundary (or very close to it), 2021-08-15

! Subroutines used:
! Tool_Normal_vector_of_3D_Tri(Tri_P1,Tri_P2,Tri_P3,Np) !Calculate the outer normal vector of a 3D
! triangle
! Tool_Intersection_of_AB_and_Triangle_3D(A,B,Tri_P1,Tri_P2,Tri_P3 !Calculate the intersection of a
! line segment and a 3D triangle
!                                 Yes_Inter,InterSection_P)
! Matrix_Delete_Near_Point_Dou.f !Used to delete points in the matrix that are close to each other
! (distance less than Dis_Tol)
! Tool_Area_Tri_3D(c_P1, c_P2, c_P3, delta_Fluid_Area)                  !Calculate the area of a subspace triangle
! Tool_Yes_Point_Near_Model_Boundary_3D(Point, Yes_Near_Boundary) !Used to check if a point is on
! the model boundary (or very close to it), 2021-08-15

! The global variable Crack_Type_Status_3D(i_C,10) is used to indicate the type and status of
! cracks.
! Column 1 (Fracture Type): =1, HF fracture; =2, natural fracture; =3, post-fracturing hydraulic
! fracture
! Note: Natural fractures and after-pressurized hydraulic fractures may potentially turn into HF
! fractures.
! Column 2 (Fracture Status): =1, HF fracturing not completed; =2, HF fracturing completed
! Column 3 (Can the crack continue to propagate): =1, yes; =0, no
! Column 4 (Whether the fracture has obtained a fluid node): =1, Yes; =0, No
! Column 5 (Did the crack propagate in the previous step?): =1, Yes; =0, No
! The global variable Cracks_Stages_Wellbores(i_WB, i_Stage, i_C) is used for the crack numbers
! corresponding to each segment of each wellbore.
      
!=============================
! Read public variable module
!=============================
use Global_Float_Type                                                                
use Global_Crack_Common
use Global_Crack_3D
use Global_Model
use Global_Elem_Area_Vol
use Global_Common
use Global_HF
use Global_Stress
use Global_XFEM_Elements
use omp_lib

use Global_Cal_Ele_Num_by_Coors_3D  
use Global_INTERFACE_Tool_ThetaX_ThetaY_ThetaZ_3D_rotation

!======================
! Variable Declaration
!======================
implicit none
integer,intent(in)::isub
integer i_C,i_P,i_E,i_E_Edge
integer i_Crack_Ele
integer Crack_Node1,Crack_Node2,Crack_Node3
real(kind=FT) Point1(3),Point2(3),Point3(3)
integer Edge_node1,Edge_node2
real(kind=FT) A(3),B(3)
logical Yes_Inter
real(kind=FT) InterSection_P(3),Unit_Norm_Vector(3)
integer c_num_CalP,c_Location,c_ele_num_CalP(Num_Elem)
logical Yes_Exist
real(kind=FT),ALLOCATABLE::tmp_ele_CalP(:,:,:)
                                               ! All elements that contain fluid nodes are necessarily enhanced elements, so here it is changed to
                                               ! num_XFEM_Elem, thereby significantly reducing memory usage.
                                               ! IMPROV2022081901.
real(kind=FT),ALLOCATABLE::tem_tmp_ele_CalP(:,:,:)
                                               
integer Uniqued_m,Uni_Mat_Count(150)
real(kind=FT) Matrix(150,3),Uniqued_Matrix(150,3)
real(kind=FT) Dis_Tol
integer i_CalP,c_count_FluidEle,i_FluidEle,i_sub_tri
logical c_Yes
real(kind=FT) c_P1(3),c_P2(3),c_P3(3)
integer P1_CalP_num,P2_CalP_num,P3_CalP_num
real(kind=FT) c_Fluid_Area,delta_Fluid_Area,c_Fluid_Vector(3)
integer c_CalP
real(kind=FT) c_Norm_2
integer,ALLOCATABLE::tem_real_CalP(:)
integer,ALLOCATABLE::Uniqued_Vec_tem_real_CalP(:) 
integer tem_real_num,tem_Uniqued_n
integer cc_location
logical c_Yes_In
real(kind=FT) c_Fluid_Centroid(3)
integer num_Fluid_nodes
real(kind=FT) Center_X,Center_Y,Center_Z
real(kind=FT) Point_Center(3),Line_O_1(2,3),Line_O_N(2,3)
integer Fluid_node1,i_FluidNode,Fluid_nodeN
real(kind=FT),ALLOCATABLE:: Angle_Lines(:)
real(kind=FT) vector_a(3),vector_b(3)
integer i_temp
integer Old_Cracks_FluidEle_CalP_3D(Max_ele_num_CalP)
integer fnode_1,fnode_2,fnode_3
integer c_count_FluEl
integer  c_Solid_El,i_added,c_node2,c_node3 
real(kind=FT) c_x_vector(3),c_y_vector(3),c_z_vector(3)
integer c_First_CalP
real(kind=FT) Coor_First_CalP(3)
real(kind=FT) ThetaX,ThetaY,ThetaZ,c_T_Matrx(3,3)
real(kind=FT) c_Point_A(3),c_Point_B(3)
integer i_Fnode,c_E_ndoes(3)
logical c_Yes_Exist       
logical c_Yes_Inter
real(kind=FT) c_Inter_Point_A(3),c_Inter_Point_B(3)
integer c_Vertex_1,c_Vertex_2
integer num_Inter_Outline
real(kind=FT) c_Four_Points(4,3)
integer Glo_location
real(kind=FT) ori_n(3),Normal_Stress
integer Solid_Elem
real(kind=FT) c_Insitu_Stress(6)
real(kind=FT) max_x_Tri,max_y_Tri,max_z_Tri
real(kind=FT) min_x_Tri,min_y_Tri,min_z_Tri
logical Logical_Yes_x,Logical_Yes_y,Logical_Yes_z
integer c_Cr_Location
integer Flag_Eles_Has_CalP(Num_Elem)
integer Count_Eles_Has_CalP
integer Eles_Has_CalP_Location(Num_Elem)
integer Eles_Has_CalP_List(Num_Elem)
integer Size_tmp_ele_CalP
integer Size_Poss_CalP
integer Ele_Num_Cache
integer c_Elem
integer c_Max_N_FluEl_3D
integer c_Max_N_CalP_3D
integer Eles_List_to_be_Checked(Num_Elem),number_Eles_to_be_Checked,c_E

logical Yes_A_on,Yes_B_on


#ifndef Silverfrost      
print *,'    Preparing calculation points of each crack...'
901  FORMAT(5X,'Caution :: crack ',I8,' has no enriched node!')       
      
do i_C = 1,num_Crack
    !########################################################################################### 
    ! If the fluid calculation point has already been computed, and the fracture did 
    ! not propagate in the previous step or is a wellbore HF-induced fracture,
    ! Then there is no need to calculate further. IMPROV2022060202.
    ! The global variable Crack_Type_Status_3D(i_C,10) is used to indicate the type and 
    ! status of cracks.
    ! Column 1 (Fracture Type): =1, HF Fracture; =2, Natural Fracture; 
    !                           =3, Post-Fracturing Hydraulic Fracture
    ! Column 2 (Fracture Status): =1, HF fracturing not completed; =2, HF fracturing completed
    ! Column 3 (Can the crack continue to propagate): =1, yes; =0, no
    ! Column 4 (Whether the fracture has obtained a fluid node): =1, Yes; =0, No
    ! Column 5 (Did the crack propagate in the previous step?): =1, Yes; =0, No
    !############################################################################################
    ! If the fluid calculation points have already been computed
    if(Crack_Type_Status_3D(i_C,4) == 1)then
         !========================
         !====   OPTION -1   ====
         !========================
            ! If the crack did not extend in the previous step, there is no need to update it.
            !        if(Crack_Type_Status_3D(i_C,5) == 0) then   
            ! Cracks_Real_CalP_Num_3D(i_C) = 0 !Number of fluid nodes participating in fluid-structure
            ! interaction calculation
            !            print*, 'Crack', i_C ,' not propagated last step'
            !            cycle
            !        endif 
            ! If it is a wellbore HF post-fracture, there is no need to calculate the fluid nodes, as they are
            ! no longer involved in the calculation.
            !        if(Crack_Type_Status_3D(i_C,1) == 3)then   
            ! Cracks_Real_CalP_Num_3D(i_C) = 0 !Number of fluid nodes involved in fluid-structure interaction
            ! calculation
            !            cycle
            !        endif   

        !=============================================
        !====   OPTION -2   ====   !BUGFIX2022062001. 
        !=============================================
        ! If it is a wellbore HF post-fracture, there is no need to calculate the fluid nodes, as they are
        ! no longer involved in the calculation.
        if(Crack_Type_Status_3D(i_C,1) == 3)then    
            Cracks_Real_CalP_Num_3D(i_C) = 0
            cycle
        endif              
        ! If the crack did not propagate in the previous step and is not an HF crack, there is no need to
        ! update it.
        if(Crack_Type_Status_3D(i_C,5) == 0 .and. Crack_Type_Status_3D(i_C,1)/=1) then    
            Cracks_Real_CalP_Num_3D(i_C) = 0
            cycle
        endif 
    endif
        
    ! The fracture has been marked as having fluid nodes calculated
    Crack_Type_Status_3D(i_C,4) = 1

    !#####################
    ! Data Initialization
    !#####################
    Cracks_FluidEle_num_3D(i_C) =0
    Cracks_CalP_Num_3D(i_C)                        = 0
    Cracks_Real_CalP_Num_3D(i_C)                   = 0
    
    ! If no memory has been allocated for the discrete fracture surfaces before this fracture,
    ! reallocate it. 2022-09-04.
    call D3_Allocate_Crack_Memory(i_C,2,0)
    call D3_Allocate_Crack_Memory(i_C,3,0)
    
    !..............................................................
    ! Allocate memory for the memory-intensive tmp_ele_CalP(:,:,:)
    ! The first dimension is allocated 1024 initially; if it goes 
    ! out of bounds later, it will be expanded by 1024 each time.
    ! The second dimension is initially allocated 16; if it 
    ! overflows later, it will be expanded by 16 each time.
    ! 2022-08-20.
    !..............................................................
    !Size_tmp_ele_CalP = int(dble(num_XFEM_Elem)/TWO)
    !Size_tmp_ele_CalP = int(dble(num_XFEM_Elem)/dble(num_Crack))
    Size_tmp_ele_CalP = 1024
    Size_Poss_CalP    = 16
    ALLOCATE(tmp_ele_CalP(Size_tmp_ele_CalP,Size_Poss_CalP,3))
        
    !############################################################################################
    ! Searching for the current crack calculation point:
    ! (1) Obtained through the intersection points of the element edges and discrete fracture
    ! surfaces.
    ! (2) Use the two points of the BaseLine of the crack tip enhancement element as calculation
    ! points.
    !############################################################################################
    c_num_CalP =0
    c_ele_num_CalP(1:Num_Elem) =0
    tmp_ele_CalP(1:Size_tmp_ele_CalP,1:Size_Poss_CalP,1:3) = ZR 
    Count_Eles_Has_CalP = 0
    Eles_Has_CalP_Location(1:num_Crack) = 0
    Eles_Has_CalP_List(1:Num_Elem)      = 0

    ! This section has data dependencies and will not use OpenMP. 2022-09-23.
    do i_Crack_Ele =1,Crack3D_Meshed_Ele_num(i_C)
      Crack_Node1 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,1)  
      Crack_Node2 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,2)  
      Crack_Node3 = Crack3D_Meshed_Ele(i_C)%row(i_Crack_Ele,3)        
      Point1(1:3) = Crack3D_Meshed_Node(i_C)%row(Crack_Node1,1:3)  
      Point2(1:3) = Crack3D_Meshed_Node(i_C)%row(Crack_Node2,1:3)  
      Point3(1:3) = Crack3D_Meshed_Node(i_C)%row(Crack_Node3,1:3)       

      max_x_Tri = maxval([Point1(1),Point2(1),Point3(1)])
      min_x_Tri = minval([Point1(1),Point2(1),Point3(1)])
      max_y_Tri = maxval([Point1(2),Point2(2),Point3(2)])
      min_y_Tri = minval([Point1(2),Point2(2),Point3(2)])
      max_z_Tri = maxval([Point1(3),Point2(3),Point3(3)])
      min_z_Tri = minval([Point1(3),Point2(3),Point3(3)])          
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! TRY 1: Obtain calculation points through the intersection 
      !        of element edges and discrete fracture surfaces.
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !Loop over 3D elements. 
      do i_E = 1,num_XFEM_Elem
        c_Elem = XFEM_Elem_List(i_E)
        !.............................................................................
        ! Check whether the coordinate ranges intersect (overlap). If x, y, or
        ! z do not overlap, an intersection is impossible.
        ! Element coordinate range and triangular fracture surface element coordinate
        ! range. 2022-06-12. IMPROV2022061202.
        !.............................................................................
        call Tool_Yes_Two_Ranges_Overlapped_Double([x_min_Elements(c_Elem),x_max_Elements(c_Elem)],  &
                                                   [min_x_Tri,max_x_Tri],Logical_Yes_x) 
        if(Logical_Yes_x .eqv. .False.)then
            cycle
        endif
        call Tool_Yes_Two_Ranges_Overlapped_Double([y_min_Elements(c_Elem),y_max_Elements(c_Elem)],  &
                                                   [min_y_Tri,max_y_Tri],Logical_Yes_y) 
        if(Logical_Yes_y .eqv. .False.)then
            cycle
        endif
        call Tool_Yes_Two_Ranges_Overlapped_Double([z_min_Elements(c_Elem),z_max_Elements(c_Elem)],  &
                                                   [min_z_Tri,max_z_Tri],Logical_Yes_z) 
        if(Logical_Yes_z .eqv. .False.)then
            cycle
        endif            
            
        do i_E_Edge = 1,12
          ! The variable Element_Edges(12,2,Num_Elem) stores the 12 edges (corresponding node numbers) of the
          ! cell.
          Edge_node1 = Element_Edges(i_E_Edge,1,c_Elem)
          Edge_node2 = Element_Edges(i_E_Edge,2,c_Elem)
          A(1:3) = Coor(Edge_node1,1:3)
          B(1:3) = Coor(Edge_node2,1:3)
          
          ! Calculate the intersection points between the current edge segments and the spatial triangles
          ! (triangles of discrete fracture surfaces).
          call Tool_Intersection_of_AB_and_Triangle_3D(A,B,Point1,Point2,Point3,Yes_Inter,InterSection_P)
          
          ! If there is an intersection.
          if(Yes_Inter .eqv. .True.) then
            ! Cumulative count of calculation points contained in the current element
            c_ele_num_CalP(c_Elem) =c_ele_num_CalP(c_Elem) +1
            ! Expand the array to prevent out-of-bounds errors.
            if(c_ele_num_CalP(c_Elem) > size(tmp_ele_CalP,2)) then
                  allocate(tem_tmp_ele_CalP(size(tmp_ele_CalP,1),size(tmp_ele_CalP,2)+16,3))
                  tem_tmp_ele_CalP(:,1:size(tmp_ele_CalP,2),:)=tmp_ele_CalP(:,1:size(tmp_ele_CalP,2),:)
                  deallocate(tmp_ele_CalP)
                  call move_alloc(tem_tmp_ele_CalP,tmp_ele_CalP)
            endif            
            
            ! The following is the algorithm adopted on 2022-08-20, which uses dynamic arrays to reduce the
            ! memory usage of tmp_ele_CalP.
            if (.not. any(Eles_Has_CalP_List(1:Count_Eles_Has_CalP) == c_Elem))then  
                 Count_Eles_Has_CalP = Count_Eles_Has_CalP + 1
                 Eles_Has_CalP_Location(c_Elem) = Count_Eles_Has_CalP
                 Eles_Has_CalP_List(Count_Eles_Has_CalP) =  c_Elem                
            endif
            if(Eles_Has_CalP_Location(c_Elem) > size(tmp_ele_CalP,1))then
                allocate(tem_tmp_ele_CalP(size(tmp_ele_CalP,1)+1024,size(tmp_ele_CalP,2),3))
                tem_tmp_ele_CalP(1:size(tmp_ele_CalP,1),:,:)=tmp_ele_CalP(1:size(tmp_ele_CalP,1),:,:)
                deallocate(tmp_ele_CalP)
                call move_alloc(tem_tmp_ele_CalP,tmp_ele_CalP)
            endif
            tmp_ele_CalP(Eles_Has_CalP_Location(c_Elem),c_ele_num_CalP(c_Elem),1:3)=InterSection_P
            
            !------------------------------------------
            ! Check if the intersection already exists
            !------------------------------------------
            call Vector_belongs_Matrix_Is_Dou(c_num_CalP,3,Cracks_CalP_Coors_3D(i_C)%row(1:c_num_CalP,1:3),&
                                              InterSection_P(1:3),c_Location,Yes_Exist)   

            !-----------------------------------------------------------------------------
            ! If this intersection does not exist, a new calculation point will be added.
            !-----------------------------------------------------------------------------
            if(Yes_Exist .eqv. .False.)then
             c_num_CalP =c_num_CalP +1       
             c_Max_N_CalP_3D = size(Cracks_CalP_Coors_3D(i_C)%row,1)
             if(c_num_CalP > c_Max_N_CalP_3D)then
                !call Warning_Message('S',Keywords_Blank) 
                ! Expand memory. NEWFTU2022110501.
                call D3_Allocate_Crack_Memory(i_C,3,1)
                !c_Max_N_CalP_3D = size(Cracks_CalP_Coors_3D(i_C)%row,1)
             endif 
             
             Cracks_CalP_Coors_3D(i_C)%row(c_num_CalP,1:3)=InterSection_P
             !2022-09-27.
             if(sum(abs(InterSection_P))<=Tol_10)then
                print *,'    ERROR-2022092704 :: illegal Cracks_CalP_Coors_3D in D3_Cal_HF_Crack_Points_Info_Linear.f'
                call Warning_Message('S',Keywords_Blank)  
             endif     
             ! Calculate the outward normal vector of the element
             call Tool_Normal_Unit_vector_of_3D_Tri(Point1,Point2,Point3,Unit_Norm_Vector)             
             Cracks_CalP_Orient_3D(i_C)%row(c_num_CalP,1:3) =Unit_Norm_Vector
             ! Calculate the discrete crack element number corresponding to the point
             Cracks_CalP_MeshedEl_3D(i_C)%row(c_num_CalP) = i_Crack_Ele
             ! Calculate the cell number and edge number of the point
             Cracks_CalP_Elem_3D(i_C)%row(c_num_CalP,1) = c_Elem
             Cracks_CalP_Elem_3D(i_C)%row(c_num_CalP,2) = i_E_Edge 
            endif
          endif
        enddo
      enddo
     enddo
     !!$omp end parallel do     
     
     
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! TRY 2: Use the two points of the BaseLine enhanced by the
     !        crack tip element as calculation points.
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     !=================================================================================
     ! Fixed the bug where tmp_ele_CalP(Num_Elem,20000,3) consumed a large amount 
     ! of memory (BUGFIX2022041601)
     ! (1) The previous algorithm had a bug. The following loop 'do i_E = 1, Num_Elem'
     !     should not be placed inside 'do i_Crack_Ele'.
     !     Otherwise, repeated detection will occur, causing too many 
     !     c_ele_num_CalP(i_E) = c_ele_num_CalP(i_E) + 1
     ! (2) Change tmp_ele_CalP(Num_Elem,20000,3) to tmp_ele_CalP(Num_Elem,50,3)
     !=================================================================================
     !Loop over 3D elements.
     ! This section has data dependencies and will not use OpenMP. 2022-09-23.
     do i_E = 1,num_XFEM_Elem
        ! If the current element is not an enriched element, skip it. Because fluid elements must be within
        ! XFEM elements. 2022-08-20. IMPROV2022081901.
        !if(Elem_XFEM_Flag(i_E)==0) cycle
        
        c_Elem = XFEM_Elem_List(i_E)
        
        call Vector_Location_Int_v2(Solid_El_Max_num_Crs,Solid_El_Crs(c_Elem,1:Solid_El_Max_num_Crs),i_C,c_Cr_Location)
        !c_Cr_Location =
        !minloc(Solid_El_Crs(i_E,1:Solid_El_Max_num_Crs),1,mask=(Solid_El_Crs(i_E,1:Solid_El_Max_num_Crs)==i_C))
        !!IMPROV2022081502.
        
        if(c_Cr_Location==0 .or. Solid_El_num_Crs(c_Elem)==0)then
            cycle
        endif
        
        !=====================================================================
        ! If it is a crack tip enrichment element, take the two points of its
        ! BaseLine as fluid nodes (BUGFIX2021081601)
        !=====================================================================
        !If do not consider fluid elements of crack tip element, then don't enter current do i_E
        !(2021-08-16).
        if (Key_Tip_Fluid_Element_3D==1)then
          if(sum(abs(Solid_El_Tip_BaseLine(c_Elem)%row(c_Cr_Location,1:2,1:3)))>Tol_10) then
              c_Point_A(1:3) = Solid_El_Tip_BaseLine(c_Elem)%row(c_Cr_Location,1,1:3)
              c_Point_B(1:3) = Solid_El_Tip_BaseLine(c_Elem)%row(c_Cr_Location,2,1:3)
              
              ! If c_Point_A(1:3) or c_Point_B(1:3) is invalid, skip this cell. 2022-10-23. IMPROV2022102302.
              if(sum(abs(c_Point_A(1:3)))<=Tol_10)then
                 print *,'    WARNING-2022102303 :: illegal Solid_El_Tip_BaseLine!'
                 print *,'                          In D3_Cal_HF_Crack_Points_Info_Linear.f90.'
                 print *,'                          c_Elem:        ',c_Elem
                 print *,'                          c_Cr_Location: ',c_Cr_Location
                 cycle
              endif               
              if(sum(abs(c_Point_B(1:3)))<=Tol_10)then
                 print *,'    WARNING-2022102304 :: illegal Solid_El_Tip_BaseLine!'
                 print *,'                          In D3_Cal_HF_Crack_Points_Info_Linear.f90.'
                 print *,'                          c_Elem:        ',c_Elem
                 print *,'                          c_Cr_Location: ',c_Cr_Location
                 cycle
              endif 
              
              !/////////////////////////////////
              ! The first point of the baseline
              !/////////////////////////////////
              ! Cumulative count of calculation points contained in the current element
              c_ele_num_CalP(c_Elem) =c_ele_num_CalP(c_Elem) +1
              ! Expand the array to prevent out-of-bounds errors.
              if(c_ele_num_CalP(c_Elem) > size(tmp_ele_CalP,2)) then
                  print *,'    Extending dim 2 of allocatable tmp_ele_CalP(:,1:extended,:) by 16...'
                  allocate(tem_tmp_ele_CalP(size(tmp_ele_CalP,1),size(tmp_ele_CalP,2)+16,3))
                  tem_tmp_ele_CalP(:,1:size(tmp_ele_CalP,2),:)=tmp_ele_CalP(:,1:size(tmp_ele_CalP,2),:)
                  deallocate(tmp_ele_CalP)
                  call move_alloc(tem_tmp_ele_CalP,tmp_ele_CalP)
              endif               
               
              if (.not. any(Eles_Has_CalP_List(1:Count_Eles_Has_CalP) == c_Elem))then  
                 Count_Eles_Has_CalP = Count_Eles_Has_CalP + 1
                 Eles_Has_CalP_Location(c_Elem) = Count_Eles_Has_CalP
                 Eles_Has_CalP_List(Count_Eles_Has_CalP) =  c_Elem                
              endif
              if(Eles_Has_CalP_Location(c_Elem) > size(tmp_ele_CalP,1))then
                  allocate(tem_tmp_ele_CalP(size(tmp_ele_CalP,1)+1024,size(tmp_ele_CalP,2),3))
                  tem_tmp_ele_CalP(1:size(tmp_ele_CalP,1),:,:)=tmp_ele_CalP(1:size(tmp_ele_CalP,1),:,:)
                  deallocate(tmp_ele_CalP)
                  call move_alloc(tem_tmp_ele_CalP,tmp_ele_CalP)
              endif              
              tmp_ele_CalP(Eles_Has_CalP_Location(c_Elem),c_ele_num_CalP(c_Elem),1:3)=c_Point_A
              
              ! Check if the intersection already exists
              call Vector_belongs_Matrix_Is_Dou(c_num_CalP,3,Cracks_CalP_Coors_3D(i_C)%row(1:c_num_CalP,1:3), &
                                                c_Point_A(1:3),c_Location,Yes_Exist)   
              if(Yes_Exist .eqv. .False.)then
                c_num_CalP =c_num_CalP +1   
                
                c_Max_N_CalP_3D = size(Cracks_CalP_Coors_3D(i_C)%row,1)
                if(c_num_CalP > c_Max_N_CalP_3D)then
                    !call Warning_Message('S',Keywords_Blank) 
                    ! Expand memory. NEWFTU2022110501.
                    call D3_Allocate_Crack_Memory(i_C,3,1)
                    !c_Max_N_CalP_3D = size(Cracks_CalP_Coors_3D(i_C)%row,1)
                endif 
                
                Cracks_CalP_Coors_3D(i_C)%row(c_num_CalP,1:3)=c_Point_A(1:3)
                ! Calculate the outward normal vector of the element
                Cracks_CalP_Orient_3D(i_C)%row(c_num_CalP,1:3) =&
                             Solid_El_Tip_BaseLine_Nor_Vec(c_Elem)%row(c_Cr_Location,1:3)
                ! Calculate the discrete crack element number corresponding to the point
                ! Cracks_CalP_MeshedEl_3D(i_C,c_num_CalP) = i_Crack_Ele !This variable is not used, so it is not
                ! assigned here
                ! Calculate the cell number and edge number of the point
                Cracks_CalP_Elem_3D(i_C)%row(c_num_CalP,1) = c_Elem
                ! Cracks_CalP_Elem_3D(i_C,c_num_CalP,2) = i_E_Edge !For fluid nodes formed by the BaseLine, this
                ! variable does not exist (in fact, this variable has no effect in the program)
              endif
              !//////////////////////////////////
              ! The second point of the baseline
              !//////////////////////////////////
              ! Cumulative count of calculation points contained in the current element
              c_ele_num_CalP(c_Elem) =c_ele_num_CalP(c_Elem) +1
              ! Expand the array to prevent out-of-bounds errors.
              if(c_ele_num_CalP(c_Elem) > size(tmp_ele_CalP,2)) then
                    print *,'    Extending dim 2 of allocatable tmp_ele_CalP(:,1:extended,:) by 16...'
                    allocate(tem_tmp_ele_CalP(size(tmp_ele_CalP,1),size(tmp_ele_CalP,2)+16,3))
                    tem_tmp_ele_CalP(:,1:size(tmp_ele_CalP,2),:)=tmp_ele_CalP(:,1:size(tmp_ele_CalP,2),:)
                    deallocate(tmp_ele_CalP)
                    call move_alloc(tem_tmp_ele_CalP,tmp_ele_CalP)
              endif               
              
              ! tmp_ele_CalP(i_E,c_ele_num_CalP(i_E),1:3)=c_Point_B !Temporary variable: coordinates of
              ! calculation points contained in the current element
              
              tmp_ele_CalP(Eles_Has_CalP_Location(c_Elem),c_ele_num_CalP(c_Elem),1:3)=c_Point_B
              
              ! Check if the intersection already exists
              call Vector_belongs_Matrix_Is_Dou(c_num_CalP,3,Cracks_CalP_Coors_3D(i_C)%row(1:c_num_CalP,1:3), &
                                                c_Point_B(1:3),c_Location,Yes_Exist)   
              if(Yes_Exist .eqv. .False.)then
                 c_num_CalP =c_num_CalP +1  
                 c_Max_N_CalP_3D = size(Cracks_CalP_Coors_3D(i_C)%row,1)
                 if(c_num_CalP > c_Max_N_CalP_3D)then
                    !call Warning_Message('S',Keywords_Blank) 
                    ! Expand memory. NEWFTU2022110501.
                    call D3_Allocate_Crack_Memory(i_C,3,1)
                    !c_Max_N_CalP_3D = size(Cracks_CalP_Coors_3D(i_C)%row,1)
                 endif  
                 
                 Cracks_CalP_Coors_3D(i_C)%row(c_num_CalP,1:3)=c_Point_B(1:3)               
                 ! Calculate the outward normal vector of the element
                 Cracks_CalP_Orient_3D(i_C)%row(c_num_CalP,1:3) =&
                               Solid_El_Tip_BaseLine_Nor_Vec(c_Elem)%row(c_Cr_Location,1:3)
                 ! Calculate the discrete crack element number corresponding to the point
                 ! Cracks_CalP_MeshedEl_3D(i_C,c_num_CalP) = i_Crack_Ele ! This variable has no effect, so it is not
                 ! assigned here
                 ! Calculate the cell number and edge number of the point
                 Cracks_CalP_Elem_3D(i_C)%row(c_num_CalP,1) = c_Elem
                 ! Cracks_CalP_Elem_3D(i_C, c_num_CalP, 2) = i_E_Edge !For fluid nodes formed by the BaseLine, this
                 ! variable does not exist (in fact, this variable has no effect in the program)
              endif
          endif
        endif
     enddo
     Cracks_CalP_Num_3D(i_C)   = c_num_CalP
         
     !######################################## 
     ! Fluid Element Calculation Point Number
     !########################################
     !..........................................................
     ! STEP 1: Delete duplicate fluid calculation points in the 
     !         current 3D fracture surface mesh elements
     !..........................................................
     !IMPROV2023071701.
     Eles_List_to_be_Checked(1:Num_Elem)   = 0
     number_Eles_to_be_Checked = 0
     do i_E = 1,Num_Elem 
        if(c_ele_num_CalP(i_E)>=2)then 
            number_Eles_to_be_Checked = number_Eles_to_be_Checked +1
            Eles_List_to_be_Checked(number_Eles_to_be_Checked) = i_E
        endif
     enddo
     !$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_E,c_E,Matrix,Uniqued_Matrix,Uniqued_m,Uni_Mat_Count) &
     !$OMP            SCHEDULE(STATIC)    
     ! OpenMP parallelization. IMPROV2022092302.
     !do i_E = 1,Num_Elem 
     do i_E = 1,number_Eles_to_be_Checked
       c_E = Eles_List_to_be_Checked(i_E)     
       !Matrix(1:c_ele_num_CalP(i_E),1:3) =tmp_ele_CalP(i_E,1:c_ele_num_CalP(i_E),1:3) 
       !Matrix(1:c_ele_num_CalP(i_E),1:3) =tmp_ele_CalP(Elem_Location(i_E,1),1:c_ele_num_CalP(i_E),1:3)
       !!IMPROV2022081901.
       Matrix(1:c_ele_num_CalP(c_E),1:3) =tmp_ele_CalP(Eles_Has_CalP_Location(c_E),1:c_ele_num_CalP(c_E),1:3) 
       call Matrix_Unique_Row_Dou(c_ele_num_CalP(c_E),3,c_ele_num_CalP(c_E),Matrix(1:c_ele_num_CalP(c_E),1:3),  &
            Uniqued_Matrix(1:c_ele_num_CalP(c_E),1:3),Uniqued_m,Uni_Mat_Count(1:c_ele_num_CalP(c_E)))   
       c_ele_num_CalP(c_E) = Uniqued_m
       tmp_ele_CalP(Eles_Has_CalP_Location(c_E),1:Uniqued_m,1:3) = Uniqued_Matrix(1:Uniqued_m,1:3)
     enddo
     !$OMP END PARALLEL DO
     
     
     !...........................................................................
     ! STEP 2: If a element contains only one or two, three, or four calculation
     !         points, further inspection is required.
     ! Whether this element is crossed by two or one crack fronts; if so,
     ! assign it to the current element
     ! Added calculation point. Added on 2022-04-22. BUGFIX2022042201.
     !...........................................................................
     do i_E = 1,Num_Elem 
       if(c_ele_num_CalP(i_E)>=1 .and. c_ele_num_CalP(i_E)<=5)then      
       !if(c_ele_num_CalP(i_E)<=5)then    
           ! Determine whether the 3D element is intersected by the outer boundary of a 3D discrete fracture,
           ! and calculate the coordinates of the intersection points
           ! Note: If it is crossed by two boundaries
           c_Yes_Inter =.False.
           call Tool_Yes_Ele_Intersected_by_3D_Crack_Outline(i_E,i_C,c_Yes_Inter,c_Inter_Point_A,c_Inter_Point_B,  &
                            c_Vertex_1,c_Vertex_2,num_Inter_Outline,c_Four_Points)

           ! Passed through by the front edges of two cracks.
           if (num_Inter_Outline==2)then
               do i_P =1,4
                 c_ele_num_CalP(i_E) = c_ele_num_CalP(i_E) +1
                 ! Expand the array to prevent out-of-bounds errors.
                 if(c_ele_num_CalP(i_E) > size(tmp_ele_CalP,2)) then
                      print *,'    Extending dim 2 of allocatable tmp_ele_CalP(:,1:extended,:) by 16...'
                      allocate(tem_tmp_ele_CalP(size(tmp_ele_CalP,1),size(tmp_ele_CalP,2)+16,3))
                      tem_tmp_ele_CalP(:,1:size(tmp_ele_CalP,2),:)=tmp_ele_CalP(:,1:size(tmp_ele_CalP,2),:)
                      deallocate(tmp_ele_CalP)
                      call move_alloc(tem_tmp_ele_CalP,tmp_ele_CalP)
                 endif 
                 ! tmp_ele_CalP(i_E, c_ele_num_CalP(i_E), 1:3) = c_Four_Points(i_P, 1:3) !Temporary variable:
                 ! coordinates of the calculation points contained in the current element
                 
                 
                 tmp_ele_CalP(Eles_Has_CalP_Location(i_E),c_ele_num_CalP(i_E),1:3)=  c_Four_Points(i_P,1:3)
                 
                 ! Check if the intersection already exists
                 call Vector_belongs_Matrix_Is_Dou(c_num_CalP,3,Cracks_CalP_Coors_3D(i_C)%row(1:c_num_CalP,1:3),&
                                                   c_Four_Points(i_P,1:3),c_Location,Yes_Exist)   
                 if(Yes_Exist .eqv. .False.)then
                   c_num_CalP =c_num_CalP +1   
                   !BUGFIX2023011001.
                   c_Max_N_CalP_3D = size(Cracks_CalP_Coors_3D(i_C)%row,1)
                   if(c_num_CalP > c_Max_N_CalP_3D)then
                        ! Expand memory. NEWFTU2022110501.
                        call D3_Allocate_Crack_Memory(i_C,3,1)
                   endif 
                   
                   Cracks_CalP_Coors_3D(i_C)%row(c_num_CalP,1:3)=c_Four_Points(i_P,1:3)      
                   !2022-09-27.
                   if(sum(abs(c_Four_Points(i_P,1:3)))<=Tol_10)then
                      print *,'    ERROR-2022092701 :: illegal Cracks_CalP_Coors_3D, in D3_Cal_HF_Crack_Points_Info_Linear.f'
                      call Warning_Message('S',Keywords_Blank)  
                   endif 
                   ! Calculate the outward normal vector of the element
                   call Tool_Normal_Unit_vector_of_3D_Tri(c_Four_Points(1,:),c_Four_Points(2,:),c_Four_Points(3,:),Unit_Norm_Vector)             
                   Cracks_CalP_Orient_3D(i_C)%row(c_num_CalP,1:3) =Unit_Norm_Vector  
                 endif
              enddo
           !oooooooooooooooooooooooooooooo
           !2023-09-20. IMPROV2023092002.
           !oooooooooooooooooooooooooooooo
           elseif(num_Inter_Outline==1)then
               do i_P =1,2
                 c_ele_num_CalP(i_E) = c_ele_num_CalP(i_E) +1
                 ! Expand the array to prevent out-of-bounds errors.
                 if(c_ele_num_CalP(i_E) > size(tmp_ele_CalP,2)) then
                      print *,'    Extending dim 2 of allocatable tmp_ele_CalP(:,1:extended,:) by 16...'
                      allocate(tem_tmp_ele_CalP(size(tmp_ele_CalP,1),size(tmp_ele_CalP,2)+16,3))
                      tem_tmp_ele_CalP(:,1:size(tmp_ele_CalP,2),:)=tmp_ele_CalP(:,1:size(tmp_ele_CalP,2),:)
                      deallocate(tmp_ele_CalP)
                      call move_alloc(tem_tmp_ele_CalP,tmp_ele_CalP)
                 endif 
                 ! tmp_ele_CalP(i_E,c_ele_num_CalP(i_E),1:3)= c_Four_Points(i_P,1:3) !Temporary variable: coordinates
                 ! of the calculation points contained in the current element
                 
                 
                 tmp_ele_CalP(Eles_Has_CalP_Location(i_E),c_ele_num_CalP(i_E),1:3)=  c_Four_Points(i_P,1:3)
                 
                 ! Check if the intersection already exists
                 call Vector_belongs_Matrix_Is_Dou(c_num_CalP,3,Cracks_CalP_Coors_3D(i_C)%row(1:c_num_CalP,1:3),&
                                                   c_Four_Points(i_P,1:3),c_Location,Yes_Exist)   
                 if(Yes_Exist .eqv. .False.)then
                   c_num_CalP =c_num_CalP +1   
                   !BUGFIX2023011001.
                   c_Max_N_CalP_3D = size(Cracks_CalP_Coors_3D(i_C)%row,1)
                   if(c_num_CalP > c_Max_N_CalP_3D)then
                        call D3_Allocate_Crack_Memory(i_C,3,1)
                   endif 
                   
                   Cracks_CalP_Coors_3D(i_C)%row(c_num_CalP,1:3)=c_Four_Points(i_P,1:3)      
                   !2022-09-27.
                   if(sum(abs(c_Four_Points(i_P,1:3)))<=Tol_10)then
                      print *,'    ERROR-2022092701 :: illegal Cracks_CalP_Coors_3D, in D3_Cal_HF_Crack_Points_Info_Linear.f'
                      call Warning_Message('S',Keywords_Blank)  
                   endif 
                   ! Calculate the outward normal vector of the element
                   call Tool_Normal_Unit_vector_of_3D_Tri(c_Four_Points(1,:),c_Four_Points(2,:),c_Four_Points(3,:),Unit_Norm_Vector)             
                   Cracks_CalP_Orient_3D(i_C)%row(c_num_CalP,1:3) =Unit_Norm_Vector  
                 endif
              enddo
           !oooooooooooo
           !2023-09-20. 
           !oooooooooooo
           elseif(num_Inter_Outline==3)then
                print *,'    WARN-2023092001 :: unforseen case in D3_Cal_HF_Crack_Points_Info_Linear.f!'
                call Warning_Message('S',Keywords_Blank)  
                !TO BE DONE.
           !oooooooooooo
           !2023-09-20. 
           !oooooooooooo
           elseif(num_Inter_Outline==4)then
                print *,'    WARN-2023092002 :: unforseen case in D3_Cal_HF_Crack_Points_Info_Linear.f!'
                call Warning_Message('S',Keywords_Blank)  
                !TO BE DONE.
           !oooooooooooo
           !2023-09-20. 
           !oooooooooooo
           elseif(num_Inter_Outline>=5)then
                print *,'    WARN-2023092003 :: unforseen case in D3_Cal_HF_Crack_Points_Info_Linear.f!'
                call Warning_Message('S',Keywords_Blank)
                !TO BE DONE.
           endif
       endif
     enddo
     Cracks_CalP_Num_3D(i_C)   = c_num_CalP

     
     !............................................................................................
     ! STEP 4: Handle the calculation point numbering and other aspects of the fluid elements in
     ! each element
     ! Store in the following variables:
     ! Cracks_FluidEle_CalP_3D(1:Max_Num_Cr_3D,1:Max_N_CalP_3D,1:9) !Calculation point numbers of
     ! fluid elements for each crack
     ! Cracks_FluidEle_num_CalP_3D(1:Max_Num_Cr_3D,1:Max_N_CalP_3D) !Number of calculation points
     ! for each fluid element in each crack
     !............................................................................................
     ! OpenMP experiment, with data dependency. 2022-10-05.
     c_count_FluidEle = 0
     do i_E = 1,Num_Elem 
       if (c_ele_num_CalP(i_E) >=3)then
          !IMPROV2023081303. 
          ! Check whether the number of fluid nodes (calculation points) exceeds 20.
          if(c_ele_num_CalP(i_E) >Max_ele_num_CalP)then
              print *,'    Error-2023051601 :: c_ele_num_CalP(i_E) > Max_ele_num_CalP (',Max_ele_num_CalP,')!'
              print *,'            In D3_Cal_HF_Crack_Points_Info_Linear.f'
              print *,'            c_ele_num_CalP(i_E):',c_ele_num_CalP(i_E)
              print *,'            i_E:',i_E
              call Warning_Message('S',Keywords_Blank)  
              !return
          endif   
          
          c_count_FluidEle = c_count_FluidEle +1

          
          ! Check if the array is out of bounds. 2022-10-04.
          c_Max_N_FluEl_3D = size(Cracks_FluidEle_CalP_3D(i_C)%row,1)
          if(c_count_FluidEle > c_Max_N_FluEl_3D)then
              ! Expand memory. NEWFTU2022110502.
              call D3_Allocate_Crack_Memory(i_C,2,1)
              !c_Max_N_FluEl_3D = size(Cracks_FluidEle_CalP_3D(i_C)%row,1)
          endif
            
          ! Number of calculation points for each fracture fluid element
          Cracks_FluidEle_num_CalP_3D(i_C)%row(c_count_FluidEle)=  c_ele_num_CalP(i_E)
          ! Obtain a calculation point number
          do i_CalP=1,c_ele_num_CalP(i_E)                              
            call Vector_belongs_Matrix_Is_Dou(Cracks_CalP_Num_3D(i_C),3,&
                                     Cracks_CalP_Coors_3D(i_C)%row(1:Cracks_CalP_Num_3D(i_C),1:3), &
                                     tmp_ele_CalP(Eles_Has_CalP_Location(i_E),i_CalP,1:3),c_Location,c_Yes)           
            Cracks_FluidEle_CalP_3D(i_C)%row(c_count_FluidEle,i_CalP) = c_Location
            Cracks_FluidEle_EleNum_3D(i_C)%row(c_count_FluidEle) = i_E
          enddo
       endif
     enddo
     Cracks_FluidEle_num_3D(i_C) = c_count_FluidEle

         
     !.........................................................................................
     ! STEP 5: Reorder the nodes of each fluid element so that the polygons formed by the node 
     ! connections are standard polygons, rather than disordered ones.
     ! Sorting steps:
     ! (1) Calculate the centroid of all coordinates;
     ! (2) Calculate the angle between each point and the midpoint;
     ! (3) Sort by angle size.
     !.........................................................................................
     ! OpenMP parallel implementation. No data dependency. IMPROV2022100501.
     !$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_FluidEle,num_Fluid_nodes,Angle_Lines,&
     !$OMP          Center_X,Center_Y,Center_Z,Point_Center,Fluid_node1,Line_O_1,Line_O_N,&
     !$OMP          i_FluidNode,Fluid_nodeN,vector_a,vector_b,Old_Cracks_FluidEle_CalP_3D,&
     !$OMP          i_temp,c_location) 
     do i_FluidEle = 1,Cracks_FluidEle_num_3D(i_C)
       num_Fluid_nodes=Cracks_FluidEle_num_CalP_3D(i_C)%row(i_FluidEle) 
       if (num_Fluid_nodes<=3) then
           cycle
       endif
       allocate(Angle_Lines(num_Fluid_nodes-1))
       !Center of all fluid ndoes of current fluid element.
       Center_X = sum(Cracks_CalP_Coors_3D(i_C)%row(&
                   Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,1:num_Fluid_nodes),1))/num_Fluid_nodes
       Center_Y = sum(Cracks_CalP_Coors_3D(i_C)%row(&
                   Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,1:num_Fluid_nodes),2))/num_Fluid_nodes
       Center_Z = sum(Cracks_CalP_Coors_3D(i_C)%row(&
                   Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,1:num_Fluid_nodes),3))/num_Fluid_nodes     
       Point_Center(1:3) = [Center_X,Center_Y,Center_Z]
       ! Fluid Node 1
       Fluid_node1 = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,1)
       
       ! Denote the line connecting the center of mass to fluid node 1 as Line_O_1
       Line_O_1(1,1:3) = Point_Center(1:3)
       Line_O_1(2,1:3) = Cracks_CalP_Coors_3D(i_C)%row(Fluid_node1,1:3)
        
       ! Calculate the angle between Line_O_1 and Line_O_N
       Line_O_N(1,1:3) = Point_Center(1:3)
       
       do i_FluidNode = 2,num_Fluid_nodes
           Fluid_nodeN =Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,i_FluidNode)
           Line_O_N(2,1:3)=Cracks_CalP_Coors_3D(i_C)%row(Fluid_nodeN,1:3)
           vector_a = Line_O_1(2,1:3) - Line_O_1(1,1:3)
           vector_b = Line_O_N(2,1:3) - Line_O_N(1,1:3)
           ! Calculate the angle and store it in Angle_Lines.
           call Tool_Angle_of_Vectors_a_and_b_3D(vector_a,vector_b,Angle_Lines(i_FluidNode-1),1)
       enddo
       
       !Old_Cracks_FluidEle_CalP_3D(1:9) = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,1:9)
       
       !2023-08-13.
       Old_Cracks_FluidEle_CalP_3D(1:Max_ele_num_CalP) = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,1:Max_ele_num_CalP)
       
       do i_temp = 2,num_Fluid_nodes
           c_location = minloc(Angle_Lines(1:num_Fluid_nodes-1),1)
           Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,i_temp)=Old_Cracks_FluidEle_CalP_3D(c_location+1)
           Angle_Lines(c_location) = Con_Big_20
       enddo
       
       deallocate(Angle_Lines)  
     enddo
     !$OMP END PARALLEL DO 
     
     
     !...............................................................................................
     ! STEP 6: Split the 3D fluid elements, ensuring that each fluid element is a triangular element
     !         (executed when Key_3D_FluEle_Triang=1).
     !...............................................................................................
     if (Key_3D_FluEle_Triang==1)then
       c_count_FluEl = Cracks_FluidEle_num_3D(i_C)
       do i_FluidEle = 1,Cracks_FluidEle_num_3D(i_C)
         num_Fluid_nodes=Cracks_FluidEle_num_CalP_3D(i_C)%row(i_FluidEle) 
         c_Solid_El = Cracks_FluidEle_EleNum_3D(i_C)%row(i_FluidEle)
         
         if (num_Fluid_nodes<=3) then
             cycle
         endif
         if (num_Fluid_nodes>=4) then 
           fnode_1 = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,1)
           fnode_2 = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,2)
           fnode_3 = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,3)
           Cracks_FluidEle_num_CalP_3D(i_C)%row(i_FluidEle) = 3            
           do i_added=1,num_Fluid_nodes-3
             c_count_FluEl = c_count_FluEl + 1
             
             ! Check if the array is out of bounds. 2022-10-04.
             c_Max_N_FluEl_3D = size(Cracks_FluidEle_CalP_3D(i_C)%row,1)
             if(c_count_FluEl>c_Max_N_FluEl_3D)then
                  ! Expand memory. NEWFTU2022110502.
                  call D3_Allocate_Crack_Memory(i_C,2,1)
                  !c_Max_N_FluEl_3D = size(Cracks_FluidEle_CalP_3D(i_C)%row,1)
             endif
          
             c_node2 = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,i_added+2)   
             c_node3 = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,i_added+3) 
             ! Add Fluid Element
             Cracks_FluidEle_CalP_3D(i_C)%row(c_count_FluEl,1)=fnode_1 
             Cracks_FluidEle_CalP_3D(i_C)%row(c_count_FluEl,2)=c_node2
             Cracks_FluidEle_CalP_3D(i_C)%row(c_count_FluEl,3)=c_node3
             
             Cracks_FluidEle_num_CalP_3D(i_C)%row(c_count_FluEl) = 3
             Cracks_FluidEle_EleNum_3D(i_C)%row(c_count_FluEl)=c_Solid_El
           enddo  
         endif 
       enddo
       ! Update the total number of fluid elements
       Cracks_FluidEle_num_3D(i_C) =  c_count_FluEl
     endif

     
1101 FORMAT(5X,'Number of fluid elements of crack ',I4,':',I8,' / ',I8)   
     c_Max_N_FluEl_3D = size(Cracks_FluidEle_CalP_3D(i_C)%row,1)    
     write(*,1101) i_C,Cracks_FluidEle_num_3D(i_C),c_Max_N_FluEl_3D
     
     ! The following algorithm has been replaced by a better one. For details, see BUGFIX2021081601,
     ! 2021-08-16.
     !.................................................................................
     ! STEP 6: Determination of the fluid elements and nodes corresponding to crack 
     ! tip enhancement Key_Tip_Fluid_Element_3D==1
     ! See algorithm: See algorithm\Goodnotes\PhiPsi algorithm and others\P2
     ! Date Added: 2020-01-21
     ! Note: This algorithm is not suitable for edge cracks and is only applicable 
     ! to cracks that are entirely within the model. For faults involving edge cracks,
     ! see the diagram.
     ! \Goodnotes\PhiPsi Algorithm and Others\P3, 2021-04-24
     !.................................................................................
     ! The code has been deleted
 
     !......................................................................
     ! STEP 7: Check whether every fluid node is being utilized, 2020-01-27
     !......................................................................
     ! OpenMP parallel implementation. IMPROV2022100502.
     !$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_Fnode,c_Yes_Exist,i_FluidEle,c_E_ndoes) 
     do i_Fnode = 1,Cracks_CalP_Num_3D(i_C)
         c_Yes_Exist = .False.
         do i_FluidEle = 1,Cracks_FluidEle_num_3D(i_C)
             c_E_ndoes(1:3) = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,1:3)
             if (any(c_E_ndoes(1:3)== i_Fnode))then
                 c_Yes_Exist = .True.
                 exit
             endif
         enddo  
         if(c_Yes_Exist .eqv. .False.)then
           print *,'    WARNING :: unemployed fluid nodes ',i_Fnode,' found!'        
           print *,'               in D3_Cal_HF_Crack_Points_Info_Linear.f, code-2020!'  
         endif
     enddo
     !$OMP END PARALLEL DO
         
     !...............................................................................................
     ! STEP 8: Calculate the area, centroid, outward normal vector, and fluid node numbering for
     ! each fluid element
     ! Store in the following variables:
     ! Cracks_FluidEle_Area_3D(1:Max_Num_Cr_3D,1:Max_N_CalP_3D) !Area of each crack's fluid element
     ! Cracks_FluidEle_Centroid_3D(1:Max_Num_Cr_3D,1:Max_N_CalP_3D) !Centroid of the fluid element
     ! for each crack
     ! Cracks_FluidEle_Vector_3D(1:Max_Num_Cr_3D,1:Max_N_CalP_3D,1:3) !Average outward normal vector
     ! of the fluid element for each crack
     !...............................................................................................
     !//////////////////////////////////////////
     ! Calculate the area of each fluid element
     !//////////////////////////////////////////
     !$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_FluidEle,P1_CalP_num,c_P1,c_Fluid_Area,i_sub_tri,P2_CalP_num,P3_CalP_num,&
     !$OMP                 c_P2, c_P3, delta_Fluid_Area) &    !BUGFIX2023071701. Missed delta_Fluid_Area. 2023-07-17.
     !$OMP            SCHEDULE(static)      
     ! OpenMP parallelization. IMPROV2022092302.
     do i_FluidEle = 1,Cracks_FluidEle_num_3D(i_C)
        ! Sub-triangle loop of fluid elements
        P1_CalP_num = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,1)     
        c_P1(1:3)   = Cracks_CalP_Coors_3D(i_C)%row(P1_CalP_num,1:3)
        c_Fluid_Area = ZR
        do i_sub_tri=1,Cracks_FluidEle_num_CalP_3D(i_C)%row(i_FluidEle)-2
            P2_CalP_num = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,i_sub_tri+1)   
            P3_CalP_num = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,i_sub_tri+2) 
            c_P2(1:3) =Cracks_CalP_Coors_3D(i_C)%row(P2_CalP_num,1:3)
            c_P3(1:3) =Cracks_CalP_Coors_3D(i_C)%row(P3_CalP_num,1:3)
            ! Calculate the area of a subspace triangle
            call Tool_Area_Tri_3D(c_P1,c_P2,c_P3,delta_Fluid_Area)
            c_Fluid_Area = c_Fluid_Area + delta_Fluid_Area
            
        enddo
        Cracks_FluidEle_Area_3D(i_C)%row(i_FluidEle) = c_Fluid_Area 
     enddo    
     !$OMP END PARALLEL DO
         
     !//////////////////////////////////////////////////////////////////////////////////////////////
     ! Calculate the outward normal vector of each fluid element (the average outward normal vector
     ! at fluid computation points)
     ! and the centroid,
     ! Local coordinate system x-axis, y-axis, z-axis (normal) of the centroid position of each
     ! fluid element in the crack
     !//////////////////////////////////////////////////////////////////////////////////////////////
     !$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_FluidEle,c_Fluid_Vector,c_Fluid_Centroid,i_CalP,c_CalP,c_Norm_2,c_Solid_El,&
     !$OMP                 c_z_vector,c_First_CalP,Coor_First_CalP,c_x_vector,c_y_vector,ThetaX,ThetaY,ThetaZ,c_T_Matrx)    &
     !$OMP            SCHEDULE(static)      
     do i_FluidEle = 1,Cracks_FluidEle_num_3D(i_C)
        !Ele_Num_Cache = 1 !2022-09-24. IMPROV2022092401. 
        ! Calculation Point Loop of Fluid Element
        c_Fluid_Vector(1:3) = ZR
        c_Fluid_Centroid(1:3) = ZR
        do i_CalP =1,Cracks_FluidEle_num_CalP_3D(i_C)%row(i_FluidEle)
            c_CalP = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,i_CalP)  
            c_Fluid_Vector = c_Fluid_Vector +Cracks_CalP_Orient_3D(i_C)%row(c_CalP,1:3)
            c_Fluid_Centroid = c_Fluid_Centroid +Cracks_CalP_Coors_3D(i_C)%row(c_CalP,1:3)     
        enddo
        c_Fluid_Vector = c_Fluid_Vector/Cracks_FluidEle_num_CalP_3D(i_C)%row(i_FluidEle)
        ! Regularization (Normalization)
        call Vector_Norm2(3,c_Fluid_Vector,c_Norm_2)   
        c_Fluid_Vector = c_Fluid_Vector/c_Norm_2
        Cracks_FluidEle_Vector_3D(i_C)%row(i_FluidEle,1:3)=c_Fluid_Vector
        ! centroid
        c_Fluid_Centroid = c_Fluid_Centroid/Cracks_FluidEle_num_CalP_3D(i_C)%row(i_FluidEle)
        Cracks_FluidEle_Centroid_3D(i_C)%row(i_FluidEle,1:3)=c_Fluid_Centroid(1:3)   
        ! If the solid element number corresponding to the fluid element has not yet been obtained, it is
        ! determined based on the element number at the centroid of the fluid element (2020-01-22; used for
        ! fluid elements corresponding to 3D crack tip enrichment elements)
        if (Cracks_FluidEle_EleNum_3D(i_C)%row(i_FluidEle)==0) then
            !call
            !Cal_Ele_Num_by_Coors_3D(c_Fluid_Centroid(1),c_Fluid_Centroid(2),c_Fluid_Centroid(3),Ele_Num_Cache,c_Solid_El)
            ! No cache.
            call Cal_Ele_Num_by_Coors_3D_v2(c_Fluid_Centroid(1),c_Fluid_Centroid(2),c_Fluid_Centroid(3),c_Solid_El)
            Cracks_FluidEle_EleNum_3D(i_C)%row(i_FluidEle)=c_Solid_El
        endif

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Local coordinate system x-y-z axes of the centroid position of each fracture fluid element
        ! (z-axis is normal)
        ! Date Added: 2020-01-19
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        c_z_vector = Cracks_FluidEle_Vector_3D(i_C)%row(i_FluidEle,1:3)
        ! The vector from the centroid to the first fluid node is defined as the x-axis
        c_First_CalP = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,1) 
        Coor_First_CalP = Cracks_CalP_Coors_3D(i_C)%row(c_First_CalP,1:3)
        c_x_vector =   Coor_First_CalP -c_Fluid_Centroid
        call Vector_Norm2(3,c_x_vector,c_Norm_2)   
        c_x_vector = c_x_vector/c_Norm_2
        ! The y-axis is obtained by the cross product of the z-axis and the x-axis (note that it is z × x,
        ! not x × z).
        call Vector_Cross_Product_3(c_z_vector,c_x_vector,c_y_vector) 
        call Vector_Norm2(3,c_y_vector,c_Norm_2)   
        c_y_vector = c_y_vector/c_Norm_2
        Cracks_FluidEle_LCS_x_3D(i_C)%row(i_FluidEle,1:3) =c_x_vector            
        Cracks_FluidEle_LCS_y_3D(i_C)%row(i_FluidEle,1:3) =c_y_vector     
        Cracks_FluidEle_LCS_z_3D(i_C)%row(i_FluidEle,1:3) =c_z_vector    
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Local coordinate system x-y-z axes (with the z-axis as the normal) of the centroid
        ! position of each fracture fluid element
        ! Transformation Matrix
        ! Date Added: 2020-01-19
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Determine ThetaX, ThetaY, ThetaZ so that all three sets of transformations above are satisfied.
        call Tool_ThetaX_ThetaY_ThetaZ_3D_rotation(1, c_x_vector,c_y_vector,c_z_vector,    &
                     [ONE,ZR,ZR],[ZR,ONE,ZR],[ZR,ZR,ONE],ThetaX,ThetaY,ThetaZ,c_T_Matrx(1:3,1:3))
        Cracks_FluidEle_LCS_T_3D(i_C)%row(i_FluidEle,1:3,1:3)=c_T_Matrx
     enddo  
     !$OMP END PARALLEL DO
  
     !///////////////////////////////////////////////////////////////////////////////////////////
     ! Calculate the actual number of fluid nodes participating in the current crack calculation
     !(store in: Cracks_Real_CalP_Num_3D(1:Max_Num_Cr_3D))
     ! Only HF fractures are counted here. IMPROV2022060201.
     !///////////////////////////////////////////////////////////////////////////////////////////
     1001 FORMAT(5X,'Number of fluid nodes of crack    ',I4,':',I8,' / ',I8)
     if (Crack_Type_Status_3D(i_C,1) == 1) then
         !tem_real_CalP(1:Max_N_CalP_3D) = 0
         c_Max_N_CalP_3D = size(Cracks_CalP_Coors_3D(i_C)%row,1)
         if(allocated(tem_real_CalP)) deallocate(tem_real_CalP)
         ALLOCATE(tem_real_CalP(10*c_Max_N_CalP_3D))   
         if(allocated(Uniqued_Vec_tem_real_CalP)) deallocate(Uniqued_Vec_tem_real_CalP)
         ALLOCATE(Uniqued_Vec_tem_real_CalP(10*c_Max_N_CalP_3D))  
         tem_real_CalP(1:10*c_Max_N_CalP_3D) = 0
         tem_real_num = 0
         do i_FluidEle = 1,Cracks_FluidEle_num_3D(i_C)
            ! Calculation Point Loop of Fluid Element
            do i_CalP =1,Cracks_FluidEle_num_CalP_3D(i_C)%row(i_FluidEle)
                c_CalP = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,i_CalP)  
                tem_real_num = tem_real_num +1
                tem_real_CalP(tem_real_num)  = c_CalP
            enddo
         enddo  
         ! Check array out-of-bounds
         if (tem_real_num > 10*c_Max_N_CalP_3D) then
              print *,'    ERROR-2023081301 :: tem_real_num > 10*c_Max_N_CalP_3D, in D3_Cal_HF_Crack_Points_Info_Linear.f!'
              print *,'                        Try to increase to 20* in D3_Cal_HF_Crack_Points_Info_Linear.f!'
              call Warning_Message('S',Keywords_Blank) 
         endif
         
         ! Delete duplicate fluid nodes
         call Vector_Unique_Int(10*c_Max_N_CalP_3D,tem_real_num,tem_real_CalP,Uniqued_Vec_tem_real_CalP,tem_Uniqued_n)  
         !Check if tem_Uniqued_n > Max_N_CalP_3D (2021-08-21)
         write(*,1001) i_C,tem_Uniqued_n,c_Max_N_CalP_3D
         
         if(tem_Uniqued_n > c_Max_N_CalP_3D)then
            ! Expand memory. NEWFTU2022110501.
            call D3_Allocate_Crack_Memory(i_C,3,1)
         endif  
         
         Cracks_Real_CalP_Num_3D(i_C) = tem_Uniqued_n
         !Cracks_Real_CalPs_3D(i_C,1:tem_Uniqued_n) = Uniqued_Vec_tem_real_CalP(1:tem_Uniqued_n)
     endif
         
     !/////////////////////////////////////////////////////////////////////////////////////////
     ! Global numbering of fluid nodes, stored in Cracks_FluidEle_Glo_CalP_3D(i_C, i_FluidEle,
     ! i_CalP)
     ! Note: For hydraulic fracturing fractures only
     !/////////////////////////////////////////////////////////////////////////////////////////
     ! If the current crack is an HF crack. 2022-06-02. IMPROV2022060201.
     if (Crack_Type_Status_3D(i_C,1) == 1) then
       ! OpenMP parallelization. IMPROV2022092302.
       do i_FluidEle = 1,Cracks_FluidEle_num_3D(i_C)
        ! Solid element number where the fluid element is located
        Solid_Elem = Cracks_FluidEle_EleNum_3D(i_C)%row(i_FluidEle)
        ! Average outward normal vector of the fluid element
        ori_n =  Cracks_FluidEle_Vector_3D(i_C)%row(i_FluidEle,1:3) 
        ! Calculation Point Loop of Fluid Element
        do i_CalP =1,Cracks_FluidEle_num_CalP_3D(i_C)%row(i_FluidEle)
            c_CalP = Cracks_FluidEle_CalP_3D(i_C)%row(i_FluidEle,i_CalP)
            call Vector_Location_Int(tem_Uniqued_n,Uniqued_Vec_tem_real_CalP(1:tem_Uniqued_n),c_CalP,cc_location,c_Yes_In) 
            if (i_C==1)then
                Cracks_FluidEle_Glo_CalP_3D(i_C)%row(i_FluidEle,i_CalP) =cc_location
                Glo_location = cc_location
            elseif(i_C >= 2)then
                !Cracks_FluidEle_Glo_CalP_3D(i_C,i_FluidEle,i_CalP)=cc_location + Cracks_Real_CalP_Num_3D(i_C-1)
                Cracks_FluidEle_Glo_CalP_3D(i_C)%row(i_FluidEle,i_CalP)= cc_location +sum(Cracks_Real_CalP_Num_3D(1:i_C-1))
                Glo_location = cc_location +sum(Cracks_Real_CalP_Num_3D(1:i_C-1))  
            endif
            ! Save the local information corresponding to the global fluid node numbers: including fracture
            ! numbers, fluid element numbers, and the fluid node numbers of the corresponding fluid elements
            Cracks_FluidEle_CalP_Glo_Info(Glo_location,1) = i_C
            Cracks_FluidEle_CalP_Glo_Info(Glo_location,2) = i_FluidEle
            Cracks_FluidEle_CalP_Glo_Info(Glo_location,3) = i_CalP
          
            ! Obtain the normal initial ground stress at the fluid node positions of hydraulic fracturing
            ! fracture surfaces under global numbering. 2022-06-04.
            !NEWFTU2022060401.
            if(Key_InSitu_Strategy == 4 )then
              ! Obtained in-situ stress tensor c_Insitu_Stress
              c_Insitu_Stress(1) =InSitu_Strs_Gaus_xx(Solid_Elem,1)
              c_Insitu_Stress(2) =InSitu_Strs_Gaus_yy(Solid_Elem,1)                    
              c_Insitu_Stress(3) =InSitu_Strs_Gaus_zz(Solid_Elem,1)                    
              c_Insitu_Stress(4) =InSitu_Strs_Gaus_xy(Solid_Elem,1)                    
              c_Insitu_Stress(5) =InSitu_Strs_Gaus_yz(Solid_Elem,1)                    
              c_Insitu_Stress(6) =InSitu_Strs_Gaus_xz(Solid_Elem,1)                    
              ! Calculate the normal stress on the fracture surface
              call Tool_Get_Normal_Stress_on_Plane_by_Stress_Tensor(c_Insitu_Stress,ori_n,Normal_Stress)
              Cracks_FluidEle_CalP_Glo_Insitu(Glo_location) =Normal_Stress 
            endif
        enddo
       enddo
     endif
     
     ! Clear temporary memory variables.
     if(allocated(tem_real_CalP)) deallocate(tem_real_CalP)
     if(allocated(Uniqued_Vec_tem_real_CalP)) deallocate(Uniqued_Vec_tem_real_CalP)
     if(allocated(tmp_ele_CalP)) deallocate(tmp_ele_CalP)
end do

!=====================================================================
! Updated the maximum values of several global variables. 2023-08-13.
!=====================================================================
Max_Max_N_FluEl_3D = maxval(Max_N_FluEl_3D)  
Max_Max_N_CalP_3D  = maxval(Max_N_CalP_3D) 
Max_Max_N_Node_3D  = maxval(Max_N_Node_3D) 


!**************************************************
! Call Create_Tri_Mesh_From_Points.exe to generate 
! output_triangulation.txt file.
!**************************************************
!call system ('Create_Tri_Mesh_From_Points.exe')   
#endif

#ifdef Silverfrost
print *,'    ERROR :: Silverfrost compiler failed to compile D3_Cal_HF_Crack_Points_Info_Linear.F90!'
print *,'             In D3_Cal_HF_Crack_Points_Info_Linear.F90.'
call Warning_Message('S',Keywords_Blank)
#endif      
return 
end SUBROUTINE D3_Cal_HF_Crack_Points_Info_Linear
