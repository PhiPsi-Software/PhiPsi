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
 
SUBROUTINE D3_Prepare_Subdivision_Integration(c_isub)
! Block-wise integral related data calculation. Abandoned (2022-07-28).
! Ref: My PhiPsi Development Notebook V1-P54-55. NEWFTU2022072701.
! 2022-07-27.

! Variable Description:
! integer Elems_Integration_Type(num_Elem, num_Crack) ! Integration type, =0, standard fixed
! integration.
! =1, Case 1, there are 3 intersections
! =2, Case 2, there are 4 intersections, which can be further divided into 2 situations
! =3, Case 3, there are 5 intersections
! =4, Case 4, there are 6 intersections
! Note: Depending on the different cracks, the integration scheme for the same element may vary.
!
! integer, ALLOCATABLE :: Elems_Integration_Type(num_Elem, num_Crack) ! Integration type.
! 2022-07-27.
! integer, ALLOCATABLE :: Elems_Num_SubEles(num_Elem, num_Crack)           ! Number of sub-elements
! integer, ALLOCATABLE :: Elems_Type_SubEles(num_Elem, num_Crack) ! Sub-element type, =1 indicates a
! hexahedral element, =2 indicates a tetrahedral element
! integer, ALLOCATABLE :: Elems_SubEles_Index(num_Elem, num_Crack) ! The index of the sub-block
! elements in the SubEles_Coors matrix.
! integer num_SubEles                                                  ! Total number of sub-block elements.
! integer, ALLOCATABLE :: SubEles_Integ_Num(num_Elem) ! Number of integration points in
! sub-elements.
! Specifically used to store the coordinates of integration points for sub-elements. 
!    1: sub-element number; 2: integration point number; 3: local coordinates and weights of the
!    integration point.

!-----------------------------
! Read public variable module
!-----------------------------
use Global_Float_Type
use Global_Common
use Global_Model
use Global_Elem_Area_Vol
use Global_Crack_Common
use Global_Crack_3D
      
      
!---------------------------
! Variable Type Declaration
!---------------------------
implicit none
integer,intent(in):: c_isub
integer i_C,i_Crack_Ele,Crack_Node1,Crack_Node2,Crack_Node3
real(kind=FT) Point1(3),Point2(3),Point3(3)
real(kind=FT) max_x_Tri,max_y_Tri,max_z_Tri
real(kind=FT) min_x_Tri,min_y_Tri,min_z_Tri
integer i_E,i_E_Edge
logical Logical_Yes_x,Logical_Yes_y,Logical_Yes_z
real(kind=FT) A(3),B(3)
logical Yes_Inter
real(kind=FT) InterSection_P(3),Unit_Norm_Vector(3)
integer c_num_InterP,c_Location,c_ele_num_InterP(Num_Elem)
logical Yes_Exist
integer Edge_node1,Edge_node2
real(kind=FT) tmp_ele_InterP(Num_Elem,50,3)
integer tmp_ele_InterP_Edge(Num_Elem,50)
integer Index_Vetor_4(4)
integer tem_Vector_4(4)
print *,'    Preparing data for subdivision integration...'

!----------------
! Initialization
!----------------
IF(ALLOCATED(Elems_Integration_Type)) DEALLOCATE(Elems_Integration_Type)    
ALLOCATE(Elems_Integration_Type(num_Elem,num_Crack))       
Elems_Integration_Type(1:num_Elem,1:num_Crack) = 0

IF(ALLOCATED(Elems_Num_SubEles)) DEALLOCATE(Elems_Num_SubEles)    
ALLOCATE(Elems_Num_SubEles(num_Elem,num_Crack))       
Elems_Num_SubEles(1:num_Elem,1:num_Crack) = 0

IF(ALLOCATED(Elems_Type_SubEles)) DEALLOCATE(Elems_Type_SubEles)    
ALLOCATE(Elems_Type_SubEles(num_Elem,num_Crack))       
Elems_Type_SubEles(1:num_Elem,1:num_Crack) = 0

IF(ALLOCATED(Elems_SubEles_Index)) DEALLOCATE(Elems_SubEles_Index)    
ALLOCATE(Elems_SubEles_Index(num_Elem,num_Crack))       
Elems_SubEles_Index(1:num_Elem,1:num_Crack) = 0

IF(ALLOCATED(SubEles_Integ_Num)) DEALLOCATE(SubEles_Integ_Num)    
ALLOCATE(SubEles_Integ_Num(num_Elem))       
SubEles_Integ_Num(1:num_Elem) = 0

IF(ALLOCATED(SubEles_Integ_Coors)) DEALLOCATE(SubEles_Integ_Coors)
ALLOCATE(SubEles_Integ_Coors(num_Elem,200,4))      
SubEles_Integ_Coors(1:num_Elem,1:200,1:4) = ZR

num_SubEles = 0

  
! ALLOCATE(SubEles_Coors(1000000,700,3)) !Memory space usage calculation;
! 1e6*700*3*64/1024/1024/1024/8 = 15.65 GB.

!---------------     
! element cycle
!---------------   
c_num_InterP =0
c_ele_num_InterP(1:Num_Elem) =0
! Temporary variable: coordinates of calculation points contained in the current cell
tmp_ele_InterP(1:Num_Elem,1:50,1:3) = ZR       
! Temporary variable: the edge number corresponding to the intersection of cracks contained in the
! current element and the element
tmp_ele_InterP_Edge(1:Num_Elem,1:50) = 0       
do i_E = 1,Num_Elem
  ! Fracture cycle
  do i_C = 1,num_Crack 
    ! If the current element is a Heaviside-enriched element relative to the current crack, further
    ! calculations are performed.
    if(Elem_Type_3D(i_E,i_C)==2)then
        !//////////////////////////////////////////////////////////////////////////////////////////
        ! Discrete fracture surface loop. Obtain the number and coordinates of intersection points
        ! between the fracture surface and the elements.
        !//////////////////////////////////////////////////////////////////////////////////////////
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
          !       of element edges and discrete fracture surfaces.
          !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          !..............................................................
          ! Check whether the coordinate ranges intersect (overlap). 
          ! If x, y, or z do not overlap, an intersection is impossible.
          ! element coordinate range and triangular fracture surface
          ! element coordinate range. 2022-06-12. IMPROV2022061202.
          !..............................................................
          call Tool_Yes_Two_Ranges_Overlapped_Double([x_min_Elements(i_E),x_max_Elements(i_E)],&
                                                     [min_x_Tri,max_x_Tri],Logical_Yes_x) 
          if(Logical_Yes_x .eqv. .False.) cycle
          call Tool_Yes_Two_Ranges_Overlapped_Double([y_min_Elements(i_E),y_max_Elements(i_E)],&
                                                     [min_y_Tri,max_y_Tri],Logical_Yes_y) 
          if(Logical_Yes_y .eqv. .False.) cycle
          call Tool_Yes_Two_Ranges_Overlapped_Double([z_min_Elements(i_E),z_max_Elements(i_E)],&
                                                     [min_z_Tri,max_z_Tri],Logical_Yes_z) 
          if(Logical_Yes_z .eqv. .False.) cycle           
                    
          do i_E_Edge = 1,12
              ! The variable Element_Edges(Num_Elem,12,2) stores the 12 edges (corresponding node numbers) of the
              ! cell.
              !Edge_node1 = Element_Edges(i_E,i_E_Edge,1)
              !Edge_node2 = Element_Edges(i_E,i_E_Edge,2)
              ! The variable Element_Edges(12,2,Num_Elem) stores the 12 edges (corresponding node numbers) of the
              ! element.
              Edge_node1 = Element_Edges(i_E_Edge,1,i_E)
              Edge_node2 = Element_Edges(i_E_Edge,2,i_E)              
              A(1:3) = Coor(Edge_node1,1:3)
              B(1:3) = Coor(Edge_node2,1:3)
              ! Calculate the intersection points between the current edge segments and spatial triangles
              ! (discrete fracture plane triangles)
              call Tool_Intersection_of_AB_and_Triangle_3D(A,B,Point1,Point2,Point3,Yes_Inter,InterSection_P)
              ! If there is an intersection
              if(Yes_Inter .eqv. .True.) then
                !------------------------------------------
                ! Check if the intersection already exists
                !------------------------------------------
                call Vector_belongs_Matrix_Is_Dou(c_ele_num_InterP(i_E),3,tmp_ele_InterP(i_E,1:c_ele_num_InterP(i_E),1:3),&                 
                                                  InterSection_P(1:3),c_Location,Yes_Exist)   
 
                !------------------------------------------------------------------ 
                ! If this intersection does not exist, add a new calculation point
                !------------------------------------------------------------------
                if(Yes_Exist .eqv. .False.)then
                    ! The intersection count contained in the current element is accumulated.
                    c_ele_num_InterP(i_E) =c_ele_num_InterP(i_E) +1
                    ! Temporary variable: Coordinates of the calculation points contained in the current cell.
                    tmp_ele_InterP(i_E,c_ele_num_InterP(i_E),1:3)=InterSection_P  
                    ! Temporary variable: the edge number corresponding to the intersection of cracks contained in the
                    ! current element and the element.
                    tmp_ele_InterP_Edge(i_E,c_ele_num_InterP(i_E)) = i_E_Edge      
                endif
              endif
          enddo
        enddo
        
        !//////////////////////////////////////////////////////////////
        ! If there are 4 intersections, handle it according to Case 2.
        !////////////////////////////////////////////////////////////// 
        if(c_ele_num_InterP(i_E)==4)then
            !Elems_Integration_Type(i_E,i_C) = 2
            call D3_Subdivision_Integration_SubEles_Case2(i_E,i_C,c_ele_num_InterP(i_E),   &
                   tmp_ele_InterP_Edge(i_E,1:c_ele_num_InterP(i_E)),tmp_ele_InterP(i_E,1:c_ele_num_InterP(i_E),1:3))
            
        endif
        
    endif
  enddo
enddo
RETURN
END SUBROUTINE D3_Prepare_Subdivision_Integration
