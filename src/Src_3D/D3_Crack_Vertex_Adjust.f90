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
 
subroutine D3_Crack_Vertex_Adjust(iter)
! Adjust the crack edge nodes (to prevent them from being on the surface of solid elements).
! Date added: 2020-01-17


!***************
! Public Module
!***************
use Global_Crack_Common
use Global_Crack_3D
use Global_Elem_Area_Vol
use Global_Model
use Global_Common
use omp_lib

!**********************
! Variable Declaration
!**********************
implicit none
integer,intent(in)::iter
integer i_V,num_Vertex,c_Mesh_Node
real(kind=FT) c_x,c_y,c_z
real(kind=FT) c_X_NODES(8),c_Y_NODES(8),c_Z_NODES(8)
real(kind=FT) mid_Ele_x,mid_Ele_y,mid_Ele_z
integer c_NN(8),Triangle(12,3),i_Tri,El 
real(kind=FT) Line_AB(2,3),new_Line_AB(2,3),delta_L
real(kind=FT) new_Point(3)
logical c_Yes_on
integer c_NN_El,i_C,i_E
real(kind=FT) c_X_NODES_El(8),c_Y_NODES_El(8),c_Z_NODES_El(8)
real(kind=FT) Point1(3),Point2(3),Point3(3)
integer c_OUT_Elem,Flag_Exit
integer c_Ele_Domain_ID,c_Ele

!*************
! Crack cycle
!*************
! OpenMP Parallelization. 2022-06-12. IMPROV2022061201.
!$OMP PARALLEL do DEFAULT(SHARED) PRIVATE(i_C,i_V,num_Vertex,          &
!$OMP&              c_Mesh_Node,c_x,c_y,c_z,Flag_Exit,i_E,i_Tri,       &
!$OMP&              c_X_NODES,c_Y_NODES,c_Z_NODES,Triangle,            &
!$OMP&              c_X_NODES_El,c_Y_NODES_El,c_Z_NODES_El,            &
!$OMP&              Point1,Point2,Point3,c_Yes_on,El,delta_L,          &
!$OMP&              mid_Ele_x,mid_Ele_y,mid_Ele_z,Line_AB,new_Point,   &
!$OMP&              new_Line_AB,c_Ele_Domain_ID,c_Ele)                 &   
!$OMP&            SCHEDULE(static)    
do i_C =1,num_Crack         
    num_Vertex = Crack3D_Meshed_Outline_num(i_C)
    do i_V =1,num_Vertex
      ! Discrete node number at the crack surface boundary
      c_Mesh_Node = Crack3D_Meshed_Outline(i_C)%row(i_V,1)    
      ! Coordinates of discrete fracture edge nodes
      c_x  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,1) 
      c_y  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,2) 
      c_z  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,3) 
      
      Flag_Exit = 0
      
      !2022-11-24. NEWFTU2022112401.
      call D3_Get_Point_Domain_Number([c_x,c_y,c_z],c_Ele_Domain_ID)
  
      !do i_E = 1,Num_Elem
      do i_E=1,Domain_Elements_Num(c_Ele_Domain_ID)
        c_Ele = Domain_Elements(c_Ele_Domain_ID,i_E)
        
        if(Flag_Exit==1)then
            exit
        endif
        !c_X_NODES = G_X_NODES(i_E,1:8)
        !c_Y_NODES = G_Y_NODES(i_E,1:8)  
        !c_Z_NODES = G_Z_NODES(i_E,1:8)  
        c_X_NODES = G_X_NODES(1:8,c_Ele)
        c_Y_NODES = G_Y_NODES(1:8,c_Ele)    
        c_Z_NODES = G_Z_NODES(1:8,c_Ele) 
        ! If the point is not within the element range, then exit. 2022-06-12. IMPROV2022061201.
        if(c_x > (maxval(c_X_NODES)+Tol_8)) cycle
        if(c_y > (maxval(c_Y_NODES)+Tol_8)) cycle
        if(c_z > (maxval(c_Z_NODES)+Tol_8)) cycle
        if(c_x < (minval(c_X_NODES)-Tol_8)) cycle
        if(c_y < (minval(c_Y_NODES)-Tol_8)) cycle
        if(c_z < (minval(c_Z_NODES)-Tol_8)) cycle
        
          ! The six faces of the element are divided into 12 triangles.
          Triangle(1,1:3) = [1,2,5]
          Triangle(2,1:3) = [2,6,5]
          Triangle(3,1:3) = [2,3,6]
          Triangle(4,1:3) = [3,7,6]
          Triangle(5,1:3) = [3,7,8]
          Triangle(6,1:3) = [3,8,4]
          Triangle(7,1:3) = [4,8,5]
          Triangle(8,1:3) = [4,5,1]
          Triangle(9,1:3) = [5,6,8]
          Triangle(10,1:3) = [6,7,8]
          Triangle(11,1:3) = [1,2,4]
          Triangle(12,1:3) = [2,3,4]    
          ! 12 Triangle Loop
          do i_Tri = 1,12
              Point1 =  [c_X_NODES(Triangle(i_Tri,1)),c_Y_NODES(Triangle(i_Tri,1)),c_Z_NODES(Triangle(i_Tri,1))]
              Point2 =  [c_X_NODES(Triangle(i_Tri,2)),c_Y_NODES(Triangle(i_Tri,2)),c_Z_NODES(Triangle(i_Tri,2))]
              Point3 =  [c_X_NODES(Triangle(i_Tri,3)),c_Y_NODES(Triangle(i_Tri,3)),c_Z_NODES(Triangle(i_Tri,3))] 
              ! Check if the checkpoint is on the spatial triangle
              call Tool_Yes_Point_on_3D_Triangle([c_x,c_y,c_z],Point1,Point2,Point3,c_Yes_on)
              if(c_Yes_on)then
                ! Vertex cell number
                El=Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(c_Mesh_Node) 
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ! The center coordinates of the cell containing the vertex
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                c_X_NODES_El = G_X_NODES(1:8,El)
                c_Y_NODES_El = G_Y_NODES(1:8,El)  
                c_Z_NODES_El = G_Z_NODES(1:8,El)  
                
                mid_Ele_x = sum(c_X_NODES_El(1:8))/EIG
                mid_Ele_y = sum(c_Y_NODES_El(1:8))/EIG
                mid_Ele_z = sum(c_Z_NODES_El(1:8))/EIG                   
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ! Move c_x, c_y, c_z toward the center of the cell containing the vertex by one
                ! ten-thousandth of the cell's average length
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                Line_AB(1,1:3) = [c_x,c_y,c_z]
                Line_AB(2,1:3) = [mid_Ele_x,mid_Ele_y,mid_Ele_z]
                delta_L = Ave_Elem_L/1.0D4
                ! Adjust the vertex position toward the element center
                call Tool_Shorten_or_Extend_Line_3D(Line_AB,-delta_L,'A',new_Line_AB,new_Point)
                ! Update coordinates of discrete crack edge nodes
                Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,1)=new_Point(1)
                Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,2)=new_Point(2) 
                Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,3)=new_Point(3)  
                !goto 100
                
                Flag_Exit = 1
                exit
              endif
          enddo
      enddo
    enddo
enddo
!$omp end parallel do       
RETURN
END SUBROUTINE D3_Crack_Vertex_Adjust