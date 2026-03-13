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
 
      subroutine Tool_Yes_Ele_Intersected_by_3D_Crack_Outline(c_E,
     &                          c_C,c_Yes_Inter,
     &                          c_Inter_Point_A,c_Inter_Point_B,
     &                          Vertex_1,Vertex_2,
     &                          num_Inter_Outline,Four_Points)
C     Determine whether a 3D element is intersected by the outer boundary of a 3D discrete crack, and calculate the coordinates of the intersection points.
c     Note: It is directly penetrated by an edge line and does not include vertices on the boundaries of discrete cracks.
c     Firstly written on 2020-01-08.
c     If there are multiple crossover outlines, then take the longer one (2021-08-21), Ref:\theory_documents\023
c     The situation and handling of a 3D crack front with a zigzag edge and two intersection lines in an element - 2021-08-22.

      !......................
      ! Variable Declaration
      !......................
      use Global_Float_Type
      use Global_Crack_3D
      use Global_Model
      
      implicit none
      integer,intent(in)::c_E,c_C
      logical,intent(out):: c_Yes_Inter
      real(KIND=FT),intent(out):: c_Inter_Point_A(3),c_Inter_Point_B(3)
      real(KIND=FT),intent(out):: Four_Points(4,3)
      integer,intent(out)::Vertex_1,Vertex_2
      integer,intent(out)::num_Inter_Outline
      integer i_Outline,num_Vertex,c_Mesh_Node,c_Mesh_Nex_Node
      real(KIND=FT) c_P_1(3),c_P_2(3)
      logical c_Outline_Yes_Inter
      integer c_num_Inter
      real(KIND=FT) c_InterSection_P(2,3)
      real(KIND=FT) Tool_Function_2Point_Dis_3D,last_Length,c_Length
      real(KIND=FT) x_min_c_P,x_max_c_P
      real(KIND=FT) y_min_c_P,y_max_c_P
      real(KIND=FT) z_min_c_P,z_max_c_P
      real(KIND=FT) x_min_Ele,x_max_Ele
      real(KIND=FT) y_min_Ele,y_max_Ele
      real(KIND=FT) z_min_Ele,z_max_Ele
      
      
      c_Yes_Inter =.False.
      num_Vertex = Crack3D_Meshed_Outline_num(c_C)  
      num_Inter_Outline = 0
      last_Length = -TEN_10
      Four_Points(1:4,1:3) = ZR
      
      x_min_Ele = minval(G_X_NODES(1:8,c_E))
      x_max_Ele = maxval(G_X_NODES(1:8,c_E))
      y_min_Ele = minval(G_Y_NODES(1:8,c_E))
      y_max_Ele = maxval(G_Y_NODES(1:8,c_E))
      z_min_Ele = minval(G_Z_NODES(1:8,c_E))
      z_max_Ele = maxval(G_Z_NODES(1:8,c_E))

      do i_Outline = 1,num_Vertex
          c_Mesh_Node = Crack3D_Meshed_Outline(c_C)%row(i_Outline,1)    
          c_Mesh_Nex_Node = Crack3D_Meshed_Outline(c_C)%row(i_Outline,2)   
          c_P_1 = Crack3D_Meshed_Node(c_C)%row(c_Mesh_Node,1:3) 
          c_P_2 = Crack3D_Meshed_Node(c_C)%row(c_Mesh_Nex_Node,1:3)
          
          x_min_c_P = min(c_P_1(1),c_P_2(1)) 
          if((x_min_c_P-Tol_8) > x_max_Ele) cycle
          
          x_max_c_P = max(c_P_1(1),c_P_2(1)) 
          if((x_max_c_P+Tol_8) < x_min_Ele) cycle
          
          y_min_c_P = min(c_P_1(2),c_P_2(2)) 
          if((y_min_c_P-Tol_8) > y_max_Ele) cycle
          
          y_max_c_P = max(c_P_1(2),c_P_2(2)) 
          if((y_max_c_P+Tol_8) < y_min_Ele) cycle
          
          z_min_c_P = min(c_P_1(3),c_P_2(3)) 
          if((z_min_c_P-Tol_8) > z_max_Ele) cycle
          
          z_max_c_P = max(c_P_1(3),c_P_2(3)) 
          if((z_max_c_P+Tol_8) < z_min_Ele) cycle
          
          call Tool_All_Intersections_of_AB_and_Brick_Ele(
     &                 c_P_1,c_P_2,c_E,
     &                 c_Outline_Yes_Inter,c_num_Inter,c_InterSection_P) 
          if (c_Outline_Yes_Inter .and. c_num_Inter==2)then
              c_Length = Tool_Function_2Point_Dis_3D(
     &                              c_InterSection_P(1,1:3),
     &                              c_InterSection_P(2,1:3))
              c_Yes_Inter =.True.
              num_Inter_Outline = num_Inter_Outline +1
              
              if(num_Inter_Outline==1) then
                  Four_Points(1,1:3) = c_InterSection_P(1,1:3)
                  Four_Points(2,1:3) = c_InterSection_P(2,1:3)
              elseif(num_Inter_Outline==2) then
                  Four_Points(3,1:3) = c_InterSection_P(1,1:3)
                  Four_Points(4,1:3) = c_InterSection_P(2,1:3)                  
              endif
              if (c_Length > last_Length)then
                  c_Inter_Point_A = c_InterSection_P(1,1:3)
                  c_Inter_Point_B = c_InterSection_P(2,1:3)
                  Vertex_1 = i_Outline
                  if(i_Outline == num_Vertex) then
                      Vertex_2 = 1
                  else
                      Vertex_2 = i_Outline+1
                  endif
              endif
              last_Length = c_Length
          endif
      enddo

      return 
      end SUBROUTINE Tool_Yes_Ele_Intersected_by_3D_Crack_Outline              
