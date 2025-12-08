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
 
      SUBROUTINE Determine_Enriched_Nodes_Pre_Check(ifra,isub)
c     (1) Correct the crack coordinates to prevent certain difficult-to-handle situations.
c     The main situation that needs to be prevented is: preventing a single element edge from having two intersections with a crack segment.
c     Corresponding notes can be found in: V5-P49
c     (2)...

c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Common   
      use Global_Model
      use Global_Elem_Area_Vol
      use Global_Crack
      use Global_Crack_Common
      use Global_HF
      
c     --------------------------
c     Variable Type Declaration
c     --------------------------
      implicit none
      integer ifra,isub,i_C,i_S,i_Side,i_E
      real(kind=FT) c_X,c_Y,c_Line_AB(2,2)
      logical c_Yes_Cross,Flag_1,Flag_2
      integer N1,N2,N3,N4,NN(4),NN_L(5),n_node_1,n_node_2,Num_Intersect
      real(kind=FT) point_side_1(2),point_side_2(2),
     &                 point_seg_1(2),point_seg_2(2),
     &                 X_NODES_L(5),Y_NODES_L(5)
      integer i_inter,i_Segment,Crack_Seg_Num(5)
      integer j_C,j_S
      logical c_Flag_1,c_Flag_2
      real(kind=FT) point_seg2_1(2),point_seg2_2(2)
      logical c_Yes_Cross2
      
      real(kind=FT) point_V(2),point_M(2),point_N(2)
      real(kind=FT) point_A(2),point_B(2),Line_VM(2,2)
      real(kind=FT) delta_L,Line_VN(2,2)
      
c     ---------------------------
c     Perform Detection Task (1)
c     ---------------------------
      do i_C=1,num_Crack
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Cycle between each element
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do i_E = 1,Num_Elem
              N1  = Elem_Node(i_E,1)                                         
              N2  = Elem_Node(i_E,2)                                             
              N3  = Elem_Node(i_E,3)                                             
              N4  = Elem_Node(i_E,4)  
              NN  = [N1,N2,N3,N4]
              NN_L= [N1,N2,N3,N4,N1]
              X_NODES_L = Coor(NN_L,1)
              Y_NODES_L = Coor(NN_L,2)
              Num_Intersect = 0
      		  Flag_1 = .False.
              Flag_2 = .False.
              !*************************************************
              ! Current cracks cycle between each crack segment
              !*************************************************
              do i_S = 1,Each_Cr_Poi_Num(i_C)-1
                  !Get the coordinate of the two points of the current segment.
                  point_seg_1=[Crack_Coor(i_C,i_S,1),
     &                         Crack_Coor(i_C,i_S,2)]
                  point_seg_2=[Crack_Coor(i_C,i_S+1,1),
     &                         Crack_Coor(i_C,i_S+1,2)] 
                  ! Adjacent crack segment
                  if ((Each_Cr_Poi_Num(i_C)-1 >= 2) .and. 
     &                (i_S  < (Each_Cr_Poi_Num(i_C)-1))) then
                      point_seg2_1=[Crack_Coor(i_C,i_S+1,1),
     &                              Crack_Coor(i_C,i_S+1,2)]
                      point_seg2_2=[Crack_Coor(i_C,i_S+2,1),
     &                              Crack_Coor(i_C,i_S+2,2)]                       
                  end if
                  ! Loop through the four edges of the current element
                  do i_Side = 1,4
                      n_node_1 = NN_L(i_Side)
                      n_node_2 = NN_L(i_Side+1)
                      !Get the coordinate of the two points of the current side.
                      point_side_1 =[Coor(n_node_1,1),Coor(n_node_1,2)]
                      point_side_2 =[Coor(n_node_2,1),Coor(n_node_2,2)]
       	              !Check if they intersects and calculate intersection point.
                      call Tool_Intersection(
     &                              point_seg_1,point_seg_2,
     &                              point_side_1,point_side_2,
     &                              point_A(1),point_A(2),c_Yes_Cross) 
                      ! If the current element boundary intersects with the current crack segment
                      if(c_Yes_Cross.eqv..True.)then
                          ! Continue checking for intersections between adjacent crack segments and the current cell
                          ! boundaries
                          if(i_S  < (Each_Cr_Poi_Num(i_C)-1))then
                              call Tool_Intersection(
     &                               point_seg2_1,point_seg2_2,
     &                               point_side_1,point_side_2,
     &                               point_B(1),point_B(2),c_Yes_Cross2)  
                              ! If there is an intersection, the coordinates of the crack need to be corrected.
                              if(c_Yes_Cross2)then
                                  ! Midpoint M of AB
                                  point_M = HLF*(point_A + point_B)
                                  ! V-line
                                  point_V = point_seg2_1
                                  ! Extend VM by a tiny amount to get point N
                                  Line_VM(1,1:2) = point_V 
                                  Line_VM(2,1:2) = point_M
                                  delta_L = 1.0D-2*Ave_Elem_L
                                  call Tool_Shorten_or_Extend_Line(
     &                                      Line_VM,delta_L,
     &                                      'B',Line_VN,point_N)
                                  ! Update crack coordinates
                                  Crack_Coor(i_C,i_S+1,1:2)=point_N
                              endif
                          endif
                      else
       					  !Check if point_seg_1 is inside the current element.
                          call Tool_Yes_In_Poly
     &                              (point_seg_1(1),point_seg_1(2),
     &                               X_NODES_L,Y_NODES_L,5,c_Flag_1)
                            !Check if point_seg_2 is inside the current element.
                          call Tool_Yes_In_Poly
     &                              (point_seg_2(1),point_seg_2(2),
     &                               X_NODES_L,Y_NODES_L,5,c_Flag_2)
 
                          ! A complete crack segment is not allowed within a single unit. Check!
      		  if (c_Flag_1 .and. c_Flag_2) then
                              print *,"    Error :: A complete crack"
     &                            // " segment found in one element."
                              print *,"             Message produced in"
     &                        // " Determine_Enriched_Nodes_Pre_Check."
                              ! Exit program
                              call Warning_Message('S',Keywords_Blank)
      		  end if
                      end if
                  end do
              end do
          end do
      end do
      
      
      RETURN
      END SUBROUTINE Determine_Enriched_Nodes_Pre_Check
