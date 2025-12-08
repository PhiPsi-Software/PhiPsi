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
 
      SUBROUTINE Determine_Enriched_Nodes_Cross_Check(ifra,isub,
     &                    Yes_Cross_Ele)
c     (1) Check for the presence of cross-shaped intersection cracks and cross cracks
c     (2)...

      ! Variable Description
      ! Ele_Cross_Point_Cr_num(Num_Elem,1:2) !The crack numbers of the two cracks that include the
      ! cross-point element
      ! Ele_Cross_Point_Cr_Seg_num(Num_Elem,1:2) !Contains the crack segment numbers corresponding to the
      ! two cracks of the element at the Cross point
      ! Ele_Cross_Point_RABCD(Num_Elem,1:10,1:2) !Element cross point coordinates, as well as the four
      ! intersection points A, B, C, D, R between the crack segment and the element
      ! Where R is the intersection point
      ! For details about this variable, see the notes on V5-PhiPsi cross-shaped fracture handling.



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
      use Global_Cross

c     --------------------------
c     Variable Type Declaration
c     --------------------------
      implicit none
      integer,intent(in)::ifra,isub
      logical,intent(out)::Yes_Cross_Ele
      integer i_C,i_S,i_Side,i_E
      integer j_C,j_S
      real(kind=FT) S1(2),S2(2)
      real(kind=FT) S3(2),S4(2)
      real(kind=FT) c_Inter_x,c_Inter_y
      logical c_Yes_Cross
      integer c_Elem
      integer N1,N2,N3,N4,NN(4),NN_L(5)
      integer n_node_1,n_node_2
      real(kind=FT) point_side_1(2),point_side_2(2)
      real(kind=FT) c_X,c_Y
      logical c_Yes_Cross_2
      integer c_num_Inters
      real(kind=FT) c_Inters(10,2)
      real(kind=FT) dis_S1_A,dis_S1_B,dis_S3_A,dis_S3_B
      real(kind=FT) Tool_Function_2Point_Dis
      real(kind=FT) Tiny_V
      logical Yes_T_Crack

      Tiny_V = 1.0D-5*Ave_Elem_L

c     --------------------
c     Data Initialization
c     --------------------
      Cross_Point_Cr_num(1:Max_Num_Cross,1:2)     = 0
      Cross_Point_Ele_num(1:Max_Num_Cross)    = 0
      !Cross_Point_Cr_Seg_num(1:Max_Num_Cross,1:2) = ZR
      Cross_Point_RABCD(1:Max_Num_Cross,1:10,1:2) = ZR

c     ---------------------------
c     Perform Detection Task (1)
c     ---------------------------
      num_Cross = 0
      Yes_Cross_Ele =.False.
      if(num_Crack>=2)then
        do i_C=1,num_Crack
          ! Current cracks are cycling between each crack segment
          do i_S = 1,Each_Cr_Poi_Num(i_C)-1
              ! Crack fragment endpoint coordinates
              S1 = [Crack_Coor(i_C,i_S,1),Crack_Coor(i_C,i_S,2)]
              S2 = [Crack_Coor(i_C,i_S+1,1), Crack_Coor(i_C,i_S+1,2)]

              do j_C =i_C+1,num_Crack
                  do j_S = 1,Each_Cr_Poi_Num(j_C)-1
                      ! Crack fragment endpoint coordinates
                      S3 = [Crack_Coor(j_C,j_S,1),Crack_Coor(j_C,j_S,2)]
                      S4 = [Crack_Coor(j_C,j_S+1,1),
     &                      Crack_Coor(j_C,j_S+1,2)]
                      ! Determine whether two crack segments intersect
                      call Tool_Intersection(
     &                            S1,S2,S3,S4,
     &                            c_Inter_x,c_Inter_y,c_Yes_Cross)
                      ! Determine whether it is a T-shaped intersection. A T-shaped intersection is not considered a
                      ! cross-shaped intersection, and the intersection point of a T-shaped intersection coincides with
                      ! one endpoint of the line segment.
                      Yes_T_Crack = .False.
                      if( (abs(c_Inter_x-S1(1)) < Tiny_V .and.
     &                     abs(c_Inter_y-S1(2)) < Tiny_V)     .or.
     &                    (abs(c_Inter_x-S2(1)) < Tiny_V .and.
     &                     abs(c_Inter_y-S2(2)) < Tiny_V)     .or.
     &                    (abs(c_Inter_x-S3(1)) < Tiny_V .and.
     &                     abs(c_Inter_y-S3(2)) < Tiny_V)     .or.
     &                    (abs(c_Inter_x-S4(1)) < Tiny_V .and.
     &                     abs(c_Inter_y-S4(2)) < Tiny_V))then
                          Yes_T_Crack = .True.
                      endif
                      ! If intersecting
                      if(c_Yes_Cross.and.(Yes_T_Crack.EQV..FALSE.))then
                          Yes_Cross_Ele =.True.

                          num_Cross = num_Cross + 1
                          !*****************************************************
                          ! Calculate the cell number of the intersection point
                          !*****************************************************
                          call Cal_Ele_Num_by_Coors(
     &                              c_Inter_x,c_Inter_y,c_Elem)
                          !********************
                          ! Save relevant data
                          !********************
                          Cross_Point_Cr_num(num_Cross,1) = i_C
                          Cross_Point_Cr_num(num_Cross,2) = j_C
                          Cross_Point_Ele_num(num_Cross)  = c_Elem
                          !Ele_Cross_Point_Cr_Seg_num(c_Elem,1)= i_S
                          !Ele_Cross_Point_Cr_Seg_num(c_Elem,2)= j_S
                          Cross_Point_RABCD(num_Cross,1,1:2) =
     &                                          [c_Inter_x,c_Inter_y]
                          !*************************************************************************
                          ! Calculate the intersection points of two fracture segments with the
                          ! element, A, B and C, D
                          ! It should be noted here that the order of the intersections must
                          ! correspond to
                          ! The crack segments are aligned in the same direction, namely: S1, A, B,
                          ! S2
                          !                    S3,C,D,S4
                          !*************************************************************************
      					  N1  = Elem_Node(c_Elem,1)
      	      N2  = Elem_Node(c_Elem,2)
      	      N3  = Elem_Node(c_Elem,3)
      	      N4  = Elem_Node(c_Elem,4)
      	      NN_L= [N1,N2,N3,N4,N1]
                          ! First, calculate the intersection points of crack i_C with the four boundaries of the element.
                          c_num_Inters = 0
                          c_Inters(1:10,1:2) = ZR
      	      do i_Side = 1,4
      	          n_node_1 = NN_L(i_Side)
                              n_node_2 = NN_L(i_Side+1)
      		      !Get the coordinate of the two points of the current element side.
      		      point_side_1 = [Coor(n_node_1,1),
     &                                        Coor(n_node_1,2)]
      		      point_side_2 = [Coor(n_node_2,1),
     &                                        Coor(n_node_2,2)]
                              !Check if they intersects and calculate intersection point.
                              call Tool_Intersection(
     &                                       S1,S2,
     &                                      point_side_1,point_side_2,
     &                                      c_X,c_Y,c_Yes_Cross_2)
                              if (c_Yes_Cross_2.eqv..True.) then
                                  c_num_Inters = c_num_Inters + 1
                                  c_Inters(c_num_Inters,1:2)= [c_X,c_Y]
                              end if
                          end do
                          ! Sequential detection and adjustment
                          dis_S1_A = Tool_Function_2Point_Dis(
     &                                       S1,c_Inters(1,1:2))
                          dis_S1_B = Tool_Function_2Point_Dis(
     &                                       S1,c_Inters(2,1:2))
                          if(dis_S1_A <=dis_S1_B) then
                              Cross_Point_RABCD(num_Cross,2,1:2) =
     &                                                 c_Inters(1,1:2)
                              Cross_Point_RABCD(num_Cross,3,1:2) =
     &                                                 c_Inters(2,1:2)
                          else
                              Cross_Point_RABCD(num_Cross,2,1:2) =
     &                                                 c_Inters(2,1:2)
                              Cross_Point_RABCD(num_Cross,3,1:2) =
     &                                                 c_Inters(1,1:2)
                          endif
                          ! First, calculate the intersection points of crack i_C with the four boundaries of the element.
                          c_num_Inters = 0
                          c_Inters(1:10,1:2) = ZR
      	      do i_Side = 1,4
      	          n_node_1 = NN_L(i_Side)
                              n_node_2 = NN_L(i_Side+1)
      		      !Get the coordinate of the two points of the current element side.
      		      point_side_1 = [Coor(n_node_1,1),
     &                                        Coor(n_node_1,2)]
      		      point_side_2 = [Coor(n_node_2,1),
     &                                        Coor(n_node_2,2)]
                              !Check if they intersects and calculate intersection point.
                              call Tool_Intersection(
     &                                       S3,S4,
     &                                      point_side_1,point_side_2,
     &                                      c_X,c_Y,c_Yes_Cross_2)
                              if (c_Yes_Cross_2.eqv..True.) then
                                  c_num_Inters = c_num_Inters + 1
                                  c_Inters(c_num_Inters,1:2)= [c_X,c_Y]
                              end if
                          end do
                          ! Sequential detection and adjustment
                          dis_S3_A = Tool_Function_2Point_Dis(
     &                                       S3,c_Inters(1,1:2))
                          dis_S3_B = Tool_Function_2Point_Dis(
     &                                       S3,c_Inters(2,1:2))
                          if(dis_S3_A <=dis_S3_B) then
                              Cross_Point_RABCD(num_Cross,4,1:2) =
     &                                                 c_Inters(1,1:2)
                              Cross_Point_RABCD(num_Cross,5,1:2) =
     &                                                 c_Inters(2,1:2)
                          else
                              Cross_Point_RABCD(num_Cross,4,1:2) =
     &                                                 c_Inters(2,1:2)
                              Cross_Point_RABCD(num_Cross,5,1:2) =
     &                                                 c_Inters(1,1:2)
                          endif
                      endif
                  enddo
              end do
          enddo
        end do
      endif

      !**********************************************************************************
      ! If a cross-shaped crack is detected, pop up a message and terminate the program.
      !**********************************************************************************
      if (Yes_Cross_Ele) then
          print *, '    Error :: crossing cracks are not allowed '
     &              // 'in 2D PhiPsi!'
          call Warning_Message('S',Keywords_Blank)
      endif

      RETURN
      END SUBROUTINE Determine_Enriched_Nodes_Cross_Check
