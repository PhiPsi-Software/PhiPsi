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
 
      SUBROUTINE Determine_Enriched_Nodes(ifra,isub)
c     This subroutine is used to determine the enhancement nodes.

c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Common   
      use Global_Filename
      use Global_Model
      use Global_Elem_Area_Vol
      use Global_Crack
      use Global_Crack_Common
      use Global_HF
      use Global_Inclusion
      use Global_Cross
      
c     --------------------------
c     Variable Type Declaration
c     --------------------------
      implicit none
      integer ifra,isub,i_C,i_S,i_O,i_Side,i_E,i_N,c_Edge_Elem
      logical Yes_in_Outline1,Yes_in_Outline2
      real(kind=FT) crack_p1(2),crack_p2(2),
     &                 Outline_p1(2),Outline_p2(2)
      integer EndByEnd_Outline(size(Outline,1)+1) 
      real(kind=FT) c_X,c_Y,c_Line_AB(2,2),shorted_Line_AB(2,2),
     &                 new_Point(2)
      logical c_Yes_Cross,Flag_1,Flag_2
      integer N1,N2,N3,N4,NN(4),NN_L(5),n_node_1,n_node_2,Num_Intersect
      real(kind=FT) point_side_1(2),point_side_2(2),
     &                 X_el_inter,Y_el_inter,crack_inter(5,2),
     &                 point_seg_1(2),point_seg_2(2),
     &                 X_NODES_L(5),Y_NODES_L(5)
      integer point_seg_inElemnt_Flag(Num_Elem)
      real(kind=FT) point_seg_inElemnt(Num_Elem,2)
      integer i_inter,i_Segment,Crack_Seg_Num(5),
     &        Uniqued_Cr_Seg_Num(5),Uniqued_Num_Cr_Seg
      real(kind=FT) x_inter,y_inter
      logical Yes_on_Segment,Yes_On_Line
      integer j_C,j_S,tem_i_1,int_Flag_1,int_Flag_2
      logical Yes_Ap_has_GP,Yes_Am_has_GP,Yes_Cutthrough
      integer DOMAIN_Outline(20,2),m_DOMAIN_Outline,
     &        Domain_El(10),n_Domain_El
      real(kind=FT) tem_01,AB(2,2),CD(2,2),tem_02(2),tem_03(2),
     &                 tem_point_seg_1(2),tem_point_seg_2(2)
      real(kind=FT) Area_Enrich
      integer iEl_contour
      integer c_NN(4),c_N,j_E,c_Adj_Ele
      ! -----------Hole-related Variables------------
      integer i_H,num_in_Nodes
      real(kind=FT) c_Hole_x,c_Hole_y,c_Hole_r,c_Dis,c_N_x,c_N_y
      logical Nodes_4_in(4)
      real(kind=FT) Tool_Function_2Point_Dis
      !-----
      integer i_Su_Ele,c_Chk_Ele,i_Chk_Node,c_Chk_N,Positive_Cr
      integer Negative_Cr
      logical Yes_Posit_Cross,Yes_Negti_Cross
      integer tem_Checked_Nodes(Num_Node)
      real(kind=FT) Clo_Supp_Domain(20,2)
      ! -----Mixed Related----------------------
      integer i_Incl,i_Edge,n_Vertex
      real(kind=FT) c_Incl_x,c_Incl_y,c_Incl_r
      ! -----Cross-Shaped Intersection Cracks Related------------
      logical Yes_Cross_Ele
      integer i_Cross,c_Cross_Ele
      integer c_Cr,c_ele
      integer  iii
      ! ------Crack and Hole Intersection Related
      real(kind=FT) x0_Hole,y0_Hole,R0_Hole,theta_Hole,a_Hole,b_Hole
      integer i_Tip
      real(kind=FT) Coor_Tip(2)
      integer c_Jun_Ele
      integer c_num_Inter,c_State,c_num_Inter_CH
      real(kind=FT) c_Inter(2,2)
      logical Yes_Hole_Cross
      integer c_Hole
      logical c_Yes_Arc_Crack
      integer c_num_Inter_Line_Arc
      real(kind=FT) c_Inter_Line_Arc(2,2)
      integer Return_Statu
      ! -------Arc-shaped Crack Related
      real(kind=FT) cc_tip_x,cc_tip_y
      integer num_points_crack
      integer i_El_Check,c_Check_El
      real(kind=FT) c_Point_1(2),c_Point_2(2)
      logical Yes_P1_ON,Yes_P2_ON
      logical Yes_P_Equal_1,Yes_P_Equal_2
      integer i_P,c_OUT_Elem
      real(kind=FT) c_tip_x,c_tip_y,check_L
      integer i_Node
      real(kind=FT) check_x_min,check_x_max,check_y_min,check_y_max
      integer c_ELe_Num 

c     ---------------------------------------------------------------------------------------
c     Initialize the matrix elements to 0. Note:
c     Some variables are only initialized in the first rupture step, while others need to be
c     initialized
c     in any rupture step.
c     ---------------------------------------------------------------------------------------
      ! All variables initialized at the point of each break step
      Crack_Tip_Coor(1:Max_Num_Cr,1:2,1:2)             = ZR
      Elem_Type(1:Num_Elem,1:Max_Num_Cr)               = 0
      Enriched_Node_Type(1:Num_Node,1:Max_Num_Cr)      = 0
      Node_Jun_elem(1:Num_Node,1:Max_Num_Cr)           = 0
      Jun_Ele_Negative_Cr_Num(1:Num_Elem,1:Max_Num_Cr) = 0
      Coors_Element_Crack(1:Num_Elem,1:Max_Num_Cr,1:4) = ZR
      Coors_Tip(1:Num_Elem,1:2)                        = ZR
      Coors_Vertex(1:Num_Elem,1:2)                     = ZR
      Coors_Junction(1:Num_Elem,1:Max_Num_Cr,1:4)      = ZR
      x_cr_tip_nodes(1:Max_Num_Cr,1:Num_Node)          = ZR
      y_cr_tip_nodes(1:Max_Num_Cr,1:Num_Node)          = ZR
      Ele_Num_Tip_Enriched_Node(1:Max_Num_Cr,1:Num_Node)   = 0
      TipEle_Adjacent_Ele(1:Num_Elem,1:Max_Num_Cr)         = 0
      !Crack_Tip_Type(1:Num_Crack,1:2)   = 0
      Edge_Disposed_Crack(1:Max_Num_Cr,1:Max_Num_Cr_P,1:2) = ZR     
      Crack_Jun_CrNum(1:Max_Num_Cr,1:2)          = 0
      Crack_Jun_HoleNum(1:Max_Num_Cr,1:2)        = 0
      Crack_Jun_Elem(1:Max_Num_Cr,1:2)           = 0  
      Node_Jun_Hole(1:Num_Node,1:Max_Num_Cr)     = 0 
      Ele_Jun_Hole(1:Num_Elem,1:Max_Num_Cr)      = 0
      ! Variables initialized only in the first iteration of the loop
      if (ifra==1) then
          Crack_Tip_Type(1:Max_Num_Cr,1:2) = 0
          ! Related to holes
          Elem_Type_Hl(1:Num_Elem,1:Max_Num_Hl) = 0
          Enriched_Node_Type_HL(1:Num_Node,1:Max_Num_Hl)=0
          ! Mixed-related initialization
          Elem_Type_Incl(1:Num_Elem,1:Max_Num_Incl) = 0
          Enriched_Node_Type_Incl(1:Num_Node,1:Max_Num_Incl)=0
      end if 
      ! Cross-related data initialization (this plan was later abandoned)
      Node_Cross_elem(1:Num_Node,1:Max_Num_Cross)  = 0
      Elem_Type_Cross(1:Num_Elem,1:Max_Num_Cross)  = 0
      Enriched_Node_Type_Cross(1:Num_Node,1:Max_Num_Cross)= 0
      ! Arc Crack Data Initialization (2017-07-18)
      Crack_Arc_Tip_A_B_C_x(1:Max_Num_Cr,1:2,1:3)  = ZR
      Crack_Arc_Tip_A_B_C_y(1:Max_Num_Cr,1:2,1:3)  = ZR

c     ---------------------------
c     Formatted output statement
c     ---------------------------
 1001 FORMAT('     Caution :: tip 1 of crack ',I3,' is edge crack.') 
 1002 FORMAT('     Caution :: tip 2 of crack ',I3,' is edge crack.')     
 
      print *,'    Determine enriched node...'
      
c     -----------------------------------------------------------------------------------------------
c     ---------------------         Step 0: Conduct Preliminary Check          ---------------------
c     -----------------------------------------------------------------------------------------------
      ! Pre-detection and coordinate correction.
      if (Key_User_Defined_2D_Crack_Path/=1)then
          call Determine_Enriched_Nodes_Pre_Check(ifra,isub)
      endif
      
      ! Detect cross-shaped cracks.
      call Determine_Enriched_Nodes_Cross_Check(ifra,isub,Yes_Cross_Ele)
     
c     -----------------------------------------------------------------------------------
c     ---------------------   Step 1: determine the element type   ---------------------
c     -----------------------------------------------------------------------------------
      Edge_Disposed_Crack = Crack_Coor
      loop_Crack: do i_C=1,num_Crack
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !Initialization of some variables
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          c_Edge_Elem = 0
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Current cracks cycle between each crack segment
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          loop_i_S: do i_S = 1,Each_Cr_Poi_Num(i_C)-1
              !*************************************
              ! Crack fragment endpoint coordinates
              !*************************************
              crack_p1 = [Crack_Coor(i_C,i_S,1),Crack_Coor(i_C,i_S,2)]
              crack_p2 = [Crack_Coor(i_C,i_S+1,1),
     &                    Crack_Coor(i_C,i_S+1,2)]
              !************************************************************************************
              ! Extract the crack tip coordinates and save them to the global crack tip coordinate
              ! variable
              !************************************************************************************
              if (i_S.eq.1) then
                  Crack_Tip_Coor(i_C,1,1:2) = crack_p1
                  ! Save to tip coordinate global variable
                  Cr_First_Tip(i_C,1:2)= crack_p1
                  ! Save to the global variable of dip angle for the fracture segment where the crack tip is located
                  Cr_First_Tip_Ori(i_C)=atan2(crack_p1(2)-crack_p2(2),
     &                                        crack_p1(1)-crack_p2(1))
              end if
              if(i_S.eq.Each_Cr_Poi_Num(i_C)-1) then
                  Crack_Tip_Coor(i_C,2,1:2) = crack_p2
                  ! Save to tip coordinate global variable
                  Cr_Second_Tip(i_C,1:2)= crack_p2
                  ! Save to the global variable of dip angle for the fracture segment where the crack tip is located
                  Cr_Second_Tip_Ori(i_C)=atan2(crack_p2(2)-crack_p1(2),
     &                                         crack_p2(1)-crack_p1(1))
              end if
              !***************************************************************************
              ! The Outline dimension is num x 2, convert it into a one-dimensional array
              ! EndByEnd_Outline
              !***************************************************************************
              EndByEnd_Outline(1:size(Outline,1)) = Outline(:,1)
              EndByEnd_Outline(size(Outline,1)+1) = Outline(1,1)        
              !********************************************************************
              !********************************************************************
              Call Tool_Yes_In_Poly
     &             (crack_p1(1),crack_p1(2),
     &              Coor(EndByEnd_Outline,1),Coor(EndByEnd_Outline,2),
     &              size(EndByEnd_Outline),Yes_in_Outline1)
              Call Tool_Yes_In_Poly
     &             (crack_p2(1),crack_p2(2),
     &              Coor(EndByEnd_Outline,1),Coor(EndByEnd_Outline,2),
     &              size(EndByEnd_Outline),Yes_in_Outline2)
              !*************************************************************************
              ! Find the intersection between the current crack segment and the outline
              !*************************************************************************
              do i_O = 1,size(Outline,1)
                 Outline_p1 = Coor(Outline(i_O,1),:)
                 Outline_p2 = Coor(Outline(i_O,2),:)   
                 call Tool_Intersection(crack_p1,crack_p2,
     &                                      Outline_p1,Outline_p2,
     &                                      c_X,c_Y,c_Yes_Cross)  
                 !**********************************************************
                 ! If the current crack segment intersects with the outline
                 !**********************************************************
                  if(c_Yes_Cross.eqv..true.) then
                      call Cal_Ele_Num_by_Coors(c_X,c_Y,c_Edge_Elem)
                      !##########################
                      ! The first crack fragment
                      !##########################
                      if ((i_S .eq. 1).and.
     &                    (Yes_in_Outline1.eqv..False.))then
                          ! Crack tip 1 is an edge crack
                          WRITE(*,1001) i_C
                          Crack_Tip_Type(i_C,1) = -2
                          Flag_Crack_Tip_Out_Mol(i_C,1) = .True.
                          ! Special treatment for edge cracks, modify coordinates
                          c_Line_AB(1,:) = [c_X,c_Y]
                          c_Line_AB(2,:) = [crack_p2]
                          ! Shorten Line_AB
                          
                          !tem_01 = -Delta_Factor_Edge*Ave_Elem_L_Enrich
                          
                          ! BUGFIX2024091412. At this point, Ave_Elem_L_Enrich has not been defined.
                          tem_01 = -Delta_Factor_Edge*Ave_Elem_L
                          
                          call Tool_Shorten_or_Extend_Line
     &                              (c_Line_AB,tem_01,'A',
     &                              shorted_Line_AB,new_Point)
                          ! Cracks for post-processing
                          Edge_Disposed_Crack(i_C,1,1:2)=new_Point
                          ! Calculate Coors_Element_Crack
                          !Get the intersection point of the shorted first crack segment 
                          !'shorted_Line_AB' and the 4 element sides.
      					  N1  = Elem_Node(c_Edge_Elem,1)                                         
      	      N2  = Elem_Node(c_Edge_Elem,2)                                             
      	      N3  = Elem_Node(c_Edge_Elem,3)                                             
      	      N4  = Elem_Node(c_Edge_Elem,4)                                            
      	      NN_L= [N1,N2,N3,N4,N1]
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
     &                        shorted_Line_AB(1,:),shorted_Line_AB(2,:),
     &                        point_side_1,point_side_2,
     &                        c_X,c_Y,c_Yes_Cross) 
                              if (c_Yes_Cross.eqv..True.) then
                                  X_el_inter = c_X
                                  Y_el_inter = c_Y
                                  exit
                              end if
                          end do  
                          Coors_Element_Crack(c_Edge_Elem,i_C,1:4) 
     &                               = [c_X,c_Y,X_el_inter,Y_el_inter]
                      !#########################
                      ! The last crack fragment
                      !#########################
                      elseif ((i_S .eq. Each_Cr_Poi_Num(i_C)-1).and.
     &                    (Yes_in_Outline2.eqv..False.))then 
                          ! Tip 2 is an edge crack
                          WRITE(*,1002) i_C
                          Crack_Tip_Type(i_C,2) = -2
                          Flag_Crack_Tip_Out_Mol(i_C,2) = .True.
                          ! Special treatment for edge cracks, modify coordinates
                          c_Line_AB(1,:) = [crack_p1]
                          c_Line_AB(2,:) = [c_X,c_Y]
                          ! Shorten Line_AB
                          !tem_01 = -Delta_Factor_Edge*Ave_Elem_L_Enrich
                          
                          ! BUGFIX2024091412. At this point, Ave_Elem_L_Enrich has not been defined yet.
                          tem_01 = -Delta_Factor_Edge*Ave_Elem_L
                          
                          call Tool_Shorten_or_Extend_Line
     &                              (c_Line_AB,tem_01,'B',
     &                              shorted_Line_AB,new_Point)
                          ! Cracks for post-processing
                          Edge_Disposed_Crack(i_C,
     &                                 Each_Cr_Poi_Num(i_C),:)=new_Point
                          ! Calculate Coors_Element_Crack
                          !Get the intersection point of the shorted first crack segment 
                          !'shorted_Line_AB' and the 4 element sides.
      					  N1  = Elem_Node(c_Edge_Elem,1)                                         
      	      N2  = Elem_Node(c_Edge_Elem,2)                                             
      	      N3  = Elem_Node(c_Edge_Elem,3)                                             
      	      N4  = Elem_Node(c_Edge_Elem,4)                                            
      	      NN_L= [N1,N2,N3,N4,N1]
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
     &                        shorted_Line_AB(1,:),shorted_Line_AB(2,:),
     &                        point_side_1,point_side_2,
     &                        c_X,c_Y,c_Yes_Cross) 
                              if (c_Yes_Cross.eqv..True.) then
                                  X_el_inter = c_X
                                  Y_el_inter = c_Y
                                  exit
                              end if
                          end do  
                          Coors_Element_Crack(c_Edge_Elem,i_C,1:4) 
     &                               = [c_X,c_Y,X_el_inter,Y_el_inter] 
                      end if
                  end if
                  
              end do
          end do loop_i_S
          !**********************************************************
          ! If there is an inline, that is, a round hole, then......
          ! Content to be completed
          !**********************************************************
          !           ToDoList, TDL
          !************************************

      
          !!!!!!!!!!!!!!!!!!!!!!!!!
          ! Cycle between elements
          !!!!!!!!!!!!!!!!!!!!!!!!!
          loop_Ele: do i_E = 1,Num_Elem
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
              !*******************************************************
              ! Current cracks are cycling between each crack segment
              !*******************************************************
              loop_i_S_2: do i_S = 1,Each_Cr_Poi_Num(i_C)-1
                  ! Is the current crack segment an arc segment (2017-07-16)
                  c_Yes_Arc_Crack = .False.
                  if(sum(abs(Arc_Crack_Coor(i_C,i_S,1:11))) > ZR)then
                      c_Yes_Arc_Crack = .True.
                      ! Change the crack tip type to 7: curved crack tip (2017-07-20)
                      if(i_S == 1)then 
                          Crack_Tip_Type(i_C,1) = 7
                      elseif(i_S == Each_Cr_Poi_Num(i_C)-1)then 
                          Crack_Tip_Type(i_C,2) = 7
                      endif
                  endif

                  
                  !Get the coordinate of the two points of the current segment.
                  point_seg_1=[Crack_Coor(i_C,i_S,1),
     &                         Crack_Coor(i_C,i_S,2)]
                  point_seg_2=[Crack_Coor(i_C,i_S+1,1),
     &                         Crack_Coor(i_C,i_S+1,2)] 
                  ! Loop through the four edges of the current element
                  do i_Side = 1,4
                      n_node_1 = NN_L(i_Side)
                      n_node_2 = NN_L(i_Side+1)
                  !Get the coordinate of the two points of the current side.
                     point_side_1 =[Coor(n_node_1,1),Coor(n_node_1,2)]
                   point_side_2 =[Coor(n_node_2,1),Coor(n_node_2,2)]            
                      !Check if they intersects and calculate intersection point.
                      if(c_Yes_Arc_Crack .EQV. .False.) then
                          call Tool_Intersection(                   
     &                                       point_seg_1,point_seg_2,
     &                                       point_side_1,point_side_2,
     &                                       c_X,c_Y,c_Yes_Cross) 
                      elseif(c_Yes_Arc_Crack .EQV. .True.) then
                          call Tool_Intersection_Line_and_Arc(
     &                                   point_side_1,point_side_2,
     &                                   Arc_Crack_Coor(i_C,i_S,1:11),
     &                                   c_num_Inter_Line_Arc,
     &                                   c_Inter_Line_Arc,
     &                                   c_Yes_Cross) 
                          ! If a element edge has two intersections with the current inclusion, it is considered to have no
                          ! intersection.
                          ! The principle behind this can be seen in the upper left diagram on page P118 of notes V5-58;
                          ! debugging on 2017-07-23
                          if(c_num_Inter_Line_Arc==2)then
                              c_Yes_Cross = .False.
                          endif
                          !if(c_Yes_Cross)then
                          !endif
                      endif
                      if(c_Yes_Cross.eqv..True.)then
                          Num_Intersect  = Num_Intersect + 1     
                          if(c_Yes_Arc_Crack .EQV. .False.) then
                              crack_inter(Num_Intersect,1:2) =[c_X,c_Y]
                          elseif(c_Yes_Arc_Crack .EQV. .True.) then 
                              crack_inter(Num_Intersect,1:2) =
     &                                   c_Inter_Line_Arc(1,1:2)
                          endif
       					  !Check if point_seg_1 is inside the current element.
                          call Tool_Yes_In_Poly
     &                              (point_seg_1(1),point_seg_1(2),
     &                               X_NODES_L,Y_NODES_L,5,Flag_1)
      	      !Check if point_seg_2 is inside the current element.
                          call Tool_Yes_In_Poly
     &                              (point_seg_2(1),point_seg_2(2),
     &                               X_NODES_L,Y_NODES_L,5,Flag_2)
      		      !Get the coordinates of the point of the crack segment 
                          !in the element, actually, it's point_seg_1 or point_seg_2. 
                          if (Flag_1.eqv..True.) then
                              int_Flag_1 = 1
                          else
                              int_Flag_1 = 0
                          end if
                          if (Flag_2.eqv..True.) then
                              int_Flag_2 = 1
                          else
                              int_Flag_2 = 0
                          end if
      	      point_seg_inElemnt(i_E,:) = 
     &                          point_seg_1*int_Flag_1 + 
     &                          point_seg_2*int_Flag_2
      	      if (Flag_1.eqv..True.) then
      	          point_seg_inElemnt_Flag(i_E) = 1
      	      elseif (Flag_2.eqv..True.) then
      	          point_seg_inElemnt_Flag(i_E) = 2
      	      end if
                      end if
                  end do
              end do loop_i_S_2
              
              !*************************************************
              ! Next, determine the type of enhancement element
              !*************************************************
              !###################################################
              !(Type 2):Fully cracked element without kink point.
              !###################################################
              if((Num_Intersect.eq.2).
     &                 and.(Flag_1.eqv..False.).
     &                 and.(Flag_2.eqv..False.))then
                  Elem_Type(i_E,i_C) = 2
                  !Check if the order of the two points of "crack_inter" need to be changed.
                  !Make sure that the order of "crack_inter" are the same as the crack.
            !Loop to find how many crack segments have "crack_inter" points, actually 1.
                  tem_i_1 = 0
                  Crack_Seg_Num(1:5) =0
                  do i_inter=1,2
                      x_inter = crack_inter(i_inter,1)
                      y_inter = crack_inter(i_inter,2)
                      do i_Segment = 1,Each_Cr_Poi_Num(i_C)-1
                          point_seg_1 = Crack_Coor(i_C,i_Segment,1:2)
                          point_seg_2 = Crack_Coor(i_C,i_Segment+1,1:2)
                          ! Is the current crack segment an arc segment (2017-07-16)
                          c_Yes_Arc_Crack = .False.
                          if(sum(abs(Arc_Crack_Coor(
     &                                 i_C,i_Segment,1:11)))> ZR)then
                              c_Yes_Arc_Crack = .True.
                          endif
                          if(c_Yes_Arc_Crack .EQV. .False.) then
                              call Tool_Yes_On_Line(x_inter,y_inter,
     &                          point_seg_1,point_seg_2,Yes_on_Segment)
                          elseif(c_Yes_Arc_Crack)then
                              call Tool_Yes_On_Arc(x_inter,y_inter,
     &                             Arc_Crack_Coor(i_C,i_Segment,1:11),
     &                             Yes_on_Segment)
                          endif
      	      if (Yes_on_Segment) then
                                tem_i_1 = tem_i_1 + 1
                                Crack_Seg_Num(tem_i_1)  = i_Segment
      	      end if
                      end do                     
                  end do
                  ! Remove duplicate elements from the first few elements in the array and store them in
                  ! Uniqued_Cr_Seg_Num.
                  ! Among them, the valid data of the newly generated array is from the first to the
                  ! Uniquer_Num_Cr_Seg element.
                  call Vector_Unique_Int(5,tem_i_1,Crack_Seg_Num,
     &                           Uniqued_Cr_Seg_Num,Uniqued_Num_Cr_Seg) 
                  point_seg_1=Crack_Coor(i_C,Uniqued_Cr_Seg_Num(1),:)
            point_seg_2=Crack_Coor(i_C,Uniqued_Cr_Seg_Num(1)+1,:)
            AB(1,:) = point_seg_1
                  AB(2,:) = point_seg_2
            CD(1,:) = [crack_inter(1,1),crack_inter(1,2)]
                  CD(2,:) = [crack_inter(2,1),crack_inter(2,2)]
                  tem_02  = AB(2,:)-AB(1,:)
                  tem_03  = CD(2,:)-CD(1,:)
                  if(dot_product(tem_02,tem_03) .GE. ZR) then
                      do iii = 1,num_Crack
                          Coors_Element_Crack(i_E,iii,1:4)=
     &                         [crack_inter(1,1),crack_inter(1,2),
     &                          crack_inter(2,1),crack_inter(2,2)]
                      enddo
                  else
                  !Change the order of "crack_inter" to make sure that the order are the same as the crack.
                      do iii = 1,num_Crack
                          Coors_Element_Crack(i_E,iii,1:4)=
     &                         [crack_inter(2,1),crack_inter(2,2),
     &                          crack_inter(1,1),crack_inter(1,2)]  
                      enddo
                  end if
              end if
              !
              !################################################
              !(Type 3):Fully cracked element with kink point.    
              !################################################
              if((Num_Intersect.eq.2).and.((Flag_1).or.(Flag_2)))then       
                  Elem_Type(i_E,i_C) = 3  
                  ! Check if the order of the two points of crack_inter need to be changed.
            ! Make sure that the order of "crack_inter" are the same as the crack.
            ! Loop to find how many crack segments have "crack_inter" points, actually 2.
                  tem_i_1 = 0
                  Crack_Seg_Num(1:5) =0
                  do i_inter=1,2
                      x_inter = crack_inter(i_inter,1)
                      y_inter = crack_inter(i_inter,2)
                      do i_Segment = 1,Each_Cr_Poi_Num(i_C)-1
                          point_seg_1 = Crack_Coor(i_C,i_Segment,:)
                          point_seg_2 = Crack_Coor(i_C,i_Segment+1,:)
                          ! Is the current crack segment an arc segment (2017-07-16)
                          c_Yes_Arc_Crack = .False.
                          if(sum(abs(Arc_Crack_Coor(
     &                                 i_C,i_Segment,1:11)))> ZR)then
                              c_Yes_Arc_Crack = .True.
                          endif
                          if(c_Yes_Arc_Crack .EQV. .False.) then
                              call Tool_Yes_On_Line(x_inter,y_inter,
     &                          point_seg_1,point_seg_2,Yes_on_Segment)
                          elseif(c_Yes_Arc_Crack)then
                              call Tool_Yes_On_Arc(x_inter,y_inter,
     &                             Arc_Crack_Coor(i_C,i_Segment,1:11),
     &                             Yes_on_Segment)
                          endif
      	      if (Yes_on_Segment) then
                                tem_i_1 = tem_i_1 + 1
                                Crack_Seg_Num(tem_i_1)  = i_Segment
      	      end if
                      end do                     
                  end do
                  ! Remove duplicate elements from the first several elements in the array and store them in
                  ! Uniqued_Cr_Seg_Num,
                  ! Among them, the valid data of the newly generated array is from the first to the
                  ! Uniquer_Num_Cr_Seg element.
                  call Vector_Unique_Int(5,tem_i_1,Crack_Seg_Num,
     &                           Uniqued_Cr_Seg_Num,Uniqued_Num_Cr_Seg) 
                  point_seg_1=Crack_Coor(i_C,Uniqued_Cr_Seg_Num(1),:)
            point_seg_2=Crack_Coor(i_C,Uniqued_Cr_Seg_Num(2)+1,:)
            AB(1,1:2) = point_seg_1
                  AB(2,1:2) = point_seg_2
            CD(1,1:2) = [crack_inter(1,1),crack_inter(1,2)]
                  CD(2,1:2) = [crack_inter(2,1),crack_inter(2,2)]
                  tem_02  = AB(2,1:2)-AB(1,1:2)
                  tem_03  = CD(2,1:2)-CD(1,1:2)
                  if(dot_product(tem_02,tem_03) .GE. ZR) then
                      do iii = 1,num_Crack
                          Coors_Element_Crack(i_E,iii,1:4)=
     &                          [crack_inter(1,1),crack_inter(1,2),
     &                           crack_inter(2,1),crack_inter(2,2)]
                      enddo
                  else
                  !Change the order of "crack_inter" to make sure that the order are the same as the crack.
                      do iii = 1,num_Crack
                          Coors_Element_Crack(i_E,iii,1:4)=
     &                         [crack_inter(2,1),crack_inter(2,2),
     &                          crack_inter(1,1),crack_inter(1,2)]  
                      enddo
                  end if   
                  Coors_Vertex(i_E,1:2)  = point_seg_inElemnt(i_E,1:2)
              end if
              !###################################################
              !(Type 1):Tip element or (Type 4):Junction element.    
              !###################################################
              if(Num_Intersect .eq. 1)then 
                  !The coordinates of the crack tip.
                  Coors_Tip(i_E,:) = point_seg_inElemnt(i_E,:)
                  ! If tip-splitting enhancement is performed
                  if(Key_TipEnrich /= 0 )then
                      Elem_Type(i_E,i_C) = 1
                      !Crack tip of tip enriched nodes.	
                      x_cr_tip_nodes(i_C,N1) = Coors_Tip(i_E,1)
                      x_cr_tip_nodes(i_C,N2) = Coors_Tip(i_E,1)
                      x_cr_tip_nodes(i_C,N3) = Coors_Tip(i_E,1)
                      x_cr_tip_nodes(i_C,N4) = Coors_Tip(i_E,1)
                      y_cr_tip_nodes(i_C,N1) = Coors_Tip(i_E,2)
                      y_cr_tip_nodes(i_C,N2) = Coors_Tip(i_E,2)
                      y_cr_tip_nodes(i_C,N3) = Coors_Tip(i_E,2)
                      y_cr_tip_nodes(i_C,N4) = Coors_Tip(i_E,2)
                  end if
              
                  ! Check Junction Enhanced Node
                  do j_C=1,num_Crack
                      if(j_C.ne.i_C) then
                          !Loop through each segment of j_C. 
                          do j_S = 1,Each_Cr_Poi_Num(j_C)-1
                              !Get the coordinate of the two points of the current segment.
                              tem_point_seg_1 = Crack_Coor(j_C,j_S,:)
                              tem_point_seg_2 = Crack_Coor(j_C,j_S+1,:)
                              !Check that if the crack tip is on the segment or not.
                              call Tool_Yes_On_Line
     &                            (Coors_Tip(i_E,1),Coors_Tip(i_E,2),
     &                             tem_point_seg_1,tem_point_seg_2,
     &                             Yes_On_Line)
                              !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                              ! If yes, then: Junction enhancement element
                              !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      		      if(Yes_On_Line)then
      		          Elem_Type(i_E,i_C) = 4
      			      Coors_Junction(i_E,i_C,1:2)  = 
     &                              crack_inter(1,:)
                                  Coors_Junction(i_E,i_C,3:4)  = 
     &                              [Coors_Tip(i_E,1),Coors_Tip(i_E,2)]
                                  Jun_Ele_Negative_Cr_Num(i_E,i_C)=j_C
      			      !Determine the type of the crack tip.
                                  if (point_seg_inElemnt_Flag(i_E).eq.1)
     &                                  then
                                        Crack_Tip_Type(i_C,1) = 1
                                        Crack_Jun_CrNum(i_C,1)= j_C
                                        Crack_Jun_Elem(i_C,1) = i_E
                                  elseif (
     &                                point_seg_inElemnt_Flag(i_E).eq.2)
     &                                  then
                                        Crack_Tip_Type(i_C,2) = 1
                                        Crack_Jun_CrNum(i_C,2)= j_C
                                        Crack_Jun_Elem(i_C,2) = i_E
                                  end if
      		      end if
                          end do
                      end if
                  end do
                  !&&&&&&&&&&&&&&&
            ! Tip splitting unit.
                  !&&&&&&&&&&&&&&&
            if ((Elem_Type(i_E,i_C) .ne. 4) .and. 
     &                (Key_TipEnrich .ne. 0)) then
                      Elem_Type(i_E,i_C) = 1
                      do iii = 1,num_Crack
                          Coors_Element_Crack(i_E,iii,1:2)=
     &                              crack_inter(1,1:2) 
                      enddo
                      cc_tip_x = Coors_Tip(i_E,1);
                      cc_tip_y = Coors_Tip(i_E,2);
                      Coors_Element_Crack(i_E,i_C,3:4)= 
     &                               [cc_tip_x,cc_tip_y]
            !If there is edge crak,then the crack tip should be changed to an edge.
                      if (i_E .eq. c_Edge_Elem) then
                          Elem_Type(c_Edge_Elem,i_C) = 2
                      end if
                      
                  end if
              end if
          end do loop_Ele
          
          
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! If there are arc-shaped cracks, special treatment is given to the crack tips of the arc-shaped
          ! cracks, 2017-07-18
          ! Calculate Crack_Arc_Tip_A_B_C data for arc-shaped cracks, 2017-07-18
          !   --------------------------------------
          ! Used for enhancing crack tip correction, specifically see page 19 of my doctoral thesis. Because
          ! the curved crack is considered to be composed of several line segments, therefore
          ! The original crack tip enhancement correction needs to be adjusted because the number of crack
          ! segments has essentially increased; the approach here is to calculate directly.
          ! Coordinates of points A, B, and C
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if (Yes_Arc_Crack)then 
              num_points_crack = Each_Cr_Poi_Num(i_C)
              ! Current cracks are cycling between each crack segment
              do i_P = 1,num_points_crack
                  crack_p1=[Crack_Coor(i_C,i_P,1),Crack_Coor(i_C,i_P,2)]
                  ! The current crack tip is crack tip number 1.
                  if(i_P==1) then
                      ! Coordinates of point A
                      Crack_Arc_Tip_A_B_C_x(i_C,1,1) =crack_p1(1)
                      Crack_Arc_Tip_A_B_C_y(i_C,1,1) =crack_p1(2)
                      ! element number of the crack tip
                      call Cal_Ele_Num_by_Coors(crack_p1(1),crack_p1(2),
     &                                          c_OUT_Elem)
                      crack_inter(1,1) = 
     &                            Coors_Element_Crack(c_OUT_Elem,i_C,1)
                      crack_inter(1,2) = 
     &                            Coors_Element_Crack(c_OUT_Elem,i_C,2)
                      ! Coordinates of point B
                      Crack_Arc_Tip_A_B_C_x(i_C,1,2) =crack_inter(1,1)
                      Crack_Arc_Tip_A_B_C_y(i_C,1,2) =crack_inter(1,2)  
                      !||||||||||||||||||||||||||||||||||
                      ! Next, the coordinates of point C
                      !||||||||||||||||||||||||||||||||||
                      ! Determine whether the current crack segment is curved; if not, there is no need to calculate it.
                      if(sum(abs(Arc_Crack_Coor(i_C,1,1:11)))>Tol_11)
     &                                                           then
                        ! First, find the cell number of the adjacent cell through which the curved crack passes.
                        do i_El_Check=1,num_Ele_Eles(c_OUT_Elem)
                          c_Check_El=Ele_Elements(c_OUT_Elem,i_El_Check)
                          ! It cannot be the current cell itself
                          if(c_Check_El/=c_OUT_Elem)then
                          ! If the current detection unit contains coordinate points that have been pierced by the arc, and
                          ! both points are on the arc, then this
                          ! A unit is the target element
                          c_Point_1 = 
     &                       [Coors_Element_Crack(c_Check_El,i_C,1),
     &                        Coors_Element_Crack(c_Check_El,i_C,2)]
                          c_Point_2 = 
     &                       [Coors_Element_Crack(c_Check_El,i_C,3),
     &                        Coors_Element_Crack(c_Check_El,i_C,4)]
                          call Tool_Yes_On_Arc(
     &                        c_Point_1(1),c_Point_1(2),
     &                        Arc_Crack_Coor(i_C,1,1:11),Yes_P1_ON)
                          call Tool_Yes_On_Arc(
     &                        c_Point_2(1),c_Point_2(2),
     &                        Arc_Crack_Coor(i_C,1,1:11),Yes_P2_ON)
                          if(Yes_P1_ON .and. Yes_P2_ON)then 
                              ! Determine point C, which is definitely either _Point_1 or _Point_2, and not the one that coincides
                              ! with point B.
                              call Tool_Yes_Point_Equal(
     &                               crack_inter(1,1:2),c_Point_1,
     &                               Yes_P_Equal_1)
                              call Tool_Yes_Point_Equal(
     &                               crack_inter(1,1:2),c_Point_2,
     &                               Yes_P_Equal_2)
                              if(Yes_P_Equal_1)then
                                  ! Coordinates of point C
                                   Crack_Arc_Tip_A_B_C_x(i_C,1,3) =
     &                                        c_Point_2(1)
                                   Crack_Arc_Tip_A_B_C_y(i_C,1,3) =
     &                                        c_Point_2(2) 
                              endif
                              if(Yes_P_Equal_2)then
                                  ! Coordinates of point C
                                   Crack_Arc_Tip_A_B_C_x(i_C,1,3) =
     &                                        c_Point_1(1)
                                   Crack_Arc_Tip_A_B_C_y(i_C,1,3) =
     &                                        c_Point_1(2) 
                              endif
                              exit
                          endif
                          endif
                        enddo
                      endif
                  ! Current crack tip is crack tip No. 2
                  elseif(i_P==num_points_crack) then
                      ! Coordinates of point A
                      Crack_Arc_Tip_A_B_C_x(i_C,2,1) =crack_p1(1)
                      Crack_Arc_Tip_A_B_C_y(i_C,2,1) =crack_p1(2)
                      ! element number of the crack tip
                      call Cal_Ele_Num_by_Coors(crack_p1(1),crack_p1(2),
     &                                          c_OUT_Elem)
                      crack_inter(1,1) = 
     &                            Coors_Element_Crack(c_OUT_Elem,i_C,1)
                      crack_inter(1,2) = 
     &                            Coors_Element_Crack(c_OUT_Elem,i_C,2)
                      ! Coordinates of point B
                      Crack_Arc_Tip_A_B_C_x(i_C,2,2) =
     &                                        crack_inter(1,1)
                      Crack_Arc_Tip_A_B_C_y(i_C,2,2) =
     &                                        crack_inter(1,2)  
                      !||||||||||||||||||||||||||||||||||
                      ! Next, the coordinates of point C
                      !||||||||||||||||||||||||||||||||||
                      ! Determine whether the current crack segment is curved; if not, there is no need to calculate it.
                      if(sum(abs(Arc_Crack_Coor(i_C,
     &                            num_points_crack-1,1:11)))>Tol_11)then
                        ! First, find the cell number of the adjacent cell through which the curved crack passes.
                        do i_El_Check=1,num_Ele_Eles(c_OUT_Elem)
                          c_Check_El=Ele_Elements(c_OUT_Elem,i_El_Check)
                          ! It cannot be the current cell itself
                          if(c_Check_El/=c_OUT_Elem)then
                          ! If the current detection unit has coordinates that have been pierced by the arc, and both points
                          ! are on the arc, then this
                          ! A unit is the target element
                          c_Point_1 = 
     &                       [Coors_Element_Crack(c_Check_El,i_C,1),
     &                        Coors_Element_Crack(c_Check_El,i_C,2)]
                          c_Point_2 = 
     &                       [Coors_Element_Crack(c_Check_El,i_C,3),
     &                        Coors_Element_Crack(c_Check_El,i_C,4)]
                          call Tool_Yes_On_Arc(
     &                        c_Point_1(1),c_Point_1(2),
     &                        Arc_Crack_Coor(i_C,num_points_crack-1,
     &                                       1:11),Yes_P1_ON)
                          call Tool_Yes_On_Arc(
     &                        c_Point_2(1),c_Point_2(2),
     &                        Arc_Crack_Coor(i_C,num_points_crack-1,
     &                                       1:11),Yes_P2_ON)
                          if(Yes_P1_ON .and. Yes_P2_ON)then 
                              ! Determine point C, which is definitely either point 1 or point 2, and is not the same as point B.
                              call Tool_Yes_Point_Equal(
     &                             crack_inter(1,1:2),
     &                             c_Point_1,Yes_P_Equal_1)
                              call Tool_Yes_Point_Equal(
     &                             crack_inter(1,1:2),
     &                             c_Point_2,Yes_P_Equal_2)
                              if(Yes_P_Equal_1)then
                                  ! Coordinates of point C
                                   Crack_Arc_Tip_A_B_C_x(i_C,2,3) =
     &                                        c_Point_2(1)
                                   Crack_Arc_Tip_A_B_C_y(i_C,2,3) =
     &                                        c_Point_2(2) 
                              endif
                              if(Yes_P_Equal_2)then
                                  ! Coordinates of point C
                                   Crack_Arc_Tip_A_B_C_x(i_C,2,3) =
     &                                        c_Point_1(1)
                                   Crack_Arc_Tip_A_B_C_y(i_C,2,3) =
     &                                        c_Point_1(2) 
                              endif
                              exit
                          endif
                          endif
                        enddo
                      endif
                   end if
               enddo
          endif


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Determine the intersection between the crack and the circular hole, and locate the
          ! enhancement
          ! type of the element at the intersection to 6.
          !2017-05-19
          ! Elliptical holes and circular holes cannot coexist.
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do i_H=1,num_Circ_Hole
              x0_Hole =  Hole_Coor(i_H,1)
              y0_Hole =  Hole_Coor(i_H,2)
              R0_Hole =  Hole_Coor(i_H,3) 
              do i_Tip = 1,2
                  Coor_Tip(1:2) = Crack_Tip_Coor(i_C,i_Tip,1:2)
                  ! Check whether the crack tip is at the boundary of the hole
                  c_Dis = Tool_Function_2Point_Dis(Coor_Tip,
     &                                            [x0_Hole,y0_Hole])
                  if(abs(c_Dis - R0_Hole) <= Tol_10) then
                      Crack_Tip_Type(i_C,i_Tip) = 1
                      Crack_Jun_HoleNum(i_C,i_Tip)= i_H
                      call Cal_Ele_Num_by_Coors(Coor_Tip(1),Coor_Tip(2),
     &                                          c_Jun_Ele)
                      Elem_Type(c_Jun_Ele,i_C)    = 6
                      Ele_Jun_Hole(c_Jun_Ele,i_C) = i_H
                      !..........................
                      ! Calculate Coors_Junction
                      !..........................
                      if(i_Tip==1)then
                          crack_p1 = [Crack_Coor(i_C,1,1),
     &                                Crack_Coor(i_C,1,2)]
                          crack_p2 = [Crack_Coor(i_C,2,1),
     &                                Crack_Coor(i_C,2,2)]
                      elseif(i_Tip==2)then
                          crack_p1 = 
     &                       [Crack_Coor(i_C,Each_Cr_Poi_Num(i_C),1),
     &                        Crack_Coor(i_C,Each_Cr_Poi_Num(i_C),2)]
                          crack_p2 = 
     &                       [Crack_Coor(i_C,Each_Cr_Poi_Num(i_C)-1,1),
     &                        Crack_Coor(i_C,Each_Cr_Poi_Num(i_C)-1,2)]
                      endif
                      ! Loop through the four edges of the current element
                      N1  = Elem_Node(c_Jun_Ele,1)                                         
                      N2  = Elem_Node(c_Jun_Ele,2)                                             
                      N3  = Elem_Node(c_Jun_Ele,3)                                             
                      N4  = Elem_Node(c_Jun_Ele,4)                                            
                      NN_L= [N1,N2,N3,N4,N1]
                      do i_Side = 1,4
                          n_node_1 = NN_L(i_Side)
                          n_node_2 = NN_L(i_Side+1)
                          !Get the coordinate of the two points of the current side.
                          point_side_1 =[Coor(n_node_1,1),
     &                                   Coor(n_node_1,2)]
                          point_side_2 =[Coor(n_node_2,1),
     &                                   Coor(n_node_2,2)]   
                          !Check if they intersects and calculate intersection point.
                          call Tool_Intersection(
     &                                       crack_p1,crack_p2,
     &                                       point_side_1,point_side_2,
     &                                       c_X,c_Y,c_Yes_Cross) 
                          if(c_Yes_Cross.eqv..True.)then
                              Coors_Junction(c_Jun_Ele,i_C,1:2) =
     &                                               [c_X,c_Y]
                              Coors_Junction(c_Jun_Ele,i_C,3:4) =
     &                                                Coor_Tip
                              exit
                          endif
                      enddo
                      !...............................
                      ! Calculate Coors_Element_Crack
                      !...............................
                      c_num_Inter =0
                      do i_Side = 1,4
                          n_node_1 = NN_L(i_Side)
                          n_node_2 = NN_L(i_Side+1)
                          !Get the coordinate of the two points of the current side.
                          point_side_1 =[Coor(n_node_1,1),
     &                                   Coor(n_node_1,2)]
                          point_side_2 =[Coor(n_node_2,1),
     &                                   Coor(n_node_2,2)]            
                          !Check if they intersects and calculate intersection point.
                          call Tool_Intersection_Line_and_Circle(
     &                                    x0_Hole,y0_Hole,R0_Hole,
     &                                    point_side_1,point_side_2,
     &                                    c_num_Inter_CH,c_State,
     &                                    c_Inter)
                          if (c_State==2 .and. c_num_Inter_CH ==1) then
                              c_num_Inter = c_num_Inter +1
                              Coors_Element_Crack(
     &                              c_Jun_Ele,i_C,c_num_Inter*2-1) = 
     &                              c_Inter(1,1)
                              Coors_Element_Crack(
     &                              c_Jun_Ele,i_C,c_num_Inter*2)   =
     &                              c_Inter(1,2)
                          endif
                      enddo
                  endif
              enddo
          enddo
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Determine the intersection between the crack and the elliptical hole, and assign element
          ! enhancement type 6 to the elements at the intersection
          !2020-08-09
          ! Elliptical holes and circular holes cannot coexist.
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do i_H=1,num_Ellip_Hole
              x0_Hole =  Ellip_Hole_Coor(i_H,1)
              y0_Hole =  Ellip_Hole_Coor(i_H,2)
              a_Hole  =  Ellip_Hole_Coor(i_H,3) 
              b_Hole  =  Ellip_Hole_Coor(i_H,4) 
              theta_Hole =  Ellip_Hole_Coor(i_H,5)
              do i_Tip = 1,2
                  !Coodinates of fracture tip.
                  Coor_Tip(1:2) = Crack_Tip_Coor(i_C,i_Tip,1:2)
                  ! Check whether the crack tip is on the boundary of the hole
                  call Tool_Yes_Point_in_Oblique_Ellipse(Coor_Tip,
     &                  x0_Hole,y0_Hole,a_Hole,b_Hole,theta_Hole,
     &                  Return_Statu,Tol_10)
                  if(Return_Statu==0) then
                      Crack_Tip_Type(i_C,i_Tip) = 1
                      Crack_Jun_HoleNum(i_C,i_Tip)= i_H
                      call Cal_Ele_Num_by_Coors(Coor_Tip(1),Coor_Tip(2),
     &                                          c_Jun_Ele)
                      Elem_Type(c_Jun_Ele,i_C)    = 6
                      Ele_Jun_Hole(c_Jun_Ele,i_C) = i_H
                      !..........................
                      ! Calculate Coors_Junction
                      !..........................
                      if(i_Tip==1)then
                          crack_p1 = [Crack_Coor(i_C,1,1),
     &                                Crack_Coor(i_C,1,2)]
                          crack_p2 = [Crack_Coor(i_C,2,1),
     &                                Crack_Coor(i_C,2,2)]
                      elseif(i_Tip==2)then
                          crack_p1 = 
     &                       [Crack_Coor(i_C,Each_Cr_Poi_Num(i_C),1),
     &                        Crack_Coor(i_C,Each_Cr_Poi_Num(i_C),2)]
                          crack_p2 = 
     &                       [Crack_Coor(i_C,Each_Cr_Poi_Num(i_C)-1,1),
     &                        Crack_Coor(i_C,Each_Cr_Poi_Num(i_C)-1,2)]
                      endif
                      ! Loop through the four edges of the current element
                      N1  = Elem_Node(c_Jun_Ele,1)                                         
                      N2  = Elem_Node(c_Jun_Ele,2)                                             
                      N3  = Elem_Node(c_Jun_Ele,3)                                             
                      N4  = Elem_Node(c_Jun_Ele,4)                                            
                      NN_L= [N1,N2,N3,N4,N1]
                      do i_Side = 1,4
                          n_node_1 = NN_L(i_Side)
                          n_node_2 = NN_L(i_Side+1)
                          !Get the coordinate of the two points of the current side.
                          point_side_1 =[Coor(n_node_1,1),
     &                                   Coor(n_node_1,2)]
                          point_side_2 =[Coor(n_node_2,1),
     &                                   Coor(n_node_2,2)]   
                          !Check if they intersects and calculate intersection point.
                          call Tool_Intersection(
     &                                       crack_p1,crack_p2,
     &                                       point_side_1,point_side_2,
     &                                       c_X,c_Y,c_Yes_Cross) 
                          if(c_Yes_Cross.eqv..True.)then
                              Coors_Junction(c_Jun_Ele,i_C,1:2) =
     &                                               [c_X,c_Y]
                              Coors_Junction(c_Jun_Ele,i_C,3:4) =
     &                                                Coor_Tip
                              exit
                          endif
                      enddo
                      !...............................
                      ! Calculate Coors_Element_Crack
                      !...............................
                      c_num_Inter =0
                      do i_Side = 1,4
                          n_node_1 = NN_L(i_Side)
                          n_node_2 = NN_L(i_Side+1)
                          !Get the coordinate of the two points of the current side.
                          point_side_1 =[Coor(n_node_1,1),
     &                                   Coor(n_node_1,2)]
                          point_side_2 =[Coor(n_node_2,1),
     &                                   Coor(n_node_2,2)]            
                          !Check if they intersects and calculate intersection point.
                          call Tool_Intersection_Line_and_Circle(
     &                                    x0_Hole,y0_Hole,R0_Hole,
     &                                    point_side_1,point_side_2,
     &                                    c_num_Inter_CH,c_State,
     &                                    c_Inter)
                          if (c_State==2 .and. c_num_Inter_CH ==1) then
                              c_num_Inter = c_num_Inter +1
                              Coors_Element_Crack(
     &                              c_Jun_Ele,i_C,c_num_Inter*2-1) = 
     &                              c_Inter(1,1)
                              Coors_Element_Crack(
     &                              c_Jun_Ele,i_C,c_num_Inter*2)   =
     &                              c_Inter(1,2)
                          endif
                      enddo
                  endif
              enddo
          enddo
      end do loop_Crack

                   
c     -----------------------------------------------------------------------------------
c     ---------------------   Step 2: determine the enriched nodes  --------------------
c     -----------------------------------------------------------------------------------
      loop 1: do i_C=1,num_Crack
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Cycle between each element
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          loop 2: do i_E = 1,Num_Elem
              N1  = Elem_Node(i_E,1)                              
              N2  = Elem_Node(i_E,2)                          
              N3  = Elem_Node(i_E,3)                          
              N4  = Elem_Node(i_E,4)  
              NN  = [N1,N2,N3,N4]
              NN_L= [N1,N2,N3,N4,N1]
              X_NODES_L = Coor(NN_L,1)
              Y_NODES_L = Coor(NN_L,2)
              !*******************************************
              !If type-1 element, i.e. tip element, then:
              !*******************************************
            if (Elem_Type(i_E,i_C) .eq. 1) then
                  if(Key_TipEnrich /= 0 )then
                    Enriched_Node_Type(NN,i_C) =1
                      Ele_Num_Tip_Enriched_Node(i_C,NN)=i_E
                  endif
              !************************************************************************
              !If type-2 element, i.e. fully cracked element without kink point, then:
              !************************************************************************
            elseif (Elem_Type(i_E,i_C) .eq. 2) then
                  !Loop through each node of the current element.
                  do i_N = 1,4
                !If the node has not been enriched yet, then:
                      if (Enriched_Node_Type(NN(i_N),i_C).eq. 0) then
                !Get the support domain of the current node.
      	      call Cal_Support_Domain_of_Node(NN(i_N),
     &                                  DOMAIN_Outline,
     &                                  m_DOMAIN_Outline,
     &                                  Domain_El,
     &                                  n_Domain_El)
      	      !Get the A+ and A- Domain.
                          call Cal_Ap_and_Am(Crack_Coor(i_C,:,:),
     &                         Each_Cr_Poi_Num(i_C),
     &                         DOMAIN_Outline(1:m_DOMAIN_Outline,:),
     &                         m_DOMAIN_Outline,
     &                         Domain_El(1:n_Domain_El),n_Domain_El,
     &                         i_E,NN(i_N),
     &                         Yes_Cutthrough,
     &                         Yes_Ap_has_GP,Yes_Am_has_GP)

      	      !------------------------------------------------------------						 
      	      !Option 1: If both A_Plus and A_Minus has gauss points, then:
      	      !------------------------------------------------------------
      	      if (Yes_Ap_has_GP .and. Yes_Am_has_GP) then
      	          Enriched_Node_Type(NN(i_N),i_C) = 2
      	      end if
                      end if
            end do  
              !*********************************************************************    
              !If type-3 element, i.e. Fully cracked element with kink point, then:
              !*********************************************************************
            elseif (Elem_Type(i_E,i_C) .eq. 3) then
                  !Loop through each node of the current element.
                  do i_N = 1,4
                !If the node has not been enriched yet, then:
                      if (Enriched_Node_Type(NN(i_N),i_C).eq. 0) then
                !Get the support domain of the current node.
      	      call Cal_Support_Domain_of_Node(NN(i_N),
     &                                  DOMAIN_Outline,
     &                                  m_DOMAIN_Outline,
     &                                  Domain_El,
     &                                  n_Domain_El)
      	      !Get the A+ and A- Domain.
                          call Cal_Ap_and_Am(Crack_Coor(i_C,:,:),
     &                         Each_Cr_Poi_Num(i_C),
     &                         DOMAIN_Outline(1:m_DOMAIN_Outline,:),
     &                         m_DOMAIN_Outline,
     &                         Domain_El(1:n_Domain_El),n_Domain_El,
     &                         i_E,NN(i_N),
     &                         Yes_Cutthrough,
     &                         Yes_Ap_has_GP,Yes_Am_has_GP)
      	      !------------------------------------------------------------						 
      	      !Option 1: If both A_Plus and A_Minus has gauss points, then:
      	      !------------------------------------------------------------
      	      if (Yes_Ap_has_GP .and. Yes_Am_has_GP) then
      	          Enriched_Node_Type(NN(i_N),i_C) = 2
      	      end if
              end if
            end do
              !************************************************
              !If type-4 element, i.e. Junction element, then:
              !************************************************
            elseif (Elem_Type(i_E,i_C) .eq. 4) then
                  do i_N = 1,4
                      Enriched_Node_Type(NN(i_N),i_C) = 3
                      Node_Jun_elem(NN(i_N),i_C) = i_E
                  end do
                  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                  ! Determine additional Junction enhancement nodes to ensure that elements in the
                  ! support domain containing two cracks are all subjected to Junction enhancement,
                  ! added on 2016-07-10
                  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                  ! See the loop for the elements around the junction element
                  tem_Checked_Nodes(1:Num_Node) = 0
                  loop 3: do i_Su_Ele = 1,num_Ele_Eles(i_E)
                      c_Chk_Ele = Ele_Elements(i_E,i_Su_Ele)
                      ! Looping through each node of the currently inspected unit
                      do i_Chk_Node = 1,4
                          c_Chk_N  = Elem_Node(c_Chk_Ele,i_Chk_Node)    
                          ! If the element has already been inspected, it will not be inspected again.
                          if(tem_Checked_Nodes(c_Chk_N)==1)then
                              exit
                          endif
                !Get the support domain of the current node.
      	      call Cal_Support_Domain_of_Node(c_Chk_N,
     &                           DOMAIN_Outline,m_DOMAIN_Outline,
     &                           Domain_El,n_Domain_El)
                          ! Convert support domain boundaries (node numbers) into a closed polygon
                          Clo_Supp_Domain(1:20,1:2) = ZR
                          Clo_Supp_Domain(1:m_DOMAIN_Outline,1:2) = 
     &                       Coor(DOMAIN_Outline(1:m_DOMAIN_Outline,1),
     &                            1:2)     
                          Clo_Supp_Domain(m_DOMAIN_Outline+1,1:2) = 
     &                       Coor(DOMAIN_Outline(1,1),1:2)  
                          ! Obtain the primary and secondary fracture numbers
                          Positive_Cr = i_C
                          Negative_Cr = Jun_Ele_Negative_Cr_Num(i_E,i_C)
                          ! Check if active cracks intersect with the support area_Yes_Posit_Cross
                          call Tool_Yes_Crack_Polygon_Intersects(
     &                         Positive_Cr,
     &                         Clo_Supp_Domain(1:m_DOMAIN_Outline+1,1:2)
     &                         ,m_DOMAIN_Outline+1,Yes_Posit_Cross)
                          ! Check if the passive cracks intersect with the support zone_Yes_Negat_Cross
                          call Tool_Yes_Crack_Polygon_Intersects(
     &                         Negative_Cr,
     &                         Clo_Supp_Domain(1:m_DOMAIN_Outline+1,1:2)
     &                         ,m_DOMAIN_Outline+1,Yes_Negti_Cross)
                          ! If it intersects both active and passive cracks at the same time, a new Junction node will be
                          ! added.
                          if(Yes_Posit_Cross .and. Yes_Negti_Cross)then
                              Enriched_Node_Type(c_Chk_N,i_C) = 3
                              Node_Jun_elem(c_Chk_N,i_C) = i_E
                          endif
                          ! Mark this node as detected
                          tem_Checked_Nodes(c_Chk_N) = 1
                      enddo
                  enddo loop 3
              !******************************************************************
              !If type-6 element, i.e. Junction element of crack and hole, then:
              !******************************************************************
            elseif (Elem_Type(i_E,i_C) .eq. 6) then
                  do i_N = 1,4
                      Enriched_Node_Type(NN(i_N),i_C) = 6
                      Node_Jun_elem(NN(i_N),i_C) = i_E
                      Node_Jun_Hole(NN(i_N),i_C) = Ele_Jun_Hole(i_E,i_C)
                  end do
                  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                  ! Determine additional junction enhancement nodes to ensure the support domain
                  ! contains the current crack
                  ! Perform Junction (type 6) enhancement on all elements including Hole
                  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                  c_Hole = Ele_Jun_Hole(i_E,i_C)
                  x0_Hole =  Hole_Coor(c_Hole,1)
                  y0_Hole =  Hole_Coor(c_Hole,2)
                  R0_Hole =  Hole_Coor(c_Hole,3) 
                  ! Inter-unit circulation around the junction element
                  tem_Checked_Nodes(1:Num_Node) = 0
                  Loop_i_Su_Ele:do i_Su_Ele = 1,num_Ele_Eles(i_E)
                      c_Chk_Ele = Ele_Elements(i_E,i_Su_Ele)
                      ! Loop through each node of the currently inspected unit
                      Loop_i_Chk_Node:do i_Chk_Node = 1,4
                          c_Chk_N  = Elem_Node(c_Chk_Ele,i_Chk_Node)    
                          ! If the element has already been inspected, it will not be inspected again.
                          if(tem_Checked_Nodes(c_Chk_N)==1)then
                              exit
                          endif
                !Get the support domain of the current node.
      	      call Cal_Support_Domain_of_Node(c_Chk_N,
     &                           DOMAIN_Outline,m_DOMAIN_Outline,
     &                           Domain_El,n_Domain_El)
                          ! Convert support domain boundaries (node numbers) into a closed polygon
                          Clo_Supp_Domain(1:20,1:2) = ZR
                          Clo_Supp_Domain(1:m_DOMAIN_Outline,1:2) = 
     &                       Coor(DOMAIN_Outline(1:m_DOMAIN_Outline,1),
     &                            1:2)     
                          Clo_Supp_Domain(m_DOMAIN_Outline+1,1:2) = 
     &                       Coor(DOMAIN_Outline(1,1),1:2)  
                          ! Obtain main fracture number
                          Positive_Cr = i_C
                          ! Check if active cracks intersect with the support area_Yes_Posit_Cross
                          call Tool_Yes_Crack_Polygon_Intersects(
     &                         Positive_Cr,
     &                         Clo_Supp_Domain(1:m_DOMAIN_Outline+1,1:2)
     &                         ,m_DOMAIN_Outline+1,Yes_Posit_Cross)
                          ! Check if there is an overlap between holes and support regions
                          call Tool_Yes_Cricle_Polygon_Intersects(
     &                         x0_Hole,y0_Hole,R0_Hole,
     &                         Clo_Supp_Domain(1:m_DOMAIN_Outline+1,1:2)
     &                         ,m_DOMAIN_Outline+1,Yes_Hole_Cross)
                          ! If they intersect at the same time, a new Junction node (Crack and hole) is added.
                          if(Yes_Posit_Cross .and. Yes_Hole_Cross)then
                              Enriched_Node_Type(c_Chk_N,i_C) = 6
                              Node_Jun_elem(c_Chk_N,i_C) = i_E
                              Node_Jun_Hole(c_Chk_N,i_C) = 
     &                                        Ele_Jun_Hole(i_E,i_C)
                          endif
                          ! Mark this node as detected
                          tem_Checked_Nodes(c_Chk_N) = 1
                      enddo Loop_i_Chk_Node
                  enddo Loop_i_Su_Ele
              end if
          end do loop 2
      end do loop 1
      
      
c     -------------------------------------------------------------------------------------------
c     ----------------------------   Step 3: Check Enhanced Node Type  -------------------------
c     -------------------------------------------------------------------------------------------
      ! A node cannot be both junction-enhanced and tip-enhanced.
      do i_N=1,Num_Node
          do i_C=1,num_Crack
            if(Enriched_Node_Type(i_N,i_C)==3) then
              do j_C=1,num_Crack
                if(j_C /= i_C)then
                  if(Enriched_Node_Type(i_N,j_C)==1) then
                    print *, '    WARNING :: Junction enrichment and'
     &                //' tip enrichment should not coexist for a node!'
                    print *, '             Ill node is',i_N
                    print *, '             Related cracks are',i_C,j_C
                    !call Warning_Message('S',Keywords_Blank)
                  end if
                endif
              enddo
            end if
          enddo
      end do
      
c     ------------------------------------------------------------------------------------------
c     ------------------------ Step 4.1: Identify circular hole reinforcement elements and
c     reinforcement nodes -------------
c     ------------------------ Circular holes and elliptical holes do not coexist -------------
c     ------------------------------------------------------------------------------------------
      ! Execute only on the first rupture step
      if (num_Circ_Hole /= 0 .and. ifra==1) then   
          ! Cycle between each element
          do i_E = 1,Num_Elem
              N1  = Elem_Node(i_E,1)                                         
              N2  = Elem_Node(i_E,2)                                             
              N3  = Elem_Node(i_E,3)                                             
              N4  = Elem_Node(i_E,4)  
              NN  = [N1,N2,N3,N4]
              NN_L= [N1,N2,N3,N4,N1]
              X_NODES_L = Coor(NN_L,1)
              Y_NODES_L = Coor(NN_L,2)
              !**************************
              ! Circulation of each hole
              !**************************
              do i_H = 1,num_Circ_Hole
                  c_Hole_x = Hole_Coor(i_H,1)
                  c_Hole_y = Hole_Coor(i_H,2)
                  c_Hole_r = Hole_Coor(i_H,3)
                  ! Loop through the 4 nodes of the current element
                  Nodes_4_in(1:4) = .False.
                  do i_N = 1,4
                      c_Dis=Tool_Function_2Point_Dis(Coor(NN(i_N),1:2),
     &                                               Hole_Coor(i_H,1:2))
                      if (c_Dis<=c_Hole_r)then
                          Nodes_4_in(i_N) = .True.
                          
                      endif
                  enddo
                  num_in_Nodes = count(Nodes_4_in(1:4))
                  if(num_in_Nodes == 4)then
                      !............................................................................
                      ! Delete the following code on 2017-05-29, reason for deletion:
                      ! (1) Because the degrees of freedom of the nodes that are completely within
                      ! the hole have already been removed in Boundary_cond.f,
                      ! Therefore, no further enhancement is needed.
                      ! (2) Strengthening these nodes can cause the stiffness matrix to have an
                      ! excessively large condition number, making it unsolvable by some solvers.
                      ! Therefore, it cannot be enhanced.
                      !............................................................................
                      ! Elem_Type_Hl(i_E,i_H) = 1  !Then this element is a Hole-enhanced element
                      !Enriched_Node_Type_HL(NN(1:4),i_H) = 1
                      
                      !.......................
                      !do nothing, 2017-05-29
                      !.......................
                  elseif(num_in_Nodes <= 3 .and. num_in_Nodes >= 1)then
                      Elem_Type_Hl(i_E,i_H) = 1
                      do i_N = 1,4
                          c_Dis=Tool_Function_2Point_Dis(
     &                        Coor(NN(i_N),1:2),Hole_Coor(i_H,1:2))
                          if (c_Dis<=c_Hole_r)then
                              Enriched_Node_Type_HL(NN(i_N),i_H) = 1
                          else
                              Enriched_Node_Type_HL(NN(i_N),i_H) = 1
                          endif
                    enddo
                  endif
                  !Enriched_Node_Type_HL(1:Num_Node,1:Max_Num_Hl)
              enddo
          enddo
      endif
c     ------------------------------------------------------------------------------------------
c     ------------------------ Step 4.2: Identify Hole Reinforcement Elements and Reinforcement
c     Nodes
c     ------------------------------------------------------------------------------------------
c     ------------------------ Circular holes and elliptical holes do not coexist -------------
c     ------------------------             2020-08-09                      -------------
c     --------------------------------------------------------------------------------------------
      ! Execute only on the first rupture step
      if (num_Ellip_Hole /= 0 .and. ifra==1) then   
          ! Cycle between elements
          do i_E = 1,Num_Elem
              N1  = Elem_Node(i_E,1)                                         
              N2  = Elem_Node(i_E,2)                                             
              N3  = Elem_Node(i_E,3)                                             
              N4  = Elem_Node(i_E,4)  
              NN  = [N1,N2,N3,N4]
              NN_L= [N1,N2,N3,N4,N1]
              X_NODES_L = Coor(NN_L,1)
              Y_NODES_L = Coor(NN_L,2)
              !*************************
              ! Cycle through each hole
              !*************************
              do i_H = 1,num_Ellip_Hole
                  x0_Hole =  Ellip_Hole_Coor(i_H,1)
                  y0_Hole =  Ellip_Hole_Coor(i_H,2)
                  a_Hole  =  Ellip_Hole_Coor(i_H,3) 
                  b_Hole  =  Ellip_Hole_Coor(i_H,4) 
                  theta_Hole =  Ellip_Hole_Coor(i_H,5) 
                  ! Loop through the 4 nodes of the current element
                  Nodes_4_in(1:4) = .False.
                  do i_N = 1,4
                      ! Check whether the crack tip is inside the hole
                      call Tool_Yes_Point_in_Oblique_Ellipse(
     &                      Coor(NN(i_N),1:2),
     &                      x0_Hole,y0_Hole,a_Hole,b_Hole,theta_Hole,
     &                      Return_Statu,Tol_10)
                      if (Return_Statu==1)then
                          Nodes_4_in(i_N) = .True.
                          
                      endif
                  enddo
                  num_in_Nodes = count(Nodes_4_in(1:4))
                  if(num_in_Nodes == 4)then
                      !............................................................................
                      ! Delete the following code on 2017-05-29, reason for deletion:
                      ! (1) Because the degrees of freedom of the nodes that are completely within
                      ! the hole have already been removed in Boundary_cond.f,
                      ! Therefore, no further enhancement is needed.
                      ! (2) Strengthening these nodes can cause the stiffness matrix to have an
                      ! excessively large condition number, making it unsolvable by some solvers.
                      ! Therefore, it cannot be enhanced.
                      !............................................................................
                      ! Elem_Type_Hl(i_E,i_H) = 1  !Then this element is a Hole-enhanced element
                      !Enriched_Node_Type_HL(NN(1:4),i_H) = 1
                      
                      !.......................
                      !do nothing, 2017-05-29
                      !.......................
                      
                  elseif(num_in_Nodes <= 3 .and. num_in_Nodes >= 1)then
                      Elem_Type_Hl(i_E,i_H) = 1
                      do i_N = 1,4
                          call Tool_Yes_Point_in_Oblique_Ellipse(
     &                        Coor(NN(i_N),1:2),
     &                        x0_Hole,y0_Hole,a_Hole,b_Hole,theta_Hole,
     &                        Return_Statu,Tol_10)
                          if (Return_Statu==1)then
                              Enriched_Node_Type_HL(NN(i_N),i_H) = 1
                          else
                              Enriched_Node_Type_HL(NN(i_N),i_H) = 1
                          endif
                    enddo
                  endif
                  !Enriched_Node_Type_HL(1:Num_Node,1:Max_Num_Hl)
              enddo
          enddo
      endif      
      
c     -----------------------------------------------------------------------------------------
c     Step 5: Determine the reinforcement elements and reinforcement nodes at the cross-shaped
c     intersection
c     -----------------------------------------------------------------------------------------
      if (num_Cross /= 0) then   
          do i_Cross = 1,num_Cross
              !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
              ! 1. Identify the basic Cross enhancement nodes
              !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
              c_Cross_Ele = Cross_Point_Ele_num(i_Cross)
              
              Elem_Type_Cross(c_Cross_Ele,i_Cross) = 1
              
              N1  = Elem_Node(c_Cross_Ele,1)                                         
              N2  = Elem_Node(c_Cross_Ele,2)                                             
              N3  = Elem_Node(c_Cross_Ele,3)                                             
              N4  = Elem_Node(c_Cross_Ele,4)  
              Enriched_Node_Type_Cross(N1,i_Cross) = 1
              Enriched_Node_Type_Cross(N2,i_Cross) = 1
              Enriched_Node_Type_Cross(N3,i_Cross) = 1
              Enriched_Node_Type_Cross(N4,i_Cross) = 1
              Node_Cross_elem(Elem_Node(c_Cross_Ele,1:4) ,i_Cross) =
     &                                                   c_Cross_Ele
              !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
              ! 2. Identify additional Cross reinforcement nodes and ensure that all elements in the
              ! support domain containing two cracks are reinforced with Cross.
              !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
              ! Elements around the Cross unit see each other cyclically
              tem_Checked_Nodes(1:Num_Node) = 0
              do i_Su_Ele = 1,num_Ele_Eles(c_Cross_Ele)
                  c_Chk_Ele = Ele_Elements(c_Cross_Ele,i_Su_Ele)
                  ! Loop through each node of the currently inspected unit
                  do i_Chk_Node = 1,4
                      c_Chk_N  = Elem_Node(c_Chk_Ele,i_Chk_Node)    
                      ! If the element has already been inspected, it will not be inspected again.
                      if(tem_Checked_Nodes(c_Chk_N)==1)then
                          exit
                      endif
                      !Get the support domain of the current node.
                      call Cal_Support_Domain_of_Node(c_Chk_N,
     &                       DOMAIN_Outline,m_DOMAIN_Outline,
     &                       Domain_El,n_Domain_El)
                      ! Convert support domain boundaries (node numbers) into a closed polygon
                      Clo_Supp_Domain(1:20,1:2) = ZR
                      Clo_Supp_Domain(1:m_DOMAIN_Outline,1:2) = 
     &                   Coor(DOMAIN_Outline(1:m_DOMAIN_Outline,1),
     &                        1:2)     
                      Clo_Supp_Domain(m_DOMAIN_Outline+1,1:2) = 
     &                   Coor(DOMAIN_Outline(1,1),1:2)  
                      ! Obtain the primary and secondary fracture numbers
                      Positive_Cr=i_C
                      Negative_Cr=int(Cross_Point_Cr_num(i_Cross,2)) 
                      ! Check if active cracks intersect with the support area_Yes_Posit_Cross
                      call Tool_Yes_Crack_Polygon_Intersects(
     &                     Positive_Cr,
     &                     Clo_Supp_Domain(1:m_DOMAIN_Outline+1,1:2),
     &                     m_DOMAIN_Outline+1,Yes_Posit_Cross)
                      ! Check if the passive cracks intersect with the support zone_Yes_Negat_Cross
                      call Tool_Yes_Crack_Polygon_Intersects(
     &                     Negative_Cr,
     &                     Clo_Supp_Domain(1:m_DOMAIN_Outline+1,1:2),
     &                     m_DOMAIN_Outline+1,Yes_Negti_Cross)
                      ! If it intersects with both active and passive cracks at the same time, a new Cross node will be
                      ! added.
                      if(Yes_Posit_Cross .and. Yes_Negti_Cross)then
                          Enriched_Node_Type_Cross(c_Chk_N,i_Cross) = 1
                          Node_Cross_elem(c_NN(i_N),i_Cross)=c_Cross_Ele
                      endif
                      ! Mark this node as detected
                      tem_Checked_Nodes(c_Chk_N) = 1
                  enddo
              enddo          
          enddo
      endif
      
   
c     -----------------------------------------------------------------------------------------
c ---------------------------- Step 6: Determine the Circular Inclusion Reinforcement Unit and
c Reinforcement Node -----------------
c     -----------------------------------------------------------------------------------------
      ! Execute only on the first rupture step
      if (num_Circ_Incl /= 0 .and. ifra==1) then   
          ! Cycle between elements
          do i_E = 1,Num_Elem
              N1  = Elem_Node(i_E,1)                                         
              N2  = Elem_Node(i_E,2)                                             
              N3  = Elem_Node(i_E,3)                                             
              N4  = Elem_Node(i_E,4)  
              NN  = [N1,N2,N3,N4]
              NN_L= [N1,N2,N3,N4,N1]
              X_NODES_L = Coor(NN_L,1)
              Y_NODES_L = Coor(NN_L,2)
              !*************************
              ! Circular Inclusion Loop
              !*************************
              do i_Incl = 1,num_Circ_Incl
                  c_Incl_x = Circ_Inclu_Coor(i_Incl,1)
                  c_Incl_y = Circ_Inclu_Coor(i_Incl,2)
                  c_Incl_r = Circ_Inclu_Coor(i_Incl,3)
                  ! Loop through the 4 nodes of the current element
                  Nodes_4_in(1:4) = .False.
                  do i_N = 1,4
                      c_Dis =Tool_Function_2Point_Dis(Coor(NN(i_N),1:2),
     &                                     Circ_Inclu_Coor(i_Incl,1:2))
                      if (c_Dis<=c_Incl_r)then
                          Nodes_4_in(i_N) = .True.
                      endif
                  enddo
                  num_in_Nodes = count(Nodes_4_in(1:4))
                  if(num_in_Nodes == 4)then
                      !no nothing!
                  elseif(num_in_Nodes <= 3 .and. num_in_Nodes >= 1)then
                      Elem_Type_Incl(i_E,i_Incl) = 1
                      do i_N = 1,4
                          Enriched_Node_Type_Incl(NN(i_N),i_Incl)=1
                      enddo
                  endif
              enddo
          enddo
      endif
c     ---------------------------------------------------------------------------------
c     Step 7: Determine the polygon inclusion enhancement element and enhancement node
c     ---------------------------------------------------------------------------------
      ! Closed-form polygon inclusion
      do i_Incl = 1,num_Poly_Incl
          ! Number of vertices (sides) of the current polygon
          n_Vertex = Poly_Inclu_Edges_Num(i_Incl)
          ! Closed-form polygon inclusion
          Poly_Incl_Coor_x_Cl(i_Incl,1:n_Vertex)=
     &                      Poly_Incl_Coor_x(i_Incl,1:n_Vertex)
          Poly_Incl_Coor_x_Cl(i_Incl,n_Vertex+1)=
     &                      Poly_Incl_Coor_x(i_Incl,1)
          Poly_Incl_Coor_y_Cl(i_Incl,1:n_Vertex)=
     &                      Poly_Incl_Coor_y(i_Incl,1:n_Vertex)
          Poly_Incl_Coor_y_Cl(i_Incl,n_Vertex+1)=
     &                      Poly_Incl_Coor_y(i_Incl,1)
      enddo
      ! Execute only on the first rupture step
      if (num_Poly_Incl /= 0 .and. ifra==1) then   
          ! Cycle between each element
          do i_E = 1,Num_Elem
              N1  = Elem_Node(i_E,1)                                         
              N2  = Elem_Node(i_E,2)                                             
              N3  = Elem_Node(i_E,3)                                             
              N4  = Elem_Node(i_E,4)  
              NN  = [N1,N2,N3,N4]
              NN_L= [N1,N2,N3,N4,N1]
              X_NODES_L = Coor(NN_L,1)
              Y_NODES_L = Coor(NN_L,2)
              !******************************************
              ! Various polygons interspersed with loops
              !******************************************
              do i_Incl = 1,num_Poly_Incl
                  ! The number of vertices (sides) of the current polygon
                  n_Vertex = Poly_Inclu_Edges_Num(i_Incl)
                  ! Loop through the 4 nodes of the current element
                  Nodes_4_in(1:4) = .False.
                  do i_N = 1,4
                      call Tool_Yes_In_Poly(
     &                        Coor(NN(i_N),1),Coor(NN(i_N),2),
     &                        Poly_Incl_Coor_x_Cl(i_Incl,1:n_Vertex+1),
     &                        Poly_Incl_Coor_y_Cl(i_Incl,1:n_Vertex+1),
     &                        n_Vertex+1,
     &                        Nodes_4_in(i_N))
                  enddo
                  num_in_Nodes = count(Nodes_4_in(1:4))
                  if(num_in_Nodes == 4)then
                      !no nothing!
                  elseif(num_in_Nodes <= 3 .and. num_in_Nodes >= 1)then
                      Elem_Type_Incl(i_E,i_Incl) = 1
                      do i_N = 1,4
                          Enriched_Node_Type_Incl(NN(i_N),i_Incl)=1
                      enddo
                  endif
              enddo
          enddo
      endif

      


c     -------------------------------------------------------------
c     Step 9: Determine the average cell area of the enhanced unit   
c     -------------------------------------------------------------
      ! Cycle between elements
      Area_Enrich = ZR
      iEl_contour = 0
      do i_E = 1,Num_Elem
          if (maxval(Elem_Type(i_E,1:num_crack)) /=0 .or.
     &        minval(Elem_Type(i_E,1:num_crack)) /=0) then
              iEl_contour  = iEl_contour  + 1
              Area_Enrich = Area_Enrich + Elem_Area(i_E)
          end if
      end do
      Ave_Elem_Area_Enrich = Area_Enrich/iEl_contour
      Ave_Elem_L_Enrich    = sqrt(Ave_Elem_Area_Enrich)
      
      ! For the local encryption process, obtain the feature length of the enhanced unit before local
      ! encryption (2021-08-22)
      if (key_local_mesh_refine>=0  .and. 
     &   (Flag_Local_Refined.eqv..false.))then
          Ave_Elem_L_Enrich_Unlocalrefined = Ave_Elem_L_Enrich 
      endif
      
c     ---------------------------------------------------------------------------------
c      Step 10: If using the multi-tip enrichment node algorithm Key_Multi_TipEnrich=1
c     It is necessary to re-determine the crack tip enhancement elements and nodes            
c     ---------------------------------------------------------------------------------
      if (Key_TipEnrich /= 0 .and. Key_Multi_TipEnrNode==1 
     &                       .and. num_Crack>=1) then
        do i_C=1,num_Crack
          do i_Tip = 1,2
              ! Only detect ordinary crack tips; detection of non-ordinary crack tips is not required.
              if (Crack_Tip_Type(i_C,i_Tip) ==0) then
                  c_tip_x= Crack_Tip_Coor(i_C,i_Tip,1)
                  c_tip_y= Crack_Tip_Coor(i_C,i_Tip,2)
                  check_L=Key_TipEnrich_Radius_Factor*Ave_Elem_L_Enrich
                  check_x_min = c_tip_x-check_L
                  check_x_max = c_tip_x+check_L
                  check_y_min = c_tip_y-check_L
                  check_y_max = c_tip_y+check_L
                  do i_Node =1,Num_Node
                      ! Preliminary judgment
                      if(Coor(i_Node,1) < check_x_min) then
                          cycle
                      endif   
                      if(Coor(i_Node,1) > check_x_max) then
                          cycle
                      endif
                      if(Coor(i_Node,2) < check_y_min) then
                          cycle
                      endif   
                      if(Coor(i_Node,2) > check_y_max) then
                          cycle
                      endif
                      c_dis =Tool_Function_2Point_Dis([c_tip_x,c_tip_y],
     &                                                 Coor(i_Node,1:2))
                      ! If it falls within the detection radius
                      if(c_dis<=check_L)then
                          ! If this node has not been enhanced before (very important)
                          if(Enriched_Node_Type(i_Node,i_C)==0)then
                              ! Add Enhanced Node
                              Enriched_Node_Type(i_Node,i_C)=1
                              ! Crack tip enhancement node corresponding to the crack tip coordinates
                              x_cr_tip_nodes(i_C,i_Node) = c_tip_x
                              y_cr_tip_nodes(i_C,i_Node) = c_tip_y
                              ! element number of the crack tip
                              call Cal_Ele_Num_by_Coors(c_tip_x,c_tip_y,
     &                                                  c_ELe_Num)
                              ! The element number corresponding to the crack tip at the crack-tip enhancement node
                              Ele_Num_Tip_Enriched_Node(i_C,i_Node)=
     &                                                  c_ELe_Num
                          endif
                      endif
                  enddo
              endif
          enddo
        enddo
      endif
      
c     ----------------------------------------------------------------------
c      The adjacent elements of each crack tip (those crossed by the crack)
c      TipEle_Adjacent_Ele(num_crack,2)
c     ----------------------------------------------------------------------
      do i_E = 1,Num_Elem
          do i_C=1,num_Crack
              if (Elem_Type(i_E,i_C)==1) then
                  c_NN    = G_NN(:,i_E)
                  do i_N=1,4
                      c_N = c_NN(i_N)
                      ! Loop through the neighboring elements of the current node of the current element
                      do j_E =1,num_Node_Elements(c_N)
                          c_Adj_Ele = Node_Elements(c_N,j_E) 
                          ! If the potential element happens to be intersected by crack i_C, the search is successful.
                          if(Elem_Type(c_Adj_Ele,i_C)==2 .or.
     &                       Elem_Type(c_Adj_Ele,i_C)==3 .or.
     &                       Elem_Type(c_Adj_Ele,i_C)==4 .or.
     &                       Elem_Type(c_Adj_Ele,i_C)==5)then
                              TipEle_Adjacent_Ele(i_E,i_C) = c_Adj_Ele 
                              goto 222
                          endif
                      enddo
                  enddo
              end if
          enddo
  222 continue
      end do
      
      !***************************************************************************************
      ! In the full hydraulic fracturing model, the injected water fracture cannot be an edge
      ! fracture.
      !***************************************************************************************
      if (Key_Symm_HF==0 .and.(Key_Analysis_Type==3 .or.
     &                        Key_Analysis_Type==4 .or.
     &                        Key_Analysis_Type==5))then
          if((Crack_Tip_Type(Inject_Crack_Num,1)==-2) .or.
     &       (Crack_Tip_Type(Inject_Crack_Num,2)==-2))    then
              print *, '    Error :: injection crack cannot be an edge'
     &                               //' crack for the full HF model!'
              call Warning_Message('S',Keywords_Blank) 
          endif
      endif
      !********************************************************************************************
      ! The symmetrical hydraulic fracturing model requires water injection fractures to have edge
      ! fractures, and they must be at fracture tip 1.
      !********************************************************************************************
      if (Key_Symm_HF==1 .and.(Key_Analysis_Type==3 .or.
     &                        Key_Analysis_Type==4 .or.
     &                        Key_Analysis_Type==5))then
          if(Crack_Tip_Type(Inject_Crack_Num,1)/=-2)    then
              print *, '    Error :: injection crack must be an edge'
     &                           //' crack for the symmetric HF model!'
              call Warning_Message('S',Keywords_Blank) 
          endif
      endif
      

      
      RETURN
      END SUBROUTINE Determine_Enriched_Nodes
