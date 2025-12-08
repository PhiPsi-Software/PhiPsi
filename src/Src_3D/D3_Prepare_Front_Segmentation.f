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
 
      SUBROUTINE D3_Prepare_Front_Segmentation(isub)
c     Preparation for 3D crack front segmentation (finding suspended points).
C     Set the suppression segment for Crack Front Segmentation.
c     If it is a suppression segment, then Crack3D_Meshed_Outline(i_C, i, 4) is set to 0, and the first vertex of this segment can only extend by a minimal step.
c     Firstly Written on 2021-08-20.
c     Midofied on 2021-08-22.

C     real(kind=FT) Crack3D_Meshed_Node(Max_Num_Cr_3D, Max_N_Node_3D, 3)
C     integer Cr3D_Meshed_Node_in_Ele_Num(Max_Num_Cr_3D, Max_N_Node_3D)
C     real(kind=FT) Cr3D_Meshed_Node_in_Ele_Local(Max_Num_Cr_3D, Max_N_Node_3D, 3)
C     integer Crack3D_Meshed_Node_num(Max_Num_Cr_3D)
C     integer Crack3D_Meshed_Ele(Max_Num_Cr_3D, Max_N_Node_3D, 3)
c     integer Crack3D_Meshed_Ele_Attri(Max_Num_Cr_3D, Max_N_Node_3D, 5)
c     real(kind=FT)  Crack3D_Meshed_Ele_Nor_Vector(Max_Num_Cr_3D, Max_N_Node_3D, 3)
c     real(kind=FT)  Crack3D_Meshed_Node_Nor_Vector(Max_Num_Cr_3D,Max_N_Node_3D,3)
C     integer Crack3D_Meshed_Ele_num(Max_Num_Cr_3D)
C     integer Crack3D_Meshed_Outline(Max_Num_Cr_3D, Max_N_Node_3D, 4)
                                                                        ! Data 1 is the first point on the boundary line of the crack front edge
                                                                        ! Data 2 is the second point on the boundary line of the crack front edge
                                                                        ! Data 3 corresponds to the discrete fracture element number
                                                                        ! Data 4 is used to mark whether the two points of the boundary line are allowed to extend,
                                                                        ! extending in very small steps (2021-08-20)
C     integer Crack3D_Meshed_Outline_num(Max_Num_Cr_3D)
C     real(kind=FT) Crack3D_Meshed_Vertex_x_Vector(Max_Num_Cr_3D, Max_N_Node_3D,3)
C     real(kind=FT) Crack3D_Meshed_Vertex_y_Vector(Max_Num_Cr_3D, Max_N_Node_3D, 3)
C     real(kind=FT) Crack3D_Meshed_Vertex_z_Vector(Max_Num_Cr_3D, Max_N_Node_3D, 3)


c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Elem_Area_Vol
      use Global_Crack_Common
      use Global_Crack_3D
      
      
      
c     --------------------------
C     Variable Type Declaration
c     --------------------------
      implicit none
      integer,intent(in)::isub
      integer i_C,Num_CrMesh_Outlines,num_Seg_Outline,i_f_S
      integer i_Out_Node,c_Mesh_Node
      real(kind=FT) c_Point(3)
      integer c_Location
      logical c_Yes
      integer i,c_Out_Line,c_mod,max_mod_inc,mod_inc_count
      
      !=================
      !  First Substep.
      !=================
      if (isub==1) then
        do i_C =1,num_Crack  
          !/////////////////////////////////////////////////
          !                  OPTION-1
          !-------------------------------------------------
          !    Automatically Generates Suspended Points
          !    according to Number_front_Segments.
          !/////////////////////////////////////////////////
          Num_CrMesh_Outlines = Crack3D_Meshed_Outline_num(i_C)
          !Check Out_num_Outline.
          if (Num_CrMesh_Outlines < 5*Number_front_Segments)then
              print *, '    Error :: Num_CrMesh_Outlines <='
     &      //' 5*Number_front_Segments!'
              print *, '             Num_CrMesh_Outlines:',
     &                               Num_CrMesh_Outlines   
              print *, '             5*Number_front_Segments:',
     &                               5*Number_front_Segments       
              print *, '             Set *key_front_segmentation=0 or'
     &               //              'reduce *Number_front_Segments!'
              print *, '             Error in D3_Prepare_Front'
     &                 //            '_Segmentation.f.'  
              call Warning_Message('S',Keywords_Blank)
          endif
          ! The number of segments corresponding to each segment
          num_Seg_Outline=ceiling(dble(Num_CrMesh_Outlines)/
     &                            dble(Number_front_Segments))
          
          if (num_Seg_Outline*Number_front_Segments > 
     &        Num_CrMesh_Outlines) then
              num_Seg_Outline = num_Seg_Outline -1
          endif
          
          ! Calculate the remainder
          c_mod = mod(Num_CrMesh_Outlines,num_Seg_Outline)
          max_mod_inc = c_mod -1
          
          !Set suspended OutLines.
          mod_inc_count = 0
          c_Out_Line    = 0
          !c_Out_Line    = c_Out_Line  + 1
          do i_f_S=1,Number_front_Segments
              c_Out_Line = c_Out_Line + num_Seg_Outline
              ! Distribute the extra max_mod_inc outlines from the last segment to the other segments
              if (c_mod >=2 .and. mod_inc_count < max_mod_inc+1) then
                  mod_inc_count = mod_inc_count + 1
                  c_Out_Line = c_Out_Line + 1
              endif
              ! Check whether c_Out_Line is reasonable
              if (c_Out_Line > Num_CrMesh_Outlines)then
                  print *, '    Error :: illegal c_Out_Line!'
                  print *, '             Error in D3_Prepare_Front'
     &                     //            '_Segmentation.f.'  
                  call Warning_Message('S',Keywords_Blank)
              endif
              Crack3D_Meshed_Outline(i_C)%row(c_Out_Line,4) = 0 
              
          end do   
          !/////////////////////////////////////////////////
          !                  OPTION-2
          !-------------------------------------------------
          !    Manually provides Suspended Points.
          !    Crack3D_Meshed_Outline(i_C,i_Out_Node,4)=0 
          !/////////////////////////////////////////////////     
          !Crack3D_Meshed_Outline(i_C,1,4)=0 
        enddo  
        
      !==========================
      !  The following Substeps.
      !==========================
      elseif (isub>=2) then
        do i_C =1,num_Crack    
          Num_CrMesh_Outlines = Crack3D_Meshed_Outline_num(i_C)
          !Initialize to 1.
          Crack3D_Meshed_Outline(i_C)%row(1:Num_CrMesh_Outlines,4) =1
          do i_Out_Node = 1,Num_CrMesh_Outlines
              c_Mesh_Node =Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,1)
              c_Point(1)  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,1) 
              c_Point(2)  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,2) 
              c_Point(3)  = Crack3D_Meshed_Node(i_C)%row(c_Mesh_Node,3)
              ! Check whether the boundary point is in Suspended_Points(num_Suspended_Point, 1:3)
              call Vector_belongs_Matrix_Is_Dou(
     &                   num_Suspended_Point,3,
     &                   Suspended_Points(1:num_Suspended_Point,1:3),
     &                   c_Point(1:3),
     &                   c_Location,c_Yes)   
              !If yes then
              if(c_Yes)then
                  Crack3D_Meshed_Outline(i_C)%row(i_Out_Node,4)=0 
              endif
          enddo
        enddo  
      endif
      
      RETURN
      END SUBROUTINE D3_Prepare_Front_Segmentation
