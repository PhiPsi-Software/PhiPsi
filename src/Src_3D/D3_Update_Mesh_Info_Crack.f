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
 
      SUBROUTINE D3_Update_Mesh_Info_Crack(isub)
c     Update mesh info of 3D crack surface. Used for local mesh refinement. 
c
c     Info need to be updated:
c          Cr3D_Meshed_Node_in_Ele_Num(Max_Num_Cr_3D,Max_N_Node_3D) 
c          Cr3D_Meshed_Node_in_Ele_Local(Max_Num_Cr_3D,Max_N_Node_3D,3)
c
c     Variables related to 3D crack mesh:
C     real(kind=FT) Crack3D_Meshed_Node(Max_Num_Cr_3D, Max_N_Node_3D, 3)
C     integer Cr3D_Meshed_Node_in_Ele_Num(Max_Num_Cr_3D, Max_N_Node_3D)
C     real(kind=FT) Cr3D_Meshed_Node_in_Ele_Local(Max_Num_Cr_3D, Max_N_Node_3D, 3)
C     integer Crack3D_Meshed_Node_num(Max_Num_Cr_3D)
C     integer Crack3D_Meshed_Ele(Max_Num_Cr_3D, Max_N_Node_3D, 3)
c     integer Crack3D_Meshed_Ele_Attri(Max_Num_Cr_3D, Max_N_Node_3D, 5)
c     real(kind=FT)  Crack3D_Meshed_Ele_Nor_Vector(Max_Num_Cr_3D, Max_N_Node_3D, 3)
c     real(kind=FT)  Crack3D_Meshed_Node_Nor_Vector(Max_Num_Cr_3D, Max_N_Node_3D, 3)
C     integer Crack3D_Meshed_Ele_num(Max_Num_Cr_3D)
C     integer Crack3D_Meshed_Outline(Max_Num_Cr_3D, Max_N_Node_3D, 4)
                                                                        ! Data 1 is the first point on the boundary line of the crack front edge
                                                                        ! Data 2 is the second point on the boundary line of the crack front edge
                                                                        ! Data 3 corresponds to the element number
                                                                        ! Data 4 is used to mark whether the two points of the boundary line are allowed to extend,
                                                                        ! extending in very small steps (2021-08-20)
C     integer Crack3D_Meshed_Outline_num(Max_Num_Cr_3D)
C     real(kind=FT) Crack3D_Meshed_Vertex_x_Vector(Max_Num_Cr_3D, Max_N_Node_3D,3)
C     real(kind=FT) Crack3D_Meshed_Vertex_y_Vector(Max_Num_Cr_3D, Max_N_Node_3D, 3)
C     real(kind=FT) Crack3D_Meshed_Vertex_z_Vector(Max_Num_Cr_3D, Max_N_Node_3D, 3)
c
c     Fristly written on 2021-08-08.

c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Elem_Area_Vol
      use Global_Crack_Common
      use Global_Crack_3D
      
      use Global_Cal_Ele_Num_by_Coors_3D  
      
c     --------------------------
C     Variable Type Declaration
c     --------------------------
      implicit none
      integer,intent(in)::isub
      real(kind=FT) c_Point(3)
      integer i_P,OUT_Elem,i_C
      real(kind=FT) c_Kesi,c_Yita,c_Zeta
      integer Ele_Num_Cache
      
      !Loop through crack.
      do i_C =1,num_Crack  
          Ele_Num_Cache = 1
          !Loop through crack mesh node.
          do i_P =1,Crack3D_Meshed_Node_num(i_C)
              c_Point(1:3) = Crack3D_Meshed_Node(i_C)%row(i_P,1:3) 
              !Get element number of crack node and save it.
              call Cal_Ele_Num_by_Coors_3D(c_Point(1),c_Point(2),
     &                                     c_Point(3),
     &                                     Ele_Num_Cache,OUT_Elem)
              Cr3D_Meshed_Node_in_Ele_Num(i_C)%row(i_P) = OUT_Elem
              !Get local coordinate system of crack node and save it.                  
              call Cal_KesiYita_by_Coor_3D(c_Point(1:3),OUT_Elem,
     &                                     c_Kesi,c_Yita,c_Zeta)
              Cr3D_Meshed_Node_in_Ele_Local(i_C)%row(i_P,1:3) = 
     &                                    [c_Kesi,c_Yita,c_Zeta]
          enddo
      enddo
 
      
      RETURN
      END SUBROUTINE D3_Update_Mesh_Info_Crack
