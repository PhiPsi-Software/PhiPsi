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
 
      SUBROUTINE D3_Get_Crack_Mesh_Outline_Clockwise(i_C,
     &           Outline_num,c_Outline,Cros_Product_Vector)
c     Determine the clockwise or counterclockwise direction of the outline 
c     using the cross product. Ref: My PhiPsi Development Notebook V1-P41.
c     2022-07-13.

c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Elem_Area_Vol
      use Global_Crack_3D
      
c     --------------------------
C     Variable Type Declaration
c     --------------------------
      implicit none
      integer,intent(in)::i_C,Outline_num
      integer,intent(in)::c_Outline(Outline_num,2)
      real(kind=FT),intent(out)::Cros_Product_Vector(3)
      integer i_outline,c_P1,c_P2,c_P3
      real(kind=FT) c_Vector_1(3),c_Vector_2(3)
      real(kind=FT) c_P1_Coor(3),c_P2_Coor(3),c_P3_Coor(3)
      real(kind=FT) c_Cros_Product(3)
      real(kind=FT) Sum_Cros_Product(3)
      real(kind=FT) c_Center(3) 

c     ---------------------
C     Main program section
c     ---------------------
      !Cros_Product_Vector = ZR 
      Sum_Cros_Product(1:3) =  ZR
      
      ! Obtain the centroid of the boundary discrete points.
      c_Center(1) = sum(Crack3D_Meshed_Node(i_C)%row(
     &                  c_Outline(1:Outline_num,1),1))/Outline_num
      c_Center(2) = sum(Crack3D_Meshed_Node(i_C)%row(
     &                  c_Outline(1:Outline_num,1),2))/Outline_num
      c_Center(3) = sum(Crack3D_Meshed_Node(i_C)%row(
     &                  c_Outline(1:Outline_num,1),3))/Outline_num   
     
      do i_outline =1,Outline_num
          c_P1 = c_Outline(i_outline,1) 
          c_P2 = c_Outline(i_outline,2) 
          !c_P3 = c_Outline(i_outline+1,2) 
          c_P1_Coor = Crack3D_Meshed_Node(i_C)%row(c_P1,1:3)
          c_P2_Coor = Crack3D_Meshed_Node(i_C)%row(c_P2,1:3)
          !c_P3_Coor = Crack3D_Meshed_Node(i_C,c_P3,1:3)
          c_Vector_1 = c_P2_Coor - c_P1_Coor
          !c_Vector_2 = c_P3_Coor - c_P2_Coor
          c_Vector_2 = c_Center - c_P1_Coor
          call Vector_Normalize(3,c_Vector_1)  
          call Vector_Normalize(3,c_Vector_2)
          call Vector_Cross_Product_3(c_Vector_1,c_Vector_2,
     &                                c_Cros_Product(1:3)) 
          call Vector_Normalize(3,c_Cros_Product)  
          Sum_Cros_Product(1:3) = Sum_Cros_Product(1:3) +c_Cros_Product
      enddo
      Cros_Product_Vector = Sum_Cros_Product(1:3)/Outline_num
      
      call Vector_Normalize(3,Cros_Product_Vector)  
      
      RETURN
      END SUBROUTINE D3_Get_Crack_Mesh_Outline_Clockwise
