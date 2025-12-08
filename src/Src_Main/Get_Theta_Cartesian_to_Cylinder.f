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
 
      SUBROUTINE Get_Theta_Cartesian_to_Cylinder
c     Calculate the cylindrical coordinate angle theta from the Cartesian coordinates of each node in each element, used for transforming stress or strain tensors.
c     Ref: \theory_documents\024.1 Stress Tensor and Strain Tensor Transformation from Cartesian Coordinate System to Cylindrical Coordinate System A.26-2021-09-10.pdf
c     \theory_documents\024.2 Stress Tensor and Strain Tensor Transformation from Cartesian Coordinate System to Cylindrical Coordinate System 2.7-2021-09-10.html
c     2021-09-10.

c     ----------------------------
c     Read public variable module
c     ----------------------------
      use Global_Float_Type
      use Global_Common
      use Global_Model
      use Global_Dynamic
      use Global_Material
      
c     --------------------------
C     Variable Type Declaration
c     --------------------------
      implicit none
      integer i_E,i_N
      real(kind=FT) Vec_x(3),Vec_y(3),Vec_z(3)
      real(kind=FT) Cylinder_Center(3)
      real(kind=FT) Node_Coor(3)
      real(kind=FT) c_Dis,c_foot_point(3),Cyl_zaxis_A(3),Cyl_zaxis_B(3)
      real(kind=FT) Tem_Vector(3),c_angle
      real(kind=FT) dot,v1(3),v2(3),n(3),det,cross_v1v2(3)
      integer c_Node
      
      print *,'    Calculating theta angle '
     &//      'from Cartesian CS to cylinder CS of each element...'
      do i_E =1,Num_Elem
          do i_N = 1,8
            ! element center.
            !El_Center = Elem_Centroid(i_E,1:3)
            
            ! Node Coordinates
            Node_Coor(1) = G_X_NODES(i_N,i_E)
            Node_Coor(2) = G_Y_NODES(i_N,i_E)
            Node_Coor(3) = G_Z_NODES(i_N,i_E)
            
            c_Node = G_NN(i_N,i_E)
            
            ! The three coordinate axes of the cylindrical coordinate system (initial Cartesian)
            Vec_x = Post_cylinder_Coor_Vector_x(1:3)
            Vec_y = Post_cylinder_Coor_Vector_y(1:3)
            Vec_z = Post_cylinder_Coor_Vector_z(1:3)      
            
            ! Vector_x crossed with Vector_y yields Vector_z
            !call Vector_Cross_Product_3(Vec_x,Vec_y,Vec_z)   
            
            call Vector_Normalize(3,Vec_x)   
            call Vector_Normalize(3,Vec_y)         
            call Vector_Normalize(3,Vec_z) 
            
            ! Center coordinates of the cylinder
            Cylinder_Center = Post_cylinder_Coor_Center(1:3)
            
            ! The foot of the perpendicular from the calculation node to the z-axis of the column coordinate
            ! system.
            Cyl_zaxis_A = Cylinder_Center
            Cyl_zaxis_B = Cylinder_Center + ONE*Vec_z
            call Tool_Dis_Point_to_Line_AB_3D(Node_Coor,
     &              Cyl_zaxis_A,Cyl_zaxis_B,c_Dis,c_foot_point,1)
     
             
            ! Obtain the vector Tem_Vector pointing from the foot to the node
            Tem_Vector = Node_Coor - c_foot_point
            call Vector_Normalize(3,Tem_Vector) 
            
            ! Calculate the angle theta
            !-----------------------------------------------------------------------------
            ! ----OPTION 1: Incorrect, because the obtained theta is between 0 and 180---
            !-----------------------------------------------------------------------------
            ! Calculating the angle between Tem_Vector and the x-axis of the cylindrical coordinate system gives
            ! the theta angle (0-180, incorrect)
C           call Tool_Angle_of_Vectors_a_and_b_3D(Tem_Vector,Vec_x,
C &                                            c_angle,1)
            !-----------------------------------------------------------------------------
            ! ----OPTION 2: Incorrect, because the obtained theta is between 0 and 180---
            !-----------------------------------------------------------------------------
            
            ! Using the GEOMETRY library to calculate the angle (0-180, incorrect).
C           call lines_exp_angle_3d (c_foot_point,Node_Coor,
C    &                     Cylinder_Center,Cylinder_Center+Vec_x*ONE,
C    &                     c_angle)
            
            !-----------------------------------------------------------------
            ! ----OPTION 3: Yes, the obtained theta is between 0-360      ---
            !Ref: \theory_documents\025 How to find an angle 
            !      in range(0, 360) between 2 vectors-2021-09-10.pdf
            !
            !-----------------------------------------------------------------
            v1 = Vec_x
            v2 = Tem_Vector
            n  = Vec_z
            dot = dot_product(v1,v2)
            
            call Vector_Cross_Product_3(v1,v2,cross_v1v2)   
            det = dot_product(n,cross_v1v2)
            
            c_angle = atan2(det, dot)
            
            Theta_Cartesian_to_Cylinder(i_E,i_N) = c_angle
            Theta_Cartesian_to_Cylinder_Node(c_Node) = c_angle
            
            
          end do
      enddo
      RETURN
      END SUBROUTINE Get_Theta_Cartesian_to_Cylinder
