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
 
      SUBROUTINE Get_Ele_Composite_Mat_Rot_Matrix
c     Calculate the coordinate rotation matrix of a 3D composite material element.
c     2020-03-19.
c     Ref: ANSYS Theory Manual. ANSYS 15.0. Formula 13.140. PDF version, page 471.

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
      integer i_E,Coor_Type,mat_num
      real(kind=FT) Vec_x(3),Vec_y(3),Vec_z(3)
      real(kind=FT) a11,a12,a13,a21,a22,a23,a31,a32,a33
      real(kind=FT) E1,E2,E3,v12,v23,v13,G12,G23,G13,temp
      real(kind=FT) T(6,6),El_Center(3),Cylinder_Center(3),Cyl_z_coor(3)
      real(kind=FT) c_Dis,c_foot_point(3),Cyl_zaxis_A(3),Cyl_zaxis_B(3)
      
      print *,'    Calculating material coordinate rotation matrix...'
      ALLOCATE(Ele_ComMat_RotMatrix(Num_Elem,6,6)) 
      Ele_ComMat_RotMatrix(1:Num_Elem,1:6,1:6) = ZR
      
c     ------------------------------
C     Main code: loop over elements  
c     ------------------------------
      do i_E =1,Num_Elem
          mat_num = Elem_Mat(i_E)
          El_Center = Elem_Centroid(i_E,1:3)
          !if the mat is composite material.
          if (Material_Type(mat_num)==5)then
            !*********************************
            !Cartesian mat coordinate system.
            !********************************* 
            if(Material_Para_Added(mat_num,11)<TWO)then
                Vec_x = Mat_Cartesian_Coor_Vector_x(mat_num,1:3)
                Vec_y = Mat_Cartesian_Coor_Vector_y(mat_num,1:3)
                call Vector_Normalize(3,Vec_x)   
                call Vector_Normalize(3,Vec_y) 
                ! Vector_x crossed with Vector_y yields Vector_z
                call Vector_Cross_Product_3(Vec_x,Vec_y,Vec_z)  
            !*******************************
            !Cylinder mat coodinate system.
            !*******************************
            elseif(Material_Para_Added(mat_num,11)>TWO)then
                Vec_z = Mat_cylinder_Coor_Vector_z(mat_num,1:3)
                Cylinder_Center = Mat_cylinder_Coor_Center(mat_num,1:3)
                !Get the foot_point of El_Center to the z-axis of the Cylinder coodinate system.
                Cyl_zaxis_A = Cylinder_Center
                Cyl_zaxis_B = Cylinder_Center + ONE*Vec_z
                call Tool_Dis_Point_to_Line_AB_3D(El_Center,
     &                  Cyl_zaxis_A,Cyl_zaxis_B,c_Dis,c_foot_point,1)
                Vec_y = El_Center - c_foot_point
                call Vector_Normalize(3,Vec_y) 
                !Get vector x by performing cross-product.
                call Vector_Cross_Product_3(Vec_y,Vec_z,Vec_x)  
            endif
            a11 = Vec_x(1);a12 = Vec_x(2);a13 = Vec_x(3)
            a21 = Vec_y(1);a22 = Vec_y(2);a23 = Vec_y(3)
            a31 = Vec_z(1);a32 = Vec_z(2);a33 = Vec_z(3)            
            T(1,1) = a11**2;  T(1,2) = a12**2;  T(1,3) = a13**2
            T(1,4) = a11*a12; T(1,5) = a12*a13; T(1,6) = a11*a13
            T(2,1) = a21**2;  T(2,2) = a22**2;  T(2,3) = a23**2
            T(2,4) = a21*a22; T(2,5) = a22*a23; T(2,6) = a21*a23
            T(3,1) = a31**2;  T(3,2) = a32**2;  T(3,3) = a33**2
            T(3,4) = a31*a32; T(3,5) = a32*a33; T(3,6) = a31*a33  
            T(4,1) = TWO*a11*a21
            T(4,2) = TWO*a12*a22
            T(4,3) = TWO*a13*a23
            T(4,4) = a11*a22+a12*a21
            T(4,5) = a12*a23+a13*a32
            T(4,6) = a11*a23+a13*a21
            T(5,1) = TWO*a21*a31
            T(5,2) = TWO*a22*a32
            T(5,3) = TWO*a23*a33
            T(5,4) = a21*a32+a22*a31
            T(5,5) = a22*a33+a23*a32
            T(5,6) = a21*a33+a23*a31
            T(6,1) = TWO*a11*a31
            T(6,2) = TWO*a12*a32
            T(6,3) = TWO*a13*a33
            T(6,4) = a11*a32+a12*a31
            T(6,5) = a12*a33+a13*a32
            T(6,6) = a11*a33+a13*a31  
            Ele_ComMat_RotMatrix(i_E,1:6,1:6)=T(1:6,1:6)            
          endif
      end do
      
      RETURN
      END SUBROUTINE Get_Ele_Composite_Mat_Rot_Matrix
