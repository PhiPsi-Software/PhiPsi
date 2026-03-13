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
 
      subroutine Tool_Offset_Point_A_to_Point_B_3D(A,B,delta_L,new_A)
      ! Offset point A towards B by delta_L.
      !2022-05-05.
      !
      !                A                       B
      !
      !                o-----------------------o
      !
      !                delta_L > 0:
      !
      !                |<---delta_L--->|
      !                                o-------o
      !              
      !              
      ! Reference: http://wenku.baidu.com/view/2613ac1a227916888486d7cb.html, Vectors in Three-Dimensional
      ! Space 2011
      
      use Global_Float_Type      
      implicit none
      real(kind=FT),intent(in) ::A(3),B(3),delta_L
      real(kind=FT),intent(out)::new_A(3)
      real(kind=FT) Vector_AB(3),norm_Vector_AB
      real(kind=FT) theta1,theta2,theta3
      
      Vector_AB = B(1:3) - A(1:3) 
      
      norm_Vector_AB = sqrt(Vector_AB(1)**2 + 
     &                      Vector_AB(2)**2 + Vector_AB(3)**2)
      theta1 = acos(Vector_AB(1)/norm_Vector_AB)
      theta2 = acos(Vector_AB(2)/norm_Vector_AB) 
      theta3 = acos(Vector_AB(3)/norm_Vector_AB) 
      
      new_A(1)= A(1)+delta_L*cos(theta1)
      new_A(2)= A(2)+delta_L*cos(theta2)
      new_A(3)= A(3)+delta_L*cos(theta3)
      
      return 
      end SUBROUTINE Tool_Offset_Point_A_to_Point_B_3D                         
