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
 
      function Tool_Function_Dis_Point_to_3D_Tri_v2(Point,
     &                            Tri_P1,Tri_P2,Tri_P3)
C     Calculating the distance from a point to a triangle in space, Reference: Jones_1995_3D Distance from a Point to a Triangle
c     Tool_Function_Dis_Point_to_3D_Tri_v2 is a simplified version of Tool_Dis_Point_to_3D_Tri, only calculating the unsigned distance.
c     2022-05-05. NEWFTU2022050502.

      !......................
      ! Variable Declaration
      !......................
      use Global_Float_Type     
      implicit none
      real(kind=FT),intent(in)::Point(3),
     &                          Tri_P1(3),Tri_P2(3),Tri_P3(3)
      real(kind=FT) P1P2(3),P1P3(3),P1P0(3),Np(3),cosa
      real(kind=FT) abs_P1P0,abs_Np
      real(kind=FT) Tool_Function_Dis_Point_to_3D_Tri_v2
      
      Tool_Function_Dis_Point_to_3D_Tri_v2   =  ZR
      
      P1P2 = Tri_P2 - Tri_P1
      P1P3 = Tri_P3 - Tri_P1
      
      
      call Vector_Cross_Product_3(P1P2,P1P3,Np)  
      
      P1P0 = Point - Tri_P1
      abs_P1P0 = sqrt(P1P0(1)**2 + P1P0(2)**2 + P1P0(3)**2)
      abs_Np   = sqrt(Np(1)**2   + Np(2)**2   + Np(3)**2)
      cosa = dot_product(P1P0,Np)/abs_P1P0/abs_Np 
      
      
      Tool_Function_Dis_Point_to_3D_Tri_v2 = abs(abs_P1P0 * cosa)
      

      end function  Tool_Function_Dis_Point_to_3D_Tri_v2                 
