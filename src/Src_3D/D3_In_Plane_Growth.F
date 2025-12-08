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
 
      subroutine D3_In_Plane_Growth(point_inout,i_C) 
c     In-plane extension. Correct the extension point at the crack front
c     vertex to the foot of the perpendicular from that point to the plane of the crack.
c     Added date: 2022-04-18, NEWFTU2022041801


c     **************
c     Public Module
c     **************
      use Global_Crack_3D
      use Global_Model
      use Global_Common
      
c     *********************
c     Variable Declaration
c     *********************
      implicit none
      real(kind=FT),intent(inout)::point_inout(3)
      integer,intent(in)::i_C
      real(kind=FT) c_Point(3),c_PER(3),c_Distance
      logical c_Yes_PER_in,c_Yes_PER_on
      real(kind=FT) c_Tri_P1(3),c_Tri_P2(3),c_Tri_P3(3)
      integer c_CalP_1,c_CalP_2,c_CalP_3
      integer Crack_Node1,Crack_Node2,Crack_Node3
      
#ifndef Silverfrost
      
      c_Point = point_inout
      

      !*****************************************
      !       OPTION - 1
      ! There is a bug. For example, when three
      ! fluid nodes are on a straight line.
      !*****************************************
C     Select the three nodes of the first fluid element of the crack as the projection surface.
C     c_CalP_1 = Cracks_FluidEle_CalP_3D(i_C,1,1) 
C     c_CalP_2 = Cracks_FluidEle_CalP_3D(i_C,1,2) 
C     c_CalP_3 = Cracks_FluidEle_CalP_3D(i_C,1,3) 
C     Coordinates of 3 fluid nodes
C     c_Tri_P1 = Cracks_CalP_Coors_3D(i_C,c_CalP_1,1:3) 
C     c_Tri_P2 = Cracks_CalP_Coors_3D(i_C,c_CalP_2,1:3) 
C     c_Tri_P3 = Cracks_CalP_Coors_3D(i_C,c_CalP_3,1:3) 

      !*****************************************
      !       OPTION - 2 
      ! OPTION - 1 has a bug. BUGFIX2022070801.
      !*****************************************
      ! Select the three nodes of the first discrete crack surface element 
      ! of the crack as the projection plane. 2022-07-08.
      Crack_Node1 = Crack3D_Meshed_Ele(i_C)%row(1,1)
      Crack_Node2 = Crack3D_Meshed_Ele(i_C)%row(1,2)
      Crack_Node3 = Crack3D_Meshed_Ele(i_C)%row(1,3)
      c_Tri_P1 = Crack3D_Meshed_Node(i_C)%row(Crack_Node1,1:3)
      c_Tri_P2 = Crack3D_Meshed_Node(i_C)%row(Crack_Node2,1:3)
      c_Tri_P3 = Crack3D_Meshed_Node(i_C)%row(Crack_Node3,1:3)
               
      
      ! Calculate the foot of the perpendicular projection
      call Tool_Dis_Point_to_3D_Tri(c_Point,
     &                        c_Tri_P1,c_Tri_P2,c_Tri_P3,
     &                  c_Distance,c_PER,c_Yes_PER_in,c_Yes_PER_on)
      ! Change the Final_Point coordinates to the foot of the perpendicular
      point_inout = c_PER
      
#endif     
      
#ifdef Silverfrost
      print *,'    ERROR :: Silverfrost compiler failed to compile '//
     &        'D3_In_Plane_Growth.F!'
      print *,'             In D3_In_Plane_Growth.F.'
      call Warning_Message('S',Keywords_Blank)
#endif  
      RETURN
      END SUBROUTINE D3_In_Plane_Growth